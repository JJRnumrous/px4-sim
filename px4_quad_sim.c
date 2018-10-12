#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>

#include "quad_parameters.h"
#include "utils.h"
#include "px4_sim_communication.h"
#include "px4_quad_sim.h"

// Function prototypes.
void testQuadDynamics();
void test6DOF();
void printSensors(Quad *quad);
void print6DOF(Quad quad);
void sighandler(int sig);

// Function definitions.
int init_sim()
{
    int i;
    double inertia[3];
    double c_d[3];

    // Initialise quad model
    inertia[0] = I_XX;
    inertia[1] = I_YY;
    inertia[2] = I_ZZ;

    for(i = 0 ; i < 3 ; i++)
        c_d[i] = C_D;

    init_quad(&quad, MASS, inertia, D_ARM, R_D, R_LD, c_d, THRUST_TC, THRUST_HOVER_NORM * THRUST_MAX_FORCE);
    init_quad_sensors(&quad, GPS_EPH, GPS_EPV, GPS_FIX, GPS_NUM_SATS, GPS_LAT_LON_NOISE, GPS_ALT_NOISE, GPS_SPEED_NOISE, IMU_ACC_NOISE, IMU_GYRO_NOISE, MAG_DECL, MAG_INCL, MAG_SCALE, MAG_NOISE, BARO_TEMP);

    // Initialise PX4 sim
    return init_px4_sim(SENSOR_FREQ, GPS_FREQ, HIL, THRUST_HOVER_NORM, THRUST_MAX_FORCE);
}

int advance_sim(uint64_t time_usec, double dt, double wind_vel_e[3])
{
    int i;
    int result;
    double q[4];
    double thrust_commands[4];

    get_att_quaternion(&quad, q);

    // Send HIL MAVLink messages
    result = send_hil_messages(time_usec, q, quad.state.euler_rates, quad.sensors.gps.lat_lon_alt, quad.sensors.gps.gps_speed, quad.sensors.gps.ground_speed, quad.sensors.gps.cog, quad.sensors.gps.eph, quad.sensors.gps.epv, quad.sensors.gps.fix, quad.sensors.gps.visible_sats, quad.sensors.imu.acc, quad.sensors.imu.gyro, quad.sensors.mag.mag_field, quad.sensors.baro.pressure, quad.sensors.baro.temperature);
    if(result < 0)
        return result;

    // Gaurd against large dt's
    // if(dt >= 2 * (double)(1.0 / SENSOR_FREQ))
    // {
    //     dt = 2.0 * (double)(1.0 / SENSOR_FREQ);
    //     printf("changed dt: %2.3f ms\n", dt);
    // }

    get_thrust_commands_force(thrust_commands);
    update_quad(&quad, thrust_commands, wind_vel_e, dt);
}

int main()
{
    /*
    int i;
    double dt;
    uint64_t cur_time;
    uint64_t prev_time;
    uint64_t sensor_update;
    double wind_vel_e[3];
    int result;
    double excess_time;

    signal(SIGINT, sighandler); // Exit simulation if Ctrl-C is pressed

    // Connect to PX4 and initialize physics sim
    printf("Connecting to PX4...\n");
    result = init_sim();
    if(result < 0)
    {
        printf("Failed! Result: %d\n", result);
        return 0;
    }
    printf("Connected!\n");

    // Initialise variables
    cur_time = get_time_usec();
    prev_time = cur_time;
    sensor_update = cur_time;

    for(i = 0 ; i < 3 ; i++)
        wind_vel_e[i] = 0.0f;

    // Loop
    while(true)
    {
        // Update the physics sim with a fixed time
        cur_time = get_time_usec();
        if((int64_t)(cur_time - sensor_update) >= 0)
        {
            dt = (double)((cur_time - prev_time) * 1e-6);
            // printf("dt: %2.3f ms\n", dt * 1e3);
            advance_sim(cur_time, dt, wind_vel_e);
            prev_time = cur_time;

            // Calculate next time the physics sim should update
            sensor_update = cur_time + (uint64_t)(1000000.0 / SENSOR_FREQ);
        }
        
        // Poll for MAVLink messages: receive actuator controls from PX4 and when in HIL - send messages from autopilo to GCS and vice-versa.
        result = pollMavlinkMessage();
        if(result < 0)
            return result;
        usleep(10); // 10 us
    }
*/
   test6DOF();
    return 0;
}

void sighandler(int sig)
{
   printf("\nExit simulation...\n");
   disconnect_sim();
   exit(1);
}

// Only for testing:
void testQuadDynamics()
{
    int i;
    double sim_time = 2;
    double wind_vel_e[3];
    double inertia[3];
    double c_d[3];
    double thrust_commands[4];


    // Initialise quad model
    inertia[0] = I_XX;
    inertia[1] = I_YY;
    inertia[2] = I_ZZ;

    for(i = 0 ; i < 3 ; i++)
        c_d[i] = C_D;

    init_quad(&quad, MASS, inertia, D_ARM, R_D, R_LD, c_d, THRUST_TC, THRUST_HOVER_NORM * THRUST_MAX_FORCE);
    init_quad_sensors(&quad, GPS_EPH, GPS_EPV, GPS_FIX, GPS_NUM_SATS, GPS_LAT_LON_NOISE, GPS_ALT_NOISE, GPS_SPEED_NOISE, IMU_ACC_NOISE, IMU_GYRO_NOISE, MAG_DECL, MAG_INCL, MAG_SCALE, MAG_NOISE, BARO_TEMP);

    // Initialise wind and thrust
    for(i = 0 ; i < 3 ; i++)
        wind_vel_e[i] = 0.0f;
    
    thrust_commands[0] = 0.3 * 0.5 * 15 * GRAVITY;
    thrust_commands[1] = 0.6 * 0.5 * 15 * GRAVITY;
    thrust_commands[2] = 0.7 * 0.5 * 15 * GRAVITY;
    thrust_commands[3] = 0.6 * 0.5 * 15 * GRAVITY;

    for(i = 0 ; i < (sim_time*SENSOR_FREQ) ; i++)
    {
        update_quad(&quad, thrust_commands, wind_vel_e, (double)(1.0 / SENSOR_FREQ));
        printSensors(&quad);
    }
}

/* *************************************************************************************************************************************************
* Function:     For testing Runga Kutta and 6DOF equations
*
*
************************************************************************************************************************************************** */
void test6DOF()
{
    int i;
    double sim_time = 1000.0/SENSOR_FREQ;
    double inertia[3];
    double c_d[3];    
    
    // Initialise quad model
    inertia[0] = I_XX;
    inertia[1] = I_YY;
    inertia[2] = I_ZZ;

    for(i = 0 ; i < 3 ; i++)
        c_d[i] = C_D;

    init_quad(&quad, MASS, inertia, D_ARM, R_D, R_LD, c_d, THRUST_TC, THRUST_HOVER_NORM * THRUST_MAX_FORCE);
   
    double forces[3], moments[3];
    forces[0] = 1;
    forces[1] = 2;
    forces[2] = 3;
    moments[0] = 4;
    moments[1] = 5;
    moments[2] = 6;

    
    for(i=1; i <= (sim_time*SENSOR_FREQ) ; i++)
    {
        six_dof((1.0/SENSOR_FREQ),&quad,forces,moments);
        //printf("Forces %d: \t %f \n",i,quad.state.omega_b[0]);
        print6DOF(quad);
    }    
    
    
}

void print6DOF(Quad quad)
{
    int i = 0;
    printf("\nVel_body:\t");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad.state.vel_b[i]);
    printf("\n");

    printf("Omega_body:\t");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad.state.omega_b[i]);
    printf("\n");
    
    printf("Acc_body:\t");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad.state.acc_b[i]);
    printf("\n");

    printf("Euler_rates:\t");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad.state.euler_rates[i]);
    printf("\n");

    printf("Euler:\t\t");
    for(i = 0 ; i < 3 ; i++)

    printf("%4.6f\t", quad.state.euler[i]*180.0/M_PI -360.0);
    printf("\n");
    printf("Pos_e:\t\t");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad.state.pos_e[i]);
    printf("\n");
   
}


void printSensors(Quad *quad)
{
    int i;

    printf("Accelerometer:\n");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad->sensors.imu.acc[i]);
    printf("\n");

    printf("Gyro:\n");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad->sensors.imu.gyro[i]);
    printf("\n");

    printf("Barometer:\n");
    printf("%4.6f\t", quad->sensors.baro.pressure);
    printf("\n");

    printf("Magnetometer:\n");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad->sensors.mag.mag_field[i]);
    printf("\n");

    printf("GPS:\n");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad->sensors.gps.lat_lon_alt[i]);
    printf("\n");

    printf("GPS Speed:\n");
    for(i = 0 ; i < 3 ; i++)
        printf("%4.6f\t", quad->sensors.gps.gps_speed[i]);
    printf("\n\n");
}