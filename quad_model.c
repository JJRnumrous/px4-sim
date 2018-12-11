#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "quad_model.h"

//Private variables
double *_Forces, *_Moments, _thrust_command;
Quad *_quad;

// Private function prototypes
void update_sensors(Quad *quad);

double D_U(double states[]);
double D_V(double states[]);
double D_W(double states[]);
double D_P(double states[]);
double D_Q(double states[]);
double D_R(double states[]);

// double D_phi(double states[]);
// double D_theta(double states[]);
// double D_psi(double states[]);
double D_q0(double state[]);
double D_q1(double state[]);
double D_q2(double state[]);
double D_q3(double state[]);
double D_N(double states[]);
double D_E(double states[]);
double D_D(double states[]);

double D_T(double states[]);

// Public function definitions
void init_quad(Quad *quad, double mass, double inertia[3], double d, double r_d, double c_d[3], double thrust_tc, double throttle_hover, double init_yaw)
{
    int i;
    int j;

    quad->mass = mass;
    quad->d = d;
    quad->r_d = r_d;    
    quad->thrust_tc = thrust_tc;
    
    for(i = 0 ; i < 4 ; i++){
        quad->thrust[i] = throttle_hover / 4;
        quad->state.quat_rate[i] = 0.0f;
        quad->state.quat[i] = 0.0f;
    }
    quad->state.quat[0] = 1.0f;
    for(i = 0 ; i < 3 ; i++)
    {
        quad->inertia[i] = inertia[i];
        quad->c_d[i] = c_d[i];

        quad->state.acc_b[i] = 0.0f;
        quad->state.vel_b[i] = 0.0f;
        quad->state.vel_e[i] = 0.0f;
        quad->state.pos_e[i] = 0.0f;
        quad->state.omega_b[i] = 0.0f;
        // quad->state.euler_rates[i] = 0.0f;
        quad->state.euler[i] = 0.0f;

        for(j = 0 ; j < 3 ; j++)
        {
            if(i == j)
            {
                quad->state.dcm_be[i][j] = 1.0f;
                continue;
            }
            quad->state.dcm_be[i][j] = 0.0f;
        }
    }
    quad->state.euler[2] = init_yaw;
    euler_to_quat(quad->state.euler, quad->state.quat);
}

void init_quad_sensors(Quad *quad, double eph, double epv, double fix, double visible_sats, double lat_lon_noise_std_dev, double alt_noise_std_dev, double speed_noise_std_dev, double acc_noise_std_dev, double gyro_noise_std_dev, double mag_decl, double mag_incl, double mag_scale, double mag_noise_std_dev, double temperature)
{
    init_gps(&(quad->sensors.gps), eph, epv, fix, visible_sats, lat_lon_noise_std_dev, alt_noise_std_dev, speed_noise_std_dev);
    init_imu(&(quad->sensors.imu), acc_noise_std_dev, gyro_noise_std_dev);
    init_mag(&(quad->sensors.mag), mag_decl, mag_incl, mag_scale, mag_noise_std_dev);
    init_baro(&(quad->sensors.baro), temperature);
}

void six_dof(double dt, Quad *quad, double forces[3], double moments[3])
{
    int i,len;
    int tlen = 13;
    double total[tlen], total_integ[tlen], quat_norm;
    // vFunctionCall total_dot[] = {D_U, D_V, D_W, D_P, D_Q, D_R,D_phi, D_theta,D_psi,D_N, D_E, D_D};
    vFunctionCall total_dot[] = {D_U, D_V, D_W, D_P, D_Q, D_R,D_q0, D_q1, D_q2, D_q3, D_N, D_E, D_D};

    QuadState *quad_state;
    quad_state = &(quad->state);

    // initialise local global for 6dof Runga Kutta
    _Forces  = forces;
    _Moments = moments;
    _quad    = quad;

    for(i=0; i<tlen;i++){
        if(i<3)             total[i] = quad_state->vel_b[i];
        if((i>=3)&&(i<6))   total[i] = quad_state->omega_b[i-3];
        if((i>=6)&&(i<10))   total[i] = quad_state->quat[i-6];
        if(i>=10)            total[i] = quad_state->pos_e[i-10];
    }

    integer_rk4(dt, total_dot,total,total_integ,tlen);    
    quat_norm = sqrt(pow(total_integ[6],2) + pow(total_integ[7],2) + pow(total_integ[8],2) + pow(total_integ[9],2));

    quad_state->vel_b[0] = total_integ[0];
    quad_state->vel_b[1] = total_integ[1];
    quad_state->vel_b[2] = total_integ[2];
    quad_state->omega_b[0] = total_integ[3];
    quad_state->omega_b[1] = total_integ[4];
    quad_state->omega_b[2] = total_integ[5];
    quad_state->quat[0]   = (total_integ[6])/quat_norm;
    quad_state->quat[1]   = (total_integ[7])/quat_norm;
    quad_state->quat[2]   = (total_integ[8])/quat_norm;
    quad_state->quat[3]   = (total_integ[9])/quat_norm;
    quad_state->pos_e[0]   = total_integ[10];
    quad_state->pos_e[1]   = total_integ[11];
    quad_state->pos_e[2]   = total_integ[12];

    
    quad_state->quat_rate[0] = D_q0(total_integ);
    quad_state->quat_rate[1] = D_q1(total_integ);
    quad_state->quat_rate[2] = D_q2(total_integ);
    quad_state->quat_rate[3] = D_q3(total_integ);    
   
   
    calc_dcm_be(quad_state->quat, quad_state->dcm_be);
    body_to_earth_rotation(quad_state->dcm_be,quad_state->vel_b, quad_state->vel_e);
    // quad_state->vel_e[0] = D_N(total_integ);
    // quad_state->vel_e[1] = D_E(total_integ);
    // quad_state->vel_e[2] = D_D(total_integ);

    // The body accel and euler rates are updated with most recent states

    // calculate and set the euler rates
    // quad_state->euler_rates[0] = D_phi(total_integ);
    // quad_state->euler_rates[1] = D_theta(total_integ);
    // quad_state->euler_rates[2] = D_psi(total_integ);
    

    //calculate and set new linear body accel for timestep
    quad_state->acc_b[0] = D_U(total_integ);
    quad_state->acc_b[1] = D_V(total_integ);
    quad_state->acc_b[2] = D_W(total_integ);    
}

/* **************************************************************************************************************************************************
* Function:    Derivative Functions for Kinetics
* Inputs:      states = u0, v1, w2, p3, q4, r5
* Output:      Value for derivitive at specific states.
**************************************************************************************************************************************************** */
double D_U(double states[])
{
    return _Forces[0]/_quad->mass + states[1]*states[5] - states[2]*states[4];    
}

double D_V(double states[])
{
    return _Forces[1]/_quad->mass - states[0]*states[5] + states[2]*states[3];
}

double D_W(double states[])
{
    return _Forces[2]/_quad->mass + states[0]*states[4] - states[1]*states[3];
}

double D_P(double states[])
{
    return (_Moments[0] - states[4]*states[5]*(_quad->inertia[2] - _quad->inertia[1]))/_quad->inertia[0];
}

double D_Q(double states[])
{
    return (_Moments[1] - states[3]*states[5]*(_quad->inertia[0] - _quad->inertia[2]))/_quad->inertia[1];
}

double D_R(double states[])
{
    return (_Moments[2] - states[3]*states[4]*(_quad->inertia[1] - _quad->inertia[0]))/_quad->inertia[2];
}
/* **************************************************************************************************************************************************
* Function:    Kinematic Derivative functions
* Inputs:      states = u0, v1, w2, p3, q4, r5, phi6, theta7, psi8, N9, E10, D11
* Output:      Value for derivitive at specific states.
****************************************************************************************************************************************************/
// double D_phi(double states[])
// {
//     //return _quad->state.omega_b[0] + sin(states[0])*tan(states[1])*_quad->state.omega_b[1] + cos(states[0])*tan(states[6])*_quad->state.omega_b[2];
//     return states[3] + sin(states[6])*tan(states[7])*states[4] + cos(states[6])*tan(states[7])*states[5];
// }
// double D_theta(double states[])
// {
//     //return cos(states[0])*_quad->state.omega_b[1] - sin(states[0])*_quad->state.omega_b[2];
//     return cos(states[6])*states[4] - sin(states[6])*states[5];
// }
// double D_psi(double states[])
// {
//     //return (sin(states[0])/cos(states[1]))*_quad->state.omega_b[1] + (cos(states[0])/cos(states[1]))*_quad->state.omega_b[2];
//     return (sin(states[6])/cos(states[7]))*states[4] + (cos(states[6])/cos(states[7]))*states[5];
// }
// double D_N(double states[])
// {
//     //return _quad->state.vel_b[0]*cos(_quad->state.euler[2])*cos(_quad->state.euler[1]) + _quad->state.vel_b[1]*(cos(_quad->state.euler[2])*sin(_quad->state.euler[1])*sin(_quad->state.euler[0]) - sin(_quad->state.euler[2])*cos(_quad->state.euler[0])) + _quad->state.vel_b[2]*(cos(_quad->state.euler[2])*sin(_quad->state.euler[1])*cos(_quad->state.euler[0]) + sin(_quad->state.euler[2])*sin(_quad->state.euler[0]));
//     return states[0]*cos(states[8])*cos(states[7]) + states[1]*(cos(states[8])*sin(states[7])*sin(states[6]) - sin(states[8])*cos(states[6])) + states[2]*(cos(states[8])*sin(states[7])*cos(states[6]) + sin(states[8])*sin(states[6]));
// }
// double D_E(double states[])
// {
//     //return _quad->state.vel_b[0]*sin(_quad->state.euler[2])*cos(_quad->state.euler[1]) + _quad->state.vel_b[1]*( sin(_quad->state.euler[2])*sin(_quad->state.euler[1])*sin(_quad->state.euler[0]) + cos(_quad->state.euler[2])*cos(_quad->state.euler[0])) + _quad->state.vel_b[2]*(sin(_quad->state.euler[2])*sin(_quad->state.euler[1])*cos(_quad->state.euler[0]) - cos(_quad->state.euler[2])*sin(_quad->state.euler[0]));
//     return states[0]*sin(states[8])*cos(states[7]) + states[1]*(sin(states[8])*sin(states[7])*sin(states[6]) + cos(states[8])*cos(states[6])) + states[2]*(sin(states[8])*sin(states[7])*cos(states[6]) - cos(states[8])*sin(states[6]));
// }
// double D_D(double states[])
// {
//     //return _quad->state.vel_b[0]*sin(_quad->state.euler[1]) + _quad->state.vel_b[1]*cos(_quad->state.euler[1])*sin(_quad->state.euler[0]) + _quad->state.vel_b[2]*cos(_quad->state.euler[1])*cos(_quad->state.euler[0]);
//     return 0.0 -states[0]*sin(states[7]) + states[1]*cos(states[7])*sin(states[6]) + states[2]*cos(states[7])*cos(states[6]);
// }

double D_q0(double states[])
{
    return 0.0f - 0.5f*(states[3]*states[7] + states[4]*states[8] + states[5]*states[9]);
}
double D_q1(double states[])
{
    return 0.5f*(states[3]*states[6] - states[4]*states[9] + states[5]*states[8]);
}
double D_q2(double states[])
{
    return 0.5f*(states[3]*states[9] + states[4]*states[6] - states[5]*states[7]);
}
double D_q3(double states[])
{
    return 0.5f*(0.0f - states[3]*states[8] + states[4]*states[7] + states[5]*states[6]);
}
double D_N(double states[])
{
    return states[0]*(pow(states[6],2)+pow(states[7],2)-pow(states[8],2)-pow(states[9],2)) + 2.0f*states[1]*(states[7]*states[8] - states[6]*states[9]) + 2.0f*states[2]*(states[7]*states[9] + states[6]*states[8]);
}
double D_E(double states[])
{
    return 2.0f*states[0]*(states[7]*states[8] + states[6]*states[9]) + states[1]*(pow(states[6],2)-pow(states[7],2)+pow(states[8],2)-pow(states[9],2)) + 2.0f*states[2]*(states[8]*states[9] - states[6]*states[7]);
}
double D_D(double states[])
{
    return 2.0f*states[0]*(states[7]*states[9] - states[6]*states[8]) + 2.0f*states[1]*(states[8]*states[9] + states[6]*states[7]) + states[2]*(pow(states[6],2)-pow(states[7],2)-pow(states[8],2)+pow(states[9],2));
}

double D_T(double states[])
{
    // states = thrust ; thrust_command ; 
    return  (_thrust_command - states[0])/_quad->thrust_tc;
}

void forces_moments_aerodynamic_model(Quad *quad, double wind_vel_e[3], double forces[], double moments[])
{
    int i;
    double scaled_wind[3];
    double wind_vel_b[3];
    double rel_vel[3];

    // Scale wind velocity in inertial frame with altitude.
    for(i = 0 ; i < 3 ; i++)
        scaled_wind[i] = pow(abs(quad->state.pos_e[2])/10.0, ETA_W) * wind_vel_e[i];

    // Convert scaled wind to body axes.
    earth_to_body_rotation(quad->state.dcm_be, scaled_wind, wind_vel_b);

    // Calculate relative velocity.
    for(i = 0 ; i < 3 ; i++)
        rel_vel[i] = -quad->state.vel_b[i] + wind_vel_b[i];

    // Calculate aerodynamic forces.
    for(i = 0 ; i < 3 ; i++)
        forces[i] = 0.5 * RHO * abs(rel_vel[i]) * rel_vel[i] * quad->c_d[i];

    for(i = 0 ; i < 3 ; i++)
        moments[i] = 0;
}

void forces_moments_thrust_model(double dt, double thrust_commands[4], Quad *quad, double forces[], double moments[])
{
    int i;  
    //double thrust_dot;    
    double thrust_integ;
    vFunctionCall thrust_dot[] = {D_T};

    _quad = quad;

    for(i = 0 ; i < 4 ; i++)
    {
        _thrust_command = thrust_commands[i];        
        integer_rk4(dt, thrust_dot, &(quad->thrust[i]), &thrust_integ, 1);
        quad->thrust[i] = thrust_integ;
    }

    forces[0] = 0.0f;
    forces[1] = 0.0f;
    forces[2] =             -(quad->thrust[0] + quad->thrust[1] + quad->thrust[2] + quad->thrust[3]);

    // Cross
    moments[0] = quad->d *   (-quad->thrust[0] + quad->thrust[1] + quad->thrust[2] - quad->thrust[3]);
    moments[1] = quad->d *   ( quad->thrust[0] - quad->thrust[1] + quad->thrust[2] - quad->thrust[3]);
    moments[2] = quad->r_d * ( quad->thrust[0] + quad->thrust[1] - quad->thrust[2] - quad->thrust[3]);

    // Plus
    // moments[0] = quad->d * (quad->thrust[3] - quad->thrust[1]);
    // moments[1] = quad->d * (quad->thrust[0] - quad->thrust[2]);
    // moments[2] = quad->r_d * (-quad->thrust[0]+quad->thrust[1]-quad->thrust[2]+quad->thrust[3]) / quad->r_ld;

}

void forces_moments_gravity_model(Quad *quad, double forces[], double moments[])
{
    int i;

    for(i = 0 ; i < 3 ; i++)
    {
        forces[i] = quad->mass * GRAVITY * quad->state.dcm_be[i][2];
        moments[i] = 0;
    }
}


void update_quad(Quad *quad, double thrust_commands[4], double wind_vel_e[3], double dt)
{
    int i;
    double forces[3], forces_aero[3], forces_thrust[3], forces_gravity[3];
    double moments[3], moments_aero[3], moments_thrust[3], moments_gravity[3];

    forces_moments_aerodynamic_model(quad, wind_vel_e, forces_aero, moments_aero);
    forces_moments_gravity_model(quad, forces_gravity, moments_gravity);
    forces_moments_thrust_model(dt, thrust_commands, quad, forces_thrust, moments_thrust);

    for(i = 0 ; i < 3 ; i++)
    {
        forces[i] = forces_thrust[i] + forces_gravity[i] + forces_aero[i];
        moments[i] = moments_thrust[i] + moments_gravity[i] + moments_aero[i];
    }

    six_dof(dt, quad, forces, moments);
    update_sensors(quad);
}

void update_sensors(Quad *quad)
{
    update_gps(&(quad->sensors.gps), quad->state.pos_e, quad->state.vel_e);
    update_imu(&(quad->sensors.imu), quad->state.acc_b, quad->state.omega_b, quad->state.dcm_be);
    update_mag(&(quad->sensors.mag), quad->state.dcm_be);
    update_baro(&(quad->sensors.baro), quad->sensors.gps.lat_lon_alt[2]);
}