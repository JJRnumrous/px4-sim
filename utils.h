#ifndef UTILS_H_
#define UTILS_H_

#include <inttypes.h>


// Constants
#define GRAVITY 9.80665 // gravity
#define Pb 101325.0 // static pressure at sea level [Pa]
#define Tb 288.15 // standard temperature at sea level [K]
#define Lb -0.0065 // standard temperature lapse rate [K/m]
#define M 0.0289644 // molar mass of Earth's air [kg/mol]
#define R 8.31432 // universal gas constant
#define ETA_W 0.2 // eta
#define RHO 1.225 // air density
#define R_EARTH 6371000.0 // radius of earth (6378137.0)

// Types and structs
typedef struct
{
    double lat_lon_noise_std_dev;
    double alt_noise_std_dev;
    double speed_noise_std_dev;

    double lat_lon_alt[3]; // deg, deg, m
    double gps_speed[3]; // m/s
    double ground_speed; // m/s
    double eph; // GpsSensor HDOP (m)
    double epv; // GpsSensor VDOP (m)
    int fix; // 0-1: No fix ; 2: 2D fix ; 3: 3D fix
    double cog; // course over ground (deg)
    int visible_sats; // number of visible satellites.
} GpsSensor;

typedef struct
{
    double acc_noise_std_dev;
    double gyro_noise_std_dev;

    double acc[3]; // accelerometer (m/(s^2))
    double gyro[3]; // gyroscope (rad/s)
} ImuSensor;

typedef struct
{
    double mag_noise_std_dev;
    double mag_scale;
    double mag_decl; // magnetic declination (deg)
    double mag_incl; // magnetic inclination (deg)

    double mag_field[3]; // magnetic field (gauss)
} MagSensor;

typedef struct
{
    double temperature; // degC
    double pressure; // barometric pressure (Pa)
    double pressure_alt; // altitude (m)
} BaroSensor;

typedef struct
{
    GpsSensor gps;
    ImuSensor imu;
    MagSensor mag;
    BaroSensor baro;
} Sensors;

// Types and structs
typedef struct
{
    double acc_b[3]; // Udot, Vdot, Wdot (m/(s^2))
    double vel_b[3]; // U, V, W (m/s)
    double vel_e[3]; // Ndot, Edot, Ndot (m/s)
    double pos_e[3]; // N, E, D (m)
    double omega_b[3]; // P, Q, R (rad/s)
    double euler_rates[3]; // phi_dot, theta_dot, psi_dot (rad/s)
    double euler[3]; // phi, theta, psi (rad)
    double dcm_be[3][3]; // direct cosine matrix
} QuadState;

typedef struct
{
    double mass;
    double inertia[3]; // Ixx, Iyy, Izz

    double d; // distance from center to motor
    double r_d; // virtual yaw moment arm
    double r_ld; // lift-to-drag ratio
    double thrust_tc; // thrust time constant

    double c_d[3]; // aerodynamic drag coefficient

    // Motor configuration:
    /*
     *     1(N)
     * 4   x(D)  2(E)
     *     3
     */
    double thrust[4]; // thrusts of motors

    QuadState state;
    Sensors sensors;
} Quad;



// Utility functions
// Prototypes
double deg2rad(double deg);
double rad2deg(double rad);
double zero_mean_noise(double std_dev);
double rand_gauss(double mean, double std_dev);

void calc_dcm_be(double euler[3], double dcm_be[3][3]);
void body_to_earth_rotation(double dcm_be[3][3], double rotate_from[3], double rotate_to[3]);
void earth_to_body_rotation(double dcm_be[3][3], double rotate_from[3], double rotate_to[3]);

double wrap_angle_2pi(double angle);
double wrap_angle_pi(double angle);

double integrate(double sum, double val, double dt);
double integrate_euler(double sum, double val, double dt);
typedef double (* vFunctionCall)(double[6], Quad*);
void RK4(double dt, vFunctionCall F1,vFunctionCall F2,vFunctionCall F3,vFunctionCall F4,vFunctionCall F5,vFunctionCall F6, double FM[6], Quad *quad,double *ret);
void state_adder(Quad *quad, double k,double l,double m,double n,double o,double p,double devider);
uint64_t get_time_usec();

#endif