#ifndef QUAD_MODEL_H_
#define QUAD_MODEL_H_


#include "sensors_model.h"
#include "utils.h"

// Types and structs
typedef struct
{
    double acc_b[3]; // Udot, Vdot, Wdot (m/(s^2))
    double vel_b[3]; // U, V, W (m/s)
    double vel_e[3]; // Ndot, Edot, Ndot (m/s)
    double pos_e[3]; // N, E, D (m)
    double omega_b[3]; // P, Q, R (rad/s)
    // double euler_rates[3]; // phi_dot, theta_dot, psi_dot (rad/s)
    double quat_rate[4];
    double quat[4];
    double euler[3]; // phi, theta, psi (rad)
    double dcm_be[3][3]; // direct cosine matrix
} QuadState;

typedef struct
{
    double mass;
    double inertia[3]; // Ixx, Iyy, Izz

    double d; // distance from center to motor
    double r_d; // virtual yaw moment arm
    // double r_ld; // lift-to-drag ratio
    double thrust_tc; // thrust time constant

    double c_d[3]; // aerodynamic drag coefficient

    // Motor configuration:
    /*
     * 3   (N)   1
     *     (D)  (E)
     * 2         4
     */
    double thrust[4]; // thrusts of motors

    QuadState state;
    Sensors sensors;
} Quad;

void init_quad(Quad *quad, double mass, double inertia[3], double d, double r_d, double c_d[3], double thrust_tc, double throttle_hover, double init_yaw);
void init_quad_sensors(Quad *quad, double eph, double epv, double fix, double visible_sats, double lat_lon_noise_std_dev, double alt_noise_std_dev, double speed_noise_std_dev, double acc_noise_std_dev, double gyro_noise_std_dev, double mag_decl, double mag_incl, double mag_scale, double mag_noise_std_dev, double temperature);

void six_dof(double dt, Quad *quad, double forces[3], double moments[3]);

void forces_moments_aerodynamic_model(Quad *quad, double wind_vel_e[3], double forces[], double moments[]);
void forces_moments_thrust_model(double dt, double thrust_commands[4], Quad *quad, double forces[], double moments[]);
void forces_moments_gravity_model(Quad *quad, double forces[], double moments[]);

void get_att_quaternion(Quad *quad, double q[4]);

void update_quad(Quad *quad, double thrust_commands[4], double wind_vel_e[3], double dt);

#endif