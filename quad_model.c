#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "quad_model.h"

// Private function prototypes
void six_dof_kinematic(double dt, QuadState *quad_state);
void update_sensors(Quad *quad);

// Public function definitions
void init_quad(Quad *quad, double mass, double inertia[3], double d, double r_d, double r_ld, double c_d[3], double thrust_tc, double throttle_hover)
{
    int i;
    int j;

    quad->mass = mass;
    quad->d = d;
    quad->r_d = r_d;
    quad->r_ld = r_ld;
    quad->thrust_tc = thrust_tc;
    
    for(i = 0 ; i < 4 ; i++)
        quad->thrust[i] = throttle_hover / 4;

    for(i = 0 ; i < 3 ; i++)
    {
        quad->inertia[i] = inertia[i];
        quad->c_d[i] = c_d[i];

        quad->state.acc_b[i] = 0.0f;
        quad->state.vel_b[i] = 0.0f;
        quad->state.vel_e[i] = 0.0f;
        quad->state.pos_e[i] = 0.0f;
        quad->state.omega_b[i] = 0.0f;
        quad->state.euler_rates[i] = 0.0f;
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
    double p_dot, q_dot, r_dot;
    QuadState *quad_state;

    quad_state = &(quad->state);

    p_dot = (moments[0] - quad_state->omega_b[1]*quad_state->omega_b[2]*(quad->inertia[2] - quad->inertia[1]))/quad->inertia[0];
    q_dot = (moments[1] - quad_state->omega_b[0]*quad_state->omega_b[2]*(quad->inertia[0] - quad->inertia[2]))/quad->inertia[1];
    r_dot = (moments[2] - quad_state->omega_b[0]*quad_state->omega_b[1]*(quad->inertia[1] - quad->inertia[0]))/quad->inertia[2];

    quad_state->acc_b[0] = forces[0]/quad->mass + quad_state->vel_b[1]*quad_state->omega_b[2] - quad_state->vel_b[2]*quad_state->omega_b[1];
    quad_state->acc_b[1] = forces[1]/quad->mass - quad_state->vel_b[0]*quad_state->omega_b[2] + quad_state->vel_b[2]*quad_state->omega_b[0];
    quad_state->acc_b[2] = forces[2]/quad->mass + quad_state->vel_b[0]*quad_state->omega_b[1] - quad_state->vel_b[1]*quad_state->omega_b[0];

    six_dof_kinematic(dt, quad_state);

    quad_state->vel_b[0] = integrate(quad_state->vel_b[0], quad_state->acc_b[0], dt);
    quad_state->vel_b[1] = integrate(quad_state->vel_b[1], quad_state->acc_b[1], dt);
    quad_state->vel_b[2] = integrate(quad_state->vel_b[2], quad_state->acc_b[2], dt);

    quad_state->omega_b[0] = integrate(quad_state->omega_b[0], p_dot, dt);
    quad_state->omega_b[1] = integrate(quad_state->omega_b[1], q_dot, dt);
    quad_state->omega_b[2] = integrate(quad_state->omega_b[2], r_dot, dt);
}

void six_dof_kinematic(double dt, QuadState *quad_state)
{
    int i, j;

    double phi, theta, psi;

    phi = quad_state->euler[0];
    theta = quad_state->euler[1];
    psi = quad_state->euler[2];

    quad_state->euler_rates[0] = quad_state->omega_b[0] + sin(phi)*tan(theta)*quad_state->omega_b[1] + cos(phi)*tan(theta)*quad_state->omega_b[2];
    quad_state->euler_rates[1] = cos(phi)*quad_state->omega_b[1] - sin(phi)*quad_state->omega_b[2];
    quad_state->euler_rates[2] = (sin(phi)/cos(theta))*quad_state->omega_b[1] + (cos(phi)/cos(theta))*quad_state->omega_b[2];

    body_to_earth_rotation(quad_state->dcm_be, quad_state->vel_b, quad_state->vel_e);

    quad_state->euler[0] = wrap_angle_pi(integrate(quad_state->euler[0], quad_state->euler_rates[0], dt));
    quad_state->euler[1] = wrap_angle_pi(integrate(quad_state->euler[1], quad_state->euler_rates[1], dt));
    quad_state->euler[2] = wrap_angle_pi(integrate(quad_state->euler[2], quad_state->euler_rates[2], dt));

    quad_state->pos_e[0] = integrate(quad_state->pos_e[0], quad_state->vel_e[0], dt);
    quad_state->pos_e[1] = integrate(quad_state->pos_e[1], quad_state->vel_e[1], dt);
    quad_state->pos_e[2] = integrate(quad_state->pos_e[2], quad_state->vel_e[2], dt);

    calc_dcm_be(quad_state->euler, quad_state->dcm_be);
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
    double thrust_dot;

    for(i = 0 ; i < 4 ; i++)
    {
        thrust_dot = (-quad->thrust[i] + thrust_commands[i])/quad->thrust_tc;
        quad->thrust[i] = integrate(quad->thrust[i], thrust_dot, dt);
    }

    forces[0] = 0.0f;
    forces[1] = 0.0f;
    forces[2] = -(quad->thrust[0]+quad->thrust[1]+quad->thrust[2]+quad->thrust[3]);

    moments[0] = quad->d * (quad->thrust[3] - quad->thrust[1]);
    moments[1] = quad->d * (quad->thrust[0] - quad->thrust[2]);
    moments[2] = quad->r_d * (-quad->thrust[0]+quad->thrust[1]-quad->thrust[2]+quad->thrust[3]) / quad->r_ld;
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

void get_att_quaternion(Quad *quad, double q[4])
{
    q[0] = 0.0f;
    q[1] = 0.0f;
    q[2] = 0.0f;
    q[3] = 0.0f;
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