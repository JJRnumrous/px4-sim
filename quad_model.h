#ifndef QUAD_MODEL_H_
#define QUAD_MODEL_H_


#include "sensors_model.h"
#include "utils.h"



void init_quad(Quad *quad, double mass, double inertia[3], double d, double r_d, double r_ld, double c_d[3], double thrust_tc, double throttle_hover);
void init_quad_sensors(Quad *quad, double eph, double epv, double fix, double visible_sats, double lat_lon_noise_std_dev, double alt_noise_std_dev, double speed_noise_std_dev, double acc_noise_std_dev, double gyro_noise_std_dev, double mag_decl, double mag_incl, double mag_scale, double mag_noise_std_dev, double temperature);

void six_dof(double dt, Quad *quad, double forces[3], double moments[3]);

double D_U(double [6], Quad *quad);
double D_V(double [6], Quad *quad);
double D_W(double [6], Quad *quad);
double D_P(double [6], Quad *quad);
double D_Q(double [6], Quad *quad);
double D_R(double [6], Quad *quad);

void forces_moments_aerodynamic_model(Quad *quad, double wind_vel_e[3], double forces[], double moments[]);
void forces_moments_thrust_model(double dt, double thrust_commands[4], Quad *quad, double forces[], double moments[]);
void forces_moments_gravity_model(Quad *quad, double forces[], double moments[]);

void get_att_quaternion(Quad *quad, double q[4]);

void update_quad(Quad *quad, double thrust_commands[4], double wind_vel_e[3], double dt);

#endif