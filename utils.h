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

typedef double (* vFunctionCall)(double[]);

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
void integer_rk4(double dt, vFunctionCall D_X[], double state[], double X[], int len);

uint64_t get_time_usec();

#endif