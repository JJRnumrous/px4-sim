#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "utils.h"
#include "quad_model.h"

double deg2rad(double deg)
{
    return deg * M_PI / 180.0f;
}

double rad2deg(double rad)
{
    return rad * 180.0f / M_PI;
}

double zero_mean_noise(double std_dev)
{
    return rand_gauss(0.0, std_dev);
}

// Gauss distribution - polar method
double rand_gauss(double mean, double std_dev)
{
    static int first = 1;
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if(first)
    {
        srand(time(NULL));
        first = 0;
    }

    if (call == 1)
    {
        call = !call;
        return (mean + std_dev * (double)X2);
    }

    do
    {
        U1 = -1 + ((double)rand() / RAND_MAX) * 2;
        U2 = -1 + ((double)rand() / RAND_MAX) * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mean + std_dev * (double)X1);
}

void quat_to_euler(double q[4], double euler[3])
{
    // roll
    euler[0] = wrap_angle_pi(atan2(2*(q[0]*q[1] + q[2]*q[3]), 1 - 2*(pow(q[1], 2) + pow(q[2], 2))));
     // pitch
    double sinp = 2*(q[0]*q[2] - q[3]*q[1]);
    if(abs(sinp) >= 1)
        euler[1] = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    else
        euler[1] = wrap_angle_pi(asin(sinp));
     // yaw
    euler[2] = wrap_angle_pi(atan2(2*(q[0]*q[3] + q[1]*q[2]), 1 - 2*(pow(q[2], 2) + pow(q[3], 2))));
}
 void euler_to_quat(double euler[3], double quat[4])
{
    quat[0] = cos(euler[0]/2)*cos(euler[1]/2)*cos(euler[2]/2) + sin(euler[0]/2)*sin(euler[1]/2)*sin(euler[2]/2);
    quat[1] = sin(euler[0]/2)*cos(euler[1]/2)*cos(euler[2]/2) - cos(euler[0]/2)*sin(euler[1]/2)*sin(euler[2]/2);
    quat[2] = cos(euler[0]/2)*sin(euler[1]/2)*cos(euler[2]/2) + sin(euler[0]/2)*cos(euler[1]/2)*sin(euler[2]/2);
    quat[3] = cos(euler[0]/2)*cos(euler[1]/2)*sin(euler[2]/2) + sin(euler[0]/2)*sin(euler[1]/2)*cos(euler[2]/2);
}
 void calc_dcm_be(double q[4], double dcm_be[3][3])
{
    dcm_be[0][0] = pow(q[0], 2) + pow(q[1], 2) - pow(q[2], 2) - pow(q[3], 2);
    dcm_be[0][1] = 2*(q[1]*q[2] + q[0]*q[3]);
    dcm_be[0][2] = 2*(q[1]*q[3] - q[0]*q[2]);
    dcm_be[1][0] = 2*(q[1]*q[2] - q[0]*q[3]);
    dcm_be[1][1] = pow(q[0], 2) - pow(q[1], 2) + pow(q[2], 2) - pow(q[3], 2);
    dcm_be[1][2] = 2*(q[2]*q[3] + q[0]*q[1]);
    dcm_be[2][0] = 2*(q[1]*q[3] + q[0]*q[2]);
    dcm_be[2][1] = 2*(q[2]*q[3] - q[0]*q[1]);
    dcm_be[2][2] = pow(q[0], 2) - pow(q[1], 2) - pow(q[2], 2) + pow(q[3], 2);
}


// TODO : not used?
void body_to_earth_rotation(double dcm_be[3][3], double rotate_from[3], double rotate_to[3])
{
    int i, j;

    for(i = 0 ; i < 3 ; i++)
    {   
        rotate_to[i] = 0;
        for(j = 0 ; j < 3 ; j++)
        {
            rotate_to[i] += (dcm_be[j][i] * rotate_from[j]);
        }
    }
}

void earth_to_body_rotation(double dcm_be[3][3], double rotate_from[3], double rotate_to[3])
{
    int i, j;

    for(i = 0 ; i < 3 ; i++)
    {   
        rotate_to[i] = 0;
        for(j = 0 ; j < 3 ; j++)
        {
            rotate_to[i] += (dcm_be[i][j] * rotate_from[j]);
        }
    }
}

double wrap_angle_2pi(double angle)
{
    double wrapped;

    if(angle >= 0 && angle <= 2*M_PI)
        return angle;

    // Wrap between 0 and 2pi.
    wrapped = angle - (floor(angle/(2*M_PI)))*2*M_PI;

    return wrapped;
}

double wrap_angle_pi(double angle)
{
    double wrapped;

    if(angle >= -M_PI && angle <= M_PI)
        return angle;

    // Wrap between 0 and 2pi.
    wrapped = wrap_angle_2pi(angle);

    // Make sure wrapped angle is between pi and -pi
    if(wrapped > M_PI)
        wrapped -= (2*M_PI);

    return wrapped;
}

double integrate(double sum, double val, double dt)
{
    return integrate_euler(sum, val, dt);
}

double integrate_euler(double sum, double val, double dt)
{
    return sum + val*dt;
}

// **************************************************************************************************************************************************
// Function:    Derivative Functions for Kinetics
// Inputs:      dt - time; D_X - Derivitive functions; states - in order of D_X; X - integrated states; len - length of states
// Output:      void
// **************************************************************************************************************************************************
void integer_rk4(double dt, vFunctionCall D_X[], double states[], double X[], int len)
{
    int i;
    double k1[len], k2[len], k3[len], k4[len], states_k2[len], states_k3[len], states_k4[len];

    for(i=0; i<len;i++)
    {
        k1[i] = dt*D_X[i](states);
        states_k2[i] = states[i] + k1[i]/2.0;        
    }
    for(i=0; i<len;i++)
    {
        k2[i] = dt*D_X[i](states_k2);
        states_k3[i] = states[i] + k2[i]/2.0;
    }
    for(i=0; i<len;i++)
    {
        k3[i] = dt*D_X[i](states_k3);
        states_k4[i] = states[i] + k3[i];
    }
    for(i=0; i<len;i++)
    {
        k4[i] = dt*D_X[i](states_k4);        
    }
    
    for(i=0; i<len;i++)
    {
        X[i] = states[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;        
    }      
}

uint64_t get_time_usec()
{
    static struct timeval _time_stamp;
    gettimeofday(&_time_stamp, NULL);
    return _time_stamp.tv_sec * 1000000 + _time_stamp.tv_usec;
}