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

void calc_dcm_be(double euler[3], double dcm_be[3][3])
{
    double c_phi, c_theta, c_psi;
    double s_phi, s_theta, s_psi;

    c_phi = cos(euler[0]);
    c_theta = cos(euler[1]);
    c_psi = cos(euler[2]);

    s_phi = sin(euler[0]);
    s_theta = sin(euler[1]);
    s_psi = sin(euler[2]);

    dcm_be[0][0] = c_psi*c_theta;
    dcm_be[0][1] = s_psi*c_theta;
    dcm_be[0][2] = -s_theta;
    dcm_be[1][0] = c_psi*s_theta*s_phi - s_psi*c_phi;
    dcm_be[1][1] = s_psi*s_theta*s_phi + c_psi*c_phi;
    dcm_be[1][2] = c_theta*s_phi;
    dcm_be[2][0] = c_psi*s_theta*c_phi + s_psi*s_phi;
    dcm_be[2][1] = s_psi*s_theta*c_phi - c_psi*s_phi;
    dcm_be[2][2] = c_theta*c_phi;
}

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
    wrapped = wrap_angle_2pi(wrapped);

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
// Inputs:      Forces and moments, quad states => mass,inertia[3], vel_b[3], omega_b[3]
// Output:      Value for derivitive at specific states.
// **************************************************************************************************************************************************
// TODO: update quad_state->a_b 
double *RK4(double dt, vFunctionCall F1,vFunctionCall F2,vFunctionCall F3,vFunctionCall F4,vFunctionCall F5,vFunctionCall F6, double FM[6], Quad *quad,double *ret)
{
    QuadState temp_quad;
    QuadState *quad_state;
    quad_state = &(quad->state);
    temp_quad = (quad->state);


    double k1 = dt*F1(FM,temp_quad);
    double l1 = dt*F2(FM,temp_quad);
    double m1 = dt*F3(FM,temp_quad);
    double n1 = dt*F4(FM,temp_quad);
    double o1 = dt*F5(FM,temp_quad);
    double p1 = dt*F6(FM,temp_quad);
    state_adder(&temp_quad, k1,l1,m1,n1,o1,p1,2.0f);
    double k2 = dt*F1(FM,temp_quad);
    double l2 = dt*F2(FM,temp_quad);
    double m2 = dt*F3(FM,temp_quad);
    double n2 = dt*F4(FM,temp_quad);
    double o2 = dt*F5(FM,temp_quad);
    double p2 = dt*F6(FM,temp_quad);
    state_adder(&temp_quad, k2,l2,m2,n2,o2,p2,2.0f);
    double k3 = dt*F1(FM,temp_quad);
    double l3 = dt*F2(FM,temp_quad);
    double m3 = dt*F3(FM,temp_quad);
    double n3 = dt*F4(FM,temp_quad);
    double o3 = dt*F5(FM,temp_quad);
    double p3 = dt*F6(FM,temp_quad);
    state_adder(&temp_quad, k3,l3,m3,n3,o3,p3,1.0f);
    double k4 = dt*F1(FM,temp_quad);
    double l4 = dt*F2(FM,temp_quad);
    double m4 = dt*F3(FM,temp_quad);
    double n4 = dt*F4(FM,temp_quad);
    double o4 = dt*F5(FM,temp_quad);
    double p4 = dt*F6(FM,temp_quad);
    temp_quad = (quad->state);
    //integrals answer
    ret[0] += (k1 + 2*k2 + 2*k3 + k4)/6;
    ret[1] += (l1 + 2*l2 + 2*l3 + l4)/6;
    ret[2] += (m1 + 2*m2 + 2*m3 + m4)/6;
    ret[3] += (n1 + 2*n2 + 2*n3 + n4)/6;
    ret[4] += (o1 + 2*o2 + 2*o3 + o4)/6;
    ret[5] += (p1 + 2*p2 + 2*p3 + p4)/6;       
}

void state_adder(QuadState *quad, double k,double l,double m,double n,double o,double p,double devider)
{
    quad->vel_b[0] += k/devider;
    quad->vel_b[1] += l/devider;
    quad->vel_b[2] += m/devider;
    quad->omega_b[0] += n/devider;
    quad->omega_b[1] += o/devider;
    quad->omega_b[2] += p/devider;
}

uint64_t get_time_usec()
{
    static struct timeval _time_stamp;
    gettimeofday(&_time_stamp, NULL);
    return _time_stamp.tv_sec * 1000000 + _time_stamp.tv_usec;
}