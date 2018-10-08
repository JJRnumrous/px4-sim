#ifndef SENSORS_MODEL_H_
#define SENSORS_MODEL_H_

#include "utils.h"



// Functions
void init_gps(GpsSensor *gps, double eph, double epv, double fix, double visible_sats, double lat_lon_noise_std_dev, double alt_noise_std_dev, double speed_noise_std_dev);
void update_gps(GpsSensor *gps, double pos_e[3], double vel_e[3]);

void init_imu(ImuSensor *imu, double acc_noise_std_dev, double gyro_noise_std_dev);
void update_imu(ImuSensor *imu, double acc_b[3], double omega_b[3], double dcm_be[3][3]);

void init_mag(MagSensor *mag, double mag_decl, double mag_incl, double mag_scale, double mag_noise_std_dev);
void update_mag(MagSensor *mag, double dcm_be[3][3]);

void init_baro(BaroSensor *baro, double temperature);
void update_baro(BaroSensor *baro, double noisy_alt);

#endif