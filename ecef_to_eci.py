# ecef_to_eci_q19.py
#
# Usage: python3 ecef_to_eci.py year month day hour minute second \
#                               ecef_x_km ecef_y_km ecef_z_km
#  Text explaining script usage
# Parameters:
#  arg1: description of argument 1
#  arg2: description of argument 2
#  ...
# Output:
#  ECI components
#
# Written by Connor Walsh
# Other contributors: Brad Denby
#

# import Python modules
import sys # argv
import math
import numpy as np


# "constants"
w = 7.292115 * 10**-5 #rad/sec

# helper functions

## function description
# def calc_something(param1, param2):
#   pass

# initialize script arguments
year = float('nan') # time in years
month = float('nan') # time in months
day = float('nan') # time in days
hour = float('nan') # time in hours
minute = float('nan') # time in minutes
second = float('nan') # time in seconds
eci_x_km = float('nan') # eci x component 
eci_y_km = float('nan') # eci y component 
eci_z_km = float('nan') # eci z component 

# parse script arguments
if len(sys.argv)==10:
    year = int(sys.argv[1]) # time in years
    month = int(sys.argv[2]) # time in months
    day = int(sys.argv[3]) # time in days
    hour = int(sys.argv[4]) # time in hours
    minute = int(sys.argv[5]) # time in minutes
    second = float(sys.argv[6]) # time in seconds
    ecef_x_km = float(sys.argv[7]) # eci x component 
    ecef_y_km = float(sys.argv[8]) # eci y component 
    ecef_z_km = float(sys.argv[9]) # eci z component 

else:
  print(\
   'Usage: python3 ecef_to_eci.py year month day hour minute second ecef_x_km ecef_y_km ecef_z_km\
  ')
  exit()

# write script below this line
if month <= 2:
    year -= 1
    month += 12

JD = (day - 32075) \
    + 1461 * (year + 4800 + (month - 14)//12)//4 \
    +367*(month-2-(month-14)//12 * 12)//12 \
    -3*((year+4900+(month-14)//12)//100)//4

#finding JD from midnight
A = year // 100
B = 2 - A + (A // 4)

JD_mnight = math.floor(365.25 * (year + 4716)) \
          + math.floor(30.6001 * (month + 1)) \
          + day + B - 1524.5

D_frac = (hour + minute / 60.0 + second / 3600.0) / 24.0

# Final Julian Date (including fractional day part)
JD_UT1 = JD_mnight + D_frac

T_UT1 = (JD_UT1 - 2451545.0)/36525

#GMST angle in seconds
theta_GMST = 67310.54841 + \
            (876600*60*60+8640184.812866)*T_UT1 \
            + 0.093104 * T_UT1**2 + -6.2*10**-6 * T_UT1**3 
theta_GMST = theta_GMST

theta_GMST_rad = (theta_GMST % 86400)*w + 2*math.pi  # Convert from seconds to radians

rECEF = [[ecef_x_km],
        [ecef_y_km],
        [ecef_z_km]
]

Rz = np.array([[math.cos(-theta_GMST_rad), -math.sin(-theta_GMST_rad), 0],
                [math.sin(-theta_GMST_rad), math.cos(-theta_GMST_rad), 0],
                [0,                                  0,                1]])

# Multiply matrices
rECI = np.linalg.inv(Rz) @ rECEF

eci_x_km = rECI[0][0]
eci_y_km = rECI[1][0] 
eci_z_km = rECI[2][0]

print(eci_x_km)
print(eci_y_km)
print(eci_z_km)