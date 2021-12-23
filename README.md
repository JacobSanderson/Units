# Units
Abstraction of a Unit, includes convertions and basic operations. 

------ EXAMPLE : Free Fall (No air resistance) -------

from units_test import Unit, Constants

g = Constants.g \n
velocity = Unit(1, "km.s^-1") \n 
time = Unit(2, "min") \n
distance = - g*time**2*(1/2) + velocity*time  \n

print( distance )                    # Unit( 4.93680e+04 m )  
print( distance.to("mi").text() )    # 105.24040 mi

--------------------------------------------------------
Notes:
 - When having a unit to the zero e.g. km^0 in an operation, it'll give a float, ussually
 - Cannot multiply by const on the left i.e. 2*Unit() gives Error, but Unit()*2 is valid
 - Try to use parenthesis, i.e. Unit()* 1/2 gives Error, but Unit* (1/2) is valid
 - Every operation translates back to metric, u have to convert back after it 

Limitations:
 - Can't compute sqrt, sin, cos, exp, or Unit()**Unit(), cuz is not implemented
 - TODO : implement temperatures units (°C,°F,K,...)
 - TODO : add computations of errors
 - TODO : implement nparrays that have the same fucntionality

Instalation:
  Download the .py file and put it in ../Python38/Lib/site-packages, or just put it in
  the same folder your working, then you can import it no problem.
  I'm using python 3.8 but it should work on versions above. No dependencies (yet).
