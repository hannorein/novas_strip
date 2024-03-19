#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "novas.h"
int main (void)
{

    double jd = 2460409.3076388887;
    
    on_surface geo_loc = {44.30, -76.47, 0.0, 10.0, 1010.0};
    observer obs = {1, geo_loc};
    
    double pos1[3], vel1[3];
    geo_posvel (jd, 0.0, 1, &obs, pos1, vel1);
    
    printf("%f %f %f\n", pos1[0]*AU_KM, pos1[1]*AU_KM, pos1[2]*AU_KM);

    return (0);
}
