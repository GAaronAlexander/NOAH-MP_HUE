# Notes for NOAH-MP HUE

1. Currently assuming that there is 10 meter winds as input, which effects the line 2283 in module_sf_noahmp.lsm, which changes the wind speed down to 2 meters
2. Currently using the exact same emissivity as other soils in urban environments. 
3. Currently, we have not changed the thermal diffusivity and the heat capacity. This effects the thermal property lsm subroutine called THERMOPROP

Groundalbedo has an option to change the surface parameters of albedo, but has not been changed yet. 



## Things that still need to be addressed: 

FAREA in the model still needs to be changed around (just divide the two adjacent values by one another!)