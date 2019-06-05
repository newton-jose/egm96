# egm96
The Earth Gravitational Model egm96 in Kotlin

This is a rewrite, in Kotlin, from egm96 FORTRAN original code called F477.f

https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/f477.f

I try to do a more cleaner program to make it ease to be developed by others.

So the program egm96.kts is a function with latitude and longitude as input and
the geoidal ondulation as return.  

All done in Kotlin :)

After download the egm96 and corrections coefficients files:

https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.z

https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/corrcoef.z

You can run the egm96.kts https://github.com/newton-jose/egm96/blob/master/egm96.kts with the command:

$ kotlinc  -script egm96.kts
