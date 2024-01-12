# Aeroemiss (Aerosol emissions for CCAM)

Aeroemiss is used to create aerosol emissions based on CMIP5 and CMIP6 forcings for the conformal cubic grid employed by the Conformal Cubic Atmospheric Model
(CCAM).

## Website

For documentation, see our website at

[https://confluence.csiro.au/display/CCAM/CCAM]

## Dependencies

Aeroemiss requires the NetCDF C library.

## Building aeroemiss

To build aeroemiss with intel, gnu and cray fortran compiler use

```
make
make GFORTRAN=yes
make CRAY=yes
```

Debugging is also enabled with

```
make TEST=yes
```
