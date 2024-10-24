# Terread (Prepare orography for CCAM)

Terread is used to create orography datatsets for the Conformal Cubic
Atmospheric Model (CCAM).

## Website

For documentation, see our website at

[https://research.csiro.au/ccam/]

## Dependencies

Terread requires the NetCDF C library.

## Building terread

To build terread with intel, gnu and cray fortran compiler use

```
make
make GFORTRAN=yes
make CRAY=yes
```

Debugging is also enabled with

```
make TEST=yes
```
