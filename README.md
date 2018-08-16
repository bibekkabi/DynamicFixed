# DynamicFixed
A Dynamic Fixed-point Arithmetic Library in C/C++ 

This is a fixed-point arithmetic library for 32-bit wordlength. 

It's purpose is to compute the data dynamic range during run-time and assign the appropriate integer wordlength for the fixed-point implementation.

We have an example of fixed-point implementation of Lanczos based tridiagonalization algorithm using the dynamic fixed-point library.

It is a part of the published work:
Pradhan, Tapan, Bibek Kabi, Ramanarayan Mohanty, and Aurobinda Routray. "Development of numerical linear algebra algorithms in dynamic fixed‐point format: a case study of Lanczos tridiagonalization." International Journal of Circuit Theory and Applications 44, no. 6 (2016): 1222-1262.

# Project Showcase
The file Fixed_dyn.cpp contains the arithmetic operations (addition, substraction, multiplication, division and square root) for 32-bit fixed-point format. 

# Requirements
  - The Compiler must be C++11 compliant.
  - The user should ideally be familar with the [Q number format](https://en.wikipedia.org/wiki/Q_(number_format)) for fixed-points.

# How to use
This is a header file #include "Fixed_dyn.h". Just include it. 
