/*
 * src/dd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 */
#include "config.h"
#include "dd_real.h"
#include "inline.h"
#include <math.h>

static const dd_real dd_2pi = { 6.283185307179586232e+00, 2.449293598294706414e-16 };
static const dd_real dd_pi = { 3.141592653589793116e+00, 1.224646799147353207e-16 };
static const dd_real dd_pi2 = { 1.570796326794896558e+00, 6.123233995736766036e-17 };
static const dd_real dd_pi4 = { 7.853981633974482790e-01, 3.061616997868383018e-17 };
static const dd_real dd_3pi4 = { 2.356194490192344837e+00, 9.1848509936051484375e-17 };
static const dd_real dd_e = { 2.718281828459045091e+00, 1.445646891729250158e-16 };
static const dd_real dd_log2 = { 6.931471805599452862e-01, 2.319046813846299558e-17 };
static const dd_real dd_log10c = { 2.302585092994045901e+00, -2.170756223382249351e-16 };
static const dd_real dd_nan = { NAN, NAN };
static const dd_real dd_inf = { INFINITY, INFINITY };
static const dd_real dd_zero = { 0, 0 };

static const double dd_eps = 4.93038065763132e-32;  // 2^-104
static const double dd_min_normalized = 2.0041683600089728e-292;  // = 2^(-1022 + 53)
static const dd_real dd_max = { 1.79769313486231570815e+308, 9.97920154767359795037e+291 };
static const dd_real dd_safe_max = { 1.7976931080746007281e+308, 9.97920154767359795037e+291 };
static const int dd_ndigits = 31;


