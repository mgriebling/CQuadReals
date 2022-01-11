/*
 * src/qd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Defines constants used in quad-double package.
 */
#include "config.h"
#include "qd_real.h"
#include "inline.h"
#include <math.h>

/* Some useful constants. */
static const qd_real qd_2pi = { 6.283185307179586232e+00, 2.449293598294706414e-16, -5.989539619436679332e-33, 2.224908441726730563e-49 };
static const qd_real qd_pi = { 3.141592653589793116e+00, 1.224646799147353207e-16, -2.994769809718339666e-33, 1.112454220863365282e-49 };
static const qd_real qd_pi2 = { 1.570796326794896558e+00, 6.123233995736766036e-17, -1.497384904859169833e-33, 5.562271104316826408e-50 };
static const qd_real qd_pi4 = { 7.853981633974482790e-01, 3.061616997868383018e-17, -7.486924524295849165e-34, 2.781135552158413204e-50 };
static const qd_real qd_3pi4 = { 2.356194490192344837e+00, 9.1848509936051484375e-17, 3.9168984647504003225e-33, -2.5867981632704860386e-49 };
static const qd_real qd_e = { 2.718281828459045091e+00, 1.445646891729250158e-16, -2.127717108038176765e-33, 1.515630159841218954e-49 };
static const qd_real qd_log2 = { 6.931471805599452862e-01, 2.319046813846299558e-17, 5.707708438416212066e-34, -3.582432210601811423e-50 };
static const qd_real qd_log10c = { 2.302585092994045901e+00, -2.170756223382249351e-16, -9.984262454465776570e-33, -4.023357454450206379e-49 };
static const qd_real qd_nan = { NAN, NAN, NAN, NAN };
static const qd_real qd_inf = { INFINITY, INFINITY, -INFINITY, INFINITY };

static const double qd_eps = 1.21543267145725e-63; // = 2^-209
static const double qd_min_normalized = 1.6259745436952323e-260; // = 2^(-1022 + 3*53)
static const qd_real qd_max = { 1.79769313486231570815e+308, 9.97920154767359795037e+291, 5.53956966280111259858e+275, 3.07507889307840487279e+259 };
static const qd_real qd_safe_max = {1.7976931080746007281e+308,  9.97920154767359795037e+291,  5.53956966280111259858e+275, 3.07507889307840487279e+259 };
static const int qd_ndigits = 62;

