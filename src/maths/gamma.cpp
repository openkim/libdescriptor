//from https://people.math.sc.edu/Burkardt/c_src/c_src.html

#include <cmath>
#include <stdexcept>
#include "gamma.hpp"

//Allan Macleod, Algorithm AS 245,
double lgamma ( double xvalue){
  double alr2pi = 0.918938533204673;
  double r1[9] = {
    -2.66685511495,
    -24.4387534237,
    -21.9698958928,
     11.1667541262,
     3.13060547623,
     0.607771387771,
     11.9400905721,
     31.4690115749,
     15.2346874070 };
  double r2[9] = {
    -78.3359299449,
    -142.046296688,
     137.519416416,
     78.6994924154,
     4.16438922228,
     47.0668766060,
     313.399215894,
     263.505074721,
     43.3400022514 };
  double r3[9] = {
    -2.12159572323E+05,
     2.30661510616E+05,
     2.74647644705E+04,
    -4.02621119975E+04,
    -2.29660729780E+03,
    -1.16328495004E+05,
    -1.46025937511E+05,
    -2.42357409629E+04,
    -5.70691009324E+02 };
  double r4[5] = {
     0.279195317918525,
     0.4917317610505968,
     0.0692910599291889,
     3.350343815022304,
     6.012459259764103 };
  double value, x, x1, x2, y;
  double xlge = 510000.0;
  double xlgst = 1.0E+30;

  x = xvalue;
  value = 0.0;
  //Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
  if ( x < 1.5 ){
    if ( x < 0.5 ){
      value = - log ( x );
      y = x + 1.0;
      //Test whether X < machine epsilon.
      if ( y == 1.0 ){
        return value;
      }
    } else {
      value = 0.0;
      y = x;
      x = ( x - 0.5 ) - 0.5;
    }

    value = value + x * ((((
        r1[4]   * y
      + r1[3] ) * y
      + r1[2] ) * y
      + r1[1] ) * y
      + r1[0] ) / ((((
                  y
      + r1[8] ) * y
      + r1[7] ) * y
      + r1[6] ) * y
      + r1[5] );

    return value;
  }
  if ( x < 4.0 ){
    y = ( x - 1.0 ) - 1.0;

    value = y * ((((
        r2[4]   * x
      + r2[3] ) * x
      + r2[2] ) * x
      + r2[1] ) * x
      + r2[0] ) / ((((
                  x
      + r2[8] ) * x
      + r2[7] ) * x
      + r2[6] ) * x
      + r2[5] );
  }   else if ( x < 12.0 ) {
    value = ((((
        r3[4]   * x
      + r3[3] ) * x
      + r3[2] ) * x
      + r3[1] ) * x
      + r3[0] ) / ((((
                  x
      + r3[8] ) * x
      + r3[7] ) * x
      + r3[6] ) * x
      + r3[5] );
  } else {
    y = log ( x );
    value = x * ( y - 1.0 ) - 0.5 * y + alr2pi;

    if ( x <= xlge ) {
      x1 = 1.0 / x;
      x2 = x1 * x1;

      value = value + x1 * ( (
             r4[2]   *
        x2 + r4[1] ) *
        x2 + r4[0] ) / ( (
        x2 + r4[4] ) *
        x2 + r4[3] );
    }
  }

  return value;
}

double gamma_inc (double p,  double x){
  double a, an, arg, b, dif, factor, g, gin, pn[6], rn, term, value;
  int i;
  double acu = 1.0E-08;
  double oflo = 1.0E+37;
  double uflo = 1.0E-37;

  if ( p <= 0.0 ){
    throw std::runtime_error("p <= 0.0");
    // return -999999999.0; //error
  }

  if ( x < 0.0 ){
    throw std::runtime_error("x < 0.0");
    // return -999999999.0; //error
  }

  if ( x == 0.0 ){
    value = 0.0;
    return value;
  }

  g = lgamma(p);

  arg = p * std::log (x) - x - g;

  if ( arg < std::log ( uflo ) ){
    throw std::runtime_error("arg < std::log ( uflo )");
    // return -999999999.0; //error
  }
  factor = std::exp(arg);
  if ( x <= 1.0 || x < p ){
    gin = 1.0;
    term = 1.0;
    rn = p;

    for ( ; ; ){
      rn = rn + 1.0;
      term = term * x / rn;
      gin = gin + term;

      if ( term <= acu ){
        break;
      }
    }

    value = gin * factor / p;
    return value;
  }

  a = 1.0 - p;
  b = a + x + 1.0;
  term = 0.0;

  pn[0] = 1.0;
  pn[1] = x;
  pn[2] = x + 1.0;
  pn[3] = x * b;

  gin = pn[2] / pn[3];

  for ( ; ; ){
    a = a + 1.0;
    b = b + 2.0;
    term = term + 1.0;
    an = a * term;
    for ( i = 0; i <= 1; i++ ){
      pn[i+4] = b * pn[i+2] - an * pn[i];
    }

    if ( pn[5] != 0.0 ){
      rn = pn[4] / pn[5];
      dif = fabs ( gin - rn );
      if ( dif <= acu ){
        if ( dif <= acu * rn ){
          value = 1.0 - factor * gin;
          break;
        }
      }
      gin = rn;
    }

    for ( i = 0; i < 4; i++ ){
      pn[i] = pn[i+2];
    }

    if ( oflo <= fabs ( pn[4] ) ){
      for ( i = 0; i < 4; i++ )
      {
        pn[i] = pn[i] / oflo;
      }
    }
  }

  return value;
}


// void gamma_log_values ( int *n_data, double *x, double *fx ){
//   double fx_vec[20] = {
//       0.1524063822430784E+01,
//       0.7966778177017837E+00,
//       0.3982338580692348E+00,
//       0.1520596783998375E+00,
//       0.0000000000000000E+00,
//      -0.4987244125983972E-01,
//      -0.8537409000331584E-01,
//      -0.1081748095078604E+00,
//      -0.1196129141723712E+00,
//      -0.1207822376352452E+00,
//      -0.1125917656967557E+00,
//      -0.9580769740706586E-01,
//      -0.7108387291437216E-01,
//      -0.3898427592308333E-01,
//      0.00000000000000000E+00,
//      0.69314718055994530E+00,
//      0.17917594692280550E+01,
//      0.12801827480081469E+02,
//      0.39339884187199494E+02,
//      0.71257038967168009E+02 };

//   double x_vec[20] = {
//       0.20E+00,
//       0.40E+00,
//       0.60E+00,
//       0.80E+00,
//       1.00E+00,
//       1.10E+00,
//       1.20E+00,
//       1.30E+00,
//       1.40E+00,
//       1.50E+00,
//       1.60E+00,
//       1.70E+00,
//       1.80E+00,
//       1.90E+00,
//       2.00E+00,
//       3.00E+00,
//       4.00E+00,
//      10.00E+00,
//      20.00E+00,
//      30.00E+00 };

//   if ( *n_data < 0 ){
//     *n_data = 0;
//   }

//   *n_data = *n_data + 1;

//   if ( 20 < *n_data ){
//     *n_data = 0;
//     *x = 0.0;
//     *fx = 0.0;
//   } else {
//     *x = x_vec[*n_data-1];
//     *fx = fx_vec[*n_data-1];
//   }

//   return;
// }


// double lanczos_lngamma ( double z, int *ier ){
//   double a[9] = {
//          0.9999999999995183,
//        676.5203681218835,
//     - 1259.139216722289,
//        771.3234287757674,
//      - 176.6150291498386,
//         12.50734324009056,
//        - 0.1385710331296526,
//          0.9934937113930748E-05,
//          0.1659470187408462E-06 };
//   int j;
//   double lnsqrt2pi = 0.9189385332046727;
//   double tmp;
//   double value;

//   if ( z <= 0.0 )
//   {
//     *ier = 1;
//     value = 0.0;
//     return value;
//   }

//   *ier = 0;

//   value = 0.0;
//   tmp = z + 7.0;
//   for ( j = 8; 1 <= j; j-- )
//   {
//     value = value + a[j] / tmp;
//     tmp = tmp - 1.0;
//   }

//   value = value + a[0];
//   value = log ( value ) + lnsqrt2pi - ( z + 6.5 )
//     + ( z - 0.5 ) * log ( z + 6.5 );

//   return value;
// }
