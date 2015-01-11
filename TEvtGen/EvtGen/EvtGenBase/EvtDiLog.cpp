// $Id: EvtDiLog.cpp,v 1.3 2009-03-16 15:52:42 robbep Exp $
// Include files

// local
#include "EvtGenBase/EvtDiLog.hh"

#include <cmath>

//-----------------------------------------------------------------------------
// Implementation file for class : EvtDiLog
//
// 2007-01-23 : Patrick Robbe
//-----------------------------------------------------------------------------

double EvtDiLog::DiLog( double x ) {

  double h , t , y , s , a , alfa , b0, b1, b2 ;
  if ( x == 1. ) h = PI6 ;
  else if ( x == -1. ) h = -PI12 ;
  else {
    t = -x ;
    if ( t <= -2. ) {
      y = -1./(1.+t) ;
      s = 1. ;
      a = -PI3 + HF * ( std::pow( log(-t) , 2 ) - 
                        std::pow( log( 1. + 1./t ) , 2 ) ) ;
    } else if ( t < -1. ) {
      y = -1. - t ;
      s = -1. ;
      a = log( -t ) ;
      a = -PI6 + a * ( a + log( 1. + 1./t ) ) ;
    } else if ( t <= -HF ) {
      y = - (1. + t ) / t ;
      s = 1. ;
      a = log( -t ) ;
      a = -PI6 + a * ( -HF * a + log( 1. + t ) ) ;
    } else if ( t < 0 ) {
      y = -t / ( 1. + t ) ;
      s = -1. ;
      a = HF * std::pow( log( 1. + t ) , 2 ) ;
    } else if ( t <= 1. ) {
      y = t ;
      s = 1. ;
      a = 0. ;
    } else {
      y = 1. / t ;
      s = -1. ;
      a = PI6 + HF * std::pow( log( t ) , 2 ) ;
    }
    
    h = y + y - 1. ;
    alfa = h + h ;
    b1 = 0. ;
    b2 = 0. ;
    for ( int i = 19 ; i >= 0 ; --i ) {
      b0 = C[ i ] + alfa * b1 - b2 ;
      b2 = b1 ;
      b1 = b0 ;
    }
    
    h = -(s * ( b0 -h * b2 ) + a ) ;
  }
  
  return h ;
}
