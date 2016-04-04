// $Id: EvtIncoherentMixing.cpp,v 1.13 2009-11-27 09:09:41 mwhitehe Exp $
// Include files 


// local
#include "EvtGenBase/EvtIncoherentMixing.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtRandom.hh"

//-----------------------------------------------------------------------------
// Implementation file for class : EvtIncoherentMixing
//
// 2003-10-09 : Patrick Robbe
//-----------------------------------------------------------------------------


bool EvtIncoherentMixing::_doB0Mixing = false ;
bool EvtIncoherentMixing::_doBsMixing = false ;
bool EvtIncoherentMixing::_enableFlip = false ;
double EvtIncoherentMixing::_dGammad = 0. ;
double EvtIncoherentMixing::_deltamd = 0.502e12 ;
// dGamma_s corresponds to DeltaGamma / Gamma = 10 %
double EvtIncoherentMixing::_dGammas = 6.852e10 ;
double EvtIncoherentMixing::_deltams = 20.e12 ;

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
EvtIncoherentMixing::EvtIncoherentMixing(  ) {
  _doB0Mixing = false ;
  _doBsMixing = false ;
  _dGammad = 0. ;
  // dGammas corresponds to DeltaGamma / Gamma = 10 %
  _dGammas = 6.852e10 ;
  _deltamd = 0.502e12 ;
  _deltams = 20.e12 ;
  _enableFlip = false ;
}
//=============================================================================
EvtIncoherentMixing::~EvtIncoherentMixing( ) 
{
}
// ============================================================================
void EvtIncoherentMixing::incoherentB0Mix( const EvtId id, double &t , 
                                           int &mix )
{
  static EvtId B0  = EvtPDL::getId( "B0" ) ;
  static EvtId B0B = EvtPDL::getId( "anti-B0" ) ;
 
  if ( ( B0 != id ) && ( B0B != id ) ) {
    report(Severity::Error,"EvtGen") << "Bad configuration in incoherentB0Mix" 
                           << std::endl ;
    ::abort() ;
  }
  
  double x = getdeltamd() * EvtPDL::getctau( B0 ) / EvtConst::c ;

  double y = getdGammad() * ( EvtPDL::getctau( B0 ) / EvtConst::c ) / 2. ;

  double fac = 1. ; // No CP violation

  double mixprob = ( x*x + y*y ) / ( x*x + y*y + ( 1./fac ) * 
                                     ( 2. + x*x - y*y ) ) ;

  int mixsign ;
  
  // decide if state is mixed
  mixsign = ( mixprob > EvtRandom::Flat( 0. , 1. ) ) ? -1 : 1 ;

  double prob ;
  
  do {
    t = -log( EvtRandom::Flat() ) * EvtPDL::getctau( B0 ) / ( 1. - y ) ;
    prob = ( 1. + exp( -2. * y * t / EvtPDL::getctau( B0 ) ) +
      mixsign * 2. * exp( -y * t / EvtPDL::getctau( B0 ) ) * 
      cos( getdeltamd() * t / EvtConst::c ) ) / 2. ;
  } while ( prob < 2. * EvtRandom::Flat() ) ;
 
  mix = 0 ;
  if ( mixsign == -1 ) mix = 1 ;
  
  return ;  
}
// ============================================================================
void EvtIncoherentMixing::incoherentBsMix( const EvtId id, double &t , 
                                           int &mix )
{
  static EvtId BS  = EvtPDL::getId( "B_s0" ) ;
  static EvtId BSB = EvtPDL::getId( "anti-B_s0" ) ;
 
  if ( ( BS != id ) && ( BSB != id ) ) {
    report(Severity::Error,"EvtGen") << "Bad configuration in incoherentBsMix" 
                           << std::endl ;
    ::abort() ;
  }
  
  double x = getdeltams() * EvtPDL::getctau( BS ) / EvtConst::c ;

  double y = getdGammas() * ( EvtPDL::getctau( BS ) / EvtConst::c ) / 2. ;

  double fac = 1. ; // No CP violation

  double mixprob = ( x*x + y*y ) / ( x*x + y*y + ( 1./fac ) * 
                                     ( 2. + x*x - y*y ) ) ;

  int mixsign ;
  
  // decide if state is mixed
  mixsign = ( mixprob > EvtRandom::Flat( 0. , 1. ) ) ? -1 : 1 ;

  double prob ;
  
  do {
    t = -log( EvtRandom::Flat() ) * EvtPDL::getctau( BS ) / ( 1. - y ) ;
    prob = ( 1. + exp( -2. * y * t / EvtPDL::getctau( BS ) ) +
      mixsign * 2. * exp( -y * t / EvtPDL::getctau( BS ) ) * 
      cos( getdeltams() * t / EvtConst::c ) ) / 2. ;
  } while ( prob < 2. * EvtRandom::Flat() ) ;
 
  mix = 0 ;
  if ( mixsign == -1 ) mix = 1 ;
  
  return ;  
}

// ========================================================================
bool EvtIncoherentMixing::isBsMixed ( EvtParticle * p ) 
{ 
  if ( ! ( p->getParent() ) ) return false ;
  
  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  
  if ( ( p->getId() != BS0 ) && ( p->getId() != BSB ) ) return false ;
  
  if ( ( p->getParent()->getId() == BS0 ) ||
       ( p->getParent()->getId() == BSB ) ) return true ;
  
  return false ;
}

// ========================================================================
bool EvtIncoherentMixing::isB0Mixed ( EvtParticle * p ) 
{ 
  if ( ! ( p->getParent() ) ) return false ;
  
  static EvtId B0 =EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");
  
  if ( ( p->getId() != B0 ) && ( p->getId() != B0B ) ) return false ;
  
  if ( ( p->getParent()->getId() == B0 ) ||
       ( p->getParent()->getId() == B0B ) ) return true ;
  
  return false ;
}
//============================================================================
// Return the tag of the event (ie the anti-flavour of the produced 
// B meson). Flip the flavour of the event with probB probability
//============================================================================
void EvtIncoherentMixing::OtherB( EvtParticle * p ,
                                  double & t ,
                                  EvtId & otherb ,
                                  double probB ) 
{
  //if(p->getId() == B0 || p->getId() == B0B) 
  //added by liming Zhang
  enableFlip();
  if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
    p->getParent()->setLifetime() ;
    t = p->getParent()->getLifetime() ;
  }
  else {
    p->setLifetime() ;
    t = p->getLifetime() ;
  }

  if ( flipIsEnabled() ) {
    //std::cout << " liming << flipIsEnabled " << std::endl;
    // Flip the flavour of the particle with probability probB
    bool isFlipped = ( EvtRandom::Flat( 0. , 1. ) < probB ) ;
    
    if ( isFlipped ) {
      if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
        p->getParent()
          ->setId( EvtPDL::chargeConj( p->getParent()->getId() ) ) ;
        p->setId( EvtPDL::chargeConj( p->getId() ) ) ;
      }
      else {
        p->setId( EvtPDL::chargeConj( p->getId() ) ) ;
      }
    }
  }
  
  if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
    // if B has mixed, tag flavour is charge conjugate of parent of B-meson
    otherb = EvtPDL::chargeConj( p->getParent()->getId() ) ;
  }
  else {
    // else it is opposite flavour than this B hadron
    otherb = EvtPDL::chargeConj( p->getId() ) ;
  }

  return ;
}
//============================================================================
// Return the tag of the event (ie the anti-flavour of the produced 
// B meson). No flip
//============================================================================
void EvtIncoherentMixing::OtherB( EvtParticle * p ,
                                  double & t ,
                                  EvtId & otherb ) 
{
  if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
    p->getParent()->setLifetime() ;
    t = p->getParent()->getLifetime() ;
  }
  else {
    p->setLifetime() ;
    t = p->getLifetime() ;
  }
  
  if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
    // if B has mixed, tag flavour is charge conjugate of parent of B-meson
    otherb = EvtPDL::chargeConj( p->getParent()->getId() ) ;
  }
  else {
    // else it is opposite flavour than this B hadron
    otherb = EvtPDL::chargeConj( p->getId() ) ;
  }

  return ;
}


// activate or desactivate the Bs mixing
void EvtIncoherentMixing::setB0Mixing()   { _doB0Mixing = true ; }
void EvtIncoherentMixing::unsetB0Mixing() { _doB0Mixing = false ; } 

// activate or desactivate the B0 mixing
void EvtIncoherentMixing::setBsMixing()   { _doBsMixing = true ; } 
void EvtIncoherentMixing::unsetBsMixing() { _doBsMixing = false ; } 

// is mixing activated ? 
bool EvtIncoherentMixing::doB0Mixing()  { return _doB0Mixing ; }
bool EvtIncoherentMixing::doBsMixing()  { return _doBsMixing ; }

// set values for the mixing
void EvtIncoherentMixing::setdGammad( double value )  { _dGammad = value ; } 
void EvtIncoherentMixing::setdeltamd( double value )  { _deltamd = value ; } 
void EvtIncoherentMixing::setdGammas( double value )  { _dGammas = value ; } 
void EvtIncoherentMixing::setdeltams( double value )  { _deltams = value ; } 

// get parameters for mixing
double EvtIncoherentMixing::getdGammad() { return _dGammad ; } 
double EvtIncoherentMixing::getdeltamd() { return _deltamd ; }
double EvtIncoherentMixing::getdGammas() { return _dGammas ; } 
double EvtIncoherentMixing::getdeltams() { return _deltams ; }

bool EvtIncoherentMixing::flipIsEnabled() { return _enableFlip ; } 
void EvtIncoherentMixing::enableFlip() { _enableFlip = true ; } 
void EvtIncoherentMixing::disableFlip() { _enableFlip = false ; } 
