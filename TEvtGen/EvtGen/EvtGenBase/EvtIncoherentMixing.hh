// $Id: EvtIncoherentMixing.hh,v 1.6 2009-11-25 13:41:59 mwhitehe Exp $
#ifndef EVTINCOHERENTMIXING_HH 
#define EVTINCOHERENTMIXING_HH 1

// Include files

#include "EvtGenBase/EvtParticle.hh"

/** @class EvtIncoherentMixing EvtIncoherentMixing.hh EvtGenBase/EvtIncoherentMixing.hh
 *  Class to contain the parameters of the incoherent B0 and B0S mixing
 *
 *  @author Patrick Robbe
 *  @date   2003-10-09
 */
class EvtIncoherentMixing {
public:
  /// Standard constructor
  EvtIncoherentMixing( );

  ~EvtIncoherentMixing( ) ; ///< Destructor

  // activate or desactivate the Bs mixing
  static void setB0Mixing() ;
  static void unsetB0Mixing() ;

  // activate or desactivate the B0 mixing
  static void setBsMixing() ;
  static void unsetBsMixing() ;

  // is mixing activated ? 
  static bool doB0Mixing() ;
  static bool doBsMixing() ;

  // set values for the mixing
  static void setdGammad( double value ) ;
  static void setdeltamd( double value ) ;
  static void setdGammas( double value ) ;
  static void setdeltams( double value ) ;

  // get parameters for mixing
  static double getdGammad() ;
  static double getdeltamd() ;
  static double getdGammas() ;
  static double getdeltams() ;

  // Functions to obtain time and type of the mixing ( 1 mix, 0 unmix )
  static void incoherentB0Mix( const EvtId id, double &t, int &mix ) ;
  static void incoherentBsMix( const EvtId id, double &t, int &mix ) ;

  // Functions to check if a B has mixed (comes from a B)
  static bool isB0Mixed( EvtParticle * ) ;
  static bool isBsMixed( EvtParticle * ) ;

  // Functions for CP models
  // returns for particle p the decay time t , the flavour of the tag
  // side (B0, anti-B0, B_s0, anti-B_s0). The flavour of p can flip
  // with probability probB
  static void OtherB( EvtParticle * p , double & t , EvtId & otherb ,
                      double probB ) ;
  static void OtherB( EvtParticle * p , double & t , EvtId & otherb ) ;

  static bool flipIsEnabled() ;
  static void enableFlip() ;
  static void disableFlip() ;

protected:

private:

  static bool _doB0Mixing ;
  static bool _doBsMixing ;
  static double _dGammad ;
  static double _deltamd ;
  static double _dGammas ;
  static double _deltams ;
  static bool _enableFlip ;
};
#endif // EVTGENBASE_EVTINCOHERENTMIXING_HH
