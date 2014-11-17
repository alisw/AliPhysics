// $Id: EvtBToDDalitzCPK.hh,v 1.2 2009-03-16 16:37:55 robbep Exp $
#ifndef EVTGENMODELS_EVTBTODDALITZCPK_HH 
#define EVTGENMODELS_EVTBTODDALITZCPK_HH 1

// Include files
#include "EvtGenBase/EvtDecayAmp.hh"

/** @class EvtBToDDalitzCPK EvtBToDDalitzCPK.hh EvtGenModels/EvtBToDDalitzCPK.hh
 *  Decay Model for B->DK, (adds the possibility to use D0->Ks pi pi to 
 *  find gamma with a Dalitz analysis
 *
 *  @author Patrick Robbe
 *  @date   2003-12-08
 */

class EvtBToDDalitzCPK : public EvtDecayAmp {
public:
  /// Standard constructor
  EvtBToDDalitzCPK( ) { }

  virtual ~EvtBToDDalitzCPK( ); ///< Destructor

  virtual std::string getName ( ) ;
  EvtDecayBase * clone ( ) ;
  
  void decay ( EvtParticle * p ) ;
  void init ( ) ;
  
  void initProbMax ( ) ;
  
protected:

private:
  int _flag ;

};
#endif // EVTGENMODELS_EVTBTODDALITZCPK_HH
