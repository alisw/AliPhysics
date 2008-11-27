// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2005

//=============================================================================
// two-particle correlation analyzer
//=============================================================================

#ifndef ALIDANALCORREL_H
#define ALIDANALCORREL_H
#include "AliDAnal.h"
#include "AliDPair.h"
class TH1D;
class TH2D;
class AliDEvent;

//=============================================================================
class AliDAnalCorrel : public AliDAnal {
   
 public:
  AliDAnalCorrel(Char_t *nam="correl", 
	      Double_t emi=-1, Double_t ema=1, 
	      Int_t pid0=0, Int_t pid1=0);     // constructor
  virtual ~AliDAnalCorrel(){}                     // destructor
  // process one (tru) or two (mix) events
  void Process(Int_t tmr, AliDEvent *ev0, AliDEvent *ev1, Double_t phirot);

 protected:
  Int_t    fPid0;                       // particle species 0
  Int_t    fPid1;                       // particle species 1
  Double_t fMass0;                      // mass 0
  Double_t fMass1;                      // mass 1
  AliDPair    fPa;                         // pair buffer for calculations

  ClassDef(AliDAnalCorrel,1)
};
//=============================================================================
#endif
