#ifndef ALIDECAYERSLEVTGEN_H
#define ALIDECAYERSLEVTGEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
///////////////////////////////////////////////////////////////////////////
//  Implementation of AliDecayer using EvtGen package. It inherits       //
//  from AliDecayer.                                                     //
//                                                                       //
//  Origin: Michal.Broz@cern.ch  				         //
///////////////////////////////////////////////////////////////////////////

#include "AliDecayer.h"

class EvtGen;
class EvtStdlibRandomEngine;
class EvtAbsRadCorr;
class EvtStdHep;
class TLorentzVector;
class TClonesArray;

class AliDecayerSLEvtGen : public AliDecayer
{
 public:
  AliDecayerSLEvtGen();
  AliDecayerSLEvtGen(const AliDecayerSLEvtGen &decayer);
  virtual ~AliDecayerSLEvtGen();
  virtual void  Init();
  virtual void  Decay(Int_t ipart, TLorentzVector *p);
  virtual void  DecayPolarized(Int_t ipart, TLorentzVector *p, Int_t alpha);
  virtual Int_t ImportParticles(TClonesArray *particles);
  Char_t*       GetDecayTablePath() {return fDecayTablePath;} 
  virtual void  ReadDecayTable();
  Bool_t SetDecayTablePath(Char_t *path);
  
  virtual void    SetForceDecay(Decay_t decay) {}
  virtual void    SetForceDecay(Int_t decay){}
  virtual void    ForceDecay(){}
  virtual Float_t GetPartialBranchingRatio(Int_t ipart){return 0;}
  virtual Float_t GetLifetime(Int_t kf){return 0;}
   
  private:
  void  Copy(TObject &decayer) const;
  AliDecayerSLEvtGen &operator=(const AliDecayerSLEvtGen &decayer) 
      {decayer.Copy(*this);return(*this);}

  protected:
  EvtStdlibRandomEngine *fRandomEngine;  //!pointer to EvtRandomEngine to generate random number   
  EvtAbsRadCorr *fRadCorrEngine;   //!pointer to EvtGenCorrEngine
  EvtGen *fGenerator;              //!pointer to EvtGen class interface 
  EvtStdHep *fEvtstdhep;           //!pointer to EvtGen common block
  Char_t *fDecayTablePath;         //!pointer to decay table path
  Char_t *fParticleTablePath;      //!pointer to particle table path	

  ClassDef(AliDecayerSLEvtGen,0)  //AliDecayer implementation using EvtGen 
};

//converts from mm/c to s
const Double_t kconv=0.001/2.999792458e8;

#endif

 
