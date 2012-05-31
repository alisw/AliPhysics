#ifndef ALIDECAYEREVTGEN_H
#define ALIDECAYEREVTGEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
///////////////////////////////////////////////////////////////////////////
//  Implementation of AliDecayer using EvtGen package. It inherits       //
//  from AliDecayer.                                                     //
//                                                                       //
//  Contact: Giuseppe.Bruno@ba.infn.it  &  Fiorella.Fionda@ba.infn.it    //
///////////////////////////////////////////////////////////////////////////


#include "AliDecayer.h"

class EvtGen;
class EvtRandomEngine;
class EvtStdHep;
class TLorentzVector;
class TClonesArray;

class AliDecayerEvtGen : public AliDecayer
{
 public:
  AliDecayerEvtGen();
  AliDecayerEvtGen(const AliDecayerEvtGen &decayer);
  virtual ~AliDecayerEvtGen();
  virtual void  Init();
  virtual void  Decay(Int_t ipart, TLorentzVector *p);
  virtual Int_t ImportParticles(TClonesArray *particles);
  virtual void    SetForceDecay(Decay_t decay) {fDecay=decay;}
  virtual void    SetForceDecay(Int_t decay){SetForceDecay((Decay_t) decay);}
  virtual void    ForceDecay();
  virtual Float_t GetPartialBranchingRatio(Int_t ipart);
  virtual Float_t GetLifetime(Int_t kf);
  Char_t*           GetDecayTablePath() {return fDecayTablePath;} 
  virtual void    ReadDecayTable();
  Bool_t SetDecayTablePath(Char_t *path);
   
  private:
  void  Copy(TObject &decayer) const;
  AliDecayerEvtGen &operator=(const AliDecayerEvtGen &decayer) 
      {decayer.Copy(*this);return(*this);}

  protected:
  EvtRandomEngine *fRandomEngine;  //!pointer to EvtRandomEngine to generate random number   
  EvtGen *fGenerator;              //!pointer to EvtGen class interface 
  EvtStdHep *fEvtstdhep;           //!pointer to EvtGen common block
  Char_t *fDecayTablePath;         //!pointer to decay table path
  Char_t *fParticleTablePath;      //!pointer to particle table path
  Decay_t     fDecay;  	           // Forced decay case		

  ClassDef(AliDecayerEvtGen,0)  //AliDecayer implementation using EvtGen 
};

//converts from mm/c to s
const Double_t kconv=0.001/2.999792458e8;

#endif

 
