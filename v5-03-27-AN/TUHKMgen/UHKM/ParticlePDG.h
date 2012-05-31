//
//  Copyright   : The FASTMC and SPHMC Collaboration
//  Author      : Ionut Cristian Arsene 
//  Affiliation : Oslo University, Norway & Institute for Space Sciences, Bucharest, Romania
//  e-mail      : i.c.arsene@fys.uio.no
//  Date        : 2007/05/30
//
//  This class is using the particle and decay lists provided by the 
//  THERMINATOR (Computer Physics Communications 174 669 (2006)) and
//  SHARE (Computer Physics Communications 167 229 (2005)) collaborations.
//

#ifndef PARTICLEPDG_H
#define PARTICLEPDG_H

#include <Rtypes.h>
#include "DecayChannel.h"

const Int_t kMaxDecayChannels = 100;

class ParticlePDG {
 public:
  ParticlePDG();
  ParticlePDG(const Char_t * const name, Int_t pdg, Double_t mass, Double_t width);
  ~ParticlePDG();
  
  void AddChannel(DecayChannel* channel);
  void SetName(const Char_t * const name) {
    for(Int_t i=0; i<9; i++)
      if(*(name+i) != '\0') fName[i] = *(name+i);
      else break;
  }
  void SetPDG(Int_t value) {fPDG = value;}
  void SetMass(Double_t value) {fMass = value;}
  void SetWidth(Double_t value) {fWidth = value;}
  void SetSpin(Double_t value) {fSpin = value;}
  void SetIsospin(Double_t value) {fIsospin = value;}
  void SetIsospinZ(Double_t value) {fIsospinZ = value;}
  void SetLightQNumber(Double_t value) {fLightQuarkNumber = value;}
  void SetLightAQNumber(Double_t value) {fAntiLightQuarkNumber = value;}
  void SetStrangeQNumber(Double_t value) {fStrangeQuarkNumber = value;}
  void SetStrangeAQNumber(Double_t value) {fAntiStrangeQuarkNumber = value;}
  void SetCharmQNumber(Double_t value) {fCharmQuarkNumber = value;}
  void SetCharmAQNumber(Double_t value) {fAntiCharmQuarkNumber = value;}
  void SetStable(Bool_t value) {fStable = value;}
  
  const Char_t* GetName() const {return fName;}
  Int_t GetPDG() const {return fPDG;}
  Double_t GetMass() const {return fMass;}
  Double_t GetWidth() const {return fWidth;}
  Int_t GetNDecayChannels() const {return fNDecayChannels;}
  Double_t GetSpin() const {return fSpin;}
  Double_t GetIsospin() const {return fIsospin;}
  Double_t GetIsospinZ() const {return fIsospinZ;}
  Double_t GetLightQNumber() const {return fLightQuarkNumber;}
  Double_t GetLightAQNumber() const {return fAntiLightQuarkNumber;}
  Double_t GetStrangeQNumber() const {return fStrangeQuarkNumber;}
  Double_t GetStrangeAQNumber() const {return fAntiStrangeQuarkNumber;}
  Double_t GetCharmQNumber() const {return fCharmQuarkNumber;}
  Double_t GetCharmAQNumber() const {return fAntiCharmQuarkNumber;}
  Double_t GetBaryonNumber() const {return (fLightQuarkNumber     + fStrangeQuarkNumber     + fCharmQuarkNumber -  
                                      fAntiLightQuarkNumber - fAntiStrangeQuarkNumber - fAntiCharmQuarkNumber)/3.;}
  Double_t GetStrangeness() const {return (fAntiStrangeQuarkNumber - fStrangeQuarkNumber);}
  Double_t GetCharmness() const {return (fCharmQuarkNumber - fAntiCharmQuarkNumber);}
  Double_t GetElectricCharge() const {return fIsospinZ + (GetBaryonNumber()+GetStrangeness()+GetCharmness())/2.;}
  Bool_t GetStableStatus() const {
    return fStable;
  }  

  Double_t GetFullBranching();
  DecayChannel* GetDecayChannel(Int_t i) const {
    if(0<=i && i<fNDecayChannels) return fDecayChannels[i];
    else return 0x0;
  }

 private:
  ParticlePDG(const ParticlePDG&);
  ParticlePDG& operator=(const ParticlePDG&);

  Char_t        fName[9];                      // particle name
  Int_t         fPDG;                          // PDG code
  Double_t      fMass;                         // mass
  Double_t      fWidth;                        // width
  Double_t      fSpin;                         // J
  Double_t      fIsospin;                      // I
  Double_t      fIsospinZ;                     // I3
  Double_t      fLightQuarkNumber;             // u, d quark number
  Double_t      fAntiLightQuarkNumber;         // u-, d- quark number
  Double_t      fStrangeQuarkNumber;           // s quark number
  Double_t      fAntiStrangeQuarkNumber;       // s- quark number
  Double_t      fCharmQuarkNumber;             // c quark number
  Double_t      fAntiCharmQuarkNumber;         // c- quark number
  Int_t         fNDecayChannels;               // number of decay channels
  Bool_t        fStable;                       // flag to turn on/off decay
  DecayChannel* fDecayChannels[kMaxDecayChannels];   // array of decay channels
};

#endif
