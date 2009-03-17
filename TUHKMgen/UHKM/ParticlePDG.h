/*
  Copyright   : The FASTMC and SPHMC Collaboration
  Author      : Ionut Cristian Arsene 
  Affiliation : Oslo University, Norway & Institute for Space Sciences, Bucharest, Romania
  e-mail      : i.c.arsene@fys.uio.no
  Date        : 2007/05/30

  This class is using the particle and decay lists provided by the 
  THERMINATOR (Computer Physics Communications 174 669 (2006)) and
  SHARE (Computer Physics Communications 167 229 (2005)) collaborations.
*/

#ifndef PARTICLE_PDG
#define PARTICLE_PDG

#include "Rtypes.h"

#ifndef DECAY_CHANNEL
#include "DecayChannel.h"
#endif

const Int_t kMaxDecayChannels = 20;

class ParticlePDG {
 public:
  ParticlePDG();
  ParticlePDG(Char_t* name, Int_t pdg, Double_t mass, Double_t width);
  ~ParticlePDG();
  
  void AddChannel(DecayChannel* channel);
  void SetName(Char_t* name) {
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
  
  Char_t* GetName() {return fName;}
  Int_t GetPDG() {return fPDG;}
  Double_t GetMass() {return fMass;}
  Double_t GetWidth() {return fWidth;}
  Int_t GetNDecayChannels() {return fNDecayChannels;}
  Double_t GetSpin() {return fSpin;}
  Double_t GetIsospin() {return fIsospin;}
  Double_t GetIsospinZ() {return fIsospinZ;}
  Double_t GetLightQNumber() {return fLightQuarkNumber;}
  Double_t GetLightAQNumber() {return fAntiLightQuarkNumber;}
  Double_t GetStrangeQNumber() {return fStrangeQuarkNumber;}
  Double_t GetStrangeAQNumber() {return fAntiStrangeQuarkNumber;}
  Double_t GetCharmQNumber() {return fCharmQuarkNumber;}
  Double_t GetCharmAQNumber() {return fAntiCharmQuarkNumber;}
  Double_t GetBaryonNumber() {return (fLightQuarkNumber     + fStrangeQuarkNumber     + fCharmQuarkNumber -  
                                      fAntiLightQuarkNumber - fAntiStrangeQuarkNumber - fAntiCharmQuarkNumber)/3.;}
  Double_t GetStrangeness() {return (fAntiStrangeQuarkNumber - fStrangeQuarkNumber);}
  Double_t GetCharmness() {return (fCharmQuarkNumber - fAntiCharmQuarkNumber);}
  Double_t GetElectricCharge() {return fIsospinZ + (GetBaryonNumber()+GetStrangeness()+GetCharmness())/2.;}
  Bool_t GetStableStatus() {
    return fStable;
  }  

  Double_t GetFullBranching();
  DecayChannel* GetDecayChannel(Int_t i) {
    if(0<=i && i<fNDecayChannels) 
      return fDecayChannels[i];
    else
      return 0x0;
  }

 private:
  ParticlePDG(const ParticlePDG&);
  ParticlePDG& operator=(const ParticlePDG&);

  Char_t        fName[9];
  Int_t         fPDG;
  Double_t      fMass;
  Double_t      fWidth;
  Double_t      fSpin;                         // J
  Double_t      fIsospin;                      // I
  Double_t      fIsospinZ;                     // I3
  Double_t      fLightQuarkNumber;             // u, d quark number
  Double_t      fAntiLightQuarkNumber;         // u-, d- quark number
  Double_t      fStrangeQuarkNumber;           // s quark number
  Double_t      fAntiStrangeQuarkNumber;       // s- quark number
  Double_t      fCharmQuarkNumber;             // c quark number
  Double_t      fAntiCharmQuarkNumber;         // c- quark number
  Int_t         fNDecayChannels;
  Bool_t        fStable;                       // flag to turn on/off decay
  DecayChannel* fDecayChannels[kMaxDecayChannels];   
};

#endif
