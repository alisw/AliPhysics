#ifndef ALIFITHits_H
#define ALIFITHits_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////
//  Manager and hits classes for set:FIT
//    Alla.Maevskaya@cern.ch    //
////////////////////////////////////////////////
 
#include "AliHit.h"
 
class AliFITHits : public AliHit {
public:

  AliFITHits();//Empty ctor
  AliFITHits(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliFITHits(){}//Empty virtual dtor

  AliFITHits& operator=(const AliFITHits&)  { return *this; }
  AliFITHits(const AliFITHits& o):AliHit(),
    fVolume(999),
    fPmt(999),
    fMCP(999),
    fParticle(0),
    fEtot(0),
    fTime(0),
    fCharge(999),
    fPx(999),
    fPy(999),
    fPz(999),
    fEloss(999),
    fTlength(999),
    fNphot(999)
    { ((AliFITHits &) o).Copy(*this);}
  

  Int_t Volume() const {return fVolume;}
  Int_t Pmt() const {return fPmt;}
  Int_t MCP() const {return fMCP;}
  Float_t Particle() const {return fParticle;} 
  Double_t Etot() const {return fEtot;}
  Float_t Time() const {return fTime;}
  Float_t Charge() const {return fCharge;};
  Float_t Px() const {return fPx;};
  Float_t Py() const {return fPy;};
  Float_t Pz() const {return fPz;};
  Float_t Eloss() const {return fEloss;};
  Float_t Tlength() const {return fTlength;};
  Int_t   Nphot() const {return fNphot;}

private:
  Int_t      fVolume;   //T0 arm mark
  Int_t      fPmt;      //PMT number on MCP 
  Int_t      fMCP;      // # MCP 
  Int_t      fParticle; //Primary particle ID
  Double_t   fEtot;     //Energy of primary particle at the entrance to radiator 
  Float_t    fTime;     //Primary particle TOF 
  Float_t    fCharge;
  Float_t    fPx;
  Float_t    fPy;
  Float_t    fPz;
  Float_t    fEloss;
  Float_t    fTlength;
  Int_t      fNphot;
   
   ClassDef(AliFITHits,1)  //Hits for detector T0
};



#endif//ALIT0hit_H
