#ifndef ALIAODRECODECAYHF_H
#define ALIAODRECODECAYHF_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliAODRecoDecayHF
// base class for AOD reconstructed heavy-flavour decays
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//***********************************************************

#include "AliAODRecoDecay.h"

class AliAODRecoDecayHF : public AliAODRecoDecay {

 public:

  AliAODRecoDecayHF();
  AliAODRecoDecayHF(AliAODVertex *vtx2,Int_t nprongs,Short_t charge,
		    Double_t *px,Double_t *py,Double_t *pz,
		    Double_t *d0,Double_t *d0err);
  AliAODRecoDecayHF(AliAODVertex *vtx2,Int_t nprongs,Short_t charge,
		    Double_t *d0,Double_t *d0err);
  AliAODRecoDecayHF(Double_t vtx1[3],Double_t vtx2[3],
		    Int_t nprongs,Short_t charge,
		    Double_t *px,Double_t *py,Double_t *pz,Double_t *d0);
  virtual ~AliAODRecoDecayHF();

  AliAODRecoDecayHF(const AliAODRecoDecayHF& source);
  AliAODRecoDecayHF& operator=(const AliAODRecoDecayHF& source); 
   

  // primary vertex
  void SetOwnPrimaryVtx(const AliAODVertex *vtx2) { UnsetOwnPrimaryVtx(); fOwnPrimaryVtx = new AliAODVertex(*vtx2);}
  void CheckOwnPrimaryVtx() const 
    {if(!fOwnPrimaryVtx) printf("fOwnPrimaryVtx not set"); return;}
  AliAODVertex* GetOwnPrimaryVtx() const {return fOwnPrimaryVtx;}
  void GetOwnPrimaryVtx(Double_t vtx[3]) const 
    {CheckOwnPrimaryVtx();fOwnPrimaryVtx->GetPosition(vtx);}
  void UnsetOwnPrimaryVtx() {if(fOwnPrimaryVtx) {delete fOwnPrimaryVtx; fOwnPrimaryVtx=0;} return;}

  // kinematics & topology
  Double_t DecayLength() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::DecayLength(fOwnPrimaryVtx);}
  Double_t DecayLengthError() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::DecayLengthError(fOwnPrimaryVtx);}
  Double_t NormalizedDecayLength() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::NormalizedDecayLength(fOwnPrimaryVtx);}
  Double_t DecayLengthXY() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::DecayLengthXY(fOwnPrimaryVtx);}
  Double_t DecayLengthXYError() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::DecayLengthXYError(fOwnPrimaryVtx);}
  Double_t NormalizedDecayLengthXY() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::NormalizedDecayLengthXY(fOwnPrimaryVtx);}
  Double_t Ct(UInt_t pdg) const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::Ct(pdg,fOwnPrimaryVtx);}
  Double_t CosPointingAngle() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::CosPointingAngle(fOwnPrimaryVtx);}
  Double_t CosPointingAngleXY() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::CosPointingAngleXY(fOwnPrimaryVtx);}
  Double_t ImpParXY() const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::ImpParXY(fOwnPrimaryVtx);}
  Double_t QtProngFlightLine(Int_t ip) const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::QtProngFlightLine(ip,fOwnPrimaryVtx);}
  Double_t QlProngFlightLine(Int_t ip) const 
    {CheckOwnPrimaryVtx();return AliAODRecoDecay::QlProngFlightLine(ip,fOwnPrimaryVtx);}

  // prongs
  Double_t Getd0errProng(Int_t ip) const {return fd0err[ip];}
  Double_t Normalizedd0Prong(Int_t ip) const 
    {return Getd0Prong(ip)/Getd0errProng(ip);}
  
  void SetProngIDs(Int_t nIDs,UShort_t *id);
  UShort_t GetProngID(Int_t ip) const 
    {if(fProngID) {return fProngID[ip];} else {return 9999;}}


 protected:

  AliAODVertex *fOwnPrimaryVtx; // primary vertex for this candidate
  Double_t     *fd0err;  //[fNProngs] error on prongs rphi impact param [cm]
  UShort_t     *fProngID;  //[fNProngs] track ID of daughters

  ClassDef(AliAODRecoDecayHF,2)  // base class for AOD reconstructed 
                                 // heavy-flavour decays
};

inline void AliAODRecoDecayHF::SetProngIDs(Int_t nIDs,UShort_t *id) 
{
  if(nIDs!=GetNProngs()) { 
    printf("Wrong number of IDs, must be nProngs\n");
    return;
  }
  if(fProngID) delete [] fProngID;
  fProngID = new UShort_t[nIDs];
  for(Int_t i=0;i<nIDs;i++) 
    fProngID[i] = id[i]; 
  return;
}

#endif
