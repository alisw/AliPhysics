#ifndef ALIAODRECODECAYHF_H
#define ALIAODRECODECAYHF_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliAODRecoDecayHF
// base class for AOD reconstructed heavy-flavour decays
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//***********************************************************

#include <TRef.h>
#include <TList.h>
#include "AliAODTrack.h"
#include "AliAODRecoDecay.h"

class AliAODEvent;
class AliRDHFCuts;
class AliKFParticle;

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
  void SetPrimaryVtxRef(TObject *vtx) { fEventPrimaryVtx = vtx; }
  AliAODVertex* GetPrimaryVtxRef() const { return (AliAODVertex*)(fEventPrimaryVtx.GetObject()); }
  void SetOwnPrimaryVtx(const AliAODVertex *vtx) { UnsetOwnPrimaryVtx(); fOwnPrimaryVtx = new AliAODVertex(*vtx);}
  void CheckOwnPrimaryVtx() const 
    {if(!fOwnPrimaryVtx) printf("fOwnPrimaryVtx not set"); return;}
  AliAODVertex* GetOwnPrimaryVtx() const {return fOwnPrimaryVtx;}
  void GetOwnPrimaryVtx(Double_t vtx[3]) const 
    {CheckOwnPrimaryVtx();fOwnPrimaryVtx->GetPosition(vtx);}
  void UnsetOwnPrimaryVtx() {if(fOwnPrimaryVtx) {delete fOwnPrimaryVtx; fOwnPrimaryVtx=0;} return;}
  AliAODVertex* GetPrimaryVtx() const { return (GetOwnPrimaryVtx() ? GetOwnPrimaryVtx() : GetPrimaryVtxRef()); }
  AliAODVertex* RemoveDaughtersFromPrimaryVtx(AliAODEvent *aod);  
  void          RecalculateImpPars(AliAODVertex *vtxAODNew,AliAODEvent *aod);


  // kinematics & topology
  Double_t DecayLength2() const 
    { return AliAODRecoDecay::DecayLength2(GetPrimaryVtx());}
  Double_t DecayLength() const 
    { return AliAODRecoDecay::DecayLength(GetPrimaryVtx());}
  Double_t DecayLengthError() const 
    { return AliAODRecoDecay::DecayLengthError(GetPrimaryVtx());}
  Double_t NormalizedDecayLength() const 
    { return AliAODRecoDecay::NormalizedDecayLength(GetPrimaryVtx());}
  Double_t NormalizedDecayLength2() const 
    { return AliAODRecoDecay::NormalizedDecayLength2(GetPrimaryVtx());}
  Double_t DecayLengthXY() const 
    { return AliAODRecoDecay::DecayLengthXY(GetPrimaryVtx());}
  Double_t DecayLengthXYError() const 
    { return AliAODRecoDecay::DecayLengthXYError(GetPrimaryVtx());}
  Double_t NormalizedDecayLengthXY() const 
    { return AliAODRecoDecay::NormalizedDecayLengthXY(GetPrimaryVtx());}
  Double_t Ct(UInt_t pdg) const 
    { return AliAODRecoDecay::Ct(pdg,GetPrimaryVtx());}
  Double_t CosPointingAngle() const 
    { return AliAODRecoDecay::CosPointingAngle(GetPrimaryVtx());}
  Double_t CosPointingAngleXY() const 
    { return AliAODRecoDecay::CosPointingAngleXY(GetPrimaryVtx());}
  Double_t ImpParXY() const 
    { return AliAODRecoDecay::ImpParXY(GetPrimaryVtx());}
  Double_t QtProngFlightLine(Int_t ip) const 
    { return AliAODRecoDecay::QtProngFlightLine(ip,GetPrimaryVtx());}
  Double_t QlProngFlightLine(Int_t ip) const 
    { return AliAODRecoDecay::QlProngFlightLine(ip,GetPrimaryVtx());}

  // prongs
  Double_t Getd0errProng(Int_t ip) const {return fd0err[ip];}
  void     Setd0errProngs(Int_t nprongs,Double_t *d0);
  Double_t Normalizedd0Prong(Int_t ip) const 
    {return Getd0Prong(ip)/Getd0errProng(ip);}
  
  void SetProngIDs(Int_t nIDs,UShort_t *id);
  UShort_t GetProngID(Int_t ip) const 
    {if(fProngID) {return fProngID[ip];} else {return 9999;}}

  // ITS clustermap for daughters
  Bool_t DaughterHasPointOnITSLayer(Int_t dg,Int_t l) const;

  // check if it is like-sign
  Bool_t IsLikeSign() const;

  // list of cuts
  void SetListOfCutsRef(TObject *obj) {fListOfCuts=obj;}
  TList *GetListOfCuts() const {return (TList*)(fListOfCuts.GetObject());}
  AliRDHFCuts *GetCuts(const char* name) const;

  // vertexing KF:
  AliKFParticle *ApplyVertexingKF(Int_t *iprongs,Int_t nprongs,Int_t *pdgs,
				 Bool_t topoCostraint,Double_t bzkG,
				 Double_t *mass) const;
  
  // misalign
  void Misalign(TString misal="null");

  // selection map
  void    SetSelectionBit(Int_t i) {SETBIT(fSelectionMap,i); return;}
  Bool_t  HasSelectionBit(Int_t i) const {return TESTBIT(fSelectionMap,i);}
  ULong_t GetSelectionMap() const {return fSelectionMap;}

  Int_t   NumberOfFakeDaughters() const;

  Bool_t  HasBadDaughters() const; // TPC+ITS tracks not passing the StandardCuts2010 with loose DCA

 protected:

  AliAODVertex *fOwnPrimaryVtx; // primary vertex for this candidate
  TRef          fEventPrimaryVtx; // ref to primary vertex of the event
  TRef          fListOfCuts;  // ref to the list of analysis cuts
  Double_t     *fd0err;  //[fNProngs] error on prongs rphi impact param [cm]
  UShort_t     *fProngID;  //[fNProngs] track ID of daughters
  ULong_t       fSelectionMap; // used to store outcome of selection in AliAnalysisVertexingHF

  ClassDef(AliAODRecoDecayHF,5)  // base class for AOD reconstructed heavy-flavour decays
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

inline Bool_t AliAODRecoDecayHF::IsLikeSign() const
{
  // check if it is like-sign

  Int_t ndg=GetNDaughters();
  if(!ndg) {
    printf("Daughters not available\n");
    return kFALSE;
  }
  Int_t chargeDg0 = ((AliAODTrack*)GetDaughter(0))->Charge();

  for(Int_t i=1; i<ndg; i++) {
    if(chargeDg0!=((AliAODTrack*)GetDaughter(i))->Charge()) return kFALSE;
  }

  return kTRUE;
}

inline AliRDHFCuts *AliAODRecoDecayHF::GetCuts(const char* name) const
{ 
  // returns the analysis cuts

  TList *list = GetListOfCuts();
  if(!list) return 0;


  return (AliRDHFCuts*)list->FindObject(name);
}

inline Bool_t AliAODRecoDecayHF::DaughterHasPointOnITSLayer(Int_t dg,Int_t l) const
{
  // ITS clustermap for daughters

  if(l<0 || l>5) {
    printf("ERROR: layer has to be in the range 0-5\n");
    return kFALSE;
  }
  AliAODTrack *t = (AliAODTrack*)GetDaughter(dg);
  if(!t) return kFALSE;

  return TESTBIT(t->GetITSClusterMap(),l);
}

inline Int_t AliAODRecoDecayHF::NumberOfFakeDaughters() const 
{
  // Count number of daughters with negative label

  Int_t nfakes=0;
  for(Int_t i=0; i<GetNDaughters(); i++) {
    AliAODTrack *track=(AliAODTrack*)GetDaughter(i);

    if(track->Charge()==0) { // this is a two prong decay
      AliAODRecoDecay *rd=(AliAODRecoDecay*)GetDaughter(i);
      for(Int_t j=0; j<rd->GetNDaughters(); j++) {
	AliAODTrack *track2=(AliAODTrack*)GetDaughter(j);
	if(track2->GetLabel()<0) nfakes++;
      }
      continue;
    }

    if(track->GetLabel()<0) nfakes++;
  }
  return nfakes;
}

inline void AliAODRecoDecayHF::Setd0errProngs(Int_t nprongs,Double_t *d0err) 
{
  if(nprongs!=GetNProngs()) { 
    printf("Wrong number of momenta, must be nProngs");
    return;
  }
  if(!fd0) {
    fd0err = new Double32_t[nprongs];
  }
  for(Int_t i=0;i<nprongs;i++) {
    fd0err[i] = d0err[i]; 
  }

  return;
}

inline Bool_t AliAODRecoDecayHF::HasBadDaughters() const {

  // checks if at least one of ITS+TPC daughters didn't pass the StandardCuts2010 with loose DCA (ie doesn't have the BIT(4) set)

  for(Int_t iDau=0; iDau<GetNDaughters(); iDau++){
    AliAODTrack* at=(AliAODTrack*)GetDaughter(iDau);
    if(at->Charge()==0) continue;
    if(at->GetTPCNcls()==0) continue;
    if(!(at->TestFilterMask(BIT(4)))) return kTRUE;
  }
  return kFALSE;
}


#endif

