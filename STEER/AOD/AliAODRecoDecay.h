#ifndef ALIAODRECODECAY_H
#define ALIAODRECODECAY_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliAODRecoDecay
// base class for AOD reconstructed decays
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//***********************************************************

#include <TMath.h>
#include <TRef.h>
#include <TClonesArray.h>
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"

class AliVVertex;

class AliAODRecoDecay : public AliVTrack {

 public:

  AliAODRecoDecay();
  AliAODRecoDecay(AliAODVertex *vtx2,Int_t nprongs,Short_t charge,
		  Double_t *px,Double_t *py,Double_t *pz,
		  Double_t *d0);
  AliAODRecoDecay(AliAODVertex *vtx2,Int_t nprongs,Short_t charge,
		  Double_t *d0);
  virtual ~AliAODRecoDecay();

  AliAODRecoDecay(const AliAODRecoDecay& source);
  AliAODRecoDecay& operator=(const AliAODRecoDecay& source); 
   

  // decay vertex
  Double_t GetSecVtxX() const {return GetSecondaryVtx()->GetX();}
  Double_t GetSecVtxY() const {return GetSecondaryVtx()->GetY();}
  Double_t GetSecVtxZ() const {return GetSecondaryVtx()->GetZ();}
  Double_t RadiusSecVtx() const;
  void     SetSecondaryVtx(AliAODVertex *vtx2) {fSecondaryVtx=vtx2;}
  AliAODVertex* GetSecondaryVtx() const { return (((AliAODVertex*)fSecondaryVtx.GetObject()) ? (AliAODVertex*)fSecondaryVtx.GetObject() : GetOwnSecondaryVtx()); }
  void     SetOwnSecondaryVtx(AliAODVertex *vtx2) {fOwnSecondaryVtx=vtx2;}
  AliAODVertex* GetOwnSecondaryVtx() const {return fOwnSecondaryVtx;}
  void     GetSecondaryVtx(Double_t vtx[3]) const;
  Double_t GetReducedChi2() const {return GetSecondaryVtx()->GetChi2perNDF();}
  Short_t  Charge() const {return fCharge;}
  Short_t  GetCharge() const {return fCharge;}
  void     SetCharge(Short_t charge=0) {fCharge=charge;}

  // Match to MC signal:
  // check if this candidate is matched to a MC signal
  // If no, return -1
  // If yes, return label (>=0) of the AliAODMCParticle
  // if ndgCk>0, checks also daughters PDGs
  Int_t    MatchToMC(Int_t pdgabs,TClonesArray *mcArray,Int_t ndgCk=0, const Int_t *pdgDg=0) const;

  // PID
  void      SetPID(Int_t nprongs,Double_t *pid);
  Double_t *GetPID() const { return fPID; }
  void      GetPIDProng(Int_t ip,Double_t *pid) const;
  virtual const Double_t *PID() const { return fPID; }

  // prong-to-prong DCAs
  void    SetDCAs(Int_t nDCA,Double_t *dca);
  void    SetDCA(Double_t dca); // 2 prong
  Double_t GetDCA(Int_t i=0) const {return fDCA[i];}

  //event and run number
  void SetEventRunNumbers(Int_t en,Int_t rn) 
    { fEventNumber=en; fRunNumber=rn; return; }
  Int_t GetEventNumber() const { return fEventNumber; }
  Int_t GetRunNumber() const { return fRunNumber; }

  // methods of AliVTrack
  virtual Int_t    GetID() const { return -1; }
  virtual UChar_t  GetITSClusterMap() const;
  virtual ULong_t  GetStatus() const;
  virtual Bool_t   GetXYZ(Double_t *p) const { return XvYvZv(p); }
  virtual Bool_t   GetCovarianceXYZPxPyPz(Double_t cv[21]) const;
  virtual Bool_t   PropagateToDCA(const AliVVertex* vtx,Double_t b,Double_t maxd,Double_t dz[2],Double_t covar[3]);

  // kinematics & topology
  Double_t Px() const; 
  Double_t Py() const;
  Double_t Pz() const;
  Double_t P2() const {return Px()*Px()+Py()*Py()+Pz()*Pz();}
  Double_t Pt2() const {return Px()*Px()+Py()*Py();}
  Double_t P() const {return TMath::Sqrt(P2());}
  Double_t Pt() const {return TMath::Sqrt(Pt2());}
  Double_t OneOverPt() const {return (Pt() ? 1./Pt() : 0.);}
  Bool_t   PxPyPz(Double_t p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE; }
  Double_t Phi() const {return TMath::Pi()+TMath::ATan2(-Py(),-Px());}
  Double_t Theta() const {return 0.5*TMath::Pi()-TMath::ATan(Pz()/(Pt()+1.e-13));}
  Double_t Eta() const {return 0.5*TMath::Log((P()+Pz())/(P()-Pz()+1.e-13));}
  Double_t Xv() const { return GetSecVtxX(); }
  Double_t Yv() const { return GetSecVtxY(); }
  Double_t Zv() const { return GetSecVtxZ(); }
  virtual Bool_t   XvYvZv(Double_t x[3]) const { x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return kTRUE; }
  Double_t E2(UInt_t pdg) const;
  Double_t E(UInt_t pdg) const {return TMath::Sqrt(E2(pdg));}
  Double_t Y(UInt_t pdg) const {Double_t e=E(pdg); return 0.5*TMath::Log((e+Pz())/(e-Pz()+1.e-13));}
  Double_t DecayLength2(Double_t point[3]) const;
  Double_t DecayLength(Double_t point[3]) const {return TMath::Sqrt(DecayLength2(point));}
  Double_t DecayLength2(AliAODVertex *vtx1) const
  {return GetSecondaryVtx()->Distance2ToVertex(vtx1);}
  Double_t DecayLength(AliAODVertex *vtx1) const
    {return TMath::Sqrt(DecayLength2(vtx1));}
  Double_t DecayLengthError2(AliAODVertex *vtx1) const
    {return GetSecondaryVtx()->Error2DistanceToVertex(vtx1);}
  Double_t DecayLengthError(AliAODVertex *vtx1) const
    {return TMath::Sqrt(DecayLengthError2(vtx1));}
  Double_t NormalizedDecayLength2(AliAODVertex *vtx1) const 
    {return DecayLength2(vtx1)/DecayLengthError2(vtx1);}
  Double_t NormalizedDecayLength(AliAODVertex *vtx1) const 
  {return TMath::Sqrt(NormalizedDecayLength2(vtx1));}
  Double_t DecayLengthXY(Double_t point[3]) const;
  Double_t DecayLengthXY(AliAODVertex *vtx1) const
    {return GetSecondaryVtx()->DistanceXYToVertex(vtx1);}
  Double_t DecayLengthXYError(AliAODVertex *vtx1) const
    {return GetSecondaryVtx()->ErrorDistanceXYToVertex(vtx1);}
  Double_t NormalizedDecayLengthXY(AliAODVertex *vtx1) const 
    {return DecayLengthXY(vtx1)/DecayLengthXYError(vtx1);}
  Double_t Ct(UInt_t pdg,Double_t point[3]) const;
  Double_t Ct(UInt_t pdg,AliAODVertex *vtx1) const;
  Double_t CosPointingAngle(Double_t point[3]) const;
  Double_t CosPointingAngle(AliAODVertex *vtx1) const;
  Double_t CosPointingAngleXY(Double_t point[3]) const;
  Double_t CosPointingAngleXY(AliAODVertex *vtx1) const;
  Double_t CosThetaStar(Int_t ip,UInt_t pdgvtx,UInt_t pdgprong0,UInt_t pdgprong1) const;
  Double_t InvMass2(Int_t npdg,UInt_t *pdg) const;
  Double_t InvMass(Int_t npdg,UInt_t *pdg) const {return TMath::Sqrt(InvMass2(npdg,pdg));}
  Double_t ImpParXY(Double_t point[3]) const;
  Double_t ImpParXY(AliAODVertex *vtx1) const;
  Bool_t   PassInvMassCut(Int_t pdgMom,Int_t npdgDg,UInt_t *pdgDg,Double_t cut) const;

  // prongs
  Int_t    GetNProngs() const {return fNProngs;}
  Int_t    GetNDaughters() const {return GetSecondaryVtx()->GetNDaughters();}
  TObject *GetDaughter(Int_t i) const {return (GetNDaughters()>i ? GetSecondaryVtx()->GetDaughter(i) : 0x0);}

  Short_t  ChargeProng(Int_t ip) const;
  Double_t Getd0Prong(Int_t ip) const {return fd0[ip];}
  Double_t Prodd0d0(Int_t ip1=0,Int_t ip2=0) const {return fd0[ip1]*fd0[ip2];} 
  void     SetPxPyPzProngs(Int_t nprongs,Double_t *px,Double_t *py,Double_t *pz);
  void     Setd0Prongs(Int_t nprongs,Double_t *d0);
  Double_t PxProng(Int_t ip) const {return fPx[ip];}
  Double_t PyProng(Int_t ip) const {return fPy[ip];}
  Double_t PzProng(Int_t ip) const {return fPz[ip];}
  Double_t PtProng(Int_t ip) const {return TMath::Sqrt(Pt2Prong(ip));}
  Double_t Pt2Prong(Int_t ip) const; 
  Double_t PProng(Int_t ip) const {return TMath::Sqrt(P2Prong(ip));}
  Double_t P2Prong(Int_t ip) const; 
  Double_t PhiProng(Int_t ip) const 
    {return TMath::ATan2(PyProng(ip),PxProng(ip));}
    Double_t ThetaProng(Int_t ip) const 
      {return 0.5*TMath::Pi()-TMath::ATan(PzProng(ip)/(PtProng(ip)+1.e-13));}
  Double_t EtaProng(Int_t ip) const 
    {return -TMath::Log(TMath::Tan(0.5*ThetaProng(ip)));}
  Double_t E2Prong(Int_t ip,UInt_t pdg) const;
  Double_t EProng(Int_t ip,UInt_t pdg) const {return TMath::Sqrt(E2Prong(ip,pdg));}
  Double_t YProng(Int_t ip,UInt_t pdg) const 
    {return 0.5*TMath::Log((EProng(ip,pdg)+PzProng(ip))/(EProng(ip,pdg)-PzProng(ip)+1.e-13));}
  Double_t Alpha() const;             // for Armenteros-Podolanski plot (V0's)
  Double_t QlProng(Int_t ip) const;
  Double_t QtProng(Int_t ip=0) const; // for Armenteros-Podolanski plot (V0's)
  Double_t QlProngFlightLine(Int_t ip,Double_t point[3]) const;
  Double_t QlProngFlightLine(Int_t ip,AliAODVertex *vtx1) const;
  Double_t QtProngFlightLine(Int_t ip,Double_t point[3]) const;
  Double_t QtProngFlightLine(Int_t ip,AliAODVertex *vtx1) const;
  Double_t InvMass2Prongs(Int_t ip1,Int_t ip2,UInt_t pdg1,UInt_t pdg2) const;
  Double_t ProngsRelAngle(Int_t ip1=0,Int_t ip2=1) const;

  // relate to other objects
  //Double_t DistanceToVertex(AliAODVertex *vtx) // distance to a AliAODVertex
  //Double_t DistanceToTrack(AliAODTrack *trk)   // distance to a AliAODTrack


  // print
  void    Print(Option_t* option = "") const;
  //void    PrintIndices() const {GetSecondaryVtx()->PrintIndices();}

  // dummy functions for inheritance from AliVParticle
  Double_t E() const 
    {printf("Dummy function; use AliAODRecoDecay::E(UInt_t pdg) instead"); return (Double_t)-999.;}
  Double_t Y() const 
    {printf("Dummy function; use AliAODRecoDecay::Y(UInt_t pdg) instead"); return (Double_t)-999.;}
  Double_t M() const 
    {printf("Dummy function"); return (Double_t)-999.;}
  Int_t GetLabel() const {return -1;}
  Int_t PdgCode()  const {return  0;}
  
 protected:

  Int_t    MatchToMC(Int_t pdgabs,TClonesArray *mcArray,Int_t dgLabels[10],Int_t ndg,Int_t ndgCk=0,const Int_t *pdgDg=0) const;
  Int_t    MatchToMC(Int_t pdgabs,TClonesArray *mcArray,Int_t dgLabels[10]) const { return MatchToMC(pdgabs,mcArray,dgLabels,GetNDaughters()); }

  TRef     fSecondaryVtx;  // decay vertex
  AliAODVertex *fOwnSecondaryVtx;  // temporary solution (to work outside AliAODEvent)
  Short_t  fCharge;  // charge, use this convention for prongs charges:
                     // if(charge== 0) even-index prongs are +
                     //                odd-index prongs are -
                     // if(charge==+1) even-index prongs are +
                     //                odd-index prongs are -
                     // if(charge==-1) even-index prongs are -
                     //                odd-index prongs are +

  // TEMPORARY, to be removed when we do analysis on AliAODEvent
  Int_t fNProngs;    // number of prongs
  Int_t fNDCA;       // number of dca's
  Int_t fNPID;       // number of PID probabilities
  Double32_t *fPx;   //[fNProngs] px of tracks at the vertex [GeV/c]
  Double32_t *fPy;   //[fNProngs] py of tracks at the vertex [GeV/c]
  Double32_t *fPz;   //[fNProngs] pz of tracks at the vertex [GeV/c]
  Double32_t *fd0;   //[fNProngs] rphi impact params w.r.t. Primary Vtx [cm]
  Double32_t *fDCA;  //[fNDCA] prong-to-prong DCA [cm]
                     // convention:fDCA[0]=p0p1,fDCA[1]=p0p2,fDCA[2]=p1p2,...
  Double32_t *fPID;  //[fNPID] combined pid
                     //  (combined detector response probabilities)
                            
  // TEMPORARY, to be removed when we do analysis on AliAODEvent
  Int_t fEventNumber;
  Int_t fRunNumber;
  // TO BE PUT IN SPECIAL MC CLASS
  //Bool_t   fSignal; // TRUE if signal, FALSE if background (for simulation)
  //Int_t  fTrkNum[2]; // numbers of the two decay tracks  
  //Int_t fPdg[2];  // PDG codes of the two tracks (for sim.)
  //Int_t fMum[2];  // PDG codes of the mothers    (for sim.)

  //

  ClassDef(AliAODRecoDecay,4)  // base class for AOD reconstructed decays
};


inline Short_t AliAODRecoDecay::ChargeProng(Int_t ip) const
{
  if(fCharge==0 || fCharge==+1) {
    if(ip%2==0) {
      return (Short_t)1;
    } else {
      return (Short_t)-1;
    }
  } else { // fCharge==-1
    if(ip%2==0) {
      return (Short_t)-1;
    } else {
      return (Short_t)1;
    }
  }
}

inline Double_t AliAODRecoDecay::RadiusSecVtx() const 
{ 
  return TMath::Sqrt(GetSecVtxX()*GetSecVtxX()+GetSecVtxY()*GetSecVtxY());
}

inline void AliAODRecoDecay::GetSecondaryVtx(Double_t vtx[3]) const 
{
  GetSecondaryVtx()->GetPosition(vtx);
  return;
}

inline Double_t AliAODRecoDecay::Px() const 
{
  Double_t px=0.; 
  for(Int_t i=0;i<GetNProngs();i++) px+=PxProng(i); 
  return px;
}

inline Double_t AliAODRecoDecay::Py() const 
{
  Double_t py=0.; 
  for(Int_t i=0;i<GetNProngs();i++) py+=PyProng(i); 
  return py;
}

inline Double_t AliAODRecoDecay::Pz() const 
{
  Double_t pz=0.; 
  for(Int_t i=0;i<GetNProngs();i++) pz+=PzProng(i); 
  return pz;
}

inline Double_t AliAODRecoDecay::Ct(UInt_t pdg,AliAODVertex *vtx1) const
{
  Double_t v[3];
  vtx1->GetPosition(v);
  return Ct(pdg,v);
}

inline Double_t AliAODRecoDecay::CosPointingAngle(AliAODVertex *vtx1) const
{
  Double_t v[3];
  vtx1->GetPosition(v);
  return CosPointingAngle(v);
}

inline Double_t AliAODRecoDecay::CosPointingAngleXY(AliAODVertex *vtx1) const
{
  Double_t v[3];
  vtx1->GetPosition(v);
  return CosPointingAngleXY(v);
}

inline Double_t AliAODRecoDecay::ImpParXY(AliAODVertex *vtx1) const
{
  Double_t v[3];
  vtx1->GetPosition(v);
  return ImpParXY(v);
}

inline Double_t AliAODRecoDecay::Pt2Prong(Int_t ip) const 
{
  return PxProng(ip)*PxProng(ip)+PyProng(ip)*PyProng(ip);
}

inline Double_t AliAODRecoDecay::P2Prong(Int_t ip) const 
{
  return Pt2Prong(ip)+PzProng(ip)*PzProng(ip);
}

inline Double_t AliAODRecoDecay::QlProngFlightLine(Int_t ip,AliAODVertex *vtx1) const
{
  Double_t v[3];
  vtx1->GetPosition(v);
  return QlProngFlightLine(ip,v);
}

inline Double_t AliAODRecoDecay::QtProngFlightLine(Int_t ip,AliAODVertex *vtx1) const
{
  Double_t v[3];
  vtx1->GetPosition(v);
  return QtProngFlightLine(ip,v);
}

inline void AliAODRecoDecay::Setd0Prongs(Int_t nprongs,Double_t *d0) 
{
  if(nprongs!=GetNProngs()) { 
    printf("Wrong number of momenta, must be nProngs");
    return;
  }
  if(!fd0) {
    fd0 = new Double32_t[nprongs];
  }
  for(Int_t i=0;i<nprongs;i++) {
    fd0[i] = d0[i]; 
  }

  return;
}

inline void AliAODRecoDecay::SetPxPyPzProngs(Int_t nprongs,Double_t *px,Double_t *py,Double_t *pz) 
{
  if(nprongs!=GetNProngs()) { 
    printf("Wrong number of momenta, must be nProngs");
    return;
  }
  if(!fPx) {
    fPx = new Double32_t[nprongs];
    fPy = new Double32_t[nprongs];
    fPz = new Double32_t[nprongs];
  }
  for(Int_t i=0;i<nprongs;i++) {
    fPx[i] = px[i]; 
    fPy[i] = py[i]; 
    fPz[i] = pz[i]; 
  }

  return;
}

inline void AliAODRecoDecay::SetDCAs(Int_t nDCA,Double_t *dca) 
{
  if(nDCA!=(GetNProngs()*(GetNProngs()-1)/2)) { 
    printf("Wrong number of DCAs, must be nProngs*(nProngs-1)/2");
    return;
  }
  if(fDCA) delete [] fDCA;
  fNDCA = nDCA;
  fDCA = new Double32_t[nDCA];
  for(Int_t i=0;i<nDCA;i++) 
    fDCA[i] = dca[i]; 
  return;
}

inline void AliAODRecoDecay::SetDCA(Double_t dca) 
{
  Double_t ddca[1]; ddca[0]=dca;
  SetDCAs(1,ddca);
  return;
}

inline void AliAODRecoDecay::SetPID(Int_t nprongs,Double_t *pid) 
{
  if(nprongs!=GetNProngs()) {
    printf("Wrong number of prongs");
    return;
  }
  if(fPID) delete [] fPID;
  fNPID = nprongs*5;
  fPID = new Double32_t[nprongs*5];
  for(Int_t i=0;i<nprongs;i++) 
    for(Int_t j=0;j<5;j++)
      fPID[i*5+j] = pid[i*5+j]; 
  return;
}

inline void AliAODRecoDecay::GetPIDProng(Int_t ip,Double_t *pid) const
{ 
  for(Int_t j=0;j<5;j++)
    pid[j] = fPID[ip*5+j];
  return;
}



#endif


