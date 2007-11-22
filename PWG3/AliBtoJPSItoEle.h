#ifndef AliBtoJPSItoEle_H
#define AliBtoJPSItoEle_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          Class AliBtoJPSItoEle
//                 Reconstructed B -> J\PSI+X --> e+ e- *X  class
//      
// Note: the two decay tracks are labelled: 0 (positive electron)
//                                          1 (negative electron)
//
//         Origin: G.E. Bruno    giuseppe.bruno@ba.infn.it                  
//          based on Class for charm golden channel (D0->Kpi) (thanks to Andrea Dainese!)
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TMath.h>

//----------------------------------------------------------------------------
//     Some constants (masses)
//
// particle masses
const Double_t kMJPsi = 3.096916;     // J/Psi mass
const Double_t kMe    = 0.000510998902;  // elec  mass




//-----------------------------------------------------------------------------
class AliBtoJPSItoEle : public TObject {
 public:
  //
  AliBtoJPSItoEle();
  AliBtoJPSItoEle(Int_t ev,Int_t trkNum[2],  
	     Double_t v1[3],Double_t v2[3],Double_t dca,
	     Double_t mom[6],Double_t d0[2]);
  virtual ~AliBtoJPSItoEle();
  AliBtoJPSItoEle(const AliBtoJPSItoEle& btoJpsi);

  Double_t Alpha() const { return (Ql(0)-Ql(1))/(Ql(0)+Ql(1)); }
  void     ApplyPID(TString pidScheme="TRDTPCparam");
  Double_t ChildrenRelAngle() const; 
  void     ComputeWgts();
  void     CorrectWgt4BR(Double_t factor);
  Double_t CosPointing() const;
  Double_t CosPointingXY() const;
  void     CosThetaStar(Double_t &ctsJPSI) const;
  Double_t Ct() const {return Length()*kMJPsi/P();}
  Double_t Energy() const { return TMath::Sqrt(P()*P()+kMJPsi*kMJPsi); }
  Double_t Eta() const;
  Double_t EtaChild(Int_t child) const;
  Int_t    EventNo() const {return TMath::Abs(fEvent);}
  Double_t GetDCA() const { return 10000.*fDCA; }
  Int_t    GetTrkNum(Int_t child) const { return fTrkNum[child]; }
  Double_t Getd0Child(Int_t child) const { return fd0[child]; }
  Int_t    GetPdgChild(Int_t child) const { return fPdg[child]; }
  Int_t    GetPdgMum(Int_t child) const {return fMum[child]; }
  Int_t    GetPdgGMum(Int_t child) const {return fGMum[child]; }
  void     GetPIDresponse(Double_t resp0[5],Double_t resp1[5]) const; 
  Double_t GetTagEl(Int_t child) const; 
  void     GetWgts(Double_t &WgtJPsi) const;
  void     GetPrimaryVtx(Double_t vtx[3]) const 
    { vtx[0]=fV1x; vtx[1]=fV1y; vtx[2]=fV1z; }
  void     GetSecondaryVtx(Double_t vtx[3]) const 
    { vtx[0]=fV2x; vtx[1]=fV2y; vtx[2]=fV2z; }

  Double_t ImpPar() const;
  void     InvMass(Double_t &mJPsi) const;
  Bool_t   IsSignal() const { if(fSignal) return kTRUE; return kFALSE; } 
  Bool_t   IsJpsiPrimary() const { if(fJpsiPrimary) return kTRUE; return kFALSE; } 
  Double_t Length() const
    { return TMath::Sqrt((fV1x-fV2x)*(fV1x-fV2x)
			 +(fV1y-fV2y)*(fV1y-fV2y)+(fV1z-fV2z)*(fV1z-fV2z)); }
  Double_t P()  const { return TMath::Sqrt(Pt()*Pt()+Pz()*Pz()); } 
  Double_t PChild(Int_t child)   const { return TMath::Sqrt(fPx[child]*fPx[child]+fPy[child]*fPy[child]+fPz[child]*fPz[child]); } 
  Double_t ProdImpParams() const { return fd0[0]*fd0[1]; } 
  Double_t Pt() const { return TMath::Sqrt(Px()*Px()+Py()*Py()); }
  Double_t PtChild(Int_t child)  const { return TMath::Sqrt(fPx[child]*fPx[child]+fPy[child]*fPy[child]); } 
  Double_t Px() const { return (fPx[0]+fPx[1]); }
  Double_t Py() const { return (fPy[0]+fPy[1]); }
  Double_t Pz() const { return (fPz[0]+fPz[1]); }
  Double_t Ql(Int_t child) const;
  Double_t Qt() const;
  Double_t Rapidity() const { return 0.5*TMath::Log((Energy()+Pz())/(Energy()-Pz()+1.e-13)); }
  Bool_t   Select(const Double_t* cuts,Int_t& okB) const;
  void     SetPrimaryVtx(Double_t vtx[3]) 
    { fV1x=vtx[0]; fV1y=vtx[1]; fV1z=vtx[2]; }
  void     SetSignal() { fSignal =  kTRUE; }
  void     SetJpsiPrimary() { fJpsiPrimary =  kTRUE; }
  void     SetTOFmasses(Double_t mass[2]) 
    { fTOFmass[0]=mass[0]; fTOFmass[1]=mass[1]; }
  void     SetPIDresponse(Double_t resp0[5],Double_t resp1[5]); 
  void     SetPdgCodes(Int_t pdg[2]) {fPdg[0]=pdg[0];fPdg[1]=pdg[1]; }
  void     SetMumPdgCodes(Int_t mum[2]) {fMum[0]=mum[0];fMum[1]=mum[1]; }
  void     SetGMumPdgCodes(Int_t gmum[2]) {fGMum[0]=gmum[0];fGMum[1]=gmum[1]; }
  Double_t TRDTPCCombinedPIDParametrization(Double_t p) const;
  //
 private:
  //
  Bool_t   fSignal; // TRUE if signal, FALSE if background (for simulation) 
                    // (background are both combinatorial and primary J/psi)
  Bool_t   fJpsiPrimary; // TRUE if the current candidate is a primary J/psi, FALSE otherway (for simulation) 
  Int_t    fEvent;  // number of the event this B comes from
                    // -1 if the B comes from ev. mixing

  Int_t fTrkNum[2]; // numbers of the two decay tracks

  Double_t fV1x; // X-position of the primary vertex of the event
  Double_t fV1y; // Y-position of the primary vertex of the event
  Double_t fV1z; // Z-position of the primary vertex of the event
  Double_t fV2x; // X-position of the reconstructed secondary vertex
  Double_t fV2y; // Y-position of the reconstructed secondary vertex
  Double_t fV2z; // Z-position of the reconstructed secondary vertex
  Double_t fDCA; // DCA of the two tracks

  Double_t fPx[2];  // X,Y,Z
  Double_t fPy[2];  // momenta of the two tracks
  Double_t fPz[2];  // at the reconstructed vertex  

  Double_t fd0[2];  //  impact parameters in the bending plane

  Int_t fPdg[2];  // PDG codes of the two tracks       (for sim.)
  Int_t fMum[2];  // PDG codes of the mothers          (for sim.)
  Int_t fGMum[2]; // PDG codes of the grand-mothers    (for sim.)

  Double_t fTagEl[2];  // probability to be tagged as electron 
  Double_t fTagPi[2];  // probability to be tagged as pion 
  Double_t fTagKa[2];  // probability to be tagged as kaon 
  Double_t fTagPr[2];  // probability to be tagged as proton 
  Double_t fTagNid[2]; // probability to be tagged as "non-identified" 

  Double_t fPIDrespEl[2]; // det. response to be electron
  Double_t fPIDrespMu[2]; // det. response to be muon
  Double_t fPIDrespPi[2]; // det. response to be pion
  Double_t fPIDrespKa[2]; // det. response to be kaon
  Double_t fPIDrespPr[2]; // det. response to be proton
  Double_t fTOFmass[2]; // mass estimated by the TOF (-1000. if track not reached TOF)
  Double_t fWgtJPsi; // weights for the pair

  ClassDef(AliBtoJPSItoEle,1)  // Reconstructed B candidate class
};

#endif







