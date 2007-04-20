#ifndef AliD0toKpi_H
#define AliD0toKpi_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          Class AliD0toKpi
//                 Reconstructed D0 -> K^- pi^+ class
//      
// Note: the two decay tracks are labelled: 0 (positive track)
//                                          1 (negative track)
//
//         Origin: A. Dainese    andrea.dainese@lnl.infn.it                  
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TMath.h>

//----------------------------------------------------------------------------
//     Some constants (masses + parameterized TOF PID)
//
// particle masses
const Double_t kMD0 = 1.8645;  // D0  mass
const Double_t kMK  = 0.49368; // K+  mass
const Double_t kMPi = 0.13957; // pi+ mass

//  --- TOF tagging probabilities --- 
//  central HIJING
//  B = 0.4 T
//  tracking errors in TPC included
//  With TRD
//
//  *** Pb-Pb dNch/dy=6000 *** 
//
//  PIONS
const Int_t kPiBinsPbPb = 10;
const Double_t kPiBinWidthPbPb = 0.250;
const Double_t kPiTagPiPbPb[kPiBinsPbPb] = {0.211421,0.652184,0.624421,0.614727,0.610777,0.628015,0.631520,0.630324,0.637551,0.575235};
const Double_t kPiTagNidPbPb[kPiBinsPbPb] = {0.788579,0.347816,0.375579,0.385273,0.389223,0.371985,0.368480,0.369676,0.362449,0.424765};
//  KAONS
const Int_t kKBinsPbPb = 10;
const Double_t kKBinWidthPbPb = 0.250;
const Double_t kKTagKPbPb[kKBinsPbPb] = {0.000000,0.101255,0.397662,0.467586,0.517008,0.555023,0.584185,0.519029,0.464117,0.247308};
const Double_t kKTagPiPbPb[kKBinsPbPb] = {0.102049,0.289930,0.101930,0.057771,0.040286,0.028567,0.053108,0.094369,0.066302,0.247308};
const Double_t kKTagNidPbPb[kKBinsPbPb] = {0.897951,0.608815,0.500408,0.474643,0.442705,0.416410,0.362707,0.386603,0.469580,0.505383};
//  PROTONS
const Int_t kPBinsPbPb = 9;
const Double_t kPBinWidthPbPb = 0.500;
const Double_t kPTagPPbPb[kPBinsPbPb] = {0.017940,0.350681,0.535286,0.583264,0.562935,0.560524,0.545992,0.598060,0.351245};
const Double_t kPTagPiPbPb[kPBinsPbPb] = {0.195955,0.094949,0.039962,0.026039,0.007556,0.016986,0.030333,0.000000,0.000000};
const Double_t kPTagNidPbPb[kPBinsPbPb] = {0.786105,0.554370,0.424751,0.390697,0.429508,0.422491,0.423675,0.401940,0.648755};
//
// pp PYTHIA 
//
// *** cuts for pp ***
//
//  PIONS
const Int_t kPiBinsPP = 10;
const Double_t kPiBinWidthPP = 0.250;
const Double_t kPiTagPiPP[kPiBinsPP] = {0.194528,0.447097,0.603364,0.646413,0.647125,0.669157,0.688139,0.682564,0.689910,0.665710};
const Double_t kPiTagNidPP[kPiBinsPP] = {0.805472,0.552903,0.396636,0.353587,0.352875,0.330843,0.311861,0.317436,0.310090,0.334290};
//  KAONS
const Int_t kKBinsPP = 10;
const Double_t kKBinWidthPP = 0.250;
const Double_t kKTagKPP[kKBinsPP] = {0.000000,0.173393,0.439690,0.519423,0.587025,0.605372,0.586021,0.650139,0.444444,0.299363};
const Double_t kKTagPiPP[kKBinsPP] = {0.000000,0.001495,0.000000,-0.000000,-0.000000,0.000000,0.032258,0.060572,0.101449,0.242038};
const Double_t kKTagNidPP[kKBinsPP] = {1.000000,0.825112,0.560310,0.480577,0.412975,0.394628,0.381720,0.289289,0.454106,0.458599};
//  PROTONS
const Int_t kPBinsPP = 9;
const Double_t kPBinWidthPP = 0.500;
const Double_t kPTagPPP[kPBinsPP] = {0.029404,0.438640,0.613710,0.665152,0.634961,0.657711,0.703704,0.685714,0.235294};
const Double_t kPTagPiPP[kPBinsPP] = {0.000000,0.000000,0.000000,-0.000000,0.000000,0.000000,-0.000000,0.014286,-0.000000};
const Double_t kPTagNidPP[kPBinsPP] = {0.970596,0.561360,0.386290,0.334848,0.365039,0.342289,0.296296,0.300000,0.764706};




//-----------------------------------------------------------------------------
class AliD0toKpi : public TObject {
 public:
  //
  AliD0toKpi();
  AliD0toKpi(Int_t ev,Int_t trkNum[2],  
	     Double_t v1[3],Double_t v2[3],Double_t dca,
	     Double_t mom[6],Double_t d0[2]);
  virtual ~AliD0toKpi();
  AliD0toKpi(const AliD0toKpi& d0toKpi);

  Double_t Alpha() const { return (Ql(0)-Ql(1))/(Ql(0)+Ql(1)); }
  void     ApplyPID(TString pidScheme="TOFparamPbPb");
  Double_t ChildrenRelAngle() const; 
  void     ComputeWgts();
  void     CorrectWgt4BR(Double_t factor);
  Double_t CosPointing() const;
  Double_t CosPointingXY() const;
  void     CosThetaStar(Double_t &ctsD0,Double_t &ctsD0bar) const;
  Double_t Ct() const {return Length()*kMD0/P();}
  Double_t Energy() const { return TMath::Sqrt(P()*P()+kMD0*kMD0); }
  Double_t Eta() const;
  Double_t EtaChild(Int_t child) const;
  Int_t    EventNo() const {return TMath::Abs(fEvent);}
  Double_t GetDCA() const { return 10000.*fDCA; }
  Int_t    GetTrkNum(Int_t child) const { return fTrkNum[child]; }
  Double_t Getd0Child(Int_t child) const { return fd0[child]; }
  Int_t    GetPdgChild(Int_t child) const { return fPdg[child]; }
  Int_t    GetPdgMum(Int_t child) const {return fMum[child]; }
  void     GetWgts(Double_t &WgtD0,Double_t &WgtD0bar,TString sample) const;
  void     GetPrimaryVtx(Double_t vtx[3]) const 
    { vtx[0]=fV1x; vtx[1]=fV1y; vtx[2]=fV1z; }
  void     GetSecondaryVtx(Double_t vtx[3]) const 
    { vtx[0]=fV2x; vtx[1]=fV2y; vtx[2]=fV2z; }

  Double_t ImpPar() const;
  void     InvMass(Double_t &mD0,Double_t &mD0bar) const;
  Bool_t   IsSignal() const { if(fSignal) return kTRUE; return kFALSE; } 
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
  Bool_t   Select(const Double_t* cuts,Int_t& okD0,Int_t& okD0bar) const;
  void     SetPrimaryVtx(Double_t vtx[3]) 
    { fV1x=vtx[0]; fV1y=vtx[1]; fV1z=vtx[2]; }
  void     SetSignal() { fSignal =  kTRUE; }
  void     SetTOFmasses(Double_t mass[2]) 
    { fTOFmass[0]=mass[0]; fTOFmass[1]=mass[1]; }
  void     SetPIDresponse(Double_t resp0[5],Double_t resp1[5]); 
  void     SetPdgCodes(Int_t pdg[2]) {fPdg[0]=pdg[0];fPdg[1]=pdg[1]; }
  void     SetMumPdgCodes(Int_t mum[2]) {fMum[0]=mum[0];fMum[1]=mum[1]; }
  Double_t LinearInterpolation(Double_t p,Int_t nBins,Double_t Bin,
			       const Double_t *values) const;
  //
 private:
  //
  Bool_t   fSignal; // TRUE if signal, FALSE if background (for simulation)
  Int_t    fEvent;  // number of the event this D0 comes from
                    // -1 if the D0 comes from ev. mixing

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

  Int_t fPdg[2];  // PDG codes of the two tracks (for sim.)
  Int_t fMum[2];  // PDG codes of the mothers    (for sim.)

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

  Double_t fWgtAD0,fWgtAD0bar; // weights for the 3 samples
  Double_t fWgtBD0,fWgtBD0bar; // weights for the 3 samples 
  Double_t fWgtCD0,fWgtCD0bar; // A: (K,Pi)+(K,?) B: (?,Pi) C: (?,?)
  Double_t fWgtDD0,fWgtDD0bar; // D: all other pairs

  ClassDef(AliD0toKpi,1)  // Reconstructed D0 candidate class
};

#endif








