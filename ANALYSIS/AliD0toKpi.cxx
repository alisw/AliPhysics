/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//----------------------------------------------------------------------------
//               Implementation of the D0toKpi class
//      for pp and PbPb interactions
// Note: the two decay tracks are labelled: 0 (positive track)
//                                          1 (negative track)
//            Origin: A. Dainese    andrea.dainese@pd.infn.it            
//----------------------------------------------------------------------------

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TVector3.h>

#include "AliD0toKpi.h"

ClassImp(AliD0toKpi)

//----------------------------------------------------------------------------
AliD0toKpi::AliD0toKpi() {
  // Default constructor
  
  fSignal = kFALSE;

  fEvent = 0;

  fTrkNum[0] = 0;
  fTrkNum[1] = 0;

  fV1x = 0.;
  fV1y = 0.;
  fV1z = 0.;
  fV2x = 0.;
  fV2y = 0.;
  fV2z = 0.;
  fDCA = 0.;

  fPx[0] = 0.;
  fPy[0] = 0.;
  fPz[0] = 0.;
  fPx[1] = 0.;
  fPy[1] = 0.;
  fPz[1] = 0.;

  fd0[0] = 0.;
  fd0[1] = 0.;

  fPdg[0] = 0;
  fPdg[1] = 0;
  fMum[0] = 0;
  fMum[1] = 0;

  fTagPi[0] = 0.;
  fTagPi[1] = 0.;
  fTagKa[0] = 0.;
  fTagKa[1] = 0.;
  fTagNid[0] = 0.;
  fTagNid[1] = 0.;

  fWgtAD0=fWgtAD0bar=fWgtBD0=fWgtBD0bar=fWgtCD0=fWgtCD0bar=fWgtDD0=fWgtDD0bar=0;

}
//----------------------------------------------------------------------------
AliD0toKpi::AliD0toKpi(Int_t ev,Int_t trkNum[2],
		       Double_t v1[3],Double_t v2[3], 
		       Double_t dca,
		       Double_t mom[6],Double_t d0[2]) {
  // Constructor

  fSignal = kFALSE;

  fEvent = ev;
  fTrkNum[0] = trkNum[0];
  fTrkNum[1] = trkNum[1];

  fV1x = v1[0];
  fV1y = v1[1];
  fV1z = v1[2];
  fV2x = v2[0];
  fV2y = v2[1];
  fV2z = v2[2];
  fDCA = dca;

  fPx[0] = mom[0];
  fPy[0] = mom[1];
  fPz[0] = mom[2];
  fPx[1] = mom[3];
  fPy[1] = mom[4];
  fPz[1] = mom[5];

  fd0[0] = d0[0];
  fd0[1] = d0[1];

  fPdg[0] = 0;
  fPdg[1] = 0;
  fMum[0] = 0;
  fMum[1] = 0;

  fTagPi[0]  = 0.;
  fTagPi[1]  = 0.;
  fTagKa[0]  = 0.;
  fTagKa[1]  = 0.;
  fTagNid[0] = 0.;
  fTagNid[1] = 0.;

  fWgtAD0=fWgtAD0bar=fWgtBD0=fWgtBD0bar=fWgtCD0=fWgtCD0bar=fWgtDD0=fWgtDD0bar=0;
}
//----------------------------------------------------------------------------
AliD0toKpi::~AliD0toKpi() {}
//____________________________________________________________________________
AliD0toKpi::AliD0toKpi( const AliD0toKpi& d0toKpi):TObject(d0toKpi) {
  // dummy copy constructor
}
//----------------------------------------------------------------------------
void AliD0toKpi::ApplyPID(const Char_t * pidScheme) {
  // Applies particle identification
  const char *tofparampbpb = strstr(pidScheme,"TOFparamPbPb");
  const char *tofparampp   = strstr(pidScheme,"TOFparamPP");

  if((tofparampbpb || tofparampp) && fPdg[0]==0) {
    printf("AliD0toKpi::ApplyPID :\n Warning: TOF parameterized PID can be used only for simulation!\n"); 
    return;
  }

  if(tofparampbpb) {
    // tagging of the positive track
    if(TMath::Abs(fPdg[0])==211 || TMath::Abs(fPdg[0])==13 
       || TMath::Abs(fPdg[0])==11) { // pion,muon,electron
      fTagPi[0]  = LinearInterpolation(PChild(0),kPiBinsPbPb,kPiBinWidthPbPb,kPiTagPiPbPb);
      fTagNid[0] = 1.-fTagPi[0];
      fTagKa[0]   = 0.;
      fTagPr[0]   = 0.;
    } 
    if(TMath::Abs(fPdg[0])==321) { // kaon
      fTagKa[0]   = LinearInterpolation(PChild(0),kKBinsPbPb,kKBinWidthPbPb,kKTagKPbPb);
      fTagNid[0] = LinearInterpolation(PChild(0),kKBinsPbPb,kKBinWidthPbPb,kKTagNidPbPb);
      if((fTagNid[0]+fTagKa[0])>1.) fTagNid[0] = 1.-fTagKa[0];
      fTagPi[0] = 1.-fTagNid[0]-fTagKa[0];
      fTagPr[0] = 0.;
    } 
    if(TMath::Abs(fPdg[0])==2212) { // proton
      fTagPr[0]  = LinearInterpolation(PChild(0),kPBinsPbPb,kPBinWidthPbPb,kPTagPPbPb);
      fTagNid[0] = LinearInterpolation(PChild(0),kPBinsPbPb,kPBinWidthPbPb,kPTagNidPbPb);
      if((fTagNid[0]+fTagPr[0])>1.) fTagNid[0] = 1.-fTagPr[0];
      fTagPi[0] = 1.-fTagNid[0]-fTagPr[0];
      fTagKa[0]   = 0.;
    } 
    // tagging of the negative track
    if(TMath::Abs(fPdg[1])==211 || TMath::Abs(fPdg[1])==13 
       || TMath::Abs(fPdg[1])==11) { // pion,muon,electron
      fTagPi[1]  = LinearInterpolation(PChild(1),kPiBinsPbPb,kPiBinWidthPbPb,kPiTagPiPbPb);
      fTagNid[1] = 1.-fTagPi[1];
      fTagKa[1]   = 0.;
      fTagPr[1]   = 0.;
    } 
    if(TMath::Abs(fPdg[1])==321) { // kaon
      fTagKa[1]   = LinearInterpolation(PChild(1),kKBinsPbPb,kKBinWidthPbPb,kKTagKPbPb);
      fTagNid[1] = LinearInterpolation(PChild(1),kKBinsPbPb,kKBinWidthPbPb,kKTagNidPbPb);
      if((fTagNid[1]+fTagKa[1])>1.) fTagNid[1] = 1.-fTagKa[1];
      fTagPi[1] = 1.-fTagNid[1]-fTagKa[1];
      fTagPr[1] = 0.;
    } 
    if(TMath::Abs(fPdg[1])==2212) { // proton
      fTagPr[1]  = LinearInterpolation(PChild(1),kPBinsPbPb,kPBinWidthPbPb,kPTagPPbPb);
      fTagNid[1] = LinearInterpolation(PChild(1),kPBinsPbPb,kPBinWidthPbPb,kPTagNidPbPb);
      if((fTagNid[1]+fTagPr[1])>1.) fTagNid[1] = 1.-fTagPr[1];
      fTagPi[1] = 1.-fTagNid[1]-fTagPr[1];
      fTagKa[1]   = 0.;
    } 
  }


  if(tofparampp) {
    // tagging of the positive track
    if(TMath::Abs(fPdg[0])==211 || TMath::Abs(fPdg[0])==13 
       || TMath::Abs(fPdg[0])==11) { // pion,muon,electron
      fTagPi[0]  = LinearInterpolation(PChild(0),kPiBinsPP,kPiBinWidthPP,kPiTagPiPP);
      fTagNid[0] = 1.-fTagPi[0];
      fTagKa[0]   = 0.;
      fTagPr[0]   = 0.;
    } 
    if(TMath::Abs(fPdg[0])==321) { // kaon
      fTagKa[0]   = LinearInterpolation(PChild(0),kKBinsPP,kKBinWidthPP,kKTagKPP);
      fTagNid[0] = LinearInterpolation(PChild(0),kKBinsPP,kKBinWidthPP,kKTagNidPP);
      if((fTagNid[0]+fTagKa[0])>1.) fTagNid[0] = 1.-fTagKa[0];
      fTagPi[0] = 1.-fTagNid[0]-fTagKa[0];
      fTagPr[0] = 0.;
    } 
    if(TMath::Abs(fPdg[0])==2212) { // proton
      fTagPr[0]  = LinearInterpolation(PChild(0),kPBinsPP,kPBinWidthPP,kPTagPPP);
      fTagNid[0] = LinearInterpolation(PChild(0),kPBinsPP,kPBinWidthPP,kPTagNidPP);
      if((fTagNid[0]+fTagPr[0])>1.) fTagNid[0] = 1.-fTagPr[0];
      fTagPi[0] = 1.-fTagNid[0]-fTagPr[0];
      fTagKa[0]   = 0.;
    } 
    // tagging of the negative track
    if(TMath::Abs(fPdg[1])==211 || TMath::Abs(fPdg[1])==13 
       || TMath::Abs(fPdg[1])==11) { // pion,muon,electron
      fTagPi[1]  = LinearInterpolation(PChild(1),kPiBinsPP,kPiBinWidthPP,kPiTagPiPP);
      fTagNid[1] = 1.-fTagPi[1];
      fTagKa[1]   = 0.;
      fTagPr[1]   = 0.;
    } 
    if(TMath::Abs(fPdg[1])==321) { // kaon
      fTagKa[1]   = LinearInterpolation(PChild(1),kKBinsPP,kKBinWidthPP,kKTagKPP);
      fTagNid[1] = LinearInterpolation(PChild(1),kKBinsPP,kKBinWidthPP,kKTagNidPP);
      if((fTagNid[1]+fTagKa[1])>1.) fTagNid[1] = 1.-fTagKa[1];
      fTagPi[1] = 1.-fTagNid[1]-fTagKa[1];
      fTagPr[1] = 0.;
    } 
    if(TMath::Abs(fPdg[1])==2212) { // proton
      fTagPr[1]  = LinearInterpolation(PChild(1),kPBinsPP,kPBinWidthPP,kPTagPPP);
      fTagNid[1] = LinearInterpolation(PChild(1),kPBinsPP,kPBinWidthPP,kPTagNidPP);
      if((fTagNid[1]+fTagPr[1])>1.) fTagNid[1] = 1.-fTagPr[1];
      fTagPi[1] = 1.-fTagNid[1]-fTagPr[1];
      fTagKa[1]   = 0.;
    } 
  }

  return;
}
//----------------------------------------------------------------------------
Double_t AliD0toKpi::ChildrenRelAngle() const {
  // relative angle between K and pi

  TVector3 mom0(fPx[0],fPy[0],fPz[0]);
  TVector3 mom1(fPx[1],fPy[1],fPz[1]);

  Double_t angle = mom0.Angle(mom1);

  return angle; 
}
//----------------------------------------------------------------------------
void AliD0toKpi::ComputeWgts() {
  // calculate the weights for PID


  // assignement of the weights from PID
  fWgtAD0    = fTagKa[1]*(fTagPi[0]+fTagNid[0]);
  fWgtAD0bar = fTagKa[0]*(fTagPi[1]+fTagNid[1]);
  fWgtBD0    = fTagPi[0]*fTagNid[1];
  fWgtBD0bar = fTagPi[1]*fTagNid[0];
  fWgtCD0    = fTagNid[0]*fTagNid[1];
  fWgtCD0bar = fTagNid[0]*fTagNid[1];
  fWgtDD0    = 1.-fWgtAD0-fWgtBD0-fWgtCD0;
  fWgtDD0bar = 1.-fWgtAD0bar-fWgtBD0bar-fWgtCD0bar;

  /*
  cerr<<fWgtAD0<<"  "<<fWgtAD0bar<<endl;
  cerr<<fWgtBD0<<"  "<<fWgtBD0bar<<endl;
  cerr<<fWgtCD0<<"  "<<fWgtCD0bar<<endl;

  if(fWgtAD0<0.) cerr<<"AliD0toKpi::ComputeWgts()  Negative weight!!!\n";
  if(fWgtAD0bar<0.) cerr<<"AliD0toKpi::ComputeWgts()  Negative weight!!!\n";
  if(fWgtBD0<0.) cerr<<"AliD0toKpi::ComputeWgts()  Negative weight!!!\n";
  if(fWgtBD0bar<0.) cerr<<"AliD0toKpi::ComputeWgts()  Negative weight!!!\n";
  if(fWgtCD0<0.) cerr<<"AliD0toKpi::ComputeWgts()  Negative weight!!!\n";
  if(fWgtCD0bar<0.) cerr<<"AliD0toKpi::ComputeWgts()  Negative weight!!!\n";
  */

  return;
}
//----------------------------------------------------------------------------
void AliD0toKpi::CorrectWgt4BR(Double_t factor) {
  // correct weights of background from charm 

  fWgtAD0    *= factor;
  fWgtAD0bar *= factor;
  fWgtBD0    *= factor;
  fWgtBD0bar *= factor;
  fWgtCD0    *= factor;
  fWgtCD0bar *= factor;
  fWgtDD0    *= factor;
  fWgtDD0bar *= factor;

  return;
}
//----------------------------------------------------------------------------
Double_t AliD0toKpi::CosPointing() const {
  // cosine of pointing angle in space

  TVector3 mom(Px(),Py(),Pz());
  TVector3 flight(fV2x-fV1x,fV2y-fV1y,fV2z-fV1z);

  Double_t pta = mom.Angle(flight);

  return TMath::Cos(pta); 
}
//----------------------------------------------------------------------------
Double_t AliD0toKpi::CosPointingXY() const {
  // cosine of pointing angle in transverse plane

  TVector3 momXY(Px(),Py(),0.);
  TVector3 flightXY(fV2x-fV1x,fV2y-fV1y,0.);

  Double_t ptaXY = momXY.Angle(flightXY);

  return TMath::Cos(ptaXY); 
}
//----------------------------------------------------------------------------
void AliD0toKpi::CosThetaStar(Double_t &ctsD0,Double_t &ctsD0bar) const {
  // cosine of decay angle in the D0 rest frame

  Double_t pStar = TMath::Sqrt(TMath::Power(kMD0*kMD0-kMK*kMK-kMPi*kMPi,2.)-4.*kMK*kMK*kMPi*kMPi)/(2.*kMD0);

  Double_t beta = P()/Energy();
  Double_t gamma = Energy()/kMD0;

  ctsD0 = (Ql(1)/gamma-beta*TMath::Sqrt(pStar*pStar+kMK*kMK))/pStar;
  // if(ctsD0 > 1.)  { cerr<<"AliD0toKpi::CosThetaStar: > 1 "<<ctsD0<<"!\n"; }
  // if(ctsD0 < -1.) { cerr<<"AliD0toKpi::CosThetaStar: < -1 "<<ctsD0<<"!\n"; }

  ctsD0bar = (Ql(0)/gamma-beta*TMath::Sqrt(pStar*pStar+kMK*kMK))/pStar;
  // if(ctsD0bar > 1.)  { cerr<<"AliD0toKpi::CosThetaStar: > 1 "<<ctsD0bar<<"!\n"; }
  // if(ctsD0bar < -1.) { cerr<<"AliD0toKpi::CosThetaStar: < -1 "<<ctsD0bar<<"!\n";} 

  return;
}
//----------------------------------------------------------------------------
Double_t AliD0toKpi::Eta() const {
  // pseudorapidity of the D0

  Double_t theta = TMath::Pi()/2.-TMath::ATan2(Pz(),Pt());
  Double_t eta = -TMath::Log(TMath::Tan(theta/2.));
  return eta;
}
//----------------------------------------------------------------------------
Double_t AliD0toKpi::EtaChild(Int_t child) const {
  // pseudorapidity of the decay tracks

  Double_t theta = TMath::Pi()/2.-TMath::ATan2(fPz[child],PtChild(child));
  Double_t eta = -TMath::Log(TMath::Tan(theta/2.));
  return eta;
}
//----------------------------------------------------------------------------
void AliD0toKpi::GetWgts(Double_t &WgtD0,Double_t &WgtD0bar,TString sample) 
  const {
  // returns the weights for pid

  const char *sampleA = strstr(sample.Data(),"A");
  const char *sampleB = strstr(sample.Data(),"B");
  const char *sampleC = strstr(sample.Data(),"C");
  const char *sampleD = strstr(sample.Data(),"D");
  const char *sampleABCD = strstr(sample.Data(),"ABCD");
  const char *sampleABC = strstr(sample.Data(),"ABC");
  const char *sampleBC = strstr(sample.Data(),"BC");

  if(sampleA) { WgtD0 = fWgtAD0;  WgtD0bar = fWgtAD0bar; }
  if(sampleB) { WgtD0 = fWgtBD0;  WgtD0bar = fWgtBD0bar; }
  if(sampleC) { WgtD0 = fWgtCD0;  WgtD0bar = fWgtCD0bar; }
  if(sampleD) { WgtD0 = fWgtDD0;  WgtD0bar = fWgtDD0bar; }
  if(sampleABCD) { 
    WgtD0    = fWgtAD0+fWgtBD0+fWgtCD0+fWgtDD0; 
    WgtD0bar = fWgtAD0bar+fWgtBD0bar+fWgtCD0bar+fWgtDD0bar; 
  }
  if(sampleABC) { 
    WgtD0    = fWgtAD0+fWgtBD0+fWgtCD0; 
    WgtD0bar = fWgtAD0bar+fWgtBD0bar+fWgtCD0bar; 
  }
  if(sampleBC) { 
    WgtD0    = fWgtBD0+fWgtCD0; 
    WgtD0bar = fWgtBD0bar+fWgtCD0bar; 
  }


  if(fSignal) {
    if(fMum[0]==421)  WgtD0bar = 0.;
    if(fMum[0]==-421) WgtD0 = 0.; 
  }

  return;
}
//----------------------------------------------------------------------------
void AliD0toKpi::InvMass(Double_t &mD0,Double_t &mD0bar) const {
  // invariant mass as D0 and as D0bar

  Double_t energy[2];

  // D0 -> K- Pi+
  energy[1] = TMath::Sqrt(kMK*kMK+PChild(1)*PChild(1));
  energy[0] = TMath::Sqrt(kMPi*kMPi+PChild(0)*PChild(0));

  mD0 = TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-P()*P());
    

  // D0bar -> K+ Pi-
  energy[0] = TMath::Sqrt(kMK*kMK+PChild(0)*PChild(0));
  energy[1] = TMath::Sqrt(kMPi*kMPi+PChild(1)*PChild(1));

  mD0bar = TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-P()*P());
    
  return;

}
//----------------------------------------------------------------------------
Double_t AliD0toKpi::Ql(Int_t child) const {
  // longitudinal momentum of decay tracks w.r.t. to D0 momentum

  Double_t qL;
  TVector3 mom(fPx[child],fPy[child],fPz[child]);
  TVector3 momD(Px(),Py(),Pz());

  qL = mom.Dot(momD)/momD.Mag();

  return qL ;
}
//----------------------------------------------------------------------------
Double_t AliD0toKpi::Qt() const {
  // transverse momentum of decay tracks w.r.t. to D0 momentum  

  TVector3 mom0(fPx[0],fPy[0],fPz[0]);
  TVector3 momD(Px(),Py(),Pz());

  return mom0.Perp(momD);
}
//----------------------------------------------------------------------------
Bool_t AliD0toKpi::Select(const Double_t* cuts,Int_t& okD0,Int_t& okD0bar) 
  const {
//
// This function compares the D0 with a set of cuts:
//
// cuts[0] = inv. mass half width [GeV]   
// cuts[1] = dca [micron]
// cuts[2] = cosThetaStar 
// cuts[3] = pTK [GeV/c]
// cuts[4] = pTPi [GeV/c]
// cuts[5] = d0K [micron]   upper limit!
// cuts[6] = d0Pi [micron]  upper limit!
// cuts[7] = d0d0 [micron^2]
// cuts[8] = cosThetaPoint
//
// If the D0/D0bar doesn't pass the cuts it sets the weights to 0
// If neither D0 nor D0bar pass the cuts return kFALSE
//
  Double_t mD0,mD0bar,ctsD0,ctsD0bar;
  okD0=1; okD0bar=1;

  if(PtChild(1) < cuts[3] || PtChild(0) < cuts[4]) okD0 = 0;
  if(PtChild(0) < cuts[3] || PtChild(1) < cuts[4]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(TMath::Abs(Getd0Child(1)) > cuts[5] || 
     TMath::Abs(Getd0Child(0)) > cuts[6]) okD0 = 0;
  if(TMath::Abs(Getd0Child(0)) > cuts[6] ||
     TMath::Abs(Getd0Child(1)) > cuts[5]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(GetDCA() > cuts[1]) { okD0 = okD0bar = 0; return kFALSE; }

  InvMass(mD0,mD0bar);
  if(TMath::Abs(mD0-kMD0)    > cuts[0]) okD0 = 0;
  if(TMath::Abs(mD0bar-kMD0) > cuts[0]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  CosThetaStar(ctsD0,ctsD0bar);
  if(TMath::Abs(ctsD0)    > cuts[2]) okD0 = 0;
  if(TMath::Abs(ctsD0bar) > cuts[2]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(ProdImpParams() > cuts[7]) { okD0 = okD0bar = 0; return kFALSE; }

  if(CosPointing()   < cuts[8]) { okD0 = okD0bar = 0; return kFALSE; }

  return kTRUE;
}
//-----------------------------------------------------------------------------
void AliD0toKpi::SetPIDresponse(Double_t resp0[5],Double_t resp1[5]) {
  // Set combined PID detector response probabilities

  fPIDrespEl[0] = resp0[0];
  fPIDrespEl[1] = resp1[0];
  fPIDrespMu[0] = resp0[1];
  fPIDrespMu[1] = resp1[1];
  fPIDrespPi[0] = resp0[2];
  fPIDrespPi[1] = resp1[2];
  fPIDrespKa[0] = resp0[3];
  fPIDrespKa[1] = resp1[3];
  fPIDrespPr[0] = resp0[4];
  fPIDrespPr[1] = resp1[4];

  return;
} 
//-----------------------------------------------------------------------------
void AliD0toKpi::DrawPIDinTOF(const Char_t * pidScheme) const {
  // Draw parameterized PID probabilities in TOF

  const char *tofparampbpb = strstr(pidScheme,"TOFparamPbPb");
  const char *tofparampp = strstr(pidScheme,"TOFparamPP");

  TH2F* framePi = new TH2F("framePi","Tag probabilities for PIONS",2,0,2.5,2,0,1);
  framePi->SetXTitle("p [GeV/c]"); 
  framePi->SetStats(0);
  TH2F* frameK = new TH2F("frameK","Tag probabilities for KAONS",2,0,2.5,2,0,1);
  frameK->SetXTitle("p [GeV/c]");
  frameK->SetStats(0);
  TH2F* frameP = new TH2F("frameP","Tag probabilities for PROTONS",2,0,4.5,2,0,1);
  frameP->SetXTitle("p [GeV/c]");
  frameP->SetStats(0);

  TH1F* hPiPi = new TH1F("hPiPi","Tag probabilities for PIONS",kPiBinsPbPb,0,2.5);
  TH1F* hPiNid = new TH1F("hPiNid","Tag probabilities for PIONS",kPiBinsPbPb,0,2.5);

  TH1F* hKK = new TH1F("hKK","Tag probabilities for KAONS",kKBinsPbPb,0,2.5);
  TH1F* hKNid = new TH1F("hKNid","Tag probabilities for KAONS",kKBinsPbPb,0,2.5);
  TH1F* hKPi = new TH1F("hKPi","Tag probabilities for KAONS",kKBinsPbPb,0,2.5);

  TH1F* hPP = new TH1F("hPP","Tag probabilities for PROTONS",kPBinsPbPb,0,4.5);
  TH1F* hPNid = new TH1F("hPNid","Tag probabilities for PROTONS",kPBinsPbPb,0,4.5);
  TH1F* hPPi = new TH1F("hPPi","Tag probabilities for PROTONS",kPBinsPbPb,0,4.5);


  if(tofparampbpb) {

    for(Int_t i=1; i<=kPiBinsPbPb; i++) {
      hPiPi->SetBinContent(i,kPiTagPiPbPb[i-1]);
      hPiNid->SetBinContent(i,kPiTagPiPbPb[i-1]+kPiTagNidPbPb[i-1]);
      
      hKK->SetBinContent(i,kKTagKPbPb[i-1]);
      hKPi->SetBinContent(i,kKTagKPbPb[i-1]+kKTagPiPbPb[i-1]);
      hKNid->SetBinContent(i,kKTagKPbPb[i-1]+kKTagPiPbPb[i-1]+kKTagNidPbPb[i-1]);
    }
    for(Int_t i=1; i<=kPBinsPbPb; i++) {    
      hPP->SetBinContent(i,kPTagPPbPb[i-1]);
      hPPi->SetBinContent(i,kPTagPPbPb[i-1]+kPTagPiPbPb[i-1]);
      hPNid->SetBinContent(i,kPTagPPbPb[i-1]+kPTagPiPbPb[i-1]+kPTagNidPbPb[i-1]);
    }

  } else if(tofparampp) {

    for(Int_t i=1; i<=kPiBinsPP; i++) {
      hPiPi->SetBinContent(i,kPiTagPiPP[i-1]);
      hPiNid->SetBinContent(i,kPiTagPiPP[i-1]+kPiTagNidPP[i-1]);
      
      hKK->SetBinContent(i,kKTagKPP[i-1]);
      hKPi->SetBinContent(i,kKTagKPP[i-1]+kKTagPiPP[i-1]);
      hKNid->SetBinContent(i,kKTagKPP[i-1]+kKTagPiPP[i-1]+kKTagNidPP[i-1]);
    }
    for(Int_t i=1; i<=kPBinsPP; i++) {    
      hPP->SetBinContent(i,kPTagPPP[i-1]);
      hPPi->SetBinContent(i,kPTagPPP[i-1]+kPTagPiPP[i-1]);
      hPNid->SetBinContent(i,kPTagPPP[i-1]+kPTagPiPP[i-1]+kPTagNidPP[i-1]);
    }

  } 


  TCanvas* c = new TCanvas("c","Parameterized PID in TOF",0,0,1000,400);
  c->Divide(3,1);
  c->cd(1);
  framePi->Draw();
  hPiNid->SetFillColor(18); hPiNid->Draw("same");
  hPiPi->SetFillColor(4); hPiPi->Draw("same");
  TPaveLabel* pav1 = new TPaveLabel(1,.2,1.4,.3,"#pi");
  pav1->SetBorderSize(0);
  pav1->Draw("same");
  TPaveLabel* pav2 = new TPaveLabel(1,.8,1.8,.9,"non-id");
  pav2->SetBorderSize(0);
  pav2->Draw("same");

  c->cd(2);
  frameK->Draw();
  hKNid->SetFillColor(18); hKNid->Draw("same");
  hKPi->SetFillColor(4); hKPi->Draw("same");
  hKK->SetFillColor(7); hKK->Draw("same");
  TPaveLabel* pav3 = new TPaveLabel(1,.2,1.5,.3,"K");
  pav3->SetBorderSize(0);
  pav3->Draw("same");
  TPaveLabel* pav4 = new TPaveLabel(1,.8,1.8,.9,"non-id");
  pav4->SetBorderSize(0);
  pav4->Draw("same");
  TPaveLabel* pav5 = new TPaveLabel(.4,.5,.8,.6,"#pi");
  pav5->SetBorderSize(0);
  pav5->Draw("same");

  c->cd(3);
  frameP->Draw();
  hPNid->SetFillColor(18); hPNid->Draw("same");
  hPPi->SetFillColor(4); hPPi->Draw("same");
  hPP->SetFillColor(3); hPP->Draw("same");
  TPaveLabel* pav6 = new TPaveLabel(1,.2,1.5,.3,"p");
  pav6->SetBorderSize(0);
  pav6->Draw("same");
  TPaveLabel* pav7 = new TPaveLabel(1,.8,2.6,.9,"non-id");
  pav7->SetBorderSize(0);
  pav7->Draw("same");
  TPaveLabel* pav8 = new TPaveLabel(.2,.5,1,.6,"#pi");
  pav8->SetBorderSize(0);
  pav8->Draw("same");


  return;
}
//----------------------------------------------------------------------------
Double_t AliD0toKpi::LinearInterpolation(Double_t p,Int_t nBins,Double_t Bin,
					 const Double_t *values) const {
  // a linear interpolation method

  Double_t value=0; 
  Double_t slope;

  if(p<0.5*Bin) {
    value = values[0];
  } else if(p>=(nBins-0.5)*Bin) {
    slope = (2*values[nBins-1]-values[nBins-2]-values[nBins-3])/Bin/2;
    value = values[nBins-2]+slope*(p-Bin*(nBins-1.5));
  } else {
    for(Int_t i=0; i<nBins; i++) {
      if(p<(i+0.5)*Bin) {
	slope = (values[i]-values[i-1])/Bin;
	value = values[i-1]+slope*(p-Bin*(i-0.5));
	break;
      }
    }
  }

  if(value<0.) value=0.;
  if(value>1.) value=1.;

  return value;
}
//----------------------------------------------------------------------------






/*
//____________________________________________________________________________
void AliD0toKpi::SetPtWgts4pp() {
  // Correct pt distribution in order to reproduce MNR pt slope
  // (for pp generated with PYTHIA min. bias)

  if(TMath::Abs(fMum[0]) != 421 && TMath::Abs(fMum[1]) != 421 &&
     TMath::Abs(fMum[0]) != 411 && TMath::Abs(fMum[1]) != 411) return;

  Double_t ptWgt = 1.;
  ptWgt = 2.05-0.47*Pt()+0.02*Pt()*Pt();
  if(Pt() >= 5.) ptWgt = 0.56*TMath::Exp(-0.12*Pt());

  fWgtAD0    *= ptWgt;
  fWgtAD0bar *= ptWgt;
  fWgtBD0    *= ptWgt;
  fWgtBD0bar *= ptWgt;
  fWgtCD0    *= ptWgt;
  fWgtCD0bar *= ptWgt;
  fWgtDD0    *= ptWgt;
  fWgtDD0bar *= ptWgt;

  return;
}
//____________________________________________________________________________
*/









