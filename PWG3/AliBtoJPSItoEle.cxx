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
//               Implementation of the BtoJPSItoEle class
//                  for pp and PbPb interactions
// Note: the two decay tracks are labelled: 0 (positive electron)
//                                          1 (negative electron)
//            Origin: G.E. Bruno    giuseppe.bruno@ba.infn.it            
//  based on Class for charm golden channel (D0->Kpi)
//----------------------------------------------------------------------------

// #include <Riostream.h> // for debugging

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TVector3.h>
#include <TString.h>

#include "AliBtoJPSItoEle.h"

ClassImp(AliBtoJPSItoEle)

//----------------------------------------------------------------------------
AliBtoJPSItoEle::AliBtoJPSItoEle() {
  // Default constructor
  
  fSignal = kFALSE;
  fJpsiPrimary = kFALSE;

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
  fGMum[0] = 0;
  fGMum[1] = 0;

  fTagPi[0] = 0.;
  fTagPi[1] = 0.;
  fTagKa[0] = 0.;
  fTagKa[1] = 0.;
  fTagNid[0] = 0.;
  fTagNid[1] = 0.;

  fWgtJPsi=0;

}
//----------------------------------------------------------------------------
AliBtoJPSItoEle::AliBtoJPSItoEle(Int_t ev,Int_t trkNum[2],
		       Double_t v1[3],Double_t v2[3], 
		       Double_t dca,
		       Double_t mom[6],Double_t d0[2]) {
  // Constructor

  fSignal = kFALSE;
  fJpsiPrimary = kFALSE;

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
  fGMum[0] = 0;
  fGMum[1] = 0;

  fTagPi[0]  = 0.;
  fTagPi[1]  = 0.;
  fTagKa[0]  = 0.;
  fTagKa[1]  = 0.;
  fTagNid[0] = 0.;
  fTagNid[1] = 0.;

  fWgtJPsi=0;
}
//----------------------------------------------------------------------------
AliBtoJPSItoEle::~AliBtoJPSItoEle() {}
//____________________________________________________________________________
AliBtoJPSItoEle::AliBtoJPSItoEle( const AliBtoJPSItoEle& btoJpsi):TObject(btoJpsi) {
  // dummy copy constructor
}
//----------------------------------------------------------------------------
void AliBtoJPSItoEle::ApplyPID(TString pidScheme) {
  // Applies particle identification

  if(!pidScheme.CompareTo("TRDTPCparam")  && fPdg[0]==0) {
    printf("AliBtoJPSItoEle::ApplyPID :\n Warning: TRD-TPC parameterized PID can be used only for simulation!\n"); 
    return;
  }

  if(!pidScheme.CompareTo("TRDTPCparam")) {
    // tagging of the positive track
    if(TMath::Abs(fPdg[0])==11) { // electron
      fTagEl[0] = 0.81;
      fTagNid[0] = 1.-fTagEl[0];
    }
    else if(TMath::Abs(fPdg[0])==211) { // pion
      fTagEl[0]   = TRDTPCCombinedPIDParametrization(PChild(0));
      fTagNid[0] = 1.-fTagEl[0];
    } 
    else { // all the others 
      fTagEl[0]  = 0.;
      fTagNid[0] = 1.;
    } 
    // tagging of the negative track
    if(TMath::Abs(fPdg[1])==11) { // electron
      fTagEl[1] = 0.81;
      fTagNid[1] = 1.-fTagEl[1];
    }
    else if(TMath::Abs(fPdg[1])==211) { // pion
      fTagEl[1]   = TRDTPCCombinedPIDParametrization(PChild(1));
      fTagNid[1] = 1.-fTagEl[1];
    }
    else { // all the others
      fTagEl[1]  = 0.;
      fTagNid[1] = 1.;
    }
  }

  if(!pidScheme.CompareTo("ESDCombinedPID")) {
    fTagEl[0]=fPIDrespEl[0];
    fTagEl[1]=fPIDrespEl[1];
    fTagNid[0] = 1.-fTagEl[0];
    fTagNid[1] = 1.-fTagEl[1];
  }
  return;
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::ChildrenRelAngle() const {
  // relative angle between K and pi

  TVector3 mom0(fPx[0],fPy[0],fPz[0]);
  TVector3 mom1(fPx[1],fPy[1],fPz[1]);

  Double_t angle = mom0.Angle(mom1);

  return angle; 
}
//----------------------------------------------------------------------------
void AliBtoJPSItoEle::ComputeWgts() {
  // calculate the weights for PID


  // assignement of the weights from PID
  fWgtJPsi    = fTagEl[0]*fTagEl[1]; // both assumed to be electrons 

  
  // if(fWgtJPsi<0.) cerr<<"AliBtoJPSItoEle::ComputeWgts()  Negative weight!!!\n";
  

  return;
}
//----------------------------------------------------------------------------
void AliBtoJPSItoEle::CorrectWgt4BR(Double_t factor) {
  // correct weights of background from charm 

  fWgtJPsi    *= factor;

  return;
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::CosPointing() const {
  // cosine of pointing angle in space

  TVector3 mom(Px(),Py(),Pz());
  TVector3 flight(fV2x-fV1x,fV2y-fV1y,fV2z-fV1z);

  Double_t pta = mom.Angle(flight);

  return TMath::Cos(pta); 
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::CosPointingXY() const {
  // cosine of pointing angle in transverse plane

  TVector3 momXY(Px(),Py(),0.);
  TVector3 flightXY(fV2x-fV1x,fV2y-fV1y,0.);

  Double_t ptaXY = momXY.Angle(flightXY);

  return TMath::Cos(ptaXY); 
}
//----------------------------------------------------------------------------
void AliBtoJPSItoEle::CosThetaStar(Double_t &ctsJPsi) const {
  // cosine of decay angle in the J/Psi rest frame (of the negative electron)

  Double_t pStar = TMath::Sqrt(TMath::Power(kMJPsi*kMJPsi-2.*kMe*kMe,2.)-4.*kMe*kMe*kMe*kMe)/(2.*kMJPsi);

  Double_t beta = P()/Energy();
  Double_t gamma = Energy()/kMJPsi;

  ctsJPsi = (Ql(1)/gamma-beta*TMath::Sqrt(pStar*pStar+kMe*kMe))/pStar;
  // if(ctsJPsi > 1.)  { cerr<<"AliBtoJPSItoEle::CosThetaStar: > 1 "<<ctsJPsi<<"!\n"; }
  // if(ctsJPsi < -1.) { cerr<<"AliBtoJPSItoEle::CosThetaStar: < -1 "<<ctsJPsi<<"!\n"; }

  return;
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::Eta() const {
  // pseudorapidity of the J/Psi

  Double_t theta = TMath::Pi()/2.-TMath::ATan2(Pz(),Pt());
  Double_t eta = -TMath::Log(TMath::Tan(theta/2.));
  return eta;
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::EtaChild(Int_t child) const {
  // pseudorapidity of the decay tracks

  Double_t theta = TMath::Pi()/2.-TMath::ATan2(fPz[child],PtChild(child));
  Double_t eta = -TMath::Log(TMath::Tan(theta/2.));
  return eta;
}
//----------------------------------------------------------------------------
void AliBtoJPSItoEle::GetWgts(Double_t &WgtJPsi) 
  const {
  // returns the weights for pid

    WgtJPsi    = fWgtJPsi; 

  return;
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::ImpPar() const {
  // J/Psi impact parameter in the bending plane
  
    Double_t k = -(fV2x-fV1x)*Px()-(fV2y-fV1y)*Py();
    k /= Pt()*Pt();
    Double_t dx = fV2x-fV1x+k*Px();
    Double_t dy = fV2y-fV1y+k*Py();
    Double_t absDD = TMath::Sqrt(dx*dx+dy*dy);
    TVector3 mom(Px(),Py(),Pz());
    TVector3 flight(fV2x-fV1x,fV2y-fV1y,fV2z-fV1z);
    TVector3 cross = mom.Cross(flight);
    return (cross.Z()>0. ? absDD : -absDD);
}
//----------------------------------------------------------------------------
void AliBtoJPSItoEle::InvMass(Double_t &mJPsi) const {
  // invariant mass as J/Psi

  Double_t energy[2];

  // J/psi -> e- e+
  energy[1] = TMath::Sqrt(kMe*kMe+PChild(1)*PChild(1));
  energy[0] = TMath::Sqrt(kMe*kMe+PChild(0)*PChild(0));

  mJPsi = TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-P()*P());
    
  return;

}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::Ql(Int_t child) const {
  // longitudinal momentum of decay tracks w.r.t. to J/Psi momentum

  Double_t qL;
  TVector3 mom(fPx[child],fPy[child],fPz[child]);
  TVector3 momJPsi(Px(),Py(),Pz());

  qL = mom.Dot(momJPsi)/momJPsi.Mag();

  return qL ;
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::Qt() const {
  // transverse momentum of decay tracks w.r.t. to JPsi momentum  

  TVector3 mom0(fPx[0],fPy[0],fPz[0]);
  TVector3 momJPsi(Px(),Py(),Pz());

  return mom0.Perp(momJPsi);
}
//----------------------------------------------------------------------------
Bool_t AliBtoJPSItoEle::Select(const Double_t* cuts,Int_t& okB) 
  const {
//
// This function compares the B candidates with a set of cuts:
//
// cuts[0] = inv. mass half width [GeV]   
// cuts[1] = dca [micron]
// cuts[2] = cosThetaStar 
// cuts[3] = pTP [GeV/c]
// cuts[4] = pTN [GeV/c]
// cuts[5] = d0P [micron]   upper limit!
// cuts[6] = d0N [micron]  upper limit!
// cuts[7] = d0d0 [micron^2]
// cuts[8] = cosThetaPoint
//
// If the candidate doesn't pass the cuts it sets the weight to 0
// and return kFALSE
//
  Double_t mJPsi,ctsJPsi;
  okB=1; 

  if(PtChild(1) < cuts[3] || PtChild(0) < cuts[4]) okB = 0;
  if(!okB) return kFALSE;

  if(TMath::Abs(Getd0Child(1)) > cuts[5] || 
     TMath::Abs(Getd0Child(0)) > cuts[6]) okB = 0;
  if(!okB) return kFALSE;

  if(GetDCA() > cuts[1]) { okB = 0; return kFALSE; }

  InvMass(mJPsi);
  if(TMath::Abs(mJPsi-kMJPsi)    > cuts[0]) okB = 0;
  if(!okB) return kFALSE;

  CosThetaStar(ctsJPsi);
  if(TMath::Abs(ctsJPsi)    > cuts[2]) okB = 0;
  if(!okB) return kFALSE;

  if(ProdImpParams() > cuts[7]) { okB = 0; return kFALSE; }

  if(CosPointing()   < cuts[8]) { okB = 0; return kFALSE; }

  return kTRUE;
}
//-----------------------------------------------------------------------------
void AliBtoJPSItoEle::SetPIDresponse(Double_t resp0[5],Double_t resp1[5]) {
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
Double_t AliBtoJPSItoEle::GetTagEl(Int_t child) const {
  // Get tag probabilities for electrons
                                                                                                                          
  if(child!=0 && child !=1) return -1;

  return fTagEl[child];
}
//-----------------------------------------------------------------------------
void AliBtoJPSItoEle::GetPIDresponse(Double_t resp0[5],Double_t resp1[5]) const {
  // Get combined PID detector response probabilities
                                                                                                                          
  resp0[0] = fPIDrespEl[0];
  resp1[0] = fPIDrespEl[1]; 
  resp0[1] = fPIDrespMu[0];
  resp1[1] = fPIDrespMu[1];
  resp0[2] = fPIDrespPi[0];
  resp1[2] = fPIDrespPi[1];
  resp0[3] = fPIDrespKa[0];
  resp1[3] = fPIDrespKa[1];
  resp0[4] = fPIDrespPr[0];
  resp1[4] = fPIDrespPr[1];
                                                                                                                          
  return;
}
//----------------------------------------------------------------------------
Double_t AliBtoJPSItoEle::TRDTPCCombinedPIDParametrization(Double_t p) const {


  // a first raw parametrization of the probability to misidentify a charged pion as electron as a 
  // function of the momentum, as given by the combined TPC and TRD response. 
  //  PID cuts are set such that the probability for correct electron id is 90% in each of the two 
  //    detectors
  
// first estimate based on parameterization in the B-> single electron analysis
  Double_t value=0;
  Double_t p1 =11.;
  Double_t p2=0.00007;
  Double_t p3=0.007;
  value=p2+p3*(1.-exp(-TMath::Power(p/p1,4.)));

  value/=0.01; // here remove from TPC+TRD the TRD contribution estimated to be 0.01
// Better estimation based on TRD test beam (as presented by Andrea at Munster)
//  if (p<10.) value*=(1.32-0.18*p+0.076*p*p-0.0037*p*p*p)/100.;
//  if (p>10.) value*=(0.48+0.287*p)/100.;
// From Silvia MASCIOCCHI (Darmstadt) at PWG3 AliceWeek October 2007
  if (p<10.) value*=(0.44-0.06*p+0.031*p*p-0.0008*p*p*p)/100.;
  if (p>10.) value*=(-0.67+0.28*p)/100.;

  return value;
}
//----------------------------------------------------------------------------




