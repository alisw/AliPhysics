#include <iostream.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include "AliD0toKpi.h"

ClassImp(AliD0toKpi)
  /*******************************************************************
   *                                                                 *
   *  Reconstructed D^0 -> K^- pi^+   candidate class                * 
   *                                                                 *
   *  origin: A. Dainese    andrea.dainese@pd.infn.it                *  
   *******************************************************************/

//________________________________________________________
AliD0toKpi::AliD0toKpi() {

  // Default constructor
  fSignal = kFALSE;

  fEvent = 0;

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
  fMum[0] =0;
  fMum[1] =0;


  fWgtAD0 = fWgtAD0bar = fWgtBD0 = fWgtBD0bar = fWgtCD0 = fWgtCD0bar = 0;

}

//________________________________________________________
AliD0toKpi::AliD0toKpi(Bool_t sgn,Int_t ev,
		       Double_t v1[3],Double_t v2[3], 
		       Double_t dca,
		       Double_t mom[6],Double_t d0[2],
		       Int_t pdg[2],Int_t mum[2]) {

  // Standard constructor
  fSignal = sgn;

  fEvent = ev;

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

  fd0[0] =  d0[0];
  fd0[1] = d0[1];

  fPdg[0] = pdg[0];
  fPdg[1] = pdg[1];
  fMum[0] = mum[0];
  fMum[1] = mum[1];


  fWgtAD0 = fWgtAD0bar = fWgtBD0 = fWgtBD0bar = fWgtCD0 = fWgtCD0bar = 0;
}

//_____________________________________________________________________
Double_t AliD0toKpi::EtaChild(Int_t child) const {

  Double_t theta = TMath::Pi()/2.-TMath::ATan2(fPz[child],PtChild(child));
  Double_t eta = -TMath::Log(TMath::Tan(theta/2.));
  return eta;
}
//_____________________________________________________________________
Double_t AliD0toKpi::Eta() const {

  Double_t theta = TMath::Pi()/2.-TMath::ATan2(Pz(),Pt());
  Double_t eta = -TMath::Log(TMath::Tan(theta/2.));
  return eta;
}

//_____________________________________________________________________
Double_t AliD0toKpi::CPta() const {

  TVector3 mom(Px(),Py(),Pz());
  TVector3 flight(fV2x-fV1x,fV2y-fV1y,fV2z-fV1z);

  Double_t pta = mom.Angle(flight);

  return TMath::Cos(pta); 
}
//_____________________________________________________________________
Double_t AliD0toKpi::CPtaXY() const {

  TVector3 momXY(Px(),Py(),0.);
  TVector3 flightXY(fV2x-fV1x,fV2y-fV1y,0.);

  Double_t ptaXY = momXY.Angle(flightXY);

  return TMath::Cos(ptaXY); 
}
//_____________________________________________________________________
Double_t AliD0toKpi::ChildrenRelAngle() const {

  TVector3 mom0(fPx[0],fPy[0],fPz[0]);
  TVector3 mom1(fPx[1],fPy[1],fPz[1]);

  Double_t angle = mom0.Angle(mom1);

  return angle; 
}
//_____________________________________________________________________
void AliD0toKpi::DrawPIDinTOF() const {

  TH2F* framePi = new TH2F("framePi","Tag probabilities for PIONS",2,0,2.5,2,0,1);
  framePi->SetXTitle("p [GeV/c]"); 
  framePi->SetStats(0);
  TH2F* frameK = new TH2F("frameK","Tag probabilities for KAONS",2,0,2.5,2,0,1);
  frameK->SetXTitle("p [GeV/c]");
  frameK->SetStats(0);
  TH2F* frameP = new TH2F("frameP","Tag probabilities for PROTONS",2,0,4.5,2,0,1);
  frameP->SetXTitle("p [GeV/c]");
  frameP->SetStats(0);

  TH1F* hPiPi = new TH1F("hPiPi","Tag probabilities for PIONS",kPiBins,0,2.5);
  TH1F* hPiNid = new TH1F("hPiNid","Tag probabilities for PIONS",kPiBins,0,2.5);

  TH1F* hKK = new TH1F("hKK","Tag probabilities for KAONS",kKBins,0,2.5);
  TH1F* hKNid = new TH1F("hKNid","Tag probabilities for KAONS",kKBins,0,2.5);
  TH1F* hKPi = new TH1F("hKPi","Tag probabilities for KAONS",kKBins,0,2.5);

  TH1F* hPP = new TH1F("hPP","Tag probabilities for PROTONS",kPBins,0,4.5);
  TH1F* hPNid = new TH1F("hPNid","Tag probabilities for PROTONS",kPBins,0,4.5);
  TH1F* hPPi = new TH1F("hPPi","Tag probabilities for PROTONS",kPBins,0,4.5);

  for(Int_t i=1; i<=kPiBins; i++) {
    hPiPi->SetBinContent(i,kPiTagPi[i-1]);
    hPiNid->SetBinContent(i,kPiTagPi[i-1]+kPiTagNid[i-1]);

    hKK->SetBinContent(i,kKTagK[i-1]);
    hKPi->SetBinContent(i,kKTagK[i-1]+kKTagPi[i-1]);
    hKNid->SetBinContent(i,kKTagK[i-1]+kKTagPi[i-1]+kKTagNid[i-1]);
  }
  for(Int_t i=1; i<=kPBins; i++) {    
    hPP->SetBinContent(i,kPTagP[i-1]);
    hPPi->SetBinContent(i,kPTagP[i-1]+kPTagPi[i-1]);
    hPNid->SetBinContent(i,kPTagP[i-1]+kPTagPi[i-1]+kPTagNid[i-1]);
  }

  TCanvas* c = new TCanvas("c","PID",0,0,1000,400);
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
//___________________________________________________________________________
Double_t AliD0toKpi::LinearInterpolation(Double_t p,Int_t nBins,Double_t Bin,const Double_t *values) const {

  Double_t value=0;
  Double_t slope;

  if(p<0.5*Bin) {
    value = values[0];
  } else if(p>=(nBins-0.5)*Bin) {
    //slope = (values[nBins-1]-values[nBins-2])/Bin;
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

  if(value<0) value=0.;
  if(value>1) value=1.;
  return value;
}
//______________________________________________________________________
void AliD0toKpi::ComputeWgts() {

  Double_t probTagPi[2]  = {0.,0.};
  Double_t probTagK[2]   = {0.,0.};
  Double_t probTagNid[2] = {0.,0.};
  Double_t probTagP[2]   = {0.,0.};

  // tagging of the positive track
  if(TMath::Abs(fPdg[0])==211 || TMath::Abs(fPdg[0])==13 
     || TMath::Abs(fPdg[0])==11) { // pion,muon,electron
    probTagPi[0]  = LinearInterpolation(PChild(0),kPiBins,kPiBinWidth,kPiTagPi);
    probTagNid[0] = 1.-probTagPi[0];
    probTagK[0]   = 0.;
    probTagP[0]   = 0.;
  } 
  if(TMath::Abs(fPdg[0])==321) { // kaon
    probTagK[0]   = LinearInterpolation(PChild(0),kKBins,kKBinWidth,kKTagK);
    probTagNid[0] = LinearInterpolation(PChild(0),kKBins,kKBinWidth,kKTagNid);
    if((probTagNid[0]+probTagK[0])>1.) probTagNid[0] = 1.-probTagK[0];
    probTagPi[0] = 1.-probTagNid[0]-probTagK[0];
    probTagP[0] = 0.;
  } 
  if(TMath::Abs(fPdg[0])==2212) { // proton
    probTagP[0]  = LinearInterpolation(PChild(0),kPBins,kPBinWidth,kPTagP);
    probTagNid[0] = LinearInterpolation(PChild(0),kPBins,kPBinWidth,kPTagNid);
    if((probTagNid[0]+probTagP[0])>1.) probTagNid[0] = 1.-probTagP[0];
    probTagPi[0] = 1.-probTagNid[0]-probTagP[0];
    probTagK[0]   = 0.;
  } 


  // tagging of the negative track
  if(TMath::Abs(fPdg[1])==211 || TMath::Abs(fPdg[1])==13 
     || TMath::Abs(fPdg[1])==11) { // pion,muon,electron
    probTagPi[1]  = LinearInterpolation(PChild(1),kPiBins,kPiBinWidth,kPiTagPi);
    probTagNid[1] = 1.-probTagPi[1];
    probTagK[1]   = 0.;
    probTagP[1]   = 0.;
  } 
  if(TMath::Abs(fPdg[1])==321) { // kaon
    probTagK[1]   = LinearInterpolation(PChild(1),kKBins,kKBinWidth,kKTagK);
    probTagNid[1] = LinearInterpolation(PChild(1),kKBins,kKBinWidth,kKTagNid);
    if((probTagNid[1]+probTagK[1])>1.) probTagNid[1] = 1.-probTagK[1];
    probTagPi[1] = 1.-probTagNid[1]-probTagK[1];
    probTagP[1] = 0.;
  } 
  if(TMath::Abs(fPdg[1])==2212) { // proton
    probTagP[1]  = LinearInterpolation(PChild(1),kPBins,kPBinWidth,kPTagP);
    probTagNid[1] = LinearInterpolation(PChild(1),kPBins,kPBinWidth,kPTagNid);
    if((probTagNid[1]+probTagP[1])>1.) probTagNid[1] = 1.-probTagP[1];
    probTagPi[1] = 1.-probTagNid[1]-probTagP[1];
    probTagK[1]   = 0.;
  } 

  // assignement of the weights from PID
  fWgtAD0    = probTagK[1]*(probTagPi[0]+probTagNid[0]);
  fWgtAD0bar = probTagK[0]*(probTagPi[1]+probTagNid[1]);
  fWgtBD0    = probTagPi[0]*probTagNid[1];
  fWgtBD0bar = probTagPi[1]*probTagNid[0];
  fWgtCD0    = probTagNid[0]*probTagNid[1];
  fWgtCD0bar = probTagNid[0]*probTagNid[1];

  // pt-dependent weight to reproduce pt distr. given by NLO QCD (MNR)
  // [only for candidates with at least one track coming from a D0]
  SetPtWgts();



  /*
  for(Int_t j=0;j<2;j++) cerr<<" PDG = "<<GetPdgChild(j)<<" p = "<<PChild(j)<<" TagPi = "<<ProbTagPi[j]<<" TagK = "<<ProbTagK[j]<<" TagNid = "<<ProbTagNid[j]<<endl;

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
//____________________________________________________________________________
void AliD0toKpi::CorrectWgt4BR(Double_t factor) {

  fWgtAD0    *= factor;
  fWgtAD0bar *= factor;
  fWgtBD0    *= factor;
  fWgtBD0bar *= factor;
  fWgtCD0    *= factor;
  fWgtCD0bar *= factor;

  return;
}
//____________________________________________________________________________
void AliD0toKpi::SetPtWgts() {

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

  return;
}
//____________________________________________________________________________
void AliD0toKpi::GetWgts(Double_t &WgtD0,Double_t &WgtD0bar,Option_t *sample) const {

  const char *sampleA = strstr(sample,"A");
  const char *sampleB = strstr(sample,"B");
  const char *sampleC = strstr(sample,"C");

  if(sampleA) { WgtD0 = fWgtAD0;  WgtD0bar = fWgtAD0bar; }
  if(sampleB) { WgtD0 = fWgtBD0;  WgtD0bar = fWgtBD0bar; }
  if(sampleC) { WgtD0 = fWgtCD0;  WgtD0bar = fWgtCD0bar; }

  if(fSignal) {
    if(fMum[0]==421)  WgtD0bar = 0.;
    if(fMum[0]==-421) WgtD0 = 0.; 
  }

  return;
}
//____________________________________________________________________________
void AliD0toKpi::InvMass(Double_t &mD0,Double_t &mD0bar) const {
    
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
//________________________________________________________________
Double_t AliD0toKpi::qT() const {

  TVector3 mom0(fPx[0],fPy[0],fPz[0]);
  TVector3 momD(Px(),Py(),Pz());

  return mom0.Perp(momD);
}
//________________________________________________________________
Double_t AliD0toKpi::qL(Int_t child) const {

  Double_t qL;
  TVector3 mom(fPx[child],fPy[child],fPz[child]);
  TVector3 momD(Px(),Py(),Pz());

  qL = mom.Dot(momD)/momD.Mag();

  return qL ;
}
//________________________________________________________________
void AliD0toKpi::CosThetaStar(Double_t &ctsD0,Double_t &ctsD0bar) const {

  Double_t pStar = TMath::Sqrt(TMath::Power(kMD0*kMD0-kMK*kMK-kMPi*kMPi,2.)-4.*kMK*kMK*kMPi*kMPi)/(2.*kMD0);

  Double_t beta = P()/Energy();
  Double_t gamma = Energy()/kMD0;

  ctsD0 = (qL(1)/gamma-beta*TMath::Sqrt(pStar*pStar+kMK*kMK))/pStar;
// if(ctsD0 > 1.)  { cerr<<"AliD0toKpi::CosThetaStar: > 1 "<<ctsD0<<"!\n"; }
// if(ctsD0 < -1.) { cerr<<"AliD0toKpi::CosThetaStar: < -1 "<<ctsD0<<"!\n"; }

  ctsD0bar = (qL(0)/gamma-beta*TMath::Sqrt(pStar*pStar+kMK*kMK))/pStar;
// if(ctsD0bar > 1.)  { cerr<<"AliD0toKpi::CosThetaStar: > 1 "<<ctsD0bar<<"!\n"; }
// if(ctsD0bar < -1.) { cerr<<"AliD0toKpi::CosThetaStar: < -1 "<<ctsD0bar<<"!\n";} 

  return;
}
//____________________________________________________________________
Bool_t AliD0toKpi::Select(const Double_t* cuts,Int_t& okD0,Int_t& okD0bar) const {
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
// If the the D0/D0bar doesn't pass the cuts it sets the weights to 0
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
  if(TMath::Abs(mD0-kMD0) > cuts[0])    okD0 = 0;
  if(TMath::Abs(mD0bar-kMD0) > cuts[0]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  CosThetaStar(ctsD0,ctsD0bar);
  if(TMath::Abs(ctsD0) > cuts[2])    okD0 = 0;
  if(TMath::Abs(ctsD0bar) > cuts[2]) okD0bar = 0;
  if(!okD0 && !okD0bar) return kFALSE;

  if(d0d0() > cuts[7]) { okD0 = okD0bar = 0; return kFALSE; }

  if(CPta() < cuts[8]) { okD0 = okD0bar = 0; return kFALSE; }

  return kTRUE;
}







