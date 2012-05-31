/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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

/* $Id: AliBalancePsi.cxx 54125 2012-01-24 21:07:41Z miweber $ */

//-----------------------------------------------------------------
//           Balance Function class
//   This is the class to deal with the Balance Function wrt Psi analysis
//   Origin: Panos Christakoglou, Nikhef, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------


//ROOT
#include <Riostream.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TString.h>

#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliTHn.h"

#include "AliBalancePsi.h"

ClassImp(AliBalancePsi)

//____________________________________________________________________//
AliBalancePsi::AliBalancePsi() :
  TObject(), 
  bShuffle(kFALSE),
  fAnalysisLevel("ESD"),
  fAnalyzedEvents(0) ,
  fCentralityId(0) ,
  fCentStart(0.),
  fCentStop(0.),
  fHistP(0),
  fHistN(0),
  fHistPN(0),
  fHistNP(0),
  fHistPP(0),
  fHistNN(0),
  fPsiInterval(15.) {
  // Default constructor
}

//____________________________________________________________________//
AliBalancePsi::AliBalancePsi(const AliBalancePsi& balance):
  TObject(balance), bShuffle(balance.bShuffle), 
  fAnalysisLevel(balance.fAnalysisLevel),
  fAnalyzedEvents(balance.fAnalyzedEvents), 
  fCentralityId(balance.fCentralityId),
  fCentStart(balance.fCentStart),
  fCentStop(balance.fCentStop),
  fHistP(balance.fHistP),
  fHistN(balance.fHistN),
  fHistPN(balance.fHistPN),
  fHistNP(balance.fHistNP),
  fHistPP(balance.fHistPP),
  fHistNN(balance.fHistNN),
  fPsiInterval(balance.fPsiInterval) {
  //copy constructor
}

//____________________________________________________________________//
AliBalancePsi::~AliBalancePsi() {
  // Destructor
  delete fHistP;
  delete fHistN;
  delete fHistPN;
  delete fHistNP;
  delete fHistPP;
  delete fHistNN;
}

//____________________________________________________________________//
void AliBalancePsi::InitHistograms() {
  // single particle histograms
  Int_t anaSteps   = 1;       // analysis steps
  Int_t iBinSingle[nTrackVariablesSingle];        // binning for track variables
  Double_t* dBinsSingle[nTrackVariablesSingle];   // bins for track variables  
  TString axisTitleSingle[nTrackVariablesSingle]; // axis titles for track variables
  
  // two particle histograms
  Int_t iBinPair[nTrackVariablesPair];         // binning for track variables
  Double_t* dBinsPair[nTrackVariablesPair];    // bins for track variables  
  TString axisTitlePair[nTrackVariablesPair];  // axis titles for track variables

  //centrality
  /*const Int_t kNCentralityBins = 9;
  Double_t centralityBins[kNCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};
  iBinSingle[0]       = kNCentralityBins;
  dBinsSingle[0]      = centralityBins;
  axisTitleSingle[0]  = "Centrality percentile [%]"; 
  iBinPair[0]       = kNCentralityBins;
  dBinsPair[0]      = centralityBins;
  axisTitlePair[0]  = "Centrality percentile [%]"; */

  //Psi_2
  const Int_t kNPsi2Bins = 3;
  Double_t psi2Bins[kNPsi2Bins+1] = {-0.5,0.5,1.5,2.5};
  iBinSingle[0]       = kNPsi2Bins;
  dBinsSingle[0]      = psi2Bins;
  axisTitleSingle[0]  = "#phi - #Psi_{2} (a.u.)";
  iBinPair[0]       = kNPsi2Bins;
  dBinsPair[0]      = psi2Bins;
  axisTitlePair[0]  = "#phi - #Psi_{2} (a.u.)"; 
  
   // delta eta
  const Int_t kNDeltaEtaBins = 80;
  Double_t deltaEtaBins[kNDeltaEtaBins+1];
  for(Int_t i = 0; i < kNDeltaEtaBins+1; i++)
    deltaEtaBins[i] = -2.0 + i * 0.05;
  iBinPair[1]       = kNDeltaEtaBins;
  dBinsPair[1]      = deltaEtaBins;
  axisTitlePair[1]  = "#Delta #eta"; 

   // delta phi
  const Int_t kNDeltaPhiBins = 72;
  Double_t deltaPhiBins[kNDeltaPhiBins+1];
  for(Int_t i = 0; i < kNDeltaPhiBins+1; i++){
    deltaPhiBins[i] = -180.0 + i * 5.;
  } 
  iBinPair[2]       = kNDeltaPhiBins;
  dBinsPair[2]      = deltaPhiBins;
  axisTitlePair[2]  = "#Delta #phi (#circ)"; 

  // pt(trigger-associated)
  const Int_t kNPtBins = 40;
  Double_t ptBins[kNPtBins+1];
  for(Int_t i = 0; i < kNPtBins+1; i++){
    ptBins[i] = 0.0 + i * 0.5;
   } 
  iBinSingle[1]       = kNPtBins;
  dBinsSingle[1]      = ptBins;
  axisTitleSingle[1]  = "p_{t}^{trig.} (GeV/c)"; 

  iBinPair[3]       = kNPtBins;
  dBinsPair[3]      = ptBins;
  axisTitlePair[3]  = "p_{t}^{trig.} (GeV/c)"; 

  iBinPair[4]       = kNPtBins;
  dBinsPair[4]      = ptBins;
  axisTitlePair[4]  = "p_{t}^{assoc.} (GeV/c)";   

  TString histName;
  //+ triggered particles
  histName = "fHistP"; 
  if(bShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistP = new AliTHn(histName.Data(),histName.Data(),anaSteps,nTrackVariablesSingle,iBinSingle);
  for (Int_t j=0; j<nTrackVariablesSingle; j++) {
    fHistP->SetBinLimits(j, dBinsSingle[j]);
    fHistP->SetVarTitle(j, axisTitleSingle[j]);
  }

  //- triggered particles
  histName = "fHistN"; 
  if(bShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistN = new AliTHn(histName.Data(),histName.Data(),anaSteps,nTrackVariablesSingle,iBinSingle);
  for (Int_t j=0; j<nTrackVariablesSingle; j++) {
    fHistN->SetBinLimits(j, dBinsSingle[j]);
    fHistN->SetVarTitle(j, axisTitleSingle[j]);
  }
  
  //+- pairs
  histName = "fHistPN";
  if(bShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistPN = new AliTHn(histName.Data(),histName.Data(),anaSteps, nTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<nTrackVariablesPair; j++) {
    fHistPN->SetBinLimits(j, dBinsPair[j]);
    fHistPN->SetVarTitle(j, axisTitlePair[j]);
  }

  //-+ pairs
  histName = "fHistNP";
  if(bShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistNP = new AliTHn(histName.Data(),histName.Data(),anaSteps, nTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<nTrackVariablesPair; j++) {
    fHistNP->SetBinLimits(j, dBinsPair[j]);
    fHistNP->SetVarTitle(j, axisTitlePair[j]);
  }

  //++ pairs
  histName = "fHistPP";
  if(bShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistPP = new AliTHn(histName.Data(),histName.Data(),anaSteps, nTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<nTrackVariablesPair; j++) {
    fHistPP->SetBinLimits(j, dBinsPair[j]);
    fHistPP->SetVarTitle(j, axisTitlePair[j]);
  }

  //-- pairs
  histName = "fHistNN";
  if(bShuffle) histName.Append("_shuffle");
  if(fCentralityId) histName += fCentralityId.Data();
  fHistNN = new AliTHn(histName.Data(),histName.Data(),anaSteps, nTrackVariablesPair, iBinPair);
  for (Int_t j=0; j<nTrackVariablesPair; j++) {
    fHistNN->SetBinLimits(j, dBinsPair[j]);
    fHistNN->SetVarTitle(j, axisTitlePair[j]);
  }
  Printf("Finished setting up the AliTHn");
}

//____________________________________________________________________//
void AliBalancePsi::CalculateBalance(Double_t gReactionPlane, 
				     vector<Double_t> **chargeVector) {
  // Calculates the balance function
  fAnalyzedEvents++;
  Int_t i = 0 , j = 0;
  
  // Initialize histograms if not done yet
  if(!fHistPN){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("This works only in local mode --> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    InitHistograms();
  }

  Double_t trackVarsSingle[nTrackVariablesSingle];
  Double_t trackVarsPair[nTrackVariablesPair];

  Int_t gNtrack = chargeVector[0]->size();
  //Printf("(AliBalancePsi) Number of tracks: %d",gNtrack);
  Double_t gPsiMinusPhi = 0.;
  Double_t dy = 0., deta = 0.;
  Double_t qLong = 0., qOut = 0., qSide = 0., qInv = 0.;
  Double_t dphi = 0.;

  Double_t charge1  = 0;
  Double_t eta1 = 0., rap1 = 0.;
  Double_t px1 = 0., py1 = 0., pz1 = 0.;
  Double_t pt1 = 0.;
  Double_t energy1 = 0.;
  Double_t phi1    = 0.;

  Double_t charge2  = 0;
  Double_t eta2 = 0., rap2 = 0.;
  Double_t px2 = 0., py2 = 0., pz2 = 0.;
  Double_t pt2 = 0.;
  Double_t energy2 = 0.;
  Double_t phi2    = 0.;
  Double_t gPsiMinusPhiBin = -10.;
  
  for(i = 0; i < gNtrack; i++) {
    gPsiMinusPhiBin = -10.;
    charge1 = chargeVector[0]->at(i);
    rap1    = chargeVector[1]->at(i);
    eta1    = chargeVector[2]->at(i);
    phi1    = chargeVector[3]->at(i);
    px1     = chargeVector[4]->at(i);
    py1     = chargeVector[5]->at(i);
    pz1     = chargeVector[6]->at(i);
    pt1     = chargeVector[7]->at(i);
    energy1 = chargeVector[8]->at(i);
    gPsiMinusPhi   = TMath::Abs(phi1 - gReactionPlane);
    if((gPsiMinusPhi <= 7.5)||
       ((172.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 187.5)))
      gPsiMinusPhiBin = 0.0;
    else if(((37.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 52.5))||
	    ((127.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 142.5))||
	    ((217.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 232.5))||
	    ((307.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 322.5)))
      gPsiMinusPhiBin = 1.0;
    else if(((82.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 97.5))||
	    ((262.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 277.5)))
      gPsiMinusPhiBin = 2.0;
    else continue;
    //cout<<"phi: "<<phi1<<" - Psi: "<<gReactionPlane<<" - phi - Psi: "<<gPsiMinusPhi<<" - Bin: "<<gPsiMinusPhiBin<<endl; 
 
    //if(gPsiMinusPhi > 180.) gPsiMinusPhi = 360. - gPsiMinusPhi;
    //if(gPsiMinusPhi < -fPsiInterval/2) gPsiMinusPhi = 360. + gPsiMinusPhi;
    
    //trackVarsSingle[0]    =  fCentrality;
    trackVarsSingle[0]    =  gPsiMinusPhiBin;
    //trackVarsSingle[2]    =  eta1;
    //trackVarsSingle[3]    =  rap1;
    //trackVarsSingle[3]    =  phi1;  
    trackVarsSingle[1]    =  pt1;  

    //Printf("Track(a) %d - phi-Psi: %lf",i+1,trackVarsSingle[1]);
    //fill single particle histograms
    if(charge1 > 0)  
      fHistP->Fill(trackVarsSingle,0,1.); 
    else if(charge1 < 0)            
      fHistN->Fill(trackVarsSingle,0,1.); 

    for(j = 0; j < i; j++) {
      charge2 = chargeVector[0]->at(j);
      rap2    = chargeVector[1]->at(j);
      eta2    = chargeVector[2]->at(j);
      phi2    = chargeVector[3]->at(j);
      px2     = chargeVector[4]->at(j);
      py2     = chargeVector[5]->at(j);
      pz2     = chargeVector[6]->at(j);
      pt2     = chargeVector[7]->at(j);
      energy2 = chargeVector[8]->at(j);
      //Printf("Track(b) %d - pt: %lf",j+1,pt2);
      
      // filling the arrays
      // RAPIDITY 
      dy = rap1 - rap2;
      
      // Eta
      deta = eta1 - eta2;
      
      //qlong
      Double_t eTot = energy1 + energy2;
      Double_t pxTot = px1 + px2;
      Double_t pyTot = py1 + py2;
      Double_t pzTot = pz1 + pz2;
      Double_t q0Tot = energy1 - energy2;
      Double_t qxTot = px1 - px2;
      Double_t qyTot = py1 - py2;
      Double_t qzTot = pz1 - pz2;
      
      Double_t eTot2 = eTot*eTot;
      Double_t pTot2 = pxTot*pxTot + pyTot*pyTot + pzTot*pzTot;
      Double_t pzTot2 = pzTot*pzTot;
      
      Double_t q0Tot2 = q0Tot*q0Tot;
      Double_t qTot2  = qxTot*qxTot + qyTot*qyTot + qzTot*qzTot;
      
      Double_t snn    = eTot2 - pTot2;
      Double_t ptTot2 = pTot2 - pzTot2 ;
      Double_t ptTot  = TMath::Sqrt( ptTot2 );
      
      qLong = TMath::Abs(eTot*qzTot - pzTot*q0Tot)/TMath::Sqrt(snn + ptTot2);
      
      //qout
      qOut = TMath::Sqrt(snn/(snn + ptTot2)) * TMath::Abs(pxTot*qxTot + pyTot*qyTot)/ptTot;
      
      //qside
      qSide = TMath::Abs(pxTot*qyTot - pyTot*qxTot)/ptTot;
      
      //qinv
      qInv = TMath::Sqrt(TMath::Abs(-q0Tot2 + qTot2 ));
      
      //phi
      dphi = phi2 - phi1;
      if(dphi < -180.) dphi += 360.;  //dphi should be between -180 and 180!
      else if(dphi > 180.) dphi -= 360.; //dphi should be between -180 and 180!
      
      //trackVarsPair[0]    =  fCentrality;             
      trackVarsPair[0]    =  gPsiMinusPhiBin;
      //trackVarsPair[2]    =  eta1;
      //trackVarsPair[3]    =  rap1;
      //trackVarsPair[4]    =  phi1;
      trackVarsPair[1]    =  deta;
      //trackVarsPair[8]    =  dy;
      trackVarsPair[2]    =  dphi;
      trackVarsPair[3]    =  pt1;
      trackVarsPair[4]    =  pt2;
      //trackVarsPair[10]   =  qSide;
      //trackVarsPair[11]   =  qOut;
      //trackVarsPair[12]   =  qLong;
      //trackVarsPair[13]   =  qInv;
      
      if( charge1 > 0 && charge2 < 0)  fHistPN->Fill(trackVarsPair,0,1.); 
      else if( charge1 < 0 && charge2 > 0)  fHistNP->Fill(trackVarsPair,0,1.); 
      else if( charge1 > 0 && charge2 > 0)  fHistPP->Fill(trackVarsPair,0,1.); 
      else if( charge1 < 0 && charge2 < 0)  fHistNN->Fill(trackVarsPair,0,1.); 
      else AliWarning("Wrong charge combination!");
      
    }//end of 2nd particle loop
  }//end of 1st particle loop
}  

//____________________________________________________________________//
TH1D *AliBalancePsi::GetBalanceFunctionHistogram(Int_t iVariableSingle,
						 Int_t iVariablePair,
						 Double_t psiMin, 
						 Double_t psiMax) {
  //Returns the BF histogram, extracted from the 6 AliTHn objects 
  //(private members) of the AliBalancePsi class.
  //iVariableSingle: 0(phi-Psi), 1(pt-trigger)
  //iVariablePair: 0(phi-Psi) 1(Delta eta), 2(Delta phi), 3(pt-trigger), 4(pt-associated
  TString gAnalysisType[ANALYSIS_TYPES] = {"y","eta","phi","qlong","qout","qside","qinv"};
  TString histName = "gHistBalanceFunctionHistogram";
  histName += gAnalysisType[iVariablePair];

  // Psi_2
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 

  //Printf("P:%lf - N:%lf - PN:%lf - NP:%lf - PP:%lf - NN:%lf",fHistP->GetEntries(0),fHistN->GetEntries(0),fHistPN->GetEntries(0),fHistNP->GetEntries(0),fHistPP->GetEntries(0),fHistNN->GetEntries(0));

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH1D* hTemp1 = (TH1D*)fHistPN->Project(0,iVariablePair);
  TH1D* hTemp2 = (TH1D*)fHistNP->Project(0,iVariablePair);
  TH1D* hTemp3 = (TH1D*)fHistPP->Project(0,iVariablePair);
  TH1D* hTemp4 = (TH1D*)fHistNN->Project(0,iVariablePair);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,iVariableSingle);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,iVariableSingle);

  TH1D *gHistBalanceFunctionHistogram = 0x0;
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {
    gHistBalanceFunctionHistogram = (TH1D*)hTemp1->Clone();
    gHistBalanceFunctionHistogram->Reset();
    
    switch(iVariablePair) {
    case 1:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta #eta");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta #eta)");
      break;
    case 2:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta #phi (deg.)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta #phi)");
      break;
    default:
      break;
    }

    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-1.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-1.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5);
  }

  return gHistBalanceFunctionHistogram;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetBalanceFunctionDeltaEtaDeltaPhi(Double_t psiMin, 
							Double_t psiMax) {
  //Returns the BF histogram in Delta eta vs Delta phi, 
  //extracted from the 6 AliTHn objects 
  //(private members) of the AliBalancePsi class.
  //iVariableSingle: 0(phi-Psi), 1(pt-trigger)
  //iVariablePair: 0(phi-Psi) 1(Delta eta), 2(Delta phi), 3(pt-trigger), 4(pt-associated
  TString histName = "gHistBalanceFunctionHistogram2D";

  // Psi_2
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 

  Printf("P:%lf - N:%lf - PN:%lf - NP:%lf - PP:%lf - NN:%lf",fHistP->GetEntries(0),fHistN->GetEntries(0),fHistPN->GetEntries(0),fHistNP->GetEntries(0),fHistPP->GetEntries(0),fHistNN->GetEntries(0));

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH2D* hTemp1 = (TH2D*)fHistPN->Project(0,1,2);
  TH2D* hTemp2 = (TH2D*)fHistNP->Project(0,1,2);
  TH2D* hTemp3 = (TH2D*)fHistPP->Project(0,1,2);
  TH2D* hTemp4 = (TH2D*)fHistNN->Project(0,1,2);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,1);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,1);

  TH2D *gHistBalanceFunctionHistogram = 0x0;
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {
    gHistBalanceFunctionHistogram = (TH2D*)hTemp1->Clone();
    gHistBalanceFunctionHistogram->Reset();
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta #eta");   
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("#Delta #phi (deg.)");
    gHistBalanceFunctionHistogram->GetZaxis()->SetTitle("B(#Delta #eta,#Delta #phi)");   
    
    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-1.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-1.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5);
  }

  return gHistBalanceFunctionHistogram;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionPN(Double_t psiMin, 
					      Double_t psiMax) {
  //Returns the 2D correlation function for (+-) pairs
  // Psi_2: axis 0
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  //fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(-0.5,2.5); 
  //fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(-0.5,2.5); 

  TH2D *gHistTest = dynamic_cast<TH2D *>(fHistP->Project(0,0,1));
  TCanvas *c1 = new TCanvas("c1","");
  gHistTest->DrawCopy("colz");

  //0:step, 2: Delta eta, 3: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistPN->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  Printf("Entries (test): %lf",(Double_t)(gHistTest->GetEntries()));
  Printf("Entries (1D): %lf",(Double_t)(fHistP->Project(0,1)->GetEntries()));
  Printf("Entries (2D): %lf",(Double_t)(fHistPN->Project(0,1,2)->GetEntries()));
  
  TCanvas *c2 = new TCanvas("c2","");
  fHistPN->Project(0,1,2)->DrawCopy("colz");

  if((Double_t)(fHistP->Project(0,1)->GetEntries())!=0)
    gHist->Scale(1./(Double_t)(fHistP->Project(0,1)->GetEntries()));
    
  return gHist;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionNP(Double_t psiMin, 
					      Double_t psiMax) {
  //Returns the 2D correlation function for (+-) pairs
  // Psi_2: axis 0
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
    
  //0:step, 2: Delta eta, 3: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistNP->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  //Printf("Entries (1D): %lf",(Double_t)(fHistN->Project(0,2)->GetEntries()));
  //Printf("Entries (2D): %lf",(Double_t)(fHistNP->Project(0,2,3)->GetEntries()));
  if((Double_t)(fHistN->Project(0,1)->GetEntries())!=0)
    gHist->Scale(1./(Double_t)(fHistN->Project(0,1)->GetEntries()));
    
  return gHist;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionPP(Double_t psiMin, 
					      Double_t psiMax) {
  //Returns the 2D correlation function for (+-) pairs
  // Psi_2: axis 0
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
      
  //0:step, 2: Delta eta, 3: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistPP->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  //Printf("Entries (1D): %lf",(Double_t)(fHistP->Project(0,2)->GetEntries()));
  //Printf("Entries (2D): %lf",(Double_t)(fHistPP->Project(0,2,3)->GetEntries()));
  if((Double_t)(fHistP->Project(0,1)->GetEntries())!=0)
    gHist->Scale(1./(Double_t)(fHistP->Project(0,1)->GetEntries()));
  
  return gHist;
}

//____________________________________________________________________//
TH2D *AliBalancePsi::GetCorrelationFunctionNN(Double_t psiMin, 
					      Double_t psiMax) {
  //Returns the 2D correlation function for (+-) pairs
  // Psi_2: axis 0
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(psiMin,psiMax); 
    
  //0:step, 2: Delta eta, 3: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistNN->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  //Printf("Entries (1D): %lf",(Double_t)(fHistN->Project(0,2)->GetEntries()));
  //Printf("Entries (2D): %lf",(Double_t)(fHistNN->Project(0,2,3)->GetEntries()));
  if((Double_t)(fHistN->Project(0,1)->GetEntries())!=0)
    gHist->Scale(1./(Double_t)(fHistN->Project(0,1)->GetEntries()));
    
  return gHist;
}

