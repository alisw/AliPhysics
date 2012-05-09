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
  const Int_t kNCentralityBins = 9;
  Double_t centralityBins[kNCentralityBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};
  //const Int_t kNCentralityBins = 200;
  //Double_t centralityBins[kNCentralityBins+1];
  //for(Int_t i = 0; i < kNCentralityBins+1; i++)
  //centralityBins[i] = 0.0 + i * 0.5;
  iBinSingle[0]       = kNCentralityBins;
  dBinsSingle[0]      = centralityBins;
  axisTitleSingle[0]  = "Centrality percentile [%]"; 
  iBinPair[0]       = kNCentralityBins;
  dBinsPair[0]      = centralityBins;
  axisTitlePair[0]  = "Centrality percentile [%]"; 

  //Psi_2
  const Int_t kNPsi2Bins = 48;
  Double_t psi2Bins[kNPsi2Bins+1];
  for(Int_t i = 0; i < kNPsi2Bins+1; i++)
    psi2Bins[i] = -15.0 + i * 7.5;
  iBinSingle[1]       = kNPsi2Bins;
  dBinsSingle[1]      = psi2Bins;
  axisTitleSingle[1]  = "#phi - #Psi_{2} (#circ)";
  iBinPair[1]       = kNPsi2Bins;
  dBinsPair[1]      = psi2Bins;
  axisTitlePair[1]  = "#phi - #Psi_{2} (#circ)"; 
  
  // eta
  const Int_t kNEtaBins = 20;
  Double_t etaBins[kNEtaBins+1];
  for(Int_t i = 0; i < kNEtaBins+1; i++)
    etaBins[i] = -1.0 + i * 0.1;
  iBinSingle[2]       = kNEtaBins;
  dBinsSingle[2]      = etaBins;
  axisTitleSingle[2]  = "#eta"; 
  
  //#eta of triggered particle
  /*iBinPair[2]       = kNEtaBins;
  dBinsPair[2]      = etaBins;
  axisTitlePair[2]  = "#eta"; */

  // delta eta
  const Int_t kNDeltaEtaBins = 20;
  Double_t deltaEtaBins[kNDeltaEtaBins+1];
  for(Int_t i = 0; i < kNDeltaEtaBins+1; i++)
    deltaEtaBins[i] = -2.0 + i * 0.2;
  iBinPair[2]       = kNDeltaEtaBins;
  dBinsPair[2]      = deltaEtaBins;
  axisTitlePair[2]  = "#Delta #eta"; 

  // y
  /*const Int_t kNYBins = 40;
  Double_t yBins[kNYBins+1];
  for(Int_t i = 0; i < kNYBins+1; i++)
    yBins[i] = -1.0 + i * 0.05;
  iBinSingle[3]       = kNYBins;
  dBinsSingle[3]      = yBins;
  axisTitleSingle[3]  = "y"; 
  
  //y of triggered particle
  iBinPair[3]       = kNYBins;
  dBinsPair[3]      = yBins;
  axisTitlePair[3]  = "y"; */

  // delta y
  /*const Int_t kNDeltaYBins = 40;
  Double_t deltaYBins[kNDeltaYBins+1];
  for(Int_t i = 0; i < kNDeltaYBins+1; i++)
    deltaYBins[i] = -2.0 + i * 0.1;
  iBinPair[8]       = kNDeltaYBins;
  dBinsPair[8]      = deltaYBins;
  axisTitlePair[8]  = "#Delta y"; */

  // phi
  const Int_t kNPhiBins = 36;
  Double_t phiBins[kNPhiBins+1];
  for(Int_t i = 0; i < kNPhiBins+1; i++){
    phiBins[i] = 0.0 + i * 10.;
  } 
  iBinSingle[3]       = kNPhiBins;
  dBinsSingle[3]      = phiBins;
  axisTitleSingle[3]  = "#phi (#circ)"; 

  /*iBinPair[4]       = kNPhiBins;
  dBinsPair[4]      = phiBins;
  axisTitlePair[4]  = "#phi (#circ)"; */
  
  // pt(trigger)
  /*const Int_t kNPtBins = 100;
  Double_t ptBins[kNPtBins+1];
  for(Int_t i = 0; i < kNPtBins+1; i++){
    ptBins[i] = 0.0 + i * 0.2;
   } 
  iBinSingle[4]       = kNPtBins;
  dBinsSingle[4]      = ptBins;
  axisTitleSingle[4]  = "p_{t}^{trig.} (GeV/c)"; */

  /*iBinPair[5]       = kNPtBins;
  dBinsPair[5]      = ptBins;
  axisTitlePair[5]  = "p_{t}^{trig.} (GeV/c)"; 

  iBinPair[6]       = kNPtBins;
  dBinsPair[6]      = ptBins;
  axisTitlePair[6]  = "p_{t}^{assoc.} (GeV/c)"; */

  // delta phi
  const Int_t kNDeltaPhiBins = 72;
  Double_t deltaPhiBins[kNDeltaPhiBins+1];
  for(Int_t i = 0; i < kNDeltaPhiBins+1; i++){
    deltaPhiBins[i] = -180.0 + i * 5.;
  } 
  iBinPair[3]       = kNDeltaPhiBins;
  dBinsPair[3]      = deltaPhiBins;
  axisTitlePair[3]  = "#Delta #phi (#circ)"; 

  // qside
  /*const Int_t kNQSideBins = 200;
  Double_t qSideBins[kNQSideBins+1];
  for(Int_t i = 0; i < kNQSideBins+1; i++)
    qSideBins[i] = 0.0 + i * 0.02;
  iBinPair[10]       = kNQSideBins;
  dBinsPair[10]      = qSideBins;
  axisTitlePair[10]  = "q_{side} (GeV/c)";

  // qout
  const Int_t kNQoutBins = 200;
  Double_t qoutBins[kNQoutBins+1];
  for(Int_t i = 0; i < kNQoutBins+1; i++)
    qoutBins[i] = 0.0 + i * 0.02;
  iBinPair[11]       = kNQoutBins;
  dBinsPair[11]      = qoutBins;
  axisTitlePair[11]  = "q_{out} (GeV/c)";

  // qlong
  const Int_t kNQlongBins = 200;
  Double_t qlongBins[kNQlongBins+1];
  for(Int_t i = 0; i < kNQlongBins+1; i++)
    qlongBins[i] = 0.0 + i * 0.02;
  iBinPair[12]       = kNQlongBins;
  dBinsPair[12]      = qlongBins;
  axisTitlePair[12]  = "q_{long} (GeV/c)";

  // qinv
  const Int_t kNQinvBins = 200;
  Double_t qinvBins[kNQinvBins+1];
  for(Int_t i = 0; i < kNQinvBins+1; i++)
    qinvBins[i] = 0.0 + i * 0.02;
  iBinPair[13]       = kNQinvBins;
  dBinsPair[13]      = qinvBins;
  axisTitlePair[13]  = "q_{inv} (GeV/c)";*/

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
void AliBalancePsi::CalculateBalance(Float_t fCentrality, 
				     Double_t gReactionPlane, 
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
  
  for(i = 0; i < gNtrack; i++) {
    charge1 = chargeVector[0]->at(i);
    rap1    = chargeVector[1]->at(i);
    eta1    = chargeVector[2]->at(i);
    phi1    = chargeVector[3]->at(i);
    px1     = chargeVector[4]->at(i);
    py1     = chargeVector[5]->at(i);
    pz1     = chargeVector[6]->at(i);
    pt1     = chargeVector[7]->at(i);
    energy1 = chargeVector[8]->at(i);
    gPsiMinusPhi   = phi1 - gReactionPlane;
    if(gPsiMinusPhi < -fPsiInterval/2) gPsiMinusPhi = 360. + gPsiMinusPhi;
    
    trackVarsSingle[0]    =  fCentrality;
    trackVarsSingle[1]    =  gPsiMinusPhi;
    trackVarsSingle[2]    =  eta1;
    //trackVarsSingle[3]    =  rap1;
    trackVarsSingle[3]    =  phi1;  
    //trackVarsSingle[5]    =  pt1;  

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
      
      trackVarsPair[0]    =  fCentrality;             
      trackVarsPair[1]    =  gPsiMinusPhi;             
      //trackVarsPair[2]    =  eta1;
      //trackVarsPair[3]    =  rap1;
      //trackVarsPair[4]    =  phi1;
      //trackVarsPair[5]    =  pt1;
      //trackVarsPair[6]    =  pt2;
      trackVarsPair[2]    =  deta;
      //trackVarsPair[8]    =  dy;
      trackVarsPair[3]    =  dphi;
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
						 Double_t centrMin, 
						 Double_t centrMax, 
						 Double_t psiMin, 
						 Double_t psiMax) {
  //Returns the BF histogram, extracted from the 6 AliTHn objects 
  //(private members) of the AliBalancePsi class.
  //iVariableSingle: 2(eta), 3(y), 4(phi) 
  //iVariablePair: 2(Delta eta), 3(Delta y), 4(Delta phi), 5(qside), 6(qout), 7(qlong) 8(qinv) 
  TString gAnalysisType[ANALYSIS_TYPES] = {"y","eta","phi","qlong","qout","qside","qinv"};
  TString histName = "gHistBalanceFunctionHistogram";
  histName += gAnalysisType[iVariablePair];

  // centrality
  fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(centrMin,centrMax); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(centrMin,centrMax); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(centrMin,centrMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(centrMin,centrMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(centrMin,centrMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRangeUser(centrMin,centrMax); 

  // Psi_2
  fHistP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(psiMin,psiMax); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(psiMin,psiMax); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(psiMin,psiMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(psiMin,psiMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(psiMin,psiMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(1)->SetRangeUser(psiMin,psiMax); 

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
    case 2:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta #eta");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta #eta)");
      break;
    case 3:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta y");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta y)");
      break;
    case 4:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta #phi (deg.)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta #phi)");
      break;
    case 5:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("q_{side} (GeV/c)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(q_{side})");
      break;
    case 6:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("q_{out} (GeV/c)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(q_{out})");
      break;
    case 7:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("q_{long} (GeV/c)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(q_{long})");
      break;
    case 8:
      gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
      gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(q_{inv})");
      break;
    default:
      break;
    }

    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-2.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-2.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5);
  }

  return gHistBalanceFunctionHistogram;
}
