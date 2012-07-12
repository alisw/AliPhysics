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
  AliInfo("Finished setting up the AliTHn");
}

//____________________________________________________________________//
void AliBalancePsi::CalculateBalance(Double_t gReactionPlane,
				     TObjArray *particles, 
				     TObjArray *particlesMixed ) {
  // Calculates the balance function
  fAnalyzedEvents++;
  
  // Initialize histograms if not done yet
  if(!fHistPN){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("This works only in local mode --> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    InitHistograms();
  }

  Double_t trackVariablesSingle[nTrackVariablesSingle];
  Double_t trackVariablesPair[nTrackVariablesPair];

  if (!particles){
    AliWarning("particles TObjArray is NULL pointer --> return");
    return;
  }
  
  // define end of particle loops
  Int_t iMax = particles->GetEntriesFast();
  Int_t jMax = iMax;
  if (particlesMixed)
    jMax = particlesMixed->GetEntriesFast();

  // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TObjArray* particlesSecond = (particlesMixed) ? particlesMixed : particles;

  TArrayF secondEta(jMax);
  TArrayF secondPhi(jMax);
  TArrayF secondPt(jMax);

  for (Int_t i=0; i<jMax; i++){
    secondEta[i] = ((AliVParticle*) particlesSecond->At(i))->Eta();
    secondPhi[i] = ((AliVParticle*) particlesSecond->At(i))->Phi();
    secondPt[i]  = ((AliVParticle*) particlesSecond->At(i))->Pt();
  }
  
  // 1st particle loop
  for (Int_t i=0; i<iMax; i++)
    {
      
      AliVParticle* firstParticle = (AliVParticle*) particles->At(i);
      
      // some optimization
      Float_t firstEta = firstParticle->Eta();
      Float_t firstPhi = firstParticle->Phi();
      Float_t firstPt  = firstParticle->Pt();

      // Event plane (determine psi bin)
      Double_t gPsiMinusPhi    =   0.;
      Double_t gPsiMinusPhiBin = -10.;
      gPsiMinusPhi   = TMath::Abs(firstPhi - gReactionPlane);
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

      
      Short_t  charge = (Short_t) firstParticle->Charge();

      trackVariablesSingle[0]    =  gPsiMinusPhiBin;
      trackVariablesSingle[1]    =  firstPt;  
      
      //fill single particle histograms
      if(charge > 0)      fHistP->Fill(trackVariablesSingle,0,1.); 
      else if(charge < 0) fHistN->Fill(trackVariablesSingle,0,1.);  
      

      
      // 2nd particle loop (only for j < i for non double counting in the same pT region)
      // --> SAME pT region for trigger and assoc: NO double counting with this
      // --> DIFF pT region for trigger and assoc: Missing assoc. particles with j > i to a trigger i 
      //                          --> can be handled afterwards by using assoc. as trigger as well ?!     
      for(Int_t j = 0; j < i; j++) {   // or go to full here (everything prepared)?
	
	if (particlesMixed && jMax < i)  // if the mixed track number is smaller than the main event one (could be done better if one loops over all tracks)
	  break;

	AliVParticle* secondParticle = (AliVParticle*) particlesSecond->At(j);

	Short_t charge2 = (Short_t) secondParticle->Charge();
	
	trackVariablesPair[0]    =  gPsiMinusPhiBin;
	trackVariablesPair[1]    =  firstEta - secondEta[j];  // delta eta
	trackVariablesPair[2]    =  firstPhi - secondPhi[j];  // delta phi
	if (trackVariablesPair[2] > 180.)   // delta phi between -180 and 180 
	  trackVariablesPair[2] -= 360.;
	if (trackVariablesPair[2] <  - 180.) 
	  trackVariablesPair[2] += 360.;
	
	trackVariablesPair[3]    =  firstPt;      // pt trigger
	trackVariablesPair[4]    =  secondPt[j];  // pt
	//	trackVariablesPair[5]    =  fCentrality;  // centrality
	
	if( charge > 0 && charge2 < 0)  fHistPN->Fill(trackVariablesPair,0,1.); 
	else if( charge < 0 && charge2 > 0)  fHistNP->Fill(trackVariablesPair,0,1.); 
	else if( charge > 0 && charge2 > 0)  fHistPP->Fill(trackVariablesPair,0,1.); 
	else if( charge < 0 && charge2 < 0)  fHistNN->Fill(trackVariablesPair,0,1.); 
	else AliWarning(Form("Wrong charge combination: charge1 = %d and charge2 = %d",charge,charge2));
	
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

  AliInfo(Form("P:%lf - N:%lf - PN:%lf - NP:%lf - PP:%lf - NN:%lf",fHistP->GetEntries(0),fHistN->GetEntries(0),fHistPN->GetEntries(0),fHistNP->GetEntries(0),fHistPP->GetEntries(0),fHistNN->GetEntries(0)));

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
  c1->cd();
  if(!gHistTest){
    AliError("Projection of fHistP = NULL");
    return gHistTest;
  }
  else{
    gHistTest->DrawCopy("colz");
  }

  //0:step, 2: Delta eta, 3: Delta phi
  TH2D *gHist = dynamic_cast<TH2D *>(fHistPN->Project(0,1,2));
  if(!gHist){
    AliError("Projection of fHistPN = NULL");
    return gHist;
  }

  AliInfo(Form("Entries (test): %lf",(Double_t)(gHistTest->GetEntries())));
  AliInfo(Form("Entries (1D): %lf",(Double_t)(fHistP->Project(0,1)->GetEntries())));
  AliInfo(Form("Entries (2D): %lf",(Double_t)(fHistPN->Project(0,1,2)->GetEntries())));
  
  TCanvas *c2 = new TCanvas("c2","");
  c2->cd();
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

