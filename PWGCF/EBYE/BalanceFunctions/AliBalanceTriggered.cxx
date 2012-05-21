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

//-----------------------------------------------------------------
//           Balance Function class
//   This is the class to deal with the Balance Function analysis
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//   Modified: Michael Weber, m.weber@cern.ch
//-----------------------------------------------------------------


//ROOT
#include <Riostream.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>


#include <AliTHn.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TString.h>

#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include "AliBalanceTriggered.h"

ClassImp(AliBalanceTriggered)

//____________________________________________________________________//
AliBalanceTriggered::AliBalanceTriggered() :
  TObject(), 
  fAnalysisLevel("AOD"),
  fHistP(0x0),
  fHistN(0x0),
  fHistPN(0x0),
  fHistNP(0x0),
  fHistPP(0x0),
  fHistNN(0x0)
{
  // Default constructor
 
}


//____________________________________________________________________//
AliBalanceTriggered::AliBalanceTriggered(const AliBalanceTriggered& balance):
  TObject(balance), 
  fAnalysisLevel(balance.fAnalysisLevel),
  fHistP(balance.fHistP),
  fHistN(balance.fHistN),
  fHistPN(balance.fHistPN),
  fHistNP(balance.fHistNP),
  fHistPP(balance.fHistPP),
  fHistNN(balance.fHistNN)
{
  //copy constructor

}

//____________________________________________________________________//
AliBalanceTriggered::~AliBalanceTriggered() {
  // Destructor
 
  delete fHistP;
  delete fHistN;
  delete fHistPN;
  delete fHistNP;
  delete fHistPP;
  delete fHistNN;
  

}

//____________________________________________________________________//
void AliBalanceTriggered::InitHistograms() {

  //Initialize the histograms

  TString title    = "";      // histogram title
  Int_t anaSteps   = 1;       // analysis steps

  // single particle histograms
  Int_t iBinSingle[nTrackVarsSingle];         // binning for track variables
  Double_t* dBinsSingle[nTrackVarsSingle];    // bins for track variables  
  TString axisTitleSingle[nTrackVarsSingle];  // axis titles for track variables
  
  // two particle histograms
  Int_t iBinPair[nTrackVarsPair];         // binning for track variables
  Double_t* dBinsPair[nTrackVarsPair];    // bins for track variables  
  TString axisTitlePair[nTrackVarsPair];  // axis titles for track variables

   
  //-----------------------------------------------------------
  // histogram settings (hard coded!)
  //-----------------------------------------------------------

  // eta
  const Int_t kNEtaBins = 40;
  Double_t etaBins[kNEtaBins+1];
  for(Int_t i = 0; i < kNEtaBins+1; i++){
    etaBins[i] = -1.0 + i * 0.05;
  } 
  iBinSingle[0]       = kNEtaBins;
  dBinsSingle[0]      = etaBins;
  axisTitleSingle[0]  = "#eta"; 
  
  // delta eta
  const Int_t kNDeltaEtaBins = 40;
  Double_t deltaEtaBins[kNDeltaEtaBins+1];
  for(Int_t i = 0; i < kNDeltaEtaBins+1; i++){
    deltaEtaBins[i] = -2.0 + i * 0.1;
  } 
  iBinPair[0]       = kNDeltaEtaBins;
  dBinsPair[0]      = deltaEtaBins;
  axisTitlePair[0]  = "#Delta #eta"; 

  // phi
  const Int_t kNPhiBins = 72;
  Double_t phiBins[kNPhiBins+1];
  for(Int_t i = 0; i < kNPhiBins+1; i++){
    phiBins[i] = 0.0 + i * 5.;
  } 
  iBinSingle[1]       = kNPhiBins;
  dBinsSingle[1]      = phiBins;
  axisTitleSingle[1]  = "#phi (#circ)"; 
  
  // delta phi
  const Int_t kNDeltaPhiBins = 72;
  Double_t deltaPhiBins[kNDeltaPhiBins+1];
  for(Int_t i = 0; i < kNDeltaPhiBins+1; i++){
    deltaPhiBins[i] = -180.0 + i * 5.;
  } 
  iBinPair[1]       = kNDeltaPhiBins;
  dBinsPair[1]      = deltaPhiBins;
  axisTitlePair[1]  = "#Delta #phi (#circ)"; 

  // pT
  const Int_t kNPtBins = 22;
  Double_t pTBins[kNPtBins+1] = {0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0};
  iBinSingle[2]       = kNPtBins;
  dBinsSingle[2]      = pTBins;
  axisTitleSingle[2]  = "p_{T,trig} (GeV/c)";
 
  iBinPair[2]       = kNPtBins;
  dBinsPair[2]      = pTBins;
  axisTitlePair[2]  = "p_{T} (GeV/c)";
 
  iBinPair[3]       = kNPtBins;
  dBinsPair[3]      = pTBins;
  axisTitlePair[3]  = "p_{T,trig} (GeV/c)"; 

  // centrality
  const Int_t kNCentBins = 9;
  Double_t centBins[kNCentBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};

  iBinSingle[3]       = kNCentBins;
  dBinsSingle[3]      = centBins;
  axisTitleSingle[3]  = "centrality"; 

  iBinPair[4]       = kNCentBins;
  dBinsPair[4]      = centBins;
  axisTitlePair[4]  = "centrality";

  //-----------------------------------------------------------
  


  // User settings (depending on the analysis settings [only one at the moment])
  if(fAnalysisLevel == "AOD"){
    title = "hdEtaVsdPhi";
  }

  //-----------------------------------------------------------
  //-----------------------------------------------------------
  // Histogram creation

  // histogram for negative particles
  fHistN = new AliTHn(Form("fHistN"), Form("%s_N",title.Data()), anaSteps, nTrackVarsSingle, iBinSingle);
  for (Int_t j=0; j<nTrackVarsSingle; j++)
    {
      fHistN->SetBinLimits(j, dBinsSingle[j]);
      fHistN->SetVarTitle(j, axisTitleSingle[j]);
    }

  // histogram for positive particles
  fHistP = new AliTHn(Form("fHistP"), Form("%s_P",title.Data()), anaSteps, nTrackVarsSingle, iBinSingle);
  for (Int_t j=0; j<nTrackVarsSingle; j++)
    {
      fHistP->SetBinLimits(j, dBinsSingle[j]);
      fHistP->SetVarTitle(j, axisTitleSingle[j]);
    }

  // histogram for +- pairs
  fHistPN = new AliTHn(Form("fHistPN"), Form("%s_PN",title.Data()), anaSteps, nTrackVarsPair, iBinPair);
  for (Int_t j=0; j<nTrackVarsPair; j++)
    {
      fHistPN->SetBinLimits(j, dBinsPair[j]);
      fHistPN->SetVarTitle(j, axisTitlePair[j]);
    }

  // histogram for -+ pairs
  fHistNP = new AliTHn(Form("fHistNP"), Form("%s_NP",title.Data()), anaSteps, nTrackVarsPair, iBinPair);
  for (Int_t j=0; j<nTrackVarsPair; j++)
    {
      fHistNP->SetBinLimits(j, dBinsPair[j]);
      fHistNP->SetVarTitle(j, axisTitlePair[j]);
    }

  // histogram for ++ pairs
  fHistPP = new AliTHn(Form("fHistPP"), Form("%s_PP",title.Data()), anaSteps, nTrackVarsPair, iBinPair);
  for (Int_t j=0; j<nTrackVarsPair; j++)
    {
      fHistPP->SetBinLimits(j, dBinsPair[j]);
      fHistPP->SetVarTitle(j, axisTitlePair[j]);
    }

  // histogram for -- pairs
  fHistNN = new AliTHn(Form("fHistNN"), Form("%s_NN",title.Data()), anaSteps, nTrackVarsPair, iBinPair);
  for (Int_t j=0; j<nTrackVarsPair; j++)
    {
      fHistNN->SetBinLimits(j, dBinsPair[j]);
      fHistNN->SetVarTitle(j, axisTitlePair[j]);
    }

  //-----------------------------------------------------------
  //-----------------------------------------------------------


}

//____________________________________________________________________//
void AliBalanceTriggered::FillBalance(Float_t fCentrality,TObjArray *particles, TObjArray *particlesMixed) {
  // Calculates the balance function

 
  // Initialize histograms if not done yet
  if(!fHistPN){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("This works only in local mode --> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    InitHistograms();
  }

  Double_t trackVarsSingle[nTrackVarsSingle];
  Double_t trackVarsPair[nTrackVarsPair];

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

      
      Short_t  charge = (Short_t) firstParticle->Charge();
      trackVarsSingle[0]    =  firstEta;    //eta
      trackVarsSingle[1]    =  firstPhi;    //phi
      trackVarsSingle[2]    =  firstPt;     //pt trigger
      trackVarsSingle[3]    =  fCentrality; //centrality 
      
      //fill single particle histograms
      if(charge > 0)  fHistP->Fill(trackVarsSingle,0,1.); 
      else            fHistN->Fill(trackVarsSingle,0,1.);  

      
      // 2nd particle loop (only for j < i for non double counting in the same pT region)
      // --> SAME pT region for trigger and assoc: NO double counting with this
      // --> DIFF pT region for trigger and assoc: Missing assoc. particles with j > i to a trigger i 
      //                          --> can be handled afterwards by using assoc. as trigger as well ?!     
      for(Int_t j = 0; j < i; j++) {   // or go to full here (everything prepared)?
	
	if (particlesMixed && jMax < i)  // if the mixed track number is smaller than the main event one (could be done better if one loops over all tracks)
	  break;

	AliVParticle* secondParticle = (AliVParticle*) particlesSecond->At(j);

	Short_t charge2 = (Short_t) secondParticle->Charge();
	trackVarsPair[0]    =  firstEta - secondEta[j];  // delta eta
	trackVarsPair[1]    =  firstPhi - secondPhi[j];  // delta phi
	if (trackVarsPair[1] > 180)   // delta phi between -180 and 180 
	  trackVarsPair[1] -= 360;
	if (trackVarsPair[1] <  - 180) 
	  trackVarsPair[1] += 360;
	
	trackVarsPair[2]    =  firstPt;      // pt trigger
	trackVarsPair[3]    =  secondPt[j];  // pt
	trackVarsPair[4]    =  fCentrality;  // centrality
	
	if( charge > 0 && charge2 < 0)  fHistPN->Fill(trackVarsPair,0,1.); 
	else if( charge < 0 && charge2 > 0)  fHistNP->Fill(trackVarsPair,0,1.); 
	else if( charge > 0 && charge2 > 0)  fHistPP->Fill(trackVarsPair,0,1.); 
	else if( charge < 0 && charge2 < 0)  fHistNN->Fill(trackVarsPair,0,1.); 
	else AliWarning("Wrong charge combination!");
	
      }//end of 2nd particle loop
    }//end of 1st particle loop
}  


TH1D* AliBalanceTriggered::GetBalanceFunctionHistogram1D(Int_t var, Double_t pTMinTrigger, Double_t pTMaxTrigger, Double_t pTMin, Double_t pTMax, Double_t centrMin, Double_t centrMax, Double_t etaGap){

  // check which variable should be analyzed
  // 0 = Delta eta
  // 1 = Delta phi

  if( var < 0 || var > 1){
    AliError("Only Variable 0 (= Delta eta) or 1 (= Delta phi) allowed");
    return NULL;
  }


  // Choose region to analyze 
  // for Single Histograms (P,N):       2 = pT,trigger; 3 = centrality
  // for Pair Histograms (PN,NP,NN,PP): 2 = pT; 3 = pT,trigger; 4 = centrality

  // pT trigger
  fHistP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 

  // pT
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 

  // centrality
  fHistP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(centrMin,centrMax); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(centrMin,centrMax); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 

  if (etaGap > 0){
    fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(20 + etaGap, 9999); 
    fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(20 + etaGap, 9999); 
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(20 + etaGap, 9999); 
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(20 + etaGap, 9999); 
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(20 + etaGap, 9999); 
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(20 + etaGap, 9999); 
  }  

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH1D* hTemp1 = (TH1D*)fHistPN->Project(0,var);
  TH1D* hTemp2 = (TH1D*)fHistNP->Project(0,var);
  TH1D* hTemp3 = (TH1D*)fHistPP->Project(0,var);
  TH1D* hTemp4 = (TH1D*)fHistNN->Project(0,var);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,var);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,var);

  if (etaGap > 0){
    fHistP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(-9999, 9999); 
    fHistN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(-9999, 9999); 
    fHistPN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(-9999, 9999); 
    fHistNP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(-9999, 9999); 
    fHistPP->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(-9999, 9999); 
    fHistNN->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(-9999, 9999); 
  }

  TH1D* gHistBalanceFunctionHistogram = NULL;

  // Calculate BF
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {
    
    gHistBalanceFunctionHistogram = (TH1D*)hTemp1->Clone();
    gHistBalanceFunctionHistogram->Reset();

    // Calculate BF
    hTemp1->Add(hTemp3,-1.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-1.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5/1.);
  }

  return gHistBalanceFunctionHistogram;
}

TH2D* AliBalanceTriggered::GetBalanceFunctionHistogram2D(Int_t var1, Int_t var2, Double_t pTMinTrigger, Double_t pTMaxTrigger, Double_t pTMin, Double_t pTMax, Double_t centrMin, Double_t centrMax){

  // check which variable should be analyzed
  // 0 = Delta eta
  // 1 = Delta phi

  if( var1 < 0 || var1 > 1 || var2 < 0 || var2 > 1){
    AliError("Only Variable 0 (= Delta eta) or 1 (= Delta phi) allowed");
    return NULL;
  }


  // Choose region to analyze 
  // for Single Histograms (P,N):       2 = pT,trigger; 3 = centrality
  // for Pair Histograms (PN,NP,NN,PP): 2 = pT; 3 = pT,trigger; 4 = centrality

  // pT trigger
  fHistP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 

  // pT
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 

  // centrality
  fHistP->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(centrMin,centrMax); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(centrMin,centrMax); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
  

  // Project into the wanted space (1st: analysis step, 2nd: axis1, 3rd: axis2)
  TH2D* hTemp1 = (TH2D*)fHistPN->Project(0,var1,var2);
  TH2D* hTemp2 = (TH2D*)fHistNP->Project(0,var1,var2);
  TH2D* hTemp3 = (TH2D*)fHistPP->Project(0,var1,var2);
  TH2D* hTemp4 = (TH2D*)fHistNN->Project(0,var1,var2);
  TH2D* hTemp5 = (TH2D*)fHistP->Project(0,var1,var2);
  TH2D* hTemp6 = (TH2D*)fHistN->Project(0,var1,var2);

  TH2D* gHistBalanceFunctionHistogram = NULL;

  // Calculate BF
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)&&(hTemp5)&&(hTemp6)) {
    
    gHistBalanceFunctionHistogram = (TH2D*)hTemp1->Clone();
    gHistBalanceFunctionHistogram->Reset();
    
    // Calculate BF
    hTemp1->Add(hTemp3,-1.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-1.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
  }

  return gHistBalanceFunctionHistogram;
}


TH1D* AliBalanceTriggered::GetHistogram1D(Int_t histo, Int_t var, Double_t pTMinTrigger, Double_t pTMaxTrigger, Double_t pTMin, Double_t pTMax, Double_t centrMin, Double_t centrMax, Double_t etaGap){

  // check which variable should be analyzed
  //
  // single histograms:
  // 0 = eta 
  // 1 = phi
  // 2 = pT, trigger
  // 3 = centrality
  //
  // pair histograms:
  // 0 = Delta eta 
  // 1 = Delta phi
  // 2 = pT
  // 3 = pT, trigger
  // 4 = centrality

  if(histo < 0 || histo > 5){
    AliError("Only 6 histograms available: 0(P), 1(N), 2(PN), 3(NP), 4(PP), 5(NN)");
    return NULL;
  }

  if( histo > 1 && (var < 0 || var > 5)){
    AliError("Only Variable 0 to 4 allowed for pair histograms (histo > 1)");
    return NULL;
  }
  if( histo < 2 && (var < 0 || var > 4)){
    AliError("Only Variable 0 to 3 allowed for single histograms (histo < 2)");
    return NULL;
  }

  // get the histogram
  AliTHn *gTHn = NULL;
  switch(histo){
 
  case 0:
    gTHn = fHistP;
    break;

  case 1:
    gTHn = fHistN;
    break;

  case 2:
    gTHn = fHistPN;
    break;

  case 3:
    gTHn = fHistNP;
    break;

  case 4:
    gTHn = fHistPP;
    break;

  case 5:
    gTHn = fHistNN;
    break;
    
  }

  if(!gTHn){
    AliError(Form("AliTHn number %d = NULL",histo));
    return NULL;
  }

  // Choose region to analyze 
  // for Single Histograms (P,N):       2 = pT,trigger; 3 = centrality
  // for Pair Histograms (PN,NP,NN,PP): 2 = pT; 3 = pT,trigger; 4 = centrality

  // pT trigger
  gTHn->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
 
  // pT
  if(histo > 1) gTHn->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
 
  // centrality
  if(histo < 2) gTHn->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(centrMin,centrMax); 
  else          gTHn->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 

  if (etaGap > 0) gTHn->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(20 + etaGap, 9999); 

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH1D* gHisto      = (TH1D*)gTHn->Project(0,var);


  gTHn->GetGrid(0)->GetGrid()->GetAxis(0)->SetRange(-9999, 9999); 


  return gHisto;
}


TH2D* AliBalanceTriggered::GetHistogram2D(Int_t histo, Int_t var1, Int_t var2, Double_t pTMinTrigger, Double_t pTMaxTrigger, Double_t pTMin, Double_t pTMax, Double_t centrMin, Double_t centrMax){

  // check which variable should be analyzed
  //
  // single histograms:
  // 0 = eta 
  // 1 = phi
  // 2 = pT, trigger
  // 3 = centrality
  //
  // pair histograms:
  // 0 = Delta eta 
  // 1 = Delta phi
  // 2 = pT
  // 3 = pT, trigger
  // 4 = centrality

  if(histo < 0 || histo > 5){
    AliError("Only 6 histograms available: 0(P), 1(N), 2(PN), 3(NP), 4(PP), 5(NN)");
    return NULL;
  }

  if( histo > 1 && (var1 < 0 || var1 > 5 || var2 < 0 || var2 > 5)){
    AliError("Only Variable 0 to 4 allowed for pair histograms (histo > 1)");
    return NULL;
  }
  if( histo < 2 && (var1 < 0 || var1 > 4 || var2 < 0 || var2 > 4)){
    AliError("Only Variable 0 to 3 allowed for single histograms (histo < 2)");
    return NULL;
  }

  // get the histogram
  AliTHn *gTHn = NULL;
  switch(histo){
 
  case 0:
    gTHn = fHistP;
    break;

  case 1:
    gTHn = fHistN;
    break;

  case 2:
    gTHn = fHistPN;
    break;

  case 3:
    gTHn = fHistNP;
    break;

  case 4:
    gTHn = fHistPP;
    break;

  case 5:
    gTHn = fHistNN;
    break;
    
  }

  if(!gTHn){
    AliError(Form("AliTHn number %d = NULL",histo));
    return NULL;
  }

  // Choose region to analyze 
  // for Single Histograms (P,N):       2 = pT,trigger; 3 = centrality
  // for Pair Histograms (PN,NP,NN,PP): 2 = pT; 3 = pT,trigger; 4 = centrality

  // pT trigger
  gTHn->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMinTrigger,pTMaxTrigger); 
 
  // pT
  if(histo > 1) gTHn->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
 
  // centrality
  if(histo < 2) gTHn->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(centrMin,centrMax); 
  else          gTHn->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
 

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH2D* gHisto = (TH2D*)gTHn->Project(0,var1,var2);

  return gHisto;
}
