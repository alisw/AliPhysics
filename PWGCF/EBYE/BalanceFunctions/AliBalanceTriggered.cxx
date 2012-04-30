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
  bShuffle(kFALSE),
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
  TObject(balance), bShuffle(balance.bShuffle), 
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
void AliBalanceTriggered::FillBalance(Float_t fCentrality,vector<Double_t> **chargeVector) {
  // Calculates the balance function

 
  // Initialize histograms if not done yet
  if(!fHistPN){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("This works only in local mode --> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    InitHistograms();
  }

  Int_t gNtrack = chargeVector[0]->size();
  Double_t trackVarsSingle[nTrackVarsSingle];
  Double_t trackVarsPair[nTrackVarsPair];

  // 1st particle loop
  for(Int_t i = 0; i < gNtrack;i++){

    Short_t  charge = (Short_t) chargeVector[0]->at(i);
    trackVarsSingle[0]    =  chargeVector[2]->at(i);  //eta
    trackVarsSingle[1]    =  chargeVector[3]->at(i);  //phi
    trackVarsSingle[2]    =  chargeVector[7]->at(i);  //pt trigger
    trackVarsSingle[3]    =  fCentrality;             //centrality (really as variable here????)

    //fill single particle histograms
    if(charge > 0)  fHistP->Fill(trackVarsSingle,0,1.); 
    else            fHistN->Fill(trackVarsSingle,0,1.); 

    // 2nd particle loop
    for(Int_t j = 0; j < i; j++) {

      // need check for single particle region!!!???
      //
      
      Short_t charge2 = (Short_t) chargeVector[0]->at(j);
      trackVarsPair[0]    =  chargeVector[2]->at(i) - chargeVector[2]->at(j) ;  //delta eta
      trackVarsPair[1]    =  chargeVector[3]->at(i) - chargeVector[3]->at(j);  //delta phi
      trackVarsPair[2]    =  chargeVector[7]->at(j);  //pt
      trackVarsPair[3]    =  chargeVector[7]->at(i);  //pt trigger
      trackVarsPair[4]    =  fCentrality;             //centrality (really as variable here????)

    if( charge > 0 && charge2 < 0)  fHistPN->Fill(trackVarsPair,0,1.); 
    else if( charge < 0 && charge2 > 0)  fHistNP->Fill(trackVarsPair,0,1.); 
    else if( charge > 0 && charge2 > 0)  fHistPP->Fill(trackVarsPair,0,1.); 
    else if( charge < 0 && charge2 < 0)  fHistNN->Fill(trackVarsPair,0,1.); 
    else AliWarning("Wrong charge combination!");

    }//end of 2nd particle loop
  }//end of 1st particle loop
}  


TH1D* AliBalanceTriggered::GetBalanceFunctionHistogram1D(Int_t var, Double_t pTMin, Double_t pTMax, Double_t centrMin, Double_t centrMax){

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
  fHistP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMin,pTMax); 
  fHistN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMin,pTMax); 
  fHistPN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMin,pTMax); 
  fHistNP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMin,pTMax); 
  fHistPP->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMin,pTMax); 
  fHistNN->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMin,pTMax); 

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
  

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH1D* hTemp1 = (TH1D*)fHistPN->Project(0,var);
  TH1D* hTemp2 = (TH1D*)fHistNP->Project(0,var);
  TH1D* hTemp3 = (TH1D*)fHistPP->Project(0,var);
  TH1D* hTemp4 = (TH1D*)fHistNN->Project(0,var);
  TH1D* hTemp5 = (TH1D*)fHistP->Project(0,var);
  TH1D* hTemp6 = (TH1D*)fHistN->Project(0,var);

  TH1D* gHistBalanceFunctionHistogram = (TH1D*)hTemp1->Clone();
  gHistBalanceFunctionHistogram->Reset();
  
  // Calculate BF
  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)) {
    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-1.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-1.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5/1.);
  }

  return gHistBalanceFunctionHistogram;
}

TH1D* AliBalanceTriggered::GetHistogram1D(Int_t histo, Int_t var, Double_t pTMin, Double_t pTMax, Double_t centrMin, Double_t centrMax){

  // check which variable should be analyzed
  //
  // pair histograms:
  // 0 = Delta eta 
  // 1 = Delta phi
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

  default:
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
  gTHn->GetGrid(0)->GetGrid()->GetAxis(2)->SetRange(pTMin,pTMax); 
 
  // pT
  if(histo > 1) gTHn->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(pTMin,pTMax); 
 
  // centrality
  if(histo < 2) gTHn->GetGrid(0)->GetGrid()->GetAxis(3)->SetRange(centrMin,centrMax); 
  else          gTHn->GetGrid(0)->GetGrid()->GetAxis(4)->SetRange(centrMin,centrMax); 
 

  // Project into the wanted space (1st: analysis step, 2nd: axis)
  TH1D* gHisto = (TH1D*)gTHn->Project(0,var);

  return gHisto;
}
