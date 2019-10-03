
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliHelperPID class
//-----------------------------------------------------------------

#include "AliHelperPID.h"
#include "AliVEvent.h"      
#include "AliAODEvent.h"      
#include "AliMCEvent.h"      
#include "AliMCEventHandler.h"      
#include "TH1F.h"             
#include "TH2F.h"             
#include "TList.h"            
#include "AliStack.h"            
#include "AliVTrack.h"      
#include "TParticle.h"
#include "AliAODMCParticle.h" 
#include "AliPIDResponse.h"   
#include "AliPIDCombined.h"   
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

using namespace AliHelperPIDNameSpace;
using namespace std;

const char * kPIDTypeName[]={"TPC","TOF","TPC-TOF"} ;
const char * kDetectorName[]={"ITS","TPC","TOF"} ;
const char * kParticleSpeciesName[]={"Pions","Kaons","Protons","Undefined"} ;

ClassImp(AliHelperPID)

AliHelperPID::AliHelperPID() : TNamed("HelperPID", "PID object"),fisMC(0), fPIDType(kNSigmaTPCTOF), fNSigmaPID(3), fBayesCut(0.8), fPIDResponse(0x0), fPIDCombined(0x0),fOutputList(0x0),fRequestTOFPID(1),fRemoveTracksT0Fill(0),fUseExclusiveNSigma(0),fPtTOFPID(.6),fHasTOFPID(0){

  // Fixing Leaks 
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  for(Int_t ipart=0;ipart<kNSpecies;ipart++)
    for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++)
      fnsigmas[ipart][ipid]=999.;
  
  for(Int_t ipart=0;ipart<kNSpecies;ipart++)fHasDoubleCounting[ipart]=kFALSE;
  
  fOutputList = new TList;
  fOutputList->SetOwner();
  fOutputList->SetName("OutputList");
  
  //nsigma plot
  for(Int_t ipart=0;ipart<kNSpecies;ipart++){
    for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
      Double_t miny=-30;
      Double_t maxy=30;
      if(ipid==kNSigmaTPCTOF){miny=0;maxy=50;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigma_%d_%d",ipart,ipid),Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //nsigmaRec plot
  for(Int_t ipart=0;ipart<kNSpecies;ipart++){
    for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
      Double_t miny=-10;
      Double_t maxy=10;
      if(ipid==kNSigmaTPCTOF){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaRec_%d_%d",ipart,ipid),
				  Form("n#sigma for reconstructed %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //BayesRec plot
  for(Int_t ipart=0;ipart<kNSpecies;ipart++){
    Double_t miny=0.;
    Double_t maxy=1;
    TH2F *fHistoBayes=new TH2F(Form("BayesRec_%d",ipart),
			       Form("probability for reconstructed %s",kParticleSpeciesName[ipart]),200,0,10,500,miny,maxy);
    fHistoBayes->GetXaxis()->SetTitle("P_{T} (GeV / c)");
    fHistoBayes->GetYaxis()->SetTitle(Form("Bayes prob %s",kParticleSpeciesName[ipart]));
    fOutputList->Add(fHistoBayes);
  }
  
  //nsigmaDC plot
  for(Int_t ipart=0;ipart<kNSpecies;ipart++){
    for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
      Double_t miny=-10;
      Double_t maxy=10;
      if(ipid==kNSigmaTPCTOF){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaDC_%d_%d",ipart,ipid),
				  Form("n#sigma for double counting %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //nsigmaMC plot
  for(Int_t ipart=0;ipart<kNSpecies;ipart++){
    for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
      Double_t miny=-30;
      Double_t maxy=30;
      if(ipid==kNSigmaTPCTOF){miny=0;maxy=50;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaMC_%d_%d",ipart,ipid),
				  Form("n#sigma for MC %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //PID signal plot
  for(Int_t idet=0;idet<kNDetectors;idet++){
    for(Int_t ipart=0;ipart<kNSpecies;ipart++){
      Double_t maxy=500;
      if(idet==kTOF)maxy=1.1;
      TH2F *fHistoPID=new TH2F(Form("PID_%d_%d",idet,ipart),Form("%s signal - %s",kDetectorName[idet],kParticleSpeciesName[ipart]),200,0,10,500,-maxy,maxy);
      fHistoPID->GetXaxis()->SetTitle("P (GeV / c)");
      fHistoPID->GetYaxis()->SetTitle(Form("%s signal",kDetectorName[idet]));
      fOutputList->Add(fHistoPID);
    }
  }
  //PID signal plot, before PID cut
  for(Int_t idet=0;idet<kNDetectors;idet++){
    Double_t maxy=500;
    if(idet==kTOF)maxy=1.1;
    TH2F *fHistoPID=new TH2F(Form("PIDAll_%d",idet),Form("%s signal",kDetectorName[idet]),200,0,10,500,-maxy,maxy);
    fHistoPID->GetXaxis()->SetTitle("P (GeV / c)");
    fHistoPID->GetYaxis()->SetTitle(Form("%s signal",kDetectorName[idet]));
    fOutputList->Add(fHistoPID);
  }
 TH1::AddDirectory(oldStatus);

}

//////////////////////////////////////////////////////////////////////////////////////////////////

TH2F* AliHelperPID::GetHistogram2D(const char * name){
  // returns histo named name
  return (TH2F*) fOutputList->FindObject(name);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliHelperPID::GetParticleSpecies(AliVTrack * trk, Bool_t FIllQAHistos){ 
  //return the specie according to the minimum nsigma value
  //no double counting, this has to be evaluated using CheckDoubleCounting()
  
  Int_t ID=kSpUndefined;
  
  //get the PID response
  if(!fPIDResponse) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  }
  if(!fPIDResponse) {
    AliFatal("Cannot get pid response");
  }
  
  //calculate nsigmas (used also by the bayesian)
  CalculateNSigmas(trk,FIllQAHistos);//fill the data member fnsigmas with the nsigmas value [ipart][iPID]
  
  //Do PID
  if(fPIDType==kBayes){//use bayesianPID
    
    if(!fPIDCombined) {
      // ------- setup PIDCombined
      fPIDCombined=new AliPIDCombined;
      fPIDCombined->SetDefaultTPCPriors();
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);  
    }
    if(!fPIDCombined) {
      AliFatal("PIDCombined object not found");
    }
    
    ID = GetIDBayes(trk,FIllQAHistos);
    
  }else{ //use nsigma PID
    
    ID=FindMinNSigma(trk,FIllQAHistos);
    
    if(fUseExclusiveNSigma){ //if one particle has double counting and exclusive nsigma is requested ID = kSpUndefined
      Bool_t *HasDC;
      HasDC=GetDoubleCounting(trk,FIllQAHistos);
      for(Int_t ipart=0;ipart<kNSpecies;ipart++){
	if(HasDC[ipart]==kTRUE)  ID = kSpUndefined;
      }
    }
  }

  if(FIllQAHistos){
    //Fill PID signal plot
    if(ID != kSpUndefined){
      for(Int_t idet=0;idet<kNDetectors;idet++){
	TH2F *h=GetHistogram2D(Form("PID_%d_%d",idet,ID));
	if(idet==kITS)h->Fill(trk->P(),trk->GetITSsignal()*trk->Charge());
	if(idet==kTPC)h->Fill(trk->P(),trk->GetTPCsignal()*trk->Charge());
	if(idet==kTOF && fHasTOFPID)h->Fill(trk->P(),TOFBetaCalc(trk)*trk->Charge());
      }
    }
    //Fill PID signal plot without cuts
    for(Int_t idet=0;idet<kNDetectors;idet++){
      TH2F *h=GetHistogram2D(Form("PIDAll_%d",idet));
      if(idet==kITS)h->Fill(trk->P(),trk->GetITSsignal()*trk->Charge());
      if(idet==kTPC)h->Fill(trk->P(),trk->GetTPCsignal()*trk->Charge());
      if(idet==kTOF && fHasTOFPID)h->Fill(trk->P(),TOFBetaCalc(trk)*trk->Charge());
    }
  }
  return ID;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliHelperPID::GetParticleSpecies(AliVParticle * part) {
  // return PID according to MC truth
  switch(TMath::Abs(part->PdgCode())){
  case 2212:
    return kSpProton;
    break;
  case 321:
    return kSpKaon;
    break;
  case 211:
    return kSpPion;
    break;
  default:
    return kSpUndefined;
  } 
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliHelperPID::GetIDBayes(AliVTrack * trk, Bool_t FIllQAHistos){ 
  
  Bool_t *IDs=GetAllCompatibleIdentitiesNSigma(trk,FIllQAHistos);
  
  Double_t probBayes[AliPID::kSPECIES];
  
  UInt_t detUsed= 0;
  CheckTOF(trk);
  if(fHasTOFPID && trk->Pt()>fPtTOFPID){//use TOF information
    detUsed = CalcPIDCombined(trk, fPIDResponse, AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC, probBayes);
    if(detUsed != (AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC))return kSpUndefined;//check that TPC and TOF are used
  }else{
    detUsed = CalcPIDCombined(trk, fPIDResponse,AliPIDResponse::kDetTPC, probBayes);
    if(detUsed != AliPIDResponse::kDetTPC)return kSpUndefined;//check that TPC is used
  }
  
  //the probability has to be normalized to one, we check it
  Double_t sump=0.;
  for(Int_t ipart=0;ipart<AliPID::kSPECIES;ipart++)sump+=probBayes[ipart];
  if(sump<.99 && sump>1.01){//FIXME precision problem in the sum, workaround
    AliFatal("Bayesian probability not normalized to one");
  }
  
  //probabilities are normalized to one, if the cut is above .5 there is no problem
  if(probBayes[AliPID::kPion]>fBayesCut && IDs[kSpPion]==1){
    TH2F *h=GetHistogram2D(Form("BayesRec_%d",kSpPion));
    h->Fill(trk->Pt(),probBayes[AliPID::kPion]);
    return kSpPion;
  }
  else if(probBayes[AliPID::kKaon]>fBayesCut && IDs[kSpKaon]==1){
    TH2F *h=GetHistogram2D(Form("BayesRec_%d",kSpKaon));
    h->Fill(trk->Pt(),probBayes[AliPID::kKaon]);
    return kSpKaon;
  }
  else if(probBayes[AliPID::kProton]>fBayesCut && IDs[kSpProton]==1){
    TH2F *h=GetHistogram2D(Form("BayesRec_%d",kSpProton));
    h->Fill(trk->Pt(),probBayes[AliPID::kProton]);
    return kSpProton;
  }
  else{
    return kSpUndefined;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

UInt_t AliHelperPID::CalcPIDCombined(const AliVTrack *track,const AliPIDResponse *PIDResponse, Int_t detMask, Double_t* prob) const{
  //
  // Bayesian PID calculation
  //
  for(Int_t i=0;i<AliPID::kSPECIES;i++)
    {
      prob[i]=0.;
    }
  fPIDCombined->SetDetectorMask(detMask);
  
  return fPIDCombined->ComputeProbabilities((AliVTrack*)track, PIDResponse, prob);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliHelperPID::CalculateNSigmas(AliVTrack * trk, Bool_t FIllQAHistos){ 
  //defines data member fnsigmas
  
  // Compute nsigma for each hypthesis
  AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(trk);
  // --- TPC
  Double_t nsigmaTPCkProton = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton);
  Double_t nsigmaTPCkKaon   = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon); 
  Double_t nsigmaTPCkPion   = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion); 
  // --- TOF
  Double_t nsigmaTOFkProton=999.,nsigmaTOFkKaon=999.,nsigmaTOFkPion=999.;
  Double_t nsigmaTPCTOFkProton=999.,nsigmaTPCTOFkKaon=999.,nsigmaTPCTOFkPion=999.;

  CheckTOF(trk);
  
  if(fHasTOFPID && trk->Pt()>fPtTOFPID){//use TOF information
    nsigmaTOFkProton = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton);
    nsigmaTOFkKaon   = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon); 
    nsigmaTOFkPion   = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion); 
    Double_t d2Proton=nsigmaTPCkProton * nsigmaTPCkProton + nsigmaTOFkProton * nsigmaTOFkProton;
    Double_t d2Kaon=nsigmaTPCkKaon * nsigmaTPCkKaon + nsigmaTOFkKaon * nsigmaTOFkKaon;
    Double_t d2Pion=nsigmaTPCkPion * nsigmaTPCkPion + nsigmaTOFkPion * nsigmaTOFkPion;
    //commented, this is done in analogy with AliESDTrackCuts, nsigma combind according to the probability
    // --- combined
    // -----------------------------------
    // How to get to a n-sigma cut?
    //
    // The accumulated statistics from 0 to d is
    //
    // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
    // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
    // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-d**2)/2)
    //
    // work around precision problem
    // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
    //if(TMath::Exp(- d2Proton / 2) > 1e-15)nsigmaTPCTOFkProton  =  TMath::Sqrt(2)*TMath::ErfInverse(1 - TMath::Exp(- d2Proton / 2));
    //if(TMath::Exp(- d2Kaon / 2) > 1e-15)nsigmaTPCTOFkKaon    = TMath::Sqrt(2)*TMath::ErfInverse(1 - TMath::Exp(- d2Kaon / 2));
    //if(TMath::Exp(- d2Pion / 2) > 1e-15)nsigmaTPCTOFkPion     = TMath::Sqrt(2)*TMath::ErfInverse(1 - TMath::Exp(- d2Pion / 2));
    
    //used for the 2PC PID paper (circular cut)
    nsigmaTPCTOFkProton  =  TMath::Sqrt(d2Proton);
    nsigmaTPCTOFkKaon    =  TMath::Sqrt(d2Kaon);
    nsigmaTPCTOFkPion    =  TMath::Sqrt(d2Pion);
  }else{
    // --- combined
    // if TOF is missing and below fPtTOFPID only the TPC information is used
    nsigmaTPCTOFkProton = TMath::Abs(nsigmaTPCkProton);
    nsigmaTPCTOFkKaon   = TMath::Abs(nsigmaTPCkKaon);
    nsigmaTPCTOFkPion   = TMath::Abs(nsigmaTPCkPion);
  }
  
  //set data member fnsigmas
  fnsigmas[kSpPion][kNSigmaTPC]=nsigmaTPCkPion;
  fnsigmas[kSpKaon][kNSigmaTPC]=nsigmaTPCkKaon;
  fnsigmas[kSpProton][kNSigmaTPC]=nsigmaTPCkProton;
  fnsigmas[kSpPion][kNSigmaTOF]=nsigmaTOFkPion;
  fnsigmas[kSpKaon][kNSigmaTOF]=nsigmaTOFkKaon;
  fnsigmas[kSpProton][kNSigmaTOF]=nsigmaTOFkProton;
  fnsigmas[kSpPion][kNSigmaTPCTOF]=nsigmaTPCTOFkPion;
  fnsigmas[kSpKaon][kNSigmaTPCTOF]=nsigmaTPCTOFkKaon;
  fnsigmas[kSpProton][kNSigmaTPCTOF]=nsigmaTPCTOFkProton;
  
  if(FIllQAHistos){
    //Fill NSigma SeparationPlot
    for(Int_t ipart=0;ipart<kNSpecies;ipart++){
      for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
	if((ipid!=kNSigmaTPC) && (!fHasTOFPID) && trk->Pt()<fPtTOFPID)continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigma_%d_%d",ipart,ipid));
	h->Fill(trk->Pt(),fnsigmas[ipart][ipid]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliHelperPID::FindMinNSigma(AliVTrack * trk,Bool_t FillQAHistos){ 
  
  CheckTOF(trk);  
  if(fRequestTOFPID && (!fHasTOFPID) && trk->Pt()>fPtTOFPID)return kSpUndefined;
  
  //get the identity of the particle with the minimum Nsigma
  Double_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  switch (fPIDType){
  case kNSigmaTPC:
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTPC]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTPC])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTPC])  ;
    break;
  case kNSigmaTOF:
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTOF])  ;
    break;
  case kNSigmaTPCTOF://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTPCTOF])  ;
    break;
  case kBayes://the nsigma in the bayesian is used to clean with a very large n-sigma value
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTPCTOF])  ;
    break;
  }
  
  // guess the particle based on the smaller nsigma (within fNSigmaPID)
  if( ( nsigmaKaon==nsigmaPion ) && ( nsigmaKaon==nsigmaProton )) return kSpUndefined;//if is the default value for the three
  
  if( ( nsigmaKaon   < nsigmaPion ) && ( nsigmaKaon < nsigmaProton ) && (nsigmaKaon   < fNSigmaPID)){
    if(FillQAHistos){
      for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
	if((ipid!=kNSigmaTPC) && (!fHasTOFPID) && (trk->Pt()<fPtTOFPID))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",kSpKaon,ipid));
	h->Fill(trk->Pt(),fnsigmas[kSpKaon][ipid]);
      }
    }
    return kSpKaon;
  }
  if( ( nsigmaPion   < nsigmaKaon ) && ( nsigmaPion < nsigmaProton ) && (nsigmaPion   < fNSigmaPID)){
    if(FillQAHistos){
      for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
	if((ipid!=kNSigmaTPC) && (!fHasTOFPID) && (trk->Pt()<fPtTOFPID))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",kSpPion,ipid));
	h->Fill(trk->Pt(),fnsigmas[kSpPion][ipid]);
      }
    }
    return kSpPion;
  }
  if( ( nsigmaProton < nsigmaKaon ) && ( nsigmaProton < nsigmaPion ) && (nsigmaProton < fNSigmaPID)){
    if(FillQAHistos){
      for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
	if((ipid!=kNSigmaTPC) && (!fHasTOFPID) && (trk->Pt()<fPtTOFPID))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",kSpProton,ipid));
	h->Fill(trk->Pt(),fnsigmas[kSpProton][ipid]);
      }
    }
    return kSpProton;
  }
  // else, return undefined
  return kSpUndefined;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t* AliHelperPID::GetDoubleCounting(AliVTrack * trk,Bool_t FIllQAHistos){ 
  //if a particle has double counting set fHasDoubleCounting[ipart]=kTRUE
  //fill DC histos
  for(Int_t ipart=0;ipart<kNSpecies;ipart++)fHasDoubleCounting[ipart]=kFALSE;//array with kTRUE for second (or third) identity of the track
  
  Int_t MinNSigma=FindMinNSigma(trk,kFALSE);//not filling the NSigmaRec histos
  
  CheckTOF(trk);
  
  if(MinNSigma==kSpUndefined)return fHasDoubleCounting;//in case of undefined no Double counting
  
  Double_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  switch (fPIDType) {
  case kNSigmaTPC:
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTPC]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTPC])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTPC])  ;
    break;
  case kNSigmaTOF:
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTOF])  ;
    break;
  case kNSigmaTPCTOF://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTPCTOF])  ;
    break;
  case kBayes://the nsigma in the bayesian is used to clean with a very large n-sigma value
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTPCTOF])  ;
    break;
  }
  if(nsigmaPion<fNSigmaPID && MinNSigma!=kSpPion)fHasDoubleCounting[kSpPion]=kTRUE;
  if(nsigmaKaon<fNSigmaPID && MinNSigma!=kSpKaon)fHasDoubleCounting[kSpKaon]=kTRUE;
  if(nsigmaProton<fNSigmaPID && MinNSigma!=kSpProton)fHasDoubleCounting[kSpProton]=kTRUE;
  
  if(FIllQAHistos){
    //fill NSigma distr for double counting
    for(Int_t ipart=0;ipart<kNSpecies;ipart++){
      if(fHasDoubleCounting[ipart]){
	for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
	  if((ipid!=kNSigmaTPC) && (!fHasTOFPID) && (trk->Pt()<fPtTOFPID))continue;//not filling TOF and combined if no TOF PID
	  TH2F *h=GetHistogram2D(Form("NSigmaDC_%d_%d",ipart,ipid));
	  h->Fill(trk->Pt(),fnsigmas[ipart][ipid]);
	}
      }
    }
  }
  return fHasDoubleCounting;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t* AliHelperPID::GetAllCompatibleIdentitiesNSigma(AliVTrack * trk,Bool_t FIllQAHistos){ 
  
  Bool_t *IDs=GetDoubleCounting(trk,FIllQAHistos);
  IDs[FindMinNSigma(trk,FIllQAHistos)]=kTRUE;
  return IDs;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliHelperPID::GetMCParticleSpecie(AliVEvent* event, AliVTrack * trk, Bool_t FillQAHistos){ 
  //return the specie according to the MC truth
  CheckTOF(trk);
  
  if(!fisMC)AliFatal("Error: AliHelperPID::GetMCParticleSpecie called on data\n");
  
  Int_t code=999;
  Bool_t isAOD=kFALSE;
  Bool_t isESD=kFALSE;
  TString nameoftrack(event->ClassName());  
  if(!nameoftrack.CompareTo("AliESDEvent"))isESD=kTRUE;
  else if(!nameoftrack.CompareTo("AliAODEvent"))isAOD=kTRUE;
  else AliFatal("Not processing AODs nor ESDs") ;
  
  if(isAOD){
    TClonesArray *arrayMC = 0;
    arrayMC = (TClonesArray*) event->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!arrayMC)AliFatal("Error: MC particles branch not found!\n");
    AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(trk->GetLabel()));
    if (!partMC){
      AliError("Cannot get MC particle");
      return kSpUndefined;
    }
    code=partMC->GetPdgCode();
  }
  if(isESD){
    AliStack* stack=0x0;
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) AliFatal("ERROR: Could not retrieve MC event handler");
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent)AliFatal("ERROR: Could not retrieve MC event");
    stack = mcEvent->Stack();
    if (!stack) AliFatal("ERROR: stack not available\n");
    TParticle *part=0;
    part = (TParticle*)stack->Particle(TMath::Abs(trk->GetLabel()));
    if (!part){
      AliError("Cannot get MC particle");
      return kSpUndefined;
    }
    code = part->GetPdgCode();
  }
  switch(TMath::Abs(code)){
  case 2212:
    if(FillQAHistos){
      for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
	if((ipid!=kNSigmaTPC) && (!fHasTOFPID) && (trk->Pt()<fPtTOFPID))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",kSpProton,ipid));
	h->Fill(trk->Pt(),fnsigmas[kSpProton][ipid]);
      }
    }
    return kSpProton;
    break;
  case 321:
    if(FillQAHistos){
      for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
	if((ipid!=kNSigmaTPC) && (!fHasTOFPID) && (trk->Pt()<fPtTOFPID))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",kSpKaon,ipid));
	h->Fill(trk->Pt(),fnsigmas[kSpKaon][ipid]);
      }
    }
    return kSpKaon;
    break;
  case 211:
    if(FillQAHistos){
      for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
	if((ipid!=kNSigmaTPC) && (!fHasTOFPID) && (trk->Pt()<fPtTOFPID))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",kSpPion,ipid));
	h->Fill(trk->Pt(),fnsigmas[kSpPion][ipid]);
      }
    }
    return kSpPion;
    break;
  default:
    return kSpUndefined;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliHelperPID::CheckTOF(AliVTrack * trk)
{
  //check if the particle has TOF Matching
  
  //get the PIDResponse
  if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,trk)==0)fHasTOFPID=kFALSE;
  else fHasTOFPID=kTRUE;
  
  //in addition to TOF status we look at the pt
  if(trk->Pt()<fPtTOFPID)fHasTOFPID=kFALSE;
  
  if(fRemoveTracksT0Fill)
    {
      Int_t startTimeMask = fPIDResponse->GetTOFResponse().GetStartTimeMask(trk->P());
      if (startTimeMask < 0)fHasTOFPID=kFALSE; 
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Double_t AliHelperPID::TOFBetaCalc(AliVTrack *track) const{
  //TOF beta calculation
  Double_t tofTime=track->GetTOFsignal();
  
  Double_t c=TMath::C()*1.E-9;// m/ns
  Float_t startTime = fPIDResponse->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P());//in ps
  Double_t length= fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c;
  tofTime -= startTime;      // subtract startTime to the signal
  Double_t tof= tofTime*1E-3; // ns, average T0 fill subtracted, no info from T0detector 	 
  tof=tof*c;
  return length/tof;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Double_t AliHelperPID::GetMass(AliHelperParticleSpecies_t id) const{
  //return Mass according to AliHelperParticleSpecies_t. If undefined return -999.
  Double_t mass=-999.;
  if (id == kSpProton) { mass=9.38271999999999995e-01; }
  if (id == kSpKaon)   { mass=4.93676999999999977e-01; }
  if (id == kSpPion)    { mass=1.39570000000000000e-01; }
  return mass;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Long64_t AliHelperPID::Merge(TCollection* list)
{
  // Merging interface.
  // Returns the number of merged objects (including this).
  Printf("Merging");
  if (!list)
    return 0;
  if (list->IsEmpty())
    return 1;
  TIterator* iter = list->MakeIterator();
  TObject* obj;
  // collections of TList with output histos
  TList collections;
  Int_t count = 0;
  while ((obj = iter->Next())) {
    AliHelperPID* entry = dynamic_cast<AliHelperPID*> (obj);
    if (entry == 0) 
      continue;
    TList * outputlist = entry->GetOutputList();      
    collections.Add(outputlist);
    count++;
  }
  fOutputList->Merge(&collections);
  delete iter;
  Printf("OK");
  return count+1;
}
