
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
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

using namespace AliHelperPIDNameSpace;
using namespace std;

const char * kPIDTypeName[]={"TPC","TOF","TPC-TOF"} ;
const char * kDetectorName[]={"TPC","TOF"} ;
const char * kParticleSpeciesName[]={"Pions","Kaons","Protons","Undefined"} ;

ClassImp(AliHelperPID)

AliHelperPID::AliHelperPID() : TNamed("HelperPID", "PID object"),fisMC(0), fPIDType(kNSigmaTPCTOF), fNSigmaPID(3), fPIDResponse(0),fOutputList(0),fRequestTOFPID(1),fPtTOFPID(.6),fHasTOFPID(0){
  
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
      TH2F *fHistoNSigma=new TH2F(Form("NSigma_%d_%d",ipart,ipid),Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,200,miny,maxy);
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
				  Form("n#sigma for reconstructed %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,100,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //nsigmaDC plot
  for(Int_t ipart=0;ipart<kNSpecies;ipart++){
    for(Int_t ipid=0;ipid<=kNSigmaPIDType;ipid++){
      Double_t miny=-10;
      Double_t maxy=10;
      if(ipid==kNSigmaTPCTOF){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaDC_%d_%d",ipart,ipid),
				  Form("n#sigma for double counting %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,100,miny,maxy);
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
				  Form("n#sigma for MC %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,100,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //PID signal plot
  for(Int_t idet=0;idet<kNDetectors;idet++){
    Double_t maxy=1000;
    if(idet==kTOF)maxy=200000;
    TH2F *fHistoPID=new TH2F(Form("PID_%d",idet),Form("%s signal",kDetectorName[idet]),200,0,10,200,-maxy,maxy);
    fHistoPID->GetXaxis()->SetTitle("P (GeV / c)");
    fHistoPID->GetYaxis()->SetTitle(Form("%s signal",kDetectorName[idet]));
    fOutputList->Add(fHistoPID);
  }
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
  //Printf("fPtTOFPID %.1f, fRequestTOFPID %d, fNSigmaPID %.1f, fPIDType %d",fPtTOFPID,fRequestTOFPID,fNSigmaPID,fPIDType);
  
  CalculateNSigmas(trk,FIllQAHistos);//fill the data member fnsigmas with the nsigmas value [ipart][iPID]
  
  return FindMinNSigma(trk,FIllQAHistos);//NSigmaRec distr filled here
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliHelperPID::CalculateNSigmas(AliVTrack * trk, Bool_t FIllQAHistos){ 
  //defines data member fnsigmas
  if(!fPIDResponse) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  }
  if(!fPIDResponse) {
    AliFatal("Cannot get pid response");
  }
  
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
    // --- combined
    nsigmaTPCTOFkProton = TMath::Sqrt((nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton)/2.);
    nsigmaTPCTOFkKaon   = TMath::Sqrt((nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon)/2.);
    nsigmaTPCTOFkPion   = TMath::Sqrt((nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion)/2.);
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
    //Fill PID signal plot
    for(Int_t idet=0;idet<kNDetectors;idet++){
      TH2F *h=GetHistogram2D(Form("PID_%d",idet));
      if(idet==kTPC)h->Fill(trk->P(),trk->GetTPCsignal()*trk->Charge());
      if(idet==kTOF && fHasTOFPID)h->Fill(trk->P(),trk->GetTOFsignal()*trk->Charge());
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
  for(Int_t ipart=0;ipart<=kNSpecies;ipart++)fHasDoubleCounting[ipart]=kFALSE;//array with kTRUE for second (or third) identity of the track
  
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
  UInt_t status;
  status=trk->GetStatus();
  if((status&AliVTrack::kTOFout)==0 || (status&AliVTrack::kTIME)==0)fHasTOFPID=kFALSE;
  else fHasTOFPID=kTRUE;
  //in addition to KTOFout and kTIME we look at the pt
  if(trk->Pt()<fPtTOFPID)fHasTOFPID=kFALSE;
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
