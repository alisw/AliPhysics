/*
 ***********************************************************
 Manager for event plane corrections framework
Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
Instructions in AddTask_EPcorrectionsExample.C
2014/12/10
 *********************************************************
 */

//#include "AliSysInfo.h"
#include <iostream>

#include <TROOT.h>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <TChain.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliESDEvent.h>
#include <TList.h>
#include "AliQnCorrectionsVarManager.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliQnCorrectionsFillEvent.h"
#include "AliQnCorrectionsHistos.h"
#include "AliQnCorrectionsCuts.h"
#include "AliAnalysisTaskQnAnalysis.h"

// make a change

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskQnAnalysis)


//_________________________________________________________________________________
AliAnalysisTaskQnAnalysis::AliAnalysisTaskQnAnalysis() :
  AliAnalysisTaskSE(),
  fEventQAList(0x0),
  fEventCuts(0x0),
  fEventPlaneHistos(0x0),
  fFillEvent(0x0)
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskQnAnalysis::AliAnalysisTaskQnAnalysis(const char* name) :
  AliAnalysisTaskSE(name),
  fEventQAList(0x0),
  fEventCuts(0x0),
  fEventPlaneHistos(0x0),
  fFillEvent(0x0)
{
  //
  // Constructor
  //

  fFillEvent = new AliQnCorrectionsFillEvent();

  fEventQAList = new TList();
  fEventQAList->SetName("EventQA");
  fEventQAList->SetOwner(kTRUE);

  fEventPlaneHistos = new AliQnCorrectionsHistos();

  for(Int_t i=0; i<kNcorrelation; i++){
    fEPflow[i][0]=-1;
    fEPflow[i][1]=-1;
    for(Int_t j=0; j<kNharmonics; j++)
      for(Int_t k=0; k<4; k++){
        fVn[i][j][k]=0x0;
  }}
  for(Int_t i=0; i<kNcorrelation*3; i++){
    fEPcorrelation[i][0]=-1;
    fEPcorrelation[i][1]=-1;
    for(Int_t j=0; j<kNharmonics; j++)
      for(Int_t k=0; k<4; k++){
        fCorrelation[i][j][k]=0x0;
  }}
  


  DefineInput(0,TChain::Class());
  DefineInput(1,TList::Class());
  DefineOutput(1, TList::Class());// Event QA histograms
}



//____________________________________________________________________________
AliAnalysisTaskQnAnalysis::~AliAnalysisTaskQnAnalysis()
{
  //
  // Constructor
  //

  //fEventQAList->Clear("");
  //fEventPlaneHistos->Clear("");
  //delete fEventQAList;
  //delete fEventPlaneHistos;

  //for(Int_t i=0; i<kNcorrelation; i++){
  //  for(Int_t j=0; j<kNharmonics; j++)
  //    for(Int_t k=0; k<4; k++){
  //      delete fVn[i][j][k];
  //}}
  //for(Int_t i=0; i<kNcorrelation*3; i++){
  //  for(Int_t j=0; j<kNharmonics; j++)
  //    for(Int_t k=0; k<4; k++){
  //      delete fCorrelation[i][j][k];
  //}}
  
}



//_________________________________________________________________________________
void AliAnalysisTaskQnAnalysis::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  fEPflow[0][0] = kTPC;
  fEPflow[0][1] = kVZEROA;

  fEPflow[1][0] = kTPC;
  fEPflow[1][1] = kVZEROC;

  fEPflow[2][0] = kTPC;
  fEPflow[2][1] = kFMDA;

  fEPflow[3][0] = kTPC;
  fEPflow[3][1] = kFMDC;

  fEPcorrelation[0][0] = kTPC;
  fEPcorrelation[0][1] = kVZEROA;

  fEPcorrelation[1][0] = kTPC;
  fEPcorrelation[1][1] = kVZEROC;

  fEPcorrelation[2][0] = kVZEROA;
  fEPcorrelation[2][1] = kVZEROC;

  fEPcorrelation[3][0] = kTPC;
  fEPcorrelation[3][1] = kFMDA;

  fEPcorrelation[4][0] = kTPC;
  fEPcorrelation[4][1] = kFMDC;

  fEPcorrelation[5][0] = kFMDA;
  fEPcorrelation[5][1] = kFMDC;

  TString names[kNqvectors] = {"TPC","V0A","V0C","FMDA","FMDC"};
  TString components[4] = {"XX","XY","YX","YY"};

  Double_t centbinning[12] = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
  for(Int_t i=0; i<kNcorrelation; i++)
    for(Int_t j=0; j<kNharmonics; j++)
      for(Int_t k=0; k<4; k++){
        fVn[i][j][k]=0x0;
        if(fEPflow[i][0]==-1) continue;
        fVn[i][j][k] = new TProfile(Form("vn_%sx%s_%s_h%d",names[fEPflow[i][0]].Data(),names[fEPflow[i][1]].Data(), components[k].Data(),j+1),Form("vn_%sx%s_%s_h%d",names[fEPflow[i][0]].Data(),names[fEPflow[i][1]].Data(), components[k].Data(),j+1), 11, centbinning); 
  }

  for(Int_t i=0; i<kNcorrelation*3; i++)
    for(Int_t j=0; j<kNharmonics; j++)
      for(Int_t k=0; k<4; k++){
        fCorrelation[i][j][k]=0x0;
        if(fEPcorrelation[i][0]==-1) continue;
        fCorrelation[i][j][k] = new TProfile(Form("cor_%sx%s_%s_h%d",names[fEPcorrelation[i][0]].Data(),names[fEPcorrelation[i][1]].Data(), components[k].Data(),j+1),Form("vn_%sx%s_%s_h%d",names[fEPcorrelation[i][0]].Data(),names[fEPcorrelation[i][1]].Data(), components[k].Data(),j+1), 11, centbinning); 
  }

  PostData(1, fEventQAList);

}



//________________________________________________________________________________________________________
void AliAnalysisTaskQnAnalysis::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //

  AliVEvent* event = InputEvent();
  TList* qnlist = dynamic_cast<TList*>(GetInputData(1));
  if(!qnlist) return;

  Float_t values[2000]={-999.};
  fFillEvent->Process((AliAnalysisTaskSE*) this, event, values);

  fEventPlaneHistos->FillHistClass("Event_NoCuts", values);
  if(!IsEventSelected(values)) return;
  fEventPlaneHistos->FillHistClass("Event_Analysis", values);

  AliQnCorrectionsQnVector* qvec[kNqvectors] = {0x0};
  qvec[kFMDA]   = GetQvector(qnlist, "FMDA"  , AliQnCorrectionsConstants::kAlignment);
  qvec[kFMDC]   = GetQvector(qnlist, "FMDC"  , AliQnCorrectionsConstants::kAlignment);
  qvec[kVZEROA] = GetQvector(qnlist, "VZEROA", AliQnCorrectionsConstants::kAlignment);
  qvec[kVZEROC] = GetQvector(qnlist, "VZEROC", AliQnCorrectionsConstants::kAlignment);
  qvec[kTPC]    = GetQvector(qnlist, "TPC"   , AliQnCorrectionsConstants::kRecentering);

  for(Int_t i=0; i<kNcorrelation; i++)
    for(Int_t j=0; j<kNharmonics; j++){
      if(fEPflow[i][0]==-1) continue;
      if(!qvec[fEPflow[i][0]]) continue;
      if(qvec[fEPflow[i][0]]->CheckQnVectorStatus(j+1,AliQnCorrectionsConstants::kUndefined)) continue;
      if(qvec[fEPflow[i][1]]->CheckQnVectorStatus(j+1,AliQnCorrectionsConstants::kUndefined)) continue;
      if(qvec[fEPflow[i][0]]->N()==0||qvec[fEPflow[i][1]]->N()==0) continue;
        fVn[i][j][0]->Fill(values[AliQnCorrectionsVarManager::kCentVZERO],  qvec[fEPflow[i][0]]->Qx(j+1)*qvec[fEPflow[i][1]]->QxNorm(j+1));
        fVn[i][j][1]->Fill(values[AliQnCorrectionsVarManager::kCentVZERO],  qvec[fEPflow[i][0]]->Qx(j+1)*qvec[fEPflow[i][1]]->QyNorm(j+1));
        fVn[i][j][2]->Fill(values[AliQnCorrectionsVarManager::kCentVZERO],  qvec[fEPflow[i][0]]->Qy(j+1)*qvec[fEPflow[i][1]]->QxNorm(j+1));
        fVn[i][j][3]->Fill(values[AliQnCorrectionsVarManager::kCentVZERO],  qvec[fEPflow[i][0]]->Qy(j+1)*qvec[fEPflow[i][1]]->QyNorm(j+1));
  }

  for(Int_t i=0; i<kNcorrelation*3; i++)
    for(Int_t j=0; j<kNharmonics; j++){
      if(fEPcorrelation[i][0]==-1) continue;
      if(!qvec[fEPcorrelation[i][0]]) continue;
      if(qvec[fEPcorrelation[i][0]]->CheckQnVectorStatus(j+1,AliQnCorrectionsConstants::kUndefined)) continue;
      if(qvec[fEPcorrelation[i][1]]->CheckQnVectorStatus(j+1,AliQnCorrectionsConstants::kUndefined)) continue;
      if(qvec[fEPcorrelation[i][0]]->N()==0||qvec[fEPcorrelation[i][1]]->N()==0) continue;
        fCorrelation[i][j][0]->Fill(values[AliQnCorrectionsVarManager::kCentVZERO],  qvec[fEPcorrelation[i][0]]->QxNorm(j+1)*qvec[fEPcorrelation[i][1]]->QxNorm(j+1));
        fCorrelation[i][j][1]->Fill(values[AliQnCorrectionsVarManager::kCentVZERO],  qvec[fEPcorrelation[i][0]]->QxNorm(j+1)*qvec[fEPcorrelation[i][1]]->QyNorm(j+1));
        fCorrelation[i][j][2]->Fill(values[AliQnCorrectionsVarManager::kCentVZERO],  qvec[fEPcorrelation[i][0]]->QyNorm(j+1)*qvec[fEPcorrelation[i][1]]->QxNorm(j+1));
        fCorrelation[i][j][3]->Fill(values[AliQnCorrectionsVarManager::kCentVZERO],  qvec[fEPcorrelation[i][0]]->QyNorm(j+1)*qvec[fEPcorrelation[i][1]]->QyNorm(j+1));
  }



}  // end loop over events


//__________________________________________________________________
void AliAnalysisTaskQnAnalysis::FinishTaskOutput()
{
  //
  // Finish Task 
  //

  THashList* hList = (THashList*) fEventPlaneHistos->HistList();
  for(Int_t i=0; i<hList->GetEntries(); ++i) {
    THashList* list = (THashList*)hList->At(i);
    fEventQAList->Add(list);
  }

  for(Int_t i=0; i<kNcorrelation; i++)
    for(Int_t j=0; j<kNharmonics; j++)
      for(Int_t k=0; k<4; k++){
        if(fVn[i][j][k]) if(fVn[i][j][k]->GetEntries()>0) fEventQAList->Add(fVn[i][j][k]);
  }

  for(Int_t i=0; i<kNcorrelation*3; i++)
    for(Int_t j=0; j<kNharmonics; j++)
      for(Int_t k=0; k<4; k++){
        if(fCorrelation[i][j][k]) if(fCorrelation[i][j][k]->GetEntries()>0) fEventQAList->Add(fCorrelation[i][j][k]);
  }


}



//__________________________________________________________________
Bool_t AliAnalysisTaskQnAnalysis::IsEventSelected(Float_t* values) {
  if(!fEventCuts) return kTRUE;
  return fEventCuts->IsSelected(values);
}




//__________________________________________________________________
AliQnCorrectionsQnVector* AliAnalysisTaskQnAnalysis::GetQvector(TList* l, TString name, Int_t maxcorrection){
  TClonesArray* c  = dynamic_cast<TClonesArray*> (l->FindObject(Form("%s_%d",name.Data(),maxcorrection)));
  AliQnCorrectionsQnVector* q = 0x0;
  if(c) q = (AliQnCorrectionsQnVector*) c->At(0);
  return q;
}

