#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliITSRecPoint.h"
#include "AliESDEvent.h"
#include "AliTrackPointArray.h"
#include "AliITSgeomTGeo.h"
#include "AliESDfriend.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliGeomManager.h"
#include <TSystem.h>
#include <TTree.h>
#include <TH1F.h>
#include <TChain.h>
#include <TGeoGlobalMagField.h>
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisTaskSDDRP.h"

ClassImp(AliAnalysisTaskSDDRP)
//______________________________________________________________________________
AliAnalysisTaskSDDRP::AliAnalysisTaskSDDRP() : AliAnalysisTask("SDD RecPoints",""), 
  fOutput(0),
  fHistNEvents(0),
  fRecPMod(0),
  fTrackPMod(0),
  fGoodAnMod(0),
  fRecPLadLay3(0),
  fRecPLadLay4(0),
  fTrackPLadLay3(0),
  fTrackPLadLay4(0),
  fGoodAnLadLay3(0),
  fGoodAnLadLay4(0),
  fESD(0),
  fESDfriend(0),
  fResp(0),
  fRunNumber(0),
  fMinITSpts(3),
  fMinPfordEdx(1.5),
  fOnlyCINT1BTrig(0)
{
  //
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskSDDRP::~AliAnalysisTaskSDDRP(){
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//___________________________________________________________________________
void AliAnalysisTaskSDDRP::ConnectInputData(Option_t *) {
  //
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if(!tree) {
    printf("ERROR: Could not read chain from input slot 0\n");
  }else{
    tree->SetBranchStatus("ESDfriend*", 1);
    tree->SetBranchAddress("ESDfriend.",&fESDfriend);
    AliESDInputHandlerRP *hand = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    fESD=hand->GetEvent();
  }
  return;
}
//___________________________________________________________________________

void AliAnalysisTaskSDDRP::CreateOutputObjects() {
  //
  AliGeomManager::LoadGeometry(fGeomFile.Data());
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(fRunNumber);
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  
  AliCDBEntry* eR=(AliCDBEntry*)man->Get("ITS/Calib/RespSDD");
  fResp=(AliITSresponseSDD*)eR->GetObject();

  AliCDBEntry* eC=(AliCDBEntry*)man->Get("ITS/Calib/CalibSDD");
  TObjArray* calsdd=(TObjArray*)eC->GetObject();
  Int_t countGood3[14];
  Int_t countGood4[22];
  Int_t countGoodMod[260];
  for(Int_t ilad=0;ilad<14;ilad++) countGood3[ilad]=0;
  for(Int_t ilad=0;ilad<22;ilad++) countGood4[ilad]=0;
  for(Int_t imod=0;imod<260;imod++) countGoodMod[imod]=0;
  for(Int_t imod=0;imod<260;imod++){
    AliITSCalibrationSDD* cal=(AliITSCalibrationSDD*)calsdd->At(imod);
    if(cal->IsBad()) continue;
    Int_t modId=imod+AliITSgeomTGeo::GetModuleIndex(3,1,1);
    Int_t lay,lad,det;
    AliITSgeomTGeo::GetModuleId(modId,lay,lad,det);
    if(!CheckModule(lay,lad,det)) continue;
    for(Int_t ian=0; ian<512; ian++){
      if(cal->IsBadChannel(ian)) continue;
      countGoodMod[imod]++;
      if(lay==3) countGood3[lad-1]++;
      else if(lay==4) countGood4[lad-1]++;
    }
  }


  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  // -- Module histos

  fRecPMod = new TH1F("hRPMod","Rec Points per Module",260,239.5,499.5);
  fRecPMod->Sumw2();
  fRecPMod->SetMinimum(0);
  fOutput->Add(fRecPMod);

  fTrackPMod = new TH1F("hTPMod","Track Points per Module",260,239.5,499.5);
  fTrackPMod->Sumw2();
  fTrackPMod->SetMinimum(0);
  fOutput->Add(fTrackPMod);

  fGoodAnMod = new TH1F("hGAMod","Good Anodes per Module",260,239.5,499.5);
  for(Int_t imod=0;imod<260;imod++) fGoodAnMod->SetBinContent(imod+1,countGoodMod[imod]);
  fGoodAnMod->SetMinimum(0);
  fOutput->Add(fGoodAnMod);

  // -- Ladder histos

  fRecPLadLay3 = new TH1F("hRPLad3","Rec Points per Ladder Layer 3",14,-0.5,13.5);
  fRecPLadLay3->Sumw2();
  fRecPLadLay3->SetMinimum(0);
  fOutput->Add(fRecPLadLay3);

  fRecPLadLay4 = new TH1F("hRPLad4","Rec Points per Ladder Layer 4",22,-0.5,21.5);
  fRecPLadLay4->Sumw2();
  fRecPLadLay4->SetMinimum(0);
  fOutput->Add(fRecPLadLay4);

  fTrackPLadLay3 = new TH1F("hTPLad3","Track Points per Ladder Layer 3",14,-0.5,13.5);
  fTrackPLadLay3->Sumw2();
  fTrackPLadLay3->SetMinimum(0);
  fOutput->Add(fTrackPLadLay3);

  fTrackPLadLay4 = new TH1F("hTPLad4","Track Points per Ladder Layer 4",22,-0.5,21.5);
  fTrackPLadLay4->Sumw2();
  fTrackPLadLay4->SetMinimum(0);
  fOutput->Add(fTrackPLadLay4);

  fGoodAnLadLay3 = new TH1F("hGALad3","Good Anodes per Ladder Layer 3",14,-0.5,13.5);
  for(Int_t ilad=0;ilad<14;ilad++) fGoodAnLadLay3->SetBinContent(ilad+1,countGood3[ilad]);
  fGoodAnLadLay3->SetMinimum(0);
  fOutput->Add(fGoodAnLadLay3);

  fGoodAnLadLay4 = new TH1F("hGALad4","Good Anodes per Ladder Layer 4",22,-0.5,21.5);
  for(Int_t ilad=0;ilad<22;ilad++) fGoodAnLadLay4->SetBinContent(ilad+1,countGood4[ilad]);
  fGoodAnLadLay4->SetMinimum(0);
  fOutput->Add(fGoodAnLadLay4);

  for(Int_t it=0; it<8; it++){
    fSignalTime[it]=new TH1F(Form("hSigTimeInt%d",it),Form("hSigTimeInt%d",it),100,0.,300.);
    fSignalTime[it]->Sumw2();
    fSignalTime[it]->SetMinimum(0);
    fOutput->Add(fSignalTime[it]);
  }
}
//______________________________________________________________________________
void AliAnalysisTaskSDDRP::Exec(Option_t *)
{
  if(!fESD) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESD\n");
    return;
  } 
  if(!fESDfriend) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESDfriend\n");
    return;
  } 
  PostData(0,fOutput);
  fESD->SetESDfriend(fESDfriend);
  fHistNEvents->Fill(0);
  if(fOnlyCINT1BTrig){
    if(!fESD->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL")) return;
    fHistNEvents->Fill(1);
  }
  const AliTrackPointArray *array = 0;
  Int_t ntracks = fESD->GetNumberOfTracks();
  for (Int_t itrack=0; itrack < ntracks; itrack++) {
    AliESDtrack * track = fESD->GetTrack(itrack);
    if (!track) continue;
    if(track->GetNcls(1)>0) continue;
    if(track->GetNcls(0) < fMinITSpts) continue;
    Double_t dedx[4];
    track->GetITSdEdxSamples(dedx);
    array = track->GetTrackPointArray();
    if(!array) continue;
    for(Int_t ipt=0; ipt<array->GetNPoints(); ipt++) {
      AliTrackPoint point;
      Int_t modId;
      array->GetPoint(point,ipt);
      Int_t volId = point.GetVolumeID();
      Int_t layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if(layerId<3 || layerId>4) continue;
      modId+=AliITSgeomTGeo::GetModuleIndex(layerId,1,1);
      Int_t lay,lad,det;
      AliITSgeomTGeo::GetModuleId(modId,lay,lad,det);
      if(!CheckModule(lay,lad,det)) continue;
      fTrackPMod->Fill(modId);
      Float_t dtime=point.GetDriftTime()-fResp->GetTimeZero(modId);
      Int_t theBin=int(dtime/6500.*8.);
      if(layerId==3){ 
	fTrackPLadLay3->Fill(lad-1);	  
	if(dedx[0]>0. && track->P()>fMinPfordEdx) fSignalTime[theBin]->Fill(dedx[0]);
      }
      if(layerId==4){
	fTrackPLadLay4->Fill(lad-1);
	if(dedx[0]>0.&& track->P()>fMinPfordEdx) fSignalTime[theBin]->Fill(dedx[1]);
      }
    }
  }

  AliESDInputHandlerRP *hand = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  TTree* tR = hand->GetTreeR("ITS");
  TClonesArray *ITSrec= new TClonesArray("AliITSRecPoint");
  TBranch *branch =tR->GetBranch("ITSRecPoints");
  branch->SetAddress(&ITSrec);
  for (Int_t modId=240; modId<500; modId++){
    Int_t lay,lad,det;
    AliITSgeomTGeo::GetModuleId(modId,lay,lad,det);
    if(!CheckModule(lay,lad,det)) continue;
    branch->GetEvent(modId);
    Int_t nrecp = ITSrec->GetEntries();	
    fRecPMod->Fill(modId,nrecp);	  
    if(lay==3) fRecPLadLay3->Fill(lad-1,nrecp);
    if(lay==4) fRecPLadLay4->Fill(lad-1,nrecp);
    
  }
  ITSrec->Delete();
  delete ITSrec;

  PostData(0,fOutput);
  
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskSDDRP::CheckModule(Int_t lay, Int_t lad, Int_t det) const{
  //
  if(lay==4){
    if(lad==3 && det==5) return kFALSE; // 1500 V
    if(lad==3 && det==6) return kFALSE; // 0 V
    if(lad==3 && det==7) return kFALSE; // 1500 V
    if(lad==4 && det==1) return kFALSE; // 0 V
    if(lad==4 && det==2) return kFALSE; // 1500 V
    if(lad==7 && det==3) return kFALSE; // noisy
    if(lad==7 && det==5) return kFALSE; // 0 MV + noisy
    if(lad==9 && det==3) return kFALSE; // 1500 V
    if(lad==9 && det==4) return kFALSE; // 0 V
    if(lad==9 && det==5) return kFALSE; // 1500 V
    if(lad==11 && det==6) return kFALSE; // 1500 V
    if(lad==11 && det==7) return kFALSE; // 0 V
    if(lad==11 && det==8) return kFALSE; // 1500 V
    if(lad==18 && det==5) return kFALSE; // 1500 V
    if(lad==18 && det==6) return kFALSE; // 0 V
    if(lad==18 && det==7) return kFALSE; // 1500 V
    if(lad==22 && det==1) return kFALSE; // 0 V
    if(lad==22 && det==2) return kFALSE; // 1500 V
  }
  if(lay==3){
    if(lad==4 && det==4) return kFALSE; // 1500 V 
    if(lad==3) return kFALSE;  // swapped in geometry
    if(lad==5 && det==1) return kFALSE; // noisy
    if(lad==5 && det==2) return kFALSE; // noisy
    if(lad==6 && det==1) return kFALSE; // noisy
  }
  return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisTaskSDDRP::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  fRecPMod= dynamic_cast<TH1F*>(fOutput->FindObject("hRPMod"));
  fTrackPMod= dynamic_cast<TH1F*>(fOutput->FindObject("hTPMod"));
  fGoodAnMod= dynamic_cast<TH1F*>(fOutput->FindObject("hGAMod"));

  fRecPLadLay3= dynamic_cast<TH1F*>(fOutput->FindObject("hRPLad3"));
  fRecPLadLay4= dynamic_cast<TH1F*>(fOutput->FindObject("hRPLad4"));
  fTrackPLadLay3= dynamic_cast<TH1F*>(fOutput->FindObject("hTPLad3"));
  fTrackPLadLay4= dynamic_cast<TH1F*>(fOutput->FindObject("hTPLad4"));
  fGoodAnLadLay3= dynamic_cast<TH1F*>(fOutput->FindObject("hGALad3"));
  fGoodAnLadLay4= dynamic_cast<TH1F*>(fOutput->FindObject("hGALad4"));

  for(Int_t it=0; it<8; it++){
    fSignalTime[it]= dynamic_cast<TH1F*>(fOutput->FindObject(Form("hSigTimeInt%d",it)));
  }

  return;
}





