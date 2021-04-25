/*************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
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

// author: Artur Furs, afurs@cern.ch
// Analysis task for conditional TZERO spectrum analysis
//

#include "AliAnalysisTaskSpectrumT0.h"
ClassImp(AliAnalysisTaskSpectrumT0)
/*******************************************************************************************************************/
AliAnalysisTaskSpectrumT0::AliAnalysisTaskSpectrumT0(const char *name): AliAnalysisTaskSE(name)
,mESD(0x0), mOutputList(0x0),mOutputListPerRun(0x0),mEventStruct{}
,mEventID{},mHistRunStats(0x0)
{
  /// Constructor

  /// Define input and output slots here
  /// Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  /// Output slot #0 id reserved by the base class for AOD
  /// Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}
/*******************************************************************************************************************/
AliAnalysisTaskSpectrumT0::~AliAnalysisTaskSpectrumT0()
{
// Destructor

  if (mOutputList)	delete mOutputList;
}
/*******************************************************************************************************************/
Bool_t AliAnalysisTaskSpectrumT0::UserNotify()	{
///
/// Calls the mother class Notify()
///
	return AliAnalysisTaskSE::UserNotify();
}
/*******************************************************************************************************************/
void AliAnalysisTaskSpectrumT0::NotifyRun()
{
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  if(mMapRun2Output.find(fCurrentRunNumber)==mMapRun2Output.end()) {
    //dummy run
    auto it = mMapRun2Output.find(0);
    mRefMapEventID2Output=it->second;
    mHistRunStats->Fill(0);
  }
  else {
    auto it = mMapRun2Output.find(fCurrentRunNumber);
    mRefMapEventID2Output=it->second;
    auto itCnt = mMapRun2BinPos.find(fCurrentRunNumber);
    if(itCnt!=mMapRun2BinPos.end())  {
      mHistRunStats->Fill(itCnt->second);
    }
  }
  //InitHists();
  #endif
}
/*******************************************************************************************************************/
void AliAnalysisTaskSpectrumT0::UserCreateOutputObjects()	{
  /// Create histograms
  /// Called once
  mOutputList = new TList();
  mOutputList->SetOwner(kTRUE); /// Will delete the histos on cleanup
  mOutputList->SetName("output");

  #if !(defined(__CINT__) || defined(__MAKECINT__))
  std::size_t bitPosition=0;
  mMapBitPos2Funct.clear();
  for(const auto& entry: sVecFunctors)  {
    mMapBitPos2Funct.insert({bitPosition,entry});
    bitPosition++;
  }
  ///
  //Temporary
  ///
  /*
  this->AddEventID("1011010010");
  this->AddEventID("10011010010");
  this->AddEventID("1");
  GetRunNumsFromROOT("runListRCT2018_pp_tzero.root");
  */
  ////////////////////////////
  //Map preparation
  mSetRunNums.insert(0);//dummy run
  std::size_t binPos=0;
  //mMapRun2BinPos.insert({0,binPos});
  mHistRunStats = new TH1F("hRunStat","hRunStat",mSetRunNums.size(),0,mSetRunNums.size());
  mOutputList->Add(mHistRunStats);
  for(const auto& runnum: mSetRunNums) {
    mMapRun2BinPos.insert({runnum,binPos});
    mHistRunStats->GetXaxis()->SetBinLabel(binPos+1,Form("%i",runnum));
    auto itInserted = mMapRun2Output.insert({runnum,{}});
    mRefMapEventID2Output = itInserted.first->second;
    mOutputListPerRun = new TList();
    mOutputListPerRun->SetName(Form("run%i",runnum));
    mOutputListPerRun->SetOwner(kTRUE);
    mOutputList->Add(mOutputListPerRun);
    for(const auto &entry: mSetEventID) {
      TString stNameSuffix = Form("_v%i_run%i",entry.to_ullong(),runnum);
      TString stConditions = GetNameEventID(entry);
      mRefMapEventID2Output.get().insert({entry,OutputData_t{mOutputListPerRun,stNameSuffix,stConditions}});
    }
    binPos++;
  }
  auto it = mMapRun2Output.find(0);
  mRefMapEventID2Output=it->second;
  #endif
  PostData(1, mOutputList);
}
/*******************************************************************************************************************/
void AliAnalysisTaskSpectrumT0::GetRunNumsFromCSV(TString filepath) {
  TTree treeCSV;
  treeCSV.ReadFile(filepath);
  int runnum;
  treeCSV.SetBranchAddress("raw_run",&runnum);
  for(int iEntry=0;iEntry<treeCSV.GetEntries();iEntry++) {
    treeCSV.GetEntry(iEntry);
    mSetRunNums.insert(static_cast<unsigned int>(runnum));
  }
}
/*******************************************************************************************************************/
bool AliAnalysisTaskSpectrumT0::GetRunNumsFromROOT(TString filepath) {
  TFile *inputFile = TFile::Open(filepath);
  if(inputFile==NULL) return false;
  if(inputFile->IsZombie()) return false;
  TTree *inputTree = dynamic_cast<TTree *>(inputFile->Get("RCT"));
  if(inputTree==NULL) return false;
  if(inputTree->GetEntries()==0) return false;
  int runnum;
  inputTree->SetBranchAddress("raw_run",&runnum);
  for(int iEntry=0;iEntry<inputTree->GetEntries();iEntry++) {
    inputTree->GetEntry(iEntry);
    mSetRunNums.insert(static_cast<unsigned int>(runnum));
  }
  delete inputTree;
  inputFile->Close();
  delete inputFile;
  return true;
}
/*******************************************************************************************************************/
void AliAnalysisTaskSpectrumT0::UserExec(Option_t *)	{
  /// Main loop
  /// Called for each event
  mESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(mESD!=NULL)	{
    mEventStruct.clear();
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    //Data accumulating

    //Global vertex
    const AliESDVertex *esdVertexGlobal = dynamic_cast<const AliESDVertex *>(mESD->GetPrimaryVertex());
    if(esdVertexGlobal!=NULL)	{
        mEventStruct.fNcontGlobal = esdVertexGlobal->GetNContributors();
        mEventStruct.fVertexGlobal_Z=esdVertexGlobal->GetZ();
        mEventStruct.fVtxGlobalResZ = esdVertexGlobal->GetZRes();
    }
    //Track vertex
    const AliESDVertex *esdVertexTrack = dynamic_cast<const AliESDVertex *>(mESD->GetPrimaryVertexTracks());
    if(esdVertexTrack!=NULL)	{
        mEventStruct.fNcontTrack = esdVertexTrack->GetNContributors();
        mEventStruct.fVertexTrack_Z=esdVertexTrack->GetZ();
        mEventStruct.fVtxTrackResZ = esdVertexTrack->GetZRes();
    }
    //V0 timing
    AliESDVZERO* esdV0 = dynamic_cast<AliESDVZERO*>(mESD->GetVZEROData());
    if(esdV0!=NULL)	{
       mEventStruct.fTimeV0C = esdV0->GetV0CTime();
       mEventStruct.fTimeV0A = esdV0->GetV0ATime();
    }
    //Pileup
    mEventStruct.fIsPileup = mESD->IsPileupFromSPDInMultBins();
    mEventStruct.fIsPileupLowMult = mESD->IsPileupFromSPD(3.,0.8,3.,2.,5.);
    //Physics selection
    if(inputHandler!=NULL) {
      mEventStruct.fIsPhysSel = inputHandler->IsEventSelected();
    }
    //CTP triggers
    mEventStruct.fTrigger.SetString(mESD->GetFiredTriggerClasses());
    mEventStruct.fTriggerMask = mESD->GetTriggerMask();
    mEventStruct.fTriggerMaskNext50 = mESD->GetTriggerMaskNext50();

    #if !(defined(__CINT__) || defined(__MAKECINT__))
    //auto it = mMapRun2Output.find(0);
    //mRefMapEventID2Output=it->second;
    this->UpdateEventID();
    this->FillHists(mESD);
    #endif
  }
  else {
    //printf("ERROR: mESD not available\n");
    return;
  }
  PostData(1, mOutputList);
}      
/*******************************************************************************************************************/
void AliAnalysisTaskSpectrumT0::Terminate(Option_t *)	{

}

#if !(defined(__CINT__) || defined(__MAKECINT__))
const std::vector<AliAnalysisTaskSpectrumT0::FunctorStruct> AliAnalysisTaskSpectrumT0::sVecFunctors ={
    {"noCuts","No cuts",
   [](const EventStruct& data) {
     return true;
   }},

  {"trkVtx","Track vertex, Ncont>1 && |Z|<10",
   [](const EventStruct& data) {
     return (data.fNcontTrack>1) && (data.fVertexTrack_Z>-10) && (data.fVertexTrack_Z<10);
   }},

  {"trkVtx2","Track vertex, Ncont>1 && |Z|<30",
   [](const EventStruct& data) {
     return (data.fNcontTrack>1) && (data.fVertexTrack_Z>-30) && (data.fVertexTrack_Z<30);
   }},

  {"glbVtx","Global vertex, Ncont>1 && |Z|<10",
   [](const EventStruct& data) {
     return (data.fNcontGlobal>0) && (data.fVertexGlobal_Z>-10) && (data.fVertexGlobal_Z<10);
   }},

  {"tV0","V0 timing cut, 11.5<sumV0<17.5 && 5.5<diffV0<11.5",
   [](const EventStruct& data) {
     return ((data.fTimeV0A+data.fTimeV0C)>11.5) && ((data.fTimeV0A+data.fTimeV0C)<17.5)
             && ((data.fTimeV0A-data.fTimeV0C)>5.5) &&  ((data.fTimeV0A-data.fTimeV0C)<11.5);
   }},

  {"tV02","V0 timing cut, 10<sumV0<18 && 4<diffV0<12",
   [](const EventStruct& data) {
     return ((data.fTimeV0A+data.fTimeV0C)>10) && ((data.fTimeV0A+data.fTimeV0C)<18)
             && ((data.fTimeV0A-data.fTimeV0C)>4) &&  ((data.fTimeV0A-data.fTimeV0C)<12);
   }},

  {"PS","Physics selection",[](const EventStruct& data) {
     return data.fIsPhysSel;
   }},

  {"PU","Pileup rejection",[](const EventStruct& data) {
     return !data.fIsPileup;
   }},

  {"PUlm","Pileup rejection(low mult)",
   [](const EventStruct& data) {
     return !data.fIsPileupLowMult;
   }},

  {"INT7","CINT7-B",
   [](const EventStruct& data) {
    TString trg = data.fTrigger.GetString();
    return trg.Contains("CINT7-B");
   }},

  {"0TVX","C0TVX-B && CINT7-B",
   [](const EventStruct& data) {
    TString trg = data.fTrigger.GetString();
    return (trg.Contains("C0TVX-B")&&trg.Contains("CINT7-B"));
   }}
};
#endif
