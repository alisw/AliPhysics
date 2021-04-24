/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// author: Artur Furs, afurs@cern.ch
// Analysis task for conditional TZERO spectrum analysis
//

#ifndef AliAnalysisTaskSpectrumT0_cxx
#define AliAnalysisTaskSpectrumT0_cxx
//STD
#include <iostream>
#include <map>
#include <utility>
#include <bitset>
#include <set>
#include <cstdio>

#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <functional>
#endif

//Custom
#include "EventStruct.h"
//ROOT
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TChain.h"
//AliRoot
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDTZERO.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliVEvent.h"
#if !(defined(__CINT__) || defined(__MAKECINT__))
template<typename HistType=TH1F,std::size_t NCHANNELS=24,std::size_t NTRGSIGNALS=3>
class OutputData {
 public:
  OutputData(TList *listOutput,TString nameSuffix, TString conditions) {
    Init(listOutput,nameSuffix,conditions);
  }
  virtual ~OutputData() {
    if(mIsReady) {
//      delete mListAll;
//      mListAll=NULL;
    }
  }
  typedef HistType Hist_t;
  constexpr static std::size_t sNCHANNELS = NCHANNELS;
  constexpr static std::size_t sNTRGSIGNALS = NTRGSIGNALS;
  enum ETRGSIGNALS {kOrA,kOrC,kTVDC};
  std::vector<Hist_t*> mAmp; //!
  std::vector<Hist_t*> mAmpNew; //!
  std::vector<Hist_t*> mTime; //!
  std::vector<Hist_t*> mTimeFull; //!
  std::vector<Hist_t*> mTimeTrg; //!
  bool mIsReady=false; //!
  TList *mListAll=NULL; //!
  void FillT0data(const AliESDEvent* esdData) {
    const AliESDTZERO* esdT0= dynamic_cast<const AliESDTZERO*> (esdData->GetESDTZERO());
    if(esdT0!=nullptr)	{
      const int nHit = 0;
      const Double32_t *time=esdT0->GetT0time();
      const  Double32_t *amp=esdT0->GetT0amplitude();
      const  Double32_t *ampNew=esdT0->GetT0NewAmplitude();
      for(int iCh=0;iCh<NCHANNELS;iCh++) {
        if(amp[iCh]>0) mAmp[iCh]->Fill(amp[iCh]);
        if(ampNew[iCh]>0) mAmpNew[iCh]->Fill(ampNew[iCh]);
        if(time[iCh]>0) mTime[iCh]->Fill(time[iCh]);
        Float_t timeFull=esdT0->GetTimeFull(iCh, nHit);
        if(timeFull>0) mTimeFull[iCh]->Fill(timeFull);
      }
      Double_t orA = esdT0->GetOrA(nHit);
      Double_t orC = esdT0->GetOrC(nHit);
      Double_t tvdc = esdT0->GetTVDC(nHit);
      if(orA!=0) mTimeTrg[kOrA]->Fill(orA);
      if(orC!=0) mTimeTrg[kOrC]->Fill(orC);
      if(tvdc!=0) mTimeTrg[kTVDC]->Fill(tvdc);
    }
  }
 private:
  void Init(TList *listOutput,TString nameSuffix, TString conditions) {
    if(listOutput==NULL) return;
    if(mListAll==NULL) {
      mListAll=new TList();
      mIsReady=true;
    }
    else {
      return;
    }
    listOutput->Add(mListAll);
    mListAll->SetName("listOutput"+nameSuffix);
    mListAll->SetOwner(kTRUE);
    mAmp.resize(NCHANNELS);
    mAmpNew.resize(NCHANNELS);
    mTime.resize(NCHANNELS);
    mTimeFull.resize(NCHANNELS);
    mTimeTrg.resize(NTRGSIGNALS);
    TList *listAmp = new TList();
    listAmp->SetName("listAmp"+nameSuffix);
    mListAll->Add(listAmp);
    TList *listAmpNew = new TList();
    listAmpNew->SetName("listAmpNew"+nameSuffix);
    mListAll->Add(listAmpNew);
    TList *listTime = new TList();
    listTime->SetName("listTime"+nameSuffix);
    mListAll->Add(listTime);
    TList *listTimeFull = new TList();
    listTimeFull->SetName("listTimeFull"+nameSuffix);
    mListAll->Add(listTimeFull);
    TList *listTimeTrg = new TList();
    listTimeTrg->SetName("listTimeTrg"+nameSuffix);
    mListAll->Add(listTimeTrg);
    for(auto entry: (*mListAll)) {
      TList *listEntry = dynamic_cast<TList *>(entry);
      if(listEntry==NULL) continue;
      listEntry->SetOwner(kTRUE);
      //listOutput->Add(listEntry);
    }
    for(int iCh=0;iCh<NCHANNELS;iCh++) {
      Hist_t *hist;
      //Amplitudes
      hist = new Hist_t(Form("hAmp_ch%i"+nameSuffix,iCh+1),Form("Amplitude ch%i, " + conditions,iCh+1),500,0,10);
      mAmp[iCh] = hist;
      listAmp->Add(hist);
      //New amplitudes
      hist = new Hist_t(Form("hAmpNew_ch%i"+nameSuffix,iCh+1),Form("Amplitude(new) ch%i, " + conditions,iCh+1),500,0,10);
      mAmpNew[iCh] = hist;
      listAmpNew->Add(hist);
      //Time
      hist = new Hist_t(Form("hTime_ch%i"+nameSuffix,iCh+1),Form("Time ch%i, " + conditions,iCh+1),1000,9000,11000);
      mTime[iCh] = hist;
      listTime->Add(hist);
      //Time full
      hist = new Hist_t(Form("hTimeFull_ch%i"+nameSuffix,iCh+1),Form("Time(full) ch%i, " + conditions,iCh+1),250,-5,5);
      mTimeFull[iCh] = hist;
      listTimeFull->Add(hist);
    }
    {
      Hist_t *hist;
      //OrA
      hist = new Hist_t("hOrA"+nameSuffix,"OrA , " + conditions,250,-5,5);
      mTimeTrg[kOrA] = hist;
      listTimeTrg->Add(hist);
      //OrC
      hist = new Hist_t("hOrC"+nameSuffix,"OrC , " + conditions,250,-5,5);
      mTimeTrg[kOrC] = hist;
      listTimeTrg->Add(hist);
      //TVDC
      hist = new Hist_t("hTVDC"+nameSuffix,"TVDC , " + conditions,250,-5,5);
      mTimeTrg[kTVDC] = hist;
      listTimeTrg->Add(hist);
    }
  }
};
#endif
class AliAnalysisTaskSpectrumT0 : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskSpectrumT0():AliAnalysisTaskSE()
    ,mESD(0x0), mOutputList(0x0),mOutputListPerRun(0x0)
    ,mHistRunStats(0x0)
    {}
    AliAnalysisTaskSpectrumT0(const char *name);
    virtual ~AliAnalysisTaskSpectrumT0();
    virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
    virtual void   NotifyRun();
		virtual void   Terminate(Option_t *);
    Bool_t UserNotify();
    //Type fields
    typedef TH1F Hist_t;
    typedef std::bitset<11> EventID_t;
    typedef std::set<unsigned int> SetRunNums_t;
    struct Comparer {
      bool operator()(const EventID_t& eventCutID1,const EventID_t eventCutID2) const {
        return eventCutID1.to_ullong()<eventCutID2.to_ullong();
      }
    };
    typedef std::set<EventID_t,Comparer> SetEventID_t;
    //
    //Only for ROOT6 !!!
    #if !(defined(__CINT__) || defined(__MAKECINT__))
    typedef OutputData<Hist_t,24,3> OutputData_t;
    typedef std::map<EventID_t,OutputData_t,Comparer> MapEventID2Output_t;
    typedef std::function<bool(const EventStruct&)> Functor_t;
    struct FunctorStruct {
      std::string mName; //!
      std::string mTitle; //!
      Functor_t mFunct; //!
    };
    static const std::vector<FunctorStruct> sVecFunctors; //!
    std::map<std::size_t, FunctorStruct> mMapBitPos2Funct; //!
    std::reference_wrapper<MapEventID2Output_t> mRefMapEventID2Output = mMapEventID2Output; //!
    void PrintCutBits() const {
      std::cout<<"\n=======CUT BIT NAMES=======";
      for(const auto &entry: mMapBitPos2Funct) {
        std::cout<<"\n|position: "<<entry.first;
        std::cout<<" |name: "<<entry.second.mName;
        std::cout<<" |title: "<<entry.second.mTitle;
      }
      std::cout<<"\n===========================\n";
    }
    void PrintEventIDs() {
      std::cout<<"\n=======EventIDs=======";
      for(const auto &entry: mSetEventID) {
        std::cout<<"\n|bitset: "<<entry.to_string();
        std::cout<<" |name: "<<GetNameEventID(entry);

      }
      std::cout<<"\n=====================\n";
    }
    TString GetNameEventID(const EventID_t& eventID) {
      if(eventID.count()==0) return TString{};
      std::string stName{};
      for(const auto &entryBits: mMapBitPos2Funct) {
        if(eventID.test(entryBits.first)) {
          stName+=entryBits.second.mName;
          stName+="_";
        }
      }
      return TString{stName.substr(0,stName.size()-1)};
    }
    void UpdateEventID() {
      mEventID.reset();
      for(const auto &entry:mMapBitPos2Funct) {
        auto &funct = entry.second.mFunct;
        bool isFired = funct(mEventStruct);
        if(isFired)mEventID.set(entry.first);
      }
    }
    void FillHists(const AliESDEvent* esdData) {
      for(auto &entry: mRefMapEventID2Output.get()) {
        auto &eventID = entry.first;
        auto &outputData = entry.second;
        if(this->IsAcceptableEventID(eventID)) {
          outputData.FillT0data(esdData);
        }
      }
    }
    #endif
    //END

    bool AddEventID(ULong64_t eventID) {return AddEventID(EventID_t{eventID});}
    bool AddEventID(const std::string &eventID) {return AddEventID(EventID_t{eventID});}
    bool AddEventID(const EventID_t& eventID) {
      auto pairInserted = mSetEventID.insert(eventID);
      return pairInserted.second;
    }
    bool IsAcceptableEventID(const EventID_t& targetEventID){
      return (mEventID&targetEventID) == targetEventID;
    }
    void GetRunNumsFromCSV(TString filepath);
    bool GetRunNumsFromROOT(TString filepath);
  private:
    AliESDEvent *mESD;    //! ESD object
    //Data container
    TList *mOutputListPerRun; //! Output list
    TList *mOutputList; //! Output list
    //Event struct
    EventStruct mEventStruct; //! Event structure
    //EventID bit set
    EventID_t mEventID; //!

    #if !(defined(__CINT__) || defined(__MAKECINT__))
    MapEventID2Output_t mMapEventID2Output; //! unused, only for reference wrapper during construction
    std::map<unsigned int,MapEventID2Output_t> mMapRun2Output; //!
    #endif
    std::map<unsigned int, unsigned int> mMapRun2BinPos; //!
    TH1F *mHistRunStats; //!
    //To store
    SetEventID_t mSetEventID; // Set of registered EventIDs
    SetRunNums_t mSetRunNums; // Set of run numbers, for output preparation
    //
    AliAnalysisTaskSpectrumT0(const AliAnalysisTaskSpectrumT0&); // not implemented
    AliAnalysisTaskSpectrumT0& operator=(const AliAnalysisTaskSpectrumT0&); // not implemented
    ClassDef(AliAnalysisTaskSpectrumT0, 1); // example of analysis
};
#endif
