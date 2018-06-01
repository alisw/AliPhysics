#ifndef ALIESDTOOLS_H
#define ALESDTOOLS_H


class AliESDtools : public TNamed
{
  public:
  AliESDtools();
  void Init(TTree* tree);
  /// caching
  Int_t  CacheTPCEventInformation();
  Int_t CalculateEventVariables();
  void TPCVertexFit(TH1F *hisVertex);
  Int_t  GetNearestTrack(const AliExternalTrackParam * trackMatch, Int_t indexSkip, AliESDEvent*event, Int_t trackType, Int_t paramType, AliExternalTrackParam & paramNearest);
  void   ProcessITSTPCmatchOut(AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, TTreeStream *pcstream);
  // static functions for querying in TTree formula
  static Int_t    SCalculateEventVariables(Int_t entry){fgInstance->fESDtree->GetEntry(entry); return fgInstance->CalculateEventVariables();}
  static Double_t GetTrackCounters(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackCounters)[index];}
  static Double_t GetTrackTPCCountersZ(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackTPCCountersZ)[index];}
  static Double_t GetTrackdEdxRatio(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackdEdxRatio)[index];}
  static Double_t GetTrackNcl(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackNcl)[index];}
  static Double_t GetTrackChi2(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackChi2)[index];}
  static Double_t GetTrackMatchEff(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackMatchEff)[index];}
  static Double_t GetMeanHisTPCVertexA(){return fgInstance->fHisTPCVertexA->GetMean();}
  static Double_t GetMeanHisTPCVertexC(){return fgInstance->fHisTPCVertexC->GetMean();}
  //
  Int_t fVerbose;
  TTree *fESDtree;
  AliESDEvent * fEvent;
  TH1F *fHisTPCVertexA;
  TH1F *fHisTPCVertexC;
  TH1F *fHisTPCVertex;
  TH1F *fHisTPCVertexACut;
  TH1F *fHisTPCVertexCCut;
  TH1F             * fHistPhiTPCcounterA;         // helper histogram phi counteres
  TH1F             * fHistPhiTPCcounterC;         // helper histogram phi counters
  TH1F             * fHistPhiTPCcounterAITS;      // helper histogram phi counters
  TH1F             * fHistPhiTPCcounterCITS;      // helper histogram phi counters
  TH1F             * fHistPhiITScounterA;         // helper histogram phi counters
  TH1F             * fHistPhiITScounterC;         // helper histogram phi counters
  TVectorF         * fCacheTrackCounters;         // track counter
  TVectorF         * fCacheTrackTPCCountersZ;     // track counter with DCA z cut
  TVectorF         * fCacheTrackdEdxRatio;        // dedx info counter
  TVectorF         * fCacheTrackNcl;              // ncl counter
  TVectorF         * fCacheTrackChi2;             // chi2 counter
  TVectorF         * fCacheTrackMatchEff;         // matchEff counter
  TGraph           * fLumiGraph;                  // graph for the interaction rate info for a run
  //
  static AliESDtools* fgInstance;                /// instance of the tool
  ClassDef(AliESDtools, 1) 
};

#endif
