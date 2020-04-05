#ifndef ALIESDTOOLS_H
#define ALIESDTOOLS_H

class AliPIDResponse;
class TTreeSRedirector;
class TTreeStream;
class TTree;
class TGraph;
class TH1F;
class AliExternalTrackParam;
class AliESDEvent;
class AliESDfriend;
class AliTriggerAnalysis;
//class TVectorF;
#include "TNamed.h"

class AliESDtools : public TNamed {
  public:
  AliESDtools();
  void Init(TTree* tree, AliESDEvent *event= nullptr);
  void SetStreamer(TTreeSRedirector *streamer){fStreamer=streamer;}
  static Double_t LoadESD(Int_t entry, Int_t verbose=0);
  /// caching
  Int_t  CacheTPCEventInformation();
  Int_t  CacheITSVertexInformation(Bool_t doReset=1, Double_t dcaCut=0.05,  Double_t dcaZcut=0.15);
  Int_t  CacheTOFEventInformation(Bool_t dumpStreamer=0);
  Int_t CalculateEventVariables();
  void TPCVertexFit(TH1F *hisVertex);
  Int_t  GetNearestTrack(const AliExternalTrackParam * trackMatch, Int_t indexSkip, AliESDEvent*event, Int_t trackType, Int_t paramType, AliExternalTrackParam & paramNearest);
  void   ProcessITSTPCmatchOut(AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, TTreeStream *pcstream);
  Double_t CachePileupVertexTPC(Int_t entry, Int_t doReset=0, Int_t verbose=0);
  //
  Int_t DumpEventVariables();
  static Int_t SDumpEventVariables(){return fgInstance->DumpEventVariables();}
  // static functions for querying cached variables in TTree formula
  static Int_t    SCalculateEventVariables(Int_t entry){LoadESD(entry,0); return fgInstance->CalculateEventVariables();}
  static Double_t GetTrackCounters(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackCounters)[index];}
  static Double_t GetTrackTPCCountersZ(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackTPCCountersZ)[index];}
  static Double_t GetTrackdEdxRatio(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackdEdxRatio)[index];}
  static Double_t GetTrackNcl(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackNcl)[index];}
  static Double_t GetTrackChi2(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackChi2)[index];}
  static Double_t GetTrackMatchEff(Int_t index, Int_t toolIndex){return (*fgInstance->fCacheTrackMatchEff)[index];}
  static Double_t GetMeanHisTPCVertexA(){return fgInstance->fHisTPCVertexA->GetMean();}
  static Double_t GetMeanHisTPCVertexC(){return fgInstance->fHisTPCVertexC->GetMean();}
  static Double_t GetVertexInfo(Int_t index){return (*fgInstance->fTPCVertexInfo)[index];}
  static Int_t SetDefaultAliases(TTree* tree);
  static Double_t SCachePileupVertexTPC(Int_t entry, Int_t doReset, Int_t verbose=0){ return fgInstance->CachePileupVertexTPC(entry, doReset, verbose);}
  static Double_t SCacheITSVertexInformation(Bool_t doReset=1, Double_t dcaCut=0.05, Double_t dcaZcut=0.15){ return fgInstance->CacheITSVertexInformation(doReset, dcaCut, dcaZcut);}
  //
  static Double_t SCacheTOFEventInformation(Bool_t dumpStreamer){return fgInstance->CacheTOFEventInformation(dumpStreamer);}
  static Double_t SGetTOFHitInfo(Int_t number, Int_t coord);
  //
  Int_t fVerbose;                                 // verbosity flag
  TTree *fESDtree;                                //! esd Tree pointer - class is not owner
  AliESDEvent * fEvent;                           //! esd event pointer - class is not owner
  AliPIDResponse   * fPIDResponse;                //! PID response object
  AliTriggerAnalysis *fTriggerAnalysis;           //! tigger analysis
  Bool_t   fTaskMode;                             // analysis task mode
  TH1F *fHisITSVertex;                            // ITS z vertex histogram
  TH1F *fHisTPCVertexA;                           // TPC z vertex A side
  TH1F *fHisTPCVertexC;                           // TPC z vertex C side
  TH1F *fHisTPCVertex;                            //
  TH1F *fHisTPCVertexACut;                        //
  TH1F *fHisTPCVertexCCut;                        //

  TVectorF         * fTPCVertexInfo;              // TPC vertex info
  TVectorF         * fITSVertexInfo;              // ITS vertex info
  TH1F             * fHistPhiTPCCounterA;         // helper histogram phi counters
  TH1F             * fHistPhiTPCCounterC;         // helper histogram phi counters
  TH1F             * fHistPhiTPCCounterAITS;      // helper histogram phi counters
  TH1F             * fHistPhiTPCCounterCITS;      // helper histogram phi counters
  TH1F             * fHistPhiITSCounterA;         // helper histogram phi counters
  TH1F             * fHistPhiITSCounterC;         // helper histogram phi counters
  TVectorF         * fCacheTrackCounters;         // track counter
  TVectorF         * fCacheTrackTPCCountersZ;     // track counter with DCA z cut
  TVectorF         * fCacheTrackdEdxRatio;        // dEdx info counter
  TVectorF         * fCacheTrackNcl;              // ncl counter
  TVectorF         * fCacheTrackChi2;             // chi2 counter
  TVectorF         * fCacheTrackMatchEff;         // matchEff counter
  TGraph           * fLumiGraph;                  // graph for the interaction rate info for a run
  //
  TTreeSRedirector * fStreamer;                  /// streamer
  static AliESDtools* fgInstance;                /// instance of the tool -needed in order to use static functions (for TTreeFormula)
  private:
  AliESDtools(AliESDtools&);
  AliESDtools &operator=(const AliESDtools&);
  ClassDef(AliESDtools, 1) 
};

#endif
