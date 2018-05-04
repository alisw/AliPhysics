#ifndef ALIESDTOOLS_H
#define ALESDTOOLS_H


class AliESDtools : public TNamed
{
  public:
  AliESDtools();
  void Init(TTree* tree);
  Int_t fVerbose;
  TTree *fESDtree;
  AliESDEvent * fEvent;
  Int_t  CacheTPCEventInformation();
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
  void CalculateEventVariables();
  void TPCVertexFit(TH1F *hisVertex);
  ClassDef(AliESDtools, 1) 
};

#endif
