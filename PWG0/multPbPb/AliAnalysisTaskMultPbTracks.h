#ifndef ALIANALYSISTASKMULTPBTRACKS_H
#define ALIANALYSISTASKMULTPBTRACKS_H

#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h" // if I don't include this, nothing compiles

//-------------------------------------------------------------------------
//                      AliAnalysisTaskMultPbTracks
// 
// 
//
//
// Author: Michele Floris, CERN
//-------------------------------------------------------------------------


class AliESDEvent;
class AliESDtrackCuts;
class AliAnalysisMultPbCentralitySelector;
class AliAnalysisMultPbTrackHistoManager;



class AliAnalysisTaskMultPbTracks : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskMultPbTracks();
  AliAnalysisTaskMultPbTracks(const char * name);
  AliAnalysisTaskMultPbTracks(const AliAnalysisTaskMultPbTracks& obj) ;
  ~AliAnalysisTaskMultPbTracks();
  //void SetCentralitySelector(AliAnalysisMultPbCentralitySelector * centr) { fCentrSelector=centr;}
  void SetTrackCuts(AliESDtrackCuts * cuts) { fTrackCuts = cuts;}
  void SetCentralityBin(Int_t bin = 0) { fCentrBin = bin; }
  void SetCentralityEstimator(const char * centr) { fCentralityEstimator = centr; }

  void SetIsMC(Bool_t flag=kTRUE) { fIsMC = flag;}
  AliAnalysisMultPbTrackHistoManager * GetHistoManager() { return fHistoManager;}
  Bool_t IsPhysicalPrimaryAndTransportBit(Int_t ipart) ;

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  

private:

  //
  AliESDEvent *  fESD;    //! ESD object  AliVEvent*     fEvent;
  //  TList * fListHisto;     // list of output object
  AliAnalysisMultPbTrackHistoManager  * fHistoManager; // wrapper for the list, takes care of merging + histo booking and getters
  //  AliAnalysisMultPbCentralitySelector * fCentrSelector; // centrality selector
  Int_t fCentrBin; // centrality bin selected (5% XS percentiles)
  TString fCentralityEstimator; // Name of the centrality estimator, for AliESDCentrality
  AliESDtrackCuts * fTrackCuts; // track cuts
  AliESDtrackCuts * fTrackCutsNoDCA; // copy of the previous one, but with no DCA cuts
  
  Bool_t fIsMC; // true if processing montecarlo

  AliAnalysisTaskMultPbTracks& operator=(const AliAnalysisTaskMultPbTracks& task);
  
  ClassDef(AliAnalysisTaskMultPbTracks, 2)


};

#endif
