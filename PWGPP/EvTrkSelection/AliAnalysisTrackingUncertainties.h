#ifndef ALIANALYSISTRACKINGUNCERTAINTIES_H
#define ALIANALYSISTRACKINGUNCERTAINTIES_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// Analysis task for the systematic study of the uncertainties related to   //
// the tracking and ITS-TPC matching efficiency for different particle      //
// species.                                                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TList;
class AliESDEvent;
class AliMCEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpid;


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "THnSparse.h"

// ITS->TPC matching constants
const int kMaxMatch=5;
const double kMaxChi2 = 200;


class AliAnalysisTrackingUncertainties : public AliAnalysisTaskSE {
 public:
  AliAnalysisTrackingUncertainties(const char *name);
  AliAnalysisTrackingUncertainties();
  virtual ~AliAnalysisTrackingUncertainties() {}
  //
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //
  void           ProcessTrackCutVariation();
  void           ProcessItsTpcMatching();
  void           Match(const AliESDtrack* tr0, const AliESDtrack* tr1,Double_t rotate=0, Int_t nmatch = 0);
  //
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  void           InitializeTrackCutHistograms();
  //


 private:
  //
  void   BinLogAxis(const THn *h, Int_t axisNumber);
  Bool_t IsVertexAccepted(AliESDEvent * esd, Float_t &vertexZ);
  //
  //
  //
  AliESDEvent * fESD;               //! ESD object
  AliESDpid   * fESDpid;            //! basic pid object
  AliAnalysisUtils * fUtils;        //! vertex and event selection classes
  Bool_t        fMCtrue;            // flag if real data or MC is processed
  //
  //
  TList           * fListHist;      //! output list for histograms
  AliESDtrackCuts * fESDtrackCuts;  // cut set which is under study
  //
  // helper variables for ITS->TPC matching
  //
  const AliESDtrack * fMatchTr[kMaxMatch];
  Double_t fMatchChi[kMaxMatch];
  //
  //
  AliAnalysisTrackingUncertainties(const AliAnalysisTrackingUncertainties&); 
  AliAnalysisTrackingUncertainties& operator=(const AliAnalysisTrackingUncertainties&); 

  ClassDef(AliAnalysisTrackingUncertainties, 1); 
};

#endif
