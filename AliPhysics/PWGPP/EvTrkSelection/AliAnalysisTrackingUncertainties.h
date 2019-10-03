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
#include "THn.h"




class AliAnalysisTrackingUncertainties : public AliAnalysisTaskSE {
 public:
  enum{
    kMaxMatch = 5,
    kNumberOfAxes = 5
  };
  // ITS->TPC matching constants
  static const double kMaxChi2;

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
  void           Match(const AliESDtrack* tr0, const AliESDtrack* tr1, Int_t& nmatch, Bool_t excludeMom = kFALSE, Double_t rotate=0);

  //
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  void           InitializeTrackCutHistograms();
  //
  void           ExcludeMomFromChi2ITSTPC(Bool_t ex = kTRUE) { fExcludeMomFromChi2ITSTPC = ex; }


 private:
  enum ESpecies_t{
      kSpecElectron = 0,
      kSpecPion = 1,
      kSpecKaon = 2,
      kSpecProton = 3,
      kUndef = 4,
      kAll = 5
  };
  //
  void   BinLogAxis(const THn *h, Int_t axisNumber);
  Bool_t IsVertexAccepted(AliESDEvent * esd, Float_t &vertexZ);
  ESpecies_t GetPid(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsElectron(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsPion(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsKaon(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  Bool_t IsProton(const AliESDtrack * const tr, Bool_t useTPCTOF = kFALSE) const;
  //
  Bool_t IsConsistentWithPid(Int_t type, const AliESDtrack * const tr);
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
  Bool_t fExcludeMomFromChi2ITSTPC; // ITS->TPC : exclude momentum from matching chi2 calculation  

  //
  //
  AliAnalysisTrackingUncertainties(const AliAnalysisTrackingUncertainties&); 
  AliAnalysisTrackingUncertainties& operator=(const AliAnalysisTrackingUncertainties&); 

  ClassDef(AliAnalysisTrackingUncertainties, 1); 
};

#endif
