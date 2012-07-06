#ifndef ALIPIDCOMBINEDTASK_H
#define ALIPIDCOMBINEDTASK_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#########################################################
//#                                                       # 
//#        Task for testing the combined PID              #
//#                                                       #
//#  Pietro Antonioli, INFN / Pietro.Antonioli@bo.infn.it #
//#  Jens Wiechula, Uni Tübingen / Jens.Wiechula@cern.ch  #
//#                                                       #
//#########################################################

#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliPIDCombined.h>

#include <AliESDtrackCuts.h>
#include <AliAnalysisFilter.h>

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2D;
class fPIDResponse;
class fPIDCombined;

class AliAnalysisTaskPIDCombined : public AliAnalysisTaskSE {
  
public:
  static const Int_t kPtBins = 6;


  AliAnalysisTaskPIDCombined();
  AliAnalysisTaskPIDCombined(const char *name);
  virtual ~AliAnalysisTaskPIDCombined(){;}
  
  virtual void  UserExec(Option_t *option);
  virtual void  UserCreateOutputObjects();
    

private:

  TList fHistList;                   //! list of histograms
  TH2D *fProbTPCnSigma[AliPID::kSPECIES];    //! probabilities vs nSigma in the TPC
  TH2D *fProbTOFnSigma[AliPID::kSPECIES];    //! probabilities vs nSigma  the TOF
  TH2D *fProbTPCTOFnSigmaTPC[AliPID::kSPECIES];    //! comb. probabilities vs nSigma TPC
  TH2D *fProbTPC[AliPID::kSPECIES];          //! probabilities vs mom in the TPC
  TH2D *fProbTOF[AliPID::kSPECIES];          //! probabilities vs mom in the TOF
  TH2D *fProbTPCTOF[AliPID::kSPECIES];       //! combined probabilities vs mom TPC-TOF
  TH1F *fPriors[AliPID::kSPECIES];           //! priors

  TH2D *fProbTPCTOFnSigTPCMom[kPtBins][AliPID::kSPECIES];  // prob. x mom. bins
  TH2D *fProbTPCnSigTPCMom[kPtBins][AliPID::kSPECIES];     // prob. x mom. bins
  TH2D *fProbTOFnSigTOFMom[kPtBins][AliPID::kSPECIES];     // prob. x mom. bins

  TH2D *fPriorsUsed[AliPID::kSPECIES];       //! priors used

  const AliPIDResponse *fPIDResponse;     //! PID response object
  AliPIDCombined       *fPIDCombined;     //! combined PID object
  AliESDtrackCuts      *fTrackCuts;            //! track selection
  AliAnalysisFilter    *fTrackFilter;         //! track filter

  TH2D *fDeDx;                              //! histo with the dedx
  TH2D *fDeDxTuned;                         //! histo to check the dedx tuning in MC


  AliAnalysisTaskPIDCombined(const AliAnalysisTaskPIDCombined &c);
  AliAnalysisTaskPIDCombined& operator= (const AliAnalysisTaskPIDCombined &c);

  void FillHistogram(const char* name, Double_t x, Double_t weight=1.);
  void FillHistogram(const char* name, Double_t x, Double_t y, Double_t weight=1.);
  Int_t GetMomBin(Float_t mom);
  static const char* fgkBinMomDesc[kPtBins];
  
  ClassDef(AliAnalysisTaskPIDCombined, 2);
};
#endif
