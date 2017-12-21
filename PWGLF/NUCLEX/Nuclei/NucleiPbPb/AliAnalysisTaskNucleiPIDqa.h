/// \class AliAnalysisTaskNucleiPIDqa
/// \author Maximiliano Puccio <maximiliano.puccio@cern.ch>, University and INFN Torino
/// \date March 2017

#ifndef __AliAnalysisTaskNucleiPIDqa__
#define __AliAnalysisTaskNucleiPIDqa__

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>
#include <TString.h>
#include "AliEventCuts.h"

class AliPIDResponse;
class TH2F;
class TList;

class AliAnalysisTaskNucleiPIDqa : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskNucleiPIDqa(TString taskname = "NucleiPIDqa");
  virtual ~AliAnalysisTaskNucleiPIDqa();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  AliEventCuts  fEventCut;

  float fNsigmaITS;  ///< Number of sigma ITS
  float fNsigmaTPC;  ///< Number of sigma TPC
  float fNsigmaTOF;  ///< Number of sigma TOF
  int   fITSsignalN; ///< Number of cluster with PID information in ITS
  int   fTPCsignalN; ///< Number of cluster with PID information in TPC
private:
  AliAnalysisTaskNucleiPIDqa (const AliAnalysisTaskNucleiPIDqa &source);
  AliAnalysisTaskNucleiPIDqa &operator=(const AliAnalysisTaskNucleiPIDqa &source);

  TList*          fList;             //!<! List of the output histograms
  AliPIDResponse* fPID;              //!<! ALICE PID framework

  // Data histograms
  TH2F *fITSsignal[2];               //!<! ITS dE/dx for positive and negative particles
  TH2F *fTPCsignal[2];               //!<! TPC dE/dx for positive and negative particles
  TH2F *fTOFsignal[2];               //!<! TOF beta for positive and negative particles
  TH2F *fITSsignalSelected[2][4][5]; //!<! ITS dE/dx for (anti)-matter for all the nuclei w/o nsigma cut on ITS, TPC or TOF pid
  TH2F *fTPCsignalSelected[2][4][5]; //!<! TPC dE/dx for (anti-)matter for all the nuclei w/o nsigma cut on ITS, TPC or TOF pid
  TH2F *fTOFsignalSelected[2][4][5]; //!<! TOF beta for (anti-)matter for all the nuclei w/o nsigma cut on ITS, TPC or TOF pid
  TH2F *fITSnSigmaSelected[2][4][5]; //!<! ITS nsigma for all the nuclei w/o nsigma cut on ITS, TPC or TOF pid
  TH2F *fTPCnSigmaSelected[2][4][5]; //!<! TPC nsigma for all the nuclei w/o nsigma cut on ITS, TPC or TOF pid
  TH2F *fTOFnSigmaSelected[2][4][5]; //!<! TOF nsigma for all the nuclei w/o nsigma cut on ITS, TPC or TOF pid
  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskNucleiPIDqa, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskNucleiPIDqa__) */
