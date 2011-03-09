#ifndef AliAnalysisTaskEMCALPi0PbPb_cxx
#define AliAnalysisTaskEMCALPi0PbPb_cxx

// $Id$

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTreeStream.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>

class TH1F;
class TH2F;
class TH3F;
class TH1I;
class AliESDEvent;
class AliAODEvent;
class AliPHOSGeometry;
class AliESDtrackCuts;
class AliESDCaloCluster;
class AliAODCaloCluster;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALPi0PbPb : public AliAnalysisTaskSE {
 public:
//  AliAnalysisTaskEMCALPi0PbPb();
  AliAnalysisTaskEMCALPi0PbPb(const char *name);
  virtual ~AliAnalysisTaskEMCALPi0PbPb(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  Double_t GetSigmaMax(AliESDCaloCluster * clust, AliESDEvent * evt); 
  Double_t GetSigmaMax(AliAODCaloCluster * clust, AliAODEvent * evt); 
 
 private:
  TList       *fOutputDataESD;  //! container of output ESD histograms
  TH1I *fhNEvents;    //!  Total analysis Event number in Histogram
  Int_t nModuleEMC;

  //Multi Dis. Vs Centrality at different module
  TH2F **  fhCentVsClusMultEMC;   //!  Cluster multi.  Vs Centrality  at EMCal Module
  TH2F **  fhCentVsCellMultEMC;   //!  Cell multi.     Vs Centrality  at EMCal Module

  TH2F * fhCentVsClusMultEMCAll;  //! Cluster multi.  Vs Centrality  at EMCal all module
  TH2F * fhCentVsCellMultEMCAll;  //! Cluster multi.  Vs Centrality  at EMCal all module

  TH3F * fhCentVsInMVsPtEMC;      //! Centrality Vs Invariant mass Vs Pt of two cluster
  TH3F * fhCentVsInMVsPhiEMC;      //! Centrality Vs Invariant mass Vs Pt of two cluster
  TH2F * fhAsyVsPt;               //! Asymmetry Vs Pt of two cluster 

 // EMCal occupancy and cell energy

  TH3F ** fhCentVsColuRowEMC;            //! Column Vs Row Vs Centrality at EMCal Module 
  TH3F ** fhCentVsColuRowEnerEMC;        //! Column Vs Row Vs Centrality at EMCal Module Vs energy weight 

 // Cluster QA and showershape
  TH3F * fhMggPtEr;            //! Invariant mass vs Pt vs Energy ratio  
  TH3F * fhPtEnSg;             //! Pt vs En vs Cluster long axis

  AliAnalysisTaskEMCALPi0PbPb& operator=(const AliAnalysisTaskEMCALPi0PbPb&); // not implemented
 
  ClassDef(AliAnalysisTaskEMCALPi0PbPb, 1); // example of analysis
};

#endif
