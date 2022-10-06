/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskMCPredictions_H
#define AliAnalysisTaskMCPredictions_H

class TTree;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliPPVsMultUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMCPredictions : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskMCPredictions();
  AliAnalysisTaskMCPredictions(const char *name, Int_t lNSmallBinning = 1000, Int_t lNLargeBinning = 2000, Int_t lRebinFactor = 1, Int_t lNBBins = 1, Int_t lNNpartBins = 1, Int_t lNEtaBins = 800 );
  virtual ~AliAnalysisTaskMCPredictions();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;
  
  //Check MC type
  Bool_t IsHijing()  const;
  Bool_t IsDPMJet()  const;
  Bool_t IsEPOSLHC() const;
  
  Double_t ComputeDeltaPhi( Double_t phi1, Double_t phi2) const;
  
  void SetSaveTree( Bool_t lOpt = kTRUE ) { fkSaveTree = lOpt; }
  void SetSelectINELgtZERO ( Bool_t lOpt ) { fkSelectINELgtZERO = lOpt; }
  void SetWideRapidityCut ( Bool_t lOpt = kTRUE) { fkWideRapiditySpeciesStudy = lOpt; }
  
  void SetDoImpactParameterStudy ( Bool_t lOpt = kTRUE) { fkDoImpactParameterStudy = lOpt; }
  void SetDoNpartStudy ( Bool_t lOpt = kTRUE) { fkDoNpartStudy = lOpt; }
  void SetDoNMPIStudy ( Bool_t lOpt = kTRUE) { fkDoNMPIStudy = lOpt; }
  void SetDoRapidityStudy ( Bool_t lOpt = kTRUE) { fkDoRapidityStudy = lOpt; }
  
  void SetMinimumMultiplicity ( Long_t lMinMult ) { fkMinimumMultiplicity = lMinMult; } ;
  void SetCheckOriginThirdArgument ( Bool_t lVal ) { fCheckOriginThirdArgument = lVal; } ;
  void SetNSpecies ( Int_t lVal ) { fkNSpecies = lVal; } ;
  void SetSpeciesSwitch ( Int_t lIndex, Bool_t lVal ) { fkSpeciesSwitch[lIndex] = lVal; } ;
  void SetSpeciesSwitch ( Bool_t *lVal ) { for (Int_t ii=0; ii<76; ii++) fkSpeciesSwitch[ii] = lVal[ii]; } ;
  void SetStandardSpecies ();
  
  //configure intervals
  void ClearEtaIntervals() { fkNIntervals = 0; };
  void AddEtaInterval( Float_t lMin, Float_t lMax, Float_t lWeight = 1.0 ) {
    fkIntervalMinEta[fkNIntervals] = lMin;
    fkIntervalMaxEta[fkNIntervals] = lMax;
    fkIntervalWeight[fkNIntervals] = lWeight;
    fkNIntervals++;
  } ;
  void PrintEtaIntervals(); 
  
  //---------------------------------------------------------------------------------------
  
private:
  // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
  // your data member object is created on the worker nodes and streaming is not needed.
  // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  TList  *fListHist;  //! List of Cascade histograms
  TTree *fTree; //! output tree (PoIs)
  Bool_t fkSaveTree;
  
  //Variables for fTree
  Float_t fPt;
  Int_t fPID;
  Int_t fPIDMother[10];
  Int_t fStatus;
  Int_t fStatusMother[10];
  Int_t fNParents;
  Int_t fV0MMultiplicity;
  
  //Histograms (Desired objects in this cross-checking task)
  TH1D *fHistEventCounter; // histogram for event counting
  TH1D *fHistChargedEta; // histogram for event counting
  
  Int_t fSmallMultRange;
  Int_t fLargeMultRange;
  Int_t fRebinFactor;
  
  Int_t fkNBBins;
  Int_t fkNNpartBins;
  Int_t fkNEtaBins;
  Int_t fkNSpecies;
  Bool_t fkSpeciesSwitch[78]; 
  
  Int_t fkNIntervals; //number of eta intervals in "SPD" estimator (def 1)
  Float_t fkIntervalMinEta[10];
  Float_t fkIntervalMaxEta[10];
  Float_t fkIntervalWeight[10]; 
  
  Bool_t fkSelectINELgtZERO;
  Bool_t fkWideRapiditySpeciesStudy; //if true, |y|<4.0 instead of |y|<0.5 for PID
  
  Bool_t fkDoImpactParameterStudy;
  Bool_t fkDoNpartStudy;
  Bool_t fkDoNMPIStudy;
  Bool_t fkDoRapidityStudy;
  Long_t fkMinimumMultiplicity;
  
  Bool_t fCheckOriginThirdArgument;
  
  //Basic Histograms for counting events as a function of V0M percentiles...
  TH1D *fHistV0MMult; //!
  TH1D *fHistSPDMult; //!
  TH2D *fHistNchVsV0MMult; //!
  TH2D *fHistNchVsSPDMult; //!
  TH1D *fHistNpart; //!
  TH2D *fHistNchVsNpart; //!
  TH1D *fHistB; //!
  TH2D *fHistNchVsB; //!
  TH1D *fHistNMPI; //!
  TH2D *fHistNchVsNMPI; //!
  
  TH1D *fHistPt[78];              //! for keeping track of base spectra
  TH1D *fHistEta[78];              //! for keeping track of base spectra
  TH2D *fHistPtVsV0MMult[78];     //! for keeping track of base spectra
  TH2D *fHistPtVsSPDMult[78];     //! for keeping track of base spectra
  TH2D *fHistEtaVsSPDMult[78];    //! for keeping track of base spectra
  TH2D *fHistYVsSPDMult[78];    //! for keeping track of base spectra
  TH2D *fHistPtVsNpart[78];       //! for keeping track of base spectra
  TH2D *fHistPtVsB[78];           //! for keeping track of base spectra
  TH2D *fHistPtVsNMPI[78];       //! for keeping track of base spectra
  
 
  //double-differential analysis
  Long_t fkDDRebin;
  Long_t fkMaxMultDDV0M;
  Long_t fkMaxMultDDSPD;
  TH2D *fHistV0MvsSPD; //! DD studies
  TH2D *fHistDDNch; //! DD studies
  TH2D *fHistDDNMPI; //! DD studies
  TH2D *fHistDDQ2; //! DD studies
  TH2D *fHistDDYield[78]; //! DD studies
  TH2D *fHistDDPt[78]; //! DD studies
  
  AliAnalysisTaskMCPredictions(const AliAnalysisTaskMCPredictions&);            // not implemented
  AliAnalysisTaskMCPredictions& operator=(const AliAnalysisTaskMCPredictions&); // not implemented
  
  ClassDef(AliAnalysisTaskMCPredictions, 1);
};

#endif

