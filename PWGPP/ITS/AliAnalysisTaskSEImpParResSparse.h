#ifndef ALIANALYSISTASKSEIMPPARRESSPARSE_H
#define ALIANALYSISTASKSEIMPPARRESSPARSE_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEImpParRes
// AliAnalysisTaskSE for the study of the track impact parameter resolution
//
// Authors: xianbao.yuan@pd.infn.it, andrea.dainese@pd.infn.it
// 
//*************************************************************************

class TList;
class TH1F;
class TH1D;
class AliTriggerConfiguration;
class AliVTrack;
class AliVVertex;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"
#include <THnSparse.h>


class AliAnalysisTaskSEImpParResSparse : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskSEImpParResSparse();
  AliAnalysisTaskSEImpParResSparse(const char *name);
  virtual ~AliAnalysisTaskSEImpParResSparse();
  
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetReadMC(Bool_t readMC) { fReadMC=readMC; return; }
  void SetIsAOD(Bool_t isAOD) { fIsAOD=isAOD; return; }
  void SetSelectedPdg(Int_t pdg) { fSelectedPdg=pdg; return; }
  void SetUseDiamond(Bool_t use=kFALSE) { fUseDiamond=use; return; }
  void SetUseRecoVertex(Bool_t use=kFALSE) { fUseRecoVertex=use; return; }
  void SetSkipTrack(Bool_t skip=kFALSE) { fSkipTrack=skip; return; }
  void SetMultiplicityRange(Int_t min,Int_t max) { fMinMult=min; fMaxMult=max; }
  void SetCheckSDDIsIn(Int_t check=0) { fCheckSDDIsIn=check; }
  void SetTriggerClass(TString tclass="") { fTriggerClass=tclass; }
  void SetTriggerMask(ULong64_t mask=0) {fTriggerMask=mask;}
  void SetOCDBPath(TString path="") { fOCDBPath=path; }
  void SetESDtrackCuts(AliESDtrackCuts *esdCuts) {fESDtrackCuts=esdCuts;}

  void SetUseCutGeoNcrNcl(Bool_t opt){fUseCutGeoNcrNcl=opt;}
  void ConfigureCutGeoNcrNcl(Double_t dz, Double_t len, Double_t onept, Double_t fncr, Double_t fncl){
     fDeadZoneWidth=dz;  fCutGeoNcrNclLength=len; fCutGeoNcrNclGeom1Pt=onept;
     fCutGeoNcrNclFractionNcr=fncr; fCutGeoNcrNclFractionNcl=fncl;
  }

  Bool_t IsSelectedCentrality(AliESDEvent *esd) const;
  void BinLogAxis(const THnSparseF *h, Int_t axisNumber);
  void BinLogPtAxis(TH1F *h);
  void SetTrackType(Int_t type=0) { fTrackType=type; }
  void SetFillSparseForExpert(Bool_t a=kFALSE) { fFillSparseForExpert=a; }
  void SetUseTriggerSelection(Bool_t a=kFALSE) {fUseTriggerSelection=a;}

  void SetUsePtWeights(Int_t opt=0, Float_t scaling=1.){ fUseptWeights=opt; fScalingFactPtWeight=scaling;}
  void ConfigurePtWeights();
  

 private:
  
  AliAnalysisTaskSEImpParResSparse(const AliAnalysisTaskSEImpParResSparse &source);
  AliAnalysisTaskSEImpParResSparse& operator=(const AliAnalysisTaskSEImpParResSparse& source);  
  
  Int_t PhiBin(Double_t phi) const;
  Int_t ClusterTypeOnITSLayer(AliESDtrack *t,Int_t layer) const;
  Bool_t IsTrackSelected(AliVTrack *t,AliVVertex *v, AliESDtrackCuts *cuts, const AliVEvent* aod) const;
  Bool_t fIsAOD;  // flag to read AOD or ESD (default is ESD)
  Bool_t fReadMC;       // flag used to switch on/off MC reading
  Int_t  fSelectedPdg;  // only for a given particle species (-1 takes all tracks)
  Bool_t fUseDiamond;   // use diamond constraint in primary vertex
  Bool_t fUseRecoVertex;   // use reco vertex also when reading MC
  Bool_t fSkipTrack;    // redo primary vertex for each track
  Int_t  fMinMult; // minimum multiplicity
  Int_t  fMaxMult; // maximum multiplicity
  Int_t  fCheckSDDIsIn; // check for ITSSDD in the trigger cluster: 0 no check; !=0 check from OCDB
  Bool_t fUseTriggerSelection;
  TString      fTriggerClass; // trigger class to be inspected
  ULong64_t           fTriggerMask;            /// trigger mask
  AliTriggerConfiguration *fTrigConfig; // trigger configuration (read from OCDB)
  TString      fOCDBPath; // to the OCDB
  AliESDtrackCuts *fESDtrackCuts; // track cuts 
  TH1F  *fNentries;   //! histogram of number of events
  TH1F  *fMultiplicity; //!
  Bool_t fUseCutGeoNcrNcl; /// flag for enabling/disabling geometrical cut on TPC track
  Double_t fDeadZoneWidth;       /// 1st parameter of GeoNcrNcl cut
  Double_t fCutGeoNcrNclLength;  /// 2nd parameter of GeoNcrNcl cut
  Double_t fCutGeoNcrNclGeom1Pt; /// 3rd parameter of GeoNcrNcl cut
  Double_t fCutGeoNcrNclFractionNcr; /// 4th parameter of GeoNcrNcl cut
  Double_t fCutGeoNcrNclFractionNcl; /// 5th parameter of GeoNcrNcl cut
  Int_t fTrackType;
  Bool_t fFillSparseForExpert;
  THnSparseF *fImpParrphiSparsePtBchargePhi; //!<! sparse
  THnSparseF *fImpParrphiSparsePtzVtxEtaPhi; //!<! sparse
  THnSparseF *fImpParrphiSparsePtEtaPhi; //!<! sparse
  THnSparseF *fImpParPullrphiSparsePtEtaPhi; //!<! sparse
  THnSparseF *fImpParPullrphiSparsePtBchargePhi; //!<! sparse
  THnSparseF *fImpParzSparsePtBchargePhi; //!<! sparse
  THnSparseF *fImpParzSparsePtzVtxEtaPhi; //!<! sparse
  THnSparseF *fImpParzSparsePtEtaPhi; //!<! sparse
  THnSparseF *fImpParPullzSparsePtEtaPhi; //!<! sparse
  THnSparseF *fImpParPullzSparsePtBchargePhi; //!<! sparse
  TH1F *fPtDistrib; //!<!
  TH1F *fhPtWeights;           // histo with pt weights
  Int_t fUseptWeights;  //0 no weights, 1 pp weights, 2 pPb weights, 3 PbPb weights
  Float_t fScalingFactPtWeight;
  TList *fOutput;  //! 

  ClassDef(AliAnalysisTaskSEImpParResSparse,3); // AliAnalysisTaskSE for the study of the impact parameter resolution
};

#endif
