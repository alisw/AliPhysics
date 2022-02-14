#ifndef AliAnalysisTaskSEImpParResSparse_H
#define AliAnalysisTaskSEImpParResSparse_H

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
#include "AliESDEvent.h"
#include "AliESDtrack.h"


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
  void SetFillSparseForExpert(Bool_t a=kFALSE) {   fFillSparseForExpert=a; 
                                                   if (!a) {
                                                      SetFillSparse_ImpParrphiSparsePtBchargePhi(kFALSE);
                                                      SetFillSparse_ImpParrphiSparsePtzVtxEtaPhi(kFALSE);
                                                      SetFillSparse_ImpParzSparsePtBchargePhi(kFALSE);
                                                      SetFillSparse_ImpParzSparsePtzVtxEtaPhi(kFALSE);
                                                      SetFillSparse_ImpParPullrphiSparsePtBchargePhi(kFALSE);
                                                      SetFillSparse_ImpParPullzSparsePtBchargePhi(kFALSE);
                                                   }
                                               }
  void SetUseTriggerSelection(Bool_t a=kFALSE) {fUseTriggerSelection=a;}

  void SetUsePtWeights(Int_t opt=0, Float_t scaling=1.){ fUseptWeights=opt; fScalingFactPtWeight=scaling;}
  void ConfigurePtWeights();
  
  void SetParticleSpecies(Int_t species) {fParticleSpecies=species;}
  Int_t GetParticleSpecies() const {return fParticleSpecies;}
    
  void SetUsePhysicalPrimaries(Bool_t opt) {fUsePhysicalPrimary=opt;}
  void SetUseGenPt(Bool_t opt) {fUseGeneratedPt=opt;}

  void SetUseFinerPhiBins(Bool_t flag) {fUseFinerPhiBins=flag;}   // mfaggin
  void SetStoreSPDmodulesInfo(Bool_t flag) {fStoreSPDmodulesInfo=flag;}  // mfaggin
  void SetFillSparse_ImpParrphiSparsePtBchargePhi(Bool_t flag) {fFillSparse_ImpParrphiSparsePtBchargePhi=flag;}   // mfaggin
  void SetFillSparse_ImpParrphiSparsePtzVtxEtaPhi(Bool_t flag)   {fFillSparse_ImpParrphiSparsePtzVtxEtaPhi=flag;}  // mfaggin
  void SetFillSparse_ImpParrphiSparsePtEtaPhi(Bool_t flag)   {fFillSparse_ImpParrphiSparsePtEtaPhi=flag;}   // mfaggin
  void SetFillSparse_ImpParPullrphiSparsePtEtaPhi(Bool_t flag) {fFillSparse_ImpParPullrphiSparsePtEtaPhi=flag;}  // mfaggin
  void SetFillSparse_ImpParPullrphiSparsePtBchargePhi(Bool_t flag)   {fFillSparse_ImpParPullrphiSparsePtBchargePhi=flag;} // mfaggin
  void SetFillSparse_ImpParzSparsePtBchargePhi(Bool_t flag) {fFillSparse_ImpParzSparsePtBchargePhi=flag;}   // mfaggin
  void SetFillSparse_ImpParzSparsePtzVtxEtaPhi(Bool_t flag) {fFillSparse_ImpParzSparsePtzVtxEtaPhi=flag;}   // mfaggin
  void SetFillSparse_ImpParzSparsePtEtaPhi(Bool_t flag)  {fFillSparse_ImpParzSparsePtEtaPhi=flag;} // mfaggin
  void SetFillSparse_ImpParPullzSparsePtEtaPhi(Bool_t flag) {fFillSparse_ImpParPullzSparsePtEtaPhi=flag;}   // mfaggin
  void SetFillSparse_ImpParPullzSparsePtBchargePhi(Bool_t flag)   {fFillSparse_ImpParPullzSparsePtBchargePhi=flag;}  // mfaggin

  // establish which tracks label an event as "high multiplicity"
  void SetEstimatorHighMultEv(Bool_t flag)   {fUsetrkITSrefit_highMev=flag;}  // mfaggin

  // establish wheter to maintain primary tracks in MC
  void SetUseOnlyPrimPartMC(Bool_t flag)  {fUseOnlyPrimPartMC=flag;}  // mfaggin

 private:
  
  AliAnalysisTaskSEImpParResSparse(const AliAnalysisTaskSEImpParResSparse &source);
  AliAnalysisTaskSEImpParResSparse& operator=(const AliAnalysisTaskSEImpParResSparse& source);  
  
  Int_t PhiBin(Double_t phi, Bool_t usefinebinsphi=kFALSE) const; // mfaggin (added a new argument)
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

  // mfaggin
  Bool_t fUseFinerPhiBins; /// flag to impose fine phi binning in fImpParrphiSparsePtEtaPhi sparse 
  Bool_t fStoreSPDmodulesInfo; /// flag to decide wheter to store info about SPD modules (layer, detector index)
  Bool_t fFillSparse_ImpParrphiSparsePtBchargePhi;       /// bool to switch on/off the fImpParrphiSparsePtBchargePhi
  Bool_t fFillSparse_ImpParrphiSparsePtzVtxEtaPhi;       /// bool to switch on/off the fImpParrphiSparsePtzVtxEtaPhi
  Bool_t fFillSparse_ImpParrphiSparsePtEtaPhi;           /// bool to switch on/off the fImpParrphiSparsePtEtaPhi
  Bool_t fFillSparse_ImpParPullrphiSparsePtEtaPhi;       /// bool to switch on/off the fImpParPullrphiSparsePtEtaPhi
  Bool_t fFillSparse_ImpParPullrphiSparsePtBchargePhi;   /// bool to switch on/off the fImpParPullrphiSparsePtBchargePhi
  Bool_t fFillSparse_ImpParzSparsePtBchargePhi;          /// bool to switch on/off the fImpParzSparsePtBchargePhi
  Bool_t fFillSparse_ImpParzSparsePtzVtxEtaPhi;          /// bool to switch on/off the fImpParzSparsePtzVtxEtaPhi
  Bool_t fFillSparse_ImpParzSparsePtEtaPhi;              /// bool to switch on/off the fImpParzSparsePtEtaPhi
  Bool_t fFillSparse_ImpParPullzSparsePtEtaPhi;          /// bool to switch on/off the fImpParPullzSparsePtEtaPhi
  Bool_t fFillSparse_ImpParPullzSparsePtBchargePhi;      /// bool to switch on/off the fImpParPullzSparsePtBchargePhi

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
  THnSparseF *fImpParrphiSparsePtEtaPhi_SPDmod; //!<! sparse   (mfaggin)
  TH1F *fPtDistrib; //!<!
  TH1F *fhPtWeights;           // histo with pt weights
  Int_t fUseptWeights;  //0 no weights, 1 pp weights, 2 pPb weights, 3 PbPb weights
  Float_t fScalingFactPtWeight;
  TList *fOutput;  //! 
  Int_t fParticleSpecies; /// particle species
  Bool_t fUsePhysicalPrimary; //
  Bool_t fUseGeneratedPt; //

  // (mfaggin)
  // establish which tracks label an event as "high multiplicity"
  Bool_t fUsetrkITSrefit_highMev;   ///
  // debug histogram for track counting 
  TH1D* fCountTracks;      //!<!
  // debug histogram for counting of number of tracks per event 
  TH1D* fNumTracksperEv;   //!<!
  // debug histogram for counting of number of contributors to the primary vertex
  TH1D* fNumContributors;  //!<!
  // establish wheter to maintain primary tracks in MC
  Bool_t fUseOnlyPrimPartMC;  ///

                                             // mfaggin (number 9 inserted)
  ClassDef(AliAnalysisTaskSEImpParResSparse,10); // AliAnalysisTaskSE for the study of the impact parameter resolution
};

#endif
