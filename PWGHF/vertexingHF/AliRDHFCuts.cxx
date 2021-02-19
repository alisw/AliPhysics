/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// Base class for cuts on AOD reconstructed heavy-flavour decay
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////
#include <Riostream.h>

#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliAODRecoDecayHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliVertexerTracks.h"
#include "AliTimeRangeMasking.h"
#include "AliEventCuts.h"
#include "AliRDHFCuts.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliAODVZERO.h"
#include "TRandom.h"
#include <TF1.h>
#include <TFile.h>
#include <TKey.h>

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCuts);
/// \endcond


//--------------------------------------------------------------------------
AliRDHFCuts::AliRDHFCuts(const Char_t* name, const Char_t* title) :
AliAnalysisCuts(name,title),
fMinVtxType(3),
fMinVtxContr(1),
fMaxVtxRedChi2(1e6),
fMaxVtxZ(10.),
fMinSPDMultiplicity(0),
fMinContrPileupMV(5),
fMaxVtxChi2PileupMV(5.),
fMinWDzPileupMV(15.),
fRejectPlpFromDiffBCMV(kFALSE),
fTriggerMask(AliVEvent::kAnyINT),
fUseOnlyOneTrigger(kFALSE),
fTrackCuts(0),
fnPtBins(1),
fnPtBinLimits(1),
fPtBinLimits(0),
fnVars(1),
fVarNames(0),
fnVarsForOpt(0),
fVarsForOpt(0),
fGlobalIndex(1),
fCutsRD(0),
fIsUpperCut(0),
fUsePID(kFALSE),
fUseAOD049(kFALSE),
fPidHF(0),
fWhyRejection(0),
fEvRejectionBits(0),
fRemoveDaughtersFromPrimary(kFALSE),
fUseMCVertex(kFALSE),
fUsePhysicsSelection(kTRUE),
fOptPileup(0),
fMinContrPileup(3),
fMinDzPileup(0.6),
fUseMultDepPileupCut(kFALSE),
fUseCentrality(0),
fMinCentrality(0.),
fMaxCentrality(100.),
fMultSelectionObjectName("MultSelection"),
fFixRefs(kFALSE),
fIsSelectedCuts(0),
fIsSelectedPID(0),
fMinPtCand(-1.),
fMaxPtCand(100000.),
fMaxRapidityCand(-999.),
fKeepSignalMC(kFALSE),
fIsCandTrackSPDFirst(kFALSE),
fMaxPtCandTrackSPDFirst(0.),
fApplySPDDeadPbPb2011(kFALSE),
fApplySPDMisalignedPP2012(kFALSE),
fApplySPDUniformAccPbPbRun2(kFALSE),
fMaxDiffTRKV0Centr(-1.),
fRemoveTrackletOutliers(kFALSE),
fCutOnzVertexSPD(3),
fKinkReject(kFALSE),
fUseTrackSelectionWithFilterBits(kTRUE),
fUseCentrFlatteningInMC(kFALSE),
fHistCentrDistr(0x0),
fCutRatioClsOverCrossRowsTPC(0),
fCutRatioSignalNOverCrossRowsTPC(0),
fCutTPCSignalN(0),
fCutMinCrossedRowsTPCPtDep(""),
f1CutMinNCrossedRowsTPCPtDep(0x0),
fUseCutGeoNcrNcl(kFALSE),
fDeadZoneWidth(3.),
fCutGeoNcrNclLength(130.),
fCutGeoNcrNclGeom1Pt(1.5),
fCutGeoNcrNclFractionNcr(0.85),
fCutGeoNcrNclFractionNcl(0.7),
fUseV0ANDSelectionOffline(kFALSE),
fUseTPCtrackCutsOnThisDaughter(kTRUE),
fApplyZcutOnSPDvtx(kFALSE),
fUsePreselect(0),
fAliEventCuts(0x0),
fApplyCentralityCorrCuts(kFALSE),
fApplyPbPbOutOfBunchPileupCuts(0),
fApplyPbPbOutOfBunchPileupCutsITSTPC(0),
fKeepOnlyPbPbOutOfBunchPileupCutsITSTPC(kFALSE),
fUseAliEventCuts(kFALSE),
fUseTimeRangeCutForPbPb2018(kTRUE),
fTimeRangeCut(),
fCurrentRun(-1),
fEnableNsigmaTPCDataCorr(kFALSE),
fSystemForNsigmaTPCDataCorr(AliAODPidHF::kNone)
{
  //
  // Default Constructor
  //
  fAliEventCuts = new AliEventCuts();
  fAliEventCuts->SetManualMode();
  fTriggerClass[0]="CINT1"; fTriggerClass[1]="";
}
//--------------------------------------------------------------------------
AliRDHFCuts::AliRDHFCuts(const AliRDHFCuts &source) :
  AliAnalysisCuts(source),
  fMinVtxType(source.fMinVtxType),
  fMinVtxContr(source.fMinVtxContr),
  fMaxVtxRedChi2(source.fMaxVtxRedChi2),
  fMaxVtxZ(source.fMaxVtxZ),
  fMinSPDMultiplicity(source.fMinSPDMultiplicity),
  fMinContrPileupMV(source.fMinContrPileupMV),
  fMaxVtxChi2PileupMV(source.fMaxVtxChi2PileupMV),
  fMinWDzPileupMV(source.fMinWDzPileupMV),
  fRejectPlpFromDiffBCMV(source.fRejectPlpFromDiffBCMV),
  fTriggerMask(source.fTriggerMask),
  fUseOnlyOneTrigger(source.fUseOnlyOneTrigger),
  fTriggerClass(),
  fTrackCuts(0),
  fnPtBins(source.fnPtBins),
  fnPtBinLimits(source.fnPtBinLimits),
  fPtBinLimits(0),
  fnVars(source.fnVars),
  fVarNames(0),
  fnVarsForOpt(source.fnVarsForOpt),
  fVarsForOpt(0),
  fGlobalIndex(source.fGlobalIndex),
  fCutsRD(0),
  fIsUpperCut(0),
  fUsePID(source.fUsePID),
  fUseAOD049(source.fUseAOD049),
  fPidHF(0),
  fWhyRejection(source.fWhyRejection),
  fEvRejectionBits(source.fEvRejectionBits),
  fRemoveDaughtersFromPrimary(source.fRemoveDaughtersFromPrimary),
  fUseMCVertex(source.fUseMCVertex),
  fUsePhysicsSelection(source.fUsePhysicsSelection),
  fOptPileup(source.fOptPileup),
  fMinContrPileup(source.fMinContrPileup),
  fMinDzPileup(source.fMinDzPileup),
  fUseMultDepPileupCut(source.fUseMultDepPileupCut),
  fUseCentrality(source.fUseCentrality),
  fMinCentrality(source.fMinCentrality),
  fMaxCentrality(source.fMaxCentrality),
  fMultSelectionObjectName(source.fMultSelectionObjectName),
  fFixRefs(source.fFixRefs),
  fIsSelectedCuts(source.fIsSelectedCuts),
  fIsSelectedPID(source.fIsSelectedPID),
  fMinPtCand(source.fMinPtCand),
  fMaxPtCand(source.fMaxPtCand),
  fMaxRapidityCand(source.fMaxRapidityCand),
  fKeepSignalMC(source.fKeepSignalMC),
  fIsCandTrackSPDFirst(source.fIsCandTrackSPDFirst),
  fMaxPtCandTrackSPDFirst(source.fMaxPtCandTrackSPDFirst),
  fApplySPDDeadPbPb2011(source.fApplySPDDeadPbPb2011),
  fApplySPDMisalignedPP2012(source.fApplySPDMisalignedPP2012),
  fApplySPDUniformAccPbPbRun2(source.fApplySPDUniformAccPbPbRun2),
  fMaxDiffTRKV0Centr(source.fMaxDiffTRKV0Centr),
  fRemoveTrackletOutliers(source.fRemoveTrackletOutliers),
  fCutOnzVertexSPD(source.fCutOnzVertexSPD),
  fKinkReject(source.fKinkReject),
  fUseTrackSelectionWithFilterBits(source.fUseTrackSelectionWithFilterBits),
  fUseCentrFlatteningInMC(source.fUseCentrFlatteningInMC),
  fHistCentrDistr(0x0),
  fCutRatioClsOverCrossRowsTPC(source.fCutRatioClsOverCrossRowsTPC),
  fCutRatioSignalNOverCrossRowsTPC(source.fCutRatioSignalNOverCrossRowsTPC),
  fCutTPCSignalN(source.fCutTPCSignalN),
  fCutMinCrossedRowsTPCPtDep(""),
  f1CutMinNCrossedRowsTPCPtDep(0x0),
  fUseCutGeoNcrNcl(source.fUseCutGeoNcrNcl),
  fDeadZoneWidth(source.fDeadZoneWidth),
  fCutGeoNcrNclLength(source.fCutGeoNcrNclLength),
  fCutGeoNcrNclGeom1Pt(source.fCutGeoNcrNclGeom1Pt),
  fCutGeoNcrNclFractionNcr(source.fCutGeoNcrNclFractionNcr),
  fCutGeoNcrNclFractionNcl(source.fCutGeoNcrNclFractionNcl),
  fUseV0ANDSelectionOffline(source.fUseV0ANDSelectionOffline),
  fUseTPCtrackCutsOnThisDaughter(source.fUseTPCtrackCutsOnThisDaughter),
  fApplyZcutOnSPDvtx(source.fApplyZcutOnSPDvtx),
  fUsePreselect(source.fUsePreselect),
  fAliEventCuts(source.fAliEventCuts),
  fApplyCentralityCorrCuts(source.fApplyCentralityCorrCuts),
  fApplyPbPbOutOfBunchPileupCuts(source.fApplyPbPbOutOfBunchPileupCuts),
  fApplyPbPbOutOfBunchPileupCutsITSTPC(source.fApplyPbPbOutOfBunchPileupCutsITSTPC),
  fKeepOnlyPbPbOutOfBunchPileupCutsITSTPC(source.fKeepOnlyPbPbOutOfBunchPileupCutsITSTPC),
  fUseAliEventCuts(source.fUseAliEventCuts),
  fUseTimeRangeCutForPbPb2018(source.fUseTimeRangeCutForPbPb2018),
  fTimeRangeCut(),
  fCurrentRun(source.fCurrentRun),
  fEnableNsigmaTPCDataCorr(source.fEnableNsigmaTPCDataCorr),
  fSystemForNsigmaTPCDataCorr(source.fSystemForNsigmaTPCDataCorr)
{
  //
  // Copy constructor
  //
  cout<<"Copy constructor"<<endl;
  fTriggerClass[0] = source.fTriggerClass[0];
  fTriggerClass[1] = source.fTriggerClass[1];
  if(source.GetTrackCuts()) AddTrackCuts(source.GetTrackCuts());
  if(source.fPtBinLimits) SetPtBins(source.fnPtBinLimits,source.fPtBinLimits);
  if(source.fVarNames) SetVarNames(source.fnVars,source.fVarNames,source.fIsUpperCut);
  if(source.fCutsRD) SetCuts(source.fGlobalIndex,source.fCutsRD);
  if(source.fVarsForOpt) SetVarsForOpt(source.fnVarsForOpt,source.fVarsForOpt);
  if(source.fPidHF) SetPidHF(source.fPidHF);
  if(source.fHistCentrDistr) fHistCentrDistr=(TH1F*)(source.fHistCentrDistr->Clone());
  if(source.fCutMinCrossedRowsTPCPtDep) fCutMinCrossedRowsTPCPtDep=source.fCutMinCrossedRowsTPCPtDep;
  if(source.f1CutMinNCrossedRowsTPCPtDep) f1CutMinNCrossedRowsTPCPtDep=new TFormula(*(source.f1CutMinNCrossedRowsTPCPtDep));
  PrintAll();

}
//--------------------------------------------------------------------------
AliRDHFCuts &AliRDHFCuts::operator=(const AliRDHFCuts &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAnalysisCuts::operator=(source);

  fMinVtxType=source.fMinVtxType;
  fMinVtxContr=source.fMinVtxContr;
  fMaxVtxRedChi2=source.fMaxVtxRedChi2;
  fMaxVtxZ=source.fMaxVtxZ;
  fMinSPDMultiplicity=source.fMinSPDMultiplicity;
  fMinContrPileupMV=source.fMinContrPileupMV;
  fMaxVtxChi2PileupMV=source.fMaxVtxChi2PileupMV;
  fMinWDzPileupMV=source.fMinWDzPileupMV;
  fRejectPlpFromDiffBCMV=source.fRejectPlpFromDiffBCMV;
  fTriggerMask=source.fTriggerMask;
  fUseOnlyOneTrigger=source.fUseOnlyOneTrigger;
  fTriggerClass[0]=source.fTriggerClass[0];
  fTriggerClass[1]=source.fTriggerClass[1];
  fnPtBins=source.fnPtBins;
  fnPtBinLimits=source.fnPtBinLimits;
  fnVars=source.fnVars;
  fGlobalIndex=source.fGlobalIndex;
  fnVarsForOpt=source.fnVarsForOpt;
  fUsePID=source.fUsePID;
  fUseAOD049=source.fUseAOD049;
  if(fPidHF) delete fPidHF;
  fPidHF=new AliAODPidHF(*(source.GetPidHF()));
  fWhyRejection=source.fWhyRejection;
  fEvRejectionBits=source.fEvRejectionBits;
  fRemoveDaughtersFromPrimary=source.fRemoveDaughtersFromPrimary;
  fUseMCVertex=source.fUseMCVertex;
  fUsePhysicsSelection=source.fUsePhysicsSelection;
  fOptPileup=source.fOptPileup;
  fMinContrPileup=source.fMinContrPileup;
  fMinDzPileup=source.fMinDzPileup;
  fUseMultDepPileupCut=source.fUseMultDepPileupCut;
  fUseCentrality=source.fUseCentrality;
  fMinCentrality=source.fMinCentrality;
  fMaxCentrality=source.fMaxCentrality;
  fMultSelectionObjectName=source.fMultSelectionObjectName;
  fFixRefs=source.fFixRefs;
  fIsSelectedCuts=source.fIsSelectedCuts;
  fIsSelectedPID=source.fIsSelectedPID;
  fMinPtCand=source.fMinPtCand;
  fMaxPtCand=source.fMaxPtCand;
  fMaxRapidityCand=source.fMaxRapidityCand;
  fKeepSignalMC=source.fKeepSignalMC;
  fIsCandTrackSPDFirst=source.fIsCandTrackSPDFirst;
  fMaxPtCandTrackSPDFirst=source.fMaxPtCandTrackSPDFirst;
  fApplySPDDeadPbPb2011=source.fApplySPDDeadPbPb2011;
  fApplySPDUniformAccPbPbRun2=source.fApplySPDUniformAccPbPbRun2;
  fApplySPDMisalignedPP2012=source.fApplySPDMisalignedPP2012;
  fMaxDiffTRKV0Centr=source.fMaxDiffTRKV0Centr;
  fRemoveTrackletOutliers=source.fRemoveTrackletOutliers;
  fCutOnzVertexSPD=source.fCutOnzVertexSPD;
  fKinkReject=source.fKinkReject;
  fUseTrackSelectionWithFilterBits=source.fUseTrackSelectionWithFilterBits;
  if(fHistCentrDistr) delete fHistCentrDistr;
  fUseCentrFlatteningInMC=source.fUseCentrFlatteningInMC;
  if(source.fHistCentrDistr)fHistCentrDistr=(TH1F*)(source.fHistCentrDistr->Clone());

  if(source.GetTrackCuts()) {delete fTrackCuts; fTrackCuts=new AliESDtrackCuts(*(source.GetTrackCuts()));}
  if(source.fPtBinLimits) SetPtBins(source.fnPtBinLimits,source.fPtBinLimits);
  if(source.fVarNames) SetVarNames(source.fnVars,source.fVarNames,source.fIsUpperCut);
  if(source.fCutsRD) SetCuts(source.fGlobalIndex,source.fCutsRD);
  if(source.fVarsForOpt) SetVarsForOpt(source.fnVarsForOpt,source.fVarsForOpt);

  if(fCutMinCrossedRowsTPCPtDep) fCutMinCrossedRowsTPCPtDep=source.fCutMinCrossedRowsTPCPtDep;
  if(f1CutMinNCrossedRowsTPCPtDep) delete f1CutMinNCrossedRowsTPCPtDep;
  if(source.f1CutMinNCrossedRowsTPCPtDep) f1CutMinNCrossedRowsTPCPtDep=new TFormula(*(source.f1CutMinNCrossedRowsTPCPtDep));
  fCutRatioClsOverCrossRowsTPC=source.fCutRatioClsOverCrossRowsTPC;
  fCutRatioSignalNOverCrossRowsTPC=source.fCutRatioSignalNOverCrossRowsTPC;
  fCutTPCSignalN=source.fCutTPCSignalN;
  fUseCutGeoNcrNcl=source.fUseCutGeoNcrNcl;
  fDeadZoneWidth=source.fDeadZoneWidth;
  fCutGeoNcrNclLength=source.fCutGeoNcrNclLength;
  fCutGeoNcrNclGeom1Pt=source.fCutGeoNcrNclGeom1Pt;
  fCutGeoNcrNclFractionNcr=source.fCutGeoNcrNclFractionNcr;
  fCutGeoNcrNclFractionNcl=source.fCutGeoNcrNclFractionNcl;
  fUseV0ANDSelectionOffline=source.fUseV0ANDSelectionOffline;
  fUseTPCtrackCutsOnThisDaughter=source.fUseTPCtrackCutsOnThisDaughter;
  fUsePreselect=source.fUsePreselect;
  fAliEventCuts=source.fAliEventCuts;
  fApplyCentralityCorrCuts=source.fApplyCentralityCorrCuts;
  fApplyPbPbOutOfBunchPileupCuts=source.fApplyPbPbOutOfBunchPileupCuts;
  fApplyPbPbOutOfBunchPileupCutsITSTPC=source.fApplyPbPbOutOfBunchPileupCutsITSTPC;
  fKeepOnlyPbPbOutOfBunchPileupCutsITSTPC=source.fKeepOnlyPbPbOutOfBunchPileupCutsITSTPC;
  fUseAliEventCuts=source.fUseAliEventCuts;
  fUseTimeRangeCutForPbPb2018=source.fUseTimeRangeCutForPbPb2018;
  fCurrentRun=source.fCurrentRun;
  fEnableNsigmaTPCDataCorr=source.fEnableNsigmaTPCDataCorr;
  fSystemForNsigmaTPCDataCorr=source.fSystemForNsigmaTPCDataCorr;

  PrintAll();

  return *this;
}
//--------------------------------------------------------------------------
AliRDHFCuts::~AliRDHFCuts() {
  //
  // Default Destructor
  //
  if(fTrackCuts) { delete fTrackCuts; fTrackCuts=0; }
  if(fPtBinLimits) {delete [] fPtBinLimits; fPtBinLimits=0;}
  if(fVarNames) {delete [] fVarNames; fVarNames=0;}
  if(fVarsForOpt) {delete [] fVarsForOpt; fVarsForOpt=0;}
  if(fCutsRD) {
    delete [] fCutsRD;
    fCutsRD=0;
  }
  if(fIsUpperCut) {delete [] fIsUpperCut; fIsUpperCut=0;}
  if(fPidHF){
    delete fPidHF;
    fPidHF=0;
  }
  if(fHistCentrDistr)delete fHistCentrDistr;

  if(f1CutMinNCrossedRowsTPCPtDep) {
    delete f1CutMinNCrossedRowsTPCPtDep;
    f1CutMinNCrossedRowsTPCPtDep = 0;
  }
  delete fAliEventCuts;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCuts::IsEventSelectedInCentrality(AliVEvent *event) {
  //
  // Centrality selection
  //
  if(fUseCentrality<kCentOff||fUseCentrality>=kCentInvalid){
    AliWarning("Centrality estimator not valid");
    return 3;
  }else{
    Float_t centvalue=GetCentrality((AliAODEvent*)event);
    if (centvalue<-998.){//-999 if no centralityP
      return 3;
    }else if(fEvRejectionBits&(1<<kMismatchOldNewCentrality)){
      return 3;
    }else{
      if (centvalue<fMinCentrality || centvalue>fMaxCentrality){
	return 2;
      }
      // selection to flatten the centrality distribution
      if(fHistCentrDistr){
	if(!IsEventSelectedForCentrFlattening(centvalue))return 4;
      }
    }
  }
  return 0;
}


//-------------------------------------------------
void AliRDHFCuts::SetHistoForCentralityFlattening(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef,Int_t switchTRand){
  // set the histo for centrality flattening
  // the centrality is flatten in the range minCentr,maxCentr
  // if centrRef is zero, the minimum in h within (minCentr,maxCentr) defines the reference
  //                positive, the value of h(centrRef) defines the reference (-> the centrality distribution might be not flat in the whole desired range)
  //                negative, h(bin with max in range)*centrRef is used to define the reference (-> defines the maximum loss of events, also in this case the distribution might be not flat)
  // switchTRand is used to set the unerflow bin of the histo: if it is < -1 in the analysis the random event selection will be done on using TRandom

  if(maxCentr<minCentr){
    AliWarning("AliRDHFCuts::Wrong centralities values while setting the histogram for centrality flattening");
  }

  if(fHistCentrDistr)delete fHistCentrDistr;
  fHistCentrDistr=(TH1F*)h->Clone("hCentralityFlat");
  fHistCentrDistr->SetTitle("Reference histo for centrality flattening");
  Int_t minbin=fHistCentrDistr->FindBin(minCentr*1.00001); // fast if fix bin width
  Int_t maxbin=fHistCentrDistr->FindBin(maxCentr*0.9999);
  fHistCentrDistr->GetXaxis()->SetRange(minbin,maxbin);
  Double_t ref=0.,bincont=0.,binrefwidth=1.;
  Int_t binref=0;
  if(TMath::Abs(centrRef)<0.0001){
    binref=fHistCentrDistr->GetMinimumBin();
    binrefwidth=fHistCentrDistr->GetBinWidth(binref);
    ref=fHistCentrDistr->GetBinContent(binref)/binrefwidth;
  }
  else if(centrRef>0.){
    binref=h->FindBin(centrRef);
    if(binref<1||binref>h->GetNbinsX()){
      AliWarning("AliRDHFCuts::Wrong centrality reference value while setting the histogram for centrality flattening");
    }
    binrefwidth=fHistCentrDistr->GetBinWidth(binref);
    ref=fHistCentrDistr->GetBinContent(binref)/binrefwidth;
  }
  else{
    if(centrRef<-1) AliWarning("AliRDHFCuts: with this centrality reference no flattening will be applied");
    binref=fHistCentrDistr->GetMaximumBin();
    binrefwidth=fHistCentrDistr->GetBinWidth(binref);
    ref=fHistCentrDistr->GetMaximum()*TMath::Abs(centrRef)/binrefwidth;
  }

  for(Int_t j=1;j<=h->GetNbinsX();j++){// Now set the "probabilities"
    if(h->GetBinLowEdge(j)*1.0001>=minCentr&&h->GetBinLowEdge(j+1)*0.9999<=maxCentr){
      bincont=h->GetBinContent(j);
      fHistCentrDistr->SetBinContent(j,ref/bincont*h->GetBinWidth(j));
      fHistCentrDistr->SetBinError(j,h->GetBinError(j)*ref/bincont);
    }
    else{
      h->SetBinContent(j,1.1);// prob > 1 to assure that events will not be rejected
    }
  }

  fHistCentrDistr->SetBinContent(0,switchTRand);
  return;

}

//-------------------------------------------------
Bool_t AliRDHFCuts::IsEventSelectedForCentrFlattening(Float_t centvalue){
  //
  //  Random event selection, based on fHistCentrDistr, to flatten the centrality distribution
  //  Can be faster if it was required that fHistCentrDistr covers
  //  exactly the desired centrality range (e.g. part of the lines below should be done during the
  // setting of the histo) and TH1::SetMinimum called
  //

  if(!fHistCentrDistr) return kTRUE;
  // Int_t maxbin=fHistCentrDistr->FindBin(fMaxCentrality*0.9999);
  //   if(maxbin>fHistCentrDistr->GetNbinsX()){
  //     AliWarning("AliRDHFCuts: The maximum centrality exceeds the x-axis limit of the histogram for centrality flattening");
  //   }

  Int_t bin=fHistCentrDistr->FindBin(centvalue); // Fast if the histo has a fix bin
  Double_t bincont=fHistCentrDistr->GetBinContent(bin);
  Double_t centDigits=centvalue-(Int_t)(centvalue*100.)/100.;// this is to extract a random number between 0 and 0.01

  if(fHistCentrDistr->GetBinContent(0)<-0.9999){
    if(gRandom->Uniform(1.)<bincont)return kTRUE;
    return kFALSE;
  }

  if(centDigits*100.<bincont)return kTRUE;
  return kFALSE;

}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetupPID(AliVEvent *event) {
  // Set the PID response object in the AliAODPidHF
  // in case of old PID sets the TPC dE/dx BB parameterization
  
  Bool_t isMC=kFALSE;
  if(fPidHF){
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    if(fPidHF->GetPidResponse()==0x0){
      AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
      fPidHF->SetPidResponse(pidResp);
    }
    if(fPidHF->GetUseCombined()) fPidHF->SetUpCombinedPID();
    if(fPidHF->GetOldPid()) {

      TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if(mcArray) {isMC=kTRUE;fUseAOD049=kFALSE;}

      // pp, from LHC10d onwards
      if((event->GetRunNumber()>121693 && event->GetRunNumber()<136851) ||
	 event->GetRunNumber()>139517) fPidHF->SetOnePad(kTRUE);
      // pp, 2011 low energy run
      if((event->GetRunNumber()>=146686 && event->GetRunNumber()<=146860)){
	fPidHF->SetppLowEn2011(kTRUE);
	fPidHF->SetOnePad(kFALSE);
      }
      // PbPb LHC10h
      if(event->GetRunNumber()>=136851 && event->GetRunNumber()<=139517) fPidHF->SetPbPb(kTRUE);
      // MC
      if(isMC) fPidHF->SetMC(kTRUE);
      if(isMC && (event->GetRunNumber()>=146686 && event->GetRunNumber()<=146860))
	fPidHF->SetMClowenpp2011(kTRUE);
      fPidHF->SetBetheBloch();
    }else{
      // check that AliPIDResponse object was properly set in case of using OADB
      if(fPidHF->GetPidResponse()==0x0) AliFatal("AliPIDResponse object not set");
    }

    // force recomputation of TOF Nsigma with tune-on-data to have latest development of tail parametrisation in old AODs
    if(isMC) {
      for(Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++) {
        AliVTrack* track=dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
        if(!track || track->GetTOFsignalTunedOnData() > 99999) continue;
        track->SetTOFsignalTunedOnData(100000);
      }
    }
    else { // apply TPC postcalibration
      if(fEnableNsigmaTPCDataCorr) {

        Bool_t isPass1 = kFALSE;
        TTree *treeAOD = inputHandler->GetTree();
        TString currentFile = treeAOD->GetCurrentFile()->GetName();
        if((currentFile.Contains("LHC18q") || currentFile.Contains("LHC18r")) && currentFile.Contains("pass1"))
          isPass1 = kTRUE;

        fPidHF->EnableNsigmaTPCDataCorr(event->GetRunNumber(),fSystemForNsigmaTPCDataCorr,isPass1);
      }
    }
  }
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::IsEventSelected(AliVEvent *event) {
  //
  // Event selection
  //
  //if(fTriggerMask && event->GetTriggerMask()!=fTriggerMask) return kFALSE;

  // commented for the time being
  if(fUseAliEventCuts) return IsEventSelectedWithAliEventCuts(event);

  fWhyRejection=0;
  fEvRejectionBits=0;
  Bool_t accept=kTRUE;


  // check if it's MC
  Bool_t isMC=kFALSE;
  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(mcArray) {isMC=kTRUE;fUseAOD049=kFALSE;}


  SetupPID(event);

  // trigger class
  TString firedTriggerClasses=((AliAODEvent*)event)->GetFiredTriggerClasses();
  // don't do for MC and for PbPb 2010 data
  if(!isMC && (event->GetRunNumber()<136851 || event->GetRunNumber()>139517)) {
    if(!firedTriggerClasses.Contains(fTriggerClass[0].Data()) &&
       (fTriggerClass[1].CompareTo("")==0 || !firedTriggerClasses.Contains(fTriggerClass[1].Data())) ) {
      fWhyRejection=5;
      fEvRejectionBits+=1<<kNotSelTrigger;
      accept=kFALSE;
    }
  }


  // TEMPORARY FIX FOR GetEvent
  Int_t nTracks=((AliAODEvent*)event)->GetNumberOfTracks();
  for(Int_t itr=0; itr<nTracks; itr++){
    AliAODTrack* tr=(AliAODTrack*)((AliAODEvent*)event)->GetTrack(itr);
    tr->SetAODEvent((AliAODEvent*)event);
  }

  // TEMPORARY FIX FOR REFERENCES
  // Fix references to daughter tracks
  //  if(fFixRefs) {
  //    AliAnalysisVertexingHF *fixer = new AliAnalysisVertexingHF();
  //    fixer->FixReferences((AliAODEvent*)event);
  //    delete fixer;
  //  }
  //

  // cuts used for run 2 Pb-Pb from AliEventCuts
  Bool_t doAliEvCuts=kFALSE;
  if(fApplyCentralityCorrCuts || fApplyPbPbOutOfBunchPileupCuts>0 || fApplyPbPbOutOfBunchPileupCutsITSTPC){
    doAliEvCuts=kTRUE;
    Int_t runNumb=event->GetRunNumber();
    if(runNumb >= 244917 && runNumb <= 246994) fAliEventCuts->SetupRun2PbPb();
    else if(runNumb >= 295369 && runNumb <= 297624) fAliEventCuts->SetupPbPb2018();
    else doAliEvCuts=kFALSE;
    if(fApplyPbPbOutOfBunchPileupCutsITSTPC)
      fAliEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(true, fApplyPbPbOutOfBunchPileupCutsITSTPC);
  }

  if(doAliEvCuts) fAliEventCuts->AcceptEvent(event);

  // physics selection requirements
  if(fUsePhysicsSelection){
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if(!isSelected) {
      if(accept) fWhyRejection=7;
      fEvRejectionBits+=1<<kPhysicsSelection;
      accept=kFALSE;
    }else{
      if(fUseOnlyOneTrigger){
	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()!=fTriggerMask){
	  if(accept) fWhyRejection=7;
	  fEvRejectionBits+=1<<kPhysicsSelection;
	  accept=kFALSE;
	}
      }
      if(fUseV0ANDSelectionOffline){
	AliAODVZERO* v0data=(AliAODVZERO*)((AliAODEvent*)event)->GetVZEROData();
	Int_t tv0a=v0data->GetV0ADecision();
	Int_t tv0c=v0data->GetV0CDecision();
	if(!(tv0a==1 && tv0c==1)){
	  if(accept) fWhyRejection=7;
	  fEvRejectionBits+=1<<kPhysicsSelection;
	  accept=kFALSE;
	}
      }
    }
  }

  if(fUseTimeRangeCutForPbPb2018){
    Int_t nrun=event->GetRunNumber();
    if(nrun!=fCurrentRun){
      fCurrentRun=nrun;
      fTimeRangeCut.InitFromRunNumber(fCurrentRun);
    }
    if(fTimeRangeCut.CutEvent((AliAODEvent*)event)){
      // use same fWhyRejection as for physics selection, to have proper counting of events for norm
      if(accept) fWhyRejection=7;
      fEvRejectionBits+=1<<kBadTimeRange;
      accept=kFALSE;
    }
  }
  
  // centrality selection
  if (fUseCentrality!=kCentOff) {
    Int_t rejection=IsEventSelectedInCentrality(event);
    Bool_t okCent=kFALSE;
    if(rejection==0) okCent=kTRUE;
    if(isMC && rejection==4 && !fUseCentrFlatteningInMC) okCent=kTRUE;
    if(!okCent){
      if(accept) fWhyRejection=rejection;
      if(fWhyRejection==4)fEvRejectionBits+=1<<kCentralityFlattening;
      else fEvRejectionBits+=1<<kOutsideCentrality;
      accept=kFALSE;
    }

  }

  // PbPb2011 outliers in tracklets vs. VZERO and centTRK vs. centV0
  if(event->GetRunNumber()>=167693 && event->GetRunNumber()<=170593){
    if(fRemoveTrackletOutliers){
      Double_t v0cent=GetCentrality((AliAODEvent*)event,kCentV0M);
      Double_t ntracklets=((AliAODEvent*)event)->GetTracklets()->GetNumberOfTracklets();
      Double_t cutval=60.-0.08*ntracklets+1./50000.*ntracklets*ntracklets;
      if(ntracklets<1000. && v0cent<cutval){
	if(accept) fWhyRejection=2;
	fEvRejectionBits+=1<<kOutsideCentrality;
	 accept=kFALSE;
      }
    }
    if(fMaxDiffTRKV0Centr>0.){
      Double_t v0cent=GetCentrality((AliAODEvent*)event,kCentV0M);
      Double_t trkcent=GetCentrality((AliAODEvent*)event,kCentTRK);
      if(TMath::Abs(trkcent-v0cent)>fMaxDiffTRKV0Centr){
	if(accept) fWhyRejection=1;
	fEvRejectionBits+=1<<kBadTrackV0Correl;
	accept=kFALSE;
      }
    }
  }

  // cuts on correlations between centrality estimators in Pb-Pb 2015 and Pb-Pb 2018
  if(fApplyCentralityCorrCuts && doAliEvCuts){
    if(!fAliEventCuts->PassedCut(AliEventCuts::kMultiplicity)){
      if(accept) fWhyRejection=1; // for fWhyRejection they are classified as pileup
      fEvRejectionBits+=1<<kBadTrackV0Correl;
      accept=kFALSE;
    }
  }

  // vertex requirements

  const AliVVertex *vertex = event->GetPrimaryVertex();

  if(!vertex){
    accept=kFALSE;
    fEvRejectionBits+=1<<kNoVertex;
  }else{
    TString title=vertex->GetTitle();
    if(title.Contains("Z") && fMinVtxType>1){
      accept=kFALSE;
      fEvRejectionBits+=1<<kNoVertex;
    }
    else if(title.Contains("3D") && fMinVtxType>2){
      accept=kFALSE;
      fEvRejectionBits+=1<<kNoVertex;
    }
    if(vertex->GetNContributors()<fMinVtxContr){
      accept=kFALSE;
      fEvRejectionBits+=1<<kTooFewVtxContrib;
    }
  }

  const AliVVertex *vSPD = ((AliAODEvent*)event)->GetPrimaryVertexSPD();
  if(fCutOnzVertexSPD>0){
    if(!vSPD || (vSPD && vSPD->GetNContributors()<fMinVtxContr)){
      accept=kFALSE;
      fEvRejectionBits+=1<<kBadSPDVertex;
    }else{
      if(fCutOnzVertexSPD==1 && TMath::Abs(vSPD->GetZ())>12.) {
        // protection for events with bad reconstructed track vertex (introduced for 2011 Pb-Pb)
        fEvRejectionBits+=1<<kZVtxSPDOutFid;
        if(accept) fWhyRejection=6;
        accept=kFALSE;
      }
      if(fCutOnzVertexSPD>=2 && vertex){
        Double_t dz = vSPD->GetZ()-vertex->GetZ();
        // cut on absolute distance between track and SPD vertex (introduced for 2011 Pb-Pb)
        if(TMath::Abs(dz)>0.5) {
          fEvRejectionBits+=1<<kBadTrackVertex;
          if(accept) fWhyRejection=0;
          accept=kFALSE;
        }
        if(accept && fCutOnzVertexSPD==3){
          // cut on nsigma distance between track and SPD vertex (for 2015 Pb-Pb)
          double covTrc[6],covSPD[6];
          vertex->GetCovarianceMatrix(covTrc);
          vSPD->GetCovarianceMatrix(covSPD);
          double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
          double errTrc = TMath::Sqrt(covTrc[5]);
          double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
          if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20){
            // reject, bad reconstructed track vertex
            fEvRejectionBits+=1<<kBadTrackVertex;
            if(accept) fWhyRejection=0;
            accept=kFALSE;
          }
        }
      }
    }
  }

  Double_t zvert = -999;
  if(!fApplyZcutOnSPDvtx && vertex) zvert=vertex->GetZ();
  else if(fApplyZcutOnSPDvtx) {
    if(!vSPD || (vSPD && vSPD->GetNContributors()<1)){
      accept=kFALSE;
      fEvRejectionBits+=1<<kBadSPDVertex;
    }
    else zvert = vSPD->GetZ();
  }

  if(TMath::Abs(zvert)>fMaxVtxZ) {
    fEvRejectionBits+=1<<kZVtxOutFid;
    if(accept) fWhyRejection=6;
    accept=kFALSE;
  }

  // pile-up rejection
  if(fOptPileup==kRejectPileupEvent){
    Bool_t isPileup=kFALSE;
    if(fUseMultDepPileupCut){
      isPileup=event->IsPileupFromSPDInMultBins();
    }else{
      Int_t cutc=(Int_t)fMinContrPileup;
      Double_t cutz=(Double_t)fMinDzPileup;
      isPileup=event->IsPileupFromSPD(cutc,cutz,3.,2.,10.);
    }
    if(isPileup){
      if(accept) fWhyRejection=1;
      fEvRejectionBits+=1<<kPileup;
      accept=kFALSE;
    }
  }
  else if(fOptPileup==kRejectMVPileupEvent){
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(fMinContrPileupMV);  // min. multiplicity of the pile-up vertex to consider
    utils.SetMaxPlpChi2MV(fMaxVtxChi2PileupMV); // max chi2 per contributor of the pile-up vertex to consider.
    utils.SetMinWDistMV(fMinWDzPileupMV);       // minimum weighted distance in Z between 2 vertices (i.e. (zv1-zv2)/sqrt(sigZv1^2+sigZv2^2) )
    utils.SetCheckPlpFromDifferentBCMV(fRejectPlpFromDiffBCMV); // vertex with |BCID|>2 will trigger pile-up (based on TOF)
    Bool_t isPUMV = utils.IsPileUpMV(event);
    if(isPUMV) {
      if(accept) fWhyRejection=1;
      fEvRejectionBits+=1<<kPileup;
      accept=kFALSE;
    }
  }

  // cut on correlations for out of bunch pileup in PbPb run2
  if(fApplyPbPbOutOfBunchPileupCuts==1 && doAliEvCuts){
    if(!fAliEventCuts->PassedCut(AliEventCuts::kCorrelations)){
      if(accept) fWhyRejection=1; // for fWhyRejection they are classified as pileup
      fEvRejectionBits+=1<<kBadCentrEstimCorrel;
      accept=kFALSE;
    }
  }else if(fApplyPbPbOutOfBunchPileupCuts==2 && event->GetRunNumber() >= 295369 && event->GetRunNumber() <= 297624){
    // Ionut cut on V0multiplicity vs. n TPC clusters (Pb-Pb 2018)
    AliAODVZERO* v0data=(AliAODVZERO*)((AliAODEvent*)event)->GetVZEROData();
    Float_t mTotV0=v0data->GetMTotV0A()+v0data->GetMTotV0C();
    Int_t nTPCcls=((AliAODEvent*)event)->GetNumberOfTPCClusters();
    Float_t mV0TPCclsCut=-2000.+(0.013*nTPCcls)+(1.25e-9*nTPCcls*nTPCcls);
    if(mTotV0<mV0TPCclsCut){
      if(accept) fWhyRejection=1;
      fEvRejectionBits+=1<<kBadTrackV0Correl;
      accept=kFALSE;
    }
  }

  // cut on ITS-TPC multiplicity correlation for OOB TPC pileup
  if(fApplyPbPbOutOfBunchPileupCutsITSTPC) {
    if(!fAliEventCuts->PassedCut(AliEventCuts::kTPCPileUp)){
      if(accept) fWhyRejection=1;
      fEvRejectionBits+=1<<kBadTPCITSCorrel;
      accept=kFALSE;      
    }
  }

  // Correcting PP2012 flag to remove tracks crossing SPD misaligned staves for periods 12def
  if(fApplySPDMisalignedPP2012 && !(event->GetRunNumber()>=195681 && event->GetRunNumber()<=197388)) fApplySPDMisalignedPP2012=false;

  // keep pileup events instead of rejecting them (only if good for all other requirementss)
  if(fKeepOnlyPbPbOutOfBunchPileupCutsITSTPC) {
    if(!(fEvRejectionBits&BIT(kBadTPCITSCorrel)))
    {
      fWhyRejection=1;
      return kFALSE; // reject if it was not rejected by pileup
    }
    if(!accept && fEvRejectionBits == BIT(kBadTPCITSCorrel))
    {
      fWhyRejection=0;
      return kTRUE; // accept if it was rejected only by OOB pileup
    }
  }

  return accept;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::IsEventSelectedWithAliEventCuts(AliVEvent *event) {
  //
  // Event selection with AliEventCuts
  //

  fWhyRejection=0;
  fEvRejectionBits=0;
  Bool_t accept=kTRUE;

  // check if it's MC
  Bool_t isMC=kFALSE;
  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(mcArray) {isMC=kTRUE;fUseAOD049=kFALSE;}


  SetupPID(event);

  // TEMPORARY FIX FOR GetEvent
  Int_t nTracks=((AliAODEvent*)event)->GetNumberOfTracks();
  for(Int_t itr=0; itr<nTracks; itr++){
    AliAODTrack* tr=(AliAODTrack*)((AliAODEvent*)event)->GetTrack(itr);
    tr->SetAODEvent((AliAODEvent*)event);
  }

  Int_t runNumb=event->GetRunNumber();
  if(runNumb >= 244917 && runNumb <= 246994) fAliEventCuts->SetupRun2PbPb();
  else if(runNumb >= 295369 && runNumb <= 297624) fAliEventCuts->SetupPbPb2018();
  else fAliEventCuts->SetManualMode(kFALSE);

  if(fUseTimeRangeCutForPbPb2018) fAliEventCuts->UseTimeRangeCut();
    
  // setup cuts
  TString selTrigClassClass="";
  if(!isMC && (event->GetRunNumber()<136851 || event->GetRunNumber()>139517)) {
    // don't do for MC and for PbPb 2010 data
    if(fTriggerClass[0].Length()>0){
      selTrigClassClass=fTriggerClass[0].Data();
      if(fTriggerClass[1].Length()>0) selTrigClassClass.Append(Form(",%s",fTriggerClass[1].Data()));
    }
  }
  fAliEventCuts->OverrideAutomaticTriggerSelection(fTriggerMask);
  if(fUseOnlyOneTrigger) fAliEventCuts->fRequireExactTriggerMask=kTRUE;
  fAliEventCuts->OverridePileUpCuts(fMinContrPileup,fMinDzPileup,3.,2.,5.);
  fAliEventCuts->fTrackletBGcut=kFALSE;
  if(fOptPileup==kRejectPileupEvent){
    fAliEventCuts->fPileUpCutMV=kFALSE;
    fAliEventCuts->fUseSPDpileUpCut=kTRUE;
    if(fUseMultDepPileupCut) fAliEventCuts->fUseMultiplicityDependentPileUpCuts=kTRUE;
  }
  else if(fOptPileup==kRejectMVPileupEvent){
    fAliEventCuts->fPileUpCutMV=kTRUE;
    fAliEventCuts->fUseSPDpileUpCut=kFALSE;
    fAliEventCuts->fUtils.SetMinPlpContribMV(fMinContrPileupMV);
    fAliEventCuts->fUtils.SetMaxPlpChi2MV(fMaxVtxChi2PileupMV);
    fAliEventCuts->fUtils.SetMinWDistMV(fMinWDzPileupMV);
    fAliEventCuts->fUtils.SetCheckPlpFromDifferentBCMV(fRejectPlpFromDiffBCMV);
  }
  fAliEventCuts->SetMaxVertexZposition(fMaxVtxZ);

  if(fApplyPbPbOutOfBunchPileupCutsITSTPC)
    fAliEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(true, fApplyPbPbOutOfBunchPileupCutsITSTPC);

  fAliEventCuts->AcceptEvent(event);

  // trigger class
  if(selTrigClassClass.Length()>0 && !fAliEventCuts->PassedCut(AliEventCuts::kTriggerClasses)){
    fWhyRejection=5;
    fEvRejectionBits+=1<<kNotSelTrigger;
    accept=kFALSE;
  }

  // physics selection requirements
  if(fUsePhysicsSelection){
    if(!fAliEventCuts->PassedCut(AliEventCuts::kTrigger)){
      if(accept) fWhyRejection=7;
      fEvRejectionBits+=1<<kPhysicsSelection;
      accept=kFALSE;
    }else{
      if(fUseV0ANDSelectionOffline){
	AliAODVZERO* v0data=(AliAODVZERO*)((AliAODEvent*)event)->GetVZEROData();
	Int_t tv0a=v0data->GetV0ADecision();
	Int_t tv0c=v0data->GetV0CDecision();
	if(!(tv0a==1 && tv0c==1)){
	  if(accept) fWhyRejection=7;
	  fEvRejectionBits+=1<<kPhysicsSelection;
	  accept=kFALSE;
	}
      }
    }
  }
  if(fUseTimeRangeCutForPbPb2018){
    if(!fAliEventCuts->PassedCut(AliEventCuts::kTimeRangeCut)){
      // use same fWhyRejection as for physics selection, to have proper counting of events for norm
      if(accept) fWhyRejection=7;
      fEvRejectionBits+=1<<kBadTimeRange;
      accept=kFALSE;
    }
  }
  
  // centrality selection
  if (fUseCentrality!=kCentOff) {
    Int_t rejection=IsEventSelectedInCentrality(event);
    Bool_t okCent=kFALSE;
    if(rejection==0) okCent=kTRUE;
    if(isMC && rejection==4 && !fUseCentrFlatteningInMC) okCent=kTRUE;
    if(!okCent){
      if(accept) fWhyRejection=rejection;
      if(fWhyRejection==4)fEvRejectionBits+=1<<kCentralityFlattening;
      else fEvRejectionBits+=1<<kOutsideCentrality;
      accept=kFALSE;
    }
  }

  // PbPb2011 outliers in tracklets vs. VZERO and centTRK vs. centV0
  if(event->GetRunNumber()>=167693 && event->GetRunNumber()<=170593){
    if(fRemoveTrackletOutliers){
      Double_t v0cent=GetCentrality((AliAODEvent*)event,kCentV0M);
      Double_t ntracklets=((AliAODEvent*)event)->GetTracklets()->GetNumberOfTracklets();
      Double_t cutval=60.-0.08*ntracklets+1./50000.*ntracklets*ntracklets;
      if(ntracklets<1000. && v0cent<cutval){
	if(accept) fWhyRejection=2;
	fEvRejectionBits+=1<<kOutsideCentrality;
	 accept=kFALSE;
      }
    }
    if(fMaxDiffTRKV0Centr>0.){
      Double_t v0cent=GetCentrality((AliAODEvent*)event,kCentV0M);
      Double_t trkcent=GetCentrality((AliAODEvent*)event,kCentTRK);
      if(TMath::Abs(trkcent-v0cent)>fMaxDiffTRKV0Centr){
	if(accept) fWhyRejection=1;
	fEvRejectionBits+=1<<kBadTrackV0Correl;
	accept=kFALSE;
      }
    }
  }

  // cuts on correlations between centrality estimators in Pb-Pb 2015 and Pb-Pb 2018
  if(fApplyCentralityCorrCuts){
    if(!fAliEventCuts->PassedCut(AliEventCuts::kMultiplicity)){
      if(accept) fWhyRejection=1; // for fWhyRejection they are classified as pileup
      fEvRejectionBits+=1<<kBadTrackV0Correl;
      accept=kFALSE;
    }
  }

  // vertex requirements

  if(fMinVtxType>2 && !fAliEventCuts->PassedCut(AliEventCuts::kVertexTracks)){
    accept=kFALSE;
    fEvRejectionBits+=1<<kNoVertex;
  }
  if(fMinVtxType<=2 && !fAliEventCuts->PassedCut(AliEventCuts::kVertexSPD)){
    accept=kFALSE;
    fEvRejectionBits+=1<<kNoVertex;
  }

  if(fCutOnzVertexSPD>0){
    if(!fAliEventCuts->PassedCut(AliEventCuts::kVertexSPD)){
      accept=kFALSE;
      fEvRejectionBits+=1<<kBadSPDVertex;
    }else{
      if(fCutOnzVertexSPD==1 && !fAliEventCuts->PassedCut(AliEventCuts::kVertexPositionSPD)){
        // protection for events with bad reconstructed track vertex (introduced for 2011 Pb-Pb)
        fEvRejectionBits+=1<<kZVtxSPDOutFid;
        if(accept) fWhyRejection=6;
        accept=kFALSE;
      }
      if(fCutOnzVertexSPD>=2 && !fAliEventCuts->PassedCut(AliEventCuts::kVertexQuality)){
	fEvRejectionBits+=1<<kBadTrackVertex;
	if(accept) fWhyRejection=0;
	accept=kFALSE;
      }
    }
  }

  if(!fAliEventCuts->PassedCut(AliEventCuts::kVertexPosition)){
    fEvRejectionBits+=1<<kZVtxOutFid;
    if(accept) fWhyRejection=6;
    accept=kFALSE;
  }

  // pile-up rejection
  if(fOptPileup==kRejectPileupEvent || fOptPileup==kRejectMVPileupEvent){
    if(!fAliEventCuts->PassedCut(AliEventCuts::kPileUp)){
      if(accept) fWhyRejection=1;
      fEvRejectionBits+=1<<kPileup;
      accept=kFALSE;
    }
  }
  // cut on correlations for out of bunch pileup in PbPb run2
  if(fApplyPbPbOutOfBunchPileupCuts==1){
    if(!fAliEventCuts->PassedCut(AliEventCuts::kCorrelations)){
      if(accept) fWhyRejection=1; // for fWhyRejection they are classified as pileup
      fEvRejectionBits+=1<<kBadCentrEstimCorrel;
      accept=kFALSE;
    }
  }else if(fApplyPbPbOutOfBunchPileupCuts==2 && event->GetRunNumber() >= 295369 && event->GetRunNumber() <= 297624){
    // Ionut cut on V0multiplicity vs. n TPC clusters (Pb-Pb 2018)
    AliAODVZERO* v0data=(AliAODVZERO*)((AliAODEvent*)event)->GetVZEROData();
    Float_t mTotV0=v0data->GetMTotV0A()+v0data->GetMTotV0C();
    Int_t nTPCcls=((AliAODEvent*)event)->GetNumberOfTPCClusters();
    Float_t mV0TPCclsCut=-2000.+(0.013*nTPCcls)+(1.25e-9*nTPCcls*nTPCcls);
    if(mTotV0<mV0TPCclsCut){
      if(accept) fWhyRejection=1;
      fEvRejectionBits+=1<<kBadTrackV0Correl;
      accept=kFALSE;
    }
  }

  // cut on ITS-TPC multiplicity correlation for OOB TPC pileup
  if(fApplyPbPbOutOfBunchPileupCutsITSTPC) {
    if(!fAliEventCuts->PassedCut(AliEventCuts::kTPCPileUp)){
      if(accept) fWhyRejection=1;
      fEvRejectionBits+=1<<kBadTPCITSCorrel;
      accept=kFALSE;      
    }
  }

  // Correcting PP2012 flag to remove tracks crossing SPD misaligned staves for periods 12def
  if(fApplySPDMisalignedPP2012 && !(event->GetRunNumber()>=195681 && event->GetRunNumber()<=197388)) fApplySPDMisalignedPP2012=false;

  return accept;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::AreDaughtersSelected(AliAODRecoDecayHF *d, const AliAODEvent* aod) const{
  //
  // Daughter tracks selection
  //
  if(!fTrackCuts) return kTRUE;

  Int_t ndaughters = d->GetNDaughters();
  AliAODVertex *vAOD = d->GetPrimaryVtx();
  Double_t pos[3],cov[6];
  vAOD->GetXYZ(pos);
  vAOD->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  Bool_t retval=kTRUE;

  for(Int_t idg=0; idg<ndaughters; idg++) {
    AliAODTrack *dgTrack = (AliAODTrack*)d->GetDaughter(idg);
    if(!dgTrack) {retval = kFALSE; continue;}
    //printf("charge %d\n",dgTrack->Charge());
    if(dgTrack->Charge()==0) continue; // it's not a track, but a V0

    if(fIsCandTrackSPDFirst && d->Pt()<fMaxPtCandTrackSPDFirst)
      { if(!dgTrack->HasPointOnITSLayer(0)) { retval = kFALSE; continue; } }

    if(!IsDaughterSelected(dgTrack,&vESD,fTrackCuts,aod)) retval = kFALSE;
  }

  return retval;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCuts::CheckMatchingAODdeltaAODevents(){
  //
  // Check if AOD and deltaAOD files are composed of the same events:
  // mismatches were observed in the merged AODs of LHC15o
  //
  // When AOD+deltaAOD are produced from ESD, mismatches can be found looking at:
  //   - the AOD trees in AliAOD.root and AliAOD.VertexingHF.root have different number of entries
  //   - the titles of the TProcessID objects do not match
  // When deltaAOD are produced from AOD, mismatches can be found looking at:
  //   - the AOD trees in AliAOD.root and AliAOD.VertexingHF.root have different number of entries
  //
  // Return values:
  //   -1: AOD and deltaAOD trees have different number of entries
  //    0: AOD and deltaAOD trees have same number of entries  +  the titles of the TProcessID objects do not match
  //    1: AOD and deltaAOD trees have same number of entries  +  the titles of the TProcessID objects match
  Bool_t okTProcessNames = kTRUE;
  AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  TTree *treeAOD      = aodHandler->GetTree();
  TTree *treeDeltaAOD = treeAOD->GetFriend("aodTree");
  TFile *mfile = treeAOD->GetCurrentFile();
  TFile *dfile = treeDeltaAOD->GetCurrentFile();
  treeAOD = (TTree*)mfile->Get("aodTree");
  treeDeltaAOD = (TTree*)dfile->Get("aodTree");
  if(!treeDeltaAOD || !treeAOD) return -1;
  if(treeAOD->GetEntries()!=treeDeltaAOD->GetEntries()){
      printf("AliRDHFCuts::CheckMatchingAODdeltaAODevents: Difference in number of entries in main and friend tree, skipping event\n");
      return -1;
  }
  TList* lm=mfile->GetListOfKeys();
  TList* ld=dfile->GetListOfKeys();
  Int_t nentm=lm->GetEntries();
  for(Int_t jm=0; jm<nentm; jm++){
    TKey* o=(TKey*)lm->At(jm);
    TString clnam=o->GetClassName();
    if(clnam=="TProcessID"){
      TString pname=o->GetName();
      TString ptit=o->GetTitle();
      if(pname.Contains("ProcessID")){
        TObject* od=(TObject*)ld->FindObject(pname.Data());
        if(od){
          TString ptit2=od->GetTitle();
          if(ptit2!=ptit){
            printf("AliRDHFCuts::CheckMatchingAODdeltaAODevents: mismatch in %s: AOD: %s  -- deltaAOD: %s\n",pname.Data(),ptit.Data(),ptit2.Data());
            okTProcessNames = kFALSE;
          }
        }
      }
    }
  }

  if (okTProcessNames) return 1;
  else return 0;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::CheckPtDepCrossedRows(TString rows,Bool_t print) const {
  //
  // Check the correctness of the string syntax
  //
  Bool_t retval=kTRUE;

  if(!rows.Contains("pt")) {
    if(print) AliError("string must contain \"pt\"");
    retval= kFALSE;
  }
  return retval;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetMinCrossedRowsTPCPtDep(const char *rows){
  //
  //Create the TFormula from TString for TPC crossed rows pT dependent cut
  //


  // setting data member that describes the TPC crossed rows pT dependent cut
  fCutMinCrossedRowsTPCPtDep = rows;

  // creating TFormula from TString
   if(f1CutMinNCrossedRowsTPCPtDep){
     delete f1CutMinNCrossedRowsTPCPtDep;
     // resetting TFormula
     f1CutMinNCrossedRowsTPCPtDep = 0;
   }
   if(!CheckPtDepCrossedRows(rows,kTRUE))return;

   TString tmp(rows);
   tmp.ReplaceAll("pt","x");
   f1CutMinNCrossedRowsTPCPtDep = new TFormula("f1CutMinNCrossedRowsTPCPtDep",tmp.Data());


}
//---------------------------------------------------------------------------
Bool_t AliRDHFCuts::IsDaughterSelected(AliAODTrack *track,const AliESDVertex *primary,AliESDtrackCuts *cuts, const AliAODEvent* aod) const{
  //
  // Convert to ESDtrack, relate to vertex and check cuts
  //
  if(!cuts) return kTRUE;

  if(cuts->GetFlagCutTOFdistance()) cuts->SetFlagCutTOFdistance(kFALSE);


  // convert to ESD track here
  AliESDtrack esdTrack(track);
  // set the TPC cluster info
  esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(track->GetTPCNclsF());
  // needed to calculate the impact parameters
  esdTrack.RelateToVertex(primary,0.,3.);

  //applying ESDtrackCut
  if(!cuts->IsSelected(&esdTrack)) return kFALSE;

  //appliyng kink rejection
  if(fKinkReject){
   AliAODVertex *maybeKink=track->GetProdVertex();
   if(maybeKink->GetType()==AliAODVertex::kKink) return kFALSE;
  }

  //appliyng TPC crossed rows pT dependent cut
  if(f1CutMinNCrossedRowsTPCPtDep && fUseTPCtrackCutsOnThisDaughter){
    Float_t nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
    if(nCrossedRowsTPC<f1CutMinNCrossedRowsTPCPtDep->Eval(esdTrack.Pt())) return kFALSE;
  }

  //appliyng NTPCcls/NTPCcrossedRows cut
  if(fCutRatioClsOverCrossRowsTPC && fUseTPCtrackCutsOnThisDaughter){
    Float_t nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
    Float_t nClustersTPC = esdTrack.GetTPCNcls();
    if(nCrossedRowsTPC!=0){
      Float_t ratio = nClustersTPC/nCrossedRowsTPC;
      if(ratio<fCutRatioClsOverCrossRowsTPC) return kFALSE;
    }
    else return kFALSE;
  }

  //appliyng TPCsignalN/NTPCcrossedRows cut
  if(fCutRatioSignalNOverCrossRowsTPC && fUseTPCtrackCutsOnThisDaughter){
    Float_t nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
    Float_t nTPCsignal = esdTrack.GetTPCsignalN();
    if(nCrossedRowsTPC!=0){
      Float_t ratio = nTPCsignal/nCrossedRowsTPC;
      if(ratio<fCutRatioSignalNOverCrossRowsTPC) return kFALSE;
    }
    else return kFALSE;
  }

  // cut on the number of TPC clusters for PID
  if(fCutTPCSignalN && fUseTPCtrackCutsOnThisDaughter){
    Float_t nTPCsignal = esdTrack.GetTPCsignalN();
    if(nTPCsignal<fCutTPCSignalN) return kFALSE;
  }
  
  // geometrical cut (note uses track at vertex instead of at TPC inner wall)
  if(fUseCutGeoNcrNcl && aod && fUseTPCtrackCutsOnThisDaughter){
    Float_t nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
    Float_t lengthInActiveZoneTPC=esdTrack.GetLengthInActiveZone(0,fDeadZoneWidth,220.,aod->GetMagneticField());
    Double_t cutGeoNcrNclLength=fCutGeoNcrNclLength-TMath::Power(TMath::Abs(esdTrack.GetSigned1Pt()),fCutGeoNcrNclGeom1Pt);
    Bool_t isOK=kTRUE;
    if (lengthInActiveZoneTPC<cutGeoNcrNclLength) isOK=kFALSE;
    if (nCrossedRowsTPC<fCutGeoNcrNclFractionNcr*cutGeoNcrNclLength) isOK=kFALSE;
    if (esdTrack.GetTPCncls()<fCutGeoNcrNclFractionNcl*cutGeoNcrNclLength) isOK=kFALSE;
    if(!isOK) return kFALSE;
  }


  if(fOptPileup==kRejectTracksFromPileupVertex){
    // to be implemented
    // we need either to have here the AOD Event,
    // or to have the pileup vertex object
  }

  if(fApplySPDDeadPbPb2011){
    Bool_t deadSPDLay1PbPb2011[20][4]={
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kFALSE,kFALSE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kFALSE,kFALSE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE}
    };
    Bool_t deadSPDLay2PbPb2011[40][4]={
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kFALSE,kFALSE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kTRUE,kTRUE,kFALSE,kFALSE},
      {kTRUE,kTRUE,kTRUE,kTRUE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE},
      {kFALSE,kFALSE,kFALSE,kFALSE}
    };
    Double_t xyz1[3],xyz2[3];
    esdTrack.GetXYZAt(3.9,0.,xyz1);
    esdTrack.GetXYZAt(7.6,0.,xyz2);
    Double_t phi1=TMath::ATan2(xyz1[1],xyz1[0]);
    if(phi1<0) phi1+=2*TMath::Pi();
    Int_t lad1=(Int_t)(phi1/(2.*TMath::Pi()/20.));
    Double_t phi2=TMath::ATan2(xyz2[1],xyz2[0]);
    if(phi2<0) phi2+=2*TMath::Pi();
    Int_t lad2=(Int_t)(phi2/(2.*TMath::Pi()/40.));
    Int_t mod1=TMath::Floor((xyz1[2]+14)/7.);
    Int_t mod2=TMath::Floor((xyz2[2]+14)/7.);
    Bool_t lay1ok=kFALSE;
    if(mod1>=0 && mod1<4 && lad1<20){
      lay1ok=deadSPDLay1PbPb2011[lad1][mod1];
    }
    Bool_t lay2ok=kFALSE;
    if(mod2>=0 && mod2<4 && lad2<40){
      lay2ok=deadSPDLay2PbPb2011[lad2][mod2];
    }
    if(!lay1ok && !lay2ok) return kFALSE;
  }

  if(fApplySPDMisalignedPP2012) {
    // Cut tracks crossing the SPD at 5.6<phi<2pi
    Double_t xyz1[3],xyz2[3];
    esdTrack.GetXYZAt(3.9,0.,xyz1);
    esdTrack.GetXYZAt(7.6,0.,xyz2);
    Double_t phi1=TMath::ATan2(xyz1[1],xyz1[0]);
    if(phi1<0) phi1+=2*TMath::Pi();
    Double_t phi2=TMath::ATan2(xyz2[1],xyz2[0]);
    if(phi2<0) phi2+=2*TMath::Pi();
    Bool_t lay1ok=kTRUE;
    if(phi1>5.6 && phi1<2.*TMath::Pi()) lay1ok=kFALSE;
    Bool_t lay2ok=kTRUE;
    if(phi2>5.6 && phi2<2.*TMath::Pi()) lay2ok=kFALSE;
    if(!lay1ok || !lay2ok) return kFALSE;
  }

  if(fApplySPDUniformAccPbPbRun2){
    // Cut tracks crossing the regions at
    // -->  0.94<phi<1.34 and 0<z<14cm (corresponding to the region of SPD HS 1A0)
    // -->  3.49<phi<3.9 and -14<z<0 cm (corresponding to the region of SPD HS 5C0)
    // Reason: HS 1A0 and 5C0 were excluded from 2018 Pb-Pb data taking, while were present in 2015 Pb-Pb data taking
    // Goal: this patch should allow to use the same SPD regions between the two periods
    Double_t xyz1[3];
    esdTrack.GetXYZAt(3.9,0.,xyz1);
    Double_t phi1=TMath::ATan2(xyz1[1],xyz1[0]);
    if((phi1<1.34 && phi1>0.94) && (xyz1[2]>0 && xyz1[2]<14)) return kFALSE; // exclude region of 1A0
    if((phi1<3.9 && phi1>3.49) && (xyz1[2]>-14 && xyz1[2]<0)) return kFALSE; //exclude region of 5C0
  }

  return kTRUE;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetPtBins(Int_t nPtBinLimits,Float_t *ptBinLimits) {
  // Set the pt bins

  if(fPtBinLimits) {
    delete [] fPtBinLimits;
    fPtBinLimits = NULL;
    printf("Changing the pt bins\n");
  }

  if(nPtBinLimits != fnPtBins+1){
    cout<<"Warning: ptBinLimits dimention "<<nPtBinLimits<<" != nPtBins+1 ("<<fnPtBins+1<<")\nSetting nPtBins to "<<nPtBinLimits-1<<endl;
    SetNPtBins(nPtBinLimits-1);
  }

  fnPtBinLimits = nPtBinLimits;
  SetGlobalIndex();
  //cout<<"Changing also Global Index -> "<<fGlobalIndex<<endl;
  fPtBinLimits = new Float_t[fnPtBinLimits];
  for(Int_t ib=0; ib<nPtBinLimits; ib++) fPtBinLimits[ib]=ptBinLimits[ib];

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetVarNames(Int_t nVars,TString *varNames,Bool_t *isUpperCut){
  // Set the variable names

  if(fVarNames) {
    delete [] fVarNames;
    fVarNames = NULL;
    //printf("Changing the variable names\n");
  }
  if(nVars!=fnVars){
    printf("Wrong number of variables: it has to be %d\n",fnVars);
    return;
  }
  //fnVars=nVars;
  fVarNames = new TString[nVars];
  fIsUpperCut = new Bool_t[nVars];
  for(Int_t iv=0; iv<nVars; iv++) {
    fVarNames[iv] = varNames[iv];
    fIsUpperCut[iv] = isUpperCut[iv];
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetVarsForOpt(Int_t nVars,Bool_t *forOpt) {
  // Set the variables to be used for cuts optimization

  if(fVarsForOpt) {
    delete [] fVarsForOpt;
    fVarsForOpt = NULL;
    //printf("Changing the variables for cut optimization\n");
  }

  if(nVars==0){//!=fnVars) {
    printf("%d not accepted as number of variables: it has to be %d\n",nVars,fnVars);
    return;
  }

  fnVarsForOpt = 0;
  fVarsForOpt = new Bool_t[fnVars];
  for(Int_t iv=0; iv<fnVars; iv++) {
    fVarsForOpt[iv]=forOpt[iv];
    if(fVarsForOpt[iv]) fnVarsForOpt++;
  }

  return;
}

//---------------------------------------------------------------------------
void AliRDHFCuts::SetUseCentrality(Int_t flag) {
  //
  // set centrality estimator
  //
  fUseCentrality=flag;
  if(fUseCentrality<kCentOff||fUseCentrality>=kCentInvalid) AliWarning("Centrality estimator not valid");

  return;
}


//---------------------------------------------------------------------------
void AliRDHFCuts::SetCuts(Int_t nVars,Int_t nPtBins,Float_t **cutsRD) {
  //
  // store the cuts
  //
  if(nVars!=fnVars) {
    printf("Wrong number of variables: it has to be %d\n",fnVars);
    AliFatal("exiting");
  }
  if(nPtBins!=fnPtBins) {
    printf("Wrong number of pt bins: it has to be %d\n",fnPtBins);
    AliFatal("exiting");
  }

  if(!fCutsRD)  fCutsRD = new Float_t[fGlobalIndex];


  for(Int_t iv=0; iv<fnVars; iv++) {

    for(Int_t ib=0; ib<fnPtBins; ib++) {

      //check
      if(GetGlobalIndex(iv,ib)>=fGlobalIndex) {
	cout<<"Overflow, exit..."<<endl;
	return;
      }

      fCutsRD[GetGlobalIndex(iv,ib)] = cutsRD[iv][ib];

    }
  }
  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::SetCuts(Int_t glIndex,Float_t* cutsRDGlob){
  //
  // store the cuts
  //
  if(glIndex != fGlobalIndex){
    cout<<"Wrong array size: it has to be "<<fGlobalIndex<<endl;
    AliFatal("exiting");
  }
  if(!fCutsRD)  fCutsRD = new Float_t[fGlobalIndex];

  for(Int_t iGl=0;iGl<fGlobalIndex;iGl++){
    fCutsRD[iGl] = cutsRDGlob[iGl];
  }
  return;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::PrintAll() const {
  //
  // print all cuts values
  //

  printf("---- Event Selecion Cuts ----\n");
  printf("Minimum vtx type %d\n",fMinVtxType);
  printf("Minimum vtx contr %d\n",fMinVtxContr);
  printf("Max vtx red chi2 %f\n",fMaxVtxRedChi2);
  printf("Min SPD mult %d\n",fMinSPDMultiplicity);
  printf("Remove daughters from vtx %d\n",(Int_t)fRemoveDaughtersFromPrimary);
  printf("Physics selection: %s\n",fUsePhysicsSelection ? "Yes" : "No");
  printf("Pileup rejection: %s\n",(fOptPileup > 0) ? "Yes" : "No");
  if(fOptPileup==1) printf(" -- Reject pileup event");
  if(fOptPileup==2) printf(" -- Reject tracks from pileup vtx");
  if(fUseCentrality>0) {
    TString estimator="";
    if(fUseCentrality==kCentV0M) estimator = "V0";
    if(fUseCentrality==kCentTRK) estimator = "Tracks";
    if(fUseCentrality==kCentTKL) estimator = "Tracklets";
    if(fUseCentrality==kCentCL1) estimator = "SPD clusters outer";
    if(fUseCentrality==kCentZNA) estimator = "ZNA";
    if(fUseCentrality==kCentZPA) estimator = "ZPA";
    if(fUseCentrality==kCentV0A) estimator = "V0A";
    if(fUseCentrality==kCentCL0) estimator = "SPD clusters inner";
    printf("Centrality class considered: %.1f-%.1f, estimated with %s\n",fMinCentrality,fMaxCentrality,estimator.Data());
  }

  printf("---- Single Track Cuts ----\n");
  TString itsSelString[8]={"kOff", "kNone", "kAny", "kFirst", "kOnlyFirst", "kSecond", "kOnlySecond", "kBoth"};
  if(fTrackCuts) {
    printf(" Require TPC refit                          = %d\n",fTrackCuts->GetRequireTPCRefit());
    printf(" Min. number of TPC Clusters                = %d\n",fTrackCuts->GetMinNClusterTPC());
    printf(" Min. number of TPC Crossed Rows            = %.0f\n",fTrackCuts->GetMinNCrossedRowsTPC());
    printf(" Min. ratio crossed rows /findable clusters = %f\n",fTrackCuts->GetMinRatioCrossedRowsOverFindableClustersTPC());
    printf(" Max. chi2/cluster TPC                      = %f\n",fTrackCuts->GetMaxChi2PerClusterTPC());
    printf(" Require ITS refit                          = %d\n",fTrackCuts->GetRequireITSRefit());
    printf(" Min. number of ITS Clusters                = %d\n",fTrackCuts->GetMinNClustersITS());
    printf(" Cluster requirement SPD                    = %s\n",itsSelString[fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)].Data());
    printf(" Cluster requirement SDD                    = %s\n",itsSelString[fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSDD)].Data());
    printf(" Cluster requirement SSD                    = %s\n",itsSelString[fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSSD)].Data());
    printf(" Max. chi2/cluster ITS                      = %f\n",fTrackCuts->GetMaxChi2PerClusterITS());
    printf(" Max. chi2 TPC constr-global (golden chi2)  = %f\n",fTrackCuts->GetMaxChi2TPCConstrainedGlobal());
    printf(" DCA to vertex (XY) Min - Max (cm)          = %f - %f\n",fTrackCuts->GetMinDCAToVertexXY(),fTrackCuts->GetMaxDCAToVertexXY());
    printf(" DCA to vertex (Z) Min - Max  (cm)          = %f - %f\n",fTrackCuts->GetMinDCAToVertexZ(),fTrackCuts->GetMaxDCAToVertexZ());
  }
  
  if(fCutRatioClsOverCrossRowsTPC) printf("N TPC Clusters > %f N TPC Crossed Rows\n", fCutRatioClsOverCrossRowsTPC);
  if(fCutRatioSignalNOverCrossRowsTPC) printf("N TPC Points for dE/dx > %f N TPC Crossed Rows\n", fCutRatioSignalNOverCrossRowsTPC);
  if(fCutTPCSignalN>0) printf("N TPC Clusters for PID for track sel > %d\n", fCutTPCSignalN);
  if(f1CutMinNCrossedRowsTPCPtDep) printf("N TPC Crossed Rows pT-dependent cut: %s\n", fCutMinCrossedRowsTPCPtDep.Data());
  if(fIsCandTrackSPDFirst) printf("Check for candidates with pt < %2.2f, that daughters fullfill kFirst criteria\n",fMaxPtCandTrackSPDFirst);

  printf("---- Candidate Cuts ----\n");
  if(fVarNames){
    cout<<"Array of variables"<<endl;
    for(Int_t iv=0;iv<fnVars;iv++){
      cout<<fVarNames[iv]<<"\t";
    }
    cout<<endl;
  }
  if(fVarsForOpt){
    cout<<"Array of optimization"<<endl;
    for(Int_t iv=0;iv<fnVars;iv++){
      cout<<fVarsForOpt[iv]<<"\t";
    }
    cout<<endl;
  }
  if(fIsUpperCut){
    cout<<"Array of upper/lower cut"<<endl;
   for(Int_t iv=0;iv<fnVars;iv++){
     cout<<fIsUpperCut[iv]<<"\t";
   }
   cout<<endl;
  }
  if(fPtBinLimits){
    cout<<"Array of ptbin limits"<<endl;
    for(Int_t ib=0;ib<fnPtBinLimits;ib++){
      cout<<fPtBinLimits[ib]<<"\t";
    }
    cout<<endl;
  }
  if(fCutsRD){
    cout<<"Matrix of cuts"<<endl;
   for(Int_t iv=0;iv<fnVars;iv++){
     for(Int_t ib=0;ib<fnPtBins;ib++){
       cout<<"fCutsRD["<<iv<<"]["<<ib<<"] = "<<fCutsRD[GetGlobalIndex(iv,ib)]<<"\t";
     }
     cout<<endl;
   }
   cout<<endl;
  }
  printf("fUsePreselect=%d \n",fUsePreselect);
  printf("---- PID Cuts ----\n");
  printf("Use PID %d  OldPid=%d\n",(Int_t)fUsePID,fPidHF ? fPidHF->GetOldPid() : -1);
  if(fPidHF) fPidHF->PrintAll();
  Printf("EnableNSigmaTPCDataCorr = %d, %d", fEnableNsigmaTPCDataCorr, fSystemForNsigmaTPCDataCorr);

  return;
}

//--------------------------------------------------------------------------
void AliRDHFCuts::PrintTrigger() const{
  // print the trigger selection

  printf("Selected trigger classes: %s %s\n",fTriggerClass[0].Data(),fTriggerClass[1].Data());

  cout<<" Trigger selection pattern: ";
  if( fTriggerMask & AliVEvent::kAny ) cout<<" kAny ";
  if( fTriggerMask & AliVEvent::kAnyINT ) cout<<" kAnyINT ";
  if( fTriggerMask & AliVEvent::kINT7 ) cout<<" kINT7 ";
  if( fTriggerMask & AliVEvent::kMB ) cout<<" kMB ";
  if( fTriggerMask & AliVEvent::kCINT5 ) cout<<" kCINT5 ";
  if( fTriggerMask & AliVEvent::kCentral ) cout<<" kCentral ";
  if( fTriggerMask & AliVEvent::kSemiCentral ) cout<<" kSemiCentral ";
  if( fTriggerMask & AliVEvent::kEMCEGA ) cout<<" kEMCEGA ";
  if( fTriggerMask & AliVEvent::kHighMult ) cout<<" kHighMult ";
  if( fTriggerMask & AliVEvent::kFastOnly ) cout<<" kFastOnly ";
  cout << endl<< endl;

}

//---------------------------------------------------------------------------
void AliRDHFCuts::GetCuts(Float_t**& cutsRD) const{
  //
  // get the cuts
  //

  //cout<<"Give back a "<<fnVars<<"x"<<fnPtBins<<" matrix."<<endl;


  Int_t iv,ib;
  if(!cutsRD) {
    //cout<<"Initialization..."<<endl;
    cutsRD=new Float_t*[fnVars];
    for(iv=0; iv<fnVars; iv++) {
      cutsRD[iv] = new Float_t[fnPtBins];
    }
  }

  for(Int_t iGlobal=0; iGlobal<fGlobalIndex; iGlobal++) {
    GetVarPtIndex(iGlobal,iv,ib);
    cutsRD[iv][ib] = fCutsRD[iGlobal];
  }

  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCuts::GetGlobalIndex(Int_t iVar,Int_t iPtBin) const{
  //
  // give the global index from variable and pt bin
  //
  return iPtBin*fnVars+iVar;
}

//---------------------------------------------------------------------------
void AliRDHFCuts::GetVarPtIndex(Int_t iGlob, Int_t& iVar, Int_t& iPtBin) const {
  //
  //give the index of the variable and of the pt bin from the global index
  //
  iPtBin=(Int_t)iGlob/fnVars;
  iVar=iGlob%fnVars;

  return;
}

//-------------------------------------------------------------------
Float_t AliRDHFCuts::GetCutValue(Int_t iVar,Int_t iPtBin) const {
  //
  // Give the value of cut set for the variable iVar and the pt bin iPtBin
  //
  if(!fCutsRD){
    cout<<"Cuts not iniziaisez yet"<<endl;
    return 0;
  }
  return fCutsRD[GetGlobalIndex(iVar,iPtBin)];
}

//-------------------------------------------------------------------
Float_t AliRDHFCuts::GetCentrality(AliAODEvent* aodEvent,AliRDHFCuts::ECentrality estimator) {

  if(aodEvent->GetRunNumber()<244824)return GetCentralityOldFramework(aodEvent,estimator);
  Double_t cent=-999;

  if(estimator==kCentOff) return -999;

  AliMultSelection *multSelection = (AliMultSelection*)aodEvent->FindListObject(fMultSelectionObjectName);
  if(!multSelection){
      AliWarning("AliMultSelection could not be found in the aod event list of objects");
      return cent;
  }

  // Compare centrality with the centrality at AOD filtering and on-the-fly
  if( fMultSelectionObjectName.CompareTo("MultSelection")!=0 ){
      AliMultSelection *defaultmultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");
      if(!defaultmultSelection){
          AliWarning("AliMultSelection default method could not be found in the aod event list of objects");
          return cent;
      }
      Float_t defaultCent = defaultmultSelection->GetMultiplicityPercentile("V0M");
      Float_t newCent = multSelection->GetMultiplicityPercentile("V0M");
      if( defaultCent<20. && newCent>20.) fEvRejectionBits+=1<<kMismatchOldNewCentrality;
      else if (defaultCent>20. && newCent<20.) fEvRejectionBits+=1<<kMismatchOldNewCentrality;
  }

  if(estimator==kCentV0M){
    cent=multSelection->GetMultiplicityPercentile("V0M");
  }else if(estimator==kCentV0A){
    cent=multSelection->GetMultiplicityPercentile("V0A");
  }else if(estimator==kCentZNA){
    cent=multSelection->GetMultiplicityPercentile("ZNA");
  }else if(estimator==kCentCL1){
    cent=multSelection->GetMultiplicityPercentile("CL1");
  }else if(estimator==kCentCL0){
    cent=multSelection->GetMultiplicityPercentile("CL0");
  }else {
    AliWarning(Form("CENTRALITY ESTIMATE WITH ESTIMATOR %d NOT YET IMPLEMENTED FOR NEW FRAMEWORK",(Int_t)estimator));
    return cent;
  }
  Int_t qual = multSelection->GetEvSelCode();
  if(qual == 199 ) cent=-999;
  return cent;
}
//-------------------------------------------------------------------
Float_t AliRDHFCuts::GetCentralityOldFramework(AliAODEvent* aodEvent,AliRDHFCuts::ECentrality estimator) {
  //
  // Get centrality percentile
  //

  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)aodEvent)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(mcArray) {fUseAOD049=kFALSE;}

  AliAODHeader *header=dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
  if(!header) AliFatal("Not a standard AOD");
  AliCentrality *centrality=header->GetCentralityP();
  Float_t cent=-999.;
  Bool_t isSelRun=kFALSE;
  Int_t selRun[5]={138364, 138826, 138828, 138836, 138871};
  if(!centrality) return cent;
  else{
    if (estimator==kCentV0M){
      cent=(Float_t)(centrality->GetCentralityPercentile("V0M"));
      if(cent<0){
	Int_t quality = centrality->GetQuality();
	if(quality<=1){ // fQuality==1 means rejected by zVertex cut that we apply a part and we want to keep separate (Giacomo)
	  cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
	}else{
	  Int_t runnum=aodEvent->GetRunNumber();
	  for(Int_t ir=0;ir<5;ir++){
	    if(runnum==selRun[ir]){
	      isSelRun=kTRUE;
	      break;
	    }
	  }
	  if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
	}
      }

      //temporary fix for AOD049 outliers
      if(fUseAOD049&&cent>=0){
	Float_t v0=0;
	AliAODVZERO* aodV0 = aodEvent->GetVZEROData();
	v0+=aodV0->GetMTotV0A();
	v0+=aodV0->GetMTotV0C();
	if(cent==0&&v0<19500)return -1;//filtering issue
	Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
	Float_t val= 1.30552 +  0.147931 * v0;
	Float_t tklSigma[101]={176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86, 120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654, 92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334, 68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224, 51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255, 37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398, 26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235, 19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504, 12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544};
	if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)cent] )return -1;//outlier
      }
    }
    else {
       if (estimator==kCentTRK) {
	cent=(Float_t)(centrality->GetCentralityPercentile("TRK"));
	if(cent<0){
	  Int_t quality = centrality->GetQuality();
	  if(quality<=1){
	    cent=(Float_t)centrality->GetCentralityPercentileUnchecked("TRK");
	  }else{
	    Int_t runnum=aodEvent->GetRunNumber();
	    for(Int_t ir=0;ir<5;ir++){
	      if(runnum==selRun[ir]){
		isSelRun=kTRUE;
		break;
	      }
	    }
	    if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("TRK");
	  }
	}
       }
      else{
	if (estimator==kCentTKL){
	  cent=(Float_t)(centrality->GetCentralityPercentile("TKL"));
	  if(cent<0){
	    Int_t quality = centrality->GetQuality();
	    if(quality<=1){
	      cent=(Float_t)centrality->GetCentralityPercentileUnchecked("TKL");
	    }else{
	      Int_t runnum=aodEvent->GetRunNumber();
	      for(Int_t ir=0;ir<5;ir++){
		if(runnum==selRun[ir]){
		  isSelRun=kTRUE;
		  break;
	    }
	      }
	      if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("TKL");
	    }
	  }
	}
	else{
	  if (estimator==kCentCL1){
	    cent=(Float_t)(centrality->GetCentralityPercentile("CL1"));
	    if(cent<0){
	      Int_t quality = centrality->GetQuality();
	      if(quality<=1){
		cent=(Float_t)centrality->GetCentralityPercentileUnchecked("CL1");
	      }else{
		Int_t runnum=aodEvent->GetRunNumber();
		for(Int_t ir=0;ir<5;ir++){
		  if(runnum==selRun[ir]){
		    isSelRun=kTRUE;
		    break;
		  }
		}
		if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("CL1");
	      }
	    }
	  }
	else{
	  if (estimator==kCentZNA){
	    cent=(Float_t)(centrality->GetCentralityPercentile("ZNA"));
	    if(cent<0){
	      Int_t quality = centrality->GetQuality();
	      if(quality<=1){
		cent=(Float_t)centrality->GetCentralityPercentileUnchecked("ZNA");
	      }else{
		Int_t runnum=aodEvent->GetRunNumber();
		for(Int_t ir=0;ir<5;ir++){
		  if(runnum==selRun[ir]){
		    isSelRun=kTRUE;
		    break;
		  }
		}
		if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("ZNA");
	      }
	    }
	  }
	else{
	  if (estimator==kCentZPA){
	    cent=(Float_t)(centrality->GetCentralityPercentile("ZPA"));
	    if(cent<0){
	      Int_t quality = centrality->GetQuality();
	      if(quality<=1){
		cent=(Float_t)centrality->GetCentralityPercentileUnchecked("ZPA");
	      }else{
		Int_t runnum=aodEvent->GetRunNumber();
		for(Int_t ir=0;ir<5;ir++){
		  if(runnum==selRun[ir]){
		    isSelRun=kTRUE;
		    break;
		  }
		}
		if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("ZPA");
	      }
	    }
	  }
	else{
	  if (estimator==kCentV0A){
	    cent=(Float_t)(centrality->GetCentralityPercentile("V0A"));
	    if(cent<0){
	      Int_t quality = centrality->GetQuality();
	      if(quality<=1){
		cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0A");
	      }else{
		Int_t runnum=aodEvent->GetRunNumber();
		for(Int_t ir=0;ir<5;ir++){
		  if(runnum==selRun[ir]){
		    isSelRun=kTRUE;
		    break;
		  }
		}
		if((quality==8||quality==9)&&isSelRun)cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0A");
	      }
	    }
	  }
	  else {
	    AliWarning("Centrality estimator not valid");

	  }
	}
    }
    }
    }
    }
    }
  }
  return cent;
}
//-------------------------------------------------------------------
Bool_t AliRDHFCuts::CompareCuts(const AliRDHFCuts *obj) const {
  //
  // Compare two cuts objects
  //

  Bool_t areEqual=kTRUE;

  if(fMinVtxType!=obj->fMinVtxType) { printf("Minimum vtx type %d  %d\n",fMinVtxType,obj->fMinVtxType); areEqual=kFALSE;}

  if(fMinVtxContr!=obj->fMinVtxContr) { printf("Minimum vtx contr %d  %d\n",fMinVtxContr,obj->fMinVtxContr); areEqual=kFALSE;}

  if(TMath::Abs(fMaxVtxRedChi2-obj->fMaxVtxRedChi2)>1.e-10) {   printf("Max vtx red chi2 %f  %f\n",fMaxVtxRedChi2,obj->fMaxVtxRedChi2);areEqual=kFALSE;}

  if(fMinSPDMultiplicity!=obj->fMinSPDMultiplicity) {  printf("Min SPD mult %d\n  %d",fMinSPDMultiplicity,obj->fMinSPDMultiplicity);areEqual=kFALSE;}

  if(fUsePID!=obj->fUsePID) { printf("Use PID %d  %d\n",(Int_t)fUsePID,(Int_t)obj->fUsePID); areEqual=kFALSE;}

  if(fRemoveDaughtersFromPrimary!=obj->fRemoveDaughtersFromPrimary) {printf("Remove daughters from vtx %d  %d\n",(Int_t)fRemoveDaughtersFromPrimary,(Int_t)obj->fRemoveDaughtersFromPrimary); areEqual=kFALSE;}
  if(fTrackCuts){
    if(fTrackCuts->GetMinNClusterTPC()!=obj->fTrackCuts->GetMinNClusterTPC()) {printf("MinNClsTPC %d  %d\n",fTrackCuts->GetMinNClusterTPC(),obj->fTrackCuts->GetMinNClusterTPC()); areEqual=kFALSE;}

    if(fTrackCuts->GetMinNClustersITS()!=obj->fTrackCuts->GetMinNClustersITS()) {printf("MinNClsITS %d  %d\n",fTrackCuts->GetMinNClustersITS(),obj->fTrackCuts->GetMinNClustersITS()); areEqual=kFALSE;}

    if(TMath::Abs(fTrackCuts->GetMaxChi2PerClusterTPC()-obj->fTrackCuts->GetMaxChi2PerClusterTPC())>1.e-10) {printf("MaxChi2ClsTPC %f  %f\n",fTrackCuts->GetMaxChi2PerClusterTPC(),obj->fTrackCuts->GetMaxChi2PerClusterTPC()); areEqual=kFALSE;}

    if(fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)!=obj->fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)) {printf("ClusterReq SPD %d  %d\n",fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD),obj->fTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD)); areEqual=kFALSE;}
  }

  if(fUsePreselect!=obj->fUsePreselect){printf("fUsePreselect: %d %d\n",fUsePreselect,obj->fUsePreselect);areEqual=kFALSE;}

  if(fCutsRD) {
   for(Int_t iv=0;iv<fnVars;iv++) {
     for(Int_t ib=0;ib<fnPtBins;ib++) {
       if(TMath::Abs(fCutsRD[GetGlobalIndex(iv,ib)]-obj->fCutsRD[GetGlobalIndex(iv,ib)])>1.e-10) {
	 cout<<"fCutsRD["<<iv<<"]["<<ib<<"] = "<<fCutsRD[GetGlobalIndex(iv,ib)]<<"   "<<obj->fCutsRD[GetGlobalIndex(iv,ib)]<<"\n";
	 areEqual=kFALSE;
       }
     }
   }
  }

  return areEqual;
}
//---------------------------------------------------------------------------
void AliRDHFCuts::MakeTable() const {
  //
  // print cuts values in table format
  //

	TString ptString = "pT range";
	if(fVarNames && fPtBinLimits && fCutsRD){
		TString firstLine(Form("*       %-15s",ptString.Data()));
		for (Int_t ivar=0; ivar<fnVars; ivar++){
			firstLine+=Form("*    %-15s  ",fVarNames[ivar].Data());
			if (ivar == fnVars){
				firstLine+="*\n";
			}
		}
		Printf("%s",firstLine.Data());

		for (Int_t ipt=0; ipt<fnPtBins; ipt++){
			TString line;
			if (ipt==fnPtBins-1){
				line=Form("*  %5.1f < pt < inf    ",fPtBinLimits[ipt]);
			}
			else{
				line=Form("*  %5.1f < pt < %4.1f   ",fPtBinLimits[ipt],fPtBinLimits[ipt+1]);
			}
			for (Int_t ivar=0; ivar<fnVars; ivar++){
				line+=Form("*     %-15f ",fCutsRD[GetGlobalIndex(ivar,ipt)]);
			}
			Printf("%s",line.Data());
		}

	}


  return;
}
//--------------------------------------------------------------------------
Bool_t AliRDHFCuts::RecalcOwnPrimaryVtx(AliAODRecoDecayHF *d,
					AliAODEvent *aod) const
{
  //
  // Recalculate primary vertex without daughters
  //

  if(!aod) {
    AliError("Can not remove daughters from vertex without AOD event");
    return 0;
  }

  AliAODVertex *recvtx=d->RemoveDaughtersFromPrimaryVtx(aod);
  if(!recvtx){
    AliDebug(2,"Removal of daughter tracks failed");
    return kFALSE;
  }


  //set recalculed primary vertex
  d->SetOwnPrimaryVtx(recvtx);
  delete recvtx;

  return kTRUE;
}
//--------------------------------------------------------------------------
Bool_t AliRDHFCuts::SetMCPrimaryVtx(AliAODRecoDecayHF *d,AliAODEvent *aod) const
{
  //
  // Recalculate primary vertex without daughters
  //

  if(!aod) {
    AliError("Can not get MC vertex without AOD event");
    return kFALSE;
  }

  // load MC header
  AliAODMCHeader *mcHeader =
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    AliError("Can not get MC vertex without AODMCHeader event");
    return kFALSE;
  }
  Double_t pos[3];
  Double_t covmatrix[6]={0.,0.,0.,0.,0.,0.};
  mcHeader->GetVertex(pos);
  AliAODVertex *recvtx=new AliAODVertex(pos,covmatrix);

  if(!recvtx){
    AliDebug(2,"Removal of daughter tracks failed");
    return kFALSE;
  }

  //set recalculed primary vertex
  d->SetOwnPrimaryVtx(recvtx);

  d->RecalculateImpPars(recvtx,aod);

  delete recvtx;

  return kTRUE;
}
//--------------------------------------------------------------------------
void AliRDHFCuts::CleanOwnPrimaryVtx(AliAODRecoDecayHF *d,
				     AliAODEvent *aod,
				     AliAODVertex *origownvtx) const
{
  //
  // Clean-up own primary vertex if needed
  //

  if(fRemoveDaughtersFromPrimary || fUseMCVertex) {
    d->UnsetOwnPrimaryVtx();
    if(origownvtx) {
      d->SetOwnPrimaryVtx(origownvtx);
      delete origownvtx; origownvtx=NULL;
    }
    d->RecalculateImpPars(d->GetPrimaryVtx(),aod);
  } else {
    if(origownvtx) {
      delete origownvtx; origownvtx=NULL;
    }
  }
  return;
}
//--------------------------------------------------------------------------
Bool_t AliRDHFCuts::IsSignalMC(AliAODRecoDecay *d,AliAODEvent *aod,Int_t pdg) const
{
  //
  // Checks if this candidate is matched to MC signal
  //

  if(!aod) return kFALSE;

  // get the MC array
  TClonesArray *mcArray = (TClonesArray*)((AliAODEvent*)aod)->GetList()->FindObject(AliAODMCParticle::StdBranchName());

  if(!mcArray) return kFALSE;

  // try to match
  Int_t label = d->MatchToMC(pdg,mcArray);

  if(label>=0) {
    //printf("MATCH!\n");
    return kTRUE;
  }

  return kFALSE;
}


//--------------------------------------------------------------------------
Bool_t AliRDHFCuts::RecomputePrimaryVertex(AliAODEvent* event) const{
  // recompute event primary vertex from AOD tracks

   AliVertexerTracks *vertexer = new AliVertexerTracks(event->GetMagneticField());
   vertexer->SetITSMode();
   vertexer->SetMinClusters(3);

   AliAODVertex* pvtx=event->GetPrimaryVertex();
   if(strstr(pvtx->GetTitle(),"VertexerTracksWithConstraint")) {
     Float_t diamondcovxy[3];
     event->GetDiamondCovXY(diamondcovxy);
     Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
     Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
     AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
     vertexer->SetVtxStart(diamond);
     delete diamond; diamond=NULL;
   }

   AliESDVertex* vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(event);
   if(!vertexESD) return kFALSE;
   if(vertexESD->GetNContributors()<=0) {
     //AliDebug(2,"vertexing failed");
     delete vertexESD; vertexESD=NULL;
     return kFALSE;
   }
   delete vertexer; vertexer=NULL;

   // convert to AliAODVertex
   Double_t pos[3],cov[6],chi2perNDF;
   vertexESD->GetXYZ(pos); // position
   vertexESD->GetCovMatrix(cov); //covariance matrix
   chi2perNDF = vertexESD->GetChi2toNDF();
   delete vertexESD; vertexESD=NULL;

   pvtx->SetPosition(pos[0],pos[1],pos[2]);
   pvtx->SetChi2perNDF(chi2perNDF);
   pvtx->SetCovMatrix(cov);

   return kTRUE;
}
