/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//  MC ONGOING CHANGES:
// - EXPLORE PID CASES: MUST MATCH THE RIGHT MASS HYPO!! CHECK
//  - CHECK WHAT HAPPENS FOR RESONANT CHANNELS (Lc ones + Xic to K*0: should be ok
//  - DEVELOP MC PART: CHECK SIGMA_C (all histos at reco level should be filled inside SigmaC loop, only if also the partner pion is found --> should be fixed by the addition of extra steps before calling SigmaCloop and specifc flagging insied); 
// FILLING FOR SPARSES AND TREE: CANDIDATE PT SHOULD BE THAT AT GEN LEVEL: should be ok now
//
//  SHOULD GenACCEPTANCE BE CHECKED WITH GEN VARIABLES AT RECO LEVEL?  YES: inefficiency not balanced by "contamination"  NO: NO MORE JUST AN EFFICIENCY (but likely feed-in and feed-out compensate). Relevant cut is eta one --> If GenAcc cut on daughter is defined at 0.9 and at reco one cuts at 0.8 we practically do not have this problem.
// 
// fhistMCSpectrumAccLc, fhistMCSpectrumAccSc, fhistMCSpectrumAccLcFromSc: thye do not contain any info about a kRecoPID step with the cut-object PID or ExplorePID. The kRecoCuts step include the filtering PID, but for adding an extra useful kRecoPID step one would need to add an axis for the PID strategy, on which one should project a posteriori. 
// IF INSTEAD ONE LOOPS OVER FILTERED CANDIDATES kRecoCuts include the cut-object PID but not the 3sigma PID applied in the task, while kRecoPID includes this. However, also in this case the info about ExplorePID options is not stored.
// --> EFFICIENCY CALCULATION: SHOULD BE DONE BY PROJECTING THE SPARSE to get the RECOPID LEVEL CORRESPONDING TO A SPECIFIC CUT & PID SELECTION! The GenAcc, LimAcc or Gen level instead must be taken from fhistMCSpectrumAccLc, fhistMCSpectrumAccSc, fhistMCSpectrumAccLcFromSc, fhistMCSpectrumAccXic
//
//  FIX SETTING OF CUT OBJECT IN CONSTRUCTOR.. also for EVENT ELECTION AND TRIGGER
//  CHECK AND FIX USAGE OF NORMALISATION COUNTER
// DELETION OF SECONDARY VTX
// SELCTION WITH CUTS TO BE FIXED FOR CASE 0: SHOULD NOT DELETE, WE MUST ADD A BIT TO THE HISTOGRAM OR DUPLICATE THE HISTO FOR INCLUDING XIC
// FIDUCIAL ACC: FOR XIC the Lc is used
// RECALCULATION OF PRIM VTX SHOULD NOT BE DEFINED BY FLAG IN TASK AS IT IS NOW
// CUT OBJECTS USED IN ANALYSIS SHOULD ALSO BE STREAMED 
// CHECK FILLING OF RECO STEP FOR EFFICIENCY CALCULATION, THERE MIGHT BE A BIAS RELATED TO THE MASS SELECTION, which is not checked to fill that histo
// ADD MASS CUT AROUND SIGMA_C MASS: IS IT CORRECT?
// EVENT SELECTION GEN STEP: CENTRALITY TO BE ADDED
// CHANGE THE LOOP STRUCTURE: SINCE CANDIDATES ARE DELETED, IT IS CONVENIENT TO INTRODUCE A ~static pointer TO AN ALIAODRECODDECAYHF3PRONG OBJECT, A "NEXT" METHOD THAT FILLS IT AND A RESET METHOD THAT CLEANS IT, AVOIDING NEW AND DELETE  
////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TArrayI.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TDatabasePDG.h>
#include <THnSparse.h>
#include "AliVertexingHFUtils.h"
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"

#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsXictopKpi.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEXicTopKpi.h"
#include "AliNormalizationCounter.h"
#include "AliPIDResponse.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEXicTopKpi);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEXicTopKpi::AliAnalysisTaskSEXicTopKpi():
  AliAnalysisTaskSE(),
  fvHF(0x0),
  fCuts(0x0),
  //fCutsLc(0x0),
  fCutsXic(0x0),
  fCounter(0),
  fPidResponse(0),
  fmcArray(0x0),
  fReadMC(kFALSE),
  fAnalysisType(0),
  fRecalPrimVtx(kFALSE),
  fNentries(0),
  fhistMonitoring(0x0),
  ftrackArraySel(0),
  ftrackArraySelSoftPi(0x0),
  fnSel(),
  fnSelSoftPi(),
  ftrackSelStatusProton(),
  ftrackSelStatusKaon(),
  ftrackSelStatusPion(),
  fSys(0),
  fAODProtection(0),
  fLikeSign(0),
  fESDtrackCutsProton(0x0),
  fESDtrackCutsKaon(0x0),
  fESDtrackCutsPion(0x0),
  fESDtrackCutsSoftPion(0x0),
  fprimVtx(0x0),
  fhistInvMassCheck(0x0),
  fhistMCSpectrumAccLc(0x0),
  fhistMCSpectrumAccLcFromSc(0x0),
  fhistMCSpectrumAccSc(0x0),
  fhistMCSpectrumAccXic(0x0),
  fhistMCSpectrumAccCdeuteron(0x0),
  fhSparseAnalysis(0x0),
  fhSparseAnalysisSigma(0x0),
  fhSparsePartReco(0x0),
  fhSparsePartGen(0x0),
  fCosPointDistrAll(0x0),
  fCosPointDistrAllFilter(0x0),
  fCosPointDistrSignal(0x0),
  fCosPointDistrSignalFilter(0x0),
  fDist12Signal(0x0),
  fDist12SignalFilter(0x0),
  fDist12All(0x0),
  fDist12AllFilter(0x0),
  fDist23Signal(0x0),
  fDist23All(0x0),
  fDist23AllFilter(0x0),
  fVtxResXPt(0x0),
  fVtxResYPt(0x0),
  fVtxResZPt(0x0),
  fVtxResXYPt(0x0),
  fVtxResXYZPt(0x0),
  fPrimVtxResXPt(0x0),
  fPrimVtxResYPt(0x0),
  fPrimVtxResZPt(0x0),
  fDecayLResXPt(0x0),
  fDecayLResYPt(0x0),
  fDecayLResZPt(0x0),
  fDecayLResXYPt(0x0),
  fDecayLResXYZPt(0x0),
  fnSigmaPIDtofProton(0x0),
  fnSigmaPIDtofPion(0x0),
  fnSigmaPIDtofKaon(0x0),
  fnSigmaPIDtpcProton(0x0),
  fnSigmaPIDtpcPion(0x0),
  fnSigmaPIDtpcKaon(0x0),
  fProtonID(0x0),
  fKaonID(0x0),
  fPionID(0x0),
  fOutput(0x0),
  fVertexerTracks(0x0),
  fSetTrackCutLcFilteringPP(kFALSE),
  fCutSelLevel(0),
  fApplykFirst(kFALSE),
  fMaxPtTrackkFirst(0.),
  fMaxVtxChi2Cut(10000.),
  //fFillTree(0),
  fFillTree(kFALSE),
  fTreeVar(0x0),
  fpT_down(-1),
  fLowpT_down(-1),
  fHighpT_down(-1),
  fminpT_treeFill(2.),
  fmaxpT_treeFill(36.),
  fCompute_dist12_dist23(kFALSE),
  fExplore_PIDstdCuts(kFALSE),
  fExplPID_BayesOnlyProt(kFALSE),
  fOnlyBayesPIDbin(kFALSE),
  fLcMassWindowForSigmaC(0.030),
  fSigmaCDeltaMassWindow(0.230),
  fSigmaCfromLcOnTheFly(kTRUE),
  fCheckOnlyTrackEfficiency(kFALSE),
  fIsCdeuteronAnalysis(kFALSE),
  fIsKeepOnlyCdeuteronSignal(kFALSE),
  fIsXicUpgradeAnalysis(kFALSE),
  fIsKeepOnlySigXicUpgradeAnalysis(kFALSE),
  fIsKeepOnlyBkgXicUpgradeAnalysis(kFALSE),
  fRejFactorBkgUpgrade(0.05),
  fNRotations(0),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fPdgFiducialYreco(4122)
  ,flowMass_treeFill(2.)
  ,fhighMass_tree_Fill(2.7)
  ,fStudyScPeakMC(kFALSE)
  ,fhsparseMC_ScPeak(0x0)
  ,fMinPtSoftPion(0.05)
  ,fZvtx_gen_within10_MC(0)
  ,fZvtx_gen_noSel10_MC(0)
  ,fZvtx_reco_noSel10_MC(0)
  ,fCandCounter(0x0)
  ,fCandCounter_onTheFly(0x0)
  ,fNSigmaPreFilterPID(3.)
  ,fApplyEvSel(kTRUE)
  ,fNoStdPIDcases(kFALSE)
  ,fPtSoftPionCand(0)
  ,fPtSoftPionCand_insideScLoop(0)
{
  /// Default constructor

}

//________________________________________________________________________
AliAnalysisTaskSEXicTopKpi::AliAnalysisTaskSEXicTopKpi(const char *name,AliRDHFCutsD0toKpi* cuts):
  AliAnalysisTaskSE(name),
  fvHF(0x0),
  fCuts(0x0),
  //fCutsLc(0x0),
  fCutsXic(0x0),
  fCounter(0),
  fPidResponse(0),
  fmcArray(0x0),
  fReadMC(kFALSE),
  fAnalysisType(0),
  fRecalPrimVtx(kFALSE),
  fNentries(0),
  fhistMonitoring(0x0),
  ftrackArraySel(0),
  ftrackArraySelSoftPi(0x0),
  fnSel(),
  fnSelSoftPi(),
  ftrackSelStatusProton(),
  ftrackSelStatusKaon(),
  ftrackSelStatusPion(),
  fSys(0),
  fAODProtection(0),
  fLikeSign(0),
  fESDtrackCutsProton(0x0),
  fESDtrackCutsKaon(0x0),
  fESDtrackCutsPion(0x0),
  fESDtrackCutsSoftPion(0x0),
  fprimVtx(0x0),
  fhistInvMassCheck(0x0),
  fhistMCSpectrumAccLc(0x0),
  fhistMCSpectrumAccLcFromSc(0x0),
  fhistMCSpectrumAccSc(0x0),
  fhistMCSpectrumAccXic(0x0),
  fhistMCSpectrumAccCdeuteron(0x0),
  fhSparseAnalysis(0x0),
  fhSparseAnalysisSigma(0x0),
  fhSparsePartReco(0x0),
  fhSparsePartGen(0x0),
  fCosPointDistrAll(0x0),
  fCosPointDistrAllFilter(0x0),
  fCosPointDistrSignal(0x0),
  fCosPointDistrSignalFilter(0x0),
  fDist12Signal(0x0),
  fDist12SignalFilter(0x0),
  fDist12All(0x0),
  fDist12AllFilter(0x0),
  fDist23Signal(0x0),
  fDist23All(0x0),
  fDist23AllFilter(0x0),
  fVtxResXPt(0x0),
  fVtxResYPt(0x0),
  fVtxResZPt(0x0),
  fVtxResXYPt(0x0),
  fVtxResXYZPt(0x0),
  fPrimVtxResXPt(0x0),
  fPrimVtxResYPt(0x0),
  fPrimVtxResZPt(0x0),
  fDecayLResXPt(0x0),
  fDecayLResYPt(0x0),
  fDecayLResZPt(0x0),
  fDecayLResXYPt(0x0),
  fDecayLResXYZPt(0x0),
  fnSigmaPIDtofProton(0x0),
  fnSigmaPIDtofPion(0x0),
  fnSigmaPIDtofKaon(0x0),
  fnSigmaPIDtpcProton(0x0),
  fnSigmaPIDtpcPion(0x0),
  fnSigmaPIDtpcKaon(0x0),
  fProtonID(0x0),
  fKaonID(0x0),
  fPionID(0x0),
  fOutput(0x0),
  fVertexerTracks(0x0),
  fSetTrackCutLcFilteringPP(kFALSE),
  fCutSelLevel(0),
  fApplykFirst(kFALSE),
  fMaxPtTrackkFirst(0.),
  fMaxVtxChi2Cut(10000.),
  fFillTree(kFALSE),
  fTreeVar(0x0),
  fpT_down(-1),
  fLowpT_down(-1),
  fHighpT_down(-1),
  fminpT_treeFill(2.),
  fmaxpT_treeFill(36.),
  fCompute_dist12_dist23(kFALSE),
  fExplore_PIDstdCuts(kFALSE),
  fExplPID_BayesOnlyProt(kFALSE),
  fOnlyBayesPIDbin(kFALSE),
  fLcMassWindowForSigmaC(0.030),
  fSigmaCDeltaMassWindow(0.230),
  fSigmaCfromLcOnTheFly(kTRUE),
  fCheckOnlyTrackEfficiency(kFALSE),
  fIsCdeuteronAnalysis(kFALSE),
  fIsKeepOnlyCdeuteronSignal(kFALSE),
  fIsXicUpgradeAnalysis(kFALSE),
  fIsKeepOnlySigXicUpgradeAnalysis(kFALSE),
  fIsKeepOnlyBkgXicUpgradeAnalysis(kFALSE),
  fRejFactorBkgUpgrade(0.05),
  fNRotations(0),
  fMinAngleForRot(5*TMath::Pi()/6),
  fMaxAngleForRot(7*TMath::Pi()/6),
  fPdgFiducialYreco(4122)
  ,flowMass_treeFill(2.)
  ,fhighMass_tree_Fill(2.7)
  ,fStudyScPeakMC(kFALSE)
  ,fhsparseMC_ScPeak(0x0)
  ,fMinPtSoftPion(0.05)
  ,fZvtx_gen_within10_MC(0)
  ,fZvtx_gen_noSel10_MC(0)
  ,fZvtx_reco_noSel10_MC(0)
  ,fCandCounter(0x0)
  ,fCandCounter_onTheFly(0x0)
  ,fNSigmaPreFilterPID(3.)
  ,fApplyEvSel(kTRUE)
  ,fNoStdPIDcases(kFALSE)
  ,fPtSoftPionCand(0)
  ,fPtSoftPionCand_insideScLoop(0)
{
  /// Default constructor


  DefineOutput(1,TH1F::Class());  //My private output
  DefineOutput(2,AliNormalizationCounter::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,TTree::Class());
  fCuts=cuts;
}

//________________________________________________________________________
AliAnalysisTaskSEXicTopKpi::~AliAnalysisTaskSEXicTopKpi()
{
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  if(fvHF){
    delete fvHF;
    fvHF=0;
  }
  //if(fCutsLc){
  //  delete fCutsLc;
  //  fCutsLc =0;
  //}
  if(fCutsXic){
    delete fCutsXic;
    fCutsXic =0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
  if(  fhistMonitoring){
    delete   fhistMonitoring;
    fhistMonitoring=0;
  }
  if(fCounter){
    delete fCounter;
    fCounter=0;
  }
  if(ftrackArraySel)delete ftrackArraySel;
  if(ftrackArraySelSoftPi)delete ftrackArraySelSoftPi;

  if(ftrackSelStatusProton)delete ftrackSelStatusProton;
  if(ftrackSelStatusKaon)delete ftrackSelStatusKaon;
  if(ftrackSelStatusPion)delete ftrackSelStatusPion;

  if(fESDtrackCutsProton)delete fESDtrackCutsProton;
  if(fESDtrackCutsKaon) delete fESDtrackCutsKaon;
  if(fESDtrackCutsPion) delete fESDtrackCutsPion;
  if(fESDtrackCutsSoftPion) delete fESDtrackCutsSoftPion;
  if(  fhistInvMassCheck) delete   fhistInvMassCheck;
  if(fhistMCSpectrumAccLc)delete fhistMCSpectrumAccLc;
  if(fhistMCSpectrumAccLcFromSc)delete fhistMCSpectrumAccLcFromSc;
  if(fhistMCSpectrumAccSc)delete fhistMCSpectrumAccSc;
  if(fhistMCSpectrumAccXic)delete fhistMCSpectrumAccXic;
  if(fhistMCSpectrumAccCdeuteron)delete fhistMCSpectrumAccCdeuteron;
  if(fhSparseAnalysis)delete fhSparseAnalysis;
  if(fhSparseAnalysisSigma)delete fhSparseAnalysisSigma;
  if(fhSparsePartReco)delete fhSparsePartReco; 
  if(fhSparsePartGen)delete fhSparsePartGen;
  if(fCosPointDistrAll)delete fCosPointDistrAll;
  if(fCosPointDistrSignal)delete fCosPointDistrSignal;
  if(fCosPointDistrAllFilter)delete fCosPointDistrAllFilter;
  if(fCosPointDistrSignalFilter)delete fCosPointDistrSignalFilter;
  if(fDist12Signal)delete fDist12Signal;
  if(fDist12SignalFilter)delete fDist12SignalFilter;
  if(fDist12All)delete fDist12All;
  if(fDist12AllFilter)delete fDist12AllFilter;
  if(fDist23Signal)delete fDist23Signal;
  if(fDist23All)delete fDist23All;
  if(fDist23AllFilter)delete fDist23AllFilter;
  if(fVtxResXPt)delete fVtxResXPt;
  if(fVtxResYPt)delete fVtxResYPt;
  if(fVtxResZPt)delete fVtxResZPt;
  if(fVtxResXYPt)delete fVtxResXYPt;
  if(fPrimVtxResXPt)delete fPrimVtxResXPt;
  if(fPrimVtxResYPt)delete fPrimVtxResYPt;
  if(fPrimVtxResZPt)delete fPrimVtxResZPt;
  if(fDecayLResXPt)delete fDecayLResXPt;
  if(fDecayLResYPt)delete fDecayLResYPt;
  if(fDecayLResZPt)delete fDecayLResZPt;
  if(fDecayLResXYPt)delete fDecayLResXYPt;
  if(fOutput)delete fOutput;
  if(fVertexerTracks)delete fVertexerTracks;
  if(fnSigmaPIDtofProton)delete fnSigmaPIDtofProton;
  if(fnSigmaPIDtofPion)delete fnSigmaPIDtofPion;
  if(fnSigmaPIDtofKaon)delete fnSigmaPIDtofKaon;
  if(fnSigmaPIDtpcProton)delete fnSigmaPIDtpcProton;
  if(fnSigmaPIDtpcPion)delete fnSigmaPIDtpcPion;
  if(fnSigmaPIDtpcKaon)delete fnSigmaPIDtpcKaon;
  if(fProtonID)delete fProtonID;
  if(fKaonID)delete fKaonID;
  if(fPionID)delete fPionID;
  if(fTreeVar)delete fTreeVar;
  if(fhsparseMC_ScPeak)delete fhsparseMC_ScPeak;
}  

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::Init()
{
  /// Initialization

  if(fDebug > 1) printf("AnalysisTaskSEXicTopKpi::Init() \n");

  
//   AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*fCuts);
//   const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
//   copyfCuts->SetName(nameoutput);

  if(!fCutsXic || fSetTrackCutLcFilteringPP){
    Printf("Setting default cuts");
    if(fCutsXic)delete fCutsXic;
    fCutsXic = new AliRDHFCutsXictopKpi("CutsXictopKpi");
    //Float_t cutsArrayLctopKpi[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};
    Float_t cutsArrayLctopKpi[13]={0.18,0.4,0.5,0.,0.,0.01,0.06,0.005,0.7,0.0,0.,0.05,0.4};// cuts in ConfigVertexingHF_pp_LambdaC.C rot 2018pp (0.18 original on mass, changed to 0.4)
  //Float_t cutsArrayLctopKpi[13]={0.18,0.3,0.3,0.,0.,0.00,0.06,0.000,0.0,-1.,0.,0.05,0.3}; // cuts in ConfigVertexingHF.C
    //Float_t cutsArrayLctopKpi[13]={0.40,0.4,0.5,0.,0.,0.00,0.06,0.000,0.0,-1.,0.,0.05,0.4}; // default filtering cuts (from ConfigVertexingHF.C) in pp collisions with enlarged mass cut (from 0.18 to 0.400) to store also candidates in Xic mass window

    /*THE VARIABLES ARE: 
      "inv. mass [GeV]",
      "pTK [GeV/c]",
      "pTP [GeV/c]",
      "d0K [cm]   lower limit!",
      "d0Pi [cm]  lower limit!",
      "dist12 (cm)",
      "sigmavert (cm)",
      "dist prim-sec (cm)",
      "pM=Max{pT1,pT2,pT3} (GeV/c)",
      "cosThetaPoint",
      "Sum d0^2 (cm^2)",
      "dca cut (cm)",
      "cut on pTpion [GeV/c]" */
    fCutsXic->SetCuts(13,cutsArrayLctopKpi);
    //    fCutsXic->SetUsePID(kFALSE); 
    fCutsXic->PrintAll();
  }
  else{
   
//     AliESDtrackCuts *esdTrCuts=new AliESDtrackCuts("AliESDtrackCutsAll","default");
//     esdTrCuts->SetRequireTPCRefit(kTRUE);
//     esdTrCuts->SetMinNClustersTPC(50);
//     esdTrCuts->SetRequireITSRefit(kTRUE);
//     //esdTrCuts->SetMinNClustersITS(4);
//     esdTrCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
// 					   AliESDtrackCuts::kAny);
//     esdTrCuts->SetMinDCAToVertexXY(0.);
//     esdTrCuts->SetPtRange(0.3,1.e10);// ORIGINAL VALUE IS 0.3 BUT WE CUT AT 0.4 in ANALYSIS
//     esdTrCuts->SetEtaRange(-0.8,+0.8);
//     fCutsXic->AddTrackCuts(esdTrCuts);
    /*
    Double_t sigmaVtxMax[8]={0.09,0.09,0.05,0.035,0.035,0.03,0.03,0.025};
    Double_t sumd02[8]={0.0003,0.0003,0.0002,0.00015,0.00015,0.0001,0.,0.};
    if(fAnalysisType==0 || fAnalysisType==4){// assure mass range is large enough
      Int_t nvars=fCutsXic->GetNVars();
      Int_t nptbins=fCutsXic->GetNPtBins();
      Float_t **cutvalues=new Float_t*[nvars];
      for(Int_t iv=0;iv<nvars;iv++){
	cutvalues[iv]=new Float_t[nptbins];
      }
      fCutsXic->GetCuts(cutvalues);
      for(Int_t km=0;km<nptbins;km++){
	//cutvalues[0][km]=0.40;
	cutvalues[6][km]=sigmaVtxMax[km]; // TO BE CHANGED IN FUTURE
	cutvalues[10][km]=sumd02[km];
      }
      fCutsXic->SetCuts(nvars,nptbins,cutvalues);
      Printf("Xic Cuts modified to assure mass window is large enough, current cuts are:");*/
      Printf("\n---cuts set by user ---");
      fCutsXic->PrintAll();
      //}
  }
  
  if(fDebug>=0 || fSetTrackCutLcFilteringPP || (!(fCutsXic->GetTrackCuts()))){// track cuts used for Lc filtering (in pp, 2018): need to set them to be sure that only tighter cuts than these are used
    AliESDtrackCuts *esdTrCuts=new AliESDtrackCuts("AliESDtrackCutsAll","default");
    esdTrCuts->SetRequireTPCRefit(kTRUE);
    esdTrCuts->SetMinNClustersTPC(50);
    esdTrCuts->SetRequireITSRefit(kTRUE);
    //esdTrCuts->SetMinNClustersITS(4);
    esdTrCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					AliESDtrackCuts::kAny);
    esdTrCuts->SetMinDCAToVertexXY(0.);
    esdTrCuts->SetPtRange(0.3,1.e10);// ORIGINAL VALUE IS 0.3 BUT WE CUT AT 0.4 in ANALYSIS
    esdTrCuts->SetEtaRange(-0.8,+0.8);
    fCutsXic->AddTrackCuts(esdTrCuts);
  }
  
  if((fAnalysisType==0 || fAnalysisType==3) && !fESDtrackCutsSoftPion){
    fESDtrackCutsSoftPion = new AliESDtrackCuts("AliESDtrackCuts","default");
    //    fESDtrackCutsSoftPion->SetRequireITSRefit(kTRUE);   
    fESDtrackCutsSoftPion->SetMinNClustersITS(3);
    fESDtrackCutsSoftPion->SetMaxDCAToVertexXY(0.065);
    fESDtrackCutsSoftPion->SetPtRange(fMinPtSoftPion,1.e10);
    fESDtrackCutsSoftPion->SetMaxDCAToVertexZ(0.15);
    fESDtrackCutsSoftPion->SetEtaRange(-0.9,+0.9);    
  } 

  
  // protection against negative values for downsampling for fTreeVar filling
  if(fFillTree){
    if(fpT_down<0 || fLowpT_down<0 || fHighpT_down<0){
      if(fpT_down<0)      fpT_down     = 4;
      if(fLowpT_down<0)   fLowpT_down  = 0.005;
      if(fHighpT_down<0)  fHighpT_down = 0.05;
      printf("\n\n=== fTreeVar filling activated without downsampling parameters ===\n");
      printf("    Using downsampling rescue values.\n    Adopted values:\n");
      printf("      fpT_down = %f\n      fLowpT_down = %f\n      fHighpT_down = %f\n\n\n",fpT_down,fLowpT_down,fHighpT_down);
    }
  }
  
  // print wheter we require the calulcation of dist12 and dist23
  if(fCompute_dist12_dist23)  printf("\n\n===== fCompute_dist12_dist23 is kTRUE ---> dist12 and dist23 will be calculated (->increase CPU time)=====\n\n");

  // Post the data
  //  PostData(4,copyfCuts);


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::UserCreateOutputObjects()
{

  /// Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEXicTopKpi::UserCreateOutputObjects() \n");
  if(fRecalPrimVtx && fSys>0)AliFatal("Cannot recalculcate prim vtx in p-Pb and Pb-Pb");
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("listOutput");

  // histogram to count candidates
  fCandCounter = new TH1F("fCandCounter","candidate counter (pre-filtered)",12,-1.5,10.5);
  fCandCounter->GetXaxis()->SetBinLabel(1,"evs. entering cand. loop");
  fCandCounter->GetXaxis()->SetBinLabel(2,"filtered cand. loop entered");
  fCandCounter->GetXaxis()->SetBinLabel(3,"Lc cuts passed (filbit & QA)");
  //fCandCounter->GetXaxis()->SetBinLabel(3,"MC matched candidate");
  fCandCounter->GetXaxis()->SetBinLabel(4,"after FillRecoCand");
  fCandCounter->GetXaxis()->SetBinLabel(5,"after IsSelectedTracks");
  //fCandCounter->GetXaxis()->SetBinLabel(5,"after FillRecoCand (signal)");
  fCandCounter->GetXaxis()->SetBinLabel(6,"in fiducial acceptance");
  fCandCounter->GetXaxis()->SetBinLabel(7,"after IsSelectedTracks and iSelPID");
  //fCandCounter->GetXaxis()->SetBinLabel(7,"in fiducial acceptance (signal)");
  //fCandCounter->GetXaxis()->SetBinLabel(8,"before IsSelected");
  //fCandCounter->GetXaxis()->SetBinLabel(9,"before IsSelected (signal)");
  fCandCounter->GetXaxis()->SetBinLabel(8,"before IsSelected");
  fCandCounter->GetXaxis()->SetBinLabel(9,"after IsSelected");
  fOutput->Add(fCandCounter);

  // histogram to count candidates on the fly
  fCandCounter_onTheFly = new TH1F("fCandCounter_onTheFly","candidate counter (on-the-fly)",10,0.5,10.5);
  fCandCounter_onTheFly->GetXaxis()->SetBinLabel(1,"evs. entering the triple nested track loop");
  fCandCounter_onTheFly->GetXaxis()->SetBinLabel(2,"on-the-fly candidate built");
  fCandCounter_onTheFly->GetXaxis()->SetBinLabel(3,"after FillRecoCand");
  fCandCounter_onTheFly->GetXaxis()->SetBinLabel(4,"in fiducial acceptance");
  fCandCounter_onTheFly->GetXaxis()->SetBinLabel(5,"after reducedChi2 cut");
  fCandCounter_onTheFly->GetXaxis()->SetBinLabel(6,"candidates after only cuts");
  fCandCounter_onTheFly->GetXaxis()->SetBinLabel(7,"candidates after cuts and PID in cutobject");
  fCandCounter_onTheFly->GetXaxis()->SetBinLabel(8,"candidates effectively used (beware fExplore_PIDstdCuts)");
  fOutput->Add(fCandCounter_onTheFly);

  // histograms to monitor the vtx_z
  fZvtx_gen_within10_MC = new TH1D("fZvtx_gen_within10_MC","vtx_{z}^{gen.} requiring to be < 10 cm",200,-20,20);
  fZvtx_gen_noSel10_MC  = new TH1D("fZvtx_gen_noSel10_MC","vtx_{z}^{gen.} after event selection",200,-20,20);
  fZvtx_reco_noSel10_MC = new TH1D("fZvtx_reco_noSel10_MC","vtx_{z}^{reco} after event selection",200,-20,20);
  fOutput->Add(fZvtx_gen_within10_MC);
  fOutput->Add(fZvtx_gen_noSel10_MC);
  fOutput->Add(fZvtx_reco_noSel10_MC);


  fhistMonitoring=new TH1F("fhistMonitoring","Overview",30,-0.5,29.5);
  fhistMonitoring->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fhistMonitoring->GetXaxis()->SetBinLabel(2,"nEventsSel");
  fhistMonitoring->GetXaxis()->SetBinLabel(3,"nEvGoodVtx");
  fhistMonitoring->GetXaxis()->SetBinLabel(4,"nTracksTot");
  fhistMonitoring->GetXaxis()->SetBinLabel(5,"nTracksSel");
  fhistMonitoring->GetXaxis()->SetBinLabel(6,"nPionComp");
  fhistMonitoring->GetXaxis()->SetBinLabel(7,"nKaonComp");
  fhistMonitoring->GetXaxis()->SetBinLabel(8,"nProtonComp");
  fhistMonitoring->GetXaxis()->SetBinLabel(9,"nTriplets");
  fhistMonitoring->GetXaxis()->SetBinLabel(10,"nCandidates");
  fhistMonitoring->GetXaxis()->SetBinLabel(11,"nCandSel");

  fOutput->Add(fhistMonitoring);

  TString strnameoutput="EntryCounter";
  fNentries=new TH1F(strnameoutput.Data(), "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of D0 selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 23,-0.5,22.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  if(fReadMC) fNentries->GetXaxis()->SetBinLabel(3,"nD0Selected");
  else fNentries->GetXaxis()->SetBinLabel(3,"Dstar<-D0");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
  if(fSys==0) fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
  //  if(fFillVarHists || fPIDCheck){
  fNentries->GetXaxis()->SetBinLabel(8,"PID=0");
  fNentries->GetXaxis()->SetBinLabel(9,"PID=1");
  fNentries->GetXaxis()->SetBinLabel(10,"PID=2");
  fNentries->GetXaxis()->SetBinLabel(11,"PID=3");
  //  }
  if(fReadMC && fSys==0){
    fNentries->GetXaxis()->SetBinLabel(12,"K");
    fNentries->GetXaxis()->SetBinLabel(13,"Lambda");
  }
  fNentries->GetXaxis()->SetBinLabel(14,"Pile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(15,"N. of 0SMH");
  if(fSys==1) fNentries->GetXaxis()->SetBinLabel(16,"Nev in centr");
  //  if(fIsRejectSDDClusters)
  fNentries->GetXaxis()->SetBinLabel(17,"SDD-Cls Rej");
  fNentries->GetXaxis()->SetBinLabel(18,"Phys.Sel.Rej");
  fNentries->GetXaxis()->SetBinLabel(19,"D0 failed to be filled");
  fNentries->GetXaxis()->SetBinLabel(20,"fisFilled is 0");
  fNentries->GetXaxis()->SetBinLabel(21,"fisFilled is 1");
  fNentries->GetXaxis()->SetBinLabel(22,"AOD/dAOD mismatch");
  fNentries->GetXaxis()->SetBinLabel(23,"AOD/dAOD #events ok");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(2)->GetContainer()->GetName()));
  fCounter->Init();


  ftrackArraySel=new TArrayI(10000);
  ftrackArraySelSoftPi=new TArrayI(10000);
  ftrackSelStatusProton=new TArrayI(10000);
  ftrackSelStatusKaon=new TArrayI(10000);
  ftrackSelStatusPion=new TArrayI(10000);

  
  //fhistInvMassCheck=new TH2F("fhistInvMassCheck","InvDistrCheck",1000,1.600,2.800,5,-0.5,4.5);
  if(fIsCdeuteronAnalysis) fhistInvMassCheck=new TH2F("fhistInvMassCheck","InvDistrCheck",1000,2.600,3.800,1010,-0.5,1010);
  else fhistInvMassCheck=new TH2F("fhistInvMassCheck","InvDistrCheck",1000,1.600,2.800,5,-0.5,15.5);
  //  fhistCheckPIDTOFTPC=new TH3F("fhistCheckPIDTOFTPC","fhistCheckPIDTOFTPC",

  //fhistMCSpectrumAccLc=new TH3F("fhistMCSpectrumAccLc","fhistMCSpectrumAccLc",250,0,50,20,-0.5,19.5,2,3.5,5.5); // 

  const Int_t nAxes_THnMC = 4;
  Int_t nBins_THnMC[nAxes_THnMC]  = {250,   20,   2,    6};
  Double_t min_THnMC[nAxes_THnMC] = {0  , -0.5, 3.5, -1.5};
  Double_t max_THnMC[nAxes_THnMC] = {50 , 19.5, 5.5,  4.5};
  fhistMCSpectrumAccLc=new THnF("fhistMCSpectrumAccLc","fhistMCSpectrumAccLc;pt;step;origin;decay_channel;",nAxes_THnMC,nBins_THnMC,min_THnMC,max_THnMC); // 

  // adding axis for Lc decay channel (MC)
  const Int_t nbinsAccLcFromSc=7;
  Int_t binsAccLcFromSc[nbinsAccLcFromSc]        = {250,  20,   2, 20, 250, 40,    6};
  Double_t lowedgesAccLcFromSc[nbinsAccLcFromSc] = {0  ,-0.5, 3.5, -1,   0, -2, -1.5};
  Double_t upedgesAccLcFromSc[nbinsAccLcFromSc]  = {50 ,19.5, 5.5,  1,  50,  2,  4.5};
  fhistMCSpectrumAccLcFromSc=new THnSparseF("fhistMCSpectrumAccLcFromSc","fhistMCSpectrumAccLcFromSc;ptLc;codeLc;Qorigin;yLc;ptSc;ySc;decay_channel;",nbinsAccLcFromSc,binsAccLcFromSc,lowedgesAccLcFromSc,upedgesAccLcFromSc); // 

  //fhistMCSpectrumAccXic=new TH3F("fhistMCSpectrumAccXic","fhistMCSpectrumAccXic",250,0,50,20,-0.5,19.5,2,3.5,5.5); // 
  fhistMCSpectrumAccXic=new THnF("fhistMCSpectrumAccXic","fhistMCSpectrumAccXic;pt;step;origin;decay_channel;",nAxes_THnMC,nBins_THnMC,min_THnMC,max_THnMC); // 
  
  //fhistMCSpectrumAccSc=new TH3F("fhistMCSpectrumAccSc","fhistMCSpectrumAccSc",250,0,50,20,-0.5,19.5,2,3.5,5.5); // 
  fhistMCSpectrumAccSc=new THnF("fhistMCSpectrumAccSc","fhistMCSpectrumAccSc;pt;step;origin;decay_channel;",nAxes_THnMC,nBins_THnMC,min_THnMC,max_THnMC); // 
  
  fhistMCSpectrumAccCdeuteron=new TH2F("fhistMCSpectrumAccCdeuteron","fhistMCSpectrumAccCdeuteron",250,0,50,20,-0.5,19.5); //

  // Sparse histos to study track reco & PID efficiency
  Int_t nbinsSparseTrack[9]={100,20,16,3,3,3,3,3,4};
  Double_t lowEdgesSparseTrack[9]={0.,-1,0,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5};
  Double_t upEdgesSparseTrack[9]={2.5,1,TMath::Pi()*2.,2.5,2.5,2.5,2.5,2.5,3.5};
  if(fReadMC && !fFillTree)  fhSparsePartReco=new THnSparseF("fhSparsePartReco","fhSparsePartReco;pt:eta:phi:piSel:softPiSel:Ksel:Psel:partType",8,nbinsSparseTrack,lowEdgesSparseTrack,upEdgesSparseTrack);
  Int_t nbinsSparsePart[4]={100,20,16,4};//same as above but at kine level
  Double_t lowEdgesSparsePart[4]={0.,-1,0,-0.5};
  Double_t upEdgesSparsePart[4]={2.5,1,TMath::Pi()*2.,3.5};
  if(fReadMC && !fFillTree)  fhSparsePartGen=new THnSparseF("fhSparsePartGen","fhSparsePartGen;pt:eta:phi:piSel:Ksel:Psel:partType",4,nbinsSparsePart,lowEdgesSparsePart,upEdgesSparsePart);


  //Int_t nbinsSparse[7]={16,125,10,16,20,10,10};
  //Double_t lowEdges[7]={0,2.15,0.,0,0.8,0,-1};
  //Double_t upEdges[7]={16,2.65,0.0500,8,1.,5,9};
  //if(!fFillTree)  fhSparseAnalysis=new THnSparseF("fhSparseAnalysis","fhSparseAnalysis;pt;mass;Lxy;nLxy;cosThatPoint;normImpParXY;seleFlag",7,nbinsSparse,lowEdges,upEdges);
  //Int_t nbinsSparseSigma[9]={16,200,10,12,10,10,10,22,20};
  //Double_t lowEdgesSigma[9]={0,0.130,0.,0,0.8,0,-1,2.262,-1};
  //Double_t upEdgesSigma[9]={16,0.330,0.0500,6.,1.,5,9,2.306,1};
  //if(!fFillTree)  fhSparseAnalysisSigma=new THnSparseF("fhSparseAnalysisSigma","fhSparseAnalysis;pt;deltamass;Lxy;nLxy;cosThetaPoint;normImpParXY;seleFlag;LcMass;CosThetaStarSoftPion",9,nbinsSparseSigma,lowEdgesSigma,upEdgesSigma);
  
  // adding PID cases study and axis for decay channel (MC)
  Int_t nbinsSparse[9]={24,125,10,16,20,10,23,11,6};
  Double_t lowEdges[9]={0,2.15,0.,0,0.8,0,-0.5,-0.5,-1.5};
  Double_t upEdges[9]={24,2.65,0.0500,8,1.,5,22.5,10.5,4.5};
  if(fIsCdeuteronAnalysis){
    lowEdges[1] = 2.95;
    upEdges[1] = 3.45;
  }
  if(fExplore_PIDstdCuts && fOnlyBayesPIDbin){
    printf("\n##############################################################################\n");
    printf("ATTENTION: bins for PID axis in reco THnSparse for Lc (Xic) reduced from 11 to 1\n");
    printf("##############################################################################\n");
    // only Bayes PID
    nbinsSparse[7]=1;
    upEdges[7]=0.5;
  }
  if(fExplore_PIDstdCuts && fExplPID_BayesOnlyProt){  // add the bin for the Bayes PID case only for proton hypothesis
    printf("\n#################################################################################\n");
    printf("ATTENTION: bins for PID axis in reco THnSparse for Lc (Xic) increased from 11 to 12\n");
    printf("#################################################################################\n");
    nbinsSparse[7]=12;
    upEdges[7]=11.5;
  }
  if(fIsXicUpgradeAnalysis){
    // bins of 0.002
    nbinsSparse[4] = 50;  // finer bins for cosThPoint
    lowEdges[4] = 0.9;
  }
  if(!fFillTree)  fhSparseAnalysis=new THnSparseF("fhSparseAnalysis","fhSparseAnalysis;pt;mass;Lxy;nLxy;cosThatPoint;normImpParXY;infoMC;PIDcase;channel",9,nbinsSparse,lowEdges,upEdges);
  
  // add also here the axis for Lc decay channel (MC)
  Int_t nbinsSparseSigma[14]={16,400,10,12,10,10,1,11,22,20,16,2,1,6};
  Double_t lowEdgesSigma[14]={0,0.130,0.,0,0.8,0,-0.5,-0.5,2.266,-1,0,3.5,0.5,-1.5};
  Double_t upEdgesSigma[14]={16,0.330,0.0500,6.,1.,5,0.5,10.5,2.306,1,16,5.5,1.5,4.5};
  if(fExplore_PIDstdCuts && fOnlyBayesPIDbin){
    printf("\n#########################################################################################\n");
    printf("ATTENTION: bins for PID axis in reco THnSparse for Lc(<-Sc) and Sc reduced from 11 to 1\n");
    printf("#########################################################################################\n");
    // only Bayes PID
    nbinsSparseSigma[7]=1;
    upEdgesSigma[7]=0.5;
  }
  if(fExplore_PIDstdCuts && fExplPID_BayesOnlyProt){  // add the bin for the Bayes PID case only for proton hypothesis
    printf("\n#################################################################################\n");
    printf("ATTENTION: bins for PID axis in reco THnSparse for SigmaC increased from 11 to 12\n");
    printf("#################################################################################\n");
    nbinsSparseSigma[7]=12;
    upEdgesSigma[7]=11.5;
  }
  if(fReadMC){
    // save the generated pT for reco particles with finer binning
    nbinsSparseSigma[0]=80;
    nbinsSparseSigma[10]=80;
  }
  if(!fFillTree)  fhSparseAnalysisSigma=new THnSparseF("fhSparseAnalysisSigma","fhSparseAnalysis;pt;deltamass;Lxy;nLxy;cosThetaPoint;normImpParXY;softPiITSrefit;PIDcase;LcMass;CosThetaStarSoftPion;ptsigmac;checkorigin;isRotated;channel",14,nbinsSparseSigma,lowEdgesSigma,upEdgesSigma);
  
  fCosPointDistrAll=new TH1F("fCosPointDistrAll","fCosPointDistrAll",200,-1.1,1.1);
  fCosPointDistrSignal=new TH1F("fCosPointDistrSignal","fCosPointDistrSignal",200,-1.1,1.1);
  fCosPointDistrAllFilter=new TH1F("fCosPointDistrAllFilter","fCosPointDistrAllFilter",200,-1.1,1.1);
  fCosPointDistrSignalFilter=new TH1F("fCosPointDistrSignalFilter","fCosPointDistrSignalFilter",200,-1.1,1.1);
  fDist12Signal=new TH1F("fDist12Signal","fDist12Signal",500,0.,1000.);
  fDist12SignalFilter=new TH1F("fDist12SignalFilter","fDist12SignalFilter",500,0.,1000.);
  fDist23Signal=new TH1F("fDist23Signal","fDist23Signal",500,0.,1000.);
  fDist12All=new TH1F("fDist12All","fDist12All",500,0.,1000.);
  fDist12AllFilter=new TH1F("fDist12AllFilter","fDist12AllFilter",500,0.,1000.);
  fDist23All=new TH1F("fDist23All","fDist23All",500,0.,1000.);
  fDist23AllFilter=new TH1F("fDist23AllFilter","fDist23AllFilter",500,0.,1000.);
  fVtxResXPt = new TH2F("fVtxResXPt","vertex resolution X vs p_{T}",1000,-0.1,0.1,200,0,20);
  fVtxResYPt = new TH2F("fVtxResYPt","vertex resolution Y vs p_{T}",1000,-0.1,0.1,200,0,20);
  fVtxResZPt = new TH2F("fVtxResZPt","vertex resolution Z vs p_{T}",1000,-0.1,0.1,200,0,20);
  fVtxResXYPt = new TH2F("fVtxResXYPt","vertex resolution XY vs p_{T}",1000,-0.1,0.1,200,0,20);
  fVtxResXYZPt = new TH2F("fVtxResXYZPt","vertex resolution XYZ vs p_{T}",1000,-0.1,0.1,200,0,20);
  fPrimVtxResXPt = new TH2F("fPrimVtxResXPt","primary vertex resolution X vs p_{T}",1000,-0.1,0.1,200,0,20);
  fPrimVtxResYPt = new TH2F("fPrimVtxResYPt","primary vertex resolution Y vs p_{T}",1000,-0.1,0.1,200,0,20);
  fPrimVtxResZPt = new TH2F("fPrimVtxResZPt","primary vertex resolution Z vs p_{T}",1000,-0.1,0.1,200,0,20);
  fDecayLResXPt = new TH2F("fDecayLResXPt","decay length resolution X vs p_{T}",1000,-0.1,0.1,200,0,20);
  fDecayLResYPt = new TH2F("fDecayLResYPt","decay length resolution Y vs p_{T}",1000,-0.1,0.1,200,0,20);
  fDecayLResZPt = new TH2F("fDecayLResZPt","decay length resolution Z vs p_{T}",1000,-0.1,0.1,200,0,20);
  fDecayLResXYPt = new TH2F("fDecayLResXYPt","decay length resolution XY vs p_{T}",1000,-0.1,0.1,200,0,20);
  fDecayLResXYZPt = new TH2F("fDecayLResXYZPt","decay length resolution XYZ vs p_{T}",1000,-0.1,0.1,200,0,20);

  fnSigmaPIDtofProton=new TH2F("fnSigmaPIDtofProton","fnSigmaPIDtofProton",100,0,20,80,-10,10);
  fnSigmaPIDtofPion=new TH2F("fnSigmaPIDtofPion","fnSigmaPIDtofPion",100,0,20,80,-10,10);
  fnSigmaPIDtofKaon=new TH2F("fnSigmaPIDtofKaon","fnSigmaPIDtofKaon",100,0,20,80,-10,10);
  fnSigmaPIDtpcProton=new TH2F("fnSigmaPIDtpcProton","fnSigmaPIDtpcProton",100,0,20,80,-10,10);
  fnSigmaPIDtpcPion=new TH2F("fnSigmaPIDtpcPion","fnSigmaPIDtpcPion",100,0,20,80,-10,10);
  fnSigmaPIDtpcKaon=new TH2F("fnSigmaPIDtpcKaon","fnSigmaPIDtpcKaon",100,0,20,80,-10,10);

  fProtonID=new TH2F("fProtonID","fProtonID",2,0,2,200,0,20);
  fProtonID->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fProtonID->GetXaxis()->SetBinLabel(1,"MC particles");
  fProtonID->GetXaxis()->SetBinLabel(2,"ID'd particles");
  fKaonID=new TH2F("fKaonID","fKaonID",2,0,2,200,0,20);
  fKaonID->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fKaonID->GetXaxis()->SetBinLabel(1,"MC particles");
  fKaonID->GetXaxis()->SetBinLabel(2,"ID'd particles");
  fPionID=new TH2F("fPionID","fPionID",2,0,2,200,0,20);
  fPionID->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fPionID->GetXaxis()->SetBinLabel(1,"MC particles");
  fPionID->GetXaxis()->SetBinLabel(2,"ID'd particles");

/*
  //  pt vs. pointing angle, lxy, nlxy, ptP,ptK,ptPi,vtxchi2,sigmaVtx,sumd02,dca1,dca2,dca3,nd01,nd02,nd03,Lc d0
  //  Float_t pt,pAngle,lxy,nlxy,ptP,ptK,ptPi,vtxchi2,sigmaVtx,sumd02,dca1,dca2,dca3,nd01,nd02,nd03,d0Lc;
  Float_t var[20];
  TString varNames[20]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","flagMC"};
  fTreeVar=new TTree("T","tree with variables");
  for(Int_t k=0;k<20;k++){
    fTreeVar->Branch(varNames[k].Data(),&var[k]);
  }
  */
  // 
  //  Rescaling of reconstructed to Lc as they were Xic added in the end of the tree
  //

  Float_t *var;

  // extra variables for c-deuteron
  if(fIsCdeuteronAnalysis){
    var = new Float_t[49];
    Short_t resp;
    TString varNames[50]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","m_pKpi","m_piKp","flagMC","w_FromLc_toXic","nSig_TPC_prot_0","nSig_TOF_prot_0","nSig_TPC_pion_0","nSig_TOF_pion_0","nSig_TPC_kaon_1","nSig_TOF_kaon_1","nSig_TPC_prot_2","nSig_TOF_prot_2","nSig_TPC_pion_2","nSig_TOF_pion_2","decayL","ndecayL","dist12","dca","pAngleXY","massHypo","cDeutMCpt","deutMCptTrk0","deutMCptTrk1","deutMCptTrk2","deutStatusTrk0","deutStatusTrk1","deutStatusTrk2","pdgMotherTrk0","pdgMotherTrk1","pdgMotherTrk2","massHypoFilt_respCuts_respPID"};
    fTreeVar=new TTree("T","tree with variables");
    for(Int_t k=0;k<49;k++){
      fTreeVar->Branch(varNames[k].Data(),&var[k]);
    }
    fTreeVar->Branch(varNames[49].Data(),&resp);
  } else{
    var = new Float_t[33];
    Short_t resp;
    TString varNames[34]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","m_pKpi","m_piKp","flagMC","w_FromLc_toXic","nSig_TPC_prot_0","nSig_TOF_prot_0","nSig_TPC_pion_0","nSig_TOF_pion_0","nSig_TPC_kaon_1","nSig_TOF_kaon_1","nSig_TPC_prot_2","nSig_TOF_prot_2","nSig_TPC_pion_2","nSig_TOF_pion_2","massHypoFilt_respCuts_respPID"};
    fTreeVar=new TTree("T","tree with variables");
    for(Int_t k=0;k<33;k++){
      fTreeVar->Branch(varNames[k].Data(),&var[k]);
    }
    fTreeVar->Branch(varNames[33].Data(),&resp);
  }
  //  fOutput->Add(fTreeVar);

  //
  //  THnSparse to study the Sc peak in MC
  //
  Int_t bin_ScPeakMC[7]        = {16,  400,  400,   500,   3,  11, 350};
  Double_t lowEdge_ScPeakMC[7] = { 0,2.250,2.440, 0.150,-0.5,-0.5,2.25};
  Double_t upEdge_ScPeakMC[7]  = {16,2.650,2.480, 0.200, 2.5,10.5,2.32};
  if(fStudyScPeakMC)  fhsparseMC_ScPeak = new THnSparseF("fhsparseMC_ScPeak","fhsparseMC_ScPeak;ptgen_Sc;recoMass_Sc;MCcalcMass_Sc;gausTerm_deltaM_Sc;charge;PIDcase;LCmass_reco;",7,bin_ScPeakMC,lowEdge_ScPeakMC,upEdge_ScPeakMC);
  //
  //  gausTerm_deltaM_Sc defined as     (massSc_reco-MCcalcMass_Sc+massSc_truefromMCparticle)-massLc_reco
  //


  // pT distribution of soft pion candidate tracks
  fPtSoftPionCand = new TH1F("fPtSoftPionCand","soft pion candidates;#it{p}_{T} (GeV/#it{c});",1000,0,0.2);
  // pT distribution of soft pion candidate tracks before SigmaC loop
  fPtSoftPionCand_insideScLoop = new TH1F("fPtSoftPionCand_insideScLoop","soft pion candidates inside SigmaC loop;#it{p}_{T} (GeV/#it{c});",1000,0,0.2);

  fOutput->Add(fDist12Signal);
  fOutput->Add(fDist12SignalFilter);
  fOutput->Add(fDist12All);
  fOutput->Add(fDist12AllFilter);
  fOutput->Add(fDist23Signal);
  fOutput->Add(fDist23All);
  fOutput->Add(fDist23AllFilter);
  fOutput->Add(fVtxResXPt);
  fOutput->Add(fVtxResYPt);
  fOutput->Add(fVtxResZPt);
  fOutput->Add(fVtxResXYPt);
  fOutput->Add(fVtxResXYZPt);
  fOutput->Add(fPrimVtxResXPt);
  fOutput->Add(fPrimVtxResYPt);
  fOutput->Add(fPrimVtxResZPt);
  fOutput->Add(fDecayLResXPt);
  fOutput->Add(fDecayLResYPt);
  fOutput->Add(fDecayLResZPt);
  fOutput->Add(fDecayLResXYPt);
  fOutput->Add(fDecayLResXYZPt);
  fOutput->Add(fCosPointDistrAll);
  fOutput->Add(fCosPointDistrAllFilter);
  fOutput->Add(fCosPointDistrSignal);
  fOutput->Add(fCosPointDistrSignalFilter);
  fOutput->Add(fhistInvMassCheck);
  fOutput->Add(fhistMCSpectrumAccLc);
  fOutput->Add(fhistMCSpectrumAccLcFromSc);
  fOutput->Add(fhistMCSpectrumAccSc);
  fOutput->Add(fhistMCSpectrumAccXic);
  fOutput->Add(fhistMCSpectrumAccCdeuteron);
  if(fhSparseAnalysis)  fOutput->Add(fhSparseAnalysis);
  if(fhSparseAnalysisSigma)  fOutput->Add(fhSparseAnalysisSigma);
  if(fReadMC && !fFillTree)  fOutput->Add(fhSparsePartReco);
  if(fReadMC && !fFillTree)  fOutput->Add(fhSparsePartGen);
  fOutput->Add(fnSigmaPIDtofProton);
  fOutput->Add(fnSigmaPIDtofPion);
  fOutput->Add(fnSigmaPIDtofKaon);
  fOutput->Add(fnSigmaPIDtpcProton);
  fOutput->Add(fnSigmaPIDtpcPion);
  fOutput->Add(fnSigmaPIDtpcKaon);
  fOutput->Add(fProtonID);
  fOutput->Add(fKaonID);
  fOutput->Add(fPionID);
  if(fStudyScPeakMC)  fOutput->Add(fhsparseMC_ScPeak);
  fOutput->Add(fPtSoftPionCand);
  fOutput->Add(fPtSoftPionCand_insideScLoop);


  // Post the data
  PostData(1,fNentries);
  PostData(2,fCounter);  
  PostData(3,fOutput);
  PostData(4,fTreeVar);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::UserExec(Option_t */*option*/)
{
  if(fDebug>=0)Printf("AliAnalysisTaskSEXicTopKpi: User Exec");
  
  /// Execute analysis for current event  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!aod) {
    printf("AliAnalysisTaskSEXicTopKpi::UserExec: input event not found!\n");
    return;
  }
  fhistMonitoring->Fill(0);
  
  if(fDebug>=0 && fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fNentries->Fill(21);
      return;
    }
    fNentries->Fill(22);
  }
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  fPidResponse=inputHandler->GetPIDResponse();
  if(!fPidResponse)return;
  
  if(fCutsXic->GetIsUsePID()){
    fCutsXic->GetPidHF()->SetPidResponse(fPidResponse);
    fCutsXic->GetPidpion()->SetPidResponse(fPidResponse);
    fCutsXic->GetPidprot()->SetPidResponse(fPidResponse);
  
    fCutsXic->GetPidHF()->SetOldPid(kFALSE);
    fCutsXic->GetPidpion()->SetOldPid(kFALSE);
    fCutsXic->GetPidprot()->SetOldPid(kFALSE);
  }
  

  //  if(!aod && AODEvent() && IsStandardAOD()) {
  //     // In case there is an AOD handler writing a standard AOD, use the AOD 
  //     // event in memory rather than the input (ESD) event.    
  //     aod = dynamic_cast<AliAODEvent*> (AODEvent());
  //     // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
  //     // have to taken from the AOD event hold by the AliAODExtension
  //     AliAODHandler* aodHandler = (AliAODHandler*) 
  //       ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());  
  //     if(aodHandler->GetExtensions()) {
  //       AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
  //       AliAODEvent* aodFromExt = ext->GetAOD();
  //       inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
  //     }
  //   } else if(aod) {
  //     inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  //   }
  
  
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001){ 
    Printf("vtx failure %p or wrong B %f",aod->GetPrimaryVertex(),aod->GetMagneticField());
    return;
  }

  AliAODMCHeader *mcHeader = 0;
  TClonesArray *lcArray = 0;
  
  if(fDebug>=0 || !fSigmaCfromLcOnTheFly){
    lcArray=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    if(!lcArray){
      Printf("Array of filtered Lc not present, delta-file not associated?");
      return;
    }
  }
  
  if(fReadMC) {
    // load MC particles    
    fmcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!fmcArray) {
      printf("AliAnalysisTaskSEXicTopKpi::UserExec: MC particles branch not found!\n");
      return;
    }
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEXicTopKpi::UserExec: MC header branch not found!\n");
      return;
    }
  }

  //
  //  Event selection in MC: require the GENERATED primary vertex to have |vtx_z|<10cm
  if(fReadMC){// EVENT SELECTION GEN STEP: CENTRALITY TO BE ADDED
    if(TMath::Abs(mcHeader->GetVtxZ())<fCutsXic->GetMaxVtxZ()){
      LoopOverGenParticles();
      fZvtx_gen_within10_MC->Fill(mcHeader->GetVtxZ()); // counter
    }
    //else return;  // apply the |vtx_z|<10cm requirement also to store reconstructed candidates
    else if(!fApplyEvSel) return;  // apply the |vtx_z|<10cm requirement also to store reconstructed candidates
  }
  
  //printf("VERTEX Z %f %f\n",vtx1->GetZ(),mcHeader->GetVtxZ());

  //------- IONUT CUT  ----------
//   Int_t nTPCout=0;
//   Float_t mTotV0=0;
//   AliAODVZERO* v0data=(AliAODVZERO*)((AliAODEvent*)aod)->GetVZEROData();
//   Float_t mTotV0A=v0data->GetMTotV0A();
//   Float_t mTotV0C=v0data->GetMTotV0C();
//   mTotV0=mTotV0A+mTotV0C;
//   Int_t ntracksEv = aod->GetNumberOfTracks();
//   for(Int_t itrack=0; itrack<ntracksEv; itrack++) { // loop on tacks
//     //    ... get the track
//     AliAODTrack * track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
//     if(!track) {AliFatal("Not a standard AOD");}
//     if(track->GetID()<0)continue;
//     if((track->GetFlags())&(AliESDtrack::kTPCout)) nTPCout++;
//     else continue;
//   }
//   if(fhMultVZEROTPCoutTrackCorrNoCut) fhMultVZEROTPCoutTrackCorrNoCut->Fill(nTPCout,mTotV0);
//   Float_t mV0Cut=-2200.+(2.5*nTPCout)+(0.000012*nTPCout*nTPCout);
//   if(fEnablePileupRejVZEROTPCout){
//     if(mTotV0<mV0Cut) return;	
//   }


  
  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fCuts,fReadMC); 
  //fCounter->StoreEvent(aod,fReadMC); 
  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);
  //   if(fReadMC && fStepMCAcc){
  //     FillMCAcceptanceHistos(fmcArray, mcHeader);
  //   }
  
  // 
  //  We need to call the IsEventSelected function
  //  using the cutobject used for the PID selections
  //  ---> it sets the Bayes PID correctly!
  //    (dirty solution, we'll clean it)
  //
  fCutsXic->IsEventSelected(aod);
  if(!fCuts->IsEventSelected(aod)) {
    if(fCuts->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
    if(fSys==1 && (fCuts->GetWhyRejection()==2 || fCuts->GetWhyRejection()==3)) fNentries->Fill(15);
    if(fCuts->GetWhyRejection()==7) fNentries->Fill(17);
    if(fApplyEvSel){     
      PostData(1,fNentries);
      PostData(2,fCounter);  
      PostData(3,fOutput);
      PostData(4,fTreeVar);
      return;
    }
  }
  fhistMonitoring->Fill(1);
  
 // Check the Nb of SDD clusters
  
  //   if (fIsRejectSDDClusters) { 
  //     Bool_t skipEvent = kFALSE;
  //     Int_t ntracks = 0;
  //     if (aod) ntracks = aod->GetNumberOfTracks();
  //     for(Int_t itrack=0; itrack<ntracks; itrack++) { // loop on tacks
  //       //    ... get the track
  //       AliAODTrack * track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
  //       if(!track) AliFatal("Not a standard AOD");
  //       if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
  // 	skipEvent=kTRUE;
  // 	fNentries->Fill(16);
  // 	break;
  //       }
  //     }
  //     if (skipEvent) return;
  //   }
  
  
  //   if(fhMultVZEROTPCoutTrackCorr)fhMultVZEROTPCoutTrackCorr->Fill(nTPCout,mTotV0);
  

  // AOD primary vertex
  fprimVtx = (AliAODVertex*)aod->GetPrimaryVertex();
  Bool_t isGoodVtx=kFALSE;
  fZvtx_reco_noSel10_MC->Fill(fprimVtx->GetZ());                // counter
  if(fReadMC) fZvtx_gen_noSel10_MC->Fill(mcHeader->GetVtxZ());  // counter
  
  TString primTitle =fprimVtx->GetTitle();
  if(primTitle.Contains("VertexerTracks") && fprimVtx->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fNentries->Fill(3);
  }
  fhistMonitoring->Fill(2);

  fvHF=new AliAnalysisVertexingHF();
 
  PrepareTracks(aod,fmcArray,mcHeader); // done always, also if only pre-filtered candidates are used: needed for SigmaC loop... or should work for a different solution
  if(fCheckOnlyTrackEfficiency){
    delete fvHF;
    PostData(1,fNentries);
    PostData(2,fCounter);  
    PostData(3,fOutput);
    PostData(4,fTreeVar);
    return;
  }

  
  //  Int_t pdgDg[3]={211,321,2212};
  
  if(fDebug>=0 || !fSigmaCfromLcOnTheFly){   // ANALYSIS WILL BE DONE ON FILTERED CANDIDATES
    LoopOverFilteredCandidates(lcArray,aod);
    if(!fSigmaCfromLcOnTheFly){
      delete fvHF;
      PostData(1,fNentries);
      PostData(2,fCounter);  
      PostData(3,fOutput);
      PostData(4,fTreeVar);
      return;
    }
  }
  
  
  /*
  // Initialize vars for TTree and set addresses if needed
  Float_t var[20];
  TString varNames[20]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","flagMC"};
  if(fFillTree){
  for(Int_t k=0;k<20;k++){
  fTreeVar->SetBranchAddress(varNames[k].Data(),&var[k]);
  }
  }*/
  // 
  //  Rescaling of reconstructed to Lc as they were Xic added in the end of the tree
  //

  Float_t *var;  
  Short_t resp;

  // c-deuteron has a few extra variables
  if(fIsCdeuteronAnalysis){
    var = new Float_t[49];
    TString varNames[50]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","m_pKpi","m_piKp","flagMC","w_FromLc_toXic","nSig_TPC_prot_0","nSig_TOF_prot_0","nSig_TPC_pion_0","nSig_TOF_pion_0","nSig_TPC_kaon_1","nSig_TOF_kaon_1","nSig_TPC_prot_2","nSig_TOF_prot_2","nSig_TPC_pion_2","nSig_TOF_pion_2","decayL","ndecayL","dist12","dca","pAngleXY","massHypo","cDeutMCpt","deutMCptTrk0","deutMCptTrk1","deutMCptTrk2","deutStatusTrk0","deutStatusTrk1","deutStatusTrk2","pdgMotherTrk0","pdgMotherTrk1","pdgMotherTrk2","massHypoFilt_respCuts_respPID"};
    if(fFillTree){
      for(Int_t k=0;k<49;k++){
        fTreeVar->SetBranchAddress(varNames[k].Data(),&var[k]);
      }
      fTreeVar->SetBranchAddress(varNames[49].Data(),&resp);
    }
  } else{
    var = new Float_t[33];
    TString varNames[34]={"pt","pAngle","lxy","nlxy","ptP","ptK","ptPi","vtxchi2","sigmaVtx","sumd02","dca1","dca2","dca3","nd01","nd02","nd03","d0Lc","cosThetaStar1","cosThetaStar2","m_pKpi","m_piKp","flagMC","w_FromLc_toXic","nSig_TPC_prot_0","nSig_TOF_prot_0","nSig_TPC_pion_0","nSig_TOF_pion_0","nSig_TPC_kaon_1","nSig_TOF_kaon_1","nSig_TPC_prot_2","nSig_TOF_prot_2","nSig_TPC_pion_2","nSig_TOF_pion_2","massHypoFilt_respCuts_respPID"};
    if(fFillTree){
      for(Int_t k=0;k<33;k++){
        fTreeVar->SetBranchAddress(varNames[k].Data(),&var[k]);
      }
      fTreeVar->SetBranchAddress(varNames[33].Data(),&resp);
    }
  } 
  
  // NOW LOOP OVER SELECTED TRACKS
  fCandCounter_onTheFly->Fill(1); // evs. entering the triple nested track loop
  for(Int_t itrack1=0;itrack1<fnSel;itrack1++){// First loop
    if(ftrackSelStatusProton->At(itrack1)<=0 && ftrackSelStatusKaon->At(itrack1)<=0 && ftrackSelStatusPion->At(itrack1)<=0)continue;//reject tracks selected only as soft-pions
    AliAODTrack *track1=(AliAODTrack*)aod->GetTrack(ftrackArraySel->At(itrack1));    
    if(!track1)continue;
    //     if(track1->Charge()<0 && !fLikeSign) continue;
    
    // HERE MAY CHECK DISPLACEMENT FOR MAKING FILTERING FASTER IN Pb-Pb
    for(Int_t itrack2=itrack1+1;itrack2<fnSel;itrack2++){// Second loop
      if(ftrackSelStatusProton->At(itrack2)<=0 && ftrackSelStatusKaon->At(itrack2)<=0 && ftrackSelStatusPion->At(itrack2)<=0)continue;//reject tracks selected only as soft-pions
      AliAODTrack *track2=(AliAODTrack*)aod->GetTrack(ftrackArraySel->At(itrack2));    
      if(!track2)continue;
      //     if(track2->Charge()>0 && !fLikeSign) continue;      
      //HERE SHOULD CHECK 2nd track and PAIR DISPLACEMENT       
      Short_t charge=track1->Charge();
      charge+=track2->Charge();
      if(charge==0){
	if(ftrackSelStatusKaon->At(itrack1)<=0 && ftrackSelStatusKaon->At(itrack2)<=0)continue; // tracks have opposite sign -> one must be compatible with being a kaon
	if( (ftrackSelStatusProton->At(itrack1)<=0 && ftrackSelStatusProton->At(itrack2)<=0) && (ftrackSelStatusPion->At(itrack1)<=0 && ftrackSelStatusPion->At(itrack2)<=0)) continue; // there should be at least one p or one pi (-> reject KK pairs)
      }
      else if(TMath::Abs(charge)==2){
	if(ftrackSelStatusProton->At(itrack1)<=0 && ftrackSelStatusProton->At(itrack2)<=0)continue;// equal-sign pair -> there must be at least one proton
	if(ftrackSelStatusPion->At(itrack1)<=0 && ftrackSelStatusPion->At(itrack2)<=0)continue;// equal-sign pair -> there must be at least one pion
      }
      else continue;// should never met this else condition since charge is checked in SelectTrack and pretended to be |charge|=1

      // THIRD LOOP
      for(Int_t itrackThird=itrack2+1;itrackThird<fnSel;itrackThird++){// Third loop
	if(itrackThird==itrack1 || itrackThird==itrack2)continue;// should never happen, line can be removed
	if(ftrackSelStatusProton->At(itrackThird)<=0 && ftrackSelStatusKaon->At(itrackThird)<=0 && ftrackSelStatusPion->At(itrackThird)<=0)continue;//needed to reject tracks selected only as soft-pions
	AliAODTrack *track3=(AliAODTrack*)aod->GetTrack(ftrackArraySel->At(itrackThird));    
	if(!track3)continue;       
	UShort_t trid[3];
	TObjArray trackArray(3);
	Int_t massHypothesis=0;// 0 = none, 1= p K pi only, 2 = pi, K, p only, 3=both
	if(charge!=0){// ++- or --+ : fill candidate anyway in the same way (opposite charge in the middle)
	  if(charge<0 && track3->Charge()<0)continue;
	  if(charge>0 && track3->Charge()>0)continue;
	  trackArray.AddAt(track1,0);
	  trackArray.AddAt(track3,1);
	  trackArray.AddAt(track2,2);
	  trid[0]=(UShort_t)track1->GetID();
	  trid[1]=(UShort_t)track3->GetID();
	  trid[2]=(UShort_t)track2->GetID();
	  // PID BASIC SELECTION HERE
	  if(ftrackSelStatusKaon->At(itrackThird)<=0)continue;// the opposite-charge particle must be compatible with the K hypo
	  if(ftrackSelStatusProton->At(itrack1)>0 && ftrackSelStatusPion->At(itrack2)>0){
	    massHypothesis+=1;// can be p K pi
	  }
	  if(ftrackSelStatusProton->At(itrack2)>0 && ftrackSelStatusPion->At(itrack1)>0){
	    massHypothesis+=2;// can be pi K p
	  }	 
	}
	else {
	  if(track3->Charge()==track1->Charge()){// +-+ or -+-
	    trackArray.AddAt(track1,0);
	    trackArray.AddAt(track2,1);
	    trackArray.AddAt(track3,2);
	    trid[0]=(UShort_t)track1->GetID();
	    trid[1]=(UShort_t)track2->GetID();
	    trid[2]=(UShort_t)track3->GetID();
	    // PID BASIC SELECTION HERE
	    if(ftrackSelStatusKaon->At(itrack2)<=0)continue;// the opposite-charge particle must be compatible with the K hypo
	    if(ftrackSelStatusProton->At(itrack1)>0 && ftrackSelStatusPion->At(itrackThird)>0){
	      massHypothesis+=1;// can be p K pi
	    }
	    if(ftrackSelStatusProton->At(itrackThird)>0 && ftrackSelStatusPion->At(itrack1)>0){
	      massHypothesis+=2;// can be pi K p
	    }	    
	  }
	  else {// -++ or +--
	    trackArray.AddAt(track2,0);
	    trackArray.AddAt(track1,1);
	    trackArray.AddAt(track3,2);
	    trid[0]=(UShort_t)track2->GetID();
	    trid[1]=(UShort_t)track1->GetID();
	    trid[2]=(UShort_t)track3->GetID();
	    // PID BASIC SELECTION HERE
	    if(ftrackSelStatusKaon->At(itrack1)<=0)continue;// the opposite-charge particle must be compatible with the K hypo
	    if(ftrackSelStatusProton->At(itrack2)>0 && ftrackSelStatusPion->At(itrackThird)>0){
	      massHypothesis+=1;// can be p K pi
	    }
	    if(ftrackSelStatusProton->At(itrackThird)>0 && ftrackSelStatusPion->At(itrack2)>0){
	      massHypothesis+=2;// can be pi K p
	    }	    
	  }
	}
	
	//      if(track3->Charge()<0 && !fLikeSign) continue;
	if(massHypothesis==0)continue; 
	fhistMonitoring->Fill(8);
	AliAODRecoDecayHF3Prong *io3Prong= new AliAODRecoDecayHF3Prong();		  
	io3Prong->SetNProngsHF(3);
	io3Prong->SetNProngs();
	io3Prong->SetProngIDs(3,trid);	
	io3Prong->SetIsFilled(0);
  fCandCounter_onTheFly->Fill(2); // on-the-fly candidate built

	if(!fvHF->FillRecoCand(aod,io3Prong)){
	  AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	  if(vtx3){delete vtx3;vtx3=0;}
	  delete io3Prong;
	  continue;
	}
  fCandCounter_onTheFly->Fill(3); // after FillRecoCand

	Bool_t unsetvtx=kFALSE;
	if(!io3Prong->GetOwnPrimaryVtx()){
	  io3Prong->SetOwnPrimaryVtx(fprimVtx);
	  unsetvtx=kTRUE;
	}
	Bool_t recPrimVtx=kFALSE;
	AliAODVertex *origownvtx=0x0;
	if(fRecalPrimVtx && fDebug<0){
	  if(io3Prong->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*io3Prong->GetOwnPrimaryVtx());
	  if(fCutsXic->RecalcOwnPrimaryVtx(io3Prong,aod))recPrimVtx=kTRUE;
	  else fCutsXic->CleanOwnPrimaryVtx(io3Prong,aod,origownvtx);
	}


	fhistMonitoring->Fill(9);
	Double_t candPt=io3Prong->Pt();
	if(!fCutsXic->IsInFiducialAcceptance(candPt,io3Prong->Y(fPdgFiducialYreco))){
	  AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	  if(vtx3){delete vtx3;vtx3=0;}
	  if(unsetvtx)io3Prong->UnsetOwnPrimaryVtx();
	  delete io3Prong;
	  continue;
	}
  fCandCounter_onTheFly->Fill(4); // in fiducial acceptance

	if(io3Prong->GetReducedChi2()>fMaxVtxChi2Cut){
	  AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	  if(vtx3){delete vtx3;vtx3=0;}
	  if(unsetvtx)io3Prong->UnsetOwnPrimaryVtx();
	  delete io3Prong;
	  continue;	  
	}
  fCandCounter_onTheFly->Fill(5); // after reducedChi2 cut

	Int_t isTrueLambdaCorXic=0,checkOrigin=-1, decay_channel=0;   // decay channel info is 0 in data
  Int_t arrayDauLabReco[3];
	AliAODMCParticle *part=0x0;
	//Double_t pointlcsc[6];
	Double_t pointlcsc[7];  // adding axis for Lc decay channel (MC)
	if(fReadMC){
	  part=MatchRecoCandtoMCAcc(io3Prong,isTrueLambdaCorXic,checkOrigin);
	  //  static Int_t CheckLcpKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);
	  //  AliVertexingHFUtils::CheckLcpKpiDecay(fmcArray,)
	  if(!part) {
	    if(fDebug>=2){
	      for(Int_t k=0;k<3;k++){
		AliAODTrack *trk = (AliAODTrack*)io3Prong->GetDaughter(k);
		Int_t dgLabel = trk->GetLabel();
		if(dgLabel>=0){
		  AliAODMCParticle *partMC=(AliAODMCParticle*)fmcArray->At(dgLabel);		
		  Int_t labMom=-1,pdgmum=-1;
		  labMom=partMC->GetMother();
		  if(labMom>=0){
		    AliAODMCParticle *partMumMC=(AliAODMCParticle*)fmcArray->At(labMom);		
		    pdgmum=partMumMC->GetPdgCode();
		  }
		  Printf("Daught %d, charge %d, label: %d, pdg: %d, mum Label: %d, pdg: %d",k,trk->Charge(),dgLabel,partMC->GetPdgCode(),labMom,pdgmum);		  
		}
		else Printf("Daught %d, chrage %d, neg label: %d",k,trk->Charge(),dgLabel);
	      }
	    }
	  }
    // if the candidate is matched, retrieve info about decay channel
    else{
      Int_t pdgcode = part->GetPdgCode();
      if(fIsXicUpgradeAnalysis) printf("===> |PDG code| = %d\n",TMath::Abs(pdgcode));
      if(TMath::Abs(pdgcode)==4122)       decay_channel = AliVertexingHFUtils::CheckLcpKpiDecay(fmcArray, part, arrayDauLabReco);
      else if(TMath::Abs(pdgcode)==4232)  decay_channel = CheckXicpKpiDecay(fmcArray, part, arrayDauLabReco);
      if((fIsCdeuteronAnalysis && (TMath::Abs(pdgcode)==2010010020)) || (fIsXicUpgradeAnalysis && TMath::Abs(pdgcode)==4232)) { // c deuteron and Xic upgrade analyses

        if((fIsXicUpgradeAnalysis && TMath::Abs(pdgcode)==4232))  printf("... Xic\n");

        // vertex information
        Double_t secVtx[3] = {0};
        Double_t secVtxMC[3] = {0};
        secVtx[0] = io3Prong->GetSecVtxX();
        secVtx[1] = io3Prong->GetSecVtxY();
        secVtx[2] = io3Prong->GetSecVtxZ();

        // get daughter for true vertex
        Int_t daughLabel = ((AliAODTrack*)io3Prong->GetDaughter(1))->GetLabel();
        if(daughLabel>=0){
          AliAODMCParticle *daugh = (AliAODMCParticle*)fmcArray->At(daughLabel);
          daugh->XvYvZv(secVtxMC);
        }
        else{
          printf("WARNING: daughLabel = %d\n",daughLabel);
        }

        fVtxResXPt->Fill(secVtx[0]-secVtxMC[0], io3Prong->Pt());
        fVtxResYPt->Fill(secVtx[1]-secVtxMC[1], io3Prong->Pt());
        fVtxResZPt->Fill(secVtx[2]-secVtxMC[2], io3Prong->Pt());
        fVtxResXYPt->Fill(TMath::Sqrt(secVtx[0]*secVtx[0] + secVtx[1]*secVtx[1]) - TMath::Sqrt(secVtxMC[0]*secVtxMC[0] + secVtxMC[1]*secVtxMC[1]), io3Prong->Pt());
        fVtxResXYZPt->Fill(TMath::Sqrt(secVtx[0]*secVtx[0] + secVtx[1]*secVtx[1] + secVtx[2]*secVtx[2]) - TMath::Sqrt(secVtxMC[0]*secVtxMC[0] + secVtxMC[1]*secVtxMC[1] + secVtxMC[2]*secVtxMC[2]), io3Prong->Pt());

        // primary vertex
        Double_t primVtx[3] = {0};
        Double_t primVtxMC[3] = {0};
        fprimVtx->GetXYZ(primVtx);

        // get origin of MC particle for true primary vertex (?)
        if(fIsXicUpgradeAnalysis) mcHeader->GetVertex(primVtxMC); // *** maybe we shall use mcHeader->GetVertex(primVtxMC); for Xic, to avoid problems with FD ***
        if(fIsCdeuteronAnalysis)  part->XvYvZv(primVtxMC);

        fPrimVtxResXPt->Fill(primVtx[0]-primVtxMC[0], io3Prong->Pt());
        fPrimVtxResYPt->Fill(primVtx[1]-primVtxMC[1], io3Prong->Pt());
        fPrimVtxResZPt->Fill(primVtx[2]-primVtxMC[2], io3Prong->Pt());

        fDecayLResXPt->Fill((secVtx[0]-primVtx[0]) - (secVtxMC[0]-primVtxMC[0]), io3Prong->Pt());
        fDecayLResYPt->Fill((secVtx[1]-primVtx[1]) - (secVtxMC[1]-primVtxMC[1]), io3Prong->Pt());
        fDecayLResZPt->Fill((secVtx[2]-primVtx[2]) - (secVtxMC[2]-primVtxMC[2]), io3Prong->Pt());
        fDecayLResXYPt->Fill(TMath::Sqrt((secVtx[0]-primVtx[0])*(secVtx[0]-primVtx[0]) +
              (secVtx[1]-primVtx[1])*(secVtx[1]-primVtx[1])) -
            TMath::Sqrt((secVtxMC[0]-primVtxMC[0])*(secVtxMC[0]-primVtxMC[0]) + 
              (secVtxMC[1]-primVtxMC[1])*(secVtxMC[1]-primVtxMC[1]))
            , io3Prong->Pt());
        fDecayLResXYZPt->Fill(TMath::Sqrt((secVtx[0]-primVtx[0])*(secVtx[0]-primVtx[0]) +
              (secVtx[1]-primVtx[1])*(secVtx[1]-primVtx[1]) +
              (secVtx[2]-primVtx[2])*(secVtx[2]-primVtx[2])) -
            TMath::Sqrt((secVtxMC[0]-primVtxMC[0])*(secVtxMC[0]-primVtxMC[0]) + 
              (secVtxMC[1]-primVtxMC[1])*(secVtxMC[1]-primVtxMC[1]) +
              (secVtxMC[2]-primVtxMC[2])*(secVtxMC[2]-primVtxMC[2]))
            , io3Prong->Pt());
      }
    }
  }

	//if(fReadMC && isTrueLambdaCorXic==1){
	// MAYBE WE'LL CAHNGE IT !!!     
	Bool_t isFromSigmaC=kFALSE;
	AliAODMCParticle *mcpartMum=0x0;
	if(fReadMC && (isTrueLambdaCorXic==10 || isTrueLambdaCorXic==20 || isTrueLambdaCorXic==40 || isTrueLambdaCorXic==50 || isTrueLambdaCorXic==80 || isTrueLambdaCorXic==100)){
	  //fhistMCSpectrumAccLc->Fill(part->Pt(),kReco,checkOrigin);
    const Double_t arr_FillkReco_Lc[4] = {part->Pt(),kReco,(Double_t)checkOrigin,(Double_t)decay_channel};
	  fhistMCSpectrumAccLc->Fill(arr_FillkReco_Lc);	  
	  // SIGMA C
	  Int_t indSc=part->GetMother();
	  if(indSc>=0){
	    mcpartMum=(AliAODMCParticle*)fmcArray->At(indSc); 
	    Int_t pdgLcMum=TMath::Abs(mcpartMum->GetPdgCode());
	    if(pdgLcMum==4112 || pdgLcMum==4222){
	      isFromSigmaC=kTRUE;
	      //fhistMCSpectrumAccSc->Fill(mcpartMum->Pt(),kRecoLc,checkOrigin);
        const Double_t arr_FillkRecoLc_Sc[4] = {mcpartMum->Pt(), kRecoLc, (Double_t)checkOrigin,(Double_t)decay_channel};	      
	      fhistMCSpectrumAccSc->Fill(arr_FillkRecoLc_Sc);	      
	      pointlcsc[0]=part->Pt();
	      pointlcsc[1]=kReco;
	      pointlcsc[2]=checkOrigin;
	      pointlcsc[3]=part->Y();
	      pointlcsc[4]=mcpartMum->Pt();
	      pointlcsc[5]=mcpartMum->Y();
        pointlcsc[6]=decay_channel;
	      fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);
	    }
	  }
	}
	
	if(fReadMC && (isTrueLambdaCorXic==30 || isTrueLambdaCorXic==60 || isTrueLambdaCorXic==120 || isTrueLambdaCorXic==240 || isTrueLambdaCorXic==150 || isTrueLambdaCorXic==300)){
    //fhistMCSpectrumAccXic->Fill(part->Pt(),kReco,checkOrigin);
    const Double_t arr_FillkReco_Xic[4] = {part->Pt(),kReco,(Double_t)checkOrigin,(Double_t)decay_channel};
	  fhistMCSpectrumAccXic->Fill(arr_FillkReco_Xic);
	}
	if(fReadMC && (isTrueLambdaCorXic==1001 || isTrueLambdaCorXic==1002)) { //c deuteron
 	  fhistMCSpectrumAccCdeuteron->Fill(part->Pt(),kReco);
	}
	if(fDebug>=0 || fCompute_dist12_dist23){
	  FillDist12and23(io3Prong,aod->GetMagneticField());
	}
	else {
	  io3Prong->SetDist12toPrim(0.05);  //needed to pass pp filtering cuts
	  io3Prong->SetDist23toPrim(0.05);	  
	}
  

	Double_t pcand[3];
	io3Prong->PxPyPz(pcand);
	AliAODTrack *trPr;
	Double_t cosThetaStarP1=-2,cosThetaStarP2=-2;
	  
	if(massHypothesis==1 || massHypothesis==3){
	  trPr=(AliAODTrack*)trackArray.At(0);
	  Double_t pprot[3];
	  trPr->PxPyPz(pprot);
	  //cosThetaStarP1=CosThetaStar(pcand,pprot,TDatabasePDG::Instance()->GetParticle(4122)->Mass(),TDatabasePDG::Instance()->GetParticle(2212)->Mass());	  
	  cosThetaStarP1=CosThetaStar(pcand,pprot,TDatabasePDG::Instance()->GetParticle(fPdgFiducialYreco)->Mass(),TDatabasePDG::Instance()->GetParticle(2212)->Mass());	  
	}
	if(massHypothesis==2 || massHypothesis==3){
	  trPr=(AliAODTrack*)trackArray.At(2);
	  Double_t pprot[3];
	  trPr->PxPyPz(pprot);
	  //cosThetaStarP2=CosThetaStar(pcand,pprot,TDatabasePDG::Instance()->GetParticle(4122)->Mass(),TDatabasePDG::Instance()->GetParticle(2212)->Mass());
	  cosThetaStarP2=CosThetaStar(pcand,pprot,TDatabasePDG::Instance()->GetParticle(fPdgFiducialYreco)->Mass(),TDatabasePDG::Instance()->GetParticle(2212)->Mass());
	}
	var[17]=cosThetaStarP1;
	var[18]=cosThetaStarP2;
	
	
	Int_t isSeleCuts=3;
	Int_t resp_onlyPID = 3;
	if(fCutsXic && (fAnalysisType==0 || fAnalysisType==2 || fAnalysisType==3)){// TO BE FIXED FOR CASE 0: SHOULD NOT DELETE, WE MUST ADD A BIT TO THE HISTOGRAM OR DUPLICATE THE HISTO FOR INCLUDING XIC
	  
	  // here cuts + PID are considered
	  isSeleCuts=fCutsXic->IsSelected(io3Prong,AliRDHFCuts::kCandidate,(AliAODEvent*)aod);
	  if((massHypothesis&isSeleCuts)>0){
	    if(fDebug>=0){
	      Printf("IsSeleCuts: %d with masshypo %d; %f , %f",isSeleCuts,massHypothesis,io3Prong->InvMassLcpKpi(),io3Prong->InvMassLcpiKp());
	      Printf("pid: 1st track: %d, %d, %d",ftrackSelStatusPion->At(itrack1),ftrackSelStatusKaon->At(itrack1),ftrackSelStatusProton->At(itrack1));
	      Printf("pid: 2nd track: %d, %d, %d",ftrackSelStatusPion->At(itrack2),ftrackSelStatusKaon->At(itrack2),ftrackSelStatusProton->At(itrack2));
	      Printf("pid: 3rd track: %d, %d, %d",ftrackSelStatusPion->At(itrackThird),ftrackSelStatusKaon->At(itrackThird),ftrackSelStatusProton->At(itrackThird));
	      PrintCandidateVariables(io3Prong,aod);
	      Printf("\n \n");
	    }
	  }

	  // store info for different selections (massHypothesis filtering, cut selection, PID selection)
	  Int_t isPIDused = fCutsXic->GetIsUsePID();
	  fCutsXic->SetUsePID(kFALSE);   // disable PID temporarly
	  Int_t resp_onlyCuts = fCutsXic->IsSelected(io3Prong,AliRDHFCuts::kCandidate,(AliAODEvent*)aod);	  
	  if(isPIDused){  // if the PID is supposed to be used, let's restore it in the cutobject
	    fCutsXic->SetUsePID(kTRUE);  // restoring PID
	    resp_onlyPID = fCutsXic->IsSelected(io3Prong,AliRDHFCuts::kPID,(AliAODEvent*)aod);
	  }
	  //
	  // fill the tree entry with a map
	  //  - first two bits: massHypothesis
	  //  - 3rd and 4th bit: resp_onlyCuts
	  //  - last 3 bits: resp_onlyPID
	  //
	  resp = SetMapCutsResponse(massHypothesis,resp_onlyCuts,resp_onlyPID);
		  
	  
	  //POSTPONED!!!
	  //
	  //massHypothesis=isSeleCuts&massHypothesis;
	  //
	  //
	  
	  //printf("massHypothesis & isSeleCuts %d\ncandPt %.2f\n==============\n",massHypothesis,candPt);
	  
	  /*
	  if(fFillTree==1){
	    FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	  }
	  else if(fFillTree==2){
	    if(isTrueLambdaCorXic){
	      if(candPt<4.){
		if(candPt*1000.-(Int_t)(candPt*1000)<0.3)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	      }
	      else FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	    }
	    else {
	      if(candPt<4.){
		if(candPt*1000.-(Int_t)(candPt*1000)<0.003)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	      }
	      else{
		if(candPt*1000.-(Int_t)(candPt*1000)<0.02)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	      }
	    }
	  }
	  else if(fFillTree==3){
	    if(isTrueLambdaCorXic){
	      if(candPt<4.){
		if(candPt*1000.-(Int_t)(candPt*1000)<0.03)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	      }
	      else FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	    }
	    else {
	      if(candPt<4.){
		if(candPt*1000.-(Int_t)(candPt*1000)<0.0003)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	      }
	      else{
		if(candPt*1000.-(Int_t)(candPt*1000)<0.002)FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	      }
	    }
	  }
	  */
	  
	  // fill the tree
	  if(fFillTree){
            if(fIsCdeuteronAnalysis){
              if(isSeleCuts) {
              // for c-deuteron: signal has no downsampling, background has downsampling
              if(fReadMC && isTrueLambdaCorXic > 1000){
                FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray,mcHeader);
	      }
              else{
                if(candPt<fpT_down){  // downsampling for low pT
	          if(candPt*1000.-(Int_t)(candPt*1000)<fLowpT_down){
	            // fill with background
		    if(fReadMC && isTrueLambdaCorXic==0)   FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray,mcHeader);
	          }
	        } 
	        else{   // downsampling for high pT
	          if(candPt*1000.-(Int_t)(candPt*1000)<fHighpT_down){
		    // fill with background
		    if(fReadMC && isTrueLambdaCorXic==0)   FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray,mcHeader);
	          }
	        }
              }
              }
            }
            else{
              if(candPt<fpT_down){  // downsampling for low pT
	        //if(candPt*1000.-(Int_t)(candPt*1000)<fLowpT_down)     FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	        if(candPt*1000.-(Int_t)(candPt*1000)<fLowpT_down){
	          // fill only with true generated particles for MC
		  if(fReadMC && isTrueLambdaCorXic>0.5)   FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray,mcHeader);
		  else if(!fReadMC)                       FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray,mcHeader);
	        }
	      } 
	      else{   // downsampling for high pT
	        //if(candPt*1000.-(Int_t)(candPt*1000)<fHighpT_down)    FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray);
	        if(candPt*1000.-(Int_t)(candPt*1000)<fHighpT_down){
		  // fill only with true generated particles for MC
		  if(fReadMC && isTrueLambdaCorXic>0.5)   FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray,mcHeader);
		  else if(!fReadMC)                       FillTree(io3Prong,massHypothesis,var,isTrueLambdaCorXic,aod,part,fmcArray,mcHeader);
	        }
	      }
            }
	  }
	  
    if(resp_onlyCuts&massHypothesis)  fCandCounter_onTheFly->Fill(6); // candidates after only cuts
    if(isSeleCuts&massHypothesis)     fCandCounter_onTheFly->Fill(7); // candidates after cuts and PID in cutobject

	  // POSTPONED HERE !!!
	  //
	  //massHypothesis=isSeleCuts&massHypothesis;
	  if(fExplore_PIDstdCuts)   massHypothesis=resp_onlyCuts&massHypothesis;  // we do not want the candidates to be filtered by Bayes PID
	  else                      massHypothesis=isSeleCuts&massHypothesis;
	  //
	  //

	  if(!massHypothesis){  // true is  massHypothesis, after PID and cuts, is null
	    if(recPrimVtx)fCutsXic->CleanOwnPrimaryVtx(io3Prong,aod,origownvtx);	    
	    if(unsetvtx)io3Prong->UnsetOwnPrimaryVtx();
	    AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	    if(vtx3){delete vtx3;  vtx3=0;}
	    delete io3Prong;
	    continue;
	  }
    fCandCounter_onTheFly->Fill(8); // candidates effectively used (beware fExplore_PIDstdCuts)
	}
      
	
	//	Int_t flagSel=FlagCandidateWithVariousCuts(io3Prong,aod,massHypothesis);
	
 	fCosPointDistrAll->Fill(io3Prong->CosPointingAngle());
 	fDist12All->Fill(io3Prong->GetDist12toPrim()*10000.);
 	fDist23All->Fill(io3Prong->GetDist23toPrim()*10000.);

	if(fReadMC){
	  if(isTrueLambdaCorXic==40 || isTrueLambdaCorXic==50 || isTrueLambdaCorXic==80 || isTrueLambdaCorXic==100 || isTrueLambdaCorXic==120 || isTrueLambdaCorXic==150 || isTrueLambdaCorXic==240 || isTrueLambdaCorXic==300 || isTrueLambdaCorXic==1001 || isTrueLambdaCorXic==1002){
	    fCosPointDistrSignal->Fill(io3Prong->CosPointingAngle());
	    fDist12Signal->Fill(io3Prong->GetDist12toPrim()*10000.);
	    //	    Printf("Dist 12, online x10000: %f",io3Prong->GetDist12toPrim()*10000.); 
	    fDist23Signal->Fill(io3Prong->GetDist23toPrim()*10000.);	    
	    //if(isTrueLambdaCorXic==1)fhistMCSpectrumAccLc->Fill(part->Pt(),kRecoCuts);
	    if(isTrueLambdaCorXic==40 || isTrueLambdaCorXic==50 || isTrueLambdaCorXic==80 || isTrueLambdaCorXic==100){
	      //fhistMCSpectrumAccLc->Fill(part->Pt(),kRecoCuts,checkOrigin);
        const Double_t arr_FillkRecoCuts_Lc[4] = {part->Pt(),kRecoCuts,(Double_t)checkOrigin,(Double_t)decay_channel};
	      fhistMCSpectrumAccLc->Fill(arr_FillkRecoCuts_Lc);
	      if(isFromSigmaC){
		//fhistMCSpectrumAccSc->Fill(mcpartMum->Pt(),kRecoLcCuts,checkOrigin);
    const Double_t arr_FillkRecoLcCuts_Sc[4] = {mcpartMum->Pt(),kRecoLcCuts,(Double_t)checkOrigin,(Double_t)decay_channel};
		fhistMCSpectrumAccSc->Fill(arr_FillkRecoLcCuts_Sc);
		pointlcsc[1]=kRecoCuts;
		fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);
	      }
	    }
	    //if(isTrueLambdaCorXic==3)fhistMCSpectrumAccXic->Fill(part->Pt(),kRecoCuts);
	    if(isTrueLambdaCorXic==120 || isTrueLambdaCorXic==150 || isTrueLambdaCorXic==240 || isTrueLambdaCorXic==300){
        //fhistMCSpectrumAccXic->Fill(part->Pt(),kRecoCuts,checkOrigin);
        const Double_t arr_FillkRecoCuts_Xic[4] = {part->Pt(),kRecoCuts,(Double_t)checkOrigin,(Double_t)decay_channel};
        fhistMCSpectrumAccXic->Fill(arr_FillkRecoCuts_Xic);
      }
      if(isTrueLambdaCorXic==1001 || isTrueLambdaCorXic==1002) fhistMCSpectrumAccCdeuteron->Fill(part->Pt(),kRecoCuts);
	  }
	}
      
	
	Double_t diffIP[3], errdiffIP[3],normIP[3],maxIP=0;
	io3Prong->Getd0MeasMinusExpProng(0,aod->GetMagneticField(),diffIP[0],errdiffIP[0]);
	io3Prong->Getd0MeasMinusExpProng(1,aod->GetMagneticField(),diffIP[1],errdiffIP[1]);
	io3Prong->Getd0MeasMinusExpProng(2,aod->GetMagneticField(),diffIP[2],errdiffIP[2]);
	normIP[0]=diffIP[0]/errdiffIP[0];
	maxIP=TMath::Abs(  normIP[0]);
	normIP[1]=diffIP[1]/errdiffIP[1];
	if(TMath::Abs(  normIP[1])>maxIP)maxIP=TMath::Abs(  normIP[1]);
	normIP[2]=diffIP[2]/errdiffIP[2];
	if(TMath::Abs(  normIP[2])>maxIP)maxIP=TMath::Abs(  normIP[2]);

	Int_t converted_isTrueLcXic = ConvertXicMCinfo(isTrueLambdaCorXic);
	//Double_t point[8]={candPt,0,io3Prong->DecayLengthXY(),io3Prong->NormalizedDecayLengthXY(),io3Prong->CosPointingAngle(),maxIP,(Double_t)converted_isTrueLcXic,0};  
	Double_t point[9]={candPt,0,io3Prong->DecayLengthXY(),io3Prong->NormalizedDecayLengthXY(),io3Prong->CosPointingAngle(),maxIP,(Double_t)converted_isTrueLcXic,0,(Double_t) decay_channel};
  
	Double_t mass1=0,mass2=0;
  int nCasesExplPID = 10;
  if(fExplPID_BayesOnlyProt)  nCasesExplPID = 11;
  const int arr_dim_PID = nCasesExplPID+1;
	Bool_t arrayPIDpkpi[arr_dim_PID],arrayPIDpikp[arr_dim_PID];
	if(massHypothesis>0 && fExplore_PIDstdCuts){
	  arrayPIDpkpi[0]=((massHypothesis&resp_onlyPID)==1 || (massHypothesis&resp_onlyPID)==3 ) ? kTRUE : kFALSE;
	  arrayPIDpikp[0]=( (massHypothesis&resp_onlyPID)==2 || (massHypothesis&resp_onlyPID)==3 ) ? kTRUE : kFALSE;
	    for(UInt_t i=1; i<=nCasesExplPID; i++){  // loop on PID cut combinations to be tested
	      arrayPIDpkpi[i]=kFALSE;
	      arrayPIDpikp[i]=kFALSE;
	      point[7] = i;
	      fCutsXic->ExplorePID(fPidResponse,io3Prong,point[7],arrayPIDpkpi[i],arrayPIDpikp[i]);	 
	    }
	}

	if(massHypothesis==1 || massHypothesis ==3) {
	  if(fIsCdeuteronAnalysis) mass1=io3Prong->InvMassCdeuterondKpi();
	  else mass1=io3Prong->InvMassLcpKpi();
	  fhistInvMassCheck->Fill(mass1,isTrueLambdaCorXic);
	  point[1]=mass1;
	  point[7]=0;
	  // this two filling fill the sparse with PID in cut object always
	  if(fhSparseAnalysis && !fExplore_PIDstdCuts){
	    if(fReadMC && (converted_isTrueLcXic==2 || converted_isTrueLcXic==6 || converted_isTrueLcXic==8 || converted_isTrueLcXic==11 || converted_isTrueLcXic==15 || converted_isTrueLcXic==17 || converted_isTrueLcXic==21 )) {
	      point[0]=part->Pt();
	      fhSparseAnalysis->Fill(point);
	    }
	    if(!fReadMC || (fReadMC&&fIsXicUpgradeAnalysis&&fIsKeepOnlyBkgXicUpgradeAnalysis) )  fhSparseAnalysis->Fill(point);
	  }
	  if(fExplore_PIDstdCuts){
	    if(fhSparseAnalysis){
	      for(UInt_t i=0; i<=nCasesExplPID; i++){  // loop on PID cut combinations to be tested
		      point[7] = i;
		      if(arrayPIDpkpi[i]){
			if(fReadMC && (converted_isTrueLcXic==2 || converted_isTrueLcXic==6 || converted_isTrueLcXic==8 || converted_isTrueLcXic==11 || converted_isTrueLcXic==15 || converted_isTrueLcXic==17)) {
			  point[0]=part->Pt();
        if(fNoStdPIDcases){
          if(i==0 || i==11) fhSparseAnalysis->Fill(point);  // avoid to fill the cases with STD PID
        }
        else  fhSparseAnalysis->Fill(point);
			}
			if(!fReadMC || (fReadMC&&fIsXicUpgradeAnalysis&&fIsKeepOnlyBkgXicUpgradeAnalysis) ){
        if(fNoStdPIDcases){
          if(i==0 || i==11) fhSparseAnalysis->Fill(point);  // avoid to fill the cases with STD PID
        }
        else  fhSparseAnalysis->Fill(point);
      }
		      }
	      }
	    }	    	   
	  }
	}
	if(massHypothesis==2 || massHypothesis ==3){
	  if(fIsCdeuteronAnalysis) mass2=io3Prong->InvMassCdeuteronpiKd();
	  else mass2=io3Prong->InvMassLcpiKp();
	  fhistInvMassCheck->Fill(mass2,isTrueLambdaCorXic);
	  point[1]=mass2;
	  point[7]=0;
	  // this two filling fill the sparse with PID in cut object always
	  if(fhSparseAnalysis && !fExplore_PIDstdCuts){
	    if(fReadMC && (converted_isTrueLcXic==3 || converted_isTrueLcXic==7 || converted_isTrueLcXic==9 || converted_isTrueLcXic==12 || converted_isTrueLcXic==16 || converted_isTrueLcXic==18 || converted_isTrueLcXic==22)) {
	      point[0]=part->Pt();
	      fhSparseAnalysis->Fill(point);
	    }
	    if(!fReadMC || (fReadMC&&fIsXicUpgradeAnalysis&&fIsKeepOnlyBkgXicUpgradeAnalysis) )  fhSparseAnalysis->Fill(point);
	  }
	  if(fExplore_PIDstdCuts){
	    if(fhSparseAnalysis){
	      for(UInt_t i=0; i<=nCasesExplPID; i++){  // loop on PID cut combinations to be tested
	        point[7] = i;
		if(arrayPIDpikp[i]){
		  if(fReadMC && (converted_isTrueLcXic==3 || converted_isTrueLcXic==7 || converted_isTrueLcXic==9 || converted_isTrueLcXic==12 || converted_isTrueLcXic==16 || converted_isTrueLcXic==18))  {
		    point[0]=part->Pt();
		    if(fNoStdPIDcases){
          if(i==0 || i==11) fhSparseAnalysis->Fill(point);  // avoid to fill the cases with STD PID
        }
        else  fhSparseAnalysis->Fill(point);
		  }
		  if(!fReadMC || (fReadMC&&fIsXicUpgradeAnalysis&&fIsKeepOnlyBkgXicUpgradeAnalysis) ){
        if(fNoStdPIDcases){
          if(i==0 || i==11) fhSparseAnalysis->Fill(point);  // avoid to fill the cases with STD PID
        }
        else  fhSparseAnalysis->Fill(point);
      }
		}
	      }
	    }	  
	  }
	}
	
	fhistMonitoring->Fill(10);
	if((fAnalysisType==0 || fAnalysisType ==3)&&fSigmaCfromLcOnTheFly){
	  if(!fReadMC){
	    if(fExplore_PIDstdCuts)SigmaCloop(io3Prong,aod,massHypothesis,mass1,mass2,point,resp_onlyPID,arrayPIDpkpi,arrayPIDpikp,itrack1,itrack2,itrackThird);
	    else SigmaCloop(io3Prong,aod,massHypothesis,mass1,mass2,point,resp_onlyPID,0x0,0x0,itrack1,itrack2,itrackThird);
	  }
	  else if(isFromSigmaC){
	    if(fExplore_PIDstdCuts)SigmaCloop(io3Prong,aod,massHypothesis,mass1,mass2,point,resp_onlyPID,arrayPIDpkpi,arrayPIDpikp,itrack1,itrack2,itrackThird,mcpartMum,(Double_t)checkOrigin,(Double_t)decay_channel);
	    else SigmaCloop(io3Prong,aod,massHypothesis,mass1,mass2,point,resp_onlyPID,0x0,0x0,itrack1,itrack2,itrackThird,mcpartMum,(Double_t)checkOrigin,(Double_t)decay_channel);
	  }
	}
	
	// NOW DELETE VTX AND CANDIDATE
	if(recPrimVtx)fCutsXic->CleanOwnPrimaryVtx(io3Prong,aod,origownvtx);	    
	if(unsetvtx)io3Prong->UnsetOwnPrimaryVtx();
	AliAODVertex *vtx3 = (AliAODVertex*)io3Prong->GetSecondaryVtx();
	if(vtx3){delete vtx3;vtx3=0;}
	delete io3Prong;
      }      
    }
  }
  delete fvHF;
  if(var) delete [] var;
  PostData(1,fNentries);
  PostData(2,fCounter);
  PostData(3,fOutput);
  PostData(4,fTreeVar);

  return;
}

//________________________________________________________
void AliAnalysisTaskSEXicTopKpi::FillArrayVariableSparse(AliAODRecoDecayHF3Prong *io3Prong,AliAODEvent *aod,Double_t *point,Int_t massHypothesis){
 
  point[0]=io3Prong->Pt();
  point[1]=0;
  point[2]=io3Prong->DecayLengthXY();
  point[3]=io3Prong->NormalizedDecayLengthXY();
  point[4]=io3Prong->CosPointingAngle();
  Double_t diffIP[3], errdiffIP[3],normIP[3],maxIP=0;
  io3Prong->Getd0MeasMinusExpProng(0,aod->GetMagneticField(),diffIP[0],errdiffIP[0]);
  io3Prong->Getd0MeasMinusExpProng(1,aod->GetMagneticField(),diffIP[1],errdiffIP[1]);
  io3Prong->Getd0MeasMinusExpProng(2,aod->GetMagneticField(),diffIP[2],errdiffIP[2]);
  normIP[0]=diffIP[0]/errdiffIP[0];
  maxIP=TMath::Abs(  normIP[0]);
  normIP[1]=diffIP[1]/errdiffIP[1];
  if(TMath::Abs(  normIP[1])>maxIP)maxIP=TMath::Abs(  normIP[1]);
  normIP[2]=diffIP[2]/errdiffIP[2];
  if(TMath::Abs(  normIP[2])>maxIP)maxIP=TMath::Abs(  normIP[2]);
  point[5]=maxIP;
  point[6]=-1;//FlagCandidateWithVariousCuts(io3Prong,aod,massHypothesis);// not needed for SigmaC! (replaced by ITSrefit for soft pion)
  point[7]=-1;
}

//________________________________________________________
void AliAnalysisTaskSEXicTopKpi::SigmaCloop(AliAODRecoDecayHF3Prong *io3Prong,AliAODEvent *aod,Int_t massHypothesis,Double_t mass1, Double_t mass2,Double_t *pointS,Int_t resp_onlyPID,Bool_t *arrayPIDselPkPi,Bool_t *arrayPIDselPikP,Int_t itrack1,Int_t itrack2,Int_t itrackThird,AliAODMCParticle *pSigmaC,Int_t checkorigin,Int_t decay_channel){
  Int_t labelSoftPi=-1;
  Double_t ptsigmacMC=-1;
  Double_t ptlambdacMC=-1;
  Double_t ysigmacMC=-9;
  Double_t ylambdacMC=-9;
  //Double_t pointlcsc[6];
  Double_t pointlcsc[7];  // adding axis for Lc decay channel (MC)
  AliAODMCParticle *mcpartLc=0x0;
  if(pSigmaC){
    ptsigmacMC=pSigmaC->Pt();
    if(pSigmaC->GetNDaughters()!=2)return;
    for(Int_t k=pSigmaC->GetDaughterLabel(0);k<=pSigmaC->GetDaughterLabel(1);k++){
      if(k>=0){
	AliAODMCParticle *mcpartScdau=(AliAODMCParticle*)fmcArray->At(k);
	if(TMath::Abs(mcpartScdau->GetPdgCode())==211){
	  labelSoftPi=k;
	}
	else if(TMath::Abs(mcpartScdau->GetPdgCode())==4122){
	  mcpartLc=mcpartScdau;
	  ptlambdacMC=mcpartLc->Pt();
	  ylambdacMC=mcpartLc->Y();
	}
      }
    }
  }    
  
  if(TMath::Abs(mass1-2.28646)>fLcMassWindowForSigmaC && TMath::Abs(mass2-2.28646)>fLcMassWindowForSigmaC)return; //Lc mass window selection  
  if(fDebug > 1){
    Printf("Good Lc candidate , will loop over %d pions",fnSelSoftPi);
  }
  //Double_t pointSigma[13];
  Double_t pointSigma[14];
  pointSigma[13] = (Double_t) decay_channel;
  pointSigma[11]=checkorigin;
  pointSigma[12]=1;
  Bool_t arrayVariableIsFilled=kFALSE;
  if(pointS){    
    for(Int_t k=0;k<8;k++){
      pointSigma[k]=pointS[k];
    }
    arrayVariableIsFilled=kTRUE;
  }
  
  // Loop over soft pions
  Double_t p2=io3Prong->P2();    
  for(Int_t isoft=0;isoft<fnSelSoftPi;isoft++){
    Int_t indsof=ftrackArraySelSoftPi->At(isoft);
    if(indsof==itrack1 || indsof==itrack2 || indsof==itrackThird)continue;
    AliAODTrack *tracksoft=(AliAODTrack*)aod->GetTrack(indsof);    		
    if(pSigmaC){
      if(TMath::Abs(tracksoft->GetLabel())!=labelSoftPi)continue;
    }

    fPtSoftPionCand_insideScLoop->Fill(tracksoft->Pt());
    
    if(itrack1==-1){// Lc from filtered candidate --> need to check ID
      Bool_t skip=kFALSE;
      for(Int_t k=0;k<3;k++){
	if((Int_t)(io3Prong->GetProngID(k))==tracksoft->GetID()){
	  //	  Printf("Skipping Lc candidate with itself");
	  skip=kTRUE;
	  break;
	}
      }
      if(skip)continue;
    }
    if(fDebug > 1){
      Printf("4plet ok, mass hypo is %d, mass1: %f, mass2: %f",massHypothesis,(massHypothesis==1||massHypothesis==3) ? mass1 : 0,(massHypothesis==2||massHypothesis==3) ? mass2:0);
    }
    
    //fhistMCSpectrumAccSc->Fill(ptsigmacMC,kReco,checkorigin);
    const Double_t arr_FillkReco_Sc[4] = {ptsigmacMC,kReco,(Double_t)checkorigin,(Double_t)decay_channel};
    fhistMCSpectrumAccSc->Fill(arr_FillkReco_Sc);
    pointlcsc[0]=ptlambdacMC;
    pointlcsc[1]=kReco;
    pointlcsc[2]=checkorigin;
    pointlcsc[3]=ylambdacMC;
    pointlcsc[4]=ptsigmacMC;
    pointlcsc[5]=ysigmacMC;
    pointlcsc[6]=decay_channel;
    fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);
    Double_t psoft[3],psoftOrig[3];
    tracksoft->PxPyPz(psoftOrig);
    psoft[0]=psoftOrig[0];
    psoft[1]=psoftOrig[1];
    psoft[2]=psoftOrig[2];
    Double_t pcand[3];
    io3Prong->PxPyPz(pcand);
    Double_t rotStep=0.;
    //printf("###### Soft pion candidate charge: %d\n",tracksoft->Charge());
    
    pointSigma[12]=1;
    if(fNRotations>1) rotStep=(fMaxAngleForRot-fMinAngleForRot)/(fNRotations-1); // -1 is to ensure that the last rotation is done with angle=fMaxAngleForRot     
    for(Int_t irot=-1; irot<fNRotations; irot++){
      // tracks are rotated to provide further background, if required
      // ASSUMPTIONS: there is no need to repeat single track selection after rotation, because we just rotate in the transverse plane (-> pt and eta does not change; the aspects related to the detector are, like potential intersection of dead modules in the roated direction are not considered)
      
      if(irot>=0){
	Double_t phirot=fMinAngleForRot+rotStep*irot;	
	psoft[0]=psoftOrig[0]*TMath::Cos(phirot)-psoftOrig[1]*TMath::Sin(phirot);
	psoft[1]=psoftOrig[0]*TMath::Sin(phirot)+psoftOrig[1]*TMath::Cos(phirot);
	pointSigma[12]=0;
      }
      
      Double_t psigma[3]={pcand[0]+psoft[0],pcand[1]+psoft[1],pcand[2]+psoft[2]};
      Double_t e1,e2;	      
      Double_t cosThetaStarSoftPi=-1.1;
      if((massHypothesis==1 || massHypothesis ==3)&&TMath::Abs(mass1-2.28646)<fLcMassWindowForSigmaC){// here we may be more restrictive and check also resp_only_pid, given that later is  done before filling the histogram
	if(fDebug > 1){
	  Printf("Ok Lc mass1 Window");
	}
	e1=TMath::Sqrt(mass1*mass1+p2);
	e2=TMath::Sqrt(0.019479785+psoft[0]*psoft[0]+psoft[1]*psoft[1]+psoft[2]*psoft[2]);// 0.019479785 =  0.13957*0.13957
	TLorentzVector lsum(psoft[0]+io3Prong->Px(),psoft[1]+io3Prong->Py(),psoft[2]+io3Prong->Pz(),e1+e2);
	//pointSigma[7]=mass1;//;
	pointSigma[8]=mass1;//;
	Double_t deltaM=lsum.M()-mass1;
	
	if(deltaM<fSigmaCDeltaMassWindow){// good candidate
	  if(fDebug > 1){
	    Printf("Ok SigmaC mass Window");
	  }	
	  // Fill array with Lc related variables if not done already
	  if(!arrayVariableIsFilled){
	    FillArrayVariableSparse(io3Prong,aod,pointSigma,massHypothesis & resp_onlyPID);// note that the last parameter here is irrelevant
	    arrayVariableIsFilled=kTRUE;	 
	  }	
	  pointSigma[6]=-1;// n.b. overwrites seleFlag (not exploited even for Lc), defined below by ITSrefit for soft pion --> this makes irrelevant which pid status is passed in the FillArrayVariableSparse above
	  // (mfaggin) now it overwrites the MC info on Lc/Xic
	  pointSigma[7]=0;
	  // now calculated remaining pair variables and fill sparse
	  if(tracksoft->TestFilterBit(AliAODTrack::kITSrefit))pointSigma[6]=0;
	  cosThetaStarSoftPi=CosThetaStar(psigma,psoft,TDatabasePDG::Instance()->GetParticle(4222)->Mass(),TDatabasePDG::Instance()->GetParticle(211)->Mass());
	  pointSigma[9]=cosThetaStarSoftPi;
	  pointSigma[1]=deltaM;	       
	  pointSigma[10]=lsum.Pt();
	  if(fhSparseAnalysisSigma && !fExplore_PIDstdCuts && (resp_onlyPID==1 || resp_onlyPID==3) )  {
	    if(!pSigmaC) fhSparseAnalysisSigma->Fill(pointSigma);
	    else {
	      AliAODTrack *trkd=(AliAODTrack*)io3Prong->GetDaughter(0);
	      AliAODMCParticle* pProt=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd->GetLabel()));
	      if(TMath::Abs(pProt->GetPdgCode())==2212){
		      pointSigma[10]=ptsigmacMC;
		      pointSigma[0]=ptlambdacMC;
		      fhSparseAnalysisSigma->Fill(pointSigma);
		      //fhistMCSpectrumAccSc->Fill(ptsigmacMC,kRecoPID,checkorigin);
          const Double_t arr_FillkRecoPID_Sc[4] = {ptsigmacMC,kRecoPID,(Double_t)checkorigin,(Double_t)decay_channel};	
          fhistMCSpectrumAccSc->Fill(arr_FillkRecoPID_Sc);      
		      pointlcsc[0]=ptlambdacMC;
		      pointlcsc[1]=kRecoPID;
		      pointlcsc[2]=checkorigin;
		      pointlcsc[3]=ylambdacMC;
		      pointlcsc[4]=ptsigmacMC;
		      pointlcsc[5]=ysigmacMC;
          pointlcsc[6]=decay_channel;
		      fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);

          //
          //  Sc peak in MC
          //
          if(fStudyScPeakMC){
            Double_t arr_ScPeakMC[7]={ptsigmacMC,-1.,-1.,-1.,(Double_t) TMath::Abs( pSigmaC->Charge()/3. ),pointSigma[7],pointSigma[8]};  // beware: AliAODMCParticle::Charge() returns the charge in unit of |e|/3 (see AliAODMCParticle--->TParticlePDG)

            // reconstructed Sc mass
            arr_ScPeakMC[1] = lsum.M();

            //////////////////////////////////////////
            // Sc mass calculated from generated pT //
            //////////////////////////////////////////
            Double_t arr_pGenP[3]      = {-1.,-1.,-1.};
            Double_t arr_pGenK[3]      = {-1.,-1.,-1.};
            Double_t arr_pGenPi[3]     = {-1.,-1.,-1.};
            Double_t arr_pGenSoftPi[3] = {-1.,-1.,-1.};
            //
            // proton (NB: here the daughter 0 is the proton)
            pProt->PxPyPz(arr_pGenP);
            Double_t m_prot  = 0.938272081;  // proton mass from PDG (GeV/c)
            Double_t E_pProt = TMath::Sqrt( m_prot*m_prot + arr_pGenP[0]*arr_pGenP[0] + arr_pGenP[1]*arr_pGenP[1] + arr_pGenP[2]*arr_pGenP[2] );
            //
            // kaon
            AliAODTrack *trkd1=(AliAODTrack*)io3Prong->GetDaughter(1);
	          AliAODMCParticle* pKaon=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd1->GetLabel()));
            if(TMath::Abs(pKaon->PdgCode())==321) pKaon->PxPyPz(arr_pGenK);
            else  printf("### WARNING: track 1 for Lc(<-Sc) NOT a kaon!\n");
            Double_t m_kaon = 0.493677; // mass kaon from PDG (GeV/c)
            Double_t E_pKaon = TMath::Sqrt( m_kaon*m_kaon + arr_pGenK[0]*arr_pGenK[0] + arr_pGenK[1]*arr_pGenK[1] + arr_pGenK[2]*arr_pGenK[2] );
            //
            // pion
            AliAODTrack *trkd2=(AliAODTrack*)io3Prong->GetDaughter(2);
	          AliAODMCParticle* pPion=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd2->GetLabel()));
            if(TMath::Abs(pPion->PdgCode())==211) pPion->PxPyPz(arr_pGenPi);
            else  printf("### WARNING: track 2 for Lc(<-Sc) NOT a pion!\n");
            Double_t m_pion = 0.13957;  // mass pion from PDG (GeV/c)
            Double_t E_pPion = TMath::Sqrt( m_pion*m_pion + arr_pGenPi[0]*arr_pGenPi[0] + arr_pGenPi[1]*arr_pGenPi[1] + arr_pGenPi[2]*arr_pGenPi[2] );
            //
            // soft pion
            AliAODMCParticle* pSoftPion=(AliAODMCParticle*)fmcArray->At(TMath::Abs(tracksoft->GetLabel()));
            if(TMath::Abs(pSoftPion->PdgCode())==211) pSoftPion->PxPyPz(arr_pGenSoftPi);
            else  printf("### WARNING: tracksoft for Sc NOT a pion!\n");
            Double_t E_pSoftPion = TMath::Sqrt( m_pion*m_pion + arr_pGenSoftPi[0]*arr_pGenSoftPi[0] + arr_pGenSoftPi[1]*arr_pGenSoftPi[1] + arr_pGenSoftPi[2]*arr_pGenSoftPi[2] );
            //
            //  create a TLorentzVector with generated p, K, pi and the soft pion
            TLorentzVector genSc_pKpi (arr_pGenP[0]+arr_pGenK[0]+arr_pGenPi[0]+arr_pGenSoftPi[0],arr_pGenP[1]+arr_pGenK[1]+arr_pGenPi[1]+arr_pGenSoftPi[1],arr_pGenP[2]+arr_pGenK[2]+arr_pGenPi[2]+arr_pGenSoftPi[2],E_pProt+E_pKaon+E_pPion+E_pSoftPion);
            //
            //  mass value computed with generated momenta
            arr_ScPeakMC[2] = genSc_pKpi.M();

            // difference between the reconstructed mass and the one computed from generated particles
            arr_ScPeakMC[3] = arr_ScPeakMC[1]-arr_ScPeakMC[2];
            // sum with the MC true mass of Sc ---> in this way, we should retrieve the gaus term of the Sc peak
            arr_ScPeakMC[3] += pSigmaC->M();
            // difference with the reconstructed Lc mass ---> in this way, we should retrieve the gaus term of the deltaM peak (since the Lc is gaussian)
            arr_ScPeakMC[3] -= pointSigma[8];

            // fill the sparse
            printf("### before filling for Sc peak: ptSc_gen=%.3f, recoMass_Sc=%.3f, calcMCmass_Sc=%3f, diff=%.6f, charge=%.3f, PIDcase=%.3f, recoMass_Lc=%.3f, PDGcode=%d\n",arr_ScPeakMC[0],arr_ScPeakMC[1],arr_ScPeakMC[2],arr_ScPeakMC[3],arr_ScPeakMC[4],arr_ScPeakMC[5],arr_ScPeakMC[6],pSigmaC->GetPdgCode());
            printf("### true mass Sc from MC: %f\n",pSigmaC->M());
            fhsparseMC_ScPeak->Fill(arr_ScPeakMC);
          }

	      }
	    }
	  }
	  if(fhSparseAnalysisSigma && fExplore_PIDstdCuts){
      int nCasesExplPID = 10;
      if(fExplPID_BayesOnlyProt)  nCasesExplPID = 11;
	    for(Int_t k=0;k<=nCasesExplPID;k++){
	      pointSigma[7]=k;
	      if(arrayPIDselPkPi[k]){
		if(!pSigmaC){
		  fhSparseAnalysisSigma->Fill(pointSigma);	    
		}
		else {
		  AliAODTrack *trkd=(AliAODTrack*)io3Prong->GetDaughter(0);
		  AliAODMCParticle* pProt=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd->GetLabel()));
		  if(TMath::Abs(pProt->GetPdgCode())==2212){
		    pointSigma[10]=ptsigmacMC;
		    pointSigma[0]=ptlambdacMC;		 
		    fhSparseAnalysisSigma->Fill(pointSigma);
		    //		  fhistMCSpectrumAccSc->Fill(ptsigmacMC,kRecoPID,checkorigin);	      
		    pointlcsc[0]=ptlambdacMC;
		    pointlcsc[1]=kRecoPID;
		    pointlcsc[2]=checkorigin;
		    pointlcsc[3]=ylambdacMC;
		    pointlcsc[4]=ptsigmacMC;
		    pointlcsc[5]=ysigmacMC;
        pointlcsc[6]=decay_channel;
		    fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);

        //
        //  Sc peak in MC
        //
        if(fStudyScPeakMC){
          Double_t arr_ScPeakMC[7]={ptsigmacMC,-1.,-1.,-1.,(Double_t) TMath::Abs( pSigmaC->Charge()/3. ),pointSigma[7],pointSigma[8]};  // beware: AliAODMCParticle::Charge() returns the charge in unit of |e|/3 (see AliAODMCParticle--->TParticlePDG)

          // reconstructed Sc mass
          arr_ScPeakMC[1] = lsum.M();

          //////////////////////////////////////////
          // Sc mass calculated from generated pT //
          //////////////////////////////////////////
          Double_t arr_pGenP[3]  = {-1.,-1.,-1.};
          Double_t arr_pGenK[3]  = {-1.,-1.,-1.};
          Double_t arr_pGenPi[3] = {-1.,-1.,-1.};
          Double_t arr_pGenSoftPi[3] = {-1.,-1.,-1.};

          //
          // proton (NB: here the daughter 0 is the proton)
          pProt->PxPyPz(arr_pGenP);
          Double_t m_prot  = 0.938272081;  // proton mass from PDG (GeV/c)
          Double_t E_pProt = TMath::Sqrt( m_prot*m_prot + arr_pGenP[0]*arr_pGenP[0] + arr_pGenP[1]*arr_pGenP[1] + arr_pGenP[2]*arr_pGenP[2] );
          //
          // kaon
          AliAODTrack *trkd1=(AliAODTrack*)io3Prong->GetDaughter(1);
	        AliAODMCParticle* pKaon=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd1->GetLabel()));
          if(TMath::Abs(pKaon->PdgCode())==321) pKaon->PxPyPz(arr_pGenK);
          else  printf("### WARNING: track 1 for Lc(<-Sc) NOT a kaon!\n");
          Double_t m_kaon = 0.493677; // mass kaon from PDG (GeV/c)
          Double_t E_pKaon = TMath::Sqrt( m_kaon*m_kaon + arr_pGenK[0]*arr_pGenK[0] + arr_pGenK[1]*arr_pGenK[1] + arr_pGenK[2]*arr_pGenK[2] );
          //
          // pion
          AliAODTrack *trkd2=(AliAODTrack*)io3Prong->GetDaughter(2);
	        AliAODMCParticle* pPion=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd2->GetLabel()));
          if(TMath::Abs(pPion->PdgCode())==211) pPion->PxPyPz(arr_pGenPi);
          else  printf("### WARNING: track 2 for Lc(<-Sc) NOT a pion!\n");
          Double_t m_pion = 0.13957;  // mass pion from PDG (GeV/c)
          Double_t E_pPion = TMath::Sqrt( m_pion*m_pion + arr_pGenPi[0]*arr_pGenPi[0] + arr_pGenPi[1]*arr_pGenPi[1] + arr_pGenPi[2]*arr_pGenPi[2] );
          //
          // soft pion
          AliAODMCParticle* pSoftPion=(AliAODMCParticle*)fmcArray->At(TMath::Abs(tracksoft->GetLabel()));
          if(TMath::Abs(pSoftPion->PdgCode())==211) pSoftPion->PxPyPz(arr_pGenSoftPi);
          else  printf("### WARNING: tracksoft for Sc NOT a pion!\n");
          Double_t E_pSoftPion = TMath::Sqrt( m_pion*m_pion + arr_pGenSoftPi[0]*arr_pGenSoftPi[0] + arr_pGenSoftPi[1]*arr_pGenSoftPi[1] + arr_pGenSoftPi[2]*arr_pGenSoftPi[2] );
          //
          //  create a TLorentzVector with generated p, K, pi and the soft pion
          TLorentzVector genSc_pKpi (arr_pGenP[0]+arr_pGenK[0]+arr_pGenPi[0]+arr_pGenSoftPi[0],arr_pGenP[1]+arr_pGenK[1]+arr_pGenPi[1]+arr_pGenSoftPi[1],arr_pGenP[2]+arr_pGenK[2]+arr_pGenPi[2]+arr_pGenSoftPi[2],E_pProt+E_pKaon+E_pPion+E_pSoftPion);
          //
          //  mass value computed with generated momenta
          arr_ScPeakMC[2] = genSc_pKpi.M();

          // difference between the reconstructed mass and the one computed from generated particles
          arr_ScPeakMC[3] = arr_ScPeakMC[1]-arr_ScPeakMC[2];
          // sum with the MC true mass of Sc ---> in this way, we should retrieve the gaus term of the Sc peak
          arr_ScPeakMC[3] += pSigmaC->M();
          // difference with the reconstructed Lc mass ---> in this way, we should retrieve the gaus term of the deltaM peak (since the Lc is gaussian)
          arr_ScPeakMC[3] -= pointSigma[8];

          // fill the sparse
          printf("### before filling for Sc peak: ptSc_gen=%.3f, recoMass_Sc=%.3f, calcMCmass_Sc=%3f, diff=%.6f, charge=%.3f, PIDcase=%.3f, recoMass_Lc=%.3f, PDGcode=%d\n",arr_ScPeakMC[0],arr_ScPeakMC[1],arr_ScPeakMC[2],arr_ScPeakMC[3],arr_ScPeakMC[4],arr_ScPeakMC[5],arr_ScPeakMC[6],pSigmaC->GetPdgCode());
          printf("### true mass Sc from MC: %f\n",pSigmaC->M());
          fhsparseMC_ScPeak->Fill(arr_ScPeakMC);
        }
		  }
		}
	      }
	    }
	  }
	}
      }
      if((massHypothesis==2 || massHypothesis ==3)&&TMath::Abs(mass2-2.28646)<fLcMassWindowForSigmaC){// here we should be more restrictive and check also resp_only_pid, given that later is done before filling the histogram
	if(fDebug > 1){
	  Printf("Ok Lc mass2 Window");
	}
	e1=TMath::Sqrt(mass2*mass2+p2);
	e2=TMath::Sqrt(0.019479785+psoft[0]*psoft[0]+psoft[1]*psoft[1]+psoft[2]*psoft[2]);// 0.019479785 =  0.13957*0.13957
	TLorentzVector lsum(psoft[0]+io3Prong->Px(),psoft[1]+io3Prong->Py(),psoft[2]+io3Prong->Pz(),e1+e2);
	pointSigma[8]=mass2;
	Double_t deltaM=lsum.M()-mass2;
	
	if(deltaM<fSigmaCDeltaMassWindow){// good candidate
	  if(fDebug > 1){
	    Printf("Ok SigmaC mass Window");
	  }
	  // Fill array with Lc related variables if not done already
	  if(!arrayVariableIsFilled){
	    FillArrayVariableSparse(io3Prong,aod,pointSigma,massHypothesis& resp_onlyPID);
	  }
	  pointSigma[6]=-1;// n.b. overwrites seleFlag (not exploited even for Lc), defined below by ITSrefit for soft pion --> this makes irrelevant which pid status is passed in the FillArrayVariableSparse above
	  pointSigma[7]=0;
	  
	  // now calculated remaining pair variables and fill sparse	
	  if(tracksoft->TestFilterBit(AliAODTrack::kITSrefit))pointSigma[6]=0;
	  cosThetaStarSoftPi=CosThetaStar(psigma,psoft,TDatabasePDG::Instance()->GetParticle(4222)->Mass(),TDatabasePDG::Instance()->GetParticle(211)->Mass());
	  pointSigma[9]=cosThetaStarSoftPi;
	  pointSigma[1]=deltaM;
	  pointSigma[10]=lsum.Pt(); // not needed
	  if(fhSparseAnalysisSigma && !fExplore_PIDstdCuts && (resp_onlyPID==2 || resp_onlyPID==3)) {
	    if(!pSigmaC)fhSparseAnalysisSigma->Fill(pointSigma);	    
	    else {
	      AliAODTrack *trkd=(AliAODTrack*)io3Prong->GetDaughter(2);
	      AliAODMCParticle* pProt=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd->GetLabel()));
	      if(TMath::Abs(pProt->GetPdgCode())==2212){
		      pointSigma[10]=ptsigmacMC;
		      pointSigma[0]=ptlambdacMC;
		      fhSparseAnalysisSigma->Fill(pointSigma);
		      //fhistMCSpectrumAccSc->Fill(ptsigmacMC,kRecoPID,checkorigin);
          const Double_t arr_FillkRecoPID_Sc[4] = {ptsigmacMC,kRecoPID,(Double_t)checkorigin,(Double_t)decay_channel};
          fhistMCSpectrumAccSc->Fill(arr_FillkRecoPID_Sc);
		      pointlcsc[0]=ptlambdacMC;
		      pointlcsc[1]=kRecoPID;
		      pointlcsc[2]=checkorigin;
		      pointlcsc[3]=ylambdacMC;
		      pointlcsc[4]=ptsigmacMC;
		      pointlcsc[5]=ysigmacMC;
          pointlcsc[6]=decay_channel;
		      fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);

          //
          //  Sc peak in MC
          //
          if(fStudyScPeakMC){
            Double_t arr_ScPeakMC[7]={ptsigmacMC,-1.,-1.,-1.,(Double_t) TMath::Abs( pSigmaC->Charge()/3. ),pointSigma[7],pointSigma[8]};  // beware: AliAODMCParticle::Charge() returns the charge in unit of |e|/3 (see AliAODMCParticle--->TParticlePDG)

            // reconstructed Sc mass
            arr_ScPeakMC[1] = lsum.M();

            //////////////////////////////////////////
            // Sc mass calculated from generated pT //
            //////////////////////////////////////////
            Double_t arr_pGenP[3]  = {-1.,-1.,-1.};
            Double_t arr_pGenK[3]  = {-1.,-1.,-1.};
            Double_t arr_pGenPi[3] = {-1.,-1.,-1.};
            Double_t arr_pGenSoftPi[3] = {-1.,-1.,-1.};
            //
            // proton (NB: here the daughter 2 is the proton)
            pProt->PxPyPz(arr_pGenP);
            Double_t m_prot  = 0.938272081;  // proton mass from PDG (GeV/c)
            Double_t E_pProt = TMath::Sqrt( m_prot*m_prot + arr_pGenP[0]*arr_pGenP[0] + arr_pGenP[1]*arr_pGenP[1] + arr_pGenP[2]*arr_pGenP[2] );
            //
            // kaon
            AliAODTrack *trkd1=(AliAODTrack*)io3Prong->GetDaughter(1);
	          AliAODMCParticle* pKaon=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd1->GetLabel()));
            if(TMath::Abs(pKaon->PdgCode())==321) pKaon->PxPyPz(arr_pGenK);
            else  printf("### WARNING: track 1 for Lc(<-Sc) NOT a kaon!\n");
            Double_t m_kaon = 0.493677; // mass kaon from PDG (GeV/c)
            Double_t E_pKaon = TMath::Sqrt( m_kaon*m_kaon + arr_pGenK[0]*arr_pGenK[0] + arr_pGenK[1]*arr_pGenK[1] + arr_pGenK[2]*arr_pGenK[2] );
            //
            // pion
            AliAODTrack *trkd0=(AliAODTrack*)io3Prong->GetDaughter(0);
	          AliAODMCParticle* pPion=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd0->GetLabel()));
            if(TMath::Abs(pPion->PdgCode())==211) pPion->PxPyPz(arr_pGenPi);
            else  printf("### WARNING: track 0 for Lc(<-Sc) NOT a pion!\n");
            Double_t m_pion = 0.13957;  // mass pion from PDG (GeV/c)
            Double_t E_pPion = TMath::Sqrt( m_pion*m_pion + arr_pGenPi[0]*arr_pGenPi[0] + arr_pGenPi[1]*arr_pGenPi[1] + arr_pGenPi[2]*arr_pGenPi[2] );
            //
            // soft pion
            AliAODMCParticle* pSoftPion=(AliAODMCParticle*)fmcArray->At(TMath::Abs(tracksoft->GetLabel()));
            if(TMath::Abs(pSoftPion->PdgCode())==211) pSoftPion->PxPyPz(arr_pGenSoftPi);
            else  printf("### WARNING: tracksoft for Sc NOT a pion!\n");
            Double_t E_pSoftPion = TMath::Sqrt( m_pion*m_pion + arr_pGenSoftPi[0]*arr_pGenSoftPi[0] + arr_pGenSoftPi[1]*arr_pGenSoftPi[1] + arr_pGenSoftPi[2]*arr_pGenSoftPi[2] );
            //
            //  create a TLorentzVector with generated p, K, pi and the soft pion
            TLorentzVector genSc_pKpi (arr_pGenP[0]+arr_pGenK[0]+arr_pGenPi[0]+arr_pGenSoftPi[0],arr_pGenP[1]+arr_pGenK[1]+arr_pGenPi[1]+arr_pGenSoftPi[1],arr_pGenP[2]+arr_pGenK[2]+arr_pGenPi[2]+arr_pGenSoftPi[2],E_pProt+E_pKaon+E_pPion+E_pSoftPion);
            //
            //  mass value computed with generated momenta
            arr_ScPeakMC[2] = genSc_pKpi.M();

            // difference between the reconstructed mass and the one computed from generated particles
            arr_ScPeakMC[3] = arr_ScPeakMC[1]-arr_ScPeakMC[2];
            // sum with the MC true mass of Sc ---> in this way, we should retrieve the gaus term of the Sc peak
            arr_ScPeakMC[3] += pSigmaC->M();
            // difference with the reconstructed Lc mass ---> in this way, we should retrieve the gaus term of the deltaM peak (since the Lc is gaussian)
            arr_ScPeakMC[3] -= pointSigma[8];

            // fill the sparse
            printf("### before filling for Sc peak: ptSc_gen=%.3f, recoMass_Sc=%.3f, calcMCmass_Sc=%3f, diff=%.6f, charge=%.3f, PIDcase=%.3f, recoMass_Lc=%.3f, PDGcode=%d\n",arr_ScPeakMC[0],arr_ScPeakMC[1],arr_ScPeakMC[2],arr_ScPeakMC[3],arr_ScPeakMC[4],arr_ScPeakMC[5],arr_ScPeakMC[6],pSigmaC->GetPdgCode());
            printf("### true mass Sc from MC: %f\n",pSigmaC->M());
            fhsparseMC_ScPeak->Fill(arr_ScPeakMC);
          }
	      }
	    }
	  }
	  if(fhSparseAnalysisSigma && fExplore_PIDstdCuts){
      int nCasesExplPID = 10;
      if(fExplPID_BayesOnlyProt)  nCasesExplPID = 11;
	    for(Int_t k=0;k<=nCasesExplPID;k++){
	      pointSigma[7]=k;
	      if(arrayPIDselPikP[k]){
		if(!pSigmaC){
		  fhSparseAnalysisSigma->Fill(pointSigma);	    
		}
		else{
		  AliAODTrack *trkd=(AliAODTrack*)io3Prong->GetDaughter(2);
		  AliAODMCParticle* pProt=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd->GetLabel()));
		  if(TMath::Abs(pProt->GetPdgCode())==2212){
		    pointSigma[10]=ptsigmacMC;
		    pointSigma[0]=ptlambdacMC;
		    fhSparseAnalysisSigma->Fill(pointSigma);
		    //		  fhistMCSpectrumAccSc->Fill(ptsigmacMC,kRecoPID,checkorigin);	      
		    pointlcsc[0]=ptlambdacMC;
		    pointlcsc[1]=kRecoPID;
		    pointlcsc[2]=checkorigin;
		    pointlcsc[3]=ylambdacMC;
		    pointlcsc[4]=ptsigmacMC;
		    pointlcsc[5]=ysigmacMC;
        pointlcsc[6]=decay_channel;
		    fhistMCSpectrumAccLcFromSc->Fill(pointlcsc);

        //
        //  Sc peak in MC
        //
        if(fStudyScPeakMC){
          Double_t arr_ScPeakMC[7]={ptsigmacMC,-1.,-1.,-1.,(Double_t) TMath::Abs( pSigmaC->Charge()/3. ),pointSigma[7],pointSigma[8]};  // beware: AliAODMCParticle::Charge() returns the charge in unit of |e|/3 (see AliAODMCParticle--->TParticlePDG)

          // reconstructed Sc mass
          arr_ScPeakMC[1] = lsum.M();

          //////////////////////////////////////////
          // Sc mass calculated from generated pT //
          //////////////////////////////////////////
          Double_t arr_pGenP[3]  = {-1.,-1.,-1.};
          Double_t arr_pGenK[3]  = {-1.,-1.,-1.};
          Double_t arr_pGenPi[3] = {-1.,-1.,-1.};
          Double_t arr_pGenSoftPi[3] = {-1.,-1.,-1.};
          //
          // proton (NB: here the daughter 2 is the proton)
          pProt->PxPyPz(arr_pGenP);
          Double_t m_prot  = 0.938272081;  // proton mass from PDG (GeV/c)
          Double_t E_pProt = TMath::Sqrt( m_prot*m_prot + arr_pGenP[0]*arr_pGenP[0] + arr_pGenP[1]*arr_pGenP[1] + arr_pGenP[2]*arr_pGenP[2] );
          //
          // kaon
          AliAODTrack *trkd1=(AliAODTrack*)io3Prong->GetDaughter(1);
	        AliAODMCParticle* pKaon=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd1->GetLabel()));
          if(TMath::Abs(pKaon->PdgCode())==321) pKaon->PxPyPz(arr_pGenK);
          else  printf("### WARNING: track 1 for Lc(<-Sc) NOT a kaon!\n");
          Double_t m_kaon = 0.493677; // mass kaon from PDG (GeV/c)
          Double_t E_pKaon = TMath::Sqrt( m_kaon*m_kaon + arr_pGenK[0]*arr_pGenK[0] + arr_pGenK[1]*arr_pGenK[1] + arr_pGenK[2]*arr_pGenK[2] );
          //
          // pion
          AliAODTrack *trkd0=(AliAODTrack*)io3Prong->GetDaughter(0);
	        AliAODMCParticle* pPion=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd0->GetLabel()));
          if(TMath::Abs(pPion->PdgCode())==211) pPion->PxPyPz(arr_pGenPi);
          else  printf("### WARNING: track 0 for Lc(<-Sc) NOT a pion!\n");
          Double_t m_pion = 0.13957;  // mass pion from PDG (GeV/c)
          Double_t E_pPion = TMath::Sqrt( m_pion*m_pion + arr_pGenPi[0]*arr_pGenPi[0] + arr_pGenPi[1]*arr_pGenPi[1] + arr_pGenPi[2]*arr_pGenPi[2] );
          //
          // soft pion
          AliAODMCParticle* pSoftPion=(AliAODMCParticle*)fmcArray->At(TMath::Abs(tracksoft->GetLabel()));
          if(TMath::Abs(pSoftPion->PdgCode())==211) pSoftPion->PxPyPz(arr_pGenSoftPi);
          else  printf("### WARNING: tracksoft for Sc NOT a pion!\n");
          Double_t E_pSoftPion = TMath::Sqrt( m_pion*m_pion + arr_pGenSoftPi[0]*arr_pGenSoftPi[0] + arr_pGenSoftPi[1]*arr_pGenSoftPi[1] + arr_pGenSoftPi[2]*arr_pGenSoftPi[2] );
          //
          //  create a TLorentzVector with generated p, K, pi and the soft pion
          TLorentzVector genSc_pKpi (arr_pGenP[0]+arr_pGenK[0]+arr_pGenPi[0]+arr_pGenSoftPi[0],arr_pGenP[1]+arr_pGenK[1]+arr_pGenPi[1]+arr_pGenSoftPi[1],arr_pGenP[2]+arr_pGenK[2]+arr_pGenPi[2]+arr_pGenSoftPi[2],E_pProt+E_pKaon+E_pPion+E_pSoftPion);
          //
          //  mass value computed with generated momenta
          arr_ScPeakMC[2] = genSc_pKpi.M();

          // difference between the reconstructed mass and the one computed from generated particles
          arr_ScPeakMC[3] = arr_ScPeakMC[1]-arr_ScPeakMC[2];
          // sum with the MC true mass of Sc ---> in this way, we should retrieve the gaus term of the Sc peak
          arr_ScPeakMC[3] += pSigmaC->M();
          // difference with the reconstructed Lc mass ---> in this way, we should retrieve the gaus term of the deltaM peak (since the Lc is gaussian)
          arr_ScPeakMC[3] -= pointSigma[8];

          // fill the sparse
          printf("### before filling for Sc peak: ptSc_gen=%.3f, recoMass_Sc=%.3f, calcMCmass_Sc=%3f, diff=%.6f, charge=%.3f, PIDcase=%.3f, recoMass_Lc=%.3f, PDGcode=%d\n",arr_ScPeakMC[0],arr_ScPeakMC[1],arr_ScPeakMC[2],arr_ScPeakMC[3],arr_ScPeakMC[4],arr_ScPeakMC[5],arr_ScPeakMC[6],pSigmaC->GetPdgCode());
          printf("### true mass Sc from MC: %f\n",pSigmaC->M());
          fhsparseMC_ScPeak->Fill(arr_ScPeakMC);
        }
		  }
		}	      
	      }
	    }
	  }
	}
      }	      
    }  
  }
}


//____________________________________________________________________________
Int_t AliAnalysisTaskSEXicTopKpi::CheckXicpKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab)const{
  /// Checks the Xic->pKpi decay channel. Returns 1 for non-resonant decays and 2 for resonant ones, -1 in other cases

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=4232) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  //Int_t labelFirstDau = mcPart->GetDaughter(0); // old
  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Int_t nProtons=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundpKpi=0;

  Int_t codeRes=-1;
  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==2212){
	nProtons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313){
	codeRes=TMath::Abs(pdgdau);
	Int_t nResDau=dau->GetNDaughters();
	if(nResDau!=2) return -1;
	//Int_t indFirstResDau=dau->GetDaughter(0); // old
	Int_t indFirstResDau=dau->GetDaughterLabel(0);
	for(Int_t resDau=0; resDau<2; resDau++){
	  Int_t indResDau=indFirstResDau+resDau;
	  if(indResDau<0) return -1;
	  AliAODMCParticle* resdau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indResDau));
	  if(!resdau) return -1;
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }else if(TMath::Abs(pdgresdau)==211){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nPions++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }
	}
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    if(nProtons!=1) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 1;
    else if(nDau==2){
      if(codeRes==313) return 2;
    }
  }
  return -1;
  
}


///--------
AliESDtrack* AliAnalysisTaskSEXicTopKpi::SelectTrack(AliAODTrack *aodtr, Int_t &isSelProton,Int_t &isSelKaon, Int_t &isSelPion,Int_t &isSelSoftPion,AliESDtrackCuts *cutsProton, AliESDtrackCuts *cutsKaon, AliESDtrackCuts *cutsPion,AliESDtrackCuts *cutsSoftPion){
  
  isSelProton=-1;
  isSelKaon=-1;
  isSelPion=-1;
  isSelSoftPion=-1;


  if(aodtr->GetID()<0)return 0x0;
  if(TMath::Abs(aodtr->Charge())!=1)return 0x0;
  Bool_t isFB4=kTRUE;
  if(!(aodtr->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA))){
    isFB4=kFALSE;
    if(fAnalysisType!=0 && fAnalysisType !=3){
      return 0x0;
    }
  }
  
  AliESDtrack *esdTrack=new AliESDtrack(aodtr);
  // set the TPC cluster info
  esdTrack->SetTPCClusterMap(aodtr->GetTPCClusterMap());
  esdTrack->SetTPCSharedMap(aodtr->GetTPCSharedMap());
  esdTrack->SetTPCPointsF(aodtr->GetTPCNclsF());
  // needed to calculate the impact parameters
  Double_t pos[3],cov[6];
  fprimVtx->GetXYZ(pos);
  fprimVtx->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  esdTrack->RelateToVertex(&vESD,0.,3.);
  if(fAnalysisType==0 || fAnalysisType ==3){
    if(cutsSoftPion){
      if(cutsSoftPion->IsSelected(esdTrack)){
	isSelSoftPion=0;
      }
    }
  }
  if(!isFB4 && isSelSoftPion<0){
    delete esdTrack;
    return 0x0;   
  }
  
  //AliESDtrackCuts *esdTrCutsAll=fCuts->GetTrackCuts();
  AliESDtrackCuts *esdTrCutsAll=fCutsXic->GetTrackCuts();
  if(!esdTrCutsAll){
    Printf("AliAnalysisTaskSEXicTopKpi :: No ESD track cuts !!!");
    delete esdTrack;
    return 0x0;
  }

  AliESDtrackCuts::ITSClusterRequirement spdreq=esdTrCutsAll->GetClusterRequirementITS(AliESDtrackCuts::kSPD);
  if(fApplykFirst){
    if(aodtr->Pt()<fMaxPtTrackkFirst){
      esdTrCutsAll->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    }
  }
  if(isFB4){
    if(esdTrCutsAll->IsSelected(esdTrack)){
      isSelProton=0;
      isSelKaon=0;
      isSelPion=0;
    }    
    else if(isSelSoftPion<0){
      if(fApplykFirst)esdTrCutsAll->SetClusterRequirementITS(AliESDtrackCuts::kSPD,spdreq);
      delete esdTrack;
      esdTrack=0x0;
    }
  }
  if(fApplykFirst)esdTrCutsAll->SetClusterRequirementITS(AliESDtrackCuts::kSPD,spdreq);

  //applying ESDtrackCut
  // if(cutsProton->IsSelected(esdTrack))isSelProton=1; 
  //   if(cutsKaon){
  //     if(cutsKaon->IsSelected(esdTrack))isSelKaon=1; 
  //   }
  //   else isSelKaon=isSelProton;
  //   if(cutsPion){
  //     if(cutsPion->IsSelected(esdTrack))isSelPion=1; 
  //   }
  //   else isSelPion=isSelProton;
  
  //   if(isSelProton==0 && isSelKaon==0 && isSelPion==0){
  //     delete esdTrack; esdTrack=0x0;
  //   }
  return esdTrack;
}



void AliAnalysisTaskSEXicTopKpi::IsSelectedPID(AliAODTrack *track,Int_t &iSelPion,Int_t &iSelKaon,Int_t &iSelProton,const Int_t iSelPionCuts,const Int_t iSelKaonCuts,const Int_t iSelProtonCuts,Bool_t fillHistos){

  iSelProton=0;
  iSelKaon=0;
  iSelPion=0;

  // TOF PID SELECTION
  AliPIDResponse::EDetPidStatus status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
  
  Double_t trpt=-1;
  if(fillHistos)trpt=track->Pt();
  
  if (status == AliPIDResponse::kDetPidOk){
    if(iSelProtonCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTOF(track,(fIsCdeuteronAnalysis?(AliPID::kDeuteron):(AliPID::kProton)));
      if(fillHistos)fnSigmaPIDtofProton->Fill(trpt,nsigma);
      //      Printf("nsigma Proton TOF: %f",nsigma);
      if(-fNSigmaPreFilterPID<=nsigma&&nsigma<=fNSigmaPreFilterPID)iSelProton++;
      else iSelProton--;
    }   
    if(iSelKaonCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
      if(fillHistos)fnSigmaPIDtofKaon->Fill(trpt,nsigma);
      if(-fNSigmaPreFilterPID<=nsigma&&nsigma<=fNSigmaPreFilterPID)iSelKaon++;
      else iSelKaon--;
    }   
    if(iSelPionCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTOF(track,AliPID::kPion);
      //	Printf("nsigma Pion TOF: %f",nsigma);
      if(fillHistos)fnSigmaPIDtofPion->Fill(trpt,nsigma);
      if(-fNSigmaPreFilterPID<=nsigma&&nsigma<=fNSigmaPreFilterPID)iSelPion++;
      else iSelPion--;
    }
  }
  else {
    if(fillHistos && trpt>0.350 && (iSelProtonCuts>=0 || iSelKaonCuts >=0 || iSelPionCuts >=0)){
      fnSigmaPIDtofPion->Fill(trpt,-30);
      fnSigmaPIDtofKaon->Fill(trpt,-30);
      fnSigmaPIDtofProton->Fill(trpt,-30);      
    }
  }
  
    //    TPC PID SELECTION
  status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,track);
  if (status == AliPIDResponse::kDetPidOk){
    if(iSelProtonCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTPC(track,(fIsCdeuteronAnalysis?(AliPID::kDeuteron):(AliPID::kProton)));
      if(fillHistos)fnSigmaPIDtpcProton->Fill(trpt,nsigma);
      //      Printf("nsigma Proton TPC: %f",nsigma);
      if(-fNSigmaPreFilterPID<=nsigma&&nsigma<=fNSigmaPreFilterPID)iSelProton++;
      else iSelProton--;
    }   
    if(iSelKaonCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
      if(fillHistos)fnSigmaPIDtpcKaon->Fill(trpt,nsigma);
      if(-fNSigmaPreFilterPID<=nsigma&&nsigma<=fNSigmaPreFilterPID)iSelKaon++;
	else iSelKaon--;
    }   
    if(iSelPionCuts>=0){
      Double_t nsigma=fPidResponse->NumberOfSigmasTPC(track,AliPID::kPion);
      if(fillHistos)fnSigmaPIDtpcPion->Fill(trpt,nsigma);
      //	Printf("nsigma Pion TPC: %f",nsigma);
      if(-fNSigmaPreFilterPID<=nsigma&&nsigma<=fNSigmaPreFilterPID)iSelPion++;
      else iSelPion--;
    }   
  }
  else {
    if(fillHistos && trpt>0.150 && (iSelProtonCuts>=0 || iSelKaonCuts >=0 || iSelPionCuts >=0)){
      fnSigmaPIDtpcPion->Fill(trpt,-30);
      fnSigmaPIDtpcKaon->Fill(trpt,-30);
      fnSigmaPIDtpcProton->Fill(trpt,-30);
    }
  }

}


//__________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::FillDist12and23(AliAODRecoDecayHF3Prong *pr,Double_t magfield){
  
  if (!fVertexerTracks) fVertexerTracks = new AliVertexerTracks(magfield);
  Double_t pos[3];
  fprimVtx->GetXYZ(pos);
//   fprimVtx->GetCovarianceMatrix(cov);
//   const AliESDVertex vESD(pos,cov,100.,100);
  
  AliESDtrack **esdTrack=new AliESDtrack*[3];
  Double_t d0z0[2],covd0z0[3];
  for(Int_t j=0;j<3;j++){
    AliAODTrack *tr=(AliAODTrack*)pr->GetDaughter(j);
    esdTrack[j]=new AliESDtrack(tr);
    // needed to calculate the impact parameters
    esdTrack[j]->PropagateToDCA(fprimVtx,magfield,kVeryBig,d0z0,covd0z0);
  }
  TObjArray *twoTrackArray=new TObjArray(2);
  twoTrackArray->AddAt(esdTrack[0],0);
  twoTrackArray->AddAt(esdTrack[1],1);

  AliESDVertex *vertexESD = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(twoTrackArray);
  //  Double_t dispersion;//,chi2perNDF;
  //  vertexESD->GetXYZ(pos); // position
  //  vertexESD->GetCovMatrix(cov); //covariance matrix
  //  chi2perNDF = vertexESD->GetChi2toNDF();
  //  dispersion = vertexESD->GetDispersion();
  Double_t dist12=TMath::Sqrt((vertexESD->GetX()-pos[0])*(vertexESD->GetX()-pos[0])+(vertexESD->GetY()-pos[1])*(vertexESD->GetY()-pos[1])+(vertexESD->GetZ()-pos[2])*(vertexESD->GetZ()-pos[2]));
  pr->SetDist12toPrim(dist12);
  delete vertexESD; vertexESD=NULL;
  
  esdTrack[1]->PropagateToDCA(fprimVtx,magfield,kVeryBig,d0z0,covd0z0);
  twoTrackArray->AddAt(esdTrack[2],0);
  twoTrackArray->AddAt(esdTrack[1],1);

  vertexESD = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(twoTrackArray);

  //  Double_t dca[3]={dcap1n1,dcap2n1,dcap1p2};
  Double_t dist23=TMath::Sqrt((vertexESD->GetX()-pos[0])*(vertexESD->GetX()-pos[0])+(vertexESD->GetY()-pos[1])*(vertexESD->GetY()-pos[1])+(vertexESD->GetZ()-pos[2])*(vertexESD->GetZ()-pos[2]));
  pr->SetDist23toPrim(dist23);
  delete vertexESD; vertexESD=NULL;
  delete twoTrackArray;

  for(Int_t j=0;j<3;j++){
    delete esdTrack[j];
  }
  delete [] esdTrack;
  return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSEXicTopKpi::FlagCandidateWithVariousCuts(AliAODRecoDecayHF3Prong *pr,AliAODEvent *aod,Int_t massHypo){
  Int_t sellevel=0;
  Double_t ptprong[3];
  ptprong[0]=pr->PtProng(0);
  ptprong[1]=pr->PtProng(1);
  ptprong[2]=pr->PtProng(2);
  
  // first set, just try harder kine selections
  if(massHypo==1||massHypo==3){
    if(ptprong[0]>1.&&ptprong[1]>0.6)sellevel|=1;
  }
  if(massHypo==2||massHypo==3){
    if(ptprong[2]>1.&&ptprong[1]>0.6)sellevel|=1;
  }


  // STRONGER PID REQUIRES DIFFERENT APPROACH
  
  // second set, improve vtx quality
  Double_t sigmaVert=pr->GetSigmaVert(aod);
  Double_t redchi2=pr->GetReducedChi2();
  Double_t maxdca2prongs=-1;
  //  Double_t maxd0prong=-1;
  Double_t sumd02=0;
  for(Int_t i=0;i<3;i++){
    Double_t d0pr=TMath::Abs(pr->Getd0Prong(i));
    sumd02+=d0pr*d0pr;
    //    if(d0pr>maxd0prong)maxd0prong=d0pr;
    Double_t dca2prongs=pr->GetDCA(i);
    if(dca2prongs>maxdca2prongs)maxdca2prongs=dca2prongs;
  }
  if(sigmaVert<0.030&&redchi2<2.&&maxdca2prongs<0.0350)sellevel|=2;
  if(sigmaVert<0.030&&sumd02>0.0001&&redchi2<1.7&&maxdca2prongs<0.0350)sellevel|=4;
  
  return sellevel;
}
//
void AliAnalysisTaskSEXicTopKpi::FillTree(AliAODRecoDecayHF3Prong *cand,Int_t massHypothesis,Float_t *varPointer, Int_t flagMC,AliAODEvent *aod, AliAODMCParticle* p=0x0, TClonesArray* array_MC=0x0, AliAODMCHeader *mcHeader=0x0){ 
  varPointer[0]=cand->Pt();
  // to save memory, let's discard all candidates with pT out of range
  if((varPointer[0]<fminpT_treeFill || varPointer[0]>fmaxpT_treeFill) && !fIsCdeuteronAnalysis)  return;
  varPointer[1]=cand->CosPointingAngle();
  varPointer[2]=cand->DecayLengthXY();
  varPointer[3]=cand->NormalizedDecayLengthXY();
    // REMOVE! put prong 0 in position 4, prong 1 in position 5, prong 2 in position 6
  /*
  UInt_t indices_prongs[3];
  if(massHypothesis<=1 || massHypothesis==3){
    varPointer[4]=cand->PtProng(0);
    varPointer[5]=cand->PtProng(1);
    varPointer[6]=cand->PtProng(2);
    indices_prongs[0]=0;
    indices_prongs[1]=1;
    indices_prongs[2]=2;
  }
  else if(massHypothesis==2){
    varPointer[4]=cand->PtProng(2);
    varPointer[5]=cand->PtProng(1);
    varPointer[6]=cand->PtProng(0);
    indices_prongs[0]=2;
    indices_prongs[1]=1;
    indices_prongs[2]=0;
  }
  */
  varPointer[4]=cand->PtProng(0);
  varPointer[5]=cand->PtProng(1);
  varPointer[6]=cand->PtProng(2);

  varPointer[7]=cand->GetReducedChi2();
  varPointer[8]=cand->GetSigmaVert(aod);
  varPointer[9]=0;
  for(Int_t i=0;i<3;i++){
    Double_t d0pr=cand->Getd0Prong(i);
    Double_t diffIP, errdiffIP;
    varPointer[10+i]=d0pr;
    varPointer[9]+=d0pr*d0pr;	    
    cand->Getd0MeasMinusExpProng(i,aod->GetMagneticField(),diffIP,errdiffIP);
    varPointer[13+i]=diffIP/errdiffIP;;
  }	     
  varPointer[16]=0.;//temporary
  varPointer[21]=flagMC;	

  // add mass (mfaggin)
    // old
  //Double_t mass_pKpi = -1.;
  //Double_t mass_piKp = -1.;
  //if(massHypothesis==3){
  //  mass_pKpi = cand->InvMassLcpKpi();
  //  mass_piKp = cand->InvMassLcpiKp();
  //}
  //else if(massHypothesis==1){ // p K pi
  //  mass_pKpi = cand->InvMassLcpKpi();
  //}
  //else if(massHypothesis==2){ // pi K p
  //  mass_piKp = cand->InvMassLcpiKp();
  //}
    // modified
  Double_t mass_pKpi, mass_piKp;
  if(fIsCdeuteronAnalysis) {
    mass_pKpi = cand->InvMassCdeuterondKpi();
    mass_piKp = cand->InvMassCdeuteronpiKd();
  }
  else{
    mass_pKpi = cand->InvMassLcpKpi();
    mass_piKp = cand->InvMassLcpiKp(); 
// to save memory, let's discard candidates with invariant mass too far from the expected one (2.467 GeV/c^2)
   //if((mass_pKpi<2.3 || mass_pKpi>2.7) && (mass_piKp<2.3 || mass_piKp>2.7)  )  return;   
   if((mass_pKpi<flowMass_treeFill || mass_pKpi>fhighMass_tree_Fill) && (mass_piKp<flowMass_treeFill || mass_piKp>fhighMass_tree_Fill)  )  return;

  }

  varPointer[19] = mass_pKpi;
  varPointer[20] = mass_piKp;

  // add info about weight to "convert" reconstructed true LC in Xic
  //if(flagMC!=1 || !p){  // generated particle associated to this reconstructed one is not a Lc or is absent
  if(flagMC<0.5 || !p){  // generated particle associated to this reconstructed one is not a Lc or is absent
    varPointer[22]=1.;
  }
  ///////////////////////////////////////////////////////////
  // flagMC==1 means that the reconstructed particle is connected to a generated Lc
  // flagMC==4(5) means that the reconstructed particle is connected to a (non-)prompt generated Lc with found quark
  // ---
  // UPDATE
  // Information about generated pKpi (*=10) or piKp (*=20), therefore the generated Lc have
  //  - flagMC = 40: prompt Lc decaying in pKpi
  //  - flagMC = 80: prompt Lc decaying in piKp
  //  - flagMC = 50: non-prompt Lc decaying in pKpi
  //  - flagMC = 100: non-prompt Lc decaying in piKp
  // ---
  ///////////////////////////////////////////////////////////
  //else if(flagMC==1 && p && array_MC){ // flagMC==1 means that the reconstructed particle is connected to a generated Lc
  else if(flagMC>0.5 /*&& flagMC<5.5*/ && p && array_MC){ 
    //Int_t index_firstProng = p->GetDaughter(0); // old
    Int_t index_firstProng = p->GetDaughterLabel(0);
    AliAODMCParticle *mc_firstProng=(AliAODMCParticle*)array_MC->At(index_firstProng);
    varPointer[22]=Weight_fromLc_toXic(p,mc_firstProng);
  }
  else  varPointer[22]=999.;

  // store the PID variables (nSigma)
  Float_t nSigma_TPC_prot_0=99, nSigma_TOF_prot_0=99, nSigma_TPC_pion_0=99, nSigma_TOF_pion_0=99;
  Float_t nSigma_TPC_kaon_1=99, nSigma_TOF_kaon_1=99; 
  Float_t nSigma_TPC_prot_2=99, nSigma_TOF_prot_2=99, nSigma_TPC_pion_2=99, nSigma_TOF_pion_2=99;
  AliPIDResponse::EDetPidStatus status_TPC_0, status_TPC_1, status_TPC_2, status_TOF_0, status_TOF_1, status_TOF_2;
  status_TPC_0 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,(AliAODTrack*)cand->GetDaughter(0));
  status_TOF_0 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,(AliAODTrack*)cand->GetDaughter(0));
  status_TPC_1 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,(AliAODTrack*)cand->GetDaughter(1));
  status_TOF_1 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,(AliAODTrack*)cand->GetDaughter(1));
  status_TPC_2 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,(AliAODTrack*)cand->GetDaughter(2));
  status_TOF_2 = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,(AliAODTrack*)cand->GetDaughter(2));
  if(status_TPC_0 == AliPIDResponse::kDetPidOk){
    nSigma_TPC_pion_0 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(0),AliPID::kPion);
    nSigma_TPC_prot_0 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(0),(fIsCdeuteronAnalysis?(AliPID::kDeuteron):(AliPID::kProton)));
  }
  if(status_TOF_0 == AliPIDResponse::kDetPidOk){
    nSigma_TOF_pion_0 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(0),AliPID::kPion);
    nSigma_TOF_prot_0 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(0),(fIsCdeuteronAnalysis?(AliPID::kDeuteron):(AliPID::kProton)));
  }
  if(status_TPC_1 == AliPIDResponse::kDetPidOk){
    nSigma_TPC_kaon_1 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(1),AliPID::kKaon);
  }
  if(status_TOF_1 == AliPIDResponse::kDetPidOk){
    nSigma_TOF_kaon_1 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(1),AliPID::kKaon);
  }
  if(status_TPC_2 == AliPIDResponse::kDetPidOk){
    nSigma_TPC_pion_2 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(2),AliPID::kPion);
    nSigma_TPC_prot_2 = fPidResponse->NumberOfSigmasTPC((AliAODTrack*)cand->GetDaughter(2),(fIsCdeuteronAnalysis?(AliPID::kDeuteron):(AliPID::kProton)));
  }
  if(status_TOF_2 == AliPIDResponse::kDetPidOk){
    nSigma_TOF_pion_2 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(2),AliPID::kPion);
    nSigma_TOF_prot_2 = fPidResponse->NumberOfSigmasTOF((AliAODTrack*)cand->GetDaughter(2),(fIsCdeuteronAnalysis?(AliPID::kDeuteron):(AliPID::kProton)));
  }
  varPointer[23] = nSigma_TPC_prot_0;
  varPointer[24] = nSigma_TOF_prot_0;
  varPointer[25] = nSigma_TPC_pion_0;
  varPointer[26] = nSigma_TOF_pion_0;
  varPointer[27] = nSigma_TPC_kaon_1;
  varPointer[28] = nSigma_TOF_kaon_1;
  varPointer[29] = nSigma_TPC_prot_2;
  varPointer[30] = nSigma_TOF_prot_2;
  varPointer[31] = nSigma_TPC_pion_2;
  varPointer[32] = nSigma_TOF_pion_2;

  // for c-deuteron fill tree with extra variables
  if(fIsCdeuteronAnalysis) {
    varPointer[33] = cand->DecayLength();
    varPointer[34] = cand->NormalizedDecayLengthXY();
    varPointer[35] = TMath::Min(cand->GetDist12toPrim(),cand->GetDist23toPrim());
    Double_t dcas[3]={0};
    cand->GetDCAs(dcas);
    varPointer[36] = TMath::Max(dcas[0],TMath::Max(dcas[1],dcas[2]));
    varPointer[37] = cand->CosPointingAngleXY();
    varPointer[38] = massHypothesis;
    varPointer[39] = -1; // MC pt c-deuteron
    varPointer[40] = -1; // MC pt background deuteron track 0
    varPointer[41] = -1; // MC pt background deuteron track 1
    varPointer[42] = -1; // MC pt background deuteron track 2
    varPointer[43] = -4; // status of deuteron - primary, secondary or other
    varPointer[44] = -4; // status of deuteron - primary, secondary or other
    varPointer[45] = -4; // status of deuteron - primary, secondary or other
    varPointer[46] = -1; // pdg of deuteron mother
    varPointer[47] = -1; // pdg of deuteron mother
    varPointer[48] = -1; // pdg of deuteron mother
    if(flagMC>1000){
      if(p && array_MC){
        Int_t mcLabel = p->GetLabel();
        if(mcLabel>=0){
          AliAODMCParticle *mcPart = (AliAODMCParticle*)array_MC->At(mcLabel);
          varPointer[39] = mcPart->Pt();
        }
      }
    }
    else{
      for(int i=0; i<3; i++){
        AliAODTrack* trk_prong = (AliAODTrack*) cand->GetDaughter(i);
        Int_t prLabel = TMath::Abs(trk_prong->GetLabel());
        if(prLabel>=0){
          AliAODMCParticle* partMC_prong = (AliAODMCParticle*) array_MC->At(prLabel);
          Int_t pdg = TMath::Abs(partMC_prong->GetPdgCode());
          if(pdg==1000010020) {
            varPointer[40+i] = partMC_prong->Pt();
            if(partMC_prong->IsPhysicalPrimary()) varPointer[43+i] = 1;
            else if(partMC_prong->IsSecondaryFromMaterial()) varPointer[43+i] = 2;
            else varPointer[43+i] = 3;
            Int_t motherLabel = partMC_prong->GetMother();
            if(motherLabel>=0) {;
              AliAODMCParticle *motherPart = (AliAODMCParticle*)array_MC->At(motherLabel);
              if(motherPart) { 
                Int_t motherPdg = TMath::Abs(motherPart->GetPdgCode());
                varPointer[46+i] = motherPdg;
              }
            }
          }
          Bool_t isTrackInjected = AliVertexingHFUtils::IsTrackInjected(trk_prong,mcHeader,array_MC);
          if(isTrackInjected) varPointer[43+i] += 3;
        }
      }
    }

  }

  fTreeVar->Fill();
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEXicTopKpi::CosThetaStar(Double_t mumVector[3],Double_t daughtVector[3],Double_t massMum,Double_t massDaught){
  
  Double_t mumP2=mumVector[0]*mumVector[0]+mumVector[1]*mumVector[1]+mumVector[2]*mumVector[2];
  Double_t mumP=TMath::Sqrt(mumP2);
  Double_t eMum=TMath::Sqrt(mumP2+massMum*massMum);
  Double_t daughtP2=daughtVector[0]*daughtVector[0]+daughtVector[1]*daughtVector[1]+daughtVector[2]*daughtVector[2];
  Double_t eDaugh=TMath::Sqrt(daughtP2+massDaught*massDaught);
  Double_t plLab=(mumVector[0]*daughtVector[0]+mumVector[1]*daughtVector[1]+mumVector[2]*daughtVector[2])/mumP;
  Double_t beta = mumP/eMum;
  Double_t gamma = eMum/massMum;
  Double_t plStar=gamma*(plLab-beta*eDaugh);
  Double_t daughtpT2=daughtP2-plLab*plLab;
  return plStar/TMath::Sqrt(plStar*plStar+daughtpT2);  
  
}

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::Terminate(Option_t */*option*/)
{
  return;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEXicTopKpi::Weight_fromLc_toXic(AliAODMCParticle* p, AliAODMCParticle* prong)
{
  // Function to calculate weight to treat reco true Lc as Xic (mfaggin)
  Double_t cTau_Lc  = 59.9E-04;  // [cTau]=cm
  Double_t cTau_Xic = 132E-04;  // [cTau]=cm  
  Double_t mass_Lc  = 2.28646;   // [m]=GeV/c^2
  Double_t mass_Xic = 2.46787;  // [m]=GeV/c^2

  TLorentzVector vecLc;
  TLorentzVector vecXic;
  vecLc.SetXYZM(p->Px(),p->Py(),p->Pz(),mass_Lc);   // [p]=GeV/c
  vecXic.SetXYZM(p->Px(),p->Py(),p->Pz(),mass_Xic);  // [p]=GeV/c

  /*
  // calc. decay length
  Double_t X_prodVtx_Lc = p->Xv();  // coordinates of vertex where the Lc is produced (primary vertex for prompt)
  Double_t Y_prodVtx_Lc = p->Yv();
  Double_t Z_prodVtx_Lc = p->Zv();
  Double_t X_prodVtx_prong = prong->Xv(); // coordinates of space point where the Lc decays
  Double_t Y_prodVtx_prong = prong->Yv();
  Double_t Z_prodVtx_prong = prong->Zv();
  Double_t decLength = TMath::Sqrt( pow(X_prodVtx_Lc-X_prodVtx_prong,2) + pow(Y_prodVtx_Lc-Y_prodVtx_prong,2) + pow(Z_prodVtx_Lc-Z_prodVtx_prong,2) );

  cout << "****** X_prodVtx_Lc: " << X_prodVtx_Lc << "   Y_prodVtx_Lc: " << Y_prodVtx_Lc << "   Z_prodVtx_Lc: " << Z_prodVtx_Lc << endl;
  cout << "****** X_prodVtx_prong: " << X_prodVtx_prong << "   Y_prodVtx_prong: " << Y_prodVtx_prong << "   Z_prodVtx_prong: " << Z_prodVtx_prong << endl; 
  Double_t ct = decLength*vecLc.E()/p->P();  // [decaylength]=cm, [m/p] = /

  // ct
  cout << "****** prong->Tv(): " << prong->Tv() << "     p->Tv(): " << p->Tv() << endl;
  */
  Double_t ct1 = 2.99792458E+10*(prong->Tv()-p->Tv());  // [c] = cm/s

  //cout << "************** ct: " << ct << "    ct1 " << ct1 << endl;

  Double_t gamma_Lc = vecLc.Gamma();
  Double_t gamma_Xic = vecXic.Gamma();

  //Double_t exp_Xic = TMath::Exp(-ct/(gamma_Xic*cTau_Xic));
  //Double_t exp_Lc  = TMath::Exp(-ct/(gamma_Lc*cTau_Lc));  
  Double_t exp_Xic = TMath::Exp(-ct1/(gamma_Xic*cTau_Xic));
  Double_t exp_Lc  = TMath::Exp(-ct1/(gamma_Lc*cTau_Lc));

  return exp_Xic/exp_Lc; 
}
//________________________________________________________________________
Short_t AliAnalysisTaskSEXicTopKpi::SetMapCutsResponse(Int_t massHypo_filtering, Int_t response_onlyCuts, Int_t response_onlyPID)
{
  //
  // Create the following map:
  //    - massHypothesis in 0 and 1 bits (00: massHypothesis = 0; 01: massHypothesis = 1; 10: massHypothesis = 2; 11: massHypothesis = 3)
  //    - resp_onlyCuts  in 2 and 3 bits (same logic, but two bits left)
  //    - resp_onlyPID   in 4 and 5 bits (same logic, but two bits left further) 
  //

  Short_t returnValue = 0;
  // massHypthesis
  if(massHypo_filtering==1)       SETBIT(returnValue,0);
  else if(massHypo_filtering==2)  SETBIT(returnValue,1);
  else if(massHypo_filtering==3){
    SETBIT(returnValue,0);
    SETBIT(returnValue,1);
  }
  // resp_onlyCuts
  if(response_onlyCuts==1)        SETBIT(returnValue,2);
  else if(response_onlyCuts==2)   SETBIT(returnValue,3);
  else if (response_onlyCuts==3)
  {
    SETBIT(returnValue,2);
    SETBIT(returnValue,3);
  }
  //resp_onlyPID
  if(response_onlyPID==1)         SETBIT(returnValue,4);
  else if(response_onlyPID==2)    SETBIT(returnValue,5);
  else if(response_onlyPID==3){
    SETBIT(returnValue,4);
    SETBIT(returnValue,5);
  }
  
  return returnValue;
}

//________________________________________________________________________
void AliAnalysisTaskSEXicTopKpi::PrepareTracks(AliAODEvent *aod,TClonesArray *mcArray, AliAODMCHeader *mcHeader){
// SELECT TRACKS and flag them: could consider to include common tracks (dca) to reduce cpu time (call once propagatetodca)
  ftrackArraySel->Reset();
  ftrackArraySelSoftPi->Reset();
  fnSel=0;
  fnSelSoftPi=0;
  
  ftrackSelStatusProton->Reset();
  ftrackSelStatusKaon->Reset();
  ftrackSelStatusPion->Reset();

  for(Int_t itrack=0;itrack<aod->GetNumberOfTracks();itrack++){
    fhistMonitoring->Fill(3);
    Int_t iSelProton=0,iSelKaon=0,iSelPion=0;
    Int_t iSelProtonCuts=-1,iSelKaonCuts=-1,iSelPionCuts=-1,iSelSoftPionCuts=-1;
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrack));
    if(!track) AliFatal("Not a standard AOD");
    if(fIsCdeuteronAnalysis && fReadMC) {
      Bool_t isTrackCD = AliVertexingHFUtils::IsTrackFromHadronDecay(2010010020, track, mcArray);
      // remove track if only keeping signal, or if it is injected and not a deuteron
      if(!isTrackCD) { 
        if(fIsKeepOnlyCdeuteronSignal) continue;
        Bool_t isBkgTrackInjected = AliVertexingHFUtils::IsTrackInjected(track,mcHeader,mcArray);
        Int_t label=TMath::Abs(track->GetLabel());
        AliAODMCParticle *mcpart=(AliAODMCParticle*)mcArray->At(label);
        Int_t pdg=TMath::Abs(mcpart->GetPdgCode());
        if(isBkgTrackInjected && pdg!=1000010020) continue; 
      }
    }
    if(fIsXicUpgradeAnalysis && fReadMC) {
      Bool_t isTrackXic = AliVertexingHFUtils::IsTrackFromHadronDecay(4232, track, mcArray);
      //std::cout << "isTrackXic: " << isTrackXic << std::endl;
      // remove track if only keeping signal, or if it is injected
      if(!isTrackXic) {
        if(fIsKeepOnlySigXicUpgradeAnalysis) continue;
        Bool_t isBkgTrackInjected = AliVertexingHFUtils::IsTrackInjected(track,mcHeader,mcArray);
        if(isBkgTrackInjected) continue;

        // if looking at bkg, keep only a small percentage of available tracks (5%)
        Double_t pt_track = track->Pt()*1000.;  // rejection from the 4th decimal digit
        if( TMath::Abs(pt_track-int(pt_track))>fRejFactorBkgUpgrade ) continue;
      }
      else{ // it means that the track comes from a true Xic
        if(fIsKeepOnlyBkgXicUpgradeAnalysis)  continue; // skip the track if we want to study pure combinatorial bkg without keeping tracks of true Xic's
      }
    }
    else if(fIsXicUpgradeAnalysis && !fReadMC && fSys==2){ // we enter here if we run on real Pb-Pb data
      Double_t pt_track = track->Pt()*1000.;  // rejection from the 4th decimal digit
      if( TMath::Abs(pt_track-int(pt_track))>fRejFactorBkgUpgrade ) continue; // if looking at bkg, keep only a small percentage of available tracks (5%)
    }

    //    Printf("selecting track");
    AliESDtrack *trackESD=SelectTrack(track,iSelProtonCuts,iSelKaonCuts,iSelPionCuts,iSelSoftPionCuts,fESDtrackCutsProton,fESDtrackCutsKaon,fESDtrackCutsPion,fESDtrackCutsSoftPion);
    
    if(!trackESD)continue;
    fhistMonitoring->Fill(4);
    //    Printf("good track");    
    
    if(iSelSoftPionCuts>=0){
      ftrackArraySelSoftPi->AddAt(itrack,fnSelSoftPi);
      fPtSoftPionCand->Fill(trackESD->Pt());
      fnSelSoftPi++;
    }    
    
    if(iSelProtonCuts < 0 && iSelKaonCuts < 0 && iSelPionCuts < 0 && iSelSoftPionCuts<0){
      delete trackESD;
      continue;
    } 
    
    ftrackArraySel->AddAt(itrack,fnSel);
  
    // PID SELECTION
    IsSelectedPID(track,iSelPion,iSelKaon,iSelProton,iSelPionCuts,iSelKaonCuts,iSelProtonCuts,kTRUE);
    //    if(itrack%50==0)Printf("Track %d, pt: %f",itrack,track->Pt());
    ftrackSelStatusProton->AddAt(iSelProton,fnSel);
    if(iSelProton>0)fhistMonitoring->Fill(7);
    
    ftrackSelStatusKaon->AddAt(iSelKaon,fnSel);
    if(iSelKaon>0)fhistMonitoring->Fill(6);
    
    ftrackSelStatusPion->AddAt(iSelPion,fnSel);
    if(iSelPion>0)fhistMonitoring->Fill(5);
    
    fnSel++;

    if(fReadMC && !fFillTree){
      Int_t label=TMath::Abs(track->GetLabel());
      AliAODMCParticle *mcpart=(AliAODMCParticle*)mcArray->At(label);
      if(mcpart->IsPhysicalPrimary()){
	Int_t pdg=TMath::Abs(mcpart->GetPdgCode());
	Int_t partType=3;
	if(pdg==211)partType=0;
	if(pdg==321)partType=1;
	if(pdg==2212)partType=2;
	Double_t point[9]={mcpart->Pt(),mcpart->Eta(),mcpart->Phi(),(Double_t)iSelPion,(Double_t)iSelSoftPionCuts,track->TestFilterBit(AliAODTrack::kITSrefit) ? (Double_t)iSelSoftPionCuts : 0,(Double_t)iSelKaon,(Double_t)iSelProton,(Double_t)partType};
	fhSparsePartReco->Fill(point);		
      }
    }

    if(fReadMC) {
      Int_t label=TMath::Abs(track->GetLabel());
      AliAODMCParticle *mcpart=(AliAODMCParticle*)mcArray->At(label);
      Int_t pdg=TMath::Abs(mcpart->GetPdgCode());
      if(iSelPion>0) {
        fPionID->Fill(1.5, mcpart->Pt());
        if(pdg==211)  fPionID->Fill(0.5, mcpart->Pt());
      }
      if(iSelKaon>0) {
        fKaonID->Fill(1.5, mcpart->Pt());
        if(pdg==321) fKaonID->Fill(0.5, mcpart->Pt());
      }
      if(iSelProton>0) {
        fProtonID->Fill(1.5, mcpart->Pt());
        if(fIsCdeuteronAnalysis?pdg==1000010020:pdg==2212) fProtonID->Fill(0.5, mcpart->Pt());
      }
    }
    
    delete trackESD;
  }
}

//______________________________________
void AliAnalysisTaskSEXicTopKpi::PrintCandidateVariables(AliAODRecoDecayHF3Prong *d,AliAODEvent *aod){
  printf("pt: %f, mass1: %f, mass 2 %f \n",d->Pt(),d->InvMassLcpKpi(),d->InvMassLcpiKp());
  printf("pt prong 0: %f, pt prong 1:%f, ptprong %f \n",d->PtProng(0),d->PtProng(1),d->PtProng(2));
  printf("d0 prong 1: %f, dist12:%f, dist23:%f \n",d->Getd0Prong(1),d->GetDist12toPrim(),d->GetDist23toPrim());
  printf("decay length: %f, sumd02 %f, CosPointAngle %f, sigmaVtx %f \n",d->DecayLength(),(d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2)),d->CosPointingAngle(),d->GetSigmaVert(aod));
  printf("dca 0: %f, dca 1: %f, dca 2: %f \n",d->GetDCA(0),d->GetDCA(1),d->GetDCA(2));
  return;
}

//______________________________________
Int_t AliAnalysisTaskSEXicTopKpi::ConvertXicMCinfo(Int_t infoMC){
  Int_t returnValue = 0.;
  
  switch (infoMC)
  {
  ////////////////////////////////////////////////////////
  //  Lc case (task version October 7th, 2019)          //
  ////////////////////////////////////////////////////////
  //    - 1:   matched Lc                     ---> 1
  //    - 10:  matched Lc, pKpi               ---> 2
  //    - 20:  matched Lc, piKp               ---> 3
  //    - 4:   matched Lc from c (prompt)     ---> 4
  //    - 5:   matched Lc from b (non-prompt) ---> 5
  //    - 40:  matched prompt Lc, pKpi        ---> 6
  //    - 80:  matched prompt Lc, piKpi       ---> 7
  //    - 50:  matched non-prompt Lc, pKpi    ---> 8
  //    - 100: matched non-prompt Lc, piKp    ---> 9
  case 1:
    returnValue = 1;
    break;
  case 10:
    returnValue = 2;
    break;
  case 20:
    returnValue = 3;
    break;
  case 4:
    returnValue = 4;
    break;
  case 5:
    returnValue = 5;
    break;
  case 40:
    returnValue = 6;
    break;
  case 80:
    returnValue = 7;
    break;
  case 50:
    returnValue = 8;
    break;
  case 100:
    returnValue = 9;
    break;

  ////////////////////////////////////////////////////////
  //  Xic case (task version October 7th, 2019)         //
  ////////////////////////////////////////////////////////
  //    - 3:   matched Xic                     ---> 10
  //    - 30:  matched Xic, pKpi               ---> 11
  //    - 60:  matched Xic, piKp               ---> 12
  //    - 12:  matched Xic from c (prompt)     ---> 13
  //    - 15:  matched Xic from b (non-prompt) ---> 14
  //    - 120: matched prompt Xic, pKpi        ---> 15
  //    - 240: matched prompt Xic, piKpi       ---> 16
  //    - 150: matched non-prompt Xic, pKpi    ---> 17
  //    - 300: matched non-prompt Xic, piKp    ---> 18
  case 3:
    returnValue = 10;
    break;
  case 30:
    returnValue = 11;
    break;
  case 60:
    returnValue = 12;
    break;
  case 12:
    returnValue = 13;
    break;
  case 15:
    returnValue = 14;
    break;
  case 120:
    returnValue = 15;
    break;
  case 240:
    returnValue = 16;
    break;
  case 150:
    returnValue = 17;
    break;
  case 300:
    returnValue = 18;
    break;
    ////////////////////////////////////////////////////////
    //  C-deuteron case (task version November 11th, 2019)//
    ////////////////////////////////////////////////////////
    //    - 1000:   matched c-deuteron, dKpi         ---> 20
    //    - 1001:   matched c-deuteron, dKpi         ---> 21
    //    - 1002:   matched c-deuteron, dKpi         ---> 22
  case 1000:
    returnValue = 20;
    break;
  case 1001:
    returnValue = 21;
    break;
  case 1002:
    returnValue = 22;
    break;
  }

  return returnValue;
}

//__________________________________________________________________________
AliAODMCParticle* AliAnalysisTaskSEXicTopKpi::MatchRecoCandtoMCAcc(AliAODRecoDecayHF3Prong *io3Prong,Int_t &isTrueLambdaCorXic,Int_t &checkOrigin){
   
  AliAODMCParticle *part=MatchRecoCandtoMC(io3Prong,isTrueLambdaCorXic,checkOrigin);
  if(!part){
    isTrueLambdaCorXic=0;
    return 0x0;
  }

  Bool_t isInAcc=kTRUE;
  // check GenAcc level
  if(fCutsXic){
    if(!fCutsXic->IsInFiducialAcceptance(part->Pt(),part->Y())){
      isInAcc=kFALSE;
      part=0x0;
      isTrueLambdaCorXic=0;
    }
  }
  else {
    if(TMath::Abs(part->Y())>0.8){
      isInAcc=kFALSE;
      part=0x0;
      isTrueLambdaCorXic=0;
    }
  }
  if(isInAcc){
    for(Int_t k=0;k<3;k++){
      AliAODTrack *trkd=(AliAODTrack*)io3Prong->GetDaughter(k);
      AliAODMCParticle *mcpartdau=(AliAODMCParticle*)fmcArray->At(TMath::Abs(trkd->GetLabel()));
      if(TMath::Abs(mcpartdau->Eta())>0.9){
	isInAcc=kFALSE;
	part=0x0;
	isTrueLambdaCorXic=0;
	break;
      }
    }
  }
  return part;
}

//__________________________________________________________________________
AliAODMCParticle* AliAnalysisTaskSEXicTopKpi::MatchRecoCandtoMC(AliAODRecoDecayHF3Prong *io3Prong,Int_t &isTrueLambdaCorXic,Int_t &checkOrigin){
  Int_t isFound=0;
  Int_t pdgDg[3]={211,321,2212};
  Int_t pdgCd[3]={211,321,1000010020};
  AliAODMCParticle *part=0x0;

  Int_t partind=io3Prong->MatchToMC(4122,fmcArray,3,pdgDg);
  if(partind>=0){
    isFound=4122;
    isTrueLambdaCorXic=1;
  }
  else {
    partind=io3Prong->MatchToMC(4232,fmcArray,3,pdgDg);
    if(partind>=0){
      isFound=4232;
      isTrueLambdaCorXic=3;
    }
    else{
      // try c deuteron
      partind=io3Prong->MatchToMC(2010010020,fmcArray,3,pdgCd);
      if(partind>=0) {
	part=(AliAODMCParticle*)fmcArray->At(partind);
	if(part){
	  // found c deuteron - label as 1000 for now
	  isTrueLambdaCorXic=1000;
	  // check daughters - label as 1001 if dKpi, 1002 if piKd
	  AliAODTrack* trk_prong = (AliAODTrack*) io3Prong->GetDaughter(0);
	  Int_t prLabel = TMath::Abs(trk_prong->GetLabel());
	  if(prLabel>0){
	    AliAODMCParticle* partMC_prong = (AliAODMCParticle*) fmcArray->At(prLabel);
	    Int_t pdg_prong = -1;
	    if(partMC_prong)  pdg_prong = TMath::Abs(partMC_prong->GetPdgCode());
	    if(pdg_prong==1000010020){      // 1st prong is a deuteron ---> dKpi
	      isTrueLambdaCorXic+=1;
	    }
	    else if(pdg_prong==211){  // 1st prong is a pion ---> piKd
	      isTrueLambdaCorXic+=2;
	    }
	  }
	}
	return part;
      }
    }  
  }

  if(partind>=0){
    part=(AliAODMCParticle*)fmcArray->At(partind);
    if(part){      
      checkOrigin= AliVertexingHFUtils::CheckOrigin(fmcArray,part,kTRUE);
      if(checkOrigin==4) isTrueLambdaCorXic*=4;      // from quark c
      else if(checkOrigin==5) isTrueLambdaCorXic*=5; // from quark b
      //
      // check if it is pKpi or piKp
      //
      AliAODTrack* trk_prong = (AliAODTrack*) io3Prong->GetDaughter(0);
      Int_t prLabel = TMath::Abs(trk_prong->GetLabel());
      if(prLabel>0){
	AliAODMCParticle* partMC_prong = (AliAODMCParticle*) fmcArray->At(prLabel);
	Int_t pdg_prong = -1;
	if(partMC_prong)  pdg_prong = TMath::Abs(partMC_prong->GetPdgCode());
	if(pdg_prong==2212){      // 1st prong is a proton ---> pKpi
	  isTrueLambdaCorXic*=10;
	}
	else if(pdg_prong==211){  // 1st prong is a pion ---> piKp
	  isTrueLambdaCorXic*=20;
	}
      }
      else{printf("---> Lc prong label %d\n",prLabel);}
    }
  }
  return part; 	  
}

//__________________________________________________________
void AliAnalysisTaskSEXicTopKpi::LoopOverGenParticles(){
  if(fDebug>=0)Printf("AliAnalysisTaskSEXicTopKpi: LoopOverGenParticless");
  // check whether lc or xic are present
  for(Int_t kmc=0;kmc<fmcArray->GetEntries();kmc++){
    AliAODMCParticle *mcpart=(AliAODMCParticle*)fmcArray->At(kmc);
    
    if(mcpart->IsPhysicalPrimary() && !fFillTree){// fill single particle histos
      Int_t pdg=TMath::Abs(mcpart->GetPdgCode());
      Int_t partType=3;
	if(pdg==211)partType=0;
	if(pdg==321)partType=1;
	if(pdg==2212)partType=2;
	Double_t point[4]={mcpart->Pt(),mcpart->Eta(),mcpart->Phi(),(Double_t)partType};
	fhSparsePartGen->Fill(point); 
      }	

      Int_t pdg=mcpart->GetPdgCode();
      Int_t arrayDauLab[3];
      if(TMath::Abs(pdg)==4122){
        Int_t decay_channel = AliVertexingHFUtils::CheckLcpKpiDecay(fmcArray, mcpart, arrayDauLab);
	if(decay_channel>=1){
	  Int_t checkOrigin=AliVertexingHFUtils::CheckOrigin(fmcArray,mcpart,kTRUE);
	  if(checkOrigin==0)continue;
	  
	  Double_t ptpart=mcpart->Pt();
	  Double_t ypart=mcpart->Y();
	  
	  // SIGMA C
    if(fDebug>=0)Printf("   AliAnalysisTaskSEXicTopKpi: LoopOverGenParticless - SigmaC stuff");
	  Bool_t isFromSigmaC=kFALSE;
	  Int_t indSc=mcpart->GetMother();
	  AliAODMCParticle *mcpartMum=0x0;
	  if(indSc>=0){
	    mcpartMum=(AliAODMCParticle*)fmcArray->At(indSc); 
	    Int_t pdgLcMum=TMath::Abs(mcpartMum->GetPdgCode());
	    if(pdgLcMum==4112 || pdgLcMum==4222)isFromSigmaC=kTRUE;
	  }
	  //Double_t pointLcSc[6];
	  Double_t pointLcSc[7];  // adding axis for Lc decay channel (MC)
	  Double_t ptpartSc;
	  Double_t ypartSc;
	  if(isFromSigmaC){
	    ptpartSc=mcpartMum->Pt();
	    ypartSc=mcpartMum->Y();
	    pointLcSc[0]=ptpart;
	    pointLcSc[1]=-1;
	    pointLcSc[2]=checkOrigin;
	    pointLcSc[3]=ypart;
	    pointLcSc[4]=ptpartSc;
	    pointLcSc[5]=ypartSc;
      pointLcSc[6]=decay_channel;
	  }
	  if(TMath::Abs(ypart)<0.5){
	    //fhistMCSpectrumAccLc->Fill(ptpart,kGenLimAcc,checkOrigin);// Gen Level
      const Double_t arr_FillkGenLimAcc_Lc[4] = {ptpart,kGenLimAcc,(Double_t)checkOrigin,(Double_t)decay_channel};
	    fhistMCSpectrumAccLc->Fill(arr_FillkGenLimAcc_Lc);// Gen Level

	    if(isFromSigmaC){
	      pointLcSc[1]=kGenLimAcc;
	      fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	    }
	  }

	  Bool_t isInAcc=kTRUE;
	  // check GenAcc level
	  if(fCutsXic){
	    if(!fCutsXic->IsInFiducialAcceptance(ptpart,ypart)){
	      isInAcc=kFALSE;
	    }
	  }
	  else {
	    if(TMath::Abs(ypart)>0.8){
	      isInAcc=kFALSE;
	    }
	  }
	  if(isInAcc){
	    //fhistMCSpectrumAccLc->Fill(ptpart,kGenAccMother,checkOrigin);// Gen Acc Mother
      const Double_t arr_FillkGenAccMother_Lc[4] = {ptpart,kGenAccMother,(Double_t)checkOrigin,(Double_t)decay_channel};
      fhistMCSpectrumAccLc->Fill(arr_FillkGenAccMother_Lc);// Gen Acc Mother
	    if(isFromSigmaC){
	      pointLcSc[1]=kGenAccMother;
	      fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	    }
	    
	    for(Int_t k=0;k<3;k++){
	      AliAODMCParticle *mcpartdau=(AliAODMCParticle*)fmcArray->At(arrayDauLab[k]);
	      if(TMath::Abs(mcpartdau->Eta())>0.9){
		isInAcc=kFALSE;
	      }	    
	    }
	    if(isInAcc){
	      //fhistMCSpectrumAccLc->Fill(mcpart->Pt(),kGenAcc,checkOrigin);// Gen Acc
        const Double_t arr_FillkGenAcc_Lc[4] = {mcpart->Pt(),kGenAcc,(Double_t)checkOrigin,(Double_t)decay_channel};
	      fhistMCSpectrumAccLc->Fill(arr_FillkGenAcc_Lc);// Gen Acc
	      if(isFromSigmaC){
		pointLcSc[1]=kGenAcc;
		fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	      }
	    }
	  }

	  if(isFromSigmaC){
	    // LimAcc level
	    if(TMath::Abs(ypartSc)<0.5){
	      //fhistMCSpectrumAccSc->Fill(ptpartSc,kGenLimAcc,checkOrigin);// Gen Level
        const Double_t arr_FillkGenLimAcc_Sc[4] = {ptpartSc,kGenLimAcc,(Double_t)checkOrigin,(Double_t)decay_channel};
	      fhistMCSpectrumAccSc->Fill(arr_FillkGenLimAcc_Sc);// Gen Level
	    }
	    // check GenAcc level
	    Bool_t isInAccSc=kTRUE;
	    if(fCutsXic){
	      if(!fCutsXic->IsInFiducialAcceptance(ptpartSc,ypartSc)){
		isInAccSc=kFALSE;
	      }
	    }
	    else {
	      if(TMath::Abs(mcpartMum->Y())>0.8){
		isInAccSc=kFALSE;
	      }
	    }
	    if(isInAccSc){
	      //fhistMCSpectrumAccSc->Fill(ptpartSc,kGenAccMother,checkOrigin);// Gen Acc Mother
        const Double_t arr_FillkGenAccMother_Sc[4] = {ptpartSc,kGenAccMother,(Double_t)checkOrigin,(Double_t)decay_channel};
	      fhistMCSpectrumAccSc->Fill(arr_FillkGenAccMother_Sc);// Gen Acc Mother
	      
	      if(isInAcc){// both Sc and Lc in fiducial acceptance + Lc daughter in Acc
		for(Int_t k=mcpartMum->GetDaughterLabel(0);k<mcpartMum->GetDaughterLabel(1);k++){
		  if(k>=0){AliAODMCParticle *mcpartMumdau=(AliAODMCParticle*)fmcArray->At(k);
		    if(TMath::Abs(mcpartMumdau->GetPdgCode()==211)&&TMath::Abs(mcpartMumdau->Eta())>0.9){
		      isInAccSc=kFALSE;
		    }	    
		  }
		}
		if(isInAccSc){
		  //fhistMCSpectrumAccSc->Fill(ptpartSc,kGenAcc,checkOrigin);// Gen Acc
      const Double_t arr_FillkGenAcc_Sc[4] = {ptpartSc,kGenAcc,(Double_t)checkOrigin,(Double_t)decay_channel};
		  fhistMCSpectrumAccSc->Fill(arr_FillkGenAcc_Sc);// Gen Acc
		}		
	      }
	    }	    
	  }
	}
      }
      else if(TMath::Abs(pdg)==4232){
        Int_t decay_channel = CheckXicpKpiDecay(fmcArray, mcpart, arrayDauLab); 
	if(decay_channel>=1){
	  Int_t checkOrigin=AliVertexingHFUtils::CheckOrigin(fmcArray,mcpart,kTRUE);
	  if(checkOrigin==0 && !fIsXicUpgradeAnalysis)continue;

	  Double_t ptpart=mcpart->Pt();
	  Double_t ypart=mcpart->Y();

	  if(TMath::Abs(ypart)<0.5){
	    //fhistMCSpectrumAccXic->Fill(ptpart,kGenLimAcc,checkOrigin);// Gen Level
      const Double_t arr_FillkGenLimAcc_Xic[4] = {ptpart,kGenLimAcc,(Double_t)checkOrigin,(Double_t)decay_channel};
	    fhistMCSpectrumAccXic->Fill(arr_FillkGenLimAcc_Xic);// Gen Level
	  }	  
	  Bool_t isInAcc=kTRUE;
	  // check GenAcc level
	  if(fCutsXic){
	    if(!fCutsXic->IsInFiducialAcceptance(ptpart,ypart)){
	      isInAcc=kFALSE;
	    }	     
	  }
	  else {
	    if(TMath::Abs(ypart)>0.8){
	      isInAcc=kFALSE;
	    }
	  }
	  if(isInAcc){
	    //fhistMCSpectrumAccXic->Fill(ptpart,kGenAccMother,checkOrigin);// Gen Acc Mother
      const Double_t arr_FillkGenAccMother_Xic[4] = {ptpart,kGenAccMother,(Double_t)checkOrigin,(Double_t)decay_channel};
	    fhistMCSpectrumAccXic->Fill(arr_FillkGenAccMother_Xic);// Gen Acc Mother
	    for(Int_t k=0;k<3;k++){
	      AliAODMCParticle *mcpartdau=(AliAODMCParticle*)fmcArray->At(arrayDauLab[k]);
	      if(TMath::Abs(mcpartdau->Eta())>0.9){
		isInAcc=kFALSE;
	      }	    
	    }
	    if(isInAcc){
	      //fhistMCSpectrumAccXic->Fill(ptpart,kGenAcc,checkOrigin);// Gen Acc
        const Double_t arr_FillkGenAcc_Xic[4] = {ptpart,kGenAcc,(Double_t)checkOrigin,(Double_t)decay_channel};
	      fhistMCSpectrumAccXic->Fill(arr_FillkGenAcc_Xic);// Gen Acc
	    }
	  }	  
	}
      }
      else if(TMath::Abs(pdg)==2010010020){
	fhistMCSpectrumAccCdeuteron->Fill(mcpart->Pt(),0);// 
	if(TMath::Abs(mcpart->Y())<0.5){
	  fhistMCSpectrumAccCdeuteron->Fill(mcpart->Pt(),kGenLimAcc);// Gen Level
	}
      }
  }
  
}

//_______________________________________
void AliAnalysisTaskSEXicTopKpi::LoopOverFilteredCandidates(TClonesArray *lcArray,AliAODEvent *aod){

  if(fDebug>=0)Printf("AliAnalysisTaskSEXicTopKpi: LoopOverFilteredCandidates");
  fCandCounter->Fill(-1);
  for(Int_t iLcFilt=0;iLcFilt<lcArray->GetEntriesFast();iLcFilt++){
    fCandCounter->Fill(0);  // loop on filtered candidates entered
    Bool_t recPrimVtx=kFALSE;
    Int_t isTrueLambdaCorXic=-1;
    Int_t checkOrigin=-1;
    Int_t decay_channel=0;  // decay channel info is 0 in data
    Int_t arrayDauLabReco[3];
    
    AliAODVertex *origownvtx=0x0;
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)lcArray->UncheckedAt(iLcFilt);      
    /*if(d->GetSelectionMap()){
      if(!d->HasSelectionBit(AliRDHFCuts::kLcCuts)&&(!d->HasSelectionBit(AliRDHFCuts::kLcPID)))continue;	
    }*/     
    fCandCounter->Fill(1);  // Lc cuts passed (filbit & QA)

    AliAODMCParticle* part=0x0;
    if( fDebug>=0 && fReadMC){
      part=MatchRecoCandtoMCAcc(d,isTrueLambdaCorXic,checkOrigin);
      // if the candidate is matched, retrieve info about decay channel
      if(part){
        Int_t pdgcode = part->GetPdgCode();
        if(TMath::Abs(pdgcode)==4122)       decay_channel =   AliVertexingHFUtils::CheckLcpKpiDecay(fmcArray, part,   arrayDauLabReco);
        else if(TMath::Abs(pdgcode)==4232)  decay_channel =   CheckXicpKpiDecay(fmcArray, part, arrayDauLabReco);
      }	  
    }
    if(!(fvHF->FillRecoCand(aod,d))) {//Fill the data members of the candidate only if they are empty.
      continue;
    }
    fCandCounter->Fill(2);  // after FillRecoCandidate

    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(fprimVtx);
      unsetvtx=kTRUE;
    }       
    if(fRecalPrimVtx && fDebug<0){// some redundancies with the above... to be fixed
      if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
      recPrimVtx=kTRUE;
      if(!fCutsXic->RecalcOwnPrimaryVtx(d,aod)){
	fCutsXic->CleanOwnPrimaryVtx(d,aod,origownvtx);
	recPrimVtx=kFALSE;
      }
    }
    Int_t iSel=3,iSelTrackCuts=3,iSelCuts=3,iSelPID=0,massHypothesis=0;
    if(!fCutsXic->IsInFiducialAcceptance(d->Pt(),d->Y(fPdgFiducialYreco)))iSel=0;
    if(!iSel){
	 if(recPrimVtx)fCutsXic->CleanOwnPrimaryVtx(d,aod,origownvtx);
	 if(unsetvtx)d->UnsetOwnPrimaryVtx();
	 continue;
    }
    fCandCounter->Fill(4);  // in fiducial acceptance

    if(d->GetReducedChi2()>fMaxVtxChi2Cut){
      if(recPrimVtx)fCutsXic->CleanOwnPrimaryVtx(d,aod,origownvtx);
	 if(unsetvtx)d->UnsetOwnPrimaryVtx();
	 continue;
    }
    Int_t iSelDebugPion[3]={-1,-1,-1},iSelDebugProton[3]={-1,-1,-1},iSelDebugKaon[3]={-1,-1,-1};
    for(Int_t itr=0;itr<3;itr++){
      AliAODTrack *track=(AliAODTrack*)d->GetDaughter(itr);
      Int_t iSelProtonCuts=-1,iSelKaonCuts=-1,iSelPionCuts=-1,iSelSoftPionCuts=-1;
      AliESDtrack *trackESD=SelectTrack(track,iSelProtonCuts,iSelKaonCuts,iSelPionCuts,iSelSoftPionCuts,fESDtrackCutsProton,fESDtrackCutsKaon,fESDtrackCutsPion,fESDtrackCutsSoftPion);// redunand because tracks were already tested and classified... one needs to build tracID[arrayIndex] maps... too heavy, could be useful only in Pb-Pb
	 if(!trackESD){
	   iSelTrackCuts=0;
	   break;
	 }
	 if(iSelProtonCuts < 0 && iSelKaonCuts < 0 && iSelPionCuts < 0 ){
	   iSelTrackCuts=0;
	   delete trackESD;
	   break;
	 }	 
	 delete trackESD;
	 
	 Int_t iSelProton=0,iSelKaon=0,iSelPion=0;//,iSelSoftPion=0;
	 IsSelectedPID(track,iSelPion,iSelKaon,iSelProton,iSelPionCuts,iSelKaonCuts,iSelProtonCuts,kFALSE);
	 iSelDebugPion[itr]=iSelPion;
	 iSelDebugProton[itr]=iSelProton;
	 iSelDebugKaon[itr]=iSelKaon;
	 if(itr==1&&iSelKaon<=0){
	   iSelPID=0;
	   break;
	 }
	 if((itr==0||itr==2)&&(iSelProton<=0&&iSelPion<=0)){	  
	   iSelPID=0;	  
	   break;
	 }
	 if(itr==0){
	   if(iSelProton>0)iSelPID++; // can be pKpi
	   if(iSelPion>0)iSelPID+=2;  // can be pi K p
	 }
	 if(itr==2){
	   if(iSelPion<=0){if(iSelPID==1||iSelPID==3)iSelPID--;}// cannot be pKpi
	   if(iSelProton<=0){if(iSelPID==2||iSelPID==3)iSelPID-=2;} //cannot be pi K p
	 }
       }

      if(iSelTrackCuts) fCandCounter->Fill(3);  // after IsSelectedTracks

       if(fDebug<0 && (iSelPID<=0 || iSelTrackCuts<=0)){
	 if(recPrimVtx)fCutsXic->CleanOwnPrimaryVtx(d,aod,origownvtx);
	 if(unsetvtx)d->UnsetOwnPrimaryVtx();	      
	 continue;
       }

       fCandCounter->Fill(5); // after IsSelectedTracks and iSelPID

       Int_t ptbin=fCutsXic->PtBin(d->Pt());
       if(ptbin!=-1)  fCandCounter->Fill(6);  // before IsSelected

       iSelCuts=fCutsXic->IsSelected(d,AliRDHFCuts::kCandidate,(AliAODEvent*)aod);// NOTE THAT PID IN CUT OBJECT IS APPLIED HERE!
       //       PrintCandidateVariables(d,aod);
       massHypothesis=iSelCuts&iSelPID;


       Bool_t isFromSigmaC=kFALSE;
       AliAODMCParticle *mcpartMum=0x0;
       if(fDebug>=0 && part){
	 
	 Double_t partpt=part->Pt(),party=part->Y();
	 //fhistMCSpectrumAccLc->Fill(partpt,kReco,checkOrigin);
   const Double_t arr_FillkReco_Lc[4] = {partpt,kReco,(Double_t)checkOrigin,(Double_t)decay_channel};
	 fhistMCSpectrumAccLc->Fill(arr_FillkReco_Lc);
	 if(iSelTrackCuts){
     //fhistMCSpectrumAccLc->Fill(partpt,kRecoCuts+2,checkOrigin);// NOT FINISHED: SHOULD CHECK REFLECTIONS + GenAcc with MC variables!!
     const Double_t arr_FillkRecoCutsPlus2_Lc[4] = {partpt,kRecoCuts+2,(Double_t)checkOrigin,(Double_t)decay_channel};
     fhistMCSpectrumAccLc->Fill(arr_FillkRecoCutsPlus2_Lc);// NOT FINISHED: SHOULD CHECK REFLECTIONS + GenAcc with MC variables!!
   }
	 if(iSelCuts){
     //fhistMCSpectrumAccLc->Fill(partpt,kRecoCuts+3,checkOrigin);
     const Double_t arr_FillkRecoCutsPlus3_Lc[4] = {partpt,kRecoCuts+3,(Double_t)checkOrigin,(Double_t)decay_channel};
     fhistMCSpectrumAccLc->Fill(arr_FillkRecoCutsPlus3_Lc);
   }
   if(iSelPID && (massHypothesis==3 || (massHypothesis ==1 && (isTrueLambdaCorXic==40 || isTrueLambdaCorXic==50)) || (massHypothesis ==2 && (isTrueLambdaCorXic==80 || isTrueLambdaCorXic==100)) )){
     //fhistMCSpectrumAccLc->Fill(partpt,kRecoCuts+4,checkOrigin);
     const Double_t arr_FillkRecoCutsPlus4_Lc[4] = {partpt,kRecoCuts+4,(Double_t)checkOrigin,(Double_t)decay_channel};
     fhistMCSpectrumAccLc->Fill(arr_FillkRecoCutsPlus4_Lc);
   }
	 if(iSelTrackCuts>0&&massHypothesis>0){
	   if(massHypothesis==3 || (massHypothesis ==1 && (isTrueLambdaCorXic==40 || isTrueLambdaCorXic==50)) || (massHypothesis ==2 && (isTrueLambdaCorXic==80 || isTrueLambdaCorXic==100)) ){
	     //fhistMCSpectrumAccLc->Fill(partpt,kRecoPID,checkOrigin);
       const Double_t arr_FillkRecoPID_Lc[4] = {partpt,kRecoPID,(Double_t)checkOrigin,(Double_t)decay_channel};
	     fhistMCSpectrumAccLc->Fill(arr_FillkRecoPID_Lc);
	   }
	 }
	 // SIGMA C
	 Int_t indSc=part->GetMother();
	 if(indSc>=0){
	   mcpartMum=(AliAODMCParticle*)fmcArray->At(indSc); 
	   Int_t pdgLcMum=TMath::Abs(mcpartMum->GetPdgCode());
	   if(pdgLcMum==4112 || pdgLcMum==4222)isFromSigmaC=kTRUE;
	 }
	 //Double_t pointLcSc[6];
	 Double_t pointLcSc[7]; // adding axis for Lc decay channel (MC)
	 Double_t ptpartSc;
	 Double_t ypartSc;
	 if(isFromSigmaC){
	   ptpartSc=mcpartMum->Pt();
	   ypartSc=mcpartMum->Y();
	   pointLcSc[0]=partpt;
	   pointLcSc[1]=kRecoLc;
	   pointLcSc[2]=checkOrigin;
	   pointLcSc[3]=party;
	   pointLcSc[4]=ptpartSc;
	   pointLcSc[5]=ypartSc;
     pointLcSc[6]=decay_channel;
	   
	   fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	   
	   if(iSelTrackCuts){
	     pointLcSc[1]=kRecoLcCuts+2;
	     fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	   }
	   
	   if(iSelCuts){
	     pointLcSc[1]=kRecoLcCuts+3;
	     fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	   }

	   if(iSelPID && (massHypothesis==3 || (massHypothesis ==1 && (isTrueLambdaCorXic==40 || isTrueLambdaCorXic==50)) || (massHypothesis ==2 && (isTrueLambdaCorXic==80 || isTrueLambdaCorXic==100)) )){
	     pointLcSc[1]=kRecoLcCuts+4;
	     fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	   }

	   if(iSelTrackCuts>0&&massHypothesis>0){
	     if(massHypothesis==3 || (massHypothesis ==1 && (isTrueLambdaCorXic==40 || isTrueLambdaCorXic==50)) || (massHypothesis ==2 && (isTrueLambdaCorXic==80 || isTrueLambdaCorXic==100)) ){
	       pointLcSc[1]=kRecoLcPID;
	       fhistMCSpectrumAccLcFromSc->Fill(pointLcSc);
	     }
	   }
	 }// end is from sigmaC
	 
	 if(iSelTrackCuts>0&&massHypothesis>0){
	   fDist12SignalFilter->Fill(d->GetDist12toPrim()*10000.);
	   fCosPointDistrSignalFilter->Fill(d->CosPointingAngle());
	 }
       }
       if(iSelTrackCuts>0 && massHypothesis>0){
	 if(fDebug>=0){
	   Printf("Lc cand after cuts with masses: %f ,  %f, IsSeleCuts %d, iSelPID %d, massHypothesis %d",d->InvMassLcpKpi(),d->InvMassLcpiKp(),iSelCuts,iSelPID,massHypothesis);
	   Printf("pid: 1st track: %d, %d, %d",iSelDebugPion[0],iSelDebugKaon[0],iSelDebugProton[0]);
	   Printf("pid: 2nd track: %d, %d, %d",iSelDebugPion[1],iSelDebugKaon[1],iSelDebugProton[1]);
	   Printf("pid: 3rd track: %d, %d, %d",iSelDebugPion[2],iSelDebugKaon[2],iSelDebugProton[2]);
	   PrintCandidateVariables(d,aod);
	   Printf("\n \n");
	 }
	 fDist12AllFilter->Fill(d->GetDist12toPrim()*10000.);
	 fDist23AllFilter->Fill(d->GetDist23toPrim()*10000.);
	 fCosPointDistrAllFilter->Fill(d->CosPointingAngle());

   if(massHypothesis) fCandCounter->Fill(7);  // after IsSelected

       }
       if(!fSigmaCfromLcOnTheFly){
	 Int_t resp_onlyPID = 3;       	
	 Int_t isPIDused = fCutsXic->GetIsUsePID();
	 fCutsXic->SetUsePID(kFALSE);   // disable PID temporarly
	 Int_t resp_onlyCuts = fCutsXic->IsSelected(d,AliRDHFCuts::kCandidate,(AliAODEvent*)aod);
	 if(isPIDused){  // if the PID is supposed to be used, let's restore it in the cutobject
	   fCutsXic->SetUsePID(kTRUE);  // restoring PID
	   resp_onlyPID = fCutsXic->IsSelected(d,AliRDHFCuts::kPID,(AliAODEvent*)aod);
	 }
	 fCutsXic->SetUsePID(isPIDused); 
       
	 if(!fExplore_PIDstdCuts && massHypothesis>0){
	   if(isFromSigmaC){
	     SigmaCloop(d,aod,massHypothesis,d->InvMassLcpKpi(),d->InvMassLcpiKp(),0x0,resp_onlyPID,0x0,0x0,-1,-1,-1,mcpartMum,(Double_t)checkOrigin,(Double_t)decay_channel);
	   }
	   else if(!fReadMC){
	     SigmaCloop(d,aod,massHypothesis,d->InvMassLcpKpi(),d->InvMassLcpiKp(),0x0,resp_onlyPID);
	   }
	 }
	 massHypothesis=resp_onlyCuts & iSelPID;
	 if(fExplore_PIDstdCuts && massHypothesis>0){
     int nCasesExplPID = 10;
     if(fExplPID_BayesOnlyProt)  nCasesExplPID = 11;
     const int arr_dim_PID = nCasesExplPID+1;	   
	   Bool_t arrayPIDpkpi[arr_dim_PID],arrayPIDpikp[arr_dim_PID];	   
	   arrayPIDpkpi[0]=((massHypothesis&resp_onlyPID)==1 || (massHypothesis&resp_onlyPID)==3 ) ? kTRUE : kFALSE;
	   arrayPIDpikp[0]=( (massHypothesis&resp_onlyPID)==2 || (massHypothesis&resp_onlyPID)==3 ) ? kTRUE : kFALSE;
	   for(UInt_t i=1; i<=nCasesExplPID; i++){  // loop on PID cut combinations to be tested
	     arrayPIDpkpi[i]=kFALSE;
	     arrayPIDpikp[i]=kFALSE;	       
	     fCutsXic->ExplorePID(fPidResponse,d,i,arrayPIDpkpi[i],arrayPIDpikp[i]);
	   }
	   if(isFromSigmaC){
	     SigmaCloop(d,aod,massHypothesis,d->InvMassLcpKpi(),d->InvMassLcpiKp(),0x0,resp_onlyPID,arrayPIDpkpi,arrayPIDpikp,-1,-1,-1,mcpartMum,(Double_t)checkOrigin,(Double_t)decay_channel);
	   }
	   else if(!fReadMC){
	     SigmaCloop(d,aod,massHypothesis,d->InvMassLcpKpi(),d->InvMassLcpiKp(),0x0,resp_onlyPID,arrayPIDpkpi,arrayPIDpikp);	   
	   }
	 } 
       }
       if(recPrimVtx)fCutsXic->CleanOwnPrimaryVtx(d,aod,origownvtx);
       if(unsetvtx)d->UnsetOwnPrimaryVtx();
  }
}// END OF LOOP OVER FILTERED CANDIDARES
