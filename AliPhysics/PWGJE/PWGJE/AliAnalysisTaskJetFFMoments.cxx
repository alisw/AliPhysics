// ***********************************************************
// Task used for Fragmentation Function Moments (FFM) analysis
// Compares (gen) and (rec) jets   
// ************************************************************

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
 
#include <TFile.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TFormula.h>
#include <TF1.h>
#include <TRandom3.h>

#include "AliAnalysisTaskJetFFMoments.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliAODJetEventBackground.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliLog.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHerwigEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliAODMCHeader.h"
#ifndef __CINT__
#include "fastjet/PseudoJet.hh"
#include "fastjet/PseudoJetStructureBase.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
// get info on how fastjet was configured
#include "fastjet/config.h"

#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/JetFFMoments.hh"
namespace cont = fastjet::contrib;
#else
namespace fastjet {
  class PseudoJet;
  class PseudoJetStructureBase;
  class ClusterSequenceArea;  
  class GhostedAreaSpec;
  class AreaType;
  class JetDefinition;
  class JetAlgorithm;
  class Strategy;
  class RecombinationScheme;
  class Subtractor;
}
namespace  cont = fastjet::contrib {
  class JetFFMoments;
}
#endif

  using std::vector;

ClassImp(AliAnalysisTaskJetFFMoments)

AliAnalysisTaskJetFFMoments::~AliAnalysisTaskJetFFMoments()
{

  //
  // Destructor
  //

  delete fRef;

  for(int iJetBranch=0; iJetBranch < fgkFFMNJetBranches; iJetBranch++){
    fListJets[iJetBranch]->Clear(ProofClearOpt());
    if (!IsProof()) { 
      fListMatchedJets[iJetBranch]->Clear();
      fHistListJets[iJetBranch]->Clear();
      delete fListMatchedJets[iJetBranch];
    }
    if(fListJets[iJetBranch]) delete fListJets[iJetBranch];
  }

  if(fTracksAODMCCharged) delete fTracksAODMCCharged;
  if(fTracksAODMCChargedSecNS) delete fTracksAODMCChargedSecNS;
  if(fTracksAODMCChargedSecS)  delete fTracksAODMCChargedSecS;
  if(fTracksRecQualityCuts)    delete fTracksRecQualityCuts;

  if (fHistList && ! IsProof() ) { delete fHistList; }

  if( !fkIsPbPb ) {
    delete fh1CentralitySelect;
    delete fh1CentralityPhySel;
  }

  if(fTCAJetsOut)fTCAJetsOut->Delete();
  delete fTCAJetsOut;
  
  if(fAODJetBackgroundOut)fAODJetBackgroundOut->Reset();
  delete fAODJetBackgroundOut;

}

//---------------------------------------------------------------
AliAnalysisTaskJetFFMoments::AliAnalysisTaskJetFFMoments(): 
  AliAnalysisTaskSE(),
  fAOD(0x0),
  fAODJets(0x0),
  fAODExtension(0x0),
  fFilterPt(0.),
  fRef(0x0),
  fkIsPbPb(kFALSE),
  fkEventSelection(kTRUE),
  fkRejectFastOnly(kFALSE),
  fkRequireVZEROAC(kFALSE),
  fkRequireTZEROvtx(kFALSE),
  fkRejectPileup(kFALSE),
  fCentCutUp(0),
  fCentCutLo(0),
  fVtxZMax(8),
  fVtxR2Max(1),
  fFilterMask(0),
  fFilterMaskBestPt(0),
  fFilterType(0),
  fkUseAODTrackInput(kFALSE),
  fkUseAODMCInput(kFALSE),
  fTrackPtMin(0.15),
  fTrackPtMax(100.),
  fExtTrackCutType(0),
  fkUseHFcuts(kFALSE),
  fkRequireITSRefit(kFALSE),
  fkApplySharedClusterCut(kFALSE),
  fTrackEtaMin(-0.9),
  fTrackEtaMax(0.9),
  fTrackPhiMin(0.),
  fTrackPhiMax(2*TMath::Pi()),
  fTCAJetsOut(0x0),
  fAODJetBackgroundOut(0x0),
  fTracksAODMCCharged(0x0),
  fTracksAODMCChargedSecNS(0x0),
  fTracksAODMCChargedSecS(0x0),
  fTracksRecQualityCuts(0x0),
  fNonStdBranch(""),
  fNonStdFile(""),
  fAnaJetType("ALLJET"),
  fJetPtMin(2.),
  fJetEtaMin(-0.5),
  fJetEtaMax(0.5),
  fJetDeltaPhiCut(TMath::Pi()/3.),
  fkDoJetMatching(kFALSE),
  fkUseClosestJetsforMatching(kFALSE),
  fkFillMismatchHisto(kFALSE),
  fJetMatchingFractionMin(0.5),
  fJetMatchedDistMax(0.3),
  fkUseJetFromInput(kFALSE),
  fNUsedJets(999),
  fJetMinLTrackPt(-1),
  fJetMaxTrackPt(-1),
  fJetMinnTracks(0),
  fkEffJetType(0),
  fkGenJetType(0),
  fTracksInJetMethod(0),
  fFFBckgMode(0),
  fFFJtValue(0),
  fFFMNMin(-0.75),
  fFFMNMax(6.0),
  fFFMMomMax(28),
  fFFMScalePower(4),
  fFFMBckgType(1),
  fFFMBckgPar1(0.9),
  fFFMBckgPar2(TMath::Pi()),
  fFFMBckgMu(25),
  fHistosLevel(1),
  fkHighResolution(kFALSE),
  fkUseTrackPtSumAsJetPt(kFALSE),
  fkDoJetReco(kFALSE),
  fkUseBackgroundCalc(kFALSE),
  fRparam(0.4),
  fAlgorithm(fastjet::antikt_algorithm),
  fBkgAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area),
  fGhostArea(0.01),
  fActiveAreaRepeats(1),
  fGhostEtaMax(1.5),
  fPrParameters(0x0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1AvgTrials(0x0),
  fh1CentralityPhySel(0x0),
  fh1CentralitySelect(0x0),
  fh1vZPhySel(0x0),
  fh1vZSelect(0x0),
  fEffi(0x0),
  fResol(0.01),
  fResolMeth(0),
  fEffivar(0),
  fResolvar(0),
  fPtHardAndPythiaJetPtFactor(0),
  fPtHardAndTrackPtFactor(0)
{
  //
  // Default Constructor
  //

  fTrackType[0] = kTrackUndef;
  fTrackType[1] = kTrackUndef;
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) { fh1GenTracks[iAxis] = 0x0; fh1RecTracks[iAxis] = 0x0; fh2PtRecVsGenPrim[iAxis] = 0x0; fh2PtRecVsGenSec[iAxis] = 0x0; }
  for(int iJetBranch=0; iJetBranch < fgkFFMNJetBranches; iJetBranch++){
    if(iJetBranch < fgkFFMNJetBranches-1){
      for( Int_t iAxis = 0; iAxis < 5; iAxis++)  fh2MatchedJets[iAxis] = 0x0;
      fhnJetMomN_Raw = 0x0;
      fhnJetMomN_Sub = 0x0;
      fhnJetMomN_Imp = 0x0;
    }
    fJetBranch[iJetBranch] ="";
    fBkgJetBranch[iJetBranch] ="";
    fListJets[iJetBranch] = new TList(); 
    fListJets[iJetBranch]->SetOwner(kTRUE);
    fListMatchedJets[iJetBranch] = new TList();
    fListMatchedJets[iJetBranch]->SetOwner(kTRUE);
    fh1Njets[iJetBranch] = 0x0;
    if( fAnaJetType.Contains("DIJET")) fh1Asy_DiJets[iJetBranch] = 0x0;
	fh2MatchedJetsRDPtVSPt[iJetBranch] = 0x0;
	fh2MatchedJetsAreaVSPt[iJetBranch] = 0x0;
        fh2MatchedJetsUE[iJetBranch] = 0x0;
    for( Int_t iAxis = 0; iAxis < 5; iAxis++) fh1JetPr_Mismatched[iJetBranch][iAxis] = 0x0;
	fh2MismatchedJetsAreaVSPt[iJetBranch] = 0x0;
	for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2TracksInJets[iJetBranch][iAxis] = 0x0;
	for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2AssociatedTracksInJets[iJetBranch][iAxis] = 0x0;
	for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2AssociatedTracksInJetsSecNS[iJetBranch][iAxis] = 0x0;
	for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2AssociatedTracksInJetsSecS[iJetBranch][iAxis] = 0x0;
	for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2AssociatedTracksInJetsSecSsc[iJetBranch][iAxis] = 0x0;
    fhnJetFFM_Raw[iJetBranch] = 0x0;
    fhnJetFFM_Sub[iJetBranch] = 0x0;
    fhnJetFFM_Imp[iJetBranch] = 0x0;
        fp2TracksInJetFFM[iJetBranch] = 0x0;
	fp2AssociatedTracksJetFFM[iJetBranch] = 0x0;
	fp2AssociatedTracksJetFFMSecNS[iJetBranch] = 0x0;
	fp2AssociatedTracksJetFFMSecS[iJetBranch] = 0x0;
	fp2AssociatedTracksJetFFMSecSsc[iJetBranch] = 0x0;
	fp2JetFFM_Raw[iJetBranch] = 0x0;
	fp2JetFFM_Sub[iJetBranch] = 0x0;
	fp2JetFFM_Imp[iJetBranch] = 0x0;
    fHistListJets[iJetBranch] = 0x0;
  }

  Double_t pi = TMath::Pi();
  Double_t halfWidth = (fFFMNMax - fFFMNMin)/(fFFMMomMax-1)/2.;
  Double_t nAxisMin = fFFMNMin - halfWidth;
  Double_t nAxisMax = fFFMNMax + halfWidth;

  // set Axis parameters, if not set in sepcial:
  // 0 - 9
  //  0 |  1  |  2  |  3   |  4  |  5 |   6   |      7      |    8    |    9    |
  // vz | ntr | ep  | epb  |  z  | xi | lnjT  |  DeltaTheta | FFM_gen | FFM_rec |
  // 10 - 24
  //  0 |  1  |  2  |    3     |      4        |||      n
  // pt | eta | phi |   area   | nconstituents |||     gen  2x5+n
  // pt | eta | phi |   area   | nconstituents |||     rec  3x5+n
  // pt | eta | phi |          |    N(order)   |||    track 4x5+n
  // 25 - 32
  //	 25	|	26	|
  //  fraction | RDeltaPt |  R for Relative
  //
  fnBinsAxis[0] = 20;          fBinMinAxis[0] = 0.;            fBinMaxAxis[0] = 100.;
  fnBinsAxis[1] = 50;          fBinMinAxis[1] = 0.;            fBinMaxAxis[1] = 50.;
  fnBinsAxis[2] = 40;          fBinMinAxis[2] = -20.;          fBinMaxAxis[2] = 20.;
  fnBinsAxis[3] = 30;          fBinMinAxis[3] = 0.;            fBinMaxAxis[3] = pi;
  fnBinsAxis[4] = 44;          fBinMinAxis[4] = 0.;            fBinMaxAxis[4] = 1.1;
  fnBinsAxis[5] = 28;          fBinMinAxis[5] = 0.;            fBinMaxAxis[5] = 7.;
  fnBinsAxis[6] = 25;          fBinMinAxis[6] = 4.-2*pi;       fBinMaxAxis[6] = 4.;
  fnBinsAxis[7] = 20;          fBinMinAxis[7] = 0.;            fBinMaxAxis[7] = 0.5;
  fnBinsAxis[8] =675;		  fBinMinAxis[8] = 0.;			fBinMaxAxis[8] = 30.;
  fnBinsAxis[9] =675;		  fBinMinAxis[9] = 0.;			fBinMaxAxis[9] = 30.;
  for( Int_t i = 0; i < 2; i ++) {
    fnBinsAxis[10+i*5] = 40;   fBinMinAxis[10+i*5] = 0.;       fBinMaxAxis[10+i*5] = 200.;
    fnBinsAxis[11+i*5] = 100;   fBinMinAxis[11+i*5] = -1.;      fBinMaxAxis[11+i*5] = 1.;
    fnBinsAxis[12+i*5] = 90;   fBinMinAxis[12+i*5] = 0.;       fBinMaxAxis[12+i*5] = 2*pi;
    fnBinsAxis[13+i*5] = 50;   fBinMinAxis[13+i*5] = 0.;	   fBinMaxAxis[13+i*5] = 1.;
    fnBinsAxis[14+i*5] = 50;   fBinMinAxis[14+i*5] = 0.5;	  fBinMaxAxis[14+i*5] = 50.5;
  }
  fnBinsAxis[20] = 400;		fBinMinAxis[20] = 0.;		   fBinMaxAxis[20] = 200.;
  fnBinsAxis[21] = 80;		 fBinMinAxis[21] = -1.;		  fBinMaxAxis[21] = 1.;
  fnBinsAxis[22] = 90;		 fBinMinAxis[22] = 0.;		   fBinMaxAxis[22] = 2*pi;
  fnBinsAxis[23] = 0;  fBinMinAxis[23] = 0.;           fBinMaxAxis[23] = 0.;
  fnBinsAxis[24] = fFFMMomMax; fBinMinAxis[24] = nAxisMin;    fBinMaxAxis[24] = nAxisMax;
  fnBinsAxis[25] = 21; fBinMinAxis[25] = 0.;		   fBinMaxAxis[25] = 1.05;
  fnBinsAxis[26] = 80; fBinMinAxis[26] = -1.;		   fBinMaxAxis[26] = 0.6;
  fnBinsAxis[27] = 400; fBinMinAxis[27] = -10.;              fBinMaxAxis[27] = 10;
  fnBinsAxis[28] = 500; fBinMinAxis[28] = 0.;              fBinMaxAxis[28] = 50;

  for( Int_t i = 29; i < 32; i ++) {
    fnBinsAxis[i] = 0; fBinMinAxis[i] = 0.;    fBinMaxAxis[i] = 0.;
  }

  SetCentralityClasses();

  fHistList = 0x0;
}

//---------------------------------------------------------------
AliAnalysisTaskJetFFMoments::AliAnalysisTaskJetFFMoments(const char* name):
  AliAnalysisTaskSE(name),
  fAOD(0x0),
  fAODJets(0x0),
  fAODExtension(0x0),
  fFilterPt(0.),
  fRef(new TRefArray),
  fkIsPbPb(kFALSE),
  fkEventSelection(kTRUE),
  fkRejectFastOnly(kFALSE),
  fkRequireVZEROAC(kFALSE),
  fkRequireTZEROvtx(kFALSE),
  fkRejectPileup(kFALSE),
  fCentCutUp(0),
  fCentCutLo(0),
  fVtxZMax(8),
  fVtxR2Max(1),
  fFilterMask(0),
  fFilterMaskBestPt(0),
  fFilterType(0),
  fkUseAODTrackInput(kFALSE),
  fkUseAODMCInput(kFALSE),
  fTrackPtMin(0.15),
  fTrackPtMax(100.),
  fExtTrackCutType(0),
  fkUseHFcuts(kFALSE),
  fkRequireITSRefit(kFALSE),
  fkApplySharedClusterCut(kFALSE),
  fTrackEtaMin(-0.9),
  fTrackEtaMax(0.9),
  fTrackPhiMin(0.),
  fTrackPhiMax(2*TMath::Pi()),
  fTCAJetsOut(0x0),
  fAODJetBackgroundOut(0x0),
  fTracksAODMCCharged(0x0),
  fTracksAODMCChargedSecNS(0x0),
  fTracksAODMCChargedSecS(0x0),
  fTracksRecQualityCuts(0x0),
  fNonStdBranch(""),
  fNonStdFile(""),
  fAnaJetType("ALLJET"),
  fJetPtMin(2.),
  fJetEtaMin(-0.5),
  fJetEtaMax(0.5),
  fJetDeltaPhiCut(TMath::Pi()/3.),
  fkDoJetMatching(kFALSE),
  fkUseClosestJetsforMatching(kFALSE),
  fkFillMismatchHisto(kFALSE),
  fJetMatchingFractionMin(0.5),
  fJetMatchedDistMax(0.3),
  fkUseJetFromInput(kFALSE),
  fNUsedJets(999),
  fJetMinLTrackPt(-1),
  fJetMaxTrackPt(-1),
  fJetMinnTracks(0),
  fkEffJetType(0),
  fkGenJetType(0),
  fTracksInJetMethod(0),
  fFFBckgMode(0),
  fFFJtValue(0),
  fFFMNMin(-0.75),
  fFFMNMax(6.0),
  fFFMMomMax(28),
  fFFMScalePower(4),
  fFFMBckgType(1),
  fFFMBckgPar1(0.9),
  fFFMBckgPar2(TMath::Pi()),
  fFFMBckgMu(25),
  fHistosLevel(1),
  fkHighResolution(kFALSE),
  fkUseTrackPtSumAsJetPt(kFALSE),
  fkDoJetReco(kFALSE),
  fkUseBackgroundCalc(kFALSE),
  fRparam(0.4),
  fAlgorithm(fastjet::antikt_algorithm),
  fBkgAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area),
  fGhostArea(0.01),
  fActiveAreaRepeats(1),
  fGhostEtaMax(1.5),
  fPrParameters(0x0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1AvgTrials(0x0),
  fh1CentralityPhySel(0x0),
  fh1CentralitySelect(0x0),
  fh1vZPhySel(0x0),
  fh1vZSelect(0x0),
  fEffi(0x0),
  fResol(0.01),
  fResolMeth(0),
  fEffivar(0),
  fResolvar(0),
  fPtHardAndPythiaJetPtFactor(0),
  fPtHardAndTrackPtFactor(0)
{
  //
  // named ctor
  //

  fTrackType[0] = kTrackUndef;
  fTrackType[1] = kTrackUndef;
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) { fh1GenTracks[iAxis] = 0x0; fh1RecTracks[iAxis] = 0x0; fh2PtRecVsGenPrim[iAxis] = 0x0; fh2PtRecVsGenSec[iAxis] = 0x0; }
  for(int iJetBranch=0; iJetBranch < fgkFFMNJetBranches; iJetBranch++){
    if(iJetBranch < fgkFFMNJetBranches-1){
      for( Int_t iAxis = 0; iAxis < 5; iAxis++)  fh2MatchedJets[iAxis] = 0x0;
      fhnJetMomN_Raw = 0x0;
      fhnJetMomN_Sub = 0x0;
      fhnJetMomN_Imp = 0x0;
    }
    fJetBranch[iJetBranch] ="";
    fBkgJetBranch[iJetBranch] ="";
    fListJets[iJetBranch] = new TList; 
    fListJets[iJetBranch]->SetOwner(kTRUE);
    fListMatchedJets[iJetBranch] = new TList();
    fListMatchedJets[iJetBranch]->SetOwner(kTRUE);  
    fh1Njets[iJetBranch] = 0x0;
    if( fAnaJetType.Contains("DIJET")) fh1Asy_DiJets[iJetBranch] = 0x0;
    fh2MatchedJetsRDPtVSPt[iJetBranch] = 0x0;
    fh2MatchedJetsAreaVSPt[iJetBranch] = 0x0;
    fh2MatchedJetsUE[iJetBranch] = 0x0;
    for( Int_t iAxis = 0; iAxis < 5; iAxis++) fh1JetPr_Mismatched[iJetBranch][iAxis] = 0x0;
    fh2MismatchedJetsAreaVSPt[iJetBranch] = 0x0;
    for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2TracksInJets[iJetBranch][iAxis] =0x0;
    for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2AssociatedTracksInJets[iJetBranch][iAxis] =0x0;
    for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2AssociatedTracksInJetsSecNS[iJetBranch][iAxis] =0x0;
    for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2AssociatedTracksInJetsSecS[iJetBranch][iAxis] =0x0;
    for( Int_t iAxis = 0; iAxis < 7; iAxis++) fh2AssociatedTracksInJetsSecSsc[iJetBranch][iAxis] =0x0;
    fhnJetFFM_Raw[iJetBranch] = 0x0;
    fhnJetFFM_Sub[iJetBranch] = 0x0;
    fhnJetFFM_Imp[iJetBranch] = 0x0;
    fp2TracksInJetFFM[iJetBranch] = 0x0;
    fp2AssociatedTracksJetFFM[iJetBranch] = 0x0;
    fp2AssociatedTracksJetFFMSecNS[iJetBranch] = 0x0;
    fp2AssociatedTracksJetFFMSecS[iJetBranch] = 0x0;
    fp2AssociatedTracksJetFFMSecSsc[iJetBranch] = 0x0;
    fp2JetFFM_Raw[iJetBranch] = 0x0;
    fp2JetFFM_Sub[iJetBranch] = 0x0;
    fp2JetFFM_Imp[iJetBranch] = 0x0;
    fHistListJets[iJetBranch] = 0x0;
  }

  Double_t pi = TMath::Pi();
  Double_t halfWidth = (fFFMNMax - fFFMNMin)/(fFFMMomMax-1)/2.;
  Double_t nAxisMin = fFFMNMin - halfWidth;
  Double_t nAxisMax = fFFMNMax + halfWidth;
  
  // set Axis parameters, if not set in sepcial:
  // 0 - 9
  //  0 |  1  |  2  |  3   |  4  |  5 |   6   |      7      |    8    |    9    |
  // vz | ntr | ep  | epb  |  z  | xi | lnjT  |  DeltaTheta | FFM_gen | FFM_rec |
  // 10 - 24
  //  0 |  1  |  2  |    3     |      4        |||      n
  // pt | eta | phi |   area   | nconstituents |||     gen  2x5+n
  // pt | eta | phi |   area   | nconstituents |||     rec  3x5+n
  // pt | eta | phi |          |    N(order)   |||    track 4x5+n
  // 25 - 32
  //	 25	|	26	|
  //  fraction | RDeltaPt |  R for Relative
  //
  fnBinsAxis[0] = 20;          fBinMinAxis[0] = 0.;            fBinMaxAxis[0] = 100.;
  fnBinsAxis[1] = 50;          fBinMinAxis[1] = 0.;            fBinMaxAxis[1] = 50.;
  fnBinsAxis[2] = 40;          fBinMinAxis[2] = -20.;          fBinMaxAxis[2] = 20.;
  fnBinsAxis[3] = 30;          fBinMinAxis[3] = 0.;            fBinMaxAxis[3] = pi;
  fnBinsAxis[4] = 44;          fBinMinAxis[4] = 0.;            fBinMaxAxis[4] = 1.1;
  fnBinsAxis[5] = 28;          fBinMinAxis[5] = 0.;            fBinMaxAxis[5] = 7.;
  fnBinsAxis[6] = 25;          fBinMinAxis[6] = 4.-2*pi;       fBinMaxAxis[6] = 4.;
  fnBinsAxis[7] = 20;          fBinMinAxis[7] = 0.;            fBinMaxAxis[7] = 0.5;
  fnBinsAxis[8] =150;          fBinMinAxis[8] = 0.;            fBinMaxAxis[8] = 30.;
  fnBinsAxis[9] =150;          fBinMinAxis[9] = 0.;            fBinMaxAxis[9] = 30.;
  for( Int_t i = 0; i < 2; i ++) {
    fnBinsAxis[10+i*5] = 40;   fBinMinAxis[10+i*5] = 0.;       fBinMaxAxis[10+i*5] = 200.;
    fnBinsAxis[11+i*5] = 100;   fBinMinAxis[11+i*5] = -1.;      fBinMaxAxis[11+i*5] = 1.;
    fnBinsAxis[12+i*5] = 90;   fBinMinAxis[12+i*5] = 0.;       fBinMaxAxis[12+i*5] = 2*pi;
    fnBinsAxis[13+i*5] = 50;   fBinMinAxis[13+i*5] = 0.;	   fBinMaxAxis[13+i*5] = 1.;
    fnBinsAxis[14+i*5] = 50;   fBinMinAxis[14+i*5] = 0.5;	  fBinMaxAxis[14+i*5] = 50.5;
  }
  fnBinsAxis[20] = 400;		fBinMinAxis[20] = 0.;		   fBinMaxAxis[20] = 200.;
  fnBinsAxis[21] = 80;		 fBinMinAxis[21] = -1.;		  fBinMaxAxis[21] = 1.;
  fnBinsAxis[22] = 90;		 fBinMinAxis[22] = 0.;		   fBinMaxAxis[22] = 2*pi;
  fnBinsAxis[23] = 0;  fBinMinAxis[23] = 0.;           fBinMaxAxis[23] = 0.;
  fnBinsAxis[24] = fFFMMomMax; fBinMinAxis[24] = nAxisMin;    fBinMaxAxis[24] = nAxisMax;
  fnBinsAxis[25] = 21;		 fBinMinAxis[25] = 0.;		   fBinMaxAxis[25] = 1.05;
  fnBinsAxis[26] = 80;		 fBinMinAxis[26] = -1.;		  fBinMaxAxis[26] = 0.6;
  fnBinsAxis[27] = 400;          fBinMinAxis[27] = -10.;              fBinMaxAxis[27] = 10;
  fnBinsAxis[28] = 500;          fBinMinAxis[28] = 0.;              fBinMaxAxis[28] = 50;

  for( Int_t i = 29; i < 32; i ++) {
    fnBinsAxis[i] = 0; fBinMinAxis[i] = 0.;    fBinMaxAxis[i] = 0.;
  }

  SetCentralityClasses();

  fHistList = 0x0;
  DefineOutput(1, TList::Class());  
}

//----------------------------------------------------------------
void AliAnalysisTaskJetFFMoments::LocalInit()
{
  if(fDebug) {
    if(fkDoJetReco) printf("Jet finding ON - ");
    printf("Type of jet analysed: %s\n ", fAnaJetType.Data());
    printf("TRACK PARAMETERS:   FilterMask: %d	 pt_min: %.3f GeV/c	eta min max: [%.2f,%.2f]  extra cut type: %d\n", fFilterMask, fTrackPtMin, fTrackEtaMin, fTrackEtaMax, fExtTrackCutType); 
    printf(" JET   PARAMETERS:   Algorithm: %d	 R: %.2f                eta min max: [%.2f,%.2f]  pt_min: %.3f \n", fAlgorithm, fRparam, fJetEtaMin, fJetEtaMax,fJetPtMin);
    if(fTrackType[0]>0) printf(" Efficiency Params: GenJetType: %d EffJetType: %d\n", fkGenJetType,fkEffJetType);
    if(fTrackType[1] < kTrackAOD) printf(" FastSim Mode ON (Toy MC - KINE): TrackResolution %.3f (method %d)  Systematics: TrackResolution %.1f %% Efficiency %.1f", fResol, fResolMeth ,fEffivar, fResolvar);
    printf(" TRACKS IN JETS:   Method: %d\n",fTracksInJetMethod);
    if(fFFBckgMode)     printf("FF Background Mode %d\n",fFFBckgMode);
    if(fkDoJetMatching) printf("Jet matching ON with parameters, distance: %.2f, energy fraction: %.3f\n", fJetMatchedDistMax, fJetMatchingFractionMin);
    printf("FFM computed with scaling power %.1f, background mode %d and parameters: %.2f %.2f %.3f\n", fFFMScalePower, fFFMBckgType, fFFMBckgPar1, fFFMBckgPar2, fFFMBckgMu);
  }
  if( fkDoJetMatching == kFALSE) {
    fJetMatchedDistMax = 0.; fJetMatchingFractionMin = 0.;
  }

  if((!fJetBranch[0].Contains("MC") && !fkDoJetReco) || (fkDoJetReco && !(fTrackType[0] == kTrackAODMCCharged || fTrackType[0] ==  kTrackAODMCChargedAcceptance  || fTrackType[0] == kTrackAODMCextra || fTrackType[0] == kTrackAODMCextraonly || fTrackType[0]< kTrackAOD))){
    AliWarning("Warning! You are running with the reconstructed data only. Both axis of the correlation plots will be filled with the reconstructed information.");
    if (fkDoJetMatching) { AliError("The Jet Matching cannot be switched on if we do not have MC data"); exit(1);}
  }
}

//---------------------------------------------------------------
Bool_t AliAnalysisTaskJetFFMoments::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 

  if(fDebug != 0) printf("AliAnalysisTaskJetFFMoments::Notify() \n");

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ntrials  = 1;

  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    if(!fh1Xsec||!fh1Trials){
      if(fDebug > 2) Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }

    if(!((!fJetBranch[0].Contains("MC") && !fkDoJetReco) || (fkDoJetReco && !(fTrackType[0] == kTrackAODMCCharged || fTrackType[0] ==  kTrackAODMCChargedAcceptance  || fTrackType[0] == kTrackAODMCextra || fTrackType[0] == kTrackAODMCextraonly)))){
      AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(), xsection, ntrials);
      fh1Xsec->Fill("<#sigma>", xsection);
      fh1Trials->Fill("#sum{ntrials}", ntrials);
      // construct a poor man average trials 
      Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
      Float_t avgTrials = -999.;
      if(ntrials>=nEntries && nEntries>0.) avgTrials = ntrials/nEntries;
      fh1AvgTrials->Fill("#sum{ntrials}", avgTrials);
    }
  }
  return kTRUE;
}

//---------------------------------------------------------------
void AliAnalysisTaskJetFFMoments::UserCreateOutputObjects()
{

  //
  // Create the output container
  //

  // Connect the AOD
  if (fDebug != 0) printf("AnalysisTaskJetFFMoments::UserCreateOutputObjects() \n");
  if(fNonStdBranch.Length()!=0)
    {
      // only create the output branch if we have a name
      // Create a new branch for jets...
      //  -> cleared in the UserExec....
      // here we can also have the case that the brnaches are written to a separate file
      
      fTCAJetsOut = new TClonesArray("AliAODJet", 0);
      fTCAJetsOut->SetName(fNonStdBranch.Data());
      AddAODBranch("TClonesArray", &fTCAJetsOut, fNonStdFile.Data());

      if(fkUseBackgroundCalc){
        if(!AODEvent()->FindListObject(Form("%s_%s", AliAODJetEventBackground::StdBranchName(), fNonStdBranch.Data()))){
          fAODJetBackgroundOut = new AliAODJetEventBackground();
          fAODJetBackgroundOut->SetName(Form("%s_%s", AliAODJetEventBackground::StdBranchName(), fNonStdBranch.Data()));
          AddAODBranch("AliAODJetEventBackground", &fAODJetBackgroundOut, fNonStdFile.Data());  
        } 
      }
    
      if(fNonStdFile.Length()!=0){
        // 
        // case that we have an AOD extension we need to fetch the jets from the extended output
        // we identify the extension aod event by looking for the branchname
        AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
        // case that we have an AOD extension we need can fetch the background maybe from the extended output                                                                  
        fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
      }
    }

  fTracksAODMCCharged = new TList();
  fTracksAODMCCharged->SetOwner(kFALSE);
  fTracksAODMCCharged->SetName("fTracksAODMCCharged");

  fTracksAODMCChargedSecNS = new TList();
  fTracksAODMCChargedSecNS->SetOwner(kFALSE);
  fTracksAODMCChargedSecNS->SetName("fTracksAODMCChargedSecNS");

  fTracksAODMCChargedSecS = new TList();
  fTracksAODMCChargedSecS->SetOwner(kFALSE);
  fTracksAODMCChargedSecS->SetName("fTracksAODMCChargedSecS");

  fTracksRecQualityCuts = new TList();
  fTracksRecQualityCuts->SetOwner(kFALSE);

  fAnaJetType.ToUpper();
  if(!fHistList)fHistList = new TList();
  fHistList->SetOwner();
  PostData(1, fHistList); // post data in any case once

  TString fAnaObjectType = "";
  // Set the name of the 2 TLists stored
  for(int iJetBranch=0; iJetBranch < fgkFFMNJetBranches; iJetBranch++){ // output list
    if( iJetBranch == 0) {
      fAnaObjectType = "Gen";
    } else {
      fAnaObjectType = "Rec";
    }
    if(!fHistListJets[iJetBranch] ) { fHistListJets[iJetBranch] = new TList(); }
    fHistListJets[iJetBranch]->SetOwner(kTRUE);
    fHistListJets[iJetBranch]->SetName(Form("Ana%s_%s", fAnaJetType.Data(), fAnaObjectType.Data()));
  }

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
 
  //
  //  Histogram
  
  //reset the Axis parameters for FFM
  Double_t halfWidth = (fFFMNMax - fFFMNMin)/(fFFMMomMax-1)/2.;
  Double_t nAxisMin = fFFMNMin - halfWidth;
  Double_t nAxisMax = fFFMNMax + halfWidth;
  fnBinsAxis[24] = fFFMMomMax; fBinMinAxis[24] = nAxisMin;    fBinMaxAxis[24] = nAxisMax;

  CreateHistos();

  // Fill histogram with analysis parameters
  fPrParameters->Fill("DoJetReco", fkDoJetReco);
  fPrParameters->Fill("FFMScalePower", fFFMScalePower);
  fPrParameters->Fill("IsPbPb", fkIsPbPb);
  fPrParameters->Fill("DoMatching", fkDoJetMatching);
  fPrParameters->Fill("Kine_Input", fJetBranch[0].Length() != 0);
  fPrParameters->Fill("Reco_Input", fJetBranch[1].Length() != 0);
  fPrParameters->Fill("NUsed(di)Jets", fNUsedJets);
  fPrParameters->Fill("MatchMaxDist", fJetMatchedDistMax);
  fPrParameters->Fill("JetMatchingFractionMin", fJetMatchingFractionMin);
  fPrParameters->Fill("FFBckgMode", fFFBckgMode);
  fPrParameters->Fill("FFMBckgType", fFFMBckgType);
  fPrParameters->Fill("FFMBckgPar1", fFFMBckgPar1);
  fPrParameters->Fill("FFMBckgPar2", fFFMBckgPar2);
  fPrParameters->Fill("FFMBckgMu/10", fFFMBckgMu/10.);
  fPrParameters->Fill("JetAlgorithm", fAlgorithm);
  fPrParameters->Fill("JetR",fRparam);
  fPrParameters->Fill("JetPtMin", fJetPtMin);
  fPrParameters->Fill("|JetEtaMin|", TMath::Abs(fJetEtaMin));
  fPrParameters->Fill("JetEtaMax", fJetEtaMax);
  fPrParameters->Fill("FilterMask/100", fFilterMask/100.);
  fPrParameters->Fill("TrackPtMin", fTrackPtMin);
  fPrParameters->Fill("|TrackEtaMin|", TMath::Abs(fTrackEtaMin));
  fPrParameters->Fill("TrackEtaMax", fTrackEtaMax);
  fPrParameters->Fill("TrackExtraCut", fExtTrackCutType);
  fPrParameters->Fill("UseTrackPtSumAsJetPt", fkUseTrackPtSumAsJetPt);
  fPrParameters->Fill("TracksInJetMethod", fTracksInJetMethod);
  fPrParameters->Fill("GenJetType", fkGenJetType);
  fPrParameters->Fill("EffJetType", fkEffJetType);
  fPrParameters->Fill("FastTrackEfficiencySyst/100", fEffivar);
  fPrParameters->Fill("FastTrackResolution", fResol);
  fPrParameters->Fill("FastTrackResolutionSyst/100", fResolvar);

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fHistList->GetEntries(); ++i) { //
    TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
  }
  for(int iJetBranch=0; iJetBranch < fgkFFMNJetBranches; iJetBranch++){
    for (Int_t i=0; i<fHistListJets[iJetBranch]->GetEntries(); ++i) {
      TH1 *h1 = dynamic_cast<TH1*>(fHistListJets[iJetBranch]->At(i));
      if (h1){
        h1->Sumw2();
        continue;
      }
    }
  }
  TH1::AddDirectory(oldStatus);

}

//---------------------------------------------------------------
void AliAnalysisTaskJetFFMoments::UserExec(Option_t */*option*/)
{

  if(fDebug != 0) Printf("AliAnalysisTaskJetFFMoments::UserExec()");

  // handle and reset the output jet branch 
  if(fTCAJetsOut)fTCAJetsOut->Delete();
  if(fAODJetBackgroundOut)fAODJetBackgroundOut->Reset();

  // track sorts
  if(fTrackType[1] >= kTrackAOD) {
  if(fkUseAODTrackInput){
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager %d", (char*)__FILE__, __LINE__, fkUseAODTrackInput);
      return;
    }
    // fetch the header
  } else {
    //  assume that the AOD is in the general output...
    fAOD  = AODEvent();
    if(!fAOD){
      Printf("%s:%d AODEvent not found in the Output", (char*)__FILE__, __LINE__);
      return;
    }
  }
  if(!fAOD){
    Printf("%s:%d AODEvent not found", (char*)__FILE__,__LINE__);
    return;
  }
  //
  // Execute analysis for current event
  //
  Bool_t selectEvent =  false;
  Bool_t physicsSelection = true;// handled by the framework(fInputHandler->IsEventSelected()&AliVEvent::kMB)==AliVEvent::kMB;

   // Trigger selection
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)
   ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
   if( fkRejectFastOnly && ((inputHandler->IsEventSelected() & AliVEvent::kFastOnly) == AliVEvent::kFastOnly)){
   if (fDebug > 10) std::cout<<fAOD->GetFiredTriggerClasses()<<std::endl;
   physicsSelection = false;
   return;
   }

  Float_t cent = 0;
  Float_t zVtx  = 0;
  Int_t cenClass = -1;
  if(fAOD) {
    const AliAODVertex *vtxAOD = fAOD->GetPrimaryVertex();
    if(!vtxAOD) {
      Printf("%s:%d vtxAOD not found\n", (char*)__FILE__, __LINE__); 
      return;
    }
    TString vtxTitle(vtxAOD->GetTitle());
    zVtx = vtxAOD->GetZ();

    cent = ((AliAODHeader*)fAOD->GetHeader())->GetCentrality();
    if(fkIsPbPb){ 
      if(cent<fCentClass[0])cenClass = 0;
      else if(cent<fCentClass[1])cenClass = 1;
      else if(cent<fCentClass[2])cenClass = 2;
      else if(cent<fCentClass[3])cenClass = 3;
      std::cout << "Centrality class: " << cenClass << std::endl;
    }

    if(physicsSelection){
      fh1CentralityPhySel->Fill(cent);
      fh1vZPhySel->Fill(zVtx);
    }

    if(fkEventSelection){
      if(vtxAOD->GetNContributors()>2&&!vtxTitle.Contains("TPCVertex")){
        Float_t yvtx = vtxAOD->GetY();
        Float_t xvtx = vtxAOD->GetX();
        Float_t r2   = yvtx*yvtx+xvtx*xvtx;  
        if(TMath::Abs(zVtx)<fVtxZMax&&r2<fVtxR2Max){ // apply vertex cut later on
          if(physicsSelection){
            selectEvent = true;
          }
        }
      }
      if(fCentCutUp>0){
        if(cent<fCentCutLo||cent>fCentCutUp){
          selectEvent = false;
        }
      }
     if(fkRejectPileup && AliAnalysisHelperJetTasks::IsPileUp()){
        selectEvent = false;
     }
    }else{
      selectEvent = true;
    }
  }

  // The selection below is used for LHC10h (PbPb)
  Bool_t T0 = false;
  Bool_t V0 = false;
  const AliAODVZERO  *vzero = fAOD->GetVZEROData();
  if(vzero){
    if((vzero->GetTriggerChargeA()>0)&&(vzero->GetTriggerChargeC()>0)){
      V0 = true;
    }
  }
  
  const AliAODTZERO  *tzero = fAOD->GetTZEROData();
  if(tzero){
    if(TMath::Abs(tzero->GetT0VertexRaw())<100){
      T0 = true;
    }
  }
  
  if(fkRequireVZEROAC&&fkRequireTZEROvtx)selectEvent = selectEvent&&V0&&T0;
  else if(fkRequireTZEROvtx)selectEvent = selectEvent&&T0;
  else if(fkRequireVZEROAC)selectEvent = selectEvent&&V0;

 if(fPtHardAndPythiaJetPtFactor || fPtHardAndTrackPtFactor)  selectEvent = selectEvent&&!IsOutlier(GetPythiaHeader());

  if(!selectEvent){
    PostData(1, fHistList);
    return;
  }

  fh1CentralitySelect->Fill(cent);  
  fh1vZSelect->Fill(zVtx);

} else  {

 // Kinematics
 AliGenPythiaEventHeader *pythiaHeader = GetPythiaHeader();
 AliGenHerwigEventHeader *HerwigHeader = GetHerwigHeader();
 if(pythiaHeader) {
   TArrayF t;
   pythiaHeader->PrimaryVertex(t);
   Float_t xvtx = t.GetAt(0);
   Float_t yvtx = t.GetAt(1);
   Float_t zvtx = t.GetAt(2);
   Float_t r2   = yvtx*yvtx+xvtx*xvtx;
   if(!(TMath::Abs(zvtx)<fVtxZMax&&r2<fVtxR2Max)) return;
   if((fPtHardAndPythiaJetPtFactor || fPtHardAndTrackPtFactor) && IsOutlier(GetPythiaHeader()))  return;
   fh1vZSelect->Fill(t.GetAt(2));
   fh1Xsec->Fill("<#sigma>",  pythiaHeader->GetXsection());
   fh1Trials->Fill("#sum{ntrials}",pythiaHeader->Trials());
  } else if(HerwigHeader) {
   TArrayF t;
   HerwigHeader->PrimaryVertex(t);
   Float_t xvtx = t.GetAt(0);
   Float_t yvtx = t.GetAt(1);
   Float_t zvtx = t.GetAt(2);
   Float_t r2   = yvtx*yvtx+xvtx*xvtx;
   if(!(TMath::Abs(zvtx)<fVtxZMax&&r2<fVtxR2Max)) return;
   fh1vZSelect->Fill(t.GetAt(2));
   fh1Xsec->Fill("<#sigma>",  HerwigHeader->Weight());
   fh1Trials->Fill("#sum{ntrials}",HerwigHeader->Trials());
  }


  else  fh1vZSelect->Fill(0);

}

  if (fDebug > 10)Printf("%s:%d vertex filled\n", (char*)__FILE__, __LINE__);

  //----------------  GET TRACKS  ---------------------------------------------------------------
  // we simply fetch the tracks/mc particles as a list of AliVParticles 
  // and it should be with typical parameters for track from Jet Branch Name, if it exsit

  vector<fastjet::PseudoJet> sorted_jets_fj[2];
  vector<fastjet::PseudoJet> full_event[2];
  //------------------------------------------------------------------------
  // we will have all the jets here(full_jets), either by algorithm, or we 
  // get them from AOD friend 
  //------------------------------------------------------------------------
  vector<fastjet::PseudoJet> full_jets;
  AliAODJet * tmp_jet;
  //------------------------------------------------------------------------
  // Setting for Jets:
  // create what we need for the clustering and the background estimation and subtraction
  //----------------------------------------------------------
  fastjet::JetDefinition jet_def(fAlgorithm, fRparam, fRecombScheme, fStrategy);
  fastjet::JetDefinition jet_def_for_rho(fBkgAlgorithm, fRparam,fRecombScheme, fStrategy);
  fastjet::GhostedAreaSpec ghostSpec(fGhostEtaMax, fActiveAreaRepeats, fGhostArea);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fAreaType,ghostSpec);
  // NB explicit ghosts do not work for moments with N < 0 with FastJet versions < 3.1  
  //  AreaDefinition area_def(fastjet::active_area_explicit_ghosts, 
  //                          fastjet::GhostedAreaSpec(fastjet::SelectorAbsRapMax(5.0)));
  //
  fastjet::Selector sel_jets = fastjet::SelectorNHardest(fNUsedJets) * fastjet::SelectorEtaRange(fJetEtaMin, fJetEtaMax);
 
  //---------------------- get all particles and jets ----------------------- 
  TList ParticleList[2];
  ParticleList[0].SetName("particleListGen");
  ParticleList[1].SetName("particleListRec");
  Int_t ifirstBr = 0;
  if((!fJetBranch[0].Contains("MC") && !fkDoJetReco) || (fkDoJetReco && !(fTrackType[0] == kTrackAODMCCharged || fTrackType[0] ==  kTrackAODMCChargedAcceptance  || fTrackType[0] == kTrackAODMCextra || fTrackType[0] == kTrackAODMCextraonly || fTrackType[0] == kTrackKineCharged || fTrackType[0] == kTrackKineChargedAcceptance || fTrackType[0] == kTrackKineAcceptance))) {
    ifirstBr = 1; // If no branch MC, loop only on data reco
  }
  if(!fJetBranch[1].Length() && !fkDoJetReco) ifirstBr = 2; // If no reco, do not do anything  
  for(Int_t iJetBranch = ifirstBr; iJetBranch < fgkFFMNJetBranches; iJetBranch++) {
    Int_t nParticles = GetListOfTracks(&ParticleList[iJetBranch], fTrackType[iJetBranch]);
    if(fDebug>2)Printf("%s:%d Number of selected tracks: %d", (char*)__FILE__, __LINE__, nParticles);
    if(!nParticles){ if(fDebug) Printf("%s:%d This event contains %d selected tracks ----> event skipped", (char*)__FILE__, __LINE__, nParticles); return;}
    if(TrackQAFilter(&ParticleList[iJetBranch], fTrackType[iJetBranch])) { if(fDebug>2) {Printf("%s:%d Event skipped in track QA filter!", (char*)__FILE__, __LINE__);}  return;}

    // Get list of particles in event
    for(Int_t i = 0; i < nParticles; i++) {
      AliVParticle *vp = (AliVParticle*)ParticleList[iJetBranch].At(i);
      Double_t fAllParticle[3] = {vp->Pt(), vp->Eta(), vp->Phi()};
      for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) { // 0: Pt, 1: Eta, 2: Phi
	if( iJetBranch == 0) fh1GenTracks[iAxis]->Fill(fAllParticle[iAxis]);
		else { //if( iJetBranch == 1) {
		  fh1RecTracks[iAxis]->Fill(fAllParticle[iAxis]);
		}
      }

      // Add particles to fastjet: used for jet finding and/or FFM
      fastjet::PseudoJet particle(vp->Px(), vp->Py(), vp->Pz(), vp->P());
      particle.set_user_index(i);
      full_event[iJetBranch].push_back(particle);

    } // End loop on list of particles
  } //End loop on list of JetBranch


  // Store array of association indexes between rec and gen tracks
  // arrays holding for each generated particle the reconstructed AOD track index & isPrimary flag, are initialized in AssociateGenRec(...) function
  TArrayI indexAODTr;
  TArrayS isGenPrim;
  // array holding for each reconstructed AOD track generated particle index, initialized in AssociateGenRec(...) function
  TArrayI indexMCTr;

  // ... and another set for secondaries from strange/non strange mothers (secondary MC tracks are stored in different lists)
  TArrayI indexAODTrSecNS;
  TArrayS isGenSecNS;
  TArrayI indexMCTrSecNS;

  TArrayI indexAODTrSecS;
  TArrayS isGenSecS;
  TArrayI indexMCTrSecS;

  Int_t  nTracksAODMCCharged = GetListOfTracks(fTracksAODMCCharged, kTrackAODMCCharged);
  if(fDebug>2)Printf("%s:%d selected AODMC tracks: %d ",(char*)__FILE__,__LINE__,nTracksAODMCCharged);

  Int_t  nTracksAODMCChargedSecNS = GetListOfTracks(fTracksAODMCChargedSecNS, kTrackAODMCChargedSecNS);
  if(fDebug>2)Printf("%s:%d selected AODMC secondary tracks NS: %d ",(char*)__FILE__,__LINE__,nTracksAODMCChargedSecNS);

  Int_t  nTracksAODMCChargedSecS = GetListOfTracks(fTracksAODMCChargedSecS, kTrackAODMCChargedSecS);
  if(fDebug>2)Printf("%s:%d selected AODMC secondary tracks S: %d ",(char*)__FILE__,__LINE__,nTracksAODMCChargedSecS);

  Int_t  nTracksRecQualityCuts = GetListOfTracks(fTracksRecQualityCuts, kTrackAODQualityCuts);
  if(fDebug>2)Printf("%s:%d selected rec tracks quality after cuts, full acceptance/pt : %d ",(char*)__FILE__,__LINE__,nTracksRecQualityCuts);

  // associate gen and rec tracks, store indices in TArrays
  AssociateGenRec(fTracksAODMCCharged, fTracksRecQualityCuts, indexAODTr, indexMCTr, isGenPrim,fh2PtRecVsGenPrim);
  AssociateGenRec(fTracksAODMCChargedSecNS,fTracksRecQualityCuts,indexAODTrSecNS,indexMCTrSecNS,isGenSecNS,fh2PtRecVsGenSec);
  AssociateGenRec(fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,indexMCTrSecS,isGenSecS,fh2PtRecVsGenSec);
  CalcSingleTrackEff(fh1TrackEffPtGen,fh2TrackEffEtaPhiGen,fh1TrackEffPtRec,fh2TrackEffEtaPhiRec,fTracksAODMCCharged,indexAODTr,isGenPrim);

 for(Int_t iJetBranch = ifirstBr; iJetBranch < fgkFFMNJetBranches; iJetBranch++) {
    
    if(fkDoJetReco) {
      // Jet finding switched on case
      if(fDebug > 3) printf("==== Jet finding (fastjet) is ON !!\n");
      // do the clustering
      fListJets[iJetBranch]->Clear();
      fastjet::ClusterSequenceArea clust_seq_full(full_event[iJetBranch], jet_def, area_def);
      full_jets = sel_jets(clust_seq_full.inclusive_jets());
      sorted_jets_fj[iJetBranch] = sorted_by_pt(full_jets);
      PseudoJetsToAODJets(sorted_jets_fj[iJetBranch], clust_seq_full, ParticleList[iJetBranch],fListJets[iJetBranch]);
    } else {
      // READ JETs from AOD (friend) trees of from output
      if(fDebug > 3) printf("==== Jets read from AOD branch !!\n");
      TClonesArray *aodJets[2] = {0x0};

      if (fkUseJetFromInput) {
	aodJets[iJetBranch] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranch[iJetBranch].Data()));
      }
      else { if(AODEvent()) aodJets[iJetBranch] = dynamic_cast<TClonesArray*>(AODEvent()->FindListObject(fJetBranch[iJetBranch].Data()));
      }
      
      if(!aodJets[iJetBranch]) {
        Printf("%s:%d no reconstructed Jet array with branch name %s found",(char*)__FILE__,__LINE__,fJetBranch[iJetBranch].Data());
        return; //continue; 
      }
      fListJets[iJetBranch]->Clear("nodelete");
      
      Int_t nJets = aodJets[iJetBranch]->GetEntries();
      if(fDebug) Printf("%s contains %d jet(s)", fJetBranch[iJetBranch].Data(), nJets);
      for(Int_t i=0; i<TMath::Min(nJets,fNUsedJets); i++){ 
        tmp_jet = dynamic_cast<AliAODJet*> ((*aodJets[iJetBranch])[i]);
		if ( SelectJet(tmp_jet, &(ParticleList[iJetBranch]))) {
          fListJets[iJetBranch]->Add(tmp_jet);
	  if(fDebug > 6) {tmp_jet->Print(Form("%d",i));}
	}
      }
    } // End else of fkDoJetReco

    if( fDebug>3 ) printf("%d jets selected after cuts in branch %d \n", fListJets[iJetBranch]->GetEntries(), iJetBranch);
    
    //Njets figs
    if(fListJets[iJetBranch]->GetEntries() > 0 ) fh1Njets[iJetBranch]->Fill(fListJets[iJetBranch]->GetEntries());

  } // End loop on iJetBranch 
  
  //---------------------------------------------------------------
  // Jet matching
  //---------------------------------------------------------------
  Int_t nGenJets = fListJets[0]->GetEntries();
  Int_t nRecJets = fListJets[1]->GetEntries();
  Int_t nUsedJets[fgkFFMNJetBranches] = {TMath::Min(fNUsedJets, nGenJets), TMath::Min(fNUsedJets, nRecJets)};
  TList* listUsedJets[fgkFFMNJetBranches];
  if(fkDoJetMatching){
    fListMatchedJets[0]->Clear("nodelete");
    fListMatchedJets[1]->Clear("nodelete");
    if(!fkUseClosestJetsforMatching) GetListOfMatchedJets(fListJets,nUsedJets,*&(fListMatchedJets[0]),*&(fListMatchedJets[1]),ifirstBr);
    else GetListOfClosestJets(fListJets,nUsedJets,*&(fListMatchedJets[0]),*&(fListMatchedJets[1]),ifirstBr);
    listUsedJets[0] = fListMatchedJets[0];
    listUsedJets[1] = fListMatchedJets[1];
  } // End kDoJetMatching
  else {
    listUsedJets[0] = fListJets[fkGenJetType];
    listUsedJets[1] = fListJets[1];
  }
  
  //--------------------------------------------------------------------
  // Do analysis  
  //--------------------------------------------------------------------
  static const Int_t nFFBins = fFFMMomMax;
  TArrayF *fFFMs_Raw[fgkFFMNJetBranches];
  TArrayF *fFFMs_Sub[fgkFFMNJetBranches];
  TArrayF *fFFMs_Imp[fgkFFMNJetBranches];

  Int_t usedDimTab[3] = {ifirstBr==1?listUsedJets[1]->GetEntries():fNUsedJets, ifirstBr==1?listUsedJets[1]->GetEntries():listUsedJets[0]->GetEntries(), listUsedJets[1]->GetEntries()};
  Int_t usedDim = TMath::MinElement(3,usedDimTab);
  for(Int_t iJetBranch = ifirstBr; iJetBranch < fgkFFMNJetBranches; iJetBranch++) {
    fFFMs_Raw[iJetBranch] = new TArrayF[nFFBins];
    fFFMs_Sub[iJetBranch] = new TArrayF[nFFBins];
    fFFMs_Imp[iJetBranch] = new TArrayF[nFFBins];
    // Initialize objects for FFM
    for(Int_t iN = 0; iN<fFFMMomMax; iN++) {
      fFFMs_Raw[iJetBranch][iN].Set(usedDim);
      fFFMs_Sub[iJetBranch][iN].Set(usedDim);
      fFFMs_Imp[iJetBranch][iN].Set(usedDim);
    }
    // Bkg tools for FFM
    fastjet::Selector rho_range = 0;
    if(fFFMBckgType==kUndef) {if(fDebug >9) Printf("FFM: Signal only mode");} // Signal only
    else if(fFFMBckgType==kDoughnut) rho_range = fastjet::SelectorDoughnut(fFFMBckgPar1,fFFMBckgPar2); // Rmin, Rmax
    else if(fFFMBckgType==kRectangle) rho_range = fastjet::SelectorRectangle(fFFMBckgPar1,fFFMBckgPar2); // HalfRapwidth, HalfPhiWidth
    else if(fFFMBckgType==kEtaRange) rho_range = fastjet::SelectorEtaRange( - fGhostEtaMax + fRparam,  fGhostEtaMax - fRparam); // EtaMin, EtaMax
    else if(fFFMBckgType==kRapPhiRange) rho_range = fastjet::SelectorRapPhiRange(- fGhostEtaMax + fRparam, fGhostEtaMax - fRparam,0, 2*TMath::Pi());
    else if(fFFMBckgType==kPerp) rho_range =  SelectorPerp(fRparam, TMath::Pi()/2.);
    else if(fFFMBckgType==kPerp2) rho_range =  SelectorPerp(fRparam, TMath::Pi()/2) || SelectorPerp(fRparam, -TMath::Pi()/2);
    else {Printf("Please check the background methods available to compute ffm !");}
    //if (fFFMBckgType==kPerp || fFFMBckgType==kPerp2) rho_range * fastjet::SelectorEtaRange(fJetEtaMin, fJetEtaMax);
    if (fFFMBckgType==kEtaRange || fFFMBckgType==kRapPhiRange) rho_range * !(fastjet::SelectorNHardest(2));
    fastjet::JetMedianBackgroundEstimator* bge = 0;
    fastjet::JetMedianBackgroundEstimator bg(rho_range, jet_def_for_rho, area_def);
    if(fFFBckgMode) {
    bge = &(bg); }
    fastjet::Subtractor subtractor(bge);
    if(fFFBckgMode) bge->set_particles(full_event[iJetBranch]);	// tell the background estimator what to use
    // Initialize JetFFMoments
    fastjet::contrib::JetFFMoments ffms_unsubtracted(fFFMNMin, fFFMNMax, nFFBins);
    fastjet::contrib::JetFFMoments ffms_subtracted(fFFMNMin, fFFMNMax, nFFBins, bge);
    fastjet::contrib::JetFFMoments ffms_improved(fFFMNMin, fFFMNMax, nFFBins, bge);
    double mu = fFFMBckgMu; // 25: typical value in the 150-200 GeV range
    if(fFFBckgMode) ffms_improved.set_improved_subtraction(mu, rho_range, full_event[iJetBranch], jet_def_for_rho, area_def);
    // With FastJet-3.1, this could be done this way:
    // ffms_improved.set_improved_subtraction(mu);

    // Clear pseudo jet list before to calculate FFM
    if(sorted_jets_fj[iJetBranch].size()) sorted_jets_fj[iJetBranch].clear();
    // jet track eff
    Double_t sumPtGenLeadingJetRecEff = 0;
    Double_t sumPtGenLeadingJetSec    = 0;
    Double_t sumPtRecLeadingJetRecEff = 0;
    Double_t sumPtGenLeadingJet = 0;
    Double_t sumPtRecLeadingJetReconstructed =0;
    for(Int_t ijet = 0; ijet < usedDim; ijet++ ) {
      AliAODJet* ujet = ((AliAODJet*)listUsedJets[iJetBranch]->At(ijet));
      if (!ujet) continue;
      if(!ClassifyJetEvent(ijet,usedDim)) break;

      fh2MatchedJetsAreaVSPt[iJetBranch]->Fill(ujet->EffectiveAreaCharged(), ujet->Pt());

       if(iJetBranch==1) {
        AliAODJet* uRecJet = (AliAODJet*)listUsedJets[iJetBranch]->At(ijet);
        AliAODJet* uGenJet = 0x0;
        AliAODJet* uEffJet = 0x0;
        if(listUsedJets[0]->GetEntries()) {uGenJet = (AliAODJet*)listUsedJets[0]->At(ijet);}
        else { uGenJet = (AliAODJet*)listUsedJets[iJetBranch]->At(ijet); }
        if(!fkEffJetType) {uEffJet = uGenJet;} else {uEffJet=uRecJet;}
      //*********************************************
      // Transverse and longitudinal jet info
      //*********************************************

          Bool_t isBadJetReconstructed     = kFALSE;
          TList* jettracklistReconstructed = new TList();
          jettracklistReconstructed->SetName("jettracklistReconstructed");
          GetTracksInJet(&ParticleList[1], jettracklistReconstructed, uRecJet, sumPtRecLeadingJetReconstructed, isBadJetReconstructed);

         if(GetJetMinNTracks()>0 && jettracklistReconstructed->GetSize()<=GetJetMinNTracks())     isBadJetReconstructed     = kTRUE;

         if(isBadJetReconstructed){
          delete jettracklistReconstructed;

          continue;
         }

        TArrayI inAOD; TArrayI inMC; TArrayS inTEST;
        CalcFFAndFFMFromTracksInJet(uRecJet, jettracklistReconstructed, 0, 0, inAOD, 0, 0 ,kFALSE, fh2TracksInJets[iJetBranch],fp2TracksInJetFFM[iJetBranch]);

        Double_t jetBkgPtRec = 0;
        Double_t jetBkgPtGen = 0;
         if(ijet == 0 && fFFBckgMode) {
           GetTracksTiltedwrpJetAxis(TMath::Pi()/2., &ParticleList[0],0x0,uGenJet,fRparam,jetBkgPtGen);
           GetTracksTiltedwrpJetAxis(TMath::Pi()/2., &ParticleList[1],0x0,uRecJet,fRparam,jetBkgPtRec);
           fh2MatchedJetsUE[0]->Fill(uGenJet->Pt(),jetBkgPtGen);
           fh2MatchedJetsUE[iJetBranch]->Fill(uRecJet->Pt(),jetBkgPtRec);
          }
      
      //*********************************************
      // MC <---> Reco Jet info
      //*********************************************
       Double_t jetEntriesMatch[10] = {
          uGenJet->Pt() - jetBkgPtGen , uGenJet->Eta(), uGenJet->Phi(), uGenJet->EffectiveAreaCharged(), (Double_t)( (TRefArray *) (uGenJet->GetRefTracks()) )->GetEntries(),
          uRecJet->Pt() - jetBkgPtRec , uRecJet->Eta(), uRecJet->Phi(), uRecJet->EffectiveAreaCharged(), (Double_t)( (TRefArray *) (uRecJet->GetRefTracks()) )->GetEntries(),
        };
        for(Int_t iAxis=0; iAxis<5; iAxis++) fh2MatchedJets[iAxis]->Fill(jetEntriesMatch[iAxis], jetEntriesMatch[iAxis+5]);

        if(listUsedJets[0]->GetEntries()) {
          //At least we need rec and gen, if no-matching, we MATCH the rec_jet and gen_jet with same number
          fh2MatchedJetsRDPtVSPt[0]->Fill((jetEntriesMatch[5]-jetEntriesMatch[0])/jetEntriesMatch[0], jetEntriesMatch[0]);
          fh2MatchedJetsRDPtVSPt[iJetBranch]->Fill((jetEntriesMatch[5]-jetEntriesMatch[0])/jetEntriesMatch[0], jetEntriesMatch[5]);

          Bool_t isBadJetGen     = kFALSE;

          TList* jettracklistGener = new TList();
          jettracklistGener->SetName("jettracklistGener");
          GetTracksInJet(&ParticleList[0], jettracklistGener, uGenJet, sumPtGenLeadingJet, isBadJetGen);

          if(GetJetMinNTracks()>0 && jettracklistGener->GetSize()<=GetJetMinNTracks())     isBadJetGen     = kTRUE;

          if(isBadJetGen){
          delete jettracklistGener;

          continue;
         }

           CalcFFAndFFMFromTracksInJet(uGenJet, jettracklistGener, 0, 0, inAOD, 0 ,0 ,kFALSE,  fh2TracksInJets[0],  fp2TracksInJetFFM[0]);

	  //efficiency:  primary tracks in jet
	  //Make list of constituents and check assosiation
	  Bool_t isBadJetGenPrim = kFALSE;
	  Bool_t isBadJetGenSec  = kFALSE;
	  Bool_t isBadJetRec     = kFALSE;

	  // for efficiency: gen tracks from pointing with gen/rec jet
	  TList* jettracklistGenPrim = new TList();
	  jettracklistGenPrim->SetName("jettracklistGenPrim");
	  GetTracksInJet(fTracksAODMCCharged, jettracklistGenPrim, uEffJet, sumPtGenLeadingJetRecEff, isBadJetGenPrim);

          TList* jettracklistGenSecNS = new TList();
	  jettracklistGenSecNS->SetName("jettracklistGenSecNS");
	  GetTracksInJet(fTracksAODMCChargedSecNS, jettracklistGenSecNS, uEffJet, sumPtGenLeadingJetSec,  isBadJetGenSec);

	  TList* jettracklistGenSecS = new TList();
	  jettracklistGenSecS->SetName("jettracklistGenSecS");
	  GetTracksInJet(fTracksAODMCChargedSecS, jettracklistGenSecS, uEffJet, sumPtGenLeadingJetSec, isBadJetGenSec);

	  // bin efficiency in jet pt bins using rec tracks
	  TList* jettracklistRec = new TList();
	  jettracklistRec->SetName("jettracklistRec");
	  GetTracksInJet(&ParticleList[1], jettracklistRec, uEffJet, sumPtRecLeadingJetRecEff, isBadJetRec);

	  if(GetJetMinNTracks()>0 && jettracklistGenPrim->GetSize()<=GetJetMinNTracks())   isBadJetGenPrim = kTRUE;
	  if(GetJetMinNTracks()>0 && jettracklistGenSecNS->GetSize()<=GetJetMinNTracks())  isBadJetGenSec  = kTRUE;
	  if(GetJetMinNTracks()>0 && jettracklistRec->GetSize()<=GetJetMinNTracks())       isBadJetRec     = kTRUE;

          if(isBadJetRec){
	    delete jettracklistGenPrim;
	    delete jettracklistGenSecNS;
	    delete jettracklistGenSecS;
	    delete jettracklistRec;

	    continue;
	  }

	    CalcFFAndFFMFromTracksInJet(uEffJet, jettracklistGenPrim, fTracksAODMCCharged, fTracksRecQualityCuts, indexAODTr, isGenPrim,jettracklistRec,kFALSE,fh2AssociatedTracksInJets[iJetBranch],fp2AssociatedTracksJetFFM[iJetBranch]);

	    // secondaries: use jet pt from primaries
	    CalcFFAndFFMFromTracksInJet(uEffJet,jettracklistGenSecNS,fTracksAODMCChargedSecNS,fTracksRecQualityCuts,indexAODTrSecNS,isGenSecNS,jettracklistRec,kFALSE,fh2AssociatedTracksInJetsSecNS[iJetBranch],fp2AssociatedTracksJetFFMSecNS[iJetBranch]);

	    CalcFFAndFFMFromTracksInJet(uEffJet,jettracklistGenSecS,fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,isGenSecS,jettracklistRec,kFALSE,fh2AssociatedTracksInJetsSecS[iJetBranch],fp2AssociatedTracksJetFFMSecS[iJetBranch]);

	    CalcFFAndFFMFromTracksInJet(uEffJet,jettracklistGenSecS,fTracksAODMCChargedSecS,fTracksRecQualityCuts,indexAODTrSecS,isGenSecS,jettracklistRec,kTRUE,fh2AssociatedTracksInJetsSecSsc[iJetBranch],fp2AssociatedTracksJetFFMSecSsc[iJetBranch]);

          delete jettracklistGener;
          delete jettracklistReconstructed;

	  delete jettracklistGenPrim;
	  delete jettracklistGenSecNS;
	  delete jettracklistGenSecS;
	  delete jettracklistRec;

	} // End 	if(listUsedJets[0]->GetEntries()
      }// End iJetBranch == 1

      //*******************************************/
      //  Dijet analysis
      //*********************************************
      if((usedDim>=2) && (ijet==1) && fAnaJetType.Contains("DIJET") ){
	Double_t leadingPhi =  ((AliAODJet*)(listUsedJets[iJetBranch]->At(0)))->Phi();
	Double_t associaPhi =  ((AliAODJet*)(listUsedJets[iJetBranch]->At(ijet)))->Phi();
	Double_t deltaPhi = TMath::Abs( TMath::Abs(associaPhi - leadingPhi));
	if( deltaPhi >= TMath::Pi()-fJetDeltaPhiCut && deltaPhi <= TMath::Pi()+fJetDeltaPhiCut ) {
	  Double_t leadingPt = ((AliAODJet*)(listUsedJets[iJetBranch]->At(0)))->Pt();
	  Double_t associaPt = ((AliAODJet*)(listUsedJets[iJetBranch]->At(ijet)))->Pt();
	  Double_t aj = (leadingPt - associaPt)/(leadingPt + associaPt);
	  fh1Asy_DiJets[iJetBranch]->Fill(aj);
	}
      }
      
       //*******************************************
       // Compute FFM
       //*******************************************
      fastjet::PseudoJet pseudoJetUsed;
      int ret = AliAODJetToPseudoJet(((AliAODJet*)listUsedJets[iJetBranch]->At(ijet)), pseudoJetUsed);
      if(ret) return;
      sorted_jets_fj[iJetBranch].push_back(pseudoJetUsed);
      if(fDebug > 9 && pseudoJetUsed.pt()<5) std::cout << "jet " << ijet << ", pt < 5: " << pseudoJetUsed.pt() << std::endl;

      if(fDebug > 9) { 
	fastjet::PseudoJet j = sorted_jets_fj[iJetBranch][ijet];
        vector<fastjet::PseudoJet> constituents = j.constituents();
        std::cout << "pseudo jet eta: " << j.eta() << ", phi: " << j.phi() << ", pt: " << j.pt() << ", area: " << j.area() << std::endl;
        fastjet::PseudoJet area4vect = j.area_4vector();
        std::cout << "a_x: " << area4vect.px() << ", a_y: " << area4vect.py() << "a_z: " << area4vect.pz() << std::endl;
        std::cout << "Number of constituents: " << constituents.size() << std::endl;
        for(UInt_t i=0; i<constituents.size(); i++){
          std::cout << "constituent " << i << ", pt: " << constituents[i].pt() << ", eta: " << constituents[i].eta() <<", phi: " << constituents[i].phi() << std::endl;      }

        if(fFFBckgMode) {
	fastjet::PseudoJet subtracted_jet = subtractor(j);
	std::cout <<"rho : " <<bge->rho(pseudoJetUsed) << " # [subtracted hard jet " << ijet+1 << "]:(pt,y,phi) = (" << subtracted_jet.pt() << ", "
	     << subtracted_jet.rap() << ", " << subtracted_jet.phi() << ")" << std::endl;
      }
     }

      // Compute FFM
      vector<double> ffm_unsubtracted = ffms_unsubtracted(pseudoJetUsed);
      vector<double> ffm_subtracted;
      vector<double> ffm_improved;

     if(fFFBckgMode) {
      ffm_subtracted = ffms_subtracted(pseudoJetUsed);
      ffm_improved = ffms_improved(pseudoJetUsed);
    } else {
      ffm_subtracted = ffms_unsubtracted(pseudoJetUsed);
      ffm_improved = ffms_unsubtracted(pseudoJetUsed);
   }
      //
      // Fill FFM histograms and ThnSparses
      //
	  Bool_t kWriteRaw = kTRUE;
	  Bool_t kWriteSub = kTRUE;
	  Bool_t kWriteImp = kTRUE;
	  for(Int_t iN = 0; iN < fFFMMomMax; iN++) {
		if(kWriteRaw && ffm_unsubtracted[iN] <= 0.) kWriteRaw = kFALSE;
		if(kWriteSub && ffm_unsubtracted[iN] <= 0.) kWriteSub = kFALSE;
		if(kWriteImp && ffm_unsubtracted[iN] <= 0.) kWriteImp = kFALSE;
	  }
	  if( !kWriteRaw || !kWriteSub || !kWriteImp ) printf("======================== one out of expected =====================\n");
      for(Int_t iN = 0; iN < fFFMMomMax; iN++) {
        Double_t A = TMath::Power( 0.5*(ffms_subtracted.N(iN) + 1), fFFMScalePower);
        fFFMs_Raw[iJetBranch][iN][ijet] = A*ffm_unsubtracted[iN];
	fFFMs_Sub[iJetBranch][iN][ijet] = A*ffm_subtracted[iN];
	fFFMs_Imp[iJetBranch][iN][ijet] = A*ffm_improved[iN];
		Double_t FFM_pt = pseudoJetUsed.pt();
        Double_t N = iN * (fFFMNMax - fFFMNMin)/(fFFMMomMax-1) + fFFMNMin;
		Double_t FFM_Raw_in_pt[3] = { fFFMs_Raw[iJetBranch][iN][ijet], FFM_pt, N };
		Double_t FFM_Sub_in_pt[3] = { fFFMs_Sub[iJetBranch][iN][ijet], FFM_pt, N };
		Double_t FFM_Imp_in_pt[3] = { fFFMs_Imp[iJetBranch][iN][ijet], FFM_pt, N };
		if(kWriteRaw) { fhnJetFFM_Raw[iJetBranch]->Fill( FFM_Raw_in_pt );
						if(FFM_pt>=5) fp2JetFFM_Raw[iJetBranch]->Fill(N, FFM_pt, fFFMs_Raw[iJetBranch][iN][ijet]);}
		if(kWriteSub) { fhnJetFFM_Sub[iJetBranch]->Fill( FFM_Sub_in_pt );
						fp2JetFFM_Sub[iJetBranch]->Fill(N, FFM_pt, fFFMs_Sub[iJetBranch][iN][ijet]);}
		if(kWriteImp) { fhnJetFFM_Imp[iJetBranch]->Fill( FFM_Imp_in_pt );
						fp2JetFFM_Imp[iJetBranch]->Fill(N, FFM_pt, fFFMs_Imp[iJetBranch][iN][ijet]);}

        if( fHistosLevel >=9 ) {
          if( iJetBranch == 1 ) {
      	    fastjet::PseudoJet ugenjet;
      	    Float_t usedFFMs_Raw_0 = 0;
      	    Float_t usedFFMs_Sub_0 = 0;
      	    Float_t usedFFMs_Imp_0 = 0;
      	    if(sorted_jets_fj[0].size()) { ugenjet = sorted_jets_fj[0][ijet]; usedFFMs_Raw_0 = fFFMs_Raw[0][iN][ijet];
	      usedFFMs_Sub_0 = fFFMs_Sub[0][iN][ijet]; usedFFMs_Imp_0 = fFFMs_Imp[0][iN][ijet]; }
      	    else { ugenjet = sorted_jets_fj[iJetBranch][ijet]; usedFFMs_Raw_0 = fFFMs_Raw[1][iN][ijet];
	      usedFFMs_Sub_0 = fFFMs_Sub[1][iN][ijet]; usedFFMs_Imp_0 = fFFMs_Imp[1][iN][ijet];  }
      	 
            Float_t usedPtFraction = ((AliAODJet*)listUsedJets[1]->At(ijet))->GetNEF();
            //for Pb-Pb, we may need centrality
            //(Double_t)cenClass, (Double_t)nInputTracks, (Double_t)rpBin,
      	    Double_t jetFFM_Raw[5] = {
      	      usedFFMs_Raw_0, fFFMs_Raw[iJetBranch][iN][ijet], ugenjet.pt(), pseudoJetUsed.pt(), usedPtFraction
      	    }; // NB: here NEF = PtFraction
	    Double_t jetFFM_Sub[5] = {
	      usedFFMs_Sub_0, fFFMs_Sub[iJetBranch][iN][ijet], ugenjet.pt(), pseudoJetUsed.pt(), usedPtFraction
	    };
	    
	    Double_t jetFFM_Imp[5] = {
	      usedFFMs_Imp_0, fFFMs_Imp[iJetBranch][iN][ijet], ugenjet.pt(), pseudoJetUsed.pt(), usedPtFraction
	    };
	    
      	    fhnJetMomN_Raw[iN]->Fill(jetFFM_Raw);
	    fhnJetMomN_Sub[iN]->Fill(jetFFM_Sub);
	    fhnJetMomN_Imp[iN]->Fill(jetFFM_Imp);

       	  } // End ijetbranch = 1
       	} // End histo level
      } // End loop on iN
      // //////////////////////////////////////// FFM done!!!!! ////////////////////////////////////////

      if(fDebug > 8) {
        std::cout << "# Fragmentation function moments for full jet "
             << ": (pt, eta, phi) = (" << pseudoJetUsed.pt() << ", "<< pseudoJetUsed.eta() << ", " << pseudoJetUsed.phi() << "), #constituents="<< pseudoJetUsed.constituents().size() << std::endl
             <<": (px, py, pz, E )  = ("<<pseudoJetUsed.px()<<", "<<pseudoJetUsed.py()<<", "<<pseudoJetUsed.pz()<<", "<<pseudoJetUsed.e()<<")"<<"  sqrt(px^2+py^2): "<<TMath::Sqrt(pseudoJetUsed.px()*pseudoJetUsed.px() + pseudoJetUsed.py()*pseudoJetUsed.py())<<std::endl;
        std::cout << "# N	M_N(unsubtracted)	M_N(subtracted)	      M_N(improved)" << std::endl;
        for (unsigned int in=0; in<ffm_subtracted.size(); in++){
	  printf("%1.2f	   %2.9f  	  %2.9f  	  %2.9f	 \n", ffms_subtracted.N(in), ffm_unsubtracted[in], ffm_subtracted[in], ffm_improved[in]);
        }
        std::cout << "# N	A*M_N(unsubtracted)	A*M_N(subtracted)     A*M_N(improved)" << std::endl;
        for (unsigned int in=0; in<ffm_subtracted.size(); in++){
	  printf("%1.2f	   %2.9f  	  %2.9f  	  %2.9f  \n", ffms_subtracted.N(in), fFFMs_Raw[iJetBranch][in][ijet], fFFMs_Sub[iJetBranch][in][ijet], fFFMs_Imp[iJetBranch][in][ijet]);
        }
      } // End fDebug
    } // End loop on ijet
    
  } // End loop on jet branches

   fTracksAODMCCharged->Clear();
   fTracksAODMCChargedSecNS->Clear();
   fTracksAODMCChargedSecS->Clear();
   fTracksRecQualityCuts->Clear();

  // Store the jet branch in the AOD
  static AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  if(fTCAJetsOut&&aodH&&fFilterPt>0){
    if(fTCAJetsOut->GetEntries()>0){
      AliAODJet *jet = (AliAODJet*)fTCAJetsOut->At(0);
      if(jet->Pt()>fFilterPt){
        aodH->SetFillAOD(kTRUE);
      }
    }
  }

  for(Int_t iJetBranch = ifirstBr; iJetBranch < fgkFFMNJetBranches; iJetBranch++) {
    delete [] fFFMs_Raw[iJetBranch];
    delete [] fFFMs_Sub[iJetBranch];
    delete [] fFFMs_Imp[iJetBranch];
  }

  PostData(1, fHistList);

}

// ---------------------------------------------------------------
void AliAnalysisTaskJetFFMoments::Terminate(Option_t */*option*/)
{
  //
  // Terminate analysis
  //
  if (fDebug !=0) printf("AnalysisJetFFMoments: Terminate() \n");

}

// ---------------------------------------------------------------
void AliAnalysisTaskJetFFMoments::GetTracksInJet(TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t& sumPt, Bool_t& isBadPt)
 {
   if(fFFBckgMode) {outputlist->AddAll(inputlist); return;}
   if(fTracksInJetMethod==0){ // Pointing for MC and Rec
     GetJetTracksPointing(inputlist, outputlist, jet, fRparam, sumPt, GetJetMinLTrackPt() , GetJetMaxTrackPt(), isBadPt);
   } else if(fTracksInJetMethod==1) { // Pointing for MC and refTracks for Rec
     if((TString(inputlist->GetName())).Contains("MC")){ // MC Case
       GetJetTracksPointing(inputlist, outputlist, jet, TMath::Abs(fRparam)+0.2, sumPt, GetJetMinLTrackPt() , GetJetMaxTrackPt(), isBadPt);
     }
     else { // Rec case
       GetJetTracksTrackrefs(outputlist, jet, GetJetMinLTrackPt() , GetJetMaxTrackPt(), isBadPt);
     }
   } else if(fTracksInJetMethod==2){ // refTracks in both MC and rec
     GetJetTracksTrackrefs(outputlist, jet, GetJetMinLTrackPt() , GetJetMaxTrackPt(), isBadPt);
   } else if(fTracksInJetMethod==3){ // Mengliang's method : simple refTracks
     TRefArray* raTrackRef = jet->GetRefTracks();
     for( Int_t icons = 0; icons<raTrackRef->GetEntries(); icons++){
     AliVParticle* track = dynamic_cast<AliVParticle*>(raTrackRef->At(icons));
     if(!track)continue;
     outputlist->Add(track);
     }
   } else  {printf("Unkown track in jet method %d",fTracksInJetMethod);}

 }

// ---------------------------------------------------------------
Int_t  AliAnalysisTaskJetFFMoments::GetListOfTracks(TList *list,Int_t type)
{

  //
  // get list of tracks/particles for different types
  //

  if(fDebug>2) Printf("%s:%d Selecting tracks with %d",(char*)__FILE__,__LINE__,type);

  Int_t iCount = 0;
  if(type==kTrackAOD || type==kTrackAODQualityCuts || type==kTrackAODextra || type==kTrackAODextraonly || type==kTrackAODMCextra || type==kTrackAODMCextraonly){

    if(type!=kTrackAODextraonly && type!=kTrackAODMCextraonly) {
      AliAODEvent *aod = 0;
      if(fkUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
      else aod = AODEvent();
      if(!aod){
	if(fDebug>2)Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
	return iCount;
      }

      for(int it = 0;it < aod->GetNumberOfTracks();++it){
	AliAODTrack *tr = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
        if(!tr) AliFatal("Not a standard AOD");
	Bool_t bGood = false;
	if(fFilterType == 0)bGood = true;
	else if(fFilterType == 1)bGood = tr->IsHybridTPCConstrainedGlobal();
	else if(fFilterType == 2)bGood = tr->IsHybridGlobalConstrainedGlobal();
	if((fFilterMask>0)&&((!tr->TestFilterBit(fFilterMask)||(!bGood)))){
	  if(fDebug>10)Printf("%s:%d Not matching filter %d/%d %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks(),fFilterMask,tr->GetFilterMap());	
	  continue;
	}

	// heavy flavor jets
	if(fFilterMask==528 && fkUseHFcuts){
          Double_t ntpcClus = tr->GetTPCNcls();
          Double_t trPt=tr->Pt();
	  TFormula NTPCClsCut("f1NClustersTPCLinearPtDep","70.+30./20.*x");
	
	  if (trPt <= 20. && (ntpcClus < NTPCClsCut.Eval(trPt))) continue;
	  else if (trPt > 20. && ntpcClus < 100) continue;

	  if (AvoidDoubleCountingHF(aod,tr)) continue;
	}
	//

        if(fkRequireITSRefit){if((tr->GetStatus()&AliESDtrack::kITSrefit)==0)continue;}
        if (fkApplySharedClusterCut) {
	  Double_t frac = Double_t(tr->GetTPCnclsS()) /Double_t(tr->GetTPCncls());
	  if (frac > 0.4) continue;
        } 

		if(type!=kTrackAODQualityCuts){
		  if( tr->Eta()>fTrackEtaMax || tr->Eta()<fTrackEtaMin){
	  if(fDebug>10)Printf("%s:%d Not matching eta %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	  continue;
	}
		  if( tr->Phi()>fTrackPhiMax || tr->Phi()<fTrackPhiMin){
		    if(fDebug>10)Printf("%s:%d Not matching phi %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());
		    continue;
		  }
	if(tr->Pt()<fTrackPtMin){
	  if(fDebug>10)Printf("%s:%d Not matching pt %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	  continue;
	}
		}
	if(fDebug>10)Printf("%s:%d MATCHED %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	list->Add(tr);
	iCount++;
      }
    }
    if(type==kTrackAODextra || type==kTrackAODextraonly) {
      AliAODEvent *aod = 0;
      if(fkUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
      else aod = AODEvent();
      
      if(!aod){
	return iCount;
      }
      TClonesArray *aodExtraTracks = dynamic_cast<TClonesArray*>(aod->FindListObject("aodExtraTracks"));
      if(!aodExtraTracks)return iCount;
      for(int it =0; it<aodExtraTracks->GetEntries(); it++) {
	AliVParticle *track = dynamic_cast<AliVParticle*> ((*aodExtraTracks)[it]);
	if (!track) continue;

	AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*> (track);
	if(!trackAOD)continue;
	Bool_t bGood = false;
	if(fFilterType == 0)bGood = true;
	else if(fFilterType == 1)bGood = trackAOD->IsHybridTPCConstrainedGlobal();
	else if(fFilterType == 2)bGood = trackAOD->IsHybridGlobalConstrainedGlobal();
	if((fFilterMask>0)&&((!trackAOD->TestFilterBit(fFilterMask)||(!bGood))))continue;
        if(fkRequireITSRefit){if((trackAOD->GetStatus()&AliESDtrack::kITSrefit)==0)continue;}
	if (fkApplySharedClusterCut) {
	  Double_t frac = Double_t(trackAOD->GetTPCnclsS()) /Double_t(trackAOD->GetTPCncls());
	  if (frac > 0.4) continue;
	}


		if( trackAOD->Eta()>fTrackEtaMax || trackAOD->Eta()<fTrackEtaMin) continue;
	if(trackAOD->Pt()<fTrackPtMin) continue;
	if(fDebug) printf("pt extra track %.2f \n", trackAOD->Pt());
	list->Add(trackAOD);
	iCount++;
      }
    }

    if(type==kTrackAODMCextra || type==kTrackAODMCextraonly) { //embed generator level particles
      AliAODEvent *aod = 0;
      if(fkUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
      else aod = AODEvent();
      if(!aod){
	return iCount;
      }
      TClonesArray *aodExtraTracks = dynamic_cast<TClonesArray*>(aod->FindListObject("aodExtraMCparticles"));
      if(!aodExtraTracks)return iCount;
      for(int it =0; it<aodExtraTracks->GetEntries(); it++) {
	AliVParticle *track = dynamic_cast<AliVParticle*> ((*aodExtraTracks)[it]);
	AliAODMCParticle *partmc = dynamic_cast<AliAODMCParticle*> ((*aodExtraTracks)[it]);
	if (!track) {
	  if(fDebug>10)  printf("track %d does not exist\n",it);
	  continue;
	}

	if(!partmc) continue;
	if(!partmc->IsPhysicalPrimary())continue;

	if (track->Pt()<fTrackPtMin) {
	  if(fDebug>10)  printf("track %d has too low pt %.2f\n",it,track->Pt());
	  continue;
	}

	if(fDebug>10) printf("pt extra track %.2f \n", track->Pt());        
	
		if( track->Eta()>fTrackEtaMax || track->Eta()<fTrackEtaMin) continue;
	if(track->Pt()<fTrackPtMin) continue;
	list->Add(track);

	iCount++;
      }
    }
    
  }
  else if (type ==  kTrackKineAll||type == kTrackKineCharged || type == kTrackKineChargedAcceptance  || type == kTrackKineChargedAcceptanceDet){
    AliMCEvent* mcEvent = MCEvent();
    if(!mcEvent)return iCount;
    // we want to have alivpartilces so use get track
    for(int it = 0;it < mcEvent->GetNumberOfTracks();++it){
      if(!mcEvent->IsPhysicalPrimary(it))continue;
      AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(it);
      if(type == kTrackKineAll){
	if(part->Pt()<fTrackPtMin)continue;
	list->Add(part);
	iCount++;
      }

      else if(type == kTrackKineCharged || type == kTrackKineChargedAcceptance  || type == kTrackKineChargedAcceptanceDet || type == kTrackKineAcceptance ||  type == kTrackKineAcceptanceDet){
        if (type == kTrackKineCharged || type == kTrackKineChargedAcceptance  || type == kTrackKineChargedAcceptanceDet) {
         if(part->Particle()->GetPDG()->Charge()==0)continue;
        }

        if((type == kTrackKineChargedAcceptance || type == kTrackKineChargedAcceptanceDet || type == kTrackKineAcceptance || type == kTrackKineAcceptanceDet) &&
           (       part->Eta() < fTrackEtaMin
                || part->Eta() > fTrackEtaMax
                || part->Phi() < fTrackPhiMin
                || part->Phi() > fTrackPhiMax
                || part->Pt()  < fTrackPtMin)) continue;

        if(type == kTrackKineChargedAcceptanceDet || type == kTrackKineAcceptanceDet)
          {

         fh1TrackEffPtGen->Fill(part->Pt());
         fh2TrackEffEtaPhiGen->Fill(part->Eta(),TVector2::Phi_0_2pi(part->Phi()));

         if (fResol) {
            TLorentzVector newpart;
           if(fResolMeth == 0) {

            Float_t p  = part->P();
            Float_t p0 = p + gRandom->Gaus(0., fResol*(1+fResolvar/100.));
            // Cinematic recomputed taking into account potential track momentum resolution
            Float_t r  = p/p0;
            Float_t px = r * part->Px();
            Float_t py = r * part->Py();
            Float_t pz = r * part->Pz();
            Float_t m = part->M();
            Float_t e  = TMath::Sqrt(px * px + py * py + pz * pz + m * m);
            newpart = TLorentzVector(px, py, pz, e);

         } else if(fResolMeth == 1) {

             Float_t r = gRandom->Gaus(1., fResol*(1+fResolvar/100.));
             // Cinematic recomputed taking into account potential track momentum resolution
             Float_t px = r * part->Px();
             Float_t py = r * part->Py();
             Float_t pz = r * part->Pz();
             Float_t m = part->M();
             Float_t e  = TMath::Sqrt(px * px + py * py + pz * pz + m * m);
             newpart = TLorentzVector(px, py, pz, e);

        } else if(fResolMeth == 2) {

            Float_t px = part->Px();
            Float_t py = part->Py();
            Float_t pz = part->Pz();
            Float_t e  = part->E();
            Float_t pt = part->Pt();
            newpart = TLorentzVector(px, py, pz,e);
            Float_t ptSmeared= pt * gRandom->Gaus(1.0,fResol*(1+fResolvar/100.));
            newpart.SetPerp(ptSmeared);

       } else if(fResolMeth == 3) {

            Float_t px = part->Px();
            Float_t py = part->Py();
            Float_t pz = part->Pz();
            Float_t e  = part->E();
            Float_t pt = part->Pt();
            newpart = TLorentzVector(px, py, pz,e);
            Float_t ptSmeared= 1./pt * gRandom->Gaus(1.0,fResol*(1+fResolvar/100.));
            newpart.SetPerp(1./ptSmeared);
       }
            fh2TrackResPt->Fill(part->Pt(),100*(TMath::Abs(part->Pt()-newpart.Pt())/part->Pt()));
            fh1TrackResPtInv->Fill(100. * (TMath::Abs(1/newpart.Pt()) - 1./part->Pt()) * part->Pt());
            part->Particle()->SetMomentum(newpart);
           }

          if (fEffi &&  fEffi->InheritsFrom("TF1") &&  (gRandom->Rndm()  > (1+fEffivar/100.)*(((TF1*) fEffi)->Eval(part->Pt())))) continue;
          if (fEffi &&  fEffi->InheritsFrom("TH1") &&  (gRandom->Rndm()  > (1+fEffivar/100.)*(((TH1*) fEffi)->GetBinContent(((TH1*)fEffi)->FindBin(part->Pt()))))) continue;
          if (fEffi &&  fEffi->InheritsFrom("TF1") &&  (((1+fEffivar/100.)*(((TF1*) fEffi)->Eval(part->Pt()))) > 1) ) continue;
          if (fEffi &&  fEffi->InheritsFrom("TH1") &&  (((1+fEffivar/100.)*(((TH1*) fEffi)->GetBinContent(((TH1*)fEffi)->FindBin(part->Pt())))) > 1)) continue;

        
        if (       part->Eta() < fTrackEtaMin
                || part->Eta() > fTrackEtaMax
                || part->Phi() < fTrackPhiMin
                || part->Phi() > fTrackPhiMax
                || part->Pt()  < fTrackPtMin) continue;


         fh1TrackEffPtRec->Fill(part->Pt());
         fh2TrackEffEtaPhiRec->Fill(part->Eta(),TVector2::Phi_0_2pi(part->Phi()));


         }

	list->Add(part);
	iCount++;
      }
    }
  }
  else if (type == kTrackAODMCCharged || type == kTrackAODMCAll || type == kTrackAODMCChargedAcceptance || type==kTrackAODMCChargedSecNS || type==kTrackAODMCChargedSecS) {
    AliAODEvent *aod = 0;
    if(fkUseAODMCInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else aod = AODEvent();
    if(!aod)return iCount;
    TClonesArray *tca = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    for(int it = 0;it < tca->GetEntriesFast();++it){
      AliAODMCParticle *part = (AliAODMCParticle*)(tca->At(it));
	  if(!part)continue;
	  if(type != kTrackAODMCChargedSecNS && type != kTrackAODMCChargedSecS  && !part->IsPhysicalPrimary())continue;
          if((type == kTrackAODMCChargedSecNS || type == kTrackAODMCChargedSecS) && part->IsPhysicalPrimary())continue;

      if(type == kTrackAODMCAll){
	if(part->Pt()<fTrackPtMin)continue;
	list->Add(part);
	iCount++;
      }
	  else if (type == kTrackAODMCCharged || type == kTrackAODMCChargedAcceptance || type==kTrackAODMCChargedSecNS || type==kTrackAODMCChargedSecS){
	if(part->Charge()==0)continue;

		if(type==kTrackAODMCChargedSecNS || type==kTrackAODMCChargedSecS){
		  Bool_t isFromStrange = kFALSE;
		  Int_t iMother = part->GetMother();

		  if(iMother < 0) continue; // throw out PYTHIA stack partons + incoming protons

		  AliAODMCParticle *partM = dynamic_cast<AliAODMCParticle*>(tca->At(iMother));
		  if(!partM) continue;

		  Int_t codeM =  TMath::Abs(partM->GetPdgCode());
		  Int_t mfl = Int_t (codeM/ TMath::Power(10, Int_t(TMath::Log10(codeM))));
		  if  (mfl == 3 && codeM != 3) isFromStrange = kTRUE;

		  if(codeM == 130) isFromStrange = kTRUE; // K0 long
		  if(part->IsSecondaryFromMaterial()) isFromStrange = kFALSE; // strange resonances from hadronic showers ?

		  if(type==kTrackAODMCChargedSecNS && isFromStrange) continue;
		  if(type==kTrackAODMCChargedSecS  && !isFromStrange) continue;

	  list->Add(part);
		  iCount++;
	}

		else if(type == kTrackAODMCCharged){
	  list->Add(part);
		  iCount++;
	}
		else if(type == kTrackAODMCChargedAcceptance){
		  if( part->Eta()>fTrackEtaMax || part->Eta()<fTrackEtaMin) continue;
		  if( part->Phi()>fTrackPhiMax || part->Phi()<fTrackPhiMin) continue;
		  if( part->Pt()<fTrackPtMin) continue;
		  list->Add(part);
	iCount++;
      }
    }
	}
  }// AODMCparticle
  else if (type == kTrackAODMCHF){
	  
    AliAODEvent *aod = 0;
    if(fkUseAODMCInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else aod = AODEvent();  
    if(!aod)return iCount;
    TClonesArray *tca = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    for(int it = 0;it < tca->GetEntriesFast();++it){
      AliAODMCParticle *part = (AliAODMCParticle*)(tca->At(it));
      if(!part) continue;
      int partpdg= part->PdgCode();
      if(!part->IsPhysicalPrimary() && !IsBMeson(partpdg) && !IsDMeson(partpdg) )continue;
		
      if (IsBMeson(partpdg) || IsDMeson(partpdg)) {
	iCount+= AddDaughters( list , part,tca);
      }
      else {
		if((part->Pt()>=fTrackPtMin) && (part->Eta()<=fTrackEtaMax) && (part->Eta()>=fTrackEtaMin) && (part->Charge()!=0))list->Add(part);
	iCount++;
      }
    }
  }
  
  list->Sort();
  return iCount;
}

// ---------------------------------------------------------------------------------
Int_t AliAnalysisTaskJetFFMoments::AddDaughters(TList * list, AliAODMCParticle *part, TClonesArray * tca){
	Int_t count=0;
	Int_t nDaugthers = part->GetNDaughters();
	for(Int_t i=0;i< nDaugthers;++i){
		AliAODMCParticle *partdaughter = (AliAODMCParticle*)(tca->At(i));
			if(!partdaughter) continue;
			if(partdaughter->Pt()<fTrackPtMin)continue;
			if(TMath::Abs(partdaughter->Eta())>fTrackEtaMax)continue;
			if(partdaughter->Charge()==0)continue;

	if(!IsDMeson(partdaughter->PdgCode()) && !IsBMeson(partdaughter->PdgCode()) ){
		list->Add(partdaughter);
		count++;
		}
		else count+=AddDaughters(list,part,tca);
	}
return count;
}

// ---------------------------------------------------------------------------------
bool AliAnalysisTaskJetFFMoments::IsBMeson(int pc)
{
  int bPdG[] = {511,521,10511,10521,513,523,10513,10523,20513,20523,20513,20523,515,525,531,
		10531,533,10533,20533,535,541,10541,543,10543,20543,545,551,10551,100551,
		110551,200551,210551,553,10553,20553,30553,100553,110553,120553,130553,200553,210553,220553,
		300553,9000533,9010553,555,10555,20555,100555,110555,120555,200555,557,100557};
  for(int i=0;i< (int)(sizeof(bPdG)/sizeof(int));++i) if(abs(pc) == bPdG[i]) return true;
  return false;
}

// ---------------------------------------------------------------------------------
bool AliAnalysisTaskJetFFMoments::IsDMeson(int pc)
{
  int bPdG[] = {411,421,10411,10421,413,423,10413,10423,20431,20423,415,
		425,431,10431,433,10433,20433,435,441,10441,100441,443,10443,20443,
		100443,30443,9000443,9010443,9020443,445,100445};
  for(int i=0;i< (int)(sizeof(bPdG)/sizeof(int));++i) if(abs(pc) == bPdG[i]) return true;
  return false;
}

// ---------------------------------------------------------------------------------
Bool_t AliAnalysisTaskJetFFMoments::AvoidDoubleCountingHF(AliAODEvent *aod, AliAODTrack *tr1){

  if(!(tr1->TestFilterBit(BIT(9)))) return kFALSE;

  Int_t idtr1 = tr1->GetID();

  for(int jt = 0;jt < aod->GetNumberOfTracks();++jt){

    AliAODTrack *tr2 = dynamic_cast<AliAODTrack*>(aod->GetTrack(jt));
    if(!tr2) AliFatal("Not a standard AOD");
    Int_t idtr2 = tr2->GetID();

    if (!(tr2->TestFilterBit(BIT(4)))) continue;
    if (idtr1==(idtr2+1)*-1.) return kTRUE;

  }
  return kFALSE;
}

// _______________________________________________________________________________
Int_t AliAnalysisTaskJetFFMoments::GetListOfMatchedJets(TList* listJets[fgkFFMNJetBranches], Int_t nUsedJets[fgkFFMNJetBranches], TList* list1, TList* list2, Int_t ifirstBr)
{

  TList* listMatchedJets[fgkFFMNJetBranches];
  listMatchedJets[0] = list1;
  listMatchedJets[1] = list2;

  TArrayI aMatchIndex;
  aMatchIndex.Set(nUsedJets[0]);
  TArrayF aPtFraction;
  aPtFraction.Set(nUsedJets[0]);
  // See comments in AliAnalysisHeplerTasks: 
  // aMatchIndex: index of the reconstructed jet which matches the generated one. If -1, no matching 
  // aPtFraction: fraction of the generated pt that ends up in the reconstructed jet  
  AliAnalysisHelperJetTasks::GetJetMatching(listJets[0], nUsedJets[0], listJets[1], nUsedJets[1],
					    aMatchIndex, aPtFraction, fDebug, fJetMatchedDistMax, fkIsPbPb?1:2);

  // Fill the mismatched histograms and remove the unmatched jets from the jet list
  for( Int_t iJetBranch = ifirstBr; iJetBranch < fgkFFMNJetBranches; iJetBranch++) {
    for( Int_t ijet = 0; ijet < nUsedJets[iJetBranch]; ijet++ ) {
      AliAODJet* tmpjet = 0x0;
      if(fDebug>2 && iJetBranch==0) printf("gen branch jets: %2d matched with rec jet: %2d with a fraction %2.4f\n", ijet, aMatchIndex[ijet], aPtFraction[ijet]); 

      // Fill the lists of matched generated and reconstructed jets
      if((iJetBranch==0)&&(aMatchIndex[ijet]!=-1)&&(aPtFraction[ijet]>=fJetMatchingFractionMin) ){
	listMatchedJets[iJetBranch]->Add( ((AliAODJet*) listJets[iJetBranch]->At(ijet)));
	((AliAODJet*)listJets[1]->At(aMatchIndex[ijet]))->SetNEF(aPtFraction[ijet]); // NB: here NEF = PtFraction
	listMatchedJets[1]->Add(((AliAODJet*) listJets[1]->At(aMatchIndex[ijet])));
      }

      if(fkFillMismatchHisto) { 
	// Fetch the mismatched generated jets
	if( ( (iJetBranch==0) && (aMatchIndex[ijet]==-1)) || 
	    ( (iJetBranch==0) && (aPtFraction[ijet] < fJetMatchingFractionMin) ) ) { 
	  tmpjet = ((AliAODJet*) listJets[iJetBranch]->At(ijet));
	}
	  
	// Fetch the mismatched reconstructed jets
	if(iJetBranch==1){
	  if(!nUsedJets[0]) { tmpjet = ((AliAODJet*) listJets[iJetBranch]->At(ijet));}
	  for(Int_t i=0; i< nUsedJets[0]; i++){
	    // Compare id of reco jet to id of match id of gen jets
	    if( (ijet==aMatchIndex[i]) && (aPtFraction[i] >= fJetMatchingFractionMin)) {break;} 
	    if(i==nUsedJets[0]-1){
	      tmpjet = ((AliAODJet*) listJets[iJetBranch]->At(ijet));
	    }
	  } // End loop on generated jets
	} // End iJetBranch == 1

	  // Mismatched histogram filling
	if(tmpjet){
	  Double_t jetEntriesMismatch[5] = {tmpjet->Pt(), tmpjet->Eta(), tmpjet->Phi(), tmpjet->EffectiveAreaCharged(), (Double_t)(( (TRefArray *) (tmpjet->GetRefTracks()))->GetEntries())};
	  for( Int_t iAxis = 0; iAxis < 5; iAxis++) fh1JetPr_Mismatched[iJetBranch][iAxis]->Fill(jetEntriesMismatch[iAxis]);
		  fh2MismatchedJetsAreaVSPt[iJetBranch]->Fill(tmpjet->EffectiveAreaCharged(), tmpjet->Pt());
	  // After keeping the mismatched information, remove the unmatched jets from the jet list
	} // End tmpjet
      } // End Fill mismatched histo
    } // End loop on jets
  } // End loop on jet branches

  return listJets[1]->GetEntries();

}
//_______________________________________________________________________________
Int_t AliAnalysisTaskJetFFMoments::GetListOfClosestJets(TList* listJets[fgkFFMNJetBranches], Int_t nUsedJets[fgkFFMNJetBranches], TList* list1, TList* list2, Int_t ifirstBr)
{

  TList* listMatchedJets[fgkFFMNJetBranches];
  listMatchedJets[0] = list1; //MC
  listMatchedJets[1] = list2;

  TArrayI aGenIndex;
  aGenIndex.Set(TMath::Max(nUsedJets[0],nUsedJets[1]));
  TArrayI aRecIndex;
  aRecIndex.Set(TMath::Max(nUsedJets[0],nUsedJets[1]));
  // See comments in AliAnalysisHeplerTasks:
  // aRecIndex: index of the reconstructed jet which matches the generated one. If -1, no matching
  // aGenIndex: index of the generated jet which matches the reconstructed one. If -1, no matching
  AliAnalysisHelperJetTasks::GetClosestJets(listJets[0], nUsedJets[0], listJets[1], nUsedJets[1],
                                            aGenIndex, aRecIndex, fDebug, fJetMatchedDistMax);




  // Fill the mismatched histograms and remove the unmatched jets from the jet list
  for( Int_t iJetBranch = ifirstBr; iJetBranch < fgkFFMNJetBranches; iJetBranch++) {
    for( Int_t ijet = 0; ijet < nUsedJets[iJetBranch]; ijet++ ) {
      if(fDebug>2 && iJetBranch==0 && aRecIndex[ijet]!=-1)  printf("gen branch jets: %2d matched with rec jet: %2d\n", ijet, aRecIndex[ijet]);

      // Fill the lists of matched generated and reconstructed jets
      if((iJetBranch==0)&&(aRecIndex[ijet]!=-1) ){
        listMatchedJets[iJetBranch]->Add( ((AliAODJet*) listJets[iJetBranch]->At(ijet)));
        listMatchedJets[1]->Add(((AliAODJet*) listJets[1]->At(aRecIndex[ijet])));
      }
   } // End loop on jets
  } // End loop on jet branches
 
 return listJets[1]->GetEntries();

}
//_______________________________________________________________________________
void AliAnalysisTaskJetFFMoments::AssociateGenRec(TList* tracksAODMCCharged, TList* tracksRec, TArrayI& indexAODTr, TArrayI& indexMCTr,  TArrayS& isRefGen, TH2D** fh2RecVsGen)
{
//   associate generated and reconstructed tracks, fill TArrays of list indices
//   if TH2D ** fh2RecVsGen not NULL, fill it (0:pt, 1:eta, 2:phi)

 Int_t nTracksRec  = tracksRec->GetSize();
 Int_t nTracksGen  = tracksAODMCCharged->GetSize();
 if(!fAOD) return; //added for kine
 TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
 if(!nTracksGen) return;
 if(!tca)return;

//  set size
 indexAODTr.Set(nTracksGen);
 indexMCTr.Set(nTracksRec);
 isRefGen.Set(nTracksGen);

 indexAODTr.Reset(-1);
 indexMCTr.Reset(-1);
  isRefGen.Reset(0);

//	loop over reconstructed tracks, get generated track
  for(Int_t iRec=0; iRec<nTracksRec; iRec++){

	AliAODTrack* rectrack = dynamic_cast<AliAODTrack*>(tracksRec->At(iRec));
	if(!rectrack)continue;
	Int_t label = TMath::Abs(rectrack->GetLabel());

	//	find MC track in our list
	AliAODMCParticle* gentrack = dynamic_cast<AliAODMCParticle*> (tca->At(label));

	if(rectrack->GetLabel() < 0 && fDebug > 5) {
	  printf("label < 0: tracks %d with charge %d asso. to MC particle (pdg = %d, rec_asso_id/genid = %d/%d) with charge %d\n", iRec, rectrack->Charge(), gentrack->GetPdgCode(), rectrack->GetLabel(), gentrack->GetLabel(),  gentrack->Charge());
	}

	Int_t listIndex = -1;
	if(gentrack) listIndex = tracksAODMCCharged->IndexOf(gentrack);
	if(listIndex>=0){
	  indexAODTr[listIndex] = iRec;
	  indexMCTr[iRec] = listIndex;
          fh2TrackResPt->Fill(gentrack->Pt(),100*(TMath::Abs(gentrack->Pt()-rectrack->Pt())/gentrack->Pt()));
          fh1TrackResPtInv->Fill(100. * (TMath::Abs(1/rectrack->Pt()) - 1./gentrack->Pt()) * gentrack->Pt());
	}
  } // End loop on rec tracks

  // define reference sample of primaries/secondaries (for reconstruction efficiency / contamination)

  for(Int_t iGen=0; iGen<nTracksGen; iGen++){

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (tracksAODMCCharged->At(iGen));
    if(!gentrack)continue;
    Int_t pdg = gentrack->GetPdgCode();

    // 211 - pi, 2212 - proton, 321 - Kaon, 11 - electron, 13 - muon
    if(TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 ||
       TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13){

      isRefGen[iGen] = kTRUE;

      Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track

      if(iRec>=0){
	AliAODTrack* vt = dynamic_cast<AliAODTrack*>(tracksRec->At(iRec));
	if(vt && fh2RecVsGen){
	  fh2RecVsGen[0]->Fill(gentrack->Pt(), vt->Pt());
	  fh2RecVsGen[1]->Fill(gentrack->Eta(), vt->Eta());
	  fh2RecVsGen[2]->Fill(gentrack->Phi(), vt->Phi());
	}
      }
    } // End if particle type
  } // End loop on gen particles
}

//_______________________________________________________________________________
void AliAnalysisTaskJetFFMoments::CalcFFAndFFMFromTracksInJet(AliAODJet* jet,  TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, TArrayI indexAODTr, const TArrayS& isRefGen, TList* jetTrackListTR, Bool_t scaleStrangeness,  TH2D** histo, TProfile2D * p2FFM)
{
  Double_t sumPt = 0.;
  Float_t alpha = TMath::Pi()/2.;
  TList* outputlist = 0x0;
  if(fFFBckgMode){
    outputlist = new TList();
  }

  Bool_t SelectTrackAsso;
  if(fTracksInJetMethod==0 || fTracksInJetMethod==1 || fTracksInJetMethod==2) SelectTrackAsso = kFALSE;
  else SelectTrackAsso = kTRUE;
  if(fTracksInJetMethod==0 || fTracksInJetMethod==1 || fTracksInJetMethod==2){
    TList* trackInJetAsso = new TList();
    trackInJetAsso->SetName(Form("%sAsso",jetTrackList->GetName()));
    TArrayD* sWeight = new TArrayD(jetTrackList->GetSize());
    if (tracksGen!=0) {
    if(!jetTrackListTR) {
      FillTrackAssoList(jet,jetTrackList,tracksGen, tracksRec, indexAODTr, isRefGen,0,scaleStrangeness, sWeight,trackInJetAsso);
    }
    else {
      FillTrackAssoList(jet,jetTrackList,tracksGen, tracksRec, indexAODTr, isRefGen,jetTrackListTR,scaleStrangeness, sWeight,trackInJetAsso);
    }

    if(fFFBckgMode)     GetTracksTiltedwrpJetAxis(alpha,trackInJetAsso,outputlist,jet,fRparam,sumPt);
    if(fFFBckgMode==2)  GetTracksTiltedwrpJetAxis(-alpha,trackInJetAsso,outputlist,jet,fRparam,sumPt);

    if(!fFFBckgMode) ComputeFFM(SelectTrackAsso,trackInJetAsso,indexAODTr,ComputeFF(SelectTrackAsso,jet,indexAODTr,trackInJetAsso,sWeight,histo),sWeight,p2FFM);
   else ComputeFFM(SelectTrackAsso,outputlist,indexAODTr,ComputeFF(SelectTrackAsso,jet,indexAODTr,outputlist,sWeight,histo),sWeight,p2FFM);
   }
   else {
    SelectTrackAsso = kFALSE;
    sWeight->Reset(1.);
    if(fFFBckgMode)     GetTracksTiltedwrpJetAxis(alpha,jetTrackList,outputlist,jet,fRparam,sumPt);
    if(fFFBckgMode==2)  GetTracksTiltedwrpJetAxis(-alpha,jetTrackList,outputlist,jet,fRparam,sumPt);

   if(!fFFBckgMode)  ComputeFFM(SelectTrackAsso,jetTrackList,indexAODTr,ComputeFF(SelectTrackAsso,jet,indexAODTr,jetTrackList,sWeight,histo),sWeight,p2FFM);
   else  ComputeFFM(SelectTrackAsso,outputlist,indexAODTr,ComputeFF(SelectTrackAsso,jet,indexAODTr,outputlist,sWeight,histo),sWeight,p2FFM);
   }
    delete trackInJetAsso;
    delete sWeight;
  }
  else if(fTracksInJetMethod==3){

    TArrayD* sWeight = new TArrayD(jetTrackList->GetSize());
    sWeight->Reset(1.);

    if(fFFBckgMode)     GetTracksTiltedwrpJetAxis(alpha,jetTrackList,outputlist,jet,fRparam,sumPt);
    if(fFFBckgMode==2)  GetTracksTiltedwrpJetAxis(-alpha,jetTrackList,outputlist,jet,fRparam,sumPt);

    //Mengliang's method
    if(fkUseTrackPtSumAsJetPt) {
     if(!fFFBckgMode) ComputeFFM(SelectTrackAsso, jetTrackList, indexAODTr, ComputeFF(SelectTrackAsso, NULL, indexAODTr,  jetTrackList, sWeight, histo), sWeight, p2FFM);
     else ComputeFFM(SelectTrackAsso, outputlist, indexAODTr, ComputeFF(SelectTrackAsso, NULL, indexAODTr, outputlist, sWeight, histo), sWeight, p2FFM);
   } else {
     if(!fFFBckgMode) ComputeFFM(SelectTrackAsso, jetTrackList, indexAODTr, ComputeFF(SelectTrackAsso, jet, indexAODTr, jetTrackList, sWeight, histo), sWeight, p2FFM);
     else ComputeFFM(SelectTrackAsso, outputlist, indexAODTr, ComputeFF(SelectTrackAsso, jet, indexAODTr, outputlist, sWeight, histo), sWeight, p2FFM);
  }

  }
  else printf("Unkown track in jet method %d",fTracksInJetMethod);

  if(fFFBckgMode) delete outputlist;
}

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFFMoments::ComputeFFM(Bool_t SelectTrackAsso, TList* tracks,  TArrayI indexTr, Float_t JetPt, TArrayD* sWeight, TProfile2D * p2FFM){

    const Int_t nFFBins = fFFMMomMax;
    Int_t nTracks =  0.;
    if(SelectTrackAsso) nTracks = indexTr.GetSize();
    else nTracks = tracks->GetSize();

    if(fkUseTrackPtSumAsJetPt) {
        Float_t SumPt = 0;
        //  loop over reconstructed tracks, get SumPt
        for(Int_t i=0; i<nTracks; i++) {
          if(SelectTrackAsso&&!fFFBckgMode){
          if(indexTr.At(i) == -1) continue;
        }
        AliVParticle* track = dynamic_cast<AliVParticle*>(tracks->At(i));
        if(!track) continue;
        SumPt += track->Pt();
       }
       if(SumPt <= 0) { printf("strange Pt_Sum <= 0."); /*return;*/ }
       else JetPt = SumPt;
    }

    for(Int_t iN = 0; iN< nFFBins; iN++) {
      Double_t N = iN * (fFFMNMax - fFFMNMin)/(fFFMMomMax-1) + fFFMNMin;
      Double_t A = TMath::Power( 0.5*(N + 1), fFFMScalePower);
      Double_t fFFMValue = 0.;
      Double_t weight = 0;
      for(Int_t i=0; i<nTracks; i++){
	if(SelectTrackAsso&&!fFFBckgMode) if(indexTr.At(i) == -1) continue;
	AliVParticle* track = 0x0;
        track = dynamic_cast<AliVParticle*>(tracks->At(i));
	if(!track) continue;
	Double_t z = track->Pt()/JetPt;
	fFFMValue += TMath::Power(z, N);
        weight += sWeight->At(i); // Sum sWeight
      }
      if(nTracks!=0) {weight /= nTracks;} else {weight = 1; }
      p2FFM->Fill(N, JetPt, A*fFFMValue,weight);
    }

}

// ________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFFMoments::GetJetTracksPointing(TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius, Double_t& sumPt, Double_t minPtL, Double_t maxPt, Bool_t& isBadPt)
{
  // fill list of tracks in cone around jet axis
  sumPt = 0;
  Bool_t isBadMaxPt = kFALSE;
  Bool_t isBadMinPt = kTRUE;

  Double_t jetMom[3];
  jet->PxPyPz(jetMom);
  TVector3 jet3mom(jetMom);

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
    if(!track)continue;
    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);

    Double_t dR = jet3mom.DeltaR(track3mom);

    if(dR<radius){

      outputlist->Add(track);
      sumPt += track->Pt();

      if(maxPt>0  && track->Pt()>maxPt)  isBadMaxPt = kTRUE;
      if(minPtL>0 && track->Pt()>minPtL) isBadMinPt = kFALSE;
    }
  }

  isBadPt = kFALSE;
  if(minPtL>0 && isBadMinPt) isBadPt = kTRUE;
  if(maxPt>0  && isBadMaxPt) isBadPt = kTRUE;

  outputlist->Sort();
}

// _________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFFMoments::GetJetTracksTrackrefs(TList* list, const AliAODJet* jet, Double_t minPtL, Double_t maxPt, Bool_t& isBadPt)
{
  // list of jet tracks from trackrefs
  Int_t nTracks = jet->GetRefTracks()->GetEntriesFast();

  Bool_t isBadMaxPt = kFALSE;
  Bool_t isBadMinPt = kTRUE;

  for(Int_t itrack=0; itrack<nTracks; itrack++) {

    AliVParticle* track = dynamic_cast<AliVParticle*>(jet->GetRefTracks()->At(itrack));
    if(!track){
      AliError("expected ref track not found ");
      continue;
    }

    if(track->Pt()  < fTrackPtMin) continue; // track refs may contain low pt cut (bug in AliFastJetInput)
    if(maxPt>0 && track->Pt()>maxPt)   isBadMaxPt = kTRUE;
    if(minPtL>0 && track->Pt()>minPtL) isBadMinPt = kFALSE;

    list->Add(track);
  }

  isBadPt = kFALSE;
  if(minPtL>0 && isBadMinPt) isBadPt = kTRUE;
  if(maxPt>0 && isBadMaxPt)  isBadPt = kTRUE;

  list->Sort();

}

// _______________________________________________________________________________
void AliAnalysisTaskJetFFMoments::CreateHistos() 
{
  if (fDebug != 0) printf("AliAnalysisTaskJetFFMoments::CreateHistos() start!\n");

  fPrParameters = new TProfile("h1Parameters","Parameters for analysis", 31, 0.5, 31.5);
  Int_t axis = 1;
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"IsPbPb");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FFMScalePower");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"Kine_Input");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"Reco_Input");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"DoJetReco");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"NUsed(di)Jets");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"DoMatching");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"MatchMaxDist");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"JetMatchingFractionMin");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FFBckgMode");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FFMBckgType");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FFMBckgPar1");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FFMBckgPar2");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FFMBckgMu/10");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"JetAlgorithm");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"JetR");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"JetPtMin");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"|JetEtaMin|");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"JetEtaMax");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FilterMask/100");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"TrackPtMin");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"|TrackEtaMin|");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"TrackEtaMax");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"TrackExtraCut");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"UseTrackPtSumAsJetPt");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"TracksInJetMethod");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"GenJetType");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"EffJetType");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FastTrackEfficiencySyst/100");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FastTrackResolution");
  fPrParameters->GetXaxis()->SetBinLabel(axis++,"FastTrackResolutionSyst/100");

  fh1Xsec = new TProfile("hpXsec","Pythia X section",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");

  fh1Trials = new TH1F("h1Trials","Pythia number of trials",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");

  fh1AvgTrials = new TH1F("h1AvgTrials","Average number of trials",1,0,1);
  fh1AvgTrials->GetXaxis()->SetBinLabel(1,"#sum{avg ntrials}");

  fh1CentralitySelect = new TH1F("h1CentralitySelect", "Centrality with physics selection bypassed; cent (%)", 102, -1.5, 100.5);
  fh1CentralityPhySel = new TH1F("h1CentralityPhySel", "Centrality after physics selection; cent (%)", 102, -1.5, 100.5);

  fh1vZSelect = new TH1F("h1vZSelect", "vZ with physics selection bypassed; zvtx", 100, -25, 25);
  fh1vZPhySel = new TH1F("h1vZPhySel", "vZ after physics selection; zvtx", 100, -25, 25);

 
  TString fAnaObjectType = "";
  TString sHisName;
  TString sAxisName[10];
  sAxisName[0]="pt";	sAxisName[1]="eta";	sAxisName[2]="phi";	sAxisName[3]="area";	sAxisName[4]="nCon";
  sAxisName[5]="z";	sAxisName[6]="xi";	sAxisName[7]="jT";	sAxisName[8]="Dtheta";	sAxisName[9]="non";

  UInt_t entries = 0; // bit coded, see GetDimParams() below
  UInt_t res = 0;  // detail high resolution
  
  // 0 - 9
  //  0 |  1  |  2  |  3   |  4  |  5 |   6   |      7      |    8    |    9    |
  // vz | ntr | ep  | epb  |  z  | xi | lnjT  |  DeltaTheta | FFM_gen | FFM_rec |
  // 10 - 24
  //  0 |  1  |  2  |    3     |      4        |||      n
  // pt | eta | phi |   area   | nconstituents |||     gen  2x5+n
  // pt | eta | phi |   area   | nconstituents |||     rec  3x5+n
  // pt | eta | phi |          |    N(order)   |||    track 4x5+n
  // 25 - 32
  //	 25	|	26	|
  //  fraction | RDeltaPt |  R for Relative
  //
  sHisName = TString("h1GenTracks");
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
    fh1GenTracks[iAxis] = CreateTH1D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, fkHighResolution);
  }
  sHisName = TString("h1RecTracks");
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
    fh1RecTracks[iAxis] = CreateTH1D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, fkHighResolution);
  }
  sHisName = TString("h2AssociatedTracks");
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
	fh2PtRecVsGenPrim[iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, 20+iAxis, fkHighResolution);
	fh2PtRecVsGenPrim[iAxis]->GetXaxis()->SetTitle(Form("%s_{gen}", fh2PtRecVsGenPrim[iAxis]->GetXaxis()->GetTitle()));
	fh2PtRecVsGenPrim[iAxis]->GetYaxis()->SetTitle(Form("%s_{rec}", fh2PtRecVsGenPrim[iAxis]->GetYaxis()->GetTitle()));
 }
    sHisName = TString("h2AssociatedSecTracks");
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
	fh2PtRecVsGenSec[iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, 20+iAxis, fkHighResolution);
	fh2PtRecVsGenSec[iAxis]->GetXaxis()->SetTitle(Form("%s_{gen}", fh2PtRecVsGenSec[iAxis]->GetXaxis()->GetTitle()));
	fh2PtRecVsGenSec[iAxis]->GetYaxis()->SetTitle(Form("%s_{rec}", fh2PtRecVsGenSec[iAxis]->GetYaxis()->GetTitle()));
  }

   sHisName = TString("h1TrackEffPtGen");
   fh1TrackEffPtGen =  CreateTH1D(sHisName.Data(), 20, fkHighResolution);
   sHisName = TString("h2TrackEffEtaPhiGen");
   fh2TrackEffEtaPhiGen =  CreateTH2D(sHisName.Data(), 21, 22, fkHighResolution);
   sHisName = TString("h1TrackEffPtRec");
   fh1TrackEffPtRec =  CreateTH1D(sHisName.Data(), 20, fkHighResolution);
   sHisName = TString("h2TrackEffEtaPhiRec");
   fh2TrackEffEtaPhiRec =  CreateTH2D(sHisName.Data(), 21, 22, fkHighResolution);

   sHisName = TString("h2TrackResPt");
   fh2TrackResPt =  CreateTH2D(sHisName.Data(), 20, 27, fkHighResolution);
   sHisName = TString("h2TrackResPtInv");
   fh1TrackResPtInv =  CreateTH1D(sHisName.Data(), 27, fkHighResolution);

  sHisName = TString("h2PtHardCutVsPythiaJetPt");
  fh2PtHardVsPtCut[0] = CreateTH2D(sHisName.Data(), 20, 20, fkHighResolution);
  sHisName = TString("h2PtHardCutVsTrackPt");
  fh2PtHardVsPtCut[1] = CreateTH2D(sHisName.Data(), 20, 20, fkHighResolution);
  sHisName = TString("h2PtHardVsPythiaJetPt");
  fh2PtHardVsPt[0] = CreateTH2D(sHisName.Data(), 20, 20, fkHighResolution);
  sHisName = TString("h2PtHardVsTrackPt");
  fh2PtHardVsPt[1] = CreateTH2D(sHisName.Data(), 20, 20, fkHighResolution);

  for(int iJetBranch=0; iJetBranch < fgkFFMNJetBranches; iJetBranch++){
    if( iJetBranch == 0) fAnaObjectType = "Gen";
    else { fAnaObjectType = "Rec"; }

    sHisName = Form("h1MisMatched_%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
    for( Int_t iAxis = 0; iAxis < 5; iAxis++) fh1JetPr_Mismatched[iJetBranch][iAxis] = CreateTH1D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), iJetBranch?15+iAxis:10+iAxis, fkHighResolution);
	fh2MismatchedJetsAreaVSPt[iJetBranch] = CreateTH2D(Form("h2MisMatched_%s_%s_Area", fAnaJetType.Data(), fAnaObjectType.Data()), 13+iJetBranch*5, 10+iJetBranch*5, fkHighResolution);

	sHisName = Form("h2TracksIn%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
	for( Int_t iAxis = 0; iAxis < 4; iAxis++ ) fh2TracksInJets[iJetBranch][iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis+5].Data()), 4+iAxis, iJetBranch?15:10, fkHighResolution);
	for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh2TracksInJets[iJetBranch][iAxis+4] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, iJetBranch?15:10, fkHighResolution);

      if( iJetBranch == 1) {
	sHisName = Form("h2AssociatedTracksIn%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
	for( Int_t iAxis = 0; iAxis < 4; iAxis++ ) fh2AssociatedTracksInJets[iJetBranch][iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis+5].Data()), 4+iAxis, iJetBranch?15:10, fkHighResolution);
	for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh2AssociatedTracksInJets[iJetBranch][iAxis+4] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, iJetBranch?15:10, fkHighResolution);

         sHisName = Form("h2AssociatedTracksSecNSIn%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
	for( Int_t iAxis = 0; iAxis < 4; iAxis++ ) fh2AssociatedTracksInJetsSecNS[iJetBranch][iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis+5].Data()), 4+iAxis, iJetBranch?15:10, fkHighResolution);
	for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh2AssociatedTracksInJetsSecNS[iJetBranch][iAxis+4] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, iJetBranch?15:10, fkHighResolution);
        sHisName = Form("h2AssociatedTracksSecSIn%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());

	for( Int_t iAxis = 0; iAxis < 4; iAxis++ ) fh2AssociatedTracksInJetsSecS[iJetBranch][iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis+5].Data()), 4+iAxis, iJetBranch?15:10, fkHighResolution);
	for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh2AssociatedTracksInJetsSecS[iJetBranch][iAxis+4] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, iJetBranch?15:10, fkHighResolution);
        sHisName = Form("h2AssociatedTracksSecSscIn%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
	for( Int_t iAxis = 0; iAxis < 4; iAxis++ ) fh2AssociatedTracksInJetsSecSsc[iJetBranch][iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis+5].Data()), 4+iAxis, iJetBranch?15:10, fkHighResolution);
	for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh2AssociatedTracksInJetsSecSsc[iJetBranch][iAxis+4] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, iJetBranch?15:10, fkHighResolution);

      sHisName = Form("h2RecJetVsGenJet_%s", fAnaJetType.Data());
      for( Int_t iAxis = 0; iAxis < 5; iAxis++) fh2MatchedJets[iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 10+iAxis, 15+iAxis, fkHighResolution);

      if( fHistosLevel >= 9 ) {
	entries   = 1<<8 | 1<<9 | 1<<10 | 1<<15 | 1<<25;
	if(fkHighResolution) res = 1<<10 | 1<<15 | 1<<25;

        const Int_t nFFBins = fFFMMomMax;
        fhnJetMomN_Raw = new THnSparse*[nFFBins];
        fhnJetMomN_Sub = new THnSparse*[nFFBins];
        fhnJetMomN_Imp = new THnSparse*[nFFBins];
        for( Int_t iN = 0; iN < fFFMMomMax; iN++) {
          Double_t fN = iN * (fFFMNMax - fFFMNMin)/(fFFMMomMax-1) + fFFMNMin;
	  fhnJetMomN_Raw[iN] =  CreateTHnSparseF(Form("hnFFM_Raw_M_{%.2f}_%s", fN, fAnaJetType.Data()), entries, res);
	  fhnJetMomN_Sub[iN] =  CreateTHnSparseF(Form("hnFFM_Sub_M_{%.2f}_%s", fN, fAnaJetType.Data()), entries, res);
	  fhnJetMomN_Imp[iN] =  CreateTHnSparseF(Form("hnFFM_Imp_M_{%.2f}_%s", fN, fAnaJetType.Data()), entries, res);
        } // End loop on iN
      } // End histo level 
    } // End iJetBranch==1

    if( iJetBranch == 0 ) {
      entries = 1<<8 | 1<<10 | 1<< 24;
      if(fkHighResolution) res = 1<<10;
    } else if ( iJetBranch == 1 ) {
      entries = 1<<9 | 1<<15 | 1<< 24;
      if(fkHighResolution) res = 1<<15;
    }
	TString FFM_suff = Form("%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
	fhnJetFFM_Raw[iJetBranch] = CreateTHnSparseF(Form("hnFFM_Raw_%s", FFM_suff.Data()), entries, res);
	fhnJetFFM_Sub[iJetBranch] = CreateTHnSparseF(Form("hnFFM_Sub_%s", FFM_suff.Data()), entries, res);
	fhnJetFFM_Imp[iJetBranch] = CreateTHnSparseF(Form("hnFFM_Imp_%s", FFM_suff.Data()), entries, res);
        fp2TracksInJetFFM[iJetBranch] = CreateTProfile2D(Form("p2FFM_Tracks_%s",FFM_suff.Data()), 24, 10+5*iJetBranch, 8, res);
	if( iJetBranch == 1) {
        fp2AssociatedTracksJetFFM[iJetBranch] = CreateTProfile2D(Form("p2FFM_AssociatedTracks_%s",FFM_suff.Data()), 24, 10+5*iJetBranch, 8, res);
        fp2AssociatedTracksJetFFMSecNS[iJetBranch] = CreateTProfile2D(Form("p2FFM_AssociatedTracksSecNS_%s",FFM_suff.Data()), 24, 10+5*iJetBranch, 8, res);
	fp2AssociatedTracksJetFFMSecS[iJetBranch] = CreateTProfile2D(Form("p2FFM_AssociatedTracksSecS_%s",FFM_suff.Data()), 24, 10+5*iJetBranch, 8, res);
	fp2AssociatedTracksJetFFMSecSsc[iJetBranch] = CreateTProfile2D(Form("p2FFM_AssociatedTracksSecSsc_%s",FFM_suff.Data()), 24, 10+5*iJetBranch, 8, res);
        }
	fp2JetFFM_Raw[iJetBranch] = CreateTProfile2D(Form("p2FFM_Raw_%s",FFM_suff.Data()), 24, 10+5*iJetBranch, 8, res);
	fp2JetFFM_Sub[iJetBranch] = CreateTProfile2D(Form("p2FFM_Sub_%s",FFM_suff.Data()), 24, 10+5*iJetBranch, 8, res);
	fp2JetFFM_Imp[iJetBranch] = CreateTProfile2D(Form("p2FFM_Imp_%s",FFM_suff.Data()), 24, 10+5*iJetBranch, 8, res);

    if( fAnaJetType.Contains("DIJET")) {
      fh1Asy_DiJets[iJetBranch]  = new TH1F(Form("h1Asy_DiJets_%s", fAnaObjectType.Data()), Form("Asymmetry DiJets_%s", fAnaObjectType.Data()), 20, 0.0 , 1.0);
      fh1Asy_DiJets[iJetBranch]->GetXaxis()->SetTitle("Asymmetry = #frac{p_{T}^1-p_{T}^2}{p_{T}^1+p_{T}}^2");
    }
	fh2MatchedJetsRDPtVSPt[iJetBranch] = CreateTH2D(Form("h2Matched_%s_%s_RelativeDeltaPt", fAnaJetType.Data(), fAnaObjectType.Data()), 26, 10+iJetBranch*5, fkHighResolution);
	fh2MatchedJetsRDPtVSPt[iJetBranch]->SetXTitle(Form("%s/%s", fh2MatchedJetsRDPtVSPt[iJetBranch]->GetXaxis()->GetTitle(), fh2MatchedJetsRDPtVSPt[iJetBranch]->GetYaxis()->GetTitle()));

        fh2MatchedJetsUE[iJetBranch] = CreateTH2D(Form("h2MatchedUE_%s_%s_pt", "LEADING", fAnaObjectType.Data()), 10+iJetBranch*5, 28 ,0);

	fh2MatchedJetsAreaVSPt[iJetBranch] = CreateTH2D(Form("h2Matched_%s_%s_Area", fAnaJetType.Data(), fAnaObjectType.Data()), 13+iJetBranch*5, 10+iJetBranch*5, fkHighResolution);

    fh1Njets[iJetBranch] = new TH1F(Form("h1Njets_%s", fAnaObjectType.Data()), "Number of jets after cut and before matching", 9, 0.5 , 9.5);
    fh1Njets[iJetBranch]->GetXaxis()->SetTitle("N_{jet}");

  } // End loop on iJetBranch

  //events
  if( fkIsPbPb ) {
    fHistList->Add(fh1CentralitySelect);
    fHistList->Add(fh1CentralityPhySel);
  }
  fHistList->Add(fh1vZSelect);
  fHistList->Add(fh1vZPhySel);
  //tracks & jets
  fHistList->Add(fPrParameters);
  fHistList->Add(fh1Xsec);
  fHistList->Add(fh1Trials);
  fHistList->Add(fh1AvgTrials);
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
    fHistList->Add(fh1GenTracks[iAxis]);
  }
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
    fHistList->Add(fh1RecTracks[iAxis]);
  }

 if(!((!fJetBranch[0].Contains("MC") && !fkDoJetReco) || (fkDoJetReco && !(fTrackType[0] == kTrackAODMCCharged || fTrackType[0] ==  kTrackAODMCChargedAcceptance  || fTrackType[0] == kTrackAODMCextra || fTrackType[0] == kTrackAODMCextraonly))) ||  fTrackType[1] < kTrackAOD ){
	for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
	  fHistList->Add(fh2PtRecVsGenPrim[iAxis]);
	  fHistList->Add(fh2PtRecVsGenSec[iAxis]);
	}

    fHistList->Add(fh1TrackEffPtGen);
    fHistList->Add(fh1TrackEffPtRec);
    fHistList->Add(fh2TrackEffEtaPhiGen);
    fHistList->Add(fh2TrackEffEtaPhiRec);
    fHistList->Add(fh2TrackResPt);
    fHistList->Add(fh1TrackResPtInv);
    fHistList->Add(fh2PtHardVsPt[0]);
    fHistList->Add(fh2PtHardVsPt[1]);
    fHistList->Add(fh2PtHardVsPtCut[0]);
    fHistList->Add(fh2PtHardVsPtCut[1]);

    if (fEffi && fTrackType[1] < kTrackAOD) fHistList->Add(fEffi);
  }

  for( Int_t iAxis = 0; iAxis < 5; iAxis++ ) fHistListJets[1]->Add(fh2MatchedJets[iAxis]);
  if( fHistosLevel >= 9 ) {
    for( Int_t iN = 0; iN < fFFMMomMax; iN++) {
      fHistListJets[1]->Add(fhnJetMomN_Raw[iN]);
      fHistListJets[1]->Add(fhnJetMomN_Sub[iN]);
      fHistListJets[1]->Add(fhnJetMomN_Imp[iN]);
    }
  }

  for(int iJetBranch=0; iJetBranch < fgkFFMNJetBranches; iJetBranch++){
    fHistList->Add(fHistListJets[iJetBranch]);
    fHistListJets[iJetBranch]->Add(fh1Njets[iJetBranch]);
    for( Int_t iAxis = 0; iAxis < 5; iAxis++ ) fHistListJets[iJetBranch]->Add(fh1JetPr_Mismatched[iJetBranch][iAxis]);
	fHistListJets[iJetBranch]->Add(fh2MismatchedJetsAreaVSPt[iJetBranch]);
    if( fAnaJetType.Contains("DIJET")) fHistListJets[iJetBranch]->Add(fh1Asy_DiJets[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fh2MatchedJetsRDPtVSPt[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fh2MatchedJetsAreaVSPt[iJetBranch]);
        fHistListJets[iJetBranch]->Add(fh2MatchedJetsUE[iJetBranch]);
	for( Int_t iAxis = 0; iAxis < 7; iAxis++) fHistListJets[iJetBranch]->Add(fh2TracksInJets[iJetBranch][iAxis]);
	if((fJetBranch[0].Contains("MC") || fkDoJetReco) && ( !fkDoJetReco ||(fTrackType[0] == kTrackAODMCCharged || fTrackType[0] ==  kTrackAODMCChargedAcceptance  || fTrackType[0] == kTrackAODMCextra || fTrackType[0] == kTrackAODMCextraonly))){
	  for( Int_t iAxis = 0; iAxis < 7; iAxis++) fHistListJets[iJetBranch]->Add(fh2AssociatedTracksInJets[iJetBranch][iAxis]);
	  for( Int_t iAxis = 0; iAxis < 7; iAxis++) fHistListJets[iJetBranch]->Add(fh2AssociatedTracksInJetsSecNS[iJetBranch][iAxis]);
	  for( Int_t iAxis = 0; iAxis < 7; iAxis++) fHistListJets[iJetBranch]->Add(fh2AssociatedTracksInJetsSecS[iJetBranch][iAxis]);
	  for( Int_t iAxis = 0; iAxis < 7; iAxis++) fHistListJets[iJetBranch]->Add(fh2AssociatedTracksInJetsSecSsc[iJetBranch][iAxis]);
	}
    fHistListJets[iJetBranch]->Add(fhnJetFFM_Raw[iJetBranch]);
    fHistListJets[iJetBranch]->Add(fhnJetFFM_Sub[iJetBranch]);
    fHistListJets[iJetBranch]->Add(fhnJetFFM_Imp[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fp2TracksInJetFFM[iJetBranch]);
        fHistListJets[iJetBranch]->Add(fp2AssociatedTracksJetFFM[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fp2AssociatedTracksJetFFMSecNS[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fp2AssociatedTracksJetFFMSecS[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fp2AssociatedTracksJetFFMSecSsc[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fp2JetFFM_Raw[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fp2JetFFM_Sub[iJetBranch]);
	fHistListJets[iJetBranch]->Add(fp2JetFFM_Imp[iJetBranch]);
  }

  if (fDebug != 0) printf("AliAnalysisTaskJetFFMoments::CreateHistos() end!\n");
}

//----------------------------------------------------------------------
fastjet::PseudoJet  AliAnalysisTaskJetFFMoments::join_with_area(const vector <fastjet::PseudoJet> & pieces,
								const fastjet::PseudoJet & area_4vector,
								const double area,
								const double area_error,
								const bool   is_pure_ghost)
{
  /// mimics the join function, supplementing it with user-specified
  /// areas
  // compute the total momentum
  //--------------------------------------------------
  fastjet::PseudoJet result;  // automatically initialised to 0
  for (unsigned int i=0; i<pieces.size(); i++)
    result += pieces[i];

  // attach a CompositeJetStructure to the result
  //--------------------------------------------------
  AliCompositeJetStructureUserArea *cj_struct =  new AliCompositeJetStructureUserArea(pieces, area_4vector, area, area_error, is_pure_ghost);

  result.set_structure_shared_ptr(fastjet::SharedPtr<fastjet::PseudoJetStructureBase>(cj_struct));

  return result;
}

//_____________________________________________________________
void AliAnalysisTaskJetFFMoments::PseudoJetsToAODJets( const vector<fastjet::PseudoJet> & fPseudojets, fastjet::ClusterSequenceArea & CluSeq, TList & tracks, TList *jets)
{
  //
  // Fill a list of AliAODJets from a table of pseudojets
  //

  Int_t nj = 0;
  for (size_t j = 0; j < fPseudojets.size(); j++) {
    Double_t area      = CluSeq.area(fPseudojets[j]);
    Double_t areaError = CluSeq.area_error(fPseudojets[j]);
    vector<fastjet::PseudoJet> constituents = CluSeq.constituents(fPseudojets[j]);
    int nCon= constituents.size();
    TArrayI ind(nCon);
    AliAODJet aodjet(fPseudojets[j].px(), fPseudojets[j].py(), fPseudojets[j].pz(), fPseudojets[j].E());
    aodjet.SetEffArea(area,areaError);
    fastjet::PseudoJet vecarea=CluSeq.area_4vector(fPseudojets[j]);
    TLorentzVector vecareab;
    vecareab.SetPxPyPzE(vecarea.px(),vecarea.py(),vecarea.pz(),vecarea.e());
    aodjet.SetVectorAreaCharged(&vecareab);
    for (int i=0; i < nCon; i++) {
      fastjet::PseudoJet mPart=constituents[i];
      ind[i]=mPart.user_index();
  
      Int_t ntracks = tracks.GetEntries();
      for(Int_t itrack=0; itrack<ntracks; itrack++) {
        if(itrack==ind[i]) {
          TObject *track = tracks.At(itrack);
	  aodjet.AddTrack(track);
        }
      }
    } // End loop on Constituents
	if ( SelectJet(&(aodjet), &tracks) ) {
      jets->Add(new AliAODJet(aodjet));
      if (fTCAJetsOut) { new ((*fTCAJetsOut)[nj++]) AliAODJet(aodjet);}
      if(fDebug > 2) {aodjet.Print("");}
    }
  }
}

//_____________________________________________________________
int AliAnalysisTaskJetFFMoments::AliAODJetToPseudoJet(AliAODJet* jet, fastjet::PseudoJet & fCurrentPseudojet)
{ 
  //
  // Fill a pseudo jet from an AliAODJet
  //

  if(jet == 0) { if( fDebug ) printf("In AliAnalysisTaskJetFFMoments::AliAODJetToPseudoJet(): Could not get the AOD jet!\n" ); return -1;}

  // Fill charged tracks as constituents
  Int_t nConstituents =  jet->GetRefTracks()->GetEntries() ;
  if( nConstituents == 0) { printf("In AliAnalysisTaskJetFFMoments::AliAODJetToPseudoJet(): no constituents in jets! The pseudojet can not be filled!!\n "); return -2;}
  vector<fastjet::PseudoJet> inputParticles;
  for(Int_t i = 0; i < nConstituents; i++)  { // loop for all input tracks
    AliAODTrack* track = (AliAODTrack *) (jet->GetTrack(i));
    fastjet::PseudoJet inputPart(track->Px(),track->Py(),track->Pz(),track->P());  // create PseudoJet object
    inputPart.set_user_index(i);      //label the particle into Fastjet algortihm
    inputParticles.push_back(inputPart); 
  } // End loop on Trk

  // Get the area and area 4vector
  Float_t area = jet->EffectiveAreaCharged();
  Float_t areaerror = jet->ErrorEffectiveAreaCharged();
  TLorentzVector* areav = jet->VectorAreaCharged();
  if(!areav) { AliError("The area 4-vector is not stored for the current jet! - Please add it..."); return -3;}
  fastjet::PseudoJet area_4vector(areav->Px(), areav->Py(), areav->Pz(), areav->E());

  // Fill the pseudo jet
  fCurrentPseudojet = join_with_area(inputParticles, area_4vector, area, areaerror);
  
  return 0;
}

//_____________________________________________________________
int AliAnalysisTaskJetFFMoments::AliAODJetToPseudoJet(AliAODJet* jet, TList* list,  fastjet::PseudoJet & fCurrentPseudojet)
{
  //
  // Fill a pseudo jet from an AliAODJet
  //

  if(jet == 0) { if( fDebug ) printf("In AliAnalysisTaskJetFFMoments::AliAODJetToPseudoJet(): Could not get the AOD jet!\n" ); return -1;}

  // Fill charged tracks as constituents
  Int_t nConstituents =  list->GetEntries() ;
  if( nConstituents == 0) { printf("In AliAnalysisTaskJetFFMoments::AliAODJetToPseudoJet(): no constituents in jets! The pseudojet can not be filled!!\n "); return -2;}
  vector<fastjet::PseudoJet> inputParticles;
  for(Int_t i = 0; i < nConstituents; i++)  { // loop for all input tracks
    AliAODTrack* track = (AliAODTrack *) (list->At(i));
    fastjet::PseudoJet inputPart(track->Px(),track->Py(),track->Pz(),track->P()); //E? // create PseudoJet object
    inputPart.set_user_index(i);      //label the particle into Fastjet algortihm
    inputParticles.push_back(inputPart);
  } // End loop on Trk

  // Get the area and area 4vector
  Float_t area = jet->EffectiveAreaCharged();
  Float_t areaerror = jet->ErrorEffectiveAreaCharged();
  TLorentzVector* areav = jet->VectorAreaCharged();
  if(!areav) { AliError("The area 4-vector is not stored for the current jet! - Please add it..."); return -3;}
  fastjet::PseudoJet area_4vector(areav->Px(), areav->Py(), areav->Pz(), areav->E());

  // Fill the pseudo jet
  fCurrentPseudojet = join_with_area(inputParticles, area_4vector, area, areaerror);

  return 0;
}

//______________________________________
Bool_t AliAnalysisTaskJetFFMoments::ClassifyJetEvent(Int_t ijet, Int_t njets)
{

  Bool_t result = kTRUE;
  if (     fAnaJetType.Contains("LEADING") ) {
    if ( njets>=1 && ijet == 1 )  { result = kFALSE;}
  } else if( fAnaJetType.Contains("DIJET") ) {
    if ( njets==1 || (njets>=2 &&  ijet == 2) ) { result = kFALSE; }
  } else if( fAnaJetType.Contains("ALLJET")) {
    { result = kTRUE; }
  } else AliError("Wrong jet analysis type!!");
  return result;
}
//______________________________________
Bool_t AliAnalysisTaskJetFFMoments::SelectJet(AliAODJet * jet, TList * tracksAfterCut)
{

  if( !jet ) {
    if(fDebug > 9) printf("no jet exist!\n");
    return kFALSE;
  }
  if( jet->Eta() > fJetEtaMax || jet->Eta() < fJetEtaMin ) {
    if(fDebug > 9) printf("eta out of range: %f\n", jet->Eta());
    return kFALSE;
  }
  if( jet->Pt()  < fJetPtMin ) {
    if(fDebug > 9) printf("pT out of range: %f\n", jet->Pt());
    return kFALSE;
  }


  //tracks in jet should be in list of tracksAftercut
  TRefArray * tracks = jet->GetRefTracks();
  Int_t ntracks = tracks -> GetEntriesFast();
  if (ntracks < fJetMinnTracks) return kFALSE;
  for(Int_t i=0; i < ntracks; i++) { 
	AliVParticle* track = (AliVParticle*) (jet->GetTrack(i));
	Int_t track_inlist = tracksAfterCut->IndexOf(track);
	if( track_inlist < 0) return kFALSE;
        if( track->Pt() > fTrackPtMax) return kFALSE;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________________________
TH1D * AliAnalysisTaskJetFFMoments::CreateTH1D(const char* name, Int_t iXAxis, Bool_t res)
{

  TString hnTitle(name);
  Int_t nbins;
  Double_t xmin;
  Double_t xmax;

  const char* label;
  GetDimParams(iXAxis, res, label, nbins, xmin, xmax);
  hnTitle += Form(";%s;",label);

  TH1D* h1tmp = new TH1D(name, hnTitle.Data(), nbins, xmin, xmax);
  h1tmp->Sumw2();

  return h1tmp;

}
//_____________________________________________________________________________________________
TProfile2D * AliAnalysisTaskJetFFMoments::CreateTProfile2D(const char* name, Int_t iXAxis, Int_t iYAxis, Int_t iZAxis, Bool_t res)
{

  TString p2Title(name);
  Int_t nbins[3];
  Double_t xmin[3];
  Double_t xmax[3];

  const char* label;
  GetDimParams(iXAxis, res, label, nbins[0], xmin[0], xmax[0]);
  p2Title += Form(";%s",label);
  GetDimParams(iYAxis, res, label, nbins[1], xmin[1], xmax[1]);
  p2Title += Form(";%s",label);
  GetDimParams(iZAxis, res, label, nbins[2], xmin[2], xmax[2]);
  p2Title += Form(";%s",label);

//  p2Title += ";";

  TProfile2D* p2tmp =  new TProfile2D(name, p2Title.Data(), nbins[0], xmin[0], xmax[0], nbins[1], xmin[1], xmax[1]);
  p2tmp->Sumw2();

  return p2tmp;
}
//_____________________________________________________________________________________________
TH2D * AliAnalysisTaskJetFFMoments::CreateTH2D(const char* name, Int_t iXAxis, Int_t iYAxis, Bool_t res)
{

  TString hnTitle(name);
  Int_t nbins[2];
  Double_t xmin[2];
  Double_t xmax[2];

  const char* label;
  GetDimParams(iXAxis, res, label, nbins[0], xmin[0], xmax[0]);
  hnTitle += Form(";%s",label);
  GetDimParams(iYAxis, res, label, nbins[1], xmin[1], xmax[1]);
  hnTitle += Form(";%s",label);

  hnTitle += ";";

  TH2D* h2tmp =  new TH2D(name, hnTitle.Data(), nbins[0], xmin[0], xmax[0], nbins[1], xmin[1], xmax[1]);
  h2tmp->Sumw2();

  return h2tmp;
}
//_______________________________________________________________________________________________
THnSparse* AliAnalysisTaskJetFFMoments::CreateTHnSparseF(const char* name, UInt_t entries, UInt_t res)
{

  // generate new THnSparseF, axes are defined in GetDimParams()
  Int_t count = 0;
  UInt_t tmp = entries;
  while(tmp!=0){
    count++;
    tmp = tmp &~ -tmp;  // clear lowest bit
  }

  TString hnTitle(name);
  const Int_t dim = count;
  Int_t nbins[dim];
  Double_t xmin[dim];
  Double_t xmax[dim];

  Int_t i=0;
  Int_t c=0;
  while(c<dim && i<32){
    if(entries&(1<<i)){
      Bool_t highres = res&(1<<i);
      const char* label;
      GetDimParams(i, highres, label, nbins[c], xmin[c], xmax[c]);
      hnTitle += Form(";%s",label);
      c++;
    }
      
    i++;
  }

  hnTitle += ";";

  THnSparseF* hntmp = new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);

  hntmp->Sumw2();

  return hntmp; 
}

//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFFMoments::GetDimParams(Int_t iEntry, Bool_t highres, const char* &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
 
  // 0 - 9
  //  0 |  1  |  2  |  3   |  4  |  5 |   6   |      7      |    8    |    9    |
  // vz | ntr | ep  | epb  |  z  | xi | lnjT  |  DeltaTheta | FFM_gen | FFM_rec |
  // 10 - 24
  //  0 |  1  |  2  |    3     |      4        |||      n
  // pt | eta | phi |   area   | nconstituents |||     gen  2x5+n
  // pt | eta | phi |   area   | nconstituents |||     rec  3x5+n
  // pt | eta | phi |          |    N(order)   |||    track 4x5+n
  // 25 - 32
  //	 25	|	26	|
  //  fraction | RDeltaPt |  R for Relative
  
  if( iEntry < 10 ) {
    switch(iEntry){
    case 0:
      label = "V0 centrality (%)";
      break;
        
    case 1:
      label = "vz";
      break;
        
    case 2:
      label = "nb. of input tracks";
      break;
        
    case 3:
      label = "event plane #psi";
      break;

    case 4:
      label = "z";
      break;
         
    case 5:
      label = "#xi";
      break;
         
    case 6: 
    if(fFFJtValue==0)  label = "ln(j_{T})";
    else if (fFFJtValue==1)  label = "1./(j_{T})";
      break;
         
    case 7:
      label = "#Delta #theta";
      break;

    case 8:
	  label = Form("(#frac{ N+1 }{2})^{%f} M_{N}^{gen}", fFFMScalePower);
      break;
         
    case 9:
	  label = Form("(#frac{ N+1 }{2})^{%f} M_{N}^{rec}", fFFMScalePower);
      break;
         
    }
  } else if( iEntry >= 10 && iEntry < 23) {
    switch(iEntry%5){
    case 0:
      label = "p_{T}";
      break;
         
    case 1:
      label = "#eta";
      break;

    case 2:
      label = "#phi";
      break;
         
    case 3:
      label = "area";
      break;
         
    case 4:
      label = "N_{constituents}";
      break;
    }
    switch(iEntry/5){
    case 2:
      label = Form("%s^{gen}",label);
      break;
    case 3:
      label = Form("%s^{rec}",label);
      break;
    case 4:
      label = Form("%s^{track}",label);
      break;
    }
  } else if( iEntry == 24) {
    label = "N";
  } else if( iEntry == 25) {
    label = "Match_fraction";
  } else if( iEntry == 26) {
	label = "(p^{Rec}_{T}-p^{Gen}_{T})";
  } else if( iEntry == 27) {
        label = "p_{T} Resolution";
  } else if (iEntry == 28) {
        label = "UE p_{T} density (GeV/c)"; 
  }
  nbins = fnBinsAxis[iEntry];
  xmin = fBinMinAxis[iEntry];
  xmax = fBinMaxAxis[iEntry];
  
  if(fkIsPbPb && iEntry == 2){
    nbins *= 5;
    xmax *= 10.;
  }

  if(highres && (!TString(label).Contains("N_{constituents}")) )  nbins *= 5;

}

//_______________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFFMoments::TrackQAFilter(TList* list, Int_t type)
{

  ////External Track Cut for chi2 in AOD
  ////reject all tracks if there is "the one" in the event
  ////see the difinition in the 7 th line after. 
  ////fExtTrackCutType = 0;// 1 for EtaChi2; 2 for PtChi2
  
  Bool_t Cut_on_track = kFALSE;
  if(type == kTrackAOD ){
    for ( Int_t i = 0 ; i< list->GetEntries(); i++) {
      AliAODTrack * fAODtrack = (AliAODTrack *) (list->At(i));
      if( fExtTrackCutType == 1) {
        Double_t fChi2AOD = fAODtrack->Chi2perNDF();
        if ( fChi2AOD > 1.8) continue;
        Double_t fEta = fAODtrack->Eta();
        // Chi2 vs Eta: Chi2 in [1.39-C, 1.89-C] or [0, 1.89-C ], C=(0.89^2 - eta^2 )^0.5
        Double_t fR = 0.90625;
        Double_t x0 = 1.90625;
        Double_t x =  x0 - TMath::Sqrt( fR*fR - fEta*fEta );
        if( fChi2AOD < x ) { Cut_on_track = kTRUE; break;}
      } else if (fExtTrackCutType == 2) {
        Double_t fChi2AOD = fAODtrack->Chi2perNDF();
        Double_t pt = fAODtrack->Pt();
        // Chi2 vs Pt: below line y/y0 + x/x0 =1 => y small than y = y0( 1- x/x0 ) => y <  y0( 1- x/x0 ); 
        Double_t x0 = 1.5;
        Double_t y0 = 20.;
        if( pt < y0 * ( 1. - fChi2AOD/x0 ) ) { Cut_on_track = kTRUE; break;}
      }
    }
    if( Cut_on_track ) { list->Clear(); } 
  }
  return Cut_on_track;
}

//______________________________________________________________________________________________________________________________________________________
void  AliAnalysisTaskJetFFMoments::FillTrackAssoList(AliAODJet* jet,TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, const TArrayI& indexAODTr, const TArrayS& isRefGen, TList* jetTrackListTR, Bool_t scaleStrangeness, TArrayD* sWeight, TList* listRecTracks)
{
  // fill objects for jet track reconstruction efficiency or secondaries contamination
  // jetTrackListTR pointer: track refs if not NULL

  Int_t nTracksInJet = jetTrackList->GetSize(); // list with AODMC tracks
  if(nTracksInJet == 0) return;
  sWeight->Reset(1.);
  listRecTracks->Clear();
  Int_t count = 0;

  for(Int_t iTr=0; iTr<nTracksInJet; iTr++){ // jet tracks loop

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (jetTrackList->At(iTr));
    if(!gentrack)continue;
    // find jet track in gen tracks list
    Int_t iGen = tracksGen->IndexOf(gentrack);
    if(iGen<0){
      if(fDebug>0) Printf("%s:%d gen jet track not found ",(char*)__FILE__,__LINE__);
      continue;
    }

    if(isRefGen[iGen] != 1) continue; // select primaries

    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // gen level acc & pt cuts - skip in case of track refs
    if(!jetTrackListTR && (etaGen < fTrackEtaMin || etaGen > fTrackEtaMax)) continue;
    if(!jetTrackListTR && (phiGen < fTrackPhiMin || phiGen > fTrackPhiMax)) continue;
    if(!jetTrackListTR &&  ptGen  < fTrackPtMin) continue;

    Int_t iRec   = indexAODTr[iGen]; // can be -1 if no good reconstructed track

    Bool_t isRec = (iRec>=0) ? kTRUE : kFALSE;

    Bool_t isJetTrack = kFALSE;
    if(!jetTrackListTR) isJetTrack = kTRUE; // skip trackRefs check for tracks in ideal cone

    if(isRec){

      AliAODTrack* rectrack = dynamic_cast<AliAODTrack*> (tracksRec->At(iRec));
      if(!rectrack) continue;

      if(jetTrackListTR){
        Int_t iRecTR = jetTrackListTR->IndexOf((AliVTrack*)rectrack);
        if(iRecTR >=0 ) isJetTrack = kTRUE; // rec tracks assigned to jet
      }

      if(isJetTrack){
        if(scaleStrangeness) sWeight->AddAt(GetMCStrangenessFactorCMS(gentrack),count);
	listRecTracks->Add(rectrack);
	count++;

      } // End isJetTrack
    } // End isRec
  } // End loop on tracks in jet
}

//_______________________________________________________________________________
TList* AliAnalysisTaskJetFFMoments::GetTrackInJetAssoList(AliAODJet* jet,  TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, TArrayI indexAODTr, const TArrayS& isRefGen, TList* jetTrackListTR, Bool_t scaleStrangeness)
{

    TList* trackInJetAsso = new TList();
    TArrayD* sWeight = new TArrayD(jetTrackList->GetSize());
    if (! TString(jetTrackList->GetName()).Contains("Gener")) {
    if(jetTrackListTR->GetSize()==0) {
      FillTrackAssoList(jet,jetTrackList,tracksGen, tracksRec, indexAODTr, isRefGen,0,scaleStrangeness, sWeight,trackInJetAsso);
    }
    else {
      FillTrackAssoList(jet,jetTrackList,tracksGen, tracksRec, indexAODTr, isRefGen,jetTrackListTR,scaleStrangeness, sWeight,trackInJetAsso);
    }
   }

  delete sWeight;
  return trackInJetAsso;
}

//_______________________________________________________________________________________________
Float_t AliAnalysisTaskJetFFMoments::ComputeFF(Bool_t SelectTrackAsso, AliAODJet * jet, TArrayI indexTr, TList* tracks, TArrayD* sWeight, TH2D** histo)
{
  Int_t nTracks = 0;
  if(SelectTrackAsso&&indexTr.GetSize()!=0) nTracks = indexTr.GetSize();
  else nTracks = tracks->GetSize();
  Float_t JetPt = -1.;
  Float_t SumPt = 0.;

   if(jet == NULL) {// no jet as input
    //	loop over reconstructed tracks, get SumPt
    for(Int_t i=0; i<nTracks; i++){
      if(SelectTrackAsso&&!fFFBckgMode){
	if(indexTr.At(i) == -1) continue;
      }
      AliVParticle* track;
      track = dynamic_cast<AliVParticle*>(tracks->At(i));
      if(!track) continue;
      SumPt += track->Pt();
    }
    if(SumPt <= 0) { printf("strange Pt_Sum <= 0."); /*return;*/ }
    else JetPt = SumPt;

    //properties of tracks
    for(Int_t i=0; i<nTracks; i++){
      AliVParticle* track = 0x0;
      if(SelectTrackAsso&&indexTr.GetSize()!=0&&!fFFBckgMode){
	if(indexTr.At(i) == -1) continue;
      }
       track = dynamic_cast<AliVParticle*>(tracks->At(i));
      if(!track) continue;
      Float_t pt = track->Pt();
      Float_t eta = track->Eta();
      Float_t phi = track->Phi();
      Float_t z = pt/JetPt - 1.e-6;
      Float_t xi = log(1./z);
      Float_t deltaTheta = -999.;
      Float_t lnjT = -999.;
      Double_t TracksInJets[7]= {  z , xi , lnjT , deltaTheta , pt , eta , phi };
      for( Int_t iAxis = 0; iAxis < 7; iAxis++) histo[iAxis]->Fill(TracksInJets[iAxis], SumPt);
    }
  } // End if(jet==NULL)
  else {
    TVector3 jetV(jet->Px(), jet->Py(), jet->Pz());
    JetPt = jet->Pt();
    TVector3 trackV;
    if(tracks){
      for(Int_t i=0; i<nTracks; i++){
	if(SelectTrackAsso&&indexTr.GetSize()!=0&&!fFFBckgMode) if(indexTr.At(i) == -1) continue;
	AliVParticle* track = 0x0;
        track = dynamic_cast<AliVParticle*>(tracks->At(i));
	if(!track) continue;
	trackV.SetXYZ(track->Px(), track->Py(), track->Pz());
	FillFF(trackV,jetV,sWeight->At(i),histo);
      }
    }
    else {
      TRefArray* trackRef = jet->GetRefTracks();
      for( Int_t icons = 0; icons<trackRef->GetEntries(); icons++){
	AliAODTrack* track = (AliAODTrack*) trackRef->At(icons);
	trackV.SetXYZ(track->Px(), track->Py(), track->Pz());
	FillFF(trackV,jetV,sWeight->At(icons),histo);
      }
    }
  } // End else
  return JetPt;
}

//_______________________________________________________________________________________________
void AliAnalysisTaskJetFFMoments::FillFF(TVector3 trackV, TVector3 jetV, Double_t sWeight, TH2D** histo)
  {
    Float_t Z, xi, lnjT, deltaTheta, pt, eta, phi;
    if(jetV.Pt() <= 0.) {
      Z = -1.; xi = -1.; lnjT = -1.; deltaTheta = -1.; pt = -1.; eta = -999.; phi = -999.;
    }
    Z = trackV.Pt()/jetV.Pt() -1e-6;
    xi = log(1./Z);
    if(fFFJtValue==0 && jetV.Mag() && (trackV.Cross(jetV).Mag())) lnjT = TMath::Log((trackV.Cross(jetV)).Mag()/jetV.Mag()); //log(jT)
    if(fFFJtValue==1 && jetV.Mag() && (trackV.Cross(jetV)).Mag()) lnjT = 1./((trackV.Cross(jetV)).Mag()/jetV.Mag()); //1/jT
    deltaTheta = trackV.Angle(jetV);
    pt = trackV.Pt();
    eta = trackV.Eta();
    phi = trackV.Phi();
    if(phi <= 0.) phi += TMath::TwoPi();
    Double_t TracksInJets[7]= {  Z , xi , lnjT , deltaTheta , pt , eta , phi };
     for( Int_t iAxis = 0; iAxis < 7; iAxis++) histo[iAxis]->Fill(TracksInJets[iAxis], jetV.Pt(),sWeight);
  }

//__________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFFMoments::GetMCStrangenessFactorCMS(AliAODMCParticle* daughter)
{
  // strangeness ratio MC/data as function of mother pt from CMS (7TeV data) in |eta|<2.0

  TClonesArray *tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!tca) return 1;

  AliAODMCParticle* currentMother   = daughter;
  AliAODMCParticle* currentDaughter = daughter;


  // find first primary mother K0s, Lambda or Xi
  while(1){

    Int_t daughterPDG   = currentDaughter->GetPdgCode();

    Int_t motherLabel   = currentDaughter->GetMother();
    if(motherLabel >= tca->GetEntriesFast()){ // protection
      currentMother = currentDaughter;
      break;
    }

    currentMother     = (AliAODMCParticle*) tca->At(motherLabel);

    if(!currentMother){
      currentMother = currentDaughter;
      break;
    }

    Int_t motherPDG   = currentMother->GetPdgCode();

    // phys. primary found ?
    if(currentMother->IsPhysicalPrimary()) break;

    if(TMath::Abs(daughterPDG) == 321){ // K+/K- e.g. from phi (ref data not feeddown corrected)
      currentMother = currentDaughter; break;
    }
    if(TMath::Abs(motherPDG) == 310 ){ // K0s e.g. from phi (ref data not feeddown corrected)
      break;
    }
    if(TMath::Abs(motherPDG) == 3212 && TMath::Abs(daughterPDG) == 3122){ // mother Sigma0, daughter Lambda (this case not included in feeddown corr.)
      currentMother = currentDaughter; break;
    }

    currentDaughter = currentMother;
  }


  Int_t motherPDG   = currentMother->GetPdgCode();
  Double_t motherPt = currentMother->Pt();

  Double_t fac = 1;

  if(TMath::Abs(motherPDG) == 310 || TMath::Abs(motherPDG)==321){ // K0s / K+ / K-

    if(0.00 <= motherPt && motherPt < 0.20) fac = 0.768049;
    else if(0.20 <= motherPt && motherPt < 0.40) fac = 0.732933;
    else if(0.40 <= motherPt && motherPt < 0.60) fac = 0.650298;
    else if(0.60 <= motherPt && motherPt < 0.80) fac = 0.571332;
    else if(0.80 <= motherPt && motherPt < 1.00) fac = 0.518734;
    else if(1.00 <= motherPt && motherPt < 1.20) fac = 0.492543;
    else if(1.20 <= motherPt && motherPt < 1.40) fac = 0.482704;
    else if(1.40 <= motherPt && motherPt < 1.60) fac = 0.488056;
    else if(1.60 <= motherPt && motherPt < 1.80) fac = 0.488861;
    else if(1.80 <= motherPt && motherPt < 2.00) fac = 0.492862;
    else if(2.00 <= motherPt && motherPt < 2.20) fac = 0.504332;
    else if(2.20 <= motherPt && motherPt < 2.40) fac = 0.501858;
    else if(2.40 <= motherPt && motherPt < 2.60) fac = 0.512970;
    else if(2.60 <= motherPt && motherPt < 2.80) fac = 0.524131;
    else if(2.80 <= motherPt && motherPt < 3.00) fac = 0.539130;
    else if(3.00 <= motherPt && motherPt < 3.20) fac = 0.554101;
    else if(3.20 <= motherPt && motherPt < 3.40) fac = 0.560348;
    else if(3.40 <= motherPt && motherPt < 3.60) fac = 0.568869;
    else if(3.60 <= motherPt && motherPt < 3.80) fac = 0.583310;
    else if(3.80 <= motherPt && motherPt < 4.00) fac = 0.604818;
    else if(4.00 <= motherPt && motherPt < 5.00) fac = 0.632630;
    else if(5.00 <= motherPt && motherPt < 6.00) fac = 0.710070;
    else if(6.00 <= motherPt && motherPt < 8.00) fac = 0.736365;
    else if(8.00 <= motherPt && motherPt < 10.00) fac = 0.835865;
  }

  if(TMath::Abs(motherPDG) == 3122){ // Lambda

    if(0.00 <= motherPt && motherPt < 0.20) fac = 0.645162;
    else if(0.20 <= motherPt && motherPt < 0.40) fac = 0.627431;
    else if(0.40 <= motherPt && motherPt < 0.60) fac = 0.457136;
    else if(0.60 <= motherPt && motherPt < 0.80) fac = 0.384369;
    else if(0.80 <= motherPt && motherPt < 1.00) fac = 0.330597;
    else if(1.00 <= motherPt && motherPt < 1.20) fac = 0.309571;
    else if(1.20 <= motherPt && motherPt < 1.40) fac = 0.293620;
    else if(1.40 <= motherPt && motherPt < 1.60) fac = 0.283709;
    else if(1.60 <= motherPt && motherPt < 1.80) fac = 0.282047;
    /////else if(1.80 <= motherPt && motherPt < 2.00) fac = 0.277261;
    else if(2.00 <= motherPt && motherPt < 2.20) fac = 0.275772;
    else if(2.20 <= motherPt && motherPt < 2.40) fac = 0.280726;
    else if(2.40 <= motherPt && motherPt < 2.60) fac = 0.288540;
    else if(2.60 <= motherPt && motherPt < 2.80) fac = 0.288315;
    else if(2.80 <= motherPt && motherPt < 3.00) fac = 0.296619;
    else if(3.00 <= motherPt && motherPt < 3.20) fac = 0.302993;
    else if(3.20 <= motherPt && motherPt < 3.40) fac = 0.338121;
   else if(3.40 <= motherPt && motherPt < 3.60) fac = 0.349800;
    else if(3.60 <= motherPt && motherPt < 3.80) fac = 0.356802;
    else if(3.80 <= motherPt && motherPt < 4.00) fac = 0.391202;
    else if(4.00 <= motherPt && motherPt < 5.00) fac = 0.422573;
    else if(5.00 <= motherPt && motherPt < 6.00) fac = 0.573815;
    else if(6.00 <= motherPt && motherPt < 8.00) fac = 0.786984;
    else if(8.00 <= motherPt && motherPt < 10.00) fac = 1.020021;
  }

  if(TMath::Abs(motherPDG) == 3312 || TMath::Abs(motherPDG) == 3322){ // xi

    if(0.00 <= motherPt && motherPt < 0.20) fac = 0.666620;
    else if(0.20 <= motherPt && motherPt < 0.40) fac = 0.575908;
    else if(0.40 <= motherPt && motherPt < 0.60) fac = 0.433198;
    else if(0.60 <= motherPt && motherPt < 0.80) fac = 0.340901;
    else if(0.80 <= motherPt && motherPt < 1.00) fac = 0.290896;
    else if(1.00 <= motherPt && motherPt < 1.20) fac = 0.236074;
    else if(1.20 <= motherPt && motherPt < 1.40) fac = 0.218681;
    else if(1.40 <= motherPt && motherPt < 1.60) fac = 0.207763;
    else if(1.60 <= motherPt && motherPt < 1.80) fac = 0.222848;
    else if(1.80 <= motherPt && motherPt < 2.00) fac = 0.208806;
    else if(2.00 <= motherPt && motherPt < 2.20) fac = 0.197275;
    else if(2.20 <= motherPt && motherPt < 2.40) fac = 0.183645;
    else if(2.40 <= motherPt && motherPt < 2.60) fac = 0.188788;
    else if(2.60 <= motherPt && motherPt < 2.80) fac = 0.188282;
    else if(2.80 <= motherPt && motherPt < 3.00) fac = 0.207442;
    else if(3.00 <= motherPt && motherPt < 3.20) fac = 0.240388;
    else if(3.20 <= motherPt && motherPt < 3.40) fac = 0.241916;
    else if(3.40 <= motherPt && motherPt < 3.60) fac = 0.208276;
    else if(3.60 <= motherPt && motherPt < 3.80) fac = 0.234550;
    else if(3.80 <= motherPt && motherPt < 4.00) fac = 0.251689;
    else if(4.00 <= motherPt && motherPt < 5.00) fac = 0.310204;
    else if(5.00 <= motherPt && motherPt < 6.00) fac = 0.343492;
  }

  Double_t weight = 1;
  if(fac > 0) weight = 1/fac;

  return weight;
}

// _____________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFFMoments::CalcSingleTrackEff(TH1D* trackEffPtGen, TH2D* trackEffEtaPhiGen, TH1D* trackEffPtRec, TH2D* trackEffEtaPhiRec,
						    TList* tracksGen, const TArrayI& indexAODTr, const TArrayS& isRefGen, Bool_t scaleStrangeness)
{

  // single track reconstruction efficiency

  Int_t nTracksGen  = tracksGen->GetSize();

  if(!nTracksGen || !fAOD) return;

  for(Int_t iGen=0; iGen<nTracksGen; iGen++){

    if(isRefGen[iGen] != 1) continue; // select primaries

    AliAODMCParticle* gentrack =  dynamic_cast<AliAODMCParticle*> (tracksGen->At(iGen));
    if(!gentrack) continue;
    Double_t ptGen  = gentrack->Pt();
    Double_t etaGen = gentrack->Eta();
    Double_t phiGen = TVector2::Phi_0_2pi(gentrack->Phi());

    // apply same acc & pt cuts as for FF

    if(etaGen < fTrackEtaMin || etaGen > fTrackEtaMax) continue;
    if(phiGen < fTrackPhiMin || phiGen > fTrackPhiMax) continue;
    if(ptGen  < fTrackPtMin) continue;

    if(trackEffPtGen) trackEffPtGen->Fill(ptGen);
    if(trackEffEtaPhiGen) trackEffEtaPhiGen->Fill(etaGen, phiGen);

    Int_t iRec = indexAODTr[iGen]; // can be -1 if no good reconstructed track

    if(iRec>=0 && trackEffPtRec && trackEffEtaPhiRec){
      if(scaleStrangeness){
        Double_t weight = GetMCStrangenessFactorCMS(gentrack);
        trackEffPtRec->Fill(ptGen, weight);
        trackEffEtaPhiRec->Fill(etaGen, phiGen, weight);
      }
      else {
       trackEffPtRec->Fill(ptGen);
       trackEffEtaPhiRec->Fill(etaGen, phiGen);
    }
    }
  }
}

// _____________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFFMoments::GetTracksTiltedwrpJetAxis(Float_t alpha, TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius,Double_t& sumPt)
{
  // List of tracks in cone perpendicular to the jet azimuthal direction

  Double_t jetMom[3];
  jet->PxPyPz(jetMom);

  TVector3 jet3mom(jetMom);
  // Rotate phi and keep eta unchanged
  Double_t etaTilted = jet3mom.Eta();
  Double_t phiTilted = TVector2::Phi_0_2pi(jet3mom.Phi()) + alpha;
  if(phiTilted > 2*TMath::Pi()) phiTilted = phiTilted - 2*TMath::Pi();
  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){

  AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
  if(!track)continue;
  Double_t trackMom[3];
  track->PxPyPz(trackMom);
  TVector3 track3mom(trackMom);

  Double_t deta = track3mom.Eta() - etaTilted;
  Double_t dphi = TMath::Abs(track3mom.Phi() - phiTilted);
  if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
  Double_t dR = TMath::Sqrt(deta * deta + dphi * dphi);

    if(dR<=radius){
     if (outputlist) outputlist->Add(track);
      sumPt += track->Pt();
    }
  }
}

// _____________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFFMoments::IsOutlier(AliGenPythiaEventHeader * const header) {
     Bool_t hasOutlier = kFALSE;
     if(!header) return hasOutlier;

    if ( fPtHardAndPythiaJetPtFactor > 0.) {
     Float_t pbuf[4];
     TLorentzVector jetvec;
     for(int ijet = 0; ijet < header->NTriggerJets(); ijet++){
      memset(pbuf, 0, sizeof(Float_t) * 4);
       header->TriggerJet(ijet, pbuf);
       jetvec.SetPxPyPzE(pbuf[0], pbuf[1], pbuf[2], pbuf[3]);
       fh2PtHardVsPt[0]->Fill(header->GetPtHard(),jetvec.Pt());
       if(TMath::Abs(jetvec.Pt()) >= fPtHardAndPythiaJetPtFactor * header->GetPtHard()){
       fh2PtHardVsPtCut[0]->Fill(header->GetPtHard(),jetvec.Pt());
        if(fDebug)  Printf("Reject : pythia jet %2.2f, factor %2.2f, ptHard %f", jetvec.Pt(), fPtHardAndPythiaJetPtFactor, header->GetPtHard());
        hasOutlier = true;
        break;
      }
    }
   }

    if ( fPtHardAndTrackPtFactor > 0.) {
       Int_t ntracks = 0;
       TClonesArray *tca = 0x0;

      if(fAOD) {
       tca = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
       ntracks = tca->GetEntriesFast();
      } else {
       AliMCEvent* mcEvent = MCEvent();
       ntracks = mcEvent->GetNumberOfTracks();
      }

      if (!ntracks) return hasOutlier;
       for(int it = 0;it < ntracks ;++it){
       AliVParticle* part = 0x0;
       if(tca) {
       part =  dynamic_cast<AliVParticle*> (tca->At(it));
       } else {
       part =  dynamic_cast<AliVParticle*> (MCEvent()->GetTrack(it));
       }
       Float_t trackpt = part->Pt();

       fh2PtHardVsPt[1]->Fill(header->GetPtHard(),trackpt);
       if (TMath::Abs(trackpt) >= (fPtHardAndTrackPtFactor * header->GetPtHard())) {
        fh2PtHardVsPtCut[1]->Fill(header->GetPtHard(),trackpt);
        if(fDebug) Printf("Reject : track %2.2f, factor %2.2f, ptHard %f", trackpt, fPtHardAndTrackPtFactor, header->GetPtHard());
        hasOutlier = true;
        break;
        }
      }
     }
     return hasOutlier;
}

// _____________________________________________________________________________________________________________________________________________________________________
AliGenPythiaEventHeader *AliAnalysisTaskJetFFMoments::GetPythiaHeader()  {
     if(!MCEvent()) return 0x0;
     AliGenPythiaEventHeader *pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
     if(fDebug>10) pythiaHeader->Dump();
     if (!pythiaHeader) {
       // Check if AOD
       if(fAOD) {
       AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::StdBranchName()));
       if (aodMCH) {
         for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
           pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
           if (pythiaHeader) break;
        }
       }
     }
    }
    return pythiaHeader;
}

// __________________________________________________________________________________________________________________________________________________________$
AliGenHerwigEventHeader *AliAnalysisTaskJetFFMoments::GetHerwigHeader()  {
     if(!MCEvent()) return 0x0;
     AliGenHerwigEventHeader *HerwigHeader = dynamic_cast<AliGenHerwigEventHeader*>(MCEvent()->GenEventHeader());
     if(fDebug>10) HerwigHeader->Dump();
     if (!HerwigHeader) {
       // Check if AOD
       if(fAOD) {
       AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::StdBranchName()));
       if (aodMCH) {
         for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
           HerwigHeader = dynamic_cast<AliGenHerwigEventHeader*>(aodMCH->GetCocktailHeader(i));
           if (HerwigHeader) break;
        }
       }
     }
    }
    return HerwigHeader;
}
