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
#include <TFormula.h>

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
  fkRequireVZEROAC(kFALSE),
  fkRequireTZEROvtx(kFALSE),
  fCentCutUp(0),
  fCentCutLo(0),
  fVtxZMax(8),
  fVtxR2Max(1),
  fFilterMask(0),
  fFilterMaskBestPt(0),
  fFilterType(0),
  fkUseAODTrackInput(kTRUE),
  fkUseAODMCInput(kFALSE),
  fTrackPtMin(0.15),
  fTrackPtMax(100.),
  fExtTrackCutType(0),
  fkUseHFcuts(kFALSE),
  fkRequireITSRefit(kFALSE),
  fkApplySharedClusterCut(kFALSE),
  fTrackEtaMin(-0.9),
  fTrackEtaMax(0.9),
  fTCAJetsOut(0x0),
  fAODJetBackgroundOut(0x0),
  fNonStdBranch(""),
  fNonStdFile(""),
  fAnaJetType("ALLJET"),
  fJetPtMin(2.),
  fJetEtaMin(-0.5),
  fJetEtaMax(0.5),
  fJetDeltaPhiCut(TMath::Pi()/3.),
  fkDoJetMatching(kFALSE),
  fkFillMismatchHisto(kFALSE),
  fJetMatchingFractionMin(0.5),
  fJetMatchedDistMax(0.3),
  fkUseJetFromInput(kFALSE),
  fNUsedJets(999),
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
  fh1vZSelect(0x0)
{
  //
  // Default Constructor
  //

  fTrackType[0] = kTrackUndef;
  fTrackType[1] = kTrackUndef;
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) { fh1GenTracks[iAxis] = 0x0; fh1RecTracks[iAxis] = 0x0;}
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
    for( Int_t iAxis = 0; iAxis < 5; iAxis++) fh1JetPr_Mismatched[iJetBranch][iAxis] = 0x0;
    for( Int_t iAxis = 0; iAxis < 4; iAxis++) fh2Constituents[iJetBranch][iAxis] =0x0;
    for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh1TracksInJets[iJetBranch][iAxis] =0x0;
    fhnJetFFM_Raw[iJetBranch] = 0x0;
    fhnJetFFM_Sub[iJetBranch] = 0x0;
    fhnJetFFM_Imp[iJetBranch] = 0x0;
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
  //     25    |
  //  fraction |
  //
  fnBinsAxis[0] = 20;          fBinMinAxis[0] = 0.;            fBinMaxAxis[0] = 100.;
  fnBinsAxis[1] = 80;          fBinMinAxis[1] = 0.;            fBinMaxAxis[1] = 400.;
  fnBinsAxis[2] = 40;          fBinMinAxis[2] = -20.;          fBinMaxAxis[2] = 20.;
  fnBinsAxis[3] = 30;          fBinMinAxis[3] = 0.;            fBinMaxAxis[3] = pi;
  fnBinsAxis[4] = 44;          fBinMinAxis[4] = 0.;            fBinMaxAxis[4] = 1.1;
  fnBinsAxis[5] = 28;          fBinMinAxis[5] = 0.;            fBinMaxAxis[5] = 7.;
  fnBinsAxis[6] = 25;          fBinMinAxis[6] = 4.-2*pi;       fBinMaxAxis[6] = 4.;
  fnBinsAxis[7] = 20;          fBinMinAxis[7] = 0.;            fBinMaxAxis[7] = 0.5;
  fnBinsAxis[8] =150;          fBinMinAxis[8] = 0.;            fBinMaxAxis[8] = 30.;
  fnBinsAxis[9] =150;          fBinMinAxis[9] = 0.;            fBinMaxAxis[9] = 30.;
  for( Int_t i = 0; i < 3; i ++) { 
    fnBinsAxis[10+i*5] = 40;   fBinMinAxis[10+i*5] = 0.;       fBinMaxAxis[10+i*5] = 200.;
    fnBinsAxis[11+i*5] = 80;   fBinMinAxis[11+i*5] = -1.;      fBinMaxAxis[11+i*5] = 1.;
    fnBinsAxis[12+i*5] = 90;   fBinMinAxis[12+i*5] = 0.;       fBinMaxAxis[12+i*5] = 2*pi;
    if( i != 2) { fnBinsAxis[13+i*5] = 20;     fBinMinAxis[13+i*5] = 0.;       fBinMaxAxis[13+i*5] = 1.;}
    if( i != 2) { fnBinsAxis[14+i*5] = 50;     fBinMinAxis[14+i*5] = 0.5;      fBinMaxAxis[14+i*5] = 50.5;}
  }
  fnBinsAxis[23] = 0;  fBinMinAxis[23] = 0.;           fBinMaxAxis[23] = 0.;
  fnBinsAxis[24] = fFFMMomMax; fBinMinAxis[24] = nAxisMin;    fBinMaxAxis[24] = nAxisMax;
  fnBinsAxis[25] = 21; fBinMinAxis[23] = 0.;           fBinMaxAxis[25] = 1.05;
  for( Int_t i = 26; i < 32; i ++) {
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
  fkRequireVZEROAC(kFALSE),
  fkRequireTZEROvtx(kFALSE),
  fCentCutUp(0),
  fCentCutLo(0),
  fVtxZMax(8),
  fVtxR2Max(1),
  fFilterMask(0),
  fFilterMaskBestPt(0),
  fFilterType(0),
  fkUseAODTrackInput(kTRUE),
  fkUseAODMCInput(kFALSE),
  fTrackPtMin(0.15),
  fTrackPtMax(100.),
  fExtTrackCutType(0),
  fkUseHFcuts(kFALSE),
  fkRequireITSRefit(kFALSE),
  fkApplySharedClusterCut(kFALSE),
  fTrackEtaMin(-0.9),
  fTrackEtaMax(0.9),
  fTCAJetsOut(0x0),
  fAODJetBackgroundOut(0x0),
  fNonStdBranch(""),
  fNonStdFile(""),
  fAnaJetType("ALLJET"),
  fJetPtMin(2.),
  fJetEtaMin(-0.5),
  fJetEtaMax(0.5),
  fJetDeltaPhiCut(TMath::Pi()/3.),
  fkDoJetMatching(kFALSE),
  fkFillMismatchHisto(kFALSE),
  fJetMatchingFractionMin(0.5),
  fJetMatchedDistMax(0.3),
  fkUseJetFromInput(kFALSE),
  fNUsedJets(999),
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
  fh1vZSelect(0x0)
{
  //
  // named ctor
  //

  fTrackType[0] = kTrackUndef;
  fTrackType[1] = kTrackUndef;
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) { fh1GenTracks[iAxis] = 0x0; fh1RecTracks[iAxis] = 0x0;}
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
    for( Int_t iAxis = 0; iAxis < 5; iAxis++) fh1JetPr_Mismatched[iJetBranch][iAxis] = 0x0;
    for( Int_t iAxis = 0; iAxis < 4; iAxis++) fh2Constituents[iJetBranch][iAxis] =0x0;
    for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh1TracksInJets[iJetBranch][iAxis] =0x0;
    fhnJetFFM_Raw[iJetBranch] = 0x0;
    fhnJetFFM_Sub[iJetBranch] = 0x0;
    fhnJetFFM_Imp[iJetBranch] = 0x0;
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
  //     25    |
  //  fraction |
  //
  fnBinsAxis[0] = 20;          fBinMinAxis[0] = 0.;            fBinMaxAxis[0] = 100.;
  fnBinsAxis[1] = 80;          fBinMinAxis[1] = 0.;            fBinMaxAxis[1] = 400.;
  fnBinsAxis[2] = 40;          fBinMinAxis[2] = -20.;          fBinMaxAxis[2] = 20.;
  fnBinsAxis[3] = 30;          fBinMinAxis[3] = 0.;            fBinMaxAxis[3] = pi;
  fnBinsAxis[4] = 44;          fBinMinAxis[4] = 0.;            fBinMaxAxis[4] = 1.1;
  fnBinsAxis[5] = 28;          fBinMinAxis[5] = 0.;            fBinMaxAxis[5] = 7.;
  fnBinsAxis[6] = 25;          fBinMinAxis[6] = 4.-2*pi;       fBinMaxAxis[6] = 4.;
  fnBinsAxis[7] = 20;          fBinMinAxis[7] = 0.;            fBinMaxAxis[7] = 0.5;
  fnBinsAxis[8] =150;          fBinMinAxis[8] = 0.;            fBinMaxAxis[8] = 30.;
  fnBinsAxis[9] =150;          fBinMinAxis[9] = 0.;            fBinMaxAxis[9] = 30.;
  for( Int_t i = 0; i < 3; i ++) {
    fnBinsAxis[10+i*5] = 40;   fBinMinAxis[10+i*5] = 0.;       fBinMaxAxis[10+i*5] = 200.;
    fnBinsAxis[11+i*5] = 80;   fBinMinAxis[11+i*5] = -1.;      fBinMaxAxis[11+i*5] = 1.;
    fnBinsAxis[12+i*5] = 90;   fBinMinAxis[12+i*5] = 0.;       fBinMaxAxis[12+i*5] = 2*pi;
    if( i != 2) { fnBinsAxis[13+i*5] = 20;     fBinMinAxis[13+i*5] = 0.;       fBinMaxAxis[13+i*5] = 1.;}
    if( i != 2) { fnBinsAxis[14+i*5] = 50;     fBinMinAxis[14+i*5] = 0.5;      fBinMaxAxis[14+i*5] = 50.5;}
  }
  fnBinsAxis[23] = 0;  fBinMinAxis[23] = 0.;           fBinMaxAxis[23] = 0.;
  fnBinsAxis[24] = fFFMMomMax; fBinMinAxis[24] = nAxisMin;    fBinMaxAxis[24] = nAxisMax;
  fnBinsAxis[25] = 21; fBinMinAxis[23] = 0.;           fBinMaxAxis[25] = 1.05;
  for( Int_t i = 26; i < 32; i ++) {
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
    printf("Type of jet analysed: %s\n", fAnaJetType.Data());
    printf("TRACK PARAMETERS:   FilterMask: %d	 pt_min: %.3f GeV/c	eta min max: [%.2f,%.2f]  extra cut type: %d\n", fFilterMask, fTrackPtMin, fTrackEtaMin, fTrackEtaMax, fExtTrackCutType); 
    printf("JET   PARAMETERS:   Algorithm: %d	 R: %.2f                eta min max: [%.2f,%.2f]  pt_min: %.3f \n", fAlgorithm, fRparam, fJetEtaMin, fJetEtaMax,fJetPtMin);
    if(fkDoJetMatching) printf("Jet matching ON with parameters, distance: %.2f, energy fraction: %.3f\n", fJetMatchedDistMax, fJetMatchingFractionMin);
    printf("FFM computed with scaling power %d, background mode %d and parameters: %f %f %.3f\n", fFFMScalePower, fFFMBckgType, fFFMBckgPar1, fFFMBckgPar2, fFFMBckgMu);
  }
  if( fkDoJetMatching == kFALSE) {
    fJetMatchedDistMax = 0.; fJetMatchingFractionMin = 0.;
  }

  if((!fJetBranch[0].Contains("MC") && !fkDoJetReco) || (fkDoJetReco && !(fTrackType[0] == kTrackAODMCCharged || fTrackType[0] ==  kTrackAODMCChargedAcceptance  || fTrackType[0] == kTrackAODMCextra || fTrackType[0] == kTrackAODMCextraonly))){
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

    if( fJetBranch[0].Length() != 0 ) {//Kine used
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
      cout << "Centrality class: " << cenClass << endl;
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


  if(!selectEvent){
    PostData(1, fHistList);
    return;
  }
  fh1CentralitySelect->Fill(cent);  
  fh1vZSelect->Fill(zVtx);

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
  fastjet::JetDefinition jet_def_for_rho(fBkgAlgorithm, fRparam);
  fastjet::GhostedAreaSpec ghostSpec(fGhostEtaMax, fActiveAreaRepeats, fGhostArea);
  fastjet::AreaDefinition area_def(fAreaType, 
				   fastjet::GhostedAreaSpec(fastjet::SelectorRapRange(-fGhostEtaMax, fGhostEtaMax)));

  // NB explicit ghosts do not work for moments with N < 0 with FastJet versions < 3.1  
  //  AreaDefinition area_def(fastjet::active_area_explicit_ghosts, 
  //                          fastjet::GhostedAreaSpec(fastjet::SelectorAbsRapMax(5.0)));
  //
  fastjet::Selector sel_jets = fastjet::SelectorNHardest(fNUsedJets) * fastjet::SelectorEtaRange(fJetEtaMin, fJetEtaMax);
 
  //---------------------- get all particles and jets ----------------------- 
  TList ParticleList[2];
  Int_t ifirstBr = 0;
  if((!fJetBranch[0].Contains("MC") && !fkDoJetReco) || (fkDoJetReco && !(fTrackType[0] == kTrackAODMCCharged || fTrackType[0] ==  kTrackAODMCChargedAcceptance  || fTrackType[0] == kTrackAODMCextra || fTrackType[0] == kTrackAODMCextraonly))) { 
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
	if( iJetBranch == 1) fh1RecTracks[iAxis]->Fill(fAllParticle[iAxis]);
      }

      // Add particles to fastjet: used for jet finding and/or FFM
      // Carefull energy is not well determined in real data, should not matter for p_T scheme?
      fastjet::PseudoJet particle(vp->Px(), vp->Py(), vp->Pz(), vp->E());
      particle.set_user_index(i);
      full_event[iJetBranch].push_back(particle);

    } // End loop on list of particles

    
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
      else {  aodJets[iJetBranch] = dynamic_cast<TClonesArray*>(AODEvent()->FindListObject(fJetBranch[iJetBranch].Data())); 
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
	if ( SelectJet(tmp_jet)) {
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
    GetListOfMatchedJets(fListJets,nUsedJets,*&(fListMatchedJets[0]),*&(fListMatchedJets[1]),ifirstBr);
    listUsedJets[0] = fListMatchedJets[0];
    listUsedJets[1] = fListMatchedJets[1];
  } // End kDoJetMatching
  else {
    listUsedJets[0] = fListJets[0];
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
    fastjet::Selector rho_range;
    if(fFFMBckgType==kDoughnut) rho_range = fastjet::SelectorDoughnut(fFFMBckgPar1,fFFMBckgPar2); // Rmin, Rmax
    else if(fFFMBckgType==kRectangle) rho_range = fastjet::SelectorRectangle(fFFMBckgPar1,fFFMBckgPar2); // HalfRapwidth, HalfPhiWidth
    else if(fFFMBckgType==kEtaRange) rho_range = fastjet::SelectorEtaRange(fFFMBckgPar1,fFFMBckgPar2); // EtaMin, EtaMax
    else {Printf("Please check the background methods available to compute ffm !");}
    fastjet::JetMedianBackgroundEstimator bge(rho_range, jet_def_for_rho, area_def);
    fastjet::Subtractor subtractor(&bge);
    bge.set_particles(full_event[iJetBranch]);	// tell the background estimator what to use

    // Initialize JetFFMoments
    fastjet::contrib::JetFFMoments ffms_unsubtracted(fFFMNMin, fFFMNMax, nFFBins);
    fastjet::contrib::JetFFMoments ffms_subtracted(fFFMNMin, fFFMNMax, nFFBins, &bge);
    fastjet::contrib::JetFFMoments ffms_improved(fFFMNMin, fFFMNMax, nFFBins, &bge);
    double mu = fFFMBckgMu; // 25: typical value in the 150-200 GeV range
    ffms_improved.set_improved_subtraction(mu, rho_range, full_event[iJetBranch], jet_def_for_rho, area_def);
    // With FastJet-3.1, this could be done this way:
    // ffms_improved.set_improved_subtraction(mu);

    // Clear pseudo jet list before to calculate FFM
    if(sorted_jets_fj[iJetBranch].size()) sorted_jets_fj[iJetBranch].clear();
    for(Int_t ijet = 0; ijet < usedDim; ijet++ ) {
      AliAODJet* ujet = ((AliAODJet*)listUsedJets[iJetBranch]->At(ijet));
      if (!ujet) continue;
      if(!ClassifyJetEvent(ijet,usedDim)) break;

      //*********************************************
      // Transverse and longitudinal jet info
      //*********************************************
      TVector3 jetV(ujet->Px(), ujet->Py(), ujet->Pz());
      TRefArray* trackRef = ujet->GetRefTracks();
      Double_t Z = -1.;
      Double_t lnjT = -1.;
      Double_t deltaTheta = -1.;
      for( Int_t icons = 0; icons<trackRef->GetEntries(); icons++){
	AliAODTrack* track = (AliAODTrack*) trackRef->At(icons);
	TVector3 trackV(track->Px(), track->Py(), track->Pz());
	if(ujet->Pt()) Z = track->Pt()/ujet->Pt();
	else {cout << "The reconstructed jet pt is NULL!! Check your jet list containt!";}
	lnjT = TMath::Log((trackV.Cross(jetV)).Mag()/jetV.Mag());
	deltaTheta = trackV.Angle(jetV);
	Double_t TracksInJets[10] = {
	  Z-1.e-6, log(1./Z), lnjT, deltaTheta, track->Pt(), track->Eta(), track->Phi(), ujet->GetNEF()
	};
	for( Int_t iAxis = 0; iAxis < 4; iAxis++) fh2Constituents[iJetBranch][iAxis]->Fill(TracksInJets[iAxis], ujet->Pt());
	for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh1TracksInJets[iJetBranch][iAxis]->Fill(TracksInJets[iAxis+4]);
      } // End loop on tracks in jets

      
      //*********************************************
      // MC <---> Reco Jet info
      //*********************************************
      if(iJetBranch==1){  
	AliAODJet* urecjet = (AliAODJet*)listUsedJets[iJetBranch]->At(ijet);
	AliAODJet* ugenjet = 0x0;
	if(listUsedJets[0]->GetEntries()) {ugenjet = (AliAODJet*)listUsedJets[0]->At(ijet);}
	else { ugenjet = (AliAODJet*)listUsedJets[iJetBranch]->At(ijet);}

	Double_t jetEntriesMatch[10] = {
	  ugenjet->Pt(), ugenjet->Eta(), ugenjet->Phi(), ugenjet->EffectiveAreaCharged(), (Double_t)( (TRefArray *) (ugenjet->GetRefTracks()))->GetEntries(),
	  urecjet->Pt(), urecjet->Eta(), urecjet->Phi(), urecjet->EffectiveAreaCharged(), (Double_t)( (TRefArray *) (urecjet->GetRefTracks()))->GetEntries(),
	};
	for(Int_t iAxis=0; iAxis<5; iAxis++) fh2MatchedJets[iAxis]->Fill(jetEntriesMatch[iAxis], jetEntriesMatch[iAxis+5]);
      } // End ijetBranch == 1 
    
      //*********************************************
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
      
      //*********************************************
      // Compute FFM
      //*********************************************
      fastjet::PseudoJet pseudoJetUsed;
      int ret = AliAODJetToPseudoJet(((AliAODJet*)listUsedJets[iJetBranch]->At(ijet)), pseudoJetUsed);
      if(ret) return;
      sorted_jets_fj[iJetBranch].push_back(pseudoJetUsed);

      if(fDebug > 9) { 
	fastjet::PseudoJet j = sorted_jets_fj[iJetBranch][ijet];
        cout << "pseudo jet eta: " << j.eta() << ", phi: " << j.phi() << ", pt: " << j.pt() << ", area: " << j.area() << endl;
        fastjet::PseudoJet area4vect = j.area_4vector();
        cout << "a_x: " << area4vect.px() << ", a_y: " << area4vect.py() << "a_z: " << area4vect.pz() << endl;

	fastjet::PseudoJet subtracted_jet = subtractor(j);
	cout << "# [subtracted hard jet " << ijet+1 << "]:(pt,y,phi) = (" << subtracted_jet.pt() << ", "
	     << subtracted_jet.rap() << ", " << subtracted_jet.phi() << ")" << endl;
      }
            
      // Compute FFM
      vector<double> ffm_unsubtracted = ffms_unsubtracted(pseudoJetUsed);
      vector<double> ffm_subtracted   = ffms_subtracted(pseudoJetUsed);
      vector<double> ffm_improved     = ffms_improved(pseudoJetUsed);
      
      //
      // Fill FFM histograms and ThnSparses
      //
      for(Int_t iN = 0; iN < fFFMMomMax; iN++) {
        Double_t A = TMath::Power( 0.5*(ffms_subtracted.N(iN) + 1), fFFMScalePower);
        fFFMs_Raw[iJetBranch][iN][ijet] = A*ffm_unsubtracted[iN];
	fFFMs_Sub[iJetBranch][iN][ijet] = A*ffm_subtracted[iN];
	fFFMs_Imp[iJetBranch][iN][ijet] = A*ffm_improved[iN];

        Double_t N = iN * (fFFMNMax - fFFMNMin)/(fFFMMomMax-1) + fFFMNMin;
        Double_t FFM_Raw_in_pt[3] = { fFFMs_Raw[iJetBranch][iN][ijet], pseudoJetUsed.pt(), N };
        fhnJetFFM_Raw[iJetBranch]->Fill( FFM_Raw_in_pt );
	Double_t FFM_Sub_in_pt[3] = { fFFMs_Sub[iJetBranch][iN][ijet], pseudoJetUsed.pt(), N };
	fhnJetFFM_Sub[iJetBranch]->Fill( FFM_Sub_in_pt );
	Double_t FFM_Imp_in_pt[3] = { fFFMs_Imp[iJetBranch][iN][ijet], pseudoJetUsed.pt(), N };
	fhnJetFFM_Imp[iJetBranch]->Fill( FFM_Imp_in_pt );

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
        cout << "# Fragmentation function moments for full jet "
             << ": (pt, eta, phi) = (" << pseudoJetUsed.pt() << ", "<< pseudoJetUsed.eta() << ", " << pseudoJetUsed.phi() << "), #constituents="<< pseudoJetUsed.constituents().size() << endl
             <<": (px, py, pz, E )  = ("<<pseudoJetUsed.px()<<", "<<pseudoJetUsed.py()<<", "<<pseudoJetUsed.pz()<<", "<<pseudoJetUsed.e()<<")"<<"  sqrt(px^2+py^2): "<<TMath::Sqrt(pseudoJetUsed.px()*pseudoJetUsed.px() + pseudoJetUsed.py()*pseudoJetUsed.py())<<endl;
        cout << "# N	M_N(unsubtracted)	M_N(subtracted)	      M_N(improved)(NOT YET)" << endl;
        for (unsigned int in=0; in<ffm_subtracted.size(); in++){
	  printf("%1.2f	   %2.9f  	  %2.9f  	  %2.9f	 \n", ffms_subtracted.N(in), ffm_unsubtracted[in], ffm_subtracted[in], ffm_improved[in]);
        }
        cout << "# N	A*M_N(unsubtracted)	A*M_N(subtracted)	      A*M_N(improved)(NOT YET)" << endl;
        for (unsigned int in=0; in<ffm_subtracted.size(); in++){
	  printf("%1.2f	   %2.9f  	  %2.9f  	  %2.9f  \n", ffms_subtracted.N(in), fFFMs_Raw[iJetBranch][in][ijet], fFFMs_Sub[iJetBranch][in][ijet], fFFMs_Imp[iJetBranch][in][ijet]);
        }
      } // End fDebug
    } // End loop on ijet
    
  } // End loop on jet branches

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
Int_t  AliAnalysisTaskJetFFMoments::GetListOfTracks(TList *list,Int_t type)
{

  //
  // get list of tracks/particles for different types
  //

  if(fDebug>2) Printf("%s:%d Selecting tracks with %d",(char*)__FILE__,__LINE__,type);

  Int_t iCount = 0;
  if(type==kTrackAOD || type==kTrackAODextra || type==kTrackAODextraonly || type==kTrackAODMCextra || type==kTrackAODMCextraonly){

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
	if(TMath::Abs(tr->Eta())>fTrackEtaMax){
	  if(fDebug>10)Printf("%s:%d Not matching eta %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	  continue;
	}
	if(tr->Pt()<fTrackPtMin){
	  if(fDebug>10)Printf("%s:%d Not matching pt %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	  continue;
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


	if(TMath::Abs(trackAOD->Eta())>fTrackEtaMax) continue;
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
	
	if(TMath::Abs(track->Eta())>fTrackEtaMax) continue;
	if(track->Pt()<fTrackPtMin) continue;
	list->Add(track);

	iCount++;
      }
    }
    
  }
  else if (type ==  kTrackKineAll||type == kTrackKineCharged){
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
      else if(type == kTrackKineCharged){
	if(part->Particle()->GetPDG()->Charge()==0)continue;
	if(part->Pt()<fTrackPtMin)continue;
	list->Add(part);
	iCount++;
      }
    }
  }
  else if (type == kTrackAODMCCharged || type == kTrackAODMCAll || type == kTrackAODMCChargedAcceptance) {
    AliAODEvent *aod = 0;
    if(fkUseAODMCInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else aod = AODEvent();
    if(!aod)return iCount;
    TClonesArray *tca = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    for(int it = 0;it < tca->GetEntriesFast();++it){
      AliAODMCParticle *part = (AliAODMCParticle*)(tca->At(it));
      if(!part->IsPhysicalPrimary())continue;
      if(type == kTrackAODMCAll){
	if(part->Pt()<fTrackPtMin)continue;
	list->Add(part);
	iCount++;
      }
      else if (type == kTrackAODMCCharged || type == kTrackAODMCChargedAcceptance ){
	if(part->Charge()==0)continue;
	if(part->Pt()<fTrackPtMin)continue;
	if(kTrackAODMCCharged){
	  list->Add(part);
	}
	else {
	  if(TMath::Abs(part->Eta())>fTrackEtaMax)continue;
	  list->Add(part);
	}
	iCount++;
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
			
	if(part->Pt()<fTrackPtMin)	continue;
	if(TMath::Abs(part->Eta())>fTrackEtaMax)	continue;
	if(part->Charge()==0)	continue;
			
	if((part->Pt()>=fTrackPtMin) && (TMath::Abs(part->Eta())<=fTrackEtaMax) && (part->Charge()!=0))list->Add(part);
	iCount++;
      }
    }
  }
  
  list->Sort();
  return iCount;
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
	  // After keeping the mismatched information, remove the unmatched jets from the jet list
	} // End tmpjet
      } // End Fill mismatched histo
    } // End loop on jets
  } // End loop on jet branches

  return listJets[1]->GetEntries();

}

// _______________________________________________________________________________
void AliAnalysisTaskJetFFMoments::CreateHistos() 
{
  if (fDebug != 0) printf("AliAnalysisTaskJetFFMoments::CreateHistos() start!\n");

  fPrParameters = new TProfile("h1Parameters","Parameters for analysis", 23, 0.5, 23.5);
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
  //     25    |
  //  fraction |
  //
  sHisName = TString("h1GenTracks");
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
    fh1GenTracks[iAxis] = CreateTH1D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, fkHighResolution);
  }
  sHisName = TString("h1RecTracks");
  for( Int_t iAxis = 0; iAxis < 3; iAxis++ ) {
    fh1RecTracks[iAxis] = CreateTH1D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, fkHighResolution);
  }

  for(int iJetBranch=0; iJetBranch < fgkFFMNJetBranches; iJetBranch++){
    if( iJetBranch == 0) fAnaObjectType = "Gen";
    else { fAnaObjectType = "Rec"; }

    sHisName = Form("h1MisMatched_%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
    for( Int_t iAxis = 0; iAxis < 5; iAxis++) fh1JetPr_Mismatched[iJetBranch][iAxis] = CreateTH1D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), iJetBranch?15+iAxis:10+iAxis, fkHighResolution);

    sHisName = Form("h2StructureOf%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
    for( Int_t iAxis = 0; iAxis < 4; iAxis++) fh2Constituents[iJetBranch][iAxis] = CreateTH2D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis+5].Data()), 4+iAxis, iJetBranch?15:10, fkHighResolution);
    sHisName = Form("h1TracksIn%s_%s", fAnaJetType.Data(), fAnaObjectType.Data());
    for( Int_t iAxis = 0; iAxis < 3; iAxis++) fh1TracksInJets[iJetBranch][iAxis] = CreateTH1D(Form("%s_%s", sHisName.Data(), sAxisName[iAxis].Data()), 20+iAxis, fkHighResolution);

    if( iJetBranch == 1) {
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
    fhnJetFFM_Raw[iJetBranch] = CreateTHnSparseF(Form("hnFFM_Raw_%s_%s", fAnaJetType.Data(), fAnaObjectType.Data()), entries, res);
    fhnJetFFM_Sub[iJetBranch] = CreateTHnSparseF(Form("hnFFM_Sub_%s_%s", fAnaJetType.Data(), fAnaObjectType.Data()), entries, res);
    fhnJetFFM_Imp[iJetBranch] = CreateTHnSparseF(Form("hnFFM_Imp_%s_%s", fAnaJetType.Data(), fAnaObjectType.Data()), entries, res);

    if( fAnaJetType.Contains("DIJET")) {
      fh1Asy_DiJets[iJetBranch]  = new TH1F(Form("h1Asy_DiJets_%s", fAnaObjectType.Data()), Form("Asymmetry DiJets_%s", fAnaObjectType.Data()), 20, 0.0 , 1.0);
      fh1Asy_DiJets[iJetBranch]->GetXaxis()->SetTitle("Asymmetry = #frac{p_{T}^1-p_{T}^2}{p_{T}^1+p_{T}}^2");
    }

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
    if( fAnaJetType.Contains("DIJET")) fHistListJets[iJetBranch]->Add(fh1Asy_DiJets[iJetBranch]);
    for( Int_t iAxis = 0; iAxis < 4; iAxis++) fHistListJets[iJetBranch]->Add(fh2Constituents[iJetBranch][iAxis]);
    for( Int_t iAxis = 0; iAxis < 3; iAxis++) fHistListJets[iJetBranch]->Add(fh1TracksInJets[iJetBranch][iAxis]);
    fHistListJets[iJetBranch]->Add(fhnJetFFM_Raw[iJetBranch]);
    fHistListJets[iJetBranch]->Add(fhnJetFFM_Sub[iJetBranch]);
    fHistListJets[iJetBranch]->Add(fhnJetFFM_Imp[iJetBranch]);
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
    if ( SelectJet(&(aodjet)) ) {
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
    fastjet::PseudoJet inputPart(track->Px(),track->Py(),track->Pz(),track->E());  // create PseudoJet object
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
Bool_t AliAnalysisTaskJetFFMoments::SelectJet(AliAODJet * jet)
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

  //tracks in jet can not lager than 100 GeV/c
  TRefArray * tracks = jet->GetRefTracks();
  Int_t ntracks = tracks -> GetEntriesFast();
  for(Int_t i=0; i < ntracks; i++) { 
    AliAODTrack * track = (AliAODTrack *) (jet->GetTrack(i));
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
  //     25    |
  //  fraction |
  
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
      label = "ln(j_{T})";
      break;
         
    case 7:
      label = "#Delta #theta";
      break;

    case 8:
      label = "FFM_M_{N}^{gen}";
      break;
         
    case 9:
      label = "FFM_M_{N}^{rec}";
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
