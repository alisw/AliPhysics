// **************************************
// Task used for the correction of determiantion of reconstructed jet spectra
// Compares input (gen) and output (rec) jets   
// *******************************************


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

 
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TRefArray.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include  "TDatabasePDG.h"
#include <TGrid.h>

#include "AliAnalysisTaskJetCluster.h"
#include "AliAnalysisManager.h"
#include "AliJetFinder.h"
#include "AliJetHeader.h"
#include "AliJetReader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetKineReaderHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
// get info on how fastjet was configured
#include "fastjet/config.h"

using std::vector;

ClassImp(AliAnalysisTaskJetCluster)

AliAnalysisTaskJetCluster::~AliAnalysisTaskJetCluster(){
  //
  // Destructor
  //

  delete fRef;
  delete fRandom;

  if(fTCAJetsOut)fTCAJetsOut->Delete();
  delete fTCAJetsOut;
  
  if(fTCAJetsOutRan)fTCAJetsOutRan->Delete();
  delete fTCAJetsOutRan;
  
  if(fTCARandomConesOut)fTCARandomConesOut->Delete();
  delete fTCARandomConesOut;
  
  if(fTCARandomConesOutRan)fTCARandomConesOutRan->Delete();
  delete fTCARandomConesOutRan;
  
  if(fAODJetBackgroundOut)fAODJetBackgroundOut->Reset();
  delete fAODJetBackgroundOut;
}

AliAnalysisTaskJetCluster::AliAnalysisTaskJetCluster(): 
  AliAnalysisTaskSE(),
  fAOD(0x0),
  fAODExtension(0x0),
  fRef(new TRefArray),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseBackgroundCalc(kFALSE),
  fEventSelection(kFALSE),     
  fRequireVZEROAC(kFALSE),     
  fRequireTZEROvtx(kFALSE),
  fFilterMask(0),
  fFilterMaskBestPt(0),
  fFilterType(0),
  fJetTypes(kJet),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),  
  fNSkipLeadingRan(0),
  fNSkipLeadingCone(0),
  fNRandomCones(0),
  fAvgTrials(1),
  fExternalWeight(1),
  fTrackEtaWindow(0.9),    
  fRecEtaWindow(0.5),
  fTrackPtCut(0.),							
  fJetOutputMinPt(0.150),
  fMaxTrackPtInJet(100.),
  fJetTriggerPtCut(0),
  fVtxZCut(8),
  fVtxR2Cut(1),
  fCentCutUp(0),
  fCentCutLo(0),
  fStoreRhoLeadingTrackCorr(kFALSE),
  fNonStdBranch(""),
  fBackgroundBranch(""),
  fNonStdFile(""),
  fMomResH1(0x0),
  fMomResH2(0x0),
  fMomResH3(0x0),
  fMomResH1Fit(0x0),
  fMomResH2Fit(0x0),
  fMomResH3Fit(0x0),
  fhEffH1(0x0),
  fhEffH2(0x0),
  fhEffH3(0x0),
  fUseTrPtResolutionSmearing(kFALSE),
  fUseDiceEfficiency(kFALSE),
  fDiceEfficiencyMinPt(-1.),
  fUseTrPtResolutionFromOADB(kFALSE),
  fUseTrEfficiencyFromOADB(kFALSE),
  fPathTrPtResolution(""),
  fPathTrEfficiency(""),
  fChangeEfficiencyFraction(0.),
  fEfficiencyFixed(1.),
  fRparam(1.0), 
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area), 
  fGhostArea(0.01),
  fActiveAreaRepeats(1),
  fGhostEtamax(1.5),
  fTCAJetsOut(0x0),
  fTCAJetsOutRan(0x0),
  fTCARandomConesOut(0x0),
  fTCARandomConesOutRan(0x0),
  fAODJetBackgroundOut(0x0),
  fRandom(0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NJetsRec(0x0),
  fh1NConstRec(0x0),
  fh1NConstLeadingRec(0x0),
  fh1PtJetsRecIn(0x0),
  fh1PtJetsLeadingRecIn(0x0),
  fh1PtJetConstRec(0x0),
  fh1PtJetConstLeadingRec(0x0),
  fh1PtTracksRecIn(0x0),
  fh1PtTracksLeadingRecIn(0x0),
  fh1NJetsRecRan(0x0),
  fh1NConstRecRan(0x0),
  fh1PtJetsLeadingRecInRan(0x0),
  fh1NConstLeadingRecRan(0x0),
  fh1PtJetsRecInRan(0x0),
  fh1PtTracksGenIn(0x0),
  fh1Nch(0x0),
  fh1CentralityPhySel(0x0), 
  fh1Centrality(0x0), 
  fh1CentralitySelect(0x0),
  fh1ZPhySel(0x0), 
  fh1Z(0x0), 
  fh1ZSelect(0x0),
  fh2NRecJetsPt(0x0),
  fh2NRecTracksPt(0x0),
  fh2NConstPt(0x0),
  fh2NConstLeadingPt(0x0),
  fh2JetPhiEta(0x0),
  fh2LeadingJetPhiEta(0x0),
  fh2JetEtaPt(0x0),
  fh2LeadingJetEtaPt(0x0),
  fh2TrackEtaPt(0x0),
  fh2LeadingTrackEtaPt(0x0),
  fh2JetsLeadingPhiEta(0x0),
  fh2JetsLeadingPhiPt(0x0),
  fh2TracksLeadingPhiEta(0x0),
  fh2TracksLeadingPhiPt(0x0),
  fh2TracksLeadingJetPhiPt(0x0),
  fh2JetsLeadingPhiPtW(0x0),
  fh2TracksLeadingPhiPtW(0x0),
  fh2TracksLeadingJetPhiPtW(0x0),
  fh2NRecJetsPtRan(0x0),
  fh2NConstPtRan(0x0),
  fh2NConstLeadingPtRan(0x0),
  fh2PtNch(0x0),
  fh2PtNchRan(0x0),
  fh2PtNchN(0x0),
  fh2PtNchNRan(0x0),
  fh2TracksLeadingJetPhiPtRan(0x0),
  fh2TracksLeadingJetPhiPtWRan(0x0),
  fh3CentvsRhoLeadingTrackPt(0x0),
  fh3CentvsSigmaLeadingTrackPt(0x0),
  fh3MultvsRhoLeadingTrackPt(0x0),
  fh3MultvsSigmaLeadingTrackPt(0x0),
  fh3CentvsRhoLeadingTrackPtQ1(0x0),
  fh3CentvsRhoLeadingTrackPtQ2(0x0),
  fh3CentvsRhoLeadingTrackPtQ3(0x0),
  fh3CentvsRhoLeadingTrackPtQ4(0x0),
  fh3CentvsSigmaLeadingTrackPtQ1(0x0),
  fh3CentvsSigmaLeadingTrackPtQ2(0x0),
  fh3CentvsSigmaLeadingTrackPtQ3(0x0),
  fh3CentvsSigmaLeadingTrackPtQ4(0x0),
  fh3MultvsRhoLeadingTrackPtQ1(0x0),
  fh3MultvsRhoLeadingTrackPtQ2(0x0),
  fh3MultvsRhoLeadingTrackPtQ3(0x0),
  fh3MultvsRhoLeadingTrackPtQ4(0x0),
  fh3MultvsSigmaLeadingTrackPtQ1(0x0),
  fh3MultvsSigmaLeadingTrackPtQ2(0x0),
  fh3MultvsSigmaLeadingTrackPtQ3(0x0),
  fh3MultvsSigmaLeadingTrackPtQ4(0x0),
  fh3CentvsDeltaRhoLeadingTrackPtQ1(0x0),
  fh3CentvsDeltaRhoLeadingTrackPtQ2(0x0),
  fh3CentvsDeltaRhoLeadingTrackPtQ3(0x0),
  fh3CentvsDeltaRhoLeadingTrackPtQ4(0x0),
  fh2PtGenPtSmeared(0x0),
  fp1Efficiency(0x0),
  fp1PtResolution(0x0),
  fHistList(0x0)  
{
  //
  // Constructor
  //

  for(int i = 0;i<3;i++){
    fh1BiARandomCones[i] = 0;
    fh1BiARandomConesRan[i] = 0;
  }
  for(int i = 0;i<kMaxCent;i++){
    fh2JetsLeadingPhiPtC[i] = 0;     
    fh2JetsLeadingPhiPtWC[i] = 0;      //! jet correlation with leading jet
    fh2TracksLeadingJetPhiPtC[i] = 0;
    fh2TracksLeadingJetPhiPtWC[i] = 0;
  }
}

AliAnalysisTaskJetCluster::AliAnalysisTaskJetCluster(const char* name):
  AliAnalysisTaskSE(name),
  fAOD(0x0),
  fAODExtension(0x0),
  fRef(new TRefArray),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseBackgroundCalc(kFALSE),
  fEventSelection(kFALSE),
  fRequireVZEROAC(kFALSE),     
  fRequireTZEROvtx(kFALSE), 
  fFilterMask(0),
  fFilterMaskBestPt(0),
  fFilterType(0),
  fJetTypes(kJet),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fNSkipLeadingRan(0),
  fNSkipLeadingCone(0),
  fNRandomCones(0),
  fAvgTrials(1),
  fExternalWeight(1),    
  fTrackEtaWindow(0.9),    
  fRecEtaWindow(0.5),
  fTrackPtCut(0.),							
  fJetOutputMinPt(0.150),
  fMaxTrackPtInJet(100.),
  fJetTriggerPtCut(0),
  fVtxZCut(8),
  fVtxR2Cut(1),
  fCentCutUp(0),
  fCentCutLo(0),
  fStoreRhoLeadingTrackCorr(kFALSE),
  fNonStdBranch(""),
  fBackgroundBranch(""),
  fNonStdFile(""),
  fMomResH1(0x0),
  fMomResH2(0x0),
  fMomResH3(0x0),
  fMomResH1Fit(0x0),
  fMomResH2Fit(0x0),
  fMomResH3Fit(0x0),
  fhEffH1(0x0),
  fhEffH2(0x0),
  fhEffH3(0x0),
  fUseTrPtResolutionSmearing(kFALSE),
  fUseDiceEfficiency(kFALSE),
  fDiceEfficiencyMinPt(-1.),
  fUseTrPtResolutionFromOADB(kFALSE),
  fUseTrEfficiencyFromOADB(kFALSE),
  fPathTrPtResolution(""),
  fPathTrEfficiency(""),
  fChangeEfficiencyFraction(0.),
  fEfficiencyFixed(1.),
  fRparam(1.0), 
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area), 
  fGhostArea(0.01),
  fActiveAreaRepeats(1),
  fGhostEtamax(1.5),
  fTCAJetsOut(0x0),
  fTCAJetsOutRan(0x0),
  fTCARandomConesOut(0x0),
  fTCARandomConesOutRan(0x0),
  fAODJetBackgroundOut(0x0),
  fRandom(0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NJetsRec(0x0),
  fh1NConstRec(0x0),
  fh1NConstLeadingRec(0x0),
  fh1PtJetsRecIn(0x0),
  fh1PtJetsLeadingRecIn(0x0),
  fh1PtJetConstRec(0x0),
  fh1PtJetConstLeadingRec(0x0),
  fh1PtTracksRecIn(0x0),
  fh1PtTracksLeadingRecIn(0x0),
  fh1NJetsRecRan(0x0),
  fh1NConstRecRan(0x0),
  fh1PtJetsLeadingRecInRan(0x0),
  fh1NConstLeadingRecRan(0x0),
  fh1PtJetsRecInRan(0x0),
  fh1PtTracksGenIn(0x0),
  fh1Nch(0x0),
  fh1CentralityPhySel(0x0), 
  fh1Centrality(0x0), 
  fh1CentralitySelect(0x0),
  fh1ZPhySel(0x0), 
  fh1Z(0x0), 
  fh1ZSelect(0x0),
  fh2NRecJetsPt(0x0),
  fh2NRecTracksPt(0x0),
  fh2NConstPt(0x0),
  fh2NConstLeadingPt(0x0),
  fh2JetPhiEta(0x0),
  fh2LeadingJetPhiEta(0x0),
  fh2JetEtaPt(0x0),
  fh2LeadingJetEtaPt(0x0),
  fh2TrackEtaPt(0x0),
  fh2LeadingTrackEtaPt(0x0),
  fh2JetsLeadingPhiEta(0x0),
  fh2JetsLeadingPhiPt(0x0),
  fh2TracksLeadingPhiEta(0x0),
  fh2TracksLeadingPhiPt(0x0),
  fh2TracksLeadingJetPhiPt(0x0),
  fh2JetsLeadingPhiPtW(0x0),
  fh2TracksLeadingPhiPtW(0x0),
  fh2TracksLeadingJetPhiPtW(0x0),
  fh2NRecJetsPtRan(0x0),
  fh2NConstPtRan(0x0),
  fh2NConstLeadingPtRan(0x0),
  fh2PtNch(0x0),
  fh2PtNchRan(0x0),
  fh2PtNchN(0x0),
  fh2PtNchNRan(0x0),
  fh2TracksLeadingJetPhiPtRan(0x0),
  fh2TracksLeadingJetPhiPtWRan(0x0),
  fh3CentvsRhoLeadingTrackPt(0x0),
  fh3CentvsSigmaLeadingTrackPt(0x0),
  fh3MultvsRhoLeadingTrackPt(0x0),
  fh3MultvsSigmaLeadingTrackPt(0x0),
  fh3CentvsRhoLeadingTrackPtQ1(0x0),
  fh3CentvsRhoLeadingTrackPtQ2(0x0),
  fh3CentvsRhoLeadingTrackPtQ3(0x0),
  fh3CentvsRhoLeadingTrackPtQ4(0x0),
  fh3CentvsSigmaLeadingTrackPtQ1(0x0),
  fh3CentvsSigmaLeadingTrackPtQ2(0x0),
  fh3CentvsSigmaLeadingTrackPtQ3(0x0),
  fh3CentvsSigmaLeadingTrackPtQ4(0x0),
  fh3MultvsRhoLeadingTrackPtQ1(0x0),
  fh3MultvsRhoLeadingTrackPtQ2(0x0),
  fh3MultvsRhoLeadingTrackPtQ3(0x0),
  fh3MultvsRhoLeadingTrackPtQ4(0x0),
  fh3MultvsSigmaLeadingTrackPtQ1(0x0),
  fh3MultvsSigmaLeadingTrackPtQ2(0x0),
  fh3MultvsSigmaLeadingTrackPtQ3(0x0),
  fh3MultvsSigmaLeadingTrackPtQ4(0x0),
  fh3CentvsDeltaRhoLeadingTrackPtQ1(0x0),
  fh3CentvsDeltaRhoLeadingTrackPtQ2(0x0),
  fh3CentvsDeltaRhoLeadingTrackPtQ3(0x0),
  fh3CentvsDeltaRhoLeadingTrackPtQ4(0x0),
  fh2PtGenPtSmeared(0x0),
  fp1Efficiency(0x0),
  fp1PtResolution(0x0),
  fHistList(0x0)
{
  //
  // named ctor
  //

  for(int i = 0;i<3;i++){
    fh1BiARandomCones[i] = 0;
    fh1BiARandomConesRan[i] = 0;
  }
  for(int i = 0;i<kMaxCent;i++){
    fh2JetsLeadingPhiPtC[i] = 0;     
    fh2JetsLeadingPhiPtWC[i] = 0;      //! jet correlation with leading jet
    fh2TracksLeadingJetPhiPtC[i] = 0;
    fh2TracksLeadingJetPhiPtWC[i] = 0;
  }
  DefineOutput(1, TList::Class());  
}



Bool_t AliAnalysisTaskJetCluster::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  return kTRUE;
}

void AliAnalysisTaskJetCluster::UserCreateOutputObjects()
{

  //
  // Create the output container
  //

  fRandom = new TRandom3(0);


  // Connect the AOD


  if (fDebug > 1) printf("AnalysisTaskJetCluster::UserCreateOutputObjects() \n");

  

  if(fNonStdBranch.Length()!=0)
    {
      // only create the output branch if we have a name
      // Create a new branch for jets...
      //  -> cleared in the UserExec....
      // here we can also have the case that the brnaches are written to a separate file
      
      if(fJetTypes&kJet){
	fTCAJetsOut = new TClonesArray("AliAODJet", 0);
	fTCAJetsOut->SetName(fNonStdBranch.Data());
	AddAODBranch("TClonesArray",&fTCAJetsOut,fNonStdFile.Data());
      }

      if(fJetTypes&kJetRan){
	fTCAJetsOutRan = new TClonesArray("AliAODJet", 0);
	fTCAJetsOutRan->SetName(Form("%s_%s",fNonStdBranch.Data(),"random"));
	if(fUseDiceEfficiency || fUseTrPtResolutionSmearing) {
	  if(  fEfficiencyFixed < 1.)
	    fTCAJetsOutRan->SetName(Form("%s_%sDetector%d%dEffFixed%d",fNonStdBranch.Data(),"random",fUseTrPtResolutionSmearing,fUseDiceEfficiency,(int)(fEfficiencyFixed*100.)));
	  else
	    fTCAJetsOutRan->SetName(Form("%s_%sDetector%d%dFr%d",fNonStdBranch.Data(),"random",fUseTrPtResolutionSmearing,fUseDiceEfficiency,(int)(fChangeEfficiencyFraction*100.)));
	}
	AddAODBranch("TClonesArray",&fTCAJetsOutRan,fNonStdFile.Data());
      }

      if(fUseBackgroundCalc){
	if(!AODEvent()->FindListObject(Form("%s_%s",AliAODJetEventBackground::StdBranchName(),fNonStdBranch.Data()))){
	  fAODJetBackgroundOut = new AliAODJetEventBackground();
	  fAODJetBackgroundOut->SetName(Form("%s_%s",AliAODJetEventBackground::StdBranchName(),fNonStdBranch.Data()));
	  if(fUseDiceEfficiency || fUseTrPtResolutionSmearing) {
	    if(  fEfficiencyFixed < 1.)
	      fAODJetBackgroundOut->SetName(Form("%s_%sDetector%d%dEffFixed%d",AliAODJetEventBackground::StdBranchName(),fNonStdBranch.Data(),fUseTrPtResolutionSmearing,fUseDiceEfficiency,(int)(fEfficiencyFixed*100.)));
	    else
	      fAODJetBackgroundOut->SetName(Form("%s_%sDetector%d%dFr%d",AliAODJetEventBackground::StdBranchName(),fNonStdBranch.Data(),fUseTrPtResolutionSmearing,fUseDiceEfficiency,(int)(fChangeEfficiencyFraction*100.)));
	  }
	  AddAODBranch("AliAODJetEventBackground",&fAODJetBackgroundOut,fNonStdFile.Data());  
	}
      }
      // create the branch for the random cones with the same R 
      TString cName = Form("%sRandomConeSkip%02d",fNonStdBranch.Data(),fNSkipLeadingCone);
      if(fUseDiceEfficiency || fUseTrPtResolutionSmearing) {
	if(  fEfficiencyFixed < 1.)
	  cName = Form("%sDetector%d%dEffFixed%d_RandomConeSkip%02d",fNonStdBranch.Data(),fUseTrPtResolutionSmearing,fUseDiceEfficiency,(int)(fEfficiencyFixed*100.),fNSkipLeadingCone);
	else
	  cName = Form("%sDetector%d%dFr%d_RandomConeSkip%02d",fNonStdBranch.Data(),fUseTrPtResolutionSmearing,fUseDiceEfficiency,(int)(fChangeEfficiencyFraction*100.),fNSkipLeadingCone);
      }
      if(fNRandomCones>0){
	if(fJetTypes&kRC){
	  if(!AODEvent()->FindListObject(cName.Data())){
	    fTCARandomConesOut = new TClonesArray("AliAODJet", 0);
	    fTCARandomConesOut->SetName(cName.Data());
	    AddAODBranch("TClonesArray",&fTCARandomConesOut,fNonStdFile.Data());
	  }
	}
	// create the branch with the random for the random cones on the random event
	if(fJetTypes&kRCRan){
	  cName = Form("%sRandomCone_random",fNonStdBranch.Data());
	  if(!AODEvent()->FindListObject(cName.Data())){
	    fTCARandomConesOutRan = new TClonesArray("AliAODJet", 0);
	    fTCARandomConesOutRan->SetName(cName.Data());
	    AddAODBranch("TClonesArray",&fTCARandomConesOutRan,fNonStdFile.Data());
	  }
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


  if(!fHistList)fHistList = new TList();
  fHistList->SetOwner();
  PostData(1, fHistList); // post data in any case once

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //
  //  Histogram
    
  const Int_t nBinPt = 100;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = 0.0;
    }
    else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.0;
    }
  }
  
  const Int_t nBinPhi = 90;
  Double_t binLimitsPhi[nBinPhi+1];
  for(Int_t iPhi = 0;iPhi<=nBinPhi;iPhi++){
    if(iPhi==0){
      binLimitsPhi[iPhi] = -1.*TMath::Pi();
    }
    else{
      binLimitsPhi[iPhi] = binLimitsPhi[iPhi-1] + 1/(Float_t)nBinPhi * TMath::Pi()*2;
    }
  }



  const Int_t nBinEta = 40;
  Double_t binLimitsEta[nBinEta+1];
  for(Int_t iEta = 0;iEta<=nBinEta;iEta++){
    if(iEta==0){
      binLimitsEta[iEta] = -2.0;
    }
    else{
      binLimitsEta[iEta] = binLimitsEta[iEta-1] + 0.1;
    }
  }

  const int nChMax = 5000;

  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");

  fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");


  fh1NJetsRec = new TH1F("fh1NJetsRec","N reconstructed jets",120,-0.5,119.5);
  fh1NJetsRecRan = new TH1F("fh1NJetsRecRan","N reconstructed jets",120,-0.5,119.5);

  fh1NConstRec = new TH1F("fh1NConstRec","# jet constituents",120,-0.5,119.5);
  fh1NConstRecRan = new TH1F("fh1NConstRecRan","# jet constituents",120,-0.5,119.5);
  fh1NConstLeadingRec = new TH1F("fh1NConstLeadingRec","jet constituents",120,-0.5,119.5);
  fh1NConstLeadingRecRan = new TH1F("fh1NConstLeadingRecRan","jet constituents",120,-0.5,119.5);


  fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);

  fh1PtJetsRecIn  = new TH1F("fh1PtJetsRecIn","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetsRecInRan  = new TH1F("fh1PtJetsRecInRan","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetsLeadingRecIn = new TH1F("fh1PtJetsLeadingRecIn","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetsLeadingRecInRan = new TH1F("fh1PtJetsLeadingRecInRan","Rec jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetConstRec = new TH1F("fh1PtJetsConstRec","Rec jets constituents P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtJetConstLeadingRec = new TH1F("fh1PtJetsConstLeadingRec","Rec jets constituents P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
  fh1PtTracksRecIn  = new TH1F("fh1PtTracksRecIn",Form("Rec tracks P_T #eta < %1.2f;p_{T} (GeV/c)",fTrackEtaWindow),nBinPt,binLimitsPt);
  fh1PtTracksLeadingRecIn  = new TH1F("fh1PtTracksLeadingRecIn",Form("Rec tracks P_T #eta < %1.2f ;p_{T} (GeV/c)",fTrackEtaWindow),nBinPt,binLimitsPt);
  fh1PtTracksGenIn  = new TH1F("fh1PtTracksGenIn",Form("gen tracks P_T #eta < %1.2f ;p_{T} (GeV/c)",fTrackEtaWindow),nBinPt,binLimitsPt);
  fh1Nch = new TH1F("fh1Nch","charged multiplicity; N_{ch}",nChMax,-0.5,nChMax-0.5);

  fh1Centrality = new TH1F("fh1Centrality",";cent (%)",111,-0.5,110.5);
  fh1CentralitySelect = new TH1F("fh1CentralitySelect",";cent (%)",111,-0.5,110.5);
  fh1CentralityPhySel = new TH1F("fh1CentralityPhySel",";cent (%)",111,-0.5,110.5);

  fh1Z = new TH1F("fh1Z",";zvtx",100,-25,25);
  fh1ZSelect = new TH1F("fh1ZSelect",";zvtx",100,-25,25);
  fh1ZPhySel = new TH1F("fh1ZPhySel",";zvtx",100,-25,25);

  fh2NRecJetsPt = new TH2F("fh2NRecJetsPt","Number of jets above threshhold;p_{T,cut} (GeV/c);N_{jets}",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NRecJetsPtRan = new TH2F("fh2NRecJetsPtRan","Number of jets above threshhold;p_{T,cut} (GeV/c);N_{jets}",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NRecTracksPt = new TH2F("fh2NRecTracksPt","Number of tracks above threshhold;p_{T,cut} (GeV/c);N_{tracks}",nBinPt,binLimitsPt,50,-0.5,49.5);
  // 


  fh2NConstPt = new TH2F("fh2NConstPt","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NConstLeadingPt = new TH2F("fh2NConstLeadingPt","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NConstPtRan = new TH2F("fh2NConstPtRan","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);
  fh2NConstLeadingPtRan = new TH2F("fh2NConstLeadingPtRan","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);

  fh2PtNch = new TH2F("fh2PtNch","p_T of cluster vs. multiplicity; N_{ch};p_{T} (GeV/c);",nChMax,-0.5,nChMax-0.5,nBinPt,binLimitsPt);
  fh2PtNchRan = new TH2F("fh2PtNchRan","p_T of cluster vs. multiplicity ran; N_{ch};p_{T} (GeV/c);",nChMax,-0.5,nChMax-0.5,nBinPt,binLimitsPt);
  fh2PtNchN = new TH2F("fh2PtNchN","p_T of cluster vs. multiplicity N weighted; N_{ch};p_{T} (GeV/c);",nChMax,-0.5,nChMax-0.5,nBinPt,binLimitsPt);
  fh2PtNchNRan = new TH2F("fh2PtNchNRan","p_T of cluster vs. multiplicity N weighted ran; N_{ch};p_{T} (GeV/c);",nChMax,-0.5,nChMax-0.5,nBinPt,binLimitsPt);



  fh2JetPhiEta  = new TH2F("fh2JetPhiEta","eta vs phi all jets;#phi;#eta",
			   nBinPhi,0.,2.*TMath::Pi(),nBinEta,binLimitsEta);
  fh2LeadingJetPhiEta  = new TH2F("fh2LeadingJetPhiEta","eta vs phi leading jets;#phi;#eta",
				  nBinPhi,0.,2.*TMath::Pi(),nBinEta,binLimitsEta);

  fh2JetEtaPt  = new TH2F("fh2JetEtaPt","pt vs eta all jets;#eta;p_{T}",
			  nBinEta,binLimitsEta,nBinPt,binLimitsPt);
  fh2LeadingJetEtaPt  = new TH2F("fh2LeadingJetEtaPt","pT vs eta leading jets;#eta;p_{T}",
				 nBinEta,binLimitsEta,nBinPt,binLimitsPt);

  fh2TrackEtaPt  = new TH2F("fh2TrackEtaPt","pt vs eta all jets;#eta;p_{T}",
			    nBinEta,binLimitsEta,nBinPt,binLimitsPt);
  fh2LeadingTrackEtaPt  = new TH2F("fh2LeadingTrackEtaPt","pT vs eta leading jets;#eta;p_{T}",
				   nBinEta,binLimitsEta,nBinPt,binLimitsPt);



  fh2JetsLeadingPhiEta = new TH2F("fh2JetsLeadingPhiEta","delta eta vs delta phi to leading jet;#Delta#phi;#Delta#eta",
				  nBinPhi,binLimitsPhi,nBinEta,binLimitsEta);
  fh2JetsLeadingPhiPt = new TH2F("fh2JetsLeadingPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				 nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingPhiEta = new TH2F("fh2TracksLeadingPhiEta","delta eta vs delta phi to leading track;#Delta#phi;#Delta#eta",
				    nBinPhi,binLimitsPhi,nBinEta,binLimitsEta);
  fh2TracksLeadingPhiPt = new TH2F("fh2TracksLeadingPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				   nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingJetPhiPt = new TH2F("fh2TracksLeadingJetPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				      nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingJetPhiPtRan = new TH2F("fh2TracksLeadingJetPhiPtRan","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
					 nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);

  fh2JetsLeadingPhiPtW      = new TH2F("fh2JetsLeadingPhiPtW","leading p_T vs delta phi p_T weigted to leading jet;#Delta#phi;p_{T} (GeV/c)",
				       nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingPhiPtW    = new TH2F("fh2TracksLeadingPhiPtW","leading p_T vs delta phi to leading jet (p_T weighted);#Delta#phi;p_{T} (GeV/c)",
				       nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);

  fh2TracksLeadingJetPhiPtW = new TH2F("fh2TracksLeadingJetPhiPtW","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
				       nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);
  fh2TracksLeadingJetPhiPtWRan = new TH2F("fh2TracksLeadingJetPhiPtWRan","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
					  nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);

  if(fStoreRhoLeadingTrackCorr) {
    fh3CentvsRhoLeadingTrackPt = new TH3F("fh3CentvsRhoLeadingTrackPt","centrality vs background density full event; centrality; #rho", 50,0.,100.,500,0.,250.,100,0.,100.);
    fh3CentvsSigmaLeadingTrackPt = new TH3F("fh3CentvsSigmaLeadingTrackPt","centrality vs sigma full event; centrality; #sigma(#rho)", 50,0.,100.,50,0.,50.,100,0.,100.);
    fh3MultvsRhoLeadingTrackPt = new TH3F("fh3MultvsRhoLeadingTrackPt","multiplicity vs background density full event; multiplicity; #rho", 100,0.,5000.,500,0.,250.,100,0.,100.);
    fh3MultvsSigmaLeadingTrackPt = new TH3F("fh3MultvsSigmaLeadingTrackPt","multiplicity vs sigma full event; multiplicity; #sigma(#rho)", 100,0.,5000.,50,0.,50.,100,0.,100.);


    fh3CentvsRhoLeadingTrackPtQ1 = new TH3F("fh3CentvsRhoLeadingTrackPtQ1","centrality vs background density Q1; centrality; #rho; leading p_{t}^{track}", 50,0.,100.,500,0.,250.,100,0.,100.);
    fh3CentvsRhoLeadingTrackPtQ2 = new TH3F("fh3CentvsRhoLeadingTrackPtQ2","centrality vs background density Q2; centrality; #rho; leading p_{t}^{track}", 50,0.,100.,500,0.,250.,100,0.,100.);
    fh3CentvsRhoLeadingTrackPtQ3 = new TH3F("fh3CentvsRhoLeadingTrackPtQ3","centrality vs background density Q3; centrality; #rho; leading p_{t}^{track}", 50,0.,100.,500,0.,250.,100,0.,100.);
    fh3CentvsRhoLeadingTrackPtQ4 = new TH3F("fh3CentvsRhoLeadingTrackPtQ4","centrality vs background density Q4; centrality; #rho; leading p_{t}^{track}", 50,0.,100.,500,0.,250.,100,0.,100.);

    fh3CentvsSigmaLeadingTrackPtQ1 = new TH3F("fh3CentvsSigmaLeadingTrackPtQ1","centrality vs background density Q1; centrality; #sigma(#rho); leading p_{t}^{track}", 50,0.,100.,50,0.,50.,100,0.,100.);
    fh3CentvsSigmaLeadingTrackPtQ2 = new TH3F("fh3CentvsSigmaLeadingTrackPtQ2","centrality vs background density Q2; centrality; #sigma(#rho); leading p_{t}^{track}", 50,0.,100.,50,0.,50.,100,0.,100.);
    fh3CentvsSigmaLeadingTrackPtQ3 = new TH3F("fh3CentvsSigmaLeadingTrackPtQ3","centrality vs background density Q3; centrality; #sigma(#rho); leading p_{t}^{track}", 50,0.,100.,50,0.,50.,100,0.,100.);
    fh3CentvsSigmaLeadingTrackPtQ4 = new TH3F("fh3CentvsSigmaLeadingTrackPtQ4","centrality vs background density Q4; centrality; #sigma(#rho); leading p_{t}^{track}", 50,0.,100.,50,0.,50.,100,0.,100.);

    fh3MultvsRhoLeadingTrackPtQ1 = new TH3F("fh3MultvsRhoLeadingTrackPtQ1","multiplicity vs background density Q1; multiplicity; #rho; leading p_{t}^{track}", 100,0.,5000.,500,0.,250.,100,0.,100.);
    fh3MultvsRhoLeadingTrackPtQ2 = new TH3F("fh3MultvsRhoLeadingTrackPtQ2","multiplicity vs background density Q2; multiplicity; #rho; leading p_{t}^{track}", 100,0.,5000.,500,0.,250.,100,0.,100.);
    fh3MultvsRhoLeadingTrackPtQ3 = new TH3F("fh3MultvsRhoLeadingTrackPtQ3","multiplicity vs background density Q3; multiplicity; #rho; leading p_{t}^{track}", 100,0.,5000.,500,0.,250.,100,0.,100.);
    fh3MultvsRhoLeadingTrackPtQ4 = new TH3F("fh3MultvsRhoLeadingTrackPtQ4","multiplicity vs background density Q4; multiplicity; #rho; leading p_{t}^{track}", 100,0.,5000.,500,0.,250.,100,0.,100.);

    fh3MultvsSigmaLeadingTrackPtQ1 = new TH3F("fh3MultvsSigmaLeadingTrackPtQ1","multiplicity vs background density Q1; multiplicity; #sigma(#rho); leading p_{t}^{track}", 100,0.,5000.,50,0.,50.,100,0.,100.);
    fh3MultvsSigmaLeadingTrackPtQ2 = new TH3F("fh3MultvsSigmaLeadingTrackPtQ2","multiplicity vs background density Q2; multiplicity; #sigma(#rho); leading p_{t}^{track}", 100,0.,5000.,50,0.,50.,100,0.,100.);
    fh3MultvsSigmaLeadingTrackPtQ3 = new TH3F("fh3MultvsSigmaLeadingTrackPtQ3","multiplicity vs background density Q3; multiplicity; #sigma(#rho); leading p_{t}^{track}", 100,0.,5000.,50,0.,50.,100,0.,100.);
    fh3MultvsSigmaLeadingTrackPtQ4 = new TH3F("fh3MultvsSigmaLeadingTrackPtQ4","multiplicity vs background density Q4; multiplicity; #sigma(#rho); leading p_{t}^{track}", 100,0.,5000.,50,0.,50.,100,0.,100.);


    fh3CentvsDeltaRhoLeadingTrackPtQ1 = new TH3F("fh3CentvsDeltaRhoLeadingTrackPtQ1","centrality vs background density Q1; centrality; #rho_{quadrant}-#rho_{full event}; leading p_{t}^{track}", 50,0.,100.,200,-10.,10.,100,0.,100.);
    fh3CentvsDeltaRhoLeadingTrackPtQ2 = new TH3F("fh3CentvsDeltaRhoLeadingTrackPtQ2","centrality vs background density Q2; centrality; #rho_{quadrant}-#rho_{full event}; leading p_{t}^{track}", 50,0.,100.,200,-10.,10.,100,0.,100.);
    fh3CentvsDeltaRhoLeadingTrackPtQ3 = new TH3F("fh3CentvsDeltaRhoLeadingTrackPtQ3","centrality vs background density Q3; centrality; #rho_{quadrant}-#rho_{full event}; leading p_{t}^{track}", 50,0.,100.,200,-10.,10.,100,0.,100.);
    fh3CentvsDeltaRhoLeadingTrackPtQ4 = new TH3F("fh3CentvsDeltaRhoLeadingTrackPtQ4","centrality vs background density Q4; centrality; #rho_{quadrant}-#rho_{full event}; leading p_{t}^{track}", 50,0.,100.,200,-10.,10.,100,0.,100.);

    fHistList->Add(fh3CentvsRhoLeadingTrackPt);
    fHistList->Add(fh3CentvsSigmaLeadingTrackPt);
    fHistList->Add(fh3MultvsRhoLeadingTrackPt);
    fHistList->Add(fh3MultvsSigmaLeadingTrackPt);

    fHistList->Add(fh3CentvsRhoLeadingTrackPtQ1);
    fHistList->Add(fh3CentvsRhoLeadingTrackPtQ2);
    fHistList->Add(fh3CentvsRhoLeadingTrackPtQ3);
    fHistList->Add(fh3CentvsRhoLeadingTrackPtQ4);

    fHistList->Add(fh3CentvsSigmaLeadingTrackPtQ1);
    fHistList->Add(fh3CentvsSigmaLeadingTrackPtQ2);
    fHistList->Add(fh3CentvsSigmaLeadingTrackPtQ3);
    fHistList->Add(fh3CentvsSigmaLeadingTrackPtQ4);

    fHistList->Add(fh3MultvsRhoLeadingTrackPtQ1);
    fHistList->Add(fh3MultvsRhoLeadingTrackPtQ2);
    fHistList->Add(fh3MultvsRhoLeadingTrackPtQ3);
    fHistList->Add(fh3MultvsRhoLeadingTrackPtQ4);

    fHistList->Add(fh3MultvsSigmaLeadingTrackPtQ1);
    fHistList->Add(fh3MultvsSigmaLeadingTrackPtQ2);
    fHistList->Add(fh3MultvsSigmaLeadingTrackPtQ3);
    fHistList->Add(fh3MultvsSigmaLeadingTrackPtQ4);

    fHistList->Add(fh3CentvsDeltaRhoLeadingTrackPtQ1);
    fHistList->Add(fh3CentvsDeltaRhoLeadingTrackPtQ2);
    fHistList->Add(fh3CentvsDeltaRhoLeadingTrackPtQ3);
    fHistList->Add(fh3CentvsDeltaRhoLeadingTrackPtQ4);

  }
  
  //Detector level effects histos
  fh2PtGenPtSmeared = new TH2F("fh2PtGenPtSmeared","fh2PtGenPtSmeared",nBinPt,binLimitsPt,nBinPt,binLimitsPt);

  fp1Efficiency = new TProfile("fp1Efficiency","fp1Efficiency",nBinPt,binLimitsPt);
  fp1PtResolution = new TProfile("fp1PtResolution","fp1PtResolution",nBinPt,binLimitsPt);

  fHistList->Add(fh2PtGenPtSmeared);
  fHistList->Add(fp1Efficiency);
  fHistList->Add(fp1PtResolution);

  if(fNRandomCones>0&&fUseBackgroundCalc){
    for(int i = 0;i<3;i++){
      fh1BiARandomCones[i] = new TH1F(Form("fh1BiARandomCones%d",i),";B_{i}^{A} (GeV/c)",200,-100,100);
      fh1BiARandomConesRan[i] =  new TH1F(Form("fh1BiARandomConesRan%d",i),";B_{i}^{A} (GeV/c)",200,-100,100);
    }
  }

  for(int i = 0;i < kMaxCent;i++){
    fh2JetsLeadingPhiPtC[i] = (TH2F*)fh2JetsLeadingPhiPt->Clone(Form("%s_C%02d",fh2JetsLeadingPhiPt->GetName(),i+1));
    fh2JetsLeadingPhiPtWC[i]= (TH2F*)fh2JetsLeadingPhiPtW->Clone(Form("%s_C%02d",fh2JetsLeadingPhiPtW->GetName(),i+1));
    fh2TracksLeadingJetPhiPtC[i] = (TH2F*)fh2TracksLeadingJetPhiPt->Clone(Form("%s_C%02d",fh2TracksLeadingJetPhiPt->GetName(),i+1));
    fh2TracksLeadingJetPhiPtWC[i] = (TH2F*)fh2TracksLeadingJetPhiPtW->Clone(Form("%s_C%02d",fh2TracksLeadingJetPhiPtW->GetName(),i+1));
  }

  const Int_t saveLevel = 3; // large save level more histos
  if(saveLevel>0){
    fHistList->Add(fh1Xsec);
    fHistList->Add(fh1Trials);

    fHistList->Add(fh1NJetsRec);
    fHistList->Add(fh1NConstRec);
    fHistList->Add(fh1NConstLeadingRec);
    fHistList->Add(fh1PtJetsRecIn);
    fHistList->Add(fh1PtJetsLeadingRecIn);
    fHistList->Add(fh1PtTracksRecIn);
    fHistList->Add(fh1PtTracksLeadingRecIn);
    fHistList->Add(fh1PtJetConstRec);
    fHistList->Add(fh1PtJetConstLeadingRec);
    fHistList->Add(fh1NJetsRecRan);
    fHistList->Add(fh1NConstRecRan);
    fHistList->Add(fh1PtJetsLeadingRecInRan);
    fHistList->Add(fh1NConstLeadingRecRan);
    fHistList->Add(fh1PtJetsRecInRan);
    fHistList->Add(fh1Nch);
    fHistList->Add(fh1Centrality);
    fHistList->Add(fh1CentralitySelect);
    fHistList->Add(fh1CentralityPhySel);
    fHistList->Add(fh1Z);
    fHistList->Add(fh1ZSelect);
    fHistList->Add(fh1ZPhySel);
    if(fNRandomCones>0&&fUseBackgroundCalc){
      for(int i = 0;i<3;i++){
	fHistList->Add(fh1BiARandomCones[i]);
	fHistList->Add(fh1BiARandomConesRan[i]);
      }
    }
    for(int i = 0;i < kMaxCent;i++){
      fHistList->Add(fh2JetsLeadingPhiPtC[i]);
      fHistList->Add(fh2JetsLeadingPhiPtWC[i]);
      fHistList->Add(fh2TracksLeadingJetPhiPtC[i]);
      fHistList->Add(fh2TracksLeadingJetPhiPtWC[i]);
    }

    fHistList->Add(fh2NRecJetsPt);
    fHistList->Add(fh2NRecTracksPt);
    fHistList->Add(fh2NConstPt);
    fHistList->Add(fh2NConstLeadingPt);
    fHistList->Add(fh2PtNch);
    fHistList->Add(fh2PtNchRan);
    fHistList->Add(fh2PtNchN);
    fHistList->Add(fh2PtNchNRan);
    fHistList->Add(fh2JetPhiEta);
    fHistList->Add(fh2LeadingJetPhiEta);
    fHistList->Add(fh2JetEtaPt);
    fHistList->Add(fh2LeadingJetEtaPt);
    fHistList->Add(fh2TrackEtaPt);
    fHistList->Add(fh2LeadingTrackEtaPt);
    fHistList->Add(fh2JetsLeadingPhiEta );
    fHistList->Add(fh2JetsLeadingPhiPt);
    fHistList->Add(fh2TracksLeadingPhiEta);
    fHistList->Add(fh2TracksLeadingPhiPt);
    fHistList->Add(fh2TracksLeadingJetPhiPt);
    fHistList->Add(fh2JetsLeadingPhiPtW);
    fHistList->Add(fh2TracksLeadingPhiPtW);
    fHistList->Add(fh2TracksLeadingJetPhiPtW);
    fHistList->Add(fh2NRecJetsPtRan);
    fHistList->Add(fh2NConstPtRan);
    fHistList->Add(fh2NConstLeadingPtRan);
    fHistList->Add(fh2TracksLeadingJetPhiPtRan);
    fHistList->Add(fh2TracksLeadingJetPhiPtWRan);
  }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fHistList->At(i));
    if(hn)hn->Sumw2();
  }
  TH1::AddDirectory(oldStatus);
}

void AliAnalysisTaskJetCluster::LocalInit()
{
  //
  // Initialization
  //

  if (fDebug > 1) printf("AnalysisTaskJetCluster::Init() \n");

  if(fUseTrPtResolutionFromOADB) LoadTrPtResolutionRootFileFromOADB();
  if(fUseTrEfficiencyFromOADB)   LoadTrEfficiencyRootFileFromOADB();


  FitMomentumResolution();

}

void AliAnalysisTaskJetCluster::UserExec(Option_t */*option*/)
{

  // handle and reset the output jet branch 

  if(fTCAJetsOut)fTCAJetsOut->Delete();
  if(fTCAJetsOutRan)fTCAJetsOutRan->Delete();
  if(fTCARandomConesOut)fTCARandomConesOut->Delete();
  if(fTCARandomConesOutRan)fTCARandomConesOutRan->Delete();
  if(fAODJetBackgroundOut)fAODJetBackgroundOut->Reset();

  AliAODJetEventBackground* externalBackground = 0;
  if(!externalBackground&&fBackgroundBranch.Length()){
    externalBackground =  (AliAODJetEventBackground*)(AODEvent()->FindListObject(fBackgroundBranch.Data()));
    if((!externalBackground)&&fAODExtension)externalBackground = (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch.Data()));
    if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
  }
  //
  // Execute analysis for current event
  //
  AliESDEvent *fESD = 0;
  if(fUseAODTrackInput){    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODTrackInput);
      return;
    }
    // fetch the header
  }
  else{
    //  assume that the AOD is in the general output...
    fAOD  = AODEvent();
    if(!fAOD){
      Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
      return;
    }
    if(fDebug>0){
      fESD = dynamic_cast<AliESDEvent*> (InputEvent());
    }
  }

  //Check if information is provided detector level effects
  if(!fMomResH1 || !fMomResH2 || !fMomResH3) fUseTrPtResolutionSmearing = kFALSE;
  if(!fhEffH1 || !fhEffH2 || !fhEffH3 )      fUseDiceEfficiency = kFALSE;
  if(  fEfficiencyFixed < 1. )               fUseDiceEfficiency = kTRUE;  

  Bool_t selectEvent =  false;
  Bool_t physicsSelection = true;// handled by the framework(fInputHandler->IsEventSelected()&AliVEvent::kMB)==AliVEvent::kMB;

  Float_t cent = 0;
  Float_t zVtx  = 0;
  Int_t cenClass = -1;
  if(fAOD){
    const AliAODVertex *vtxAOD = fAOD->GetPrimaryVertex();
    TString vtxTitle(vtxAOD->GetTitle());
    zVtx = vtxAOD->GetZ();

    cent = fAOD->GetHeader()->GetCentrality();
    if(cent<10)cenClass = 0;
    else if(cent<30)cenClass = 1;
    else if(cent<50)cenClass = 2;
    else if(cent<80)cenClass = 3;
    if(physicsSelection){
      fh1CentralityPhySel->Fill(cent);
      fh1ZPhySel->Fill(zVtx);
    }

    if(fEventSelection){
      if(vtxAOD->GetNContributors()>2&&!vtxTitle.Contains("TPCVertex")){
	Float_t yvtx = vtxAOD->GetY();
	Float_t xvtx = vtxAOD->GetX();
	Float_t r2   = yvtx*yvtx+xvtx*xvtx;  
	if(TMath::Abs(zVtx)<fVtxZCut&&r2<fVtxR2Cut){ // apply vertex cut later on
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
  
  if(fRequireVZEROAC&&fRequireTZEROvtx)selectEvent = selectEvent&&V0&&T0;
  else if(fRequireTZEROvtx)selectEvent = selectEvent&&T0;
  else if(fRequireVZEROAC)selectEvent = selectEvent&&V0;


  if(!selectEvent){
    PostData(1, fHistList);
    return;
  }
  fh1Centrality->Fill(cent);  
  fh1Z->Fill(zVtx);
  fh1Trials->Fill("#sum{ntrials}",1);
  

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);

  // ==== General variables needed



  // we simply fetch the tracks/mc particles as a list of AliVParticles

  TList recParticles;
  TList genParticles;

  Int_t nT = GetListOfTracks(&recParticles,fTrackTypeRec);
  Float_t nCh = recParticles.GetEntries(); 
  fh1Nch->Fill(nCh);
  if(fDebug>2)Printf("%s:%d Selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,recParticles.GetEntries());
  nT = GetListOfTracks(&genParticles,fTrackTypeGen);
  if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nT,genParticles.GetEntries());

  // find the jets....

  vector<fastjet::PseudoJet> inputParticlesRec;
  vector<fastjet::PseudoJet> inputParticlesRecRan;
  
  // Generate the random cones
  
  AliAODJet vTmpRan(1,0,0,1);
  for(int i = 0; i < recParticles.GetEntries(); i++){
    AliVParticle *vp = (AliVParticle*)recParticles.At(i);

    // Carefull energy is not well determined in real data, should not matter for p_T scheme?
    // we take total momentum here

    if((!fUseTrPtResolutionSmearing) && (!fUseDiceEfficiency)) {
      //Add particles to fastjet in case we are not running toy model
      fastjet::PseudoJet jInp(vp->Px(),vp->Py(),vp->Pz(),vp->P());
      jInp.set_user_index(i);
      inputParticlesRec.push_back(jInp);
    }
    else if(fUseDiceEfficiency) {

      // Dice to decide if particle is kept or not - toy  model for efficiency
      //
      Double_t sumEff = 0.;
      Double_t pT = 0.;
      Double_t eff[3] = {0.};
      Int_t cat[3] = {0};

      Double_t rnd = fRandom->Uniform(1.);
      if(  fEfficiencyFixed<1. ) {
	sumEff = fEfficiencyFixed;
      } else {

	pT = vp->Pt();
	Double_t pTtmp = pT;
	if(pT>10.) pTtmp = 10.;
	Double_t eff1 = fhEffH1->GetBinContent(fhEffH1->FindBin(pTtmp));
	Double_t eff2 = fhEffH2->GetBinContent(fhEffH2->FindBin(pTtmp));
	Double_t eff3 = fhEffH3->GetBinContent(fhEffH3->FindBin(pTtmp));

	//Sort efficiencies from large to small
	if(eff1>eff2 && eff1>eff3) { 
	  eff[0] = eff1; 
	  cat[0] = 1;
	  if(eff2>eff3) {
	    eff[1] = eff2;
	    eff[2] = eff3; 
	    cat[1] = 2;
	    cat[2] = 3;
	  } else {
	    eff[1] = eff3;
	    eff[2] = eff2; 
	    cat[1] = 3;
	    cat[2] = 2;
	  }
	}
	else if(eff2>eff1 && eff2>eff3) {
	  eff[0] = eff2;
	  cat[0] = 2;
	  if(eff1>eff3) {
	    eff[1] = eff1;
	    eff[2] = eff3; 
	    cat[1] = 1;
	    cat[2] = 3;
	  } else {
	    eff[1] = eff3;
	    eff[2] = eff1; 
	    cat[1] = 3;
	    cat[2] = 1;
	  }
	}
	else if(eff3>eff1 && eff3>eff2) {
	  eff[0] = eff3;
	  cat[0] = 3;
	  if(eff1>eff2) {
	    eff[1] = eff1;
	    eff[2] = eff2; 
	    cat[1] = 1;
	    cat[2] = 2;
	  } else {
	    eff[1] = eff2;
	    eff[2] = eff1; 
	    cat[1] = 2;
	    cat[2] = 1;
	  }
	}
	
	sumEff = eff[0]+eff[1]+eff[2];
      }
      fp1Efficiency->Fill(vp->Pt(),sumEff);
      if(rnd>sumEff && pT > fDiceEfficiencyMinPt) continue;

      if(fUseTrPtResolutionSmearing) {
	//Smear momentum of generated particle
	Double_t smear = 1.;
	//Select hybrid track category
	if(rnd<=eff[2]) 
	  smear = GetMomentumSmearing(cat[2],pT);
	else if(rnd<=(eff[2]+eff[1])) 
	  smear = GetMomentumSmearing(cat[1],pT);
	else 
	  smear = GetMomentumSmearing(cat[0],pT);

	fp1PtResolution->Fill(vp->Pt(),smear);

	Double_t sigma = vp->Pt()*smear;
	Double_t pTrec = fRandom->Gaus(vp->Pt(),sigma);
	
	Double_t phi   = vp->Phi();
	Double_t theta = 2.*TMath::ATan(TMath::Exp(-1.*(vp->Eta())));
	Double_t pX    = pTrec * TMath::Cos(phi);
	Double_t pY    = pTrec * TMath::Sin(phi);
	Double_t pZ    = pTrec/TMath::Tan(theta);
	Double_t p=TMath::Sqrt(pTrec*pTrec+pZ*pZ);

	fh2PtGenPtSmeared->Fill(vp->Pt(),pTrec);

	fastjet::PseudoJet jInp(pX,pY,pZ,p);
	jInp.set_user_index(i);
	inputParticlesRec.push_back(jInp);

      }
      else {
	fastjet::PseudoJet jInp(vp->Px(),vp->Py(),vp->Pz(),vp->P());
	jInp.set_user_index(i);
	inputParticlesRec.push_back(jInp);

      }

    }

    // the randomized input changes eta and phi, but keeps the p_T
    if(i>=fNSkipLeadingRan){// eventually skip the leading particles
      Double_t pT = vp->Pt();
      Double_t eta = 2.*fTrackEtaWindow * fRandom->Rndm() - fTrackEtaWindow;
      Double_t phi = 2.* TMath::Pi() * fRandom->Rndm();
      
      Double_t theta = 2.*TMath::ATan(TMath::Exp(-eta));  
      Double_t pZ = pT/TMath::Tan(theta);

      Double_t pX = pT * TMath::Cos(phi);
      Double_t pY = pT * TMath::Sin(phi);
      Double_t p  = TMath::Sqrt(pT*pT+pZ*pZ); 
      fastjet::PseudoJet jInpRan(pX,pY,pZ,p);

      jInpRan.set_user_index(i);
      inputParticlesRecRan.push_back(jInpRan);
      vTmpRan.SetPxPyPzE(pX,pY,pZ,p);
    }

    // fill the tref array, only needed when we write out jets
    if(fTCAJetsOut){
      if(i == 0){
	fRef->Delete(); // make sure to delete before placement new...
	if((!fUseTrPtResolutionSmearing) && (!fUseDiceEfficiency)) {
	  new(fRef) TRefArray(TProcessID::GetProcessWithUID(vp)); //TRefArray does not work with toy model ...
	} 
      }
      if((!fUseTrPtResolutionSmearing) && (!fUseDiceEfficiency)) fRef->Add(vp);  //TRefArray does not work with toy model ...
    }
  }// recparticles

  if(inputParticlesRec.size()==0){
    if(fDebug)Printf("%s:%d No input particles found, skipping event",(char*)__FILE__,__LINE__);
    PostData(1, fHistList);
    return;
  }
  
  // run fast jet
  // employ setters for these...

 
  // now create the object that holds info about ghosts                        
  /*
    if(!fUseBackgroundCalc&& fNonStdBranch.Length()==0){
    // reduce CPU time...
    fGhostArea = 0.5; 
    fActiveAreaRepeats = 0; 
    }
  */

  fastjet::GhostedAreaSpec ghostSpec(fGhostEtamax, fActiveAreaRepeats, fGhostArea);
  fastjet::AreaType areaType =   fastjet::active_area;
  fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
  fastjet::JetDefinition jetDef(fAlgorithm, fRparam, fRecombScheme, fStrategy);
  fastjet::ClusterSequenceArea clustSeq(inputParticlesRec, jetDef,areaDef);
  
  //range where to compute background
  Double_t phiMin = 0, phiMax = 0, rapMin = 0, rapMax = 0;
  phiMin = 0;
  phiMax = 2*TMath::Pi();
  rapMax = fGhostEtamax - fRparam;
  rapMin = - fGhostEtamax + fRparam;
  fastjet::RangeDefinition range(rapMin,rapMax, phiMin, phiMax);
 

  const vector <fastjet::PseudoJet> &inclusiveJets = clustSeq.inclusive_jets();
  const vector <fastjet::PseudoJet> &sortedJets = sorted_by_pt(inclusiveJets);

 
  fh1NJetsRec->Fill(sortedJets.size());

  // loop over all jets an fill information, first on is the leading jet

  Int_t nRecOver = inclusiveJets.size();
  Int_t nRec     = inclusiveJets.size();
  if(inclusiveJets.size()>0){
    AliAODJet leadingJet (sortedJets[0].px(), sortedJets[0].py(), sortedJets[0].pz(), sortedJets[0].E());
    Double_t area = clustSeq.area(sortedJets[0]);
    leadingJet.SetEffArea(area,0);
    Float_t pt = leadingJet.Pt();
    Int_t nAodOutJets = 0;
    Int_t nAodOutTracks = 0;
    AliAODJet *aodOutJet = 0;

    Int_t iCount = 0;
    for(int i = 1;i <= fh2NRecJetsPt->GetNbinsX();i++){
      Float_t ptCut = fh2NRecJetsPt->GetXaxis()->GetBinCenter(i);
      while(pt<ptCut&&iCount<nRec){
	nRecOver--;
	iCount++;
	if(iCount<nRec){
	  pt = sortedJets[iCount].perp();
	}
      }
      if(nRecOver<=0)break;
      fh2NRecJetsPt->Fill(ptCut,nRecOver);
    }
    Float_t phi = leadingJet.Phi();
    if(phi<0)phi+=TMath::Pi()*2.;    
    Float_t eta = leadingJet.Eta();
    Float_t pTback = 0;
    if(externalBackground){
      // carefull has to be filled in a task before
      // todo, ReArrange to the botom 
      pTback = externalBackground->GetBackground(2)*leadingJet.EffectiveAreaCharged();
    }
    pt = leadingJet.Pt() - pTback;
    // correlation of leading jet with tracks
    TIterator *recIter = recParticles.MakeIterator();
    recIter->Reset();
    AliVParticle *tmpRecTrack = 0;
    while((tmpRecTrack = (AliVParticle*)(recIter->Next()))){
      Float_t tmpPt = tmpRecTrack->Pt();
      // correlation
      Float_t tmpPhi =  tmpRecTrack->Phi();     
      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    
      Float_t dPhi = phi - tmpPhi;
      if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
      if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
      fh2TracksLeadingJetPhiPt->Fill(dPhi,pt);
      fh2TracksLeadingJetPhiPtW->Fill(dPhi,pt,tmpPt);
      if(cenClass>=0){
	fh2TracksLeadingJetPhiPtC[cenClass]->Fill(dPhi,pt);
	fh2TracksLeadingJetPhiPtWC[cenClass]->Fill(dPhi,pt,tmpPt);
      }

    }  
    
   
    TLorentzVector vecareab;
    for(int j = 0; j < nRec;j++){
      AliAODJet tmpRec (sortedJets[j].px(), sortedJets[j].py(), sortedJets[j].pz(), sortedJets[j].E());
      aodOutJet = 0;
      nAodOutTracks = 0;
      Float_t tmpPt = tmpRec.Pt();  
      
      if(tmpPt>fJetOutputMinPt&&fTCAJetsOut){// cut on the non-background subtracted...
	aodOutJet =  new ((*fTCAJetsOut)[nAodOutJets++]) AliAODJet(tmpRec);
	aodOutJet->GetRefTracks()->Clear();
	Double_t area1 = clustSeq.area(sortedJets[j]);
	aodOutJet->SetEffArea(area1,0);
        fastjet::PseudoJet vecarea=clustSeq.area_4vector(sortedJets[j]);  
        vecareab.SetPxPyPzE(vecarea.px(),vecarea.py(),vecarea.pz(),vecarea.e());     
	aodOutJet->SetVectorAreaCharged(&vecareab);
      }


      Float_t tmpPtBack = 0;
      if(externalBackground){
	// carefull has to be filled in a task before
	// todo, ReArrange to the botom
	tmpPtBack = externalBackground->GetBackground(2)*tmpRec.EffectiveAreaCharged();
      }
      tmpPt = tmpPt - tmpPtBack;
      if(tmpPt<0)tmpPt = 0; // avoid negative weights...
      
      fh1PtJetsRecIn->Fill(tmpPt);
      // Fill Spectra with constituentsemacs
      const vector<fastjet::PseudoJet> &constituents = clustSeq.constituents(sortedJets[j]);

      fh1NConstRec->Fill(constituents.size());
      fh2PtNch->Fill(nCh,tmpPt);
      fh2PtNchN->Fill(nCh,tmpPt,constituents.size());
      fh2NConstPt->Fill(tmpPt,constituents.size());
      // loop over constiutents and fill spectrum
   
      AliVParticle *partLead = 0;
      Float_t ptLead = -1;

      for(unsigned int ic = 0; ic < constituents.size();ic++){
	AliVParticle *part = (AliVParticle*)recParticles.At(constituents[ic].user_index());
	if(!part) continue;
	fh1PtJetConstRec->Fill(part->Pt());
	if(aodOutJet){
	  if((!fUseTrPtResolutionSmearing) && (!fUseDiceEfficiency)) aodOutJet->AddTrack(fRef->At(constituents[ic].user_index()));
	  if(part->Pt()>fMaxTrackPtInJet){
	    aodOutJet->SetTrigger(AliAODJet::kHighTrackPtTriggered);
	  }
	}
	if(part->Pt()>ptLead){
	  ptLead = part->Pt();
	  partLead = part;
	}
	if(j==0)fh1PtJetConstLeadingRec->Fill(part->Pt());
      }

      AliAODTrack *aodT = 0;
      if(partLead){
	if(aodOutJet){
	  //set pT of leading constituent of jet
	  aodOutJet->SetPtLeading(partLead->Pt());
	  aodT = dynamic_cast<AliAODTrack*>(partLead);
	  if(aodT){
	    if(aodT->TestFilterBit(fFilterMaskBestPt)){
	      aodOutJet->SetTrigger(AliAODJet::kHighTrackPtBest);
	    }
	  }
	}
      }
    
      // correlation
      Float_t tmpPhi =  tmpRec.Phi();
      Float_t tmpEta =  tmpRec.Eta();
      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;        
      if(j==0){
	fh1PtJetsLeadingRecIn->Fill(tmpPt);
	fh2LeadingJetPhiEta->Fill(tmpPhi,tmpEta);
	fh2LeadingJetEtaPt->Fill(tmpEta,tmpPt);
	fh1NConstLeadingRec->Fill(constituents.size());
	fh2NConstLeadingPt->Fill(tmpPt,constituents.size());
	continue;
      }
      fh2JetPhiEta->Fill(tmpRec.Phi(),tmpEta);
      fh2JetEtaPt->Fill(tmpEta,tmpPt);
      Float_t dPhi = phi - tmpPhi;
      if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
      if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
      Float_t dEta = eta - tmpRec.Eta();
      fh2JetsLeadingPhiEta->Fill(dPhi,dEta);
      fh2JetsLeadingPhiPt->Fill(dPhi,pt);
      if(cenClass>=0){
	fh2JetsLeadingPhiPtC[cenClass]->Fill(dPhi,pt);
	fh2JetsLeadingPhiPtWC[cenClass]->Fill(dPhi,pt,tmpPt);
      }
      fh2JetsLeadingPhiPtW->Fill(dPhi,pt,tmpPt);
    }// loop over reconstructed jets
    delete recIter;



    // Add the random cones...
    if(fNRandomCones>0&&fTCARandomConesOut){       
      // create a random jet within the acceptance
      Double_t etaMax = fTrackEtaWindow - fRparam;
      Int_t nCone = 0;
      Int_t nConeRan = 0;
      Double_t pTC = 1; // small number
      for(int ir = 0;ir < fNRandomCones;ir++){
	Double_t etaC = etaMax*2.*(fRandom->Rndm()-0.5); // +- etamax
	Double_t phiC = fRandom->Rndm()*2.*TMath::Pi(); // 0 - 2pi
	// massless jet
	Double_t thetaC = 2.*TMath::ATan(TMath::Exp(-etaC));  
	Double_t pZC = pTC/TMath::Tan(thetaC);
	Double_t pXC = pTC * TMath::Cos(phiC);
	Double_t pYC = pTC * TMath::Sin(phiC);
	Double_t pC  = TMath::Sqrt(pTC*pTC+pZC*pZC); 
	AliAODJet tmpRecC (pXC,pYC,pZC, pC); 
	bool skip = false;
	for(int jj = 0; jj < TMath::Min(nRec,fNSkipLeadingCone);jj++){// test for overlap with leading jets
	  AliAODJet jet (sortedJets[jj].px(), sortedJets[jj].py(), sortedJets[jj].pz(), sortedJets[jj].E());
	  if(jet.DeltaR(& tmpRecC)<2.*fRparam+0.2){
	    skip = true;
	    break;
	  }
	}
	// test for overlap with previous cones to avoid double counting
	for(int iic = 0;iic<ir;iic++){
	  AliAODJet *iicone = (AliAODJet*)fTCARandomConesOut->At(iic);
	  if(iicone){
	    if(iicone->DeltaR(&tmpRecC)<2.*fRparam){
	      skip = true;
	      break;
	    }
	  }
	}
	if(skip)continue;
	tmpRecC.SetBgEnergy(0,0); // this is use as temporary storage of the summed p_T below
	tmpRecC.SetPtLeading(-1.);
	if(fTCARandomConesOut)new ((*fTCARandomConesOut)[nCone++]) AliAODJet(tmpRecC);
	if(fTCARandomConesOutRan)new ((*fTCARandomConesOutRan)[nConeRan++]) AliAODJet(tmpRecC);
      }// loop over random cones creation

  
      // loop over the reconstructed particles and add up the pT in the random cones
      // maybe better to loop over randomized particles not in the real jets...
      // but this by definition brings dow average energy in the whole  event
      AliAODJet vTmpRanR(1,0,0,1);
      for(int i = 0; i < recParticles.GetEntries(); i++){
	AliVParticle *vp = (AliVParticle*)recParticles.At(i);
	if(fTCARandomConesOut){
	  for(int ir = 0;ir < fNRandomCones;ir++){
	    AliAODJet *jC = (AliAODJet*)fTCARandomConesOut->At(ir);  
	    if(jC&&jC->DeltaR(vp)<fRparam){
	      if(vp->Pt()>fMaxTrackPtInJet)jC->SetTrigger(AliAODJet::kHighTrackPtTriggered);
	      jC->SetBgEnergy(jC->ChargedBgEnergy()+vp->Pt(),0);
	      if(vp->Pt() > jC->GetPtLeading()) jC->SetPtLeading(vp->Pt());
	    }
	  }  
	}// add up energy in cone

	// the randomized input changes eta and phi, but keeps the p_T
	if(i>=fNSkipLeadingRan){// eventually skip the leading particles
	  Double_t pTR = vp->Pt();
	  Double_t etaR = 2.*fTrackEtaWindow* fRandom->Rndm() - fTrackEtaWindow;
	  Double_t phiR = 2.* TMath::Pi() * fRandom->Rndm();
	 
	  Double_t thetaR = 2.*TMath::ATan(TMath::Exp(-etaR));  
	  Double_t pZR = pTR/TMath::Tan(thetaR);
	 
	  Double_t pXR = pTR * TMath::Cos(phiR);
	  Double_t pYR = pTR * TMath::Sin(phiR);
	  Double_t pR  = TMath::Sqrt(pTR*pTR+pZR*pZR); 
	  vTmpRanR.SetPxPyPzE(pXR,pYR,pZR,pR);
	  if(fTCARandomConesOutRan){
	    for(int ir = 0;ir < fTCARandomConesOutRan->GetEntriesFast();ir++){
	      AliAODJet *jC = (AliAODJet*)fTCARandomConesOutRan->At(ir);  
	      if(jC&&jC->DeltaR(&vTmpRanR)<fRparam){
		if(vTmpRanR.Pt()>fMaxTrackPtInJet)jC->SetTrigger(AliAODJet::kHighTrackPtTriggered);
		jC->SetBgEnergy(jC->ChargedBgEnergy()+vTmpRanR.Pt(),0);
		if(vTmpRanR.Pt() > jC->GetPtLeading()) jC->SetPtLeading(vTmpRanR.Pt());
	      }
	    }  
	  }
	}
      }// loop over recparticles
    
      Float_t jetArea = fRparam*fRparam*TMath::Pi();
      if(fTCARandomConesOut){
	for(int ir = 0;ir < fTCARandomConesOut->GetEntriesFast();ir++){
	  // rescale the momentum vectors for the random cones
	 
	  AliAODJet *rC = (AliAODJet*)fTCARandomConesOut->At(ir);
	  if(rC){
	    Double_t etaC = rC->Eta();
	    Double_t phiC = rC->Phi();
	    // massless jet, unit vector
	    pTC = rC->ChargedBgEnergy();
	    if(pTC<=0)pTC = 0.001; // for almost empty events
	    Double_t thetaC = 2.*TMath::ATan(TMath::Exp(-etaC));  
	    Double_t pZC = pTC/TMath::Tan(thetaC);
	    Double_t pXC = pTC * TMath::Cos(phiC);
	    Double_t pYC = pTC * TMath::Sin(phiC);
	    Double_t pC  = TMath::Sqrt(pTC*pTC+pZC*pZC); 
	    rC->SetPxPyPzE(pXC,pYC,pZC, pC); 
	    rC->SetBgEnergy(0,0);
	    rC->SetEffArea(jetArea,0);
	  }
	}
      }
      if(fTCARandomConesOutRan){
	for(int ir = 0;ir < fTCARandomConesOutRan->GetEntriesFast();ir++){
	  AliAODJet* rC = (AliAODJet*)fTCARandomConesOutRan->At(ir);
	  // same wit random
	  if(rC){
	    Double_t etaC = rC->Eta();
	    Double_t phiC = rC->Phi();
	    // massless jet, unit vector
	    pTC = rC->ChargedBgEnergy();
	    if(pTC<=0)pTC = 0.001;// for almost empty events
	    Double_t thetaC = 2.*TMath::ATan(TMath::Exp(-etaC));  
	    Double_t pZC = pTC/TMath::Tan(thetaC);
	    Double_t pXC = pTC * TMath::Cos(phiC);
	    Double_t pYC = pTC * TMath::Sin(phiC);
	    Double_t pC  = TMath::Sqrt(pTC*pTC+pZC*pZC); 
	    rC->SetPxPyPzE(pXC,pYC,pZC, pC); 
	    rC->SetBgEnergy(0,0);
	    rC->SetEffArea(jetArea,0);
	  }
	}
      }
    }// if(fNRandomCones
  
    //background estimates:all bckg jets(0) & wo the 2 hardest(1)

    if(fAODJetBackgroundOut){
      vector<fastjet::PseudoJet> jets2=sortedJets;
      if(jets2.size()>2) jets2.erase(jets2.begin(),jets2.begin()+2); 

      Double_t bkg1=0;
      Double_t sigma1=0.;
      Double_t meanarea1=0.;
      Double_t bkg2=0;
      Double_t sigma2=0.;
      Double_t meanarea2=0.;

      Double_t bkg4=0;
      Double_t sigma4=0.;
      Double_t meanarea4=0.;

      clustSeq.get_median_rho_and_sigma(jets2, range, true, bkg1, sigma1, meanarea1, true);
      fAODJetBackgroundOut->SetBackground(0,bkg1,sigma1,meanarea1);

      //     fh1BiARandomCones[0]->Fill(omCone-(bkg1*areaRandomCone));    
      //  fh1BiARandomConesRan[0]->Fill(ptRandomConeRan-(bkg1*areaRandomCone));      
      clustSeq.get_median_rho_and_sigma(sortedJets, range, true, bkg4, sigma4, meanarea4, true);
      fAODJetBackgroundOut->SetBackground(3,bkg4,sigma4,meanarea4);

      //////////////////////
     
      clustSeq.get_median_rho_and_sigma(jets2, range, false, bkg2, sigma2, meanarea2, true);
      fAODJetBackgroundOut->SetBackground(1,bkg2,sigma2,meanarea2);
      //  fh1BiARandomCones[1]->Fill(ptRandomCone-(bkg2*areaRandomCone));    
      //   fh1BiARandomConesRan[1]->Fill(ptRandomConeRan-(bkg2*areaRandomCone));    

    }
  }

  if(fStoreRhoLeadingTrackCorr) {
    vector<fastjet::PseudoJet> jets3=sortedJets;
    if(jets3.size()>2) jets3.erase(jets3.begin(),jets3.begin()+2); 

    Double_t bkg2=0;
    Double_t sigma2=0.;
    Double_t meanarea2=0.;

    clustSeq.get_median_rho_and_sigma(jets3, range, false, bkg2, sigma2, meanarea2, true);
    fAODJetBackgroundOut->SetBackground(1,bkg2,sigma2,meanarea2);

    //Get leading track in event
    AliVParticle *leading = (AliVParticle*)recParticles.At(0);

    fh3CentvsRhoLeadingTrackPt->Fill(cent,bkg2,leading->Pt());
    fh3CentvsSigmaLeadingTrackPt->Fill(cent,sigma2,leading->Pt());
    fh3MultvsRhoLeadingTrackPt->Fill(nCh,bkg2,leading->Pt());
    fh3MultvsSigmaLeadingTrackPt->Fill(nCh,sigma2,leading->Pt());
     

    //Calculate rho with 4-vector area method for different quadrants with respect to the leading track in the event
    //exclude 2 leading jets from event
    //Quadrant 1: |phi_leadingTrack - phi_bkgCluster| < pi/2./2. - R  (Near side to leading track)
    //Quadrant 2: pi/2 - (pi/2./2. - R) < |phi_leadingTrack - phi_bkgCluster| < pi/2 + (pi/2./2. - R) (Orthogonal to leading track)
    //Quadrant 3: pi - (pi/2/2 - R) < |phi_leadingTrack - phi_bkgCluster| < pi + (pi/2./2. - R)  (Away side to leading track)
    //Quadrant 4: 3/2*pi - (pi/2./2. - R) < |phi_leadingTrack - phi_bkgCluster| < 3/2*pi + (pi/2./2. - R)   (Orthogonal to leading track)

    //Assuming R=0.2 for background clusters

    Double_t bkg2Q[4]      = {0.};
    Double_t sigma2Q[4]    = {0.};
    Double_t meanarea2Q[4] = {0.};

    //Get phi of leading track in event
    Float_t phiLeadingTrack = leading->Phi();
    Float_t phiStep = (TMath::Pi()/2./2. - 0.2);

    //Quadrant1 - near side
    phiMin = phiLeadingTrack - phiStep;
    phiMax = phiLeadingTrack + phiStep;
    fastjet::RangeDefinition rangeQ1(rapMin,rapMax, phiMin, phiMax);
    clustSeq.get_median_rho_and_sigma(jets3, rangeQ1, false, bkg2Q[0], sigma2Q[0], meanarea2Q[0], true);

    //Quadrant2 - orthogonal
    Double_t phiQ2 = phiLeadingTrack + TMath::Pi()/2.;
    if(phiQ2>TMath::TwoPi()) phiQ2 = phiQ2 - TMath::TwoPi();
    phiMin = phiQ2 - phiStep;
    phiMax = phiQ2 + phiStep;
    fastjet::RangeDefinition rangeQ2(rapMin,rapMax, phiMin, phiMax);
    clustSeq.get_median_rho_and_sigma(jets3, rangeQ2, false, bkg2Q[1], sigma2Q[1], meanarea2Q[1], true);

    //Quadrant3 - away side
    Double_t phiQ3 = phiLeadingTrack + TMath::Pi();
    if(phiQ3>TMath::TwoPi()) phiQ3 = phiQ3 - TMath::TwoPi();
    phiMin = phiQ3 - phiStep;
    phiMax = phiQ3 + phiStep;
    fastjet::RangeDefinition rangeQ3(rapMin,rapMax, phiMin, phiMax);
    clustSeq.get_median_rho_and_sigma(jets3, rangeQ3, false, bkg2Q[2], sigma2Q[2], meanarea2Q[2], true);

    //Quadrant4 - othogonal
    Double_t phiQ4 = phiLeadingTrack + 3./2.*TMath::Pi();
    if(phiQ4>TMath::TwoPi()) phiQ4 = phiQ4 - TMath::TwoPi();
    phiMin = phiQ4 - phiStep;
    phiMax = phiQ4 + phiStep;
    fastjet::RangeDefinition rangeQ4(rapMin,rapMax, phiMin, phiMax);
    clustSeq.get_median_rho_and_sigma(jets3, rangeQ4, false, bkg2Q[3], sigma2Q[3], meanarea2Q[3], true);

    fh3CentvsRhoLeadingTrackPtQ1->Fill(cent,bkg2Q[0],leading->Pt());
    fh3CentvsRhoLeadingTrackPtQ2->Fill(cent,bkg2Q[1],leading->Pt());
    fh3CentvsRhoLeadingTrackPtQ3->Fill(cent,bkg2Q[2],leading->Pt());
    fh3CentvsRhoLeadingTrackPtQ4->Fill(cent,bkg2Q[3],leading->Pt());

    fh3CentvsSigmaLeadingTrackPtQ1->Fill(cent,sigma2Q[0],leading->Pt());
    fh3CentvsSigmaLeadingTrackPtQ2->Fill(cent,sigma2Q[1],leading->Pt());
    fh3CentvsSigmaLeadingTrackPtQ3->Fill(cent,sigma2Q[2],leading->Pt());
    fh3CentvsSigmaLeadingTrackPtQ4->Fill(cent,sigma2Q[3],leading->Pt());

    fh3MultvsRhoLeadingTrackPtQ1->Fill(nCh,bkg2Q[0],leading->Pt());
    fh3MultvsRhoLeadingTrackPtQ2->Fill(nCh,bkg2Q[1],leading->Pt());
    fh3MultvsRhoLeadingTrackPtQ3->Fill(nCh,bkg2Q[2],leading->Pt());
    fh3MultvsRhoLeadingTrackPtQ4->Fill(nCh,bkg2Q[3],leading->Pt());

    fh3MultvsSigmaLeadingTrackPtQ1->Fill(nCh,sigma2Q[0],leading->Pt());
    fh3MultvsSigmaLeadingTrackPtQ2->Fill(nCh,sigma2Q[1],leading->Pt());
    fh3MultvsSigmaLeadingTrackPtQ3->Fill(nCh,sigma2Q[2],leading->Pt());
    fh3MultvsSigmaLeadingTrackPtQ4->Fill(nCh,sigma2Q[3],leading->Pt());

    fh3CentvsDeltaRhoLeadingTrackPtQ1->Fill(cent,bkg2Q[0]-bkg2,leading->Pt());
    fh3CentvsDeltaRhoLeadingTrackPtQ2->Fill(cent,bkg2Q[1]-bkg2,leading->Pt());
    fh3CentvsDeltaRhoLeadingTrackPtQ3->Fill(cent,bkg2Q[2]-bkg2,leading->Pt());
    fh3CentvsDeltaRhoLeadingTrackPtQ4->Fill(cent,bkg2Q[3]-bkg2,leading->Pt());

  }
   

  
  
 
  // fill track information
  Int_t nTrackOver = recParticles.GetSize();
  // do the same for tracks and jets

  if(nTrackOver>0){
    TIterator *recIter = recParticles.MakeIterator();
    AliVParticle *tmpRec = (AliVParticle*)(recIter->Next());  
    Float_t pt = tmpRec->Pt();

    //    Printf("Leading track p_t %3.3E",pt);
    for(int i = 1;i <= fh2NRecTracksPt->GetNbinsX();i++){
      Float_t ptCut = fh2NRecTracksPt->GetXaxis()->GetBinCenter(i);
      while(pt<ptCut&&tmpRec){
	nTrackOver--;
	tmpRec = (AliVParticle*)(recIter->Next()); 
	if(tmpRec){
	  pt = tmpRec->Pt();
	}
      }
      if(nTrackOver<=0)break;
      fh2NRecTracksPt->Fill(ptCut,nTrackOver);
    }
   
    recIter->Reset();
    AliVParticle *leading = (AliVParticle*)recParticles.At(0);
    Float_t phi = leading->Phi();
    if(phi<0)phi+=TMath::Pi()*2.;    
    Float_t eta = leading->Eta();
    pt  = leading->Pt();
    while((tmpRec = (AliVParticle*)(recIter->Next()))){
      Float_t tmpPt = tmpRec->Pt();
      Float_t tmpEta = tmpRec->Eta();
      fh1PtTracksRecIn->Fill(tmpPt);
      fh2TrackEtaPt->Fill(tmpEta,tmpPt);
      if(tmpRec==leading){
	fh1PtTracksLeadingRecIn->Fill(tmpPt);
	fh2LeadingTrackEtaPt->Fill(tmpEta,tmpPt);
	continue;
      }
      // correlation
      Float_t tmpPhi =  tmpRec->Phi();
     
      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    
      Float_t dPhi = phi - tmpPhi;
      if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
      if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
      Float_t dEta = eta - tmpRec->Eta();
      fh2TracksLeadingPhiEta->Fill(dPhi,dEta);
      fh2TracksLeadingPhiPt->Fill(dPhi,pt);
      fh2TracksLeadingPhiPtW->Fill(dPhi,pt,tmpPt);
    }  
    delete recIter;
  }

  // find the random jets

  fastjet::ClusterSequenceArea clustSeqRan(inputParticlesRecRan, jetDef, areaDef);
  
  // fill the jet information from random track
  const vector <fastjet::PseudoJet> &inclusiveJetsRan = clustSeqRan.inclusive_jets();
  const vector <fastjet::PseudoJet> &sortedJetsRan    = sorted_by_pt(inclusiveJetsRan);

  fh1NJetsRecRan->Fill(sortedJetsRan.size());

  // loop over all jets an fill information, first on is the leading jet

  Int_t nRecOverRan = inclusiveJetsRan.size();
  Int_t nRecRan     = inclusiveJetsRan.size();

  if(inclusiveJetsRan.size()>0){
    AliAODJet leadingJet (sortedJetsRan[0].px(), sortedJetsRan[0].py(), sortedJetsRan[0].pz(), sortedJetsRan[0].E());
    Float_t pt = leadingJet.Pt();
   
    Int_t iCount = 0;
    TLorentzVector vecarearanb;

    for(int i = 1;i <= fh2NRecJetsPtRan->GetNbinsX();i++){
      Float_t ptCut = fh2NRecJetsPtRan->GetXaxis()->GetBinCenter(i);
      while(pt<ptCut&&iCount<nRecRan){
	nRecOverRan--;
	iCount++;
	if(iCount<nRecRan){
	  pt = sortedJetsRan[iCount].perp();
	}
      }
      if(nRecOverRan<=0)break;
      fh2NRecJetsPtRan->Fill(ptCut,nRecOverRan);
    }
    Float_t phi = leadingJet.Phi();
    if(phi<0)phi+=TMath::Pi()*2.;    
    pt = leadingJet.Pt();

    // correlation of leading jet with random tracks

    for(unsigned int ip = 0; ip < inputParticlesRecRan.size();ip++)
      { 
	Float_t tmpPt = inputParticlesRecRan[ip].perp();
	// correlation
	Float_t tmpPhi =  inputParticlesRecRan[ip].phi();
	if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    
	Float_t dPhi = phi - tmpPhi;
	if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
	if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();      
	fh2TracksLeadingJetPhiPtRan->Fill(dPhi,pt);
	fh2TracksLeadingJetPhiPtWRan->Fill(dPhi,pt,tmpPt);
      }  

    Int_t nAodOutJetsRan = 0;
    AliAODJet *aodOutJetRan = 0;
    for(int j = 0; j < nRecRan;j++){
      AliAODJet tmpRec (sortedJetsRan[j].px(), sortedJetsRan[j].py(), sortedJetsRan[j].pz(), sortedJetsRan[j].E());
      Float_t tmpPt = tmpRec.Pt();
      fh1PtJetsRecInRan->Fill(tmpPt);
      // Fill Spectra with constituents
      const vector<fastjet::PseudoJet> &constituents = clustSeqRan.constituents(sortedJetsRan[j]);
      fh1NConstRecRan->Fill(constituents.size());
      fh2NConstPtRan->Fill(tmpPt,constituents.size());
      fh2PtNchRan->Fill(nCh,tmpPt);
      fh2PtNchNRan->Fill(nCh,tmpPt,constituents.size());


      if(tmpPt>fJetOutputMinPt&&fTCAJetsOutRan){
	aodOutJetRan =  new ((*fTCAJetsOutRan)[nAodOutJetsRan++]) AliAODJet(tmpRec);
	Double_t arearan=clustSeqRan.area(sortedJetsRan[j]);
	aodOutJetRan->GetRefTracks()->Clear();
	aodOutJetRan->SetEffArea(arearan,0);
	fastjet::PseudoJet vecarearan=clustSeqRan.area_4vector(sortedJetsRan[j]);  
	vecarearanb.SetPxPyPzE(vecarearan.px(),vecarearan.py(),vecarearan.pz(),vecarearan.e());
	aodOutJetRan->SetVectorAreaCharged(&vecarearanb);

      }

      // correlation
      Float_t tmpPhi =  tmpRec.Phi();
      if(tmpPhi<0)tmpPhi+=TMath::Pi()*2.;    

      if(j==0){
	fh1PtJetsLeadingRecInRan->Fill(tmpPt);
	fh1NConstLeadingRecRan->Fill(constituents.size());
	fh2NConstLeadingPtRan->Fill(tmpPt,constituents.size());
	continue;
      }
    }  

     
    if(fAODJetBackgroundOut){
      Double_t bkg3=0.;
      Double_t sigma3=0.;
      Double_t meanarea3=0.;
      clustSeqRan.get_median_rho_and_sigma(sortedJetsRan ,range, false, bkg3, sigma3, meanarea3, true);
      fAODJetBackgroundOut->SetBackground(2,bkg3,sigma3,meanarea3);
      //     float areaRandomCone = rRandomCone2 *TMath::Pi();         
      /*
	fh1BiARandomCones[2]->Fill(ptRandomCone-(bkg3*areaRandomCone));    
	fh1BiARandomConesRan[2]->Fill(ptRandomConeRan-(bkg3*areaRandomCone));    
      */
    }



  }


  // do the event selection if activated
  if(fJetTriggerPtCut>0){
    bool select = false;
    Float_t minPt = fJetTriggerPtCut;
    /*
    // hard coded for now ...
    // 54.50 44.5 29.5 18.5 for anti-kt rejection 1E-3
    if(cent<10)minPt = 50;
    else if(cent<30)minPt = 42;
    else if(cent<50)minPt = 28;
    else if(cent<80)minPt = 18;
    */
    float rho = 0;
    if(externalBackground)rho = externalBackground->GetBackground(2);
    if(fTCAJetsOut){
      for(int i = 0;i < fTCAJetsOut->GetEntriesFast();i++){
	AliAODJet *jet = (AliAODJet*)fTCAJetsOut->At(i);
	Float_t ptSub = jet->Pt() - rho *jet->EffectiveAreaCharged();
	if(ptSub>=minPt){
	  select = true;
	  break;
	}
      }
    }   
 
    if(select){
      static AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
      fh1CentralitySelect->Fill(cent);
      fh1ZSelect->Fill(zVtx);
      aodH->SetFillAOD(kTRUE);
    }
  }
  if (fDebug > 2){
    if(fTCAJetsOut)Printf("%s:%d Rec Jets %d",(char*)__FILE__,__LINE__,fTCAJetsOut->GetEntriesFast());
    if(fTCAJetsOutRan)Printf("%s:%d Rec Jets Ran %d",(char*)__FILE__,__LINE__,fTCAJetsOutRan->GetEntriesFast());
    if(fTCARandomConesOut)Printf("%s:%d RC %d",(char*)__FILE__,__LINE__,fTCARandomConesOut->GetEntriesFast());
    if(fTCARandomConesOutRan)Printf("%s:%d RC Ran %d",(char*)__FILE__,__LINE__,fTCARandomConesOutRan->GetEntriesFast());
  }
  PostData(1, fHistList);
}

void AliAnalysisTaskJetCluster::Terminate(Option_t */*option*/)
{
  //
  // Terminate analysis
  //
  if (fDebug > 1) printf("AnalysisJetCluster: Terminate() \n");

  if(fMomResH1Fit) delete fMomResH1Fit;
  if(fMomResH2Fit) delete fMomResH2Fit;
  if(fMomResH3Fit) delete fMomResH3Fit;
    
}


Int_t  AliAnalysisTaskJetCluster::GetListOfTracks(TList *list,Int_t type){

  //
  // get list of tracks/particles for different types
  //

  if(fDebug>2)Printf("%s:%d Selecting tracks with %d",(char*)__FILE__,__LINE__,type);

  Int_t iCount = 0;
  if(type==kTrackAOD || type==kTrackAODextra || type==kTrackAODextraonly || type==kTrackAODMCextra || type==kTrackAODMCextraonly){

    if(type!=kTrackAODextraonly && type!=kTrackAODMCextraonly) {
      AliAODEvent *aod = 0;
      if(fUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
      else aod = AODEvent();
      if(!aod){
	if(fDebug>2)Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
	return iCount;
      }

      for(int it = 0;it < aod->GetNumberOfTracks();++it){
	AliAODTrack *tr = aod->GetTrack(it);
	Bool_t bGood = false;
	if(fFilterType == 0)bGood = true;
	else if(fFilterType == 1)bGood = tr->IsHybridTPCConstrainedGlobal();
	else if(fFilterType == 2)bGood = tr->IsHybridGlobalConstrainedGlobal();
	if((fFilterMask>0)&&((!tr->TestFilterBit(fFilterMask)||(!bGood)))){
	  if(fDebug>10)Printf("%s:%d Not matching filter %d/%d %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks(),fFilterMask,tr->GetFilterMap());	
	  continue;
	}
	if(TMath::Abs(tr->Eta())>fTrackEtaWindow){
	  if(fDebug>10)Printf("%s:%d Not matching eta %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	  continue;
	}
	if(tr->Pt()<fTrackPtCut){
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
      if(fUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
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
	if(TMath::Abs(trackAOD->Eta())>fTrackEtaWindow) continue;
	if(trackAOD->Pt()<fTrackPtCut) continue;
	list->Add(trackAOD);
	iCount++;
      }
    }

    if(type==kTrackAODMCextra || type==kTrackAODMCextraonly) { //embed generator level particles
      AliAODEvent *aod = 0;
      if(fUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
      else aod = AODEvent();
      if(!aod){
	return iCount;
      }
      TClonesArray *aodExtraTracks = dynamic_cast<TClonesArray*>(aod->FindListObject("aodExtraMCparticles"));
      if(!aodExtraTracks)return iCount;
      for(int it =0; it<aodExtraTracks->GetEntries(); it++) {
	AliVParticle *track = dynamic_cast<AliVParticle*> ((*aodExtraTracks)[it]);
	if (!track) continue;

	AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*> (track);
	if(!trackAOD)continue;
        
	if(TMath::Abs(trackAOD->Eta())>fTrackEtaWindow) continue;
	if(trackAOD->Pt()<fTrackPtCut) continue;
	list->Add(trackAOD);
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
	if(part->Pt()<fTrackPtCut)continue;
	list->Add(part);
	iCount++;
      }
      else if(type == kTrackKineCharged){
	if(part->Particle()->GetPDG()->Charge()==0)continue;
	if(part->Pt()<fTrackPtCut)continue;
	list->Add(part);
	iCount++;
      }
    }
  }
  else if (type == kTrackAODMCCharged || type == kTrackAODMCAll || type == kTrackAODMCChargedAcceptance) {
    AliAODEvent *aod = 0;
    if(fUseAODMCInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else aod = AODEvent();
    if(!aod)return iCount;
    TClonesArray *tca = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    for(int it = 0;it < tca->GetEntriesFast();++it){
      AliAODMCParticle *part = (AliAODMCParticle*)(tca->At(it));
      if(!part->IsPhysicalPrimary())continue;
      if(type == kTrackAODMCAll){
	if(part->Pt()<fTrackPtCut)continue;
	list->Add(part);
	iCount++;
      }
      else if (type == kTrackAODMCCharged || type == kTrackAODMCChargedAcceptance ){
	if(part->Charge()==0)continue;
	if(part->Pt()<fTrackPtCut)continue;
	if(kTrackAODMCCharged){
	  list->Add(part);
	}
	else {
	  if(TMath::Abs(part->Eta())>fTrackEtaWindow)continue;
	  list->Add(part);
	}
	iCount++;
      }
    }
  }// AODMCparticle
  list->Sort();
  return iCount;
}

void AliAnalysisTaskJetCluster::LoadTrPtResolutionRootFileFromOADB() {

   if (!gGrid) {
     AliInfo("Trying to connect to AliEn ...");
     TGrid::Connect("alien://");
   }

  TFile *f = TFile::Open(fPathTrPtResolution.Data());

  if(!f)return;

  TProfile *fProfPtPtSigma1PtGlobSt     = (TProfile*)f->Get("fProfPtPtSigma1PtGlobSt");
  TProfile *fProfPtPtSigma1PtGlobCnoITS = (TProfile*)f->Get("fProfPtPtSigma1PtGlobCnoITS");
  TProfile *fProfPtPtSigma1PtGlobCnoSPD = (TProfile*)f->Get("fProfPtPtSigma1PtGlobCnoSPD");

  SetSmearResolution(kTRUE);
  SetMomentumResolutionHybrid(fProfPtPtSigma1PtGlobSt,fProfPtPtSigma1PtGlobCnoITS,fProfPtPtSigma1PtGlobCnoSPD);
  

}

void AliAnalysisTaskJetCluster::LoadTrEfficiencyRootFileFromOADB() {

   if (!gGrid) {
     AliInfo("Trying to connect to AliEn ...");
     TGrid::Connect("alien://");
   }

  TFile *f = TFile::Open(fPathTrEfficiency.Data());
  if(!f)return;

  TH1D *hEffPosGlobSt = (TH1D*)f->Get("hEffPosGlobSt");
  TH1D *hEffPosGlobCnoITS = (TH1D*)f->Get("hEffPosGlobCnoITS");
  TH1D *hEffPosGlobCnoSPD = (TH1D*)f->Get("hEffPosGlobCnoSPD");

  SetDiceEfficiency(kTRUE);

  if(fChangeEfficiencyFraction>0.) {

    TH1D *hEffPosGlobStMin = (TH1D*)hEffPosGlobSt->Clone("hEffPosGlobStMin");
    
    for(int bin=1; bin<=hEffPosGlobSt->GetNbinsX(); bin++) {
      Double_t content = hEffPosGlobSt->GetBinContent(bin);
      hEffPosGlobStMin->SetBinContent(bin,content-fChangeEfficiencyFraction);
    }

    SetEfficiencyHybrid(hEffPosGlobStMin,hEffPosGlobCnoITS,hEffPosGlobCnoSPD);

  }
  else
    SetEfficiencyHybrid(hEffPosGlobSt,hEffPosGlobCnoITS,hEffPosGlobCnoSPD);

}

void AliAnalysisTaskJetCluster::SetMomentumResolutionHybrid(TProfile *p1, TProfile *p2, TProfile *p3) {

  //
  // set mom res profiles
  //

  if(fMomResH1) delete fMomResH1;
  if(fMomResH2) delete fMomResH2;
  if(fMomResH3) delete fMomResH3;

  fMomResH1 = new TProfile(*p1);//(TProfile*)p1->Clone("fMomResH1");
  fMomResH2 = new TProfile(*p2);//(TProfile*)p2->Clone("fMomResH2");
  fMomResH3 = new TProfile(*p3);//(TProfile*)p3->Clone("fMomResH3");

}

void AliAnalysisTaskJetCluster:: SetEfficiencyHybrid(TH1 *h1, TH1 *h2, TH1 *h3) {
  //
  // set tracking efficiency histos
  //

  fhEffH1 = (TH1*)h1->Clone("fhEffH1");
  fhEffH2 = (TH1*)h2->Clone("fhEffH2");
  fhEffH3 = (TH1*)h3->Clone("fhEffH3");
}

Double_t AliAnalysisTaskJetCluster::GetMomentumSmearing(Int_t cat, Double_t pt) {

  //
  // Get smearing on generated momentum
  //

  //printf("GetMomentumSmearing for cat %d and pt = %f \n",cat,pt);

  TProfile *fMomRes = 0x0;
  if(cat==1) fMomRes = (TProfile*)fMomResH1->Clone("fMomRes");
  if(cat==2) fMomRes = (TProfile*)fMomResH2->Clone("fMomRes");
  if(cat==3) fMomRes = (TProfile*)fMomResH3->Clone("fMomRes");

  if(!fMomRes) {
    return 0.;
  }


  Double_t smear = 0.;

  if(pt>20.) {
    if(cat==1 && fMomResH1Fit) smear = fMomResH1Fit->Eval(pt);
    if(cat==2 && fMomResH2Fit) smear = fMomResH2Fit->Eval(pt);
    if(cat==3 && fMomResH3Fit) smear = fMomResH3Fit->Eval(pt);
  }
  else {

    Int_t bin = fMomRes->FindBin(pt);

    smear = fRandom->Gaus(fMomRes->GetBinContent(bin),fMomRes->GetBinError(bin));

  }

  if(fMomRes) delete fMomRes;
  
  return smear;
}

void AliAnalysisTaskJetCluster::FitMomentumResolution() {
  //
  // Fit linear function on momentum resolution at high pT
  //

  if(!fMomResH1Fit && fMomResH1) {
    fMomResH1Fit = new TF1("fMomResH1Fit","[0]+[1]*x",0.,200.);
    fMomResH1->Fit(fMomResH1Fit,"LL V0","",5.,30.);
    fMomResH1Fit ->SetRange(5.,100.);
  }

  if(!fMomResH2Fit && fMomResH2) {
    fMomResH2Fit = new TF1("fMomResH2Fit","[0]+[1]*x",0.,200.);
    fMomResH2->Fit(fMomResH2Fit,"LL V0","",5.,30.);
    fMomResH2Fit ->SetRange(5.,100.);
  }

  if(!fMomResH3Fit && fMomResH3) {
    fMomResH3Fit = new TF1("fMomResH3Fit","[0]+[1]*x",0.,200.);
    fMomResH3->Fit(fMomResH3Fit,"LL V0","",5.,30.);
    fMomResH3Fit ->SetRange(5.,100.);
  }

}

/*
  Int_t AliAnalysisTaskJetCluster::AddParticlesFastJet(TList &particles,vector<fastjet::PseudoJet> &inputParticles){
  for(int i = 0; i < particles.GetEntries(); i++){
  AliVParticle *vp = (AliVParticle*)particles.At(i);
  // Carefull energy is not well determined in real data, should not matter for p_T scheme?
  fastjet::PseudoJet jInp(vp->Px(),vp->Py(),vp->Pz(),vp->E());
  jInp.set_user_index(i);
  inputParticles.push_back(jInp);
  }

  return 0;

  }
*/
