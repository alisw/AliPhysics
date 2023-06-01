/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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
// AliAnalysisTaskSE for Lambdac candidates (V0 candidates)
// and hadrons correlations
//
// Authors:
// Antonio Palasciano,  antonio.palasciano@ba.infn.it 
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
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
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSELcTopK0sCorrelations.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliHFOfflineCorrelator.h"
#include "AliMultSelection.h"
#include "AliAODVZERO.h"
#include "AliESDUtils.h"

// add includes for KFParticle
#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "KFParticleBase.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSELcTopK0sCorrelations)


//________________________________________________________________________
AliAnalysisTaskSELcTopK0sCorrelations::AliAnalysisTaskSELcTopK0sCorrelations():
AliAnalysisTaskSE(),
  fNPtBinsCorr(0), 
  fBinLimsCorr(),
  fPtThreshLow(),
  fPtThreshUp(), 
  fLSBLowLim(), 
  fLSBUppLim(), 
  fRSBLowLim(), 
  fRSBUppLim(),
  fSignLowLim(),
  fSignUppLim(),
  fDaughTrackID(),
  fDaughTrigNum(),
  fEvents(0),
  fAlreadyFilled(kFALSE),
  fNtrigD(0),
  fOutputMass(0),
  fOutputCorr(0),
  fOutputStudy(0),
  fNentries(0), 
  fCutsLambdac(0),
  fCutsTracks(0),
  fCorrelatorTr(0),
  fCorrelatorKc(0),
  fCorrelatorK0(0),
  fReadMC(0),
  fRecoTr(kTRUE),
  fRecoLambdac(kTRUE),
  fSelEvType(kFALSE),
  fMixing(kFALSE),
  fCounter(0),
  fNPtBins(1),
  fFillOnlyLambdacLambdacbar(0),
  fIsSelectedCandidate(0),
  fSys(0),
  fEtaForCorrel(0),
  fIsRejectSDDClusters(0),
  fFillGlobal(kFALSE),
  fMultEv(0.),
  fMultEvOrig(0.),
  fMultEvV0M(0.),
  fMultEvV0MEqual(0.),
  fCentEvV0M(0.),
  fzVtx(0.),
  fMEAxisThresh(kFALSE),
  fSoftPiCut(kFALSE),
  fKaonCorr(kFALSE),
  fSignLeft_LowPt(0),
  fSignRight_LowPt(0),
  fSignLeft_HighPt(0),
  fSignRight_HighPt(0),
  fPoolNum(0),
  fSpeed(kOneBinSB),
  fMergePools(kFALSE),
  fUseLceff(kTRUE),
  fUseTrackeff(kTRUE),
  fPtAssocLimit(1.),
  fMinDPt(2.),
  fV0CentMin(0.),
  fV0CentMax(0.),
  fTrkMultMin(0.),  
  fTrkMultMax(0.),  
  fVsMultAnalysis(kFALSE),
  fFillTrees(kNoTrees),
  fFractAccME(100),
  fAODProtection(0),
  fPurityStudies(kFALSE),
  fUseNtrklWeight(kFALSE),
  fHistNtrklWeight(0x0),
  fWeight(1.),
  fEqualizeTracklets(kFALSE),
  fRefMult(0.),
  fBranchD(),
  fBranchTr(),
  fBranchDCutVars(),
  fTreeD(0x0),
  fTreeTr(0x0),
  fTrackArray(0x0),
  fTrackArrayFilled(kFALSE)     
{
  // Default constructor
  for(Int_t i=0; i<33; i++) fTrackletProfiles[i]=0;
}

//________________________________________________________________________
AliAnalysisTaskSELcTopK0sCorrelations::AliAnalysisTaskSELcTopK0sCorrelations(const char *name,AliRDHFCutsLctoV0* cutsLambdac, AliHFAssociatedTrackCuts* cutsTrk):
  AliAnalysisTaskSE(name),
  fNPtBinsCorr(0),  
  fBinLimsCorr(),
  fPtThreshLow(),
  fPtThreshUp(), 
  fLSBLowLim(), 
  fLSBUppLim(), 
  fRSBLowLim(), 
  fRSBUppLim(),
  fSignLowLim(),
  fSignUppLim(),
  fDaughTrackID(),
  fDaughTrigNum(),
  fEvents(0),
  fAlreadyFilled(kFALSE),
  fNtrigD(0),
  fOutputMass(0),
  fOutputCorr(0),
  fOutputStudy(0),
  fNentries(0),
  fCutsLambdac(0),
  fCutsTracks(cutsTrk),
  fCorrelatorTr(0),
  fCorrelatorKc(0),
  fCorrelatorK0(0),
  fReadMC(0),
  fRecoTr(kTRUE),
  fRecoLambdac(kTRUE),
  fSelEvType(kFALSE),
  fMixing(kFALSE),
  fCounter(0),
  fNPtBins(1),
  fFillOnlyLambdacLambdacbar(0),
  fIsSelectedCandidate(0),
  fSys(0),
  fEtaForCorrel(0),
  fIsRejectSDDClusters(0),
  fFillGlobal(kFALSE),
  fMultEv(0.),
  fMultEvOrig(0.),
  fMultEvV0M(0.),
  fMultEvV0MEqual(0.),
  fCentEvV0M(0.),
  fzVtx(0.),
  fMEAxisThresh(kFALSE),
  fSoftPiCut(kFALSE),
  fKaonCorr(kFALSE),
  fSignLeft_LowPt(0),
  fSignRight_LowPt(0),
  fSignLeft_HighPt(0),
  fSignRight_HighPt(0),
  fPoolNum(0),
  fSpeed(kOneBinSB),
  fMergePools(kFALSE),
  fUseLceff(kTRUE),
  fUseTrackeff(kTRUE),
  fPtAssocLimit(1.),
  fMinDPt(2.),
  fV0CentMin(0.),
  fV0CentMax(0.),
  fTrkMultMin(0.),  
  fTrkMultMax(0.),  
  fVsMultAnalysis(kFALSE),
  fFillTrees(kNoTrees),
  fFractAccME(100),
  fAODProtection(0),
  fPurityStudies(kFALSE),
  fUseNtrklWeight(kFALSE),
  fHistNtrklWeight(0x0),
  fWeight(1.),
  fEqualizeTracklets(kFALSE),
  fRefMult(0.),
  fBranchD(),
  fBranchTr(),
  fBranchDCutVars(),
  fTreeD(0x0),
  fTreeTr(0x0),
  fTrackArray(0x0),
  fTrackArrayFilled(kFALSE)        
{
  // Default constructor

  fNPtBins=cutsLambdac->GetNPtBins();
    
  fCutsLambdac=cutsLambdac;

  for(Int_t i=0; i<33; i++) fTrackletProfiles[i]=0;

  // Output slot #1 writes into a TList container (mass with cuts)
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TH1F container (number of events)
  DefineOutput(2,TH1F::Class());  //My private output
  // Output slot #3 writes into a AliRDHFCutsLctoV0 container (cuts)
  DefineOutput(3,AliRDHFCutsLctoV0::Class());  //My private output
  // Output slot #4 writes Normalization Counter 
  DefineOutput(4,AliNormalizationCounter::Class());
  // Output slot #5 writes into a TList container (correl output)
  DefineOutput(5,TList::Class());  //My private output
  // Output slot #6 writes into a TList container (correl advanced)
  DefineOutput(6,TList::Class());  //My private output
  // Output slot #7 writes into a AliHFAssociatedTrackCuts container (cuts)
  DefineOutput(7,AliHFAssociatedTrackCuts::Class());  //My private output
  // Output slot #8 writes into a TTree (Lambdac)
  DefineOutput(8,TTree::Class());  //My private output
  // Output slot #9 writes into a TTree (Tracks)
  DefineOutput(9,TTree::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSELcTopK0sCorrelations::AliAnalysisTaskSELcTopK0sCorrelations(const AliAnalysisTaskSELcTopK0sCorrelations &source):
  AliAnalysisTaskSE(source),
  fNPtBinsCorr(source.fNPtBinsCorr), 
  fBinLimsCorr(source.fBinLimsCorr),
  fPtThreshLow(source.fPtThreshLow),
  fPtThreshUp(source.fPtThreshUp), 
  fLSBLowLim(source.fLSBLowLim), 
  fLSBUppLim(source.fLSBUppLim), 
  fRSBLowLim(source.fRSBLowLim), 
  fRSBUppLim(source.fRSBUppLim),
  fSignLowLim(source.fSignLowLim),
  fSignUppLim(source.fSignUppLim),
  fDaughTrackID(source.fDaughTrackID),
  fDaughTrigNum(source.fDaughTrigNum),
  fEvents(source.fEvents),
  fAlreadyFilled(source.fAlreadyFilled),
  fNtrigD(source.fNtrigD),
  fOutputMass(source.fOutputMass),
  fOutputCorr(source.fOutputCorr),
  fOutputStudy(source.fOutputStudy),
  fNentries(source.fNentries), 
  fCutsLambdac(source.fCutsLambdac),
  fCutsTracks(source.fCutsTracks),
  fCorrelatorTr(source.fCorrelatorTr),
  fCorrelatorKc(source.fCorrelatorKc),
  fCorrelatorK0(source.fCorrelatorK0),
  fReadMC(source.fReadMC),
  fRecoTr(source.fRecoTr),
  fRecoLambdac(source.fRecoLambdac),
  fSelEvType(source.fSelEvType),
  fMixing(source.fMixing),
  fCounter(source.fCounter),
  fNPtBins(source.fNPtBins),
  fFillOnlyLambdacLambdacbar(source.fFillOnlyLambdacLambdacbar),
  fIsSelectedCandidate(source.fIsSelectedCandidate),
  fSys(source.fSys),
  fEtaForCorrel(source.fEtaForCorrel),
  fIsRejectSDDClusters(source.fIsRejectSDDClusters),
  fFillGlobal(source.fFillGlobal),
  fMultEv(source.fMultEv),
  fMultEvOrig(source.fMultEvOrig),
  fMultEvV0M(source.fMultEvV0M),
  fMultEvV0MEqual(source.fMultEvV0MEqual),
  fCentEvV0M(source.fCentEvV0M),
  fzVtx(source.fzVtx),
  fSoftPiCut(source.fSoftPiCut),
  fMEAxisThresh(source.fMEAxisThresh),
  fKaonCorr(source.fKaonCorr),
  fSignLeft_LowPt(source.fSignLeft_LowPt),
  fSignRight_LowPt(source.fSignRight_LowPt),
  fSignLeft_HighPt(source.fSignLeft_HighPt),
  fSignRight_HighPt(source.fSignRight_HighPt),
  fPoolNum(source.fPoolNum),
  fSpeed(source.fSpeed),
  fMergePools(source.fMergePools),
  fUseLceff(source.fUseLceff),
  fUseTrackeff(source.fUseTrackeff),
  fPtAssocLimit(source.fPtAssocLimit),
  fMinDPt(source.fMinDPt),
  fV0CentMin(source.fV0CentMin),
  fV0CentMax(source.fV0CentMax),
  fTrkMultMin(source.fTrkMultMin),
  fTrkMultMax(source.fTrkMultMax),
  fVsMultAnalysis(source.fVsMultAnalysis),
  fFillTrees(source.fFillTrees),
  fFractAccME(source.fFractAccME),
  fAODProtection(source.fAODProtection),
  fPurityStudies(source.fPurityStudies),
  fUseNtrklWeight(source.fUseNtrklWeight),
  fHistNtrklWeight(source.fHistNtrklWeight),
  fWeight(source.fWeight),
  fEqualizeTracklets(source.fEqualizeTracklets),
  fRefMult(source.fRefMult),
  fBranchD(source.fBranchD),
  fBranchTr(source.fBranchTr),
  fBranchDCutVars(source.fBranchDCutVars), 
  fTreeD(source.fTreeD),
  fTreeTr(source.fTreeTr),  
  fTrackArray(source.fTrackArray),
  fTrackArrayFilled(source.fTrackArrayFilled)   
{
  // Copy constructor
  for(Int_t i=0; i<33; i++) fTrackletProfiles[i]=source.fTrackletProfiles[i];
}

//________________________________________________________________________
AliAnalysisTaskSELcTopK0sCorrelations::~AliAnalysisTaskSELcTopK0sCorrelations()
{
  if (fOutputMass) {
    delete fOutputMass;
    fOutputMass = 0;
  }
  if (fOutputCorr) {
    delete fOutputCorr;
    fOutputCorr = 0;
  }
  if (fOutputStudy) {
    delete fOutputStudy;
    fOutputStudy = 0;
  }
  if (fCutsLambdac) {
    delete fCutsLambdac;
    fCutsLambdac = 0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
  if (fCorrelatorTr) {
    delete fCorrelatorTr;
    fCorrelatorTr = 0;
  }
  if (fCorrelatorKc) {
    delete fCorrelatorKc;
    fCorrelatorKc = 0;
  }
  if (fCorrelatorK0) {
    delete fCorrelatorK0;
    fCorrelatorK0 = 0;
  }
  if (fCounter){
    delete fCounter;
    fCounter=0;
  }
  for(Int_t i=0; i<33; i++) {
    if (fTrackletProfiles[i]) delete fTrackletProfiles[i];
  }
  if (fApplyML && fMLResponse)
    delete fMLResponse;

  if (fCreateMLtree && fMLhandler)
    delete fMLhandler;
}  

//______________________________________________________________________________
AliAnalysisTaskSELcTopK0sCorrelations& AliAnalysisTaskSELcTopK0sCorrelations::operator=(const AliAnalysisTaskSELcTopK0sCorrelations& orig)
{
// Assignment
  if (&orig == this) return *this; //if address is the same (same object), returns itself

  AliAnalysisTaskSE::operator=(orig); //Uses the AliAnalysisTaskSE operator to assign the inherited part of the class
  fNPtBinsCorr = orig.fNPtBinsCorr; 
  fBinLimsCorr = orig.fBinLimsCorr;
  fPtThreshLow = orig.fPtThreshLow;
  fPtThreshUp = orig.fPtThreshUp; 
  fLSBLowLim = orig.fLSBLowLim; 
  fLSBUppLim = orig.fLSBUppLim; 
  fRSBLowLim = orig.fRSBLowLim;  
  fRSBUppLim = orig.fRSBUppLim; 
  fSignLowLim = orig.fSignLowLim;
  fSignUppLim = orig.fSignUppLim;
  fDaughTrackID = orig.fDaughTrackID;
  fDaughTrigNum = orig.fDaughTrigNum;
  fEvents = orig.fEvents;
  fAlreadyFilled = orig.fAlreadyFilled;
  fNtrigD = orig.fNtrigD;
  fOutputMass = orig.fOutputMass;
  fOutputCorr = orig.fOutputCorr;
  fOutputStudy = orig.fOutputStudy;
  fNentries = orig.fNentries; 
  fCutsLambdac = orig.fCutsLambdac;
  fCutsTracks = orig.fCutsTracks;
  fCorrelatorTr = orig.fCorrelatorTr;
  fCorrelatorKc = orig.fCorrelatorKc;
  fCorrelatorK0 = orig.fCorrelatorK0;
  fReadMC = orig.fReadMC;
  fRecoTr = orig.fRecoTr;
  fRecoLambdac = orig.fRecoLambdac;
  fSelEvType = orig.fSelEvType;
  fMixing = orig.fMixing;
  fCounter = orig.fCounter;
  fNPtBins = orig.fNPtBins;
  fFillOnlyLambdacLambdacbar = orig.fFillOnlyLambdacLambdacbar;
  fIsSelectedCandidate = orig.fIsSelectedCandidate;
  fSys = orig.fSys;
  fEtaForCorrel = orig.fEtaForCorrel;
  fIsRejectSDDClusters = orig.fIsRejectSDDClusters;
  fFillGlobal = orig.fFillGlobal;
  fMultEv = orig.fMultEv;
  fMultEvOrig = orig.fMultEvOrig;
  fMultEvV0M = orig.fMultEvV0M;
  fMultEvV0MEqual = orig.fMultEvV0MEqual;
  fCentEvV0M = orig.fCentEvV0M;
  fzVtx = orig.fzVtx;
  fSoftPiCut = orig.fSoftPiCut;
  fMEAxisThresh = orig.fMEAxisThresh;
  fKaonCorr = orig.fKaonCorr;
  fSignLeft_LowPt = orig.fSignLeft_LowPt;
  fSignRight_LowPt = orig.fSignRight_LowPt;
  fSignLeft_HighPt = orig.fSignLeft_HighPt;
  fSignRight_HighPt = orig.fSignRight_HighPt;
  fPoolNum = orig.fPoolNum;
  fSpeed = orig.fSpeed;   
  fMergePools = orig.fMergePools;
  fUseLceff = orig.fUseLceff;
  fUseTrackeff = orig.fUseTrackeff;
  fPtAssocLimit = orig.fPtAssocLimit;
  fMinDPt = orig.fMinDPt;
  fV0CentMin = orig.fV0CentMin;
  fV0CentMax = orig.fV0CentMax;
  fTrkMultMin = orig.fTrkMultMin;
  fTrkMultMax = orig.fTrkMultMax;
  fVsMultAnalysis = orig.fVsMultAnalysis;
  fFillTrees = orig.fFillTrees;
  fFractAccME = orig.fFractAccME; 
  fAODProtection = orig.fAODProtection;
  fPurityStudies = orig.fPurityStudies;
  fUseNtrklWeight = orig.fUseNtrklWeight;
  fHistNtrklWeight = orig.fHistNtrklWeight;
  fWeight = orig.fWeight;
  fEqualizeTracklets = orig.fEqualizeTracklets;
  fRefMult = orig.fRefMult;
  fBranchD = orig.fBranchD;
  fBranchTr = orig.fBranchTr;
  fBranchDCutVars = orig.fBranchDCutVars;
  fTreeD = orig.fTreeD;
  fTreeTr = orig.fTreeTr;
  fTrackArray = orig.fTrackArray;      
  fTrackArrayFilled = orig.fTrackArrayFilled;
  
  for(Int_t i=0; i<33; i++) {
    fTrackletProfiles[i] = orig.fTrackletProfiles[i];
  }

  return *this; //returns pointer of the class
}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::Init()
{
  // Initialization

  if(fDebug > 1) printf("AliAnalysisTaskSELcTopK0sCorrelations::Init() \n");
  
  if(fUseNtrklWeight && !fReadMC){ AliFatal("Nch weights can only be used in MC mode"); return; }

  //Copy of cuts objects
  AliRDHFCutsLctoV0* copyfCutsLambdac = new AliRDHFCutsLctoV0(*fCutsLambdac);
  const char* nameoutput=GetOutputSlot(3)->GetContainer()->GetName();
  copyfCutsLambdac->SetName(nameoutput);

  fTrackArray = new TObjArray();
  
  //needer to clear completely the objects inside with Clear() method
  // Post the data
  PostData(3,copyfCutsLambdac);
  PostData(7,fCutsTracks);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::UserCreateOutputObjects()
{

  if(fFillTrees==kFillTrees) {

    fBranchD = new AliHFCorrelationBranchD();
    fBranchTr = new AliHFCorrelationBranchTr();

    fTreeD = new TTree("fTreeD","TTree for Lambdac mesons");
    fTreeD->Branch("branchD",&fBranchD);

    fTreeTr = new TTree("fTreeTr","TTree for Associated Tracks");
    fTreeTr->Branch("branchTr",&fBranchTr);

    PostData(8,fTreeD);
    PostData(9,fTreeTr);
  }
  
  if(fFillTrees==kFillCutOptTree) {

    fBranchDCutVars = new AliD0hCutOptim();

    fTreeD = new TTree("fTreeD","TTree for Lambdac mesons - Vars for Cut Optimization");
    fTreeD->Branch("branchD",&fBranchDCutVars);

    PostData(8,fTreeD);
  }  

  // Create the output container
  //
  if(fDebug > 1) {
    printf("AliAnalysisTaskSELcTopK0sCorrelations::UserCreateOutputObjects() \n");
  }
  //HFCorrelator creation and definition
  fCorrelatorTr = new AliHFCorrelator("CorrelatorTr",fCutsTracks,fSys,fCutsLambdac);//fSys=0 use multiplicity, =1 use centrality


  fCorrelatorTr->SetDeltaPhiInterval(-TMath::Pi()/2,3*TMath::Pi()/2);// set the Delta Phi Interval you want (in this case -0.5Pi to 1.5 Pi)


  fCorrelatorTr->SetEventMixing(fMixing);// sets the analysis on a single event (kFALSE) or mixed events (kTRUE)


  fCorrelatorTr->SetAssociatedParticleType(1);// set 1 for correlations with hadrons, 2 with kaons, 3 with KZeros


  fCorrelatorTr->SetApplyDisplacementCut(2); //0: don't calculate Lambdac; 1: return Lambdac; 2: return Lambdac/Lambdacerr


  fCorrelatorTr->SetUseMC(fReadMC);// sets Montecarlo flag


  fCorrelatorTr->SetUseReco(fRecoTr);// sets (if MC analysis) wheter to analyze Reco or Kinem tracks


  if(fMixing && fSoftPiCut) {
    fCorrelatorTr->SetStoreInfoSoftPiME(kTRUE);
  }
 
  Bool_t pooldefTr = fCorrelatorTr->DefineEventPool();// method that defines the properties ot the event mixing (zVtx and Multipl. bins)


  if(!pooldefTr) AliInfo("Warning:: Event pool not defined properly");


  // Several histograms are more conveniently managed in a TList
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("listMass");

  fOutputCorr = new TList();
  fOutputCorr->SetOwner();
  fOutputCorr->SetName("correlationslist");

  fOutputStudy = new TList();
  fOutputStudy->SetOwner();
  fOutputStudy->SetName("controlplots");

  TString nameMass=" ",nameSgn=" ", nameBkg=" ", nameRfl=" ",nameMassWg=" ",nameSgnWg=" ", nameBkgWg=" ", nameRflWg=" ";

//for origin c case (or for data)
  for(Int_t i=0;i<fCutsLambdac->GetNPtBins();i++){

    nameMass="histMass_"; if(fReadMC) nameMass+="c_";
    nameMass+=i;
    nameMassWg="histMass_WeigLambdacEff_"; if(fReadMC) nameMassWg+="c_";
    nameMassWg+=i;
    nameSgn="histSgn_"; if(fReadMC) nameSgn+="c_";
    nameSgn+=i;
    nameSgnWg="histSgn_WeigLambdacEff_"; if(fReadMC) nameSgnWg+="c_";
    nameSgnWg+=i;
    nameBkg="histBkg_"; if(fReadMC) nameBkg+="c_";
    nameBkg+=i;
    nameBkgWg="histBkg_WeigLambdacEff_"; if(fReadMC) nameBkgWg+="c_";
    nameBkgWg+=i;
    nameRfl="histRfl_"; if(fReadMC) nameRfl+="c_";
    nameRfl+=i;
    nameRflWg="histRfl_WeigLambdacEff_"; if(fReadMC) nameRflWg+="c_";
    nameRflWg+=i;

    //histograms of invariant mass distributions

    //MC signal
    if(fReadMC){
      TH1F* tmpSt = new TH1F(nameSgn.Data(), "#Lambda_{c} invariant mass c - MC; M [GeV]; Entries",150,1.9864,2.5864);
      TH1F* tmpStWg = new TH1F(nameSgnWg.Data(), "#Lambda_{c} invariant mass c - MC; M [GeV] - weight 1/Lambdaceff; Entries",150,1.9864,2.5864);
      tmpSt->Sumw2();
      tmpStWg->Sumw2();

      //Reflection: histo filled with LambdacMass which pass the cut (also) as Lambdacbar and with Lambdacbar which pass (also) the cut as Lambdac
      TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass c - MC; M [GeV]; Entries",150,1.9864,2.5864);
      TH1F* tmpRtWg = new TH1F(nameRflWg.Data(), "Reflected signal invariant mass c - MC - weight 1/Lambdaceff; M [GeV]; Entries",150,1.9864,2.5864);
      TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass c - MC; M [GeV]; Entries",150,1.9864,2.5864);
      TH1F* tmpBtWg = new TH1F(nameBkgWg.Data(), "Background invariant mass c - MC - weight 1/Lambdaceff; M [GeV]; Entries",150,1.9864,2.5864);
      tmpBt->Sumw2();
      tmpBtWg->Sumw2();
      tmpRt->Sumw2();
      tmpRtWg->Sumw2();
      fOutputMass->Add(tmpSt);
      fOutputMass->Add(tmpStWg);
      fOutputMass->Add(tmpRt);
      fOutputMass->Add(tmpRtWg);
      fOutputMass->Add(tmpBt);
      fOutputMass->Add(tmpBtWg);
    }

    //mass
    TH1F* tmpMt = new TH1F(nameMass.Data(),"#Lambda_{c} invariant mass c; M [GeV]; Entries",150,1.9864,2.5864);
    tmpMt->Sumw2();
    fOutputMass->Add(tmpMt);
    //mass weighted by 1/Lambdaceff
    TH1F* tmpMtwg = new TH1F(nameMassWg.Data(),"#Lambda_{c} invariant mass c - weight 1/Lambdaceff; M [GeV]; Entries",150,1.9864,2.5864);
    tmpMtwg->Sumw2();
    fOutputMass->Add(tmpMtwg);
    
    if(fFillTrees>0) { //multi-histo for mass, pT, centrality for offline code (use in place of TH1F to select centrality and change pT bins offline)
      nameMass="histMass2D_";  nameMass+=i;
      TH2F* hMass2D = new TH2F(nameMass.Data(),"Mass histogram vs centrality; Entries",150,1.9864,2.5864,100,0.,100.);
      hMass2D->Sumw2();
      fOutputMass->Add(hMass2D);

      nameMass="histMass2D_WeigLambdacEff_";  nameMass+=i;      
      TH2F* hMass2DW = new TH2F(nameMass.Data(),"Mass histogram vs centrality - weight 1/Lambdaceff; Entries",150,1.9864,2.5864,100,0.,100.);
      hMass2DW->Sumw2();
      fOutputMass->Add(hMass2DW);      
    }
    
  }

//for origin b case (no Bkg and Mass histos, here for weights you should use c+b efficiencies, while on data (on MC they're useless))
  for(Int_t i=0;i<fCutsLambdac->GetNPtBins();i++){

    nameSgn="histSgn_b_";
    nameSgn+=i;
    nameSgnWg="histSgn_WeigLambdacEff_b_";
    nameSgnWg+=i;
    nameRfl="histRfl_b_";
    nameRfl+=i;
    nameRflWg="histRfl_WeigLambdacEff_b_";
    nameRflWg+=i;

    //histograms of invariant mass distributions

    //MC signal
    if(fReadMC){
      TH1F* tmpSt = new TH1F(nameSgn.Data(), "#Lambda_{c} invariant mass b - MC; M [GeV]; Entries",150,1.9864,2.5864);
      TH1F* tmpStWg = new TH1F(nameSgnWg.Data(), "#Lambda_{c} invariant mass b - MC; M [GeV] - weight 1/Lambdaceff; Entries",150,1.9864,2.5864);
      tmpSt->Sumw2();
      tmpStWg->Sumw2();

      //Reflection: histo filled with LambdacMass which pass the cut (also) as Lambdacbar and with Lambdacbar which pass (also) the cut as Lambdac
      TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass b - MC; M [GeV]; Entries",150,1.9864,2.5864);
      TH1F* tmpRtWg = new TH1F(nameRflWg.Data(), "Reflected signal invariant mass b - MC - weight 1/Lambdaceff; M [GeV]; Entries",150,1.9864,2.5864);
      tmpRt->Sumw2();
      tmpRtWg->Sumw2();
      fOutputMass->Add(tmpSt);
      fOutputMass->Add(tmpStWg);
      fOutputMass->Add(tmpRt);
      fOutputMass->Add(tmpRtWg);
    }
  }

  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();

  fNentries=new TH1F(nameoutput, "Control plot", 22,-0.5,21.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nEventsSelected");
  fNentries->GetXaxis()->SetBinLabel(3,"nEventsGoodVtxPrim");
  fNentries->GetXaxis()->SetBinLabel(4,"mismatch AOD/dAOD");
  fNentries->GetXaxis()->SetBinLabel(5,"REJ: no prim vtx");
  fNentries->GetXaxis()->SetBinLabel(6,"REJ: Pile-up");  
  fNentries->GetXaxis()->SetBinLabel(7,"REJ: centrality");
  fNentries->GetXaxis()->SetBinLabel(8,"REJ: cent wrong ext");
  fNentries->GetXaxis()->SetBinLabel(9,"REJ: cent flatten");
  fNentries->GetXaxis()->SetBinLabel(10,"REJ: trigger class");
  fNentries->GetXaxis()->SetBinLabel(11,"REJ: trigger mask");
  fNentries->GetXaxis()->SetBinLabel(12,"REJ: zVtx>10cm");
  fNentries->GetXaxis()->SetBinLabel(13,"N. of 0SMH");    
  if(fIsRejectSDDClusters) fNentries->GetXaxis()->SetBinLabel(14,"SDD-Cls Rej");
  if(fReadMC) fNentries->GetXaxis()->SetBinLabel(15,"nEvsWithProdMech");
  fNentries->GetXaxis()->SetBinLabel(16,"Lambdac failed to be filled");
  fReadMC ? fNentries->GetXaxis()->SetBinLabel(17,"nTrueLambdacSelected(MC)") : fNentries->GetXaxis()->SetBinLabel(17,"SigmaC<-Lambdac");
  fNentries->GetXaxis()->SetBinLabel(18,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(19,"Preselect Rejection");
  if(fSys==0) fNentries->GetXaxis()->SetBinLabel(20,"nCandSel(QualTr)");
  fNentries->GetXaxis()->SetBinLabel(21,"nCandSel(Cuts)");
  fNentries->GetXaxis()->SetBinLabel(22,"REJ: V0fly");

  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(4)->GetContainer()->GetName()));
  fCounter->Init();

//Loading of ML models
  if (fApplyML)
  {

    fMLResponse = new AliHFMLResponseLambdactopK0s("LambdactopK0sMLResponse", "LambdactopK0sMLResponse", fConfigPath.Data());
    fMLResponse->MLResponseInit();
  }

  CreateCorrelationsObjs(); //creates histos for correlations analysis

  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fNentries);
  PostData(4,fCounter);
  PostData(5,fOutputCorr);
  PostData(6,fOutputStudy);
  PostData(8,fTreeD);
  PostData(9,fTreeTr);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  //cout<<"I'm in UserExec"<<endl;


  //cuts order
  //     printf("    |M-MLambdac| [GeV]    < %f\n",fLambdactopK0s[0]);
  //     printf("    dca    [cm]  < %f\n",fLambdactopK0s[1]);
  //     printf("    cosThetaStar     < %f\n",fLambdactopK0s[2]);
  //     printf("    pTK     [GeV/c]    > %f\n",fLambdactopK0s[3]);
  //     printf("    pTpi    [GeV/c]    > %f\n",fLambdactopK0s[4]);
  //     printf("    |LambdacK|  [cm]  < %f\n",fLambdactopK0s[5]);
  //     printf("    |Lambdacpi| [cm]  < %f\n",fLambdactopK0s[6]);
  //     printf("    LambdacLambdac  [cm^2] < %f\n",fLambdactopK0s[7]);
  //     printf("    cosThetaPoint    > %f\n",fLambdactopK0s[8]);
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  fEvents++;

  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fNentries->Fill(3);
      return;
    }
  }

  TString bname="CascadesHF";

  TClonesArray *inputArray=0;

  fMultEv = 0.; //reset event multiplicity
  fMultEvOrig = 0.; //reset event multiplicity
  fMultEvV0M = 0.;
  fCentEvV0M = 0.;
  fzVtx = 0.; //reset event multiplicity
  fPoolNum = 0; //reset event pool

  fDaughTrackID.clear(); //removes daugher IDs from previous event
  fDaughTrigNum.clear(); //removes daugher trigger matchings from previous event
  fTrackArray->Clear(); //removes associated tracks selected from previous event
  fTrackArrayFilled = kFALSE; //associated track array is now not filled
  
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
    }
  } else if(aod) {
    inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  }

  if(!inputArray || !aod) {
    printf("AliAnalysisTaskSELcTopK0sCorrelations::UserExec: input branch not found!\n");
    return;
  }

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSELcTopK0sCorrelations::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSELcTopK0sCorrelations::UserExec: MC header branch not found!\n");
      return;
    }
  }

  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fCutsLambdac,fReadMC); 

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(12);

  //Call IsEventSelected only for Reco! (and Data of course)
  if(fRecoLambdac && !fCutsLambdac->IsEventSelected(aod)) {
    if(fCutsLambdac->GetWhyRejection()==0) fNentries->Fill(4); // no prim vertex
    if(fCutsLambdac->GetWhyRejection()==1) fNentries->Fill(5); // rejected for pileup
    if(fCutsLambdac->GetWhyRejection()==2) fNentries->Fill(6); // rejected for centrality
    if(fCutsLambdac->GetWhyRejection()==3) fNentries->Fill(7); // rejected for centrality wrong extim
    if(fCutsLambdac->GetWhyRejection()==4) fNentries->Fill(8); // rejected for centrality flattening
    if(fCutsLambdac->GetWhyRejection()==5) fNentries->Fill(9); // rejected for trigger class
    if(fCutsLambdac->GetWhyRejection()==6) fNentries->Fill(11); // rejected for vtx outside 10
    if(fCutsLambdac->GetWhyRejection()==7) fNentries->Fill(10); // rejected for trigger mask
    return;
  }

  //On Kine, instead of IsEventSelected just select on zVtx and trigger mask in pPb
  if(!fRecoLambdac) {

    Double_t zVtxMC = mcHeader->GetVtxZ();

    if(TMath::Abs(zVtxMC)>10) return;

    if(aod->GetTriggerMask()==0 && (aod->GetRunNumber()>=195344 && aod->GetRunNumber()<=195677)) return;
  }

  fNentries->Fill(1); //event selected after selection

  //Setting PIDResponse for associated tracks
  fCorrelatorTr->SetPidAssociated();

  //Selection on production type (MC)
  if(fReadMC && fSelEvType){ 

    Bool_t isMCeventgood = kFALSE;
            
    Int_t eventType = mcHeader->GetEventType();
    Int_t NMCevents = fCutsTracks->GetNofMCEventType();
               
    for(Int_t k=0; k<NMCevents; k++){
      Int_t * MCEventType = fCutsTracks->GetMCEventType();
          
      if(eventType == MCEventType[k]) isMCeventgood= kTRUE;
      ((TH1D*)fOutputStudy->FindObject("EventTypeMC"))->Fill(eventType);
    }
                
    if(NMCevents && !isMCeventgood){
      if(fDebug > 2) std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
      return; 
    }
    fNentries->Fill(14); //event with particular production type                
  
  } //end of selection


  // Check the Nb of SDD clusters
  if (fIsRejectSDDClusters) { 
    Bool_t skipEvent = kFALSE;
    Int_t ntracks = 0;
    if (aod) ntracks = aod->GetNumberOfTracks();
    for(Int_t itrack=0; itrack<ntracks; itrack++) { // loop on tacks
      //    ... get the track
      AliAODTrack * track = (AliAODTrack*)aod->GetTrack(itrack);
      if(!track) {
	AliWarning("Error in casting to AOD track. Not a standard AOD?");
        continue;
      }
      if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
	skipEvent=kTRUE;
	fNentries->Fill(13);
	break;
      }
    }
    if (skipEvent) return;
  }

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  //Fill Event Multiplicity/Centrality
  if(fVsMultAnalysis) { //for v2 analysis (V0M estimator)
    AliMultSelection *multSel = (AliMultSelection*)aod->FindListObject("MultSelection");
    if(!multSel) AliWarning("AliMultSelection object not found!");
    else fCentEvV0M = multSel->GetMultiplicityPercentile("V0M"); //this is the line for the V0M centrality as from DPG!!
    Int_t vzeroMultA=0, vzeroMultC=0; //all the following part comes from D2H vs mult codes
    Int_t vzeroMultAEq=0, vzeroMultCEq=0;
    AliAODVZERO *vzeroAOD = (AliAODVZERO*)aod->GetVZEROData();
    if(vzeroAOD) {
      vzeroMultA = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
      vzeroMultC = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
      vzeroMultAEq = static_cast<Int_t>(AliVertexingHFUtils::GetVZEROAEqualizedMultiplicity(aod));
      vzeroMultCEq = static_cast<Int_t>(AliVertexingHFUtils::GetVZEROCEqualizedMultiplicity(aod));
      fMultEvV0M = vzeroMultA + vzeroMultC;
      fMultEvV0MEqual = vzeroMultAEq + vzeroMultCEq;
    }
    fMultEvOrig = (Double_t)(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.)); //still the uncorrected value of tracklets
  } else { //std behaviour
    if(fSys==0) fMultEvOrig = (Double_t)(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.)); //pp (or pPb)
    else fMultEvOrig = fCutsLambdac->GetCentrality(aod); //PbPb
  }

  //tracklet equalisation, if requested
  if(fEqualizeTracklets) {
    TProfile* estimatorAvg = GetEstimatorHistogram(aod);
    if(!estimatorAvg) {
      printf("Exiting due to issues in tracklet mult correction\n");
      return;
    }
    fMultEv = AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,fMultEvOrig,vtx1->GetZ(),fRefMult);
    if(fDebug > 1) printf("Original tracklet multiplicity is %.2f --> corrected to %.2f\n",fMultEvOrig,fMultEv);
  } else { //don't equalize tracklets
    fMultEv = fMultEvOrig;
  }
  if (fSys!=0) fMultEv = fMultEvOrig; //this variable is the centrality, for PbPb!

  //Fill control plots for event centrality and zVtx
  if(fCutsLambdac->GetUseCentrality()) ((TH1F*)fOutputStudy->FindObject("hCentralEvts"))->Fill(fCutsLambdac->GetCentrality(aod));
  ((TH1F*)fOutputStudy->FindObject("hCentEvV0M"))->Fill(fCentEvV0M);
  ((TH1F*)fOutputStudy->FindObject("hMultEvV0M"))->Fill(fMultEvV0M);
  ((TH1F*)fOutputStudy->FindObject("hMultEvV0MEqual"))->Fill(fMultEvV0MEqual);
  ((TH1F*)fOutputStudy->FindObject("hMultEvTrkl1"))->Fill(fMultEvOrig);
  ((TH1F*)fOutputStudy->FindObject("hZvtxEvts"))->Fill(vtx1->GetZ());
  if(fEqualizeTracklets) {
    ((TH1F*)fOutputStudy->FindObject("hMultEvTrkl1Equal"))->Fill(fMultEv);
    ((TH2F*)fOutputStudy->FindObject("hNtrVsZvtx"))->Fill(vtx1->GetZ(),fMultEvOrig);   
    ((TH2F*)fOutputStudy->FindObject("hNtrCorrVsZvtx"))->Fill(vtx1->GetZ(),fMultEv); 
  }
  
  //Select Centrality range for V2 in pp analysis
  if(fVsMultAnalysis) {
    if((fV0CentMin!=0 || fV0CentMax!=0) && (fCentEvV0M < fV0CentMin || fCentEvV0M >= fV0CentMax)) { //rejected by V0M centrality out of range
    	fNentries->Fill(6); 
    	return;
    } else if((fTrkMultMin!=0 || fTrkMultMax!=0) && (fMultEv < fTrkMultMin || fMultEv >= fTrkMultMax)) { //rejected by SPD tracklet centrality out of range
      fNentries->Fill(6); 
      return;
    } else { //accepted events. Fill histograms and continue running
      ((TH1F*)fOutputStudy->FindObject("hCentEvV0MSelEvents"))->Fill(fCentEvV0M);
      ((TH1F*)fOutputStudy->FindObject("hMultEvV0MSelEvents"))->Fill(fMultEvV0M);
      ((TH1F*)fOutputStudy->FindObject("hMultEvV0MEqualSelEvents"))->Fill(fMultEvV0MEqual);
      ((TH1F*)fOutputStudy->FindObject("hMultEvTrkl1SelEvents"))->Fill(fMultEvOrig);
      if(fEqualizeTracklets) {
        ((TH1F*)fOutputStudy->FindObject("hMultEvTrkl1EqualSelEvents"))->Fill(fMultEv);    
        ((TH2F*)fOutputStudy->FindObject("hNtrVsZvtxSelEvents"))->Fill(vtx1->GetZ(),fMultEvOrig);   
        ((TH2F*)fOutputStudy->FindObject("hNtrCorrVsZvtxSelEvents"))->Fill(vtx1->GetZ(),fMultEv);   
      }
    }
  }
 
  //HFCorrelators initialization (for this event)
  fCorrelatorTr->SetAODEvent(aod); // set the AOD event from which you are processing

  Bool_t correlatorONTr = fCorrelatorTr->Initialize(); // initialize the pool for event mixing

  if(!correlatorONTr) {AliInfo("AliHFCorrelator (tracks) didn't initialize the pool correctly or processed a bad event"); return;}
  if(fReadMC) {
    fCorrelatorTr->SetMCArray(mcArray); // set the TClonesArray *fmcArray for analysis on monte carlo


  }

  //Pool definition
  Double_t MultipOrCent = fCorrelatorTr->GetCentrality();
  Double_t zVtxPosition = vtx1->GetZ();
  fzVtx = zVtxPosition;
  if(!fMergePools) fPoolNum = fCutsTracks->GetPoolBin(MultipOrCent, zVtxPosition);

  //vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    fNentries->Fill(2);
  }

  //Reset flag for tracks distributions fill and counter of Lambdac triggers and of Soft pions
  fAlreadyFilled=kFALSE;
  fNtrigD=0;

  //Reset (and, in case, evaluate), Ntrkl event weights
  fWeight=1.;
  if(fReadMC && fUseNtrklWeight) {
    Int_t nTracklets = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.));
    fWeight *= GetNtrklWeight(nTracklets);
    if(fDebug > 1) printf("Using Ntrkl weights, tracklets=%d, Weight=%f\n",nTracklets,fWeight);
  }

  //***** Loop over Lambdac candidates *****
  Int_t nInLambdactopK0s = inputArray->GetEntriesFast();
  if(fDebug>2) printf("Number of Lambdac->pK0s: %d\n",nInLambdactopK0s);

  Int_t nSelectedloose=0,nSelectedtight=0;  

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

  //RecoLambdac case ************************************************
  if(fRecoLambdac) {

    for (Int_t iLambdactopK0s = 0; iLambdactopK0s < nInLambdactopK0s; iLambdactopK0s++) {
      AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF*)inputArray->UncheckedAt(iLambdactopK0s);

      //new preselection check, from PbPb 2018 (to spped up analysis)
      TObjArray arrTracks(2);
      for(Int_t ipr=0;ipr<2;ipr++){
        AliAODTrack *tr=vHF->GetProng(aod,d,ipr);
        arrTracks.AddAt(tr,ipr);
      }

      if(!fCutsLambdac->PreSelect(arrTracks)){
        fNentries->Fill(18);
        continue;
      }

      if(!(vHF->FillRecoCand(aod,d))) {//Fill the data members of the candidate only if they are empty.   
        fNentries->Fill(15); //monitor how often this fails 
        continue;
      }

      if(d->Pt() < fMinDPt) continue; //to save time and merging memory...

      if(d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kLctoV0Cuts)){
  	    fNentries->Fill(16);
  	    continue;      
      }

      if(fCutsLambdac->IsInFiducialAcceptance(d->Pt(),d->Y(4122))) {
        nSelectedloose++;
        nSelectedtight++;      
        if(fSys==0){
  	     if(fCutsLambdac->IsSelected(d,AliRDHFCuts::kTracks,aod)) fNentries->Fill(19);       
        }  
        Int_t ptbin=fCutsLambdac->PtBin(d->Pt());
        if(ptbin==-1) {fNentries->Fill(17); continue;} //out of bounds

        fIsSelectedCandidate=fCutsLambdac->IsSelected(d,AliRDHFCuts::kAll,aod); //Lambdac selected
        if (fIsSelectedCandidate % 2 == 0){
          continue; // v0s: only K0s hypothesis (remove Lambda)
        }

        AliAODv0 *v0part = dynamic_cast<AliAODv0*>(dynamic_cast<AliAODRecoCascadeHF *>(d)->Getv0());

        if (v0part->GetOnFlyStatus()) {
            fNentries->Fill(22);
            continue;
        }

        //Lambdac
        if(!ReconstructKFLc(d)) continue;
        /* Recompute the main properties of the Lc according to KF particle
        printf("Check values with previous\n");
        printf("d->Phi()=%2f, d->Pt()=%.2f, d->Eta()=%.2f, d->Y()=%.2f, d->InvMassLctoK0sP()=%.2f\n", d->Phi(), d->Pt(), d->Eta(), d->Y(4122), d->InvMassLctoK0sP());
        printf("fKFLcPhi=%2f, fKFLcPt=%.2f, fKFLcEta=%.2f, fKFLcRapidity=%.2f, fKFLcInvMass=%.2f\n",fKFLcPhi, fKFLcPt, fKFLcEta, fKFLcRapidity, fKFLcInvMass);
        */
        //Define pT bin checking the KFParticle properties
        ptbin=fCutsLambdac->PtBin(fKFLcPt);
        if(ptbin==-1) {fNentries->Fill(17); continue;} //out of bounds

        Double_t phiLambdac = fCorrelatorTr->SetCorrectPhiRange(fKFLcPhi);
        fCorrelatorTr->SetTriggerParticleProperties(fKFLcPt,phiLambdac,fKFLcEta); // sets the parameters of the trigger particles that are needed

        //variables for ML application
        AliAODPidHF *Pid_HF = nullptr;
        std::vector<Double_t> modelPred{};
        bool isMLsel = true;

        if (fApplyML)
        {
          Pid_HF = fCutsLambdac->GetPidHF();
          isMLsel = fMLResponse->IsSelectedMultiClass(modelPred, d, aod->GetMagneticField(), Pid_HF);
        }

        if(isMLsel)
        {
          if(modelPred.size() == 0)
          {
            modelPred.push_back(-1.);
            modelPred.push_back(-1.);
            modelPred.push_back(-1.);
          }
          else if(modelPred.size() == 1)
          {
            modelPred.push_back(-1.);
            modelPred.push_back(-1.);
          }
          
          FillMassHists(d,mcArray,fCutsLambdac,fOutputMass,aod);
          //Set special Properties for Lc (i.e. soft pion subtraction)
          fCorrelatorTr->SetLcpK0sProperties(d,fIsSelectedCandidate); //sets special properties for Lc

           if(!fReadMC) {
            if (TMath::Abs(d->Eta())<fEtaForCorrel) {
	            if(!fAlreadyFilled && !fFillTrees) ((TH1F*)fOutputStudy->FindObject(Form("hEvtsPerPool_%d",ptbin)))->Fill(fPoolNum+0.5);			
              if(!fMixing && !fAlreadyFilled) {
 	              ((TH1F*)fOutputStudy->FindObject("hZvtx"))->Fill(vtx1->GetZ());
	              ((TH1F*)fOutputStudy->FindObject(Form("hMultiplEvt_Bin%d",ptbin)))->Fill(fMultEv);
              }
   	          if(fFillTrees==kNoTrees) CalculateCorrelations(d); //correlations on real data
	          }
          } else { //correlations on MC -> association of selected Lambdac to MCinfo with MCtruth
            if (TMath::Abs(d->Eta())<fEtaForCorrel) {
              Int_t pdgLctopK0s[2] = {2212, 310};
              Int_t pdgDgK0stoDaughters[2] = {211, 211};
                TClonesArray *arrayMC;
              //Int_t oldlabLambdac = d->MatchToMC(4122,mcArray,3,pdgDgLambdactopK0s); //return MC particle label if the array corresponds to a Lambdac, -1 if not
              //Int_t labLambdac = d->MatchToMC(4122,mcArray,0); //Checks if it is a lambdac and return a labelParticle
              Int_t labLambdac = d->MatchToMC(4122, 310, pdgLctopK0s, pdgDgK0stoDaughters, arrayMC, true);
              if (labLambdac>-1) {
  	            if(!fAlreadyFilled && !fFillTrees) ((TH1F*)fOutputStudy->FindObject(Form("hEvtsPerPool_%d",ptbin)))->Fill(fPoolNum+0.5);
                if(!fMixing && !fAlreadyFilled) {
		            ((TH1F*)fOutputStudy->FindObject("hZvtx"))->Fill(vtx1->GetZ());
                ((TH1F*)fOutputStudy->FindObject(Form("hMultiplEvt_Bin%d",ptbin)))->Fill(fMultEv); //Fill multiplicity histo
                }
	              if(fFillTrees==kNoTrees) CalculateCorrelations(d,labLambdac,mcArray);
	            }
            }
          }
        } //isMLsel loop ends here
      }
    }
  }
  //End RecoLambdac case ************************************************

  //MCKineLambdac case ************************************************
  if(fReadMC && !fRecoLambdac) {

    for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { //Loop over all the tracks of MCArray
      AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
      if (!mcPart) {
        AliWarning("Particle not found in tree, skipping"); 
        continue;
      } 
      if(TMath::Abs(mcPart->GetPdgCode()) == 4122){  // THIS IS A Lambdac
        if (fCutsLambdac->IsInFiducialAcceptance(mcPart->Pt(),mcPart->Y()) ) {
          nSelectedloose++;
          nSelectedtight++; 
          TClonesArray *arrayMC;
          /* if( MatchLcToMC(mcPart, mcArray) == -1){
            AliAODMCParticle* mcDau1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughterLabel(0)));
	          //AliAODMCParticle* mcDau2 = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughterLabel(1)-1));
	          AliAODMCParticle* mcDau3 = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughterLabel(1)));
            cout<<"daughter levels ------"<<mcDau1<<"  "<<mcDau3<<"  "<<endl;
            if(!mcDau1 || !mcDau3 ) continue;
	              cout << "get pdg     "<< mcDau1->GetPdgCode() <<"  "<<mcDau3->GetPdgCode() <<"  "<<endl;
          } */
          Int_t labDau[3] = {-1, -1, -1}; 
          Int_t deca = AliVertexingHFUtils::CheckLcV0bachelorDecay(arrayMC, mcPart, labDau);
          Bool_t daughInAcc= CheckDaugAcc(arrayMC, 3, labDau);
          //if( MatchLcToMC(mcPart, mcArray) > -1){
            if( deca == 1 && daughInAcc == kTRUE){
            //Removal of cases in which Lambdac decay is not in pK0s!
	          /*  if(mcPart->GetNDaughters()!=3) continue;
	          AliAODMCParticle* mcDau1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughterLabel(0)));
	          AliAODMCParticle* mcDau2 = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughterLabel(1)-1));
	          AliAODMCParticle* mcDau3 = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughterLabel(1)));
	          cout<<"daughter levels ------"<<mcDau1<<"  "<<mcDau2<<"  "<<mcDau3<<endl;
	          if(!mcDau1 || !mcDau2 || !mcDau3) continue;
	              cout<<"get pdg     "<<mcDau1->GetPdgCode()<<"  "<<mcDau2->GetPdgCode()<<"  "<<mcDau3->GetPdgCode()<<endl;
	          Int_t pdg1 = TMath::Abs(mcDau1->GetPdgCode());
	          Int_t pdg2 = TMath::Abs(mcDau2->GetPdgCode());
	          Int_t pdg3 = TMath::Abs(mcDau3->GetPdgCode());
            if(!((pdg1 == 2212 && pdg2 == 321 && pdg3 == 211) || (pdg1 == 211 && pdg2 == 321 && pdg3 == 2212))) continue;
            if(TMath::Abs(mcDau1->Eta())>0.8||TMath::Abs(mcDau2->Eta())>0.8) continue;
            //Check momentum conservation (to exclude 4-prong decays with tracks outside y=1.5)
            Double_t p1[3]  = {mcDau1->Px(),mcDau1->Py(),mcDau1->Pz()};
            Double_t p2[3]  = {mcDau2->Px(),mcDau2->Py(),mcDau2->Pz()};
            Double_t p3[3]  = {mcDau3->Px(),mcDau3->Py(),mcDau3->Pz()};
            Double_t pLambdac[3] = {mcPart->Px(),mcPart->Py(),mcPart->Pz()};
            if(TMath::Abs( (p1[0]+p2[0]+p3[0]-pLambdac[0])*(p1[0]+p2[0]+p3[0]-pLambdac[0]) + (p1[1]+p2[1]+p3[1]-pLambdac[1])*(p1[1]+p2[1]+p3[1]-pLambdac[1]) + (p1[2]+p2[2]+p3[2]-pLambdac[2])*(p1[2]+p2[2]+p3[2]-pLambdac[2]) )>0.1) continue; */
            if(fSys==0) fNentries->Fill(19);
            Int_t ptbin=fCutsLambdac->PtBin(mcPart->Pt());
            if(ptbin==-1) {fNentries->Fill(17); continue;} //out of bounds  
  
            //Lambdac infos
            Double_t phiLambdac = fCorrelatorTr->SetCorrectPhiRange(mcPart->Phi());
            fCorrelatorTr->SetTriggerParticleProperties(mcPart->Pt(),phiLambdac,mcPart->Eta()); // sets the parameters of the trigger particles that are needed

            if (TMath::Abs(mcPart->Eta())<fEtaForCorrel) {
              if (mcPart->GetPdgCode()==4122) fIsSelectedCandidate = 1;
    	        else fIsSelectedCandidate = 2;
	            TString fillthis="histSgn_"; 
	            if(CheckLambdacOrigin(mcArray,mcPart)==4) fillthis+="c_";
	            else if(CheckLambdacOrigin(mcArray,mcPart)==5) fillthis+="b_";
              else continue;
              fillthis+=ptbin;
	            ((TH1F*)(fOutputMass->FindObject(fillthis)))->Fill(2.286);
              CalculateCorrelationsMCKine(mcPart,mcArray);
              if(!fMixing) ((TH1F*)fOutputStudy->FindObject(Form("hMultiplEvt_Bin%d",ptbin)))->Fill(fMultEv); //Fill multiplicity histo
            }
          }
        }
      }
    }
  }
  //End MCKineLambdac case ************************************************

  if(fMixing && fFillTrees!=kFillTrees /* && fAlreadyFilled*/) { // update the pool for Event Mixing, if: enabled,  event is ok, at least a SelLambdac found! (fAlreadyFilled's role!)
    Bool_t updatedTr = fCorrelatorTr->PoolUpdate();
    if(!updatedTr) AliInfo("Pool was not updated");
  }
  if(fFillTrees==kFillTrees && fAlreadyFilled) FillTreeTracks(aod);
  
  ((TH1F*)fOutputStudy->FindObject("hNtrUnCorrEvSel"))->Fill(fMultEvOrig); //Fill multiplicity histo
  if(fEqualizeTracklets) ((TH1F*)fOutputStudy->FindObject("hNtrCorrEvSel"))->Fill(fMultEv); //Fill multiplicity histo
  if(fAlreadyFilled) { //there's a selected D candidate in the event
    ((TH1F*)fOutputStudy->FindObject("hNtrUnCorrEvWithCand"))->Fill(fMultEvOrig); //Fill multiplicity histo
    if(fEqualizeTracklets) ((TH1F*)fOutputStudy->FindObject("hNtrCorrEvWithCand"))->Fill(fMultEv); //Fill multiplicity histo
  }
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);  
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);  
  delete vHF;

  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fNentries);
  PostData(4,fCounter);
  PostData(5,fOutputCorr);
  PostData(6,fOutputStudy);
  if(fFillTrees!=kNoTrees) PostData(8,fTreeD); //fill in case kFillTrees or kFillCutOptTree
  if(fFillTrees==kFillTrees) PostData(9,fTreeTr); //fill only in case kFillTrees
  
  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::FillMassHists(AliAODRecoCascadeHF *part, TClonesArray *arrMC, AliRDHFCutsLctoV0* cuts, TList *listout, AliAODEvent *aod) {
  //
  // function used in UserExec to fill mass histograms:
  //

  if(fDebug>2)  cout << "Candidate selected" << endl;

  //Double_t invmassLambdac = part->InvMassLctoK0sP(), invmassLambdacbar = part->InvMassLctoK0sP();
  Double_t invmassLambdac = fKFLcInvMass;
  Int_t ptbin = cuts->PtBin(fKFLcPt);

  TString fillthis="";
  Int_t pdgLctopK0s[2] = {2212, 310};
  Int_t pdgDgK0stoDaughters[2] = {211, 211};
  Int_t labLambdac=-1;
  TClonesArray *arrayMC;
  if (fReadMC) labLambdac = part->MatchToMC(4122, 310, pdgLctopK0s, pdgDgK0stoDaughters, arrayMC, true); //return MC particle label if the array corresponds to a Lambdac, -1 if not (cf. AliAODRecoDecay.cxx)

  //count candidates selected by cuts
  fNentries->Fill(20);
  //count true Lambdac selected by cuts
  if (fReadMC && labLambdac>=0) fNentries->Fill(16);

  //if ((fIsSelectedCandidate==1 || fIsSelectedCandidate==3) && fFillOnlyLambdacLambdacbar<2) { //Lambdac
  if ((fIsSelectedCandidate > 0)) { //Lambdac

    if(fReadMC){ //on MC
      if(labLambdac>=0 && CheckLambdacOrigin(arrMC,(AliAODMCParticle*)arrMC->At(labLambdac))==4) {
  	    AliAODMCParticle *partLambdac = (AliAODMCParticle*)arrMC->At(labLambdac);
	      Int_t pdgLambdac = partLambdac->GetPdgCode();
        if (pdgLambdac==4122){ //Lambdac
          fillthis="histSgn_c_";
          fillthis+=ptbin;
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac);
          fillthis="histSgn_WeigLambdacEff_c_";
          fillthis+=ptbin;
          Double_t effLambdac = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
          if(!fUseLceff || !effLambdac) effLambdac=1.;
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac,1./effLambdac);
        } else{ //it was a Lambdacbar
          fillthis="histRfl_c_";
          fillthis+=ptbin;
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac);
          fillthis="histRfl_WeigLambdacEff_c_";
          fillthis+=ptbin;
          Double_t effLambdac = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
          if(!fUseLceff || !effLambdac) effLambdac=1.;          
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac,1./effLambdac);
        }
      } else if(labLambdac>=0 && CheckLambdacOrigin(arrMC,(AliAODMCParticle*)arrMC->At(labLambdac))==5) {
        AliAODMCParticle *partLambdac = (AliAODMCParticle*)arrMC->At(labLambdac);
        Int_t pdgLambdac = partLambdac->GetPdgCode();
        if (pdgLambdac==4122){ //Lambdac
          fillthis="histSgn_b_";
          fillthis+=ptbin;
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac);
          fillthis="histSgn_WeigLambdacEff_b_";
          fillthis+=ptbin;
          Double_t effLambdac = fCutsTracks->GetTrigWeightB(part->Pt(),fMultEv);
          if(!fUseLceff || !effLambdac) effLambdac=1.;            
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac,1./effLambdac);
  	    } else{ //it was a Lambdacbar
	        fillthis="histRfl_b_";
	        fillthis+=ptbin;
	        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac);
          fillthis="histRfl_WeigLambdacEff_b_";
          fillthis+=ptbin;
          Double_t effLambdac = fCutsTracks->GetTrigWeightB(part->Pt(),fMultEv);
          if(!fUseLceff || !effLambdac) effLambdac=1.;          
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac,1./effLambdac);
	      }
      } else {//background
  	    fillthis="histBkg_c_";
	      fillthis+=ptbin;
	      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac);
        fillthis="histBkg_WeigLambdacEff_c_";
        fillthis+=ptbin;
        Double_t effLambdac = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
        if(!fUseLceff || !effLambdac) effLambdac=1.; 
        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac,1./effLambdac);
      }
    }else{ //on data
      fillthis="histMass_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac);
      fillthis="histMass_WeigLambdacEff_";
      fillthis+=ptbin;
      Double_t effLambdac = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
       
      if(!fUseLceff || !effLambdac) effLambdac=1.; 
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdac,1./effLambdac);
      if(fFillTrees>0) {
	      Double_t centFill = 0.;
	      if(fCutsLambdac->GetUseCentrality()) centFill = fCutsLambdac->GetCentrality(aod);
        ((TH2F*)(listout->FindObject(Form("histMass2D_%d",ptbin))))->Fill(invmassLambdac,centFill);
        ((TH2F*)(listout->FindObject(Form("histMass2D_WeigLambdacEff_%d",ptbin))))->Fill(invmassLambdac,centFill,1./effLambdac);
      }
      
    }
     
  }
  /* if (fIsSelectedCandidate>1 && (fFillOnlyLambdacLambdacbar==0 || fFillOnlyLambdacLambdacbar==2)) { //Lambdacbar

    if(fReadMC){ //on MC
      if(labLambdac>=0 && CheckLambdacOrigin(arrMC,(AliAODMCParticle*)arrMC->At(labLambdac))==4) {
  	    AliAODMCParticle *partLambdac = (AliAODMCParticle*)arrMC->At(labLambdac);
	      Int_t pdgLambdac = partLambdac->GetPdgCode();
	      if (pdgLambdac==-4122){ //Lambdac
	        fillthis="histSgn_c_";
	        fillthis+=ptbin;
	        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar);
          fillthis="histSgn_WeigLambdacEff_c_";
          fillthis+=ptbin;
          Double_t effLambdac = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
          if(!fUseLceff || !effLambdac) effLambdac=1.; 
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar,1./effLambdac);
  	    } else{ //it was a Lambdacbar
	        fillthis="histRfl_c_";
	        fillthis+=ptbin;
	        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar);
          fillthis="histRfl_WeigLambdacEff_c_";
          fillthis+=ptbin;
          Double_t effLambdac = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
          if(!fUseLceff || !effLambdac) effLambdac=1.; 
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar,1./effLambdac);
  	    }
      } else if(labLambdac>=0 && CheckLambdacOrigin(arrMC,(AliAODMCParticle*)arrMC->At(labLambdac))==5) {
  	    AliAODMCParticle *partLambdac = (AliAODMCParticle*)arrMC->At(labLambdac);
	      Int_t pdgLambdac = partLambdac->GetPdgCode();
	      if (pdgLambdac==-4122){ //Lambdac
	        fillthis="histSgn_b_";
	        fillthis+=ptbin;
	        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar);
          fillthis="histSgn_WeigLambdacEff_b_";
          fillthis+=ptbin;
          Double_t effLambdac = fCutsTracks->GetTrigWeightB(part->Pt(),fMultEv);
          if(!fUseLceff || !effLambdac) effLambdac=1.; 
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar,1./effLambdac);
  	    } else{ //it was a Lambdacbar
	        fillthis="histRfl_b_";
	        fillthis+=ptbin;
	        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar);
          fillthis="histRfl_WeigLambdacEff_b_";
          fillthis+=ptbin;
          Double_t effLambdac = fCutsTracks->GetTrigWeightB(part->Pt(),fMultEv);
          if(!fUseLceff || !effLambdac) effLambdac=1.; 
          ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar,1./effLambdac);
	      }
      } else {//background
  	    fillthis="histBkg_c_";
	      fillthis+=ptbin;
	      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar);
        fillthis="histBkg_WeigLambdacEff_c_";
        fillthis+=ptbin;
        Double_t effLambdac = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
        if(!fUseLceff || !effLambdac) effLambdac=1.; 
        ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar,1./effLambdac);
      }
    }else{ //on data
      fillthis="histMass_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar);
      fillthis="histMass_WeigLambdacEff_";
      fillthis+=ptbin;
      Double_t effLambdac = fCutsTracks->GetTrigWeight(part->Pt(),fMultEv);
      if(!fUseLceff || !effLambdac) effLambdac=1.; 
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassLambdacbar,1./effLambdac);
      if(fFillTrees>0) {
	      Double_t centFill = 0.;
	      if(fCutsLambdac->GetUseCentrality()) centFill = fCutsLambdac->GetCentrality(aod);
        ((TH2F*)(listout->FindObject(Form("histMass2D_%d",ptbin))))->Fill(invmassLambdacbar,centFill);
        ((TH2F*)(listout->FindObject(Form("histMass2D_WeigLambdacEff_%d",ptbin))))->Fill(invmassLambdacbar,centFill,1./effLambdac);
      }      
    }
  } */

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSELambdacCorrelations: Terminate() \n");

  fOutputMass = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputMass) {     
    printf("ERROR: fOutputMass not available\n");
    return;
  }

  fNentries = dynamic_cast<TH1F*>(GetOutputData(2));
  
  if(!fNentries){
    printf("ERROR: fNEntries not available\n");
    return;
  }

  fCutsLambdac = dynamic_cast<AliRDHFCutsLctoV0*>(GetOutputData(3));
  if(!fCutsLambdac){
    printf("ERROR: fCuts not available\n");
    return;
  }

  fCounter = dynamic_cast<AliNormalizationCounter*>(GetOutputData(4));    
  if (!fCounter) {
    printf("ERROR: fCounter not available\n");
    return;
  }
  fOutputCorr = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutputCorr) {     
    printf("ERROR: fOutputCorr not available\n");
    return;
  }
  fOutputStudy = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputStudy) {     
    printf("ERROR: fOutputStudy not available\n");
    return;
  }
  fCutsTracks = dynamic_cast<AliHFAssociatedTrackCuts*>(GetOutputData(7));
  if(!fCutsTracks){
    printf("ERROR: fCutsTracks not available\n");
    return;
  }

  return;
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSELcTopK0sCorrelations::CheckLambdacOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
  if(fDebug > 2) printf("AliAnalysisTaskSELcTopK0sCorrelations::CheckLambdacOrigin() \n");
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;

  while (mother > 0){
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcMoth){
      pdgGranma = mcMoth->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcMoth->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if(isQuarkFound) {
    if(isFromB) return 5;
    else return 4;
  }
  else return 1;
}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::CreateCorrelationsObjs() {
//

  TString namePlot = "";

  //These for limits in THnSparse (one per bin, same limits). 
  //Vars: DeltaPhi, InvMass, PtTrack, Displacement, DeltaEta --> Last bin for pTassoc is to avoid the overflow!
  Int_t nBinsPhi[5] = {32,150,(int)(2*fPtAssocLimit+1),3,16}; 
  Double_t binMinPhi[5] = {-TMath::Pi()/2.,1.9864,0.,0.,-1.6};  //is the minimum for all the bins
  Double_t binMaxPhi[5] = {3.*TMath::Pi()/2.,2.5864,fPtAssocLimit+0.5,3.,1.6};  //is the maximum for all the bins

  //Vars: DeltaPhi, InvMass, DeltaEta
  Int_t nBinsMix[5] = {32,150,16,(int)(2*fPtAssocLimit+1),2};
  Double_t binMinMix[5] = {-TMath::Pi()/2.,1.9864,-1.6,0.,-0.5};  //is the minimum for all the bins
  Double_t binMaxMix[5] = {3.*TMath::Pi()/2.,2.5864,1.6,fPtAssocLimit+0.5,1.5};  //is the maximum for all the bins

  Int_t nPoolForHistos=1;
  if(!fMergePools) nPoolForHistos= fCutsTracks->GetNZvtxPoolBins()*fCutsTracks->GetNCentPoolBins(); //multeplicity of histos in case of correct pools treatment: sum(SE_i/ME_i)
 
  for(Int_t i=0;i<fNPtBinsCorr;i++) {

    //Modify n of bins with fast speed: in the "for" loop since bins can depend on pT (e.g. mass bin)
    //setting of mass bin is done at the end of the loop!
    if(fSpeed==kOneBinSB) { //these with fast speed, only 1 SBL and 1 SBR bins
      if(i>=9) {nBinsPhi[0] = 32; nBinsPhi[1] = 67; nBinsPhi[3] = 1; nBinsPhi[4] = 16;}
      else {nBinsPhi[0] = 32; nBinsPhi[1] = 43; nBinsPhi[3] = 1; nBinsPhi[4] = 16;}
      binMinPhi[0] = -TMath::Pi()/2.; binMinPhi[1] = 1.9864; binMinPhi[3] = 0.; binMinPhi[4] = -1.6;
      binMaxPhi[0] = 3.*TMath::Pi()/2.; binMaxPhi[1] = 2.5864; binMaxPhi[3] = 3.; binMaxPhi[4] = 1.6;
    
      if(i>=9) {nBinsMix[0] = 32; nBinsMix[1] = 67; nBinsMix[2] = 16;}
      else {nBinsMix[0] = 32; nBinsMix[1] = 43; nBinsMix[2] = 16;} 
      binMinMix[0] = -TMath::Pi()/2.; binMinMix[1] = 1.9864; binMinMix[2] = -1.6;
      binMaxMix[0] = 3.*TMath::Pi()/2.; binMaxMix[1] = 2.5864; binMaxMix[2] = 1.6;
    }  	  
    if(fSpeed==kOneBinSBandS) { //these with fast speed, only 1 SBL+SBR bin and 1 Sgin bin = total of 2 bins in mass axis!
      if(i>=9) {nBinsPhi[0] = 32; nBinsPhi[1] = 2; nBinsPhi[3] = 1; nBinsPhi[4] = 16;}
      else {nBinsPhi[0] = 32; nBinsPhi[1] = 2; nBinsPhi[3] = 1; nBinsPhi[4] = 16;}
      binMinPhi[0] = -TMath::Pi()/2.; binMinPhi[1] = 1.9864; binMinPhi[3] = 0.; binMinPhi[4] = -1.6;
      binMaxPhi[0] = 3.*TMath::Pi()/2.; binMaxPhi[1] = 2.5864; binMaxPhi[3] = 3.; binMaxPhi[4] = 1.6;
    
      if(i>=9) {nBinsMix[0] = 32; nBinsMix[1] = 2; nBinsMix[2] = 16;}
      else {nBinsMix[0] = 32; nBinsMix[1] = 2; nBinsMix[2] = 16;} 
      binMinMix[0] = -TMath::Pi()/2.; binMinMix[1] = 1.9864; binMinMix[2] = -1.6;
      binMaxMix[0] = 3.*TMath::Pi()/2.; binMaxMix[1] = 2.5864; binMaxMix[2] = 1.6;
    }   


    if(!fMixing) {
     if(!fFillTrees) {
      for(Int_t k=0; k<nPoolForHistos; k++) {    	    
    	    
        //THnSparse plots: correlations for various invariant mass (MC and data)
    
        namePlot="hPhi_Charg_Bin";
        namePlot+=i; namePlot+="_p"; namePlot+=k;

        THnSparseF *hPhiC = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiC->Sumw2();
        fOutputCorr->Add(hPhiC);
  
        //histos for c/b origin for Lambdac (MC only)
        if (fReadMC) {

          namePlot="hPhi_Charg_From_c_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiC_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_c->Sumw2();
          fOutputCorr->Add(hPhiC_c);
 

          namePlot="hPhi_Charg_From_b_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiC_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_b->Sumw2();
          fOutputCorr->Add(hPhiC_b);


          namePlot="hPhi_Charg_HF_From_c_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiC_HF_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_HF_c->Sumw2();
          fOutputCorr->Add(hPhiC_HF_c);


          namePlot="hPhi_Charg_HF_From_b_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;
     
          THnSparseF *hPhiC_HF_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_HF_b->Sumw2();
          fOutputCorr->Add(hPhiC_HF_b);

          namePlot="hPhi_Charg_NonHF_Bin";
          namePlot+=i; namePlot+="_p"; namePlot+=k;

          THnSparseF *hPhiC_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
          hPhiC_Non->Sumw2();
          fOutputCorr->Add(hPhiC_Non);
        } //end of MC
  
        //modify here the mass axis of THnSparse! 
        if(fSpeed==kOneBinSB) {
      	  Int_t nBins; Double_t mBin;      
      	  if(i>=9) { //signal range is xxx to xxx, plus 1 bin L and R for sidebands
      	    nBins = 67;
      	    mBin = 2.1704;
      	  }
          else { //signal range is xxx to xxx, plus 1 bin L and R for sidebands
            nBins = 43;
      	    mBin = 2.2224;
          }
     	
          Double_t varBins[nBins+1];
          varBins[0] = 1.9864;
          varBins[1] = 1.9904;
      	  for(int j = 2; j<nBins-1; j++) {varBins[j]=mBin; mBin+=0.004;}
      	  varBins[nBins-1] = 2.5824;
      	  varBins[nBins] = 2.5864;
        


      	  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
          if (fReadMC) {


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_HF_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_HF_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_NonHF_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
          }
        } //end of fSpeed==1

                //modify here the mass axis of THnSparse! 
        if(fSpeed==kOneBinSBandS) {
          Int_t nBins = 2;
          Double_t varBins[3];
          varBins[0] = -0.5;
          varBins[1] = 0.5;
          varBins[2] = 1.5;
        
       
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
          if (fReadMC) {


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_HF_From_c_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_HF_From_b_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);


            ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_NonHF_Bin%d_p%d",i,k)))->GetAxis(1)->Set(nBins, varBins);
          }
        } //end of fSpeed==2

      } //end of pool multiplicity      
   
      //Resume the definition of histos
      if(!fSpeed) {
        //leading hadron correlations
        namePlot="hPhi_Lead_Bin";
        namePlot+=i;

        THnSparseF *hCorrLead = new THnSparseF(namePlot.Data(), "Leading particle correlations; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
        hCorrLead->Sumw2();
        fOutputCorr->Add(hCorrLead);

        if (fReadMC) {
          namePlot="hPhi_Lead_From_c_Bin";
          namePlot+=i;

          THnSparseF *hCorrLead_c = new THnSparseF(namePlot.Data(), "Leading particle correlations - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_c->Sumw2();
          fOutputCorr->Add(hCorrLead_c);
  
          namePlot="hPhi_Lead_From_b_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrLead_b = new THnSparseF(namePlot.Data(), "Leading particle correlations - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_b->Sumw2();
          fOutputCorr->Add(hCorrLead_b);
  
          namePlot="hPhi_Lead_HF_From_c_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrLead_HF_c = new THnSparseF(namePlot.Data(), "Leading particle correlations HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_HF_c->Sumw2();
          fOutputCorr->Add(hCorrLead_HF_c);
  
          namePlot="hPhi_Lead_HF_From_b_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrLead_HF_b = new THnSparseF(namePlot.Data(), "Leading particle correlations HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_HF_b->Sumw2();
          fOutputCorr->Add(hCorrLead_HF_b);

          namePlot="hPhi_Lead_NonHF_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrLead_Non = new THnSparseF(namePlot.Data(), "Leading particle correlations - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          hCorrLead_Non->Sumw2();
          fOutputCorr->Add(hCorrLead_Non);
        }
      
        //pT weighted correlations
        namePlot="hPhi_Weig_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrWeig = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted); #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
        fOutputCorr->Add(hCorrWeig);
  
        if (fReadMC) {
          namePlot="hPhi_Weig_From_c_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_c = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_c);
  
          namePlot="hPhi_Weig_From_b_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_b = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_b);
  
          namePlot="hPhi_Weig_HF_From_c_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_HF_c = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_HF_c);
  
          namePlot="hPhi_Weig_HF_From_b_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_HF_b = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_HF_b);

          namePlot="hPhi_Weig_NonHF_Bin";
          namePlot+=i;
  
          THnSparseF *hCorrWeig_Non = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
          fOutputCorr->Add(hCorrWeig_Non);
        }
      } //end of fSpeed
     } //end of !fFillTrees

     //pT distribution histos
     namePlot = "hist_Pt_Charg_Bin"; namePlot+=i;
     TH1F *hPtC = new TH1F(namePlot.Data(), "Charged track pT (in Lambdac evs); p_{T} (GeV/c)",240,0.,12.);
     hPtC->SetMinimum(0);
     fOutputStudy->Add(hPtC);

     //Events multiplicity
     namePlot = "hMultiplEvt_Bin"; namePlot+=i;
     TH1F *hMultEv = new TH1F(namePlot.Data(), "Event multiplicity",1500,0.,6000.);
     hMultEv->SetMinimum(0);
     fOutputStudy->Add(hMultEv);

    } //end of !fMixing

    if(fMixing && !fFillTrees) {
      for(Int_t k=0; k<nPoolForHistos; k++) {     	    
        //THnSparse plots for event mixing!
  
        namePlot="hPhi_Charg_Bin";
        namePlot+=i; namePlot+="_p"; namePlot+=k; namePlot+="_EvMix";

        THnSparseF *hPhiC_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsMix,binMinMix,binMaxMix);
        hPhiC_EvMix->Sumw2();
        fOutputCorr->Add(hPhiC_EvMix);  

        //modify here the mass axis of THnSparse! 
        if(fSpeed==kOneBinSB) {
      	  Int_t nBins; Double_t mBin;      
      	  if(i>=9) { //signal range is xxx to xxx, plus 1 bin L and R for sidebands
      	    nBins = 67;
      	    mBin = 2.1704;
      	  }
          else { //signal range is xxx to xxx, plus 1 bin L and R for sidebands
            nBins = 43;
      	    mBin = 2.2224;
          }
     	
          Double_t varBins[nBins+1];
          varBins[0] = 1.9864;
          varBins[1] = 1.9904;
      	  for(int j = 2; j<nBins-1; j++) {varBins[j]=mBin; mBin+=0.004;}
      	  varBins[nBins-1] = 2.5824;
      	  varBins[nBins] = 2.5864;
        
      	 
      	  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d_EvMix",i,k)))->GetAxis(1)->Set(nBins, varBins);
      	  
        } //end of fSpeed==1

        if(fSpeed==kOneBinSBandS) {
          Int_t nBins = 2;
          Double_t varBins[3];          
          varBins[0] = -0.5;
          varBins[1] = 0.5;
          varBins[2] = 1.5;
        
        
          ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d_EvMix",i,k)))->GetAxis(1)->Set(nBins, varBins);
          
        } //end of fSpeed==2          
      
      } //end of Mult pools
    } //end of Mix
 
    //both for SE and for ME
    
    //Sigmac feeddown pions rejection histos
    namePlot = "hSigmacPionsVsDmass_Bin"; namePlot+=i;
    TH2F *hSigmacPions = new TH2F(namePlot.Data(), "Tracks rejected for Sigmac inv.mass cut vs Lc inv mass; # Tracks",2,0.,2.,150,1.9864,2.5864);
    hSigmacPions->GetXaxis()->SetBinLabel(1,"Not rejected");
    hSigmacPions->GetXaxis()->SetBinLabel(2,"Rejected");
    hSigmacPions->SetMinimum(0);
    fOutputStudy->Add(hSigmacPions); 

    namePlot = "hSigmacPionsVsdeltaPhi_Bin"; namePlot+=i;
    TH2F *hSigmacPions2 = new TH2F(namePlot.Data(), "Tracks rejected for Sigmac inv.mass cut vs deltaPhi; # Tracks",2,0.,2.,64,-TMath::Pi()/2.,3.*TMath::Pi()/2.);
    hSigmacPions2->GetXaxis()->SetBinLabel(1,"Not rejected");
    hSigmacPions2->GetXaxis()->SetBinLabel(2,"Rejected");
    hSigmacPions2->SetMinimum(0);
    fOutputStudy->Add(hSigmacPions2); 
    
    if(!fFillTrees) {
      //ME filling control plots
      namePlot="hEvtsPerPool_"; namePlot+=i;
      TH1F *hEvPerPool = new TH1F(namePlot.Data(), "Events With selLambdac in ME pools",nPoolForHistos,0.,nPoolForHistos);
      hEvPerPool->SetMinimum(0);
      fOutputStudy->Add(hEvPerPool);
    }

  } //end of bin loop

  //out of bin loop
  TH1F *hCountC = new TH1F("hist_Count_Charg", "Charged track counter; # Tracks",6000,0.,6000.);
  hCountC->SetMinimum(0);
  fOutputStudy->Add(hCountC);

  TH1F *hZvtx = new TH1F("hZvtx", "z of Primary vtx (for events with selected D); z (cm); # Events",48,-12.,12.);
  hZvtx->SetMinimum(0);
  fOutputStudy->Add(hZvtx);
  
  TH1F *hZvtxEvts = new TH1F("hZvtxEvts", "z of Primary vtx (for selected events); z (cm); # Events",120,-30.,30.);
  hZvtxEvts->SetMinimum(0);
  fOutputStudy->Add(hZvtxEvts);
  
  TH1F *hCentralEvts = new TH1F("hCentralEvts","Centrality of events (std approach); centrality; # Events",10002,-0.01,100.01);
  hCentralEvts->SetMinimum(0);
  fOutputStudy->Add(hCentralEvts);  

  TH1F *hCentEvV0M = new TH1F("hCentEvV0M","Centrality of events (v2 pp analysis) in V0M percentiles; V0M perc.; # Events",10002,-0.01,100.01);
  hCentEvV0M->SetMinimum(0);
  fOutputStudy->Add(hCentEvV0M);  

  TH1F *hMultEvV0M = new TH1F("hMultEvV0M","Multiplicity of events (v2 pp analysis) in V0M amplitude; V0M amplitude; # Events",1000,0,1000);
  hMultEvV0M->SetMinimum(0);
  fOutputStudy->Add(hMultEvV0M);  

  TH1F *hMultEvV0MEqual = new TH1F("hMultEvV0MEqual","Multiplicity of events (v2 pp analysis) in V0M equalized amplitude; V0M amplitude (eq.); # Events",1000,0,1000);
  hMultEvV0MEqual->SetMinimum(0);
  fOutputStudy->Add(hMultEvV0MEqual);  

  TH1F *hMultEvTrkl1 = new TH1F("hMultEvTrkl1","Multiplicity of events (v2 pp analysis) in Tracklets <1; SPD tracklets in |eta|<1; # Events",200,0,200);
  hMultEvTrkl1->SetMinimum(0);
  fOutputStudy->Add(hMultEvTrkl1);   

  TH1F *hNtrUnCorrEvSel = new TH1F("hNtrUnCorrEvSel","Uncorrected Trkl multiplicity for selected events; Trkl ; Entries",200,-0.5,199.5);
  hNtrUnCorrEvSel->SetMinimum(0);
  hNtrUnCorrEvSel->Sumw2();
  fOutputStudy->Add(hNtrUnCorrEvSel); 

  TH1F *hNtrUnCorrEvWithCand = new TH1F("hNtrUnCorrEvWithCand","Uncorrected Trkl multiplicity for events with D candidates; Trkl ; Entries",200,-0.5,199.5);// Total multiplicity
  hNtrUnCorrEvWithCand->SetMinimum(0);
  hNtrUnCorrEvWithCand->Sumw2();
  fOutputStudy->Add(hNtrUnCorrEvWithCand); 

  if(fEqualizeTracklets) {
    TH1F *hMultEvTrkl1Equal = new TH1F("hMultEvTrkl1Equal","Multiplicity of events (v2 pp analysis) in Tracklets <1, EQUALIZED; SPD tracklets in |eta|<1; # Events",200,0,200);
    hMultEvTrkl1Equal->SetMinimum(0);
    fOutputStudy->Add(hMultEvTrkl1Equal);  

    TH2F *hNtrVsZvtx = new TH2F("hNtrVsZvtx","Ntracklets vs VtxZ; VtxZ;N_{trkl};",300,-15,15,150,-0.5,149.5); //
    hNtrVsZvtx->SetMinimum(0);
    fOutputStudy->Add(hNtrVsZvtx);       

    TH2F *hNtrCorrVsZvtx = new TH2F("hNtrCorrVsZvtx","Ntracklets (corrected) vs VtxZ; VtxZ;N_{trkl};",300,-15,15,150,-0.5,149.5); //
    hNtrCorrVsZvtx->SetMinimum(0);
    fOutputStudy->Add(hNtrCorrVsZvtx); 

    TH1F *hNtrCorrEvSel = new TH1F("hNtrCorrEvSel","Corrected Trkl multiplicity for selected events; Trkl ; Entries",200,-0.5,199.5);
    hNtrCorrEvSel->SetMinimum(0);
    hNtrCorrEvSel->Sumw2();
    fOutputStudy->Add(hNtrCorrEvSel); 

    TH1F *hNtrCorrEvWithCand = new TH1F("hNtrCorrEvWithCand","Corrected Trkl multiplicity for events with D candidates; Trkl ; Entries",200,-0.5,199.5);// Total multiplicity
    hNtrCorrEvWithCand->SetMinimum(0);
    hNtrCorrEvWithCand->Sumw2();
    fOutputStudy->Add(hNtrCorrEvWithCand); 
  }

  if(fVsMultAnalysis) {
    TH1F *hCentEvV0MSelEvents = new TH1F("hCentEvV0MSelEvents","Centrality of events (for selected events, v2 pp analysis) in V0M percentiles; V0M perc.; # Events",10002,-0.01,100.01);
    hCentEvV0MSelEvents->SetMinimum(0);
    fOutputStudy->Add(hCentEvV0MSelEvents);  
    
    TH1F *hMultEvV0MSelEvents = new TH1F("hMultEvV0MSelEvents","Multiplicity of events (for selected events, v2 pp analysis) in V0M amplitude; V0M amplitude; # Events",1000,0,1000);
    hMultEvV0MSelEvents->SetMinimum(0);
    fOutputStudy->Add(hMultEvV0MSelEvents);  

    TH1F *hMultEvV0MEqualSelEvents = new TH1F("hMultEvV0MEqualSelEvents","Multiplicity of events (for selected events, v2 pp analysis) in V0M equalized amplitude; V0M amplitude (eq.); # Events",1000,0,1000);
    hMultEvV0MEqualSelEvents->SetMinimum(0);
    fOutputStudy->Add(hMultEvV0MEqualSelEvents);  

    TH1F *hMultEvTrkl1SelEvents = new TH1F("hMultEvTrkl1SelEvents","Multiplicity of events (for selected events, v2 pp analysis) in Tracklets <1; SPD tracklets in |eta|<1; # Events",200,0,200);
    hMultEvTrkl1SelEvents->SetMinimum(0);
    fOutputStudy->Add(hMultEvTrkl1SelEvents);       

    if(fEqualizeTracklets) {
      TH1F *hMultEvTrkl1EqualSelEvents = new TH1F("hMultEvTrkl1EqualSelEvents","Multiplicity of events (for selected events, v2 pp analysis) in Tracklets <1, EQUALIZED; SPD tracklets in |eta|<1; # Events",200,0,200);
      hMultEvTrkl1EqualSelEvents->SetMinimum(0);
      fOutputStudy->Add(hMultEvTrkl1EqualSelEvents);       

      TH2F *hNtrVsZvtxSelEvents = new TH2F("hNtrVsZvtxSelEvents","Ntracklets vs VtxZ; VtxZ;N_{trkl};",300,-15,15,150,-0.5,149.5); //
      hNtrVsZvtxSelEvents->SetMinimum(0);
      fOutputStudy->Add(hNtrVsZvtxSelEvents);       

      TH2F *hNtrCorrVsZvtxSelEvents = new TH2F("hNtrCorrVsZvtxSelEvents","Ntracklets (corrected) vs VtxZ; VtxZ;N_{trkl};",300,-15,15,150,-0.5,149.5); //
      hNtrCorrVsZvtxSelEvents->SetMinimum(0);
      fOutputStudy->Add(hNtrCorrVsZvtxSelEvents);       
    }
  }

  TH1F *hZeroEff = new TH1F("hZeroEff","Efficiency debug plot (0-eff cases); # Events",4,0.,4.);
  hZeroEff->GetXaxis()->SetBinLabel(1,"Lambdac ok eff");
  hZeroEff->GetXaxis()->SetBinLabel(2,"Lambdac 0 eff");
  hZeroEff->GetXaxis()->SetBinLabel(3,"Track ok eff");
  hZeroEff->GetXaxis()->SetBinLabel(4,"Track 0 eff");
  hZeroEff->SetMinimum(0);
  fOutputStudy->Add(hZeroEff);  

  if (fReadMC) {
    TH1D *hEventTypeMC = new TH1D("EventTypeMC","EventTypeMC",100,-0.5,99.5);
    fOutputStudy->Add(hEventTypeMC); 
  }

  if (fFillGlobal) { //all-events plots
    //pt distributions
    TH1F *hPtCAll = new TH1F("hist_Pt_Charg_AllEv", "Charged track pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtCAll->SetMinimum(0);
    fOutputStudy->Add(hPtCAll);

    TH1F *hPtHAll = new TH1F("hist_Pt_Kcharg_AllEv", "Kaons pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtHAll->SetMinimum(0);
    fOutputStudy->Add(hPtHAll);

    TH1F *hPtKAll = new TH1F("hist_Pt_K0_AllEv", "K0 pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtKAll->SetMinimum(0);
    fOutputStudy->Add(hPtKAll);

    //K0 Invariant Mass plots
    TH2F *hK0MassInv = new TH2F("hK0MassInv", "K0 invariant mass; Invariant mass (MeV/c^{2}); pT (GeV/c)",200,0.4,0.6,100,0.,10.);
    hK0MassInv->SetMinimum(0);
    fOutputStudy->Add(hK0MassInv);
  }

  if(!fMixing) {
    //phi distributions
    TH1F *hPhiDistCAll = new TH1F("hist_PhiDistr_Charg", "Charged track phi distr. (All); #varphi (rad)",64,0,6.283);
    hPhiDistCAll->SetMinimum(0);
    fOutputStudy->Add(hPhiDistCAll);

    TH1F *hPhiDistHAll = new TH1F("hist_PhiDistr_Kcharg", "Kaons phi distr. (All); #varphi (rad)",64,0,6.283);
    hPhiDistHAll->SetMinimum(0);
    fOutputStudy->Add(hPhiDistHAll);

    TH1F *hPhiDistKAll = new TH1F("hist_PhiDistr_K0", "K0 phi distr. (All); #varphi (rad)",64,0,6.283);
    hPhiDistKAll->SetMinimum(0);
    fOutputStudy->Add(hPhiDistKAll);

    TH1F *hPhiDistDAll = new TH1F("hist_PhiDistr_Lambdac", "#Lambda_{c} phi distr. (All); #varphi (rad)",64,0,6.283);
    hPhiDistDAll->SetMinimum(0);
    fOutputStudy->Add(hPhiDistDAll);

    //eta distributions
    TH1F *hEtaDistCAll = new TH1F("hist_EtaDistr_Charg", "Charged track eta distr. (All); #eta (rad)",40,-1,1);
    hEtaDistCAll->SetMinimum(0);
    fOutputStudy->Add(hEtaDistCAll);

    TH1F *hEtaDistHAll = new TH1F("hist_EtaDistr_Kcharg", "Kaons eta distr. (All); #eta (rad)",40,-1,1);
    hEtaDistHAll->SetMinimum(0);
    fOutputStudy->Add(hEtaDistHAll);

    TH1F *hEtaDistKAll = new TH1F("hist_EtaDistr_K0", "K0 eta distr. (All); #eta (rad)",40,-1,1);
    hEtaDistKAll->SetMinimum(0);
    fOutputStudy->Add(hEtaDistKAll);

    TH1F *hEtaDistDAll = new TH1F("hist_EtaDistr_Lambdac", "#Lambda_{c} eta distr. (All); #eta (rad)",40,-1,1);
    hEtaDistDAll->SetMinimum(0);
    fOutputStudy->Add(hEtaDistDAll);
    
    //phivsEta
    TH2F *hPhiVsEtaDistCAll = new TH2F("hist_PhiVsEtaDistr_Charg", "Phi vs Eta distribution - Charged tracks",64,0,6.283,40,-1,1);
    hPhiVsEtaDistCAll->SetMinimum(0);
    fOutputStudy->Add(hPhiVsEtaDistCAll);     
    
    TH2F *hPhiVsEtaDistHAll = new TH2F("hist_PhiVsEtaDistr_Kcharg", "Phi vs Eta distribution - Charged Kaons",64,0,6.283,40,-1,1);
    hPhiVsEtaDistHAll->SetMinimum(0);
    fOutputStudy->Add(hPhiVsEtaDistHAll);  
    
    TH2F *hPhiVsEtaDistKAll = new TH2F("hist_PhiVsEtaDistr_K0", "Phi vs Eta distribution - K^{0}",64,0,6.283,40,-1,1);
    hPhiVsEtaDistKAll->SetMinimum(0);
    fOutputStudy->Add(hPhiVsEtaDistKAll);  
    
    TH2F *hPhiVsEtaDistDAll = new TH2F("hist_PhiVsEtaDistr_Lambdac", "Phi vs Eta distribution - #Lambda_{c}",64,0,6.283,40,-1,1);
    hPhiVsEtaDistDAll->SetMinimum(0);
    fOutputStudy->Add(hPhiVsEtaDistDAll);  
    }

  //for MC analysis only
  for(Int_t i=0;i<fNPtBinsCorr;i++) {

    if (fReadMC && !fMixing) {

      namePlot="histDispl_Charg_Bin"; namePlot+=i;
      TH1F *hDisplCharg = new TH1F(namePlot.Data(), "Charged tracks Displacement; DCA",150,0.,0.15);
      hDisplCharg->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg);
  
      namePlot="histDispl_Charg_HF_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplCharg_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF);

      namePlot="histDispl_Charg_From_c_Bin"; namePlot+=i;
      TH1F *hDisplCharg_c = new TH1F(namePlot.Data(), "Charged tracks Displacement - c origin; DCA",150,0.,0.15);
      hDisplCharg_c->Sumw2();
      hDisplCharg_c->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_c);
  
      namePlot="histDispl_Charg_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF_c = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplCharg_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF_c);

      namePlot="histDispl_Charg_From_b_Bin"; namePlot+=i;
      TH1F *hDisplCharg_b = new TH1F(namePlot.Data(), "Charged tracks Displacement - b origin; DCA",150,0.,0.15);
      hDisplCharg_b->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_b);
  
      namePlot="histDispl_Charg_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF_b = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplCharg_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF_b);

      //origin of tracks histos
      namePlot="histOrig_Charg_Bin";  namePlot+=i;
      TH1F *hOrigin_Charm = new TH1F(namePlot.Data(), "Origin of charged tracks",9,0.,9.);
      hOrigin_Charm->SetMinimum(0);
      hOrigin_Charm->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_Charm->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(4,"c hadr.");
      hOrigin_Charm->GetXaxis()->SetBinLabel(5,"B->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(6,"B->X-># (X!=D)");
      hOrigin_Charm->GetXaxis()->SetBinLabel(7,"B->D->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(8,"B->D->X->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_Charm);
    }

    if (fReadMC) {
      //origin of Lambdac histos
      namePlot="histOrig_Lambdac_Bin";  namePlot+=i;
      TH1F *hOrigin_Lambdac = new TH1F(namePlot.Data(), "Origin of Lambdac",2,0.,2.);
      hOrigin_Lambdac->SetMinimum(0);
      hOrigin_Lambdac->GetXaxis()->SetBinLabel(1,"From c");
      hOrigin_Lambdac->GetXaxis()->SetBinLabel(2,"From b");
      fOutputStudy->Add(hOrigin_Lambdac);

      //primary tracks (Kine & Reco)
      namePlot="hPhysPrim_Bin";  namePlot+=i;
      TH1F *hPhysPrim = new TH1F(namePlot.Data(), "Origin of hadrons",2,0.,2.);
      hPhysPrim->SetMinimum(0);
      hPhysPrim->GetXaxis()->SetBinLabel(1,"OK");
      hPhysPrim->GetXaxis()->SetBinLabel(2,"NO");
      fOutputStudy->Add(hPhysPrim);
    }
  } //end of for on pTbins

  if(fPurityStudies) {
    
    TString namebinD[6] = {"2to3","3to5","5to8","8to16","16to24","24to36"};
    TString namebinAss[7] = {"03to99","03to1","1to99","1to3","1to2","2to3","3to99"};

    for(int i=0; i<6; i++) { //pTD
      for(int j=0; j<7; j++) { //pTass
	namePlot=Form("hPurityCount_PrimAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data());
        TH1F *hpurity_prim = new TH1F(namePlot.Data(), "Prim accepted",1,-0.5,0.5);
        hpurity_prim->SetMinimum(0);
        hpurity_prim->Sumw2();
        hpurity_prim->GetXaxis()->SetBinLabel(1,"Accepted");
        fOutputStudy->Add(hpurity_prim);

	namePlot=Form("hPurityCount_SecAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data());
        TH1F *hpurity_sec = new TH1F(namePlot.Data(), "Sec accepted",1,-0.5,0.5);
        hpurity_sec->SetMinimum(0);
        hpurity_sec->Sumw2();
        hpurity_sec->GetXaxis()->SetBinLabel(1,"Accepted");
        fOutputStudy->Add(hpurity_sec);

	namePlot=Form("hPurityCount_CharmAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data());
        TH1F *hpurity_c = new TH1F(namePlot.Data(), "Charm accepted",1,-0.5,0.5);
        hpurity_c->SetMinimum(0);
        hpurity_c->Sumw2();
        hpurity_c->GetXaxis()->SetBinLabel(1,"Accepted");
        fOutputStudy->Add(hpurity_c);

	namePlot=Form("hPurityCount_BeautyAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data());
        TH1F *hpurity_b = new TH1F(namePlot.Data(), "Beauty accepted",1,-0.5,0.5);
        hpurity_b->SetMinimum(0);
        hpurity_b->Sumw2();
        hpurity_b->GetXaxis()->SetBinLabel(1,"Accepted");
        fOutputStudy->Add(hpurity_b);

	namePlot=Form("hPuritydPhi_PrimAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data());
        TH1F *hpuritydphi_prim = new TH1F(namePlot.Data(), "Prim accepted vs dPhi",32,-TMath::Pi()/2.,3*TMath::Pi()/2.);
        hpuritydphi_prim->SetMinimum(0);
        hpuritydphi_prim->Sumw2();
        fOutputStudy->Add(hpuritydphi_prim);

	namePlot=Form("hPuritydPhi_SecAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data());
        TH1F *hpuritydphi_sec = new TH1F(namePlot.Data(), "Sec accepted vs dPhi",32,-TMath::Pi()/2.,3*TMath::Pi()/2.);
        hpuritydphi_sec->SetMinimum(0);
        hpuritydphi_sec->Sumw2();
        fOutputStudy->Add(hpuritydphi_sec);

	namePlot=Form("hPuritydPhi_CharmAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data());
        TH1F *hpuritydphi_c = new TH1F(namePlot.Data(), "Charm accepted vs dPhi",32,-TMath::Pi()/2.,3*TMath::Pi()/2.);
        hpuritydphi_c->SetMinimum(0);
        hpuritydphi_c->Sumw2();
        fOutputStudy->Add(hpuritydphi_c);

	namePlot=Form("hPuritydPhi_BeautyAccepted_pTD%s_pTass%s",namebinD[i].Data(),namebinAss[j].Data());
        TH1F *hpuritydphi_b = new TH1F(namePlot.Data(), "Beauty accepted vs dPhi",32,-TMath::Pi()/2.,3*TMath::Pi()/2.);
        hpuritydphi_b->SetMinimum(0);
        hpuritydphi_b->Sumw2();
        fOutputStudy->Add(hpuritydphi_b);
      }
    }

  } //end of purity studies

}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::CalculateCorrelations(AliAODRecoCascadeHF* d, Int_t labLambdac, TClonesArray* mcArray) {
//
// Method for correlations Lambdac-hadrons study
//
  Int_t N_Charg = 0, N_KCharg = 0, N_Kaons = 0;
//  Double_t mLambdac, mLambdacbar;
  Double_t mLambdac;//, mLambdacbar;
  Int_t origLambdac = 0, PDGLambdac = 0, ptbin = 0;
  //mLambdac = d->InvMassLctoK0sP();
  mLambdac = fKFLcInvMass;
  //mLambdacbar = d->InvMassLctoK0sP();
  //Double_t mInv[2] = {mLambdac, mLambdacbar};
  Double_t mInv[] = {mLambdac};
  ptbin = PtBinCorr(fKFLcPt);

  if(ptbin < 0) return;

  //Fill of Lambdac phi distribution
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_PhiDistr_Lambdac"))->Fill(fKFLcPhi);  
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_EtaDistr_Lambdac"))->Fill(fKFLcEta); 
  if (!fMixing) ((TH2F*)fOutputStudy->FindObject("hist_PhiVsEtaDistr_Lambdac"))->Fill(fKFLcPhi,fKFLcEta); 

  //Origin of Lambdac
  TString orig="";
  if(fReadMC) {
    origLambdac=CheckLambdacOrigin(mcArray,(AliAODMCParticle*)mcArray->At(labLambdac));
    PDGLambdac = ((AliAODMCParticle*)mcArray->At(labLambdac))->GetPdgCode();
    switch (CheckLambdacOrigin(mcArray,(AliAODMCParticle*)mcArray->At(labLambdac))) {
      case 4:
        orig = "_From_c";
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_Lambdac_Bin%d",ptbin)))->Fill(0.);
        break;
      case 5:
        orig = "_From_b";
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_Lambdac_Bin%d",ptbin)))->Fill(1.);
        break;
      default:
        return;
    }
  }

  Double_t highPt = 0; Double_t lead[4] = {0,0,0,1};  //infos for leading particle (pt,deltaphi)

  //loop over the tracks in the pool 
  Bool_t execPoolTr = fCorrelatorTr->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
 
		
  Int_t NofEventsinPool = 1;
  if(fMixing) {
    NofEventsinPool = fCorrelatorTr->GetNofEventsInPool();
    if(!execPoolTr) {
      AliInfo("Mixed event analysis: track pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //Charged tracks
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksTr = fCorrelatorTr->ProcessAssociatedTracks(jMix);// process all the tracks in the aodEvent, by applying the selection cuts
    if(!analyzetracksTr) {
      AliInfo("AliHFCorrelator::Cannot process the track array");
      continue;
    }
	
    for(Int_t iTrack = 0; iTrack<fCorrelatorTr->GetNofTracks(); iTrack++){ // looping on track candidates

      Bool_t runcorrelation = fCorrelatorTr->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* track = fCorrelatorTr->GetAssociatedParticle();

      if(!fMixing) {
        Int_t idDaughs[3] = {((AliVTrack*)d->GetBachelor())->GetID(),((AliVTrack*)d->Getv0PositiveTrack())->GetID(),((AliVTrack*)d->Getv0NegativeTrack())->GetID()}; //IDs of daughters to be skipped
        if(track->GetID() == idDaughs[0] || track->GetID() == idDaughs[1] || track->GetID() == idDaughs[2]) continue; //discards daughters of candidate
      }
      if(track->Pt() < fPtThreshLow.at(ptbin) || track->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K

      if(fReadMC) {
        AliAODMCParticle* trkKine = (AliAODMCParticle*)mcArray->At(track->GetLabel());
        if (!trkKine) continue;
        //remove secondary tracks (for MC closure test, but obviously not for purity studies)
        if (!trkKine->IsPhysicalPrimary()) {
 	        ((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(1.);  
  	      if(!fPurityStudies) continue; //reject the Reco track if correspondent Kine track is not primary
        } else ((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(0.);
        //remove tracks not being pi/K/p/e/mu (for MC closure studies, but not for purity studies
        // --> the sense is that in the purity correction we also remove primary reco particles not being pi/K/p/e/mu
        // --> (we don't want them in the MC closure instead, since there we don't apply purity)
        if(!fPurityStudies) {
          Int_t pdg = TMath::Abs(trkKine->GetPdgCode());
          if(!((pdg==321)||(pdg==211)||(pdg==2212)||(pdg==13)||(pdg==11))) continue;
        }
      }

      Double_t effTr = track->GetWeight(); //extract track efficiency
      Double_t effLambdac = 1.;
      if(fReadMC) {
        if(origLambdac==4) effLambdac = fCutsTracks->GetTrigWeight(fKFLcPt,fMultEv);
        if(origLambdac==5) effLambdac = fCutsTracks->GetTrigWeightB(fKFLcPt,fMultEv);
      } else effLambdac = fCutsTracks->GetTrigWeight(fKFLcPt,fMultEv);
      if(!fUseLceff) effLambdac=1.; 
      if(!fUseTrackeff) effTr=1.; 
      Double_t eff = effTr*effLambdac;
      if(!eff) eff = 1; //safety check

      //debug histogram
      if(!effLambdac) ((TH1F*)fOutputStudy->FindObject("hZeroEff"))->Fill(1.5); 
        else ((TH1F*)fOutputStudy->FindObject("hZeroEff"))->Fill(0.5);
      if(!effTr) ((TH1F*)fOutputStudy->FindObject("hZeroEff"))->Fill(3.5); 
        else ((TH1F*)fOutputStudy->FindObject("hZeroEff"))->Fill(2.5);

	    if(!fMixing) {
        if(fSoftPiCut && !track->CheckSoftPi()) { //removal of soft pions
          /* if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mLambdac);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mLambdacbar);
           */
          ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mLambdac);
          ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(1.,fCorrelatorTr->GetDeltaPhi());
    	  continue; //in SE events, just reject the soft pion
        } else { //not a soft pion
          /* if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mLambdac);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mLambdacbar); */
          ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mLambdac);
          ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(0.,fCorrelatorTr->GetDeltaPhi());
        }
      }
      if(fMixing) { 
        if(fSoftPiCut && !fCutsTracks->InvMassSigmacRejection(d,track,fKFLcInvMass)) { //removal of soft pions
          /* if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mLambdac);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mLambdacbar); */
          ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(1.,mLambdac);
          ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(1.,fCorrelatorTr->GetDeltaPhi());
          if(fMixing) FillSparsePlots(mcArray,mInv,origLambdac,PDGLambdac,track,ptbin,kTrack,1,1./eff); //in ME events, fill the THnSparse under the softpi hypothesis
    	  continue; 
        } else { //not a soft pion
          /* if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mLambdac);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mLambdacbar); */
          ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsDmass_Bin%d",ptbin)))->Fill(0.,mLambdac);
          ((TH2F*)fOutputStudy->FindObject(Form("hSigmacPionsVsdeltaPhi_Bin%d",ptbin)))->Fill(0.,fCorrelatorTr->GetDeltaPhi());
        }
      } 


      FillSparsePlots(mcArray,mInv,origLambdac,PDGLambdac,track,ptbin,kTrack,0,1./eff); //fills for charged tracks

      if(!fMixing) N_Charg++;

      //retrieving leading info...
      if(track->Pt() > highPt) {
        if(fReadMC && track->GetLabel()<1) continue;
        if(fReadMC && !(AliAODMCParticle*)mcArray->At(track->GetLabel())) continue;
        lead[0] = fCorrelatorTr->GetDeltaPhi();
        lead[1] = fCorrelatorTr->GetDeltaEta();
        lead[2] = fReadMC ? CheckTrackOrigin(mcArray,(AliAODMCParticle*)mcArray->At(track->GetLabel())) : 0;
        if(fReadMC) {
  	      if(origLambdac==4) lead[3] = 1./(track->GetWeight()*fCutsTracks->GetTrigWeight(fKFLcPt,fMultEv)); //weight is 1./efficiency
	        if(origLambdac==5) lead[3] = 1./(track->GetWeight()*fCutsTracks->GetTrigWeightB(fKFLcPt,fMultEv)); //weight is 1./efficiency
	      } else lead[3] = 1./(track->GetWeight()*fCutsTracks->GetTrigWeight(fKFLcPt,fMultEv));
        highPt = track->Pt();
      }

    } // end of tracks loop
  } //end of event loop for fCorrelatorTr

  Double_t fillSpLeadLambdac[4] = {lead[0],mLambdac,lead[1],0.4}; //dummy value for threshold of leading!
  //Double_t fillSpLeadLambdacbar[4] = {lead[0],mLambdacbar,lead[1],0.4};

  //leading track correlations fill
  if(!fMixing && !fSpeed) {
    if(fReadMC) {
      /* if(((AliAODMCParticle*)mcArray->At(labLambdac))->GetPdgCode()==4122  && (fIsSelectedCandidate==1||fIsSelectedCandidate==3)) { //Lambdac
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadLambdac,lead[3]); //c and b Lambdac
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdac,lead[3]); //c or b Lambdac
        if(origLambdac==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdac,lead[3]);  
        if(origLambdac==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdac,lead[3]);  
        if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadLambdac,lead[3]);  //non HF  
      }
      if(((AliAODMCParticle*)mcArray->At(labLambdac))->GetPdgCode()==-4122 && fIsSelectedCandidate>1 ) { //Lambdacbar
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadLambdacbar,lead[3]);
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdacbar,lead[3]); //c or b Lambdac
        if(origLambdac==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdacbar,lead[3]);  
        if(origLambdac==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdacbar,lead[3]); 
        if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadLambdacbar,lead[3]);  //non HF  
      } */
      if(std::abs(((AliAODMCParticle*)mcArray->At(labLambdac))->GetPdgCode())==4122  && (fIsSelectedCandidate>0)) { //Lambdac
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadLambdac,lead[3]); //c and b Lambdac
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdac,lead[3]); //c or b Lambdac
        if(origLambdac==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdac,lead[3]);  
        if(origLambdac==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadLambdac,lead[3]);  
        if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadLambdac,lead[3]);  //non HF  
      }
    } else {
        /* if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadLambdac,lead[3]); 
        if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadLambdacbar,lead[3]); */
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadLambdac,lead[3]); 
    }
  }
    //Fill of count histograms
  if (!fAlreadyFilled && !fMixing) { 
    ((TH1F*)fOutputStudy->FindObject("hist_Count_Charg"))->Fill(N_Charg);
  }


  fAlreadyFilled=kTRUE; //at least a Lambdac analyzed in the event; distribution plots already filled

}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::CalculateCorrelationsMCKine(AliAODMCParticle* d, TClonesArray* mcArray) {
//
// Method for correlations Lambdac-hadrons study
//
  Int_t N_Charg = 0, N_KCharg = 0, N_Kaons = 0;
  Double_t mLambdac = 2.286, mLambdacbar = 2.286;
  Double_t mInv[2] = {mLambdac, mLambdacbar};
  Int_t origLambdac = 0, PDGLambdac = 0;
  Int_t ptbin = PtBinCorr(d->Pt());

  if(ptbin < 0) return;

  //Fill of Lambdac phi distribution
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_PhiDistr_Lambdac"))->Fill(d->Phi()); 
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_EtaDistr_Lambdac"))->Fill(d->Phi()); 
  if (!fMixing) ((TH2F*)fOutputStudy->FindObject("hist_PhiVsEtaDistr_Lambdac"))->Fill(d->Phi(),d->Eta()); 
  
  //Origin of Lambdac
  TString orig="";
  origLambdac=CheckLambdacOrigin(mcArray,d);
  PDGLambdac = d->GetPdgCode();
  switch (CheckLambdacOrigin(mcArray,d)) {
    case 4:
      orig = "_From_c";
      ((TH1F*)fOutputStudy->FindObject(Form("histOrig_Lambdac_Bin%d",ptbin)))->Fill(0.);
      break;
    case 5:
      orig = "_From_b";
      ((TH1F*)fOutputStudy->FindObject(Form("histOrig_Lambdac_Bin%d",ptbin)))->Fill(1.);
      break;
    default:
      return;
  }

  Double_t highPt = 0; Double_t lead[3] = {0,0,0};  //infos for leading particle (pt,deltaphi)

  //loop over the tracks in the pool 
  Bool_t execPoolTr = fCorrelatorTr->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
 // Bool_t execPoolKc = fCorrelatorKc->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
//  Bool_t execPoolK0 = fCorrelatorK0->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
		
  Int_t NofEventsinPool = 1;
  if(fMixing) {
    NofEventsinPool = fCorrelatorTr->GetNofEventsInPool(); 
    if(!execPoolTr) {
      AliInfo("Mixed event analysis: track pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //Charged tracks
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)

    Bool_t analyzetracksTr = fCorrelatorTr->ProcessAssociatedTracks(jMix);// process all the tracks in the aodEvent, by applying the selection cuts
    if(!analyzetracksTr) {
      AliInfo("AliHFCorrelator::Cannot process the track array");
      continue;
    }
	
    for(Int_t iTrack = 0; iTrack<fCorrelatorTr->GetNofTracks(); iTrack++){ // looping on track candidates

      Bool_t runcorrelation = fCorrelatorTr->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* track = fCorrelatorTr->GetAssociatedParticle();
      if(track->GetLabel()<0) continue;
      if(track->Pt() < fPtThreshLow.at(ptbin) || track->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
      if(track->Pt() < 0.3 || TMath::Abs(track->Eta())>0.8) continue; //discard tracks outside barrel (since it's kinematic MC and produces tracks all over rapidity region
      if(!fMixing) N_Charg++;

      AliAODMCParticle *trkMC = (AliAODMCParticle*)mcArray->At(track->GetLabel());
      if(!trkMC) continue;

      if (!trkMC->IsPhysicalPrimary()) {  //reject material budget, or other fake tracks
 	((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(1.);  
  	continue;
      } else ((TH1F*)fOutputStudy->FindObject(Form("hPhysPrim_Bin%d",ptbin)))->Fill(0.);

      if (IsDDaughter(d,trkMC,mcArray)) continue;
      //if (fSoftPiCut && IsSoftPion_MCKine(d,trkMC,mcArray)) continue; //remove soft pions (if requestes, e.g. for templates)

      FillSparsePlots(mcArray,mInv,origLambdac,PDGLambdac,track,ptbin,kTrack,0); //fills for charged tracks

      //retrieving leading info...
      if(track->Pt() > highPt) {
        lead[0] = fCorrelatorTr->GetDeltaPhi();
        lead[1] = fCorrelatorTr->GetDeltaEta();
        lead[2] = fReadMC ? CheckTrackOrigin(mcArray,trkMC) : 0;
        highPt = track->Pt();
      }

    } // end of tracks loop
  } //end of event loop for fCorrelatorTr



  Double_t fillSpLeadMC[4] = {lead[0],mLambdac,lead[1],0.4}; //mLambdac = mLambdacbar = 2.286

  //leading track correlations fill
  if(!fMixing && !fSpeed) {
    /* if(d->GetPdgCode()==4122 && (fIsSelectedCandidate==1||fIsSelectedCandidate==3)) { //Lambdac
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadMC); //c and b Lambdac
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); //c or b Lambdac
      if(origLambdac==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if(origLambdac==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadMC);  //non HF
    }
    if(d->GetPdgCode()==-4122 && fIsSelectedCandidate>1) { //Lambdacbar
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadMC);
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); //c or b Lambdac
      if(origLambdac==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if(origLambdac==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); 
      if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadMC);  //non HF
    } */
    if(d->GetPdgCode()==4122 && (fIsSelectedCandidate>0)) { //Lambdac
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadMC); //c and b Lambdac
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); //c or b Lambdac
      if(origLambdac==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if(origLambdac==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadMC);  //non HF
    }
    if(d->GetPdgCode()==-4122 && fIsSelectedCandidate>0) { //Lambdacbar
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadMC);
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); //c or b Lambdac
      if(origLambdac==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if(origLambdac==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); 
      if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadMC);  //non HF
    }
  }
    //Fill of count histograms
  if (!fAlreadyFilled && !fMixing) { 
    ((TH1F*)fOutputStudy->FindObject("hist_Count_Charg"))->Fill(N_Charg);
  }

  fAlreadyFilled=kTRUE; //at least a Lambdac analyzed in the event; distribution plots already filled

}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::FillSparsePlots(TClonesArray* mcArray, Double_t mInv[], Int_t origLambdac, Int_t PdgLambdac, AliReducedParticle* track, Int_t ptbin, Int_t type, Int_t softpiME, Double_t wg) {
  //
  //fills the THnSparse for correlations, calculating the variables
  //

  //for MC, in case Ntrkl reweight is active, add the event weights to THnSparse
  if(fReadMC && fUseNtrklWeight) wg*=fWeight;

  //Initialization of variables
  Double_t mLambdac, /* mLambdacbar, */ deltaphi = 0., deltaeta = 0.;
  mLambdac = mInv[0];
  //mLambdacbar = mInv[1];

  if (fReadMC && track->GetLabel()<1) return;
  if (fReadMC && !(AliAODMCParticle*)mcArray->At(track->GetLabel())) return;
  Double_t ptTrack = track->Pt();
  Double_t LambdacTrack = type!=kK0 ? track->GetImpPar() : 0.;
  Double_t phiTr = track->Phi();
  Double_t etaTr = track->Eta();
  Double_t origTr = fReadMC ? CheckTrackOrigin(mcArray,(AliAODMCParticle*)mcArray->At(track->GetLabel())) : 0;

  TString part = "", orig = "";

  switch (type) {
    case(kTrack): {
      part = "Charg";
      deltaphi = fCorrelatorTr->GetDeltaPhi();
      deltaeta = fCorrelatorTr->GetDeltaEta();
      break;
    }
   
  }
  
  if(fMixing == kSE) {

    //Fixes limits; needed to include overflow into THnSparse projections!
    Double_t pTorig = track->Pt();
    Double_t Lambdacorig = track->GetImpPar();
    Double_t ptLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",ptbin,fPoolNum)))->GetAxis(2)->GetXmax(); //all plots have same axes...
    Double_t displLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",ptbin,fPoolNum)))->GetAxis(3)->GetXmax();
    Double_t EtaLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d",ptbin,fPoolNum)))->GetAxis(4)->GetXmax();
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01;
    if(LambdacTrack > displLim_Sparse) LambdacTrack = (displLim_Sparse-0.001);
    if(deltaeta > EtaLim_Sparse) deltaeta = EtaLim_Sparse-0.01;
    if(deltaeta < -EtaLim_Sparse) deltaeta = -EtaLim_Sparse+0.01;
  
    //variables for filling histos
    Double_t fillSpPhiLambdac[5] = {deltaphi,mLambdac,ptTrack,LambdacTrack,deltaeta};
    //Double_t fillSpPhiLambdacbar[5] = {deltaphi,mLambdacbar,ptTrack,LambdacTrack,deltaeta};
    Double_t fillSpWeigLambdac[5] = {deltaphi,mLambdac,deltaeta,ptTrack};
    //Double_t fillSpWeigLambdacbar[5] = {deltaphi,mLambdacbar,deltaeta,ptTrack};

    Bool_t allowLambdac = 0;
    //Bool_t allowLambdacbar = 0;
    if(fSpeed==kOneBinSB) { //filling of sidebands in speed mode: 1 bin for LSB, 1 for RSB, no filling outside signal region and SB
      if(ptbin<PtBinCorr(8.01)) {	    
        if(mLambdac > fSignLeft_LowPt && mLambdac < fSignRight_LowPt) allowLambdac = 1;
        //if(mLambdacbar > fSignLeft_LowPt && mLambdacbar < fSignRight_LowPt) allowLambdacbar = 1;
      } else {
        if(mLambdac > fSignLeft_HighPt && mLambdac < fSignRight_HighPt) allowLambdac = 1;
        //if(mLambdacbar > fSignLeft_HighPt && mLambdacbar < fSignRight_HighPt) allowLambdacbar = 1;
      }
      if(mLambdac > fLSBLowLim.at(ptbin) && mLambdac < fLSBUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 1.99; fillSpWeigLambdac[1] = 1.99;} //in LSB bin!
      //if(mLambdacbar > fLSBLowLim.at(ptbin) && mLambdacbar < fLSBUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 1.99; fillSpWeigLambdacbar[1] = 1.99;} //in LSB bin!
      if(mLambdac > fRSBLowLim.at(ptbin) && mLambdac < fRSBUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 2.585; fillSpWeigLambdac[1] = 2.585;} //in RSB bin!
      //if(mLambdacbar > fRSBLowLim.at(ptbin) && mLambdacbar < fRSBUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 2.585; fillSpWeigLambdacbar[1] = 2.585;} //in RSB bin!
    } //in this way if sidebands overlap with signal range in Mass axis, those overlapping bins will be void. But this creates no problems...
    if(fSpeed==kOneBinSBandS) { //filling of sidebands in speed mode: 1 bin for LSB, 1 for RSB, no filling outside signal region and SB
      if(mLambdac > fSignLowLim.at(ptbin) && mLambdac < fSignUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 0; fillSpWeigLambdac[1] = 0;} //in Signal region-->Bin 1 (-0.5->0.5)!
      //if(mLambdacbar > fSignLowLim.at(ptbin) && mLambdacbar < fSignUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 0; fillSpWeigLambdacbar[1] = 0;} //in Signal region-->Bin 1 (-0.5->0.5)!
      if(mLambdac > fLSBLowLim.at(ptbin) && mLambdac < fLSBUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 1; fillSpWeigLambdac[1] = 1;} //in LSB!-->Bin 2 (0.5->1.5)!
      //if(mLambdacbar > fLSBLowLim.at(ptbin) && mLambdacbar < fLSBUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 1; fillSpWeigLambdacbar[1] = 1;} //in LSB!-->Bin 2 (0.5->1.5)!
      if(mLambdac > fRSBLowLim.at(ptbin) && mLambdac < fRSBUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 1; fillSpWeigLambdac[1] = 1;} //in RSB!-->Bin 2 (0.5->1.5)!
      //if(mLambdacbar > fRSBLowLim.at(ptbin) && mLambdacbar < fRSBUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 1; fillSpWeigLambdacbar[1] = 1;} //in RSB!-->Bin 2 (0.5->1.5)!
    }    
    if(!fSpeed) { // Full Minv range in THnSparse!
      /* if((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3)) allowLambdac = 1;   
      if((fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3)) allowLambdacbar = 1; */
      if(fIsSelectedCandidate > 0) allowLambdac = 1;
    }
    
    if(fReadMC == 0) {
      //sparse fill for data (tracks, K+-, K0) + weighted

      /* if((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) && allowLambdac) { //Lambdac //
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
        if(!fSpeed) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
      }
      if((fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) && allowLambdacbar) { //Lambdacbar
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdacbar,wg);
        if(!fSpeed) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigLambdacbar,pTorig*wg);
      } */
      if(fIsSelectedCandidate > 0 && allowLambdac) { //Lambdac //
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
        if(!fSpeed) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
      }
      if(!fAlreadyFilled) {
 	      ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_%s_Bin%d",part.Data(),ptbin)))->Fill(pTorig);
 	      ((TH1F*)fOutputStudy->FindObject(Form("hist_PhiDistr_%s",part.Data())))->Fill(phiTr);
 	      ((TH1F*)fOutputStudy->FindObject(Form("hist_EtaDistr_%s",part.Data())))->Fill(etaTr);
 	      ((TH2F*)fOutputStudy->FindObject(Form("hist_PhiVsEtaDistr_%s",part.Data())))->Fill(phiTr,etaTr);
      }
    }

    if(fReadMC) {

      if(origLambdac==4) {orig = "_From_c";} else {orig = "_From_b";}

      //sparse fill for data (tracks, K+-, K0) + weighted
      /* if(PdgLambdac==4122 && (fIsSelectedCandidate==1||fIsSelectedCandidate==3)) { //Lambdac (from MCTruth)
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
         if(origLambdac==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
         if(origLambdac==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_NonHF_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
         if(!fSpeed) {
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
           if(origLambdac==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
           if(origLambdac==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
           if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_NonHF_Bin%d",ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
         }  
      }
      if(PdgLambdac==-4122 && fIsSelectedCandidate>1) { //Lambdacbar
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdacbar,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdacbar,wg);
         if(origLambdac==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdacbar,wg);
         if(origLambdac==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdacbar,wg); 
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_NonHF_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdacbar,wg);
         if(!fSpeed) {
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigLambdacbar,pTorig*wg);
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdacbar,pTorig*wg);
           if(origLambdac==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdacbar,pTorig*wg);
           if(origLambdac==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdacbar,pTorig*wg);
           if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_NonHF_Bin%d",ptbin)))->Fill(fillSpWeigLambdacbar,pTorig*wg);
         }
      }  */
      if(std::abs(PdgLambdac)==4122 && fIsSelectedCandidate>0) {
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpWeigLambdac,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpWeigLambdac,wg);
         if(origLambdac==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpWeigLambdac,wg);
         if(origLambdac==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d_p%d",part.Data(),orig.Data(),ptbin,fPoolNum)))->Fill(fillSpWeigLambdac,wg); 
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_NonHF_Bin%d_p%d",part.Data(),ptbin,fPoolNum)))->Fill(fillSpWeigLambdac,wg);
         if(!fSpeed) {
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
           ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
           if(origLambdac==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
           if(origLambdac==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
           if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_NonHF_Bin%d",ptbin)))->Fill(fillSpWeigLambdac,pTorig*wg);
         }
      }
      if(!fAlreadyFilled) {
	      ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_Bin%d",part.Data(),ptbin)))->Fill(Lambdacorig); //Fills displacement histos
        if (origTr>=1&&origTr<=8) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_HF_Bin%d",part.Data(),ptbin)))->Fill(Lambdacorig);
        if (origTr>=1&&origTr<=8) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_HF%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(Lambdacorig);
        ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(Lambdacorig); //Fills displacement histos
        ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_%s_Bin%d",part.Data(),ptbin)))->Fill(pTorig);
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_%s_Bin%d",part.Data(),ptbin)))->Fill(origTr);
 	      ((TH1F*)fOutputStudy->FindObject(Form("hist_PhiDistr_%s",part.Data())))->Fill(phiTr);
 	      ((TH1F*)fOutputStudy->FindObject(Form("hist_PhiDistr_%s",part.Data())))->Fill(phiTr);
 	      ((TH1F*)fOutputStudy->FindObject(Form("hist_EtaDistr_%s",part.Data())))->Fill(etaTr);
 	      ((TH2F*)fOutputStudy->FindObject(Form("hist_PhiVsEtaDistr_%s",part.Data())))->Fill(phiTr,etaTr);
      }

      if(fPurityStudies) FillPurityPlots(mcArray,track,ptbin,deltaphi);
    }//end MC case

  } //end of SE fill

  if(fMixing == kME) {

    //Fixes limits; needed to include overflow into THnSparse projections!
    Double_t EtaLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d_EvMix",ptbin,fPoolNum)))->GetAxis(2)->GetXmax();
    if(deltaeta > EtaLim_Sparse) deltaeta = EtaLim_Sparse-0.01;
    if(deltaeta < -EtaLim_Sparse) deltaeta = -EtaLim_Sparse+0.01;
    Double_t ptLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_p%d_EvMix",ptbin,fPoolNum)))->GetAxis(3)->GetXmax(); //all plots have same axes...
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01;

    //variables for filling histos
    Double_t fillSpPhiLambdac[5] = {deltaphi,mLambdac,deltaeta,0.4,0}; //dummy for ME threshold! unless explicitly set by flag...
    //Double_t fillSpPhiLambdacbar[5] = {deltaphi,mLambdacbar,deltaeta,0.4,0};
    if(fMEAxisThresh) {
      fillSpPhiLambdac[3] = ptTrack;
      //fillSpPhiLambdacbar[3] = ptTrack;
    }
   
    if(softpiME==1) { //it's a softPi in the ME analysis! Fill it in the dedicated slice of ME THnSparse
          fillSpPhiLambdac[4] = 1;
          //fillSpPhiLambdacbar[4] = 1;
    	}
    	
    Bool_t allowLambdac = 0;
    //Bool_t allowLambdacbar = 0;
    if(fSpeed==kOneBinSB) { //filling of sidebands in speed mode: 1 bin for LSB, 1 for RSB, no filling outside signal region and SB
      if(ptbin<PtBinCorr(8.01)) {	    
        if(mLambdac > fSignLeft_LowPt && mLambdac < fSignRight_LowPt) allowLambdac = 1;
        //if(mLambdacbar > fSignLeft_LowPt && mLambdacbar < fSignRight_LowPt) allowLambdacbar = 1;
      } else {
        if(mLambdac > fSignLeft_HighPt && mLambdac < fSignRight_HighPt) allowLambdac = 1;
        //if(mLambdacbar > fSignLeft_HighPt && mLambdacbar < fSignRight_HighPt) allowLambdacbar = 1;
      }
      if(mLambdac > fLSBLowLim.at(ptbin) && mLambdac < fLSBUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 1.99;} //in LSB bin!
      //if(mLambdacbar > fLSBLowLim.at(ptbin) && mLambdacbar < fLSBUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 1.99;} //in LSB bin!
      if(mLambdac > fRSBLowLim.at(ptbin) && mLambdac < fRSBUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 2.585;} //in RSB bin!
      //if(mLambdacbar > fRSBLowLim.at(ptbin) && mLambdacbar < fRSBUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 2.585;} //in RSB bin!
    } //in this way if sidebands overlap with signal range in Mass axis, those overlapping bins will be void. But this creates no problems...
    if(fSpeed==kOneBinSBandS) { //filling of sidebands in speed mode: 1 bin for LSB, 1 for RSB, no filling outside signal region and SB
      if(mLambdac > fSignLowLim.at(ptbin) && mLambdac < fSignUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 0;} //in Signal region-->Bin 1 (-0.5->0.5)!
      //if(mLambdacbar > fSignLowLim.at(ptbin) && mLambdacbar < fSignUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 0;} //in Signal region-->Bin 1 (-0.5->0.5)!
      if(mLambdac > fLSBLowLim.at(ptbin) && mLambdac < fLSBUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 1;} //in LSB!-->Bin 2 (0.5->1.5)!
      //if(mLambdacbar > fLSBLowLim.at(ptbin) && mLambdacbar < fLSBUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 1;} //in LSB!-->Bin 2 (0.5->1.5)!
      if(mLambdac > fRSBLowLim.at(ptbin) && mLambdac < fRSBUppLim.at(ptbin)) {allowLambdac = 1; fillSpPhiLambdac[1] = 1;} //in RSB!-->Bin 2 (0.5->1.5)!
      //if(mLambdacbar > fRSBLowLim.at(ptbin) && mLambdacbar < fRSBUppLim.at(ptbin)) {allowLambdacbar = 1; fillSpPhiLambdacbar[1] = 1;} //in RSB!-->Bin 2 (0.5->1.5)!
    }      
    if(!fSpeed) { // Full Minv range in THnSparse!
      /* if((fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3)) allowLambdac = 1;   
      if((fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3)) allowLambdacbar = 1; */
      allowLambdac = 1;
    }    
    
    if(fReadMC == 0) {
      //sparse fill for data (tracks, K+-, K0)
      /* if((fIsSelectedCandidate == 1||fIsSelectedCandidate == 3) && allowLambdac) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
      if((fIsSelectedCandidate == 2||fIsSelectedCandidate == 3) && allowLambdacbar) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdacbar,wg); */
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
    }
    if(fReadMC == 1) {
      //sparse fill for data (tracks, K+-, K0)
      /* if(PdgLambdac==4122 && (fIsSelectedCandidate==1||fIsSelectedCandidate==3))  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
      if(PdgLambdac==-4122 && fIsSelectedCandidate>1) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdacbar,wg); */
      if(std::abs(PdgLambdac)==4122 && fIsSelectedCandidate>0)((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_p%d_EvMix",part.Data(),ptbin,fPoolNum)))->Fill(fillSpPhiLambdac,wg);
    }//end MC case

  } //end of ME fill
  
  return;
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSELcTopK0sCorrelations::CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checks on particle (#) origin:
  // 0) Not HF
  // 1) D->#
  // 2) D->X->#
  // 3) c hadronization
  // 4) B->#
  // 5) B->X-># (X!=D)
  // 6) B->D->#
  // 7) B->D->X->#
  // 8) b hadronization
  //
  if(fDebug>2) printf("AliAnalysisTaskSELcTopK0sCorrelations::CheckTrkOrigin() \n");
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isDdaugh=kFALSE;
  Bool_t isDchaindaugh=kFALSE;
  Bool_t isBdaugh=kFALSE;
  Bool_t isBchaindaugh=kFALSE;
  Bool_t isQuarkFound=kFALSE;

  if (mother<0) return -1;
  while (mother >= 0){
    istep++;
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcMoth){
      pdgGranma = mcMoth->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isBchaindaugh=kTRUE;
        if(istep==1) isBdaugh=kTRUE;
      }
      if ((abspdgGranma > 400 && abspdgGranma < 500) || (abspdgGranma > 4000 && abspdgGranma < 5000)){
	isDchaindaugh=kTRUE;
        if(istep==1) isDdaugh=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) {isQuarkFound=kTRUE; if(abspdgGranma==5) isFromB = kTRUE;}
      mother = mcMoth->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      return -1;
    }
  }

  //decides what to return based on the flag status
  if(isQuarkFound) {
    if(!isFromB) {  //charm
      if(isDdaugh) return 1; //charm immediate
      else if(isDchaindaugh) return 2; //charm chain
      else return 3; //charm hadronization
    }
    else { //beauty
      if(isBdaugh) return 4; //b immediate
      else if(isBchaindaugh) { //b chain
        if(isDchaindaugh) {
          if(isDdaugh) return 6; //d immediate
          return 7; //d chain
          }
        else return 5; //b, not d
      }
      else return 8; //b hadronization
    }
  }
  else if(!isDdaugh && !isDchaindaugh && !isBdaugh && !isBchaindaugh) return 0; //no HF decay 
     //in this case, it's !isQuarkFound, but not in 100% cases it's a non HF particle!
     //rarely you can find a D/B meson which comes from a -1! It isn't a Non-HF, in that case! And I'll return -1...

  return -1; //some problem spotted
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSELcTopK0sCorrelations::IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* mcArray) const {

  //Daughter removal in MCKine case
  Bool_t isDaughter = kFALSE;
  Int_t labelLambdac = d->GetLabel();

  Int_t mother = track->GetMother();

  //Loop on the mothers to find the Lambdac label (it must be the trigger Lambdac, not a generic Lambdac!)
  while (mother > 0){
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother)); //it's the mother of the track!
    if (mcMoth){
      if (mcMoth->GetLabel() == labelLambdac) isDaughter = kTRUE;
      mother = mcMoth->GetMother(); //goes back by one
    } else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  return isDaughter;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSELcTopK0sCorrelations::PtBinCorr(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin=-1;
  if(pt<fBinLimsCorr.at(0)) return ptbin; //out of bounds
  
  Int_t i = 0;
  while(pt>fBinLimsCorr.at(i)) {ptbin=i; i++;}
  
  return ptbin;
}


//____________________________________________________________________________
TProfile* AliAnalysisTaskSELcTopK0sCorrelations::GetEstimatorHistogram(const AliVEvent* event){
  /// Get Estimator Histogram from period event->GetRunNumber();
  ///
  /// If you select SPD tracklets in |eta|<1 you should use type == 1
  ///
    
  Int_t runNo  = event->GetRunNumber();
  Int_t period = -1;
  
  //2016
    if(runNo>=252235 && runNo<=252375)period = 0;//16d
    if(runNo>=252603 && runNo<=253591)period = 1;//16e
    if(runNo>=254124 && runNo<=254332)period = 2;//16g
    if(runNo>=254378 && runNo<=255469)period = 3;//16h
    if(runNo>=256146 && runNo<=256420)period = 4;//16j
    if(runNo>=256504 && runNo<=258537)period = 5;//16k
    if(runNo>=258883 && runNo<=260187)period = 6;//16l
    if(runNo>=262395 && runNo<=264035)period = 7;//16o
    if(runNo>=264076 && runNo<=264347)period = 8;//16p
  //2017
    if(runNo>=270531 && runNo<=270667)period = 9;//17c
    if(runNo>=270822 && runNo<=270830)period = 10;//17e
    if(runNo>=270854 && runNo<=270865)period = 11;//17f
    if(runNo>=271868 && runNo<=273103)period = 12;//17h
    if(runNo>=273591 && runNo<=274442)period = 13;//17i
    if(runNo>=274593 && runNo<=274671)period = 14;//17j 
    if(runNo>=274690 && runNo<=276508)period = 15;//17k
    if(runNo>=276551 && runNo<=278216)period = 16;//17l
    if(runNo>=278914 && runNo<=280140)period = 17;//17m
    if(runNo>=280282 && runNo<=281961)period = 18;//17o
    if(runNo>=282504 && runNo<=282704)period = 19;//17r
  //2018
    if(runNo>=284706 && runNo<=285447)period = 20;//18b
    if(runNo>=285978 && runNo<=286350)period = 21;//18d
    if(runNo>=286380 && runNo<=286937)period = 22;//18e
    if(runNo>=287000 && runNo<=287977)period = 23;//18f
    if(runNo>=288619 && runNo<=288750)period = 24;//18g
    if(runNo>=288804 && runNo<=288806)period = 25;//18h
    if(runNo>=288861 && runNo<=288909)period = 26;//18i
    if(runNo>=289165 && runNo<=289201)period = 27;//18k
    if(runNo>=289240 && runNo<=289971)period = 28;//18l
    if(runNo>=290222 && runNo<=292839)period = 29;//18m
    if(runNo>=293357 && runNo<=293359)period = 30;//18n
    if(runNo>=293368 && runNo<=293898)period = 31;//18o
    if(runNo>=294009 && runNo<=294925)period = 32;//18p  
  
  if(period==-1) {
     printf("Error! No corresponding profile for tracklets for this run (%d)! Skipping the correction...\n",runNo);
     return 0x0;
  }

  return fTrackletProfiles[period];
}

//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::FillPurityPlots(TClonesArray* mcArray, AliReducedParticle* track, Int_t ptbin, Double_t deltaphi) {

  //Purity studies (only in MC reco mode)

  if(!fReadMC || !fRecoLambdac || !fRecoTr) return;

  TString namebinD[6] = {"2to3","3to5","5to8","8to16","16to24","24to36"};
  TString namebinAss[7] = {"03to99","03to1","1to99","1to3","1to2","2to3","3to99"};

  AliAODMCParticle* trkKine = (AliAODMCParticle*)mcArray->At(track->GetLabel());
  if (!trkKine) return;
  Double_t origTr = CheckTrackOrigin(mcArray,trkKine);  
  Bool_t primTrack = trkKine->IsPhysicalPrimary();
  Double_t pTtr = track->Pt();

  Bool_t fillAssocRange[7] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
  TString stringpTD = "";
  Bool_t okpTD = kFALSE;
  if(fBinLimsCorr.at(ptbin) >= 2 && fBinLimsCorr.at(ptbin) < 3)   {stringpTD = namebinD[0]; okpTD = kTRUE;}
  if(fBinLimsCorr.at(ptbin) >= 3 && fBinLimsCorr.at(ptbin) < 5)   {stringpTD = namebinD[1]; okpTD = kTRUE;}
  if(fBinLimsCorr.at(ptbin) >= 5 && fBinLimsCorr.at(ptbin) < 8)   {stringpTD = namebinD[2]; okpTD = kTRUE;}
  if(fBinLimsCorr.at(ptbin) >= 8 && fBinLimsCorr.at(ptbin) < 16)  {stringpTD = namebinD[3]; okpTD = kTRUE;}
  if(fBinLimsCorr.at(ptbin) >= 16 && fBinLimsCorr.at(ptbin) < 24) {stringpTD = namebinD[4]; okpTD = kTRUE;}
  if(fBinLimsCorr.at(ptbin) >= 24 && fBinLimsCorr.at(ptbin) < 36) {stringpTD = namebinD[5]; okpTD = kTRUE;}

  if(pTtr >= 0.3) fillAssocRange[0] = kTRUE;
  if(pTtr >= 0.3 && pTtr < 1) fillAssocRange[1] = kTRUE;
  if(pTtr >= 1) fillAssocRange[2] = kTRUE;
  if(pTtr >= 1 && pTtr < 3) fillAssocRange[3] = kTRUE;
  if(pTtr >= 1 && pTtr < 2) fillAssocRange[4] = kTRUE;
  if(pTtr >= 2 && pTtr < 3) fillAssocRange[5] = kTRUE;
  if(pTtr >= 3) fillAssocRange[6] = kTRUE;

  if(!okpTD) return;
  for(int j=0; j<7; j++) {
    if(fillAssocRange[j]==kTRUE) {
      if(primTrack) {
        ((TH1F*)fOutputStudy->FindObject(Form("hPurityCount_PrimAccepted_pTD%s_pTass%s",stringpTD.Data(),namebinAss[j].Data())))->Fill(0.,fWeight); 
        ((TH1F*)fOutputStudy->FindObject(Form("hPuritydPhi_PrimAccepted_pTD%s_pTass%s",stringpTD.Data(),namebinAss[j].Data())))->Fill(deltaphi,fWeight); 
      }
      if(!primTrack) {
        ((TH1F*)fOutputStudy->FindObject(Form("hPurityCount_SecAccepted_pTD%s_pTass%s",stringpTD.Data(),namebinAss[j].Data())))->Fill(0.,fWeight); 
        ((TH1F*)fOutputStudy->FindObject(Form("hPuritydPhi_SecAccepted_pTD%s_pTass%s",stringpTD.Data(),namebinAss[j].Data())))->Fill(deltaphi,fWeight); 
      }
      if(primTrack && origTr>=1&&origTr<=3) {  //fill for acccepted primary charm tracks
        ((TH1F*)fOutputStudy->FindObject(Form("hPurityCount_CharmAccepted_pTD%s_pTass%s",stringpTD.Data(),namebinAss[j].Data())))->Fill(0.,fWeight); 
        ((TH1F*)fOutputStudy->FindObject(Form("hPuritydPhi_CharmAccepted_pTD%s_pTass%s",stringpTD.Data(),namebinAss[j].Data())))->Fill(deltaphi,fWeight); 
      }
      if(primTrack && origTr>=4&&origTr<=8) {  //fill for accepted primary beauty tracks
        ((TH1F*)fOutputStudy->FindObject(Form("hPurityCount_BeautyAccepted_pTD%s_pTass%s",stringpTD.Data(),namebinAss[j].Data())))->Fill(0.,fWeight); 
        ((TH1F*)fOutputStudy->FindObject(Form("hPuritydPhi_BeautyAccepted_pTD%s_pTass%s",stringpTD.Data(),namebinAss[j].Data())))->Fill(deltaphi,fWeight); 
      }
    }
  }
  
  return;
}

//__________________________________________________________________________________________________
Double_t AliAnalysisTaskSELcTopK0sCorrelations::GetNtrklWeight(Int_t ntrkl){
  //
  //  extracts the Ntrkl weight using the histo with data/MC Ntracklets ratio
  //
  if(ntrkl<=0) return 1.;
  if(!fHistNtrklWeight) { AliError("Input histogram to evaluate Ntrkl weight missing"); return 0.; }
  Double_t histweight=fHistNtrklWeight->GetBinContent(fHistNtrklWeight->FindBin(ntrkl));
  Double_t weight = histweight>0 ? histweight : 0.;
  return weight;
}


//________________________________________________________________________
void AliAnalysisTaskSELcTopK0sCorrelations::PrintBinsAndLimits() {

  cout << "--------------------------\n";
  cout << "PtBins = " << fNPtBinsCorr << "\n";
  cout << "PtBin limits--------------\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << "Bin "<<i+1<<" = "<<fBinLimsCorr.at(i)<<" to "<<fBinLimsCorr.at(i+1)<<"\n";
  }
  cout << "\n--------------------------\n";
  cout << "PtBin tresh. tracks low---\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fPtThreshLow.at(i) << ", ";
  }
  cout << "\nPtBin tresh. tracks up----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fPtThreshUp.at(i) << ", ";
  }
  cout << "\nSB limits (LSBLow)----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fLSBLowLim.at(i) << ", ";
  }
  cout << "\nSB limits (LSBUpp)----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fLSBUppLim.at(i) << ", ";
  }
  cout << "\nSB limits (RSBLow)----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fRSBLowLim.at(i) << ", ";
  }
  cout << "\nSB limits (RSBUpp)----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fRSBUppLim.at(i) << ", ";
  }
  if(fSpeed==kOneBinSBandS) {
    cout << "\nSignal limits (Low)----\n";
    for (int i=0; i<fNPtBinsCorr; i++) {
      cout << fSignLowLim.at(i) << ", ";
    }
    cout << "\nSignal limits (Upp)----\n";
    for (int i=0; i<fNPtBinsCorr; i++) {
      cout << fSignUppLim.at(i) << ", ";
    }        
  }

  cout << "\n--------------------------\n";
  cout << "Lambdac Eta cut for Correl = "<<fEtaForCorrel<<"\n";
  cout << "--------------------------\n";
  cout << "MC Truth = "<<fReadMC<<" - MC Reco: Trk = "<<fRecoTr<<", Lambdac = "<<fRecoLambdac<<"\n";
  cout << "--------------------------\n";
  cout << "Sel of Event tpye (PP,GS,FE,...)= "<<fSelEvType<<"\n";
  cout << "--------------------------\n";
  cout << "Ev Mixing = "<<fMixing<<"\n";
  cout << "--------------------------\n";
  cout << "ME thresh axis = "<<fMEAxisThresh<<"\n";
  cout << "--------------------------\n";
  cout << "Speed (1 SBL/SBR and eventually Sign bin) = "<<fSpeed<<"\n";
  cout << "--------------------------\n";
  cout << "All entries in Pool0 = "<<fMergePools<<"\n";
  cout << "--------------------------\n";  
  cout << "PtBin associated maximum edge = "<<fPtAssocLimit<<"\n";
  cout << "--------------------------\n";    
  cout << "Minimum D-meson pT = "<<fMinDPt<<"\n";
  cout << "--------------------------\n";  
  cout << "TTree filling = "<<fFillTrees<<"\n";
  cout << "--------------------------\n";  
  cout << "Purity studies (for MC) = "<<fPurityStudies<<"\n";
  cout << "--------------------------\n";  
  cout << "Ntrkl reweighting (for MC) = "<<fUseNtrklWeight<<"\n";
  cout << "--------------------------\n";  
  cout << "Analysis vs Mult = "<<fVsMultAnalysis<<"\n";
  cout << "--------------------------\n";  
  if(fVsMultAnalysis) cout << "V0M lims = "<<fV0CentMin<<"-"<<fV0CentMax<<"  -  SPD trkl lims = "<<fTrkMultMin<<"-"<<fTrkMultMax<<"\n";
  if(fVsMultAnalysis) cout << "--------------------------\n";  

  if(fPurityStudies) {
  cout << "---------------------------------------------------------------------------------------------------------------------------------------------------\n";  
  cout << "WARNING! Task launched in 'purity studies' mode! If on MC, secondary tracks at reco will NOT be discarded, in the THnSparse! On data, no effects...\n";
  cout << "---------------------------------------------------------------------------------------------------------------------------------------------------\n"; 
  getchar();
  }

}

bool AliAnalysisTaskSELcTopK0sCorrelations::CheckDaugAcc(TClonesArray *arrayMC, int nProng, int *labDau)
{
    /// check if the decay products are in the good eta and pt range

    for (int iProng = 0; iProng < nProng; iProng++)
    {
        AliAODMCParticle *mcPartDaughter = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labDau[iProng]));
        if (!mcPartDaughter)
            return false;

        double eta = mcPartDaughter->Eta();

        if (TMath::Abs(eta) > 0.8)
            return false;
    }
    return true;
}

bool AliAnalysisTaskSELcTopK0sCorrelations::ReconstructKFLc(AliAODRecoCascadeHF *cand){
// topological variables
    AliAODRecoCascadeHF* candCasc = (AliAODRecoCascadeHF*)cand;
    AliAODv0 *v0part = candCasc->Getv0();
    AliAODTrack *bachPart = dynamic_cast<AliAODTrack*>(candCasc->GetBachelor());
    AliAODVertex *primVert = dynamic_cast<AliAODVertex*>(candCasc->GetPrimaryVtx());
    float massK0s = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    float massL = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    //Use KF Particle
        // convert primary vertex
        KFPVertex pVertex;
        double pos[3], cov[6];
        primVert->GetXYZ(pos);
        primVert->GetCovarianceMatrix(cov);
        float posF[3], covF[6];
        for(int iEl = 0; iEl < 3; iEl++)
            posF[iEl] = (float)pos[iEl];
        for(int iEl = 0; iEl < 6; iEl++)
            covF[iEl] = (float)cov[iEl];

        pVertex.SetXYZ((float)pos[0], (float)pos[1], (float)pos[2]);
        pVertex.SetCovarianceMatrix(covF);
        pVertex.SetChi2(primVert->GetChi2());
        pVertex.SetNDF(primVert->GetNDF());
        pVertex.SetNContributors(primVert->GetNContributors());
        KFParticle primVertKF(pVertex);

        double covP[21], covN[21], covB[21];

        int pdgBach = -1;
        int pdgV0dau[2] = {-1, -1};
        /* if(massHypo == AliHFMLVarHandlerNonPromptLc2V0bachelor::kpK0s)
        { */
            pdgBach = 2212;
            pdgV0dau[0] = 211;
            pdgV0dau[1] = 211;

        AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(candCasc->Getv0PositiveTrack());
        AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(candCasc->Getv0NegativeTrack());
        // check charge of the first daughter, if negative, define it as the second one
        if (v0Pos->Charge()<0) {
            v0Pos = dynamic_cast<AliAODTrack*>(candCasc->Getv0NegativeTrack());
            v0Neg = dynamic_cast<AliAODTrack*>(candCasc->Getv0PositiveTrack());
        }

        if ( !bachPart->GetCovarianceXYZPxPyPz(covB) || !v0Pos->GetCovarianceXYZPxPyPz(covP) || !v0Neg->GetCovarianceXYZPxPyPz(covN) )
            return false;
        if ( !AliVertexingHFUtils::CheckAODtrackCov(bachPart) || !AliVertexingHFUtils::CheckAODtrackCov(v0Pos) || !AliVertexingHFUtils::CheckAODtrackCov(v0Neg) )
            return false;

        KFParticle KFBachelor;
        if(bachPart->Charge() > 0) KFBachelor = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, pdgBach);
        if(bachPart->Charge() < 0) KFBachelor = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, -pdgBach);

        KFParticle KFV0DauPlus, KFV0DauMinus;

            KFV0DauPlus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, pdgV0dau[0]);
            KFV0DauMinus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -pdgV0dau[1]);

        KFParticle KFV0;
        const KFParticle *Ks0Daughters[2] = {&KFV0DauPlus, &KFV0DauMinus};
        KFV0.Construct(Ks0Daughters, 2);
        // check V0 covariance matrix
        if(!AliVertexingHFUtils::CheckKFParticleCov(KFV0))
            return false;

        float dcaPointV0[8], dcaPointV0Cov[36];
        KFV0.GetParametersAtPoint(posF, covF, dcaPointV0, dcaPointV0Cov);
        KFV0.SetNonlinearMassConstraint(massK0s);

        // reconstruct Lc
        KFParticle KFLc;
        const KFParticle *LcDaughters[2] = {&KFBachelor, &KFV0};
        KFLc.Construct(LcDaughters, 2);

        // check Lc covariance matrix
        if (!AliVertexingHFUtils::CheckKFParticleCov(KFLc))
            return false;

        fKFLcInvMass = KFLc.GetMass();
        fKFLcPt = KFLc.GetPt();
        fKFLcEta = KFLc.GetEta();  
        fKFLcPhi = KFLc.GetPhi();  
        fKFLcRapidity = KFLc.GetRapidity();

        return true;

}