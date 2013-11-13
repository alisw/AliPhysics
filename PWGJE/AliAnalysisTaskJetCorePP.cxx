
// ******************************************
// This task computes several jet observables like 
// the fraction of energy in inner and outer coronnas,
// jet-track correlations,triggered jet shapes and 
// correlation strength distribution of particles inside jets.    
// Author: lcunquei@cern.ch
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


#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TArrayI.h" 
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TRandom3.h"

#include "AliLog.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODJet.h"

#include "AliAnalysisTaskJetCorePP.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetCorePP)

//Filip Krizek 1st March 2013

//---------------------------------------------------------------------
AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP() :
AliAnalysisTaskSE(),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
fAODExtension(0x0),
fJetBranchName(""),
fJetBranchNameChargMC(""),
fJetBranchNameFullMC(""),
fJetBranchNameBg(""),
fJetBranchNameBgChargMC(""),
fListJets(0x0),
fListJetsGen(0x0),
fListJetsGenFull(0x0),
fListJetsBg(0x0),
fListJetsBgGen(0x0),
fNonStdFile(""),
fSystem(0), //pp=0  pPb=1
fJetParamR(0.4),
fBgJetParamR(0.3),
fBgMaxJetPt(8.0),
fBgConeR(0.4),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.0),
fVtxZMax(10.0),
fFilterMask(0),
fCentMin(0.0),
fCentMax(100.0),
fJetEtaMin(-0.5),
fJetEtaMax(0.5),
fTriggerEtaCut(0.9),
fTrackEtaCut(0.9),
fTrackLowPtCut(0.15),
fOutputList(0x0),
fHistEvtSelection(0x0),
fh2Ntriggers(0x0),
fHJetSpec(0x0),
fHJetSpecSubUeMedian(0x0),
fHJetSpecSubUeCone(0x0),
fHJetUeMedian(0x0),
fHJetUeCone(0x0),
fHRhoUeMedianVsCone(0x0),
//fHJetDensity(0x0),
//fHJetDensityA4(0x0),
fhJetPhi(0x0),
fhTriggerPhi(0x0),
fhJetEta(0x0),
fhTriggerEta(0x0),
fhVertexZ(0x0),
fhVertexZAccept(0x0),
fhContribVtx(0x0),
fhContribVtxAccept(0x0),
fhDphiTriggerJet(0x0),
fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0),
fhCentralityAccept(0x0),
//fHJetPtRaw(0x0),
//fHLeadingJetPtRaw(0x0), 
//fHDphiVsJetPtAll(0x0), 
fhJetPtGenVsJetPtRec(0x0),
fhJetPtGenVsJetPtRecSubUeMedian(0x0),
fhJetPtGenVsJetPtRecSubUeCone(0x0),
fhJetPtGen(0x0), 
fhJetPtSubUeMedianGen(0x0), 
fhJetPtSubUeConeGen(0x0), 
fhJetPtGenChargVsJetPtGenFull(0x0),
fhJetPtGenFull(0x0),
fh2NtriggersGen(0x0),
fHJetSpecGen(0x0),
fHJetSpecSubUeMedianGen(0x0),
fHJetSpecSubUeConeGen(0x0),
fHJetUeMedianGen(0x0),
fHJetUeConeGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkSecOrFakeRec(0x0),
fHRhoUeMedianVsConeGen(0x0),
fIsChargedMC(0),
fIsFullMC(0),
faGenIndex(0),
faRecIndex(0),
fkAcceptance(2.0*TMath::Pi()*1.8),
fkDeltaPhiCut(TMath::Pi()-0.8),
fh1Xsec(0x0),
fh1Trials(0x0),
fh1AvgTrials(0x0),
fh1PtHard(0x0),
fh1PtHardNoW(0x0),  
fh1PtHardTrials(0x0),
fAvgTrials(1),
fHardest(0),
fEventNumberRangeLow(0),
fEventNumberRangeHigh(99),
fTriggerPtRangeLow(0.0),
fTriggerPtRangeHigh(10000.0),
fFillRespMx(0),
fRandom(0x0),
fnTrials(1000),
fnPhi(9),
fnEta(2),
fEtaSize(0.7),
fPhiSize(2*TMath::Pi()/fnPhi),
fCellArea(fPhiSize*fEtaSize),
fSafetyMargin(1.1)
{
   // default Constructor
}

//---------------------------------------------------------------------

AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
fAODExtension(0x0),
fJetBranchName(""),
fJetBranchNameChargMC(""),
fJetBranchNameFullMC(""),
fJetBranchNameBg(""),
fJetBranchNameBgChargMC(""),
fListJets(0x0),
fListJetsGen(0x0),
fListJetsGenFull(0x0),
fListJetsBg(0x0),
fListJetsBgGen(0x0),
fNonStdFile(""),
fSystem(0),  //pp=0   pPb=1
fJetParamR(0.4),
fBgJetParamR(0.3),
fBgMaxJetPt(8.0),
fBgConeR(0.4),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.0),
fVtxZMax(10.0),
fFilterMask(0),
fCentMin(0.0),
fCentMax(100.0),
fJetEtaMin(-0.5),
fJetEtaMax(0.5),
fTriggerEtaCut(0.9),
fTrackEtaCut(0.9),
fTrackLowPtCut(0.15),
fOutputList(0x0),
fHistEvtSelection(0x0),
fh2Ntriggers(0x0),
fHJetSpec(0x0),
fHJetSpecSubUeMedian(0x0),
fHJetSpecSubUeCone(0x0),
fHJetUeMedian(0x0),
fHJetUeCone(0x0),
fHRhoUeMedianVsCone(0x0),
//fHJetDensity(0x0),
//fHJetDensityA4(0x0),
fhJetPhi(0x0),
fhTriggerPhi(0x0),
fhJetEta(0x0),
fhTriggerEta(0x0),
fhVertexZ(0x0),
fhVertexZAccept(0x0),
fhContribVtx(0x0),
fhContribVtxAccept(0x0),
fhDphiTriggerJet(0x0),
fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0),
fhCentralityAccept(0x0),
//fHJetPtRaw(0x0),
//fHLeadingJetPtRaw(0x0), 
//fHDphiVsJetPtAll(0x0), 
fhJetPtGenVsJetPtRec(0x0),
fhJetPtGenVsJetPtRecSubUeMedian(0x0),
fhJetPtGenVsJetPtRecSubUeCone(0x0),
fhJetPtGen(0x0),
fhJetPtSubUeMedianGen(0x0), 
fhJetPtSubUeConeGen(0x0), 
fhJetPtGenChargVsJetPtGenFull(0x0),
fhJetPtGenFull(0x0),
fh2NtriggersGen(0x0),
fHJetSpecGen(0x0),
fHJetSpecSubUeMedianGen(0x0),
fHJetSpecSubUeConeGen(0x0),
fHJetUeMedianGen(0x0),
fHJetUeConeGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkSecOrFakeRec(0x0),
fHRhoUeMedianVsConeGen(0x0),
fIsChargedMC(0),
fIsFullMC(0),
faGenIndex(0),
faRecIndex(0),
fkAcceptance(2.0*TMath::Pi()*1.8),
fkDeltaPhiCut(TMath::Pi()-0.8),
fh1Xsec(0x0),
fh1Trials(0x0), 
fh1AvgTrials(0x0),
fh1PtHard(0x0),
fh1PtHardNoW(0x0),  
fh1PtHardTrials(0x0),
fAvgTrials(1),
fHardest(0),
fEventNumberRangeLow(0),
fEventNumberRangeHigh(99),
fTriggerPtRangeLow(0.0),
fTriggerPtRangeHigh(10000.0),
fFillRespMx(0),
fRandom(0x0),
fnTrials(1000),
fnPhi(9),
fnEta(2),
fEtaSize(0.7),
fPhiSize(2*TMath::Pi()/fnPhi),
fCellArea(fPhiSize*fEtaSize),
fSafetyMargin(1.1)
{
// Constructor

   DefineOutput(1, TList::Class());
}

//--------------------------------------------------------------
AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP(const AliAnalysisTaskJetCorePP& a):
AliAnalysisTaskSE(a.GetName()),
fESD(a.fESD),
fAODIn(a.fAODIn),
fAODOut(a.fAODOut),
fAODExtension(a.fAODExtension),
fJetBranchName(a.fJetBranchName),
fJetBranchNameChargMC(a.fJetBranchNameChargMC),
fJetBranchNameFullMC(a.fJetBranchNameFullMC),
fJetBranchNameBg(a.fJetBranchNameBg),
fJetBranchNameBgChargMC(a.fJetBranchNameBgChargMC),
fListJets(a.fListJets),
fListJetsGen(a.fListJetsGen),
fListJetsGenFull(a.fListJetsGenFull),
fListJetsBg(a.fListJetsBg),
fListJetsBgGen(a.fListJetsBgGen),
fNonStdFile(a.fNonStdFile),
fSystem(a.fSystem),  
fJetParamR(a.fJetParamR),
fBgJetParamR(a.fBgJetParamR),
fBgMaxJetPt(a.fBgMaxJetPt),
fBgConeR(a.fBgConeR),
fOfflineTrgMask(a.fOfflineTrgMask),
fMinContribVtx(a.fMinContribVtx),
fVtxZMin(a.fVtxZMin),
fVtxZMax(a.fVtxZMax),
fFilterMask(a.fFilterMask),
fCentMin(a.fCentMin),
fCentMax(a.fCentMax),
fJetEtaMin(a.fJetEtaMin),
fJetEtaMax(a.fJetEtaMax),
fTriggerEtaCut(a.fTriggerEtaCut),
fTrackEtaCut(a.fTrackEtaCut),
fTrackLowPtCut(a.fTrackLowPtCut),
fOutputList(a.fOutputList),
fHistEvtSelection(a.fHistEvtSelection),
fh2Ntriggers(a.fh2Ntriggers),
fHJetSpec(a.fHJetSpec),
fHJetSpecSubUeMedian(a.fHJetSpecSubUeMedian),
fHJetSpecSubUeCone(a.fHJetSpecSubUeCone),
fHJetUeMedian(a.fHJetUeMedian),
fHJetUeCone(a.fHJetUeCone),
fHRhoUeMedianVsCone(a.fHRhoUeMedianVsCone), 
//fHJetDensity(a.fHJetDensity),
//fHJetDensityA4(a.fHJetDensityA4),
fhJetPhi(a.fhJetPhi),
fhTriggerPhi(a.fhTriggerPhi),
fhJetEta(a.fhJetEta),
fhTriggerEta(a.fhTriggerEta),
fhVertexZ(a.fhVertexZ),
fhVertexZAccept(a.fhVertexZAccept),
fhContribVtx(a.fhContribVtx),
fhContribVtxAccept(a.fhContribVtxAccept),
fhDphiTriggerJet(a.fhDphiTriggerJet),
fhDphiTriggerJetAccept(a.fhDphiTriggerJetAccept),
fhCentrality(a.fhCentrality),
fhCentralityAccept(a.fhCentralityAccept),
//fHJetPtRaw(a.fHJetPtRaw),
//fHLeadingJetPtRaw(a.fHLeadingJetPtRaw),
//fHDphiVsJetPtAll(a.fHDphiVsJetPtAll),
fhJetPtGenVsJetPtRec(a.fhJetPtGenVsJetPtRec),
fhJetPtGenVsJetPtRecSubUeMedian(a.fhJetPtGenVsJetPtRecSubUeMedian),
fhJetPtGenVsJetPtRecSubUeCone(a.fhJetPtGenVsJetPtRecSubUeCone),
fhJetPtGen(a.fhJetPtGen),
fhJetPtSubUeMedianGen(a.fhJetPtSubUeMedianGen), 
fhJetPtSubUeConeGen(a.fhJetPtSubUeConeGen), 
fhJetPtGenChargVsJetPtGenFull(a.fhJetPtGenChargVsJetPtGenFull),
fhJetPtGenFull(a.fhJetPtGenFull),
fh2NtriggersGen(a.fh2NtriggersGen),
fHJetSpecGen(a.fHJetSpecGen),
fHJetSpecSubUeMedianGen(a.fHJetSpecSubUeMedianGen),
fHJetSpecSubUeConeGen(a.fHJetSpecSubUeConeGen),
fHJetUeMedianGen(a.fHJetUeMedianGen),
fHJetUeConeGen(a.fHJetUeConeGen),
fhPtTrkTruePrimRec(a.fhPtTrkTruePrimRec),
fhPtTrkTruePrimGen(a.fhPtTrkTruePrimGen),
fhPtTrkSecOrFakeRec(a.fhPtTrkSecOrFakeRec),
fHRhoUeMedianVsConeGen(a.fHRhoUeMedianVsConeGen),
fIsChargedMC(a.fIsChargedMC),
fIsFullMC(a.fIsFullMC),
faGenIndex(a.faGenIndex),
faRecIndex(a.faRecIndex),
fkAcceptance(a.fkAcceptance),
fkDeltaPhiCut(a.fkDeltaPhiCut),
fh1Xsec(a.fh1Xsec),
fh1Trials(a.fh1Trials),
fh1AvgTrials(a.fh1AvgTrials),
fh1PtHard(a.fh1PtHard),
fh1PtHardNoW(a.fh1PtHardNoW),  
fh1PtHardTrials(a.fh1PtHardTrials),
fAvgTrials(a.fAvgTrials),
fHardest(a.fHardest),
fEventNumberRangeLow(a.fEventNumberRangeLow),
fEventNumberRangeHigh(a.fEventNumberRangeHigh),
fTriggerPtRangeLow(a.fTriggerPtRangeLow),
fTriggerPtRangeHigh(a.fTriggerPtRangeHigh),
fFillRespMx(a.fFillRespMx),
fRandom(a.fRandom),
fnTrials(a.fnTrials),
fnPhi(a.fnPhi),
fnEta(a.fnEta),
fEtaSize(a.fEtaSize),
fPhiSize(a.fPhiSize),
fCellArea(a.fCellArea),
fSafetyMargin(a.fSafetyMargin)
{
   //Copy Constructor
   fRandom->SetSeed(0);
}
//--------------------------------------------------------------

AliAnalysisTaskJetCorePP& AliAnalysisTaskJetCorePP::operator = (const AliAnalysisTaskJetCorePP& a){
  // assignment operator
  this->~AliAnalysisTaskJetCorePP();
  new(this) AliAnalysisTaskJetCorePP(a);
  return *this;
}
//--------------------------------------------------------------

AliAnalysisTaskJetCorePP::~AliAnalysisTaskJetCorePP()
{
   //Destructor 
   delete fListJets;
   delete fListJetsGen;
   delete fListJetsGenFull;
   delete fListJetsBg;
   delete fListJetsBgGen;
   delete fOutputList; // ????
   delete fRandom;
}

//--------------------------------------------------------------


Bool_t AliAnalysisTaskJetCorePP::Notify()
{
   //Implemented Notify() to read the cross sections
   //and number of trials from pyxsec.root
   //inspired by AliAnalysisTaskJetSpectrum2::Notify()
   if(!fIsChargedMC) return kFALSE; 

   fESD = dynamic_cast<AliESDEvent*>(InputEvent());
   if(!fESD){
      if(fDebug>1) AliError("ESD not available");
      fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
   } 
 
   fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());


   if(fNonStdFile.Length()!=0){
      // case that we have an AOD extension we can fetch the jets from the extended output
      AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
      fAODExtension = aodH ? aodH->GetExtension(fNonStdFile.Data()) : 0;
      if(!fAODExtension){
         if(fDebug>1) Printf("AODExtension found for %s",fNonStdFile.Data());
      } 
   }
 
   TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
   Float_t xsection = 0;
   Float_t ftrials  = 1;

   fAvgTrials = 1;
   if(tree){
      TFile *curfile = tree->GetCurrentFile();
      if(!curfile) {
         Error("Notify","No current file");
         return kFALSE;
      }
      if(!fh1Xsec || !fh1Trials){
         Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
         return kFALSE;
      }
      AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
      fh1Xsec->Fill("<#sigma>",xsection);
      // construct a poor man average trials
      Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
      if(ftrials>=nEntries && nEntries>0.) fAvgTrials = ftrials/nEntries;
      fh1Trials->Fill("#sum{ntrials}",ftrials);
   }  

   if(fDebug)Printf("Reading File %s",fInputHandler->GetTree()->GetCurrentFile()->GetName());

   return kTRUE;
}
//--------------------------------------------------------------

void AliAnalysisTaskJetCorePP::Init()
{
   // check for jet branches
   if(!strlen(fJetBranchName.Data())){
      AliError("Jet branch name not set.");
   }
}

//--------------------------------------------------------------

void AliAnalysisTaskJetCorePP::UserCreateOutputObjects()
{
   // Create histograms and initilize variables
   fSafetyMargin = fBgConeR*fBgConeR /(fBgJetParamR*fBgJetParamR);
  
   // Called once
   fListJets   = new TList();  //reconstructed level antikt jets
   fListJetsBg = new TList();  //reconstructed jets to be removed from bg

   fIsChargedMC = (fJetBranchNameChargMC.Length()>0) ? kTRUE : kFALSE;
   fIsFullMC    = (fJetBranchNameFullMC.Length()>0)  ? kTRUE : kFALSE;

   fRandom = new TRandom3(0);

   if(fIsChargedMC){
      fListJetsGen   = new TList(); //generator level charged antikt jets
      fListJetsBgGen = new TList(); //generator level jets to be removed from bg 

      if(fIsFullMC)
         fListJetsGenFull = new TList(); //generator level full jets
   }
   OpenFile(1);
   if(!fOutputList) fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"event number (rejected)");
   
   fOutputList->Add(fHistEvtSelection);

   Int_t nBinsCentrality = (fSystem==0) ? 1 : 10; // pp=1 else 10
   //trigger pt spectrum (reconstructed) 
   fh2Ntriggers = new TH2F("fh2Ntriggers","# of triggers",
                             nBinsCentrality,0.0,100.0,50,0.0,50.0);
   fOutputList->Add(fh2Ntriggers);

   //Centrality, A, pTjet, pTtrigg, dphi
   const Int_t dimSpec   = 5;
   const Int_t nBinsSpec[dimSpec]     = {nBinsCentrality, 100,   220,  50, TMath::Nint(10*(TMath::Pi()-fkDeltaPhiCut))};
   const Double_t lowBinSpec[dimSpec] = {0.0,             0.0,   -20, 0.0, fkDeltaPhiCut};
   const Double_t hiBinSpec[dimSpec]  = {100.0,           1.0, 200.0,50.0, TMath::Pi()};
   fHJetSpec = new THnSparseF("fHJetSpec",
                   "Recoil jet spectrum [cent,A,pTjet,pTtrig,dphi]",
                   dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
   fOutputList->Add(fHJetSpec);  

   //background estimated as  median of kT jets 
   fHJetSpecSubUeMedian = (THnSparseF*) fHJetSpec->Clone("fHJetSpecSubUeMedian");
   fHJetSpecSubUeMedian->SetTitle("Recoil jet spectrum [cent,A,pTjet-pTUe,pTtrig,dphi]");
   fOutputList->Add(fHJetSpecSubUeMedian); 
   //background estimated as weighted  median of kT jets  ala Cone
   fHJetSpecSubUeCone = (THnSparseF*) fHJetSpec->Clone("fHJetSpecSubUeCone");
   fHJetSpecSubUeCone->SetTitle("Recoil jet spectrum [cent,A,pTjet-pTUe,pTtrig,dphi]");
   fOutputList->Add(fHJetSpecSubUeCone); 



   //------------------- HISTOS FOR DIAGNOSTIC ----------------------
   //A, pTjet, pTjet-pTUe, pTUe, rhoUe     bg estimate from kT median
   const Int_t    dimSpecMed   = 5;
   const Int_t    nBinsSpecMed[dimSpecMed]  = {50,     50,  120,     50,   50};
   const Double_t lowBinSpecMed[dimSpecMed] = {0.0,   0.0, -20.0,   0.0,  0.0};
   const Double_t hiBinSpecMed[dimSpecMed]  = {1.0, 100.0, 100.0,  20.0, 20.0};
   fHJetUeMedian = new THnSparseF("fHJetUeMedian",
                   "Recoil jet spectrum [A,pTjet,pTjet-pTUe, pTUe, rhoUe]",
                   dimSpecMed, nBinsSpecMed, lowBinSpecMed, hiBinSpecMed);
   fOutputList->Add(fHJetUeMedian);  
   
   //A, pTjet, pTjet-pTUe, pTUe, rhoUe     bg estimate from kT median Cone
   fHJetUeCone = (THnSparseF*) fHJetUeMedian->Clone("fHJetUeCone");
   fHJetUeCone->SetTitle("Recoil jet spectrum [A,pTjet,pTjet-pTUe, pTUe, rhoUe]");
   fOutputList->Add(fHJetUeCone); 

   //rho bacground reconstructed data
   const Int_t    dimRho   = 2;
   const Int_t    nBinsRho[dimRho]  = {200  , 200};
   const Double_t lowBinRho[dimRho] = {0.0  , 0.0};
   const Double_t hiBinRho[dimRho]  = {20.0 , 20.0};

   fHRhoUeMedianVsCone = new THnSparseF("hRhoUeMedianVsCone","[Rho Cone, Rho Median]",  
                                      dimRho, nBinsRho, lowBinRho, hiBinRho);
   fOutputList->Add(fHRhoUeMedianVsCone);

   //Jet number density histos [Trk Mult, jet density, pT trigger]
   /*const Int_t    dimJetDens   = 3;
   const Int_t    nBinsJetDens[dimJetDens]  = {100,   100,   10};
   const Double_t lowBinJetDens[dimJetDens] = {0.0,   0.0,  0.0};
   const Double_t hiBinJetDens[dimJetDens]  = {500.0, 5.0, 50.0};

   fHJetDensity = new THnSparseF("fHJetDensity","Jet dens vs trk mult A>0.07",
                                   dimJetDens,nBinsJetDens,lowBinJetDens,hiBinJetDens);

   fHJetDensityA4 =new THnSparseF("fHJetDensityA4","Jet dens vs trk mult A>0.4",
                                   dimJetDens,nBinsJetDens,lowBinJetDens,hiBinJetDens);

   fOutputList->Add(fHJetDensity);
   fOutputList->Add(fHJetDensityA4);
   */      

   //inclusive azimuthal and pseudorapidity histograms
   fhJetPhi = new TH2D("fhJetPhi","Azim dist jets vs pTjet",
                        50, 0, 100, 50,-TMath::Pi(),TMath::Pi());
   fhTriggerPhi= new TH2D("fhTriggerPhi","azim dist trig had vs pTtrigg",
                        50, 0, 50, 50,-TMath::Pi(),TMath::Pi());
   fhJetEta = new TH2D("fhJetEta","Eta dist jets vs pTjet",
                        50,0, 100, 40,-0.9,0.9);
   fhTriggerEta = new TH2D("fhTriggerEta","Eta dist trig had vs pTtrigg",
                        50, 0, 50, 40,-0.9,0.9);

   fhVertexZ = new TH1D("fhVertexZ","z vertex",40,-20,20);  
   fhVertexZAccept = new TH1D("fhVertexZAccept","z vertex after cut",40,-20,20);  
   fhContribVtx = new TH1D("fhContribVtx","contrib to vtx",200,0,200);   
   fhContribVtxAccept = new TH1D("fhContribVtxAccept","contrib to vtx after cut",200,0,200);   
   fhDphiTriggerJet = new TH1D("fhDphiTriggerJet","Deltaphi trig-jet",50, -TMath::Pi(),TMath::Pi()); 
   fhDphiTriggerJetAccept = new TH1D("fhDphiTriggerJetAccept","Deltaphi trig-jet after cut",50, -TMath::Pi(),TMath::Pi()); 
   fhCentrality = new TH1D("fhCentrality","Centrality",20,0,100);
   fhCentralityAccept = new TH1D("fhCentralityAccept","Centrality after cut",20,0,100);

   fOutputList->Add(fhJetPhi);
   fOutputList->Add(fhTriggerPhi);
   fOutputList->Add(fhJetEta);
   fOutputList->Add(fhTriggerEta);
   fOutputList->Add(fhVertexZ);    
   fOutputList->Add(fhVertexZAccept);    
   fOutputList->Add(fhContribVtx); 
   fOutputList->Add(fhContribVtxAccept); 
   fOutputList->Add(fhDphiTriggerJet);
   fOutputList->Add(fhDphiTriggerJetAccept);
   fOutputList->Add(fhCentrality); 
   fOutputList->Add(fhCentralityAccept);

   // raw spectra of INCLUSIVE jets  
   //Centrality, pTjet, A
   /*const Int_t dimRaw   = 3;
   const Int_t nBinsRaw[dimRaw]     = {nBinsCentrality,  50,   100};
   const Double_t lowBinRaw[dimRaw] = {0.0,             0.0,   0.0};
   const Double_t hiBinRaw[dimRaw]  = {100.0,           100,   1.0};
   fHJetPtRaw = new THnSparseF("fHJetPtRaw",
                                "Incl. jet spectrum [cent,pTjet,A]",
                                dimRaw,nBinsRaw,lowBinRaw,hiBinRaw);
   fOutputList->Add(fHJetPtRaw);  

   // raw spectra of LEADING jets  
   //Centrality, pTjet, A
   fHLeadingJetPtRaw = new THnSparseF("fHLeadingJetPtRaw",
                                "Leading jet spectrum [cent,pTjet,A]",
                                dimRaw,nBinsRaw,lowBinRaw,hiBinRaw);
   fOutputList->Add(fHLeadingJetPtRaw);  

   // Dphi versus pT jet 
   //Centrality, Dphi=phiTrig-phiJet, pTjet, pTtrigg 
   const Int_t dimDp   = 4;
   const Int_t nBinsDp[dimDp]     = {nBinsCentrality,  50,     50,    50};
   const Double_t lowBinDp[dimDp] = {0.0,       -TMath::Pi(),   0.0,   0.0};
   const Double_t hiBinDp[dimDp]  = {100.0,      TMath::Pi(), 100.0, 100.0};
   fHDphiVsJetPtAll = new THnSparseF("fHDphiVsJetPtAll",
                                "Dphi vs jet pT [cent,Dphi,pTjet,pTtrigg]",
                                dimDp,nBinsDp,lowBinDp,hiBinDp);
   fOutputList->Add(fHDphiVsJetPtAll);  
   */

   //analyze MC generator level 
   if(fIsChargedMC){   
      if(fFillRespMx){
         //Fill response matrix only once 
         fhJetPtGenVsJetPtRec = new TH2D("fhJetPtGenVsJetPtRec","JetPtGenVsJetPtRec", 200,0,200, 200,0,200); 
         fOutputList->Add(fhJetPtGenVsJetPtRec); //gen MC charg jet pt spectrum versus rec charged jet pt spectrum
         //....
         fhJetPtGenVsJetPtRecSubUeMedian = new TH2D("fhJetPtGenVsJetPtRecSubUeMedian","fhJetPtGenVsJetPtRecSubUeMedian", 220,-20,200, 220,-20,200); 
         fOutputList->Add(fhJetPtGenVsJetPtRecSubUeMedian); // with kT median bg subtr 
         //....
         fhJetPtGenVsJetPtRecSubUeCone=(TH2D*)fhJetPtGenVsJetPtRecSubUeMedian->Clone("fhJetPtGenVsJetPtRecSubUeCone");
         fhJetPtGenVsJetPtRecSubUeCone->SetTitle("fhJetPtGenVsJetPtRecSubUeCone");
         fOutputList->Add(fhJetPtGenVsJetPtRecSubUeCone); // with weighted kT median bg subtr 
         //....
         fhJetPtGen = new TH1D("fhJetPtGen","Jet Pt (MC Gen)",200,0,200); //MC generator charged jet pt spectrum
         fOutputList->Add(fhJetPtGen);  
         //....
         fhJetPtSubUeMedianGen = new TH1D("fhJetPtSubUeMedianGen","Jet Pt - UE pT (MC Gen)",220,-20,200); 
         fOutputList->Add(fhJetPtSubUeMedianGen);  // with kT median bg subtr
         //....
         fhJetPtSubUeConeGen = (TH1D*) fhJetPtSubUeMedianGen->Clone("fhJetPtSubUeConeGen");
         fOutputList->Add(fhJetPtSubUeConeGen); // with weighted kT median bg subtr
      
         //....
         if(fIsFullMC){
            fhJetPtGenChargVsJetPtGenFull = new TH2D("fhJetPtGenChargVsJetPtGenFull","fhJetPtGenChargVsJetPtGenFull", 200,0,200, 200,0,200);
            fOutputList->Add(fhJetPtGenChargVsJetPtGenFull); //gen full MC jet pt versus gen charged jet MC pt
         //....
            fhJetPtGenFull = new TH1D("fhJetPtGenFull","Jet Pt (MC Full jets Gen)",200,0,200); //MC generator full jet pt spectrum
            fOutputList->Add(fhJetPtGenFull); 
         }
      }
      //....
      fh2NtriggersGen = (TH2F*) fh2Ntriggers->Clone("fh2NtriggersGen");
      fh2NtriggersGen->SetTitle(Form("%s Gen MC",fh2Ntriggers->GetTitle()));
      fOutputList->Add(fh2NtriggersGen);

      //Centrality, A, pT, pTtrigg
      fHJetSpecGen = (THnSparseF*) fHJetSpec->Clone("fHJetSpecGen");
      fHJetSpecGen->SetTitle(Form("%s Gen MC",fHJetSpec->GetTitle()));
      fOutputList->Add(fHJetSpecGen); 

      fHJetSpecSubUeMedianGen = (THnSparseF*)  fHJetSpecSubUeMedian->Clone("fHJetSpecSubUeMedianGen");
      fHJetSpecSubUeMedianGen->SetTitle(Form("%s Gen MC",fHJetSpecSubUeMedian->GetTitle()));
      fOutputList->Add(fHJetSpecSubUeMedianGen);  

      fHJetSpecSubUeConeGen =  (THnSparseF*) fHJetSpecSubUeCone->Clone("fHJetSpecSubUeConeGen"); 
      fHJetSpecSubUeConeGen->SetTitle(Form("%s Gen MC",fHJetSpecSubUeCone->GetTitle()));
      fOutputList->Add(fHJetSpecSubUeConeGen);
      //---
      fHJetUeMedianGen  =  (THnSparseF*) fHJetUeMedian->Clone("fHJetUeMedianGen");  
      fHJetUeMedianGen->SetTitle(Form("%s Gen MC", fHJetUeMedian->GetTitle()));
      fOutputList->Add(fHJetUeMedianGen);

      fHJetUeConeGen =  (THnSparseF*) fHJetUeCone->Clone("fHJetUeConeGen"); 
      fHJetUeConeGen->SetTitle(Form("%s Gen MC", fHJetUeCone->GetTitle())); 
      fOutputList->Add(fHJetUeConeGen);

      if(fFillRespMx){
         //track efficiency/contamination histograms  eta versus pT
         fhPtTrkTruePrimRec = new TH2D("fhPtTrkTruePrimRec","PtTrkTruePrimRec",100,0,20,18,-0.9,0.9);
         fOutputList->Add(fhPtTrkTruePrimRec); 
      
         fhPtTrkTruePrimGen = (TH2D*) fhPtTrkTruePrimRec->Clone("fhPtTrkTruePrimGen");
         fhPtTrkTruePrimGen->SetTitle("PtTrkTruePrimGen");    
         fOutputList->Add(fhPtTrkTruePrimGen);
      
         fhPtTrkSecOrFakeRec = (TH2D*) fhPtTrkTruePrimRec->Clone("fhPtTrkSecOrFakeRec");    
         fhPtTrkSecOrFakeRec->SetTitle("PtTrkSecOrFakeRec");    
         fOutputList->Add(fhPtTrkSecOrFakeRec);
      }

      fHRhoUeMedianVsConeGen = (THnSparseF*) fHRhoUeMedianVsCone->Clone("hRhoUeMedianVsConeGen");
      fHRhoUeMedianVsConeGen->SetTitle(Form("%s Gen MC", fHRhoUeMedianVsCone->GetTitle())); 
      fOutputList->Add(fHRhoUeMedianVsConeGen);
   }
   //-------------------------------------
   //     pythia histograms
   const Int_t nBinPt = 150;
   Double_t binLimitsPt[nBinPt+1];
   for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
      if(iPt == 0){
         binLimitsPt[iPt] = -50.0;
      }else{// 1.0
         binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.;
      }
   }
   
   fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
   fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
   fOutputList->Add(fh1Xsec);
   fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
   fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
   fOutputList->Add(fh1Trials);
   fh1AvgTrials = new TH1F("fh1AvgTrials","trials root file",1,0,1);
   fh1AvgTrials->GetXaxis()->SetBinLabel(1,"#sum{avg ntrials}");
   fOutputList->Add(fh1AvgTrials);
   fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
   fOutputList->Add(fh1PtHard);
   fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
   fOutputList->Add(fh1PtHardNoW);
   fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);
   fOutputList->Add(fh1PtHardTrials);      
   

   // =========== Switch on Sumw2 for all histos ===========
   for(Int_t i=0; i<fOutputList->GetEntries(); i++){
      TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
      if(h1){
         h1->Sumw2();
         continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fOutputList->At(i));
      if(hn){
         hn->Sumw2();
      }	  
   }
   TH1::AddDirectory(oldStatus);

   PostData(1, fOutputList);
}

//--------------------------------------------------------------------

void AliAnalysisTaskJetCorePP::UserExec(Option_t *)
{
   //Event loop
   Double_t eventW  = 1.0;
   Double_t ptHard  = 0.0;
   Double_t nTrials = 1.0; // Trials for MC trigger
   if(fIsChargedMC) fh1AvgTrials->Fill("#sum{avg ntrials}",fAvgTrials); 

   if(TMath::Abs((Float_t) fJetParamR)<0.00001){
      AliError("ANTIKT Cone radius is set to zero.");  
      return;
   }
   if(TMath::Abs((Float_t) fBgJetParamR)<0.00001){
      AliError("KT Cone radius is set to zero.");  
      return;
   }
   if(!strlen(fJetBranchName.Data())){
      AliError("Jet branch name not set.");
      return;
   }
   if(!strlen(fJetBranchNameBg.Data())){
      AliError("Jet bg branch name not set.");
      return;
   }


   fESD = dynamic_cast<AliESDEvent*>(InputEvent());
   if(!fESD){
      if(fDebug>1) AliError("ESD not available");
      fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
   } 
 
   fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());
   AliAODEvent* aod = NULL;
   // take all other information from the aod we take the tracks from
   if(!aod){
      if(!fESD) aod = fAODIn;
      else      aod = fAODOut;
   }

   if(fNonStdFile.Length()!=0){
      // case that we have an AOD extension we can fetch the jets from the extended output
      AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
      fAODExtension = aodH ? aodH->GetExtension(fNonStdFile.Data()) : 0;
      if(!fAODExtension){
         if(fDebug>1) Printf("AODExtension found for %s",fNonStdFile.Data());
      } 
   }
    
   // ----------------- EVENT SELECTION  --------------------------
   fHistEvtSelection->Fill(1); // number of events before event selection

   // physics selection
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)
           ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());

   if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
      if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
      fHistEvtSelection->Fill(2);
      PostData(1, fOutputList);
      return;
   }
  
   //check AOD pointer
   if(!aod){
      if(fDebug) Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
      fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
      return;
   }
   
   // vertex selection
   AliAODVertex* primVtx = aod->GetPrimaryVertex();

   if(!primVtx){
      if(fDebug) Printf("%s:%d No primVtx",(char*)__FILE__,__LINE__);
      fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
      return;
   }

   Int_t nTracksPrim = primVtx->GetNContributors();
   Float_t vtxz = primVtx->GetZ();
   //Input events
   fhContribVtx->Fill(nTracksPrim);
   if( nTracksPrim > 0 ) fhVertexZ->Fill(vtxz);

   if((nTracksPrim < fMinContribVtx) ||
      (vtxz < fVtxZMin) ||
      (vtxz > fVtxZMax)){
      if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",
                         (char*)__FILE__,__LINE__,vtxz);
      fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
      return;
   }else{
      //Accepted events
      fhContribVtxAccept->Fill(nTracksPrim);
      fhVertexZAccept->Fill(vtxz);
   }
   //FK// No event class selection imposed
   // event class selection (from jet helper task)
   //Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
   //if(fDebug) Printf("Event class %d", eventClass);
   //if(eventClass < fEvtClassMin || eventClass > fEvtClassMax){
   //   fHistEvtSelection->Fill(4);
   //   PostData(1, fOutputList);
   //   return;
   //}

   //------------------ CENTRALITY SELECTION ---------------
   AliCentrality *cent = 0x0;
   Double_t centValue  = 0.0; 
   if(fSystem){  //fSystem=0 for pp,   fSystem=1 for pPb
      if(fESD){
         cent = fESD->GetCentrality();
         if(cent) centValue = cent->GetCentralityPercentile("V0M");
      }else{
         centValue = aod->GetHeader()->GetCentrality();
      }   
      if(fDebug) printf("centrality: %f\n", centValue);
      //Input events
      fhCentrality->Fill(centValue); 

      if(centValue < fCentMin || centValue > fCentMax){
         fHistEvtSelection->Fill(4);
         PostData(1, fOutputList);
         return;
      }else{
         //Accepted events
         fhCentralityAccept->Fill( centValue );
      }
   }
 
   //-----------------select disjunct event subsamples ----------------
   Int_t eventnum  = aod->GetHeader()->GetEventNumberESDFile();
   Int_t lastdigit = eventnum % 10;
   if(!(fEventNumberRangeLow<=lastdigit && lastdigit<=fEventNumberRangeHigh)){
      fHistEvtSelection->Fill(5);
      PostData(1, fOutputList);
      return;
   } 

   if(fDebug) std::cout<<" ACCEPTED EVENT "<<endl;
   fHistEvtSelection->Fill(0); // accepted events 
   // ==================== end event selection ============================ 
 
   Double_t tmpArrayFive[5];
 
   //=============== Reconstructed antikt jets =============== 
   ReadTClonesArray(fJetBranchName.Data()  , fListJets); 
   ReadTClonesArray(fJetBranchNameBg.Data(), fListJetsBg); 

   //============ Estimate background in reconstructed events ===========
   Double_t rhoFromCellMedian=0.0,    rhoCone=0.0;

   //Find Hadron trigger and Estimate rho from cone
   TList particleList; //list of tracks
   Int_t indexTrigg = GetListOfTracks(&particleList); //index of trigger hadron in Particle list
   EstimateBgCone(fListJets, &particleList, rhoCone);

   //Estimate rho from cell median minus jets
   EstimateBgRhoMedian(fListJetsBg, &particleList, rhoFromCellMedian);
   //fListJetsBg->Clear(); //this list is further not needed

   //==============  analyze generator level MC  ================ 
   TList particleListGen; //list of tracks in MC

   if(fIsChargedMC){       
      //================= generated charged antikt jets ================
      ReadTClonesArray(fJetBranchNameChargMC.Data(),   fListJetsGen); 
      ReadTClonesArray(fJetBranchNameBgChargMC.Data(), fListJetsBgGen); 

      if(fIsFullMC){ //generated full jets
         ReadTClonesArray(fJetBranchNameFullMC.Data(),  fListJetsGenFull); 
      }

      //========================================================
      //serarch for charged trigger at the MC generator level

      Int_t    indexTriggGen = -1;
      Double_t ptTriggGen    = -1;
      Int_t    iCounterGen   =  0; //number of entries in particleListGen array
      Int_t    triggersMC[200];//list of trigger candidates
      Int_t    ntriggersMC   = 0; //index in triggers array

      if(fESD){//ESD input

         AliMCEvent* mcEvent = MCEvent();
         if(!mcEvent){
            PostData(1, fOutputList);
            return;
         }
         
         AliGenPythiaEventHeader *pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
         if(pythiaGenHeader){
            nTrials = pythiaGenHeader->Trials();
            ptHard  = pythiaGenHeader->GetPtHard();
            fh1PtHard->Fill(ptHard,eventW);
            fh1PtHardNoW->Fill(ptHard,1);
            fh1PtHardTrials->Fill(ptHard,nTrials);
         }

         for(Int_t it = 0; it < mcEvent->GetNumberOfTracks(); it++){
            if(!mcEvent->IsPhysicalPrimary(it)) continue;
            AliMCParticle* part = (AliMCParticle*) mcEvent->GetTrack(it);
            if(!part) continue;  
            if(SelectMCGenTracks((AliVParticle*) part, &particleListGen, ptTriggGen, indexTriggGen, iCounterGen)){ 

               if(fHardest==0 && ntriggersMC<200){//single inclusive trigger
                  if(indexTriggGen > -1){//trigger candidate was found
                     triggersMC[ntriggersMC] = indexTriggGen;
                     ntriggersMC++; 
                  }
               }

               iCounterGen++;//index in particleListGen array 
            }
         }
      }else{  //AOD input
         if(!fAODIn){
            PostData(1, fOutputList);
            return;
         }
         TClonesArray *tca = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(AliAODMCParticle::StdBranchName()));
         if(!tca){ 
            PostData(1, fOutputList);
            return;
         }
         for(Int_t it = 0; it < tca->GetEntriesFast(); it++){
            AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(it));
            if(!part) continue;  
            if(!part->IsPhysicalPrimary()) continue;
            if(SelectMCGenTracks((AliVParticle*) part, &particleListGen, ptTriggGen, indexTriggGen, iCounterGen)){

               if(fHardest==0 && ntriggersMC<200){//single inclusive trigger
                  if(indexTriggGen > -1){ //trigger candidater was found
                     triggersMC[ntriggersMC] = indexTriggGen;
                     ntriggersMC++; 
                  }
               }

               iCounterGen++;//number of entries in particleListGen array
            }
         }
      }

      //==============  Estimate bg in generated events ==============
      Double_t rhoFromCellMedianGen=0.0, rhoConeGen=0.0;

      //Estimate rho from cone
      EstimateBgCone(fListJetsGen, &particleListGen, rhoConeGen);

      //Estimate rho from cell median minus jets
      EstimateBgRhoMedian(fListJetsBgGen, &particleListGen, rhoFromCellMedianGen);
      //fListJetsBgGen->Clear();

      //============  Generator trigger+jet ==================
      if(fHardest==0){ //single inclusive trigger
         if(ntriggersMC>0){ //there is at least one trigger 
            Int_t rnd     = fRandom->Integer(ntriggersMC); //0 to ntriggers-1
            indexTriggGen = triggersMC[rnd];
         }else{
            indexTriggGen = -1; //trigger not found
         }
      }

      //-----------
      Int_t ilowGen  = (fHardest==0 || fHardest==1) ? indexTriggGen : 0;
      Int_t ihighGen = (fHardest==0 || fHardest==1) ? indexTriggGen+1 : particleListGen.GetEntries();
      Bool_t fillOnceGen = kTRUE;
      //-----------
      
      for(Int_t igen= ilowGen; igen< ihighGen; igen++){ //loop over possible trigger
         indexTriggGen = igen; //trigger hadron 

         if(indexTriggGen == -1) continue;  
         AliVParticle* triggerGen=(AliVParticle*) particleListGen.At(indexTriggGen);  
         if(!triggerGen) continue;
 
         if(fHardest >= 2){ 
            if(triggerGen->Pt() < 6.0) continue; //all hadrons pt>6  
         }
         if(TMath::Abs((Float_t) triggerGen->Eta()) > fTriggerEtaCut) continue;

         ptTriggGen = triggerGen->Pt(); //count triggers
         fh2NtriggersGen->Fill(centValue, ptTriggGen);

         //Count jets and trigger-jet pairs at MC  generator level
         for(Int_t ij=0; ij<fListJetsGen->GetEntries(); ij++){
            AliAODJet* jet = (AliAODJet*)(fListJetsGen->At(ij));
            if(!jet) continue;
            Double_t etaJetGen = jet->Eta();
            Double_t areaJetGen = jet->EffectiveAreaCharged();
 
            if((fJetEtaMin<=etaJetGen) && (etaJetGen<=fJetEtaMax)){ 

               Double_t dphi = RelativePhi(triggerGen->Phi(), jet->Phi()); 
               if(TMath::Abs((Double_t) dphi) < fkDeltaPhiCut) continue;

               //Centrality, A, pT, pTtrigg
               tmpArrayFive[0] = centValue;
               tmpArrayFive[1] = areaJetGen;
               tmpArrayFive[2] = jet->Pt();
               tmpArrayFive[3] = ptTriggGen;
               tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
               fHJetSpecGen->Fill(tmpArrayFive);


               Double_t ptUeFromCellMedianGen = rhoFromCellMedianGen*areaJetGen;
               Double_t ptUeConeGen         = rhoConeGen*areaJetGen;

               //Centrality, A, pTjet-pTbgCellMedian, pTtrigg, dphi
               tmpArrayFive[0] = centValue;
               tmpArrayFive[1] = areaJetGen;
               tmpArrayFive[2] = jet->Pt() - ptUeFromCellMedianGen;
               tmpArrayFive[3] = ptTriggGen;
               tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
               fHJetSpecSubUeMedianGen->Fill(tmpArrayFive);

               //Centrality, A, pTjet-pTbgCone, pTtrigg, dphi
               tmpArrayFive[0] = centValue;
               tmpArrayFive[1] = areaJetGen;
               tmpArrayFive[2] = jet->Pt() - ptUeConeGen;
               tmpArrayFive[3] = ptTriggGen;
               tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
               fHJetSpecSubUeConeGen->Fill(tmpArrayFive);

               //Ue diagnostics  "[A,pTjet,pTjet-pTUe, pTUe, rhoUe]", kT median
               tmpArrayFive[0] = areaJetGen;
               tmpArrayFive[1] = jet->Pt();
               tmpArrayFive[2] = jet->Pt() - ptUeFromCellMedianGen;
               tmpArrayFive[3] = ptUeFromCellMedianGen;
               tmpArrayFive[4] = rhoFromCellMedianGen;
               fHJetUeMedianGen->Fill(tmpArrayFive); 

               //Ue diagnostics  "[A,pTjet,pTjet-pTUe, pTUe, rhoUe]", perpendicular Cone
               tmpArrayFive[0] = areaJetGen;
               tmpArrayFive[1] = jet->Pt();
               tmpArrayFive[2] = jet->Pt() - ptUeConeGen;
               tmpArrayFive[3] = ptUeConeGen;
               tmpArrayFive[4] = rhoConeGen;
               fHJetUeConeGen->Fill(tmpArrayFive);

               if(fillOnceGen){
                  Double_t fillRhoGen[] = { rhoConeGen, rhoFromCellMedianGen};
                  fHRhoUeMedianVsConeGen->Fill(fillRhoGen); 
                  fillOnceGen = kFALSE;
               }
            }//back to back jet-trigger pair
         }//jet loop
      }//trigger loop

      if(fFillRespMx){
         //================ RESPONSE MATRIX ===============
         //Count jets and trigger-jet pairs at MC  generator level
         for(Int_t ij=0; ij<fListJetsGen->GetEntries(); ij++){
            AliAODJet* jet = (AliAODJet*)(fListJetsGen->At(ij));
            if(!jet) continue;
            Double_t etaJetGen = jet->Eta();
            Double_t ptJetGen  = jet->Pt();
       
            if((fJetEtaMin<=etaJetGen) && (etaJetGen<=fJetEtaMax)){ 
               fhJetPtGen->Fill(ptJetGen); // gen level pt spectum of jets response mx normalization

               Double_t areaJetGen = jet->EffectiveAreaCharged();
               Double_t ptUeFromCellMedianGen = rhoFromCellMedianGen*areaJetGen;
               Double_t ptUeConeGen         = rhoConeGen*areaJetGen;

               fhJetPtSubUeMedianGen->Fill(ptJetGen - ptUeFromCellMedianGen); 
               fhJetPtSubUeConeGen->Fill(ptJetGen   - ptUeConeGen);    
            }
         }
         if(fListJets->GetEntries()>0 && fListJetsGen->GetEntries()>0){ //at least some reconstructed jets
         
            Int_t ng = (Int_t) fListJetsGen->GetEntries();
            Int_t nr = (Int_t) fListJets->GetEntries();

            //Find closest MC generator - reconstructed jets
            if(faGenIndex.GetSize()<nr) faGenIndex.Set(nr); //idx of gen jet assoc to rec jet
            if(faRecIndex.GetSize()<ng) faRecIndex.Set(ng); //idx of rec jet assoc to gen jet

            if(fDebug){
               Printf("New Rec List %d gen index Array %d",nr,faGenIndex.GetSize());
               Printf("New Gen List %d rec index Array %d",ng,faRecIndex.GetSize());
            }
            //matching of MC genrator level and reconstructed jets
            AliAnalysisHelperJetTasks::GetClosestJets(fListJetsGen,ng,fListJets,nr,faGenIndex,faRecIndex,fDebug); 

            // Fill response matrix
            for(Int_t ir = 0; ir < nr; ir++){
               AliAODJet *recJet   = (AliAODJet*) fListJets->At(ir);
               Double_t etaJetRec  = recJet->Eta();
               Double_t ptJetRec   = recJet->Pt();
               Double_t areaJetRec = recJet->EffectiveAreaCharged();
               //fill response matrix if generator and reconstructed jets are within |eta|<0.9-fiduc

               if((fJetEtaMin <= etaJetRec) && (etaJetRec <= fJetEtaMax)){ 
                  Int_t ig = faGenIndex[ir]; //associated generator level jet
                  if(ig >= 0 && ig < ng){
                     if(fDebug > 10) Printf("%s:%d ig = %d ir = %d",(char*)__FILE__,__LINE__,ig,ir);
                     AliAODJet *genJet  = (AliAODJet*) fListJetsGen->At(ig);
                     Double_t ptJetGen  = genJet->Pt();
                     Double_t etaJetGen = genJet->Eta();

                     //fill response matrix if generator and reconstructed jets are within |eta|<0.9-fiduc
                     if((fJetEtaMin <= etaJetGen) && (etaJetGen <= fJetEtaMax)){
                        fhJetPtGenVsJetPtRec->Fill(ptJetRec, ptJetGen);

                        Double_t areaJetGen  = genJet->EffectiveAreaCharged();
                        Double_t ptUeFromCellMedianGen = rhoFromCellMedianGen*areaJetGen;
                        Double_t ptUeConeGen         = rhoConeGen*areaJetGen;
                        Double_t ptUeFromCellMedianRec = rhoFromCellMedian*areaJetRec;
                        Double_t ptUeConeRec         = rhoCone*areaJetRec;
                        fhJetPtGenVsJetPtRecSubUeMedian->Fill(ptJetRec-ptUeFromCellMedianRec, 
                                                           ptJetGen-ptUeFromCellMedianGen);
                        fhJetPtGenVsJetPtRecSubUeCone->Fill(ptJetRec-ptUeConeRec, ptJetGen-ptUeConeGen);
                     }
                  }//ig>=0
               }//rec jet in eta acceptance
            }//loop over reconstructed jets
         }// # of  rec jets >0
      
         //=========================== Full jet vs charged jet matrix  ==========================
         if(fIsFullMC){
            //Count full jets and charged-jet pairs at MC  generator level
            for(Int_t ij=0; ij<fListJetsGenFull->GetEntries(); ij++){
               AliAODJet* jetFull = (AliAODJet*)(fListJetsGenFull->At(ij));
               if(!jetFull) continue;
               Double_t etaJetFull = jetFull->Eta();
               Double_t ptJetFull  = jetFull->Pt();

               if((fJetEtaMin<=etaJetFull) && (etaJetFull<=fJetEtaMax)){
                  fhJetPtGenFull->Fill(ptJetFull); // generator level pt spectum of full jets 
               }
            }
            if(fListJetsGen->GetEntries()>0 && fListJetsGenFull->GetEntries()>0){ //at least some reconstructed jets
               Int_t nful = (Int_t) fListJetsGenFull->GetEntries();
               Int_t nchr = (Int_t) fListJetsGen->GetEntries();
 
               //Find closest MC generator full - charged jet
               if(faGenIndex.GetSize()<nchr) faGenIndex.Set(nchr); //idx of gen FULL jet assoc to gen CHARGED jet
               if(faRecIndex.GetSize()<nful) faRecIndex.Set(nful); //idx of gen CHARGED jet assoc to gen FULL jet
 
               if(fDebug){
                  Printf("New Charg List %d Full index Array %d",nchr,faGenIndex.GetSize());
                  Printf("New Full List %d Charg index Array %d",nful,faRecIndex.GetSize());
               }
               //matching of MC genrator level and reconstructed jets
               AliAnalysisHelperJetTasks::GetClosestJets(fListJetsGenFull,nful,fListJetsGen,nchr,faGenIndex,faRecIndex,fDebug);

              // Fill response matrix
               for(Int_t ichr = 0; ichr < nchr; ichr++){ //charged jet loop
                  AliAODJet *chJet  = (AliAODJet*) fListJetsGen->At(ichr);
                  Double_t etaJetCh = chJet->Eta();
                  Double_t ptJetCh  = chJet->Pt();
                  //fill response matrix if generator and reconstructed jets are within |eta|<0.9-fiduc
    
                  if((fJetEtaMin <= etaJetCh) && (etaJetCh <= fJetEtaMax)){
                     Int_t iful = faGenIndex[ichr]; //associated generator level jet
                     if(iful >= 0 && iful < nful){
                        if(fDebug > 10) Printf("%s:%d iful = %d ichr = %d",(char*)__FILE__,__LINE__,iful,ichr);
                        AliAODJet *genJetFull  = (AliAODJet*) fListJetsGenFull->At(iful);
                        Double_t ptJetFull  = genJetFull->Pt();
                        Double_t etaJetFull = genJetFull->Eta();

                        //fill response matrix if generator and reconstructed jets are within |eta|<0.9-fiduc
                        if((fJetEtaMin <= etaJetFull) && (etaJetFull <= fJetEtaMax)){
                           fhJetPtGenChargVsJetPtGenFull->Fill(ptJetFull,ptJetCh);
                        }
                     }//iful>=0
                  }//rec jet in eta acceptance
               }//loop over reconstructed jets
            }// # of  rec jets >0
         }//pointer MC generator jets
      } //fill resp mx only for bin 
   }//analyze generator level MC
    
    
    
   //=============  RECONSTRUCTED INCLUSIVE JETS ===============

   Double_t etaJet  = 0.0;
   Double_t pTJet   = 0.0;
   Double_t areaJet = 0.0;
   Double_t phiJet  = 0.0;
   //Int_t indexLeadingJet     = -1;
   //Double_t pTLeadingJet     = -10.0; 
   //Double_t areaLeadingJet   = -10.0;
  
   for(Int_t ij=0; ij<fListJets->GetEntries(); ij++){
      AliAODJet* jet = (AliAODJet*)(fListJets->At(ij));
      if(!jet) continue;
      etaJet  = jet->Eta();
      phiJet  = jet->Phi();
      pTJet   = jet->Pt();
      if(pTJet==0) continue; 
     
      if((etaJet<fJetEtaMin) || (etaJet>fJetEtaMax)) continue;
      /*areaJet = jet->EffectiveAreaCharged();*/

      //Jet Diagnostics---------------------------------
      fhJetPhi->Fill(pTJet, RelativePhi(phiJet,0.0)); //phi -pi,pi
      fhJetEta->Fill(pTJet, etaJet);
      //search for leading jet
      /*if(pTJet > pTLeadingJet){
         indexLeadingJet  = ij; 
         pTLeadingJet     = pTJet; 
         areaLeadingJet   = areaJet; 
      } 
 
      // raw spectra of INCLUSIVE jets  
      //Centrality, pTjet, A
      Double_t fillraw[] = { centValue,
                             pTJet,
                             areaJet
                           };
      fHJetPtRaw->Fill(fillraw);*/
   }
   /*
   if(indexLeadingJet > -1){ 
      // raw spectra of LEADING jets  
      //Centrality, pTjet,  A
      Double_t fillleading[] = { centValue,
                                 pTLeadingJet,
                                 areaLeadingJet
                               };
      fHLeadingJetPtRaw->Fill(fillleading);
   } 
   */
 
   // ===============  RECONSTRUCTED TRIGGER-JET PAIRS ================
   if(fIsChargedMC && fFillRespMx){
     FillEffHistos(&particleList, &particleListGen); //Fill efficiency histos
   }
   Bool_t filledOnce = kTRUE; //fill rho histogram only once per event
   //set ranges of the trigger loop
   Int_t ilow  = (fHardest==0 || fHardest==1) ? indexTrigg : 0;
   Int_t ihigh = (fHardest==0 || fHardest==1) ? indexTrigg+1 : particleList.GetEntries();

   for(Int_t itrk= ilow; itrk< ihigh; itrk++){ //loop over possible trigger
      indexTrigg = itrk; //trigger hadron with pT >10 GeV 
 
      if(indexTrigg < 0) continue;

      AliVParticle *triggerHadron = (AliVParticle*) particleList.At(indexTrigg);     
      if(!triggerHadron){  
         PostData(1, fOutputList);
         return;
      }
      if(fHardest >= 2){ 
         if(triggerHadron->Pt() < 6.0) continue; //all hadrons pt>6  
      }
      if(TMath::Abs((Float_t) triggerHadron->Eta())> fTriggerEtaCut) continue;
 
      //Fill trigger histograms
      fh2Ntriggers->Fill(centValue,triggerHadron->Pt()); //trigger pT
      fhTriggerPhi->Fill(triggerHadron->Pt(),RelativePhi(triggerHadron->Phi(),0.0)); //phi -pi,pi
      fhTriggerEta->Fill(triggerHadron->Pt(),triggerHadron->Eta());

  
      //---------- make trigger-jet pairs ---------
      //Int_t injet4     = 0;
      //Int_t injet      = 0; 

      for(Int_t ij=0; ij<fListJets->GetEntries(); ij++){
         AliAODJet* jet = (AliAODJet*)(fListJets->At(ij));
         if(!jet) continue;
         etaJet  = jet->Eta();
         phiJet  = jet->Phi();
         pTJet   = jet->Pt();
         if(pTJet==0) continue; 
     
         if((etaJet<fJetEtaMin) || (etaJet>fJetEtaMax)) continue;
         areaJet = jet->EffectiveAreaCharged();
         //if(areaJet >= 0.07) injet++; 
         //if(areaJet >= 0.4)  injet4++;

         Double_t dphi = RelativePhi(triggerHadron->Phi(), phiJet); 
         fhDphiTriggerJet->Fill(dphi); //Input

         //Dphi versus jet pT   
         //Centrality, Dphi=phiTrig-phiJet, pTjet, pTtrigg 
         /*Double_t filldp[] = { centValue,
                               dphi,
                               pTJet,
                               triggerHadron->Pt()
                             };
         fHDphiVsJetPtAll->Fill(filldp);
         */ 
         // Select back to back trigger - jet pairs
         if(TMath::Abs((Double_t) dphi) < fkDeltaPhiCut) continue;
         fhDphiTriggerJetAccept->Fill(dphi); //Accepted
          
         //NO bg subtraction
         //Centrality, A, pTjet, pTtrigg, dphi
         tmpArrayFive[0] = centValue;
         tmpArrayFive[1] = areaJet;
         tmpArrayFive[2] = pTJet;
         tmpArrayFive[3] = triggerHadron->Pt();
         tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
         fHJetSpec->Fill(tmpArrayFive);

         //subtract bg using estinates base on median of kT jets
         Double_t ptUeFromCellMedian = rhoFromCellMedian*areaJet;
         Double_t ptUeCone         = rhoCone*areaJet;

         //Centrality, A, pTjet-pTbgCellMedian, pTtrigg, dphi
         tmpArrayFive[0] = centValue;
         tmpArrayFive[1] = areaJet;
         tmpArrayFive[2] = pTJet-ptUeFromCellMedian;
         tmpArrayFive[3] = triggerHadron->Pt();
         tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
         fHJetSpecSubUeMedian->Fill(tmpArrayFive);

         //Centrality, A, pTjet-pTbgCone, pTtrigg, dphi
         tmpArrayFive[0] = centValue;
         tmpArrayFive[1] = areaJet;
         tmpArrayFive[2] = pTJet - ptUeCone;
         tmpArrayFive[3] = triggerHadron->Pt();
         tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
         fHJetSpecSubUeCone->Fill(tmpArrayFive);

         //Ue diagnostics  "[A,pTjet,pTjet-pTUe, pTUe, rhoUe]", kT median
         tmpArrayFive[0] = areaJet;
         tmpArrayFive[1] = pTJet;
         tmpArrayFive[2] = pTJet - ptUeFromCellMedian;
         tmpArrayFive[3] = ptUeFromCellMedian;
         tmpArrayFive[4] = rhoFromCellMedian;
         fHJetUeMedian->Fill(tmpArrayFive);
 
         //Ue diagnostics  "[A,pTjet,pTjet-pTUe, pTUe, rhoUe]", Cone median
         tmpArrayFive[0] = areaJet;
         tmpArrayFive[1] = pTJet;
         tmpArrayFive[2] = pTJet - ptUeCone;
         tmpArrayFive[3] = ptUeCone;
         tmpArrayFive[4] = rhoCone;
         fHJetUeCone->Fill(tmpArrayFive);

         if(filledOnce){ //fill for each event only once
            Double_t fillRho[] = { rhoCone,rhoFromCellMedian};
            fHRhoUeMedianVsCone->Fill(fillRho);
            filledOnce = kFALSE;
         }
      }//jet loop
 
      //Fill Jet Density In the Event A>0.07
      /*if(injet>0){
         Double_t filldens[]={ (Double_t) particleList.GetEntries(),
                            injet/fkAcceptance,
                            triggerHadron->Pt()
                          };
         fHJetDensity->Fill(filldens);
      } 

      //Fill Jet Density In the Event A>0.4
      if(injet4>0){ 
         Double_t filldens4[]={ (Double_t) particleList.GetEntries(), 
                                injet4/fkAcceptance,
                                triggerHadron->Pt()
                              };
         fHJetDensityA4->Fill(filldens4);
      }*/
   }


   PostData(1, fOutputList);
}

//----------------------------------------------------------------------------
void AliAnalysisTaskJetCorePP::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if(fDebug) printf("AliAnalysisTaskJetCorePP DONE\n");
   if(!GetOutputData(1)) return;
}


//----------------------------------------------------------------------------
Int_t  AliAnalysisTaskJetCorePP::GetListOfTracks(TList *list){
   //Fill the list of accepted tracks (passed track cut)
   //return consecutive index of the hardest ch hadron in the list
   Int_t iCount        = 0;
   AliAODEvent *aodevt = NULL;

   if(!fESD) aodevt = fAODIn;
   else      aodevt = fAODOut;   

   if(!aodevt) return -1;

   Int_t    index = -1; //index of the highest particle in the list
   Double_t ptmax = -10;
   Int_t    triggers[200];
   Int_t    ntriggers = 0; //index in triggers array

   for(Int_t it = 0; it < aodevt->GetNumberOfTracks(); it++){
      AliAODTrack *tr = aodevt->GetTrack(it);
      
      //if((fFilterMask > 0) && !(tr->TestFilterBit(fFilterMask))) continue;
      if((fFilterMask > 0) && !(tr->IsHybridGlobalConstrainedGlobal())) continue;
      if(TMath::Abs((Float_t) tr->Eta()) > fTrackEtaCut) continue;
      if(tr->Pt() < fTrackLowPtCut) continue;
      list->Add(tr);

      //Search trigger:
      if(fHardest>0){ //leading particle 
         if(tr->Pt()>ptmax){ 
            ptmax = tr->Pt();	
            index = iCount;
         }
      }
    
      if(fHardest==0 && ntriggers<200){ //single inclusive 
         if(fTriggerPtRangeLow <= tr->Pt() && 
            tr->Pt() < fTriggerPtRangeHigh){ 
            triggers[ntriggers] = iCount;
            ntriggers++;
         }
      }

      iCount++;
   }

   if(fHardest==0 && ntriggers>0){ //select random inclusive trigger
      Int_t rnd = fRandom->Integer(ntriggers); //0 to ntriggers-1
      index     = triggers[rnd];
   }

   return index;
}

//----------------------------------------------------------------------------
/*
Double_t AliAnalysisTaskJetCorePP::GetBackgroundInPerpCone(Float_t jetR, Double_t jetPhi, Double_t jetEta, TList* trkList){
   //calculate sum of track pT in the cone perpendicular in phi to the jet 
   //jetR = cone radius
   // jetPhi, jetEta = direction of the jet 
   Int_t numberOfTrks = trkList->GetEntries();
   Double_t pTsum = 0.0;
   Double_t perpConePhi = jetPhi + TMath::Pi()/2;//perp cone w.r.t. jet in phi
   for(Int_t it=0; it<numberOfTrks; it++){
      AliVParticle *trk = (AliVParticle*) trkList->At(it); 
      Double_t dphi = RelativePhi(perpConePhi,trk->Phi());     
      Double_t deta = trk->Eta()-jetEta;     

      if( (dphi*dphi + deta*deta)< (jetR*jetR)){
         pTsum += trk->Pt();
      } 
   }

   return pTsum;
}
*/
//----------------------------------------------------------------------------

Double_t AliAnalysisTaskJetCorePP::RelativePhi(Double_t mphi,Double_t vphi){
   //Get relative azimuthal angle of two particles -pi to pi
   if      (vphi < -TMath::Pi()) vphi += TMath::TwoPi();
   else if (vphi > TMath::Pi())  vphi -= TMath::TwoPi();

   if      (mphi < -TMath::Pi()) mphi += TMath::TwoPi();
   else if (mphi > TMath::Pi())  mphi -= TMath::TwoPi();

   Double_t dphi = mphi - vphi;
   if      (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
   else if (dphi > TMath::Pi())  dphi -= TMath::TwoPi();

   return dphi;//dphi in [-Pi, Pi]
}


//----------------------------------------------------------------------------
Bool_t AliAnalysisTaskJetCorePP::SelectMCGenTracks(AliVParticle *trk, TList *trkList, Double_t &ptLeading, Int_t &index, Int_t counter){
   //fill track efficiency denominator
   if(trk->Pt() < fTrackLowPtCut) return kFALSE;
   if(trk->Charge()==0)        return kFALSE;
   if(TMath::Abs((Float_t) trk->Eta()) > fTrackEtaCut) return kFALSE;
   trkList->Add(trk);
   if(fFillRespMx) fhPtTrkTruePrimGen->Fill(trk->Pt(),trk->Eta()); //Efficiency denominator

   //Search trigger:
   if(fHardest>0){ //leading particle 
      if(ptLeading < trk->Pt()){
         index      = counter;
         ptLeading  = trk->Pt();
      }
   }

   if(fHardest==0){ //single inclusive 
      index = -1;  
      if(fTriggerPtRangeLow <= trk->Pt() && 
            trk->Pt() < fTriggerPtRangeHigh){
         index      = counter;
      }
   }

   return kTRUE;
} 

//----------------------------------------------------------------------------
void AliAnalysisTaskJetCorePP::FillEffHistos(TList *recList, TList *genList){
   //fill track efficiency numerator 
   if(!recList) return;
   if(!genList) return;
   Int_t nRec = recList->GetEntries();
   Int_t nGen = genList->GetEntries();
   for(Int_t ir=0; ir<nRec; ir++){ 
      AliAODTrack *trkRec = (AliAODTrack*) recList->At(ir);
      if(!trkRec) continue; 
      Int_t recLabel = TMath::Abs(trkRec->GetLabel());
      Bool_t matched = kFALSE;
 
      for(Int_t ig=0; ig<nGen; ig++){
	AliVParticle *trkGen = (AliVParticle*) genList->At(ig);  
        if(!trkGen) continue;
        Int_t genLabel = TMath::Abs(trkGen->GetLabel()); 
        if(recLabel==genLabel){
          fhPtTrkTruePrimRec->Fill(trkGen->Pt(),trkGen->Eta());
          matched = kTRUE;
          break;
        }
      }
      if(!matched) fhPtTrkSecOrFakeRec->Fill(trkRec->Pt(),trkRec->Eta()); 
   }

   return;
}
//________________________________________________________________
void AliAnalysisTaskJetCorePP::EstimateBgRhoMedian(TList *listJet, TList* listPart,  Double_t &rhoMedian){
   //Estimate background rho by means of integrating track pT outside identified jet cones

   rhoMedian = 0;

   //identify leading jet in the full track acceptance
   Int_t idxLeadingJet    = -1;
   Double_t pTleading     = 0.0; 
   AliAODJet* jet = NULL;
 
   for(Int_t ij = 0; ij < listJet->GetEntries(); ij++){ //loop over bg kT jets
      jet = (AliAODJet*)(listJet->At(ij));
      if(!jet) continue;
      if(TMath::Abs(jet->Eta()) > fTrackEtaCut) continue; //track acceptance cut
      if(jet->Pt() > pTleading){
         idxLeadingJet = ij; 
         pTleading     = jet->Pt();
      }
   } 

   Int_t  njets = 0;
   static Double_t jphi[100];
   static Double_t jeta[100];
   static Double_t jRsquared[100];

   for(Int_t ij = 0; ij< listJet->GetEntries(); ij++){ //loop over bg kT jets

      jet = (AliAODJet*)(listJet->At(ij));
      if(!jet) continue;
      if(TMath::Abs(jet->Eta()) > fTrackEtaCut) continue; 

      if((jet->Pt() > fBgMaxJetPt || ij== idxLeadingJet) && njets<100){
         jphi[njets]      = RelativePhi(jet->Phi(),0.0); //-pi,pi
         jeta[njets]      = jet->Eta();
         jRsquared[njets] = fSafetyMargin*jet->EffectiveAreaCharged()/TMath::Pi(); //scale jet radius
         njets++;
      }
   }

   
   static Double_t nOutCone[10][4];
   static Double_t sumPtOutOfCone[10][4];

   for(Int_t ie=0; ie<fnEta; ie++){
      for(Int_t ip=0; ip<fnPhi; ip++){
         nOutCone[ip][ie]       = 0.0;     //initialize counter
         sumPtOutOfCone[ip][ie] = 0.0;
      }
   }

   Double_t rndphi, rndeta;
   Double_t rndphishift, rndetashift;
   Double_t dphi, deta; 
   Bool_t   bIsInCone;

   for(Int_t it=0; it<fnTrials; it++){

      rndphi = fRandom->Uniform(0, fPhiSize);
      rndeta = fRandom->Uniform(0, fEtaSize);
                              
      for(Int_t ip=0; ip<fnPhi; ip++){  //move radom position to each cell
         rndphishift = rndphi + ip*fPhiSize - TMath::Pi();
         for(Int_t ie=0; ie<fnEta; ie++){
            rndetashift = rndeta + ie*fEtaSize - fEtaSize;

            bIsInCone = 0; //tag if trial is in the jet cone
            for(Int_t ij=0; ij<njets; ij++){
               deta = jeta[ij] - rndetashift;
               dphi = RelativePhi(rndphishift,jphi[ij]);
               if((dphi*dphi + deta*deta) < jRsquared[ij]){
                  bIsInCone = 1;
                  break;
               }
            }
            if(!bIsInCone) nOutCone[ip][ie]++;
         }
      }
   }
  

   //Sum particle pT outside jets in cells
   Int_t npart = listPart->GetEntries();
   Int_t phicell,etacell; 
   AliVParticle *part = NULL; 
   for(Int_t ip=0; ip < npart; ip++){ 
         
      part = (AliVParticle*) listPart->At(ip);
      if(!part) continue;
      if( TMath::Abs(part->Eta()) > fEtaSize) continue; //skip part outside +-0.7 in eta
 
      bIsInCone = 0; //init
      for(Int_t ij=0; ij<njets; ij++){
         dphi = RelativePhi(jphi[ij], part->Phi());
         deta = jeta[ij] - part->Eta();
         if((dphi*dphi + deta*deta) < jRsquared[ij]){
            bIsInCone = 1;
            break;
         }
      } 
      if(!bIsInCone){
         phicell = TMath::Nint(TMath::Floor((RelativePhi(part->Phi(),0.0) + TMath::Pi())/fPhiSize));
         etacell = TMath::Nint(TMath::Floor((part->Eta()+fEtaSize)/fEtaSize));
         sumPtOutOfCone[phicell][etacell]+= part->Pt();
      }
   }

   static Double_t rhoInCells[20];
   Double_t  relativeArea;
   Int_t  nCells=0;

   for(Int_t ip=0; ip<fnPhi; ip++){
      for(Int_t ie=0; ie<fnEta; ie++){
         relativeArea = nOutCone[ip][ie]/fnTrials;
         
         if(relativeArea < 0.5) continue; 
         rhoInCells[nCells] =  sumPtOutOfCone[ip][ie]/(relativeArea*fCellArea);
         nCells++;
      }
   }

   if(nCells>0)
     rhoMedian = TMath::Median(nCells, rhoInCells); 

}

//__________________________________________________________________
void AliAnalysisTaskJetCorePP::EstimateBgCone(TList *listJet, TList* listPart,  Double_t &rhoPerpCone){

   rhoPerpCone = 0.0; //init

   if(!listJet) return;
   Int_t njet  = listJet->GetEntries();
   Int_t npart = listPart->GetEntries();
   Double_t pTleading  = 0.0;
   Double_t phiLeading = 1000.;
   Double_t etaLeading = 1000.;
   
   for(Int_t jsig=0; jsig < njet; jsig++){ 
      AliAODJet* signaljet = (AliAODJet*)(listJet->At(jsig));
      if(!signaljet) continue;
      if((signaljet->Eta()<fJetEtaMin) && (fJetEtaMax<signaljet->Eta())) continue; //acceptance cut
      if(signaljet->Pt() >= pTleading){ //replace leading and subleading jet
         pTleading  = signaljet->Pt();
         phiLeading = signaljet->Phi();
         etaLeading = signaljet->Eta();
      }
   }   
  
   Double_t phileftcone  = phiLeading + TMath::Pi()/2; 
   Double_t phirightcone = phiLeading - TMath::Pi()/2; 
   Double_t dp, de;

   for(Int_t ip=0; ip < npart; ip++){ 
         
      AliVParticle *part = (AliVParticle*) listPart->At(ip);
      if(!part){
         continue;
      }
      
      dp = RelativePhi(phileftcone, part->Phi());
      de = etaLeading - part->Eta();
      if( sqrt(dp*dp + de*de)< fBgConeR) rhoPerpCone += part->Pt();

      dp = RelativePhi(phirightcone, part->Phi());
      if( sqrt(dp*dp + de*de)< fBgConeR) rhoPerpCone += part->Pt();
 
   }

   //normalize total pT by two times cone are 
   rhoPerpCone = rhoPerpCone/(2*TMath::Pi()*fBgConeR*fBgConeR);

}
//__________________________________________________________________

void AliAnalysisTaskJetCorePP::ReadTClonesArray(TString bname, TList *list){
   //Convert TClones array of jets to TList 

   if(!strlen(bname.Data())){
      AliError(Form("Jet branch %s not set.", bname.Data()));
      return;
   }

   TClonesArray *array=0x0;

   if(fAODOut&&!array){
      array = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(bname.Data()));
   }
   if(fAODExtension&&!array){
      array = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(bname.Data()));
   }
   if(fAODIn&&!array){
      array = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(bname.Data()));
   }

   list->Clear(); //clear list beforehand

   if(array){
      if(fDebug){
         Printf("########## %s: %d jets",bname.Data(), array->GetEntriesFast());
      }
      for(Int_t iJet = 0; iJet < array->GetEntriesFast(); iJet++) {
         AliAODJet *jet = dynamic_cast<AliAODJet*>((*array)[iJet]);
         if (jet) list->Add(jet);
      }
   }

   return;
}
