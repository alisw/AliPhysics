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
// AliAnalysisTaskSE for D0 candidates (2Prongs)
// and hadrons correlations
//
// Authors: Chiara Bianchin, chiara.bianchin@pd.infn.it (invariant mass)
// Fabio Colamaria, fabio.colamaria@ba.infn.it (correlations)
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
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
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSED0Correlations.h"
#include "AliNormalizationCounter.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSED0Correlations)


//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations():
AliAnalysisTaskSE(),
  fNPtBinsCorr(0), 
  fBinLimsCorr(),
  fPtThreshLow(),
  fPtThreshUp(), 
  fEvents(0),
  fAlreadyFilled(kFALSE),
  fOutputMass(0),
  fOutputCorr(0),
  fOutputStudy(0),
  fNentries(0), 
  fCutsD0(0),
  fCutsTracks(0),
  fCorrelatorTr(0),
  fCorrelatorKc(0),
  fCorrelatorK0(0),
  fReadMC(0),
  fRecoTr(kTRUE),
  fRecoD0(kTRUE),
  fSelEvType(kFALSE),
  fMixing(kFALSE),
  fCounter(0),
  fNPtBins(1),
  fFillOnlyD0D0bar(0),
  fIsSelectedCandidate(0),
  fSys(0),
  fEtaForCorrel(0),
  fIsRejectSDDClusters(0),
  fFillGlobal(kTRUE)
{
  // Default constructor

}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations(const char *name,AliRDHFCutsD0toKpi* cutsD0, AliHFAssociatedTrackCuts* cutsTrk):
  AliAnalysisTaskSE(name),
  fNPtBinsCorr(0),  
  fBinLimsCorr(),
  fPtThreshLow(),
  fPtThreshUp(),
  fEvents(0),
  fAlreadyFilled(kFALSE),
  fOutputMass(0),
  fOutputCorr(0),
  fOutputStudy(0),
  fNentries(0),
  fCutsD0(0),
  fCutsTracks(cutsTrk),
  fCorrelatorTr(0),
  fCorrelatorKc(0),
  fCorrelatorK0(0),
  fReadMC(0),
  fRecoTr(kTRUE),
  fRecoD0(kTRUE),
  fSelEvType(kFALSE),
  fMixing(kFALSE),
  fCounter(0),
  fNPtBins(1),
  fFillOnlyD0D0bar(0),
  fIsSelectedCandidate(0),
  fSys(0),
  fEtaForCorrel(0),
  fIsRejectSDDClusters(0),
  fFillGlobal(kTRUE)
{
  // Default constructor

  fNPtBins=cutsD0->GetNPtBins();
    
  fCutsD0=cutsD0;

  // Output slot #1 writes into a TList container (mass with cuts)
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TH1F container (number of events)
  DefineOutput(2,TH1F::Class());  //My private output
  // Output slot #3 writes into a AliRDHFD0toKpi container (cuts)
  DefineOutput(3,AliRDHFCutsD0toKpi::Class());  //My private output
  // Output slot #4 writes Normalization Counter 
  DefineOutput(4,AliNormalizationCounter::Class());
  // Output slot #5 writes into a TList container (correl output)
  DefineOutput(5,TList::Class());  //My private output
  // Output slot #6 writes into a TList container (correl advanced)
  DefineOutput(6,TList::Class());  //My private output
  // Output slot #7 writes into a AliHFAssociatedTrackCuts container (cuts)
  DefineOutput(7,AliHFAssociatedTrackCuts::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations(const AliAnalysisTaskSED0Correlations &source):
  AliAnalysisTaskSE(source),
  fNPtBinsCorr(source.fNPtBinsCorr), 
  fBinLimsCorr(source.fBinLimsCorr),
  fPtThreshLow(source.fPtThreshLow),
  fPtThreshUp(source.fPtThreshUp),
  fEvents(source.fEvents),
  fAlreadyFilled(source.fAlreadyFilled),
  fOutputMass(source.fOutputMass),
  fOutputCorr(source.fOutputCorr),
  fOutputStudy(source.fOutputStudy),
  fNentries(source.fNentries), 
  fCutsD0(source.fCutsD0),
  fCutsTracks(source.fCutsTracks),
  fCorrelatorTr(source.fCorrelatorTr),
  fCorrelatorKc(source.fCorrelatorKc),
  fCorrelatorK0(source.fCorrelatorK0),
  fReadMC(source.fReadMC),
  fRecoTr(source.fRecoTr),
  fRecoD0(source.fRecoD0),
  fSelEvType(source.fSelEvType),
  fMixing(source.fMixing),
  fCounter(source.fCounter),
  fNPtBins(source.fNPtBins),
  fFillOnlyD0D0bar(source.fFillOnlyD0D0bar),
  fIsSelectedCandidate(source.fIsSelectedCandidate),
  fSys(source.fSys),
  fEtaForCorrel(source.fEtaForCorrel),
  fIsRejectSDDClusters(source.fIsRejectSDDClusters),
  fFillGlobal(source.fFillGlobal)
{
  // Copy constructor
}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::~AliAnalysisTaskSED0Correlations()
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
  if (fCutsD0) {
    delete fCutsD0;
    fCutsD0 = 0;
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
}  

//______________________________________________________________________________
AliAnalysisTaskSED0Correlations& AliAnalysisTaskSED0Correlations::operator=(const AliAnalysisTaskSED0Correlations& orig)
{
// Assignment
  if (&orig == this) return *this; //if address is the same (same object), returns itself

  AliAnalysisTaskSE::operator=(orig); //Uses the AliAnalysisTaskSE operator to assign the inherited part of the class
  fNPtBinsCorr = orig.fNPtBinsCorr; 
  fBinLimsCorr = orig.fBinLimsCorr;
  fPtThreshLow = orig.fPtThreshLow;
  fPtThreshUp = orig.fPtThreshUp;
  fEvents = orig.fEvents;
  fAlreadyFilled = orig.fAlreadyFilled;
  fOutputMass = orig.fOutputMass;
  fOutputCorr = orig.fOutputCorr;
  fOutputStudy = orig.fOutputStudy;
  fNentries = orig.fNentries; 
  fCutsD0 = orig.fCutsD0;
  fCutsTracks = orig.fCutsTracks;
  fCorrelatorTr = orig.fCorrelatorTr;
  fCorrelatorKc = orig.fCorrelatorKc;
  fCorrelatorK0 = orig.fCorrelatorK0;
  fReadMC = orig.fReadMC;
  fRecoTr = orig.fRecoTr;
  fRecoD0 = orig.fRecoD0;
  fSelEvType = orig.fSelEvType;
  fMixing = orig.fMixing;
  fCounter = orig.fCounter;
  fNPtBins = orig.fNPtBins;
  fFillOnlyD0D0bar = orig.fFillOnlyD0D0bar;
  fIsSelectedCandidate = orig.fIsSelectedCandidate;
  fSys = orig.fSys;
  fEtaForCorrel = orig.fEtaForCorrel;
  fIsRejectSDDClusters = orig.fIsRejectSDDClusters;
  fFillGlobal = orig.fFillGlobal;

  return *this; //returns pointer of the class
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Correlations::Init() \n");
  
  //Copy of cuts objects
  AliRDHFCutsD0toKpi* copyfCutsD0 = new AliRDHFCutsD0toKpi(*fCutsD0);
  const char* nameoutput=GetOutputSlot(3)->GetContainer()->GetName();
  copyfCutsD0->SetName(nameoutput);

  //needer to clear completely the objects inside with Clear() method
  // Post the data
  PostData(3,copyfCutsD0);
  PostData(7,fCutsTracks);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::UserCreateOutputObjects()
{

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Correlations::UserCreateOutputObjects() \n");

  //HFCorrelator creation and definition
  fCorrelatorTr = new AliHFCorrelator("CorrelatorTr",fCutsTracks,fSys);
  fCorrelatorKc = new AliHFCorrelator("CorrelatorKc",fCutsTracks,fSys);
  fCorrelatorK0 = new AliHFCorrelator("CorrelatorK0",fCutsTracks,fSys);
  fCorrelatorTr->SetDeltaPhiInterval(-TMath::Pi()/2,3*TMath::Pi()/2);// set the Delta Phi Interval you want (in this case -0.5Pi to 1.5 Pi)
  fCorrelatorKc->SetDeltaPhiInterval(-TMath::Pi()/2,3*TMath::Pi()/2);
  fCorrelatorK0->SetDeltaPhiInterval(-TMath::Pi()/2,3*TMath::Pi()/2);
  fCorrelatorTr->SetEventMixing(fMixing);// sets the analysis on a single event (kFALSE) or mixed events (kTRUE)
  fCorrelatorKc->SetEventMixing(fMixing);
  fCorrelatorK0->SetEventMixing(fMixing);
  fCorrelatorTr->SetAssociatedParticleType(1);// set 1 for correlations with hadrons, 2 with kaons, 3 with KZeros
  fCorrelatorKc->SetAssociatedParticleType(2);// set 1 for correlations with hadrons, 2 with kaons, 3 with KZeros
  fCorrelatorK0->SetAssociatedParticleType(3);// set 1 for correlations with hadrons, 2 with kaons, 3 with KZeros
  fCorrelatorTr->SetApplyDisplacementCut(2); //0: don't calculate d0; 1: return d0; 2: return d0/d0err
  fCorrelatorKc->SetApplyDisplacementCut(2);
  fCorrelatorK0->SetApplyDisplacementCut(0);
  fCorrelatorTr->SetUseMC(fReadMC);// sets Montecarlo flag
  fCorrelatorKc->SetUseMC(fReadMC);
  fCorrelatorK0->SetUseMC(fReadMC);
  fCorrelatorTr->SetUseReco(fRecoTr);// sets (if MC analysis) wheter to analyze Reco or Kinem tracks
  fCorrelatorKc->SetUseReco(fRecoTr);
  fCorrelatorK0->SetUseReco(fRecoTr);
  fCorrelatorKc->SetPIDmode(2); //switch for K+/- PID option
  Bool_t pooldefTr = fCorrelatorTr->DefineEventPool();// method that defines the properties ot the event mixing (zVtx and Multipl. bins)
  Bool_t pooldefKc = fCorrelatorKc->DefineEventPool();// method that defines the properties ot the event mixing (zVtx and Multipl. bins)
  Bool_t pooldefK0 = fCorrelatorK0->DefineEventPool();// method that defines the properties ot the event mixing (zVtx and Multipl. bins)
  if(!pooldefTr) AliInfo("Warning:: Event pool not defined properly");
  if(!pooldefKc) AliInfo("Warning:: Event pool not defined properly");
  if(!pooldefK0) AliInfo("Warning:: Event pool not defined properly");

  // Several histograms are more conveniently managed in a TList
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("listMass");

  fOutputCorr = new TList();
  fOutputCorr->SetOwner();
  fOutputCorr->SetName("correlationslist");

  fOutputStudy = new TList();
  fOutputStudy->SetOwner();
  fOutputStudy->SetName("MCstudyplots");

  TString nameMass=" ",nameSgn=" ", nameBkg=" ", nameRfl=" ";

  for(Int_t i=0;i<fCutsD0->GetNPtBins();i++){

    nameMass="histMass_";
    nameMass+=i;
    nameSgn="histSgn_";
    nameSgn+=i;
    nameBkg="histBkg_";
    nameBkg+=i;
    nameRfl="histRfl_";
    nameRfl+=i;

    //histograms of invariant mass distributions

    //MC signal
    if(fReadMC){
      TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",120,1.5648,2.1648);
      tmpSt->Sumw2();

      //Reflection: histo filled with D0Mass which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
      TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries",120,1.5648,2.1648);
      TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",120,1.5648,2.1648);
      tmpBt->Sumw2();
      tmpRt->Sumw2();
      fOutputMass->Add(tmpSt);
      fOutputMass->Add(tmpRt);
      fOutputMass->Add(tmpBt);
    }

    //mass
    TH1F* tmpMt = new TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",120,1.5648,2.1648);
    tmpMt->Sumw2();
    fOutputMass->Add(tmpMt);
  }

  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();

  fNentries=new TH1F(nameoutput, "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of D0 selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 20,-0.5,19.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  fReadMC ? fNentries->GetXaxis()->SetBinLabel(3,"nD0Selected") : fNentries->GetXaxis()->SetBinLabel(3,"Dstar<-D0");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
  if(fSys==0) fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
  if(fReadMC && fSys==0){
    fNentries->GetXaxis()->SetBinLabel(12,"K");
    fNentries->GetXaxis()->SetBinLabel(13,"Lambda");
  }
  fNentries->GetXaxis()->SetBinLabel(14,"Pile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(15,"N. of 0SMH");
  if(fSys==1) fNentries->GetXaxis()->SetBinLabel(16,"Nev in centr");
  if(fIsRejectSDDClusters) fNentries->GetXaxis()->SetBinLabel(17,"SDD-Cls Rej");
  fNentries->GetXaxis()->SetBinLabel(18,"Phys.Sel.Rej");
  fNentries->GetXaxis()->SetBinLabel(19,"nEventsSelected");
  if(fReadMC) fNentries->GetXaxis()->SetBinLabel(20,"nEvsWithProdMech");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(4)->GetContainer()->GetName()));
  fCounter->Init();

  CreateCorrelationsObjs(); //creates histos for correlations analysis

  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fNentries);
  PostData(4,fCounter);
  PostData(5,fOutputCorr);
  PostData(6,fOutputStudy);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  //cout<<"I'm in UserExec"<<endl;


  //cuts order
  //     printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
  //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
  //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
  //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
  //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
  //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
  //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
  //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
  //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  fEvents++;

  TString bname="D0toKpi";

  TClonesArray *inputArray=0;

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
    printf("AliAnalysisTaskSED0Correlations::UserExec: input branch not found!\n");
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
      printf("AliAnalysisTaskSED0Correlations::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSED0Correlations::UserExec: MC header branch not found!\n");
      return;
    }
  }

  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fCutsD0,fReadMC); 

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);

  if(!fCutsD0->IsEventSelected(aod)) {
    if(fCutsD0->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
    if(fSys==1 && (fCutsD0->GetWhyRejection()==2 || fCutsD0->GetWhyRejection()==3)) fNentries->Fill(15);
    if(fCutsD0->GetWhyRejection()==7) fNentries->Fill(17);
    return;
  }

  fNentries->Fill(18); //event selected after selection

  //Setting PIDResponse for associated tracks
  fCorrelatorTr->SetPidAssociated();
  fCorrelatorKc->SetPidAssociated();
  fCorrelatorK0->SetPidAssociated();

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
      std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
      return; 
    }
    fNentries->Fill(19); //event with particular production type                
  
  } //end of selection


  // Check the Nb of SDD clusters
  if (fIsRejectSDDClusters) { 
    Bool_t skipEvent = kFALSE;
    Int_t ntracks = 0;
    if (aod) ntracks = aod->GetNTracks();
    for(Int_t itrack=0; itrack<ntracks; itrack++) { // loop on tacks
      //    ... get the track
      AliAODTrack * track = aod->GetTrack(itrack);
      if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
	skipEvent=kTRUE;
	fNentries->Fill(16);
	break;
      }
    }
    if (skipEvent) return;
  }
  
  //HFCorrelators initialization (for this event)
  fCorrelatorTr->SetAODEvent(aod); // set the AOD event from which you are processing
  fCorrelatorKc->SetAODEvent(aod);
  fCorrelatorK0->SetAODEvent(aod);
  Bool_t correlatorONTr = fCorrelatorTr->Initialize(); // initialize the pool for event mixing
  Bool_t correlatorONKc = fCorrelatorKc->Initialize();
  Bool_t correlatorONK0 = fCorrelatorK0->Initialize();
  if(!correlatorONTr) {AliInfo("AliHFCorrelator (tracks) didn't initialize the pool correctly or processed a bad event"); return;}
  if(!correlatorONKc) {AliInfo("AliHFCorrelator (charged K) didn't initialize the pool correctly or processed a bad event"); return;}
  if(!correlatorONK0) {AliInfo("AliHFCorrelator (K0) didn't initialize the pool correctly or processed a bad event"); return;}

  if(fReadMC) {
    fCorrelatorTr->SetMCArray(mcArray); // set the TClonesArray *fmcArray for analysis on monte carlo
    fCorrelatorKc->SetMCArray(mcArray);
    fCorrelatorK0->SetMCArray(mcArray);
  }

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  //vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    fNentries->Fill(3);
  }

  //Reset flag for tracks distributions fill
  fAlreadyFilled=kFALSE;

  //***** Loop over D0 candidates *****
  Int_t nInD0toKpi = inputArray->GetEntriesFast();
  if(fDebug>2) printf("Number of D0->Kpi: %d\n",nInD0toKpi);

  if(fFillGlobal) { //loop on V0 and tracks for each event, to fill Pt distr. and InvMass distr.

    TClonesArray *v0array = (TClonesArray*)aod->GetList()->FindObject("v0s");
    Int_t pdgCodes[2] = {211,211};
    Int_t idArrayV0[v0array->GetEntriesFast()][2];
    for(int iV0=0; iV0<v0array->GetEntriesFast(); iV0++) {
      for(int j=0; j<2; j++) {idArrayV0[iV0][j]=-2;}
      AliAODv0 *v0 = (AliAODv0*)v0array->UncheckedAt(iV0);
      if(SelectV0(v0,vtx1,2,idArrayV0)) { //option 2 = for mass inv plots only
        if(fReadMC && fRecoTr && (v0->MatchToMC(310,mcArray,2,pdgCodes)<0)) continue; //310 = K0s, 311 = K0 generico!!
        ((TH2F*)fOutputStudy->FindObject("hK0MassInv"))->Fill(v0->MassK0Short(),v0->Pt()); //invariant mass plot
        ((TH1F*)fOutputStudy->FindObject("hist_Pt_K0_AllEv"))->Fill(v0->Pt()); //pT distribution (in all events), K0 case
      }
    }

    for(Int_t itrack=0; itrack<aod->GetNTracks(); itrack++) { // loop on tacks
      AliAODTrack * track = aod->GetTrack(itrack);
      //rejection of tracks
      if(track->GetID() < 0) continue; //discard negative ID tracks
      if(track->Pt() < fPtThreshLow.at(0) || track->Pt() > fPtThreshUp.at(0)) continue; //discard tracks outside pt range for hadrons/K
      if(!fCutsTracks->IsHadronSelected(track) || !fCutsTracks->CheckHadronKinematic(track->Pt(),0.1)) continue; //0.1 = dummy (d0 value, no cut on it for me)
      //pT distribution (in all events), charg and hadr cases
      ((TH1F*)fOutputStudy->FindObject("hist_Pt_Charg_AllEv"))->Fill(track->Pt()); 
      if(fCutsTracks->CheckKaonCompatibility(track,kFALSE,0,2)) ((TH1F*)fOutputStudy->FindObject("hist_Pt_Kcharg_AllEv"))->Fill(track->Pt());
    }
 
  } //end of loops for global plot fill

  Int_t nSelectedloose=0,nSelectedtight=0;  
  
  //RecoD0 case ************************************************
  if(fRecoD0) {

    for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
      AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iD0toKpi);
 
      if(d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
  	fNentries->Fill(2);
  	continue; //skip the D0 from Dstar  
      }
    
      if (fCutsD0->IsInFiducialAcceptance(d->Pt(),d->Y(421)) ) {
        nSelectedloose++;
        nSelectedtight++;      
        if(fSys==0){
  	  if(fCutsD0->IsSelected(d,AliRDHFCuts::kTracks,aod))fNentries->Fill(6);       
        }  
        Int_t ptbin=fCutsD0->PtBin(d->Pt());
        if(ptbin==-1) {fNentries->Fill(4); continue;} //out of bounds

        fIsSelectedCandidate=fCutsD0->IsSelected(d,AliRDHFCuts::kAll,aod); //D0 selected
        if(!fIsSelectedCandidate) continue;

        //D0 infos
        Double_t phiD0 = fCorrelatorTr->SetCorrectPhiRange(d->Phi());
                 phiD0 = fCorrelatorKc->SetCorrectPhiRange(d->Phi());  //bad usage, but returns a Double_t...
                 phiD0 = fCorrelatorK0->SetCorrectPhiRange(d->Phi());
        fCorrelatorTr->SetTriggerParticleProperties(d->Pt(),phiD0,d->Eta()); // sets the parameters of the trigger particles that are needed
        fCorrelatorKc->SetTriggerParticleProperties(d->Pt(),phiD0,d->Eta());
        fCorrelatorK0->SetTriggerParticleProperties(d->Pt(),phiD0,d->Eta());
        fCorrelatorTr->SetD0Properties(d,fIsSelectedCandidate); //sets special properties for D0
        fCorrelatorKc->SetD0Properties(d,fIsSelectedCandidate);
        fCorrelatorK0->SetD0Properties(d,fIsSelectedCandidate);

        if(!fReadMC) {
          if (TMath::Abs(d->Eta())<fEtaForCorrel) CalculateCorrelations(d); //correlations on real data
        } else { //correlations on MC -> association of selected D0 to MCinfo with MCtruth
          if (TMath::Abs(d->Eta())<fEtaForCorrel) {
            Int_t pdgDgD0toKpi[2]={321,211};
    	    Int_t labD0 = d->MatchToMC(421,mcArray,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not
            if (labD0>-1) CalculateCorrelations(d,labD0,mcArray);
          }
        }

        FillMassHists(d,mcArray,fCutsD0,fOutputMass);
      }
    }
  }
  //End RecoD0 case ************************************************
  
  //MCKineD0 case ************************************************
  if(fReadMC && !fRecoD0) {

    for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { //Loop over all the tracks of MCArray
      AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
      if (!mcPart) {
        AliWarning("Particle not found in tree, skipping"); 
        continue;
      } 
  
      if(TMath::Abs(mcPart->GetPdgCode()) == 421){  // THIS IS A D0
        if (fCutsD0->IsInFiducialAcceptance(mcPart->Pt(),mcPart->Y()) ) {
          nSelectedloose++;
          nSelectedtight++;      
          if(fSys==0) fNentries->Fill(6);
          Int_t ptbin=fCutsD0->PtBin(mcPart->Pt());
          if(ptbin==-1) {fNentries->Fill(4); continue;} //out of bounds  
  
          //D0 infos
          Double_t phiD0 = fCorrelatorTr->SetCorrectPhiRange(mcPart->Phi());
                   phiD0 = fCorrelatorKc->SetCorrectPhiRange(mcPart->Phi());  //bad usage, but returns a Double_t...
                   phiD0 = fCorrelatorK0->SetCorrectPhiRange(mcPart->Phi());
          fCorrelatorTr->SetTriggerParticleProperties(mcPart->Pt(),phiD0,mcPart->Eta()); // sets the parameters of the trigger particles that are needed
          fCorrelatorKc->SetTriggerParticleProperties(mcPart->Pt(),phiD0,mcPart->Eta());
          fCorrelatorK0->SetTriggerParticleProperties(mcPart->Pt(),phiD0,mcPart->Eta());
          //fCorrelatorTr->SetD0Properties(mcPart,fIsSelectedCandidate); //needed for D* soft pions rejection, useless in MCKine
          //fCorrelatorKc->SetD0Properties(mcPart,fIsSelectedCandidate);
          //fCorrelatorK0->SetD0Properties(mcPart,fIsSelectedCandidate);
  
          if (TMath::Abs(mcPart->Eta())<fEtaForCorrel) {
  
            //Removal of D0 from D* feeddown! This solves also the problem of soft pions, now excluded
         /*   Int_t mother = mcPart->GetMother();
  	    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother));
            if(!mcMoth) continue;
	    if(TMath::Abs(mcMoth->GetPdgCode())==413) continue;
	*/
            if (mcPart->GetPdgCode()==421) fIsSelectedCandidate = 1;
    	    else fIsSelectedCandidate = 2;

	    TString fillthis="histSgn_"; fillthis+=ptbin;
	    ((TH1F*)(fOutputMass->FindObject(fillthis)))->Fill(1.864);
	  
            CalculateCorrelationsMCKine(mcPart,mcArray);
          }
        }
      }
    }
  }
  //End MCKineD0 case ************************************************

  if(fMixing /* && fAlreadyFilled*/) { // update the pool for Event Mixing, if: enabled,  event is ok, at least a SelD0 found! (fAlreadyFilled's role!)
    Bool_t updatedTr = fCorrelatorTr->PoolUpdate();
    Bool_t updatedKc = fCorrelatorKc->PoolUpdate();
    Bool_t updatedK0 = fCorrelatorK0->PoolUpdate();
    if(!updatedTr || !updatedKc || !updatedK0) AliInfo("Pool was not updated");
  }

  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);  
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);  

  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fNentries);
  PostData(4,fCounter);
  PostData(5,fOutputCorr);
  PostData(6,fOutputStudy);

  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0Correlations::FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi* cuts, TList *listout){
  //
  // function used in UserExec to fill mass histograms:
  //
  if (!fIsSelectedCandidate) return;

  if(fDebug>2)  cout<<"Candidate selected"<<endl;

  Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
  Int_t ptbin = cuts->PtBin(part->Pt());
  
  TString fillthis="";
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t labD0=-1;
  if (fReadMC) labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

  //count candidates selected by cuts
  fNentries->Fill(1);
  //count true D0 selected by cuts
  if (fReadMC && labD0>=0) fNentries->Fill(2);

  if ((fIsSelectedCandidate==1 || fIsSelectedCandidate==3) && fFillOnlyD0D0bar<2) { //D0

    if(fReadMC){
      if(labD0>=0) {
	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	if (pdgD0==421){ //D0
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	} else{ //it was a D0bar
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	}
      } else {//background
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      }
    }else{
      fillthis="histMass_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
    }
     
  }
  if (fIsSelectedCandidate>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)) { //D0bar

    if(fReadMC){
      if(labD0>=0) {
	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();

	if (pdgD0==-421){ //D0bar
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	} else{
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	}
      } else {//background or LS
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
      }
    }else{
      fillthis="histMass_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(invmassD0bar);
    }
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Correlations: Terminate() \n");

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

  fCutsD0 = dynamic_cast<AliRDHFCutsD0toKpi*>(GetOutputData(3));
  if(!fCutsD0){
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
Int_t AliAnalysisTaskSED0Correlations::CheckD0Origin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
  printf("AliAnalysisTaskSED0Correlations::CheckD0Origin() \n");
	
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
void AliAnalysisTaskSED0Correlations::CreateCorrelationsObjs() {
//

  TString namePlot = "";

  //These for limits in THnSparse (one per bin, same limits). 
  //Vars: DeltaPhi, InvMass, PtTrack, Displacement, DeltaEta
  Int_t nBinsPhi[5] = {32,150,6,3,16};
  Double_t binMinPhi[5] = {-TMath::Pi()/2,1.6,0.,0.,-1.6};  //is the minimum for all the bins
  Double_t binMaxPhi[5] = {3*TMath::Pi()/2,2.2,3.0,3.,1.6};  //is the maximum for all the bins

  //Vars: DeltaPhi, InvMass, DeltaEta
  Int_t nBinsMix[3] = {32,150,16};
  Double_t binMinMix[3] = {-TMath::Pi()/2,1.6,-1.6};  //is the minimum for all the bins
  Double_t binMaxMix[3] = {3*TMath::Pi()/2,2.2,1.6};  //is the maximum for all the bins

  for(Int_t i=0;i<fNPtBinsCorr;i++){

    if(!fMixing) {
      //THnSparse plots: correlations for various invariant mass (MC and data)
      namePlot="hPhi_K0_Bin";
      namePlot+=i;

      THnSparseF *hPhiK = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiK->Sumw2();
      fOutputCorr->Add(hPhiK);

      namePlot="hPhi_Kcharg_Bin";
      namePlot+=i;

      THnSparseF *hPhiH = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiH->Sumw2();
      fOutputCorr->Add(hPhiH);

      namePlot="hPhi_Charg_Bin";
      namePlot+=i;

      THnSparseF *hPhiC = new THnSparseF(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiC->Sumw2();
      fOutputCorr->Add(hPhiC);
  
      //histos for c/b origin for D0 (MC only)
      if (fReadMC) {

        //generic origin for tracks
        namePlot="hPhi_K0_From_c_Bin";
        namePlot+=i;

        THnSparseF *hPhiK_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK_c->Sumw2();
        fOutputCorr->Add(hPhiK_c);

        namePlot="hPhi_Kcharg_From_c_Bin";
        namePlot+=i;

        THnSparseF *hPhiH_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH_c->Sumw2();
        fOutputCorr->Add(hPhiH_c);

        namePlot="hPhi_Charg_From_c_Bin";
        namePlot+=i;

        THnSparseF *hPhiC_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiC_c->Sumw2();
        fOutputCorr->Add(hPhiC_c);
  
        namePlot="hPhi_K0_From_b_Bin";
        namePlot+=i;

        THnSparseF *hPhiK_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK_b->Sumw2();
        fOutputCorr->Add(hPhiK_b);

        namePlot="hPhi_Kcharg_From_b_Bin";
        namePlot+=i;

        THnSparseF *hPhiH_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH_b->Sumw2();
        fOutputCorr->Add(hPhiH_b);

        namePlot="hPhi_Charg_From_b_Bin";
        namePlot+=i;

        THnSparseF *hPhiC_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiC_b->Sumw2();
        fOutputCorr->Add(hPhiC_b);

        //HF-only tracks (c for c->D0, b for b->D0)
        namePlot="hPhi_K0_HF_From_c_Bin";
        namePlot+=i;

        THnSparseF *hPhiK_HF_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK_HF_c->Sumw2();
        fOutputCorr->Add(hPhiK_HF_c);

        namePlot="hPhi_Kcharg_HF_From_c_Bin";
        namePlot+=i;

        THnSparseF *hPhiH_HF_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH_HF_c->Sumw2();
        fOutputCorr->Add(hPhiH_HF_c);

        namePlot="hPhi_Charg_HF_From_c_Bin";
        namePlot+=i;

        THnSparseF *hPhiC_HF_c = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiC_HF_c->Sumw2();
        fOutputCorr->Add(hPhiC_HF_c);

        namePlot="hPhi_K0_HF_From_b_Bin";
        namePlot+=i;

        THnSparseF *hPhiK_HF_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK_HF_b->Sumw2();
        fOutputCorr->Add(hPhiK_HF_b);

        namePlot="hPhi_Kcharg_HF_From_b_Bin";
        namePlot+=i;

        THnSparseF *hPhiH_HF_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH_HF_b->Sumw2();
        fOutputCorr->Add(hPhiH_HF_b);

        namePlot="hPhi_Charg_HF_From_b_Bin";
        namePlot+=i;

        THnSparseF *hPhiC_HF_b = new THnSparseF(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiC_HF_b->Sumw2();
        fOutputCorr->Add(hPhiC_HF_b);

        namePlot="hPhi_K0_NonHF_Bin";
        namePlot+=i;

        THnSparseF *hPhiK_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiK_Non->Sumw2();
        fOutputCorr->Add(hPhiK_Non);

        namePlot="hPhi_Kcharg_NonHF_Bin";
        namePlot+=i;

        THnSparseF *hPhiH_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiH_Non->Sumw2();
        fOutputCorr->Add(hPhiH_Non);

        namePlot="hPhi_Charg_NonHF_Bin";
        namePlot+=i;

        THnSparseF *hPhiC_Non = new THnSparseF(namePlot.Data(), "Azimuthal correlation - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",5,nBinsPhi,binMinPhi,binMaxPhi);
        hPhiC_Non->Sumw2();
        fOutputCorr->Add(hPhiC_Non);
      }

      //leading hadron correlations
      namePlot="hPhi_Lead_Bin";
      namePlot+=i;

      THnSparseF *hCorrLead = new THnSparseF(namePlot.Data(), "Leading particle correlations; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      hCorrLead->Sumw2();
      fOutputCorr->Add(hCorrLead);

      if (fReadMC) {
        namePlot="hPhi_Lead_From_c_Bin";
        namePlot+=i;

        THnSparseF *hCorrLead_c = new THnSparseF(namePlot.Data(), "Leading particle correlations - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        hCorrLead_c->Sumw2();
        fOutputCorr->Add(hCorrLead_c);
  
        namePlot="hPhi_Lead_From_b_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrLead_b = new THnSparseF(namePlot.Data(), "Leading particle correlations - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        hCorrLead_b->Sumw2();
        fOutputCorr->Add(hCorrLead_b);
  
        namePlot="hPhi_Lead_HF_From_c_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrLead_HF_c = new THnSparseF(namePlot.Data(), "Leading particle correlations HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        hCorrLead_HF_c->Sumw2();
        fOutputCorr->Add(hCorrLead_HF_c);
  
        namePlot="hPhi_Lead_HF_From_b_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrLead_HF_b = new THnSparseF(namePlot.Data(), "Leading particle correlations HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        hCorrLead_HF_b->Sumw2();
        fOutputCorr->Add(hCorrLead_HF_b);

        namePlot="hPhi_Lead_NonHF_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrLead_Non = new THnSparseF(namePlot.Data(), "Leading particle correlations - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        hCorrLead_Non->Sumw2();
        fOutputCorr->Add(hCorrLead_Non);
      }
      
      //pT weighted correlations
      namePlot="hPhi_Weig_Bin";
      namePlot+=i;
  
      THnSparseF *hCorrWeig = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted); #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      fOutputCorr->Add(hCorrWeig);
  
      if (fReadMC) {
        namePlot="hPhi_Weig_From_c_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrWeig_c = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        fOutputCorr->Add(hCorrWeig_c);
  
        namePlot="hPhi_Weig_From_b_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrWeig_b = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        fOutputCorr->Add(hCorrWeig_b);
  
        namePlot="hPhi_Weig_HF_From_c_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrWeig_HF_c = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        fOutputCorr->Add(hCorrWeig_HF_c);
  
        namePlot="hPhi_Weig_HF_From_b_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrWeig_HF_b = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        fOutputCorr->Add(hCorrWeig_HF_b);

        namePlot="hPhi_Weig_NonHF_Bin";
        namePlot+=i;
  
        THnSparseF *hCorrWeig_Non = new THnSparseF(namePlot.Data(), "Charged particle correlations (pT weighted) - Non HF; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
        fOutputCorr->Add(hCorrWeig_Non);
      }

    //pT distribution histos
    namePlot = "hist_Pt_Charg_Bin"; namePlot+=i;
    TH1F *hPtC = new TH1F(namePlot.Data(), "Charged track pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
    hPtC->SetMinimum(0);
    fOutputStudy->Add(hPtC);

    namePlot = "hist_Pt_Kcharg_Bin"; namePlot+=i;
    TH1F *hPtH = new TH1F(namePlot.Data(), "Hadrons pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
    hPtH->SetMinimum(0);
    fOutputStudy->Add(hPtH);

    namePlot = "hist_Pt_K0_Bin"; namePlot+=i;
    TH1F *hPtK = new TH1F(namePlot.Data(), "Kaons pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
    hPtK->SetMinimum(0);
    fOutputStudy->Add(hPtK);

    //D* feeddown pions rejection histos
    namePlot = "hDstarPions_Bin"; namePlot+=i;
    TH2F *hDstarPions = new TH2F(namePlot.Data(), "Tracks rejected for D* inv.mass cut; # Tracks",2,0.,2.,150,1.6,2.2);
    hDstarPions->GetXaxis()->SetBinLabel(1,"Not rejected");
    hDstarPions->GetXaxis()->SetBinLabel(2,"Rejected");
    hDstarPions->SetMinimum(0);
    fOutputStudy->Add(hDstarPions); 
 
    }

    if(fMixing) {
      //THnSparse plots for event mixing!
      namePlot="hPhi_K0_Bin";
      namePlot+=i;namePlot+="_EvMix";

      THnSparseF *hPhiK_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      hPhiK_EvMix->Sumw2();
      fOutputCorr->Add(hPhiK_EvMix);

      namePlot="hPhi_Kcharg_Bin";
      namePlot+=i;namePlot+="_EvMix";
  
      THnSparseF *hPhiH_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      hPhiH_EvMix->Sumw2();
      fOutputCorr->Add(hPhiH_EvMix);

      namePlot="hPhi_Charg_Bin";
      namePlot+=i;namePlot+="_EvMix";

      THnSparseF *hPhiC_EvMix = new THnSparseF(namePlot.Data(), "Az. corr. EvMix; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",3,nBinsMix,binMinMix,binMaxMix);
      hPhiC_EvMix->Sumw2();
      fOutputCorr->Add(hPhiC_EvMix);
    }
  }

  //out of bin loop
  if(!fMixing) {
    TH1F *hCountC = new TH1F("hist_Count_Charg", "Charged track counter; # Tracks",100,0.,100.);
    hCountC->SetMinimum(0);
    fOutputStudy->Add(hCountC);

    TH1F *hCountH = new TH1F("hist_Count_Kcharg", "Hadrons counter; # Tracks",20,0.,20.);
    hCountH->SetMinimum(0);
    fOutputStudy->Add(hCountH);

    TH1F *hCountK = new TH1F("hist_Count_K0", "Kaons counter; # Tracks",20,0.,20.);
    hCountK->SetMinimum(0);
    fOutputStudy->Add(hCountK);
  }

  if (fReadMC) {
    TH1D *hEventTypeMC = new TH1D("EventTypeMC","EventTypeMC",100,-0.5,99.5);
    fOutputStudy->Add(hEventTypeMC); 
  }

  if (fFillGlobal) { //all-events plots
    //pt distributions
    TH1F *hPtCAll = new TH1F("hist_Pt_Charg_AllEv", "Charged track pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtCAll->SetMinimum(0);
    fOutputStudy->Add(hPtCAll);

    TH1F *hPtHAll = new TH1F("hist_Pt_Kcharg_AllEv", "Hadrons pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtHAll->SetMinimum(0);
    fOutputStudy->Add(hPtHAll);

    TH1F *hPtKAll = new TH1F("hist_Pt_K0_AllEv", "Kaons pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtKAll->SetMinimum(0);
    fOutputStudy->Add(hPtKAll);

    if(!fMixing) {
      //phi distributions
      TH1F *hPhiDistCAll = new TH1F("hist_PhiDistr_Charg", "Charged track phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
      hPhiDistCAll->SetMinimum(0);
      fOutputStudy->Add(hPhiDistCAll);

      TH1F *hPhiDistHAll = new TH1F("hist_PhiDistr_Kcharg", "Hadrons phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
      hPhiDistHAll->SetMinimum(0);
      fOutputStudy->Add(hPhiDistHAll);

      TH1F *hPhiDistKAll = new TH1F("hist_PhiDistr_K0", "Kaons phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
      hPhiDistKAll->SetMinimum(0);
      fOutputStudy->Add(hPhiDistKAll);

      TH1F *hPhiDistDAll = new TH1F("hist_PhiDistr_D0", "D^{0} phi distr. (All); p_{T} (GeV/c)",64,0,6.283);
      hPhiDistDAll->SetMinimum(0);
      fOutputStudy->Add(hPhiDistDAll);
    }

    //K0 Invariant Mass plots
    TH2F *hK0MassInv = new TH2F("hK0MassInv", "K0 invariant mass; Invariant mass (MeV/c^{2}); pT (GeV/c)",200,0.4,0.6,100,0.,10.);
    hK0MassInv->SetMinimum(0);
    fOutputStudy->Add(hK0MassInv);
  }

  //for MC analysis only
  for(Int_t i=0;i<fNPtBinsCorr;i++) {

    if (fReadMC && !fMixing) {

      //displacement histos
      namePlot="histDispl_K0_Bin"; namePlot+=i;
      TH1F *hDisplK = new TH1F(namePlot.Data(), "Kaons Displacement; DCA",150,0.,0.15);
      hDisplK->SetMinimum(0);
      fOutputStudy->Add(hDisplK);
  
      namePlot="histDispl_K0_HF_Bin";  namePlot+=i;
      TH1F *hDisplK_HF = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplK_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF);

      namePlot="histDispl_Kcharg_Bin"; namePlot+=i;
      TH1F *hDisplHadr = new TH1F(namePlot.Data(), "Hadrons Displacement; DCA",150,0.,0.15);
      hDisplHadr->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr);
  
      namePlot="histDispl_Kcharg_HF_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplHadr_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF);

      namePlot="histDispl_Charg_Bin"; namePlot+=i;
      TH1F *hDisplCharg = new TH1F(namePlot.Data(), "Charged tracks Displacement; DCA",150,0.,0.15);
      hDisplCharg->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg);
  
      namePlot="histDispl_Charg_HF_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplCharg_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF);

      namePlot="histDispl_K0_From_c_Bin"; namePlot+=i;
      TH1F *hDisplK_c = new TH1F(namePlot.Data(), "Kaons Displacement - c origin; DCA",150,0.,0.15);
      hDisplK_c->SetMinimum(0);
      fOutputStudy->Add(hDisplK_c);
  
      namePlot="histDispl_K0_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplK_HF_c = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplK_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF_c);

      namePlot="histDispl_Kcharg_From_c_Bin"; namePlot+=i;
      TH1F *hDisplHadr_c = new TH1F(namePlot.Data(), "Hadrons Displacement - c origin; DCA",150,0.,0.15);
      hDisplHadr_c->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_c);
  
      namePlot="histDispl_Kcharg_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF_c = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplHadr_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF_c);

      namePlot="histDispl_Charg_From_c_Bin"; namePlot+=i;
      TH1F *hDisplCharg_c = new TH1F(namePlot.Data(), "Charged tracks Displacement - c origin; DCA",150,0.,0.15);
      hDisplCharg_c->Sumw2();
      hDisplCharg_c->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_c);
  
      namePlot="histDispl_Charg_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF_c = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplCharg_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF_c);

      namePlot="histDispl_K0_From_b_Bin"; namePlot+=i;
      TH1F *hDisplK_b = new TH1F(namePlot.Data(), "Kaons Displacement - b origin; DCA",150,0.,0.15);
      hDisplK_b->SetMinimum(0);
      fOutputStudy->Add(hDisplK_b);
  
      namePlot="histDispl_K0_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplK_HF_b = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplK_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF_b);

      namePlot="histDispl_Kcharg_From_b_Bin"; namePlot+=i;
      TH1F *hDisplHadr_b = new TH1F(namePlot.Data(), "Hadrons Displacement - b origin; DCA",150,0.,0.15);
      hDisplHadr_b->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_b);

      namePlot="histDispl_Kcharg_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF_b = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplHadr_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF_b);

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

      namePlot="histOrig_Kcharg_Bin";  namePlot+=i;
      TH1F *hOrigin_Kcharg = new TH1F(namePlot.Data(), "Origin of hadrons",9,0.,9.);
      hOrigin_Kcharg->SetMinimum(0);
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(4,"c hadr.");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(5,"B->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(6,"B->X-># (X!=D)");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(7,"B->D->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(8,"B->D->X->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_Kcharg);

      namePlot="histOrig_K0_Bin";  namePlot+=i;
      TH1F *hOrigin_K = new TH1F(namePlot.Data(), "Origin of kaons",9,0.,9.);
      hOrigin_K->SetMinimum(0);
      hOrigin_K->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_K->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_K->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_K->GetXaxis()->SetBinLabel(4,"c hadr.");
      hOrigin_K->GetXaxis()->SetBinLabel(5,"B->#");
      hOrigin_K->GetXaxis()->SetBinLabel(6,"B->X-># (X!=D)");
      hOrigin_K->GetXaxis()->SetBinLabel(7,"B->D->#");
      hOrigin_K->GetXaxis()->SetBinLabel(8,"B->D->X->#");
      hOrigin_K->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_K);
    }

    if (fReadMC) {
      //origin of D0 histos
      namePlot="histOrig_D0_Bin";  namePlot+=i;
      TH1F *hOrigin_D0 = new TH1F(namePlot.Data(), "Origin of D0",2,0.,2.);
      hOrigin_D0->SetMinimum(0);
      hOrigin_D0->GetXaxis()->SetBinLabel(1,"From c");
      hOrigin_D0->GetXaxis()->SetBinLabel(2,"From b");
      fOutputStudy->Add(hOrigin_D0);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::CalculateCorrelations(AliAODRecoDecayHF2Prong* d, Int_t labD0, TClonesArray* mcArray) {
//
// Method for correlations D0-hadrons study
//
  Int_t N_Charg = 0, N_KCharg = 0, N_Kaons = 0;
  Double_t mD0, mD0bar;
  Int_t origD0 = 0, PDGD0 = 0, ptbin = 0;
  d->InvMassD0(mD0,mD0bar);
  Double_t mInv[2] = {mD0, mD0bar};
  ptbin = PtBinCorr(d->Pt());

  if(ptbin < 0) return;

  //Fill of D0 phi distribution
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_PhiDistr_D0"))->Fill(d->Phi());  

  //Origin of D0
  TString orig="";
  if(fReadMC) {
    origD0=CheckD0Origin(mcArray,(AliAODMCParticle*)mcArray->At(labD0));
    PDGD0 = ((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode();
    switch (CheckD0Origin(mcArray,(AliAODMCParticle*)mcArray->At(labD0))) {
      case 4:
        orig = "_From_c";
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(0.);
        break;
      case 5:
        orig = "_From_b";
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(1.);
        break;
      default:
        return;
    }
  }

  Double_t highPt = 0; Double_t lead[4] = {0,0,0,1};  //infos for leading particle (pt,deltaphi)

  //loop over the tracks in the pool 
  Bool_t execPoolTr = fCorrelatorTr->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
  Bool_t execPoolKc = fCorrelatorKc->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
  Bool_t execPoolK0 = fCorrelatorK0->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
		
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
	/*if(!track->CheckSoftPi()) { //removal of soft pions
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(1.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(1.,mD0bar);
          continue;
        } else { //not a soft pion
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(0.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(0.,mD0bar);
        }*/
        Int_t idDaughs[2] = {((AliVTrack*)d->GetDaughter(0))->GetID(),((AliVTrack*)d->GetDaughter(1))->GetID()}; //IDs of daughters to be skipped
        if(track->GetID() == idDaughs[0] || track->GetID() == idDaughs[1]) continue; //discards daughters of candidate
      }
      if(track->Pt() < fPtThreshLow.at(ptbin) || track->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K

      Double_t eff = track->GetWeight(); //extract track efficiency

      FillSparsePlots(mcArray,mInv,origD0,PDGD0,track,ptbin,kTrack,1./eff); //fills for charged tracks, weight = 1./eff

      if(!fMixing) N_Charg++;

      //retrieving leading info...
      if(track->Pt() > highPt) {
        if(fReadMC && track->GetLabel()<1) continue;
        if(fReadMC && !(AliAODMCParticle*)mcArray->At(track->GetLabel())) continue;
        lead[0] = fCorrelatorTr->GetDeltaPhi();
        lead[1] = fCorrelatorTr->GetDeltaEta();
        lead[2] = fReadMC ? CheckTrackOrigin(mcArray,(AliAODMCParticle*)mcArray->At(track->GetLabel())) : 0;
	lead[3] = 1./track->GetWeight(); //weight is 1./efficiency
        highPt = track->Pt();
      }

    } // end of tracks loop
  } //end of event loop for fCorrelatorTr

  if(fMixing) {
    NofEventsinPool = fCorrelatorKc->GetNofEventsInPool(); 
    if(!execPoolKc) {
      AliInfo("Mixed event analysis: K+/- pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //Charged Kaons loop
  for (Int_t jMix = 0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksKc = fCorrelatorKc->ProcessAssociatedTracks(jMix);
    if(!analyzetracksKc) {
      AliInfo("AliHFCorrelator::Cannot process the K+/- array");
      continue;
    }  

    for(Int_t iTrack = 0; iTrack<fCorrelatorKc->GetNofTracks(); iTrack++){ // looping on charged kaons candidates

      Bool_t runcorrelation = fCorrelatorKc->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* kCharg = fCorrelatorKc->GetAssociatedParticle();

      if(!fMixing) {      
  	if(!kCharg->CheckSoftPi()) { //removal of soft pions
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(1.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(1.,mD0bar);
          continue;
        } else {
          if (fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(0.,mD0);
          if (fIsSelectedCandidate >= 2) ((TH2F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(0.,mD0bar);
        }
        Int_t idDaughs[2] = {((AliVTrack*)d->GetDaughter(0))->GetID(),((AliVTrack*)d->GetDaughter(1))->GetID()}; //IDs of daughters to be skipped
        if(kCharg->GetID() == idDaughs[0] || kCharg->GetID() == idDaughs[1]) continue; //discards daughters of candidate
      }
      if(kCharg->Pt() < fPtThreshLow.at(ptbin) || kCharg->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K

      FillSparsePlots(mcArray,mInv,origD0,PDGD0,kCharg,ptbin,kKCharg); //fills for charged tracks

      if(!fMixing) N_KCharg++;

    } // end of charged kaons loop
  } //end of event loop for fCorrelatorKc

  if(fMixing) {
    NofEventsinPool = fCorrelatorK0->GetNofEventsInPool(); 
    if(!execPoolK0) {
      AliInfo("Mixed event analysis: K0 pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //K0 loop
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksK0 = fCorrelatorK0->ProcessAssociatedTracks(jMix);
    if(!analyzetracksK0) {
      AliInfo("AliHFCorrelator::Cannot process the K0 array");
      continue;
    }  

    for(Int_t iTrack = 0; iTrack<fCorrelatorK0->GetNofTracks(); iTrack++){ // looping on k0 candidates

      Bool_t runcorrelation = fCorrelatorK0->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* k0 = fCorrelatorK0->GetAssociatedParticle();

      if(k0->Pt() < fPtThreshLow.at(ptbin) || k0->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
  
      FillSparsePlots(mcArray,mInv,origD0,PDGD0,k0,ptbin,kK0); //fills for charged tracks

      if(!fMixing) N_Kaons++;

    } // end of charged kaons loop
  } //end of event loop for fCorrelatorK0

  Double_t fillSpLeadD0[3] = {lead[0],mD0,lead[1]};
  Double_t fillSpLeadD0bar[3] = {lead[0],mD0bar,lead[1]};

  //leading track correlations fill
  if(!fMixing) {
    if(fReadMC) {
      if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==421) { //D0
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadD0,lead[3]); //c and b D0
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0,lead[3]); //c or b D0
        if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0,lead[3]);  
        if(origD0==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0,lead[3]);  
        if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadD0,lead[3]);  //non HF  
      }
      if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==-421) { //D0bar
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadD0bar,lead[3]);
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0bar,lead[3]); //c or b D0
        if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0bar,lead[3]);  
        if(origD0==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadD0bar,lead[3]); 
        if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadD0bar,lead[3]);  //non HF  
      }
    } else {
        if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadD0,lead[3]); 
        if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadD0bar,lead[3]);
    }

    //Fill of count histograms
    if (!fAlreadyFilled) { 
      ((TH1F*)fOutputStudy->FindObject("hist_Count_Charg"))->Fill(N_Charg);
      ((TH1F*)fOutputStudy->FindObject("hist_Count_Kcharg"))->Fill(N_KCharg);
      ((TH1F*)fOutputStudy->FindObject("hist_Count_K0"))->Fill(N_Kaons);
    }
  }

  fAlreadyFilled=kTRUE; //at least a D0 analyzed in the event; distribution plots already filled

}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::CalculateCorrelationsMCKine(AliAODMCParticle* d, TClonesArray* mcArray) {
//
// Method for correlations D0-hadrons study
//
  Int_t N_Charg = 0, N_KCharg = 0, N_Kaons = 0;
  Double_t mD0 = 1.864, mD0bar = 1.864;
  Double_t mInv[2] = {mD0, mD0bar};
  Int_t origD0 = 0, PDGD0 = 0;
  Int_t ptbin = PtBinCorr(d->Pt());

  if(ptbin < 0) return;

  //Fill of D0 phi distribution
  if (!fMixing) ((TH1F*)fOutputStudy->FindObject("hist_PhiDistr_D0"))->Fill(d->Phi()); 
  
  //Origin of D0
  TString orig="";
  origD0=CheckD0Origin(mcArray,d);
  PDGD0 = d->GetPdgCode();
  switch (CheckD0Origin(mcArray,d)) {
    case 4:
      orig = "_From_c";
      ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(0.);
      break;
    case 5:
      orig = "_From_b";
      ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(1.);
      break;
    default:
      return;
  }

  Double_t highPt = 0; Double_t lead[3] = {0,0,0};  //infos for leading particle (pt,deltaphi)

  //loop over the tracks in the pool 
  Bool_t execPoolTr = fCorrelatorTr->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
  Bool_t execPoolKc = fCorrelatorKc->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
  Bool_t execPoolK0 = fCorrelatorK0->ProcessEventPool(); //pool is ready? (only in ME, in SE returns kFALSE)
		
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

      if (!trkMC->IsPhysicalPrimary()) continue; //reject material budget, or other fake tracks
      if (IsDDaughter(d,trkMC,mcArray)) continue;
      FillSparsePlots(mcArray,mInv,origD0,PDGD0,track,ptbin,kTrack); //fills for charged tracks

      //retrieving leading info...
      if(track->Pt() > highPt) {
        lead[0] = fCorrelatorTr->GetDeltaPhi();
        lead[1] = fCorrelatorTr->GetDeltaEta();
        lead[2] = fReadMC ? CheckTrackOrigin(mcArray,trkMC) : 0;
        highPt = track->Pt();
      }

    } // end of tracks loop
  } //end of event loop for fCorrelatorTr

  if(fMixing) {
    NofEventsinPool = fCorrelatorKc->GetNofEventsInPool(); 
    if(!execPoolKc) {
      AliInfo("Mixed event analysis: K+/- pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //Charged Kaons loop
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksKc = fCorrelatorKc->ProcessAssociatedTracks(jMix);
    if(!analyzetracksKc) {
      AliInfo("AliHFCorrelator::Cannot process the K+/- array");
      continue;
    }  

    for(Int_t iTrack = 0; iTrack<fCorrelatorKc->GetNofTracks(); iTrack++){ // looping on charged kaons candidates

      Bool_t runcorrelation = fCorrelatorKc->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* kCharg = fCorrelatorKc->GetAssociatedParticle();
      if(kCharg->GetLabel()<1) continue;
      if(kCharg->Pt() < fPtThreshLow.at(ptbin) || kCharg->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
      if(TMath::Abs(kCharg->Eta())>0.8) continue; //discard tracks outside barrel (since it's kinematic MC and produces tracks all over rapidity region
      if(!fMixing) N_KCharg++;

      AliAODMCParticle *kChargMC = (AliAODMCParticle*)mcArray->At(kCharg->GetLabel());
      if(!kChargMC) continue;

      if (IsDDaughter(d,kChargMC,mcArray)) continue;
      FillSparsePlots(mcArray,mInv,origD0,PDGD0,kCharg,ptbin,kKCharg); //fills for charged tracks

    } // end of charged kaons loop
  } //end of event loop for fCorrelatorKc

  if(fMixing) {
    NofEventsinPool = fCorrelatorK0->GetNofEventsInPool(); 
    if(!execPoolK0) {
      AliInfo("Mixed event analysis: K0 pool is not ready");
      NofEventsinPool = 0;
    }
  }

  //K0 loop
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++) {// loop on events in the pool; if it is SE analysis, stops at one (index not needed there)
    Bool_t analyzetracksK0 = fCorrelatorK0->ProcessAssociatedTracks(jMix);
    if(!analyzetracksK0) {
      AliInfo("AliHFCorrelator::Cannot process the K0 array");
      continue;
    }  

    for(Int_t iTrack = 0; iTrack<fCorrelatorK0->GetNofTracks(); iTrack++){ // looping on k0 candidates

      Bool_t runcorrelation = fCorrelatorK0->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* k0 = fCorrelatorK0->GetAssociatedParticle();
      if(k0->GetLabel()<1) continue;
      if(k0->Pt() < fPtThreshLow.at(ptbin) || k0->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
      if(TMath::Abs(k0->Eta())>0.8) continue; //discard tracks outside barrel (since it's kinematic MC and produces tracks all over rapidity region
  
      AliAODMCParticle *k0MC = (AliAODMCParticle*)mcArray->At(k0->GetLabel());
      if(!k0MC) continue;

      if (IsDDaughter(d,k0MC,mcArray)) continue;
      FillSparsePlots(mcArray,mInv,origD0,PDGD0,k0,ptbin,kK0); //fills for charged tracks

      if(!fMixing) N_Kaons++;

    } // end of charged kaons loop
  } //end of event loop for fCorrelatorK0

  Double_t fillSpLeadMC[3] = {lead[0],mD0,lead[1]}; //mD0 = mD0bar = 1.864

  //leading track correlations fill
  if(!fMixing) {
    if(d->GetPdgCode()==421) { //D0
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadMC); //c and b D0
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); //c or b D0
      if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if(origD0==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadMC);  //non HF
    }
    if(d->GetPdgCode()==-421) { //D0bar
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(fillSpLeadMC);
      ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); //c or b D0
      if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC);  
      if(origD0==5&&(int)lead[2]>=4&&(int)lead[2]<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpLeadMC); 
      if((int)lead[2]==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Lead_NonHF_Bin%d",ptbin)))->Fill(fillSpLeadMC);  //non HF
    }

    //Fill of count histograms
    if (!fAlreadyFilled) { 
      ((TH1F*)fOutputStudy->FindObject("hist_Count_Charg"))->Fill(N_Charg);
      ((TH1F*)fOutputStudy->FindObject("hist_Count_Kcharg"))->Fill(N_KCharg);
      ((TH1F*)fOutputStudy->FindObject("hist_Count_K0"))->Fill(N_Kaons);
    }

    fAlreadyFilled=kTRUE; //at least a D0 analyzed in the event; distribution plots already filled
  }

}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::FillSparsePlots(TClonesArray* mcArray, Double_t mInv[], Int_t origD0, Int_t PdgD0, AliReducedParticle* track, Int_t ptbin, Int_t type, Double_t wg) {
  //
  //fills the THnSparse for correlations, calculating the variables
  //

  //Initialization of variables
  Double_t mD0, mD0bar, deltaphi = 0., deltaeta = 0.;
  mD0 = mInv[0];
  mD0bar = mInv[1];

  if (fReadMC && track->GetLabel()<1) return;
  if (fReadMC && !(AliAODMCParticle*)mcArray->At(track->GetLabel())) return;
  Double_t ptTrack = track->Pt();
  Double_t d0Track = type!=kK0 ? track->GetImpPar() : 0.;
  Double_t phiTr = track->Phi();
  Double_t origTr = fReadMC ? CheckTrackOrigin(mcArray,(AliAODMCParticle*)mcArray->At(track->GetLabel())) : 0;

  TString part = "", orig = "";

  switch (type) {
    case(kTrack): {
      part = "Charg";
      deltaphi = fCorrelatorTr->GetDeltaPhi();
      deltaeta = fCorrelatorTr->GetDeltaEta();
      break;
    }
    case(kKCharg): {
      part = "Kcharg";
      deltaphi = fCorrelatorKc->GetDeltaPhi();
      deltaeta = fCorrelatorKc->GetDeltaEta();
      break;
    }
    case(kK0): {
      part = "K0";
      deltaphi = fCorrelatorK0->GetDeltaPhi();
      deltaeta = fCorrelatorK0->GetDeltaEta();
      break;
    }
  }
  
  if(fMixing == kSE) {

    //Fixes limits; needed to include overflow into THnSparse projections!
    Double_t pTorig = track->Pt();
    Double_t d0orig = track->GetImpPar();
    Double_t ptLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->GetAxis(2)->GetXmax(); //all plots have same axes...
    Double_t displLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->GetAxis(3)->GetXmax();
    Double_t EtaLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->GetAxis(4)->GetXmax();
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01;
    if(d0Track > displLim_Sparse) d0Track = (displLim_Sparse-0.001);
    if(deltaeta > EtaLim_Sparse) deltaeta = EtaLim_Sparse-0.01;
    if(deltaeta < -EtaLim_Sparse) deltaeta = -EtaLim_Sparse+0.01;
  
    //variables for filling histos
    Double_t fillSpPhiD0[5] = {deltaphi,mD0,ptTrack,d0Track,deltaeta};
    Double_t fillSpPhiD0bar[5] = {deltaphi,mD0bar,ptTrack,d0Track,deltaeta};
    Double_t fillSpWeigD0[3] = {deltaphi,mD0,deltaeta};
    Double_t fillSpWeigD0bar[3] = {deltaphi,mD0bar,deltaeta};

    if(fReadMC == 0) {
      //sparse fill for data (tracks, K+-, K0) + weighted
      if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) { //D0
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d",part.Data(),ptbin)))->Fill(fillSpPhiD0,wg);
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
      }
      if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) { //D0bar
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d",part.Data(),ptbin)))->Fill(fillSpPhiD0bar,wg);
        ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
      }
      if(!fAlreadyFilled) {
 	((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_%s_Bin%d",part.Data(),ptbin)))->Fill(pTorig);
 	((TH1F*)fOutputStudy->FindObject(Form("hist_PhiDistr_%s",part.Data())))->Fill(phiTr);
      }
    }

    if(fReadMC) {

      if(origD0==4) {orig = "_From_c";} else {orig = "_From_b";}

      //sparse fill for data (tracks, K+-, K0) + weighted
      if(PdgD0==421) { //D0 (from MCTruth)
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d",part.Data(),ptbin)))->Fill(fillSpPhiD0,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(fillSpPhiD0,wg);
         if(origD0==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(fillSpPhiD0,wg);
         if(origD0==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(fillSpPhiD0,wg);
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_NonHF_Bin%d",part.Data(),ptbin)))->Fill(fillSpPhiD0,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
         if(origD0==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
         if(origD0==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_NonHF_Bin%d",ptbin)))->Fill(fillSpWeigD0,pTorig*wg);
      }
      if(PdgD0==-421) { //D0bar
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d",part.Data(),ptbin)))->Fill(fillSpPhiD0bar,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(fillSpPhiD0bar,wg);
         if(origD0==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(fillSpPhiD0bar,wg);
         if(origD0==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_HF%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(fillSpPhiD0bar,wg); 
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_NonHF_Bin%d",part.Data(),ptbin)))->Fill(fillSpPhiD0bar,wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
         ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
         if(origD0==4&&origTr>=1&&origTr<=3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
         if(origD0==5&&origTr>=4&&origTr<=8) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
         if(origTr==0) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Weig_NonHF_Bin%d",ptbin)))->Fill(fillSpWeigD0bar,pTorig*wg);
      } 
      if(!fAlreadyFilled) {
	((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_Bin%d",part.Data(),ptbin)))->Fill(d0orig); //Fills displacement histos
        if (origTr>=1&&origTr<=8) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_HF_Bin%d",part.Data(),ptbin)))->Fill(d0orig);
        if (origTr>=1&&origTr<=8) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s_HF%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(d0orig);
	((TH1F*)fOutputStudy->FindObject(Form("histDispl_%s%s_Bin%d",part.Data(),orig.Data(),ptbin)))->Fill(d0orig); //Fills displacement histos
        ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_%s_Bin%d",part.Data(),ptbin)))->Fill(pTorig);
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_%s_Bin%d",part.Data(),ptbin)))->Fill(origTr);
 	((TH1F*)fOutputStudy->FindObject(Form("hist_PhiDistr_%s",part.Data())))->Fill(phiTr);
      }
    }//end MC case

  } //end of SE fill

  if(fMixing == kME) {

    //Fixes limits; needed to include overflow into THnSparse projections!
    Double_t EtaLim_Sparse = ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d_EvMix",ptbin)))->GetAxis(2)->GetXmax();
    if(deltaeta > EtaLim_Sparse) deltaeta = EtaLim_Sparse-0.01;
    if(deltaeta < -EtaLim_Sparse) deltaeta = -EtaLim_Sparse+0.01;

    //variables for filling histos
    Double_t fillSpPhiD0[3] = {deltaphi,mD0,deltaeta};
    Double_t fillSpPhiD0bar[3] = {deltaphi,mD0bar,deltaeta};

    if(fReadMC == 0) {
      //sparse fill for data (tracks, K+-, K0)
      if(fIsSelectedCandidate == 1||fIsSelectedCandidate == 3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_EvMix",part.Data(),ptbin)))->Fill(fillSpPhiD0,wg);
      if(fIsSelectedCandidate == 2||fIsSelectedCandidate == 3) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_EvMix",part.Data(),ptbin)))->Fill(fillSpPhiD0bar,wg);
    }
    if(fReadMC == 1) {
      //sparse fill for data (tracks, K+-, K0)
      if(PdgD0==421)  ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_EvMix",part.Data(),ptbin)))->Fill(fillSpPhiD0,wg);
      if(PdgD0==-421) ((THnSparseF*)fOutputCorr->FindObject(Form("hPhi_%s_Bin%d_EvMix",part.Data(),ptbin)))->Fill(fillSpPhiD0bar,wg);
    }//end MC case

  } //end of ME fill
  
  return;
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSED0Correlations::CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
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
  if(fDebug>2) printf("AliAnalysisTaskSED0Correlations::CheckTrkOrigin() \n");
	
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
Bool_t AliAnalysisTaskSED0Correlations::IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* mcArray) const {

  //Daughter removal in MCKine case
  Bool_t isDaughter = kFALSE;
  Int_t labelD0 = d->GetLabel();

  Int_t mother = track->GetMother();

  //Loop on the mothers to find the D0 label (it must be the trigger D0, not a generic D0!)
  while (mother > 0){
    AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother)); //it's the mother of the track!
    if (mcMoth){
      if (mcMoth->GetLabel() == labelD0) isDaughter = kTRUE;
      mother = mcMoth->GetMother(); //goes back by one
    } else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  return isDaughter;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSED0Correlations::PtBinCorr(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin=-1;
  if(pt<fBinLimsCorr.at(0)) return ptbin; //out of bounds
  
  Int_t i = 0;
  while(pt>fBinLimsCorr.at(i)) {ptbin=i; i++;}
  
  return ptbin;
}

//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskSED0Correlations::SelectV0(AliAODv0* v0, AliAODVertex *vtx, Int_t opt, Int_t idArrayV0[][2]) const
{
  //
  // Selection for K0 hypotheses
  // options: 1 = selects mass invariant about 3 sigma inside the peak + threshold of 0.3 GeV
  // 	      2 = no previous selections

  if(!fCutsTracks->IsKZeroSelected(v0,vtx)) return kFALSE;

  AliAODTrack *v0Daug1 = (AliAODTrack*)v0->GetDaughter(0);
  AliAODTrack *v0Daug2 = (AliAODTrack*)v0->GetDaughter(1);

  if(opt==1) { //additional cuts for correlations (V0 has to be closer than 3 sigma from K0 mass)
    if(TMath::Abs(v0->MassK0Short()-0.4976) > 3*0.004) return kFALSE;
  }

  //This part removes double counting for swapped tracks!
  Int_t i = 0;  //while loop (until the last-written entry pair of ID!
  while(idArrayV0[i][0]!=-2 && idArrayV0[i][1]!=-2) {
    if((v0Daug1->GetID()==idArrayV0[i][0] && v0Daug2->GetID()==idArrayV0[i][1])||
       (v0Daug1->GetID()==idArrayV0[i][1] && v0Daug2->GetID()==idArrayV0[i][0])) return kFALSE;
    i++;
  }
  idArrayV0[i][0]=v0Daug1->GetID();
  idArrayV0[i][1]=v0Daug2->GetID();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::PrintBinsAndLimits() {

  cout << "--------------------------\n";
  cout << "PtBins = " << fNPtBinsCorr << "\n";
  cout << "PtBin limits--------------\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << "Bin "<<i+1<<" = "<<fBinLimsCorr.at(i)<<" to "<<fBinLimsCorr.at(i+1)<<"\n";
  }
  cout << "\n--------------------------\n";
  cout << "PtBin tresh. tracks low---\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fPtThreshLow.at(i) << "; ";
  }
  cout << "PtBin tresh. tracks up----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fPtThreshUp.at(i) << "; ";
  }
  cout << "\n--------------------------\n";
  cout << "D0 Eta cut for Correl = "<<fEtaForCorrel<<"\n";
  cout << "--------------------------\n";
  cout << "MC Truth = "<<fReadMC<<" - MC Reco: Trk = "<<fRecoTr<<", D0 = "<<fRecoD0<<"\n";
  cout << "--------------------------\n";
  cout << "Sel of Event tpye (PP,GS,FE,...)= "<<fSelEvType<<"\n";
  cout << "--------------------------\n";
  cout << "Ev Mixing = "<<fMixing<<"\n";
  cout << "--------------------------\n";
}

