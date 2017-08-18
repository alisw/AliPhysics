/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Author: Friederike Bock, Mike Sas, Oton Vazquez Doce                    *
* Version 1.0                                                             *
*                                                                         *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on electromagnetic cocktail output
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisTaskLMeeCocktailMC.h"
#include "AliVParticle.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskLMeeCocktailMC)

//________________________________________________________________________
AliAnalysisTaskLMeeCocktailMC::AliAnalysisTaskLMeeCocktailMC(): AliAnalysisTaskSE(),
  fArr(0x0),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fOutputContainer(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fHistNEvents(NULL),
  fhwEffpT(NULL),
  fhwMultpT(NULL),
  fhwMultmT(NULL),
  fhwMultpT2(NULL),
  fhwMultmT2(NULL),
  fmee(NULL),
  fmee_wmult(NULL),
  fmee_orig(NULL),
  fmee_orig_wmult(NULL),
  fphi(NULL),
  frap(NULL),
  fpteevsmee(NULL),
  fpteevsmee_wmult(NULL),
  fpteevsmee_orig(NULL),
  fpteevsmee_orig_wmult(NULL),
  fd1origpt(-1),
  fd1origp(-1),
  fd1origeta(-1),
  fd1origphi(-1),
  fd2origpt(-1),
  fd2origp(-1),
  fd2origeta(-1),
  fd2origphi(-1),
  fd1pt(-1),
  fd1p(-1),
  fd1eta(-1),
  fd1phi(-1),
  fd2pt(-1),
  fd2p(-1),
  fd2eta(-1),
  fd2phi(-1),
  feeorigpt(-1),
  feeorigp(-1),
  feeorigm(-1),
  feeorigeta(-1),
  feeorigphi(-1),
  feeorigphiv(-1),
  feept(-1),
  feemt(-1),
  feep(-1),
  feem(-1),
  feeeta(-1),
  feephi(-1),
  feephiv(-1),
  fmotherpt(-1),
  fmothermt(-1),
  fmotherp(-1),
  fmotherm(-1),
  fmothereta(-1),
  fmotherphi(-1),
  fID(-1),
  fweight(-1),
  fwEffpT(-1),
  fwMultpT(-1),
  fwMultmT(-1),
  fwMultpT2(-1),
  fwMultmT2(-1),
  fpass(-1),
  fFileName(0),
  fFile(0),
  fFileNameW(0),
  fFileW(0),
  teeTTree(NULL),
  fIsMC(1),
  fMaxY(2),
  fMinPt(2),
  fWriteTTree(2),
  fcollisionSystem(2)
{
  
}

//________________________________________________________________________
AliAnalysisTaskLMeeCocktailMC::AliAnalysisTaskLMeeCocktailMC(const char *name):
  AliAnalysisTaskSE(name),
  fArr(0x0),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fOutputContainer(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fHistNEvents(NULL),
  fhwEffpT(NULL),
  fhwMultpT(NULL),
  fhwMultmT(NULL),
  fhwMultpT2(NULL),
  fhwMultmT2(NULL),
  fmee(NULL),
  fmee_wmult(NULL),
  fmee_orig(NULL),
  fmee_orig_wmult(NULL),
  fphi(NULL),
  frap(NULL),
  fpteevsmee(NULL),
  fpteevsmee_wmult(NULL),
  fpteevsmee_orig(NULL),
  fpteevsmee_orig_wmult(NULL),
  fd1origpt(-1),
  fd1origp(-1),
  fd1origeta(-1),
  fd1origphi(-1),
  fd2origpt(-1),
  fd2origp(-1),
  fd2origeta(-1),
  fd2origphi(-1),
  fd1pt(-1),
  fd1p(-1),
  fd1eta(-1),
  fd1phi(-1),
  fd2pt(-1),
  fd2p(-1),
  fd2eta(-1),
  fd2phi(-1),
  feeorigpt(-1),
  feeorigp(-1),
  feeorigm(-1),
  feeorigeta(-1),
  feeorigphi(-1),
  feeorigphiv(-1),
  feept(-1),
  feemt(-1),
  feep(-1),
  feem(-1),
  feeeta(-1),
  feephi(-1),
  feephiv(-1),
  fmotherpt(-1),
  fmothermt(-1),
  fmotherp(-1),
  fmotherm(-1),
  fmothereta(-1),
  fmotherphi(-1),
  fID(-1),
  fweight(-1),
  fwEffpT(-1),
  fwMultpT(-1),
  fwMultmT(-1),
  fwMultpT2(-1),
  fwMultmT2(-1),
  fpass(-1),
  fFileName(0),
  fFile(0),
  fFileNameW(0),
  fFileW(0),
  teeTTree(NULL),
  fIsMC(1),
  fMaxY(2),
  fMinPt(2),
  fWriteTTree(2),
  fcollisionSystem(2)

{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskLMeeCocktailMC::~AliAnalysisTaskLMeeCocktailMC()
{

}

//________________________________________________________________________
void AliAnalysisTaskLMeeCocktailMC::UserCreateOutputObjects(){
  
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer          = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  fHistNEvents                = new TH1F("NEvents", "NEvents", 1, 0, 1);
  fHistNEvents->Sumw2();
  fOutputContainer->Add(fHistNEvents);


  // Get Multiplicity weight file:
  fFileNameW = "$ALICE_PHYSICS/PWGDQ/dielectron/files/MultiplicityWeight_mT.root";
  fFileW = TFile::Open(fFileNameW.Data());
  if(!fFileW->IsOpen()){
   AliError(Form("Could not open Multiplicity weight file %s",fFileNameW.Data() ));
  }
  fhwEffpT = (TH1F*) fFileW->Get("fhwEffpT"); // histo: eff weight in function of pT.
  fhwMultpT = (TH1F*) fFileW->Get("fhwMultpT"); // histo: multiplicity weight in function of pT.
  fhwMultmT = (TH1F*) fFileW->Get("fhwMultmT"); // histo: multiplicity weight in function of mT.
  fhwMultpT2 = (TH1F*) fFileW->Get("fhwMultpT_upperlimit"); // histo: multiplicity weight in function of pT.
  fhwMultmT2 = (TH1F*) fFileW->Get("fhwMultmT_upperlimit"); // histo: multiplicity weight in function of mT.
  fOutputContainer->Add(fhwEffpT);
  fOutputContainer->Add(fhwMultpT);
  fOutputContainer->Add(fhwMultmT);
  fOutputContainer->Add(fhwMultpT2);
  fOutputContainer->Add(fhwMultmT2);

  // Prepare resolution file
  if(fcollisionSystem<=200) fFileName = "$ALICE_PHYSICS/PWGDQ/dielectron/files/histos_smear_brems_pp_LongP_fix.root";
  if(fcollisionSystem==300) fFileName = "$ALICE_PHYSICS/PWGDQ/dielectron/files/histos_smear_brems_pPb_LongP_fix.root";
  if(fcollisionSystem==400) fFileName = "$ALICE_PHYSICS/PWGDQ/dielectron/files/histos_smear_brems_PbPb_LongP_fix.root";
  fFile = TFile::Open(fFileName.Data());
  if(!fFile->IsOpen()){
   AliError(Form("Could not open file %s",fFileName.Data() ));
  }
  TObjArray* arr=0x0;
  arr = (TObjArray*) fFile->Get("ptSlices");
  if (!arr) printf("no resolution array set! using internal parameterization. \n");
  else      printf("using resolution array: \n");
  fArr=arr;

  // Define the tree
  teeTTree = new TTree("eeTTree","a simple TTree");
  teeTTree->SetDirectory(0); // This is to force a memory-resident Tree, and avoid errors. // ????? is this necessary? does it create memory problems?
  teeTTree->Branch("d1origpt",&fd1origpt,"fd1origpt/F");
  teeTTree->Branch("d1origp",&fd1origp,"fd1origp/F");
  teeTTree->Branch("d1origeta",&fd1origeta,"fd1origeta/F");
  teeTTree->Branch("d1origphi",&fd1origphi,"fd1origphi/F");
  teeTTree->Branch("d2origpt",&fd2origpt,"fd2origpt/F");
  teeTTree->Branch("d2origp",&fd2origp,"fd2origp/F");
  teeTTree->Branch("d2origeta",&fd2origeta,"fd2origeta/F");
  teeTTree->Branch("d2origphi",&fd2origphi,"fd2origphi/F");
  teeTTree->Branch("d1pt",&fd1pt,"fd1pt/F");
  teeTTree->Branch("d1p",&fd1p,"fd1p/F");
  teeTTree->Branch("d1eta",&fd1eta,"fd1eta/F");
  teeTTree->Branch("d1phi",&fd1phi,"fd1phi/F");
  teeTTree->Branch("d2pt",&fd2pt,"fd2pt/F");
  teeTTree->Branch("d2p",&fd2p,"fd2p/F");
  teeTTree->Branch("d2eta",&fd2eta,"fd2eta/F");
  teeTTree->Branch("d2phi",&fd2phi,"fd2phi/F");
  teeTTree->Branch("eeorigpt",&feeorigpt,"feeorigpt/F");
  teeTTree->Branch("eeorigp",&feeorigp,"feeorigp/F");
  teeTTree->Branch("eeorigm",&feeorigm,"feeorigm/F");
  teeTTree->Branch("eeorigeta",&feeorigeta,"feeorigeta/F");
  teeTTree->Branch("eeorigphi",&feeorigphi,"feeorigphi/F");
  teeTTree->Branch("eeorigphiv",&feeorigphiv,"feeorigphiv/F");
  teeTTree->Branch("eept",&feept,"feept/F");
  teeTTree->Branch("eemt",&feemt,"feemt/F");
  teeTTree->Branch("eep",&feep,"feep/F");
  teeTTree->Branch("eem",&feem,"feem/F");
  teeTTree->Branch("eeeta",&feeeta,"feeeta/F");
  teeTTree->Branch("eephi",&feephi,"feephi/F");
  teeTTree->Branch("eephiv",&feephiv,"feephiv/F");
  teeTTree->Branch("motherpt",&fmotherpt,"fmotherpt/F");
  teeTTree->Branch("mothermt",&fmothermt,"fmothermt/F");
  teeTTree->Branch("motherp",&fmotherp,"fmotherp/F");
  teeTTree->Branch("motherm",&fmotherm,"fmotherm/F");
  teeTTree->Branch("mothereta",&fmothereta,"fmothereta/F");
  teeTTree->Branch("motherphi",&fmotherphi,"fmotherphi/F");
  teeTTree->Branch("ID",&fID,"fID/I");
  teeTTree->Branch("weight",&fweight,"fweight/D");
  teeTTree->Branch("wEffpT",&fwEffpT,"fwEffpT/D");
  teeTTree->Branch("wMultpT",&fwMultpT,"fwMultpT/D");
  teeTTree->Branch("wMultmT",&fwMultmT,"fwMultmT/D");
  teeTTree->Branch("wMultpT2",&fwMultpT2,"fwMultpT2/D");
  teeTTree->Branch("wMultmT2",&fwMultmT2,"fwMultmT2/D");
  teeTTree->Branch("pass",&fpass,"fpass/O");

  fOutputContainer->Add(teeTTree);
  /*
   //to read this tree in a root session something like this has to be done:
   TDirectory *thedir = _file0->GetDirectory("LMeeCocktailMC");
   TList *thelist = (TList*)thedir->Get("LMeeCocktailMC_0.80");
   TTree *eeTTree = (TTree*)thelist->FindObject("eeTTree");
   //or simply
   TList* list =(TList*)_file0->Get("LMeeCocktailMC/LMeeCocktailMC_0.80");TTree* eeTTree=(TTree*)list.FindObject("eeTTree");
  */

  // Define the histograms
  // ---------------------
  //   "22" Gamma
  //   "111" Pi0
  //   "221" Eta
  //   "331" EtaPrim
  //   "223" Omega
  //   "211" Pi+
  //   "-211" Pi-
  //   "113" rho0
  //   "213" rho+
  //   "-213" rho-
  //   "333" phi
  //   "443" J/psi
  //   "1114" Delta-
  //   "2114" Delta0
  //   "2214" Delta+
  //   "2224" Delta++
  //   "3212" Sigma0
  const Int_t nInputParticles = 14;
  Int_t   fParticleList_local[] = {111,221,331,223,113,213,-213,333,443,1114,2114,2214,2224,3212};
  TString fParticleListNames_local[] = {"Pi0","Eta","EtaPrim","omega","rho0","rho+","rho-","phi","J/psi","Delta-","Delta0","Delta+","Delta++","Sigma0"};  
  fParticleList       = fParticleList_local;
  fParticleListNames  = fParticleListNames_local;
  fmee = new TH1F*[nInputParticles+1];
  fmee_wmult = new TH1F*[nInputParticles+1];
  fmee_orig = new TH1F*[nInputParticles+1];
  fmee_orig_wmult = new TH1F*[nInputParticles+1];
  fphi = new TH1F*[nInputParticles+1];
  frap = new TH1F*[nInputParticles+1];
  fpteevsmee = new TH2F*[nInputParticles];
  fpteevsmee_wmult = new TH2F*[nInputParticles];
  fpteevsmee_orig = new TH2F*[nInputParticles];
  fpteevsmee_orig_wmult = new TH2F*[nInputParticles];
 
  //Int_t   histBinM  = 2000;
  Int_t   histBinM  = 600; //400
  Float_t histMinM  = 0.;
  Float_t histMaxM  = 6.; //10.
  Int_t   histBinPt = 160; //80
  Float_t histMinPt = 0.;
  Float_t histMaxPt = 8.; //8.
  Int_t   histBinPhi = 240; //320;
  Float_t histMinPhi = 0.;
  Float_t histMaxPhi = TMath::TwoPi(); //3.2; 
  Int_t   histBinRap  = 240;
  Float_t histMinRap  = -1.2;
  Float_t histMaxRap  =  1.2;
  fmee[nInputParticles] = new TH1F("mee","mee",histBinM,histMinM,histMaxM);
  fmee[nInputParticles]->Sumw2();
  fOutputContainer->Add(fmee[nInputParticles]);
  fmee_wmult[nInputParticles] = new TH1F("mee_wmult","mee_wmult",histBinM,histMinM,histMaxM);
  fmee_wmult[nInputParticles]->Sumw2();
  fOutputContainer->Add(fmee_wmult[nInputParticles]);
  fmee_orig[nInputParticles] = new TH1F("mee_orig","mee_orig",histBinM,histMinM,histMaxM);
  fmee_orig[nInputParticles]->Sumw2();
  fOutputContainer->Add(fmee_orig[nInputParticles]);
  fmee_orig_wmult[nInputParticles] = new TH1F("mee_orig_wmult","mee_orig_wmult",histBinM,histMinM,histMaxM);
  fmee_orig_wmult[nInputParticles]->Sumw2();
  fOutputContainer->Add(fmee_orig_wmult[nInputParticles]);
  fphi[nInputParticles] = new TH1F("phi","phi",histBinPhi,histMinPhi,histMaxPhi);
  fphi[nInputParticles]->Sumw2();
  fOutputContainer->Add(fphi[nInputParticles]);
  frap[nInputParticles] = new TH1F("rap","rap",histBinRap,histMinRap,histMaxRap);
  frap[nInputParticles]->Sumw2();
  fOutputContainer->Add(frap[nInputParticles]);
  fpteevsmee[nInputParticles] = new TH2F("pteevsmee","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fpteevsmee[nInputParticles]->Sumw2();
  fOutputContainer->Add(fpteevsmee[nInputParticles]);
  fpteevsmee_wmult[nInputParticles] = new TH2F("pteevsmee_wmult","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fpteevsmee_wmult[nInputParticles]->Sumw2();
  fOutputContainer->Add(fpteevsmee_wmult[nInputParticles]);
  fpteevsmee_orig[nInputParticles] = new TH2F("pteevsmee_orig","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fpteevsmee_orig[nInputParticles]->Sumw2();
  fOutputContainer->Add(fpteevsmee_orig[nInputParticles]);
  fpteevsmee_orig_wmult[nInputParticles] = new TH2F("pteevsmee_orig_wmult","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fpteevsmee_orig_wmult[nInputParticles]->Sumw2();
  fOutputContainer->Add(fpteevsmee_orig_wmult[nInputParticles]);

  for(Int_t i=0; i<nInputParticles; i++){
   fmee[i] = new TH1F(Form("mee_%s",fParticleListNames[i].Data()),Form("mee_%s",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM);
   fmee[i]->Sumw2();
   fOutputContainer->Add(fmee[i]);
   fmee_wmult[i] = new TH1F(Form("mee_%s_wmult",fParticleListNames[i].Data()),Form("mee_%s_wmult",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM);
   fmee_wmult[i]->Sumw2();
   fOutputContainer->Add(fmee_wmult[i]);
   fmee_orig[i] = new TH1F(Form("mee_%s_orig",fParticleListNames[i].Data()),Form("mee_%s_orig",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM);
   fmee_orig[i]->Sumw2();
   fOutputContainer->Add(fmee_orig[i]);
   fmee_orig_wmult[i] = new TH1F(Form("mee_%s_orig_wmult",fParticleListNames[i].Data()),Form("mee_%s_orig_wmult",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM);
   fmee_orig_wmult[i]->Sumw2();
   fOutputContainer->Add(fmee_orig_wmult[i]);
   fphi[i] = new TH1F(Form("phi_%s",fParticleListNames[i].Data()),Form("phi_%s",fParticleListNames[i].Data()),histBinPhi,histMinPhi,histMaxPhi);
   fphi[i]->Sumw2();
   fOutputContainer->Add(fphi[i]);
   frap[i] = new TH1F(Form("rap_%s",fParticleListNames[i].Data()),Form("rap_%s",fParticleListNames[i].Data()),histBinRap,histMinRap,histMaxRap);
   frap[i]->Sumw2();
   fOutputContainer->Add(frap[i]);
   fpteevsmee[i] = new TH2F(Form("pteevsmee_%s",fParticleListNames[i].Data()),Form("%s;#it{m}_{ee};#it{p}_{T,ee}",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
   fpteevsmee[i]->Sumw2();
   fOutputContainer->Add(fpteevsmee[i]);
   fpteevsmee_wmult[i] = new TH2F(Form("pteevsmee_%s_wmult",fParticleListNames[i].Data()),Form("%s;#it{m}_{ee};#it{p}_{T,ee}",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
   fpteevsmee_wmult[i]->Sumw2();
   fOutputContainer->Add(fpteevsmee_wmult[i]);
   fpteevsmee_orig[i] = new TH2F(Form("pteevsmee_%s_orig",fParticleListNames[i].Data()),Form("%s;#it{m}_{ee};#it{p}_{T,ee}",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
   fpteevsmee_orig[i]->Sumw2();
   fOutputContainer->Add(fpteevsmee_orig[i]);
   fpteevsmee_orig_wmult[i] = new TH2F(Form("pteevsmee_%s_orig_wmult",fParticleListNames[i].Data()),Form("%s;#it{m}_{ee};#it{p}_{T,ee}",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
   fpteevsmee_orig_wmult[i]->Sumw2();
   fOutputContainer->Add(fpteevsmee_orig_wmult[i]);
  }
  
  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskLMeeCocktailMC::UserExec(Option_t *)
{

  fInputEvent = InputEvent();
  //cout << "I found an Event" << endl;
  
  fMCEvent = MCEvent();
  if(fMCEvent == NULL) fIsMC = 0;
  
  if (fIsMC==0) return;
  //cout << "I found an MC header" << endl;
    
  fMCStack = fMCEvent->Stack();
  if(fMCStack == NULL) fIsMC = 0;
  if (fIsMC==0) return;
  
  fHistNEvents->Fill(0.5);
  //   cout << "the stack is intact" << endl;
  ProcessMCParticles();

  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliAnalysisTaskLMeeCocktailMC::ProcessMCParticles(){

  Int_t Skip2ndLeg=0;
  // Loop over all primary MC particle  
  for(UInt_t i = 0; i < fMCStack->GetNtrack(); i++) {

   if(i!=0&&i==Skip2ndLeg) continue; //skip if maked as second electron

    //get the particle
    TParticle* particle         = NULL;
    TParticle* particle2         = NULL;
    particle                    = (TParticle *)fMCStack->Particle(i);
    if (!particle) continue;

    //get the mother
    Bool_t hasMother            = kFALSE;
    Bool_t particleIsPrimary    = kTRUE;
    //cout << i << "\t"<< particle->GetMother(0) << endl;
    if (particle->GetMother(0)>-1){
      hasMother = kTRUE;
      particleIsPrimary = kFALSE;
    }
    TParticle* motherParticle   = NULL;
    if( hasMother ) motherParticle = (TParticle *)fMCStack->Particle(particle->GetMother(0));
    if (motherParticle){
      hasMother                 = kTRUE;
    }else{
      hasMother                 = kFALSE;
    }
    Bool_t motherIsPrimary    = kFALSE;

    if(hasMother){
     //if(motherParticle->GetMother(0)>-1)motherIsPrimary = kTRUE;
     if(motherParticle->GetMother(0)==-1)motherIsPrimary = kTRUE;
     
     //skip for the moment other particles rather than pi0, eta, etaprime, omega, rho, phi.
     switch(motherParticle->GetPdgCode()){
      case 111:
       break;
      case 221:
       break;
      case 331:
       break;
      case 223:
       break;
      case 113:
       break;
      case 333:
       break;
      default:
       continue;
      }
    }

    // Not sure about this cut. From GammaConv group
    if (!(fabs(particle->Energy()-particle->Pz())>0.)) continue;

    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
    if (yPre == 0.) continue;
    Double_t y = 0.5*TMath::Log(yPre); 
    //if (fabs(y) > fMaxY) continue; // this cut done later in acceptance cuts
    
    // We have an electron with a mother. Check mother that mother is primary and number of daughters
    if(abs(particle->GetPdgCode())==11 && hasMother==kTRUE){
     UInt_t dectyp = 0;
     dectyp=motherParticle->GetDaughter(1)-motherParticle->GetDaughter(0)+1;
     if(dectyp > 4) continue;  // exclude five or more particles decay
     if(motherIsPrimary){

     //check second leg
     if(i<fMCStack->GetNtrack()-1) {
      particle2  = (TParticle *)fMCStack->Particle(i+1);
      if(particle2->GetMother(0)==particle->GetMother(0)&&particle->GetMother(1)==-1&&particle2->GetMother(1)==-1){
       Skip2ndLeg=i+1;

       TLorentzVector dau1,dau2,ee;
       dau1.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->Energy());
       dau2.SetPxPyPzE(particle2->Px(),particle2->Py(),particle2->Pz(),particle2->Energy());
  
        //create dielectron:
       ee=dau1+dau2;
 
       //get index for histograms
       Int_t hindex;
       switch(motherParticle->GetPdgCode()){
        case 111:
         hindex=0;
         break;
        case 221:
         hindex=1;
         break;
        case 331:
         hindex=2;
         break;
        case 223:
         hindex=3;
         break;
        case 113:
         hindex=4;
         break;
        case 213:
         hindex=5;
         break;
        case -213:
         hindex=6;
         break;
        case 333:
         hindex=7;
         break;
        case 443:
         hindex=8;
         break;
        case 1114:
         hindex=9;
         break;
        case 2114:
         hindex=10;
         break;
        case 2214:
         hindex=11;
         break;
        case 2224:
         hindex=12;
         break;
        case 3212:
         hindex=13;
         break;
        }

        // Fill tree words before resolution/acceptance
        fd1origpt=dau1.Pt();
        fd1origp=dau1.P();
        fd1origeta=dau1.Eta();
        fd1origphi=dau1.Phi();
        fd2origpt=dau2.Pt();
        fd2origp=dau2.P();
        fd2origeta=dau2.Eta();
        fd2origphi=dau2.Phi();
        feeorigpt=ee.Pt();
        feeorigp=ee.P();
        feeorigm=ee.M();
        feeorigeta=ee.Eta();
        feeorigphi=ee.Phi();
        if(particle->GetPdgCode()>0) { feeorigphiv=PhiV(dau1,dau2); }else{ feeorigphiv=PhiV(dau2,dau1);};
        //get the efficiency weight
        Int_t effbin=fhwEffpT->FindBin(dau1.Pt());
        fwEffpT=fhwEffpT->GetBinContent(effbin);
        effbin=fhwEffpT->FindBin(dau2.Pt());
        fwEffpT=fwEffpT*fhwEffpT->GetBinContent(effbin);

        //Resolution and acceptance  
        dau1=ApplyResolution(dau1);
        dau2=ApplyResolution(dau2);
        fpass=kTRUE;
        if(dau1.Pt()<fMinPt||dau2.Pt()<fMinPt) fpass=kFALSE; //leg pT cut
        if(TMath::Abs(dau1.Rapidity())>fMaxY||TMath::Abs(dau2.Rapidity())>fMaxY) fpass=kFALSE;

        //Fill words after resolution/acceptance
        ee=dau1+dau2;
        fd1pt=dau1.Pt();
        fd1p=dau1.P();
        fd1eta=dau1.Eta();
        fd1phi=dau1.Phi();
        fd2pt=dau2.Pt();
        fd2p=dau2.P();
        fd2eta=dau2.Eta();
        fd2phi=dau2.Phi();
        feept=ee.Pt();
        feemt=ee.Mt();
        feep=ee.P();
        feem=ee.M();
        feeeta=ee.Eta();
        feephi=ee.Phi();
        if(particle->GetPdgCode()>0) { feephiv=PhiV(dau1,dau2); }else{ feephiv=PhiV(dau2,dau1);};
        fmotherpt=motherParticle->Pt();
        fmothermt=sqrt(pow(motherParticle->GetCalcMass(),2)+pow(motherParticle->Pt(),2));
        fmotherp=motherParticle->P();
        fmotherm=motherParticle->GetCalcMass();
        fmothereta=motherParticle->Eta();
        fmotherphi=motherParticle->Phi();
        fID=1000*dectyp+motherParticle->GetPdgCode();
        fweight=particle->GetWeight(); //get particle weight from generator
        //now get multiplicity based weight:
        int iwbin=fhwMultpT->FindBin(fmotherpt);
        fwMultpT=fhwMultpT->GetBinContent(iwbin); //pT weight
        fwMultpT2=fhwMultpT2->GetBinContent(iwbin); //pT weight
        double min_mT=fhwMultmT->GetBinLowEdge(1); // consider as minimum valid mT value the edge of the weight histo.
        if(fmothermt>min_mT){
         iwbin=fhwMultmT->FindBin(fmothermt);
         fwMultmT = fhwMultmT->GetBinContent(iwbin); //mT weight
         fwMultmT2 = fhwMultmT2->GetBinContent(iwbin); //mT weight
        }else{
         printf("AliAnalysisTaskLMeeCocktailMC ERROR = Generated particle with mT < Pion mass cannot be weighted \n");
         fwMultmT = 0.;
         fwMultmT2 = 0.;
        }
 
        //Fill the tree
        if(fWriteTTree) teeTTree->Fill();

        //Fill the histograms
        if(dectyp<4){
         fmee_orig[hindex]->Fill(ee.M(), particle->GetWeight()*fwEffpT);
         fmee_orig_wmult[hindex]->Fill(ee.M(), particle->GetWeight()*fwEffpT*fwMultmT);
         fpteevsmee_orig[hindex]->Fill(ee.M(),ee.Pt(), particle->GetWeight()*fwEffpT);
         fpteevsmee_orig_wmult[hindex]->Fill(ee.M(),ee.Pt(), particle->GetWeight()*fwEffpT*fwMultmT);
         const Int_t nInputParticles = 14;
         fmee_orig[nInputParticles]->Fill(ee.M(), particle->GetWeight()*fwEffpT);
         fmee_orig_wmult[nInputParticles]->Fill(ee.M(), particle->GetWeight()*fwEffpT*fwMultmT);
         fpteevsmee_orig[nInputParticles]->Fill(ee.M(),ee.Pt(), particle->GetWeight()*fwEffpT);
         fpteevsmee_orig_wmult[nInputParticles]->Fill(ee.M(),ee.Pt(), particle->GetWeight()*fwEffpT*fwMultmT);
         if(fpass&&dectyp<4){
          fmee[hindex]->Fill(ee.M(), particle->GetWeight()*fwEffpT);
          fmee_wmult[hindex]->Fill(ee.M(), particle->GetWeight()*fwEffpT*fwMultmT);
          fphi[hindex]->Fill(ee.Phi(), particle->GetWeight()*fwEffpT);
          frap[hindex]->Fill(ee.Phi(), particle->GetWeight()*fwEffpT);
          fpteevsmee[hindex]->Fill(ee.M(),ee.Pt(), particle->GetWeight()*fwEffpT);
          fpteevsmee_wmult[hindex]->Fill(ee.M(),ee.Pt(), particle->GetWeight()*fwEffpT*fwMultmT);
          //const Int_t nInputParticles = 14;
          fmee[nInputParticles]->Fill(ee.M(), particle->GetWeight()*fwEffpT);
          fmee_wmult[nInputParticles]->Fill(ee.M(), particle->GetWeight()*fwEffpT*fwMultmT);
          fphi[nInputParticles]->Fill(ee.Phi(), particle->GetWeight()*fwEffpT);
          frap[nInputParticles]->Fill(ee.Phi(), particle->GetWeight()*fwEffpT);
          fpteevsmee[nInputParticles]->Fill(ee.M(),ee.Pt(), particle->GetWeight()*fwEffpT);
          fpteevsmee_wmult[nInputParticles]->Fill(ee.M(),ee.Pt(), particle->GetWeight()*fwEffpT*fwMultmT);
         }
        }

       } // legs coincide
      } // check the 2nd leg

     } // mother is primary
   } // pdgid==11 and HasMother
    
  }//MC particles loop
}

//________________________________________________________________________
void AliAnalysisTaskLMeeCocktailMC::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}


//_______________________________________________________________________________________________
TLorentzVector AliAnalysisTaskLMeeCocktailMC::ApplyResolution(TLorentzVector vec)
{
// from Theos LightFlavorGenerator, modified

  Double_t theta, phi, pt, p, px, py, pz, E, mass;
  TLorentzVector resvec;

  mass  = 0.51099906e-3;
  pt    = vec.Pt();
  p     = vec.P();
  theta = vec.Theta();
  phi   = vec.Phi();

  TH1D *hisSlice(0x0);
  if(fArr){
    TH2D *hDeltaPtvsPt = static_cast<TH2D*> (fArr->At(0));
    //std::cout << "2D momentum reso hist = " << hDeltaPtvsPt->GetTitle() << std::endl;
    // Get the momentum slice histogram for sampling of the smearing
    // (in some input file versions this also contains a langau fit )
    // since histogram bins start at 1, we can use the bin number directly for the array index (first slice stored in position 1).
    Int_t histIndex = hDeltaPtvsPt->GetXaxis()->FindBin(pt);
    if (histIndex<1) histIndex=1; // in case some track is below the first p-bin (which currently starts at 100 MeV).
    if (histIndex>fArr->GetLast()) {
      printf(" track with too high momentum (%f GeV), using last slice.\n", pt); // add some printouts when testing new input files!
      histIndex=fArr->GetLast();
    }
    hisSlice = static_cast<TH1D*> (fArr->At(histIndex));
  }
  // get smear parameter via random selection from the p slices retreived from the deltaP/P plot
  Double_t SmearingPPercent(0.);
  if(hisSlice) {
    //std::cout << "using momentum resolution from histogram." << std::endl;
    SmearingPPercent = hisSlice->GetRandom();
  } else {
    //std::cout << "using parameterized resolution." << std::endl;
    SmearingPPercent = gRandom->Gaus(0,p*sqrt(0.004*0.004 + (0.012*p)*(0.012*p))) / p ; // for B=0.5 Tesla
  }

  Double_t SmearingP = p * SmearingPPercent;
  p   -= SmearingP;
  px   = p*sin(theta)*cos(phi);
  py   = p*sin(theta)*sin(phi);
  pz   = p*cos(theta);
  E    = sqrt(p*p + mass*mass);

//  vec.SetPxPyPzE(px,py,pz,E);
  resvec.SetPxPyPzE(px,py,pz,E);
 // printf("         RESVEC P %f \n",resvec.P());
  return resvec;
}


Double_t AliAnalysisTaskLMeeCocktailMC::PhiV(TLorentzVector e1, TLorentzVector e2)
{
  Double_t outPhiV;
  TVector3 p1=e1.Vect();
  TVector3 p2=e2.Vect();
  TVector3 p12=p1+p2;
  TVector3 u=p12.Unit();
  TVector3 p1u=p1.Unit();
  TVector3 p2u=p2.Unit();
  TVector3 v=p1u.Cross(p2u);
  TVector3 w=u.Cross(v);
  TVector3 zu(0,0,1);
  TVector3 wc=u.Cross(zu);
  outPhiV=w.Angle(wc);
  return outPhiV;
}


