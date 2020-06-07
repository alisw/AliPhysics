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
#include "TF1.h"
#include "TH2F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TDatabasePDG.h"
#include "TGenPhaseSpace.h"
#include "TSystem.h"
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
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fHistNEvents(NULL),
  fhwEffpT(NULL),
  fhwMultpT(NULL),
  fhwMultmT(NULL),
  fhwMultpT2(NULL),
  fhwMultmT2(NULL),
  fh_DCAtemplates(NULL),
  fArr(0x0),
  fArrResoPt(0x0),
  fArrResoEta(0x0),
  fArrResoPhi_Pos(0x0),
  fArrResoPhi_Neg(0x0),
  ffVPHpT(0x0),
  fhKW(NULL),
  fmee_orig(NULL),
  fpteevsmee_orig(NULL),
  fmotherpT_orig(NULL),
  fphi_orig(NULL),
  frap_orig(NULL),
  fmee_orig_wALT(NULL),
  fpteevsmee_orig_wALT(NULL),
  fmotherpT_orig_wALT(NULL),
  fmee(NULL),
  fpteevsmee(NULL),
  fphi(NULL),
  frap(NULL),
  fDCAeevsmee(NULL),
  fDCAeevsptee(NULL),
  fmee_wALT(NULL),
  fpteevsmee_wALT(NULL),
  fULS_orig(NULL),
  fLSpp_orig(NULL),
  fLSmm_orig(NULL),
  fULS(NULL),
  fLSpp(NULL),
  fLSmm(NULL),
  fd1DCA(-1),
  fd2DCA(-1),
  fpairDCA(-1),
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
  fdectyp(-1),
  fdau3pdg(-1),
  fweight(-1),
  fwEffpT(-1),
  fwMultpT(-1),
  fwMultmT(-1),
  fwMultpT2(-1),
  fwMultmT2(-1),
  fpass(-1),
  fFileName(0),
  fFile(0),
  fFileNameDCA(0),
  fFileDCA(0),
  fFileNameEff(0),
  fFileEff(0),
  fFileNameVPH(0),
  fFileVPH(0),
  fResolDataSetName(""),
  teeTTree(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fIsMC(1),
  fMaxEta(2),
  fMinPt(0),
  fMaxPt(1000),
  fWriteTTree(2),
  fcollisionSystem(2),
  fResolType(2),
  fALTweightType(2)
{

}

//________________________________________________________________________
AliAnalysisTaskLMeeCocktailMC::AliAnalysisTaskLMeeCocktailMC(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fHistNEvents(NULL),
  fhwEffpT(NULL),
  fhwMultpT(NULL),
  fhwMultmT(NULL),
  fhwMultpT2(NULL),
  fhwMultmT2(NULL),
  fh_DCAtemplates(NULL),
  fArr(0x0),
  fArrResoPt(0x0),
  fArrResoEta(0x0),
  fArrResoPhi_Pos(0x0),
  fArrResoPhi_Neg(0x0),
  ffVPHpT(0x0),
  fhKW(NULL),
  fmee_orig(NULL),
  fpteevsmee_orig(NULL),
  fmotherpT_orig(NULL),
  fphi_orig(NULL),
  frap_orig(NULL),
  fmee_orig_wALT(NULL),
  fpteevsmee_orig_wALT(NULL),
  fmotherpT_orig_wALT(NULL),
  fmee(NULL),
  fpteevsmee(NULL),
  fphi(NULL),
  frap(NULL),
  fDCAeevsmee(NULL),
  fDCAeevsptee(NULL),
  fmee_wALT(NULL),
  fpteevsmee_wALT(NULL),
  fULS_orig(NULL),
  fLSpp_orig(NULL),
  fLSmm_orig(NULL),
  fULS(NULL),
  fLSpp(NULL),
  fLSmm(NULL),
  fd1DCA(-1),
  fd2DCA(-1),
  fpairDCA(-1),
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
  fdectyp(-1),
  fdau3pdg(-1),
  fweight(-1),
  fwEffpT(-1),
  fwMultpT(-1),
  fwMultmT(-1),
  fwMultpT2(-1),
  fwMultmT2(-1),
  fpass(-1),
  fFileName(0),
  fFile(0),
  fFileNameDCA(0),
  fFileDCA(0),
  fFileNameEff(0),
  fFileEff(0),
  fFileNameVPH(0),
  fFileVPH(0),
  fResolDataSetName(""),
  teeTTree(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fIsMC(1),
  fMaxEta(2),
  fMinPt(0),
  fMaxPt(1000),
  fWriteTTree(2),
  fcollisionSystem(2),
  fResolType(2),
  fALTweightType(2)
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

  // Get Efficiency (and Multiplicity) weight file:
  fFileNameEff = "$ALICE_PHYSICS/PWGDQ/dielectron/files/LMeeCocktailInputs_EffMult.root";
  fFileEff = TFile::Open(fFileNameEff.Data());
  if(!fFileEff->IsOpen()){
   AliError(Form("Could not open Efficiency and Multiplicity weight file %s",fFileNameEff.Data() ));
  }
  fhwEffpT = (TH1F*) fFileEff->Get("fhwEffpT"); // histo: eff weight in function of pT.
  fhwMultpT = (TH1F*) fFileEff->Get("fhwMultpT"); // histo: multiplicity weight in function of pT.
  fhwMultmT = (TH1F*) fFileEff->Get("fhwMultmT"); // histo: multiplicity weight in function of mT.
  fhwMultpT2 = (TH1F*) fFileEff->Get("fhwMultpT_upperlimit"); // histo: multiplicity weight in function of pT.
  fhwMultmT2 = (TH1F*) fFileEff->Get("fhwMultmT_upperlimit"); // histo: multiplicity weight in function of mT.
  // store those weights in the output file:
  //fOutputContainer->Add(fhwEffpT);
  //fOutputContainer->Add(fhwMultpT);
  //fOutputContainer->Add(fhwMultmT);
  //fOutputContainer->Add(fhwMultpT2);
  //fOutputContainer->Add(fhwMultmT2);

  // Get DCA templates:
  fFileNameDCA = "$ALICE_PHYSICS/PWGDQ/dielectron/files/LMeeCocktailInputs_DCA.root";
  fFileDCA = TFile::Open(fFileNameDCA.Data());
  if(!fFileDCA->IsOpen()){
   AliError(Form("Could not open DCA templates file %s",fFileNameDCA.Data() ));
  }
  fh_DCAtemplates = new TH1F*[6];
  fh_DCAtemplates[0] = (TH1F*) fFileDCA->Get("fh_DCAtemplate1"); // histo: DCA template 0.2 < pT < 0.3 GeV/c.
  fh_DCAtemplates[1] = (TH1F*) fFileDCA->Get("fh_DCAtemplate2"); // histo: DCA template 0.3 < pT < 0.4 GeV/c.
  fh_DCAtemplates[2] = (TH1F*) fFileDCA->Get("fh_DCAtemplate3"); // histo: DCA template 0.4 < pT < 0.6 GeV/c.
  fh_DCAtemplates[3] = (TH1F*) fFileDCA->Get("fh_DCAtemplate4"); // histo: DCA template 0.6 < pT < 1. GeV/c.
  fh_DCAtemplates[4] = (TH1F*) fFileDCA->Get("fh_DCAtemplate5"); // histo: DCA template 1. < pT < 2. GeV/c.
  fh_DCAtemplates[5] = (TH1F*) fFileDCA->Get("fh_DCAtemplate6"); // histo: DCA template pT > 2. GeV/c.
  //for(int ii=0;ii<6;ii++) fOutputContainer->Add(fh_DCAtemplates[ii]);

  //get the template for virtual photons pT (from now use 7TeV pi0 pT parametrization):
  fFileNameVPH = "$ALICE_PHYSICS/PWG/Cocktail/parametrisations/pp_7TeV.root";
  fFileVPH = TFile::Open(fFileNameVPH.Data());
  if(!fFileVPH->IsOpen()){
   AliError(Form("Could not open Virtual Photon templates file %s",fFileNameVPH.Data() ));
  }
  ffVPHpT = (TF1*)fFileVPH->GetDirectory("7TeV_Comb")->Get("111_pt");
  //Build Kroll-wada for virtual photon mass parametrization:
  Double_t KWmass = 0.;
  Double_t  emass = (TDatabasePDG::Instance()->GetParticle(11))->Mass();
  Int_t KWnbins = 10000;
  Float_t KWmin   = 2.*emass;
  Float_t KWmax         = 1.1;
  Double_t KWbinwidth     = (KWmax - KWmin) / (Double_t)KWnbins;
  fhKW = new TH1F("fhKW","fhKW",KWnbins,KWmin,KWmax);
  for(Int_t ibin = 1; ibin <= KWnbins; ibin++ ){
    KWmass     = KWmin + (Double_t)(ibin - 1) * KWbinwidth + KWbinwidth / 2.0;
    fhKW->AddBinContent(ibin,2.*(1./137.03599911)/3./3.14159265359/KWmass
       *sqrt(1.-4.*emass*emass/KWmass/KWmass)*(1.+2.*emass*emass/KWmass/KWmass));
  }

  // Prepare resolution file
  //RUN1
  if(fResolType == 1) {
   if(fcollisionSystem<=200) fFileName = "$ALICE_PHYSICS/PWGDQ/dielectron/files/LMeeCocktailInputs_Respp.root";
   if(fcollisionSystem==300) fFileName = "$ALICE_PHYSICS/PWGDQ/dielectron/files/LMeeCocktailInputs_RespPb.root";
   if(fcollisionSystem==400) fFileName = "$ALICE_PHYSICS/PWGDQ/dielectron/files/LMeeCocktailInputs_ResPbPb.root";
   if(fcollisionSystem<=200||fcollisionSystem==300||fcollisionSystem==400){
    fFile = TFile::Open(fFileName.Data());
    if(!fFile->IsOpen()){
     AliError(Form("Could not open file %s",fFileName.Data() ));
    }
    TObjArray* arr=0x0;
    arr = (TObjArray*) fFile->Get("ptSlices");
    if (!arr) printf("no resolution array set! using internal parameterization. \n");
    else      printf("using resolution array: \n");
    fArr=arr;
   }
  }
  //RUN2
  if(fResolType == 2) {
    if(fResolDataSetName.Contains("alien")){
      // file is copied from alien path to local directory
      gSystem->Exec(Form("alien_cp %s .", fResolDataSetName.Data()));
      
      // obtain ROOT file name only and local directory
      TObjArray* Strings = fResolDataSetName.Tokenize("/");
      fFileName = Form("%s/%s",gSystem->pwd(),Strings->At(Strings->GetEntriesFast()-1)->GetName());
      
      Printf("Set resolution file name to %s (copied from %s)",fFileName.Data(),fResolDataSetName.Data());
    }
    else{
      if(fcollisionSystem==200){ //pp 13TeV
	fFileName = "$ALICE_PHYSICS/PWGDQ/dielectron/files/LMeeCocktailInputs_Respp13TeV.root";
      }
      else{
	fFileName = "$ALICE_PHYSICS/PWGDQ/dielectron/files/"+ fResolDataSetName;
      }
    }
   fFile = TFile::Open(fFileName.Data());
   if(!fFile->IsOpen()){
     AliError(Form("Could not open file %s",fFileName.Data() ));
   }
    TObjArray* ArrResoPt=0x0;
    ArrResoPt = (TObjArray*) fFile->Get("RelPtResArrCocktail");
    TObjArray* ArrResoEta=0x0;
    ArrResoEta = (TObjArray*) fFile->Get("EtaResArr");
    TObjArray* ArrResoPhi_Pos=0x0;
    ArrResoPhi_Pos = (TObjArray*) fFile->Get("PhiPosResArr");
    TObjArray* ArrResoPhi_Neg=0x0;
    ArrResoPhi_Neg = (TObjArray*) fFile->Get("PhiEleResArr");
    fArrResoPt=ArrResoPt;
    fArrResoEta=ArrResoEta;
    fArrResoPhi_Pos=ArrResoPhi_Pos;
    fArrResoPhi_Neg=ArrResoPhi_Neg;
  } 

  // Define the output tree
  teeTTree = new TTree("eeTTree","a simple TTree");
  teeTTree->SetDirectory(0); // This is to force a memory-resident Tree, and avoid errors. // ????? is this necessary? does it create memory problems?
  teeTTree->Branch("d1DCA",&fd1DCA,"fd1DCA/F");
  teeTTree->Branch("d2DCA",&fd2DCA,"fd2DCA/F");
  teeTTree->Branch("pairDCA",&fpairDCA,"fpairDCA/F");
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
  teeTTree->Branch("dectyp",&fdectyp,"fdectyp/I");
  teeTTree->Branch("dau3pdg",&fdau3pdg,"fdau3pdg/I");
  teeTTree->Branch("weight",&fweight,"fweight/D");
  teeTTree->Branch("wEffpT",&fwEffpT,"fwEffpT/D");
  teeTTree->Branch("wMultpT",&fwMultpT,"fwMultpT/D");
  teeTTree->Branch("wMultmT",&fwMultmT,"fwMultmT/D");
  teeTTree->Branch("wMultpT2",&fwMultpT2,"fwMultpT2/D");
  teeTTree->Branch("wMultmT2",&fwMultmT2,"fwMultmT2/D");
  teeTTree->Branch("pass",&fpass,"fpass/O");

  fOutputContainer->Add(teeTTree);
  //to read such tree in a root session do:
  //TList* list =(TList*)_file0->Get("LMeeCocktailMC/LMeeCocktailMC_0.80");TTree* eeTTree=(TTree*)list.FindObject("eeTTree");


  // Define the histograms
  // ---------------------
  //  "111" Pi0 (0)
  //  "221" Eta (1)
  //XX      3221 // Eta_dalitz 
  //XX    114221 // eta 4-body -> e+e- e+e- 
  //XX   2114221 // eta 4-body -> e+e- pi+pi- 
  //  "331" EtaP (2)
  //      223331 // EtaP_dalitz_photon (3)
  //     2233331 // EtaP_dalitz_omega  (4)
  //XX      4331 // eta' 4-body 
  //  "113" Rho (5)
  //  "223" Omega (6)
  //       2223 // Omega_2body (7)
  //       3223 // Omega_dalitz (8)
  //  "333" Phi (9)
  //       2333 // Phi_2body (10)
  //    2213333 // Phi_dalitz_eta (11)
  //    1113333 // Phi_dalitz_pi0 (12)
  //XX  "443" Jpsi 
  // 

  //THE TOTAL NUMBER OF HISTOS is defined IN THE HEADER: const Int_t nInputParticles = 14; (#of particles)
  Int_t fParticleList_local[] = {111, 221, 331, 223331, 2233331, 113, 223, 2223, 3223, 333, 2333, 2213333, 1113333, 000};
  TString fParticleListNames_local[] = {"Pi0","Eta","EtaP","EtaP_dalitz_photon","EtaP_dalitz_omega","Rho","Omega","Omega_2body","Omega_dalitz","Phi","Phi_2body","Phi_dalitz_eta","Phi_dalitz_pi0","Virtual_Photon"};
  fParticleList       = fParticleList_local;
  fParticleListNames  = fParticleListNames_local;

  //booking
  Int_t   histBinM  = 1200; //600
  Float_t histMinM  = 0.;
  Float_t histMaxM  = 6.; //10.
  Int_t   histBinPt = 400; //160//80
  Float_t histMinPt = 0.;
  Float_t histMaxPt = 10.; //8.
  Int_t   histBinPhi = 240; //320;
  Float_t histMinPhi = 0.;
  Float_t histMaxPhi = TMath::TwoPi(); //3.2; 
  Int_t   histBinRap  = 240;
  Float_t histMinRap  = -1.2;
  Float_t histMaxRap  =  1.2;
  const Int_t nbm = 8;
  Double_t mbins[nbm+1] = {0.,0.08,0.14,0.2,1.1,2.7,2.8,3.2,5.0};
  const Int_t nbDCA=11;
  Double_t DCAbins[nbDCA+1]= {0.,0.4,0.8,1.2,1.6,2.0,2.4,3.,4.,5.,7.,10.};
  const Int_t nbpt=16;
  Double_t ptbins[nbpt+1]= {0.,0.5,1,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.};

  fmee_orig = new TH1F*[nInputParticles+1];
  fmee_orig[nInputParticles] = new TH1F("mee_orig","mee_orig",histBinM,histMinM,histMaxM);
  fmee_orig[nInputParticles]->Sumw2();
  fOutputContainer->Add(fmee_orig[nInputParticles]);
  fpteevsmee_orig = new TH2F*[nInputParticles+1];
  fpteevsmee_orig[nInputParticles] = new TH2F("pteevsmee_orig","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fpteevsmee_orig[nInputParticles]->Sumw2();
  fOutputContainer->Add(fpteevsmee_orig[nInputParticles]);
  fmotherpT_orig = new TH1F*[nInputParticles+1];
  fmotherpT_orig[nInputParticles] = new TH1F("motherpT_orig","motherpT_orig",histBinPt,histMinPt,histMaxPt);
  fmotherpT_orig[nInputParticles]->Sumw2();
  fOutputContainer->Add(fmotherpT_orig[nInputParticles]);
  fphi_orig = new TH1F*[nInputParticles+1];
  fphi_orig[nInputParticles] = new TH1F("phi_orig","phi_orig",histBinPhi,histMinPhi,histMaxPhi);
  fphi_orig[nInputParticles]->Sumw2();
  fOutputContainer->Add(fphi_orig[nInputParticles]);
  frap_orig = new TH1F*[nInputParticles+1];
  frap_orig[nInputParticles] = new TH1F("rap_orig","rap_orig",histBinRap,histMinRap,histMaxRap);
  frap_orig[nInputParticles]->Sumw2();
  fOutputContainer->Add(frap_orig[nInputParticles]);
  for(Int_t i=0; i<nInputParticles; i++){
   fmee_orig[i] = new TH1F(Form("mee_orig_%s",fParticleListNames[i].Data()),Form("mee_orig_%s",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM);
   fmee_orig[i]->Sumw2();
   fOutputContainer->Add(fmee_orig[i]);
   fphi_orig[i] = new TH1F(Form("phi_orig_%s",fParticleListNames[i].Data()),Form("phi_orig_%s",fParticleListNames[i].Data()),histBinPhi,histMinPhi,histMaxPhi);
   fphi_orig[i]->Sumw2();
   fOutputContainer->Add(fphi_orig[i]);
   frap_orig[i] = new TH1F(Form("rap_orig_%s",fParticleListNames[i].Data()),Form("rap_orig_%s",fParticleListNames[i].Data()),histBinRap,histMinRap,histMaxRap);
   frap_orig[i]->Sumw2();
   fOutputContainer->Add(frap_orig[i]);
   fpteevsmee_orig[i] = new TH2F(Form("pteevsmee_orig_%s",fParticleListNames[i].Data()),Form("%s;#it{m}_{ee};#it{p}_{T,ee}",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
   fpteevsmee_orig[i]->Sumw2();
   fOutputContainer->Add(fpteevsmee_orig[i]);
   fmotherpT_orig[i] = new TH1F(Form("motherpT_orig_%s",fParticleListNames[i].Data()),Form("motherpT_orig_%s",fParticleListNames[i].Data()),histBinPt,histMinPt,histMaxPt);
   fmotherpT_orig[i]->Sumw2();
   fOutputContainer->Add(fmotherpT_orig[i]);
  }

  fmee = new TH1F*[nInputParticles+1];
  fmee[nInputParticles] = new TH1F("mee","mee",histBinM,histMinM,histMaxM);
  fmee[nInputParticles]->Sumw2();
  fOutputContainer->Add(fmee[nInputParticles]);
  fpteevsmee = new TH2F*[nInputParticles+1];
  fpteevsmee[nInputParticles] = new TH2F("pteevsmee","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fpteevsmee[nInputParticles]->Sumw2();
  fOutputContainer->Add(fpteevsmee[nInputParticles]);
  fphi = new TH1F*[nInputParticles+1];
  fphi[nInputParticles] = new TH1F("phi","phi",histBinPhi,histMinPhi,histMaxPhi);
  fphi[nInputParticles]->Sumw2();
  fOutputContainer->Add(fphi[nInputParticles]);
  frap = new TH1F*[nInputParticles+1];
  frap[nInputParticles] = new TH1F("rap","rap",histBinRap,histMinRap,histMaxRap);
  frap[nInputParticles]->Sumw2();
  fOutputContainer->Add(frap[nInputParticles]);
  for(Int_t i=0; i<nInputParticles; i++){
   fmee[i] = new TH1F(Form("mee_%s",fParticleListNames[i].Data()),Form("mee_%s",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM);
   fmee[i]->Sumw2();
   fOutputContainer->Add(fmee[i]);
   fphi[i] = new TH1F(Form("phi_%s",fParticleListNames[i].Data()),Form("phi_%s",fParticleListNames[i].Data()),histBinPhi,histMinPhi,histMaxPhi);
   fphi[i]->Sumw2();
   fOutputContainer->Add(fphi[i]);
   frap[i] = new TH1F(Form("rap_%s",fParticleListNames[i].Data()),Form("rap_%s",fParticleListNames[i].Data()),histBinRap,histMinRap,histMaxRap);
   frap[i]->Sumw2();
   fOutputContainer->Add(frap[i]);
   fpteevsmee[i] = new TH2F(Form("pteevsmee_%s",fParticleListNames[i].Data()),Form("%s;#it{m}_{ee};#it{p}_{T,ee}",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
   fpteevsmee[i]->Sumw2();
   fOutputContainer->Add(fpteevsmee[i]);
  }
  fDCAeevsmee  = new TH2F("DCAeevsmee","DCAvsmee;#it{m}_{ee};DCA_{xy}^{ee}",nbm,mbins,nbDCA,DCAbins);
  fDCAeevsptee = new TH2F("DCAeevsptee","DCAvsptee;#it{p}_{T}^{ee};DCA_{xy}^{ee}",nbm,mbins,nbpt,ptbins);
  fDCAeevsmee->Sumw2();
  fDCAeevsptee->Sumw2();
  fOutputContainer->Add(fDCAeevsmee);
  fOutputContainer->Add(fDCAeevsptee);

  fmee_orig_wALT = new TH1F*[nInputParticles+1];
  fmee_orig_wALT[nInputParticles] = new TH1F("mee_orig_wALT","mee_orig_wALT",histBinM,histMinM,histMaxM);
  fmee_orig_wALT[nInputParticles]->Sumw2();
  if(fALTweightType>0)fOutputContainer->Add(fmee_orig_wALT[nInputParticles]);
  fpteevsmee_orig_wALT = new TH2F*[nInputParticles+1];
  fpteevsmee_orig_wALT[nInputParticles] = new TH2F("pteevsmee_orig_wALT","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fpteevsmee_orig_wALT[nInputParticles]->Sumw2();
  if(fALTweightType>0)fOutputContainer->Add(fpteevsmee_orig_wALT[nInputParticles]);
  fmotherpT_orig_wALT = new TH1F*[nInputParticles+1];
  fmotherpT_orig_wALT[nInputParticles] = new TH1F("motherpT_orig_wALT","motherpT_orig_wALT",histBinPt,histMinPt,histMaxPt);
  fmotherpT_orig_wALT[nInputParticles]->Sumw2();
  if(fALTweightType>0)fOutputContainer->Add(fmotherpT_orig_wALT[nInputParticles]);
  for(Int_t i=0; i<nInputParticles; i++){
   fmee_orig_wALT[i] = new TH1F(Form("mee_orig_wALT_%s",fParticleListNames[i].Data()),Form("mee_orig_wALT_%s",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM);
   fmee_orig_wALT[i]->Sumw2();
   if(fALTweightType>0)fOutputContainer->Add(fmee_orig_wALT[i]);
   fpteevsmee_orig_wALT[i] = new TH2F(Form("pteevsmee_orig_wALT_%s",fParticleListNames[i].Data()),Form("%s;#it{m}_{ee};#it{p}_{T,ee}",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
   fpteevsmee_orig_wALT[i]->Sumw2();
   if(fALTweightType>0)fOutputContainer->Add(fpteevsmee_orig_wALT[i]);
   fmotherpT_orig_wALT[i] = new TH1F(Form("motherpT_orig_wALT_%s",fParticleListNames[i].Data()),Form("motherpT_orig_wALT_%s",fParticleListNames[i].Data()),histBinPt,histMinPt,histMaxPt);
   fmotherpT_orig_wALT[i]->Sumw2();
   if(fALTweightType>0)fOutputContainer->Add(fmotherpT_orig_wALT[i]);
  }
  fmee_wALT = new TH1F*[nInputParticles+1];
  fmee_wALT[nInputParticles] = new TH1F("mee_wALT","mee_wALT",histBinM,histMinM,histMaxM);
  fmee_wALT[nInputParticles]->Sumw2();
  if(fALTweightType>0)fOutputContainer->Add(fmee_wALT[nInputParticles]);
  fpteevsmee_wALT = new TH2F*[nInputParticles+1];
  fpteevsmee_wALT[nInputParticles] = new TH2F("pteevsmee_wALT","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fpteevsmee_wALT[nInputParticles]->Sumw2();
  if(fALTweightType>0)fOutputContainer->Add(fpteevsmee_wALT[nInputParticles]);
  for(Int_t i=0; i<nInputParticles; i++){
   fmee_wALT[i] = new TH1F(Form("mee_wALT_%s",fParticleListNames[i].Data()),Form("mee_wALT_%s",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM);
   fmee_wALT[i]->Sumw2();
   if(fALTweightType>0)fOutputContainer->Add(fmee_wALT[i]);
   fpteevsmee_wALT[i] = new TH2F(Form("pteevsmee_wALT_%s",fParticleListNames[i].Data()),Form("%s;#it{m}_{ee};#it{p}_{T,ee}",fParticleListNames[i].Data()),histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
   fpteevsmee_wALT[i]->Sumw2();
   if(fALTweightType>0)fOutputContainer->Add(fpteevsmee_wALT[i]);
  }

  fULS_orig = new TH2F("ULS_orig","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fLSpp_orig = new TH2F("LSpp_orig","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fLSmm_orig = new TH2F("LSmm_orig","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fULS_orig->Sumw2();
  fLSpp_orig->Sumw2();
  fLSmm_orig->Sumw2();
  fOutputContainer->Add(fULS_orig);
  fOutputContainer->Add(fLSpp_orig);
  fOutputContainer->Add(fLSmm_orig);
  fULS = new TH2F("ULS","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fLSpp = new TH2F("LSpp","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fLSmm = new TH2F("LSmm","ptvsmee;#it{m}_{ee};#it{p}_{T,ee}",histBinM,histMinM,histMaxM,histBinPt,histMinPt,histMaxPt);
  fULS->Sumw2();
  fLSpp->Sumw2();
  fLSmm->Sumw2();
  fOutputContainer->Add(fULS);
  fOutputContainer->Add(fLSpp);
  fOutputContainer->Add(fLSmm);

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
    
  fHistNEvents->Fill(0.5);

  ProcessMCParticles();

  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliAnalysisTaskLMeeCocktailMC::ProcessMCParticles(){

  std::vector<TLorentzVector> eBuff;
  std::vector<Char_t> echBuff;
  std::vector<Double_t> eweightBuff;

  Int_t Skip2ndLeg=0;

  // Loop over all primary MC particle  
  for(UInt_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {

   //LS and ULS spectra
   //------------------
   if(abs(fMCEvent->GetTrack(i)->PdgCode())==11){
    //get the electron
    //---------------
    TLorentzVector e,dielectron;
    Char_t ech,dielectron_ch;
    Double_t eweight,dielectron_weight;
    e.SetPxPyPzE(fMCEvent->GetTrack(i)->Px(),fMCEvent->GetTrack(i)->Py(),fMCEvent->GetTrack(i)->Pz(),fMCEvent->GetTrack(i)->E());
    if(fMCEvent->GetTrack(i)->PdgCode()>0) { ech=1.; } else { ech=-1.; };
    eweight=fMCEvent->Particle(i)->GetWeight();
    //put in the buffer
    //-----------------
    eBuff.push_back(e);
    echBuff.push_back(ech);
    eweightBuff.push_back(eweight);
    //loop the buffer and pair
    //------------------------
    for(Int_t jj=eBuff.size()-2;jj>=0;jj--){
     dielectron=eBuff.at(jj)+e;
     dielectron_ch=(echBuff.at(jj)+ech)/2;
     dielectron_weight=eweightBuff.at(jj)*eweight;
     if(dielectron_ch==0) fULS_orig->Fill(dielectron.M(),dielectron.Pt(),dielectron_weight);
     if(dielectron_ch>0) fLSpp_orig->Fill(dielectron.M(),dielectron.Pt(),dielectron_weight);
     if(dielectron_ch<0) fLSmm_orig->Fill(dielectron.M(),dielectron.Pt(),dielectron_weight);
     if(e.Pt()>fMinPt&&eBuff.at(jj).Pt()>fMinPt&&e.Pt()<fMaxPt&&eBuff.at(jj).Pt()<fMaxPt&&TMath::Abs(e.Eta())<fMaxEta&&TMath::Abs(eBuff.at(jj).Eta())<fMaxEta){
      if(dielectron_ch==0) fULS->Fill(dielectron.M(),dielectron.Pt(),dielectron_weight);
      if(dielectron_ch>0) fLSpp->Fill(dielectron.M(),dielectron.Pt(),dielectron_weight);
      if(dielectron_ch<0) fLSmm->Fill(dielectron.M(),dielectron.Pt(),dielectron_weight);
     }
    }
   }


   if(i!=0&&i==Skip2ndLeg) continue; //skip if maked as second electron

    //get the particle
    AliVParticle* particle       = NULL;
    AliVParticle* particle2      = NULL;
    particle                    = (AliVParticle *)fMCEvent->GetTrack(i);
    if (!particle) continue;

    //get the mother
    Bool_t hasMother            = kFALSE;
    Bool_t particleIsPrimary    = kTRUE;
    //cout << i << "\t"<< particle->Particle(0) << endl;
    if (particle->Particle()->GetMother(0)>-1){
      hasMother = kTRUE;
      particleIsPrimary = kFALSE;
    }
    AliVParticle* motherParticle   = NULL;
    AliVParticle* dau3Particle   = NULL;
    fdau3pdg   = 0;
    if( hasMother ) motherParticle = (AliVParticle *)fMCEvent->GetTrack(particle->Particle()->GetMother(0));
    if (motherParticle){
      hasMother                 = kTRUE;
    }else{
      hasMother                 = kFALSE;
    }
    Bool_t motherIsPrimary    = kFALSE;

    if(hasMother){
     //if(motherParticle->Particle()->GetMother(0)>-1)motherIsPrimary = kTRUE;
     if(motherParticle->Particle()->GetMother(0)==-1)motherIsPrimary = kTRUE;
     
     //skip for the moment other particles rather than pi0, eta, etaprime, omega, rho, phi.
     switch(motherParticle->PdgCode()){
      case 111:
       break;
      case 221:
       break;
      case 331:
       break;
      case 113:
       break;
      case 223:
       break;
      case 333:
       break;
      default:
       continue;
      }
    }

    // Not sure about this cut. From GammaConv group. Harmless a priori.
    if (!(fabs(particle->E()-particle->Pz())>0.)) continue;

    Double_t yPre = (particle->E()+particle->Pz())/(particle->E()-particle->Pz());
    if (yPre == 0.) continue;
    
    // We have an electron with a mother. Check that mother is primary and number of daughters
    if(abs(particle->PdgCode())==11 && hasMother==kTRUE){
     fdectyp = 0; // fdectyp: decay type (based on number of daughters).
     fdectyp=motherParticle->GetDaughterLabel(1)-motherParticle->GetDaughterLabel(0)+1;
     if(fdectyp > 4) continue;  // exclude five or more particles decay
     if(motherIsPrimary){

     //check second leg
     if(i<fMCEvent->GetNumberOfTracks()-1) {
      particle2  = (AliVParticle *)fMCEvent->GetTrack(i+1);
      if(particle2->Particle()->GetMother(0)==particle->Particle()->GetMother(0)&&particle->Particle()->GetMother(1)==-1&&particle2->Particle()->GetMother(1)==-1){
       Skip2ndLeg=i+1;

       TLorentzVector dau1,dau2,ee,ee_orig;
       dau1.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->E());
       dau2.SetPxPyPzE(particle2->Px(),particle2->Py(),particle2->Pz(),particle2->E());  
       //create dielectron before resolution effects:
       ee=dau1+dau2;
       ee_orig=ee;

       // get info of the other particles in the decay: //////////////////////////////////////
       for(Int_t jj=motherParticle->GetDaughterLabel(0);jj<=motherParticle->GetDaughterLabel(1);jj++){
        if(jj==i||jj==Skip2ndLeg) {
        continue;
       }
       dau3Particle = (AliVParticle *)fMCEvent->GetTrack(jj);
       fdau3pdg= abs(dau3Particle->PdgCode());
      }

      ///////////////////////////////////////////////////////////////////////////////////////
      //FOR THE TIME BEING SKIP eta' -> omega e+ e- !!!!!!!!!!! ////////////////////////////
      //  skip as well phi -> pi0 e+ e- //
      ///////////////////////////////////////////////////////////////////////////////////////
      //if(fdectyp==3&&motherParticle->PdgCode()==331&&fdau3pdg==223) continue;
      //if(fdectyp==3&&motherParticle->PdgCode()==333&&fdau3pdg==111) continue; // skip as well phi -> pi0 e+ e-
      ///////////////////////////////////////////////////////////////////////////////////////

       //get index for histograms
       Int_t hindex[3];
       for(Int_t jj=0;jj<3;jj++){hindex[jj]=-1;};
       switch(motherParticle->PdgCode()){
        case 111:
         hindex[0]=0;
         break;
        case 221:
         hindex[0]=1;
         break;
        case 331:
         hindex[0]=2;
         if(fdectyp==3&&fdau3pdg==22) hindex[1]=3;
         if(fdectyp==3&&fdau3pdg==223) hindex[1]=4;
         break;
        case 113:
         hindex[0]=5;
         break;
        case 223:
         hindex[0]=6;
         if(fdectyp==2) hindex[1]=7;
         if(fdectyp==3&&fdau3pdg==111) hindex[1]=8;
         break;
        case 333:
         hindex[0]=9;
         if(fdectyp==2) hindex[1]=10;
         if(fdectyp==3&&fdau3pdg==221) hindex[1]=11;
         if(fdectyp==3&&fdau3pdg==111) hindex[1]=12;
         break;
        }

        hindex[2]=nInputParticles;

        if(hindex[0]<0) {
         printf("Error LMeeCocktail hindex[0]<0 \n");
         continue;
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
        if(particle->PdgCode()>0) { feeorigphiv=PhiV(dau1,dau2); }else{ feeorigphiv=PhiV(dau2,dau1);};

        //get the efficiency weight
        Int_t effbin=fhwEffpT->FindBin(dau1.Pt());
        fwEffpT=fhwEffpT->GetBinContent(effbin);
        effbin=fhwEffpT->FindBin(dau2.Pt());
        fwEffpT=fwEffpT*fhwEffpT->GetBinContent(effbin);

        //Resolution and acceptance  
        //-------------------------
        if(particle->PdgCode()>0) { dau1=ApplyResolution(dau1,-1,fResolType); }else{ dau1=ApplyResolution(dau1,1,fResolType);};
        if(particle2->PdgCode()>0) { dau2=ApplyResolution(dau2,-1,fResolType); }else{ dau2=ApplyResolution(dau2,1,fResolType);};
        fpass=kTRUE;
        if(dau1.Pt()<fMinPt||dau2.Pt()<fMinPt) fpass=kFALSE; //leg pT cut
        if(dau1.Pt()>fMaxPt||dau2.Pt()>fMaxPt) fpass=kFALSE; //leg pT cut
        if(TMath::Abs(dau1.Eta())>fMaxEta||TMath::Abs(dau2.Eta())>fMaxEta) fpass=kFALSE;

        //get the pair DCA (based in smeared pT)
        Float_t DCAtemplateLowEdge[] = {0., .3, .4, .6, 1., 2. };
        Float_t DCAtemplateUpEdge[]  =     {.3, .4, .6, 1., 2., 100000.}  ;
        for(int jj=0;jj<nbDCAtemplate;jj++){ //loop over DCA templates
         if(dau1.Pt()>=DCAtemplateLowEdge[jj]&&dau1.Pt()<DCAtemplateUpEdge[jj]){fd1DCA=fh_DCAtemplates[jj]->GetRandom();};
         if(dau2.Pt()>=DCAtemplateLowEdge[jj]&&dau2.Pt()<DCAtemplateUpEdge[jj]){fd2DCA=fh_DCAtemplates[jj]->GetRandom();};
        }
        fpairDCA=sqrt((pow(fd1DCA,2)+pow(fd2DCA,2))/2);

        //Fill tree words after resolution/acceptance
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
        if(particle->PdgCode()>0) { feephiv=PhiV(dau1,dau2); }else{ feephiv=PhiV(dau2,dau1);};
        fmotherpt=motherParticle->Pt();
        fmothermt=sqrt(pow(motherParticle->Particle()->GetCalcMass(),2)+pow(motherParticle->Pt(),2));
        fmotherp=motherParticle->P();
        fmotherm=motherParticle->Particle()->GetCalcMass();
        fmothereta=motherParticle->Eta();
        fmotherphi=motherParticle->Phi();
        //fID=1000*fdectyp+motherParticle->PdgCode();
        fID=motherParticle->PdgCode();
        fweight=particle->Particle()->GetWeight(); //get particle weight from generator

        //get multiplicity based weight:
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

        //Which ALT weight to use?:
        Double_t fwALT = fwEffpT; //by default use pt efficiency weight
        if(fALTweightType == 1) fwALT = fwMultmT;  //mT multiplicity weight
        if(fALTweightType == 11) fwALT = fwMultmT2;  //mT multiplicity weight, higher mult
        if(fALTweightType == 2) fwALT = fwMultpT;  //pT multiplicity weight
        if(fALTweightType == 22) fwALT = fwMultpT2;  //pT multiplicity weight, higher mult
 
        //Fill the tree
        if(fWriteTTree) teeTTree->Fill();

	
        //Fill the histograms
        if(fdectyp<4){ //skip for the moment 4-particle decays
         for(Int_t jj=0;jj<3;jj++){ // fill the different hindex -> particles
          if(hindex[jj]>-1){
           fmee_orig[hindex[jj]]->Fill(ee_orig.M(), fweight);
           if(fALTweightType == 1||fALTweightType == 11) {fmotherpT_orig[hindex[jj]]->Fill(fmothermt,fweight);
           }else if(fALTweightType == 2||fALTweightType == 22) {fmotherpT_orig[hindex[jj]]->Fill(fmotherpt,fweight);}
           fpteevsmee_orig[hindex[jj]]->Fill(ee_orig.M(),ee.Pt(), fweight);
           fphi_orig[hindex[jj]]->Fill(ee_orig.Phi(), fweight);
           frap_orig[hindex[jj]]->Fill(ee_orig.Rapidity(), fweight);
           fmee_orig_wALT[hindex[jj]]->Fill(ee_orig.M(), fweight*fwALT);
           if(fALTweightType == 1||fALTweightType == 11) {fmotherpT_orig_wALT[hindex[jj]]->Fill(fmothermt,fweight*fwALT);
           }else if(fALTweightType == 2||fALTweightType == 22) {fmotherpT_orig_wALT[hindex[jj]]->Fill(fmotherpt,fweight*fwALT);}
           fpteevsmee_orig_wALT[hindex[jj]]->Fill(ee_orig.M(),ee.Pt(), fweight*fwALT);
           if(fpass){
            fmee[hindex[jj]]->Fill(ee.M(), fweight);
            fpteevsmee[hindex[jj]]->Fill(ee.M(),ee.Pt(), fweight);
            fphi[hindex[jj]]->Fill(ee.Phi(), fweight);
            frap[hindex[jj]]->Fill(ee.Rapidity(), fweight);
            fDCAeevsmee->Fill(ee.M(),fpairDCA, fweight);
            fDCAeevsptee->Fill(ee.Pt(),fpairDCA, fweight);
            fmee_wALT[hindex[jj]]->Fill(ee.M(), fweight*fwALT);
            fpteevsmee_wALT[hindex[jj]]->Fill(ee.M(),ee.Pt(), fweight*fwALT);
           }
          }
         }
        }
  
        //Virtual photon generation
        //-------------------------
        //We will generate one virtual photon per histogrammed pion
        if(motherParticle->PdgCode()==111){
         //get mass and pt from histos and flat eta and phi
         Double_t VPHpT = ffVPHpT->GetRandom();
         Double_t VPHmass = fhKW->GetRandom();
         Double_t VPHeta = -1.+gRandom->Rndm()*2.;
         Double_t VPHphi = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
         TLorentzVector beam;
         beam.SetPtEtaPhiM(VPHpT,VPHeta,VPHphi,VPHmass);
         Double_t decaymasses[2] = {(TDatabasePDG::Instance()->GetParticle(11))->Mass(),
                                    (TDatabasePDG::Instance()->GetParticle(11))->Mass()};
         TGenPhaseSpace VPHgen;
         Bool_t SetDecay;
         SetDecay = VPHgen.SetDecay(beam, 2, decaymasses);
         if(SetDecay==0) printf(" ERROR: decay not permitted by kinematics \n");
         Double_t VPHweight = VPHgen.Generate();
         //get electrons from the decay
         TLorentzVector *decay1,*decay2;
         decay1 = VPHgen.GetDecay(0);
         decay2 = VPHgen.GetDecay(1);
         dau1.SetPxPyPzE(decay1->Px(),decay1->Py(),decay1->Pz(),decay1->E());
         dau2.SetPxPyPzE(decay2->Px(),decay2->Py(),decay2->Pz(),decay2->E());

         //create dielectron before resolution effects:
         ee=dau1+dau2;
         ee_orig=ee;

         //get index for histograms
         hindex[0]=nInputParticles-1;
         hindex[1]=-1;
         hindex[2]=-1;

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
        feeorigphiv=PhiV(dau1,dau2); 

        //get the efficiency weight
        Int_t effbin=fhwEffpT->FindBin(dau1.Pt());
        fwEffpT=fhwEffpT->GetBinContent(effbin);
        effbin=fhwEffpT->FindBin(dau2.Pt());
        fwEffpT=fwEffpT*fhwEffpT->GetBinContent(effbin);

        //Resolution and acceptance  
        //-------------------------
        dau1=ApplyResolution(dau1,1,fResolType);
        dau2=ApplyResolution(dau2,-1,fResolType);
        fpass=kTRUE;
        if(dau1.Pt()<fMinPt||dau2.Pt()<fMinPt) fpass=kFALSE; //leg pT cut
        if(dau1.Pt()>fMaxPt||dau2.Pt()>fMaxPt) fpass=kFALSE; //leg pT cut
        if(TMath::Abs(dau1.Eta())>fMaxEta||TMath::Abs(dau2.Eta())>fMaxEta) fpass=kFALSE;

        //get the pair DCA (based in smeared pT) -> no DCA for virtual photon for the moment
        //Float_t DCAtemplateLowEdge[nbDCAtemplate] = {0., .3, .4, .6, 1., 2. };
        //Float_t DCAtemplateUpEdge[nbDCAtemplate]  =     {.3, .4, .6, 1., 2., 100000.}  ;
        //for(int jj=0;jj<nbDCAtemplate;jj++){ //loop over DCA templates
        // if(dau1.Pt()>=DCAtemplateLowEdge[jj]&&dau1.Pt()<DCAtemplateUpEdge[jj]){fd1DCA=fh_DCAtemplates[jj]->GetRandom();};
        // if(dau2.Pt()>=DCAtemplateLowEdge[jj]&&dau2.Pt()<DCAtemplateUpEdge[jj]){fd2DCA=fh_DCAtemplates[jj]->GetRandom();};
        //}
        fpairDCA=10000.;

        //Fill tree words after resolution/acceptance
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
        feephiv=PhiV(dau1,dau2);
        fmotherpt=beam.Pt();
        fmothermt=sqrt(pow(beam.M(),2)+pow(beam.Pt(),2));
        fmotherp=beam.P();
        fmotherm=beam.M();
        fmothereta=beam.Eta();
        fmotherphi=beam.Phi();
        fID=0; // set ID to Zero for VPH
        fweight=VPHweight;
        //get multiplicity based weight:
        fwMultmT = 1; //no weight for photons so far

        //Fill the tree
        if(fWriteTTree) teeTTree->Fill();

        //Fill the histograms
         for(Int_t jj=0;jj<3;jj++){ // fill the different hindex -> particles
          if(hindex[jj]>-1){
           fmee_orig[hindex[jj]]->Fill(ee_orig.M(), VPHweight);
           fpteevsmee_orig[hindex[jj]]->Fill(ee_orig.M(),ee.Pt(), VPHweight);
           fphi_orig[hindex[jj]]->Fill(ee_orig.Phi(), VPHweight);
           frap_orig[hindex[jj]]->Fill(ee_orig.Rapidity(), VPHweight);
           if(fpass){
            fmee[hindex[jj]]->Fill(ee.M(), VPHweight);
            fpteevsmee[hindex[jj]]->Fill(ee.M(),ee.Pt(), VPHweight);
            fphi[hindex[jj]]->Fill(ee.Phi(), VPHweight);
            frap[hindex[jj]]->Fill(ee.Rapidity(), VPHweight);
           }
          }
         }
        }//----------- create one v.ph. per pi0
  

       } // legs coincide
      } // check the 2nd leg

     } // mother is primary
   } // pdgid==11 and HasMother
    
  }//MC particles loop
  
  //Clear buffers
  eBuff.clear();
  echBuff.clear();
  eweightBuff.clear();

}

//________________________________________________________________________
void AliAnalysisTaskLMeeCocktailMC::Terminate(const Option_t *)
{
  //fOutputContainer->Print(); // Will crash on GRID
}


//_______________________________________________________________________________________________
TLorentzVector AliAnalysisTaskLMeeCocktailMC::ApplyResolution(TLorentzVector vec, Char_t ch = 0, Int_t Run = 2 )
{
// from Theos LightFlavorGenerator, modified

  Double_t theta, phi, pt, p, px, py, pz, E, mass, eta;
  TLorentzVector resvec;

  mass  = 0.51099906e-3;
  pt    = vec.Pt();
  p     = vec.P();
  theta = vec.Theta();
  phi   = vec.Phi();
  eta   = vec.Eta();

  // Run == 0 --> no resolution is applied
  if(Run == 0){
    px   = p*sin(theta)*cos(phi);
    py   = p*sin(theta)*sin(phi);
    pz   = p*cos(theta);
    E    = sqrt(p*p + mass*mass);
    resvec.SetPxPyPzE(px,py,pz,E);
  }
  else if(Run == 1){
  
   TH1D *hisSlice(0x0);
   if(fArr){
     TH2D *hDeltaPtvsPt = static_cast<TH2D*> (fArr->At(0));
     //std::cout << "2D momentum reso hist = " << hDeltaPtvsPt->GetTitle() << std::endl;
     // Get the momentum slice histogram for sampling of the smearing
     // (in some input file versions this also contains a landau fit )
     // since histogram bins start at 1, we can use the bin number directly for the array index (first slice stored in position 1).
     Int_t histIndex = hDeltaPtvsPt->GetXaxis()->FindBin(pt);
     if (histIndex<1) histIndex=1; // in case some track is below the first p-bin (which currently starts at 100 MeV).
     if (histIndex>fArr->GetLast()) {
       //printf(" track with too high momentum (%f GeV), using last slice.\n", pt); // add some printouts when testing new input files!
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


  } else if(Run == 2){

   //smear pt
   Int_t ptbin=((TH2D*)(fArrResoPt->At(0)))->GetXaxis()->FindBin(pt);
   if(ptbin<1) ptbin=1;
   if(ptbin > fArrResoPt->GetLast()) ptbin = fArrResoPt->GetLast();
   Double_t smearing = ((TH1D*)(fArrResoPt->At(ptbin)))->GetRandom() * pt;
   Double_t sPt = pt - smearing;

   //smear eta
   ptbin=((TH2D*)(fArrResoEta->At(0)))->GetXaxis()->FindBin(pt);
   if(ptbin<1) ptbin=1;
   if(ptbin > fArrResoEta->GetLast()) ptbin = fArrResoEta->GetLast();
   smearing = ((TH1D*)(fArrResoEta->At(ptbin)))->GetRandom();
   Double_t sEta = eta - smearing;

   //smear phi
   ptbin=((TH2D*)(fArrResoPhi_Pos->At(0)))->GetXaxis()->FindBin(pt);
   if(ptbin<1) ptbin=1;
   if(ptbin > fArrResoPhi_Pos->GetLast()) ptbin = fArrResoPhi_Pos->GetLast();
   if(ch>0){
    smearing = ((TH1D*)(fArrResoPhi_Pos->At(ptbin)))->GetRandom();
   }else if(ch<0){
    smearing = ((TH1D*)(fArrResoPhi_Neg->At(ptbin)))->GetRandom();
   }
   Double_t sPhi = phi - smearing;

   //printf(" Original Pt = %f Phi %f Eta %f -> final pt = %f Phi %f Eta %f \n",pt,phi,eta,sPt,sPhi,sEta);

   Double_t sPx = sPt * cos(sPhi);
   Double_t sPy = sPt * sin(sPhi);
   Double_t sPz = sPt * sinh(sEta);
   Double_t sP = sPt * cosh(sEta);
   Double_t   sE    = sqrt(sP*sP+ mass*mass);

   resvec.SetPxPyPzE(sPx,sPy,sPz,sE);
  }

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


