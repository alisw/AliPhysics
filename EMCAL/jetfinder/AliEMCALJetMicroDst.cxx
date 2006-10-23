/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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

//----------------------------------------
//
// JetMicroDst to store information about
// jetfinding for offline analysis
//
//*-- Authors: Aleksei Pavlinov (WSU) 
//
//-----------------------------------------

//*

#include <TBrowser.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>

#include "AliEMCALJetFinder.h"
#include "AliEMCALJetMicroDst.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliHeader.h"
#include "AliMC.h"
#include "AliRun.h"

ClassImp(AliEMCALJetMicroDst)

TString gAliNameTree("jetMDST"); // 7-feb-2002

//TH1F*  fhPtPart, *fhNJet, *fhPtJet;
//TH2F*  fhEtaPhiPart, *fhEtaPhiJet;
//TH1F*  fhNcell, *fhCellId, *fhCellEt, *fhSumEt;
//TH1F*  fhNgrid, *fhGridId, *fhGridEt, *fhSumEtGrForJF;

extern "C" void sgpdge_(Int_t &i, Int_t &pdggea);

AliEMCALJetMicroDst::AliEMCALJetMicroDst(const char *name, const char *tit) : 
  TNamed(name,tit),
  fDebug(0),fFile(0),fTree(0),fListHist(0),fFileName(0),fdecone(0.),fptcone(0.),
  fnpart(0),fnjet(0),fncell(0),fngrid(0),fnchp(0),fhPtPart(0),fhNJet(0),fhPtJet(0),
  fhEtaPhiPart(0),fhEtaPhiJet(0),fhNcell(0),fhCellId(0),fhCellEt(0),fhSumEt(0),fhNgrid(0),
  fhGridId(0),fhGridEt(0),fhSumEtGrForJF(0)
{
	//constructor
  fFile  = 0;
  fTree  = 0;
  fDebug = 0;
    
//  Don't add histos to the current directory
//  TH1::AddDirectory(0);
  gROOT->cd();
  fhPtPart     = new TH1F("fhPtPart","P_{T} for partons", 300, 0., 300.);
  // 16-jan-2002 - new limit fo phi 
  fhEtaPhiPart = new TH2F("fhEtaPhiPart","#eta #phi distr.for partons after HSc", 
		       28, -0.7, 0.7, 21, 0.0, (2.0/3.0)*TMath::Pi());

  fhNJet  = new TH1F("fhNJet","number of jets", 11, -0.5, 10.5);
  fhPtJet = new TH1F("fhPtJet","P_{T} for jets", 500, 0., 500.);
  fhEtaPhiJet = new TH2F("fhEtaPhiJet","#eta #phi distr.for jets (W)", 
		       28, -0.7, 0.7, 21, 0.0, (2.0/3.0)*TMath::Pi());

  fhNcell  = new TH1F("fhNcell","#cell with de>0.0 for EMCAL", 1400, 0.0, 14000.);
  fhCellId = new TH1F("fhCellId","cell ID with de>0.0 for EMCAL", 1400, 0.0, 14000.);
  fhCellEt = new TH1F("fhCellEt","cell Et for EMCAL", 1000, 0.0, 10.);
  fhSumEt  = new TH1F("fhSumEt","sum Et for EMCAL", 1000, 0.0, 1000.);

  fhNgrid  = new TH1F("fhNgrid","#cell with de>0.0 in EMCAL grid for JF", 1400, 0.0, 14000.);
  fhGridId = new TH1F("fhGridId","cell ID with de>0.0 in EMCAL grid for JF", 1400, 0.0, 14000.);
  fhGridEt = new TH1F("fhGridEt","cell Et in EMCAL grid for JF", 1000, 0.0, 10.);
  fhSumEtGrForJF  = new TH1F("fhSumEtGrForJF","sum Et in EMCAL grid for JF", 1000, 0.0, 1000.);

  fListHist = MoveHistsToList("Hist For AliEMCALJetMicroDst", kFALSE); 
}


AliEMCALJetMicroDst::AliEMCALJetMicroDst(const AliEMCALJetMicroDst& jet) : 
  TNamed(jet.GetName(),jet.GetTitle()),
  fDebug(jet.fDebug),
  fFile(jet.fFile),
  fTree(jet.fTree),
  fListHist(jet.fListHist),
  fFileName(jet.fFileName),
  fdecone(jet.fdecone),
  fptcone(jet.fptcone),
  fnpart(jet.fnpart),
  fnjet(jet.fnjet),
  fncell(jet.fncell),
  fngrid(jet.fngrid),
  fnchp(jet.fnchp),
  fhPtPart(jet.fhPtPart),
  fhNJet(jet.fhNJet),
  fhPtJet(jet.fhPtJet),
  fhEtaPhiPart(jet.fhEtaPhiPart),
  fhEtaPhiJet(jet.fhEtaPhiJet),
  fhNcell(jet.fhNcell),
  fhCellId(jet.fhCellId),
  fhCellEt(jet.fhCellEt),
  fhSumEt(jet.fhSumEt),
  fhNgrid(jet.fhNgrid),
  fhGridId(jet.fhGridId),
  fhGridEt(jet.fhGridEt),
  fhSumEtGrForJF(jet.fhSumEtGrForJF)
{
  //copy constructor
}

AliEMCALJetMicroDst::~AliEMCALJetMicroDst()
{
	//destructor
  if(fFile) fFile->Close();
}  

Bool_t AliEMCALJetMicroDst::Create(TFile *file)
{
  // Creates the DST file
  if(!file) {
    Error("Create", "define TFile for output\n");
    return kFALSE;
  }
  fFile     = file;
  fFileName = fFile->GetName();
  fFile->cd();
  fTree = new TTree(gAliNameTree.Data(),"Temporary micro DST for jets analysis");
  // for jet calibration - 4-mar-2003
  fTree->Branch("fdecone", &fdecone, "fdecone/F");
  fTree->Branch("fptcone", &fptcone, "fptcone/F");
  // partons 
  fTree->Branch("fnpart", &fnpart, "fnpart/I");
  fTree->Branch("fxpt",  fxpt,  "fxpt[fnpart]/F");
  fTree->Branch("fxeta", fxeta, "fxeta[fnpart]/F");
  fTree->Branch("fxphi", fxphi, "fxphi[fnpart]/F");
  // jets 
  fTree->Branch("fnjet", &fnjet, "fnjet/I");
  fTree->Branch("fjet",  fjet,  "fjet[fnjet]/F");
  fTree->Branch("fjetaw", fjetaw, "fjetaw[fnjet]/F");
  fTree->Branch("fjphiw", fjphiw, "fjphiw[fnjet]/F");
  fTree->Branch("fjetal", fjetal, "fjetal[fnjet]/F");
  fTree->Branch("fjphil", fjphil, "fjphil[fnjet]/F");

  // Et in EMCAL itself 
  fTree->Branch("fncell", &fncell, "fncell/I");
  fTree->Branch("fidcell", fidcell, "fidcell[fncell]/I");
  fTree->Branch("fetcell", fetcell, "fetcell[fncell]/F");

  // Et in EMCAL grid for JF 
  fTree->Branch("fngrid", &fngrid,  "fngrid/I");
  fTree->Branch("fidgrid", fidgrid, "fidgrid[fngrid]/I");
  fTree->Branch("fetgrid", fetgrid, "fetgrid[fngrid]/F");

  // charge particle which hit to EMCAL
  fTree->Branch("fnchp", &fnchp, "fnchp/I");
  fTree->Branch("fpid", fpid, "fpid[fnchp]/I");
  fTree->Branch("fppt", fppt, "fppt[fnchp]/F");
  fTree->Branch("fpeta", fpeta, "fpeta[fnchp]/F");
  fTree->Branch("fpphi", fpphi, "fpphi[fnchp]/F");

  return kTRUE;
}

Bool_t AliEMCALJetMicroDst::Create(const char *fname)
{
	// Create member
  TFile *file = new TFile(fname, "RECREATE");
  if(file) {
    //    fNameFile = fname;
    return Create(file);
  } else   return kFALSE;
}

Bool_t AliEMCALJetMicroDst::Open(const char *fname)
{
	//Open member
  if(fFile && fFile->IsOpen()) fFile->Close();
  if(strlen(fname)) fName = fname;
  TFile *file = new TFile(fName.Data(), "READ");
  if(file) {
    Bool_t ini =  Initialize(file);
    Info("Open", "open file %s : initialize TTree %i",fName.Data(), Int_t(ini)); 
    return ini;
  } else {
    Error("Open", "can not open file %s",fName.Data()); 
    return kFALSE;
  }
}

const Char_t* AliEMCALJetMicroDst::DefineName(Int_t mode)
{
	//DefineName member
  static TString dir, name;
  //  dir = "jetDST/"; // 24-jan-2003
  dir = "/auto/alice/pavlinov/jet/microDST/"; // 24-jan-2003
  switch (mode) {
  case 1: // for characteristic of BG
    name = dir + "mDst1_1.root"; // Bg 2000 - first version 
    SetTitle("Bg2000");
    break;
  case 2: // for characteristic of BG
    name = dir + "mDst2_1.root"; // Bg 4000 - first version 
    SetTitle("Bg4000");
    break;
  case 3: // for characteristic of BG
    name = dir + "mDst3_1.root"; // Bg 8000 - first version 
    SetTitle("Bg8000");
    break;
  case 4: // Central Hijing - 18-mar-2003
    name  = dir  + "march/";
    name += "jF_R0.50MinCell1.0PtCut0.0EtSeed4.0MinEt40.0BGSubtr0SF11.6Smear0Eff0HijingCentral.root";
    SetTitle("HijingCentral");
    break;
  case 5: // Para Hijing (Dn/Dy=8000) - 21-mar-2003
    name  = dir  + "march/";
    name += "jF_R0.50MinCell1.0PtCut0.0EtSeed4.0MinEt40.0BGSubtr0SF11.6Smear0Eff0ParaHijing8000.root";
    SetTitle("HIJINGparaDnDy8000");
    break;
  case 6: // Para Hijing (Dn/Dy=4000) - 21-mar-2003
    name  = dir  + "march/";
    name += "jF_R0.50MinCell1.0PtCut0.0EtSeed4.0MinEt40.0BGSubtr0SF11.6Smear0Eff0ParaHijing4000.root";
    SetTitle("HIJINGparaDnDy4000");
    break;
  case 7: // Para Hijing (Dn/Dy=2000) - 21-mar-2003
    name  = dir  + "march/";
    name += "jF_R0.50MinCell1.0PtCut0.0EtSeed4.0MinEt40.0BGSubtr0SF11.6Smear0Eff0ParaHijing2000.root";
    SetTitle("HIJINGparaDnDy2000");
    break;
  case 11: // pure PYTHIA with default value of parameters
    name = dir + "jF_R0.50MinCell0.0PtCut0.0EtSeed8.0MinEt40.0BGSubtr0SF11.6.root";
    break;
  case 12: // 0 + background 
    name = dir + "jF_R0.50MinCell0.0PtCut0.0EtSeed8.0MinEt40.0BGSubtr0SF11.6kBackground2000.root";
    break;
    // calibration case 
  case 101:
    name  = dir  + "march/";
    name += "Pythia100_1.root";
    SetTitle("Pythia100_1");
    break;
  case 102: // 2-apr-2003
    name  = dir  + "march/";
    name += "Pythia50_1.root";
    SetTitle("Pythia50_1");
    //    name = "microDst3th.root"; // 101 + (smearing and eff) - 14-mar-2003
    break;
  case 103:// 4-apr-2003
    name  = dir  + "march/";
    name += "Pythia200_1.root";
    SetTitle("Pythia200_1");
    //    name = "microDst4th.root"; // 102 + MinCell= 1.0
    break;
  default:
    Fatal("DefineName", "NO  D E F A U L T : mode %i\n", mode);
  }
  Info("DefineName", "mode %5i file : %s : Title %s\n", mode, name.Data(), GetTitle());
  return name.Data();
}

Bool_t AliEMCALJetMicroDst::Initialize(TFile *file)
{
	// Initialize method
  if(file) fFile = file; 
  fFile->cd();
  fTree = (TTree*)fFile->Get(gAliNameTree.Data());
  if(!fTree) return kFALSE;
  // for jet calibration - 4-mar-2003
  fTree->SetBranchAddress("fdecone",&fdecone);
  fTree->SetBranchAddress("fptcone",&fptcone);
  // partons
  fTree->SetBranchAddress("fnpart",&fnpart);
  fTree->SetBranchAddress("fxpt",   fxpt);
  fTree->SetBranchAddress("fxeta",  fxeta);
  fTree->SetBranchAddress("fxphi",  fxphi);
  // jets 
  fTree->SetBranchAddress("fnjet", &fnjet);
  fTree->SetBranchAddress("fjet",   fjet);
  fTree->SetBranchAddress("fjetaw", fjetaw);
  fTree->SetBranchAddress("fjphiw", fjphiw);
  fTree->SetBranchAddress("fjetal", fjetal);
  fTree->SetBranchAddress("fjphil", fjphil);
  // eT in EMCAL
  fTree->SetBranchAddress("fncell", &fncell);
  fTree->SetBranchAddress("fidcell", fidcell);
  fTree->SetBranchAddress("fetcell", fetcell);
  // eT in EMCAL grid for JF
  fTree->SetBranchAddress("fngrid", &fngrid);
  fTree->SetBranchAddress("fidgrid", fidgrid);
  fTree->SetBranchAddress("fetgrid", fetgrid);
  // 28-jan-2003
  fTree->SetBranchAddress("fnchp", &fnchp);
  fTree->SetBranchAddress("fpid", fpid);
  fTree->SetBranchAddress("fppt", fppt);
  fTree->SetBranchAddress("fpeta", fpeta);
  fTree->SetBranchAddress("fpphi", fpphi);

  return kTRUE;
}

void AliEMCALJetMicroDst::Print(Option_t* option) const
{
	// 
  if(option);
  if(fFile) {
    fFile->Print();
    if(fTree) fTree->Print();
    else Info("Print", "TRee is zero\n");
  } else {
     Info("Print", "File with TRee is closed \n Name of file %s(from fFileName", fFileName.Data());
  }

  Info("Print", "******* Current(last) event *****");
  printf("#partons %2i \n", fnpart);
  for(Int_t i=0; i<fnpart; i++){
    printf("     %1i Pt %7.1f eta  %7.4f phi  %7.4f \n",
    i, fxpt[i],fxeta[i],fxphi[i]);  
  }
  printf("#jets    %2i \n", fnjet);
  for(Int_t i=0; i<fnjet; i++){
    printf("     %1i Et %7.1f etaw %7.4f phiw %7.4f \n",
    i,fjet[i],fjetaw[i],fjphiw[i]);  
  }
  printf(" Title %s \n", GetTitle());
}

void AliEMCALJetMicroDst::Fill(AliRun *run, AliEMCALJetFinder* jetFinder, Int_t modeFilling)
{
  //  modeFilling >=1 - fill info for EMCAL grid
  if(!run) run = gAlice;
  AliGenEventHeader* evHeader = run->GetHeader()->GenEventHeader();
  TString tmp(evHeader->ClassName());
  if(tmp.Contains("Hijing")) {
     AliGenHijingEventHeader *hijEvHeader = (AliGenHijingEventHeader*)evHeader;
     FillPartons(hijEvHeader);
  } else if(tmp.Contains("Pythia")) {
     FillPartons();     
  } else {
    Error("Fill", "Wrong type of generator -> %s \n Info about partons will be absent",tmp.Data()); 
  }

  FillJets(jetFinder);

  if(modeFilling >= 1) {
     FillEtForEMCAL(jetFinder);
     FillEtForGrid(jetFinder);
     FillChargeParticles(jetFinder);
  } else {
     fncell = 0; // 27-jan-2003
     fngrid = 0; // 27-jan-2003
     fnchp  = 0;
     // negative - no signal 
     fdecone = -1.;
     fptcone = -1.;
  }

  FillJetsControl();  //24-jan-2003

  fTree->Fill();
}

void AliEMCALJetMicroDst::FillPartons(AliGenHijingEventHeader *header)
{
  //Make partons arrays
  TLorentzVector parton[4];
  header->GetJets(parton[0], parton[1], parton[2], parton[3]);

  fnpart = 4; // 
  for(Int_t i=0; i<4; i++){
     fxpt[i]  = parton[i].Pt();
     fxeta[i] = parton[i].Eta();
     fxphi[i] = parton[i].Phi();
  }
}

void AliEMCALJetMicroDst::FillPartons()
{
  // for case of Pythia -> get info from full event record

  fnpart = 2;
  TParticle *mPart;
  Int_t ind;
  for(Int_t i=6; i<8; i++){
     mPart = gAlice->GetMCApp()->Particle(i);
     ind   = i-6;
     fxpt[ind]  = mPart->Pt();
     fxeta[ind] = mPart->Eta();
     fxphi[ind] = mPart->Phi();
  }
}

void AliEMCALJetMicroDst::FillJets(AliEMCALJetFinder* jetFinder)
{
  // Fill Jets
  fnjet = 0;
  if(fDebug>1) Info("FillJets", "Debug"); 
  if(!jetFinder) {
    if(fDebug>1) Info("FillJets", "jetFinder is zero"); 
    return;
  }
  fnjet = jetFinder->Njets();
  if(fnjet>10) {
    if(fDebug>1) Warning("FillJets", "wrong value of jetFinder->Njets() %i ", fnjet); 
    fnjet = 10;
  }
  //  fhNJet->Fill(njet);
  if(fDebug>1) Info("FillJets", "njet %i", fnjet); 
  if(fnjet){
    for(Int_t i=0; i<fnjet; i++){
      fjet[i]   = jetFinder->JetEnergy(i);
      fjetaw[i] = jetFinder->JetEtaW(i);
      fjphiw[i] = jetFinder->JetPhiW(i);
      fjetal[i] = jetFinder->JetEtaL(i);
      fjphil[i] = jetFinder->JetPhiL(i);
    }
  }
}

void AliEMCALJetMicroDst::FillEtForEMCAL(AliEMCALJetFinder* jetFinder)
{
  // Fill Et for EMCAL
   fncell = 0;
   TH2F *hid = jetFinder->GetLegoEMCAL();
   if(!hid) return;

   Double_t de = 0.;
   Int_t neta = hid->GetNbinsX(), nphi = hid->GetNbinsY();
   for(Int_t ieta=1; ieta<=neta; ieta++) {
      for(Int_t iphi=1; iphi<=nphi; iphi++) {
	 de = hid->GetBinContent(ieta,iphi);
         if(de > 0.0) {
	   fetcell[fncell] = Float_t(de);
           fidcell[fncell] = nphi*(ieta-1) + iphi;
	   fncell++;
           if(fncell >= 13824) break; 
	   // Info("FillEtForEMCAL", " ncell %i6 id %i6 de %f \n", ncell, idcell[ncell], etcell[ncell]); 
         }
      }
   }
   if(fnjet == 1) {
     // jet energy calculate around LP direction !!! - 10-mar-2003 
      fdecone = jetFinder->EMCALConeEnergy(fjetal[0],fjphil[0]); 
      fptcone = jetFinder->TrackConeEnergy(fjetal[0],fjphil[0]); // get from lego plot fo ch.part
      Info("FillEtForEMCAL", " njet %i Emcal in cone %f pt ch.part in cone %f\n", fnjet, fdecone, fptcone); 
      Info("FillEtForEMCAL", " jet - decone - ptcone : %9.2f\n", fjet[0]-fdecone-fptcone);
   } else {
      fdecone = -1.;
      fptcone = -1.;
   }

   Info("FillEtForEMCAL", "neta %3i nphi %3i # array size %i Sum.Et %f\n", 
   neta,nphi, fncell, hid->Integral());
}

void AliEMCALJetMicroDst::FillEtForGrid(AliEMCALJetFinder* jetFinder)
{
  // Fill ET for Grid
   TH2F *hid = jetFinder->GetLego();
   if(!hid) return;

   FillArrays(hid, fngrid, fidgrid, fetgrid);
}

void AliEMCALJetMicroDst::FillArrays(TH2* hid, Int_t &n, Int_t *id, Float_t *et)
{
  // Fill arays
   n = 0;
   Double_t de = 0.;
   Int_t neta = hid->GetNbinsX(), nphi = hid->GetNbinsY();
   for(Int_t ieta=1; ieta<=neta; ieta++) {
      for(Int_t iphi=1; iphi<=nphi; iphi++) {
	 de = hid->GetBinContent(ieta,iphi);
         if(de > 0.0) {
	   et[n] = Float_t(de);
           id[n] = nphi*(ieta-1) + iphi;
	   n++;
           if(n >= 13824) break; 
         }
      }
   }
   Info("FillArrays", "neta %3i nphi %3i # array size %i Sum.Et %f\n", 
   neta, nphi, n, hid->Integral());
}

void AliEMCALJetMicroDst::FillChargeParticles(AliEMCALJetFinder* jetFinder)
{
  // 28-jan-2003 for fullness ; 18-mar - sometimes 
  fnchp = 0;
  Int_t gid=0;
  for(Int_t i=0; i<jetFinder->fNt; i++) {
    //     fPdgT[i];
     if(jetFinder->fTrackList[i] >= 1) {
        sgpdge_(jetFinder->fPdgT[i], gid);
        fpid[fnchp]  = gid; 
        fppt[fnchp]  = jetFinder->fPtT[i];
        fpeta[fnchp] = jetFinder->fEtaT[i];
        fpphi[fnchp] = jetFinder->fPhiT[i];
        fnchp++;
     }
  }
  Info("FillChargedParticles", "fNtS %i : nchp %i -> %i\n", jetFinder->fNtS, fnchp, jetFinder->fNtS - fnchp);
}

void AliEMCALJetMicroDst::FillJetsControl()
{  
  // see FillJets(AliEMCALJetFinder* jetFinder) and FillPartons
  fhNJet->Fill(fnjet);
  for(Int_t i=0; i<fnjet; i++){
     fhPtJet->Fill(fjet[i]);
     fhEtaPhiJet->Fill(fjetaw[i],fjphiw[i]); 
  }

  for(Int_t i=0; i < fnpart; i++){
     fhEtaPhiPart->Fill(fxeta[i], fxphi[i]);
     fhPtPart->Fill(fxpt[i]);
  }

  Double_t sum = 0.0;
  fhNcell->Fill(fncell);
  for(Int_t i=0; i < fncell; i++){
    fhCellId->Fill(fidcell[i]);
    fhCellEt->Fill(fetcell[i]);
    sum += Double_t(fetcell[i]);
  }
  fhSumEt->Fill(sum);

  sum = 0.0;
  fhNgrid->Fill(fngrid);
  for(Int_t i=0; i < fngrid; i++){
    fhGridId->Fill(fidgrid[i]);
    fhGridEt->Fill(fetgrid[i]);
    sum += Double_t(fetgrid[i]);
  }
  fhSumEtGrForJF->Fill(sum);
}

Int_t AliEMCALJetMicroDst::GetEntry(Int_t entry)
{ 
  // Read contents of entry.
   if (!fTree) {
      Error("GetEntry", "define TTree");
      return -1;
   }
   return fTree->GetEntry(entry);
}

Bool_t AliEMCALJetMicroDst::GetParton(Int_t i, Float_t& pt, Float_t& eta, Float_t& phi) const
{
  // Get parton
  if(i>=0 && i<fnpart) {
    pt  = fxpt[i];
    eta = fxeta[i];
    phi = fxphi[i];
    return kTRUE;
  } else return kFALSE; 
}

Bool_t AliEMCALJetMicroDst::GetParton(Int_t i, TVector3& vec) const
{
  // Get Parton
  static Float_t pt, eta, phi;

  if(!GetParton(i, pt, eta, phi)) return kFALSE;

  FillVector(pt, eta, phi, vec);
  return kTRUE;
}

Bool_t AliEMCALJetMicroDst::GetJet(Int_t i,Int_t mode,Float_t& pt, Float_t& eta, Float_t& phi) const
{
  // mode=1(W) mode=any(L)
  if(i>=0 && i<fnjet) {
    pt  = fjet[i];
    if(mode==1) {
      eta = fjetaw[i];
      phi = fjphiw[i];
    } else {
      eta = fjetal[i];
      phi = fjphil[i];
    }
    return kTRUE;
  } else return kFALSE; 
}

Bool_t AliEMCALJetMicroDst::GetJet(Int_t i, Int_t mode, TVector3& vec) const 
{
  // Get Jet
  static Float_t pt, eta, phi;

  if(!GetJet(i, mode, pt, eta, phi)) return kFALSE;
  FillVector(pt, eta, phi, vec);
  return kTRUE;
}

void AliEMCALJetMicroDst::Test()
{
  // Test
  if(!fFile || !fTree ) {
    Info("Test", "define file with proper TTree !");
    return;
  }
  Int_t nbytes=0, nb=0, nentries=Int_t(fTree->GetEntries());
  for(Int_t i=0; i<nentries; i++){
    nb = fTree->GetEntry(i);  
    nbytes += nb;
    for(Int_t j=0; j<fnpart; j++){
      fhEtaPhiPart->Fill(fxeta[j], fxphi[j]);
      fhPtPart->Fill(fxpt[j]);
    }

    fhNJet->Fill(fnjet);
    if(fnjet){
      for(Int_t j=0; j<fnjet; j++) {
        fhPtJet->Fill(fjet[j]);
        fhEtaPhiJet->Fill(fjetaw[j],fjphiw[j]); 
      }
    }
  }
  Info("Test", "Entries %5i Bytes %10i\n", nentries, nbytes);
}

void AliEMCALJetMicroDst::FillVector(Float_t pt, Float_t eta, Float_t phi, TVector3& vec)
{ 
  // service function 
  static Float_t px, py, pz;

  px = pt*TMath::Cos(phi);
  py = pt*TMath::Sin(phi);
  pz = pt*TMath::SinH(eta);  // sinh(eta) = cot(theta)

  vec.SetXYZ(px, py, pz);
}

void AliEMCALJetMicroDst::GetEtaPhi(Int_t id, Double_t &eta, Double_t &phi) const
{ 
  // see AliEMCALGeometry 
  static Int_t ieta, iphi, nphi=144, neta=96;
  static Double_t phiMax=(2.0/3.0)*TMath::Pi(), phiMin=0.0;
  static Double_t phiStep=(phiMax-phiMin)/nphi, phiBeg = phiMin + phiStep/2.; 
  static Double_t etaMax=0.7, etaMin=-etaMax;
  static Double_t etaStep=(etaMax-etaMin)/neta, etaBeg = etaMin + etaStep/2.;

  ieta = (id-1)/nphi + 1; // id = nphi*(ieta-1) + iphi
  iphi = id - nphi*(ieta-1);
  if(ieta<=0 || ieta>neta) {
    Fatal("GetEtaPhi", "wrong id %i(ieta %i,iphi %i) : nphi %i neta %i", id,iphi,ieta, nphi,neta);

  }

  eta  = etaBeg + etaStep*(ieta-1);
  phi  = phiBeg + phiStep*(iphi-1);
}

TVector3& AliEMCALJetMicroDst::GetCellVector(Int_t i) const 
{
  // Get cell vector
  static Double_t eta,phi;
  static TVector3 vec;
  vec.SetXYZ(0.,0.,0.);
  if(i>=0 && i<fncell) {
     GetEtaPhi(fidcell[i], eta, phi);
     vec.SetPtEtaPhi(Double_t(fetcell[i]),eta,phi);
  }
  return vec;
}

TVector3& AliEMCALJetMicroDst::GetGridVector(Int_t i) const 
{
  // Get grid vector
  static Double_t eta,phi;
  static TVector3 vec;
  vec.SetXYZ(0.,0.,0.);
  if(i>=0 && i<fngrid) {
     GetEtaPhi(fidgrid[i], eta, phi);
     vec.SetPtEtaPhi(Double_t(fetgrid[i]),eta,phi);
  }
  return vec;
}

Double_t AliEMCALJetMicroDst::GetSumInCone(TVector3 &jet,Int_t nc, Float_t *et,Float_t *eta,Float_t *phi, Double_t cellEtCut, Double_t rJet) const
{
  // Get Sum in cone
  static Double_t sum=0.;
  static TVector3 cell(0., 0., 0.);
  if(nc<=0 || et==0 || eta==0 || phi==0) {
    Error("GetSumInCone", "nc %d %f %f %f ", nc, et, eta, phi);
    return -1.;
  }

  sum=0.;
  //  jet.SetPtEtaPhi(jet[0],jetaw[0],jphiw[0]); // must be one jet !!
  Info("GetSumInCone", "jet.Mag() %f : njet %i\n", jet.Mag(), fnjet);
  for(Int_t i=0; i<nc; i++){
    if(et[i] < cellEtCut)    continue;
    cell.SetPtEtaPhi(et[i], eta[i], phi[i]);
    if(jet.DeltaR(cell) > rJet) continue;
    sum += et[i];
  }
  Info("GetSumCone", "Sum %f \n", sum);
  return sum;
}

Double_t AliEMCALJetMicroDst::GetEmcalEtInCone(TVector3 &jet, Double_t cellEtCut, Double_t rJet)
{
  // Get EMCAL Et in cone
  Int_t nc = fncell;
  if(nc<=0) return 0.;

  Float_t *et=fetcell, *eta=new Float_t[nc], *phi=new Float_t[nc];
  Double_t etaCell=0., phiCell=0., eTotal=0;

  for(Int_t i=0; i<nc; i++) {
    GetEtaPhi(fidcell[i], etaCell, phiCell);
    eta[i] = etaCell;
    phi[i] = phiCell;
  }

  eTotal = GetSumInCone(jet, nc, et,eta,phi, cellEtCut,rJet);
  delete [] eta;
  delete [] phi;

  return eTotal;
}

Double_t AliEMCALJetMicroDst::GetTpcPtInCone(TVector3 &jet,Double_t cellEtCut, Double_t rJet)
{
  // Get TPC PT in cone
  if(fnchp<=0) return 0.;
  return GetSumInCone(jet, fnchp, fppt,fpeta,fpphi, cellEtCut,rJet);
}

Double_t AliEMCALJetMicroDst::GetSum(Int_t n, Float_t *ar, Double_t cut) const 
{ 
  // 25-apr-2003
  Double_t sum=0.0;
  if(n<=0 || ar==0) return sum;
  for(Int_t i=0; i<n; i++) {if(ar[i]>=cut) sum += Double_t(ar[i]);}
  return sum;
}

void AliEMCALJetMicroDst::Close()
{
  // Close
  fFile->Write(); 
  fTree->Print(); 
  fFile->Close();
  fFile = 0;
  fTree = 0;
}

void AliEMCALJetMicroDst::Browse(TBrowser* b) const
{
  // Browse
   if(fTree)      b->Add((TObject*)fTree);
   if(fListHist)  b->Add((TObject*)fListHist);
   //   TObject::Browse(b);
}

Bool_t  AliEMCALJetMicroDst::IsPythiaDst() const 
{
  // Is Pythia DST
  TString st(GetTitle());
  if(st.Contains("py",TString::kIgnoreCase)||!st.Contains("hijing", TString::kIgnoreCase)) return kTRUE;
  else return kFALSE;
}

Bool_t  AliEMCALJetMicroDst::IsFolder() const
{
  // Is folder
  if(fTree || fListHist) return kTRUE;
  else                   return kFALSE;
}

TList* AliEMCALJetMicroDst::MoveHistsToList(const char* name, Bool_t putToBrowser)
{
  // Move HIST to list
  gROOT->cd();
  TIter nextHist(gDirectory->GetList());
  TList *list = new TList;
  list->SetName(name);
  TObject *objHist;
  while((objHist=nextHist())){
    if (!objHist->InheritsFrom("TH1")) continue;
    ((TH1*)objHist)->SetDirectory(0); // Remove from gROOT
    list->Add(objHist);
  }
  if(putToBrowser) gROOT->GetListOfBrowsables()->Add((TObject*)list);
  return list;
}

void AliEMCALJetMicroDst::FillH1(TList *l, Int_t ind, Double_t x, Double_t w)
{
  //Fill 1d histogram with input data

  static TH1* hid=0;
  if(l == 0) return;
  if(ind < l->GetSize()){
    hid = (TH1*)l->At(ind);
    hid->Fill(x,w);
  }
}

void AliEMCALJetMicroDst::FillH2(TList *l, Int_t ind, Double_t x, Double_t y, Double_t w)
{
  //Fill 2d histogram with input data

  static TH2* hid=0;
  if(l == 0) return;
  if(ind < l->GetSize()){
    hid = (TH2*)l->At(ind);
    hid->Fill(x,y,w);
  }
}

int AliEMCALJetMicroDst::SaveListOfHists(TList *list,const char* name,Bool_t kSingleKey,const char* opt)
{
  //Save histograms to file

  printf(" Name of out file |%s|\n", name); 
  int save = 0;
  if(list && list->GetSize() && strlen(name)){
    TString nf(name); 
    if(nf.Contains(".root") == kFALSE) nf += ".root";
    TFile file(nf.Data(),opt);
    TIter nextHist(list);
    TObject* objHist=0;
    int nh=0;
    if(kSingleKey) {
       file.cd();
       list->Write(list->GetName(),TObject::kSingleKey);
       list->ls();
       save = 1;
    } else {
      while((objHist=nextHist())) { // loop over list 
        if(objHist->InheritsFrom("TH1")) {
          TH1* hid = (TH1*)objHist;
          file.cd();
          hid->Write();
          nh++;
          printf("Save hist. %s \n",hid ->GetName());
        }
      }
      printf("%i hists. save to file -> %s\n", nh,file.GetName());
      if(nh>0) save = 1;
    }
    file.Close();
  } else {
    printf("TAliasPAI::saveListOfHists : N O  S A V I N G \n");
    if(list==0) printf("List of object 0 : %p \n", (void*)list);
    else printf("Size of list %i \n", list->GetSize());
  }
  return save;
}

void AliEMCALJetMicroDst::Sgpdge(Int_t pdgId, Int_t &gId)
{ // 8-nov-05 
  sgpdge_(pdgId, gId);
}
