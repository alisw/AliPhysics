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

//*-- Authors: Aleksei Pavlinov (WSU) 

//*

#include "AliEMCALJetMicroDst.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliEMCALJetFinder.h"
#include <TVector3.h>
#include <TROOT.h>
#include <TBrowser.h>
#include <TString.h>
#include <TParticle.h>

ClassImp(AliEMCALJetMicroDst)

TString nameTree("jetMDST"); // 7-feb-2002

TH1F*  hPtPart, *hNJet, *hPtJet;
TH2F*  hEtaPhiPart, *hEtaPhiJet;
TH1F*  hNcell, *hCellId, *hCellEt, *hSumEt;
TH1F*  hNgrid, *hGridId, *hGridEt, *hSumEtGrForJF;

extern "C" void sgpdge_(Int_t &i, Int_t &pdggea);

AliEMCALJetMicroDst::AliEMCALJetMicroDst(char *name, char *tit) : TNamed(name,tit)
{
  fFile  = 0;
  fTree  = 0;
  fDebug = 0;
    
//  Don't add histos to the current directory
//  TH1::AddDirectory(0);
  gROOT->cd();
  hPtPart     = new TH1F("hPtPart","P_{T} for partons", 300, 0., 300.);
  // 16-jan-2002 - new limit fo phi 
  hEtaPhiPart = new TH2F("hEtaPhiPart","#eta #phi distr.for partons after HSc", 
		       28, -0.7, 0.7, 21, 0.0, (2.0/3.0)*TMath::Pi());

  hNJet  = new TH1F("hNJet","number of jets", 11, -0.5, 10.5);
  hPtJet = new TH1F("hPtJet","P_{T} for jets", 500, 0., 500.);
  hEtaPhiJet = new TH2F("hEtaPhiJet","#eta #phi distr.for jets (W)", 
		       28, -0.7, 0.7, 21, 0.0, (2.0/3.0)*TMath::Pi());

  hNcell  = new TH1F("hNcell","#cell with de>0.0 for EMCAL", 1400, 0.0, 14000.);
  hCellId = new TH1F("hCellId","cell ID with de>0.0 for EMCAL", 1400, 0.0, 14000.);
  hCellEt = new TH1F("hCellEt","cell Et for EMCAL", 1000, 0.0, 10.);
  hSumEt  = new TH1F("hSumEt","sum Et for EMCAL", 1000, 0.0, 1000.);

  hNgrid  = new TH1F("hNgrid","#cell with de>0.0 in EMCAL grid for JF", 1400, 0.0, 14000.);
  hGridId = new TH1F("hGridId","cell ID with de>0.0 in EMCAL grid for JF", 1400, 0.0, 14000.);
  hGridEt = new TH1F("hGridEt","cell Et in EMCAL grid for JF", 1000, 0.0, 10.);
  hSumEtGrForJF  = new TH1F("hSumEtGrForJF","sum Et in EMCAL grid for JF", 1000, 0.0, 1000.);

  fListHist = MoveHistsToList("Hist For AliEMCALJetMicroDst", kFALSE); 
}

AliEMCALJetMicroDst::~AliEMCALJetMicroDst()
{
  if(fFile) fFile->Close();
}  

Bool_t AliEMCALJetMicroDst::Create(TFile *file)
{
  if(!file) {
    Error("Create", "define TFile for output\n");
    return kFALSE;
  }
  fFile     = file;
  fFileName = fFile->GetName();
  fFile->cd();
  fTree = new TTree(nameTree.Data(),"Temporary micro DST for jets analysis");
  // for jet calibration - 4-mar-2003
  fTree->Branch("decone", &decone, "decone/F");
  fTree->Branch("ptcone", &ptcone, "ptcone/F");
  // partons 
  fTree->Branch("npart", &npart, "npart/I");
  fTree->Branch("xpt",  xpt,  "xpt[npart]/F");
  fTree->Branch("xeta", xeta, "xeta[npart]/F");
  fTree->Branch("xphi", xphi, "xphi[npart]/F");
  // jets 
  fTree->Branch("njet", &njet, "njet/I");
  fTree->Branch("jet",  jet,  "jet[njet]/F");
  fTree->Branch("jetaw", jetaw, "jetaw[njet]/F");
  fTree->Branch("jphiw", jphiw, "jphiw[njet]/F");
  fTree->Branch("jetal", jetal, "jetal[njet]/F");
  fTree->Branch("jphil", jphil, "jphil[njet]/F");

  // Et in EMCAL itself 
  fTree->Branch("ncell", &ncell, "ncell/I");
  fTree->Branch("idcell", idcell, "idcell[ncell]/I");
  fTree->Branch("etcell", etcell, "etcell[ncell]/F");

  // Et in EMCAL grid for JF 
  fTree->Branch("ngrid", &ngrid,  "ngrid/I");
  fTree->Branch("idgrid", idgrid, "idgrid[ngrid]/I");
  fTree->Branch("etgrid", etgrid, "etgrid[ngrid]/F");

  // charge particle which hit to EMCAL
  fTree->Branch("nchp", &nchp, "nchp/I");
  fTree->Branch("pid", pid, "pid[nchp]/I");
  fTree->Branch("ppt", ppt, "ppt[nchp]/F");
  fTree->Branch("peta", peta, "peta[nchp]/F");
  fTree->Branch("pphi", pphi, "pphi[nchp]/F");

  return kTRUE;
}

Bool_t AliEMCALJetMicroDst::Create(const char *fname)
{
  TFile *file = new TFile(fname, "RECREATE");
  if(file) {
    //    fNameFile = fname;
    return Create(file);
  } else   return kFALSE;
}

Bool_t AliEMCALJetMicroDst::Open(const char *fname)
{
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

const Char_t* AliEMCALJetMicroDst::DefineName(const Int_t mode)
{
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
  if(file) fFile = file; 
  fFile->cd();
  fTree = (TTree*)fFile->Get(nameTree.Data());
  if(!fTree) return kFALSE;
  // for jet calibration - 4-mar-2003
  fTree->SetBranchAddress("decone",&decone);
  fTree->SetBranchAddress("ptcone",&ptcone);
  // partons
  fTree->SetBranchAddress("npart",&npart);
  fTree->SetBranchAddress("xpt",   xpt);
  fTree->SetBranchAddress("xeta",  xeta);
  fTree->SetBranchAddress("xphi",  xphi);
  // jets 
  fTree->SetBranchAddress("njet", &njet);
  fTree->SetBranchAddress("jet",   jet);
  fTree->SetBranchAddress("jetaw", jetaw);
  fTree->SetBranchAddress("jphiw", jphiw);
  fTree->SetBranchAddress("jetal", jetal);
  fTree->SetBranchAddress("jphil", jphil);
  // eT in EMCAL
  fTree->SetBranchAddress("ncell", &ncell);
  fTree->SetBranchAddress("idcell", idcell);
  fTree->SetBranchAddress("etcell", etcell);
  // eT in EMCAL grid for JF
  fTree->SetBranchAddress("ngrid", &ngrid);
  fTree->SetBranchAddress("idgrid", idgrid);
  fTree->SetBranchAddress("etgrid", etgrid);
  // 28-jan-2003
  fTree->SetBranchAddress("nchp", &nchp);
  fTree->SetBranchAddress("pid", pid);
  fTree->SetBranchAddress("ppt", ppt);
  fTree->SetBranchAddress("peta", peta);
  fTree->SetBranchAddress("pphi", pphi);

  return kTRUE;
}

void AliEMCALJetMicroDst::Print(Option_t* option) const
{
  if(option);
  if(fFile) {
    fFile->Print();
    if(fTree) fTree->Print();
    else Info("Print", "TRee is zero\n");
  } else {
     Info("Print", "File with TRee is closed \n Name of file %s(from fFileName", fFileName.Data());
  }

  Info("Print", "******* Current(last) event *****");
  printf("#partons %2i \n", npart);
  for(Int_t i=0; i<npart; i++){
    printf("     %1i Pt %7.1f eta  %7.4f phi  %7.4f \n",
    i, xpt[i],xeta[i],xphi[i]);  
  }
  printf("#jets    %2i \n", njet);
  for(Int_t i=0; i<njet; i++){
    printf("     %1i Et %7.1f etaw %7.4f phiw %7.4f \n",
    i,jet[i],jetaw[i],jphiw[i]);  
  }
  printf(" Title %s \n", GetTitle());
}

void AliEMCALJetMicroDst::Fill(AliRun *run, AliEMCALJetFinder* jetFinder, Int_t modeFilling)
{//  modeFilling >=1 - fill info for EMCAL grid
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
     ncell = 0; // 27-jan-2003
     ngrid = 0; // 27-jan-2003
     nchp  = 0;
     // negative - no signal 
     decone = -1.;
     ptcone = -1.;
  }

  FillJetsControl();  //24-jan-2003

  fTree->Fill();
}

void AliEMCALJetMicroDst::FillPartons(AliGenHijingEventHeader *header)
{
  TLorentzVector parton[4];
  header->GetJets(parton[0], parton[1], parton[2], parton[3]);

  npart = 4; // 
  for(Int_t i=0; i<4; i++){
     xpt[i]  = parton[i].Pt();
     xeta[i] = parton[i].Eta();
     xphi[i] = parton[i].Phi();
  }
}

void AliEMCALJetMicroDst::FillPartons()
{// for case of Pythia -> get info from full event record

  npart = 2;
  TParticle *MPart;
  Int_t ind;
  for(Int_t i=6; i<8; i++){
     MPart = gAlice->Particle(i);
     ind   = i-6;
     xpt[ind]  = MPart->Pt();
     xeta[ind] = MPart->Eta();
     xphi[ind] = MPart->Phi();
  }
}

void AliEMCALJetMicroDst::FillJets(AliEMCALJetFinder* jetFinder)
{
  njet = 0;
  if(fDebug>1) Info("FillJets", "Debug"); 
  if(!jetFinder) {
    if(fDebug>1) Info("FillJets", "jetFinder is zero"); 
    return;
  }
  njet = jetFinder->Njets();
  if(njet>10) {
    if(fDebug>1) Warning("FillJets", "wrong value of jetFinder->Njets() %i ", njet); 
    njet = 10;
  }
  //  hNJet->Fill(njet);
  if(fDebug>1) Info("FillJets", "njet %i", njet); 
  if(njet){
    for(Int_t i=0; i<njet; i++){
      jet[i]   = jetFinder->JetEnergy(i);
      jetaw[i] = jetFinder->JetEtaW(i);
      jphiw[i] = jetFinder->JetPhiW(i);
      jetal[i] = jetFinder->JetEtaL(i);
      jphil[i] = jetFinder->JetPhiL(i);
    }
  }
}

void AliEMCALJetMicroDst::FillEtForEMCAL(AliEMCALJetFinder* jetFinder)
{
   ncell = 0;
   TH2F *hid = jetFinder->GetLegoEMCAL();
   if(!hid) return;

   Double_t de = 0.;
   Int_t neta = hid->GetNbinsX(), nphi = hid->GetNbinsY();
   for(Int_t ieta=1; ieta<=neta; ieta++) {
      for(Int_t iphi=1; iphi<=nphi; iphi++) {
	 de = hid->GetBinContent(ieta,iphi);
         if(de > 0.0) {
	   etcell[ncell] = Float_t(de);
           idcell[ncell] = nphi*(ieta-1) + iphi;
	   ncell++;
           if(ncell >= 13824) break; 
	   // Info("FillEtForEMCAL", " ncell %i6 id %i6 de %f \n", ncell, idcell[ncell], etcell[ncell]); 
         }
      }
   }
   if(njet == 1) {
     // jet energy calculate around LP direction !!! - 10-mar-2003 
      decone = jetFinder->EMCALConeEnergy(jetal[0],jphil[0]); 
      ptcone = jetFinder->TrackConeEnergy(jetal[0],jphil[0]); // get from lego plot fo ch.part
      Info("FillEtForEMCAL", " njet %i Emcal in cone %f pt ch.part in cone %f\n", njet, decone, ptcone); 
      Info("FillEtForEMCAL", " jet - decone - ptcone : %9.2f\n", jet[0]-decone-ptcone);
   } else {
      decone = -1.;
      ptcone = -1.;
   }

   Info("FillEtForEMCAL", "neta %3i nphi %3i # array size %i Sum.Et %f\n", 
   neta,nphi, ncell, hid->Integral());
}

void AliEMCALJetMicroDst::FillEtForGrid(AliEMCALJetFinder* jetFinder)
{
   TH2F *hid = jetFinder->GetLego();
   if(!hid) return;

   FillArrays(hid, ngrid, idgrid, etgrid);
}

void AliEMCALJetMicroDst::FillArrays(TH2* hid, Int_t &n, Int_t *id, Float_t *et)
{
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
{// 28-jan-2003 for fullness ; 18-mar - sometimes 
  nchp = 0;
  Int_t gid=0;
  for(Int_t i=0; i<jetFinder->fNt; i++) {
    //     fPdgT[i];
     if(jetFinder->fTrackList[i] >= 1) {
        sgpdge_(jetFinder->fPdgT[i], gid);
        pid[nchp]  = gid; 
        ppt[nchp]  = jetFinder->fPtT[i];
        peta[nchp] = jetFinder->fEtaT[i];
        pphi[nchp] = jetFinder->fPhiT[i];
        nchp++;
     }
  }
  Info("FillChargedParticles", "fNtS %i : nchp %i -> %i\n", jetFinder->fNtS, nchp, jetFinder->fNtS - nchp);
}

void AliEMCALJetMicroDst::FillJetsControl()
{  // see FillJets(AliEMCALJetFinder* jetFinder) and FillPartons
  hNJet->Fill(njet);
  for(Int_t i=0; i<njet; i++){
     hPtJet->Fill(jet[i]);
     hEtaPhiJet->Fill(jetaw[i],jphiw[i]); 
  }

  for(Int_t i=0; i < npart; i++){
     hEtaPhiPart->Fill(xeta[i], xphi[i]);
     hPtPart->Fill(xpt[i]);
  }

  Double_t sum = 0.0;
  hNcell->Fill(ncell);
  for(Int_t i=0; i < ncell; i++){
    hCellId->Fill(idcell[i]);
    hCellEt->Fill(etcell[i]);
    sum += Double_t(etcell[i]);
  }
  hSumEt->Fill(sum);

  sum = 0.0;
  hNgrid->Fill(ngrid);
  for(Int_t i=0; i < ngrid; i++){
    hGridId->Fill(idgrid[i]);
    hGridEt->Fill(etgrid[i]);
    sum += Double_t(etgrid[i]);
  }
  hSumEtGrForJF->Fill(sum);
}

Int_t AliEMCALJetMicroDst::GetEntry(Int_t entry)
{ // Read contents of entry.
   if (!fTree) {
      Error("GetEntry", "define TTree");
      return -1;
   }
   return fTree->GetEntry(entry);
}

Bool_t AliEMCALJetMicroDst::GetParton(Int_t i, Float_t& pt, Float_t& eta, Float_t& phi)
{
  if(i>=0 && i<npart) {
    pt  = xpt[i];
    eta = xeta[i];
    phi = xphi[i];
    return kTRUE;
  } else return kFALSE; 
}

Bool_t AliEMCALJetMicroDst::GetParton(Int_t i, TVector3& vec)
{
  static Float_t pt, eta, phi;

  if(!GetParton(i, pt, eta, phi)) return kFALSE;

  FillVector(pt, eta, phi, vec);
  return kTRUE;
}

Bool_t AliEMCALJetMicroDst::GetJet(Int_t i,Int_t mode,Float_t& pt, Float_t& eta, Float_t& phi)
{// mode=1(W) mode=any(L)
  if(i>=0 && i<njet) {
    pt  = jet[i];
    if(mode==1) {
      eta = jetaw[i];
      phi = jphiw[i];
    } else {
      eta = jetal[i];
      phi = jphil[i];
    }
    return kTRUE;
  } else return kFALSE; 
}

Bool_t AliEMCALJetMicroDst::GetJet(Int_t i, Int_t mode, TVector3& vec)
{
  static Float_t pt, eta, phi;

  if(!GetJet(i, mode, pt, eta, phi)) return kFALSE;
  FillVector(pt, eta, phi, vec);
  return kTRUE;
}

void AliEMCALJetMicroDst::Test()
{
  if(!fFile || !fTree ) {
    Info("Test", "define file with proper TTree !");
    return;
  }
  Int_t nbytes=0, nb=0, nentries=Int_t(fTree->GetEntries());
  for(Int_t i=0; i<nentries; i++){
    nb = fTree->GetEntry(i);  
    nbytes += nb;
    for(Int_t j=0; j<npart; j++){
      hEtaPhiPart->Fill(xeta[j], xphi[j]);
      hPtPart->Fill(xpt[j]);
    }

    hNJet->Fill(njet);
    if(njet){
      for(Int_t j=0; j<njet; j++) {
        hPtJet->Fill(jet[j]);
        hEtaPhiJet->Fill(jetaw[j],jphiw[j]); 
      }
    }
  }
  Info("Test", "Entries %5i Bytes %10i\n", nentries, nbytes);
}

void AliEMCALJetMicroDst::FillVector(Float_t pt, Float_t eta, Float_t phi, TVector3& vec)
{ // service function 
  static Float_t px, py, pz;

  px = pt*TMath::Cos(phi);
  py = pt*TMath::Sin(phi);
  pz = pt*TMath::SinH(eta);  // sinh(eta) = cot(theta)

  vec.SetXYZ(px, py, pz);
}

void AliEMCALJetMicroDst::GetEtaPhi(Int_t id, Double_t &eta, Double_t &phi)
{ // see AliEMCALGeometry 
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

TVector3& AliEMCALJetMicroDst::GetCellVector(Int_t i)
{
  static Double_t eta,phi;
  static TVector3 vec;
  vec.SetXYZ(0.,0.,0.);
  if(i>=0 && i<ncell) {
     GetEtaPhi(idcell[i], eta, phi);
     vec.SetPtEtaPhi(Double_t(etcell[i]),eta,phi);
  }
  return vec;
}

TVector3& AliEMCALJetMicroDst::GetGridVector(Int_t i)
{
  static Double_t eta,phi;
  static TVector3 vec;
  vec.SetXYZ(0.,0.,0.);
  if(i>=0 && i<ngrid) {
     GetEtaPhi(idgrid[i], eta, phi);
     vec.SetPtEtaPhi(Double_t(etgrid[i]),eta,phi);
  }
  return vec;
}

Double_t AliEMCALJetMicroDst::GetSumInCone(TVector3 &jet,Int_t nc, Float_t *et,Float_t *eta,Float_t *phi, Double_t cellEtCut, Double_t rJet)
{
  static Double_t sum=0.;
  static TVector3 cell(0., 0., 0.);
  if(nc<=0 || et==0 || eta==0 || phi==0) {
    Error("GetSumInCone", "nc %d %f %f %f ", nc, et, eta, phi);
    return -1.;
  }

  sum=0.;
  //  jet.SetPtEtaPhi(jet[0],jetaw[0],jphiw[0]); // must be one jet !!
  Info("GetSumInCone", "jet.Mag() %f : njet %i\n", jet.Mag(), njet);
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
  Int_t nc = ncell;
  if(nc<=0) return 0.;

  Float_t *et=etcell, *eta=new Float_t[nc], *phi=new Float_t[nc];
  Double_t etaCell=0., phiCell=0., ET=0;

  for(Int_t i=0; i<nc; i++) {
    GetEtaPhi(idcell[i], etaCell, phiCell);
    eta[i] = etaCell;
    phi[i] = phiCell;
  }

  ET = GetSumInCone(jet, nc, et,eta,phi, cellEtCut,rJet);
  delete [] eta;
  delete [] phi;

  return ET;
}

Double_t AliEMCALJetMicroDst::GetTpcPtInCone(TVector3 &jet,Double_t cellEtCut, Double_t rJet)
{
  if(nchp<=0) return 0.;
  return GetSumInCone(jet, nchp, ppt,peta,pphi, cellEtCut,rJet);
}

Double_t AliEMCALJetMicroDst::GetSum(Int_t n, Float_t *ar, Double_t cut)
{ // 25-apr-2003
  Double_t sum=0.0;
  if(n<=0 || ar==0) return sum;
  for(Int_t i=0; i<n; i++) {if(ar[i]>=cut) sum += Double_t(ar[i]);}
  return sum;
}

void AliEMCALJetMicroDst::Close()
{
  fFile->Write(); 
  fTree->Print(); 
  fFile->Close();
  fFile = 0;
  fTree = 0;
}

void AliEMCALJetMicroDst::Browse(TBrowser* b)
{
   if(fTree)      b->Add((TObject*)fTree);
   if(fListHist)  b->Add((TObject*)fListHist);
   //   TObject::Browse(b);
}

Bool_t  AliEMCALJetMicroDst::IsPythiaDst()
{
  TString st(GetTitle());
  if(st.Contains("py",TString::kIgnoreCase)||!st.Contains("hijing", TString::kIgnoreCase)) return kTRUE;
  else return kFALSE;
}

Bool_t  AliEMCALJetMicroDst::IsFolder() const
{
  if(fTree || fListHist) return kTRUE;
  else                   return kFALSE;
}

TList* AliEMCALJetMicroDst::MoveHistsToList(char* name, Bool_t putToBrowser)
{
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
