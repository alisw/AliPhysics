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

AliEMCALJetMicroDst::AliEMCALJetMicroDst(char *name, char *tit) : TNamed(name,tit)
{
  fFile  = 0;
  fTree  = 0;
  fDebug = 0;
  fName  = "dstTest.root"; // default name
    
//  Don't add histos to the current directory
//  TH1::AddDirectory(0);
  gROOT->cd();
  hPtPart     = new TH1F("hPtPart","P_{T} for partons", 300, 0., 300.);
  hEtaPhiPart = new TH2F("hEtaPhiPart","#eta #phi distr.for partons after HSc", 
		       28, -0.7, 0.7, 21, 0.0, TMath::Pi()*2./3.);

  hNJet  = new TH1F("hNJet","number of jets", 11, -0.5, 10.5);
  hPtJet = new TH1F("hPtJet","P_{T} for jets", 300, 0., 300.);
  hEtaPhiJet = new TH2F("hEtaPhiJet","#eta #phi distr.for jets (W)", 
		       28, -0.7, 0.7, 21, 0.0, TMath::Pi()*2./3.);
  fListHist = MoveHistsToList("Hist For AliEMCALJetMicroDst", kFALSE); 
}

AliEMCALJetMicroDst::~AliEMCALJetMicroDst()
{
  if(fFile) fFile->Close();
}  

Bool_t AliEMCALJetMicroDst::Create(TFile *file)
{
  if(!file) {
    printf("<E> AliEMCALJetMicroDst::create -> define TFile for output\n");
    return kFALSE;
  }
  fFile = file;
  fFile->cd();
  fTree = new TTree(nameTree.Data(),"Temporary micro DST for jets analysis");
  // partons 
  fTree->Branch("npart", &fNpart, "npart/I");
  fTree->Branch("xpt",  fXpt,  "xpt[npart]/F");
  fTree->Branch("xeta", fXeta, "xeta[npart]/F");
  fTree->Branch("xphi", fXphi, "xphi[npart]/F");
  // jets 
  fTree->Branch("njet", &fNjet, "njet/I");
  fTree->Branch("jet",  fJet,  "jet[njet]/F");
  fTree->Branch("jetaw", fJetaW, "jetaw[njet]/F");
  fTree->Branch("jphiw", fJphiW, "jphiw[njet]/F");
  fTree->Branch("jetal", fJetaL, "jetal[njet]/F");
  fTree->Branch("jphil", fJphiL, "jphil[njet]/F");
  return kTRUE;
}

Bool_t AliEMCALJetMicroDst::Create(const char *fname)
{
  TFile *file = new TFile(fname, "RECREATE");
  if(file) {
    fName = fname;
    return Create(file);
  } else     return kFALSE;
}

Bool_t AliEMCALJetMicroDst::Open(const char *fname)
{
  if(fFile && fFile->IsOpen()) fFile->Close();
  if(strlen(fname)) fName = fname;
  TFile *file = new TFile(fName.Data(), "READ");
  if(file) {
    printf("<I> open file %s \n",fName.Data()); 
    return Initialize(file);
  } else {
    printf("<E> can not open file %s \n",fName.Data()); 
    return kFALSE;
  }
}

Bool_t AliEMCALJetMicroDst::Initialize(TFile *file)
{
  if(file) fFile = file; 
  fFile->cd();
  fTree = (TTree*)fFile->Get(nameTree.Data());
  if(!fTree) return kFALSE;
  // partons
  fTree->SetBranchAddress("npart",&fNpart);
  fTree->SetBranchAddress("xpt",   fXpt);
  fTree->SetBranchAddress("xeta",  fXeta);
  fTree->SetBranchAddress("xphi",  fXphi);
  // jets 
  fTree->SetBranchAddress("njet", &fNjet);
  fTree->SetBranchAddress("jet",   fJet);
  fTree->SetBranchAddress("jetaw", fJetaW);
  fTree->SetBranchAddress("jphiw", fJphiW);
  fTree->SetBranchAddress("jetal", fJetaL);
  fTree->SetBranchAddress("jphil", fJphiL);

  return kTRUE;
}

void AliEMCALJetMicroDst::Print(Option_t* option) const
{
  if(option);
  if(fFile) {
    fFile->Print();
    if(fTree) fTree->Print();
    else printf("<I> TRee is zero\n");
  } else printf("<I> File with TRee is zero\n");
  printf("Default name of file %s\n", fName.Data());

  printf("\n #partons %2i ", fNpart);
  for(Int_t i=0; i<fNpart; i++){
    printf("\n     %1i Pt %7.1f eta  %f7.4 phi  %f7.4 ",
    i, fXpt[i],fXeta[i],fXphi[i]);  
  }

  printf("\n #jets    %2i ", fNjet);
  for(Int_t i=0; i<fNjet; i++){
    printf("\n     %1i Et %7.1f etaW %f7.4 phiW %f7.4 ",
    i,fJet[i],fJetaW[i],fJphiW[i]);  
  }
}

void AliEMCALJetMicroDst::Fill(AliRun *run, AliEMCALJetFinder* jetFinder)
{
  if(!run) run = gAlice;
  AliGenEventHeader* evHeader = run->GetHeader()->GenEventHeader();
  TString tmp(evHeader->ClassName());
  if(tmp.Contains("Hijing")) {
     AliGenHijingEventHeader *hijEvHeader = (AliGenHijingEventHeader*)evHeader;
     FillPartons(hijEvHeader);
  } else if(tmp.Contains("Pythia")) {
     FillPartons();     
  } else {
    printf("\n <E> Wrong type of generator -> %s ",tmp.Data()); 
    printf("\n     Inof about partons will be absent "); 
  }

  FillJets(jetFinder);
  fTree->Fill();
}

void AliEMCALJetMicroDst::FillPartons(AliGenHijingEventHeader *header)
{
  TLorentzVector parton[4];
  header->GetJets(parton[0], parton[1], parton[2], parton[3]);

  fNpart = 4; // 
  for(Int_t i=0; i<4; i++){
     fXpt[i]  = parton[i].Pt();
     fXeta[i] = parton[i].Eta();
     fXphi[i] = parton[i].Phi();
     if(i>1) {
       hEtaPhiPart->Fill(parton[i].Eta(), parton[i].Phi());
       hPtPart->Fill(parton[i].Pt());
     }
  }
}

void AliEMCALJetMicroDst::FillPartons()
{// for case of Pythia -> get info from full event record

  fNpart = 2;
  TParticle *MPart;
  Int_t ind;
  for(Int_t i=6; i<8; i++){
     MPart = gAlice->Particle(i);
     ind   = i-6;
     fXpt[ind]  = MPart->Pt();
     fXeta[ind] = MPart->Eta();
     fXphi[ind] = MPart->Phi();

     hEtaPhiPart->Fill(fXeta[ind], fXphi[ind]);
     hPtPart->Fill(fXpt[ind]);
  }
}

void AliEMCALJetMicroDst::FillJets(AliEMCALJetFinder* jetFinder)
{
  fNjet = 0;
  if(fDebug>1) printf("\n<I> AliEMCALJetMicroDst::FillJets"); 
  if(!jetFinder) {
    if(fDebug>1) printf("\n : jetFinder is zero"); 
    return;
  }
  fNjet = jetFinder->Njets();
  if(fNjet>10) {
    if(fDebug>1) printf("\n <W> wrong value of jetFinder->Njets() %i ", fNjet); 
    fNjet = 10;
  }
  hNJet->Fill(fNjet);
  if(fDebug>1) printf("\n <I> fNjet %i", fNjet); 
  if(fNjet){
    for(Int_t i=0; i<fNjet; i++){
      fJet[i]   = jetFinder->JetEnergy(i);
      fJetaW[i] = jetFinder->JetEtaW(i);
      fJphiW[i] = jetFinder->JetPhiW(i);
      fJetaL[i] = jetFinder->JetEtaL(i);
      fJphiL[i] = jetFinder->JetPhiL(i);
      hPtJet->Fill(fJet[i]);
      hEtaPhiJet->Fill(fJetaW[i],fJphiW[i]); 
    }
  }
}

Int_t AliEMCALJetMicroDst::GetEntry(Int_t entry)
{ // Read contents of entry.
   if (!fTree) {
      printf("\n<E> AliEMCALJetMicroDst::GetEntry() -> define TTree \n");
      return -1;
   }
   return fTree->GetEntry(entry);
}

Bool_t AliEMCALJetMicroDst::GetParton(Int_t i, Float_t& pt, Float_t& eta, Float_t& phi)
{
  if(i>=0 && i<fNpart) {
    pt  = fXpt[i];
    eta = fXeta[i];
    phi = fXphi[i];
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
  if(i>=0 && i<fNjet) {
    pt  = fJet[i];
    if(mode==1) {
      eta = fJetaW[i];
      phi = fJphiW[i];
    } else {
      eta = fJetaL[i];
      phi = fJphiL[i];
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
    printf("\n<I> AliEMCALJetMicroDst::Test() -> define file with proper TTree !");
    return;
  }
  Int_t nbytes=0, nb=0, nentries=Int_t(fTree->GetEntries());
  for(Int_t i=0; i<nentries; i++){
    nb = fTree->GetEntry(i);  
    nbytes += nb;
    for(Int_t j=0; j<fNpart; j++){
      hEtaPhiPart->Fill(fXeta[j], fXphi[j]);
      hPtPart->Fill(fXpt[j]);
    }

    hNJet->Fill(fNjet);
    if(fNjet){
      for(Int_t j=0; j<fNjet; j++) {
        hPtJet->Fill(fJet[j]);
        hEtaPhiJet->Fill(fJetaW[j],fJphiW[j]); 
      }
    }
  }
  printf("\n<I> AliEMCALJetMicroDst::Test() -> Entries %5i Bytes %10i\n", nentries, nb);
}

void AliEMCALJetMicroDst::FillVector(Float_t pt, Float_t eta, Float_t phi, TVector3& vec)
{ // service function 
  static Float_t px, py, pz;

  px = pt*TMath::Cos(phi);
  py = pt*TMath::Sin(phi);
  pz = pt*TMath::SinH(eta);  // sinh(eta) = cot(theta)

  vec.SetXYZ(px, py, pz);
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
