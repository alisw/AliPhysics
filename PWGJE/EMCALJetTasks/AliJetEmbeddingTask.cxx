// $Id$
//
// Jet embedding task.
//
// Author: S.Aiola, C.Loizides

#include <TRandom3.h>
#include <TFile.h>
#include <AliLog.h>
#include <TGrid.h>
#include <THnSparse.h>
#include <AliPicoTrack.h>

#include "AliJetEmbeddingTask.h"

ClassImp(AliJetEmbeddingTask)

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask() : 
  AliJetModelBaseTask("AliJetEmbeddingTask", kTRUE),
  fMassless(kFALSE),
  fNeutralFraction(0),
  fNeutralMass(0.135),
  fMass(0.1396),
  fPathMinputFile(""),
  fPathpTinputFile(""),
  fMinputName(""),
  fpTinputName(""),
  fFromTree(0),
  fPathTreeinputFile(""),
  fTreeinputName("fTreeJet"),
  fBranchJDetName("fJetDet."),
  fBranchJParName("fJetPart."),
  fTreeJet4Vect(0),
  fCurrentEntry(0), 
  fInput(0),
  fRandomEntry(0),
  fNBins(10),
  fDownscale(0),
  fPtRanges(0),
  fNevPerBin(0),
  fCount(0),
  fCurrentBin(0),
  fGoBack(0),
  fhPartJet(0),
  fhEtaPart(0),
  fhPhiPart(0),
  fhTreeEntriesUsed(0),
  fNTreeEntries(0),
  fOldCutOnDetPt(0)
  
{
  // Default constructor.
  SetSuffix("Embedded");
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask(const char *name) : 
  AliJetModelBaseTask(name, kTRUE),
  fMassless(kFALSE),
  fNeutralFraction(0),
  fNeutralMass(0.135),
  fMass(0.1396),
  fPathMinputFile(""),
  fPathpTinputFile(""),
  fMinputName(""),
  fpTinputName(""),
  fFromTree(0),
  fPathTreeinputFile(""),
  fTreeinputName("fTreeJet"),
  fBranchJDetName("fJetDet."),
  fBranchJParName("fJetPart."),
  fTreeJet4Vect(0),
  fCurrentEntry(0),
  fInput(0),
  fRandomEntry(0),
  fNBins(10),
  fDownscale(0),
  fPtRanges(0),
  fNevPerBin(0),
  fCount(0),
  fCurrentBin(0),
  fGoBack(0),
  fhPartJet(0),
  fhEtaPart(0),
  fhPhiPart(0),
  fhTreeEntriesUsed(0),
  fNTreeEntries(0),
  fOldCutOnDetPt(0)
{
  // Standard constructor.
  SetSuffix("Embedded");
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
AliJetEmbeddingTask::~AliJetEmbeddingTask()
{
  // Destructor
  delete[] fDownscale; fDownscale = 0;
  delete[] fPtRanges; fPtRanges = 0;
}

//________________________________________________________________________

void AliJetEmbeddingTask::UserCreateOutputObjects(){
   
   AliJetModelBaseTask::UserCreateOutputObjects();
   
   OpenFile(1);
   fInput = new TList();
   fInput->SetOwner();
   
   delete gRandom;
   gRandom = new TRandom3(0); //random seed
   
   if(!fPathTreeinputFile.IsNull()){
      SetTreeFromFile(fPathTreeinputFile, fTreeinputName);
      if(!fTreeJet4Vect) AliFatal("Something went wrong in setting the tree");
      Printf("Emb task");
      TLorentzVector *detjet = 0;
      fTreeJet4Vect->SetBranchAddress(fBranchJDetName, &detjet);
      fTreeJet4Vect->GetEntry(0);
      fTreeJet4Vect->Show();
      fNTreeEntries = fTreeJet4Vect->GetEntries();
      fhTreeEntriesUsed = new TH1F("fhTreeEntriesUsed", "Entries;Entry in TTree", fNTreeEntries, 0, fNTreeEntries-1);
      fOutput->Add(fhTreeEntriesUsed);
      
      fCurrentEntry = gRandom->Integer(fNTreeEntries); //in each worker it starts from a different entry
      Printf("Entries %lld, start from %d", fNTreeEntries, fCurrentEntry);
   }
   
   if(!fPathMinputFile.IsNull() && fPathpTinputFile.IsNull()){
      SetMassDistributionFromFile(fPathMinputFile,fMinputName);
      if(!fHMassDistrib) AliFatal("Something went wrong in setting the M distribution");
      fInput->Add(fHMassDistrib);
   }
   
   if(fPathMinputFile.IsNull() && !fPathpTinputFile.IsNull()){
      SetpTDistributionFromFile(fPathpTinputFile, fpTinputName);
      if(!fPtSpectrum) AliFatal("Something went wrong in setting the pT distribution");
       fInput->Add(fPtSpectrum);
   }
   
   if(!fPathMinputFile.IsNull() && !fPathpTinputFile.IsNull()){
      SetMassAndPtDistributionFromFile(fPathMinputFile, fPathpTinputFile, fMinputName, fpTinputName);
      if(!fHMassDistrib) AliFatal("Something went wrong in setting the M distribution");
      if(!fPtSpectrum) AliFatal("Something went wrong in setting the pT distribution");
      fInput->Add(fHMassDistrib);
      fInput->Add(fPtSpectrum);
   }
   PostData(2, fInput);
   
   const Int_t ndim = 4;
   Int_t nbins = 80, mind = -30, maxd = 30, minm = 0, maxm = 80, minpt = 0, maxpt = 120;
   const Int_t nbinshnsp[ndim] = {nbins, nbins, nbins, nbins};  
   const Double_t minhnsp[ndim] = {(Double_t)mind, (Double_t)mind, (Double_t)minm, (Double_t)minpt};
   const Double_t maxhnsp[ndim] = {(Double_t)maxd, (Double_t)maxd, (Double_t)maxm, (Double_t)maxpt};
   TString title = "Check part level;#it{M}_{det} - #it{M}_{par} (GeV); #it{p}_{T, det} - #it{p}_{T, par} (GeV/#it{c}); #it{M} (GeV); #it{p}_{T} (GeV/#it{c})";
   
   fhPartJet = new THnSparseF("fhPartJet", title.Data(), ndim, nbinshnsp, minhnsp, maxhnsp);
   fhPartJet->Sumw2();
   fOutput->Add(fhPartJet);
   
   fhEtaPart = new TH1F("fhEtaPart" , "#eta distribution part; #eta", 100, -1, 1);
   fhEtaPart->Sumw2();
   fOutput->Add(fhEtaPart);
   fhPhiPart = new TH1F("fhPhiPart" , "#varphi distribution; #varphi", 100, (-1)*TMath::Pi(), TMath::Pi());
   fhPhiPart->Sumw2();
   fOutput->Add(fhPhiPart);
   
}

//________________________________________________________________________
void AliJetEmbeddingTask::Run() 
{
  // Embed particles.
  if (fNClusters > 0 && fOutClusters) {
    if (fCopyArray) 
      CopyClusters();
    for (Int_t i = 0; i < fNClusters; ++i) {
      AddCluster();
    }
  }
 
  if (fNTracks > 0 && fOutTracks) {
    if (fCopyArray) 
      CopyTracks();
    for (Int_t i = 0; i < fNTracks; ++i) {

       Short_t charge = 1;
       if(fNeutralFraction>0.) {
       	  Double_t rnd = gRandom->Rndm();
       	  if(rnd<fNeutralFraction) {
       	     charge = 0;
       	     //mass = fNeutralMass;
       	  }
       }
       // Add track from tree of 4-vectors (jet reco) and save the particle level somewhere
       if(fFromTree){
       	  
       	  if(!fTreeJet4Vect || fBranchJDetName.IsNull()) {
       	     AliFatal(Form("Tree or branch name not found"));
       	  }
       	  TLorentzVector *jetDet = 0;
       	  TLorentzVector *jetPar = 0;
       	  fTreeJet4Vect->ResetBranchAddresses();
       	  //fTreeJet4Vect->SetBranchAddress(fBranchJDetName.Data(), &jetDet, &bDet);
       	  fTreeJet4Vect->SetBranchAddress(fBranchJDetName.Data(), &jetDet);
       	  fTreeJet4Vect->SetBranchAddress(fBranchJParName.Data(), &jetPar);
       	  
       	  Double_t pTemb = -1, pTpar = -1;
       	  if(fCurrentEntry < fNTreeEntries) fTreeJet4Vect->GetEntry(fCurrentEntry);
       	  else {
       	     fCurrentEntry = 0;
       	     AliWarning("Starting from first entry again");
       	     fTreeJet4Vect->GetEntry(fCurrentEntry);
       	  }
       	  
       	  pTpar = jetPar->Pt(); //cut on the particle level pT
       	  if(fOldCutOnDetPt) pTpar = jetDet->Pt(); //cut on the detector level pT (as it was before, keeping the functionality)
       	  // selected pT range 
       	  if((fPtMin != 0 && fPtMax != 0) && !fRandomEntry) {
       	     Int_t countWhile = 0;
       	     while(!(pTpar > fPtMin && pTpar < fPtMax) && countWhile < fNTreeEntries){
       	     	fCurrentEntry++;
       	     	if(fCurrentEntry < fNTreeEntries) fTreeJet4Vect->GetEntry(fCurrentEntry);
       	     	else {
       	     	   fCurrentEntry = 0;
       	     	   AliWarning("Starting from first entry again");
       	     	   fTreeJet4Vect->GetEntry(fCurrentEntry);
       	     	}
       	     	pTpar = jetPar->Pt();  //cut on the particle level pT
       	     	if(fOldCutOnDetPt) pTpar = jetDet->Pt(); //cut on the detector level pT (as it was before, keeping the functionality)
       	     	countWhile++;
       	     }
       	  }
       	  // exclude a fraction of the entries -- doesn't work very well
       	  //if(fRandomEntry){
          //
     	  //   if(fCurrentEntry < nentries) fTreeJet4Vect->GetEntry(fCurrentEntry);
       	  //   else {
       	  //   	fCurrentEntry = 0;
       	  //   	AliWarning("Starting from first entry again");
       	  //   	fTreeJet4Vect->GetEntry(fCurrentEntry);
       	  //   }
       	  //   pTemb = jetDet->Pt();
       	  //   
       	  //   Float_t downscl = GetDownscalinigFactor();
       	  //   Float_t random = gRandom->Rndm();
       	  //   
       	  //   while (random > downscl){
       	  //   	fCurrentEntry++;
       	  //   	random = gRandom->Rndm();
       	  //   	if(fCurrentEntry < nentries) fTreeJet4Vect->GetEntry(fCurrentEntry);
       	  //   	else {
       	  //   	   fCurrentEntry = 0;
       	  //   	   AliWarning("Starting from first entry again");
       	  //   	   fTreeJet4Vect->GetEntry(fCurrentEntry);
       	  //   	}
       	  //   	pTemb = jetDet->Pt();
       	  //   	if(pTemb < fPtRanges[fCurrentBin]) {
       	  //   	   random = gRandom->Rndm();
       	  //   	   
       	  //   	}
       	  //   	   
       	  //   }
       	  //   
       	  //}

       	  fhTreeEntriesUsed->Fill(fCurrentEntry);
       	  // Add the track that complies with the settings 
       	  AddTrack(jetDet->Pt(), jetDet->Eta(), jetDet->Phi(),0,0,0,0,kFALSE,  fCurrentEntry, charge, jetDet->M());
       	  //Printf("Embedded det %.2f, part %.2f", jetDet->Pt(), jetPar->Pt());
       	  
       	  Double_t x[4] = {jetDet->M() - jetPar->M(), jetDet->Pt() - jetPar->Pt(), jetPar->M(), jetPar->Pt()};
       	  fhPartJet->Fill(x);
       	  
       	  fhEtaPart->Fill(jetPar->Eta());
       	  fhPhiPart->Fill(jetPar->Phi());
       	  fCount++; // count the number of track embedded in the current pT range
       	  fCurrentEntry++; //increase for next iteration
       	  
       } else { //other inputs
       	  
       	  Double_t mass = fMass;
       	  Short_t charge = 1;
       	  if(fNeutralFraction>0.) {
       	     Double_t rnd = gRandom->Rndm();
       	     if(rnd<fNeutralFraction) {
       	     	charge = 0;
       	     	mass = fNeutralMass;
       	     }
       	  }
       	  if(fMassless) mass = 0.;
       	  if(fMassFromDistr) mass = -999;
       	  
       	  AddTrack(-999,-999,-999,0,0,0,0,kFALSE,0,charge,mass);
       }
    }
  }
  
  AliJetModelBaseTask::FillHistograms();
}

//________________________________________________________________________

Float_t AliJetEmbeddingTask::GetDownscalinigFactor(){
   
   if(fCount > fNevPerBin) {
      //Printf("%d = fCount, Increasing fCurrentBin %d -> %d", fCount, fCurrentBin, fCurrentBin+1);
      fCurrentBin++;
      fCount = 0;
      
   }

   if (fCurrentBin >= fNBins) {
      AliError(Form("Bin %d out of bound %d, set to fNBins - 1 = %d", fCurrentBin, fNBins, fNBins - 1));
      fCurrentBin = fNBins - fGoBack;
      if (fCurrentBin < 0) fCurrentBin = fNBins - 1;
      fGoBack++;
   }
   return fDownscale[fCurrentBin];
}


//________________________________________________________________________
void AliJetEmbeddingTask::SetTree(TTree *tree)  {
   if(!tree){
      AliError("Null tree");
      return;
   }
   fFromTree = kTRUE; 
   fTreeJet4Vect = (TTree*)tree->Clone(Form("%sCpEmb", tree->GetName()));
   AliInfo(Form("Input tree set %d (%p -> %p)", fTreeJet4Vect->GetNbranches(), tree, fTreeJet4Vect));
   
   return;
}

//________________________________________________________________________

void AliJetEmbeddingTask::SetTreeFromFile(TString filename, TString treename){
   
   if(filename.Contains("alien")) {
      TGrid::Connect("alien://");
   }
   TFile *f = TFile::Open(filename);
   if(!f){
      Printf("File %s not found, cannot SetTree", filename.Data());
      return;
   }
   
   TTree *tree = dynamic_cast<TTree*>(f->Get(treename));
   if(!tree){
      Printf("Tree %s not found!!!", treename.Data());
      f->ls();
      return;
   }
   SetTree(tree);
   
   //f->Close();
   //delete f;

   return;
}

//________________________________________________________________________
  
void AliJetEmbeddingTask::SetRejection(Float_t* rej) {
   
   if(fNBins == 0){
      AliError("Set number of bins first! SetNBinsEmbedding(nbins);");
      return;
   }
   Printf("Creating array of factors with size %d", fNBins);
   fDownscale = new Float_t[fNBins];
   for(Int_t i = 0; i<fNBins; i++) {
      fDownscale[i] = rej[i];
      Printf("Bin %d -> Factor = %e", i, fDownscale[i]);
   }
   return;
}


void AliJetEmbeddingTask::Terminate(){
   Printf("fGoBack = %d", fGoBack);
}

