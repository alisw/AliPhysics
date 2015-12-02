// $Id$
//
// Jet embedding task.
//
// Author: S.Aiola, C.Loizides

#include <TRandom3.h>
#include <TFile.h>
#include <AliLog.h>
#include <TGrid.h>

#include "AliJetEmbeddingTask.h"

ClassImp(AliJetEmbeddingTask)

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask() : 
  AliJetModelBaseTask("AliJetEmbeddingTask", kTRUE),
  fMassless(kFALSE),
  fMassFromDistr(kFALSE),
  fNeutralFraction(0),
  fNeutralMass(0.135),
  fMass(0.1396),
  fHMassDistrib(0),
  fPathMinputFile(""),
  fPathpTinputFile(""),
  fMinputName(""),
  fpTinputName(""),
  fFromTree(0),
  fPathTreeinputFile(""),
  fTreeinputName("fTreeJet"),
  fBranchJDetName("fJetDet"),
  fTreeJet4Vect(0),
  fCurrentEntry(0), 
  fInput(0),
  fRandomEntry(0),
  fMaxMemory(2000000000)
{
  // Default constructor.
  SetSuffix("Embedded");
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask(const char *name) : 
  AliJetModelBaseTask(name, kTRUE),
  fMassless(kFALSE),
  fMassFromDistr(kFALSE),
  fNeutralFraction(0),
  fNeutralMass(0.135),
  fMass(0.1396),
  fHMassDistrib(0),
  fPathMinputFile(""),
  fPathpTinputFile(""),
  fMinputName(""),
  fpTinputName(""),
  fFromTree(0),
  fPathTreeinputFile(""),
  fTreeinputName("fTreeJet"),
  fBranchJDetName("fJetDet"),
  fTreeJet4Vect(0),
  fCurrentEntry(0),
  fInput(0),
  fRandomEntry(0),
  fMaxMemory(2000000000)
{
  // Standard constructor.
  SetSuffix("Embedded");
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
AliJetEmbeddingTask::~AliJetEmbeddingTask()
{
  // Destructor
}

//________________________________________________________________________

void AliJetEmbeddingTask::UserCreateOutputObjects(){
   
   AliJetModelBaseTask::UserCreateOutputObjects();
   
   OpenFile(1);
   fInput = new TList();
   fInput->SetOwner();
   
   if(!fPathTreeinputFile.IsNull()){
      SetTreeFromFile(fPathTreeinputFile, fTreeinputName);
      if(!fTreeJet4Vect) AliFatal("Something went wrong in setting the tree");
      Printf("Emb task");
      TLorentzVector *detjet = 0;
      fTreeJet4Vect->SetBranchAddress(fBranchJDetName, &detjet);
      fTreeJet4Vect->GetEntry(0);
      fTreeJet4Vect->Show();
      if(fRandomEntry) fTreeJet4Vect->LoadBaskets(fMaxMemory); //default 2GB (2000000000)
      else{
      	 if(fPtMin != 0 && fPtMax != 0){
      	    AliInfo(Form("Using range %.2f - %.2f GeV/c", fPtMin, fPtMax));
      	    for(Int_t i = 0; i<fTreeJet4Vect->GetEntries(); i++){
      	       fTreeJet4Vect->GetEntry(i);
      	       if((detjet->Pt()> fPtMin) && (detjet->Pt()> fPtMax)) {
      	       	  fCurrentEntry = i;
      	       	  AliInfo(Form("Setting fCurrentEntry to %d: will loop from there", fCurrentEntry));
      	       	  
      	       	  break;
      	       }
      	    }
      	 }
      }
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
       
       // Add track from tree of 4-vectors (jet reco) and save the particle level somewhere
       if(fFromTree){
       	  if(!fTreeJet4Vect || fBranchJDetName.IsNull()) {
       	     AliFatal(Form("Tree or branch name not found"));
       	  }
       	  TLorentzVector *jetDet = 0;
       	  TBranch *bDet = 0;
       	  fTreeJet4Vect->ResetBranchAddresses();
       	  fTreeJet4Vect->SetBranchAddress(fBranchJDetName.Data(), &jetDet, &bDet);
       	  Int_t nentries = fTreeJet4Vect->GetEntries();
       	  Double_t pTemb = -1;
       	  while(pTemb < fMinPtEmb){
       	     if(fRandomEntry) fCurrentEntry = gRandom->Integer(nentries);
       	     if(fCurrentEntry < nentries) bDet->GetEntry(fCurrentEntry);
       	     else {
       	     	fCurrentEntry = 0;
       	     	AliWarning("Starting from first entry again");
       	     	bDet->GetEntry(fCurrentEntry);
       	     }
       	     pTemb = jetDet->Pt();
       	     //Printf("Embedding entry %d -> Det Lev %.2f, %.2f, %.2f, %.2f", fCurrentEntry, jetDet->Pt(), jetDet->Phi(), jetDet->Eta(), jetDet->M());

       	     if ((fRandomEntry && (pTemb >= fMinPtEmb)) || (!fRandomEntry && ((fPtMin == 0 && fPtMax == 0) || (pTemb > fPtMin && pTemb < fPtMax)))) {
       	     	
       	     	AddTrack(jetDet->Pt(), jetDet->Eta(), jetDet->Phi(), 0,0,0,0, kFALSE, fCurrentEntry, charge, jetDet->M());
       	     	//Printf("Embedded (pT = %.2f)!!!!", jetDet->Pt());
       	     }
       	     if(!fRandomEntry) fCurrentEntry++;
       	  }
       } else {
       	  
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
       	  if(fMassFromDistr) {
       	     if(fHMassDistrib)
       	     	mass = fHMassDistrib->GetRandom();
       	     else {
       	     	AliError(Form("Template distribution for mass of track embedding not found, use %f", fMass));
       	     	mass = fMass;
       	     }
       	  }
       	  AddTrack(-1,-999,-1,0,0,0,0,kFALSE,0,charge,mass);
       }
    }
  }
  
  AliJetModelBaseTask::FillHistograms();
}

//________________________________________________________________________

void AliJetEmbeddingTask::SetMassDistribution(TH1F *hM)  {
   if(!hM){
      AliError("Null histogram for mass distribution");
      return;
   }
   fMassFromDistr = kTRUE; 
   fHMassDistrib = hM;
   AliInfo("Input mass distribution set");
   
   return;
}
//________________________________________________________________________
void AliJetEmbeddingTask::SetMassDistributionFromFile(TString filename, TString histoname){
   
   if(filename.Contains("alien")) {
      TGrid::Connect("alien://");
   }
   TFile *f = TFile::Open(filename);
   if(!f){
      AliFatal(Form("File %s not found, cannot SetMassDistribution", filename.Data()));
      return;
   }
   
   TH1F* h = dynamic_cast<TH1F*> (f->Get(histoname));
   if(!h) {
      AliError("Input file for Mass not found");
      f->ls();
   }
   SetMassDistribution(h);
   
   //f->Close();
   //delete f;
   
   return;

}

//________________________________________________________________________

void AliJetEmbeddingTask::SetMassAndPtDistributionFromFile(TString filenameM, TString filenamepT, TString histonameM, TString histonamepT){
   SetMassDistributionFromFile(filenameM, histonameM);
   SetpTDistributionFromFile(filenamepT, histonamepT);
   return;
}

//________________________________________________________________________
void AliJetEmbeddingTask::SetpTDistributionFromFile(TString filename, TString histoname){
   
   if(filename.Contains("alien")) {
      TGrid::Connect("alien://");
   }
   TFile *f = TFile::Open(filename);
   if(!f){
      AliFatal(Form("File %s not found, cannot SetpTDistribution", filename.Data()));
      return;
   }

   TH1F* h = dynamic_cast<TH1F*> (f->Get(histoname));
   if(!h) {
      AliError("Input file for pT not found");
      f->ls();
   }

   AliJetModelBaseTask::SetPtSpectrum(h);

   //f->Close();
   //delete f;

   return;

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

