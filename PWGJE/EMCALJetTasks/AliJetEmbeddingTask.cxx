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
  AliJetModelBaseTask("AliJetEmbeddingTask", kFALSE),
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
  fPathTreeinputFile("TreesJet4VectorForEmbedding.root"),
  fTreeinputName("fTreeJetA"),
  fBranchJDetName("fJetGenSub"),
  fTreeJet4Vect(0),
  fCurrentEntry(0)
{
  // Default constructor.
  SetSuffix("Embedded");
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask(const char *name) : 
  AliJetModelBaseTask(name, kFALSE),
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
  fPathTreeinputFile("TreesJet4VectorForEmbedding.root"),
  fTreeinputName("fTreeJetA"),
  fBranchJDetName("fJetGenSub"),
  fTreeJet4Vect(0),
  fCurrentEntry(0)
{
  // Standard constructor.
  SetSuffix("Embedded");
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliJetEmbeddingTask::~AliJetEmbeddingTask()
{
  // Destructor
}

//________________________________________________________________________

void AliJetEmbeddingTask::UserCreateOutputObjects(){
   
   AliJetModelBaseTask::UserCreateOutputObjects();
   
   fOutput = new TList();
   fOutput->SetOwner();
   
   if(!fPathTreeinputFile.IsNull()){
      SetTreeFromFile(fPathTreeinputFile, fTreeinputName);
      if(!fTreeJet4Vect) AliFatal("Something went wrong in setting the tree");
      fOutput->Add(fTreeJet4Vect);
   }
   
   if(!fPathMinputFile.IsNull() && fPathpTinputFile.IsNull()){
      SetMassDistributionFromFile(fPathMinputFile,fMinputName);
      if(!fHMassDistrib) AliFatal("Something went wrong in setting the M distribution");
      fOutput->Add(fHMassDistrib);
   }
   
   if(fPathMinputFile.IsNull() && !fPathpTinputFile.IsNull()){
      SetpTDistributionFromFile(fPathpTinputFile, fpTinputName);
      if(!fPtSpectrum) AliFatal("Something went wrong in setting the pT distribution");
       fOutput->Add(fPtSpectrum);
   }
   
   if(!fPathMinputFile.IsNull() && !fPathpTinputFile.IsNull()){
      SetMassAndPtDistributionFromFile(fPathMinputFile, fPathpTinputFile, fMinputName, fpTinputName);
      if(!fHMassDistrib) AliFatal("Something went wrong in setting the M distribution");
      if(!fPtSpectrum) AliFatal("Something went wrong in setting the pT distribution");
      fOutput->Add(fHMassDistrib);
      fOutput->Add(fPtSpectrum);
   }
   
   PostData(1, fOutput);
   
   
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
       	  Printf("Adding track from tree");
       	  TLorentzVector *jetDet = 0;
       	  TBranch *bDet = 0;
       	  Int_t nbranches = fTreeJet4Vect->GetNbranches();
       	  Printf("Input: %p, %d , %s, %p %p", fTreeJet4Vect,nbranches, fBranchJDetName.Data(), jetDet,bDet );
       	  fTreeJet4Vect->SetBranchAddress(fBranchJDetName.Data(), &jetDet, &bDet);
       	  Printf("Branch set, %p",bDet);
       	  Int_t nentries = fTreeJet4Vect->GetEntries();
       	  if(fCurrentEntry < nentries) bDet->GetEntry(fCurrentEntry);
       	  else {
       	     fCurrentEntry = 0;
       	     AliWarning("Starting from first entry again");
       	     bDet->GetEntry(fCurrentEntry);
       	  }
       	  fCurrentEntry++;
       	  
       	  AddTrack(jetDet->Pt(), jetDet->Eta(), jetDet->Phi(), 0,0,0,0, kFALSE, 0, charge, jetDet->M());
       	  
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
   fTreeJet4Vect = (TTree*)tree->Clone(Form("%sCp", tree->GetName()));
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
      Printf("Tree not found!!!");
      f->ls();
      return;
   }
   SetTree(tree);
   
   //f->Close();
   //delete f;

   return;
}
