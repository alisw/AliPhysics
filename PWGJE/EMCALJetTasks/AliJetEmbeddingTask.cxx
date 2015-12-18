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
  fhDeltapT(0),
  fhDeltaM(0),
  fhpTPart(0),
  fhMPart(0),
  fhEtaPart(0),
  fhPhiPart(0)
  
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
  fhDeltapT(0),
  fhDeltaM(0),
  fhpTPart(0),
  fhMPart(0),
  fhEtaPart(0),
  fhPhiPart(0)
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
   
   if(!fPathTreeinputFile.IsNull()){
      SetTreeFromFile(fPathTreeinputFile, fTreeinputName);
      if(!fTreeJet4Vect) AliFatal("Something went wrong in setting the tree");
      Printf("Emb task");
      TLorentzVector *detjet = 0;
      fTreeJet4Vect->SetBranchAddress(fBranchJDetName, &detjet);
      fTreeJet4Vect->GetEntry(0);
      fTreeJet4Vect->Show();
      if(!fRandomEntry && fPtMin != 0 && fPtMax != 0){
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
   
   

  
  fhDeltapT = new TH1F("fhDeltapT", "Delta #it{p}_{T}; #Delta #it{p}_{T} (GeV/#it{c})", 60, -30., 30.);
  fhDeltapT->Sumw2();
  fOutput->Add(fhDeltapT);
  
  fhDeltaM = new TH1F("fhDeltaM", "Delta M; #Delta M (GeV)", 60, -30., 30.);
  fhDeltaM->Sumw2();
  fOutput->Add(fhDeltaM);

  fhpTPart  = new TH1F("fhpTPart"  , "#it{p}_{T, part};#it{p}_{T} (GeV/#it{c})", 60, -10, 100);
  fhpTPart->Sumw2();
  fOutput->Add(fhpTPart);
  fhMPart   = new TH1F("fhMPart"   , "#it{M}_{part};#it{M} (GeV)", 60, -10, 20);
  fhMPart->Sumw2();
  fOutput->Add(fhMPart);
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
       	  Int_t nentries = fTreeJet4Vect->GetEntries();
       	  Double_t pTemb = -1;
       	  if(fCurrentEntry < nentries) fTreeJet4Vect->GetEntry(fCurrentEntry);
       	  else {
       	     fCurrentEntry = 0;
       	     AliWarning("Starting from first entry again");
       	     fTreeJet4Vect->GetEntry(fCurrentEntry);
       	  }
       	  pTemb = jetDet->Pt();
       	  
       	  // selected pT range 
       	  if((fPtMin != 0 && fPtMax != 0) && !fRandomEntry) {
       	     while(!(pTemb > fPtMin && pTemb < fPtMax)){
       	     	fCurrentEntry++;
       	     	if(fCurrentEntry < nentries) fTreeJet4Vect->GetEntry(fCurrentEntry);
       	     	else {
       	     	   fCurrentEntry = 0;
       	     	   AliWarning("Starting from first entry again");
       	     	   fTreeJet4Vect->GetEntry(fCurrentEntry);
       	     	}
       	     	pTemb = jetDet->Pt();
       	     }
       	  }
       	  // exclude a fraction of the entries -- doesn't work very well
       	  if(fRandomEntry){
  
     	     if(fCurrentEntry < nentries) fTreeJet4Vect->GetEntry(fCurrentEntry);
       	     else {
       	     	fCurrentEntry = 0;
       	     	AliWarning("Starting from first entry again");
       	     	fTreeJet4Vect->GetEntry(fCurrentEntry);
       	     }
       	     pTemb = jetDet->Pt();
       	     
       	     Float_t downscl = GetDownscalinigFactor();
       	     Float_t random = gRandom->Rndm();
       	     
       	     while (random > downscl){
       	     	fCurrentEntry++;
       	     	random = gRandom->Rndm();
       	     	if(fCurrentEntry < nentries) fTreeJet4Vect->GetEntry(fCurrentEntry);
       	     	else {
       	     	   fCurrentEntry = 0;
       	     	   AliWarning("Starting from first entry again");
       	     	   fTreeJet4Vect->GetEntry(fCurrentEntry);
       	     	}
       	     	pTemb = jetDet->Pt();
       	     	if(pTemb < fPtRanges[fCurrentBin]) {
       	     	   random = gRandom->Rndm();
       	     	   
       	     	}
       	     	   
       	     }
       	     
       	  }

       	  // Add the track that complies with the settings 
       	  AddTrack(jetDet->Pt(), jetDet->Eta(), jetDet->Phi(),0,0,0,0,kFALSE,  fCurrentEntry, charge, jetDet->M());
       	  //Printf("Embedded det %.2f, part %.2f", jetDet->Pt(), jetPar->Pt());
       	  fhDeltapT->Fill(jetPar->Pt() - jetDet->Pt());
       	  fhDeltaM ->Fill(jetPar->M() - jetDet->M());
       	  fhpTPart ->Fill(jetPar->Pt());
       	  fhMPart  ->Fill(jetPar->M());
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
       	  if(fMassFromDistr) {
       	     if(fHMassDistrib)
       	     	mass = fHMassDistrib->GetRandom();
       	     else {
       	     	AliError(Form("Template distribution for mass of track embedding not found, use %f", fMass));
       	     	mass = fMass;
       	     }
       	  }
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

//________________________________________________________________________

void AliJetEmbeddingTask::SetPtRangesEmb(Float_t* ptlims) {
   
   if(fNBins == 0){
      AliError("Set number of bins first! SetNBinsEmbedding(nbins);");
      return;
   }
   Printf("Creating array of pt limits with size %d", fNBins);
   fPtRanges = new Float_t[fNBins];
   for(Int_t i = 0; i<fNBins; i++){
      fPtRanges[i] = ptlims[i];
      if(fPtRanges[i] < fMinPtEmb){
      	 fPtRanges[i] = fMinPtEmb;
      	 AliError(Form("Minimum pT set to %.2f", fMinPtEmb));
      }
      Printf("Bin %d -> PtRange = %e", i, fPtRanges[i]);
   }
   return;
 }
 
void AliJetEmbeddingTask::Terminate(){
   Printf("fGoBack = %d", fGoBack);
}

