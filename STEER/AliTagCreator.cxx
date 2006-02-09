/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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

//-----------------------------------------------------------------
//           AliTagCreator class
//   This is the class to deal with the tag creation (post process)
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

//ROOT
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>

//ROOT-AliEn
#include <TGrid.h>
#include <TGridResult.h>

//AliRoot
#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliLog.h"


#include "AliTagCreator.h"


ClassImp(AliTagCreator)


//______________________________________________________________________________
AliTagCreator::AliTagCreator() //local mode
{
  //==============Default constructor for a AliTagCreator==================
  fgridpath = "";
  fSE = "";   
  fStorage = 0; 
}

//______________________________________________________________________________
AliTagCreator::~AliTagCreator()
{
//================Default destructor for a AliTagCreator=======================
}

//______________________________________________________________________________
void AliTagCreator::SetStorage(Int_t storage)
{
  // Sets correctly the storage: 0 for local, 1 for GRID
  fStorage = storage;
  if(fStorage == 0)
    AliInfo(Form("Tags will be stored locally...."));
  if(fStorage == 1)
    AliInfo(Form("Tags will be stored in the grid...."));
  if((fStorage != 0)&&(fStorage != 1))
    {
      AliInfo(Form("Storage was not properly set!!!"));
      abort();
    }  
}


//______________________________________________________________________________
Bool_t AliTagCreator::ReadESDCollection(TGridResult *fresult)
{
  // Reads the entry of the TGridResult and creates the tags
  Int_t nEntries = fresult->GetEntries();

  TString alienUrl;
  const char *guid;
  const char *md5;
  const char *turl;
  Long64_t size = -1;

  Int_t counter = 0;
  for(Int_t i = 0; i < nEntries; i++)
    {
      alienUrl = fresult->GetKey(i,"turl");
      guid = fresult->GetKey(i,"guid");
      if(fresult->GetKey(i,"size"))
	size = atol (fresult->GetKey(i,"size"));
      md5 = fresult->GetKey(i,"md5");
      turl = fresult->GetKey(i,"turl");
      if(md5 && !strlen(guid))
	md5 = 0;
      if(guid && !strlen(guid))
	guid = 0;

      TFile *f = TFile::Open(alienUrl,"READ");
      CreateTag(f,guid,md5,turl,size,counter);
      f->Close();
      delete f;	 
      counter += 1;
    }//grid result loop

  return kTRUE;
}

//__________________________________________________________________________
Bool_t AliTagCreator::MergeTags()
{
  //Merges the tags and stores the merged tag file 
  //locally if fStorage=0 or in the grid if fStorage=1
  AliInfo(Form("Merging tags....."));
  TChain *fgChain = new TChain("T");

  if(fStorage == 0) {
    const char * tagPattern = "tag";
    // Open the working directory
    void * dirp = gSystem->OpenDirectory(gSystem->pwd());
    const char * name = 0x0;
    // Add all files matching *pattern* to the chain
    while((name = gSystem->GetDirEntry(dirp))) {
      if (strstr(name,tagPattern))       
	fgChain->Add(name);  
    }//directory loop
    AliInfo(Form("Chained tag files: %d",fgChain->GetEntries()));
  }//local mode

  else if(fStorage == 1) {
    TString alienLocation = gGrid->Pwd();
    alienLocation += fgridpath.Data();
    alienLocation += "/";

    TGridResult *tagresult = gGrid->Query(alienLocation,"*tag.root","","");
    Int_t nEntries = tagresult->GetEntries();
    for(Int_t i = 0; i < nEntries; i++) {
      TString alienUrl = tagresult->GetKey(i,"turl");
      fgChain->Add(alienUrl);
    }//grid result loop      
    AliInfo(Form("Chained tag files: %d",fgChain->GetEntries()));
  }//grid mode
 
  AliRunTag *tag = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fgChain->SetBranchAddress("AliTAG",&tag);
   
  //Defining new tag objects
  AliRunTag *newTag = new AliRunTag();
  TTree ttag("T","A Tree with event tags");
  TBranch * btag = ttag.Branch("AliTAG", &newTag);
  btag->SetCompressionLevel(9);
  for(Int_t iTagFiles = 0; iTagFiles < fgChain->GetEntries(); iTagFiles++) {
    fgChain->GetEntry(iTagFiles);
    newTag->SetRunId(tag->GetRunId());
    const TClonesArray *tagList = tag->GetEventTags();
    for(Int_t j = 0; j < tagList->GetEntries(); j++) {
      evTag = (AliEventTag *) tagList->At(j);
      newTag->AddEventTag(*evTag);
    }
    ttag.Fill();
    newTag->Clear();
  }//tag file loop 
  
  TString localFileName = "Run"; localFileName += tag->GetRunId(); 
  localFileName += ".Merged"; localFileName += ".ESD.tag.root";
     
  TString alienFileName = "/alien";
  alienFileName += gGrid->Pwd();
  alienFileName += fgridpath.Data();
  alienFileName += "/";
  alienFileName +=  localFileName;
  alienFileName += "?se=";
  alienFileName += fSE.Data();

  TString filename = 0x0;
  
  if(fStorage == 0) {
    filename = localFileName.Data();      
    AliInfo(Form("Writing merged tags to local file: %s",filename.Data()));
  } 
  if(fStorage == 1) {
    filename = alienFileName.Data();
    AliInfo(Form("Writing merged tags to grid file: %s",filename.Data()));     
  }
  
  TFile* ftag = TFile::Open(filename, "recreate");
  ftag->cd();
  ttag.Write();
  ftag->Close();

  delete tag;
  delete evTag;
  delete newTag;

  return kTRUE;
}

//_____________________________________________________________________________
void AliTagCreator::CreateTag(TFile* file, const char *guid, const char *md5, const char *turl, Long64_t size, Int_t Counter)
{
  /////////////
  //muon code//
  ////////////
  Double_t MUON_MASS = 0.105658369;
  //Variables
  Double_t X,Y,Z ;
  Double_t ThetaX, ThetaY, Pyz,Chisquare;
  Double_t PxRec,PyRec, PzRec, Energy;
  Int_t Charge;
  TLorentzVector EPvector;

  Float_t ZVertexCut = 10.0; 
  Float_t RhoVertexCut = 2.0; 

  Float_t LowPtCut = 1.0;
  Float_t HighPtCut = 3.0;
  Float_t VeryHighPtCut = 10.0;
  ////////////

  Double_t partFrac[5] = {0.01, 0.01, 0.85, 0.10, 0.05};

  // Creates the tags for all the events in a given ESD file
  Int_t ntrack;
  Int_t nProtons, nKaons, nPions, nMuons, nElectrons;
  Int_t nPos, nNeg, nNeutr;
  Int_t nK0s, nNeutrons, nPi0s, nGamas;
  Int_t nCh1GeV, nCh3GeV, nCh10GeV;
  Int_t nMu1GeV, nMu3GeV, nMu10GeV;
  Int_t nEl1GeV, nEl3GeV, nEl10GeV;
  Float_t maxPt = .0, meanPt = .0, totalP = .0;

  AliRunTag *tag = new AliRunTag();
  AliEventTag *evTag = new AliEventTag();
  TTree ttag("T","A Tree with event tags");
  TBranch * btag = ttag.Branch("AliTAG", &tag);
  btag->SetCompressionLevel(9);
  
  AliInfo(Form("Creating the tags......."));	
  
  Int_t firstEvent = 0,lastEvent = 0;
  TTree *t = (TTree*) file->Get("esdTree");
  TBranch * b = t->GetBranch("ESD");
  AliESD *esd = 0;
  b->SetAddress(&esd);
  
  tag->SetRunId(esd->GetRunNumber());
  
  Int_t iNumberOfEvents = b->GetEntries();
  for (Int_t iEventNumber = 0; iEventNumber < iNumberOfEvents; iEventNumber++) {
    ntrack = 0;
    nPos = 0;
    nNeg = 0;
    nNeutr =0;
    nK0s = 0;
    nNeutrons = 0;
    nPi0s = 0;
    nGamas = 0;
    nProtons = 0;
    nKaons = 0;
    nPions = 0;
    nMuons = 0;
    nElectrons = 0;	  
    nCh1GeV = 0;
    nCh3GeV = 0;
    nCh10GeV = 0;
    nMu1GeV = 0;
    nMu3GeV = 0;
    nMu10GeV = 0;
    nEl1GeV = 0;
    nEl3GeV = 0;
    nEl10GeV = 0;
    maxPt = .0;
    meanPt = .0;
    totalP = .0;
    
    b->GetEntry(iEventNumber);
    const AliESDVertex * vertexIn = esd->GetVertex();
    
    for (Int_t iTrackNumber = 0; iTrackNumber < esd->GetNumberOfTracks(); iTrackNumber++) {
      AliESDtrack * esdTrack = esd->GetTrack(iTrackNumber);
      UInt_t status = esdTrack->GetStatus();
      
      //select only tracks with ITS refit
      if ((status&AliESDtrack::kITSrefit)==0) continue;
      //select only tracks with TPC refit
      if ((status&AliESDtrack::kTPCrefit)==0) continue;
      
      //select only tracks with the "combined PID"
      if ((status&AliESDtrack::kESDpid)==0) continue;
      Double_t p[3];
      esdTrack->GetPxPyPz(p);
      Double_t momentum = sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2));
      Double_t fPt = sqrt(pow(p[0],2) + pow(p[1],2));
      totalP += momentum;
      meanPt += fPt;
      if(fPt > maxPt) maxPt = fPt;
      
      if(esdTrack->GetSign() > 0) {
	nPos++;
	if(fPt > LowPtCut) nCh1GeV++;
	if(fPt > HighPtCut) nCh3GeV++;
	if(fPt > VeryHighPtCut) nCh10GeV++;
      }
      if(esdTrack->GetSign() < 0) {
	nNeg++;
	if(fPt > LowPtCut) nCh1GeV++;
	if(fPt > HighPtCut) nCh3GeV++;
	if(fPt > VeryHighPtCut) nCh10GeV++;
      }
      if(esdTrack->GetSign() == 0) nNeutr++;
      
      //PID
      Double_t prob[5];
      esdTrack->GetESDpid(prob);
      
      Double_t rcc = 0.0;
      for(Int_t i = 0; i < AliPID::kSPECIES; i++) rcc += prob[i]*partFrac[i];
      if(rcc == 0.0) continue;
      //Bayes' formula
      Double_t w[5];
      for(Int_t i = 0; i < AliPID::kSPECIES; i++) w[i] = prob[i]*partFrac[i]/rcc;
      
      //protons
      if ((w[4]>w[3])&&(w[4]>w[2])&&(w[4]>w[1])&&(w[4]>w[0])) nProtons++;
      //kaons
      if ((w[3]>w[4])&&(w[3]>w[2])&&(w[3]>w[1])&&(w[3]>w[0])) nKaons++;
      //pions
      if ((w[2]>w[4])&&(w[2]>w[3])&&(w[2]>w[1])&&(w[2]>w[0])) nPions++; 
      //electrons
      if ((w[0]>w[4])&&(w[0]>w[3])&&(w[0]>w[2])&&(w[0]>w[1])) {
	nElectrons++;
	if(fPt > LowPtCut) nEl1GeV++;
	if(fPt > HighPtCut) nEl3GeV++;
	if(fPt > VeryHighPtCut) nEl10GeV++;
      }	  
      ntrack++;
    }//esd track loop
    
    /////////////
    //muon code//
    ////////////
    Int_t nMuonTracks = esd->GetNumberOfMuonTracks();
    // loop over all reconstructed tracks (also first track of combination)
    for (Int_t iTrack = 0; iTrack <  nMuonTracks;  iTrack++) {
      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);
      if (muonTrack == 0x0) continue;
      
      // Coordinates at vertex
      Z = muonTrack->GetZ(); 
      Y = muonTrack->GetBendingCoor();
      X = muonTrack->GetNonBendingCoor(); 
      
      ThetaX = muonTrack->GetThetaX();
      ThetaY = muonTrack->GetThetaY();
      
      Pyz = 1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
      PzRec = - Pyz / TMath::Sqrt(1.0 + TMath::Tan(ThetaY)*TMath::Tan(ThetaY));
      PxRec = PzRec * TMath::Tan(ThetaX);
      PyRec = PzRec * TMath::Tan(ThetaY);
      Charge = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));
      
      //ChiSquare of the track if needed
      Chisquare = muonTrack->GetChi2()/(2.0 * muonTrack->GetNHit() - 5);
      Energy = TMath::Sqrt(MUON_MASS * MUON_MASS + PxRec * PxRec + PyRec * PyRec + PzRec * PzRec);
      EPvector.SetPxPyPzE(PxRec, PyRec, PzRec, Energy);
      
      // total number of muons inside a vertex cut 
      if((TMath::Abs(Z)<ZVertexCut) && (TMath::Sqrt(Y*Y+X*X)<RhoVertexCut)) {
	nMuons++;
	if(EPvector.Pt() > LowPtCut) {
	  nMu1GeV++; 
	  if(EPvector.Pt() > HighPtCut) {
	    nMu3GeV++; 
	    if (EPvector.Pt() > VeryHighPtCut) {
	      nMu10GeV++;
	    }
	  }
	}
      }
    }//muon track loop
    
    // Fill the event tags 
    if(ntrack != 0)
      meanPt = meanPt/ntrack;
    
    evTag->SetEventId(iEventNumber+1);
    evTag->SetGUID(guid);
    evTag->SetMD5(md5);
    evTag->SetTURL(turl);
    evTag->SetSize(size);
    evTag->SetVertexX(vertexIn->GetXv());
    evTag->SetVertexY(vertexIn->GetYv());
    evTag->SetVertexZ(vertexIn->GetZv());
    
    evTag->SetT0VertexZ(esd->GetT0zVertex());
    
    evTag->SetTrigger(esd->GetTrigger());
    
    evTag->SetZDCNeutronEnergy(esd->GetZDCNEnergy());
    evTag->SetZDCProtonEnergy(esd->GetZDCPEnergy());
    evTag->SetZDCEMEnergy(esd->GetZDCEMEnergy());
    evTag->SetNumOfParticipants(esd->GetZDCParticipants());
    
    
    evTag->SetNumOfTracks(esd->GetNumberOfTracks());
    evTag->SetNumOfPosTracks(nPos);
    evTag->SetNumOfNegTracks(nNeg);
    evTag->SetNumOfNeutrTracks(nNeutr);
    
    evTag->SetNumOfV0s(esd->GetNumberOfV0s());
    evTag->SetNumOfCascades(esd->GetNumberOfCascades());
    evTag->SetNumOfKinks(esd->GetNumberOfKinks());
    evTag->SetNumOfPMDTracks(esd->GetNumberOfPmdTracks());
    
    evTag->SetNumOfProtons(nProtons);
    evTag->SetNumOfKaons(nKaons);
    evTag->SetNumOfPions(nPions);
    evTag->SetNumOfMuons(nMuons);
    evTag->SetNumOfElectrons(nElectrons);
    evTag->SetNumOfPhotons(nGamas);
    evTag->SetNumOfPi0s(nPi0s);
    evTag->SetNumOfNeutrons(nNeutrons);
    evTag->SetNumOfKaon0s(nK0s);
    
    evTag->SetNumOfChargedAbove1GeV(nCh1GeV);
    evTag->SetNumOfChargedAbove3GeV(nCh3GeV);
    evTag->SetNumOfChargedAbove10GeV(nCh10GeV);
    evTag->SetNumOfMuonsAbove1GeV(nMu1GeV);
    evTag->SetNumOfMuonsAbove3GeV(nMu3GeV);
    evTag->SetNumOfMuonsAbove10GeV(nMu10GeV);
    evTag->SetNumOfElectronsAbove1GeV(nEl1GeV);
    evTag->SetNumOfElectronsAbove3GeV(nEl3GeV);
    evTag->SetNumOfElectronsAbove10GeV(nEl10GeV);
    
    evTag->SetNumOfPHOSTracks(esd->GetNumberOfPHOSParticles());
    evTag->SetNumOfEMCALTracks(esd->GetNumberOfEMCALParticles());
    
    evTag->SetTotalMomentum(totalP);
    evTag->SetMeanPt(meanPt);
    evTag->SetMaxPt(maxPt);
    
    tag->AddEventTag(*evTag);
  }//event loop
  lastEvent = iNumberOfEvents;
  
  t->Delete("");
  
  ttag.Fill();
  tag->Clear();

  TString localFileName = "Run"; localFileName += tag->GetRunId(); 
  localFileName += ".Event"; localFileName += firstEvent; localFileName += "_"; localFileName += lastEvent; localFileName += "."; localFileName += Counter;
  localFileName += ".ESD.tag.root";

  TString alienLocation = "/alien";
  alienLocation += gGrid->Pwd();
  alienLocation += fgridpath.Data();
  alienLocation += "/";
  alienLocation +=  localFileName;
  alienLocation += "?se=";
  alienLocation += fSE.Data();

  TString fileName;
  
  if(fStorage == 0)
    {
      fileName = localFileName.Data();      
      AliInfo(Form("Writing tags to local file: %s",fileName.Data()));
    }
  if(fStorage == 1)
    {
      fileName = alienLocation.Data();
      AliInfo(Form("Writing tags to grid file: %s",fileName.Data()));
    }

  TFile* ftag = TFile::Open(fileName, "recreate");
  ftag->cd();
  ttag.Write();
  ftag->Close();

  delete ftag;
  delete esd;

  delete tag;
  delete evTag;
}

