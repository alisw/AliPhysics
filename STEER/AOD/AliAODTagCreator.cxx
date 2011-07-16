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
//           AliAODTagCreator class
//   This is the class to deal with the tag creation (post process)
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

//ROOT
#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TRefArray.h>

//ROOT-AliEn
#include <TGrid.h>
#include <TGridResult.h>

//AliRoot
#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliFileTag.h"
#include "AliPID.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliLog.h"

#include "AliAODTagCreator.h"


ClassImp(AliAODTagCreator)


//______________________________________________________________________________
  AliAODTagCreator::AliAODTagCreator() :
    AliTagCreator(), 
    fChain(0), 
    fAODEvent(0), 
    fTreeT(0), 
    fRunTag(0),
    fTreeTEsd(0), 
    fRunTagEsd(0)  
{
  //==============Default constructor for a AliAODTagCreator================
}

//______________________________________________________________________________
AliAODTagCreator::~AliAODTagCreator() {
//================Default destructor for a AliAODTagCreator===================
  delete fChain;
  delete fAODEvent;
}

//______________________________________________________________________________
Bool_t AliAODTagCreator::ReadGridCollection(TGridResult *fresult) {
  // Reads the entry of the TGridResult and creates the tags
  Int_t nEntries = fresult->GetEntries();

  TString alienUrl;
  const char* guid;
  const char* md5;
  //  const char* turl;
  Long64_t size = -1;

  fChain = new TChain("aodTree");
 
  for(Int_t i = 0; i < nEntries; i++) {
    alienUrl = fresult->GetKey(i,"turl");
    guid = fresult->GetKey(i,"guid");
    if(fresult->GetKey(i,"size")) size = atol (fresult->GetKey(i,"size"));
    md5 = fresult->GetKey(i,"md5");
    //    turl = fresult->GetKey(i,"turl");
    if(md5 && !strlen(guid)) md5 = 0;
    if(guid && !strlen(guid)) guid = 0;
    
    fChain->Add(alienUrl);
  }//grid result loop
  
  AliInfo(Form("AOD chain created......."));	
  AliInfo(Form("Chain entries: %lld",fChain->GetEntries()));	

  CreateTag(fChain, "grid");
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODTagCreator::ReadLocalCollection(const char *localpath, const char* pattern) {
  // Checks the different subdirs of the given local path and in the
  // case where it finds an AliAOD.root file it creates the tags
  
  void *dira =  gSystem->OpenDirectory(localpath);
  Char_t fPath[512];
  const char * dirname  = 0x0;
  const char * filename = 0x0;

  fChain = new TChain("aodTree");

  while((dirname = gSystem->GetDirEntry(dira))) {
    snprintf(fPath,512,"%s/%s",localpath,dirname);
    void *dirb =  gSystem->OpenDirectory(fPath);
    while((filename = gSystem->GetDirEntry(dirb))) {
	TString bstr = dirname;
	if(bstr.Contains("..")) continue;
	if(strstr(filename,pattern)) {
	    TString aodFileName;
	    aodFileName = fPath;
	    aodFileName += "/";
	    aodFileName += pattern;
	    fChain->Add(aodFileName);
	} //pattern check
    } //child directory's entry loop
  } //parent directory's entry loop
  
  AliInfo(Form("AOD chain created......."));	
  AliInfo(Form("Chain entries: %lld",fChain->GetEntries()));	

  CreateTag(fChain, "local");

  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODTagCreator::ReadCAFCollection(const char *filename) {
  // Temporary solution for CAF: Takes as an input the ascii file that
  // lists the AODs stored in the SE of the CAF and creates the tags.

  // Open the input stream
  ifstream in;
  in.open(filename);

  TString aodfile;

  fChain = new TChain("aodTree");

  // Read the input list of files and add them to the chain
  while(in.good()) {
    in >> aodfile;
    if (!aodfile.Contains("root")) continue; // protection
    fChain->Add(aodfile);
  }

  AliInfo(Form("AOD chain created......."));	
  AliInfo(Form("Chain entries: %lld",fChain->GetEntries()));	

  CreateTag(fChain, "proof");

  return kTRUE;
}

//__________________________________________________________________________
void AliAODTagCreator::CreateAODTags(Int_t fFirstEvent, Int_t fLastEvent, TList */*grpList*/) {
  // Creates tag files for AODs

    AliInfo(Form("Creating the AOD tags......."));	
    
    TFile *file = TFile::Open("AliAOD.root");
    if (!file || !file->IsOpen()) {
	AliError(Form("opening failed"));
	delete file;
	return ;
    }

    fChain = new TChain("aodTree");
    fChain->Add("AliAOD.root");
    
    fAODEvent = new AliAODEvent();
    fAODEvent->ReadFromTree(fChain);
    
    Int_t lastEvent = 0;  
    if(fLastEvent == -1) lastEvent = (Int_t)fChain->GetEntries();
    else lastEvent = fLastEvent;

    char fileName[256];
    snprintf(fileName, 256, "Run%d.Event%d_%d.AOD.tag.root", 
	    fAODEvent->GetRunNumber(), fFirstEvent, lastEvent );
    AliInfo(Form("writing tags to file %s", fileName));
    AliDebug(1, Form("writing tags to file %s", fileName));
 
    TFile* ftag = TFile::Open(fileName, "recreate");

    fRunTag = new AliRunTag();
    fTreeT = new TTree("T","A Tree with event tags");
    TBranch * btag = fTreeT->Branch("AliTAG", &fRunTag);
    btag->SetCompressionLevel(9);

    CreateTags();

    ftag->cd();
    fTreeT->Fill();
    fRunTag->Clear();
    fTreeT->Write();
    ftag->Close();
    file->cd();
    file->Close();
}

//_____________________________________________________________________________
void AliAODTagCreator::CreateTag(TChain* chain, const char *type) {
    
    // Private method that creates tag files
    //

    //reading the esd tag file 
    fTreeTEsd = new TChain("T");
    const char * tagPattern = "ESD.tag";
    // Open the working directory
    void * dirp = gSystem->OpenDirectory(gSystem->pwd());
    const char * name = 0x0;
    // Add all files matching *pattern* to the chain
    while((name = gSystem->GetDirEntry(dirp))) {
	if (strstr(name,tagPattern)) fTreeTEsd->Add(name);
    }//directory loop
    AliInfo(Form("Chained tag files: %lld",fTreeTEsd->GetEntries()));
      
    fChain = chain;
    
    TString fSession = type;
    TString fguid, fmd5, fturl;
    fAODEvent = new AliAODEvent();
    fAODEvent->ReadFromTree(fChain);

    Int_t firstEvent = 0;
    
    TString localFileName = "Run"; 
    localFileName += fAODEvent->GetRunNumber(); 
    localFileName += ".Event"; 
    localFileName += firstEvent; 
    localFileName += "_"; 
    localFileName += chain->GetEntries(); //localFileName += "."; localFileName += Counter;
    localFileName += ".AOD.tag.root";
    
    TString fileName;
    
    if(fStorage == 0) {
	fileName = localFileName.Data();      
	AliInfo(Form("Writing tags to local file: %s",fileName.Data()));
    }
    else if(fStorage == 1) {
	TString alienLocation = "/alien";
	alienLocation += gGrid->Pwd();
	alienLocation += fgridpath.Data();
	alienLocation += "/";
	alienLocation +=  localFileName;
	alienLocation += "?se=";
	alienLocation += fSE.Data();
	fileName = alienLocation.Data();
	AliInfo(Form("Writing tags to grid file: %s",fileName.Data()));
    }
    
    TFile* ftag = TFile::Open(fileName, "recreate");
    
    fRunTag = new AliRunTag();
    fTreeT  = new TTree("T", "A Tree with event tags");
    TBranch * btag = fTreeT->Branch("AliTAG", &fRunTag);
    btag->SetCompressionLevel(9);
    
    // Access information from esd tag
    fRunTagEsd = new AliRunTag();
    fTreeTEsd->SetBranchAddress("AliTAG",&fRunTagEsd);
    // Creating new information of aod
    AliInfo(Form("Creating the AOD tags......."));
    CreateTags(type);
    ftag->cd();
    fRunTag->Clear();
    fTreeT->Write();
    ftag->Close();
}

void AliAODTagCreator::CreateTags(const char* /*type*/)
{
    // Event loop for tag creation
    TString fturl;
    TString fguid;
    Int_t   oldRun = -1;
    fChain->GetEntry(0);
    TFile *f = fChain->GetFile();
    TString ftempGuid = f->GetUUID().AsString();
    // Loop over events 
    Int_t nEvents = fChain->GetEntries();
    Int_t ntags    = 0;
    Int_t tagentry = 0;
    //    const TClonesArray *evTagList = 0;
    TString foldguid = "";

    for (Int_t iEventNumber = 0; iEventNumber < nEvents; iEventNumber++) {
	// Copy old tag information
	if (iEventNumber >= ntags) {
	    fTreeTEsd->GetEntry(tagentry++);
	    fRunTag->CopyStandardContent(fRunTagEsd);
// 	    evTagList = fRunTagEsd->GetEventTags();
// 	    ntags += evTagList->GetEntries();
	    ntags = fRunTagEsd->GetNEvents();
	}

	// Create a new Tag
	AliEventTag* evTag = new AliEventTag();
	// Read event
	fChain->GetEntry(iEventNumber);
	if (iEventNumber == 0) oldRun = fAODEvent->GetRunNumber();
	// Reference to the input file
	TFile *file = fChain->GetFile();
	//	const TUrl *url = file->GetEndpointUrl();
	fguid = file->GetUUID().AsString();

// 	if (!strcmp(type,"grid")) {
// 	    TString fturltemp = "alien://"; fturltemp += url->GetFile();
// 	    fturl = fturltemp(0,fturltemp.Index(".root",5,0,TString::kExact)+5);
// 	} else {
// 	    fturl = url->GetFile();
// 	}
	fturl = file->GetName();
	
	fAODEvent->GetStdContent();
	
	// Fill the event tag from the aod informatiom
	FillEventTag(fAODEvent, evTag);
	// Set the event and input file references
	//evTag->SetEventId(iEventNumber+1);
	
	// **** FIXME ****
// 	evTag->SetGUID(fguid);
// 	if(!strcmp(type,"grid")) {
// 	    evTag->SetMD5("");
// 	    evTag->SetTURL(fturl);
// 	    evTag->SetSize(0);
// 	    }
// 	else evTag->SetPath(fturl);
	//  **** FIXME ****

	// Check if a new run has to be created
	// File has changed
	if(fguid != ftempGuid) {
	    ftempGuid = fguid;
	    fTreeT->Fill();
	    fRunTag->Clear("");

	    AliFileTag *nftag = new AliFileTag();
	    
	    // if(fSession == "grid") {
	      nftag->SetMD5("");
	      nftag->SetTURL(fturl);
	      nftag->SetSize(0);
	    // }
	    // else {
	    //   nftag->SetPath(fturl);
	    //   nftag->SetSize(0);
	    //   nftag->SetMD5("");
	    //   nftag->SetTURL(fturl);
	    // }
      
	    if (fRunTag->GetFileId(fguid) > -1)
	      AliFatal("Adding a file which is already in the RunTag.");
	    
	    fRunTag->AddFileTag(nftag);
	    
	}

	// Run# has changed
	if (oldRun != (fAODEvent->GetRunNumber()))
	{
	    oldRun = fAODEvent->GetRunNumber();

	    fTreeT->Fill();
	    fRunTag->Clear("");
	    ftempGuid = fguid;
	    fTreeT->Fill();
	    fRunTag->Clear("");

	    AliFileTag *nftag = new AliFileTag();
	    
	    // if(fSession == "grid") {
	      nftag->SetMD5("");
	      nftag->SetTURL(fturl);
	      nftag->SetSize(0);
	    // }
	    // else {
	    //   nftag->SetPath(fturl);
	    //   nftag->SetSize(0);
	    //   nftag->SetMD5("");
	    //   nftag->SetTURL(fturl);
	    // }
      
	    if (fRunTag->GetFileId(fguid) > -1)
	      AliFatal("Adding a file which is already in the RunTag.");
	    
	    fRunTag->AddFileTag(nftag);
	    
	}
	
	// Add the event tag
	fRunTag->AddEventTag(*evTag);
	delete evTag;
	// Last event
	if(iEventNumber+1 == fChain->GetEntries()) {
	    fTreeT->Fill();
	    fRunTag->Clear("");
	}      
    }//event loop
}



void AliAODTagCreator::FillEventTag(AliAODEvent* aod, AliEventTag* evTag)
{
// 
// Fill the event tag information	
    //
    fAODEvent = aod;
    
    //
    Float_t fLowPtCut      =  1.0;
    Float_t fHighPtCut     =  3.0;
    Float_t fVeryHighPtCut = 10.0;
    ////////////                            
    Double_t partFrac[10] = {0.01, 0.01, 0.85, 0.10, 0.05, 0., 0., 0., 0., 0.};
    
    // Creates the tags for all the events in a given AOD file
    Int_t   ntrack = 0;
    Int_t   nPos = 0, nNeg = 0, nNeutr =0;
    Int_t   nKinks = 0, nV0s = 0, nCascades = 0;
    Int_t   nK0s = 0, nNeutrons = 0, nPi0s = 0, nGamas = 0;
    Int_t   nProtons = 0,  nKaons = 0, nPions = 0, nMuons = 0, nElectrons = 0, nFWMuons = 0;
    Int_t   nCh1GeV = 0, nCh3GeV = 0, nCh10GeV = 0;
    Int_t   nMu1GeV = 0, nMu3GeV = 0, nMu10GeV = 0;
    Int_t   nEl1GeV = 0, nEl3GeV = 0, nEl10GeV = 0;
    Float_t maxPt =  .0, etamaxPt = -999., phimaxPt = -999., meanPt = .0, totalP =  .0;

    TRefArray tmp;


    // Primary Vertex
    AliAODVertex *pVertex = fAODEvent->GetPrimaryVertex();
    if (pVertex) {
	evTag->SetVertexX(pVertex->GetX());
	evTag->SetVertexY(pVertex->GetY());
	evTag->SetVertexZ(pVertex->GetZ());
	Double_t covmatrix[6];
	pVertex->GetCovarianceMatrix(covmatrix);
	evTag->SetVertexZError(sqrt(covmatrix[5]));
    }
    // loop over vertices 
    Int_t nVtxs = fAODEvent->GetNumberOfVertices();
    for (Int_t nVtx = 0; nVtx < nVtxs; nVtx++) {
	AliAODVertex *vertex = fAODEvent->GetVertex(nVtx);
	if(vertex->GetType() == 1) nKinks    += 1;
	if(vertex->GetType() == 2) nV0s      += 1;
	if(vertex->GetType() == 3) nCascades += 1;
    }
    Int_t nTracks = fAODEvent->GetNTracks();
    for (Int_t nTr = 0; nTr < nTracks; nTr++) {
	AliAODTrack *track = fAODEvent->GetTrack(nTr);
	
	Double_t fPt = track->Pt();
	if(fPt > maxPt) {
	    maxPt = fPt;
	    etamaxPt = track->Eta();
	    phimaxPt = track->Phi();
	}
	
	if(track->Charge() > 0) {
	    nPos++;
	    if(fPt > fLowPtCut)      nCh1GeV++;
	    if(fPt > fHighPtCut)     nCh3GeV++;
	    if(fPt > fVeryHighPtCut) nCh10GeV++;
	}
	if(track->Charge() < 0) {
	    nNeg++;
	    if(fPt > fLowPtCut)      nCh1GeV++;
	    if(fPt > fHighPtCut)     nCh3GeV++;
	    if(fPt > fVeryHighPtCut) nCh10GeV++;
	    }
	if(track->Charge() == 0) nNeutr++;
	//PID
	const Double32_t *prob = track->PID();
	Double_t rcc = 0.0;
	for(Int_t i = 0; i < AliPID::kSPECIES; i++) rcc += prob[i]*partFrac[i];
	if(rcc == 0.0) continue;
	//Bayes' formula  
	Double_t w[10];
	for(Int_t i = 0; i < AliPID::kSPECIES; i++) w[i] = prob[i]*partFrac[i]/rcc;
	
	//protons
	if ((w[4]>w[3])&&(w[4]>w[2])&&(w[4]>w[1])&&(w[4]>w[0])) nProtons++;
	//kaons
	if ((w[3]>w[4])&&(w[3]>w[2])&&(w[3]>w[1])&&(w[3]>w[0])) nKaons++;
	//pions
	if ((w[2]>w[4])&&(w[2]>w[3])&&(w[2]>w[1])&&(w[2]>w[0])) nPions++;
	//muons 
	if ((w[1]>w[4])&&(w[1]>w[3])&&(w[1]>w[2])&&(w[1]>w[0])) {
	    nMuons++;
	    if(fPt > fLowPtCut)      nMu1GeV++;
	    if(fPt > fHighPtCut)     nMu3GeV++;
	    if(fPt > fVeryHighPtCut) nMu10GeV++;
	}
	//electrons  
	if ((w[0]>w[4])&&(w[0]>w[3])&&(w[0]>w[2])&&(w[0]>w[1])) {
	    nElectrons++;
	    if(fPt > fLowPtCut) nEl1GeV++;
	    if(fPt > fHighPtCut) nEl3GeV++;
	    if(fPt > fVeryHighPtCut) nEl10GeV++;
	}
	totalP += track->P();
	meanPt += fPt;
	ntrack++;
	// forward muons (in the dimuon spectrometer)
	if(track->IsMuonTrack()) nFWMuons++;   
  	                          
    }//track loop
    //
    // Fill the event tags  
    if(ntrack != 0)
	meanPt = meanPt/ntrack;
    
    
    evTag->SetNumOfTracks(nTracks);
    evTag->SetNumOfPosTracks(nPos);
    evTag->SetNumOfNegTracks(nNeg);
    evTag->SetNumOfNeutrTracks(nNeutr);
    
    evTag->SetNumOfV0s(nV0s);
    evTag->SetNumOfCascades(nCascades);
    evTag->SetNumOfKinks(nKinks);
    
    evTag->SetNumOfProtons(nProtons);
    evTag->SetNumOfKaons(nKaons);
    evTag->SetNumOfPions(nPions);
    evTag->SetNumOfMuons(nMuons);
    evTag->SetNumOfFWMuons(nFWMuons);
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

    tmp.Clear();
    evTag->SetNumOfPHOSClusters(fAODEvent->GetPHOSClusters(&tmp));
    tmp.Clear();
    evTag->SetNumOfEMCALClusters(fAODEvent->GetEMCALClusters(&tmp));
	
    evTag->SetTotalMomentum(totalP);
    evTag->SetMeanPt(meanPt);
    evTag->SetMaxPt(maxPt);
    evTag->SetEtaMaxPt(etamaxPt);
    evTag->SetPhiMaxPt(phimaxPt);
}
