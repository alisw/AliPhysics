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
#include <TObject.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <THashTable.h>

//ROOT-AliEn
#include <TGrid.h>
//#include <TAlienCollection.h>
#include <TGridResult.h>
#include <TFileMerger.h>
#include <TMap.h>
#include <TXMLParser.h>
#include <TThread.h>
#include <TTreePlayer.h>
#include <TProof.h>

//AliRoot
#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliLog.h"


#include "AliTagCreator.h"


ClassImp(AliTagCreator)


//______________________________________________________________________________
AliTagCreator::AliTagCreator() //local mode
{
  //==============Default constructor for a AliTagCreator==================
  fUser = "";
  fPasswd = "";  
  fSE = "";   
  fHost = "";
  fPort = 0; 
  fresult = 0;
}

//______________________________________________________________________________
AliTagCreator::AliTagCreator(const char *host, Int_t port, const char *username)
{
  //==============Default constructor for a AliTagCreator==================
  fresult = 0;
  fHost = host;
  fUser = username;
  fPort = port;

  TString grid = "alien://";
  grid += fHost;
  grid += ":";
  grid += port;
  
  //connect to the grid
  TGrid::Connect(grid.Data(),fUser.Data());
  if(!gGrid)
    {
      AliError("Connection failed!!!");
      return;
    }
}


//______________________________________________________________________________
AliTagCreator::AliTagCreator(const char *host, Int_t port, const char *username, const char *passwd)
{
  //==============Default constructor for a AliTagCreator==================
  fresult = 0;
  fHost = host;
  fUser = username;
  fPasswd = passwd;
  fPort = port;

  TString grid = "alien://";
  grid += fHost;
  grid += ":";
  grid += port;
  
  //connect to the grid
  TGrid::Connect(grid.Data(),fUser.Data(),fPasswd.Data());
   if(!gGrid)
    {
      AliError("Connection failed!!!");
      return;
    }
}


//______________________________________________________________________________
AliTagCreator::~AliTagCreator()
{
//================Default destructor for a AliTagCreator=======================
}

//______________________________________________________________________________
void AliTagCreator::SetSE(const char *se)
{
  fSE = se;
}

//______________________________________________________________________________
Bool_t AliTagCreator::ConnectToGrid(const char *host, Int_t port, const char *username)
{
  fHost = host;
  fUser = username;
  fPort = port;

  TString grid = "alien://";
  grid += fHost;
  grid += ":";
  grid += port;
  
  //connect to the grid
  TGrid::Connect(grid.Data(),fUser.Data());
  if(!gGrid)
    {
      AliError("Connection failed!!!");
      return kFALSE;
    } //connect to the grid

  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTagCreator::ReadESDCollection(TGridResult *res)
{
  fresult = res;
  Int_t NEntries = fresult->GetEntries();

  TString gridname = "alien://";
  TString AlienUrl;
  const char *guid;
 
  Int_t counter = 0;
  for(Int_t i = 0; i < NEntries; i++)
    {
      AlienUrl = gridname;
      AlienUrl += fresult->GetKey(i,"lfn");
      guid = fresult->GetKey(i,"guid");
      TFile *f = TFile::Open(AlienUrl,"READ");
      CreateTag(f,guid,counter);
      f->Close();
      delete f;	 
      counter += 1;
    }//grid result loop

  return kTRUE;
}


//_____________________________________________________________________________
void AliTagCreator::CreateTag(TFile* file, const char *guid, Int_t Counter)
{
  Int_t ntrack;
  Int_t NProtons, NKaons, NPions, NMuons, NElectrons;
  Int_t Npos, Nneg, Nneutr;
  Int_t NK0s, Nneutrons, Npi0s, Ngamas;
  Int_t Nch1GeV, Nch3GeV, Nch10GeV;
  Int_t Nmu1GeV, Nmu3GeV, Nmu10GeV;
  Int_t Nel1GeV, Nel3GeV, Nel10GeV;
  Float_t MaxPt = .0, MeanPt = .0, TotalP = .0;

  AliRunTag *tag = new AliRunTag();
  AliDetectorTag *detTag = new AliDetectorTag();
  AliEventTag *evTag = new AliEventTag();
  TTree ttag("T","A Tree with event tags");
  TBranch * btag = ttag.Branch("AliTAG", "AliRunTag", &tag);
  btag->SetCompressionLevel(9);
  
  AliInfo(Form("Creating the tags......."));	
  
  Int_t firstEvent = 0,lastEvent = 0;
  TTree *t = (TTree*) file->Get("esdTree");
  TBranch * b = t->GetBranch("ESD");
  AliESD *esd = 0;
  b->SetAddress(&esd);
  
  tag->SetRunId(esd->GetRunNumber());
  
  Int_t i_NumberOfEvents = b->GetEntries();
  for (Int_t i_EventNumber = 0; i_EventNumber < i_NumberOfEvents; i_EventNumber++)
    {
      ntrack = 0;
      Npos = 0;
      Nneg = 0;
      Nneutr =0;
      NK0s = 0;
      Nneutrons = 0;
      Npi0s = 0;
      Ngamas = 0;
      NProtons = 0;
      NKaons = 0;
      NPions = 0;
      NMuons = 0;
      NElectrons = 0;	  
      Nch1GeV = 0;
      Nch3GeV = 0;
      Nch10GeV = 0;
      Nmu1GeV = 0;
      Nmu3GeV = 0;
      Nmu10GeV = 0;
      Nel1GeV = 0;
      Nel3GeV = 0;
      Nel10GeV = 0;
      MaxPt = .0;
      MeanPt = .0;
      TotalP = .0;
      
      b->GetEntry(i_EventNumber);
      const AliESDVertex * VertexIn = esd->GetVertex();
      
      for (Int_t i_TrackNumber = 0; i_TrackNumber < esd->GetNumberOfTracks(); i_TrackNumber++)
	{
	  AliESDtrack * ESDTrack = esd->GetTrack(i_TrackNumber);
	  UInt_t status = ESDTrack->GetStatus();
	  
	  //select only tracks with ITS refit
	  if ((status&AliESDtrack::kITSrefit)==0) continue;
	  
	  //select only tracks with TPC refit-->remove extremely high Pt tracks
	  if ((status&AliESDtrack::kTPCrefit)==0) continue;
	  
	  //select only tracks with the "combined PID"
	  if ((status&AliESDtrack::kESDpid)==0) continue;
	  Double_t p[3];
	  ESDTrack->GetPxPyPz(p);
	  Double_t P = sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2));
	  Double_t fPt = sqrt(pow(p[0],2) + pow(p[1],2));
	  TotalP += P;
	  MeanPt += fPt;
	  if(fPt > MaxPt)
	    MaxPt = fPt;
	  
	  if(ESDTrack->GetSign() > 0)
	    {
	      Npos++;
	      if(fPt > 1.0)
		Nch1GeV++;
	      if(fPt > 3.0)
		Nch3GeV++;
	      if(fPt > 10.0)
		Nch10GeV++;
	    }
	  if(ESDTrack->GetSign() < 0)
	    {
	      Nneg++;
	      if(fPt > 1.0)
		Nch1GeV++;
	      if(fPt > 3.0)
		Nch3GeV++;
	      if(fPt > 10.0)
		Nch10GeV++;
	    }
	  if(ESDTrack->GetSign() == 0)
	    Nneutr++;
	  
	  //PID
	  Double_t prob[10];
	  ESDTrack->GetESDpid(prob);
	  
	  //K0s
	  if ((prob[8]>prob[7])&&(prob[8]>prob[6])&&(prob[8]>prob[5])&&(prob[8]>prob[4])&&(prob[8]>prob[3])&&(prob[8]>prob[2])&&(prob[8]>prob[1])&&(prob[8]>prob[0]))
	    NK0s++;
	  //neutrons
	  if ((prob[7]>prob[8])&&(prob[7]>prob[6])&&(prob[7]>prob[5])&&(prob[7]>prob[4])&&(prob[7]>prob[3])&&(prob[7]>prob[2])&&(prob[7]>prob[1])&&(prob[7]>prob[0]))
	    Nneutrons++; 
	  //pi0s
	  if ((prob[6]>prob[8])&&(prob[6]>prob[7])&&(prob[6]>prob[5])&&(prob[6]>prob[4])&&(prob[6]>prob[3])&&(prob[6]>prob[2])&&(prob[6]>prob[1])&&(prob[6]>prob[0]))
	    Npi0s++;
	  //gamas
	  if ((prob[5]>prob[8])&&(prob[5]>prob[7])&&(prob[5]>prob[6])&&(prob[5]>prob[4])&&(prob[5]>prob[3])&&(prob[5]>prob[2])&&(prob[5]>prob[1])&&(prob[5]>prob[0]))
	    Ngamas++;
	  //protons
	  if ((prob[4]>prob[8])&&(prob[4]>prob[7])&&(prob[4]>prob[6])&&(prob[4]>prob[5])&&(prob[4]>prob[3])&&(prob[4]>prob[2])&&(prob[4]>prob[1])&&(prob[4]>prob[0]))
	    NProtons++;
	  //kaons
	  if ((prob[3]>prob[8])&&(prob[3]>prob[7])&&(prob[3]>prob[6])&&(prob[3]>prob[5])&&(prob[3]>prob[4])&&(prob[3]>prob[2])&&(prob[3]>prob[1])&&(prob[3]>prob[0]))
	    NKaons++;
	  //pions
	  if ((prob[2]>prob[8])&&(prob[2]>prob[7])&&(prob[2]>prob[6])&&(prob[2]>prob[5])&&(prob[2]>prob[4])&&(prob[2]>prob[3])&&(prob[2]>prob[1])&&(prob[2]>prob[0]))
	    NPions++; 
	  //muons
	  if ((prob[1]>prob[8])&&(prob[1]>prob[7])&&(prob[1]>prob[6])&&(prob[1]>prob[5])&&(prob[1]>prob[4])&&(prob[1]>prob[3])&&(prob[1]>prob[2])&&(prob[1]>prob[0]))
	    {
	      NMuons++;
	      if(fPt > 1.0)
		Nmu1GeV++;
	      if(fPt > 3.0)
		Nmu3GeV++;
	      if(fPt > 10.0)
		Nmu10GeV++;
	    }
	  //electrons
	  if ((prob[0]>prob[8])&&(prob[0]>prob[7])&&(prob[0]>prob[6])&&(prob[0]>prob[5])&&(prob[0]>prob[4])&&(prob[0]>prob[3])&&(prob[0]>prob[2])&&(prob[0]>prob[1]))
	    {
	      NElectrons++;
	      if(fPt > 1.0)
		Nel1GeV++;
	      if(fPt > 3.0)
		Nel3GeV++;
	      if(fPt > 10.0)
		Nel10GeV++;
	    }
	  
	  
	  
	  ntrack++;
	}//track loop
      // Fill the event tags 
      if(ntrack != 0)
	MeanPt = MeanPt/ntrack;
      
      evTag->SetEventId(i_EventNumber+1);
      evTag->SetGUID(guid);
      evTag->SetVertexX(VertexIn->GetXv());
      evTag->SetVertexY(VertexIn->GetYv());
      evTag->SetVertexZ(VertexIn->GetZv());
      
      evTag->SetT0VertexZ(esd->GetT0zVertex());
      
      evTag->SetTrigger(esd->GetTrigger());
      
      evTag->SetZDCNeutronEnergy(esd->GetZDCNEnergy());
      evTag->SetZDCProtonEnergy(esd->GetZDCPEnergy());
      evTag->SetZDCEMEnergy(esd->GetZDCEMEnergy());
      evTag->SetNumOfParticipants(esd->GetZDCParticipants());
      
      
      evTag->SetNumOfTracks(esd->GetNumberOfTracks());
      evTag->SetNumOfPosTracks(Npos);
      evTag->SetNumOfNegTracks(Nneg);
      evTag->SetNumOfNeutrTracks(Nneutr);
      
      evTag->SetNumOfV0s(esd->GetNumberOfV0s());
      evTag->SetNumOfCascades(esd->GetNumberOfCascades());
      evTag->SetNumOfKinks(esd->GetNumberOfKinks());
      evTag->SetNumOfPMDTracks(esd->GetNumberOfPmdTracks());
      
      evTag->SetNumOfProtons(NProtons);
      evTag->SetNumOfKaons(NKaons);
      evTag->SetNumOfPions(NPions);
      evTag->SetNumOfMuons(NMuons);
      evTag->SetNumOfElectrons(NElectrons);
      evTag->SetNumOfPhotons(Ngamas);
      evTag->SetNumOfPi0s(Npi0s);
      evTag->SetNumOfNeutrons(Nneutrons);
      evTag->SetNumOfKaon0s(NK0s);
      
      evTag->SetNumOfChargedAbove1GeV(Nch1GeV);
      evTag->SetNumOfChargedAbove3GeV(Nch3GeV);
      evTag->SetNumOfChargedAbove10GeV(Nch10GeV);
      evTag->SetNumOfMuonsAbove1GeV(Nmu1GeV);
      evTag->SetNumOfMuonsAbove3GeV(Nmu3GeV);
      evTag->SetNumOfMuonsAbove10GeV(Nmu10GeV);
      evTag->SetNumOfElectronsAbove1GeV(Nel1GeV);
      evTag->SetNumOfElectronsAbove3GeV(Nel3GeV);
      evTag->SetNumOfElectronsAbove10GeV(Nel10GeV);
      
      evTag->SetNumOfPHOSTracks(esd->GetNumberOfPHOSParticles());
      evTag->SetNumOfEMCALTracks(esd->GetNumberOfEMCALParticles());
      
      evTag->SetTotalMomentum(TotalP);
      evTag->SetMeanPt(MeanPt);
      evTag->SetMaxPt(MaxPt);
  
      tag->AddEventTag(evTag);
    }
  lastEvent = i_NumberOfEvents;
  
  t->Delete("");
	
  ttag.Fill();
  tag->Clear();

  TString LocalfileName = "Run"; LocalfileName += tag->GetRunId(); 
  LocalfileName += ".Event"; LocalfileName += firstEvent; LocalfileName += "_"; LocalfileName += lastEvent; LocalfileName += "."; LocalfileName += Counter;
  LocalfileName += ".ESD.tag.root";

  cout<<"Writing tags to local file: "<<LocalfileName<<endl;

  TFile* ftag = TFile::Open(LocalfileName, "recreate");
  ftag->cd();
  ttag.Write();
  ftag->Close();

  delete ftag;
  delete esd;

  delete tag;
  delete detTag;
  delete evTag;
}

//_____________________________________________________________________________
Bool_t AliTagCreator::StoreGridTagFile(const char *localpath, const char *gridpath)
{
  gSystem->Load("libThread.so");
  gSystem->Load("libTreePlayer.so");
  gSystem->Load("libProof.so");
  cout<<"Storing tag file to alien's file catalog..."<<endl;
  TFileMerger merger;
  const char * pattern = "tag.root";
  // Open the working directory
  void * dirp = gSystem->OpenDirectory(localpath);
  const char * name = 0x0;
  TString AlienLocation;
  Char_t LocalLocation[256];
 
  // Add all files matching *pattern* to the chain
  while((name = gSystem->GetDirEntry(dirp)))
    {
      if (strstr(name,pattern))
	{
	  sprintf(LocalLocation,"file:%s/%s",localpath,name);
	  AlienLocation = "/alien";
	  AlienLocation += gGrid->Pwd();
	  AlienLocation += gridpath;
	  AlienLocation += "/";
	  AlienLocation += name;
	  AlienLocation += "?se=";
	  AlienLocation += fSE.Data();
	  merger.Cp(LocalLocation,AlienLocation);
	    
	}
    }	
  gSystem->FreeDirectory(dirp);
   
  return kTRUE;
}
