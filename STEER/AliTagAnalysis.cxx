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
//           AliTagAnalysis class
//   This is the class to deal with the tag analysis
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

//ROOT
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TEventList.h>
#include <TEntryList.h>

//ROOT-AliEn
#include <TGridResult.h>

#include "AliLog.h"

#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliTagAnalysis.h"
#include "AliEventTagCuts.h"
#include "AliRunTagCuts.h"
#include "AliXMLCollection.h"

class TTree;

ClassImp(AliTagAnalysis)

TChain *AliTagAnalysis::fgChain = 0;

//___________________________________________________________________________
AliTagAnalysis::AliTagAnalysis(): 
  TObject(),
  ftagresult(0x0),
  fTagDirName(),
  fChain(0x0)
{
  //Default constructor for a AliTagAnalysis
}

//___________________________________________________________________________
AliTagAnalysis::~AliTagAnalysis() {
//Default destructor for a AliTagAnalysis
}

//___________________________________________________________________________
void AliTagAnalysis::ChainLocalTags(const char *dirname) {
  //Searches the entries of the provided direcory
  //Chains the tags that are stored locally
  fTagDirName = dirname;
  TString fTagFilename;
  
  TChain *fgChain = new TChain("T");
  fChain = fgChain;
  
  const char * tagPattern = "tag.root";
  // Open the working directory
  void * dirp = gSystem->OpenDirectory(fTagDirName);
  const char * name = 0x0;
  // Add all files matching *pattern* to the chain
  while((name = gSystem->GetDirEntry(dirp))) {
    if (strstr(name,tagPattern)) { 
      fTagFilename = fTagDirName;
      fTagFilename += "/";
      fTagFilename += name;
	  	
      fChain->Add(fTagFilename);  
    }//pattern check
  }//directory loop
  AliInfo(Form("Chained tag files: %d ",fChain->GetEntries()));
}


//___________________________________________________________________________
void AliTagAnalysis::ChainGridTags(TGridResult *res) {
  //Loops overs the entries of the TGridResult
  //Chains the tags that are stored in the GRID
  ftagresult = res;
  Int_t nEntries = ftagresult->GetEntries();
 
  TChain *fgChain = new TChain("T");
  fChain = fgChain;

  TString gridname = "alien://";
  TString alienUrl;
 
  for(Int_t i = 0; i < nEntries; i++) {
    alienUrl = ftagresult->GetKey(i,"turl");
    fChain->Add(alienUrl);
  }//grid result loop  
}


//___________________________________________________________________________
TChain *AliTagAnalysis::QueryTags(AliRunTagCuts *RunTagCuts, AliEventTagCuts *EvTagCuts) {
  //Queries the tag chain using the defined 
  //event tag cuts from the AliEventTagCuts object
  //and returns a TChain along with the associated TEventList
  AliInfo(Form("Querying the tags........"));
  
  //ESD file chain
  TChain *fESDchain = new TChain("esdTree");
  //Event list
  TEventList *fEventList = new TEventList();
  
  //Defining tag objects
  AliRunTag *tag = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  TString guid = 0;
  TString turl = 0;
  TString path = 0;

  Int_t iAccepted = 0;
  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    fChain->GetEntry(iTagFiles);
    if(RunTagCuts->IsAccepted(tag)) {
      Int_t iEvents = tag->GetNEvents();
      const TClonesArray *tagList = tag->GetEventTags();
      for(Int_t i = 0; i < iEvents; i++) {
	evTag = (AliEventTag *) tagList->At(i);
	guid = evTag->GetGUID(); 
	turl = evTag->GetTURL(); 
	path = evTag->GetPath();
	if(EvTagCuts->IsAccepted(evTag)) fEventList->Enter(iAccepted+i);
      }//event loop
      iAccepted += iEvents;
    
      if(path != "") fESDchain->AddFile(path);
      else if(turl != "") fESDchain->AddFile(turl);
    }//run tags cut
  }//tag file loop
  AliInfo(Form("Accepted events: %d",fEventList->GetN()));
  fESDchain->SetEventList(fEventList);
   
  return fESDchain;
}

//___________________________________________________________________________
Bool_t AliTagAnalysis::CreateXMLCollection(const char* name, AliRunTagCuts *RunTagCuts, AliEventTagCuts *EvTagCuts) {
  //Queries the tag chain using the defined 
  //event tag cuts from the AliEventTagCuts object
  //and returns a XML collection
  AliInfo(Form("Creating the collection........"));

  AliXMLCollection *collection = new AliXMLCollection();
  collection->SetCollectionName(name);
  collection->WriteHeader();

  //Event list
  TEventList *fEventList = new TEventList();
  TString guid = 0x0;
  TString turl = 0x0;
  TString lfn = 0x0;
  
  //Defining tag objects
  AliRunTag *tag = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    fChain->GetEntry(iTagFiles);
    if(RunTagCuts->IsAccepted(tag)) {
      Int_t iEvents = tag->GetNEvents();
      const TClonesArray *tagList = tag->GetEventTags();
      for(Int_t i = 0; i < iEvents; i++) {
	evTag = (AliEventTag *) tagList->At(i);
	guid = evTag->GetGUID(); 
	turl = evTag->GetTURL(); 
	lfn = turl(8,turl.Length());

	if(EvTagCuts->IsAccepted(evTag)) fEventList->Enter(i);
      }//event loop
      collection->WriteBody(iTagFiles+1,guid,lfn,turl,fEventList);
      fEventList->Clear();
    }//run tag cuts
  }//tag file loop
  collection->Export();

  return kTRUE;
}

//___________________________________________________________________________
TChain *AliTagAnalysis::GetInputChain(const char* system, const char*global, const char *wn) {
  //returns the chain+event list - used in batch sessions
  TString fsystem = system;
  Int_t iAccepted = 0;
  //ESD file chain
  TChain *fESDchain = new TChain("esdTree");
  //Event list
  TEventList *fEventList = new TEventList();
  AliXMLCollection *collection1 = AliXMLCollection::Open(global);
  AliXMLCollection *collection2 = AliXMLCollection::Open(wn);
  collection1->OverlapCollection(collection2);
  collection1->Reset();
  while (collection1->Next()) {
    AliInfo(Form("Adding: %s",collection1->GetTURL("")));
    fESDchain->Add(collection1->GetTURL(""));
    TEntryList *list = (TEntryList *)collection1->GetEventList("");
    for(Int_t i = 0; i < list->GetN(); i++) fEventList->Enter(iAccepted+list->GetEntry(i));

    if(fsystem == "pp") iAccepted += 100;
    else if(fsystem == "PbPb") iAccepted += 1;
  }

  fESDchain->SetEventList(fEventList);
  
  AliInfo(Form("Number of selected events: %d",fEventList->GetN()));

  return fESDchain;
}
