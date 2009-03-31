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
#include <Riostream.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TEventList.h>
#include <TEntryList.h>
#include <TTreeFormula.h>

//ROOT-AliEn
#include <TGridResult.h>

#include "AliLog.h"

#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliTagAnalysis.h"
#include "AliEventTagCuts.h"
#include "AliDetectorTagCuts.h"
#include "AliLHCTagCuts.h"
#include "AliRunTagCuts.h"
#include "AliXMLCollection.h"

class TTree;

ClassImp(AliTagAnalysis)

//___________________________________________________________________________
AliTagAnalysis::AliTagAnalysis(): 
  TObject(),
  ftagresult(0x0),
  fTagDirName(),
  fChain(0x0),
  fAnalysisType(),
  fGlobalList(0) {
  //Default constructor for a AliTagAnalysis
}

//___________________________________________________________________________
AliTagAnalysis::AliTagAnalysis(const char* type): 
  TObject(),
  ftagresult(0x0),
  fTagDirName(),
  fChain(0x0),
  fAnalysisType(type),
  fGlobalList(0) {
  //constructor for a AliTagAnalysis
}

//___________________________________________________________________________
AliTagAnalysis::~AliTagAnalysis() {
  //Default destructor for a AliTagAnalysis
  if(ftagresult) delete ftagresult;
  if(fChain) delete fChain;
  if(fGlobalList) delete fGlobalList;
}

//___________________________________________________________________________
Bool_t  AliTagAnalysis::AddTagsFile(const char *alienUrl) {
  // Add a single tags file to the chain

  Bool_t rv = kTRUE ;

  if (! fChain) fChain = new TChain("T");

  TFile *f = TFile::Open(alienUrl,"READ");
  fChain->Add(alienUrl);
  AliInfo(Form("Chained tag files: %d ",fChain->GetEntries()));
  delete f;

  if (fChain->GetEntries() == 0 )
    rv = kFALSE ;

  return rv ;
}

//___________________________________________________________________________
void AliTagAnalysis::ChainLocalTags(const char *dirname) {
  //Searches the entries of the provided direcory
  //Chains the tags that are stored locally
  fTagDirName = dirname;
  TString fTagFilename;
  
  if (! fChain)  fChain = new TChain("T");
  const char * tagPattern = 0x0;
  if(fAnalysisType == "ESD") tagPattern = "ESD.tag.root";
  else if(fAnalysisType == "AOD") tagPattern = "AOD.tag.root";
  else AliFatal("Only ESD and AOD type is implemented!!!");

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
TChain * AliTagAnalysis::ChainGridTags(TGridResult *res) {
  //Loops overs the entries of the TGridResult
  //Chains the tags that are stored in the GRID
  ftagresult = res;
  Int_t nEntries = ftagresult->GetEntries();
 
  if (! fChain)  fChain = new TChain("T");

  TString gridname = "alien://";
  TString alienUrl;
 
  for(Int_t i = 0; i < nEntries; i++) {
    alienUrl = ftagresult->GetKey(i,"turl");
    fChain->Add(alienUrl);
  }//grid result loop
  return fChain;
}


//___________________________________________________________________________
TChain *AliTagAnalysis::QueryTags(AliRunTagCuts *runTagCuts, 
				  AliLHCTagCuts *lhcTagCuts, 
				  AliDetectorTagCuts *detTagCuts, 
				  AliEventTagCuts *evTagCuts) {
  //Queries the tag chain using the defined 
  //event tag cuts from the AliEventTagCuts object
  //and returns a TChain along with the associated TEventList
  AliInfo(Form("Querying the tags........"));

  TString fAliceFile;
  if(fAnalysisType == "ESD") fAliceFile = "esdTree";
  else if(fAnalysisType == "AOD") fAliceFile = "aodTree";
  else AliFatal("Only ESD and AOD type is implemented!!!");

  //ESD file chain
  TChain *fESDchain = new TChain(fAliceFile.Data());
  //global entry list
  fGlobalList = new TEntryList();
  
  //Defining tag objects
  AliRunTag   *tag     = new AliRunTag;
  AliEventTag *evTag   = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  TString guid = 0;
  TString turl = 0;
  TString path = 0;

  Int_t iAccepted = 0;
  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    fChain->GetEntry(iTagFiles);
    if(runTagCuts->IsAccepted(tag)) {
      if(lhcTagCuts->IsAccepted(tag->GetLHCTag())) {
	if(detTagCuts->IsAccepted(tag->GetDetectorTags())) {
	  TEntryList *fLocalList = new TEntryList();
	  Int_t iEvents = tag->GetNEvents();
	  const TClonesArray *tagList = tag->GetEventTags();
	  for(Int_t i = 0; i < iEvents; i++) {
	    evTag = (AliEventTag *) tagList->At(i);
	    guid = evTag->GetGUID(); 
	    turl = evTag->GetTURL(); 
	    path = evTag->GetPath();
	    fLocalList->SetTreeName(fAliceFile.Data());
	    if(turl!="") fLocalList->SetFileName(turl.Data());
	    else fLocalList->SetFileName(path.Data());
	    if(evTagCuts->IsAccepted(evTag)) fLocalList->Enter(i);
	  }//event loop
	  if(path != "") fESDchain->AddFile(path);
	  else if(turl != "") fESDchain->AddFile(turl);
	  fGlobalList->Add(fLocalList);
	  iAccepted += fLocalList->GetN();
	}//detector tag cuts
      }//lhc tag cuts
    }//run tags cut
    tag->Clear();
  }//tag file loop
  AliInfo(Form("Accepted events: %d",iAccepted));
  fESDchain->SetEntryList(fGlobalList,"ne");
   
  return fESDchain;
}

//___________________________________________________________________________
TChain *AliTagAnalysis::QueryTags(const char *fRunCut, 
				  const char *fLHCCut, 
				  const char *fDetectorCut, 
				  const char *fEventCut) { 	 
  //Queries the tag chain using the defined 	 
  //event tag cuts from the AliEventTagCuts object 	 
  //and returns a TChain along with the associated TEventList 	 
  AliInfo(Form("Querying the tags........")); 	 
  
  TString fAliceFile;
  if(fAnalysisType == "ESD") fAliceFile = "esdTree";
  else if(fAnalysisType == "AOD") fAliceFile = "aodTree";
  else AliFatal("Only ESD and AOD type is implemented!!!");

  //ESD file chain
  TChain *fESDchain = new TChain(fAliceFile.Data());
  //global entry list
  fGlobalList = new TEntryList();
  
  //Defining tag objects 	 
  AliRunTag *tag = new AliRunTag; 	 
  AliEventTag *evTag = new AliEventTag; 	 
  fChain->SetBranchAddress("AliTAG",&tag); 	 
  
  TString guid = 0; 	 
  TString turl = 0; 	 
  TString path = 0; 	 
  
  TTreeFormula *fRunFormula = new TTreeFormula("fRun",fRunCut,fChain); 	 
  TTreeFormula *fLHCFormula = new TTreeFormula("fLHC",fLHCCut,fChain); 	 
  TTreeFormula *fDetectorFormula = new TTreeFormula("fDetector",fDetectorCut,fChain);
  TTreeFormula *fEventFormula = new TTreeFormula("fEvent",fEventCut,fChain);
  
  Int_t current = -1; 	 
  Int_t iAccepted = 0; 	 
  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    fChain->GetEntry(iTagFiles); 	 
    if (current != fChain->GetTreeNumber()) { 	 
      fRunFormula->UpdateFormulaLeaves(); 	 
      fLHCFormula->UpdateFormulaLeaves(); 	 
      fDetectorFormula->UpdateFormulaLeaves(); 	 
      fEventFormula->UpdateFormulaLeaves(); 	 
      current = fChain->GetTreeNumber(); 	 
    } 	 
    if(fRunFormula->EvalInstance(iTagFiles) == 1) { 	 
      if(fLHCFormula->EvalInstance(iTagFiles) == 1) { 	 
	if(fDetectorFormula->EvalInstance(iTagFiles) == 1) { 	 
	  TEntryList *fLocalList = new TEntryList();
	  Int_t iEvents = fEventFormula->GetNdata(); 	 
	  const TClonesArray *tagList = tag->GetEventTags(); 	 
	  for(Int_t i = 0; i < iEvents; i++) { 	 
	    evTag = (AliEventTag *) tagList->At(i); 	 
	    guid = evTag->GetGUID(); 	 
	    turl = evTag->GetTURL(); 	 
	    path = evTag->GetPath(); 	 
	    fLocalList->SetTreeName(fAliceFile.Data());
	    fLocalList->SetFileName(turl.Data());
	    if(fEventFormula->EvalInstance(i) == 1) fLocalList->Enter(i);
	  }//event loop 	 
	  iAccepted += fLocalList->GetN(); 	 
	  
	  if(path != "") fESDchain->AddFile(path); 	 
	  else if(turl != "") fESDchain->AddFile(turl); 	 
	  fGlobalList->Add(fLocalList);
	  iAccepted += fLocalList->GetN();
	}//detector tag cuts
      }//lhc tag cuts
    }//run tag cut 	 
  }//tag file loop 	 
  AliInfo(Form("Accepted events: %d",iAccepted)); 	 
  fESDchain->SetEntryList(fGlobalList,"ne"); 	 
  
  return fESDchain; 	 
}

//___________________________________________________________________________
Bool_t AliTagAnalysis::CreateXMLCollection(const char* name, 
					   AliRunTagCuts *runTagCuts, 
					   AliLHCTagCuts *lhcTagCuts, 
					   AliDetectorTagCuts *detTagCuts, 
					   AliEventTagCuts *evTagCuts) {
  //Queries the tag chain using the defined 
  //event tag cuts from the AliEventTagCuts object
  //and returns a XML collection
  AliInfo(Form("Creating the collection........"));

  AliXMLCollection *collection = new AliXMLCollection();
  collection->SetCollectionName(name);
  collection->WriteHeader();

  TString guid = 0x0;
  TString turl = 0x0;
  TString lfn = 0x0;
  
  //Defining tag objects
  AliRunTag *tag = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    //Event list
    TEntryList *fList = new TEntryList();
    fChain->GetEntry(iTagFiles);
    if(runTagCuts->IsAccepted(tag)) {
      if(lhcTagCuts->IsAccepted(tag->GetLHCTag())) {
	if(detTagCuts->IsAccepted(tag->GetDetectorTags())) {
	  Int_t iEvents = tag->GetNEvents();
	  const TClonesArray *tagList = tag->GetEventTags();
	  for(Int_t i = 0; i < iEvents; i++) {
	    evTag = (AliEventTag *) tagList->At(i);
	    guid = evTag->GetGUID(); 
	    turl = evTag->GetTURL(); 
	    lfn = turl(8,turl.Length());
	    if(evTagCuts->IsAccepted(evTag)) fList->Enter(i);
	  }//event loop
	  collection->WriteBody(iTagFiles+1,guid,lfn,turl,fList);
	}//detector tag cuts
      }//lhc tag cuts 
    }//run tag cuts
    tag->Clear();
  }//tag file loop
  collection->Export();

  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliTagAnalysis::CreateXMLCollection(const char* name, 
					   const char *fRunCut, 
					   const char *fLHCCut, 
					   const char *fDetectorCut, 
					   const char *fEventCut) {
  //Queries the tag chain using the defined 
  //event tag cuts from the AliEventTagCuts object
  //and returns a XML collection
  AliInfo(Form("Creating the collection........"));

  AliXMLCollection *collection = new AliXMLCollection();
  collection->SetCollectionName(name);
  collection->WriteHeader();

  TString guid = 0x0;
  TString turl = 0x0;
  TString lfn = 0x0;
  
  //Defining tag objects
  AliRunTag *tag = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  TTreeFormula *fRunFormula = new TTreeFormula("fRun",fRunCut,fChain);
  TTreeFormula *fLHCFormula = new TTreeFormula("fLHC",fLHCCut,fChain); 	 
  TTreeFormula *fDetectorFormula = new TTreeFormula("fDetector",fDetectorCut,fChain);
  TTreeFormula *fEventFormula = new TTreeFormula("fEvent",fEventCut,fChain);

  Int_t current = -1;
  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    //Event list
    TEntryList *fList = new TEntryList();
    fChain->GetEntry(iTagFiles);
    if (current != fChain->GetTreeNumber()) {
      fRunFormula->UpdateFormulaLeaves();
      fLHCFormula->UpdateFormulaLeaves();
      fDetectorFormula->UpdateFormulaLeaves();
      fEventFormula->UpdateFormulaLeaves();
      current = fChain->GetTreeNumber();
    }
    if(fRunFormula->EvalInstance(iTagFiles) == 1) {
      if(fLHCFormula->EvalInstance(iTagFiles) == 1) { 	 
	if(fDetectorFormula->EvalInstance(iTagFiles) == 1) { 	 
	  Int_t iEvents = fEventFormula->GetNdata();
	  const TClonesArray *tagList = tag->GetEventTags();
	  for(Int_t i = 0; i < iEvents; i++) {
	    evTag = (AliEventTag *) tagList->At(i);
	    guid = evTag->GetGUID(); 
	    turl = evTag->GetTURL(); 
	    lfn = turl(8,turl.Length());
	    if(fEventFormula->EvalInstance(i) == 1) fList->Enter(i);
	  }//event loop
	  collection->WriteBody(iTagFiles+1,guid,lfn,turl,fList);
	}//detector tag cuts
      }//lhc tag cuts 
    }//run tag cuts
  }//tag file loop
  collection->Export();

  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliTagAnalysis::CreateAsciiCollection(const char* name, 
					     AliRunTagCuts *runTagCuts, 
					     AliLHCTagCuts *lhcTagCuts, 
					     AliDetectorTagCuts *detTagCuts, 
					     AliEventTagCuts *evTagCuts) {
  //Queries the tag chain using the defined 
  //event tag cuts from the AliEventTagCuts object
  //and returns a XML collection
  AliInfo(Form("Creating the collection........"));

  ofstream fout;
  fout.open(name);

  TString guid = 0x0;
  TString turl = 0x0;
  TString lfn = 0x0;

  TString line0 = 0;

  //Defining tag objects
  AliRunTag *tag = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    //Event list
    TEntryList *fList = new TEntryList();
    fChain->GetEntry(iTagFiles);
    if(runTagCuts->IsAccepted(tag)) {
      if(lhcTagCuts->IsAccepted(tag->GetLHCTag())) {
	if(detTagCuts->IsAccepted(tag->GetDetectorTags())) {
	  Int_t iEvents = tag->GetNEvents();
	  const TClonesArray *tagList = tag->GetEventTags();
	  for(Int_t i = 0; i < iEvents; i++) {
	    evTag = (AliEventTag *) tagList->At(i);
	    guid = evTag->GetGUID(); 
	    turl = evTag->GetTURL(); 
	    lfn = turl(8,turl.Length());
	    if(evTagCuts->IsAccepted(evTag)) fList->Enter(i);
	  }//event loop
	  line0 = guid; line0 += " "; line0 += turl; line0 += " ";
	  for(Int_t i = 0; i < fList->GetN(); i++) {
	    line0 += fList->GetEntry(i); 
	    line0 += " ";
	  }  
	  fout<<line0<<"\n";
	}//detector tag cuts
      }//lhc tag cuts 
    }//run tag cuts
    tag->Clear();
  }//tag file loop

  fout.close();

  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliTagAnalysis::CreateAsciiCollection(const char* name, 
					     const char *fRunCut, 
					     const char *fLHCCut, 
					     const char *fDetectorCut, 
					     const char *fEventCut) {
  //Queries the tag chain using the defined 
  //event tag cuts from the AliEventTagCuts object
  //and returns a XML collection
  AliInfo(Form("Creating the collection........"));

  ofstream fout;
  fout.open(name);

  TString guid = 0x0;
  TString turl = 0x0;
  TString lfn = 0x0;

  TString line0 = 0;
  
  //Defining tag objects
  AliRunTag *tag = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  TTreeFormula *fRunFormula = new TTreeFormula("fRun",fRunCut,fChain);
  TTreeFormula *fLHCFormula = new TTreeFormula("fLHC",fLHCCut,fChain); 	 
  TTreeFormula *fDetectorFormula = new TTreeFormula("fDetector",fDetectorCut,fChain);
  TTreeFormula *fEventFormula = new TTreeFormula("fEvent",fEventCut,fChain);

  Int_t current = -1;
  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    //Event list
    TEntryList *fList = new TEntryList();
    fChain->GetEntry(iTagFiles);
    if (current != fChain->GetTreeNumber()) {
      fRunFormula->UpdateFormulaLeaves();
      fLHCFormula->UpdateFormulaLeaves();
      fDetectorFormula->UpdateFormulaLeaves();
      fEventFormula->UpdateFormulaLeaves();
      current = fChain->GetTreeNumber();
    }
    if(fRunFormula->EvalInstance(iTagFiles) == 1) {
      if(fLHCFormula->EvalInstance(iTagFiles) == 1) { 	 
	if(fDetectorFormula->EvalInstance(iTagFiles) == 1) { 	 
	  Int_t iEvents = fEventFormula->GetNdata();
	  const TClonesArray *tagList = tag->GetEventTags();
	  for(Int_t i = 0; i < iEvents; i++) {
	    evTag = (AliEventTag *) tagList->At(i);
	    guid = evTag->GetGUID(); 
	    turl = evTag->GetTURL(); 
	    lfn = turl(8,turl.Length());
	    if(fEventFormula->EvalInstance(i) == 1) fList->Enter(i);
	  }//event loop
	  line0 = guid; line0 += " "; line0 += turl; line0 += " ";
	  for(Int_t i = 0; i < fList->GetN(); i++) {
	    line0 += fList->GetEntry(i); 
	    line0 += " ";
	  }  
	  fout<<line0<<"\n";
	}//detector tag cuts
      }//lhc tag cuts 
    }//run tag cuts
  }//tag file loop

  fout.close();

  return kTRUE;
}

//___________________________________________________________________________
TChain *AliTagAnalysis::GetInputChain(const char* system, const char *wn) {
  //returns the chain+event list - used in batch sessions
  // this function will be removed once the new root 
  // improvements are committed
  TString fsystem = system;
  Int_t iAccepted = 0;

  TChain *fAnalysisChain = 0;
  if(fAnalysisType == "ESD") fAnalysisChain = new TChain("esdTree");
  else if(fAnalysisType == "AOD") fAnalysisChain = new TChain("aodTree");
  else AliFatal("Only ESD and AOD type is implemented!!!");
  
  //Event list
  TEventList *fEventList = new TEventList();
  AliXMLCollection *collection = AliXMLCollection::Open(wn);

  collection->Reset();
  while (collection->Next()) {
    AliInfo(Form("Adding: %s",collection->GetTURL("")));
    fAnalysisChain->Add(collection->GetTURL(""));
    TEntryList *list = (TEntryList *)collection->GetEventList("");
    for(Int_t i = 0; i < list->GetN(); i++) fEventList->Enter(iAccepted+list->GetEntry(i));

    if(fsystem == "pp") iAccepted += 100;
    else if(fsystem == "PbPb") iAccepted += 1;
  }

  fAnalysisChain->SetEventList(fEventList);
  
  AliInfo(Form("Number of selected events: %d",fEventList->GetN()));

  return fAnalysisChain;
}

//___________________________________________________________________________
TChain *AliTagAnalysis::GetChainFromCollection(const char* collectionname, 
					       const char* treename) {
  //returns the TChain+TEntryList object- used in batch sessions
  TString fAliceFile = treename;
  Int_t iAccepted = 0;
  TChain *fAnalysisChain = 0;
  if(fAliceFile == "esdTree") fAnalysisChain = new TChain("esdTree");
  else if(fAliceFile == "aodTree") fAnalysisChain = new TChain("aodTree");
  else AliFatal("Inconsistent tree name - use esdTree or aodTree!");

  //Event list
  fGlobalList = new TEntryList();
  AliXMLCollection *collection = AliXMLCollection::Open(collectionname);

  collection->Reset();
  while (collection->Next()) {
    AliInfo(Form("Adding: %s",collection->GetTURL("")));
    fAnalysisChain->Add(collection->GetTURL(""));
    TEntryList *list = (TEntryList *)collection->GetEventList("");
    list->SetTreeName(fAliceFile.Data());
    list->SetFileName(collection->GetTURL(""));
    fGlobalList->Add(list);
    iAccepted += list->GetN();
  }

  fAnalysisChain->SetEntryList(fGlobalList,"ne");
  
  AliInfo(Form("Number of selected events: %d",iAccepted));

  return fAnalysisChain;
}
