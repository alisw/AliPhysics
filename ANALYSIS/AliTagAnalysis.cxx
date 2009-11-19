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
#include <TMap.h>

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
      printf("Tag file %s\n", fTagFilename.Data());
      
    }//pattern check
  }//directory loop
  //AliInfo(Form("Chained tag files: %d ",fChain->GetEntries()));
 // AliDebug(Form("Chained tag files: %d ",fChain->GetEntries()));
  fChain->ls();
  
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

  Bool_t aod = kFALSE;
  TString aliceFile;
  if(fAnalysisType == "ESD") aliceFile = "esdTree";
  else if(fAnalysisType == "AOD") {
      aliceFile = "aodTree";
      aod = kTRUE;
  }
  else AliFatal("Only ESD and AOD type is implemented!!!");

  //ESD file chain
  TChain *esdChain = new TChain(aliceFile.Data());
  //global entry list
  fGlobalList = new TEntryList();
  
  //Defining tag objects
  AliRunTag   *tag     = new AliRunTag;
  AliEventTag *evTag   = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  TString guid;
  TString turl;
  TString path;

  TTree*      cTree     = 0; 
  TEntryList* localList = new TEntryList();

  Int_t iAccepted = 0;
  Int_t iev       = 0;
  Int_t ientry    = 0;
  Int_t cEntries  = 0;
  
  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    fChain->GetEntry(iTagFiles);
    TTree* tree = fChain->GetTree();
    if (cTree != tree) {
	// Fix for aod tags: for each tree in the chain, merge the entries
	cTree    = tree;
	cEntries = tree->GetEntries();
	iev      = 0;
	ientry   = 0;
    }

    if(runTagCuts->IsAccepted(tag)) {
      if(lhcTagCuts->IsAccepted(tag->GetLHCTag())) {
	if(detTagCuts->IsAccepted(tag->GetDetectorTags())) {
	  if ((iev == 0) || !aod) localList->Reset();
	  Int_t iEvents = tag->GetNEvents();
	  const TClonesArray *tagList = tag->GetEventTags();
	  for(Int_t i = 0; i < iEvents; i++) {
	    evTag = (AliEventTag *) tagList->At(i);
	    guid = evTag->GetGUID(); 
	    turl = evTag->GetTURL(); 
	    path = evTag->GetPath();
	    localList->SetTreeName(aliceFile.Data());
	    if(turl!="") localList->SetFileName(turl.Data());
	    else localList->SetFileName(path.Data());

	    if(evTagCuts->IsAccepted(evTag)) {
		if(aod) localList->Enter(iev);
		else localList->Enter(i);
	    }
	    iev++;
	  }//event loop
	  if ((ientry == cEntries-1) || !aod) {
	      iAccepted += localList->GetN();
	      if(path != "") esdChain->AddFile(path);
	      else if(turl != "") esdChain->AddFile(turl);
	      fGlobalList->Add(localList);
	  }
	}//detector tag cuts
      }//lhc tag cuts
    }//run tags cut
    tag->Clear();
    ientry++;
  }//tag file loop
  AliInfo(Form("Accepted events: %d",iAccepted));
  esdChain->ls();
  esdChain->SetEntryList(fGlobalList,"ne");
  delete tag;
  delete localList;
  
  return esdChain;
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
  
  Bool_t aod = kFALSE;
  TString aliceFile;
  if(fAnalysisType == "ESD") aliceFile = "esdTree";
  else if(fAnalysisType == "AOD") {
      aliceFile = "aodTree";
      aod = kTRUE;
  }
  else AliFatal("Only ESD and AOD type is implemented!!!");

  //ESD file chain
  TChain *esdChain = new TChain(aliceFile.Data());
  //global entry list
  fGlobalList = new TEntryList();
  
  //Defining tag objects 	 
  AliRunTag   *tag   = new AliRunTag; 	 
  AliEventTag *evTag = new AliEventTag; 	 
  fChain->SetBranchAddress("AliTAG",&tag); 	 
  
  TString guid; 	 
  TString turl; 	 
  TString path; 	 
  
  TTreeFormula *fRunFormula = new TTreeFormula("fRun",fRunCut,fChain); 	 
  TTreeFormula *fLHCFormula = new TTreeFormula("fLHC",fLHCCut,fChain); 	 
  TTreeFormula *fDetectorFormula = new TTreeFormula("fDetector",fDetectorCut,fChain);
  TTreeFormula *fEventFormula = new TTreeFormula("fEvent",fEventCut,fChain);
  
  TEntryList* localList = new TEntryList();

  Int_t iev       = 0;
  Int_t ientry    = 0;
  Int_t cEntries  = 0;
  Int_t current   = -1; 	 
  Int_t iAccepted = 0; 	 

  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {
    fChain->GetEntry(iTagFiles); 	 
    if (current != fChain->GetTreeNumber()) { 	 
      fRunFormula->UpdateFormulaLeaves(); 	 
      fLHCFormula->UpdateFormulaLeaves(); 	 
      fDetectorFormula->UpdateFormulaLeaves(); 	 
      fEventFormula->UpdateFormulaLeaves(); 	 
      // Fix for aod tags: for each tree in the chain, merge the entries
      cEntries = (fChain->GetTree())->GetEntries();
      iev      = 0;
      ientry   = 0;
      //
      current = fChain->GetTreeNumber(); 	 
    } 	 
    if(fRunFormula->EvalInstance(iTagFiles) == 1) { 	 
      if(fLHCFormula->EvalInstance(iTagFiles) == 1) { 	 
	if(fDetectorFormula->EvalInstance(iTagFiles) == 1) {
          if ((iev == 0) || !aod) localList->Reset(); 	 
	  Int_t iEvents = fEventFormula->GetNdata(); 	 
	  const TClonesArray *tagList = tag->GetEventTags(); 	 
	  for(Int_t i = 0; i < iEvents; i++) { 	 
	    evTag = (AliEventTag *) tagList->At(i); 	 
	    guid = evTag->GetGUID(); 	 
	    turl = evTag->GetTURL(); 	 
	    path = evTag->GetPath(); 	 
	    localList->SetTreeName(aliceFile.Data());
	    localList->SetFileName(turl.Data());
	    if(fEventFormula->EvalInstance(i) == 1) {
		if(aod) localList->Enter(iev);
		else localList->Enter(i);
	    }
	    iev++;
	  }//event loop 	 

	  if ((ientry == cEntries-1) || !aod) {  
	      if(path != "") esdChain->AddFile(path); 	 
	      else if(turl != "") esdChain->AddFile(turl); 	 
	      fGlobalList->Add(localList);
	      iAccepted += localList->GetN();
	  }
	}//detector tag cuts
      }//lhc tag cuts
    }//run tag cut 	 
    tag->Clear();
    ientry++;
  }//tag file loop 	 
  AliInfo(Form("Accepted events: %d",iAccepted)); 	 
  esdChain->SetEntryList(fGlobalList,"ne"); 	 

  delete tag;
  delete localList;
  return esdChain; 	 
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

  Bool_t aod = kFALSE;
  if(fAnalysisType == "AOD") aod = kTRUE;


  AliXMLCollection *collection = new AliXMLCollection();
  collection->SetCollectionName(name);
  collection->WriteHeader();

  TString guid;
  TString turl;
  TString lfn;
  
  TTree*      cTree = 0; 
  TEntryList* localList = new TEntryList();
  Int_t iAccepted = 0;
  Int_t iev       = 0;
  Int_t ientry    = 0;
  Int_t cEntries  = 0;

  Int_t iRejectedRun = 0;
  Int_t iRejectedLHC = 0;
  Int_t iRejectedDet = 0;
  Int_t iRejectedEvt = 0;

  Int_t iTotalEvents = 0;

  Int_t iAcceptedEvtInFile = 0;
  Int_t iRejectedEvtInFile = 0;

  //Defining tag objects
  AliRunTag   *tag   = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {

    fChain->GetEntry(iTagFiles);
    TTree* tree = fChain->GetTree();
    if (cTree != tree) {
	// Fix for aod tags: for each tree in the chain, merge the entries
	cTree    = tree;
	cEntries = tree->GetEntries();
	iev      = 0;
	ientry   = 0;
    }
    //Event list
    iTotalEvents += tag->GetNEvents();
    if ((iev == 0) || !aod) localList->Reset();
    if(runTagCuts->IsAccepted(tag)) {
      if(lhcTagCuts->IsAccepted(tag->GetLHCTag())) {
	if(detTagCuts->IsAccepted(tag->GetDetectorTags())) {
	  Int_t iEvents = tag->GetNEvents();
	  const TClonesArray *tagList = tag->GetEventTags();
	  iRejectedEvtInFile = 0;
	  iAcceptedEvtInFile = 0;
	  for(Int_t i = 0; i < iEvents; i++) {
	    evTag = (AliEventTag *) tagList->At(i);
	    guid = evTag->GetGUID(); 
	    turl = evTag->GetTURL(); 
	    lfn = turl(8,turl.Length());
	    if(evTagCuts->IsAccepted(evTag)) {
		if(aod) localList->Enter(iev);
		else localList->Enter(i);
		iAcceptedEvtInFile++;
	    }
	    else {
	      iRejectedEvt++;
	      iRejectedEvtInFile++;
	    }
	    iev++;
	  }//event loop
	  if ((ientry == cEntries-1) || !aod) {
	    iAccepted += localList->GetN();
	    collection->WriteBody(iTagFiles+1,guid,lfn,turl,localList,iAcceptedEvtInFile,iRejectedEvtInFile);
	  }
	}//detector tag cuts
	else {
	  iRejectedDet += tag->GetNEvents();
	}
      }//lhc tag cuts 
      else {
	iRejectedLHC += tag->GetNEvents();
      }
    }//run tag cuts
    else {
      iRejectedRun += tag->GetNEvents();
    }
    tag->Clear();
    ientry++;
  }//tag file loop
  collection->WriteSummary(iTotalEvents, iAccepted, iRejectedRun, iRejectedLHC, iRejectedDet, iRejectedEvt);
  collection->Export();
  
  delete tag;
  delete localList;
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

  Bool_t aod = kFALSE;
  if(fAnalysisType == "AOD") aod = kTRUE;

  AliXMLCollection *collection = new AliXMLCollection();
  collection->SetCollectionName(name);
  collection->WriteHeader();

  TString guid;
  TString turl;
  TString lfn;
  TEntryList* localList = new TEntryList();
  
  Int_t iAccepted = 0;
  Int_t iev       = 0;
  Int_t ientry    = 0;
  Int_t cEntries  = 0;

  Int_t iRejectedRun = 0;
  Int_t iRejectedLHC = 0;
  Int_t iRejectedDet = 0;
  Int_t iRejectedEvt = 0;

  Int_t iTotalEvents = 0;

  Int_t iAcceptedEvtInFile = 0;
  Int_t iRejectedEvtInFile = 0;

  //Defining tag objects
  AliRunTag *tag     = new AliRunTag;
  AliEventTag *evTag = new AliEventTag;
  fChain->SetBranchAddress("AliTAG",&tag);

  TTreeFormula *fRunFormula = new TTreeFormula("fRun",fRunCut,fChain);
  TTreeFormula *fLHCFormula = new TTreeFormula("fLHC",fLHCCut,fChain); 	 
  TTreeFormula *fDetectorFormula = new TTreeFormula("fDetector",fDetectorCut,fChain);
  TTreeFormula *fEventFormula = new TTreeFormula("fEvent",fEventCut,fChain);

  Int_t current = -1;
  for(Int_t iTagFiles = 0; iTagFiles < fChain->GetEntries(); iTagFiles++) {

    fChain->GetEntry(iTagFiles);
    if (current != fChain->GetTreeNumber()) {
      fRunFormula->UpdateFormulaLeaves();
      fLHCFormula->UpdateFormulaLeaves();
      fDetectorFormula->UpdateFormulaLeaves();
      fEventFormula->UpdateFormulaLeaves();
      // Fix for aod tags: for each tree in the chain, merge the entries
      cEntries = (fChain->GetTree())->GetEntries();
      iev      = 0;
      ientry   = 0;
      //
      current = fChain->GetTreeNumber();
    }
    //Event list
    iTotalEvents += tag->GetNEvents();
    if ((iev == 0) || !aod) localList->Reset();
    if(fRunFormula->EvalInstance(iTagFiles) == 1) {
      if(fLHCFormula->EvalInstance(iTagFiles) == 1) { 	 
	if(fDetectorFormula->EvalInstance(iTagFiles) == 1) { 	 
	  Int_t iEvents = fEventFormula->GetNdata();
	  const TClonesArray *tagList = tag->GetEventTags();
	  iRejectedEvtInFile = 0;
	  iAcceptedEvtInFile = 0;
	  for(Int_t i = 0; i < iEvents; i++) {
	    evTag = (AliEventTag *) tagList->At(i);
	    guid = evTag->GetGUID(); 
	    turl = evTag->GetTURL(); 
	    lfn = turl(8,turl.Length());
	    if(fEventFormula->EvalInstance(i) == 1) {
		if(aod) localList->Enter(iev);
		else localList->Enter(i);
		iAcceptedEvtInFile++;
	    }
	    else {
	      iRejectedEvt++;
	      iRejectedEvtInFile++;
	    }
	    iev++;
	  }//event loop
	  if ((ientry == cEntries-1) || !aod) {
	    collection->WriteBody(iTagFiles+1,guid,lfn,turl,localList,iAcceptedEvtInFile, iRejectedEvtInFile);
	    iAccepted += localList->GetN();
	  }
	}//detector tag cuts
	else {
	  iRejectedDet += tag->GetNEvents();
	}
      }//lhc tag cuts 
      else {
	iRejectedLHC += tag->GetNEvents();
      }
    }//run tag cuts
    else {
      iRejectedRun += tag->GetNEvents();
    }
    ientry++;
  }//tag file loop
  collection->WriteSummary(iTotalEvents, iAccepted, iRejectedRun, iRejectedLHC, iRejectedDet, iRejectedEvt);
  collection->Export();

  delete tag;
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
  TString aliceFile = treename;
  Int_t iAccepted = 0;
  TChain *fAnalysisChain = 0;
  if(aliceFile == "esdTree") fAnalysisChain = new TChain("esdTree");
  else if(aliceFile == "aodTree") fAnalysisChain = new TChain("aodTree");
  else AliFatal("Inconsistent tree name - use esdTree or aodTree!");

  //Event list
  fGlobalList = new TEntryList();
  AliXMLCollection *collection = AliXMLCollection::Open(collectionname);

  // Tag selection summary per file
  TMap *tagCutSummary = new TMap();
  tagCutSummary->SetName("TagCutSumm");

  collection->Reset();
  while (collection->Next()) {
    AliInfo(Form("Adding: %s",collection->GetTURL("")));
    fAnalysisChain->Add(collection->GetTURL(""));
    TEntryList *list = (TEntryList *)collection->GetEventList("");
    list->SetTreeName(aliceFile.Data());
    list->SetFileName(collection->GetTURL(""));
    fGlobalList->Add(list);
    iAccepted += list->GetN();
    if (collection->GetCutSumm())
      tagCutSummary->Add(new TObjString(collection->GetTURL("")), new TObjString(collection->GetCutSumm()));
  }

  fAnalysisChain->SetEntryList(fGlobalList,"ne");
  
  AliInfo(Form("Number of selected events: %d",iAccepted));

  TList *aUserInfo = fAnalysisChain->GetUserInfo();
  aUserInfo->Add(tagCutSummary);

  Int_t iAccEv;
  Int_t iTotalEvents;
  Int_t iRejRun;
  Int_t iRejLHC;
  Int_t iRejDet;
  Int_t iRejEvt;

  collection->GetCollectionSummary(&iTotalEvents, &iAccEv, &iRejRun, &iRejLHC, &iRejDet, &iRejEvt);
 
  char nstr[2000];

  sprintf(nstr, "TotalEvents=%i", iTotalEvents);
  TObjString *iTotStr = new TObjString(nstr);
  aUserInfo->Add(iTotStr);

  sprintf(nstr, "AcceptedEvents=%i", iAccepted);
  TObjString *iAccStr = new TObjString(nstr);
  aUserInfo->Add(iAccStr);

  sprintf(nstr, "RejectedRun=%i", iRejRun);
  TObjString *iRejRunStr = new TObjString(nstr);
  aUserInfo->Add(iRejRunStr);

  sprintf(nstr, "RejectedLHC=%i", iRejLHC);
  TObjString *iRejLHCStr = new TObjString(nstr);
  aUserInfo->Add(iRejLHCStr);

  sprintf(nstr, "RejectedDet=%i", iRejDet);
  TObjString *iRejDetStr = new TObjString(nstr);
  aUserInfo->Add(iRejDetStr);

  sprintf(nstr, "RejectedEvt=%i", iRejEvt);
  TObjString *iRejEvtStr = new TObjString(nstr);
  aUserInfo->Add(iRejEvtStr);

  return fAnalysisChain;
}
