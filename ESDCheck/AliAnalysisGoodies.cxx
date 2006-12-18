/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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
//_________________________________________________________________________
// Various utilities usefull for analysis
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisDataContainer.h" 
#include "AliTagAnalysis.h" 
#include "AliEventTagCuts.h" 
#include "AliRunTagCuts.h" 
#include "AliXMLCollection.h" 
#include "AliAnalysisGoodies.h" 
#include "AliAnalysisManager.h" 
#include "AliAnalysisTask.h" 
  //#include "AliPhotonAnalysisTask.h" 
#include "AliLog.h" 

#include <Riostream.h>
#include <TAlienCollection.h>
#include <TChain.h>
#include <TFileMerger.h>
#include <TGrid.h>
#include <TROOT.h> 
#include <TSystem.h>


//______________________________________________________________________________
AliAnalysisGoodies::AliAnalysisGoodies() :
  fESDTreeName("esdTree"), 
  fnumberOfTasks(0),
  fTaskList(0),
  fTaskInType(0), 
  fTaskOuType(0)
{
  fTimer.Reset() ; 
   
  TString token = gSystem->Getenv("GRID_TOKEN") ; 
  
  if ( token == "OK" ) 
    TGrid::Connect("alien://");
  else 
    AliInfo("You are not connected to the GRID") ; 
}

//______________________________________________________________________________
void AliAnalysisGoodies::Help() const  
{
  AliInfo("Analysis utilities:\n") ; 
  printf("                ***  Copy  : copy files ESD files listed in an xml collection from AliEn catalog to local storage and creates a local xml collection  \n") ; 
  printf("                                        usage: Copy(in, out)\n") ; 
  printf("                                                in: a xml esd collection file name    \n") ;  
  printf("                                                ou: the local directory where to save the esd root files   \n") ;  
  printf("                ***  Make  : makes esd collection from tags  \n") ; 
  printf("                                        usage: Make(tags, esds)\n") ; 
  printf("                                                tags: is either a tag root file or an xml tag collection   \n") ;  
  printf("                                                esds: is an esd collection     \n") ;  
  printf("                ***  Merge  : merges files listed in a xml collection \n") ; 
  printf("                                        usage Merge(collection, outputDile)\n") ; 
  printf("                                               collection: is a xml collection \n") ;  
  printf("                ***  Process : process the events with an Analysis Task \n") ;
  printf("                                        usage: Process(esdFile, tagCuts) \n") ;
  printf("                                                esdFile: can be a root file with the ESD Tree ( ex: esd?AliESDs.root) \n") ;
  printf("                                                        or a root file with the Tag Tree ( ex: tag?Run100.Event0_100.ESD.tag.root) \n") ;
  printf("                                                        or a local or alien xml file with the ESD collection ( ex: esd?esdCollection.xml) \n") ;
  printf("                                                        or a local or alien xml file with the TAG collection ( ex: tag?tagCollection.xml) \n") ;
  printf("                                                        or a TChain of esd TTrees \n") ;
  printf("                                               tagCuts: is the AliEventTagCuts (needed only for tag? cases \n") ;
  printf("                ***  Register: register files already stored in a MSS into the AliEn catalog\n") ;
  printf("                                        usage: Register(lfndir, pfndir, pfnFileName) \n") ; 
  printf("                                                lfndir : AliEn directory ( ex:  /alice/data/2006/LHC06c/PHOS_TestBeam/\n") ;  
  printf("                                                pfndir : MSS directory   ( ex: /castor/cern.ch/alice/testbeam/phos/2006 \n") ;
  printf("                                                file   : text file with a list of the file names to be registered\n ") ; 

}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::Alien2Local(const TString collectionNameIn, const TString localDir)
{
  // copy files ESD files listed in an xml collection from AliEn catalog to local storage and creates a local xml collection
  // usage: Alien2Local(in, out)
  //        in: a xml esd collection file name 
  //        ou: the local directory where to save the esd root files          

  Bool_t rv = kTRUE ; 

  fTimer.Start() ; 

  AliXMLCollection * collectionIn = AliXMLCollection::Open(collectionNameIn) ;
  collectionIn->Reset() ; 

  AliXMLCollection * collectionOu = new AliXMLCollection() ; 
  TString collectionNameOu(collectionIn->GetCollectionName()) ; 
  collectionNameOu.Append("Local") ; 
  collectionOu->SetCollectionName(collectionNameOu) ; 
  collectionOu->WriteHeader() ; 

  TFileMerger merger ; 
  
  const char* ocwd = gSystem->WorkingDirectory();

  Int_t counter = 0 ;  
  while ( collectionIn->Next() ) {
    gSystem->ChangeDirectory(localDir) ; 
    TString fileTURL = collectionIn->GetTURL("") ; 

    TString tempo(fileTURL) ; 
    tempo.Remove(tempo.Last('/'), tempo.Length()) ; 
    TString evtsNumber = tempo(tempo.Last('/')+1, tempo.Length())+"/";
    tempo.Remove(tempo.Last('/'), tempo.Length()) ; 
    TString runNumber = tempo(tempo.Last('/')+1, tempo.Length())+"/" ; 
    TString dir = localDir + runNumber ; 
    gSystem->MakeDirectory(dir) ; 
    gSystem->ChangeDirectory(dir) ; 
    dir += evtsNumber + "/"; 
    gSystem->MakeDirectory(dir) ; 
    gSystem->ChangeDirectory(dir) ; 
    dir += collectionIn->GetCollectionName() ; 
    TEntryList * list = collectionIn->GetEventList("") ; 
    
    collectionOu->WriteBody(counter, collectionIn->GetGUID(""), collectionIn->GetLFN(""), collectionIn->GetTURL(""), list) ;
    counter++ ; 
    printf("Copying %s to %s\n", fileTURL.Data(), dir.Data()) ;  
    merger.Cp(fileTURL, dir) ;
  }
  collectionOu->Export() ;
  gSystem->ChangeDirectory(ocwd) ; 
  
  fTimer.Stop();
  fTimer.Print();

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::Make(AliRunTagCuts *runCuts, AliEventTagCuts *evtCuts, const char * in, const char * out) const  
{
  // makes esd collection from tags 
  // usage Make(tags, esds)
  //              tags: is either a tag root file or an xml tag collection  
  //              esds: is an esd collection  

  Bool_t rv = kTRUE ; 

  if ( !evtCuts && !runCuts ) {
    AliError("No Tag cuts provided") ; 
    return kFALSE ; 
  }
 
  TString file(in) ; 
  if ( file.Contains(".root") ) 
    rv = MakeEsdCollectionFromTagFile(runCuts, evtCuts, file.Data(), out) ; 
  else  if ( file.Contains(".xml") ) 
    rv = MakeEsdCollectionFromTagCollection(runCuts, evtCuts, file.Data(), out) ;
  else {
    AliError(Form("%s is not a valid file format", in)) ; 
    rv = kFALSE ; 
  }

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::MakeEsdCollectionFromTagFile(AliRunTagCuts *runCuts, AliEventTagCuts *evtCuts, const char * in, const char * out) const 
{
  // Makes an esd collection from a root tag file 
  Bool_t rv = kTRUE ; 
    // Open the file collection 
  printf("*** Create Collection       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  file   = |%s|             \n",in);              	
 
  AliTagAnalysis * tagAna = new AliTagAnalysis(); 
  rv = tagAna->AddTagsFile(in);
  if ( ! rv ) 
    return rv ; 
 
  tagAna->CreateXMLCollection(out, runCuts, evtCuts) ;
 
  return rv ; 

}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::MakeEsdCollectionFromTagCollection(AliRunTagCuts * runCuts, AliEventTagCuts * evtCuts, const char * in, const char * out) const 
{
  // Makes an esd collection from a xml tag collection 
  Bool_t rv = kTRUE ; 
   // Open the file collection 
  printf("*** Create Collection       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",in);              	
  
  TAlienCollection * collection = TAlienCollection::Open(in);
  TGridResult* result = collection->GetGridResult("");
  AliTagAnalysis * tagAna = new AliTagAnalysis(); 
  tagAna->ChainGridTags(result);

  tagAna->CreateXMLCollection(out, runCuts, evtCuts) ;

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::MakeEsdCollectionFromTagCollection(const char * runCuts, const char * evtCuts, const char * in, const char * out) const 
{
  // Makes an esd collection from a xml tag collection 
  
  Bool_t rv = kTRUE ; 
 
  // Open the file collection 
  printf("*** Create Collection       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",in);              	
  
  TAlienCollection * collection = TAlienCollection::Open(in);
  TGridResult* result = collection->GetGridResult("");
  AliTagAnalysis * tagAna = new AliTagAnalysis(); 
  tagAna->ChainGridTags(result);
  
  tagAna->CreateXMLCollection(out, runCuts, evtCuts) ;

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::Merge(const char * collectionFile, const char * subFile, const char * outFile) 
{
  // merges files listed in a xml collection 
  // usage Merge(collection, outputFile))
  //              collection: is a xml collection  
  
  Bool_t rv = kFALSE ; 

  if ( strstr(collectionFile, ".xml") == 0 ) {
    AliError("Input collection file must be an \".xml\" file\n") ; 
    return kFALSE ; 
  }

  fTimer.Start() ;

  // Open the file collection 
  printf("*** Create Collection       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",collectionFile);              	
  
  TAlienCollection * collection = TAlienCollection::Open(collectionFile);
  TGridResult* result = collection->GetGridResult("");
  
  Int_t index = 0  ;
  const char * turl ;
  TFileMerger merger ; 
  if (!outFile) {
    TString tempo(collectionFile) ; 
    if ( subFile) 
      tempo.ReplaceAll(".xml", subFile) ; 
    else 
      tempo.ReplaceAll(".xml", "_Merged.root") ; 
    outFile = tempo.Data() ; 
  }
  merger.OutputFile(outFile) ; 

  while ( (turl = result->GetKey(index, "turl")) ) {
    char file[2048] ;
    if ( subFile )
      sprintf(file, "%s#%s", turl, subFile) ; 
    else 
      sprintf(file, "%s", turl) ; 
      
    printf("%s\n", file) ; 
    merger.AddFile(file) ; 
    index++ ;  
  }

  if (index) 
    merger.Merge() ; 
  
  AliInfo(Form("Files merged into %s\n", outFile)) ;
 
  fTimer.Stop();
  fTimer.Print();
  
  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::Process(TChain * chain) 
{
  // process events starting from a chain of esd Trees
  Bool_t rv = kFALSE ; 

  fTimer.Start() ;

  rv = ProcessChain(chain) ; 

  fTimer.Stop();
  fTimer.Print();

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::Process(const char * inFile) 
{ 
  // process the events with an Analysis Task 
  // usage Process(esdFile)
  //              esdFile: is of the form opt?file_lfn 
  Bool_t rv = kFALSE ; 
  AliRunTagCuts   * runCuts = 0x0 ; 
  AliEventTagCuts * evtCuts = 0x0 ;

  rv = Process(inFile, runCuts, evtCuts) ; 

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::Process(const char * inFile, AliRunTagCuts *runCuts, AliEventTagCuts * evtCuts ) 
{
  // process the events with an Analysis Task 
  // usage Process(esdFile, runtagCuts, evtTagCuts)
  //              esdFile: is of the form opt?file_lfn 
  
  Bool_t rv = kFALSE ; 

  fTimer.Start() ;

  TString file(inFile) ; 
  if ( file.Contains("esd?") && file.Contains(".root") ) {
    file.ReplaceAll("esd?", "") ; 
    rv = ProcessEsdFile(file.Data()) ; 

  } else if ( file.Contains("esd?") && file.Contains(".xml") ) { 
    file.ReplaceAll("esd?", "") ; 
    rv = ProcessEsdXmlCollection(file.Data()) ; 

  } else if (file.Contains("tag?") && file.Contains(".root") ) {
    file.ReplaceAll("tag?", "") ; 
    rv = ProcessTagFile(file.Data(), runCuts, evtCuts) ; 

  } else if (file.Contains("tag?") && file.Contains(".xml") ) {
    file.ReplaceAll("tag?", "") ; 
    rv = ProcessTagXmlCollection(file.Data(), runCuts, evtCuts) ; 

  } else { 
    AliError(Form("%s is not a valid file format", inFile)) ; 
    rv = kFALSE ;
  }
  
  fTimer.Stop();
  fTimer.Print();

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::Process(const char * inFile, const char * runCuts, const char * evtCuts) 
{
  // process the events with an Analysis Task 
  // usage Process(esdFile, runtagCuts, evtTagCuts)
  //              esdFile: is of the form opt?file_lfn 
  
  Bool_t rv = kFALSE ; 

  fTimer.Start() ;

  TString file(inFile) ; 
  if ( file.Contains("esd?") && file.Contains(".root") ) {
    file.ReplaceAll("esd?", "") ; 
    rv = ProcessEsdFile(file.Data()) ; 

  } else if ( file.Contains("esd?") && file.Contains(".xml") ) { 
    file.ReplaceAll("esd?", "") ; 
    rv = ProcessEsdXmlCollection(file.Data()) ; 

  } else if (file.Contains("tag?") && file.Contains(".root") ) {
    file.ReplaceAll("tag?", "") ; 
    rv = ProcessTagFile(file.Data(), runCuts, evtCuts) ; 

  } else if (file.Contains("tag?") && file.Contains(".xml") ) {
    file.ReplaceAll("tag?", "") ; 
    rv = ProcessTagXmlCollection(file.Data(), runCuts, evtCuts) ; 

  } else { 
    AliError(Form("%s is not a valid file format", inFile)) ; 
    rv = kFALSE ;
  }
  
  fTimer.Stop();
  fTimer.Print();

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::ProcessChain(TChain * chain) const
{
  // Procees a TChain. 

  Bool_t rv = kTRUE ;

  if (! fTaskList ) {
    AliError("No tasks defined") ; 
    return kFALSE ;
  }

  // Make the analysis manager
  AliAnalysisManager * mgr = new AliAnalysisManager() ;

  // Make tasks 
  // The top input must be common to all top tasks
  TClass * classIn = fTaskInType[0] ; 
  AliAnalysisDataContainer * taskInput  = mgr->CreateContainer("Input  Container", classIn, AliAnalysisManager::kInputContainer) ;
  Int_t index ; 
  for (index = 0; index < fnumberOfTasks; index++) {
    AliAnalysisTask * task = fTaskList[index] ;
    mgr->AddTask(task) ;
  
    // Create containers for input/output
    TClass * classOu = fTaskOuType[index] ; 
    AliAnalysisDataContainer * taskOutput = mgr->CreateContainer("Output Container", classOu, AliAnalysisManager::kOutputContainer) ;
    mgr->ConnectInput (task, 0, taskInput);
    mgr->ConnectOutput(task, 0, taskOutput);
  }
  
  // Open data
  taskInput->SetData(chain);

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    chain->Process(mgr);
  } else 
    rv = kFALSE ; 
  
  return rv ; 
}
 
//______________________________________________________________________
const Bool_t AliAnalysisGoodies::ProcessEsdFile(const char * esdFile) const   
{
  // process the events in a single ESD file with an Analysis Task 
  // usage ProcessLocalEsdFile(esdFile)
  //              esdFile: is the root file (local or in alien) with the ESD Tree ( ex: AliESDs.root) 
 
  Bool_t rv = kTRUE ;  
  
  printf("*** Process       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",esdFile);              	

  // Makes the ESD chain 
  printf("*** Getting the Chain       ***\n");
  TChain* analysisChain = new TChain(fESDTreeName) ;
  analysisChain->AddFile(esdFile);
 
  // Process the events
  rv = ProcessChain(analysisChain) ; 

  return rv;
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::ProcessTagFile(const char * tagFile, AliRunTagCuts *runCuts, AliEventTagCuts *evtCuts) const   
{
  // process the events in a single Tag file with an Analysis Task 
  // usage ProcessLocalEsdFile(tagFile)
  //              tagFile: is the root file (local or in alien) with the Tag Tree (ex: Run102.Event0_100.ESD.tag.root) 
 
  Bool_t rv = kTRUE ;  
  
  if ( !evtCuts && !runCuts ) {
    AliError("No Tag cuts provided") ; 
    return kFALSE ; 
  }
  
  printf("*** Process       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",tagFile);              	

  AliTagAnalysis * tagAna = new AliTagAnalysis(); 
  rv = tagAna->AddTagsFile(tagFile);
  if ( ! rv ) 
    return rv ; 

  // Query the tag file and make the analysis chain
  TChain * analysisChain = new TChain(fESDTreeName)  ;
  analysisChain = tagAna->QueryTags(runCuts, evtCuts);
  
  // Process the events
  rv = ProcessChain(analysisChain) ; 

  return rv;
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::ProcessTagFile(const char * tagFile, const char * runCuts, const char * evtCuts) const   
{
  // process the events in a single Tag file with an Analysis Task 
  // usage ProcessLocalEsdFile(tagFile)
  //              tagFile: is the root file (local or in alien) with the Tag Tree (ex: Run102.Event0_100.ESD.tag.root) 
 
  Bool_t rv = kTRUE ;  
  

  if ( !evtCuts && !runCuts ) {
    AliError("No Tag cuts provided") ; 
    return kFALSE ; 
  }
  
  printf("*** Process       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",tagFile);              	

  AliTagAnalysis * tagAna = new AliTagAnalysis(); 
  rv = tagAna->AddTagsFile(tagFile);
  if ( ! rv ) 
    return rv ; 

  // Query the tag file and make the analysis chain
  TChain * analysisChain = new TChain(fESDTreeName)  ;
  analysisChain = tagAna->QueryTags(runCuts, evtCuts);
  
  // Process the events
 rv = ProcessChain(analysisChain) ; 

  return rv;
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::ProcessEsdXmlCollection(const char * xmlFile) const   
{
  // process the events in a xml ESD collection  with an Analysis Task 
  // usage ProcessLocalEsdFile(xmlFile)
  //              xmlFile: is the local xml file with the ESD collection ( ex: esdCollection.xml) 
 
  Bool_t rv = kTRUE ;  
  
  printf("*** Process       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",xmlFile);              	

  TAlienCollection * collection = TAlienCollection::Open(xmlFile) ; 
  if (! collection) {
    AliError(Form("%s not found", xmlFile)) ; 
    return kFALSE ; 
  }

  TGridResult* result = collection->GetGridResult("");
  TList* analysisfilelist = result->GetFileInfoList();
  
  // Makes the ESD chain 
  printf("*** Getting the Chain       ***\n");
  TChain* analysisChain = new TChain(fESDTreeName);
  analysisChain->AddFileInfoList(analysisfilelist);
 
  // Process the events
  rv = ProcessChain(analysisChain) ; 

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::ProcessTagXmlCollection(const char * xmlFile, AliRunTagCuts *runCuts, AliEventTagCuts * evtCuts) const   
{
  // process the events in a xml ESD collection  with an Analysis Task 
  // usage ProcessLocalEsdFile(xmlFile)
  //              xmlFile: is the local xml file with the tag collection ( ex: tagCollection.xml) 
 
  Bool_t rv = kTRUE ;  
  
  if ( !evtCuts && !runCuts ) {
    AliError("No Tag cuts provided") ; 
    return kFALSE ; 
  }

  printf("*** Process       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",xmlFile);              	
 
  // check if file is local or alien
  if ( gSystem->AccessPathName(xmlFile) ) 
    TGrid::Connect("alien://"); 

  TAlienCollection * collection = TAlienCollection::Open(xmlFile) ; 
  if (! collection) {
    AliError(Form("%s not found", xmlFile)) ; 
    return kFALSE ; 
  }

  TGridResult* result = collection->GetGridResult("");
  AliTagAnalysis * tagAna = new AliTagAnalysis(); 
  tagAna->ChainGridTags(result);
  
  // Query the tag file and make the analysis chain
  TChain * analysisChain = new TChain(fESDTreeName)  ;
  analysisChain = tagAna->QueryTags(runCuts, evtCuts);

  // Process the events
  rv = ProcessChain(analysisChain) ; 

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::ProcessTagXmlCollection(const char * xmlFile, const char * runCuts, const char * evtCuts) const   
{
  // process the events in a xml ESD collection  with an Analysis Task 
  // usage ProcessLocalEsdFile(xmlFile)
  //              xmlFile: is the local xml file with the tag collection ( ex: tagCollection.xml) 
 
  Bool_t rv = kTRUE ;  

 if ( !evtCuts && !runCuts ) {
    AliError("No Tag cuts provided") ; 
    return kFALSE ; 
  }

  printf("*** Process       ***\n");
  printf("***  Wk-Dir = |%s|             \n",gSystem->WorkingDirectory());
  printf("***  Coll   = |%s|             \n",xmlFile);              	
 
  // check if file is local or alien
  if ( gSystem->AccessPathName(xmlFile) ) 
    TGrid::Connect("alien://"); 

  TAlienCollection * collection = TAlienCollection::Open(xmlFile) ; 
  if (! collection) {
    AliError(Form("%s not found", xmlFile)) ; 
    return kFALSE ; 
  }

  TGridResult* result = collection->GetGridResult("");
  AliTagAnalysis * tagAna = new AliTagAnalysis(); 
  tagAna->ChainGridTags(result);
  
  // Query the tag file and make the analysis chain
  TChain * analysisChain = new TChain(fESDTreeName)  ;
  analysisChain = tagAna->QueryTags(runCuts, evtCuts);

  // Process the events
  rv = ProcessChain(analysisChain) ; 

  return rv ; 
}

//______________________________________________________________________
const Bool_t AliAnalysisGoodies::Register( const char * lfndir, const char * pfndir, const char * file) 
{
  // register files already stored in a MSS into the AliEn catalog
  // usage: Register(lfndir, pfndir, pfnFileName)
  //         lfndir : AliEn directory ( ex:  /alice/data/2006/LHC06c/PHOS_TestBeam/ ) 
  //         pfndir : MSS directory   ( ex: /castor/cern.ch/alice/testbeam/phos/2006 )
  //         file   : text file with a list of the file names to be registered

  Bool_t rv = kTRUE ;  
  fTimer.Start() ; 

  ifstream in;
  in.open(file);
  if ( in.bad() ) {
    AliError(Form("Cannot open file %s\n", file)) ; 
    return kFALSE ; 
  }

  TGrid::Connect("alien://");

  char fileName[1024] ;

  while (1) {
    in >> fileName ;
    if (!in.good()) 
      break;
    char lfn[1024] ; 
    
    sprintf(lfn, "%s/%s", lfndir, fileName) ; 
    
    char pfn[1024] ; 
    
    sprintf(pfn, "castor://Alice::CERN::Castor2/%s/%s", pfndir, fileName) ;  
        
    printf("Register %s as %s\n", pfn, lfn) ; 
    
    gGrid->Register(lfn, pfn) ;
  }
  
  fTimer.Stop();
  fTimer.Print();
  
  return rv;
}
 
//______________________________________________________________________
void AliAnalysisGoodies::SetTasks(Int_t nb, AliAnalysisTask ** taskList, TClass ** inputType, TClass ** outputType)
{
  // define a task with its output and input type

  
  fnumberOfTasks= nb; 
  fTaskList   = taskList ;
  fTaskInType = inputType ; 
  fTaskOuType = outputType ; 
}
