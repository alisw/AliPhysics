////////////////////////////////////////////////////////////////////////
//
//  Class to steer several steps in merging:
//       - extraction of vertex from bgr event
//       - event simulation
//       - creation of sdigits
//       - creation of digits
//       - analysis
//                  
//  Author: Jiri Chudoba (CERN), 2002
//
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>

// --- ROOT system ---

#include "TTask.h"
#include "TString.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TArrayF.h"
#include "TSystem.h"
#include "TGeometry.h"
#include "TInterpreter.h"

// --- AliRoot header files ---

#include "STEER/AliRun.h"
#include "STEER/AliHeader.h"
#include "STEER/AliGenEventHeader.h"
#endif

class AliTRDdigitizer;
class AliPHOSSDigitizer;
class AliTOFSDigitizer;

class AliMergeSteer: public TTask {

public:
  AliMergeSteer(const Text_t* name="AliMergeSteer",
		  const Text_t* title="AliMergeSteer");
//  AliMergeSteer(Option_t *option);
  virtual ~AliMergeSteer();
  Bool_t SetDetectorFlag(Option_t* option, Int_t flag);
  Int_t GetDetectorFlag(Option_t* option);  
//  void Print();
  virtual void Exec(const Option_t *option);
  Bool_t ImportgAlice(TFile *file);

  Bool_t ExtractVertex(TString fn, Int_t eventNr);
  void PrintVertex(TArrayF &vertex);
  void ExportVertex(TArrayF &vertex);

  Bool_t Simulate();
  Bool_t CreateSDigits();
  Bool_t Merge();
  Bool_t NoMerge();
  Bool_t ITSFastPoints(const char *outputFile, const char *inputFile);
  Bool_t RecoMerged();
  Bool_t RecoSignalOnly();
  Bool_t CmpMerged();
  Bool_t CmpSignalOnly();
  
  Bool_t AliCopy(TFile *inputFile, TFile *outputFile);
  Bool_t AliCopy(TString inputFileName, TString outputFileName);
  
  
  
// file names setter/getters  
  void SetFileNameHits(TString fileName) {fFileNameHits = fileName;}
  void SetFileNameSDigits(TString fileName) {fFileNameSDigits = fileName;}
  void SetFileNameBgrHits(TString fileName) {fFileNameBgrHits = fileName;}
  void SetFileNameBgrSDigits(TString fileName) {fFileNameBgrSDigits = fileName;}
  void SetFileNameDigitsMerged(TString fileName) {fFileNameDigitsMerged = fileName;}
  void SetFileNameDigitsSignalOnly(TString fileName) {fFileNameDigitsSignalOnly = fileName;}
  TString FileNameHits() {return fFileNameHits;}
  TString FileNameSDigits() {return fFileNameSDigits;}
  TString FileNameBgr() {return fFileNameBgrSDigits;}
  TString FileNameDigitsMerged() {return fFileNameDigitsMerged;}
  TString FileNameDigitsSignalOnly() {return fFileNameDigitsSignalOnly;}

// flags
  void SetFlagSim(Bool_t flag) {fFlagSim = flag;}
  void SetFlagSDigits(Bool_t flag) {fFlagSDigits = flag;}
  void SetFlagMerge(Bool_t flag) {fFlagMerge = flag;}
  void SetFlagDigitsSignalOnly(Bool_t flag) {fFlagDigitsSignalOnly = flag;}
  void SetFlagRecoMerged(Bool_t flag) {fFlagRecoMerged = flag;}
  void SetFlagRecoSignalOnly(Bool_t flag) {fFlagRecoSignalOnly = flag;}
  void SetFlagCmpMerged(Bool_t flag) {fFlagCmpMerged = flag;}
  void SetFlagCmpSignalOnly(Bool_t flag) {fFlagCmpSignalOnly = flag;}
  
// event numbers
  void SetNEvents(Int_t nEvents) {fNEvents = nEvents;}
  void SetBgrEventNr(Int_t i) {fBgrEventNr = i;}
  
private:  

  TString fFileNameHits;
  TString fFileNameSDigits;
  TString fFileNameBgrHits;
  TString fFileNameBgrSDigits;
  TString fFileNameDigitsMerged;
  TString fFileNameDigitsSignalOnly;
  TString fFileNameCmpMerged;
  TString fFileNameCmpSignalOnly;
  
  Bool_t fFlagSim;
  Bool_t fFlagSDigits;
  Bool_t fFlagMerge;
  Bool_t fFlagDigitsSignalOnly;
  Bool_t fFlagRecoMerged;
  Bool_t fFlagRecoSignalOnly;
  Bool_t fFlagCmpMerged;
  Bool_t fFlagCmpSignalOnly;
  Bool_t fCopy2S;
  Bool_t fCopy2D;
  
  Int_t fNEvents;
  Int_t fFirstEvent;
  Int_t fBgrEventNr;
  
  Int_t fDEBUG;
  
// flags for detectors - determine which detector will be used
// during merging.
//  0 - do not use
//  1 - process the whole detector
//  2 - process only active region (region of interest on)  
  
  Int_t   fFMD;
  Int_t   fITS;
  Int_t   fMUON;
  Int_t   fPHOS;
  Int_t   fPMD;
  Int_t   fHMPID;
  Int_t   fT0;
  Int_t   fTOF;
  Int_t   fTPC;
  Int_t   fTRD;
  Int_t   fZDC;
  Int_t   fEMCAL;
  
};


////////////////////////////////////////////////////////////////////////
//
// AliMergeSteer.cxx
//
//
////////////////////////////////////////////////////////////////////////

AliMergeSteer::AliMergeSteer(const Text_t *name, const Text_t* title) : TTask(name,title)
{
//
// default ctor
//
  fFileNameHits="galice.root";
  fFileNameSDigits = "sdigits.root";
  fFileNameBgrSDigits = "bgr.sdigits.root";
  fFileNameBgrHits = "bgr.hits.root";
  fFileNameDigitsMerged = "digits.root";
  fFileNameDigitsSignalOnly = "digits.signal.root";
  fFileNameCmpMerged = "CmpGaRS_merged.root";
  fFileNameCmpSignalOnly = "CmpGaRS_signal.root";

  fFlagSim = kFALSE;
  fFlagSDigits = kFALSE;
  fFlagMerge = kFALSE;
  fFlagDigitsSignalOnly = kFALSE;
  fFlagRecoMerged = kFALSE;
  fFlagRecoSignalOnly = kFALSE;
  fFlagCmpMerged = kFALSE;
  fFlagCmpSignalOnly = kFALSE;    
  fCopy2S = kTRUE;
  fCopy2D = kTRUE;
  
  fNEvents = 1;
  fFirstEvent = 0;
  fBgrEventNr = 0;

  fDEBUG = 1;

  fFMD = 0;
  fITS = 0;
  fMUON = 0;
  fPHOS = 0;
  fPMD = 0;
  fHMPID = 0;
  fT0 = 0;
  fTOF = 0;
  fTPC = 0;
  fTRD = 0;
  fZDC = 0;
  fEMCAL = 0;
}

////////////////////////////////////////////////////////////////////////
AliMergeSteer::~AliMergeSteer()
{
// default dtor

}
////////////////////////////////////////////////////////////////////////
void AliMergeSteer::Exec(Option_t* option)
{
  Bool_t rc = kTRUE;
  if (gAlice) delete gAlice;
  gAlice = 0;

  if (fFlagSim) rc = Simulate();
  if (!rc) {Error("Exec","Simulate wrong"); return;}

  if (fFlagSDigits) rc = CreateSDigits();
  if (!rc) {Error("Exec","CreateSDigits wrong"); return;}

  if (fFlagSDigits && fCopy2S) rc = AliCopy(fFileNameHits,fFileNameSDigits);
  if (!rc) {Error("Exec","AliCopy to SD wrong"); return;}

  if (fFlagMerge && fCopy2D) rc = AliCopy(fFileNameHits,fFileNameDigitsMerged);
  if (!rc) {Error("Exec","AliCopy to DigitsMerged wrong"); return;}

  if (fFlagMerge) rc = Merge();
  if (!rc) {Error("Exec","Merge wrong"); return;}

  if (fFlagDigitsSignalOnly) rc = NoMerge();
  if (!rc) {Error("Exec","NoMerge wrong"); return;}
  
  if (fFlagRecoMerged) rc = RecoMerged();
  if (!rc) {Error("Exec","RecoMerged wrong"); return;}
  
  if (fFlagRecoSignalOnly) rc = RecoSignalOnly();
  if (!rc) {Error("Exec","RecoSignalOnly wrong"); return;}

  if (fFlagCmpMerged) rc = CmpMerged();
  if (!rc) {Error("Exec","CmpMerged wrong"); return;}
  
  if (fFlagCmpSignalOnly) rc = CmpSignalOnly();
  if (!rc) {Error("Exec","CmpSignalOnly wrong"); return;}
  
 


}
////////////////////////////////////////////////////////////////////////
Int_t AliMergeSteer::GetDetectorFlag(Option_t* option)
{
//
// return current flag value for a given detector
//
  if (strstr(option,"FMD")) {
    return fFMD;
  } else if (strstr(option,"ITS")) {
    return fITS;
  } else if (strstr(option,"MUON")) {
    return fMUON;
  } else if (strstr(option,"PHOS")) {
    return fPHOS;
  } else if (strstr(option,"PMD")) {
    return fPMD;
  } else if (strstr(option,"HMPID")) {
    return fHMPID;
  } else if (strstr(option,"T0")) {
    return fT0;
  } else if (strstr(option,"TOF")) {
    return fTOF;
  } else if (strstr(option,"TPC")) {
    return fTPC;
  } else if (strstr(option,"TRD")) {
    return fTRD;
  } else if (strstr(option,"ZDC")) {
    return fZDC;
  } else if (strstr(option,"EMCAL")) {
    return fEMCAL;
  } else {
    cerr<<"Unknown detector required."<<endl;
    return -1;
  }
}
////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::SetDetectorFlag(Option_t* option, Int_t flag)
{
//
// return current flag value for a given detector
//
  if (strstr(option,"FMD")) {
    fFMD = flag;
  } else if (strstr(option,"ITS")) {
    fITS = flag;
  } else if (strstr(option,"MUON")) {
    fMUON = flag;
  } else if (strstr(option,"PHOS")) {
    fPHOS = flag;
  } else if (strstr(option,"PMD")) {
    fPMD = flag;
  } else if (strstr(option,"HMPID")) {
    fHMPID = flag;
  } else if (strstr(option,"T0")) {
    fT0 = flag;
  } else if (strstr(option,"TOF")) {
    fTOF = flag;
  } else if (strstr(option,"TPC")) {
    fTPC = flag;
  } else if (strstr(option,"TRD")) {
    fTRD = flag;
  } else if (strstr(option,"ZDC")) {
    fZDC = flag;
  } else if (strstr(option,"EMCAL")) {
    fEMCAL = flag;
  } else {
    cerr<<"Unknown detector required."<<endl;
    return kFALSE;
  }
  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::ImportgAlice(TFile *file) {
// read in gAlice object from the file
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice)  return kFALSE;
  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::ExtractVertex(TString fn, Int_t eventNr)
{
// Open file with TPC geom and digits
  TFile *file=TFile::Open(fn);
  if (!file->IsOpen()) {cerr<<"Cannnot open "<<fn.Data()<<" !\n"; return kFALSE;}
  if (!(gAlice=(AliRun*)file->Get("gAlice"))) {
    cerr<<"gAlice was not found in "<<fn.Data()<<" !\n";
    return kFALSE;
  }
  
  AliHeader *header = gAlice->GetHeader();
  if (!header) {
    cerr<<"header was not found in "<<fn.Data()<<" !\n";
    return kFALSE;
  } 
  AliGenEventHeader* genEventHeader = header->GenEventHeader();
  if (!genEventHeader) {
    cerr<<"GenEventHeader was not found in "<<fn.Data()<<" !\n";
    return kFALSE;
  } 

  TArrayF primaryVertex(3);
  genEventHeader->PrimaryVertex(primaryVertex);
  PrintVertex(primaryVertex);
  ExportVertex(primaryVertex);
//  delete header;

// Following two lines should be there, but ....
//  delete genEventHeader;
//  delete gAlice;
  gAlice = 0;
  file->Close();
    
  return kTRUE;

}
////////////////////////////////////////////////////////////////////////
void AliMergeSteer::PrintVertex(TArrayF &primaryVertex) 
{
  cout <<"CONFIG_VERTEX: "
       <<primaryVertex[0]<<" "
       <<primaryVertex[1]<<" "
       <<primaryVertex[2]<<" "<<endl;
  return;
} 

////////////////////////////////////////////////////////////////////////
void AliMergeSteer::ExportVertex(TArrayF &primaryVertex) 
{
  char vertexAsString[30];
  sprintf(vertexAsString,"%f",primaryVertex[0]);
  gSystem->Setenv("CONFIG_VERTEX_X",vertexAsString);
  sprintf(vertexAsString,"%f",primaryVertex[1]);
  gSystem->Setenv("CONFIG_VERTEX_Y",vertexAsString);
  sprintf(vertexAsString,"%f",primaryVertex[2]);
  gSystem->Setenv("CONFIG_VERTEX_Z",vertexAsString);
  return;
}
////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::Simulate()
{
  char cmd[100];
  TFile *f;

  TDirectory *saveDir = gDirectory;
  cout<< "The original directory is: ";
  saveDir->pwd();
  if (!ExtractVertex(fFileNameBgrHits,fBgrEventNr)) {
    cerr<<" ExtractVertexAndSimulateSignal:  Error in ExtractVertex"<<endl;
    return kFALSE;
  }
  saveDir->cd();
  new AliRun("gAlice","Signal for Merging");
  gAlice->Init("Config.C");
//    gSystem->Setenv("CONFIG_FILE",fileNameSigHits);
  TStopwatch timer;
  timer.Start();
  gAlice->Run(fNEvents);
  timer.Stop();
  cout<<"Simulation of "<<fNEvents<<" of signal took: "<<endl;
  timer.Print();
  delete gAlice;
  gAlice = 0;
  f = static_cast<TFile *>(gROOT->FindObject("galice.root"));
  if (f) f->Close();
  f = 0;    
  sprintf(cmd,"mv galice.root %s",fFileNameHits.Data());
  gSystem->Exec(cmd);
  return kTRUE;
 
}

////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::CreateSDigits()
{
  char macroName[200];
  char funcName[200]; 
  sprintf(macroName,"AliHits2SDigits.C");
  sprintf(funcName,
	  "AliHits2SDigits(\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d);",
	  fFileNameSDigits.Data(),fFileNameHits.Data(),
	  fNEvents,0,fITS,fTPC,fTRD,fPHOS,fTOF,0);
  cerr<<"I'll do: "<<funcName<<endl;
  gROOT->LoadMacro(macroName);
  if (fDEBUG) cerr<<"I'll do: "<<funcName<<endl;
  gInterpreter->ProcessLine(funcName);
  if (fDEBUG) cerr<<"SDigits created"<<endl;
  return kTRUE;
}


////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::Merge()
{
  char macroName[200];
  char funcName[200]; 
  cerr<<"Start merging"<<endl;
  sprintf(macroName,"MergeV1.C");
  sprintf(funcName,
	  "Merge(\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d);",
	  fFileNameDigitsMerged.Data(),fFileNameSDigits.Data(),
	  fFileNameBgrSDigits.Data(),
	  fNEvents,fITS,fTPC,fTRD,fPHOS,fMUON,fHMPID,0);
  cerr<<"I'll do: "<<funcName<<endl;
  gROOT->LoadMacro(macroName);
  if (fDEBUG) cerr<<"I'll do: "<<funcName<<endl;
  gInterpreter->ProcessLine(funcName);  
  if (fDEBUG) cerr<<"Merging done"<<endl;

//  return kTRUE;
// add ITS fast points, no merging yet
//  return ITSFastPoints(fFileNameDigitsMerged.Data(),
//		       fFileNameHits.Data());

}
////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::NoMerge()
{
  char macroName[200];
  char funcName[200]; 
  cerr<<"Start NoMerging"<<endl;
  sprintf(macroName,"AliSDigits2Digits.C");
  sprintf(funcName,
	  "AliSDigits2Digits(\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d);",
	  fFileNameDigitsSignalOnly.Data(),fFileNameSDigits.Data(),
	  fNEvents,fITS,fTPC,fTRD,fPHOS,fMUON,fHMPID,0);
  if (fDEBUG) cerr<<"I'll do: "<<funcName<<endl;
  gROOT->LoadMacro(macroName);
  gInterpreter->ProcessLine(funcName);  
  if (fDEBUG) cerr<<"NoMerging done"<<endl;
//  return kTRUE;
// add ITS fast points, no merging
//  return ITSFastPoints(fFileNameDigitsSignalOnly.Data(),
//		       fFileNameHits.Data());

}
////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::ITSFastPoints(const char *outputFile, const char *inputFile) {

  char macroName[200];
  char funcName[200]; 
  sprintf(macroName,"AliITSHits2SDR.C");
  sprintf(funcName,"AliITSH2FR2files(\"%s\",\"%s\");",
	  inputFile,  outputFile);
  if (fDEBUG) cerr<<"I'll do: "<<funcName<<endl;
  gROOT->LoadMacro(macroName);
  gInterpreter->ProcessLine(funcName);  
  if (fDEBUG) cerr<<"ITSFastPoints done"<<endl;

  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::RecoMerged()
{
//
//
//
  char macroName[200];
  char funcName[200]; 
  cerr<<"Start RecoMerged"<<endl;
  sprintf(macroName,"AliBarrelRecoV3.C");
  sprintf(funcName,"AliBarrelRecoMerged(%d);",fNEvents);
  gROOT->LoadMacro(macroName);
  if (fDEBUG) cerr<<"I'll do: "<<funcName<<endl;
  gInterpreter->ProcessLine(funcName);
  if (fDEBUG) cerr<<"RecoMerged done"<<endl;
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::RecoSignalOnly()
{
//
//
//
  char macroName[200];
  char funcName[200]; 
  cerr<<"Start RecoSignalOnly"<<endl;
  sprintf(macroName,"AliBarrelRecoNoITSClass.C");
  sprintf(funcName,"AliBarrelReco(%d);",fNEvents);
  gROOT->LoadMacro(macroName);
  if (fDEBUG) cerr<<"I'll do: "<<funcName<<endl;
  gInterpreter->ProcessLine(funcName);  
  if (fDEBUG) cerr<<"RecoSignalOnly done"<<endl;
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::CmpMerged()
{
//
//
//
  char macroName[200];
  char funcName[200]; 
  cerr<<"Start CmpMerged"<<endl;
  sprintf(macroName,"CmpGaRS.C");
  sprintf(funcName,
	  "CmpGaRS(%d,%d,\"%s\",\"AliTPCtracks_merged.root\",\"%s\");",
	  fNEvents, fFirstEvent, fFileNameHits.Data(),
	  fFileNameCmpMerged.Data());
  gROOT->LoadMacro(macroName);
  if (fDEBUG) cerr<<"I'll do: "<<funcName<<endl;
  gInterpreter->ProcessLine(funcName);  
  if (fDEBUG) cerr<<"CmpMerged done"<<endl;
  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::CmpSignalOnly()
{
//
//
//
  char macroName[200];
  char funcName[200]; 
  cerr<<"Start CmpSignalOnly"<<endl;
  sprintf(macroName,"CmpGaRS.C");
  sprintf(funcName,
	  "CmpGaRS(%d,%d,\"%s\",\"AliTPCtracks.root\",\"%s\");",
	  fNEvents, fFirstEvent, fFileNameHits.Data(),
	  fFileNameCmpSignalOnly.Data());
  gROOT->LoadMacro(macroName);
  if (fDEBUG) cerr<<"I'll do: "<<funcName<<endl;
  gInterpreter->ProcessLine(funcName);  
  if (fDEBUG) cerr<<"CmpSignalOnly done"<<endl;
  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::AliCopy(TFile *inputFile, TFile *outputFile)
{
//
// copy gAlice object, AliceGeom and TreeE
//

// copy gAlice
  if (fDEBUG) cout<<"Copy gAlice: ";
  outputFile->cd();
  if (gAlice) {
    Error("AliCopy",
	  "gAlice must be deleted before AliCopy is called.");
    return kFALSE;
  }
  if (!ImportgAlice(inputFile)) return kFALSE;
  gAlice->Write();
  if (fDEBUG) cout<<"done"<<endl;


// copy TreeE
  TTree *treeE  = gAlice->TreeE();
  if (!treeE) {
    cerr<<"No TreeE found "<<endl;
    return kFALSE;
  }      
  if (fDEBUG) cout<<"Copy TreeE: ";
  AliHeader *header = new AliHeader();
  treeE->SetBranchAddress("Header", &header);
  treeE->SetBranchStatus("*",1);
  TTree *treeENew =  treeE->CloneTree();
  treeENew->Write();
  if (fDEBUG) cout<<"done"<<endl;

// copy AliceGeom
  if (fDEBUG) cout<<"Copy AliceGeom: ";
  TGeometry *AliceGeom = static_cast<TGeometry*>(inputFile->Get("AliceGeom"));
  if (!AliceGeom) {
    cerr<<"AliceGeom was not found in the input file "<<endl;
    return kFALSE;
  }
  AliceGeom->Write();
  if (fDEBUG) cout<<"done"<<endl;

  delete gAlice;
  gAlice = 0;

  return kTRUE;

}

////////////////////////////////////////////////////////////////////////
Bool_t AliMergeSteer::AliCopy(TString inputFileName, TString outputFileName)
{
//
// open iput and output files, 
// ask to copy gAlice object, AliceGeom and TreeE
// close input and ouput files
//
  if (fDEBUG) {
    cout<<"AliCopy: will copy gAlice from "<<inputFileName.Data()<<" to "
	<<outputFileName.Data()<<endl;
  }

  TFile *inputFile = TFile::Open(inputFileName.Data());
  if (!inputFile->IsOpen()) {
    cerr<<"Can't open "<<inputFileName.Data()<<" !\n";
    return kFALSE;
  }

  TFile *outputFile = TFile::Open(outputFileName.Data(),"UPDATE");
  if (!outputFile->IsOpen()) {
    cerr<<"Can't open "<<outputFileName.Data()<<" !\n";
    return kFALSE;
  }

  AliCopy(inputFile, outputFile);

  inputFile->Close();
  outputFile->Close();

  if (fDEBUG) {
    cout<<"AliCopy copied gAlice from "<<inputFileName.Data()<<" to "
	<<outputFileName.Data()<<endl;
  }

  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
