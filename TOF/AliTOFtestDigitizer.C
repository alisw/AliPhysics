////////////////////////////////////////////////////////////////////////
//
// name: AliTOFtestDigitizer
// date: 11-VI-2002
// last update: 11-VI-2002
// author: F. Pierella | pierella@bo.infn.it
// version: 1.0
//
// description: 
//       creates digits from sdigits for TOF detector
//       stores sdigits in separate file (or in the source file
//       with sdigits). Stores gAlice object and copies TE to the
//       file with digits
//
// input:
//       TString fileNameSDigits ... input file with sdigits
//       TString fileNameDigits ... output file with digits
//       Int_t nEvents  ... how many events to process
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"
#include "AliTOFDigitizer.h"
#include "../STEER/AliRunDigitizer.h"
#include "../STEER/AliDigitizer.h"
#include "TStopwatch.h"
#endif

//#include "AliHits2SDigits.C"

TFile* OpenFile(TString fileName);
Bool_t ImportgAlice(TFile *file);
void AliCopyN(TString inputFile, TString outputFile);
void AliCopy(TFile *inputFile, TFile *outputFile);
Int_t gDEBUG = 1;


Int_t AliTOFtestDigitizer(TString fileNameDigits="digits.root", 
			TString fileNameSDigits="rfio:sdigits.root", 
			Int_t nEvents = 1, Int_t iTOF = 1, Int_t iCopy = 1)
{
// delete the current gAlice object, the one from input file
//  will be used

  if(gAlice){
    delete gAlice;
    gAlice = 0;
  } // end if gAlice
  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  manager->SetInputStream(0,fileNameSDigits.Data());
  if (fileNameDigits != "") {
    if (iCopy) {
      AliCopyN(fileNameSDigits,fileNameDigits);
    }
    manager->SetOutputFile(fileNameDigits);
  }
  manager->SetNrOfEventsToWrite(nEvents);
  if (iTOF) AliTOFDigitizer *dTOF  = new AliTOFDigitizer(manager);
  TStopwatch timer;
  timer.Start();
  manager->Exec("deb all");
  timer.Stop(); 
  timer.Print();
  delete manager;
}


////////////////////////////////////////////////////////////////////////
void AliCopyN(TString inputFileName, TString outputFileName) {
// copy some objects

  TFile *inputFile = OpenFile(inputFileName);
  if (!inputFile) return;

  TFile *outputFile = TFile::Open(outputFileName.Data(),"update");
  if (!outputFile->IsOpen()) {
    cerr<<"Can't open "<<outputFileName.Data()<<" !\n";
    return;
  }
  if (!ImportgAlice(inputFile)) return;
  AliCopy(inputFile, outputFile);
  inputFile->Close();
  delete inputFile;
  outputFile->Close();
  delete outputFile;
}

////////////////////////////////////////////////////////////////////////
TFile* OpenFile(TString fileName) {
  // open file fileName
  TFile *file = TFile::Open(fileName.Data());
  if (!file->IsOpen()) {
    cerr<<"Can't open "<<fileName.Data()<<" !\n";
    return 0;
  }
  return file;
}

////////////////////////////////////////////////////////////////////////
Bool_t ImportgAlice(TFile *file) {
  // read in gAlice object from the file
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice)  return kFALSE;
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
void AliCopy(TFile *inputFile, TFile *outputFile) {
// copy some objects
  
// copy gAlice
  if (gDEBUG) cout<<"Copy gAlice: ";
  outputFile->cd();
  gAlice->Write();
  if (gDEBUG) cout<<"done"<<endl;
 
  TTree *treeE  = gAlice->TreeE();
  if (!treeE) {
    cerr<<"No TreeE found "<<endl;
    return;
  }
  
// copy TreeE
  if (gDEBUG) cout<<"Copy TreeE: ";
  AliHeader *header = new AliHeader();
  treeE->SetBranchAddress("Header", &header);
  treeE->SetBranchStatus("*",1);  
  TTree *treeENew =  treeE->CloneTree();
  treeENew->Write();
  if (gDEBUG) cout<<"done"<<endl;
   
// copy AliceGeom
  if (gDEBUG) cout<<"Copy AliceGeom: ";
  TGeometry *AliceGeom = static_cast<TGeometry*>(inputFile->Get("AliceGeom"));
  if (!AliceGeom) {
    cerr<<"AliceGeom was not found in the input file "<<endl;
    return;
  }
  AliceGeom->Write();
  if (gDEBUG) cout<<"done"<<endl;
   
}

