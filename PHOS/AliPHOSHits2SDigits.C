////////////////////////////////////////////////////////////////////////
//
// name: AliPHOSHits2SDigits
// date: 4.4.2002
// last update: 4.4.2002
// author: Jiri Chudoba
// version: 1.0
//
// description: 
//       creates sdigits for digits for PHOS.
//       stores sdigits in separate file (or in the source file
//       with hits). Stores gAlice object and copies TE to the
//       file with sdigits.
//
// input:
//        TString fileNameHits ... input file with hits
//        TString fileNameSDigits ... output file with sdigits
//
// History:
//
// 04.04.02 - first version
// 
////////////////////////////////////////////////////////////////////////


Int_t AliPHOSHits2SDigits(TString
fileNameSDigits="PHOS.sdigits.root",TString fileNameHits="rfio:galice.root")

{
  TFile *fileSDigits;
  fileSDigits = Init(fileNameSDigits, fileNameHits);
  if (!fileSDigits) return 1;

  AliPHOSSDigitizer *sdPHOS = new AliPHOSSDigitizer(fileNameHits.Data(),"PHOS");

  TStopwatch timer;
  timer.Start();

  gAlice->MakeTree("S",fileSDigits);
  sdPHOS->ExecuteTask("deb all");

  timer.Stop(); 
  timer.Print();

  fileSDigits->Close();
  delete fileSDigits;

}
 

////////////////////////////////////////////////////////////////////////
TFile* Init(TString fileNameSDigits, TString fileNameHits) {
// open input file, read in gAlice, prepare output file
  if (gAlice) delete gAlice;
  gAlice = 0;

  Bool_t sameFiles = kFALSE;
  if (fileNameSDigits == fileNameHits || fileNameSDigits == "") sameFiles = kTRUE;

  TString fileMode = "read";
  if (sameFiles) fileMode = "update";

  TFile *fileHits =  OpenFile(fileNameHits.Data(),fileMode.Data());
  if (!fileHits) return 0;
  if (!ImportgAlice(fileHits)) return 0;
  if (!sameFiles) return gAlice->InitTreeFile("S",fileNameSDigits.Data());
  return fileHits;

}

////////////////////////////////////////////////////////////////////////
TFile* OpenFile(TString fileName, TString fileMode) {
// open file fileName
  TFile *file = TFile::Open(fileName.Data(),fileMode.Data());
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
