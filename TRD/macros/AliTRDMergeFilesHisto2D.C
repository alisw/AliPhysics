#if !defined( __CINT__) || defined(__MAKECINT__)

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2.h>
#include <AliCDBManager.h>
#include <../TRD/AliTRDCalibra.h>

#endif


void AliTRDMergeFilesHisto2D(const char* variablecali, const char* nome, const char* dire){
  //
  // After having simulated and reconstructed events in subrepertories 000%d of dire
  // this macro searchs in the subdirectories the file TRD.calibration.root
  // takes the 2D histos and merge them
  // nome is the name of the file where the sum will be written
  // variablecali can be: CH2d, PH2d, PRF2d
  //



  //The final sum histo
  TH2 *histosum = 0x0;
  Int_t j = 0;


  //The patterns
  const char *name="TRD.calibration.root";
  const char *patterndir="0";
  const char *namesubdir = 0x0;
  char dir[200];

  //Open the current directory
  void * dircu = gSystem->OpenDirectory(dire);
  while((namesubdir = gSystem->GetDirEntry(dircu))) {
    if(strstr(namesubdir,patterndir)) {
      sprintf(dir,"%s/%s",dire,namesubdir);
      printf("process subdirectories: %s\n",dir);
  
      //full name of the file and tree
      char fullname[255] = "";
      
      sprintf(fullname,"%s/%s",dir,name);
      printf("Process file: %s\n",fullname);
      TFile *file = (TFile *) TFile::Open(fullname,"READ");
      if(!file) continue;
      TH2 *histo = (TH2 *) file->Get(variablecali);
      histo->SetDirectory(0);
      if(!histo) continue;
      if(j == 0) histosum = histo;
      else histosum->Add(histo,1);
      j++;
      delete file;
    }// if patterndir
  }//loop in the current directory


  //Open a file to put the histosum
  TFile *fout = TFile::Open((const char*) nome ,"UPDATE");
  fout->WriteTObject(histosum,histosum->GetName(),(Option_t *) "kOverWrite");
  fout->Close();

}
