#if !defined( __CINT__) || defined(__MAKECINT__)

#include <vector>
#include <TChain.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TString.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TH1F.h>
#include <AliCDBManager.h>
#include "../TRD/AliTRDCalibra.h"

#endif


Bool_t AliTRDCreateVector2DSum(const char* variablecali, const char* noime, const char* dire, const char* namefile){
  //
  // After having simulated and reconstructed events in subrepertories 000%d of dire
  // this macro searchs in the subdirectories the file TRD.calibration.root
  // takes the vectors and merge them
  // variablecali can be treeCH2d, treePH2d or treePRF2d
  // noime is the name of the "sum" tree
  // namefile is the name of the file where noime will be written
  //


  // The tree
  TChain *treeChain = new TChain(variablecali);
   

  //Variables
  TObjArray *vectorplace = new TObjArray();
  TObjArray *whereinthechain = new TObjArray();


  //TRDCalibra
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  AliTRDCalibra *calibra = AliTRDCalibra::Instance();

  //The patterns
  const char *name="TRD.calibration.root";
  const char *patterndir="0";
  const char *namesubdir = 0x0;
  char dir[200];

  //Open the directory dire
  void * dircu = gSystem->OpenDirectory(dire);
  while((namesubdir = gSystem->GetDirEntry(dircu))) {
    if(strstr(namesubdir,patterndir)) {
      sprintf(dir,"%s/%s",dire,namesubdir);
      printf("process subdirectories: %s\n",dir);
      
      //full name of the file and tree
      char fullname[255] = "";
      
      sprintf(fullname,"%s/%s",dir,name);
      printf("Process file: %s\n",fullname);
      TFile *file = new TFile(fullname,"READ");
      if(!file) continue;
      TTree *treecurrent = (TTree *) file->Get(variablecali);
      if(!treecurrent) continue;
      gDirectory = gROOT;
      TObjArray *vectorcurrent = calibra->ConvertTreeVector(treecurrent);
      printf("Size of the current tree: %d\n",(Int_t) vectorcurrent->GetEntriesFast());
      printf("Size of whereinthechain: %d\n",(Int_t) whereinthechain->GetEntriesFast());
      printf("Size of the chain: %d\n", (Int_t) treeChain->GetEntries());
      Int_t j = (Int_t) treeChain->GetEntries();
	for(Int_t jui = 0; jui < (Int_t) vectorcurrent->GetEntriesFast(); jui++){
	  //Search if already found
	  Int_t place = calibra->SearchInTreeVector(vectorplace,((AliTRDCalibra::AliTRDPlace *)vectorcurrent->At(jui))->GetPlace());
	  //Create a new element in the two std vectors
	  if(place == -1){
	    AliTRDCalibra::AliTRDPlace *placejui = new AliTRDCalibra::AliTRDPlace();
	    placejui->SetPlace(j+jui);
	    TObjArray *chainplace = new TObjArray();
	    chainplace->Add((TObject *) placejui);
	    vectorplace->Add((TObject *) (vectorcurrent->At(jui)));
	    whereinthechain->Add((TObject *) chainplace);
	  }
	  //Update the element at the place "place" in the std vector whereinthechain
	  else {
	    AliTRDCalibra::AliTRDPlace *placejui = new AliTRDCalibra::AliTRDPlace();
	    placejui->SetPlace((j+jui));
	    TObjArray *chainplace = ((TObjArray *) whereinthechain->At(place));
	    chainplace->Add((TObject *) placejui);
	    whereinthechain->AddAt((TObject *)chainplace, place);
	  }
	} 
	treeChain->AddFile(fullname);
	delete file;
    }//if pattern
  }// in the sub
  
  //Take care of the profile
  const char* pattern = "P";
  TTree *tree;

  if(!strstr(variablecali,pattern)){
    //Ready to read the chain
    TH1F *his = 0x0;
    treeChain->SetBranchAddress("histo",&his);
    
    //Initialise the final tree
    Int_t group = -1;
    TH1F *histsum = 0x0;
    
    tree = new TTree(noime,noime);
    tree->Branch("groupnumber",&group,"groupnumber/I");
    tree->Branch("histo","TH1F",&histsum,32000,0);
    
    //Init histsum
    if(treeChain->GetEntries() < 1) return kFALSE; 
    
    printf("FINAL Size of the chain: %d\n", (Int_t) treeChain->GetEntriesFast());
    printf("FINAL Size of vectorplace: %d\n", (Int_t) vectorplace->GetEntriesFast());
    for(Int_t h = 0; h < (Int_t) vectorplace->GetEntriesFast(); h++){
      group = ((AliTRDCalibra::AliTRDPlace *)vectorplace->At(h))->GetPlace();
      TObjArray *chainplace = ((TObjArray *)whereinthechain->At(h));
      treeChain->GetEntry(((AliTRDCalibra::AliTRDPlace *)chainplace->At(0))->GetPlace());
      //Init for the first time
      if(h == 0)  {
	histsum = new TH1F("","",his->GetXaxis()->GetNbins(),his->GetXaxis()->GetBinLowEdge(1),his->GetXaxis()->GetBinUpEdge(his->GetXaxis()->GetNbins()));
	histsum->Sumw2();
      }
      //Reset for each new group
      histsum->SetEntries(0);
      for(Int_t l = 0; l <= histsum->GetXaxis()->GetNbins(); l++){
	histsum->SetBinContent(l,0.0);
	histsum->SetBinError(l,0.0);
      }
      histsum->Add(his,1);
      if((Int_t) chainplace->GetEntriesFast() > 1){
	for(Int_t s = 1; s < (Int_t)chainplace->GetEntriesFast(); s++){
	  treeChain->GetEntry(((AliTRDCalibra::AliTRDPlace *)chainplace->At(s))->GetPlace());
	  histsum->Add(his,1);
	}
      }
      tree->Fill();
    }
  }
  else {
    //Ready to read the chain
    TGraphErrors *his = 0x0;
    treeChain->SetBranchAddress("histo",&his);
    
    //Initialise the final tree
    Int_t group = -1;
    TGraphErrors *histsum = 0x0;
    Double_t *xref = 0x0;
   
    
    tree = new TTree(noime,noime);
    tree->Branch("groupnumber",&group,"groupnumber/I");
    tree->Branch("histo","TGraphErrors",&histsum,32000,0);
    
    //Init histsum
    if(treeChain->GetEntries() < 1) return kFALSE; 
    
    printf("FINAL Size of the chain: %d\n", (Int_t) treeChain->GetEntriesFast());
    printf("FINAL Size of vectorplace: %d\n", (Int_t) vectorplace->GetEntriesFast());
    for(Int_t h = 0; h < (Int_t) vectorplace->GetEntriesFast(); h++){
      group = ((AliTRDCalibra::AliTRDPlace *)vectorplace->At(h))->GetPlace();
      TObjArray *chainplace = ((TObjArray *)whereinthechain->At(h));
      treeChain->GetEntry(((AliTRDCalibra::AliTRDPlace *)chainplace->At(0))->GetPlace());
      //Init for the fisrt time
      Int_t nbins = his->GetN();
      Double_t *x;
      x = new Double_t[nbins];
      xref = his->GetX();
      Double_t *ex;
      ex = new Double_t[nbins];
      Double_t *y;
      y = new Double_t[nbins];
      Double_t *ey;
      ey = new Double_t[nbins];
      
      for(Int_t lo = 0; lo < nbins; lo++){
	x[lo] = xref[lo];
	ex[lo] = 0.0;
	y[lo]  = 0.0;
	ey[lo] = 0.0;
      }
      delete histsum;
      histsum = new TGraphErrors(nbins,x,y,ex,ey);
      //Reset for each group
      Double_t *xp;
      xp = histsum->GetX();
      Double_t *exp;
      exp = histsum->GetEX();
      Double_t *yp;
      yp = histsum->GetY();
      Double_t *eyp;
      eyp = histsum->GetEY();
     
      //Add the first
      histsum = calibra->AddProfiles(his,histsum);
      if((Int_t) chainplace->GetEntriesFast() > 1){
	for(Int_t s = 1; s < (Int_t)chainplace->GetEntriesFast(); s++){
	  treeChain->GetEntry(((AliTRDCalibra::AliTRDPlace *)chainplace->At(s))->GetPlace());
	  histsum = calibra->AddProfiles(his,histsum);
	}
      }
      tree->Fill();
    }
    
  }


  //Write in a file
  TFile *fout = TFile::Open((const char*) namefile,"UPDATE");
  fout->WriteTObject(tree,tree->GetName(),(Option_t *) "kOverWrite");
  fout->Close();

  return kTRUE;
}
