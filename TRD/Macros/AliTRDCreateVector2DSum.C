#if !defined( __CINT__) || defined(__MAKECINT__)

#include <vector>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TH1F.h>
#include <AliCDBManager.h>
#include <AliTRDCalibra.h>

#endif


Bool_t AliTRDCreateVector2DSum(const char* variablecali, const char* noime, const char* dire, const char* namefile){
  //
  // After having simulated and reconstructed events in subrepertories 0%d of dire
  // this macro searchs in the subdirectories the file TRD.calibration.root
  // takes the vectors and merge them
  // variablecali can be treeCH2d, treePH2d or treePRF2d
  // noime is the name of the "sum" tree
  // namefile is the name of the file where noime will be written
  //


  // The tree
  TChain *treeChain = new TChain(variablecali);
   

  //Variables
  std::vector<Int_t> vectorplace;
  std::vector<std::vector<Int_t> > whereinthechain;


  //TRDCalibra
  AliCDBManager *man = AliCDBManager::Instance();
  man->GetStorage("local://$ALICE_ROOT");
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
      TFile *file = (TFile *) TFile::Open(fullname,"READ");
      if(!file) continue;
      TTree *treecurrent = (TTree *) file->Get(variablecali);
      if(!treecurrent) continue;
      std::vector<Int_t> vectorcurrent = calibra->ConvertTreeVector(treecurrent);
      printf("Size of the current tree: %d\n",(Int_t) vectorcurrent.size());
      printf("Size of whereinthechain: %d\n",(Int_t) whereinthechain.size());
      printf("Size of the chain: %d\n", (Int_t) treeChain->GetEntries());
      Int_t j = (Int_t) treeChain->GetEntries();
	for(Int_t jui = 0; jui < (Int_t) vectorcurrent.size(); jui++){
	  //Search if already found
	  Int_t place = calibra->SearchInTreeVector(vectorplace,vectorcurrent[jui]);
	  //Create a new element in the two std vectors
	  if(place == -1){
	    std::vector<Int_t> chainplace;
	    chainplace.push_back((Int_t) (j+jui));
	    for(Int_t i = 0; i < (Int_t) chainplace.size(); i++){
	      //cout << "i: " << i << "value: " << chainplace[i] << endl;
	    }
	    vectorplace.push_back((Int_t) (vectorcurrent[jui]));
	    whereinthechain.push_back(chainplace);
	   
	  }
	  //Update the element at the place "place" in the std vector whereinthechain
	  else {
	    std::vector<Int_t> chainplace = whereinthechain[place];
	    chainplace.push_back((Int_t) (j+jui));
	    //cout << "size of chainplace place > -1: " << (Int_t) chainplace.size() << endl;
	    for(Int_t i = 0; i < (Int_t) chainplace.size(); i++){
	      //cout << "i: " << i << "value: " << chainplace[i] << endl;
	    }
	    std::vector<std::vector<Int_t> >::iterator it = whereinthechain.begin()+place;
	    whereinthechain.erase(it);
	    whereinthechain.insert(it,chainplace);
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
    
    printf("FINAL Size of the chain: %d\n", (Int_t) treeChain->GetEntries());
    printf("FINAL Size of vectorplace: %d\n", (Int_t) vectorplace.size());
    for(Int_t h = 0; h < (Int_t) vectorplace.size(); h++){
      group = vectorplace[h];
      std::vector<Int_t> chainplace = whereinthechain[h];
      treeChain->GetEntry(chainplace[0]);
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
      if((Int_t) chainplace.size() > 1){
	for(Int_t s = 1; s < (Int_t)chainplace.size(); s++){
	  treeChain->GetEntry(chainplace[s]);
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
    
    printf("FINAL Size of the chain: %d\n", (Int_t) treeChain->GetEntries());
    printf("FINAL Size of vectorplace: %d\n", (Int_t) vectorplace.size());
    for(Int_t h = 0; h < (Int_t) vectorplace.size(); h++){
      group = vectorplace[h];
      std::vector<Int_t> chainplace = whereinthechain[h];
      treeChain->GetEntry(chainplace[0]);
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
      if((Int_t) chainplace.size() > 1){
	for(Int_t s = 1; s < (Int_t)chainplace.size(); s++){
	  treeChain->GetEntry(chainplace[s]);
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
