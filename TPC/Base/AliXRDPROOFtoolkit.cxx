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

/*
  marian.ivanov@cern.ch
  
  //TOOLKIT for chain manipulation:  
  //Example usage:
  //
  // 1. Check the list of files
  AliXRDPROOFtoolkit toolkit;
  AliXRDPROOFtoolkit::FilterList("pp.txt","AliESDs.root esdTree AliESDfriends.root * Kinematics.root *",0)
  AliXRDPROOFtoolkit::FilterList("pp.txt","AliESDs.root esdTree AliESDfriends.root * Kinematics.root *",1)
  //
  //2. make a chain or random chain form the list of files
  TChain * chain = toolkit.MakeChain("esd.txt","esdTree",0,10)
  TChain * chainRandom = toolkit.MakeChainrandom("esd.txt","esdTree",0,10)
  chain->Draw("fTPCnclsF");
  

*/
#include <TTree.h>
#include <TEnv.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TDSet.h>
#include <TH1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TPad.h>
#include <exception>
#include <fstream>
#include <TRandom.h>
#include <TTimeStamp.h>
#include <AliXRDPROOFtoolkit.h>


ClassImp(AliXRDPROOFtoolkit)



//______________________________________________________________________________
  AliXRDPROOFtoolkit::AliXRDPROOFtoolkit () : 
    TObject () ,
    fVerbose(kFALSE),          // verbso mode  - print command 
    fUserName(""),         // user name
    fUserGroup(0)        // user group info
{
  //
  // 
  //
  fUserGroup = gSystem->GetUserInfo();
  fUserName  = fUserGroup->fUser;       
  fVerbose=1;
}





TChain* AliXRDPROOFtoolkit::MakeChain(const char*fileIn, const char * treeName, const char *fName, Int_t maxFiles, Int_t startFile)
{
  //
  // Create the chain
  //
  TChain* chain = new TChain(treeName);

  // Open the input stream
  ifstream in;
  in.open(fileIn);

  // Read the input list of files and add them to the chain
  TString currentFile;
  Int_t counter=0;
  while(in.good()) {
    in >> currentFile;
    if (fName) {
      currentFile+="#";
      currentFile+=fName;
    }
    if (!currentFile.Contains("root")) continue; // protection
    counter++;
    if (counter<startFile) continue;
    if (counter>maxFiles+startFile) break;
    TFile * f = TFile::Open(currentFile.Data());
    if (f){
      chain->Add(currentFile.Data());
    }
    delete f;
  }

  in.close();

  return chain;
}

TChain* AliXRDPROOFtoolkit::MakeChainRandom(const char*fileIn, const char * treeName,const char *fName, Int_t maxFiles, Int_t startFile)
{
  //
  // Create a TDSet - files are in random order
  //
  // filein    - input list text file
  // treename  - containg tree 
  // maxFiles  - maximum number of files included

  TObjArray array(10000);
  
  TChain* chain = new TChain(treeName);

  // Open the input stream
  ifstream in;
  in.open(fileIn);

  // Read the input list of files and add them to the chain
  TString currentFile;
  Int_t counter=0;
  while(in.good()) {
    in >> currentFile;
    if (fName) {
      currentFile+="#";
      currentFile+=fName;
    }
    if (!currentFile.Contains("root")) continue; // protection
    counter++;
    //    chain->Add(currentFile.Data());
    array.AddLast(new TObjString(currentFile));
  }
  in.close();
  Int_t entries = array.GetEntries();
  printf("Number of entries\t%d\n",entries);
  //
  //
  //
  Double_t *randomI = new Double_t[entries];
  Int_t *indexes = new Int_t[entries];
  for (Int_t i=0;i<entries; i++) randomI[i]=gRandom->Rndm();
  TMath::Sort(entries,randomI,indexes); 
  
  for (Int_t i=startFile; (i<startFile+maxFiles) && (i<entries); i++){
    Int_t ifile = indexes[i];
    if (ifile<entries && (array.At(ifile)) &&  array.At(ifile)->TestBit(TObject::kCannotPick)==kFALSE){ 
      printf("%d\t%d\t%s\n",i, ifile, array.At(ifile)->GetName());
      chain->Add(array.At(ifile)->GetName());
      array.At(ifile)->SetBit(TObject::kCannotPick);
    }
  }
  return chain;
}



TDSet* AliXRDPROOFtoolkit::MakeSet(const char*fileIn, const char * treeName, const char *fName, Int_t maxFiles)
{
  //
  // Create the TDSet out of list
  // filein    - input list text file
  // treename  - containg tree 
  // maxFiles  - maximum number of files included

  TDSet* chain = new TDSet(treeName);

  // Open the input stream
  ifstream in;
  in.open(fileIn);

  // Read the input list of files and add them to the chain
  TString currentFile;
  Int_t counter=0;
  while(in.good()) {
    in >> currentFile;
    if (fName) {
      currentFile+="#";
      currentFile+=fName;
    }
    if (!currentFile.Contains("root")) continue; // protection
    counter++;
    if (maxFiles>0 && counter>maxFiles) break;
    chain->Add(currentFile.Data());
  }

  in.close();
  chain->Validate();
  return chain;
}


TDSet* AliXRDPROOFtoolkit::MakeSetRandom(const char*fileIn, const char * treeName, const char *fName, Int_t maxFiles)
{
  //
  // Create a TDSet - files are in random order
  //
  // filein    - input list text file
  // treename  - containg tree 
  // maxFiles  - maximum number of files included

  TObjArray array(10000);
  
  TDSet* chain = new TDSet(treeName);

  // Open the input stream
  ifstream in;
  in.open(fileIn);

  // Read the input list of files and add them to the chain
  TString currentFile;
  Int_t counter=0;
  while(in.good()) {
    in >> currentFile;
    if (fName) {
      currentFile+="#";
      currentFile+=fName;
    }
    if (!currentFile.Contains("root")) continue; // protection
    counter++;
    //    chain->Add(currentFile.Data());
    array.AddLast(new TObjString(currentFile));
  }
  in.close();
  Int_t entries = array.GetEntries();
  printf("Number of entries\t%d",entries);
  if (maxFiles<0) maxFiles=entries;
  if (maxFiles>entries) maxFiles=entries;
  for (Int_t i=0; i<maxFiles; i++){
    Int_t ifile = TMath::Nint(gRandom->Rndm()*Float_t(entries));
    if (ifile<entries && (array.At(ifile)) &&  array.At(ifile)->TestBit(TObject::kCannotPick)==kFALSE){ 
      printf("%d\t%d\t%s\n",i, ifile, array.At(ifile)->GetName());
      chain->Add(array.At(ifile)->GetName());
      array.At(ifile)->SetBit(TObject::kCannotPick);
    }
  }


  chain->Validate();
  return chain;
}







Int_t  AliXRDPROOFtoolkit::CheckTreeInFile(const char*fileName,const char*treeName, Int_t debugLevel, const char *branchName){
  //
  // Check the tree in file 
  // fileName   - the name of the file with tree
  // treeName   - the name of file
  // debugLevel - 0 check the existance of the file -  1 make loop over entries
  // branchName - if debugLevel>0 the branch is chcecked
  //              if brnachName =0 the content of full tree is chcecked
  // return value = 0 - Check things  OK
  //               -1 - file not exist or not accesible
  //               -2 - file is zombie
  //		   -3 - tree not present
  //               -4 - branch not present
  TFile * file = TFile::Open(fileName);
  if (!file) { return -1;}
  if (file->IsZombie()) {file->Close(); delete file; return -2;};

  TString TrName(treeName);
  if (TrName=="*") {
    //cout <<"        treename ==== *"<<endl;;
    file->Close(); delete file; 
    return 0;
  }
  TTree * tree = (TTree*)file->Get(treeName);
  if (!tree) {file->Close(); delete file; return -3;}
  TBranch * branch = 0;
  if (branchName) {
    branch = tree->GetBranch(branchName);
    if (!branch) {file->Close(); delete file; return -4;}
  }
  //
  if (debugLevel==1 &&  tree->GetEntries()==0 ) return 1; //empty 

  tree->SetBranchStatus("*",1);
  try {
    if (debugLevel>1){
      Int_t entries = tree->GetEntries();
      for (Int_t i=0;i<entries; i++){
	if (branch) branch->GetEntry(i);
	else tree->GetEntry();      
      }
    }
  }catch ( ... ) {
    printf("PROBLEM\n");  
    // never catched  - as there is no exception in the ROOT IO
    file->Close(); delete file;
    return 1 ;
  }

  file->Close(); delete file;
  return 0;
}


Bool_t  AliXRDPROOFtoolkit::FilterList(const char*inputList, const char*fileList, Int_t checkLevel){
  //
  // Filter the list  
  // inputList - list of original file names
  // fileList  - list of file to be checked
  //           - 0 - fileName
  //           - 1 - treeName (if * not checked)
  //           - 2 - fileName 
  //                 ....
  // checkLevel - 0 - check only existance of the files and tree's + 
  //                  simple file corruption
  //            > 1 - check the content of the tree - 
  //                  (can crash as there do not exest exception handling in ROOT)
  // Output -  two streams are created - file with good entries
  // "%s.Good a,d file with bad entries %s.Bad
  //EXAMPLE:
  // AliXRDPROOFtoolkit::FilterList("ppgrid2.txt","AliESDs.root esdTree AliESDfriends.root * Kinematics.root *",1) 
  gEnv->SetValue("TFile.Recover", 0);
  //
  fstream finput;
  finput.open(inputList, ios_base::in);
  fstream focGood;
  fstream focBad;
  focGood.open(Form("%s.Good",inputList), ios_base::out|ios_base::trunc);
  focBad.open(Form("%s.Bad",inputList), ios_base::out|ios_base::trunc);
  //
  if(!finput.is_open()) {
    cout<<"Can't open file "<<inputList<<endl;
    return kFALSE;
  }
  //
  // Read the input list of files and add them to the chain
  //
  TObjArray *array = (TString(fileList)).Tokenize(" ");
  TString currentFile;
  Int_t counter=0;
  while(finput.good()) {
    finput >> currentFile;
    if (!currentFile.Contains("root")) continue; // protection
    if (currentFile.Contains("alien://")){
      focGood<<currentFile<<endl;
      continue;
    }
    Bool_t isZip = currentFile.Contains("#");
    const char * dirname = gSystem->DirName(currentFile.Data());
    Int_t status = 0;
    //
    for (Int_t i=0; i<array->GetEntries(); i+=2){
      char fname[1000];
      if (!isZip){
	snprintf(fname,1000, "%s/%s",dirname,array->At(i)->GetName());
        if (((TObjString*)array->At(i))->String().Contains("*")){
	  snprintf(fname,1000, "%s", currentFile.Data());
	}
      }
      if (isZip) {
	const char * fileName   =  gSystem->BaseName(currentFile.Data());
	TString fstring=fileName;
	fstring[fstring.First("#")]=0;
	snprintf(fname,1000, "%s/%s#%s",dirname,fstring.Data(),array->At(i)->GetName());
	printf(fname, "To check %s%s#%s\n",dirname,fstring.Data(),array->At(i)->GetName());
      }

      printf("\nFile to be checked %s\n",fname);
      //cout <<"\n arguments: "<< array->At(i+1)->GetName()<<" "<<checkLevel<<endl;
      Int_t cstatus = CheckTreeInFile(fname, array->At(i+1)->GetName(), checkLevel,0);
      //printf("  CheckTreeInFile returns %d",cstatus);
      if (cstatus!=0) {
	status = cstatus; 
	break;
      }
    }
    if (status==0){
      focGood<<currentFile<<endl;
    }else{
      focBad<<currentFile<<endl;
    }
    counter++;    
  }
  finput.close();
  delete array;
  return kTRUE;
}


Bool_t  AliXRDPROOFtoolkit::FilterListZip(const char*inputList, const char*fileList, Int_t checkLevel){
  //
  // Filter the list  
  // inputList - list of original file names
  // fileList  - list of file to be checked
  //           - 0 - fileName
  //           - 1 - treeName (if * not checked)
  //           - 2 - fileName 
  //                 ....
  // checkLevel - 0 - check only existance of the files and tree's + 
  //                  simple file corruption
  //            > 1 - check the content of the tree - 
  //                  (can crash as there do not exest exception handling in ROOT)
  // Output -  two streams are created - file with good entries
  // "%s.Good a,d file with bad entries %s.Bad
  //EXAMPLE:
  // AliXRDPROOFtoolkit::FilterList("ppgrid2.txt","AliESDs.root esdTree AliESDfriends.root * Kinematics.root *",1) 

  fstream finput;
  finput.open(inputList, ios_base::in);
  fstream focGood;
  fstream focBad;
  focGood.open(Form("%s.Good",inputList), ios_base::out|ios_base::trunc);
  focBad.open(Form("%s.Bad",inputList), ios_base::out|ios_base::trunc);
  //
  if(!finput.is_open()) {
    cout<<"Can't open file "<<inputList<<endl;
    return kFALSE;
  }
  //
  // Read the input list of files and add them to the chain
  //
  TObjArray *array = (TString(fileList)).Tokenize(" ");
  TString currentFile;
  Int_t counter=0;
  while(finput.good()) {
    finput >> currentFile;
    if (!currentFile.Contains("root")) continue; // protection
    if (currentFile.Contains("alien://")){
      focGood<<currentFile<<endl;
      continue;
    }
    //Bool_t isZip = currentFile.Contains("#");
    const char * dirname = gSystem->DirName(currentFile.Data());
    const char * fileName   =  gSystem->BaseName(currentFile.Data());
    TString fstring=fileName;
    fstring[fstring.First("#")]=0;
    Int_t status = 0;
    for (Int_t i=0; i<array->GetEntries(); i+=2){
      char fname[1000];
      //if (isZip) sprintf(fname,
      snprintf(fname,1000, "%s/%s#%s",dirname,fstring.Data(),array->At(i)->GetName());
      printf(fname, "To check %s%s#%s\n",dirname,fstring.Data(),array->At(i)->GetName());
      //cout <<"\n arguments: "<< array->At(i+1)->GetName()<<" "<<checkLevel<<endl;
      Int_t cstatus = CheckTreeInFile(fname, array->At(i+1)->GetName(), checkLevel,0);
      //printf("  CheckTreeInFile returns %d",cstatus);
      if (cstatus!=0) {
	status = cstatus; 
	break;
      }
    }
    if (status==0){
      focGood<<currentFile<<endl;
    }else{
      focBad<<currentFile<<endl;
    }
    counter++;    
  }
  finput.close();
  return kTRUE;
}





Bool_t  AliXRDPROOFtoolkit::XRDCopyDir(const char * idir, const char * files, const char *odir, Bool_t /*zip*/){
  //
  // idir  - input directory
  // odir  - output directory
  // files - the list of files to be coppied
  // zip   - not supported yet
  //
  // Example :									
  //
  // idir ="root://gsiaf.gsi.de:1094//sma/sim/v4-05-Rev-03/pp/0000";
  // odir ="root://lxgrid2.gsi.de:1094//miranov/test/pp/0000"; 
  // char *files="AliESDs.root AliESDfriend.root Kinematics.root";
  TString str(files);
  TObjArray * array = str.Tokenize(" "); 
  Int_t nfiles = array->GetEntries();
  char infile[1000];
  char outfile[1000];
  Bool_t succes=kTRUE;
  for (Int_t ifile =0; ifile<nfiles; ifile++){
    snprintf(infile,1000,"%s/%s", idir, array->At(ifile)->GetName());
    snprintf(outfile,1000,"%s/%s", odir, array->At(ifile)->GetName());
    printf("%s - %s\n",infile, outfile);
    Bool_t result = TFile::Cp(infile,outfile); 
    succes &= result;
  }
  delete array;
  return succes;
}



void AliXRDPROOFtoolkit::JoinTreesIndex(const char * outputFile, const char * outputTree, const char *indexName, const char *inputTrees, Int_t debugLevel){
  //
  // Join together several tree according to the index
  // 
  // Parameters:
  // Output:
  //     outputFile : name of the output file
  //     outputTree : name of the output Tree
  //     indexName  : name of the branch to be used as an index
  // Input:
  //     inputTrees : decription of the input trees setup
  /*
    Example usage:
    
    AliXRDPROOFtoolkit::JoinTreesIndex("outAll.root","joinAll","run","1#CPass1#run#tpcQA#TPCCPass1.root+1#VPass1#run#tpcQA#TPCVPass1.root+1#Pass1#run#tpcQA#TPCPass1.root+0#DAQ#run#joinTree#fproductionJoin.root+0#C#run#dcs#OCDBscan.root+0#CE#run#Fits#CEtrend.root");
    ==>
    Combine information form the Cpass1,VPass, and Pass1QA, calibration tree, DAQ information, trigger information 
    Make a File "outAll.root",  with tree "joinAll", index of tree with name "run"
    //
    // Input tree configuration string:
    //
    const char *inputTrees="1#CPass1#run#tpcQA#TPCCPass1.root+1#VPass1#run#tpcQA#TPCVPass1.root+1#Pass1#run#tpcQA#TPCPass1.root+0#DAQ#run#joinTree#/home/miranov/test/dbQueries/fproductionJoin.root+0#C#run#dcs#OCDBscan.root+0#CE#run#Fits#CEtrend.root"
    Describe 6 trees to be merged (string separated be +):
      TObjArray *arrayInput = TString(inputTrees).Tokenize("+");
      TObjString = 1#CPass1#run#tpcQA#TPCCPass1.root
      TObjString = 1#VPass1#run#tpcQA#TPCVPass1.root
      TObjString = 1#Pass1#run#tpcQA#TPCPass1.root
      TObjString = 0#DAQ#run#joinTree#/home/miranov/test/dbQueries/fproductionJoin.root
      TObjString = 0#C#run#dcs#OCDBscan.root
      TObjString = 0#CE#run#Fits#CEtrend.root
    //  
    Each tree is characterize by 5 parameters - separate by #
       description="1#CPass1#run#tpcQA#TPCCPass1.root"
       TString(description)->Tokenize("#").Print()
       Collection name='TObjArray', class='TObjArray', size=16
       TObjString = 1                    ==> (0/1) index is used 
       TObjString = CPass1               ==> name of output branch  in output tree
       TObjString = run                  ==> name of the index
       TObjString = tpcQA                ==> name of the input tree in the input file
       TObjString = TPCCPass1.root       ==> name of the input file
  */
  //
  //
  //
                
  TFile * fout = new TFile(outputFile,"recreate");
  fout->cd();
  TTree *joinTree=new TTree(outputTree,outputTree);
  //
  // 1. Define setup. parse definition string
  //
  TObjArray *arrayInput = TString(inputTrees).Tokenize("+");
  Int_t nTrees = arrayInput->GetEntries();
  TObjArray * arrayFile  = new TObjArray(nTrees);    // array of TFiles with trees
  TObjArray * arrayTrees = new TObjArray(nTrees);    // array of trees 
  TObjArray * arrayNames = new TObjArray(nTrees);    // name of tree
  TObjArray * arrayRunID = new TObjArray(nTrees);    // name of tree
  TArrayI arrayEnableTree(nTrees);  
  for (Int_t i=0; i<2; i++) printf("\n");
  printf("Joing query\n");
  arrayInput->Print();
  for (Int_t i=0; i<2; i++) printf("\n");
  {for (Int_t itree=0; itree<nTrees; itree++){
      //
      TObjArray *description = TString(arrayInput->At(itree)->GetName()).Tokenize("#");
      if (description->GetEntries()<4) {
	printf("Fatal: Invalid description:  %s\n", arrayInput->At(itree)->GetName());
	continue;
      }
      TFile * f = TFile::Open(description->At(4)->GetName());
      if (!f){
	printf("Fatal: Invalid description: fileName %s\n", description->At(4)->GetName());
	delete arrayInput;
	return;
      }
      arrayFile->AddAt(f,itree);
      TTree * tree = (TTree*)f->Get(description->At(3)->GetName());
      if (!tree){
	printf("Fatal: Invalid description. Tree name\t%s\n", description->At(3)->GetName());
	delete arrayInput;
	return;
      }
      tree->SetCacheSize(400000000);
      //    
      arrayTrees->AddAt(tree,itree);
      //
      arrayRunID->AddAt(new TObjString(description->At(2)->GetName()),itree);
      arrayNames->AddAt(new TObjString(description->At(1)->GetName()),itree);
      arrayEnableTree[itree]=atoi(description->At(0)->GetName());    

    }}
  //  
  delete arrayInput;
  // 2. Make the run list
  //
  //
  map<int, int> runMap;
  map<int, int> *runMapTree = new map<int, int>[nTrees];
  //map<int, int> runMapTree[nTrees];
  {for (Int_t itree=0; itree<nTrees; itree++){
      TTree * tree = (TTree*)arrayTrees->At(itree);
      Int_t entries=tree->GetEntries();
      char query[2000];
      snprintf(query,2000,"%s:Entry$", arrayRunID->At(itree)->GetName());
      entries = tree->Draw(query,"","goff");      
      for (Int_t ientry=0;ientry<entries; ientry++){
	Int_t irun=Int_t(tree->GetV1()[ientry]);
	//	Int_t entryNr=Int_t(tree->GetV2()[ientry]);
 	if (arrayEnableTree[itree]>0) runMap[irun]+=1;
 	runMapTree[itree][irun]=ientry;
 	if (debugLevel>0) printf("%s\t%d\t%d\n",tree->GetName(), irun, 	runMapTree[itree][irun]);
      }
    }
  }
  //
  // 3. Make join tree
  //
  Int_t jrun=0;
  fout->cd();
  joinTree->Branch(indexName, &jrun,Form("%s/I",indexName));
  Int_t *status=new Int_t[nTrees];
  char *brName = new char[10000];
  char *brTitle= new char[10000];
  //
  
  {for (Int_t itree=0; itree<nTrees; itree++){
      TTree * tree = (TTree*)arrayTrees->At(itree);      
      tree->GetEntry(1);
      TString treeName=arrayNames->At(itree)->GetName();
      if (treeName.Length()>0){
	joinTree->Branch(Form("%s.status",treeName.Data()), &status[itree],Form("%s.status/I",treeName.Data()));
      }else{
	joinTree->Branch("status", &status[itree],"status/I");
      }
      //
      Int_t nbranches= tree->GetListOfBranches()->GetEntries();
      for (Int_t ibr=0; ibr<nbranches; ibr++){
	TBranch * br = (TBranch*)(tree->GetListOfBranches()->At(ibr));
	if (treeName.Length()>0){
	  sprintf(brName,"%s.%s",treeName.Data(), br->GetName());
	  sprintf(brTitle,"%s.%s",treeName.Data(), br->GetTitle());
	}else{
	  sprintf(brName,"%s",br->GetName());
	  sprintf(brTitle,"%s",br->GetTitle());
	}
	void* addr = 0;
	TString className=br->GetClassName();
	if (className.Length()==0){
	  TString str(br->GetTitle());
	  if (str[str.Length()-1]=='I') addr=new Int_t;
	  if (str[str.Length()-1]=='F') addr=new Float_t;
	  if (str[str.Length()-1]=='D') addr=new Double_t;
	  if (str[str.Length()-1]=='C') addr=new Char_t[10000];
	  if (addr) joinTree->Branch(brName, addr, brTitle);
	  br->SetAddress(addr);
	}else{
	  TClass cclass(className);
	  TObject **addrClass =  new TObject *;
	  (*addrClass)=0;
	  printf("%s\t%s\n",br->GetName(), className.Data());
	  br->SetAddress(addrClass);	  
	  br->GetEntry(0);	  
	  joinTree->Branch(brName,addrClass);	  	  
	}	
      }
    }
  }
  joinTree->Write();
  //
  // 4. Fill the trees
  //
  map<int, int>::iterator riter;   
  {for (riter=runMap.begin(); riter != runMap.end(); ++riter){
      printf("%d\t%d\t", riter->first, riter->second);
      jrun=riter->first;
      for (Int_t itree=0; itree<nTrees; itree++){
	TTree * tree = (TTree*)arrayTrees->At(itree); 
	Int_t entry= runMapTree[itree][jrun];
	status[itree]=(entry>0)?1:0;
	if (entry>=0) tree->GetEntry(entry);
	printf("%d\t",entry);
	//
      }
      joinTree->Fill();
      printf("\n");
    }}
  fout->cd();
  joinTree->Write(outputTree);
  fout->Close();

}



void AliXRDPROOFtoolkit::CacheFileList(const char * fileIn, const char* cachePrefix){
  //
  // cache the list of file locally, cache valeus are stored in the location
  // specified by optional argumen prefix 
  // 2 new files are created 
  //       <fileIn>.cache    - file with the location of cahe files
  //       <fileIn>.cacheLog - log file +list of files which can not be cached 
  //      
  /*
    fileIn = "TPCCPass1.list";
    cachePrefix = "";
  */
  ifstream fin;
  fin.open(fileIn);  
  ofstream fout;
  fout.open(Form("%s.cache",fileIn));  
  ofstream foutLog;
  foutLog.open(Form("%s.cacheLog",fileIn));  
  // Read the input list of files and add them to the chain
  TString currentFile;
  TString cacheFile;
  Int_t counter=0;
  {while(fin.good()) {
      TTimeStamp s;
      TString fname;
      fin >> currentFile;
      fname=currentFile;
      fname.ReplaceAll("-","_");
      fname.ReplaceAll("/","_");
      fname.ReplaceAll(":","_");
      cacheFile=cachePrefix;
      cacheFile+=fname;
      printf("%s\t%s\n",currentFile.Data(),cacheFile.Data());
      if (TFile::Cp(currentFile.Data(),cacheFile.Data())){
	fout<<cacheFile.Data()<<"\n";
	foutLog<<s.AsString();
	foutLog<<cacheFile.Data()<<"n";
      }else{
	foutLog<<"Copy failed"<<currentFile.Data()<<cacheFile.Data()<<"\n";
      }
    }}
  fout.close();
  foutLog.close();
}



void   AliXRDPROOFtoolkit::MakeTreeFromList(const char *fout, const char * treeOut, const char * treeIn, const char * flist, Bool_t debug){
  //
  // join trees from the list  and make a common tree - stored in the file 
  // 
  /*
    Example:
    const char * fout="TPCCpass1.root";
    const char *treeOut="tpcQA"
    const char *treeIn="tpcQA"
    const char * flist="TPCCPass1.list"    
  */
  if (debug>0){
    printf("MakeTreeFromList\n");
    printf("fout=%s\n",fout);
    printf("treeOut=%s\n",treeOut);
    printf("treeIn=%s\n",treeIn);
    printf("fileList=%s\n",flist);
  } 
  ifstream fin;
  fin.open(flist);  
  ofstream foutLog;
  foutLog.open(Form("%s.chainLog",flist));  
  // Read the input list of files and add them to the chain
  TString currentFile;
  Int_t counter=0;
  Int_t nbranches=0;
  {while(fin.good()) {
      fin >> currentFile;
      TFile * f = TFile::Open(currentFile.Data());
      foutLog<<"Opening file"<<currentFile.Data();
      if (!f) {
	foutLog<<"Error opening file\t"<<currentFile<<"\n";
	cout<<"Error opening file\t"<<currentFile<<"\n";
	continue;
      }
      TTree * tree = (TTree*)f->Get(treeIn);
      if (!tree) {
	foutLog<<"Error opening tree\t"<<currentFile<<treeIn<<"\n";
	cout<<"Error opening tree\t"<<currentFile<<treeIn<<"\n";
	f->ls();	
	continue;
      }
      if (tree->GetListOfBranches()==0){
	foutLog<<"Error opening tree\t"<<currentFile<<treeIn<<"\n";
	cout<<"Error opening tree\t"<<currentFile<<treeIn<<"\n";
	continue;
      }
      Int_t nbranchesCurrent = tree->GetListOfBranches()->GetEntries();
      if ( nbranches ==0 ) nbranches=nbranchesCurrent;
      if ( nbranches!=nbranchesCurrent){
	foutLog<<"Error  tree layout\t"<<currentFile<<treeIn<<nbranches<<nbranchesCurrent<<"\n";
	cout<<"Error tree layout\t"   <<currentFile<<treeIn<<nbranches<<nbranchesCurrent<<"\n";
      }     
      counter++;
    }
  }
  foutLog<<"Number of files"<<counter<<"\n";
  cout<<   "Number of files"<<counter<<"\n";
  //

  TChain * chain = AliXRDPROOFtoolkit::MakeChain(flist,treeIn,0,1000000000,0);
  Bool_t status=kTRUE;
  if (!chain) status=kFALSE;
  if (chain->GetEntries()==0) status=kFALSE;
  if (!status){
    printf("Incorrect list (%s) or trees (%s)", flist,treeIn);
    return;
  }
  TFile *fileOut= TFile::Open(fout, "recreate");
  TTree * tree = chain->CopyTree("1");
  fileOut->cd();
  tree->Write(treeOut);
  fileOut->Close(); 
}
