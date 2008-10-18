/*
  TTOOLKIT to acces files on XRD 
  WORKING ONLY AT GSI
//

  Example usage:
  Load toolkit
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros")
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")

  1. 
  Retrieve the list of files : POSSIBLE ONLY FOR USER WITH LOGIN ACCESS TO THE XRD MACHINES 

  AliXRDPROOFtoolkit tool;
  tool.ListOfFiles("all.txt","/data.local2/sma/sim/v4-06-Rev-03/pp*", "AliESDs.root", kFALSE);
  tool.ListOfFiles("pp.txt","/data.local2/sma/sim/v4-06-Rev-03/pp", "AliESDs.root", kFALSE);

  tool.ListOfFiles("pp64.txt","/data.local2/sma/sim/v4-06-Rev-03/pp64", "AliESDs.root", kFALSE);
  tool.ListOfFiles("pp64m.txt","/data.local2/sma/sim/v4-06-Rev-03/pp64", "AliESDs.root", kTRUE);
  
  // check the list
  AliXRDPROOFtoolkit::FilterList("pp.txt","AliESDs.root esdTree AliESDfriends.root * Kinematics.root *",0)  
  AliXRDPROOFtoolkit::FilterList("pp.txt","AliESDs.root esdTree AliESDfriends.root * Kinematics.root *",1) 

  2. 
  
  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("pp.txt","esdTree",0,10)
  chain->Draw("fTPCnclsF");

  
  3. Process logs - ONLY priority users - with ssh acces
  AliXRDPROOFtoolkit tool;
  tool.FilterSegFault();
  TTree * treeSys  = tool.DumpSys(kTRUE)

*/


#include <TTree.h>
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


#include <AliXRDPROOFtoolkit.h>


ClassImp(AliXRDPROOFtoolkit);



//______________________________________________________________________________
AliXRDPROOFtoolkit::AliXRDPROOFtoolkit () : TObject () 
{
  //
  // Default  -GSI specific setup
  //
  for(Int_t i=255;i<300;i++)
  //for(Int_t i=259;i<260;i++)
    listeMachine.push_back(new TString(Form("lxb%d.gsi.de", i)));
  fUserGroup = gSystem->GetUserInfo();
  fUserName = fUserGroup->fUser;       
  fVerbose=1;
}

void AliXRDPROOFtoolkit::Print(Option_t* /*option*/) const{
  //
  // Print current setup of toolkit
  //
  printf("User name\t%s\n", fUserName.Data());
  printf("List of slaves\n");

  for(UInt_t i=0;i<listeMachine.size();i++)
    printf("%s\n",listeMachine[i]->Data());
}


//______________________________________________________________________________
void AliXRDPROOFtoolkit::AddMachine (const char*name)
{
  //
  // add computer to the list
  //
  listeMachine.push_back( new TString(name));
}


//______________________________________________________________________________
Bool_t AliXRDPROOFtoolkit::ListOfFiles(const char*fileName, const char*path, const char*filter, Bool_t displayMachine)
{
  //
  // Get the list of files on "registerd slaves"
  //
  // fileName - Resultinfg file with list
  // path     - starting path on slave e.g /data.local2/sma/
  // filter   - filter expression e.g. AliESDs.root
  // display machine - not used yet

  char filterXRD[100];
  char command[1000];
  //
  // FILTER expression - to get REDIRECTOR NAME
  //
  //what we want for sed
  //  sprintf(filterXRD,"sed s/\\/data.local2/root:\\/\\/gsiaf.gsi.de:1094\\//"); 
  //what we should write
  sprintf(filterXRD,"sed s/\\\\/data.local2/root:\\\\/\\\\/gsiaf.gsi.de:1094\\\\//");
  //
  // 
  //
  gSystem->Exec(Form("echo  >%s",fileName));
  for(UInt_t i=0; i<listeMachine.size(); i++){
    if (displayMachine){
      sprintf(filterXRD,"sed s/\\\\/data.local2/root:\\\\/\\\\/%s:1094\\\\//",listeMachine[i]->Data());
    }
    cout<<"Inspecting "<<listeMachine[i]->Data()<<" ..."<<endl;
    sprintf(command,"lsrun -m %s   find %s | grep %s | %s >> %s", listeMachine[i]->Data(), path, filter, filterXRD,fileName);
    printf(command);
    gSystem->Exec(command);    
  }
  return kTRUE;
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
    chain->Add(currentFile.Data());
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
  printf("Number of entries\t%d",entries);
  if (maxFiles<0) maxFiles=entries;
  if (maxFiles>entries) maxFiles=entries;
  //
  //
  //
  for (Int_t i=0; i<maxFiles; i++){
    Int_t ifile = TMath::Nint(gRandom->Rndm()*Float_t(entries));
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

//_______________________________________________________________________________


void   AliXRDPROOFtoolkit::FilterSegFault(const char *filter){
  //
  // Print the list of slaves finished with seg fault
  //  
  printf("List of the PROOF slaves  with seg fault (Last session) \n");
  for(UInt_t i=0; i<listeMachine.size(); i++){
    if (HasSegFault(listeMachine[i]->Data(), filter)){
      printf("%s\n",listeMachine[i]->Data());
    }
  }
}

//_______________________________________________________________________________

Bool_t    AliXRDPROOFtoolkit::HasSegFault(const char * machine, const char *filter){
  //
  // check if segmentation fault on the node
  //
  char command0[1000];
  char command[1000];
  char *path  =  "segfault.log";
  Long_t id,size,flags,modtime;
  
  // first check if we have a log file
  sprintf(command0, Form("lsrun -m %s find /data.local2/proof/%s/  | grep %s | grep .log >filelist.list", machine,fUserName.Data(), filter));
  if (fVerbose>1) {
    printf("%s\n",command0);
  }
  gSystem->Exec(command0);
  gSystem->GetPathInfo("filelist.list", &id, &size, &flags, &modtime);
  if (size<1) return kFALSE;
  
  //
  // now check the content
  //
  sprintf(command, Form("lsrun -m %s  cat `lsrun -m %s find /data.local2/proof/%s/  | grep %s|grep log ` | grep segmentation  >segfault.log", machine, machine,fUserName.Data(), filter));
  
  if (fVerbose) {
    printf("%s\n",command);
  }
  if (fVerbose>1) gSystem->Exec(command);
  gSystem->GetPathInfo(path, &id, &size, &flags, &modtime);
  return (size>10);
}


//--------------------------------------------------------------------

TTree *    AliXRDPROOFtoolkit::DumpSys(Bool_t verbose){
  //
  // Dump memory usage on worker nodes
  gSystem->Exec("echo > syswatch.log");
  for(UInt_t i=0; i<listeMachine.size(); i++){
    if (verbose) gSystem->Exec(Form("ssh %s cat /data.local2/proof/miranov/last-worker*/syswatch.log ",listeMachine[i]->Data()));
    gSystem->Exec(Form("ssh %s cat /data.local2/proof/miranov/last-worker*/syswatch.log >> syswatch.log",listeMachine[i]->Data()));
  } 
  
  TTree * tree = new TTree("asci","asci"); 
  tree->SetDirectory(0);
  //  tree->ReadFile("syswatch.log","hname/C:fname/C:entry/I:time/I:mem0/F:mem1/F:mem2/F:mem3/F");
  tree->ReadFile("syswatch.log","hname/C:fname/C:entry/I:time/I:mem0/F:mem1/F:mem2/F:mem3/F:cpu0/F:cpu1/F:cpu2/F:cpu3/F");
  return tree;
}

// //--------------------------------------------------------------------

TTree * AliXRDPROOFtoolkit::DumpSys2(Bool_t verbose){
  //
  // Dump memory usage on worker nodes
  //
  gSystem->Exec("echo > syswatch.log");
  for(UInt_t i=0; i<listeMachine.size(); i++){
    if (verbose) 
      gSystem->Exec(Form("lsrun -m %s cat /data.local2/proof/miranov/last-worker*/syswatch.log ",listeMachine[i]->Data()));
    gSystem->Exec(Form("lsrun -m  %s cat /data.local2/proof/miranov/last-workersession/syswatch.log >> syswatch.log",listeMachine[i]->Data()));
  } 
  
  TTree * tree = new TTree("asci","asci"); 
  tree->SetDirectory(0);
  //  tree->ReadFile("syswatch.log","hname/C:fname/C:entry/I:time/I:mem0/F:mem1/F:mem2/F:mem3/F");
  tree->ReadFile("syswatch.log","hname/C:fname/C:entry/I:time/I:mem0/F:mem1/F:mem2/F:mem3/F:cpu0/F:cpu1/F:cpu2/F:cpu3/F");
  return tree;
}




//


TTree*    AliXRDPROOFtoolkit::DumpFiles(Bool_t verbose){
  //
  //
  //
  gSystem->Exec("touch filewatch.log");
  for(UInt_t i=0; i<listeMachine.size(); i++){ 
    if (verbose) gSystem->Exec(Form("ssh %s cat /data.local2/proof/miranov/last-worker*/worker*/filewatch.log >>filewatch.log",listeMachine[i]->Data()));
    gSystem->Exec(Form("ssh %s cat /data.local2/proof/miranov/last-worker*/worker*/filewatch.log >>filewatch.log",listeMachine[i]->Data()));
  } 
 TTree * tree = new TTree("asci","asci"); 
 tree->SetDirectory(0);
 tree->ReadFile("filewatch.log","hname/C:fname/C:entry/I:status/I");
 return tree;  
}






//______________________________________________________________________________
void AliXRDPROOFtoolkit::CheckFiles (const char*fileIn, UInt_t checkLevel, const char*treeToRetrieve, const char*varexp, const char*selection)
{
  //
  // check the files
  //  
  fstream fic;
  fic.open(fileIn, ios_base::in);
  fstream focGood;
  fstream focBad;
  focGood.open(Form("%s.Good",fileIn), ios_base::out|ios_base::trunc);
  focBad.open(Form("%s.Bad",fileIn), ios_base::out|ios_base::trunc);
  //
  if(!fic.is_open()) {
    cout<<"Can't open file "<<fileIn<<endl;
    return;
  }
  //
  //
  //
  Long64_t size;
  char buffer[256];
  char * line;
  char * machine;
  Int_t level=0;
  Float_t compressionFactor;
  Int_t nkey;
  Int_t version;
  TObjString * fileName=0x0;
  TObjString * machineName=0x0;
  //  TChain chain ("check_chain");
  
  TFile fout(Form("%s.root",fileIn),"recreate");
  TTree * tree=new TTree("stats", "stats of AliXRDPROOFtoolkit::CheckFiles function");
  
  tree->Branch ("machine", "TObjString", &machineName, 256000,0);
  tree->Branch ("file", "TObjString", &fileName, 256000,0);
  tree->Branch ("level", &level, "level/I");
  tree->Branch ("size", &size, "size/L");
  tree->Branch ("nkey", &nkey, "nkey/I");
  tree->Branch ("compress", &compressionFactor, "compress/F");
  tree->Branch ("version", &version, "version/I");

  // start loop over the files
  fic.getline(buffer, sizeof(buffer));
  while (fic.good()){

    // get the machine name if present in the file
    machine=strchr(buffer, '\t');
    if(machine!=0x0){
      machine[0]='\0';
      line=machine+1;
      machine=buffer;
      machineName=new TObjString(machine);
    }else {
      machineName=new TObjString("x");
      line=buffer;
    }
    cout<<"Inspecting file :"<<line<<endl;  

 
    TTree * getTree=0x0;
    fileName=new TObjString(line);
    
    TFile * f=TFile::Open(line);
    //    chain.AddFile(line);
    level=0;
    size=-1;
    nkey=-1;
    compressionFactor=-1;
    version=-1;
    
    if (fileName->String().Contains("AliESDs.root")){
      //
      // check  friend file 
      //
      char  fnameFriend[1000];
      sprintf(fnameFriend,"%s/AliESDfriends.root",gSystem->DirName(fileName->String().Data()));
      cout<<"Inspecting file :"<<fnameFriend<<endl;  
      TFile * ffriend=TFile::Open(fnameFriend);
      if (!ffriend) level=-4;
      else{
	if (ffriend->IsZombie()) level=-3;
	if (ffriend->TestBit(TFile::kRecovered) || ffriend->TestBit(TFile::kWriteError)){
	  level=-2;
	}
	ffriend->Close();
	delete ffriend;
      }
    }
    

    if(level>=-1 && f!=0x0){
      level=1;
      size=f->GetSize();
      nkey=f->GetNkeys();
      compressionFactor=f->GetCompressionFactor();
      version=f->GetVersion();

      if(checkLevel>0 && !f->IsZombie()){
	level=2;
	if(checkLevel>1 && treeToRetrieve!="" && (getTree=(TTree*)f->Get(treeToRetrieve))!=0x0){
	  level=3;
	  Int_t tentries = getTree->GetEntries();
	  if (tentries>=0) level=4;
	  cout<<"Number of entries :"<<getTree->GetEntries()<<endl;  
	  
	  if(checkLevel>3 && varexp!="" &&tentries>0) {
	    getTree->SetBranchStatus("*",1);	    
	    try{
	      TH1F his("his","his",100,-1,1);
	      Int_t selected = getTree->Draw(Form("%s>>his",varexp), selection, "goff", 1, tentries-1);
	      cout<<"Selected\t"<<selected<<endl;
	      if(selected>-1){
		level=5;
		
		//try to remove the created histogrames ...
// 		TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp"); // 1D
// 		TGraph *graph = (TGraph*)gPad->GetPrimitive("Graph"); // 2D
// 		if(htemp!=0x0) {cout<<"removing TH1D"<<endl; delete htemp;}
// 		else if(graph!=0x0) {cout<<"removing TGraph"<<endl; delete graph;}
// 		else cout<<"nothing to remove : memory leak ..."<<endl;

	      }
	    }catch(std::bad_alloc){
	      cout<<"Warning : Draw option send std::badalloc"<<endl;
	    }
	  }
	  delete getTree;
	}
      }
      f->Close();
      delete f;
    }
    if (level>checkLevel){
      focGood<<line<<endl;
    }
    else{
      focBad<<line<<"\t"<<level<<endl;
    }

      
    tree->Fill();
    fic.getline(buffer, sizeof(buffer)); 
  }

  //now use chain to check
  //chain.Lookup(kTRUE);


  // Save the tree
  tree->Write();
  fout.Close();
  focGood.close();
  focBad.close();
  
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
    //    Bool_t isZip = currentFile.Contains("#");
    const char * dirname = gSystem->DirName(currentFile.Data());
    Int_t status = 0;
    for (Int_t i=0; i<array->GetEntries(); i+=2){
      char fname[1000];
      //if (isZip) sprintf(fname,
      sprintf(fname, "%s/%s",dirname,array->At(i)->GetName());
      printf("\nFile to be checked %s/%s",dirname,array->At(i)->GetName());
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
}





Bool_t  AliXRDPROOFtoolkit::XRDCopyDir(const char * idir, const char * files, const char *odir, Bool_t zip){
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
    sprintf(infile,"%s/%s", idir, array->At(ifile)->GetName());
    sprintf(outfile,"%s/%s", odir, array->At(ifile)->GetName());
    printf("%s - %s\n",infile, outfile);
    Bool_t result = TFile::Cp(infile,outfile); 
    succes &= result;
  }
  return succes;
}




