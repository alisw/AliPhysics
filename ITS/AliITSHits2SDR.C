#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TFile.h>

#include "AliRun.h"
#include "AliHeader.h"
#include "AliITS.h"
#include "AliITSDetType.h"
#include "AliITSsimulationFastPoints.h"
#include "AliITSresponse.h"
#include "AliITSresponseSDD.h"
#include "AliITSreconstruction.h"
#include "AliRunDigitizer.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSDigitizer.h"

#endif

TFile* AccessInpFile(TString inFile);
void writeAR(TFile * fin, TFile *fou);
void AliITSSD2D(TString inFile = "galice.root", TString outFile = "");
void AliITSH2FR2files (TString inFile, TString outfile);
void AliITSH2SD2Files(TFile *file2);
void copyTK(TString inFile, TString outFile);

void AliITSHits2SDR(TString inFile = "galice.root", TString File1 = "galiceS.root",TString File2 = "galiceD.root", TString File3 = "galiceR.root"){
  ////////////////////////////////////////////////////////////////////
  // This macro takes hits from inFile file and 
  // produce fast RecPoints, Sdigits, Digits and slow Rec Points 
  // for the ITS using the standard detector simulations. 
  // The first argument is the input file which contains the tree with hits
  // The second argument is the output file name onto which TreeS will  
  // be written. 
  // The third and fourth arguments are the file names onto which TreeD  
  // and TreeR+TreeC will be written (respectively).
  // The AliRun object and the KINE tree will be saved onto 
  // these files as well. 
  // This macro will process all of the events on input the root file.
  // This macro can be compiled and it should be launched from an aliroot 
  // session.
  // WARNING: Fast Rec Points are written on a branch named ITSRecPointsF
  // this macro should be used with aliroot (it does not load the libraries) 
  ////////////////////////////////////////////////////////////////////
 
 
 
  // Produce Fast Rec Points  (on File3)
  cout<<"Producing Fast rec points...\n";
  AliITSH2FR2files(inFile,File3);

  // TreeK is copied to File3 (for comparison purposes)
 
  copyTK(inFile,File3);
 
  // Produce SDigits
  TFile *file = AccessInpFile(inFile);
  if(!file){
    cerr<<"Problems in accessing the input file\n";
    return;
  }
  TFile *file2 = gAlice->InitTreeFile("S",File1);
  file2->cd();
  cout<<"Producing Summable digits ....\n";
  AliITSH2SD2Files(file2);
  // write the AliRun object to the output file
  writeAR(file,file2);
  delete gAlice;
  gAlice = 0;
  file->Close();
  delete file;
  file2 = 0;
  
  // Produce Digits 
  cout<<"Producing Digits ...\n";
  AliITSSD2D(File1,File2);

  // Produce Slow Rec. Points 
  TFile *fin = new TFile(File2);
  gAlice = (AliRun*) fin->Get("gAlice");
  const char *nulptr=0;
  AliITSreconstruction *itsr = new AliITSreconstruction(nulptr);
  itsr->SetOutputFile(File3);
  itsr->Init();
  itsr->Exec();
  TFile *fou = gAlice->GetTreeRFile();
  writeAR(fin,fou);
  delete itsr;
  
}

void AliITSH2SD2Files(TFile *file2){

  // Function used to produce s-digits
 
  AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
  if (!ITS) {
	cerr<<"AliITS object not found on file" << endl;
	return;
  }
  if(!(ITS->GetITSgeom())){
	cerr << " AliITSgeom not found. Can't digitize without it." << endl;
	return ;
  }

  // for old files
  AliITSresponseSDD *resp1 = (AliITSresponseSDD*)ITS->DetType(1)->GetResponseModel();
  TDatime *ct0= new TDatime(2002,04,26,00,00,00),ct=file2->GetCreationDate();
  if(ct0->GetDate()<ct.GetDate()){
    resp1 = new AliITSresponseSDD();
    ITS->SetResponseModel(1,resp1);
    cout << "Changed response class for SDD: \n";
    resp1->Print();
  } // end if
  // end mods for old files (<26/4/2002)
  for (Int_t nevent=0; nevent<gAlice->TreeE()->GetEntries(); nevent++) {
    gAlice->GetEvent(nevent);
    gAlice->MakeTree("S",file2);
    ITS->MakeBranch("S");
    ITS->SetTreeAddress();   
    ITS->Hits2SDigits();
  }
}


TFile * AccessInpFile(TString inFile){

  // Function used to open the input file and fetch the AliRun object

  TFile *retfil = 0;
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
  if (file) {file->Close(); delete file;}
  file = new TFile(inFile);
  if (!file->IsOpen()) {
	cerr<<"Can't open "<<inFile<<" !" << endl;
	return retfil;
  } 

  // Get AliRun object from file or return if not on file
  if (gAlice) delete gAlice; gAlice = 0;
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
	cerr << "AliRun object not found on file"<< endl;
	return retfil;
  } 
  return file;
}

void writeAR(TFile * fin, TFile *fou) {
  TDirectory *current = gDirectory;
  TTree *Te;
  TTree *TeNew;
  AliHeader *alhe = new AliHeader();
  Te = (TTree*)fin->Get("TE");
  Te->SetBranchAddress("Header",&alhe);
  Te->SetBranchStatus("*",1);
  fou->cd();
  TeNew = Te->CloneTree();
  TeNew->Write(0,TObject::kOverwrite);
  gAlice->Write(0,TObject::kOverwrite);
  current->cd();
  delete alhe;
  cout<<"AliRun object written to file\n";
}


void AliITSSD2D(TString inFile, TString outFile){
  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  char ftmp[50];
  sprintf(ftmp,"%s",inFile.Data());
  manager->SetInputStream(0,ftmp);
  if(outFile != "")manager->SetOutputFile(outFile);
  AliITSDigitizer *dITS  = new AliITSDigitizer(manager);
  manager->Exec("");
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
  TFile *file2= new TFile(outFile,"update");
  writeAR(file,file2);
  delete dITS;
  delete manager;
  file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
  if(file){
    file->Close();
    delete file;
  }
  file2->Close();
  delete file2;
}

void AliITSH2FR2files (TString inFile, TString outfile)
{

  TFile * file = AccessInpFile(inFile);
  if(!file){
    cerr<<"*** AliITSH2FR2files: Problems in accessing the input file\n";
    return;
  }
 
  // open output file and create TreeR on it
  TFile * file2 = gAlice->InitTreeFile("R",outfile);
  file2->cd();
  // get ITS  
  AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
  if (!ITS) {
	cerr<<"ITSH2FR2files.C : AliITS object not found on file" << endl;
	return;
  }  // end if !ITS

  // Set the simulation model
 
  for (Int_t i=0;i<3;i++) {
    ITS->SetSimulationModel(i,new AliITSsimulationFastPoints());
  }
 

  //
  // Event Loop
  //
  Int_t nsignal=25;
  Int_t size=-1;
  TStopwatch timer;
  for (Int_t nevent=0; nevent<gAlice->TreeE()->GetEntries(); nevent++) {
    gAlice->GetEvent(nevent);
    if(!gAlice->TreeR())gAlice->MakeTree("R",file2);
    ITS->MakeBranch("RF"); 
    ITS->SetTreeAddress();  
    Int_t bgr_ev=nevent/nsignal;
    ITS->HitsToFastRecPoints(nevent,bgr_ev,size," ","All"," ");
  }

  timer.Stop(); timer.Print();

  delete gAlice;
  gAlice = 0;
  file->Close();
}

void copyTK(TString inFile, TString outFile) {

  TDirectory *current = gDirectory;
  TParticle *part = new TParticle();
  TTree *Tk;
  TTree *TkNew;
  TFile *fin = AccessInpFile(inFile);
  TFile *fou = new TFile(outFile,"update");
  char tname[30];
  for(Int_t event= 0; event < gAlice->TreeE()->GetEntries(); event++){
    current->cd();
    sprintf(tname,"TreeK%d",event);
    Tk = (TTree*)fin->Get(tname);
    if(!Tk){
      cerr<<"Trying to access a non-existent TTree : "<<tname<<endl;
    }
    else {
      Tk->SetBranchAddress("Particles",&part);
      Tk->SetBranchStatus("*",1);
      fou->cd();
      TkNew = Tk->CloneTree();
      TkNew->Write();
    }
  }
  delete gAlice;
  gAlice = 0;
  fin->Close();
  fou->Close();
  delete fin;
  delete fou;
  current->cd();

}
