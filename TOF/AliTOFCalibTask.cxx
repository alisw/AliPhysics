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
$Log
*/

// task to perform TOF calibration
// optimized for pp runs
// expect a maximum of 100 entries per channel
// different ways to calibrate are implemented
// C. Zampolli 

#include <TChain.h>
#include <TObject.h> 
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH1D.h>
#include <TH2F.h>
  
#include "AliTOFCalibTask.h" 
#include "AliTOFChannelTask.h" 
#include "AliESDEvent.h" 
#include "AliESDtrack.h" 
#include "AliLog.h"
//______________________________________________________________________________
AliTOFCalibTask::AliTOFCalibTask(const char *name) : 
  AliAnalysisTask(name,""),  
  fdir(0),
  fChain(0),
  fESD(0),
  fMinEntries(0),
  fIch(0),
  fToT(0),
  fTime(0),
  fExpTimePi(0),
  fExpTimeKa(0),
  fExpTimePr(0),
  ftree(0x0),
  fMinTime(0),
  fbigarray(0x0),
  findexarray(0x0),
  fnESD(0),
  fnESDselected(0),
  fnESDkTOFout(0),
  fnESDkTIME(0),
  fnESDassTOFcl(0),
  fnESDTIMEcut(0),
  fnESDTRDcut(0),
  fhToT(0),
  fhTime(0),
  fhExpTimePi(0),
  fhExpTimeKa(0),
  fhExpTimePr(0),
  fhPID(0),
  fhch(0),
  fOutputContainer(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple

  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
  fdir = gDirectory->GetPath();
  fbigarray = new Float_t*[TOFCHANNELS];
  findexarray = new Int_t[TOFCHANNELS];

  for (Int_t i =0;i<TOFCHANNELS;i++){
    fbigarray[i] = new Float_t[CHENTRIES];
    findexarray[i]=0;
    for (Int_t j =0;j<CHENTRIES;j++){
      fbigarray[i][j]=-1;
    }
  }
}

//______________________________________________________________________________
AliTOFCalibTask::AliTOFCalibTask(const AliTOFCalibTask &calibtask) : 
  AliAnalysisTask("AliTOFCalibTask",""),  
  fdir(0),
  fChain(0),
  fESD(0),
  fMinEntries(0),
  fIch(0),
  fToT(0),
  fTime(0),
  fExpTimePi(0),
  fExpTimeKa(0),
  fExpTimePr(0),
  ftree(0x0),
  fMinTime(0),
  fbigarray(0x0),
  findexarray(0x0),
  fnESD(0),
  fnESDselected(0),
  fnESDkTOFout(0),
  fnESDkTIME(0),
  fnESDassTOFcl(0),
  fnESDTIMEcut(0),
  fnESDTRDcut(0),
  fhToT(0),
  fhTime(0),
  fhExpTimePi(0),
  fhExpTimeKa(0),
  fhExpTimePr(0),
  fhPID(0),
  fhch(0),
  fOutputContainer(0)
{
  // Copy Constructor.
  fdir=calibtask.fdir;
  fChain=calibtask.fChain;
  fESD=calibtask.fESD;
  fMinEntries=calibtask.fMinEntries;
  fIch=calibtask.fIch;
  fToT=calibtask.fToT;
  fTime=calibtask.fTime;
  fExpTimePi=calibtask.fExpTimePi;
  fExpTimeKa=calibtask.fExpTimeKa;
  fExpTimePr=calibtask.fExpTimePr;
  ftree=calibtask.ftree;
  fMinTime=calibtask.fMinTime;
  fnESD=calibtask.fnESD;
  fnESDselected=calibtask.fnESDselected;
  fnESDkTOFout=calibtask.fnESDkTOFout;
  fnESDkTIME=calibtask.fnESDkTIME;
  fnESDassTOFcl=calibtask.fnESDassTOFcl;
  fnESDTIMEcut=calibtask.fnESDTIMEcut;
  fnESDTRDcut=calibtask.fnESDTRDcut;
  fhToT=calibtask.fhToT;
  fhTime=calibtask.fhTime;
  fhExpTimePi=calibtask.fhExpTimePi;
  fhExpTimeKa=calibtask.fhExpTimeKa;
  fhExpTimePr=calibtask.fhExpTimePr;
  fhPID=calibtask.fhPID;
  fhch=calibtask.fhch;
  fOutputContainer=calibtask.fOutputContainer; 

  fbigarray = new Float_t*[TOFCHANNELS];
  findexarray = new Int_t[TOFCHANNELS];

  for (Int_t i =0;i<TOFCHANNELS;i++){
    fbigarray[i] = new Float_t[CHENTRIES];
    findexarray[i]=calibtask.findexarray[i];
    for (Int_t j =0;j<CHENTRIES;j++){
      fbigarray[i][j]=calibtask.fbigarray[i][j];
    }
  }

}
//______________________________________________________________________________
AliTOFCalibTask:: ~AliTOFCalibTask() 
{
  // destructor

  AliInfo("TOF Calib Task: Deleting");
  delete[] fbigarray;
  delete findexarray;
  delete fOutputContainer;
  delete fhToT;
  delete fhTime;
  delete fhExpTimePi; 
  delete fhExpTimeKa; 
  delete fhExpTimePr; 
  delete fhPID; 
  delete fhch;
}
//______________________________________________________________________________
AliTOFCalibTask& AliTOFCalibTask::operator=(const AliTOFCalibTask &calibtask)  
{ 
   //assignment operator
  this->fdir=calibtask.fdir;
  this->fChain=calibtask.fChain;
  this->fESD=calibtask.fESD;
  this->fMinEntries=calibtask.fMinEntries;
  this->fIch=calibtask.fIch;
  this->fToT=calibtask.fToT;
  this->fTime=calibtask.fTime;
  this->fExpTimePi=calibtask.fExpTimePi;
  this->fExpTimeKa=calibtask.fExpTimeKa;
  this->fExpTimePr=calibtask.fExpTimePr;
  this->ftree=calibtask.ftree;
  this->fMinTime=calibtask.fMinTime;
  this->fnESD=calibtask.fnESD;
  this->fnESDselected=calibtask.fnESDselected;
  this->fnESDkTOFout=calibtask.fnESDkTOFout;
  this->fnESDkTIME=calibtask.fnESDkTIME;
  this->fnESDassTOFcl=calibtask.fnESDassTOFcl;
  this->fnESDTIMEcut=calibtask.fnESDTIMEcut;
  this->fnESDTRDcut=calibtask.fnESDTRDcut;
  this->fhToT=calibtask.fhToT;
  this->fhTime=calibtask.fhTime;
  this->fhExpTimePi=calibtask.fhExpTimePi;
  this->fhExpTimeKa=calibtask.fhExpTimeKa;
  this->fhExpTimePr=calibtask.fhExpTimePr;
  this->fOutputContainer=calibtask.fOutputContainer; 
  this->fhPID=calibtask.fhPID;
  this->fhch=calibtask.fhch;

  this->fbigarray = new Float_t*[TOFCHANNELS];
  this->findexarray = new Int_t[TOFCHANNELS];
  for (Int_t i =0;i<TOFCHANNELS;i++){
    this->fbigarray[i] = new Float_t[CHENTRIES];
    this->findexarray[i]=calibtask.findexarray[i];
    for (Int_t j =0;j<CHENTRIES;j++){
      this->fbigarray[i][j]=calibtask.fbigarray[i][j];
    }
  }
  return *this;
}
//--------------------------------------------------------------------------
void AliTOFCalibTask::BookHistos(){

  //booking histos

  AliInfo(Form("*** Booking Histograms %s", GetName())) ; 

  fhToT=
    new TH1F("hToT", " ToT distribution (ns) ", 400, 0, 40);
  fhTime=
    new TH1F("hTime", " Time distribution (ns) ", 400, 0, 40);
  fhExpTimePi=
    new TH1F("hExpTimePi", " Exp Time distribution, Pi (ns) ", 400, 0, 40);
  fhExpTimeKa=
    new TH1F("hExpTimeKa", " Exp Time distribution, Ka (ns) ", 400, 0, 40);
  fhExpTimePr=
    new TH1F("hExpTimePr", " Exp Time distribution, Pr (ns) ", 400, 0, 40);
  fhPID=
    new TH1I("hPID", " Combinatorial PID identity ", 3, 0, 3);
  fhch=
    new TH1D("hch", " TOF channel ", TOFCHANNELS, 0, TOFCHANNELS);

  //create the putput container
  fOutputContainer = new TObjArray(7) ; 
  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fhToT,             0) ; 
  fOutputContainer->AddAt(fhTime,            1) ; 
  fOutputContainer->AddAt(fhExpTimePi,       2) ; 
  fOutputContainer->AddAt(fhExpTimeKa,       3) ; 
  fOutputContainer->AddAt(fhExpTimePr,       4) ; 
  fOutputContainer->AddAt(fhPID,             5) ; 
  fOutputContainer->AddAt(fhch,              6) ; 

}
//----------------------------------------------------------------------------
void AliTOFCalibTask::DrawHistos(){

  // drawing output histos

  AliInfo(Form("*** Drawing Histograms %s", GetName())) ; 

  TCanvas * canvasToTTime = new TCanvas("canvasToTTime", " ToT and Time ",400, 30, 550, 630) ;
  canvasToTTime->Divide(1,2);
  canvasToTTime->cd(1);
  fhToT->SetLineColor(4);
  fhToT->GetXaxis()->SetTitle("ToT (ns)");
  fhToT->Draw("hist");
  canvasToTTime->cd(2);
  fhTime->SetLineColor(4);
  fhTime->GetXaxis()->SetTitle("Time (ns)");
  fhTime->Draw("hist");
  canvasToTTime->Update();
  canvasToTTime->Print("ToTTime.gif");

  TCanvas * canvasExpTime = new TCanvas("canvasExpTime", " Expected Times ",400, 30, 550, 630) ;
  canvasExpTime->Divide(1,3);
  canvasExpTime->cd(1);
  fhExpTimePi->SetLineColor(4);
  fhExpTimePi->GetXaxis()->SetTitle("Exp Time (ns), #pi");
  fhExpTimePi->Draw("hist");
  canvasExpTime->cd(2);
  fhExpTimeKa->SetLineColor(4);
  fhExpTimeKa->GetXaxis()->SetTitle("Exp Time (ns), K");
  fhExpTimeKa->Draw("hist");
  canvasExpTime->cd(3);
  fhExpTimePr->SetLineColor(4);
  fhExpTimePr->GetXaxis()->SetTitle("Exp Time (ns), p");
  fhExpTimePr->Draw("hist");

  canvasExpTime->Print("ExpTime.gif");

  TCanvas * canvasPID = new TCanvas("canvasPID", " Combinatorial PID ",400, 30, 550, 400);
  fhPID->GetXaxis()->SetTitle("Comb PID");
  fhPID->GetXaxis()->SetBinLabel(1,"#pi");
  fhPID->GetXaxis()->SetBinLabel(2,"K");
  fhPID->GetXaxis()->SetBinLabel(3,"p");
  fhPID->Draw("hist");

  canvasPID->Print("PID.gif");

  TCanvas * canvasrndch = new TCanvas("canvasrndch", " TOF channel ",400, 30, 550, 400);
  fhch->GetXaxis()->SetTitle("TOF ch");
  fhch->Draw("hist");
  Float_t meanTOFch = 0;
  for (Int_t ibin=0;ibin<TOFCHANNELS;ibin++){
    meanTOFch+=(Float_t)fhch->GetBinContent(ibin+1);
  }

  meanTOFch/=TOFCHANNELS;
  AliDebug(1,Form(" Mean number of tracks/channel = %f ",meanTOFch));

  canvasrndch->Print("rndch.gif");

  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.gif", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.gif"); 
  gROOT->ProcessLine(line);
  AliInfo(Form("*** TOF Calib Task: plots saved in %s.tar.gz...\n", GetName())) ;
}

//______________________________________________________________________________
void AliTOFCalibTask::ConnectInputData(const Option_t*)
{
  // Initialization of branch container and histograms 
    
  //  AliLog::SetClassDebugLevel("AliTOFCalibTask",1);
  AliInfo(Form("*** Initialization of %s", GetName())) ; 
  
  // Get input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return ;
  }
  
  // One should first check if the branch address was taken by some other task
  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESDEvent*)(*address);
  } else {
    fESD = new AliESDEvent();
    AliInfo(" qui ok ");
    //    fESD = (AliESDEvent*)fChain->GetTree()->GetUserInfo()->FindObject("AliESDEvent");
    fESD->ReadFromTree(fChain) ;  
  }

  BookHistos();

}
//-----------------------------------------------------------------------
Bool_t AliTOFCalibTask::Notify()
{
  // Initialisation of branch container and histograms 
    
  AliInfo(Form("*** We are in  %s::Notify()", GetName())) ; 
  
  // Get input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return kFALSE;
  }
  
  // One should first check if the branch address was taken by some other task
  char ** address = (char **)GetBranchAddress(0, "ESD") ;
  if (address) 
    fESD = (AliESDEvent *)(*address) ; 
  else {
    fESD = new AliESDEvent();
    //    fESD = (AliESDEvent*)fChain->GetTree()->GetUserInfo()->FindObject("AliESDEvent");
    fESD->ReadFromTree(fChain) ;  
    //SetBranchAddress(0,"ESD",&fESD);
  }

  return kTRUE;
}

//________________________________________________________________________
void AliTOFCalibTask::CreateOutputObjects()
{
// Create histograms
}

//______________________________________________________________________________
void AliTOFCalibTask::Exec(Option_t * opt) 
{

  // main 

  AliInfo(Form("*** Executing %s", GetName())) ; 

//******* The loop over events -----------------------------------------------

  // Processing of one event
  Long64_t entry = fChain->GetReadEntry() ;  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }

  if ( !((entry-1)%100) ) 
    AliDebug(1,Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  // ************************  TOF *************************************

  fMinTime=22E3;   //ns; not used
  Int_t ntrk = fESD->GetNumberOfTracks() ;
  fnESD+=ntrk;
  Int_t nselected = 0;
  Int_t itr = -1;
  while ( ntrk-- ) {
    itr++;
    AliESDtrack * t = fESD->GetTrack(ntrk) ;
    //selecting only good quality tracks
    if (!Select(t)) continue;
    nselected++;
    Int_t ich = Int_t(t->GetTOFCalChannel()); 
    fhch->Fill(ich);
    AliDebug(2,Form(" ESD in channel %i, filling array %i",t->GetTOFCalChannel(),ich));
    Float_t tot = t->GetTOFsignalToT();
    Float_t time = t->GetTOFsignalRaw();
    AliDebug(2,Form(" track # %i in schannel %i, time = %f \n",ntrk,ich,time)); 
    Double_t expTime[10]; 
    t->GetIntegratedTimes(expTime);
    Float_t expTimePi = expTime[2]*1.E-3;
    Float_t expTimeKa = expTime[3]*1.E-3;
    Float_t expTimePr = expTime[4]*1.E-3;
    if (findexarray[ich]==(Int_t)(CHENTRIES/NIDX)) {
      AliInfo(Form("too many tracks in channel %i, not storing any more...",ich));
      continue;
    } 
    findexarray[ich]++;
    AliDebug(2,Form("tracks in channel %i = %i, storing... ",ich, findexarray[ich] ));
    Int_t ientry=(findexarray[ich]-1)*NIDX;
    fbigarray[ich][ientry+DELTAIDXTOT]=tot;  //in ns
    fbigarray[ich][ientry+DELTAIDXTIME]=time*1E-3; // in ns
    fbigarray[ich][ientry+DELTAIDXEXTIMEPI]=expTimePi;
    fbigarray[ich][ientry+DELTAIDXEXTIMEKA]=expTimeKa;
    fbigarray[ich][ientry+DELTAIDXEXTIMEPR]=expTimePr;
    fhToT->Fill(fbigarray[ich][ientry+DELTAIDXTOT]);
    fhTime->Fill(fbigarray[ich][ientry+DELTAIDXTIME]);
    fhExpTimePi->Fill(fbigarray[ich][ientry+DELTAIDXEXTIMEPI]);
    fhExpTimeKa->Fill(fbigarray[ich][ientry+DELTAIDXEXTIMEKA]);
    fhExpTimePr->Fill(fbigarray[ich][ientry+DELTAIDXEXTIMEPR]);
    AliDebug(2, Form("track = %i, tot = %f, time = %f, and Exp time in TOF: pi = %f, K = %f, p = %f",itr, fbigarray[ich][ientry+DELTAIDXTOT], fbigarray[ich][ientry+DELTAIDXTIME], expTimePi,expTimeKa,expTimePr));
  }
  fnESDselected+=nselected;

  PostData(0, fOutputContainer);  
}

//_____________________________________________________________________________
void AliTOFCalibTask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  
  // some plots

  TH1::AddDirectory(0);
  AliInfo("TOF Calib Task: End of events loop");
  AliInfo(Form(" Number of analyzed ESD tracks: %i\n",fnESD));
  AliInfo(Form(" Number of selected ESD tracks: %i\n",fnESDselected));
  AliInfo(Form(" Number of ESD tracks with kTOFout: %i\n",fnESDkTOFout));
  AliInfo(Form(" Number of ESD tracks with kTIME: %i\n",fnESDkTIME));
  AliInfo(Form(" Number of ESD tracks with TRDcut: %i\n",fnESDTRDcut));
  AliInfo(Form(" Number of ESD tracks with TIMEcut: %i\n",fnESDTIMEcut));
  AliInfo(Form(" Number of ESD tracks with assTOFcl: %i\n",fnESDassTOFcl));

  for (Int_t i = 0;i<TOFCHANNELS;i++){
    Int_t size=findexarray[i]*NIDX;
    AliDebug(2, Form(" entries %i in channel %i ",findexarray[i],i));
    if (findexarray[i]<6) {
      AliDebug(1, Form(" not enough statistics for combined PID for channel %i, putting all the tracks as if they were pions",i));
      continue;
    }
    if (!CombPID(&fbigarray[i][0], size)) AliError("ERROR!!!!ERROR!!!");
  }
  
  DrawHistos();

  // saving data in a tree
  
  TDirectory *dir = gDirectory;
  TFile *outFile = 0;
  TTree * tree = 0x0;
  Bool_t isThere=kFALSE;
  const char *dirname = "./";
  TString filename= "TOFCalib.root";

  if((gSystem->FindFile(dirname,filename))!=NULL){
    isThere=kTRUE;
  }
  
  if(!isThere){
    AliInfo("File with tree for Calibration not there, creating it");
    outFile = new TFile( "TOFCalib.root","RECREATE");
    outFile->cd();
    tree = new TTree("T", "Tree for TOF calibration");
    Float_t p[CHENTRIESSMALL];
    Int_t nentries;
    tree->Branch("nentries",&nentries,"nentries/I");
    tree->Branch("TOFentries",p,"TOFentries[nentries]/F");
    for (Int_t i=0;i<TOFCHANNELS;i++){
      nentries=findexarray[i]*(NIDXSMALL); // when filling small array, 
                                           // only first 3 floats taken 
                                           // into account
      for (Int_t j=0; j<findexarray[i];j++){
	for (Int_t k=0; k<NIDXSMALL;k++){
	  Int_t index1= j*NIDXSMALL+k;   // index in small array
	  Int_t index2=j*NIDX+k;         // index in big array
	  p[index1]=fbigarray[i][index2];
	}
      }
      tree->Fill();
    }
    tree->Write("",TObject::kOverwrite);
    outFile->Close();
    //delete outFile;
    //outFile=0x0;
    //delete tree;
    //tree=0x0;
  }
  else{
    AliInfo("File with tree for Calibration already there, updating it");
    outFile = new TFile( "TOFCalib.root","READ");
    outFile->cd();
    tree=(TTree*)outFile->Get("T");
    Float_t p[CHENTRIESSMALL];
    Float_t ptemp[MAXCHENTRIESSMALL];
    Int_t nentries;
    Int_t nentries1;
    Int_t nentriestemp;
    tree->SetBranchAddress("nentries",&nentries);
    tree->SetBranchAddress("TOFentries",p);
    TFile * outFile1 = new TFile( "TOFCalib.root","RECREATE");
    outFile1->cd();
    TTree *treeTemp = new TTree("T","Tree for TOF calibration");
    treeTemp->Branch("nentries",&nentriestemp,"nentries/I");
    treeTemp->Branch("TOFentries",ptemp,"TOFentries[nentries]/F");
    for (Int_t i=0;i<TOFCHANNELS;i++){
      tree->GetEntry(i);
      nentries1=nentries;
      if (nentries+findexarray[i]*(NIDXSMALL) > MAXCHENTRIESSMALL){
	AliDebug(2, Form("findexarray[%i] = %i, nentries = %i",i,findexarray[i],nentries)); 
	AliInfo(Form("Too many entries in channel %i, stopping at %i  entries, no storing any more...",i,MAXCHENTRIESSMALL/NIDXSMALL));
	findexarray[i]=(MAXCHENTRIESSMALL-nentries)/NIDXSMALL;
      }
      for (Int_t kk = 0;kk<nentries;kk++){
	ptemp[kk]=p[kk];
      }
      nentries+=findexarray[i]*(NIDXSMALL); // when filling small array, 
	                                    // only first 3 floats taken 
                                            // into account 
      nentriestemp=nentries; 
      for (Int_t j=0; j<findexarray[i];j++){
	for (Int_t k=0; k<NIDXSMALL;k++){
	  Int_t index1= j*NIDXSMALL+k;   // index in small array
	  Int_t index2=j*NIDX+k;   // index in big array
	  ptemp[index1+nentries1]=fbigarray[i][index2];
	}
      }
      treeTemp->Fill();
    }
    treeTemp->SetName("T");
    treeTemp->SetTitle("Tree for TOF calibration");
    outFile->Close();
    treeTemp->Write("",TObject::kOverwrite);
    outFile1->Close();
  }
  dir->cd();
  Int_t calibrationStatus = 0;
  //  Int_t ch[2]={3,1003};
  //  calibrationStatus = Calibrate(2,ch,"save");
  //calibrationStatus = Calibrate(3,4,"save");
  //calibrationStatus = CalibrateFromProfile(3,"save");
  //  calibrationStatus = Calibrate(3,"save");
  calibrationStatus = Calibrate("");
  AliInfo(Form("Calibration Status = %i, bye.... \n",calibrationStatus));
  
}
//----------------------------------------------------------------------------
Int_t AliTOFCalibTask::Calibrate(Int_t ichmin, Int_t ichmax, Option_t *optionSave, Option_t *optionFit){

  // calibrating summing more than one channels
  // computing calibration parameters
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays
  
  TH1::AddDirectory(0);

  AliInfo(Form("*** Calibrating Histograms %s, summing more channels, from channel %i, to channel %i, storing Calib Pars in channel %i", GetName(),ichmin,ichmax,ichmin)) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 
  const char *dirname = "./";
  TString filename= "TOFCalib.root";
  
  if((gSystem->FindFile(dirname,filename))==NULL){
    AliInfo("No file with tree for calibration found! exiting....");
    return 1;
  }

  TFile *file = new TFile("TOFCalib.root","READ");
  TTree *tree = (TTree*)file->Get("T");
  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  tree->SetBranchAddress("nentries",&nentries);
  tree->SetBranchAddress("TOFentries",p);

  Float_t nentriesmean =0;
  for (Int_t i=ichmin; i<ichmax; i++){
    tree->GetEntry(i);
    Int_t ntracks=nentries/3;
    nentriesmean+=ntracks;
  }

  //  nentriesmean/=(ichmax-ichmin);
  if (nentriesmean < MEANENTRIES) {
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",nentriesmean));
    return 2;
  }

  //filling ToT and Time arrays

  Int_t nbinToT = 100;  // ToT bin width in Profile = 48.8 ps 
  Float_t minToT = 0;   // ns
  Float_t maxToT = 4.88;  // ns
  TObjArray * arrayCal = new TObjArray(TOFCHANNELS);

  TH1F *hToT = new TH1F("htot","htot",nbinToT, minToT, maxToT);
  TH1F *hdeltaTime = new TH1F("hdeltaTime","hdeltaTime",200,2,4);
  Int_t ntracksTot = 0;
  Int_t ntracks = 0;
  Double_t binsProfile[101]; // sized larger than necessary, the correct 
                             // dim being set in the booking of the profile
  Int_t nusefulbins=0;
  Float_t meantime=0;
  for (Int_t i = ichmin;i<ichmax;i++){
    tree->GetEntry(i);
    ntracksTot+=nentries/3;
    ntracks=nentries/3;
    AliDebug(2,Form("channel %i, nentries = %i, ntracks = %i",i,nentries, ntracks));
    for (Int_t j=0;j<ntracks;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      hdeltaTime->Fill(p[idxexTime]-p[idxexExTime]);
      meantime+=p[idxexTime]-p[idxexExTime];
      hToT->Fill(tot);
    }
  }
  nusefulbins = FindBins(hToT,&binsProfile[0]);
  meantime/=ntracksTot;
  AliDebug(2, Form("meantime = %f",meantime));
  
  for (Int_t j=1;j<=nusefulbins;j++) {
    AliDebug(2,Form(" summing channels from %i to %i, nusefulbins = %i, bin %i = %f",ichmin,ichmax,nusefulbins,j,binsProfile[j])); 
  }

  TProfile* hSlewingProf = new TProfile("hSlewingProf", "hSlewingProf",nusefulbins, binsProfile, "G");  // CHECK THE BUILD OPTION, PLEASE!!!!!!
  TH2F * htimetot = new TH2F("htimetot","htimetot",nbinToT, minToT, maxToT,600,-5,10);
  for (Int_t i=ichmin; i<ichmax; i++){
    tree->GetEntry(i);
    ntracks=nentries/3;
    for (Int_t j=0;j<ntracks;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      Float_t time = p[idxexTime]-p[idxexExTime];
      AliDebug (2, Form("track = %i, time = %f, tot = %f, time-meantime = %f",j,time, tot, time-meantime));
      hSlewingProf->Fill(tot,time);  // if meantime is not used, the fill may be moved in the loop above
      htimetot->Fill(tot,time-meantime);  // if meantime is not used, the fill may be moved in the loop above
    }
  }
  hSlewingProf->Fit("pol5",optionFit,"",0,4);
  TF1 * calibfunc = (TF1*)hSlewingProf->GetFunction("pol5");
  Float_t par[6];    
  for(Int_t kk=0;kk<6;kk++){
    par[kk]=calibfunc->GetParameter(kk);
    AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
  }

  if(strstr(optionSave,"save")){
    TFile * fileProf = new TFile("TOFCalibSave.root","recreate");
    fileProf->cd(); 
    TString profName=Form("Profile%06i_%06i",ichmin,ichmax);
    TString timeTotName=Form("TimeTot%06i_%06i",ichmin,ichmax);
    TString totName=Form("Tot%06i_%06i",ichmin,ichmax);
    TString deltaName=Form("Delta%06i_%06i",ichmin,ichmax);
    hSlewingProf->Write(profName);
    htimetot->Write(timeTotName);
    hToT->Write(totName);
    hdeltaTime->Write(deltaName);
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }

  delete hToT;
  hToT=0x0;
  delete hSlewingProf;
  hSlewingProf=0x0;
  delete htimetot;
  htimetot=0x0;
  delete hdeltaTime;
  hdeltaTime=0x0;
  file->Close();
  delete file;
  file=0x0;

  AliTOFChannelTask * calChannel = new AliTOFChannelTask();
  calChannel->SetSlewPar(par);
  // saving parameters in chmin
  arrayCal->AddAt(calChannel,ichmin);
  TFile *filecalib = new TFile("outCalArray.root","RECREATE");
  arrayCal->Write("array",TObject::kSingleKey);
  filecalib->Close();
  delete filecalib;
  filecalib = 0x0;
  delete arrayCal;
  arrayCal = 0x0;
  delete calChannel;
  calChannel =0x0;
  return 0;
}
//----------------------------------------------------------------------------
Int_t AliTOFCalibTask::Calibrate(Int_t i, Option_t *optionSave, Option_t *optionFit){

  // computing calibration parameters for channel i
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays

  TH1::AddDirectory(0);
  
  AliInfo(Form("*** Calibrating Histograms (one channel) %s", GetName())) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 
  const char *dirname = "./";
  TString filename= "TOFCalib.root";
  
  if((gSystem->FindFile(dirname,filename))==NULL){
    AliInfo("No file with tree for calibration found! exiting....");
    return 1;
  }

  TFile *file = new TFile("TOFCalib.root","READ");
  TTree *tree = (TTree*)file->Get("T");
  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  tree->SetBranchAddress("nentries",&nentries);
  tree->SetBranchAddress("TOFentries",p);

  tree->GetEntry(i);
  Int_t nentriesmean=nentries/3;

  if (nentriesmean < MEANENTRIES) {  
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",nentriesmean));
    return 2;
  }

  //filling ToT and Time arrays

  Int_t nbinToT = 100;  // ToT bin width in Profile = 48.8 ps 
  Float_t minToT = 0;   // ns
  Float_t maxToT = 4.88;  // ns
  TObjArray * arrayCal = new TObjArray(TOFCHANNELS);

  TH1F *hToT = new TH1F("htot","htot",nbinToT, minToT, maxToT);
  TH1F *hdeltaTime = new TH1F("hdeltaTime","hdeltaTime",200,2,4);
  Int_t ntracks = 0;
  Double_t binsProfile[101]; // sized larger than necessary, the correct 
                             // dim being set in the booking of the profile
  Int_t nusefulbins=0;
  Float_t meantime=0;
  ntracks=nentries/3;
  AliDebug(2,Form("channel %i, nentries = %i, ntracks = %i",i ,nentries, ntracks));
  for (Int_t j=0;j<ntracks;j++){
    Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
    Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
    Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
    Float_t tot = p[idxexToT];
    meantime+=p[idxexTime]-p[idxexExTime];
    hdeltaTime->Fill(p[idxexTime]-p[idxexExTime]);
    hToT->Fill(tot);
  }
  
  nusefulbins = FindBins(hToT,&binsProfile[0]);
  meantime/=ntracks;
  AliDebug(2,Form("meantime = %f",meantime));
  
  for (Int_t j=1;j<=nusefulbins;j++) {
    AliDebug(2,Form(" channel %i, nusefulbins = %i, bin %i = %f",i,nusefulbins,j,binsProfile[j])); 
  }

  TProfile* hSlewingProf = new TProfile("hSlewingProf", "hSlewingProf",nusefulbins, binsProfile, "G");  // CHECK THE BUILD OPTION, PLEASE!!!!!!
  TH2F * htimetot = new TH2F("htimetot","htimetot",nbinToT, minToT, maxToT,600,-5,10);
  tree->GetEntry(i);
  ntracks=nentries/3;
  for (Int_t j=0;j<ntracks;j++){
    Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
    Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
    Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
    Float_t tot = p[idxexToT];
    Float_t time = p[idxexTime]-p[idxexExTime];
    AliDebug (2,Form("track = %i, time = %f, tot = %f, time-meantime = %f",j,time, tot, time-meantime));
    hSlewingProf->Fill(tot,time);  // if meantime is not used, the fill may be moved in the loop above
    htimetot->Fill(tot,time-meantime);  // if meantime is not used, the fill may be moved in the loop above
  }
  hSlewingProf->Fit("pol5",optionFit,"",0,4);
  TF1 * calibfunc = (TF1*)hSlewingProf->GetFunction("pol5");
  Float_t par[6];    
  for(Int_t kk=0;kk<6;kk++){
    par[kk]=calibfunc->GetParameter(kk);
    AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
  }

  file->Close();
  delete file;
  file=0x0;

  if(strstr(optionSave,"save")){
    TFile * fileProf = new TFile("TOFCalibSave.root","recreate");
    fileProf->cd();   
    TString profName=Form("Profile%06i",i);
    TString timeTotName=Form("TimeTot%06i",i);
    TString totName=Form("Tot%06i",i);
    TString deltaName=Form("Delta%06i",i);
    hSlewingProf->Write(profName);
    htimetot->Write(timeTotName);
    hToT->Write(totName);
    hdeltaTime->Write(deltaName);
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }

  delete hToT;
  hToT=0x0; 
  delete hSlewingProf;
  hSlewingProf=0x0;
  delete htimetot;
  htimetot=0x0;
  delete hdeltaTime;
  hdeltaTime=0x0;

  AliTOFChannelTask * calChannel = new AliTOFChannelTask();
  calChannel->SetSlewPar(par);
  arrayCal->AddAt(calChannel,i);
  TFile *filecalib = new TFile("outCalArray.root","RECREATE");
  arrayCal->Write("array",TObject::kSingleKey);
  filecalib->Close();
  delete filecalib;
  filecalib=0x0;
  delete arrayCal;
  arrayCal = 0x0;
  delete calChannel;
  calChannel =0x0;
  
  return 0;
}
//----------------------------------------------------------------------------
Int_t AliTOFCalibTask::Calibrate(Int_t nch, Int_t *ch, Option_t *optionSave, Option_t *optionFit){

  // calibrating an array of channels
  // computing calibration parameters
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays
  
  TH1::AddDirectory(0);

  AliInfo(Form("*** Calibrating Histograms %s, number of channels = %i", GetName(),nch)) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 
  for (Int_t ich=0; ich<nch; ich++){
    Int_t i = ch[ich];
    AliInfo(Form("Calibrating channel = %i",i )) ; 
  }
  const char *dirname = "./";
  TString filename= "TOFCalib.root";
  
  if((gSystem->FindFile(dirname,filename))==NULL){
    AliInfo("No file with tree for calibration found! exiting....");
    return 1;
  }

  TFile *file = new TFile("TOFCalib.root","READ");
  TTree *tree = (TTree*)file->Get("T");
  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  tree->SetBranchAddress("nentries",&nentries);
  tree->SetBranchAddress("TOFentries",p);

  Float_t nentriesmean =0;
  for (Int_t ich=0; ich<nch; ich++){
    Int_t i = ch[ich];
    tree->GetEntry(i);
    Int_t ntracks=nentries/3;
    nentriesmean+=ntracks;
  }

  nentriesmean/=nch;
  if (nentriesmean < MEANENTRIES) { 
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",nentriesmean));
    return 2;
  }

  //filling ToT and Time arrays

  Int_t nbinToT = 100;  // ToT bin width in Profile = 48.8 ps 
  Float_t minToT = 0;   // ns
  Float_t maxToT = 4.88;  // ns
  TObjArray * arrayCal = new TObjArray(TOFCHANNELS);
  arrayCal->SetOwner();
  TFile * fileProf=0x0;
  if(strstr(optionSave,"save")){
    fileProf = new TFile("TOFCalibSave.root","recreate");
  }
  for (Int_t ich=0; ich<nch; ich++) {
    TH1F *hToT = new TH1F("htot","htot",nbinToT, minToT, maxToT);
    TH1F *hdeltaTime = new TH1F("hdeltaTime","hdeltaTime",200,2,4);
    Double_t binsProfile[101]; // sized larger than necessary, the correct 
                              // dim being set in the booking of the profile
    TH2F * htimetot = new TH2F("htimetot","htimetot",nbinToT, minToT, maxToT,600,-5,10);
    Int_t ntracksTot = 0;
    Int_t ntracks = 0;
    Int_t nusefulbins=0;
    Float_t meantime=0;
    Int_t i = ch[ich];
    AliDebug(2,Form("Calibrating channel %i",i));
    tree->GetEntry(i);
    ntracksTot+=nentries/3;
    ntracks=nentries/3;
    if (ntracks < MEANENTRIES) {
      AliInfo(Form(" Too small mean number of entires in channel %i (number of tracks = %f), not calibrating channel and continuing.....",i,ntracks));
      continue;
    }
    for (Int_t j=0;j<ntracks;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      hdeltaTime->Fill(p[idxexTime]-p[idxexExTime]);
      meantime+=p[idxexTime]-p[idxexExTime];
      hToT->Fill(tot);
    }

    nusefulbins = FindBins(hToT,&binsProfile[0]);
    meantime/=ntracksTot;
    for (Int_t j=1;j<=nusefulbins;j++) {
      AliDebug(2,Form(" channel %i, nusefulbins = %i, bin %i = %f",i,nusefulbins,j,binsProfile[j])); 
    }

    TProfile* hSlewingProf = new TProfile("hSlewingProf", "hSlewingProf",nusefulbins, binsProfile, "G");  // CHECK THE BUILD OPTION, PLEASE!!!!!!
    tree->GetEntry(i);
    ntracks=nentries/3;
    for (Int_t j=0;j<ntracks;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      Float_t time = p[idxexTime]-p[idxexExTime];
      AliDebug(2,Form("track = %i, time = %f, tot = %f, time-meantime = %f",j,time, tot, time-meantime));
      hSlewingProf->Fill(tot,time);  // if meantime is not used, the fill may be moved in the loop above
      htimetot->Fill(tot,time-meantime);  // if meantime is not used, the fill may be moved in the loop above
    }
    
    hSlewingProf->Fit("pol5",optionFit,"",1,4);
    TF1 * calibfunc = (TF1*)hSlewingProf->GetFunction("pol5");
    Float_t par[6];    
    for(Int_t kk=0;kk<6;kk++){
      par[kk]=calibfunc->GetParameter(kk);
      AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
    }
    
    if(strstr(optionSave,"save") && fileProf){
      TString profName=Form("Profile%06i",i);
      TString timeTotName=Form("TimeTot%06i",i);
      TString totName=Form("Tot%06i",i);
      TString deltaName=Form("Delta%06i",i);
      fileProf->cd();
      hSlewingProf->Write(profName);
      htimetot->Write(timeTotName);
      hToT->Write(totName);
      hdeltaTime->Write(deltaName);
    }

    AliTOFChannelTask * calChannel = new AliTOFChannelTask();
    calChannel->SetSlewPar(par);
    arrayCal->AddAt(calChannel,i);
    delete hToT;
    hToT=0x0;
    delete hSlewingProf;
    hSlewingProf=0x0;
    delete htimetot;
    htimetot=0x0;
    delete hdeltaTime;
    hdeltaTime=0x0;
  }

  if(strstr(optionSave,"save") && fileProf){
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }
  file->Close();
  delete file;
  file=0x0;
  TFile *filecalib = new TFile("outCalArray.root","RECREATE");
  arrayCal->Write("array",TObject::kSingleKey);
  filecalib->Close();
  delete filecalib;
  filecalib=0x0;
  arrayCal->Clear();
  delete arrayCal;
  arrayCal=0x0;

  return 0;
}
//----------------------------------------------------------------------------
Int_t AliTOFCalibTask::CalibrateFromProfile(Int_t i, Option_t *optionSave, Option_t *optionFit){

  // computing calibration parameters using the old profiling algo
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays

  TH1::AddDirectory(0);

  TObjArray * arrayCal = new TObjArray(TOFCHANNELS);
  AliInfo(Form("*** Calibrating Histograms From Profile %s", GetName())) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 
  const char *dirname = "./";
  TString filename= "TOFCalib.root";

  if((gSystem->FindFile(dirname,filename))==NULL){
    AliInfo("No file with tree for calibration found! exiting....");
    return 1;
  }

  TFile *file = new TFile("TOFCalib.root","READ");
  TTree *tree = (TTree*)file->Get("T");
  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  tree->SetBranchAddress("nentries",&nentries);
  tree->SetBranchAddress("TOFentries",p);

  tree->GetEntry(i);
  Int_t ntracks=nentries/3;

  if (ntracks < MEANENTRIES) {  
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",ntracks));
    return 2;
  }

  TH1F * hProf = new TH1F();
  hProf = Profile(i);
  hProf->Fit("pol5",optionFit,"",0,4);
  TF1 * calibfunc = (TF1*)hProf->GetFunction("pol5");
  Float_t par[6];    
  for(Int_t kk=0;kk<6;kk++){
    par[kk]=calibfunc->GetParameter(kk);
    AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
  }

  if(strstr(optionSave,"save")){
    TFile * fileProf = new TFile("TOFCalibSave.root","recreate");
    fileProf->cd(); 
    TString profName=Form("Profile%06i",i);
    hProf->Write(profName);
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }

  delete hProf;
  hProf=0x0;
  file->Close();
  delete file;
  file=0x0;
  AliTOFChannelTask * calChannel = new AliTOFChannelTask();
  calChannel->SetSlewPar(par);
  arrayCal->AddAt(calChannel,i);
  TFile *filecalib = new TFile("outCalArray.root","RECREATE");
  arrayCal->Write("array",TObject::kSingleKey);
  filecalib->Close();
  delete filecalib;
  filecalib = 0x0;
  delete calChannel;
  delete arrayCal;
  return 0;
}
//----------------------------------------------------------------------------
Int_t AliTOFCalibTask::Calibrate(Option_t *optionSave, Option_t *optionFit){

  // calibrating the whole TOF
  // computing calibration parameters
  // Returning codes:
  // 0 -> everything was ok
  // 1 -> no tree for calibration found
  // 2 -> not enough statistics to perform calibration
  // 3 -> problems with arrays

  TH1::AddDirectory(0);

  AliInfo(Form("*** Calibrating Histograms %s, all channels", GetName())) ; 
  AliInfo(Form("Option for Saving histos = %s",optionSave )) ; 
  AliInfo(Form("Option for Fitting Profile histos = %s",optionFit )) ; 
  const char *dirname = "./";
  TString filename= "TOFCalib.root";

  if((gSystem->FindFile(dirname,filename))==NULL){
    AliInfo("No file with tree for calibration found! exiting....");
    return 1;
  }

  TFile * fileProf=0x0;
  if(strstr(optionSave,"save")){
    fileProf = new TFile("TOFCalibSave.root","recreate");
  }
  TFile *file = new TFile("TOFCalib.root","READ");
  TTree *tree = (TTree*)file->Get("T");
  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  tree->SetBranchAddress("nentries",&nentries);
  tree->SetBranchAddress("TOFentries",p);

  Float_t nentriesmean =0;
  for (Int_t i=0; i<TOFCHANNELS; i++){
    tree->GetEntry(i);
    Int_t ntracks=nentries/3;
    nentriesmean+=ntracks;
  }

  nentriesmean/=TOFCHANNELS;
  if (nentriesmean < MEANENTRIES) {
    AliInfo(Form(" Too small mean number of entires per channel (mean number = %f) not calibrating and exiting.....",nentriesmean));
    return 2;
  }

  //filling ToT and Time arrays

  Int_t nbinToT = 100;  // ToT bin width in Profile = 50.0 ps 
  Float_t minToT = 0;   // ns
  Float_t maxToT = 4.88;// ns
  TObjArray * arrayCal = new TObjArray(TOFCHANNELS);
  arrayCal->SetOwner();
  for (Int_t ii=0; ii<TOFCHANNELS; ii++) {
    TH1F *hToT = new TH1F("htot","htot",nbinToT, minToT, maxToT);
    TH1F *hdeltaTime = new TH1F("hdeltaTime","hdeltaTime",200,2,4);
    TH2F * htimetot = new TH2F("htimetot","htimetot",nbinToT, minToT, maxToT,600,-5,10);
    if (ii%1000 == 0) AliDebug(1,Form("Calibrating channel %i ",ii));
    Int_t i = 3;
    Int_t nusefulbins=0;
    Double_t binsProfile[101]; // sized larger than necessary, the correct 
                              // dim being set in the booking of the profile
    tree->GetEntry(i);
    Int_t ntracks=nentries/3;
    if (ntracks < MEANENTRIES) {
      AliInfo(Form(" Too small mean number of entires in channel %i (number of tracks = %f), not calibrating channel and continuing.....",i,ntracks));
      continue;
    }
    Float_t meantime=0;
    for (Int_t j=0;j<ntracks;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      hdeltaTime->Fill(p[idxexTime]-p[idxexExTime]);
      meantime+=p[idxexTime]-p[idxexExTime];
      hToT->Fill(tot);
    }
    nusefulbins = FindBins(hToT,&binsProfile[0]);
    meantime/=ntracks;
    for (Int_t j=0;j<nusefulbins;j++) {
      AliDebug(2,Form(" channel %i, usefulbin = %i, low edge = %f",i,j,binsProfile[j])); 
    }
    TProfile* hSlewingProf = new TProfile("hSlewingProf", "hSlewingProf",nusefulbins, binsProfile, "G");  // CHECK THE BUILD OPTION, PLEASE!!!!!!
    for (Int_t j=0;j<ntracks;j++){
      Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
      Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
      Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
      Float_t tot = p[idxexToT];
      Float_t time = p[idxexTime]-p[idxexExTime];
      AliDebug (2,Form("track = %i, time = %f, tot = %f, time-meantime = %f",j,time, tot, time-meantime));
      hSlewingProf->Fill(tot,time);  // if meantime is not used, the fill may be moved in the loop above
      htimetot->Fill(tot,time-meantime);  // if meantime is not used, the fill may be moved in the loop above
    }

    //    hSlewingProf->Fit("pol5",optionFit,"",1,4);
    //TF1 * calibfunc = (TF1*)hSlewingProf->GetFunction("pol5");
    Float_t par[6];    
    for(Int_t kk=0;kk<6;kk++){
      //      par[kk]=calibfunc->GetParameter(kk);
      par[kk]=kk;
      AliDebug(2,Form("parameter %i = %f",kk,par[kk]));
    }

    if(strstr(optionSave,"save") && fileProf){
      TString profName=Form("Profile%06i",ii);
      TString timeTotName=Form("TimeTot%06i",ii);
      TString totName=Form("Tot%06i",ii);
      TString deltaName=Form("Delta%06i",ii);
      fileProf->cd();
      hSlewingProf->Write(profName);
      htimetot->Write(timeTotName);
      hToT->Write(totName);
      hdeltaTime->Write(deltaName);
    }
    AliTOFChannelTask * calChannel = new AliTOFChannelTask();
    calChannel->SetSlewPar(par);
    arrayCal->AddAt(calChannel,ii);

    delete hToT;
    hToT=0x0;
    delete hSlewingProf;
    hSlewingProf=0x0;
    delete htimetot;
    htimetot=0x0;
    delete hdeltaTime;
    hdeltaTime=0x0;
  }

  if(strstr(optionSave,"save")){
    fileProf->Close();
    delete fileProf;
    fileProf=0x0;
  }

  file->Close();
  delete file;
  file=0x0;
  TFile *filecalib = new TFile("outCalArray.root","RECREATE");
  arrayCal->Write("array", TObject::kSingleKey);
  filecalib->Close();
  delete filecalib;
  filecalib=0x0;
  arrayCal->Clear();
  delete arrayCal;
  arrayCal=0x0;

  return 0;
}

//_____________________________________________________________________________

Bool_t AliTOFCalibTask::Select(AliESDtrack *t){

  // selecting good quality tracks
  //track selection: reconstrution to TOF:
  if ((t->GetStatus()&AliESDtrack::kTOFout)==0) {
    return 0;
  }
  fnESDkTOFout++;
  //IsStartedTimeIntegral
  if ((t->GetStatus()&AliESDtrack::kTIME)==0) {
    return 0;
  }
  fnESDkTIME++;
  if (t->GetStatus() & AliESDtrack::kTRDbackup) { 
    Float_t xout= t->GetOuterParam()->GetX();
    if (xout<364.25 &&  xout > 300.) return 0;
  }
  fnESDTRDcut++;
  Double_t time=t->GetTOFsignal();	
  time*=1.E-3; // tof given in nanoseconds
  if(time >= fMinTime){
    return 0;
  }
  fnESDTIMEcut++;
  
  Double_t mom=t->GetP();
  if (!(mom<=UPPERMOMBOUND && mom>=LOWERMOMBOUND)){
  //  return 0;  // skipping momentum cut
  } 
   
  UInt_t assignedTOFcluster=t->GetTOFcluster();//index of the assigned TOF cluster, >0 ?
  if(assignedTOFcluster==0){ // not matched
    return 0;
  }
  fnESDassTOFcl++;
  return 1;
}
//_____________________________________________________________________________

Bool_t AliTOFCalibTask::CombPID(Float_t *smallarray, Int_t size){

  //track Combinatorial PID for calibration

  //cout << " In comb PID " << endl;
  Float_t t0offset=0;
  const Int_t kntracksinset=6;
  if (kntracksinset != 6) {
    AliDebug(1,"Number of tracks in set smaller than expected one, identifying every particle as if it was a pion!");
    return 0;
  }
  Float_t exptof[kntracksinset][3];
  Int_t   assparticle[kntracksinset];
  Float_t timeofflight[kntracksinset];
  Float_t texp[kntracksinset];
  Float_t timezero[kntracksinset];
  Float_t weightedtimezero[kntracksinset];
  Float_t besttimezero[kntracksinset];
  Float_t bestchisquare[kntracksinset];
  Float_t bestweightedtimezero[kntracksinset];

  for (Int_t i=0;i<kntracksinset;i++){
    assparticle[i]=3;
    timeofflight[i]=0;
    texp[i]=0;
    timezero[i]=0;
    weightedtimezero[i]=0;
    besttimezero[i]=0;
    bestchisquare[i]=0;
    bestweightedtimezero[i]=0;
  }

  Int_t nset= (Int_t)(size/kntracksinset/NIDX);
  for (Int_t i=0; i< nset; i++) {   
    for (Int_t j=0; j<kntracksinset; j++) {
      Int_t idxtime = ((kntracksinset*i+j)*NIDX)+DELTAIDXTIME; 
      Int_t idxextimePi = ((kntracksinset*i+j)*NIDX)+DELTAIDXEXTIMEPI; 
      Int_t idxextimeKa = ((kntracksinset*i+j)*NIDX)+DELTAIDXEXTIMEKA; 
      Int_t idxextimePr = ((kntracksinset*i+j)*NIDX)+DELTAIDXEXTIMEPR; 
      Double_t time=smallarray[idxtime]; // TOF time in ns
      timeofflight[j]=time+t0offset;
      exptof[j][0]=smallarray[idxextimePi];
      exptof[j][1]=smallarray[idxextimeKa];
      exptof[j][2]=smallarray[idxextimePr];
      AliDebug(2,Form("j = %i, Time = %f, and Exp time in PID: pi = %f, K = %f, p = %f",j,timeofflight[j],exptof[j][0],exptof[j][1],exptof[j][2]));
    }
    Float_t t0best=999.;
    Float_t et0best=999.;
    Float_t chisquarebest=999.;
    for (Int_t i1=0; i1<3;i1++) {
      texp[0]=exptof[0][i1];
      for (Int_t i2=0; i2<3;i2++) { 
      	texp[1]=exptof[1][i2];
	for (Int_t i3=0; i3<3;i3++) {
	  texp[2]=exptof[2][i3];
	  for (Int_t i4=0; i4<3;i4++) {
	    texp[3]=exptof[3][i4];
	    for (Int_t i5=0; i5<3;i5++) {
	      texp[4]=exptof[4][i5];
	      for (Int_t i6=0; i6<3;i6++) {
		texp[5]=exptof[5][i6];
					
		Float_t sumAllweights=0.;
		Float_t meantzero=0.;
		Float_t emeantzero=0.;
		
		for (Int_t itz=0; itz<kntracksinset;itz++) {
		  timezero[itz]=texp[itz]-timeofflight[itz];		    
		  weightedtimezero[itz]=timezero[itz]/TRACKERROR;
		  sumAllweights+=1./TRACKERROR;
		  meantzero+=weightedtimezero[itz];
		  
		} // end loop for (Int_t itz=0; itz<15;itz++)
		
		meantzero=meantzero/sumAllweights; // it is given in [ns]
		emeantzero=sqrt(1./sumAllweights); // it is given in [ns]
		
		// calculate chisquare
		
		Float_t chisquare=0.;		
		for (Int_t icsq=0; icsq<kntracksinset;icsq++) {
		  chisquare+=(timezero[icsq]-meantzero)*(timezero[icsq]-meantzero)/TRACKERROR;
		} // end loop for (Int_t icsq=0; icsq<15;icsq++) 
		
		Int_t npion=0;
		if(i1==0)npion++;
		if(i2==0)npion++;
		if(i3==0)npion++;
		if(i4==0)npion++;
		if(i5==0)npion++;
		if(i6==0)npion++;
		
	     	if(chisquare<=chisquarebest  && ((Float_t) npion/ ((Float_t) kntracksinset)>0.3)){
		  for(Int_t iqsq = 0; iqsq<kntracksinset; iqsq++) {
		    besttimezero[iqsq]=timezero[iqsq]; 
		    bestweightedtimezero[iqsq]=weightedtimezero[iqsq]; 
		    bestchisquare[iqsq]=(timezero[iqsq]-meantzero)*(timezero[iqsq]-meantzero)/TRACKERROR; 
		  }
		  
		  assparticle[0]=i1;
		  assparticle[1]=i2;
		  assparticle[2]=i3;
		  assparticle[3]=i4;
		  assparticle[4]=i5;
		  assparticle[5]=i6;
		  
		  chisquarebest=chisquare;
     		  t0best=meantzero;
		  et0best=emeantzero;
		} // close if(dummychisquare<=chisquare)
	      } // end loop on i6
	    } // end loop on i5
	  } // end loop on i4
	} // end loop on i3
      } // end loop on i2
    } // end loop on i1
		
    Float_t confLevel=999;
    if(chisquarebest<999.){
      Double_t dblechisquare=(Double_t)chisquarebest;
      confLevel=(Float_t)TMath::Prob(dblechisquare,kntracksinset-1); 
    }
    // assume they are all pions for fake sets
    if(confLevel<0.01 || confLevel==999. ){
      for (Int_t itrk=0; itrk<kntracksinset; itrk++)assparticle[itrk]=0;
    }

    AliDebug(2,Form(" Best Assignment, selection %i %i %i %i %i %i ",assparticle[0],assparticle[1],assparticle[2],assparticle[3],assparticle[4],assparticle[5]));

    for (Int_t kk=0;kk<kntracksinset;kk++){
      Int_t idxextimePi = ((kntracksinset*i+kk)*NIDX)+DELTAIDXEXTIMEPI; 
      Int_t idxextimeKa = ((kntracksinset*i+kk)*NIDX)+DELTAIDXEXTIMEKA; 
      Int_t idxextimePr = ((kntracksinset*i+kk)*NIDX)+DELTAIDXEXTIMEPR; 
      // storing in third slot of smallarray the assigned expected time
      fhPID->Fill(assparticle[kk]);
      //      assparticle[kk]=0;  //assuming all particles are pions
      if (assparticle[kk]==0){       //assigned to be a Pi
	smallarray[idxextimeKa]=-1;
	smallarray[idxextimePr]=-1;
      }
      else if (assparticle[kk]==1){  //assigned to be a Ka
	smallarray[idxextimePi]=smallarray[idxextimeKa];
	smallarray[idxextimeKa]=-1;
	smallarray[idxextimePr]=-1;
      }
      else if (assparticle[kk]==2){  //assigned to be a Pr
	smallarray[idxextimePi]=smallarray[idxextimePr];
	smallarray[idxextimeKa]=-1;
	smallarray[idxextimePr]=-1;
      }
    }
  }
  return 1;
}
//-----------------------------------------------------------------------
TH1F* AliTOFCalibTask::Profile(Int_t ich)
{
  // profiling algo

  TFile *file = new TFile("TOFCalib.root","READ");
  TTree *tree = (TTree*)file->Get("T");
  Float_t p[MAXCHENTRIESSMALL];
  Int_t nentries;
  tree->SetBranchAddress("nentries",&nentries);
  tree->SetBranchAddress("TOFentries",p);

  //Prepare histograms for Slewing Correction
  const Int_t knbinToT = 100;
  Int_t nbinTime = 200;
  Float_t minTime = -5.5; //ns
  Float_t maxTime = 5.5; //ns
  Float_t minToT = 0; //ns
  Float_t maxToT = 5.; //ns
  Float_t deltaToT = (maxToT-minToT)/knbinToT;
  Double_t mTime[knbinToT+1],mToT[knbinToT+1],meanTime[knbinToT+1], meanTime2[knbinToT+1],vToT[knbinToT+1], vToT2[knbinToT+1],meanToT[knbinToT+1],meanToT2[knbinToT+1],vTime[knbinToT+1],vTime2[knbinToT+1],xlow[knbinToT+1],sigmaTime[knbinToT+1];
  Int_t n[knbinToT+1], nentrx[knbinToT+1];
  Double_t sigmaToT[knbinToT+1];
  for (Int_t i = 0; i < knbinToT+1 ; i++){
    mTime[i]=0;
    mToT[i]=0;
    n[i]=0;
    meanTime[i]=0;
    meanTime2[i]=0;
    vToT[i]=0;
    vToT2[i]=0;
    meanToT[i]=0;
    meanToT2[i]=0;
    vTime[i]=0;
    vTime2[i]=0;
    xlow[i]=0;
    sigmaTime[i]=0;
    sigmaToT[i]=0;
    n[i]=0;
    nentrx[i]=0;
  }
  TH2F* hSlewing = new TH2F("hSlewing", "hSlewing", knbinToT, minToT, maxToT, nbinTime, minTime, maxTime);
  Int_t ntracks = 0;
  TH1F *histo = new TH1F("histo", "1D Time vs ToT", knbinToT, minToT, maxToT);
  tree->GetEntry(ich);
  ntracks=nentries/3;
  for (Int_t j=0;j<ntracks;j++){
    Int_t idxexToT = (j* NIDXSMALL)+DELTAIDXTOT; 
    Int_t idxexTime = (j* NIDXSMALL)+DELTAIDXTIME; 
    Int_t idxexExTime = (j* NIDXSMALL)+DELTAIDXPID; 
    Float_t tot = p[idxexToT];
    Float_t time = p[idxexTime]-p[idxexExTime];
    Int_t nx = (Int_t)((tot-minToT)/deltaToT)+1;
    if ((tot != 0) && ( time!= 0)){
      vTime[nx]+=time;
      vTime2[nx]+=time*time;
      vToT[nx]+=tot;
      vToT2[nx]+=tot*tot;
      nentrx[nx]++;
      hSlewing->Fill(tot,time);
    }
  }
  Int_t nbinsToT=hSlewing->GetNbinsX();
  if (nbinsToT != knbinToT) {
    AliError("Profile :: incompatible numbers of bins");
    return 0x0;
  }
  
  Int_t usefulBins=0;
  for (Int_t i=1;i<=nbinsToT;i++){
    if (nentrx[i]!=0){
      n[usefulBins]+=nentrx[i];
      if (n[usefulBins]==0 && i == nbinsToT) {
	break;
      }
      meanTime[usefulBins]+=vTime[i];
      meanTime2[usefulBins]+=vTime2[i];
      meanToT[usefulBins]+=vToT[i];
      meanToT2[usefulBins]+=vToT2[i];
      if (n[usefulBins]<10 && i!=nbinsToT) continue; 
      mTime[usefulBins]=meanTime[usefulBins]/n[usefulBins];
      mToT[usefulBins]=meanToT[usefulBins]/n[usefulBins];
      sigmaTime[usefulBins]=TMath::Sqrt(1./n[usefulBins]/n[usefulBins]
			    *(meanTime2[usefulBins]-meanTime[usefulBins]
			    *meanTime[usefulBins]/n[usefulBins]));
      if ((1./n[usefulBins]/n[usefulBins]
	   *(meanToT2[usefulBins]-meanToT[usefulBins]
	     *meanToT[usefulBins]/n[usefulBins]))< 0) {
	AliError(" too small radical" );
	sigmaToT[usefulBins]=0;
      }
      else{       
	sigmaToT[usefulBins]=TMath::Sqrt(1./n[usefulBins]/n[usefulBins]
			     *(meanToT2[usefulBins]-meanToT[usefulBins]
			     *meanToT[usefulBins]/n[usefulBins]));
      }
      usefulBins++;
    }
  }
  for (Int_t i=0;i<usefulBins;i++){
    Int_t binN = (Int_t)((mToT[i]-minToT)/deltaToT)+1;
    histo->Fill(mToT[i],mTime[i]);
    histo->SetBinError(binN,sigmaTime[i]);
  } 
  delete file;
  file=0x0;
  delete hSlewing;
  hSlewing=0x0;

  return histo;
}
//----------------------------------------------------------------------------
Int_t AliTOFCalibTask::FindBins(TH1F* h, Double_t *binsProfile) const{

  // to determine the bins for ToT histo

  Int_t cont = 0;
  Int_t startBin = 1;
  Int_t nbin = h->GetNbinsX();
  Int_t nentries = (Int_t)h->GetEntries();
  Float_t max = h->GetBinLowEdge(nbin);
  Int_t nusefulbins=0;
  Int_t maxcont=0;
  // setting maxvalue of entries per bin
  if (nentries <= 60) maxcont = 2;
  else  if (nentries <= 100) maxcont = 5;
  else  if (nentries <= 500) maxcont = 10;
  else  maxcont = 20;
  for (Int_t j=1;j<=nbin;j++) {
    cont += (Int_t)h->GetBinContent(j);
    if (j<nbin){
      if (cont>=maxcont){
	nusefulbins++;
	binsProfile[nusefulbins-1]=h->GetBinLowEdge(startBin);
	cont=0;
	startBin=j+1;
	continue;
      }
    }
    else{
      if (cont>=maxcont){
	nusefulbins++;
	binsProfile[nusefulbins-1]=h->GetBinLowEdge(startBin);
	binsProfile[nusefulbins]=max;
      }
      else {
	binsProfile[nusefulbins]=h->GetBinLowEdge(startBin);
      }
    }
  }
  return nusefulbins;
}
