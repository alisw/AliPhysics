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
$Log$
Revision 1.5  2007/10/05 10:38:10  zampolli
Oversight fixed

Revision 1.4  2007/10/04 15:36:44  zampolli
Updates to new TOF offline calibration schema

Revision 1.3  2007/07/31 07:26:16  zampolli
Bug fixed in the $ field

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
#include <TDirectory.h>
#include <TGrid.h>
  
#include "AliTOFCalibTask.h" 
#include "AliESDEvent.h" 
#include "AliESDtrack.h" 
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBStorage.h"

//______________________________________________________________________________
AliTOFCalibTask::AliTOFCalibTask(const char *name) : 
  AliAnalysisTask(name,""),  
  fdir(0),
  fChain(0),
  fESD(0),
  fToT(0),
  fTime(0),
  fExpTimePi(0),
  fExpTimeKa(0),
  fExpTimePr(0),
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
  fOutputContainer(0),
  frun(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple

  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
  fdir = gDirectory->GetPath();
  fbigarray = new Float_t*[TOFCHANNELS];
  findexarray = new Int_t[TOFCHANNELS];

  for (Int_t i=0;i<11;i++){
    fassparticle[i]=-1;
  } 

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
  fToT(0),
  fTime(0),
  fExpTimePi(0),
  fExpTimeKa(0),
  fExpTimePr(0),
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
  fOutputContainer(0),
  frun(0)
{
  // Copy Constructor.
  fdir=calibtask.fdir;
  fChain=calibtask.fChain;
  fESD=calibtask.fESD;
  fToT=calibtask.fToT;
  fTime=calibtask.fTime;
  fExpTimePi=calibtask.fExpTimePi;
  fExpTimeKa=calibtask.fExpTimeKa;
  fExpTimePr=calibtask.fExpTimePr;
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
  frun=calibtask.frun;
  fOutputContainer=calibtask.fOutputContainer; 

  fbigarray = new Float_t*[TOFCHANNELS];
  findexarray = new Int_t[TOFCHANNELS];

  for (Int_t i=0;i<11;i++){
    fassparticle[i]=calibtask.fassparticle[i];
  } 

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
  this->fToT=calibtask.fToT;
  this->fTime=calibtask.fTime;
  this->fExpTimePi=calibtask.fExpTimePi;
  this->fExpTimeKa=calibtask.fExpTimeKa;
  this->fExpTimePr=calibtask.fExpTimePr;
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

  for (Int_t i=0;i<11;i++){
    this->fassparticle[i]=calibtask.fassparticle[i];
  } 
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
  }
  fESD->ReadFromTree(fChain) ;  

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
  
  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESDEvent*)(*address);
  } else {
    fESD = new AliESDEvent();
  }
  fESD->ReadFromTree(fChain) ;  

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
    //    ich=3; //only for debug purpose
    AliDebug(2,Form(" ESD in channel %i, filling array %i",t->GetTOFCalChannel(),ich));
    Float_t tot = t->GetTOFsignalToT();
    Float_t time = t->GetTOFsignalRaw();
    AliDebug(2,Form(" track # %i in channel %i, time = %f \n",ntrk,ich,time)); 
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
  TDirectory *dir = gDirectory;
  AliInfo("TOF Calib Task: End of events loop");
  for (Int_t ich = 0;ich<TOFCHANNELS;ich++){
    if (findexarray[ich]>0){
      Int_t ncuttime = SelectOnTime(&fbigarray[ich][0],findexarray[ich],ich);
      fnESDselected-=ncuttime;
    }
  }

  AliInfo(Form(" Number of analyzed ESD tracks: %i\n",fnESD));
  AliInfo(Form(" Number of selected ESD tracks: %i\n",fnESDselected));
  AliInfo(Form(" Number of ESD tracks with kTOFout: %i\n",fnESDkTOFout));
  AliInfo(Form(" Number of ESD tracks with kTIME: %i\n",fnESDkTIME));
  AliInfo(Form(" Number of ESD tracks with TRDcut: %i\n",fnESDTRDcut));
  AliInfo(Form(" Number of ESD tracks with TIMEcut: %i\n",fnESDTIMEcut));
  AliInfo(Form(" Number of ESD tracks with assTOFcl: %i\n",fnESDassTOFcl));

  for (Int_t i = 0;i<TOFCHANNELS;i++){
    Int_t size=findexarray[i]*NIDX;
    if (i==3) AliInfo(Form(" entries %i in channel %i ",findexarray[i],i));
    AliDebug(2, Form(" entries %i in channel %i ",findexarray[i],i));
    if (findexarray[i]<=2) {
      AliDebug(1, Form(" not enough statistics for combined PID for channel %i, putting all the tracks as if they were pions",i));
      continue;
    }
    if (!CombPID(&fbigarray[i][0], size)) AliError("ERROR!!!!ERROR!!!");
  }
  
  
  DrawHistos();

  // saving data in a tree
  
  AliInfo("Building tree for Calibration");
  TTree * tree = new TTree("T", "Tree for TOF calibration");
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
  
  AliInfo("Putting tree for calibration in Reference data");

  // grid file option
  Char_t filename[100];
  sprintf(filename,"alien:///alice/cern.ch/user/c/czampolli/TOFCalibReference_%i.root",frun);
  TGrid::Connect("alien://");
  TFile *filegrid = TFile::Open(filename,"CREATE");
  tree->Write();
  dir->cd();
  
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

Int_t AliTOFCalibTask::SelectOnTime(Float_t *charray, Int_t ntracks, Int_t ich){

  // discarding tracks with time-mintime < MINTIME
  
  Int_t ndeleted=0;
  Float_t mintime = 1E6;
  for (Int_t itr = 0;itr<ntracks;itr++){
    Int_t ientry=itr*NIDX;
    Float_t time = charray[ientry+DELTAIDXTIME];// in ns
      if (time<mintime) mintime = time;
  }
  AliDebug(1,Form("Mintime for channel %i = %f",ich, mintime));
  for (Int_t itr = 0;itr<ntracks;itr++){
    Int_t ientry=itr*NIDX;
    Float_t time = charray[ientry+DELTAIDXTIME];// in ns
    if ((time-mintime)>MINTIME) {
      ndeleted++;
      AliDebug(1,Form("Deleting %i track from channel %i, time = %f",ndeleted, ich, time));
      findexarray[ich]--;
      for (Int_t j=itr+1;j<ntracks;j++){
	Int_t ientry=j*NIDX;
	Int_t idxtot = ientry+DELTAIDXTOT; 
	Int_t idxtime = ientry+DELTAIDXTIME; 
	Int_t idxextimePi = ientry+DELTAIDXEXTIMEPI; 
	Int_t idxextimeKa = ientry+DELTAIDXEXTIMEKA; 
	Int_t idxextimePr = ientry+DELTAIDXEXTIMEPR;
	Int_t ientrydel=(j-1)*NIDX;
	Int_t idxtotdel = ientrydel+DELTAIDXTOT; 
	Int_t idxtimedel = ientrydel+DELTAIDXTIME; 
	Int_t idxextimePidel = ientrydel+DELTAIDXEXTIMEPI; 
	Int_t idxextimeKadel = ientrydel+DELTAIDXEXTIMEKA; 
	Int_t idxextimePrdel = ientrydel+DELTAIDXEXTIMEPR;
	charray[idxtotdel]=charray[idxtot]; 
	charray[idxtimedel]=charray[idxtime]; 
	charray[idxextimePidel]=charray[idxextimePi]; 
	charray[idxextimeKadel]=charray[idxextimeKa]; 
	charray[idxextimePrdel]=charray[idxextimePr];
      }
    } 
  }
  return ndeleted;
}
//_____________________________________________________________________________

Bool_t AliTOFCalibTask::CombPID(Float_t *smallarray, Int_t size){

  // track Combinatorial PID for calibration, only when more than 2 particles
  // fall in channel

  Float_t t0offset=0;

  if (size/NIDX<=2){
    AliDebug(1,"Number of tracks in channel smaller than 2, identifying every particle as if it was a pion!");
    return 0;
  }

  Int_t ntracksinchannel = size/NIDX;
  Int_t ntrkinset=-1;
  Int_t nset = ntracksinchannel/6;

  if (nset ==0) {
    nset=1;
    if (ntracksinchannel < 6) {
      AliDebug(2,"Number of tracks in set smaller than 6, Combinatorial PID still applied to this set.");
    }
  }

  AliInfo(Form("number of sets = %i", nset));
  // loop on sets
  for (Int_t i=0; i< nset; i++) {   
    if (i<nset-1)ntrkinset=6;
    else {
      if (ntracksinchannel<6){
	ntrkinset=ntracksinchannel;
      }
      else{
	ntrkinset=6+Int_t((ntracksinchannel)%6);
      }
    }
    AliInfo(Form("set = %i, number of tracks in set = %i", i,ntrkinset));
      
    for (Int_t ii=0;ii<ntrkinset;ii++){
      fassparticle[ii]=-1;
    }
    Float_t **exptof;
    Float_t* timeofflight;
    Float_t* texp;
    Int_t* index;

    exptof = new Float_t*[ntrkinset]; 
    timeofflight = new Float_t[ntrkinset];
    texp = new Float_t[ntrkinset];
    index = new Int_t[ntrkinset];
    for (Int_t ii=0;ii<ntrkinset;ii++){
      exptof[ii] = new Float_t[3];
      timeofflight[ii]=0;
      texp[ii]=0;
      index[ii]=-1;
    }

    for (Int_t j=0; j<ntrkinset; j++) {
      Int_t idxtime = ((6*i+j)*NIDX)+DELTAIDXTIME; 
      Int_t idxextimePi = ((6*i+j)*NIDX)+DELTAIDXEXTIMEPI; 
      Int_t idxextimeKa = ((6*i+j)*NIDX)+DELTAIDXEXTIMEKA; 
      Int_t idxextimePr = ((6*i+j)*NIDX)+DELTAIDXEXTIMEPR; 
      AliDebug(2,Form("idxtime = %i, idxextimePi = %i, idxextimeKa = %i, idxextimePr = %i", idxtime, idxextimePi, idxextimeKa, idxextimePr));
      Double_t time=smallarray[idxtime]; // TOF time in ns
      timeofflight[j]=time+t0offset;
      exptof[j][0]=smallarray[idxextimePi];
      exptof[j][1]=smallarray[idxextimeKa];
      exptof[j][2]=smallarray[idxextimePr];
      AliDebug(2,Form("j = %i, Time = %f, and Exp time in PID: pi = %f, K = %f, p = %f",j,timeofflight[j],exptof[j][0],exptof[j][1],exptof[j][2]));
    }
    
    
    Float_t chisquarebest=999.;
    AliInfo(Form(" Set = %i with %i tracks ", i,ntrkinset));
    chisquarebest = LoopCombPID(ntrkinset, ntrkinset,exptof,&texp[0],&timeofflight[0], &index[0],chisquarebest); 
    
    Float_t confLevel=999;
    if(chisquarebest<999.){
      Double_t dblechisquare=(Double_t)chisquarebest;
      confLevel=(Float_t)TMath::Prob(dblechisquare,ntrkinset-1); 
    }

    // assume they are all pions for fake sets
    if(confLevel<0.01 || confLevel==999. ){
      for (Int_t itrk=0; itrk<ntrkinset; itrk++)fassparticle[itrk]=0;
    }
    
    AliDebug(2,Form(" Best Assignment, selection %i %i %i %i %i %i %i %i %i %i",fassparticle[0],fassparticle[1],fassparticle[2],fassparticle[3],fassparticle[4],fassparticle[5],fassparticle[6],fassparticle[7],fassparticle[8],fassparticle[9]));

    for (Int_t kk=0;kk<ntrkinset;kk++){
      Int_t idxextimePi = ((6*i+kk)*NIDX)+DELTAIDXEXTIMEPI; 
      Int_t idxextimeKa = ((6*i+kk)*NIDX)+DELTAIDXEXTIMEKA; 
      Int_t idxextimePr = ((6*i+kk)*NIDX)+DELTAIDXEXTIMEPR; 
      // storing in third slot of smallarray the assigned expected time
      fhPID->Fill(fassparticle[kk]);
      //      fassparticle[kk]=0;  //assuming all particles are pions
      if (fassparticle[kk]==0){       //assigned to be a Pi
	smallarray[idxextimeKa]=-1;
	smallarray[idxextimePr]=-1;
      }
      else if (fassparticle[kk]==1){  //assigned to be a Ka
	smallarray[idxextimePi]=smallarray[idxextimeKa];
	smallarray[idxextimeKa]=-1;
	smallarray[idxextimePr]=-1;
      }
      else if (fassparticle[kk]==2){  //assigned to be a Pr
	smallarray[idxextimePi]=smallarray[idxextimePr];
	smallarray[idxextimeKa]=-1;
	smallarray[idxextimePr]=-1;
      }
    }
    
  }
  return 1;
}
//---------------------------------------------------------------------
Float_t  AliTOFCalibTask::LoopCombPID(Int_t itrkinset, Int_t ntrkinset, Float_t **exptof, Float_t *texp, Float_t *timeofflight, Int_t *index, Float_t chisquarebest){

  Int_t indextr = ntrkinset-itrkinset;
 
  for (index[indextr]=0;index[indextr]<3;index[indextr]++){
   Int_t ii = index[indextr];
   texp[indextr]=exptof[indextr][ii];
    if (indextr<ntrkinset-1){
      chisquarebest = LoopCombPID(ntrkinset-indextr-1,ntrkinset,exptof,&texp[0],&timeofflight[0],&index[0], chisquarebest);
    }
    
    else {
      Float_t *besttimezero = new Float_t[ntrkinset];
      Float_t *bestchisquare = new Float_t[ntrkinset];
      Float_t *bestweightedtimezero = new Float_t[ntrkinset];
      Float_t sumAllweights=0.;
      Float_t meantzero=0.;
      Float_t *weightedtimezero = new Float_t[ntrkinset];
      Float_t *timezero = new Float_t[ntrkinset];
      
      AliDebug(2,Form("Configuration = %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, so far chisquarebest = %f ",index[0],index[1],index[2],index[3],index[4],index[5],index[6],index[7],index[8],index[9],chisquarebest)); 
      for (Int_t itz=0; itz<ntrkinset;itz++) {
	timezero[itz]=texp[itz]-timeofflight[itz];		    
	weightedtimezero[itz]=timezero[itz]/TRACKERROR;
	sumAllweights+=1./TRACKERROR;
	meantzero+=weightedtimezero[itz];
      } // end loop for (Int_t itz=0; itz<ntrkinset;itz++)
      
      meantzero=meantzero/sumAllweights; // it is given in [ns]
      
      // calculate chisquare
      
      Float_t chisquare=0.;		
      for (Int_t icsq=0; icsq<ntrkinset;icsq++) {
	chisquare+=(timezero[icsq]-meantzero)*(timezero[icsq]-meantzero)/TRACKERROR;
      } // end loop for (Int_t icsq=0; icsq<ntrkinset;icsq++) 
      
      Int_t npion=0;
      for (Int_t j=0;j<ntrkinset;j++){
	if(index[j]==0)npion++;
      }
      
      if(chisquare<=chisquarebest  && ((Float_t) npion/ ((Float_t) ntrkinset)>0.3)){
	for(Int_t iqsq = 0; iqsq<ntrkinset; iqsq++) {
	  besttimezero[iqsq]=timezero[iqsq]; 
	  bestweightedtimezero[iqsq]=weightedtimezero[iqsq]; 
	  bestchisquare[iqsq]=(timezero[iqsq]-meantzero)*(timezero[iqsq]-meantzero)/TRACKERROR; 
	}
	
	for (Int_t j=0;j<ntrkinset;j++){
	  fassparticle[j]=index[j];
	}		  
	
	chisquarebest=chisquare;
      }
    }
  }
  return chisquarebest;
}
