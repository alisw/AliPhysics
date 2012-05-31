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
Revision 1.6  2007/10/05 10:44:11  zampolli
Oversight for debugging purpose fixed

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
#include <TGrid.h>
#include <TList.h>
#include <TArrayF.h>
#include <TBenchmark.h>
  
#include "AliTOFCalibTask.h" 
#include "AliESDEvent.h" 
#include "AliESDtrack.h" 
#include "AliLog.h"
#include "AliTOFArray.h"

//_________________________________________________________________________
AliTOFCalibTask::AliTOFCalibTask(const char *name) : 
  AliAnalysisTaskSE(name),  
  fESD(0),
  fToT(0),
  fTime(0),
  fExpTimePi(0),
  fExpTimeKa(0),
  fExpTimePr(0),
  fMinTime(0),
  fTOFArray(0x0),
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
  fhESD(0),
  fhESDselected(0),
  fhESDkTOFout(0),
  fhESDkTIME(0),
  fhESDassTOFcl(0),
  fhESDTIMEcut(0),
  fhESDTRDcut(0),
  fListOfHistos(0x0),
  fListArray(0x0)
{
  // Constructor.

  DefineOutput(1,TList::Class()); 
  DefineOutput(2,TList::Class()); 

  for (Int_t i=0;i<11;i++){
    fassparticle[i]=-1;
  } 

}

//______________________________________________________________________________
AliTOFCalibTask::AliTOFCalibTask(const AliTOFCalibTask &calibtask) : 
  AliAnalysisTaskSE("AliTOFCalibTask"),  
  fESD(calibtask.fESD),
  fToT(calibtask.fToT),
  fTime(calibtask.fTime),
  fExpTimePi(calibtask.fExpTimePi),
  fExpTimeKa(calibtask.fExpTimeKa),
  fExpTimePr(calibtask.fExpTimePr),
  fMinTime(calibtask.fMinTime),
  fTOFArray(calibtask.fTOFArray),
  fnESD(calibtask.fnESD),
  fnESDselected(calibtask.fnESDselected),
  fnESDkTOFout(calibtask.fnESDkTOFout),
  fnESDkTIME(calibtask.fnESDkTIME),
  fnESDassTOFcl(calibtask.fnESDassTOFcl),
  fnESDTIMEcutcalibtask.fnESDTIMEcut(),
  fnESDTRDcutcalibtask.fnESDTRDcut(),
  fhToT(calibtask.fhToT),
  fhTime(calibtask.fhTime),
  fhExpTimePi(calibtask.fhExpTimePi),
  fhExpTimeKa(calibtask.fhExpTimeKa),
  fhExpTimePr(calibtask.fhExpTimePr),
  fhPID(calibtask.fhPID),
  fhch(calibtask.fhch),
  fhESD(calibtask.fhESD),
  fhESDselected(calibtask.fhESDselected),
  fhESDkTOFout(calibtask.fhESDkTOFout),
  fhESDkTIME(calibtask.fhESDkTIME),
  fhESDassTOFcl(calibtask.fhESDassTOFcl),
  fhESDTIMEcut(calibtask.fhESDTIMEcut),
  fhESDTRDcut(calibtask.fhESDTRDcut),
  fListOfHistos(calibtask.fListOfHistos),
  fListArray(calibtask.fListArray)

{
  // Copy Constructor.

  for (Int_t i=0;i<11;i++){
    fassparticle[i]=calibtask.fassparticle[i];
  } 

}
//______________________________________________________________________________
AliTOFCalibTask:: ~AliTOFCalibTask() 
{
  // destructor

  AliInfo("TOF Calib Task: Deleting");
  
  delete fTOFArray;  
  delete fhToT;
  delete fhTime;
  delete fhExpTimePi; 
  delete fhExpTimeKa; 
  delete fhExpTimePr; 
  delete fhPID; 
  delete fhch;
  delete fhESD;
  delete fhESDselected;
  delete fhESDkTOFout;
  delete fhESDkTIME;
  delete fhESDTRDcut;
  delete fhESDTIMEcut;
  delete fhESDassTOFcl;
  //  delete fListOfHistos;
}
//______________________________________________________________________________
AliTOFCalibTask& AliTOFCalibTask::operator=(const AliTOFCalibTask &calibtask)  
{ 
   //assignment operator
  fESD=calibtask.fESD;
  fToT=calibtask.fToT;
  fTime=calibtask.fTime;
  fExpTimePi=calibtask.fExpTimePi;
  fExpTimeKa=calibtask.fExpTimeKa;
  fExpTimePr=calibtask.fExpTimePr;
  fMinTime=calibtask.fMinTime;
  fTOFArray=calibtask.fTOFArray;
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
  fhESD=calibtask.fhESD;
  fhESDselected=calibtask.fhESDselected;
  fhESDkTOFout=calibtask.fhESDkTOFout;
  fhESDkTIME=calibtask.fhESDkTIME;
  fhESDassTOFcl=calibtask.fhESDassTOFcl;
  fhESDTIMEcut=calibtask.fhESDTIMEcut;
  fhESDTRDcut=calibtask.fhESDTRDcut;
  fListOfHistos=calibtask.fListOfHistos;
  fListArray=calibtask.fListArray;

  for (Int_t i=0;i<11;i++){
    fassparticle[i]=calibtask.fassparticle[i];
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


  // create the output list of histos
  if (!fListOfHistos) fListOfHistos = new TList();
  fListOfHistos->SetOwner();

  fListOfHistos->AddAt(fhToT,             0) ; 
  fListOfHistos->AddAt(fhTime,            1) ; 
  fListOfHistos->AddAt(fhExpTimePi,       2) ; 
  fListOfHistos->AddAt(fhExpTimeKa,       3) ; 
  fListOfHistos->AddAt(fhExpTimePr,       4) ; 
  fListOfHistos->AddAt(fhPID,             5) ; 
  fListOfHistos->AddAt(fhch,              6) ; 

  fhESD=
    new TH1I("hESD","Number of analyzed ESDs",1,0,1);
  fhESDselected=
    new TH1I("hESDselected","Number of selected ESDs",1,0,1);
  fhESDkTOFout=
    new TH1I("hESDkTOFout","Number of ESDs with kTOFout",1,0,1);
  fhESDkTIME=
    new TH1I("hESDkTIME","Number of ESDs with kTime",1,0,1);
  fhESDTRDcut=
    new TH1I("hESDTRDcut","Number of ESDs with TRDcut",1,0,1);
  fhESDTIMEcut=
    new TH1I("hESDTIMEcut","Number of ESDs with TIMEcut",1,0,1);
  fhESDassTOFcl=
    new TH1I("hESDassTOFcl","Number of ESDs with assTOFcl",1,0,1);

  fListOfHistos->AddAt(fhESD,             7) ; 
  fListOfHistos->AddAt(fhESDselected,     8) ; 
  fListOfHistos->AddAt(fhESDkTOFout,      9) ; 
  fListOfHistos->AddAt(fhESDkTIME,        10) ; 
  fListOfHistos->AddAt(fhESDTRDcut,       11) ; 
  fListOfHistos->AddAt(fhESDTIMEcut,      12) ; 
  fListOfHistos->AddAt(fhESDassTOFcl,     13) ; 

  return;
  
}
//----------------------------------------------------------------------------
void AliTOFCalibTask::DrawHistos(){

  // drawing output histos

  AliInfo(Form("*** Drawing Histograms %s", GetName())) ; 

  if (!gROOT->IsBatch()){

	  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
	  
	  fhToT = (TH1F*)fListOfHistos->At(0);
	  fhTime = (TH1F*)fListOfHistos->At(1);
	  fhExpTimePi = (TH1F*)fListOfHistos->At(2);
	  fhExpTimeKa = (TH1F*)fListOfHistos->At(3);
	  fhExpTimePr = (TH1F*)fListOfHistos->At(4);
	  fhPID = (TH1I*)fListOfHistos->At(5);
	  fhch = (TH1D*)fListOfHistos->At(6);
	  
	  fhESD = (TH1I*)fListOfHistos->At(7);
	  fhESDselected = (TH1I*)fListOfHistos->At(8);
	  fhESDkTOFout = (TH1I*)fListOfHistos->At(9);
	  fhESDkTIME = (TH1I*)fListOfHistos->At(10);
	  fhESDassTOFcl = (TH1I*)fListOfHistos->At(11);
	  fhESDTIMEcut = (TH1I*)fListOfHistos->At(12);
	  fhESDTRDcut = (TH1I*)fListOfHistos->At(13);
	  
	  TCanvas * canvasToTTime = new TCanvas("canvasToTTime", " ToT and Time ",400, 30, 550, 630) ;
	  canvasToTTime->Divide(1,2);
	  canvasToTTime->cd(1);
	  fhToT->SetLineColor(4);
	  fhToT->GetXaxis()->SetTitle("ToT (ns)");
	  fhToT->DrawCopy("hist");
	  canvasToTTime->cd(2);
	  fhTime->SetLineColor(4);
	  fhTime->GetXaxis()->SetTitle("Time (ns)");
	  fhTime->DrawCopy("hist");
	  canvasToTTime->Update();
	  //canvasToTTime->Print("ToTTime.gif");
	  
	  TCanvas * canvasExpTime = new TCanvas("canvasExpTime", " Expected Times ",400, 30, 550, 630) ;
	  canvasExpTime->Divide(1,3);
	  canvasExpTime->cd(1);
	  fhExpTimePi->SetLineColor(4);
	  fhExpTimePi->GetXaxis()->SetTitle("Exp Time (ns), #pi");
	  fhExpTimePi->DrawCopy("hist");
	  canvasExpTime->cd(2);
	  fhExpTimeKa->SetLineColor(4);
	  fhExpTimeKa->GetXaxis()->SetTitle("Exp Time (ns), K");
	  fhExpTimeKa->DrawCopy("hist");
	  canvasExpTime->cd(3);
	  fhExpTimePr->SetLineColor(4);
	  fhExpTimePr->GetXaxis()->SetTitle("Exp Time (ns), p");
	  fhExpTimePr->DrawCopy("hist");
	  
	  //canvasExpTime->Print("ExpTime.gif");
	  
	  TCanvas * canvasPID = new TCanvas("canvasPID", " Combinatorial PID ",400, 30, 550, 400);
	  canvasPID->cd();
	  fhPID->GetXaxis()->SetTitle("Comb PID");
	  fhPID->GetXaxis()->SetBinLabel(1,"#pi");
	  fhPID->GetXaxis()->SetBinLabel(2,"K");
	  fhPID->GetXaxis()->SetBinLabel(3,"p");
	  fhPID->DrawCopy("hist");
	  
	  //canvasPID->Print("PID.gif");
	  
	  TCanvas * canvasrndch = new TCanvas("canvasrndch", " TOF channel ",400, 30, 550, 400);
	  canvasrndch->cd();
	  fhch->GetXaxis()->SetTitle("TOF ch");
	  fhch->Draw("hist");
	  Float_t meanTOFch = 0;
	  for (Int_t ibin=0;ibin<TOFCHANNELS;ibin++){
		  meanTOFch+=(Float_t)fhch->GetBinContent(ibin+1);
	  }
	  
	  meanTOFch/=TOFCHANNELS;
	  AliDebug(1,Form(" Mean number of tracks/channel = %f ",meanTOFch));
	  
	  //canvasrndch->Print("rndch.gif");

	  /*
	  char line[1024] ; 
	  sprintf(line, ".!tar -zcvf %s.tar.gz *.gif", GetName()) ; 
	  gROOT->ProcessLine(line);
	  sprintf(line, ".!rm -fR *.gif"); 
	  gROOT->ProcessLine(line);
	  AliInfo(Form("*** TOF Calib Task: plots saved in %s.tar.gz...\n", GetName())) ;
	  */
	  
	  AliInfo(Form(" Number of analyzed ESD tracks: %i",(Int_t)fhESD->GetEntries()));
	  AliInfo(Form(" Number of selected ESD tracks: %i",(Int_t)fhESDselected->GetEntries()));
	  AliInfo(Form(" Number of ESD tracks with kTOFout: %i",(Int_t)fhESDkTOFout->GetEntries()));
	  AliInfo(Form(" Number of ESD tracks with kTIME: %i",(Int_t)fhESDkTIME->GetEntries()));
	  AliInfo(Form(" Number of ESD tracks with TRDcut: %i",(Int_t)fhESDTRDcut->GetEntries()));
	  AliInfo(Form(" Number of ESD tracks with TIMEcut: %i",(Int_t)fhESDTIMEcut->GetEntries()));
	  AliInfo(Form(" Number of ESD tracks with assTOFcl: %i",(Int_t)fhESDassTOFcl->GetEntries()));
  }
  return;
}

//________________________________________________________________________
void AliTOFCalibTask::UserCreateOutputObjects()
{
// Create histograms and AliTOFArray

	AliInfo("UserCreateObjects");
	OpenFile(1);
	BookHistos();

	if (!fTOFArray) fTOFArray = new AliTOFArray(TOFCHANNELS);
	for (Int_t i =0;i<TOFCHANNELS;i++){
		fTOFArray->SetArray(i);
	}
	if (!fListArray) fListArray = new TList();
	fListArray->SetOwner();
	fListArray->Add(fTOFArray);
	return;
}

//______________________________________________________________________________
void AliTOFCalibTask::UserExec(Option_t * /*opt*/) 
{

  // main 


//******* The loop over events -----------------------------------------------

  AliVEvent*    fESD = fInputEvent ;
  if (!fESD) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }


  fMinTime=22E3;   //ns; not used
  Int_t ntrk = fESD->GetNumberOfTracks() ;
  fnESD+=ntrk;
  Int_t nselected = 0;
  Int_t itr = -1;
  while ( ntrk-- ) {
    fhESD->Fill(0);
    itr++;
    AliESDtrack * t = (AliESDtrack*)fESD->GetTrack(ntrk) ;
    //selecting only good quality tracks
    if (!Select(t)) continue;
    nselected++;
    fhESDselected->Fill(0);
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
    Int_t currentSize = fTOFArray->GetArraySize(ich);
    if (currentSize==CHENTRIES){ 
	    AliDebug(2,Form("too many tracks in channel %i, not storing any more...",ich));
      continue;
    } 
    AliDebug(2,Form("tracks in channel %i = %i, storing... ",ich, currentSize/NIDX ));

    fTOFArray->ReSetArraySize(ich,currentSize+5);
    fTOFArray->SetAt(ich,currentSize+DELTAIDXTOT,tot);
    fTOFArray->SetAt(ich,currentSize+DELTAIDXTIME,time*1E-3);
    fTOFArray->SetAt(ich,currentSize+DELTAIDXEXTIMEPI,expTimePi);
    fTOFArray->SetAt(ich,currentSize+DELTAIDXEXTIMEKA,expTimeKa);
    fTOFArray->SetAt(ich,currentSize+DELTAIDXEXTIMEPR,expTimePr);
    fhToT->Fill(fTOFArray->GetArrayAt(ich,currentSize+DELTAIDXTOT));
    fhTime->Fill(fTOFArray->GetArrayAt(ich,currentSize+DELTAIDXTIME));
    fhExpTimePi->Fill(fTOFArray->GetArrayAt(ich,currentSize+DELTAIDXEXTIMEPI));
    fhExpTimeKa->Fill(fTOFArray->GetArrayAt(ich,currentSize+DELTAIDXEXTIMEKA));
    fhExpTimePr->Fill(fTOFArray->GetArrayAt(ich,currentSize+DELTAIDXEXTIMEPR));
    AliDebug(2,Form("track = %i, tot = %f, time = %f, and Exp time in TOF: pi = %f, K = %f, p = %f",itr, fTOFArray->GetArrayAt(ich,currentSize+DELTAIDXTOT), fTOFArray->GetArrayAt(ich,currentSize+DELTAIDXTIME), expTimePi,expTimeKa,expTimePr));


  }
  fnESDselected+=nselected;
  PostData(1, fListOfHistos);  
  PostData(2, fListArray);  

}

//_____________________________________________________________________________
void AliTOFCalibTask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  
  // some plots

  TH1::AddDirectory(0);
  AliInfo("TOF Calib Task: End of events loop");
  DrawHistos();
  fListArray = dynamic_cast<TList*>(GetOutputData(2));
  fTOFArray = dynamic_cast<AliTOFArray*>(fListArray->At(0));

  for (Int_t ich = 0;ich<TOFCHANNELS;ich++){
	  Int_t entriesChannel = fTOFArray->GetArraySize(ich)/NIDX;
	  if (entriesChannel>0){
		  Int_t ncuttime = SelectOnTime(fTOFArray->GetArray(ich),entriesChannel,ich);
		  fnESDselected-=ncuttime;
	  }
  }

  TBenchmark bench;
  bench.Start("CombPID");
  for (Int_t i = 0;i<TOFCHANNELS;i++){
    Int_t size=fTOFArray->GetArraySize(i);
    AliDebug(2, Form(" entries %i in channel %i ",size/NIDX,i));
    if (size/NIDX<=2) {
      AliDebug(1, Form(" not enough statistics for combined PID for channel %i, putting all the tracks as if they were pions",i));
      continue;
    }
    if (i%1000 == 0) AliInfo(Form("At channel %d",i));
    if (!CombPID(fTOFArray->GetArray(i), size)) AliError("ERROR!!!!ERROR!!!");
  }
  bench.Stop("CombPID");
  bench.Print("CombPID");

  // saving data in a tree --> obsolete code; keeping for backup with new structure
  // using AliTOFArray
  /*  
  AliInfo("Building tree for Calibration");
  TTree * tree = new TTree("T", "Tree for TOF calibration");
  AliTOFArray* tempArray = new AliTOFArray(TOFCHANNELS);
  tree->Branch("TOFentries","AliTOFArray",&tempArray,32000,0);
  for (Int_t i = 0;i<TOFCHANNELS;i++){
	  Int_t sizeChannel = fTOFArray->GetArraySize(i)/NIDX;
      	  if (i==3) AliDebug(2,Form("Entries in channel %d = %d ",i,sizeChannel));  // just for debug 
	  if (sizeChannel!=0){
		  tempArray->SetArray(i,NIDXSMALL*sizeChannel);
	  }
	  for (Int_t j =0; j<sizeChannel;j++){
		  for (Int_t k=0; k<NIDXSMALL;k++){
			  tempArray->SetAt(i,j*NIDXSMALL+k,fTOFArray->GetArrayAt(i,j*NIDX+k));
		  }
	  }
  }
  tree->Fill();
  */	  
  
  AliInfo("Putting object for calibration in Reference data");

  TFile *fileout = TFile::Open("outputTOFCalibration.root","RECREATE");
  tempArray->Write();
  fileout->Close();
  delete tempArray;
}
//_____________________________________________________________________________

Bool_t AliTOFCalibTask::Select(AliESDtrack *t){

  // selecting good quality tracks
  //track selection: reconstrution to TOF:
  if ((t->GetStatus()&AliESDtrack::kTOFout)==0) {
    return 0;
  }
  fnESDkTOFout++;
  fhESDkTOFout->Fill(0);
  //IsStartedTimeIntegral
  if ((t->GetStatus()&AliESDtrack::kTIME)==0) {
    return 0;
  }
  fnESDkTIME++;
  fhESDkTIME->Fill(0);
  if (t->GetStatus() & AliESDtrack::kTRDbackup) { 
    Float_t xout= t->GetOuterParam()->GetX();
    if (xout<364.25 &&  xout > 300.) return 0;
  }
  fnESDTRDcut++;
  fhESDTRDcut->Fill(0);
  Double_t time=t->GetTOFsignal();	
  time*=1.E-3; // tof given in nanoseconds
  if(time >= fMinTime){
    return 0;
  }
  fnESDTIMEcut++;
  fhESDTIMEcut->Fill(0);
  
  Double_t mom=t->GetP();
  if (!(mom<=UPPERMOMBOUND && mom>=LOWERMOMBOUND)){
  //  return 0;  // skipping momentum cut
  } 
   
  Int_t assignedTOFcluster=t->GetTOFcluster();//index of the assigned TOF cluster, >0 ?
  if(assignedTOFcluster==-1){ // not matched
    return 0;
  }
  fnESDassTOFcl++;
  fhESDassTOFcl->Fill(0);
  AliDebug(2,"selecting the track");
  return 1;
}
//_____________________________________________________________________________

Int_t AliTOFCalibTask::SelectOnTime(Float_t *charray, Int_t ntracks, Int_t ich){

	// to be re-implemented with new object

  // discarding tracks with time-mintime < MINTIME
  
  Int_t ndeleted=0;
  Float_t mintime = 1E6;
  for (Int_t itr = 0;itr<ntracks;itr++){
    Int_t ientry=itr*NIDX;
    Float_t time = charray[ientry+DELTAIDXTIME];// in ns
      if (time<mintime) mintime = time;
  }
  AliDebug(1,Form("Mintime for channel %i = %f",ich, mintime));
  /*
  Int_t nTracksInChannel = fTOFArray->GetArray(ich)->GetSize();
  for (Int_t itr = 0;itr<ntracks;itr++){
    Int_t ientry=itr*NIDX;
    Float_t time = charray[ientry+DELTAIDXTIME];// in ns
    if ((time-mintime)>MINTIME) {
      ndeleted++;
      AliDebug(1,Form("Deleting %i track from channel %i, time = %f",ndeleted, ich, time));
      nTracksInChannel;
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
  */
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

  AliDebug(2,Form("number of sets = %i", nset));
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
    AliDebug(2,Form("set = %i, number of tracks in set = %i", i,ntrkinset));
      
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
    AliDebug(2,Form(" Set = %i with %i tracks ", i,ntrkinset));
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
    delete[] exptof; 
    delete[] timeofflight;
    delete[] texp;
    delete[] index;
    
  }
  return 1;
}
//---------------------------------------------------------------------
Float_t  AliTOFCalibTask::LoopCombPID(Int_t itrkinset, Int_t ntrkinset, Float_t **exptof, Float_t *texp, Float_t *timeofflight, Int_t *index, Float_t chisquarebest){

	// performing combinatorial PID in recursive way

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
      delete[] besttimezero;
      delete[] bestchisquare;
      delete[] bestweightedtimezero;
      delete[] weightedtimezero;
      delete[] timezero;
    }
  }
  return chisquarebest;
}
