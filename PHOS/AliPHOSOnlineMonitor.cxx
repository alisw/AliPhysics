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

/* $Id$ */

//_________________________________________________________________________
// Class intended to perform online monitoring of PHOS beamtests
// Being constructed, produces menu with list of available histograms to fill
// Once histograms are selected, button "Go" should be pressed to start scan of data.
// Prepared histograms will be periodically updated during scan of the data.
// Note:
// 1. To plot most of histograms, a "Connection Table", relating ADC signal index and AbsId 
//    of PHOS crystal, should be created beforehand. To do this, call macro
//       $ALICE_ROOT/PHOS/macros/BeamTest/MakeConTableDB.C 
//    with apropriate number of raws and columns of prototype.
// 2. To perform reconstruction and e.g. invariant mass analysis, a "Calibration Database"
//    should be created beforehand. To do this, call macro 
//       $ALICE_ROOT/PHOS/macros/BeamTest/MakeConTableDB.C 
//    to read calibration parameters from file or use AliPHOSCalibrator to calculate 
//    pedestal and gains.
// 3. Once histograms are filled with "Go" method, they can be written to file
//    with WriteHisto("Filename.root") method. 
//
//*-- Author : D.Peressounko (RRC KI) after A.V. Kuryakin, (Sarov) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TROOT.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGroupButton.h"
#include "TFile.h"
#include "TSystem.h"


// --- Standard library ---
#include "TBenchmark.h"
#include "Riostream.h"

// --- AliRoot header files ---
#include "AliPHOSOnlineMonitor.h"
#include "AliPHOSConTableDB.h"
#include "AliPHOSGeometry.h"
#include "AliRawReaderDateV3.h"
#include "AliRawEventHeaderBase.h"
#include "AliPHOSRawStream2004.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGetterLight.h"  
#include "AliPHOSClusterizerv1.h"  
#include "AliPHOSTrackSegmentMakerv1.h"  
#include "AliPHOSPIDv1.h"  
#include "AliPHOSCalibrManager.h" 
#include "AliPHOSCalibrationDB.h"

ClassImp(AliPHOSOnlineMonitor)


//____________________________________________________________________________ 
AliPHOSOnlineMonitor::AliPHOSOnlineMonitor(): 
  TDialogCanvas("PHOS","PHOS",150,300),
  fScanPed(kFALSE),
  fScanSig(kFALSE),
  fReconstruct(kFALSE),
  fNevents(0),
  fNUpdate(1000),
  fCanvasList(new TList),
  fHistosList(new TList),
  fInputFile(),
  fGeom(0),
  fcdb(0)
{
  MakeButtons() ;
  Modified(kTRUE);
  Update();
  SetEditable(kFALSE);
  
  //add this TFitPanel to the list of cleanups such that in case
  //the referenced object is deleted, its pointer be reset
  gROOT->GetListOfCleanups()->Add(this);
  
  fRefObject = this;
  fRefPad    = (TPad*)gROOT->GetSelectedPad();

  fGeom = AliPHOSGeometry::GetInstance("IHEP","") ;
}
//____________________________________________________________________________ 
AliPHOSOnlineMonitor::AliPHOSOnlineMonitor(const char * inputfile): 
  TDialogCanvas("PHOS","PHOS",150,300),
  fScanPed(kFALSE),
  fScanSig(kFALSE),
  fReconstruct(kFALSE),
  fNevents(0),
  fNUpdate(1000),
  fCanvasList(new TList),
  fHistosList(new TList),
  fInputFile(inputfile),
  fGeom(0),
  fcdb(0)
{
  MakeButtons() ;
  Modified(kTRUE);
  Update();
  SetEditable(kFALSE);
  
  //add this TFitPanel to the list of cleanups such that in case
  //the referenced object is deleted, its pointer be reset
  gROOT->GetListOfCleanups()->Add(this);
  
  fRefObject = this;
  fRefPad    = (TPad*)gROOT->GetSelectedPad();

  fGeom = AliPHOSGeometry::GetInstance("IHEP","") ;
}

//____________________________________________________________________________ 
AliPHOSOnlineMonitor::AliPHOSOnlineMonitor(const AliPHOSOnlineMonitor & /*rhs*/):
  TDialogCanvas(),
  fScanPed(kFALSE),
  fScanSig(kFALSE),
  fReconstruct(kFALSE),
  fNevents(0),
  fNUpdate(0),
  fCanvasList(0),
  fHistosList(0),
  fInputFile(),
  fGeom(0),
  fcdb(0)
{
  Fatal("AliPHOSOnlineMonitor", "not implemented");
}

//____________________________________________________________________________ 
AliPHOSOnlineMonitor & AliPHOSOnlineMonitor::operator = (const AliPHOSOnlineMonitor &)
{
  Fatal("operator = ", "not implemented");
  return *this;
}

//____________________________________________________________________________ 
AliPHOSOnlineMonitor::~AliPHOSOnlineMonitor()
{
  //Obvious, but unevoidable comment for destructor: cleans up everething.
  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  while((c=(TCanvas*)nextCanvas()))
    delete c ;
  delete fCanvasList ;
  
  TIter nextHisto(fHistosList);
  TH1D * h ;
  while((h=(TH1D*)nextHisto()))
    delete h ;
  delete fHistosList ; 

  if(fcdb)
    delete fcdb ;

}
//____________________________________________________________________________ 
void AliPHOSOnlineMonitor::MakeButtons(void){
  //Make buttons on graphical menu
  Int_t nButtons = 16;
  TGroupButton * b ;
  Float_t xmin = 0.0;
  Float_t ymin = 0.01;
  Float_t xmax = 0.99;
  Float_t ymax = 0.99;
  Float_t dy = (ymax-ymin)/nButtons ;

  Float_t y2=ymax ;
  Float_t y1=y2-dy ;
  b = new TGroupButton("APPLY","Triggers","",xmin,y1,xmax,y2);
  b->Draw();
  y2=y1 ;
  y1=y1-dy ;
  b = new TGroupButton("APPLY","Pedestals","",xmin,y1,xmax,y2);
  b->Draw();
  y2=y1 ;
  y1=y1-dy ;
  b = new TGroupButton("APPLY","Spectrum all","",xmin,y1,xmax,y2);
  b->Draw();
  y2=y1 ;
  y1=y1-dy ;
  b = new TGroupButton("APPLY","Spectrum g","",xmin,y1,xmax,y2);
  b->Draw();
  y2=y1 ;
  y1=y1-dy ;
  b = new TGroupButton("APPLY","Inv Mass","",xmin,y1,xmax,y2);
  b->Draw();
  for(Int_t i=1; i<=5; i++){
    y2=y1 ;
    y1=y1-dy ;
    char name[10] ;
    sprintf(name,"Edep(ADC) %d",i) ;
    b = new TGroupButton("APPLY",name,"",xmin,y1,xmax,y2);
    b->Draw();
  }
  for(Int_t i=1; i<=5; i++){
    y2=y1 ;
    y1=y1-dy ;
    char name[10] ;
    sprintf(name,"Edep(Cal) %d",i) ;
    b = new TGroupButton("APPLY",name,"",xmin,y1,xmax,y2);
    b->Draw();
  }
  y2=y1 ;
  y1=y1-dy ;
  b = new TGroupButton("APPLY","Go","",xmin,y1,xmax,y2);
  b->SetTextColor(2);
  b->Draw();
}
//____________________________________________________________________________ 
void AliPHOSOnlineMonitor::Apply(const char *action){
  //Function to handle button actions

  TDialogCanvas::Apply() ;

  TObject *obj;
  TGroupButton *button;
  TIter next(fPrimitives);
  
  if (!strcmp(action,"Triggers")) {
    DrawTriggers() ;
  }

  if (!strcmp(action,"Pedestals")) {
    DrawPedestals() ;
  }
  if (!strcmp(action,"Spectrum all")) {
    DrawSpectrum("all") ;
  }
  if (!strcmp(action,"Spectrum g")) {
    DrawSpectrum("gamma") ;
  }
  if (!strcmp(action,"Inv Mass")) {
    DrawMinv() ;
  }
  if(strstr(action,"Edep")){
    Int_t n ;
    char tmp[10] ;
    sscanf(action,"%s %d",tmp,&n) ;
    char opt[5]="" ;
    if(strstr(action,"Cal"))
      sprintf(opt,"Cal") ;
    DrawEdep(n,opt) ;
  }
  if (!strcmp(action,"Go")) {
    Go() ;
  }

  //Mark button as pressed
  if(strcmp(action,"Go")){ //Do not mark "Go" button
    while ((obj = next())) {
      if (obj->InheritsFrom(TGroupButton::Class())) {
	button = (TGroupButton*)obj;
	if(!strcmp(button->GetTitle(),action)){
	  if (button->GetBorderMode() > 0){
	    button->SetBorderMode(-1) ;
	    button->Modified(kTRUE);
	  }
	}
      }
    }
  }
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::SetInputFile(const char * filename){
  //close previously opened
  
  fInputFile = filename ;
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::DrawPedestals(){
  //Prepare canvas and histograms for drawing pedestals

  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  Bool_t exists = kFALSE ;
  while((c=(TCanvas*)nextCanvas())){
    if(!strcmp(c->GetName(),"Pedestals")){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    c = new TDialogCanvas("Pedestals","Pedestals",300,200) ;
    fCanvasList->AddLast(c) ;
  }
  
  TIter nextHisto(fHistosList);
  TH1D * h ;
  exists = kFALSE ;
  while((h=(TH1D*)nextHisto())){
    if(!strcmp(h->GetName(),"hPedestals")){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    h = new TH1D("hPedestals","Pedestals per event",fGeom->GetNModules()*fGeom->GetNCristalsInModule(),0.,
		 1.*fGeom->GetNModules()*fGeom->GetNCristalsInModule()) ;
    fHistosList->AddLast(h) ;
  }
  
  c->cd() ;
  h->SetStats(0) ;
  h->Draw() ;
  fScanPed = kTRUE ; //We will scan pedestals

}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::DrawTriggers(){
  //Prepare canvas and histogram for drawing triggers

  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  Bool_t exists = kFALSE ;
  while((c=(TCanvas*)nextCanvas())){
    if(!strcmp(c->GetName(),"Triggers")){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    c = new TDialogCanvas("Triggers","Triggers",200,200) ;
    fCanvasList->AddLast(c) ;
  }
  
  TIter nextHisto(fHistosList);
  TH1D * h ;
  exists = kFALSE ;
  while((h=(TH1D*)nextHisto())){
    if(!strcmp(h->GetName(),"hTriggers")){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    h = new TH1D("hTriggers","Triggers",2,0.,2.) ;
    fHistosList->AddLast(h) ;
  }
  //Make Labels
  h->SetBit(TH1::kCanRebin);  
  h->Fill("LED",0.0000001) ; 
  h->Fill("PUL",0.0000001) ; 
  h->Fill("PED",0.0000001) ; 
  h->Fill("NEL",0.0000001) ; 
  h->Fill("WEL",0.0000001) ; 
  h->Fill("SOB",0.0000001) ; 
  h->Fill("EOB",0.0000001) ; 
  h->Fill("wrong",0.0000001) ;
  h->LabelsOption("h");
  h->LabelsDeflate();
  h->SetStats(0) ;
  c->cd() ;
  h->Draw() ;
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::DrawSpectrum(const char * opt){
  //Prepare canvas and histograms for drawing spectra of all reconstructed particles or photons

  TString name("Spectrum") ;
  name+=opt ;
  
  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  Bool_t exists = kFALSE ;
  while((c=(TCanvas*)nextCanvas())){
    if(!strcmp(c->GetName(),name.Data())){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    c = new TDialogCanvas(name,name,250,300) ;
    fCanvasList->AddLast(c) ;
  }
  
  TIter nextHisto(fHistosList);
  TH1D * h ;
  exists = kFALSE ;
  name.Prepend("h") ;
  while((h=(TH1D*)nextHisto())){
    if(!strcmp(h->GetName(),name.Data())){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    h = new TH1D(name,name,100,0.,100.) ;
    fHistosList->AddLast(h) ;
  }
  
  h->SetStats(0) ;
  c->cd() ;
  h->Draw() ;
  fReconstruct = kTRUE ;
  fScanSig = kTRUE ; //We will scan pedestals

}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::DrawMinv(){
  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  Bool_t exists = kFALSE ;
  while((c=(TCanvas*)nextCanvas())){
    if(!strcmp(c->GetName(),"InvMass")){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    c = new TDialogCanvas("InvMass","Invariant mass",300,200) ;
    fCanvasList->AddLast(c) ;
  }
  
  TIter nextHisto(fHistosList);
  TH1D * h ;
  exists = kFALSE ;
  while((h=(TH1D*)nextHisto())){
    if(!strcmp(h->GetName(),"hInvMass")){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    h = new TH1D("hInvMass","hInvMass",1000,0.,1.0) ;
    fHistosList->AddLast(h) ;
  }
  
  c->cd() ;
  h->Draw() ;
  h->SetStats(0) ;
  fReconstruct = kTRUE ;
  fScanSig = kTRUE ; //We will scan pedestals
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::DrawEdep(Int_t mod,const char * opt){
  char name[15] ;
  sprintf(name,"Edep%s %d",opt,mod) ;

  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  Bool_t exists = kFALSE ;
  while((c=(TCanvas*)nextCanvas())){
    if(!strcmp(c->GetName(),name)){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    c = new TDialogCanvas(name,name,300,200) ;
    fCanvasList->AddLast(c) ;
  }
  
  TIter nextHisto(fHistosList);
  TH2D * h ;
  exists = kFALSE ;
  sprintf(name,"hEdep%s%d",opt,mod) ;
  while((h=(TH2D*)nextHisto())){
    if(!strcmp(h->GetName(),name)){
      exists = kTRUE ;
      break;
    }
  }
  if(!exists){
    h = new TH2D(name,name,fGeom->GetNPhi(),0.,1.*fGeom->GetNPhi(),fGeom->GetNZ(),0.,1.*fGeom->GetNZ()) ;
    fHistosList->AddLast(h) ;
  }
  
  c->cd() ;
  h->Draw("col") ;
  h->SetStats(0) ;
  fScanSig = kTRUE ; //We will scan signal events
  if(strstr(opt,"Cal"))
    fReconstruct = kTRUE ;
 
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::ScanPedestals(TClonesArray * digits){
  //This method is called for events with PED trigger
  //We fill bins with ADC values

  TH1D * h = (TH1D*)gROOT->FindObjectAny("hPedestals");
  if(!h){
    Error("ScanPedestals","Can not fild histogram hPedestals") ;
    return ;
  }
  for(Int_t i=0; i<digits->GetEntriesFast(); i++){
    AliPHOSDigit * dig = static_cast<AliPHOSDigit*>(digits->At(i)) ;
    h->AddBinContent(dig->GetId(),dig->GetAmp()) ;
  }
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::ScanEdep(TClonesArray * digits){
  //Fill 2D distribution of ADC values in NEL and WEL events
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ;
  AliPHOSCalibrationDB *cdb = 0 ;
  if(gime)
    cdb = gime->CalibrationDB() ;
  Int_t mod = 0 ;
  char name[15] ;
  TH2D * h = 0 ;
  TH2D * hCal = 0 ;
  for(Int_t i=0; i<digits->GetEntriesFast(); i++){
    AliPHOSDigit * dig = static_cast<AliPHOSDigit*>(digits->At(i)) ;
    Int_t relId[4] ;
    fGeom->AbsToRelNumbering(dig->GetId(),relId) ;
    if(mod != relId[0]){ //new module, look for histograms
      mod = relId[0] ;
      sprintf(name,"hEdep%d",mod) ;
      h = (TH2D*)gROOT->FindObjectAny(name);
      sprintf(name,"hEdepCal%d",mod) ;
      hCal = (TH2D*)gROOT->FindObjectAny(name);
    }
    if(h)
      h->Fill(relId[2]-0.1,relId[3]-0.1,1.*dig->GetAmp()) ;
    if(hCal)
      hCal->Fill(relId[2]-0.1,relId[3]-0.1,cdb->Calibrate(dig->GetAmp(),dig->GetId())) ;
  }
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::ScanRecon(TClonesArray * recParticles){
  if(!recParticles || recParticles->GetEntries()==0) return ;

  TH1D* hSpectr   = (TH1D*)gROOT->FindObjectAny("hSpectrumall");
  TH1D* hSpectrGam= (TH1D*)gROOT->FindObjectAny("hSpectrumgamma");
  TH1D* hInvMass  = (TH1D*)gROOT->FindObjectAny("hInvMass");
  for(Int_t i=0; i<recParticles->GetEntriesFast() ; i++){
    AliPHOSRecParticle * p = (AliPHOSRecParticle *)recParticles->At(i) ;
    if(hSpectr)hSpectr->Fill(p->Energy()) ;
    if(hSpectrGam && p->IsPhoton())hSpectrGam->Fill(p->Energy()) ;
    if(hInvMass){
      for(Int_t j=i+1; j<recParticles->GetEntriesFast() ; j++){
	AliPHOSRecParticle * p2 = (AliPHOSRecParticle *)recParticles->At(j) ;
	Double_t e = p->Energy() + p2->Energy() ;
	Double_t x = p->Px() + p2->Px() ;
	Double_t y = p->Py() + p2->Py() ;
	Double_t z = p->Pz() + p2->Pz() ;
	Double_t m = e*e-x*x-y*y-z*z ;
	hInvMass->Fill(m>0?TMath::Sqrt(m): 0. ) ;
      }
    }
  }
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::ScanTrigger(Int_t trig){
  //Fills trigger distribution

  TH1D * h = (TH1D*)gROOT->FindObjectAny("hTriggers");
  if(!h) return ;
  switch(trig){
  case AliPHOSRawStream2004::kLED : h->Fill("LED",1.) ; break ;
  case AliPHOSRawStream2004::kPUL : h->Fill("PUL",1.) ; break ;
  case AliPHOSRawStream2004::kPED : h->Fill("PED",1.) ; break ;
  case AliPHOSRawStream2004::kNEL : h->Fill("NEL",1.) ; break ;
  case AliPHOSRawStream2004::kWEL : h->Fill("WEL",1.) ; break ;
  case AliPHOSRawStream2004::kSOB : h->Fill("SOB",1.) ; break ;
  case AliPHOSRawStream2004::kEOB : h->Fill("EOB",1.) ; break ;
  default : h->Fill("wrong",1.) ;
  }
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::SetConTableDB(const char * filename){
  //Read ConnectionTableDB from file
  TFile * file = new TFile(filename) ;
  AliPHOSConTableDB * tmp = (AliPHOSConTableDB*)file->Get("AliPHOSConTableDB") ;
  fcdb = new AliPHOSConTableDB(*tmp) ;
  file->Close() ; 
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::Go(){
  //Perform scan of curent event
  gBenchmark->Start("PHOSOnlineMon"); 

  //First test if we need "Connection table" then open it
  if(!fcdb){
    SetConTableDB() ;
    if(fcdb){
      Info("Go","Read Connection table from file \"ConTableDB.root\"") ;
    }else{
      Error("Go","Please, set connection table with SetConTableDB() method") ;
      return ;
    }
  }

  AliPHOSGetterLight * gime = AliPHOSGetterLight::Instance("PHOS","On Flight") ;

  //Configure CalibrManager to read data from file
  //Create calibration database and read it
  AliPHOSCalibrationDB * calibDB = 0 ;
  if(fScanSig || fReconstruct){ //We will ned calibration parameters
    AliPHOSCalibrManager::GetInstance("CalibrDB.root","root") ;
    //If we configured manager to read from ASCII file, 
    //give him connection table. OK, it will not harm in any case.
    AliPHOSCalibrManager::GetInstance()->SetConTable(fcdb) ;
    
    calibDB = new AliPHOSCalibrationDB("OnLine") ;
    calibDB->GetParameters() ; //Read parameters using Manager
    gime->SetCalibrationDB(calibDB) ;
  }
  
  //Now open data file
  AliRawReaderDateV3 *rawReader = new AliRawReaderDateV3(fInputFile) ; 
  rawReader->RequireHeader(kFALSE);
  AliPHOSRawStream2004     *rawStream = new AliPHOSRawStream2004(rawReader) ;
  rawStream->SetConTableDB(fcdb) ;
  
  TClonesArray * digits = gime->Digits() ;
  TClonesArray * recParticles = gime->RecParticles() ;
  AliPHOSClusterizerv1* clu = 0 ;
  AliPHOSTrackSegmentMakerv1 * tsm = 0 ;
  AliPHOSPIDv1 * pid = 0 ;
  if(fReconstruct){ //We will need calibation parameters
    clu = new AliPHOSClusterizerv1(gime->PHOSGeometry()) ;
    clu->SetWriting(0) ; //Do not write to file
    clu->SetEmcMinE(0.05) ;  //Minimal energy of the digit
    clu->SetEmcLocalMaxCut(0.05) ; //Height of local maximum over environment
    clu->SetEmcClusteringThreshold(0.2) ; //Minimal energy to start cluster
//    clu->SetUnfolding(kFALSE) ; //Do not unfold
    tsm = new AliPHOSTrackSegmentMakerv1(gime->PHOSGeometry());
    tsm->SetWriting(0) ; //Do not write to file
    pid = new AliPHOSPIDv1(gime->PHOSGeometry()) ;
    pid->SetWriting(0) ; //Do not write to file    
  }
  
  fNevents=0 ;
  //Scan all event in file
  printf("      ") ;
  while(rawReader->NextEvent()){
    //Is it PHYSICAL event
    if(rawReader->GetType() == AliRawEventHeaderBase::kPhysicsEvent){
      fNevents++ ;
      if(fNevents%100 ==0){
	printf("\b\b\b\b\b\b%6d",fNevents) ;
      }
      if(rawStream->ReadDigits(digits)){
	
	//Test trigger
	//Pedestal Event
	ScanTrigger(rawStream->GetTrigger()) ;
	if(rawStream->IsPEDevent() && fScanPed){
	  ScanPedestals(digits) ;
	}
	if((rawStream->IsNELevent() || rawStream->IsWELevent()) && fScanSig){
	  ScanEdep(digits) ;
	  if(fReconstruct){
	    // PLEASE FIX IT !!!
	    //	    gime->Clusterizer()->Exec("") ;
	    //	    gime->TrackSegmentMaker()->Exec("") ;
	    //	    gime->PID()->Exec("") ;
	    ScanRecon(recParticles) ;
	  }
	}
      }
      
      if(fNevents%fNUpdate == 0 ){ //upate all histograms	
	TIter nextCanvas(fCanvasList);
	TCanvas * c ;
	while((c=(TCanvas*)nextCanvas())){
	  c->Modified() ;
	  c->Update() ;
	}
      }
      gSystem->ProcessEvents(); 
    }
    //    if(fNevents>=200)break ;
  }
  printf("\n") ;
  gBenchmark->Stop("PHOSOnlineMon");
  Float_t time = gBenchmark->GetCpuTime("PHOSOnlineMon") ;
  printf("took %f seconds for scanning, i.e. %f seconds per event %d  \n",
	 time,time/fNevents,fNevents) ; 
  
  //Update canvas with histograms at the end
  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  while((c=(TCanvas*)nextCanvas())){
    c->Modified(kTRUE) ;
  }
  
  if(clu)delete clu ;
  if(tsm)delete tsm ;
  if(pid)delete pid ;
  printf("delete 1 \n") ;
  if(calibDB) delete calibDB ;
  delete rawStream ;
  delete rawReader ;
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::Clean(){
  //Cleans content of all histograms

  TIter nextHisto(fHistosList);
  TH1D * h ;
  while((h=(TH1D*)nextHisto())){
    h->Reset("ISE") ;
  }
  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  while((c=(TCanvas*)nextCanvas())){
    c->Modified() ;
  }
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::Reset(){
  //delets all canvas and histograms,
  //marks buttons as unpressed

  TIter nextHisto(fHistosList);
  TH1D * h ;
  while((h=(TH1D*)nextHisto())){
    fHistosList->Remove(h) ;
    delete h ;
  }
  TIter nextCanvas(fCanvasList);
  TCanvas * c ;
  while((c=(TCanvas*)nextCanvas())){
    fCanvasList->Remove(c) ;
    delete c ;
  }
  TObject *obj;
  TGroupButton *button;
  TIter next(fPrimitives);
  
  //Mark buttons as anpressed
  while ((obj = next())) {
    if (obj->InheritsFrom(TGroupButton::Class())) {
      button = (TGroupButton*)obj;
      if (button->GetBorderMode() < 0){
	button->SetBorderMode(1) ;
	button->Modified(kTRUE);
      }
    }
  }
  
}
//____________________________________________________________________________ 
void  AliPHOSOnlineMonitor::WriteHistograms(const char * filename){
  //Write filled histograms to file
  TFile * file = new TFile(filename,"Update") ;
  file->cd() ;
  TIter nextHisto(fHistosList);
  TH1 * h ;
  while((h=(TH1*)nextHisto())){
    h->Write(0,TObject::kOverwrite) ;
  }
  file->Close() ; 
}
