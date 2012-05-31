
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
  // Load libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  
  
  .x ~/NimStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");

  TFile f("CalibObjects.root");
  AliTPCcalibTrigger *calibTrigger = (AliTPCcalibTrigger *)f->Get("TPCCalib")->FindObject("calibTrigger");


  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chainTrack = tool.MakeChain("trigger.txt","Track",0,10200);
  chainTrack->Lookup();
  TChain * chainEvent = tool.MakeChain("trigger.txt","Event",0,10200);
  chainEvent->Lookup();


*/

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"

#include "AliTracker.h"
#include "AliMagF.h"
#include "AliTPCCalROC.h"

#include "AliLog.h"

#include "AliTPCcalibTrigger.h"

#include "TTreeStream.h"
#include "AliTPCTracklet.h"
#include "TTimeStamp.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibLaser.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"

ClassImp(AliTPCcalibTrigger)

AliTPCcalibTrigger::AliTPCcalibTrigger():
  AliTPCcalibBase("calibTrigger","calibTrigger"),
  fHisMap(0)
{

}

AliTPCcalibTrigger::AliTPCcalibTrigger(const char * name, const char * title):
  AliTPCcalibBase(name,title),
  fHisMap(0)
{
  //
  //
  //
  fHisMap = new TMap;
}

AliTPCcalibTrigger::~AliTPCcalibTrigger(){
  //
  // delete histograms
  // class is owner of all histograms
  //
  if (!fHisMap) return;
  fHisMap->SetOwner(kTRUE);
  fHisMap->DeleteAll();
  delete fHisMap;
}


Long64_t AliTPCcalibTrigger::Merge(TCollection *li) {
  //
  // Merge histograms
  //
  TIterator* iter = li->MakeIterator();
  AliTPCcalibTrigger* cal = 0;

  while ((cal = (AliTPCcalibTrigger*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibTrigger::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    TMap * addMap=(cal->fHisMap);
    if(!addMap) return 0;
    TIterator* iterator = addMap->MakeIterator();
    iterator->Reset();
    //    TPair* addPair=0;
    TObject *object=0;
    //
    while((object=iterator->Next())){
      THnSparse* his1 = dynamic_cast<THnSparseF*>(cal->fHisMap->GetValue(object->GetName()));
      if (!his1) continue;      
      his1->Print();
      THnSparse* his0 = dynamic_cast<THnSparseF*>(fHisMap->GetValue(object->GetName()));

      if(!his0){
        his0=MakeHisto(object->GetName());
        fHisMap->Add(new TObjString(object->GetName()),his0);
      }
      his0->Add(his1);
    }
  }
  return 0;
}



void AliTPCcalibTrigger::Process(AliESDEvent *event){
  //
  //
  //
  if (!event) return;
  const TString &trigger = event->GetFiredTriggerClasses();
  TTreeSRedirector * cstream =  GetDebugStreamer();
  //
  TObjString str(event->GetFiredTriggerClasses());
  Bool_t hasPIXEL=HasPIXEL(&str);
  Int_t  hasTRD=HasTRD(&str);
  Bool_t hasTOF=HasTOF(&str);
  Bool_t hasACORDE=HasACORDE(&str);
  //
  if (!GetHisto(trigger.Data())){
    AddHisto(trigger.Data(),MakeHisto(trigger.Data()));
  }
  if (!GetHisto("all")){
    AddHisto("all",MakeHisto("all"));
  }

  THnSparse *histoAll = GetHisto("all");
  THnSparse *histo = GetHisto(trigger.Data());
  Double_t xcont[9]={0,0,0,0,0,0,0,0,0};
  
  Int_t ntracks = event->GetNumberOfTracks();
  xcont[0] = ntracks;
  xcont[8] = 1;
  //
  // GetLongest track
  //  
  AliESDtrack * longest=0;
  Int_t nclmax=0;
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    AliESDtrack *track=event->GetTrack(itrack);
    if (!track) continue;
    if (track->GetTPCNcls()<=nclmax) continue;
    nclmax = track->GetTPCNcls();
    longest= track;
  }
  xcont[1]  =nclmax;
  histoAll->Fill(xcont);
  histo->Fill(xcont);
  if (cstream) {
    (*cstream) << "Event" <<
      "run="<<fRun<<
      "time="<<fTime<<
      "tname.="<<&str<<
      "pixel="<<hasPIXEL<<
      "trd="<<hasTRD<<
      "tof="<<hasTOF<<
      "acorde="<<hasACORDE<<
      "ntracks="<<ntracks<<
      "\n";
  }
  //
  xcont[8] = -1.;
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    AliESDtrack *track=event->GetTrack(itrack);
    if (!track) continue;
    Float_t dca[2];
    Double_t pxyz[3];
    track->GetDZ(0.,0.,0.,event->GetMagneticField(),dca);
    Bool_t status = track->GetPxPyPz(pxyz);
    Double_t alpha = TMath::ATan2(pxyz[1],pxyz[0]);
    xcont[1]=track->GetTPCNcls();
    xcont[2]=dca[0];
    xcont[3]=dca[1];
    xcont[4]=alpha;
    xcont[5]=track->GetParameter()[3];
    xcont[6]=track->Pt();
    xcont[7]=track->GetTPCsignal();    
    histoAll->Fill(xcont);
    histo->Fill(xcont);
    //
    //
    if (cstream) {
      Double_t mpt = track->GetParameter()[4];
      Int_t kokot[1000];
      Int_t nclITS=track->GetITSclusters(kokot);
      Int_t nclTPC=track->GetTPCNcls();
      Int_t nclTRD=track->GetTRDclusters(kokot);
      Int_t ntlTRD=track->GetTRDntracklets();
      ULong_t tstatus = track->GetStatus();
      (*cstream) << "Track" <<
	"run="<<fRun<<
	"time="<<fTime<<
	"tname.="<<&str<<
	"status="<<status<<	
	"tstatus="<<tstatus<<	
	//
	"ntracks="<<ntracks<<
	"tstatus="<<status<<
	"nclITS="<<nclITS<<
	"nclTPC="<<nclTPC<<
	"nclTRD="<<nclTRD<<
  "ntlTRD="<<ntlTRD<<
	//
	"pixel="<<hasPIXEL<<
	"trd="<<hasTRD<<
	"tof="<<hasTOF<<
	"acorde="<<hasACORDE<<
	"ncl="<<xcont[1]<<
	"dcaR="<<xcont[2]<<
	"dcaZ="<<xcont[3]<<
	"alpha="<<xcont[4]<<
	"theta="<<xcont[5]<<
	"pt="<<xcont[6]<<
	"dEdx="<<xcont[7]<<	
	"mpt="<<mpt<<
	"\n";
    }
  }
}

THnSparse *AliTPCcalibTrigger::MakeHisto(const char* trigger){
  //
  // Make event/track histograms
  // trigger histo 
  //
  //                 ntracks  nclMax dcaR dcaZ alpha   theta pt dEdx ev
  Int_t    bins[9] = {50,     40,    20,  20,  18,     25,   25, 25, 2 };
  //Int_t    bins[9] = {50*   20*    25*  25*  18*     25*   25* 25 };
  Double_t xmin[9] = {0.,     0,      0,  -250,   -3.14, -1.5,  0, 0, -1.};
  Double_t xmax[9] = {50,     160,  150,  250, 3.14,  1.5,   100, 100, 1.};
  TString  axisName[9]={
    "ntracks",
    "ncl",
    "dcaR",
    "dcaZ",
    "alpha",
    "theta",
    "pt",
    "dEdx",
    "ev"
  };
  TString  axisTitle[9]={
    "Number of tracks",
    "N_{cl}",
    "dca_{R} (cm)",
    "dca_{z} (cm)",
    "alpha (mrad)",
    "theta",
    "p_{t} (GeV/c)",
    "dEdx (a.u.)",
    "ev"
  };

  
  THnSparse *sparse = new THnSparseF(Form("his_%s",trigger), Form("his_%s",trigger), 9, bins, xmin, xmax);
  for (Int_t iaxis=0; iaxis<9; iaxis++){
    sparse->GetAxis(iaxis)->SetName(axisName[iaxis]);
    sparse->GetAxis(iaxis)->SetTitle(axisTitle[iaxis]);
  }
  return sparse;
}

THnSparse * AliTPCcalibTrigger::GetHisto(const char *trigger) { 
  //
  // return histogram for given class
  if (!fHisMap) fHisMap=new TMap;
  return (THnSparse*) fHisMap->GetValue(trigger);
}

void   AliTPCcalibTrigger::AddHisto(const char *trigger, THnSparse *his) { 
  if (!GetHisto(trigger)) {
    TObjString *pstr = new TObjString(trigger);
    fHisMap->Add(pstr,his);
  }
}

TTree * AliTPCcalibTrigger::MakeTree(const char * fname){
  //
  //
  //
  TTreeSRedirector * sred = new TTreeSRedirector(fname);
  TTreeStream &pcstream = (*sred)<<"Trigger";
  //
  //
  TIterator* iterator = fHisMap->MakeIterator();
  TObject * object=0;
  //
  while((object=iterator->Next())){
    MakeTree(pcstream, object->GetName());
  }
  delete sred;
  TFile *f = new TFile(fname);
  TTree *tree = (TTree*)f->Get("Trigger");
  return tree;
}


void AliTPCcalibTrigger::MakeTree(TTreeStream &pcstream, const char *tname){
  //
  //  TTreeSRedirector * sred = new TTreeSRedirector("trigger.root");
  //  TTreeStream &pcstream = (*sred)<<"Trigger";
  //
  //AliTPCcalibTrigger *calibTrigger = this;
  Double_t value;
  THnSparse * his = GetHisto(tname);
  if (!his) return;
  //
  Int_t bins[1000];
  Int_t ndim = his->GetNdimensions();
  Double_t position[10];
  //
  TObjString str(tname);
  Bool_t isAll  = str.String().Contains("all");
  Bool_t hasPIXEL=HasPIXEL(&str);
  Int_t  hasTRD=HasTRD(&str);
  Bool_t hasTOF=HasTOF(&str);
  Bool_t hasACORDE=HasACORDE(&str);
  for (Long64_t i = 0; i < his->GetNbins(); ++i) {
    value = his->GetBinContent(i, bins);
    pcstream<<"val="<<value;
    pcstream<<"tname.="<<&str;
    //
    pcstream<<"all="<<isAll;
    pcstream<<"pixel="<<hasPIXEL;
    pcstream<<"trd="<<hasTRD;
    pcstream<<"tof="<<hasTOF;
    pcstream<<"acorde="<<hasACORDE;
    //
    for (Int_t idim = 0; idim < ndim; idim++) {
      position[idim] = his->GetAxis(idim)->GetBinCenter(bins[idim]);
      pcstream<<Form("%s=",his->GetAxis(idim)->GetName())<<position[idim];
    }
    pcstream<<"\n";
  }
}


Bool_t AliTPCcalibTrigger::HasTOF(TObjString *tname){
  //
  Bool_t result = kFALSE;
  result|=(tname->String().Contains("0OB")>0);
  result|=(tname->String().Contains("0OC")>0);
  return result;
}

Bool_t AliTPCcalibTrigger::HasACORDE(TObjString *tname){
  Bool_t result = kFALSE;
  result|=(tname->String().Contains("0ASL")>0);
  result|=(tname->String().Contains("0AMU")>0);
  result|=(tname->String().Contains("0ASC")>0);
  return result;
}

Bool_t AliTPCcalibTrigger::HasPIXEL(TObjString *tname){
  return (tname->String().Contains("0SCO")>0);
}

Int_t AliTPCcalibTrigger::HasTRD(TObjString *tname){
  //
  // Returns a mask containing TRD trigger information
  // 0: No TRD trigger fired
  // 1: TRD L1 fired
  // 2: TRD L0 (krypton trigger) fired
  //
  Int_t result = 0;
  if(tname->String().Contains("TRD")) result = 1;     // Normal TRD L1 name
  if(tname->String().Contains("0HPT1")) result = 1;   // Old TRD L1 name
  if(tname->String().Contains("0HWU") && !tname->String().Contains("TRD")) result = 2;  // pretrigger always input for L1
  return result;
}

