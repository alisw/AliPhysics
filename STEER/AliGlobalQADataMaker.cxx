/*
 The class for calculating the global (not detector specific) quality assurance.
 It reuses the following TLists from its base class 
    AliQADataMaker::fRecPointsQAList (for keeping the track residuals)
    AliQADataMaker::fESDsQAList      (for keeping global ESD QA data)
*/

#include <TH1F.h>

#include "AliGlobalQADataMaker.h"
#include "AliGeomManager.h"
#include "AliESDEvent.h"

ClassImp(AliGlobalQADataMaker)
 
void AliGlobalQADataMaker::InitRecPoints() {
  //------------------------------------------------------
  // This function books the histograms of *track*residuals*
  // as a part of global QA
  //------------------------------------------------------
  Char_t *name[]={
    "SPD1 residuals Y","SPD1 residuals Z",
    "SPD2 residuals Y","SPD2 residuals Z",
    "SDD1 residuals Y","SDD1 residuals Z",
    "SDD2 residuals Y","SDD2 residuals Z",
    "SSD1 residuals Y","SSD1 residuals Z",
    "SSD2 residuals Y","SSD2 residuals Z",

    "TPC1 residuals Y","TPC1 residuals Z",
    "TPC2 residuals Y","TPC2 residuals Z",

    "TRD1 residuals Y","TRD1 residuals Z",
    "TRD2 residuals Y","TRD2 residuals Z",
    "TRD3 residuals Y","TRD3 residuals Z",
    "TRD4 residuals Y","TRD4 residuals Z",
    "TRD5 residuals Y","TRD5 residuals Z",
    "TRD6 residuals Y","TRD6 residuals Z",

    "TOF residuals Y","TOF residuals Z",

    "PHOS1 residuals Y","PHOS1 residuals Z",
    "PHOS2 residuals Y","PHOS2 residuals Z",

    "HMPID residuals Y","HMPID residuals Z",

    "MUON residuals Y","MUON residuals Z",

    "EMCAL residuals Y","EMCAL residuals Z"
  };

  for (Int_t m=1; m<AliGeomManager::kLastLayer; m++) {
    Int_t i=2*m-2;
    TH1F *h=new TH1F(name[i],name[i],100,-5.,5.);
    Add2RecPointsList(h,i);    
    h=new TH1F(name[i+1],name[i+1],100,-5.,5.);
    Add2RecPointsList(h,i+1);    
  }
}

void AliGlobalQADataMaker::InitESDs() {
  //------------------------------------------------------
  // This function books the ESD QA histograms
  // as a part of global QA
  //------------------------------------------------------
  Char_t *name[]={
    "Fraction of the assigned clusters in ITS",
    "Fraction of the assigned clusters in TPC",
    "Fraction of the assigned clusters in TRD"
  };
  Add2ESDsList(new TH1F(name[0],name[0],100,0.,2.),0);
  Add2ESDsList(new TH1F(name[1],name[1],100,0.,2.),1);
  Add2ESDsList(new TH1F(name[2],name[2],100,0.,2.),2);
}

void AliGlobalQADataMaker::MakeESDs(AliESDEvent * esd) {
  //-----------------------------------------------------------
  // This function fills the ESD QA histograms
  // as a part of global QA
  //-----------------------------------------------------------
  Int_t ntrk=esd->GetNumberOfTracks() ; 
  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *track=esd->GetTrack(i);

    if (track->IsOn(AliESDtrack::kITSrefit)) {
      Int_t n=track->GetITSclusters(0);
      GetESDsData(0)->Fill(Float_t(n)/6.); //6 is the number of ITS layers
    }

    if (track->IsOn(AliESDtrack::kTPCrefit)) {
      Int_t n =track->GetTPCNcls();
      Int_t nf=track->GetTPCNclsF();      // number of crossed TPC pad rows
      GetESDsData(1)->Fill(Float_t(n)/nf);
    }

    if (track->IsOn(AliESDtrack::kTRDrefit)) {
      Int_t n=track->GetTRDclusters(0);
      GetESDsData(2)->Fill(Float_t(n)/120.);//120 is the number of TRD time bins
    }

  }
}
