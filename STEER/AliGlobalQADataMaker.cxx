/*
 The class for calculating the global (not detector specific) quality assurance.
 It reuses the following TLists from its base class 
    AliQADataMaker::fRecPointsQAList (for keeping the track residuals)
    AliQADataMaker::fESDsQAList      (for keeping global ESD QA data)
*/

#include <TH1F.h>

#include "AliGlobalQADataMaker.h"
#include "AliGeomManager.h"

ClassImp(AliGlobalQADataMaker)
 
void AliGlobalQADataMaker::InitRecPoints() {
  //------------------------------------------------------
  // This function fills the histograms of *track*residuals*
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
