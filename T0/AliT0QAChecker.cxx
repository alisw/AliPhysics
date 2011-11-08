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
//...
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for T0
//---------------------------------------------
//checkig without reference data:
//for RAW QA all histograms should have approximatly the same 
//number of entries as RefPoint
//for Rec Points checks 
//  - amplitude measured by 2 methos
//  - online and offline T0 measurements
// for ESD quality of reconstruction ( and measurements):
// RMS of vertex and T0 less than 75ps
//
// Alla.Maevskaya@cern.ch   
//...

// --- ROOT system ---
#include <Riostream.h>
#include <TClass.h>
#include <TH1F.h> 
#include <TF1.h> 
#include <TFitResultPtr.h>
#include <TH2.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TMath.h>
#include <TString.h>
#include <TPaveText.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"
#include "AliT0QAChecker.h"

ClassImp(AliT0QAChecker)
//____________________________________________________________________________
AliT0QAChecker::AliT0QAChecker() :
AliQACheckerBase("T0","T0 Quality Assurance Checker")

{
  // Standard constructor

}

//____________________________________________________________________________
AliT0QAChecker::AliT0QAChecker(const AliT0QAChecker& qac):
  AliQACheckerBase(qac.GetName(), qac.GetTitle()) 
{
  // copy constructor
  AliError("Copy should not be used with this class\n");
}
//____________________________________________________________________________
AliT0QAChecker& AliT0QAChecker::operator=(const AliT0QAChecker& qac){
  // assignment operator
  this->~AliT0QAChecker();
  new(this)AliT0QAChecker(qac);
  return *this;
}


//____________________________________________________________________________
AliT0QAChecker::~AliT0QAChecker(){
  // destructor

}

//__________________________________________________________________
void AliT0QAChecker::Check(Double_t *  test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/)
{

  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!
    
 char * detOCDBDir = Form("T0/%s/%s", AliQAv1::GetRefOCDBDirName(), AliQAv1::GetRefDataDirName()) ; 

  AliCDBEntry *QARefRec = AliQAManager::QAManager()->Get(detOCDBDir);
  //  QARefRec->Dump();
  if( !QARefRec){
    AliInfo("QA reference data NOT retrieved for Reconstruction check. No T0 reference distribution");
  }
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    test[specie]    = 10.0 ; 

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    //  TString dataType = AliQAv1::GetAliTaskName(index);
    if (!(AliQAv1::Instance()->IsEventSpecieSet(specie) && list[specie]) || list[specie]->GetEntries() == 0) {
      test[specie] = 1. ; // nothing to check
      continue;
    }
    if (index == AliQAv1::kRAW && AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib)
      //      if (index == AliQAv1::kRAW )
      {
	test[specie] = CheckRaw(list[specie]);
      }
    if (index == AliQAv1::kESD && AliRecoParam::Convert(specie) != AliRecoParam::kCalib)
      test[specie] = CheckESD(list[specie]);
    //
  }
}

//--------------------------------------------------------------------------
Double_t AliT0QAChecker::CheckRaw(TObjArray *listrec) const
{
   Float_t checkr = 0;
   
  TString hname[10];
   TH1*hdata; 
   const char *cname;
   TH1 *fhRawEff[10];
   Int_t nh=0;
 

   //   Int_t nnn[4] = { 420, 458, 459, 460};
   Int_t nnn[4] = { 169, 207, 208, 209};
   for (Int_t ir=0; ir<4; ir++)
     {
	 hdata = (TH1*) listrec->UncheckedAt(nnn[ir]);
	 if(hdata) {
	   cname = hdata->GetName();
	   hname[ir] = cname;
	   fhRawEff[nh] = hdata;
	   nh++;
	 }
    }
     
   TLine linelowyellow(0, 0.7, 24, 0.7);    
   linelowyellow.SetLineColor(5);
   linelowyellow.SetLineStyle(3);
   linelowyellow.SetLineWidth(4);
   TLine linelowred(0, 0.2, 24, 0.2);    
   linelowred.SetLineColor(2);
   linelowred.SetLineStyle(3);
   linelowred.SetLineWidth(4);

   Float_t thryell = 0.7;
   //   TPaveText text(0.30,0.50,0.99,0.99,"NDC");    
   Float_t thrred = 0.2;
   Float_t chcont =0;
    for (Int_t ih= 0; ih<4; ih++)
     { 
       // clean objects added at previous checks
       TList* lstF = fhRawEff[ih]->GetListOfFunctions();
       if (lstF) {
	 TObject *stats = lstF->FindObject("stats");
	 lstF->Remove(stats);
	 TObject *obj;
	 while ((obj = lstF->First())) {
	   while(lstF->Remove(obj)) { }
	   delete obj;
	 }
	 if (stats) lstF->Add(stats);
       } 
       fhRawEff[ih]->SetLineWidth(2);
       fhRawEff[ih]->SetMaximum(2.);
       fhRawEff[ih]->SetMinimum(0.);
       fhRawEff[ih]->GetListOfFunctions()->Add((TLine*)linelowyellow.Clone());
       fhRawEff[ih]->GetListOfFunctions()->Add((TLine*)linelowred.Clone());

       Int_t nbins= fhRawEff[ih]->GetNbinsX();
       Bool_t yell = kFALSE;
       Bool_t red = kFALSE;
       for (Int_t in=1; in<nbins-1; in++)
	 {
	   if(ih==0 && in==5) continue;
	   chcont=fhRawEff[ih]->GetBinContent(in);
	   if (chcont < thryell  ) yell = kTRUE;
	   if (chcont < thrred  ) red = kTRUE;
	 }
       
       if (! yell && !red) {
	 AliDebug(AliQAv1::GetQADebugLevel(), Form(" efficiency in all channes %s  is good", fhRawEff[ih]->GetName() ));
	 checkr=1.;
	 //	 text.AddText(Form("T0 RUN %d ",AliCDBManager::Instance()->GetRun()));
	 //	 text.AddText(Form(" No problems "));
	 //	 text.SetFillColor(3);
       }
       
       if(red ) {
	 checkr = 0.;
	 AliDebug(AliQAv1::GetQADebugLevel(), Form(" efficiency in all channes %s  is not so good", fhRawEff[ih]->GetName() ));
	 //	 text.AddText(Form("T0 RUN %d ",AliCDBManager::Instance()->GetRun()));
	 //	 text.AddText(Form("Very serious problem, call expert "));
	 //	 text.SetFillColor(2);
       }
       
       if ( yell && !red) {
	 AliDebug(AliQAv1::GetQADebugLevel(), Form(" efficiency in all channes %s  is not so good", fhRawEff[ih]->GetName() ));
	 checkr=0.75;
	 //	 text.AddText(Form("T0 RUN %d ",AliCDBManager::Instance()->GetRun()));
	 //	 text.AddText(Form("Some problems "));
	 //	 text.SetFillColor(5);
      }
       // fhRawEff[ih]->GetListOfFunctions()->Add((TPaveText*)text.Clone());	       
     }
    return checkr;
  
}


//--------------------------------------------------------------------------
Double_t AliT0QAChecker::CheckESD(TObjArray *listrec ) const
{
  Float_t checkr = 0;
  TH1 *fhESD;
 
  fhESD  = (TH1*) listrec->UncheckedAt(2);
  if(fhESD){
    AliDebug(AliQAv1::GetQADebugLevel(), Form("count %s ", fhESD->GetName()) );
     TF1 *f1 = new TF1("f1","gaus",-1,1);
    fhESD->Fit("f1","R","Q", -1,1);
    Double_t par[3];
    f1->GetParameters(&par[0]);
    
    TPaveText text(0.30,0.50,0.99,0.99,"NDC");
    
    text.AddText(Form("T0 RUN %d ",AliCDBManager::Instance()->GetRun()));
    
    AliDebug(AliQAv1::GetQADebugLevel(), Form("numentries %d mean %f  #sigma %f", (int)fhESD->GetEntries(),par[1], par[2]));
    
    
    if (par[2] > 0.07 && par[2] < 1.) {
      checkr=0.5;
       text.AddText(Form("not good resolution :\n %f ns\n", par[2] ));
       text.SetFillColor(5);
       printf("T0 detector resolution is not good enouph: %f ns\n",par[2] );
    }
    if(TMath::Abs(par[1])>0.05) {
      checkr = 0.5;
      text.AddText(Form(" Check clock shift on %f ns", par[1]));
      text.SetFillColor(5);
    }
    if (par[2] >  1. || TMath::Abs(par[1])>0.1) {
      checkr = 0.25;
      text.AddText(Form(" Bad resolution:\n mean %f ns sigma %f ns", par[1], par[2]));
      text.SetFillColor(2);
      { // RS Clean previous additions
	TList* lstF = fhESD->GetListOfFunctions();
	if (lstF) {
	  TObject *stats = lstF->FindObject("stats");
	  lstF->Remove(stats);
	  TObject *obj;
	  while ((obj = lstF->First())) {
	    while(lstF->Remove(obj)) { }
	    delete obj;
	  }
	  if (stats) lstF->Add(stats);
	} 
      }
      fhESD->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
      AliDebug(AliQAv1::GetQADebugLevel(),
	       Form("Please, check calibration: shift= %f resolution %f test=%f\n",
		    par[1], par[2], checkr) ) ; 
    }
  }
  else
    {
      AliDebug(AliQAv1::GetQADebugLevel(),
	       Form("No ESD QA histogram found, nothing to check"));
	checkr=0;
    }
  
  
  return checkr;
}
