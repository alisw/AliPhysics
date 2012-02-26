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

    
  for(Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++){ 
    test[specie]    = 1.0; //FK//  initiate qa flag for the whole set of histograms as good 
  }


  for(Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if(!(AliQAv1::Instance()->IsEventSpecieSet(specie) && list[specie]) || list[specie]->GetEntries() == 0) {
      continue;
    }
    if(index == AliQAv1::kRAW){

      if(AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib){//      if (index == AliQAv1::kRAW )
        //check laser data efficiencies   
        Double_t qaFlag = CheckLaser(list[specie]);
        if(qaFlag < test[specie]) test[specie] = qaFlag;
      }

      if(AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib   ||
         AliRecoParam::ConvertIndex(specie) == AliRecoParam::kDefault){ 
         
        //check BCID   
        Double_t qaFlag = CheckBCID(list[specie]);
        if(qaFlag < test[specie]) test[specie] = qaFlag;
      }

      if(AliRecoParam::ConvertIndex(specie) == AliRecoParam::kDefault){ 
        //check physics 
        Double_t qaFlag = CheckRaw(list[specie]);
        if(qaFlag < test[specie]) test[specie] = qaFlag;
      }
    }

    if(index == AliQAv1::kESD && AliRecoParam::Convert(specie) != AliRecoParam::kCalib){
      test[specie] = CheckESD(list[specie]);
    } 
  }
}

//--------------------------------------------------------------------------
Double_t AliT0QAChecker::CheckLaser(TObjArray *listrec) const {
   
  TH1 *hdata; 
  TH1 *fhRawEff[10];
  Int_t nEffHistos=0;

  //thresholds for warning and error on efficiencies
  Float_t thrWarning = 0.5; //FK//  warning level
  Float_t thrError   = 0.2; //FK//  error level

  const int kNumberOfHistos = 3; 
  Int_t consecutiveHistoNumber[kNumberOfHistos] = { 207, 208, 209};  //Checked histos   fhCDFeff, hEffLED, hEffQTC
  Int_t qualityFlag[kNumberOfHistos]; //quality flag for a given histogram

  for(Int_t ir=0; ir<kNumberOfHistos; ir++){
    qualityFlag[ir] = kT0Info; //init quality flag for a given histogram

    hdata = (TH1*) listrec->UncheckedAt(consecutiveHistoNumber[ir]);
    if(hdata){
      fhRawEff[nEffHistos] = hdata;
      nEffHistos++;
    }
  }
     
  TLine linelowyellow(0, thrWarning, 24, thrWarning);    
  linelowyellow.SetLineColor(5);
  linelowyellow.SetLineStyle(3);
  linelowyellow.SetLineWidth(4);
  TLine linelowred(0, thrError, 24, thrError);    
  linelowred.SetLineColor(2);
  linelowred.SetLineStyle(3);
  linelowred.SetLineWidth(4);

  Bool_t bEffHistosNotEmpty = kFALSE; //check if all histograms have some counts
 
  for(Int_t ih = 0; ih < nEffHistos; ih++){
     
    EraseOldMessages((TH1*) fhRawEff[ih]);// clean objects added at previous checks
 
    fhRawEff[ih]->SetLineWidth(2);
    fhRawEff[ih]->SetMaximum(2.);
    fhRawEff[ih]->SetMinimum(0.);
    fhRawEff[ih]->GetListOfFunctions()->Add((TLine*)linelowyellow.Clone());
    fhRawEff[ih]->GetListOfFunctions()->Add((TLine*)linelowred.Clone());

    if(fhRawEff[ih]->Integral()>0) bEffHistosNotEmpty = kTRUE; //this histo does have some counts in it

    Int_t nbins= fhRawEff[ih]->GetNbinsX();
    for(Int_t ib=1; ib<=nbins; ib++){ //loop over bins and check if the efficiency is above level
        
      Float_t chcont = fhRawEff[ih]->GetBinContent(ib);
      if(chcont < thrWarning && qualityFlag[ih] > kT0Error  ) qualityFlag[ih] = kT0Warning;//Warning level
      if(chcont < thrError)  qualityFlag[ih] = kT0Error;//Error level
    }
   
    if(qualityFlag[ih] == kT0Info ){
      AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 efficiency  %s  is good", fhRawEff[ih]->GetName() ));
    }else if(qualityFlag[ih] == kT0Warning){ 
      AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 efficiency  %s  is not so good", fhRawEff[ih]->GetName() ));
    }else if(qualityFlag[ih] == kT0Error){
      AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 efficiency  %s  is not good", fhRawEff[ih]->GetName() ));
    }
  }

  //executive summary
  int lowestQualityFlag = (int) kT0Info;
  for(Int_t ih = 0; ih < nEffHistos; ih++){

    if(!bEffHistosNotEmpty){ //all laser efficiency plots are empty
      TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
      text.AddText(Form("1) T0 is in BEAMTUNIG: empty plots are ok")); 
      text.AddText(Form("2) T0 is in READY: check calibriation trigger")); 
      text.AddText(Form("if also physics data are empty report"));
      text.AddText(Form("readout problem to the T0 on-call expert")); 
      fhRawEff[ih]->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
    }
 
    if( qualityFlag[ih] <lowestQualityFlag )  lowestQualityFlag = qualityFlag[ih];
  }
   
  return ConvertQualityFlagToDouble(lowestQualityFlag); 
}
//--------------------------------------------------------------------------
Double_t AliT0QAChecker::CheckBCID(TObjArray *listrec) const {
   
  Int_t qualityFlagBCID = kT0Info; //init quality flag for a given histogram; 

  TH2F *hBCID = (TH2F*) listrec->UncheckedAt(224); //BCID versus TRM  BCID

     
  // clean objects added at previous checks
  EraseOldMessages((TH1*)hBCID);

  if(hBCID->Integral()>0){
    //BCID does have some counts in it

    Int_t nbinsX = hBCID->GetNbinsX();
    Int_t nbinsY = hBCID->GetNbinsY();
    double entriesOnDiagonal  = 0; //count diagonal and off diagonal entries
    double entriesOffDiagonal = 0;

    for(Int_t ix=1; ix<=nbinsX; ix++){ 
      for(Int_t iy=1; iy<=nbinsY; iy++){ 
        if(TMath::Abs(ix-iy)<6) entriesOnDiagonal  += hBCID->GetBinContent(ix,iy); //On  Diagonal
        else       entriesOffDiagonal += hBCID->GetBinContent(ix,iy); //Off Diagonal
      }
    }
    if(entriesOnDiagonal<1 || entriesOffDiagonal>0){
      qualityFlagBCID = kT0Error; //no entries on diagonal
      AliDebug(AliQAv1::GetQADebugLevel(), Form("T0   %s is not diagonal", hBCID->GetName() ));

      TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
      text.AddText(Form("Check if entries are on a diagonal.")); 
      text.AddText(Form("Report readout problem to the T0 on-call expert")); 
      hBCID->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
    }
  }else{ //BCID empty

    qualityFlagBCID = kT0Error;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 :  %s has NO entries", hBCID->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("NO ENTRIES!!!")); 
    text.AddText(Form("If T0 is READY report")); 
    text.AddText(Form("readout problem to the T0 on-call expert")); 
    hBCID->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }

  //executive summary
  int lowestQualityFlag = (int) qualityFlagBCID;

  return ConvertQualityFlagToDouble(lowestQualityFlag); 
  
}

//--------------------------------------------------------------------------
Double_t AliT0QAChecker::CheckRaw(TObjArray *listrec) const {
   

  Int_t qualityFlagTrigger = kT0Info; //init quality flag for a given histogram; 

  TH1F *hTrigger = (TH1F*) listrec->UncheckedAt(169);//hRawTrigger 

     
  // clean objects added at previous checks
  EraseOldMessages((TH1*) hTrigger); 

  if(hTrigger->Integral()>0){
    //trigger plot does have some counts in it
    //are Mean, ORA and ORC not empty?  
    if( hTrigger->GetBinContent(1)<0.001 || hTrigger->GetBinContent(3)<0.001 || hTrigger->GetBinContent(4)<0.001){
      qualityFlagTrigger = kT0Error; //no entries on diagonal
      AliDebug(AliQAv1::GetQADebugLevel(), Form("T0: too little ORA and ORC in  %s", hTrigger->GetName() ));

      TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
      text.AddText(Form("Check ORA and ORC")); 
      text.AddText(Form("Report problem to the T0 on-call expert")); 
      hTrigger->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
    }
  }else{ //Trigger histo empty

    qualityFlagTrigger = kT0Error;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 histogram  %s has NO entries", hTrigger->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("NO ENTRIES!!!")); 
    text.AddText(Form("If T0 is READY report")); 
    text.AddText(Form("readout problem to the T0 on-call expert")); 
    hTrigger->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }

  //executive summary
  int lowestQualityFlag = (int) qualityFlagTrigger;
   

  return ConvertQualityFlagToDouble(lowestQualityFlag); 
  
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


//--------------------------------------------------------------------------
void AliT0QAChecker::EraseOldMessages(TH1* h) const 
{
  //erase the old captions 
  TList* lstF = h->GetListOfFunctions();
  if(lstF){  
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
//--------------------------------------------------------------------------
Double_t AliT0QAChecker::ConvertQualityFlagToDouble(int qualityFlag) const 
{
  //covert quality flag to double
  Double_t checkr=1.0;

  switch ( qualityFlag ){
    case kT0Info:
        checkr = 1.0; break;
    case kT0Warning:
        checkr = 0.75; break;
    case kT0Error:
          checkr = 0.25; break;
    case kT0Fatal:
        checkr = -1.0; break;
    default:
         AliError("Invalid ecc value. FIXME !");
         checkr = 0.25; break;
  };

  return checkr; 
} 
