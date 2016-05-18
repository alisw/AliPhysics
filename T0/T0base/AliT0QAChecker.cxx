/**************************************************************************
 * Coyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
#include <TLegend.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"
#include "AliT0QAChecker.h"
#include "AliQAThresholds.h"
#include "AliDAQ.h"

ClassImp(AliT0QAChecker)
//____________________________________________________________________________
AliT0QAChecker::AliT0QAChecker() :
  AliQACheckerBase("T0","T0 Quality Assurance Checker"),
  fCFDErrorThreshold(0),
  fLEDErrorThreshold(0),
  fQTCErrorThreshold(0),
  fRatioCFDEffLEDEffErrorThreshold(1),
  fQTCEfficiencyErrorThreshold(0),
  fBCIDPeriodParam(3564),
  fBCIDOffsetParam(37),
  fBCIDBandWidthParam(10),
  fTZeroAPlusCErrorThreshold(2000.0),
  fTZeroAMinusCErrorThreshold(2000.0)
{
  // Standard constructor
  for(Int_t i=0; i<24; i++){ 
    fMeanCFDFromGoodRunParam[i]=0; 
    fMeanLEDFromGoodRunParam[i]=0; 
    fMeanQTCFromGoodRunParam[i]=0; 
  }

}

//____________________________________________________________________________
AliT0QAChecker::AliT0QAChecker(const AliT0QAChecker& qac):
  AliQACheckerBase(qac.GetName(), qac.GetTitle()),
  fCFDErrorThreshold(qac.fCFDErrorThreshold),
  fLEDErrorThreshold(qac.fLEDErrorThreshold),
  fQTCErrorThreshold(qac.fQTCErrorThreshold),
  fRatioCFDEffLEDEffErrorThreshold(qac.fRatioCFDEffLEDEffErrorThreshold),
  fQTCEfficiencyErrorThreshold(qac.fQTCEfficiencyErrorThreshold),
  fBCIDPeriodParam(qac.fBCIDPeriodParam),
  fBCIDOffsetParam(qac.fBCIDOffsetParam),
  fBCIDBandWidthParam(qac.fBCIDBandWidthParam),
  fTZeroAPlusCErrorThreshold(qac.fTZeroAPlusCErrorThreshold),
  fTZeroAMinusCErrorThreshold(qac.fTZeroAMinusCErrorThreshold)
{
  // copy constructor
  AliError("Copy should not be used with this class\n");
  for(Int_t i=0; i<24; i++){ 
    fMeanCFDFromGoodRunParam[i]=qac.fMeanCFDFromGoodRunParam[i]; 
    fMeanLEDFromGoodRunParam[i]=qac.fMeanLEDFromGoodRunParam[i]; 
    fMeanQTCFromGoodRunParam[i]=qac.fMeanQTCFromGoodRunParam[i]; 
  }

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

  AliCDBManager* man = AliCDBManager::Instance();
  //man->SetDefaultStorage(gSystem->Getenv("AMORE_CDB_URI"));
  if(!man) return; 
  AliCDBEntry* entry = man->Get("GRP/Calib/QAThresholds");
  if(!entry) return; 
  TObjArray* t0branch = (TObjArray*) entry->GetObject();
  if(!list) return;
  AliQAThresholds*  thresholds = (AliQAThresholds*) t0branch->FindObject("T00");
  // here you should test that you got a non-null pointer

  if(!thresholds) return;
  if(AliDAQ::DetectorID("T0")!= thresholds->GetDetectorId()){
    AliInfo(Form("DETECTOR ID %d DOES NOT MATCH TO TZERO",thresholds->GetDetectorId()));
    return;
  }

  int iparam; 
  for(int ipmt=0; ipmt<24;ipmt++){ 
    iparam = ipmt + 1; //current consecutive number of parameter
    if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean CFD from a good run 
      fMeanCFDFromGoodRunParam[ipmt] = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
    }

    iparam = ipmt + 25;
    if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean LED from a good run 
      fMeanLEDFromGoodRunParam[ipmt] = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
    } 
    iparam = ipmt + 49;
    if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean QTC from a good run 
      fMeanQTCFromGoodRunParam[ipmt] = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
    } 
  }
  iparam = 73; //CFD threshold
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ 
    fCFDErrorThreshold = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  }
  iparam = 74; //LED threshold
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ 
    fLEDErrorThreshold = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  }
  iparam = 75; //QTC threshold
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ 
    fQTCErrorThreshold = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  }
 
  iparam = 82; //Error level threshold on CFD efficiency/LED efficiency ratio
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ 
    fRatioCFDEffLEDEffErrorThreshold = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  }
  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!

  iparam = 83; //Error level threshold on QTC efficiency 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ 
    fQTCEfficiencyErrorThreshold = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  }

  iparam = 84; 
  if((TParameter<int>*) thresholds->GetThreshold(iparam)){ 
    fBCIDPeriodParam = ((TParameter<int>*) thresholds->GetThreshold(iparam))->GetVal();
  }

  iparam = 85; 
  if((TParameter<int>*) thresholds->GetThreshold(iparam)){ 
    fBCIDOffsetParam = ((TParameter<int>*) thresholds->GetThreshold(iparam))->GetVal();
  }

  iparam = 86; 
  if((TParameter<int>*) thresholds->GetThreshold(iparam)){ 
    fBCIDBandWidthParam = ((TParameter<int>*) thresholds->GetThreshold(iparam))->GetVal();
  }

  iparam = 87; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ 
    fTZeroAPlusCErrorThreshold = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  }

  iparam = 88; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ 
    fTZeroAMinusCErrorThreshold = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  }


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

      //if(AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib){//      if (index == AliQAv1::kRAW )
        //check laser data efficiencies   
      //  Double_t qaFlag = CheckLaser(list[specie]);
      //  if(qaFlag < test[specie]) test[specie] = qaFlag;
      //}

      //if(//AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib   ||
       //  AliRecoParam::ConvertIndex(specie) == AliRecoParam::kHighMult ||
      //   AliRecoParam::ConvertIndex(specie) == AliRecoParam::kLowMult){ 
         //AliRecoParam::ConvertIndex(specie) == AliRecoParam::kDefault ||
         
        //check BCID   
      //  Double_t qaFlag = CheckBCID(list[specie]);
      //  if(qaFlag < test[specie]) test[specie] = qaFlag;
      //}

      if(AliRecoParam::ConvertIndex(specie) == AliRecoParam::kHighMult ||
        AliRecoParam::ConvertIndex(specie) == AliRecoParam::kLowMult){ 
        //AliRecoParam::ConvertIndex(specie) == AliRecoParam::kDefault ||
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
//Double_t AliT0QAChecker::CheckLaser(TObjArray *listrec) const {
//   
//  return 1.0; 
//}

//--------------------------------------------------------------------------
Double_t AliT0QAChecker::CheckRaw(TObjArray *listrec) const {
   
  //Fk Set drawing options for LED and CFD efficiencies from the raw data
  TH1F *hCFDEffData = (TH1F*) listrec->UncheckedAt(207);//hRawTrigger 
  TH1F *hLEDEffData = (TH1F*) listrec->UncheckedAt(208);//hRawTrigger 

  //clean objects added at previous checks
  EraseOldMessages((TH1*) hCFDEffData); 
  hCFDEffData->GetListOfFunctions()->Add((TH1D*) hLEDEffData->Clone());	      

  TLegend leg(0.12,0.76,0.9,0.92," ","brNDC");
  leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.04); leg.SetNColumns(2);
  leg.AddEntry((TH1D*) hCFDEffData,"CFD","p");
  leg.AddEntry((TH1D*) hLEDEffData,"LED","p");
  hCFDEffData->GetListOfFunctions()->Add((TLegend*) leg.Clone());	      


  //Fk Draw CFD-mean for each PMT
  TH2F* fhCFD = (TH2F*) listrec->UncheckedAt(210);
  TH1F* fhCFDSubtrMean = (TH1F*) listrec->UncheckedAt(231);

  EraseOldMessages((TH1*) fhCFDSubtrMean); 
  for(int ipmt=0; ipmt<24;ipmt++){ 
    TH1F*  hProjDummy = (TH1F*) fhCFD->ProjectionY("dummy",ipmt+1,ipmt+1);
    Float_t mean=0.0, rms=0.0;
    GetMeanAndRmsAroundMainMaximum(mean, rms,  hProjDummy,0);

    Float_t deviation = mean - fMeanCFDFromGoodRunParam[ipmt]; 

    fhCFDSubtrMean->SetBinContent(ipmt+1,deviation);
    fhCFDSubtrMean->SetBinError(ipmt+1,rms);
      
    delete hProjDummy;
  }
  TLine linelowredCFD(0, fCFDErrorThreshold, 24, fCFDErrorThreshold);    
  linelowredCFD.SetLineColor(2);
  linelowredCFD.SetLineStyle(3);
  linelowredCFD.SetLineWidth(4);
  TLine linehighredCFD(0, -fCFDErrorThreshold, 24, -fCFDErrorThreshold);    
  linehighredCFD.SetLineColor(2);
  linehighredCFD.SetLineStyle(3);
  linehighredCFD.SetLineWidth(4);
  fhCFDSubtrMean->GetListOfFunctions()->Add((TLine*) linelowredCFD.Clone());	      
  fhCFDSubtrMean->GetListOfFunctions()->Add((TLine*) linehighredCFD.Clone());	      

 
 
  //Fk Draw LED-mean for each PMT
  TH2F* fhLED = (TH2F*) listrec->UncheckedAt(211);
  TH1F* fhLEDSubtrMean = (TH1F*) listrec->UncheckedAt(232);
 
  EraseOldMessages((TH1*) fhLEDSubtrMean); 
  for(int ipmt=0; ipmt<24;ipmt++){ 
    TH1F*  hProjDummy = (TH1F*) fhLED->ProjectionY("dummy",ipmt+1,ipmt+1);
    Float_t mean=0.0, rms=0.0;
    GetMeanAndRmsAroundMainMaximum(mean, rms,  hProjDummy,1);
    Float_t deviation = mean - fMeanLEDFromGoodRunParam[ipmt]; 

    fhLEDSubtrMean->SetBinContent(ipmt+1,deviation);
    fhLEDSubtrMean->SetBinError(ipmt+1,rms);
      
    delete hProjDummy;
  }
  TLine linelowredLED(0, fLEDErrorThreshold, 24, fLEDErrorThreshold);    
  linelowredLED.SetLineColor(2);
  linelowredLED.SetLineStyle(3);
  linelowredLED.SetLineWidth(4);
  TLine linehighredLED(0, -fLEDErrorThreshold, 24, -fLEDErrorThreshold);    
  linehighredLED.SetLineColor(2);
  linehighredLED.SetLineStyle(3);
  linehighredLED.SetLineWidth(4);
  fhLEDSubtrMean->GetListOfFunctions()->Add((TLine*) linelowredLED.Clone());	      
  fhLEDSubtrMean->GetListOfFunctions()->Add((TLine*) linehighredLED.Clone());	      

     
  //Fk Draw QTC-mean for each PMT
  TH2F* fhQTC = (TH2F*) listrec->UncheckedAt(212);
  TH1F* fhQTCSubtrMean = (TH1F*) listrec->UncheckedAt(233);
   
  EraseOldMessages((TH1*) fhQTCSubtrMean); 
  for(int ipmt=0; ipmt<24;ipmt++){ 
    TH1F*  hProjDummy = (TH1F*) fhQTC->ProjectionY("dummy",ipmt+1,ipmt+1);
    Float_t mean=0.0, rms=0.0;
    GetMeanAndRmsAroundMainMaximum(mean, rms,  hProjDummy,2);
    Float_t deviation = mean - fMeanQTCFromGoodRunParam[ipmt]; 

    fhQTCSubtrMean->SetBinContent(ipmt+1,deviation);
    fhQTCSubtrMean->SetBinError(ipmt+1,rms);
      
    delete hProjDummy;
  }
  TLine linelowredQTC(0, fQTCErrorThreshold, 24, fQTCErrorThreshold);    
  linelowredQTC.SetLineColor(2);
  linelowredQTC.SetLineStyle(3);
  linelowredQTC.SetLineWidth(4);
  TLine linehighredQTC(0, -fQTCErrorThreshold, 24, -fQTCErrorThreshold);    
  linehighredQTC.SetLineColor(2);
  linehighredQTC.SetLineStyle(3);
  linehighredQTC.SetLineWidth(4);
  fhQTCSubtrMean->GetListOfFunctions()->Add((TLine*) linelowredQTC.Clone());	      
  fhQTCSubtrMean->GetListOfFunctions()->Add((TLine*) linehighredQTC.Clone());	      

  //CFD and LED efficiency in range ~2000- ~3000 
  TH1F* hCFDeffSubRange = (TH1F*) listrec->UncheckedAt(237);
  TH1F* hEffLEDSubRange = (TH1F*) listrec->UncheckedAt(238);
  // ratio CDF eff /LEF eff in subragne 
  TH1F* hRatioCFDLEDeff = (TH1F*) listrec->UncheckedAt(239);//FK   
  EraseOldMessages((TH1*) hRatioCFDLEDeff); 
  int npmt = hRatioCFDLEDeff->GetNbinsX();
  for(int ipmt=1;ipmt<=npmt;ipmt++){
    Float_t c0 = hCFDeffSubRange->GetBinContent(ipmt); 
    Float_t c1 = hEffLEDSubRange->GetBinContent(ipmt);
    if(c1){
      hRatioCFDLEDeff->SetBinContent(ipmt,c0/c1);  
    }else{
      hRatioCFDLEDeff->SetBinContent(ipmt,0);  
    }  
  }

  TLine linelowredRatioCFDLEDeff(0, 1+fRatioCFDEffLEDEffErrorThreshold, 24, 1+fRatioCFDEffLEDEffErrorThreshold);    
  linelowredRatioCFDLEDeff.SetLineColor(2);
  linelowredRatioCFDLEDeff.SetLineStyle(3);
  linelowredRatioCFDLEDeff.SetLineWidth(4);
  TLine linehighredRatioCFDLEDeff(0, 1-fRatioCFDEffLEDEffErrorThreshold, 24, 1-fRatioCFDEffLEDEffErrorThreshold);    
  linehighredRatioCFDLEDeff.SetLineColor(2);
  linehighredRatioCFDLEDeff.SetLineStyle(3);
  linehighredRatioCFDLEDeff.SetLineWidth(4);
  hRatioCFDLEDeff->GetListOfFunctions()->Add((TLine*) linelowredRatioCFDLEDeff.Clone());	      
  hRatioCFDLEDeff->GetListOfFunctions()->Add((TLine*) linehighredRatioCFDLEDeff.Clone());	      

  //        PERFROM CHECKS on HISTOGRAMS

  //-------- triggers -----------
  Int_t qualityFlagTrigger = kT0Info; //init quality flag for a given histogram; 

  TH1F *hTrigger = (TH1F*) listrec->UncheckedAt(169);//hRawTrigger 

  // clean objects added at previous checks
  EraseOldMessages((TH1*) hTrigger); 

  if(hTrigger->Integral()>0){
    //trigger plot does have some counts in it
    //are Mean, ORA and ORC not empty?  
    if(/* hTrigger->GetBinContent(1)<0.001 ||*/ hTrigger->GetBinContent(3)<0.001 || hTrigger->GetBinContent(4)<0.001){
      qualityFlagTrigger = kT0Error; //no entries on diagonal
      AliDebug(AliQAv1::GetQADebugLevel(), Form("T0: too little ORA and ORC in  %s", hTrigger->GetName() ));

      TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
      text.AddText(Form("Check ORA and ORC")); 
      text.AddText(Form("Report problem to the T0 on-call expert")); 
      text.SetBorderSize(0);
      text.SetFillStyle(0);
      hTrigger->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
    }

  }else{ //Trigger histo empty

    qualityFlagTrigger = kT0Error;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 histogram  %s has NO entries", hTrigger->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("NO ENTRIES!!!")); 
    text.AddText(Form("If T0 is READY report")); 
    text.AddText(Form("readout problem to the T0 on-call expert")); 
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hTrigger->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      

  }
  //---------- CFD eff/LED eff within subrange  -----------
  Int_t qualityFlagRatioCFDeffLEDeff = kT0Info; //init quality flag for a given histogram;
  int nPMTs = hRatioCFDLEDeff->GetNbinsX(); 

  for(int ipmt=1; ipmt<=nPMTs; ipmt++){
    if(TMath::Abs( hRatioCFDLEDeff->GetBinContent(ipmt) -1) > fRatioCFDEffLEDEffErrorThreshold){ //mean is expected to be around 1
      qualityFlagRatioCFDeffLEDeff = kT0Error;
    }
  }
  if(qualityFlagRatioCFDeffLEDeff == kT0Error){
    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("Problem with efficiency ratio CFD/LED !!!")); 
    text.AddText(Form("If T0 is READY and beam is on report")); 
    text.AddText(Form("the problem to the T0 on-call expert"));
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hRatioCFDLEDeff->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }
  

  //---------- CFD -  mean CFD  -----------
  Int_t qualityFlagCFDSubtr = kT0Info; //init quality flag for a given histogram;
  nPMTs = fhCFDSubtrMean->GetNbinsX(); 

  for(int ipmt=1; ipmt<=nPMTs; ipmt++){
    if(TMath::Abs( fhCFDSubtrMean->GetBinContent(ipmt)) > fCFDErrorThreshold){
      qualityFlagCFDSubtr = kT0Error;

   }
  }
  if(qualityFlagCFDSubtr == kT0Error){
    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("Problem with CFD timing  !!!")); 
    text.AddText(Form("If T0 is READY and beam is on report")); 
    text.AddText(Form("the problem to the T0 on-call expert"));
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    fhCFDSubtrMean->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }
  //--------- QTC efficiency ---------------

  Int_t qualityFlagEffQTC = kT0Info; //init quality flag for a given histogram;
  TH1F *hEffQTC = (TH1F*) listrec->UncheckedAt(209);//hRawTrigger 
  
  EraseOldMessages((TH1*) hEffQTC); // clean objects added at previous checks

  nPMTs = hEffQTC->GetNbinsX(); 
  for(int ipmt=1; ipmt<=nPMTs; ipmt++){
    if(TMath::Abs( hEffQTC->GetBinContent(ipmt)) < fQTCEfficiencyErrorThreshold){
      qualityFlagEffQTC = kT0Error;
    }
  }
  if( qualityFlagEffQTC == kT0Error){
    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("Problem with QTC efficiency !!!")); 
    text.AddText(Form("If T0 is READY and beam is on report")); 
    text.AddText(Form("the problem to the T0 on-call expert"));
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hEffQTC->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }
  TLine linehighredQTCeff(0, fQTCEfficiencyErrorThreshold, 24, fQTCEfficiencyErrorThreshold);    
  linehighredQTCeff.SetLineColor(2);
  linehighredQTCeff.SetLineStyle(3);
  linehighredQTCeff.SetLineWidth(4);
  hEffQTC->GetListOfFunctions()->Add((TLine*) linehighredQTCeff.Clone());	      
 
  //---------- BCID --------------------
  Int_t qualityFlagBCID = kT0Info; //init quality flag for a given histogram; 
  TH2F *hBCID = (TH2F*) listrec->UncheckedAt(236); //BCID versus TRM  BCID
     
  // clean objects added at previous checks
  EraseOldMessages((TH1*)hBCID);

  if(hBCID->Integral()>0){
    //BCID does have some counts in it

    Int_t nbinsX = hBCID->GetNbinsX();
    Int_t startX = hBCID->GetXaxis()->FindBin(100); //skip region close to orbit period
    Int_t nbinsY = hBCID->GetNbinsY();
    Int_t binWidthX = (Int_t) hBCID->GetXaxis()->GetBinWidth(1);
    Int_t binWidthY = (Int_t) hBCID->GetYaxis()->GetBinWidth(1);
    double entriesOnDiagonal  = 0; //count diagonal and off diagonal entries
    double entriesOffDiagonal = 0;

    for(Int_t itrm=startX; itrm<=nbinsX; itrm++){ //BCID TRM
      for(Int_t ibcid=1; ibcid<=nbinsY; ibcid++){ //BCID
        if(TMath::Abs( (itrm*binWidthX - fBCIDOffsetParam) % fBCIDPeriodParam - binWidthY*ibcid) < fBCIDBandWidthParam ){ 
          entriesOnDiagonal  += hBCID->GetBinContent(itrm,ibcid); //On  Diagonal
          //hBCID->Fill(itrm*binWidthX,ibcid*binWidthY,0.001); // visualize the diagonal belt
        }else{
          entriesOffDiagonal += hBCID->GetBinContent(itrm,ibcid); //Off Diagonal
        }
      }
    }
    if(entriesOnDiagonal<1 || entriesOffDiagonal>20){
      qualityFlagBCID = kT0Error; //no entries on diagonal
      AliDebug(AliQAv1::GetQADebugLevel(), Form("T0   %s is not diagonal", hBCID->GetName() ));

      TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
      text.AddText(Form("Check if entries are on diagonal")); 
      text.AddText(Form("If NOT - report readout problem to the T0 on-call expert")); 
      text.SetBorderSize(0);
      text.SetFillStyle(0);
      hBCID->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
    }
  }else{ //BCID empty

    qualityFlagBCID = kT0Error;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 :  %s has NO entries", hBCID->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("NO ENTRIES!!!")); 
    text.AddText(Form("If T0 is READY report")); 
    text.AddText(Form("readout problem to the T0 on-call expert"));
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hBCID->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }

  //--------- mean versus vertex 1st ---------
  Int_t qualityFlagMeanVersusVertex1st = kT0Info; //init quality flag for a given histogram;
  TH2F *hMeanVersusVertex1st  = (TH2F*) listrec->UncheckedAt(220);  //fhMeanBest  time
  EraseOldMessages((TH1*) hMeanVersusVertex1st); 
 
  if(hMeanVersusVertex1st->Integral()<1){
    qualityFlagMeanVersusVertex1st = kT0Error;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 histogram  %s has NO entries", hMeanVersusVertex1st->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("NO ENTRIES!!!")); 
    text.AddText(Form("If T0 is READY and beam is on report")); 
    text.AddText(Form("the problem to the T0 on-call expert"));
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hMeanVersusVertex1st->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }


  //--------- vertex TVDC on ---------
  Int_t qualityFlagVertex1stTVDCon = kT0Info; //init quality flag for a given histogram;
  TH1F *hVertex1stTVDCon  = (TH1F*) listrec->UncheckedAt(223);  //fhMeanBest  time
  EraseOldMessages((TH1*) hVertex1stTVDCon); 
 
  if(hVertex1stTVDCon->Integral()<1){
    qualityFlagVertex1stTVDCon = kT0Error;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 histogram  %s has NO entries", hVertex1stTVDCon->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("NO ENTRIES!!!")); 
    text.AddText(Form("If T0 is READY and beam is on report")); 
    text.AddText(Form("the problem to the T0 on-call expert"));
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hVertex1stTVDCon->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }else{
    if(TMath::Abs(hVertex1stTVDCon->GetMean()) > fTZeroAMinusCErrorThreshold){ 
      qualityFlagVertex1stTVDCon = kT0Error;

      TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
      text.AddText(Form("Displaced vertex")); 
      text.AddText(Form("If T0 is READY and beam is on report")); 
      text.AddText(Form("the problem to the T0 on-call expert"));
      text.SetBorderSize(0);
      text.SetFillStyle(0);
      hVertex1stTVDCon->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
    }
  }

  //--------- vertex TVDC off ---------
  Int_t qualityFlagVertex1stTVDCoff = kT0Info; //init quality flag for a given histogram;
  TH1F *hVertex1stTVDCoff  = (TH1F*) listrec->UncheckedAt(225);  //fhMeanBest  time
  EraseOldMessages((TH1*) hVertex1stTVDCoff); 
 
  if(hVertex1stTVDCoff->Integral()<1){
    qualityFlagVertex1stTVDCoff = kT0Warning;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 histogram  %s has NO entries", hVertex1stTVDCoff->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("Warning: NO ENTRIES")); 
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hVertex1stTVDCoff->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }



  //----------- time TVDC on ---------
    
  Int_t qualityFlagMean1stTVDCon = kT0Info; //init quality flag for a given histogram;
  TH1F *hMean1stTVDCon  = (TH1F*) listrec->UncheckedAt(226);  //fhMeanBest  time
  EraseOldMessages((TH1*) hMean1stTVDCon); 
 
  if(hMean1stTVDCon->Integral()<1){
    qualityFlagMean1stTVDCon = kT0Error;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 histogram  %s has NO entries", hMean1stTVDCon->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("NO ENTRIES!!!")); 
    text.AddText(Form("If T0 is READY and beam is on report")); 
    text.AddText(Form("the problem to the T0 on-call expert"));
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hMean1stTVDCon->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }else{
    //cout<<"Mean: "<<TMath::Abs(hMean1stTVDCon->GetMean())<<" threshold "<<fTZeroAPlusCErrorThreshold<<endl;
    if(TMath::Abs(hMean1stTVDCon->GetMean()) > fTZeroAPlusCErrorThreshold){ 
      qualityFlagMean1stTVDCon = kT0Error;

      TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
      text.AddText(Form("Shift of mean time")); 
      text.AddText(Form("If T0 is READY and beam is on report")); 
      text.AddText(Form("the problem to the T0 on-call expert"));
      text.SetBorderSize(0);
      text.SetFillStyle(0);
      hMean1stTVDCon->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
    }
  }

  //----------- time TVDC off ---------
    
  Int_t qualityFlagMean1stTVDCoff = kT0Info; //init quality flag for a given histogram;
  TH1F *hMean1stTVDCoff  = (TH1F*) listrec->UncheckedAt(227);  //fhMeanBest  time
  EraseOldMessages((TH1*) hMean1stTVDCoff); 
 
  if(hMean1stTVDCoff->Integral()<1){
    qualityFlagMean1stTVDCoff = kT0Warning;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("T0 histogram  %s has NO entries", hMean1stTVDCoff->GetName() ));

    TPaveText text(0.20,0.50,0.99,0.99,"NDC");   
    text.AddText(Form("Warning: NO ENTRIES")); 
    text.SetBorderSize(0);
    text.SetFillStyle(0);
    hMean1stTVDCoff->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
  }
 


  //----------------------- executive summary ---------------------
  int lowestQualityFlag = (int) qualityFlagTrigger;
  if(qualityFlagRatioCFDeffLEDeff   < lowestQualityFlag)  lowestQualityFlag = qualityFlagRatioCFDeffLEDeff;
  if(qualityFlagCFDSubtr            < lowestQualityFlag)  lowestQualityFlag = qualityFlagCFDSubtr;
  if(qualityFlagEffQTC              < lowestQualityFlag)  lowestQualityFlag = qualityFlagEffQTC;
  if(qualityFlagMeanVersusVertex1st < lowestQualityFlag)  lowestQualityFlag = qualityFlagMeanVersusVertex1st;
  if(qualityFlagBCID                < lowestQualityFlag)  lowestQualityFlag = qualityFlagBCID;
  if(qualityFlagVertex1stTVDCon     < lowestQualityFlag)  lowestQualityFlag = qualityFlagVertex1stTVDCon;
  if(qualityFlagVertex1stTVDCoff    < lowestQualityFlag)  lowestQualityFlag = qualityFlagVertex1stTVDCoff;
  if(qualityFlagMean1stTVDCon       < lowestQualityFlag)  lowestQualityFlag = qualityFlagMean1stTVDCon;
  if(qualityFlagMean1stTVDCoff      < lowestQualityFlag)  lowestQualityFlag = qualityFlagMean1stTVDCoff;
  

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
 
//--------------------------------------------------------------------------
Float_t AliT0QAChecker::GetMeanAboveThreshold(TH1F* hV, Float_t thr) const{
  //caculate mean value of histo bins above threshold
  Int_t nBins = hV->GetNbinsX();
  Int_t nBinsAboveThr = 0; 
  Float_t sum = 0;

  for(Int_t ib=1;ib<=nBins;ib++){
    Float_t val = hV->GetBinContent(ib);
    if(val<thr){
      sum+=val;
      nBinsAboveThr++; 
    }
  }

  if(nBinsAboveThr>0) return sum/nBinsAboveThr;
  else return hV->GetMean();
}
    
//--------------------------------------------------------------------------
void AliT0QAChecker::GetMeanAndRmsAroundMainMaximum(Float_t &meanHisto,Float_t &rmsHisto, TH1F *histo, int type) const{

  if(!histo){
    meanHisto=0.0;
    rmsHisto=0.0;
    return;
  }
  if(histo->Integral()<0.00001){
    meanHisto=0.0;
    rmsHisto=0.0;
    return;
  }

  double nSigma      = 3.0; //n sigma window around mean and main maximum 
  double expectedRMS = 13.0; // expected rms of the main peak in case of CFD
  if(type == 1) expectedRMS = 34.0; //LED
  if(type == 2) expectedRMS = 34.0; //QTC

  //0) approx of mean is global maximum
  Int_t   nb     =  histo->GetNbinsX();
  Int_t   bmax   =  histo->GetMaximumBin();
  Float_t xmax   =  histo->GetBinCenter(bmax);

  double window = expectedRMS * nSigma;
  //1) estimate a mean in the nSigma window around main max
  int nlow = histo->FindBin( xmax - window);
  int nhigh = histo->FindBin( xmax + window);
  if(nlow<1)   nlow = 1; 
  if(nhigh>nb) nhigh = nb;

  Float_t sum=0.0, sumWeight=0.0, sumWeightSqr=0.0;
  for(int ii=nlow;ii<=nhigh;ii++){
    sum          += histo->GetBinContent(ii);  
    sumWeight    += histo->GetBinContent(ii)*histo->GetBinCenter(ii);
  }  
  if(sum>0.0){
    meanHisto = sumWeight/sum; //mean in 1st itteration
  }else{
    meanHisto = 0.0;
    rmsHisto=0.0;
    return;
  }
  
  //2) recalculte mean and rms in the nSigma window around mean 1)
  nlow  = histo->FindBin( meanHisto - window);
  nhigh = histo->FindBin( meanHisto + window);
  if(nlow<1)   nlow = 1; 
  if(nhigh>nb) nhigh = nb;

  sum=0.0; sumWeight=0.0; sumWeightSqr=0.0;
  for(int ii=nlow;ii<=nhigh;ii++){
    sum          += histo->GetBinContent(ii);  
    sumWeight    += histo->GetBinContent(ii)*histo->GetBinCenter(ii);
    sumWeightSqr += histo->GetBinContent(ii)*histo->GetBinCenter(ii)*histo->GetBinCenter(ii);
  }
  if(sum>0.0){
    meanHisto = sumWeight/sum; //mean in 2nd itteration
    rmsHisto  = sumWeightSqr/sum - meanHisto*meanHisto;//rms square
    if(rmsHisto > 0.0) rmsHisto = sqrt(rmsHisto); //rms
    else  rmsHisto = 0.0;
    return;
  }else{
    meanHisto = 0.0;
    rmsHisto  = 0.0;
    return;
  }
}
