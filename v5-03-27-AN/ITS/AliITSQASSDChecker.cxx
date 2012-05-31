/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

// *****************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  P. Cerello Apr 2008
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
//#include "Riostream.h"

// --- AliRoot header files ---
#include "AliITSQASSDChecker.h"
#include "AliITSQADataMakerRec.h"
#include "AliLog.h"

ClassImp(AliITSQASSDChecker)
//__________________________________________________________________
AliITSQASSDChecker& AliITSQASSDChecker::operator = (const AliITSQASSDChecker& qac ) 
{
  // Equal operator.
  this->~AliITSQASSDChecker();
  new(this) AliITSQASSDChecker(qac);
  return *this;
}

void AliITSQASSDChecker::CheckRaws(TH1* histo) {  
  // checker for RAWS
  Double_t minSSDDataSize = 0;
  Double_t maxSSDDataSize = 200;
  Double_t minDDLDataSize = 0;
  Double_t maxDDLDataSize = 50;
  Double_t minLDCDataSize = 0;
  Double_t maxLDCDataSize = 100;
  Double_t minMeanDDLDataSize = 0;
  Double_t maxMeanDDLDataSize = 50;
  Double_t minMeanLDCDataSize = 0;
  Double_t maxMeanLDCDataSize = 100;
  //  Double_t maxOccupancy = 5;

  TString histname = histo->GetName();

  if (histname.EndsWith("SSDEventType")) {
    if (histo->GetEntries()==0) {
      AliWarning("Event type histogram is empty");
    }
    else if (histo->GetBinContent(histo->FindBin(7))==0) AliWarning("No type 7 (physics) events in EventType");
  }

  if (histname.EndsWith("SSDDataSize")) {
    if (histo->GetEntries()==0) AliWarning("SSD data size histogram is empty");
    if (histo->GetMean()>maxSSDDataSize||histo->GetMean()<minSSDDataSize) AliWarning(Form("SSD mean data size is %-.2g kB", histo->GetMean()));
  }

  if (histname.EndsWith("SSDDataSizePerDDL")) {
    if (histo->GetEntries()==0) {
      AliWarning("Data size per DDL histogram is empty");
    }
    else {
      for(Int_t i = 512; i < 528; i++) {
        if(histo->GetBinContent(histo->FindBin(i))==0) {
           AliWarning(Form("Data size / DDL histogram: bin for DDL %i is empty",i));
        }
        else if(histo->GetBinContent(histo->FindBin(i))<minDDLDataSize||histo->GetBinContent(histo->FindBin(i))>maxDDLDataSize) AliWarning(Form("Data size DDL %i is %-.2g kB",i,histo->GetBinContent(histo->FindBin(i))));
     }
    }
  }

  if (histname.EndsWith("SSDDataSizePerLDC")) {
    if (histo->GetEntries()==0) {
      AliWarning("Data size per LDC histogram is empty");
    }    
    else {
      AliInfo(Form("Data size per LDC histogram has %f entries",histo->GetEntries()));
      for(Int_t i = 170; i < 178; i++) {
        if(histo->GetBinContent(histo->FindBin(i))==0) {
          AliWarning(Form("Data size / LDC histogram: bin for LDC %i is empty",i));
        }
	else if(AliITSQADataMakerRec::AreEqual(histo->GetBinContent(histo->FindBin(i)),minLDCDataSize) ||histo->GetBinContent(histo->FindBin(i))>maxLDCDataSize) AliWarning(Form("Data size LDC %i is %-.2g kB",i,histo->GetBinContent(i)));
      }
    }
  }

  if (histname.EndsWith("SSDLDCId")) {
    if (histo->GetEntries()==0) {
      AliWarning("LDC ID histogram is empty");
    }    
    else {
      for(Int_t i = 170; i < 177; i++) {
        if(histo->GetBinContent(histo->FindBin(i))==0) {
          AliWarning(Form("LDC ID histogram: No entries for LDC %i",i));
        }
        else if(histo->GetBinContent(histo->FindBin(i))!=histo->GetBinContent(histo->FindBin(i+1))) {
          AliWarning("LDC Id distribution is not uniform");
          i=176;
        }
      }
    }
  }

  if (histname.EndsWith("SSDDDLId")) {
    if (histo->GetEntries()==0) {
      AliWarning("DDL ID histogram is empty");
    }
    else {
      for(Int_t i = 512; i < 527; i++) {
        if(histo->GetBinContent(histo->FindBin(i))==0) {
          AliWarning(Form("DDL ID histogram: No entries for DDL %i",i));
        }
        else if(histo->GetBinContent(histo->FindBin(i))!=histo->GetBinContent(histo->FindBin(i+1))) {
          AliWarning("DDL Id distribution is not uniform");
          i=526;
        }
      }
    }
  }

  if (histname.Contains("SSDDataSizeLDC")) {
    if (histo->GetEntries()==0) {
      AliWarning(Form("LDC %s data size distribution is empty", histname(histname.Length()-3,3).Data()));
    }
    else if (histo->GetMean()<minMeanLDCDataSize||histo->GetMean()>maxMeanLDCDataSize) AliWarning(Form("Mean data size of LDC %s is %-.2g kB",histname(histname.Length()-3,3).Data(), histo->GetMean()));
  }

  if (histname.Contains("SSDDataSizeDDL")) {
    if (histo->GetEntries()==0) {
      AliWarning(Form("DDL %s data size distribution is empty", histname(histname.Length()-3,3).Data()));
    } 
    else if (histo->GetMean()<minMeanDDLDataSize||histo->GetMean()>maxMeanDDLDataSize) AliWarning(Form("Mean data size of DDL %s is %-.2g kB",histname(histname.Length()-3,3).Data(), histo->GetMean()));
  }

  if (histname.Contains("SSDAverageOccupancy")) {
 
    const char* side = "";
    int ladder = 0;
    int layernr = 0;

    if (histname.EndsWith("5")) layernr = 499;
    if (histname.EndsWith("6")) layernr = 599;

    for (Int_t i = 1; i < histo->GetNbinsY() + 1; i++) { //ladder/side loop
      if(i==3.*int(i/3.)){
        ladder=int(i/3.)+layernr;
        side="P side";
      }
      else if(i==3.*int(i+1/3.)){
        ladder=int((i+1)/3.)+layernr;
        side="N side";
      }

      for (Int_t j = 1; j < histo->GetNbinsX() + 1; j++) { //module loop
        //if(histo->GetBinContent(j,i)>maxOccupancy)
          // AliWarning(Form("Occupancy ladder %i, module %i, %s is %-.2f %%",ladder,j,side, histo->GetBinContent(j,i)));
      }//module loop
    }//ladder loop
  }

}


//__________________________________________________________________
Bool_t  AliITSQASSDChecker::MakeSSDImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode) {
  //create the image for raws and recpoints. In the other case, the default methodof CheckerBase class will be used  //
  gStyle->SetPalette(1,0);
  Bool_t rval=kFALSE;
  fImage=(TCanvas**)AliQAChecker::Instance()->GetDetQAChecker(0)->GetImage();

  if(task == AliQAv1::kRAWS) 
    rval = MakeSSDRawsImage(list, task,mode);
  else rval=kFALSE;
  
  return rval;
}

//_______________________________________________________________________
Bool_t AliITSQASSDChecker::MakeSSDRawsImage(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode ) {
  // MakeSSDRawsImage: raw data QA plots
  for (Int_t esIndex = 0 ; esIndex < AliRecoParam::kNSpecies ; esIndex++) {
    //printf("-------------------------> %i \n", esIndex);
    if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(esIndex)) || list[esIndex]->GetEntries() == 0) continue;
    else {
      const Char_t * title = Form("QA_%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(esIndex)) ; 
      if ( !fImage[esIndex] ) fImage[esIndex] = new TCanvas(title, title,1280,980) ;
	
      fImage[esIndex]->Clear() ; 
      fImage[esIndex]->SetTitle(title) ; 
      fImage[esIndex]->cd();
 
      //TPaveText someText(0.015, 0.015, 0.98, 0.98);
      //someText.AddText(title);
      //someText.Draw(); 
      fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
      fImage[esIndex]->Clear() ; 
      Int_t nx =2; //TMath::Nint(TMath::Sqrt(nImages));
      Int_t ny =2; // nx  ; 
	
      fImage[esIndex]->Divide(nx, ny) ; 
      TIter nexthist(list[esIndex]) ; 
      TH1* hist = NULL ;
      Int_t npad = 1 ; 
      fImage[esIndex]->cd(npad); 
      fImage[esIndex]->cd(npad)->SetBorderMode(0) ;
      while ( (hist=static_cast<TH1*>(nexthist())) ) {
	//gPad=fImage[esIndex]->cd(npad)->GetPad(npad);
	TString cln(hist->ClassName()) ; 
	if ( ! cln.Contains("TH") ) continue ;
	
	if(hist->TestBit(AliQAv1::GetImageBit())) {
	  //Printf("Histo name: %s - Class: %s",hist->GetName(),hist->ClassName());
	  hist->GetXaxis()->SetTitleSize(0.02);
	  hist->GetYaxis()->SetTitleSize(0.02);
	  hist->GetXaxis()->SetLabelSize(0.02);
	  hist->GetYaxis()->SetLabelSize(0.02);
	  if(cln.Contains("TH2")) {
	    gPad->SetRightMargin(0.15);
	    gPad->SetLeftMargin(0.05);
	    hist->SetStats(0);
	    hist->SetOption("colz") ;
	    //hist->GetListOfFunctions()->FindObject("palette")->SetLabelSize(0.025);
	    //gPad->Update();
	  }
	  hist->DrawCopy() ; 
	  fImage[esIndex]->cd(++npad) ; 
	  fImage[esIndex]->cd(npad)->SetBorderMode(0) ; 
	}
      }
      fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
    }
  }
  
  return kTRUE;
}

//__________________________________________________________________
Double_t AliITSQASSDChecker::Check(AliQAv1::ALITASK_t /*index*/, const TObjArray * list, const AliDetectorRecoParam * /*recoParam*/) { 
  // main checker method 
  AliDebug(AliQAv1::GetQADebugLevel(),Form("AliITSQASSDChecker called with offset: %d\n", fSubDetOffset));

  AliInfo(Form("AliITSQASSDChecker called with offset: %d\n", fSubDetOffset) );
  //cout<<"(AliITSQASSDChecker::Check): List name "<<list->GetName()<<endl;
  Double_t test = 0.0  ;
  Int_t count = 0 ;
  TString listname = list->GetName();

  if (list->GetEntries() == 0){
    test = 1. ; // nothing to check
  }
  else {

    TIter next(list) ;
    TH1 * hdata ;
    count = 0 ;
    while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
      if (hdata) {
	TString histname = hdata->GetName();
	if(!histname.Contains("fHistSSD")) continue;
        Double_t rv = 0.;
        if(hdata->GetEntries()>0) {
           rv = 1;

	   //if(histname.Contains("PerDDL")) cout << "(AliITSQASSDChecker::Check) " << histname << " has " << hdata->GetEntries() << " entries. Mean: " << hdata->GetMean() << endl;
       
       //    if(hdata->GetMean()>0&&!histname.Contains("_Ladder")) cout << "(AliITSQASSDChecker::Check) " << histname << " not empty! " << hdata->GetEntries() << " entries. Mean: " << hdata->GetMean() << endl;
        }

    //    if (listname.Contains("Raws")) CheckRaws(hdata);
   //     if (listname.Contains("RecPoints")) CheckRecPoints(hdata);

        //AliDebug(AliQAv1::GetQADebugLevel(), Form("%s -> %f", hdata->GetName(), rv)) ;
	//cout<<hdata->GetName()<<" - "<<rv<<endl;
        count++ ;
        test += rv ;
      }
      else{
        AliError("Data type cannot be processed") ;
      }
    }
    if (count != 0) {
      if (AliITSQADataMakerRec::AreEqual(test,0.)) {
        AliWarning("Histograms are there, but they are all empty: setting flag to kWARNING");
        test = 0.5;  //upper limit value to set kWARNING flag for a task
      }
      else {
        test /= count ;
      }
    }
  }
  
  //AliDebug(AliQAv1::GetQADebugLevel(), Form("Test Result = %f", test)) ;
  //cout<<"Test result: "<<test<<endl;

  return test ;

  //return 0.;


}

//__________________________________________________________________
void AliITSQASSDChecker::SetTaskOffset(Int_t TaskOffset){
  // defines offset for SSD
  fSubDetOffset = TaskOffset;
}

//__________________________________________________________________
void AliITSQASSDChecker::SetStepBit(const Double_t *steprange) {
  // defines step range
  fStepBitSSD = new Double_t[AliQAv1::kNBIT];
  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fStepBitSSD[bit]=steprange[bit];
    }
}

//__________________________________________________________________
void  AliITSQASSDChecker::SetSSDLimits(const Float_t *lowvalue, const Float_t * highvalue){
  // defines 
  fLowSSDValue = new Float_t[AliQAv1::kNBIT];
  fHighSSDValue= new Float_t[AliQAv1::kNBIT];

  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fLowSSDValue[bit]=lowvalue[bit];
      fHighSSDValue[bit]= highvalue[bit];
    }

}
