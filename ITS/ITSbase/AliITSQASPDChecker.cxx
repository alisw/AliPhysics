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

/* $Id $ */

// *****************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  P. Cerello Apr 2008
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"
#include "TH2F.h"
#include "TString.h"
#include "TList.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TFrame.h"

// --- AliRoot header files ---
#include "AliQACheckerBase.h"
#include "AliITSQASPDChecker.h"
#include "AliITSQADataMakerRec.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSTriggerConditions.h"

ClassImp(AliITSQASPDChecker)
//__________________________________________________________________
AliITSQASPDChecker::AliITSQASPDChecker() : 
  TObject(),
  fSubDetOffset(0), 
  fStepBitSPD(NULL),
  fLowSPDValue(NULL),
  fHighSPDValue(NULL),
  fImage(NULL)
{
  // default contructor
}
//__________________________________________________________________
AliITSQASPDChecker& AliITSQASPDChecker::operator = (const AliITSQASPDChecker& qac ) 
{
  // Equal operator.
  this->~AliITSQASPDChecker();
  new(this) AliITSQASPDChecker(qac);
  return *this;
}
//__________________________________________________________________
AliITSQASPDChecker::~AliITSQASPDChecker() {
  // destructor
  if(fStepBitSPD) delete[] fStepBitSPD ;
  if(fLowSPDValue)delete[]fLowSPDValue;
  if(fHighSPDValue) delete[]fHighSPDValue;
  if(fImage) delete[]fImage;
} 

//__________________________________________________________________
Double_t AliITSQASPDChecker::Check(AliQAv1::ALITASK_t index, TObjArray * list, const AliDetectorRecoParam * /*recoParam*/)
{
  //
  // General methods for SPD Cheks to be used in RAWS and REC ALITASK_t
  //

  AliDebug(2, Form("AliITSQASPDChecker called with offset: %d\n", fSubDetOffset));

  Double_t test = 0.0;
  Int_t count = 0;
  // Checks for ALITASK_t AliQAv1::kRAW
  if(index == AliQAv1::kRAW) {
    return CheckRawData(list);
  } else {
    if (list->GetEntries() == 0) {
      test = 1.; // nothing to check
    }
    else {
      TIter next(list);
      TH1 * hdata;
      count = 0;
      while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
	if (hdata) {
	  TString histName = hdata->GetName();
	  if (!histName.Contains("_SPD")) continue;
	  Double_t rv = 0.;
	  if (hdata->GetEntries()>0) rv = 1;
	  if (histName.Contains("LayPattern")) {
	    if (hdata->GetBinContent(1)) {
	      Double_t ratio=hdata->GetBinContent(2)/hdata->GetBinContent(1);
	      AliDebug(2, Form("%s: ratio RecPoints lay2 / lay1 = %f", hdata->GetName(), ratio));
	    }
	    else
	      AliDebug(AliQAv1::GetQADebugLevel(), "No RecPoints in lay1");
	  }
	  else if(histName.Contains("ModPattern")) {
	    Int_t ndead=0;
	    for(Int_t ibin=0;ibin<hdata->GetNbinsX();ibin++) {
	      if(histName.Contains("SPD1") && ibin<80 && hdata->GetBinContent(ibin+1)>0) ndead++;
	      if(histName.Contains("SPD2") && ibin>79 && hdata->GetBinContent(ibin+1)>0) ndead++;
	    }
	    AliDebug(2, Form("%s: Entries = %d  number of empty modules = %d", 
			     hdata->GetName(),(Int_t)hdata->GetEntries(),ndead));
	  }
	  else if(histName.Contains("SizeYvsZ")) {
	    Double_t meanz=hdata->GetMean(1);
	    Double_t meany=hdata->GetMean(2);
	    Double_t rmsz=hdata->GetRMS(1);
	    Double_t rmsy=hdata->GetRMS(2);
	    AliDebug(AliQAv1::GetQADebugLevel(), Form("%s: Cluster sizeY mean = %f  rms = %f", hdata->GetName(),meany,rmsy));
	    AliDebug(AliQAv1::GetQADebugLevel(), Form("%s: Cluster sizeZ mean = %f  rms = %f", hdata->GetName(),meanz,rmsz));
	  }
	  else if(histName.Contains("SPDMultiplicity")) {
	    AliDebug(2, Form("%s: Events = %d  mean = %f  rms = %f",
			     hdata->GetName(),(Int_t)hdata->GetEntries(),hdata->GetMean(),hdata->GetRMS()));}

	  // else AliDebug(AliQAv1::GetQADebugLevel(), Form("%s -> %f", hdata->GetName(), rv));
	  count++;
	  test += rv;
	}
	else {
	  AliError("Data type cannot be processed") ;
	}
      }

      if (count != 0) {
	if (AliITSQADataMakerRec::AreEqual(test,0)) {
	  AliWarning("Histograms are there, but they are all empty: setting flag to kWARNING");
	  test = fHighSPDValue[AliQAv1::kWARNING];  //upper limit value to set kWARNING flag for a task
	}
	else {
	  test /= count;
	}
      }
    }
  }
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Test Result = %f", test));
  return test ;

}
//__________________________________________________________________
Double_t AliITSQASPDChecker::CheckRawData(const TObjArray * list) {
  //
  // Checks on the raw data histograms [ preliminary version ]
  // The output of this method is the fraction of SPD histograms which are processed by the checker. 
  // The methods returns fHighSPDValue[AliQAv1::kFATAL] in case of data format errors or MEB errors
  // 
  // A. Mastroserio

  Double_t test =0;

  // basic checks on input data
  if(!list) {
    AliError("NO histogram list for RAWS");
    return test;
  }

  if(list->GetEntries() == 0) {
    AliWarning("No histograms in RAW list \n");
    return test;
  }

  // loop over the raw data histograms
  TIter next(list);
  TH1 * hdata;
  Double_t totalHistos = 0;
  Double_t goodHistos = 0; // number of histograms which passed the checks
  Double_t response =0;
  Bool_t fatalProblem = kFALSE;

  while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
    if (hdata) {
      TString histName = hdata->GetName();
      if(!histName.Contains("SPD")) continue;
      totalHistos++;
      // data format error
      if(histName.Contains("SPDErrorsAll")){
	TPaveText *p = ((TPaveText*)hdata->GetListOfFunctions()->FindObject("TPave"));	    
	if(hdata->Integral(0,hdata->GetNbinsX())>0){
	  Bool_t isHighMult = kFALSE;
	  Bool_t isDataCorrupted=kFALSE;
	  for(Int_t ieq=0; ieq<20; ieq++){
	    if(hdata->GetBinContent(ieq+1,17+1)>0 && hdata->GetBinContent(ieq+1,20+1)>0) isHighMult = kTRUE;
	    for(Int_t iErr=1; iErr<20; iErr++){
	      if(iErr==20 || iErr==17) continue;
	      if(hdata->GetBinContent(ieq+1,iErr+1)>0) isDataCorrupted=kTRUE; 
	    }
	  }
	  if(isHighMult && !isDataCorrupted) {
	    response = fHighSPDValue[AliQAv1::kWARNING];
	    if(p) {
	      p->Clear();
	      p->SetFillColor(kOrange);
	      p->AddText("High occupancy in a chip detected (only error types 17,20 and 0 are present) ");
	      p->AddText("ONLY if OTHER error types are present CALL the expert");
	    }
	  } else if(isDataCorrupted) {
	    response = fHighSPDValue[AliQAv1::kFATAL];
	    fatalProblem=kTRUE;
	    if(p) {
	      p->Clear();
	      p->SetFillColor(kRed);
	      p->AddText("Data Format NOT OK. Please call the expert!");		
	    }
	  }
	  continue;
	} // if errors 
	else {
	  if(p) {
	    p->Clear();
	    p->SetFillColor(kGreen);
	    p->AddText("OK");
	  }
	}
      } // data format error

      // MEB error
      else if(histName.Contains("MEB")){
	TPaveText *pmeb = ((TPaveText*)hdata->GetListOfFunctions()->FindObject("TPave"));
	if(hdata->GetEntries()>0){
	  response = fHighSPDValue[AliQAv1::kFATAL];
	  fatalProblem=kTRUE;
	  if(pmeb) {
	    pmeb->Clear();
	    pmeb->SetFillColor(kRed);
	    pmeb->AddText("MEB problem could be present. Please check if SPD is in READY state.");
	    pmeb->AddText("If SPD is in -READY- state, please notify it to the expert."); 
	  }
	  continue;

	} else {
	  if(pmeb) {
	    pmeb->Clear();
	    pmeb->SetFillColor(kGreen);
	    pmeb->AddText("OK");
	  }
	}   

      }
      goodHistos++;
    }
  }
  if(!fatalProblem) response = goodHistos/totalHistos;
  // printf("n histos %f - good ones %f ----> ratio %f , fatal response %i\n",totalHistos,goodHistos,goodHistos/totalHistos,(Int_t)fatalProblem);
  return response;
}

//__________________________________________________________________
void AliITSQASPDChecker::SetTaskOffset(Int_t TaskOffset)
{
  // Offset for SPD within ITS QA
  fSubDetOffset = TaskOffset;
}

//__________________________________________________________________
void AliITSQASPDChecker::SetStepBit(const Double_t *steprange) 
{
  // Step bit for SPD within ITS QA
  fStepBitSPD = new Double_t[AliQAv1::kNBIT];
  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fStepBitSPD[bit]=steprange[bit];
    }
}

//__________________________________________________________________
void  AliITSQASPDChecker::SetSPDLimits(const Float_t *lowvalue, const Float_t * highvalue)
{
  // SPD limints for QA bit within general ITS QA
  fLowSPDValue = new Float_t[AliQAv1::kNBIT];
  fHighSPDValue= new Float_t[AliQAv1::kNBIT];

  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      fLowSPDValue[bit]=lowvalue[bit];
      fHighSPDValue[bit]= highvalue[bit];
    }

}

//__________________________________________________________________
void AliITSQASPDChecker::MaskFastOrChips(TH2F *hFOmap, Int_t iLayer)
{
  // Empty masked FastOr Chips from OCDB
  AliCDBEntry *pitEntry = 0x0;       
  AliITSTriggerConditions *triCd = 0x0;
  Int_t sect = 999;
  Int_t chipbin = 999;
  //
  if(!hFOmap) {
    AliWarning("Unable to retrieve histogram!");
    return;
  }
  
  // Read PIT trigger conditions
  AliCDBManager *man = AliCDBManager::Instance();
  if(man) {
    man->SetDefaultStorage(gSystem->Getenv("AMORE_CDB_URI"));
    pitEntry = man->Get("TRIGGER/SPD/PITConditions", AliQAChecker::Instance()->GetRunNumber());
  }
  if(!pitEntry) {
    AliWarning("Trigger conditions retrieval from OCDB failed!");
    return;
  }
  triCd = (AliITSTriggerConditions*)pitEntry->GetObject();
  
  //exclude masked chips
  for (UInt_t ieq=0; ieq<20; ieq++) {
    for (UInt_t ihs=2*iLayer; ihs<(2+4*iLayer); ihs++) {
      for (UInt_t ichip=0; ichip<10; ichip++) {
	if(!triCd->IsChipActive(ieq,ihs,ichip)) {
	  if(ieq<10) {sect = ieq; chipbin = (19-ichip)+1;} //side A
	  else {sect = ieq-10; chipbin = ichip+1;} //side C
	  if(iLayer==0)   hFOmap->SetBinContent(chipbin,1+sect*2+ihs,0.); //inner layer
	  if(iLayer==1)   hFOmap->SetBinContent(chipbin,1+sect*4+(ihs-2),0.); //outer layer
	}
      }
    }
  }
  
}

//__________________________________________________________________
Bool_t  AliITSQASPDChecker::MakeSPDImage( TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode)
{
  //create the image for raws and recpoints. In the other case, the default methodof CheckerBase class will be used

  Bool_t val=kFALSE;

  fImage=(TCanvas**)AliQAChecker::Instance()->GetDetQAChecker(0)->GetImage();

  switch(task)
    {
    case AliQAv1::kRAWS:{
      val = MakeSPDRawsImage(list, task,mode);
    }
      break;
    case AliQAv1::kRECPOINTS:;
    case AliQAv1::kHITS:; 
    case AliQAv1::kESDS:; 
    case AliQAv1::kDIGITS:;
    case AliQAv1::kDIGITSR:;
    case AliQAv1::kSDIGITS:;
    case AliQAv1::kTRACKSEGMENTS:;
    case AliQAv1::kRECPARTICLES:; 
    default:
      {
	//AliQAChecker::Instance()->GetDetQAChecker(0)->MakeImage(list,task,mode);
	val = kFALSE;
      }
    break;
    case AliQAv1::kNULLTASKINDEX:; case  AliQAv1::kNTASKINDEX: 
      {AliWarning(Form("No histograms for these tasks ( %s ) \n", AliQAv1::GetTaskName(task).Data())); val = kFALSE;}
      break;
    }
  return val; 
}
//_______________________________________________________________________
Bool_t AliITSQASPDChecker::MakeSPDRawsImage(TObjArray ** list, AliQAv1::TASKINDEX_t task, AliQAv1::MODE_t mode )
{
  //
  // create layout of the histograms used in the DQM
  //

  // some style settings
  gStyle->SetOptStat(10);  
  gStyle->SetStatY(0.97);
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.165);
  gStyle->SetTitleX(0.45);
  gStyle->SetTitleOffset(1.2);

  for (Int_t esIndex = 0 ; esIndex < AliRecoParam::kNSpecies ; esIndex++) {
    //printf("-------------------------> %i \n", esIndex);
    if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(esIndex)) || list[esIndex]->GetEntries() == 0) 
      {printf ("Nothing for %s \n", AliRecoParam::GetEventSpecieName(esIndex)); continue;}
    else{
      const Char_t * title = Form("QA_%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(esIndex)) ; 
      if ( !fImage[esIndex] ) {
	fImage[esIndex] = new TCanvas(title, title,6000,3200) ;
      }

      fImage[esIndex]->Clear() ; 
      fImage[esIndex]->SetTitle(title) ; 
      fImage[esIndex]->cd();

      TPaveText someText(0.015, 0.015, 0.98, 0.98);
      someText.AddText(title);
      someText.Draw(); 
      fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
      fImage[esIndex]->Clear() ; 
      Int_t nx =3; 
      Int_t ny =2; 

      fImage[esIndex]->Divide(nx, ny) ; 

      TH1* hist = NULL ;
      Int_t npad = 1 ; 
      fImage[esIndex]->cd(npad); 
      fImage[esIndex]->cd(npad)->SetBorderMode(0) ;

      TIter next(list[esIndex]);

      while ( (hist=static_cast<TH1*>(next())) ) {
	//gPad=fImage[esIndex]->cd(npad)->GetPad(npad);
	if(!hist->TestBit(AliQAv1::GetImageBit())) continue;
	TString name(hist->GetName());
	if(name.Contains("SPDErrorsAll")) {
	  fImage[esIndex]->cd(1) ;    
	  TPaveText *pImErrs = ((TPaveText*)hist->GetListOfFunctions()->FindObject("TPave"));
	  if(pImErrs) {
	    if(pImErrs->GetFillColor()!=kGreen){
	      gPad->SetFillColor(pImErrs->GetFillColor());
	      gPad->SetFillStyle(3001);
	      gPad->SetFrameFillColor(19);
	    }    
	  }
	  gPad->SetBorderMode(0) ;  
	  gPad->SetRightMargin(0.05);
	  gPad->SetLeftMargin(0.1);
	  gPad->SetTopMargin(0.25);
	  gPad->SetGridx();
	  gPad->SetGridy();
	  hist->SetOption("text");
	  hist->GetYaxis()->SetLabelSize(0.030);
	  hist->GetYaxis()->SetNdivisions(211);
	  hist->GetXaxis()->SetLabelSize(0.030);
	  hist->GetXaxis()->SetNdivisions(20);
	  hist->DrawCopy();
	  gPad->Update();
	}     
	if(name.Contains("MEB")) {
	  fImage[esIndex]->cd(4) ; 
	  TPaveText *pImMeb = ((TPaveText*)hist->GetListOfFunctions()->FindObject("TPave"));
	  if(pImMeb) {
	    if(pImMeb->GetFillColor()!=kGreen){
	      gPad->SetFillColor(pImMeb->GetFillColor());
	      gPad->SetFillStyle(3001);
	      gPad->SetFrameFillColor(19);
	    }
	  }
	  gPad->SetBorderMode(0) ;  
	  gPad->SetGridy();
	  gPad->SetGridx();
	  gPad->SetRightMargin(0.05);
	  gPad->SetLeftMargin(0.1);
	  gPad->SetBottomMargin(0.1);
	  gPad->SetTopMargin(0.27);
	  hist->GetXaxis()->SetLabelSize(0.033);
	  hist->SetOption("text") ;
	  hist->DrawCopy();  
	}     
	if(name.Contains("SPDFastOrCorrelation")){
	  fImage[esIndex]->cd(7) ; 
	  gPad->SetBorderMode(0) ;  
	  hist->SetOption("") ;
	  hist->DrawCopy();               
	}
	if(name.Contains("SPDHitMapStaveChipInner")){
	  fImage[esIndex]->cd(2) ;
	  gStyle->SetStatX(0.81);
	  gStyle->SetStatY(0.942);
	  gPad->SetBorderMode(0) ;  
	  gPad->SetRightMargin(0.25);
	  gPad->SetGridx();
	  gPad->SetGridy();
	  hist->SetObjectStat(0);
	  hist->SetOption("colz") ;
	  hist->DrawCopy();               
	}
	if(name.Contains("SPDHitMapStaveChipOuter")){
	  fImage[esIndex]->cd(3) ;
	  gStyle->SetStatX(0.81);
	  gStyle->SetStatY(0.942);
	  gPad->SetBorderMode(0) ;  
	  gPad->SetRightMargin(0.25);
	  gPad->SetGridx();
	  gPad->SetGridy();
	  hist->SetObjectStat(0);
	  hist->SetOption("colz") ;
	  hist->DrawCopy();   
	  gPad->Update();
	}
	if(name.Contains("SPDFastOrMapStaveChip")){
	  fImage[esIndex]->cd(8) ; 
	  gPad->SetBorderMode(0) ;  
	  gPad->SetRightMargin(0.1);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetBottomMargin(0.15);
	  gPad->SetTopMargin(-0.2);
	  gPad->SetGridy();
	  gPad->SetGridx();
	  hist->SetObjectStat(0);
	  hist->SetOption("colz") ;
	  hist->DrawCopy("colz");   
	  TH2F *h2 =  (TH2F*)(gPad->GetListOfPrimitives()->At(0)); 
	  for(Int_t i=0; i<h2->GetListOfFunctions()->GetEntries(); i++){
	    TString cname = h2->GetListOfFunctions()->At(i)->ClassName();
	    if(cname.Contains("TPaletteAxis")){
	      TPaletteAxis *palette = (TPaletteAxis*)(h2->GetListOfFunctions()->At(i));
	      if(palette) palette->SetLabelSize(0.02);
	    }
	  }
	  hist->DrawCopy("colz");   
	  //gPad->Update();
	}
	if(name.Contains("SPDFastOrMapStaveInner")){
	  fImage[esIndex]->cd(5) ;
	  gStyle->SetStatX(0.81);
	  gStyle->SetStatY(0.942);
	  gPad->SetBorderMode(0) ;  
	  gPad->SetRightMargin(0.25);
	  gPad->SetGridx();
	  gPad->SetGridy();
	  hist->SetObjectStat(0);
	  hist->SetOption("colz") ;
	  MaskFastOrChips(dynamic_cast<TH2F*>(hist),0);
	  hist->DrawCopy();
	}
	if(name.Contains("SPDFastOrMapStaveOuter")){
	  fImage[esIndex]->cd(6) ; 
	  gStyle->SetStatX(0.81);
	  gStyle->SetStatY(0.942);
	  gPad->SetBorderMode(0) ;  
	  gPad->SetRightMargin(0.25);
	  gPad->SetGridx();
	  gPad->SetGridy();
	  hist->SetObjectStat(0);
	  hist->SetOption("colz") ;
	  MaskFastOrChips(dynamic_cast<TH2F*>(hist),1);        
	  hist->DrawCopy();
	}

      }

      fImage[esIndex]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), AliQAv1::GetModeName(mode), AliQAChecker::Instance()->GetRunNumber(), AliQAv1::GetImageFileFormat()), "ps") ; 
      //fImage[esIndex]->SaveAs("image.png");
    }
  }
  return kTRUE;
}

