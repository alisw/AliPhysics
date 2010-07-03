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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Example: fill pedestal with Gaussian noise                             //
//                                                                        //
//  AliTRDCalibPadStatus ped;                                             //
//  ped.TestEvent(numberofevent);                                         //
//                                                                        //
//  // Method without histo                                               //
//  ped.Analyse();                                                        //
//                                                                        //
//  // Create the histo of the AliTRDCalROC                               //
//  TH2F * histo2dm = ped.GetCalRocMean(0,kFALSE)->MakeHisto2D();         //
//  histo2dm->Scale(10.0);                                                //
//  TH1F * histo1dm = ped.GetCalRocMean(0,kFALSE)->MakeHisto1D();         //
//  histo1dm->Scale(10.0);                                                //
//  TH2F * histo2ds = ped.GetCalRocSquares(0,kFALSE)->MakeHisto2D();      //
//  histo2ds->Scale(10.0);                                                //
//  TH1F * histo1ds = ped.GetCalRocSquares(0,kFALSE)->MakeHisto1D();      //
//  histo1ds->Scale(10.0)                                                 //
//                                                                        //
//  // Draw output                                                        //
//  TCanvas* c1 = new TCanvas;                                            //
//  c1->Divide(2,2);                                                      //
//  c1->cd(1);                                                            //
//  histo2dm->Draw("colz");                                               //
//  c1->cd(2);                                                            //
//  histo1dm->Draw();                                                     //
//  c1->cd(3);                                                            //
//  histo2ds->Draw("colz");                                               //
//  c1->cd(4);                                                            //
//  histo1ds->Draw();                                                     //
//                                                                        //
//  // Method with histo                                                  //
//  ped.AnalyseHisto();                                                   //
//                                                                        //
//  // Take the histo                                                     //
//  TH1F *histo = ped.GetHisto(31);                                       //
//  histo->SetEntries(1);                                                 //
//  histo->Draw();                                                        //
//
// Authors:
//   R. Bailhache (R.Bailhache@gsi.de, rbailhache@ikf.uni-frankfurt.de)
//   J. Book (jbook@ikf.uni-frankfurt.de)
//                                                                                                    //
////////////////////////////////////////////////////////////////////////////


//Root includes
#include <TObjArray.h>
#include <TH2F.h>
#include <TString.h>
#include <TMath.h>
#include <TRandom.h>

//#include <TRandom.h>
#include <TDirectory.h>
#include <TFile.h>

//AliRoot includes
#include <AliMathBase.h>
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"

//header file
#include "AliLog.h"
#include "AliTRDCalibPadStatus.h"
#include "AliTRDrawStreamBase.h"
#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"
#include "./Cal/AliTRDCalROC.h"
#include "./Cal/AliTRDCalPadStatus.h"
#include "./Cal/AliTRDCalDet.h"
#include "./Cal/AliTRDCalPad.h"
#include "./Cal/AliTRDCalSingleChamberStatus.h"

#include "AliTRDrawFastStream.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDdigitsParam.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDarraySignal.h"
#include "AliTRDarrayADC.h"
#include "AliTRDfeeParam.h"

#include "AliTRDrawStream.h"

#ifdef ALI_DATE
#include "event.h"
#endif

ClassImp(AliTRDCalibPadStatus) /*FOLD00*/

//_____________________________________________________________________
AliTRDCalibPadStatus::AliTRDCalibPadStatus() : /*FOLD00*/
  TObject(),
  fGeo(0),
  fAdcMin(0),
  fAdcMax(21),
  fDetector(-1),
  fNumberOfTimeBins(0),
  fCalRocArrayMean(540),
  fCalRocArrayRMS(540),
  fCalRocArrayMeand(540),
  fCalRocArrayRMSd(540),
  fHistoArray(540)
{
    //
    // default constructor
    //

  fGeo = new AliTRDgeometry();

}

//_____________________________________________________________________
AliTRDCalibPadStatus::AliTRDCalibPadStatus(const AliTRDCalibPadStatus &ped) : /*FOLD00*/
  TObject(ped),
  fGeo(0),
  fAdcMin(ped.GetAdcMin()),
  fAdcMax(ped.GetAdcMax()),
  fDetector(ped.fDetector),
  fNumberOfTimeBins(ped.fNumberOfTimeBins),
  fCalRocArrayMean(540),
  fCalRocArrayRMS(540),
  fCalRocArrayMeand(540),
  fCalRocArrayRMSd(540),
  fHistoArray(540)
{
    //
    // copy constructor
    //
    for (Int_t idet = 0; idet < 540; idet++){
	const AliTRDCalROC *calRocMean  = (AliTRDCalROC*)ped.fCalRocArrayMean.UncheckedAt(idet);
	const AliTRDCalROC *calRocRMS   = (AliTRDCalROC*)ped.fCalRocArrayRMS.UncheckedAt(idet);
	const AliTRDCalROC *calRocMeand = (AliTRDCalROC*)ped.fCalRocArrayMeand.UncheckedAt(idet);
	const AliTRDCalROC *calRocRMSd  = (AliTRDCalROC*)ped.fCalRocArrayRMSd.UncheckedAt(idet);
	const TH2F         *hped        = (TH2F*)ped.fHistoArray.UncheckedAt(idet);
    
	if ( calRocMean != 0x0 ) fCalRocArrayMean.AddAt(new AliTRDCalROC(*calRocMean), idet);
	if ( calRocRMS != 0x0 )  fCalRocArrayRMS.AddAt(new AliTRDCalROC(*calRocRMS), idet);

	if ( calRocMeand != 0x0 ) fCalRocArrayMeand.AddAt(new AliTRDCalROC(*calRocMeand), idet);
	if ( calRocRMSd != 0x0 )  fCalRocArrayRMSd.AddAt(new AliTRDCalROC(*calRocRMSd), idet);

	if ( hped != 0x0 ){
	  TH2F *hNew = new TH2F(*hped);
	  hNew->SetDirectory(0);
	  fHistoArray.AddAt(hNew,idet);
	}
	
    }
    if (fGeo) {
      delete fGeo;
    }
    fGeo = new AliTRDgeometry();
}

//_____________________________________________________________________
AliTRDCalibPadStatus& AliTRDCalibPadStatus::operator = (const  AliTRDCalibPadStatus &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTRDCalibPadStatus(source);

  return *this;
}
//_____________________________________________________________________
AliTRDCalibPadStatus::~AliTRDCalibPadStatus() /*FOLD00*/
{
  //
  // destructor
  //
  fCalRocArrayMean.Delete();
  fCalRocArrayRMS.Delete();
  fCalRocArrayMeand.Delete();
  fCalRocArrayRMSd.Delete();
  fHistoArray.Delete();
  if (fGeo) {
    delete fGeo;
  }
}
//_____________________________________________________________________
void AliTRDCalibPadStatus::Destroy()
{
  //
  // Destroy
  //
  fCalRocArrayMean.Delete();
  fCalRocArrayRMS.Delete();
  fCalRocArrayMeand.Delete();
  fCalRocArrayRMSd.Delete();
  fHistoArray.Delete();
}
//_____________________________________________________________________
Int_t AliTRDCalibPadStatus::UpdateHisto(const Int_t icdet, /*FOLD00*/
					const Int_t icRow,
					const Int_t icCol,
					const Int_t csignal,
					const Int_t crowMax,
					const Int_t ccold,
					const Int_t icMcm)
{
  //
  // Signal filling methode 
  //
  Int_t nbchannel = icRow+icCol*crowMax;

  // now the case of double read channel
  if(ccold > 0){
    nbchannel = (((ccold-1)*8+ icMcm)*crowMax+icRow)+144*crowMax;
    //printf("nbchannel %d, ccold %d, icMcm %d, crowMax %d, icRow %d\n",nbchannel,ccold,icMcm,crowMax,icRow);
  }
  
  // fast filling methode.
  // Attention: the entry counter of the histogram is not increased
  //            this means that e.g. the colz draw option gives an empty plot
  
  Int_t bin = 0;
  
  if ( !(((Int_t)csignal>=fAdcMax ) || ((Int_t)csignal<fAdcMin)) )
    bin = (nbchannel+1)*(fAdcMax-fAdcMin+2)+((Int_t)csignal-fAdcMin+1);
  
  //GetHisto(icdet,kTRUE)->Fill(csignal,nbchannel);
  
  GetHisto(icdet,kTRUE)->GetArray()[bin]++;
  
  return 0;
}
//_____________________________________________________________________
Int_t AliTRDCalibPadStatus::UpdateHisto2(const Int_t icdet, /*FOLD00*/
					 const Int_t icRow,
					 const Int_t icCol,
					 const Int_t csignal,
					 const Int_t crowMax,
					 const Int_t ccold,
					 const Int_t icMcm,
					 const Int_t icRob
					 )
{
  //
  // Signal filling methode 
  //
  Int_t nbchannel = icRow+icCol*crowMax;
  Int_t mCm = icMcm%4;
  Int_t rOb = icRob%2;

  // now the case of double read channel
  if(ccold > 0){
    nbchannel = (((ccold-1)*8+ (mCm+rOb*4))*crowMax+icRow)+144*crowMax;
    //printf("nbchannel %d, ccold %d, icMcm %d, crowMax %d, icRow %d\n",nbchannel,ccold,icMcm,crowMax,icRow);
  }
  
  // fast filling methode.
  // Attention: the entry counter of the histogram is not increased
  //            this means that e.g. the colz draw option gives an empty plot
  
  Int_t bin = 0;
  
  if ( !(((Int_t)csignal>=fAdcMax ) || ((Int_t)csignal<fAdcMin)) )
    bin = (nbchannel+1)*(fAdcMax-fAdcMin+2)+((Int_t)csignal-fAdcMin+1);

  //GetHisto(icdet,kTRUE)->Fill(csignal,nbchannel);

  GetHisto(icdet,kTRUE)->GetArray()[bin]++;
  
  return 0;
}
//_____________________________________________________________________
Int_t AliTRDCalibPadStatus::ProcessEvent(AliTRDrawStreamBase *rawStream, Bool_t nocheck)
{
  //
  // Event Processing loop - AliTRDRawStreamCosmic
  // 0 time bin problem or zero suppression
  // 1 no input
  // 2 input
  // 

  //
  // Raw version number: 
  // [3,31] non zero suppressed
  // 2,4 and [32,63] zero suppressed 
  //

  Int_t withInput = 1;

  rawStream->SetSharedPadReadout(kTRUE);

  if(!nocheck) {

    // Check the raw version and if all have the same number of timebins. 

    while (rawStream->Next()) {

      Int_t rawversion = rawStream->GetRawVersion();                     //  current raw version
      //printf("Raw version is %d\n",rawversion);

      // Could eventually change, have to check with time    
      if((rawversion < 3) || (rawversion > 31)) {
	AliInfo(Form("this is not no-zero-suppressed data, the version is %d",rawversion));
	return 0;
      }
      Int_t idetector  = rawStream->GetDet();                            //  current detector
      Int_t iRow       = rawStream->GetRow();                            //  current row
      Int_t iRowMax    = rawStream->GetMaxRow();                         //  current rowmax
      Int_t iCol       = rawStream->GetCol();                            //  current col
      Int_t iADC       = 21-rawStream->GetADC();                         //  current ADC
      
      // It goes in the opposite direction
      Int_t col        = 0;
      if(iADC == 1) col = 1;
      else {
	col = TMath::Max(0,(Int_t)(iADC-19));
	if(col > 0) col++;
      }
      Int_t mcm        = (Int_t)(iCol/18);                               //  current group of 18 col pads
      if(col > 1) mcm -= 1;      
      if(col ==1) mcm += 1;

      // printf to check
      //Bool_t shared = rawStream->IsCurrentPadShared();                  
      //printf("ADC %d, iCol %d, col %d, mcm %d, shared %d\n",iADC,iCol,col,mcm,(Int_t)shared);

      // Take the signal
      Int_t *signal    = rawStream->GetSignals();                        //  current ADC signal
      Int_t nbtimebin  = rawStream->GetNumberOfTimeBins();               //  number of time bins read from data

      if((fDetector != -1) && (nbtimebin != fNumberOfTimeBins)) {
	AliInfo(Form("the number of time bins is %d, is different from the previous one %d",nbtimebin,fNumberOfTimeBins));
      	return 0;
      }
      fNumberOfTimeBins = nbtimebin;
      fDetector         = idetector;      

      for(Int_t k = 0; k < fNumberOfTimeBins; k++){
	if(signal[k]>0 && iCol != -1) UpdateHisto(idetector,iRow,iCol,signal[k],iRowMax,col,mcm);
      }
      
      withInput = 2;
    }
  }
  else {

    while (rawStream->Next()) {
    
      Int_t idetector  = rawStream->GetDet();                            //  current detector
      Int_t iRow       = rawStream->GetRow();                            //  current row
      Int_t iRowMax    = rawStream->GetMaxRow();                         //  current rowmax
      Int_t iCol       = rawStream->GetCol();                            //  current col
      Int_t iADC       = 21-rawStream->GetADC();                            //  current ADC

      // It goes in the opposite direction      
      Int_t col        = 0;
      if(iADC == 1) col = 1;
      else {
	col = TMath::Max(0,(Int_t)(iADC-19));
	if(col > 0) col++;
      }
      Int_t mcm        = (Int_t)(iCol/18);                               //  current group of 18 col pads
      if(col > 1) mcm -= 1;      
      if(col ==1) mcm += 1;

      // Take the signal
      Int_t *signal    = rawStream->GetSignals();                        //  current ADC signal
      Int_t nbtimebin = rawStream->GetNumberOfTimeBins();               //  number of time bins read from data
      
      
      //printf("det %d, row %d, signal[0] %d, signal[1] %d, signal [2] %d\n", idetector, iRow, signal[0], signal[1], signal[2]);
 
      for(Int_t k = 0; k < nbtimebin; k++){
	if(signal[k]>0 && iCol != -1) {
	  UpdateHisto(idetector,iRow,iCol,signal[k],iRowMax,col,mcm);
	  //printf("Update with det %d, row %d, col %d, signal %d, rowmax %d, col %d, mcm %d\n",idetector,iRow,iCol,signal[n],iRowMax,col,mcm);
	}
      }
      
      withInput = 2;
    }
  }
  
  return withInput;
}
//_____________________________________________________________________
Int_t AliTRDCalibPadStatus::ProcessEvent(AliRawReader *rawReader, Bool_t nocheck)
{
  //
  //  Event processing loop - AliRawReader
  //

  Int_t result;
  
  rawReader->Select("TRD");
  
  AliTRDrawStreamBase *pstream = AliTRDrawStreamBase::GetRawStream(rawReader);
 
  result = ProcessEvent(pstream, nocheck);

  delete pstream;

  return result;
}

//_________________________________________________________________________
Int_t AliTRDCalibPadStatus::ProcessEvent(
#ifdef ALI_DATE
					  const eventHeaderStruct *event,
					  Bool_t nocheck
#else
					  const eventHeaderStruct* /*event*/,
					  Bool_t /*nocheck*/
	    
#endif 
					  )
{
  //
  //  process date event
  //
#ifdef ALI_DATE
    AliRawReader *rawReader = new AliRawReaderDate((void*)event);
    Bool_t result=ProcessEvent(rawReader, nocheck);
    delete rawReader;
    return result;
#else
    Fatal("AliTRDCalibPadStatus", "this class was compiled without DATE");
    return 0;
#endif

}

//_____________________________________________________________________
Int_t AliTRDCalibPadStatus::ProcessEvent2(AliRawReader *rawReader)
{
  //
  // Event Processing loop - AliTRDRawStreamCosmic
  // 0 time bin problem or zero suppression
  // 1 no input
  // 2 input
  // Raw version number: 
  // [3,31] non zero suppressed
  // 2,4 and [32,63] zero suppressed 
  //
  
  Int_t withInput = 1;

  AliTRDrawFastStream *rawStream = new AliTRDrawFastStream(rawReader);
  rawStream->SetNoErrorWarning();
  rawStream->SetSharedPadReadout(kTRUE);

  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager(kTRUE);
  digitsManager->CreateArrays();
  
  AliTRDfeeParam *feeParam = AliTRDfeeParam::Instance();

  Int_t det    = 0;
  while ((det = rawStream->NextChamber(digitsManager, NULL, NULL)) >= 0) { //idetector
    if (digitsManager->GetIndexes(det)->HasEntry()) {//QA
      //	printf("there is ADC data on this chamber!\n");
      
       AliTRDarrayADC *digits = (AliTRDarrayADC *) digitsManager->GetDigits(det); //mod
      if (digits->HasData()) { //array
	
	AliTRDSignalIndex   *indexes = digitsManager->GetIndexes(det);
	if (indexes->IsAllocated() == kFALSE) {
	  AliError("Indexes do not exist!");
	  break;
	}
	Int_t iRow  = 0;
	Int_t iCol  = 0;
	indexes->ResetCounters();

	while (indexes->NextRCIndex(iRow, iCol)) { //column,row
	  

	  AliTRDdigitsParam *digitParam = (AliTRDdigitsParam *)digitsManager->GetDigitsParam();
	  
	  Int_t mcm          = 0;     // MCM from AliTRDfeeParam
	  Int_t rob          = 0;     // ROB from AliTRDfeeParam
	  Int_t extCol       = 0;     // extended column from  AliTRDfeeParam  
	  mcm = feeParam->GetMCMfromPad(iRow,iCol);
	  rob = feeParam->GetROBfromPad(iRow,iCol);

	  Int_t idetector  = det;                            //  current detector
	  Int_t iRowMax    = rawStream->GetMaxRow();         //  current rowmax

	  Int_t adc        = 20 - (iCol%18) -1;                 //  current adc
	  Int_t col        = 0;                              //  col!=0 ->Shared Pad
	  extCol = feeParam->GetExtendedPadColFromADC(rob,mcm,adc);
	  //printf("  iCol %d  iRow %d  iRowMax %d  rob %d  mcm %d  adc %d  extCol %d\n",iCol,iRow,iRowMax,rob,mcm,adc,extCol);	  
	  
	  // Signal for regular pads
	  Int_t nbtimebin  = digitParam->GetNTimeBins(idetector);  //  number of time bins read from data	  
	  for(Int_t k = 0; k < nbtimebin; k++){
	    Short_t signal = 0;
	    signal = digits->GetData(iRow,iCol,k);
	    
	    if(signal>0) {
	      UpdateHisto2(idetector,iRow,iCol,signal,iRowMax,col,mcm,rob);
	    }
	  }

	  
	  
	  if((adc==3-1 || adc==20-1 || adc==19-1) && (iCol > 1 && iCol <142) /* && fSharedPadsOn*/ ) { //SHARED PADS
	    
	      switch(adc) {
	      case 2:  
		adc = 20;                                       //shared Pad adc 
		mcm = feeParam->GetMCMfromSharedPad(iRow,iCol); //shared Pad mcm 
		col =  1;
		break;
	      case 19:  
		adc = 1;                                        //shared Pad adc  
		mcm = feeParam->GetMCMfromSharedPad(iRow,iCol); //shared Pad mcm  
		col =  2;
		break;
	      case 18: 
		adc =  0;                                       //shared Pad adc  
		mcm = feeParam->GetMCMfromSharedPad(iRow,iCol); //shared Pad mcm 
		col =  3;
		break;
	      }
	      rob = feeParam->GetROBfromSharedPad(iRow,iCol);     //shared Pad rob 
	      
	    
	    extCol = feeParam->GetExtendedPadColFromADC(rob,mcm,adc);     //extended pad col via the shared pad rob,mcm and adc
	    
	    //printf("SHARED PAD ---  iCol %d  iRow %d  rob %d  mcm %d  adc %d  extCol %d  col %d\n",iCol,iRow,rob,mcm,adc,extCol,col);
	    for(Int_t k = 0; k < nbtimebin; k++){
	      Short_t signal = 0;
	      signal = digits->GetDataByAdcCol(iRow,extCol,k);
	      
	      if(signal>0) {
		UpdateHisto2(idetector,iRow,iCol,signal,iRowMax,col,mcm,rob);
	      }
	    }
	  } //shared pads end
	    
	  
	  withInput = 2;
	}//column,row

      }//array
    }//QA
    digitsManager->ClearArrays(det);
  }//idetector
  delete digitsManager;
  delete rawStream;
  return withInput;
}

//_____________________________________________________________________

Int_t AliTRDCalibPadStatus::ProcessEvent3(AliRawReader *rawReader)
{
  //
  // RawReader = AliTRDrawStream (Jochen Klein) 
  //
  // Event Processing loop - AliTRDRawStreamCosmic
  // 0 time bin problem or zero suppression
  // 1 no input
  // 2 input
  // Raw version number: 
  // [3,31] non zero suppressed
  // 2,4 and [32,63] zero suppressed 
  //
  
 

  Int_t withInput = 1;

  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager(kTRUE);
  digitsManager->CreateArrays();

  AliTRDrawStream *rawStream = new AliTRDrawStream(rawReader);
  rawStream->SetDigitsManager(digitsManager);
  //rawStream->SetNoErrorWarning();
  //rawStream->SetSharedPadReadout(kTRUE);
  
  AliTRDfeeParam *feeParam = AliTRDfeeParam::Instance();

  Int_t det    = 0;
  while ((det = rawStream->NextChamber(digitsManager, NULL, NULL)) >= 0) { //idetector
    if (digitsManager->GetIndexes(det)->HasEntry()) {//QA
      //      printf("there is ADC data on this chamber!\n");
      
       AliTRDarrayADC *digits = (AliTRDarrayADC *) digitsManager->GetDigits(det); //mod
      if (digits->HasData()) { //array
	
	AliTRDSignalIndex   *indexes = digitsManager->GetIndexes(det);
	if (indexes->IsAllocated() == kFALSE) {
	  AliError("Indexes do not exist!");
	  break;
	}
	Int_t iRow  = 0;
	Int_t iCol  = 0;
	indexes->ResetCounters();
	
	while (indexes->NextRCIndex(iRow, iCol)) { //column,row
	
	  AliTRDdigitsParam *digitParam = (AliTRDdigitsParam *)digitsManager->GetDigitsParam();
	  
	  Int_t mcm          = 0;     // MCM from AliTRDfeeParam
	  Int_t rob          = 0;     // ROB from AliTRDfeeParam
	  Int_t extCol       = 0;     // extended column from AliTRDfeeParam  
	  mcm = feeParam->GetMCMfromPad(iRow,iCol);
	  rob = feeParam->GetROBfromPad(iRow,iCol);
	  
	  Int_t idetector  = det;                            //  current detector
	  Int_t iRowMax    = 16;                              //  current rowmax
	  if(GetStack(det) == 2) iRowMax = 12;
	  	  
	  Int_t adc        = 20 - (iCol%18) -1;                 //  current adc
	  Int_t col        = 0;                              //  col!=0 ->Shared Pad
	  extCol = feeParam->GetExtendedPadColFromADC(rob,mcm,adc);
	  //printf("  iCol %d  iRow %d  iRowMax %d  rob %d  mcm %d  adc %d  extCol %d\n",iCol,iRow,iRowMax,rob,mcm,adc,extCol);	  
	  
	  // Signal for regular pads
	  Int_t nbtimebin  = digitParam->GetNTimeBins(idetector);  //  number of time bins read from data	  
	  for(Int_t k = 0; k < nbtimebin; k++){
	    Short_t signal = 0;
	    signal = digits->GetData(iRow,iCol,k);

	    if(signal>0) {
	      UpdateHisto2(idetector,iRow,iCol,signal,iRowMax,col,mcm,rob);
	    }
	  }
	  
	  
	  
	  if((adc==3-1 || adc==20-1 || adc==19-1) && (iCol > 1 && iCol <142)  ) { //SHARED PADS
	    
	    switch(adc) {
	    case 2:  
	      adc = 20;                                       //shared Pad adc 
	      mcm = feeParam->GetMCMfromSharedPad(iRow,iCol); //shared Pad mcm 
	      col =  1;
	      break;
	    case 19:  
	      adc = 1;                                        //shared Pad adc  
	      mcm = feeParam->GetMCMfromSharedPad(iRow,iCol); //shared Pad mcm  
	      col =  2;
	      break;
	    case 18: 
	      adc =  0;                                       //shared Pad adc  
	      mcm = feeParam->GetMCMfromSharedPad(iRow,iCol); //shared Pad mcm 
	      col =  3;
	      break;
	    }
	    rob = feeParam->GetROBfromSharedPad(iRow,iCol);     //shared Pad rob 
	    
	    
	    extCol = feeParam->GetExtendedPadColFromADC(rob,mcm,adc);     //extended pad col via the shared pad rob,mcm and adc
	    
	    //printf("SHARED PAD ---  iCol %d  iRow %d  rob %d  mcm %d  adc %d  extCol %d  col %d\n",iCol,iRow,rob,mcm,adc,extCol,col);
	    for(Int_t k = 0; k < nbtimebin; k++){
	      Short_t signal = 0;
	      signal = digits->GetDataByAdcCol(iRow,extCol,k);
	      
	      if(signal>0) {
		UpdateHisto2(idetector,iRow,iCol,signal,iRowMax,col,mcm,rob);
	      }
	    }
	  } //shared pads end
	  
	  
	  withInput = 2;
	}//column,row

      }//array
    }//QA
    digitsManager->ClearArrays(det);
  }//idetector
  delete digitsManager;
  delete rawStream;
  return withInput;
  
}


//_____________________________________________________________________
Bool_t AliTRDCalibPadStatus::TestEventHisto(Int_t nevent, Int_t sm, Int_t ch) /*FOLD00*/
{
  //
  //  Test event loop
  // fill one oroc and one iroc with random gaus
  //

  gRandom->SetSeed(0);

    for (Int_t ism=sm; ism<sm+1; ism++){
       	for (Int_t ich=ch; ich < ch+1; ich++){
	    for (Int_t ipl=0; ipl < 6; ipl++){
	      for(Int_t irow = 0; irow < fGeo->GetRowMax(ipl,ich,ism); irow++){
		for(Int_t icol = 0; icol < fGeo->GetColMax(ipl); icol++){
		  for (Int_t iTimeBin=0; iTimeBin<(30*nevent); iTimeBin++){
		    Int_t signal=TMath::Nint(gRandom->Gaus(10,1.5));
		    if ( signal>0 )UpdateHisto((ipl+ich*6+ism*6*5),irow,icol,signal,fGeo->GetRowMax(ipl,ich,ism),0,0);
		  }
		}
	      }
	    }
	}
    }
    return kTRUE;
}

//_____________________________________________________________________
TH2F* AliTRDCalibPadStatus::GetHisto(Int_t det, TObjArray *arr, /*FOLD00*/
				  Int_t nbinsY, Float_t ymin, Float_t ymax,
				  const Char_t *type, Bool_t force)
{
    //
    // return pointer to histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(det) )
	return (TH2F*)arr->UncheckedAt(det);

    // if we are forced and histogram doesn't yes exist create it
    Char_t name[255], title[255];

    sprintf(name,"hCalib%s%.3d",type,det);
    sprintf(title,"%s calibration histogram detector %.2d;ADC channel;Channel (pad)",type,det);

   
    Int_t nbchannels = fGeo->GetRowMax(GetLayer(det),GetStack(det),GetSector(det))*fGeo->GetColMax(GetLayer(det));
    
    // we will add 3*8*rowMax channels at the end for the double counted
    nbchannels += 3*8*(fGeo->GetRowMax(GetLayer(det),GetStack(det),GetSector(det)));


    // new histogram with calib information. One value for each pad!
    TH2F* hist = new TH2F(name,title,
			  nbinsY, ymin, ymax,
			  nbchannels,0,nbchannels
			  );
    hist->SetDirectory(0);
    arr->AddAt(hist,det);
    return hist;
}

//_____________________________________________________________________
TH2F* AliTRDCalibPadStatus::GetHisto(Int_t det, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fHistoArray;
    return GetHisto(det, arr, fAdcMax-fAdcMin, fAdcMin-0.5, fAdcMax-0.5, "Pedestal", force);
}

//_____________________________________________________________________
AliTRDCalROC* AliTRDCalibPadStatus::GetCalRoc(Int_t det, TObjArray* arr, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to ROC Calibration
    // if force is true create a new AliTRDCalROC if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(det) )
	return (AliTRDCalROC*)arr->UncheckedAt(det);

    // if we are forced and histogram doesn't yes exist create it

    // new AliTRDCalROC. One value for each pad!
    AliTRDCalROC *croc = new AliTRDCalROC(GetLayer(det),GetStack(det));
    arr->AddAt(croc,det);
    return croc;
}
//_____________________________________________________________________
AliTRDCalROC* AliTRDCalibPadStatus::GetCalRocMean(Int_t det, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayMean;
    return GetCalRoc(det, arr, force);
}

//_____________________________________________________________________
AliTRDCalROC* AliTRDCalibPadStatus::GetCalRocRMS(Int_t det, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayRMS;
    return GetCalRoc(det, arr, force);
}
//_____________________________________________________________________
AliTRDCalROC* AliTRDCalibPadStatus::GetCalRocMeand(Int_t det, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayMeand;
    return GetCalRoc(det, arr, force);
}

//_____________________________________________________________________
AliTRDCalROC* AliTRDCalibPadStatus::GetCalRocRMSd(Int_t det, Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayRMSd;
    return GetCalRoc(det, arr, force);
}

//_____________________________________________________________________
void AliTRDCalibPadStatus::AnalyseHisto() /*FOLD00*/
{
    //
    //  Calculate calibration constants
    //

    Int_t nbinsAdc = fAdcMax-fAdcMin;

    TVectorD param(4);
    TMatrixD dummy(3,3);

    Float_t *arrayHP=0;


    for (Int_t idet=0; idet<540; idet++){
	TH2F *hP = GetHisto(idet);
        if ( !hP ) {
	  continue;
	}

        //printf("Entries for %d\n",idet);

	AliTRDCalROC *rocMean     = GetCalRocMean(idet,kTRUE);
	AliTRDCalROC *rocRMS      = GetCalRocRMS(idet,kTRUE);

	arrayHP = hP->GetArray();
        Int_t nChannels = rocMean->GetNchannels();

	for (Int_t iChannel=0; iChannel<nChannels; iChannel++){
            Int_t offset = (nbinsAdc+2)*(iChannel+1)+1;
	    Double_t ret = AliMathBase::FitGaus(arrayHP+offset,nbinsAdc,fAdcMin-0.5,fAdcMax-0.5,&param,&dummy);
	    // if the fitting failed set noise and pedestal to 0
	    if ((ret==-4) || (ret==-1) || (ret==-2)) {
		param[1]=0.0;
		param[2]=0.0;
	    }
	    if((param[1]/10.0) > 65534.0) param[1] = 0.0;
	    if((param[2]/10.0) > 65534.0) param[2] = 0.0;
	    rocMean->SetValue(iChannel,param[1]/10.0);
            rocRMS->SetValue(iChannel,param[2]/10.0);
	}

	// here we analyse doubled read channels
	
	AliTRDCalROC *rocMeand     = GetCalRocMeand(idet,kTRUE);
	AliTRDCalROC *rocRMSd      = GetCalRocRMSd(idet,kTRUE);

	Int_t nrows = rocMeand->GetNrows();
	Int_t shift = 144*nrows;
	Int_t total = shift+3*8*nrows; 

      	for (Int_t iChannel=shift; iChannel<total; iChannel++){
            Int_t offset = (nbinsAdc+2)*(iChannel+1)+1;
	    Double_t ret = AliMathBase::FitGaus(arrayHP+offset,nbinsAdc,fAdcMin-0.5,fAdcMax-0.5,&param,&dummy);
            // if the fitting failed set noise and pedestal to 0
	    if ((ret==-4) || (ret==-1) || (ret==-2)) {
		param[1]=0.0;
		param[2]=0.0;
	    }
	    if((param[1]/10.0) > 65534.0) param[1] = 0.0;
	    if((param[2]/10.0) > 65534.0) param[2] = 0.0;
	    
	    // here we have to recalculate backward
	    Int_t nb   = iChannel-shift;
	    Int_t row  = nb%nrows;
	    Int_t j    = (Int_t)(nb/nrows);
	    Int_t imcm = j%8;
	    Int_t icol = (Int_t)(j/8);
	    
	    Int_t finalcol = 18*imcm;
	    if(icol > 0) icol += 17;
	    else icol = -1;
	    finalcol += icol;

	    Int_t channel = row+finalcol*nrows;

	    //printf("iChannel %d, nrows %d, finalcol %d, row %d, channel %d\n",iChannel,nrows,finalcol,row,channel);
	    if((finalcol < 0) || (finalcol >= 144)) continue;

	    rocMeand->SetValue(channel,param[1]/10.0);
            rocRMSd->SetValue(channel,param[2]/10.0);
	}

    }
      
}

//_______________________________________________________________________________________
AliTRDCalPadStatus* AliTRDCalibPadStatus::CreateCalPadStatus()
{
  //
  // Create Pad Status out of Mean and RMS values
  // The chamber without data are masked, this is the corrected in the preprocessor
  //

  AliTRDCalPadStatus* obj = new AliTRDCalPadStatus("padstatus", "padstatus");
  
  for (Int_t idet=0; idet<540; ++idet)
    {
      AliTRDCalSingleChamberStatus *calROC = obj->GetCalROC(idet);


      if ( !GetCalRocMean(idet)) {
	for(Int_t k = 0; k < calROC->GetNchannels(); k++){
	  calROC->SetStatus(k,AliTRDCalPadStatus::kNotConnected);
	}
	continue;
      }
      

      //Take the stuff
      AliTRDCalROC *calRocMean    = new AliTRDCalROC(*( (AliTRDCalROC *) GetCalRocMean(idet)));
      AliTRDCalROC *calRocRMS     = new AliTRDCalROC(*( (AliTRDCalROC *) GetCalRocRMS(idet)));

      //Take the stuff second chance
      AliTRDCalROC *calRocMeand    = new AliTRDCalROC(*( (AliTRDCalROC *) GetCalRocMeand(idet)));
      AliTRDCalROC *calRocRMSd     = new AliTRDCalROC(*( (AliTRDCalROC *) GetCalRocRMSd(idet)));

      calRocRMS->Unfold();
      calRocRMSd->Unfold();

     
      //Range
      Int_t row      = calROC->GetNrows();
      Int_t col      = calROC->GetNcols();
      
      Double_t rmsmean       = calRocMean->GetRMS()*10.0;
      Double_t meanmean      = calRocMean->GetMean()*10.0;
      Double_t meansquares   = calRocRMS->GetMean();

      
      for(Int_t irow = 0; irow < row; irow++){
	
	// for bridged pads
	Float_t meanprevious = 0.0;
	Float_t rmsprevious  = 0.0; 
	Float_t mean         = 0.0;
	Float_t rms          = 0.0;
	
	for(Int_t icol = 0; icol < col; icol++){
	  
	  mean     = calRocMean->GetValue(icol,irow)*10.0;
	  rms      = calRocRMS->GetValue(icol,irow);

	  if(icol > 0) {
	    meanprevious     = calRocMean->GetValue((icol -1),irow)*10.0;
	    rmsprevious      = calRocRMS->GetValue((icol - 1),irow);
	  }
	  
	  Bool_t pb = kFALSE;
	  // masked if two noisy
	  if((rms <= 0.0001) || (TMath::Abs(mean-meanmean)>(5*rmsmean)) || (TMath::Abs(rms)>(5.0*TMath::Abs(meansquares)))) {
	    
	    pb = kTRUE;
	    // look at second chance
	    Float_t meand     = calRocMeand->GetValue(icol,irow)*10.0;
	    Float_t rmsd      = calRocRMSd->GetValue(icol,irow);
	    
	    if((rmsd <= 0.0001) || (TMath::Abs(meand-meanmean)>(5*rmsmean)) || (TMath::Abs(rmsd)>(5.0*TMath::Abs(meansquares)))) {
	      if((rmsd <= 0.0001) && (rms <= 0.0001)) {
		calROC->SetStatus(icol,irow,AliTRDCalPadStatus::kNotConnected);
	      }
	      else {
		calROC->SetStatus(icol, irow, AliTRDCalPadStatus::kMasked);
	      }
	    }
	    else {
	      calROC->SetStatus(icol, irow, AliTRDCalPadStatus::kReadSecond);
	    }
	  }


	  // bridge if previous pad found something
	  if(!pb) {
	    if((meanprevious == mean) && (rmsprevious == rms) && (mean > 0.0001)) {
	      //printf("mean previous %f, mean %f, rms %f, rmsprevious %f, col %d\n",meanprevious,mean,rms,rmsprevious,icol);
	      calROC->SetStatus(icol -1 ,irow, AliTRDCalPadStatus::kPadBridgedRight);
	      calROC->SetStatus(icol ,irow, AliTRDCalPadStatus::kPadBridgedLeft);
	    }	    
	  }

	}
      }

      delete calRocMean;
      delete calRocRMS;
      delete calRocMeand;
      delete calRocRMSd;


    }
  
  return obj;
  
}
//_______________________________________________________________________________________
AliTRDCalPad* AliTRDCalibPadStatus::CreateCalPad()
{
  //
  // Create Pad Noise out of RMS values
  //

  AliTRDCalPad* obj = new AliTRDCalPad("PadNoise", "PadNoise");
  
  
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)  {
    
    AliTRDCalROC *calROC22 = obj->GetCalROC(det);

    AliTRDCalROC *calRocRMS     = ((AliTRDCalROC *)GetCalRocRMS(det,kTRUE));
   
    for(Int_t k = 0; k < calROC22->GetNchannels(); k++){
      calROC22->SetValue(k,calRocRMS->GetValue(k));
    }

  }
  
  return obj;
  
}

//_______________________________________________________________________________________
AliTRDCalDet* AliTRDCalibPadStatus::CreateCalDet() const
{
  //
  // Create Det Noise correction factor
  //

  AliTRDCalDet* obj = new AliTRDCalDet("DetNoise", "DetNoise (correction factor)");

  for(Int_t l = 0; l < 540; l++){
    obj->SetValue(l,10.0);
  }
  
  return obj;
  
}

//_____________________________________________________________________
void AliTRDCalibPadStatus::DumpToFile(const Char_t *filename, const Char_t *dir, Bool_t append) /*FOLD00*/
{
    //
    //  Write class to file
    //

    TString sDir(dir);
    TString option;

    if ( append )
	option = "update";
    else
        option = "recreate";

    TDirectory *backup = gDirectory;
    TFile f(filename,option.Data());
    f.cd();
    if ( !sDir.IsNull() ){
	f.mkdir(sDir.Data());
	f.cd(sDir);
    }
    this->Write();
    f.Close();

    if ( backup ) backup->cd();
}

//_____________________________________________________________________
void AliTRDCalibPadStatus::SetCalRocMean(AliTRDCalROC *mean, Int_t det) /*FOLD00*/
{
    //
    //  Put the AliTRDCalROC in the array fCalRocArrayMean
    //


  AliTRDCalROC *rocMean = GetCalRocMean(det,kTRUE);
  
  Int_t nChannels = rocMean->GetNchannels();
  
  for (Int_t iChannel=0; iChannel<nChannels; iChannel++){
    
    rocMean->SetValue(iChannel,mean->GetValue(iChannel));
    
  }
  
}

//_____________________________________________________________________
void AliTRDCalibPadStatus::SetCalRocRMS(AliTRDCalROC *rms, Int_t det) /*FOLD00*/
{
    //
    //  Put the AliTRDCalROC in the array fCalRocArrayRMS
    //


  AliTRDCalROC *rocRms = GetCalRocRMS(det,kTRUE);
  
  Int_t nChannels = rocRms->GetNchannels();
  
  for (Int_t iChannel=0; iChannel<nChannels; iChannel++){
    
    rocRms->SetValue(iChannel,rms->GetValue(iChannel));
    
  }
  
}
//_____________________________________________________________________
void AliTRDCalibPadStatus::SetCalRocMeand(AliTRDCalROC *mean, Int_t det) /*FOLD00*/
{
    //
    //  Put the AliTRDCalROC in the array fCalRocArrayMean
    //


  AliTRDCalROC *rocMean = GetCalRocMeand(det,kTRUE);
  
  Int_t nChannels = rocMean->GetNchannels();
  
  for (Int_t iChannel=0; iChannel<nChannels; iChannel++){
    
    rocMean->SetValue(iChannel,mean->GetValue(iChannel));
    
  }
  
}

//_____________________________________________________________________
void AliTRDCalibPadStatus::SetCalRocRMSd(AliTRDCalROC *rms, Int_t det) /*FOLD00*/
{
    //
    //  Put the AliTRDCalROC in the array fCalRocArrayRMS
    //


  AliTRDCalROC *rocRms = GetCalRocRMSd(det,kTRUE);
  
  Int_t nChannels = rocRms->GetNchannels();
  
  for (Int_t iChannel=0; iChannel<nChannels; iChannel++){
    
    rocRms->SetValue(iChannel,rms->GetValue(iChannel));
    
  }
  
}
//_____________________________________________________________________________
Int_t AliTRDCalibPadStatus::GetLayer(Int_t d) const
{
  //
  // Reconstruct the layer number from the detector number
  //

  return ((Int_t) (d % 6));

}

//_____________________________________________________________________________
Int_t AliTRDCalibPadStatus::GetStack(Int_t d) const
{
  //
  // Reconstruct the chamber number from the detector number
  //

  return ((Int_t) (d % 30) / 6);

}

//_____________________________________________________________________________
Int_t AliTRDCalibPadStatus::GetSector(Int_t d) const
{
  //
  // Reconstruct the sector number from the detector number
  //

  return ((Int_t) (d / 30));

}
