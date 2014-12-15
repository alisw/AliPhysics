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

/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
//                  Implementation of the TPC pulser calibration                       //
//                                                                                     //
//   Origin: Jens Wiechula, Marian Ivanov   J.Wiechula@gsi.de, Marian.Ivanov@cern.ch   //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *                      Class Description                                  *
 ***************************************************************************

 The AliTPCCalibPulser class is used to get calibration data concerning the FEE using
 runs performed with the calibration pulser.

 The information retrieved is
 - Time0 differences
 - Signal width differences
 - Amplification variations

 the seen differences arise from the manufacturing tolerances of the PASAs and are very small within
 one chip and somewhat large between different chips.

 Histograms:
   For each ROC three TH2S histos 'Reference Histograms'  (ROC channel vs. [Time0, signal width, Q sum]) is created when
   it is filled for the first time (GetHisto[T0,RMS,Q](ROC,kTRUE)). The histos are stored in the
   TObjArrays fHistoT0Array, fHistoRMSArray and fHistoQArray.


 Working principle:
 ------------------
 Raw calibration pulser data is processed by calling one of the ProcessEvent(...) functions
 (see below). These in the end call the Update(...) function.

 - the Update(...) function:
   In this function the array fPadSignal is filled with the adc signals between the specified range
   fFirstTimeBin and fLastTimeBin for the current pad.
   before going to the next pad the ProcessPad() function is called, which analyses the data for one pad
   stored in fPadSignal.

   - the ProcessPad() function:
     Find Pedestal and Noise information
     - use database information which has to be set by calling
       SetPedestalDatabase(AliTPCCalPad *pedestalTPC, AliTPCCalPad *padNoiseTPC)
     - if no information from the pedestal data base
       is available the informaion is calculated on the fly ( see FindPedestal() function )

     Find the Pulser signal information
     - calculate  mean = T0, RMS = signal width and Q sum in a range of -2+7 timebins around Q max
       the Q sum is scaled by pad area
       (see FindPulserSignal(...) function)

     Fill a temprary array for the T0 information (GetPadTimesEvent(fCurrentSector,kTRUE)) (why see below)
     Fill the Q sum and RMS values in the histograms (GetHisto[RMS,Q](ROC,kTRUE)),

 At the end of each event the EndEvent() function is called

 - the EndEvent() function:
   calculate the mean T0 for each ROC and fill the Time0 histogram with Time0-<Time0 for ROC>
   This is done to overcome syncronisation problems between the trigger and the fec clock.

 After accumulating the desired statistics the Analyse() function has to be called.
 - the Analyse() function
   Whithin this function the mean values of T0, RMS, Q are calculated for each pad, using
   the AliMathBase::GetCOG(...) function, and the calibration
   storage classes (AliTPCCalROC) are filled for each ROC.
   The calibration information is stored in the TObjArrays fCalRocArrayT0, fCalRocArrayRMS and
   fCalRocArrayQ;



 User interface for filling data:
 --------------------------------

 To Fill information one of the following functions can be used:

 Bool_t ProcessEvent(eventHeaderStruct *event);
   - process Date event
   - use AliTPCRawReaderDate and call ProcessEvent(AliRawReader *rawReader)

 Bool_t ProcessEvent(AliRawReader *rawReader);
   - process AliRawReader event
   - use AliTPCRawStreamV3 to loop over data and call ProcessEvent(AliTPCRawStreamV3 *rawStream)

 Bool_t ProcessEvent(AliTPCRawStreamV3 *rawStream);
   - process event from AliTPCRawStreamV3
   - call Update function for signal filling

 Int_t Update(const Int_t isector, const Int_t iRow, const Int_t
              iPad,  const Int_t iTimeBin, const Float_t signal);
   - directly  fill signal information (sector, row, pad, time bin, pad)
     to the reference histograms

 It is also possible to merge two independently taken calibrations using the function

 void Merge(AliTPCCalibPulser *sig)
   - copy histograms in 'sig' if the do not exist in this instance
   - Add histograms in 'sig' to the histograms in this instance if the allready exist
   - After merging call Analyse again!



 -- example: filling data using root raw data:
 void fillSignal(Char_t *filename)
 {
    rawReader = new AliRawReaderRoot(fileName);
    if ( !rawReader ) return;
    AliTPCCalibPulser *calib = new AliTPCCalibPulser;
    while (rawReader->NextEvent()){
      calib->ProcessEvent(rawReader);
    }
    calib->Analyse();
    calib->DumpToFile("SignalData.root");
    delete rawReader;
    delete calib;
 }


 What kind of information is stored and how to retrieve them:
 ------------------------------------------------------------

 - Accessing the 'Reference Histograms' (Time0, signal width and Q sum information pad by pad):

   TH2F *GetHistoT0(Int_t sector);
   TH2F *GetHistoRMS(Int_t sector);
   TH2F *GetHistoQ(Int_t sector);

 - Accessing the calibration storage objects:

   AliTPCCalROC *GetCalRocT0(Int_t sector);   // for the Time0 values
   AliTPCCalROC *GetCalRocRMS(Int_t sector);  // for the signal width values
   AliTPCCalROC *GetCalRocQ(Int_t sector);    // for the Q sum values

   example for visualisation:
   if the file "SignalData.root" was created using the above example one could do the following:

   TFile fileSignal("SignalData.root")
   AliTPCCalibPulser *sig = (AliTPCCalibPulser*)fileSignal->Get("AliTPCCalibPulser");
   sig->GetCalRocT0(0)->Draw("colz");
   sig->GetCalRocRMS(0)->Draw("colz");

   or use the AliTPCCalPad functionality:
   AliTPCCalPad padT0(ped->GetCalPadT0());
   AliTPCCalPad padSigWidth(ped->GetCalPadRMS());
   padT0->MakeHisto2D()->Draw("colz");       //Draw A-Side Time0 Information
   padSigWidth->MakeHisto2D()->Draw("colz"); //Draw A-Side signal width Information
*/


//Root includes
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2S.h>
#include <TString.h>
#include <TVectorF.h>
#include <TMath.h>
#include <TMap.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>

//AliRoot includes
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCROC.h"
#include "AliTPCParam.h"
#include "AliTPCCalibPulser.h"
#include "AliTPCcalibDB.h"
#include "AliMathBase.h"
#include "AliLog.h"
#include "TTreeStream.h"

//date
#include "event.h"




ClassImp(AliTPCCalibPulser)

AliTPCCalibPulser::AliTPCCalibPulser() :
AliTPCCalibRawBase(),
fNbinsT0(200),
fXminT0(-2),
fXmaxT0(2),
fNbinsQ(200),
fXminQ(10),
fXmaxQ(40),
fNbinsRMS(100),
fXminRMS(0.1),
fXmaxRMS(5.1),
fPeakIntMinus(2),
fPeakIntPlus(2),
fIsZeroSuppressed(kFALSE),
fLastSector(-1),
fParam(new AliTPCParam),
fPedestalTPC(0x0),
fPadNoiseTPC(0x0),
fOutliers(0x0),
fPedestalROC(0x0),
fPadNoiseROC(0x0),
fCalRocArrayT0(72),
fCalRocArrayQ(72),
fCalRocArrayRMS(72),
fCalRocArrayOutliers(72),
fHistoQArray(72),
fHistoT0Array(72),
fHistoRMSArray(72),
fHMeanTimeSector(0x0),
fVMeanTimeSector(72),
fPadTimesArrayEvent(72),
fPadQArrayEvent(72),
fPadRMSArrayEvent(72),
fPadPedestalArrayEvent(72),
fCurrentChannel(-1),
fCurrentSector(-1),
fCurrentRow(-1),
fCurrentPad(-1),
fMaxPadSignal(-1),
fMaxTimeBin(-1),
fPadSignal(1024),
fPadPedestal(0),
fPadNoise(0),
fVTime0Offset(72),
fVTime0OffsetCounter(72)
{
  //
  // AliTPCSignal default constructor
  //
  SetNameTitle("AliTPCCalibPulser","AliTPCCalibPulser");
  fFirstTimeBin=60;
  fLastTimeBin=900;
  fParam->Update();
}
//_____________________________________________________________________
AliTPCCalibPulser::AliTPCCalibPulser(const AliTPCCalibPulser &sig) :
AliTPCCalibRawBase(sig),
fNbinsT0(sig.fNbinsT0),
fXminT0(sig.fXminT0),
fXmaxT0(sig.fXmaxT0),
fNbinsQ(sig.fNbinsQ),
fXminQ(sig.fXminQ),
fXmaxQ(sig.fXmaxQ),
fNbinsRMS(sig.fNbinsRMS),
fXminRMS(sig.fXminRMS),
fXmaxRMS(sig.fXmaxRMS),
fPeakIntMinus(sig.fPeakIntMinus),
fPeakIntPlus(sig.fPeakIntPlus),
fIsZeroSuppressed(sig.fIsZeroSuppressed),
fLastSector(-1),
fParam(new AliTPCParam),
fPedestalTPC(0x0),
fPadNoiseTPC(0x0),
fOutliers(0x0),
fPedestalROC(0x0),
fPadNoiseROC(0x0),
fCalRocArrayT0(72),
fCalRocArrayQ(72),
fCalRocArrayRMS(72),
fCalRocArrayOutliers(72),
fHistoQArray(72),
fHistoT0Array(72),
fHistoRMSArray(72),
fHMeanTimeSector(0x0),
fVMeanTimeSector(72),
fPadTimesArrayEvent(72),
fPadQArrayEvent(72),
fPadRMSArrayEvent(72),
fPadPedestalArrayEvent(72),
fCurrentChannel(-1),
fCurrentSector(-1),
fCurrentRow(-1),
fCurrentPad(-1),
fMaxPadSignal(-1),
fMaxTimeBin(-1),
fPadSignal(1024),
fPadPedestal(0),
fPadNoise(0),
fVTime0Offset(72),
fVTime0OffsetCounter(72)
{
  //
  // AliTPCSignal default constructor
  //
  
  for (Int_t iSec = 0; iSec < 72; ++iSec){
    const AliTPCCalROC *calQ   = (AliTPCCalROC*)sig.fCalRocArrayQ.UncheckedAt(iSec);
    const AliTPCCalROC *calT0  = (AliTPCCalROC*)sig.fCalRocArrayT0.UncheckedAt(iSec);
    const AliTPCCalROC *calRMS = (AliTPCCalROC*)sig.fCalRocArrayRMS.UncheckedAt(iSec);
    const AliTPCCalROC *calOut = (AliTPCCalROC*)sig.fCalRocArrayOutliers.UncheckedAt(iSec);
    
    const TH2S *hQ   = (TH2S*)sig.fHistoQArray.UncheckedAt(iSec);
    const TH2S *hT0  = (TH2S*)sig.fHistoT0Array.UncheckedAt(iSec);
    const TH2S *hRMS = (TH2S*)sig.fHistoRMSArray.UncheckedAt(iSec);
    
    if ( calQ   != 0x0 ) fCalRocArrayQ.AddAt(new AliTPCCalROC(*calQ), iSec);
    if ( calT0  != 0x0 ) fCalRocArrayT0.AddAt(new AliTPCCalROC(*calT0), iSec);
    if ( calRMS != 0x0 ) fCalRocArrayRMS.AddAt(new AliTPCCalROC(*calRMS), iSec);
    if ( calOut != 0x0 ) fCalRocArrayOutliers.AddAt(new AliTPCCalROC(*calOut), iSec);
    
    if ( hQ != 0x0 ){
      TH2S *hNew = new TH2S(*hQ);
      hNew->SetDirectory(0);
      fHistoQArray.AddAt(hNew,iSec);
    }
    if ( hT0 != 0x0 ){
      TH2S *hNew = new TH2S(*hT0);
      hNew->SetDirectory(0);
      fHistoQArray.AddAt(hNew,iSec);
    }
    if ( hRMS != 0x0 ){
      TH2S *hNew = new TH2S(*hRMS);
      hNew->SetDirectory(0);
      fHistoQArray.AddAt(hNew,iSec);
    }
    fVMeanTimeSector[iSec]=sig.fVMeanTimeSector[iSec];
  }
  
  if ( sig.fHMeanTimeSector ) fHMeanTimeSector=(TH2F*)sig.fHMeanTimeSector->Clone();
  fParam->Update();
}
AliTPCCalibPulser::AliTPCCalibPulser(const TMap *config) :
AliTPCCalibRawBase(),
fNbinsT0(200),
fXminT0(-2),
fXmaxT0(2),
fNbinsQ(200),
fXminQ(10),
fXmaxQ(40),
fNbinsRMS(100),
fXminRMS(0.1),
fXmaxRMS(5.1),
fPeakIntMinus(2),
fPeakIntPlus(2),
fIsZeroSuppressed(kFALSE),
fLastSector(-1),
fParam(new  AliTPCParam),
fPedestalTPC(0x0),
fPadNoiseTPC(0x0),
fOutliers(0x0),
fPedestalROC(0x0),
fPadNoiseROC(0x0),
fCalRocArrayT0(72),
fCalRocArrayQ(72),
fCalRocArrayRMS(72),
fCalRocArrayOutliers(72),
fHistoQArray(72),
fHistoT0Array(72),
fHistoRMSArray(72),
fHMeanTimeSector(0x0),
fVMeanTimeSector(72),
fPadTimesArrayEvent(72),
fPadQArrayEvent(72),
fPadRMSArrayEvent(72),
fPadPedestalArrayEvent(72),
fCurrentChannel(-1),
fCurrentSector(-1),
fCurrentRow(-1),
fCurrentPad(-1),
fMaxPadSignal(-1),
fMaxTimeBin(-1),
fPadSignal(1024),
fPadPedestal(0),
fPadNoise(0),
fVTime0Offset(72),
fVTime0OffsetCounter(72)
{
  //
  // This constructor uses a TMap for setting some parametes
  //
  SetNameTitle("AliTPCCalibPulser","AliTPCCalibPulser");
  fFirstTimeBin=60;
  fLastTimeBin=900;
  if (config->GetValue("FirstTimeBin")) fFirstTimeBin = ((TObjString*)config->GetValue("FirstTimeBin"))->GetString().Atoi();
  if (config->GetValue("LastTimeBin")) fLastTimeBin = ((TObjString*)config->GetValue("LastTimeBin"))->GetString().Atoi();
  if (config->GetValue("NbinsT0")) fNbinsT0 = ((TObjString*)config->GetValue("NbinsT0"))->GetString().Atoi();
  if (config->GetValue("XminT0")) fXminT0 = ((TObjString*)config->GetValue("XminT0"))->GetString().Atof();
  if (config->GetValue("XmaxT0")) fXmaxT0 = ((TObjString*)config->GetValue("XmaxT0"))->GetString().Atof();
  if (config->GetValue("NbinsQ")) fNbinsQ = ((TObjString*)config->GetValue("NbinsQ"))->GetString().Atoi();
  if (config->GetValue("XminQ")) fXminQ = ((TObjString*)config->GetValue("XminQ"))->GetString().Atof();
  if (config->GetValue("XmaxQ")) fXmaxQ = ((TObjString*)config->GetValue("XmaxQ"))->GetString().Atof();
  if (config->GetValue("NbinsRMS")) fNbinsRMS = ((TObjString*)config->GetValue("NbinsRMS"))->GetString().Atoi();
  if (config->GetValue("XminRMS")) fXminRMS = ((TObjString*)config->GetValue("XminRMS"))->GetString().Atof();
  if (config->GetValue("XmaxRMS")) fXmaxRMS = ((TObjString*)config->GetValue("XmaxRMS"))->GetString().Atof();
  if (config->GetValue("PeakIntMinus")) fPeakIntMinus = (Int_t)((TObjString*)config->GetValue("PeakIntMinus"))->GetString().Atof();
  if (config->GetValue("PeakIntPlus")) fPeakIntPlus = (Int_t)((TObjString*)config->GetValue("PeakIntPlus"))->GetString().Atof();
  if (config->GetValue("IsZeroSuppressed")) fIsZeroSuppressed = (Bool_t)((TObjString*)config->GetValue("IsZeroSuppressed"))->GetString().Atoi();
  
  fParam->Update();
}
//_____________________________________________________________________
AliTPCCalibPulser& AliTPCCalibPulser::operator = (const  AliTPCCalibPulser &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCCalibPulser(source);
  
  return *this;
}
//_____________________________________________________________________
AliTPCCalibPulser::~AliTPCCalibPulser()
{
  //
  // destructor
  //
  
  Reset();
  delete fParam;
}
void AliTPCCalibPulser::Reset()
{
  //
  // Delete all information: Arrays, Histograms, CalRoc objects
  //
  fCalRocArrayT0.Delete();
  fCalRocArrayQ.Delete();
  fCalRocArrayRMS.Delete();
  fCalRocArrayOutliers.Delete();
  
  fHistoQArray.Delete();
  fHistoT0Array.Delete();
  fHistoRMSArray.Delete();
  
  fPadTimesArrayEvent.Delete();
  fPadQArrayEvent.Delete();
  fPadRMSArrayEvent.Delete();
  fPadPedestalArrayEvent.Delete();
  
  if (fHMeanTimeSector) delete fHMeanTimeSector;
}
//_____________________________________________________________________
Int_t AliTPCCalibPulser::Update(const Int_t icsector,
                                const Int_t icRow,
                                const Int_t icPad,
                                const Int_t icTimeBin,
                                const Float_t csignal)
{
    //
    // Signal filling methode on the fly pedestal and time offset correction if necessary.
    // no extra analysis necessary. Assumes knowledge of the signal shape!
    // assumes that it is looped over consecutive time bins of one pad
    //
  
  if (icRow<0) return 0;
  if (icPad<0) return 0;
  if (icTimeBin<0) return 0;
  if ( (icTimeBin>fLastTimeBin) || (icTimeBin<fFirstTimeBin)   ) return 0;
  
  if ( icRow<0 || icPad<0 ){
    AliWarning("Wrong Pad or Row number, skipping!");
    return 0;
  }
  
  Int_t iChannel  = fROC->GetRowIndexes(icsector)[icRow]+icPad; //  global pad position in sector
  
    //init first pad and sector in this event
  if ( fCurrentChannel == -1 ) {
    fCurrentChannel = iChannel;
    fCurrentSector  = icsector;
    fCurrentRow     = icRow;
    fCurrentPad     = icPad;
  }
  
    //process last pad if we change to a new one
  if ( iChannel != fCurrentChannel ){
    ProcessPad();
    fLastSector=fCurrentSector;
    fCurrentChannel = iChannel;
    fCurrentSector  = icsector;
    fCurrentRow     = icRow;
    fCurrentPad     = icPad;
  }
  
    //fill signals for current pad
  fPadSignal[icTimeBin]=csignal;
  if ( csignal > fMaxPadSignal ){
    fMaxPadSignal = csignal;
    fMaxTimeBin   = icTimeBin;
  }
  return 0;
}
//_____________________________________________________________________
void AliTPCCalibPulser::FindPedestal(Float_t part)
{
    //
    // find pedestal and noise for the current pad. Use either database or
    // truncated mean with part*100%
    //
  Bool_t noPedestal = kTRUE;;
  if (fPedestalTPC&&fPadNoiseTPC){
        //use pedestal database
        //only load new pedestals if the sector has changed
    if ( fCurrentSector!=fLastSector ){
      fPedestalROC = fPedestalTPC->GetCalROC(fCurrentSector);
      fPadNoiseROC = fPadNoiseTPC->GetCalROC(fCurrentSector);
    }
    
    if ( fPedestalROC&&fPadNoiseROC ){
      fPadPedestal = fPedestalROC->GetValue(fCurrentChannel);
      fPadNoise    = fPadNoiseROC->GetValue(fCurrentChannel);
      noPedestal   = kFALSE;
    }
    
  }
  
    //if we are not running with pedestal database, or for the current sector there is no information
    //available, calculate the pedestal and noise on the fly
  if ( noPedestal ) {
    const Int_t kPedMax = 100;  //maximum pedestal value
    Float_t  max    =  0;
    Int_t    median =  -1;
    Int_t    count0 =  0;
    Int_t    count1 =  0;
  //
    Float_t padSignal=0;
        //
    UShort_t histo[kPedMax];
    memset(histo,0,kPedMax*sizeof(UShort_t));
    
    for (Int_t i=fFirstTimeBin; i<=fLastTimeBin; ++i){
      padSignal = fPadSignal.GetMatrixArray()[i];
      if (padSignal<=0) continue;
      if (padSignal>max && i>10) {
        max = padSignal;
      }
      if (padSignal>kPedMax-1) continue;
      histo[Int_t(padSignal+0.5)]++;
      count0++;
    }
      //
    for (Int_t i=1; i<kPedMax; ++i){
      if (count1<count0*0.5) median=i;
      count1+=histo[i];
    }
  // truncated mean
  //
        // what if by chance histo[median] == 0 ?!?
    Float_t count=histo[median] ,mean=histo[median]*median,  rms=histo[median]*median*median ;
  //
    for (Int_t idelta=1; idelta<10; ++idelta){
      if (median-idelta<=0) continue;
      if (median+idelta>kPedMax) continue;
      if (count<part*count1){
        count+=histo[median-idelta];
        mean +=histo[median-idelta]*(median-idelta);
        rms  +=histo[median-idelta]*(median-idelta)*(median-idelta);
        count+=histo[median+idelta];
        mean +=histo[median+idelta]*(median+idelta);
        rms  +=histo[median+idelta]*(median+idelta)*(median+idelta);
      }
    }
    fPadPedestal = 0;
    fPadNoise    = 0;
    if ( count > 0 ) {
      mean/=count;
      rms    = TMath::Sqrt(TMath::Abs(rms/count-mean*mean));
      fPadPedestal = mean;
      fPadNoise    = rms;
    }
  }
  fPadPedestal*=(Float_t)(!fIsZeroSuppressed);
}
//_____________________________________________________________________
void AliTPCCalibPulser::FindPulserSignal(TVectorD &param, Float_t &qSum)
{
//
    //  Find position, signal width and height of the CE signal (last signal)
    //  param[0] = Qmax, param[1] = mean time, param[2] = rms;
    //  maxima: array of local maxima of the pad signal use the one closest to the mean CE position
    //
  
  Float_t ceQmax  =0, ceQsum=0, ceTime=0, ceRMS=0;
  Int_t   cemaxpos       = fMaxTimeBin;
  Float_t ceSumThreshold = 10.*TMath::Max(fPadNoise,Float_t(1.));  // threshold for the signal sum
  Float_t ceMaxThreshold = 5.*TMath::Max(fPadNoise,Float_t(1.));  // threshold for the signal max
  const Int_t    kCemin  = fPeakIntMinus;             // range for the analysis of the ce signal +- channels from the peak
  const Int_t    kCemax  = fPeakIntPlus;
  param[0] = ceQmax;
  param[1] = ceTime;
  param[2] = ceRMS;
  qSum     = ceQsum;
  
  if (cemaxpos>0){
    ceQmax = fPadSignal.GetMatrixArray()[cemaxpos]-fPadPedestal;
    if ( ceQmax<ceMaxThreshold ) return;
    for (Int_t i=cemaxpos-kCemin; i<=cemaxpos+kCemax; ++i){
      Float_t signal = fPadSignal.GetMatrixArray()[i]-fPadPedestal;
      if ( (i>fFirstTimeBin) && (i<fLastTimeBin) && (signal>0) ){
        ceTime+=signal*(i+0.5);
        ceRMS +=signal*(i+0.5)*(i+0.5);
        ceQsum+=signal;
      }
    }
  }
  if (ceQsum>ceSumThreshold) {
    ceTime/=ceQsum;
    ceRMS  = TMath::Sqrt(TMath::Abs(ceRMS/ceQsum-ceTime*ceTime));
    ceTime-=GetL1PhaseTB();
        //only fill the Time0Offset if pad was not marked as an outlier!
    if ( !fOutliers ){
      //skip edge pads for calculating the mean time
      if ( !IsEdgePad(fCurrentSector,fCurrentRow,fCurrentPad) ){
        fVTime0Offset.GetMatrixArray()[fCurrentSector]+=ceTime;   // mean time for each sector
        fVTime0OffsetCounter.GetMatrixArray()[fCurrentSector]++;
        GetHistoTSec()->Fill(ceTime,fCurrentSector);
      }
    } else {
      if ( !(fOutliers->GetCalROC(fCurrentSector)->GetValue(fCurrentChannel)) ){
        fVTime0Offset.GetMatrixArray()[fCurrentSector]+=ceTime;   // mean time for each sector
        fVTime0OffsetCounter.GetMatrixArray()[fCurrentSector]++;
      }
    }
    
  //Normalise Q to the 'cell-size': The wire density is the same in the IROC and OROC, therefore the
  //                                the pick-up signal should scale with the pad area. In addition
  //                                the signal should decrease with the wire distance (4mm in IROC, 6mm in OROC),
  //                                ratio 2/3. The pad area we express in mm2 (factor 100). We normalise the signal
  //                                to the OROC signal (factor 2/3 for the IROCs).
    Float_t norm = fParam->GetPadPitchWidth(fCurrentSector)*fParam->GetPadPitchLength(fCurrentSector,fCurrentRow)*100;
    if ( fCurrentSector<fParam->GetNInnerSector() ) norm*=3./2.;
    
    ceQsum/=norm;
  } else {
    ceQmax=0;
    ceTime=0;
    ceRMS =0;
    ceQsum=0;
  }
  param[0] = ceQmax;
  param[1] = ceTime;
  param[2] = ceRMS;
  qSum     = ceQsum;
}
//_____________________________________________________________________
void AliTPCCalibPulser::ProcessPad()
{
  //
  //  Process data of current pad
  //
  
  FindPedestal();
  TVectorD param(3);
  Float_t  qSum;
  FindPulserSignal(param, qSum);
  
  Double_t meanT  = param[1];
  Double_t sigmaT = param[2];
  
  
  //Fill Event T0 counter
  (*GetPadTimesEvent(fCurrentSector,kTRUE)).GetMatrixArray()[fCurrentChannel] = meanT;
  
  //Fill Q histogram
//   GetHistoQ(fCurrentSector,kTRUE)->Fill( TMath::Sqrt(qSum), fCurrentChannel ); //use linear scale, needs also a change in the Analyse funciton.
  GetHistoQ(fCurrentSector,kTRUE)->Fill( qSum, fCurrentChannel );
  
  //Fill RMS histogram
  GetHistoRMS(fCurrentSector,kTRUE)->Fill( sigmaT, fCurrentChannel );
  
  
    //Fill debugging info
  if ( GetStreamLevel()>0 ){
    TTreeSRedirector *streamer=GetDebugStreamer();
    if ( GetStreamLevel() == 1 ){
      if ( streamer ) {
        Int_t padc = fCurrentPad-(fROC->GetNPads(fCurrentSector,fCurrentRow)/2);
        (*streamer) << "PadSignals" <<
          "Event="  <<fNevents <<
          "Sector=" <<fCurrentSector<<
          "Row="    <<fCurrentRow<<
          "Pad="    <<fCurrentPad<<
          "PadC="   <<padc<<
          "Channel="<<fCurrentChannel<<
          "Sum="    <<qSum<<
          "params.="<<&param<<
          "signal.=" <<&fPadSignal<<
          "\n";
      }
    } else { //debug > 1
      (*GetPadPedestalEvent(fCurrentSector,kTRUE))[fCurrentChannel]=fPadPedestal;
      (*GetPadRMSEvent(fCurrentSector,kTRUE))[fCurrentChannel]=sigmaT;
      (*GetPadQEvent(fCurrentSector,kTRUE))[fCurrentChannel]=qSum;
    }
  }
  ResetPad();
}
//_____________________________________________________________________
void AliTPCCalibPulser::EndEvent()
{
    //
    //  Process data of current event
    //
  
    //check if last pad has allready been processed, if not do so
  if ( fMaxTimeBin>-1 ) ProcessPad();
  
    //loop over all ROCs, fill Time0 histogram corrected for the mean Time0 of each ROC
  for ( Int_t iSec = 0; iSec<72; ++iSec ){
    TVectorF *vTimes = GetPadTimesEvent(iSec);
    if ( !vTimes || fVTime0OffsetCounter[iSec]==0 ) continue;
    Float_t time0 = fVTime0Offset[iSec]/fVTime0OffsetCounter[iSec];
    for ( UInt_t iChannel=0; iChannel<fROC->GetNChannels(iSec); ++iChannel ){
      Float_t time  = (*vTimes).GetMatrixArray()[iChannel];
      
      GetHistoT0(iSec,kTRUE)->Fill( time-time0,iChannel );
            //GetHistoT0(iSec,kTRUE)->Fill( time,iChannel );
      
      
      //Debug start
      if ( GetStreamLevel()>1 ){
        TTreeSRedirector *streamer=GetDebugStreamer();
        if ( streamer ) {
          Int_t row=0;
          Int_t pad=0;
          Int_t padc=0;
          
          Float_t q   = (*GetPadQEvent(iSec)).GetMatrixArray()[iChannel];
          Float_t rms = (*GetPadRMSEvent(iSec)).GetMatrixArray()[iChannel];
          
          UInt_t channel=iChannel;
          Int_t sector=iSec;
          
          while ( channel > (fROC->GetRowIndexes(sector)[row]+fROC->GetNPads(sector,row)-1) ) row++;
          pad = channel-fROC->GetRowIndexes(sector)[row];
          padc = pad-(fROC->GetNPads(sector,row)/2);
          
          (*streamer) << "DataPad" <<
    		    "Event=" << fNevents <<
            "Sector="<< sector <<
            "Row="   << row<<
            "Pad="   << pad <<
            "PadC="  << padc <<
            "PadSec="<< channel <<
            "Time0="  << time0 <<
            "Time="  << time <<
            "RMS="   << rms <<
            "Sum="   << q <<
            "\n";
        }
      }
      //Debug end
    }
  }
  ++fNevents;
}
//_____________________________________________________________________
TH2S* AliTPCCalibPulser::GetHisto(Int_t sector, TObjArray *arr,
                                  Int_t nbinsY, Float_t ymin, Float_t ymax,
                                  const Char_t *type, Bool_t force)
{
    //
    // return pointer to Q histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
  if ( !force || arr->UncheckedAt(sector) )
    return (TH2S*)arr->UncheckedAt(sector);
  
  // if we are forced and histogram doesn't yes exist create it
  // new histogram with Q calib information. One value for each pad!
  TH2S* hist = new TH2S(Form("hCalib%s%.2d",type,sector),Form("%s calibration histogram sector %.2d",type,sector),
                        nbinsY, ymin, ymax,
                        fROC->GetNChannels(sector),0,fROC->GetNChannels(sector));
  hist->SetDirectory(0);
  arr->AddAt(hist,sector);
  return hist;
}
//_____________________________________________________________________
TH2S* AliTPCCalibPulser::GetHistoT0(Int_t sector, Bool_t force)
{
    //
    // return pointer to T0 histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
  TObjArray *arr = &fHistoT0Array;
  return GetHisto(sector, arr, fNbinsT0, fXminT0, fXmaxT0, "T0", force);
}
//_____________________________________________________________________
TH2S* AliTPCCalibPulser::GetHistoQ(Int_t sector, Bool_t force)
{
    //
    // return pointer to Q histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
  TObjArray *arr = &fHistoQArray;
  return GetHisto(sector, arr, fNbinsQ, fXminQ, fXmaxQ, "Q", force);
}
//_____________________________________________________________________
TH2S* AliTPCCalibPulser::GetHistoRMS(Int_t sector, Bool_t force)
{
    //
    // return pointer to Q histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
  TObjArray *arr = &fHistoRMSArray;
  return GetHisto(sector, arr, fNbinsRMS, fXminRMS, fXmaxRMS, "RMS", force);
}
//_____________________________________________________________________
TH2F* AliTPCCalibPulser::GetHistoTSec()
{
    //
    // return the pointer to the abs time distribution per sector
    // create it if it does not exist
    //
  if ( !fHMeanTimeSector )   //!!!if you change the binning here, you should also change it in the Analyse Function!!
    fHMeanTimeSector = new TH2F("fHMeanTimeSector","Abs mean time per sector",
                                20*(fLastTimeBin-fFirstTimeBin), fFirstTimeBin, fLastTimeBin,
                                72,0,72);
  return fHMeanTimeSector;
}
//_____________________________________________________________________
TVectorF* AliTPCCalibPulser::GetPadInfoEvent(Int_t sector, TObjArray *arr, Bool_t force)
{
    //
    // return pointer to Pad Info from 'arr' for the current event and sector
    // if force is true create it if it doesn't exist allready
    //
  if ( !force || arr->UncheckedAt(sector) )
    return (TVectorF*)arr->UncheckedAt(sector);
  
  TVectorF *vect = new TVectorF(fROC->GetNChannels(sector));
  arr->AddAt(vect,sector);
  return vect;
}
//_____________________________________________________________________
TVectorF* AliTPCCalibPulser::GetPadTimesEvent(Int_t sector, Bool_t force)
{
    //
    // return pointer to Pad Times Array for the current event and sector
    // if force is true create it if it doesn't exist allready
    //
  TObjArray *arr = &fPadTimesArrayEvent;
  return GetPadInfoEvent(sector,arr,force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibPulser::GetPadQEvent(Int_t sector, Bool_t force)
{
    //
    // return pointer to Pad Q Array for the current event and sector
    // if force is true create it if it doesn't exist allready
    // for debugging purposes only
    //
  
  TObjArray *arr = &fPadQArrayEvent;
  return GetPadInfoEvent(sector,arr,force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibPulser::GetPadRMSEvent(Int_t sector, Bool_t force)
{
    //
    // return pointer to Pad RMS Array for the current event and sector
    // if force is true create it if it doesn't exist allready
    // for debugging purposes only
    //
  TObjArray *arr = &fPadRMSArrayEvent;
  return GetPadInfoEvent(sector,arr,force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibPulser::GetPadPedestalEvent(Int_t sector, Bool_t force)
{
    //
    // return pointer to Pad RMS Array for the current event and sector
    // if force is true create it if it doesn't exist allready
    // for debugging purposes only
    //
  TObjArray *arr = &fPadPedestalArrayEvent;
  return GetPadInfoEvent(sector,arr,force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPulser::GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force) const
{
    //
    // return pointer to ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
  if ( !force || arr->UncheckedAt(sector) )
    return (AliTPCCalROC*)arr->UncheckedAt(sector);
  
    // if we are forced and histogram doesn't yes exist create it
  
    // new AliTPCCalROC for T0 information. One value for each pad!
  AliTPCCalROC *croc = new AliTPCCalROC(sector);
  arr->AddAt(croc,sector);
  return croc;
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPulser::GetCalRocT0(Int_t sector, Bool_t force)
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
  TObjArray *arr = &fCalRocArrayT0;
  return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPulser::GetCalRocQ(Int_t sector, Bool_t force)
{
    //
    // return pointer to T0 ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
  TObjArray *arr = &fCalRocArrayQ;
  return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPulser::GetCalRocRMS(Int_t sector, Bool_t force)
{
    //
    // return pointer to signal width ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
  TObjArray *arr = &fCalRocArrayRMS;
  return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibPulser::GetCalRocOutliers(Int_t sector, Bool_t force)
{
    //
    // return pointer to Outliers
    // if force is true create a new histogram if it doesn't exist allready
    //
  TObjArray *arr = &fCalRocArrayOutliers;
  return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
void AliTPCCalibPulser::ResetEvent()
{
    //
    //  Reset global counters  -- Should be called before each event is processed
    //
  fLastSector=-1;
  fCurrentSector=-1;
  fCurrentRow=-1;
  fCurrentPad=-1;
  fCurrentChannel=-1;
  
  ResetPad();
  
  fPadTimesArrayEvent.Delete();
  
  fPadQArrayEvent.Delete();
  fPadRMSArrayEvent.Delete();
  fPadPedestalArrayEvent.Delete();
  
  for ( Int_t i=0; i<72; ++i ){
    fVTime0Offset[i]=0;
    fVTime0OffsetCounter[i]=0;
  }
}
//_____________________________________________________________________
void AliTPCCalibPulser::ResetPad()
{
    //
    //  Reset pad infos -- Should be called after a pad has been processed
    //
  for (Int_t i=fFirstTimeBin; i<fLastTimeBin+1; ++i)
    fPadSignal[i] = 0;
  fMaxTimeBin = -1;
  fMaxPadSignal = -1;
  fPadPedestal  = -1;
  fPadNoise     = -1;
}
//_____________________________________________________________________
Bool_t AliTPCCalibPulser::IsEdgePad(Int_t sector, Int_t row, Int_t pad)
{
    //
    // return true if pad is on the edge of a row
    //
  Int_t edge1   = 0;
  Int_t edge2   = fROC->GetNPads(sector,row)-1;
  if ( pad == edge1 || pad == edge2 ) return kTRUE;
  
  return kFALSE;
}
//_____________________________________________________________________
void AliTPCCalibPulser::Merge(AliTPCCalibPulser * const sig)
{
  //
  //  Merge reference histograms of sig to the current AliTPCCalibPulser
  //

  MergeBase(sig);
  //merge histograms
  for (Int_t iSec=0; iSec<72; ++iSec){
    TH2S *hRefQmerge   = sig->GetHistoQ(iSec);
    TH2S *hRefT0merge  = sig->GetHistoT0(iSec);
    TH2S *hRefRMSmerge = sig->GetHistoRMS(iSec);
    
    
    if ( hRefQmerge ){
      TDirectory *dir = hRefQmerge->GetDirectory(); hRefQmerge->SetDirectory(0);
      TH2S *hRefQ   = GetHistoQ(iSec);
      if ( hRefQ ) hRefQ->Add(hRefQmerge);
      else {
        TH2S *hist = new TH2S(*hRefQmerge);
        hist->SetDirectory(0);
        fHistoQArray.AddAt(hist, iSec);
      }
      hRefQmerge->SetDirectory(dir);
    }
    if ( hRefT0merge ){
      TDirectory *dir = hRefT0merge->GetDirectory(); hRefT0merge->SetDirectory(0);
      TH2S *hRefT0  = GetHistoT0(iSec);
      if ( hRefT0 ) hRefT0->Add(hRefT0merge);
      else {
        TH2S *hist = new TH2S(*hRefT0merge);
        hist->SetDirectory(0);
        fHistoT0Array.AddAt(hist, iSec);
      }
      hRefT0merge->SetDirectory(dir);
    }
    if ( hRefRMSmerge ){
      TDirectory *dir = hRefRMSmerge->GetDirectory(); hRefRMSmerge->SetDirectory(0);
      TH2S *hRefRMS = GetHistoRMS(iSec);
      if ( hRefRMS ) hRefRMS->Add(hRefRMSmerge);
      else {
        TH2S *hist = new TH2S(*hRefRMSmerge);
        hist->SetDirectory(0);
        fHistoRMSArray.AddAt(hist, iSec);
      }
      hRefRMSmerge->SetDirectory(dir);
    }
    
  }
  if ( sig->fHMeanTimeSector ){
    TDirectory *dir = sig->fHMeanTimeSector->GetDirectory(); sig->fHMeanTimeSector->SetDirectory(0);
    if ( fHMeanTimeSector ) fHMeanTimeSector->Add(sig->fHMeanTimeSector);
    else {
      fHMeanTimeSector = new TH2F(*sig->fHMeanTimeSector);
      fHMeanTimeSector->SetDirectory(0);
    }
    sig->fHMeanTimeSector->SetDirectory(dir);
  }
}


//_____________________________________________________________________
Long64_t AliTPCCalibPulser::Merge(TCollection * const list)
{
  //
  // Merge all objects of this type in list
  //
  
  Long64_t nmerged=1;
  
  TIter next(list);
  AliTPCCalibPulser *ce=0;
  TObject *o=0;
  
  while ( (o=next()) ){
    ce=dynamic_cast<AliTPCCalibPulser*>(o);
    if (ce){
      Merge(ce);
      ++nmerged;
    }
  }
  
  return nmerged;
}

//_____________________________________________________________________
void AliTPCCalibPulser::Analyse()
{
  //
  //  Calculate calibration constants
  //
  
  TVectorD paramQ(3);
  TVectorD paramT0(3);
  TVectorD paramRMS(3);
  TMatrixD dummy(3,3);
  //calculate mean time for each sector and mean time for each side
  TH1F hMeanTsec("hMeanTsec","hMeanTsec",20*(fLastTimeBin-fFirstTimeBin),fFirstTimeBin,fLastTimeBin);
  fVMeanTimeSector.Zero();
  
  for (Int_t iSec=0; iSec<72; ++iSec){
    TH2S *hT0 = GetHistoT0(iSec);
    if (!hT0 ) continue;
    //calculate sector mean T
    if ( fHMeanTimeSector ){
      Int_t nbinsT = fHMeanTimeSector->GetNbinsX();
      Int_t offset = (nbinsT+2)*(iSec+1);
      Float_t *arrP=fHMeanTimeSector->GetArray()+offset;
      Int_t entries=0;
      for ( Int_t i=0; i<nbinsT; i++ ) entries+=(Int_t)arrP[i+1];
      hMeanTsec.Set(nbinsT+2,arrP);
      hMeanTsec.SetEntries(entries);
      paramT0.Zero();
      // truncated mean: remove lower 5% and upper 5%
      if ( entries>0 ) AliMathBase::TruncatedMean(&hMeanTsec,&paramT0,0.05,.95);
      fVMeanTimeSector[iSec]=paramT0[1];
    }
    
    AliTPCCalROC *rocQ   = GetCalRocQ  (iSec,kTRUE);
    AliTPCCalROC *rocT0  = GetCalRocT0 (iSec,kTRUE);
    AliTPCCalROC *rocRMS = GetCalRocRMS(iSec,kTRUE);
    AliTPCCalROC *rocOut = GetCalRocOutliers(iSec,kTRUE);
    
    TH2S *hQ   = GetHistoQ(iSec);
    TH2S *hRMS = GetHistoRMS(iSec);
    
    Short_t *arrayhQ   = hQ->GetArray();
    Short_t *arrayhT0  = hT0->GetArray();
    Short_t *arrayhRMS = hRMS->GetArray();
    
    UInt_t nChannels = fROC->GetNChannels(iSec);
    Float_t meanTsec = fVMeanTimeSector[iSec];
    
  //debug
    Int_t row=0;
    Int_t pad=0;
    Int_t padc=0;
  //! debug
    
    for (UInt_t iChannel=0; iChannel<nChannels; ++iChannel){
      
      Float_t cogTime0 = -1000;
      Float_t cogQ     = -1000;
      Float_t cogRMS   = -1000;
      Float_t cogOut   = 0;
      
      Int_t offsetQ = (fNbinsQ+2)*(iChannel+1)+1;
      Int_t offsetT0 = (fNbinsT0+2)*(iChannel+1)+1;
      Int_t offsetRMS = (fNbinsRMS+2)*(iChannel+1)+1;
/*
      AliMathBase::FitGaus(arrayhQ+offsetQ,fNbinsQ,fXminQ,fXmaxQ,&paramQ,&dummy);
      AliMathBase::FitGaus(arrayhT0+offsetT0,fNbinsT0,fXminT0,fXmaxT0,&paramT0,&dummy);
      AliMathBase::FitGaus(arrayhRMS+offsetRMS,fNbinsRMS,fXminRMS,fXmaxRMS,&paramRMS,&dummy);
      cogQ     = paramQ[1];
      cogTime0 = paramT0[1];
      cogRMS   = paramRMS[1];
*/
      cogQ     = AliMathBase::GetCOG(arrayhQ+offsetQ,fNbinsQ,fXminQ,fXmaxQ);
      cogTime0 = AliMathBase::GetCOG(arrayhT0+offsetT0,fNbinsT0,fXminT0,fXmaxT0);
      cogRMS   = AliMathBase::GetCOG(arrayhRMS+offsetRMS,fNbinsRMS,fXminRMS,fXmaxRMS);
      
      /*
      if ( (cogQ < ??) && (cogTime0 > ??) && (cogTime0<??) && ( cogRMS>??) ){
    cogOut = 1;
    cogTime0 = 0;
    cogQ     = 0;
    cogRMS   = 0;
      }
*/
//       rocQ->SetValue(iChannel, cogQ*cogQ); // changed to linear scale again
      rocQ->SetValue(iChannel, cogQ);
      rocT0->SetValue(iChannel, cogTime0+meanTsec); //offset by mean time of the sector
      rocRMS->SetValue(iChannel, cogRMS);
      rocOut->SetValue(iChannel, cogOut);

      // in case a channel has no data set the value to 0
      if (TMath::Abs(cogTime0-fXminT0)<1e-10){
        rocQ->SetValue(iChannel, 0);
        rocT0->SetValue(iChannel, 0); //offset by mean time of the sector
        rocRMS->SetValue(iChannel, 0);
      }
      
      //debug
      if ( GetStreamLevel() > 2 ){
        TTreeSRedirector *streamer=GetDebugStreamer();
        if ( streamer ) {
          while ( iChannel > (fROC->GetRowIndexes(iSec)[row]+fROC->GetNPads(iSec,row)-1) ) row++;
          pad = iChannel-fROC->GetRowIndexes(iSec)[row];
          padc = pad-(fROC->GetNPads(iSec,row)/2);
          
          (*streamer) << "DataEnd" <<
            "Sector="  << iSec      <<
            "Pad="     << pad       <<
            "PadC="    << padc      <<
            "Row="     << row       <<
            "PadSec="  << iChannel   <<
            "Q="       << cogQ      <<
            "T0="      << cogTime0  <<
            "RMS="     << cogRMS    <<
            "\n";
        }
      }
      //! debug
    }
    
    
  }
}
//_____________________________________________________________________
//_________________________  Test Functions ___________________________
//_____________________________________________________________________
TObjArray* AliTPCCalibPulser::TestBinning()
{
  //
  //  Function to test the binning of the reference histograms
  //  type: T0, Q or RMS
  //  mode: 0 - number of filled bins per channel
  //        1 - number of empty bins between filled bins in one ROC
  //  returns TObjArray with the test histograms type*2+mode:
  //  position 0 = T0,0 ; 1 = T0,1 ; 2 = Q,0 ...
  
  
  TObjArray *histArray = new TObjArray(6);
  const Char_t *type[] = {"T0","Q","RMS"};
  Int_t fNbins[3] = {fNbinsT0,fNbinsQ,fNbinsRMS};
  
  for (Int_t itype = 0; itype<3; ++itype){
    for (Int_t imode=0; imode<2; ++imode){
      Int_t icount = itype*2+imode;
      histArray->AddAt(new TH1F(Form("hTestBinning%s%d",type[itype],imode),
                                Form("Test Binning of '%s', mode - %d",type[itype],imode),
                                72,0,72),
                       icount);
    }
  }
  
  
  TH2S *hRef=0x0;
  Short_t *array=0x0;
  for (Int_t itype = 0; itype<3; ++itype){
    for (Int_t iSec=0; iSec<72; ++iSec){
      if ( itype == 0 ) hRef = GetHistoT0(iSec);
      if ( itype == 1 ) hRef = GetHistoQ(iSec);
      if ( itype == 2 ) hRef = GetHistoRMS(iSec);
      if ( hRef == 0x0 ) continue;
      array = (hRef->GetArray());
      UInt_t nChannels = fROC->GetNChannels(iSec);
      
      Int_t nempty=0;
      for (UInt_t iChannel=0; iChannel<nChannels; ++iChannel){
        Int_t nfilled=0;
        Int_t offset = (fNbins[itype]+2)*(iChannel+1)+1;
        Int_t c1 = 0;
        Int_t c2 = 0;
        for (Int_t iBin=0; iBin<fNbins[itype]; ++iBin){
          if ( array[offset+iBin]>0 ) {
            nfilled++;
            if ( c1 && c2 ) nempty++;
            else c1 = 1;
          }
          else if ( c1 ) c2 = 1;
          
        }
        ((TH1F*)histArray->At(itype*2))->Fill(nfilled);
      }
      ((TH1F*)histArray->At(itype*2+1))->Fill(iSec,nempty);
    }
  }
  return histArray;
}
