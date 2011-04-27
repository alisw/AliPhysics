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

////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//             Implementation of the TPC Central Electrode calibration                //
//                                                                                    //
//   Origin: Jens Wiechula, Marian Ivanov   J.Wiechula@gsi.de, Marian.Ivanov@cern.ch  //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////
//
//
// *************************************************************************************
// *                                Class Description                                  *
// *************************************************************************************
//
/* BEGIN_HTML
 <h4>The AliTPCCalibCE class is used to get calibration data from the Central Electrode
 using laser runs.</h4>

 The information retrieved is
 <ul style="list-style-type: square;">
   <li>Time arrival from the CE</li>
   <li>Signal width</li>
   <li>Signal sum</li>
 </ul>

<h4>Overview:</h4>
 <ol style="list-style-type: upper-roman;">
   <li><a href="#working">Working principle</a></li>
   <li><a href="#user">User interface for filling data</a></li>
   <li><a href="#info">Stored information</a></li>
 </ol>

 <h3><a name="working">I. Working principle</a></h3>

 <h4>Raw laser data is processed by calling one of the ProcessEvent(...) functions
 (see below). These in the end call the Update(...) function.</h4>

 <ul style="list-style-type: square;">
   <li>the Update(...) function:<br />
       In this function the array fPadSignal is filled with the adc signals between the specified range
       fFirstTimeBin and fLastTimeBin for the current pad.
       before going to the next pad the ProcessPad() function is called, which analyses the data for one pad
       stored in fPadSignal.
   </li>
   <ul style="list-style-type: square;">
   <li>the ProcessPad() function:</li>
   <ol style="list-style-type: decimal;">
     <li>Find Pedestal and Noise information</li>
     <ul style="list-style-type: square;">
       <li>use database information which has to be set by calling<br />
           SetPedestalDatabase(AliTPCCalPad *pedestalTPC, AliTPCCalPad *padNoiseTPC)</li>
       <li>if no information from the pedestal data base
           is available the informaion is calculated on the fly
           ( see FindPedestal() function )</li>
     </ul>
     <li>Find local maxima of the pad signal</li>
     <ul style="list-style-type: square;">
       <li>maxima arise from the laser tracks, the CE and also periodic postpeaks after the CE signal have
       have been observed ( see FindLocalMaxima(...) )</li>
     </ul>
     <li>Find the CE signal information</li>
     <ul style="list-style-type: square;">
       <li>to find the position of the CE signal the Tmean information from the previos event is used
           as the CE signal the local maximum closest to this Tmean is identified</li>
       <li>calculate  mean = T0, RMS = signal width and Q sum in a range of -4+7 timebins around Q max position
           the Q sum is scaled by pad area (see FindPulserSignal(...) function)</li>
     </ul>
     <li>Fill a temprary array for the T0 information (GetPadTimesEvent(fCurrentSector,kTRUE)) (why see below)</li>
     <li>Fill the Q sum and RMS values in the histograms (GetHisto[RMS,Q](ROC,kTRUE))</li>
     </ol>
   </ul>
 </ul>

 <h4>At the end of each event the EndEvent() function is called</h4>

 <ul style="list-style-type: square;">
   <li>the EndEvent() function:</li>
   <ul style="list-style-type: square;">
     <li>calculate the mean T0 for side A and side C. Fill T0 histogram with Time0-<Time0 for side[A,C]>
         This is done to overcome syncronisation problems between the trigger and the fec clock.</li>
     <li>calculate Mean T for each ROC using the COG aroud the median of the LocalMaxima distribution in one sector</li>
     <li>calculate Mean Q</li>
     <li>calculate Global fit parameters for Pol1 and Pol2 fits</li>
   </ul>
 </ul>

 <h4>After accumulating the desired statistics the Analyse() function has to be called.</h4>
  <ul style="list-style-type: square;">
  <li>the Analyse() function:</li>
    <ul style="list-style-type: square;">
      <li>calculate the mean values of T0, RMS, Q for each pad, using
          the AliMathBase::GetCOG(...) function</li>
      <li>fill the calibration storage classes (AliTPCCalROC) for each ROC</li>
         (The calibration information is stored in the TObjArrays fCalRocArrayT0, fCalRocArrayRMS and fCalRocArrayQ</li>
    </ul>
  </ul>

 <h3><a name="user">II. User interface for filling data</a></h3>

 <h4>To Fill information one of the following functions can be used:</h4>

 <ul style="list-style-type: none;">
  <li> Bool_t ProcessEvent(eventHeaderStruct *event);</li>
    <ul style="list-style-type: square;">
      <li>process Date event</li>
      <li>use AliTPCRawReaderDate and call ProcessEvent(AliRawReader *rawReader)</li>
    </ul>
    <br />

  <li> Bool_t ProcessEvent(AliRawReader *rawReader);</li>
    <ul style="list-style-type: square;">
      <li>process AliRawReader event</li>
      <li>use AliTPCRawStream to loop over data and call ProcessEvent(AliTPCRawStream *rawStream)</li>
    </ul>
    <br />

  <li> Bool_t ProcessEvent(AliTPCRawStream *rawStream);</li>
    <ul style="list-style-type: square;">
      <li>process event from AliTPCRawStream</li>
      <li>call Update function for signal filling</li>
    </ul>
    <br />

  <li> Int_t Update(const Int_t isector, const Int_t iRow, const Int_t
              iPad,  const Int_t iTimeBin, const Float_t signal);</li>
    <ul style="list-style-type: square;">
      <li>directly  fill signal information (sector, row, pad, time bin, pad)
          to the reference histograms</li>
    </ul>
 </ul>

 <h4>It is also possible to merge two independently taken calibrations using the function</h4>

 <ul style="list-style-type: none;">
  <li> void Merge(AliTPCCalibSignal *sig)</li>
    <ul style="list-style-type: square;">
      <li>copy histograms in 'sig' if they do not exist in this instance</li>
      <li>Add histograms in 'sig' to the histograms in this instance if the allready exist</li>
      <li>After merging call Analyse again!</li>
    </ul>
 </ul>


 <h4>example: filling data using root raw data:</h4>
 <pre> 
 void fillCE(Char_t *filename)
 {
    rawReader = new AliRawReaderRoot(fileName);
    if ( !rawReader ) return;
    AliTPCCalibCE *calib = new AliTPCCalibCE;
    while (rawReader->NextEvent()){
      calib->ProcessEvent(rawReader);
    }
    calib->Analyse();
    calib->DumpToFile("CEData.root");
    delete rawReader;
    delete calib;
 }
 </pre>

 <h3><a name="info">III. What kind of information is stored and how to retrieve it</a></h4>

 <h4><a name="info:stored">III.1 Stored information</a></h4>
 <ul style="list-style-type: none;">
  <li>Histograms:</li>
  <ul style="list-style-type: none;">
    <li>For each ROC three TH2S histos 'Reference Histograms'  (ROC channel vs. [Time0, signal width, Q sum])
        is created when it is filled for the first time (GetHisto[T0,RMS,Q](ROC,kTRUE)). The histos are
        stored in the TObjArrays fHistoT0Array, fHistoRMSArray and fHistoQArray.</li>
  </ul>
  <br />

 <li>Calibration Data:</li>
 <ul style="list-style-type: none;">
      <li>For each ROC three types of calibration data (AliTPCCalROC) is stored: for the mean arrival Time,
          the signal width and the signal Sum. The AliTPCCalROC objects are stored in the TObjArrays
          fCalRocArrayT0, fCalRocArrayRMS , fCalRocArrayQ. The object for each roc is created the first time it
          is accessed (GetCalRoc[T0,RMS,Q](ROC,kTRUE));</li>
 </ul>
 <br />

 <li>For each event the following information is stored:</li>
   
 <ul style="list-style-type: square;">
   <li>event time ( TVectorD  fVEventTime )</li>
   <li>event id   ( TVectorD  fVEventNumber )</li>
   <br />
   <li>mean arrival time for each ROC                ( TObjArray fTMeanArrayEvent )</li>
   <li>mean Q for each ROC                           ( TObjArray fQMeanArrayEvent )</li>
   <li>parameters of a plane fit for each ROC        ( TObjArray fParamArrayEventPol1 )</li>
   <li>parameters of a 2D parabola fit for each ROC  ( TObjArray fParamArrayEventPol2 )</li>
  </ul>
 </ul>

 <h4><a name="info:retrieve">III.2 Retrieving information</a></h4>
 <ul style="list-style-type: none;">
  <li>Accessing the 'Reference Histograms' (Time0, signal width and Q sum information pad by pad):</li>
    <ul style="list-style-type: square;">
      <li>TH2F *GetHistoT0(Int_t sector);</li>
      <li>TH2F *GetHistoRMS(Int_t sector);</li>
      <li>TH2F *GetHistoQ(Int_t sector);</li>
    </ul>
    <br />
    
  <li>Accessing the calibration storage objects:</li>
    <ul style="list-style-type: square;">
      <li>AliTPCCalROC *GetCalRocT0(Int_t sector);   // for the Time0 values</li>
      <li>AliTPCCalROC *GetCalRocRMS(Int_t sector);  // for the signal width values</li>
      <li>AliTPCCalROC *GetCalRocQ(Int_t sector);    // for the Q sum values</li>
    </ul>
    <br />

  <li>Accessin the event by event information:</li>
    <ul style="list-style-type: square;">
      <li>The event by event information can be displayed using the</li>
      <li>MakeGraphTimeCE(Int_t sector, Int_t xVariable, Int_t fitType, Int_t fitParameter)</li>
      <li>which creates a graph from the specified variables</li>
    </ul>
  </ul>
  
  <h4>example for visualisation:</h4>
  <pre>
  //if the file "CEData.root" was created using the above example one could do the following:
  TFile fileCE("CEData.root")
  AliTPCCalibCE *ce = (AliTPCCalibCE*)fileCE->Get("AliTPCCalibCE");
  ce->GetCalRocT0(0)->Draw("colz");
  ce->GetCalRocRMS(0)->Draw("colz");

  //or use the AliTPCCalPad functionality:
  AliTPCCalPad padT0(ped->GetCalPadT0());
  AliTPCCalPad padSigWidth(ped->GetCalPadRMS());
  padT0->MakeHisto2D()->Draw("colz");       //Draw A-Side Time0 Information
  padSigWidth->MakeHisto2D()->Draw("colz"); //Draw A-Side signal width Information

  //display event by event information:
  //Draw mean arrival time as a function of the event time for oroc sector A00
  ce->MakeGraphTimeCE(36, 0, 2)->Draw("alp");
  //Draw first derivative in local x from a plane fit as a function of the event time for oroc sector A00
  ce->MakeGraphTimeCE(36, 0, 0, 1)->Draw("alp");  
  </pre>
END_HTML */
//////////////////////////////////////////////////////////////////////////////////////


//Root includes
#include <TObjArray.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2S.h>
#include <TF1.h>
#include <TString.h>
#include <TVectorF.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TMap.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCollection.h>
#include <TTimeStamp.h>
#include <TList.h>
#include <TKey.h>

//AliRoot includes
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliRawEventHeaderBase.h"
#include "AliTPCRawStream.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCROC.h"
#include "AliTPCParam.h"
#include "AliTPCCalibCE.h"
#include "AliMathBase.h"
#include "AliTPCTransform.h"
#include "AliTPCLaserTrack.h"
#include "TTreeStream.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
//date
#include "event.h"
ClassImp(AliTPCCalibCE)


AliTPCCalibCE::AliTPCCalibCE() :
  AliTPCCalibRawBase(),
  fNbinsT0(200),
  fXminT0(-5),
  fXmaxT0(5),
  fNbinsQ(200),
  fXminQ(1),
  fXmaxQ(40),
  fNbinsRMS(100),
  fXminRMS(0.1),
  fXmaxRMS(5.1),
  fPeakDetMinus(2),
  fPeakDetPlus(3),
  fPeakIntMinus(2),
  fPeakIntPlus(2),
  fNoiseThresholdMax(5.),
  fNoiseThresholdSum(8.),
  fIsZeroSuppressed(kFALSE),
  fLastSector(-1),
  fSecRejectRatio(.4),
  fParam(new AliTPCParam),
  fPedestalTPC(0x0),
  fPadNoiseTPC(0x0),
  fPedestalROC(0x0),
  fPadNoiseROC(0x0),
  fCalRocArrayT0(72),
  fCalRocArrayT0Err(72),
  fCalRocArrayQ(72),
  fCalRocArrayRMS(72),
  fCalRocArrayOutliers(72),
  fHistoQArray(72),
  fHistoT0Array(72),
  fHistoRMSArray(72),
  fMeanT0rms(0),
  fMeanQrms(0),
  fMeanRMSrms(0),
  fHistoTmean(72),
  fParamArrayEventPol1(72),
  fParamArrayEventPol2(72),
  fTMeanArrayEvent(72),
  fQMeanArrayEvent(72),
  fVEventTime(1000),
  fVEventNumber(1000),
  fVTime0SideA(1000),
  fVTime0SideC(1000),
  fEventId(-1),
  fOldRunNumber(0),
  fPadTimesArrayEvent(72),
  fPadQArrayEvent(72),
  fPadRMSArrayEvent(72),
  fPadPedestalArrayEvent(72),
  fCurrentChannel(-1),
  fCurrentSector(-1),
  fCurrentRow(-1),
  fMaxPadSignal(-1),
  fMaxTimeBin(-1),
//   fPadSignal(1024),
  fPadPedestal(0),
  fPadNoise(0),
  fVTime0Offset(72),
  fVTime0OffsetCounter(72),
  fVMeanQ(72),
  fVMeanQCounter(72),
  fCurrentCETimeRef(0),
  fProcessOld(kTRUE),
  fProcessNew(kFALSE),
  fAnalyseNew(kTRUE),
  fHnDrift(0x0),
  fArrHnDrift(100),
  fTimeBursts(100),
  fArrFitGraphs(0x0)
{
  //
  // AliTPCSignal default constructor
  //
  SetNameTitle("AliTPCCalibCE","AliTPCCalibCE");
  fFirstTimeBin=650;
  fLastTimeBin=1030;
  fParam->Update();
  for (Int_t i=0;i<1024;++i) fPadSignal[i]=0;
  for (Int_t i=0;i<5;++i){
    fPeaks[i]=0;
    fPeakWidths[i]=0;
  }
  for (Int_t i=0; i<100; ++i) fBinsLastAna[i]=0;
}
//_____________________________________________________________________
AliTPCCalibCE::AliTPCCalibCE(const AliTPCCalibCE &sig) :
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
  fPeakDetMinus(sig.fPeakDetMinus),
  fPeakDetPlus(sig.fPeakDetPlus),
  fPeakIntMinus(sig.fPeakIntMinus),
  fPeakIntPlus(sig.fPeakIntPlus),
  fNoiseThresholdMax(sig.fNoiseThresholdMax),
  fNoiseThresholdSum(sig.fNoiseThresholdSum),
  fIsZeroSuppressed(sig.fIsZeroSuppressed),
  fLastSector(-1),
  fSecRejectRatio(.4),
  fParam(new AliTPCParam),
  fPedestalTPC(0x0),
  fPadNoiseTPC(0x0),
  fPedestalROC(0x0),
  fPadNoiseROC(0x0),
  fCalRocArrayT0(72),
  fCalRocArrayT0Err(72),
  fCalRocArrayQ(72),
  fCalRocArrayRMS(72),
  fCalRocArrayOutliers(72),
  fHistoQArray(72),
  fHistoT0Array(72),
  fHistoRMSArray(72),
  fMeanT0rms(sig.fMeanT0rms),
  fMeanQrms(sig.fMeanQrms),
  fMeanRMSrms(sig.fMeanRMSrms),
  fHistoTmean(72),
  fParamArrayEventPol1(72),
  fParamArrayEventPol2(72),
  fTMeanArrayEvent(72),
  fQMeanArrayEvent(72),
  fVEventTime(sig.fVEventTime),
  fVEventNumber(sig.fVEventNumber),
  fVTime0SideA(sig.fVTime0SideA),
  fVTime0SideC(sig.fVTime0SideC),
  fEventId(-1),
  fOldRunNumber(0),
  fPadTimesArrayEvent(72),
  fPadQArrayEvent(72),
  fPadRMSArrayEvent(72),
  fPadPedestalArrayEvent(72),
  fCurrentChannel(-1),
  fCurrentSector(-1),
  fCurrentRow(-1),
  fMaxPadSignal(-1),
  fMaxTimeBin(-1),
//   fPadSignal(1024),
  fPadPedestal(0),
  fPadNoise(0),
  fVTime0Offset(72),
  fVTime0OffsetCounter(72),
  fVMeanQ(72),
  fVMeanQCounter(72),
  fCurrentCETimeRef(0),
  fProcessOld(sig.fProcessOld),
  fProcessNew(sig.fProcessNew),
  fAnalyseNew(sig.fAnalyseNew),
  fHnDrift(0x0),
  fArrHnDrift(100),
  fTimeBursts(100),
  fArrFitGraphs(0x0)
{
  //
  // AliTPCSignal copy constructor
  //
  for (Int_t i=0;i<1024;++i) fPadSignal[i]=0;
  
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
      fHistoT0Array.AddAt(hNew,iSec);
    }
    if ( hRMS != 0x0 ){
      TH2S *hNew = new TH2S(*hRMS);
      hNew->SetDirectory(0);
      fHistoRMSArray.AddAt(hNew,iSec);
    }
  }

  //copy fit parameters event by event
  TObjArray *arr=0x0;
  for (Int_t iSec=0; iSec<72; ++iSec){
    arr = (TObjArray*)sig.fParamArrayEventPol1.UncheckedAt(iSec);
    if ( arr ){
      TObjArray *arrEvents = new TObjArray(arr->GetSize());
      fParamArrayEventPol1.AddAt(arrEvents, iSec);
      for (Int_t iEvent=0; iEvent<arr->GetSize(); ++iEvent)
        if ( TVectorD *vec=(TVectorD*)arr->UncheckedAt(iEvent) )
          arrEvents->AddAt(new TVectorD(*vec),iEvent);
    }

    arr = (TObjArray*)sig.fParamArrayEventPol2.UncheckedAt(iSec);
    if ( arr ){
      TObjArray *arrEvents = new TObjArray(arr->GetSize());
      fParamArrayEventPol2.AddAt(arrEvents, iSec);
      for (Int_t iEvent=0; iEvent<arr->GetSize(); ++iEvent)
        if ( TVectorD *vec=(TVectorD*)arr->UncheckedAt(iEvent) )
          arrEvents->AddAt(new TVectorD(*vec),iEvent);
    }

    TVectorF *vMeanTime = (TVectorF*)sig.fTMeanArrayEvent.UncheckedAt(iSec);
    TVectorF *vMeanQ    = (TVectorF*)sig.fQMeanArrayEvent.UncheckedAt(iSec);
    if ( vMeanTime )
      fTMeanArrayEvent.AddAt(new TVectorF(*vMeanTime), iSec);
    if ( vMeanQ )
      fQMeanArrayEvent.AddAt(new TVectorF(*vMeanQ), iSec);
  }


  fVEventTime.ResizeTo(sig.fVEventTime);
  fVEventNumber.ResizeTo(sig.fVEventNumber);
  fVEventTime.SetElements(sig.fVEventTime.GetMatrixArray());
  fVEventNumber.SetElements(sig.fVEventNumber.GetMatrixArray());

  fParam->Update();

  for (Int_t i=0; i<sig.fArrHnDrift.GetEntries();++i){
    TObject *o=sig.fArrHnDrift.UncheckedAt(i);
    if (o){
      TObject *newo=o->Clone("fHnDrift");
      fArrHnDrift.AddAt(newo,i);
      if (sig.fHnDrift && o==sig.fHnDrift) fHnDrift=(THnSparseI*)newo;
    }
  }

  for (Int_t i=0;i<sig.fTimeBursts.GetNrows();++i){
    fTimeBursts[i]=sig.fTimeBursts[i];
  }
  
  for (Int_t i=0;i<5;++i){
    fPeaks[i]=sig.fPeaks[i];
    fPeakWidths[i]=sig.fPeakWidths[i];
  }
  if (sig.fArrFitGraphs) {
    fArrFitGraphs=(TObjArray*)sig.fArrFitGraphs->Clone();
    fArrFitGraphs->SetOwner();
  }

  for (Int_t i=0; i<100; ++i) fBinsLastAna[i]=sig.fBinsLastAna[i];
  
}
//_____________________________________________________________________
AliTPCCalibCE::AliTPCCalibCE(const TMap *config) :
  AliTPCCalibRawBase(),
  fNbinsT0(200),
  fXminT0(-5),
  fXmaxT0(5),
  fNbinsQ(200),
  fXminQ(1),
  fXmaxQ(40),
  fNbinsRMS(100),
  fXminRMS(0.1),
  fXmaxRMS(5.1),
  fPeakDetMinus(2),
  fPeakDetPlus(3),
  fPeakIntMinus(2),
  fPeakIntPlus(2),
  fNoiseThresholdMax(5.),
  fNoiseThresholdSum(8.),
  fIsZeroSuppressed(kFALSE),
  fLastSector(-1),
  fSecRejectRatio(.4),
  fParam(new  AliTPCParam),
  fPedestalTPC(0x0),
  fPadNoiseTPC(0x0),
  fPedestalROC(0x0),
  fPadNoiseROC(0x0),
  fCalRocArrayT0(72),
  fCalRocArrayT0Err(72),
  fCalRocArrayQ(72),
  fCalRocArrayRMS(72),
  fCalRocArrayOutliers(72),
  fHistoQArray(72),
  fHistoT0Array(72),
  fHistoRMSArray(72),
  fMeanT0rms(0),
  fMeanQrms(0),
  fMeanRMSrms(0),
  fHistoTmean(72),
  fParamArrayEventPol1(72),
  fParamArrayEventPol2(72),
  fTMeanArrayEvent(72),
  fQMeanArrayEvent(72),
  fVEventTime(1000),
  fVEventNumber(1000),
  fVTime0SideA(1000),
  fVTime0SideC(1000),
  fEventId(-1),
  fOldRunNumber(0),
  fPadTimesArrayEvent(72),
  fPadQArrayEvent(72),
  fPadRMSArrayEvent(72),
  fPadPedestalArrayEvent(72),
  fCurrentChannel(-1),
  fCurrentSector(-1),
  fCurrentRow(-1),
  fMaxPadSignal(-1),
  fMaxTimeBin(-1),
//   fPadSignal(1024),
  fPadPedestal(0),
  fPadNoise(0),
  fVTime0Offset(72),
  fVTime0OffsetCounter(72),
  fVMeanQ(72),
  fVMeanQCounter(72),
  fCurrentCETimeRef(0),
  fProcessOld(kTRUE),
  fProcessNew(kFALSE),
  fAnalyseNew(kTRUE),
  fHnDrift(0x0),
  fArrHnDrift(100),
  fTimeBursts(100),
  fArrFitGraphs(0x0)
{
  //
  // constructor which uses a tmap as input to set some specific parameters
  //
  SetNameTitle("AliTPCCalibCE","AliTPCCalibCE");
  fFirstTimeBin=650;
  fLastTimeBin=1030;
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
  if (config->GetValue("PeakDetMinus")) fPeakDetMinus = ((TObjString*)config->GetValue("PeakDetMinus"))->GetString().Atoi();
  if (config->GetValue("PeakDetPlus")) fPeakDetPlus = ((TObjString*)config->GetValue("PeakDetPlus"))->GetString().Atoi();
  if (config->GetValue("PeakIntMinus")) fPeakIntMinus = ((TObjString*)config->GetValue("PeakIntMinus"))->GetString().Atoi();
  if (config->GetValue("PeakIntPlus")) fPeakIntPlus = ((TObjString*)config->GetValue("PeakIntPlus"))->GetString().Atoi();
  if (config->GetValue("NoiseThresholdMax")) fNoiseThresholdMax = ((TObjString*)config->GetValue("NoiseThresholdMax"))->GetString().Atof();
  if (config->GetValue("NoiseThresholdSum")) fNoiseThresholdSum = ((TObjString*)config->GetValue("NoiseThresholdSum"))->GetString().Atof();
  if (config->GetValue("IsZeroSuppressed")) fIsZeroSuppressed = (Bool_t)((TObjString*)config->GetValue("IsZeroSuppressed"))->GetString().Atoi();
  if (config->GetValue("UseL1Phase")) fUseL1Phase = (Bool_t)((TObjString*)config->GetValue("UseL1Phase"))->GetString().Atoi();
  if (config->GetValue("SecRejectRatio")) fSecRejectRatio = ((TObjString*)config->GetValue("SecRejectRatio"))->GetString().Atof();

  if (config->GetValue("ProcessOld")) fProcessOld = (Bool_t)((TObjString*)config->GetValue("ProcessOld"))->GetString().Atoi();
  if (config->GetValue("ProcessNew")) fProcessNew = (Bool_t)((TObjString*)config->GetValue("ProcessNew"))->GetString().Atoi();
  if (config->GetValue("AnalyseNew")) fAnalyseNew = (Bool_t)((TObjString*)config->GetValue("AnalyseNew"))->GetString().Atoi();
  
  for (Int_t i=0;i<1024;++i) fPadSignal[i]=0;
  for (Int_t i=0;i<5;++i){
    fPeaks[i]=0;
    fPeakWidths[i]=0;
  }
  
  fParam->Update();
  for (Int_t i=0; i<100; ++i) fBinsLastAna[i]=0;
}

//_____________________________________________________________________
AliTPCCalibCE& AliTPCCalibCE::operator = (const  AliTPCCalibCE &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCCalibCE(source);

  return *this;
}
//_____________________________________________________________________
AliTPCCalibCE::~AliTPCCalibCE()
{
  //
  // destructor
  //
  
  fCalRocArrayT0.Delete();
  fCalRocArrayT0Err.Delete();
  fCalRocArrayQ.Delete();
  fCalRocArrayRMS.Delete();
  fCalRocArrayOutliers.Delete();
  
  fHistoQArray.Delete();
  fHistoT0Array.Delete();
  fHistoRMSArray.Delete();
  
  fHistoTmean.Delete();
  
  fParamArrayEventPol1.Delete();
  fParamArrayEventPol2.Delete();
  fTMeanArrayEvent.Delete();
  fQMeanArrayEvent.Delete();
  
  fPadTimesArrayEvent.Delete();
  fPadQArrayEvent.Delete();
  fPadRMSArrayEvent.Delete();
  fPadPedestalArrayEvent.Delete();
  
  fArrHnDrift.SetOwner();
  fArrHnDrift.Delete();
  
  if (fArrFitGraphs){
    fArrFitGraphs->SetOwner();
    delete fArrFitGraphs;
  }
}
//_____________________________________________________________________
Int_t AliTPCCalibCE::Update(const Int_t icsector,
				const Int_t icRow,
				const Int_t icPad,
				const Int_t icTimeBin,
				const Float_t csignal)
{
  //
  // Signal filling methode on the fly pedestal and Time offset correction if necessary.
  // no extra analysis necessary. Assumes knowledge of the signal shape!
  // assumes that it is looped over consecutive time bins of one pad
  //

  if (!fProcessOld) return 0;
  //temp

  if (icRow<0) return 0;
  if (icPad<0) return 0;
  if (icTimeBin<0) return 0;
  if ( (icTimeBin>fLastTimeBin) || (icTimeBin<fFirstTimeBin)   ) return 0;

  Int_t iChannel  = fROC->GetRowIndexes(icsector)[icRow]+icPad; //  global pad position in sector

  //init first pad and sector in this event
  if ( fCurrentChannel == -1 ) {
    fLastSector=-1;
    fCurrentChannel = iChannel;
    fCurrentSector  = icsector;
    fCurrentRow     = icRow;
  }

  //process last pad if we change to a new one
  if ( iChannel != fCurrentChannel ){
    ProcessPad();
    fLastSector=fCurrentSector;
    fCurrentChannel = iChannel;
    fCurrentSector  = icsector;
    fCurrentRow     = icRow;
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
void AliTPCCalibCE::ProcessBunch(const Int_t sector, const Int_t row, const Int_t pad,
                  const Int_t length, const UInt_t startTimeBin, const UShort_t* signal)
{
  //
  // new filling method to fill the THnSparse histogram
  //

  //only in new processing mode
  if (!fProcessNew) return;
  //don't use the IROCs
  if (sector<36) return;
  //only bunches with reasonable length
  if (length<3||length>10) return;
  
  UShort_t  timeBin    = (UShort_t)startTimeBin;
  //skip first laser layer
  if (timeBin<200) return;
  Double_t timeBurst=SetBurstHnDrift();
  
  //after 1 event setup peak ranges
  if (fNevents==1&&fPeaks[4]==0) {
    // set time range
    fHnDrift->GetAxis(4)->SetRangeUser(timeBurst-2*60,timeBurst+2*60);
    FindLaserLayers();
    // set time range
    fHnDrift->GetAxis(4)->SetRangeUser(fHnDrift->GetAxis(4)->GetXmin(),fHnDrift->GetAxis(4)->GetXmax());
  }
  
  // After the first event only fill every 5th  bin in a row with the CE information
  Int_t padFill=pad;
  if (fPeaks[4]>100&&TMath::Abs((Short_t)fPeaks[4]-(Short_t)timeBin)<(Short_t)fPeakWidths[4]){
    Int_t mod=5;
    Int_t n=pad/mod;
    padFill=mod*n+mod/2;
  }

  //noise removal
  if (!IsPeakInRange(timeBin+length/2)) return;
  
  Double_t x[kHnBinsDV]={(Double_t)sector,(Double_t)row,
      (Double_t)padFill,(Double_t)timeBin,timeBurst};
  
  for (Int_t iTimeBin = 0; iTimeBin<length; iTimeBin++){
    Float_t sig=(Float_t)signal[iTimeBin];
    x[3]=timeBin;
    fHnDrift->Fill(x,sig);
    --timeBin;
  }
}
//_____________________________________________________________________
void AliTPCCalibCE::FindLaserLayers()
{
  //
  // Find the laser layer positoins
  //


  //find CE signal position and width
  TH1D *hproj=fHnDrift->Projection(3);
  hproj->GetXaxis()->SetRangeUser(700,1030);
  Int_t maxbin=hproj->GetMaximumBin();
  Double_t binc=hproj->GetBinCenter(maxbin);
  hproj->GetXaxis()->SetRangeUser(binc-10,binc+10);
  
  fPeaks[4]=(UShort_t)TMath::Nint(hproj->GetMean());
//   fPeakWidths[4]=(UShort_t)TMath::Nint(4*hproj->GetRMS()+.5);
  fPeakWidths[4]=(UShort_t)TMath::Nint(10.);
  
  //other peaks
  Int_t timepos=fPeaks[4]-2*fPeakWidths[4];
  Int_t width=100;
  
  for (Int_t i=3; i>=0; --i){
    hproj->GetXaxis()->SetRangeUser(timepos-width,timepos);
    fPeaks[i]=(UShort_t)TMath::Nint(hproj->GetMean());
    fPeakWidths[i]=(UShort_t)TMath::Nint(4*hproj->GetRMS()+.5);
    width=250;
    timepos=fPeaks[i]-width/2;
  }
  
  //check width and reset peak if >100
  for (Int_t i=0; i<5; ++i){
    if (fPeakWidths[i]>100) {
      fPeaks[i]=0;
      fPeakWidths[i]=0;
    }
  }

  delete hproj;
}

//_____________________________________________________________________
void AliTPCCalibCE::FindPedestal(Float_t part)
{
  //
    // find pedestal and noise for the current pad. Use either database or
    // truncated mean with part*100%
  //
  Bool_t noPedestal = kTRUE;

    //use pedestal database if set
  if (fPedestalTPC&&fPadNoiseTPC){
        //only load new pedestals if the sector has changed
    if ( fCurrentSector!=fLastSector ){
      fPedestalROC = fPedestalTPC->GetCalROC(fCurrentSector);
      fPadNoiseROC = fPadNoiseTPC->GetCalROC(fCurrentSector);
    }

    if ( fPedestalROC&&fPadNoiseROC ){
      fPadPedestal = fPedestalROC->GetValue(fCurrentChannel)*(Float_t)(!fIsZeroSuppressed);
      fPadNoise    = fPadNoiseROC->GetValue(fCurrentChannel);
      noPedestal   = kFALSE;
    }

  }

    //if we are not running with pedestal database, or for the current sector there is no information
    //available, calculate the pedestal and noise on the fly
  if ( noPedestal ) {
    fPadPedestal = 0;
    fPadNoise    = 0;
    if ( fIsZeroSuppressed ) return;
    const Int_t kPedMax = 100;  //maximum pedestal value
    Float_t  max    =  0;
    Float_t  maxPos =  0;
    Int_t    median =  -1;
    Int_t    count0 =  0;
    Int_t    count1 =  0;
    //
    Float_t padSignal=0;
    //
    UShort_t histo[kPedMax];
    memset(histo,0,kPedMax*sizeof(UShort_t));

        //fill pedestal histogram
    for (Int_t i=fFirstTimeBin; i<=fLastTimeBin; ++i){
      padSignal = fPadSignal[i];
      if (padSignal<=0) continue;
      if (padSignal>max && i>10) {
        max = padSignal;
        maxPos = i;
      }
      if (padSignal>kPedMax-1) continue;
      histo[int(padSignal+0.5)]++;
      count0++;
    }
	//find median
    for (Int_t i=1; i<kPedMax; ++i){
      if (count1<count0*0.5) median=i;
      count1+=histo[i];
    }
	// truncated mean
    //
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
    if ( count > 0 ) {
      mean/=count;
      rms    = TMath::Sqrt(TMath::Abs(rms/count-mean*mean));
      fPadPedestal = mean;
      fPadNoise    = rms;
    }
  }
}
//_____________________________________________________________________
void AliTPCCalibCE::UpdateCETimeRef()
{
  // Find the time reference of the last valid CE signal in sector
  // for irocs of the A-Side the reference of the corresponging OROC is returned
  // the reason are the non reflective bands on the A-Side, which make the reference very uncertain
  if ( fLastSector == fCurrentSector ) return;
  Int_t sector=fCurrentSector; 
  if ( sector < 18 ) sector+=36;
  fCurrentCETimeRef=0;
  TVectorF *vtRef = GetTMeanEvents(sector);
  if ( !vtRef ) return; 
  Int_t vtRefSize= vtRef->GetNrows();
  if ( vtRefSize < fNevents+1 ) vtRef->ResizeTo(vtRefSize+100);
  else vtRefSize=fNevents; 
  while ( (*vtRef)[vtRefSize]==0 && vtRefSize>=0 ) --vtRefSize;
  fCurrentCETimeRef=(*vtRef)[vtRefSize];
  AliDebug(3,Form("Sector: %02d - T0 ref: %.2f",fCurrentSector,fCurrentCETimeRef)); 
} 
//_____________________________________________________________________
void AliTPCCalibCE::FindCESignal(TVectorD &param, Float_t &qSum, const TVectorF maxima)
{
  //
    //  Find position, signal width and height of the CE signal (last signal)
    //  param[0] = Qmax, param[1] = mean time, param[2] = rms;
    //  maxima: array of local maxima of the pad signal use the one closest to the mean CE position
  //

  Float_t ceQmax  =0, ceQsum=0, ceTime=0, ceRMS=0;
  Int_t   cemaxpos       = 0;
  Float_t ceSumThreshold = fNoiseThresholdSum*fPadNoise;  // threshold for the signal sum
  const Int_t    kCemin  = fPeakIntMinus;             // range for the analysis of the ce signal +- channels from the peak
  const Int_t    kCemax  = fPeakIntPlus;

  Float_t minDist  = 25;  //initial minimum distance betweek roc mean ce signal and pad ce signal

    // find maximum closest to the sector mean from the last event
  for ( Int_t imax=0; imax<maxima.GetNrows(); ++imax){
        // get sector mean of last event
    Float_t tmean = fCurrentCETimeRef;
    if ( TMath::Abs( tmean-maxima[imax] ) < minDist ) {
      minDist  = tmean-maxima[imax];
      cemaxpos = (Int_t)maxima[imax];
    }
  }
//   printf("L1 phase TB: %f\n",GetL1PhaseTB());
  if (cemaxpos!=0){
    ceQmax = fPadSignal[cemaxpos]-fPadPedestal;
    for (Int_t i=cemaxpos-kCemin; i<=cemaxpos+kCemax; ++i){
      if ( (i>fFirstTimeBin) && (i<fLastTimeBin) ){
        Float_t signal = fPadSignal[i]-fPadPedestal;
        if (signal>0) {
          ceTime+=signal*(i+0.5);
          ceRMS +=signal*(i+0.5)*(i+0.5);
          ceQsum+=signal;
        }
      }
    }
  }
  if (ceQmax&&ceQsum>ceSumThreshold) {
    ceTime/=ceQsum;
    ceRMS  = TMath::Sqrt(TMath::Abs(ceRMS/ceQsum-ceTime*ceTime));
    ceTime-=GetL1PhaseTB();
    fVTime0Offset.GetMatrixArray()[fCurrentSector]+=ceTime;   // mean time for each sector
    fVTime0OffsetCounter.GetMatrixArray()[fCurrentSector]++;

  //Normalise Q to the 'cell-size': The wire density is the same in the IROC and OROC, therefore the
  //                                the pick-up signal should scale with the pad area. In addition
  //                                the signal should decrease with the wire distance (4mm in IROC, 6mm in OROC),
  //                                ratio 2/3. The pad area we express in cm2. We normalise the signal
  //                                to the OROC signal (factor 2/3 for the IROCs).  
    Float_t norm = fParam->GetPadPitchWidth(fCurrentSector)*fParam->GetPadPitchLength(fCurrentSector,fCurrentRow);
    if ( fCurrentSector<fParam->GetNInnerSector() ) norm*=3./2.;

    ceQsum/=norm;
    fVMeanQ.GetMatrixArray()[fCurrentSector]+=ceQsum;
    fVMeanQCounter.GetMatrixArray()[fCurrentSector]++;
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
Bool_t AliTPCCalibCE::IsPeak(Int_t pos, Int_t tminus, Int_t tplus) const
{
  //
  // Check if 'pos' is a Maximum. Consider 'tminus' timebins before
  // and 'tplus' timebins after 'pos'
  //
  if ( (pos-tminus)<fFirstTimeBin || (pos+tplus)>fLastTimeBin ) return kFALSE;
  for (Int_t iTime = pos; iTime>pos-tminus; --iTime)
    if ( fPadSignal[iTime-1] >= fPadSignal[iTime] ) return kFALSE;
  for (Int_t iTime = pos, iTime2=pos; iTime<pos+tplus; ++iTime, ++iTime2){
    if ( (iTime==pos) && (fPadSignal[iTime+1]==fPadSignal[iTime]) ) // allow two timebins with same adc value
      ++iTime2;
    if ( fPadSignal[iTime2+1] >= fPadSignal[iTime2] ) return kFALSE;
  }
  return kTRUE;
}
//_____________________________________________________________________
void AliTPCCalibCE::FindLocalMaxima(TVectorF &maxima)
{
  //
    // Find local maxima on the pad signal and Histogram them
  //
  Float_t ceThreshold = fNoiseThresholdMax*TMath::Max(fPadNoise,Float_t(1.));  // threshold for the signal
  Int_t   count       = 0;
  
  for (Int_t i=fLastTimeBin-fPeakDetPlus+1; i>=fFirstTimeBin+fPeakDetMinus; --i){
    if ( (fPadSignal[i]-fPadPedestal)<ceThreshold ) continue;
    if (IsPeak(i,fPeakDetMinus,fPeakDetPlus) ){
      if (count<maxima.GetNrows()){
        maxima.GetMatrixArray()[count++]=i;
        GetHistoTmean(fCurrentSector,kTRUE)->Fill(i);
        i-=(fPeakDetMinus+fPeakDetPlus-1); // next peak cannot be at bin  fPeakDetMinus+fPeakDetPlus-1
      }
    }
  }
}
//_____________________________________________________________________
void AliTPCCalibCE::ProcessPad()
{
  //
  //  Process data of current pad
  //
  FindPedestal();
  
  TVectorF maxima(15);     // the expected maximum number of maxima in the complete TPC should be 8 laser beam layers
                             // + central electrode and possibly post peaks from the CE signal
                             // however if we are on a high noise pad a lot more peaks due to the noise might occur
  FindLocalMaxima(maxima);
  if ( (fNevents == 0) || (fOldRunNumber!=fRunNumber) ) return;  // return because we don't have Time0 info for the CE yet
  
  UpdateCETimeRef();                       // update the time refenrence for the current sector
  if ( fCurrentCETimeRef<1e-30 ) return;      //return if we don't have time 0 info, eg if only one side has laser
  TVectorD param(3);
  Float_t  qSum;
  FindCESignal(param, qSum, maxima);
  
  Double_t meanT  = param[1];
  Double_t sigmaT = param[2];
  
    //Fill Event T0 counter
  (*GetPadTimesEvent(fCurrentSector,kTRUE)).GetMatrixArray()[fCurrentChannel] = meanT;
  
    //Fill Q histogram
  GetHistoQ(fCurrentSector,kTRUE)->Fill( TMath::Sqrt(qSum), fCurrentChannel );
  
    //Fill RMS histogram
  GetHistoRMS(fCurrentSector,kTRUE)->Fill( sigmaT, fCurrentChannel );
  
  
    //Fill debugging info
  if ( GetStreamLevel()>0 ){
    (*GetPadPedestalEvent(fCurrentSector,kTRUE)).GetMatrixArray()[fCurrentChannel]=fPadPedestal;
    (*GetPadRMSEvent(fCurrentSector,kTRUE)).GetMatrixArray()[fCurrentChannel]=sigmaT;
    (*GetPadQEvent(fCurrentSector,kTRUE)).GetMatrixArray()[fCurrentChannel]=qSum;
  }
  
  ResetPad();
}
//_____________________________________________________________________
void AliTPCCalibCE::EndEvent()
{
  //  Process data of current pad
  //  The Functions 'SetTimeStamp' and 'SetRunNumber'  should be called
  //  before the EndEvent function to set the event timestamp and number!!!
  //  This is automatically done if the ProcessEvent(AliRawReader *rawReader)
  //  function was called
  if (!fProcessOld) {
    if (fProcessNew) ++fNevents;
    return;
  }
  
  //check if last pad has allready been processed, if not do so
  if ( fMaxTimeBin>-1 ) ProcessPad();

  AliDebug(3, Form("EndEvent() - Start; Event: %05d", fNevents));

  TVectorD param(3);
  TMatrixD dummy(3,3);
//    TVectorF vMeanTime(72);
//    TVectorF vMeanQ(72);
  AliTPCCalROC *calIroc=new AliTPCCalROC(0);
  AliTPCCalROC *calOroc=new AliTPCCalROC(36);

  //find mean time0 offset for side A and C
  //use only orocs due to the better statistics
  Double_t time0Side[2];       //time0 for side A:0 and C:1
  Double_t time0SideCount[2];  //time0 counter for side A:0 and C:1
  time0Side[0]=0;time0Side[1]=0;time0SideCount[0]=0;time0SideCount[1]=0;
  for ( Int_t iSec = 36; iSec<72; ++iSec ){
    time0Side[(iSec/18)%2] += fVTime0Offset.GetMatrixArray()[iSec];
    time0SideCount[(iSec/18)%2] += fVTime0OffsetCounter.GetMatrixArray()[iSec];
  }
  if ( time0SideCount[0] >0  )
    time0Side[0]/=time0SideCount[0];
  if ( time0SideCount[1] >0 )
    time0Side[1]/=time0SideCount[1];
    // end find time0 offset
  AliDebug(3,Form("time0Side/time0SideCount: A=%.2f/%.2f, C=%.2f/%.2f",time0Side[0],time0SideCount[0],time0Side[1],time0SideCount[1]));
  Int_t nSecMeanT=0;
  //loop over all ROCs, fill CE Time histogram corrected for the mean Time0 of each ROC
  for ( Int_t iSec = 0; iSec<72; ++iSec ){
    AliDebug(4,Form("Processing sector '%02d'\n",iSec));
    //find median and then calculate the mean around it
    TH1S *hMeanT    = GetHistoTmean(iSec); //histogram with local maxima position information
    if ( !hMeanT ) continue;
    //continue if not enough data is filled in the meanT histogram. This is the case if we do not have a laser event.
    if ( hMeanT->GetEffectiveEntries() < fROC->GetNChannels(iSec)*fSecRejectRatio ){
      hMeanT->Reset();
      AliDebug(3,Form("Skipping sec. '%02d': Not enough statistics\n",iSec));
      continue;
    }
    
    Double_t entries = hMeanT->GetEffectiveEntries();
    Double_t sum     = 0;
    Short_t *arr     = hMeanT->GetArray()+1;
    Int_t ibin=0;
    for ( ibin=0; ibin<hMeanT->GetNbinsX(); ++ibin){
      sum+=arr[ibin];
      if ( sum>=(entries/2.) ) break;
    }
    Int_t delta = 4;
    Int_t firstBin = fFirstTimeBin+ibin-delta;
    Int_t lastBin  = fFirstTimeBin+ibin+delta;
    if ( firstBin<fFirstTimeBin ) firstBin=fFirstTimeBin;
    if ( lastBin>fLastTimeBin   ) lastBin =fLastTimeBin;
    Float_t median =AliMathBase::GetCOG(arr+ibin-delta,2*delta,firstBin,lastBin);
    
	// check boundaries for ebye info of mean time
    TVectorF *vMeanTime=GetTMeanEvents(iSec,kTRUE);
    Int_t vSize=vMeanTime->GetNrows();
    if ( vSize < fNevents+1 ){
      vMeanTime->ResizeTo(vSize+100);
    }

    // store mean time for the readout sides
    vSize=fVTime0SideA.GetNrows();
    if ( vSize < fNevents+1 ){
      fVTime0SideA.ResizeTo(vSize+100);
      fVTime0SideC.ResizeTo(vSize+100);
    }
    fVTime0SideA.GetMatrixArray()[fNevents]=time0Side[0];
    fVTime0SideC.GetMatrixArray()[fNevents]=time0Side[1];
    
    vMeanTime->GetMatrixArray()[fNevents]=median;
    nSecMeanT++;
    // end find median
    
    TVectorF *vTimes = GetPadTimesEvent(iSec);
    if ( !vTimes ) continue;                     //continue if no time information for this sector is available
    
    AliTPCCalROC calIrocOutliers(0);
    AliTPCCalROC calOrocOutliers(36);
    
    // calculate mean Q of the sector
    TVectorF *vMeanQ=GetQMeanEvents(iSec,kTRUE);
    vSize=vMeanQ->GetNrows();
    if ( vSize < fNevents+1 ){
      vMeanQ->ResizeTo(vSize+100);
    }   
    Float_t meanQ = 0;
    if ( fVMeanQCounter.GetMatrixArray()[iSec]>0 ) meanQ=fVMeanQ.GetMatrixArray()[iSec]/fVMeanQCounter.GetMatrixArray()[iSec];
    vMeanQ->GetMatrixArray()[fNevents]=meanQ;
   
    for ( UInt_t iChannel=0; iChannel<fROC->GetNChannels(iSec); ++iChannel ){
      Float_t time  = (*vTimes).GetMatrixArray()[iChannel];

	    //set values for temporary roc calibration class
      if ( iSec < 36 ) {
        calIroc->SetValue(iChannel, time);
        if ( TMath::Abs(time) < 1e-30 ) calIrocOutliers.SetValue(iChannel,1);

      } else {
        calOroc->SetValue(iChannel, time);
        if ( TMath::Abs(time) < 1e-30 ) calOrocOutliers.SetValue(iChannel,1);
      }

      if ( (fNevents>0) && (fOldRunNumber==fRunNumber) )
        // test that we really found the CE signal reliably 
        if ( TMath::Abs(fVTime0SideA.GetMatrixArray()[fNevents-1]-time0Side[0])<.05)
          GetHistoT0(iSec,kTRUE)->Fill( time-time0Side[(iSec/18)%2],iChannel );



	    //-------------------------------  Debug start  ------------------------------
      if ( GetStreamLevel()>0 ){
        TTreeSRedirector *streamer=GetDebugStreamer();
        if (streamer){
          Int_t row=0;
          Int_t pad=0;
          Int_t padc=0;
          
          Float_t q   = (*GetPadQEvent(iSec))[iChannel];
          Float_t rms = (*GetPadRMSEvent(iSec))[iChannel];
          
          UInt_t channel=iChannel;
          Int_t sector=iSec;
          
          while ( channel > (fROC->GetRowIndexes(sector)[row]+fROC->GetNPads(sector,row)-1) ) row++;
          pad = channel-fROC->GetRowIndexes(sector)[row];
          padc = pad-(fROC->GetNPads(sector,row)/2);
          
//		TH1F *h1 = new TH1F(Form("hSignalD%d.%d.%d",sector,row,pad),
//				    Form("hSignalD%d.%d.%d",sector,row,pad),
//				    fLastTimeBin-fFirstTimeBin,
//				    fFirstTimeBin,fLastTimeBin);
//		h1->SetDirectory(0);
        //
//		for (Int_t i=fFirstTimeBin; i<fLastTimeBin+1; ++i)
//		    h1->Fill(i,fPadSignal(i));
          
          Double_t t0Sec = 0;
          if (fVTime0OffsetCounter.GetMatrixArray()[iSec]>0)
            t0Sec = fVTime0Offset.GetMatrixArray()[iSec]/fVTime0OffsetCounter.GetMatrixArray()[iSec];
          Double_t t0Side = time0Side[(iSec/18)%2];
          (*streamer) << "DataPad" <<
            "Event=" << fNevents <<
            "RunNumber=" << fRunNumber <<
            "TimeStamp="   << fTimeStamp <<
            "Sector="<< sector <<
            "Row="   << row<<
            "Pad="   << pad <<
            "PadC="  << padc <<
            "PadSec="<< channel <<
            "Time0Sec="  << t0Sec <<
            "Time0Side=" << t0Side <<
            "Time="  << time <<
            "RMS="   << rms <<
            "Sum="   << q <<
            "MeanQ=" << meanQ <<
        //		    "hist.=" << h1 <<
            "\n";
          
    //		delete h1;
        }
      }
      //-----------------------------  Debug end  ------------------------------
    }// end channel loop


    //do fitting now only in debug mode
    if (GetDebugLevel()>0){
      TVectorD paramPol1(3);
      TVectorD paramPol2(6);
      TMatrixD matPol1(3,3);
      TMatrixD matPol2(6,6);
      Float_t  chi2Pol1=0;
      Float_t  chi2Pol2=0;
      
      if ( (fNevents>0) && (fOldRunNumber==fRunNumber) ){
        if ( iSec < 36 ){
          calIroc->GlobalFit(&calIrocOutliers,0,paramPol1,matPol1,chi2Pol1,0);
          calIroc->GlobalFit(&calIrocOutliers,0,paramPol2,matPol2,chi2Pol2,1);
        } else {
          calOroc->GlobalFit(&calOrocOutliers,0,paramPol1,matPol1,chi2Pol1,0);
          calOroc->GlobalFit(&calOrocOutliers,0,paramPol2,matPol2,chi2Pol2,1);
        }
        
        GetParamArrayPol1(iSec,kTRUE)->AddAtAndExpand(new TVectorD(paramPol1), fNevents);
        GetParamArrayPol2(iSec,kTRUE)->AddAtAndExpand(new TVectorD(paramPol2), fNevents);
      }
      
  //-------------------------------  Debug start  ------------------------------
      if ( GetStreamLevel()>0 ){
        TTreeSRedirector *streamer=GetDebugStreamer();
        if ( streamer ) {
          (*streamer) << "DataRoc" <<
//    "Event=" << fEvent <<
            "RunNumber=" << fRunNumber <<
            "TimeStamp="   << fTimeStamp <<
            "Sector="<< iSec <<
            "hMeanT.=" << hMeanT <<
            "median=" << median <<
            "paramPol1.=" << &paramPol1 <<
            "paramPol2.=" << &paramPol2 <<
            "matPol1.="   << &matPol1 <<
            "matPol2.="   << &matPol2 <<
            "chi2Pol1="   << chi2Pol1 <<
            "chi2Pol2="   << chi2Pol2 <<
            "\n";
        }
      }
    }
	//-------------------------------  Debug end  ------------------------------
    hMeanT->Reset();
  }// end sector loop
    //return if no sector has a valid mean time
  if ( nSecMeanT == 0 ) return;
    
    
//    fTMeanArrayEvent.AddAtAndExpand(new TVectorF(vMeanTime),fNevents);
//    fQMeanArrayEvent.AddAtAndExpand(new TVectorF(vMeanQ),fNevents);
  if ( fVEventTime.GetNrows() < fNevents+1 ) {
    fVEventTime.ResizeTo((Int_t)(fVEventTime.GetNrows()+100));
    fVEventNumber.ResizeTo((Int_t)(fVEventNumber.GetNrows()+100));
  }
  fVEventTime.GetMatrixArray()[fNevents] = fTimeStamp;
  fVEventNumber.GetMatrixArray()[fNevents] = fEventId;

  ++fNevents;
  fOldRunNumber = fRunNumber;

  delete calIroc;
  delete calOroc;
  AliDebug(3, Form("EndEvent() - End; Event: %05d", fNevents));
}
//_____________________________________________________________________
TH2S* AliTPCCalibCE::GetHisto(Int_t sector, TObjArray *arr,
				  Int_t nbinsY, Float_t ymin, Float_t ymax,
				  const Char_t *type, Bool_t force)
{
    //
    // return pointer to TH2S histogram of 'type'
    // if force is true create a new histogram if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(sector) )
      return (TH2S*)arr->UncheckedAt(sector);

    // if we are forced and histogram doesn't exist yet create it
    // new histogram with Q calib information. One value for each pad!
    TH2S* hist = new TH2S(Form("hCalib%s%.2d",type,sector),Form("%s calibration histogram sector %.2d",type,sector),
			  nbinsY, ymin, ymax,
			  fROC->GetNChannels(sector),0,fROC->GetNChannels(sector));
    hist->SetDirectory(0);
    arr->AddAt(hist,sector);
    return hist;
}
//_____________________________________________________________________
TH2S* AliTPCCalibCE::GetHistoT0(Int_t sector, Bool_t force)
{
    //
    // return pointer to T0 histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fHistoT0Array;
    return GetHisto(sector, arr, fNbinsT0, fXminT0, fXmaxT0, "T0", force);
}
//_____________________________________________________________________
TH2S* AliTPCCalibCE::GetHistoQ(Int_t sector, Bool_t force)
{
    //
    // return pointer to Q histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fHistoQArray;
    return GetHisto(sector, arr, fNbinsQ, fXminQ, fXmaxQ, "Q", force);
}
//_____________________________________________________________________
TH2S* AliTPCCalibCE::GetHistoRMS(Int_t sector, Bool_t force)
{
    //
    // return pointer to Q histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fHistoRMSArray;
    return GetHisto(sector, arr, fNbinsRMS, fXminRMS, fXmaxRMS, "RMS", force);
}
//_____________________________________________________________________
TH1S* AliTPCCalibCE::GetHisto(Int_t sector, TObjArray *arr,
			      const Char_t *type, Bool_t force)
{
    //
    // return pointer to TH1S histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(sector) )
      return (TH1S*)arr->UncheckedAt(sector);

    // if we are forced and histogram doesn't yes exist create it
    // new histogram with calib information. One value for each pad!
    TH1S* hist = new TH1S(Form("hCalib%s%.2d",type,sector),Form("%s calibration histogram sector %.2d",type,sector),
			  fLastTimeBin-fFirstTimeBin,fFirstTimeBin,fLastTimeBin);
    hist->SetDirectory(0);
    arr->AddAt(hist,sector);
    return hist;
}
//_____________________________________________________________________
TH1S* AliTPCCalibCE::GetHistoTmean(Int_t sector, Bool_t force)
{
    //
    // return pointer to Q histogram
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fHistoTmean;
    return GetHisto(sector, arr, "LastTmean", force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibCE::GetVectSector(Int_t sector, TObjArray *arr, UInt_t size, Bool_t force) const
{
    //
    // return pointer to Pad Info from 'arr' for the current event and sector
    // if force is true create it if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(sector) )
	return (TVectorF*)arr->UncheckedAt(sector);

    TVectorF *vect = new TVectorF(size);
    arr->AddAt(vect,sector);
    return vect;
}
//_____________________________________________________________________
TVectorF* AliTPCCalibCE::GetPadTimesEvent(Int_t sector, Bool_t force)
{
    //
    // return pointer to Pad Times Array for the current event and sector
    // if force is true create it if it doesn't exist allready
    //
    TObjArray *arr = &fPadTimesArrayEvent;
    return GetVectSector(sector,arr,fROC->GetNChannels(sector),force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibCE::GetPadQEvent(Int_t sector, Bool_t force)
{
    //
    // return pointer to Pad Q Array for the current event and sector
    // if force is true create it if it doesn't exist allready
    // for debugging purposes only
    //

    TObjArray *arr = &fPadQArrayEvent;
    return GetVectSector(sector,arr,fROC->GetNChannels(sector),force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibCE::GetPadRMSEvent(Int_t sector, Bool_t force)
{
    //
    // return pointer to Pad RMS Array for the current event and sector
    // if force is true create it if it doesn't exist allready
    // for debugging purposes only
    //
    TObjArray *arr = &fPadRMSArrayEvent;
    return GetVectSector(sector,arr,fROC->GetNChannels(sector),force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibCE::GetPadPedestalEvent(Int_t sector, Bool_t force)
{
    //
    // return pointer to Pad RMS Array for the current event and sector
    // if force is true create it if it doesn't exist allready
    // for debugging purposes only
    //
    TObjArray *arr = &fPadPedestalArrayEvent;
    return GetVectSector(sector,arr,fROC->GetNChannels(sector),force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibCE::GetTMeanEvents(Int_t sector, Bool_t force)
{
    //
    // return pointer to the EbyE info of the mean arrival time for 'sector'
    // if force is true create it if it doesn't exist allready
    //
    TObjArray *arr = &fTMeanArrayEvent;
    return GetVectSector(sector,arr,100,force);
}
//_____________________________________________________________________
TVectorF* AliTPCCalibCE::GetQMeanEvents(Int_t sector, Bool_t force)
{
    //
    // return pointer to the EbyE info of the mean arrival time for 'sector'
    // if force is true create it if it doesn't exist allready
    //
    TObjArray *arr = &fQMeanArrayEvent;
    return GetVectSector(sector,arr,100,force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibCE::GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force) const
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
AliTPCCalROC* AliTPCCalibCE::GetCalRocT0(Int_t sector, Bool_t force)
{
    //
    // return pointer to Time 0 ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayT0;
    return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibCE::GetCalRocT0Err(Int_t sector, Bool_t force)
{
    //
    // return pointer to the error of Time 0 ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayT0Err;
    return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibCE::GetCalRocQ(Int_t sector, Bool_t force)
{
    //
    // return pointer to T0 ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayQ;
    return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibCE::GetCalRocRMS(Int_t sector, Bool_t force)
{
    //
    // return pointer to signal width ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayRMS;
    return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
AliTPCCalROC* AliTPCCalibCE::GetCalRocOutliers(Int_t sector, Bool_t force)
{
    //
    // return pointer to Outliers
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fCalRocArrayOutliers;
    return GetCalRoc(sector, arr, force);
}
//_____________________________________________________________________
TObjArray* AliTPCCalibCE::GetParamArray(Int_t sector, TObjArray* arr, Bool_t force) const
{
    //
    // return pointer to TObjArray of fit parameters
    // if force is true create a new histogram if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(sector) )
	return (TObjArray*)arr->UncheckedAt(sector);

    // if we are forced and array doesn't yes exist create it

    // new TObjArray for parameters
    TObjArray *newArr = new TObjArray;
    arr->AddAt(newArr,sector);
    return newArr;
}
//_____________________________________________________________________
TObjArray* AliTPCCalibCE::GetParamArrayPol1(Int_t sector, Bool_t force)
{
    //
    // return pointer to TObjArray of fit parameters from plane fit
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fParamArrayEventPol1;
    return GetParamArray(sector, arr, force);
}
//_____________________________________________________________________
TObjArray* AliTPCCalibCE::GetParamArrayPol2(Int_t sector, Bool_t force)
{
    //
    // return pointer to TObjArray of fit parameters from parabola fit
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fParamArrayEventPol2;
    return GetParamArray(sector, arr, force);
}

//_____________________________________________________________________
void AliTPCCalibCE::CreateDVhist()
{
  //
  // Setup the THnSparse for the drift velocity determination
  //

  //HnSparse bins
  //roc, row, pad, timebin, timestamp
  TTimeStamp begin(2010,01,01,0,0,0);
  TTimeStamp end(2035,01,01,0,0,0);
  Int_t nbinsTime=(end.GetSec()-begin.GetSec())/60; //Minutes resolution
  
  Int_t    bins[kHnBinsDV] = { 72,  96,  140,  1030, nbinsTime};
  Double_t xmin[kHnBinsDV] = { 0.,  0.,   0.,    0., (Double_t)begin.GetSec()};
  Double_t xmax[kHnBinsDV] = {72., 96., 140., 1030., (Double_t)end.GetSec()};
  
  fHnDrift=new THnSparseI("fHnDrift","Laser digits",kHnBinsDV, bins, xmin, xmax);
  fHnDrift->GetAxis(0)->SetNameTitle("ROC","Read-out chamber number");
  fHnDrift->GetAxis(1)->SetNameTitle("Row","Row number");
  fHnDrift->GetAxis(2)->SetNameTitle("Pad","Pad number");
  fHnDrift->GetAxis(3)->SetNameTitle("Timebin","Time bin [x100ns]");
  fHnDrift->GetAxis(4)->SetNameTitle("EventTime","Event time");
  
}

//_____________________________________________________________________
void AliTPCCalibCE::ResetEvent()
{
    //
    //  Reset global counters  -- Should be called before each event is processed
    //
    fLastSector=-1;
    fCurrentSector=-1;
    fCurrentRow=-1;
    fCurrentChannel=-1;

    ResetPad();

    fPadTimesArrayEvent.Delete();
    fPadQArrayEvent.Delete();
    fPadRMSArrayEvent.Delete();
    fPadPedestalArrayEvent.Delete();

    for ( Int_t i=0; i<72; ++i ){
	fVTime0Offset.GetMatrixArray()[i]=0;
	fVTime0OffsetCounter.GetMatrixArray()[i]=0;
	fVMeanQ.GetMatrixArray()[i]=0;
        fVMeanQCounter.GetMatrixArray()[i]=0;
    }
}
//_____________________________________________________________________
void AliTPCCalibCE::ResetPad()
{
    //
    //  Reset pad infos -- Should be called after a pad has been processed
    //
    for (Int_t i=fFirstTimeBin; i<fLastTimeBin+1; ++i)
      fPadSignal[i] = 0;
    fMaxTimeBin   = -1;
    fMaxPadSignal = -1;
    fPadPedestal  = -1;
    fPadNoise     = -1;
}
//_____________________________________________________________________
void AliTPCCalibCE::Merge(AliTPCCalibCE * const ce)
{
  //
  //  Merge ce to the current AliTPCCalibCE
  //
  MergeBase(ce);
  Int_t nCEevents = ce->GetNeventsProcessed();
  
  if (fProcessOld&&ce->fProcessOld){
  //merge histograms
    for (Int_t iSec=0; iSec<72; ++iSec){
      TH2S *hRefQmerge   = ce->GetHistoQ(iSec);
      TH2S *hRefT0merge  = ce->GetHistoT0(iSec);
      TH2S *hRefRMSmerge = ce->GetHistoRMS(iSec);
      
      
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
    
    // merge time information
    
    
    for (Int_t iSec=0; iSec<72; ++iSec){
      TObjArray *arrPol1CE  = ce->GetParamArrayPol1(iSec);
      TObjArray *arrPol2CE  = ce->GetParamArrayPol2(iSec);
      TVectorF *vMeanTimeCE = ce->GetTMeanEvents(iSec);
      TVectorF *vMeanQCE    = ce->GetQMeanEvents(iSec);
      
      TObjArray *arrPol1  = 0x0;
      TObjArray *arrPol2  = 0x0;
      TVectorF *vMeanTime = 0x0;
      TVectorF *vMeanQ    = 0x0;
      
  //resize arrays
      if ( arrPol1CE && arrPol2CE ){
        arrPol1 = GetParamArrayPol1(iSec,kTRUE);
        arrPol2 = GetParamArrayPol2(iSec,kTRUE);
        arrPol1->Expand(fNevents+nCEevents);
        arrPol2->Expand(fNevents+nCEevents);
      }
      if ( vMeanTimeCE && vMeanQCE ){
        vMeanTime = GetTMeanEvents(iSec,kTRUE);
        vMeanQ    = GetQMeanEvents(iSec,kTRUE);
        vMeanTime->ResizeTo(fNevents+nCEevents);
        vMeanQ->ResizeTo(fNevents+nCEevents);
      }
      
      for (Int_t iEvent=0; iEvent<nCEevents; ++iEvent){
        if ( arrPol1CE && arrPol2CE ){
          TVectorD *paramPol1 = (TVectorD*)(arrPol1CE->UncheckedAt(iEvent));
          TVectorD *paramPol2 = (TVectorD*)(arrPol2CE->UncheckedAt(iEvent));
          if ( paramPol1 && paramPol2 ){
            GetParamArrayPol1(iSec,kTRUE)->AddAt(new TVectorD(*paramPol1), fNevents+iEvent);
            GetParamArrayPol2(iSec,kTRUE)->AddAt(new TVectorD(*paramPol2), fNevents+iEvent);
          }
        }
        if ( vMeanTimeCE && vMeanQCE ){
          vMeanTime->GetMatrixArray()[fNevents+iEvent]=vMeanTimeCE->GetMatrixArray()[iEvent];
          vMeanQ->GetMatrixArray()[fNevents+iEvent]=vMeanQCE->GetMatrixArray()[iEvent];
        }
      }
    }
    
    
    
    const TVectorD&  eventTimes  = ce->fVEventTime;
    const TVectorD&  eventIds    = ce->fVEventNumber;
    const TVectorF&  time0SideA  = ce->fVTime0SideA;
    const TVectorF&  time0SideC  = ce->fVTime0SideC;
    fVEventTime.ResizeTo(fNevents+nCEevents);
    fVEventNumber.ResizeTo(fNevents+nCEevents);
    fVTime0SideA.ResizeTo(fNevents+nCEevents);
    fVTime0SideC.ResizeTo(fNevents+nCEevents);

    for (Int_t iEvent=0; iEvent<nCEevents; ++iEvent){
      Double_t evTime     = eventTimes.GetMatrixArray()[iEvent];
      Double_t evId       = eventIds.GetMatrixArray()[iEvent];
      Float_t  t0SideA    = time0SideA.GetMatrixArray()[iEvent];
      Float_t  t0SideC    = time0SideC.GetMatrixArray()[iEvent];
      
      fVEventTime.GetMatrixArray()[fNevents+iEvent]   = evTime;
      fVEventNumber.GetMatrixArray()[fNevents+iEvent] = evId;
      fVTime0SideA.GetMatrixArray()[fNevents+iEvent]  = t0SideA;
      fVTime0SideC.GetMatrixArray()[fNevents+iEvent]  = t0SideC;
    }
  }

  if (fProcessNew&&ce->fProcessNew) {
    if (fArrHnDrift.GetEntries() != ce->fArrHnDrift.GetEntries() ){
      AliError("Number of bursts in the instances to merge are different. No merging done!");
    } else {
      for (Int_t i=0;i<fArrHnDrift.GetEntries();++i){
        THnSparseI *h=(THnSparseI*)fArrHnDrift.UncheckedAt(i);
        THnSparseI *hce=(THnSparseI*)ce->fArrHnDrift.UncheckedAt(i);
        if (h && hce) h->Add(hce);
        else AliError(Form("AliTPCCalibCE::Merge - one THnSparse missing in burst %d",i));
      }
      //TODO: What to do with fTimeBursts???
    }
  }
  
  fNevents+=nCEevents; //increase event counter
}

//_____________________________________________________________________
Long64_t AliTPCCalibCE::Merge(TCollection * const list)
{
  //
  // Merge all objects of this type in list
  //

  Long64_t nmerged=1;

  TIter next(list);
  AliTPCCalibCE *ce=0;
  TObject *o=0;

  while ( (o=next()) ){
    ce=dynamic_cast<AliTPCCalibCE*>(o);
    if (ce){
      Merge(ce);
      ++nmerged;
    }
  }

  return nmerged;
}

//_____________________________________________________________________
TGraph *AliTPCCalibCE::MakeGraphTimeCE(Int_t sector, Int_t xVariable, Int_t fitType, Int_t fitParameter)
{
  //
  // Make graph from fit parameters of pol1 fit, pol2 fit, mean arrival time or mean Q for ROC 'sector'
  // or side (-1: A-Side, -2: C-Side)
  // xVariable:    0-event time, 1-event id, 2-internal event counter
  // fitType:      0-pol1 fit, 1-pol2 fit, 2-mean time, 3-mean Q
  // fitParameter: fit parameter ( 0-2 for pol1 ([0]+[1]*x+[2]*y),
  //                               0-5 for pol2 ([0]+[1]*x+[2]*y+[3]*x*x+[4]*y*y+[5]*x*y),
  //                               not used for mean time and mean Q )
  // for an example see class description at the beginning
  //

  TVectorD *xVar = 0x0;
  TObjArray *aType = 0x0;
  Int_t npoints=0;

  // sanity checks
  if ( (sector<-2)   || (sector>71)   ) return 0x0;  //sector outside valid range
  if ( (xVariable<0) || (xVariable>2) ) return 0x0;  //invalid x-variable
  if ( (fitType<0)   || (fitType>3)   ) return 0x0;  //invalid fit type
  if ( sector>=0 && fitType==2 && !GetTMeanEvents(sector) ) return 0x0; //no mean time information available
  if ( sector>=0 && fitType==3 && !GetQMeanEvents(sector) ) return 0x0; //no mean charge information available
  if ( sector<0 && fitType!=2) return 0x0;  //for side wise information only mean time is available

  if (sector>=0){
    if ( fitType==0 ){
      if ( (fitParameter<0) || (fitParameter>2) ) return 0x0;
      aType = &fParamArrayEventPol1;
      if ( aType->At(sector)==0x0 ) return 0x0;
    }
    else if ( fitType==1 ){
      if ( (fitParameter<0) || (fitParameter>5) ) return 0x0;
      aType = &fParamArrayEventPol2;
      if ( aType->At(sector)==0x0 ) return 0x0;
    }

  }
  if ( xVariable == 0 ) xVar = &fVEventTime;
  if ( xVariable == 1 ) xVar = &fVEventNumber;
  if ( xVariable == 2 ) {
    xVar = new TVectorD(fNevents);
    for ( Int_t i=0;i<fNevents; ++i) (*xVar)[i]=i;
  }
  
  Double_t *x = new Double_t[fNevents];
  Double_t *y = new Double_t[fNevents];
  
  for (Int_t ievent =0; ievent<fNevents; ++ievent){
    if ( fitType<2 ){
      TObjArray *events = (TObjArray*)(aType->At(sector));
      if ( events->GetSize()<=ievent ) break;
      TVectorD *v = (TVectorD*)(events->At(ievent));
      if ( (v!=0x0) && ((*xVar)[ievent]>0) ) { x[npoints]=(*xVar)[ievent]; y[npoints]=(*v)[fitParameter]; npoints++;}
    } else if (fitType == 2) {
      Double_t xValue=(*xVar)[ievent];
      Double_t yValue=0;
      if (sector>=0) yValue = (*GetTMeanEvents(sector))[ievent];
      else if (sector==-1) yValue=fVTime0SideA(ievent);
      else if (sector==-2) yValue=fVTime0SideC(ievent);
      if ( yValue>0 && xValue>0 ) { x[npoints]=xValue; y[npoints]=yValue;npoints++;}
    }else if (fitType == 3) {
      Double_t xValue=(*xVar)[ievent];
      Double_t yValue=(*GetQMeanEvents(sector))[ievent];
      if ( yValue>0 && xValue>0 ) { x[npoints]=xValue; y[npoints]=yValue;npoints++;}
    }
  }

  TGraph *gr = new TGraph(npoints);
    //sort xVariable increasing
  Int_t    *sortIndex = new Int_t[npoints];
  TMath::Sort(npoints,x,sortIndex, kFALSE);
  for (Int_t i=0;i<npoints;++i){
    gr->SetPoint(i,x[sortIndex[i]],y[sortIndex[i]]);
  }


  if ( xVariable == 2 ) delete xVar;
  delete [] x;
  delete [] y;
  delete [] sortIndex;
  return gr;
}
//_____________________________________________________________________
void AliTPCCalibCE::Analyse()
{
  //
  //  Calculate calibration constants
  //

  if (fProcessOld){
    TVectorD paramQ(3);
    TVectorD paramT0(3);
    TVectorD paramRMS(3);
    TMatrixD dummy(3,3);
    
    Float_t channelCounter=0;
    fMeanT0rms=0;
    fMeanQrms=0;
    fMeanRMSrms=0;
    
    for (Int_t iSec=0; iSec<72; ++iSec){
      TH2S *hT0 = GetHistoT0(iSec);
      if (!hT0 ) continue;
      
      AliTPCCalROC *rocQ     = GetCalRocQ  (iSec,kTRUE);
      AliTPCCalROC *rocT0    = GetCalRocT0 (iSec,kTRUE);
      AliTPCCalROC *rocT0Err = GetCalRocT0Err (iSec,kTRUE);
      AliTPCCalROC *rocRMS   = GetCalRocRMS(iSec,kTRUE);
      AliTPCCalROC *rocOut   = GetCalRocOutliers(iSec,kTRUE);
      
      TH2S *hQ   = GetHistoQ(iSec);
      TH2S *hRMS = GetHistoRMS(iSec);
      
      Short_t *arrayhQ   = hQ->GetArray();
      Short_t *arrayhT0  = hT0->GetArray();
      Short_t *arrayhRMS = hRMS->GetArray();
      
      UInt_t nChannels = fROC->GetNChannels(iSec);
      
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
        Float_t rms      = 0;
        Float_t rmsT0    = 0;
        
        
        Int_t offsetQ = (fNbinsQ+2)*(iChannel+1)+1;
        Int_t offsetT0 = (fNbinsT0+2)*(iChannel+1)+1;
        Int_t offsetRMS = (fNbinsRMS+2)*(iChannel+1)+1;
        
        cogQ     = AliMathBase::GetCOG(arrayhQ+offsetQ,fNbinsQ,fXminQ,fXmaxQ,&rms);
        fMeanQrms+=rms;
        cogTime0 = AliMathBase::GetCOG(arrayhT0+offsetT0,fNbinsT0,fXminT0,fXmaxT0,&rmsT0);
        fMeanT0rms+=rmsT0;
        cogRMS   = AliMathBase::GetCOG(arrayhRMS+offsetRMS,fNbinsRMS,fXminRMS,fXmaxRMS,&rms);
        fMeanRMSrms+=rms;
        channelCounter++;
        
      /*
             //outlier specifications
         if ( (cogQ < ??) && (cogTime0 > ??) && (cogTime0<??) && ( cogRMS>??) ){
         cogOut = 1;
          cogTime0 = 0;
          cogQ     = 0;
          cogRMS   = 0;
      }
      */
        rocQ->SetValue(iChannel, cogQ*cogQ);
        rocT0->SetValue(iChannel, cogTime0);
        rocT0Err->SetValue(iChannel, rmsT0);
        rocRMS->SetValue(iChannel, cogRMS);
        rocOut->SetValue(iChannel, cogOut);
        
        
      //debug
        if ( GetStreamLevel() > 0 ){
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
    if ( channelCounter>0 ){
      fMeanT0rms/=channelCounter;
      fMeanQrms/=channelCounter;
      fMeanRMSrms/=channelCounter;
    }
    //   if ( fDebugStreamer ) fDebugStreamer->GetFile()->Write();
    //    delete fDebugStreamer;
    //    fDebugStreamer = 0x0;
    fVEventTime.ResizeTo(fNevents);
    fVEventNumber.ResizeTo(fNevents);
    fVTime0SideA.ResizeTo(fNevents);
    fVTime0SideC.ResizeTo(fNevents);
  }

  if (fProcessNew&&fAnalyseNew){
    AnalyseTrack();
    for (Int_t iburst=0; iburst<fArrHnDrift.GetEntries(); ++iburst){
      THnSparseI *h=(THnSparseI*)fArrHnDrift.UncheckedAt(iburst);
      h->GetAxis(4)->SetRangeUser(fTimeBursts[iburst]-60*5,fTimeBursts[iburst]+60*5);
    }
  }
}




//
// New functions that also use the laser tracks
//



//_____________________________________________________________________
void AliTPCCalibCE::FindLocalMaxima(TObjArray * const arrObj, Double_t timestamp, Int_t burst)
{
  //
  //Find the local maximums for the projections to each axis
  //
  
  //find laser layer positoins
  fHnDrift->GetAxis(4)->SetRangeUser(timestamp-2*60,timestamp+2*60);
  FindLaserLayers();
  THnSparse *hProj=fHnDrift;
  Double_t posCE[4]={0.,0.,0.,0.};
  Double_t widthCE[4]={0.,0.,0.,0.};
  
//   if(fPeaks[4]!=0){
  // find central electrode position once more, separately for IROC, OROC, A-, C-Side
  
  for (Int_t i=0; i<4; ++i){
    hProj->GetAxis(0)->SetRangeUser(i*18,(i+1)*18-1);
    TH1 *h=fHnDrift->Projection(3);
    h->GetXaxis()->SetRangeUser(fPeaks[4]-fPeakWidths[4],fPeaks[4]+fPeakWidths[4]);
    Int_t nbinMax=h->GetMaximumBin();
    Double_t maximum=h->GetMaximum();
//     Double_t maxExpected=fNevents/fArrHnDrift->GetEntries()*556568./5./10.;
//     if (nbinMax<700||maximum<maxExpected) continue;
    Double_t xbinMax=h->GetBinCenter(nbinMax);
    TF1 fgaus("gaus","gaus",xbinMax-10,xbinMax+10);
    fgaus.SetParameters(maximum,xbinMax,2);
    fgaus.SetParLimits(1,xbinMax-5.,xbinMax+5.);
    fgaus.SetParLimits(2,0.2,4.);
    h->Fit(&fgaus,"RQN");
//     Double_t deltaX=4*fgaus.GetParameter(2);
//     xbinMax=fgaus.GetParameter(1);
    delete h;
    posCE[i]=fgaus.GetParameter(1);
    widthCE[i]=4*fgaus.GetParameter(2);
    hProj->GetAxis(0)->SetRangeUser(0,72);
  }
//   }
  //Current drift velocity
  Float_t vdrift=2.61301900000000000e+06;//fParam->GetDriftV();
//   cout<<"vdrift="<<vdrift<<endl;
  
  AliDebug(5,Form("Timestamp %f - default drift velocity %f",timestamp,vdrift));
  //loop over all entries in the histogram
  Int_t coord[5];
  for(Long64_t ichunk=0;ichunk<hProj->GetNbins();ichunk++){
    //get entry position and content
    Double_t adc=hProj->GetBinContent(ichunk,coord);
    
    
    Int_t sector = coord[0]-1;
    Int_t row    = coord[1]-1;
    Int_t pad    = coord[2]-1;
    Int_t timeBin= coord[3]-1;
    Double_t time   = fHnDrift->GetAxis(4)->GetBinCenter(coord[4]);
    Int_t side   = (sector/18)%2;
//     return;
//     fPeaks[4]=(UInt_t)posCE[sector/18];
//     fPeakWidths[4]=(UInt_t)widthCE[sector/18];
    
    //cuts
    if (time<timestamp-120||time>timestamp+120) continue; //window of +- 2min
    if (adc < 5 ) continue;
    if (IsEdgePad(sector,row,pad)) continue;
//     if (!IsPeakInRange(timeBin)) continue;
//     if (isCE&&((row%2)||(row%2)||(sector%2))) continue;
//     if (isCE&&(sector!=0)) continue;
    
    Int_t padmin=-2, padmax=2;
    Int_t timemin=-2, timemax=2;
    Int_t minsumperpad=2;
    //CE or laser tracks
    Bool_t isCE=kFALSE;
    if (TMath::Abs((Short_t)timeBin-(Short_t)posCE[sector/18])<(Short_t)widthCE[sector/18]) {
      isCE=kTRUE;
      padmin=0;
      padmax=0;
      timemin=-3;
      timemax=7;
    }
    
    //
    // Find local maximum and cogs
    //
    Bool_t isMaximum=kTRUE;
    Float_t totalmass=0, tbcm=0, padcm=0, rmstb=0, rmspad=0;
    Double_t cogY=0, rmsY=0;
    Int_t npart=0;
    
    // for position calculation use
    for(Int_t ipad=padmin;ipad<=padmax;++ipad){
      Float_t lxyz[3];
      fROC->GetPositionLocal(sector,row,pad+ipad,lxyz);
      
      for(Int_t itime=timemin;itime<=timemax;++itime){
        
        Int_t a[5]={coord[0],coord[1],coord[2]+ipad,coord[3]+itime,coord[4]};
        Double_t val=hProj->GetBinContent(a);
        totalmass+=val;
        
        tbcm +=(timeBin+itime)*val;
        padcm+=(pad+ipad)*val;
        cogY +=lxyz[1]*val;
        
        rmstb +=(timeBin+itime)*(timeBin+itime)*val;
        rmspad+=(pad+ipad)*(pad+ipad)*val;
        rmsY  +=lxyz[1]*lxyz[1]*val;
        
        if (val>0) ++npart;
        if (val>adc) {
          isMaximum=kFALSE;
          break;
        }
      }
      if (!isMaximum)  break;
    }
    
    if (!isMaximum||!npart)  continue;
    if (totalmass<npart*minsumperpad) continue;
    if (!isCE&&rmspad<.1) continue; //most probably noise, since signal only in one pad,
                                    //for CE we want only one pad by construction
    
    tbcm/=totalmass;
    padcm/=totalmass;
    cogY/=totalmass;
    
    rmstb/=totalmass;
    rmspad/=totalmass;
    rmsY/=totalmass;
    
    rmstb=TMath::Sqrt(TMath::Abs(tbcm*tbcm-rmstb));
    rmspad=TMath::Sqrt(TMath::Abs(padcm*padcm-rmspad));
    rmsY=TMath::Sqrt(TMath::Abs(cogY*cogY-rmsY));
    
    Int_t cog=TMath::Nint(padcm);
    
    // timebin --> z position
    Float_t zlength=fParam->GetZLength(side);
//     Float_t timePos=tbcm+(1000-fPeaks[4])
    // drift velocity is in m/s we would like to have cm/100ns, so 100cm/(10^7*100ns)
    Double_t gz=(zlength-(tbcm*vdrift*1.e-7))*TMath::Power(-1,side);
    
    // local to global transformation--> x and y positions
    Float_t padlxyz[3];
    fROC->GetPositionLocal(sector,row,pad,padlxyz);
    
    Double_t gxyz[3]={padlxyz[0],cogY,gz};
    Double_t lxyz[3]={padlxyz[0],cogY,gz};
    Double_t igxyz[3]={0,0,0};
    AliTPCTransform t1;
    t1.RotatedGlobal2Global(sector,gxyz);
    
    Double_t mindist=0;
    Int_t trackID=-1;
    Int_t trackID2=-1;
    
    //find track id and index of the position in the track (row)
    Int_t index=0;
    if (!isCE){
      index=row+(sector>35)*fROC->GetNRows(0);
      trackID=FindLaserTrackID(sector,index,gxyz,mindist,lxyz,trackID2);
    } else {
      trackID=336+((sector/18)%2);
      index= fROC->GetRowIndexes(sector)[row]+pad; //  global pad position in sector
      if (sector<36) {
        index+=(sector%18)*fROC->GetNChannels(sector);
      } else {
        index+=18*fROC->GetNChannels(0);
        index+=(sector%18)*fROC->GetNChannels(sector);
      }
      //TODO: find out about the multiple peaks in the CE
//       mindist=TMath::Abs(fPeaks[4]-tbcm);
      mindist=1.;
    }
    
    // fill track vectors
    if (trackID>0){
      AliTPCLaserTrack *ltr=(AliTPCLaserTrack*)arrObj->UncheckedAt(trackID);
      Double_t oldMinDist=ltr->fVecPhi->GetMatrixArray()[index];
      
      
//      travel time effect of light includes
      
      Double_t raylength=ltr->GetRayLength();
      Double_t globmir[3]={ltr->Xv(),ltr->Yv(),ltr->Zv()};
      ltr->GetXYZ(globmir);
      if(trackID<336){
        if(side==0){
          gxyz[2]=gxyz[2]-(TMath::Sqrt((gxyz[0]-globmir[0])*(gxyz[0]-globmir[0])
                                       +(gxyz[1]-globmir[1])*(gxyz[1]-globmir[1])
                                       +(gxyz[2]-globmir[2])*(gxyz[2]-globmir[2])+raylength))*vdrift*TMath::Power(10.,-6.)/30000;
        }
        else {
          gxyz[2]=gxyz[2]-(TMath::Sqrt((gxyz[0]-globmir[0])*(gxyz[0]-globmir[0])
                                       +(gxyz[1]-globmir[1])*(gxyz[1]-globmir[1])
                                       +(gxyz[2]-globmir[2])*(gxyz[2]-globmir[2])+raylength))*vdrift*TMath::Power(10.,-6.)/30000;
        }
      }
      
      if (TMath::Abs(oldMinDist)<1.e-20||oldMinDist>mindist){
        ltr->fVecSec->GetMatrixArray()[index]=sector;
        ltr->fVecP2->GetMatrixArray()[index]=0;
        ltr->fVecPhi->GetMatrixArray()[index]=mindist;
        ltr->fVecGX->GetMatrixArray()[index]=gxyz[0];
        ltr->fVecGY->GetMatrixArray()[index]=gxyz[1];
        ltr->fVecGZ->GetMatrixArray()[index]=gxyz[2];
        ltr->fVecLX->GetMatrixArray()[index]=lxyz[0];
        ltr->fVecLY->GetMatrixArray()[index]=lxyz[1];
        ltr->fVecLZ->GetMatrixArray()[index]=lxyz[2];
//         ltr->SetUniqueID((UInt_t)(mindist*10000)); //distance in um
      }
      TObjArray *arr=AliTPCLaserTrack::GetTracks();
      ltr=(AliTPCLaserTrack*)arr->UncheckedAt(trackID);
      igxyz[0]=ltr->fVecGX->GetMatrixArray()[row];
      igxyz[1]=ltr->fVecGY->GetMatrixArray()[row];
      igxyz[2]=ltr->fVecGZ->GetMatrixArray()[row];
    }
    
    
    if (fStreamLevel>4){
      (*GetDebugStreamer()) << "clusters" <<
        "run="   << fRunNumber <<
        "timestamp=" << timestamp <<
        "burst="     << burst     <<
        "side="      << side      <<
        "sec="       << sector    <<
        "row="       << row       <<
        "pad="       << pad       <<
        "padCog="    << cog       <<
        "timebin="   << timeBin   <<
        "cogCE="     << posCE[sector/18] <<
        "withCE="    << widthCE[sector/18] <<
        "index="     << index     <<
        
        "padcm="     << padcm     <<
        "rmspad="    << rmspad    <<
        
        "cogtb="     << tbcm      <<
        "rmstb="     << rmstb     <<
        
        "npad="      << npart     <<
        
        "lx="        << padlxyz[0]<<
        "ly="        << cogY      <<
        "lypad="     << padlxyz[1]<<
        "rmsY="      << rmsY      <<
        
        "gx="        << gxyz[0]   <<
        "gy="        << gxyz[1]   <<
        "gz="        << gxyz[2]   <<
        
        "igx="        << igxyz[0] <<
        "igy="        << igxyz[1] <<
        "igz="        << igxyz[2] <<
        
        "mind="      << mindist   <<
        "max="       << adc       <<
        "trackid="   << trackID   <<
        "trackid2="   << trackID2   <<
        "npart="     << npart     <<
        "\n";
    } // end stream levelmgz.fElements
    
  }
  
}

//_____________________________________________________________________
void AliTPCCalibCE::AnalyseTrack()
{
  //
  //  Analyse the tracks
  //
  
  
  AliTPCLaserTrack::LoadTracks();
//   AliTPCParam *param=0x0;
//   //cdb run number
//   AliCDBManager *man=AliCDBManager::Instance();
//   if (man->GetDefaultStorage()){
//     AliCDBEntry *entry=man->Get("TPC/Calib/Parameters",fRunNumber);
//     if (entry){
//       entry->SetOwner(kTRUE);
//       param = (AliTPCParam*)(entry->GetObject()->Clone());
//     }
//   }
//   if (param){
//     if (fParam) delete fParam;
//     fParam=param;
//   } else {
//     AliError("Could not get updated AliTPCParam from OCDB!!!");
//   }
  
  //Measured and ideal laser tracks
  TObjArray* arrMeasured = SetupMeasured();
  TObjArray* arrIdeal    = AliTPCLaserTrack::GetTracks();
  AddCEtoIdeal(arrIdeal);
  
  //find bursts and loop over them
  for (Int_t iburst=0; iburst<fArrHnDrift.GetEntries();++iburst){
    Double_t timestamp=fTimeBursts[iburst];
    AliDebug(5,Form("Burst: %d (%f)",iburst,timestamp));
    fHnDrift=(THnSparseI*)fArrHnDrift.UncheckedAt(iburst);
    if (!fHnDrift) continue;
    UInt_t entries=(UInt_t)fHnDrift->GetEntries();
    if (fBinsLastAna[iburst]>=entries) continue; //already analysed!!!
    fBinsLastAna[iburst]=entries;
    
    for (Int_t iDim=0; iDim<fHnDrift->GetNdimensions(); ++iDim) fHnDrift->GetAxis(iDim)->SetRange(0,0);
//     if (iburst==0) FindLaserLayers();
    
    //reset laser tracks
    ResetMeasured(arrMeasured);
    
    //find clusters and associate to the tracks
    FindLocalMaxima(arrMeasured, timestamp, iburst);
    
    //calculate drift velocity
    CalculateDV(arrIdeal,arrMeasured,iburst);
    
    //Dump information to file if requested
    if (fStreamLevel>2){
      printf("make tree\n");
      //laser track information
      
      for (Int_t itrack=0; itrack<338; ++itrack){
        TObject *iltr=arrIdeal->UncheckedAt(itrack);
        TObject *mltr=arrMeasured->UncheckedAt(itrack);
        (*GetDebugStreamer()) << "tracks" <<
          "run="   << fRunNumber <<
          "time=" << timestamp <<
          "burst="<< iburst <<
          "iltr.=" << iltr <<
          "mltr.=" << mltr <<
          "\n";
      }
    }
  }
  if (fStreamLevel>0) GetDebugStreamer()->GetFile()->Write();
}

//_____________________________________________________________________
Int_t AliTPCCalibCE::FindLaserTrackID(Int_t sector,Int_t row, const Double_t *peakpos,Double_t &mindist,
                                      const Double_t *peakposloc, Int_t &itrackMin2)
{
  //
  //  Find the tracks, which are closest to the ideal tracks, from clusters closest to the ideal tracks
  //
  
  
  TObjArray *arr=AliTPCLaserTrack::GetTracks();
  TVector3 vP(peakpos[0],peakpos[1],peakpos[2]);
  TVector3 vDir;
  TVector3 vSt;
  
  Int_t firstbeam=0;
  Int_t lastbeam=336/2;
  if ( (sector/18)%2 ) {
    firstbeam=336/2;
    lastbeam=336;
  }
  
  mindist=1000000;
  Int_t itrackMin=-1;
  for (Int_t itrack=firstbeam; itrack<lastbeam; ++itrack){
    AliTPCLaserTrack *ltr=(AliTPCLaserTrack*)arr->At(itrack);  //get the track
//     if (ltr->GetVecSec()->GetMatrixArray()[row]!=sector) continue;
    vSt.SetXYZ(ltr->GetX(),ltr->GetY(),ltr->GetZ());
    Double_t deltaZ=ltr->GetZ()-peakpos[2];
    if (TMath::Abs(deltaZ)>40) continue;
    vDir.SetMagThetaPhi(1,ltr->Theta(),TMath::ASin(ltr->GetSnp()));
    vSt.RotateZ(ltr->GetAlpha());
    vDir.RotateZ(ltr->GetAlpha());
    
    Double_t dist=(vDir.Cross(vSt-vP)).Mag()/vDir.Mag();
    
    if (dist<mindist){
      mindist=dist;
      itrackMin=itrack;
    }
    
  }
  itrackMin2=-1;
  Float_t mindist2=10;
  for (Int_t itrack=firstbeam; itrack<lastbeam; ++itrack){
    AliTPCLaserTrack *ltr=(AliTPCLaserTrack*)arr->At(itrack);  //get the track
    if ((ltr->fVecSec->GetMatrixArray())[row]!=sector) continue;
    
    Double_t deltaZ=ltr->GetZ()-peakpos[2];
    if (TMath::Abs(deltaZ)>40) continue;
    
    Double_t dist=TMath::Abs((ltr->fVecLY->GetMatrixArray())[row]-peakposloc[1]);
    if (dist>1) continue;
    
    if (dist<mindist2){
      mindist2=dist;
      itrackMin2=itrack;
    }
  }
  mindist=mindist2;
  return itrackMin2;
  
}

//_____________________________________________________________________
Bool_t AliTPCCalibCE::IsEdgePad(Int_t sector, Int_t row, Int_t pad) const
{
  //
  // return true if pad is on the edge of a row
  //
  Int_t edge1   = 0;
  if ( pad == edge1 ) return kTRUE;
  Int_t edge2   = fROC->GetNPads(sector,row)-1;
  if ( pad == edge2 ) return kTRUE;
  
  return kFALSE;
}

//_____________________________________________________________________
TObjArray* AliTPCCalibCE::SetupMeasured()
{
  //
  // setup array of measured laser tracks and CE
  //
  
  TObjArray *arrIdeal    = AliTPCLaserTrack::GetTracks();
  TObjArray *arrMeasured = new TObjArray(338);
  arrMeasured->SetOwner();
  for(Int_t itrack=0;itrack<336;itrack++){
    arrMeasured->AddAt(new AliTPCLaserTrack(*((AliTPCLaserTrack*)arrIdeal->At(itrack))),itrack);
  }
  
  //"tracks" for CE
  AliTPCLaserTrack *ltrce=new AliTPCLaserTrack;
  ltrce->SetId(336);
  ltrce->SetSide(0);
  ltrce->fVecSec=new TVectorD(557568/2);
  ltrce->fVecP2=new TVectorD(557568/2);
  ltrce->fVecPhi=new TVectorD(557568/2);
  ltrce->fVecGX=new TVectorD(557568/2);
  ltrce->fVecGY=new TVectorD(557568/2);
  ltrce->fVecGZ=new TVectorD(557568/2);
  ltrce->fVecLX=new TVectorD(557568/2);
  ltrce->fVecLY=new TVectorD(557568/2);
  ltrce->fVecLZ=new TVectorD(557568/2);
  
  arrMeasured->AddAt(ltrce,336); //CE A-Side
  
  ltrce=new AliTPCLaserTrack;
  ltrce->SetId(337);
  ltrce->SetSide(1);
  ltrce->fVecSec=new TVectorD(557568/2);
  ltrce->fVecP2=new TVectorD(557568/2);
  ltrce->fVecPhi=new TVectorD(557568/2);
  ltrce->fVecGX=new TVectorD(557568/2);
  ltrce->fVecGY=new TVectorD(557568/2);
  ltrce->fVecGZ=new TVectorD(557568/2);
  ltrce->fVecLX=new TVectorD(557568/2);
  ltrce->fVecLY=new TVectorD(557568/2);
  ltrce->fVecLZ=new TVectorD(557568/2);
  arrMeasured->AddAt(ltrce,337); //CE C-Side
  
  return arrMeasured;
}

//_____________________________________________________________________
void AliTPCCalibCE::ResetMeasured(TObjArray * const arr)
{
  //
  // reset array of measured laser tracks and CE
  //
  for(Int_t itrack=0;itrack<338;itrack++){
    AliTPCLaserTrack *ltr=(AliTPCLaserTrack*)arr->UncheckedAt(itrack);
    ltr->fVecSec->Zero();
    ltr->fVecP2->Zero();
    ltr->fVecPhi->Zero();
    ltr->fVecGX->Zero();
    ltr->fVecGY->Zero();
    ltr->fVecGZ->Zero();
    ltr->fVecLX->Zero();
    ltr->fVecLY->Zero();
    ltr->fVecLZ->Zero();
  }
}

//_____________________________________________________________________
void AliTPCCalibCE::AddCEtoIdeal(TObjArray *arr)
{
  //
  // Add ideal CE positions to the ideal track data
  //
  
  arr->Expand(338);
  //"tracks" for CE
  AliTPCLaserTrack *ltrceA=new AliTPCLaserTrack;
  ltrceA->SetId(336);
  ltrceA->SetSide(0);
  ltrceA->fVecSec=new TVectorD(557568/2);
  ltrceA->fVecP2=new TVectorD(557568/2);
  ltrceA->fVecPhi=new TVectorD(557568/2);
  ltrceA->fVecGX=new TVectorD(557568/2);
  ltrceA->fVecGY=new TVectorD(557568/2);
  ltrceA->fVecGZ=new TVectorD(557568/2);
  ltrceA->fVecLX=new TVectorD(557568/2);
  ltrceA->fVecLY=new TVectorD(557568/2);
  ltrceA->fVecLZ=new TVectorD(557568/2);
  arr->AddAt(ltrceA,336); //CE A-Side
  
  AliTPCLaserTrack *ltrceC=new AliTPCLaserTrack;
  ltrceC->SetId(337);
  ltrceC->SetSide(1);
  ltrceC->fVecSec=new TVectorD(557568/2);
  ltrceC->fVecP2=new TVectorD(557568/2);
  ltrceC->fVecPhi=new TVectorD(557568/2);
  ltrceC->fVecGX=new TVectorD(557568/2);
  ltrceC->fVecGY=new TVectorD(557568/2);
  ltrceC->fVecGZ=new TVectorD(557568/2);
  ltrceC->fVecLX=new TVectorD(557568/2);
  ltrceC->fVecLY=new TVectorD(557568/2);
  ltrceC->fVecLZ=new TVectorD(557568/2);
  arr->AddAt(ltrceC,337); //CE C-Side
  
  //Calculate ideal positoins
  Float_t gxyz[3];
  Float_t lxyz[3];
  Int_t channelSideA=0;
  Int_t channelSideC=0;
  Int_t channelSide=0;
  AliTPCLaserTrack *ltrce=0x0;
  
  for (Int_t isector=0; isector<72; ++isector){
    Int_t side=((isector/18)%2);
    for (UInt_t irow=0;irow<fROC->GetNRows(isector);++irow){
      for (UInt_t ipad=0;ipad<fROC->GetNPads(isector,irow);++ipad){
        fROC->GetPositionGlobal(isector,irow,ipad,gxyz);
        fROC->GetPositionLocal(isector,irow,ipad,lxyz);
        if (side==0) {
          ltrce=ltrceA;
          channelSide=channelSideA;
        } else {
          ltrce=ltrceC;
          channelSide=channelSideC;
        }
        
        ltrce->fVecSec->GetMatrixArray()[channelSide]=isector;
        ltrce->fVecP2->GetMatrixArray()[channelSide]=0;
        ltrce->fVecPhi->GetMatrixArray()[channelSide]=0;
        ltrce->fVecGX->GetMatrixArray()[channelSide]=gxyz[0];
        ltrce->fVecGY->GetMatrixArray()[channelSide]=gxyz[1];
//         ltrce->fVecGZ->GetMatrixArray()[channelSide]=-1;
        ltrce->fVecLX->GetMatrixArray()[channelSide]=lxyz[0];
        ltrce->fVecLY->GetMatrixArray()[channelSide]=lxyz[1];
//         ltrce->fVecLZ->GetMatrixArray()[channelSide]=-1;
        
        if (side==0){
          ltrce->fVecGZ->GetMatrixArray()[channelSide]=-0.335;
          ltrce->fVecLZ->GetMatrixArray()[channelSide]=-0.335;
          ++channelSideA;
        }
        else {
          ltrce->fVecGZ->GetMatrixArray()[channelSide]=0.15;
          ltrce->fVecLZ->GetMatrixArray()[channelSide]=0.15;
          ++channelSideC;
        }
      }
    }
  }
  
  
}

//_____________________________________________________________________
void AliTPCCalibCE::CalculateDV(TObjArray * const arrIdeal, TObjArray * const arrMeasured, Int_t burst)
{
  //
  // calculate the drift velocity from the reconstructed clusters associated
  // to the ideal laser tracks
  // use two different fit scenarios: Separate fit for A- and C-Side
  //                                  Common fit for A- and C-Side
  //
  
  if (!fArrFitGraphs){
    fArrFitGraphs=new TObjArray;
  }
  
//   static TLinearFitter fdriftA(5,"hyp4");
//   static TLinearFitter fdriftC(5,"hyp4");
//   static TLinearFitter fdriftAC(6,"hyp5");
  Double_t timestamp=fTimeBursts[burst];
  
  static TLinearFitter fdriftA(4,"hyp3");
  static TLinearFitter fdriftC(4,"hyp3");
  static TLinearFitter fdriftAC(5,"hyp4");
  TVectorD fitA(7),fitC(7),fitAC(8); //fit values+chi2+npoints
  
  Float_t chi2A = 10;
  Float_t chi2C = 10;
  Float_t chi2AC = 10;
  Int_t npointsA=0;
  Int_t npointsC=0;
  Int_t npointsAC=0;
  
  Double_t minres[3]={20.,1,0.8};
  //----
  for(Int_t i=0;i<3;i++){
    
    fdriftA.ClearPoints();
    fdriftC.ClearPoints();
    fdriftAC.ClearPoints();
    
    chi2A = 10;
    chi2C = 10;
    chi2AC = 10;
    npointsA=0;
    npointsC=0;
    npointsAC=0;
    
    for (Int_t itrack=0; itrack<338; ++itrack){
      AliTPCLaserTrack *iltr=(AliTPCLaserTrack*)arrIdeal->UncheckedAt(itrack);
      AliTPCLaserTrack *mltr=(AliTPCLaserTrack*)arrMeasured->UncheckedAt(itrack);

      //-- Exclude the tracks which has the biggest inclanation angle
      if ((itrack%7==0||itrack%7==6)&&itrack<336) continue;
      Int_t clustercounter=0;
      Int_t indexMax=159;
      
      //-- exclude the low intensity tracks
      
      for (Int_t index=0; index<indexMax; ++index){
        
        Double_t mGx=mltr->fVecGX->GetMatrixArray()[index];
        Double_t mGy=mltr->fVecGY->GetMatrixArray()[index];
        Double_t mGz=mltr->fVecGZ->GetMatrixArray()[index];
        
        if (TMath::Abs(mGz)<1e-20 && TMath::Abs(mGy)<1e-20 && TMath::Abs(mGx)<1e-20) clustercounter++;
      }
      if (clustercounter>130&&itrack<336) continue; // don't accept tracks with <= 159-130=29 clusters
      clustercounter=0;

      
      //-- drift length
      Double_t zlength = (iltr->GetSide()==0)? fParam->GetZLength(36): fParam->GetZLength(71);
      
      if (itrack>335) indexMax=557568/2;
      for (Int_t index=0; index<indexMax; ++index){
        Double_t iGx=iltr->fVecGX->GetMatrixArray()[index];
        Double_t iGy=iltr->fVecGY->GetMatrixArray()[index];
        Double_t iGz=iltr->fVecGZ->GetMatrixArray()[index];
        Double_t iR=TMath::Sqrt(iGx*iGx+iGy*iGy);
        
        Double_t mGx=mltr->fVecGX->GetMatrixArray()[index];
        Double_t mGy=mltr->fVecGY->GetMatrixArray()[index];
        Double_t mGz=mltr->fVecGZ->GetMatrixArray()[index];
        Double_t mR=TMath::Sqrt(mGx*mGx+mGy*mGy);
        
      //cut if no track info available
        if (iltr->GetBundle()==0) continue;
        if (iR<133||mR<133) continue;
        if(mltr->fVecP2->GetMatrixArray()[index]>minres[i]) continue;
        
        Double_t ldrift  = (iltr->GetSide()==0)?zlength-iGz:iGz+zlength;
        Double_t mdrift  = (iltr->GetSide()==0)?zlength-mGz:mGz+zlength;
        
        //Double_t xxx[4] = {ldrift,iGy*ldrift/(zlength*250.), 250.-mR, iltr->fVecSec->GetMatrixArray()[index]>35};
        Double_t xxx[3] = {ldrift,iGy*ldrift/(zlength*250.), 250.-mR};
        
        if (iltr->GetSide()==0){
          fdriftA.AddPoint(xxx,mdrift,1);
        }else{
          fdriftC.AddPoint(xxx,mdrift,1);
        }
//         Double_t xxx2[4] = { ldrift,iGy*ldrift/(zlength*250.), 250.-mR, iltr->fVecSec->GetMatrixArray()[index]>35, iltr->GetSide()};
        Double_t xxx2[4] = { ldrift,iGy*ldrift/(zlength*250.), 250.-mR, iltr->GetSide()};
        fdriftAC.AddPoint(xxx2,mdrift,1);
        
      }//end index loop
    }//end laser track loop
    
  //perform fit
    fdriftA.Eval();
    fdriftC.Eval();
    fdriftAC.Eval();
    
    
    
  //get fit values
    fdriftA.GetParameters(fitA);
    fdriftC.GetParameters(fitC);
    fdriftAC.GetParameters(fitAC);
    
  //Parameters:  0 linear offset
  //             1 mean drift velocity correction factor
  //             2 relative global y gradient
  //             3 radial deformation
  //             4 IROC/OROC offset
    
//      FindResiduals(arrMeasured,arrIdeal,fitA,fitC);
    
    for (Int_t itrack=0; itrack<338; ++itrack){
      AliTPCLaserTrack *iltr=(AliTPCLaserTrack*)arrIdeal->UncheckedAt(itrack);
      AliTPCLaserTrack *mltr=(AliTPCLaserTrack*)arrMeasured->UncheckedAt(itrack);

      //-- Exclude the tracks which has the biggest inclanation angle
      if ((itrack%7==0||itrack%7==6)&&itrack<336) continue;
      Int_t clustercounter=0;
      Int_t indexMax=159;

      //-- exclude the low intensity tracks

      for (Int_t index=0; index<indexMax; ++index){
        Double_t mGx=mltr->fVecGX->GetMatrixArray()[index];
        Double_t mGy=mltr->fVecGY->GetMatrixArray()[index];
        Double_t mGz=mltr->fVecGZ->GetMatrixArray()[index];
        if (TMath::Abs(mGz)<1e-20 && TMath::Abs(mGy)<1e-20 && TMath::Abs(mGx)<1e-20) clustercounter++;
      }
      if (clustercounter>130&&itrack<336) continue; // don't accept tracks with <= 159-130=29 clusters
      clustercounter=0;

      //-- drift length
      Double_t zlength = (iltr->GetSide()==0)? fParam->GetZLength(36): fParam->GetZLength(71);

      if (itrack>335) indexMax=557568/2;
      for (Int_t index=0; index<indexMax; ++index){
        Double_t iGx=iltr->fVecGX->GetMatrixArray()[index];
        Double_t iGy=iltr->fVecGY->GetMatrixArray()[index];
        Double_t iGz=iltr->fVecGZ->GetMatrixArray()[index];
        Double_t iR=TMath::Sqrt(iGx*iGx+iGy*iGy);
        
        Double_t mGx=mltr->fVecGX->GetMatrixArray()[index];
        Double_t mGy=mltr->fVecGY->GetMatrixArray()[index];
        Double_t mGz=mltr->fVecGZ->GetMatrixArray()[index];
        Double_t mR=TMath::Sqrt(mGx*mGx+mGy*mGy);
        
      //cut if no track info available
        if (iR<60||mR<60) continue;
        
        Double_t ldrift  = (iltr->GetSide()==0)?zlength-iGz:iGz+zlength;
        Double_t mdrift  = (iltr->GetSide()==0)?zlength-mGz:mGz+zlength;
        
        TVectorD *v=&fitA;
        if (iltr->GetSide()==1) v=&fitC;
//         Double_t iCorr=(*v)[0]+(*v)[1]*ldrift+(*v)[2]*iGy*ldrift/(zlength*250.)+(*v)[3]*(250.-mR)+(*v)[4]*( iltr->fVecSec->GetMatrixArray()[index]>35);
        Double_t iCorr=(*v)[0]+(*v)[1]*ldrift+(*v)[2]*iGy*ldrift/(zlength*250.)+(*v)[3]*(250.-mR);
        
        mltr->fVecP2->GetMatrixArray()[index]=mdrift-iCorr;
        
      }
    }
    
    fitA.ResizeTo(7);
    fitC.ResizeTo(7);
    fitAC.ResizeTo(8);
    
//set statistics values

    npointsA= fdriftA.GetNpoints();
    if (npointsA>0) chi2A = fdriftA.GetChisquare()/fdriftA.GetNpoints();
    fitA[5]=npointsA;
    fitA[6]=chi2A;
    
    npointsC= fdriftC.GetNpoints();
    if (npointsC>0) chi2C = fdriftC.GetChisquare()/fdriftC.GetNpoints();
    fitC[5]=npointsC;
    fitC[6]=chi2C;
    
    npointsAC= fdriftAC.GetNpoints();
    if (npointsAC>0) chi2AC = fdriftAC.GetChisquare()/fdriftAC.GetNpoints();
    fitAC[5]=npointsAC;
    fitAC[6]=chi2AC;
    
    if (fStreamLevel>2){
    //laser track information
      (*GetDebugStreamer()) << "DriftV" <<
        "iter="   << i <<
        "run="    << fRunNumber <<
        "time="   << timestamp <<
        "fitA.="  << &fitA <<
        "fitC.="  << &fitC <<
        "fitAC.=" << &fitAC <<
        "\n";
      
      
    }
    
  }
//-----
  
  
  //Parameters:  0 linear offset (global)
  //             1 mean drift velocity correction factor
  //             2 relative global y gradient
  //             3 radial deformation
  //             4 IROC/OROC offset
  //             5 linear offset relative A-C
  
  //get graphs
  TGraphErrors *grA[7];
  TGraphErrors *grC[7];
  TGraphErrors *grAC[8];
  TString names("GRAPH_MEAN_DELAY_LASER_ALL_;GRAPH_MEAN_DRIFT_LASER_ALL_;GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_;GRAPH_MEAN_RGRADIENT_LASER_ALL_;GRAPH_MEAN_IROCOROCOFFSET_LASER_ALL_;GRAPH_MEAN_NPOINTS_LASER_ALL_;GRAPH_MEAN_CHI2_LASER_ALL_");
  TString namesAC("GRAPH_MEAN_DELAY_LASER_ALL_;GRAPH_MEAN_DRIFT_LASER_ALL_;GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_;GRAPH_MEAN_RGRADIENT_LASER_ALL_;GRAPH_MEAN_IROCOROCOFFSET_LASER_ALL_;GRAPH_MEAN_NPOINTS_LASER_ALL_;GRAPH_MEAN_CHI2_LASER_ALL_;GRAPH_MEAN_DELAYC_LASER_ALL_");
  
  TObjArray *arrNames=names.Tokenize(";");
  TObjArray *arrNamesAC=namesAC.Tokenize(";");
  
  //A-Side graphs
  for (Int_t i=0; i<7; ++i){
    TString grName=arrNames->UncheckedAt(i)->GetName();
    grName+="A";
    grA[i]=(TGraphErrors*)fArrFitGraphs->FindObject(grName.Data());
    if (!grA[i]){
      grA[i]=new TGraphErrors;
      grA[i]->SetName(grName.Data());
      grA[i]->SetTitle(grName.ReplaceAll("_"," ").Data());
      fArrFitGraphs->Add(grA[i]);
    }
//     Int_t ipoint=grA[i]->GetN();
    Int_t ipoint=burst;
    grA[i]->SetPoint(ipoint,timestamp,fitA(i));
    grA[i]->SetPointError(ipoint,60,.0001);
    if (i<4) grA[i]->SetPointError(ipoint,60,fdriftA.GetCovarianceMatrixElement(i,i));
  }
  
  //C-Side graphs
  for (Int_t i=0; i<7; ++i){
    TString grName=arrNames->UncheckedAt(i)->GetName();
    grName+="C";
    grC[i]=(TGraphErrors*)fArrFitGraphs->FindObject(grName.Data());
    if (!grC[i]){
      grC[i]=new TGraphErrors;
      grC[i]->SetName(grName.Data());
      grC[i]->SetTitle(grName.ReplaceAll("_"," ").Data());
      fArrFitGraphs->Add(grC[i]);
    }
//     Int_t ipoint=grC[i]->GetN();
    Int_t ipoint=burst;
    grC[i]->SetPoint(ipoint,timestamp,fitC(i));
    grC[i]->SetPointError(ipoint,60,.0001);
    if (i<4) grC[i]->SetPointError(ipoint,60,fdriftA.GetCovarianceMatrixElement(i,i));
  }
  
  //AC-Side graphs
  for (Int_t i=0; i<8; ++i){
    TString grName=arrNamesAC->UncheckedAt(i)->GetName();
    grName+="AC";
    grAC[i]=(TGraphErrors*)fArrFitGraphs->FindObject(grName.Data());
    if (!grAC[i]){
      grAC[i]=new TGraphErrors;
      grAC[i]->SetName(grName.Data());
      grAC[i]->SetTitle(grName.ReplaceAll("_"," ").Data());
      fArrFitGraphs->Add(grAC[i]);
    }
//     Int_t ipoint=grAC[i]->GetN();
    Int_t ipoint=burst;
    grAC[i]->SetPoint(ipoint,timestamp,fitAC(i));
    grAC[i]->SetPointError(ipoint,60,.0001);
    if (i<5) grAC[i]->SetPointError(ipoint,60,fdriftA.GetCovarianceMatrixElement(i,i));
  }
  
  if (fDebugLevel>10){
    printf("A side fit parameters:\n");
    fitA.Print();
    printf("\nC side fit parameters:\n");
    fitC.Print();
    printf("\nAC side fit parameters:\n");
    fitAC.Print();
  }
  delete arrNames;
  delete arrNamesAC;
}

//_____________________________________________________________________
Double_t AliTPCCalibCE::SetBurstHnDrift()
{
  //
  // Create a new THnSparse for the current burst
  // return the time of the current burst
  //
  Int_t i=0;
  for(i=0; i<fTimeBursts.GetNrows(); ++i){
    if(fTimeBursts.GetMatrixArray()[i]<1.e-20) break;
    if(TMath::Abs(fTimeBursts.GetMatrixArray()[i]-fTimeStamp)<300){
      fHnDrift=(THnSparseI*)fArrHnDrift.UncheckedAt(i);
      return fTimeBursts(i);
    }
  }
  
  CreateDVhist();
  fArrHnDrift.AddAt(fHnDrift,i);
  fTimeBursts.GetMatrixArray()[i]=fTimeStamp;
//   fPeaks[4]=0;
  return fTimeStamp;
}

//_____________________________________________________________________
void AliTPCCalibCE::DumpToFile(const Char_t *filename, const Char_t *dir, Bool_t /*append*/)
{
  //
  //  Write class to file
  //  option can be specified in the dir option:
  //  options:
  //    name=<objname>: the name of the calibration object in file will be <objname>
  //    type=<type>:    the saving type:
  //                    0 - write the complte object
  //                    1 - Store the histogram arrays separately to make the streamed object smaller, Analyse to be called
  //                    2 - like 2, but in addition delete objects that will most probably not be used for calibration
  //                    3 - store only calibration output, don't store the reference histograms
  //                        and THnSparse (requires Analyse called before)
  //
  //  NOTE: to read the object back, the ReadFromFile function should be used
  //
  
  TString sDir(dir);
  TString objName=GetName();
  Int_t type=0;
  
  //get options
  TObjArray *arr=sDir.Tokenize(",");
  TIter next(arr);
  TObjString *s;
  while ( (s=(TObjString*)next()) ){
    TString optString=s->GetString();
    optString.Remove(TString::kBoth,' ');
    if (optString.BeginsWith("name=")){
      objName=optString.ReplaceAll("name=","");
    }
    if (optString.BeginsWith("type=")){
      optString.ReplaceAll("type=","");
      type=optString.Atoi();
    }
  }

  if (type==1||type==2) {
    //delete most probably not needed stuff
    //This requires Analyse to be called after reading the object from file
    fCalRocArrayT0.Delete();
    fCalRocArrayT0Err.Delete();
    fCalRocArrayQ.Delete();
    fCalRocArrayRMS.Delete();
    fCalRocArrayOutliers.Delete();
  }
  if (type==2||type==3){
    fParamArrayEventPol1.Delete();
    fParamArrayEventPol2.Delete();
  }
  
  TObjArray histoQArray(72);
  TObjArray histoT0Array(72);
  TObjArray histoRMSArray(72);
  TObjArray arrHnDrift(fArrHnDrift.GetEntries());

  //save all large 2D histograms in separte pointers
  //to have a smaller memory print when saving the object
  if (type==1||type==2||type==3){
    for (Int_t i=0; i<72; ++i){
      histoQArray.AddAt(fHistoQArray.UncheckedAt(i),i);
      histoT0Array.AddAt(fHistoT0Array.UncheckedAt(i),i);
      histoRMSArray.AddAt(fHistoRMSArray.UncheckedAt(i),i);
    }
    fHistoQArray.SetOwner(kFALSE);
    fHistoT0Array.SetOwner(kFALSE);
    fHistoRMSArray.SetOwner(kFALSE);
    fHistoQArray.Clear();
    fHistoT0Array.Clear();
    fHistoRMSArray.Clear();
    
    for (Int_t i=0;i<fArrHnDrift.GetEntries();++i){
      arrHnDrift.AddAt(fArrHnDrift.UncheckedAt(i),i);
    }
    fArrHnDrift.SetOwner(kFALSE);
    fArrHnDrift.Clear();
  }
  
  
  TDirectory *backup = gDirectory;
  
  TFile f(filename,"recreate");
  this->Write(objName.Data());
  if (type==1||type==2) {
    histoQArray.Write("histoQArray",TObject::kSingleKey);
    histoT0Array.Write("histoT0Array",TObject::kSingleKey);
    histoRMSArray.Write("histoRMSArray",TObject::kSingleKey);
    arrHnDrift.Write("arrHnDrift",TObject::kSingleKey);
  }

  f.Save();
  f.Close();

  //move histograms back to the object
  if (type==1||type==2){
    for (Int_t i=0; i<72; ++i){
      fHistoQArray.AddAt(histoQArray.UncheckedAt(i),i);
      fHistoT0Array.AddAt(histoT0Array.UncheckedAt(i),i);
      fHistoRMSArray.AddAt(histoRMSArray.UncheckedAt(i),i);
    }
    fHistoQArray.SetOwner(kTRUE);
    fHistoT0Array.SetOwner(kTRUE);
    fHistoRMSArray.SetOwner(kTRUE);

    for (Int_t i=0;i<arrHnDrift.GetEntries();++i){
      fArrHnDrift.AddAt(arrHnDrift.UncheckedAt(i),i);
    }
    fArrHnDrift.SetOwner(kTRUE);
  }
  
  if ( backup ) backup->cd();
}
//_____________________________________________________________________
AliTPCCalibCE* AliTPCCalibCE::ReadFromFile(const Char_t *filename)
{
  //
  // Read object from file
  // Handle properly if the histogram arrays were stored separately
  // call Analyse to make sure to have the calibration relevant information in the object
  //
  
  TFile f(filename);
  if (!f.IsOpen() || f.IsZombie() ) return 0x0;
  TList *l=f.GetListOfKeys();
  TIter next(l);
  TKey *key=0x0;
  TObject *o=0x0;

  AliTPCCalibCE *ce=0x0;
  TObjArray *histoQArray=0x0;
  TObjArray *histoT0Array=0x0;
  TObjArray *histoRMSArray=0x0;
  TObjArray *arrHnDrift=0x0;
  
  while ( (key=(TKey*)next()) ){
    o=key->ReadObj();
    if ( o->IsA()==AliTPCCalibCE::Class() ){
      ce=(AliTPCCalibCE*)o;
    } else if ( o->IsA()==TObjArray::Class() ){
      TString name=key->GetName();
      if ( name=="histoQArray") histoQArray=(TObjArray*)o;
      if ( name=="histoT0Array") histoT0Array=(TObjArray*)o;
      if ( name=="histoRMSArray") histoRMSArray=(TObjArray*)o;
      if ( name=="arrHnDrift") arrHnDrift=(TObjArray*)o;
    }
  }

  if (ce){
  //move histograms back to the object
    TH1* hist=0x0;
    if (histoQArray){
      for (Int_t i=0; i<72; ++i){
        hist=(TH1*)histoQArray->UncheckedAt(i);
        if (hist) hist->SetDirectory(0x0);
        ce->fHistoQArray.AddAt(hist,i);
      }
      ce->fHistoQArray.SetOwner(kTRUE);
    }
    
    if (histoT0Array) {
      for (Int_t i=0; i<72; ++i){
        hist=(TH1*)histoT0Array->UncheckedAt(i);
        if (hist) hist->SetDirectory(0x0);
        ce->fHistoT0Array.AddAt(hist,i);
      }
      ce->fHistoT0Array.SetOwner(kTRUE);
    }
    
    if (histoRMSArray){
      for (Int_t i=0; i<72; ++i){
        hist=(TH1*)histoRMSArray->UncheckedAt(i);
        if (hist) hist->SetDirectory(0x0);
        ce->fHistoRMSArray.AddAt(hist,i);
      }
      ce->fHistoRMSArray.SetOwner(kTRUE);
    }
    
    if (arrHnDrift){
      for (Int_t i=0; i<arrHnDrift->GetEntries(); ++i){
        THnSparseI *hSparse=(THnSparseI*)arrHnDrift->UncheckedAt(i);
        ce->fArrHnDrift.AddAt(hSparse,i);
      }
    }
    
    ce->Analyse();
  }
  f.Close();
  
  return ce;
}
