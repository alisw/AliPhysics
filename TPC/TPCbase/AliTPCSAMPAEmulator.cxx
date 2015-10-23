
/////////////////////////////////////////////////////////////////////////////////////////////////////
//     Class for emulation of the SAMPA chip (SAMPA digital Chain) in C++                          //
/////////////////////////////////////////////////////////////////////////////////////////////////////


/**		@file AliTPCSAMPAEmulator.h
 *      @brief This the header File for the SAMPA class
 
 	author: marian.ivanov@cern.ch
                mesut.arslandok@cern.ch  

For the moment very limitted subsut of the functionality implemented.
Only part realted to the Digiital filtering  - SB3 and SB2

In current interface all actions are done in once DigitalFilterFloat(Int_t npoints, Double_t *dataArray,  Double_t &baseline);
Class control variables:
   Int_t fDigitFilterType;    //  to select appropraite filter BC3(0), MAFMI(1), BC3MI(2)
   Int_t fgBaselineExportType; //  to export corrected signal (0) or pedestal estimators  (1)

Interface will change soon - implementing pipe of the filters
 e.g: 
  Signal->BC3->Rounding->ZerroSupression
    https://alice.its.cern.ch/jira/browse/ATO-129
Using  that the code will be more modular.

*/



#include <AliTPCSAMPAEmulator.h>
#include <TH1F.h>
#include <TMath.h>
#include <TSystem.h>
#include <AliDAQ.h>
#include <AliRawReader.h>
#include <AliRawVEvent.h>
#include <AliRawData.h>
#include <AliRawVEquipment.h>
#include <AliRawEquipmentHeader.h>
#include <AliTPCRawStreamV3.h>
#include <TCanvas.h>
#include <AliRawDataHeader.h>
#include <AliRawDataHeaderV3.h>
#include <cassert>

/**	@brief Constroctor of SAMPA Class
 *
 *	Consturctor of SAMPA Class, some variables are set.\n
 *	The input Data is altered, so after running the complete emulation you have the
 *	SAMPA Processed Data in the Channel Pointer.\n
 */


ClassImp(AliTPCSAMPAEmulator)

Int_t AliTPCSAMPAEmulator::fgBaselineExportType=0;

AliTPCSAMPAEmulator::AliTPCSAMPAEmulator() : 
  TNamed(),
  fDigitFilterType(0),     // type of the digital filter - BC3 default
  //
  fBC3SlopeDown(0.2),    // BC3 slope down parameter
  fBC3SlopeUp(0.1),      // BC3 slope up   parameter
  fBC3Round(-1),         // Rounding error of BC3 filter 
  fBC3DiffCutMI(2.5),      // BC3 cut on the signal difference - for MI implementation
  //
  fMAFMIKernelWidth(10), // kernel width for MAF filtering 
  fMAFMIDiffCut(-1)      // cut on the diff to skip "signal"  
{
  //
  // Constructor of SAMPA Class
  //


}



/**	@brief Destructor of SAMPA Class
 *
 *	Destructor of SAMPA Class\n
 */
AliTPCSAMPAEmulator::~AliTPCSAMPAEmulator() {
  //
  // Destructor of SAMPA Class
  // 
}


Bool_t   AliTPCSAMPAEmulator::DigitalFilterFloat(Int_t npoints, Double_t *dataArray,  Double_t &baseline){
  //
  // Run digittal filter - as configured using fDigitFilterType
  // In future more sophisticated configuration of filters, pipe  of filters should be implemented
  //    to be agreed  on interface ( e.g BC3 -> Rounding -> ZerroSpuression )
  // Current implementation used to make fast studies of the indeusced noise
  //  
  // 
  if (fDigitFilterType==0) return BC3SlopeFilterFloat(npoints,dataArray,baseline);
  if (fDigitFilterType==1) return MovingAverageFilter(npoints,dataArray,baseline);
  //
  if (fDigitFilterType==2) return  BC3SlopeFilterMI(npoints,dataArray,baseline);
  //
  assert(false);
}

void AliTPCSAMPAEmulator::SetBC3Parameters(Double_t slopeDown, Double_t slopeUp, Double_t round){
  //
  // 
  //
  fBC3SlopeDown=slopeDown;
  fBC3SlopeUp=slopeUp;
  fBC3Round=round;
}


Bool_t   AliTPCSAMPAEmulator::BC3SlopeFilterFloat(Int_t npoints, Double_t *dataArray, Double_t &baseline){
  //
  // BC3 filter
  //
  return  AliTPCSAMPAEmulator::BC3SlopeFilterFloat(npoints,dataArray,  fBC3SlopeDown,fBC3SlopeUp,  fBC3Round,baseline); 
}

Bool_t   AliTPCSAMPAEmulator::BC3SlopeFilterMI(Int_t npoints, Double_t *dataArray, Double_t &baseline){
  //
  // BC3 filter
  //
  return  AliTPCSAMPAEmulator::BC3SlopeFilterMI(npoints,dataArray,  fBC3SlopeDown,fBC3SlopeUp,  fBC3Round,baseline,fBC3DiffCutMI); 
}

Bool_t  AliTPCSAMPAEmulator::BC3SlopeFilterFloat(Int_t npoints, Double_t *dataArray, Double_t slopeDown, Double_t slopeUp, Double_t round, Double_t &slopeBaseline) {
  //
  //
  //
  // BC2 filter as should be implemented in SAMPA
  // from Konstantin description
  // https://alice.its.cern.ch/jira/browse/ATO-129
  //    discussed here the root/C code snippet for the slope based filter. The filter is applied as
  //
  //     data_filtered=data-baseline(data);
  //
  //    The filter is floating point. To simulate the integer behavior of the hardware based filter, define INTFILTER and choose the slopes to be binary compatible (1.0, 0.5, 0.25 etc.).
  // #define SLOPEUP		1.0
  // #define SLOPEDOWN	2.0

  if (npoints<=1) return kFALSE;
  for (Int_t iTimeBin=0; iTimeBin<npoints; iTimeBin++){
    Double_t data=dataArray[iTimeBin];
    if (data>slopeBaseline) {
      slopeBaseline+=slopeUp;
      if (slopeBaseline>data) slopeBaseline=data;
    } else if (data<slopeBaseline) {
      slopeBaseline-=slopeDown;
      if (slopeBaseline<data) slopeBaseline=data;
    };
    if (round>0){
      //    return round(slopeBaseline);
      slopeBaseline=TMath::Nint(slopeBaseline*round)/round;
    }
    if (fgBaselineExportType==0){
      dataArray[iTimeBin]-=slopeBaseline;
      dataArray[iTimeBin]=TMath::Nint(dataArray[iTimeBin]);
    }
    if (fgBaselineExportType==1){
      dataArray[iTimeBin]=slopeBaseline;
    }
  }
  return kTRUE;
};


Bool_t  AliTPCSAMPAEmulator::BC3SlopeFilterMI(Int_t npoints, Double_t *dataArray, Double_t slopeDown, Double_t slopeUp, Double_t round, Double_t &slopeBaseline, Double_t diffCutMI) {
  //
  //
  //
  // BC2 filter as should be implemented in SAMPA
  // from Konstantin description
  // https://alice.its.cern.ch/jira/browse/ATO-129
  //    discussed here the root/C code snippet for the slope based filter. The filter is applied as
  //
  //     data_filtered=data-baseline(data);
  //
  //    The filter is floating point. To simulate the integer behavior of the hardware based filter, define INTFILTER and choose the slopes to be binary compatible (1.0, 0.5, 0.25 etc.).
  // #define SLOPEUP		1.0
  // #define SLOPEDOWN	2.0
  
  Double_t *tmpArray= new Double_t[npoints];      
  memcpy(tmpArray, dataArray, npoints*sizeof(Double_t));

  if (npoints<=1) return kFALSE;
  for (Int_t iTimeBin=0; iTimeBin<npoints; iTimeBin++){
    Double_t data=tmpArray[iTimeBin];
    Double_t maxDiff=TMath::Max( TMath::Abs(tmpArray[(iTimeBin+npoints+1)%npoints]-data),  TMath::Abs(tmpArray[(iTimeBin+npoints-1)%npoints]-data));
    if (maxDiff<=diffCutMI){
      if (data>slopeBaseline) {
	slopeBaseline+=slopeUp;
	if (slopeBaseline>data) slopeBaseline=data;
      } else if (data<slopeBaseline) {
	slopeBaseline-=slopeDown;
	if (slopeBaseline<data) slopeBaseline=data;
      };    
      if (round>0){
	//    return round(slopeBaseline);
	slopeBaseline=TMath::Nint(slopeBaseline*round)/round;
      }
    }
    if (fgBaselineExportType==0){
      dataArray[iTimeBin]-=slopeBaseline;
      dataArray[iTimeBin]=TMath::Nint(dataArray[iTimeBin]);
    }
    if (fgBaselineExportType==1){
      dataArray[iTimeBin]=slopeBaseline;
    }
  }
  delete [] tmpArray;
  return kTRUE;
};


Bool_t  AliTPCSAMPAEmulator::MovingAverageFilter(Int_t npoints, Double_t *dataArray, Double_t &baseline){
  //
  //
  //
  return MovingAverageFilter(npoints,dataArray, fMAFMIKernelWidth, fMAFMIDiffCut,fMAFMIOnlyMinima, baseline);
}

Bool_t  AliTPCSAMPAEmulator::MovingAverageFilter(Int_t npoints, Double_t *dataArray, Double_t length, Double_t skipDiff,  Bool_t onlyMinima, Double_t &baseline){
  //
  // Baseline update formula:
  //     baseline:=(baseline*length+data)/(length+1);
  //
  if (npoints<=1) return kFALSE;
  for (Int_t iTimeBin=0; iTimeBin<npoints; iTimeBin++){
    Double_t data=dataArray[iTimeBin];
    if (skipDiff>0){ // skip data point in case it is assumed to be signal ()too big diff)
      if ( (TMath::Abs(dataArray[(iTimeBin+npoints-1)%npoints]+baseline-data)>skipDiff) ||
	   (TMath::Abs(dataArray[(iTimeBin+npoints+1)%npoints]-data)>skipDiff)) {
	continue;
      }
      if ( !(dataArray[(iTimeBin+npoints-1)%npoints]+baseline>data &&TMath::Abs(dataArray[(iTimeBin+npoints+1)%npoints]>data))){
	continue;  //use only local minima
      }
    }
    baseline=(baseline*length+data)/(length+1);
    dataArray[iTimeBin]-=baseline;
    dataArray[iTimeBin]=TMath::Nint(dataArray[iTimeBin]);
    if (fgBaselineExportType==1){
      dataArray[iTimeBin]=baseline;
    }
  }
  assert(false);
}


void AliTPCSAMPAEmulator::SetMAFMIParameters(Double_t  MAFMIKernelWidth,  Double_t  MAFMIDiffCut, Bool_t onlyMinima){
  //
  //
  //
  fMAFMIKernelWidth=MAFMIKernelWidth;   // kernel width for MAF filtering 
  fMAFMIDiffCut=MAFMIDiffCut;           // cut on the diff to skip "signal"   
  fMAFMIOnlyMinima=onlyMinima;          //   
}


Bool_t  AliTPCSAMPAEmulator::ZeroSuppression(Int_t npoints, Double_t *dataArray, Double_t threshold) {
  //
  // Zero suppression - sets all values below threshold to zero.
  // Useful for evaluating the effect of zero suppression on the results.
  // pre- and post-samples may be implemented as well.
  //
  if (npoints<=1) return kFALSE;
  for (Int_t iTimeBin=0; iTimeBin<npoints; iTimeBin++){
    Double_t data=dataArray[iTimeBin];
    if (data<threshold) dataArray[iTimeBin]=0;
  }
  return kTRUE;
}; 
