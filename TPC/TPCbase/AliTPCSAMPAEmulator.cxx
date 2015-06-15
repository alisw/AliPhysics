
/////////////////////////////////////////////////////////////////////////////////////////////////////
//     Class for emulation of the ALTRO chip (SAMPA digital Chain) in C++                          //
//     Author: Roland Bramm                                                                        //
/////////////////////////////////////////////////////////////////////////////////////////////////////


/**		@file AliTPCSAMPAEmulator.h
 *      @brief This the header File for the SAMPA class
 
 	author: marian.ivanov@cern.ch
                mesut.arslandok@cern.ch  
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

/**	@brief Constroctor of SAMPA Class
 *
 *	Consturctor of SAMPA Class, some variables are set.\n
 *	The input Data is altered, so after running the complete emulation you have the
 *	SAMPA Processed Data in the Channel Pointer.\n
 *
 *	@param timebins an <tt> int </tt> sets the length of the input Data (Channel)
 *	@param Channel an <tt> short* </tt> Pointer to a 1d Short_tArray with the input Data
 */


ClassImp(AliTPCSAMPAEmulator)

AliTPCSAMPAEmulator::AliTPCSAMPAEmulator() : 
  TNamed(),
  fDigitFilterType(0),   // type of the digital filter
  //
  fBC3SlopeDown(0.2),    // BC3 slope down parameter
  fBC3SlopeUp(0.1),      // BC3 slope up   parameter
  fBC3Round(-1),         // Rounding error of BC3 filter
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
  if (fDigitFilterType==0) return BC3SlopeFilterFloat(npoints,dataArray,baseline);
  if (fDigitFilterType==1) return MovingAverageFilter(npoints,dataArray,baseline);
  //
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
    dataArray[iTimeBin]-=slopeBaseline;
    dataArray[iTimeBin]=TMath::Nint(dataArray[iTimeBin]);
  }
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
  }
}


void AliTPCSAMPAEmulator::SetMAFMIParameters(Double_t  MAFMIKernelWidth,  Double_t  MAFMIDiffCut, Bool_t onlyMinima){
  //
  //
  //
  fMAFMIKernelWidth=MAFMIKernelWidth;   // kernel width for MAF filtering 
  fMAFMIDiffCut=MAFMIDiffCut;           // cut on the diff to skip "signal"   
  fMAFMIOnlyMinima=onlyMinima;          //   
}
