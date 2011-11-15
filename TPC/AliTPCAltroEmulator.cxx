/**		@file Altro.C
 *	@brief The Altro class implements the Altro digital Chain in C++
 *
 *	This Class represents a C++ version of the ALTRO. For a complete Documentation of the Altro
 *	Look at : http://ep-ed-alice-tpc.web.cern.ch/ep-ed-alice-tpc/altro_chip.htm\n
 *	Due to the fact that the real ALTRO constantly samples in between the recorded events,
 *	it has the knowledge on what happened in the period. This affects the BSL1, TCF and BSL2 module.
 *	In the BSL1 the ALTRO follows slow baseline drifts e.g. temperature change, the TCF has a infinite
 *	(IIR Filter) memory of "old samples" i.e. a cluster at the start of a readout cycle will be treated
 *	differently, and the BSL2 has a 8 step pipeline. The ALTRO Class can't emulate this behavior,
 *	since the data is not recorded.\n
 *
 *	@author Roland Bramm
 *	@version $LastChangedRevision: 688 $
 *	@date    $LastChangedDate: 2005-12-16 14:07:11 +0100 (Fri, 16 Dec 2005) $
 *
 *	\verbinclude Altro/Altro.C.log
 *
 */

/////////////////////////////////////////////////////////////////////////////////////////////////////
//     Class for emulation of the ALTRO chip (Altro digital Chain) in C++                          //
//     Author: Roland Bramm                                                                        //
//                                                                                                 //
//     NOTE: This class has been modified to be conform with the coding conventions of the         //
//           ALICE Off-line Project. Keywords for setting the mode of BSC1 were modified           //
//           and are shown in the header file ...                                                  //
//                           Stefan Rossegger, 8th february 2008                                   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

#include <AliTPCAltroEmulator.h>
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


/**	@brief Consturctor of Altro Class
 *
 *	Consturctor of Altro Class, some variables are set.\n
 *	The input Data is altered, so after running the complete emulation you have the
 *	Altro Processed Data in the Channel Pointer.\n
 *
 *	@param timebins an <tt> int </tt> sets the length of the input Data (Channel)
 *	@param Channel an <tt> short* </tt> Pointer to a 1d Short_tArray with the input Data
 */


ClassImp(AliTPCAltroEmulator)

AliTPCAltroEmulator::AliTPCAltroEmulator(Int_t timebins, short* Channel) : 
  TNamed(),
  ftimebins(timebins),
//  fChannelIn(Channel),
  fChannelShort(Channel), 
  fADCkeep(0),     
  fOnBSL1(0), 
  fOnTCF(0),  
  fOnBSL2(0), 
  fOnClip(0), 
  fOnZSU(0),  

  fConfiguredAltro(0),   // ConfiguredAltro
  fConfiguredBSL1(0),    // ConfiguredBSL1
  fConfiguredTCF(0),     // ConfiguredTCF
  fConfiguredTCFraw(0),     // ConfiguredTCF
  fConfiguredBSL2(0),    // ConfiguredBSL2
  fConfiguredZSU(0),     // ConfiguredZSU
  fBSL1mode(0),          // BSL1mode
  fBSL1ValuePeDestal(0), // BSL1ValuePeDestal
  fBSL1PedestalMem(0),   // BSL1PedestalMem
  fBSL1polarity(0),      // BSL1polarity

  fTCFK1(0), // K1
  fTCFK2(0), // K2
  fTCFK3(0), // K3
  fTCFL1(0), // L1
  fTCFL2(0), // L2
  fTCFL3(0), // L3

  fTCFK1Int(0), // K1Int
  fTCFK2Int(0), // K2Int
  fTCFK3Int(0), // K3Int
  fTCFL1Int(0), // L1Int
  fTCFL2Int(0), // L2Int
  fTCFL3Int(0), // L3Int

  fBSL2HighThreshold(0), // BSL2HighThreshold
  fBSL2LowThreshold(0),  // BSL2LowThreshold
  fBSL2Offset(0),        // BSL2Offset
  fBSL2Presamples(0),    // BSL2Presamples(0),
  fBSL2Postsamples(0),   // BSL2Postsamples

  fZSUThreshold(0),      // ZSUThreshold

  fZSUMinSamplesaboveThreshold(0), // ZSUMinSamplesaboveThreshold
  fZSUPresamples(0),     // ZSUPresamples
  fZSUPostsamples(0),     // ZSUPostsamples
  
  fReader(0),           // for Altro Emulation on Raw Reader
  fDecoder(0),
  fRunNumber(0), 
  fDDLFolderName("./"),
  fOutputDateFileName("./tmpRaw.date"),
  fOutputRootFileName("./tmpRaw.root"),
  fIsRandom(kTRUE),
  fChannels(0),
  fCDHs(0),
  fADCs(0),
  fTrailers(0),
  fRawData(0) {
  //
  // Constructor of Altro Class
  //

  fADCkeep = new Short_t[1024]; 

  fTCFK1IntROC[0]=0; fTCFK1IntROC[1]=0; // dummy defaults
  fTCFK2IntROC[0]=0; fTCFK2IntROC[1]=0; // dummy defaults
  fTCFK3IntROC[0]=0; fTCFK3IntROC[1]=0; // dummy defaults
  fTCFL1IntROC[0]=0; fTCFL1IntROC[1]=0; // dummy defaults
  fTCFL2IntROC[0]=0; fTCFL2IntROC[1]=0; // dummy defaults
  fTCFL3IntROC[0]=0; fTCFL3IntROC[1]=0; // dummy defaults

}



/**	@brief Destructor of Altro Class
 *
 *	Destructor of Altro Class\n
 */
AliTPCAltroEmulator::~AliTPCAltroEmulator() {
  //
  // Destructor of Altro Class
  //

  //  if(fConfiguredZSU == 1)
  delete[] fADCkeep;

  delete[] fChannels;
  delete[] fCDHs    ;
  delete[] fADCs    ;
  delete[] fTrailers;
  delete[] fRawData ;
  delete fDecoder   ;

}


/**  @brief Configures which modules of the Altro should be on.
 *
 *	Configures which modules of the Altro should be on. Each of the modules
 *	which are configured to be on, have to be configured later before running
 *	the emulation!\n
 *
 *	@param ONBaselineCorrection1 an <tt> Int_t </tt> Switch (0,1) to turn on the Base Line Correction 1 (BSL1) Module
 *	@param ONTailcancellation an <tt> Int_t </tt> Switch (0,1) to turn on the Tail Cancellation Filter (TCF) Module
 *	@param ONBaselineCorrection2 an <tt> Int_t </tt> Switch (0,1) to turn on the Moving Average Filter (BSL2) Module
 *	@param ONClipping an <tt> Int_t </tt> Switch (0,1) to turn on the Clipping Module. This is not possible in the real Altro, there it is always on.
 *	@param ONZerosuppression an <tt> Int_t </tt> Switch (0,1) to turn on the Zero Suppression (ZSU) Module
 *	@param ONDataFormatting an <tt> Int_t </tt> Switch (0,1) to turn on the Data Formatting on (not implemented)
 */
void AliTPCAltroEmulator::ConfigAltro(Int_t ONBaselineCorrection1, Int_t ONTailcancellation, Int_t ONBaselineCorrection2, Int_t ONClipping, Int_t ONZerosuppression, Int_t ONDataFormatting){
  //
  // Configures which modules of the Altro should be on 
  //
  fOnBSL1 = InRange(ONBaselineCorrection1,0,1,"AliTPCAltroEmulator::ConfigAltro","ONBaselineCorrection1");
  fOnTCF  = InRange(ONTailcancellation,0,1,"AliTPCAltroEmulator::ConfigAltro","ONTailcancellation");
  fOnBSL2 = InRange(ONBaselineCorrection2,0,1,"AliTPCAltroEmulator::ConfigAltro","ONBaselineCorrection2");
  fOnClip = InRange(ONClipping,0,1,"AliTPCAltroEmulator::ConfigAltro","ONClipping");
  fOnZSU = InRange(ONZerosuppression,0,1,"AliTPCAltroEmulator::ConfigAltro","ONZerosuppression");
  fConfiguredAltro = 1;
  if (!fConfiguredAltro) { //dummy code to avoid warning
    printf("%d\n",ONDataFormatting); // does not have to be checked
  }
}

/**  @brief Configures the Base Line Correction 1 (BSL1) Module
 *
 *	Configures the Base Line Correction 1 (BSL1) Module. You dont have to build a proper pedestalMemory
 *	array, a pointer of the correct type is enough, of course you are not allowed to use Basline
 *	Correction Modes which need then the array ...\n
 *	All configurable values are "Range checked" and if out of the Range set to the nearest extreme.
 *	So the Emulation will work, but the result is maybe not the expected one.
 *
 *	@param mode an <tt> Int_t </tt> sets the mode of the Baseline Correction. See the Altro manual for a description
 *	@param ValuePeDestal an <tt> Int_t </tt> this is the baseline of the Channel.
 *	@param PedestalMem an <tt> *Int_t </tt> Pointer to a 1d Short_t Array with the pedestal memory Data
 *	@param polarity an <tt> Int_t </tt> Switch (0,1) for the polarity
 */
void AliTPCAltroEmulator::ConfigBaselineCorrection1(Int_t mode, Int_t ValuePeDestal, Int_t *PedestalMem, Int_t polarity){
  //
  // Configures the Base Line Correction 1 (BSL1) Module
  //
  fBSL1mode          = InRange(mode,0,16,"AliTPCAltroEmulator::ConfigBaselineCorrection1","mode");
  fBSL1ValuePeDestal = InRange(ValuePeDestal,0,1023,"AliTPCAltroEmulator::BaselineCorrection1","ValuePeDestal");
  fBSL1PedestalMem = PedestalMem;
  fBSL1polarity = InRange(polarity,0,1,"AliTPCAltroEmulator::BaselineCorrection1","polarity");
  fConfiguredBSL1 = 1;
}

/**  @brief Configures the Tail Cancellation Filter (TCF) Module
 *
 *	Configures the Tail Cancellation Filter (TCF) Module. You have to set the coefficients in the
 *	Integer version.\n
 *	To convert from Int_t to Float_t use (int)*(pow(2,-16)-1)
 *	To convert from Float_t to Int_t usw (Float_t)*(pow(2,16)-1)
 *	All configurable values are "Range checked" and if out of the Range set to the nearest extreme.
 *	So the Emulation will work, but the result is maybe not the expected one.
 *
 *	@param K1 an <tt> Int_t </tt> sets the K1 coeeficient of the TCF
 *	@param K2 an <tt> Int_t </tt> sets the K2 coeeficient of the TCF
 *	@param K3 an <tt> Int_t </tt> sets the K3 coeeficient of the TCF
 *	@param L1 an <tt> Int_t </tt> sets the L1 coeeficient of the TCF
 *	@param L2 an <tt> Int_t </tt> sets the L2 coeeficient of the TCF
 *	@param L3 an <tt> Int_t </tt> sets the L3 coeeficient of the TCF
 */
void AliTPCAltroEmulator::ConfigTailCancellationFilter(Int_t K1, Int_t K2, Int_t K3, Int_t L1, Int_t L2, Int_t L3){
  //
  // Configures the Tail Cancellation Filter (TCF) Module
  //
  // conf from Int_t to fp:  (int)*(pow(2,-16)-1)
  //             backway:  (Float_t)*(pow(2,16)-1)
  fTCFK1Int = InRange(K1,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K1");
  fTCFK2Int = InRange(K2,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K2");
  fTCFK3Int = InRange(K3,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K3");
  
  fTCFL1Int = InRange(L1,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L1");
  fTCFL2Int = InRange(L2,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L2");
  fTCFL3Int = InRange(L3,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L3");
  fConfiguredTCF = 1;
}
void AliTPCAltroEmulator::ConfigTailCancellationFilterForRAWfiles(const Int_t *K1, const Int_t *K2, const Int_t *K3, 
								  const Int_t *L1, const Int_t *L2, const Int_t *L3){
  //
  // Configures the Tail Cancellation Filter (TCF) Module - Different settings for IROC and OROC
  //
  // conf from Int_t to fp:  (int)*(pow(2,-16)-1)
  //             backway:  (Float_t)*(pow(2,16)-1)

  // IROC
  fTCFK1IntROC[0] = InRange(K1[0],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K1[0]");
  fTCFK2IntROC[0] = InRange(K2[0],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K2[0]");
  fTCFK3IntROC[0] = InRange(K3[0],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K3[0]");
  fTCFL1IntROC[0] = InRange(L1[0],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L1[0]");
  fTCFL2IntROC[0] = InRange(L2[0],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L2[0]");
  fTCFL3IntROC[0] = InRange(L3[0],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L3[0]");
  // OROC
  fTCFK1IntROC[1] = InRange(K1[1],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K1[1]");
  fTCFK2IntROC[1] = InRange(K2[1],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K2[1]");
  fTCFK3IntROC[1] = InRange(K3[1],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K3[1]");
  fTCFL1IntROC[1] = InRange(L1[1],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L1[1]");
  fTCFL2IntROC[1] = InRange(L2[1],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L2[1]");
  fTCFL3IntROC[1] = InRange(L3[1],0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L3[1]");


  fConfiguredTCFraw = 1;
}


/**  @brief Configures the Moving Average Filter (BSL2) Module
 *
 *	Configures the Moving Average Filter (BSL2) Module.
 *	All configurable values are "Range checked" and if out of the Range set to the nearest extreme.
 *	So the Emulation will work, but the result is maybe not the expected one.
 *
 *	@param HighThreshold an <tt> Int_t </tt> sets the high Threshold
 *	@param LowThreshold an <tt> Int_t </tt> sets the low Theshold
 *	@param Offset an <tt> Int_t </tt> sets the the offset which is added to the Signal
 *	@param Presamples an <tt> Int_t </tt> sets the number of pre samples excluded from the moving average caclulation
 *	@param Postsamples an <tt> Int_t </tt> sets the number of post samples excluded from the moving average caclulation
 */
void AliTPCAltroEmulator::ConfigBaselineCorrection2(Int_t HighThreshold, Int_t LowThreshold, Int_t Offset, Int_t Presamples, Int_t Postsamples){
  //
  // Configures the Moving Average Filter (BSL2) Module
  //
  fBSL2HighThreshold = InRange(HighThreshold,0,1023,"AliTPCAltroEmulator::ConfigBaselineCorrection2","HighThreshold");
  fBSL2LowThreshold  = InRange(LowThreshold,0,1023,"AliTPCAltroEmulator::ConfigBaselineCorrection2","LowThreshold");
  fBSL2Offset        = InRange(Offset,0,1023,"AliTPCAltroEmulator::ConfigBaselineCorrection2","Offset");
  fBSL2Presamples    = InRange(Presamples,0,3,"AliTPCAltroEmulator::ConfigBaselineCorrection2","Presamples");
  fBSL2Postsamples   = InRange(Postsamples,0,15,"AliTPCAltroEmulator::ConfigBaselineCorrection2","Postsamples");
  fConfiguredBSL2 = 1;
}

/**  @brief Configures the Zero Suppression Module (ZSU)
 *
 *	Configures the Zero Suppression Module (ZSU).
 *	All configurable values are "Range checked" and if out of the Range set to the nearest extreme.
 *	So the Emulation will work, but the result is maybe not the expected one.
 *
 *	@param Threshold an <tt> Int_t </tt> sets the Threshold
 *	@param MinSamplesaboveThreshold an <tt> Int_t </tt> sets the minimum number of samples which have to be greater than the threshold
 *	@param Presamples an <tt> Int_t </tt> sets the number of pre samples which are kept
 *	@param Postsamples an <tt> Int_t </tt> sets the number of post samples which are kept
 */
void AliTPCAltroEmulator::ConfigZerosuppression(Int_t Threshold, Int_t MinSamplesaboveThreshold, Int_t Presamples, Int_t Postsamples){
  //
  // Configures the Zero Suppression Module (ZSU)
  //
  fZSUThreshold                = InRange(Threshold,0,1023,"AliTPCAltroEmulator::BaselineCorrection1","Threshold");
  fZSUMinSamplesaboveThreshold = InRange(MinSamplesaboveThreshold,1,3,"AliTPCAltroEmulator::BaselineCorrection1","MinSamplesaboveThreshold");
  fZSUPresamples               = InRange(Presamples,0,3,"AliTPCAltroEmulator::BaselineCorrection1","Presamples");
  fZSUPostsamples              = InRange(Postsamples,0,7,"AliTPCAltroEmulator::BaselineCorrection1","Postsamples");
  
  for(Int_t i = 0; i < ftimebins; i++){
    fADCkeep[i] = 0;
  }
  fConfiguredZSU = 1;
}

/**  @brief Prints the set Parameters, if module is configured
 *
 *	Prints the set Parameters, if module is configured.
 */
void AliTPCAltroEmulator::PrintParameters(){
  //
  // Prints the set Parameters, if module is configured
  //
  cout << "+-------------------------------------------+" << endl;
  cout << "| Configured Parameters of the Altro Module |" << endl;
  cout << "+-------------------------------------------+" << endl << endl;
  
  cout << "Parameters set in the Altro Modules:" << endl << endl;
  cout << "ONBaselineCorrection1: " << fOnBSL1 << endl;
  cout << "ONTailcancellation   : " << fOnTCF << endl;
  cout << "ONBaselineCorrection2: " << fOnBSL2 << endl;
  cout << "ONClipping           : " << fOnClip << endl;
  cout << "ONZerosuppression    : " << fOnZSU << endl << endl << endl;
  if(fConfiguredBSL1 == 1){
    cout << "Parameters set in the BSL1 (Baseline Correction 1) Module:" << endl << endl;
    cout << "mode                 : " << fBSL1mode << endl;
    cout << "ValuePeDestal        : " << fBSL1ValuePeDestal << endl;
    cout << "polarity             : " << fBSL1ValuePeDestal << endl << endl << endl;
  }else{
    cout << "BSL1 (Baseline Correction 1) Module not configured!" << endl << endl << endl;
  }
  if(fConfiguredTCF == 1){
    cout << "Parameters set in the TCF (TailCancellation Filter) Module:" << endl << endl;
    cout << "K1       (int|Float_t) : " << fTCFK1Int << " | " << fTCFK1Int/(Float_t)((1<<16)-1) << endl;
    cout << "K2       (int|Float_t) : " << fTCFK2Int << " | " << fTCFK2Int/(Float_t)((1<<16)-1) << endl;
    cout << "K3       (int|Float_t) : " << fTCFK3Int << " | " << fTCFK3Int/(Float_t)((1<<16)-1) << endl;
    cout << "L1       (int|Float_t) : " << fTCFL1Int << " | " << fTCFL1Int/(Float_t)((1<<16)-1) << endl;
    cout << "L2       (int|Float_t) : " << fTCFL2Int << " | " << fTCFL2Int/(Float_t)((1<<16)-1) << endl;
    cout << "L3       (int|Float_t) : " << fTCFL3Int << " | " << fTCFL3Int/(Float_t)((1<<16)-1) << endl << endl << endl;
  }else{
    cout << "TCF (TailCancellation Filter) Module not configured!" << endl << endl << endl;
  }
  if(fConfiguredBSL2 == 1){
    cout << "Parameters set in the BSL2 (Baseline Correction 2) Module:" << endl << endl;
    cout << "HighThreshold        : " << fBSL2HighThreshold << endl;
    cout << "LowThreshold         : " << fBSL2LowThreshold << endl;
    cout << "Offset               : " << fBSL2Offset << endl;
    cout << "Presamples           : " << fBSL2Presamples << endl;
    cout << "Postsamples          : " << fBSL2Postsamples << endl << endl << endl;
  }else{
    cout << "BSL2 (Baseline Correction 2) Module not configured!" << endl << endl << endl;
  }
  if(fConfiguredZSU == 1){
    cout << "Parameters set in the ZSU (Zero Suppression Unit) Module:" << endl << endl;
    cout << "Threshold            : " << fZSUThreshold << endl;
    cout << "MinSampaboveThreshold: " << fZSUMinSamplesaboveThreshold << endl;
    cout << "Presamples           : " << fZSUPresamples << endl;
    cout << "Postsamples          : " << fZSUPostsamples << endl << endl << endl;
  }else{
    cout << "ZSU (Zero Suppression Unit) Module not configured!" << endl << endl << endl;
  }
}


void AliTPCAltroEmulator::SetChannelData(Int_t timebins, Short_t* channelData) {
  //
  // Set channel data, for example a new channel
  //

  ftimebins = timebins;
  fChannelShort = channelData;

}
 

/**  @brief Runs the emulation of all configured Modules.
 *
 *	Runs the emulation of all configured Modules. This changes then the content of the
 *	input Array
 */
void AliTPCAltroEmulator::RunEmulation(Int_t roc){
  //
  // Runs the emulation of all configured Modules.
  //

  if (!fChannelShort) {
    printf("ERROR cant run Altro Emulation: Channel input not set.\nUse for example: SetChannelData(Int_t timebins, Short_t* Channel)\n");
    return;
  }

  //cout << "AliTPCAltroEmulator::RunEmulation | start" << endl;
  if(fConfiguredAltro == 0){
    cout << "ERROR cant run Altro Emulation because not configured" << endl;
    return;
  }
  
  //cout << "AliTPCAltroEmulator::RunEmulation | start BSL1 on: " << fOnBSL1 << " configures: " << fConfiguredBSL1 << endl;
  if(fOnBSL1 == 1){
    if(fConfiguredBSL1 == 1){
      BaselineCorrection1(fBSL1mode, fBSL1ValuePeDestal, fBSL1PedestalMem, fBSL1polarity);
    }else{
      cout << "ERROR cant run Baseline Correction 1 because not configured" << endl;
      return;
    }
  }
  
  //cout << "AliTPCAltroEmulator::RunEmulation | start TCF on: " << fOnTCF << " configures: " << fConfiguredTCF << endl;
  if(fOnTCF == 1){
    if (roc==-1) { // use one set of TCF params
      if(fConfiguredTCF == 1){
	TailCancellationFilterFixedPoint(fTCFK1Int, fTCFK2Int, fTCFK3Int, fTCFL1Int, fTCFL2Int, fTCFL3Int);
      }else{
	cout << "ERROR cant run Tail Cancellation Filter because not configured" << endl;
	return;
      }
    } else { // use different TCF params for IROC and OROC
      if(fConfiguredTCFraw == 1){
	if (roc==0)      //IROC
	  TailCancellationFilterFixedPoint(fTCFK1IntROC[0], fTCFK2IntROC[0], fTCFK3IntROC[0], 
					   fTCFL1IntROC[0], fTCFL2IntROC[0], fTCFL3IntROC[0]);
	else if (roc==1) // OROC
	  TailCancellationFilterFixedPoint(fTCFK1IntROC[1], fTCFK2IntROC[1], fTCFK3IntROC[1], 
					   fTCFL1IntROC[1], fTCFL2IntROC[1], fTCFL3IntROC[1]);
	else
	  cout << "ERROR cant run Tail Cancellation Filter because TCF settings for ROC not found" << endl;
      } else {
	cout << "ERROR cant run Tail Cancellation Filter because not configured (for RAW data files!)" << endl;
	return;
      }

    }
  }
  
  //cout << "AliTPCAltroEmulator::RunEmulation | start BSL2 on: " << fOnBSL2 << " configures: " << fConfiguredBSL2 << endl;
  if(fOnBSL2 == 1){
    if(fConfiguredBSL2 == 1){
      BaselineCorrection2RTL(fBSL2HighThreshold, fBSL2LowThreshold, fBSL2Offset, fBSL2Presamples, fBSL2Postsamples);
    }else{
      cout << "ERROR cant run Baseline Correction 2 because not configured" << endl;
      return;
    }
  }
  //cout << "AliTPCAltroEmulator::RunEmulation | start CLIP on: " << fOnClip << endl;
  if(fOnClip == 1){
    Clipping();
  }
  //cout << "AliTPCAltroEmulator::RunEmulation | start ZSU on: " << fOnZSU << " configures: " << fConfiguredZSU << endl;
  if(fOnZSU == 1){
    if(fConfiguredZSU == 1){
      Zerosuppression(fZSUThreshold,fZSUMinSamplesaboveThreshold,fZSUPresamples,fZSUPostsamples);
    }else{
      cout << "ERROR cant run Zero Suppression Unit because not configured" << endl;
      return;
    }
  }

  

}

void AliTPCAltroEmulator::BaselineCorrection1(Int_t mode, Int_t ValuePeDestal, Int_t *PedestalMem, Int_t polarity){
  //
  // BaselineCorrection1
  //

  //VPD == 0 !!
  Int_t fixedPeDestal = 0;

  // take first and last bins to calculate a mean pedestal value
  Int_t window = 3;
  Int_t meanPeDestal = 0;
  if (mode == kDINxMPD && ftimebins>=6) {
    for(Int_t i = 0; i < window; i++) {
      meanPeDestal += fChannelShort[i];
      meanPeDestal += fChannelShort[ftimebins-1-i];
    }
    meanPeDestal /= (window*2);
  }
    
  if(polarity ==1){
    for(Int_t i = 0; i < ftimebins; i++){
      fChannelShort[i]  = 1023 - fChannelShort[i];
    }
  }
  
  switch(mode) {
  case kDINxFPD:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - fixedPeDestal;
    break;
  case kDINxFT:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - PedestalMem[i];
    break;
  case kDINxFDIN:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - PedestalMem[ fChannelShort[i] ];
    break;
  case kDINxFDINxVPD:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - PedestalMem[ fChannelShort[i] - ValuePeDestal];
    break;
  case kDINxVPDxFPD:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - ValuePeDestal - fixedPeDestal;
    break;
  case kDINxVPDxFT:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - ValuePeDestal - PedestalMem[i];
    break;
  case kDINxVPDxFDIN:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - ValuePeDestal - PedestalMem[ fChannelShort[i] ];
    break;
  case kDINxVPDxFDINxVPD:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - ValuePeDestal - PedestalMem[ fChannelShort[i] - ValuePeDestal ];
    break;
  case kFDINxFPD:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = PedestalMem[ fChannelShort[i] ] - fixedPeDestal;
    break;
  case kFDINxVPDxFPD:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = PedestalMem[ fChannelShort[i] - ValuePeDestal ] - fixedPeDestal;
    break;
  case kFTxFPD:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = PedestalMem[i] - fixedPeDestal;
    break;
  case kDINxMPD:
    for(Int_t i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i] - meanPeDestal;
    break;
  }
}

Int_t AliTPCAltroEmulator::Multiply36(Int_t P, Int_t N){
  //
  // multiply function to emulate the 36 bit fixed point multiplication of the Altro.
  //
  long long retval =0;
  long long temp = 0;
  long long vAX = 0;
  temp = (long long)P*(long long)N;
  vAX = (( Mask(temp,35,18) + ((long long)(-P)<<18) ) + Mask(temp,17,0));
  if ( Maskandshift(N,17,17) == 1){
    retval = ((Maskandshift(vAX,35,35)<<17) + Maskandshift(vAX,32,16));
  }else{
    retval = Maskandshift(temp,32,16);
  }
  return retval;
}
long long AliTPCAltroEmulator::Mask(long long in, Int_t left, Int_t right){
  //
  // mask
  //
  long long retval;
  long long pattern;
  long long length = abs(left - right)+1;
  pattern = ((((long long)1)<<length)-1)<<right;
  retval = in&pattern;
  return retval;
}

long long AliTPCAltroEmulator::Maskandshift(long long in, Int_t left, Int_t right){
  //
  // maskandshift
  //
  long long retval;
  long long pattern;
  long long length = abs(left - right)+1;
  pattern = ((((long long)1)<<length)-1);
  retval = (in>>right)&pattern;
  return retval;
}

void AliTPCAltroEmulator::TailCancellationFilterFixedPoint(Int_t K1, Int_t K2, Int_t K3, Int_t L1, Int_t L2, Int_t L3){
  //
  // TailCancellationFilterFixedPoint
  //
  Int_t c1n = 0, c2n = 0, c3n = 0;
  Int_t c1o = 0, c2o = 0, c3o = 0;
  Int_t d1  = 0, d2  = 0;
  Int_t dout = 0;
  Int_t din = 0;
  Int_t bit = 0;

  //  printf("%5d %5d %5d %5d %5d %5d\n",K1,K2,K3,L1,L2,L3);

  for(Int_t i = 0; i < ftimebins; i++){
    din = fChannelShort[i];
    
    din = (din<<2);
    c1n             = Mask( (Mask(din,17,0) + Multiply36(K1,Mask(c1o,17,0)) ) ,17,0);
    d1              = Mask( (Mask(c1n,17,0) - Multiply36(L1,Mask(c1o,17,0)) ) ,17,0);
    //d1              = Mask( (Mask(c1n,17,0) + Mask(~Multiply36(L1,Mask(c1o,17,0))+1,17,0) ) ,17,0);
    
    c2n             = Mask( (Mask(d1 ,17,0) + Multiply36(K2,Mask(c2o,17,0)) ) ,17,0);
    d2              = Mask( (Mask(c2n,17,0) - Multiply36(L2,Mask(c2o,17,0)) ) ,17,0);
    //d2              = Mask( (Mask(c2n,17,0) + Mask(~Multiply36(L2,Mask(c2o,17,0))+1,17,0) ) ,17,0);
    
    c3n             = Mask( (Mask(d2 ,17,0) + Multiply36(K3,Mask(c3o,17,0)) ) ,17,0);
    dout            = Mask( (Mask(c3n,17,0) - Multiply36(L3,Mask(c3o,17,0)) ) ,17,0);
    //dout            = Mask( (Mask(c3n,17,0) + Mask(~Multiply36(L3,Mask(c3o,17,0))+1,17,0) ) ,17,0);
    
    if( (Maskandshift(dout,2,2) == 1) || (Maskandshift(dout,1,1) == 1)){
      bit = 1;
    }else{
      bit = 0;
    }
    
    dout = ((dout>>3)<<1) + bit;
    if(Maskandshift(dout,15,15) == 1){
      //is needed to get the correct coding when getting negative results
      dout = -Mask((-Mask(dout,9,0)),9,0);
    }else{
      dout = Mask(dout,9,0);
    }
    
    fChannelShort[i] = (short) dout;
    c1o = c1n;
    c2o = c2n;
    c3o = c3n;
  }
}

void AliTPCAltroEmulator::BaselineCorrection2RTL(Int_t HighThreshold, Int_t LowThreshold, Int_t Offset, Int_t Presamples, Int_t Postsamples){
  //
  // BaselineCorrection2RTL
  //

  //cout << "Altro::BaselineCorrection2RTL | HighThreshold: " << HighThreshold << " LowThreshold: " << LowThreshold << " Offset: " << Offset << " Presamples: " << Presamples << " Postsamples: " << Postsamples << endl;
  //more or less direct "translation" of the hdl code.
  //Input signals
  Int_t din;
  Int_t dout;
  Int_t edges[6]; // = Postsamples*4 + Presamples;
  Int_t offset = Offset;
  Int_t thrlo = LowThreshold;//called thr_mau[19] ...
  Int_t thrhi = HighThreshold;
  
  // Variables
  Int_t fOld[4]; //flag pipe
  Int_t fNew[4]; //flag pipe
  Int_t dOld[4]; //data pipe
  Int_t dNew[4]; //data pipe
  Int_t dxOld;
  Int_t dxNew;
  Int_t pstscnt; // Counter for Postsamples
  Int_t zOld[9]; // Filter stages
  Int_t zNew[9]; // Filter stages
  Int_t zxOld; //Accumulator stage
  Int_t zxNew; //Accumulator stage
  Int_t valcntOld; //Valid sample counter
  Int_t valcntNew = 0; //Valid sample counter
  
  Int_t valid; //Valid flag
  Int_t fx; //postsample flag
  //Int_t s07; // differentiator result
  Int_t s8; // Acc + Diff result
  Int_t flag;
  //Int_t bsth; //baseline threshold
  //Int_t din_p; //Data input strictly positive
  Int_t bsl;
  //Int_t dx_bsls; // dx -bsl
  //Int_t dx_clip; // dxbsl clipped
  //Int_t bsl_of = 0;
  
  //initialisation
  for(Int_t i = 0; i < 9 ; i++)
    zOld[i] = 0;
  for(Int_t i = 0; i < 4 ; i++){
    fOld[i] = 0;
    dOld[i] = 0;
  }
  dxOld= 0;
  pstscnt = 0;
  zxOld = 0;
  valcntOld = 0;
  valid = 0;
  for(Int_t i = 0; i < 2 ; i++){
    edges[i] = (Presamples&(1<<i))>>i;
  }
  for(Int_t i = 0; i < 4 ; i++){
    edges[(3-i)+2] = (Postsamples&(1<<i))>>i;
  }
  /*cout << "edges :";
    for(Int_t i = 0; i < 6 ; i++)
    cout << edges[i] << ":";
    cout << " Presamples: " << Presamples << " Postsamples: " << Postsamples << endl;*/
  
  //Loop
  //cout << "AliTPCAltroEmulator::BaselineCorrection2_RTL | starting Loop" << endl;
  for(Int_t timebin = -12; timebin < ftimebins+10; timebin++){
    //cout << "AliTPCAltroEmulator::BaselineCorrection2_RTL | in Loop timebin: " << timebin << endl;
    din = GetElement(fChannelShort,timebin);
    
    s8 = zxOld + (zOld[8] - zOld[0]);
    
    if(valid == 1)
      bsl = s8>>3;// ...
    else
      bsl = 0;
    
    //assign flag = (din_p > thrhi) | (thrlo > din_p);	// Signal samples between thresholds
    if( (din <= (bsl + thrhi)) && (din >= (bsl - thrlo)) )
      flag = 0;
    else
      flag = 1;
    
    if(pstscnt == 0)
      fx = 0;
    else
      fx = 1;
    
    if(valcntOld >= 12)
      valid = 1;
    else
      valid = 0;
    
    fNew[3] = flag;
    
    if( (fOld[3] == 1) || ( (flag == 1) && ( (edges[0] == 1) || (edges[1] == 1) ) ) ) //f[2] =  f[3] | (flag&(edges[0]|edges[1]));
      fNew[2] = 1;
    else
      fNew[2] = 0;
    
    if( (fOld[2] == 1) || ( (edges[1] == 1) && (flag == 1) ) ) //		f[1] =  f[2] | (edges[1] & flag);
      fNew[1] = 1;
    else
      fNew[1] = 0;
    
    if( ( (fOld[1] == 1) || ( (flag == 1) && (edges[0] == 1) && (edges[1] == 1) )  || (fx==1) ) && (valid==1) ) //		f[0] = (f[1] | (edges[1] & edges[0] & flag) | fx) & valid;
      fNew[0] = 1;
    else
      fNew[0] = 0;
    
    dxNew = dOld[0];
    for(Int_t i = 0; i < 3; i++)
      dNew[i] = dOld[i+1];
    dNew[3] = din;
    
    if( (fOld[1]==1) && (fOld[2]==0) )
      pstscnt = Postsamples;
    else if(fx == 1)
      pstscnt--;
    
    if(fOld[0] == 0){
      if(valid == 0)
	valcntNew =  ++valcntOld;
      
      zxNew = s8;
      for(Int_t i = 0; i < 8; i++)
	zNew[i] = zOld[i+1];
      zNew[8] = dOld[0];
    }else{
      zxNew = zxOld;
      for(Int_t i = 0; i < 9; i++)
	zNew[i] = zOld[i];
    }
    dout = dxOld - (bsl - offset);
    //if(dout <0)
    //	dout = 0;
    
    SetElement(fChannelShort,timebin-5,(short)dout);
    //sim clockschange
    for(Int_t i = 0; i < 9 ; i++)
      zOld[i] = zNew[i];
    zxOld = zxNew;
    for(Int_t i = 0; i < 4 ; i++){
      fOld[i] = fNew[i];
      dOld[i] = dNew[i];
    }
    dxOld = dxNew;
    valcntOld = valcntNew;
  }
}

void AliTPCAltroEmulator::Clipping(){ 
  //
  // implement if no BC2 clipping has to run
  //
  for(Int_t i = 0; i < ftimebins; i++){
    if(fChannelShort[i] < -1)
      fChannelShort[i] = -1;
  }
}

void AliTPCAltroEmulator::Zerosuppression(Int_t Threshold, Int_t MinSamplesaboveThreshold, Int_t Presamples, Int_t Postsamples){
  //
  // add again altro feature
  //

  //TODO: Implement "Altro zsu merging"
  //Int_t Postsamplecounter = 0;
  //Int_t setPostsample = 0;

  for(Int_t i = 0; i < ftimebins; i++){
    if(fChannelShort[i] >= Threshold)
      fADCkeep[i] = 1;
    else
      fADCkeep[i] = 0;
  }

  Int_t startofclustersequence = -1;
  Int_t endofClustersInSequence = -1;
  
  for(Int_t i = 0; i < ftimebins; i++){
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i-1) == 0) ){
      startofclustersequence = i;
    }
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i+1) == 0) ){
      endofClustersInSequence = i;
    }
    //cout << i << " startofclustersequence: " << startofclustersequence << " endofClustersInSequence: " << endofClustersInSequence;
    if( (startofclustersequence != -1) && (endofClustersInSequence != -1) ){
      //cout << " found! " <<  (endofClustersInSequence - startofclustersequence + 1);
      if ( (endofClustersInSequence - startofclustersequence + 1) < MinSamplesaboveThreshold ){
	for(Int_t j = startofclustersequence; j <= endofClustersInSequence ; j++){
	  fADCkeep[j] = 0;
	}
      }
      startofclustersequence = -1;
      endofClustersInSequence = -1;
    }
    //cout << endl;
  }
  
  /*for(Int_t i = 0; i < ftimebins; i++){
    if( (GetElement(fADCkeep,i-1) == 1) && (GetElement(fADCkeep,i) == 0) && (GetElement(fADCkeep,i+1) == 1) ){
    SetElement(fADCkeep,i,1);
    }
    }*/
  
  for(Int_t i = 0; i < ftimebins; i++){
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i-1) == 0) ){
      for(Int_t j = i-Presamples ; j <= i; j++){
	SetElement(fADCkeep,j,1);
      }
    }
  }
  for(Int_t i = ftimebins; i >= 0; i--){
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i+1) == 0) ){
      for(Int_t j = i ; j <= i+Postsamples; j++){
	SetElement(fADCkeep,j,1);
      }
    }
  }
  /*cout << " Postsamplecounter: " << Postsamplecounter;
    for(Int_t j = i+1 ; j <= i+Postsamples; j++){
    SetElement(fADCkeep,j,1);
    i+=Postsamples;
    }
    cout << endl;
    }
    cout << i << " ADCK: " << GetElement(fADCkeep,i);
    cout << " Postsam: " << Postsamplecounter << " ADCK: " << GetElement(fADCkeep,i);*/
  
  for(Int_t i = 0; i < ftimebins; i++){
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i+1) == 0) && ( (GetElement(fADCkeep,i+3) == 1) || (GetElement(fADCkeep,i+2) == 1) ) ){
      SetElement(fADCkeep,i+1,1);
      SetElement(fADCkeep,i+2,1);
    }
  }

  for(Int_t i = 0; i < ftimebins; i++){
    if( !GetKeepChannel(i) ) {
      SetElement(fChannelShort,i,-1); // set non relevant data to -1
    }
  }


}

/**  @brief formats the data like the ALTRO. Result is a 64 bit array
 *
 *	formats the data like the ALTRO. Result is a 64 bit array
 *
 */

void AliTPCAltroEmulator::DataFormater(){
  //
  // formats the data like the ALTRO. Result is a 64 bit array
  //


}


/**  @brief calculates the compression out of the bitmask
 *
 *	calculates the compression out of the bitmask with the set adc values
 *
 * 	@return \c Float_t consisting of the compression factor
 */
Float_t AliTPCAltroEmulator::CalculateCompression() {
  //
  // calculates the compression out of the bitmask
  //

  // calculation is based on altro 10 bit words ..
  Int_t sample = 0;
  Int_t cluster = 0;
  Int_t data = 0;
  Float_t retval = 0.0;
  
  for(Int_t i = 0; i < ftimebins; i++){
    if(fADCkeep[i] == 1){
      sample++;
    }
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i+1) == 0) ){
      cluster++;
    }
  }
  data = sample + cluster*2;
  data = data + data%4 + 4;
  if(data >0){
    retval = ftimebins / (Float_t)data;//num of timebins is equal to max number of samples
  }else{
    retval = 1.0;
  }
  return retval;
}

Short_t AliTPCAltroEmulator::GetElement(short* Array,Int_t index){
  //
  // GetElement of array
  //
  if (index < 0)
    return 0;
  else if(index >= ftimebins)
    return 0;
  else
    return Array[index];
}

void AliTPCAltroEmulator::SetElement(short* Array,Int_t index,Short_t value){
  //
  // SetElement of array
  //
  if (index < 0)
    return;
  else if(index >= ftimebins)
    return;
  else
    Array[index] = value;
}

Int_t AliTPCAltroEmulator::InBand(Int_t ADC,Int_t bsl, Int_t LowThreshold, Int_t HighThreshold){
  //
  // check if it's within the band of search
  //
  Int_t fLow = bsl - LowThreshold;
  Int_t fHigh = bsl + HighThreshold;
  if( (ADC <= fHigh) && (ADC >= fLow) )
    return 1;
  else
    return 0;
}

Int_t AliTPCAltroEmulator::InRange(Int_t parameter,Int_t Low,Int_t High,const char *Module,const char *ParameterName){
  //
  // checks it it's within the range
  //

  char out[255];
  Int_t retval;
  if(parameter > High){
    snprintf(out,255,"Error | %s | Parameter %s is to big, has to be %d <= %s <= %d, is %d, now set to %d",Module,ParameterName,Low,ParameterName,High,parameter,High);
    cout << out << endl;
    retval = High;
  }else if(parameter < Low){
    snprintf(out,255,"Error | %s | Parameter %s is to small, has to be %d <= %s <= %d, is %d, now set to %d",Module,ParameterName,Low,ParameterName,High,parameter,Low);
    cout << out << endl;
    retval = Low;
  }else{
    retval = parameter;
  }
  return retval;
}

Short_t AliTPCAltroEmulator::GetShortChannel(Int_t i){
  //
  // GetElement of channel
  //
  return GetElement(fChannelShort,i);
}

Short_t AliTPCAltroEmulator::GetKeepChannel(Int_t i){
  //
  // GetElement of Keep Channel
  //
  return GetElement(fADCkeep,i);
}

Bool_t AliTPCAltroEmulator::WriteEvent(Int_t ievent) {
  //
  // Write event to the DDL data folders
  //

  for (Int_t ddlID=0;ddlID<216;++ddlID) {
    Bool_t  *channelsDDL=fChannels+ddlID*4096     ;
    Short_t *adcsDDL    =fADCs    +ddlID*4096*1024;
    UInt_t  *cdhDDL     =fCDHs    +ddlID*8        ;
    UInt_t  *trailerDDL =fTrailers+ddlID*9        ;
    
    FILE *file=fopen(Form("%s/raw%d/TPC_%03d.ddl",
			  fDDLFolderName.Data(),ievent,768+ddlID),
		     "wb");
    if (!file) {
      return kFALSE;
    } else { 
      Int_t i32;
      // write CDH (first word to be altered later)
      for (i32=0;i32<8;++i32)
	fRawData[i32]=cdhDDL[i32];
      
      // process payload
      for (Int_t hwaddr=0;hwaddr<4096;++hwaddr) if (channelsDDL[hwaddr]) {
	Short_t *adcsChannel=adcsDDL+hwaddr*1024;
	// merge custers
	// TODO: acqusition window
	for (Int_t it=0;it<1024-3;++it) {
	  if (adcsChannel[it]>=0&&adcsChannel[it+3]>=0) {
	    if (adcsChannel[it+1]<0) {
	      //	    printf("merge");
	      adcsChannel[it+1]=0;
	    }
	    if (adcsChannel[it+2]<0) {
	      //	    printf("merge");
	      adcsChannel[it+2]=0;
	    }
	  }
	}
	Int_t i10=3;
	Int_t icw=0;
	Int_t its=1;
	Int_t cw =0;
	Int_t ts =0;
	for (Int_t it=1023;it>=0;--it) {
	  Short_t w10=adcsChannel[it];
	  if (w10>=0) {
	    if (cw<0) {
	      icw=i10++;
	      its=i10++;
	      cw =0    ;
	      ts=it    ;
	    }
	    fRawData[i32+i10/3]|=w10<<(10*(2-i10%3));
	    ++i10;
	    ++cw;
	  }
	  else {
	    if (cw>=0) {
	      cw+=2;
	      fRawData[i32+icw/3]|=cw <<(10*(2-icw%3));
	      fRawData[i32+its/3]|=ts <<(10*(2-its%3));
	      cw=-1;
	    }
	  }
	}
	fRawData[i32]=0x1<<30|(i10-3)<<16|hwaddr;
	i32+=(i10+2)/3;
	
	// clean up
	for (Int_t i=0;i<1024;++i) adcsChannel[i]=-1;
	channelsDDL[hwaddr]=kFALSE;
      }
      
      // write RCU trailer
      fRawData[i32]=0x2<<30|(i32-8);i32++;
      for (Int_t i=0;i<8;++i)
	fRawData[i32++]=trailerDDL[i];
      
      // write first word of CDH
      fRawData[0]=i32*4;
      
      Int_t nwritten=fwrite(fRawData,sizeof(UInt_t),i32,file);
      fclose(file);

      if (nwritten!=i32) return kFALSE;
      
      // clean up
      do {fRawData[--i32]=0;} while (i32>0);
    }
  }
  return kTRUE;
}

Bool_t AliTPCAltroEmulator::GDC2DDLs(AliRawVEvent *gdc,Int_t ievent) {
  //
  // Converte GDC data to DDL format
  //
  for(Int_t iLDC=0;iLDC<gdc->GetNSubEvents();++iLDC) {
    AliRawVEvent *ldc=gdc->GetSubEvent(iLDC);
    for(Int_t iEq=0;iEq<ldc->GetNEquipments();++iEq) {
      AliRawVEquipment *eq=ldc->GetEquipment(iEq);
      AliRawEquipmentHeader *eqHeader=eq->GetEquipmentHeader();
      Int_t eqSize=eqHeader->GetEquipmentSize();
      if (eqSize>0) {
	Int_t ddlIndex;
	Int_t detId=AliDAQ::DetectorIDFromDdlID(eqHeader->GetId(),ddlIndex);
	Int_t nwritten=0;
	FILE *ddlFile=fopen(Form("%s/raw%d/%s",
				 fDDLFolderName.Data(),
				 ievent,
				 AliDAQ::DdlFileName(detId,ddlIndex)),
			    "wb");
	AliRawData *rawData=eq->GetRawData();
	if (ddlFile) {
	  nwritten=fwrite(rawData->GetBuffer(),1,rawData->GetSize(),ddlFile);
	  fclose(ddlFile);
	}
	if (nwritten<rawData->GetSize()) return kFALSE;
      }
    }
  }
  return kTRUE;
}


Bool_t AliTPCAltroEmulator::ConvertRawFilesToDate(Int_t nevents) {
  //
  //  Convertes Raw files to Date format
  //
  
  // from $ALICE_ROOT/STEER/AliSimulation.cxx
  
  char command[100];
  FILE *pipe;

  printf(" RAW to DATE CONVERSION: Run Number %d\n",fRunNumber);

  if (fRunNumber>0)
    pipe=gSystem->OpenPipe(Form("dateStream -c -s -D -o %s -C -# %d -run %d", 
				fOutputDateFileName.Data(),
				nevents,
				fRunNumber),
			   "w");
  else
    pipe=gSystem->OpenPipe(Form("dateStream -c -s -D -o %s -C -# %d", 
				fOutputDateFileName.Data(),
				nevents),
			   "w");
  if (!pipe) {
    fprintf(stderr,"error: cannot execute command: %s",command);
    return kFALSE;
  }

  for (Int_t ievent=0;ievent<nevents;++ievent) {
    UInt_t detectorPattern = 0xFFFFFFFF;
    fprintf(pipe, "GDC DetectorPattern %u\n", detectorPattern);

    Float_t ldc = 0;
    Int_t prevLDC = -1;
  
    // loop over detectors and DDLs
    for (Int_t iDet = 0; iDet < AliDAQ::kNDetectors; iDet++) {
      if (!(iDet<=5 || iDet==17 )) continue;
      for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls(iDet); iDDL++) {
	//	printf("iDet=%d, iDDL=%d, filenmae=%s\n",iDet,iDDL,filename);
	Int_t ddlID = AliDAQ::DdlID(iDet,iDDL);
	Int_t ldcID = Int_t(ldc + 0.0001);
	ldc += AliDAQ::NumberOfLdcs(iDet) / AliDAQ::NumberOfDdls(iDet);

	// check existence and size of raw data file
	FILE* file = fopen(Form("%s/raw%d/%s",
				fDDLFolderName.Data(),
				ievent,
				AliDAQ::DdlFileName(iDet,iDDL)),
			   "rb");
	if (!file) continue;
	fseek(file, 0, SEEK_END);
	unsigned long size = ftell(file);
	fclose(file);
	if (!size) continue;

	if (ldcID != prevLDC) {
	  fprintf(pipe, " LDC Id %d\n", ldcID);
	  prevLDC = ldcID;
   	}
	fprintf(pipe,Form("  Equipment Id %d Payload %s/raw%d/%s\n",
			  ddlID,
			  fDDLFolderName.Data(),
			  ievent,
			  (char*)AliDAQ::DdlFileName(iDet,iDDL))
		);
      }
    }
  }
  Int_t result = gSystem->ClosePipe(pipe);
  return (result == 0);
}

void AliTPCAltroEmulator::InitBuffers() {
  // 
  // Initialization of the Buffers
  //
  if (!fChannels) fChannels=new Bool_t [216*4096     ];
  if (!fCDHs    ) fCDHs    =new UInt_t [216*8        ];
  if (!fADCs    ) fADCs    =new Short_t[216*4096*1024];
  if (!fTrailers) fTrailers=new UInt_t [216*9        ];
  if (!fRawData ) fRawData =new UInt_t [  1*4096*1024]; // be save...

  for (Int_t i=0;i<216*4096     ;++i) fChannels[i]=kFALSE;
  // no need to init CDHs
  for (Int_t i=0;i<216*4096*1024;++i) fADCs    [i]=-1    ;
  // no need to init trailers
  for (Int_t i=0;i<  1*4096*1024;++i) fRawData [i]= 0    ;
}

Bool_t AliTPCAltroEmulator::ConvertDateToRoot() {
  //
  // convert a DATE file to a root file with the program "alimdc"
  //

  // from $ALICE_ROOT/STEER/AliSimulation.cxx

  printf(" DATE to ROOT CONVERSION: Run Number %d\n",fRunNumber);

  // ALIMDC setup
  const Int_t kDBSize    = 2000000000; //2GB
  const Int_t kTagDBSize = 1000000000;
  const Bool_t kFilter = kFALSE;
  const Int_t kCompression = 1;

  //  AliInfo(Form("converting DATE file %s to root file %s", 
  //               dateFileName, rootFileName));

  const char* rawDBFS[2] = { "/tmp/mdc1", "/tmp/mdc2" };
  const char* tagDBFS    = "/tmp/mdc1/tags";

  // User defined file system locations
  if (gSystem->Getenv("ALIMDC_RAWDB1")) 
    rawDBFS[0] = gSystem->Getenv("ALIMDC_RAWDB1");
  if (gSystem->Getenv("ALIMDC_RAWDB2")) 
    rawDBFS[1] = gSystem->Getenv("ALIMDC_RAWDB2");
  if (gSystem->Getenv("ALIMDC_TAGDB")) 
    tagDBFS = gSystem->Getenv("ALIMDC_TAGDB");

  gSystem->Exec(Form("rm -rf %s",rawDBFS[0]));
  gSystem->Exec(Form("rm -rf %s",rawDBFS[1]));
  gSystem->Exec(Form("rm -rf %s",tagDBFS));

  gSystem->Exec(Form("mkdir %s",rawDBFS[0]));
  gSystem->Exec(Form("mkdir %s",rawDBFS[1]));
  gSystem->Exec(Form("mkdir %s",tagDBFS));

  Int_t result = gSystem->Exec(Form("alimdc %d %d %d %d %s", 
				    kDBSize, kTagDBSize, kFilter, kCompression, fOutputDateFileName.Data()));
  gSystem->Exec(Form("mv %s/*.root %s", rawDBFS[0],fOutputRootFileName.Data()));

  gSystem->Exec(Form("rm -rf %s",rawDBFS[0]));
  gSystem->Exec(Form("rm -rf %s",rawDBFS[1]));
  gSystem->Exec(Form("rm -rf %s",tagDBFS));

  return (result == 0);
}

void AliTPCAltroEmulator::RunEmulationOnRAWdata(AliRawReader *reader, Int_t plotFlag) {
  //
  // Run the Altro Emulation on a full Raw data set (AliRawReader)
  // plus write the outcome in a RAW data format
  //

  if (!reader) {
    printf("ERROR cant run Altro Emulation: AliRawReader is zero. No RAW data file.\n");
    return;
  }

  fReader=reader;
  if (fDecoder) delete fDecoder;
  fDecoder=new AliTPCRawStreamV3(reader);

  InitBuffers();
  
  Int_t chanCount=0;
  TH1F hisO("DINO","DINO",1024,0,1024); 
  TH1F his("DIN","DIN",1024,0,1024); 
  his.GetYaxis()->SetRangeUser(-20,90);
  Short_t *data = new Short_t[1024]; 
  TCanvas *c1 =0;
  if (plotFlag) {
    c1 = new TCanvas("c1","c1");
    c1->SetGridx(); c1->SetGridy();
  }

  // event loop
  Int_t ievent=0;
  while (fReader->NextEvent()) {
    
    if (fReader->GetRunNumber()>0) fRunNumber=fReader->GetRunNumber();
    gSystem->Exec(Form("mkdir -p %s/raw%d/",fDDLFolderName.Data(),ievent));
    GDC2DDLs(const_cast<AliRawVEvent*>(fReader->GetEvent()),ievent);
    
    Int_t ddlC =0;
    while (fDecoder->NextDDL()) {
      Int_t ddlID=fDecoder->GetDDLNumber();
      printf("ddl: %d (%d/216)\n",ddlID,++ddlC);
     
      Bool_t  *channelsDDL=fChannels+ddlID*4096      ;
      Short_t *adcsDDL    =fADCs    +ddlID*4096*1024 ;
      UInt_t  *cdhDDL     =fCDHs    +ddlID*8         ;
      UInt_t  *trailerDDL =fTrailers+ddlID*9         ;
      
      // CDH 
      for (Int_t i=0;i<8;++i)
      	// just to show how ugly it is...
      	cdhDDL[i]=reinterpret_cast<UInt_t*>(const_cast<AliRawDataHeader*>(fReader->GetDataHeader()))[i]; 
      
      // PAYLOAD
      while (fDecoder->NextChannel()) {
	Int_t hwaddr=fDecoder->GetHWAddress();
	Int_t sector=fDecoder->GetSector();
	Int_t row=fDecoder->GetRow();
	Int_t pad=fDecoder->GetPad();
	Short_t *adcsChannel=adcsDDL+hwaddr*1024;
	while (fDecoder->NextBunch()) {
	  UInt_t          ts     =fDecoder->GetStartTimeBin();
	  Int_t           cw     =fDecoder->GetBunchLength() ;
	  const UShort_t *signals=fDecoder->GetSignals()     ;
	  for (Int_t ci=0;ci<cw;++ci) {
	    Short_t s=signals[ci];
	    Int_t   t=ts-ci;
	    // TODO aqcuisition window
	    if (0<=t&&t<(1024-3)) {
	      channelsDDL[hwaddr]=kTRUE;
	      if (adcsChannel[t]<0)
		adcsChannel[t]=s;
	      else
		adcsChannel[t]+=s;
	      if (adcsChannel[t]>0x3ff)
		adcsChannel[t]=0x3ff;
	    }
	  }
	}

	// search start of aquisition
	Int_t t0 = 0; 	 while (adcsChannel[t0]==-1) t0++;
	// search end of aquisition
	Int_t tE = 1024; while (adcsChannel[tE]==-1) tE--;
	
	// channel is complete - Perform Altro Emulation
	if (plotFlag!=0 && !(chanCount%plotFlag) ) {
	  for (Int_t t=0; t<1024; t++) {
	    his.SetBinContent(t+1,adcsChannel[t]);
	  }
	  his.SetTitle(Form("sig_sec%d_row%d_pad%d",sector,row,pad));
	  his.SetStats(0); his.GetXaxis()->SetTitle("timebin"); 
	  his.DrawCopy();
	}

	// FEED THE ALTRO EMULATOR WITH CLEAN SIGNAL (no aquisition-window ghosts)
	Int_t timebins = tE-t0+1;
	for (Int_t t=t0;t<(t0+timebins);t++)
	  data[t-t0]=adcsChannel[t];
	SetChannelData(timebins,data);
	
	Int_t roc = (sector%36)>=18; // 0 for IROC, 1 for OROC
	RunEmulation(roc); 
	for (Int_t t=t0;t<(t0+timebins);t++) 
	  adcsChannel[t]=data[t-t0];
	  
	if (plotFlag!=0 && !(chanCount%plotFlag) ) {
	  for (Int_t t=0; t<1024; t++)
	    hisO.SetBinContent(t+1,adcsChannel[t]);
	  hisO.SetStats(0); hisO.SetLineColor(2);
	  hisO.DrawCopy("same");
	    
	  c1->SaveAs(Form("%s/sig_sec%02d_row%02d_pad%03d_%s_%d%d%d%d%d.png",
			  fDDLFolderName.Data(),sector,row,pad,
			  fOutputRootFileName.Data(),
			  fOnBSL1,fOnTCF,fOnBSL2,fOnClip,fOnZSU));
	    
	  his.Reset();
	  hisO.Reset();

	}
	chanCount++;

      }
      
      // TRAILER 
      UChar_t *rcuTrailer;
      fDecoder->GetRCUTrailerData(rcuTrailer);
      for (Int_t i=0;i<= /* (!) */ fDecoder->GetRCUTrailerSize()/4;++i)
      	trailerDDL[i]=reinterpret_cast<UInt_t*>(rcuTrailer)[i]; // again: UGLY!
      
    }
            
    if (ddlC>0) WriteEvent(ievent++);
  }

  delete[] data; // free space

  // convert to date and back
  ConvertRawFilesToDate(ievent);
  ConvertDateToRoot();


}
