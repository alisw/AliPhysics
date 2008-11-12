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

#include "AliTPCAltroEmulator.h"

/**	@brief Consturctor of Altro Class
 *
 *	Consturctor of Altro Class, some variables are set.\n
 *	The input Data is altered, so after running the complete emulation you have the
 *	Altro Processed Data in the Channel Pointer.\n
 *
 *	@param timebins an <tt> int </tt> sets the length of the input Data (Channel)
 *	@param Channel an <tt> short* </tt> Pointer to a 1d short Array with the input Data
 */


ClassImp(AliTPCAltroEmulator)

AliTPCAltroEmulator::AliTPCAltroEmulator(int timebins, short* Channel) : 
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
  fConfiguredBSL2(0),    // ConfiguredBSL2
  fConfiguredZSU(0),     // ConfiguredZSU
  fBSL1mode(0),          // BSL1mode
  fBSL1ValuePeDestal(0), // BSL1ValuePeDestal
  fBSL1PedestalMem(0),  // BSL1PedestalMem
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
  fZSUPostsamples(0)     // ZSUPostsamples

{
  //
  // Constructor of Altro Class
  //
  /*
  ftimebins = timebins;
  
  fChannelShort = Channel;
  
  fOnBSL1 = 0;
  fOnTCF = 0;
  fOnBSL2 = 0;
  fOnClip = 0;
  fOnZSU = 0;
  
  fConfiguredAltro = 0;
  fConfiguredBSL1 = 0;
  fConfiguredTCF = 0;
  fConfiguredBSL2 = 0;
  fConfiguredZSU = 0;
  */
}


AliTPCAltroEmulator::AliTPCAltroEmulator(const AliTPCAltroEmulator &altro):
  TNamed(),
  ftimebins(altro.ftimebins),
//  fChannelIn(Channel),
  fChannelShort(altro.fChannelShort), 
  fADCkeep(0),     
  fOnBSL1(0), 
  fOnTCF(0),  
  fOnBSL2(0), 
  fOnClip(0), 
  fOnZSU(0),  

  fConfiguredAltro(0),   // ConfiguredAltro
  fConfiguredBSL1(0),    // ConfiguredBSL1
  fConfiguredTCF(0),     // ConfiguredTCF
  fConfiguredBSL2(0),    // ConfiguredBSL2
  fConfiguredZSU(0),     // ConfiguredZSU
  fBSL1mode(0),          // BSL1mode
  fBSL1ValuePeDestal(0), // BSL1ValuePeDestal
  fBSL1PedestalMem(0),  // BSL1PedestalMem
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
  fZSUPostsamples(0)     // ZSUPostsamples

{
  //
  // copy constructor of Altro Class
  //

}


/**	@brief Destructor of Altro Class
 *
 *	Destructor of Altro Class\n
 */
AliTPCAltroEmulator::~AliTPCAltroEmulator(){
  //
  // Destructor of Altro Class
  //

  if(fConfiguredZSU == 1)
    delete fADCkeep;
}

//_____________________________________________________________________________
AliTPCAltroEmulator& AliTPCAltroEmulator::operator = (const AliTPCAltroEmulator &source)
{
  //
  // AliTPCAltroEmulator assignment operator
  //

  if (&source == this) return *this;
  new (this) AliTPCAltroEmulator(source);

  return *this;

}



/**  @brief Configures which modules of the Altro should be on.
 *
 *	Configures which modules of the Altro should be on. Each of the modules
 *	which are configured to be on, have to be configured later before running
 *	the emulation!\n
 *
 *	@param ONBaselineCorrection1 an <tt> int </tt> Switch (0,1) to turn on the Base Line Correction 1 (BSL1) Module
 *	@param ONTailcancellation an <tt> int </tt> Switch (0,1) to turn on the Tail Cancellation Filter (TCF) Module
 *	@param ONBaselineCorrection2 an <tt> int </tt> Switch (0,1) to turn on the Moving Average Filter (BSL2) Module
 *	@param ONClipping an <tt> int </tt> Switch (0,1) to turn on the Clipping Module. This is not possible in the real Altro, there it is always on.
 *	@param ONZerosuppression an <tt> int </tt> Switch (0,1) to turn on the Zero Suppression (ZSU) Module
 *	@param ONDataFormatting an <tt> int </tt> Switch (0,1) to turn on the Data Formatting on (not implemented)
 */
void AliTPCAltroEmulator::ConfigAltro(int ONBaselineCorrection1, int ONTailcancellation, int ONBaselineCorrection2, int ONClipping, int ONZerosuppression, int ONDataFormatting){
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
 *	@param mode an <tt> int </tt> sets the mode of the Baseline Correction. See the Altro manual for a description
 *	@param ValuePeDestal an <tt> int </tt> this is the baseline of the Channel.
 *	@param PedestalMem an <tt> *int </tt> Pointer to a 1d short Array with the pedestal memory Data
 *	@param polarity an <tt> int </tt> Switch (0,1) for the polarity
 */
void AliTPCAltroEmulator::ConfigBaselineCorrection1(int mode, int ValuePeDestal, int *PedestalMem, int polarity){
  //
  // Configures the Base Line Correction 1 (BSL1) Module
  //
  fBSL1mode          = InRange(mode,0,10,"AliTPCAltroEmulator::ConfigBaselineCorrection1","mode");
  fBSL1ValuePeDestal = InRange(ValuePeDestal,0,1023,"AliTPCAltroEmulator::BaselineCorrection1","ValuePeDestal");
  fBSL1PedestalMem = PedestalMem;
  fBSL1polarity = InRange(polarity,0,1,"AliTPCAltroEmulator::BaselineCorrection1","polarity");
  fConfiguredBSL1 = 1;
}

/**  @brief Configures the Tail Cancellation Filter (TCF) Module
 *
 *	Configures the Tail Cancellation Filter (TCF) Module. You have to set the coefficients in the
 *	Integer version.\n
 *	To convert from int to float use (int)*(pow(2,-16)-1)
 *	To convert from float to int usw (float)*(pow(2,16)-1)
 *	All configurable values are "Range checked" and if out of the Range set to the nearest extreme.
 *	So the Emulation will work, but the result is maybe not the expected one.
 *
 *	@param K1 an <tt> int </tt> sets the K1 coeeficient of the TCF
 *	@param K2 an <tt> int </tt> sets the K2 coeeficient of the TCF
 *	@param K3 an <tt> int </tt> sets the K3 coeeficient of the TCF
 *	@param L1 an <tt> int </tt> sets the L1 coeeficient of the TCF
 *	@param L2 an <tt> int </tt> sets the L2 coeeficient of the TCF
 *	@param L3 an <tt> int </tt> sets the L3 coeeficient of the TCF
 */
void AliTPCAltroEmulator::ConfigTailCancellationFilter(int K1, int K2, int K3, int L1, int L2, int L3){
  //
  // Configures the Tail Cancellation Filter (TCF) Module
  //
  // conf from int to fp:  (int)*(pow(2,-16)-1)
  //             backway:  (float)*(pow(2,16)-1)
  fTCFK1Int = InRange(K1,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K1");
  fTCFK2Int = InRange(K2,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K2");
  fTCFK3Int = InRange(K3,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","K3");
  
  fTCFL1Int = InRange(L1,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L1");
  fTCFL2Int = InRange(L2,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L2");
  fTCFL3Int = InRange(L3,0,65535,"AliTPCAltroEmulator::ConfigTailCancellationFilter","L3");
  fConfiguredTCF = 1;
}

/**  @brief Configures the Moving Average Filter (BSL2) Module
 *
 *	Configures the Moving Average Filter (BSL2) Module.
 *	All configurable values are "Range checked" and if out of the Range set to the nearest extreme.
 *	So the Emulation will work, but the result is maybe not the expected one.
 *
 *	@param HighThreshold an <tt> int </tt> sets the high Threshold
 *	@param LowThreshold an <tt> int </tt> sets the low Theshold
 *	@param Offset an <tt> int </tt> sets the the offset which is added to the Signal
 *	@param Presamples an <tt> int </tt> sets the number of pre samples excluded from the moving average caclulation
 *	@param Postsamples an <tt> int </tt> sets the number of post samples excluded from the moving average caclulation
 */
void AliTPCAltroEmulator::ConfigBaselineCorrection2(int HighThreshold, int LowThreshold, int Offset, int Presamples, int Postsamples){
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
 *	@param Threshold an <tt> int </tt> sets the Threshold
 *	@param MinSamplesaboveThreshold an <tt> int </tt> sets the minimum number of samples which have to be greater than the threshold
 *	@param Presamples an <tt> int </tt> sets the number of pre samples which are kept
 *	@param Postsamples an <tt> int </tt> sets the number of post samples which are kept
 */
void AliTPCAltroEmulator::ConfigZerosuppression(int Threshold, int MinSamplesaboveThreshold, int Presamples, int Postsamples){
  //
  // Configures the Zero Suppression Module (ZSU)
  //
  fZSUThreshold                = InRange(Threshold,0,1023,"AliTPCAltroEmulator::BaselineCorrection1","Threshold");
  fZSUMinSamplesaboveThreshold = InRange(MinSamplesaboveThreshold,1,3,"AliTPCAltroEmulator::BaselineCorrection1","MinSamplesaboveThreshold");
  fZSUPresamples               = InRange(Presamples,0,3,"AliTPCAltroEmulator::BaselineCorrection1","Presamples");
  fZSUPostsamples              = InRange(Postsamples,0,7,"AliTPCAltroEmulator::BaselineCorrection1","Postsamples");
  fADCkeep = (short *)calloc(sizeof(short),ftimebins);
  
  for(int i = 0; i < ftimebins; i++){
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
    cout << "K1       (int|float) : " << fTCFK1Int << " | " << fTCFK1Int/(float)((1<<16)-1) << endl;
    cout << "K2       (int|float) : " << fTCFK2Int << " | " << fTCFK2Int/(float)((1<<16)-1) << endl;
    cout << "K3       (int|float) : " << fTCFK3Int << " | " << fTCFK3Int/(float)((1<<16)-1) << endl;
    cout << "L1       (int|float) : " << fTCFL1Int << " | " << fTCFL1Int/(float)((1<<16)-1) << endl;
    cout << "L2       (int|float) : " << fTCFL2Int << " | " << fTCFL2Int/(float)((1<<16)-1) << endl;
    cout << "L3       (int|float) : " << fTCFL3Int << " | " << fTCFL3Int/(float)((1<<16)-1) << endl << endl << endl;
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

/**  @brief Runs the emulation of all configured Modules.
 *
 *	Runs the emulation of all configured Modules. This changes then the content of the
 *	input Array
 */
void AliTPCAltroEmulator::RunEmulation(){
  //
  // Runs the emulation of all configured Modules.
  //

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
    if(fConfiguredTCF == 1){
      TailCancellationFilterFixedPoint(fTCFK1Int, fTCFK2Int, fTCFK3Int, fTCFL1Int, fTCFL2Int, fTCFL3Int);
    }else{
      cout << "ERROR cant run Tail Cancellation Filter because not configured" << endl;
      return;
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

void AliTPCAltroEmulator::BaselineCorrection1(int mode, int ValuePeDestal, int *PedestalMem, int polarity){
  //
  // BaselineCorrection1
  //

  //VPD == 0 !!
  int fixedPeDestal = 0;
  
  if(polarity ==1){
    for(int i = 0; i < ftimebins; i++){
      fChannelShort[i]  = 1023 - fChannelShort[i];
    }
  }
  
  switch(mode) {
  case kDINxFPD:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - fixedPeDestal;
    break;
  case kDINxFT:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - PedestalMem[i];
    break;
  case kDINxFDIN:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - PedestalMem[ fChannelShort[i] ];
    break;
  case kDINxFDINxVPD:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - PedestalMem[ fChannelShort[i] - ValuePeDestal];
    break;
  case kDINxVPDxFPD:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - ValuePeDestal - fixedPeDestal;
    break;
  case kDINxVPDxFT:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - ValuePeDestal - PedestalMem[i];
    break;
  case kDINxVPDxFDIN:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - ValuePeDestal - PedestalMem[ fChannelShort[i] ];
    break;
  case kDINxVPDxFDINxVPD:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = fChannelShort[i]  - ValuePeDestal - PedestalMem[ fChannelShort[i] - ValuePeDestal ];
    break;
  case kFDINxFPD:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = PedestalMem[ fChannelShort[i] ] - fixedPeDestal;
    break;
  case kFDINxVPDxFPD:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = PedestalMem[ fChannelShort[i] - ValuePeDestal ] - fixedPeDestal;
    break;
  case kFTxFPD:
    for(int i = 0; i < ftimebins; i++)
      fChannelShort[i]  = PedestalMem[i] - fixedPeDestal;
    break;
  }
}

int AliTPCAltroEmulator::Multiply36(int P, int N){
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
long long AliTPCAltroEmulator::Mask(long long in, int left, int right){
  //
  //
  //
  long long retval;
  long long pattern;
  long long length = abs(left - right)+1;
  pattern = ((1<<length)-1)<<right;
  retval = in&pattern;
  return retval;
}

long long AliTPCAltroEmulator::Maskandshift(long long in, int left, int right){
  //
  //
  //
  long long retval;
  long long pattern;
  long long length = abs(left - right)+1;
  pattern = ((1<<length)-1);
  retval = (in>>right)&pattern;
  return retval;
}

void AliTPCAltroEmulator::TailCancellationFilterFixedPoint(int K1, int K2, int K3, int L1, int L2, int L3){
  //
  // TailCancellationFilterFixedPoint
  //
  int c1n = 0, c2n = 0, c3n = 0;
  int c1o = 0, c2o = 0, c3o = 0;
  int d1  = 0, d2  = 0;
  int dout = 0;
  int din = 0;
  int bit = 0;
  for(int i = 0; i < ftimebins; i++){
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

void AliTPCAltroEmulator::BaselineCorrection2RTL(int HighThreshold, int LowThreshold, int Offset, int Presamples, int Postsamples){
  //
  // BaselineCorrection2RTL
  //

  //cout << "Altro::BaselineCorrection2RTL | HighThreshold: " << HighThreshold << " LowThreshold: " << LowThreshold << " Offset: " << Offset << " Presamples: " << Presamples << " Postsamples: " << Postsamples << endl;
  //more or less direct "translation" of the hdl code.
  //Input signals
  int din;
  int dout;
  int edges[6]; // = Postsamples*4 + Presamples;
  int offset = Offset;
  int thrlo = LowThreshold;//called thr_mau[19] ...
  int thrhi = HighThreshold;
  
  // Variables
  int fOld[4]; //flag pipe
  int fNew[4]; //flag pipe
  int dOld[4]; //data pipe
  int dNew[4]; //data pipe
  int dxOld;
  int dxNew;
  int pstscnt; // Counter for Postsamples
  int zOld[9]; // Filter stages
  int zNew[9]; // Filter stages
  int zxOld; //Accumulator stage
  int zxNew; //Accumulator stage
  int valcntOld; //Valid sample counter
  int valcntNew = 0; //Valid sample counter
  
  int valid; //Valid flag
  int fx; //postsample flag
  //int s07; // differentiator result
  int s8; // Acc + Diff result
  int flag;
  //int bsth; //baseline threshold
  //int din_p; //Data input strictly positive
  int bsl;
  //int dx_bsls; // dx -bsl
  //int dx_clip; // dxbsl clipped
  //int bsl_of = 0;
  
  //initialisation
  for(int i = 0; i < 9 ; i++)
    zOld[i] = 0;
  for(int i = 0; i < 4 ; i++){
    fOld[i] = 0;
    dOld[i] = 0;
  }
  dxOld= 0;
  pstscnt = 0;
  zxOld = 0;
  valcntOld = 0;
  valid = 0;
  for(int i = 0; i < 2 ; i++){
    edges[i] = (Presamples&(1<<i))>>i;
  }
  for(int i = 0; i < 4 ; i++){
    edges[(3-i)+2] = (Postsamples&(1<<i))>>i;
  }
  /*cout << "edges :";
    for(int i = 0; i < 6 ; i++)
    cout << edges[i] << ":";
    cout << " Presamples: " << Presamples << " Postsamples: " << Postsamples << endl;*/
  
  //Loop
  //cout << "AliTPCAltroEmulator::BaselineCorrection2_RTL | starting Loop" << endl;
  for(int timebin = -12; timebin < ftimebins+10; timebin++){
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
    for(int i = 0; i < 3; i++)
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
      for(int i = 0; i < 8; i++)
	zNew[i] = zOld[i+1];
      zNew[8] = dOld[0];
    }else{
      zxNew = zxOld;
      for(int i = 0; i < 9; i++)
	zNew[i] = zOld[i];
    }
    dout = dxOld - (bsl - offset);
    //if(dout <0)
    //	dout = 0;
    
    SetElement(fChannelShort,timebin-5,(short)dout);
    //sim clockschange
    for(int i = 0; i < 9 ; i++)
      zOld[i] = zNew[i];
    zxOld = zxNew;
    for(int i = 0; i < 4 ; i++){
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
  for(int i = 0; i < ftimebins; i++){
    if(fChannelShort[i] < 0)
      fChannelShort[i] = 0;
  }
}

void AliTPCAltroEmulator::Zerosuppression(int Threshold, int MinSamplesaboveThreshold, int Presamples, int Postsamples){
  //
  // add again altro feature
  //

  //TODO: Implement "Altro zsu merging"
  //int Postsamplecounter = 0;
  //int setPostsample = 0;
  for(int i = 0; i < ftimebins; i++){
    if(fChannelShort[i] >= Threshold){
      fADCkeep[i] = 1;
    }
  }
  
  int startofclustersequence = -1;
  int endofClustersInSequence = -1;
  
  for(int i = 0; i < ftimebins; i++){
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
	for(int j = startofclustersequence; j <= endofClustersInSequence ; j++){
	  fADCkeep[j] = 0;
	}
      }
      startofclustersequence = -1;
      endofClustersInSequence = -1;
    }
    //cout << endl;
  }
  
  /*for(int i = 0; i < ftimebins; i++){
    if( (GetElement(fADCkeep,i-1) == 1) && (GetElement(fADCkeep,i) == 0) && (GetElement(fADCkeep,i+1) == 1) ){
    SetElement(fADCkeep,i,1);
    }
    }*/
  
  for(int i = 0; i < ftimebins; i++){
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i-1) == 0) ){
      for(int j = i-Presamples ; j <= i; j++){
	SetElement(fADCkeep,j,1);
      }
    }
  }
  for(int i = ftimebins; i >= 0; i--){
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i+1) == 0) ){
      for(int j = i ; j <= i+Postsamples; j++){
	SetElement(fADCkeep,j,1);
      }
    }
  }
  /*cout << " Postsamplecounter: " << Postsamplecounter;
    for(int j = i+1 ; j <= i+Postsamples; j++){
    SetElement(fADCkeep,j,1);
    i+=Postsamples;
    }
    cout << endl;
    }
    cout << i << " ADCK: " << GetElement(fADCkeep,i);
    cout << " Postsam: " << Postsamplecounter << " ADCK: " << GetElement(fADCkeep,i);*/
  
  for(int i = 0; i < ftimebins; i++){
    if( (fADCkeep[i] == 1) && (GetElement(fADCkeep,i+1) == 0) && ( (GetElement(fADCkeep,i+3) == 1) || (GetElement(fADCkeep,i+2) == 1) ) ){
      SetElement(fADCkeep,i+1,1);
      SetElement(fADCkeep,i+2,1);
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
 * 	@return \c float consisting of the compression factor
 */
float AliTPCAltroEmulator::CalculateCompression(){
  //
  // calculates the compression out of the bitmask
  //

  // calculation is based on altro 10 bit words ..
  int sample = 0;
  int cluster = 0;
  int data = 0;
  float retval = 0.0;
  
  for(int i = 0; i < ftimebins; i++){
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
    retval = ftimebins / (float)data;//num of timebins is equal to max number of samples
  }else{
    retval = 1.0;
  }
  return retval;
}

short AliTPCAltroEmulator::GetElement(short* Array,int index){
  //
  //
  //
  if (index < 0)
    return 0;
  else if(index >= ftimebins)
    return 0;
  else
    return Array[index];
}

void AliTPCAltroEmulator::SetElement(short* Array,int index,short value){
  //
  //
  //
  if (index < 0)
    return;
  else if(index >= ftimebins)
    return;
  else
    Array[index] = value;
}

int AliTPCAltroEmulator::InBand(int ADC,int bsl, int LowThreshold, int HighThreshold){
  //
  //
  //
  int fLow = bsl - LowThreshold;
  int fHigh = bsl + HighThreshold;
  if( (ADC <= fHigh) && (ADC >= fLow) )
    return 1;
  else
    return 0;
}

int AliTPCAltroEmulator::InRange(int parameter,int Low,int High,const char *Module,const char *ParameterName){
  //
  //
  //

  char out[255];
  int retval;
  if(parameter > High){
    sprintf(out,"Error | %s | Parameter %s is to big, has to be %d <= %s <= %d, is %d, now set to %d",Module,ParameterName,Low,ParameterName,High,parameter,High);
    cout << out << endl;
    retval = High;
  }else if(parameter < Low){
    sprintf(out,"Error | %s | Parameter %s is to small, has to be %d <= %s <= %d, is %d, now set to %d",Module,ParameterName,Low,ParameterName,High,parameter,Low);
    cout << out << endl;
    retval = Low;
  }else{
    retval = parameter;
  }
  return retval;
}

short AliTPCAltroEmulator::GetShortChannel(int i){
  //
  //
  //
  return GetElement(fChannelShort,i);
}

short AliTPCAltroEmulator::GetKeepChannel(int i){
  //
  //
  //
  return GetElement(fADCkeep,i);
}
