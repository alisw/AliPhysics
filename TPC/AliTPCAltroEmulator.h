#ifndef ALI_TPC_ALTRO_EMULATOR_H
#define ALI_TPC_ALTRO_EMULATOR_H


/**		@file AliTPCAltroEmulator.h
	*	@brief This the header File for the Altro class
	*
	*	@author Roland Bramm
	*	@version $LastChangedRevision: 688 $
	*	@date    $LastChangedDate: 2005-12-16 14:07:11 +0100 (Fri, 16 Dec 2005) $
	*
	*	\verbinclude Altro/Altro.h.log
*/


///////////////////////////////////////////////////////////////////////////////
//                        Class AliTPCAltroEmulator                          //
//  Class for emulation of the ALTRO chip (Altro digital Chain) in C++       //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"


#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace std;

class AliTPCAltroEmulator : public TNamed {

public:
	AliTPCAltroEmulator(int timebins, short* Channel);
	AliTPCAltroEmulator(const AliTPCAltroEmulator &sig);
	~AliTPCAltroEmulator();
	AliTPCAltroEmulator& operator = (const  AliTPCAltroEmulator &source);

	void ConfigAltro(int ONBaselineCorrection1, int ONTailcancellation, int ONBaselineCorrection2, int ONClipping, int ONZerosuppression, int ONDataFormatting);
	void ConfigBaselineCorrection1(int mode, int ValuePeDestal, int *PedestalMem, int polarity);
	void ConfigTailCancellationFilter(int K1, int K2, int K3, int L1, int L2, int L3);
	void ConfigBaselineCorrection2(int HighThreshold, int LowThreshold, int Offset, int Presamples, int Postsamples);
	void ConfigZerosuppression(int Threshold, int MinSamplesaboveThreshold, int Presamples, int Postsamples);
	void PrintParameters();
	void RunEmulation();
	float CalculateCompression();

	enum {
		/**din - fpd*/				kDINxFPD,
		/**din - f(t)*/				kDINxFT,
		/**din - f(din)*/			kDINxFDIN,
		/**din - f(din-vpd)*/			kDINxFDINxVPD,
		/**din - vpd - fpd*/			kDINxVPDxFPD,
		/**din - vpd - f(t)*/			kDINxVPDxFT,
		/**din - vpd - f(din)*/			kDINxVPDxFDIN,
		/**din - vpd - f(din - vpd)*/	        kDINxVPDxFDINxVPD,
		/**f(din) - fpd*/			kFDINxFPD,
		/**f(din - vpd) - fpd*/			kFDINxVPDxFPD,
		/**f(t) - fpd*/				kFTxFPD,
		/**f(t) - f(t)*/			kFTxFT,
		/**f(din) - f(din)*/			kFDINxFDIN,
		/**f(din - vpd) - f(din - vpd)*/        kFDINxVPDxFDINxVPD,
		/**din - fpd*/				kDINxFPD1,
		/**din - fpd*/				kDINxFPD2
	};

	private:
	int ftimebins;          // timebins

	//	short *fChannelIn;      // ChannelIn
	short *fChannelShort;   // incoming signal in short format
	short *fADCkeep;        // ADCkeep

	int fOnBSL1;            // Baseline correction and substraction 1 on
	int fOnTCF;             // Tail Cancelation Filter on
	int fOnBSL2;            // Baseline correction and substraction 2 (MAF) on
	int fOnClip;            // Clipping on (to reverse the signal for ZSU if BSL2 is on)
	int fOnZSU;             // Zero Suppression on

	int fConfiguredAltro;   // ConfiguredAltro
	int fConfiguredBSL1;    // ConfiguredBSL1
	int fConfiguredTCF;     // ConfiguredTCF
	int fConfiguredBSL2;    // ConfiguredBSL2
	int fConfiguredZSU;     // ConfiguredZSU

	int fBSL1mode;          // BSL1mode
	int fBSL1ValuePeDestal; // BSL1ValuePeDestal
	int* fBSL1PedestalMem;  // BSL1PedestalMem
	int fBSL1polarity;      // BSL1polarity

	float fTCFK1; // K1
	float fTCFK2; // K2
	float fTCFK3; // K3
	float fTCFL1; // L1
	float fTCFL2; // L2
	float fTCFL3; // L3

	int fTCFK1Int; // K1Int
	int fTCFK2Int; // K2Int
	int fTCFK3Int; // K3Int
	int fTCFL1Int; // L1Int
	int fTCFL2Int; // L2Int
	int fTCFL3Int; // L3Int

	int fBSL2HighThreshold; // BSL2HighThreshold
	int fBSL2LowThreshold;  // BSL2LowThreshold
	int fBSL2Offset;        // BSL2Offset
	int fBSL2Presamples;    // BSL2Presamples;
	int fBSL2Postsamples;   // BSL2Postsamples

	int fZSUThreshold;      // ZSUThreshold
	int fZSUMinSamplesaboveThreshold; // ZSUMinSamplesaboveThreshold
	int fZSUPresamples;     // ZSUPresamples
	int fZSUPostsamples;    // ZSUPostsamples

	void BaselineCorrection1(int mode, int FixedPeDestal, int *PedestalMem, int polarity);
	void TailCancellationFilterFixedPoint(int K1, int K2, int K3, int L1, int L2, int L3);
	void BaselineCorrection2RTL(int HighThreshold, int LowThreshold, int Offset, int Presamples, int Postsamples);
	void Clipping();
	void Zerosuppression(int Threshold, int MinSamplesaboveThreshold, int Presamples, int Postsamples);
	void DataFormater();

	short GetElement(short* Array,int index);
	void SetElement(short* Array,int index,short value);

	int InBand(int ADC,int bsl, int LowThreshold, int HighThreshold);
	int InRange(int parameter,int Low,int High,const char *Module,const char *ParameterName);
	short GetShortChannel(int i);
	short GetKeepChannel(int i);
	int Multiply36(int P, int N);
	long long Mask(long long in, int left, int right);
	long long Maskandshift(long long in, int left, int right);
  ClassDef(AliTPCAltroEmulator,0)
};
#endif
