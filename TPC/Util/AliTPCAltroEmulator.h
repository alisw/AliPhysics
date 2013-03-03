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
#include <AliRawReader.h>
#include <AliTPCRawStreamV3.h>


using namespace std;

class AliTPCAltroEmulator : public TNamed {

 public:
  AliTPCAltroEmulator(Int_t timebins=0, Short_t* Channel=0);

  ~AliTPCAltroEmulator();

  void ConfigAltro(Int_t ONBaselineCorrection1, Int_t ONTailcancellation, Int_t ONBaselineCorrection2, Int_t ONClipping, Int_t ONZerosuppression, Int_t ONDataFormatting);
  void ConfigBaselineCorrection1(Int_t mode, Int_t ValuePeDestal, Int_t *PedestalMem, Int_t polarity);
  void ConfigTailCancellationFilter(Int_t K1, Int_t K2, Int_t K3, Int_t L1, Int_t L2, Int_t L3);
  void ConfigTailCancellationFilterForRAWfiles(const Int_t* K1, const Int_t* K2, const Int_t* K3, 
					       const Int_t* L1, const Int_t* L2, const Int_t* L3);
  void ConfigBaselineCorrection2(Int_t HighThreshold, Int_t LowThreshold, Int_t Offset, Int_t Presamples, Int_t Postsamples);
  void ConfigZerosuppression(Int_t Threshold, Int_t MinSamplesaboveThreshold, Int_t Presamples, Int_t Postsamples);

  void SetChannelData(Int_t timebins, Short_t* Channel);
  void PrintParameters();
  void RunEmulation(Int_t roc=-1); // if -1, the standard "single" TCF is used
  Float_t CalculateCompression();

  // perform altro emulation on raw-reader level

  void RunEmulationOnRAWdata(AliRawReader *reader, Int_t plotFlag=0);

  TString GetDDLFolderName()      const {return fDDLFolderName     ;}
  TString GetOutputDateFileName() const {return fOutputDateFileName;}
  TString GetOutputRootFileName() const {return fOutputRootFileName;}
  void SetDDLFolderName     (const TString &name) {fDDLFolderName     =name;}
  void SetOutputDateFileName(const TString &name) {fOutputDateFileName=name;}
  void SetOutputRootFileName(const TString &name) {fOutputRootFileName=name;}




  enum {
    /**din - fpd*/			kDINxFPD,
    /**din - f(t)*/			kDINxFT,
    /**din - f(din)*/			kDINxFDIN,
    /**din - f(din-vpd)*/		kDINxFDINxVPD,
    /**din - vpd - fpd*/		kDINxVPDxFPD,
    /**din - vpd - f(t)*/		kDINxVPDxFT,
    /**din - vpd - f(din)*/		kDINxVPDxFDIN,
    /**din - vpd - f(din - vpd)*/	kDINxVPDxFDINxVPD,
    /**f(din) - fpd*/			kFDINxFPD,
    /**f(din - vpd) - fpd*/		kFDINxVPDxFPD,
    /**f(t) - fpd*/			kFTxFPD,
    /**f(t) - f(t)*/			kFTxFT,
    /**f(din) - f(din)*/		kFDINxFDIN,
    /**f(din - vpd) - f(din - vpd)*/    kFDINxVPDxFDINxVPD,
    /**din - fpd*/			kDINxFPD1,
    /**din - fpd*/			kDINxFPD2,
    /** 16. din-mean*/                  kDINxMPD
  };

 private:

  AliTPCAltroEmulator(const AliTPCAltroEmulator &sig);
  AliTPCAltroEmulator& operator = (const  AliTPCAltroEmulator &source);

  Int_t ftimebins;          // timebins

  //	Short_t *fChannelIn;      // ChannelIn
  Short_t *fChannelShort;   // incoming signal in Short_t format
  Short_t *fADCkeep;        // ADCkeep

  Int_t fOnBSL1;            // Baseline correction and substraction 1 on
  Int_t fOnTCF;             // Tail Cancelation Filter on
  Int_t fOnBSL2;            // Baseline correction and substraction 2 (MAF) on
  Int_t fOnClip;            // Clipping on (to reverse the signal for ZSU if BSL2 is on)
  Int_t fOnZSU;             // Zero Suppression on

  Int_t fConfiguredAltro;   // ConfiguredAltro
  Int_t fConfiguredBSL1;    // ConfiguredBSL1
  Int_t fConfiguredTCF;     // ConfiguredTCF
  Int_t fConfiguredTCFraw;  // ConfiguredTCF for RAW data files
  Int_t fConfiguredBSL2;    // ConfiguredBSL2
  Int_t fConfiguredZSU;     // ConfiguredZSU

  Int_t fBSL1mode;          // BSL1mode
  Int_t fBSL1ValuePeDestal; // BSL1ValuePeDestal
  Int_t* fBSL1PedestalMem;  // BSL1PedestalMem
  Int_t fBSL1polarity;      // BSL1polarity

  Float_t fTCFK1; // K1
  Float_t fTCFK2; // K2
  Float_t fTCFK3; // K3
  Float_t fTCFL1; // L1
  Float_t fTCFL2; // L2
  Float_t fTCFL3; // L3

  Int_t fTCFK1Int; // K1Int
  Int_t fTCFK2Int; // K2Int
  Int_t fTCFK3Int; // K3Int
  Int_t fTCFL1Int; // L1Int
  Int_t fTCFL2Int; // L2Int
  Int_t fTCFL3Int; // L3Int

  Int_t fTCFK1IntROC[2]; // K1Int (IROC/OROC)
  Int_t fTCFK2IntROC[2]; // K2Int (IROC/OROC)
  Int_t fTCFK3IntROC[2]; // K3Int (IROC/OROC)
  Int_t fTCFL1IntROC[2]; // L1Int (IROC/OROC)
  Int_t fTCFL2IntROC[2]; // L2Int (IROC/OROC)
  Int_t fTCFL3IntROC[2]; // L3Int (IROC/OROC)

  Int_t fBSL2HighThreshold; // BSL2HighThreshold
  Int_t fBSL2LowThreshold;  // BSL2LowThreshold
  Int_t fBSL2Offset;        // BSL2Offset
  Int_t fBSL2Presamples;    // BSL2Presamples;
  Int_t fBSL2Postsamples;   // BSL2Postsamples

  Int_t fZSUThreshold;      // ZSUThreshold
  Int_t fZSUMinSamplesaboveThreshold; // ZSUMinSamplesaboveThreshold
  Int_t fZSUPresamples;     // ZSUPresamples
  Int_t fZSUPostsamples;    // ZSUPostsamples

  void BaselineCorrection1(Int_t mode, Int_t FixedPeDestal, Int_t *PedestalMem, Int_t polarity);
  void TailCancellationFilterFixedPoint(Int_t K1, Int_t K2, Int_t K3, Int_t L1, Int_t L2, Int_t L3);
  void BaselineCorrection2RTL(Int_t HighThreshold, Int_t LowThreshold, Int_t Offset, Int_t Presamples, Int_t Postsamples);
  void Clipping();
  void Zerosuppression(Int_t Threshold, Int_t MinSamplesaboveThreshold, Int_t Presamples, Int_t Postsamples);
  void DataFormater();

  Short_t GetElement(short* Array,Int_t index);
  void SetElement(short* Array,Int_t index,Short_t value);

  Int_t InBand(Int_t ADC,Int_t bsl, Int_t LowThreshold, Int_t HighThreshold);
  Int_t InRange(Int_t parameter,Int_t Low,Int_t High,const char *Module,const char *ParameterName);
  Short_t GetShortChannel(Int_t i);
  Short_t GetKeepChannel(Int_t i);
  Int_t Multiply36(Int_t P, Int_t N);
  long long Mask(long long in, Int_t left, Int_t right);
  long long Maskandshift(long long in, Int_t left, Int_t right);

  
  
  void InitBuffers();
  Bool_t AddEvent(Int_t dt,Bool_t isFirst);
  Bool_t CreateEvent(Int_t ievent);
  Bool_t GDC2DDLs(AliRawVEvent *gdc,Int_t ievent);
  Bool_t ConvertRawFilesToDate(Int_t nevents);
  Bool_t ConvertDateToRoot();
  Bool_t WriteEvent(Int_t ievent);

  AliRawReader      *fReader ; // RAW reader
  AliTPCRawStreamV3 *fDecoder; // ALTRO decoder

  Int_t fRunNumber;            // Run Number

  TString fDDLFolderName;      // folder name for ddl files
  TString fOutputDateFileName; // filename for date output
  TString fOutputRootFileName; // filename for root output

  //  Float_t  fP[2047] ; // Interaction probabilities for times (T-1023,...T,...T+1023)
  Bool_t   fIsRandom; // Indicates if fP are treated as probabilities (in terms of Possionian statistics), or fixed numbers
  Bool_t  *fChannels; //! field of active channels
  UInt_t  *fCDHs    ; //! CDHs
  Short_t *fADCs    ; //! field of ADC counts
  UInt_t  *fTrailers; //! RCU trailers
  UInt_t  *fRawData ; //! Raw Data


  ClassDef(AliTPCAltroEmulator,1);
};
#endif
