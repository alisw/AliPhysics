#ifndef ALICALOCALIBPEDESTAL_H
#define ALICALOCALIBPEDESTAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// \file AliCaloCalibPedestal.h
//   \brief Description:
//   A help class for monitoring and calibration tools: MOOD, AMORE etc.,
//   that can process events from a standard AliCaloRawStream,
//   most usually from LED/pulser runs. It stores signal info as
//   typical (highest) amplitude, and pedestal info in geometrically-binned
//   2D profiles of the detectors (EMCAL and PHOS).
//   Comparisons (ratios and differences) can be done with references.

//   \author: Timo Alho (Jyvaskyla), original version. 
//   [Consultant: D. Silvermyr (ORNL)]
//   Partly based on AliTPCCalibPedestal.
//   
//   \version $Revision$
//   \date $Date$

#include "TProfile2D.h"
#include "TH2.h"
#include "TObjArray.h"
class AliCaloRawStream;
class AliCaloAltroMapping;
class AliRawReader;

class AliCaloCalibPedestal : public TObject {
  
 public:

  enum kDetType {kPhos, kEmCal, kNone};//The detector types
  enum kDeadMapEntry{kAlive = 0, kDead, kResurrected, kRecentlyDeceased, kNumDeadMapStates};//The entries being put to the deadmap
  
  AliCaloCalibPedestal(kDetType detectorType = kPhos);
  virtual ~AliCaloCalibPedestal();

  // copy ctor, and '=' operator, are not fully tested/debugged yet
  // at least for now; the reference info is not copied from one to the other
  AliCaloCalibPedestal(const AliCaloCalibPedestal &ped); 
  AliCaloCalibPedestal& operator = (const  AliCaloCalibPedestal &source);
  
  //Functions to ask for the constants (in case a GUI needs them, for an example
  int GetSampleMax() const {return fgkSampleMax;};
  int GetSampleMin() const {return fgkSampleMin;};

  // Event processing methods:  
  Bool_t ProcessEvent(AliRawReader *rawReader);
  Bool_t ProcessEvent(AliCaloRawStream    *in);
  
  // Mapping handling
  AliCaloAltroMapping **GetAltroMapping() { return fMapping; };
  void  SetAltroMapping(AliCaloAltroMapping **mapp) { fMapping = mapp; };

  ////////////////////////////
  //Simple getters
  // Main profiles:
  TProfile2D * GetPedProfileLowGain(int i) const {return (TProfile2D*)fPedestalLowGain[i];};	// Return a pointer to the low-gain pedestal profile
  TProfile2D * GetPedProfileHighGain(int i) const {return (TProfile2D*)fPedestalHighGain[i];};	// Return a pointer to the high-gain pedestal profile
  TProfile2D * GetPedRMSProfileLowGain(int i) const {return (TProfile2D*)fPedestalRMSLowGain[i];};	// Return a pointer to the low-gain rms profile 
  TProfile2D * GetPedRMSProfileHighGain(int i) const {return (TProfile2D*)fPedestalRMSHighGain[i];};	// Return a pointer to the high-gain rms profile 
  TProfile2D * GetPeakProfileLowGain(int i) const {return (TProfile2D*)fPeakMinusPedLowGain[i];};	// Return a pointer to the low-gain peak-pedestal profile
  TProfile2D * GetPeakProfileHighGain(int i) const {return (TProfile2D*)fPeakMinusPedHighGain[i];};	// Return a pointer to the high-gain peak-pedestal profile
  
  // Differences to references:
  TProfile2D * GetPedProfileLowGainDiff(int i){ValidateComparisonProfiles(); return (TProfile2D*)fPedestalLowGainDiff[i];};	// Return a pointer to the low-gain pedestal profile difference
  TProfile2D * GetPedProfileHighGainDiff(int i){ValidateComparisonProfiles(); return (TProfile2D*)fPedestalHighGainDiff[i];};	// Return a pointer to the high-gain pedestal profile difference
  TProfile2D * GetPeakProfileLowGainDiff(int i){ValidateComparisonProfiles(); return (TProfile2D*)fPeakMinusPedLowGainDiff[i];};	// Return a pointer to the low-gain peak-pedestal profile difference
  TProfile2D * GetPeakProfileHighGainDiff(int i){ValidateComparisonProfiles(); return (TProfile2D*)fPeakMinusPedHighGainDiff[i];};	// Return a pointer to the high-gain peak-pedestal profile difference
  
  // Ratio to references:
  TProfile2D * GetPedProfileLowGainRatio(int i){ValidateComparisonProfiles(); return (TProfile2D*)fPedestalLowGainRatio[i];};	// Return a pointer to the low-gain pedestal profile ratio
  TProfile2D * GetPedProfileHighGainRatio(int i){ValidateComparisonProfiles(); return (TProfile2D*)fPedestalHighGainRatio[i];};	// Return a pointer to the high-gain pedestal profile ratio
  TProfile2D * GetPeakProfileLowGainRatio(int i){ValidateComparisonProfiles(); return (TProfile2D*)fPeakMinusPedLowGainRatio[i];};	// Return a pointer to the low-gain peak-pedestal profile ratio
  TProfile2D * GetPeakProfileHighGainRatio(int i){ValidateComparisonProfiles(); return (TProfile2D*)fPeakMinusPedHighGainRatio[i];};	// Return a pointer to the high-gain peak-pedestal profile ratio
  
  TH2D * GetDeadMap(int i) const {return (TH2D*)fDeadMap[i];};

  // Basic info: getters  
  kDetType GetDetectorType() const {return fDetType;};//Returns if this is a PHOS or EMCAL object
  TString GetCaloString() const {return fCaloString;}; //Returns if this is a PHOS or EMCAL object
  
  int GetColumns() const {return fColumns;}; //The number of columns per module
  int GetRows() const {return fRows;}; //The number of rows per module
  int GetModules() const {return fModules;}; //The number of modules
  int GetRowMin() const {return fRowMin;}; //for histo def.
  int GetRowMax() const {return fRowMax;}; //for histo def.
  int GetRowMultiplier() const {return fRowMultiplier;}; //for histo filling

  // RunNumbers : setters and getters
  void SetRunNumber(int runNo) {fRunNumber = runNo;};
  int GetRunNumber() const {return fRunNumber;};
  int GetRefRunNumber() const {if (fReference) return fReference->GetRunNumber(); else return -1;};

  // Possibility to select only some samples for the pedestal calculation
  void SetSelectPedestalSamples(Bool_t flag = kFALSE) {fSelectPedestalSamples = flag;} // select to to use only some range of samples for pedestal calc.
  Bool_t GetSelectPedestalSamples() const {return fSelectPedestalSamples;} // select to to use only some range of samples for pedestal calc.
  void SetFirstPedestalSample(int i) {fFirstPedestalSample = i;} // first sample to use
  void SetLastPedestalSample(int i) {fLastPedestalSample = i;} // last sample to use
  int GetFirstPedestalSample() const {return fFirstPedestalSample;}; // first sample to use
  int GetLastPedestalSample() const {return fLastPedestalSample;}; // last sample to use

  // Basic counters
  int GetNEvents() const {return fNEvents;};
  int GetNChanFills() const {return fNChanFills;};
  
  /////////////////////////////
  //Analysis functions
  int GetDeadTowerCount() const {return fDeadTowers;};//Returns the number of dead towers, by counting the bins in peak-pedestal smaller than threshold
  double GetDeadTowerRatio() const {return fDeadTowers/(double)(fRows*fColumns);}; //returns the percentage of dead towers, relative to a full module
  int GetDeadTowerNew() const {return fNewDeadTowers;}; //return the new dead towers compared to the reference
  int GetDeadTowerResurrected() const {return fResurrectedTowers;}; //The the towers resurrected since the reference run

  void Reset();//Resets the whole class.
  Bool_t AddInfo(const AliCaloCalibPedestal *ped);//picks up new info from supplied argument
  
  //////////////////////////////////////////////////////
  //Functions related to comparing this with another (reference) run.
  Bool_t LoadReferenceCalib(TString fileName, TString objectName); //Loads another AliCaloCalibPedestal by name "objectName" from the file "fileName", for reference
  void ComputeDiffAndRatio();//Actually computes the difference and ratio into the histo's in memory
  AliCaloCalibPedestal * GetReference() const {return fReference;}; //Get the reference object. Needed for debug, will probably be removed later
  void ComputeDeadTowers(int threshold = 5, const char * deadMapFile = 0);//Computes the dead tower values
  
  //Saving functions
  Bool_t SaveHistograms(TString fileName, Bool_t saveEmptyHistos = kFALSE); //Saves the histograms to a .root file
  
 private:
  
  void ValidateComparisonProfiles(); //Makes sure that fPe..Diff and fPe..Ratio profiles exist
  
  //The histograms. We use a TObjArray instead of a simple array,because this gives automatic streaming properties for the
  //class. A TClonesArray would be more efficient, but it's a bit more difficult to use and it doesn't matter too much
  //since we have only around 12 objects (maximum) in the array anyway.
  TObjArray fPedestalLowGain; // pedestal info for low gain
  TObjArray fPedestalHighGain; // pedestal info for high gain
  TObjArray fPedestalRMSLowGain; // pedestal rms info for low gain
  TObjArray fPedestalRMSHighGain; // pedestal rms info for high gain
  TObjArray fPeakMinusPedLowGain; // (peak-pedestal) info for low gain
  TObjArray fPeakMinusPedHighGain; // (peak-pedestal) info for high gain
  
  //The difference of profiles between this and the reference object
  TObjArray fPedestalLowGainDiff; //!
  TObjArray fPedestalHighGainDiff; //!
  TObjArray fPeakMinusPedLowGainDiff; //!
  TObjArray fPeakMinusPedHighGainDiff; //!
  
  //The ratio of profiles between this and the reference object
  TObjArray fPedestalLowGainRatio; //!
  TObjArray fPedestalHighGainRatio; //!
  TObjArray fPeakMinusPedLowGainRatio; //!
  TObjArray fPeakMinusPedHighGainRatio; //!
  
  TObjArray fDeadMap;//The deadmap

  // status counters
  int fNEvents; //# total events processed, 
  int fNChanFills; //# total channel fills (NChan * NEvents if not zero-suppressed)

  //The dead tower counts
  int fDeadTowers; //!
  int fNewDeadTowers; //! Towers that have died since the reference run
  int fResurrectedTowers; //! Towers that have been resurrected from the dead, compared to the reference
  
  AliCaloCalibPedestal * fReference; //! A reference object, for comparing the accumulated results to a previous run
  
  kDetType fDetType; //The detector type for this object
  int fColumns;	//The number of columns per module
  int fRows;	//The number of rows per module
  int fModules;	//The number of modules
  int fRowMin; // Mimimum Row number
  int fRowMax; // Maximum now number
  int fRowMultiplier; // Multiplication factor to get proper row range between PHOS and EMCAL
  TString fCaloString; // id for which detector type we have 
  AliCaloAltroMapping **fMapping;    //! Altro Mapping object
  int fRunNumber; //The run number. Needs to be set by the user.
  Bool_t fSelectPedestalSamples; // select to to use only some range of samples for pedestal calc.
  int fFirstPedestalSample; // first sample to use
  int fLastPedestalSample; // last sample to use

  //Constants needed by the class
  static const int fgkSampleMax = 1023; // highest possible sample value (10-bit = 0x3ff)
  static const int fgkSampleMin = 0; // lowest possible sample value 
  
  static const int fgkPhosRows = 64; // number of rows per module for PHOS
  static const int fgkPhosCols = 56; // number of columns per module for PHOS
  static const int fgkPhosModules = 5; // number of modules for PHOS
  
  static const int fgkEmCalRows = 24; // number of rows per module for EMCAL
  static const int fgkEmCalCols = 48; // number of columns per module for EMCAL
  static const int fgkEmCalModules = 12; // number of modules for EMCAL
  
  ClassDef(AliCaloCalibPedestal,3)

};
    
#endif
