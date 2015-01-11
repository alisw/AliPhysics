#ifndef ALICALOCALIBSIGNAL_H
#define ALICALOCALIBSIGNAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliCaloCalibSignal.h  $ */

// \file AliCaloCalibSignal.h
//   \brief Description:
//   A help class for monitoring and calibration tools: MOOD, AMORE etc.,
//   that can process events from a standard AliCaloRawStreamV3,
//   most usually from LED/pulser runs. It stores signal info as
//   typical (highest) amplitude vs time in TGraphs (one per channel)
//   or TProfiles if we decide to just store the averages (and not all points) 
//   for the detectors (EMCAL and PHOS).

//   \author: Josh Hamblen (UTenn), original version. 
//   [Consultant: D. Silvermyr (ORNL)]
//   Partly based on AliCaloCalibPedestal.
//   
//   \version $Revision:  $
//   \date $Date: $

#include "TString.h"
#include "TTree.h"
#include "AliEMCALGeoParams.h"
class AliCaloRawStreamV3;
class AliCaloAltroMapping;
class AliRawReader;
class AliCaloRawAnalyzer;

class AliCaloCalibSignal : public TObject {
  
 public:

  enum kDetType {kPhos, kEmCal, kNone};//The detector types
  
  AliCaloCalibSignal(kDetType detectorType = kPhos); //ctor
  virtual ~AliCaloCalibSignal(); //dtor

private:
  //Just declare them, avoid compilation warning
  AliCaloCalibSignal(const AliCaloCalibSignal & /*sig*/); // copy ctor
  AliCaloCalibSignal& operator = (const  AliCaloCalibSignal &/*source*/); // assing operator
  
public:
  // Event processing methods:
  Bool_t ProcessEvent(AliRawReader *rawReader);
  Bool_t ProcessEvent(AliCaloRawStreamV3 *in, UInt_t Timestamp); // added header for time info
  Bool_t CheckFractionAboveAmp(const int *AmpVal, int resultArray[]) const; // check fraction of signals to check for LED events
  Bool_t CheckLEDRefAboveAmp(const int *AmpVal, int resultArray[]) const; // check if LED Ref is also above cut

  // Mapping handling
  AliCaloAltroMapping **GetAltroMapping() const { return fMapping; };
  void  SetAltroMapping(AliCaloAltroMapping **mapp) { fMapping = mapp; };

  // Fitter / Analyzer
  Int_t    GetFittingAlgorithm()  const {return fFittingAlgorithm; }
  void SetFittingAlgorithm(Int_t val) ;         
  AliCaloRawAnalyzer *GetRawAnalyzer()  const { return fRawAnalyzer;}  

  // Parameter/cut handling
  void SetParametersFromFile(const char *parameterFile);
  void WriteParametersToFile(const char *parameterFile);

  ////////////////////////////
  //Simple getters
  // for TTree
  TTree * GetTreeAmpVsTime() const { return fTreeAmpVsTime; } //!
  TTree * GetTreeAvgAmpVsTime() const {return fTreeAvgAmpVsTime; } //!
  TTree * GetTreeLEDAmpVsTime() const {return fTreeLEDAmpVsTime; } //!
  TTree * GetTreeLEDAvgAmpVsTime() const {return fTreeLEDAvgAmpVsTime; } //!

  // how many points do we have for each tower&gain
  int GetNHighGain(int imod, int icol, int irow) const //!
    { int towId = GetTowerNum(imod, icol, irow); return fNHighGain[towId];};	//!
  int GetNLowGain(int imod, int icol, int irow) const //!
    { int towId = GetTowerNum(imod, icol, irow); return fNLowGain[towId];};	//!
  int GetNHighGain(int towId) const { return fNHighGain[towId];};	//!
  int GetNLowGain(int towId) const { return fNLowGain[towId];};	//!

  // also for LED reference
  int GetNRef(const int imod, const int istripMod, const int igain) const //!
    { int refId = GetRefNum(imod, istripMod, igain); return fNRef[refId];}; //!
  int GetNRef(int refId) const { return fNRef[refId];}; //!

  // Basic info: getters  
  kDetType GetDetectorType() const {return fDetType;};//Returns if this is a PHOS or EMCAL object
  TString GetCaloString() const {return fCaloString;}; //Returns if this is a PHOS or EMCAL object  

  int GetColumns() const {return fColumns;}; //The number of columns per module
  int GetRows() const {return fRows;}; //The number of rows per module
  int GetLEDRefs() const {return fLEDRefs;}; //The number of LED references/monitors per module
  int GetModules() const {return fModules;}; //The number of modules

  int GetTowerNum(const int imod, const int icol, const int irow) const { return (imod*fColumns*fRows + icol*fRows + irow);}; // help index

  int GetChannelNum(const int imod, const int icol, const int irow, const int igain) const { return (igain*fModules*fColumns*fRows + imod*fColumns*fRows + icol*fRows + irow);}; // channel number with gain included

  Bool_t DecodeChannelNum(const int chanId, 
			  int *imod, int *icol, int *irow, int *igain) const; // return the module, column, row, and gain for a given channel number

  // LED reference indexing
  int GetRefNum(const int imod, const int istripMod, const int igain) const { return (igain*fModules*fLEDRefs + imod*fLEDRefs + istripMod);}; // channel number with gain included

  Bool_t DecodeRefNum(const int refId, 
		      int *imod, int *istripMod, int *igain) const; // return the module, stripModule, and gain for a given reference number

  // Basic Counters
  int GetNEvents() const {return fNEvents;};
  int GetNAcceptedEvents() const {return fNAcceptedEvents;};

  ///////////////////////////////
  //  Get and Set Cuts
  // Section for if we should help with the event selection of what is likely LED events
  void SetAmpCut(double d) { fAmpCut = d; } //!
  double GetAmpCut() const { return fAmpCut; }; //!
  void SetReqFractionAboveAmpCutVal(double d) { fReqFractionAboveAmpCutVal = d; } //!
  double GetReqFractionAboveAmpCutVal() const { return fReqFractionAboveAmpCutVal; }; //!
  void SetReqFractionAboveAmp(bool b) { fReqFractionAboveAmp = b; } //!
  bool GetReqFractionAboveAmp() const { return fReqFractionAboveAmp; }; //!
  // also for LED Reference/Mon channels
  void SetAmpCutLEDRef(double d) { fAmpCutLEDRef = d; } //!
  double GetAmpCutLEDRef() const { return fAmpCutLEDRef; }; //!
  void SetReqLEDRefAboveAmpCutVal(bool b) { fReqLEDRefAboveAmpCutVal = b; } //!
  bool GetReqLEDRefAboveAmpCutVal() const { return fReqLEDRefAboveAmpCutVal; }; //!

  // We may select to get averaged info
  void SetUseAverage(bool b) { fUseAverage = b; } //!
  bool GetUseAverage() const { return fUseAverage; }; //!
  void SetSecInAverage(int secInAverage) {fSecInAverage = secInAverage;}; // length of the interval that should be used for the average calculation (determines number of bins in TProfile)
  int GetSecInAverage() const {return fSecInAverage;}; //!

  void SetDownscale(int i) {fDownscale = i;}; //!
  int GetDownscale() const {return fDownscale;}; //!

  // Info on time since start of run
  double GetHour() const { return fHour; }; // time info for current event
  double GetCurrentHour() const { return fHour; }; // time info for current event (same as GetHour(), just more explicitly named)
  double GetLatestHour() const { return fLatestHour; }; // the latest time encountered
  // These times are typically the same, but not necessarily if the events do not come in order 
  void SetLatestHour(double d) { fLatestHour = d; }; // could be useful when we know the length of the run (i.e. after it is over), e.g. for PreProcessor

  // RunNumbers : setters and getters
  void SetRunNumber(int runNo) {fRunNumber = runNo;}; //!
  int GetRunNumber() const {return fRunNumber;};  //!
  
  // Start-of-run timestamp : set and get
  void SetStartTime(int startTime) {fStartTime = startTime;}; //!
  int GetStartTime() const {return fStartTime;}; //!

  /////////////////////////////
  //Analysis functions
  void ResetInfo();// trees and counters.
  Bool_t AddInfo(const AliCaloCalibSignal *sig);//picks up new info from supplied argument  

  //Saving functions
  Bool_t Save(TString fileName); //Saves the objects to a .root file
  Bool_t Analyze(); // makes average tree and summary tree 

 private:
 
  void DeleteTrees(); // delete old objects and set pointers
  void Zero(); // set all counters to 0
  void CreateTrees(); //! create/setup the TTrees
    
 private:

  kDetType fDetType; //The detector type for this object
  int fColumns;	//The number of columns per module
  int fRows;	//The number of rows per module
  int fLEDRefs;	//The number of LED references/monitors per module
  int fModules;	//The number of modules
  TString fCaloString; // id for which detector type we have 
  AliCaloAltroMapping **fMapping;    //! Altro Mapping object
  Int_t   fFittingAlgorithm;            // select the fitting algorithm
  AliCaloRawAnalyzer *fRawAnalyzer;     //! e.g. for sample selection for fits
  int fRunNumber; //The run number. Needs to be set by the user.
  int fStartTime;  // Time of first event

  double fAmpCut; // amplitude cut value
  double fReqFractionAboveAmpCutVal; // required fraction that should be above cut
  bool fReqFractionAboveAmp; // flag to select if we should do some event selection based on amplitudes

  double fAmpCutLEDRef; // amplitude cut value for LED reference
  bool fReqLEDRefAboveAmpCutVal; // flag to select if we should require that signal is also seen in LED Reference/Monitoring channel

  double fHour; // fraction of hour since beginning of run, for amp vs. time graphs, for current event
  double fLatestHour; // largest fraction of hour since beginning of run, for amp vs. time graphs
  bool fUseAverage; // flag to average graph points into over a time interval
  int fSecInAverage; // time interval for the graph averaging

  int fDownscale; // to select 1 out every N (fDownscale) events

  // status counters
  int fNEvents; // # events processed
  int fNAcceptedEvents; // # events accepted

  //Constants needed by the class: EMCAL ones are kept in AliEMCALGeoParams.h
  static const int fgkPhosRows = 64; // number of rows per module for PHOS
  static const int fgkPhosCols = 56; // number of columns per module for PHOS
  static const int fgkPhosLEDRefs = 0; // no LED monitor channels for PHOS
  static const int fgkPhosModules = 5; // number of modules for PHOS
  
  // From numbers above: PHOS has more possible towers (17920) than EMCAL (13824) 
  // so use PHOS numbers to set max. array sizes
  static const int fgkMaxTowers = 17920; // fgkPhosModules * fgkPhosCols * fgkPhosRows; 
  // for LED references; maximum from EMCAL
  static const int fgkMaxRefs = 288; // AliEMCALGeoParams::fgkEMCALModules * AliEMCALGeoParams::fgkEMCALLEDRefs

  static const int fgkNumSecInHr = 3600;  // number of seconds in an hour, for the fractional hour conversion on the time graph
  
  // trees
  TTree *fTreeAmpVsTime; // stores channel, gain, amp, and time info
  TTree *fTreeAvgAmpVsTime; // same, for averages
  TTree *fTreeLEDAmpVsTime; // same, for LED reference
  TTree *fTreeLEDAvgAmpVsTime; // same, for LED reference - averages

  // counters
  int fNHighGain[fgkMaxTowers]; // Number of Amp. vs. Time readings per tower
  int fNLowGain[fgkMaxTowers]; // same, for low gain
  int fNRef[fgkMaxRefs * 2]; // same, for LED refs; *2 for both gains
  
  ClassDef(AliCaloCalibSignal, 8) // don't forget to change version if you change class member list..
    
};
    
#endif
