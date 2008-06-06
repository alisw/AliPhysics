#ifndef ALICALOCALIBSIGNAL_H
#define ALICALOCALIBSIGNAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliCaloCalibSignal.h  $ */

// \file AliCaloCalibSignal.h
//   \brief Description:
//   A help class for monitoring and calibration tools: MOOD, AMORE etc.,
//   that can process events from a standard AliCaloRawStream,
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

#include "TGraph.h"
#include "TProfile.h"
class AliCaloRawStream;
class AliCaloAltroMapping;
class AliRawReader;
class AliRawEventHeaderBase;

class AliCaloCalibSignal : public TObject {
  
 public:

  enum kDetType {kPhos, kEmCal, kNone};//The detector types
  
  AliCaloCalibSignal(kDetType detectorType = kPhos); //ctor
  virtual ~AliCaloCalibSignal(); //dtor

  // copy ctor, and '=' operator, are not fully tested/debugged yet
  AliCaloCalibSignal(const AliCaloCalibSignal &sig); // copy ctor
  AliCaloCalibSignal& operator = (const  AliCaloCalibSignal &source); //!
  
  // Event processing methods:
  Bool_t ProcessEvent(AliRawReader *rawReader);
  Bool_t ProcessEvent(AliCaloRawStream *in, AliRawEventHeaderBase *aliHeader); // added header for time info
  Bool_t CheckFractionAboveAmp(int *AmpVal, int nTotChan); // check fraction of signals to check for LED events

  // Mapping handling
  AliCaloAltroMapping **GetAltroMapping() { return fMapping; };
  void  SetAltroMapping(AliCaloAltroMapping **mapp) { fMapping = mapp; };

  ////////////////////////////
  //Simple getters
  // need public access to the TGraphs.. 
  TGraph * GetGraphAmpVsTimeHighGain(int imod, int icol, int irow) const // Return a pointer to the high gain graph
    { int towId = GetTowerNum(imod, icol, irow); return fGraphAmpVsTimeHighGain[towId];}; //!
  TGraph * GetGraphAmpVsTimeLowGain(int imod, int icol, int irow) const // Return a pointer to the low gain graph
    {int towId = GetTowerNum(imod, icol, irow); return fGraphAmpVsTimeLowGain[towId];};	
  TGraph * GetGraphAmpVsTimeHighGain(int towId) const // Return a pointer to the high gain graph
    { return fGraphAmpVsTimeHighGain[towId];}; //!
  TGraph * GetGraphAmpVsTimeLowGain(int towId) const // Return a pointer to the low gain graph
    { return fGraphAmpVsTimeLowGain[towId];}; //!	

  // and similarly for the TProfiles
  TProfile * GetProfAmpVsTimeHighGain(int imod, int icol, int irow) const // Return a pointer to the high gain profile
    { int towId = GetTowerNum(imod, icol, irow); return fProfAmpVsTimeHighGain[towId];}; //!	
  TProfile * GetProfAmpVsTimeLowGain(int imod, int icol, int irow) const // Return a pointer to the low gain profile
    { int towId = GetTowerNum(imod, icol, irow); return fProfAmpVsTimeLowGain[towId];}; //!	
  TProfile * GetProfAmpVsTimeHighGain(int towId) const // Return a pointer to the high gain profile
    { return fProfAmpVsTimeHighGain[towId];}; //!	
  TProfile * GetProfAmpVsTimeLowGain(int towId) const // Return a pointer to the low gain profile
    { return fProfAmpVsTimeLowGain[towId];}; //!	

  // how many points do we have in each TGraph
  int GetNHighGain(int imod, int icol, int irow) const //!
    { int towId = GetTowerNum(imod, icol, irow); return fNHighGain[towId];};	//!
  int GetNLowGain(int imod, int icol, int irow) const //!
    { int towId = GetTowerNum(imod, icol, irow); return fNLowGain[towId];};	//!
  int GetNHighGain(int towId) const { return fNHighGain[towId];};	//!
  int GetNLowGain(int towId) const { return fNLowGain[towId];};	//!

  // Basic info: getters  
  kDetType GetDetectorType() const {return fDetType;};//Returns if this is a PHOS or EMCAL object
  TString GetCaloString() const {return fCaloString;}; //Returns if this is a PHOS or EMCAL object  

  int GetColumns() const {return fColumns;}; //The number of columns per module
  int GetRows() const {return fRows;}; //The number of rows per module
  int GetModules() const {return fModules;}; //The number of modules
  int GetTowerNum(int imod, int icol, int irow) const { return imod*fColumns*fRows + icol*fRows + irow;}; // help index

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
  double GetReqFractionAboveAmp() const { return fReqFractionAboveAmp; }; //!

  // We may select to only use the averaged info in the TProfiles rather than the
  // the full in the TGraphs
  void SetUseAverage(bool b) { fUseAverage = b; } //!
  double GetUseAverage() const { return fUseAverage; }; //!
  void SetSecInAverage(int secInAverage) {fSecInAverage = secInAverage;}; // length of the interval that should be used for the average calculation (determines number of bins in TProfile)
  int GetSecInAverage() const {return fSecInAverage;}; //!

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
  void Reset();//Resets the whole class.
  Bool_t AddInfo(const AliCaloCalibSignal *sig);//picks up new info from supplied argument  

  //Saving functions
  Bool_t Save(TString fileName, Bool_t saveEmptyGraphs = kFALSE); //Saves the TGraphs to a .root file

 private:

  void ClearObjects(); // delete old objects and set pointers
  void Zero(); // set all counters to 0
  void CreateGraphs(); //! create/setup the TGraphs
  void CreateProfile(int imod, int icol, int irow, int towerId, int gain,
		     int nbins, double min, double max); //! create/setup a TProfile
    
 private:

  kDetType fDetType; //The detector type for this object
  int fColumns;	//The number of columns per module
  int fRows;	//The number of rows per module
  int fModules;	//The number of modules
  TString fCaloString; // id for which detector type we have 
  AliCaloAltroMapping **fMapping;    //! Altro Mapping object
  int fRunNumber; //The run number. Needs to be set by the user.
  int fStartTime;  // Time of first event

  double fAmpCut; // amplitude cut value
  double fReqFractionAboveAmpCutVal; // required fraction that should be above cut
  bool fReqFractionAboveAmp; // flag to select if we should do some event selection based on amplitudes

  double fHour; // fraction of hour since beginning of run, for amp vs. time graphs, for current event
  double fLatestHour; // largest fraction of hour since beginning of run, for amp vs. time graphs
  bool fUseAverage; // flag to average graph points into over a time interval
  int fSecInAverage; // time interval for the graph averaging

  // status counters
  int fNEvents; // # events processed
  int fNAcceptedEvents; // # events accepted

  //Constants needed by the class
  static const int fgkSampleMax = 1023; // highest possible sample value (10-bit = 0x3ff)
  static const int fgkSampleMin = 0; // lowest possible sample value 
  
  static const int fgkPhosRows = 64; // number of rows per module for PHOS
  static const int fgkPhosCols = 56; // number of columns per module for PHOS
  static const int fgkPhosModules = 5; // number of modules for PHOS
  
  static const int fgkEmCalRows = 24; // number of rows per module for EMCAL
  static const int fgkEmCalCols = 48; // number of columns per module for EMCAL
  static const int fgkEmCalModules = 12; // number of modules for EMCAL

  // From numbers above: PHOS has more possible towers (17920) than EMCAL (13824) 
  // so use PHOS numbers to set max. array sizes
  static const int fgkMaxTowers = 17920; // fgkPhosModules * fgkPhosCols * fgkPhosRows; 
  
  static const int fgkNumSecInHr = 3600;  // number of seconds in an hour, for the fractional hour conversion on the time graph
  
  TGraph *fGraphAmpVsTimeHighGain[fgkMaxTowers]; // Amplitude vs. Time Graph for each high gain channel
  TGraph *fGraphAmpVsTimeLowGain[fgkMaxTowers]; // Amplitude vs. Time Graph for each low gain channel
  TProfile *fProfAmpVsTimeHighGain[fgkMaxTowers]; // Amplitude vs. Time Profile for each high gain channel
  TProfile *fProfAmpVsTimeLowGain[fgkMaxTowers]; // Amplitude vs. Time Profile for each low gain channel
  
  int fNHighGain[fgkMaxTowers]; // Number of points for each Amp. vs. Time graph
  int fNLowGain[fgkMaxTowers]; // Number of points for each Amp. vs. Time graph
  
  ClassDef(AliCaloCalibSignal,1)
    
};
    
#endif
