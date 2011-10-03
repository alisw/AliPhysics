#ifndef ALIEMCALQADataMakerRec_H
#define ALIEMCALQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.

  Based on PHOS code written by
  Y. Schutz CERN July 2007
 
 Created one histogram for QA shifter;-- Yaxian Mao: 11/2009
 The idea:average counts for all the towers should be flat 
 Change all existing histograms as experts
 
 Change histograms for DQM shifter: --  Yaxian Mao 04/2010
 Calcuate the amplitude ratio from current run and the LED reference, for QAChecker use
 Also calculate the ratio of amplitude from LED Monitor system (current/Reference), to check LED system  
 
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ;
class TH2F ;
class TH2 ; 
class TLine ;
class TText ;
class TProfile ;
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"
class AliCaloRawAnalyzer;
class AliCaloRawAnalyzerKStandard;
class AliEMCALGeoParams;

class AliEMCALQADataMakerRec: public AliQADataMakerRec {

public:
  //Histograms for Raw data control
  enum HRawType_t { 
    // first normal Low Gain and High Gain info
    kNsmodLG,kNsmodHG,kTimeLG,kTimeHG,
    kNtotLG,kNtotHG,kSigHG,kSigLG,
    kPedLG,kPedHG,
    k2DRatioAmp,kRatioDist, kLEDMonRatio, kLEDMonRatioDist,
    // then TRU info
    kNsmodTRU,
    kSigTRU,kNtotTRU,
    kNL0TRU, kTimeL0TRU,
		kNL0FirstTRU, kTimeL0FirstTRU,
    // and also LED Mon info
    kNsmodLGLEDMon,kNsmodHGLEDMon,kTimeLGLEDMon,kTimeHGLEDMon,
    kSigLGLEDMon,kSigHGLEDMon,kNtotLGLEDMon,kNtotHGLEDMon,
    kPedLGLEDMon,kPedHGLEDMon
  } ;

  //Histograms for RecPoints  control
  enum HRPType_t {kRecPE,kRecPM,kRecPDigM};

  //Histograms for ESDs  control
  enum HESDType_t {kESDCaloClusE,kESDCaloClusM,kESDCaloCellA,kESDCaloCellM} ;
                 

public:
  enum fitAlgorithm {kFastFit=1, kNeuralNet = 2, kLMS = 4, kPeakFinder = 5, kCrude = 6};
  AliEMCALQADataMakerRec(fitAlgorithm fitAlgo = kNeuralNet) ;          // ctor
 
  AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) ;   
  AliEMCALQADataMakerRec& operator = (const AliEMCALQADataMakerRec& qadm) ;
  virtual ~AliEMCALQADataMakerRec() {;} // dtor

  Int_t GetFittingAlgorithm() const {return fFittingAlgorithm; }
  void SetFittingAlgorithm(Int_t val);
  AliCaloRawAnalyzer *GetRawAnalyzer() const { return fRawAnalyzer;}
  AliCaloRawAnalyzerKStandard *GetRawAnalyzerTRU() const { return fRawAnalyzerTRU;}

  void SetSuperModules(int i) {fSuperModules = i;}; //The number of SuperModules
  int GetSuperModules() const {return fSuperModules;}; //The number of SuperModules

  // for pedestal calculation with raw data
  void SetFirstPedestalSample(int i) {fFirstPedestalSample = i;}; // first sample 
  int GetFirstPedestalSample() const {return fFirstPedestalSample;}; // first sample 
  void SetLastPedestalSample(int i) {fLastPedestalSample = i;}; // last sample 
  int GetLastPedestalSample() const {return fLastPedestalSample;}; // last sample 
  void SetFirstPedestalSampleTRU(int i) {fFirstPedestalSampleTRU = i;}; // first sample, TRU 
  int GetFirstPedestalSampleTRU() const {return fFirstPedestalSampleTRU;}; // first sample, TRU 
  void SetLastPedestalSampleTRU(int i) {fLastPedestalSampleTRU = i;}; // last sample, TRU 
  int GetLastPedestalSampleTRU() const {return fLastPedestalSampleTRU;}; // last sample, TRU 
  
  // for selection of interesting signal (max-min) range 
  // Low Gain channels
  void SetMinSignalLG(int i) {fMinSignalLG = i;}; 
  int GetMinSignalLG() const {return fMinSignalLG;}; 
  void SetMaxSignalLG(int i) {fMaxSignalLG = i;}; 
  int GetMaxSignalLG() const {return fMaxSignalLG;}; 
  // High Gain channels
  void SetMinSignalHG(int i) {fMinSignalHG = i;}; 
  int GetMinSignalHG() const {return fMinSignalHG;}; 
  void SetMaxSignalHG(int i) {fMaxSignalHG = i;}; 
  int GetMaxSignalHG() const {return fMaxSignalHG;}; 
  // TRU channels
  void SetMinSignalTRU(int i) {fMinSignalTRU = i;}; 
  int GetMinSignalTRU() const {return fMinSignalTRU;}; 
  void SetMaxSignalTRU(int i) {fMaxSignalTRU = i;}; 
  int GetMaxSignalTRU() const {return fMaxSignalTRU;}; 
  // LEDMon channels
  void SetMinSignalLGLEDMon(int i) {fMinSignalLGLEDMon = i;}; 
  int GetMinSignalLGLEDMon() const {return fMinSignalLGLEDMon;}; 
  void SetMaxSignalLGLEDMon(int i) {fMaxSignalLGLEDMon = i;}; 
  int GetMaxSignalLGLEDMon() const {return fMaxSignalLGLEDMon;}; 
  void SetMinSignalHGLEDMon(int i) {fMinSignalHGLEDMon = i;}; 
  int GetMinSignalHGLEDMon() const {return fMinSignalHGLEDMon;}; 
  void SetMaxSignalHGLEDMon(int i) {fMaxSignalHGLEDMon = i;}; 
  int GetMaxSignalHGLEDMon() const {return fMaxSignalHGLEDMon;}; 

  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list) ;
  void           GetCalibRefFromOCDB() ;
  void		 GetTruChannelPosition( Int_t &globRow, Int_t &globColumn, Int_t module, Int_t ddl, Int_t branch, Int_t column ) ;
  virtual void   InitESDs() ; 
  virtual void   InitDigits() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeDigits() ;
  virtual void   MakeDigits(TTree * digTree) ; 
  virtual void   MakeRecPoints(TTree * recpoTree) ; 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   StartOfDetectorCycle() ; 

private:
 void ConvertProfile2H(TProfile * p, TH2 * histo) ; //change the profile plot to a 2D histogram
  
 Int_t fFittingAlgorithm;             // select the fitting algorithm

 AliCaloRawAnalyzer *fRawAnalyzer;    // for signal fitting
 AliCaloRawAnalyzerKStandard *fRawAnalyzerTRU;    // for signal fitting, for TRU

  int fSuperModules; //The number of SuperModules activated
  int fFirstPedestalSample; // first sample for pedestal calculation, in bunch
  int fLastPedestalSample; // last sample for pedestal calculation, in bunch
  int fFirstPedestalSampleTRU; // first sample for pedestal calculation, in bunch
  int fLastPedestalSampleTRU; // last sample for pedestal calculation, in bunch
  int fMinSignalLG; // minimum signal, for Low Gain channels
  int fMaxSignalLG; // maximum signal, for Low Gain channels
  int fMinSignalHG; // minimum signal, for High Gain channels
  int fMaxSignalHG; // maximum signal, for High Gain channels
  int fMinSignalTRU; // minimum signal, for TRU channels
  int fMaxSignalTRU; // maximum signal, for TRU channels
  int fMinSignalLGLEDMon; // minimum signal, for LEDMon channels, low gain
  int fMaxSignalLGLEDMon; // maximum signal, for LEDMon channels, low gain
  int fMinSignalHGLEDMon; // minimum signal, for LEDMon channels, high gain
  int fMaxSignalHGLEDMon; // maximum signal, for LEDMon channels, high gain

  TProfile * fCalibRefHistoPro ; // Profile reference histogram from LED run
  TH2F     * fCalibRefHistoH2F ; // H2F reference histogram from LED run
  TProfile * fLEDMonRefHistoPro ; // Profile reference histogram from LED monitor
  TH2F     * fHighEmcHistoH2F ; // H2F reference histogram from LED run
//  TText **    fTextSM        ; //! Text info for each SM  
//  TLine *     fLineCol       ; //! line to distinguish the different SM side: A side and C side
//  TLine *     fLineRow       ; //! line to distinguish the different SM sector 0 and 1 

  ClassDef(AliEMCALQADataMakerRec,5)  // description 

};

#endif // AliEMCALQADataMakerRec_H
