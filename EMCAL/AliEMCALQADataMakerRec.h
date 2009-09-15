#ifndef ALIEMCALQADataMakerRec_H
#define ALIEMCALQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.

  Based on PHOS code written by
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"

class AliEMCALQADataMakerRec: public AliQADataMakerRec {

public:
  //Histograms for Raw data control
  enum HRawType_t { 
    // first normal Low Gain and High Gain info
    kNsmodLG,kNsmodHG,kTimeLG,kTimeHG,
    kSigLG,kSigHG,kNtotLG,kNtotHG,
    kPedLG,kPedHG,
    kPedRMSLG,kPedRMSHG,
    // then TRU info
    kNsmodTRU,kTimeTRU,
    kSigTRU,kNtotTRU,
    kPedTRU,kPedRMSTRU,
    // and also LED Mon info
    kNsmodLGLEDMon,kNsmodHGLEDMon,kTimeLGLEDMon,kTimeHGLEDMon,
    kSigLGLEDMon,kSigHGLEDMon,kNtotLGLEDMon,kNtotHGLEDMon,
    kPedLGLEDMon,kPedHGLEDMon,
    kPedRMSLGLEDMon,kPedRMSHGLEDMon
  } ;

  //Histograms for RecPoints  control
  enum HRPType_t {kRecPE,kRecPM,kRecPDigM};

  //Histograms for ESDs  control
  enum HESDType_t {kESDCaloClusE,kESDCaloClusM,kESDCaloCellA,kESDCaloCellM} ;
                 

public:
  AliEMCALQADataMakerRec() ;          // ctor
  AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) ;   
  AliEMCALQADataMakerRec& operator = (const AliEMCALQADataMakerRec& qadm) ;
  virtual ~AliEMCALQADataMakerRec() {;} // dtor

  void SetSuperModules(int i) {fSuperModules = i;}; //The number of SuperModules
  int GetSuperModules() const {return fSuperModules;}; //The number of SuperModules

  // for pedestal calculation with raw data
  void SetFirstPedestalSample(int i) {fFirstPedestalSample = i;}; // first sample 
  int GetFirstPedestalSample() const {return fFirstPedestalSample;}; // first sample 
  void SetLastPedestalSample(int i) {fLastPedestalSample = i;}; // last sample 
  int GetLastPedestalSample() const {return fLastPedestalSample;}; // last sample 
  // for selection of interesting signal (max-min) range for High Gain channels
  // (useful for MIP/cosmic studies) 
  void SetMinSignalHG(int i) {fMinSignalHG = i;}; // minimum signal
  int GetMinSignalHG() const {return fMinSignalHG;}; // minimum signal
  void SetMaxSignalHG(int i) {fMaxSignalHG = i;}; // maximum signal
  int GetMaxSignalHG() const {return fMaxSignalHG;}; // maximum signal

  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list) ;
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
  int fSuperModules; //The number of SuperModules activated
  int fFirstPedestalSample; // first sample for pedestal calculation
  int fLastPedestalSample; // last sample for pedestal calculation
  int fMinSignalHG; // minimum signal, for High Gain channels
  int fMaxSignalHG; // maximum signal, for High Gain channels

  ClassDef(AliEMCALQADataMakerRec,4)  // description 

};

#endif // AliEMCALQADataMakerRec_H
