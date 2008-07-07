#ifndef ALITOFCALIB_H
#define ALITOFCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//  class for TOF calibration:: simulation of uncalibrated data //
//////////////////////////////////////////////////////////////////

#define CHENTRIESSMALL       300   // number of entries per TOF channel per run
                                   // (to be divided by 3 to get the 
                                   // real number of entries), smallarray 
#define MAXCHENTRIESSMALL  10000   // max number of entries per TOF channel 
                                   // (to be divided by 3 to get the 
                                   // real number of entries), smallarray 
				   
#define NIDXSMALL              3   // number of values stored 
                                   // after Comb PID per ESD track
#define DELTAIDXTOT            0   // index for ToT in bigarray
#define DELTAIDXTIME           1   // index for TOF time in big/smallarray
#define DELTAIDXPID            2   // index for Exp Time after Comb PID 
                                   // in smallarray
#define MEANENTRIES           15   // Mean number of entries per channel 
                                   // to perform calibration

#include "TTask.h"

class TArrayF;
class TF1;
class TH1F;
class TObjArray;
class TTree;
class TChain;
class TMap;

//class AliESD;

class AliTOFCal;
class AliTOFRecoParam;
class AliTOFChannelOnlineStatusArray;
class AliTOFChannelOnlineArray;

class AliTOFcalib:public TTask{
public:
  AliTOFcalib();          // ctor
  AliTOFcalib(const AliTOFcalib & calib); // copy constructor
  AliTOFcalib& operator=(const AliTOFcalib & calib); // assignment operator
  virtual ~AliTOFcalib() ; // dtor
  void CreateCalArrays();
  void CreateCalObjects();
  TObjArray * GetTOFCalArrayOnline() const {return fTOFCalOnline;}
  AliTOFChannelOnlineArray * GetTOFOnlineDelay() const {return fCal;}
  AliTOFChannelOnlineStatusArray * GetTOFOnlineStatus() const {return fStatus;}
  TObjArray * GetTOFCalArrayOnlinePulser() const {return fTOFCalOnlinePulser;}
  TObjArray * GetTOFCalArrayOnlineNoise() const {return fTOFCalOnlineNoise;}
  TObjArray * GetTOFCalArrayOnlineHW() const {return fTOFCalOnlineHW;}
  TObjArray * GetTOFCalArrayOffline() const {return fTOFCalOffline;}
  TMap * GetConfigMap() const {return fConfigMap;}
  TH1F * GetTOFSimToT() const {return fTOFSimToT;}
  TTree * GetTOFCalibTree() const {return fTree;}
  TChain * GetTOFCalibChain() const {return fChain;}
  const char * GetOfflineValidity() const {return fkValidity;}
  void SetOfflineValidity(const char* validity) {fkValidity = validity;}
  Int_t NChannels()const{return fNChannels;}

  // Methods to retrieve/write parameters from/on CDB
  // writing

  void WriteSimHistoOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun, TH1F *histo);
  void WriteConfigMapOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteConfigMapOnCDB(Char_t *sel);
  // new calib objs
  void WriteParOnlineDelayOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlineStatusOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlineDelayOnCDB(Char_t *sel);
  void WriteParOnlineStatusOnCDB(Char_t *sel);
  // old calib objs
  void WriteParOnlineOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlinePulserOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlineNoiseOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlineHWOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOfflineOnCDB(Char_t *sel, const Char_t *validity, Int_t minrun, Int_t maxrun);
  void WriteParOnlineOnCDB(Char_t *sel);
  void WriteParOnlinePulserOnCDB(Char_t *sel);  // old, before unification of status info
  void WriteParOnlineNoiseOnCDB(Char_t *sel);   // old, before unification of status info
  void WriteParOnlineHWOnCDB(Char_t *sel);      // old, before unification of status info
  void WriteParOfflineOnCDB(Char_t *sel, const Char_t *validity);

  // reading
  Bool_t ReadSimHistoFromCDB(Char_t *sel, Int_t nrun);
  Bool_t ReadConfigMapFromCDB(Char_t *sel, Int_t nrun);
  // new objs
  Bool_t ReadParOnlineDelayFromCDB(Char_t *sel, Int_t nrun);
  Bool_t ReadParOnlineStatusFromCDB(Char_t *sel, Int_t nrun);
  // old objs
  Bool_t ReadParOnlineFromCDB(Char_t *sel, Int_t nrun);
  Bool_t ReadParOnlinePulserFromCDB(Char_t *sel, Int_t nrun);  // old, before unification of status info
  Bool_t ReadParOnlineNoiseFromCDB(Char_t *sel, Int_t nrun);   // old, before unification of status info
  Bool_t ReadParOnlineHWFromCDB(Char_t *sel, Int_t nrun);      // old, before unification of status info
  Bool_t ReadParOfflineFromCDB(Char_t *sel, Int_t nrun);
  void WriteRecParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun, AliTOFRecoParam *param);
  AliTOFRecoParam * ReadRecParFromCDB(Char_t *sel, Int_t nrun);
  void CreateTreeFromCDB(Int_t minrun, Int_t maxrun);
  void CreateTreeFromFile(Int_t minrun, Int_t maxrun);
  void CreateTreeFromGrid(Int_t minrun, Int_t maxrun);
  void CreateChainFromGrid(Int_t minrun, Int_t maxrun);
  Int_t Calibrate(Option_t *optionSave="", Option_t *optionFit="RQ");
  Int_t Calibrate(Int_t nch,Int_t *ich, Option_t *optionSave="", Option_t *optionFit="RQ");
  Int_t Calibrate(Int_t ichmin, Int_t ichmax, Option_t *optionSave="", Option_t *optionFit="RQ");
  Int_t Calibrate(Int_t ich, Option_t *optionSave="", Option_t *optionFit="RQ");
  Int_t CalibrateFromProfile(Int_t ich, Option_t *optionSave="", Option_t *optionFit="RQ");
  TH1F* Profile(Int_t i);
  Int_t FindBins (TH1F* h, Double_t *bins) const;
  void SetNruns(Int_t nruns) {fNruns=nruns;}
  Int_t GetNruns() const {return fNruns;}
  void SetFirstRun(Int_t firstRun) {fFirstRun=firstRun;}
  Int_t GetFirstRun() const {return fFirstRun;}
  void SetLastRun(Int_t lastRun) {fLastRun=lastRun;}
  Int_t GetLastRun() const {return fLastRun;}

private:
  Int_t fNChannels; // number of TOF channels

  // old calibration objects
  TObjArray *fTOFCalOnline;       // array of AliTOFChannels storing calib parameters
  TObjArray *fTOFCalOnlinePulser; // array of AliTOFChannels storing calib status from pulser   // old, before unification of status info
  TObjArray *fTOFCalOnlineNoise;  // array of AliTOFChannels storing calib status from noise    // old, before unification of status info
  TObjArray *fTOFCalOnlineHW;  // array of AliTOFChannels storing calib status from hardware    // old, before unification of status info
  TObjArray *fTOFCalOffline;       // array of AliTOFChannels storing calib parameters

  // new calibration objects
  AliTOFChannelOnlineArray *fCal; // object with delay array for TOF channels
  AliTOFChannelOnlineStatusArray *fStatus; // object with status array for TOF channels

  TH1F *fTOFSimToT;        // histo with realistic ToT signal from TB Data
  const char *fkValidity;  // validity for offline calibration object
  TTree *fTree;            // tree for TOF calibration
  TChain *fChain;          // chain for TOF calibration
  Int_t fNruns;            // number of runs to be processed
  Int_t fFirstRun;            // first run for calibration obj validity
  Int_t fLastRun;            // last run for calib obj validity
  TMap* fConfigMap;          // map holding configuration obj

  ClassDef(AliTOFcalib,7);
};

#endif // AliTOFcalib_H

