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
class TH1C;
class TObjArray;
class TTree;
class TChain;
class TMap;

//class AliESD;

class AliTOFCal;
class AliTOFRecoParam;
class AliTOFChannelOnlineStatusArray;
class AliTOFChannelOnlineArray;
class AliTOFDeltaBCOffset;
class AliTOFCTPLatency;
class AliTOFT0Fill;
class AliTOFRunParams;
class AliTOFResponseParams;
class AliESDEvent;
class AliLHCClockPhase;

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

  void CreateDeltaBCOffset();
  void CreateCTPLatency();
  void CreateT0Fill();
  void CreateRunParams();
  AliTOFDeltaBCOffset *GetDeltaBCOffset() const {return fDeltaBCOffset;};
  AliTOFCTPLatency *GetCTPLatency() const {return fCTPLatency;};
  AliTOFT0Fill *GetT0Fill() const {return fT0Fill;};
  AliTOFRunParams *GetRunParams() const {return fRunParams;};
  AliTOFResponseParams *GetResponseParams() const {return fResponseParams;};

  // Methods to retrieve/write parameters from/on CDB
  // writing

  void WriteSimHistoOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun, TH1F *histo);
  void WriteConfigMapOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteConfigMapOnCDB(const Char_t *sel);
  // new calib objs
  void WriteParOnlineDelayOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlineStatusOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlineDelayOnCDB(const Char_t *sel);
  void WriteParOnlineStatusOnCDB(const Char_t *sel);
  // old calib objs
  void WriteParOnlineOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlinePulserOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlineNoiseOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOnlineHWOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteParOfflineOnCDB(const Char_t *sel, const Char_t *validity, Int_t minrun, Int_t maxrun);
  void WriteParOnlineOnCDB(const Char_t *sel);
  void WriteParOnlinePulserOnCDB(const Char_t *sel);  // old, before unification of status info
  void WriteParOnlineNoiseOnCDB(const Char_t *sel);   // old, before unification of status info
  void WriteParOnlineHWOnCDB(const Char_t *sel);      // old, before unification of status info
  void WriteParOfflineOnCDB(const Char_t *sel, const Char_t *validity);

  void WriteDeltaBCOffsetOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteCTPLatencyOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteT0FillOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteRunParamsOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteReadoutEfficiencyOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteProblematicOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun);

  // reading
  Bool_t ReadSimHistoFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadConfigMapFromCDB(const Char_t *sel, Int_t nrun);
  // new objs
  Bool_t ReadParOnlineDelayFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadParOnlineStatusFromCDB(const Char_t *sel, Int_t nrun);
  // old objs
  Bool_t ReadParOnlineFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadParOnlinePulserFromCDB(const Char_t *sel, Int_t nrun);  // old, before unification of status info
  Bool_t ReadParOnlineNoiseFromCDB(const Char_t *sel, Int_t nrun);   // old, before unification of status info
  Bool_t ReadParOnlineHWFromCDB(const Char_t *sel, Int_t nrun);      // old, before unification of status info
  Bool_t ReadParOfflineFromCDB(const Char_t *sel, Int_t nrun);
  void WriteRecParOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun, AliTOFRecoParam *param);
  void WriteRecParOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun, TObjArray *arr);
  AliTOFRecoParam * ReadRecParFromCDB(const Char_t *sel, Int_t nrun, Int_t eventType=0);
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

  Bool_t ReadDeltaBCOffsetFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadCTPLatencyFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadT0FillFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadRunParamsFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadLHCClockPhaseFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadReadoutEfficiencyFromCDB(const Char_t *sel, Int_t nrun);
  Bool_t ReadProblematicFromCDB(const Char_t *sel, Int_t nrun);

  Bool_t Init(Int_t run = -1); // init
  Double_t GetTimeCorrection(Int_t index, Double_t tot, Int_t deltaBC, Int_t l0l1, UInt_t timestamp); // get time correction
  void CalibrateESD(AliESDEvent *event); // calibrate ESD
  void CalibrateTExp(AliESDEvent *event) const; // calibrate TExp
  void SetRemoveMeanT0(Bool_t value) {fRemoveMeanT0 = value;}; // setter
  void SetUseLHCClockPhase(Bool_t value) {fUseLHCClockPhase = value;}; // setter
  Bool_t GetUseLHCClockPhase() const {return fUseLHCClockPhase;}; // getter
  void SetCalibrateTOFsignal(Bool_t value) {fCalibrateTOFsignal = value;}; // setter
  void SetCorrectTExp(Bool_t value) {fCorrectTExp = value;}; // setter
  Bool_t IsChannelEnabled(Int_t index, Bool_t checkEfficiency = kTRUE, Bool_t checkProblematic = kTRUE); // is channel enabled
  Bool_t IsChannelEfficient(Int_t index); // is channel efficient
  Bool_t IsChannelProblematic(Int_t index); // is channel problematic
  Double_t TuneForMC(AliESDEvent *event, Double_t resolution); // tune for MC

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

  AliTOFDeltaBCOffset *fDeltaBCOffset; // deltaBC offset
  AliTOFCTPLatency *fCTPLatency; // CTP latency
  AliTOFT0Fill *fT0Fill; // T0 fill
  AliTOFRunParams *fRunParams; // run params
  AliLHCClockPhase *fLHCClockPhase; // LHC clock-phase
  AliTOFResponseParams *fResponseParams; // run params
  TH1F *fReadoutEfficiency; // readout efficiency
  TH1C *fProblematic; // problematic
  
  Bool_t fInitFlag; // init flag
  Bool_t fRemoveMeanT0; // remove mean T0
  Bool_t fUseLHCClockPhase; // use LHC clock-phase
  Bool_t fCalibrateTOFsignal; // calibrate TOF signal
  Bool_t fCorrectTExp; // correct expected time

  ClassDef(AliTOFcalib,11);
};

#endif // AliTOFcalib_H

