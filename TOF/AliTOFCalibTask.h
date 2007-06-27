#ifndef ALITOFCALIBTASK_H
#define ALITOFCALIBTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////
//    task for TOF calibration            //
//         C.Zampolli                     //
////////////////////////////////////////////

/* $Id$ */

#define TOFCHANNELS       160000   // number TOF channels 
#define CHENTRIES            500   // number of entries per TOF channel per run
                                   // (to be divided by 5 to get the 
                                   // real number of entries), bigarray 
#define CHENTRIESSMALL       300   // number of entries per TOF channel per run
                                   // (to be divided by 3 to get the 
                                   // real number of entries), smallarray 
#define MAXCHENTRIESSMALL  10000   // max number of entries per TOF channel 
                                   // (to be divided by 3 to get the 
                                   // real number of entries), smallarray 
				   
#define LOWERMOMBOUND        1.0   // [GeV/c] default value Pb-Pb
#define UPPERMOMBOUND        1.8   // [GeV/c] default value Pb-Pb
#define MINTIME              5E3   // min time of flight value (ns)
#define NIDX                   5   // number of values stored 
                                   // in big/smallarray per ESD track
#define NIDXSMALL              3   // number of values stored 
                                   // after Comb PID per ESD track
#define DELTAIDXTOT            0   // index for ToT in bigarray
#define DELTAIDXTIME           1   // index for TOF time in big/smallarray
#define DELTAIDXEXTIMEPI       2   // index for Exp Time Pi in bigarray
#define DELTAIDXEXTIMEKA       3   // index for Exp Time Ka in bigarray
#define DELTAIDXEXTIMEPR       4   // index for Exp Time Pr in bigarray
#define DELTAIDXPID            2   // index for Exp Time after Comb PID 
                                   // in smallarray
#define TRACKERROR       90*1E-3   // error on the tracks for 
                                   // Combinatorial PID (ns)
#define MEANENTRIES           15   // Mean number of entries per channel 
                                   // to perform calibration

#include "AliAnalysisTask.h"  

class TTree;
class AliESDtrack ;  
class TFile; 
class TH1F ;
class TH1I ;
class TH1D ;
class TH2F ;
class AliESD ; 

class AliTOFCalibTask : public AliAnalysisTask {

public:
  AliTOFCalibTask(const char *name) ; //ctor
  AliTOFCalibTask(const AliTOFCalibTask & calibtask); // copy constructor
  AliTOFCalibTask& operator=(const AliTOFCalibTask & calibtask); // assignment operator
  virtual ~AliTOFCalibTask(); //dtor
  virtual void Exec(Option_t * opt="") ;
  virtual void ConnectInputData(Option_t *) ;
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

 public: 
  class AliTOFArrayTask : public TObject {
  public:
    AliTOFArrayTask(): TObject(),fSize(0),fArray(0x0){}
    AliTOFArrayTask(Int_t size) :
      TObject(),
      fSize(size),
      fArray(new TArrayF*[size]) {
    } 
    AliTOFArrayTask(const AliTOFArrayTask & source):
      TObject(),fSize(0),fArray(0x0){ // copy constructor
      this->fSize= source.fSize;
      this->fArray= source.fArray;
    };

    AliTOFArrayTask& operator=(const AliTOFArrayTask & source) { // assignment operator
      this->fSize= source.fSize;
      this->fArray= source.fArray;
      return *this;
    }

    Int_t GetSize() const {return fSize;}
    void AddArray(Int_t pos, TArrayF * parr) {
      if (pos>-1 && pos < fSize)fArray[pos] = parr;}
    TArrayF *  GetArray(Int_t pos) const {
      TArrayF * parr = 0x0; 
      if  (pos>-1 && pos < fSize)parr = fArray[pos];
      return parr;
    }
    virtual ~AliTOFArrayTask() {
      delete [] fArray;
    }
    
  private:
    
    Int_t fSize;       // Size of the array of TArrayFs
    TArrayF ** fArray; //[fSize]};
    
  };

private:
  Bool_t Select(AliESDtrack *t);
  Bool_t CombPID(Float_t *smallarray, Int_t size);
  void BookHistos();
  void DrawHistos();
  Int_t Calibrate(Option_t *optionSave="", Option_t *optionFit="RQ");
  Int_t Calibrate(Int_t nch,Int_t *ich, Option_t *optionSave="", Option_t *optionFit="RQ");
  Int_t Calibrate(Int_t ichmin, Int_t ichmax, Option_t *optionSave="", Option_t *optionFit="RQ");
  Int_t Calibrate(Int_t ich, Option_t *optionSave="", Option_t *optionFit="RQ");
  Int_t CalibrateFromProfile(Int_t ich, Option_t *optionSave="", Option_t *optionFit="RQ");
  TH1F* Profile(Int_t i);
  Int_t FindBins (TH1F* h, Double_t *bins) const;

  const Char_t *fdir;              // initial directory
  TTree* fChain ;            //!pointer to the analyzed TTree or TChain
  //leaf types
  AliESD* fESD ;              //! Declaration of leave types
  Int_t fMinEntries;            // minimum number of entries to steer the task
  Int_t fIch;                   // TOF channel number
  Float_t fToT;                 // Time over Threshold, ns
  Float_t fTime;                // TOF time, ns
  Float_t fExpTimePi;           // exp time, Pions, ns
  Float_t fExpTimeKa;           // exp time, Kaons, ns
  Float_t fExpTimePr;           // exp time, Protons, ns
  TTree* ftree;                 // tree for calibration
  Float_t fMinTime;             // min TOF time for track selection; not used
  Float_t** fbigarray;          // big array for calibration
  Int_t* findexarray;           // array for entry index in each channel 
  Int_t fnESD;                  // number of analyzed ESD tracks 
  Int_t fnESDselected;          // number of selected ESD tracks 
  Int_t fnESDkTOFout;           // number of ESD tracks with kTOFout
  Int_t fnESDkTIME;             // number of ESD tracks with kTIME
  Int_t fnESDassTOFcl;          // number of ESD tracks with assigned TOF cluster
  Int_t fnESDTIMEcut;           // number of ESD tracks with TOF time < 17 ns
  Int_t fnESDTRDcut;            // number of ESD tracks with TRD ok

  // Histos
  TH1F* fhToT;                  // ToT histo
  TH1F* fhTime;                 // Time histo
  TH1F* fhExpTimePi;            // Exp Time Pi histo                  
  TH1F* fhExpTimeKa;            // Exp Time Ka histo                  
  TH1F* fhExpTimePr;            // Exp Time Pr histo                  
  TH1I* fhPID;                  // PID histo                  
  TH1D* fhch;                   // TOF channel histo               

  TObjArray * fOutputContainer ; //! output data container

  ClassDef(AliTOFCalibTask, 1); //  TOF Calib task 
};
#endif // ALITOFCALIBTASK_H

