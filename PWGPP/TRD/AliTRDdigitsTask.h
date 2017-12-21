#ifndef ALITRDDIGITSTASK_H
#define ALITRDDIGITSTASK_H

// example of an analysis task to analyse TRD digits 
// Authors: Tom Dietel
// based on the AliAnalysisTaskPt

class AliTRDdigitsManager;

#include "AliAnalysisTaskSE.h"

class AliTRDdigitsTask : public AliAnalysisTaskSE {
public:
  AliTRDdigitsTask()
    : AliAnalysisTaskSE(),
      fDigitsInputFileName("TRD.Digits.root"), fDigitsInputFile(0),
      fDigitsOutputFileName(""), fDigitsOutputFile(0),
      fDigMan(0),fGeo(0), fEventNoInFile(-2), fDigitsLoadedFlag(kFALSE)
  {}
  AliTRDdigitsTask(const char *name);
  //virtual ~AliTRDdigitsTask() {}
  
  //virtual void   UserCreateOutputObjects();
  virtual Bool_t UserNotify();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetDigitsInputFilename(TString x) {fDigitsInputFileName=x;}
  void SetDigitsOutputFilename(TString x) {fDigitsOutputFileName=x;}
  
protected:

  Bool_t NextEvent(Bool_t preload=kFALSE);
  Bool_t ReadDigits();
  Bool_t WriteDigits();
  
  AliTRDtrackV1* FindTRDtrackV1(AliESDfriendTrack* friendtrack);

  Int_t FindDigitsTrkl(const AliTRDtrackV1* trdTrack, Int_t layer,
		       Int_t* det, Int_t* row, Int_t* col,
		       Float_t* x, Float_t* y, Float_t* z);
  
  Int_t FindDigits(const AliExternalTrackParam* param,
		   Float_t bfield, Int_t layer,
		   Int_t* det, Int_t* row, Int_t* col);

private:

  TFile* OpenDigitsFile(TString inputfile, TString digfile, TString opt);

  TString fDigitsInputFileName;         //! Name of digits file for reading
  TFile*  fDigitsInputFile;             //! Digits file for reading
  TString fDigitsOutputFileName;        //! Name of digits file for writing
  TFile*  fDigitsOutputFile;            //! Digits file for writing

  AliTRDdigitsManager* fDigMan; //! digits manager
  AliTRDgeometry* fGeo; //! TRD geometry
  
  Int_t fEventNoInFile;
  Int_t fDigitsLoadedFlag;

  
  AliTRDdigitsTask(const AliTRDdigitsTask&); // not implemented
  AliTRDdigitsTask& operator=(const AliTRDdigitsTask&); // not implemented
  
  ClassDef(AliTRDdigitsTask, 1); 
};

#endif
