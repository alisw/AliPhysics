#ifndef ALITOFCALIB_H
#define ALITOFCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//  class for TOF calibration:: simulation of uncalibrated data //
//////////////////////////////////////////////////////////////////

#include "TTask.h"
#include "TH1F.h"
#include "AliTOFChannel.h"
#include "TClonesArray.h"
#include "TList.h"
#include "AliTOFCal.h"
#include "AliTOFGeometry.h"
#include "AliESD.h"
#include "AliRunLoader.h"

class AliTOFcalib:public TTask{
public:
  AliTOFcalib();          // ctor
  AliTOFcalib(char* headerFile, Int_t nEvents=0) ; 
  AliTOFcalib(const AliTOFcalib & calib);
  void Init();
  virtual ~AliTOFcalib() ; // dtor
  Int_t NSector()const {return fNSector;}
  Int_t NPlate()const {return fNPlate;}
  Int_t NStripA()const {return fNStripA;}
  Int_t NStripB()const {return fNStripB;}
  Int_t NStripC()const {return fNStripC;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  TClonesArray * DecalibrateDigits(TClonesArray *digits);
  void SelectESD(AliESD *event, AliRunLoader * rl); 
  void CombESDId();
  void CalibrateESD();
  TH1F* Profile(Int_t i); 
  Int_t Size()const{return fsize;}
  void SetFitFunctions();
  TF1* SetFitFunctions(TH1F* histo);
  TList* GetFitFunctions() {return flistFunc;}
  TH1F** GetHistosToT() {return fhToT;}
  void SetHistos();
  void CorrectESDTime();
  void CorrectESDTime(AliESD *event, AliRunLoader *rl);
  void WriteOnCDB();
  void ReadFromCDB(Char_t *sel, Int_t nrun);
  Int_t GetIndex(Int_t *detId);

 public: 
  class AliTOFArray : public TObject {
  public:
    AliTOFArray(): TObject(),fSize(0),fArray(0x0){}
    AliTOFArray(Int_t size) :
      TObject(),
      fSize(size),
      fArray(new TArrayF*[size]) {
    } 
    Int_t GetSize() const {return fSize;}
    void AddArray(Int_t pos, TArrayF * parr) {
      if (pos>-1 && pos < fSize)
	fArray[pos] = parr;
      //else
      //AliError("Index  out of range");  
    }
    TArrayF * GetArray(Int_t pos) {
      TArrayF * parr = 0x0;
      if  (pos>-1 && pos < fSize)
	parr = fArray[pos];
      //else
        //AliError("Index out of range");  
      return parr;
    }
    virtual ~AliTOFArray() {
      delete [] fArray;
    }
    
  private:
    
    Int_t fSize;       // Size of the array of TArrayFs
    TArrayF ** fArray; //[fSize]};
    
  };


private:
  static const Int_t fgkchannel;  // max number of entries per channel 
  static const Char_t * ffile[6]; // spectra
  Int_t fsize;     // number of channels
  Int_t fNSector;  // number of TOF sectors
  Int_t fNPlate;   // number of TOF plates
  Int_t fNStripA;  // number of TOF strips A
  Int_t fNStripB;  // number of TOF strips B
  Int_t fNStripC;  // number of TOF strips C
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X
  TObjArray * fESDsel;   // selected ESD tracks for calibration
  TList *flistFunc;      // functions for simulated Time Slewing spectra
  TH1F* fhToT[6];           // simulated ToT distributions
  Float_t fMaxToT[6];       // max simulated ToT
  Float_t fMinToT[6];       // min simulated ToT
  Float_t fMaxToTDistr[6];  // max value in the ToT distributions
  AliTOFArray *fArrayToT;   // array for ToT values
  AliTOFArray *fArrayTime;  // array for Time values
  Int_t fNevents;           // number of events
  AliTOFCal *fTOFCal;       // array of AliTOFChannels storing calib parameters
  TString fT0File ;         // output file;
  TString fHeadersFile;     // input file
  AliTOFGeometry *fGeom;    // AliTOFgeometry pointer
  ClassDef(AliTOFcalib,1);
};

#endif // AliTOFcalib_H

