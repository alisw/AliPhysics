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

class AliTOFcalib:public TTask{
public:
  AliTOFcalib();          // ctor
  AliTOFcalib(AliTOFGeometry *geom);
  AliTOFcalib(const AliTOFcalib & calib); // copy constructor
  AliTOFcalib& operator=(const AliTOFcalib & calib); // assignment operator
  virtual ~AliTOFcalib() ; // dtor
  Int_t NSector()const {return fNSector;}
  Int_t NPlate()const {return fNPlate;}
  Int_t NStripA()const {return fNStripA;}
  Int_t NStripB()const {return fNStripB;}
  Int_t NStripC()const {return fNStripC;}
  Int_t NpadZ()const {return fNpadZ;}
  Int_t NpadX()const {return fNpadX;}
  AliTOFCal * GetTOFCalArray() const {return fTOFCal;}
  AliTOFCal * GetTOFCalSimArray() const {return fTOFSimCal;}
  TH1F * GetTOFSimToT() const {return fTOFSimToT;}
  void SelectESD(AliESD *event);   
  void CombESDId();
  void CalibrateESD();
  TH1F* Profile(Int_t i); 
  Int_t NChannels()const{return fNChannels;}
  TF1* SetFitFunctions(TH1F* histo);
  void CorrectESDTime();// useless method, kept to make Chiara happy
  void CorrectESDTime(AliESD *event);
  // Methods to retrieve/write parameters from/on CDB
  void WriteSimParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void WriteSimParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun, AliTOFCal *cal, TH1F *histo);
  void ReadSimParFromCDB(Char_t *sel, Int_t nrun);
  void WriteParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun, AliTOFCal *cal);
  void WriteParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  void ReadParFromCDB(Char_t *sel, Int_t nrun);
  Int_t GetIndex(Int_t *detId); // Get channel index for Calibration 

 public: 
  class AliTOFArray : public TObject {
  public:
    AliTOFArray(): TObject(),fSize(0),fArray(0x0){}
    AliTOFArray(Int_t size) :
      TObject(),
      fSize(size),
      fArray(new TArrayF*[size]) {
    } 
    AliTOFArray(const AliTOFArray & source):
      TObject(){ // copy constructor
      this->fSize= source.fSize;
      this->fArray= source.fArray;
    };

    AliTOFArray& operator=(const AliTOFArray & source) { // assignment operator
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
    virtual ~AliTOFArray() {
      delete [] fArray;
    }
    
  private:
    
    Int_t fSize;       // Size of the array of TArrayFs
    TArrayF ** fArray; //[fSize]};
    
  };


private:
  static const Int_t fgkchannel;  // max number of entries per channel 
  Int_t fNChannels; // number of TOF channels
  Int_t fNSector;  // number of TOF sectors
  Int_t fNPlate;   // number of TOF plates
  Int_t fNStripA;  // number of TOF strips A
  Int_t fNStripB;  // number of TOF strips B
  Int_t fNStripC;  // number of TOF strips C
  Int_t fNpadZ;    // number of TOF pads Z
  Int_t fNpadX;    // number of TOF pads X
  Int_t fNevents;           // number of events
  TObjArray * fESDsel;   // selected ESD tracks for calibration
  AliTOFArray *fArrayToT;   // array for ToT values
  AliTOFArray *fArrayTime;  // array for Time values
  AliTOFCal *fTOFCal;       // array of AliTOFChannels storing calib parameters
  AliTOFCal *fTOFSimCal;       // array of AliTOFChannels storing calib parameters
  TH1F *fTOFSimToT;        // histo with realistic ToT signal from TB Data
  ClassDef(AliTOFcalib,1);
};

#endif // AliTOFcalib_H

