#ifndef AliTRDtestBeam_h
#define AliTRDtestBeam_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/*
The class to read the test beam 2007 data
*/

#define MAX_SI 2000

#include "TObject.h"

class AliTRDrawStreamV2;
class AliTRDRawStreamTB;
using namespace std;

class AliTRDtestBeam: public TObject {

public:
  AliTRDtestBeam();       // ctor
  AliTRDtestBeam(const char *filename); // constructor
  //AliTRDtestBeam(const AliTRDtestBeam &tb);  
  //AliTRDtestBeam& operator = (const AliTRDtestBeam& tb) {return *this}
  virtual ~AliTRDtestBeam() {;} // dtor

  Int_t NextEvent();
  
  AliTRDRawStreamTB *GetTRDrawStream(); // needs RawStreamTB

  // silicon
  Short_t GetNSi1() {return fNSi1;}
  Short_t GetNSi2() {return fNSi2;}
  
  Int_t GetSi1Address(Int_t i) {return (i<fNSi1)? fSi1Address[i] : -1;};
  Int_t GetSi2Address(Int_t i) {return (i<fNSi2)? fSi2Address[i] : -1;};
  
  Int_t GetSi1Charge(Int_t i) {return (i<fNSi1)? fSi1Charge[i] : -1;};
  Int_t GetSi2Charge(Int_t i) {return (i<fNSi2)? fSi1Charge[i] : -1;};
  
  Double_t GetX(Int_t n) {return (n<2)? fX[n] : -1;}
  Double_t GetY(Int_t n) {return (n<2)? fY[n] : -1;}
  Double_t GetQx(Int_t n) {return (n<2)? fQx[n] : -1;}
  Double_t GetQy(Int_t n) {return (n<2)? fQy[n] : -1;}

  // calo
  Double_t GetCher() {return fCher;}
  Double_t GetPb() {return fPb;}

protected:
  
  ifstream *fDataStream;
  
  Bool_t fHeaderIsRead;
  Int_t fEventCount;
  Int_t fLimit; // = 4
  Int_t fCurrent;

  Int_t fDdlOff;
  Int_t fSiOff;
  Int_t fQdcOff;
  Int_t fDdlSize;

  Char_t *fFileHeader;
  Char_t *fEventHeader;
  Char_t *fEventData;

  // silicon data
  
  Short_t fNSi1;
  Short_t fNSi2;
  
  Int_t fSi1Address[MAX_SI];
  Int_t fSi2Address[MAX_SI];
  
  Int_t fSi1Charge[MAX_SI];
  Int_t fSi2Charge[MAX_SI];
  
  // reconstructed Silicon data 

  Double_t fX[2];
  Double_t fY[2];
  Double_t fQx[2];
  Double_t fQy[2];

  // cherenkov glass
  Double_t fCher;
  Double_t fPb;
  

  // data reading
  
  Int_t Int(Int_t i, Char_t *start);
  Int_t DecodeSi();

  //
  static const Long_t file_head_size; //= 544; // ?
  static const Long_t event_head_size; // = 68; //?
  static const Long_t ldc_head_size; // = 68; //?
  static const Long_t equip_head_size; // = 28; //
  static const Int_t vme_in; //=1; //VME event in
  static const Int_t sim_in; //=1; //Si-strips in
  
  //typedef char byte;
  
  //offsets in bytes
  static const Int_t pos_run; // = 20; //run nr. (in file and event header)
  static const Int_t pos_length; // = 0; //event/equip. length
  static const Int_t pos_eqid; // = 8;  //equipment id.
  static const Int_t pos_sioff; // = 12;  //Si data size offset (3 extra words!!!)

  ClassDef(AliTRDtestBeam,1)  // description 
};

#endif // AliTRDQADatamaker_H

