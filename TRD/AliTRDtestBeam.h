#ifndef AliTRDtestBeam_h
#define AliTRDtestBeam_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/*
The class to read the test beam 2007 data
*/

//#define MAX_SI 2000

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
  virtual ~AliTRDtestBeam(); // dtor

  Int_t NextEvent();
  
  AliTRDRawStreamTB *GetTRDrawStream(); // needs RawStreamTB

  // silicon
  Short_t GetNSi1() const {return fNSi1;}
  Short_t GetNSi2() const {return fNSi2;}
  
  Int_t GetSi1Address(Int_t i) const {return (i<fNSi1)? fSi1Address[i] : -1;};
  Int_t GetSi2Address(Int_t i) const {return (i<fNSi2)? fSi2Address[i] : -1;};
  
  Int_t GetSi1Charge(Int_t i) const {return (i<fNSi1)? fSi1Charge[i] : -1;};
  Int_t GetSi2Charge(Int_t i) const {return (i<fNSi2)? fSi1Charge[i] : -1;};
  
  Double_t GetX(Int_t n)   const {return (n<2)? fX[n] : -1;}
  Double_t GetY(Int_t n)   const {return (n<2)? fY[n] : -1;}
  Double_t GetQx(Int_t n)  const {return (n<2)? fQx[n] : -1;}
  Double_t GetQy(Int_t n)  const {return (n<2)? fQy[n] : -1;}

  // calo
  Double_t GetCher() const {return fCher;}
  Double_t GetPb() const {return fPb;}

protected:
  
  ifstream *fDataStream;      // input data stream
  
  Bool_t fHeaderIsRead;       // do we read Header ?
  Int_t fEventCount;          // number of events
  Int_t fLimit;               // = 4
  Int_t fCurrent;             // ...

  Int_t fDdlOff;              // offset to DDL data
  Int_t fSiOff;               // offset to Silicon data
  Int_t fQdcOff;              // offset to Cherenkov and LeadGlass data
  Int_t fDdlSize;             // size of DDL data
 
  Char_t *fFileHeader;        // file header data
  Char_t *fEventHeader;       // event header data 
  Char_t *fEventData;         // actual event data

  // silicon data
  
  Short_t fNSi1;              // number of fired silicon pads from Si1
  Short_t fNSi2;              // fired pads in Si2 
  
  Int_t fSi1Address[1270];  // addresses of fires silicon pads Si1
  Int_t fSi2Address[1270];  // addresses if fires silicon pads Si2
  
  Int_t fSi1Charge[1270];   // charge collected on Si1
  Int_t fSi2Charge[1270];   // charge collected on Si2
  
  // reconstructed Silicon data 

  Double_t fX[2];             // x position for Si1 and Si2
  Double_t fY[2];             // y position
  Double_t fQx[2];            // charge on X
  Double_t fQy[2];            // charge on y 

  // cherenkov glass
  Double_t fCher;             // cherenkov signal
  Double_t fPb;               // lead glass signal
  

  // data reading
  
  Int_t Int(Int_t i, Char_t *start);
  Int_t DecodeSi();

  //
  static const Long_t fgkFileHeadSize; //= 544; // ?
  static const Long_t fgkEventHeadSize; // = 68; //?
  static const Long_t fgkLdcHeadSize; // = 68; //?
  static const Long_t fgkEquipHeadSize; // = 28; //
  static const Int_t fgkVmeIn; //=1; //VME event in
  static const Int_t fgkSimIn; //=1; //Si-strips in
  
  //typedef char byte;
  
  //offsets in bytes
  static const Int_t fgkPosRun; // = 20; //run nr. (in file and event header)
  static const Int_t fgkPosLength; // = 0; //event/equip. length
  static const Int_t fgkEqId; // = 8;  //equipment id.
  static const Int_t fgkPosSiOff; // = 12;  //Si data size offset (3 extra words!!!)

  ClassDef(AliTRDtestBeam,1)  // description 

};

#endif 

