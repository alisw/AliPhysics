#ifndef ALIEVE_MUONChamberData_H
#define ALIEVE_MUONChamberData_H

#include <Reve/Reve.h>

#include <TObject.h>

class AliMUONSegmentation;
class AliMUONGeometryTransformer;

namespace Alieve {

class MUONChamberData : public TObject
{

  MUONChamberData(const MUONChamberData&);            // Not implemented
  MUONChamberData& operator=(const MUONChamberData&); // Not implemented

private:

  static AliMUONSegmentation*        fgSegmentation;  // detector segmentation
  static AliMUONGeometryTransformer* fgTransformer;   // geometry transformer

protected:

  Int_t   fChamberID;                 // number of the chamber, 0 to 13
  Float_t fFrameCoord[26][5];         // detector elements frames
  Int_t   fNDetElem;                  // number of detector elements
  Int_t   fNDigits;                   // number of found digits (times 7)
  Float_t fDigitBuffer[7*4096];       // digits coordinates, etc.
  Float_t fChamberBox[6];             // chamber envelope box

public:

  MUONChamberData(Int_t chamber);
  virtual ~MUONChamberData();

  void DropData();

  void     Init(Int_t chamber);
  void     RegisterDigit(Int_t detElemId, Int_t cathode, Int_t ix, Int_t iy, Int_t charge);
  Float_t* GetFrameCoord(Int_t detElemId) { return fFrameCoord[detElemId]; };
  Int_t    GetNDetElem() { return fNDetElem; };
  Int_t    GetNDigits() { return fNDigits/7; };
  Float_t* GetDigitBuffer(Int_t pos) { return &fDigitBuffer[7*pos]; };
  Float_t* GetChamberBox() { return &fChamberBox[0]; };

  ClassDef(MUONChamberData,1);     // class with data for one chamber

};

}

#endif
