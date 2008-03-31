// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveMUONChamberData_H
#define AliEveMUONChamberData_H

#include <TObject.h>

class AliMUONGeometryTransformer;

class AliEveMUONChamberData : public TObject
{
public:

  AliEveMUONChamberData(Int_t chamber);
  virtual ~AliEveMUONChamberData();

  void     DropData();

  void     Init(Int_t chamber);

  void     RegisterDigit(Int_t detElemId, Int_t cathode, Int_t ix, Int_t iy, Int_t charge);
  void     RegisterCluster(Int_t detElemId, Int_t cathode, Float_t x, Float_t y, Float_t z, Float_t charge);
  void     RegisterHit(Int_t detElemId, Float_t x, Float_t y, Float_t z);

  Float_t* GetFrameCoord(Int_t detElemId) { return fFrameCoord[detElemId]; };

  Int_t    GetNDetElem()  const { return fNDetElem;    };
  Int_t    GetNDigits()   const { return fNDigits/7;   };
  Int_t    GetNClusters() const { return fNClusters/5; };
  Int_t    GetNHits()     const { return fNHits/3;     };

  Float_t* GetDigitBuffer(Int_t pos)   { return &fDigitBuffer[7*pos];   };
  Float_t* GetClusterBuffer(Int_t pos) { return &fClusterBuffer[5*pos]; };
  Float_t* GetHitBuffer(Int_t pos)     { return &fHitBuffer[3*pos];     };

  Float_t* GetChamberBox() { return &fChamberBox[0]; };


protected:

  Int_t   fChamberID;                 // number of the chamber, 0 to 13
  Float_t fFrameCoord[26][5];         // detector elements frames
  Int_t   fNDetElem;                  // number of detector elements
  Int_t   fNDigits;                   // number of found digits (times 7)
  Int_t   fNClusters;                 // number of found rec points
  Int_t   fNHits;                     // number of simulation hits
  Float_t fDigitBuffer[7*4096];       // digits coordinates, etc.
  Float_t fClusterBuffer[5*256];      // cluster coordinates, etc.
  Float_t fHitBuffer[3*256];          // hits coordinates
  Float_t fChamberBox[6];             // chamber envelope box


private:

  static AliMUONGeometryTransformer* fgTransformer;   // geometry transformer

  AliEveMUONChamberData(const AliEveMUONChamberData&);            // Not implemented
  AliEveMUONChamberData& operator=(const AliEveMUONChamberData&); // Not implemented

  ClassDef(AliEveMUONChamberData, 0);     // class with data for one chamber
};

#endif
