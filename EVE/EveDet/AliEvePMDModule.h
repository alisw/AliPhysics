// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEvePMDModule_H
#define AliEvePMDModule_H

#include <TEveQuadSet.h>

class TH1F;
class TTree;

class AliEvePMDModule : public TEveQuadSet
{
public:
  AliEvePMDModule();
  virtual ~AliEvePMDModule();

  static Int_t GetPRETotPads() { return fgPreTotPads; }
  static Int_t GetCPVTotPads() { return fgCpvTotPads; }
  static Int_t GetPRETotAdc()  { return fgPreTotAdc;  }
  static Int_t GetCPVTotAdc()  { return fgCpvTotAdc;  }

  Int_t GetNPads()      const { return fNPads; }
  Int_t GetAdc()        const { return fAdc;   }
  TH1F *GetHisto()      const { return fH1;    }

  void DisplayInit(Int_t ism);
  void DisplayDigitsData(Int_t ism, TTree *pmdt);
  void DisplayRawData(Int_t ism,   TObjArray *ddlcont);
  void SetPosition(Float_t x, Float_t y, Float_t z);

protected:
  TH1F*        fH1;         // histogram
  Float_t      fX, fY, fZ;  // coordinates
  Int_t        fNPads;      // number of pads
  Int_t        fAdc;        // ad count

  static const Float_t fgkRad;      // fooconst
  static const Float_t fgkSqRoot3;  // fooconst
  static const Float_t fgkZpos;     // position of PMD

  static Int_t fgPreTotPads; // total pre pads
  static Int_t fgCpvTotPads; // total cpv pads
  static Int_t fgPreTotAdc;  // total pre signal
  static Int_t fgCpvTotAdc;  // total cpv signal

private:
  void RectGeomCellPos(Int_t ism, Int_t irow, Int_t icol,
		       Float_t &xpos, Float_t &ypos);
  void GenerateBox(Int_t ism, Float_t &xism, Float_t &yism,
		   Float_t &dxism, Float_t &dyism);

  AliEvePMDModule(const AliEvePMDModule&);            // Not implemented
  AliEvePMDModule& operator=(const AliEvePMDModule&); // Not implemented

  ClassDef(AliEvePMDModule, 0);
}; // endclass AliEvePMDModule

#endif
