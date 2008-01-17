// $Header$

#ifndef ALIEVE_PMDModule_H
#define ALIEVE_PMDModule_H

#include <TEveUtil.h>
#include <TEveQuadSet.h>

#include <TObject.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TH1F.h>

namespace Alieve {

class PMDModule : public TEveQuadSet
{
private:
  PMDModule(const PMDModule&);            // Not implemented
  PMDModule& operator=(const PMDModule&); // Not implemented

  void RectGeomCellPos(Int_t ism, Int_t irow, Int_t icol,
		       Float_t &xpos, Float_t &ypos);
  void GenerateBox(Int_t ism, Float_t &xism, Float_t &yism,
		   Float_t &dxism, Float_t &dyism);

protected:
  TH1F*                fH1;
  Float_t              fX, fY, fZ;
  Int_t                fNPads;
  Int_t                fAdc;

  static const Float_t fgkRad;
  static const Float_t fgkSqRoot3;
  static const Float_t fgkZpos;

  static Int_t         fPreTotPads;
  static Int_t         fCpvTotPads;
  static Int_t         fPreTotAdc;
  static Int_t         fCpvTotAdc;


public:
  PMDModule();
  virtual ~PMDModule() { delete fH1; }

  Int_t GetPRETotPads() const { return fPreTotPads; }
  Int_t GetCPVTotPads() const { return fCpvTotPads; }
  Int_t GetNPads()      const { return fNPads; }
  Int_t GetPRETotAdc()  const { return fPreTotAdc; }
  Int_t GetCPVTotAdc()  const { return fCpvTotAdc; }
  Int_t GetAdc()        const { return fAdc; }
  TH1F *GetHisto()      const { return fH1;}

  void DisplayInit(Int_t ism);
  void DisplayDigitsData(Int_t ism, TTree *pmdt);
  void DisplayRawData(Int_t ism,   TObjArray *ddlcont);
  void SetPosition(Float_t x, Float_t y, Float_t z);

  ClassDef(PMDModule, 1);
}; // endclass PMDModule

}

#endif
