// $Header$

#ifndef ALIEVE_PMDModule_H
#define ALIEVE_PMDModule_H

#include <Reve/Reve.h>
#include <Reve/QuadSet.h>

#include <TObject.h>
#include <TObjArray.h>
#include <TTree.h>

namespace Alieve {

class PMDModule : public Reve::QuadSet
{
private:
  PMDModule(const PMDModule&);            // Not implemented
  PMDModule& operator=(const PMDModule&); // Not implemented

  void RectGeomCellPos(Int_t ism, Int_t irow, Int_t icol,
		       Float_t &xpos, Float_t &ypos);
  void GenerateBox(Int_t ism, Float_t &xism, Float_t &yism,
		   Float_t &dxism, Float_t &dyism);

protected:
  Float_t       fX, fY, fZ;
  Int_t         fNPads;

  static const Float_t fgkRad;
  static const Float_t fgkSqRoot3;
  static const Float_t fgkZpos;


public:
  PMDModule();
  virtual ~PMDModule() {}

  Int_t GetNPads() const { return fNPads; }

  void DisplayInit(Int_t ism);
  void DisplayDigitsData(Int_t ism, TTree *pmdt);
  void DisplayRawData(Int_t ism,   TObjArray *ddlcont);
  void SetPosition(Float_t x, Float_t y, Float_t z);

  ClassDef(PMDModule, 1);
}; // endclass PMDModule

}

#endif
