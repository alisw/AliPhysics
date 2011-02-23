#ifndef AliOADBCentrality_H
#define AliOADBCentrality_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     OADB class for run dependent centrality scaling factors and 
//     data for centrality determination
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>


class AliOADBCentrality : public TNamed {

 public :
  AliOADBCentrality();
  AliOADBCentrality(char* name);
  virtual ~AliOADBCentrality();
  Float_t V0MScaleFactor() const  {return fV0MScaleFactor;}
  Float_t SPDScaleFactor() const  {return fSPDScaleFactor;}
  Float_t TPCScaleFactor() const  {return fTPCScaleFactor;}
  TH1F*   V0hist()         const  {return ((TH1F*) (Hists1D()->FindObject("hmultV0_percentile")));}
  TH1F*   TPChist()        const  {return ((TH1F*) (Hists1D()->FindObject("hNtracks_percentile")));}
  TH1F*   SPDhist()        const  {return ((TH1F*) (Hists1D()->FindObject("hNclusters1_percentile")));}
  TH2F*   ZEMvsZDChist()   const  {return ((TH2F*) (Hists2D()->FindObject("hEzemvsEzdc_all_percentile")));}
  TList*  Hists1D()        const  {return f1DHistos;}
  TList*  Hists2D()        const  {return f2DHistos;}
  void    SetScaleFactors(Float_t v0m, Float_t spd, Float_t tpc)
      {fV0MScaleFactor = v0m; fSPDScaleFactor = spd; fTPCScaleFactor = tpc;}
  void    SetHistReferences(TList* l1, TList* l2)
      {f1DHistos = l1; f2DHistos = l2;}
 private:
  AliOADBCentrality(const AliOADBCentrality& cont); 
  AliOADBCentrality& operator=(const AliOADBCentrality& cont);

 private:
  Float_t fV0MScaleFactor; // V0  scale factor
  Float_t fSPDScaleFactor; // SPD scale factor
  Float_t fTPCScaleFactor; // TPC scale factor
  TList*    f1DHistos; // Reference to list of 1D Centrality histos 
  TList*    f2DHistos; // Reference to list of 2D Centrality histos
  ClassDef(AliOADBCentrality, 1);
};

#endif
