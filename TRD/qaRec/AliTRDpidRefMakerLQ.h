#ifndef ALITRDPIDREFMAKERLQ_H
#define ALITRDPIDREFMAKERLQ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDpidRefMakerLQ.h 34125 2009-08-06 09:35:40Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for building reference data for PID                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TH1;
class TH2;
class TH3;
class TPrincipal;
class TLinearFitter;

class AliTRDpidRefMakerLQ : public AliTRDrecoTask {

public:
  AliTRDpidRefMakerLQ();
  ~AliTRDpidRefMakerLQ();

  Bool_t   BuildLQReferences(const Char_t *filename="TRDpidLQ.root", const Char_t *dir =".");
  Bool_t   BuildNNReferences(const Char_t *filename="TRDpidNN.root", const Char_t *dir =".") const;
 
  static Double_t Estimate2D2(TH2 * const h, Float_t &x, Float_t &y);
  static Double_t Estimate2D1(TH2 * const h, Float_t &x, Float_t &y, const Float_t &dCT
                            , const Float_t &rmin, const Float_t &rmax);
  //     Double_t Estimate3D2(TH3 * const h, Float_t &x, Float_t &y, Float_t &z);

private:
  AliTRDpidRefMakerLQ(const AliTRDpidRefMakerLQ &ref);
  AliTRDpidRefMakerLQ& operator=(const AliTRDpidRefMakerLQ &ref);
  Int_t    CheckProdDirTree(const Char_t *dir=".");
  void     Prepare2D();
  void     Reset();
  void     SaveReferences(const Int_t mom, const char *fn);
 
private:
          TPrincipal    *fPrinc[5];   // Used for principal component analysis
  static TLinearFitter *fgFitter2D2; // Working object for linear fitter
  static TLinearFitter *fgFitter2D1; // Working object for linear fitter
          TH2           *fH2dEdx[5];  // dE/dx data holders
          TH1           *fH1TB[2];    // Max time bin data holders

  ClassDef(AliTRDpidRefMakerLQ, 3)  // Reference histograms builder

};

#endif

