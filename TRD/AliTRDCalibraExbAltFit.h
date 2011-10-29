#ifndef ALITRDCALIBRAEXBALTFIT_H
#define ALITRDCALIBRAEXBALTFIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalibraExbAltFit.h 46327 2011-01-10 13:29:56Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for online calibration                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//#include "TObjArray.h"
#include "TObject.h"

//class TVectorD;
class TObjArray;
class TH2S;
class TTreeSRedirector;

class AliTRDCalibraExbAltFit : public TObject {

 public:

  AliTRDCalibraExbAltFit();
  AliTRDCalibraExbAltFit(const AliTRDCalibraExbAltFit &ped);
  AliTRDCalibraExbAltFit(const TObjArray &obja);
  virtual ~AliTRDCalibraExbAltFit();
  virtual Long64_t Merge(const TCollection* list);
  virtual void Copy(TObject &c) const;

  AliTRDCalibraExbAltFit& operator = (const  AliTRDCalibraExbAltFit &source);

  void            Update(Int_t detector, Float_t tnp, Float_t pars1);
  void            FillPEArray();
  void            Add(const AliTRDCalibraExbAltFit *ped);
  TH2S           *GetFitterHisto(Int_t detector, Bool_t force=kFALSE);
  TH2S           *GetFitterHistoForce(Int_t detector);
  TH2S           *GetFitterHistoNoForce(Int_t detector) const   { return (TH2S*)fFitterHistoArray.UncheckedAt(detector);};
  Bool_t          GetParam(Int_t detector, TVectorD *param);
  Bool_t          GetError(Int_t detector, TVectorD *error);

  TObjArray      *GetPArray()                    { return &fFitterPArray;       };
  TObjArray      *GetEArray()                    { return &fFitterEArray;       };
  TObjArray       GetHistoArray() const          { return fFitterHistoArray;    };

 private:
   
  Int_t           fVersion;                 // Version of the object

  TObjArray       fFitterHistoArray;  // TObjArray of histo2D for debugging  Fitters
  TObjArray       fFitterPArray;      // Array of result parameters from  fitters for the detectors
  TObjArray       fFitterEArray;      // Array of result errors from  fitters for the detectors

  
  ClassDef(AliTRDCalibraExbAltFit,1)  // Online ExB Calibration

};



#endif

