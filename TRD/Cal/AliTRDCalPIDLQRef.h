#ifndef ALITRDCALPIDLQREF_H
#define ALITRDCALPIDLQREF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for 2-dim reference histograms                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class TH2;
class TH3;
class TPrincipal;
class TLinearFitter;

class AliTRDCalPIDLQRef : public TObject {

 public:

  AliTRDCalPIDLQRef();
  AliTRDCalPIDLQRef(const AliTRDCalPIDLQRef &ref);
  ~AliTRDCalPIDLQRef();
  AliTRDCalPIDLQRef& operator=(const AliTRDCalPIDLQRef &ref);

  TPrincipal*	GetPrincipal(const Int_t spec){return (spec>=0 && spec<5) ? fPrinc[spec] : 0x0;}
  Double_t	Estimate2D2(TH2 *h, Float_t &x, Float_t &y);
  Double_t	Estimate2D1(TH2 *h, Float_t &x, Float_t &y, Float_t &dCT, Float_t &rmin, Float_t &rmax);
  //Double_t	  Estimate3D2(TH3 *h, Float_t &x, Float_t &y, Float_t &z);
  void		Prepare2DReferences();
  void		Reset();
  void      	SaveReferences(const Int_t mom, const char *fn);
	
 private:

  TPrincipal	*fPrinc[5];      // Used for principal component analysis
  TLinearFitter	*fFitter2D2;     // Working object for linear fitter
  TLinearFitter	*fFitter2D1;     // Working object for linear fitter

  TH2 		*h2dEdx[5];      // Data holders

  ClassDef(AliTRDCalPIDLQRef, 1) // Reference histograms builder

};

#endif

