#ifndef ALITRDTRACKERFITTER_H
#define ALITRDTRACKERFITTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Track fitter                                                             //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef AliRIEMAN_H
#include <AliRieman.h>
#endif

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

// implementations in the future
class TObject;
class TLinearFitter;
class TTreeSRedirector;

class AliTRDcluster;
class AliTRDseedV1;

class AliTRDtrackerFitter : public TObject {
 
 public:

	AliTRDtrackerFitter();
	AliTRDtrackerFitter(const AliTRDtrackerFitter &);
	AliTRDtrackerFitter& operator=(const AliTRDtrackerFitter &);
	virtual ~AliTRDtrackerFitter();

	void            Copy(TObject &f) const;
	AliRieman      *GetRiemanFitter() const { return fRieman1; }
	void            GetHyperplaneFitChi2(Double_t *chisquares) const;
	void            GetHyperplaneFitResults(Double_t *params) const;
	Double_t        GetRiemanCurvature() const { return fRieman1->GetC(); }

	void            FitRieman(AliTRDcluster **cl, Int_t nLayers);
	void            FitRieman(AliTRDseedV1 *cseed, Int_t *planes=0x0);
	Double_t        FitHyperplane(AliTRDseedV1 *cseed, Double_t chi2ZF, Double_t zval);
	void            Reset();
	void            SetDebugStream(TTreeSRedirector *debug){fDebugStream = debug;}
	void            SetLayers(const Int_t nlayers) {fNlayers = nlayers;}

 private:

	// Linear fitters in planes
	TLinearFitter    *fFitterTC;       // Fitting with tilting pads - kz fixed - kz= Z/x, + vertex const
	TLinearFitter    *fFitterT2;       // Fitting with tilting pads - kz not fixed
	AliRieman        *fRieman1;        // Rieman fitter
	AliRieman        *fRieman2;        // Rieman fitter
	// Linear Fitter chisquares
	Double_t          fChi2TR;         // Chisquared
	Double_t          fChi2TC;         // Chisquared
	// Hyperplane Fit results
	Double_t          fCR;             // Track radius
	Double_t          fCC;             // Track curvature
	Double_t          fDca;            // DCA       
	// Helping Parameters
	Double_t          fDzmf;           // Something
	Double_t          fZmf;            // Something else
	Int_t             fNlayers;        // Some layers

	// Debug
	TTreeSRedirector *fDebugStream;    //! debugger

	ClassDef(AliTRDtrackerFitter,1)    // TRD track fitter
};
#endif // ALITRDTRACKERFITTER_H
