#ifndef ALIITSALIGNMENTMODULE_H
#define ALIITSALIGNMENTMODULE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
/* $Author$ */
/* $Date$ */
/* $Name$ */
/* $Header$ */
/*
  $Log$
  Revision 1.1.2.2  2000/03/04 23:39:56  nilsen
  Fixed the logs???

  Revision 1.1.2.1  2000/03/02 20:12:57  nilsen
  A new class useful for ITS detector Alignment studdies
*/
/* $Revision$ */

// Standard C & C++ libraries
#include <TObject.h>

// Standard Root Libraries
#include "TParticle.h"

// ITS libraries
#include "AliITSgeom.h"

class AliITSAlignmentModule : public TObject{
///////////////////////////////////////////////////////////////////////////
//      A track class define exclusively for the use in doing ITS detector
// alignment studdies. Not intended for general use.
// Author: B. S. Nilsen
// Date:   January 17 2000
///////////////////////////////////////////////////////////////////////////

 protected:

    Int_t     findex,flay,flad,fdet;
    TObjArray *ftrksM;
    Double_t  fChi2;
    Double_t  fx0[3],fM[3][3],fangles[3];

 public:

    AliITSAlignmentModule();
    AliITSAlignmentModule(Int_t index,AliITSgeom *gm,
			  Int_t  ntrk,AliITSAlignmentTrack *trk);
    virtual ~AliITSAlignmentModule() {;}; // default destructor OK
    Int_t GetIndex(){return findex;}
    void GetId(Int_t &lay,Int_t &lad,Int_t &det){lay=flay;lad=flad;det=fdet;}
    Double_t * GetTranslationVector(){return fx0;}
    Double_t * GetRotationAngles(){return fangles;}
    Double_t ComputeChi2();
    Double_t GetChi2(){return fChi2;}
    void AlignModule();
    void lnsrch(Int_t npar,Double_t *xold,Double_t fold,Double_t *g,
		Double_t *p,Double_t *x,Double_t &f,Double_t stpmax,
		Int_t &check);
    void MRVMminimization(Int_t npar,Double_t *p,Double_t &fret,Double_t gtol,
			  Int_t &iter);
    void SetByAngles(Double_t *th);

 private:
    void dfMdthx(Double_t dfMx[3][3]);
    void dfMdthy(Double_t dfMx[3][3]);
    void dfMdthz(Double_t dfMx[3][3]);
    void LtoG(Double_t xl[],Double_t xg[]);
    void GtoL(Double_t xg[],Double_t xl[]);
    Double_t Chi2(Double_t p[]);
    void dChi2(Double_t p[],Double_t dChi2[]);
    

    ClassDef(AliITSAlignmentModule,1) // Module class for ITS Alignment


};
#endif
