#ifndef ALIMUONTRACKRECONSTRUCTOR_H
#define ALIMUONTRACKRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class AliMUONTrackReconstructor : public TObject {
 public:
    //
    // Track Reconstruction
    AliMUONTrackReconstructor();
    virtual       ~AliMUONTrackReconstructor(){;}
    void     Init(Double_t &, Double_t &, Double_t &);
    void     Reconst(Int_t &,Int_t &,Int_t,Int_t &,Int_t&,Int_t&, Option_t *option,Text_t *filename);
    void     Reconst2(Int_t &,Int_t &,Int_t &);
    void     FinishEvent();
    void     Close();
    void     SetCutPxz(Double_t p) {fSPxzCut=p;}
    void     SetSigmaCut(Double_t p) {fSSigmaCut=p;}
    void     SetXPrec(Double_t p) {fSXPrec=p;}
    void     SetYPrec(Double_t p) {fSYPrec=p;}
    Double_t GetCutPxz() {return fSPxzCut;}
    Double_t GetSigmaCut() {return fSSigmaCut;}
    Double_t GetXPrec() {return fSXPrec;}
    Double_t GetYPrec() {return fSYPrec;}
 private:
//  Parameters for reconstruction program
    Double_t fSPxzCut;        // Pxz cut  (GeV/c) to begin the track finding
    Double_t fSSigmaCut;      // Number of sig. delimiting the searching areas
    Double_t fSXPrec;         // Chamber precision in X (cm) 
    Double_t fSYPrec;         // Chamber precision in Y (cm)
    Text_t *fFileName;        //! ?????????
    ClassDef(AliMUONTrackReconstructor,1)  // Interface to muon tracking code
	};
	

#endif
