#ifndef ALITRDV1_H
#define ALITRDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Manager and hits classes for set: TRD version 1                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

// Energy spectrum of the delta-rays 
Double_t Ermilova(Double_t *x, Double_t *par);
Double_t IntSpecGeant(Double_t *x, Double_t *par);
 
#include "AliTRD.h"

class TF1;
class TTree;
class TFile;

class AliTRDsim;

//_____________________________________________________________________________
class AliTRDv1 : public AliTRD {

 public:

  AliTRDv1();
  AliTRDv1(const char *name, const char *title);
  virtual ~AliTRDv1();

  virtual void     Init();
  virtual Int_t    IsVersion() const          { return 1;      }

  virtual void     AddAlignableVolumes() const;
  virtual void     CreateGeometry();
  virtual void     CreateMaterials();
  virtual void     CreateTRhit(Int_t det);

  virtual void     StepManager();
          void     StepManagerErmilova();
          void     StepManagerGeant();
          void     StepManagerFixedStep();
          void     SelectStepManager(Int_t t);

          void     SetStepSize(Double_t s)    { fStepSize = s; }
          void     SetTR(Bool_t kTRUE)        { fTRon = kTRUE; }

          Bool_t   GetTR() const              { return fTRon;  }
  AliTRDsim       *GetTRDsim() const          { return fTR;    }

 protected:

          void    *StepManagerEntity();

          Bool_t   fTRon;               //  Switch for TR simulation
  AliTRDsim       *fTR;                 //  TR simulator

          Int_t    fTypeOfStepManager;  //  Type of Step Manager.
          Double_t fStepSize;           //  Used for the fixed step size

 private:
  AliTRDv1(const AliTRDv1 &trd);
  AliTRDv1 &operator=(const AliTRDv1 &trd);

          Double_t BetheBloch(Double_t bg);
          Double_t BetheBlochGeant(Double_t bg);
  
          TF1     *fDeltaE;             //  Energy distribution of the delta-electrons (Ermilova)
          TF1     *fDeltaG;             //  Energy distribution of the

          Float_t  fTrackLength0;       //  Save the track length at chamber entrance  
          Int_t	   fPrimaryTrackPid;    //  Save the id of the primary track  

  ClassDef(AliTRDv1,5)                  //  Transition Radiation Detector version 1 (slow simulator)

};

#endif
