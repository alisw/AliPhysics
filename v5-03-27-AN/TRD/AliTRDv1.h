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

class AliTRDsimTR;

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

          void     SetStepSize(Double_t s)    { fStepSize = s; }
          void     SetTR(Bool_t tr)           { fTRon = tr;    }

          Bool_t   GetTR() const              { return fTRon;  }
  AliTRDsimTR     *GetTRDsim() const          { return fTR;    }

 protected:

          Bool_t   fTRon;               //  Switch for TR simulation
  AliTRDsimTR     *fTR;                 //  TR simulator

          Double_t fStepSize;           //  Used for the fixed step size
          Float_t  fWion;               //  Ionization potential

 private:

  AliTRDv1(const AliTRDv1 &trd);
  AliTRDv1 &operator=(const AliTRDv1 &trd);

  ClassDef(AliTRDv1,8)                  //  Transition Radiation Detector version 1 (slow simulator)

};

#endif
