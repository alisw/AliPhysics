#ifndef ALITRDTESTG4_H
#define ALITRDTESTG4_H
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
class AliTRDtestG4 : public AliTRD {

 public:

  AliTRDtestG4();
  AliTRDtestG4(const char *name, const char *title);
  virtual ~AliTRDtestG4();

  virtual void     Init();
  virtual Int_t    IsVersion() const          { return 1;      }

  virtual void     AddAlignableVolumes() const;
  virtual void     CreateGeometry();
  virtual void     CreateMaterials();
  virtual void     CreateTRhit(Int_t det);

  virtual void     StepManager();

          void     SetStepSize(Double_t s)    { fStepSize = s; }
          void     SetTR(Bool_t tr)           { fTRon = tr;    }
          void     SetScaleG4(Float_t f)      { fScaleG4 = f;  }

          Bool_t   GetTR() const              { return fTRon;  }
  AliTRDsimTR     *GetTRDsim() const          { return fTR;    }

 protected:

          Bool_t   fTRon;               //  Switch for TR simulation
  AliTRDsimTR     *fTR;                 //  TR simulator

          Double_t fStepSize;           //  Used for the fixed step size
          Float_t  fWion;               //  Ionization potential
	  Float_t  fScaleG4;            //  G4 scaling factor of de/dx rel. to g3 

 private:

  AliTRDtestG4(const AliTRDtestG4 &trd);
  AliTRDtestG4 &operator=(const AliTRDtestG4 &trd);

  ClassDef(AliTRDtestG4,2)                  //  Transition Radiation Detector (test version for G4 simulations)

};

#endif
