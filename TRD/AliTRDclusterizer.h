#ifndef ALITRDCLUSTERIZER_H
#define ALITRDCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>

class TFile;
class AliRunLoader;
class AliTRDparameter;

///////////////////////////////////////////////////////
//  Finds and handles cluster                        //
///////////////////////////////////////////////////////

class AliTRDclusterizer : public TNamed {

 public:

  AliTRDclusterizer();
  AliTRDclusterizer(const Text_t* name, const Text_t* title);
  AliTRDclusterizer(const AliTRDclusterizer &c);
  virtual ~AliTRDclusterizer();
  AliTRDclusterizer &operator=(const AliTRDclusterizer &c);

  virtual void    Copy(TObject &c);
  virtual Bool_t  Open(const Char_t *name, Int_t nEvent = 0);
  
  virtual Bool_t  OpenInput(Int_t nEvent = 0);
  virtual Bool_t  OpenOutput();
  virtual Bool_t  MakeClusters() = 0;
  virtual Bool_t  WriteClusters(Int_t det);
  virtual void     SetParameter(AliTRDparameter *par)      { fPar           = par; };
  void     SetVerbose(Int_t v = 1)                 { fVerbose       = v;   };

  AliTRDparameter *GetParameter()                    const { return fPar;          };

 protected:

  AliRunLoader * fRunLoader;       //! Run Loader
  
  TTree           *fClusterTree;   //! Tree with the cluster
  AliTRD          *fTRD;           //! The TRD object
  AliTRDparameter *fPar;           //  TRD digitization parameter object

  Int_t            fEvent;         //! Event number
  Int_t            fVerbose;       //  Sets the verbose level

  ClassDef(AliTRDclusterizer,3)    //  TRD-Cluster manager base class

};

#endif
