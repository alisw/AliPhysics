#ifndef ALITRDCLUSTERIZER_H
#define ALITRDCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TObjArray.h>

class TFile;
class TTree;
class AliRunLoader;
class AliTRDparameter;
class AliTRD;

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

  virtual void    Copy(TObject &c) const;
  virtual Bool_t  Open(const Char_t *name, Int_t nEvent = 0);
  
  virtual Bool_t  OpenInput(Int_t nEvent = 0);
  virtual Bool_t  OpenOutput();
  virtual Bool_t  MakeClusters() = 0;
  virtual Bool_t  WriteClusters(Int_t det);
  virtual void     SetParameter(AliTRDparameter *par)      { fPar           = par; };
  void     SetVerbose(Int_t v = 1)                 { fVerbose       = v;   };

  AliTRDparameter *GetParameter()                    const { return fPar;          };

  TObjArray*      RecPoints() {if (!fRecPoints) fRecPoints = new TObjArray(400); return fRecPoints;}
  virtual void    AddCluster(Float_t *pos, Int_t det, Float_t amp, Int_t *tracks
			     , Float_t *sig, Int_t iType);
  void            ResetRecPoints() {if (fRecPoints) fRecPoints->Delete();}


 protected:

  AliRunLoader * fRunLoader;       //! Run Loader
  
  TTree           *fClusterTree;   //! Tree with the cluster
  AliTRDparameter *fPar;           //  TRD digitization parameter object

  TObjArray*       fRecPoints;     //! Array of clusters
  Int_t            fVerbose;       //  Sets the verbose level

  ClassDef(AliTRDclusterizer,3)    //  TRD-Cluster manager base class

};

#endif
