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
class AliTRDcluster;
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
  virtual void    SetParameter(AliTRDparameter *par)      { fPar           = par; };
  void            SetVerbose(Int_t v = 1)                 { fVerbose       = v;   };

  AliTRDparameter *GetParameter()                   const { return fPar;          };

  TObjArray*      RecPoints() {if (!fRecPoints) fRecPoints = new TObjArray(400); return fRecPoints;}
  virtual AliTRDcluster  * AddCluster(Double_t *pos, Int_t timebin, Int_t det, Double_t amp, Int_t *tracks
			     , Double_t *sig, Int_t iType, Float_t center = 0);
  void            ResetRecPoints() {if (fRecPoints) fRecPoints->Delete();}

 protected:

   Double_t CalcXposFromTimebin(Float_t timebin, Float_t vdrift);
       
  AliRunLoader    *fRunLoader;     //! Run Loader
  
  TTree           *fClusterTree;   //! Tree with the cluster
  AliTRDparameter *fPar;           //  TRD digitization parameter object

  TObjArray*       fRecPoints;     //! Array of clusters
  Int_t            fVerbose;       //  Sets the verbose level

  ClassDef(AliTRDclusterizer,3)    //  TRD-Cluster manager base class

};

#endif
