#ifndef ALITRDCLUSTERIZERMI_H
#define ALITRDCLUSTERIZERMI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDclusterizerV1.h"

///////////////////////////////////////////////////////
//  Finds and handles cluster (slow simulation)      //
///////////////////////////////////////////////////////

class AliTRDdigitsManager;
class AliTRDparameter;
class AliTRDclusterMI;
class AliTRDclusterizerMI : public AliTRDclusterizerV1 {

 public:

  AliTRDclusterizerMI();
  AliTRDclusterizerMI(const Text_t* name, const Text_t* title);
  virtual ~AliTRDclusterizerMI();
  virtual Bool_t   MakeClusters(); 
  void MakeCluster(Float_t * padSignal, Float_t * pos, Float_t &sigma, Float_t & relpad);
  virtual void AddCluster(Float_t*, int, float, Int_t*, Float_t*, int) {};
  AliTRDclusterMI *  AddCluster();
  void SetCluster(AliTRDclusterMI * cl, Float_t *pos, Int_t det, Float_t amp
		  , Int_t *tracks, Float_t *sig, Int_t iType, Float_t sigmay,Float_t relpos);
 protected:

 private:
  //  AliTRDclusterizerMI &operator=(const AliTRDclusterizerMI &c){;}
  //AliTRDclusterizerMI(const AliTRDclusterizerMI &c){;}
  TObjArray * fClusterContainer;
  virtual Float_t  Unfold(Float_t eps, Int_t plane, Float_t *padSignal);

  ClassDef(AliTRDclusterizerMI,1)           // TRD-Cluster finder, slow simulator
};

#endif
