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
class AliTRDcluster;
class AliTRDclusterizerMI : public AliTRDclusterizerV1 {

 public:

  AliTRDclusterizerMI();
  AliTRDclusterizerMI(const Text_t* name, const Text_t* title);
  virtual ~AliTRDclusterizerMI();
  virtual Bool_t   MakeClusters(); 
  void MakeCluster(Double_t * padSignal, Double_t * pos, Double_t &sigma, Double_t & relpad);
  virtual AliTRDcluster  * AddCluster(Double_t *, Int_t , Double_t , Int_t *
				      , Double_t *, Int_t , Float_t = 0){return 0;}
  AliTRDclusterMI *  AddCluster();
  void SetCluster(AliTRDclusterMI * cl, Double_t *pos, Int_t det, Double_t amp
		  , Int_t *tracks, Double_t *sig, Int_t iType, Double_t sigmay,Double_t relpos);
 protected:

 private:
  //  AliTRDclusterizerMI &operator=(const AliTRDclusterizerMI &c){;}
  //AliTRDclusterizerMI(const AliTRDclusterizerMI &c){;}
  TObjArray * fClusterContainer;
  virtual Double_t Unfold(Double_t eps, Int_t plane, Double_t *padSignal);

  ClassDef(AliTRDclusterizerMI,1)           // TRD-Cluster finder, slow simulator

};

#endif
