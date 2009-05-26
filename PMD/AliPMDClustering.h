#ifndef ALIPMDCLUSTERING_H
#define ALIPMDCLUSTERING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Header File : PMDClustering.h, Version 00          //
//                                                     //
//  Date   : September 26 2002                         //
//                                                     //
//  clustering code for alice pmd                      //
//                                                     //
//-----------------------------------------------------//
/*-----------------------------------------------------------------------*/
#include "Rtypes.h"

class TNtuple;
class TObjArray;
class AliPMDClustering: public TObject
{

 public:
  AliPMDClustering(){};
  virtual ~AliPMDClustering(){};

  virtual void DoClust(Int_t idet, Int_t ismn, Int_t celltrack[][96],
		       Int_t cellpid[][96], Double_t celladc[][96],
		       TObjArray *pmdisocell, TObjArray *pmdcont) = 0;

  virtual void SetEdepCut(Float_t decut) = 0;

  ClassDef(AliPMDClustering,7) // Does clustering for PMD
};
#endif
