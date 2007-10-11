#ifndef ALIITSCLUSTERPARAM_H
#define ALIITSCLUSTERPARAM_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////
//                                                //
//  ITS cluster error and shape parameterization  //
//  andrea.dainese@lnl.infn.it                    //
////////////////////////////////////////////////////


#include <TObject.h>
#include "AliITSRecPoint.h"

//class TTree;

//_____________________________________________________________________________
class AliITSClusterParam : public TObject {
 public:
  static AliITSClusterParam* Instance();
  AliITSClusterParam(){}
  virtual           ~AliITSClusterParam(){;}
  virtual void	Print(Option_t* option = "") const;
  void SetInstance(AliITSClusterParam *param){fgInstance = param;}
  static void GetNTeor(Int_t layer,const AliITSRecPoint* cl,
		Float_t theta,Float_t phi,
		Float_t &ny,Float_t &nz);
  static Int_t GetError(Int_t layer,const AliITSRecPoint*cl,
		 Float_t theta,Float_t phi,Float_t expQ,
		 Float_t &erry,Float_t &errz);

  //void FitData(TTree * tree);
  //
 protected:
  static AliITSClusterParam*   fgInstance; //! Instance of this class (singleton implementation)
  ClassDef(AliITSClusterParam,1)    //  ITS cluster parametrization class
};

#endif
