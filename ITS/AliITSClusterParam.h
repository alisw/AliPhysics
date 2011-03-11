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
//#include "AliITSRecPoint.h"

class AliITSRecPoint;

//_____________________________________________________________________________
class AliITSClusterParam : public TObject {
 public:
  static AliITSClusterParam* Instance();
  AliITSClusterParam(){}
  virtual           ~AliITSClusterParam(){;}
  virtual void	Print(Option_t* option = "") const;
  void SetInstance(AliITSClusterParam *param){fgInstance = param;}
  static void GetNTeor(Int_t layer,const AliITSRecPoint* cl,
		       Float_t tgl,Float_t tgphitr,
		       Float_t &ny,Float_t &nz);
  static Int_t GetError(Int_t layer,const AliITSRecPoint*cl,
			Float_t tgl,Float_t tgphitr,Float_t expQ,
			Float_t &erry,Float_t &errz,Float_t &covyz,
			Bool_t addMisalErr=kTRUE);
  static Int_t GetError(Int_t layer,const AliITSRecPoint*cl,
			Float_t tgl,Float_t tgphitr,Float_t expQ,
			Float_t &erry,Float_t &errz,
			Bool_t addMisalErr=kTRUE) {
                                  Float_t covyz;
    return GetError(layer,cl,tgl,tgphitr,expQ,erry,errz,covyz,addMisalErr);
  }

  //void FitData(TTree * tree);
  //
 protected:
  static AliITSClusterParam*   fgInstance; //! Instance of this class (singleton implementation)
  static Int_t GetErrorOrigRecPoint(const AliITSRecPoint*cl,
				    Float_t &erry,Float_t &errz,Float_t &covyz);
  static Int_t GetErrorParamMI(Int_t layer,const AliITSRecPoint*cl,
			       Float_t tgl,Float_t tgphitr,Float_t expQ,
			       Float_t &erry,Float_t &errz);
  static Int_t GetErrorParamAngle(Int_t layer,const AliITSRecPoint*cl,
				  Float_t tgl,Float_t tgphitr,
				  Float_t &erry,Float_t &errz);
  static Int_t GetErrorParamAngleOld(Int_t layer,const AliITSRecPoint*cl,
                                  Float_t tgl,Float_t tgphitr,
                                  Float_t &erry,Float_t &errz);

  ClassDef(AliITSClusterParam,1)    //  ITS cluster parametrization class
};

#endif
