#ifndef ALIRECOPARAM_H
#define ALIRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Base Class for Detector reconstruction parameters                         //
// Revision: cvetan.cheshkov@cern.ch 12/06/2008                              //
// Its structure has been revised and it is interfaced to AliEventInfo.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TNamed.h"
class AliDetectorRecoParam;
class AliEventInfo;

class AliRecoParam : public TNamed
{

 public: 
  AliRecoParam();
  AliRecoParam(const char *detector);
  virtual ~AliRecoParam();  
  //
  virtual void                  Print(Option_t *option="") const;
  TObjArray                    *GetAllRecoParams() const { return fRecoParamArray; }
  virtual AliDetectorRecoParam *GetRecoParam(const AliEventInfo &/*evInfo*/) const;
  void                          AddRecoParam(AliDetectorRecoParam* param);

protected:

  TObjArray *fRecoParamArray;   //array with reconstruction-parameter objects

private:

  AliRecoParam(const AliRecoParam&); // Not implemented
  AliRecoParam& operator=(const AliRecoParam&); // Not implemented

  ClassDef(AliRecoParam, 2)
};


#endif
