#ifndef ALIMUONVQACHECKER_H
#define ALIMUONVQACHECKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVQAChecker
/// \brief Base class for a MUON QA checker
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TObjArray;
class AliMUONRecoParam;

class AliMUONVQAChecker : public TObject
{
public:
  enum ECheckCode {
    kFatal=-1,
    kError=0,
    kWarning=1,
    kInfo=2
  };
  
  AliMUONVQAChecker();
  virtual ~AliMUONVQAChecker();
  
  virtual ECheckCode * CheckRaws(TObjArray** list, AliMUONRecoParam* recoParam) = 0;
  virtual ECheckCode * CheckRecPoints(TObjArray** list, AliMUONRecoParam* recoParam) = 0;
  virtual ECheckCode * CheckESD(TObjArray** list, AliMUONRecoParam* recoParam) = 0;
  
  ClassDef(AliMUONVQAChecker,1) // Interface for a MUON QA checker
};

#endif
