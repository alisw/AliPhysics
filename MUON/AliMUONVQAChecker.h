#ifndef ALIMUONVQACHECKER_H
#define ALIMUONVQACHECKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVQAChecker
/// \brief Base class for a MUON QA checker
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TObjArray;
class AliMUONRecoParam;

class AliMUONVQAChecker : public TObject
{
public:
  /// Classification of errors severity
  enum ECheckCode {
    kFatal=-1,  ///< error is really serious
    kError=0,   ///< normal error, i.e. something is wrong
    kWarning=1, ///< warning, i.e. might become an error later on
    kInfo=2     ///< just so you know...
  };
  
  enum EColor {
    kInfoColor=kSpring-8, ///< color for information (online convention)
    kWarningColor=kOrange, ///< color for warning (online convention)
    kErrorColor=kRed, ///< color for normal error (online convention)
    kFatalColor=kMagenta+1 ///< color for fatal error (online convention)
  };
  
  AliMUONVQAChecker();
  virtual ~AliMUONVQAChecker();
  
  /// Check the QA object(s) for the raw data
  virtual ECheckCode * CheckRaws(TObjArray** list, const AliMUONRecoParam* recoParam) = 0;
  
  /// Check the QA object(s) for the RecPoints
  virtual ECheckCode * CheckRecPoints(TObjArray** list, const AliMUONRecoParam* recoParam) = 0;
  
  /// Check the QA object(s) for the ESD
  virtual ECheckCode * CheckESD(TObjArray** list, const AliMUONRecoParam* recoParam) = 0;
  
  ClassDef(AliMUONVQAChecker,1) // Interface for a MUON QA checker
};

#endif
