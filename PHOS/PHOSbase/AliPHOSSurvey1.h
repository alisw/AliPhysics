#ifndef ALIPHOSSURVEY1_H
#define ALIPHOSSURVEY1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1  2007/07/10 12:41:38  kharlov
 * Added a new class AliPHOSSurvet1 which read survey data from EDMS files
 *
 */

// A survey "reader" class, based on AliSurveyObj class.
// The source of input data is a text file in a standartized format
// downloaded from EDMS

#include "AliPHOSSurvey.h"

//  AliPHOSSurvey1 class is survey "reader" class
//  based on AliSurveyObj class.

class AliPHOSSurvey1 : public AliPHOSSurvey {
public:
  AliPHOSSurvey1(const TString &surveyFileName, const TString &namePrefix);
  virtual ~AliPHOSSurvey1();

private:
  enum EHardcoded {kNumberOfPoints = 452, kStartingPoint = 4};

  AliPHOSSurvey1(const AliPHOSSurvey1 &rhs);
  AliPHOSSurvey1 &operator = (const AliPHOSSurvey1 &rhs);

  ClassDef(AliPHOSSurvey1, 1) //Data reader, based on AliSurveyObj
};

#endif
