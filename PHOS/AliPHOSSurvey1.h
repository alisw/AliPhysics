#ifndef ALIPHOSSURVEY1_H
#define ALIPHOSSURVEY1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

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
