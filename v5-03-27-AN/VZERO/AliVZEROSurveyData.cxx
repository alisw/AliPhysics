/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliVZEROSurveyData.cxx,                                            */

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// class for VZERO survey points management                                //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliVZEROSurveyData.h"

ClassImp(AliVZEROSurveyData)

//________________________________________________________________
AliVZEROSurveyData::AliVZEROSurveyData()
{
  
}

//________________________________________________________________
void AliVZEROSurveyData::Reset()
{
  
}

//________________________________________________________________
AliVZEROSurveyData::AliVZEROSurveyData(const char* name)
{
  TString namst = "Survey_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

}

//________________________________________________________________
AliVZEROSurveyData::AliVZEROSurveyData(const AliVZEROSurveyData& surveyda) :
  TNamed(surveyda)
{
// copy constructor

  SetName(surveyda.GetName());
  SetTitle(surveyda.GetName());
  
  for(int t=0; t<3; t++) { 
      fngA[t] = surveyda.GetPointA(t);
      fngB[t] = surveyda.GetPointB(t);
      fngC[t] = surveyda.GetPointC(t);
      fngD[t] = surveyda.GetPointD(t);
  }   
  
}

//________________________________________________________________
AliVZEROSurveyData &AliVZEROSurveyData::operator =(const AliVZEROSurveyData& surveyda)
{
// assignment operator

  SetName(surveyda.GetName());
  SetTitle(surveyda.GetName());
  
  for(int t=0; t<3; t++) { 
      fngA[t] = surveyda.GetPointA(t);
      fngB[t] = surveyda.GetPointB(t);
      fngC[t] = surveyda.GetPointC(t);
      fngD[t] = surveyda.GetPointD(t);
  }   
    
  return *this;
  
}

//________________________________________________________________
AliVZEROSurveyData::~AliVZEROSurveyData()
{
  
}                                                                                  

//________________________________________________________________
void AliVZEROSurveyData::SetPointA(Float_t* ngA)
{
  if(ngA) for(int t=0; t<3; t++) fngA[t] = ngA[t];
  else for(int t=0; t<3; t++) fngA[t] = 0.0;
}

//________________________________________________________________
void AliVZEROSurveyData::SetPointB(Float_t* ngB)
{
  if(ngB) for(int t=0; t<3; t++) fngB[t] = ngB[t];
  else for(int t=0; t<3; t++) fngB[t] = 0.0;
}
//________________________________________________________________
void AliVZEROSurveyData::SetPointC(Float_t* ngC)
{
  if(ngC) for(int t=0; t<3; t++) fngC[t] = ngC[t];
  else for(int t=0; t<3; t++) fngC[t] = 0.0;
}

//________________________________________________________________
void AliVZEROSurveyData::SetPointD(Float_t* ngD)
{
  if(ngD) for(int t=0; t<3; t++) fngD[t] = ngD[t];
  else for(int t=0; t<3; t++) fngD[t] = 0.0;
}

