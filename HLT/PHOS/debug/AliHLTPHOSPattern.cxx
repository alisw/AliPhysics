// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Per Thomas Hille for the ALICE                                *
 * offline/HLT Project. Contributors are mentioned in the code where      *
 * appropriate.                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSPattern.h"
#include <iostream>
#include "AliHLTPHOSUtilities.h" 

using namespace std;


AliHLTPHOSPattern::AliHLTPHOSPattern(const int *pattern, const int length) : AliHLTPHOSBase(),
									     fPatternLength(0), 
									     fN(0), 
									     fCnt(0), 
									     fPattern(0),
									     fUtilitiesPtr(0)
{
  fUtilitiesPtr = new  AliHLTPHOSUtilities(); 
  SetPattern(pattern, length);
}



AliHLTPHOSPattern::~AliHLTPHOSPattern()
{
  //  cout << "AliHLTPHOSPattern::~AliHLTPHOSPattern()"  << endl;
  delete  fPattern;
  fPattern = 0;
}



/*
 * Checks if the input pattern is equal to fPattern
 * @param inputPattern the pattern to be compared to @ref fPattern
 * @param length the number of samples to compare i.e starting at zero and comparing "length samples"
 * @retrun the number of samples that mismatch, zero fo a complete macth, or -1 if the comparison
 * could not be done due to an incorrect input (wrong number of samples and/or presamples)
*/
int 
AliHLTPHOSPattern::ValidatePattern(const int *readbackpattern, const int nSamples, const int nPresamples) const
{
  int iRet = 0;

  if(CheckPatternLength(nSamples, nPresamples) == true)
    {
      // The presamples must be exluded from the comparison since they are not present in the altro pattern memeory
      iRet =  DoComparePattern(&readbackpattern[nPresamples], fVal, nSamples);
    }
  else
    {
      iRet = -1;
    }
  
  return iRet;
}


int 
AliHLTPHOSPattern::AddPattern(const int *readbackpattern,  const int nSamples, const int nPresamples)
{
  int iRet = ValidatePattern(readbackpattern, nSamples, nPresamples);
  if(iRet == 0)
    {
      fN ++;
    }
  else if(iRet > 0)
    {
      if(fPattern != 0)
	{
	  //	  fPattern->ValidatePattern(readbackpattern, nSamples, nPresamples);
	  fPattern->AddPattern(readbackpattern, nSamples, nPresamples);
	}
      else
	{
	  //	  cout << "AliHLTPHOSPattern::AddPattern, creating new pattern"  << endl;
	  fPattern = new AliHLTPHOSPattern(readbackpattern, nSamples);
	}
    }
  return iRet;
}


bool 
AliHLTPHOSPattern::CheckDoExistPattern(const int *readbackpattern,  const int nSamples, const int nPresamples)
{
  bool iRet = false;
  
  if(ValidatePattern(readbackpattern, nSamples, nPresamples) == 0)
    {
      iRet = true;
    }
  else if(fPattern !=0)
    {
      iRet = fPattern->CheckDoExistPattern(readbackpattern, nSamples, nPresamples);
    }
      
  return iRet;
}


int 
AliHLTPHOSPattern::DoComparePattern(const int *pattern1, const int *pattern2, const int length) const
{
  int nomatch = 0;
  
  for(int i=0; i< length; i++)
    {
      if(pattern1[i] != pattern2[i])
	{
	  nomatch ++;
	}
    }
  return nomatch;
}


bool 
AliHLTPHOSPattern::CheckPatternLength(const int nSamples, const int nPresamples) const
{
  bool iRet = true;

  if( nSamples  > ALTROMAXSAMPLES)
    {
      cout << "Warning: attemp to set pattern array of length " << nSamples << " wich is out of range" <<endl;
      cout <<"Max length of pattern array is " << ALTROMAXSAMPLES  <<endl;
      iRet = false;
    }
  else if(nPresamples >  ALTROMAXPRESAMPLES)
    {
      cout << "ERROR: attemp to set teh number of  " << nPresamples << " wich is out of range" <<endl;
      cout <<"Max length of pattern array is " << ALTROMAXPRESAMPLES  <<endl;
      iRet = false;
    }

  return iRet;
}


void 
AliHLTPHOSPattern::SetPattern(const int *pattern, const int length)
{
  for(int i=0; i< ALTROMAXSAMPLES; i++)
     {
       fVal[i] = 0;
     }

  for(int i = 0; i < length; i++)
     {
       fVal[i] = pattern[i];
     }

  if(length > ALTROMAXSAMPLES)
    {
      fPatternLength = ALTROMAXSAMPLES;
    }
  else
    {
      fPatternLength = length;
    }
}


/*
 * Returns the pattern integer array
 * @param *pattern the pattern will be filleded in the array pointed to by this pointer
 * @param maxlength the maximum number of samples to fill 
 * @return the number of samples actually filled, retruns the smalles value of maxlength and
 * and the length of the pattern array
 */
int  
AliHLTPHOSPattern::GetPattern(int *pattern,  const int maxlength) const
{
  int tmplength = 0;

  if(maxlength > fPatternLength )
    {
      tmplength = fPatternLength;
    }
  else
    {
      tmplength = maxlength;
    }

  for(int i=0; i < tmplength; i++)
    {
      pattern[i] = fVal[i];
    }
  return  tmplength;
}



void 
AliHLTPHOSPattern::PrintPattern(const int nPerLine)
{
  fUtilitiesPtr->DumpData(fVal, fPatternLength  , nPerLine);
}

// void  DumpData(T *array, int N, int nPerLine)
