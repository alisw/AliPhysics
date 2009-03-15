// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Boris Polichtchouk & Per Thomas Hille for the ALICE           *
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

#include "AliHLTPHOSRcuAltroPatternTest.h"
#include <iostream>
#include "AliHLTPHOSPattern.h"
#include "AliHLTPHOSConstants.h"

using  namespace PhosHLTConst;
using  namespace std;

//  fPatternTestPtr = new  AliHLTPHOSRcuAltroPatternTest(fModuleID, fRcuX, fRcuZ, tmpPattern, ALTRO_MAX_SAMPLES );
AliHLTPHOSRcuAltroPatternTest::AliHLTPHOSRcuAltroPatternTest(const AliHLTUInt8_t moduleID, const AliHLTUInt8_t rcuX, 
							     const AliHLTUInt8_t rcuZ, const int *pattern, const  int length): AliHLTPHOSBase(),
															       fModuleID(moduleID), 
															       fRcuX(rcuX), 
															       fRcuZ(rcuZ),
															       fReferenceAltroPattern(0),
															       fCnt(0)
{
  fReferenceAltroPattern = new  AliHLTPHOSPattern(pattern, length);

  for(int z=0; z < NZROWSRCU; z++)
    {
      for(int x=0; x < NXCOLUMNSRCU; x++)
	{
	  for(int gain = 0; gain < NGAINS; gain ++)
	    {
	      fNEqual[z][x][gain] = 0;
	      fNNotEqual[z][x][gain] = 0;
	      fPerChannelPatterns[z][x][gain] = 0;
	    }
	}  
    } 
}



 AliHLTPHOSRcuAltroPatternTest::~  AliHLTPHOSRcuAltroPatternTest()
{
  //Destructor
}


int 
AliHLTPHOSRcuAltroPatternTest::ValidateAltroPattern(const int *inputPattern, const int samples, const int presamples) const
{
  if(fReferenceAltroPattern !=0)
    {
      return fReferenceAltroPattern->ValidatePattern(inputPattern,  samples,  presamples); 
    }
  else
    {
      return -99;
    }
}


int 
AliHLTPHOSRcuAltroPatternTest::AddPattern(const int *inputPattern,  const int z, const int x, const int gain, const int nSamples, const int nPresamples)
{
 
  if(fPerChannelPatterns[z][x][gain] == 0)
    {
      //    cout << "AliHLTPHOSRcuAltroPatternTest::AddPattern creating new pattern z = "<< z <<" x =" << x <<"  gain = "<< gain <<endl;
      fPerChannelPatterns[z][x][gain] = new AliHLTPHOSPattern(&inputPattern[nPresamples], nSamples);
      return 1;
    }
  else
    {
      //    cout << "AliHLTPHOSRcuAltroPatternTest::AddPattern adding new pattern to z = "<< z <<" x =" << x <<"  gain = "<< gain <<endl;
      fPerChannelPatterns[z][x][gain]->AddPattern(inputPattern,  nSamples, nPresamples);
      return 0;
    }

  
}

/*
 *Conts the number of linkes (patterns) in the linked list)
 */  
int 
AliHLTPHOSRcuAltroPatternTest::countPatterns(const AliHLTPHOSPattern *pattern) const
{
  int tmp = 0; 

  const AliHLTPHOSPattern *tmpPattern = pattern;

  while(tmpPattern !=0)
    {
      tmp ++;
      //    cout <<"tmp =" << tmp <<endl;
      tmpPattern = tmpPattern->GetNextPtr(); 
    }
  return tmp;
}


//   const int  GetPattern(int *pattern,  const int maxlengths =  ALTROMAXSAMPLES) const;
//   const int GetPatternLength() const {return  fPatternLength;}; 

/* Counts the total number of differen patterns dtetected across all channels
 * for one RCU.
 * @param length the only the number of samples spcified by length will be compared 
 * starting from zero. The default value is the maximum number of altro samples (1008, presamples exluded)
 * tyically the electronics is configured to read out (playback) less than the maximum number of samples so length can
 * be either equal or smaller than the lenght of the pattern.
 * @return the number of different patterns detected across all channels
 */

//   const int  GetPattern(int *pattern,  const int maxlengths =  ALTROMAXSAMPLES) const;
//   const int GetPatternLength() const {return  fPatternLength;}; 
int 
AliHLTPHOSRcuAltroPatternTest::countAllPatterns(const int length, const bool /*doprintpattern*/)
{
  fCnt ++;

  int tmpPatternArray[ALTROMAXSAMPLES];
  int tmplength = fReferenceAltroPattern->GetPattern(tmpPatternArray, length);
  AliHLTPHOSPattern *tmpPattern =  new AliHLTPHOSPattern(tmpPatternArray, tmplength); 
  
  for(int z=0; z < NZROWSRCU; z++)
    {
      for(int x=0; x < NXCOLUMNSRCU; x++)
	{
	  for(int gain = 0; gain < NGAINS; gain ++)
	    {
	      if(fPerChannelPatterns[z][x][gain] != 0)
		{
		  const AliHLTPHOSPattern *tmpLinkedPattern = 0; 
		  fPerChannelPatterns[z][x][gain]->GetPattern(tmpPatternArray);
		  tmpPattern->AddPattern(tmpPatternArray, tmplength);
		  tmpLinkedPattern = fPerChannelPatterns[z][x][gain]->GetNextPtr(); 
		  
		  while(tmpLinkedPattern !=0)
		    {
		      tmpLinkedPattern->GetPattern(tmpPatternArray);
		      tmpPattern->AddPattern(tmpPatternArray, tmplength);
		      tmpLinkedPattern = tmpLinkedPattern->GetNextPtr();    
		    }
		}
	    }
	}  
    } 
  
  int tmpcnt = countPatterns(tmpPattern);
  cout << "AliHLTPHOSRcuAltroPatternTest::countAllPatterns the total number of patterns found is" << tmpcnt << endl;
  PrintPatterns(tmpPattern);
  
  delete tmpPattern;
  
  return tmpcnt;

}


void                          
AliHLTPHOSRcuAltroPatternTest::PrintStatistics() const
{
  for(int z=0; z < NZROWSRCU; z++)
    {
      for(int x=0; x < NXCOLUMNSRCU; x++)
	{
	  for(int gain = 0; gain < NGAINS; gain ++)
	    {
	      int tmp = countPatterns(fPerChannelPatterns[z][x][gain]);
	      if(tmp > 2)
		{
		  cout << "z = "<< z <<" x " << x << " gain "<< gain << "  Has   "<< tmp  <<" patterns " <<endl;
		}
	    }
	}  
    } 

}


//TODO, both argument and fucntion should be constant, but coannot right
//because of some problem with the template Dumpdata function of AliHLTPHOSBase
void 
AliHLTPHOSRcuAltroPatternTest::PrintPatterns(AliHLTPHOSPattern *pattern) 
{
  //  const AliHLTPHOSPattern *tmpPattern = pattern;
  const AliHLTPHOSPattern *tmpPattern = pattern;

  cout << endl;
  cout <<"***********************************************************" << endl;
  cout <<"***************PRINTING PATTERNS **************************" << endl;
  cout <<"***********************************************************" << endl;

  while(tmpPattern != 0)
    {
      const_cast<AliHLTPHOSPattern *>(tmpPattern)->PrintPattern();
      tmpPattern = tmpPattern->GetNextPtr();
   }

  cout <<endl;

}

