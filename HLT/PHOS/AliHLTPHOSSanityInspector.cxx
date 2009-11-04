// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** 
 * @file   AliHLTPHOSSanityInspector.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Sanity inspector for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSSanityInspector.h"
#include "Rtypes.h"

ClassImp(AliHLTPHOSSanityInspector);


AliHLTPHOSSanityInspector::AliHLTPHOSSanityInspector() : 
  //  AliHLTPHOSBase(),
  fMaxDifference(120)
{
  //See header file for documentation
}


AliHLTPHOSSanityInspector::~AliHLTPHOSSanityInspector()
{
  //See header file for documentation
}



Int_t  
AliHLTPHOSSanityInspector::CheckInsanity(const UShort_t* data, const Int_t N) const
{
   //See header file for documentation

  for(Int_t i = 1; i < N; i++)
    {
      if((((Int_t)data[i] - (Int_t)data[i-1]) > fMaxDifference) || (((Int_t)data[i-1] - (Int_t)data[i]) > fMaxDifference))
	return 1;
    }
  return 0;
}

/*
Int_t
AliHLTPHOSSanityInspector::CheckInsanity(Int_t* data, Int_t N)
{
   //See header file for documentation
  for(Int_t i = 1; i < N; i++)
  {
    if((((Int_t)data[i] - (Int_t)data[i-1]) > fMaxDifference) || (((Int_t)data[i-1] - (Int_t)data[i]) > fMaxDifference))
      return 1;
  }
  return 0;
}
*/


Int_t 
AliHLTPHOSSanityInspector::CheckAndHealInsanity(UShort_t* data, Int_t N)
{
   //See header file for documentation

  Int_t crazyness = 0;

  if(N > 3)
    {
      if((((Short_t)data[0] - (Short_t)data[1]) > fMaxDifference) || (((Short_t)data[1] - (Short_t)data[0]) > fMaxDifference))
	return -1;
      if((((Short_t)data[1] - (Short_t)data[2]) > fMaxDifference) || (((Short_t)data[2] - (Short_t)data[1]) > fMaxDifference))
	return -1;


      for(Short_t i = 2; i < N - 3; i++)
	{
	  if((((Short_t)data[i] - (Short_t)data[i+1]) > fMaxDifference) || (((Short_t)data[i+1] - (Short_t)data[i]) > fMaxDifference))
	    {
	      i++;
	      if((((Short_t)data[i] -(Short_t)data[i+1]) > fMaxDifference) || (((Short_t)data[i+1] - (Short_t)data[i]) > fMaxDifference))
		{
		  i++;
		  if((((Short_t)data[i] - (Short_t)data[i+1]) > fMaxDifference) || (((Short_t)data[i+1] - (Short_t)data[i]) > fMaxDifference))
		    {
		      return -2;  //Too crazy
		    }
		  data[i-1] = ((Short_t)data[i] + (Short_t)data[i-2])/2;
		  crazyness++;
		}
	      else 
		return -3;    //Two spikes in a row? 
	    }
	}
      
      
      
      if((((Short_t)data[N - 3] -(Short_t) data[N - 2]) > fMaxDifference) || 
	 (((Short_t)data[N - 2] - (Short_t)data[N - 3]) > fMaxDifference))
	{
	  if((((Short_t)data[N - 2] - (Short_t)data[N - 1]) > fMaxDifference) || 
	     (((Short_t)data[N - 1] - (Short_t)data[N - 2]) > fMaxDifference))
	    {
	      data[N - 2] = ((Short_t)data[N - 3] +  (Short_t)data[N - 1])/2;
	      return crazyness++;
	    }
	  return -4;

	}
      
      if((((Short_t)data[N - 2] - (Short_t)data[N - 1]) > fMaxDifference) || 
	 (((Short_t)data[N - 1] - (Short_t)data[N - 2]) > fMaxDifference))
	{
	  //	  (Short_t)data[N - 3] = (Short_t)data[N - 4] -(Short_t) data[N - 5] + (Short_t)data[N-4];
	  data[N - 1] = data[N - 2];
	  return crazyness++;
	}
      
    }
  
  return crazyness;
  
}
