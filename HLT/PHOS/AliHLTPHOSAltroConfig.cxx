// $Id$

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTPHOSAltroConfig.h"
#include <stdio.h>
//#include <unistd.h>
#ifndef __SUNPRO_CC
#include <cstdlib>
#else
#include <stdlib.h>
#endif


AliHLTPHOSAltroConfig::AliHLTPHOSAltroConfig() : fNPresamples(900), 
						 fNSamples(15),
						 fNTotalSamples(fNPresamples + fNSamples),
						 fIsAltroZeroSupressed(false),
						 fIsAltroBaselineSubtraction(false)


{

}


AliHLTPHOSAltroConfig::~AliHLTPHOSAltroConfig()
{
  //comment
}

void 
AliHLTPHOSAltroConfig:: PrintAltroDefaultValues() const
{
  //comment
  printf("\n AliHLTPHOSAltroConfig Default  Values\n");
  printf("Presamples = %d\n", fNPresamples);
  printf("NSamples = %d\n", fNSamples);

  if(fIsAltroZeroSupressed == true)
    {
      printf("fIsAltroZeroSupressed = true\n");
    }
  else
    {
      printf("fIsAltroZeroSupressed = false\n");
    }

  
  if(fIsAltroBaselineSubtraction == true)
    {
      printf("fIsAltroBaselineSubtraction = true\n");
    }
  else
    {
      printf("fIsAltroBaselineSubtraction = false\n");
    }
  //    `fIsSoftwareBaselinesubtraction 
}


void 
AliHLTPHOSAltroConfig::SetNPresSamples(const int presamples)
{
  //comment
  fNPresamples =  presamples;
}

void 
AliHLTPHOSAltroConfig::SetNSamples(const int samples)
{
  //comment
  fNSamples = samples;
}

void 
AliHLTPHOSAltroConfig::SetAltroZeroSupression(const bool isZeroSupressed)
{
  //comment
  fIsAltroZeroSupressed = isZeroSupressed;
}


void 
AliHLTPHOSAltroConfig::SetAltroBaselineSubtraction(const bool isAltroBaselineSubtraction)
{
  //comment
  fIsAltroBaselineSubtraction = isAltroBaselineSubtraction;
}


