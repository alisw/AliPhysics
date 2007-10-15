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


AliHLTPHOSAltroConfig::AliHLTPHOSAltroConfig() : fNPresamples(55), 
						 fNSamples(15), 
						 fIsAltroZeroSupressed(false),
						 fIsAltroBaselineSubtraction(false)
{
  fNTotalSamples =  fNPresamples + fNSamples ; 

  FILE *fp = fopen("hltAltroConfig.txt", "r");
  if(fp == 0)
    {
     // printf("\nAliHLTPHOSConfig::AliHLTPHOSConfig(): WARNING: Could not fine file \"hltConfig.txt\" \n");
     // printf("Default values will be used\n");
     // PrintAltroDefaultValues();
    }
  else
    {
      printf("Reading PHOS HLT configurations from file is no yest implemented\n");
      printf("You can use Setter functions of  AliHLTPHOSConfig to set the appropriate parameters\n");
      printf("See  AliHLTPHOSConfig.h for details\n");
      printf("Using default values for the moment\n");
      PrintAltroDefaultValues(); 
    }

}


AliHLTPHOSAltroConfig::~AliHLTPHOSAltroConfig()
{

}

void 
AliHLTPHOSAltroConfig:: PrintAltroDefaultValues()
{
  printf("\n AliHLTPHOSAltroConfig Default  Values\n");
  printf("Presamples = \n", fNPresamples);
  printf("NSamples = \n", fNSamples);

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
AliHLTPHOSAltroConfig::SetNPresSamples(int presamples)
{
  fNPresamples =  presamples;
}

void 
AliHLTPHOSAltroConfig::SetNSamples(int samples)
{
  fNSamples = samples;
}

void 
AliHLTPHOSAltroConfig::SetAltroZeroSupression(bool isZeroSupressed)
{
  fIsAltroZeroSupressed = isZeroSupressed;
}

void 
AliHLTPHOSAltroConfig::SetAltroBaselineSubtraction(bool isAltroBaselineSubtraction)
{
  fIsAltroBaselineSubtraction = isAltroBaselineSubtraction;
}


