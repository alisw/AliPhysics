/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
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

#include "AliHLTPHOSRawAnalyzerPeakFinderComponent.h"
#include "AliHLTPHOSRawAnalyzerPeakFinder.h"
#include <cstdlib>
#include "AliHLTPHOSCommonDefs.h"


AliHLTPHOSRawAnalyzerPeakFinderComponent gAliHLTPHOSRawAnalyzerPeakFinderComponent;

AliHLTPHOSRawAnalyzerPeakFinderComponent::AliHLTPHOSRawAnalyzerPeakFinderComponent():AliHLTPHOSRawAnalyzerComponent()
{
  char tmpPFPath[PF_MAX_PATH_LENGTH];
  //  cout <<"ALICE_ROOT ="<<getenv("ALICE_ROOT") << endl;

  Double_t tmpAVector[PF_DEFAULT_N_SAMPLES];
  Double_t tmpTVector[PF_DEFAULT_N_SAMPLES]; 
  analyzerPtr = new AliHLTPHOSRawAnalyzerPeakFinder();
  analyzerPtr->SetStartIndex(PF_DEFAULT_STARTINDEX);

  sprintf(tmpPFPath,"%s%s/start%dN%dtau%dfs%d.txt", getenv("ALICE_ROOT"), PF_VECTOR_DIR, PF_DEFAULT_STARTINDEX,  PF_DEFAULT_N_SAMPLES, DEFAULT_TAU, DEFAULT_FS);

  cout <<"PF PATH =" << tmpPFPath << endl;

  FILE *fp;

  fp = fopen(tmpPFPath, "r");
  
  if(fp != 0)
    {
      for(int i=0; i <  PF_DEFAULT_N_SAMPLES; i++)
	{
	  fscanf(fp, "%lf", &tmpAVector[i]);
	}


      fscanf(fp, "\n");

      for(int i=0; i < PF_DEFAULT_N_SAMPLES; i++)
	{
	  	  fscanf(fp, "%lf", &tmpTVector[i]);
	}

      analyzerPtr->SetAVector(tmpAVector,  PF_DEFAULT_N_SAMPLES);
      analyzerPtr->SetTVector(tmpTVector,  PF_DEFAULT_N_SAMPLES);

      fclose(fp);

    }
  
  else
    {
      HLTFatal("ERROR: could not  open PF vector file");
    }
  
  
} 

AliHLTPHOSRawAnalyzerPeakFinderComponent::~AliHLTPHOSRawAnalyzerPeakFinderComponent()
{

}



AliHLTPHOSRawAnalyzerPeakFinderComponent::AliHLTPHOSRawAnalyzerPeakFinderComponent(const AliHLTPHOSRawAnalyzerPeakFinderComponent & ) : AliHLTPHOSRawAnalyzerComponent()
{

}


const char* 
AliHLTPHOSRawAnalyzerPeakFinderComponent::GetComponentID()
{
  return "PhosRawPeakFinder";
}


AliHLTComponent*
AliHLTPHOSRawAnalyzerPeakFinderComponent::Spawn()
{
  return new AliHLTPHOSRawAnalyzerPeakFinderComponent;
}

