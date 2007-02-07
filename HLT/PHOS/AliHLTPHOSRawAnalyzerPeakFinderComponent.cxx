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

//ClassImp(AliHLTPHOSRawAnalyzerPeakFinderComponent) 
AliHLTPHOSRawAnalyzerPeakFinderComponent gAliHLTPHOSRawAnalyzerPeakFinderComponent;

AliHLTPHOSRawAnalyzerPeakFinderComponent::AliHLTPHOSRawAnalyzerPeakFinderComponent():AliHLTPHOSRawAnalyzerComponent()
{
  /*
  Double_t tmpAVector[70];
  Double_t tmpTVector[70]; 
  analyzerPtr = new AliHLTPHOSRawAnalyzerPeakFinder();
  analyzerPtr->SetStartIndex(0);
  FILE *fp;
  fp = fopen("/home/perthi/cern/aliroot/AliRoot_head/HLT/PHOS/PFVectors/start0N70tau2fs10.txt", "r");
  
  if(fp != 0)
    {
      for(int i=0; i < 70; i++)
	{
	  fscanf(fp, "%lf", &tmpAVector[i]);
	}

      fscanf(fp, "\n");

      for(int i=0; i < 70; i++)
	{
	  	  fscanf(fp, "%lf", &tmpTVector[i]);
	}

      analyzerPtr->SetAVector(tmpAVector, 70);
      analyzerPtr->SetTVector(tmpTVector, 70);

      fclose(fp);

    }
  
  else
    {
      //   cout <<"AliHLTPHOSRawAnalyzerPeakFinderComponent, ERROR: could not  open PF vector file" << endl;
    }
  
  */
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
  //cout << "AliHLTPHOSRawAnalyzerPeakFinderComponent returning ID" << endl;
  return "PhosRawPeakFinder";
}


AliHLTComponent*
AliHLTPHOSRawAnalyzerPeakFinderComponent::Spawn()
{
  //  cout << "AliHLTPHOSRawAnalyzerPeakFinderComponent spawning new instance" << endl;
  return new AliHLTPHOSRawAnalyzerPeakFinderComponent;
}

