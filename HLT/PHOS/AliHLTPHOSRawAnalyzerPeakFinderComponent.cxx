// $Id$

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
//#include <cstdlib>
//#include "AliHLTPHOSCommonDefs.h"


AliHLTPHOSRawAnalyzerPeakFinderComponent gAliHLTPHOSRawAnalyzerPeakFinderComponent;

//___________________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerPeakFinderComponent::AliHLTPHOSRawAnalyzerPeakFinderComponent():AliHLTPHOSRawAnalyzerComponent()
{
  fAnalyzerPtr = new AliHLTPHOSRawAnalyzerPeakFinder();

  if(LoadPFVector() == kFALSE)
    {
      //      cout << "Warning, could not load PF vectors" << endl;
    }
  else 
    {
      //    cout << "Loaded PF vectors" << endl;
    }
} 


//___________________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerPeakFinderComponent::~AliHLTPHOSRawAnalyzerPeakFinderComponent()
{

  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
}


//___________________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerPeakFinderComponent::AliHLTPHOSRawAnalyzerPeakFinderComponent(const AliHLTPHOSRawAnalyzerPeakFinderComponent & ) : AliHLTPHOSRawAnalyzerComponent()
{

}

//-----------------------------------------------------------------------------------------------------------
int
AliHLTPHOSRawAnalyzerPeakFinderComponent::Deinit()
{
  
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerCrudeComponent Deinit");
  return 0;
}

//___________________________________________________________________________________________________________
const char* 
AliHLTPHOSRawAnalyzerPeakFinderComponent::GetComponentID()
{
  return "PhosRawPeakFinder";
}

//___________________________________________________________________________________________________________
Bool_t 
AliHLTPHOSRawAnalyzerPeakFinderComponent::LoadPFVector()
{
  return LoadPFVector(PFDEFAULTSTARTINDEX,  PFDEFAULTNSAMPLES, DEFAULTTAU, DEFAULTFS );
}


//___________________________________________________________________________________________________________
Bool_t 
AliHLTPHOSRawAnalyzerPeakFinderComponent::LoadPFVector(int startIndex, int nSamples, int tau, int fs)
{
  char tmpPFPath[PFMAXPATHLENGTH];
  Double_t * tmpAVector = new Double_t[nSamples];
  Double_t * tmpTVector = new Double_t[nSamples]; 
  sprintf(tmpPFPath,"%s%s/start%dN%dtau%dfs%d.txt", getenv("ALICE_ROOT"), PFVECTORDIR, startIndex, nSamples, tau, fs);
  FILE *fp;
  fp = fopen(tmpPFPath, "r");
  
  Int_t res = 0; //OD to get rid of warnings
  if(fp != 0)
    {
      for(int i=0; i <  nSamples; i++)
	{
	  res = fscanf(fp, "%lf", &tmpAVector[i]);
	}

      res = fscanf(fp, "\n");

      for(int i=0; i < nSamples; i++)
	{
	  res = fscanf(fp, "%lf", &tmpTVector[i]);
	}
      fAnalyzerPtr->SetAVector(tmpAVector,  nSamples);
      fAnalyzerPtr->SetTVector(tmpTVector,  nSamples);
      fclose(fp);
      delete [] tmpAVector;
      delete [] tmpTVector;
      return kTRUE;
    }
  
  else
    {
      delete [] tmpAVector;
      delete [] tmpTVector;
      HLTFatal("ERROR: could not  open PF vector file");
      return kFALSE;
    }
}


//___________________________________________________________________________________________________________
AliHLTComponent*
AliHLTPHOSRawAnalyzerPeakFinderComponent::Spawn()
{
  return new AliHLTPHOSRawAnalyzerPeakFinderComponent;
}

