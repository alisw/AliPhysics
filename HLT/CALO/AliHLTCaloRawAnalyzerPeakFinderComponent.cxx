// $Id: AliHLTPHOSRawAnalyzerPeakFinderComponent.cxx 31490 2009-03-15 16:27:11Z odjuvsla $

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

#include "AliHLTCaloRawAnalyzerPeakFinderComponent.h"
#include "AliHLTCaloRawAnalyzerPeakFinder.h"
//#include <cstdlib>
//#include "AliHLTCaloCommonDefs.h"



//AliHLTCaloRawAnalyzerPeakFinderComponent gAliHLTCaloRawAnalyzerPeakFinderComponent;

//___________________________________________________________________________________________________________
AliHLTCaloRawAnalyzerPeakFinderComponent::AliHLTCaloRawAnalyzerPeakFinderComponent():AliHLTCaloRawAnalyzerComponentv3()
{
  fAnalyzerPtr = new AliHLTCaloRawAnalyzerPeakFinder();

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
AliHLTCaloRawAnalyzerPeakFinderComponent::~AliHLTCaloRawAnalyzerPeakFinderComponent()
{

  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
}


//___________________________________________________________________________________________________________
AliHLTCaloRawAnalyzerPeakFinderComponent::AliHLTCaloRawAnalyzerPeakFinderComponent(const AliHLTCaloRawAnalyzerPeakFinderComponent & ) : AliHLTCaloRawAnalyzerComponentv3()
{

}

//-----------------------------------------------------------------------------------------------------------
int
AliHLTCaloRawAnalyzerPeakFinderComponent::Deinit()
{
  
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTCaloRawAnalyzerCrudeComponent Deinit");
  return 0;
}

//___________________________________________________________________________________________________________


/*
const char* 
AliHLTCaloRawAnalyzerPeakFinderComponent::GetComponentID()
{
  return "PhosRawPeakFinder";
}
*/


 /*
//___________________________________________________________________________________________________________
Bool_t 
AliHLTCaloRawAnalyzerPeakFinderComponent::LoadPFVector()
{
  return LoadPFVector(PFDEFAULTSTARTINDEX,  PFDEFAULTNSAMPLES, DEFAULTTAU, DEFAULTFS );
}


//___________________________________________________________________________________________________________
Bool_t 
AliHLTCaloRawAnalyzerPeakFinderComponent::LoadPFVector(int startIndex, int nSamples, int tau, int fs)
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
 */


//___________________________________________________________________________________________________________

/*
AliHLTComponent*
AliHLTCaloRawAnalyzerPeakFinderComponent::Spawn()
{
  return new AliHLTCaloRawAnalyzerPeakFinderComponent;
}
*/
