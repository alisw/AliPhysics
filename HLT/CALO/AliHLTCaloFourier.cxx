// $Id: AliHLTPHOSFourier.cxx 34951 2009-09-23 14:35:38Z phille $

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
#include "AliHLTCaloFourier.h"

//#include "AliHLTCaloRcuFFTDataStruct.h"

#include  "AliHLTCaloConstants.h"

ClassImp(AliHLTCaloFourier);



AliHLTCaloFourier::AliHLTCaloFourier() :  fFFT_own(0),
					  fFFTInputArray(0),
					  fFFTOutputArray(0),
					  fIsFirstChannel(true),
					  fFixedDataSize(0), 
					  fFFTOupuStruct(),
					  fCurrentEvent(0)
{
  

}




AliHLTCaloFourier::AliHLTCaloFourier(const AliHLTCaloFourier&) : fFFT_own(0),
								 fFFTInputArray(0),
								 fFFTOutputArray(0),
								 fIsFirstChannel(true),
								 fFixedDataSize(0),
								 fFFTOupuStruct(), 
								 fCurrentEvent(0)

{
  
}



AliHLTCaloFourier::~AliHLTCaloFourier()
{

}


 
AliHLTCaloRcuFFTDataStruct 
AliHLTCaloFourier::GetPSD()
{
  return fFFTOupuStruct;
}



void 
AliHLTCaloFourier::ProcessFourier(const Int_t *data, const int length, const int /*z*/, const int /*x*/, const int gain, const int event)
{
  Double_t  re = 0;
  Double_t  im = 0;

  if( (event > 0 ) && (event != fCurrentEvent ))
    {
      fCurrentEvent = event;
      ResetEventPSD(gain);
    }

  if(fIsFirstChannel == true)
    {
      fCurrentEvent = event;
      fIsFirstChannel = false;
      fFixedDataSize = length;
      Int_t n_size = fFixedDataSize +1;
      fFFT_own = TVirtualFFT::FFT(1, &n_size, "R2C ES K"); 
      Init();
  }

  if( CheckSignal(data, length) == true)
    {
      Int2Double(data, fFFTInputArray,  fFixedDataSize );
      fFFT_own->SetPoints( fFFTInputArray );	 
      fFFT_own->Transform();

      for(int j=0; j < length; j++)
	{
	  fFFT_own->GetPointComplex(j,  re,  im);
	  fFFTOupuStruct.fGlobalAccumulatedPSD[gain][j] +=  EvaluateMagnitude(re, im);
	  fFFTOupuStruct.fGlobalLastPSD[gain][j] +=  EvaluateMagnitude(re, im);

	}
    }
  //  printf("AliHLTCaloFourier::ProcessFourier;  (z, x, gain)  =  (%d, %d, %d),  length = %d",  z,  x,  gain, length);
}

void 
AliHLTCaloFourier::ResetEventPSD(const int gain)
{
  cout << " AliHLTCaloFourier::ResetEventPS, resetting event PSD "<< endl;
  for(int i = 0;  i < fFixedDataSize; i++ )
    {
      fFFTOupuStruct.fGlobalLastPSD[gain][i] = 0;
    }
}


bool 
//AliHLTCaloFourier::CheckSignal(const UInt_t *data, const int length)
AliHLTCaloFourier::CheckSignal(const Int_t *data, const int length)
{
  //  UInt_t tmpMax =  Max(  const_cast< UInt_t *>(data), length);
  //  UInt_t tmpMin =  Min(  const_cast< UInt_t *>(data), length);


  Int_t tmpMax =  Max(  const_cast< Int_t *>(data), length);
  Int_t tmpMin =  Min(  const_cast< Int_t *>(data), length);


  // if( (tmpMax -tmpMin) > 200)
  if( (tmpMax -tmpMin) > 100)
   {
      cout << "FourierAna::CheckSignal min = "<< tmpMin << "  max =  " << tmpMax << endl;
    }
  
  if( (tmpMax >= AliHLTCaloConstants::GetMAXBINVALUE() ) || tmpMin < 1 )
    {
      cout << "ERROR, FourierAna::CheckSignal failed, signal out of range, min= "<< tmpMin << "max = " << tmpMax << endl;
      return false;
    }
  else
    {
      return true;
    }
  
  //return true;
}


void 
AliHLTCaloFourier::Init()
{
  fFFTInputArray = new double[fFixedDataSize];
  fFFTOutputArray = new double[fFixedDataSize];
  
  for(int gain = 0; gain <  AliHLTCaloConstants::GetNGAINS(); gain ++)
    {
      fFFTOupuStruct.fDataLength = fFixedDataSize;
 
      for(int k= 0; k <fFixedDataSize; k++ )
	{
	  fFFTInputArray[k] = 0;
	  fFFTOutputArray[k] = 0;

	}
      for(int i=0; i <  AliHLTCaloConstants::GetALTROMAXSAMPLES()  ; i++)
	{
	  fFFTOupuStruct.fGlobalAccumulatedPSD[gain][i] = 0;
	  fFFTOupuStruct.fGlobalLastPSD[gain][i] = 0;
	}
    }
}


double 
AliHLTCaloFourier::EvaluateMagnitude(const double re, const double im)
{
  return re*re + im*im;
}


void 
AliHLTCaloFourier::Int2Double(const Int_t *inputarray, double *outputarray, const int size)
{
  for(int i=0; i< size; i++)
    {
      outputarray[i] = (double)inputarray[i];
    }
}
