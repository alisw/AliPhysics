/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthomas.hille@yale.edu>                    *
 * for the ALICE HLT Project.                                             * 
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

// Evaluation of peak position
// and amplitude using Neural Networks (NN)
// ------------------
// ------------------
// ------------------


#include "AliCaloRawAnalyzerNN.h"
#include "AliCaloNeuralFit.h"
#include "AliCaloFitResults.h"
#include "AliCaloBunchInfo.h"

#include <iostream>

using namespace std;

ClassImp( AliCaloRawAnalyzerNN )

AliCaloRawAnalyzerNN::AliCaloRawAnalyzerNN() : AliCaloRawAnalyzer("Neural Network", "NN"), fNeuralNet(0)
{
  // Comment

  fNeuralNet = new AliCaloNeuralFit();

  for(int i=0; i < 5 ; i++)
    {
      fNNInput[i]  = 0;
    }

}


AliCaloRawAnalyzerNN::~AliCaloRawAnalyzerNN()
{
  delete fNeuralNet;
}


AliCaloFitResults 
AliCaloRawAnalyzerNN::Evaluate( const vector<AliCaloBunchInfo> &bunchvector, 
				       const UInt_t altrocfg1,  const UInt_t altrocfg2 )
{
  // The eveluation of  Peak position and amplitude using the Neural Network
  if( bunchvector.size()  <=  0 )
    {
      return AliCaloFitResults(AliCaloFitResults::kInvalid, AliCaloFitResults::kInvalid, AliCaloFitResults::kInvalid, AliCaloFitResults::kInvalid , AliCaloFitResults::kInvalid, AliCaloFitResults::kInvalid, AliCaloFitResults::kInvalid );
    } 
 
  short maxampindex;
  short maxamp;

  int index = SelectBunch( bunchvector, &maxampindex , &maxamp ) ;
  
  if( index   < 0 )
    {
      return AliCaloFitResults(AliCaloFitResults::kInvalid, AliCaloFitResults::kInvalid);
    }
  
  int first = 0;
  int last = 0;
 
  Float_t ped = ReverseAndSubtractPed( &(bunchvector.at( index ) )  ,  altrocfg1, altrocfg2, fReversed  );
  
  short maxrev = maxampindex  -  bunchvector.at(index).GetStartBin();
  short timebinOffset = maxampindex - (bunchvector.at(index).GetLength()-1);
  double maxf =  maxamp - ped;

  SelectSubarray( fReversed,  bunchvector.at(index).GetLength(),  maxrev , &first, &last);

  if(maxrev  < 1000 )
    {
      if (  ( maxrev   - first) < 2  &&  (last -   maxrev ) < 2)
	{
	  return AliCaloFitResults( maxamp, ped, AliCaloFitResults::kNoFit, maxf, maxrev+timebinOffset, AliCaloFitResults::kNoFit, AliCaloFitResults::kNoFit,
				    AliCaloFitResults::kNoFit, AliCaloFitSubarray(index, maxrev, first, last) ); 
	}
      else
	{

	  for(int i=0; i < 5 ; i++)
	    {
	      fNNInput[i]  = fReversed[maxrev-2 +i]/(maxamp -ped);
	    } 

	  	  
	  double amp = (maxamp - ped)*fNeuralNet->Value( 0,  fNNInput[0],  fNNInput[1], fNNInput[2], fNNInput[3], fNNInput[4]);
	  double tof = (fNeuralNet->Value( 1,  fNNInput[0],  fNNInput[1], fNNInput[2], fNNInput[3], fNNInput[4]) + timebinOffset ) ;

	  return AliCaloFitResults( maxamp, ped , AliCaloFitResults::kDummy, amp , tof, AliCaloFitResults::kDummy, AliCaloFitResults::kDummy,
				    AliCaloFitResults::kDummy, AliCaloFitSubarray(index, maxrev, first, last) );

	}
    }
  return AliCaloFitResults( maxamp, ped, AliCaloFitResults::kNoFit, maxf, maxrev+timebinOffset, AliCaloFitResults::kNoFit, AliCaloFitResults::kNoFit,
			    AliCaloFitResults::kNoFit, AliCaloFitSubarray(index, maxrev, first, last) ); 
}


