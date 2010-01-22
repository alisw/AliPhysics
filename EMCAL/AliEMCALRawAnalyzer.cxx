/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <p.t.hille@fys.uio.no>                *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to p.t.hille@fys.uio.no                             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


// Base class for extraction 
// of signal amplitude and peak position
// From EMCAL Calorimeter RAW data (from the RCU)
// Contains some utilities for preparing / selecting
// Signals suitable for signal extraction
// By derived classes

#include "AliLog.h"
#include "AliEMCALRawAnalyzer.h"
#include "AliEMCALBunchInfo.h"
#include "AliEMCALFitResults.h"
#include <iostream>
using namespace std;


AliEMCALRawAnalyzer::AliEMCALRawAnalyzer() :  TObject(),
					      fMinTimeIndex(-1),
					      fMaxTimeIndex(-1),
					      fFitArrayCut(5),
					      fAmpCut(4),
					      fNsampleCut(5),
					      fIsZerosupressed( false ),
					      fVerbose( false )
{
  //Comment
  for(int i=0; i < MAXSAMPLES; i++ )
    {
      fReversed[i] = 0;
    }
}

AliEMCALRawAnalyzer::~AliEMCALRawAnalyzer()
{

}


void 
AliEMCALRawAnalyzer::SetTimeConstraint(const int min, const int max ) 
{
  //Require that the bin if the maximum ADC value is between min and max (timebin)
  if(  ( min > max ) || min > MAXSAMPLES  || max > MAXSAMPLES  )
    {
      AliWarning( Form( "Attempt to set Invalid time bin range (Min , Max) = (%d, %d), Ingored",  min, max ) ); 
    }
  else
    {
      fMinTimeIndex  = min;
      fMaxTimeIndex  = max;
    }
}


UShort_t 
AliEMCALRawAnalyzer::Max(const UShort_t *data, const int length ) const
{
  //------------
  UShort_t tmpmax  = data[0];

  for(int i=0; i < length; i++)
    {
      if( tmpmax  <  data[i] )
	{
	  tmpmax = data[i];
	}
    }
  return tmpmax;
}


void 
AliEMCALRawAnalyzer::SelectSubarray( const Double_t *fData, const int length, const short maxindex, int *const first,  int *const last ) const
{
  //Selection of subset of data from one bunch that will be used for fitting or
  //Peak finding.  Go to the left and right of index of the maximum time bin
  //Untile the ADC value is less that fFitArrayCut
  int tmpfirst =  maxindex;
  int tmplast  =  maxindex;
  
  while((  tmpfirst  ) > 0  &&  ( fData[tmpfirst] >  fFitArrayCut   ))  
    {
      tmpfirst -- ;
    }
  
  while(( tmplast ) <  length   && ( fData [tmplast] >  fFitArrayCut ))
    {
      tmplast ++;
    }
  
  *first = tmpfirst +1;
  *last =  tmplast -1;
}



Float_t 
AliEMCALRawAnalyzer::ReverseAndSubtractPed( const AliEMCALBunchInfo *bunch, const UInt_t /*altrocfg1*/,  const UInt_t /*altrocfg2*/, double *outarray ) const
{
  //Time sample comes in reversed order, revers them back
  //Subtract the baseline based on content of altrocfg1 and altrocfg2.
  Int_t length =  bunch->GetLength(); 
  const UShort_t *sig = bunch->GetData();

  double ped = 0;
  double tmp = 0;

  if( fIsZerosupressed == false ) 
    {
      for(int i=0; i < 5; i++ )
	{
	  tmp += sig[i];
	}
    }
  ped = tmp / 5;
  for( int i=0; i < length; i++ )
    {
      outarray[i] = sig[length -i -1] - ped;
    }

  return ped;
}


short  
AliEMCALRawAnalyzer::Max( const AliEMCALBunchInfo *const bunch , int *const maxindex ) const
{
  //comment
  short tmpmax = -1;
  int tmpindex = -1;
  const UShort_t *sig = bunch->GetData();

  for(int i=0; i < bunch->GetLength(); i++ )
    {
      if( sig[i] > tmpmax )
	{
	  tmpmax   =  sig[i];
	  tmpindex =  i; 
	}
    }
  
  if(maxindex != 0 )
    {
      *maxindex =  bunch->GetLength() -1 - tmpindex + bunch->GetStartBin(); 
    }
  
  return  tmpmax;
}


int 
AliEMCALRawAnalyzer::SelectBunch( const vector<AliEMCALBunchInfo> &bunchvector,short *const maxampbin, short *const maxamplitude ) const
{
  //We select the bunch with the highest amplitude unless any time constraints is set
  short max = -1;
  short bunchindex = -1;
  short maxall = -1;
  int indx = -1;

  for(unsigned int i=0; i < bunchvector.size(); i++ )
    {
      max = Max(  &bunchvector.at(i), &indx );  
      if( IsInTimeRange( indx) )
	{
	  if( max > maxall )
	    {
	      maxall = max;
	      bunchindex = i;
	    }
	}
    }
 
  *maxampbin     = indx;
  *maxamplitude  = max;
  return  bunchindex;
}


bool 
AliEMCALRawAnalyzer::IsInTimeRange( const int maxindex  ) const
{
  // Ckeck if the index of the max ADC vaue is consistent with trigger.
  if( ( fMinTimeIndex  < 0 && fMaxTimeIndex  < 0) ||fMaxTimeIndex  < 0 )
    {
      return true; 
    }
  return ( maxindex < fMaxTimeIndex ) && ( maxindex > fMinTimeIndex  ) ? true : false;
}


void 
AliEMCALRawAnalyzer::PrintBunches( const vector<AliEMCALBunchInfo> &bvctr ) const
{
  //comment
  cout << __FILE__ << __LINE__<< "*************** Printing Bunches *******************" << endl;
  cout << __FILE__ << __LINE__<< "*** There are " << bvctr.size() << ", bunches" << endl;

  for(unsigned int i=0; i < bvctr.size() ; i++ )
    {
      PrintBunch( bvctr.at(i) );
      cout << " bunch = "  <<  i  << endl;
    }
  cout << __FILE__ << __LINE__<< "*************** Done ... *******************" << endl;
}


void 
AliEMCALRawAnalyzer::PrintBunch( const AliEMCALBunchInfo &bunch ) const
{
  //comment
  cout << __FILE__ << ":" << __LINE__ << endl;
  cout << __FILE__ << __LINE__   << ", startimebin = " << bunch.GetStartBin() << ", length =" <<  bunch.GetLength() << endl;
  const UShort_t *sig =  bunch.GetData();  
  
  for ( Int_t j = 0; j <  bunch.GetLength();  j++) 
    {
      printf("%d\t", sig[j] );
    }
  cout << endl; 
}


