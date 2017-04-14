// -*- mode: c++ -*-
/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <perthomas.hille@yale.edu>            *
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

// AliRoot/EMCal system
#include "AliCaloRawAnalyzer.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliLog.h"

// ROOT sytem
#include "TMath.h"

// Standard libraries
#include <iostream>
using namespace std;

/// \cond CLASSIMP
ClassImp(AliCaloRawAnalyzer) ;
/// \endcond
  
///
/// Constructor.
//_______________________________________________________________________
AliCaloRawAnalyzer::AliCaloRawAnalyzer(const char *name, const char *nameshort) :  TObject(),
  fMinTimeIndex(-1),
  fMaxTimeIndex(-1),
  fFitArrayCut(5),
  fAmpCut(4),
  fNsampleCut(5),
  fOverflowCut(950),
  fNsamplePed(3),
  fIsZerosupressed( false ),
  fVerbose( false ),
  fAlgo(Algo::kNONE), 
  fL1Phase(0),
  fAmp(0),
  fTof(0),
  fTau( EMCAL::TAU ),
  fFixTau( true )
{
  snprintf(fName,     256, "%s", name);
  snprintf(fNameShort,256, "%s", nameshort);
  
  for(int i=0; i < ALTROMAXSAMPLES; i++ )
  {
    fReversed[i] = 0;
  }
}

///
/// Require that the bin if the maximum ADC value is between min and max (timebin)
//_______________________________________________________________________
void 
AliCaloRawAnalyzer::SetTimeConstraint(const int min, const int max ) 
{
  
  if(  ( min > max ) || min > ALTROMAXSAMPLES  || max > ALTROMAXSAMPLES  )
    {
      AliWarning( Form( "Attempt to set Invalid time bin range (Min , Max) = (%d, %d), Ingored",  min, max ) ); 
    }
  else
    {
      fMinTimeIndex  = min;
      fMaxTimeIndex  = max;
    }
}

///
/// Get maximum of array
//_______________________________________________________________________
UShort_t 
AliCaloRawAnalyzer::Max(const UShort_t *data, const int length ) const
{  
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

///
/// Selection of subset of data from one bunch that will be used for fitting or
/// Peak finding.  Go to the left and right of index of the maximum time bin
/// Until the ADC value is less that fFitArrayCut, or derivative changes sign (data jump)
//_______________________________________________________________________
void 
AliCaloRawAnalyzer::SelectSubarray( const Double_t *data, int length, short maxindex,
                                    int * first, int * last, int cut) const
{
  int tmpfirst =  maxindex;
  int tmplast  =  maxindex;
  Double_t prevFirst =  data[maxindex];
  Double_t prevLast  =  data[maxindex];  
  bool firstJump = false;
  bool lastJump = false;
  
  while( (tmpfirst >= 0) && (data[tmpfirst] >= cut ) && (!firstJump) ) 
  {
    // jump check:
    if (tmpfirst != maxindex) { // neighbor to maxindex can share peak with maxindex
      if ( data[tmpfirst] >= prevFirst) {
        firstJump = true;
      }
    }
    prevFirst = data[tmpfirst];
    tmpfirst -- ;
  }
  
  while( (tmplast < length) && (data[tmplast] >=  cut ) && (!lastJump) ) 
  {
    // jump check:
    if (tmplast != maxindex) { // neighbor to maxindex can share peak with maxindex
      if ( data[tmplast] >= prevLast) {
        lastJump = true;
      }
    }
    prevLast = data[tmplast];
    tmplast ++;
  }
  
  // we keep one pre- or post- sample if we can (as in online)
  // check though if we ended up on a 'jump', or out of bounds: if so, back up
  if (firstJump || tmpfirst<0) tmpfirst ++;
  if (lastJump || tmplast>=length) tmplast --;
  
  *first = tmpfirst;
  *last =  tmplast;
  
  return;
}


///
/// Time sample comes in reversed order, revers them back
/// Subtract the baseline based on content of altrocfg1 and altrocfg2.
//_______________________________________________________________________
Float_t 
AliCaloRawAnalyzer::ReverseAndSubtractPed( const AliCaloBunchInfo *bunch,
                                           UInt_t /*altrocfg1*/, UInt_t /*altrocfg2*/,
                                           double *outarray ) const
{  
  Int_t length =  bunch->GetLength(); 
  const UShort_t *sig = bunch->GetData();
  
  double ped = EvaluatePedestal( sig , length);
  
  for( int i=0; i < length; i++ )
  {
    outarray[i] = sig[length -i -1] - ped;
  }
  
  return ped;
}

///
/// Pedestal evaluation if not zero suppressed
//_______________________________________________________________________
Float_t 
AliCaloRawAnalyzer::EvaluatePedestal(const UShort_t * data, int /*length*/ ) const
{  
  double tmp = 0;
  
  if( fIsZerosupressed == false ) 
  {
    for(int i=0; i < fNsamplePed; i++ )
    {
      tmp += data[i];
    }
  }
  
  return  tmp/fNsamplePed;
}

///
///  Get maximum in bunch array.
//_______________________________________________________________________
short  
AliCaloRawAnalyzer::Max( const AliCaloBunchInfo * bunch , int * maxindex ) const
{  
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
    //   *maxindex =  bunch->GetLength() -1 - tmpindex + bunch->GetStartBin(); 
    *maxindex =  bunch->GetLength() -1 - tmpindex + bunch->GetStartBin(); 
  }
  
  return  tmpmax;
}

///
/// A bunch is considered invalid if the maximum is in the first or last time-bin.
//_______________________________________________________________________
bool  
AliCaloRawAnalyzer::CheckBunchEdgesForMax( const AliCaloBunchInfo *const bunch ) const
{
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
  
  bool bunchOK = true;
  if (tmpindex == 0 || tmpindex == (bunch->GetLength()-1) )
  {
    bunchOK = false;
  }
  
  return  bunchOK;
}

///
/// We select the bunch with the highest amplitude unless any time constraints is set.
//_______________________________________________________________________
int 
AliCaloRawAnalyzer::SelectBunch( const vector<AliCaloBunchInfo> &bunchvector,
                                 short * maxampbin, short * maxamplitude )
{
  short max = -1;
  short bunchindex = -1;
  short maxall = -1;
  int indx = -1;
  
  for(unsigned int i=0; i < bunchvector.size(); i++ )
  {
    max = Max(  &bunchvector.at(i), &indx ); // CRAP PTH, bug fix, trouble if more than one bunches  
    if( IsInTimeRange( indx, fMaxTimeIndex, fMinTimeIndex) )
    {
      if( max > maxall )
      {
        maxall = max;
        bunchindex = i;
        *maxampbin     = indx;
        *maxamplitude  = max;
      }
    }
  }
  
  if (bunchindex >= 0) 
  {
    bool bunchOK = CheckBunchEdgesForMax( &bunchvector.at(bunchindex) );
    if (! bunchOK)
    { 
      bunchindex = -1; 
    }
  }
  
  return  bunchindex;
}

///
/// Check if the index of the max ADC vaue is consistent with trigger.
//_______________________________________________________________________
bool 
AliCaloRawAnalyzer::IsInTimeRange( const int maxindex, const int maxtindx, const int mintindx ) const
{
  if( ( mintindx  < 0 && maxtindx   < 0) ||maxtindx  < 0 )
  {
    return true; 
  }
  
  return ( maxindex < maxtindx ) && ( maxindex > mintindx  ) ? true : false;
}

///
/// Print bunch vector infomation.
//_______________________________________________________________________
void 
AliCaloRawAnalyzer::PrintBunches( const vector<AliCaloBunchInfo> &bvctr ) 
{
  cout << __FILE__ << __LINE__<< "*************** Printing Bunches *******************" << endl;
  cout << __FILE__ << __LINE__<< "*** There are " << bvctr.size() << ", bunches" << endl;
  
  for(unsigned int i=0; i < bvctr.size() ; i++ )
  {
    PrintBunch( bvctr.at(i) );
    cout << " bunch = "  <<  i  << endl;
  }
  cout << __FILE__ << __LINE__<< "*************** Done ... *******************" << endl;
}

///
/// Print bunch information.
//_______________________________________________________________________
void 
AliCaloRawAnalyzer::PrintBunch( const AliCaloBunchInfo &bunch )
{
  cout << __FILE__ << ":" << __LINE__ << endl;
  cout << __FILE__ << __LINE__   << ", startimebin = " << bunch.GetStartBin() << ", length =" <<  bunch.GetLength() << endl;
  const UShort_t *sig =  bunch.GetData();  
  
  for ( Int_t j = 0; j <  bunch.GetLength();  j++) 
  {
    printf("%d\t", sig[j] );
  }
  cout << endl; 
}

///
///   \param amp    - max amplitude;
///   \param time   - time of max amplitude; 
///   \param first  - sample array indices to be used
///   \param last   - sample array indices to be used
///   \param adcErr - nominal error of amplitude measurement (one value for all channels) if adcErr<0 that mean adcErr=1.
///   \param tau    - filter time response (in timebin units)
///
///   \return chi2  
//_______________________________________________________________________
Double_t
AliCaloRawAnalyzer::CalculateChi2( Double_t amp,    Double_t time,
                                   Int_t    first,  Int_t    last,
                                   Double_t adcErr, Double_t tau) const
{
  if (first == last || first<0 ) { // signal consists of single sample, chi2 estimate (0) not too well defined.. 
                                   // or, first is negative, the indices are not valid
    return Ret::kDummy;
  }
  
  int nsamples =  last - first + 1;
  // printf(" AliCaloRawAnalyzer::CalculateChi2 : first %i last %i : nsamples %i : amp %3.2f time %3.2f \n", first, last, nsamples, amp, time); 
  
  Int_t x = 0;
  Double_t chi2 = 0;
  Double_t dy = 0.0, xx = 0.0, f=0.0;
  
  for (Int_t i=0; i<nsamples; i++) 
  {
    x     = first + i; // timebin
    xx    = (x - time + tau) / tau; // help variable
    f     = 0;
    
    if (xx > 0) 
    {
      f = amp * xx*xx * TMath::Exp(2 * (1 - xx )) ;
    }
    
    dy    = fReversed[x] - f; 
    chi2 += dy*dy;
    // printf(" AliCaloRawAnalyzer::CalculateChi2 : %i : y %f -> f %f : dy %f \n", i, fReversed[first+i], f, dy); 
  }
  
  if (adcErr>0.0) 
  { // weight chi2
    chi2 /= (adcErr*adcErr);
  }
  
  return chi2;
}

///
/// * Input:
///   \param first: sample array indices to be used
///   \param last: sample array indices to be used 
/// * Output:
///   \param mean: of samples
///   \param rms: of samples 
///
/// To possibly be used to differentiate good signals from bad before fitting.
//_______________________________________________________________________
void
AliCaloRawAnalyzer::CalculateMeanAndRMS(Int_t first, Int_t last,
                                        Double_t & mean, Double_t & rms)
{
  mean = Ret::kDummy;
  rms =  Ret::kDummy;
  
  // signal consists of single sample, chi2 estimate (0) not too well defined.. 
  // or, first is negative, the indices are not valid
  if (first == last || first<0 ) 
    return;
  
  int nsamples =  last - first + 1;
  //  printf(" AliCaloRawAnalyzer::CalculateMeanAndRMS : first %i last %i : nsamples %i \n", first, last, nsamples); 
  
  int x = 0;
  Double_t sampleSum = 0; // sum of samples
  Double_t squaredSampleSum = 0; // sum of samples squared
  
  for (Int_t i=0; i<nsamples; i++) 
  {
    x = first + i;
    sampleSum += fReversed[x];
    squaredSampleSum += (fReversed[x] * fReversed[x]);
  }
  
  mean = sampleSum / nsamples; 	 
  Double_t squaredMean = squaredSampleSum / nsamples; 	 
  // The variance (rms squared) is equal to the mean of the squares minus the square of the mean.. 	 
  rms = sqrt(squaredMean - mean*mean);
  
  return;
}


///
/// Method to do the selection of what should possibly be fitted.
//_______________________________________________________________________
int
AliCaloRawAnalyzer::PreFitEvaluateSamples( const vector<AliCaloBunchInfo>  &bunchvector,
                                          UInt_t altrocfg1, UInt_t altrocfg2, Int_t & index,
                                          Float_t & maxf, short & maxamp,
                                          short & maxrev, Float_t & ped,
                                          int & first, int & last, int acut )
{  
  int nsamples  = 0;
  short maxampindex = 0;
  
  // select the bunch with the highest amplitude unless any time constraints is set
  index = SelectBunch( bunchvector,  &maxampindex,  &maxamp ); 
  
  // something valid was found, and non-zero amplitude
  if( index >= 0 && maxamp >= acut ) 
  {
    // use more convenient numbering and possibly subtract pedestal
    ped  = ReverseAndSubtractPed( &(bunchvector.at(index)),  altrocfg1, altrocfg2, fReversed  );
    maxf = TMath::MaxElement( bunchvector.at(index).GetLength(),  fReversed );
    
    if ( maxf >= acut  ) // possibly significant signal
    {
      // select array around max to possibly be used in fit
      maxrev = maxampindex - bunchvector.at(index).GetStartBin(); 
      SelectSubarray( fReversed,  bunchvector.at(index).GetLength(),  maxrev, &first, &last, acut );
      
      // sanity check: maximum should not be in first or last bin
      // if we should do a fit
      if (first!=maxrev && last!=maxrev) 
      {
        // calculate how many samples we have 
        nsamples =  last - first + 1;	  
      }
    }
  }
  
  return nsamples;
}

