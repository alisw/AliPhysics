/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with FIT reconstruction parameters                                  //
//    fMeanAmplitude -
//           for low flux time-amplitude correction equalize time to amplitude 1 MIP; 
//           for high flux - to 15MIP   
//    To have nice time spectra after reconstruction we need to know 
//    reference point to write t(i) - RefPoint. 
//    It can be apparatus RefPoint or one of PMT                       //  
//    fRefPoint - number of channel with RF
//
//       Alla.Maevskaya@cern.ch
/////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
 
#include "AliFITRecoParam.h"
#include "Riostream.h"

using std::cout;
using std::endl;
ClassImp(AliFITRecoParam)




//_____________________________________________________________________________
AliFITRecoParam::AliFITRecoParam():
  AliDetectorRecoParam()

{
  //
  // constructor
  SetName("FIT");
  SetTitle("FIT");

  for (Int_t i=0; i<240; i++) fBadChannels[i]=-1;
  
 
}

//_____________________________________________________________________________
AliFITRecoParam::~AliFITRecoParam() 
{
  //
  // destructor
  //  
}

//_____________________________________________________________________________

AliFITRecoParam::AliFITRecoParam(const AliFITRecoParam &p):
  AliDetectorRecoParam(p)      
{ 
 //copy constructor
  for (Int_t i=0; i<240; i++)  
    fBadChannels[i] = p.fBadChannels[i];
  
}
//_____________________________________________________________________________

AliFITRecoParam& AliFITRecoParam:: operator=(const AliFITRecoParam &p)
{
  //
  // assign. operator
  //

  if (this == &p)
    return *this;
  
  AliDetectorRecoParam::operator=(p);
  for (Int_t i=0; i<240; i++)  fBadChannels[i] = p.fBadChannels[i];
  
 
 return *this;

}
//_____________________________________________________________________________
 
AliFITRecoParam *AliFITRecoParam::GetLowFluxParam()
{
  //
  // make default reconstruction  parameters for low  flux env.
  //
  AliFITRecoParam *param = new AliFITRecoParam();
   return param;
}

//_____________________________________________________________________________

AliFITRecoParam *AliFITRecoParam::GetHighFluxParam()
{
  //
  // make reco parameters for high flux env.
  //

  AliFITRecoParam *param = new AliFITRecoParam();

   param->SetName("High Flux");
  param->SetTitle("High Flux");
  return param;
}


//_____________________________________________________________________________

AliFITRecoParam *AliFITRecoParam::GetLaserTestParam()
{
  //
  // special setting for laser
  //
  AliFITRecoParam *param = new AliFITRecoParam();
   param->SetName("Laser Flux");
  param->SetTitle("Laser Flux");
  return param;
}
//_____________________________________________________________________________

void AliFITRecoParam::PrintParameters() const
{
  //
  // Printing of the used FIT reconstruction parameters
  //
   cout<<" AliFITRecoParam::PrintParameters() "<<endl;

}
