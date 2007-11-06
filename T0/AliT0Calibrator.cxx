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

/* $Id$ */

#include <Riostream.h>

#include "AliLog.h"
#include "AliT0Parameters.h"
#include "AliT0Calibrator.h"
#include <TArrayI.h>
#include <TGraph.h>
#include <TH1F.h>

ClassImp(AliT0Calibrator)

//____________________________________________________________________
  AliT0Calibrator::AliT0Calibrator():TNamed(),
  fChannelWidth(0),  
  fWalk(0)
{
  AliT0Parameters* param = AliT0Parameters::Instance();
  param->Init();
  
  fChannelWidth = param->GetChannelWidth() ;  
 
  for (Int_t i=0; i<24; i++){
    fTimeDelayCFD[i] = Int_t (param->GetTimeDelayCFD(i));
    TGraph* fu = param ->GetWalk(i);
    fWalk.AddAtAndExpand(fu,i);

  }
  
  //
}
//_____________________________________________________________________________

AliT0Calibrator::AliT0Calibrator(const AliT0Calibrator &r): TNamed(),
  fChannelWidth(0),  
  fWalk(0)

{
  //
  // AliT0calibartor copy constructor
  //

  ((AliT0Calibrator &) r).Copy(*this);

}

//_____________________________________________________________________________
AliT0Calibrator &AliT0Calibrator::operator=(const AliT0Calibrator &r)
{
  //
  // Assignment operator
  //

  if (this != &r) ((AliT0Calibrator &) r).Copy(*this);
  return *this;

}



//____________________________________________________________________
Int_t  AliT0Calibrator::WalkCorrection(Int_t ipmt, Int_t qt, Int_t time) 
{
  //slewing correcion and equalizing channels

  Int_t timeEq=0, timeWalk=0;  
  TGraph *fu1=(TGraph*) fWalk.At(ipmt);
  Float_t walk=fu1->Eval(Float_t(qt));
  TH1F*hr=fu1->GetHistogram();
  Float_t maxValue=hr->GetMaximum(50);
  timeWalk = time + Int_t((maxValue-walk)/fChannelWidth) ;
  
  timeEq= timeWalk - (fTimeDelayCFD[ipmt]-fTimeDelayCFD[0]);
  AliDebug(10,Form(" time before %i timeWalk %i ,  qt %i timeEq %i \n ",
		  time,timeWalk, qt, timeEq ));
  return timeEq;
}



