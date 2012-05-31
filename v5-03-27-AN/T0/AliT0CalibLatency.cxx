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

/* $Id: AliT0CalibLatency.cxx 39062 2010-02-22 09:17:47Z alla $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for T0 calibration                       TM-AC-AM_6-02-2006  
// equalize time shift for each time CFD channel
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0CalibLatency.h"
#include "AliLog.h"
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TProfile.h>
#include <iostream>

ClassImp(AliT0CalibLatency)

//________________________________________________________________
  AliT0CalibLatency::AliT0CalibLatency():TNamed(),
  fLatencyL1(0),
  fLatencyL1A(0),
  fLatencyL1C(0),
  fLatencyHPTDC(9000)

{
  //

}

//________________________________________________________________
AliT0CalibLatency::AliT0CalibLatency(const char* name):TNamed(),
  fLatencyL1(0),
  fLatencyL1A(0),
  fLatencyL1C(0),
  fLatencyHPTDC(9000)


{
  //constructor

  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

 }

//________________________________________________________________
AliT0CalibLatency::AliT0CalibLatency(const AliT0CalibLatency& calibda):TNamed(calibda),		
  fLatencyL1(0),
  fLatencyL1A(0),
  fLatencyL1C(0),
  fLatencyHPTDC(9000)
{
  // copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());


}
//________________________________________________________________
AliT0CalibLatency &AliT0CalibLatency::operator =(const AliT0CalibLatency& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
 
  return *this;
}

//________________________________________________________________
AliT0CalibLatency::~AliT0CalibLatency()
{
  //
  // destrictor
}

//________________________________________________________________
void  AliT0CalibLatency::Print(Option_t*) const
{
  // print time values

  printf("\n	----	Latencies	----\n\n");
  printf(" Latency HPTDC %f A&C %f A %f C %f\n", fLatencyHPTDC,fLatencyL1, fLatencyL1A, fLatencyL1C );
} 

