
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for FIT calibration                       TM-AC-AM_6-02-2006  
// equalize time shift for each time CFD channel
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliFITCalibTimeEq.h"
#include "AliLog.h"
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TProfile.h>
#include <iostream>

ClassImp(AliFITCalibTimeEq)

//________________________________________________________________
AliFITCalibTimeEq::AliFITCalibTimeEq():TNamed()          
{
  //
  for(Int_t i=0; i<200; i++) {
    fTimeEq[i] = 0;	      // Time Equalized for OCDB	 
    fCFDvalue[i] = 0;
  }
}

//________________________________________________________________
AliFITCalibTimeEq::AliFITCalibTimeEq(const char* name):TNamed()
{
  //constructor
  
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  for(Int_t i=0; i<200; i++) {
    fTimeEq[i] = 0;	      // Time Equalized for OCDB	 
    fCFDvalue[i] = 0;
  }
}

//________________________________________________________________
AliFITCalibTimeEq::AliFITCalibTimeEq(const AliFITCalibTimeEq& calibda):TNamed(calibda)     
{
  // copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  ((AliFITCalibTimeEq &) calibda).Copy(*this);
  
}

//________________________________________________________________
AliFITCalibTimeEq &AliFITCalibTimeEq::operator =(const AliFITCalibTimeEq& calibda)
{
  // assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  if (this != &calibda) (( AliFITCalibTimeEq &) calibda).Copy(*this);
  return *this;
}

//________________________________________________________________
AliFITCalibTimeEq::~AliFITCalibTimeEq()
{
  //
  // destrictor
}
//________________________________________________________________
void AliFITCalibTimeEq::Reset()
{
  //reset values
  
  memset(fCFDvalue,0,200*sizeof(Float_t));
  memset(fTimeEq,1,200*sizeof(Float_t));
}


//________________________________________________________________
void  AliFITCalibTimeEq::Print(Option_t*) const
{
  // print time values
  
  printf("\n	----	PM Arrays	----\n\n");
  printf(" Time delay CFD \n");
  for (Int_t i=0; i<200; i++) 
    printf(" CFD  %f diff %f  \n",fCFDvalue[i],fTimeEq[i]);
} 


 

