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

/* $Header$ */

//-------------------------------------------------------------------------
//                Implementation of the AliKalmanTrack class
//
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "AliKalmanTrack.h"

ClassImp(AliKalmanTrack)

Double_t AliKalmanTrack::fgConvConst;

//_______________________________________________________________________
AliKalmanTrack::AliKalmanTrack():
  fLab(-3141593),
  fChi2(0),
  fMass(0.13957),
  fN(0)
{
  //
  // Default constructor
  //
    if (fgConvConst==0) 
      Fatal("AliKalmanTrack()","The magnetic field has not been set !\n"); 
}

//_______________________________________________________________________
AliKalmanTrack::AliKalmanTrack(const AliKalmanTrack &t):
  TObject(t),
  fLab(t.fLab),
  fChi2(t.fChi2),
  fMass(t.fMass),
  fN(t.fN)
{
  //
  // Copy constructor
  //
  if (fgConvConst==0) 
    Fatal("AliKalmanTrack(const AliKalmanTrack&)",
          "The magnetic field has not been set !\n"); 
}

