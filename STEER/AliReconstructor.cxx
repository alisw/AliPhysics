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
// base class for reconstruction algorithms                                  //
//                                                                           //
// Derived classes should implement a default constructor and                //
// the virtual methods                                                       //
// - Reconstruct : to perform the local reconstruction for all events        //
// - FillESD     : to fill the ESD for the current event                     //
//                                                                           //
// The reconstructor classes for the barrel detectors should in addition     //
// implement the method                                                      //
// - CreateTracker : to create a tracker object for the barrel detector      //
//                                                                           //
// The ITS reconstructor should in addition implement the method             //
// - CreateVertexer : to create an object for the vertex finding             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliReconstructor.h"
#include <TString.h>


ClassImp(AliReconstructor)


//_____________________________________________________________________________
const char* AliReconstructor::GetDetectorName() const
{
// get the name of the detector

  static TString detName;
  detName = GetName();
  detName.Remove(0, 3);
  detName.Remove(detName.Index("Reconstructor"));
  return detName.Data();
}
