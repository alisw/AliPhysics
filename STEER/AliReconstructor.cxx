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


#include "AliLog.h"
#include "AliReconstructor.h"
#include <TString.h>


ClassImp(AliReconstructor)


//_____________________________________________________________________________
void AliReconstructor::ConvertDigits(AliRawReader* /*rawReader*/, 
				     TTree* /*digitsTree*/) const
{
// convert raw data digits into digit objects in a root tree

  AliError("conversion of raw data digits into digit objects not implemented");
}


//_____________________________________________________________________________
void AliReconstructor::Reconstruct(TTree* /*digitsTree*/,
				   TTree* /*clustersTree*/) const
{
// run the local reconstruction

  AliError("local event reconstruction not implemented");
}

//_____________________________________________________________________________
void AliReconstructor::Reconstruct(AliRawReader* /*rawReader*/, 
				   TTree* /*clustersTree*/) const
{
// run the local reconstruction with raw data input

  AliError("local event reconstruction not implemented for raw data input");
}

//_____________________________________________________________________________
void AliReconstructor::FillESD(TTree* /*digitsTree*/, TTree* /*clustersTree*/,
			       AliESDEvent* /*esd*/) const
{
// fill the ESD.
// by default nothing is done

}

//_____________________________________________________________________________
void AliReconstructor::FillESD(AliRawReader* /*rawReader*/, 
			       TTree* clustersTree, AliESDEvent* esd) const
{
// fill the ESD in case of raw data input.
// by default the FillESD method for MC is called

  FillESD((TTree*)NULL, clustersTree, esd);
}

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
