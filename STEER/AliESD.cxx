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

/*
$Log:
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"

#include "AliESD.h"

ClassImp(AliESD)

//__________________________________________________________________________
AliESD::AliESD() 
{
  cout << "ESD def ctor" << endl;
}

ClassImp(AliESDVertex)

//__________________________________________________________________________
AliESDVertex::AliESDVertex() 
{
  cout << "ESDVertex def ctor" << endl;
  fCoordinates.Set(3);
  fErrorMatrix.Set(6);
}

ClassImp(AliESDTrack)

//__________________________________________________________________________
AliESDTrack::AliESDTrack() 
{
  cout << "ESDTrack def ctor" << endl;
  fPVertex.Set(5);
  fPEVertex.Set(15);
  fPFMeasPoint.Set(6);
  fPFMeasPointErr.Set(15);
  fPLMeasPoint.Set(6);
  fPLMeasPointErr.Set(15);
}
