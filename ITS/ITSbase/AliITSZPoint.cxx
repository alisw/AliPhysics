/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliITSZPoint.h"
#include <Riostream.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used by AliITSVertexerZ                                 //
// Contains Z coordinates with their error                       //
// is sortable                                                   //
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliITSZPoint)

//______________________________________________________________________
AliITSZPoint::AliITSZPoint():TObject(),
fZ(0.),
fErrz(0.){
  // Default constructor
}

//______________________________________________________________________
AliITSZPoint::AliITSZPoint(Float_t z, Float_t ez):TObject(),
fZ(z),
fErrz(ez){
  // Standard Constructor
}

//______________________________________________________________________
  AliITSZPoint::~AliITSZPoint(){
  // Destructor
  fZ=0.;
  fErrz=0.;
}

//______________________________________________________________________
void AliITSZPoint::Print(Option_t* /*option */) const {
  // print data members
  printf("Z projection %f , error %f \n",fZ,fErrz);
}
