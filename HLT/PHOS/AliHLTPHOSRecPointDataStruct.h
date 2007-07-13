#ifndef ALIHLTPHOSRECPOINTDATASTRUCT_H
#define ALIHLTPHOSRECPOINTDATASTRUCT_H

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Øystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

struct AliHLTPHOSRecPointDataStruct
{
  AliHLTUInt8_t fPHOSModule;
  AliHLTUInt8_t fMultiplicity;
  AliHLTUInt8_t fCoordinatesPtr[2]; 
  Float_t fX;
  Float_t fZ;
  Float_t fM2x;
  Float_t fM2z;
  Float_t fM3x;
  Float_t fM4z;
  Float_t fPhixe;
  Float_t fDistanceToBadChannel;
  Float_t* fEnergiesListPtr;

  void New()
  {
    fEnergiesListPtr = new Float_t[fMultiplicity];
  }

  void Del()
  {
    if(fEnergiesListPtr)
      {
	delete [] fEnergiesListPtr;
	fEnergiesListPtr = 0;
      }
  }
};

#endif
