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

///////////////////////////////////////////////////////////////////////////
//    AliITSUGeomTGeo is a simple interface class to TGeoManager       //
//    It is used in the simulation and reconstruction in order to        //
//    query the TGeo ITS geometry                                        //
//                                                                       //
//    author - cvetan.cheshkov@cern.ch                                   //
//    15/02/2007                                                         //
//    adapted to ITSupg 18/07/2012 - ruben.shahoyan@cern.ch              //
//                                                                       //
//    ATTENTION: In opposite to ols AliITSgeomTGeo, all indices start    //
//    from 0, not from 1!!!                                              //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#include <TDatime.h>
#include <TMath.h>
#include <TSystem.h>

#include "AliITSUGeomTGeo.h"
#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliITSMFTSegmentationPix.h"
using namespace TMath;

ClassImp(AliITSUGeomTGeo)

