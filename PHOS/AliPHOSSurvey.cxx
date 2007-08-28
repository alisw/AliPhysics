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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.5  2007/07/10 12:41:38  kharlov
 * Added a new class AliPHOSSurvet1 which read survey data from EDMS files
 *
 * Revision 1.4  2007/05/17 17:34:54  kharlov
 * Merging differences if v1.2 and 1.3
 *
 * Revision 1.3  2007/05/17 17:13:32  kharlov
 * Coding convensions satisfied (T.Pocheptsov)
 *
 * Revision 1.1  2007/04/19 15:47:20  kharlov
 * Add misalignment of strip units with AliPHOSSurvey class
 *
 */

// Objects of this class read txt file with survey (photogrammetry) data
// and convert the data into AliAlignObjParams of alignable PHOS volumes.
// It can be used as a base class, you need to override GetStripTransformation.
// AliPHOSSurvey inherits TObject only to use AliLog "functions".
// Author: Timur Pocheptsov (JINR)

#include <fstream>

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TString.h>
#include <TMath.h>

#include "AliSurveyObj.h"

#include "AliPHOSEMCAGeometry.h"
#include "AliAlignObjParams.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSSurvey.h"
#include "AliLog.h"

ClassImp(AliPHOSSurvey)

//____________________________________________________________________________
AliPHOSSurvey::AliPHOSSurvey()
                : fStrNum(0),
                  fStripData(0)
{
  //Default constructor.
}

namespace {

  struct AliPHOSStripCoords {
    Double_t fX1; //x coordinate of the first strip point
    Double_t fZ1; //z coordinate of the first strip point
    Double_t fX2; //x coordinate of the second strip point
    Double_t fZ2; //z coordinate of the second strip point
  };

}

//____________________________________________________________________________
AliPHOSSurvey::AliPHOSSurvey(const TString &txtFileName)
                : fStrNum(0),
                  fStripData(0)
{
  //Read survey data from txt file.
  const AliPHOSGeometry *phosGeom = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  if (!phosGeom) {
    AliError("Cannot obtain AliPHOSGeometry instance.");
    return;
  }

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    AliError(("Cannot open the survey file " + txtFileName).Data());
    return;
  }

  AliPHOSEMCAGeometry * emcaGeom = phosGeom->GetEMCAGeometry();
  fStrNum = emcaGeom->GetNStripX() * emcaGeom->GetNStripZ();

  Int_t dummyInt = 0;
  Double_t dummyY = 0.;
  Double_t *xReal = new Double_t[fStrNum * 2];//2
  Double_t *zReal = new Double_t[fStrNum * 2];//3

  for (Int_t i = 0; i < fStrNum * 2; ++i) {
    if (!inputFile) {
      AliError("Error while reading input file.");
      delete [] xReal;
      delete [] zReal;
      return;
    }
    inputFile>>dummyInt>>xReal[i]>>dummyY>>zReal[i];
    xReal[i] *= 0.1;
    zReal[i] *= 0.1;
  }

  InitStripData(xReal, zReal);

  delete [] zReal;
  delete [] xReal;
}

//____________________________________________________________________________
AliPHOSSurvey::~AliPHOSSurvey()
{
  delete [] fStripData;
}

//____________________________________________________________________________
void AliPHOSSurvey::CreateAliAlignObjParams(TClonesArray &array)
{
  //Create AliAlignObjParams.
  const AliPHOSGeometry * phosGeom = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  if (!phosGeom) {
    AliError("Cannot obtain AliPHOSGeometry instance.");
    return;
  }

  if (!gGeoManager) {
    AliWarning("Cannot create local transformations for strip units - gGeoManager does not exist.");
    AliInfo("Null shifts and rotations will be created instead.");
    return CreateNullObjects(array, phosGeom);
  }

  AliPHOSEMCAGeometry * emcaGeom = phosGeom->GetEMCAGeometry();
  Int_t arrayInd = array.GetEntries(), iIndex = 0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  for (Int_t module = 1; module <= phosGeom->GetNModules(); ++module) {
    for (Int_t i = 0, stripNum = 0; i < emcaGeom->GetNStripX(); ++i) {
      for (Int_t j = 0; j < emcaGeom->GetNStripZ(); ++j) {
        TString stripName(TString::Format("PHOS/Module%d/Strip_%d_%d", module, i, j));
        AliPHOSStripDelta t(GetStripTransformation(stripNum++, module));
        new(array[arrayInd])
          AliAlignObjParams(
                            stripName.Data(), volid, 
                            t.fXShift, t.fYShift, t.fZShift, 
                            -t.fPsi, -t.fTheta, -t.fPhi, 
                            false
                           );
        ++arrayInd;
      }
    }
  }
}

//____________________________________________________________________________
void AliPHOSSurvey::CreateNullObjects(TClonesArray &array, const AliPHOSGeometry *phosGeom)const
{
  //Create null shifts and rotations.
  const AliPHOSEMCAGeometry * emcaGeom = phosGeom->GetEMCAGeometry();
  Int_t arrayInd = array.GetEntries(), iIndex = 0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  for (Int_t module = 1; module <= phosGeom->GetNModules(); ++module)
    for (Int_t i = 0; i < emcaGeom->GetNStripX(); ++i)
      for (Int_t j = 0; j < emcaGeom->GetNStripZ(); ++j) {
        TString stripName(TString::Format("PHOS/Module%d/Strip_%d_%d", module, i, j));
        new(array[arrayInd]) AliAlignObjParams(stripName.Data(), volid, 0., 0., 0., 0., 0., 0., true);
        ++arrayInd;
      }
}

//____________________________________________________________________________
AliPHOSSurvey::AliPHOSStripDelta AliPHOSSurvey::GetStripTransformation(Int_t stripIndex, Int_t module)const
{
  //Strip 'stripIndex' transformation.
  AliPHOSStripDelta t = {0., 0., 0., 0., 0., 0.};
  if (module != 3 || !fStripData)
    return t;
  return fStripData[stripIndex];
}

//____________________________________________________________________________
void AliPHOSSurvey::InitStripData(const Double_t *xReal, const Double_t *zReal)
{
  const AliPHOSGeometry *phosGeom = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  AliPHOSEMCAGeometry * emcaGeom = phosGeom->GetEMCAGeometry();
  const Float_t *strip = emcaGeom->GetStripHalfSize();
  const Float_t *cell = emcaGeom->GetSteelCellHalfSize(); 

  AliPHOSStripCoords *idealStrips = new AliPHOSStripCoords[fStrNum];//1
  for (Int_t ix = 0, stripNumber = 0; ix < emcaGeom->GetNStripX(); ++ix) {
    for (Int_t iz = 0; iz < emcaGeom->GetNStripZ(); ++iz) {
      AliPHOSStripCoords &str = idealStrips[stripNumber++];
      str.fX1 = ix * 2 * strip[0];
      str.fX2 = str.fX1 + 14 * cell[0];
      str.fZ1 = iz * 2 * strip[2];
      str.fZ2 = str.fZ1 + 2 * cell[2];
    }
  }

  AliPHOSStripCoords *realStrips = new AliPHOSStripCoords[fStrNum];//4
  for (Int_t j = 0, stripNumber = 0; j < emcaGeom->GetNStripX() * 2; j += 2) {
    for (Int_t i = 0; i < emcaGeom->GetNStripZ(); ++i) {
      AliPHOSStripCoords &str = realStrips[stripNumber++];
      str.fX1 = xReal[i + j * emcaGeom->GetNStripZ()];
      str.fZ1 = zReal[i + j * emcaGeom->GetNStripZ()];
      str.fX2 = xReal[i + (j + 1) * emcaGeom->GetNStripZ()];
      str.fZ2 = zReal[i + (j + 1) * emcaGeom->GetNStripZ()];
    }
  }

  fStripData = new AliPHOSStripDelta[fStrNum];
  
  for (Int_t i = 0; i < fStrNum; ++i) {
    const AliPHOSStripCoords &real = realStrips[i];
    const AliPHOSStripCoords &ideal = idealStrips[i];
    AliPHOSStripDelta &t = fStripData[i];
    t.fTheta = TMath::ATan((real.fZ2 - real.fZ1)  / (real.fX2 - real.fX1)) - 
               TMath::ATan((ideal.fZ2 - ideal.fZ1) / (ideal.fX2 - ideal.fX1));
    t.fTheta *= TMath::RadToDeg();
    t.fXShift = (real.fX1 + real.fX2) / 2 - (ideal.fX1 + ideal.fX2) / 2;
    t.fZShift = (real.fZ1 + real.fZ2) / 2 - (ideal.fZ1 + ideal.fZ2) / 2;
    t.fYShift = 0., t.fPsi = 0., t.fPhi = 0.;
  }

  delete [] realStrips;
  delete [] idealStrips;
}
