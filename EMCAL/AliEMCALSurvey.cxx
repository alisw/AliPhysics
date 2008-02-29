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

/* $Id: $ */

// Objects of this class read txt file with survey data
// and convert the data into AliAlignObjParams of alignable EMCAL volumes.
// AliEMCALSurvey inherits TObject only to use AliLog "functions".
//
// This whole thing is just a dummy class with dummy variables until
// the actual survey points are determined.  For now I assumed that
// the x,y center and starting z value of each supermodule is the only
// surveyed point
//
// J.L. Klay - Cal Poly
// 28-Feb-2008
//

#include <fstream>

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TString.h>
#include <TMath.h>

#include "AliSurveyObj.h"

#include "AliAlignObjParams.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALSurvey.h"
#include "AliLog.h"

ClassImp(AliEMCALSurvey)

//____________________________________________________________________________
AliEMCALSurvey::AliEMCALSurvey()
                : fNSuperModule(0),
                  fSuperModuleData(0)
{
  //Default constructor.
}

namespace {

  //for now, measurements assumed to be taken at x,y center and
  //nominal starting z value of each supermodule

  struct AliEMCALSuperModuleCoords {
    Double_t fX1; //x coordinate of the first supermodule point
    Double_t fY1; //y coordinate of the first supermodule point
    Double_t fZ1; //z coordinate of the first supermodule point
  };

}

//____________________________________________________________________________
AliEMCALSurvey::AliEMCALSurvey(const TString &txtFileName)
                : fNSuperModule(0),
                  fSuperModuleData(0)
{
  //Read survey data from txt file.
  const AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom) {
    AliError("Cannot obtain AliEMCALGeometry instance.");
    return;
  }

  std::ifstream inputFile(txtFileName.Data());
  if (!inputFile) {
    AliError(("Cannot open the survey file " + txtFileName).Data());
    return;
  }

  fNSuperModule = geom->GetNumberOfSuperModules();

  Int_t dummyInt = 0;
  Double_t *xReal = new Double_t[fNSuperModule];
  Double_t *yReal = new Double_t[fNSuperModule];
  Double_t *zReal = new Double_t[fNSuperModule];

  for (Int_t i = 0; i < fNSuperModule; ++i) {
    if (!inputFile) {
      AliError("Error while reading input file.");
      delete [] xReal;
      delete [] yReal;
      delete [] zReal;
      return;
    }
    inputFile>>dummyInt>>xReal[i]>>yReal[i]>>zReal[i];
  }

  InitSuperModuleData(xReal, yReal, zReal);

  delete [] xReal;
  delete [] yReal;
  delete [] zReal;
}

//____________________________________________________________________________
AliEMCALSurvey::~AliEMCALSurvey()
{
  delete [] fSuperModuleData;
}

//____________________________________________________________________________
void AliEMCALSurvey::CreateAliAlignObjParams(TClonesArray &array)
{
  //Create AliAlignObjParams.
  const AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  if (!geom) {
    AliError("Cannot obtain AliEMCALGeometry instance.");
    return;
  }

  if (!gGeoManager) {
    AliWarning("Cannot create local transformations for supermodules - gGeoManager does not exist.");
    AliInfo("Null shifts and rotations will be created instead.");
    return CreateNullObjects(array, geom);
  }

  Int_t arrayInd = array.GetEntries(), iIndex = 0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  for (Int_t smodnum = 0; smodnum < geom->GetNumberOfSuperModules(); ++smodnum) {
    TString smodName(TString::Format("EMCAL/FullSupermodule%d", smodnum+1));
    if(geom->GetKey110DEG() && smodnum >= 10) {
      smodName = "EMCAL/HalfSupermodule";
      smodName += (smodnum-10+1);
    }    
    AliEMCALSuperModuleDelta t(GetSuperModuleTransformation(smodnum));
    new(array[arrayInd])
      AliAlignObjParams(
			smodName.Data(), volid, 
			t.fXShift, t.fYShift, t.fZShift, 
			-t.fPsi, -t.fTheta, -t.fPhi, 
			false
			);
    ++arrayInd;
  }

}

//____________________________________________________________________________
void AliEMCALSurvey::CreateNullObjects(TClonesArray &array, const AliEMCALGeometry *geom)const
{
  //Create null shifts and rotations.
  Int_t arrayInd = array.GetEntries(), iIndex = 0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  for (Int_t smodnum = 0; smodnum < geom->GetNumberOfSuperModules(); ++smodnum) {
    TString smodName(TString::Format("EMCAL/FullSupermodule%d", smodnum+1));
    if(geom->GetKey110DEG() && smodnum >= 10) {
      smodName = "EMCAL/HalfSupermodule";
      smodName += (smodnum-10+1);
    }
    new(array[arrayInd]) AliAlignObjParams(smodName.Data(), volid, 0., 0., 0., 0., 0., 0., true);
    ++arrayInd;
  }
}

//____________________________________________________________________________
AliEMCALSurvey::AliEMCALSuperModuleDelta AliEMCALSurvey::GetSuperModuleTransformation(Int_t supModIndex)const
{
  //Supermodule transformation.
  AliEMCALSuperModuleDelta t = {0., 0., 0., 0., 0., 0.};
  if (!fSuperModuleData)
    return t;

  return fSuperModuleData[supModIndex];
}

//____________________________________________________________________________
void AliEMCALSurvey::InitSuperModuleData(const Double_t *xReal, const Double_t *yReal, const Double_t *zReal)
{

  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  //Center of supermodules
  Float_t *pars = geom->GetSuperModulesPars();
  Double_t rpos = (geom->GetEnvelop(0) + geom->GetEnvelop(1))/2.;
  Double_t phi, phiRad, xpos, ypos, zpos;

  AliEMCALSuperModuleCoords *idealSM = new AliEMCALSuperModuleCoords[fNSuperModule];
  for (Int_t smodnum = 0; smodnum < geom->GetNumberOfSuperModules(); ++smodnum) {
    AliEMCALSuperModuleCoords &smc = idealSM[smodnum];
    phiRad = geom->GetPhiCenterOfSM(smodnum); //comes in radians
    phi = phiRad*180./TMath::Pi(); //need degrees for AliAlignObjParams
    xpos = rpos * TMath::Cos(phiRad);
    ypos = rpos * TMath::Sin(phiRad);
    zpos = pars[2];
    if(geom->GetKey110DEG() && smodnum >= 10) {
      xpos += (pars[1]/2. * TMath::Sin(phiRad));
      ypos -= (pars[1]/2. * TMath::Cos(phiRad));
    }
    smc.fX1 = xpos;
    smc.fY1 = ypos;
    if(smodnum%2==0) {
      smc.fZ1 = zpos;
    } else {
      smc.fZ1 = -zpos;
    }

  }
  
  AliEMCALSuperModuleCoords *realSM = new AliEMCALSuperModuleCoords[fNSuperModule];
  for (Int_t smodnum = 0; smodnum < geom->GetNumberOfSuperModules(); ++smodnum) {
    AliEMCALSuperModuleCoords &smc = realSM[smodnum];
    smc.fX1 = xReal[smodnum];  //x and y match
    smc.fY1 = yReal[smodnum];  //x and y match
    smc.fZ1 = zReal[smodnum]/2.;  //z measured is along end, need to
			          //convert to middle of SM
  }
  
  fSuperModuleData = new AliEMCALSuperModuleDelta[fNSuperModule];
  
  for (Int_t i = 0; i < fNSuperModule; ++i) {
    const AliEMCALSuperModuleCoords &real = realSM[i];
    const AliEMCALSuperModuleCoords &ideal = idealSM[i];
    AliEMCALSuperModuleDelta &t = fSuperModuleData[i];
    t.fTheta = TMath::ATan(real.fZ1  / real.fX1) - 
               TMath::ATan(ideal.fZ1 / ideal.fX1);
    t.fTheta *= TMath::RadToDeg();
    t.fPsi = 0.;
    t.fPhi = TMath::ATan(real.fY1 / real.fX1) - TMath::ATan(ideal.fY1 / ideal.fX1);
    t.fPhi *= TMath::RadToDeg();
    t.fXShift = real.fX1 - ideal.fX1;
    t.fYShift = real.fY1 - ideal.fY1;
    t.fZShift = real.fZ1 - ideal.fZ1;

  }

  delete [] realSM;
  delete [] idealSM;
}
