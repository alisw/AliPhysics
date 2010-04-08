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
// Dummy functions originally written before EMCAL installation and
// survey are kept for backward compatibility, but now they are not
// used.
// Surveyed points on the face of each module are read in from file
// and converted to position of the center and roll, pitch, yaw angles
// of each surveyed SM.
//
// J.L. Klay - Cal Poly
// 07-Apr-2010
//

#include <fstream>

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TString.h>
#include <TMath.h>

#include "AliSurveyObj.h"
#include "AliSurveyPoint.h"

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

  //Coordinates for each SM described in survey reports

  struct AliEMCALSuperModuleCoords {
    Double_t fX1; //x coordinate of the center of supermodule
    Double_t fY1; //y coordinate of the center of supermodule
    Double_t fZ1; //z coordinate of the center of supermodule
    Double_t fPhi;  //roll angle (phi) of supermodule
    Double_t fTheta; //pitch (theta) of supermodule
    Double_t fPsi;   //yaw (psi) of supermodule
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

  AliSurveyObj *s1 = new AliSurveyObj();
  s1->FillFromLocalFile(txtFileName);

  TObjArray* points = s1->GetData();

  fNSuperModule = geom->GetNumberOfSuperModules();

  InitSuperModuleData(points);

  //////////////////////////////////
  //Old way with dummy survey file
  //////////////////////////////////
  //std::ifstream inputFile(txtFileName.Data());
  //if (!inputFile) {
  //  AliError(("Cannot open the survey file " + txtFileName).Data());
  //  return;
  //}
  //
  //Int_t dummyInt = 0;
  //Double_t *xReal = new Double_t[fNSuperModule];
  //Double_t *yReal = new Double_t[fNSuperModule];
  //Double_t *zReal = new Double_t[fNSuperModule];
  //
  //for (Int_t i = 0; i < fNSuperModule; ++i) {
  //  if (!inputFile) {
  //    AliError("Error while reading input file.");
  //    delete [] xReal;
  //    delete [] yReal;
  //    delete [] zReal;
  //    return;
  //  }
  //  inputFile>>dummyInt>>xReal[i]>>yReal[i]>>zReal[i];
  //}
  //
  //InitSuperModuleData(xReal, yReal, zReal);
  //
  //delete [] xReal;
  //delete [] yReal;
  //delete [] zReal;

}

//____________________________________________________________________________
AliEMCALSurvey::~AliEMCALSurvey()
{
  //destructor
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
void AliEMCALSurvey::InitSuperModuleData(const TObjArray *svypts)
{
  //This method uses the data points from the actual EMCAL survey to
  //create the alignment matrices.  Only valid for measured(installed)
  //SM, others will have null objects
  
  
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
    smc.fPhi = phi; //degrees
    smc.fTheta = 0.; //degrees
    smc.fPsi = 0.; //degrees
    if(smodnum%2==0) {
      smc.fZ1 = zpos;
    } else {
      smc.fZ1 = -zpos;
    }

    printf("PHI OF IDEAL SM = %.2f\n",smc.fPhi);

  }

  //Real coordinates of center and rotation angles need to be computed
  //from survey points

  char substr[100];
  AliEMCALSuperModuleCoords *realSM = new AliEMCALSuperModuleCoords[fNSuperModule];
  for (Int_t smodnum = 0; smodnum < geom->GetNumberOfSuperModules(); ++smodnum) {
    AliEMCALSuperModuleCoords &smc = realSM[smodnum];
    zpos = pars[2]; //center of SM in z from geometry

    sprintf(substr,"4096%d",smodnum);
    //retrieve components of four face points and determine average position of center
    //in x,y,z

    std::vector<Double_t> xval;
    std::vector<Double_t> yval;
    std::vector<Double_t> zval;
    
    for(Int_t i = 0; i < svypts->GetEntries(); i++) {
      AliSurveyPoint* pt = (AliSurveyPoint*)svypts->At(i);
      TString ptname = pt->GetPointName();
      if(ptname.Contains(substr)) {
	xval.push_back(pt->GetX()*100.); //convert m to cm
	yval.push_back(pt->GetY()*100.); 
	zval.push_back(pt->GetZ()*100.); 
      }
    }

    //Take average of all relevant pairs of points
    Double_t xReal = ((xval[1]+xval[7])/2.      //4X02 and 4X08
                    + (xval[2]+xval[6])/2.      //4X03 and 4X07
                    + (xval[3]+xval[5])/2.      //4X04 and 4X06
		    + (xval[2]+xval[4])/2.)/4.; //4X03 and 4X05
    Double_t yReal = ((yval[1]+yval[7])/2.      //4X02 and 4X08
                    + (yval[2]+yval[6])/2.      //4X03 and 4X07
                    + (yval[3]+yval[5])/2.      //4X04 and 4X06
		    + (yval[2]+yval[4])/2.)/4.; //4X03 and 4X05
    smc.fX1 = xReal;
    smc.fY1 = yReal;

    //Find average value of z for front face of SM
    Double_t zReal = 0.;
    Int_t nPoints = zval.size();
    for(Int_t iz = 0; iz < nPoints; iz++) {
      zReal += zval[iz];
    }
    if(nPoints > 0) zReal = zReal/nPoints;

    if(smodnum%2==0) {
      smc.fZ1 = zReal-zpos;  //z measured is along end,
    } else {                 //convert to middle of SM
      smc.fZ1 = zReal+zpos;
    }

    Double_t roll = ( TMath::ATan((yval[5]-yval[3])/(xval[5]-xval[3])) //4X04 and 4X06
		      + TMath::ATan((yval[4]-yval[2])/(xval[4]-xval[2])) )/2.; //4X05 and 4X03
    smc.fPhi = 90. + roll*TMath::RadToDeg();

    //Note pitch calc only uses first 8 values, even though 10 are
    //measured on the topmost modules
    Double_t pitch = ( TMath::ATan((zval[0]-zval[1])/(yval[1]-yval[0])) //4X01 and 4X02
		       + TMath::ATan((zval[2]-zval[3])/(yval[3]-yval[2])) //4X03 and 4X04
		       + TMath::ATan((zval[4]-zval[5])/(yval[5]-yval[4])) //4X05 and 4X06
		       + TMath::ATan((zval[6]-zval[7])/(yval[7]-yval[6])) )/4.; //4X07 and 4X08
    smc.fTheta = 0. + pitch*TMath::RadToDeg();

    Double_t yaw = ( TMath::ATan((zval[3]-zval[5])/(xval[5]-xval[3])) //4X04 and 4X06
		     + TMath::ATan((zval[2]-zval[4])/(xval[4]-xval[2])) //4X03 and 4X05
		     + TMath::ATan((zval[1]-zval[7])/(xval[7]-xval[1])) //4X02 and 4X08
		     + TMath::ATan((zval[0]-zval[6])/(xval[6]-xval[0])) )/4.; //4X01 and 4X07
    smc.fPsi = 0. + yaw*TMath::RadToDeg();

  }//loop over supermodules
  
  fSuperModuleData = new AliEMCALSuperModuleDelta[fNSuperModule];
  
  for (Int_t i = 0; i < fNSuperModule; ++i) {
    const AliEMCALSuperModuleCoords &real = realSM[i];
    const AliEMCALSuperModuleCoords &ideal = idealSM[i];
    AliEMCALSuperModuleDelta &t = fSuperModuleData[i];
    t.fXShift = real.fX1 - ideal.fX1;
    t.fYShift = real.fY1 - ideal.fY1;
    t.fZShift = real.fZ1 - ideal.fZ1;
    t.fPhi = real.fPhi - ideal.fPhi;
    t.fTheta = real.fTheta - ideal.fTheta;
    t.fPsi = real.fPsi - ideal.fPsi;

    printf("===================== SM %d =======================\n",i);
    printf("real x (%.2f) - ideal x (%.2f) = shift in x (%.2f)\n",real.fX1,ideal.fX1,t.fXShift);
    printf("real y (%.2f) - ideal y (%.2f) = shift in y (%.2f)\n",real.fY1,ideal.fY1,t.fYShift);
    printf("real z (%.2f) - ideal z (%.2f) = shift in z (%.2f)\n",real.fZ1,ideal.fZ1,t.fZShift);
    printf("real theta (%.2f) - ideal theta (%.2f) = shift in theta %.2f\n",real.fTheta,ideal.fTheta,t.fTheta);
    printf("real psi (%.2f) - ideal psi (%.2f) = shift in psi %.2f\n",real.fPsi,ideal.fPsi,t.fPsi);
    printf("real phi (%.2f) - ideal phi (%.2f) = shift in phi %.2f\n",real.fPhi,ideal.fPhi,t.fPhi);
    printf("===================================================\n");    
  }

  delete [] realSM;
  delete [] idealSM;
}


//____________________________________________________________________________
void AliEMCALSurvey::InitSuperModuleData(const Double_t *xReal, const Double_t *yReal, const Double_t *zReal)
{
  ///////////////////////////////////////
  //Old dummy file way of doing it
  //////////////////////////////////////
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
    zpos = pars[2];
    smc.fX1 = xReal[smodnum];  //x and y match
    smc.fY1 = yReal[smodnum];  //x and y match
    if(smodnum%2==0) {
      smc.fZ1 = zReal[smodnum]-zpos;  //z measured is along end,
    } else {                          //convert to middle of SM
      smc.fZ1 = zReal[smodnum]+zpos;
    }
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

    printf("===================== SM %d =======================\n",i);
    printf("real x (%.2f) - ideal x (%.2f) = shift in x (%.2f)\n",real.fX1,ideal.fX1,t.fXShift);
    printf("real y (%.2f) - ideal y (%.2f) = shift in y (%.2f)\n",real.fY1,ideal.fY1,t.fYShift);
    printf("real z (%.2f) - ideal z (%.2f) = shift in z (%.2f)\n",real.fZ1,ideal.fZ1,t.fZShift);
    printf("theta %.2f\n",t.fTheta);
    printf("psi %.2f\n",t.fPsi);
    printf("phi %.2f\n",t.fPhi);
    printf("===================================================\n");    
  }

  delete [] realSM;
  delete [] idealSM;
}

