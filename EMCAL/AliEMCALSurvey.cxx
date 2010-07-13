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
//
// Surveyed points on the EMCAL support rails were used with the CATIA
// 3D graphics program to determine the positions of the bottom
// corners of the active area for each supermodule.  These numbers are
// read in from file and converted to position of the center and roll, 
// pitch, yaw angles of each installed SM.
//
// J.L. Klay - Cal Poly
// 21-May-2010
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
   fSuperModuleData(0),
   fDataType(kSurvey)
{
  //Default constructor.
}

namespace {

  //Coordinates for each SM described in survey reports

  struct AliEMCALSuperModuleCoords {
    Double_t fX1; //x coordinate of the center of supermodule
    Double_t fY1; //y coordinate of the center of supermodule
    Double_t fZ1; //z coordinate of the center of supermodule
    Double_t fPsi;   //yaw (psi) of supermodule
    Double_t fTheta; //pitch (theta) of supermodule
    Double_t fPhi;  //roll angle (phi) of supermodule

  };

}

//____________________________________________________________________________
AliEMCALSurvey::AliEMCALSurvey(const TString &txtFileName,const SurveyDataType_t type)
  : fNSuperModule(0),
    fSuperModuleData(0),
    fDataType(type)
{
  //Get the geometry object and then attempt to
  //read survey data from a file, depending on which
  //method (kSurvey or kDummy) is selected.

  const AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom) {
    AliError("Cannot obtain AliEMCALGeometry instance.");
    return;
  }

  fNSuperModule = geom->GetNumberOfSuperModules();

  if(fDataType == kSurvey) {

    AliSurveyObj *s1 = new AliSurveyObj();
    s1->FillFromLocalFile(txtFileName);
    TObjArray* points = s1->GetData();
    InitSuperModuleData(points);

  } else {

    //Use a dummy file that stores x,y,z of the center of each SM
    //useful for testing...
    std::ifstream inputFile(txtFileName.Data());
    if (!inputFile) {
      AliError(("Cannot open the survey file " + txtFileName).Data());
      return;
    }
    Int_t dummyInt = 0;
    Double_t *xReal = new Double_t[fNSuperModule];
    Double_t *yReal = new Double_t[fNSuperModule];
    Double_t *zReal = new Double_t[fNSuperModule];
    Double_t *psiReal = new Double_t[fNSuperModule];
    Double_t *thetaReal = new Double_t[fNSuperModule];
    Double_t *phiReal = new Double_t[fNSuperModule];

    for (Int_t i = 0; i < fNSuperModule; ++i) {
      if (!inputFile) {
	AliError("Error while reading input file.");
	delete [] xReal;
	delete [] yReal;
	delete [] zReal;
	delete [] psiReal;
	delete [] thetaReal;
	delete [] phiReal;
	return;
      }
      inputFile>>dummyInt>>xReal[i]>>yReal[i]>>zReal[i]>>psiReal[i]>>thetaReal[i]>>phiReal[i];
    }

    InitSuperModuleData(xReal, yReal, zReal, psiReal, thetaReal, phiReal);

    delete [] xReal;
    delete [] yReal;
    delete [] zReal;
    delete [] psiReal;
    delete [] thetaReal;
    delete [] phiReal;

  } //kDummy way of doing it

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

    ///////////////////////////////
    // JLK 13-July-2010
    //
    // VERY IMPORTANT!!!!
    //
    // All numbers were calculated in ALICE global c.s., which means
    // that the last argument in the creation of AliAlignObjParams
    // MUST BE set to true
    //////////////////////////////
    new(array[arrayInd])
      AliAlignObjParams(
			smodName.Data(), volid, 
			t.fXShift, t.fYShift, t.fZShift, 
			-t.fPsi, -t.fTheta, -t.fPhi, 
			true
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
  //This method uses the data points from the EMCAL survey and CATIA program to
  //create the alignment matrices.  Only valid for (installed)
  //SM, others will have null objects

  /*--------------------------------------
    The bottom edges of the strip modules
    define the active area of the EMCAL, but
    in software we have a box to hold them which
    is longer than that.  We need to convert
    info about the position of the corners of the
    bottom of the active area to the center of
    the software box that contains the strip
    modules.

    View from beam axis up to EMCAL
              Ai                Ci

        0,1         0,0 1,0          1,1
    xxxxxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxxxxxx
    x   x             x x              x   x
    x   x    % *      x x      * %     x   x
    x   x             x x              x   x
    xxxxxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxxxxxx
        1,1         1,0 0,0          0,1
    <--> = added length                 <--> = added length

    * represents the center of the active area
    % represents the center of the full box (with added length)

         View from side of topmost SM

              Ai                Ci

        0,1         0,0 1,0          1,1
    xxxxxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxxxxxx
    x   x    % *      x x      % *     x   x
    xxxxxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxxxxxx
        1,1         1,0 0,0          0,1
    <--> = added length                 <--> = added length

    * represents the center of the active area
    % represents the center of the full box (with added length)

    -------------------------------------*/

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

    //printf("PHI OF IDEAL SM = %.2f\n",smc.fPhi);

  }

  //Real coordinates of center and rotation angles need to be computed
  //from the survey/CATIA points

  char substr[100];
  AliEMCALSuperModuleCoords *realSM = new AliEMCALSuperModuleCoords[fNSuperModule];
  for (Int_t smodnum = 0; smodnum < geom->GetNumberOfSuperModules(); ++smodnum) {
    AliEMCALSuperModuleCoords &smc = realSM[smodnum];
    Double_t zLength = pars[2]*2.; //length of SM in z from software
    Double_t halfHeight = pars[0]; //half the height of the SM in y direction

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
	//Note: order of values is 00, 01, 10, 11
	xval.push_back(pt->GetX()*100.); //convert m to cm
	yval.push_back(pt->GetY()*100.); 
	zval.push_back(pt->GetZ()*100.); 
      }
    }

    //compute center of active area of each SM on bottome face from survey points
    Double_t activeX = ((xval[0] + (xval[2] - xval[0])/2.)        //x00 and x10
			+(xval[1] + (xval[3] - xval[1])/2.) ) /2.; //x01 and x11
    
//    Double_t activeY = ((yval[0] + (yval[2] - yval[0])/2.)
//			+(yval[1] + (yval[3] - yval[1])/2.) ) /2.;
//    
//    Double_t activeZ = ((zval[0] + (zval[2] - zval[0])/2.)
//			+(zval[1] + (zval[3] - zval[1])/2.) ) /2.;
    
    //printf("Bottom Center of active area of SM %s: %.2f, %.2f, %.2f\n",substr,activeX,activeY,activeZ);
    
    //compute angles for each SM
    //rotation about each axis
    //phi = angle in x-y plane
    
    Double_t realphi = 0.;
    //Note: this is phi wrt y axis.  To get phi wrt to x, add pi/2
    if(smodnum%2 == 0) {
      realphi = (TMath::ATan((yval[2] - yval[0])/(xval[2] - xval[0])) 
		 +TMath::ATan((yval[3] - yval[1])/(xval[3] - xval[1])) )/2.;
    } else {
      realphi = (TMath::ATan((yval[0] - yval[2])/(xval[0] - xval[2]))
		 +TMath::ATan((yval[1] - yval[3])/(xval[1] - xval[3])) )/2.;
    }

    //NOTE: Psi angle is always zero because the two z values being
    //subtracted are exactly the same, but just in case that could change...
    //psi = angle in x-z plane
    Double_t realpsi = (TMath::ATan((zval[0] - zval[2])/(xval[2] - xval[0]))
			+TMath::ATan((zval[1] - zval[3])/(xval[3] - xval[1])) )/2.;
    
    //theta = angle in y-z plane
    Double_t realtheta = TMath::Pi()/2. 
                       - (TMath::ATan((zval[2] - zval[3])/(yval[3] - yval[2]))
			  +TMath::ATan((zval[0] - zval[1])/(yval[1] - yval[0])) )/2.;

    //printf("Old edge of %s 01: %.2f, %.2f, %.2f\n",substr,xval[1],yval[1],zval[1]);
    //printf("Old edge of %s 11: %.2f, %.2f, %.2f\n",substr,xval[3],yval[3],zval[3]);

    //Now calculate the center of the box in z with length added to the 01
    //and 11 corners, corrected by the theta angle
    Double_t activeLength = TMath::Abs(((zval[1] - zval[0]) + (zval[3] - zval[2]))/2.);
    //printf("ACTIVE LENGTH = %.2f\n",activeLength);
    if(smodnum%2 == 0) {
      yval[1] += (zLength - activeLength)*sin(realtheta);
      yval[3] += (zLength - activeLength)*sin(realtheta);
      zval[1] += (zLength - activeLength)*cos(realtheta);
      zval[3] += (zLength - activeLength)*cos(realtheta);
    } else {
      yval[1] -= (zLength - activeLength)*sin(realtheta);
      yval[3] -= (zLength - activeLength)*sin(realtheta);
      zval[1] -= (zLength - activeLength)*cos(realtheta);
      zval[3] -= (zLength - activeLength)*cos(realtheta);
    }

    //printf("New extended edge of %s 01: %.2f, %.2f, %.2f\n",substr,xval[1],yval[1],zval[1]);            
    //printf("New extended edge of %s 11: %.2f, %.2f, %.2f\n",substr,xval[3],yval[3],zval[3]);

    //Compute the center of the bottom of the box in x,y,z
    Double_t realX = activeX;    
    Double_t realY = ((yval[0] + (yval[2] - yval[0])/2.)
			+(yval[1] + (yval[3] - yval[1])/2.) ) /2.;    
    Double_t realZ = ((zval[0] + (zval[2] - zval[0])/2.)
			+(zval[1] + (zval[3] - zval[1])/2.) ) /2.;


    //printf("Bottom Center of SM %s Box: %.2f, %.2f, %.2f\n",substr,realX,realY,realZ);

    //correct the SM centers so that we have the center of the box in
    //x,y using the phi,theta angles                   
    realX += halfHeight*TMath::Cos(TMath::Pi()/2+realphi);
    realY += halfHeight*(TMath::Sin(TMath::Pi()/2+realphi) + TMath::Sin(realtheta));
    realZ += halfHeight*TMath::Cos(TMath::Pi()/2-realtheta);

    //printf("Rotation angles of SM %s: %.4f, %.4f, %.4f\n",substr,realphi*TMath::RadToDeg(),realpsi*TMath::RadToDeg(),realtheta*TMath::RadToDeg());
    //printf("Middle of SM %s: %.2f, %.2f, %.2f\n\n",substr,realX,realY,realZ);

    smc.fX1 = realX;
    smc.fY1 = realY;
    smc.fZ1 = realZ;

    smc.fPhi = 90. + realphi*TMath::RadToDeg();
    smc.fTheta = 0. + realtheta*TMath::RadToDeg();
    smc.fPsi = 0. + realpsi*TMath::RadToDeg();

  }//loop over supermodules
  
  fSuperModuleData = new AliEMCALSuperModuleDelta[fNSuperModule];
  
  for (Int_t i = 0; i < fNSuperModule; ++i) {
    const AliEMCALSuperModuleCoords &real = realSM[i];
    const AliEMCALSuperModuleCoords &ideal = idealSM[i];
    AliEMCALSuperModuleDelta &t = fSuperModuleData[i];
    t.fXShift = real.fX1 - ideal.fX1;
    t.fYShift = real.fY1 - ideal.fY1;
    t.fZShift = ideal.fZ1 - real.fZ1; //due to z flip for C side
    if(i%2==0) {
      t.fZShift *= -1.0;  //correct shift for C side
    }
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
void AliEMCALSurvey::InitSuperModuleData(const Double_t *xReal, const Double_t *yReal, 
					 const Double_t *zReal, const Double_t *psiReal,
					 const Double_t *thetaReal, const Double_t *phiReal)
{
  ///////////////////////////////////////
  //Dummy method just takes the inputted values and applies them
  //
  //Useful for testing small changes
  //////////////////////////////////////
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  //Center of supermodules
  Float_t *pars = geom->GetSuperModulesPars();
  Double_t rpos = (geom->GetEnvelop(0) + geom->GetEnvelop(1))/2.;
  Double_t phi, phiRad, xpos, ypos, zpos;

  zpos = pars[2];

  AliEMCALSuperModuleCoords *idealSM = new AliEMCALSuperModuleCoords[fNSuperModule];
  for (Int_t smodnum = 0; smodnum < geom->GetNumberOfSuperModules(); ++smodnum) {
    AliEMCALSuperModuleCoords &smc = idealSM[smodnum];
    phiRad = geom->GetPhiCenterOfSM(smodnum); //comes in radians
    phi = phiRad*180./TMath::Pi(); //need degrees for AliAlignObjParams
    xpos = rpos * TMath::Cos(phiRad);
    ypos = rpos * TMath::Sin(phiRad);
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

  }

  AliEMCALSuperModuleCoords *realSM = new AliEMCALSuperModuleCoords[fNSuperModule];
  for (Int_t smodnum = 0; smodnum < geom->GetNumberOfSuperModules(); ++smodnum) {
    AliEMCALSuperModuleCoords &smc = realSM[smodnum];
    smc.fX1    = xReal[smodnum];  
    smc.fY1    = yReal[smodnum];  
    smc.fZ1    = zReal[smodnum];  
    smc.fTheta = thetaReal[smodnum];
    smc.fPsi   = psiReal[smodnum];
    smc.fPhi   = phiReal[smodnum];
  }
  
  fSuperModuleData = new AliEMCALSuperModuleDelta[fNSuperModule];
  
  for (Int_t i = 0; i < fNSuperModule; ++i) {
    const AliEMCALSuperModuleCoords &real = realSM[i];
    const AliEMCALSuperModuleCoords &ideal = idealSM[i];
    AliEMCALSuperModuleDelta &t = fSuperModuleData[i];
    t.fTheta = real.fTheta - ideal.fTheta;
    t.fPsi = 0.;
    t.fPhi = real.fPhi - ideal.fPhi;
    t.fXShift = real.fX1 - ideal.fX1;
    t.fYShift = real.fY1 - ideal.fY1;
    t.fZShift = real.fZ1 - ideal.fZ1;

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

