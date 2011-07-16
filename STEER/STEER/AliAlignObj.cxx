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

//-----------------------------------------------------------------
//  Implementation of the alignment object class, holding the alignment
//  constants for a single volume, through the abstract class AliAlignObj.
//  From it two derived concrete representation of alignment object class
//  (AliAlignObjParams, AliAlignObjMatrix) are derived in separate files.
//-----------------------------------------------------------------

#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <TGeoOverlap.h>
#include <TMath.h>

#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
 
ClassImp(AliAlignObj)

//_____________________________________________________________________________
AliAlignObj::AliAlignObj():
  TObject(),
  fVolPath(),
  fVolUID(0)
{
  // default constructor
  for(Int_t i=0; i<6; i++) fDiag[i]=-999.;
  for(Int_t i=0; i<15; i++) fODia[i]=-999.;
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const char* symname, UShort_t voluid) :
  TObject(),
  fVolPath(symname),
  fVolUID(voluid)
{
  // standard constructor
  //
  for(Int_t i=0; i<6; i++) fDiag[i]=-999.;
  for(Int_t i=0; i<15; i++) fODia[i]=-999.;
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const char* symname, UShort_t voluid, Double_t* cmat) :
  TObject(),
  fVolPath(symname),
  fVolUID(voluid)
{
  // standard constructor
  //
  SetCorrMatrix(cmat);
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const AliAlignObj& theAlignObj) :
  TObject(theAlignObj),
  fVolPath(theAlignObj.GetSymName()),
  fVolUID(theAlignObj.GetVolUID())
{
  //copy constructor
  for(Int_t i=0; i<6; i++) fDiag[i]=theAlignObj.fDiag[i];
  for(Int_t i=0; i<15; i++) fODia[i]=theAlignObj.fODia[i];
}

//_____________________________________________________________________________
AliAlignObj &AliAlignObj::operator =(const AliAlignObj& theAlignObj)
{
  // assignment operator
  if(this==&theAlignObj) return *this;
  fVolPath = theAlignObj.GetSymName();
  fVolUID = theAlignObj.GetVolUID();
  for(Int_t i=0; i<6; i++) fDiag[i]=theAlignObj.fDiag[i];
  for(Int_t i=0; i<15; i++) fODia[i]=theAlignObj.fODia[i];
  return *this;
}

//_____________________________________________________________________________
AliAlignObj &AliAlignObj::operator*=(const AliAlignObj& theAlignObj)
{
  // multiplication operator
  // The operator can be used to 'combine'
  // two alignment objects
  TGeoHMatrix m1;
  GetMatrix(m1);
  TGeoHMatrix m2;
  theAlignObj.GetMatrix(m2);
  m1.MultiplyLeft(&m2);
  SetMatrix(m1);
  // temporary solution: the covariance matrix of the resulting combined object
  // is set equal to the covariance matrix of the right operand
  // (not to be used for combining alignment objects for different levels)
  for(Int_t i=0; i<6; i++)  fDiag[i] = theAlignObj.fDiag[i];
  for(Int_t i=0; i<15; i++)  fODia[i] = theAlignObj.fODia[i];  
  return *this;
}

//_____________________________________________________________________________
AliAlignObj::~AliAlignObj()
{
  // dummy destructor
}

//_____________________________________________________________________________
void AliAlignObj::SetVolUID(AliGeomManager::ELayerID detId, Int_t modId)
{
  // From detector name and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for detID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  fVolUID = AliGeomManager::LayerToVolUID(detId,modId);
}

//_____________________________________________________________________________
void AliAlignObj::GetVolUID(AliGeomManager::ELayerID &layerId, Int_t &modId) const
{
  // From the fVolUID, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), sets
  // the argument layerId to the identity of the layer to which that volume
  // belongs and sets the argument modId to the identity of that volume
  // internally to the layer.
  //
  layerId = AliGeomManager::VolUIDToLayer(fVolUID,modId);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetPars(Double_t tr[], Double_t angles[]) const
{
  GetTranslation(tr);
  return GetAngles(angles);
}

//_____________________________________________________________________________
Int_t AliAlignObj::GetLevel() const
{
  // Return the geometry level of the alignable volume to which
  // the alignment object is associated; this is the number of
  // slashes in the corresponding volume path
  //
  if(!gGeoManager){
    AliWarning("gGeoManager doesn't exist or it is still open: unable to return meaningful level value.");
    return (-1);
  }
  const char* symname = GetSymName();
  const char* path;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    path = pne->GetTitle();
  }else{
    path = symname;
  }

  TString pathStr = path;
  if(pathStr[0]!='/') pathStr.Prepend('/');
  return pathStr.CountChar('/');
}

//_____________________________________________________________________________
Int_t AliAlignObj::Compare(const TObject *obj) const
{
  // Compare the levels of two
  // alignment objects
  // Used in the sorting during
  // the application of alignment
  // objects to the geometry
  //
  Int_t level = GetLevel();
  Int_t level2 = ((AliAlignObj *)obj)->GetLevel();
  if (level == level2)
    return 0;
  else
    return ((level > level2) ? 1 : -1);
}

//______________________________________________________________________________
void AliAlignObj::GetCovMatrix(Double_t *cmat) const
{
  // Fills the cmat argument with the coefficients of the external cov matrix (21 elements)
  // calculating them from the correlation matrix data member
  //

  for(Int_t i=0; i<6; ++i) {
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
      cmat[i*(i+1)/2+j] = (fDiag[j] >= 0. && fDiag[i] >= 0.) ? fODia[(i-1)*i/2+j]*fDiag[j]*fDiag[i]: -999.;
    }

    // Diagonal elements
    cmat[i*(i+1)/2+i] = (fDiag[i] >= 0.) ? fDiag[i]*fDiag[i] : -999.;
  }

  return;
}

//______________________________________________________________________________
void AliAlignObj::GetCovMatrix(TMatrixDSym& mcov) const
{
  // Fills the matrix m passed as argument as the covariance matrix calculated
  // from the coefficients of the reduced covariance matrix data members
  //

  for(Int_t i=0; i<6; ++i) {
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
      mcov(j,i) = mcov(i,j) = (fDiag[j] >= 0. && fDiag[i] >= 0.) ? fODia[(i-1)*i/2+j]*fDiag[j]*fDiag[i]: -999.;
    }

    // Diagonal elements
    mcov(i,i) = (fDiag[i] >= 0.) ? fDiag[i]*fDiag[i] : -999.;
  }

}

//______________________________________________________________________________
Bool_t AliAlignObj::GetLocalCovMatrix(TMatrixDSym& lCov) const
{
  // Calculates the covariance matrix (6x6) associated to the six parameters
  // defining the current alignment in the global coordinates system (and sets
  // in the internal data members) from the covariance matrix (6x6) for the six
  // parameters defining the alignment transformation in the local coordinates
  // system, passed as an argument.
  //
  TMatrixD mJ(6,6);// the jacobian of the transformation from local to global parameters
  if(!GetJacobian(mJ)) return kFALSE;
  
  TMatrixDSym gCov(6);
  GetCovMatrix(gCov);
  
  // Compute the local covariance matrix lcov = mJ^T gcov mJ
  TMatrixD gcovJ(gCov,TMatrixD::kMult,mJ);
  TMatrixD lCovM(mJ,TMatrixD::kTransposeMult,gcovJ);
  // To be done: somehow check that lCovM is close enough to be symmetric
  for(Int_t i=0; i<6; i++)
  {
    lCov(i,i) = lCovM(i,i);
    for(Int_t j=i+1; j<6; j++)
    {
      lCov(i,j)=lCovM(i,j);
      lCov(j,i)=lCovM(i,j);
    }
  }
  
  return kTRUE;
  
}

//______________________________________________________________________________
Bool_t AliAlignObj::GetLocalCovMatrix(Double_t *lCov) const
{
  // Calculates the covariance matrix (6x6) associated to the six parameters
  // defining the current alignment in the global coordinates system (and sets
  // in the internal data members) from the covariance matrix (6x6) for the six
  // parameters defining the alignment transformation in the local coordinates
  // system, passed as an argument.
  //
  TMatrixDSym lCovMatrix(6);
  GetLocalCovMatrix(lCovMatrix);
  
  Int_t k=0;
  for(Int_t i=0; i<6; i++)
    for(Int_t j=i; j<6; j++)
    {
       lCov[k++] = lCovMatrix(i,j);
    }
			
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAlignObj::GetJacobian(TMatrixD& mJ) const
{
  // Compute the jacobian J of the transformation of the six local to the six global delta parameters
  //
  // R00 R01 R02 | (R01Rk2 - R02Rk1)Tk  (R02Rk0 - R00Rk2)Tk  (R00Rk1 - R01Rk0)Tk
  // R00 R01 R02 | (R11Rk2 - R12Rk1)Tk  (R12Rk0 - R10Rk2)Tk  (R10Rk1 - R11Rk0)Tk
  // R00 R01 R02 | (R21Rk2 - R22Rk1)Tk  (R22Rk0 - R20Rk2)Tk  (R20Rk1 - R21Rk0)Tk
  //  -  -   -   -   -   -   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  //  0   0   0  |   R11R22 - R12R21      R12R20 - R10R22      R10R21 - R11R20
  //  0   0   0  |   R21R02 - R22R01      R22R00 - R20R02      R20R01 - R21R00
  //  0   0   0  |   R01R12 - R02R11      R02R10 - R00R12      R00R11 - R01R10
  //
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't compute the global covariance matrix from the local one without an open geometry!");
    return kFALSE;
  }

  const char* symname = GetSymName();
  TGeoPhysicalNode* node;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    if(!pne->GetPhysicalNode()){
      node = gGeoManager->MakeAlignablePN(pne);
    }else{
      node = pne->GetPhysicalNode();
    }
  }else{
    AliWarning(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as volume path!",symname));
    node = (TGeoPhysicalNode*) gGeoManager->MakePhysicalNode(symname);
  }

  if (!node) {
    AliError(Form("Volume name or path %s not valid!",symname));
    return kFALSE;
  }

  TGeoHMatrix gm; //global matrix
  gm = *node->GetMatrix();
  Double_t *tr  = gm.GetTranslation();
  Double_t *rot = gm.GetRotationMatrix();
  
  TGeoHMatrix m; // global delta transformation matrix
  GetMatrix(m);
  // We should probably check that it's sufficinetly close to identity
  // if it's not return because the "small angles" approximation cannot hold

  // 3x3 upper left part (global shifts derived w.r.t. local shifts)
  for(Int_t i=0; i<3; i++)
  {
    for(Int_t j=0; j<3; j++)
    {
      mJ(i,j) = rot[i+3*j];
    }
  }
  
  // 3x3 lower left part (global angles derived w.r.t. local shifts)
  for(Int_t i=0; i<3; i++)
  {
    for(Int_t j=0; j<3; j++)
    {
      mJ(i+3,j) = 0.;
    }
  }
  
  // 3x3 upper right part (global shifts derived w.r.t. local angles)
  for(Int_t i=0; i<3; i++)
  {
    for(Int_t j=0; j<3; j++)
    {
      Double_t mEl = 0.;
      Int_t b = (j+1)%3;
      Int_t d = (j+2)%3;
      for(Int_t k=0; k<3; k++)
      {
	mEl += (rot[3*i+b]*rot[3*k+d])*tr[k]-(rot[3*i+d]*rot[3*k+b])*tr[k];
      }
      mJ(i,j+3) = mEl;
    }
  }
  
  // 3x3 lower right part (global angles derived w.r.t. local angles)
  for(Int_t i=0; i<3; i++)
    for(Int_t j=0; j<3; j++)
    {
      Int_t a = (i+1)%3;
      Int_t b = (j+1)%3;
      Int_t c = (i+2)%3;
      Int_t d = (j+2)%3;
      mJ(i+3,j+3) = rot[3*a+b]*rot[3*c+d]-rot[3*a+d]*rot[3*c+b];
    }

  return kTRUE;

}

//______________________________________________________________________________
Bool_t AliAlignObj::SetFromLocalCov(TMatrixDSym& lCov)
{
  // Calculates the covariance matrix (6x6) associated to the six parameters
  // defining the current alignment in the global coordinates system (and sets
  // in the internal data members) from the covariance matrix (6x6) for the six
  // parameters defining the alignment transformation in the local coordinates
  // system, passed as an argument.
  //
  TMatrixD mJ(6,6);// the jacobian of the transformation from local to global parameters
  if(!GetJacobian(mJ)) return kFALSE;
  
  // Compute the global covariance matrix gcov = mJ lcov mJ'
  TMatrixD trJ(TMatrixD::kTransposed, mJ);
  TMatrixD lcovTrJ(lCov,TMatrixD::kMult,trJ);
  TMatrixD gCovM(mJ,TMatrixD::kMult,lcovTrJ);
  // To be done: somehow check that gCovM is close enough to be symmetric
  TMatrixDSym gCov(6);
  for(Int_t i=0; i<6; i++)
  {
    gCov(i,i) = gCovM(i,i);
    for(Int_t j=i+1; j<6; j++)
    {
      gCov(i,j)=gCovM(i,j);
      gCov(j,i)=gCovM(i,j);
    }
  }
  SetCorrMatrix(gCov);

  return kTRUE;
  
}

//______________________________________________________________________________
Bool_t AliAlignObj::SetFromLocalCov(Double_t *lCov)
{
  // Calculates the covariance matrix (6x6) associated to the six parameters
  // defining the current alignment in the global coordinates system, and sets
  // in the internal data members, from the 21 coefficients, passed as argument,
  // of the covariance matrix (6x6) for the six parameters defining the
  // alignment transformation in the local coordinates system.
  //
  TMatrixDSym lCovMatrix(6);
  
  Int_t k=0;
  for(Int_t i=0; i<6; i++)
    for(Int_t j=i; j<6; j++)
    {
      lCovMatrix(i,j) = lCov[k++];
      if(j!=i) lCovMatrix(j,i) = lCovMatrix(i,j);
    }
			
  return SetFromLocalCov(lCovMatrix);

}
  

//______________________________________________________________________________
void AliAlignObj::SetCorrMatrix(Double_t *cmat)
{
  // Sets the correlation matrix data member from the coefficients of the external covariance
  // matrix (21 elements passed as argument). 
  //
  if(cmat) {

    // Diagonal elements first
    for(Int_t i=0; i<6; ++i) {
      fDiag[i] = (cmat[i*(i+1)/2+i] >= 0.) ? TMath::Sqrt(cmat[i*(i+1)/2+i]) : -999.;
    }

    // ... then the ones off diagonal
    for(Int_t i=0; i<6; ++i)
      // Off diagonal elements
      for(Int_t j=0; j<i; ++j) {
	fODia[(i-1)*i/2+j] = (fDiag[i] > 0. && fDiag[j] > 0.) ? cmat[i*(i+1)/2+j]/(fDiag[j]*fDiag[i]) : 0.;       // check for division by zero (due to diagonal element of 0) and for fDiag != -999. (due to negative input diagonal element).
	if (fODia[(i-1)*i/2+j]>1.)  fODia[(i-1)*i/2+j] =  1.; // check upper boundary
	if (fODia[(i-1)*i/2+j]<-1.) fODia[(i-1)*i/2+j] = -1.; // check lower boundary
      }
  } else {
    for(Int_t i=0; i< 6; ++i) fDiag[i]=-999.;
    for(Int_t i=0; i< 6*(6-1)/2; ++i) fODia[i]=0.;
  }

  return;
}

//______________________________________________________________________________
void AliAlignObj::SetCorrMatrix(TMatrixDSym& mcov)
{
  // Sets the correlation matrix data member from the covariance matrix mcov passed
  // passed as argument. 
  //
  if(mcov.IsValid()) {

    // Diagonal elements first
    for(Int_t i=0; i<6; ++i) {
      fDiag[i] = (mcov(i,i) >= 0.) ? TMath::Sqrt(mcov(i,i)) : -999.;
    }

    // ... then the ones off diagonal
    for(Int_t i=0; i<6; ++i)
      // Off diagonal elements
      for(Int_t j=0; j<i; ++j) {
	fODia[(i-1)*i/2+j] = (fDiag[i] > 0. && fDiag[j] > 0.) ? mcov(i,j)/(fDiag[j]*fDiag[i]) : 0.;       // check for division by zero (due to diagonal element of 0) and for fDiag != -999. (due to negative input diagonal element).
	if (fODia[(i-1)*i/2+j]>1.)  fODia[(i-1)*i/2+j] =  1.; // check upper boundary
	if (fODia[(i-1)*i/2+j]<-1.) fODia[(i-1)*i/2+j] = -1.; // check lower boundary
      }
  } else {
    for(Int_t i=0; i< 6; ++i) fDiag[i]=-999.;
    for(Int_t i=0; i< 6*(6-1)/2; ++i) fODia[i]=0.;
  }

  return;
}

//_____________________________________________________________________________
void AliAlignObj::AnglesToMatrix(const Double_t *angles, Double_t *rot) const
{
  // Calculates the rotation matrix using the 
  // Euler angles in "x y z" notation
  //
  Double_t degrad = TMath::DegToRad();
  Double_t sinpsi = TMath::Sin(degrad*angles[0]);
  Double_t cospsi = TMath::Cos(degrad*angles[0]);
  Double_t sinthe = TMath::Sin(degrad*angles[1]);
  Double_t costhe = TMath::Cos(degrad*angles[1]);
  Double_t sinphi = TMath::Sin(degrad*angles[2]);
  Double_t cosphi = TMath::Cos(degrad*angles[2]);

  rot[0] =  costhe*cosphi;
  rot[1] = -costhe*sinphi;
  rot[2] =  sinthe;
  rot[3] =  sinpsi*sinthe*cosphi + cospsi*sinphi;
  rot[4] = -sinpsi*sinthe*sinphi + cospsi*cosphi;
  rot[5] = -costhe*sinpsi;
  rot[6] = -cospsi*sinthe*cosphi + sinpsi*sinphi;
  rot[7] =  cospsi*sinthe*sinphi + sinpsi*cosphi;
  rot[8] =  costhe*cospsi;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::MatrixToAngles(const Double_t *rot, Double_t *angles) const
{
  // Calculates the Euler angles in "x y z" notation
  // using the rotation matrix
  // Returns false in case the rotation angles can not be
  // extracted from the matrix
  //
  if(TMath::Abs(rot[0])<1e-7 || TMath::Abs(rot[8])<1e-7) {
    AliError("Failed to extract roll-pitch-yall angles!");
    return kFALSE;
  }
  Double_t raddeg = TMath::RadToDeg();
  angles[0]=raddeg*TMath::ATan2(-rot[5],rot[8]);
  angles[1]=raddeg*TMath::ASin(rot[2]);
  angles[2]=raddeg*TMath::ATan2(-rot[1],rot[0]);
  return kTRUE;
}

//______________________________________________________________________________
void AliAlignObj::Transform(AliTrackPoint &p, Bool_t copycov) const
{
  // The method transforms the space-point coordinates using the
  // transformation matrix provided by the AliAlignObj
  // In case the copycov flag is set to kTRUE, the covariance matrix 
  // of the alignment object is copied into the space-point
  //
  if (fVolUID != p.GetVolumeID())
    AliWarning(Form("Alignment object ID is not equal to the space-point ID (%d != %d)",fVolUID,p.GetVolumeID())); 

  TGeoHMatrix m;
  GetMatrix(m);
  Double_t *rot = m.GetRotationMatrix();
  Double_t *tr  = m.GetTranslation();

  Float_t xyzin[3],xyzout[3];
  p.GetXYZ(xyzin);
  for (Int_t i = 0; i < 3; i++)
    xyzout[i] = tr[i]+
                xyzin[0]*rot[3*i]+
                xyzin[1]*rot[3*i+1]+
                xyzin[2]*rot[3*i+2];
  p.SetXYZ(xyzout);

  if(copycov){
    TMatrixDSym covmat(6);
    GetCovMatrix(covmat); 
    p.SetAlignCovMatrix(covmat);
  }
  
}

//_____________________________________________________________________________
void AliAlignObj::Transform(AliTrackPointArray &array) const
{
  // This method is used to transform all the track points
  // from the input AliTrackPointArray
  // 
  AliTrackPoint p;
  for (Int_t i = 0; i < array.GetNPoints(); i++) {
    array.GetPoint(p,i);
    Transform(p);
    array.AddPoint(i,&p);
  }
}

//_____________________________________________________________________________
void AliAlignObj::Print(Option_t *) const
{
  // Print the contents of the
  // alignment object in angles and
  // matrix representations
  //
  Double_t tr[3];
  GetTranslation(tr);
  Double_t angles[3];
  GetAngles(angles);
  TGeoHMatrix m;
  GetMatrix(m);
  const Double_t *rot = m.GetRotationMatrix();

  printf("Volume=%s\n",GetSymName());
  if (GetVolUID() != 0) {
    AliGeomManager::ELayerID layerId;
    Int_t modId;
    GetVolUID(layerId,modId);
    printf("VolumeID=%d LayerID=%d ( %s ) ModuleID=%d\n", GetVolUID(),layerId,AliGeomManager::LayerName(layerId),modId);
  }
  printf("%12.8f%12.8f%12.8f    Tx = %12.8f    Psi   = %12.8f\n", rot[0], rot[1], rot[2], tr[0], angles[0]);
  printf("%12.8f%12.8f%12.8f    Ty = %12.8f    Theta = %12.8f\n", rot[3], rot[4], rot[5], tr[1], angles[1]);
  printf("%12.8f%12.8f%12.8f    Tz = %12.8f    Phi   = %12.8f\n", rot[6], rot[7], rot[8], tr[2], angles[2]);

}

//_____________________________________________________________________________
void AliAlignObj::SetPars(Double_t x, Double_t y, Double_t z,
			  Double_t psi, Double_t theta, Double_t phi)
{
  // Set the global delta transformation by passing 3 angles (expressed in
  // degrees) and 3 shifts (in centimeters)
  // 
  SetTranslation(x,y,z);
  SetRotation(psi,theta,phi);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalPars(Double_t x, Double_t y, Double_t z,
				 Double_t psi, Double_t theta, Double_t phi)
{
  // Set the global delta transformation by passing the parameters
  // for the local delta transformation (3 shifts and 3 angles).
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix m;
  Double_t tr[3] = {x, y, z};
  m.SetTranslation(tr);
  Double_t angles[3] = {psi, theta, phi};
  Double_t rot[9];
  AnglesToMatrix(angles,rot);
  m.SetRotation(rot);

  return SetLocalMatrix(m);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalTranslation(Double_t x, Double_t y, Double_t z)
{
  // Set the global delta transformation by passing the three shifts giving
  // the translation in the local reference system of the alignable
  // volume (known by TGeo geometry).
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix m;
  Double_t tr[3] = {x, y, z};
  m.SetTranslation(tr);

  return SetLocalMatrix(m);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalTranslation(const TGeoMatrix& m)
{
  // Set the global delta transformation by passing the matrix of
  // the local delta transformation and taking its translational part
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  const Double_t* tr = m.GetTranslation();
  TGeoHMatrix mtr;
  mtr.SetTranslation(tr);

  return SetLocalMatrix(mtr);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalRotation(Double_t psi, Double_t theta, Double_t phi)
{
  // Set the global delta transformation by passing the three angles giving
  // the rotation in the local reference system of the alignable
  // volume (known by TGeo geometry).
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix m;
  Double_t angles[3] = {psi, theta, phi};
  Double_t rot[9];
  AnglesToMatrix(angles,rot);
  m.SetRotation(rot);

  return SetLocalMatrix(m);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalRotation(const TGeoMatrix& m)
{
  // Set the global delta transformation by passing the matrix of
  // the local delta transformation and taking its rotational part
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix rotm;
  const Double_t* rot = m.GetRotationMatrix();
  rotm.SetRotation(rot);

  return SetLocalMatrix(rotm);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalMatrix(const TGeoMatrix& m)
{
  // Set the global delta transformation by passing the TGeo matrix
  // for the local delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't set the local alignment object parameters! gGeoManager doesn't exist or it is still open!");
    return kFALSE;
  }

  const char* symname = GetSymName();
  TGeoHMatrix gprime,gprimeinv;
  TGeoPhysicalNode* pn = 0;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne)
  {
    pn = pne->GetPhysicalNode();
    if(pn){
      if (pn->IsAligned())
	AliWarning(Form("Volume %s has been misaligned already!",symname));
      gprime = *pn->GetMatrix();
    }else{
      gprime = pne->GetGlobalOrig();
    }
  }else{
    AliWarning(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as volume path!",symname));
    if(!gGeoManager->cd(symname)) {
      AliError(Form("Volume name or path %s not valid!",symname));
      return kFALSE;
    }
    gprime = *gGeoManager->GetCurrentMatrix();
  }

  TGeoHMatrix m1; // the TGeoHMatrix copy of the local delta "m"
  const Double_t *tr = m.GetTranslation();
  m1.SetTranslation(tr);
  const Double_t* rot = m.GetRotationMatrix();
  m1.SetRotation(rot);

  gprimeinv = gprime.Inverse();
  m1.Multiply(&gprimeinv);
  m1.MultiplyLeft(&gprime);

  return SetMatrix(m1);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetMatrix(const TGeoMatrix& m)
{
  // Set the global delta transformation by passing the TGeoMatrix
  // for it
  //
  SetTranslation(m);
  return SetRotation(m);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetLocalPars(Double_t transl[], Double_t angles[]) const
{
  // Get the translations and angles (in degrees) expressing the
  // local delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  if(!GetLocalTranslation(transl)) return kFALSE;
  return GetLocalAngles(angles);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetLocalTranslation(Double_t* tr) const
{
  // Get the 3 shifts giving the translational part of the local
  // delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix ml;
  if(!GetLocalMatrix(ml)) return kFALSE;
  const Double_t* transl;
  transl = ml.GetTranslation();
  tr[0]=transl[0];
  tr[1]=transl[1];
  tr[2]=transl[2];
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetLocalAngles(Double_t* angles) const
{
  // Get the 3 angles giving the rotational part of the local
  // delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix ml;
  if(!GetLocalMatrix(ml)) return kFALSE;
  const Double_t *rot = ml.GetRotationMatrix();
  return MatrixToAngles(rot,angles);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetLocalMatrix(TGeoHMatrix& m) const
{
  // Get the matrix for the local delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't get the local alignment object parameters! gGeoManager doesn't exist or it is still open!");
    return kFALSE;
  }

  const char* symname = GetSymName();
  TGeoPhysicalNode* node;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    if(!pne->GetPhysicalNode()){
      node = gGeoManager->MakeAlignablePN(pne);
    }else{
      node = pne->GetPhysicalNode();
    }
  }else{
    AliWarning(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as volume path!",symname));
    node = (TGeoPhysicalNode*) gGeoManager->MakePhysicalNode(symname);
  }

  if (!node) {
    AliError(Form("Volume name or path %s not valid!",symname));
    return kFALSE;
  }
//  if (node->IsAligned())
//    AliWarning(Form("Volume %s has been misaligned already!",symname));

  GetMatrix(m);
  TGeoHMatrix gprime,gprimeinv;
  gprime = *node->GetMatrix();
  gprimeinv = gprime.Inverse();
  m.Multiply(&gprime);
  m.MultiplyLeft(&gprimeinv);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::ApplyToGeometry(Bool_t ovlpcheck)
{
  // Apply the current alignment object to the TGeo geometry
  // This method returns FALSE if the symname of the object was not
  // valid neither to get a TGeoPEntry nor as a volume path
  //
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't apply the alignment object! gGeoManager doesn't exist or it is still open!");
    return kFALSE;
  }

  if (gGeoManager->IsLocked()){
    AliError("Can't apply the alignment object! Geometry is locked!");
    return kFALSE;
  }

  const char* symname = GetSymName();
  const char* path;
  TGeoPhysicalNode* node;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    path = pne->GetTitle();
    node = gGeoManager->MakeAlignablePN(pne);
  }else{
    AliDebug(1,Form("The symbolic volume name %s does not correspond to a physical entry. Using it as a volume path!",symname));
    path=symname;
    if (!gGeoManager->CheckPath(path)) {
      AliDebug(1,Form("Volume path %s not valid!",path));
      return kFALSE;
    }
    if (gGeoManager->GetListOfPhysicalNodes()->FindObject(path)) {
      AliError(Form("Volume %s has been misaligned already!",path));
      return kFALSE;
    }
    node = (TGeoPhysicalNode*) gGeoManager->MakePhysicalNode(path);
  }

  if (!node) {
    AliError(Form("Volume path %s not valid!",path));
    return kFALSE;
  }

  //  Double_t threshold = 0.001;
  
  TGeoHMatrix align,gprime;
  gprime = *node->GetMatrix();
  GetMatrix(align);
  gprime.MultiplyLeft(&align);
  TGeoHMatrix *ginv = new TGeoHMatrix;
  TGeoHMatrix *g = node->GetMatrix(node->GetLevel()-1);
  *ginv = g->Inverse();
  *ginv *= gprime;
  AliGeomManager::ELayerID layerId; // unique identity for layer in the alobj
  Int_t modId; // unique identity for volume inside layer in the alobj
  GetVolUID(layerId, modId);
  AliDebug(2,Form("Aligning volume %s of detector layer %d with local ID %d",symname,layerId,modId));
  if(ovlpcheck){
    node->Align(ginv,0,kTRUE); //(trunk of root takes threshold as additional argument)
  }else{
    node->Align(ginv,0,kFALSE);
  }
  if(ovlpcheck)
  {
    TObjArray* ovlpArray =  gGeoManager->GetListOfOverlaps();
    Int_t nOvlp = ovlpArray->GetEntriesFast();
    if(nOvlp)
    {
      AliInfo(Form("Misalignment of node %s generated the following overlaps/extrusions:",node->GetName()));
      for(Int_t i=0; i<nOvlp; i++)
	((TGeoOverlap*)ovlpArray->UncheckedAt(i))->PrintInfo();
    }
  }
      
  return kTRUE;
}


