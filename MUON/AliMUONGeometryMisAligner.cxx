/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *      SigmaEffect_thetadegrees                                          * 
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
//
//-----------------------------------------------------------------------------
/// \class AliMUONGeometryMisAligner
///
/// This performs the misalignment on an existing muon arm geometry
/// based on the standard definition of the detector elements in 
/// $ALICE_ROOT/MUON/data
///
/// --> User has to specify the magnitude of the alignments, in the Cartesian 
/// co-ordiantes (which are used to apply translation misalignments) and in the
/// spherical co-ordinates (which are used to apply angular displacements)
///
/// --> If the constructor is used with no arguments, user has to set 
/// misalignment ranges by hand using the methods : 
/// SetApplyMisAlig, SetMaxCartMisAlig, SetMaxAngMisAlig, SetXYAngMisAligFactor
/// (last method takes account of the fact that the misalingment is greatest in 
/// the XY plane, since the detection elements are fixed to a support structure
/// in this plane. Misalignments in the XZ and YZ plane will be very small 
/// compared to those in the XY plane, which are small already - of the order 
/// of microns)
///
/// Note : If the detection elements are allowed to be misaligned in all
/// directions, this has consequences for the alignment algorithm
/// (AliMUONAlignment), which needs to know the number of free parameters. 
/// Eric only allowed 3 :  x,y,theta_xy, but in principle z and the other 
/// two angles are alignable as well.
///
/// \author Bruce Becker, Javier Castillo
//-----------------------------------------------------------------------------

#include "AliMUONGeometryMisAligner.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"

#include "AliAlignObjMatrix.h"
#include "AliMathBase.h"
#include "AliLog.h"

#include <TClonesArray.h>
#include <TGeoMatrix.h>
#include <TMatrixDSym.h>
#include <TMath.h>
#include <TRandom.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONGeometryMisAligner)
/// \endcond

//______________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner(Double_t cartXMisAligM, Double_t cartXMisAligW, Double_t cartYMisAligM, Double_t cartYMisAligW, Double_t angMisAligM, Double_t angMisAligW)
  : TObject(),
    fUseUni(kFALSE),
    fUseGaus(kTRUE),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (Int_t i=0; i<6; i++){
    for (Int_t j=0; j<2; j++){
      fDetElemMisAlig[i][j] = 0.0;
      fModuleMisAlig[i][j] = 0.0;
    }
  }
  fDetElemMisAlig[0][0] = cartXMisAligM; 
  fDetElemMisAlig[0][1] = cartXMisAligW; 
  fDetElemMisAlig[1][0] = cartYMisAligM; 
  fDetElemMisAlig[1][1] = cartYMisAligW; 
  fDetElemMisAlig[5][0] = angMisAligM; 
  fDetElemMisAlig[5][1] = angMisAligW;

}

//______________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner(Double_t cartMisAligM, Double_t cartMisAligW, Double_t angMisAligM, Double_t angMisAligW)
  : TObject(), 
    fUseUni(kFALSE),
    fUseGaus(kTRUE),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (Int_t i=0; i<6; i++){
    for (Int_t j=0; j<2; j++){
      fDetElemMisAlig[i][j] = 0.0;
      fModuleMisAlig[i][j] = 0.0;
    }
  }
  fDetElemMisAlig[0][0] = cartMisAligM; 
  fDetElemMisAlig[0][1] = cartMisAligW; 
  fDetElemMisAlig[1][0] = cartMisAligM; 
  fDetElemMisAlig[1][1] = cartMisAligW; 
  fDetElemMisAlig[5][0] = angMisAligM; 
  fDetElemMisAlig[5][1] = angMisAligW;

}

//______________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner(Double_t cartMisAlig, Double_t angMisAlig)
  : TObject(), 
    fUseUni(kTRUE),
    fUseGaus(kFALSE),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (Int_t i=0; i<6; i++){
    for (Int_t j=0; j<2; j++){
      fDetElemMisAlig[i][j] = 0.0;
      fModuleMisAlig[i][j] = 0.0;
    }
  }
  fDetElemMisAlig[0][1] = cartMisAlig; 
  fDetElemMisAlig[1][1] = cartMisAlig; 
  fDetElemMisAlig[5][1] = angMisAlig;

}

//_____________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner()
  : TObject(), 
    fUseUni(kTRUE),
    fUseGaus(kFALSE),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0)
{
  /// Default constructor
  for (Int_t i=0; i<6; i++){
    for (Int_t j=0; j<2; j++){
      fDetElemMisAlig[i][j] = 0.0;
      fModuleMisAlig[i][j] = 0.0;
    }
  }
}

//______________________________________________________________________________
AliMUONGeometryMisAligner::~AliMUONGeometryMisAligner()
{
/// Destructor

}

//_________________________________________________________________________
void
AliMUONGeometryMisAligner::SetXYAngMisAligFactor(Double_t factor)
{
  /// Set XY angular misalign factor 

  if (TMath::Abs(factor) > 1.0 && factor > 0.){
    fXYAngMisAligFactor = factor;
    fDetElemMisAlig[3][0] = fDetElemMisAlig[5][0]*factor; // These lines were 
    fDetElemMisAlig[3][1] = fDetElemMisAlig[5][1]*factor; // added to keep
    fDetElemMisAlig[4][0] = fDetElemMisAlig[5][0]*factor; // backward 
    fDetElemMisAlig[4][1] = fDetElemMisAlig[5][1]*factor; // compatibility 
  }
  else
    AliError(Form("Invalid XY angular misalign factor, %d", factor));
}

//_________________________________________________________________________
void AliMUONGeometryMisAligner::SetZCartMisAligFactor(Double_t factor) 
{
  /// Set XY angular misalign factor 
  if (TMath::Abs(factor)<1.0 && factor>0.) {
    fZCartMisAligFactor = factor;
    fDetElemMisAlig[2][0] = fDetElemMisAlig[0][0];        // These lines were added to 
    fDetElemMisAlig[2][1] = fDetElemMisAlig[0][1]*factor; // keep backward compatibility
  }
  else
    AliError(Form("Invalid Z cartesian misalign factor, %d", factor));
}

//_________________________________________________________________________
void AliMUONGeometryMisAligner::GetUniMisAlign(Double_t cartMisAlig[3], Double_t angMisAlig[3], const Double_t lParMisAlig[6][2]) const
{
  /// Misalign using uniform distribution
  /**
    misalign the centre of the local transformation
    rotation axes : 
    fAngMisAlig[1,2,3] = [x,y,z]
    Assume that misalignment about the x and y axes (misalignment of z plane)
    is much smaller, since the entire detection plane has to be moved (the 
    detection elements are on a support structure), while rotation of the x-y
    plane is more free.
  */
  cartMisAlig[0] = gRandom->Uniform(-lParMisAlig[0][1]+lParMisAlig[0][0], lParMisAlig[0][0]+lParMisAlig[0][1]);
  cartMisAlig[1] = gRandom->Uniform(-lParMisAlig[1][1]+lParMisAlig[1][0], lParMisAlig[1][0]+lParMisAlig[1][1]);
  cartMisAlig[2] = gRandom->Uniform(-lParMisAlig[2][1]+lParMisAlig[2][0], lParMisAlig[2][0]+lParMisAlig[2][1]);  
 
  angMisAlig[0] = gRandom->Uniform(-lParMisAlig[3][1]+lParMisAlig[3][0], lParMisAlig[3][0]+lParMisAlig[3][1]);
  angMisAlig[1] = gRandom->Uniform(-lParMisAlig[4][1]+lParMisAlig[4][0], lParMisAlig[4][0]+lParMisAlig[4][1]);
  angMisAlig[2] = gRandom->Uniform(-lParMisAlig[5][1]+lParMisAlig[5][0], lParMisAlig[5][0]+lParMisAlig[5][1]);	// degrees
}

//_________________________________________________________________________
void AliMUONGeometryMisAligner::GetGausMisAlign(Double_t cartMisAlig[3], Double_t angMisAlig[3], const Double_t lParMisAlig[6][2]) const
{
  /// Misalign using gaussian distribution
  /**
    misalign the centre of the local transformation
    rotation axes : 
    fAngMisAlig[1,2,3] = [x,y,z]
    Assume that misalignment about the x and y axes (misalignment of z plane)
    is much smaller, since the entire detection plane has to be moved (the 
    detection elements are on a support structure), while rotation of the x-y
    plane is more free.
  */
  cartMisAlig[0] = AliMathBase::TruncatedGaus(lParMisAlig[0][0], lParMisAlig[0][1], 3.*lParMisAlig[0][1]);
  cartMisAlig[1] = AliMathBase::TruncatedGaus(lParMisAlig[1][0], lParMisAlig[1][1], 3.*lParMisAlig[1][1]);
  cartMisAlig[2] = AliMathBase::TruncatedGaus(lParMisAlig[2][0], lParMisAlig[2][1], 3.*lParMisAlig[2][1]);
 
  angMisAlig[0] = AliMathBase::TruncatedGaus(lParMisAlig[3][0], lParMisAlig[3][1], 3.*lParMisAlig[3][1]);
  angMisAlig[1] = AliMathBase::TruncatedGaus(lParMisAlig[4][0], lParMisAlig[4][1], 3.*lParMisAlig[4][1]);
  angMisAlig[2] = AliMathBase::TruncatedGaus(lParMisAlig[5][0], lParMisAlig[5][1], 3.*lParMisAlig[5][1]);	// degrees
}

//_________________________________________________________________________
TGeoCombiTrans AliMUONGeometryMisAligner::MisAlignDetElem(const TGeoCombiTrans & transform) const
{
  /// Misalign given transformation and return the misaligned transformation. 
  /// Use misalignment parameters for detection elements.
  /// Note that applied misalignments are small deltas with respect to the detection 
  /// element own ideal local reference frame. Thus deltaTransf represents 
  /// the transformation to go from the misaligned d.e. local coordinates to the 
  /// ideal d.e. local coordinates. 
  /// Also note that this -is not- what is in the ALICE alignment framework known 
  /// as local nor global (see AliMUONGeometryMisAligner::MisAlign) 
  
  Double_t cartMisAlig[3] = {0,0,0};
  Double_t angMisAlig[3] = {0,0,0};

  if (fUseUni) { 
    GetUniMisAlign(cartMisAlig,angMisAlig,fDetElemMisAlig);
  }
  else { 
    if (!fUseGaus) {
      AliWarning("Neither uniform nor gausian distribution is set! Will use gausian...");
    } 
    GetGausMisAlign(cartMisAlig,angMisAlig,fDetElemMisAlig);
  }

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans,deltaRot);
  TGeoHMatrix newTransfMat = transform * deltaTransf;
    
  AliInfo(Form("Rotated DE by %f about Z axis.", angMisAlig[2]));

  return TGeoCombiTrans(newTransfMat);
}

//_________________________________________________________________________
TGeoCombiTrans AliMUONGeometryMisAligner::MisAlignModule(const TGeoCombiTrans & transform) const
{
  /// Misalign given transformation and return the misaligned transformation. 
  /// Use misalignment parameters for modules.
  /// Note that applied misalignments are small deltas with respect to the module 
  /// own ideal local reference frame. Thus deltaTransf represents 
  /// the transformation to go from the misaligned module local coordinates to the 
  /// ideal module local coordinates. 
  /// Also note that this -is not- what is in the ALICE alignment framework known 
  /// as local nor global (see AliMUONGeometryMisAligner::MisAlign) 
  
  Double_t cartMisAlig[3] = {0,0,0};
  Double_t angMisAlig[3] = {0,0,0};

  if (fUseUni) { 
    GetUniMisAlign(cartMisAlig,angMisAlig,fModuleMisAlig);
  }
  else { 
    if (!fUseGaus) {
      AliWarning("Neither uniform nor gausian distribution is set! Will use gausian...");
    } 
    GetGausMisAlign(cartMisAlig,angMisAlig,fModuleMisAlig);
  }

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans,deltaRot);
  TGeoHMatrix newTransfMat = transform * deltaTransf;

  AliInfo(Form("Rotated Module by %f about Z axis.", angMisAlig[2]));

  return TGeoCombiTrans(newTransfMat);
}

//______________________________________________________________________
AliMUONGeometryTransformer *
AliMUONGeometryMisAligner::MisAlign(const AliMUONGeometryTransformer *
				    transformer, Bool_t verbose)
{
  /// Takes the internal geometry module transformers, copies them to
  /// new geometry module transformers. 
  /// Calculates  module misalignment parameters and applies these
  /// to the new module transformer.
  /// Calculates the module misalignment delta transformation in the 
  /// Alice Alignment Framework newTransf = delta * oldTransf.
  /// Add a module misalignment to the new geometry transformer.
  /// Gets the Detection Elements from the module transformer. 
  /// Calculates misalignment parameters and applies these
  /// to the local transformation of the Detection Element.
  /// Obtains the new global transformation by multiplying the new 
  /// module transformer transformation with the new local transformation. 
  /// Applies the new global transform to a new detection element.
  /// Adds the new detection element to a new module transformer.
  /// Calculates the d.e. misalignment delta transformation in the 
  /// Alice Alignment Framework (newGlobalTransf = delta * oldGlobalTransf).
  /// Add a d.e. misalignment to the new geometry transformer.
  /// Adds the new module transformer to a new geometry transformer.
  /// Returns the new geometry transformer.


  AliMUONGeometryTransformer *newGeometryTransformer =
    new AliMUONGeometryTransformer();
  for (Int_t iMt = 0; iMt < transformer->GetNofModuleTransformers(); iMt++)
    {				// module transformers
      const AliMUONGeometryModuleTransformer *kModuleTransformer =
	transformer->GetModuleTransformer(iMt, true);
      
      AliMUONGeometryModuleTransformer *newModuleTransformer =
	new AliMUONGeometryModuleTransformer(iMt);
      newGeometryTransformer->AddModuleTransformer(newModuleTransformer);

      TGeoCombiTrans moduleTransform =
	TGeoCombiTrans(*kModuleTransformer->GetTransformation());
      // New module transformation
      TGeoCombiTrans newModuleTransform = MisAlignModule(moduleTransform);
      newModuleTransformer->SetTransformation(newModuleTransform);

      // Get delta transformation: 
      // Tdelta = Tnew * Told.inverse
      TGeoHMatrix deltaModuleTransform = 
	AliMUONGeometryBuilder::Multiply(
	  newModuleTransform, 
	  kModuleTransformer->GetTransformation()->Inverse());

      // Create module mis alignment matrix
      newGeometryTransformer
	->AddMisAlignModule(kModuleTransformer->GetModuleId(), deltaModuleTransform);

      AliMpExMap *detElements = kModuleTransformer->GetDetElementStore();

      if (verbose)
	AliInfo(Form("%i DEs in old GeometryStore  %i",detElements->GetSize(), iMt));

      TIter next(detElements->CreateIterator());
      AliMUONGeometryDetElement *detElement;
      
      while ( ( detElement = static_cast<AliMUONGeometryDetElement*>(next()) ) )
      {
	  /// make a new detection element
	  AliMUONGeometryDetElement *newDetElement =
	    new AliMUONGeometryDetElement(detElement->GetId(),
					  detElement->GetVolumePath());

	  // local transformation of this detection element.
          TGeoCombiTrans localTransform
	    = TGeoCombiTrans(*detElement->GetLocalTransformation());
	  TGeoCombiTrans newLocalTransform = MisAlignDetElem(localTransform);
          newDetElement->SetLocalTransformation(newLocalTransform);


	  // global transformation
	  TGeoHMatrix newGlobalTransform =
	    AliMUONGeometryBuilder::Multiply(newModuleTransform,
					     newLocalTransform);
	  newDetElement->SetGlobalTransformation(newGlobalTransform);
	  
	  // add this det element to module
	  newModuleTransformer->GetDetElementStore()->Add(newDetElement->GetId(),
							  newDetElement);

	  // In the Alice Alignment Framework misalignment objects store
	  // global delta transformation
	  // Get detection "intermediate" global transformation
	  TGeoHMatrix newOldGlobalTransform = newModuleTransform * localTransform;
          // Get detection element global delta transformation: 
	  // Tdelta = Tnew * Told.inverse
	  TGeoHMatrix  deltaGlobalTransform
	    = AliMUONGeometryBuilder::Multiply(
	        newGlobalTransform, 
		newOldGlobalTransform.Inverse());
	  
	  // Create mis alignment matrix
	  newGeometryTransformer
	    ->AddMisAlignDetElement(detElement->GetId(), deltaGlobalTransform);
	}

      
      if (verbose)
	AliInfo(Form("Added module transformer %i to the transformer", iMt));
      newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
    }
  return newGeometryTransformer;
}


void AliMUONGeometryMisAligner::SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t rChId, Double_t rChResX, Double_t rChResY, Double_t rDeResX, Double_t rDeResY){

  Int_t chIdMin = (rChId<0)? 0 : rChId;
  Int_t chIdMax = (rChId<0)? 9 : rChId;
  Double_t chResX = (rChResX<0)? fModuleMisAlig[0][1] : rChResX;
  Double_t chResY = (rChResY<0)? fModuleMisAlig[1][1] : rChResY;
  Double_t deResX = (rDeResX<0)? fDetElemMisAlig[0][1] : rDeResX;
  Double_t deResY = (rDeResY<0)? fDetElemMisAlig[1][1] : rDeResY;

  TMatrixDSym mChCorrMatrix(6);
  mChCorrMatrix[0][0]=chResX*chResX;
  mChCorrMatrix[1][1]=chResY*chResY;
  //  mChCorrMatrix.Print();

  TMatrixDSym mDECorrMatrix(6);
  mDECorrMatrix[0][0]=deResX*deResX;
  mDECorrMatrix[1][1]=deResY*deResY;
  //  mDECorrMatrix.Print();

  AliAlignObjMatrix *alignMat = 0x0;

  for(Int_t chId=chIdMin; chId<=chIdMax; chId++) {
    TString chName1;
    TString chName2;
    if (chId<4){
      chName1 = Form("GM%d",chId);
      chName2 = Form("GM%d",chId);
    } else {
      chName1 = Form("GM%d",4+(chId-4)*2);
      chName2 = Form("GM%d",4+(chId-4)*2+1);
    }
    
    for (int i=0; i<misAlignArray->GetEntries(); i++) {
      alignMat = (AliAlignObjMatrix*)misAlignArray->At(i);
      TString volName(alignMat->GetSymName());
      if((volName.Contains(chName1)&&
	  ((volName.Last('/')==volName.Index(chName1)+chName1.Length())||
	   (volName.Length()==volName.Index(chName1)+chName1.Length())))||
	 (volName.Contains(chName2)&&
	  ((volName.Last('/')==volName.Index(chName2)+chName2.Length())||
	   (volName.Length()==volName.Index(chName2)+chName2.Length())))){
	volName.Remove(0,volName.Last('/')+1);
	if (volName.Contains("GM")) {
	  //	alignMat->Print("NULL");
	  alignMat->SetCorrMatrix(mChCorrMatrix);
	} else if (volName.Contains("DE")) {
	  //	alignMat->Print("NULL");
	  alignMat->SetCorrMatrix(mDECorrMatrix);
	}
      }
    }
  }
}



