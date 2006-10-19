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
//__________________________________________________________________
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
/// \authors Bruce Becker, Javier Castillo

#include "AliMUONGeometryMisAligner.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"

#include "AliMpExMap.h"

#include "AliLog.h"

#include <TGeoMatrix.h>
#include <TMath.h>
#include <TRandom.h>

/// \cond CLASSIMP
ClassImp(AliMUONGeometryMisAligner)
/// \endcond

//______________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner(Double_t cartXMisAligM, Double_t cartXMisAligW, Double_t cartYMisAligM, Double_t cartYMisAligW, Double_t angMisAligM, Double_t angMisAligW)
  : TObject(),
    fUseUni(kFALSE),
    fUseGaus(kTRUE),
    fCartXMisAligM(cartXMisAligM),
    fCartXMisAligW(cartXMisAligW), // 0.5 mm. Perhaps this should go into AliMUONConstants.h ? 
    fCartYMisAligM(cartYMisAligM),
    fCartYMisAligW(cartYMisAligW), // 0.5 mm. Perhaps this should go into AliMUONConstants.h ? 
    fAngMisAligM(angMisAligM),
    fAngMisAligW(angMisAligW),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0),
    fDisplacementGenerator(0)
{
  /// Standard constructor
  fDisplacementGenerator = new TRandom(0);
}

//______________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner(Double_t cartMisAligM, Double_t cartMisAligW, Double_t angMisAligM, Double_t angMisAligW)
  : TObject(), 
    fUseUni(kFALSE),
    fUseGaus(kTRUE),
    fCartXMisAligM(cartMisAligM),
    fCartXMisAligW(cartMisAligW), // 0.5 mm. Perhaps this should go into AliMUONConstants.h ? 
    fCartYMisAligM(cartMisAligM),
    fCartYMisAligW(cartMisAligW), // 0.5 mm. Perhaps this should go into AliMUONConstants.h ? 
    fAngMisAligM(angMisAligM),
    fAngMisAligW(angMisAligW),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0),
    fDisplacementGenerator(0)
{
  /// Standard constructor
  fDisplacementGenerator = new TRandom(0);
}

//______________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner(Double_t cartMisAlig, Double_t angMisAlig)
  : TObject(), 
    fUseUni(kTRUE),
    fUseGaus(kFALSE),
    fCartXMisAligM(0.),
    fCartXMisAligW(cartMisAlig), // 0.5 mm. Perhaps this should go into AliMUONConstants.h ? 
    fCartYMisAligM(0.),
    fCartYMisAligW(cartMisAlig), // 0.5 mm. Perhaps this should go into AliMUONConstants.h ? 
    fAngMisAligM(0.),
    fAngMisAligW(angMisAlig),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0),
    fDisplacementGenerator(0)
{
  /// Standard constructor
  fDisplacementGenerator = new TRandom(0);
}

//_____________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner()
  : TObject(), 
    fUseUni(kTRUE),
    fUseGaus(kFALSE),
    fCartXMisAligM(0.),
    fCartXMisAligW(0.),
    fCartYMisAligM(0.),
    fCartYMisAligW(0.),
    fAngMisAligM(0.),
    fAngMisAligW(0.),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0),
    fDisplacementGenerator(0)
{
  /// Default constructor
}

//______________________________________________________________________________
AliMUONGeometryMisAligner::~AliMUONGeometryMisAligner()
{
/// Destructor

  if (fDisplacementGenerator) delete fDisplacementGenerator;
}

//_________________________________________________________________________
void
AliMUONGeometryMisAligner::SetXYAngMisAligFactor(Double_t factor)
{
  /// Set XY angular misalign factor 

  if (TMath::Abs(factor) > 1.0 && factor > 0.)
    fXYAngMisAligFactor = factor;
  else
    AliError(Form("Invalid XY angular misalign factor, %d", factor));
}

//_________________________________________________________________________
void AliMUONGeometryMisAligner::SetZCartMisAligFactor(Double_t factor) 
{
  /// Set XY angular misalign factor 
  if (TMath::Abs(factor)<1.0 && factor>0.)
    fZCartMisAligFactor = factor;
  else
    AliError(Form("Invalid Z cartesian misalign factor, %d", factor));
}

//_________________________________________________________________________
void AliMUONGeometryMisAligner::GetUniMisAlign(Double_t cartMisAlig[3], Double_t angMisAlig[3]) const
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
  cartMisAlig[0] = fDisplacementGenerator->Uniform(-fCartXMisAligW+fCartXMisAligM, fCartXMisAligM+fCartXMisAligW);
  cartMisAlig[1] = fDisplacementGenerator->Uniform(-fCartYMisAligW+fCartYMisAligM, fCartYMisAligM+fCartYMisAligW);
  cartMisAlig[2] = fDisplacementGenerator->Uniform(-fZCartMisAligFactor*(fCartXMisAligW+fCartXMisAligM), fZCartMisAligFactor*(fCartXMisAligM+fCartXMisAligW));  
 
  angMisAlig[0] = fDisplacementGenerator->Uniform(-fXYAngMisAligFactor*(fAngMisAligW+fAngMisAligM), fXYAngMisAligFactor*(fAngMisAligM+fAngMisAligW));
  angMisAlig[1] = fDisplacementGenerator->Uniform(-fXYAngMisAligFactor*(fAngMisAligW+fAngMisAligM), fXYAngMisAligFactor*(fAngMisAligM+fAngMisAligW));
  angMisAlig[2] = fDisplacementGenerator->Uniform(-fAngMisAligW+fAngMisAligM, fAngMisAligM+fAngMisAligW);	// degrees
}

//_________________________________________________________________________
void AliMUONGeometryMisAligner::GetGausMisAlign(Double_t cartMisAlig[3], Double_t angMisAlig[3]) const
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
  cartMisAlig[0] = fDisplacementGenerator->Gaus(fCartXMisAligM, fCartXMisAligW);
  cartMisAlig[1] = fDisplacementGenerator->Gaus(fCartYMisAligM, fCartYMisAligW);
  cartMisAlig[2] = fDisplacementGenerator->Gaus(fCartXMisAligM, fZCartMisAligFactor*fCartXMisAligW);
 
  angMisAlig[0] = fDisplacementGenerator->Gaus(fAngMisAligM, fXYAngMisAligFactor*fAngMisAligW);
  angMisAlig[1] = fDisplacementGenerator->Gaus(fAngMisAligM, fXYAngMisAligFactor*fAngMisAligW);
  angMisAlig[2] = fDisplacementGenerator->Gaus(fAngMisAligM, fAngMisAligW);	// degrees
}

//_________________________________________________________________________
TGeoCombiTrans AliMUONGeometryMisAligner::MisAlign(const TGeoCombiTrans & transform) const
{
  /// Misalign given transformation and return the misaligned transformation
  
  Double_t cartMisAlig[3] = {0,0,0};
  Double_t angMisAlig[3] = {0,0,0};
  const Double_t *trans = transform.GetTranslation();
  TGeoRotation *rot;
  // check if the rotation we obtain is not NULL
  if (transform.GetRotation())
    {
      rot = transform.GetRotation();
    }
  else
    {
      rot = new TGeoRotation("rot");
    }			// default constructor.

  if (fUseUni) { 
    GetUniMisAlign(cartMisAlig,angMisAlig);
  }
  else { 
    if (!fUseGaus) {
      AliWarning("Neither uniform nor gausian distribution is set! Will use gausian...");
    } 
    GetGausMisAlign(cartMisAlig,angMisAlig);
  }

  TGeoTranslation newTrans(cartMisAlig[0] + trans[0], cartMisAlig[1] + trans[1], cartMisAlig[2] + trans[2]);
  
  AliInfo(Form("Rotated by %f about Z axis.", angMisAlig[2]));
  rot->RotateX(angMisAlig[0]);
  rot->RotateY(angMisAlig[1]);
  rot->RotateZ(angMisAlig[2]);

  return TGeoCombiTrans(newTrans, *rot);
}

//______________________________________________________________________
AliMUONGeometryTransformer *
AliMUONGeometryMisAligner::MisAlign(const AliMUONGeometryTransformer *
				 transformer, Bool_t verbose)
{
  /// Takes the internal geometry module transformers, copies them
  /// and gets the Detection Elements from them.
  /// Calculates misalignment parameters and applies these
  /// to the local transform of the Detection Element
  /// Obtains the global transform by multiplying the module transformer
  /// transformation with the local transformation 
  /// Applies the global transform to a new detection element
  /// Adds the new detection element to a new module transformer
  /// Adds the new module transformer to a new geometry transformer
  /// Returns the new geometry transformer


  AliMUONGeometryTransformer *newGeometryTransformer =
    new AliMUONGeometryTransformer(kTRUE);
  for (Int_t iMt = 0; iMt < transformer->GetNofModuleTransformers(); iMt++)
    {				// module transformers
      
      const AliMUONGeometryModuleTransformer *kModuleTransformer =
	transformer->GetModuleTransformer(iMt, true);
      
      AliMUONGeometryModuleTransformer *newModuleTransformer =
	new AliMUONGeometryModuleTransformer(iMt);
      newGeometryTransformer->AddModuleTransformer(newModuleTransformer);

      TGeoCombiTrans moduleTransform =
	TGeoCombiTrans(*kModuleTransformer->GetTransformation());
      TGeoCombiTrans *newModuleTransform = new TGeoCombiTrans(moduleTransform);	
              // same module transform as the previous one 
	      // no mis align object created
      newModuleTransformer->SetTransformation(moduleTransform);

      AliMpExMap *detElements = kModuleTransformer->GetDetElementStore();

      if (verbose)
	AliInfo(Form
		("%i DEs in old GeometryStore  %i",
		 detElements->GetSize(), iMt));

      for (Int_t iDe = 0; iDe < detElements->GetSize(); iDe++)
	{			// detection elements.
	  AliMUONGeometryDetElement *detElement =
	    (AliMUONGeometryDetElement *) detElements->GetObject(iDe);
	  if (!detElement)
	    AliFatal("Detection element not found.");

	  /// make a new detection element
	  AliMUONGeometryDetElement *newDetElement =
	    new AliMUONGeometryDetElement(detElement->GetId(),
					  detElement->GetVolumePath());

	  // local transformation of this detection element.
          TGeoCombiTrans localTransform
	    = TGeoCombiTrans(*detElement->GetLocalTransformation());
	  TGeoCombiTrans newLocalTransform = MisAlign(localTransform);
          newDetElement->SetLocalTransformation(newLocalTransform);					  

	  // global transformation
	  TGeoHMatrix newGlobalTransform =
	    AliMUONGeometryBuilder::Multiply(*newModuleTransform,
	                                      newLocalTransform);
	  newDetElement->SetGlobalTransformation(newGlobalTransform);
	  
	  // add this det element to module
	  newModuleTransformer->GetDetElementStore()->Add(newDetElement->GetId(),
							  newDetElement);
          // Get delta transformation: 
	  // Tdelta = Tnew * Told.inverse
	  TGeoHMatrix  deltaTransform
	    = AliMUONGeometryBuilder::Multiply(
	        newGlobalTransform, 
		detElement->GetGlobalTransformation()->Inverse());
	  
	  // Create mis alignment matrix
	  newGeometryTransformer
	    ->AddMisAlignDetElement(detElement->GetId(), deltaTransform);
	}
      if (verbose)
	AliInfo(Form("Added module transformer %i to the transformer", iMt));
      newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
    }
  return newGeometryTransformer;
}




