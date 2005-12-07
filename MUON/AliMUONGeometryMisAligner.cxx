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
// Class AliMUONGeometryMisAligner 
// ----------------------------
// Class for misalignment of geometry transformations
//
// Author:Bruce Becker

//__________________________________________________________________
//
/////////////////////////////////////////////////////////////////////
//This performs the misalignment on an existing muon arm geometry
//  based on the standard definition of the detector elements in 
//  $ALICE_ROOT/MUON/data
//
//  --> User has to specify the magnitude of the alignments, in the Cartesian 
//  co-ordiantes (which are used to apply translation misalignments) and in the
//  spherical co-ordinates (which are used to apply angular displacements)
//  --> If the constructor is used with no arguments, user has to set 
//  misalignment ranges by hand using the methods : 
//  SetApplyMisAlig, SetMaxCartMisAlig, SetMaxAngMisAlig, SetXYAngMisAligFactor
//  (last method takes account of the fact that the misalingment is greatest in 
//  the XY plane, since the detection elements are fixed to a support structure
//  in this plane. Misalignments in the XZ and YZ plane will be very small 
//  compared to those in the XY plane, which are small already - of the order 
//  of microns)

//  Note : If the detection elements are allowed to be misaligned in all
//  directions, this has consequences for the alignment algorithm
//  (AliMUONAlignment), which needs to know the number of free parameters. 
//  Eric only allowed 3 :  x,y,theta_xy, but in principle z and the other 
//  two angles are alignable as well.



#include "AliMUONGeometryMisAligner.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONGeometryBuilder.h"
#include <TGeoManager.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>

#include "AliLog.h"

#include <sstream>

ClassImp(AliMUONGeometryMisAligner)
//______________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner(Double_t cartMisAlig, Double_t angMisAlig)
:TObject(), fDisplacementGenerator(0)
{
  /// Standard constructor
  fMaxCartMisAlig = cartMisAlig;	// 0.5 mm. Perhaps this should go into AliMUONConstants.h ? 
  fMaxAngMisAlig = angMisAlig;
  fXYAngMisAligFactor = 1.0;
  fDisplacementGenerator = new TRandom(0);
}

//_____________________________________________________________________________
AliMUONGeometryMisAligner::AliMUONGeometryMisAligner()
:TObject(), fDisplacementGenerator(0)
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONGeometryMisAligner::
AliMUONGeometryMisAligner(const AliMUONGeometryMisAligner & right):
TObject(right)
{
  /// Copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//______________________________________________________________________________
AliMUONGeometryMisAligner::~AliMUONGeometryMisAligner()
{
/// Destructor

  delete fDisplacementGenerator;
}

//______________________________________________________________________________
AliMUONGeometryMisAligner & AliMUONGeometryMisAligner::
operator=(const AliMUONGeometryMisAligner & right)
{
  /// Assignement operator (not implemented)

  // check assignement to self
  if (this == &right)
    return *this;

  AliFatal("Assignement operator not provided.");

  return *this;
}

void
AliMUONGeometryMisAligner::SetXYAngMisAligFactor(Double_t factor)
{
  if (TMath::Abs(factor) > 1.0 && factor > 0.)
    fXYAngMisAligFactor = factor;
  else
    AliError(Form("Invalid factor, %d", factor));
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
  
  cartMisAlig[0] = fDisplacementGenerator->Uniform(-1. * fMaxCartMisAlig, fMaxCartMisAlig);
  cartMisAlig[1] = fDisplacementGenerator->Uniform(-1. * fMaxCartMisAlig, fMaxCartMisAlig);
  cartMisAlig[2] = fDisplacementGenerator->Uniform(-1. * fMaxCartMisAlig, fMaxCartMisAlig);
  
  TGeoTranslation newTrans(cartMisAlig[0] + trans[0], cartMisAlig[1] + trans[1], cartMisAlig[2] + trans[2]);
  
  /*
    misalign the centre of the local transformation
    rotation axes : 
    fAngMisAlig[1,2,3] = [x,y,z]
    Assume that misalignment about the x and y axes (misalignment of z plane)
    is much smaller, since the entire detection plane has to be moved (the 
    detection elements are on a support structure), while rotation of the x-y
    plane is more free.
  */
  
  angMisAlig[0] = fDisplacementGenerator->Uniform(fXYAngMisAligFactor * fMaxAngMisAlig,  1.0 * fMaxAngMisAlig);
  angMisAlig[1] =    fDisplacementGenerator->Uniform(fXYAngMisAligFactor * fMaxAngMisAlig, 1.0 * fMaxAngMisAlig);
  angMisAlig[2] = fDisplacementGenerator->Uniform(-1. * fMaxAngMisAlig, fMaxAngMisAlig);	// degrees
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
  /////////////////////////////////////////////////////////////////////
  //   Takes the internal geometry module transformers, copies them
  // and gets the Detection Elements from them.
  // Calculates misalignment parameters and applies these
  // to the local transform of the Detection Element
  // Obtains the global transform by multiplying the module transformer
  // transformation with the local transformation 
  // Applies the global transform to a new detection element
  // Adds the new detection element to a new module transformer
  // Adds the new module transformer to a new geometry transformer
  // Returns the new geometry transformer


  AliMUONGeometryTransformer *newGeometryTransformer =
    new AliMUONGeometryTransformer(kTRUE);
  for (Int_t iMt = 0; iMt < transformer->GetNofModuleTransformers(); iMt++)
    {				// module transformers
      
      const AliMUONGeometryModuleTransformer *kModuleTransformer =
	transformer->GetModuleTransformer(iMt, true);
      
      AliMUONGeometryModuleTransformer *newModuleTransformer =
	new AliMUONGeometryModuleTransformer(iMt);
      newGeometryTransformer->AddModuleTransformer(newModuleTransformer);

      const TGeoCombiTrans *kModuleTransform =
	kModuleTransformer->GetTransformation();
      TGeoCombiTrans *newModuleTransform = new TGeoCombiTrans(*kModuleTransform);	// same module transform as the previous one 
      newModuleTransformer->SetTransformation(*kModuleTransform);

      AliMUONGeometryStore *detElements =
	kModuleTransformer->GetDetElementStore();

      if (verbose)
	AliInfo(Form
		("%i DEs in old GeometryStore  %i",
		 detElements->GetNofEntries(), iMt));

      for (Int_t iDe = 0; iDe < detElements->GetNofEntries(); iDe++)
	{			// detection elements.
	  AliMUONGeometryDetElement *detElement =
	    (AliMUONGeometryDetElement *) detElements->GetEntry(iDe);
	  if (!detElement)
	    AliFatal("Detection element not found.");

	  // local transformation of this detection element.
	  const TGeoCombiTrans *kLocalTransform =
	    detElement->GetLocalTransformation();
	  TGeoCombiTrans newLocalTransform = MisAlign(*kLocalTransform);

	  /// make a new detection element with this local transform
	  AliMUONGeometryDetElement *newDetElement =
	    new AliMUONGeometryDetElement(detElement->GetId(),
					  detElement->GetAlignedVolume(),
					  newLocalTransform);

	  TGeoHMatrix newGlobalTransform =
	    AliMUONGeometryBuilder::Multiply(newLocalTransform,
					     *newModuleTransform);

	  newDetElement->SetGlobalTransformation(newGlobalTransform);
	  if (verbose)
	    AliInfo("GlobalTransforms:");
	  newModuleTransformer->GetDetElementStore()->Add(
							  newDetElement->GetId(),
							  newDetElement);
	}
      if (verbose)
	AliInfo(Form("Added module transformer %i to the transformer", iMt));
      newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
    }
  return newGeometryTransformer;
}




