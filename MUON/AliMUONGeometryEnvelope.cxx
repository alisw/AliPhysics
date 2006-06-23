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

// $Id$
//
// Class AliMUONGeometryEnvelope
// -----------------------------
// Helper class for definititon an assembly of volumes.
// Author: Ivana Hrivnacova, IPN Orsay
// 23/01/2004

#include <TGeoMatrix.h>
#include <TString.h>
#include <TObjArray.h>

#include "AliMUONGeometryEnvelope.h"
#include "AliMUONGeometryConstituent.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONGeometryEnvelope)
/// \endcond

//______________________________________________________________________________
AliMUONGeometryEnvelope::AliMUONGeometryEnvelope(const TString& name, 
                                                 Int_t id, 
                                                 Bool_t isVirtual,
						 const char* only)
 : TNamed(name, name),
   fIsVirtual(isVirtual),
   fIsMANY(false),
   fCopyNo(1),
   fTransformation(0),
   fConstituents(0)
{
/// Standard constructor

  if (TString(only) == TString("MANY")) fIsMANY = true;

  // Create the envelope transformation
  fTransformation = new TGeoCombiTrans("");
  fConstituents = new TObjArray(20);
  
  // Set id
  SetUniqueID(id);
}


//______________________________________________________________________________
AliMUONGeometryEnvelope::AliMUONGeometryEnvelope(const TString& name,   
                                                 Int_t id,
                                                 Int_t copyNo, 
						 const char* only)
 : TNamed(name, name),
   fIsVirtual(false),
   fIsMANY(false),
   fCopyNo(copyNo),
   fTransformation(0),
   fConstituents(0)
{
/// Standard constructor for a non virtual enevelope with a specified copy 
/// number

  if (TString(only) == TString("MANY")) fIsMANY = true;

  // Create the envelope transformation
  fTransformation = new TGeoCombiTrans("");
  fConstituents = new TObjArray(20);
  
  // Set id
  SetUniqueID(id);
}


//______________________________________________________________________________
AliMUONGeometryEnvelope::AliMUONGeometryEnvelope()
 : TNamed(),
   fIsVirtual(0),
   fIsMANY(false),
   fCopyNo(0),
   fTransformation(0),
   fConstituents(0)
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONGeometryEnvelope::~AliMUONGeometryEnvelope() 
{
/// Destructor

  // Add deleting rotation matrices 
  
  delete fTransformation;
  
  if (fConstituents) {
    fConstituents->Delete();
    delete fConstituents;
  }  
}

//
// public methods
//

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituent(const TString& name, Int_t copyNo) 
{
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  fConstituents->Add(new AliMUONGeometryConstituent(name, copyNo, 0, 0));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituent(const TString& name, Int_t copyNo, 
                                          const TGeoTranslation& translation)
{
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  fConstituents
    ->Add(new AliMUONGeometryConstituent(name, copyNo, translation, 0, 0));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituent(const TString& name, Int_t copyNo,
                                          const TGeoTranslation& translation,
		                          const TGeoRotation& rotation)
{
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  fConstituents
    ->Add(new AliMUONGeometryConstituent(
                     name, copyNo, translation, rotation, 0, 0));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituent(const TString& name, Int_t copyNo,
                                          const TGeoCombiTrans& transform )
{
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  fConstituents
    ->Add(new AliMUONGeometryConstituent(
                     name, copyNo, transform, 0, 0));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituentParam(const TString& name, 
                                  Int_t copyNo, Int_t npar, Double_t* param) 
{
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  fConstituents
    ->Add(new AliMUONGeometryConstituent(name, copyNo, npar, param));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituentParam(const TString& name, 
                                  Int_t copyNo, const TGeoTranslation& translation,
				  Int_t npar, Double_t* param)
{
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  fConstituents
    ->Add(new AliMUONGeometryConstituent(
                     name, copyNo, translation, npar, param));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituentParam(const TString& name, 
                                  Int_t copyNo, const TGeoTranslation& translation,
		                  const TGeoRotation& rotation, 
				  Int_t npar, Double_t* param)
{
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  fConstituents
    ->Add(new AliMUONGeometryConstituent(
                     name, copyNo, translation, rotation, npar, param));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituentParam(const TString& name, 
                                  Int_t copyNo, 
				  const TGeoCombiTrans& transform, 
				  Int_t npar, Double_t* param)
{
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  fConstituents
    ->Add(new AliMUONGeometryConstituent(
                     name, copyNo, transform, npar, param));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::SetTranslation(const TGeoTranslation& translation)
{
/// Set the envelope position

  fTransformation
    ->SetTranslation(const_cast<Double_t*>(translation.GetTranslation()));
}  

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::SetRotation(const TGeoRotation& rotation)
{
/// Set the envelope rotation

  TGeoRotation* rot = new TGeoRotation();
  rot->SetMatrix(const_cast<Double_t*>(rotation.GetRotationMatrix()));

  fTransformation->SetRotation(rot);
}  

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::SetTransform(const TGeoCombiTrans& transform)
{
/// Set the envelope transformation

  fTransformation
    ->SetTranslation(const_cast<Double_t*>(transform.GetTranslation()));

  TGeoRotation* rot = new TGeoRotation();
  rot->SetMatrix(const_cast<Double_t*>(transform.GetRotationMatrix()));

  fTransformation->SetRotation(rot);
}  
