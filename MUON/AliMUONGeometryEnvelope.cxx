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

ClassImp(AliMUONGeometryEnvelope)

//______________________________________________________________________________
AliMUONGeometryEnvelope::AliMUONGeometryEnvelope(const TString& name, 
                                                 Bool_t isVirtual,
						 const char* only)
 : TNamed(name, name),
   fIsVirtual(isVirtual),
   fIsMANY(false),
   fCopyNo(0),
   fTransformation(0),
   fConstituents(0)
{
// Standard constructor

  if (TString(only) == TString("MANY")) fIsMANY = true;

  // Create the envelope transformation
  fTransformation = new TGeoCombiTrans("");
  fConstituents = new TObjArray(20);
}


//______________________________________________________________________________
AliMUONGeometryEnvelope::AliMUONGeometryEnvelope(const TString& name, 
                                                 Int_t copyNo, 
						 const char* only)
 : TNamed(name, name),
   fIsVirtual(false),
   fIsMANY(false),
   fCopyNo(copyNo),
   fTransformation(0),
   fConstituents(0)
{
// Standard constructor

  if (TString(only) == TString("MANY")) fIsMANY = true;

  // Create the envelope transformation
  fTransformation = new TGeoCombiTrans("");
  fConstituents = new TObjArray(20);
}


//______________________________________________________________________________
AliMUONGeometryEnvelope::AliMUONGeometryEnvelope()
 : TNamed(),
   fIsVirtual(0),
   fTransformation(0),
   fConstituents(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONGeometryEnvelope::AliMUONGeometryEnvelope(
                                        const AliMUONGeometryEnvelope& rhs)
  : TNamed(rhs)
{
  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONGeometryEnvelope::~AliMUONGeometryEnvelope() 
{
//
  // Add deleting rotation matrices 
  
  delete fTransformation;
  
  if (fConstituents) {
    fConstituents->Delete();
    delete fConstituents;
  }  
}

//______________________________________________________________________________
AliMUONGeometryEnvelope& 
AliMUONGeometryEnvelope::operator = (const AliMUONGeometryEnvelope& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  Fatal("operator=", 
        "Assignment operator is not implemented.");
    
  return *this;  
}

//
// public methods
//

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituent(const TString& name, Int_t copyNo) 
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  fConstituents->Add(new AliMUONGeometryConstituent(name, copyNo, 0, 0));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituent(const TString& name, Int_t copyNo, 
                                          const TGeoTranslation& translation)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  fConstituents
    ->Add(new AliMUONGeometryConstituent(name, copyNo, translation, 0, 0));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituent(const TString& name, Int_t copyNo,
                                          const TGeoTranslation& translation,
		                          const TGeoRotation& rotation)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  fConstituents
    ->Add(new AliMUONGeometryConstituent(
                     name, copyNo, translation, rotation, 0, 0));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituentParam(const TString& name, 
                                  Int_t copyNo, Int_t npar, Double_t* param) 
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  fConstituents
    ->Add(new AliMUONGeometryConstituent(name, copyNo, npar, param));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::AddConstituentParam(const TString& name, 
                                  Int_t copyNo, const TGeoTranslation& translation,
				  Int_t npar, Double_t* param)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

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
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  fConstituents
    ->Add(new AliMUONGeometryConstituent(
                     name, copyNo, translation, rotation, npar, param));
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::SetTranslation(const TGeoTranslation& translation)
{
// Sets the chamber position wrt ALIC.
// ---

  fTransformation
    ->SetTranslation(const_cast<Double_t*>(translation.GetTranslation()));
}  

//______________________________________________________________________________
void  AliMUONGeometryEnvelope::SetRotation(const TGeoRotation& rotation)
{
// Sets the chamber rotation wrt ALIC.
// ---

  TGeoRotation* rot = new TGeoRotation();
  rot->SetMatrix(const_cast<Double_t*>(rotation.GetRotationMatrix()));

  fTransformation->SetRotation(rot);
}  
