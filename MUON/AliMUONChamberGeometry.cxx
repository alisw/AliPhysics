// $Id$
//
// Class AliMUONChamberGeometry
// -----------------------------
// Class for definititon of the MUON chamber positions in ALIC
// Author: Ivana Hrivnacova, IPN Orsay
// 23/01/2004

#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <TObjArray.h>
#include <TArrayI.h>
#include <Riostream.h>

#include "AliMUONChamberGeometry.h"
#include "AliMUONGeometryEnvelope.h"

ClassImp(AliMUONChamberGeometry)

//______________________________________________________________________________
AliMUONChamberGeometry::AliMUONChamberGeometry(Int_t chamberId)
 : TObject(),
   fChamberId(chamberId),
   fMotherVolume("ALIC"),
   fTransformation(0),
   fEnvelopes(0),
   fNofSensVolumeIds(0),
   fSensVolumeIds(0),
   fDebug(kFALSE)
{
// Standard constructor

  // Create the chamber transformation
  fTransformation = new TGeoCombiTrans("");
  fEnvelopes = new TObjArray(20);
  fSensVolumeIds = new TArrayI(20);
}


//______________________________________________________________________________
AliMUONChamberGeometry::AliMUONChamberGeometry()
 : TObject(),
   fChamberId(0),
   fMotherVolume(),
   fTransformation(0),
   fEnvelopes(0),
   fNofSensVolumeIds(0),
   fSensVolumeIds(0),
   fDebug(kFALSE)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONChamberGeometry::AliMUONChamberGeometry(const AliMUONChamberGeometry& rhs)
  : TObject(rhs)
{
  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONChamberGeometry::~AliMUONChamberGeometry() {
//

  // Add deleting rotation matrices 
  
  delete fTransformation;
  
  if (fEnvelopes) {
    fEnvelopes->Delete();
    delete fEnvelopes;
  }  
}

//______________________________________________________________________________
AliMUONChamberGeometry& 
AliMUONChamberGeometry::operator = (const AliMUONChamberGeometry& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  Fatal("operator=", 
        "Assignment operator is not implemented.");
    
  return *this;  
}

//
// private methods
//

//______________________________________________________________________________
AliMUONGeometryEnvelope* 
AliMUONChamberGeometry::FindEnvelope(const TString& name) const
{
// Finds the envelope specified by name.
// ---

  for (Int_t i=0; i<fEnvelopes->GetEntriesFast(); i++) {
    AliMUONGeometryEnvelope* envelope 
      = (AliMUONGeometryEnvelope*)fEnvelopes->At(i);
    
    if (envelope->GetName() == name) return envelope;
  }
  
  return 0;    
}  

//
// public methods
//

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelope(const TString& name, Bool_t isVirtual,
                                          const char* only) 
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding ";
    if (!isVirtual) cout << " non-";
    cout << "virtual envelope " << name << endl;
  }  

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, isVirtual, only);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelope(const TString& name, Bool_t isVirtual,
                                          const TGeoTranslation& translation,
					  const char* only)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding ";
    if (!isVirtual) cout << " non-";
    cout << "virtual envelope " << name 
         << " with translation" << endl;
  }  

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, isVirtual, only);
  envelope->SetTranslation(translation);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelope(const TString& name, Bool_t isVirtual, 
                                          const TGeoTranslation& translation,
		                          const TGeoRotation& rotation,
					  const char* only)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding ";
    if (!isVirtual) cout << " non-";
    cout << "virtual envelope " << name 
         << " with translation and rotation" << endl;
  }  

  // fEnvelopes->Add(new TGeoCombiTrans(name, translation, rotation));
           // would be nice to be so simple 

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, isVirtual, only);
  envelope->SetRotation(rotation);
  envelope->SetTranslation(translation);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelope(const TString& name, Int_t copyNo,
                                          const char* only) 
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding "
         << " non-virtual envelope " << name 
         << " with copyNo " << copyNo << endl;
   }  

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, copyNo, only);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelope(const TString& name, Int_t copyNo,
                                          const TGeoTranslation& translation,
					  const char* only)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding "
         << " non-virtual envelope " << name 
         << " with copyNo " << copyNo
         << " with translation " << endl;
  }  

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, copyNo, only);
  envelope->SetTranslation(translation);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelope(const TString& name, Int_t copyNo, 
                                          const TGeoTranslation& translation,
		                          const TGeoRotation& rotation,
					  const char* only)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding "
         << " non-virtual envelope " << name 
         << " with copyNo " << copyNo
         << " with translation and rotation" << endl;
  }  

  // fEnvelopes->Add(new TGeoCombiTrans(name, translation, rotation));
           // would be nice to be so simple 

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, copyNo, only);
  envelope->SetRotation(rotation);
  envelope->SetTranslation(translation);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelopeConstituent(const TString& name, 
                                         const TString& envName, Int_t copyNo) 
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding constituent " << name
         << " to envelope " << envName 
         << " with copyNo " << copyNo << endl;
  }  

  AliMUONGeometryEnvelope* envelope = FindEnvelope(envName);
  
  if (!envelope) {
    // add warning
    return;
  }  
   
  envelope->AddConstituent(name, copyNo);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelopeConstituent(const TString& name, 
                                          const TString& envName, Int_t copyNo,
                                          const TGeoTranslation& translation)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding constituent " << name
         << " to envelope " << envName 
         << " with copyNo " << copyNo
         << " with translation" << endl;
  }  

  AliMUONGeometryEnvelope* envelope = FindEnvelope(envName);
  
  if (!envelope) {
    // add warning
    return;
  }  
   
  envelope->AddConstituent(name, copyNo, translation);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelopeConstituent(const TString& name, 
                                          const TString& envName, Int_t copyNo, 
                                          const TGeoTranslation& translation,
		                          const TGeoRotation& rotation)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding constituent " << name
         << " to envelope " << envName 
         << " with copyNo " << copyNo
         << " with translation and rotation" << endl;
  }  

  AliMUONGeometryEnvelope* envelope = FindEnvelope(envName);
  
  if (!envelope) {
    // add warning
    return;
  }  
   
  envelope->AddConstituent(name, copyNo, translation, rotation);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelopeConstituentParam(const TString& name, 
                                         const TString& envName, Int_t copyNo,
					 Int_t npar, Double_t* param) 
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding parameterised constituent " << name
         << " to envelope " << envName 
         << " with copyNo " << copyNo << endl;
  }  

  AliMUONGeometryEnvelope* envelope = FindEnvelope(envName);
  
  if (!envelope) {
    // add warning
    return;
  }  
   
  envelope->AddConstituentParam(name, copyNo, npar, param);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelopeConstituentParam(const TString& name, 
                                          const TString& envName, Int_t copyNo,
                                          const TGeoTranslation& translation,
					  Int_t npar, Double_t* param)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding parameterised constituent " << name
         << " to envelope " << envName 
         << " with copyNo " << copyNo
         << " with translation" << endl;
  }  

  AliMUONGeometryEnvelope* envelope = FindEnvelope(envName);
  
  if (!envelope) {
    // add warning
    return;
  }  
   
  envelope->AddConstituentParam(name, copyNo, translation, npar, param);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::AddEnvelopeConstituentParam(const TString& name, 
                                          const TString& envName, Int_t copyNo, 
                                          const TGeoTranslation& translation,
		                          const TGeoRotation& rotation,
					  Int_t npar, Double_t* param)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding parameterised constituent " << name
         << " to envelope " << envName 
         << " with copyNo " << copyNo
         << " with translation and rotation" << endl;
  }  

  AliMUONGeometryEnvelope* envelope = FindEnvelope(envName);
  
  if (!envelope) {
    // add warning
    return;
  }  
   
  envelope->AddConstituentParam(name, copyNo, translation, rotation, npar, param);
}

//______________________________________________________________________________
void  AliMUONChamberGeometry::SetTranslation(const TGeoTranslation& translation)
{
// Sets the chamber position wrt ALIC.
// ---

  fTransformation
    ->SetTranslation(const_cast<Double_t*>(translation.GetTranslation()));
}  

//______________________________________________________________________________
void  AliMUONChamberGeometry::SetRotation(const TGeoRotation& rotation)
{
// Sets the chamber rotation wrt ALIC.
// ---

  TGeoRotation* rot = new TGeoRotation();
  rot->SetMatrix(const_cast<Double_t*>(rotation.GetRotationMatrix()));

  fTransformation->SetRotation(rot);
}  

//______________________________________________________________________________
void  AliMUONChamberGeometry::SetSensitiveVolume(Int_t volId)
{
// Adds the volume specified by volId to the list of sensitive
// volumes
  
  fSensVolumeIds->AddAt(volId,fNofSensVolumeIds++);
}      

//______________________________________________________________________________
void  AliMUONChamberGeometry::SetSensitiveVolume(const TString& volName)
{
// Adds the volume specified by volId to the list of sensitive
// volumes

  fSensVolumeIds->AddAt(gMC->VolId(volName),fNofSensVolumeIds++);
}      

//______________________________________________________________________________
Bool_t AliMUONChamberGeometry::IsSensitiveVolume(Int_t volId) const
{
// Checks if the volume specified by volId is present in the list
// of sensitive volumes.

  for (Int_t i=0; i<fNofSensVolumeIds; i++) {
      if (fSensVolumeIds->At(i) == volId) return kTRUE;
  }
  return kFALSE;
}
