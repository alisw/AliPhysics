// $Id$
//
// Class AliMUONGeometryEnvelopeStore
// ----------------------------------
// Class for definititon of the temporary volume envelopes
// used in geometry construction
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <TObjArray.h>
#include <TArrayI.h>
#include <Riostream.h>

#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryTransformStore.h"
#include "AliMUONGeometryEnvelope.h"
#include "AliMUONConstants.h"

ClassImp(AliMUONGeometryEnvelopeStore)

//______________________________________________________________________________
AliMUONGeometryEnvelopeStore::AliMUONGeometryEnvelopeStore(
                                 AliMUONGeometryTransformStore* transforms)
 : TObject(),
   fDETransforms(transforms),
   fEnvelopes(0),
   fDebug(false),
   fAlign(false)
{
// Standard constructor

  fEnvelopes = new TObjArray(100);
}


//______________________________________________________________________________
AliMUONGeometryEnvelopeStore::AliMUONGeometryEnvelopeStore()
 : TObject(),
   fDETransforms(0),
   fEnvelopes(0),
   fDebug(false),
   fAlign(false)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONGeometryEnvelopeStore::AliMUONGeometryEnvelopeStore(const AliMUONGeometryEnvelopeStore& rhs)
  : TObject(rhs)
{
  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONGeometryEnvelopeStore::~AliMUONGeometryEnvelopeStore() 
{
//

  // Add deleting rotation matrices 
  
  if (fEnvelopes) {
    fEnvelopes->Delete();
    delete fEnvelopes;
  }  
}

//______________________________________________________________________________
AliMUONGeometryEnvelopeStore& 
AliMUONGeometryEnvelopeStore::operator = (const AliMUONGeometryEnvelopeStore& rhs) 
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
AliMUONGeometryEnvelopeStore::FindEnvelope(const TString& name) const
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

//______________________________________________________________________________
Bool_t AliMUONGeometryEnvelopeStore::AlignEnvelope(
                                          AliMUONGeometryEnvelope* envelope) const
{
// Find transformation by the detection element	Id (if not 0)
// (= unique ID of enevelope) and set it to the envelope.
// Return true if transformation is applied, false otherwise.
// ---

  Int_t detElemId = envelope->GetUniqueID();
  if (detElemId == 0) return false;
  
  const TGeoCombiTrans* kTransform = fDETransforms->Get(detElemId);
  if (!kTransform) {
    Warning("AlignEnvelope", "Transformation not found.");
    return false;
  };

  envelope->SetTransform(*kTransform);
  return true;
}  

//
// public methods
//

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelope(const TString& name, 
                                          Int_t id, 
                                          Bool_t isVirtual,
                                          const char* only) 
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding ";
    if (!isVirtual) cout << " non-";
    cout << "virtual envelope " << name 
         << "  id " << id << endl;
  }  

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, id, isVirtual, only);
    
  if (fAlign) AlignEnvelope(envelope); 

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelope(const TString& name, 
                                          Int_t id, 
                                          Bool_t isVirtual,
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
         << "  id " << id
         << " with translation" << endl;
  }  

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, id, isVirtual, only);  
    
  Bool_t aligned = false;
  if (fAlign) aligned = AlignEnvelope(envelope); 

  if  (!aligned)
    envelope->SetTranslation(translation);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelope(const TString& name, 
                                          Int_t id, 
                                          Bool_t isVirtual, 
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
         << "  id " << id
         << " with translation and rotation" << endl;
  }  

/*
  cout << "Adding env...  name: " << name;
   
   const Double_t* xyz = translation.GetTranslation();
   cout << "  translation: " << xyz[0] << ", " << xyz[1] << ", " << xyz[2]
        << "  rotation: ";
	   
   Double_t a1, a2, a3, a4, a5, a6;
   rotation.GetAngles(a1, a2, a3, a4, a5, a6);
   cout << a1 << ", " << a2 << ", " << a3 << ", " << a4 << ", " << a5 << ", " << a6 << endl; 	   	 
*/
  // fEnvelopes->Add(new TGeoCombiTrans(name, translation, rotation));
           // would be nice to be so simple 

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, id, isVirtual, only);

  Bool_t aligned = false;
  if (fAlign) aligned = AlignEnvelope(envelope); 

  if  (!aligned) {
    envelope->SetRotation(rotation);
    envelope->SetTranslation(translation);
  }  

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelope(const TString& name, 
                                          Int_t id, 
                                          Bool_t isVirtual, 
                                          const TGeoCombiTrans& transform,
					  const char* only)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding ";
    if (!isVirtual) cout << " non-";
    cout << "virtual envelope " << name 
         << "  id " << id
         << " with transformation" << endl;
  }  

  // fEnvelopes->Add(new TGeoCombiTrans(name, translation, rotation));
           // would be nice to be so simple 

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, id, isVirtual, only);

  Bool_t aligned = false;
  if (fAlign) aligned = AlignEnvelope(envelope); 

  if  (!aligned) 
    envelope->SetTransform(transform);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelope(const TString& name, 
                                          Int_t id, 
                                          Int_t copyNo,
                                          const char* only) 
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding "
         << " non-virtual envelope " << name 
         << "  id " << id
         << " with copyNo " << copyNo << endl;
   }  

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, id, copyNo, only);

  if (fAlign) AlignEnvelope(envelope); 

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelope(const TString& name, 
                                          Int_t id, 
                                          Int_t copyNo,
                                          const TGeoTranslation& translation,
					  const char* only)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding "
         << " non-virtual envelope " << name 
         << "  id " << id
         << " with copyNo " << copyNo
         << " with translation " << endl;
  }  

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, id, copyNo, only);

  Bool_t aligned = false;
  if (fAlign) aligned = AlignEnvelope(envelope); 

  if  (!aligned) 
    envelope->SetTranslation(translation);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelope(const TString& name, 
                                          Int_t id, 
                                          Int_t copyNo, 
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
         << "  id " << id
         << " with copyNo " << copyNo
         << " with translation and rotation" << endl;
  }  

  // fEnvelopes->Add(new TGeoCombiTrans(name, translation, rotation));
           // would be nice to be so simple 

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, id, copyNo, only);

  Bool_t aligned = false;
  if (fAlign) aligned = AlignEnvelope(envelope); 

  if  (!aligned) {
    envelope->SetRotation(rotation);
    envelope->SetTranslation(translation);
  }  

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelope(const TString& name, 
                                          Int_t id, 
                                          Int_t copyNo, 
                                          const TGeoCombiTrans& transform,
					  const char* only)
{
// Adds the volume with the specified name and transformation
// to the list of envelopes.
// ---  					   

  if (fDebug) {
    cout << "... Adding "
         << " non-virtual envelope " << name 
         << "  id " << id
         << " with copyNo " << copyNo
         << " with translation and rotation" << endl;
  }  

  // fEnvelopes->Add(new TGeoCombiTrans(name, translation, rotation));
           // would be nice to be so simple 

  AliMUONGeometryEnvelope* envelope 
    = new AliMUONGeometryEnvelope(name, id, copyNo, only);

  Bool_t aligned = false;
  if (fAlign) aligned = AlignEnvelope(envelope); 

  if  (!aligned)
    envelope->SetTransform(transform);

  fEnvelopes->Add(envelope);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelopeConstituent(const TString& name, 
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
void  AliMUONGeometryEnvelopeStore::AddEnvelopeConstituent(const TString& name, 
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
void  AliMUONGeometryEnvelopeStore::AddEnvelopeConstituent(const TString& name, 
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
void  AliMUONGeometryEnvelopeStore::AddEnvelopeConstituent(const TString& name, 
                                          const TString& envName, Int_t copyNo, 
                                          const TGeoCombiTrans& transform)
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
   
  envelope->AddConstituent(name, copyNo, transform);
}

//______________________________________________________________________________
void  AliMUONGeometryEnvelopeStore::AddEnvelopeConstituentParam(const TString& name, 
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
void  AliMUONGeometryEnvelopeStore::AddEnvelopeConstituentParam(const TString& name, 
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
void  AliMUONGeometryEnvelopeStore::AddEnvelopeConstituentParam(const TString& name, 
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
void  AliMUONGeometryEnvelopeStore::AddEnvelopeConstituentParam(const TString& name, 
                                          const TString& envName, Int_t copyNo, 
                                          const TGeoCombiTrans& transform,
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
   
  envelope->AddConstituentParam(name, copyNo, transform, npar, param);
}

