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

// ----------------------------------
// Class AliMUONGeometryEnvelopeStore
// ----------------------------------
// Class for definititon of the temporary volume envelopes
// used in geometry construction
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryEnvelope.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"

#include "AliMpExMap.h"

#include "AliLog.h"

#include <TGeoMatrix.h>
#include <TObjArray.h>
#include <Riostream.h>
#include <TString.h>

/// \cond CLASSIMP
ClassImp(AliMUONGeometryEnvelopeStore)
/// \endcond

//______________________________________________________________________________
AliMUONGeometryEnvelopeStore::AliMUONGeometryEnvelopeStore(
                                    AliMpExMap* detElements)
 : TObject(),
   fEnvelopes(0),
   fDetElements(detElements),
   fReferenceFrame(),
   fDebug(false),
   fAlign(false)
{
/// Standard constructor

  fEnvelopes = new TObjArray(100);
}


//______________________________________________________________________________
AliMUONGeometryEnvelopeStore::AliMUONGeometryEnvelopeStore()
 : TObject(),
   fEnvelopes(0),
   fDetElements(0),
   fReferenceFrame(),
   fDebug(false),
   fAlign(false)
{
/// Default constructor
}


//______________________________________________________________________________
AliMUONGeometryEnvelopeStore::~AliMUONGeometryEnvelopeStore() 
{
/// Destructor

  // Add deleting rotation matrices 
  
  if (fEnvelopes) {
    fEnvelopes->Delete();
    delete fEnvelopes;
  }  
}

//
// private methods
//

//______________________________________________________________________________
TGeoHMatrix 
AliMUONGeometryEnvelopeStore::ConvertDETransform(const TGeoHMatrix& transform) const
{
/// Convert transformation into the reference frame

  if ( fReferenceFrame.IsIdentity() )
    return transform;
  else  {
    return AliMUONGeometryBuilder::Multiply( fReferenceFrame.Inverse(),
  				  	     transform );  
  }			    
}

//______________________________________________________________________________
AliMUONGeometryEnvelope* 
AliMUONGeometryEnvelopeStore::FindEnvelope(const TString& name) const
{
/// Find the envelope specified by name.

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
/// Find transformation by the detection element	Id (if not 0)
/// (= unique ID of enevelope) and set it to the envelope.
/// Return true if transformation is applied, false otherwise.

  Int_t detElemId = envelope->GetUniqueID();
  if (detElemId == 0) return false;
  
  AliMUONGeometryDetElement* detElement 
    = (AliMUONGeometryDetElement*) fDetElements->GetValue(detElemId);
  if (!detElement) {
    AliWarning("Transformation not found.");
    return false;
  };

  // Apply frame transform
  TGeoHMatrix newTransform 
    = ConvertDETransform(*(detElement->GetLocalTransformation()));

  envelope->SetTransform(newTransform);
  
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
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

  if (!isVirtual) AliDebug(1,Form("Adding non-virtual envelope %s id %d",name.Data(),id));
//  else AliDebug(1,Form("Adding virtual envelope %s id %d",name.Data(),id));

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
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

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
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

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
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

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
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

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
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

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
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

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
/// Add the volume with the specified name and transformation
/// to the list of envelopes.

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
/// Add the volume with the specified name and transformation
/// as a constituent of the envelope envName.

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
/// Add the volume with the specified name and transformation
/// as a constituent of the envelope envName.

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
/// Add the volume with the specified name and transformation
/// as a constituent of the envelope envName.

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
/// Add the volume with the specified name and transformation
/// as a constituent of the envelope envName.

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
/// Add the volume with the specified name and transformation
/// as a constituent of the envelope envName.

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
/// Add the volume with the specified name and transformation
/// as a constituent of the envelope envName.

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
/// Add the volume with the specified name and transformation
/// as a constituent of the envelope envName.

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
/// Add the volume with the specified name and transformation
/// as a constituent of the envelope envName.

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

//______________________________________________________________________________
Int_t AliMUONGeometryEnvelopeStore::GetNofDetElements() const
{
/// Return the number od envelopes with detElemId>0.

  Int_t nofDetElems = 0;
  
  for(Int_t i=0; i<fEnvelopes->GetEntriesFast(); i++) 
    if ( fEnvelopes->At(i)->GetUniqueID() > 0 ) nofDetElems++;
  
  return nofDetElems;
}
