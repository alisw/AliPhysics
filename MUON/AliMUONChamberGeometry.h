/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

//
// Class AliMUONChamberGeometry
// -----------------------------
// Class for definititon of the MUON chamber positions in ALIC.
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_CHAMBER_GEOMETRY_H
#define ALI_MUON_CHAMBER_GEOMETRY_H

#include <TObject.h>
#include <TString.h>

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;
class TObjArray;
class TArrayI;

class AliMUONChamber;
class AliMUONGeometryEnvelope;
class AliMUONGeometryEnvelopeStore;
class AliMUONGeometryTransformStore;
class AliMUONGeometrySVMap;

class AliMUONChamberGeometry : public TObject
{
  public:
    AliMUONChamberGeometry(Int_t chamberId);
    AliMUONChamberGeometry();
    virtual ~AliMUONChamberGeometry();

    // methods
    void  SetMotherVolume(const TString& motherVolumeName);
    void  SetTranslation(const TGeoTranslation& translation);
    void  SetRotation(const TGeoRotation& rotation);
    
    void  SetSensitiveVolume(Int_t volId);
    void  SetSensitiveVolume(const TString& name);
    void  SetAlign(Bool_t align);
 
    // get methods
    TString                GetMotherVolume() const;
    const TGeoCombiTrans*  GetTransformation() const;    
    AliMUONGeometryEnvelopeStore*  GetEnvelopeStore() const;
    AliMUONGeometryTransformStore* GetTransformStore() const;
    AliMUONGeometrySVMap*          GetSVMap() const;
    Bool_t IsSensitiveVolume(Int_t volId) const; 
    Bool_t IsSensitiveVolume(const TString& volName) const; 
//
//Int_t  GetDEVolId(Int_t svVolId) const; 

  protected:
    AliMUONChamberGeometry(const AliMUONChamberGeometry& rhs);
    // operators  
    AliMUONChamberGeometry& operator = (const AliMUONChamberGeometry& rhs);

  private:
    // methods
    Int_t  GetSVIndex(Int_t svVolId) const; 
  
    // data members
    Int_t            fChamberId;     // the chamber Id
    TString          fMotherVolume;  // mother volume name
    Int_t            fNofSVs;        // number of sensitive volumes   
    TArrayI*         fSVVolumeIds;   // densitive volumes IDs  
    TGeoCombiTrans*  fTransformation;// the chamber transformation wrt to mother
                                     // volume
    AliMUONGeometryTransformStore* fDETransforms; // det elements transformations
    AliMUONGeometryEnvelopeStore*  fEnvelopes;    // envelopes                                 
    AliMUONGeometrySVMap*          fSVMap;        // sensitive volumes map
 
  ClassDef(AliMUONChamberGeometry,2) // MUON chamber geometry base class
};

// inline functions

inline void  
AliMUONChamberGeometry::SetMotherVolume(const TString& motherVolumeName)
{ fMotherVolume = motherVolumeName; }

inline TString  AliMUONChamberGeometry::GetMotherVolume() const
{ return fMotherVolume; }

inline const TGeoCombiTrans* AliMUONChamberGeometry::GetTransformation() const 
{ return fTransformation; }

inline  AliMUONGeometryEnvelopeStore*  
AliMUONChamberGeometry::GetEnvelopeStore() const
{ return fEnvelopes; }

inline AliMUONGeometryTransformStore*  
AliMUONChamberGeometry::GetTransformStore() const
{ return fDETransforms; }

inline AliMUONGeometrySVMap* AliMUONChamberGeometry::GetSVMap() const
{ return fSVMap; }

#endif //ALI_MUON_CHAMBER_GEOMETRY_H
