/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

//
// Class AliMUONGeometryModule
// -----------------------------
// Class for definition of the detector module parameters
// (the transformations of detection elements, mapping between
//  sensitive volumes and detection elements).
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_MODULE_GEOMETRY_H
#define ALI_MUON_MODULE_GEOMETRY_H

#include <TObject.h>
#include <TString.h>

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;
class TObjArray;
class TArrayI;

class AliMUONGeometryEnvelope;
class AliMUONGeometryEnvelopeStore;
class AliMUONGeometryDetElement;
class AliMUONGeometryStore;
class AliMUONGeometrySVMap;
class AliMUONVGeometryDEIndexing;

class AliMUONGeometryModule : public TObject
{
  public:
    AliMUONGeometryModule(Int_t moduleId);
    AliMUONGeometryModule();
    virtual ~AliMUONGeometryModule();

    // methods
    void Global2Local(Int_t detElemId,
                 Float_t xg, Float_t yg, Float_t zg, 
                 Float_t& xl, Float_t& yl, Float_t& zl) const;
    void Global2Local(Int_t detElemId,
                 Double_t xg, Double_t yg, Double_t zg, 
                 Double_t& xl, Double_t& yl, Double_t& zl) const;

    void Local2Global(Int_t detElemId,
                 Float_t xl, Float_t yl, Float_t zl, 
                 Float_t& xg, Float_t& yg, Float_t& zg) const;
    void Local2Global(Int_t detElemId,
                 Double_t xl, Double_t yl, Double_t zl, 
                 Double_t& xg, Double_t& yg, Double_t& zg) const;

    // set methods
    void  SetMotherVolume(const TString& motherVolumeName);
    void  SetVolume(const TString& volumeName);
    void  SetTranslation(const TGeoTranslation& translation);
    void  SetRotation(const TGeoRotation& rotation);
    
    void  SetSensitiveVolume(Int_t volId);
    void  SetSensitiveVolume(const TString& name);
    void  SetAlign(Bool_t align);
 
    // get methods
    Bool_t                 IsVirtual() const;  
    Int_t                  GetModuleId() const;
    TString                GetMotherVolume() const;
    TString                GetVolume() const;
    const TGeoCombiTrans*  GetTransformation() const;    
    AliMUONGeometryDetElement* FindBySensitiveVolume(
                                         const TString& volumePath) const;
    AliMUONVGeometryDEIndexing*    GetDEIndexing() const;
    AliMUONGeometryEnvelopeStore*  GetEnvelopeStore() const;
    AliMUONGeometryStore*          GetDetElementStore() const;
    AliMUONGeometryDetElement*     GetDetElement(Int_t detElemId) const;    
    AliMUONGeometrySVMap*          GetSVMap() const;
    Bool_t IsSensitiveVolume(Int_t volId) const; 
    Bool_t IsSensitiveVolume(const TString& volName) const; 

  protected:
    AliMUONGeometryModule(const AliMUONGeometryModule& rhs);
    // operators  
    AliMUONGeometryModule& operator = (const AliMUONGeometryModule& rhs);

  private:
    // methods
    Int_t  GetSVIndex(Int_t svVolId) const; 
  
    // data members
    Bool_t           fIsVirtual;     // true if module is not represented
                                     // by a real volume
    Int_t            fModuleId;      // the module Id
    TString          fMotherVolume;  // mother volume name
    TString          fVolume;        // the volume name if not virtual
    Int_t            fNofSVs;        // number of sensitive volumes   
    TArrayI*         fSVVolumeIds;   // densitive volumes IDs  
    TGeoCombiTrans*  fTransformation;// the module transformation wrt to mother
                                     // volume
    AliMUONGeometryEnvelopeStore* fEnvelopes;  // envelopes                                 
    AliMUONVGeometryDEIndexing*   fDEIndexing; // DE indexing
    AliMUONGeometryStore*         fDetElements;// detection elements
    AliMUONGeometrySVMap*         fSVMap;      // sensitive volumes map
 
  ClassDef(AliMUONGeometryModule,2) // MUON geometry module class
};

// inline functions

inline void  
AliMUONGeometryModule::SetMotherVolume(const TString& motherVolumeName)
{ fMotherVolume = motherVolumeName; }

inline Bool_t AliMUONGeometryModule::IsVirtual() const
{ return fIsVirtual; }  

inline Int_t  AliMUONGeometryModule::GetModuleId() const
{ return fModuleId; }

inline TString  AliMUONGeometryModule::GetMotherVolume() const
{ return fMotherVolume; }

inline TString  AliMUONGeometryModule::GetVolume() const
{ return fVolume; }

inline const TGeoCombiTrans* AliMUONGeometryModule::GetTransformation() const 
{ return fTransformation; }

inline  AliMUONGeometryEnvelopeStore* 
AliMUONGeometryModule::GetEnvelopeStore() const
{ return fEnvelopes; }

inline  AliMUONVGeometryDEIndexing* 
AliMUONGeometryModule::GetDEIndexing() const 
{ return fDEIndexing; }

inline  AliMUONGeometryStore* AliMUONGeometryModule::GetDetElementStore() const
{ return fDetElements; }

inline AliMUONGeometrySVMap* AliMUONGeometryModule::GetSVMap() const
{ return fSVMap; }

#endif //ALI_MUON_MODULE_GEOMETRY_H
