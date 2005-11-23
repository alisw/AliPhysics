/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup geometry
/// \class AliMUONGeometryModule
/// \brief Geometry parameters for detector module
///
/// Class for definition of the detector module parameters
/// (the transformations of detection elements, mapping between
///  sensitive volumes and detection elements).
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_MODULE_H
#define ALI_MUON_GEOMETRY_MODULE_H

#include <TObject.h>
#include <TString.h>

#include "AliMUONGeometryModuleTransformer.h"

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

class AliMUONGeometryModule : public TObject
{
  public:
    AliMUONGeometryModule(Int_t moduleId);
    AliMUONGeometryModule();
    virtual ~AliMUONGeometryModule();

    // set methods
    //
    void  SetMotherVolume(const TString& motherVolumeName);
    void  SetVolume(const TString& volumeName);
    void  SetTransformation(const TGeoCombiTrans& transform);
    
    void  SetSensitiveVolume(Int_t volId);
    void  SetSensitiveVolume(const TString& name);
    void  SetAlign(Bool_t align);
 
    // get methods
    //
    Bool_t   IsVirtual() const;  
    Int_t    GetModuleId() const;
    TString  GetMotherVolume() const;
    TString  GetVolume() const;
    
    AliMUONGeometryDetElement* FindBySensitiveVolume(
                                         const TString& volumePath) const;
    Bool_t IsSensitiveVolume(Int_t volId) const; 
    Bool_t IsSensitiveVolume(const TString& volName) const; 

    AliMUONGeometryEnvelopeStore*     GetEnvelopeStore() const;
    AliMUONGeometrySVMap*             GetSVMap() const;
    AliMUONGeometryModuleTransformer* GetTransformer() const;

  protected:
    AliMUONGeometryModule(const AliMUONGeometryModule& rhs);
     AliMUONGeometryModule& operator = (const AliMUONGeometryModule& rhs);

  private:
    // methods
    Int_t  GetSVIndex(Int_t svVolId) const; 
  
    // data members
    Bool_t           fIsVirtual;     // true if module is not represented
                                     // by a real volume
    TString          fMotherVolume;  // mother volume name
    TString          fVolume;        // the volume name if not virtual
    Int_t            fNofSVs;        // number of sensitive volumes   
    TArrayI*         fSVVolumeIds;   // densitive volumes IDs  

    AliMUONGeometryEnvelopeStore*     fEnvelopes;  // envelopes                                 
    AliMUONGeometrySVMap*             fSVMap;      // sensitive volumes map
    AliMUONGeometryModuleTransformer* fTransformer;// geometry transformations
 
  ClassDef(AliMUONGeometryModule,3) // MUON geometry module class
};

// inline functions

inline 
void  AliMUONGeometryModule::SetMotherVolume(const TString& motherVolumeName)
{ fMotherVolume = motherVolumeName; }

inline Bool_t AliMUONGeometryModule::IsVirtual() const
{ return fIsVirtual; }  

inline Int_t  AliMUONGeometryModule::GetModuleId() const
{ return fTransformer->GetModuleId(); }

inline TString  AliMUONGeometryModule::GetMotherVolume() const
{ return fMotherVolume; }

inline TString  AliMUONGeometryModule::GetVolume() const
{ return fVolume; }

inline  
AliMUONGeometryEnvelopeStore* AliMUONGeometryModule::GetEnvelopeStore() const
{ return fEnvelopes; }

inline 
AliMUONGeometrySVMap* AliMUONGeometryModule::GetSVMap() const
{ return fSVMap; }

inline 
AliMUONGeometryModuleTransformer* AliMUONGeometryModule::GetTransformer() const
{ return fTransformer; }

#endif //ALI_MUON_GEOMETRY_MODULE_H
