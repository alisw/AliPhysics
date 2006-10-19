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
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_MODULE_H
#define ALI_MUON_GEOMETRY_MODULE_H

#include "AliMUONGeometryModuleTransformer.h"

#include <TObject.h>
#include <TString.h>

class AliMUONGeometryEnvelope;
class AliMUONGeometryEnvelopeStore;
class AliMUONGeometryDetElement;
class AliMUONStringIntMap;

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;
class TObjArray;
class TArrayI;

class AliMUONGeometryModule : public TObject
{
  public:
    AliMUONGeometryModule(Int_t moduleId);
    AliMUONGeometryModule();
    virtual ~AliMUONGeometryModule();

    // set methods
    //
    void  SetTransformation(const TGeoCombiTrans& transform);
    void  SetVolumePath(const TString& volumePath);
    void  SetIsVirtual(Bool_t isVirtual);
    
    void  SetSensitiveVolume(Int_t volId);
    void  SetSensitiveVolume(const TString& name);
    void  SetAlign(Bool_t align);
 
    // get methods
    //
    Bool_t   IsVirtual() const;  
    Int_t    GetModuleId() const;
    TString  GetVolumePath() const;
    
    AliMUONGeometryDetElement* FindBySensitiveVolume(
                                         const TString& volumePath) const;
    Bool_t IsSensitiveVolume(Int_t volId) const; 
    Bool_t IsSensitiveVolume(const TString& volName) const; 

    AliMUONGeometryEnvelopeStore*     GetEnvelopeStore() const;
    AliMUONStringIntMap*              GetSVMap() const;
    AliMUONGeometryModuleTransformer* GetTransformer() const;

  protected:
    AliMUONGeometryModule(const AliMUONGeometryModule& rhs);
    AliMUONGeometryModule& operator = (const AliMUONGeometryModule& rhs);

  private:
    // methods
    Int_t  GetSVIndex(Int_t svVolId) const; 
  
    // data members
    Bool_t           fIsVirtual;     ///< \brief true if module is not represented
                                     /// by a real volume
    Int_t            fNofSVs;        ///< number of sensitive volumes   
    TArrayI*         fSVVolumeIds;   ///< sensitive volumes IDs  

    AliMUONGeometryEnvelopeStore*     fEnvelopes;  ///< envelopes                                 
    AliMUONStringIntMap*              fSVMap;      ///< sensitive volumes map
    AliMUONGeometryModuleTransformer* fTransformer;///< geometry transformations
 
  ClassDef(AliMUONGeometryModule,4) // MUON geometry module class
};

// inline functions

/// Set virtuality (true if module is not represented by a real volume)
inline void  AliMUONGeometryModule::SetIsVirtual(Bool_t isVirtual)
{ fIsVirtual = isVirtual; }

/// Return true if module is not represented by a real volume
inline Bool_t AliMUONGeometryModule::IsVirtual() const
{ return fIsVirtual; }  

/// Return module ID
inline Int_t  AliMUONGeometryModule::GetModuleId() const
{ return fTransformer->GetModuleId(); }

/// Return the full path of aligned module volume or envelope in geometry
inline TString AliMUONGeometryModule::GetVolumePath() const
{ return fTransformer->GetVolumePath(); }

/// Return envelopes associated with this module
inline  
AliMUONGeometryEnvelopeStore* AliMUONGeometryModule::GetEnvelopeStore() const
{ return fEnvelopes; }

/// Return sensitive volume map
inline 
AliMUONStringIntMap* AliMUONGeometryModule::GetSVMap() const
{ return fSVMap; }

/// Return transformer
inline 
AliMUONGeometryModuleTransformer* AliMUONGeometryModule::GetTransformer() const
{ return fTransformer; }

#endif //ALI_MUON_GEOMETRY_MODULE_H
