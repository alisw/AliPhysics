/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004

/// \ingroup geometry
/// \class AliMUONVGeometryBuilder
/// \brief Abstract base class for geometry construction per module(s)
///
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_V_GEOMETRY_BUILDER_H
#define ALI_MUON_V_GEOMETRY_BUILDER_H

#include <fstream>

#include <TObject.h>
#include <TObjArray.h>
#include <TGeoMatrix.h>

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;

class AliMUONGeometryModule;
class AliMUONGeometryEnvelopeStore;
class AliMUONStringIntMap;

class AliMUONVGeometryBuilder : public TObject
{
  public:
    AliMUONVGeometryBuilder(Int_t firstModuleId, Int_t nofModules);
    AliMUONVGeometryBuilder();
    virtual ~AliMUONVGeometryBuilder();
  
    // methods
    void  SetReferenceFrame(const TGeoCombiTrans& referenceFrame);
    void  RebuildSVMaps(Bool_t withEnvelopes = true) const;
    void  CreateDetElements() const;

                  /// Function to be overriden in a concrete chamber/station
		  /// geometry builder class.
		  /// Only materials that are not defined in the common
		  /// functions should be defined here.
    virtual void CreateMaterials() {}  // make = 0; ?

                  /// Function to be overriden in a concrete chamber/station
		  /// geometry builder class. \n
		  /// The geometry built there should not be placed
		  /// in ALIC; but all volumes going to ALIC
		  /// have to be added as envelopes to the chamber
		  /// geometries
		  /// (They will be then placed automatically 
		  /// usind the provided transformation.
    virtual void CreateGeometry() = 0;

                  /// Function to be overriden in a concrete chamber/station
		  /// geometry class. \n
		  /// The transformation of each chamber(s) wrt ALICE
		  /// should be defined and set to its geometry class. 
    virtual void SetTransformations() = 0;

                  /// Function to be overriden in a concrete chamber/station
		  /// geometry class. \n
		  /// The sensitive volumes Ids for each chamber
		  /// should be defined and set to its geometry class. 
    virtual void SetSensitiveVolumes() = 0;

                  /// Function to be overriden (and return false) 
		  /// in the concrete geometry builder classes 
		  /// which are already defined in the new ALICE
		  /// coordinate frame
    virtual bool ApplyGlobalTransformation() { return true; }

    // access to module geometries
    Int_t  NofGeometries() const;
    AliMUONGeometryModule* Geometry(Int_t i) const;
                  // In difference from protected GetGeometry()
		  // this function access geometry via index and not
		  // via moduleId

  protected:
    // methods
    AliMUONGeometryModule*         GetGeometry(Int_t moduleId) const;
    AliMUONGeometryEnvelopeStore*  GetEnvelopes(Int_t moduleId) const;
    AliMUONStringIntMap*           GetSVMap(Int_t moduleId) const;
    Int_t                          GetModuleId(const TString& envName) const;
    
    // set module transformation
    void SetTranslation(Int_t moduleId, 
                        const TGeoTranslation& translation);
    void SetTransformation(Int_t moduleId, 
                        const TGeoTranslation& translation,
			const TGeoRotation& rotation);
			
    // set volumes 
    void SetVolume(Int_t moduleId, const TString& volumeName, 
                   Bool_t isVirtual = false);			
    void SetMotherVolume(Int_t moduleId, const TString& volumeName);			
    
  private:
    //methods
    AliMUONVGeometryBuilder(const AliMUONVGeometryBuilder& rhs);
    AliMUONVGeometryBuilder& operator = (const AliMUONVGeometryBuilder& rhs);

    TGeoHMatrix ConvertTransform(const TGeoHMatrix& transform) const;
    TGeoHMatrix ConvertDETransform(const TGeoHMatrix& transform) const;
    TString     ComposePath(const TString& volName, Int_t copyNo) const; 
    void        MapSV(const TString& path0, 
                      const TString& volName, Int_t detElemId) const;

    // data members
    TObjArray*  fGeometryModules;   ///< \brief the modules geometries that will be built
                                    /// by this builder				    
    TGeoCombiTrans fReferenceFrame; ///< \brief the transformation from the builder 
                                    /// reference frame to that of the transform 
				    /// data files
				        
  ClassDef(AliMUONVGeometryBuilder,4) // MUON chamber geometry base class
};

// inline functions

/// Return the number of geometry modules
inline Int_t  AliMUONVGeometryBuilder::NofGeometries() const
{ return fGeometryModules->GetEntriesFast(); }

/// Return the \a i th geometry module
inline AliMUONGeometryModule* AliMUONVGeometryBuilder::Geometry(Int_t i) const
{ return (AliMUONGeometryModule*)fGeometryModules->At(i); }

#endif //ALI_MUON_V_GEOMETRY_BUILDER_H
