/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004

/// \ingroup geometry
/// \class AliMUONVGeometryBuilder
/// \brief Abstract base class for geometry construction per module(s)
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_V_GEOMETRY_BUILDER_H
#define ALI_MUON_V_GEOMETRY_BUILDER_H

#include <fstream>

#include <TObject.h>
#include <TObjArray.h>

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;

class AliMUONGeometryModule;
class AliMUONGeometryEnvelopeStore;
class AliMUONGeometryStore;
class AliMUONGeometrySVMap;

class AliMUONVGeometryBuilder : public TObject
{
  public:
    AliMUONVGeometryBuilder(const TString& fileName,
                            AliMUONGeometryModule* mg1,
                            AliMUONGeometryModule* mg2 = 0,
                            AliMUONGeometryModule* mg3 = 0,
                            AliMUONGeometryModule* mg4 = 0,
                            AliMUONGeometryModule* mg5 = 0,
                            AliMUONGeometryModule* mg6 = 0);
    AliMUONVGeometryBuilder();
    virtual ~AliMUONVGeometryBuilder();
  
    // methods
    virtual void  FillTransformations() const;
    virtual void  RebuildSVMaps() const;
    virtual Bool_t  ReadTransformations() const;
    virtual Bool_t  ReadSVMap() const;
    virtual Bool_t  WriteTransformations() const;
    virtual Bool_t  WriteSVMap(Bool_t rebuild) const;
    
    virtual void CreateMaterials() {}  // make = 0; ?
                  // Function to be overriden in a concrete chamber/station
		  // geometry builder class.
		  // Only materials that are not defined in the common
		  // functions should be defined here.
    virtual void CreateGeometry() = 0;
                  // Function to be overriden in a concrete chamber/station
		  // geometry builder class.
		  // The geometry built there should not be placed
		  // in ALIC; but all volumes going to ALIC
		  // have to be added as envelopes to the chamber
		  // geometries
		  // (They will be then placed automatically 
		  // usind the provided transformation.
    virtual void SetTransformations() = 0;
                  // Function to be overriden in a concrete chamber/station
		  // geometry class.
		  // The transformation of each chamber(s) wrt ALICE
		  // should be defined and set to its geometry class. 
    virtual void SetSensitiveVolumes() = 0;
                  // Function to be overriden in a concrete chamber/station
		  // geometry class.
		  // The sensitive volumes Ids for each chamber
		  // should be defined and set to its geometry class. 

    // access to module geometries
    Int_t  NofGeometries() const;
    AliMUONGeometryModule* Geometry(Int_t i) const;
                  // In difference from protected GetGeometry()
		  // this function access geometry via index and not
		  // via moduleId


  protected:
    AliMUONVGeometryBuilder(const AliMUONVGeometryBuilder& rhs);

    // operators  
    AliMUONVGeometryBuilder& operator = (const AliMUONVGeometryBuilder& rhs);

    // methods
    AliMUONGeometryModule*         GetGeometry(Int_t moduleId) const;
    AliMUONGeometryEnvelopeStore*  GetEnvelopes(Int_t moduleId) const;
    AliMUONGeometryStore*          GetDetElements(Int_t moduleId) const;
    AliMUONGeometrySVMap*          GetSVMap(Int_t moduleId) const;
    
  private:
    //methods
    TString  ComposePath(const TString& volName, Int_t copyNo) const; 
    void     MapSV(const TString&path, const TString& volName, 
                  Int_t detElemId) const;

    void FillData(Int_t moduleId, Int_t nofDetElements, 
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6) const; 
    void FillData(Int_t id, const TString& volName, Int_t copyNo,
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6) const;
    void FillData(const TString& sensVolumePath, Int_t detElemId) const;		   

    TString ReadData1(ifstream& in) const;
    TString ReadData2(ifstream& in) const;
    TString ReadData3(ifstream& in) const;

    void WriteTransform(ofstream& out, const TGeoCombiTrans* transform) const;
    void WriteData1(ofstream& out) const;
    void WriteData2(ofstream& out) const;
    void WriteData3(ofstream& out) const;

    // static data members
    static const TString fgkTransformFileNamePrefix; // the prefix for the name 
                                                 // of file with transformations
    static const TString fgkSVMapFileNamePrefix; // the prefix for the name of file 
                                                 // with sensitive volume map
    static const TString fgkOutFileNameSuffix;   // the suffix for the name of 
                                                 // generated files
    
    // data members
    TString     fTransformFileName; // the name file with transformations 
    TString     fSVMapFileName;     // the name file with sensitive volume map 
    TObjArray*  fModuleGeometries;  // the modules geometries that will be built
                                    // by this builder
    
  ClassDef(AliMUONVGeometryBuilder,3) // MUON chamber geometry base class
};

// inline functions

inline Int_t  AliMUONVGeometryBuilder::NofGeometries() const
{ return fModuleGeometries->GetEntriesFast(); }

inline AliMUONGeometryModule* AliMUONVGeometryBuilder::Geometry(Int_t i) const
{ return (AliMUONGeometryModule*)fModuleGeometries->At(i); }

#endif //ALI_MUON_V_GEOMETRY_BUILDER_H
