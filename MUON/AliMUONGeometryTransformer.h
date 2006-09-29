/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometryTransformer
/// \brief Top container class for geometry transformations
///
/// Geometry transformations can be filled in these ways:
/// - by geometry builder when geometry is built via builders
/// - from Root geometry file (*.root) or Root geometry manager
/// - from ASCII file (*.dat)
/// If geometry is loaded from a file, the list of aligned volume paths
/// has to be read first from volpaths.dat file.
/// 
/// \author Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_TRANSFORMER_H
#define ALI_MUON_GEOMETRY_TRANSFORMER_H

#include <TObject.h>
#include <TObjArray.h>
#include <TGeoMatrix.h>

class AliMUONGeometryModuleTransformer;
class AliMUONGeometryDetElement;

class TGeoManager;
class TClonesArray;

class AliMUONGeometryTransformer : public TObject
{
  public:
    AliMUONGeometryTransformer(Bool_t isOwner, const TString& detName = "MUON");
    AliMUONGeometryTransformer();
    virtual  ~AliMUONGeometryTransformer();
    
    // methods
    void  AddModuleTransformer(AliMUONGeometryModuleTransformer* transformer);
    void  AddMisAlignModule(Int_t moduleId, const TGeoHMatrix& matrix);
    void  AddMisAlignDetElement(Int_t detElemId, const TGeoHMatrix& matrix);

    void  AddAlignableVolumes() const; 
    TClonesArray* CreateZeroAlignmentData() const;			       

    // IO
    //
    Bool_t  ReadGeometryData(const TString& volPathFileName,
                             const TString& transformFileName);
    Bool_t  ReadGeometryData(const TString& volPathFileName,
                             TGeoManager* geoManager);

    Bool_t  WriteGeometryData(const TString& volPathFileName,
                             const TString& transformFileName,
			     const TString& misalignFileName = "") const;
   
    Bool_t  WriteVolumePaths(const TString& fileName) const;
    Bool_t  WriteTransformations(const TString& fileName) const;
    Bool_t  WriteMisAlignmentData(const TString& fileName) const;

    // Transformation methods 
    //
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

    // Get methods
    //
    Int_t GetNofModuleTransformers() const;
    const AliMUONGeometryModuleTransformer* GetModuleTransformer(
                               Int_t index, Bool_t warn = true) const;

    const AliMUONGeometryModuleTransformer* GetModuleTransformerByDEId(
                               Int_t detElemId, Bool_t warn = true) const;

    const AliMUONGeometryDetElement* GetDetElement(
                               Int_t detElemId, Bool_t warn = true) const;

    const TClonesArray* GetMisAlignmentData() const;
    
    Bool_t  HasDE(Int_t detElemId) const;

  protected:
    AliMUONGeometryTransformer(const AliMUONGeometryTransformer& right);
    AliMUONGeometryTransformer&  operator = (const AliMUONGeometryTransformer& right);
 
  private:
    // methods
    AliMUONGeometryModuleTransformer* GetModuleTransformerNonConst(
                                    Int_t index, Bool_t warn = true) const;

    TGeoHMatrix GetTransform(
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6) const;

    void FillModuleVolPath(Int_t moduleId, const TString& volPath); 
    void FillDetElemVolPath(Int_t detElemId, const TString& volPath); 

    void FillModuleTransform(Int_t moduleId,
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6); 
    void FillDetElemTransform(Int_t id, 
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6);

    Bool_t  ReadVolPaths(ifstream& in);
    TString ReadModuleTransforms(ifstream& in);
    TString ReadDetElemTransforms(ifstream& in);
    Bool_t  LoadTransforms(TGeoManager* tgeoManager); 

    Bool_t  ReadVolPaths(const TString& fileName);
    Bool_t  ReadTransformations(const TString& fileName);
    Bool_t  ReadTransformations2(const TString& fileName);

    void    WriteTransform(ofstream& out, const TGeoMatrix* transform) const;
    void    WriteModuleVolPaths(ofstream& out) const;
    void    WriteDetElemVolPaths(ofstream& out) const;
    void    WriteModuleTransforms(ofstream& out) const;
    void    WriteDetElemTransforms(ofstream& out) const;
    
    TString GetModuleSymName(Int_t moduleId) const;
    TString GetDESymName(Int_t detElemId) const;

    // data members
    TString        fDetectorName;       /// < Detector name
    TObjArray*     fModuleTransformers; ///< array of module transformers
    TClonesArray*  fMisAlignArray;      ///< array of misalignment data

  ClassDef(AliMUONGeometryTransformer,2)  // Geometry parametrisation
};

// inline methods

/// Return the number of contained module transformers
inline Int_t AliMUONGeometryTransformer::GetNofModuleTransformers() const
{ return fModuleTransformers->GetEntriesFast(); }

/// Return the array of misalignment data
inline const TClonesArray* AliMUONGeometryTransformer::GetMisAlignmentData() const	
{ return fMisAlignArray; }		       

#endif //ALI_MUON_GEOMETRY_TRANSFORMER_H







