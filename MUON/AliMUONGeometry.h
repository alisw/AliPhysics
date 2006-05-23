/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup geometry
/// \class AliMUONGeometry
/// \brief Container class for geometry modules
///
/// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_H
#define ALI_MUON_GEOMETRY_H

#include <TObject.h>
#include <TGeoMatrix.h>

class AliMUONGeometryModule;
class AliMUONGeometryTransformer;

class TObjArray;

class AliMUONGeometry : public TObject
{
  public:
    AliMUONGeometry(Bool_t isOwner);
    AliMUONGeometry();
    virtual  ~AliMUONGeometry();
    
    // methods
    void    AddModule(AliMUONGeometryModule* module);

    Bool_t  ReadSVMap(const TString& fileName);
    Bool_t  WriteSVMap(const TString& fileName) const;

    // get methods
    const AliMUONGeometryModule* GetModule(
                                    Int_t index, Bool_t warn = true) const;

    const AliMUONGeometryModule* GetModuleByDEId(
                                    Int_t detElemId, Bool_t warn = true) const;

    AliMUONGeometryTransformer* GetTransformer() const;


  protected:
    AliMUONGeometry(const AliMUONGeometry& right);
    AliMUONGeometry&  operator = (const AliMUONGeometry& right);
 
  private:
    //methods
    TString  ComposePath(const TString& volName, Int_t copyNo) const; 

    void    FillData3(const TString& sensVolumePath, Int_t detElemId);		   
    TString ReadData3(ifstream& in);
    void    WriteData3(ofstream& out) const;

    // data members
    TObjArray*                  fModules;     ///< Array of geometry modules
    AliMUONGeometryTransformer* fTransformer; ///< Geometry transformer

  ClassDef(AliMUONGeometry,1)  // Geometry parametrisation
};

inline AliMUONGeometryTransformer* AliMUONGeometry::GetTransformer() const
{ return fTransformer; }

#endif //ALI_MUON_GEOMETRY_H







