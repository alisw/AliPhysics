/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONGeometrySVMap
// --------------------------
// As the detection element frame is different from the
// frame of the sensitive volume(s) defined in Geant,
// the sensitive volumes have to be mapped to the detection 
// elements. In the map, fSVMap, the sensitive voolumes are specified
// by the full path in the volume hierarchy, defined as:
//  /volname.copyNo/volName.copyNo1/...
//
// The array of global positions of sensitive volumes fSVPositions
// is included to make easier the verification of the assignements 
// in the fSVMap.
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_SV_MAP_H
#define ALI_MUON_GEOMETRY_SV_MAP_H

#include <Riostream.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TArrayI.h>

class TGeoCombiTrans;
class TGeoTranslation;

// Substitutes map <string, int>
// which ALICE does not allow to use 

class AliMUONStringIntMap : public TObject
{
  public:
    AliMUONStringIntMap();
    virtual ~AliMUONStringIntMap();
    
    // methods
    Bool_t  Add(const TString& first, Int_t second);
    Int_t   Get(const TString& first) const;
    Int_t   GetNofItems() const;
    void    Clear();
    virtual void Print(const char* /*option*/ = "") const;
    void Print(const TString& key, ofstream& out) const;
    
  protected:
    AliMUONStringIntMap(const AliMUONStringIntMap& rhs);

    // operators  
    AliMUONStringIntMap& operator = (const AliMUONStringIntMap& rhs);
 
  private:
    // data members
    Int_t      fNofItems;    // number of items
    TObjArray  fFirstArray;  // first item array
    TArrayI    fSecondArray; // second item array
 
  ClassDef(AliMUONStringIntMap,1)  // motif map
};    


class AliMUONGeometrySVMap : public TObject
{
  public:
    AliMUONGeometrySVMap(Int_t initSize);
    AliMUONGeometrySVMap();
    virtual ~AliMUONGeometrySVMap();

    // methods
    void Add(const TString& volumePath, 
             Int_t detElemId);  
    void AddPosition(const TString& volumePath, 
             const TGeoTranslation& globalPosition);  

    void Clear();
    void ClearPositions();
    void SortPositions();
    virtual void Print(Option_t* option) const;
    void PrintPositions() const;
    void WriteMap(ofstream& out) const;
    
    // get methods
    Int_t  GetDetElemId(const TString& volumePath) const;

  protected:
    AliMUONGeometrySVMap(const AliMUONGeometrySVMap& rhs);

    // operators  
    AliMUONGeometrySVMap& operator 
      = (const AliMUONGeometrySVMap& rhs);
  
  private:
    const TGeoCombiTrans* FindByName(const TString& name) const;

    // data members
    AliMUONStringIntMap  fSVMap;       // Map of sensitive volume paths
                                       // and detector element id 
    TObjArray            fSVPositions; // The array of transformations

  ClassDef(AliMUONGeometrySVMap,1) // MUON sensitive volume map
};

#endif //ALI_MUON_GEOMETRY_TRANSFORM_STORE_H
