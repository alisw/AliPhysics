/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// Revision of includes 07/05/2004
//
// Class AliMUONVGeometryBuilder
// -----------------------------
// Abstract base class for geometry construction per chamber(s).
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_V_GEOMETRY_BUILDER_H
#define ALI_MUON_V_GEOMETRY_BUILDER_H

#include <TObject.h>

class TGeoTranslation;
class TGeoRotation;
class TGeoCombiTrans;
class TObjArray;

class AliMUONChamber;
class AliMUONChamberGeometry;

class AliMUONVGeometryBuilder : public TObject
{
  public:
    AliMUONVGeometryBuilder(AliMUONChamber* ch1,
                            AliMUONChamber* ch2 = 0,
                            AliMUONChamber* ch3 = 0,
                            AliMUONChamber* ch4 = 0,
                            AliMUONChamber* ch5 = 0,
                            AliMUONChamber* ch6 = 0);
    AliMUONVGeometryBuilder();
    virtual ~AliMUONVGeometryBuilder();
  
    // methods
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
		   
  protected:
    AliMUONVGeometryBuilder(const AliMUONVGeometryBuilder& rhs);

    // operators  
    AliMUONVGeometryBuilder& operator = (const AliMUONVGeometryBuilder& rhs);

    // methods
    AliMUONChamber* GetChamber(Int_t chamberId) const;
    
  private:
    // data members
    TObjArray*  fChambers; // the chambers which geometry will be built
                           // by this builder
    
  ClassDef(AliMUONVGeometryBuilder,1) // MUON chamber geometry base class
};

#endif //ALI_MUON_V_GEOMETRY_BUILDER_H
