// $Id$
// Category: geometry
// by I. Hrivnacova, 27.07.2000 
//
// The interface for XML convertor that
// converts G4 basic geometry objects to XML. 

#ifndef TG4_V_XML_CONVERTOR_H
#define TG4_V_XML_CONVERTOR_H

#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>

class G4Material;
class G4VSolid;
class G4LogicalVolume;

class TG4VXMLConvertor
{
  public:
    TG4VXMLConvertor();
    virtual ~TG4VXMLConvertor();

    // methods
    virtual void OpenMaterials(const G4String& version, const G4String& date, 
            const G4String& author, const G4String dtdVersion) = 0; 
    virtual void OpenSection(const G4String& name, const G4String& version,
	    const G4String& date, const G4String& author,
            const G4String& topVolume) = 0;
    virtual void OpenComposition(const G4String& name) = 0;
    virtual void CloseMaterials() = 0;
    virtual void CloseSection() = 0;
    virtual void CloseComposition() = 0;

    virtual void WriteMaterial(const G4Material* material) = 0; 
    virtual void WriteSolid(const G4VSolid* solid, G4String materialName) = 0; 
    virtual void WriteRotation(const G4RotationMatrix* rotation) = 0; 
    virtual void WritePosition(G4String solidName, G4ThreeVector position) = 0; 
    virtual void WritePositionWithRotation(
                               G4String solidName, G4ThreeVector position,
                               const G4RotationMatrix* rotation) = 0; 
    virtual void WriteEmptyLine() = 0;
    virtual void IncreaseIndention() = 0;
    virtual void DecreaseIndention() = 0;
};

#endif //TG4_V_XML_CONVERTOR_H

