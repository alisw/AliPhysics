// $Id$
// Category: geometry
// by I. Hrivnacova, 27.07.2000 
//
// XML convertor that converts G4 basic geometry objects 
// to XML defined by AGDD.dtd
// (ATLAS Generic Detector Description)

#ifndef TG4_XML_CONVERTOR_H
#define TG4_XML_CONVERTOR_H

#include "TG4VXMLConvertor.h"
#include "TG4Globals.h"

#include <globals.hh>
#include <g4std/fstream>

class G4Material;
class G4VSolid;
class G4LogicalVolume;
class G4PVReplica;
class G4Box;
class G4Tubs;
class G4Cons;
class G4Trd;
class G4Trap;
class G4Polycone;
class G4Polyhedra;

class TG4XMLConvertor : public TG4VXMLConvertor
{
  public:
    TG4XMLConvertor(G4std::ofstream& outFile);
    virtual ~TG4XMLConvertor();

    // methods
    virtual void OpenMaterials(const G4String& version, const G4String& date, 
            const G4String& author, const G4String dtdVersion);
    virtual void OpenSection(const G4String& name, const G4String& version,
	    const G4String& date, const G4String& author,
            const G4String& topVolume);
    virtual void OpenComposition(const G4String& name);
    virtual void CloseMaterials();
    virtual void CloseSection();
    virtual void CloseComposition();

    virtual void WriteMaterial(const G4Material* material); 
    virtual void WriteSolid(G4String lvName, const G4VSolid* solid, 
                            G4String materialName); 
    virtual void WriteRotation(const G4RotationMatrix* rotation); 
    virtual void WritePosition(G4String lvName, G4ThreeVector position); 
    virtual void WritePositionWithRotation(
                               G4String lvName, G4ThreeVector position,
   			       const G4RotationMatrix* rotation); 
    virtual void WriteReplica(G4String lvName, G4PVReplica* pvr);			       
    virtual void WriteEmptyLine();
    virtual void IncreaseIndention();
    virtual void DecreaseIndention();

  private:
    //methods
    void CutName(G4String& name) const;
    void CutName(G4String& name, G4int size) const;
    void PutName(G4String& element, G4String name, G4String templ) const;
    
         // writing solids
    void WriteBox (G4String lvName, const G4Box*  box,  G4String materialName); 
    void WriteTubs(G4String lvName, const G4Tubs* tubs, G4String materialName); 
    void WriteCons(G4String lvName, const G4Cons* cons, G4String materialName); 
    void WriteTrd (G4String lvName, const G4Trd*  trd,  G4String materialName); 
    void WriteTrap(G4String lvName, const G4Trap* trap, G4String materialName); 
    void WritePolycone(G4String lvName, const G4Polycone* polycone, 
                   G4String materialName); 
    void WritePolyhedra(G4String lvName, const G4Polyhedra* polyhedra, 
                   G4String materialName); 
  
    // static data members
    static const G4int fgkMaxVolumeNameLength;
    static const G4int fgkMaxMaterialNameLength;

    // data members
    G4std::ofstream&  fOutFile;          //output file
    const G4String    fBasicIndention;   //basic indention 
    G4String          fIndention;        //indention string
    G4int             fRotationCounter;  //counter of rotations
    TG4RotationMatrixVector  fRotations; // vector of rot matrices
};

#endif //TG4_XML_CONVERTOR_H

