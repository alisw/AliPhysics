// $Id$
// Category: geometry
// by I. Hrivnacova, 27.07.2000 
//
// Singleton class for generation of geometry data files in XML,
// the XML format is independent from G4 geometry
// object model. 

#ifndef TG4_XML_GEOMETRY_GENERATOR_H
#define TG4_XML_GEOMETRY_GENERATOR_H

#include <globals.hh>
#include <g4std/fstream>

class TG4VXMLConvertor;

class G4LogicalVolume;

class TG4XMLGeometryGenerator
{
  public:
    TG4XMLGeometryGenerator();
    // --> protected
    // TG4XMLGeometryGenerator(const TG4XMLGeometryGenerator& right);
    virtual ~TG4XMLGeometryGenerator();

    // static access method
    static TG4XMLGeometryGenerator* Instance();

    // methods
    void GenerateMaterials(const G4String& version, const G4String& date,
            const G4String& author,  const G4String dtdVersion,
	    G4LogicalVolume* lv);
    void GenerateSection(const G4String& name, const G4String& version,
	    const G4String& date, const G4String& author,
            const G4String& topVolume, G4LogicalVolume* lv);
    void OpenFile(G4String filePath);
    void CloseFile();

  protected:
    TG4XMLGeometryGenerator(const TG4XMLGeometryGenerator& right);

    // operators
    TG4XMLGeometryGenerator& operator=(const TG4XMLGeometryGenerator& right);

  private:
    // methods
    void ProcessLogicalVolume(G4LogicalVolume* lv); 
    void ProcessMaterials(G4LogicalVolume* lv); 
    void ProcessSolids(G4LogicalVolume* lv); 
    void ProcessRotations(G4LogicalVolume* lv); 

    // static data members
    static TG4XMLGeometryGenerator*  fgInstance;     //this instance

    // data members
    TG4VXMLConvertor*  fConvertor; //interface to XML convertor 
    G4std::ofstream    fOutFile;   //output file
};


// inline methods
inline TG4XMLGeometryGenerator* TG4XMLGeometryGenerator::Instance()
{ return fgInstance; }

#endif //TG4_GEOMETRY_XML_MANAGER_H

