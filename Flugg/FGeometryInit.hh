
// Flugg tag 

// modified 10/IX/99 for including preStepPoint
// modified 28/IX/99 for delating allocated memory
// modified 4/X/99 function FreeMemory
// modified 2/III/00 base class G4TransportationManager included   
// modified 20/III/00 PrintHistories() included

 
#ifndef FGeometryInit_h
#define FGeometryInit_h 1

//#include "g4std/fstream"
//#include "g4std/iomanip"
//#include "g4rw/cstring.h"
#include "globals.hh"

#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Isotope.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4TouchableHistory.hh"
#include "G4GeometryManager.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"


#include <map>

class FluggNavigator;
class FlukaMaterial;
class FlukaCompound;

class FGeometryInit : public G4TransportationManager {

public:
  ~FGeometryInit();	//destructor
  static FGeometryInit *GetInstance();
  inline FluggNavigator *getNavigatorForTracking();
  inline G4FieldManager * getFieldManager();
  inline void setDetConstruction(G4VUserDetectorConstruction* detector);
  inline void setDetector();
  inline void setMotherVolume();
  void createFlukaMatFile();
  void closeGeometry();
  
  void PrintHistories();
  void InitHistories();
  void DeleteHistories();
  void UpdateHistories(const G4NavigationHistory *, G4int);
  inline G4TouchableHistory * GetTouchableHistory();
  inline G4TouchableHistory * GetOldNavHist();
  inline G4TouchableHistory * GetTempNavHist();
  
  void InitHistArray();
  inline void DelHistArray();
  inline G4int * GetHistArray();
  
  void InitJrLtGeantArray();
  inline G4int * GetJrLtGeantArray();
  inline G4int GetLttcFlagGeant();
  void SetLttcFlagGeant(G4int);
  void PrintJrLtGeant();
  void Init();

  //Map access methods
  void  BuildMediaMap();
  void  SetMediumFromName(const char* volName, int med, int volid);
    //G4int GetRegionFromName(const char* volName) const;
  G4int GetMedium(int) const;
  int CurrentVolID(int ir, int& copyNo);
  int CurrentVolOffID(int ir, int off, int& copyNo);
  void Gmtod(double* xm, double* xd, int iflag);
  void Gdtom(double* xd, double* xm, int iflag);    

protected:
  void BuildRegionsMap();
  void PrintRegionsMap(std::ostream& os);
  void BuildMaterialTables();
  FlukaMaterial* BuildFlukaMaterialFromElement(const G4Element* element, 
					       G4double matDensity);
  FlukaMaterial* BuildFlukaMaterialFromIsotope(const G4Isotope* isotope,
					       G4double matDensity);
  FlukaCompound* BuildFlukaCompoundFromMaterial(const G4Material* material); 
  FlukaCompound* BuildFlukaCompoundFromElement(const G4Element* element,
					       G4double matDensity);
  void PrintMaterialTables(std::ostream& os);
  void PrintAssignmat(std::ostream& os);
  void PrintMagneticField(std::ostream& os);

private:    
  FGeometryInit();	//costructor

private:
  G4VUserDetectorConstruction * fDetector;
  G4FieldManager * fFieldManager; 
  G4TransportationManager * fTransportationManager;
  static FGeometryInit *flagInstance;
  G4VPhysicalVolume * myTopNode;
  G4GeometryManager * ptrGeoMan;
  G4int * ptrArray;
  G4TouchableHistory * ptrTouchHist;
  G4TouchableHistory * ptrOldNavHist;
  G4TouchableHistory * ptrTempNavHist;
  G4int * ptrJrLtGeant;
  G4int flagLttcGeant;
  G4int  fNRegions;
  int* fRegionMediumMap;
    
  std::map<G4VPhysicalVolume*, int, std::less<G4VPhysicalVolume*> > fRegionVolumeMap;
  std::map<G4VPhysicalVolume*, int, std::less<G4VPhysicalVolume*> > fMediumVolumeMap;
  std::map<G4VPhysicalVolume*, int, std::less<G4VPhysicalVolume*> > fVolIdVolumeMap;

  std::map<G4Material*, FlukaMaterial*, std::less<G4Material*> > G4FlukaMaterialMap;
  std::map<G4Material*, FlukaCompound*, std::less<G4Material*> > G4FlukaCompoundMap;
  //G4int NOfMaterials;
};

typedef  std::map<G4VPhysicalVolume*, int, std::less<G4VPhysicalVolume*> >::const_iterator RegionIterator;
typedef  std::vector<G4Material*>::const_iterator MatTableIterator;


//Include the file with the inline methods
#include "FGeometryInit.icc"

#endif  
