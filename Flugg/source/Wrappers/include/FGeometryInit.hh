// $Id$
// Flugg tag $Name$

// modified 10/IX/99 for including preStepPoint
// modified 28/IX/99 for delating allocated memory
// modified 4/X/99 function FreeMemory
// modified 2/III/00 base class G4TransportationManager included   
// modified 20/III/00 PrintHistories() included

 
#ifndef FGeometryInit_h
#define FGeometryInit_h 1

#include "g4std/fstream"
#include "g4std/iomanip"
#include "globals.hh"

#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Isotope.hh"
#include "G4Navigator.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4TouchableHistory.hh"
#include "G4GeometryManager.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"


class FGeometryInit : public G4TransportationManager
{
   public:
	~FGeometryInit();	//destructor
	static FGeometryInit *GetInstance();
   	G4Navigator *getNavigatorForTracking();
        G4FieldManager * getFieldManager();
	void setDetConstruction(G4VUserDetectorConstruction* detector);
	void setDetector();
	void setMotherVolume();
	void createFlukaMatFile();
        void closeGeometry();

        void PrintHistories();
	void InitHistories();
        void DeleteHistories();
        void UpdateHistories(const G4NavigationHistory *, const G4int &);
	G4TouchableHistory * GetTouchableHistory();
	G4TouchableHistory * GetOldNavHist();
	G4TouchableHistory * GetTempNavHist();

	void InitHistArray();
	void DelHistArray();
	G4int * GetHistArray();

	void InitJrLtGeantArray();
	G4int * GetJrLtGeantArray();
        G4int GetLttcFlagGeant();
        void SetLttcFlagGeant(G4int);
        void PrintJrLtGeant(); 

   private:    
	FGeometryInit();	//costructor
 	G4VUserDetectorConstruction * fDetector;
	G4FieldManager * fFieldManager; 
        G4TransportationManager * fTransportationManager;
	G4Navigator *fNavigatorForTracking;     
	static FGeometryInit *flagInstance;
        G4VPhysicalVolume * myTopNode;
	G4GeometryManager * ptrGeoMan;
	G4int * ptrArray;
        G4TouchableHistory * ptrTouchHist;
 	G4TouchableHistory * ptrOldNavHist;
 	G4TouchableHistory * ptrTempNavHist;
        G4int * ptrJrLtGeant;
        G4int flagLttcGeant;
};

#endif  
