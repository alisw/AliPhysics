
// Flugg tag 
// modified 17/06/02: by I. Gonzalez. STL migration

//#include <stdio.h>
//#include <iomanip.h>
#include "FGeometryInit.hh"
#include "FluggNavigator.hh"
#include "WrapUtils.hh"
#include "FlukaMaterial.hh"
#include "FlukaCompound.hh"

FGeometryInit * FGeometryInit::flagInstance=0;

FGeometryInit* FGeometryInit::GetInstance()  {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg::FGeometryInit::GetInstance(), instance=" 
	 << flagInstance << G4endl;
#endif
  if (!flagInstance) 
    flagInstance = new FGeometryInit();
  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg::FGeometryInit::GetInstance(), instance=" 
	 << flagInstance << G4endl;
#endif
  return flagInstance;
}  


FGeometryInit::FGeometryInit():
  fDetector(0),
  fFieldManager(0),
  fTransportationManager(G4TransportationManager::GetTransportationManager()),
  myTopNode(0),
  ptrGeoMan(0),
  ptrArray(0),
  ptrTouchHist(0),
  ptrOldNavHist(0),
  ptrTempNavHist(0),
  ptrJrLtGeant(0)
{

#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::FGeometryInit()" << G4endl;
  G4cout << "\t+ Changing the G4Navigator for FluggNavigator..." << G4endl;
#endif
  G4Navigator* actualnav = fTransportationManager->GetNavigatorForTracking();
  if (actualnav) {
    FluggNavigator* newnav = new FluggNavigator();
    fTransportationManager->SetNavigatorForTracking(newnav);
  }
  else {
    G4cerr << "ERROR: Could not find the actual G4Navigator" << G4endl;
    abort();
  }

  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::FGeometryInit()" << G4endl;
#endif

}


FGeometryInit::~FGeometryInit() {
  G4cout << "==> Flugg FGeometryInit::~FGeometryInit()" << G4endl;
  DeleteHistories();
  ptrGeoMan->OpenGeometry();  
  if (fTransportationManager)
    delete fTransportationManager;
  if (ptrJrLtGeant)
    delete[] ptrJrLtGeant;
  DelHistArray();
  
  //keep ATTENTION: never delete a pointer twice!
  G4cout << "<== Flugg FGeometryInit::FGeometryInit()" << G4endl;
}

void FGeometryInit::Init() {
// Build and initialize G4 geometry
   setDetector();
   setMotherVolume(); 
   closeGeometry();
   InitHistories();   
   InitJrLtGeantArray(); 
   InitHistArray();
   createFlukaMatFile();
}   
   
     
void FGeometryInit::closeGeometry() {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::closeGeometry()" << G4endl;
#endif

  ptrGeoMan = G4GeometryManager::GetInstance();
  if (ptrGeoMan) {
    ptrGeoMan->OpenGeometry();
    
    //true argoment allows voxel construction; if false voxels are built 
    //only for replicated volumes  
    ptrGeoMan->CloseGeometry(true);
  }
  else {
    G4cerr << "ERROR in FLUGG: Could not get G4GeometryManager instance" 
	   << G4endl;
    G4cerr << "                in FGeometryInit::closeGeometry. Exiting!!!"
	   << G4endl;
  }

#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::closeGeometry()" << G4endl;
#endif
}

//*************************************************************************

void FGeometryInit::InitHistArray() {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::InitHistArray()" << G4endl;
#endif
  ptrArray = new G4int[1000000];
  for(G4int i=0;i<1000000;i++) 
    ptrArray[i]=0;
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::InitHistArray()" << G4endl;
#endif
}



//*************************************************************************
//jrLtGeant stores all crossed lattice volume histories.

void FGeometryInit::InitJrLtGeantArray() {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::InitJrLtGeantArray()" << G4endl;
  G4cout << "Initializing JrLtGeant array" << G4endl;
#endif
  ptrJrLtGeant = new G4int[10000];
  for(G4int x=0;x<10000;x++) 
    ptrJrLtGeant[x]=-1;
  flagLttcGeant = -1;
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::InitJrLtGeantArray()" << G4endl;
#endif
}


void FGeometryInit::SetLttcFlagGeant(G4int newFlagLttc) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::SetLttcFlagGeant()" << G4endl;
#endif
  // Added by A.Solodkov
  if (newFlagLttc >= 10000) {
    G4cout << "Problems in FGeometryInit::SetLttcFlagGeant" << G4endl;
    G4cout << "Index newFlagLttc=" << newFlagLttc << " is outside array bounds"
	   << G4endl;
    G4cout << "Better to stop immediately !" << G4endl;
    exit(1);
  }
  flagLttcGeant = newFlagLttc;
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::SetLttcFlagGeant()" << G4endl;
#endif
}

void FGeometryInit::PrintJrLtGeant() {
#ifdef G4GEOMETRY_DEBUG
  //G4cout << "jrLtGeant:" << G4endl;
  //for(G4int y=0;y<=flagLttcGeant;y++)
  //
  //	 G4cout << "jrLtGeant[" << y << "]=" << ptrJrLtGeant[y] << G4endl;
#endif
}

//**************************************************************************

void FGeometryInit::PrintHistories() {
  /*
    #ifdef G4GEOMETRY_DEBUG
    G4cout << "Touch hist:" << G4endl;
    G4cout << *(ptrTouchHist->GetHistory()) << G4endl;
    G4cout << "Tmp hist:" << G4endl;
    G4cout << *(ptrTempNavHist->GetHistory()) << G4endl;
    G4cout << "Old hist:" << G4endl;
    G4cout << *(ptrOldNavHist->GetHistory()) << G4endl;
    #endif
 */
}



void FGeometryInit::InitHistories() {  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::InitHistories()" << G4endl;
#endif
  //init utility histories with navigator history
  ptrTouchHist = 
    fTransportationManager->GetNavigatorForTracking()->CreateTouchableHistory();
  ptrTempNavHist = 
    fTransportationManager->GetNavigatorForTracking()->CreateTouchableHistory();   
  ptrOldNavHist = new G4TouchableHistory();
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::InitHistories()" << G4endl;
#endif
}

void FGeometryInit::DeleteHistories() {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::DeleteHistories()" << G4endl;
#endif

  delete ptrTouchHist;
  delete ptrOldNavHist;
  delete ptrTempNavHist;
  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "Deleting step-history objects at end of run!" << G4endl;
  G4cout << "<== Flugg FGeometryInit::DeleteHistories()" << G4endl;
#endif
}


void FGeometryInit::UpdateHistories(const G4NavigationHistory * history,
				    G4int flagHist) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::UpdateHistories()" << G4endl;
#endif
  PrintHistories();
  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "...updating histories!" << G4endl;
#endif
  
  G4VPhysicalVolume * pPhysVol = history->GetTopVolume();
  
  switch (flagHist) {
  case 0: {
    //this is the case when a new history is given to the 
    //navigator and old history has to be resetted
    //touchable history has not been updated jet, so:
    
    ptrTouchHist->UpdateYourself(pPhysVol,history);
    ptrTempNavHist->UpdateYourself(pPhysVol,history);
    G4NavigationHistory * ptrOldNavHistNotConst = 
      const_cast<G4NavigationHistory * >(ptrOldNavHist->GetHistory()); 
    ptrOldNavHistNotConst->Reset();
    ptrOldNavHistNotConst->Clear();
    PrintHistories();
    break; 
  } //case 0

  case 1: {
    //this is the case when a new history is given to the 
    //navigator but old history has to be kept (e.g. LOOKZ
    //is call during an event);
    //touchable history has not been updated jet, so:
    
    ptrTouchHist->UpdateYourself(pPhysVol,history);
    ptrTempNavHist->UpdateYourself(pPhysVol,history);
    PrintHistories();
    break;
  } //case 1

  case 2: {
    //this is the case when the touchable history has been 
    //updated by a LocateGlobalPointAndUpdateTouchable call
    
    G4VPhysicalVolume * pPhysVolTemp = ptrTempNavHist->GetVolume();
    ptrOldNavHist->UpdateYourself(pPhysVolTemp,
				  ptrTempNavHist->GetHistory());
    
    ptrTempNavHist->UpdateYourself(pPhysVol,history);
    PrintHistories();
    break;
  } //case 2

  default: {
    G4cout <<" ERROR in updating step-histories!" << G4endl;
    break;
  } //default
  } //switch
  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::UpdateHistories()" << G4endl;
#endif
}

//*****************************************************************************

void FGeometryInit::createFlukaMatFile() {
  // last modification Sara Vanini 1/III/99
  // NAMES OF ELEMENTS AND COMPOUNDS: the names must be written in upper case,
  // according to the fluka standard. In addition,. they must be equal to the
  // names of the fluka materials - see fluka manual - in order that the 
  // program load the right cross sections, and equal to the names included in
  // the .pemf. Otherwise the user must define the LOW-MAT CARDS, and make his
  // own .pemf, in order to get the right cross sections loaded in memory.

#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::createFlukaMatFile()" << G4endl;
  G4cout << "================== FILEWR =================" << G4endl;
#endif 


  //Regions map
  BuildRegionsMap();
  std::ofstream vos("Volumes_index.inp");
  PrintRegionsMap(vos);
  vos.close();

  //Materials and compounds
  BuildMaterialTables();
  std::ofstream fos("flukaMat.inp");  
  PrintMaterialTables(fos);
  PrintAssignmat(fos);
  PrintMagneticField(fos);
  fos.close();

#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::createFlukaMatFile()" << G4endl;
#endif
}

////////////////////////////////////////////////////////////////////////
// 
void FGeometryInit::BuildRegionsMap() {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::BuildRegionsMap()" << G4endl;
#endif

  //Find number of Volumes in physical volume store
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  unsigned int numVol = pVolStore->size();

  G4cout << "\t* G4PhysicalVolumeStore (" << pVolStore 
	 << ") has " << numVol << " volumes. Iterating..." 
	 << G4endl;

  for(unsigned int l=0; l < numVol; l++) {
    //Get each of the physical volumes
    G4VPhysicalVolume * physicalvolume = (*pVolStore)[l];
    G4int iFlukaRegion = l+1;
    fRegionVolumeMap[physicalvolume] = iFlukaRegion;
  }



#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::BuildRegionsMap()" << G4endl;
#endif
}

void FGeometryInit::PrintRegionsMap(std::ostream& os) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::PrintRegionsMap()" << G4endl;
#endif

  //Print some header
  PrintHeader(os, "GEANT4 VOLUMES");

  //Iterate over all volumes in the map
  for (RegionIterator i = fRegionVolumeMap.begin(); 
       i != fRegionVolumeMap.end(); 
       i++) {
    
    //Get info in the map
    G4VPhysicalVolume* ptrVol = (*i).first;
    int index = (*i).second;

    //Print index and region name in some fixed format
    os.setf(std::ios::left, std::ios::adjustfield);
    os << setw10 << index;
    os << std::setw(20) << ptrVol->GetName() << std::setw(20) << "";

    //If volume is a replica... print some more stuff
    if(ptrVol->IsReplicated()) {
      EAxis axis;
      G4int nRep = -1;
      G4double width = -1;
      G4double offset = -1;
      G4bool consum = false;
      ptrVol->GetReplicationData(axis, nRep, width, offset, consum);
      os.setf(std::ios::left, std::ios::adjustfield);
      os << setw10 << "Repetion Nb: " << std::setw(3) << nRep;
    }
    os << G4endl;
    
  }
  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::PrintRegionsMap()" << G4endl;
#endif
}

////////////////////////////////////////////////////////////////////////
//
void    FGeometryInit::BuildMediaMap()
{
    fRegionMediumMap = new int[fNRegions+1];
    for (RegionIterator i = fRegionVolumeMap.begin(); 
	 i != fRegionVolumeMap.end(); 
	 i++) {
	//Get info in the map
	G4VPhysicalVolume* ptrVol = (*i).first;
	int region = (*i).second;
	G4int imed = fMediumVolumeMap[ptrVol];
	fRegionMediumMap[region] = imed;
	printf("BuildMediaMap %s %d %d\n",(ptrVol->GetName()).data(), region, imed);
	
    }
}

G4int FGeometryInit::GetMedium(int region) const
{
    return fRegionMediumMap[region];
}


void  FGeometryInit::SetMediumFromName(const char* volName, int medium, int volid) 
 {
  char name4[5];
  char tmp[5];
  strncpy(tmp, volName, 4);
  tmp[4]='\0';
  fNRegions = 0;
  
  for (RegionIterator i = fRegionVolumeMap.begin(); 
       i != fRegionVolumeMap.end(); 
       i++) {
    fNRegions++;
    //Get info in the map
    G4VPhysicalVolume* ptrVol = (*i).first;
    strncpy(name4, (ptrVol->GetName()).data(), 4);
    name4[4]='\0';
    for (int j = 0; j < 4; j++) {
	if (name4[j] == '\0') {
	    for (int k = j; k < 4; k++) {
		name4[k] = ' ';
	    }
	    break;
	}
    }
    if (! strncmp(name4, tmp, 4)) {
	fMediumVolumeMap[ptrVol] = medium;
	fVolIdVolumeMap[ptrVol]  = volid;
    }
  }
}



////////////////////////////////////////////////////////////////////////
// 
void FGeometryInit::BuildMaterialTables() {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::BuildMaterialTables()" << G4endl;
#endif

  //some terminal printout also
  G4cout << "\t* Storing information..." << G4endl;

  //The logic is the folloing:
  //Get the Material Table and:
  // 1) For materials with density <= 1.00e-10*g/cm3 assign vacuum
  // 2) For each single element material build a material equivalent
  // 3) For the rest:
  //   3.a) Build materials for each not already known element
  //   3.b) Build the compound out of them

  //Get the Material Table and iterate
  const G4MaterialTable* matTable = G4Material::GetMaterialTable();
  for (MatTableIterator i = matTable->begin(); i != matTable->end(); i++) {

    //Get some basic material information
    G4Material* material = (*i);
    G4String matName = material->GetName();
    const G4double matDensity = material->GetDensity();
    const G4int nMatElements  = material->GetNumberOfElements();

    G4cout << "\t\t+ " << matName 
	   << ": dens. = " << matDensity/(g/cm3) << "g/cm3"
	   << ", nElem = " << nMatElements << G4endl;

    // 1) For materials with density <= 1.00e-10*g/cm3 assign vacuum
    //    FlukaMaterial* is 0 in that case
    if (matDensity <= 1.00e-10*g/cm3) {
      G4FlukaMaterialMap[material] = 0;
      G4cout << "\t\t  Stored as vacuum" << G4endl;
    }
    // 2) For each single element material build a material equivalent
    else if (nMatElements == 1) {
      
      FlukaMaterial *flukamat = 
	BuildFlukaMaterialFromElement(material->GetElement(0),
				      matDensity);
      
      G4FlukaMaterialMap[material] = flukamat;
      G4cout << "\t\t  Stored as " << flukamat->GetRealName() << G4endl;
      
    } //else if (material->GetNumberOfElements() == 1)
    
    // 3) For the rest:
    //   3.a) Build materials for each not already known element
    //   3.b) Build the compound out of them
    else {
      FlukaCompound* flukacomp = 
	BuildFlukaCompoundFromMaterial(material);
      G4FlukaCompoundMap[material] = flukacomp;
      G4cout << "\t\t  Stored as " << flukacomp->GetRealName() << G4endl;
    } //else for case 3)
  } //for (materials)
  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::BuildMaterialTables()" << G4endl;
#endif
}

FlukaMaterial* 
FGeometryInit::BuildFlukaMaterialFromElement(const G4Element* element,
					     G4double matDensity) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::BuildFlukaMaterialFromElement()" 
	 << G4endl;
#endif

  //Get element and its properties
  G4String elemName(ToFlukaString(element->GetName()));

  FlukaMaterial* flukamat = FlukaMaterial::GetFlukaMaterial(elemName);
  if (matDensity != 0 || (matDensity == 0 && flukamat == 0)) {
    //Check for isotopes
    G4int nIsotopes = element->GetNumberOfIsotopes();
    if (nIsotopes == 0) {
      G4double elemA = element->GetA()/g;
      G4double elemZ = element->GetZ();
      
      if (elemA != G4int(elemA) && elemZ != G4int(elemZ)) {
	G4cout << "WARNING: Element \'" << elemName 
	       << "\' has non integer Z (" << elemZ << ") or A (" 
	       << elemA << ")"
	       << G4endl;
      }
    
      flukamat = new FlukaMaterial(elemName,
				   G4int(elemZ),
				   elemA,
				   matDensity/(g/cm3));
    }
    else if (nIsotopes == 1) {
      const G4Isotope* isotope = element->GetIsotope(0);
      flukamat = BuildFlukaMaterialFromIsotope(isotope, matDensity);
    }
    else {
      FlukaCompound *flucomp = BuildFlukaCompoundFromElement(element,
							     matDensity);
      flukamat = flucomp->GetFlukaMaterial();      
    }
  }
#ifdef G4GEOMETRY_DEBUG
  else {
    G4cout << "INFO: Element \'" << elemName 
	   << "\' already exists in the DB. It will not be recreated."
	   << G4endl;
  }
#endif

  return flukamat;
  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::BuildFlukaMaterialFromElement()" 
	 << G4endl;
#endif
}

FlukaMaterial* 
FGeometryInit::BuildFlukaMaterialFromIsotope(const G4Isotope* isotope,
					     G4double matDensity) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::BuildFlukaMaterialFromIsotope()" 
	 << G4endl;
#endif
  G4String isoName(ToFlukaString(isotope->GetName()));
  FlukaMaterial* flukamat = FlukaMaterial::GetFlukaMaterial(isoName);
  if (matDensity != 0 || (matDensity == 0 && flukamat == 0)) {
    G4int isoZ = isotope->GetZ();
    G4double isoA = (isotope->GetA())/(g);
    G4int isoN = isotope->GetN();
    flukamat = new FlukaMaterial(isoName,
				 isoZ,
				 isoA,
				 matDensity/(g/cm3),
				 isoN);
  }

  return flukamat;

#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::BuildFlukaMaterialFromIsotope()" 
	 << G4endl;
#endif
}

FlukaCompound* 
FGeometryInit::BuildFlukaCompoundFromMaterial(const G4Material* material) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::BuildFlukaCompoundFromMaterial()" 
	 << G4endl;
#endif
  //Material properties
  const G4double* elemFractions = material->GetFractionVector();
  const G4int nMatElements  = material->GetNumberOfElements();
  const G4double matDensity = material->GetDensity();
  G4String matName(ToFlukaString(material->GetName()));
  FlukaCompound* flukacomp = new FlukaCompound(matName, matDensity/(g/cm3),
					       nMatElements);
  for (G4int i = 0; i < nMatElements; i++) {
    FlukaMaterial *flukamat = 
      BuildFlukaMaterialFromElement(material->GetElement(i), 0.0);
    
    flukacomp->AddElement(flukamat->GetIndex(), -elemFractions[i]);
    
  } //for (elements)

  return flukacomp;

#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::BuildFlukaCompoundFromMaterial()" 
	 << G4endl;
#endif
}

FlukaCompound* 
FGeometryInit::BuildFlukaCompoundFromElement(const G4Element* element,
					     G4double matDensity) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::BuildFlukaCompoundFromElement()" 
	 << G4endl;
#endif
  G4int nIsotopes = element->GetNumberOfIsotopes();
  //fraction of nb of atomes per volume (= volume fraction?)
  const G4double* isoAbundance =  element->GetRelativeAbundanceVector();
  G4String elemName(ToFlukaString(element->GetName()));

  //Material properties
  FlukaCompound* flukacomp = new FlukaCompound(elemName, matDensity/(g/cm3),
					       nIsotopes);
  for (G4int i = 0; i < nIsotopes; i++) {
    FlukaMaterial *flukamat = 
      BuildFlukaMaterialFromIsotope(element->GetIsotope(i), 0.0);
    
    flukacomp->AddElement(flukamat->GetIndex(), isoAbundance[i]);
    
  } //for (elements)

  return flukacomp;

#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::BuildFlukaCompoundFromElement()" 
	 << G4endl;
#endif
}

void FGeometryInit::PrintMaterialTables(std::ostream& os) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::PrintMaterialTables()" << G4endl;
#endif
  //Print Header
  PrintHeader(os, "GEANT4 MATERIALS AND COMPOUNDS");
  
  //And some more stuff  
  size_t nIsotopes = G4Isotope::GetNumberOfIsotopes();
  size_t nElements = G4Element::GetNumberOfElements();
  size_t nMaterials = G4Material::GetNumberOfMaterials();

  os << "* In Geant4 there are " << nMaterials << " materials" << G4endl;
  os << "* In Geant4 there are " << nElements  << " elements"  << G4endl;
  os << "* In Geant4 there are " << nIsotopes  << " isotopes"  << G4endl;

  //Materials
  G4cout << "\t* Printing FLUKA materials..." << G4endl;
  FlukaMaterial::PrintMaterialsByIndex(os);
  //FlukaMaterial::PrintMaterialsByName(os);

  //Compounds
  G4cout << "\t* Printing FLUKA compounds..." << G4endl;
  FlukaCompound::PrintCompounds(os);

#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::PrintMaterialTables()" << G4endl;
#endif
}

////////////////////////////////////////////////////////////////////////
// 
void FGeometryInit::PrintAssignmat(std::ostream& os) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::PrintAssignmat()" << G4endl;
#endif

  //Find number of Volumes in physical volume store
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  unsigned int numVol = pVolStore->size();

  G4cout << "\t* G4PhysicalVolumeStore (" << pVolStore 
	 << ") has " << numVol << " volumes. " << G4endl;

  G4cout << "\t* Printing ASSIGNMAT..." << G4endl;


  PrintHeader(os,"GEANT4 MATERIAL ASSIGNMENTS");
  for(unsigned int l=0; l < numVol; l++) {

    //Get each of the physical volumes
    G4VPhysicalVolume * physicalvol = (*pVolStore)[l];

    //Get index for that volume
    G4int iFlukaRegion = fRegionVolumeMap[physicalvol];

    //Find G4 material and navigate to its fluka compound/material
    G4LogicalVolume * logicalVol = physicalvol->GetLogicalVolume();
    G4Material* material = logicalVol->GetMaterial();
    G4int matIndex = 2;
    if (G4FlukaCompoundMap[material])
      matIndex = G4FlukaCompoundMap[material]->GetIndex();
    if (G4FlukaMaterialMap[material])
      matIndex = G4FlukaMaterialMap[material]->GetIndex();

    //Find if there is a magnetic field in the region
    //check if Magnetic Field is present in the region
    G4double flagField = 0.0;
    G4FieldManager * pMagFieldMan = logicalVol->GetFieldManager();
    if(pMagFieldMan && pMagFieldMan->GetDetectorField())
      flagField = 1.0;
    
    //Print card
    os << setw10 << "ASSIGNMAT ";
    os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
    os << setw10 << setfixed << G4double(matIndex);
    os << setw10 << setfixed << G4double(iFlukaRegion);
    os << setw10 << "0.0";
    os << setw10 << setfixed << flagField;
    os << G4endl;
  }



#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::PrintAssignmat()" << G4endl;
#endif
}


void FGeometryInit::PrintMagneticField(std::ostream& os) {
#ifdef G4GEOMETRY_DEBUG
  G4cout << "==> Flugg FGeometryInit::PrintMagneticField()" << G4endl;
#endif

  G4cout << "\t* Printing Magnetic Field..." << G4endl;

  if(fTransportationManager->GetFieldManager()->DoesFieldExist()) {
    
    //get magnetic field pointer
    const G4Field * pMagField = 
      fTransportationManager->GetFieldManager()->GetDetectorField();     
    
    
    if(pMagField) {
      //Check if it can be made a uniform magnetic field
      const G4UniformMagField *pUnifMagField = 
	dynamic_cast<const G4UniformMagField*>(pMagField);
      if(pUnifMagField) {
	G4double B[3];
	G4double point[4]; //it is not really used
	pUnifMagField->GetFieldValue(point,B);

	//write MGNFIELD card 
	PrintHeader(os,"GEANT4 MAGNETIC FIELD");
	os << setw10 << "MGNFIELD  ";
	os << setw10 << "";
	os << setw10 << "";
	os << setw10 << "";
	os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
	os << setw10 << setfixed
	   << std::setprecision(4) << B[0]
	   << setw10 << B[1]
	   << setw10 << B[2]
	   << G4endl;
      }
      else {
	G4cout << "WARNING: No Uniform Magnetic Field found." << G4endl;
	G4cout << "         Manual intervention might be needed." << G4endl;
      }
    }
    else
      G4cout << "\t  No detector field found... " << G4endl;
  } // end if magnetic field
  else
    G4cout << "\t  No field found... " << G4endl;

#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::PrintMagneticField()" << G4endl;
#endif
}

int FGeometryInit::CurrentVolID(int ir, int& copyNo)
{
    if (ir == 0) 
    {
	copyNo = -1;
	return -1;
    }
    
    G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
    G4VPhysicalVolume   * physicalvol = (*pVolStore)[ir- 1];
    
    if (physicalvol) {
	copyNo =  physicalvol->GetCopyNo();
    } else {
	copyNo = -1;
	return -1;
    }
    
    
    int id = fVolIdVolumeMap[physicalvol];
    return id;
}

int FGeometryInit::CurrentVolOffID(int ir, int off, int& copyNo)
{
    if (off == 0) return CurrentVolID(ir, copyNo);
    
    G4PhysicalVolumeStore* pVolStore = G4PhysicalVolumeStore::GetInstance();
    G4VPhysicalVolume*     physicalvol = (*pVolStore)[ir- 1];
    G4VPhysicalVolume*     mother = physicalvol; 

    int index;
//============================================================================
    if (mother) {
  // Check touchable depth
  //
       if (ptrTouchHist->GetHistoryDepth() < off) {
          mother = 0;
       } else {                                                                                                                                                             
          // Get the off-th mother
          index = ptrTouchHist->GetHistoryDepth() - off;
          // in the touchable history volumes are ordered
          // from top volume up to mother volume;
          // the touchable volume is not in the history	                                                                                
          mother = ptrTouchHist->GetHistory()->GetVolume(index);
       }
    }	
//============================================================================
    
    int id;
    
    if (!mother) {
	G4cout << "Flugg FGeometryInit::CurrentVolOffID mother not found" << G4endl;
	id = -1;
	copyNo = -1;
    } else {
	copyNo =  mother ->GetCopyNo();
	id =  fVolIdVolumeMap[mother];
    }
    return id;
}

void  FGeometryInit::Gmtod(double* xm, double* xd, int iflag)
{
// Transforms a position from the world reference frame
// to the current volume reference frame.
//
//  Geant3 desription:
//  ==================
//       Computes coordinates XD (in DRS) 
//       from known coordinates XM in MRS 
//       The local reference system can be initialized by
//         - the tracking routines and GMTOD used in GUSTEP
//         - a call to GMEDIA(XM,NUMED)
//         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER) 
//             (inverse routine is GDTOM) 
//
//        If IFLAG=1  convert coordinates 
//           IFLAG=2  convert direction cosinus
//
// ---
    FluggNavigator        * ptrNavig     = getNavigatorForTracking();
    //setting variables (and dimension: Fluka uses cm.!)
    G4ThreeVector pGlob(xm[0],xm[1],xm[2]);
        G4ThreeVector pLoc;
    
    if (iflag == 1) {
	pGlob *= 10.0; // in mm
// change because of geant 4 6.0
//	pLoc = ptrNavig->ComputeLocalPoint(pGlob);
        pLoc = ptrNavig->GetGlobalToLocalTransform().TransformPoint(pGlob);

	pLoc /= 10.0;  // in cm
    } else if (iflag == 2) {
	pLoc = 
	    ptrNavig->ComputeLocalAxis(pGlob);	
    } else {
	G4cout << "Flugg FGeometryInit::Gmtod called with undefined flag" << G4endl;
    }
    xd[0] = pLoc[0]; xd[1] = pLoc[1]; xd[2] = pLoc[2];
}

void  FGeometryInit::Gdtom(double* xd, double* xm, int iflag)
{
// Transforms a position from the current volume reference frame
// to the world reference frame.
//
//  Geant3 desription:
//  ==================
//  Computes coordinates XM (Master Reference System
//  knowing the coordinates XD (Detector Ref System)
//  The local reference system can be initialized by
//    - the tracking routines and GDTOM used in GUSTEP
//    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
//        (inverse routine is GMTOD)
// 
//   If IFLAG=1  convert coordinates
//      IFLAG=2  convert direction cosinus
//
// ---

    FluggNavigator        * ptrNavig     = getNavigatorForTracking();
    G4ThreeVector pLoc(xd[0],xd[1],xd[2]);
        
    G4ThreeVector pGlob;
     if (iflag == 1) {
	 pLoc  *= 10.0; // in mm
	 pGlob = ptrNavig->GetLocalToGlobalTransform().
	     TransformPoint(pLoc);
	 pGlob /= 10.0; // in cm
     } else if (iflag == 2) {
	 pGlob = ptrNavig->GetLocalToGlobalTransform().
	     TransformAxis(pLoc);
     } else {
	 G4cout << "Flugg FGeometryInit::Gdtom called with undefined flag" << G4endl;
     }
     
     xm[0] = pGlob[0]; xm[1] = pGlob[1]; xm[2] = pGlob[2];
}
