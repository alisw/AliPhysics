
// Flugg tag 
// modified 17/06/02: by I. Gonzalez. STL migration

//#include <stdio.h>
//#include <iomanip.h>
#include "FGeometryInit.hh"
#include "FluggNavigator.hh"
#include "WrapUtils.hh"

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
  
  //open file for output
  G4std::ofstream fout("flukaMat.inp");  
  
  //PhysicalVolumeStore, Volume and MaterialTable pointers
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  unsigned int  numVol = pVolStore->size();
  
  G4int* indexMatFluka = 0; 
  static G4Material * ptrMat = 0;
  unsigned int x = 0;

  //Find a volume with a material associated
  while(!ptrMat && x<numVol) {
    G4VPhysicalVolume * ptrVol = (*pVolStore)[x];
    G4LogicalVolume * ptrLogVol = ptrVol->GetLogicalVolume();
    ptrMat = ptrLogVol->GetMaterial();
    x++;
  }
  
  if(ptrMat) {
    static const G4MaterialTable * ptrMatTab = G4Material::GetMaterialTable();
    
    //number of materials, elements, variable initialisations
    static size_t totNumMat = G4Material::GetNumberOfMaterials();
#ifdef G4GEOMETRY_DEBUG
    G4cout << "Number of materials: " << totNumMat << G4endl;
#endif 

    static size_t totNumElem = G4Element::GetNumberOfElements();
#ifdef G4GEOMETRY_DEBUG
    G4cout << "Number of elements: " << totNumMat << G4endl;
#endif 


    G4int * elemIndexInMATcard = new G4int[totNumElem];
    for(unsigned int t=0; t<totNumElem; t++) 
      elemIndexInMATcard[t] = 0;
    static const G4IsotopeTable * ptrIsotTab = 0;
    G4int initIsot = 0;
    static size_t totNumIsot;
    G4int* isotIndexInMATcard = 0;
    G4int isotPresence = 0; 
    
    
    // title
    PrintHeader(fout,"GEANT4 ELEMENTS AND COMPOUNDS");
    
    // *** loop over G4Materials to assign Fluka index 
    G4int indexCount = 3;
    indexMatFluka = new G4int[totNumMat];
    for(unsigned int i=0; i<totNumMat; i++) {
      //pointer material, state 
      ptrMat = (*ptrMatTab)[i];
      G4double denMat = ptrMat->GetDensity();
      G4String nameMat = ptrMat->GetName();
      // Fluka index: bh=1, vacuum=2, others=3,4..
      if(denMat<=1.00e-10*g/cm3)
	//N.B. fluka density limit decided on XI-`98  
	indexMatFluka[i]=2;	
      else {  
	indexMatFluka[i]=indexCount;
	indexCount++;
      }
      
      // *** write single-element material MATERIAL card
      size_t numElem = ptrMat->GetNumberOfElements();
      if(numElem==1) {
	G4int index = indexMatFluka[i];
	const G4Element * ptrElem = ptrMat->GetElement(0);
	size_t indElemTab = ptrElem->GetIndex();
	size_t numIsot = ptrElem->GetNumberOfIsotopes();
	G4double A = (ptrElem->GetA())/(g);
	if(!numIsot) {
	  if(index!=2 && !elemIndexInMATcard[indElemTab]) {
	    G4double Z = ptrElem->GetZ();
	    elemIndexInMATcard[indElemTab] = index;
	    G4String nameEl = ptrElem->GetName();
	    nameEl.toUpper();

	    //write on file MATERIAL card of element
	    PrintMaterial(fout, "MATERIAL  ",
			  Z, A,
			  denMat/(g/cm3),
			  index,
			  -1,
			  nameEl);
	  }
	}
	if(numIsot==1) {
	  // G4Isotope pointer 
	  const G4Isotope* ptrIsot = ptrElem->GetIsotope(0);
	  size_t indIsotTab = ptrIsot->GetIndex();
	  //initialize variables
	  if(!initIsot) {
	    totNumIsot = ptrIsot->GetNumberOfIsotopes();
	    ptrIsotTab = ptrIsot->GetIsotopeTable(); 
	    isotIndexInMATcard = new G4int[totNumIsot];
	    for(unsigned int t=0; t<totNumIsot; t++) 
	      isotIndexInMATcard[t] = 0;
	    initIsot = 1;
	  }
	  if(!isotIndexInMATcard[indIsotTab]) {
	    //compute physical data and counters
	    G4int ZIs = ptrIsot->GetZ();
	    G4double AIs = (ptrIsot->GetA())/(g);
	    G4int NIs = ptrIsot->GetN();
	    G4String nameIsot = ptrIsot->GetName();
	    nameIsot.toUpper();
	    isotIndexInMATcard[indIsotTab] = index;
	    
	    //write on file MATERIAL card of isotope
	    PrintMaterial(fout, "MATERIAL  ",
			  G4double(ZIs),
			  AIs,
			  denMat/(g/cm3),
			  G4double(index),
			  G4double(NIs),
			  nameIsot);
	  }
	}// end if(numIsot==1)
      }// end if(numElem==1)
    }// end for loop
    
    
    // *** material definitions: elements, compound made of G4Elements
    // or made of G4Materials
    for(unsigned int j=0; j<totNumMat; j++) {
      //pointer material, and material data 
      ptrMat = (*ptrMatTab)[j];
      size_t numElem = ptrMat->GetNumberOfElements();
      G4double densityMat = (ptrMat->GetDensity())/(g/cm3);
      G4String nameMat = ptrMat->GetName();
      nameMat.toUpper();
      isotPresence = 0;
      
      //fraction vector of compounds of the material
      const G4double* ptrFracVect = ptrMat->GetFractionVector();
      
      //loop on elements of the material
      for(unsigned int el=0; el<numElem; el++) {
	//compute physical data, initialize variables
	const G4Element * ptrElem = ptrMat->GetElement(el);
	size_t indElemTab = ptrElem->GetIndex();
	G4String nameElem = ptrElem->GetName();
	nameElem.toUpper();
	size_t numIsot = ptrElem->GetNumberOfIsotopes();
	G4double A = (ptrElem->GetA())/(g);
	
	if(!numIsot) {
	  if(!elemIndexInMATcard[indElemTab]) {
	    G4double Z = ptrElem->GetZ();
	    G4double density = ptrFracVect[el]*densityMat;
	    elemIndexInMATcard[indElemTab] = indexCount;
	    
	    //write on file MATERIAL card of element
	    PrintMaterial(fout, "MATERIAL  ",
			  Z, A,
			  density,
			  G4double(indexCount),
			  -1,
			  nameElem);
	  }
	}
	
	else {
	  if(numIsot>=2) 
	    isotPresence = 1;
				//loop on isotopes
	  for(unsigned int nis=0; nis<numIsot; nis++) {
	    // G4Isotope pointer 
	    const G4Isotope* ptrIsot = ptrElem->GetIsotope(nis);
	    size_t indIsotTab = ptrIsot->GetIndex();
	    //initialize variables
	    if(!initIsot) {
	      totNumIsot = ptrIsot->GetNumberOfIsotopes();
	      ptrIsotTab = ptrIsot->GetIsotopeTable(); 
	      isotIndexInMATcard = new G4int[totNumIsot];
	      for(unsigned int t=0; t<totNumIsot; t++) 
		isotIndexInMATcard[t] = 0;
	      initIsot = 1;
	    }
	    if(!isotIndexInMATcard[indIsotTab]) {
	      //compute physical data and counters
	      G4int ZIs = ptrIsot->GetZ();
	      G4double AIs = (ptrIsot->GetA())/(g);
	      G4int NIs = ptrIsot->GetN();
	      G4String nameIsot = ptrIsot->GetName();
	      nameIsot.toUpper();
	      G4double* ptrRelAbVect = ptrElem->GetRelativeAbundanceVector();
	      G4double density = 
		ptrFracVect[el]*densityMat*ptrRelAbVect[nis]*AIs/A;
	      G4int index = indexCount;
	      isotIndexInMATcard[indIsotTab] = indexCount;
	      indexCount+=1;
	      
	      //write on file MATERIAL card of isotope
	      PrintMaterial(fout, "MATERIAL  ",
			    G4double(ZIs), AIs,
			    density,
			    G4double(index),
			    NIs,
			    nameIsot);
	    }
	  }
	}
	
      }
      
      if(numElem>1 || isotPresence==1) { 
	// write MATERIAL+COMPOUND card specifing the compound
	
	//flags for writing COMPOUND card
	G4int treCount=0;
	
	//make MATERIAL card for compound, start COMPOUND card
	fout <<"*"<<G4endl;
	fout <<"*   Define GEANT4 compound " << nameMat << G4endl;
	PrintMaterial(fout, "MATERIAL  ",
		      -1, -1,
		      densityMat,
		      G4double(indexMatFluka[j]),
		      -1,
		      nameMat);
  	fout << setw10 << "COMPOUND  ";
	
	
	//write elements in COMPOUND card
	for(unsigned int h=0; h<numElem; h++) {
	  const G4Element * ptrElemMat = ptrMat->GetElement(h);
	  size_t indexElemMat = ptrElemMat->GetIndex();
	  size_t numIsotElem = ptrElemMat->GetNumberOfIsotopes();
	  if(!numIsotElem) {
	    PrintCompound(fout, "COMPOUND  ",
			  treCount,
			  nameMat,
			  -ptrFracVect[h],
			  G4double(elemIndexInMATcard[indexElemMat]));

	    if(treCount==3)
	      treCount=0;
	    treCount+=1;
	  }
	  else {
	    G4double * ptrIsotAbbVect = 
	      ptrElemMat->GetRelativeAbundanceVector();
	    
	    for(unsigned int iso=0; iso<numIsotElem; iso++) {
	      const G4Isotope * ptrIsotElem =ptrElemMat->GetIsotope(iso);
	      size_t indexIsotMat = ptrIsotElem->GetIndex();
	      G4double isotAbundPerVol = 
		ptrIsotAbbVect[iso]*Avogadro*densityMat*
		ptrFracVect[h]/(ptrElemMat->GetA()/(g));
	      
	      PrintCompound(fout, "COMPOUND  ",
			    treCount,
			    nameMat,
			    isotAbundPerVol,
			    G4double(isotIndexInMATcard[indexIsotMat]));
	      if(treCount==3)
		treCount=0;
	      treCount+=1;
	    }
	  }
	}
	
	//end COMPOUND card
	if(treCount==1) 
	  fout << setw10 << " " << setw10 << " "
	       << setw10 << " " << setw10 << " "
	       << nameMat << G4endl;
	if(treCount==2) 
	  fout << setw10 << " " << setw10 << " "
	       << nameMat << G4endl;
	if(treCount==3) 
	  fout << nameMat << G4endl;
	fout << "*" << G4endl;
      }
      
      
    } // end for loop
    delete elemIndexInMATcard;
    if(initIsot) 
      delete isotIndexInMATcard;

  } // end if (ptrMat)
  
  // *** material-volume correspondence
  PrintHeader(fout,"GEANT4 MATERIAL ASSIGNMENTS");
  
  //initializations
  G4int indexMatOld = 0;
  G4int indexRegFlukaFrom = 0;
  G4int indexRegFlukaTo = 0;
  G4int existTo = 0;
  G4int flagField = 0;
  G4int lastFlagField = 0;
  
  //open file for volume-index correspondence
  G4std::ofstream vout("Volumes_index.inp");
  
  //... and write title
  vout << "*" << G4endl;
  vout << "********************  GEANT4 VOLUMES *******************\n";
  vout << "*" << G4endl;
  
  //loop over all volumes...
  for(unsigned int l=0;l<numVol;l++) {
    //index of the region
    G4VPhysicalVolume * ptrVol = (*pVolStore)[l];
    G4LogicalVolume * ptrLogVol = ptrVol->GetLogicalVolume();
    G4int indexRegFluka = l+1;
    
    
    //write index volume and name on file Volumes_index.inp
    vout.setf(G4std::ios::left,G4std::ios::adjustfield);
    vout << setw10 << indexRegFluka;
    vout << G4std::setw(20) << ptrVol->GetName() << G4std::setw(20) << "";
    if(ptrVol->IsReplicated()) {
      EAxis axis;
      G4int nRep;
      G4double width;
      G4double offset;
      G4bool consum;
      ptrVol->GetReplicationData(axis,nRep,width,offset,consum);
      vout.setf(G4std::ios::left,G4std::ios::adjustfield);
      vout << setw10 << "Repetion Nb: " << G4std::setw(3) << nRep;
    }
    vout << G4endl;
    
    //check if Magnetic Field is present in the region
    G4FieldManager * pMagFieldMan = ptrLogVol->GetFieldManager();
    const G4Field * pMagField = 0;
    if(pMagFieldMan) 
      pMagField = pMagFieldMan->GetDetectorField();
    if(pMagField)	  
      flagField = 1;
    else 	  
      flagField = 0;


    //index of material in the region
    G4Material * ptrMat = ptrLogVol->GetMaterial();
    G4int indexMat; 
    if(ptrMat) {
      size_t indexMatGeant = ptrMat->GetIndex();
      indexMat = indexMatFluka[indexMatGeant];
    }
    else 
      indexMat = 2;
	
    //if materials are repeated
    if(indexMat==indexMatOld && flagField==lastFlagField) {
      indexRegFlukaTo=indexRegFluka;
      existTo=1; 
      if((l+1)==numVol) {	
	fout << setw10 << G4double(indexRegFlukaTo);
	fout << setw10 << "0.0";
	fout << setw10 << G4double(flagField);
	fout << setw10 << "0.0" << G4endl;
      }
      
    }
    else {
      //write on file ASSIGNMAT card 
      
      //first complete last material card
      if(!existTo) {
	if(l) {
	  fout << setw10 << "0.0";
	  fout << setw10 << "0.0";
	  fout << setw10 << G4double(lastFlagField);
	  fout << setw10 << "0.0" << G4endl;
	}
      }
      else {
	fout << setw10 <<  G4double(indexRegFlukaTo);
	fout << setw10 << "0.0";
	fout << setw10 << G4double(lastFlagField);
	fout << setw10 << "0.0" << G4endl;
      }
      
      // begin material card		
      fout << setw10 << "ASSIGNMAT ";
      fout.setf(static_cast<G4std::ios::fmtflags>(0),G4std::ios::floatfield);	
      fout << setw10 << setfixed
	   << G4std::setprecision(1) << G4double(indexMat);
      fout << setw10 << G4double(indexRegFluka);
      
      existTo=0;
      indexRegFlukaFrom=indexRegFluka;
      
      if((l+1)==numVol) {	
	fout << setw10 << "0.0";
	fout << setw10 << "0.0";
	fout << setw10 << G4double(flagField);
	fout << setw10 << "0.0" << G4endl;
      }
    }
    lastFlagField = flagField;
    indexMatOld = indexMat; 
  } // end of loop 
  
  //assign material 1 to black-hole=n.vol+1
  fout << setw10 << "ASSIGNMAT ";
  fout << setw10 << "1.0";
  fout << setw10 << G4double(numVol+1);
  fout << setw10 << "0.0";
  fout << setw10 << "0.0";
  fout << setw10 << "0.0";
  fout << setw10 << "0.0" << G4endl;
  
  // *** magnetic field ***
  if(fTransportationManager->GetFieldManager()->DoesFieldExist()) {
    PrintHeader(fout,"GEANT4 MAGNETIC FIELD");
    
    //get magnetic field pointer
    const G4Field * pMagField = 
      fTransportationManager->GetFieldManager()->GetDetectorField();     
    
    //if uniform magnetic field, get value
    G4double Bx=0.0;
    G4double By=0.0;
    G4double Bz=0.0;
    
    if(pMagField) {
      //G4ThreeVector* Field = new G4ThreeVector(1.,2.,3.);
      const G4UniformMagField *pUnifMagField = 
	dynamic_cast<const G4UniformMagField*>(pMagField);
      if(pUnifMagField) {
	G4double *pFieldValue = 0;
	G4double *point = new G4double[3];
	point[0] = 0.;
	point[1] = 0.;
	point[2] = 0.;
	pUnifMagField->GetFieldValue(point,pFieldValue);
	//non capisco perche' l'instruzione seguente non fa linkare. Indaga!!
	//per ora posso usare GetFieldValue: va bene lo stesso. 
	//G4ThreeVector FieldValue = pUnifMagField->GetConstantFieldValue();
	Bx = pFieldValue[0];
	By = pFieldValue[1];
	Bz = pFieldValue[2];
      }
      
    }
    
    //write MGNFIELD card 
    fout << setw10 << "MGNFIELD  ";
    fout << setw10 << "";
    fout << setw10 << "";
    fout << setw10 << "";
    fout.setf(static_cast<G4std::ios::fmtflags>(0),G4std::ios::floatfield);
    fout << setw10 << setfixed
	 << G4std::setprecision(4) << Bx
	 << setw10 << By
	 << setw10 << Bz 
	 << G4endl;
  } // end if magnetic field
  
  vout.close();
  fout.close();
  delete [] indexMatFluka;
#ifdef G4GEOMETRY_DEBUG
  G4cout << "<== Flugg FGeometryInit::createFlukaMatFile()" << G4endl;
#endif
}
