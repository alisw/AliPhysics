// $Id$
// Flugg tag $Name$

#include "FGeometryInit.hh"
#include <stdio.h>

FGeometryInit * FGeometryInit::flagInstance=0;

FGeometryInit* FGeometryInit::GetInstance() 
{
  if (!flagInstance) new FGeometryInit();
  
  return flagInstance;
}  


FGeometryInit::FGeometryInit()
  {
   flagInstance = this;
   fTransportationManager = G4TransportationManager::GetTransportationManager();
  }


FGeometryInit::~FGeometryInit()
  {
  DeleteHistories();
  ptrGeoMan->OpenGeometry();  
  delete fTransportationManager;
  delete ptrJrLtGeant;
  DelHistArray();
  
  //keep ATTENTION: never delete a pointer twice!
  }


G4Navigator* FGeometryInit::getNavigatorForTracking()
{
   return fTransportationManager->GetNavigatorForTracking();
} 

void FGeometryInit::setDetConstruction(G4VUserDetectorConstruction* detector)
{
  fDetector = detector;;
}

void FGeometryInit::setDetector()
{
  myTopNode = fDetector->Construct(); 
}

void FGeometryInit::setMotherVolume()
{
    fTransportationManager->GetNavigatorForTracking()->SetWorldVolume(myTopNode);
}

void FGeometryInit::closeGeometry()
{
   ptrGeoMan = G4GeometryManager::GetInstance();
   ptrGeoMan->OpenGeometry();

   //true argoment allows voxel construction; if false voxels are built 
   //only for replicated volumes  
   ptrGeoMan->CloseGeometry(true);
}
 
G4FieldManager * FGeometryInit::getFieldManager()
{
  return fTransportationManager->GetFieldManager();
}

//*************************************************************************

void FGeometryInit::InitHistArray()
{
  ptrArray = new G4int[1000000];
  for(G4int i=0;i<1000000;i++) ptrArray[i]=0;
}

void FGeometryInit::DelHistArray()
{
  delete  ptrArray;
}

G4int * FGeometryInit::GetHistArray()
{
  return ptrArray;
}



//*************************************************************************
//jrLtGeant stores all crossed lattice volume histories.

void FGeometryInit::InitJrLtGeantArray()
{
#ifdef G4GEOMETRY_DEBUG
  G4cout<<"Initializing JrLtGeant array"<<G4endl;
#endif
  ptrJrLtGeant = new G4int[10000];
  for(G4int x=0;x<10000;x++) ptrJrLtGeant[x]=-1;
  flagLttcGeant = -1;
}


G4int * FGeometryInit::GetJrLtGeantArray()
{
  return ptrJrLtGeant;
}


G4int FGeometryInit::GetLttcFlagGeant()
{
  return flagLttcGeant;
}

void FGeometryInit::SetLttcFlagGeant(G4int newFlagLttc)
{
  // Added by A.Solodkov
  if (newFlagLttc >= 10000) {
      G4cout<<"Problems in FGeometryInit::SetLttcFlagGeant"<<G4endl;
      G4cout<<"Index newFlagLttc="<<newFlagLttc<<" is outside array bounds"<<G4endl;
      G4cout<<"Better to stop immediately !"<<G4endl;
      exit(1);
  }
  flagLttcGeant = newFlagLttc;
}
 
void FGeometryInit::PrintJrLtGeant()
{
#ifdef G4GEOMETRY_DEBUG
   //G4cout<<"jrLtGeant:"<<G4endl;
   //for(G4int y=0;y<=flagLttcGeant;y++)
   //
   //	 G4cout<<"jrLtGeant["<<y<<"]="<<ptrJrLtGeant[y]<<G4endl;
#endif
}
	 
//**************************************************************************

void FGeometryInit::PrintHistories()
{
/*
#ifdef G4GEOMETRY_DEBUG
  G4cout<<"Touch hist:"<<G4endl;
  G4cout<<*(ptrTouchHist->GetHistory())<<G4endl;
  G4cout<<"Tmp hist:"<<G4endl;
  G4cout<<*(ptrTempNavHist->GetHistory())<<G4endl;
  G4cout<<"Old hist:"<<G4endl;
  G4cout<<*(ptrOldNavHist->GetHistory())<<G4endl;
#endif
*/
}




void FGeometryInit::InitHistories()
{  
#ifdef G4GEOMETRY_DEBUG
	G4cout <<" InitHistories start" << G4endl;
#endif
  //init utility histories with navigator history

    G4cout << fTransportationManager<< G4endl;
    G4cout << fTransportationManager->GetNavigatorForTracking()<< G4endl;

    ptrTouchHist = fTransportationManager->GetNavigatorForTracking()->CreateTouchableHistory();
    G4cout << "Touchable history " << G4endl;
    ptrTempNavHist = fTransportationManager->GetNavigatorForTracking()->CreateTouchableHistory();   
    ptrOldNavHist = new G4TouchableHistory();
#ifdef G4GEOMETRY_DEBUG
	G4cout <<" InitHistories end" << G4endl;
#endif
}

void FGeometryInit::DeleteHistories()
  {
    delete ptrTouchHist;
    delete ptrOldNavHist;
    delete ptrTempNavHist;

#ifdef G4GEOMETRY_DEBUG
  G4cout<<"Deleting step-history objects at end of run!"<<G4endl;
#endif
  }

G4TouchableHistory * FGeometryInit::GetTouchableHistory()
{
  return ptrTouchHist;
}

G4TouchableHistory * FGeometryInit::GetOldNavHist()
{
  return ptrOldNavHist;
}

G4TouchableHistory * FGeometryInit::GetTempNavHist()
{
  return ptrTempNavHist;
}


void FGeometryInit::UpdateHistories(const G4NavigationHistory * history,
				    const G4int & flagHist)
{
  PrintHistories();

#ifdef G4GEOMETRY_DEBUG
  G4cout<<"...updating histories!"<<G4endl;
#endif

  G4VPhysicalVolume * pPhysVol = history->GetTopVolume();
  
  switch (flagHist)
    {
    case 0:
      {
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
      }

    case 1:
      {
	//this is the case when a new history is given to the 
	//navigator but old history has to be kept (e.g. LOOKZ
	//is call during an event);
	//touchable history has not been updated jet, so:

	ptrTouchHist->UpdateYourself(pPhysVol,history);
	ptrTempNavHist->UpdateYourself(pPhysVol,history);
	PrintHistories();
	break;
      }

    case 2:
      {
	//this is the case when the touchable history has been 
	//updated by a LocateGlobalPointAndUpdateTouchable call

	G4VPhysicalVolume * pPhysVolTemp = ptrTempNavHist->GetVolume();
	ptrOldNavHist->UpdateYourself(pPhysVolTemp,
				      ptrTempNavHist->GetHistory());

	ptrTempNavHist->UpdateYourself(pPhysVol,history);
	PrintHistories();
	break;
      }
    default:
      {
	G4cout<<"ERROR in updating step-histories!"<<G4endl;
	break;
      }
    }

}

//*****************************************************************************


void FGeometryInit::createFlukaMatFile()
{
// ultima modifica Sara Vanini 1/III/99
// NOMI DI ELEMENTI E COMPOSTI: i nomi devono essere scritti maiuscolo,
// secondo lo standard di fluka. Devono inoltre essere uguali ai nomi dei
// materiali fluka - vedere il manuale di fluka - perche` il programma
// carichi le giuste sezioni d`urto, e uguali ai nomi inclusi nel .pemf.
// Altrimenti l`utente deve definirsi le CARDS LOW-MAT, e costruirsi il
// .pemf, per avere le giuste sezioni d`urto caricate in memoria.

     //flag
#ifdef G4GEOMETRY_DEBUG
     G4cout<<"================== FILEWR ================="<<G4endl;
#endif 

      //open file for output
     ofstream fout("flukaMat.inp");  

     //PhysicalVolumeStore, Volume and MaterialTable pointers
     G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
     G4int numVol = G4int(pVolStore->size());
   
     G4int* indexMatFluka = 0; 
     static G4Material * ptrMat = 0;
     G4int x = 0;
     while(!ptrMat && x<numVol)
       {
         G4VPhysicalVolume * ptrVol = (*pVolStore)[x];
         G4LogicalVolume * ptrLogVol = ptrVol->GetLogicalVolume();
         ptrMat = ptrLogVol->GetMaterial();
         x+=1;
       }

     if(ptrMat)
       {
	static const G4MaterialTable * ptrMatTab = G4Material::GetMaterialTable();

	//number of materials, elements, variable initialisations
	static size_t totNumMat = G4Material::GetNumberOfMaterials();
#ifdef G4GEOMETRY_DEBUG
	G4cout<<"Number of materials: "<<totNumMat<<G4endl;
#endif 
	const G4Element * ptrElem = ptrMat->GetElement(0);
	static size_t totNumElem = ptrElem->GetNumberOfElements();
	static const G4ElementTable * ptrElemTab = ptrElem->GetElementTable();
	G4int * elemIndexInMATcard = new G4int[totNumElem];
	for(G4int t=0; t<totNumElem; t++) 
	  elemIndexInMATcard[t] = 0;
	static const G4IsotopeTable * ptrIsotTab;
	G4int initIsot = 0;
	static size_t totNumIsot;
        G4int* isotIndexInMATcard = 0;
       	G4int isotPresence = 0; 


	// title
	fout<<"*\n"<<"*\n"<<"*\n";
        fout<<"*********************  GEANT4 ELEMENTS AND COMPOUNDS *********************\n"<<"*\n";
	fout<<("*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...")<<G4endl;
	fout<<("*")<<G4endl;

	// *** loop over G4Materials to assign Fluka index 
	G4int indexCount=3;
	indexMatFluka = new G4int[totNumMat];
	for(G4int i=0; i<totNumMat; i++)
		{
		//pointer material, state 
		ptrMat = (*ptrMatTab)[i];
		G4double denMat = ptrMat->GetDensity();
		G4String nameMat = ptrMat->GetName();
		// Fluka index: bh=1, vacuum=2, others=3,4..
		if(denMat<=1.00e-10*g/cm3)
		//N.B. fluka density limit decided on XI-`98  
			{
			indexMatFluka[i]=2;	
			}
		else
			{  
			indexMatFluka[i]=indexCount;
			indexCount+=1;
			}

		// *** write single-element material MATERIAL card
		size_t numElem = ptrMat->GetNumberOfElements();
		if(numElem==1)
		  {
		    G4int index = indexMatFluka[i];
		    const G4Element * ptrElem = ptrMat->GetElement(0);
		    size_t indElemTab = ptrElem->GetIndex();
		    size_t numIsot = ptrElem->GetNumberOfIsotopes();
		    G4double A = (ptrElem->GetA())/(g);
		    if(!numIsot)
		      {
			if(index!=2 && !elemIndexInMATcard[indElemTab])
			  {
			    G4double Z = ptrElem->GetZ();
			    elemIndexInMATcard[indElemTab] = index;
			    G4String nameEl = ptrElem->GetName();
			    nameEl.toUpper();

			    //write on file MATERIAL card of element
			    fout<<G4std::setw(10)<<"MATERIAL  ";
			    fout.setf(0,G4std::ios::floatfield);
			    fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)
			      		<<G4std::setprecision(1)<<Z;
			    fout<<G4std::setw(10)<<G4std::setprecision(3)<<A;
			    fout.setf(0,G4std::ios::floatfield);
			    fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::scientific)
			       		<<G4std::setprecision(3)<<denMat/(g/cm3);
			    fout.setf(0,G4std::ios::floatfield);
			    fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)<<
			       		G4std::setprecision(1)<<G4double(index);
			    fout<<G4std::setw(10)<<" ";
			    fout<<G4std::setw(10)<<" ";
			    fout<<nameEl<<G4endl;
			  }
		       }
		     if(numIsot==1)
		       {
		          // G4Isotope pointer 
			  const G4Isotope* ptrIsot = ptrElem->GetIsotope(0);
			  size_t indIsotTab = ptrIsot->GetIndex();
			  //initialize variables
		          if(!initIsot)
			    {
	        	      totNumIsot = ptrIsot->GetNumberOfIsotopes();
			      ptrIsotTab = ptrIsot->GetIsotopeTable(); 
			      isotIndexInMATcard = new G4int[totNumIsot];
			      for(G4int t=0; t<totNumIsot; t++) isotIndexInMATcard[t] = 0;
			      initIsot = 1;
			    }
		       	  if(!isotIndexInMATcard[indIsotTab])
			    {//compute physical data and counters
		              G4int ZIs = ptrIsot->GetZ();
			      G4double AIs = (ptrIsot->GetA())/(g);
			      G4int NIs = ptrIsot->GetN();
			      G4String nameIsot = ptrIsot->GetName();
			      nameIsot.toUpper();
			      G4double* ptrRelAbVect = ptrElem->GetRelativeAbundanceVector();
			      isotIndexInMATcard[indIsotTab] = index;
			
			      //write on file MATERIAL card of isotope
			      fout<<G4std::setw(10)<<"MATERIAL  ";
			      fout.setf(0,G4std::ios::floatfield);
			      fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)
				<<G4std::setprecision(1)<<G4double(ZIs);
			      fout<<G4std::setw(10)<<G4std::setprecision(3)<<AIs;
			      fout.setf(0,ios::floatfield);
			      fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::scientific)
				<<G4std::setprecision(3)<<denMat/(g/cm3);
			      fout.setf(0,G4std::ios::floatfield);
			      fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)
				<<G4std::setprecision(1)<<G4double(index);
			      fout<<G4std::setw(10)<<" ";
      			      fout<<G4std::setw(10)<<G4double(NIs);
			      fout<<nameIsot<<G4endl;
			    }
		        }// end if(numIsot==1)
		    }// end if(numElem==1)
		}// end for loop
		 

	// *** material definitions: elements, compound made of G4Elements
	// or made of G4Materials
	for(G4int j=0; j<totNumMat; j++)
		{
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
		for(G4int el=0; el<numElem; el++)
		      	{
			//compute physical data, initialize variables
			const G4Element * ptrElem = ptrMat->GetElement(el);
			size_t indElemTab = ptrElem->GetIndex();
			G4String nameElem = ptrElem->GetName();
			nameElem.toUpper();
			size_t numIsot = ptrElem->GetNumberOfIsotopes();
			G4double A = (ptrElem->GetA())/(g);

			if(!numIsot)
				{
				if(!elemIndexInMATcard[indElemTab])
					{
					G4double Z = ptrElem->GetZ();
					G4double density = ptrFracVect[el]*densityMat;
					elemIndexInMATcard[indElemTab] = indexCount;

					//write on file MATERIAL card of element
					fout<<G4std::setw(10)<<"MATERIAL  ";
					fout.setf(0,G4std::ios::floatfield);
					fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)
						<<G4std::setprecision(1)<<Z;
					fout<<G4std::setw(10)<<G4std::setprecision(3)<<A;
					fout.setf(0,G4std::ios::floatfield);
					fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::scientific)
						<<G4std::setprecision(3)<<density;
					fout.setf(0,G4std::ios::floatfield);
					fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)<<
						G4std::setprecision(1)<<G4double(indexCount);
					fout<<G4std::setw(10)<<" ";
					fout<<G4std::setw(10)<<" ";
					fout<<nameElem<<G4endl;
					indexCount+=1;
					}
				}

			else
				{
				if(numIsot>=2) isotPresence = 1;
				//loop on isotopes
				for(G4int nis=0; nis<numIsot; nis++)
					{
					// G4Isotope pointer 
					const G4Isotope* ptrIsot = ptrElem->GetIsotope(nis);
					size_t indIsotTab = ptrIsot->GetIndex();
					//initialize variables
					if(!initIsot)
						{
						totNumIsot = ptrIsot->GetNumberOfIsotopes();
						ptrIsotTab = ptrIsot->GetIsotopeTable(); 
						isotIndexInMATcard = new G4int[totNumIsot];
						for(G4int t=0; t<totNumIsot; t++) isotIndexInMATcard[t] = 0;
						initIsot = 1;
						}
					if(!isotIndexInMATcard[indIsotTab])
						{//compute physical data and counters
						G4int ZIs = ptrIsot->GetZ();
						G4double AIs = (ptrIsot->GetA())/(g);
						G4int NIs = ptrIsot->GetN();
						G4String nameIsot = ptrIsot->GetName();
						nameIsot.toUpper();
						G4double* ptrRelAbVect = ptrElem->GetRelativeAbundanceVector();
						G4double density = ptrFracVect[el]*densityMat*
							ptrRelAbVect[nis]*AIs/A;
						G4int index = indexCount;
						isotIndexInMATcard[indIsotTab] = indexCount;
						indexCount+=1;
			
			  			//write on file MATERIAL card of isotope
						fout<<G4std::setw(10)<<"MATERIAL  ";
						fout.setf(0,G4std::ios::floatfield);
						fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)
							<<G4std::setprecision(1)<<G4double(ZIs);
						fout<<G4std::setw(10)<<G4std::setprecision(3)<<AIs;
						fout.setf(0,G4std::ios::floatfield);
						fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::scientific)
							<<G4std::setprecision(3)<<density;
						fout.setf(0,G4std::ios::floatfield);
						fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)
							<<G4std::setprecision(1)<<G4double(index);
						fout<<G4std::setw(10)<<" ";
						fout<<G4std::setw(10)<<G4double(NIs);
						fout<<nameIsot<<G4endl;
						}
					}
				}
	
			}

		if(numElem>1 || isotPresence==1)
			{ 
			// write MATERIAL+COMPOUND card specifing the compound

			//flags for writing COMPOUND card
			G4int treCount=0;

			//make MATERIAL card for compound, start COMPOUND card
			fout<<"*"<<G4endl;
			fout<<"*   Define GEANT4 compound "<<nameMat<<G4endl;
			fout<<G4std::setw(10)<<"MATERIAL  ";
			fout.setf(0,G4std::ios::floatfield);
			fout<<G4std::setw(10)<<" "<<G4std::setw(10)<<" ";
			fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::scientific)
				<<G4std::setprecision(3)<<densityMat;
			fout.setf(0,G4std::ios::floatfield);
			fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)<<G4std::setprecision(1)
				<<G4double(indexMatFluka[j]);
			fout<<G4std::setw(10)<<" ";
			fout<<G4std::setw(10)<<" ";
			fout<<nameMat<<G4endl;
			fout<<G4std::setw(10)<<"COMPOUND  ";


 			//write elements in COMPOUND card
			for(G4int h=0; h<numElem; h++)
				{
				const G4Element * ptrElemMat = ptrMat->GetElement(h);
				size_t indexElemMat = ptrElemMat->GetIndex();
				size_t numIsotElem = ptrElemMat->GetNumberOfIsotopes();
				if(!numIsotElem)
					{		
					if(treCount==3)
					 	{
 						treCount=0;
	 					fout<<nameMat<<G4endl;
						fout<<G4std::setw(10)<<"COMPOUND  ";
						}

					fout.setf(0,G4std::ios::floatfield);
					fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::scientific)
						<<G4std::setprecision(2)<<-ptrFracVect[h];
					fout.setf(0,G4std::ios::floatfield);
					fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)<<
						G4std::setprecision(1)<<G4double(elemIndexInMATcard[indexElemMat]);
					treCount+=1;
					}
				else
					{
					G4double * ptrIsotAbbVect = ptrElemMat->GetRelativeAbundanceVector();
  
					for(G4int iso=0; iso<numIsotElem; iso++)
						{
					        const G4Isotope * ptrIsotElem =ptrElemMat->GetIsotope(iso);
						size_t indexIsotMat = ptrIsotElem->GetIndex();
						G4double isotAbundPerVol = 
						   ptrIsotAbbVect[iso]*Avogadro*densityMat*
						   ptrFracVect[h]/(ptrElemMat->GetA()/(g));
						
						if(treCount==3)
							{
							treCount=0;
							fout<<nameMat<<G4endl;
							fout<<G4std::setw(10)<<"COMPOUND  ";
							}
						fout.setf(0,G4std::ios::floatfield);
						fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::scientific)
							<<G4std::setprecision(2)<<isotAbundPerVol;
						fout.setf(0,G4std::ios::floatfield);
						fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)<<
						G4std::setprecision(1)<<G4double(isotIndexInMATcard[indexIsotMat]);
						treCount+=1;
						}
					}
				}

			//end COMPOUND card
			if(treCount==1) fout<<G4std::setw(10)<<" "<<G4std::setw(10)<<" "<<
			  G4std::setw(10)<<" "<<G4std::setw(10)<<" "<<nameMat<<G4endl;
			if(treCount==2) fout<<G4std::setw(10)<<" "<<G4std::setw(10)<<
			  " "<<nameMat<<G4endl;
			if(treCount==3) fout<<nameMat<<G4endl;
			fout<<"*"<<G4endl;
			}


		} // end for loop
	delete elemIndexInMATcard;
       	if(initIsot) delete isotIndexInMATcard;

       } // end if (ptrMat)

     // *** material-volume correspondence
     fout<<"*\n"<<"*\n"<<"*\n";
     fout<<"********************  GEANT4 MATERIAL ASSIGNMENTS *******************\n"<<"*\n";
     fout<<("*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...")<<G4endl;
     fout<<("*")<<G4endl;

     //initializations
     G4int indexMatOld = 0;
     G4int indexRegFlukaFrom = 0;
     G4int indexRegFlukaTo = 0;
     G4int existTo = 0;
     G4int flagField = 0;
     G4int lastFlagField = 0;

     //open file for volume-index correspondence
     ofstream vout("Volumes_index.inp");

     //... and write title
     vout<<"*"<<G4endl;
     vout<<"********************  GEANT4 VOLUMES *******************\n";
     vout<<("*")<<G4endl;

     //loop su tutti i volumi
     for(G4int l=0;l<numVol;l++)
       {
       	//index of the region
       	G4VPhysicalVolume * ptrVol = (*pVolStore)[l];
       	G4LogicalVolume * ptrLogVol = ptrVol->GetLogicalVolume();
       	G4int indexRegFluka = l+1;


        //write index volume and name on file Volumes_index.inp
	vout.setf(G4std::ios::left,G4std::ios::adjustfield);
        vout<<G4std::setw(10)<<indexRegFluka;
	vout<<G4std::setw(20)<<ptrVol->GetName()<<G4std::setw(20)<<"";
        if(ptrVol->IsReplicated())
                {
                EAxis axis;
                G4int nRep;
                G4double width;
                G4double offset;
                G4bool consum;
                ptrVol->GetReplicationData(axis,nRep,width,offset,consum);
		vout.setf(G4std::ios::left,G4std::ios::adjustfield);
		vout<<G4std::setw(10)<<"Repetion Nb: "<<G4std::setw(3)<<nRep;
                }
        vout<<G4endl;

       	//check if Magnetic Field is present in the region
       	G4FieldManager * pMagFieldMan = ptrLogVol->GetFieldManager();
	const G4Field * pMagField = 0;
        if(pMagFieldMan) pMagField = pMagFieldMan->GetDetectorField();
	if(pMagField)	  flagField = 1;
       	else 	  flagField = 0;


       	//index of material in the region
       	G4Material * ptrMat = ptrLogVol->GetMaterial();
       	G4int indexMat; 
       	if(ptrMat)
       		{
       		size_t indexMatGeant = ptrMat->GetIndex();
       		indexMat = indexMatFluka[indexMatGeant];
       		}
       	else indexMat = 2;
	
       	//if materials are repeated
       	if(indexMat==indexMatOld && flagField==lastFlagField)
       		{
       		indexRegFlukaTo=indexRegFluka;
       		existTo=1; 
       		if(l==(numVol-1))
       			{	
       			fout<<G4std::setw(10)<<G4double(indexRegFlukaTo);
       			fout<<G4std::setw(10)<<"0.0";
       			fout<<G4std::setw(10)<<G4double(flagField);
       			fout<<G4std::setw(10)<<"0.0"<<G4endl;
       			}

	       	}
	       else
	       	{
	       	//write on file ASSIGNMAT card 

	       	//first complete last material card
	       	if(!existTo) 
	       		{
	       		if(l) 
	       			{
	       			fout<<G4std::setw(10)<<"0.0";
	       			fout<<G4std::setw(10)<<"0.0";
	       			fout<<G4std::setw(10)<<G4double(lastFlagField);
	       			fout<<G4std::setw(10)<<"0.0"<<G4endl;
	       			}
	       		}
	       	else 	
	       		{
	       		fout<<G4std::setw(10)<<G4double(indexRegFlukaTo);
	       		fout<<G4std::setw(10)<<"0.0";
	       		fout<<G4std::setw(10)<<G4double(lastFlagField);
	       		fout<<G4std::setw(10)<<"0.0"<<G4endl;
	       		}

	       	// begin material card		
	       	fout<<G4std::setw(10)<<"ASSIGNMAT ";
	       	fout.setf(0,G4std::ios::floatfield);	
	       	fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)<<
	  			G4std::setprecision(1)<<G4double(indexMat);
	       	fout<<G4std::setw(10)<<G4double(indexRegFluka);
	       	
	       	existTo=0;
	       	indexRegFlukaFrom=indexRegFluka;

	       	if(l==(numVol-1))
	       		{	
	       		fout<<G4std::setw(10)<<"0.0";
	       		fout<<G4std::setw(10)<<"0.0";
	       		fout<<G4std::setw(10)<<G4double(flagField);
	       		fout<<G4std::setw(10)<<"0.0"<<G4endl;
	       		}
	       	}
	       lastFlagField = flagField;
	       indexMatOld = indexMat; 
              } // end of loop ??

	//assign material 1 to black-hole=n.vol+1
	fout<<G4std::setw(10)<<"ASSIGNMAT ";
	fout<<G4std::setw(10)<<"1.0";
	fout<<G4std::setw(10)<<G4double(numVol+1);
	fout<<G4std::setw(10)<<"0.0";
	fout<<G4std::setw(10)<<"0.0";
	fout<<G4std::setw(10)<<"0.0";
	fout<<G4std::setw(10)<<"0.0"<<G4endl;

        // *** magnetic field ***
	if(fTransportationManager->GetFieldManager()->DoesFieldExist())
	  {
	    fout<<"*\n"<<"*\n"<<"*\n";
	    fout<<"***********************  GEANT4 MAGNETIC FIELD ************************\n"<<"*\n";
	    fout<<("*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...")<<G4endl;
	    fout<<("*")<<G4endl;

	    //get magnetic field pointer
	    const G4Field * pMagField = fTransportationManager->GetFieldManager()->GetDetectorField();     

	    //if uniform magnetic field, get value
	    G4double Bx=0.0;
	    G4double By=0.0;
	    G4double Bz=0.0;

	    if(pMagField) 
	      {
		G4ThreeVector*Field = new G4ThreeVector(1.,2.,3.);
		const G4UniformMagField *pUnifMagField = 
		  dynamic_cast<const G4UniformMagField*>(pMagField);
		if(pUnifMagField)
		  {
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
            fout<<G4std::setw(10)<<"MGNFIELD  ";
            fout<<G4std::setw(10)<<"";
            fout<<G4std::setw(10)<<"";
            fout<<G4std::setw(10)<<"";
	    fout.setf(0,G4std::ios::floatfield);
	    fout<<G4std::setw(10)<<G4std::setiosflags(G4std::ios::fixed)<<G4std::setprecision(4)<<Bx;
            fout<<G4std::setw(10)<<By;
            fout<<G4std::setw(10)<<Bz<<G4endl;
          } // end if magnetic field

        vout.close();
	fout.close();
     	delete [] indexMatFluka;
}
