
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapIncrHist.hh - Sara Vanini
//
// Wrapper for updating  secondary particles history counter. 
// If counter=0 the history is deleted. 
//
// modified 14/I/9999
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
//////////////////////////////////////////////////////////////////


#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "NavHistWithCount.hh"
#include "globals.hh"


void conhwr(G4int& intHist, G4int* incrCount)
{
//flag
#ifdef G4GEOMETRY_DEBUG
  G4cout << "============= CONHWR ==============" << G4endl;    
  G4cout << "Ptr History = " << intHist << G4endl;
#endif 
  
  //get NavHistWithCount pointer
  if(intHist!=-1) {
    NavHistWithCount* ptrNavHistCount=
      reinterpret_cast<NavHistWithCount*>(intHist); 
 
    //for debugging...
#ifdef G4GEOMETRY_DEBUG
    G4cout << "Secondary counter=" << ptrNavHistCount->GetCount();
    if(*incrCount>0) 
      G4cout << "+" << *incrCount << G4endl;
    if(*incrCount<0) 
      G4cout<< *incrCount << G4endl;   
    if(*incrCount==0) 
      G4cout << G4endl; 
#endif 
    
    //update secondary particles counter
    ptrNavHistCount->UpdateCount(*incrCount);
    
    //delete history if counter=0 or if counter=-1
    G4int counter = ptrNavHistCount->GetCount();
#ifdef G4GEOMETRY_DEBUG
    G4cout << "Counter = " << counter << G4endl;
#endif 
    if(!counter || counter==-1) {
#ifdef G4GEOMETRY_DEBUG
      G4cout << "Deleting Nav Hist object..." << G4endl;
#endif 
      /*
	//for history checking....
	G4int index = ptrNavHistCount->GetCheckInd();
	static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
	G4int * ptrArray = ptrGeoInit->GetHistArray();
	ptrArray[index]=0; 
      */	  
      
      delete ptrNavHistCount;
#ifdef G4GEOMETRY_DEBUG
      G4cout << "end delete" << G4endl;
#endif 
      intHist=-1;
    }
  }
#ifdef G4GEOMETRY_DEBUG
  G4cout << "============= Out of CONHWR ==============" << G4endl;    
#endif 
 }










