
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapSavHist.hh - Sara Vanini
//
// Wrapper for saving current navigation history (fCheck=default) 
// and returning its pointer. If fCheck=-1 copy of history pointed 
// by intHist is made in NavHistWithCount object, and its pointer 
// is returned. fCheck=1 and fCheck=2 cases are only in debugging 
// version: an array is created by means of FGeometryInit functions
// (but could be a static int * ptrArray = new int[10000] with 
// file scope as well) that stores a flag for deleted/undeleted 
// histories and at the end of event is checked to verify that 
// all saved history objects have been deleted.
//
// modified 6/III/99: history check array implemented
// modified 14/IV/00: fCheck=-1 case modified
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
//////////////////////////////////////////////////////////////////


#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "NavHistWithCount.hh"
#include "G4TouchableHistory.hh"
#include "globals.hh"

G4int isvhwr(const G4int& fCheck, const G4int& intHist)
{
  //flag
#ifdef G4GEOMETRY_DEBUG
  G4cout << "============= ISVHWR ==============" << G4endl;    
  G4cout << "fCheck=" << fCheck << G4endl;
  if(fCheck==-1) 
    G4cout << "intHist=" << intHist  << G4endl;
#endif 
  
  //Geoinit pointer
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  static G4int j=0;
  
  switch (fCheck) {
  case 1:
    {
#ifdef G4GEOMETRY_DEBUG
      G4cout << "Start event." << G4endl;
#endif
      
      return 0;
    }
    
  case 2:
    {
#ifdef G4GEOMETRY_DEBUG
      G4cout << "End event. Check histories in debug-array..." << G4endl;
#endif
      
      //check that all fluka-banked histories have been delated
      // commented by A.Solodkov
      /*
	G4int * ptrArray = ptrGeoInit->GetHistArray();
	
	for(G4int k=0;k<j;k++) {
	NavHistWithCount* ptrNavHistCount=reinterpret_cast
	<NavHistWithCount*>(ptrArray[k]);
	
	if(ptrArray[k] && !ptrNavHistCount->GetDelateFlag()) 
	G4cout << "WARNING! History pointed by " <<ptrArray[k]<<
	" has not been deleted at end of event" << G4endl;
	}
	
	//reinitialise debug histories array
	for(G4int i=0;i<1000000;i++) ptrArray[i]=0;
      */
      j=0;
      
      return 0;
    }
    
  case -1:
    {
      //get history object from intHist
      NavHistWithCount* ptrNavHistCountCopy = 
	reinterpret_cast<NavHistWithCount*>(intHist);
      G4NavigationHistory* ptrNavHistCopy = 
	ptrNavHistCountCopy->GetNavHistPtr();
      
      //copy history in another NavHistWithCount object
      //and save index of check-array
      NavHistWithCount * ptrNavHistCount = 
	new NavHistWithCount(*(ptrNavHistCopy));
      ptrNavHistCount->SaveCheckInd(j);
      G4int intHistCopy = G4int(ptrNavHistCount);
      
      //store ptr in array
      // commented by PS
      //G4int * ptrArray = ptrGeoInit->GetHistArray();
      //ptrArray[j]=intHistCopy;
      j+=1;
      
#ifdef G4GEOMETRY_DEBUG
      G4cout << "Copying history..." << G4endl;
      G4cout << "Ptr History-copy =" <<intHistCopy<< G4endl;
      G4cout<<*ptrNavHistCopy<< G4endl;
#endif 
      
      return intHistCopy;
    }
    
  default:
    {
      //save G4NavigationHistory and index of check-array
      NavHistWithCount *ptrNavHistCount = 
	new NavHistWithCount(*(ptrGeoInit->GetTempNavHist()->GetHistory()));
      
      ptrNavHistCount->SaveCheckInd(j);
      G4int histInt=G4int(ptrNavHistCount);
      
      //store ptr in array
      // comented by PS
      //G4int * ptrArray = ptrGeoInit->GetHistArray();
      // ptrArray[j]=histInt;
      j+=1;
      
#ifdef G4GEOMETRY_DEBUG
      //TouchableHistory 
      G4TouchableHistory * ptrTouchHist = ptrGeoInit->GetTouchableHistory();
      G4cout << "Saving history..." << G4endl;
      G4cout << "Ptr saved History=" << histInt << G4endl;
      G4cout << *(ptrTouchHist->GetHistory()) << G4endl;
#endif 
      
      return histInt;
    }
  }
}






