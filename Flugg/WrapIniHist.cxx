
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapIniHist.hh - Sara Vanini
//
// Wrapper for reinitialization of FluggNavigator history.
//
// modified 14/I/99
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
///////////////////////////////////////////////////////////////////


#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "NavHistWithCount.hh"
#include "G4NavigationHistory.hh"
#include "FluggNavigator.hh"
#include "globals.hh"

void inihwr(G4int& intHist)
{
//flag
#ifdef G4GEOMETRY_DEBUG
  G4cout << "============= INIHWR ==============" << G4endl;    
  G4cout << "Ptr History=" <<intHist<< G4endl;
#endif 
  if(intHist==-1) {
    G4cout << "ERROR! This history has been deleted!" << G4endl;
    return;
  }
  else {
    //get NavHistWithCount,G4NavigationHistory,FGeometryInit,
    //FluggNavigator pointers
    NavHistWithCount* ptrNavHistCount=
      reinterpret_cast<NavHistWithCount*>(intHist);
    G4NavigationHistory* ptrNavHist=ptrNavHistCount->GetNavHistPtr();
    static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
    FluggNavigator* ptrNavig = ptrGeoInit->getNavigatorForTracking();
    
    //reinitialize navigator history 
    ptrNavig->UpdateNavigatorHistory(ptrNavHist);
    
    //update utility histories: touch,temp, and reset old history
    ptrGeoInit->UpdateHistories(ptrNavHist,0);
    
    //save new history in jrLtGeant if not present
    G4int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
    G4int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
    
    G4bool intHistInJrLtGeant = false;
    for(G4int h=0; h<=LttcFlagGeant; h++)
      if(jrLtGeant[h]==intHist) 
	intHistInJrLtGeant = true;
    
    if(!intHistInJrLtGeant) {  
      LttcFlagGeant += 1; 
      ptrGeoInit->SetLttcFlagGeant(LttcFlagGeant);
      
      jrLtGeant[LttcFlagGeant]=intHist;
      
#ifdef G4GEOMETRY_DEBUG
      G4cout << "* CONHWR call to increment counter" << G4endl;
#endif 
      G4int incrCount=1;
      conhwr(jrLtGeant[LttcFlagGeant],&incrCount);
    }    
    
    //print history....
#ifdef G4GEOMETRY_DEBUG
    G4cout << "History reinitialized in:" << G4endl;
    G4cout << *ptrNavHist << G4endl;
    ptrGeoInit->PrintJrLtGeant();
#endif
  }
}





