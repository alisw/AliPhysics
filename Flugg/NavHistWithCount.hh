
// Flugg tag 

// 
// NavHistWithCount.hh - Sara Vanini 
// last modified 2/II/'99
// This class stores a G4NavigationHistory and it adds to it a 
// counter for secondary particles. 
//
//


#ifndef NAVHISTWITHCOUNT_HH
#define NAVHISTWITHCOUNT_HH

#include "G4NavigationHistory.hh"


class NavHistWithCount
{
public:
   NavHistWithCount(const G4NavigationHistory &history);
   ~NavHistWithCount();
   void UpdateCount(G4int incrCount);
   G4NavigationHistory* GetNavHistPtr();
   G4int GetCount();
   G4int GetDelateFlag();
   G4int GetCheckInd();
   void SaveCheckInd(G4int);

private:
   G4int count;
   G4NavigationHistory * fhistory;
   G4int fDelate;
   G4int fCheck;
};

#include "NavHistWithCount.icc"

#endif



