// $Id$
// Category: visualization
//
// G4VisManager derived class that takes care of registering
// graphics syystem and provides 
// Geant4 implementation of the MonteCarlo interface methods                    
// for visualization.

#ifndef TG4_VIS_MANAGER_H
#define TG4_VIS_MANAGER_H
#ifdef G4VIS_USE

#include "TG4G3Attribute.h"

#include <G4VisManager.hh>
#include <g4rw/tpordvec.h>

#include <Rtypes.h>

class TG4VisManager: public G4VisManager 
{
  public:
    TG4VisManager(G4int verboseLevel = 0);  
      // Controls initial verbose level of VisManager and VisMessenger.
      // Can be changed by /vis/set/verbose.
    // --> protected  
    // TG4VisManager(const TG4VisManager& right);
    virtual ~TG4VisManager();  
    
    //----------------------  
    // functions for drawing
    //----------------------
    void DrawOneSpec(const char* name); // not implemented
    
    // see TG4VisManager.cxx for detailed description of Gsatt(), Gdraw()
    void Gsatt(const char* name, const char* att, Int_t val);
    void Gdraw(const char* name, Float_t theta , Float_t phi, 
               Float_t psi, Float_t u0, Float_t v0, Float_t ul, 
	       Float_t vl);
    void SetColors();

  protected:
    TG4VisManager(const TG4VisManager& right);

    // operators
    TG4VisManager& operator=(const TG4VisManager& right);

  private:
    // methods
    //--------
    void RegisterGraphicsSystems();
    G4bool NeedSetColours();
    void SetColourFlag(G4bool value);
    
    // methods used by Gsatt(), Gdraw()
    //---------------------------------
    
    // Get the logical volume list corresponding to NAME
    // 	Either a logical or physical volume name can be supplied
    // Clones of G3VOLUME_NUMBER will be atached to the list	
    G4RWTPtrOrderedVector<G4LogicalVolume> GetLVList(G4String name);

    // Get the physical volume list corresponding to NAME
    G4RWTPtrOrderedVector<G4VPhysicalVolume> GetPVList(G4String name);
    
    // Case insensitive string comparison
    G4bool CaseInsensitiveEqual(const G4String string1,
				const G4String string2);
    
    // Return true if the vis. attributes pointer corresponding to the 
    //	selected volume is shared by others. In this case, duplication
    //  of those is mandatory
    G4bool IsSharedVisAttributes(const G4LogicalVolume* pLV);
    
    // Set an attribute to a specific volume
    void SetG4Attribute(G4LogicalVolume* const lv, const TG4G3Attribute att,
			const G4int val);
    // Set an attribute to the tree coresponding to a volume			
    void SetAtt4Daughters(G4LogicalVolume* const lv, const TG4G3Attribute att,
			const G4int val);			    

    //data members
    //------------
    G4bool fColourFlag; //colour flag
};

// inline methods

inline G4bool TG4VisManager::NeedSetColours()
{ return fColourFlag; }

inline void TG4VisManager::SetColourFlag(G4bool value)
{ fColourFlag = value; }


#endif //G4VIS_USE
#endif //TG4_VIS_MANAGER_H
