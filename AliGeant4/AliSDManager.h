// $Id$
// Category: geometry
//
// Class providing methods for creating sensitive detectors
// and switching between lego sensitive detectors and 
// the standard ones.

#ifndef ALI_SD_MANAGER_H
#define ALI_SD_MANAGER_H

#include <globals.hh>

//class AliSDMessenger;
class AliLego;
class AliModule;

class G4VPhysicalVolume;
class G4LogicalVolume;

class AliSDManager 
{  
  public:
    // --> protected
    // AliSDManager();
    // AliSDManager(const AliSDManager& right);
    virtual ~AliSDManager();
    
    // static methods
    static AliSDManager* Instance();
    
    // methods
    void CreateSD(G4LogicalVolume* lv, AliModule* module) const;
    AliModule* FindAliModule(G4LogicalVolume* lv) const;
    void SetLego(AliLego* lego) const;
    void UnsetLego() const;

    // set/get methods 
    inline void SetNofLVWithSD(G4int nofLV);
    inline G4int GetNofLVWithSD() const;
    
  protected:
    AliSDManager();
    AliSDManager(const AliSDManager& right);

    // operators
    AliSDManager& operator=(const AliSDManager& right);
  
  private:
    // methods
    void CreateLegoSD(G4LogicalVolume* lv, AliLego* lego) const;
    void UnsetLegoSD(G4LogicalVolume* lv) const;

    // static data members
    static AliSDManager*    fgInstance;   //this instance

    // data members
    G4int                   fNofLVWithSD; //counter of logical volumes
                                          //with sensitive detectors
    //AliSDManagerMessenger*  fMessenger;
};             

// inline methods

inline G4int AliSDManager::GetNofLVWithSD() const
{ return fNofLVWithSD; }

inline void AliSDManager::SetNofLVWithSD(G4int nofLV)
{ fNofLVWithSD = nofLV; }

#endif //ALI_SD_MANAGER_H
