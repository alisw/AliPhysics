// $Id$
// Category: digits+hits
//
// Geant4 implementation of the MonteCarlo interface methods                    
// for access to Geant4 geometry related with sensitive detectors.

#ifndef TG4_SD_MANAGER_H
#define TG4_SD_MANAGER_H

#include <globals.hh>

#include <Rtypes.h>

class TG4SDServices;
class TG4VSDConstruction;

class TG4SDManager
{
  public:
    TG4SDManager(TG4VSDConstruction* sdConstruction);
    // --> protected
    // TG4SDManager();
    // TG4SDManager(const TG4SDManager& right);
    virtual ~TG4SDManager();

    // static methods
    static TG4SDManager* Instance();

    // methods
    void Initialize();
    
    // AliMC methods
    Int_t VolId(const Text_t* volName) const;                
    const char* VolName(Int_t id) const;
    Int_t NofVolumes() const; 
    Int_t VolId2Mate(Int_t volumeId) const;

    // get methods
    TG4VSDConstruction* GetSDConstruction() const;

  protected:
    TG4SDManager();
    TG4SDManager(const TG4SDManager& right);

    // operators
    TG4SDManager& operator=(const TG4SDManager& right);

    // static data members
    static TG4SDManager* fgInstance;   //this instance
    
    // data members
    TG4VSDConstruction* fSDConstruction; //sensitive detectors construction
    TG4SDServices*      fSDServices;     //services related with sensitive
                                         //detectors

};

// inline methods

inline TG4SDManager* TG4SDManager::Instance() 
{ return fgInstance; }

inline TG4VSDConstruction* TG4SDManager::GetSDConstruction() const
{ return fSDConstruction; }

#endif //TG4_SD_MANAGER_H

