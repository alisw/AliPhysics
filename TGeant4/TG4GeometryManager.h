// $Id$
// Category: geometry
//
// Geant4 implementation of the MonteCarlo interface methods                    
// for building Geant4 geometry and access to it

#ifndef TG4_GEOMETRY_MANAGER_H
#define TG4_GEOMETRY_MANAGER_H

#include "TG4NameMap.h"
#include "TG4Globals.h"
#include "TG3Cut.h"
#include "TG3Flag.h"

#include <globals.hh>
#include <G3SensVolVector.hh>

#include <Rtypes.h>

#include "g4std/fstream"
#include "g4std/vector"

class TG4CutVector;
class TG4FlagVector;
class TG4GeometryOutputManager;

class G4Material;
class G4VPhysicalVolume;

class TG4GeometryManager
{
  public:
    TG4GeometryManager();
    // --> protected
    // TG4GeometryManager(const TG4GeometryManager& right);
    virtual ~TG4GeometryManager();

    // static access method
    static TG4GeometryManager* Instance();

    //
    // methods (from the base class)
    
    // detector composition
    void  Material(Int_t& kmat, const char* name, Float_t a, 
                     Float_t z, Float_t dens, Float_t radl, Float_t absl,
                     Float_t* buf, Int_t nwbuf);
    void  Mixture(Int_t& kmat, const char *name, Float_t *a, 
                     Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat);
    void  Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
                     Float_t stemax, Float_t deemax, Float_t epsil, 
		     Float_t stmin, Float_t* ubuf, Int_t nbuf);
    void  Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
                     Float_t thetaY, Float_t phiY, Float_t thetaZ, 
		     Float_t phiZ);

    // NEW - for G4 only
    G4Material* MixMaterials(G4String name, G4double density,
                     TG4StringVector* matNames, TG4doubleVector* matWeights);

    // functions from GBASE 
    void  Ggclos(); 

    // functions from GCONS 
    void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		         Float_t &dens, Float_t &radl, Float_t &absl,
		         Float_t* ubuf, Int_t& nbuf); 
    void  Gstpar(Int_t itmed, const char *param, Float_t parval); 
    void  Gsckov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			 Float_t *absco, Float_t *effic, Float_t *rindex); 

    // functions from GGEOM 
    Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np); 
    void  Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                        Int_t iaxis);
 
    void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Float_t c0i, Int_t numed); 
    void  Gsdvt(const char *name, const char *mother, Float_t step,
                         Int_t iaxis, Int_t numed, Int_t ndvmx); 
    void  Gsdvt2(const char *name, const char *mother, Float_t step, 
                         Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx); 

    void  Gsord(const char *name, Int_t iax); 
    void  Gspos(const char *name, Int_t nr, const char *mother,  
                        Float_t x, Float_t y, Float_t z, Int_t irot, 
                        const char *konly); 
    void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np); 
        
    // Euclid		       
    void WriteEuclid(const char*, const char*, Int_t, Int_t); //new
		               
    // get methods
    Int_t VolId(const Text_t* volName) const;                
    const char* VolName(Int_t id) const; //new
    Int_t NofVolumes() const; 
    
    // end of methods
    // 

    //
    // methods for Geant4 only 
 
    G4VPhysicalVolume* CreateG4Geometry();
    void ReadG3Geometry(G4String filePath);
    void UseG3TrackingMediaLimits();
    void ClearG3Tables();       
    void ClearG3TablesFinal();
    void OpenOutFile(G4String filePath);
    void PrintNameMap();
    
    // set methods
    void SetWriteGeometry(G4bool writeGeometry);
    void SetMapSecond(const G4String& name);

    // get methods
          // volumes
    Int_t NofG3Volumes() const; 
    Int_t NofG4LogicalVolumes() const; 
    Int_t NofG4PhysicalVolumes() const; 
    Int_t NofSensitiveDetectors() const; 
    G4bool IsG3Volume(G4String lvName) const;
    void G4ToG3VolumeName(G4String& name) const;
    const G4String& GetMapSecond(const G4String& name);

          // sensitive volumes
    G3SensVolVector GetG3SensVolVector() const;

          // materials
    G4int GetMediumId(G4Material* material) const;    
    G4double GetEffA(G4Material* material) const;
    G4double GetEffZ(G4Material* material) const;

    // end of methods for Geant4 only 
    //
    
  protected:
    TG4GeometryManager(const TG4GeometryManager& right);

    // operators
    TG4GeometryManager& operator=(const TG4GeometryManager& right);

  private:
    // methods
    G4double* CreateG4doubleArray(Float_t* array, G4int size) const;
    G4String  CutName(const char* name) const;
    void GstparCut(G4int itmed, TG3Cut par, G4double parval);
    void GstparFlag(G4int itmed, TG3Flag par, G4double parval);
    void FillMediumIdVector();
        
    // static data members
    static TG4GeometryManager*  fgInstance;     //this instance

    // data members
    TG4GeometryOutputManager*   fOutputManager; //output manager 
    TG4NameMap       fNameMap;         //map of volumes names to modules names
    vector<G4int>    fMediumIdVector;  //vector of second indexes for materials
    G4int            fMediumCounter;   //global medium counter
    G4int            fMaterialCounter; //global material counter
    G4int            fMatrixCounter;   //global matrix counter
    G4bool           fUseG3TMLimits;   //if true: G3 limits are passed to G4 
                                       //(in development)
    G4bool           fWriteGeometry;   //if true: geometry parameters are written
                                       //in a file (ASCII)  
};

// inline methods
inline TG4GeometryManager* TG4GeometryManager::Instance()
{ return fgInstance; }

inline G3SensVolVector TG4GeometryManager::GetG3SensVolVector() const
{ return G3SensVol; }

#endif //TG4_GEOMETRY_MANAGER_H

