// $Id$
// Category: geometry
//
// Geant4 implementation of the MonteCarlo interface methods                    
// for building Geant4 geometry and access to it

#ifndef TG4_GEOMETRY_MANAGER_H
#define TG4_GEOMETRY_MANAGER_H

#include "TG4NameMap.h"
#include "TG4Globals.h"
#include "TG4G3Cut.h"
#include "TG4G3Control.h"

#include <globals.hh>
#include <G3SensVolVector.hh>

#include <Rtypes.h>

#include <g4std/fstream>
#include <g4std/vector>

class TG4CutVector;
class TG4FlagVector;
class TG4GeometryOutputManager;
class TG4GeometryServices;

class G4Material;
class G4VPhysicalVolume;
class G4LogicalVolume;

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

    // functions from GBASE 
    void  Ggclos(); 

    // functions from GCONS 
    void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		         Float_t &dens, Float_t &radl, Float_t &absl,
		         Float_t* ubuf, Int_t& nbuf); 
    void  Gstpar(Int_t itmed, const char *param, Float_t parval); 
    void  SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
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
    void WriteEuclid(const char* fileName, const char* topVolName, 
                         Int_t number, Int_t nlevel); //new
		               
    // get methods
    Int_t VolId(const Text_t* volName) const;                
    const char* VolName(Int_t id) const;
    Int_t NofVolumes() const; 
    Int_t VolId2Mate(Int_t volumeId) const;
    
    // end of methods
    // 

    //
    // methods for Geant4 only 
 
    G4VPhysicalVolume* CreateG4Geometry();
    void ReadG3Geometry(G4String filePath);
    void UseG3TrackingMediaLimits();
    void FillMediumIdVector();
    void ClearG3Tables();       
    void ClearG3TablesFinal();
    void OpenOutFile(G4String filePath);
    void CloseOutFile();
    void PrintNameMap();
    
    // set methods
    void SetWriteGeometry(G4bool writeGeometry);
    void SetMapSecond(const G4String& name);

    // get methods 
    G3SensVolVector GetG3SensVolVector() const;
     
  protected:
    TG4GeometryManager(const TG4GeometryManager& right);

    // operators
    TG4GeometryManager& operator=(const TG4GeometryManager& right);

  private:
    // methods
    void GstparCut(G4int itmed, TG4G3Cut par, G4double parval);
    void GstparControl(G4int itmed, TG4G3Control control, G4double parval);
        
    // static data members
    static TG4GeometryManager*  fgInstance;     //this instance

    // data members
    TG4GeometryOutputManager*   fOutputManager;   //output manager 
    TG4GeometryServices*        fGeometryServices;//geometry services
    TG4intVector  fMediumIdVector;  //vector of second indexes for materials
    TG4NameMap    fNameMap;         //map of volumes names to modules names
    G4int         fMediumCounter;   //global medium counter
    G4int         fMaterialCounter; //global material counter
    G4int         fMatrixCounter;   //global matrix counter
    G4bool        fUseG3TMLimits;   //if true: G3 limits are passed to G4 
                                    //(in development)
    G4bool        fWriteGeometry;   //if true: geometry parameters are written
                                    //in a file (ASCII)  
};

// inline methods
inline TG4GeometryManager* TG4GeometryManager::Instance()
{ return fgInstance; }

inline G3SensVolVector TG4GeometryManager::GetG3SensVolVector() const
{ return G3SensVol; }

#endif //TG4_GEOMETRY_MANAGER_H

