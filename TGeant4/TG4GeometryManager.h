// $Id$
// Category: geometry
//
// Author: V. Berejnoi, I. Hrivnacova
//
// Class TG4GeometryManager
// ------------------------
// Geant4 implementation of the MonteCarlo interface methods                    
// for building Geant4 geometry and access to it.

#ifndef TG4_GEOMETRY_MANAGER_H
#define TG4_GEOMETRY_MANAGER_H

#include "TG4NameMap.h"
#include "TG4IntMap.h"
#include "TG4Globals.h"

#include <globals.hh>

#include <Rtypes.h>


class TG4GeometryOutputManager;
class TG4GeometryServices;
class TG4G3CutVector;
class TG4G3ControlVector;

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
    void  Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, 
                     Double_t thetaY, Double_t phiY, Double_t thetaZ, 
		     Double_t phiZ);
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
                         Double_t *upar, Int_t np); 
    Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np); 
    void  Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                        Int_t iaxis);
    void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Double_t c0i, Int_t numed); 
    void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Float_t c0i, Int_t numed); 
    void  Gsdvt(const char *name, const char *mother, Double_t step,
                         Int_t iaxis, Int_t numed, Int_t ndvmx); 
    void  Gsdvt(const char *name, const char *mother, Float_t step,
                         Int_t iaxis, Int_t numed, Int_t ndvmx); 
    void  Gsdvt2(const char *name, const char *mother, Double_t step, 
                         Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx); 
    void  Gsdvt2(const char *name, const char *mother, Float_t step, 
                         Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx); 

    void  Gsord(const char *name, Int_t iax); 
    void  Gspos(const char *name, Int_t nr, const char *mother,  
                        Double_t x, Double_t y, Double_t z, Int_t irot, 
                        const char *konly); 
    void  Gspos(const char *name, Int_t nr, const char *mother,  
                        Float_t x, Float_t y, Float_t z, Int_t irot, 
                        const char *konly); 
    void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Double_t x, Double_t y, Double_t z, Int_t irot,
                         const char *konly, Double_t *upar, Int_t np); 
    void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np); 
        
    // Euclid		       
    void WriteEuclid(const char* fileName, const char* topVolName, 
                         Int_t number, Int_t nlevel); //new
		               
    // end of methods
    // 

    //
    // methods for Geant4 only 
 
    G4VPhysicalVolume* CreateG4Geometry();
    void SetUserLimits(const TG4G3CutVector& cuts,
                       const TG4G3ControlVector& controls) const;
    void ReadG3Geometry(G4String filePath);
    void UseG3TrackingMediaLimits();
    void ClearG3Tables();       
    void ClearG3TablesFinal();
    void OpenOutFile(G4String filePath);
    void CloseOutFile();
    
    // set methods
    void SetWriteGeometry(G4bool writeGeometry);
    void SetMapSecond(const G4String& name);
     
  protected:
    TG4GeometryManager(const TG4GeometryManager& right);

    // operators
    TG4GeometryManager& operator=(const TG4GeometryManager& right);

  private:
    // methods
    void FillMediumMap();
        
    // static data members
    static TG4GeometryManager*  fgInstance;     //this instance

    // data members
    TG4GeometryOutputManager*   fOutputManager;   //output manager 
    TG4GeometryServices*        fGeometryServices;//geometry services
    TG4IntMap        fMediumMap;       //map of volumes names to medias IDs
    TG4NameMap       fNameMap;         //map of volumes names to modules names
    TG4StringVector  fMaterialNameVector; // vector of material names sorted in the
                                          // the order of materials in G3Mat
    TG4StringVector  fMediumNameVector;   // vector of material names sorted in the
                                          // the order of medias in G3Med
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

#endif //TG4_GEOMETRY_MANAGER_H

