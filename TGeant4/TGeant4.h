// $Id$
// Category: run
//
// Geant4 implementation of the MonteCarlo interface                      

#ifndef TGEANT4_H
#define TGEANT4_H

#include "AliMC.h"
#include "AliMCProcess.h"

class TG4VRunConfiguration;
class TG4GeometryManager;
class TG4PhysicsManager;
class TG4StepManager;
class TG4VisManager;
class TG4RunManager;
class TG4Messenger;

class AliDecayer;

class TGeant4: public AliMC
{
  public:
    TGeant4(const char* name, const char* title,
            TG4VRunConfiguration* configuration, int argc, char** argv);
    TGeant4(const char* name, const char* title,
            TG4VRunConfiguration* configuration);
    // --> protected
    // TGeant4();
    // TGeant4(const TGeant4& right);
    virtual ~TGeant4();

    //
    // methods for building/management of geometry
    // ------------------------------------------------
    //

    // functions from GBASE 
    virtual void  FinishGeometry(); 
                  //Ggclos(); 

    // functions from GCONS 
    virtual void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		         Float_t &dens, Float_t &radl, Float_t &absl,
		         Float_t* ubuf, Int_t& nbuf); 

    // detector composition
    virtual void  Material(Int_t& kmat, const char* name, Float_t a, 
                     Float_t z, Float_t dens, Float_t radl, Float_t absl,
                     Float_t* buf, Int_t nwbuf);
    virtual void  Mixture(Int_t& kmat, const char *name, Float_t *a, 
                     Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat);
    virtual void  Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
                     Float_t stemax, Float_t deemax, Float_t epsil, 
		     Float_t stmin, Float_t* ubuf, Int_t nbuf);
    virtual void  Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
                     Float_t thetaY, Float_t phiY, Float_t thetaZ, 
		     Float_t phiZ);
    virtual void  Gstpar(Int_t itmed, const char *param, Float_t parval); 
    virtual void  Gsckov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			 Float_t *absco, Float_t *effic, Float_t *rindex); 

    // functions from GGEOM 
    virtual Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np); 
    virtual void  Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                        Int_t iaxis); 
    virtual void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Float_t c0i, Int_t numed); 
    virtual void  Gsdvt(const char *name, const char *mother, Float_t step, 
                        Int_t iaxis, Int_t numed, Int_t ndvmx); 
    virtual void  Gsdvt2(const char *name, const char *mother, Float_t step, 
                         Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx); 
    virtual void  Gsord(const char *name, Int_t iax); 
    virtual void  Gspos(const char *name, Int_t nr, const char *mother,  
                        Float_t x, Float_t y, Float_t z, Int_t irot, 
                        const char *konly); 
    virtual void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np); 
    
    // Euclid		       
    virtual void WriteEuclid(const char* fileName, const char* topVol, 
                             Int_t number, Int_t nlevel);
		               
    // get methods
    virtual Int_t VolId(const Text_t* volName) const;
    virtual const char* VolName(Int_t id) const;
    virtual Int_t NofVolumes() const;

    //
    // methods for physics management
    // ------------------------------------------------
    //
 
    virtual void BuildPhysics();

    // set methods
    virtual void SetCut(const char* cutName, Float_t cutValue);
    virtual void SetProcess(const char* flagName, Int_t flagValue);
    virtual Float_t Xsec(char* reac, Float_t energy, Int_t part, Int_t mate);
    virtual void SetExternalDecayer(AliDecayer* decayer); //NEW

    // get methods
    virtual AliDecayer* Decayer() const; //NEW
 
        // particle table usage         
    virtual Int_t IdFromPDG(Int_t pdgID) const;
    virtual Int_t PDGFromId(Int_t mcID) const;
    virtual void  DefineParticles();       

    //
    // methods for step management
    // ------------------------------------------------
    //

    // action methods
    virtual void StopTrack();
    virtual void StopEvent();   

    // set methods
    virtual void SetMaxStep(Float_t);
    virtual void SetMaxNStep(Int_t);
    virtual void SetUserDecay(Int_t);  //NEW

    // get methods
         // tracking volume(s) 
    virtual Int_t CurrentVolID(Int_t& copyNo) const;
    virtual Int_t CurrentVolOffID(Int_t off, Int_t& copyNo) const;
    virtual const char* CurrentVolName() const;
    virtual const char* CurrentVolOffName(Int_t off) const;
    virtual Int_t CurrentMaterial(Float_t &a, Float_t &z, 
                    Float_t &dens, Float_t &radl, Float_t &absl) const;  
    virtual void  Gmtod(Float_t* xm, Float_t* xd, Int_t iflag);
    virtual void  Gdtom(Float_t* xd, Float_t* xm, Int_t iflag);
    virtual Float_t MaxStep() const;
    virtual Int_t GetMaxNStep() const;
    virtual Int_t GetMedium() const;

        // tracking particle 
        // dynamic properties
    virtual void    TrackPosition(TLorentzVector& position) const;
    virtual void    TrackMomentum(TLorentzVector& momentum) const;
    virtual void    TrackVertexPosition(TLorentzVector& position) const;
    virtual void    TrackVertexMomentum(TLorentzVector& momentum) const;
    virtual Float_t TrackStep() const;
    virtual Float_t TrackLength() const; 
    virtual Float_t TrackTime() const;
    virtual Float_t Edep() const;
        // static properties
    virtual Int_t   TrackPid() const;
    virtual Float_t TrackCharge() const;
    virtual Float_t TrackMass() const;
    virtual Float_t Etot() const;

        // track status
    virtual Bool_t  IsTrackInside() const;
    virtual Bool_t  IsTrackEntering() const;
    virtual Bool_t  IsTrackExiting() const;
    virtual Bool_t  IsTrackOut() const;
    virtual Bool_t  IsTrackDisappeared() const;
    virtual Bool_t  IsTrackStop() const;
    virtual Bool_t  IsTrackAlive() const;
    virtual Bool_t  IsNewTrack() const;

        // secondaries
    virtual Int_t NSecondaries() const;
    virtual void  GetSecondary(Int_t isec, Int_t& particleId, 
                    TLorentzVector& position, TLorentzVector& momentum);
    virtual AliMCProcess ProdProcess() const; 

	// random number generator	    
    virtual void Rndm(Float_t* array, const Int_t size) const;    
  
    //
    // methods for visualization
    // ------------------------------------------------
    //
    // functions for drawing
    virtual void  DrawOneSpec(const char* name);
    virtual void  Gsatt(const char* name, const char* att, Int_t val);
    virtual void  Gdraw(const char* name, Float_t theta, Float_t phi,
		        Float_t psi, Float_t u0, Float_t v0,
		        Float_t ul, Float_t vl);

    //
    // NEW
    // Geant3 specific methods
    // !!! need to be transformed to common interface
    //
    virtual void Gdopt(const char* name , const char* value);
    virtual void SetClipBox(const char *name, Float_t xmin, Float_t xmax,
		       Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax);
    virtual void DefaultRange();
    virtual void Gdhead(Int_t isel, const char* name, Float_t chrsiz);   
    virtual void Gdman(Float_t u, Float_t v, const char* type);
    virtual void SetColors();
    virtual void Gtreve();
    virtual void GtreveRoot();
    virtual void Gckmat(Int_t itmed, char* natmed);
    virtual void InitLego();
    virtual void Gfpart(Int_t ipart, char *name, Int_t& itrtyp,  
		       Float_t& amass, Float_t& charge, Float_t& tlife);
    virtual void Gspart(Int_t ipart, const char *name, Int_t itrtyp,  
		       Float_t amass, Float_t charge, Float_t tlife); 

    //
    // methods for run control
    // ------------------------------------------------
    //

    virtual void Init();
    virtual void ProcessEvent();
    virtual void ProcessRun(Int_t nofEvents);

        // UI control methods
    void StartGeantUI();	
    void StartRootUI();	
    void ProcessGeantMacro(const char* macroName);
    void ProcessGeantCommand(const char* commandPath);

        // get methods
    virtual Int_t CurrentEvent() const; 

  protected:
    TGeant4();
    TGeant4(const TGeant4& right);

    // operators
    TGeant4& operator=(const TGeant4& right);

  private:
    // data members
    TG4GeometryManager*  fGeometryManager; //geometry manager
    TG4PhysicsManager*   fPhysicsManager;  //physics manager
    TG4StepManager*      fStepManager;     //step manager
    TG4VisManager*       fVisManager;      //visualization manager
    TG4RunManager*       fRunManager;      //run manager
    TG4Messenger*        fMessenger;       //messenger
};

#ifndef __CINT__

// inline methods
#include "TGeant4.icc"

#endif
#endif // TGEANT4_H

