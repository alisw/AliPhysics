// @(#):$Name$:$Id$
// Author: Andrei Gheata 10/07/2003

#ifndef ROOT_TFlukaMCGeometry
#define ROOT_TFlukaMCGeometry 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//
// Class TFlukaMCGeometry
// --------------------
// Implementation of the TVirtualMCGeometry interface
// for defining and using TGeo geometry with FLUKA.
//  

#ifndef ROOT_Riostream
#include "Riostream.h"
#endif

#ifndef ROOT_TVirtualMCGeometry
#include "TVirtualMCGeometry.h"
#endif

class TFlukaMCGeometry : public TVirtualMCGeometry {

  public:
    TFlukaMCGeometry();
    TFlukaMCGeometry(const char* name, const char* title);
    virtual ~TFlukaMCGeometry();
  
  // functions from GCONS 
    virtual void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
		       Float_t &dens, Float_t &radl, Float_t &absl,
		       Float_t* ubuf, Int_t& nbuf);
    virtual void  Gfmate(Int_t imat, char *name, Double_t &a, Double_t &z,  
		       Double_t &dens, Double_t &radl, Double_t &absl,
		       Double_t* ubuf, Int_t& nbuf);
    // detector composition
    virtual void  Material(Int_t& kmat, const char* name, Double_t a, 
                     Double_t z, Double_t dens, Double_t radl, Double_t absl,
                     Float_t* buf, Int_t nwbuf);
    virtual void  Material(Int_t& kmat, const char* name, Double_t a, 
                     Double_t z, Double_t dens, Double_t radl, Double_t absl,
                     Double_t* buf, Int_t nwbuf);
    virtual void  Mixture(Int_t& kmat, const char *name, Float_t *a, 
                     Float_t *z, Double_t dens, Int_t nlmat, Float_t *wmat);
    virtual void  Mixture(Int_t& kmat, const char *name, Double_t *a, 
                     Double_t *z, Double_t dens, Int_t nlmat, Double_t *wmat);
    virtual void  Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
                     Double_t stemax, Double_t deemax, Double_t epsil, 
		     Double_t stmin, Float_t* ubuf, Int_t nbuf);
    virtual void  Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
                     Double_t stemax, Double_t deemax, Double_t epsil, 
		     Double_t stmin, Double_t* ubuf, Int_t nbuf);
    virtual void  Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, 
                     Double_t thetaY, Double_t phiY, Double_t thetaZ, 
		     Double_t phiZ);

    // functions from GGEOM 
    void     Gmtod(Float_t* xm, Float_t* xd, Int_t iflag);
  
    void     Gmtod(Double_t* xm, Double_t* xd, Int_t iflag);
  
    void     Gdtom(Float_t* xd, Float_t* xm, Int_t iflag);
  
    void     Gdtom(Double_t* xd, Double_t* xm, Int_t iflag);
    virtual Int_t  Gsvolu(const char *name, const char *shape, Int_t nmed,  
                          Float_t *upar, Int_t np); 
    virtual Int_t  Gsvolu(const char *name, const char *shape, Int_t nmed,  
                          Double_t *upar, Int_t np); 
    virtual void  Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis); 
    virtual void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Double_t c0i, Int_t numed); 
    virtual void  Gsdvt(const char *name, const char *mother, Double_t step, 
                         Int_t iaxis, Int_t numed, Int_t ndvmx); 
    virtual void  Gsdvt2(const char *name, const char *mother, Double_t step, 
                         Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx); 
    virtual void  Gsord(const char *name, Int_t iax); 
    virtual void  Gspos(const char *name, Int_t nr, const char *mother,  
                         Double_t x, Double_t y, Double_t z, Int_t irot, 
                         const char *konly); 
    virtual void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Double_t x, Double_t y, Double_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np);
    virtual void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Double_t x, Double_t y, Double_t z, Int_t irot,
                         const char *konly, Double_t *upar, Int_t np);
    virtual void  Gsbool(const char* /*onlyVolName*/, const char* /*manyVolName*/) {}
  
    
    // functions for drawing
    //virtual void  DrawOneSpec(const char* name);
    void  Gsatt(const char* name, const char* att, Int_t val);
    //virtual void  Gdraw(const char*,Double_t theta, Double_t phi,
    //		        Double_t psi, Double_t u0, Double_t v0,
    //		        Double_t ul, Double_t vl);

    // Euclid
    //virtual void  WriteEuclid(const char*, const char*, Int_t, Int_t);
		               
    // get methods
    Int_t         CurrentVolID(Int_t& copyNo) const;
    Int_t         CurrentVolOffID(Int_t off, Int_t& copyNo) const;
    const char*   CurrentVolName() const;
    const char*   CurrentVolOffName(Int_t off) const;
    Int_t         GetMedium() const;
    Int_t        *GetRegionList(Int_t imed, Int_t &nreg);
    Int_t        *GetMaterialList(Int_t imat, Int_t &nreg);
    Int_t         GetFlukaMaterial(Int_t imed) const;
    Int_t         GetLastMaterialIndex() const {return fLastMaterial;}
    virtual Int_t VolId(const Text_t* volName) const;
    virtual const char* VolName(Int_t id) const;
    virtual Int_t NofVolumes() const;
    virtual Int_t VolId2Mate(Int_t id) const;

   // FLUKA specific methods
    void          CreateFlukaMatFile(const char *fname=0);
    void          PrintHeader(ofstream &out, const char *text) const;
    void          SetMreg(Int_t mreg);
    void          SetNextRegion(Int_t mreg, Int_t latt);
    Int_t         RegionId() const; 
    void          ToFlukaString(TString &str) const;

  private:
    Double_t* CreateDoubleArray(Float_t* array, Int_t size) const;
    void     Vname(const char *name, char *vname) const;
   
    TFlukaMCGeometry(const TFlukaMCGeometry& rhs);
    TFlukaMCGeometry& operator=(const TFlukaMCGeometry& rhs) {return (*this);}

    static TFlukaMCGeometry*  fgInstance; // singleton instance
    Int_t        fLastMaterial;           // last FLUKA material index
    Int_t        fNextRegion;             // next region number
    Int_t        fNextLattice;            // next lattice history
    Int_t       *fRegionList;             //! region list matching a given medium number
  ClassDef(TFlukaMCGeometry,1)  //Virtual MonteCarlo Interface
};

#endif //ROOT_TFlukaMCGeometry
