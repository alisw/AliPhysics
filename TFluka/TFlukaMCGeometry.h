#ifndef TFLUKAMCGEOMETRY_H
#define TFLUKAMCGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// Class TFlukaMCGeometry
// --------------------
// Implementation of the TVirtualMCGeometry interface
// for defining and using TGeo geometry with FLUKA.
//
// Author: Andrei Gheata 10/07/2003

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TGeoMaterial;


class TFlukaMCGeometry :public TNamed {

  public:
    enum EFlukaLatticeTypes {
       kLttcOutside = 999999999,
       kLttcVirtual = 1000000000
    };   
  
    TFlukaMCGeometry();
    TFlukaMCGeometry(const char* name, const char* title);
    virtual ~TFlukaMCGeometry();
  
    // get methods
    Int_t         GetNstep();  // to be removed
    Int_t         GetMedium() const;
    Int_t        *GetRegionList(Int_t imed, Int_t &nreg);
    Int_t        *GetMaterialList(Int_t imat, Int_t &nreg);
    Int_t         GetFlukaMaterial(Int_t imed) const;
    Int_t         GetLastMaterialIndex() const {return fLastMaterial;}
    virtual Int_t NofVolumes() const;
   // FLUKA specific methods
    void          CreateFlukaMatFile(const char *fname=0);
    void          CreatePemfFile();
    void          PrintHeader(ofstream &out, const char *text) const;
    Bool_t        IsDebugging() const {return fDebug;}
    void          SetDebugMode(Bool_t flag=kTRUE) {fDebug = flag;}
    void          SetMreg(Int_t mreg, Int_t lttc);
    void          SetCurrentRegion(Int_t mreg, Int_t latt);
    void          GetCurrentRegion(Int_t &mreg, Int_t &latt) const {mreg=fCurrentRegion; latt=fCurrentLattice;}
    Int_t         GetCurrentRegion() const {return fCurrentRegion;}
    Int_t         GetDummyRegion() const {return fDummyRegion;}
    Int_t         GetDummyLattice() const {return kLttcVirtual;}
    void          SetNextRegion(Int_t mreg, Int_t latt);
    void          GetNextRegion(Int_t &mreg, Int_t &latt) const {mreg=fNextRegion; latt=fNextLattice;}
    TGeoMaterial *GetMakeWrongMaterial(Double_t z);
    TObjArray    *GetMatList() {return fMatList;}
    TObjArray    *GetMatNames() {return fMatNames;}
    Int_t         GetElementIndex(Int_t z) const;
    Int_t         RegionId() const; 
    void          ToFlukaString(TString &str) const;
    void          FlukaMatName(TString &str) const;
    void          WritePegFile(Int_t imat, Int_t *NoStern, Int_t *ElemError,
                       Int_t *MixError, Int_t *countGas) const;
    Double_t *    GetISSB(Double_t rho, Int_t nElem, Double_t *zelem, Double_t *welem ) const;

    Double_t* CreateDoubleArray(Float_t* array, Int_t size) const;
    void     Vname(const char *name, char *vname) const;
   
  private:
    // Copy constructor and operator= declared but not implemented (-Weff++ flag)
    TFlukaMCGeometry(const TFlukaMCGeometry& rhs);
    TFlukaMCGeometry& operator=(const TFlukaMCGeometry& /*rhs*/); // {return (*this);}

    static TFlukaMCGeometry*  fgInstance; // singleton instance
    Bool_t       fDebug;                  // debug flag
    Int_t        fLastMaterial;           // last FLUKA material index
    Int_t        fDummyRegion;            // index of dummy region
    Int_t        fCurrentRegion;          // current region number
    Int_t        fCurrentLattice;         // current lattice history
    Int_t        fNextRegion;             // next region number
    Int_t        fNextLattice;            // next lattice history
    Int_t       *fRegionList;             //! region list matching a given medium number
    Int_t        fIndmat;                 // material index where pemf file creation starts
    TObjArray   *fMatList;                //! material list as known by FLUKA
    TObjArray   *fMatNames;               //! list of FLUKA material names
    
  ClassDef(TFlukaMCGeometry,1)  //Virtual MonteCarlo Interface
};

#endif //ROOT_TFlukaMCGeometry
