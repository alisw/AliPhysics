#ifndef FLUKAMATERIAL_HH
#define FLUKAMATERIAL_HH 1

#include "G4String.hh"

class FlukaLowMat;

#include <map>

class FlukaMaterial;
typedef std::map<G4String, FlukaMaterial*, std::less<G4String> > FlukaMaterialsTable;
typedef std::map<G4String, FlukaMaterial*, std::less<G4String> >::const_iterator FlukaMaterialsIterator;
typedef std::map<G4int, FlukaMaterial*, std::less<G4int> > FlukaMaterialsIndexTable;
typedef std::map<G4int, FlukaMaterial*, std::less<G4int> >::const_iterator FlukaMaterialsIndexIterator;

class FlukaMaterial {
public:

  //Constructor
  FlukaMaterial(const G4String& name, 
		G4int Z, G4double A, 
		G4double density,
		G4int N=0);
  virtual ~FlukaMaterial();

  //Getters
  G4int    GetZ() const {return fZ;}
  G4double GetA() const {return fA;}
  G4double GetDensity() const {return fDensity;}
  G4int    GetN() const {return fN;}
  G4String GetName() const {return fName;}
  G4String GetRealName() const;
  G4int    GetIndex() const {return fIndex;}
  const FlukaLowMat* GetLowMat() const {return fFlukaLowMat;}


  //Setters
  void SetZ(G4int Z) {fZ = Z;}
  void SetA(G4double A) {fA = A;}
  void SetDensity(G4double d) {fDensity = d;}
  void SetN(G4int N) {fN = N;}
  void SetName(const G4String& n) {fName = n;}

  //Static
  static const FlukaMaterialsTable* GetMaterialTable() { 
    return &fFlukaMaterials;}
  static FlukaMaterial* GetFlukaMaterial(const G4String& name) { 
    return fFlukaMaterials[name];}
  static std::ostream& PrintMaterialsByName(std::ostream& os);
  static std::ostream& PrintMaterialsByIndex(std::ostream& os);

protected:
  //Other
  void AddLowMat(const G4String& name);

public:
  G4String fName;
  G4int    fZ;
  G4double fA;
  G4double fDensity;
  G4int    fN; //Isotope number
  FlukaLowMat* fFlukaLowMat;

  size_t   fIndex;

  static   FlukaMaterialsTable fFlukaMaterials; //Map of name/material
  static   FlukaMaterialsIndexTable fFlukaIndexMaterials; //Map index/material

};

std::ostream& operator<<(std::ostream& os, const FlukaMaterial& material);


#endif
