#ifndef FLUKACOMPOUND_HH
#define FLUKACOMPOUND_HH 1

#include "G4String.hh"
#include "FlukaMaterial.hh"

#include <map>

class FlukaCompound;
typedef std::map<G4String, FlukaCompound*, std::less<G4String> > FlukaCompoundsTable;
typedef std::map<G4String, FlukaCompound*, std::less<G4String> >::const_iterator FlukaCompoundsIterator;

class FlukaCompound {
public:

  //Constructor
  FlukaCompound(const G4String& name, G4double density, G4int nElem);
  virtual ~FlukaCompound();

  //Getters
  // * Fluka name
  G4String GetName() const {return fFlukaMaterial->GetName();}
  // * Real name if the fluka name is duplicated
  G4String GetRealName() const {return fFlukaMaterial->GetRealName();}
  // * Density
  G4double GetDensity() const {return fFlukaMaterial->GetDensity();}
  // * Index (it comes from the material)
  G4int    GetIndex() const {return fFlukaMaterial->GetIndex();}
  // * Number of materials, and index and fraction for each material
  G4int    GetNMaterials() const {return fNMaterials;}
  G4int    GetMaterialIndex(G4int i) const {return fElIndex[i];}
  G4double GetMaterialFraction(G4int i) const {return fFraction[i];}

  // * Associated material
  const FlukaMaterial* GetFlukaMaterial() const {return fFlukaMaterial;}
  FlukaMaterial* GetFlukaMaterial() {return fFlukaMaterial;}

  //Setters
  void SetName(const G4String& n) {fFlukaMaterial->SetName(n);}
  void SetDensity(G4double d) {fFlukaMaterial->SetDensity(d);}

  //Other
  void AddElement(G4int index, G4double fraction);

  //Static
  static inline const FlukaCompoundsTable* GetCompoundTable();
  static inline const FlukaCompound* GetFlukaCompound(const G4String& name);
  static std::ostream& PrintCompounds(std::ostream& os);

public:
  G4int     fNMaterials; //Number of elements in total
  G4int     fNAdded;     //Number of elements added
  G4int*    fElIndex;    //Array of element indices
  G4double* fFraction;   //Array of element fracions

  FlukaMaterial* fFlukaMaterial; //Each compound has a "dummy" mat associated

  static   FlukaCompoundsTable fFlukaCompounds;

};

inline const FlukaCompoundsTable* FlukaCompound::GetCompoundTable() {
  return &fFlukaCompounds;
}

inline const FlukaCompound* FlukaCompound::GetFlukaCompound(const G4String& name) { 
  return fFlukaCompounds[name];
}


//Ostream operator
std::ostream& operator<<(std::ostream& os, const FlukaCompound& flucomp);

#endif

