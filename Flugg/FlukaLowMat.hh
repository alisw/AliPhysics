#ifndef FLUKALOWMAT_HH
#define FLUKALOWMAT_HH 1

#include "G4String.hh"
#include "FlukaMaterial.hh"

class FlukaLowMat {
public:

  //Constructor
  FlukaLowMat(const G4String& name, FlukaMaterial*);

  //Getters
  G4String GetName() const {return fName;}
  inline G4int GetIndex() const;
  const FlukaMaterial* GetFlukaMaterial() const {return fFlukaMaterial;}

  //Setters
  void SetName(const G4String& n) {fName = n;}

public:
  G4String fName;
  FlukaMaterial* fFlukaMaterial;

};


inline G4int FlukaLowMat::GetIndex() const {return fFlukaMaterial->GetIndex();}

std::ostream& operator<<(std::ostream& os, const FlukaLowMat& flowmat);

#endif

