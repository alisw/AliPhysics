#include "FlukaCompound.hh"
#include "G4ios.hh"
#include "WrapUtils.hh"

FlukaCompoundsTable FlukaCompound::fFlukaCompounds;

FlukaCompound::FlukaCompound(const G4String& name, 
			     G4double density, 
			     G4int nelem):
  fNMaterials(nelem),
  fNAdded(0) {
  //Initialise arrays for element indices and fractions
  fElIndex = new G4int[nelem];
  fFraction = new G4double[nelem];
  for (G4int i = 0; i < nelem; i++) {
    fElIndex[i] = 0;
    fFraction[i] = 0;
  }

  //Build the associated material
  fFlukaMaterial = new FlukaMaterial(name, 0, 0, density);

  //If it already exists it will have got a new name
  G4String testname(fFlukaMaterial->GetRealName());
  if (testname != name) {
    G4cerr << "INFO: Found two materials with the same name! ("
	   << name << ")" << G4endl;
    G4cerr << "      Renaming to " << testname << G4endl;
  }
  fFlukaCompounds[testname] = this;
}

FlukaCompound::~FlukaCompound() {
  delete[] fElIndex;
  delete[] fFraction;
}


void FlukaCompound::AddElement(G4int index, G4double fraction) {
  if (fNAdded < fNMaterials) {
    fElIndex[fNAdded] = index;
    fFraction[fNAdded] = fraction;
    fNAdded++;
  }
  else {
    G4cerr << "ERROR: Trying to add to many elements to compound \'" 
	   << GetName() << "\'!" << G4endl;
    G4cerr << "       Last element not added!!!" << G4endl;
  }
}

std::ostream& FlukaCompound::PrintCompounds(std::ostream& os) {
  PrintHeader(os, "COMPOUNDS");  
  
  for (FlukaCompoundsIterator i = fFlukaCompounds.begin(); 
       i != fFlukaCompounds.end(); 
       i++) {
    FlukaCompound* flucomp     = (*i).second;
    os << *flucomp;
  }

  return os;
}

std::ostream& operator<<(std::ostream& os, const FlukaCompound& flucomp) {
  G4int nmats = flucomp.GetNMaterials();
  G4String matName = flucomp.GetName();
  G4String matRealName = flucomp.GetRealName().substr(0,8);
  //Some comment
  os << "* " << matName << " COMPOUND (" << nmats << ")" << G4endl;
  
  //Material card
  //os << *(flucomp.GetFlukaMaterial());

  //The card
  G4int counttothree = 0;
  os << setw10 <<"COMPOUND  ";
  for (G4int i = 0; i < nmats; i++) {
    os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
    os << setw10
       << setfixed
       << std::setprecision(6)
       << flucomp.GetMaterialFraction(i);
    os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
    os << setw10
       << setfixed
       << std::setprecision(1)
       << G4double(flucomp.GetMaterialIndex(i));
    counttothree++;
    if (counttothree == 3 ) {
      os << matRealName;
      os << G4endl;
      if ( (i+1) != nmats)
	os << setw10 <<"COMPOUND  ";
      counttothree = 0;
    }
  }

  if ((nmats % 3) != 0) {
    //Unless we have 3, 6, 9... submaterials we need to put some empty
    //space and the compound name
    for (G4int i = 0; i < (3 - (nmats % 3)); i++)
      os << setw10 << " " << setw10 << " ";
    os << matRealName;
    os << G4endl;
  }
  
  return os;
}
