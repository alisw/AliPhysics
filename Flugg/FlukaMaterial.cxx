#include "FlukaMaterial.hh"
#include "FlukaLowMat.hh"
#include "WrapUtils.hh"
#include "G4ios.hh"

FlukaMaterialsTable FlukaMaterial::fFlukaMaterials;
FlukaMaterialsIndexTable FlukaMaterial::fFlukaIndexMaterials;

FlukaMaterial::FlukaMaterial(const G4String& name, 
			     G4int Z, G4double A, 
			     G4double density,
			     G4int N):
  fName(name),
  fZ(Z),
  fA(A),
  fDensity(density),
  fN(N),
  fFlukaLowMat(0) {

  G4String testname(name);
  G4int matrep = 1;
  while (fFlukaMaterials[testname] && matrep < 100) {
    matrep++;
    char smatrep[3];
    sprintf(smatrep,"%.2d",matrep);

    testname = name;
    if (testname.length() <= 6)
      testname += smatrep;
    else
      testname.replace(6,testname.length()-6, smatrep, 2);

#ifdef G4GEOMETRY_DEBUG
    G4cout << "INFO: We found material \'" << name << " previously defined."
	   << G4endl;
    G4cout << "         Checking if \'" << testname << "\' exists." << G4endl;
#endif
  }

  if (matrep > 99) {
    G4cerr << "ERROR: Too many materials with the same name. Exiting!"
	  << G4endl;
    abort();
  }
  
  fFlukaMaterials[testname] = this;
  if (name != testname)
    AddLowMat(testname);
  fIndex = fFlukaMaterials.size() + 2;
  fFlukaIndexMaterials[fIndex] = this;
}

FlukaMaterial::~FlukaMaterial() {
  delete fFlukaLowMat;
}

void FlukaMaterial::AddLowMat(const G4String& name) {
  fFlukaLowMat = new FlukaLowMat(name, this);
}


G4String FlukaMaterial::GetRealName() const {
  if (fFlukaLowMat)
    return fFlukaLowMat->GetName();
  return GetName();
}

std::ostream& FlukaMaterial::PrintMaterialsByName(std::ostream& os) {
  PrintHeader(os, "MATERIALS");
  for (FlukaMaterialsIterator i = fFlukaMaterials.begin(); 
       i != fFlukaMaterials.end(); 
       i++) {
    FlukaMaterial* flumat     = (*i).second;
    
    //if (flumat->GetZ()) //Skip materials that describe only compounds
    os << *flumat;
  }
  return os;
}

std::ostream& FlukaMaterial::PrintMaterialsByIndex(std::ostream& os) {
  PrintHeader(os, "MATERIALS");
  for (FlukaMaterialsIndexIterator i = fFlukaIndexMaterials.begin(); 
       i != fFlukaIndexMaterials.end(); 
       i++) {
    FlukaMaterial* flumat     = (*i).second;
    
    //if (flumat->GetZ()) //Skip materials that describe only compounds
    os << *flumat;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const FlukaMaterial& material){
  os << setw10 << "MATERIAL  ";

  os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
  G4double Z = G4double(material.GetZ());
  if (Z <= 0)
    os << setw10 << " ";
  else
    os << setw10 
       << setfixed
       << std::setprecision(1) 
       << Z;
  
  G4double A = material.GetA();
  if (A <= 0)
    os << setw10 << " ";
  else
    os << setw10 << std::setprecision(3)
       << A;

  G4double density = material.GetDensity();
  if (density <=0)
    density = 0.999;
  os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
  os << setw10 
     << setscientific
     << std::setprecision(3) 
     << density;

  os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
  os << setw10 
     << setfixed
     << std::setprecision(1) 
     << G4double(material.GetIndex());

  os << setw10 << " ";
  if (material.GetN())
    os << setw10 << G4double(material.GetN());
  else
    os << setw10 << " ";


  os << material.GetRealName().substr(0,8) << G4endl;

  if (material.GetLowMat() && material.GetZ() != 0)
    os << *(material.GetLowMat());

  return os;
}
