#include "FlukaLowMat.hh"
#include "WrapUtils.hh"

FlukaLowMat::FlukaLowMat(const G4String& name, FlukaMaterial* fmat):
  fName(name),
  fFlukaMaterial(fmat){
}

G4std::ostream& operator<<(G4std::ostream& os, const FlukaLowMat& flowmat){
  os << setw10 << "LOW-MAT   ";
  os.setf(static_cast<G4std::ios::fmtflags>(0),G4std::ios::floatfield);
  os << setw10 << setfixed << G4std::setprecision(1) 
     << G4double(flowmat.GetIndex());
  os << setw10 << " " 
     << setw10 << " " 
     << setw10 << " " 
     << setw10 << " " 
     << setw10 << " ";
  os << flowmat.GetFlukaMaterial()->GetName() << G4endl;

  return os;
}
