#include "FlukaLowMat.hh"
#include "WrapUtils.hh"

FlukaLowMat::FlukaLowMat(const G4String& name, FlukaMaterial* fmat):
  fName(name),
  fFlukaMaterial(fmat){
}

std::ostream& operator<<(std::ostream& os, const FlukaLowMat& flowmat){
  os << setw10 << "LOW-MAT   ";
  os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
  os << setw10 << setfixed << std::setprecision(1) 
     << G4double(flowmat.GetIndex());
  os << setw10 << " " 
     << setw10 << " " 
     << setw10 << " " 
     << setw10 << " " 
     << setw10 << " ";
  os << flowmat.GetFlukaMaterial()->GetName() << G4endl;

  return os;
}
