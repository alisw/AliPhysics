
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapMag.hh - Sara Vanini
//
// Wrapper for geometry tracking in magnetic field: returns magnetic 
// field values in a given position.
//
// modified 26/X/1998
// modified 18/XI/1999
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
/////////////////////////////////////////////////////////////////


#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "globals.hh"

void magfld(const G4double& pX, const G4double& pY, const G4double& pZ,
            G4double& cosBx, G4double& cosBy, G4double& cosBz, 
            G4double& Bmag, G4int& reg, G4int& idiscflag)

{
  //flag
#ifdef G4GEOMETRY_DEBUG
  G4cout<<"================== MAGFLD ================="<<G4endl;
#endif 
  
  //Geoinit Pointer
  FGeometryInit * ptrGeoInit=FGeometryInit::GetInstance();
  
  //get FieldManager, Field pointers for magnetic field handling
  G4FieldManager * pFieldMgr = ptrGeoInit->getFieldManager();
  const G4Field * ptrField = pFieldMgr->GetDetectorField();
  
  //compute field
  G4double point[3];
  point[0] = pX*10.;
  point[1] = pY*10.;
  point[2] = pZ*10.;
  G4double B[3];
  ptrField->GetFieldValue(point,B);
  G4double Bnor = sqrt(sqr(B[0])+sqr(B[1])+sqr(B[2]));
  if(Bnor) {
    cosBx = B[0]/Bnor;
    cosBy = B[1]/Bnor;
    cosBz = B[2]/Bnor;
  }
  else {
    cosBx = 0;
    cosBy = 0;
    cosBz = 1;
  }
  
  Bmag = Bnor/tesla;
  idiscflag = 0;
}




