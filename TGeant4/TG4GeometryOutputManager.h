// $Id$
// Category: geometry
//
// This class provides methods for writing
// tokens (ASCII form of Geant3 routines calls)
// into the specified output file

#ifndef TG4_GEOMETRY_OUTPUT_MANAGER_H
#define TG4_GEOMETRY_OUTPUT_MANAGER_H

#include "TG4Globals.h"

#include <globals.hh>
#include <g4std/fstream>

class TG4GeometryOutputManager
{
  public:
    TG4GeometryOutputManager();
    virtual ~TG4GeometryOutputManager();

    // methods
    void OpenFile(G4String filePath);
    void CloseFile();

    void WriteGsvolu(G4String name, G4String shape, G4int nmed, G4double* par,
              G4int npar);

    void WriteGspos(G4String name, G4int num, G4String moth, G4double x,
              G4double y, G4double z, G4int irot, G4String only);

    void WriteGsposp(G4String name, G4int num, G4String moth, G4double x,
              G4double y, G4double z, G4int irot, G4String only,
              G4double* Rpar, G4int npar);

    void WriteGsrotm(G4int irot, G4double theta1, G4double phi1,
              G4double theta2, G4double phi2, G4double theta3, G4double phi3);

    //void WriteGsatt(G4String name, G4String attr, G4int ival);

    void WriteGsdvn(G4String vname, G4String vmoth, G4int ndiv, G4int iaxis);

    void WriteGsdvt(G4String name, G4String moth, G4double Step, G4int iaxis,
              G4int numed, G4int ndvmx);

    void WriteGsdvx(G4String name, G4String moth, G4int ndiv, G4int iaxis,
              G4double Step, G4double c0, G4int numed, G4int ndvmx);

    void WriteGsdvn2(G4String name, G4String moth, G4int ndiv, G4int iaxis,
              G4double c0, G4int numed);

    void WriteGsdvt2(G4String name, G4String moth, G4double Step, G4int iaxis,
              G4double c0, G4int numed, G4int ndvmx);

    void WriteGsmate(G4int imate, G4String name, G4double a, G4double z,
              G4double dens, G4double radl, G4int nwbf, G4double* ubuf);

    void WriteGsmixt(G4int imate, G4String name, G4double* a, G4double* z,
              G4double dens, G4int nlmat, G4double* wmat);

    void WriteGstmed(G4int itmed, G4String name, G4int nmat, G4int isvol,
              G4int ifield, G4double fieldm, G4double tmaxfd,
              G4double stemax, G4double deemax, G4double epsil,
              G4double stmin, G4double* par, G4int npar);

    void WriteGstpar(G4int itmed, G4String par, G4double parval);

    void WriteGgclos();

  private:
    // data members
    G4std::ofstream  fOutFile; //output file
};

#endif //TG4_GEOMETRY_OUTPUT_MANAGER_H

