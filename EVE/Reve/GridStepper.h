// $Header$

#ifndef REVE_GridStepper_H
#define REVE_GridStepper_H

#include <Reve/Reve.h>

#include <TObject.h>

namespace Reve {

class ZTrans;

class GridStepper : public TObject
{
  Int_t *ls[3], *ns[3];
private:
  GridStepper(const GridStepper&);            // Not implemented
  GridStepper& operator=(const GridStepper&); // Not implemented

public: 
  enum StepMode_e { SM_XYZ, SM_YXZ, SM_XZY };
  StepMode_e Mode; 

  Int_t   nx, ny, nz;
  Int_t   Nx, Ny, Nz;
  Float_t Dx, Dy, Dz;
  Float_t Ox, Oy, Oz;

  GridStepper(Int_t sm=SM_XYZ);
  virtual ~GridStepper() {}

  void Reset();
  void Subtract(GridStepper& s);
  void SetNs(Int_t nx, Int_t ny, Int_t nz=1)
  { Nx = nx; Ny = ny; Nz = nz; }
  void SetDs(Float_t dx, Float_t dy, Float_t dz=0)
  { Dx = dx; Dy = dy; Dz = dz; }
  void SetOs(Float_t ox, Float_t oy, Float_t oz=0)
  { Ox = ox; Oy = oy; Oz = oz; }
  
  bool Step();

  void GetPosition(Float_t* p);

  void SetTrans(ZTrans* mx);
  void SetTransAdvance(ZTrans* mx);

  ClassDef(GridStepper, 1);
}; // end class GridStepper

}  // namespace Reve

#endif
