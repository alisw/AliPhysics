// $Header$

#ifndef REVE_GridStepper_H
#define REVE_GridStepper_H

#include <Reve/Reve.h>

#include <TObject.h>

namespace Reve {

class ZTrans;

class GridStepper : public TObject
{
private:
  Int_t *ls[3], *ns[3]; //! Internal traversal variables.

  GridStepper(const GridStepper&);            // Not implemented
  GridStepper& operator=(const GridStepper&); // Not implemented

public: 
  enum StepMode_e { SM_XYZ, SM_YXZ, SM_XZY };
  StepMode_e Mode;      // Stepping mode, order of filling.

  Int_t   nx, ny, nz;   // Current positions during filling / traversal.
  Int_t   Nx, Ny, Nz;   // Number of slots in eaxh direction.
  Float_t Dx, Dy, Dz;   // Step size in each direction.
  Float_t Ox, Oy, Oz;   // Initial offset for each direction.

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
  
  Bool_t Step();

  void GetPosition(Float_t* p);

  void SetTrans(ZTrans* mx);
  void SetTransAdvance(ZTrans* mx);

  ClassDef(GridStepper, 1);
}; // end class GridStepper

}  // namespace Reve

#endif
