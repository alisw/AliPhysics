#ifndef ALISIGNAL_H
#define ALISIGNAL_H
///////////////////////////////////////////////////////////////////////////
// Class AliSignal
// Handling of ALICE (extrapolated) signals.
//
// Note :
// ------
// Signal positions (r) and reference frames (f) are specified via
// SetPosition(r,f) under the following conventions :
//
// f="car" ==> r is Cartesian   (x,y,z)
// f="sph" ==> r is Spherical   (r,theta,phi)
// f="cyl" ==> r is Cylindrical (rho,phi,z)
//
// All angles are in radians.
//
// Example :
// ---------
//
// AliSignal s;
// Float_t pos[3]={-1,25,7};
// Float_t signal=120.8;
// s.SetPosition(pos,"car");
// s.SetSignal(signal);
// Float_t loc[3];
// s.GetPosition(loc,"sph");
// Float_t adc=s.GetSignal();
//
//--- NvE 23-jan-1999 UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliPosition.h"

class AliSignal : public TObject,public AliPosition
{
 public:
  AliSignal();                         // Default constructor
  ~AliSignal();                        // Destructor
  virtual void SetSignal(Float_t sig); // Store signal value
  virtual Float_t GetSignal();         // Provide signal value
  virtual void Reset();                // Reset all values to 0

 protected:
  Float_t fSignal; // Signal value

 ClassDef(AliSignal,1) // Class definition to enable ROOT I/O
};
#endif
