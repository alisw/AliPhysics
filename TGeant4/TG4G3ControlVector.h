// $Id$
// Category: global
//
// Vector of control process values
// with convenient set/get methods

#ifndef TG4_G3_CONTROL_VECTOR_H
#define TG4_G3_CONTROL_VECTOR_H

#include "TG4Globals.h"
#include "TG4G3Control.h"

class G4VProcess;

class TG4G3ControlVector
{
  typedef G4RWTValOrderedVector<TG4G3ControlValue> TG4ControlValueVector;

  public:
    TG4G3ControlVector();
    TG4G3ControlVector(const TG4G3ControlVector& right);
    virtual ~TG4G3ControlVector();
    
    // operators
    TG4G3ControlVector& operator=(const TG4G3ControlVector& right);
    G4double operator[](G4int index) const;

    // set methods
    void SetG3Control(TG4G3Control control, G4double controlValue);
    void SetG3Defaults();
    
    // get methods
    G4int GetControl(G4VProcess* process) const; 

  private:
    // data members
    TG4ControlValueVector*  fControlVector; //vector of control process values
};

#endif //TG4_G3_CONTROL_VECTOR_H



