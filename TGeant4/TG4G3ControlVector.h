// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4G3ControlVector
// ------------------------
// Vector of control process values
// with convenient set/get methods.

#ifndef TG4_G3_CONTROL_VECTOR_H
#define TG4_G3_CONTROL_VECTOR_H

#include "TG4Globals.h"
#include "TG4G3Control.h"

#include <g4std/vector>

class TG4G3CutVector;

class G4VProcess;

class TG4G3ControlVector
{
  typedef G4std::vector<TG4G3ControlValue> TG4ControlValueVector;

  public:
    TG4G3ControlVector();
    TG4G3ControlVector(const TG4G3ControlVector& right);
    virtual ~TG4G3ControlVector();
    
    // operators
    TG4G3ControlVector& operator=(const TG4G3ControlVector& right);
    TG4G3ControlValue operator[](G4int index) const;

    // static methods 
    static TG4G3Control      GetControl(const G4String& controlName);
    static const G4String&   GetControlName(TG4G3Control control);
    static TG4G3ControlValue GetControlValue(G4int value, 
                                             TG4G3Control control);
    static TG4G3ControlValue GetControlValue(G4double value, 
                                             TG4G3Control control);

    // set methods
    G4bool SetControl(TG4G3Control control, TG4G3ControlValue controlValue,
                    TG4G3CutVector& cuts);
    void   SetG3Defaults();
    G4bool Update(const TG4G3ControlVector& vector);
    
    // methods
    G4String Format() const;
    void Print() const;
    
    // get methods
    TG4G3ControlValue GetControlValue(G4VProcess* process) const; 
    TG4G3ControlValue GetControlValue(TG4G3Control control) const; 
    G4bool IsControl() const;

  private:
    // static methods 
    static void FillControlNameVector();
  
    // static data members
    static TG4StringVector  fgControlNameVector; //vector of control parameters
                                                 //names

    // data members
    TG4ControlValueVector   fControlVector; //vector of control process values
};

#endif //TG4_G3_CONTROL_VECTOR_H



