// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4G3Defaults
// -------------------
// Class stores the default G3 values of the kinetic energy cuts
// for particles and the control process parameters.

#ifndef TG4_G3_DEFAULTS_H
#define TG4_G3_DEFAULTS_H

#include "TG4G3CutVector.h"
#include "TG4G3ControlVector.h"
#include "TG4G3Control.h"
#include "TG4G3Cut.h"

#include <globals.hh>

class TG4G3Defaults
{
  public:
    TG4G3Defaults();      
    // --> protected
    // TG4G3Defaults(const TG4G3Defaults& right);
    virtual ~TG4G3Defaults();

    // static access methods
    static TG4G3Defaults* Instance();

    // methods
    G4bool IsDefaultCut(TG4G3Cut cut, G4double value) const; 
    G4bool IsDefaultControl(TG4G3Control control, TG4G3ControlValue value) const; 

    // get methods
    G4double CutValue(G4int cut) const;   
    TG4G3ControlValue ControlValue(G4int control) const; 
      
  protected:
    TG4G3Defaults(const TG4G3Defaults& right);      

    // operators
    TG4G3Defaults& operator=(const TG4G3Defaults& right);

  private:
    // static data members  
    static TG4G3Defaults*  fgInstance;      //this instance
    
    // data members
    TG4G3CutVector     fCutVector;    // vector of default cut values       
    TG4G3ControlVector fControlVector;// vector of default control values       
};     

// inline methods

inline TG4G3Defaults* TG4G3Defaults::Instance() 
{ return fgInstance; }

#endif //ALI_G3_DEFAULTS_H
