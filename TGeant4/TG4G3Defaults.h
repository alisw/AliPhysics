// $Id$
// Category: physics
//
// Class stores the default G3 values of the kinetic energy cuts
// for particles and the control process parameters

#ifndef TG4_G3_DEFAULTS_H
#define TG4_G3_DEFAULTS_H

#include "TG4G3Cut.h"
#include "TG4G3Control.h"

#include <globals.hh>

class TG4G3Defaults
{
  public:
    // --> protected
    // TG4G3Defaults();
    virtual ~TG4G3Defaults();

    // static methods
    static G4bool IsDefaultCut(TG4G3Cut cut, G4double value); 
    static G4bool IsDefaultControl(TG4G3Control control, G4double value); 

    // static get methods
       // precision tolerance
    static G4double CutTolerance();
    static G4double CutValue(G4int cut);   
    static TG4G3ControlValue ControlValue(G4int control); 
      
  protected:
    TG4G3Defaults();      
      // only static data members and methods

  private:
    // static data members  
      // precision tolerance
    static const G4double fgkCutTolerance; //tolerance for checking
                                          //cut values with default
       // kinetic energy cuts
    static const G4double fgkCUTGAM;   //CUTGAM default value 
    static const G4double fgkCUTELE;   //CUTELE default value
    static const G4double fgkCUTNEU;   //CUTNEU default value
    static const G4double fgkCUTHAD;   //CUTHAD default value
    static const G4double fgkCUTMUO;   //CUTMUO default value
    static const G4double fgkBCUTE;    //BCUTE default value
    static const G4double fgkBCUTM;    //BCUTM default value
    static const G4double fgkDCUTE;    //DCUTE default value
    static const G4double fgkDCUTM;    //DCUTM default value
    static const G4double fgkPPCUTM;   //PPCUTM default value
       // physics processes controls
    static const TG4G3ControlValue fgkPAIR; //PAIR default value   
    static const TG4G3ControlValue fgkCOMP; //COMP default value   
    static const TG4G3ControlValue fgkPHOT; //PHOT default value   
    static const TG4G3ControlValue fgkPFIS; //PFIS default value   
    static const TG4G3ControlValue fgkDRAY; //DRAY default value  
    static const TG4G3ControlValue fgkANNI; //ANNI default value   
    static const TG4G3ControlValue fgkBREM; //BREM default value  
    static const TG4G3ControlValue fgkHADR; //HADR default value   
    static const TG4G3ControlValue fgkMUNU; //MUNU default value   
    static const TG4G3ControlValue fgkDCAY; //DCAY default value  
    static const TG4G3ControlValue fgkLOSS; //LOSS default value  
    static const TG4G3ControlValue fgkMULS; //MULS default value  
};     

// inline methods

inline G4double TG4G3Defaults::CutTolerance() { return fgkCutTolerance; }

#endif //ALI_G3_DEFAULTS_H
