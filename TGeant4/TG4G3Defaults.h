// $Id$
// Category: physics
//
// Class stores the default G3 values of the kinetic energy cuts
// for particles and the control process flags parameters

#ifndef TG4_G3DEFAULTS_H
#define TG4_G3DEFAULTS_H

#include "TG3Cut.h"
#include "TG3Flag.h"

#include <globals.hh>

class TG4G3Defaults
{
  public:
    // --> protected
    // TG4G3Defaults();
    virtual ~TG4G3Defaults();

    // static methods
    static G4bool IsDefaultCut(TG3Cut g3Cut, G4double value); 
    static G4bool IsDefaultFlag(TG3Flag g3Flag, G4double value); 

    // static get methods
       // precision tolerance
    static G4double CutTolerance();
    static G4double CutValue(G4int g3Cut);   
    static TG3FlagValue FlagValue(G4int g3Flag); 
      
  protected:
    TG4G3Defaults();      
      // only static data members and methods

  private:
    // static data members  
      // precision tolerance
    static const G4double fgCutTolerance; //tolerance for checking
                                          //cut values with default
       // kinetic energy cuts
    static const G4double fgCUTGAM;   //CUTGAM default value 
    static const G4double fgCUTELE;   //CUTELE default value
    static const G4double fgCUTNEU;   //CUTNEU default value
    static const G4double fgCUTHAD;   //CUTHAD default value
    static const G4double fgCUTMUO;   //CUTMUO default value
    static const G4double fgBCUTE;    //BCUTE default value
    static const G4double fgBCUTM;    //BCUTM default value
    static const G4double fgDCUTE;    //DCUTE default value
    static const G4double fgDCUTM;    //DCUTM default value
    static const G4double fgPPCUTM;   //PPCUTM default value
       // physics processes flags
    static const TG3FlagValue fgPAIR; //PAIR default value   
    static const TG3FlagValue fgCOMP; //COMP default value   
    static const TG3FlagValue fgPHOT; //PHOT default value   
    static const TG3FlagValue fgPFIS; //PFIS default value   
    static const TG3FlagValue fgDRAY; //DRAY default value  
    static const TG3FlagValue fgANNI; //ANNI default value   
    static const TG3FlagValue fgBREM; //BREM default value  
    static const TG3FlagValue fgHADR; //HADR default value   
    static const TG3FlagValue fgMUNU; //MUNU default value   
    static const TG3FlagValue fgDCAY; //DCAY default value  
    static const TG3FlagValue fgLOSS; //LOSS default value  
    static const TG3FlagValue fgMULS; //MULS default value  
};     

// inline methods

inline G4double TG4G3Defaults::CutTolerance() { return fgCutTolerance; }

#endif //ALI_G3DEFAULTS_H
