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
       // physics processes flags
    static const TG3FlagValue fgkPAIR; //PAIR default value   
    static const TG3FlagValue fgkCOMP; //COMP default value   
    static const TG3FlagValue fgkPHOT; //PHOT default value   
    static const TG3FlagValue fgkPFIS; //PFIS default value   
    static const TG3FlagValue fgkDRAY; //DRAY default value  
    static const TG3FlagValue fgkANNI; //ANNI default value   
    static const TG3FlagValue fgkBREM; //BREM default value  
    static const TG3FlagValue fgkHADR; //HADR default value   
    static const TG3FlagValue fgkMUNU; //MUNU default value   
    static const TG3FlagValue fgkDCAY; //DCAY default value  
    static const TG3FlagValue fgkLOSS; //LOSS default value  
    static const TG3FlagValue fgkMULS; //MULS default value  
};     

// inline methods

inline G4double TG4G3Defaults::CutTolerance() { return fgkCutTolerance; }

#endif //ALI_G3DEFAULTS_H
