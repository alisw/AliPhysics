// $ Id:$
// Category: global

// Concrete class of G4UIcommand. 
// The command defined by this class takes up to three string values.
// General information of G4UIcommand is given in G4UIcommand.hh.

#ifndef TG4_UI_CMD_WITH_A_COMPLEX_STRING_H
#define TG4_UI_CMD_WITH_A_COMPLEX_STRING_H

#include <G4UIcommand.hh>

class TG4UICmdWithAComplexString : public G4UIcommand
{
  public:
    TG4UICmdWithAComplexString(G4String commandPath, G4UImessenger* messenger);
    virtual ~TG4UICmdWithAComplexString();			       

    // set methods
    void SetParameterName(G4String name, G4bool omittable);
    void SetDefaultValue(G4String defVal);

    // get methods
    G4String GetNewStringValue(G4String paramString);
};

#endif //TG4_UI_CMD_WITH_A_COMPLEX_STRING_H
