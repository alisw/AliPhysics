//  $Id$
// Category: geometry
//
// Uniform magnetic field.
// According to:
// Id: ExN02MagneticField.hh,v 1.1 1999/01/07 16:05:47 gunter Exp 
// GEANT4 tag Name: geant4-00-01 

#ifndef ALI_MAGNETIC_FIELD_H
#define ALI_MAGNETIC_FIELD_H

#include <G4UniformMagField.hh>

class G4FieldManager;

class AliMagneticField: public G4UniformMagField
{
  public:
    AliMagneticField();                   //  A zero field
    AliMagneticField(G4ThreeVector fied); //  The value of the field
    AliMagneticField(const AliMagneticField& right);
    virtual ~AliMagneticField();  
      
    // operators
    AliMagneticField& operator=(const AliMagneticField& right);

    // set methods
    void SetFieldValue(G4ThreeVector fieldVector);
    void SetFieldValue(G4double      fieldValue);
              // Set the field to (0, 0, fieldValue)

    // get methods
    G4ThreeVector GetConstantFieldValue();

  protected:
    // Find the global Field Manager
    G4FieldManager* GetGlobalFieldManager(); 
};

#endif //ALI_MAGNETIC_FIELD_H
