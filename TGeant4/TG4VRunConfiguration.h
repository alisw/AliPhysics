// $Id$
// Category: run
//
// Abstract class that takes care of creating  all user defined classes 
// that will be initialized and managed by Geant4 kernel (G4RunManager).
// It has one pure virtual method CreateUserConfiguration()
// that has to be be implemented by derived class.

#ifndef TG4V_RUN_CONFIGURATION_H
#define TG4V_RUN_CONFIGURATION_H

class G4VUserDetectorConstruction;
class G4VModularPhysicsList;
class G4VUserPrimaryGeneratorAction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4UserStackingAction;
class G4RunManager;

class TG4VRunConfiguration
{
  public:
    TG4VRunConfiguration();
    // --> protected
    // TG4VRunConfiguration(const TG4VRunConfiguration& right);
    virtual ~TG4VRunConfiguration();

    // methods
    void ConfigureRunManager(G4RunManager* runManager);

    // get methods
    G4VModularPhysicsList* GetPhysicsList() const;

  protected:
    TG4VRunConfiguration(const TG4VRunConfiguration& right);

    // operators
    TG4VRunConfiguration& operator=(const TG4VRunConfiguration& right);

    // methods
    virtual void CreateUserConfiguration() = 0;

    // data members
    G4VUserDetectorConstruction*    fDetectorConstruction; //det construction
    G4VModularPhysicsList*          fPhysicsList;          //physics list
    G4VUserPrimaryGeneratorAction*  fPrimaryGenerator;     //primary generator
    G4UserRunAction*                fRunAction;            //run action
    G4UserEventAction*              fEventAction;          //event action
    G4UserTrackingAction*           fTrackingAction;       //tracking action
    G4UserSteppingAction*           fSteppingAction;       //stepping action
    G4UserStackingAction*           fStackingAction;       //stacking action
};

#endif //TG4V_RUN_CONFIGURATION_H

