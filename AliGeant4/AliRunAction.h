// $Id$
// Category: run
//
// Author: I. Hrivnacova
//
// Class AliRunAction
// ------------------
// Class that defines actions at the beginning and the end of run.

#ifndef ALI_RUN_ACTION_H
#define ALI_RUN_ACTION_H

#include "AliVerbose.h"

#include <G4UserRunAction.hh>
#include <globals.hh>

class G4Timer;
    // in order to avoid the odd dependency for the
    // times system function this declaration must be the first

class AliSDConstruction;
class G4Run;

class AliRunAction : public G4UserRunAction,
                     public AliVerbose
{
  public:
    AliRunAction();
    // --> protected
    // AliRunAction(const AliRunAction& right);
    virtual ~AliRunAction();

    // methods
    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);

  protected:
    AliRunAction(const AliRunAction& right);

    // operators
    AliRunAction& operator=(const AliRunAction& right);

  private:
    // methods
    AliSDConstruction* GetSDConstruction() const;

    // data members
    G4Timer*  fTimer; //G4Timer
    G4int     fRunID; //run ID
};

#endif //ALI_RUN_ACTION_H

