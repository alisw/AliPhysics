// $Id$
// Category: run
//
// Class that defines actions at the beginning and the end of run.

#ifndef ALI_RUN_ACTION_H
#define ALI_RUN_ACTION_H

#include <globals.hh>
#include <G4UserRunAction.hh>

class G4Timer;
    // in order to avoid the odd dependency for the
    // times system function this declaration must be the first

class AliRunActionMessenger;
class G4Run;

class AliRunAction : public G4UserRunAction
{
  public:
    AliRunAction();
    // --> protected
    // AliRunAction(const AliRunAction& right);
    virtual ~AliRunAction();

    // methods
    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);

    // set methods
    void SetVerboseLevel(G4int level);

    // get methods
    G4int GetVerboseLevel() const;

  protected:
    AliRunAction(const AliRunAction& right);

    // operators
    AliRunAction& operator=(const AliRunAction& right);

  private:
    // data members
    AliRunActionMessenger*  fMessenger;    //messenger 
    G4Timer*                fTimer;        //G4Timer
    G4int                   fRunID;        //run ID
    G4int                   fVerboseLevel; //verbose level
};

// inline methods

inline void AliRunAction::SetVerboseLevel(G4int level)
{ fVerboseLevel = level; }

inline G4int AliRunAction::GetVerboseLevel() const
{ return fVerboseLevel; }

#endif //ALI_RUN_ACTION_H

