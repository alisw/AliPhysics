// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliEventAction
// --------------------
// Class that defines actions at the beginning and the end of event.

#ifndef ALI_EVENT_ACTION_H
#define ALI_EVENT_ACTION_H 

#include "AliVerbose.h"
#include "AliEventActionMessenger.h"

#include <G4UserEventAction.hh>
#include <globals.hh>

class G4Timer;
    // in order to avoid the odd dependency for the
    // times system function this declaration must be the first
class G4Event;

class AliEventAction : public G4UserEventAction,
                       public AliVerbose
{
  public:
    AliEventAction();
    // --> protected
    // AliEventAction(const AliEventAction& right);
    virtual ~AliEventAction();
    
    // methods
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    
    // set methods
    void SetDrawFlag(G4String drawFlag);
    
    // get methods
    G4String GetDrawFlag() const;
    
  protected:
    AliEventAction(const AliEventAction& right);

    // operators
    AliEventAction& operator=(const AliEventAction& right);

  private:
    // methods 
    void DisplayEvent(const G4Event* event) const;
  
    // data members
    AliEventActionMessenger   fMessenger;    //messenger
    G4Timer*                  fTimer;        //G4Timer
    G4String                  fDrawFlag;     //control drawing of the event
};

// inline methods

inline void AliEventAction::SetDrawFlag(G4String drawFlag)
{ fDrawFlag = drawFlag; }

inline G4String AliEventAction::GetDrawFlag() const
{ return fDrawFlag; }

#endif //ALI_EVENT_ACTION_H

    
