// $Id$
// Category: event
//
// Class that ensures calling sensitive detector
// before track starts stepping.
// It also takes care of setting step status (kVertex)
// and passing G4Track to TG4StepManager.

#ifndef TG4_TRACKING_ACTION_H
#define TG4_TRACKING_ACTION_H

#include <G4UserTrackingAction.hh>

class G4Track;

class TG4TrackingAction : public G4UserTrackingAction 
{
  public:
    TG4TrackingAction();
    // --> protected
    // TG4TrackingAction(const TG4TrackingAction& right);
    virtual ~TG4TrackingAction();
   
    // methods
    virtual void PreTrackingAction(const G4Track* aTrack) {;}
    virtual void PostTrackingAction(const G4Track* aTrack) {;}
                  // the following methods should not
		  // be overwritten in a derived class
    virtual void PreUserTrackingAction(const G4Track* aTrack);
    virtual void PostUserTrackingAction(const G4Track* aTrack);


  protected:
    TG4TrackingAction(const TG4TrackingAction& right);

    // operators
    TG4TrackingAction& operator=(const TG4TrackingAction& right);
};

#endif //TG4_TRACKING_ACTION_H
