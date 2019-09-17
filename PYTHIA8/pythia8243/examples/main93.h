#include "Pythia8/Pythia.h"

class RootTrack  {
public:

  bool init(Pythia8::Particle& p) {
    if (p.isFinal()) {
      phi = p.phi(), eta = p.eta(), y = p.y();
      pT = p.pT(), pid = p.id();
      return true;
    }
    return false;
  }

  double phi, eta, y, pT;
  int pid;
};

class RootEvent {
public:

  bool init(Pythia8::Info* infoPtr) {
    tracks.clear();
    // An event level cut on eg. impact parameter, number of
    // MPIs etc. can be implemented here.
    // if () return false;
    weight = infoPtr->weight();
    return true;
  }

  double weight;
  std::vector<RootTrack> tracks;

};
