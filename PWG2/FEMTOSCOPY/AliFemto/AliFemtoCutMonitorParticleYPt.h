////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticleYPt - the cut monitor for particles to study  ///
/// the difference between reconstructed and true momentum    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticleYPt_hh
#define AliFemtoCutMonitorParticleYPt_hh

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
class TH1D;
class TH2D;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorParticleYPt : public AliFemtoCutMonitor{
  
public:
  AliFemtoCutMonitorParticleYPt();
  AliFemtoCutMonitorParticleYPt(const char *aName, float aMass);
  AliFemtoCutMonitorParticleYPt(const AliFemtoCutMonitorParticleYPt &aCut);
  virtual ~AliFemtoCutMonitorParticleYPt();

  AliFemtoCutMonitorParticleYPt& operator=(const AliFemtoCutMonitorParticleYPt& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoTrack* aTrack);
  void Write();

  virtual TList *GetOutputList();

private:
  TH2D *fYPt;    // Rapidity vs. Pt monitor
  float fMass;   // Mass hypothesis
};

#endif
