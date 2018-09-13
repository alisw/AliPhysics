////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticlePIDBeta - the cut monitor for particles to study   ///
/// various aspects of the PID determination                                 ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticlePIDBeta_hh
#define AliFemtoCutMonitorParticlePIDBeta_hh

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
#include "AliFemtoCutMonitorParticlePID.h"
#include "AliFemtoCutMonitor.h"
class AliFemtoCutMonitorParticlePIDBeta : public AliFemtoCutMonitorParticlePID{

public:
  AliFemtoCutMonitorParticlePIDBeta();
  AliFemtoCutMonitorParticlePIDBeta(const char *aName, Int_t aTOFParticle, Double_t yTOFTimeMin= -4000.0, Double_t yTOFTimeMax=4000.0);
  AliFemtoCutMonitorParticlePIDBeta(const AliFemtoCutMonitorParticlePIDBeta &aCut);
  virtual ~AliFemtoCutMonitorParticlePIDBeta();
  AliFemtoCutMonitorParticlePIDBeta& operator=(const AliFemtoCutMonitorParticlePIDBeta& aCut);

  virtual void Fill(const AliFemtoTrack* aTrack);



  void Write();
  virtual TList *GetOutputList();

protected:

  TH2D *fBeta;           ///< Beta
  TH2D *fMass;           ///< Mass
  TH1D *fDifference;           ///< Mass

};

#endif
