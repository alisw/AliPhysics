
/// \brief A cut monitor for particles to store the distribution of kT and PID
/// tested only for MC data
///

#ifndef AliFemtoCutMonitorPairKT_hh
#define AliFemtoCutMonitorPairKT_hh

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

class AliFemtoCutMonitorPairKT : public AliFemtoCutMonitor {
public:
  AliFemtoCutMonitorPairKT();
  AliFemtoCutMonitorPairKT(const char *aName, float lowLimit, float highLimit);
  AliFemtoCutMonitorPairKT(const AliFemtoCutMonitorPairKT &aCut);
  virtual ~AliFemtoCutMonitorPairKT();

  AliFemtoCutMonitorPairKT& operator=(const AliFemtoCutMonitorPairKT& aCut);

  virtual AliFemtoString Report();
  void Write();
  virtual TList *GetOutputList();
  virtual void Fill(const AliFemtoEvent* aEvent) { AliFemtoCutMonitor::Fill(aEvent); }
  virtual void Fill(const AliFemtoTrack* aTrack) { AliFemtoCutMonitor::Fill(aTrack); }
  virtual void Fill(const AliFemtoV0* aV0) { AliFemtoCutMonitor::Fill(aV0); }
  virtual void Fill(const AliFemtoXi* aXi) {AliFemtoCutMonitor::Fill(aXi);}
  virtual void Fill(const AliFemtoKink* aKink) { AliFemtoCutMonitor::Fill(aKink); }
  virtual void Fill(const AliFemtoPair* aPair);
  virtual void Fill(const AliFemtoParticleCollection* aCollection) { AliFemtoCutMonitor::Fill(aCollection); }
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection) { AliFemtoCutMonitor::Fill(aEvent, aCollection); }
  virtual void Fill(const AliFemtoParticleCollection* aCollection1, const AliFemtoParticleCollection* aCollection2) { AliFemtoCutMonitor::Fill(aCollection1, aCollection2); }


private:
  TH1D *fkT;
  TH1D *fpTsum;  
  TH1D *fpTdiff;
  TH1D *fQinv;   
  TH1D *fMinv;    
  
};

#endif
