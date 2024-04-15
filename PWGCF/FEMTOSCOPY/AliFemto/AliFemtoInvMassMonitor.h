///
/// \file AliFemtoInvMassMonitor.h
///
/// \class AliFemtoInvMassMonitor
/// \brief A cut monitor for particles storing the difference between
///        reconstructed and true momentum
///

#ifndef AliFemtoInvMassMonitor_hh
#define AliFemtoInvMassMonitor_hh

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
//class AliFemtoPair; // Gael 12/04/02
class TH1D;
class TH2D;
class TList;

#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"
#include "AliFemtoPair.h"


class AliFemtoInvMassMonitor : public AliFemtoCutMonitor {
public:
  AliFemtoInvMassMonitor();
  AliFemtoInvMassMonitor(const char *aName, double m1, double m2);
  AliFemtoInvMassMonitor(const AliFemtoInvMassMonitor &aCut);
  virtual ~AliFemtoInvMassMonitor();

  AliFemtoInvMassMonitor& operator=(const AliFemtoInvMassMonitor& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack) {AliFemtoCutMonitor::Fill(aTrack);}
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoXi* aXi) {AliFemtoCutMonitor::Fill(aXi);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* pair);
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2) {AliFemtoCutMonitor::Fill(aCollection1, aCollection2);}

  void Write();

  virtual TList *GetOutputList();

private:
  TH1D *fInvMass;     // Rapidity vs. Pt monitor
  double fM1;
  double fM2;
};

#endif
