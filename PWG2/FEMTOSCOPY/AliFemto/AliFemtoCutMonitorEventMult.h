////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorEventMult - the cut monitor for events to study        ///
/// the multiplicity distribution of events                                  ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorEventMult_hh
#define AliFemtoCutMonitorEventMult_hh

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair; 
class TH1D;
class TH2D;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorEventMult : public AliFemtoCutMonitor{
  
 public:
  AliFemtoCutMonitorEventMult();
  AliFemtoCutMonitorEventMult(const char *aName);
  AliFemtoCutMonitorEventMult(const AliFemtoCutMonitorEventMult &aCut);
  virtual ~AliFemtoCutMonitorEventMult();

  AliFemtoCutMonitorEventMult& operator=(const AliFemtoCutMonitorEventMult& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent);
  virtual void Fill(const AliFemtoTrack* aTrack) {AliFemtoCutMonitor::Fill(aTrack);}
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}

  void SetReadMC(Bool_t mc);
  void AdditionalMultHistsOn(Bool_t addhists);
  void Write();

  virtual TList *GetOutputList();

 private:
  TH1D *fEvMult;     // Multiplicity distribution
  TH1D *fNormEvMult; // Normalized event multiplicity distribution
  TH1D *fSPDMult;    // SPD tracklet multiplicity
  TH2D *fMultSumPt;  // Event total pT vs. multiplicity

  Bool_t freadMC;     // If true - add only one histogram to the output
  Bool_t faddhists;   // If true - add only additional multiplicity histograms
  
  TH1D *fEstimateITSTPC;     // Multiplicity estimate ITS+TPC
  TH1D *fEstimateTracklets;  // Multiplicity estimate Tracklets
  TH1D *fEstimateITSPure;    // Multiplicity estimate ITS Pure

  TH2D *fEst1Est2;           // ITS+TPC vs Tracklets
  TH2D *fEst1Est3;           // ITS+TPC vs ITS Pure
  TH2D *fEst2Est3;           // Tracklets vs ITS Pure
  TH2D *fEst1Norm;           // ITS+TPC vs Normalized
  TH2D *fEst2Norm;           // Tracklets vs Normalized
  TH2D *fEst3Norm;           // ITS Pure vs Normalized

};

#endif
