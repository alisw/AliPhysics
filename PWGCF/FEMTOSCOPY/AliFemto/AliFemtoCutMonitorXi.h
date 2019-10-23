/// \class AliFemtoCutMonitorV0
/// \brief AliFemtoCutMonitorXi -

#ifndef AliFemtoCutMonitorXi_H
#define AliFemtoCutMonitorXi_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoXi;
class AliFemtoKink;
class AliFemtoPair; // Gael 12/04/02
class TH1F;
class TH1D;
class TH2D;
class TList;
#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorXi : public AliFemtoCutMonitor{

public:
  AliFemtoCutMonitorXi();
  AliFemtoCutMonitorXi(const char *aName);
  AliFemtoCutMonitorXi(const AliFemtoCutMonitorXi &aCut);
  virtual ~AliFemtoCutMonitorXi();

  AliFemtoCutMonitorXi& operator=(const AliFemtoCutMonitorXi& aCut);

  virtual AliFemtoString Report();
  virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
  virtual void Fill(const AliFemtoTrack* aTrack){AliFemtoCutMonitor::Fill(aTrack);}
  virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
  virtual void Fill(const AliFemtoXi* axi);
  virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
  virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
  virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
  {AliFemtoCutMonitor::Fill(aEvent, aCollection);}
  virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2) {AliFemtoCutMonitor::Fill(aCollection1, aCollection2);}
  void Write();

  virtual TList *GetOutputList();

private:
  TH1F *fLambdaMass;     ///< Mass assuming lambda hypothesis
  TH1F *fAntiLambdaMass; ///< Mass assuming antilambda hypothesis
  TH1F *fK0ShortMass;    ///< Mass assuming k-short hypothesis
  TH1F *fDcaDaughters;   ///< DCA of v0 daughters at Decay vertex
  TH1F *fDcaV0ToPrimVertex;///< DCA of v0 to primary vertex
  TH1F *fDcaPosToPrimVertex;
  TH1F *fDcaNegToPrimVertex;
  TH1F *fCosPointingAngle;
  TH1F *fDecayLength;
  TH1F *fEtaV0;
  TH1F *fPtV0;
  TH1F *fPtPosDaughter;
  TH1F *fPtNegDaughter;

  TH2D *fdEdxPosDaughter;
  TH2D *fdEdxNegDaughter;
  TH2D *fTOFtimePosDaughter;
  TH2D *fTOFtimeNegDaughter;

  TH2D *fMINVvsPt;

  TH1D *fnsigmaPosL;
  TH1D *fnsigmaNegL;
  TH1D *fnsigmaPosAL;
  TH1D *fnsigmaNegAL;

  TH1D *fParticleOrigin; ///< particle origin from MC
  TH1D *fParticleId;     ///< true particle identification from MC


  TH1F *fXiMass;
  TH1F *fOmegaMass;
  TH1F *fDcaXiDaughters;
  TH1F *fDcaXiToPrimVertex;
  TH1F *fDcaBacToPrimVertex;
  TH1F *fCosPointingAngleXi;
  TH1F *fDecayLengthXi;
  TH1F *fEtaXi;
  TH1F *fPtXi;
  TH1F *fPtBac;
  TH2D *fdEdxBac;
  TH2D *fTOFtimeBacPi;
  TH1D *fnsigmaBacPi;
  TH2D *fTOFtimeBacK;
  TH1D *fnsigmaBacK;
};

#endif
