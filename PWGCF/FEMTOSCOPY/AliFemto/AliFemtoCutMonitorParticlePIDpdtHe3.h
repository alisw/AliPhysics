////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoCutMonitorParticlePIDpdtHe3 - the cut monitor for particles to study   ///
/// various aspects of the PID determination                                 ///
/// dowang 2023.10.25                                                                         ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoCutMonitorParticlePIDpdtHe3_hh
#define AliFemtoCutMonitorParticlePIDpdtHe3_hh

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

class AliFemtoCutMonitorParticlePIDpdtHe3 : public AliFemtoCutMonitor{

    public:
        AliFemtoCutMonitorParticlePIDpdtHe3();
        AliFemtoCutMonitorParticlePIDpdtHe3(const char *aName, Int_t aTOFParticle, int MassBin, float LowMass, float UpMass);
	AliFemtoCutMonitorParticlePIDpdtHe3(const AliFemtoCutMonitorParticlePIDpdtHe3 &aCut);
        virtual ~AliFemtoCutMonitorParticlePIDpdtHe3();

        AliFemtoCutMonitorParticlePIDpdtHe3& operator=(const AliFemtoCutMonitorParticlePIDpdtHe3& aCut);

        virtual AliFemtoString Report();
        virtual void Fill(const AliFemtoEvent* aEvent) {AliFemtoCutMonitor::Fill(aEvent);}
        virtual void Fill(const AliFemtoTrack* aTrack);
        virtual void Fill(const AliFemtoV0* aV0) {AliFemtoCutMonitor::Fill(aV0);}
        virtual void Fill(const AliFemtoXi* aXi) {AliFemtoCutMonitor::Fill(aXi);}
        virtual void Fill(const AliFemtoKink* aKink) {AliFemtoCutMonitor::Fill(aKink);}
        virtual void Fill(const AliFemtoPair* aPair) {AliFemtoCutMonitor::Fill(aPair);}
        virtual void Fill(const AliFemtoParticleCollection* aCollection) {AliFemtoCutMonitor::Fill(aCollection);}
        virtual void Fill(const AliFemtoEvent* aEvent,const AliFemtoParticleCollection* aCollection)
        {AliFemtoCutMonitor::Fill(aEvent, aCollection);}
        virtual void Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2) {AliFemtoCutMonitor::Fill(aCollection1, aCollection2);}
        /* void SetTOFParticle(Int_t ipart); */

        void Write();
        void SetUsePt(Bool_t usept){fIfUsePt=usept;}
        virtual TList *GetOutputList();

    protected:
        Int_t fTOFParticle; ///< Select TOF time hypothesis, 0-pion, 1-kaon, 2-proton
        Bool_t fIfUsePt;    ///< Plot pT instead of p in all momentum histograms

        TH2D *fTPCdEdx;        ///< TPC dEdx information
        TH2D *fTOFNSigma;      ///< TOF NSigma values vs mom
        TH2D *fTPCNSigma;      ///< TPC NSigma values vs mom
        TH2D *fMass;

};

#endif
