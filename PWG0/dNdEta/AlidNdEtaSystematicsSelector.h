/* $Id$ */

#ifndef ALIDNDETASYSTEMATICSSELECTOR_H
#define ALIDNDETASYSTEMATICSSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class AlidNdEtaCorrection;
class TParticlePDG;

class TH2F;
class TH1F;

class AlidNdEtaSystematicsSelector : public AliSelectorRL {
  public:
    AlidNdEtaSystematicsSelector();
    virtual ~AlidNdEtaSystematicsSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    Bool_t SignOK(TParticlePDG* particle);

    void ReadUserObjects(TTree* tree);

    void FillCorrectionMaps(TObjArray* listOfTracks);
    void FillSecondaries();
    void FillSigmaVertex();

    TH2F* fSecondaries; // (Nprim/Nsec for the cases: all/above3GeV/reconstructed tracks/accepted tracks) vs (particle count)

    AlidNdEtaCorrection* fdNdEtaCorrectionSpecies[4];      // correction for different particle species: here pi, K, p, others

    TH1F* fSigmaVertex; // (accepted tracks) vs (n of sigma to vertex cut)

    AliESDtrackCuts* fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

    TH1F* fPIDParticles; // pid of primary particles
    TH1F* fPIDTracks; // pid of reconstructed tracks

    AlidNdEtaCorrection* fdNdEtaCorrectionVertexReco[3]; // correction for vertex reco eff

    AlidNdEtaCorrection* fdNdEtaCorrectionTriggerBias[3]; // correction for trigger bias

    Int_t fSignMode; // 1 = only positive particles are counted, -1 = only negative, 0 = both (default)
    Int_t fMultiplicityMode; // 1 = only events with low multiplicity, 2 = high multiplicity, 0 = all (default)

 private:
    AlidNdEtaSystematicsSelector(const AlidNdEtaSystematicsSelector&);
    AlidNdEtaSystematicsSelector& operator=(const AlidNdEtaSystematicsSelector&);

  ClassDef(AlidNdEtaSystematicsSelector, 0);
};

#endif
