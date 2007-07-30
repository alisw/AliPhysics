/* $Id$ */

#ifndef ALIMULTIPLICITYMCSELECTOR_H
#define ALIMULTIPLICITYMCSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class AliMultiplicityCorrection;
class AliCorrection;
class TNtuple;
class TH1;

class AliMultiplicityMCSelector : public AliSelectorRL {
  public:
    AliMultiplicityMCSelector();
    virtual ~AliMultiplicityMCSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    AliMultiplicityCorrection* fMultiplicity; // object containing the extracted data
    AliESDtrackCuts* fEsdTrackCuts;           // Object containing the parameters of the esd track cuts

    Bool_t fSystSkipParticles;     // if true skips particles (systematic study)
    AliCorrection* fParticleCorrection[4]; // correction from measured to generated particles for trigger, vertex sample in |eta| < 2;
                                           // for each of the species: pi, k, p, other; for systematic study of pt cut off
    Int_t fSelectProcessType;        // 0 = all (default), 1 = ND, 2 = SD, 3 = DD (for systematic study)
    TNtuple *fParticleSpecies;       // per event: vtx_mc, (pi, k, p, rest (in |eta| < 2)) X (true, recon) + (nolabel,
                                     // doubleTracks, doublePrimaries) [doubleTracks + doublePrimaries are already part of
                                     // rec. particles!)

    TH1* fPtSpectrum;                // function that modifies the pt spectrum (syst. study)

 private:
    AliMultiplicityMCSelector(const AliMultiplicityMCSelector&);
    AliMultiplicityMCSelector& operator=(const AliMultiplicityMCSelector&);

  ClassDef(AliMultiplicityMCSelector, 0);
};

#endif
