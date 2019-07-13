/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
//______________________________________________________________________________
// Analysis task for providing various dijet informations
// author: O. Saarimaki, D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

// used in local and grid execution

#ifndef ALIJCDIJETANA_H
#define ALIJCDIJETANA_H

#include <TH1D.h>
#include <TClonesArray.h>
#include <AliJBaseTrack.h>
#include "AliJCDijetHistos.h"

// Fastjet includes
#include <FJ_includes.h>

class AliJCDijetHistos;

class AliJCDijetAna : public TObject
{
    public:
        AliJCDijetAna(); // Default contructor
        virtual ~AliJCDijetAna(); // Destructor
        AliJCDijetAna(const AliJCDijetAna& obj); // Copy constructor
        AliJCDijetAna& operator=(const AliJCDijetAna& obj); // Equal sign operator

        static void CalculateJetsDijets(TClonesArray *inList,
                                        int    lDebug,
                                        int    lCBin,
                                        double lParticleEtaCut,
                                        double lParticlePtCut,
                                        double lJetCone,
                                        double lktJetCone,
                                        int    lktScheme,
                                        bool   lusePionMassInkt,
                                        bool   luseDeltaPhiBGSubtr,
                                        double lConstituentCut,
                                        double lLeadingJetCut,
                                        double lSubleadingJetCut,
                                        double lDeltaPhiCut,
                                        AliJCDijetHistos *fhistos);

    private:

        ClassDef(AliJCDijetAna, 1); // ClassDef needed if inheriting from TObject

};

#endif
