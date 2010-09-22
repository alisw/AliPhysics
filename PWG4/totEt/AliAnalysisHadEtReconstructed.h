//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for ESD analysis
//  - reconstruction output
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________
#ifndef ALIANALYSISHADETRECONSTRUCTED_H
#define ALIANALYSISHADETRECONSTRUCTED_H

#include "AliAnalysisHadEt.h"

class AliVParticle;

class AliAnalysisHadEtReconstructed : public AliAnalysisHadEt
{

public:
   
    AliAnalysisHadEtReconstructed();
    virtual ~AliAnalysisHadEtReconstructed();
   
    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();
    
protected:

    bool CheckGoodVertex(AliVParticle *track);
    //virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField) = 0;

    ClassDef(AliAnalysisHadEtReconstructed, 1);
};

#endif // ALIANALYSISHADETRECONSTRUCTED_H
