//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
#ifndef ALIANALYSISHADETRECONSTRUCTED_H
#define ALIANALYSISHADETRECONSTRUCTED_H

#include "AliAnalysisHadEt.h"

class AliVParticle;

class AliAnalysisHadEtReconstructed : public AliAnalysisHadEt
{

public:
   
    AliAnalysisHadEtReconstructed();
   
    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();
    
protected:

    bool CheckGoodVertex(AliVParticle *track);
    //virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField) = 0;

    Int_t fNTpcClustersCut;
    Int_t fNItsClustersCut;
   
    Double_t fTrackDistanceCut;
    
};

#endif // ALIANALYSISHADETRECONSTRUCTED_H
