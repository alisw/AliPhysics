#ifndef ALIANALYSISETRECONSTRUCTED_H
#define ALIANALYSISETRECONSTRUCTED_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis
//  - reconstruction output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEt.h"

class AliVParticle;
class AliESDEvent;

class AliAnalysisEtReconstructed : public AliAnalysisEt
{

public:
   
    AliAnalysisEtReconstructed();
    virtual ~AliAnalysisEtReconstructed();

    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();

    /** Fill the objects you want to output, classes which add new histograms should overload this. */
    virtual void FillOutputList(TList *list);

    /** Create the histograms, must be overloaded if you want to add your own */
    virtual void CreateHistograms();
    
protected:

    bool CheckGoodVertex(AliVParticle *track);
    virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);

    Double_t fTrackDistanceCut; // cut on track distance    
    Double_t fPidCut; // cut on the pid probability
    
    Char_t fClusterType; // selection on cluster type
        
    TH2F *fHistChargedPionEnergyDeposit; /** Energy deposited in calorimeter by charged pions */    
    TH2F *fHistProtonEnergyDeposit; /** Energy deposited in calorimeter by protons */    
    TH2F *fHistAntiProtonEnergyDeposit; /** Energy deposited in calorimeter by anti-protons */    
    TH2F *fHistChargedKaonEnergyDeposit; /** Energy deposited in calorimeter by charged kaons */    
    TH2F *fHistMuonEnergyDeposit; /** Energy deposited in calorimeter by muons */
    
  private:
    
    AliAnalysisEtReconstructed(const AliAnalysisEtReconstructed& g);
    AliAnalysisEtReconstructed & operator=(const AliAnalysisEtReconstructed&);
    
    Double_t CalcTrackClusterDistance(const Float_t pos[3],Int_t *trkMatchId, const AliESDEvent *event);

    ClassDef(AliAnalysisEtReconstructed, 1);
};

#endif // ALIANALYSISETRECONSTRUCTED_H
