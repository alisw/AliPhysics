#ifndef _AliNuclexEventCuts_h_
#define _AliNuclexEventCuts_h_

#include <TNamed.h>
#include <string>
using std::string;

#include "AliVEvent.h"

class AliNuclexEventCuts : public TNamed {
  public:
    AliNuclexEventCuts();
    virtual ~AliNuclexEventCuts();
    AliNuclexEventCuts(const AliNuclexEventCuts& ev);
    AliNuclexEventCuts& operator=(const AliNuclexEventCuts& from);
 
    bool  AcceptEvent (AliVEvent *ev);
    float GetCentrality (unsigned int estimator = 0) { return fCentPercentiles[estimator]; }
    
    void  SetupLHC15o();

    bool          fRequireTrackVertex;
    float         fMinVtz;
    float         fMaxVtz;
    float         fMaxDeltaSpdTrackAbsolute;
    float         fMaxDeltaSpdTrackNsigmaSPD;
    float         fMaxDeltaSpdTrackNsigmaTrack;

    bool          fRejectDAQincomplete;         ///< Reject events that have incomplete information

    int           fRejectPileupSPD;             ///< Reject all the events with SPD pile-up vertices with more than fRejectPileupSPD contributors

    int           fCentralityFramework;         ///< 0: skip centrality checks, 1: legacy centrality framework, 2: multiplicity framework
    float         fMinCentrality;               ///< Minimum centrality to be analised
    float         fMaxCentrality;               ///< Maximum centrality to be analised
    string        fCentEstimators[2];           ///< Centrality estimators: the first is used as main estimators, that is correlated with the second to monitor spurious events.
    float         fCentPercentiles[2];          ///< Centrality percentiles 
    float         fMaxDeltaEstimators;          ///< Maximum difference between two centrality estimators

    AliVEvent::EOfflineTriggerTypes fTriggerMask; ///< Trigger mask 

    TH1I*         fCutStats;                    //-><! Event selection statistics

    /// Monitoring plots: first element of the array is before the selection, the second element is after.
    TH1D*         fVtz[2];                      //-><! Vertex z distribution
    TH1D*         fDeltaTrackSPDvtz[2];         //-><! Difference between the vertex computed using SPD and the track vertex

    TH1D*         fCentrality[2];               //-><! Centrality percentile distribution
    TH2D*         fEstimCorrelation[2];         //-><! Correlation between centrality estimators
    TH2D*         fMultCentCorrelation[2];      //-><! Correlation between main centrality estimator and multiplicity
  private:

    ClassDef(AliNuclexEventCuts,1)
};

#endif

