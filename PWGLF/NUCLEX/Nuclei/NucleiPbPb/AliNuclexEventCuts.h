#ifndef _AliNuclexEventCuts_h_
#define _AliNuclexEventCuts_h_

#include <TNamed.h>
#include <string>
using std::string;

#include "AliVEvent.h"

class TList;

class AliNuclexEventCuts : public TNamed {
  public:
    AliNuclexEventCuts();
 
    bool   AcceptEvent (AliVEvent *ev, TList* qaList = 0x0);
    void AddQAplotsToList(TList *qaList);
    
    void  SetupLHC15o();

    /// While the general philosophy here is to avoid setters and getters (this is not an API we expose)
    /// for some variables (like the max z vertex position) standard the cuts usually follow some patterns
    /// (e.g. the max vertex z position is always symmetric wrt the nominal beam collision point) thus
    /// is convenient having special setters in this case.
    float             GetCentrality (unsigned int estimator = 0) const;
    string            GetCentralityEstimator (unsigned int estimator = 0) const;
    const AliVVertex* GetPrimaryVertex() const { return fPrimaryVertex; }

    void          SetCentralityEstimators (string first = "V0M", string second = "CL0") { fCentEstimators[0] = first; fCentEstimators[1] = second; }
    void          SetCentralityRange (float min, float max) { fMinCentrality = min; fMaxCentrality = max; }
    void          SetMaxVertexZposition (float max) { fMinVtz = -fabs(max); fMaxVtz = fabs(max); }

    bool          fRequireTrackVertex;            ///< if true all the events with only the SPD vertex are rejected
    float         fMinVtz;                        ///< Min z position for the primary vertex
    float         fMaxVtz;                        ///< Max z position for the primary vertex
    float         fMaxDeltaSpdTrackAbsolute;      ///< 
    float         fMaxDeltaSpdTrackNsigmaSPD;     ///<
    float         fMaxDeltaSpdTrackNsigmaTrack;   ///<
//    float         fMaxDispersionSPDvertex;
    float         fMaxResolutionSPDvertex;        ///<

    bool          fRejectDAQincomplete;           ///< Reject events that have incomplete information

    int           fSPDpileupMinContributors;      ///< Reject all the events with SPD pile-up vertices with more than fRejectPileupSPD contributors
    double        fSPDpileupMinZdist;             ///<
		double        fSPDpileupNsigmaZdist;          ///<
		double        fSPDpileupNsigmaDiamXY;         ///<
		double        fSPDpileupNsigmaDiamZ;          ///<

    int           fCentralityFramework;           ///< 0: skip centrality checks, 1: legacy centrality framework, 2: multiplicity framework
    float         fMinCentrality;                 ///< Minimum centrality to be analised
    float         fMaxCentrality;                 ///< Maximum centrality to be analised
    float         fMaxDeltaEstimators;            ///< Maximum difference between two centrality estimators

    bool          fRequireExactTriggerMask;       ///< If true the event selection mask is required to be equal to fTriggerMask
    AliVEvent::EOfflineTriggerTypes fTriggerMask; ///< Trigger mask 

    static const string fgkLabels[2];
  private:
    string        fCentEstimators[2];             ///< Centrality estimators: the first is used as main estimators, that is correlated with the second to monitor spurious events.
    float         fCentPercentiles[2];            ///< Centrality percentiles 
    AliVVertex   *fPrimaryVertex;                 //!<! Primary vertex pointer

    ClassDef(AliNuclexEventCuts,2)
};

#endif

