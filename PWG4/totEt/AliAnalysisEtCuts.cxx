//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD & MC analysis
//  - cuts for reconstruction and MonteCarlo 
// implementation file
//  
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtCuts.h"

#include <iostream>
#include "TNamed.h"

using namespace std;

ClassImp(AliAnalysisEtCuts);


AliAnalysisEtCuts::AliAnalysisEtCuts() : 
  TNamed()
				    //
  ,fCommonEtaCut(0.8)
  ,fCommonClusterEnergyCut(0.0)
  ,fCommonTrackPtCut(0.0)
  ,fCommonSingleCell(1)
				    //
  ,fGeometryPhosEtaAccCut(0.12)
  ,fGeometryPhosPhiAccMinCut(260.0)
  ,fGeometryPhosPhiAccMaxCut(320.0)
  ,fGeometryPhosDetectorRadius(460.0)
				    //
  ,fGeometryEmcalEtaAccCut(0.7)
  ,fGeometryEmcalPhiAccMinCut(80.0) // rad 1.4
  ,fGeometryEmcalPhiAccMaxCut(120.0) // rad 2.1
  ,fGeometryEmcalDetectorRadius(440.0)
				    //
  ,fReconstructedVertexXCut(0.5)
  ,fReconstructedVertexYCut(0.5)
  ,fReconstructedVertexZCut(12.0)
  ,fReconstructedIPxyCut(1.5)
  ,fReconstructedIPzCut(1.5)
  ,fReconstructedNTpcClustersCut(30)
  ,fReconstructedNItsClustersCut(3)
				    //
  ,fReconstructedPhosClusterType(-1)
  ,fReconstructedPhosClusterEnergyCut(0.0)
  ,fReconstructedPhosSingleCellEnergyCut(0.5)
  ,fReconstructedPhosTrackDistanceCut(15.0)
				    //
  ,fReconstructedEmcalClusterType(1)
  ,fReconstructedEmcalClusterEnergyCut(0.1) // GeV
  ,fReconstructedEmcalSingleCellEnergyCut(0.5)
  ,fReconstructedEmcalTrackDistanceCut(15.0)
  
  ,fMonteCarloSingleChargedParticle(3)
  ,fMonteCarloNeutralParticle(0)
{ // ctor
}

AliAnalysisEtCuts::~AliAnalysisEtCuts()
{ // dtor
}


