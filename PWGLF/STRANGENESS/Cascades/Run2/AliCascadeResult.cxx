//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "TList.h"
#include "TH3F.h"
#include "TProfile.h"
#include "AliVWeakResult.h"
#include "AliCascadeResult.h"
#include "AliLog.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliCascadeResult);
//________________________________________________________________
AliCascadeResult::AliCascadeResult() :
AliVWeakResult(),
fMassHypo(AliCascadeResult::kXiMinus),
fhNCentBounds(-1),
fhCentBins(0x0),
fhNPtBounds(-1),
fhPtBins(0x0),
fhNMassBins(-1),
fhMinMass(-1),
fhMaxMass(-1),
fHisto(0x0),
fProtonProfile(0x0),
//Acceptance Cuts
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
//V0 Cuts
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutV0Radius(5.0),
//Cascade Cuts
fCutDCAV0ToPV(0.05),
fCutV0Mass(0.010),
fCutV0MassSigma(1000),
fCutDCABachToPV(0.03),
fCutDCACascDaughters(2.0),
fCutCascCosPA(0.95),
fCutCascRadius(0.4),
fCutDCABachToBaryon(-1),
fCutBachBaryonCosPA(+2),
fCutMinV0Lifetime(-2),
fCutMaxV0Lifetime(1e+6),
//Miscellaneous
fCutProperLifetime(1000),
fCutTPCdEdx(4.0),
fCutXiRejection(0.008),
fCutMCPhysicalPrimary(kTRUE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fCutMCSelectBump(kFALSE),
fSwapBachCharge(kFALSE),
fSwapBaryon(kFALSE),
fSwapV0MesonCharge(kFALSE),
fSwapV0BaryonCharge(kFALSE),
fCutUse276TeVV0CosPA(kFALSE),
fCutUseTOFUnchecked(kFALSE),
fCutUseITSRefitTracks(kFALSE),
fCutUseITSRefitNegative(kFALSE),
fCutUseITSRefitPositive(kFALSE),
fCutUseITSRefitBachelor(kFALSE),
fCutLeastNumberOfClusters(70),
fCutMinEtaTracks(-0.8),
fCutMaxEtaTracks(+0.8),
fCutMaxChi2PerCluster(1e+5),
fCutMinTrackLength(-1),
fCutUseParametricLength(kFALSE),
fCutLeastNumberOfCrossedRows(-1),
fCutMinCrossedRowsOverLength(-1),

fCutUseVariableCascCosPA(kFALSE),
fCutVarCascCosPA_Exp0Const(0),
fCutVarCascCosPA_Exp0Slope(0),
fCutVarCascCosPA_Exp1Const(0),
fCutVarCascCosPA_Exp1Slope(0),
fCutVarCascCosPA_Const(1),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fCutUseVariableBBCosPA(kFALSE),
fCutVarBBCosPA_Exp0Const(0),
fCutVarBBCosPA_Exp0Slope(0),
fCutVarBBCosPA_Exp1Const(0),
fCutVarBBCosPA_Exp1Slope(0),
fCutVarBBCosPA_Const(1),
fCutUseVariableDCACascDau(kFALSE),
fCutVarDCACascDau_Exp0Const(0),
fCutVarDCACascDau_Exp0Slope(0),
fCutVarDCACascDau_Exp1Const(0),
fCutVarDCACascDau_Exp1Slope(0),
fCutVarDCACascDau_Const(1),
fCutDCACascadeToPV(1e+3),
fCutDCANegToPVWeighted(-1),
fCutDCAPosToPVWeighted(-1),
fCutDCABachToPVWeighted(-1),
fCutAtLeastOneTOF(kFALSE),
fCutITSorTOF(kFALSE),
fCutIsCowboy(0),
fCutIsCascadeCowboy(0)
{
    // Dummy Constructor - not to be used!
    fhNCentBounds = 21;
    fhCentBins = new Double_t[fhNCentBounds];
    for(Long_t ibin=0; ibin<fhNCentBounds; ibin++)
        fhCentBins[ibin] = ibin*100./(fhNCentBounds-1);
    
    //momentum binning assignment
    fhNPtBounds = 201;
    fhPtBins = new Double_t[fhNPtBounds];
    for(Long_t ibin=0; ibin<fhNPtBounds; ibin++)
        fhPtBins[ibin] = ibin*100./(fhNPtBounds-1);
    
    //Invariant mass assignment
    fhNMassBins = 400;
    fhMaxMass = GetMass() + 0.1;
    fhMinMass = GetMass() - 0.1;
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(const char * name, AliCascadeResult::EMassHypo lMassHypo, const char * title):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fhNCentBounds(-1),
fhCentBins(0x0),
fhNPtBounds(-1),
fhPtBins(0x0),
fhNMassBins(-1),
fhMinMass(-1),
fhMaxMass(-1),
fHisto(0x0),
fProtonProfile(0x0),
//Acceptance Cuts
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
//V0 Cuts
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutV0Radius(5.0),
//Cascade Cuts
fCutDCAV0ToPV(0.05),
fCutV0Mass(0.010),
fCutV0MassSigma(1000),
fCutDCABachToPV(0.03),
fCutDCACascDaughters(2.0),
fCutCascCosPA(0.95),
fCutCascRadius(0.4),
fCutDCABachToBaryon(-1),
fCutBachBaryonCosPA(+2),
fCutMinV0Lifetime(-2),
fCutMaxV0Lifetime(1e+6),
//Miscellaneous
fCutProperLifetime(1000),
fCutTPCdEdx(4.0),
fCutXiRejection(0.008),
fCutMCPhysicalPrimary(kTRUE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fCutMCSelectBump(kFALSE),
fSwapBachCharge(kFALSE),
fSwapBaryon(kFALSE),
fSwapV0MesonCharge(kFALSE),
fSwapV0BaryonCharge(kFALSE),
fCutUse276TeVV0CosPA(kFALSE),
fCutUseTOFUnchecked(kFALSE),
fCutUseITSRefitTracks(kFALSE),
fCutUseITSRefitNegative(kFALSE),
fCutUseITSRefitPositive(kFALSE),
fCutUseITSRefitBachelor(kFALSE),
fCutLeastNumberOfClusters(70),
fCutMinEtaTracks(-0.8),
fCutMaxEtaTracks(+0.8),
fCutMaxChi2PerCluster(1e+5),
fCutMinTrackLength(-1),
fCutUseParametricLength(kFALSE),
fCutLeastNumberOfCrossedRows(-1),
fCutMinCrossedRowsOverLength(-1),
fCutUseVariableCascCosPA(kFALSE),
fCutVarCascCosPA_Exp0Const(0),
fCutVarCascCosPA_Exp0Slope(0),
fCutVarCascCosPA_Exp1Const(0),
fCutVarCascCosPA_Exp1Slope(0),
fCutVarCascCosPA_Const(1),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fCutUseVariableBBCosPA(kFALSE),
fCutVarBBCosPA_Exp0Const(0),
fCutVarBBCosPA_Exp0Slope(0),
fCutVarBBCosPA_Exp1Const(0),
fCutVarBBCosPA_Exp1Slope(0),
fCutVarBBCosPA_Const(1),
fCutUseVariableDCACascDau(kFALSE),
fCutVarDCACascDau_Exp0Const(0),
fCutVarDCACascDau_Exp0Slope(0),
fCutVarDCACascDau_Exp1Const(0),
fCutVarDCACascDau_Exp1Slope(0),
fCutVarDCACascDau_Const(1),
fCutDCACascadeToPV(1e+3),
fCutDCANegToPVWeighted(-1),
fCutDCAPosToPVWeighted(-1),
fCutDCABachToPVWeighted(-1),
fCutAtLeastOneTOF(kFALSE),
fCutITSorTOF(kFALSE),
fCutIsCowboy(0),
fCutIsCascadeCowboy(0)
{
    // Named constructor
    fhNCentBounds = 21;
    fhCentBins = new Double_t[fhNCentBounds];
    for(Long_t ibin=0; ibin<fhNCentBounds; ibin++)
        fhCentBins[ibin] = ibin*100./(fhNCentBounds-1);
    
    //momentum binning assignment
    fhNPtBounds = 201;
    fhPtBins = new Double_t[fhNPtBounds];
    for(Long_t ibin=0; ibin<fhNPtBounds; ibin++)
        fhPtBins[ibin] = ibin*100./(fhNPtBounds-1);
    
    //Invariant mass assignment
    fhNMassBins = 400;
    fhMaxMass = GetMass() + 0.1;
    fhMinMass = GetMass() - 0.1;
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(const char * name, AliCascadeResult::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fhNCentBounds(-1),
fhCentBins(0x0),
fhNPtBounds(-1),
fhPtBins(0x0),
fhNMassBins(-1),
fhMinMass(-1),
fhMaxMass(-1),
fHisto(0x0),
fProtonProfile(0x0),
//Acceptance Cuts
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
//V0 Cuts
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutV0Radius(5.0),
//Cascade Cuts
fCutDCAV0ToPV(0.05),
fCutV0Mass(0.010),
fCutV0MassSigma(1000),
fCutDCABachToPV(0.03),
fCutDCACascDaughters(2.0),
fCutCascCosPA(0.95),
fCutCascRadius(0.4),
fCutDCABachToBaryon(-1),
fCutBachBaryonCosPA(+2),
fCutMinV0Lifetime(-2),
fCutMaxV0Lifetime(1e+6),
//Miscellaneous
fCutProperLifetime(1000),
fCutTPCdEdx(4.0),
fCutXiRejection(0.008),
fCutMCPhysicalPrimary(kTRUE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fCutMCSelectBump(kFALSE),
fSwapBachCharge(kFALSE),
fSwapBaryon(kFALSE),
fSwapV0MesonCharge(kFALSE),
fSwapV0BaryonCharge(kFALSE),
fCutUse276TeVV0CosPA(kFALSE),
fCutUseTOFUnchecked(kFALSE),
fCutUseITSRefitTracks(kFALSE),
fCutUseITSRefitNegative(kFALSE),
fCutUseITSRefitPositive(kFALSE),
fCutUseITSRefitBachelor(kFALSE),
fCutLeastNumberOfClusters(70),
fCutMinEtaTracks(-0.8),
fCutMaxEtaTracks(+0.8),
fCutMaxChi2PerCluster(1e+5),
fCutMinTrackLength(-1),
fCutUseParametricLength(kFALSE),
fCutLeastNumberOfCrossedRows(-1),
fCutMinCrossedRowsOverLength(-1),
fCutUseVariableCascCosPA(kFALSE),
fCutVarCascCosPA_Exp0Const(0),
fCutVarCascCosPA_Exp0Slope(0),
fCutVarCascCosPA_Exp1Const(0),
fCutVarCascCosPA_Exp1Slope(0),
fCutVarCascCosPA_Const(1),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fCutUseVariableBBCosPA(kFALSE),
fCutVarBBCosPA_Exp0Const(0),
fCutVarBBCosPA_Exp0Slope(0),
fCutVarBBCosPA_Exp1Const(0),
fCutVarBBCosPA_Exp1Slope(0),
fCutVarBBCosPA_Const(1),
fCutUseVariableDCACascDau(kFALSE),
fCutVarDCACascDau_Exp0Const(0),
fCutVarDCACascDau_Exp0Slope(0),
fCutVarDCACascDau_Exp1Const(0),
fCutVarDCACascDau_Exp1Slope(0),
fCutVarDCACascDau_Const(1),
fCutDCACascadeToPV(1e+3),
fCutDCANegToPVWeighted(-1),
fCutDCAPosToPVWeighted(-1),
fCutDCABachToPVWeighted(-1),
fCutAtLeastOneTOF(kFALSE),
fCutITSorTOF(kFALSE),
fCutIsCowboy(0),
fCutIsCascadeCowboy(0)
{
    //centrality binning assignment
    fhNCentBounds = lNCentBins+1;
    fhCentBins = new Double_t[fhNCentBounds];
    for(Long_t ibin=0; ibin<fhNCentBounds; ibin++)
        fhCentBins[ibin] = lCentBins[ibin];
    
    //momentum binning assignment
    fhNPtBounds = lNPtBins+1;
    fhPtBins = new Double_t[fhNPtBounds];
    for(Long_t ibin=0; ibin<fhNPtBounds; ibin++)
        fhPtBins[ibin] = lPtBins[ibin];
    
    //Invariant mass assignment
    fhNMassBins = 400;
    fhMaxMass = GetMass() + 0.1;
    fhMinMass = GetMass() - 0.1;
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(const char * name, AliCascadeResult::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins, Long_t lNMassBins, Double_t lMinMass, Double_t lMaxMass):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fhNCentBounds(-1),
fhCentBins(0x0),
fhNPtBounds(-1),
fhPtBins(0x0),
fhNMassBins(-1),
fhMinMass(-1),
fhMaxMass(-1),
fHisto(0x0),
fProtonProfile(0x0),
//Acceptance Cuts
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
//V0 Cuts
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutV0Radius(5.0),
//Cascade Cuts
fCutDCAV0ToPV(0.05),
fCutV0Mass(0.010),
fCutV0MassSigma(1000),
fCutDCABachToPV(0.03),
fCutDCACascDaughters(2.0),
fCutCascCosPA(0.95),
fCutCascRadius(0.4),
fCutDCABachToBaryon(-1),
fCutBachBaryonCosPA(+2),
fCutMinV0Lifetime(-2),
fCutMaxV0Lifetime(1e+6),
//Miscellaneous
fCutProperLifetime(1000),
fCutTPCdEdx(4.0),
fCutXiRejection(0.008),
fCutMCPhysicalPrimary(kTRUE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fCutMCSelectBump(kFALSE),
fSwapBachCharge(kFALSE),
fSwapBaryon(kFALSE),
fSwapV0MesonCharge(kFALSE),
fSwapV0BaryonCharge(kFALSE),
fCutUse276TeVV0CosPA(kFALSE),
fCutUseTOFUnchecked(kFALSE),
fCutUseITSRefitTracks(kFALSE),
fCutUseITSRefitNegative(kFALSE),
fCutUseITSRefitPositive(kFALSE),
fCutUseITSRefitBachelor(kFALSE),
fCutLeastNumberOfClusters(70),
fCutMinEtaTracks(-0.8),
fCutMaxEtaTracks(+0.8),
fCutMaxChi2PerCluster(1e+5),
fCutMinTrackLength(-1),
fCutUseParametricLength(kFALSE),
fCutLeastNumberOfCrossedRows(-1),
fCutMinCrossedRowsOverLength(-1),
fCutUseVariableCascCosPA(kFALSE),
fCutVarCascCosPA_Exp0Const(0),
fCutVarCascCosPA_Exp0Slope(0),
fCutVarCascCosPA_Exp1Const(0),
fCutVarCascCosPA_Exp1Slope(0),
fCutVarCascCosPA_Const(1),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fCutUseVariableBBCosPA(kFALSE),
fCutVarBBCosPA_Exp0Const(0),
fCutVarBBCosPA_Exp0Slope(0),
fCutVarBBCosPA_Exp1Const(0),
fCutVarBBCosPA_Exp1Slope(0),
fCutVarBBCosPA_Const(1),
fCutUseVariableDCACascDau(kFALSE),
fCutVarDCACascDau_Exp0Const(0),
fCutVarDCACascDau_Exp0Slope(0),
fCutVarDCACascDau_Exp1Const(0),
fCutVarDCACascDau_Exp1Slope(0),
fCutVarDCACascDau_Const(1),
fCutDCACascadeToPV(1e+3),
fCutDCANegToPVWeighted(-1),
fCutDCAPosToPVWeighted(-1),
fCutDCABachToPVWeighted(-1),
fCutAtLeastOneTOF(kFALSE),
fCutITSorTOF(kFALSE),
fCutIsCowboy(0),
fCutIsCascadeCowboy(0)
{
    //centrality binning assignment
    fhNCentBounds = lNCentBins+1;
    fhCentBins = new Double_t[fhNCentBounds];
    for(Long_t ibin=0; ibin<fhNCentBounds; ibin++)
        fhCentBins[ibin] = lCentBins[ibin];
    
    //momentum binning assignment
    fhNPtBounds = lNPtBins+1;
    fhPtBins = new Double_t[fhNPtBounds];
    for(Long_t ibin=0; ibin<fhNPtBounds; ibin++)
        fhPtBins[ibin] = lPtBins[ibin];
    
    //Invariant mass assignment
    fhNMassBins = lNMassBins;
    fhMaxMass = lMaxMass;
    fhMinMass = lMinMass;
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(const AliCascadeResult& lCopyMe, TString lNewName)
: AliVWeakResult(lCopyMe),
fMassHypo(lCopyMe.fMassHypo),
//Binning matters
fhNCentBounds( lCopyMe.fhNCentBounds),
fhNPtBounds( lCopyMe.fhNPtBounds),
fhNMassBins( lCopyMe.fhNMassBins),
fhMinMass(lCopyMe.fhMinMass),
fhMaxMass(lCopyMe.fhMaxMass),
//Acceptance Cuts
fCutMinRapidity(lCopyMe.fCutMinRapidity),
fCutMaxRapidity(lCopyMe.fCutMaxRapidity),
//V0 Cuts
fCutDCANegToPV(lCopyMe.fCutDCANegToPV),
fCutDCAPosToPV(lCopyMe.fCutDCAPosToPV),
fCutDCAV0Daughters(lCopyMe.fCutDCAV0Daughters),
fCutV0CosPA(lCopyMe.fCutV0CosPA),
fCutV0Radius(lCopyMe.fCutV0Radius),
//Cascade Cuts
fCutDCAV0ToPV(lCopyMe.fCutDCAV0ToPV),
fCutV0Mass(lCopyMe.fCutV0Mass),
fCutV0MassSigma(lCopyMe.fCutV0MassSigma),
fCutDCABachToPV(lCopyMe.fCutDCABachToPV),
fCutDCACascDaughters(lCopyMe.fCutDCACascDaughters),
fCutCascCosPA(lCopyMe.fCutCascCosPA),
fCutCascRadius(lCopyMe.fCutCascRadius),
fCutDCABachToBaryon(lCopyMe.fCutDCABachToBaryon),
fCutBachBaryonCosPA(lCopyMe.fCutBachBaryonCosPA),
fCutMinV0Lifetime(lCopyMe.fCutMinV0Lifetime),
fCutMaxV0Lifetime(lCopyMe.fCutMaxV0Lifetime),
//Miscellaneous
fCutProperLifetime(lCopyMe.fCutProperLifetime),
fCutTPCdEdx(lCopyMe.fCutTPCdEdx),
fCutXiRejection(lCopyMe.fCutXiRejection),
//MC specific
fCutMCPhysicalPrimary(lCopyMe.fCutMCPhysicalPrimary),
fCutMCPDGCodeAssociation(lCopyMe.fCutMCPDGCodeAssociation),
fCutMCUseMCProperties(lCopyMe.fCutMCUseMCProperties),
fCutMCSelectBump(lCopyMe.fCutMCSelectBump),
//Rsn-like bg subtraction
fSwapBachCharge(lCopyMe.fSwapBachCharge),
fSwapBaryon(lCopyMe.fSwapBaryon),
fSwapV0MesonCharge(lCopyMe.fSwapV0MesonCharge),
fSwapV0BaryonCharge(lCopyMe.fSwapV0BaryonCharge),
//276 Reanalysis
fCutUse276TeVV0CosPA(lCopyMe.fCutUse276TeVV0CosPA),
fCutUseTOFUnchecked(lCopyMe.fCutUseTOFUnchecked),
//Track selections
fCutUseITSRefitTracks(lCopyMe.fCutUseITSRefitTracks),
fCutUseITSRefitNegative(lCopyMe.fCutUseITSRefitNegative),
fCutUseITSRefitPositive(lCopyMe.fCutUseITSRefitPositive),
fCutUseITSRefitBachelor(lCopyMe.fCutUseITSRefitBachelor),
fCutLeastNumberOfClusters(lCopyMe.fCutLeastNumberOfClusters),
fCutMinEtaTracks(lCopyMe.fCutMinEtaTracks),
fCutMaxEtaTracks(lCopyMe.fCutMaxEtaTracks),
fCutMaxChi2PerCluster(lCopyMe.fCutMaxChi2PerCluster),
fCutMinTrackLength(lCopyMe.fCutMinTrackLength),
fCutUseParametricLength(lCopyMe.fCutUseParametricLength),
fCutLeastNumberOfCrossedRows(lCopyMe.fCutLeastNumberOfCrossedRows),
fCutMinCrossedRowsOverLength(lCopyMe.fCutMinCrossedRowsOverLength),

fCutUseVariableCascCosPA(lCopyMe.fCutUseVariableCascCosPA),
fCutVarCascCosPA_Exp0Const(lCopyMe.fCutVarCascCosPA_Exp0Const),
fCutVarCascCosPA_Exp0Slope(lCopyMe.fCutVarCascCosPA_Exp0Slope),
fCutVarCascCosPA_Exp1Const(lCopyMe.fCutVarCascCosPA_Exp1Const),
fCutVarCascCosPA_Exp1Slope(lCopyMe.fCutVarCascCosPA_Exp1Slope),
fCutVarCascCosPA_Const(lCopyMe.fCutVarCascCosPA_Const),
fCutUseVariableV0CosPA(lCopyMe.fCutUseVariableV0CosPA),
fCutVarV0CosPA_Exp0Const(lCopyMe.fCutVarV0CosPA_Exp0Const),
fCutVarV0CosPA_Exp0Slope(lCopyMe.fCutVarV0CosPA_Exp0Slope),
fCutVarV0CosPA_Exp1Const(lCopyMe.fCutVarV0CosPA_Exp1Const),
fCutVarV0CosPA_Exp1Slope(lCopyMe.fCutVarV0CosPA_Exp1Slope),
fCutVarV0CosPA_Const(lCopyMe.fCutVarV0CosPA_Const),
fCutUseVariableBBCosPA(lCopyMe.fCutUseVariableBBCosPA),
fCutVarBBCosPA_Exp0Const(lCopyMe.fCutVarBBCosPA_Exp0Const),
fCutVarBBCosPA_Exp0Slope(lCopyMe.fCutVarBBCosPA_Exp0Slope),
fCutVarBBCosPA_Exp1Const(lCopyMe.fCutVarBBCosPA_Exp1Const),
fCutVarBBCosPA_Exp1Slope(lCopyMe.fCutVarBBCosPA_Exp1Slope),
fCutVarBBCosPA_Const(lCopyMe.fCutVarBBCosPA_Const),
fCutUseVariableDCACascDau(lCopyMe.fCutUseVariableDCACascDau),
fCutVarDCACascDau_Exp0Const(lCopyMe.fCutVarDCACascDau_Exp0Const),
fCutVarDCACascDau_Exp0Slope(lCopyMe.fCutVarDCACascDau_Exp0Slope),
fCutVarDCACascDau_Exp1Const(lCopyMe.fCutVarDCACascDau_Exp1Const),
fCutVarDCACascDau_Exp1Slope(lCopyMe.fCutVarDCACascDau_Exp1Slope),
fCutVarDCACascDau_Const(lCopyMe.fCutVarDCACascDau_Const),
fCutDCACascadeToPV(lCopyMe.fCutDCACascadeToPV),
fCutDCANegToPVWeighted(lCopyMe.fCutDCANegToPVWeighted),
fCutDCAPosToPVWeighted(lCopyMe.fCutDCAPosToPVWeighted),
fCutDCABachToPVWeighted(lCopyMe.fCutDCABachToPVWeighted),
fCutAtLeastOneTOF(lCopyMe.fCutAtLeastOneTOF),
fCutITSorTOF(lCopyMe.fCutITSorTOF),
fCutIsCowboy(lCopyMe.fCutIsCowboy),
fCutIsCascadeCowboy(lCopyMe.fCutIsCascadeCowboy)

{
    SetName( lNewName.Data() );
    
    // Constructor
    Double_t lThisMass = GetMass();
    
    //Clone objects, if they exist
    if ( lCopyMe.fhCentBins ){
        //centrality binning assignment
        fhCentBins = new Double_t[fhNCentBounds];
        for(Long_t ibin=0; ibin<fhNCentBounds; ibin++)
            fhCentBins[ibin] = lCopyMe.fhCentBins[ibin];
    }
    if ( lCopyMe.fhPtBins){
        //momentum binning assignment
        fhPtBins = new Double_t[fhNPtBounds];
        for(Long_t ibin=0; ibin<fhNPtBounds; ibin++)
            fhPtBins[ibin] = lCopyMe.fhPtBins[ibin];
    }
    if (lCopyMe.GetHistogramToCopy() ){
        //Main output histogram
        fHisto = (TH3F*) lCopyMe.GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    }
    if( lCopyMe.GetProtonProfileToCopy() ){
        //Proton profile (G3/F correction)
        fProtonProfile = (TProfile*) lCopyMe.GetProtonProfileToCopy()->Clone(Form("fProtonProfile_%s",GetName()));
    }
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(AliCascadeResult *lCopyMe, TString lNewName)
: AliVWeakResult(*lCopyMe),
fHisto(0), fProtonProfile(0)
{
    SetName(lNewName.Data());
    fMassHypo = lCopyMe->GetMassHypothesis();
    //Binning matters
    fhNCentBounds = lCopyMe->GetNCentBins()+1;
    fhNPtBounds = lCopyMe->GetNPtBins()+1;
    fhNMassBins = lCopyMe->GetNMassBins();
    fhMinMass = lCopyMe->GetMinMass();
    fhMaxMass = lCopyMe->GetMaxMass();
    //Acceptance Cuts
    fCutMinRapidity     = lCopyMe->GetCutMinRapidity();
    fCutMaxRapidity     = lCopyMe->GetCutMaxRapidity();
    //V0 Cuts
    fCutDCANegToPV     = lCopyMe->GetCutDCANegToPV();
    fCutDCAPosToPV     = lCopyMe->GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe->GetCutDCAV0Daughters();
    fCutV0CosPA        = lCopyMe->GetCutV0CosPA();
    fCutV0Radius       = lCopyMe->GetCutV0Radius();
    //Cascade Cuts
    fCutDCAV0ToPV = lCopyMe -> GetCutDCAV0ToPV();
    fCutV0Mass    = lCopyMe -> GetCutV0Mass();
    fCutV0MassSigma    = lCopyMe -> GetCutV0MassSigma();
    fCutDCABachToPV  = lCopyMe -> GetCutDCABachToPV();
    fCutDCACascDaughters = lCopyMe -> GetCutDCACascDaughters();
    fCutCascCosPA  = lCopyMe -> GetCutCascCosPA();
    fCutCascRadius = lCopyMe -> GetCutCascRadius();
    fCutDCABachToBaryon = lCopyMe -> GetCutDCABachToBaryon();
    fCutBachBaryonCosPA = lCopyMe -> GetCutBachBaryonCosPA();
    fCutMinV0Lifetime = lCopyMe -> GetCutMinV0Lifetime();
    fCutMaxV0Lifetime = lCopyMe -> GetCutMaxV0Lifetime();
    
    //Miscellaneous
    fCutProperLifetime = lCopyMe->GetCutProperLifetime();
    fCutTPCdEdx = lCopyMe->GetCutTPCdEdx();
    fCutXiRejection = lCopyMe->GetCutXiRejection();
    
    //MC specific
    fCutMCPhysicalPrimary    = lCopyMe -> GetCutMCPhysicalPrimary();
    fCutMCPDGCodeAssociation = lCopyMe -> GetCutMCPDGCodeAssociation();
    fCutMCUseMCProperties    = lCopyMe -> GetCutMCUseMCProperties();
    fCutMCSelectBump         = lCopyMe -> GetCutMCSelectBump();
    
    //Rsn-like bg subtraction (experimental)
    fSwapBachCharge = lCopyMe -> GetSwapBachelorCharge();
    fSwapBaryon     = lCopyMe -> GetSwapBaryon();
    fSwapV0MesonCharge     = lCopyMe -> GetSwapV0MesonCharge();
    fSwapV0BaryonCharge     = lCopyMe -> GetSwapV0BaryonCharge();
    
    //2.76 TeV reanalysis
    fCutUse276TeVV0CosPA = lCopyMe ->GetCutUse276TeVV0CosPA();
    
    fCutUseTOFUnchecked = lCopyMe ->GetCutUseTOFUnchecked(),
    
    //Track cuts
    fCutUseITSRefitTracks    = lCopyMe -> GetCutUseITSRefitTracks();
    fCutUseITSRefitNegative  = lCopyMe -> GetCutUseITSRefitNegative();
    fCutUseITSRefitPositive  = lCopyMe -> GetCutUseITSRefitPositive();
    fCutUseITSRefitBachelor  = lCopyMe -> GetCutUseITSRefitBachelor();
    fCutLeastNumberOfClusters= lCopyMe -> GetCutLeastNumberOfClusters();
    fCutMinEtaTracks = lCopyMe -> GetCutMinEtaTracks();
    fCutMaxEtaTracks = lCopyMe -> GetCutMaxEtaTracks();
    fCutMaxChi2PerCluster = lCopyMe -> GetCutMaxChi2PerCluster();
    fCutMinTrackLength = lCopyMe -> GetCutMinTrackLength();
    fCutUseParametricLength = lCopyMe -> GetCutUseParametricLength();
    fCutLeastNumberOfCrossedRows = lCopyMe->GetCutLeastNumberOfCrossedRows();
    fCutMinCrossedRowsOverLength = lCopyMe->GetCutMinCrossedRowsOverLength();
    
    //Variable CascCosPA
    fCutUseVariableCascCosPA = lCopyMe -> GetCutUseVarCascCosPA();
    fCutVarCascCosPA_Exp0Const = lCopyMe -> GetCutVarCascCosPAExp0Const();
    fCutVarCascCosPA_Exp0Slope = lCopyMe -> GetCutVarCascCosPAExp0Slope();
    fCutVarCascCosPA_Exp1Const = lCopyMe -> GetCutVarCascCosPAExp1Const();
    fCutVarCascCosPA_Exp1Slope = lCopyMe -> GetCutVarCascCosPAExp1Slope();
    fCutVarCascCosPA_Const = lCopyMe -> GetCutVarCascCosPAConst();
    
    //Variable V0CosPA
    fCutUseVariableV0CosPA = lCopyMe -> GetCutUseVarV0CosPA();
    fCutVarV0CosPA_Exp0Const = lCopyMe -> GetCutVarV0CosPAExp0Const();
    fCutVarV0CosPA_Exp0Slope = lCopyMe -> GetCutVarV0CosPAExp0Slope();
    fCutVarV0CosPA_Exp1Const = lCopyMe -> GetCutVarV0CosPAExp1Const();
    fCutVarV0CosPA_Exp1Slope = lCopyMe -> GetCutVarV0CosPAExp1Slope();
    fCutVarV0CosPA_Const = lCopyMe -> GetCutVarV0CosPAConst();
    
    //Variable BBCosPA
    fCutUseVariableBBCosPA = lCopyMe -> GetCutUseVarBBCosPA();
    fCutVarBBCosPA_Exp0Const = lCopyMe -> GetCutVarBBCosPAExp0Const();
    fCutVarBBCosPA_Exp0Slope = lCopyMe -> GetCutVarBBCosPAExp0Slope();
    fCutVarBBCosPA_Exp1Const = lCopyMe -> GetCutVarBBCosPAExp1Const();
    fCutVarBBCosPA_Exp1Slope = lCopyMe -> GetCutVarBBCosPAExp1Slope();
    fCutVarBBCosPA_Const = lCopyMe -> GetCutVarBBCosPAConst();
    
    //Variable DCACascDau
    fCutUseVariableDCACascDau = lCopyMe -> GetCutUseVarDCACascDau();
    fCutVarDCACascDau_Exp0Const = lCopyMe -> GetCutVarDCACascDauExp0Const();
    fCutVarDCACascDau_Exp0Slope = lCopyMe -> GetCutVarDCACascDauExp0Slope();
    fCutVarDCACascDau_Exp1Const = lCopyMe -> GetCutVarDCACascDauExp1Const();
    fCutVarDCACascDau_Exp1Slope = lCopyMe -> GetCutVarDCACascDauExp1Slope();
    fCutVarDCACascDau_Const = lCopyMe -> GetCutVarDCACascDauConst();
    
    fCutDCACascadeToPV      = lCopyMe -> GetCutDCACascadeToPV();
    fCutDCANegToPVWeighted  = lCopyMe -> GetCutDCANegToPVWeighted();
    fCutDCAPosToPVWeighted  = lCopyMe -> GetCutDCAPosToPVWeighted();
    fCutDCABachToPVWeighted = lCopyMe -> GetCutDCABachToPVWeighted();
    
    fCutAtLeastOneTOF = lCopyMe -> GetCutAtLeastOneTOF();
    fCutITSorTOF = lCopyMe -> GetCutITSorTOF(); 
    
    fCutIsCowboy = lCopyMe -> GetCutIsCowboy();
    fCutIsCascadeCowboy = lCopyMe -> GetCutIsCascadeCowboy();
    
    // Constructor
    Double_t lThisMass = GetMass();
    
    //Clone objects, if they exist
    if ( lCopyMe->GetCentBins() ){
        //centrality binning assignment
        fhCentBins = new Double_t[fhNCentBounds];
        for(Long_t ibin=0; ibin<fhNCentBounds; ibin++)
            fhCentBins[ibin] = (lCopyMe->GetCentBins())[ibin];
    }
    if ( lCopyMe->GetPtBins()){
        //momentum binning assignment
        fhPtBins = new Double_t[fhNPtBounds];
        for(Long_t ibin=0; ibin<fhNPtBounds; ibin++)
            fhPtBins[ibin] = (lCopyMe->GetPtBins())[ibin];
    }
    if (lCopyMe->GetHistogramToCopy() ){
        //Main output histogram
        fHisto = (TH3F*) lCopyMe->GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    }
    if( lCopyMe->GetProtonProfileToCopy() ){
        fProtonProfile = (TProfile*) lCopyMe->GetProtonProfileToCopy()->Clone(Form("fProtonProfile_%s",GetName()));
    }
}
//________________________________________________________________
AliCascadeResult::~AliCascadeResult(){
    // Proper destructor: delete all pointer data members
    if (fhCentBins) {
        delete [] fhCentBins;
    }
    if (fhPtBins) {
        delete [] fhPtBins;
    }
    if (fHisto) {
        delete fHisto;
        fHisto = 0x0;
    }
    if (fProtonProfile) {
        delete fProtonProfile;
        fProtonProfile = 0x0;
    }
}

//________________________________________________________________
AliCascadeResult& AliCascadeResult::operator=(const AliCascadeResult& lCopyMe)
{
    if (&lCopyMe == this) return *this;
    SetName(lCopyMe.GetName());
    SetTitle(lCopyMe.GetTitle());
    
    fMassHypo = lCopyMe.GetMassHypothesis();
    //Binning matters
    fhNCentBounds = lCopyMe.GetNCentBins()+1;
    fhNPtBounds = lCopyMe.GetNPtBins()+1;
    fhNMassBins = lCopyMe.GetNMassBins();
    fhMinMass = lCopyMe.GetMinMass();
    fhMaxMass = lCopyMe.GetMaxMass();
    //Acceptance cuts
    fCutMinRapidity = lCopyMe.GetCutMinRapidity();
    fCutMaxRapidity = lCopyMe.GetCutMaxRapidity();
    //V0 Cuts
    fCutDCANegToPV = lCopyMe.GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe.GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe.GetCutDCAV0Daughters();
    fCutV0CosPA = lCopyMe.GetCutV0CosPA();
    fCutV0Radius = lCopyMe.GetCutV0Radius();
    //Cascade Cuts
    fCutDCAV0ToPV = lCopyMe.GetCutDCAV0ToPV();
    fCutV0Mass    = lCopyMe.GetCutV0Mass();
    fCutV0MassSigma  = lCopyMe.GetCutV0MassSigma();
    fCutDCABachToPV  = lCopyMe.GetCutDCABachToPV();
    fCutDCACascDaughters = lCopyMe.GetCutDCACascDaughters();
    fCutCascCosPA  = lCopyMe.GetCutCascCosPA();
    fCutCascRadius = lCopyMe.GetCutCascRadius();
    fCutDCABachToBaryon = lCopyMe.GetCutDCABachToBaryon();
    fCutBachBaryonCosPA = lCopyMe.GetCutBachBaryonCosPA();
    fCutMinV0Lifetime = lCopyMe.GetCutMinV0Lifetime();
    fCutMaxV0Lifetime = lCopyMe.GetCutMaxV0Lifetime();
    
    //Miscellaneous
    fCutProperLifetime = lCopyMe.GetCutProperLifetime();
    fCutTPCdEdx = lCopyMe.GetCutTPCdEdx();
    fCutXiRejection = lCopyMe.GetCutXiRejection();
    
    //MC specific
    fCutMCPhysicalPrimary = lCopyMe.GetCutMCPhysicalPrimary();
    fCutMCPDGCodeAssociation = lCopyMe.GetCutMCPDGCodeAssociation();
    fCutMCUseMCProperties = lCopyMe.GetCutMCUseMCProperties();
    fCutMCSelectBump      = lCopyMe.GetCutMCSelectBump();
    
    //Rsn-like bg subtraction (experimental)
    fSwapBachCharge = lCopyMe.GetSwapBachelorCharge();
    fSwapBaryon          = lCopyMe.GetSwapBaryon();
    fSwapV0MesonCharge     = lCopyMe.GetSwapV0MesonCharge();
    fSwapV0BaryonCharge     = lCopyMe.GetSwapV0BaryonCharge();
    
    //2.76 TeV reanalysis
    fCutUse276TeVV0CosPA = lCopyMe.GetCutUse276TeVV0CosPA();
    
    fCutUseTOFUnchecked = lCopyMe.GetCutUseTOFUnchecked(),
    
    //Track cuts
    fCutUseITSRefitTracks = lCopyMe.GetCutUseITSRefitTracks();
    fCutUseITSRefitNegative  = lCopyMe.GetCutUseITSRefitNegative();
    fCutUseITSRefitPositive  = lCopyMe.GetCutUseITSRefitPositive();
    fCutUseITSRefitBachelor  = lCopyMe.GetCutUseITSRefitBachelor();
    fCutLeastNumberOfClusters = lCopyMe.GetCutLeastNumberOfClusters();
    fCutMinEtaTracks = lCopyMe.GetCutMinEtaTracks();
    fCutMaxEtaTracks = lCopyMe.GetCutMaxEtaTracks();
    fCutMaxChi2PerCluster = lCopyMe.GetCutMaxChi2PerCluster();
    fCutMinTrackLength = lCopyMe.GetCutMinTrackLength();
    fCutUseParametricLength = lCopyMe.GetCutUseParametricLength();
    fCutLeastNumberOfCrossedRows = lCopyMe.GetCutLeastNumberOfCrossedRows();
    fCutMinCrossedRowsOverLength = lCopyMe.GetCutMinCrossedRowsOverLength();
    
    //Variable CascCosPA
    fCutUseVariableCascCosPA = lCopyMe.GetCutUseVarCascCosPA();
    fCutVarCascCosPA_Exp0Const = lCopyMe.GetCutVarCascCosPAExp0Const();
    fCutVarCascCosPA_Exp0Slope = lCopyMe.GetCutVarCascCosPAExp0Slope();
    fCutVarCascCosPA_Exp1Const = lCopyMe.GetCutVarCascCosPAExp1Const();
    fCutVarCascCosPA_Exp1Slope = lCopyMe.GetCutVarCascCosPAExp1Slope();
    fCutVarCascCosPA_Const = lCopyMe.GetCutVarCascCosPAConst();
    
    //Variable V0CosPA
    fCutUseVariableV0CosPA = lCopyMe.GetCutUseVarV0CosPA();
    fCutVarV0CosPA_Exp0Const = lCopyMe.GetCutVarV0CosPAExp0Const();
    fCutVarV0CosPA_Exp0Slope = lCopyMe.GetCutVarV0CosPAExp0Slope();
    fCutVarV0CosPA_Exp1Const = lCopyMe.GetCutVarV0CosPAExp1Const();
    fCutVarV0CosPA_Exp1Slope = lCopyMe.GetCutVarV0CosPAExp1Slope();
    fCutVarV0CosPA_Const = lCopyMe.GetCutVarV0CosPAConst();
    
    //Variable BBCosPA
    fCutUseVariableBBCosPA = lCopyMe.GetCutUseVarBBCosPA();
    fCutVarBBCosPA_Exp0Const = lCopyMe.GetCutVarBBCosPAExp0Const();
    fCutVarBBCosPA_Exp0Slope = lCopyMe.GetCutVarBBCosPAExp0Slope();
    fCutVarBBCosPA_Exp1Const = lCopyMe.GetCutVarBBCosPAExp1Const();
    fCutVarBBCosPA_Exp1Slope = lCopyMe.GetCutVarBBCosPAExp1Slope();
    fCutVarBBCosPA_Const = lCopyMe.GetCutVarBBCosPAConst();
    
    //Variable DCACascDau
    fCutUseVariableDCACascDau = lCopyMe.GetCutUseVarDCACascDau();
    fCutVarDCACascDau_Exp0Const = lCopyMe.GetCutVarDCACascDauExp0Const();
    fCutVarDCACascDau_Exp0Slope = lCopyMe.GetCutVarDCACascDauExp0Slope();
    fCutVarDCACascDau_Exp1Const = lCopyMe.GetCutVarDCACascDauExp1Const();
    fCutVarDCACascDau_Exp1Slope = lCopyMe.GetCutVarDCACascDauExp1Slope();
    fCutVarDCACascDau_Const = lCopyMe.GetCutVarDCACascDauConst();
    
    fCutDCACascadeToPV      = lCopyMe.GetCutDCACascadeToPV();
    fCutDCANegToPVWeighted  = lCopyMe.GetCutDCANegToPVWeighted();
    fCutDCAPosToPVWeighted  = lCopyMe.GetCutDCAPosToPVWeighted();
    fCutDCABachToPVWeighted = lCopyMe.GetCutDCABachToPVWeighted();
    
    fCutAtLeastOneTOF = lCopyMe.GetCutAtLeastOneTOF();
    fCutITSorTOF = lCopyMe.GetCutITSorTOF(); 
    
    fCutIsCowboy = lCopyMe.GetCutIsCowboy();
    fCutIsCascadeCowboy = lCopyMe.GetCutIsCascadeCowboy();
    
    if (fhCentBins) {
        delete [] fhCentBins;
    }
    if (fhPtBins) {
        delete [] fhPtBins;
    }
    if (fHisto) {
        delete fHisto;
        fHisto = 0;
    }
    if (fProtonProfile) {
        delete fProtonProfile;
        fProtonProfile = 0;
    }
    // Constructor
    Double_t lThisMass = GetMass();
    
    //Clone objects, if they exist
    if ( lCopyMe.GetCentBins() ){
        //centrality binning assignment
        fhCentBins = new Double_t[fhNCentBounds];
        for(Long_t ibin=0; ibin<fhNCentBounds; ibin++)
            fhCentBins[ibin] = (lCopyMe.GetCentBins())[ibin];
    }
    if ( lCopyMe.GetPtBins()){
        //momentum binning assignment
        fhPtBins = new Double_t[fhNPtBounds];
        for(Long_t ibin=0; ibin<fhNPtBounds; ibin++)
            fhPtBins[ibin] = (lCopyMe.GetPtBins())[ibin];
    }
    if (lCopyMe.GetHistogramToCopy() ){
        //Main output histogram
        fHisto = (TH3F*) lCopyMe.GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    }
    if( lCopyMe.GetProtonProfileToCopy() ){
        fProtonProfile = (TProfile*) lCopyMe.GetProtonProfileToCopy()->Clone(Form("fProtonProfile_%s",GetName()));
    }
    
    return *this;
}

//________________________________________________________________
Long64_t AliCascadeResult::Merge(TCollection *hlist)
//Merging function to allow for usage as analysis output
{
    
    if (!hlist) return 0;
    if (hlist->IsEmpty()) return (Long64_t) (GetHistogram() ? GetHistogram()->GetEntries() : 0 );
    
    if (hlist) {
        AliCascadeResult *xh = 0;
        TIter nxh(hlist);
        while ((xh = (AliCascadeResult *) nxh())) {
            // Check if you're not committing a crime
            if( ! HasSameCuts( xh ) ){
                
                AliFatal(Form("FATAL: you're trying to sum output that was obtained with different selections! Offending object: %s",GetName()));
                
            }
            //... if all fine, add this histogram
            if ( GetHistogram() && xh->GetHistogram() )
                GetHistogram()->Add(xh->GetHistogram());
            
            //... if proton profiles are both defined, merge that as well, please
            if ( fProtonProfile && xh->GetProtonProfileToCopy() )
                GetProtonProfile()->Add(xh->GetProtonProfile());
        }
    }
    return (Long64_t) (GetHistogram() ? GetHistogram()->GetEntries() : 0 );
}
//________________________________________________________________
Bool_t AliCascadeResult::HasSameCuts(AliVWeakResult *lCompare, Bool_t lCheckdEdx )
//Function to compare the cuts contained in this result with another
//Returns kTRUE if all selection cuts are identical within 1e-6
//WARNING: Does not check MC association flags 
{
    Bool_t lReturnValue = kTRUE;
    
    if( !lCompare->InheritsFrom(AliCascadeResult::Class() ) ){
        //Apples and oranges! Return kFALSE
        return kFALSE;
    }
    
    AliCascadeResult *lCompareCascade = (AliCascadeResult*) lCompare;
    
    if( fMassHypo != lCompareCascade->GetMassHypothesis() ) lReturnValue = kFALSE;
    
    //Acceptance
    if( TMath::Abs( fCutMinRapidity - lCompareCascade->GetCutMinRapidity() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMaxRapidity - lCompareCascade->GetCutMaxRapidity() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //V0 Selection Criteria
    if( TMath::Abs( fCutDCANegToPV - lCompareCascade->GetCutDCANegToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCAPosToPV - lCompareCascade->GetCutDCAPosToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCAV0Daughters - lCompareCascade->GetCutDCAV0Daughters() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0CosPA - lCompareCascade->GetCutV0CosPA() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0Radius - lCompareCascade->GetCutV0Radius() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //Cascade Selection Criteria
    if( TMath::Abs( fCutDCAV0ToPV - lCompareCascade->GetCutDCAV0ToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0Mass - lCompareCascade->GetCutV0Mass() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0MassSigma - lCompareCascade->GetCutV0MassSigma() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCABachToPV - lCompareCascade->GetCutDCABachToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCACascDaughters - lCompareCascade->GetCutDCACascDaughters() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutCascCosPA - lCompareCascade->GetCutCascCosPA() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutCascRadius - lCompareCascade->GetCutCascRadius() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCABachToBaryon - lCompareCascade->GetCutDCABachToBaryon() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutBachBaryonCosPA - lCompareCascade->GetCutBachBaryonCosPA() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMinV0Lifetime - lCompareCascade->GetCutMinV0Lifetime() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMaxV0Lifetime - lCompareCascade->GetCutMaxV0Lifetime() ) > 1e-6 ) lReturnValue = kFALSE;
    
    if( TMath::Abs( fCutProperLifetime - lCompareCascade->GetCutProperLifetime() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutTPCdEdx - lCompareCascade->GetCutTPCdEdx() ) > 1e-6 && lCheckdEdx ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutXiRejection - lCompareCascade->GetCutXiRejection() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //Track cuts
    if( TMath::Abs( fCutLeastNumberOfClusters - lCompareCascade->GetCutLeastNumberOfClusters() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutUseITSRefitTracks - lCompareCascade->GetCutUseITSRefitTracks() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutUseITSRefitNegative - lCompareCascade->GetCutUseITSRefitNegative() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutUseITSRefitPositive - lCompareCascade->GetCutUseITSRefitPositive() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutUseITSRefitBachelor - lCompareCascade->GetCutUseITSRefitBachelor() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //Check if parametric V0 CosPA (as in 2.76 analysis) used
    if( fCutUse276TeVV0CosPA != lCompareCascade->GetCutUse276TeVV0CosPA() ) lReturnValue = kFALSE;
    
    if( fCutUseTOFUnchecked != lCompareCascade->GetCutUseTOFUnchecked() ) lReturnValue = kFALSE;
    
    if( TMath::Abs( fCutMinEtaTracks - lCompareCascade->GetCutMinEtaTracks() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMaxEtaTracks - lCompareCascade->GetCutMaxEtaTracks() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMaxChi2PerCluster - lCompareCascade->GetCutMaxChi2PerCluster() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMinTrackLength - lCompareCascade->GetCutMinTrackLength() ) > 1e-6 ) lReturnValue = kFALSE;
    
    if( TMath::Abs( fCutUseParametricLength - lCompareCascade->GetCutUseParametricLength() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutLeastNumberOfCrossedRows - lCompareCascade->GetCutLeastNumberOfCrossedRows() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMinCrossedRowsOverLength - lCompareCascade->GetCutMinCrossedRowsOverLength() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //Variable CascCosPA
    if ( TMath::Abs(fCutUseVariableCascCosPA - lCompareCascade->GetCutUseVarCascCosPA()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarCascCosPA_Exp0Const - lCompareCascade->GetCutVarCascCosPAExp0Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarCascCosPA_Exp0Slope - lCompareCascade->GetCutVarCascCosPAExp0Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarCascCosPA_Exp1Const - lCompareCascade->GetCutVarCascCosPAExp1Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarCascCosPA_Exp1Slope - lCompareCascade->GetCutVarCascCosPAExp1Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarCascCosPA_Const  - lCompareCascade->GetCutVarCascCosPAConst()) > 1e-6 ) lReturnValue = kFALSE;
    
    //Variable V0CosPA
    if ( TMath::Abs(fCutUseVariableV0CosPA - lCompareCascade->GetCutUseVarV0CosPA()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp0Const - lCompareCascade->GetCutVarV0CosPAExp0Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp0Slope - lCompareCascade->GetCutVarV0CosPAExp0Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp1Const - lCompareCascade->GetCutVarV0CosPAExp1Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp1Slope - lCompareCascade->GetCutVarV0CosPAExp1Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Const  - lCompareCascade->GetCutVarV0CosPAConst()) > 1e-6 ) lReturnValue = kFALSE;
    
    //Variable BBCosPA
    if ( TMath::Abs(fCutUseVariableBBCosPA - lCompareCascade->GetCutUseVarBBCosPA()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarBBCosPA_Exp0Const - lCompareCascade->GetCutVarBBCosPAExp0Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarBBCosPA_Exp0Slope - lCompareCascade->GetCutVarBBCosPAExp0Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarBBCosPA_Exp1Const - lCompareCascade->GetCutVarBBCosPAExp1Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarBBCosPA_Exp1Slope - lCompareCascade->GetCutVarBBCosPAExp1Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarBBCosPA_Const  - lCompareCascade->GetCutVarBBCosPAConst()) > 1e-6 ) lReturnValue = kFALSE;
    
    //Variable DCACascDau
    if ( TMath::Abs(fCutUseVariableDCACascDau - lCompareCascade->GetCutUseVarDCACascDau()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarDCACascDau_Exp0Const - lCompareCascade->GetCutVarDCACascDauExp0Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarDCACascDau_Exp0Slope - lCompareCascade->GetCutVarDCACascDauExp0Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarDCACascDau_Exp1Const - lCompareCascade->GetCutVarDCACascDauExp1Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarDCACascDau_Exp1Slope - lCompareCascade->GetCutVarDCACascDauExp1Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarDCACascDau_Const  - lCompareCascade->GetCutVarDCACascDauConst()) > 1e-6 ) lReturnValue = kFALSE;
    
    if ( TMath::Abs(fCutDCACascadeToPV       - lCompareCascade->GetCutDCACascadeToPV()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutDCANegToPVWeighted   - lCompareCascade->GetCutDCANegToPVWeighted()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutDCAPosToPVWeighted   - lCompareCascade->GetCutDCAPosToPVWeighted()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutDCABachToPVWeighted  - lCompareCascade->GetCutDCABachToPVWeighted()) > 1e-6 ) lReturnValue = kFALSE;
    
    if ( TMath::Abs(fCutAtLeastOneTOF - lCompareCascade->GetCutAtLeastOneTOF()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutITSorTOF - lCompareCascade->GetCutITSorTOF()) > 1e-6 ) lReturnValue = kFALSE;
    
    if ( TMath::Abs(fCutIsCowboy - lCompareCascade->GetCutIsCowboy()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutIsCascadeCowboy - lCompareCascade->GetCutIsCascadeCowboy()) > 1e-6 ) lReturnValue = kFALSE;
    
    return lReturnValue;
}
//________________________________________________________________
void AliCascadeResult::Print()
//Function to compare the cuts contained in this result with another
//Returns kTRUE if all selection cuts are identical within 1e-6
//WARNING: Does not check MC association flags
{
    cout<<"========================================"<<endl;
    cout<<"    AliCascadeResult Configuration      "<<endl;
    cout<<"========================================"<<endl;
    cout<<" Object Name........: "<<this->GetName()<<endl;
    if( fHisto)
        cout<<" Histogram Name.....: "<<fHisto->GetName()<<endl;
    if( fProtonProfile )
        cout<<" Proton profile.....: will be saved"<<endl;
    if( fMassHypo == AliCascadeResult::kXiMinus      ) cout<<" Mass Hypothesis....: XiMinus"<<endl;
    if( fMassHypo == AliCascadeResult::kXiPlus       ) cout<<" Mass Hypothesis....: XiPlus"<<endl;
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) cout<<" Mass Hypothesis....: OmegaMinus"<<endl;
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) cout<<" Mass Hypothesis....: OmegaPlus"<<endl;
    cout<<" Expected mass......: "<<GetMass()<<endl;
    cout<<" Min y..............: "<<fCutMinRapidity<<endl;
    cout<<" Max y..............: "<<fCutMaxRapidity<<endl;
    
    cout<<" DCA Neg to PV......: "<<fCutDCANegToPV<<endl;
    cout<<" DCA Pos to PV......: "<<fCutDCAPosToPV<<endl;
    cout<<" DCA V0 Daughters...: "<<fCutDCAV0Daughters<<endl;
    cout<<" V0 Cos PA..........: "<<fCutV0CosPA<<endl;
    cout<<" Use Var V0CosPA..: "<<fCutUseVariableV0CosPA<<endl;
    if( fCutUseVariableV0CosPA ){
        cout<<" ^--Exp. 0 Const....: "<<fCutVarV0CosPA_Exp0Const<<endl;
        cout<<" ^--Exp. 0 Slope....: "<<fCutVarV0CosPA_Exp0Slope<<endl;
        cout<<" ^--Exp. 1 Const....: "<<fCutVarV0CosPA_Exp1Const<<endl;
        cout<<" ^--Exp. 1 Slope....: "<<fCutVarV0CosPA_Exp1Slope<<endl;
        cout<<" ^--Constant........: "<<fCutVarV0CosPA_Const<<endl;
    }
    cout<<" V0 2D Radius.......: "<<fCutV0Radius<<endl;
    
    cout<<" DCA V0 to PV.......: "<<fCutDCAV0ToPV<<endl;
    cout<<" V0 Mass............: "<<fCutV0Mass<<endl;
    cout<<" V0 Mass (in sigma).: "<<fCutV0MassSigma<<endl;
    cout<<" DCA Bach to PV.....: "<<fCutDCABachToPV<<endl;
    cout<<" DCA V0 Daughters.: "<<fCutDCACascDaughters<<endl;
    cout<<" Casc Cos PA........: "<<fCutCascCosPA<<endl;
    cout<<" Use Var CascCosPA..: "<<fCutUseVariableCascCosPA<<endl;
    if( fCutUseVariableCascCosPA ){
        cout<<" ^--Exp. 0 Const....: "<<fCutVarCascCosPA_Exp0Const<<endl;
        cout<<" ^--Exp. 0 Slope....: "<<fCutVarCascCosPA_Exp0Slope<<endl;
        cout<<" ^--Exp. 1 Const....: "<<fCutVarCascCosPA_Exp1Const<<endl;
        cout<<" ^--Exp. 1 Slope....: "<<fCutVarCascCosPA_Exp1Slope<<endl;
        cout<<" ^--Constant........: "<<fCutVarCascCosPA_Const<<endl;
    }
    cout<<" Casc 2D Radius.....: "<<fCutCascRadius<<endl;
    cout<<" DCA Bach to Baryon.: "<<fCutDCABachToBaryon<<endl;
    cout<<" Bach Baryon CosPA..: "<<fCutBachBaryonCosPA<<endl;
    cout<<" Use Var BBCosPA....: "<<fCutUseVariableBBCosPA<<endl;
    if( fCutUseVariableBBCosPA ){
        cout<<" ^--Exp. 0 Const....: "<<fCutVarBBCosPA_Exp0Const<<endl;
        cout<<" ^--Exp. 0 Slope....: "<<fCutVarBBCosPA_Exp0Slope<<endl;
        cout<<" ^--Exp. 1 Const....: "<<fCutVarBBCosPA_Exp1Const<<endl;
        cout<<" ^--Exp. 1 Slope....: "<<fCutVarBBCosPA_Exp1Slope<<endl;
        cout<<" ^--Constant........: "<<fCutVarBBCosPA_Const<<endl;
    }
    
    cout<<" Use Var DCACascDau.: "<<fCutUseVariableDCACascDau<<endl;
    if( fCutUseVariableDCACascDau ){
        cout<<" ^--Exp. 0 Const....: "<<fCutVarDCACascDau_Exp0Const<<endl;
        cout<<" ^--Exp. 0 Slope....: "<<fCutVarDCACascDau_Exp0Slope<<endl;
        cout<<" ^--Exp. 1 Const....: "<<fCutVarDCACascDau_Exp1Const<<endl;
        cout<<" ^--Exp. 1 Slope....: "<<fCutVarDCACascDau_Exp1Slope<<endl;
        cout<<" ^--Constant........: "<<fCutVarDCACascDau_Const<<endl;
    }
    cout<<" Min V0 Lifetime....: "<<fCutMinV0Lifetime<<endl;
    cout<<" Max V0 Lifetime....: "<<fCutMaxV0Lifetime<<endl;
    
    cout<<" Proper Lifetime....: "<<fCutProperLifetime<<endl;
    cout<<" TPC dEdx (sigmas)..: "<<fCutTPCdEdx<<endl;
    cout<<" Xi Rej (for Omega).: "<<fCutXiRejection<<endl;
    cout<<" Use 2.76TeV v0cospa: "<<fCutUse276TeVV0CosPA<<endl;
    cout<<" Use ITSref tracks..: "<<fCutUseITSRefitTracks<<endl;
    cout<<" Use ITSref neg.....: "<<fCutUseITSRefitNegative<<endl;
    cout<<" Use ITSref pos.....: "<<fCutUseITSRefitPositive<<endl;
    cout<<" Use ITSref bach.....: "<<fCutUseITSRefitBachelor<<endl;
    cout<<" Nbr Clusters.......: "<<fCutLeastNumberOfClusters<<endl;
    cout<<" Min track eta......: "<<fCutMinEtaTracks<<endl;
    cout<<" Max track eta......: "<<fCutMaxEtaTracks<<endl;
    cout<<" Max chi2/cluster...: "<<fCutMaxChi2PerCluster<<endl;
    cout<<" Min Track Length...: "<<fCutMinTrackLength<<endl;
    cout<<" Param Track Length.: "<<fCutUseParametricLength<<endl;
    cout<<" N(crossed rows)....: "<<fCutLeastNumberOfCrossedRows<<endl;
    cout<<" N(cr)/length.......: "<<fCutMinCrossedRowsOverLength<<endl;
    
    cout<<" MC PDG Association.: "<<fCutMCPDGCodeAssociation<<endl;
    cout<<" MC Phys Primary....: "<<fCutMCPhysicalPrimary<<endl;
    cout<<" MC Use MC pT, y....: "<<fCutMCUseMCProperties<<endl;
    cout<<" MC Select bump.....: "<<fCutMCSelectBump<<endl;
    
    cout<<" Swap bach charge...: "<<fSwapBachCharge<<endl;
    cout<<" Swap lam/lambar....: "<<fSwapBaryon<<endl;
    cout<<" Swap v0 mes charge.: "<<fSwapV0MesonCharge<<endl;
    cout<<" Swap v0 bar charge.: "<<fSwapV0BaryonCharge<<endl;
    
    cout<<" DCA cascade to PV..: "<<fCutDCACascadeToPV<<endl;
    cout<<" wDCA Neg to PV.....: "<<fCutDCANegToPVWeighted<<endl;
    cout<<" wDCA Pos to PV.....: "<<fCutDCAPosToPVWeighted<<endl;
    cout<<" wDCA Bach to PV....: "<<fCutDCABachToPVWeighted<<endl;
    cout<<" At least 1 tof.....: "<<fCutAtLeastOneTOF<<endl;
    cout<<" ITS||TOF...........: "<<fCutITSorTOF<<endl; 
    
    cout<<" Is cowboy..........: "<<fCutIsCowboy<<endl;
    cout<<" Is cascade cowboy..: "<<fCutIsCascadeCowboy<<endl;
    cout<<"========================================"<<endl;
    return;
}

//________________________________________________________________
Double_t AliCascadeResult::GetMass () const
//Get Mass under expected hypothesis 
//N.B. masses are rounded within 1MeV/c^2 just to simplify binning
{
    Double_t lReturnValue = 0;
    //if( fMassHypo == AliCascadeResult::kXiMinus      ) lThisMass = 1.32171;
    //if( fMassHypo == AliCascadeResult::kXiPlus       ) lThisMass = 1.32171;
    //if( fMassHypo == AliCascadeResult::kOmegaMinus   ) lThisMass = 1.67245;
    //if( fMassHypo == AliCascadeResult::kOmegaPlus    ) lThisMass = 1.67245;
    if( fMassHypo == AliCascadeResult::kXiMinus      ) lReturnValue = 1.322;
    if( fMassHypo == AliCascadeResult::kXiPlus       ) lReturnValue = 1.322;
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) lReturnValue = 1.672;
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) lReturnValue = 1.672;
    return lReturnValue;
}

//________________________________________________________________
TString AliCascadeResult::GetParticleName () const
//Get particle name
{
    TString lName = "";
    if( fMassHypo == AliCascadeResult::kXiMinus      ) lName = "XiMinus";
    if( fMassHypo == AliCascadeResult::kXiPlus       ) lName = "XiPlus";
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) lName = "OmegaMinus";
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) lName = "OmegaPlus";
    return lName;
}

//________________________________________________________________
void AliCascadeResult::InitializeHisto ()
//Initialize the main output histogram
{
    if( fHisto ) return; // do nothing
    //AliWarning( Form("Initializing output histogram named %s", GetName() ) ) ;
    //=============================================================
    //Nothing determined
    if ( fhNCentBounds < 0 && fhNPtBounds < 0 && fhNMassBins < 0 ){
        //Main output histogram: Centrality, pt, mass
        Double_t lThisMass = GetMass();
        fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,0,20, 400,lThisMass-0.1,lThisMass+0.1 );
        fHisto->Sumw2();
    }
    //=============================================================
    //centrality, pT determined, but mass not
    if ( fhNCentBounds > 0 && fhNPtBounds > 0 && fhNMassBins < 0 ){
        Double_t lThisMass = GetMass();
        
        //Construct binning in invariant mass as standard: 400 bins from lThisMass-0.1 to lThisMass+1
        Long_t lNMassBins = 400;
        
        Double_t lMassWindow = 0.1 ;
        Double_t lMassDelta = (lMassWindow * 2.) / lNMassBins;
        Double_t lMassBins[lNMassBins+1];
        
        for( Long_t ibound = 0; ibound<lNMassBins+1; ibound++) lMassBins[ibound] = (lThisMass-0.1) + ( ( (Double_t) ibound )*lMassDelta );
        
        //Main output histogram: Centrality, mass, transverse momentum: Variable binning
        fHisto = new TH3F(Form("fHisto_%s",GetName()),"", fhNCentBounds-1, fhCentBins, fhNPtBounds-1, fhPtBins, lNMassBins, lMassBins );
        fHisto->Sumw2();
    }
    //=============================================================
    //Fully custom: centrality + pT + mass
    if ( fhNCentBounds > 0 && fhNPtBounds > 0 && fhNMassBins > 0 ){
        const Long_t lNMassBinsConst = fhNMassBins;
        
        Double_t lMassWindow = (fhMaxMass - fhMinMass)/2.0 ;
        Double_t lMassDelta = (lMassWindow * 2.) / fhNMassBins;
        Double_t lMassBins[lNMassBinsConst+1];
        
        for( Long_t ibound = 0; ibound<fhNMassBins+1; ibound++) lMassBins[ibound] = fhMinMass + ( ( (Double_t) ibound )*lMassDelta );
        /*
         cout<<"Centrality binning: "<<flush;
         for( Long_t ibound = 0; ibound<fhNCentBounds; ibound++) cout<<fhCentBins[ibound]<<" "<<flush;
         cout<<endl;
         cout<<"Pt binning: "<<flush;
         for( Long_t ibound = 0; ibound<fhNPtBounds; ibound++) cout<<fhPtBins[ibound]<<" "<<flush;
         cout<<endl;
         cout<<"Mass binning: "<<flush;
         for( Long_t ibound = 0; ibound<fhNMassBins+1; ibound++) cout<<lMassBins[ibound]<<" "<<flush;
         cout<<endl;
         */
        //Main output histogram: Centrality, mass, transverse momentum: Variable binning
        fHisto = new TH3F(Form("fHisto_%s",GetName()),"", fhNCentBounds-1, fhCentBins, fhNPtBounds-1, fhPtBins, fhNMassBins, lMassBins );
        fHisto->Sumw2();
    }
    //=============================================================
}

//________________________________________________________________
void AliCascadeResult::InitializeProtonProfile () 
//Initialize TProfile to do bookkeeping of proton momenta
{
    if(!fProtonProfile) fProtonProfile = new TProfile( Form("fProtonProfile_%s",GetName()), "", fhNPtBounds-1, fhPtBins);
}

//________________________________________________________________
void AliCascadeResult::InitializeProtonProfile ( Long_t lNPtBins, Double_t *lPtBins ) //kept for compatibility
//Initialize TProfile to do bookkeeping of proton momenta
{
    if(!fProtonProfile) fProtonProfile = new TProfile( Form("fProtonProfile_%s",GetName()), "", lNPtBins, lPtBins);
}



