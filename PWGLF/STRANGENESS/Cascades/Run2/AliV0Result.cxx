//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram 
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "TList.h"
#include "TH3F.h"
#include "AliVWeakResult.h"
#include "AliV0Result.h"
#include "AliLog.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliV0Result);
//________________________________________________________________
AliV0Result::AliV0Result() :
AliVWeakResult(),
fMassHypo(AliV0Result::kK0Short),
fhNPtBoundsFeeddown(-1),
fhPtBinsFeeddown(0x0),
fProtonProfile(0x0),
fHisto(0x0),
fHistoFeeddown(0x0),
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
fCutV0Radius(5.0),
fCutMaxV0Radius(200.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
fCutArmenterosParameter(0.2),
fCutTPCdEdx(3.0),
fCutMinBaryonMomentum(-1),
fCutMCPhysicalPrimary(kTRUE),
fCutMCLambdaFromPrimaryXi(kFALSE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fCutUseITSRefitTracks(kFALSE),
fCutLeastNumberOfCrossedRows(70),
fCutLeastNumberOfCrossedRowsOverFindable(0.8),
fCutMinEtaTracks(-0.8),
fCutMaxEtaTracks(+0.8),
fCutMaxChi2PerCluster(1e+5),
fCutMinTrackLength(-1),
fCutUseParametricLength(kFALSE),
fCutMinCrossedRowsOverLength(-1),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fUseOnTheFly(kFALSE),
fCut276TeVLikedEdx(kFALSE),
fCutAtLeastOneTOF(kFALSE),
fCutITSorTOF(kFALSE), 
fCutIsCowboy(0)
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
AliV0Result::AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fhNPtBoundsFeeddown(-1),
fhPtBinsFeeddown(0x0),
fProtonProfile(0x0),
fHisto(0x0),
fHistoFeeddown(0x0),
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
fCutV0Radius(5.0),
fCutMaxV0Radius(200.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
fCutArmenterosParameter(0.2),
fCutTPCdEdx(3.0),
fCutMinBaryonMomentum(-1),
fCutMCPhysicalPrimary(kTRUE),
fCutMCLambdaFromPrimaryXi(kFALSE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fCutUseITSRefitTracks(kFALSE),
fCutLeastNumberOfCrossedRows(70),
fCutLeastNumberOfCrossedRowsOverFindable(0.8),
fCutMinEtaTracks(-0.8),
fCutMaxEtaTracks(+0.8),
fCutMaxChi2PerCluster(1e+5),
fCutMinTrackLength(-1),
fCutUseParametricLength(kFALSE),
fCutMinCrossedRowsOverLength(-1),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fUseOnTheFly(kFALSE),
fCut276TeVLikedEdx(kFALSE),
fCutAtLeastOneTOF(kFALSE),
fCutITSorTOF(kFALSE), 
fCutIsCowboy(0)
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
AliV0Result::AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fhNPtBoundsFeeddown(-1),
fhPtBinsFeeddown(0x0),
fProtonProfile(0x0),
fHisto(0x0),
fHistoFeeddown(0x0),
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
fCutV0Radius(5.0),
fCutMaxV0Radius(200.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
fCutArmenterosParameter(0.2),
fCutTPCdEdx(3.0),
fCutMinBaryonMomentum(-1),
fCutMCPhysicalPrimary(kTRUE),
fCutMCLambdaFromPrimaryXi(kFALSE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fCutUseITSRefitTracks(kFALSE),
fCutLeastNumberOfCrossedRows(70),
fCutLeastNumberOfCrossedRowsOverFindable(0.8),
fCutMinEtaTracks(-0.8),
fCutMaxEtaTracks(+0.8),
fCutMaxChi2PerCluster(1e+5),
fCutMinTrackLength(-1),
fCutUseParametricLength(kFALSE),
fCutMinCrossedRowsOverLength(-1),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fUseOnTheFly(kFALSE),
fCut276TeVLikedEdx(kFALSE),
fCutAtLeastOneTOF(kFALSE),
fCutITSorTOF(kFALSE), 
fCutIsCowboy(0)
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
AliV0Result::AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins, Long_t lNMassBins, Double_t lMinMass, Double_t lMaxMass):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fhNPtBoundsFeeddown(-1),
fhPtBinsFeeddown(0x0),
fProtonProfile(0x0),
fHisto(0x0),
fHistoFeeddown(0x0),
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
fCutV0Radius(5.0),
fCutMaxV0Radius(200.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
fCutArmenterosParameter(0.2),
fCutTPCdEdx(3.0),
fCutMinBaryonMomentum(-1),
fCutMCPhysicalPrimary(kTRUE),
fCutMCLambdaFromPrimaryXi(kFALSE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fCutUseITSRefitTracks(kFALSE),
fCutLeastNumberOfCrossedRows(70),
fCutLeastNumberOfCrossedRowsOverFindable(0.8),
fCutMinEtaTracks(-0.8),
fCutMaxEtaTracks(+0.8),
fCutMaxChi2PerCluster(1e+5),
fCutMinTrackLength(-1),
fCutUseParametricLength(kFALSE),
fCutMinCrossedRowsOverLength(-1),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fUseOnTheFly(kFALSE),
fCut276TeVLikedEdx(kFALSE),
fCutAtLeastOneTOF(kFALSE),
fCutITSorTOF(kFALSE), 
fCutIsCowboy(0)
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
AliV0Result::AliV0Result(const AliV0Result& lCopyMe, TString lNewName)
: AliVWeakResult(lCopyMe),
fMassHypo(lCopyMe.fMassHypo),
//Binning matters
fhNCentBounds( lCopyMe.fhNCentBounds),
fhNPtBounds( lCopyMe.fhNPtBounds),
fhNPtBoundsFeeddown( lCopyMe.fhNPtBoundsFeeddown),
fhNMassBins( lCopyMe.fhNMassBins),
fhMinMass(lCopyMe.fhMinMass),
fhMaxMass(lCopyMe.fhMaxMass),
//Acceptance Cuts
fCutMinRapidity(lCopyMe.fCutMinRapidity),
fCutMaxRapidity(lCopyMe.fCutMaxRapidity),
//Topological
fCutV0Radius(lCopyMe.fCutV0Radius),
fCutMaxV0Radius(lCopyMe.fCutMaxV0Radius),
fCutDCANegToPV(lCopyMe.fCutDCANegToPV),
fCutDCAPosToPV(lCopyMe.fCutDCAPosToPV),
fCutDCAV0Daughters(lCopyMe.fCutDCAV0Daughters),
fCutV0CosPA(lCopyMe.fCutV0CosPA),
fCutProperLifetime(lCopyMe.fCutProperLifetime),
fCutCompetingV0Rejection(lCopyMe.fCutCompetingV0Rejection),
fCutArmenteros(lCopyMe.fCutArmenteros),
fCutArmenterosParameter(lCopyMe.fCutArmenterosParameter),
fCutTPCdEdx(lCopyMe.fCutTPCdEdx),
fCutMCPhysicalPrimary(lCopyMe.fCutMCPhysicalPrimary),
fCutMCLambdaFromPrimaryXi(lCopyMe.fCutMCLambdaFromPrimaryXi),
fCutMCPDGCodeAssociation(lCopyMe.fCutMCPDGCodeAssociation),
fCutMCUseMCProperties(lCopyMe.fCutMCUseMCProperties),

fCutUseITSRefitTracks(lCopyMe.fCutUseITSRefitTracks),
fCutLeastNumberOfCrossedRows(lCopyMe.fCutLeastNumberOfCrossedRows),
fCutLeastNumberOfCrossedRowsOverFindable(lCopyMe.fCutLeastNumberOfCrossedRowsOverFindable),
fCutMinEtaTracks(lCopyMe.fCutMinEtaTracks),
fCutMaxEtaTracks(lCopyMe.fCutMaxEtaTracks),
fCutMaxChi2PerCluster(lCopyMe.fCutMaxChi2PerCluster),
fCutMinTrackLength(lCopyMe.fCutMinTrackLength),
fCutUseParametricLength(lCopyMe.fCutUseParametricLength),
fCutMinCrossedRowsOverLength(lCopyMe.fCutMinCrossedRowsOverLength),

fCutUseVariableV0CosPA(lCopyMe.fCutUseVariableV0CosPA),
fCutVarV0CosPA_Exp0Const(lCopyMe.fCutVarV0CosPA_Exp0Const),
fCutVarV0CosPA_Exp0Slope(lCopyMe.fCutVarV0CosPA_Exp0Slope),
fCutVarV0CosPA_Exp1Const(lCopyMe.fCutVarV0CosPA_Exp1Const),
fCutVarV0CosPA_Exp1Slope(lCopyMe.fCutVarV0CosPA_Exp1Slope),
fCutVarV0CosPA_Const(lCopyMe.fCutVarV0CosPA_Const),
fUseOnTheFly(lCopyMe.fUseOnTheFly),
fCut276TeVLikedEdx(lCopyMe.fCut276TeVLikedEdx),
fCutAtLeastOneTOF(lCopyMe.fCutAtLeastOneTOF),
fCutITSorTOF(lCopyMe.fCutITSorTOF), 
fCutIsCowboy(lCopyMe.fCutIsCowboy)
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
    if ( lCopyMe.fhPtBinsFeeddown){
        //momentum binning assignment
        fhPtBinsFeeddown = new Double_t[fhNPtBoundsFeeddown];
        for(Long_t ibin=0; ibin<fhNPtBoundsFeeddown; ibin++)
            fhPtBinsFeeddown[ibin] = lCopyMe.fhPtBinsFeeddown[ibin];
    }
    fHisto=0x0;
    if (lCopyMe.GetHistogramToCopy() ){
        //Main output histogram
        fHisto = (TH3F*) lCopyMe.GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    }
    //Copy feeddown matrix, if it exists
    fHistoFeeddown = 0x0;
    if( lCopyMe.GetHistogramFeeddownToCopy() )
        fHistoFeeddown = (TH3F*) lCopyMe.GetHistogramFeeddownToCopy()->Clone(Form("fHistoFeeddown_%s",GetName()));
    fProtonProfile = 0x0;
    //Copy proton profile, if it exists
    if( lCopyMe.GetProtonProfileToCopy() ){
        fProtonProfile = (TProfile*) lCopyMe.GetProtonProfileToCopy()->Clone(Form("fProtonProfile_%s",GetName()));
    }
}
//________________________________________________________________
AliV0Result::AliV0Result(AliV0Result *lCopyMe, TString lNewName)
: AliVWeakResult(*lCopyMe),
fHisto(0)
{
    SetName(lNewName.Data());
    fMassHypo = lCopyMe->GetMassHypothesis();
    
    //Binning matters
    fhNCentBounds = lCopyMe->GetNCentBins()+1;
    fhNPtBounds = lCopyMe->GetNPtBins()+1;
    fhNPtBoundsFeeddown = lCopyMe->GetNPtBinsFeeddown()+1;
    fhNMassBins = lCopyMe->GetNMassBins();
    fhMinMass = lCopyMe->GetMinMass();
    fhMaxMass = lCopyMe->GetMaxMass();
    
    //Acceptance Cuts
    fCutMinRapidity     = lCopyMe->GetCutMinRapidity();
    fCutMaxRapidity     = lCopyMe->GetCutMaxRapidity();
    
    fCutV0Radius = lCopyMe->GetCutV0Radius();
    fCutMaxV0Radius = lCopyMe->GetCutMaxV0Radius();
    fCutDCANegToPV = lCopyMe->GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe->GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe->GetCutDCAV0Daughters();
    fCutV0CosPA = lCopyMe->GetCutV0CosPA();
    fCutProperLifetime = lCopyMe->GetCutProperLifetime();
    fCutCompetingV0Rejection = lCopyMe->GetCutCompetingV0Rejection();
    fCutArmenteros = lCopyMe->GetCutArmenteros();
    fCutArmenterosParameter = lCopyMe->GetCutArmenterosParameter();
    fCutTPCdEdx = lCopyMe->GetCutTPCdEdx();
    fCutMinBaryonMomentum = lCopyMe->GetCutMinBaryonMomentum();
    
    fCutMCPhysicalPrimary = lCopyMe->GetCutMCPhysicalPrimary();
    fCutMCLambdaFromPrimaryXi = lCopyMe->GetCutMCLambdaFromPrimaryXi();
    fCutMCPDGCodeAssociation = lCopyMe->GetCutMCPDGCodeAssociation();
    fCutMCUseMCProperties    = lCopyMe -> GetCutMCUseMCProperties();
    
    //Track cuts
    fCutUseITSRefitTracks    = lCopyMe -> GetCutUseITSRefitTracks();
    fCutLeastNumberOfCrossedRows = lCopyMe->GetCutLeastNumberOfCrossedRows();
    fCutLeastNumberOfCrossedRowsOverFindable = lCopyMe->GetCutLeastNumberOfCrossedRowsOverFindable();
    fCutMinEtaTracks = lCopyMe -> GetCutMinEtaTracks();
    fCutMaxEtaTracks = lCopyMe -> GetCutMaxEtaTracks();
    fCutMaxChi2PerCluster = lCopyMe -> GetCutMaxChi2PerCluster();
    fCutMinTrackLength = lCopyMe -> GetCutMinTrackLength();
    fCutUseParametricLength = lCopyMe -> GetCutUseParametricLength();
    fCutMinCrossedRowsOverLength = lCopyMe->GetCutMinCrossedRowsOverLength();
    
    //Variable V0CosPA
    fCutUseVariableV0CosPA = lCopyMe -> GetCutUseVarV0CosPA();
    fCutVarV0CosPA_Exp0Const = lCopyMe -> GetCutVarV0CosPAExp0Const();
    fCutVarV0CosPA_Exp0Slope = lCopyMe -> GetCutVarV0CosPAExp0Slope();
    fCutVarV0CosPA_Exp1Const = lCopyMe -> GetCutVarV0CosPAExp1Const();
    fCutVarV0CosPA_Exp1Slope = lCopyMe -> GetCutVarV0CosPAExp1Slope();
    fCutVarV0CosPA_Const = lCopyMe -> GetCutVarV0CosPAConst();
    
    //OTF use
    fUseOnTheFly = lCopyMe -> GetUseOnTheFly();
    
    //special dedx
    fCut276TeVLikedEdx = lCopyMe -> GetCut276TeVLikedEdx();
    
    fCutAtLeastOneTOF = lCopyMe -> GetCutAtLeastOneTOF();
    fCutITSorTOF = lCopyMe -> GetCutITSorTOF(); 
    
    fCutIsCowboy = lCopyMe -> GetCutIsCowboy();
    
    // Constructor
    Double_t lThisMass = GetMass();
    
    //Clone objects, if they exist
    if ( lCopyMe->GetCentBins() ){
        //centrality binning assignment
        fhCentBins = new Double_t[fhNCentBounds];
        for(Long_t ibin=0; ibin<fhNCentBounds; ibin++)
            fhCentBins[ibin] = lCopyMe->fhCentBins[ibin];
    }
    if ( lCopyMe->GetPtBins()){
        //momentum binning assignment
        fhPtBins = new Double_t[fhNPtBounds];
        for(Long_t ibin=0; ibin<fhNPtBounds; ibin++)
            fhPtBins[ibin] = lCopyMe->fhPtBins[ibin];
    }
    if ( lCopyMe->GetPtBinsFeeddown()){
        //momentum binning assignment
        fhPtBinsFeeddown = new Double_t[fhNPtBoundsFeeddown];
        for(Long_t ibin=0; ibin<fhNPtBoundsFeeddown; ibin++)
            fhPtBinsFeeddown[ibin] = lCopyMe->fhPtBinsFeeddown[ibin];
    }
    fHisto=0x0;
    if (lCopyMe->GetHistogramToCopy() ){
        //Main output histogram
        fHisto = (TH3F*) lCopyMe->GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    }
    //Copy feeddown matrix, if it exists
    fHistoFeeddown = 0x0;
    if( lCopyMe->GetHistogramFeeddownToCopy() )
        fHistoFeeddown = (TH3F*) lCopyMe->GetHistogramFeeddownToCopy()->Clone(Form("fHistoFeeddown_%s",GetName()));
    fProtonProfile = 0x0;
    //Copy proton profile, if it exists
    if( lCopyMe->GetProtonProfileToCopy() ){
        fProtonProfile = (TProfile*) lCopyMe->GetProtonProfileToCopy()->Clone(Form("fProtonProfile_%s",GetName()));
    }
}

//________________________________________________________________
AliV0Result::~AliV0Result(){
    // Proper destructor: delete pointer data member
    if (fHisto) {
        delete fHisto;
        fHisto = 0x0;
    }
    if (fHistoFeeddown) {
        delete fHistoFeeddown;
        fHistoFeeddown = 0x0;
    }
    if (fProtonProfile) {
        delete fProtonProfile;
        fProtonProfile = 0x0;
    }
}

//________________________________________________________________
AliV0Result& AliV0Result::operator=(const AliV0Result& lCopyMe)
{
    if (&lCopyMe == this) return *this;
    //Careful with names
    SetName(lCopyMe.GetName());
    SetTitle(lCopyMe.GetTitle());
    
    fMassHypo = lCopyMe.GetMassHypothesis();
    //Binning matters
    fhNCentBounds = lCopyMe.GetNCentBins()+1;
    fhNPtBounds = lCopyMe.GetNPtBins()+1;
    fhNPtBoundsFeeddown = lCopyMe.GetNPtBinsFeeddown()+1;
    fhNMassBins = lCopyMe.GetNMassBins();
    fhMinMass = lCopyMe.GetMinMass();
    fhMaxMass = lCopyMe.GetMaxMass();
    
    //Acceptance cuts
    fCutMinRapidity = lCopyMe.GetCutMinRapidity();
    fCutMaxRapidity = lCopyMe.GetCutMaxRapidity(),
    
    fCutV0Radius = lCopyMe.GetCutV0Radius();
    fCutMaxV0Radius = lCopyMe.GetCutMaxV0Radius();
    fCutDCANegToPV = lCopyMe.GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe.GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe.GetCutDCAV0Daughters();
    fCutV0CosPA = lCopyMe.GetCutV0CosPA();
    fCutProperLifetime = lCopyMe.GetCutProperLifetime();
    fCutCompetingV0Rejection = lCopyMe.GetCutCompetingV0Rejection();
    fCutArmenteros = lCopyMe.GetCutArmenteros();
    fCutArmenterosParameter = lCopyMe.GetCutArmenterosParameter();
    fCutTPCdEdx = lCopyMe.GetCutTPCdEdx();
    fCutMinBaryonMomentum = lCopyMe.GetCutMinBaryonMomentum();
    
    fCutMCPhysicalPrimary = lCopyMe.GetCutMCPhysicalPrimary();
    fCutMCLambdaFromPrimaryXi = lCopyMe.GetCutMCLambdaFromPrimaryXi();
    fCutMCPDGCodeAssociation = lCopyMe.GetCutMCPDGCodeAssociation();
    fCutMCUseMCProperties = lCopyMe.GetCutMCUseMCProperties();
    
    //Track Cuts
    fCutUseITSRefitTracks = lCopyMe.GetCutUseITSRefitTracks();
    fCutLeastNumberOfCrossedRows = lCopyMe.GetCutLeastNumberOfCrossedRows();
    fCutLeastNumberOfCrossedRowsOverFindable = lCopyMe.GetCutLeastNumberOfCrossedRowsOverFindable();
    fCutMinEtaTracks = lCopyMe.GetCutMinEtaTracks();
    fCutMaxEtaTracks = lCopyMe.GetCutMaxEtaTracks();
    fCutMaxChi2PerCluster = lCopyMe.GetCutMaxChi2PerCluster();
    fCutMinTrackLength = lCopyMe.GetCutMinTrackLength();
    fCutUseParametricLength = lCopyMe.GetCutUseParametricLength();
    fCutMinCrossedRowsOverLength = lCopyMe.GetCutMinCrossedRowsOverLength();
    
    //Variable V0CosPA
    fCutUseVariableV0CosPA = lCopyMe.GetCutUseVarV0CosPA();
    fCutVarV0CosPA_Exp0Const = lCopyMe.GetCutVarV0CosPAExp0Const();
    fCutVarV0CosPA_Exp0Slope = lCopyMe.GetCutVarV0CosPAExp0Slope();
    fCutVarV0CosPA_Exp1Const = lCopyMe.GetCutVarV0CosPAExp1Const();
    fCutVarV0CosPA_Exp1Slope = lCopyMe.GetCutVarV0CosPAExp1Slope();
    fCutVarV0CosPA_Const = lCopyMe.GetCutVarV0CosPAConst();
    
    //OTF use
    fUseOnTheFly = lCopyMe.GetUseOnTheFly();
    
    //special dedx
    fCut276TeVLikedEdx = lCopyMe.GetCut276TeVLikedEdx();
    
    fCutAtLeastOneTOF = lCopyMe.GetCutAtLeastOneTOF();
    fCutITSorTOF = lCopyMe.GetCutITSorTOF();
    
    fCutIsCowboy = lCopyMe.GetCutIsCowboy();
    
    if (fHisto) {
        delete fHisto;
        fHisto = 0;
    }
    if (fHistoFeeddown) {
        delete fHistoFeeddown;
        fHistoFeeddown = 0;
    }
    if (fProtonProfile) {
        delete fProtonProfile;
        fProtonProfile = 0;
    }
    
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
    if ( lCopyMe.fhPtBinsFeeddown){
        //momentum binning assignment
        fhPtBinsFeeddown = new Double_t[fhNPtBoundsFeeddown];
        for(Long_t ibin=0; ibin<fhNPtBoundsFeeddown; ibin++)
            fhPtBinsFeeddown[ibin] = lCopyMe.fhPtBinsFeeddown[ibin];
    }
    fHisto=0x0;
    if (lCopyMe.GetHistogramToCopy() ){
        //Main output histogram
        fHisto = (TH3F*) lCopyMe.GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    }
    //Copy feeddown matrix, if it exists
    fHistoFeeddown = 0x0;
    if( lCopyMe.GetHistogramFeeddownToCopy() )
        fHistoFeeddown = (TH3F*) lCopyMe.GetHistogramFeeddownToCopy()->Clone(Form("fHistoFeeddown_%s",GetName()));
    fProtonProfile = 0x0;
    //Copy proton profile, if it exists
    if( lCopyMe.GetProtonProfileToCopy() ){
        fProtonProfile = (TProfile*) lCopyMe.GetProtonProfileToCopy()->Clone(Form("fProtonProfile_%s",GetName()));
    }
    
    return *this;
}

//________________________________________________________________
Long64_t AliV0Result::Merge(TCollection *hlist)
//Merging function to allow for usage as analysis output
{
    
    if (!hlist) return 0;
    if (hlist->IsEmpty()) return (Long64_t) GetHistogram()->GetEntries();
    
    if (hlist) {
        AliV0Result *xh = 0;
        TIter nxh(hlist);
        while ((xh = (AliV0Result *) nxh())) {
            // Check if you're not committing a crime
            if( ! HasSameCuts( xh ) ){
                AliFatal("FATAL: you're trying to sum output that was obtained with different selections!");
            }
            //... if all fine, add this histogram
            GetHistogram()->Add(xh->GetHistogram());
            
            //... if feeddown matrices are both defined, merge that as well, please
            if ( fHistoFeeddown && xh->GetHistogramFeeddown() )
                GetHistogramFeeddown()->Add(xh->GetHistogramFeeddown());
            
            //... if proton profiles are both defined, merge that as well, please
            if ( fProtonProfile && xh->GetProtonProfileToCopy() )
                GetProtonProfile()->Add(xh->GetProtonProfile());
        }
    }
    return (Long64_t) GetHistogram()->GetEntries();
}
//________________________________________________________________
Bool_t AliV0Result::HasSameCuts(AliVWeakResult *lCompare, Bool_t lCheckdEdx )
//Function to compare the cuts contained in this result with another
//Returns kTRUE if all selection cuts are identical within 1e-6
//WARNING: Does not check MC association flags
{
    Bool_t lReturnValue = kTRUE;
    
    if( !lCompare->InheritsFrom(AliV0Result::Class() ) ){
        //Apples and oranges! Return kFALSE
        return kFALSE;
    }
    
    AliV0Result *lCompareV0 = (AliV0Result*) lCompare;
    
    if( fMassHypo != lCompareV0->GetMassHypothesis() ) lReturnValue = kFALSE;
    
    //Acceptance
    if( TMath::Abs( fCutMinRapidity - lCompareV0->GetCutMinRapidity() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMaxRapidity - lCompareV0->GetCutMaxRapidity() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //V0 Selection Criteria
    if( TMath::Abs( fCutDCANegToPV - lCompareV0->GetCutDCANegToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCAPosToPV - lCompareV0->GetCutDCAPosToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCAV0Daughters - lCompareV0->GetCutDCAV0Daughters() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0CosPA - lCompareV0->GetCutV0CosPA() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0Radius - lCompareV0->GetCutV0Radius() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMaxV0Radius - lCompareV0->GetCutMaxV0Radius() ) > 1e-6 ) lReturnValue = kFALSE;
    
    if( TMath::Abs( fCutProperLifetime - lCompareV0->GetCutProperLifetime() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //if( fCutCompetingV0Rejection != lCompareV0->GetCutCompetingV0Rejection() ) lReturnValue = kFALSE;
    if( fCutArmenteros != lCompareV0->GetCutArmenteros() ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutArmenterosParameter - lCompareV0->GetCutArmenterosParameter() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutTPCdEdx - lCompareV0->GetCutTPCdEdx() ) > 1e-6 && lCheckdEdx ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMinBaryonMomentum - lCompareV0->GetCutMinBaryonMomentum() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //Track Selections
    if( TMath::Abs( fCutUseITSRefitTracks - lCompareV0->GetCutUseITSRefitTracks() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutLeastNumberOfCrossedRows - lCompareV0->GetCutLeastNumberOfCrossedRows() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutLeastNumberOfCrossedRowsOverFindable - lCompareV0->GetCutLeastNumberOfCrossedRowsOverFindable() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMinEtaTracks - lCompareV0->GetCutMinEtaTracks() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMaxEtaTracks - lCompareV0->GetCutMaxEtaTracks() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMaxChi2PerCluster - lCompareV0->GetCutMaxChi2PerCluster() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMinTrackLength - lCompareV0->GetCutMinTrackLength() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutUseParametricLength - lCompareV0->GetCutUseParametricLength() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMinCrossedRowsOverLength - lCompareV0->GetCutMinCrossedRowsOverLength() ) > 1e-6 ) lReturnValue = kFALSE;
    
    //Variable V0CosPA
    if ( TMath::Abs(fCutUseVariableV0CosPA - lCompareV0->GetCutUseVarV0CosPA()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp0Const - lCompareV0->GetCutVarV0CosPAExp0Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp0Slope - lCompareV0->GetCutVarV0CosPAExp0Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp1Const - lCompareV0->GetCutVarV0CosPAExp1Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp1Slope - lCompareV0->GetCutVarV0CosPAExp1Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Const  - lCompareV0->GetCutVarV0CosPAConst()) > 1e-6 ) lReturnValue = kFALSE;
    
    //Use OTF
    if( fUseOnTheFly != lCompareV0->GetUseOnTheFly() ) lReturnValue = kFALSE;
    
    //special dedx
    if( fCut276TeVLikedEdx != lCompareV0->GetCut276TeVLikedEdx() ) lReturnValue = kFALSE;
    
    if( fCutAtLeastOneTOF != lCompareV0->GetCutAtLeastOneTOF() ) lReturnValue = kFALSE;
    if( fCutITSorTOF != lCompareV0->GetCutITSorTOF() ) lReturnValue = kFALSE;
    
    if( fCutIsCowboy != lCompareV0->GetCutIsCowboy() ) lReturnValue = kFALSE;
    
    return lReturnValue;
}
//________________________________________________________________
void AliV0Result::Print()
//Function to compare the cuts contained in this result with another
//Returns kTRUE if all selection cuts are identical within 1e-6
//WARNING: Does not check MC association flags
{
    cout<<"========================================"<<endl;
    cout<<"    AliV0Result Configuration      "<<endl;
    cout<<"========================================"<<endl;
    cout<<" Object Name........: "<<this->GetName()<<endl;
    cout<<" Use OTF V0s........: "<<fUseOnTheFly<<endl;
    cout<<" 2.76TeV-like dE/dx.: "<<fCut276TeVLikedEdx<<endl;
    if( fMassHypo == AliV0Result::kK0Short      ) cout<<" Mass Hypothesis....: K0Short"<<endl;
    if( fMassHypo == AliV0Result::kLambda       ) cout<<" Mass Hypothesis....: Lambda"<<endl;
    if( fMassHypo == AliV0Result::kAntiLambda   ) cout<<" Mass Hypothesis....: AntiLambda"<<endl;
    cout<<" Expected Mass......: "<<GetMass()<<endl;
    cout<<" Min y..............: "<<fCutMinRapidity<<endl;
    cout<<" Max y..............: "<<fCutMaxRapidity<<endl;
    cout<<" DCA Neg to PV......: "<<fCutDCANegToPV<<endl;
    cout<<" DCA Pos to PV......: "<<fCutDCAPosToPV<<endl;
    cout<<" DCA V0 Daughters...: "<<fCutDCAV0Daughters<<endl;
    cout<<" V0 Cos PA..........: "<<fCutV0CosPA<<endl;
    cout<<" Use Var V0CosPA....: "<<fCutUseVariableV0CosPA<<endl;
    if( fCutUseVariableV0CosPA ){
        cout<<" ^--Exp. 0 Const....: "<<fCutVarV0CosPA_Exp0Const<<endl;
        cout<<" ^--Exp. 0 Slope....: "<<fCutVarV0CosPA_Exp0Slope<<endl;
        cout<<" ^--Exp. 1 Const....: "<<fCutVarV0CosPA_Exp1Const<<endl;
        cout<<" ^--Exp. 1 Slope....: "<<fCutVarV0CosPA_Exp1Slope<<endl;
        cout<<" ^--Constant........: "<<fCutVarV0CosPA_Const<<endl;
    }
    cout<<" V0 2D Radius.......: "<<fCutV0Radius<<endl;
    cout<<" V0 2D Radius (max).: "<<fCutMaxV0Radius<<endl;
    
    cout<<" Proper Lifetime....: "<<fCutProperLifetime<<endl;
    cout<<" Armenteros (for K0): "<<fCutArmenteros<<endl;
    cout<<" Armenteros param...: "<<fCutArmenterosParameter<<endl;
    cout<<" TPC dEdx (sigmas)..: "<<fCutTPCdEdx<<endl;
    cout<<" Min baryon momentum: "<<fCutMinBaryonMomentum<<endl;
    
    cout<<" MC PDG Association.: "<<fCutMCPDGCodeAssociation<<endl;
    cout<<" MC Phys Primary....: "<<fCutMCPhysicalPrimary<<endl;
    cout<<" Lambda from Xi.....: "<<fCutMCLambdaFromPrimaryXi<<endl;
    cout<<" MC Use MC pT, y....: "<<fCutMCUseMCProperties<<endl;
    cout<<" Use ITSref tracks..: "<<fCutUseITSRefitTracks<<endl;
    cout<<" Nbr Crossed Rows...: "<<fCutLeastNumberOfCrossedRows<<endl;
    cout<<" Nbr Crsd Rows / Fdb: "<<fCutLeastNumberOfCrossedRowsOverFindable<<endl;
    cout<<" Min track eta......: "<<fCutMinEtaTracks<<endl;
    cout<<" Max track eta......: "<<fCutMaxEtaTracks<<endl;
    cout<<" Max chi2/clusters..: "<<fCutMaxChi2PerCluster<<endl;
    cout<<" Min Track Length...: "<<fCutMinTrackLength<<endl;
    cout<<" At least 1 tof.....: "<<fCutAtLeastOneTOF<<endl;
    cout<<" ITS||TOF...........: "<<fCutITSorTOF<<endl; 
    cout<<" Is cowboy..........: "<<fCutIsCowboy<<endl;
    cout<<"========================================"<<endl;
    return;
}

//________________________________________________________________
Double_t AliV0Result::GetMass () const
//Get Mass under expected hypothesis
//N.B. masses are rounded within 1MeV/c^2 just to simplify binning
{
    Double_t lReturnValue = 0;
    
    if( fMassHypo == AliV0Result::kK0Short    ) lReturnValue = 0.498;
    if( fMassHypo == AliV0Result::kLambda     ) lReturnValue = 1.116;
    if( fMassHypo == AliV0Result::kAntiLambda ) lReturnValue = 1.116;
    return lReturnValue;
}

//________________________________________________________________
TString AliV0Result::GetParticleName () const
//Get particle name
{
    TString lName = "";
    if( fMassHypo == AliV0Result::kK0Short    ) lName = "K0Short";
    if( fMassHypo == AliV0Result::kLambda     ) lName = "Lambda";
    if( fMassHypo == AliV0Result::kAntiLambda ) lName = "AntiLambda";
    return lName;
}

//________________________________________________________________
void AliV0Result::InitializeHisto ()
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
void AliV0Result::SetupFeeddownMatrix (Long_t lNXiPtPins, Double_t *lXiPtPins)
//Setup arrays containing FD matrix range for cascades (Xi->Lambda feeddown correction) 
{
    fhNPtBoundsFeeddown = lNXiPtPins+1;
    fhPtBinsFeeddown = new Double_t[fhNPtBoundsFeeddown];
    for(Int_t i=0; i<fhNPtBoundsFeeddown; i++) fhPtBinsFeeddown[i] = lXiPtPins[i];
}

//________________________________________________________________
void AliV0Result::InitializeFeeddownMatrix(Long_t lNLambdaPtBins, Double_t *lLambdaPtBins,
                                           Long_t lNXiPtPins, Double_t *lXiPtPins,
                                           Long_t lNCentBins, Double_t *lCentBins)
//Initialize feeddown matrix
{
    if( fMassHypo == AliV0Result::kK0Short){
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"    Cannot set up feeddown matrix for K0Short, exiting!"<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        return;
    }
    
    if( fHistoFeeddown ){
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"       Feeddown matrix already exists, please check! "<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        return;
    }
    
    //Initialize
    fHistoFeeddown = new TH3F( Form("fHistoFeeddown_%s",GetName()), "", lNLambdaPtBins, lLambdaPtBins, lNXiPtPins, lXiPtPins, lNCentBins, lCentBins );
}

//________________________________________________________________
void AliV0Result::InitializeFeeddownMatrix ()
//Initialize TProfile to do bookkeeping of proton momenta
{
    if( fMassHypo == AliV0Result::kK0Short){
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"    Cannot set up feeddown matrix for K0Short, exiting!"<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        return;
    }
    
    if( fHistoFeeddown ){
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"       Feeddown matrix already exists, please check! "<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        return;
    }
    if( fhNPtBoundsFeeddown<0 ){
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cout<<"       Can't initialize FD matrix without Xi binning! "<<endl;
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        return;
    }
    
    //Initialize
    fHistoFeeddown = new TH3F( Form("fHistoFeeddown_%s",GetName()), "", fhNPtBounds-1, fhPtBins, fhNPtBoundsFeeddown-1, fhPtBinsFeeddown, fhNCentBounds-1, fhCentBins );
}

//________________________________________________________________
void AliV0Result::InitializeProtonProfile ()
//Initialize TProfile to do bookkeeping of proton momenta
{
    if(!fProtonProfile) fProtonProfile = new TProfile( Form("fProtonProfile_%s",GetName()), "", fhNPtBounds-1, fhPtBins);
}

//________________________________________________________________
void AliV0Result::InitializeProtonProfile (Long_t lNPtBins, Double_t *lPtBins) //kept for compatibility
//Initialize TProfile to do bookkeeping of proton momenta
{
    if(!fProtonProfile) fProtonProfile = new TProfile( Form("fProtonProfile_%s",GetName()), "", lNPtBins, lPtBins);
}






