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
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
fCutV0Radius(5.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
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
fHistoFeeddown(0),
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fUseOnTheFly(kFALSE)
{
    // Dummy Constructor - not to be used! 
    //Main output histogram: Centrality, mass, transverse momentum
    //Warning: This has super-fine binning in all dimensions!
    //It may be quite costly to use, memory-wise!
    fHisto = new TH3F("fHisto","", 100,0,100, 200,0,20, 400,0,2);
    fHisto->Sumw2();
}
//________________________________________________________________
AliV0Result::AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
fCutV0Radius(5.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
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
fHistoFeeddown(0), //do not initialize by default
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fUseOnTheFly(kFALSE)
{
    // Constructor
    Double_t lThisMass = GetMass();
    Double_t lMassWindow = 0.1; // Default : good for Lambdas, not good for K0
    
    if( lMassHypo == AliV0Result::kK0Short      ){
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    
    //Main output histogram: Centrality, mass, transverse momentum
    //Warning: This has super-fine binning in all dimensions!
    //It may be quite costly to use, memory-wise!
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 100,0,100, 200,0,20, 400,lThisMass-lMassWindow,lThisMass+lMassWindow);
    fHisto->Sumw2();
}
//________________________________________________________________
AliV0Result::AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
fCutV0Radius(5.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
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
fHistoFeeddown(0), //do not initialize by default
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fUseOnTheFly(kFALSE)
{
    // Constructor
    Double_t lThisMass = GetMass();
    Double_t lMassWindow = 0.1 ;
    
    if( lMassHypo == AliV0Result::kK0Short      ){
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    
    //Construct binning in invariant mass as standard: 400 bins from lThisMass-0.1 to lThisMass+1
    const Long_t lNMassBins = 400;
    Double_t lMassDelta = (lMassWindow * 2.) / lNMassBins;
    Double_t lMassBins[lNMassBins+1];
    
    for( Long_t ibound = 0; ibound<lNMassBins+1; ibound++) lMassBins[ibound] = (lThisMass-lMassWindow) + ( ( (Double_t) ibound )*lMassDelta );
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", lNCentBins, lCentBins, lNPtBins, lPtBins, lNMassBins, lMassBins );
    fHisto->Sumw2();
}
//________________________________________________________________
AliV0Result::AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins, Long_t lNMassBins, Double_t lMinMass, Double_t lMaxMass):
AliVWeakResult(name,title),
fMassHypo(lMassHypo),
fCutMinRapidity(-0.5),
fCutMaxRapidity(+0.5),
fCutV0Radius(5.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
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
fHistoFeeddown(0), //do not initialize by default
fCutUseVariableV0CosPA(kFALSE),
fCutVarV0CosPA_Exp0Const(0),
fCutVarV0CosPA_Exp0Slope(0),
fCutVarV0CosPA_Exp1Const(0),
fCutVarV0CosPA_Exp1Slope(0),
fCutVarV0CosPA_Const(1),
fUseOnTheFly(kFALSE)
{
    // Constructor
    Double_t lMassWindow = (lMaxMass-lMinMass)/2.0 ;
    
    //Construct binning in invariant mass as standard: 400 bins from lThisMass-0.1 to lThisMass+1
    const Long_t lNMassBinsConst = lNMassBins;
    Double_t lMassDelta = (lMassWindow * 2.) / lNMassBins;
    Double_t lMassBins[lNMassBinsConst+1];
    
    for( Long_t ibound = 0; ibound<lNMassBinsConst+1; ibound++) lMassBins[ibound] = lMinMass + ( ( (Double_t) ibound )*lMassDelta );
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", lNCentBins, lCentBins, lNPtBins, lPtBins, lNMassBins, lMassBins );
    fHisto->Sumw2();
}
//________________________________________________________________
AliV0Result::AliV0Result(const AliV0Result& lCopyMe, TString lNewName)
: AliVWeakResult(lCopyMe),
fMassHypo(lCopyMe.fMassHypo),
//Acceptance Cuts
fCutMinRapidity(lCopyMe.fCutMinRapidity),
fCutMaxRapidity(lCopyMe.fCutMaxRapidity),
//Topological
fCutV0Radius(lCopyMe.fCutV0Radius),
fCutDCANegToPV(lCopyMe.fCutDCANegToPV),
fCutDCAPosToPV(lCopyMe.fCutDCAPosToPV),
fCutDCAV0Daughters(lCopyMe.fCutDCAV0Daughters),
fCutV0CosPA(lCopyMe.fCutV0CosPA),
fCutProperLifetime(lCopyMe.fCutProperLifetime),
fCutCompetingV0Rejection(lCopyMe.fCutCompetingV0Rejection),
fCutArmenteros(lCopyMe.fCutArmenteros),
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

fCutUseVariableV0CosPA(lCopyMe.fCutUseVariableV0CosPA),
fCutVarV0CosPA_Exp0Const(lCopyMe.fCutVarV0CosPA_Exp0Const),
fCutVarV0CosPA_Exp0Slope(lCopyMe.fCutVarV0CosPA_Exp0Slope),
fCutVarV0CosPA_Exp1Const(lCopyMe.fCutVarV0CosPA_Exp1Const),
fCutVarV0CosPA_Exp1Slope(lCopyMe.fCutVarV0CosPA_Exp1Slope),
fCutVarV0CosPA_Const(lCopyMe.fCutVarV0CosPA_Const),
fUseOnTheFly(lCopyMe.fUseOnTheFly)
{
    SetName( lNewName.Data() ); 
    
    // Constructor
    Double_t lThisMass = GetMass();
    Double_t lMassWindow = 0.1 ;
    
    if( fMassHypo == AliV0Result::kK0Short      ){
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    
    //Main output histogram: Centrality, mass, transverse momentum: Clone from copied object
    fHisto = (TH3F*) lCopyMe.GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    
    //Copy feeddown matrix, if it exists
    fHistoFeeddown = 0x0;
    if( lCopyMe.GetHistogramFeeddownToCopy() )
        fHistoFeeddown = (TH3F*) lCopyMe.GetHistogramFeeddownToCopy()->Clone(Form("fHistoFeeddown_%s",GetName()));
}
//________________________________________________________________
AliV0Result::AliV0Result(AliV0Result *lCopyMe, TString lNewName)
    : AliVWeakResult(*lCopyMe),
      fHisto(0)
{
    SetName(lNewName.Data()); 
    fMassHypo = lCopyMe->GetMassHypothesis();
    
    //Acceptance Cuts
    fCutMinRapidity     = lCopyMe->GetCutMinRapidity();
    fCutMaxRapidity     = lCopyMe->GetCutMaxRapidity();
    
    fCutV0Radius = lCopyMe->GetCutV0Radius();
    fCutDCANegToPV = lCopyMe->GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe->GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe->GetCutDCAV0Daughters();
    fCutV0CosPA = lCopyMe->GetCutV0CosPA();
    fCutProperLifetime = lCopyMe->GetCutProperLifetime();
    fCutCompetingV0Rejection = lCopyMe->GetCutCompetingV0Rejection();
    fCutArmenteros = lCopyMe->GetCutArmenteros();
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
    
    //Variable V0CosPA
    fCutUseVariableV0CosPA = lCopyMe -> GetCutUseVarV0CosPA();
    fCutVarV0CosPA_Exp0Const = lCopyMe -> GetCutVarV0CosPAExp0Const();
    fCutVarV0CosPA_Exp0Slope = lCopyMe -> GetCutVarV0CosPAExp0Slope();
    fCutVarV0CosPA_Exp1Const = lCopyMe -> GetCutVarV0CosPAExp1Const();
    fCutVarV0CosPA_Exp1Slope = lCopyMe -> GetCutVarV0CosPAExp1Slope();
    fCutVarV0CosPA_Const = lCopyMe -> GetCutVarV0CosPAConst();
    
    //OTF use
    fUseOnTheFly = lCopyMe -> GetUseOnTheFly();
    
    // Constructor
    Double_t lThisMass = GetMass();
    Double_t lMassWindow = 0.1 ;
    
    if( fMassHypo == AliV0Result::kK0Short      ){
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    
    //Main output histogram: Centrality, mass, transverse momentum: Clone from copied object
    fHisto = (TH3F*) lCopyMe->GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    
    //Copy feeddown matrix, if it exists
    fHistoFeeddown = 0x0;
    if( lCopyMe->GetHistogramFeeddownToCopy() )
        fHistoFeeddown = (TH3F*) lCopyMe->GetHistogramFeeddownToCopy()->Clone(Form("fHistoFeeddown_%s",GetName()));
}
//________________________________________________________________
AliV0Result::~AliV0Result(){
    // Proper destructor: delete pointer data member
    if (fHisto) {
        delete fHisto;
        fHisto = 0x0;
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
    
    //Acceptance cuts
    fCutMinRapidity = lCopyMe.GetCutMinRapidity();
    fCutMaxRapidity = lCopyMe.GetCutMaxRapidity(),
    
    fCutV0Radius = lCopyMe.GetCutV0Radius();
    fCutDCANegToPV = lCopyMe.GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe.GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe.GetCutDCAV0Daughters();
    fCutV0CosPA = lCopyMe.GetCutV0CosPA();
    fCutProperLifetime = lCopyMe.GetCutProperLifetime();
    fCutCompetingV0Rejection = lCopyMe.GetCutCompetingV0Rejection();
    fCutArmenteros = lCopyMe.GetCutArmenteros();
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
    
    //Variable V0CosPA
    fCutUseVariableV0CosPA = lCopyMe.GetCutUseVarV0CosPA();
    fCutVarV0CosPA_Exp0Const = lCopyMe.GetCutVarV0CosPAExp0Const();
    fCutVarV0CosPA_Exp0Slope = lCopyMe.GetCutVarV0CosPAExp0Slope();
    fCutVarV0CosPA_Exp1Const = lCopyMe.GetCutVarV0CosPAExp1Const();
    fCutVarV0CosPA_Exp1Slope = lCopyMe.GetCutVarV0CosPAExp1Slope();
    fCutVarV0CosPA_Const = lCopyMe.GetCutVarV0CosPAConst();
    
    //OTF use
    fUseOnTheFly = lCopyMe.GetUseOnTheFly();
    
    if (fHisto) {
        delete fHisto;
        fHisto = 0;
    }
    // Constructor
    Double_t lThisMass = GetMass();
    Double_t lMassWindow = 0.1 ;
    
    if( fMassHypo == AliV0Result::kK0Short      ){
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    
    //Main output histogram: Centrality, mass, transverse momentum: Clone from copied object
    fHisto = (TH3F*) lCopyMe.GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    
    //Copy feeddown matrix, if it exists
    fHistoFeeddown = 0x0;
    if( lCopyMe.GetHistogramFeeddownToCopy() )
        fHistoFeeddown = (TH3F*) lCopyMe.GetHistogramFeeddownToCopy()->Clone(Form("fHistoFeeddown_%s",GetName()));
    
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
    
    if( TMath::Abs( fCutProperLifetime - lCompareV0->GetCutProperLifetime() ) > 1e-6 ) lReturnValue = kFALSE;

    //if( fCutCompetingV0Rejection != lCompareV0->GetCutCompetingV0Rejection() ) lReturnValue = kFALSE;
    if( fCutArmenteros != lCompareV0->GetCutArmenteros() ) lReturnValue = kFALSE;
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
    
    //Variable V0CosPA
    if ( TMath::Abs(fCutUseVariableV0CosPA - lCompareV0->GetCutUseVarV0CosPA()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp0Const - lCompareV0->GetCutVarV0CosPAExp0Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp0Slope - lCompareV0->GetCutVarV0CosPAExp0Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp1Const - lCompareV0->GetCutVarV0CosPAExp1Const()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Exp1Slope - lCompareV0->GetCutVarV0CosPAExp1Slope()) > 1e-6 ) lReturnValue = kFALSE;
    if ( TMath::Abs(fCutVarV0CosPA_Const  - lCompareV0->GetCutVarV0CosPAConst()) > 1e-6 ) lReturnValue = kFALSE;
    
    //Use OTF
    if( fUseOnTheFly != lCompareV0->GetUseOnTheFly() ) lReturnValue = kFALSE;
    
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
    cout<<" Histogram Name.....: "<<fHisto->GetName()<<endl;
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
    
    cout<<" Proper Lifetime....: "<<fCutProperLifetime<<endl;
    cout<<" Armenteros (for K0): "<<fCutArmenteros<<endl;
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
    cout<<"========================================"<<endl;
    return;
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





