//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "TList.h"
#include "TH3F.h"
#include "AliV0Result.h"
#include "AliLog.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliV0Result);
//________________________________________________________________
AliV0Result::AliV0Result() :
  TNamed(),
fMassHypo(AliV0Result::kK0Short),
fCutV0Radius(5.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutLeastNumberOfCrossedRows(70),
fCutLeastNumberOfCrossedRowsOverFindable(0.8),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
fCutTPCdEdx(3.0),
fCutMinBaryonMomentum(-1),
fCutMCPhysicalPrimary(kTRUE),
fCutMCLambdaFromPrimaryXi(kFALSE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fHistoFeeddown(0)
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
TNamed(name,title),
fMassHypo(lMassHypo),
fCutV0Radius(5.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutLeastNumberOfCrossedRows(70),
fCutLeastNumberOfCrossedRowsOverFindable(0.8),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
fCutTPCdEdx(3.0),
fCutMinBaryonMomentum(-1),
fCutMCPhysicalPrimary(kTRUE),
fCutMCLambdaFromPrimaryXi(kFALSE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fHistoFeeddown(0) //do not initialize by default
{
    // Constructor
    Double_t lThisMass = 0;
    Double_t lMassWindow = 0.1; // Default : good for Lambdas, not good for K0
    
    if( lMassHypo == AliV0Result::kK0Short      ){
        lThisMass = 0.497;
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    if( lMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( lMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
    //Main output histogram: Centrality, mass, transverse momentum
    //Warning: This has super-fine binning in all dimensions!
    //It may be quite costly to use, memory-wise!
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 100,0,100, 200,0,20, 400,lThisMass-lMassWindow,lThisMass+lMassWindow);
    fHisto->Sumw2();
}
//________________________________________________________________
AliV0Result::AliV0Result(const char * name, AliV0Result::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins):
TNamed(name,title),
fMassHypo(lMassHypo),
fCutV0Radius(5.0),
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutProperLifetime(10),
fCutLeastNumberOfCrossedRows(70),
fCutLeastNumberOfCrossedRowsOverFindable(0.8),
fCutCompetingV0Rejection(-1),
fCutArmenteros(kTRUE),
fCutTPCdEdx(3.0),
fCutMinBaryonMomentum(-1),
fCutMCPhysicalPrimary(kTRUE),
fCutMCLambdaFromPrimaryXi(kFALSE),
fCutMCPDGCodeAssociation(kTRUE),
fCutMCUseMCProperties(kTRUE),
fHistoFeeddown(0) //do not initialize by default
{
    // Constructor
    Double_t lThisMass = 0;
    Double_t lMassWindow = 0.1 ;
    
    if( lMassHypo == AliV0Result::kK0Short      ){
        lThisMass = 0.497;
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    if( lMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( lMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
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
AliV0Result::AliV0Result(const AliV0Result& lCopyMe, TString lNewName)
: TNamed(lCopyMe),
fMassHypo(lCopyMe.fMassHypo),
fCutV0Radius(lCopyMe.fCutV0Radius),
fCutDCANegToPV(lCopyMe.fCutDCANegToPV),
fCutDCAPosToPV(lCopyMe.fCutDCAPosToPV),
fCutDCAV0Daughters(lCopyMe.fCutDCAV0Daughters),
fCutV0CosPA(lCopyMe.fCutV0CosPA),
fCutProperLifetime(lCopyMe.fCutProperLifetime),
fCutLeastNumberOfCrossedRows(lCopyMe.fCutLeastNumberOfCrossedRows),
fCutLeastNumberOfCrossedRowsOverFindable(lCopyMe.fCutLeastNumberOfCrossedRowsOverFindable),
fCutCompetingV0Rejection(lCopyMe.fCutCompetingV0Rejection),
fCutArmenteros(lCopyMe.fCutArmenteros),
fCutTPCdEdx(lCopyMe.fCutTPCdEdx),
fCutMCPhysicalPrimary(lCopyMe.fCutMCPhysicalPrimary),
fCutMCLambdaFromPrimaryXi(lCopyMe.fCutMCLambdaFromPrimaryXi),
fCutMCPDGCodeAssociation(lCopyMe.fCutMCPDGCodeAssociation),
fCutMCUseMCProperties(lCopyMe.fCutMCUseMCProperties)
{
    SetName( lNewName.Data() ); 
    
    // Constructor
    Double_t lThisMass = 0;
    Double_t lMassWindow = 0.1 ;
    
    if( fMassHypo == AliV0Result::kK0Short      ){
        lThisMass = 0.497;
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    if( fMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( fMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
    //Main output histogram: Centrality, mass, transverse momentum: Clone from copied object
    fHisto = (TH3F*) lCopyMe.GetHistogramToCopy()->Clone(Form("fHisto_%s",GetName()));
    
    //Copy feeddown matrix, if it exists
    fHistoFeeddown = 0x0;
    if( lCopyMe.GetHistogramFeeddownToCopy() )
        fHistoFeeddown = (TH3F*) lCopyMe.GetHistogramFeeddownToCopy()->Clone(Form("fHistoFeeddown_%s",GetName()));
}
//________________________________________________________________
AliV0Result::AliV0Result(AliV0Result *lCopyMe, TString lNewName)
    : TNamed(*lCopyMe),
      fHisto(0)
{
    SetName(lNewName.Data()); 
    fMassHypo = lCopyMe->GetMassHypothesis();
    fCutV0Radius = lCopyMe->GetCutV0Radius();
    fCutDCANegToPV = lCopyMe->GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe->GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe->GetCutDCAV0Daughters();
    fCutV0CosPA = lCopyMe->GetCutV0CosPA();
    fCutProperLifetime = lCopyMe->GetCutProperLifetime();
    fCutLeastNumberOfCrossedRows = lCopyMe->GetCutLeastNumberOfCrossedRows();
    fCutLeastNumberOfCrossedRowsOverFindable = lCopyMe->GetCutLeastNumberOfCrossedRowsOverFindable();
    fCutCompetingV0Rejection = lCopyMe->GetCutCompetingV0Rejection();
    fCutArmenteros = lCopyMe->GetCutArmenteros();
    fCutTPCdEdx = lCopyMe->GetCutTPCdEdx();
    fCutMinBaryonMomentum = lCopyMe->GetCutMinBaryonMomentum();
    
    fCutMCPhysicalPrimary = lCopyMe->GetCutMCPhysicalPrimary();
    fCutMCLambdaFromPrimaryXi = lCopyMe->GetCutMCLambdaFromPrimaryXi();
    fCutMCPDGCodeAssociation = lCopyMe->GetCutMCPDGCodeAssociation();
    fCutMCUseMCProperties    = lCopyMe -> GetCutMCUseMCProperties();
    
    // Constructor
    Double_t lThisMass = 0;
    Double_t lMassWindow = 0.1 ;
    
    if( fMassHypo == AliV0Result::kK0Short      ){
        lThisMass = 0.497;
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    if( fMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( fMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
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
    fCutV0Radius = lCopyMe.GetCutV0Radius();
    fCutDCANegToPV = lCopyMe.GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe.GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe.GetCutDCAV0Daughters();
    fCutV0CosPA = lCopyMe.GetCutV0CosPA();
    fCutProperLifetime = lCopyMe.GetCutProperLifetime();
    fCutLeastNumberOfCrossedRows = lCopyMe.GetCutLeastNumberOfCrossedRows();
    fCutLeastNumberOfCrossedRowsOverFindable = lCopyMe.GetCutLeastNumberOfCrossedRowsOverFindable(),
    fCutCompetingV0Rejection = lCopyMe.GetCutCompetingV0Rejection();
    fCutArmenteros = lCopyMe.GetCutArmenteros();
    fCutTPCdEdx = lCopyMe.GetCutTPCdEdx();
    fCutMinBaryonMomentum = lCopyMe.GetCutMinBaryonMomentum();
    
    fCutMCPhysicalPrimary = lCopyMe.GetCutMCPhysicalPrimary();
    fCutMCLambdaFromPrimaryXi = lCopyMe.GetCutMCLambdaFromPrimaryXi();
    fCutMCPDGCodeAssociation = lCopyMe.GetCutMCPDGCodeAssociation();
    fCutMCUseMCProperties = lCopyMe.GetCutMCUseMCProperties();
    
    if (fHisto) {
        delete fHisto;
        fHisto = 0;
    }
    // Constructor
    Double_t lThisMass = 0;
    Double_t lMassWindow = 0.1 ;
    
    if( fMassHypo == AliV0Result::kK0Short      ){
        lThisMass = 0.497;
        lMassWindow = 0.15; // will be 300 MeV/c^2 wide
    }
    if( fMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( fMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
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
Bool_t AliV0Result::HasSameCuts(AliV0Result *lCompare)
//Function to compare the cuts contained in this result with another
//Returns kTRUE if all selection cuts are identical within 1e-6
//WARNING: Does not check MC association flags
{
    Bool_t lReturnValue = kTRUE;
    
    if( fMassHypo != lCompare->GetMassHypothesis() ) lReturnValue = kFALSE;
    
    //V0 Selection Criteria
    if( TMath::Abs( fCutDCANegToPV - lCompare->GetCutDCANegToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCAPosToPV - lCompare->GetCutDCAPosToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCAV0Daughters - lCompare->GetCutDCAV0Daughters() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0CosPA - lCompare->GetCutV0CosPA() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0Radius - lCompare->GetCutV0Radius() ) > 1e-6 ) lReturnValue = kFALSE;
    
    if( TMath::Abs( fCutProperLifetime - lCompare->GetCutProperLifetime() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutLeastNumberOfCrossedRows - lCompare->GetCutLeastNumberOfCrossedRows() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutLeastNumberOfCrossedRowsOverFindable - lCompare->GetCutLeastNumberOfCrossedRowsOverFindable() ) > 1e-6 ) lReturnValue = kFALSE;
    //if( fCutCompetingV0Rejection != lCompare->GetCutCompetingV0Rejection() ) lReturnValue = kFALSE;
    if( fCutArmenteros != lCompare->GetCutArmenteros() ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutTPCdEdx - lCompare->GetCutTPCdEdx() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutMinBaryonMomentum - lCompare->GetCutMinBaryonMomentum() ) > 1e-6 ) lReturnValue = kFALSE;
    
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
    cout<<" Histogram Name.....: "<<fHisto->GetName()<<endl;
    if( fMassHypo == AliV0Result::kK0Short      ) cout<<" Mass Hypothesis....: K0Short"<<endl;
    if( fMassHypo == AliV0Result::kLambda       ) cout<<" Mass Hypothesis....: Lambda"<<endl;
    if( fMassHypo == AliV0Result::kAntiLambda   ) cout<<" Mass Hypothesis....: AntiLambda"<<endl;
    
    cout<<" DCA Neg to PV......: "<<fCutDCANegToPV<<endl;
    cout<<" DCA Pos to PV......: "<<fCutDCAPosToPV<<endl;
    cout<<" DCA V0 Daughters...: "<<fCutDCAV0Daughters<<endl;
    cout<<" V0 Cos PA..........: "<<fCutV0CosPA<<endl;
    cout<<" V0 2D Radius.......: "<<fCutV0Radius<<endl;
    
    cout<<" Proper Lifetime....: "<<fCutProperLifetime<<endl;
    cout<<" Nbr Crossed Rows...: "<<fCutLeastNumberOfCrossedRows<<endl;
    cout<<" Nbr Crsd Rows / Fdb: "<<fCutLeastNumberOfCrossedRowsOverFindable<<endl;
    cout<<" Armenteros (for K0): "<<fCutArmenteros<<endl;
    cout<<" TPC dEdx (sigmas)..: "<<fCutTPCdEdx<<endl;
    cout<<" Min baryon momentum: "<<fCutMinBaryonMomentum<<endl;
    
    cout<<" MC PDG Association.: "<<fCutMCPDGCodeAssociation<<endl;
    cout<<" MC Phys Primary....: "<<fCutMCPhysicalPrimary<<endl;
    cout<<" Lambda from Xi.....: "<<fCutMCLambdaFromPrimaryXi<<endl;
    cout<<" MC Use MC pT, y....: "<<fCutMCUseMCProperties<<endl;
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

