//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "TList.h"
#include "TH3F.h"
#include "AliCascadeResult.h"
#include "AliLog.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliCascadeResult);
//________________________________________________________________
AliCascadeResult::AliCascadeResult() :
  TNamed(),
fMassHypo(AliCascadeResult::kXiMinus),
//V0 Cuts
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutV0Radius(5.0),
//Cascade Cuts
fCutDCAV0ToPV(0.05),
fCutV0Mass(0.010),
fCutDCABachToPV(0.03),
fCutDCACascDaughters(2.0),
fCutCascCosPA(0.95),
fCutCascRadius(0.4),
//Miscellaneous
fCutProperLifetime(1000),
fCutLeastNumberOfClusters(70),
fCutTPCdEdx(4.0),
fCutXiRejection(0.008),
fCutMCPhysicalPrimary(kTRUE),
fCutMCPDGCodeAssociation(kTRUE)
{
    // Dummy Constructor - not to be used!
    //Main output histogram: Centrality, pt, mass
    fHisto = new TH3F("fHisto","", 20,0,100, 200,0,20, 400,0,2 );
    fHisto->Sumw2();
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(const char * name, AliCascadeResult::EMassHypo lMassHypo, const char * title):
TNamed(name,title),
fMassHypo(lMassHypo),
//V0 Cuts
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutV0Radius(5.0),
//Cascade Cuts
fCutDCAV0ToPV(0.05),
fCutV0Mass(0.010),
fCutDCABachToPV(0.03),
fCutDCACascDaughters(2.0),
fCutCascCosPA(0.95),
fCutCascRadius(0.4),
//Miscellaneous
fCutProperLifetime(1000),
fCutLeastNumberOfClusters(70),
fCutTPCdEdx(4.0),
fCutXiRejection(0.008),
fCutMCPhysicalPrimary(kTRUE),
fCutMCPDGCodeAssociation(kTRUE)
{
    // Constructor
    Double_t lThisMass = 0;
    if( fMassHypo == AliCascadeResult::kXiMinus      ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kXiPlus       ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) lThisMass = 1.67245;
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) lThisMass = 1.67245;
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,0,20, 400,lThisMass-0.1,lThisMass+0.1);
    fHisto->Sumw2();
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(const char * name, AliCascadeResult::EMassHypo lMassHypo, const char * title, Long_t lNCentBins, Double_t *lCentBins, Long_t lNPtBins, Double_t *lPtBins):
TNamed(name,title),
fMassHypo(lMassHypo),
//V0 Cuts
fCutDCANegToPV(0.1),
fCutDCAPosToPV(0.1),
fCutDCAV0Daughters(1.0),
fCutV0CosPA(0.998),
fCutV0Radius(5.0),
//Cascade Cuts
fCutDCAV0ToPV(0.05),
fCutV0Mass(0.010),
fCutDCABachToPV(0.03),
fCutDCACascDaughters(2.0),
fCutCascCosPA(0.95),
fCutCascRadius(0.4),
//Miscellaneous
fCutProperLifetime(1000),
fCutLeastNumberOfClusters(70),
fCutTPCdEdx(4.0),
fCutXiRejection(0.008),
fCutMCPhysicalPrimary(kTRUE),
fCutMCPDGCodeAssociation(kTRUE)
{
    // Constructor
    Double_t lThisMass = 0;
    if( fMassHypo == AliCascadeResult::kXiMinus      ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kXiPlus       ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) lThisMass = 1.67245;
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) lThisMass = 1.67245;
    
    //Construct binning in invariant mass as standard: 400 bins from lThisMass-0.1 to lThisMass+1
    const Long_t lNMassBins = 400;
    
    Double_t lMassWindow = 0.1 ;
    Double_t lMassDelta = (lMassWindow * 2.) / lNMassBins;
    Double_t lMassBins[lNMassBins+1];
    
    for( Long_t ibound = 0; ibound<lNMassBins+1; ibound++) lMassBins[ibound] = (lThisMass-0.1) + ( ( (Double_t) ibound )*lMassDelta );
    
    //Main output histogram: Centrality, mass, transverse momentum: Variable binning
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", lNCentBins, lCentBins, lNPtBins, lPtBins, lNMassBins, lMassBins );
    fHisto->Sumw2();
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(const AliCascadeResult& lCopyMe)
: TNamed(lCopyMe),
fMassHypo(lCopyMe.fMassHypo),
//V0 Cuts
fCutDCANegToPV(lCopyMe.fCutDCANegToPV),
fCutDCAPosToPV(lCopyMe.fCutDCAPosToPV),
fCutDCAV0Daughters(lCopyMe.fCutDCAV0Daughters),
fCutV0CosPA(lCopyMe.fCutV0CosPA),
fCutV0Radius(lCopyMe.fCutV0Radius),
//Cascade Cuts
fCutDCAV0ToPV(lCopyMe.fCutDCAV0ToPV),
fCutV0Mass(lCopyMe.fCutV0Mass),
fCutDCABachToPV(lCopyMe.fCutDCABachToPV),
fCutDCACascDaughters(lCopyMe.fCutDCACascDaughters),
fCutCascCosPA(lCopyMe.fCutCascCosPA),
fCutCascRadius(lCopyMe.fCutCascRadius),
//Miscellaneous
fCutProperLifetime(lCopyMe.fCutProperLifetime),
fCutLeastNumberOfClusters(lCopyMe.fCutLeastNumberOfClusters),
fCutTPCdEdx(lCopyMe.fCutTPCdEdx),
fCutXiRejection(0.008),
//MC specific
fCutMCPhysicalPrimary(lCopyMe.fCutMCPhysicalPrimary),
fCutMCPDGCodeAssociation(lCopyMe.fCutMCPDGCodeAssociation)
{
    // Constructor
    Double_t lThisMass = 0;
    if( fMassHypo == AliCascadeResult::kXiMinus      ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kXiPlus       ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) lThisMass = 1.67245;
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) lThisMass = 1.67245;
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,0,20, 400,lThisMass-0.1,lThisMass+0.1);
    fHisto->Sumw2();
}
//________________________________________________________________
AliCascadeResult::AliCascadeResult(AliCascadeResult *lCopyMe)
    : TNamed(*lCopyMe),
      fHisto(0)
{
    fMassHypo = lCopyMe->GetMassHypothesis();
    //V0 Cuts
    fCutDCANegToPV     = lCopyMe->GetCutDCANegToPV();
    fCutDCAPosToPV     = lCopyMe->GetCutDCAPosToPV();
    fCutDCAV0Daughters = lCopyMe->GetCutDCAV0Daughters();
    fCutV0CosPA        = lCopyMe->GetCutV0CosPA();
    fCutV0Radius       = lCopyMe->GetCutV0Radius();
    //Cascade Cuts
    fCutDCAV0ToPV = lCopyMe -> GetCutDCAV0ToPV();
    fCutV0Mass    = lCopyMe -> GetCutV0Mass();
    fCutDCABachToPV  = lCopyMe -> GetCutDCABachToPV();
    fCutDCACascDaughters = lCopyMe -> GetCutDCACascDaughters();
    fCutCascCosPA  = lCopyMe -> GetCutCascCosPA();
    fCutCascRadius = lCopyMe -> GetCutCascRadius();
    
    //Miscellaneous
    fCutProperLifetime = lCopyMe->GetCutProperLifetime();
    fCutLeastNumberOfClusters = lCopyMe->GetCutLeastNumberOfClusters();
    fCutTPCdEdx = lCopyMe->GetCutTPCdEdx();
    fCutXiRejection = lCopyMe->GetCutXiRejection();
    
    //MC specific
    fCutMCPhysicalPrimary    = lCopyMe -> GetCutMCPhysicalPrimary();
    fCutMCPDGCodeAssociation = lCopyMe -> GetCutMCPDGCodeAssociation();
    
    // Constructor
    Double_t lThisMass = 0;
    if( fMassHypo == AliCascadeResult::kXiMinus      ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kXiPlus       ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) lThisMass = 1.67245;
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) lThisMass = 1.67245;
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,0,20, 400,lThisMass-0.1,lThisMass+0.1);
    fHisto->Sumw2();
}
//________________________________________________________________
AliCascadeResult::~AliCascadeResult(){
    // Proper destructor: delete pointer data member
    if (fHisto) {
        delete fHisto;
        fHisto = 0x0;
    }
}

//________________________________________________________________
AliCascadeResult& AliCascadeResult::operator=(const AliCascadeResult& lCopyMe)
{
    if (&lCopyMe == this) return *this;
    SetName(lCopyMe.GetName());
    SetTitle(lCopyMe.GetTitle());

    fMassHypo = lCopyMe.GetMassHypothesis();
    //V0 Cuts
    fCutDCANegToPV = lCopyMe.GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe.GetCutDCAPosToPV(),
    fCutDCAV0Daughters = lCopyMe.GetCutDCAV0Daughters();
    fCutV0CosPA = lCopyMe.GetCutV0CosPA();
    fCutV0Radius = lCopyMe.GetCutV0Radius();
    //Cascade Cuts
    fCutDCAV0ToPV = lCopyMe.GetCutDCAV0ToPV();
    fCutV0Mass    = lCopyMe.GetCutV0Mass();
    fCutDCABachToPV  = lCopyMe.GetCutDCABachToPV();
    fCutDCACascDaughters = lCopyMe.GetCutDCACascDaughters();
    fCutCascCosPA  = lCopyMe.GetCutCascCosPA();
    fCutCascRadius = lCopyMe.GetCutCascRadius();
    
    //Miscellaneous
    fCutProperLifetime = lCopyMe.GetCutProperLifetime();
    fCutLeastNumberOfClusters = lCopyMe.GetCutLeastNumberOfClusters();
    fCutTPCdEdx = lCopyMe.GetCutTPCdEdx();
    fCutXiRejection = lCopyMe.GetCutXiRejection();
    
    //MC specific
    fCutMCPhysicalPrimary = lCopyMe.GetCutMCPhysicalPrimary();
    fCutMCPDGCodeAssociation = lCopyMe.GetCutMCPDGCodeAssociation();
    
    if (fHisto) {
        delete fHisto;
        fHisto = 0;
    }
    // Constructor
    Double_t lThisMass = 0;
    if( fMassHypo == AliCascadeResult::kXiMinus      ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kXiPlus       ) lThisMass = 1.32171;
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) lThisMass = 1.67245;
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) lThisMass = 1.67245;
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,0,20, 400,lThisMass-0.1,lThisMass+0.1);
    fHisto->Sumw2();
    
    return *this;
}

//________________________________________________________________
Long64_t AliCascadeResult::Merge(TCollection *hlist)
//Merging function to allow for usage as analysis output
{
    
    if (!hlist) return 0;
    if (hlist->IsEmpty()) return (Long64_t) GetHistogram()->GetEntries();
    
    if (hlist) {
        AliCascadeResult *xh = 0;
        TIter nxh(hlist);
        while ((xh = (AliCascadeResult *) nxh())) {
            // Check if you're not committing a crime
            if( ! HasSameCuts( xh ) ){
                AliFatal("FATAL: you're trying to sum output that was obtained with different selections!");
            }
            //... if all fine, add this histogram
            GetHistogram()->Add(xh->GetHistogram());
        }
    }
    return (Long64_t) GetHistogram()->GetEntries();
}
//________________________________________________________________
Bool_t AliCascadeResult::HasSameCuts(AliCascadeResult *lCompare)
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
    
    //Cascade Selection Criteria
    if( TMath::Abs( fCutDCAV0ToPV - lCompare->GetCutDCAV0ToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutV0Mass - lCompare->GetCutV0Mass() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCABachToPV - lCompare->GetCutDCABachToPV() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutDCACascDaughters - lCompare->GetCutDCACascDaughters() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutCascCosPA - lCompare->GetCutCascCosPA() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutCascRadius - lCompare->GetCutCascRadius() ) > 1e-6 ) lReturnValue = kFALSE;
    
    if( TMath::Abs( fCutProperLifetime - lCompare->GetCutProperLifetime() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutLeastNumberOfClusters - lCompare->GetCutLeastNumberOfClusters() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutTPCdEdx - lCompare->GetCutTPCdEdx() ) > 1e-6 ) lReturnValue = kFALSE;
    if( TMath::Abs( fCutXiRejection - lCompare->GetCutXiRejection() ) > 1e-6 ) lReturnValue = kFALSE;
    
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
    if( fMassHypo == AliCascadeResult::kXiMinus      ) cout<<" Mass Hypothesis....: XiMinus"<<endl;
    if( fMassHypo == AliCascadeResult::kXiPlus       ) cout<<" Mass Hypothesis....: XiPlus"<<endl;
    if( fMassHypo == AliCascadeResult::kOmegaMinus   ) cout<<" Mass Hypothesis....: OmegaMinus"<<endl;
    if( fMassHypo == AliCascadeResult::kOmegaPlus    ) cout<<" Mass Hypothesis....: OmegaPlus"<<endl;
    
    cout<<" DCA Neg to PV......: "<<fCutDCANegToPV<<endl;
    cout<<" DCA Pos to PV......: "<<fCutDCAPosToPV<<endl;
    cout<<" DCA V0 Daughters...: "<<fCutDCAV0Daughters<<endl;
    cout<<" V0 Cos PA..........: "<<fCutV0CosPA<<endl;
    cout<<" V0 2D Radius.......: "<<fCutV0Radius<<endl;
    
    cout<<" DCA V0 to PV.......: "<<fCutDCAV0ToPV<<endl;
    cout<<" V0 Mass............: "<<fCutV0Mass<<endl;
    cout<<" DCA Bach to PV.....: "<<fCutDCABachToPV<<endl;
    cout<<" DCA Casc Daughters.: "<<fCutDCACascDaughters<<endl;
    cout<<" Casc Cos PA........: "<<fCutCascCosPA<<endl;
    cout<<" Casc 2D Radius.....: "<<fCutCascRadius<<endl;
    
    cout<<" Proper Lifetime....: "<<fCutProperLifetime<<endl;
    cout<<" Nbr Clusters.......: "<<fCutLeastNumberOfClusters<<endl;
    cout<<" TPC dEdx (sigmas)..: "<<fCutTPCdEdx<<endl;
    cout<<" Xi Rej (for Omega).: "<<fCutXiRejection<<endl;
    
    cout<<" MC PDG Association.: "<<fCutMCPDGCodeAssociation<<endl;
    cout<<" MC Phys Primary....: "<<fCutMCPhysicalPrimary<<endl;
    cout<<"========================================"<<endl;
    return;
}