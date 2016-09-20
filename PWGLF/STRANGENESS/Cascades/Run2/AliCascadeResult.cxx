//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "TList.h"
#include "TH3F.h"
#include "AliCascadeResult.h"
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
fCutMCPhysicalPrimary(kTRUE),
fCutMCPDGCodeAssociation(kTRUE)
{
    // Dummy Constructor - not to be used!
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F("fHisto","", 20,0,100, 400,0,2, 100,0,10);
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
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 400,lThisMass-0.1,lThisMass+0.1, 100,0,10);
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
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 400,lThisMass-0.1,lThisMass+0.1, 100,0,10);
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
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 400,lThisMass-0.1,lThisMass+0.1, 100,0,10);
    fHisto->Sumw2();
}
//________________________________________________________________
AliCascadeResult::~AliCascadeResult(){
    // destructor: clean stuff up
    
    // Actual deletion of the objects causes corruption of event object
    // - no idea why - on Proof(Lite). Hence it is disabled here.
    //
    //delete fEstimatorList;
    //fEstimatorList=0x0;
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
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 400,lThisMass-0.1,lThisMass+0.1, 100,0,10);
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
            // Add this histogram
            GetHistogram()->Add(xh->GetHistogram());
        }
    }
    return (Long64_t) GetHistogram()->GetEntries();
}
