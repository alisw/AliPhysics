//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold V0 configuration + results histogram
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include "TList.h"
#include "TH3F.h"
#include "AliV0Result.h"
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
fCutTPCdEdx(3.0)
{
    // Dummy Constructor - not to be used!
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F("fHisto","", 20,0,100, 200,0,2, 100,0,10);
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
fCutTPCdEdx(3.0)
{
    // Constructor
    Double_t lThisMass = 0;
    if( lMassHypo == AliV0Result::kK0Short      ) lThisMass = 0.497;
    if( lMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( lMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,lThisMass-0.1,lThisMass+0.1, 100,0,10);
}
//________________________________________________________________
AliV0Result::AliV0Result(const AliV0Result& lCopyMe)
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
fCutTPCdEdx(lCopyMe.fCutTPCdEdx)
{
    // Constructor
    Double_t lThisMass = 0;
    if( fMassHypo == AliV0Result::kK0Short      ) lThisMass = 0.497;
    if( fMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( fMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,lThisMass-0.1,lThisMass+0.1, 100,0,10);
}
//________________________________________________________________
AliV0Result::AliV0Result(AliV0Result *lCopyMe)
    : TNamed(*lCopyMe),
      fHisto(0)
{
    fMassHypo = lCopyMe->GetMassHypothesis();
    fCutV0Radius = lCopyMe->GetCutV0Radius();
    fCutDCANegToPV = lCopyMe->GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe->GetCutDCAPosToPV(),
    fCutDCAV0Daughters = lCopyMe->GetCutDCAV0Daughters(),
    fCutV0CosPA = lCopyMe->GetCutV0CosPA(),
    fCutProperLifetime = lCopyMe->GetCutProperLifetime(),
    fCutLeastNumberOfCrossedRows = lCopyMe->GetCutLeastNumberOfCrossedRows(),
    fCutLeastNumberOfCrossedRowsOverFindable = lCopyMe->GetCutLeastNumberOfCrossedRowsOverFindable(),
    fCutCompetingV0Rejection = lCopyMe->GetCutCompetingV0Rejection(),
    fCutArmenteros = lCopyMe->GetCutArmenteros();
    fCutTPCdEdx = lCopyMe->GetCutTPCdEdx();
    
    // Constructor
    Double_t lThisMass = 0;
    if( fMassHypo == AliV0Result::kK0Short      ) lThisMass = 0.497;
    if( fMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( fMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,lThisMass-0.1,lThisMass+0.1, 100,0,10);
}
//________________________________________________________________
AliV0Result::~AliV0Result(){
    // destructor: clean stuff up
    
    // Actual deletion of the objects causes corruption of event object
    // - no idea why - on Proof(Lite). Hence it is disabled here.
    //
    //delete fEstimatorList;
    //fEstimatorList=0x0;
}

//________________________________________________________________
AliV0Result& AliV0Result::operator=(const AliV0Result& lCopyMe)
{
    if (&lCopyMe == this) return *this;
    SetName(lCopyMe.GetName());
    SetTitle(lCopyMe.GetTitle());

    fMassHypo = lCopyMe.GetMassHypothesis();
    fCutV0Radius = lCopyMe.GetCutV0Radius();
    fCutDCANegToPV = lCopyMe.GetCutDCANegToPV();
    fCutDCAPosToPV = lCopyMe.GetCutDCAPosToPV(),
    fCutDCAV0Daughters = lCopyMe.GetCutDCAV0Daughters(),
    fCutV0CosPA = lCopyMe.GetCutV0CosPA(),
    fCutProperLifetime = lCopyMe.GetCutProperLifetime(),
    fCutLeastNumberOfCrossedRows = lCopyMe.GetCutLeastNumberOfCrossedRows(),
    fCutLeastNumberOfCrossedRowsOverFindable = lCopyMe.GetCutLeastNumberOfCrossedRowsOverFindable(),
    fCutCompetingV0Rejection = lCopyMe.GetCutCompetingV0Rejection(),
    fCutArmenteros = lCopyMe.GetCutArmenteros();
    fCutTPCdEdx = lCopyMe.GetCutTPCdEdx();
    
    if (fHisto) {
        delete fHisto;
        fHisto = 0;
    }
    // Constructor
    Double_t lThisMass = 0;
    if( fMassHypo == AliV0Result::kK0Short      ) lThisMass = 0.497;
    if( fMassHypo == AliV0Result::kLambda       ) lThisMass = 1.116;
    if( fMassHypo == AliV0Result::kAntiLambda   ) lThisMass = 1.116;
    
    //Main output histogram: Centrality, mass, transverse momentum
    fHisto = new TH3F(Form("fHisto_%s",GetName()),"", 20,0,100, 200,lThisMass-0.1,lThisMass+0.1, 100,0,10);
    
    return *this;
}