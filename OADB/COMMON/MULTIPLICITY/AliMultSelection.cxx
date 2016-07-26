/**********************************************
 *
 * Class designed to store all multiplicity
 * estimators in a TClonesArray 
 * 
 * Instancing an empty AliMultSelection class
 * will allow you to, among other things, 
 * initialize a list of standard estimators. 
 * This class will be used as a general-
 * purpose data holding class, and will also 
 * be attached to user tasks.
 *
 **********************************************/

#include "TList.h"
#include "TFormula.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCuts.h"
#include "AliMultEstimator.h"
#include <iostream>
#include <TROOT.h>
using namespace std;

ClassImp(AliMultSelection);
//________________________________________________________________
AliMultSelection::AliMultSelection() :
  TNamed(), fNEsts(0), fEvSelCode(0), fEstimatorList(0x0),
fThisEvent_VtxZCut(0),
fThisEvent_IsNotPileup(0),
fThisEvent_IsNotPileupMV(0),
fThisEvent_IsNotPileupInMultBins(0),
fThisEvent_Triggered(0),
fThisEvent_INELgtZERO(0),
fThisEvent_HasNoInconsistentVertices(0),
fThisEvent_PassesTrackletVsCluster(0),
fThisEvent_IsNotAsymmetricInVZERO(0),
fThisEvent_IsNotIncompleteDAQ(0)
{
  // Constructor
    fEstimatorList = new TList();
}
//________________________________________________________________
AliMultSelection::AliMultSelection(const char * name, const char * title):
TNamed(name,title), fNEsts(0), fEvSelCode(0), fEstimatorList(0x0),
fThisEvent_VtxZCut(0),
fThisEvent_IsNotPileup(0),
fThisEvent_IsNotPileupMV(0),
fThisEvent_IsNotPileupInMultBins(0),
fThisEvent_Triggered(0),
fThisEvent_INELgtZERO(0),
fThisEvent_HasNoInconsistentVertices(0),
fThisEvent_PassesTrackletVsCluster(0),
fThisEvent_IsNotAsymmetricInVZERO(0),
fThisEvent_IsNotIncompleteDAQ(0)
{
  // Constructor
    fEstimatorList = new TList();
    fEstimatorList->SetOwner(kTRUE);
}
//________________________________________________________________
AliMultSelection::AliMultSelection(const AliMultSelection& lCopyMe)
: TNamed(lCopyMe),
fNEsts(0),
fEvSelCode(lCopyMe.fEvSelCode),
fEstimatorList(0), 
fThisEvent_VtxZCut(lCopyMe.fThisEvent_VtxZCut),
fThisEvent_IsNotPileup(lCopyMe.fThisEvent_IsNotPileup),
fThisEvent_IsNotPileupMV(lCopyMe.fThisEvent_IsNotPileupMV),
fThisEvent_IsNotPileupInMultBins(lCopyMe.fThisEvent_IsNotPileupInMultBins),
fThisEvent_Triggered(lCopyMe.fThisEvent_Triggered),
fThisEvent_INELgtZERO(lCopyMe.fThisEvent_INELgtZERO),
fThisEvent_HasNoInconsistentVertices(lCopyMe.fThisEvent_HasNoInconsistentVertices),
fThisEvent_PassesTrackletVsCluster(lCopyMe.fThisEvent_PassesTrackletVsCluster),
      fThisEvent_IsNotAsymmetricInVZERO(lCopyMe.fThisEvent_IsNotAsymmetricInVZERO),
      fThisEvent_IsNotIncompleteDAQ(lCopyMe.fThisEvent_IsNotIncompleteDAQ)
{
    TIter next(lCopyMe.fEstimatorList);
    AliMultEstimator* est = 0;
    while ((est = static_cast<AliMultEstimator*>(next())))
        AddEstimator(new AliMultEstimator(*est));
}
//________________________________________________________________
AliMultSelection::AliMultSelection(AliMultSelection *lCopyMe)
    : TNamed(*lCopyMe),
      fNEsts(0),
      fEstimatorList(0)
{
    fEvSelCode = lCopyMe->GetEvSelCode();

    fThisEvent_VtxZCut = lCopyMe->GetThisEventVtxZCut();
    fThisEvent_IsNotPileup = lCopyMe->GetThisEventIsNotPileup();
    fThisEvent_IsNotPileupMV = lCopyMe->GetThisEventIsNotPileupMV();
    fThisEvent_IsNotPileupInMultBins = lCopyMe->GetThisEventIsNotPileupInMultBins();
    fThisEvent_Triggered = lCopyMe->GetThisEventTriggered();
    fThisEvent_INELgtZERO = lCopyMe->GetThisEventINELgtZERO();
    fThisEvent_HasNoInconsistentVertices = lCopyMe->GetThisEventHasNoInconsistentVertices();
    fThisEvent_PassesTrackletVsCluster = lCopyMe->GetThisEventPassesTrackletVsCluster();
    fThisEvent_IsNotAsymmetricInVZERO = lCopyMe->GetThisEventIsNotAsymmetricInVZERO();
    fThisEvent_IsNotIncompleteDAQ = lCopyMe->GetThisEventIsNotIncompleteDAQ();

    TIter next(lCopyMe->fEstimatorList);
    AliMultEstimator* est = 0;
    while ((est = static_cast<AliMultEstimator*>(next())))
        AddEstimator(new AliMultEstimator(*est));
}
//________________________________________________________________
void AliMultSelection::Set(AliMultSelection* s)
{
    // Printf("Setting %s from %s", GetName(), s->GetName());
    fEvSelCode = s->GetEvSelCode();

    fThisEvent_VtxZCut = s->GetThisEventVtxZCut();
    fThisEvent_IsNotPileup = s->GetThisEventIsNotPileup();
    fThisEvent_IsNotPileupMV = s->GetThisEventIsNotPileupMV();
    fThisEvent_IsNotPileupInMultBins = s->GetThisEventIsNotPileupInMultBins();
    fThisEvent_Triggered = s->GetThisEventTriggered();
    fThisEvent_INELgtZERO = s->GetThisEventINELgtZERO();
    fThisEvent_HasNoInconsistentVertices = s->GetThisEventHasNoInconsistentVertices();
    fThisEvent_PassesTrackletVsCluster = s->GetThisEventPassesTrackletVsCluster();
    fThisEvent_IsNotAsymmetricInVZERO = s->GetThisEventIsNotAsymmetricInVZERO();
    fThisEvent_IsNotIncompleteDAQ = s->GetThisEventIsNotIncompleteDAQ();
    TIter next(s->fEstimatorList);
    AliMultEstimator* e = 0;
    while ((e = static_cast<AliMultEstimator*>(next()))) {
        AliMultEstimator* ee = GetEstimator(e->GetName());
        if (!ee) {
            // Printf("Adding estimator %s", e->GetName());
            ee = new AliMultEstimator(e->GetName(),e->GetTitle(),
                                      e->GetDefinition());
            AddEstimator(ee);
        }
        ee->Set(e);
    }
}
//________________________________________________________________
AliMultSelection::~AliMultSelection(){
    // destructor: clean stuff up
    
    // Actual deletion of the objects causes corruption of event object
    // - no idea why - on Proof(Lite). Hence it is disabled here.
    //
    //delete fEstimatorList;
    //fEstimatorList=0x0;
}
//________________________________________________________________
void AliMultSelection::CleanUp()
{
    if (fEstimatorList) delete fEstimatorList;
    fEstimatorList = 0;
    fNEsts = 0;
    fEvSelCode = 0;
}
//________________________________________________________________
AliMultSelection& AliMultSelection::operator=(const AliMultSelection& lCopyMe)
{
    if (&lCopyMe == this) return *this;
    SetName(lCopyMe.GetName());
    SetTitle(lCopyMe.GetTitle());
    fNEsts = 0;
    fEvSelCode = lCopyMe.fEvSelCode; 
    if (fEstimatorList) {
        delete fEstimatorList;
        fEstimatorList = 0;
    }
    TIter next(lCopyMe.fEstimatorList);
    AliMultEstimator* est = 0;
    while ((est = static_cast<AliMultEstimator*>(next())))
        AddEstimator(new AliMultEstimator(*est));
    
    return *this;
}
//________________________________________________________________
void AliMultSelection::AddEstimator(AliMultEstimator *lEst)
{
    if (!lEst) return;
    if (!fEstimatorList) {
        fEstimatorList = new TList;
        fEstimatorList->SetOwner();
        fNEsts = 0;
    }
    fEstimatorList->Add(lEst);
    fNEsts++;
}
//________________________________________________________________
AliMultEstimator* AliMultSelection::GetEstimator (const TString& lName) const
{
    if (!fEstimatorList) return 0;
    return static_cast<AliMultEstimator*>(fEstimatorList->FindObject(lName));
}
//________________________________________________________________
AliMultEstimator* AliMultSelection::GetEstimator (Long_t lEstIdx) const
{
    if (!fEstimatorList) return 0;
    if (lEstIdx < 0 || lEstIdx >= fNEsts) return 0;
    return static_cast<AliMultEstimator*>(fEstimatorList->At(lEstIdx));
}
//________________________________________________________________
void AliMultSelection::PrintInfo()
{
    cout<<"AliMultSelection Name..: "<<GetName()<<endl;
    cout<<"Estimators defined.....: "<<fNEsts<<endl;
    cout<<"EvSel Code.............: "<<fEvSelCode<<endl;
    for( Long_t iEst=0; iEst<fNEsts; iEst++){
        cout<<"N: "<<GetEstimator(iEst)->GetName()<<", def: "<<GetEstimator(iEst)->GetDefinition();
        cout<<", at: "<<GetEstimator(iEst)->GetValue()<<", Perc: "<<GetEstimator(iEst)->GetPercentile()<<", <> = "<<GetEstimator(iEst)->GetMean()<<", AnchPoint = "<<GetEstimator(iEst)->GetAnchorPoint()<<", AnchPerc = "<<GetEstimator(iEst)->GetAnchorPercentile()<<endl;
    }
}
//________________________________________________________________
Float_t AliMultSelection::GetMultiplicityPercentile(TString lName, Bool_t lEmbedEvSel)
{
    Float_t lReturnValue = AliMultSelectionCuts::kNoCalib;
    AliMultEstimator *lThis = GetEstimator(lName.Data());
    if( lThis ){
        lReturnValue = lThis->GetPercentile();
        //Bypass event selection if requested to do so 
        if ( fEvSelCode > 0 && lEmbedEvSel ) lReturnValue = fEvSelCode;
    }
    return lReturnValue;
}

//________________________________________________________________
Bool_t AliMultSelection::IsEventSelected()
//Function to determine if event has been selected: simple boolean output 
{
    Bool_t lReturnValue = kFALSE ; 
    if ( fEvSelCode == 0 ) lReturnValue = kTRUE; 
    return lReturnValue ;
}

//________________________________________________________________
void AliMultSelection::Evaluate( AliMultInput *lInput )
//Master function to evaluate all existing estimators based on
//a set of input variables. Error handling to be done with care...
{
    //Loop over estimators defined in the acquired list
    AliMultEstimator* estimator = 0;
    TIter             next(fEstimatorList);
    while ((estimator = static_cast<AliMultEstimator*>(next())))
        estimator->Evaluate(lInput);

//deprecated evaluation
#if 0
    AliMultEstimator *fThisEstimator = 0;
    for( Long_t iEst = 0; iEst<fNEsts; iEst++)
    {
        //Acquiring estimator from list
        fThisEstimator = GetEstimator(iEst);
        //cout<<"Estimator found: "<<fThisEstimator->GetName()<<", defined as "<<fThisEstimator->GetDefinition()<<endl;
        TString lEvaluateMe = fThisEstimator->GetDefinition();
        
        //FIXME: replace mean values if so desired
        
        //Looping over all defined inputs
        for(Long_t iVar = 0; iVar<lInput->GetNVariables(); iVar++){
            AliMultVariable *lVar1 = (AliMultVariable*) lInput->GetVariable(iVar);
            TString lReplacement = lVar1->GetName();
            lReplacement.Append(")");
            lReplacement.Prepend("(");
            if( !lVar1->IsInteger() ){
		lEvaluateMe.ReplaceAll( lReplacement, Form("%.10f",lVar1->GetValue() ) );
	    }else{
		lEvaluateMe.ReplaceAll( lReplacement, Form("%i",lVar1->GetValueInteger() ) );
	    }
        }
        //cout<<"String to evaluate: "<<lEvaluateMe<<endl;
        //FIXME: Error Handling needs to improve here
        TFormula lEvaluateMeForm("lEvaluateMeForm", lEvaluateMe );
        
        fThisEstimator -> SetValue ( lEvaluateMeForm.Eval(0) );
    }
#endif
}
//________________________________________________________________
void AliMultSelection::Setup(const AliMultInput* inp)
{
    AliMultEstimator* estimator = 0;
    TIter             next(fEstimatorList);
    
    while ((estimator = static_cast<AliMultEstimator*>(next())))
        estimator->SetupFormula(inp);
}
