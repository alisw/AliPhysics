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
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

ClassImp(AliMultSelection);

AliMultSelection::AliMultSelection() :
  TNamed(), fNEsts(0), fEstimatorList(0x0)
{
  // Constructor
    fEstimatorList = new TList();
}
AliMultSelection::AliMultSelection(const char * name, const char * title):
TNamed(name,title), fNEsts(0), fEstimatorList(0x0)
{
  // Constructor
    fEstimatorList = new TList();
    fEstimatorList->SetOwner(kTRUE);
}
AliMultSelection::AliMultSelection(AliMultSelection *lCopyMe){
    fNEsts = lCopyMe->GetNEstimators();
    
    fEstimatorList = new TList;
    fEstimatorList->SetOwner(kTRUE);
    
    for(Long_t iEst=0; iEst<lCopyMe->GetNEstimators(); iEst++){
        AliMultEstimator *lTemp = new AliMultEstimator( *lCopyMe->GetEstimator (iEst) );
        fEstimatorList->Add( lTemp );
    }
}
AliMultSelection::~AliMultSelection(){
    // destructor: clean stuff up 
    delete fEstimatorList;
    fEstimatorList=0x0;
}

void AliMultSelection::PrintInfo()
{
    cout<<"AliMultSelection Name..: "<<GetName()<<endl;
    cout<<"Estimators defined.....: "<<fNEsts<<endl;
    for( Long_t iEst=0; iEst<fNEsts; iEst++){
        cout<<"N: "<<GetEstimator(iEst)->GetName()<<", def: "<<GetEstimator(iEst)->GetDefinition();
        cout<<", at: "<<GetEstimator(iEst)->GetValue()<<", Perc: "<<GetEstimator(iEst)->GetPercentile()<<", <> = "<<GetEstimator(iEst)->GetMean()<<endl;
    }
}

Float_t AliMultSelection::GetMultiplicityPercentile(TString lName)
{
    Float_t lReturnValue = 200;
    AliMultEstimator *lThis = GetEstimator(lName.Data());
    if( lThis ) lReturnValue = lThis->GetPercentile(); 
    return lReturnValue;
}

void AliMultSelection::Evaluate( AliMultInput *lInput )
//Master function to evaluate all existing estimators based on
//a set of input variables. Error handling to be done with care...
{
    //Loop over estimators defined in the acquired list

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
}
