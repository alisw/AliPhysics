/**********************************************
 *
 * This is the basic estimator definition class
 *
 *  Estimator definitions have to contain 
 *  the following: 
 *
 *  --- Estimator Definition (logic)
 *  --- Estimator Value
 *  --- Estimator Mean Value
 *  --- Estimator Percentile
 *
 **********************************************/

#include "AliMultEstimator.h"
#include "TFolder.h"
#include "TObjString.h"
#include "TBrowser.h"

ClassImp(AliMultEstimator);

AliMultEstimator::AliMultEstimator() :
  TNamed(), fDefinition(""), fIsInteger(kFALSE), fValue(0), fMean(0), fPercentile(0)
{
  // Constructor
  
}
AliMultEstimator::AliMultEstimator(const char * name, const char * title, TString lInitDef):
TNamed(name,title), fDefinition(""), fIsInteger(kFALSE), fValue(0), fMean(0), fPercentile(0)
{
    //Named, titled, definition constructor
    fDefinition=lInitDef;
}

AliMultEstimator::~AliMultEstimator(){
  // destructor
  
}

Float_t AliMultEstimator::GetZ(){
    //Don't die for zero
    Double_t lReturnVal = 0;
    if ( TMath::Abs( fValue ) > 1e-6 ){
        lReturnVal = fValue/fMean;
    }
    return lReturnVal; 
}