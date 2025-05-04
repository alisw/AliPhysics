#include "utils_TH1.h"
#include "TH1_ExponentialInterpolation.h"


//_________________________________________________________________________________________________
TF1 *utils_TH1::GlobalPieceWiseExponentialInterpolationTF1(std::string const &theNewName, 
                                                           TH1         const &theTH1, 
                                                           bool               theIntegrate /* = false*/,
                                                           bool               theUseXtimesExp /* = false*/) // fits a function of the f(x) = x * exp([0] + [1]*x)
{
    printf("utils_TH1::GlobalPieceWiseExponentialInterpolationTF1(): called with theNewName: %s, theTH1: %s, theIntegrate = %d, theUseXtimesExp = %d\n",
           theNewName.data(), theTH1.GetName(), theIntegrate, theUseXtimesExp);
    
    printf("line14\n");
    TH1_ExponentialInterpolation_static *lInstance_ptr = 
        new TH1_ExponentialInterpolation_static(theNewName, 
                                                 theTH1, 
                                                 theIntegrate, 
                                                 theUseXtimesExp);
    printf("line19\n");

    return lInstance_ptr 
        ?   lInstance_ptr->GetInterpolationTF1(theTH1, theIntegrate, theUseXtimesExp)
        :   static_cast<TF1*>(nullptr);
}
