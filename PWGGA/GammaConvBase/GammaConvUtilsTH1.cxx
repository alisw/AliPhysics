#include "GammaConvUtilsTH1.h"
#include "TF1.h"
#include "TH1.h"

//_________________________________________________________________________________________________
utils_TH1::utils_TH1(std::string const &theId /*= "utils_TH1_defConstructor"*/) 
    :   id{theId},
        fTH1_ExponentialInterpolation_static_instance_ptr{nullptr}
{
    printf("INFO: utils_TH1::utils_TH1(std::string const &theId): created instance %s\n",
           id.data());
}

//_________________________________________________________________________________________________
TF1 *utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1(std::string const &theNewName, 
                                                               TH1         const &theTH1, 
                                                               bool               theIntegrate /* = false*/,
                                                               bool               theUseXtimesExp /* = false*/) // fits a function of the f(x) = x * exp([0] + [1]*x)
{
    printf("utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1(): called with theNewName: %s, theTH1: %s, theIntegrate = %d, theUseXtimesExp = %d\n",
           theNewName.data(), theTH1.GetName(), theIntegrate, theUseXtimesExp);
    
    // create on heap so it survice longer
    // todo: check if that actually works.
    utils_TH1::TH1_ExponentialInterpolation_static *lStaticInstancePtr = 
        new utils_TH1::TH1_ExponentialInterpolation_static(theNewName, 
                                                           theTH1, 
                                                           theIntegrate, 
                                                           theUseXtimesExp);
    
        if (!lStaticInstancePtr){
            printf("FATAL: utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1():\n"
                    "Could not initialize TH1_ExponentialInterpolation object with above"
                    " parameters. Returning nullptr.\n");
            return nullptr;
        }
        
        if (!lStaticInstancePtr->IsInitialized()){
        printf("FATAL: utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1():\n"
                "Could not initialize TH1_ExponentialInterpolation object %s with provided"
                " parameters. Returning nullptr.\n",
               lStaticInstancePtr->GetId().data());
               return nullptr;
    }
            
    fTH1_ExponentialInterpolation_static_instance_ptr = lStaticInstancePtr;
    printf ("INFO: utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1(): instance: %s\n"
            "set fTH1_ExponentialInterpolation_static_instance_ptr to %p\n",
            id.data(),
            fTH1_ExponentialInterpolation_static_instance_ptr);
    
    TF1 *lResult = fTH1_ExponentialInterpolation_static_instance_ptr->GetInterpolationTF1(
        theTH1, 
        theIntegrate, 
        theUseXtimesExp);  
    
    printf("utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1(): %sturning %s%s %s.\n",
           lResult ? "Re" : "Not re",
           lResult ? lResult->GetName() : "TF1* for histo",
           lResult ? "." : theTH1.GetName(),
           lResult ? "" : ". Returning nullptr.\n");
    return lResult;
}
