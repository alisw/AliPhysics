#include "GammaConvUtilsTH1.h"
#include "TF1.h"
#include "TH1.h"

//_________________________________________________________________________________________________
utils_TH1::utils_TH1(std::string const &theId /*= "utils_TH1_defConstructor"*/) 
    :   id{theId},
        fTH1_ExponentialInterpolation_static_instance("utils_TH1")
{
    printf("INFO: utils_TH1::utils_TH1(std::string const &theId): created instance %s\n",
           id.data());
}

//_________________________________________________________________________________________________
utils_TH1::utils_TH1(utils_TH1 const &theRef)
:   id{theRef.id},
    fTH1_ExponentialInterpolation_static_instance{theRef.fTH1_ExponentialInterpolation_static_instance}
{
    printf("INFO: utils_TH1::utils_TH1(utils_TH1 const &theRef: created instance %s\n"
           "\tfrom &theRef = %p\n",
           id.data(),
           &theRef);
}


//_________________________________________________________________________________________________
TF1 *utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1(std::string const &theNewName, 
                                                               TH1         const &theTH1, 
                                                               bool               theIntegrate /* = false*/,
                                                               bool               theUseXtimesExp /* = false*/) // fits a function of the f(x) = x * exp([0] + [1]*x)
{
    printf("utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1(): called with theNewName: %s, theTH1: %s, theIntegrate = %d, theUseXtimesExp = %d\n",
           theNewName.data(), theTH1.GetName(), theIntegrate, theUseXtimesExp);
    
    TF1 *lTF1_result = 
        fTH1_ExponentialInterpolation_static_instance.InitializeWithHistoAndInsertInMapTF1(
            theTH1, 
            theIntegrate, 
            theUseXtimesExp);
    
    if (!lTF1_result){
        printf("FATAL: utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1():\n"
                "Could not initialize TH1_ExponentialInterpolation object with above"
                " parameters. Returning nullptr.\n");
        return nullptr;
    }
     
    printf("utils_TH1::InitGlobalPieceWiseExponentialInterpolationTF1(): %sturning %s%s %s.\n",
           lTF1_result ? "Re" : "Not re",
           lTF1_result ? lTF1_result->GetName() : "TF1* for histo",
           lTF1_result ? "." : theTH1.GetName(),
           lTF1_result ? "" : ". Returning nullptr.\n");
    return lTF1_result;
}
