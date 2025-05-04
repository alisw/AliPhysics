#pragma once
#include "TH1.h"


#include <string>



class utils_TH1
{
    public:

        utils_TH1();
        /* 
            get globally defined TF1 which is either fitted to the bin centers 
            OR such that the integrals in each in bin agree with their content. 
            
            if theUseXtimesExp=true: use f(x) = x * exp([0] + [1]*x). 
            else:                        f(x) =     exp([0] + [1]*x) 
        */
        // _________________________________________________________________________________________________
        static TF1 
            *GlobalPieceWiseExponentialInterpolationTF1(std::string const &theNewName, 
                                                        TH1         const &theTH1, 
                                                        bool               theIntegrate = false,
                                                        bool               theUseXtimesExp = false); // x*f(x) instead f(x)
};