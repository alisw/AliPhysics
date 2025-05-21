#include "GammaConvUtilsTH1.h"
#include "TF1.h"
#include "TH1D.h"
#include "TString.h"
#include <iostream>

// ============================= class utils_TH1::TH1_ExponentialInterpolation ===============================
// copy constructor
//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(
    utils_TH1::TH1_ExponentialInterpolation const &theRef
)
:   id{theRef.id + "_cp-ass"},
    fStaticParent{theRef.fStaticParent},
    fStaticParentId{theRef.GetParentId()},    
    fTH1{theRef.fTH1},
    fIntegrate{theRef.fIntegrate},
    fUseXtimesExp{theRef.fUseXtimesExp},
    fVector_tf1_local{theRef.fVector_tf1_local},
    fTF1_global{theRef.fTF1_global}
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(\n"
             "\tutils_TH1::TH1_ExponentialInterpolation const &theRef): instance: %s\n"
             "\tCopy created instance from theRef = %p\n",
           id.data(),
           &theRef);

    // reinitialize fVector_tf1_local
    fTF1_global = createNewInitializedGlobalTF1(
        Form("%s_TF1_global_recreation_after_copy_assignment\n", id.data()));
   
    printf("%s: utils_TH1::TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(\n"
              "\tutils_TH1::TH1_ExponentialInterpolation const &theRef): instance %s\n."
              "\tRecreation of global TF1 from copied theRef %s\n",
           fTF1_global 
               ?    "INFO"
               :    "FATAL",
           id.data(),
           fTF1_global
               ? "was successfull."
               : "failed!");
}

// the normal constructor
//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(
    std::string const                               &_id,
    utils_TH1::TH1_ExponentialInterpolation_static  &_parentRef,
    TH1 const                                       &_th1,
    bool                                             _integrate,
    bool                                             _useXtimesExp,
    bool                                             _verbose /* = false */)
:   id(_id.size() ? _id.data() : Form("TH1_ExponentialInterpolation_%s", _th1.GetName())),
    fStaticParent(&_parentRef),
    fStaticParentId{_parentRef.GetId()},    
    fTH1{ *dynamic_cast<TH1 *>(_th1.Clone()) },
    fIntegrate{_integrate},
    fUseXtimesExp{_useXtimesExp},
    fVector_tf1_local{},
    fTF1_global{nullptr}
{
    fTF1_global = produceNewNativeUninitializedGlobalTF1("");
    printf("utils_TH1::TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(): Created instance: %s\n"
            "\tfTF1_global = %s.\n",
           id.data(),
           fTF1_global 
               ?   fTF1_global->GetName()
               :   "nullptr");
    
    // This fills the instances Vector with local interpolations
    fTF1_global = initGlobalFunctionObject(*fTF1_global, fTH1);
    bool isFullyInitialized = static_cast<bool>(fTF1_global);

    if (_verbose){
        printf("utils_TH1::TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(): "
               "instance %s has properties:\n"
               "\tid: %s, fStaticParent: %p, fStaticParentId: %s, fTH1: %s, fIntegrate = %d, fUseXtimesExp = %d, fTF1_global: %s, isFullyInitialized = %d\n",
               id.data(),
               id.data(),  
               &fStaticParent,
               fStaticParentId.data(),
               fTH1.GetName(), 
               fIntegrate, 
               fUseXtimesExp, 
               fTF1_global 
                  ?   fTF1_global->GetName() 
                  :   "nullptr",
               isFullyInitialized);
    }
}

//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation::~TH1_ExponentialInterpolation()
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation::~TH1_ExponentialInterpolation(): id: %s\n"
    "\tDestructor called.\n",
    id.data());
    
    if (fTF1_global){
        delete fTF1_global;
        fTF1_global = nullptr;
    }
}

// private
//_________________________________________________________________________________________________
TF1 *utils_TH1::TH1_ExponentialInterpolation::createNewInitializedGlobalTF1(
    std::string const &theNewName)
{
    TF1 *lResult = produceNewNativeUninitializedGlobalTF1(theNewName);
    if (!lResult){
        printf("FATAL: utils_TH1::TH1_ExponentialInterpolation::createNewInitializedGlobalTF1(): instance %s\n"
                "\tCall to produceNewNativeUninitializedGlobalTF1(std::string const &) retrieved nullptr.\n"
                "Returning nullptr.\n",
               id.data());
        return lResult;
    }
    return initGlobalFunctionObject(*lResult, fTH1);
}

//_________________________________________________________________________________________________
TF1 *utils_TH1::TH1_ExponentialInterpolation::GetNewLocalExponentialTF1(TH1    &theTH1, 
                                                                        double  theX, 
                                                                        bool    theIntegrate, 
                                                                        bool    theUseXtimesExp)
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation::GetNewLocalExponentialTF1(): id: %s\n"
            "\tCalled with theX = %f, theIntegrate = %d, theUseXtimesExp = %d\n",
           id.data(),
           theX,
           theIntegrate,
           theUseXtimesExp);

    // return zero function for under/overflow bins
    TAxis const &lAxis = *theTH1.GetXaxis();   
    int const iBinHisto = lAxis.FindBin(theX);

    int const iLeftBinMin = 1;
    int const iRightBinMax = lAxis.GetNbins();

    bool isEdgeLeft  = (iBinHisto == iLeftBinMin);
    bool isEdgeRight = (iBinHisto >= iRightBinMax);

    int const iLeftBin = isEdgeRight
        ?    iRightBinMax - 1                  // upper edge case
        :    std::max(iBinHisto, iLeftBinMin); // normal case
    int const iRightBin = iLeftBin + 1;

    double xmin = lAxis.GetBinLowEdge(iLeftBin);
    double xmax = lAxis.GetBinUpEdge(iRightBin);
    
    if (!iBinHisto || (iBinHisto == iRightBinMax + 1)){
        double width_n = lAxis.GetBinWidth(iBinHisto 
            ?    iRightBinMax 
            :    1);
        xmin = iBinHisto 
            ?    lAxis.GetBinUpEdge(iRightBinMax)
            :    lAxis.GetBinLowEdge(1) - width_n;
        xmax = iBinHisto
            ?    lAxis.GetBinUpEdge(iRightBinMax) + width_n
            :    lAxis.GetBinLowEdge(1);

        return ZeroFunctionTF1("", 
                               xmin,
                               xmax);
    }
    
    printf("lFound_in_map bins: iBinHisto = %d, iLeftBin = %d, iRightBin = %d\n",
           iBinHisto, iLeftBin, iRightBin);
    
    std::string lFunctionName(Form("TF1_%s_localExponential%s%s_bins_%d-%d-%d%s", 
                                   theTH1.GetName(),
                                   theUseXtimesExp 
                                       ?    "_*x" 
                                       :    "",
                                   theIntegrate 
                                       ?    "fitted_w/_int_cond" 
                                       :    "calc_analyt_through_bin_centers",
                                   iBinHisto, 
                                   iLeftBin,
                                   iRightBin,
                                   isEdgeLeft 
                                     ?  "_el"
                                     : isEdgeRight
                                       ? "_er"
                                       : ""));

    printf("Will create new TF1 with name = %s in range [ %f - %f |\n",
           lFunctionName.data(),
           xmin,
           xmax);
    
    TF1 *lResult = new TF1(lFunctionName.data(),
                           Form("%sexpo(0)", 
                                theUseXtimesExp ? "x*" : ""),  // = [x*] exp([0] + [1]*x) 
                           xmin,
                           xmax);
                           
    if (!lResult){
        // printf("no lResult, returning nullptr.\n");
        // return nullptr;
    }                  

    std::string lFitOptions("QFMN0"); // Q = minimum printing
    if (theIntegrate){
        theTH1.Fit(lResult, 
                    lFitOptions.append(theIntegrate ? "I" : "").data(), 
                    "" /* global fit options I believe */, 
                    xmin, 
                    xmax);
        lResult->SetRange(xmin, xmax);
    } 
    
    // dont integrate, use bin contents
    else { // dont integrate, use bin contents
        double centerLeft  = lAxis.GetBinCenter(iLeftBin);
        double centerRight = lAxis.GetBinCenter(iRightBin);
        double x1 = centerLeft;
        double x2 = centerRight;
        double y1 = theTH1.GetBinContent(iLeftBin);
        double y2 = theTH1.GetBinContent(iRightBin);

        if (!(y1&&y2)) {
            printf("INFO: utils_TH1::TH1_ExponentialInterpolation::GetNewLocalExponentialTF1(): instance %s\n"
                    "\tOne of y1 = %f and y2 %f are zero. Bin numbers: [ %d | %d |\n"
                    "\tCan only return zero function.\n",
                    id.data(),
                    y1,
                    y2,
                    iLeftBin,
                    iRightBin);

            lResult = new TF1(Form("%s_zero_function_bins_%d-%d",
                                    id.data(),
                                    iLeftBin,
                                    iRightBin),
                              "0",
                              xmin,
                              iRightBin);
        // success
        } else {
            printf("INFO: utils_TH1::TH1_ExponentialInterpolation::GetNewLocalExponentialTF1(): instance %s\n"
                    "\t using bins [ %d | %d |   with centers [ %f | %f |  and contents: [ %f | %f |. \n", 
                id.data(),
                iLeftBin, 
                iRightBin,
                x1, 
                x2,
                y1,
                y2);
            double b = TMath::Log(y2/y1) / (x2-x1);
            double a = TMath::Log(y1) - b*x1;
            lResult->SetParameters(a,b);
            lResult->SetRange(isEdgeLeft 
                                  ?    lAxis.GetBinLowEdge(iLeftBin) 
                                  :    lAxis.GetBinCenter(iLeftBin),
                              isEdgeRight 
                                  ?    lAxis.GetBinUpEdge(iRightBin)
                                  :    lAxis.GetBinCenter(iRightBin));
        }
    } // end don't integrate, use bin contents
    return lResult;
}

//_________________________________________________________________________________________________
double 
    utils_TH1::TH1_ExponentialInterpolation::EvaluateLocalExponentialInterpolate(
        TH1    &theTH1, 
        double  theX, 
        bool    theIntegrate /*= false*/,
        bool    theUseXtimesExp /*= false*/)
{
    TF1 *f = GetNewLocalExponentialTF1(theTH1, theX, theIntegrate, theUseXtimesExp);
    return f ? f->Eval(theX) : 0.;
}

// _________________________________________________________________________________________________
// private
// returns true if fMap_local_tf1 is fully filled, else nullptr
TF1* utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject(TF1 &theGlobalTF1, TH1 &theTH1)
{
    printf("utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject() instance id: %s\n\t"
            "Start initialization of a global function object based on TF1 %s and TH1 %s\n\t\t"
            "fIntegrate: %i, fUseXtimesExp: %i\n",
           id.data(), 
           theGlobalTF1.GetName(), 
           theTH1.GetName(), 
           fIntegrate, 
           fUseXtimesExp);

    size_t const nBinsX = theTH1.GetNbinsX();
    size_t const nBinsXplus2 = nBinsX + 2;

    // will initialize this expInters instance with this theGlobalTF1
    fVector_tf1_local.clear();
    for (size_t iBin = 0; iBin <= nBinsX + 1; ++iBin){
        double x = theTH1.GetBinCenter(iBin);

        // this creates one local function per bin and stores it in fVector
        double y = theGlobalTF1.Eval(x);
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject() instance id: %s\n"
                "\tiBin = %zu, nBinsX = %zu, x = %f, y = %f, size vector = %zu\n",
               id.data(),
               iBin,
               nBinsX,
               x, 
               y,
               fVector_tf1_local.size());
    }

    size_t const lVectorSize_after_loop = fVector_tf1_local.size();
    bool isGood = static_cast<int>(lVectorSize_after_loop) == nBinsXplus2; 
    
    printf("utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject() instance id: %s\n"
            "\tlVectorSize_after_loop = %zu, theTH1.GetNbinsX() = %d, that means it is %s fully initialized.\n"
            "\tReturning global function object with name: %s, and address = %p\n",
           id.data(), 
           lVectorSize_after_loop,
           theTH1.GetNbinsX(),
           isGood
              ?   "now"
              :   "still not",
           theGlobalTF1.GetName(),
           &theGlobalTF1);   

    return &theGlobalTF1;
}

//_________________________________________________________________________________________________
double 
    utils_TH1::TH1_ExponentialInterpolation::Evaluate(double *x, double *)
{
    double lResultValue = 0.;
    TF1 *lTF1_local_good = nullptr;
    bool wasObtainedFromCache = false;
    bool shouldInsertInside = false ;

    // check if there is already a local interpol for this x
    int const lBinHisto = fTH1.FindBin(*x);

    // not sure if this explicit cast is necessary
    int const lCurrentVectorSize = static_cast<int>(fVector_tf1_local.size());

    // return 0 here if outside the range bins 1..nBinsX of fTH1
    bool isUnderOverFlow = !lBinHisto || (lBinHisto == (fTH1.GetNbinsX() + 1)); 
    if (isUnderOverFlow){
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation::Evaluate(): instance %s\n"
                "\tcalled for x = %f, the %sflow bin of %s. This will later return 0.\n",
               id.data(),
               *x,
               lBinHisto ? "over" : "under",
               fTH1.GetName());
    }
    
    bool canInsertInside = lBinHisto < lCurrentVectorSize;
    bool isBinAlreadyInVector = !shouldInsertInside &&  lBinHisto < lCurrentVectorSize-1;
    bool canInsertAtBack = !isBinAlreadyInVector && (lCurrentVectorSize == lBinHisto);
    bool shouldProceed = isBinAlreadyInVector || canInsertInside || canInsertAtBack;
    
    // this means we can't even insert a zero function
    if (!shouldProceed){
        // uncomment for per Evaluate call print statement.
        /*
            printf("INFO: utils_TH1::TH1_ExponentialInterpolation::Evaluate(): instance %s\n"
                    "\tcalled for x = %f, lBinHisto = %d\n. This is outside the histo's %s range: [ %f | .. | %f |.\n"
                    "\tThere is no reason to perceed, we will return %f already here.\n",
                id.data(),
                *x,
                lBinHisto,
                fTH1.GetName(),
                fTH1.GetXaxis()->GetXmin(),
                fTH1.GetXaxis()->GetXmax(),
                lResultValue);
        */
        return lResultValue;
    }

    // this means there is alrady an existing local TF1.
    if (isBinAlreadyInVector ){
        lTF1_local_good = &fVector_tf1_local.at(static_cast<size_t>(lBinHisto));
        if (lTF1_local_good){
            shouldInsertInside = false;
            wasObtainedFromCache = true;
        }

    // if the vector length is such that it can be inserted at the back, do it
    } else if (shouldInsertInside || canInsertAtBack){

        // this should call copy constructor on the new TF1 returned by GetNewLocalExponentialTF1
        TF1 *lTF1_temp = GetNewLocalExponentialTF1(fTH1, 
                                                   *x, 
                                                   fIntegrate, 
                                                   fUseXtimesExp);
        if (lTF1_temp){
            auto lIt_vec_bin_i = fVector_tf1_local.begin() + lBinHisto; 
            if (lIt_vec_bin_i != fVector_tf1_local.end()){
                fVector_tf1_local.erase(lIt_vec_bin_i);
            }
            lTF1_local_good = &*fVector_tf1_local.emplace(lIt_vec_bin_i,  *lTF1_temp);
            if (lTF1_local_good != lTF1_temp){
                printf("INFO: Copy constructor was called as assumend. Deleting lTF1_temp.\n");
                delete lTF1_temp;
            } else {
                printf("INFO: Strange: Copy constructor was NOT called as assumend. NOT deleting lTF1_temp.\n");
            }           
        } else {
            printf("FATAL: utils_TH1::TH1_ExponentialInterpolation::Evaluate(): instance %s\n"
                    "\tGetNewLocalExponentialTF1() returned nullptr.\n"
                    "\tThis call of Evaluate at *x = %f will return 0. later.\n",
                   id.data(),
                   *x);
        }
    } else {
        // todo: I think this branch never gets called because it will exit already at the first possibility
        printf("// line 349: todo: I think this branch never gets called because it will exit already at the first possibility.\n"
                "\tINFO: utils_TH1::TH1_ExponentialInterpolation::Evaluate(): instance %s\n"
                "\t\tThe requested x = %f corresponds to the overflow bin number or even higher.\n"
                "\t\tfTH1.GetNbinsX() = %d, fVector_tf1_local.size() = %zu.\n"
                "\t\tWon't insert new function. This will return 0 later.\n",
              id.data(),
              *x,
              fTH1.GetNbinsX(),
              fVector_tf1_local.size());
    }

    if (lTF1_local_good){

        // uncomment for a printf statement for every function value returned
        /* 
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation::Evaluate(): id: %s\n"
                "\t x = %f, lBinHisto = %d, lTF1_local_good = %p, was obtained from cache ? %s.\n", 
            id.data(),
            *x,
            lBinHisto,
            lTF1_local_good,
            wasObtainedFromCache ? "yes" : "no"); 
        */
        lResultValue = lTF1_local_good->Eval(*x);
    }

    if (!lResultValue){
        // printf("INFO: utils_TH1::TH1_ExponentialInterpolation::Evaluate(): id: %s\n"
        //         "\t returning lResultValue = %f\n",
        //     id.data(),
        //     lResultValue);
    }
    return lResultValue;
}

//_________________________________________________________________________________________________
TF1 *utils_TH1::TH1_ExponentialInterpolation::GetTF1_global()
{
    return fTF1_global;
}

//_________________________________________________________________________________________________
TF1 const
*utils_TH1::TH1_ExponentialInterpolation::GetTF1_global_const() const
{
    return static_cast<TF1 const *>(fTF1_global);
}

// to be the future global TF1 
//_________________________________________________________________________________________________
TF1 *utils_TH1::TH1_ExponentialInterpolation::produceNewNativeUninitializedGlobalTF1(std::string const &theNewName)
{
    std::string lNewName(theNewName.size() 
        ?    theNewName.data() 
        :    Form("tf1_global_expInter_%s", id.data()));

    printf("utils_TH1::TH1_ExponentialInterpolation::produceNewNativeUninitializedGlobalTF1(): id: %s returning new TF1 with name '%s'.\n",
           id.data(), 
           lNewName.data());

    TAxis const &lXaxis = *fTH1.GetXaxis();
    return new TF1(lNewName.data(),
                   this,
                   &utils_TH1::TH1_ExponentialInterpolation::Evaluate,
                   //   == -1 for a histo h =  TH1D("h", "h", 1, 1, 2)
                   lXaxis.GetBinLowEdge(0) - lXaxis.GetBinWidth(1),  
                   lXaxis.GetBinUpEdge(lXaxis.GetNbins() + 1) 
                       + lXaxis.GetBinWidth(lXaxis.GetNbins()),                      //
                   0, // nPar
                   "TH1_ExponentialInterpolation",
                   "Evaluate");
}

//_________________________________________________________________________________________________
bool utils_TH1::TH1_ExponentialInterpolation::TF1GoodForX(TF1 const  &theTF1, 
                                                          double      theX, 
                                                          bool        thePrintInfoIfNot /* = false */) const
{
    bool lResult = (theX >= theTF1.GetXmin()) && (theX <= theTF1.GetXmax());
    if (lResult){
        return true;
    }

    if (thePrintInfoIfNot){
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation::Evaluate(double *x, double *)::hasTF1GoodRangeForTheX():\n"
                "\t\ttheTF1 = %s's range is %f - %f. That %scludes theX = %f\n",
                 theTF1.GetName(), theTF1.GetXmin(), theTF1.GetXmax(), lResult ? "in" : "ex", theX);
        }
    return false;
}

///////////////////// end utils_TH1::TH1_ExponentialInterpolation ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/*
Static base class to maintain a map over all existing utils_TH1::TH1_ExponentialInterpolations
*/
//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static()
:    id(Form("TH1_ExponentialInterpolation_static_default")),
     fMap_TH1_ExponentialInterpolation()
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::(): Default constructor called.\n"
            "\tCreated Instance. %s\n",
            id.data());
}

//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(
    std::string const &theIdSuffix
)
:    id(Form("TH1_ExponentialInterpolation_static_%s", theIdSuffix.data())),
     fMap_TH1_ExponentialInterpolation()
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static"
            "(std::string const &theIdSuffix): instance: %s\n"
            "\tDefault constructor called.\n",
            id.data());
}

//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(
    utils_TH1::TH1_ExponentialInterpolation_static const &theRef
)
:   
    id(theRef.id),
    fMap_TH1_ExponentialInterpolation(theRef.fMap_TH1_ExponentialInterpolation)    
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static(TH1_ExponentialInterpolation_static const &theRef): instance: %s\n"
            "\tCalled copy constructor. New Element has properties:\n"
            "\t id: %s, fMap_TH1_ExponentialInterpolation = %p\n",
            id.data(),
            id.data(),
            &fMap_TH1_ExponentialInterpolation);
}

// //_________________________________________________________________________________________________
// utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(
//     utils_TH1::TH1_ExponentialInterpolation_static const &&theOther
// )
// :   
//     id{std::move(theOther.id)},
//     fMap_TH1_ExponentialInterpolation{std::move(theOther.fMap_TH1_ExponentialInterpolation)}    
// {
//     printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static(TH1_ExponentialInterpolation_static const &&theOther): instance: %s\n"
//             "\tCalled move constructor. New Element has properties:\n"
//             "\tfMap_TH1_ExponentialInterpolation.size() = %zu\n",
//            id.data(),
//            fMap_TH1_ExponentialInterpolation.size());
// }

/* 
    The only function that fills the the base member fMap_TH1_ExponentialInterpolation.
    That means it creates new TH1_ExponentialInterpolation on heap and fills the pair <&TH1, TH1_ExponentialInterpolation> 
    into the map
*/

// checks first whether element is contained already
//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation&
    utils_TH1::TH1_ExponentialInterpolation_static::insertNewExpInterInstance(TH1 const  &_th1,
                                                                              bool        _integrate,
                                                                              bool        _useXtimesExp)
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(): instance %s:\n"
            "\tWill insert a new <&_TH1, TH1_ExponentialInterpolation> , into map fMap_TH1_ExponentialInterpolation if key"
            "\tnot already in map.\n",
           id.data());
    
    auto lIt_found = fMap_TH1_ExponentialInterpolation.find(&_th1);
    bool lFoundInMap = lIt_found != fMap_TH1_ExponentialInterpolation.end();
    if (lFoundInMap){
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::insertNewExpInterInstance(): instance %s"
                "Found element for TH1 %s in map. Will return an iterator to the element.\n",
                id.data(),
               _th1.GetName());
        return lIt_found->second;
    }
    
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::insertNewExpInterInstance(): id: %s\n"
            "\tAbout to emplace an instance in fMap_TH1_ExponentialInterpolation with target properties:\n"
            "\tparent = %p, th1: %s, integrate = %d, useXtimesExp = %d.\n",
           id.data(),
           this,
           _th1.GetName(),
           _integrate,
           _useXtimesExp);

    auto const &lPair_ref_emplace_success = fMap_TH1_ExponentialInterpolation.emplace(
            static_cast<TH1 const *>(&_th1), 
            utils_TH1::TH1_ExponentialInterpolation(
                id.data(),
                *this,
                _th1,
                _integrate,
                _useXtimesExp,
                true /* _verbose */)
    );

    if (!lPair_ref_emplace_success.second){
        printf("WARNING: utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(): id: %s\n"
                "\tInsertion did not take place, that means the element is already in the map and stays as is.\n"
                "Normally this should not happen here because it was just searched for before!\n"
                "\thisto name and address: %s %p .\n",
               id.data(),
               _th1.GetName(),
               &_th1);
    }
    auto &lIt_mapElement = lPair_ref_emplace_success.first;

    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static(): instance: %s\n"
             "\tInserted lTH1_ExpInter_Instance into fMap_TH1_ExponentialInterpolation.\n"
             "\tInserted instance has properties:\n"
             "\tIsInitialized = %d, _th1: %s, TF1 name: %s, fTF1_global = %p\n"
           "Constructor done.\n",
           id.data(),
           lIt_mapElement->second.IsInitialized(),
           _th1.GetName(),
           lIt_mapElement->second.GetTF1_global()->GetName(),
           lIt_mapElement->second.GetTF1_global());

    return lIt_mapElement->second;
}

//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation_static::~TH1_ExponentialInterpolation_static()
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::~TH1_ExponentialInterpolation_static(): instance %s\n"
            "\tDestructor called.\n",
           id.data());
}

// private:
//_________________________________________________________________________________________________
// the only function that returns TF1s with write access
TF1 *utils_TH1::TH1_ExponentialInterpolation_static::createNew_TH1_ExponentialInterpolation(
    TH1 const  &_th1,
    bool        _integrate,
    bool        _useXtimesExp)
{   
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::createNew_TH1_ExponentialInterpolation() called.\n"
           "\tparams: _th1: %s, _integrate = %d, _useXtimesExp = %d\n",
            _th1.GetName(), 
            _integrate, 
            _useXtimesExp);

    utils_TH1::TH1_ExponentialInterpolation &lExpInter_ref = 
        insertNewExpInterInstance(_th1, _integrate, _useXtimesExp);

    TF1 *lTF1_ptr = lExpInter_ref.GetTF1_global();
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::createNew_TH1_ExponentialInterpolation(): id: %s\n"
            "\tReturning global TF1: %s for histo: %s\n", 
          id.data(), 
          lTF1_ptr->GetName(), 
          _th1.GetName());

    return lTF1_ptr;
} // last utils_TH1::TH1_ExponentialInterpolation_static methods

//_________________________________________________________________________________________________
TF1 *utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1_opt_createNewIfNecessary(
    TH1 const  &theTH1,
    bool       theIntegrate,
    bool       theUseXtimesExp,
    bool       theCreateNewIfNecessary /* = true*/ )
{
printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1_opt_createNewIfNecessary() called.\n"
        "\tparams: _th1: %s, _integrate = %d, _useXtimesExp = %d, theCreateNewIfNecessary = %d.\n",
        theTH1.GetName(), 
        theIntegrate, 
        theUseXtimesExp,
        theCreateNewIfNecessary);

    // first check whether maps contains requested element
    auto lIt_map   = fMap_TH1_ExponentialInterpolation.find(&theTH1);
    bool lFound_in_map = (lIt_map != fMap_TH1_ExponentialInterpolation.end());

    utils_TH1::TH1_ExponentialInterpolation *lExpInter = lFound_in_map 
        ?    &lIt_map->second
        :     static_cast<TH1_ExponentialInterpolation*>(nullptr);
    
    if (!lExpInter){
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1_opt_createNewIfNecessary(): id: %s\n"
                "\tRetrived nullptr from fMap_TH1_ExponentialInterpolation! This should never happen.\n"
                "\tWill reinitialize if wanted.\n",
              id.data());
    }
    
    bool lExpInterIsInitialized = lExpInter && lExpInter->IsInitialized();
    bool lWillReinitialize = (!lExpInterIsInitialized && lExpInter) && 
                              theCreateNewIfNecessary;
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1_opt_createNewIfNecessary(): id: %s\n"
            "\tInfos: theCreateNewIfNecessary = %d, lFoundInMap = %d, lExpInterDefined = %d, lExpInter = %p, lExpInterIsInitialized = %d, lWillReinitialize = %d.\n",
           id.data(),
           theCreateNewIfNecessary,
           lFound_in_map,
           static_cast<bool>(lExpInter),
           lExpInter,
           lExpInterIsInitialized,
           lWillReinitialize);

    return lWillReinitialize
        ?   InitializeWithHistoAndInsertInMapTF1(theTH1, theIntegrate, theUseXtimesExp)
        :   static_cast<TF1*>(nullptr);
}

//_________________________________________________________________________________________________
TF1 const
    *utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1_const(TH1 const &theTH1) const
{
    auto const &lMap_constref = 
        static_cast<std::map<TH1 const*, 
                    utils_TH1::TH1_ExponentialInterpolation> const>(fMap_TH1_ExponentialInterpolation);
    
    auto const lConstIt = lMap_constref.find(&theTH1);
    bool lFound_in_map = lConstIt != lMap_constref.cend();
    
    utils_TH1::TH1_ExponentialInterpolation const *lTH1_expInter = lFound_in_map 
        ?    static_cast<utils_TH1::TH1_ExponentialInterpolation const*>(&lConstIt->second) 
        :    static_cast<utils_TH1::TH1_ExponentialInterpolation const*>( nullptr); 

    TF1 const *lTF1_global_const = lTH1_expInter && lTH1_expInter->IsInitialized() 
        ?   lTH1_expInter->GetTF1_global_const()
        :   static_cast<TF1 const *>(nullptr);

    if (!lTF1_global_const){
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1_const(): instance id: %s:\n"
                "\tReturning nullptr.\n", 
             id.data());
    }
    return lTF1_global_const;
} // last utils_TH1::TH1_ExponentialInterpolation_static methods

//_________________________________________________________________________________________________
TF1 *utils_TH1::TH1_ExponentialInterpolation_static::InitializeWithHistoAndInsertInMapTF1(
    TH1 const  &theTH1, 
    bool        theIntegrate, 
    bool        theUseXtimesExp)
{
    return insertNewExpInterInstance(theTH1, 
                                     theIntegrate, 
                                     theUseXtimesExp).GetTF1_global();
}

//_________________________________________________________________________________________________
bool utils_TH1::TH1_ExponentialInterpolation_static::IsInitialized() const
{ 
    size_t const lSize = fMap_TH1_ExponentialInterpolation.size();
    bool lResult = static_cast<bool>(lSize);
    printf("utils_TH1::IsInitialized(): IsInitialized ? %s . size = %zu\n", 
           lResult 
               ?    "yes" 
               :    "no", 
           lSize); 
    return lResult; 
}
