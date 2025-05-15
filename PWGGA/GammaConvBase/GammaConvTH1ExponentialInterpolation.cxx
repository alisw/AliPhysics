#include "GammaConvUtilsTH1.h"
#include "TF1.h"
#include "TH1D.h"
#include "TString.h"
#include <iostream>

// ============================= class utils_TH1::TH1_ExponentialInterpolation ===============================
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
    // fill up 0th element so vectors indices will be same as histo bin numbers
    fTF1_global = produceNewNativeGlobalTF1("");
    printf("utils_TH1::TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(): Created instance: %s\n"
            "\tfTF1_global = %s.\n",
           id.data(),
           fTF1_global 
               ?   fTF1_global->GetName()
               :   "nullptr");
    
    // This fills the instances Vector with local interpolations
    bool isFullyInitialized = fTF1_global && initGlobalFunctionObject(*fTF1_global, fTH1);

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

    auto deletePointer = [](TObject *p){ if (p){ delete p; p = nullptr; } };
    
    deletePointer(fTF1_global);

    if (fVector_tf1_local.size()){
        for (auto *p : fVector_tf1_local){
            deletePointer(p);
        }
    }
}

//_________________________________________________________________________________________________
TF1 *utils_TH1::TH1_ExponentialInterpolation::GetLocalExponentialTF1(TH1    &theTH1, 
                                                                     double  theX, 
                                                                     bool    theIntegrate, 
                                                                     bool    theUseXtimesExp)
{
    TAxis const &lAxis = *theTH1.GetXaxis();    
    auto findFunctionI = [&](){
        int histoBinI = lAxis.FindBin(theX);
        return (theX >= lAxis.GetBinCenter(histoBinI)) 
        ? histoBinI 
        : std::max(histoBinI-1, 0);
    };
    int const iLeftBin = findFunctionI();
    int const iRightBin = iLeftBin+1;

    std::pair<double, double> lRange_edgeToEdge( 
        { lAxis.GetBinLowEdge(iLeftBin), lAxis.GetBinUpEdge(iRightBin) });    

    std::pair<double, double> lRange_centerToCenter( 
        { lAxis.GetBinCenter(iLeftBin), lAxis.GetBinCenter(iRightBin) });

    std::pair<double, double> const &lFunctionDefineRange = theIntegrate 
        ?   lRange_edgeToEdge
        :   lRange_centerToCenter;

    std::string lFunctionName(Form("%s_localExponential%s%s", 
                                   theTH1.GetName(), 
                                   theUseXtimesExp 
                                       ?    "_*x" 
                                       :    "",
                                   theIntegrate 
                                       ?    "fitted_w/_int_cond" 
                                       :    "calc_analyt_through_bin_centers"));

    TF1 *lResult = new TF1(lFunctionName.data(),
                           Form("%sexpo(0)", 
                           theUseXtimesExp ? "x*" : ""),  // = [x*] exp([0] + [1]*x) 
                           lFunctionDefineRange.first,
                           lFunctionDefineRange.second);

    std::string lFitOptions("QFMN0"); // Q = minimum printing
    if (theIntegrate){
        theTH1.Fit(lResult, 
                    lFitOptions.append(theIntegrate ? "I" : "").data(), 
                    "" /* global fit options I believe */, 
                    lFunctionDefineRange.first, 
                    lFunctionDefineRange.second);
        lResult->SetRange(lRange_centerToCenter.first, lRange_centerToCenter.second);
    } 
     // dont integrate, use bin contents
    else { // dont integrate, use bin contents
        double x1 = lRange_centerToCenter.first;
        double x2 = lRange_centerToCenter.second;
        double y1 = theTH1.GetBinContent(iLeftBin);
        double y2 = theTH1.GetBinContent(iRightBin);

        if (!(y1&&y2)) { 
            return nullptr;
        }
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation::GetLocalExponentialTF1(): instance %s\n"
                "\t using bins [ %d | %d |   with centers [ %f | %f |  and contents: [ %f | %f | \n", 
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
    } // end don't integrate, use bin contents
    return lResult;
}

//_________________________________________________________________________________________________
double utils_TH1::TH1_ExponentialInterpolation::EvaluateLocalExponentialInterpolate(
    TH1    &theTH1, 
    double  theX, 
    bool    theIntegrate /*= false*/,
    bool    theUseXtimesExp /*= false*/)
{
    TF1* f = GetLocalExponentialTF1(theTH1, theX, theIntegrate, theUseXtimesExp);
    return f ? f->Eval(theX) : 0.;
}

// _________________________________________________________________________________________________
// public
// returns true if fMap_local_tf1 is fully filled
bool utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject(TF1 &theGlobalTF1, TH1 &theTH1)
{
    printf("utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject() instance id: %s\n\t"
            "Start initialization of a global function object based on TF1 %s and TH1 %s\n\t\t"
            "fIntegrate: %i, fUseXtimesExp: %i\n",
           id.data(), 
           theGlobalTF1.GetName(), 
           theTH1.GetName(), 
           fIntegrate, 
           fUseXtimesExp);

    size_t nBinsX = theTH1.GetNbinsX();

    // enter 0 function for underflow bin
    fVector_tf1_local.push_back(new TF1("TF1bin0", 
                                        "0",  
                                        theTH1.GetXaxis()->GetBinLowEdge(1) 
                                            - theTH1.GetXaxis()->GetBinWidth(1),
                                        theTH1.GetXaxis()->GetBinLowEdge(1),
                                        0)
                                );
    
    size_t lNumberOfInsertions = 0;
    for (size_t iBin = 1; iBin <= nBinsX; ++iBin){
        double x = theTH1.GetBinCenter(iBin);
        theGlobalTF1.Eval(x); // this creates one local function per bin and stores it in fVector

        // TF1 *lTF1_candidate = fVector_tf1_local[iBin];
        TF1 *lTF1_candidate = (iBin <= fVector_tf1_local.size()) 
            ?   fVector_tf1_local[iBin]
            :   utils_TH1::TH1_ExponentialInterpolation::GetLocalExponentialTF1(theTH1, x, fIntegrate, fUseXtimesExp);

        TF1 *lTF1_local_good = (lTF1_candidate && TF1GoodForX(*lTF1_candidate, x))
            ?   lTF1_candidate
            :   static_cast<TF1*>(nullptr);
        
        if (lTF1_local_good){
            // all good
            ++lNumberOfInsertions;
        } else {    
            // try here explicitly
            printf("WARNING: utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject(): instance id: %s\n"
                    "\tFor some reason the automatic initialization through Eval() did not work!\n",
                    id.data());
            lTF1_local_good = GetLocalExponentialTF1(fTH1, 
                                                     x,
                                                     fIntegrate,
                                                     fUseXtimesExp);
            if (!lTF1_local_good){
                printf("WARNING: utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject(): instance id: %s\n"
                       "\tThrough GetLocalExponentialTF1() did return nullptr. Returning nullptr.\n",
                       id.data());
                return false;
            }

            printf("INFO: TF1 *utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject():\n"
                    "\tiBin: %zu x = %f: Added %s to map.\n",
                   iBin, 
                   x, 
                   lTF1_local_good->GetName());

            fVector_tf1_local[iBin] = lTF1_local_good;
            ++lNumberOfInsertions;                        
        }
    }

    // enter 0 function for overflow bin
    std::string lName(Form("TF1bin%zu", nBinsX+1));
    fVector_tf1_local.push_back(new TF1(lName.data(), 
                                        "0", 
                                        theTH1.GetXaxis()->GetBinLowEdge(1) 
                                            - theTH1.GetXaxis()->GetBinWidth(1),
                                        theTH1.GetXaxis()->GetBinLowEdge(1),
                                        0)
                                );


    fTF1_global = &theGlobalTF1;
    bool isFullyInitializedNow = lNumberOfInsertions == nBinsX;

    bool wasAlreadyFullyInitialized = !lNumberOfInsertions;
    if (lNumberOfInsertions){
        printf("utils_TH1::TH1_ExponentialInterpolation::initGlobalFunctionObject() instance id: %s\n"
                 "\tInserted new local TF1s for %zu out of %zu bins. Instance was %s initialized.\n",
               id.data(), 
               lNumberOfInsertions, 
               nBinsX, 
               wasAlreadyFullyInitialized
                   ?   "already fully"
                   :   "not");   
    }
    return isFullyInitializedNow;
}

//_________________________________________________________________________________________________
double utils_TH1::TH1_ExponentialInterpolation::Evaluate(double *x, double *)
{
    // check if there is already a local interpol for this x
    size_t const lBin = static_cast<size_t>(fTH1.FindBin(*x));

    // return 0 here if outside the range bins 1..nBinsX of fTH1
    if (!lBin || (lBin == (fTH1.GetNbinsX()+1))){
        printf("INFO: utils_TH1::TH1_ExponentialInterpolation::Evaluate(): instance %s\n"
                "\tcalled for x = %f, the %sflow bin of %s. Returning 0.\n",
               id.data(),
               *x,
               lBin ? "over" : "under",
               fTH1.GetName());
        return 0;
    }

    // try to get local tf1 from vector
    auto  *lTF1_local_good = (lBin <= fVector_tf1_local.size()) 
        ?   fVector_tf1_local[lBin]
        :   static_cast<TF1*>(nullptr);
    
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation::Evaluate(): id: %s\n"
            "\t x = %f, lBin = %zu, lTF1_local_good = %p, found existing lTF1_local_good ? %s.\n", 
           id.data(),
           *x,
           lBin,
           lTF1_local_good,
           lTF1_local_good ? "yes" : "no");

    // (re)insert if necessary
    if (!lTF1_local_good){   
        // this creates a new TF1 on HEAP!     
        lTF1_local_good = GetLocalExponentialTF1(fTH1, 
                                                 *x,
                                                 fIntegrate,
                                                 fUseXtimesExp);  

        fVector_tf1_local[lBin] = lTF1_local_good;
        printf("%s: utils_TH1::TH1_ExponentialInterpolation::Evaluate():"
                "\n\t\t%s local TF1 with name = %s"
                "\n\t\tfor bin %zu at x = %f into fVector_tf1_local.\n",
               lTF1_local_good ? "INFO" : "WARNING",
               lTF1_local_good ? "Inserted new" : "Could not store", 
               lTF1_local_good->GetName(), 
               lBin, 
               *x);
    } // done checking all conditions and pointers

    return lTF1_local_good ? lTF1_local_good->Eval(*x) : 0.;
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
TF1 *utils_TH1::TH1_ExponentialInterpolation::produceNewNativeGlobalTF1(std::string const &theNewName)
{
    std::string lNewName(theNewName.size() 
        ?    theNewName.data() 
        :    Form("tf1_global_expInter_%s", id.data()));

    printf("utils_TH1::TH1_ExponentialInterpolation::produceNewNativeGlobalTF1(): id: %s returning new TF1 with name '%s'.\n",
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
utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(std::string const &_id,
                                                                                    TH1         const &_th1,
                                                                                    bool               _integrate,
                                                                                    bool               _useXtimesExp)
    : id(Form("TH1_ExponentialInterpolation_static_%s_%s", _id.data(), _th1.GetName())),
      fMap_TH1_ExponentialInterpolation()
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(): instance %s:\n",
           id.data());

   utils_TH1::TH1_ExponentialInterpolation *lTH1_ExpInter_Instance_ptr = new utils_TH1::TH1_ExponentialInterpolation(
        id,
        *this,
        _th1,
        _integrate,
        _useXtimesExp,
        true /* _verbose */);

    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static():\n"
           "\tCreated instance %s with options:\n"
           "\t\t parent: %p, th1: %s, integrate = %d, useXtimesExp = %d.\n",
           id.data(),
           this,
           _th1.GetName(),
           _integrate,
           _useXtimesExp);

    if (!lTH1_ExpInter_Instance_ptr){
        printf("FATAL: utils_TH1::TH1_ExponentialInterpolation_static():\n"
               "\tTH1_ExponentialInterpolation instance could not be created. Returning.\n");
               return;
    }
    auto const &lPair = fMap_TH1_ExponentialInterpolation.insert(
        std::pair{ &_th1, lTH1_ExpInter_Instance_ptr }
    );

    TH1_ExponentialInterpolation &lTH1_ExpInter_Instance = *lTH1_ExpInter_Instance_ptr;
    if (!lPair.second){
        printf("FATAL: utils_TH1::TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(): id: %s\n"
                "\tCould not insert created lTH1_ExpInter_Instance into fMap_TH1_ExponentialInterpolation!\n"
                "\thisto name: %s. Returning.\n",
               id.data(),
               _th1.GetName());
        return;
    }

    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static(): instance: %s\n"
           "\tInserted lTH1_ExpInter_Instance into fMap_TH1_ExponentialInterpolation.\n"
           "\tInserted instance has properties:\n"
           "\tIsInitialized = %d, _th1: %s, fTF1_global = %p, TF1 name: %s\n"
           "Constructor done.\n",
           id.data(),
           lTH1_ExpInter_Instance.IsInitialized(),
           _th1.GetName(),
           lTH1_ExpInter_Instance.GetTF1_global_const(),
           lTH1_ExpInter_Instance.IsInitialized()
               ?    lTH1_ExpInter_Instance.GetTF1_global_const()->GetName()
               :    "nullptr");
}

//_________________________________________________________________________________________________
utils_TH1::TH1_ExponentialInterpolation_static::~TH1_ExponentialInterpolation_static()
{
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::~TH1_ExponentialInterpolation_static(): instance %s\n"
            "\tDestructor called.\n",
           id.data());
    
    for (auto &iPairRef : fMap_TH1_ExponentialInterpolation){
        if (iPairRef.second){
            delete iPairRef.second;
            iPairRef.second = nullptr;
        }
    }
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

    // create new instance and take its global TF1
    utils_TH1::TH1_ExponentialInterpolation *lExpInter_ptr = new utils_TH1::TH1_ExponentialInterpolation(
         "",
        *this,
         _th1,
         _integrate,
         _useXtimesExp);
    
    TF1 *lTF1_result_ptr = lExpInter_ptr
        ?    lExpInter_ptr->IsInitialized()
            ?    lExpInter_ptr->GetTF1_global()
            :    static_cast<TF1*>(nullptr) 
        :    static_cast<TF1*>(nullptr);

    if (!lTF1_result_ptr){
        printf("WARNING: utils_TH1::TH1_ExponentialInterpolation_static::createNew_TH1_ExponentialInterpolation(): id: %s\n"
               "\tNullptr found in map. This should never happen. Returning nullptr from this function.\n",
               id.data());
        return static_cast<TF1*>(nullptr);
    }

    // todo check that this does not trigger copy constructor
    auto const &lPair_it_success = fMap_TH1_ExponentialInterpolation.insert(std::pair{ &_th1,  lExpInter_ptr });
    bool isNew = lPair_it_success.second;
    printf("%s: utils_TH1::TH1_ExponentialInterpolation_static::createNew_TH1_ExponentialInterpolation(): id: %s\n"
            "\tReturning %s TF1: %s for histo: %s\n"
            "\t%s\n\n",
          isNew 
              ?    "INFO" 
              :    "\n\n\nWARNING", 
          id.data(), 
          isNew 
              ?    "newly inserted" 
              :    "from map", 
          lTF1_result_ptr->GetName(), 
          _th1.GetName(),
          isNew 
                ? "" 
                : "This is maybe alarming since it means the function was called in wrong believe,"
                   "namely that an instance for the given fTH1 already existed in the map.\n"
    );

    return lTF1_result_ptr;
} // last utils_TH1::TH1_ExponentialInterpolation_static methods

// creates new TF1 on heap if non existing and theCreateNewIfNecessary == true
//_________________________________________________________________________________________________
TF1 *utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1(
    TH1 const  &theTH1,
    bool       theIntegrate,
    bool       theUseXtimesExp,
    bool       theCreateNewIfNecessary /* = true*/ 
)
{
printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1() called.\n"
        "\tparams: _th1: %s, _integrate = %d, _useXtimesExp = %d, theCreateNewIfNecessary = %d.\n",
        theTH1.GetName(), 
        theIntegrate, 
        theUseXtimesExp,
        theCreateNewIfNecessary);

    auto lIt   = fMap_TH1_ExponentialInterpolation.find(&theTH1);
    bool found = (lIt != fMap_TH1_ExponentialInterpolation.end());

    utils_TH1::TH1_ExponentialInterpolation *lExpInter = found 
        ?    lIt->second
        :    static_cast<TH1_ExponentialInterpolation*>(nullptr);
    
    bool lExpInterIsInitialized = lExpInter && lExpInter->IsInitialized();
    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1(): id: %s\n"
            "\tfound = %d, lExpInter = %p, lExpInterIsInitialized = %d\n",
           id.data(),
           found,
           lExpInter,
           lExpInterIsInitialized);
    
    TF1 *lTF1result_ptr = lExpInterIsInitialized 
        ?   lExpInter->GetTF1_global() 
        :   theCreateNewIfNecessary
            ?   createNew_TH1_ExponentialInterpolation(theTH1, theIntegrate, theUseXtimesExp)
            :   nullptr;

    printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1():"
            " Found %s in map.\n"
            "\t%s%s.\n\n",
            found 
                ? lTF1result_ptr->GetName() 
                : "no TF1",
            !found
                ?  lTF1result_ptr
                    ?  "Created new one with name "
                    :  "Creation failed"
                : "",
            lTF1result_ptr
                ?  lTF1result_ptr->GetName()
                :  "");
    
    if (!lTF1result_ptr){
        printf("\n\n\nFATAL: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1(): instance: %s\n"
               "\tCould not (re)create a TF1 for histo %s for you. Apologies.\nReturning nullptr.\n\n\n",
               id.data(),
               theTH1.GetName());
        return nullptr;
    }
    
    auto const &lPair = fMap_TH1_ExponentialInterpolation.insert(std::pair{ &theTH1, lExpInter });
    if (!lPair.second){
        printf("WARNING/FATAL: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1(): instance: %s\n"
                "\tNewly created global TF1 %s was not inserted. That is very strange since it can only mean, it was already in the map.\n",
               id.data(),
               lTF1result_ptr->GetName());

        auto *lExpInter = lPair.first->second;
        if (lExpInter){
            if (lExpInter->IsInitialized()){
                lTF1result_ptr = lExpInter->GetTF1_global();
                return lTF1result_ptr;
            } else {
                printf("WARNING: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1(): instance: %s\n"
                        "\tElement that prevented insertion is not initialized. Initialize now",
                       id.data());
                
                

                // lTF1result_ptr = lExpInter->produceNewNativeGlobalTF1("GetInterpolationTF1()_cornercase");
                lTF1result_ptr = createNew_TH1_ExponentialInterpolation(theTH1,  theIntegrate, theUseXtimesExp );

                if (!lTF1result_ptr){
                    printf("FATAL: Could not create global TF1 from cornercase for instance %s.\n",
                        id.data());
                    return nullptr;
                }
                printf("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1(): instance: %s\n"
                        "\tCreation in cornercase was successfull:)\n",
                        id.data());
                return lTF1result_ptr;
            }
        } else {
            printf("FATAL: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1(): instance: %s\n"
                   "Nullptr in map. Erasing it.\n",
                   id.data());
            fMap_TH1_ExponentialInterpolation.erase(&theTH1);
        }
    }
    return lTF1result_ptr;
}

//_________________________________________________________________________________________________
TF1 const
    *utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1_const(TH1 const &theTH1) const
{
    auto const &lMap_constref 
        = static_cast<std::map<TH1 const*, utils_TH1::TH1_ExponentialInterpolation*>>(
            fMap_TH1_ExponentialInterpolation);
    
    auto const lConstIt = lMap_constref.find(&theTH1);
    bool found = lConstIt != lMap_constref.cend();
    
    utils_TH1::TH1_ExponentialInterpolation const *lTH1_expInter = found 
        ?   lConstIt->second 
        :    static_cast<TH1_ExponentialInterpolation*>(nullptr); 

    TF1 const *lTF1_global_const = lTH1_expInter 
        ?   lTH1_expInter->GetTF1_global_const()
        :   static_cast<TF1 const *>(nullptr);

    printf("%s\n", lTF1_global_const ? "" : Form("INFO: utils_TH1::TH1_ExponentialInterpolation_static::GetInterpolationTF1_const(): instance id: %s:\n"
           "\treturning nullptr.\n", id.data()));

    return lTF1_global_const;
} // last utils_TH1::TH1_ExponentialInterpolation_static methods

//_________________________________________________________________________________________________
bool utils_TH1::TH1_ExponentialInterpolation_static::IsInitialized() const
{ 
    size_t lSize = fMap_TH1_ExponentialInterpolation.size();
    printf("utils_TH1::IsInitialized(): IsInitialized ? %s . size = %zu\n", 
           lSize ? "yes" : "no", 
           lSize); 
    return static_cast<bool>(lSize); 
}
