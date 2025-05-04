#include "TH1_ExponentialInterpolation.h"


#include "TF1.h"
#include "TH1.h"
#include <iostream>


/*
Static base class to maintain a map over all existing TH1_ExponentialInterpolations
*/
TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(std::string const &_id,
                                                                         TH1 const &_th1,
                                                                         bool _integrate,
                                                                         bool _useXtimesExp)
    : id(Form("TH1_ExponentialInterpolation_static_%s_%s", _id.data(), _th1.GetName())),
      fMap_TH1_ExponentialInterpolation()
{
    printf("INFO: TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static(): instance %s:\n"
           "line 20\n",
           id.data());

    TH1_ExponentialInterpolation *lInstance = new TH1_ExponentialInterpolation(
        id,
        *this,
        _th1,
        _integrate,
        _useXtimesExp,
        true /* _verbose */);

    printf("INFO: TH1_ExponentialInterpolation_static::TH1_ExponentialInterpolation_static():\n"
           "\tCreated instance %s with options:\n"
           "\t\t parent: %p, th1: %s, integrate = %d, useXtimesExp = %d.\n",
           id.data(),
           this,
           _integrate,
           _useXtimesExp);

    if (!lInstance){
        printf("FATAL: TH1_ExponentialInterpolation_static():\n"
               "\tTH1_ExponentialInterpolation instance could not be created. Returning.\n");
               return;
    }

    auto const &lPair = fMap_TH1_ExponentialInterpolation.insert(
        { &_th1, lInstance }
    );

    fIsInitialized = lPair.second;
    printf("%s instance for histo %s%s.\n",
           fIsInitialized ? "Inserted" : "Did not insert",
           fIsInitialized ? ""         : " because it was nullptr",
           _th1.GetName());
}
 
//_________________________________________________________________________________________________
TF1 *TH1_ExponentialInterpolation_static::GetInterpolationTF1(TH1 const &theTH1,
                                                              bool       theIntegrate,
                                                              bool       theUseXtimesExp,
                                                              bool       theCreateNewIfNecessary /* = false*/ )
{
printf("INFO: TH1_ExponentialInterpolation_static::GetInterpolationTF1() called.\n"
           "\tparams: _th1: %s, _integrate = %d, _useXtimesExp = %d, theCreateNewIfNecessary = %d.\n",
            theTH1.GetName(), 
            theIntegrate, 
            theUseXtimesExp,
            theCreateNewIfNecessary);

    auto const   &lIt  = fMap_TH1_ExponentialInterpolation.find(&theTH1);
    bool found =  lIt != fMap_TH1_ExponentialInterpolation.end();

    TF1 *lResult = found 
        ?   lIt->second->GetTF1_global() 
        :   nullptr;  

    printf("INFO: TH1_ExponentialInterpolation_static::GetInterpolationTF1(): Found %s in map. %s%s.\n", 
            lResult ? lResult->GetName() : "no TF1",
           !lResult 
                ? theCreateNewIfNecessary 
                    ?   "Will create new one for histo "
                    :   "Returning nullptr. Set theCreateNewIfNecessary to true if you want creation of a new one.\n\n\n",
                : "");
    
    if (!lResult && theCreateNewIfNecessary){
        lResult = CreateNewInterpolation(theTH1, theIntegrate, theUseXtimesExp);  
    }
    
    if (!lResult){
        printf("\n\n\nFATAL: TH1_ExponentialInterpolation_static::GetInterpolationTF1(): instance: %s\n"
               "\tCould not (re)create a TF1 for histo %s for you. Apologies.\nReturning nullptr.\n\n\n",
               theTH1.GetName());
    }
    return lResult;
}


//_________________________________________________________________________________________________
TF1 const
    *TH1_ExponentialInterpolation_static::GetInterpolationTF1_const(TH1 const &theTH1) const
{
    auto const &lMap_constref 
        = static_cast<std::map<TH1 const*, TH1_ExponentialInterpolation*>>(
            fMap_TH1_ExponentialInterpolation);
    
    auto const lConstIt = lMap_constref.find(&theTH1);
    bool found = lConstIt != lMap_constref.cend();
    
    TH1_ExponentialInterpolation const *lTH1_expInter = found 
        ?   lConstIt->second 
        :    static_cast<TH1_ExponentialInterpolation*>(nullptr); 

    TF1 const *lTF1_global_const = lTH1_expInter 
        ?   lTH1_expInter->GetTF1_global_const()
        :   static_cast<TF1 const *>(nullptr);

    printf("%s\n", lTF1_global_const ? "" : Form("INFO: TH1_ExponentialInterpolation_static::GetInterpolationTF1_const(): instance id: %s:\n"
           "\treturning nullptr.\n", id.data()));

    return lTF1_global_const;
}

//_________________________________________________________________________________________________
// the only function that returns TF1s with write access
TF1 *TH1_ExponentialInterpolation_static::CreateNewInterpolation(TH1 const  &_th1,
                                                                 bool        _integrate,
                                                                 bool        _useXtimesExp)
{   
    printf("INFO: TH1_ExponentialInterpolation_static::CreateNewInterpolation() called.\n"
           "\tparams: _th1: %s, _integrate = %d, _useXtimesExp = %d\n",
            _th1.GetName(), 
            _integrate, 
            _useXtimesExp);

    TH1_ExponentialInterpolation *lExpInter = new TH1_ExponentialInterpolation(
         "",
        *this,
         _th1,
         _integrate,
         _useXtimesExp);

    printf("INFO: line 84.\n");
        
    // create new instance and take its global TF1
    auto &lPair_key_value = *new std::pair<TH1 const*, TH1_ExponentialInterpolation*>({ &_th1, lExpInter }); 

    printf("INFO: line 91\n");
 
    TF1 *lResult = lPair_key_value.second->GetTF1_global();
    if (!lResult){
        printf("WARNING: TH1_ExponentialInterpolation_static::CreateNewInterpolation(): id: %s\n"
               "\tNullptr found in map. This should never happen. Returning nullptr from this function.\n",
               id.data());
        return nullptr;
    }
printf("INFO: line 100\n");
    // todo check that this does not trigger copy constructor
    auto const &lPair_it_success = fMap_TH1_ExponentialInterpolation.insert(lPair_key_value);
printf("INFO: line 103\n");
    bool isNew = lPair_it_success.second;
    printf("%s: TH1_ExponentialInterpolation_static::CreateNewInterpolation(): id: %s\n"
          "\tReturning %s TF1: %s for histo: %s\n"
          "\t%s\n\n",
          isNew ? "INFO" : "\n\n\nWARNING", 
          id.data(), 
          isNew ? "newly inserted" : "from map", 
          lResult->GetName(), 
          _th1.GetName(),
          isNew ? "" : "This is maybe alarming since it means the function was called in wrong believe,"
                       "namely that an instance for the given fTH1 already existed in the map.\n");

    return lResult;
} // last TH1_ExponentialInterpolation_static methods

// ============================= class TH1_ExponentialInterpolation ===============================
// the normal constructor
//_________________________________________________________________________________________________
TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(std::string const                   &_id,
                                                           TH1_ExponentialInterpolation_static &theParentRef,
                                                           TH1 const                           &_th1,
                                                           bool                                 _integrate,
                                                           bool                                 _useXtimesExp,
                                                           bool                                 _verbose /* = false */)
:   id(_id.size() ? _id.data() : Form("TH1_ExponentialInterpolation_%s", _th1.GetName())),
    fStaticParent{theParentRef},
    fTH1{ * dynamic_cast<TH1 *>(_th1.Clone()) },
    fIntegrate{_integrate},
    fUseXtimesExp{_useXtimesExp},
    fVector_tf1_local{ std::vector<TF1*>( (1 + fTH1.GetNbinsX()), nullptr) }, // fill up 0th element so vectors indices will be same as histo bin numbers
    
    // std::vector<TF1*> v(n, nullptr);

    fTF1_global{nullptr}
{
    fTF1_global = ProduceNewNativeTF1("");
    bool isFullyInitialized = initGlobalFunctionObject(*fTF1_global, fTH1);
    
    if (_verbose){
        printf("TH1_ExponentialInterpolation::TH1_ExponentialInterpolation(): created instance with properties:\n\t"
                "id: %s, fTH1: %s, fIntegrate = %d, fUseXtimesExp = %d, fTF1_global: %s, isFullyInitialized = %d\n",
                id.data(),  
                fTH1.GetName(), 
                fIntegrate, 
                fUseXtimesExp, 
                fTF1_global 
                    ?   fTF1_global->GetName() 
                    :   "nullptr",
                isFullyInitialized);
    }
    
}

TH1_ExponentialInterpolation::~TH1_ExponentialInterpolation()
{

}

//_________________________________________________________________________________________________
TF1 *TH1_ExponentialInterpolation::GetLocalExponentialTF1(TH1   &theTH1, 
                                                          double theX, 
                                                          bool   theIntegrate, 
                                                          bool   theUseXtimesExp)
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
        ? lRange_edgeToEdge
        : lRange_centerToCenter;
    
    std::string lFunctionName(Form("%s_localExponential%s%s", 
                                   theTH1.GetName(), 
                                   theUseXtimesExp 
                                     ? "_*x" 
                                     : "",
                                   theIntegrate 
                                     ? "fitted_w/_int_cond" 
                                     : "calc_analyt_through_bin_centers"));
 
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
    else {
        double x1 = lRange_centerToCenter.first;
        double x2 = lRange_centerToCenter.second;
        double y1 = theTH1.GetBinContent(iLeftBin);
        double y2 = theTH1.GetBinContent(iRightBin);

        if (!(y1&&y2)) { 
            return nullptr;
        }
        // printf("using bin %d %f %f and bin %d %f %f\n", bin, x1, y1, adjBin, x2, y2);
        double b = TMath::Log(y2/y1) / (x2-x1);
        double a = TMath::Log(y1) - b*x1;
        lResult->SetParameters(a,b);
    }
    return lResult;
}

//_________________________________________________________________________________________________
double TH1_ExponentialInterpolation::LocalExponentialInterpolate(TH1   &theTH1, 
                                                                 double theX, 
                                                                 bool   theIntegrate /*= false*/,
                                                                 bool   theUseXtimesExp /*= false*/)
{
    TF1* f = GetLocalExponentialTF1(theTH1, theX, theIntegrate, theUseXtimesExp);
    return f ? f->Eval(theX) : 0.;
}


// _________________________________________________________________________________________________
// public
// returns true if fMap_local_tf1 is fully filled
    bool TH1_ExponentialInterpolation::initGlobalFunctionObject(TF1 &theGlobalTF1, TH1 &theTH1)
{
    printf("TH1_ExponentialInterpolation::initGlobalFunctionObject() instance id: %s\n\t"
           "Start initialization of a global function object based on TF1 %s and TH1 %s\n\t\t"
           "fIntegrate: %i, fUseXtimesExp: %i\n",
           id.data(), theGlobalTF1.GetName(), theTH1.GetName(), fIntegrate, fUseXtimesExp);

    size_t nBinsX = theTH1.GetNbinsX();
    size_t lNinsertions = 0;
    for (size_t iBin = 1; iBin <= nBinsX; ++iBin){
    
        double x = theTH1.GetBinCenter(iBin);
        theGlobalTF1.Eval(x); // this creates one local function per bin and stores it in fVector
        
        // TF1 *lTF1_candidate = fVector_tf1_local[iBin];
        TF1 *lTF1_candidate = (iBin <= fVector_tf1_local.size()) 
            ?   fVector_tf1_local[iBin]
            :   TH1_ExponentialInterpolation::GetLocalExponentialTF1(theTH1, x, fIntegrate, fUseXtimesExp);
        
            TF1 *lTF1_local_good = (lTF1_candidate && utils_TF1::IsInRange(*lTF1_candidate, x))
                ?   lTF1_candidate
                :   static_cast<TF1*>(nullptr);

        if (lTF1_local_good){
            // all good
            ++lNinsertions;
        } else {    
            // try here explicitly
            printf("WARNING: TH1_ExponentialInterpolation::initGlobalFunctionObject(): instance id: %s\n"
                   "\tFor some reason the automatic initialization through Eval() did not work!\n",
                   id.data());
            lTF1_local_good = GetLocalExponentialTF1(fTH1, 
                                                     x,
                                                     fIntegrate,
                                                     fUseXtimesExp);
            if (!lTF1_local_good){
                printf("WARNING: TH1_ExponentialInterpolation::initGlobalFunctionObject(): instance id: %s\n"
                       "\tThrough GetLocalExponentialTF1() did return nullptr. Returning nullptr.\n",
                       id.data());

                return false;
            }  
            printf("INFO: TF1 *TH1_ExponentialInterpolation::initGlobalFunctionObject():\n"
                   "\tiBin: %zu x = %f: Added %s to map.\n",
                   iBin, 
                   x, 
                   lTF1_local_good->GetName());
            fVector_tf1_local[iBin] = lTF1_local_good;
            ++lNinsertions;                        
        }
    }

    bool wasAlreadyFullyInitialized = !lNinsertions;
    if (lNinsertions){
        printf("TH1_ExponentialInterpolation::initGlobalFunctionObject() instance id: %s\n"
               "\tInserted new local TF1s for %zu out %zu bins. Instance was %s initialized.\n",
               id.data(), lNinsertions, nBinsX, wasAlreadyFullyInitialized
                   ?  "already fully"
                   :  "partly");   
    }
    bool wasFreshlyFullyInitialized = lNinsertions == nBinsX;
    bool isFullyInitializedNow = wasAlreadyFullyInitialized || wasFreshlyFullyInitialized;

    return isFullyInitializedNow;
}

//_________________________________________________________________________________________________
double TH1_ExponentialInterpolation::Evaluate(double *x, double *)
{
    auto hasTF1GoodRangeForTheX = [](TF1 const  &theTF1, 
                                      double     theX, 
                                      bool       thePrintInfoIfNot = false){
        bool lResult = (theX >= theTF1.GetXmin()) && (theX < theTF1.GetXmax());
        if (lResult){
            return true;
        }
        if (thePrintInfoIfNot){
            printf("INFO: TH1_ExponentialInterpolation::Evaluate(double *x, double *)::hasTF1GoodRangeForTheX():\n"
                   "\t\ttheTF1 = %s's range is %f - %f. That %scludes theX = %f\n",
            theTF1.GetName(), theTF1.GetXmin(), theTF1.GetXmax(), lResult ? "in" : "ex", theX);
        }
        return false;
    };
    
    // check if there is already a local interpol for this x
    size_t const lBin = static_cast<size_t>(fTH1.FindBin(*x));

    // try to get local tf1 from vector
    auto  *lTF1_local_good = (lBin <= fVector_tf1_local.size()) 
        ? fVector_tf1_local[lBin]
        : static_cast<TF1*>(nullptr);

    // (re)insert if necessary
    if (!lTF1_local_good){        
        lTF1_local_good = GetLocalExponentialTF1(fTH1, 
                                                *x,
                                                fIntegrate,
                                                fUseXtimesExp);  
             
        fVector_tf1_local[lBin] = lTF1_local_good;
        
        printf("%s: TH1_ExponentialInterpolation::Evaluate():"
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

TF1 *TH1_ExponentialInterpolation::GetTF1_global()
{
    return fTF1_global;
}
        

TF1 const
    *TH1_ExponentialInterpolation::GetTF1_global_const() const
{
    return static_cast<TF1 const *>(fTF1_global);
}

//_________________________________________________________________________________________________
TF1 *TH1_ExponentialInterpolation::ProduceNewNativeTF1(std::string const &theNewName)
{
    std::string lNewName(theNewName.size() 
        ? theNewName.data() 
        : Form("tf1_global_expInter_%s", id.data()));
    
    printf("TH1_ExponentialInterpolation::ProduceNewNativeTF1(): id: %s returning new TF1 with name '%s'.\n",
           id.data(), lNewName.data());

    TAxis const &lXaxis = *fTH1.GetXaxis();
    return new TF1(theNewName.data(),
                   this,
                   &TH1_ExponentialInterpolation::Evaluate,
                   lXaxis.GetXmin(), 
                   lXaxis.GetXmax(),
                   0, // nPar
                   "TH1_ExponentialInterpolation",
                   "Evaluate");
}


///////////////////// end TH1_ExponentialInterpolation ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
