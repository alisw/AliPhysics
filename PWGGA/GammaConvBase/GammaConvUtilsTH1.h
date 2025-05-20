#pragma once

class TF1;
class TH1;

#include <map>
#include <vector>
#include <string>

/*
    In utils_TH1 there live 
    utils_TH1::TH1_ExponentialInterpolation_static:     a private helper class 
    utils_TH1::TH1_ExponentialInterpolation:            another private helper class 
*/
    
class utils_TH1
{ 
    private:
    // ===================== class utils_TH1::TH1_ExponentialInterpolation ===================================
    /*
        Main class for exponential interpolations of TH1 histograms
        It creates TF1 function objects that are driven from utils_TH1::TH1_ExponentialInterpolation instances,
        which need to survive the lifetime of the result TF1.    
    */
        class TH1_ExponentialInterpolation_static;

        class TH1_ExponentialInterpolation
        {   
            public:
                TH1_ExponentialInterpolation() = delete;    
                TH1_ExponentialInterpolation(utils_TH1::TH1_ExponentialInterpolation const &theRef);

                TH1_ExponentialInterpolation(utils_TH1::TH1_ExponentialInterpolation const &&theRef) = delete;

                TH1_ExponentialInterpolation(std::string const                               &_id,
                                             utils_TH1::TH1_ExponentialInterpolation_static  &_parentRef,
                                             TH1 const                                       &_th1,
                                             bool                                             _integrate,
                                             bool                                             _useXtimesExp,
                                             bool                                             _verbose = true);
                
                ~TH1_ExponentialInterpolation();

                double 
                    Evaluate(double *x, double *);
                    
                bool TF1GoodForX(TF1 const  &theTF1, 
                                double      theX, 
                                bool        thePrintInfoIfNot = false) const;
            
                // get a TF1 exponential defined by:
                /* 
                case theIntegrate = false: 
                    - the next two neares bin centers to theX
                    - the two coefficients are calculated analytically from the two bin content values
                case theIntegrate = true: 
                    - the lower and upper edges of the two nearest bins to theX 
                    - the two coefficients are determined from a fit using the two bins  
                */
                TF1 *GetNewLocalExponentialTF1(TH1    &theTH1, 
                                            double  theX, 
                                            bool    theIntegrate, 
                                            bool    theUseXtimesExp);
                
                double 
                    EvaluateLocalExponentialInterpolate(TH1    &theTH1, 
                                                        double  theX, 
                                                        bool    theIntegrate = false,
                                                        bool    theUseXtimesExp = false);

                // getters 
                std::string const
                    &GetParentId() const         { return fStaticParentId; }

                std::string const
                    &GetId() const               { return id; } 

                TH1 const 
                    &GetTH1() const              { return fTH1; }
                
                // the only write accessor to fTF1_global
                TF1 *GetTF1_global();
                
                TF1 const
                    *GetTF1_global_const() const;
                
                bool IsInitialized() const { return fTF1_global; }
                    
                utils_TH1::TH1_ExponentialInterpolation 
                    &operator=(utils_TH1::TH1_ExponentialInterpolation &theRef) = delete;

                utils_TH1::TH1_ExponentialInterpolation 
                    &operator=(utils_TH1::TH1_ExponentialInterpolation const &theConstRef) = delete;
            
            private:
                
                // fill the vector of an instance returned by produceNewNativeUninitializedGlobalTF1
                TF1 *initGlobalFunctionObject(TF1 &theGlobalTF1, TH1  &theTH1);
                
                // produce a new TF1 with full definedness except its fVector_tf1_local is empty
                TF1 *produceNewNativeUninitializedGlobalTF1(std::string const &theNewName);
                
                // this calls initGlobalFunctionObject() on produceNewNativeUninitializedGlobalTF1()
                TF1 *createNewInitializedGlobalTF1(std::string const &theNewName);
                                
                std::string                          id;
                TH1_ExponentialInterpolation_static *fStaticParent;
                std::string const                    fStaticParentId;
                TH1                                 &fTH1;
                
                bool                                 fIntegrate;    /* whether the fit minimizes: fIntegrate = false:  f(x)    -   h.ExpInter(h.GetBinCenter(h.FindBin(x))) true:   int(LeftEdge(i), RightEdge(i+1), f(x))  -   h.GetBinContent(i)) */
                bool                                 fUseXtimesExp; /* whether      false: th1.BinContent(i) -> th1.BinContents(i) * th1.BinCenter(i)      true: f(x) = exp([0] + [1]*x) -> x * f(x) */
        
                // holds for every bin in fTH1 a TF1 pointer to the local interpolations
                std::vector<TF1>                    fVector_tf1_local; // lifetime must exceed lifetime of fTF1_global
                TF1                                *fTF1_global; 
       }; // end class utils_TH1::TH1_ExponentialInterpolation

    
    // =================== utils_TH1::TH1_ExponentialInterpolation_static ======================================== 
        
        //  Static base class to maintain a map over all existing utils_TH1::TH1_ExponentialInterpolations
        class TH1_ExponentialInterpolation_static {
            public:
                TH1_ExponentialInterpolation_static();
                
                TH1_ExponentialInterpolation_static(TH1_ExponentialInterpolation_static const &theRef);
                
                TH1_ExponentialInterpolation_static(TH1_ExponentialInterpolation_static const &&theRef);
                
                TH1_ExponentialInterpolation_static(std::string const &theIdSuffix);
                
                // the only function that creates new TH1_ExponentialInterpolation on heap and stores its pointer in a member map.
                TH1_ExponentialInterpolation_static(std::string const &_id,
                                                    TH1 const         &_th1,
                                                    bool               _integrate,
                                                    bool               _useXtimesExp);
                
                ~TH1_ExponentialInterpolation_static();

                // returns a valid TF1 interpolation from stack or create new one if necessary
                TF1 *GetInterpolationTF1_opt_createNewIfNecessary(TH1 const  &theTH1,
                                                                  bool       theIntegrate,
                                                                  bool       theUseXtimesExp,
                                                                  bool       theCreateNewIfNecessary = true );        
                TF1 const
                    *GetInterpolationTF1_const(TH1 const &theTH1) const;
                
                std::string const 
                    GetId() const                { return id; } 
                
                TF1 *InitializeWithHistoAndInsertInMapTF1(TH1 const  &theTH1, 
                                                          bool        theIntegrate, 
                                                          bool        theUseXtimesExp);
                
                bool IsInitialized() const;
        
            private:
        
                // one of the two ways to create new TF1_globals with write access
                // creates new TF1 on heap and inserts into map
                //_________________________________________________________________________________________________
                TF1 *createNew_TH1_ExponentialInterpolation(TH1 const  &_th1,
                                                            bool        _integrate,
                                                            bool        _useXtimesExp);
                
                // looks for existence of an element <_th1, TH1_ExponentialInterpolation*> 
                // if (found): returns fGlobal_TF1 from found element
                // else: inserts newly created
                // return value: in all cases a TF1* except creation fails
                //_________________________________________________________________________________________________
                TH1_ExponentialInterpolation
                    &insertNewExpInterInstance(TH1   const &_th1,
                                               bool         _integrate,
                                               bool         _useXtimesExp);
                
                // utils_TH1::TH1_ExponentialInterpolation_static data members
                std::string                                                    id;
                
                std::map<TH1 const*, utils_TH1::TH1_ExponentialInterpolation>  fMap_TH1_ExponentialInterpolation;
        }; // end class utils_TH1::TH1_ExponentialInterpolation_static {
    
    public:
    // ===================== class utils_TH1 ================================================================
        /*
            The public interface of class utils_sstiefel/utils_TH1
        */
        utils_TH1(std::string const &theId = "utils_TH1_defConstructor");

        utils_TH1(utils_TH1 const &theRef);  

        ~utils_TH1();  
        
        /* 
            get globally defined TF1 which is either fitted to the bin centers 
            OR such that the integrals in each in bin agree with their content. 
            
            if theUseXtimesExp=true: use f(x) = x * exp([0] + [1]*x). 
            else:                        f(x) =     exp([0] + [1]*x) 
        */
        // _________________________________________________________________________________________________ 
        TF1 *InitGlobalPieceWiseExponentialInterpolationTF1(std::string const &theNewName, 
                                                            TH1         const &theTH1, 
                                                            bool               theIntegrate = false,
                                                            bool               theUseXtimesExp = false); // x*f(x) instead f(x)
        
        // _________________________________________________________________________________________________ 
        bool IsGlobalPieceWiseExponentialInterpolationInitialized() const  
        {   
            return static_cast<bool>(fTH1_ExponentialInterpolation_static_instance.IsInitialized());
        }

        // _________________________________________________________________________________________________ 
        static TF1* 
        ZeroFunctionTF1(std::string theNewName = "", 
                        double      theXmin    = 0., 
                        double      theXmax    = 1.);

    // namespace utils_TH1::
    private:       
        std::string                          id;

        // this can hold all exponential interpolations for this instance of utils_TH1
        TH1_ExponentialInterpolation_static  fTH1_ExponentialInterpolation_static_instance;  

}; //  end class utils_TH1
 