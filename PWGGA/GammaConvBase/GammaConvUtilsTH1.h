#pragma once

class TF1;
class TH1;

#include <map>
#include <vector>
#include <string>

class utils_TH1
{ 
    /*
        The public interface of class utils_sstiefel/utils_TH1
    */
    public:
        utils_TH1();
        /* 
            get globally defined TF1 which is either fitted to the bin centers 
            OR such that the integrals in each in bin agree with their content. 
            
            if theUseXtimesExp=true: use f(x) = x * exp([0] + [1]*x). 
            else:                        f(x) =     exp([0] + [1]*x) 
        */
        // _________________________________________________________________________________________________
        static 
        TF1 *GlobalPieceWiseExponentialInterpolationTF1(std::string const &theNewName, 
                                                        TH1         const &theTH1, 
                                                        bool               theIntegrate = false,
                                                        bool               theUseXtimesExp = false); // x*f(x) instead f(x)


    /*
        In utils_TH1 there live 
            utils_TH1::TH1_ExponentialInterpolation_static:    the public utils_TH1::TH1_ExponentialInterpolation interface
            utils_TH1::TH1_ExponentialInterpolation:            a private helper class
    */
    
    class TH1_ExponentialInterpolation_static;    
    private:       
        // ===================== class utils_TH1::TH1_ExponentialInterpolation ===================================
        /*
        Main class for exponential interpolations of TH1 histograms
        It creates TF1 function objects that are driven from utils_TH1::TH1_ExponentialInterpolation instances,
        which need to survive the lifetime of the result TF1.
    
        data members:
          std::string   id;
          TH1&          fTH1;
          bool fIntegrate; 
          bool fUseXtimesExp;
      
          // get created during construction 
          TF1   &tf1_global;           
        */
        class TH1_ExponentialInterpolation
        {   
            public:
                TH1_ExponentialInterpolation() = delete;    
                TH1_ExponentialInterpolation(utils_TH1::TH1_ExponentialInterpolation const &) = delete;    

                TH1_ExponentialInterpolation(std::string const                   &_id,
                                             utils_TH1::TH1_ExponentialInterpolation_static &theParentRef,
                                             TH1 const                           &_th1,
                                             bool                                 _integrate,
                                             bool                                 _useXtimesExp,
                                             bool                                 _verbose = true);

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
                TF1 *GetLocalExponentialTF1(TH1    &theTH1, 
                                            double  theX, 
                                            bool    theIntegrate, 
                                            bool    theUseXtimesExp);
                
                double 
                    LocalExponentialInterpolate(TH1    &theTH1, 
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

                TF1 *ProduceNewNativeTF1(std::string const &theNewName);
                    
                ~TH1_ExponentialInterpolation();
        
            private:
                bool initGlobalFunctionObject(TF1 &theGlobalTF1, TH1  &theTH1);
                
                std::string                          id;
                TH1_ExponentialInterpolation_static &fStaticParent;
                std::string const                    fStaticParentId;
                TH1                                 &fTH1;
                
                bool                                 fIntegrate;    /* whether the fit minimizes: fIntegrate = false:  f(x)    -   h.ExpInter(h.GetBinCenter(h.FindBin(x))) true:   int(LeftEdge(i), RightEdge(i+1), f(x))  -   h.GetBinContent(i)) */
                bool                                 fUseXtimesExp; /* whether      false: th1.BinContent(i) -> th1.BinContents(i) * th1.BinCenter(i)      true: f(x) = exp([0] + [1]*x) -> x * f(x) */
        
                std::vector<TF1*>                    fVector_tf1_local;
                TF1                                 *fTF1_global;

    }; // end class utils_TH1::TH1_ExponentialInterpolation

    public:
// =================== utils_TH1::TH1_ExponentialInterpolation_static ======================================== 
    /*
        Static base class to maintain a map over all existing utils_TH1::TH1_ExponentialInterpolations
        
        data members:
            std::string id;
            std::map<TH1* const, utils_TH1::TH1_ExponentialInterpolation> fMap_utils_TH1::TH1_ExponentialInterpolation;
    */
    class TH1_ExponentialInterpolation_static {
        public:
            TH1_ExponentialInterpolation_static() = delete;
            TH1_ExponentialInterpolation_static(TH1_ExponentialInterpolation_static const &) = delete;
    
            TH1_ExponentialInterpolation_static(std::string const &_id,
                                                TH1 const         &_th1,
                                                bool               _integrate,
                                                bool               _useXtimesExp);
            
            // one of the two ways to create new TF1_globals with write access
            TF1 *CreateNewInterpolation(TH1 const  &_th1,
                                        bool        _integrate,
                                        bool        _useXtimesExp);
            
            // returns a valid TF1 interpolation from stack or create new one if necessary
            TF1 *GetInterpolationTF1(TH1 const  &theTH1,
                                     bool        theIntegrate,
                                     bool        theUseXtimesExp,
                                     bool        theCreateNewIfNecessary = true);
    
            TF1 const
                *GetInterpolationTF1_const(TH1 const &theTH1) const;
            
            std::string const 
                GetId() const                { return id; } 
            
            bool IsInitialized() const       { printf("in bool IsInitialized() const (): line 62\n"); 
                                               return fMap_TH1_ExponentialInterpolation.size(); }

        private:        
            // utils_TH1::TH1_ExponentialInterpolation_static data members
            std::string                                         id;
            std::map<TH1 const*, utils_TH1::TH1_ExponentialInterpolation*> fMap_TH1_ExponentialInterpolation;

    }; // end class utils_TH1::TH1_ExponentialInterpolation_static {
}; //  end class utils_TH1       
   //  
 
