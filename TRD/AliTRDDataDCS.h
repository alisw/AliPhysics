#ifndef ALITRDDATADCS_H
#define ALITRDDATADCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Extracts the DCS information                                           //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class TObjArray;
class TString;

class TGraph;
class AliSplineFit;
class TMap;

class AliTRDDataDCS : public TNamed
{
  
 public :
    
  AliTRDDataDCS ();
  ~AliTRDDataDCS ();
  
  Bool_t ExtractDCS (TMap * aliDCS);
  Bool_t PerformFit ();
  void   ClearFits ();
  void   ClearGraphs ();
  
  UInt_t 	GetNAlias () const {return fNAlias;}
  TString       GetAmandaStr (UInt_t iAlias) const;
  UInt_t 	GetNChannel (UInt_t iAlias) const; 
  
  TGraph * 		GetGraph (UInt_t iAlias, UInt_t iChannel = 0) const;
  AliSplineFit * 	GetFit (UInt_t iAlias, UInt_t iChannel = 0) const;
  
  // Get TGraph
  TGraph * GetGraphChamberByteStatus (UInt_t iSensor) const	        {return GetGraph (kChamberByteStatus, iSensor);}
  TGraph * GetGraphPreTrigger () const 					{return GetGraph (kPreTrigger);}
  TGraph * GetGraphGoofyHv () const			        	{return GetGraph (kGoofyHv);}
  TGraph * GetGraphGoofyPeakPos (UInt_t iSensor) const	                {return GetGraph (kGoofyPeakPos, iSensor);}
  TGraph * GetGraphGoofyPeakArea (UInt_t iSensor) const 	        {return GetGraph (kGoofyPeakArea, iSensor);}
  TGraph * GetGraphGoofyTemp (UInt_t iSensor) const		        {return GetGraph (kGoofyTemp, iSensor);}
  TGraph * GetGraphGoofyPressure () const 				{return GetGraph (kGoofyPressure);}
  TGraph * GetGraphGoofyVelocity () const 				{return GetGraph (kGoofyVelocity);}
  TGraph * GetGraphGoofyGain (UInt_t iSensor) const		        {return GetGraph (kGoofyGain, iSensor);}
  TGraph * GetGraphGoofyCO2 ()  const					{return GetGraph (kGoofyCO2);}
  TGraph * GetGraphGoofyN2 () const 					{return GetGraph (kGoofyN2);}
  TGraph * GetGraphGasO2 () const 					{return GetGraph (kGasO2);}
  TGraph * GetGraphGasOverpressure () const 				{return GetGraph (kGasOverpressure);}
  TGraph * GetGraphEnvTemp (UInt_t iSensor) const			{return GetGraph (kEnvTemp, iSensor);}
  TGraph * GetGraphHvAnodeImon (UInt_t iSensor) const		        {return GetGraph (kHvAnodeImon, iSensor);}
  TGraph * GetGraphHvDriftImon (UInt_t iSensor) const		        {return GetGraph (kHvDriftImon, iSensor);}
  TGraph * GetGraphHvAnodeUmon (UInt_t iSensor) const		        {return GetGraph (kHvAnodeUmon, iSensor);}
  TGraph * GetGraphHvDriftUmon (UInt_t iSensor) const		        {return GetGraph (kHvDriftUmon, iSensor);}
  TGraph * GetGraphAdcClkPhase () const 				{return GetGraph (kAdcClkPhase);}
  TGraph * GetGraphAtmPressure () const 				{return GetGraph (kAtmPressure);}
  TGraph * GetGraphLuminosity () const 					{return GetGraph (kLuminosity);}
  TGraph * GetGraphMagneticField () const				{return GetGraph (kMagneticField);}
  
  AliSplineFit * GetFitChamberByteStatus (UInt_t iSensor) const	        {return GetFit (kChamberByteStatus, iSensor);}
  AliSplineFit * GetFitPreTrigger () const 				{return GetFit (kPreTrigger);}
  AliSplineFit * GetFitGoofyHv () const					{return GetFit (kGoofyHv);}
  AliSplineFit * GetFitGoofyPeakPos (UInt_t iSensor) const	        {return GetFit (kGoofyPeakPos, iSensor);}
  AliSplineFit * GetFitGoofyPeakArea (UInt_t iSensor) const 	        {return GetFit (kGoofyPeakArea, iSensor);}
  AliSplineFit * GetFitGoofyTemp (UInt_t iSensor) const		        {return GetFit (kGoofyTemp, iSensor);}
  AliSplineFit * GetFitGoofyPressure () const 				{return GetFit (kGoofyPressure);}
  AliSplineFit * GetFitGoofyVelocity () const 				{return GetFit (kGoofyVelocity);}
  AliSplineFit * GetFitGoofyGain (UInt_t iSensor) const		        {return GetFit (kGoofyGain, iSensor);}
  AliSplineFit * GetFitGoofyCO2 ()  const				{return GetFit (kGoofyCO2);}
  AliSplineFit * GetFitGoofyN2 () const 				{return GetFit (kGoofyN2);}
  AliSplineFit * GetFitGasO2 () const 					{return GetFit (kGasO2);}
  AliSplineFit * GetFitGasOverpressure () const 			{return GetFit (kGasOverpressure);}
  AliSplineFit * GetFitEnvTemp (UInt_t iSensor) const			{return GetFit (kEnvTemp, iSensor);}
  AliSplineFit * GetFitHvAnodeImon (UInt_t iSensor) const		{return GetFit (kHvAnodeImon, iSensor);}
  AliSplineFit * GetFitHvDriftImon (UInt_t iSensor) const		{return GetFit (kHvDriftImon, iSensor);}
  AliSplineFit * GetFitHvAnodeUmon (UInt_t iSensor) const		{return GetFit (kHvAnodeUmon, iSensor);}
  AliSplineFit * GetFitHvDriftUmon (UInt_t iSensor) const		{return GetFit (kHvDriftUmon, iSensor);}
  AliSplineFit * GetFitAdcClkPhase () const 				{return GetFit (kAdcClkPhase);}
  AliSplineFit * GetFitAtmPressure () const 				{return GetFit (kAtmPressure);}
  AliSplineFit * GetFitLuminosity () const 				{return GetFit (kLuminosity);}
  AliSplineFit * GetFitMagneticField () const				{return GetFit (kMagneticField);}
  
  void Print (const Option_t * const option = "") const;
    
 protected :
          
    TGraph       * FindAndMakeGraph (TMap * const dcsMap
                                   , const char * amandaStr
				   , char dataType);
    AliSplineFit * Fit (const TGraph * const graph,
			Int_t  kMinPoints, Int_t  kIter, 
			Double_t  kMaxDelta, Int_t  kFitReq);
    
    void Init ();
    void InitFits ();
    void InitGraphs ();
    
    void SetConf (UInt_t iAlias, const char * amanda, char dataType, UInt_t nChannel, 
		  Bool_t enableGraph, Bool_t enableFit, Int_t kMinPoints, 
		  Int_t kIter, Double_t kMaxDelta, Int_t kFitReq);
    
 private :
    
    enum { kChamberByteStatus = 0
         , kPreTrigger        = 1
         , kGoofyHv           = 2
         , kGoofyPeakPos      = 3
         , kGoofyPeakArea     = 4
         , kGoofyTemp         = 5
         , kGoofyPressure     = 6
         , kGoofyVelocity     = 7
         , kGoofyGain         = 8
         , kGoofyCO2          = 9
         , kGoofyN2           = 10
         , kGasO2             = 11
         , kGasOverpressure   = 12
         , kEnvTemp           = 13
         , kHvAnodeImon       = 14
         , kHvDriftImon       = 15
         , kHvAnodeUmon       = 16
         , kHvDriftUmon       = 17
         , kAdcClkPhase       = 18
         , kAtmPressure       = 19
         , kLuminosity        = 20
         , kMagneticField     = 21
    };
              
  Bool_t fGraphsAreIni;              // Check whether graphs are initialized
  Bool_t fFitsAreIni;                // Check whether firs are initialized
  UInt_t fNAlias;                    // Number of aliases
      
    class AliTRDDataDCSdata {
     public:
      AliTRDDataDCSdata()
	:fFit(0x0)
        ,fGraph(0x0) { };
      virtual ~AliTRDDataDCSdata() { };
      TObjArray GetFit() const          { return fFit;      }
      TObjArray GetGraph() const        { return fGraph;    }
      TObject*  GetFit(Int_t i) const   { return fFit[i];   }
      TObject*  GetGraph(Int_t i) const { return fGraph[i]; }
     protected:
      TObjArray fFit;			// array of AliSplineFit
      TObjArray fGraph;		        // array of TGraph
    };      

    class AliTRDDataDCSconf {
     public:
      AliTRDDataDCSconf()
	:fAmanda(0)
	,fDataType(0)
	,fNChannel(0)
	,fEnableGraph(0)
	,fEnableFit(0)
	,fMinPoints(0)
	,fIter(0)
	,fMaxDelta(0)
	,fFitReq(0) { };
      virtual ~AliTRDDataDCSconf() { };
      TString  GetAmanda() const        { return fAmanda;      }
      Char_t   GetDataType() const      { return fDataType;    }
      UInt_t   GetNChannel() const      { return fNChannel;    }
      Bool_t   GetEnableGraph() const   { return fEnableGraph; }
      Bool_t   GetEnableFit() const     { return fEnableFit;   }
      Int_t    GetMinPoints() const     { return fMinPoints;   }
      Int_t    GetIter() const          { return fIter;        }
      Double_t GetMaxDelta() const      { return fMaxDelta;    }
      Int_t    GetFitReq() const        { return fFitReq;      }
      void     SetAmanda(TString s)     { fAmanda      = s;    }
      void     SetDataType(Char_t d)    { fDataType    = d;    }
      void     SetNChannel(UInt_t n)    { fNChannel    = n;    }
      void     SetEnableGraph(Bool_t e) { fEnableGraph = e;    }
      void     SetEnableFit(Bool_t e)   { fEnableFit   = e;    }
      void     SetMinPoints(Int_t m)    { fMinPoints   = m;    }
      void     SetIter(Int_t i)         { fIter        = i;    }
      void     SetMaxDelta(Double_t m)  { fMaxDelta    = m;    }
      void     SetFitReq(Int_t f)       { fFitReq      = f;    }
     protected:
        TString   fAmanda;	      	// amanda string
        Char_t    fDataType;		// 'c' for char, 'f' for float
        UInt_t    fNChannel;		// number of channel
        Bool_t    fEnableGraph;		// will be converted in TGraph
        Bool_t    fEnableFit;	        // will be converted in AliSplineFit
        Int_t     fMinPoints;		// minimum number of points per knot in fit
        Int_t     fIter;		// number of iterations for spline fit
        Double_t  fMaxDelta;	        // precision parameter for spline fit
        Int_t     fFitReq;	        // fit requirement, 2 = continuous 2nd derivative
      };      

      AliTRDDataDCSdata fDatas [22];	// configurations	
      AliTRDDataDCSconf fConfs [22];	// data arrays		
      
      ClassDef(AliTRDDataDCS,1)         // TRD calibration class
	
};

#endif
