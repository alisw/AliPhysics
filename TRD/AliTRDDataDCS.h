#ifndef AliTRDDataDCS_H
#define AliTRDDAtaDCS_H

#include <TNamed.h>
#include <TObjArray.h>
#include <TString.h>

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
  TGraph * GetGraphChamberByteStatus (UInt_t iSensor) const	        {return GetGraph (chamberByteStatus, iSensor);}
  TGraph * GetGraphPreTrigger () const 					{return GetGraph (preTrigger);}
  TGraph * GetGraphGoofyHv () const			        	{return GetGraph (goofyHv);}
  TGraph * GetGraphGoofyPeakPos (UInt_t iSensor) const	                {return GetGraph (goofyPeakPos, iSensor);}
  TGraph * GetGraphGoofyPeakArea (UInt_t iSensor) const 	        {return GetGraph (goofyPeakArea, iSensor);}
  TGraph * GetGraphGoofyTemp (UInt_t iSensor) const		        {return GetGraph (goofyTemp, iSensor);}
  TGraph * GetGraphGoofyPressure () const 				{return GetGraph (goofyPressure);}
  TGraph * GetGraphGoofyVelocity () const 				{return GetGraph (goofyVelocity);}
  TGraph * GetGraphGoofyGain (UInt_t iSensor) const		        {return GetGraph (goofyGain, iSensor);}
  TGraph * GetGraphGoofyCO2 ()  const					{return GetGraph (goofyCO2);}
  TGraph * GetGraphGoofyN2 () const 					{return GetGraph (goofyN2);}
  TGraph * GetGraphGasO2 () const 					{return GetGraph (gasO2);}
  TGraph * GetGraphGasOverpressure () const 				{return GetGraph (gasOverpressure);}
  TGraph * GetGraphEnvTemp (UInt_t iSensor) const			{return GetGraph (envTemp, iSensor);}
  TGraph * GetGraphHvAnodeImon (UInt_t iSensor) const		        {return GetGraph (hvAnodeImon, iSensor);}
  TGraph * GetGraphHvDriftImon (UInt_t iSensor) const		        {return GetGraph (hvDriftImon, iSensor);}
  TGraph * GetGraphHvAnodeUmon (UInt_t iSensor) const		        {return GetGraph (hvAnodeUmon, iSensor);}
  TGraph * GetGraphHvDriftUmon (UInt_t iSensor) const		        {return GetGraph (hvDriftUmon, iSensor);}
  TGraph * GetGraphAdcClkPhase () const 				{return GetGraph (adcClkPhase);}
  TGraph * GetGraphAtmPressure () const 				{return GetGraph (atmPressure);}
  TGraph * GetGraphLuminosity () const 					{return GetGraph (luminosity);}
  TGraph * GetGraphMagneticField () const				{return GetGraph (magneticField);}
  
  AliSplineFit * GetFitChamberByteStatus (UInt_t iSensor) const	        {return GetFit (chamberByteStatus, iSensor);}
  AliSplineFit * GetFitPreTrigger () const 				{return GetFit (preTrigger);}
  AliSplineFit * GetFitGoofyHv () const					{return GetFit (goofyHv);}
  AliSplineFit * GetFitGoofyPeakPos (UInt_t iSensor) const	        {return GetFit (goofyPeakPos, iSensor);}
  AliSplineFit * GetFitGoofyPeakArea (UInt_t iSensor) const 	        {return GetFit (goofyPeakArea, iSensor);}
  AliSplineFit * GetFitGoofyTemp (UInt_t iSensor) const		        {return GetFit (goofyTemp, iSensor);}
  AliSplineFit * GetFitGoofyPressure () const 				{return GetFit (goofyPressure);}
  AliSplineFit * GetFitGoofyVelocity () const 				{return GetFit (goofyVelocity);}
  AliSplineFit * GetFitGoofyGain (UInt_t iSensor) const		        {return GetFit (goofyGain, iSensor);}
  AliSplineFit * GetFitGoofyCO2 ()  const				{return GetFit (goofyCO2);}
  AliSplineFit * GetFitGoofyN2 () const 				{return GetFit (goofyN2);}
  AliSplineFit * GetFitGasO2 () const 					{return GetFit (gasO2);}
  AliSplineFit * GetFitGasOverpressure () const 			{return GetFit (gasOverpressure);}
  AliSplineFit * GetFitEnvTemp (UInt_t iSensor) const			{return GetFit (envTemp, iSensor);}
  AliSplineFit * GetFitHvAnodeImon (UInt_t iSensor) const		{return GetFit (hvAnodeImon, iSensor);}
  AliSplineFit * GetFitHvDriftImon (UInt_t iSensor) const		{return GetFit (hvDriftImon, iSensor);}
  AliSplineFit * GetFitHvAnodeUmon (UInt_t iSensor) const		{return GetFit (hvAnodeUmon, iSensor);}
  AliSplineFit * GetFitHvDriftUmon (UInt_t iSensor) const		{return GetFit (hvDriftUmon, iSensor);}
  AliSplineFit * GetFitAdcClkPhase () const 				{return GetFit (adcClkPhase);}
  AliSplineFit * GetFitAtmPressure () const 				{return GetFit (atmPressure);}
  AliSplineFit * GetFitLuminosity () const 				{return GetFit (luminosity);}
  AliSplineFit * GetFitMagneticField () const				{return GetFit (magneticField);}
  
  void Print (Option_t* option = "") const;
  
  
  public :
    
    const UInt_t  chamberByteStatus;
    const UInt_t  preTrigger;
    const UInt_t  goofyHv;
    const UInt_t  goofyPeakPos;
    const UInt_t  goofyPeakArea;
    const UInt_t  goofyTemp;
    const UInt_t  goofyPressure;
    const UInt_t  goofyVelocity;
    const UInt_t  goofyGain;
    const UInt_t  goofyCO2;
    const UInt_t  goofyN2;
    const UInt_t  gasO2;
    const UInt_t  gasOverpressure;
    const UInt_t  envTemp ;
    const UInt_t  hvAnodeImon;
    const UInt_t  hvDriftImon;
    const UInt_t  hvAnodeUmon;
    const UInt_t  hvDriftUmon;
    const UInt_t  adcClkPhase;
    const UInt_t  atmPressure;
    const UInt_t  luminosity;
    const UInt_t  magneticField;
    const UInt_t  fNAlias;
    
    
    
    protected :
      
      private :
      
      TGraph * FindAndMakeGraph (TMap * dcsMap, const char * amandaStr,
				 char dataType);
    AliSplineFit * Fit (TGraph * graph,
			Int_t  kMinPoints, Int_t  kIter, 
			Double_t  kMaxDelta, Int_t  kFitReq);
    
    void Init ();
    void InitFits ();
    void InitGraphs ();
    
    void SetConf (UInt_t iAlias, const char * amanda, char dataType, UInt_t nChannel, 
		  Bool_t enableGraph, Bool_t enableFit, Int_t kMinPoints, 
		  Int_t kIter, Double_t kMaxDelta, Int_t kFitReq);
    
    private :
      
      Bool_t graphsAreIni;
      Bool_t fitsAreIni;
      
      struct data {
	TObjArray fit;			// array of AliSplineFit
	TObjArray graph;		// array of TGraph
      };
      
      struct conf {
	TString amanda;			// amanda string
	char dataType;			// 'c' for char, 'f' for float
	UInt_t nChannel;		// number of channel
	Bool_t enableGraph;		// will be converted in TGraph
	Bool_t enableFit;		// will be converted in AliSplineFit
	Int_t  kMinPoints;		// minimum number of points per knot in fit
	Int_t  kIter;			// number of iterations for spline fit
	Double_t  kMaxDelta;	        // precision parameter for spline fit
	Int_t  kFitReq;			// fit requirement, 2 = continuous 2nd derivative
      };
      
      data fDatas [22];		
      conf fConfs [22];			
      
      ClassDef(AliTRDDataDCS,1)         //  TRD calibration class
	
	};

#endif
