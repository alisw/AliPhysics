#ifndef AliMultiplictyLoaderTask_H
#define AliMultiplictyLoaderTask_H
class TH1F;
class TH2F;
class TH3D;	
class AliESDEvent;
class  AliESDtrackCuts;
class AliPPVsMultUtils;
template<class T> class TParameter;
//class TParameter; 
#include "AliAnalysisTaskSE.h"



class AliMultiplictyLoaderTask : public AliAnalysisTaskSE 
{
 	public:
 		 AliMultiplictyLoaderTask(const char *name = "AliMultiplictyLoaderTask");
 		 virtual ~AliMultiplictyLoaderTask();
  
 		 //virtual void   ConnectInputData(Option_t *);
  		 virtual void   UserCreateOutputObjects();
 		 virtual void   UserExec(Option_t *option);
 		 virtual void   Terminate(Option_t *){} 
  		//virtual void   LocalInit();
 		void SetCentEstimator(TString cent = "V0M")    {fCentEstimator = cent; }
 		void SetUseAliPPVsMultUtils(Bool_t flag)	{fUseAliPPVsMultUtils=flag;}
		void SetDonotusetrackelts(Bool_t flag)        {fDonotusetrackelts=flag;}
 
	private:
		AliESDEvent *fESD;    //! ESD object
 	        AliPPVsMultUtils* fAliPPVsMultUtils; // tool to get V0M multiplicty/centrailty 
 		TString fCentEstimator; // type of the centrailty estimator  
 		Bool_t fUseAliPPVsMultUtils; // if true uses the centrality from AliPPVsMultUtils
 		TParameter<Double_t>* fcentvalue; // value of centrailty 
	        TParameter<Int_t>* fncharged05value; // value of Nch for |eta|<0.5
		TParameter<Int_t>* fncharged08value; // value of Nch for |eta|<0.8
		Bool_t fFirstEvent; // first Event Flag
		Bool_t fDonotusetrackelts; // if true the events with only spd veretx will have <0 value oe ref multiplicty 
		AliMultiplictyLoaderTask(const AliMultiplictyLoaderTask&); // private copy const
 		AliMultiplictyLoaderTask& operator=(const AliMultiplictyLoaderTask&); // private = operator

	
	 ClassDef(AliMultiplictyLoaderTask, 2); 
};

#endif
