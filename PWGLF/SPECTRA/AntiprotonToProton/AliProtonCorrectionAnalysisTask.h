#ifndef ALIPROTONCORRECTIONTASK_H
#define ALIPROTONCORRECTIONTASK_H
#include "AliAnalysisTask.h"

class TList;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliProtonFeedDownAnalysis;
class AliProtonAbsorptionCorrection;
class AliProtonSpectraCorrection; 
class AliProtonAnalysisBase;



class AliProtonCorrectionAnalysisTask : public AliAnalysisTask {
 public:
	AliProtonCorrectionAnalysisTask();
	AliProtonCorrectionAnalysisTask(const char *name);
	virtual ~AliProtonCorrectionAnalysisTask() {}
	
	virtual void   ConnectInputData(Option_t *);
	virtual void   CreateOutputObjects();
	virtual void   Exec(Option_t *option);
	virtual void   Terminate(Option_t *);
	
	void SetAnalysisObjectAbsorptionCorrection(AliProtonAbsorptionCorrection* const analysis) ;
	void SetAnalysisObjectFeedDown(AliProtonFeedDownAnalysis* const analysis); 
	void SetAnalysisObjectSpectraCorrection(AliProtonSpectraCorrection* const analysis);
    	void SetBaseAnalysis(AliProtonAnalysisBase* const baseAnalysis) { fProtonAnalysisBase = baseAnalysis;}
 private:
	AliESDEvent *fESD;    //ESD object 
	AliAODEvent *fAOD;    //AOD object
	AliMCEvent  *fMC;     //MC object 
	
	TList  *fList; //TList output object 
	
	AliProtonAnalysisBase *fProtonAnalysisBase; //base analysis object
	AliProtonAbsorptionCorrection *fProtonAbsorptionCorrection;//analysis object 
	AliProtonFeedDownAnalysis *fProtonFeedDownAnalysis; //analysis object 
	AliProtonSpectraCorrection *fProtonSpectraCorrection;//analysis object 
	
	TH1F *fStatHist;
	
	
	Bool_t fIsOn_AliProtonAbsorptionCorrection;
	Bool_t fIsOn_AliProtonFeedDownAnalysis;
	Bool_t fIsOn_AliProtonSpectraCorrection;
	
	AliProtonCorrectionAnalysisTask(const AliProtonCorrectionAnalysisTask&); // not implemented
	AliProtonCorrectionAnalysisTask& operator=(const AliProtonCorrectionAnalysisTask&); // not implemented
	
ClassDef(AliProtonCorrectionAnalysisTask, 1); // example of analysis
};

#endif



