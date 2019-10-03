#ifndef ALIPROTONFEEDDOWNANALYSISTASK_H
#define ALIPROTONFEEDDOWNANALYSISTASK_H
#include "AliAnalysisTask.h"

class TList;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliProtonFeedDownAnalysis;



class AliProtonFeedDownAnalysisTask : public AliAnalysisTask {
 public:
	AliProtonFeedDownAnalysisTask();
	AliProtonFeedDownAnalysisTask(const char *name);
	virtual ~AliProtonFeedDownAnalysisTask() {}
	
	virtual void   ConnectInputData(Option_t *);
	virtual void   CreateOutputObjects();
	virtual void   Exec(Option_t *option);
	virtual void   Terminate(Option_t *);
	
	void SetAnalysisObject(AliProtonFeedDownAnalysis *const analysis) {fProtonAnalysis = analysis;}
  
 private:
	AliESDEvent *fESD;    //ESD object 
	AliAODEvent *fAOD;    //AOD object
	AliMCEvent  *fMC;     //MC object 
	
	TList  *fList; //TList output object 
	
	AliProtonFeedDownAnalysis *fProtonAnalysis; //analysis object 
	
	TH1F *fStatHist;
	
	AliProtonFeedDownAnalysisTask(const AliProtonFeedDownAnalysisTask&); // not implemented
	AliProtonFeedDownAnalysisTask& operator=(const AliProtonFeedDownAnalysisTask&); // not implemented
	
ClassDef(AliProtonFeedDownAnalysisTask, 1); // example of analysis
};

#endif



