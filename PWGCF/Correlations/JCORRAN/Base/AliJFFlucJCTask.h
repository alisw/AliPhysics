#ifndef ALIJFFLUCJCTASK_H
#define ALIJFFLUCJCTASK_H

class AliJFFlucJCTask : public AliAnalysisTaskSE{
public:
	AliJFFlucJCTask();
	AliJFFlucJCTask(const char *);
	AliJFFlucJCTask(const AliJFFlucJCTask &);
	//AliJFFlucJCTask & operator=(const AliJFFlucJCTask &);
	virtual ~AliJFFlucJCTask();
	//
	void UserCreateOutputObjects(); 
	void UserExec(Option_t *);

  void SetJCatalystTaskName(TString name) {fJCatalystTaskName = name;}	/// New

	enum{
		FLUC_PHI_CORRECTION  = 0x1,
		FLUC_EBE_WEIGHTING   = 0x2,
	};
	inline void AddFlags(UInt_t _flags){
		flags |= _flags;
	}
	enum SUBEVENT{
		SUBEVENT_A = 0x1,
		SUBEVENT_B = 0x2
	};
	inline void SelectSubevents(UInt_t _subeventMask){
		subeventMask = _subeventMask;
	}

private:
	AliJCatalystTask *fJCatalystTask;
	TString fTaskName;
	TString fJCatalystTaskName;
	AliJFFlucAnalysis *pfa;
	UInt_t subeventMask;
	UInt_t flags;

	ClassDef(AliJFFlucJCTask,2);
};

#endif

