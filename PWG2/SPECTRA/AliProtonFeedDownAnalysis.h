#ifndef ALIPROTONFEEDDOWNANALYSIS_H
#define ALIPROTONFEEDDOWNANALYSIS_H

#include "TObject.h"
#include "TH1I.h"
#include "AliCFContainer.h"
class TF1;
class TH2D;
class TH1F;
class TList;

//#include "AliPID.h"
class AliPID;

class AliCFDataGrid;
class AliAODEvent;
class AliAODtrack;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliStack;
class AliESDVertex;
class AliProtonAnalysisBase;

class AliProtonFeedDownAnalysis : public TObject 
{
 	public:
	enum 
		{
    			kAll      = 0, //without protons reject the secondaries from hadronic inter. 
    			kPrimary = 1,
    			kFromLambda        = 2,
   			kSelected = 3 ,
    			kSelectedfromLambda = 4,
			kNSteps=5
		};
		AliProtonFeedDownAnalysis();
		//AliProtonFeedDownAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
		virtual ~AliProtonFeedDownAnalysis();
		
		void SetBaseAnalysis(AliProtonAnalysisBase * const baseAnalysis) {fProtonAnalysisBase = baseAnalysis;}
		AliProtonAnalysisBase *GetProtonAnalysisBaseObject() const {return fProtonAnalysisBase;}
		
		void InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY, Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
		void Analyze(AliESDEvent *fESD, const AliESDVertex *vertex,AliStack *stack);
		void Analyze(AliAODEvent *fAOD);
		void Analyze(AliStack *stack);
		
		AliCFContainer *GetProtonContainer() const {return fProtonContainer;}
		AliCFContainer *GetAntiProtonContainer() const {return fAntiProtonContainer;}
		
		void SetWeightFunction(TF1* weightfunction){fweightfunction=weightfunction;}
		TF1* GetWeightFunction() const{return fweightfunction;}	
		TH2F* GetLambdaHist() const{return fLambda;}
		TH2F* GetLambdaweightedHist() const{return fLambdaweighted;}
		TH2F* GetAntiLambdaHist() const{return fAntiLambda;}
		TH2F* GetAntiLambdaweightedHist() const{return fAntiLambdaweighted;}
 	private:
		AliProtonFeedDownAnalysis(const AliProtonFeedDownAnalysis&); // Not implemented
		AliProtonFeedDownAnalysis& operator=(const AliProtonFeedDownAnalysis&); // Not implemented
		
		AliProtonAnalysisBase *fProtonAnalysisBase;//base analysis object
		Int_t LambdaIsMother(Int_t numbe, AliStack *stack); 
		Float_t GetWeightforProton(Int_t number, AliStack *stack);
		Float_t GetWeightforLambda(Float_t pt);
		
		Int_t fNBinsY; //number of bins in y or eta
		Double_t fMinY, fMaxY; //min & max value of y or eta
		Int_t fNBinsPt;  //number of bins in pT
		Double_t fMinPt, fMaxPt; //min & max value of pT
		
		//Analysis containers
		AliCFContainer *fProtonContainer; //container for protons
		AliCFContainer *fAntiProtonContainer; //container for antiprotons
		
		TF1*fweightfunction; 
		TH2F *fLambda;
		TH2F *fLambdaweighted;
		TH2F *fAntiLambda;
		TH2F *fAntiLambdaweighted;
		
  ClassDef(AliProtonFeedDownAnalysis,2);
};

#endif

