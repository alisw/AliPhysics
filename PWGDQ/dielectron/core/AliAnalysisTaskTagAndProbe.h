#ifndef AliAnalysisTaskTagAndProbe_cxx
#define AliAnalysisTaskTagAndProbe_cxx

//Author: Daiki Sekihata (Center for Nuclear Study, the University of Tokyo)
//daiki.sekihata@cern.ch

#include "AliAnalysisTaskSE.h"
#include "AliESDv0KineCuts.h"
#include "AliAODv0KineCuts.h"

class AliAnalysisTaskTagAndProbe : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskTagAndProbe();
    AliAnalysisTaskTagAndProbe(const char *name);
    virtual ~AliAnalysisTaskTagAndProbe();
    void SetCentralityMin(Float_t min) {fCentralityMin = min;}
    void SetCentralityMax(Float_t max) {fCentralityMax = max;}
    void SetDepthNMixed(Int_t Nmix)    {fNMixed        = Nmix;}
    void SetCentralityEstimator(TString estimator) {fEstimator = estimator;}
		void SetTriggerMask(UInt_t trigger){fTriggerMask = trigger;}
		void SetPhivCutRange(Float_t Mmax, Float_t PhiVmin){fMmax = Mmax, fPhiVmin = PhiVmin;}

		//void SetEventFilter(AliAnalysisFilter *filter){fEventFilter = filter;}
		//void SetTagFilter(AliAnalysisFilter *filter)  {fTagFilter   = filter;}
		//void SetProbeFilter(AliAnalysisFilter *filter){fProbeFilter = filter;}
		//void SetPassingProbeFilter(AliAnalysisFilter *filter){fPassingProbeFilter = filter;}
		//void SetPIDFilter(AliAnalysisFilter *filter){fPIDFilter = filter;}

		AliAnalysisFilter *GetEventFilter()     {return fEventFilter;}
		AliAnalysisFilter *GetTagFilter()       {return fTagFilter  ;}
		AliAnalysisFilter *GetProbeFilter(){return fProbeFilter;}
		AliAnalysisFilter *GetPassingProbeFilter()  {return fPassingProbeFilter;}
		AliAnalysisFilter *GetPIDFilter()  {return fPIDFilter;}

    void SetPIDCalibMode(Bool_t flag) {fPIDCalibMode = flag;}
		void SetPIDCaibinPU(Bool_t flag) {fPIDCalibinPU = flag;}

		void SetCentroidCorrFunctionPU(UInt_t detID, UInt_t parID, THnBase *fun, UInt_t var0, UInt_t var1, UInt_t var2, UInt_t var3, UInt_t var4) {

			// clone temporare histogram, otherwise it will not be streamed to file!
			TString key = Form("cntrd%d%d%d%d%d_%d%d",var0,var1,var2,var3,var4,detID,parID);
			printf("key = %s\n",key.Data());
			fPostPIDCntrdCorrPU[detID][parID] = (THnBase*)fun->Clone(key.Data());
			// check for corrections and add their variables to the fill map
			printf("detID = %u , parID = %u, POST PID CORRECTION in PU added for centroids:  ",detID,parID);
			switch(fun->GetNdimensions()) {
				case 5: printf(" %s, ",fun->GetAxis(4)->GetName());
				case 4: printf(" %s, ",fun->GetAxis(3)->GetName());
				case 3: printf(" %s, ",fun->GetAxis(2)->GetName());
				case 2: printf(" %s, ",fun->GetAxis(1)->GetName());
				case 1: printf(" %s " ,fun->GetAxis(0)->GetName());
			}
			printf("\n");
		}

		void SetWidthCorrFunctionPU(UInt_t detID, UInt_t parID, THnBase *fun, UInt_t var0, UInt_t var1, UInt_t var2, UInt_t var3, UInt_t var4) {

			// clone temporare histogram, otherwise it will not be streamed to file!
			TString key = Form("wdth%d%d%d%d%d_%d%d",var0,var1,var2,var3,var4,detID,parID);

			fPostPIDWdthCorrPU[detID][parID] = (THnBase*)fun->Clone(key.Data());
			// check for corrections and add their variables to the fill map
			printf("detID = %u , parID = %u, POST PID CORRECTION IN PU added for widths:  ",detID,parID);
			switch(fun->GetNdimensions()) {
				case 5: printf(" %s, ",fun->GetAxis(4)->GetName());
				case 4: printf(" %s, ",fun->GetAxis(3)->GetName());
				case 3: printf(" %s, ",fun->GetAxis(2)->GetName());
				case 2: printf(" %s, ",fun->GetAxis(1)->GetName());
				case 1: printf(" %s " ,fun->GetAxis(0)->GetName());
			}
			printf("\n");
		}

  protected:
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    void CutEfficiency();
    void TrackQA();
    void FillV0InfoESD();
    void FillV0InfoAOD();
		Double_t PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2);
    Bool_t HasConversionPointOnSPD(AliAODv0 *v0, AliAODTrack *legPos, AliAODTrack *legNeg){
      //if(v0->RadiusV0() < 4.5 && ( (legPos->HasSharedPointOnITSLayer(0) && legNeg->HasSharedPointOnITSLayer(0)) || (legPos->HasSharedPointOnITSLayer(1) && legNeg->HasSharedPointOnITSLayer(1)) )) return kTRUE;//SPD0
      if(v0->RadiusV0() < 4.5 && ( legPos->HasSharedPointOnITSLayer(0) && legNeg->HasSharedPointOnITSLayer(0) )) return kTRUE;//SPD0
      else if((6.5 < v0->RadiusV0() && v0->RadiusV0() < 8.0) && legPos->HasSharedPointOnITSLayer(1) && legNeg->HasSharedPointOnITSLayer(1)) return kTRUE;//SPD1
      else return kFALSE;//none of above
    }

    void FillHistogramTH1(TList *list, const Char_t *name, Double_t x, Double_t w=1., Option_t *opt = "") const ;
    void FillHistogramTH2(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t w=1., Option_t *opt = "") const ;
    void FillHistogramTH3(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t z, Double_t w=1., Option_t *opt = "") const ;
    void FillSparse(TList *list, const Char_t *name, Double_t *x, Double_t w=1.) const;

  protected:
    THashList *fOutputContainer;
		UInt_t fTriggerMask;
    AliVEvent *fEvent;
    AliESDEvent *fESDEvent;
    AliAODEvent *fAODEvent;
    TString fEstimator;//V0[M|A|C], ZN[A|C], CL[0|1]
    AliMultSelection *fMultSelection;
    Float_t fCentrality;
    Float_t fCentralityMin;
    Float_t fCentralityMax;
    Int_t fNMixed;
    Double_t fVertex[3];
    Int_t fZvtxBin;
		AliPIDResponse *fPIDResponse;
    TBits *fUsedVars;
		Bool_t fPIDCalibinPU;
		THnBase *fPostPIDCntrdCorrPU[15][15];   // post pid correction object //multi-dimension for pileup, 3 for TPC/ITS/TOF, 5 for e/mu/pi/k/p
		THnBase *fPostPIDWdthCorrPU[15][15];    // post pid correction object //multi-dimension for pileup, 3 for TPC/ITS/TOF, 5 for e/mu/pi/k/p
		TObjArray *fTagTrackArray;//ele or pos
		TObjArray *fProbeTrackArray;//ele or pos
		TObjArray *fPassingProbeTrackArray;//ele or pos
		TList *fEventList[2][10];//p/pp x zvtx
		AliAnalysisFilter *fEventFilter;
		AliAnalysisFilter *fTagFilter;
		AliAnalysisFilter *fProbeFilter;
		AliAnalysisFilter *fPassingProbeFilter;
		AliAnalysisFilter *fPIDFilter;
		Float_t fMmax;
		Float_t fPhiVmin;
    AliESDv0KineCuts *fESDv0KineCuts;
    AliAODv0KineCuts *fAODv0KineCuts;
    Bool_t fPIDCalibMode;

  private:
    AliAnalysisTaskTagAndProbe(const AliAnalysisTaskTagAndProbe&);
    AliAnalysisTaskTagAndProbe& operator=(const AliAnalysisTaskTagAndProbe&);

    ClassDef(AliAnalysisTaskTagAndProbe, 4);
};


#endif
