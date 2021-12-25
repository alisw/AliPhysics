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
    void SetMC(Bool_t flag) {fIsMC = flag;}
    void SetCentralityMin(Float_t min) {fCentralityMin = min;}
    void SetCentralityMax(Float_t max) {fCentralityMax = max;}
    void SetDepthNMixed(Int_t Nmix)    {fNMixed        = Nmix;}
    void SetCentralityEstimator(TString estimator) {fEstimator = estimator;}
    void SetTriggerMask(UInt_t trigger){fTriggerMask = trigger;}
    void SetPhivCutRange(Float_t Mmax, Float_t PhiVmin){fMmax = Mmax, fPhiVmin = PhiVmin;}

    //void SetEventFilter(AliAnalysisFilter *filter){fEventFilter = filter;}
    //void SetTagFilter(AliAnalysisFilter *filter)  {fTagFilter   = filter;}
    //void SetProbeFilter(AliAnalysisFilter *filter){fProbeFilter = filter;}

    AliAnalysisFilter *GetEventFilter()     {return fEventFilter;}
    AliAnalysisFilter *GetTagFilter()       {return fTagFilter  ;}
    AliAnalysisFilter *GetProbeFilter(){return fProbeFilter;}
    void AddPassingProbeFilter(AliAnalysisFilter *filter) {fListPassingProbeFilters->Add(filter);}

    void SetPIDCalibMode(Bool_t flag) {fPIDCalibMode = flag;}
    void SetPIDCalibinPU(Bool_t flag) {fPIDCalibinPU = flag;}

    void SetCentroidCorrFunctionPU(UInt_t detID, UInt_t parID, THnBase *fun, UInt_t var0, UInt_t var1, UInt_t var2, UInt_t var3, UInt_t var4) {

      // clone temporare histogram, otherwise it will not be streamed to file!
      TString key = Form("cntrd%d%d%d%d%d_%d%d",var0,var1,var2,var3,var4,detID,parID);
      //printf("key = %s\n",key.Data());
      fun->GetAxis(4)->SetUniqueID(var4);
      fun->GetAxis(3)->SetUniqueID(var3);
      fun->GetAxis(2)->SetUniqueID(var2);
      fun->GetAxis(1)->SetUniqueID(var1);
      fun->GetAxis(0)->SetUniqueID(var0);
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
      fUsedVars->SetBitNumber(var0, kTRUE);
      fUsedVars->SetBitNumber(var1, kTRUE);
      fUsedVars->SetBitNumber(var2, kTRUE);
      fUsedVars->SetBitNumber(var3, kTRUE);
      fUsedVars->SetBitNumber(var4, kTRUE);
    }

    void SetWidthCorrFunctionPU(UInt_t detID, UInt_t parID, THnBase *fun, UInt_t var0, UInt_t var1, UInt_t var2, UInt_t var3, UInt_t var4) {

      // clone temporare histogram, otherwise it will not be streamed to file!
      TString key = Form("wdth%d%d%d%d%d_%d%d",var0,var1,var2,var3,var4,detID,parID);
      fun->GetAxis(4)->SetUniqueID(var4);
      fun->GetAxis(3)->SetUniqueID(var3);
      fun->GetAxis(2)->SetUniqueID(var2);
      fun->GetAxis(1)->SetUniqueID(var1);
      fun->GetAxis(0)->SetUniqueID(var0);
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
      fUsedVars->SetBitNumber(var0, kTRUE);
      fUsedVars->SetBitNumber(var1, kTRUE);
      fUsedVars->SetBitNumber(var2, kTRUE);
      fUsedVars->SetBitNumber(var3, kTRUE);
      fUsedVars->SetBitNumber(var4, kTRUE);
    }

    void SetCentroidCorrFunctionTPC(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz) {
      //UInt_t valType[20] = {0};
      //valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
      // clone temporare histogram, otherwise it will not be streamed to file!
      TString key = Form("cntrd%d%d%d",varx,vary,varz);
      fun->GetZaxis()->SetUniqueID(varz);
      fun->GetYaxis()->SetUniqueID(vary);
      fun->GetXaxis()->SetUniqueID(varx);
      fPostPIDCntrdCorrTPC = (TH1*)fun->Clone(key.Data());
      // check for corrections and add their variables to the fill map
      if(fPostPIDCntrdCorrTPC)  {
        printf("POST TPC PID CORRECTION added for centroids:  ");
        switch(fPostPIDCntrdCorrTPC->GetDimension()) {
          case 3: printf(" %s, ",fPostPIDCntrdCorrTPC->GetZaxis()->GetName());
          case 2: printf(" %s, ",fPostPIDCntrdCorrTPC->GetYaxis()->GetName());
          case 1: printf(" %s ",fPostPIDCntrdCorrTPC->GetXaxis()->GetName());
        }
        printf("\n");
        fUsedVars->SetBitNumber(varx, kTRUE);
        fUsedVars->SetBitNumber(vary, kTRUE);
        fUsedVars->SetBitNumber(varz, kTRUE);
      }
    } 

    void SetWidthCorrFunctionTPC(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz) {
      //UInt_t valType[20] = {0};
      //valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
      // clone temporare histogram, otherwise it will not be streamed to file!
      TString key = Form("wdth%d%d%d",varx,vary,varz);
      fun->GetZaxis()->SetUniqueID(varz);
      fun->GetYaxis()->SetUniqueID(vary);
      fun->GetXaxis()->SetUniqueID(varx);
      fPostPIDWdthCorrTPC = (TH1*)fun->Clone(key.Data());
      // check for corrections and add their variables to the fill map
      if(fPostPIDWdthCorrTPC)  {
        printf("POST TPC PID CORRECTION added for widths:  ");
        switch(fPostPIDWdthCorrTPC->GetDimension()) {
          case 3: printf(" %s, ",fPostPIDWdthCorrTPC->GetZaxis()->GetName());
          case 2: printf(" %s, ",fPostPIDWdthCorrTPC->GetYaxis()->GetName());
          case 1: printf(" %s ",fPostPIDWdthCorrTPC->GetXaxis()->GetName());
        }
        printf("\n");
        fUsedVars->SetBitNumber(varx, kTRUE);
        fUsedVars->SetBitNumber(vary, kTRUE);
        fUsedVars->SetBitNumber(varz, kTRUE);
      }
    }

    void SetCentroidCorrFunctionITS(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz) {
      //UInt_t valType[20] = {0};
      //valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
      // clone temporare histogram, otherwise it will not be streamed to file!
      TString key = Form("cntrd%d%d%d",varx,vary,varz);
      fun->GetZaxis()->SetUniqueID(varz);
      fun->GetYaxis()->SetUniqueID(vary);
      fun->GetXaxis()->SetUniqueID(varx);
      fPostPIDCntrdCorrITS = (TH1*)fun->Clone(key.Data());
      // check for corrections and add their variables to the fill map
      if(fPostPIDCntrdCorrITS)  {
        printf("POST ITS PID CORRECTION added for centroids:  ");
        switch(fPostPIDCntrdCorrITS->GetDimension()) {
          case 3: printf(" %s, ",fPostPIDCntrdCorrITS->GetZaxis()->GetName());
          case 2: printf(" %s, ",fPostPIDCntrdCorrITS->GetYaxis()->GetName());
          case 1: printf(" %s ",fPostPIDCntrdCorrITS->GetXaxis()->GetName());
        }
        printf("\n");
        fUsedVars->SetBitNumber(varx, kTRUE);
        fUsedVars->SetBitNumber(vary, kTRUE);
        fUsedVars->SetBitNumber(varz, kTRUE);
      }
    } 

    void SetWidthCorrFunctionITS(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz) {
      //UInt_t valType[20] = {0};
      //valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
      // clone temporare histogram, otherwise it will not be streamed to file!
      TString key = Form("wdth%d%d%d",varx,vary,varz);
      fun->GetZaxis()->SetUniqueID(varz);
      fun->GetYaxis()->SetUniqueID(vary);
      fun->GetXaxis()->SetUniqueID(varx);
      fPostPIDWdthCorrITS = (TH1*)fun->Clone(key.Data());
      // check for corrections and add their variables to the fill map
      if(fPostPIDWdthCorrITS)  {
        printf("POST ITS PID CORRECTION added for widths:  ");
        switch(fPostPIDWdthCorrITS->GetDimension()) {
          case 3: printf(" %s, ",fPostPIDWdthCorrITS->GetZaxis()->GetName());
          case 2: printf(" %s, ",fPostPIDWdthCorrITS->GetYaxis()->GetName());
          case 1: printf(" %s ",fPostPIDWdthCorrITS->GetXaxis()->GetName());
        }
        printf("\n");
        fUsedVars->SetBitNumber(varx, kTRUE);
        fUsedVars->SetBitNumber(vary, kTRUE);
        fUsedVars->SetBitNumber(varz, kTRUE);
      }
    }

    void SetCentroidCorrFunctionTOF(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz) {
      //UInt_t valType[20] = {0};
      //valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
      // clone temporare histogram, otherwise it will not be streamed to file!
      TString key = Form("cntrd%d%d%d",varx,vary,varz);
      fun->GetZaxis()->SetUniqueID(varz);
      fun->GetYaxis()->SetUniqueID(vary);
      fun->GetXaxis()->SetUniqueID(varx);
      fPostPIDCntrdCorrTOF = (TH1*)fun->Clone(key.Data());
      // check for corrections and add their variables to the fill map
      if(fPostPIDCntrdCorrTOF)  {
        printf("POST TOF PID CORRECTION added for centroids:  ");
        switch(fPostPIDCntrdCorrTOF->GetDimension()) {
          case 3: printf(" %s, ",fPostPIDCntrdCorrTOF->GetZaxis()->GetName());
          case 2: printf(" %s, ",fPostPIDCntrdCorrTOF->GetYaxis()->GetName());
          case 1: printf(" %s ",fPostPIDCntrdCorrTOF->GetXaxis()->GetName());
        }
        printf("\n");
        fUsedVars->SetBitNumber(varx, kTRUE);
        fUsedVars->SetBitNumber(vary, kTRUE);
        fUsedVars->SetBitNumber(varz, kTRUE);
      }
    } 

    void SetWidthCorrFunctionTOF(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz) {
      //UInt_t valType[20] = {0};
      //valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
      // clone temporare histogram, otherwise it will not be streamed to file!
      TString key = Form("wdth%d%d%d",varx,vary,varz);
      fun->GetZaxis()->SetUniqueID(varz);
      fun->GetYaxis()->SetUniqueID(vary);
      fun->GetXaxis()->SetUniqueID(varx);
      fPostPIDWdthCorrTOF = (TH1*)fun->Clone(key.Data());
      // check for corrections and add their variables to the fill map
      if(fPostPIDWdthCorrTOF)  {
        printf("POST TOF PID CORRECTION added for widths:  ");
        switch(fPostPIDWdthCorrTOF->GetDimension()) {
          case 3: printf(" %s, ",fPostPIDWdthCorrTOF->GetZaxis()->GetName());
          case 2: printf(" %s, ",fPostPIDWdthCorrTOF->GetYaxis()->GetName());
          case 1: printf(" %s ",fPostPIDWdthCorrTOF->GetXaxis()->GetName());
        }
        printf("\n");
        fUsedVars->SetBitNumber(varx, kTRUE);
        fUsedVars->SetBitNumber(vary, kTRUE);
        fUsedVars->SetBitNumber(varz, kTRUE);
      }
    }

  protected:
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void ProcessMC(Option_t *option);
    void CutEfficiency(TObjArray *arr1, TObjArray *arr2, const TString str);
    void TrackQA();
    void FillV0InfoESD();
    void FillV0InfoAOD();
    void GetMCInfoESD();
    void GetMCInfoAOD();
    Double_t PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2);

    Bool_t HasConversionPointOnSPD(AliESDv0 *v0, AliESDtrack *legPos, AliESDtrack *legNeg){
           if( (3.5 < v0->GetRr() && v0->GetRr() < 4.3) && (legPos->HasSharedPointOnITSLayer(0) && legNeg->HasSharedPointOnITSLayer(0)) ) return kTRUE;//SPD0
      else if( (6.9 < v0->GetRr() && v0->GetRr() < 7.7) && (legPos->HasSharedPointOnITSLayer(1) && legNeg->HasSharedPointOnITSLayer(1)) ) return kTRUE;//SPD1
      else return kFALSE;//none of above
    }

    Bool_t HasConversionPointOnSPD(AliAODv0 *v0, AliAODTrack *legPos, AliAODTrack *legNeg){
           if( (3.5 < v0->RadiusV0() && v0->RadiusV0() < 4.3) && (legPos->HasSharedPointOnITSLayer(0) && legNeg->HasSharedPointOnITSLayer(0)) ) return kTRUE;//SPD0
      else if( (6.9 < v0->RadiusV0() && v0->RadiusV0() < 7.7) && (legPos->HasSharedPointOnITSLayer(1) && legNeg->HasSharedPointOnITSLayer(1)) ) return kTRUE;//SPD1

      //     if( (3.5 < v0->RadiusV0() && v0->RadiusV0() < 4.3) && ( (legPos->HasSharedPointOnITSLayer(0) && legNeg->HasSharedPointOnITSLayer(0)) || (legPos->HasSharedPointOnITSLayer(1) && legNeg->HasSharedPointOnITSLayer(1)) ) ) return kTRUE;//SPD0
      //else if( (6.9 < v0->RadiusV0() && v0->RadiusV0() < 7.7) && ( (legPos->HasSharedPointOnITSLayer(1) && legNeg->HasSharedPointOnITSLayer(1)) || (legPos->HasSharedPointOnITSLayer(2) && legNeg->HasSharedPointOnITSLayer(2)) ) ) return kTRUE;//SPD1

      //if(3.5 < v0->RadiusV0() && v0->RadiusV0() < 4.3 )       return kTRUE;//SPD0
      //else if(6.9  < v0->RadiusV0() && v0->RadiusV0() < 7.7 ) return kTRUE;//SPD1
      //else if(14.  < v0->RadiusV0() && v0->RadiusV0() < 16. ) return kTRUE;//SDD0
      //else if(23.3 < v0->RadiusV0() && v0->RadiusV0() < 24.7) return kTRUE;//SDD1
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
    AliMCEvent *fMCEvent;
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
    TH1 *fPostPIDCntrdCorrTPC;
    TH1 *fPostPIDWdthCorrTPC;
    TH1 *fPostPIDCntrdCorrITS;
    TH1 *fPostPIDWdthCorrITS;
    TH1 *fPostPIDCntrdCorrTOF;
    TH1 *fPostPIDWdthCorrTOF;
    TObjArray *fTrackArrayPos;
    TObjArray *fTrackArrayNeg;
    TList *fEventList[2][10];//ele/pos x zvtx
    AliAnalysisFilter *fEventFilter;
    AliAnalysisFilter *fTagFilter;
    AliAnalysisFilter *fProbeFilter;
    TList *fListPassingProbeFilters;
    Float_t fMmax;
    Float_t fPhiVmin;
    AliESDtrackCuts *fESDtrackCutsGlobalNoDCA;
    AliESDv0KineCuts *fESDv0KineCuts;
    AliAODv0KineCuts *fAODv0KineCuts;
    AliStack *fMCArrayESD;     //MC particles array in ESD
    TClonesArray *fMCArrayAOD; //MC particles array in AOD
    Bool_t fPIDCalibMode;
    Bool_t fIsMC;

  private:
    AliAnalysisTaskTagAndProbe(const AliAnalysisTaskTagAndProbe&);
    AliAnalysisTaskTagAndProbe& operator=(const AliAnalysisTaskTagAndProbe&);

    ClassDef(AliAnalysisTaskTagAndProbe, 11);
};


#endif
