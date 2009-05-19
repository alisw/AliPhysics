#ifndef ALIFMDANALYSISTASKGENERATEBACKGROUND_H
#define ALIFMDANALYSISTASKGENERATEBACKGROUND_H

#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "AliFMDFloatMap.h"
#include "TH1F.h"

class AliFMDAnaCalibBackgroundCorrection;

class AliFMDAnalysisTaskGenerateBackground : public AliAnalysisTaskSE
{
 public:
  AliFMDAnalysisTaskGenerateBackground();
    AliFMDAnalysisTaskGenerateBackground(const char* name);
    ~AliFMDAnalysisTaskGenerateBackground() {;}
 AliFMDAnalysisTaskGenerateBackground(const AliFMDAnalysisTaskGenerateBackground& o) : AliAnalysisTaskSE(), 
      fListOfHits(), 
      fListOfPrimaries(),
      fListOfCorrection(),
      fVertexBins(o.fVertexBins),
      fLastTrackByStrip(o.fLastTrackByStrip),
      fZvtxCut(o.fZvtxCut),
      fNvtxBins(o.fNvtxBins),
      fNbinsEta(o.fNbinsEta),
      fBackground(o.fBackground)
      {}
    AliFMDAnalysisTaskGenerateBackground& operator=(const AliFMDAnalysisTaskGenerateBackground&) { return *this; }
    
    virtual void Init();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* /*option*/);
    void  Terminate(Option_t */*option*/);
    void SetZvtxCut(Float_t vtxcut) {fZvtxCut = vtxcut;}
    void SetNvtxBins(Int_t nvtxbins) {fNvtxBins = nvtxbins;}
    void SetNbinsEta(Int_t netabins) {fNbinsEta = netabins;}
    void ReadFromFile(const Char_t* filename = "background.root", Bool_t storeInOCDB = kFALSE, Int_t runNo=0);
 private:
    
    void GenerateCorrection();
    
    TList fListOfHits;
    TList fListOfPrimaries;
    TList fListOfCorrection;
    TH1F  fVertexBins;
    AliFMDFloatMap fLastTrackByStrip;
    AliFMDFloatMap fHitsByStrip;
    Float_t fZvtxCut;
    Int_t fNvtxBins;
    Int_t fNbinsEta;
    AliFMDAnaCalibBackgroundCorrection* fBackground;
    ClassDef(AliFMDAnalysisTaskGenerateBackground, 1);

};
#endif
