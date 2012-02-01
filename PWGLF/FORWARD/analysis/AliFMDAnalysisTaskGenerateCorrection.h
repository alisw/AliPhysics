#ifndef ALIFMDANALYSISTASKGENERATECORRECTION_H
#define ALIFMDANALYSISTASKGENERATECORRECTION_H

#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "AliFMDFloatMap.h"
#include "TH1F.h"

class AliFMDAnaCalibBackgroundCorrection;
class AliFMDAnaCalibEventSelectionEfficiency;

class AliFMDAnalysisTaskGenerateCorrection : public AliAnalysisTaskSE
{
 public:
  AliFMDAnalysisTaskGenerateCorrection();
    AliFMDAnalysisTaskGenerateCorrection(const char* name);
    ~AliFMDAnalysisTaskGenerateCorrection() {;}
 AliFMDAnalysisTaskGenerateCorrection(const AliFMDAnalysisTaskGenerateCorrection& o) : AliAnalysisTaskSE(), 
      fListOfHits(), 
      fListOfPrimaries(),
      fListOfCorrection(),
      fVertexBins(o.fVertexBins),
      fLastTrackByStrip(o.fLastTrackByStrip),
      fHitsByStrip(o.fHitsByStrip),
      fZvtxCut(o.fZvtxCut),
      fNvtxBins(o.fNvtxBins),
      fNbinsEta(o.fNbinsEta),
      fBackground(o.fBackground),
      fEventSelectionEff(o.fEventSelectionEff)
      {}
    AliFMDAnalysisTaskGenerateCorrection& operator=(const AliFMDAnalysisTaskGenerateCorrection&) { return *this; }
    
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
    AliFMDAnaCalibEventSelectionEfficiency*     fEventSelectionEff;
    
    ClassDef(AliFMDAnalysisTaskGenerateCorrection, 1);

};
#endif
