#ifndef ALIFMDANACALIBEVENTSELECTIONEFFICIENCY_H
#define ALIFMDANACALIBEVENTSELECTIONEFFICIENCY_H

#include <TObject.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TList.h>
#include <TH2F.h>
#include <TString.h>
class TBrowser;

class AliFMDAnaCalibEventSelectionEfficiency : public TObject
{
  
 public:
  
  AliFMDAnaCalibEventSelectionEfficiency();
  AliFMDAnaCalibEventSelectionEfficiency(const AliFMDAnaCalibEventSelectionEfficiency& o);
  AliFMDAnaCalibEventSelectionEfficiency& operator=(const AliFMDAnaCalibEventSelectionEfficiency& o);
  
  void    Init();
  Bool_t  IsFolder() const { return kTRUE; }
  void    Browse(TBrowser* b);
  void    SetCorrection(TH1F* hCorrection);
  Float_t GetCorrection(Int_t vtxbin);
  void    SetCorrection(TString trig, Int_t vtxbin, Char_t ring, TH2F* hCorrection);
  TH2F*   GetCorrection(TString trig, Int_t vtxbin, Char_t ring);
  Float_t GetVtxToTriggerRatio() {return fVtxEff; }
  void    SetVtxToTriggerRatio(Float_t vtxEff) {fVtxEff = vtxEff;}
 protected:
  
  TH1F fCorrection;
  TList fCorrectionList;
  Bool_t fIsInit;
  Float_t fVtxEff;
  ClassDef(AliFMDAnaCalibEventSelectionEfficiency,3);
};

#endif
