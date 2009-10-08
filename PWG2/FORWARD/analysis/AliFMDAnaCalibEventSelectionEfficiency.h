#ifndef ALIFMDANACALIBEVENTSELECTIONEFFICIENCY_H
#define ALIFMDANACALIBEVENTSELECTIONEFFICIENCY_H

#include <TObject.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TList.h>

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
 protected:
  
  TH1F fCorrection;
  Bool_t fIsInit;
  ClassDef(AliFMDAnaCalibEventSelectionEfficiency,1);
};

#endif
