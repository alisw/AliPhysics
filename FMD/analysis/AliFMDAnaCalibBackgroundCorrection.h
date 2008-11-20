#ifndef ALIFMDANACALIBBACKGROUNDCORRECTION_H
#define ALIFMDANACALIBBACKGROUNDCORRECTION_H

#include <TObject.h>
#include <TObjArray.h>
// #include <TH2F.h>
#include <TAxis.h>
class TH2F;
class TBrowser;

class AliFMDAnaCalibBackgroundCorrection : public TObject
{
  
 public:
  
  AliFMDAnaCalibBackgroundCorrection();
  AliFMDAnaCalibBackgroundCorrection(const AliFMDAnaCalibBackgroundCorrection& o);
  AliFMDAnaCalibBackgroundCorrection& operator=(const AliFMDAnaCalibBackgroundCorrection& o);
  
  TH2F*   GetBgCorrection(Int_t det, Char_t ring, Int_t vtxbin);
  void    SetBgCorrection(Int_t det, Char_t ring, Int_t vtxbin, TH2F* hCorrection);
  void    SetRefAxis(TAxis* axis);
  Int_t   GetNvtxBins();
  Float_t GetVtxCutZ();
  void    Init();
  Bool_t  IsFolder() const { return kTRUE; }
  void    Browse(TBrowser* b);
 protected:
  
  TObjArray* GetRingArray(Int_t det, Char_t ring);
  TObjArray  fArray;
  TAxis      fAxis;
  Bool_t     fIsInit;
  
  ClassDef(AliFMDAnaCalibBackgroundCorrection,1);
};

#endif
