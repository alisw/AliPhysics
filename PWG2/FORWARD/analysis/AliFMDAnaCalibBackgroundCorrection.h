#ifndef ALIFMDANACALIBBACKGROUNDCORRECTION_H
#define ALIFMDANACALIBBACKGROUNDCORRECTION_H

#include <TObject.h>
#include <TObjArray.h>
// #include <TH2F.h>
#include <TAxis.h>
#include <TList.h>
class TH2F;
class TH1F;
class TBrowser;

/**
 * @ingroup FMD_ana
 * @brief Do the background correction
 * 
 */
class AliFMDAnaCalibBackgroundCorrection : public TObject
{
  
 public:
  
  AliFMDAnaCalibBackgroundCorrection();
  AliFMDAnaCalibBackgroundCorrection(const AliFMDAnaCalibBackgroundCorrection& o);
  AliFMDAnaCalibBackgroundCorrection& operator=(const AliFMDAnaCalibBackgroundCorrection& o);
  
  TH2F*   GetBgCorrection(Int_t det, Char_t ring, Int_t vtxbin);
  void    SetBgCorrection(Int_t det, Char_t ring, Int_t vtxbin, TH2F* hCorrection);
  TH1F*   GetDoubleHitCorrection(Int_t det, Char_t ring);
  void    SetDoubleHitCorrection(Int_t det, Char_t ring, TH1F* hCorrection);
  TH1F*   GetSPDDeadCorrection(Int_t vtxbin);
  void    SetSPDDeadCorrection(Int_t vtxbin, TH1F* hCorrection);
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
  TList      fListOfDoubleHitCorrection;
  ClassDef(AliFMDAnaCalibBackgroundCorrection,2);
};

#endif
// Local Variables:
//   mode: C++
// End:
