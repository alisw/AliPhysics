#ifndef ALIFMDANACALIBSHARINGEFFICIENCY_H
#define ALIFMDANACALIBSHARINGEFFICIENCY_H

#include <TObject.h>
#include <TAxis.h>
#include <TObjArray.h>
class TH2F;
class TH1F;
class TBrowser;

/**
 * @ingroup FMD_ana
 * @brief Do the background correction
 * 
 */
class AliFMDAnaCalibSharingEfficiency : public TObject
{
  
 public:
  
  AliFMDAnaCalibSharingEfficiency();
  AliFMDAnaCalibSharingEfficiency(const AliFMDAnaCalibSharingEfficiency& o);
  AliFMDAnaCalibSharingEfficiency& operator=(const AliFMDAnaCalibSharingEfficiency& o);

  void    SetSharingEff(Int_t det, Char_t ring, Int_t vtxbin, TH1F* hCorrection);
  TH1F*   GetSharingEff(Int_t det, Char_t ring, Int_t vtxbin);
  void    SetSharingEffTrVtx(Int_t det, Char_t ring, Int_t vtxbin, TH1F* hCorrection);
  TH1F*   GetSharingEffTrVtx(Int_t det, Char_t ring, Int_t vtxbin);
  void    Init();
  Bool_t  IsFolder() const { return kTRUE; }
  void    Browse(TBrowser* b);
 protected:

  TObjArray* GetRingArray(Int_t det, Char_t ring);
  TObjArray* GetRingArrayTrVtx(Int_t det, Char_t ring);
  
  TObjArray  fArray;
  TObjArray  fArrayTrVtx;
  Bool_t     fIsInit;

  ClassDef(AliFMDAnaCalibSharingEfficiency,1);
};

#endif
// Local Variables:
//   mode: C++
// End:
