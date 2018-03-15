/// \class AliAODTrdTracklet

#ifndef ALIAODTRDTRACKLET_H
#define ALIAODTRDTRACKLET_H

#include "AliVTrdTracklet.h"

class AliAODTrdTracklet : public AliVTrdTracklet {
 public:

  AliAODTrdTracklet();
  AliAODTrdTracklet(const AliVTrdTracklet &rhs);
  AliAODTrdTracklet(UInt_t trackletWord, Short_t hcid, Int_t label = -1);
  virtual ~AliAODTrdTracklet() {};
  AliAODTrdTracklet(const AliAODTrdTracklet& track);
  AliAODTrdTracklet& operator=(const AliAODTrdTracklet& track);
  virtual void Copy(TObject &obj) const;

  void SetTrackletWord(UInt_t trklWord) { fTrackletWord = trklWord; }
  void SetHCId(Short_t hcid) { fHCId = hcid; }
  void SetLabel(Int_t label) { fLabel = label; }

  virtual UInt_t GetTrackletWord() const { return fTrackletWord; }
  virtual Int_t  GetBinY()  const;
  virtual Int_t  GetBinDy() const;
  virtual Int_t  GetBinZ()  const { return ((fTrackletWord >> 20) & 0xf);  }
  virtual Int_t  GetPID()   const { return ((fTrackletWord >> 24) & 0xff); }

  virtual Int_t GetHCId() const { return fHCId; }
  virtual Int_t GetDetector() const { return fHCId / 2; }

  virtual Int_t GetLabel() const { return fLabel; }

 protected:
  Short_t fHCId;		///< half-chamber ID
  UInt_t fTrackletWord;		///< tracklet word (as from FEE)
				// pppp : pppp : zzzz : dddd : dddy : yyyy : yyyy : yyyy
  Int_t  fLabel;		///< MC label

  ClassDef(AliAODTrdTracklet,1)
};

#endif
