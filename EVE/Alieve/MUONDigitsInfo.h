#ifndef ALIEVE_MUONDigitsInfo_H
#define ALIEVE_MUONDigitsInfo_H

#include <TObject.h>

#include <Reve/VSD.h>

namespace Alieve {

class MUONDigitsInfo : public TObject
{
  MUONDigitsInfo(const MUONDigitsInfo&);            // Not implemented
  MUONDigitsInfo& operator=(const MUONDigitsInfo&); // Not implemented

 public:

  MUONDigitsInfo() : TObject(), fDTree(0), fRTree(0), fTTree(0) {}
  virtual ~MUONDigitsInfo() {}

  void SetDTree(TTree* tree);
  void SetRTree(TTree* tree);
  void SetTTree(TTree* tree);

  TClonesArray* GetDigits(Int_t chamber);
  TClonesArray* GetClusters(Int_t chamber);
  TClonesArray* GetTracks();

  TTree* fDTree;
  TTree* fRTree;
  TTree* fTTree;
  
 private:

  void CreateColors();

 protected:

 ClassDef(MUONDigitsInfo,1);

};  // end class MUONDigitsInfo

}

#endif
