#ifndef ALIEVE_MUONDigitsInfo_H
#define ALIEVE_MUONDigitsInfo_H

#include <TObject.h>

#include <Reve/VSD.h>

namespace Alieve {

class MUONDigitsInfo : public TObject
{

 public:

  MUONDigitsInfo(const Text_t* /*n*/="MUONDigitsInfo", const Text_t* /*t*/=0) :
      TObject()
      { Init(); } 
 
  virtual ~MUONDigitsInfo();

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

  void Init();
  void CreateColors();

 protected:

 ClassDef(MUONDigitsInfo,1);

};  // end class MUONDigitsInfo

}

#endif
