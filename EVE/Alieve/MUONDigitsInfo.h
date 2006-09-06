/* HEAD11Jul06 */
#ifndef ALIEVE_MUONDigitsInfo_H
#define ALIEVE_MUONDigitsInfo_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Access interface to the trees with digits, clusters, tracks          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TObject.h>

#include <Reve/VSD.h>

namespace Alieve {

class MUONDigitsInfo : public TObject
{

 public:

  MUONDigitsInfo() : TObject(), fDTree(0), fRTree(0), fTTree(0) {}
  virtual ~MUONDigitsInfo() {}

  MUONDigitsInfo(const MUONDigitsInfo&);            
  MUONDigitsInfo& operator=(const MUONDigitsInfo&);

  void SetDTree(TTree* tree);
  void SetRTree(TTree* tree);
  void SetTTree(TTree* tree);

  TClonesArray* GetDigits(Int_t chamber);
  TClonesArray* GetClusters(Int_t chamber);
  TClonesArray* GetTracks();

 private:

  void CreateColors();

  TTree* fDTree;          // Tree with digits
  TTree* fRTree;          // Tree with clusters
  TTree* fTTree;          // Tree with tracks
  
 protected:

 ClassDef(MUONDigitsInfo,1);

};

}

#endif
