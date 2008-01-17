#ifndef ALIEVE_TOFDigitsInfo_H
#define ALIEVE_TOFDigitsInfo_H

#include <TEveVSD.h>

//#include <map>

#include <TObject.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <AliTOF.h>
#include <AliTOFGeometry.h>
#include <AliTOFDigitMap.h>

namespace Alieve {
  
class TOFDigitsInfo : public TObject, public TEveRefCnt
  {
    TOFDigitsInfo(const TOFDigitsInfo&);            // Not implemented
    TOFDigitsInfo& operator=(const TOFDigitsInfo&); // Not implemented
    
  private:

  protected:

    void        SetTOFSegmentation();

  public:
    TTree*                fTree;
    TTree*                fNewTree;
    AliTOFGeometry*       fGeom;
    AliTOFDigitMap*       fTOFdigitMap;

    TOFDigitsInfo();
    virtual ~TOFDigitsInfo();
    
    void SetTree(TTree* tree);
    void LoadDigits();

    //TClonesArray* GetDigits(Int_t nSector,
    void GetDigits(Int_t nSector, Int_t nPlate, Int_t nStrip,
		   Int_t nPadZ, Int_t nPadX,
		   Int_t indexDigit[3]);
    TClonesArray* GetDigits(Int_t nSector, Int_t nPlate, Int_t nStrip);
    TClonesArray* GetDigits(Int_t nSector);
    void GetDigits();
  
    ClassDef(TOFDigitsInfo, 1);
  }; // endclass TOFDigitsInfo
  
}
#endif
