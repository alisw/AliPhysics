#ifndef AliESDAD_H
#define AliESDAD_H

#include <TObject.h>
// Send comments to:
// Mario Rodriguez <mrodriguez@fis.cinvestav.mx>

class AliESDAD : public TObject
{

 public:
  AliESDAD();
  AliESDAD(const AliESDAD&);
  AliESDAD(Bool_t *ADBitCell); 
  virtual ~AliESDAD() {};
  virtual void Copy(TObject &) const;

 // We define the "setters" for AD
	// fake bit pattern, but enought to MC studies
  void SetADBitCell(Bool_t ADBitCell[16]){for (Int_t i=0;i<16;i++){fADCellID[i]=ADBitCell[i];}}

 // Getters  	
  Bool_t GetADCell(Int_t i) const;
  AliESDAD &operator=(const AliESDAD& source);
  
 protected:

  Bool_t	fADCellID[16];  // Array with the AD's bitpattern

  ClassDef(AliESDAD, 1)

};

#endif
