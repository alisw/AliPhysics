#ifndef AliPHOSClusterizerv2_H
#define AliPHOSClusterizerv2_H

// --- AliRoot header files ---
#include "AliPHOSClusterizerv1.h"

class AliPHOSClusterizerv2 : public AliPHOSClusterizerv1 {

public:

  AliPHOSClusterizerv2();
  AliPHOSClusterizerv2(const char * headerFile, const char * name = "Default", const Bool_t toSplit=kFALSE) ;
  ~AliPHOSClusterizerv2() {}
  
  Int_t AreNeighbours(AliPHOSDigit* d1, AliPHOSDigit* d2) const ;
  void GetNumberOfClustersFound(int* numb ) const;

  void Exec(Option_t* option);

  ClassDef(AliPHOSClusterizerv2,1)
};

#endif
