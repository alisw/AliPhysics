#ifndef AliPHOSClusterizerv2_H
#define AliPHOSClusterizerv2_H

// --- AliRoot header files ---
#include "AliPHOSClusterizerv1.h"

class AliPHOSClusterizerv2 : public AliPHOSClusterizerv1 {

public:

  AliPHOSClusterizerv2();
  AliPHOSClusterizerv2(const char * headerFile, const char * name = "Default") ;
  AliPHOSClusterizerv2(const AliPHOSClusterizerv2 & clu) ;
  ~AliPHOSClusterizerv2() {}
  
  Int_t AreNeighbours(AliPHOSDigit* d1, AliPHOSDigit* d2) const ;
  void GetNumberOfClustersFound(int* numb ) const;

  void Exec(Option_t* option);
  AliPHOSClusterizerv2 & operator = (const AliPHOSClusterizerv2 & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }

  ClassDef(AliPHOSClusterizerv2,1)
};

#endif
