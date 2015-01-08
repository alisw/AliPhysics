#ifndef AliPHOSClusterizerv2_H
#define AliPHOSClusterizerv2_H

// --- AliRoot header files ---
#include "AliPHOSClusterizerv1.h"

class AliPHOSClusterizerv2 : public AliPHOSClusterizerv1 {

public:

  AliPHOSClusterizerv2();
  AliPHOSClusterizerv2(AliPHOSGeometry *geom);
  ~AliPHOSClusterizerv2() {}
  
  Int_t AreNeighbours(AliPHOSDigit* d1, AliPHOSDigit* d2) const ;
  void GetNumberOfClustersFound(int* numb ) const;

  virtual void Digits2Clusters(Option_t* option);

private:
  AliPHOSClusterizerv2(const AliPHOSClusterizerv2 & clu) ;
  AliPHOSClusterizerv2 & operator = (const AliPHOSClusterizerv2 &rvalue);

  ClassDef(AliPHOSClusterizerv2,2)
};

#endif
