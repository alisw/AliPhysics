#ifndef ALIPHOSCLUSTERIZER_H
#define ALIPHOSCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Algorithme class for the clusterization   //
//  interface class                           //
//  Version SUBATECH                          //
//  Author Yves Schutz     SUBATECH           //
//                                            //  
//   pABC                                     //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h" 
#include "TClonesArray.h"

// --- Standard library ---

// --- AliRoot header files ---


typedef TClonesArray    RecPointsList ; // a cluster has a variable size (see ROOT FAQ)  
typedef TClonesArray    DigitsList ; //for digits saved on disk

class AliPHOSClusterizer : public TObject {

public:

  AliPHOSClusterizer() ;          // ctor            
  virtual ~AliPHOSClusterizer() ; // dtor

  virtual Float_t Calibrate(Int_t Amp) = 0 ; 
  virtual void  GetNumberOfClustersFound(Int_t * numb) = 0 ; 
  virtual void  MakeClusters(const DigitsList * dl, RecPointsList * emccl, RecPointsList * ppsdl) = 0 ; 

  ClassDef(AliPHOSClusterizer,1)  // clusterization interface, version 1

} ;

#endif // AliPHOSCLUSTERIZER_H
