#ifndef PMDrecpoint_H
#define PMDrecpoint_H
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store reconstructed points  for PMD                //
//                                                     //
//-----------------------------------------------------//

#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TClonesArray.h"

class AliPMDrecpoint : public TObject
{
  
 protected:

  Float_t fClusData[7];
  /*
    fClusData[0] : Detector Number,  fClusData[1] : SuperModule Number
    fClusData[2] : Cluster x      ,  fClusData[3] : Cluster y
    fClusData[4] : Cluster adc    ,  fClusData[5] : Cluster Cells
    fClusData[6] : Cluster radius
  */

 public:
  AliPMDrecpoint();
  AliPMDrecpoint(Float_t * /* clusdata */);
  AliPMDrecpoint(AliPMDrecpoint *pmdrecpoint) {*this = *pmdrecpoint;}
  
  virtual ~AliPMDrecpoint();

  Float_t GetDetector() const;
  Float_t GetSMNumber() const;
  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusRadius() const;
  
  ClassDef(AliPMDrecpoint,1)
};

#endif
