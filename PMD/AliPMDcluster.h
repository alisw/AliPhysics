#ifndef PMDcluster_H
#define PMDcluster_H
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cluster informations for PMD                 //
//                                                     //
//-----------------------------------------------------//

#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TClonesArray.h"

class AliPMDcluster : public TObject
{
  
 protected:

  Float_t fClusData[5];
  /*
    fClusData[0] : Cluster x      ,  fClusData[1] : Cluster y
    fClusData[2] : Cluster adc    ,  fClusData[3] : Cluster Cells
    fClusData[4] : Cluster radius
  */

 public:
  AliPMDcluster();
  AliPMDcluster(Float_t * /* clusdata */);
  AliPMDcluster(AliPMDcluster *pmdcluster) {*this = *pmdcluster;}
  
  virtual ~AliPMDcluster();

  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusRadius() const;
  
  ClassDef(AliPMDcluster,1)
};

#endif
