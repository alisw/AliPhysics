#ifndef ALIPMDCLUSTER_H
#define ALIPMDCLUSTER_H
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cluster informations for PMD                 //
//                                                     //
//-----------------------------------------------------//

#include "Rtypes.h"
#include "TObject.h"
class TClonesArray;

class AliPMDcluster : public TObject
{
 public:
  AliPMDcluster();
  AliPMDcluster(Float_t * /* clusdata */);
  AliPMDcluster(AliPMDcluster *pmdcluster) {*this = *pmdcluster;}
  AliPMDcluster (const AliPMDcluster &pmdcluster);  // copy constructor
  AliPMDcluster &operator=(const AliPMDcluster &pmdcluster); // assignment op
  
  virtual ~AliPMDcluster();

  Float_t GetClusX() const;
  Float_t GetClusY() const;
  Float_t GetClusADC() const;
  Float_t GetClusCells() const;
  Float_t GetClusRadius() const;

 protected:

  Float_t fClusData[5];  // Array containing cluster information
  /*
    fClusData[0] : Cluster x      ,  fClusData[1] : Cluster y
    fClusData[2] : Cluster adc    ,  fClusData[3] : Cluster Cells
    fClusData[4] : Cluster radius
  */
  
  ClassDef(AliPMDcluster,1) // Keep Cluster information
};

#endif
