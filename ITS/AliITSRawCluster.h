#ifndef ALIITSRAWCLUSTER_H
#define ALIITSRAWCLUSTER_H

////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
////////////////////////////////////////////////////

#include <TObject.h>
// this class is subject to changes !!! - info used for declustering
// and  eventual debugging

class AliITSRawCluster : public TObject {  
 public:
  AliITSRawCluster();
    virtual ~AliITSRawCluster() {// destructor
    }
    virtual Bool_t IsSortable() const {// is sortable
	return kTRUE;}
    virtual void SetMultiplicity(Int_t m) {fMultiplicity=m;}
 protected:
    Int_t       fMultiplicity;        // cluster multiplicity
  
    ClassDef(AliITSRawCluster,1)  // AliITSRawCluster class
};

#endif
