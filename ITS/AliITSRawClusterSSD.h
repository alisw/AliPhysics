#ifndef ALIITSRAWCLUSTERSSD_H
#define ALIITSRAWCLUSTERSSD_H

#include "AliITSRawCluster.h"

////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
//  Raw Clusters for SSD                          //
//                                                //
////////////////////////////////////////////////////

class AliITSRawClusterSSD : public AliITSRawCluster {
 public:
  AliITSRawClusterSSD();
    AliITSRawClusterSSD(Float_t Prob,Int_t Sp,Int_t Sn);
    virtual ~AliITSRawClusterSSD() {// destructor
    }
    void SetMultN(Int_t m){fMultiplicityN = m;}
    void SetQErr(Float_t qe){fQErr = qe;}
    void SetSignalP(Float_t sp){fSignalP = sp;}
    void SetSignalN(Float_t sn){fSignalN = sn;}
    void SetFlag(Int_t st){fStatus = st;}
    void SetNTrack(Int_t nt){fNtracks = nt;}
    Int_t  GetStatus() const {// get status
	return fStatus;}
    void   SetStatus(Int_t status) {// set status
	fStatus=status;}
 protected:
    Int_t   fMultiplicityN;  // The number of N side strips involved 
                             // in this point calculations
    Float_t fQErr;           // Total charge error
    Float_t fSignalP;        // Signal of P side cluster
    Float_t fSignalN;        // Signal of N side cluster
    Int_t   fStatus;         // Flag status : 0 - real point
                             //               1 - ghost 
                             //               2 - EIC ? 
                             //               3 - single side 
    Int_t fNtracks;          // Number of tracks created the cluster
    ClassDef(AliITSRawClusterSSD,1)  // AliITSRawCluster class for SSD
};
#endif
