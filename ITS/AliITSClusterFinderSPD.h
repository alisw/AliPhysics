#ifndef ALIITSCLUSTERFINDERSPD_H
#define ALIITSCLUSTERFINDERSPD_H

////////////////////////////////////////////////
//  ITS Cluster Finder Class                 //
////////////////////////////////////////////////

#include "AliITSClusterFinder.h"

class AliITSMapA1;
class AliITSsegmentation;
class AliITSsegmentationSPD;
class AliITSresponse;
class AliITSresponseSPD;
class TClonesArray;

class AliITSClusterFinderSPD : public AliITSClusterFinder{
 public:
    AliITSClusterFinderSPD();
    AliITSClusterFinderSPD(AliITSsegmentation *segmentation,
                           AliITSresponse *res);
    AliITSClusterFinderSPD(AliITSsegmentation *segmentation,
			   TClonesArray *digits,TClonesArray *recpoints);
    virtual ~AliITSClusterFinderSPD(){}// destructor

    virtual AliITSresponseSPD* GetResp()const{
        return (AliITSresponseSPD*) GetResponse();}//Return Response
    //Returns fSegmentation
    virtual AliITSsegmentationSPD* GetSeg()const{
        return (AliITSsegmentationSPD*)GetSegmentation();}  
    virtual void SetDx(Double_t dx=1.) {fDx=dx;}// set dx
    virtual void SetDz(Double_t dz=0.) {fDz=dz;}// set dz
    // Search for clusters
    virtual void FindRawClusters(Int_t module); 
    void  DigitToPoint(Int_t nclus, Double_t *xcenter, Double_t *zcenter,
                       Double_t *errxcenter,Double_t *errzcenter,
                       Int_t *tr1clus, Int_t *tr2clus, Int_t *tr3clus);
    void  ClusterFinder(Int_t ndigits,Int_t digx[],Int_t digz[],
                        Int_t digtr1[],Int_t digtr2[],Int_t digtr3[],
                        Int_t digtr4[],
                        Int_t &nclus,
                        Double_t xcenter[],Double_t zcenter[],
                        Double_t errxcenter[],Double_t errzcenter[],  
                        Int_t tr1clus[],Int_t tr2clus[], Int_t tr3clus[]);
 protected:
    // copy constructor
    AliITSClusterFinderSPD(const AliITSClusterFinderSPD &source);
    // assignment operator
    AliITSClusterFinderSPD& operator=(const AliITSClusterFinderSPD &source);

    Double_t             fDz;            // dz
    Double_t             fDx;            // dx
    Int_t               fMinNCells;     // min num of cells in the cluster
  
    ClassDef(AliITSClusterFinderSPD,2)  // SPD clustering
};
#endif
