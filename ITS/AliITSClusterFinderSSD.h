#ifndef ALIITSCLUSTERFINDERSSD_H
#define ALIITSCLUSTERFINDERSSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#include <TMath.h>
#include "AliITSClusterFinder.h"
//#include "AliITSsegmentationSSD.h"

class TArrayI;
class AliITSclusterSSD;
class AliITSpackageSSD;
class AliITSsegmentation;
class AliITSsegmentationSSD;
class AliITSresponse;
class AliITSresponseSSD;

class AliITSClusterFinderSSD: public AliITSClusterFinder{
  public:
    AliITSClusterFinderSSD();
    AliITSClusterFinderSSD(AliITSsegmentation *seg,AliITSresponse *resp);
    AliITSClusterFinderSSD(AliITSsegmentation *seg, TClonesArray *digits);
    virtual ~AliITSClusterFinderSSD();
    void FindRawClusters(Int_t module);

  protected:
    virtual AliITSresponseSSD* GetResp()const{
        return (AliITSresponseSSD*) GetResponse();}//Return Response
    //Returns fSegmentation
    virtual AliITSsegmentationSSD* GetSeg()const{
        return (AliITSsegmentationSSD*)GetSegmentation();} 
    void      InitReconstruction();
    Bool_t    CreateNewRecPoint(Double_t P,Double_t dP,Double_t N,Double_t dN,
                                Double_t Sig,Double_t dSig,
                                AliITSclusterSSD *clusterP,
                                AliITSclusterSSD *clusterN,Stat_t prob);
    AliITSclusterSSD* GetPSideCluster(Int_t idx);
    AliITSclusterSSD* GetNSideCluster(Int_t idx);
    AliITSclusterSSD* GetCluster(Int_t idx, Bool_t side){
    return (side) ? GetPSideCluster(idx) : GetNSideCluster(idx);};
    void   FindNeighbouringDigits();
    void   SeparateOverlappedClusters();
    void   SplitCluster(TArrayI *list,Int_t nsplits,Int_t indx,Bool_t side);
    Int_t  SortDigitsP(Int_t start, Int_t end);
    Int_t  SortDigitsN(Int_t start, Int_t end);
    void   FillDigitsIndex();
    void   SortDigits();
    void   FillClIndexArrays(Int_t* arrayP, Int_t *arrayN) const;
    void   SortClusters(Int_t* arrayP, Int_t *arrayN);
    Int_t  SortClustersP(Int_t start, Int_t end,Int_t *array);
    Int_t  SortClustersN(Int_t start, Int_t end,Int_t *array);
    void   ClustersToPackages();
    Int_t  GetDiff(Double_t */*retx*/, Double_t */*rety*/) const {return 0;}
    void   CalcStepFactor(Double_t Psteo, Double_t Nsteo );
    Bool_t GetCrossing(Double_t &x, Double_t &z); //x, y of strips crossing
    //x, y of strips crossing err.
    void   GetCrossingError(Double_t& dp, Double_t& dn);

    // Data memebers
    AliITS          *fITS;           //!Pointer to AliITS object
    TClonesArray    *fClusterP;      //!
    Int_t            fNClusterP;     //!Number of P side clusters in the array
    TClonesArray    *fClusterN;      //!Number of N side clusters in the array
    Int_t            fNClusterN;     //!
    TClonesArray    *fPackages;      //!packages  
    Int_t            fNPackages;     //!
    TArrayI         *fDigitsIndexP;  //!Digits on P side
    Int_t            fNDigitsP;      //!Number of Digits on P side
    TArrayI         *fDigitsIndexN;  //!Digits on N side
    Int_t            fNDigitsN;      //!Number of Digits on N side

    Double_t          fPitch;         //!Strip pitch
    Double_t          fTanP;          //!Pside stereo angle tangent
    Double_t          fTanN;          //!Nside stereo angle tangent

/*************************************************/
/**  parameters for reconstruction            ****/
/**  to be tune when slow simulation raliable ****/
/*************************************************/

    //Double_t fAlpha1;         //!
    //Double_t fAlpha2;         //!
    //Double_t fAlpha3;         //!
    Double_t fPNsignalRatio;  //!
    
    static const Bool_t fgkSIDEP;  //!
    static const Bool_t fgkSIDEN;  //!

    Int_t fSFF;              //!forward stepping factor 
    Int_t fSFB;              //!backward stepping factor 

    ClassDef(AliITSClusterFinderSSD,1) //Class for clustering and reconstruction of space points in SSDs 

};


#endif
