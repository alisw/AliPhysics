#ifndef ALIITSCLUSTERFINDERSSD_H
#define ALIITSCLUSTERFINDERSSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <TObject.h>
#include <TArrayS.h>
#include <TClonesArray.h>
#include <TRandom.h>
#include <TMath.h>

#include "AliITS.h"
#include "AliITSclusterSSD.h"
#include "AliITSpackageSSD.h"
#include "AliITSClusterFinder.h"


class AliITSClusterFinderSSD: public AliITSClusterFinder 
{
    
public:       
    
  AliITSClusterFinderSSD(AliITSsegmentation *seg, TClonesArray *digits, TClonesArray *recp);
		
  virtual ~AliITSClusterFinderSSD();


  void FindRawClusters();

    
  void SetAlpha1(Float_t a) {falpha1 =a;}
  void SetAlpha2(Float_t a) {falpha2 =a;}
  void SetAlpha3(Float_t a) {falpha3 =a;}



 protected:
   
  void      InitReconstruction();
  Bool_t    CreateNewRecPoint(Float_t P, Float_t dP, Float_t N, Float_t dN,
                              Float_t Sig,Float_t dSig,
                              AliITSclusterSSD *clusterP, AliITSclusterSSD *clusterN,
                              Stat_t prob);
  Bool_t   CreateNewRecPoint(AliITSclusterSSD *clusterP, AliITSclusterSSD *clusterN, Stat_t prob);
  
  AliITSclusterSSD* GetPSideCluster(Int_t idx);
  AliITSclusterSSD* GetNSideCluster(Int_t idx);
  AliITSclusterSSD* GetCluster(Int_t idx, Bool_t side);


  void      FindNeighbouringDigits();
  void      SeparateOverlappedClusters();
  void      SplitCluster(TArrayI *list, Int_t nsplits, Int_t index, Bool_t side);
  Int_t     SortDigitsP(Int_t start, Int_t end);
  Int_t     SortDigitsN(Int_t start, Int_t end);
  void      FillDigitsIndex();
  void      SortDigits();
  void      FillClIndexArrays(Int_t* arrayP, Int_t *arrayN);
  void      SortClusters(Int_t* arrayP, Int_t *arrayN);
  Int_t     SortClustersP(Int_t start, Int_t end,Int_t *array);
  Int_t     SortClustersN(Int_t start, Int_t end,Int_t *array);
  void      ConsumeClusters();
  void      ClustersToPackages();
  void      PackagesToPoints();
  void      ReconstructNotConsumedClusters();
  Bool_t    Strip2Local( Float_t stripP, Float_t stripN, Float_t &Z,Float_t &X);
  Float_t   GetClusterZ(AliITSclusterSSD* clust);
  Bool_t    IsCrossing(AliITSclusterSSD* p, AliITSclusterSSD* n);
  //returns index of best combination in "comb"
  Int_t     GetBestComb(Int_t** comb,Int_t Ncomb, Int_t Ncl, AliITSpackageSSD * pkg);

  //get point that have best signal ratio
  void      GetBestMatchingPoint(Int_t & ip, Int_t & in, AliITSpackageSSD* pkg );					
  
  //calculates Distance To Perfect Matching Line
  Float_t   DistToPML(Float_t psig, Float_t nsig){ return  (TMath::Abs( (7.0*nsig - 8.0*psig )/10.630146) );}
  
  
  Int_t     GetDiff(Float_t *retx, Float_t *rety) {return 0;}
  
  void      CalcStepFactor(Float_t Psteo, Float_t Nsteo );
  
/*************************************************/
/**  methods for resolving packages           ****/
/*************************************************/ 
//names may not be meaningful for all, see implementations for descriptions

  void      ResolveSimplePackage(AliITSpackageSSD *pkg);
  void      ResolvePackageWithOnePSideCluster(AliITSpackageSSD *pkg);
  void      ResolvePackageWithOneNSideCluster(AliITSpackageSSD *pkg);
  void      ResolveTwoForTwoPackage(AliITSpackageSSD *pkg);

  void      ResolveClusterWithOneCross(AliITSpackageSSD *pkg,
                                  Int_t clusterIndex, Bool_t clusterSide);

  void      ResolvePClusterWithOneCross(AliITSpackageSSD *pkg, Int_t clusterIndex);
  void      ResolveNClusterWithOneCross(AliITSpackageSSD *pkg, Int_t clusterIndex);
  Bool_t    ResolvePackageBestCombin(AliITSpackageSSD *pkg);
  void      ResolveOneBestMatchingPoint(AliITSpackageSSD *pkg);

  Bool_t GetCrossing(Float_t &x, Float_t &z);     //x, y of strips crossing
  void   GetCrossingError(Float_t&, Float_t&);    //x, y of strips crossing errors

  // Data memebers

  AliITS             *fITS;         //Pointer to AliITS object
  TClonesArray       *fDigits;      //Pointer to TClonesArray of digits

  TClonesArray       *fRecPoints;   //Pointer to TClonesArray of rec points
  
	
  TClonesArray    *fClusterP;    //
  Int_t            fNClusterP;   //Number of P side clusters in the array
		
  TClonesArray    *fClusterN;    //Number of N side clusters in the array
  Int_t            fNClusterN; 
    
  TClonesArray    *fPackages;    //packages  
  Int_t            fNPackages;
    
  TArrayI         *fDigitsIndexP;       //Digits on P side
  Int_t            fNDigitsP;           //Number of Digits on P side
		
  TArrayI         *fDigitsIndexN;       //Digits on N side
  Int_t            fNDigitsN;           //Number of Digits on N side


  Float_t          fPitch;              //Strip pitch
  Float_t          fTanP;               //Pside stereo angle tangent
  Float_t          fTanN;               //Nside stereo angle tangent

/*************************************************/
/**  parameters for reconstruction            ****/
/**  to be tune when slow simulation raliable ****/
/*************************************************/ 
  
  Float_t falpha1; 
  Float_t falpha2;
  Float_t falpha3;

    
  static const Bool_t SIDEP=kTRUE;
  static const Bool_t SIDEN=kFALSE;
  static const Float_t PNsignalRatio = 7./8.;
  Int_t fSFF;              //forward stepping factor 
  Int_t fSFB;              //backward stepping factor 

public:
    ClassDef(AliITSClusterFinderSSD, 1) //Class for clustering and reconstruction of space points in SSDs 

};


#endif
