#ifndef AliRICHClusterFinder_h
#define AliRICHClusterFinder_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"
#include "AliRICH.h"

class AliHitMap;
class TF1;
class TClonesArray;
class AliSegmentation;
class AliRICHResponse;
class TClonesArray;

class AliRICHClusterFinder : public TObject
{
public:
  AliRICHClusterFinder(AliRICH *pRICH);
  virtual ~AliRICHClusterFinder()                                                   {;}
  
  void   Exec();                                                                     //Loop on events and chambers  
  void   FindRawClusters();                                                          //Search for clusters  
  void   AddDigit2Cluster(Int_t i, Int_t j, AliRICHRawCluster &c);
  void   Decluster(AliRICHRawCluster *cluster);                                      //Decluster
  void   CalibrateCOG();                                                             //Self Calibration of COG 
  void   SinoidalFit(Float_t x, Float_t y, TF1 *func);
  void   CorrectCOG(){;}      
  Bool_t Centered(AliRICHRawCluster *cluster);
  void   SplitByLocalMaxima(AliRICHRawCluster *cluster);
  void   FillCluster(AliRICHRawCluster *cluster, Int_t flag);
  
  void   AddRawCluster(const AliRICHRawCluster c)  {c.Print("");Rich()->AddClusterOld(fChamber,c);fNRawClusters++;}
  void   FillCluster(AliRICHRawCluster *cluster)   {FillCluster(cluster,1);}
  void   SetNperMax(Int_t npermax=5)               {fNperMax = npermax;}     //Set max. Number of pads per local cluster
  void   SetDeclusterFlag(Int_t flag=1)            {fDeclusterFlag =flag;}   //Decluster ?
  void   SetClusterSize(Int_t clsize=10)           {fClusterSize = clsize;}  //Max. cluster size; bigger clusters will be rejected
  
  AliRICH * Rich() {return fRICH;}
//protected:
  AliRICH                *fRICH;  
  AliSegmentation        *fSegmentation;                 //Segmentation model
  AliRICHResponse*        fResponse;                     //Response model
  AliHitMap              *fHitMap;                       //Hit Map with digit positions
  TF1*                    fCogCorr;                      //Correction for center of gravity
  TClonesArray*           fDigits;                       //List of digits
  Int_t                   fNdigits;                      //Number of digits
  Int_t                   fChamber;                      //Chamber number
  Int_t                   fNRawClusters;                 //Number of raw clusters
  Int_t                   fNperMax;                      //Number of pad hits per local maximum
  Int_t                   fDeclusterFlag;                //Split clusters flag
  Int_t                   fClusterSize;                  //Size of cluster 
  Int_t                   fNPeaks;                       //Number of maxima in the cluster
  ClassDef(AliRICHClusterFinder,0) //Class for clustering and reconstruction of space points    
};
#endif
