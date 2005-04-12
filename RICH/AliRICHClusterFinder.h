#ifndef AliRICHClusterFinder_h
#define AliRICHClusterFinder_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TTask.h"

#include "AliRICH.h"
class AliHitMap;

class AliRICHClusterFinder : public TTask
{
public:    
           AliRICHClusterFinder(AliRunLoader *pRunLoader);
  virtual ~AliRICHClusterFinder()                                          {;}
  
  AliRICH        *R()                                              {return fRICH;}   //returns pointer to RICH  
  void            Exec(const Option_t *option="");                                   //loop on events and chambers  
  void            FindClusters(Int_t iChamber);                                      //find all clusters for a given chamber
  void            FindClusterContribs(AliRICHCluster *pCluster);                     //find CFM for the current cluster
  void            FormRawCluster(Int_t i, Int_t j);                                  //form a raw cluster
  void            FindLocalMaxima();                                                 //find local maxima in a cluster
  void            FitCoG();                                                          //evaluate the CoG as the best 
  void            WriteRawCluster();                                                 //write in the list of cluster  
  void            WriteResolvedCluster();                                            //write in the list of cluster  
  AliRICHCluster *GetRawCluster()                           {return &fRawCluster;}   //returns pointer to the current raw cluster
protected:
  AliRICH                *fRICH;                         //pointer to RICH
  AliHitMap              *fDigitMap;                     //map of digits positions
  AliRICHCluster         fRawCluster;                    //current raw cluster before deconvolution
  AliRICHCluster         fResolvedCluster;               //current cluster after deconvolution
  Int_t                  fNlocals;                       //number of local maxima
  Double_t               fLocalX[100],fLocalY[100];      //list of locals X,Y
  Double_t               fLocalQ[100];                   //list of locals charge Q
  Int_t                  fLocalC[100];                   //list of locals CombiPid
  ClassDef(AliRICHClusterFinder,0) //finds raw clusters, trasfers them to resolved clusters through declustering.
};
#endif
