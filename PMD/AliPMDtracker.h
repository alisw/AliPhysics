#ifndef ALIPMDTRACKER_H
#define ALIPMDTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//  This reads the file PMD.digits.root(TreeD),        //
//  calls the Clustering algorithm and stores the      //
//  clustering output in PMD.RecPoints.root(TreeR)     // 
//                                                     //
//-----------------------------------------------------//

class TClonesArray;
class TObjArray;
class TTree;

class AliPMDcluster;
class AliPMDclupid;
class AliPMDrecpoint1;
class AliPMDUtility;
class AliPMDDiscriminator;

class AliESD;

class AliPMDtracker:public TObject
{

 public:

  AliPMDtracker();
  virtual ~AliPMDtracker();

  void LoadClusters(TTree *treein);
  void Clusters2Tracks(AliESD *event);
  void SetVertex(Double_t vtx[3], Double_t evtx[3]);
  void ResetClusters();

 protected:

  TTree        *fTreeR;     // Reconstructed points
  TClonesArray *fRecpoints; // List of reconstructed points
  TObjArray    *fPMDcontin;
  TObjArray    *fPMDcontout;

  AliPMDDiscriminator *fPMDdiscriminator;
  AliPMDUtility       *fPMDutil;
  AliPMDrecpoint1     *fPMDrecpoint;
  AliPMDcluster       *fPMDclin;
  AliPMDclupid        *fPMDclout;

  Double_t fXvertex;        // X-vertex position
  Double_t fYvertex;        // Y-vertex position
  Double_t fZvertex;        // Z-vertex position
  Double_t fSigmaX;         // X-vertex error
  Double_t fSigmaY;         // Y-vertex error
  Double_t fSigmaZ;         // Z-vertex error

  ClassDef(AliPMDtracker,2) // To run PMD clustering
};
#endif

