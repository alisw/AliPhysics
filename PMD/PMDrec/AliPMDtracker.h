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
class AliPMDrecdata;
class AliPMDrechit;
class AliPMDUtility;
class AliPMDDiscriminator;

class AliESDEvent;

class AliPMDtracker:public TObject
{

 public:

  AliPMDtracker();
  AliPMDtracker(const AliPMDtracker &tracker);  // copy constructor
  AliPMDtracker &operator=(const AliPMDtracker &tracker); // assignment op

  virtual ~AliPMDtracker();

  void LoadClusters(TTree *treein);
  void Clusters2Tracks(AliESDEvent *event);
  void AssignTrPidToCluster(Int_t nentry, Int_t *itra, Int_t *ipid,
			    Float_t *cadc, Int_t &trackno, Int_t &trackpid);

  void SetVertex(Double_t vtx[3], Double_t evtx[3]);
  void ResetClusters();

 protected:

  TTree        *fTreeR;     // Reconstructed points
  TClonesArray *fRecpoints; // List of reconstructed points
  TClonesArray *fRechits;   // List of cells associated with a cluster
  TObjArray    *fPMDcontin;
  TObjArray    *fPMDcontout;

  AliPMDUtility       *fPMDutil;
  AliPMDrecpoint1     *fPMDrecpoint;
  AliPMDrecdata       *fPMDclin;
  AliPMDclupid        *fPMDclout;

  Double_t fXvertex;        // X-vertex position
  Double_t fYvertex;        // Y-vertex position
  Double_t fZvertex;        // Z-vertex position
  Double_t fSigmaX;         // X-vertex error
  Double_t fSigmaY;         // Y-vertex error
  Double_t fSigmaZ;         // Z-vertex error

  ClassDef(AliPMDtracker,5) // To run PMD clustering
};
#endif

