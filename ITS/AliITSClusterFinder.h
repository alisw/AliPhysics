#ifndef ALIITSCLUSTERFINDER_H
#define ALIITSCLUSTERFINDER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */

////////////////////////////////////////////////
//  ITS Cluster Finder Class                  //
////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>

class AliITSMap;
class AliITSresponse;
class AliITSsegmentation;
class AliITSRawCluster;
class AliITS;
class AliITSdigit;
class AliITSRecPoint;

//----------------------------------------------------------------------
class AliITSClusterFinder :public TObject{
 public:
    AliITSClusterFinder();
    AliITSClusterFinder(AliITSsegmentation *seg, AliITSresponse *resp,
			TClonesArray *digits);
    virtual ~AliITSClusterFinder();
    AliITSClusterFinder(const AliITSClusterFinder &source); // copy constructor
    // assignment operator
    AliITSClusterFinder& operator=(const AliITSClusterFinder &source);
    virtual void SetResponse(AliITSresponse *response) {
	// set response
	fResponse=response;
    }
    virtual void SetSegmentation(AliITSsegmentation *segmentation) {
	// set segmentation
	fSegmentation=segmentation;
    }
    virtual void SetDigits(TClonesArray *ITSdigits) {
	// set digits
	fDigits=ITSdigits;
	fNdigits = fDigits->GetEntriesFast();
    }
    virtual AliITSdigit* GetDigit(Int_t i){
	return (AliITSdigit*) fDigits->UncheckedAt(i);
    }
    virtual TClonesArray* Digits(){
	return fDigits;
    }
    virtual Int_t   NDigits() const {
	// Get Number of Digits
	return fNdigits;
    }
    AliITSMap   *Map() {
	// map
	return fMap;
    }
    //
    virtual void AddCluster(Int_t branch, AliITSRawCluster *c);
    virtual void AddCluster(Int_t branch, AliITSRawCluster *c,
			    AliITSRecPoint &rp);
    virtual void FindRawClusters(Int_t mod=0); // Finds cluster of digits.
    // Determins if digit i has a neighbor and if so that neighor index is j.
    virtual Bool_t IsNeighbor(TObjArray *digs,Int_t i,Int_t j[]) const;
    // Given a cluster of digits, creates the nessesary RecPoint. May also
    // do some peak separation.
    virtual void CreateRecPoints(TObjArray *cluster,Int_t mod){};
    virtual void FindCluster(Int_t i, Int_t j, AliITSRawCluster *c) {
	// find cluster
    }
    virtual void Decluster(AliITSRawCluster *cluster) {
	// Decluster
    }
    virtual void SetNperMax(Int_t npermax=3) {
	// Set max. Number of cells per local cluster
	fNperMax = npermax;
    }
    virtual void SetDeclusterFlag(Int_t flag=1) {
	// Decluster ?
	fDeclusterFlag =flag;
    }
    virtual void SetClusterSize(Int_t clsize=3) {
	// Set max. cluster size ; bigger clusters will be rejected
	fClusterSize = clsize;
    }
    virtual void CalibrateCOG() {
	// Self Calibration of COG 
    }
    virtual void CorrectCOG(){
	// correct COG
    }
    virtual Bool_t Centered(AliITSRawCluster *cluster) const {
	// cluster
	return kTRUE;
    }
    virtual void   SplitByLocalMaxima(AliITSRawCluster *cluster) {
	// split by local maxima
    }
    virtual void   FillCluster(AliITSRawCluster *cluster, Int_t) {
	// fiil cluster
    }
    virtual void   FillCluster(AliITSRawCluster *cluster) {
	// fill cluster
	FillCluster(cluster,1);
    }
    // set the fitting methods in the derived classes
    // data members

 
 protected:
    TClonesArray       *fDigits;       //! digits
    Int_t              fNdigits;       //! num of digits 
    AliITSresponse     *fResponse;     //! response
    AliITSsegmentation *fSegmentation; //!segmentation
    Int_t              fNRawClusters;  //! in case we split the cluster
                                       // and want to keep track of 
                                       // the cluster which was splitted
    AliITSMap          *fMap;          //! map
    Int_t              fNperMax;       //! NperMax
    Int_t              fDeclusterFlag; //! DeclusterFlag
    Int_t              fClusterSize;   //! ClusterSize
    Int_t              fNPeaks;        //! NPeaks  
  
    ClassDef(AliITSClusterFinder,2) //Class for clustering and reconstruction of space points
};
#endif
