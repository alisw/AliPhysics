#ifndef AliMUONClusterFinder_H
#define AliMUONClusterFinder_H
////////////////////////////////////////////////
//  MUON Cluster Finder Class                 //
////////////////////////////////////////////////
#include "AliMUONHitMap.h"
#include "TF1.h"
class AliMUONClusterFinder :
 public TObject
{
public:
    TClonesArray*           fDigits;
    Int_t                   fNdigits;
protected:
    AliMUONsegmentation*    fSegmentation;
    AliMUONresponse*        fResponse;
    TClonesArray*           fRawClusters;
    Int_t                   fChamber;
    Int_t                   fNRawClusters;
    AliMUONHitMapA1*        fHitMap;
    TF1*                    fCogCorr;
    Int_t                   fNperMax;
    Int_t                   fDeclusterFlag;
    Int_t                   fClusterSize;
    Int_t                   fNPeaks; 
 public:
    AliMUONClusterFinder
	(AliMUONsegmentation *segmentation,
	 AliMUONresponse *response, TClonesArray *digits, Int_t chamber);
    AliMUONClusterFinder();
    ~AliMUONClusterFinder(){delete fRawClusters;}
    virtual void SetSegmentation(
	AliMUONsegmentation *segmentation){
	fSegmentation=segmentation;
    }
    virtual void SetResponse(AliMUONresponse *response) {
	fResponse=response;
    }

    virtual void SetDigits(TClonesArray *MUONdigits) {
	fDigits=MUONdigits;
	fNdigits = fDigits->GetEntriesFast();
    }
    
    virtual void SetChamber(Int_t ich){
	fChamber=ich;
    }
    
    virtual void AddRawCluster(const AliMUONRawCluster);
    // Search for raw clusters
    virtual void FindRawClusters();
    virtual void  FindCluster(Int_t i, Int_t j, AliMUONRawCluster &c);
    // Decluster
    virtual void Decluster(AliMUONRawCluster *cluster);
    // Set max. Number of pads per local cluster
    virtual void SetNperMax(Int_t npermax=5) {fNperMax = npermax;}
    // Decluster ?
    virtual void SetDeclusterFlag(Int_t flag=1) {fDeclusterFlag =flag;}
    // Set max. cluster size ; bigger clusters will be rejected
    virtual void SetClusterSize(Int_t clsize=5) {fClusterSize = clsize;}
    // Self Calibration of COG 
    virtual void CalibrateCOG();
    virtual void SinoidalFit(Float_t x, Float_t y, TF1 &func);
    //
    virtual void CorrectCOG(){;}
    
    //
    virtual Bool_t Centered(AliMUONRawCluster *cluster);
    virtual void   SplitByLocalMaxima(AliMUONRawCluster *cluster);
    virtual void   FillCluster(AliMUONRawCluster *cluster, Int_t);
    virtual void   FillCluster(AliMUONRawCluster *cluster) {
	FillCluster(cluster,1);}
    TClonesArray* RawClusters(){return fRawClusters;}
    ClassDef(AliMUONClusterFinder,1) //Class for clustering and reconstruction of space points
};
#endif







