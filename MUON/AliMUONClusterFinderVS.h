#ifndef ALIMUONCLUSTERFINDERVS_H
#define ALIMUONCLUSTERFINDERVS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONClusterFinderVS
/// \brief Class for clustering and reconstruction of space points

#include <TObject.h>

class TClonesArray;

class AliMUONClusterInput;
class AliMUONDigitMapA1;
class AliMUONGeometrySegmentation;
class AliMUONRawCluster;
class AliMUONDigit;


class AliMUONClusterFinderVS : public TObject 
{
 public:
    AliMUONClusterFinderVS();
    virtual ~AliMUONClusterFinderVS();
// Decluster ?
    virtual void SetDeclusterFlag(Int_t flag=1) {fDeclusterFlag =flag;}
// Set max. cluster size ; bigger clusters will deconvoluted
    virtual void SetClusterSize(Int_t clsize=5) {fClusterSize = clsize;}
// Set max. number of pads per local cluster
    virtual void SetNperMax(Int_t npermax=5) {fNperMax = npermax;}
// Search for raw clusters
    virtual void  FindRawClusters();
// Find cluster    
    virtual void  FindCluster(Int_t i, Int_t j, Int_t cath, AliMUONRawCluster &c);
// Decluster
    virtual void Decluster(AliMUONRawCluster *cluster);
//  Perform split by local maxima  
    virtual void   SplitByLocalMaxima(AliMUONRawCluster *cluster);
    virtual void   FindLocalMaxima(AliMUONRawCluster *cluster);
    virtual void   Split(AliMUONRawCluster * cluster);
//  Perform Double Mathieson Fit
    Bool_t  DoubleMathiesonFit(AliMUONRawCluster *c, Int_t cath);
    Float_t CombiDoubleMathiesonFit(AliMUONRawCluster *c);
    Float_t SingleMathiesonFit(AliMUONRawCluster *c, Int_t cath);
    Float_t CombiSingleMathiesonFit(AliMUONRawCluster *c);    
//  Build up full cluster information    
    virtual void   FillCluster(AliMUONRawCluster *cluster, Int_t flag, Int_t cath);
    virtual void   FillCluster(AliMUONRawCluster *cluster, Int_t cath);
    virtual void   FillCluster(AliMUONRawCluster *cluster) {FillCluster(cluster,1,0);}
// Add a new raw cluster    
    virtual void AddRawCluster(AliMUONRawCluster& cluster);
//  Set tracks for debugging    
    virtual void SetTracks(Int_t t1, Int_t t2) {fTrack[0]=t1; fTrack[1]=t2;}
    void SetGhostChi2Cut(Float_t cut) {fGhostChi2Cut = cut;}
// get raw cluster pointer 
    TClonesArray*  GetRawClusters() {return fRawClusters;}
// reset raw clusters
    void ResetRawClusters();

 protected:

    AliMUONClusterInput*          fInput;        //!< AliMUONClusterInput instance
    AliMUONDigitMapA1*            fDigitMap[2];  ///< Hit Maps for cathode 1 and 2
    AliMUONGeometrySegmentation*  fSeg2[2];      ///< New Segmentations for cathode 1 and 2

    
// Configuration    
    Int_t                   fDeclusterFlag;      ///< flag for declusterin
    Int_t                   fClusterSize;        ///< cluster size 
    Int_t                   fNperMax;            ///< Maximum number of pads per peak
    Float_t                 fGhostChi2Cut;       ///< \brief Cut in charge matching chi2
	                                         /// (2 degrees of freedom)
                                                 /// Used by ghost removal
    // Current decluster result    
    Int_t                   fMul[2];             ///< current multiplicity
    Int_t                   fNPeaks;             ///< number of local maxima
    Int_t                   fNRawClusters;       ///< Number of Raw Clusters
    TClonesArray*           fRawClusters;        ///< array of cluster per ch.

    // Local data store    
    AliMUONDigit*           fDig[100][2];        ///< current list of digits 
    Int_t                   fIx[100][2];         ///< current list of x-pad-coord.
    Int_t                   fIy[100][2];         ///< current list of y-pad-coord.
    Float_t                 fX[100][2];          ///< current list of x-coord.
    Float_t                 fY[100][2];          ///< current list of y-coord.
    Float_t                 fZ[100][2];          ///< current list of z-coord.
    Int_t                   fIndLocal[100][2];   ///< indices of local maxima
    Int_t                   fNLocal[2];          ///< Number of local maxima
    Float_t                   fQ[100][2];          ///< current list of charges
    Float_t                 fZPlane;             ///< currenz z-plane position
    Int_t                   fSector;             ///< current sector
    
    // Current Fit
    Double_t                 fXFit[2];         ///< x-coordinate
    Double_t                 fYFit[2];         ///< y-coordinate
    Double_t                 fQrFit[2];        ///< charge ratio
    Float_t                  fChi2[2];         ///< chi2 of fit
    Float_t                  fXInit[2];        ///< start values
    Float_t                  fYInit[2];        ///< start values
    Float_t                  fQrInit[2];       ///< start values
    Int_t                    fFitStat;         ///< status of fit
    
    // Selected track for debugging
    Int_t                    fTrack[2];        ///< Only digits with main contributions from these tracks are
    // considered 
    
 private:
    AliMUONClusterFinderVS(const AliMUONClusterFinderVS& clusterFinder);
//  Assignment operator
    AliMUONClusterFinderVS & operator = (const AliMUONClusterFinderVS& rhs);

    ClassDef(AliMUONClusterFinderVS,3) //Class for clustering and reconstruction of space points
      };
#endif















