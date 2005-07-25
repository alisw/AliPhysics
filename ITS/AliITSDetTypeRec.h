#ifndef ALIITSDETTYPEREC_H
#define ALIITSDETTYPEREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$ 
*/

////////////////////////////////////////////////////////////////////////
// This class contains all of the "external" information needed to do //
// detector specific reconstruction for the ITS.                      //
////////////////////////////////////////////////////////////////////////
#include <TObject.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliITSLoader.h"
#include "AliRunLoader.h"

class AliITSsegmentation;
class AliITSresponse;
class AliITSClusterFinder;
class AliITSRawCluster;
class AliITSRecPoint;
class AliITSclusterV2;
class AliRawReader;
class AliITSgeom;

class AliITSDetTypeRec : public TObject {
  public:
    AliITSDetTypeRec(); // Default constructor
    AliITSDetTypeRec(const AliITSDetTypeRec& rec);
    AliITSDetTypeRec& operator=(const AliITSDetTypeRec &source);

    virtual ~AliITSDetTypeRec(); // Proper Destructor
    AliITSgeom* GetITSgeom()const{return fGeom;}
    void SetITSgeom(AliITSgeom *geom) {fGeom=geom;}

    virtual void SetLoader(AliITSLoader* loader) {fLoader=loader;}
    AliITSLoader* GetLoader() const {return fLoader;}
    virtual void SetDefaults();
    virtual void SetDefaultClusterFinders();
    virtual void SetDefaultClusterFindersV2(Bool_t rawdata=kFALSE);
    virtual void MakeBranch(Option_t *opt);
    virtual void SetTreeAddress();
    virtual void SetTreeAddressD(TTree* treeD);

    virtual void SetSegmentationModel(Int_t dettype, AliITSsegmentation *seg);
    virtual void SetCalibrationModel(Int_t dettype, AliITSresponse *cal);
    virtual void SetReconstructionModel(Int_t dettype, AliITSClusterFinder *rec);
    virtual AliITSsegmentation* GetSegmentationModel(Int_t dettype);
    virtual AliITSresponse* GetCalibrationModel(Int_t dettype);
    virtual AliITSClusterFinder* GetReconstructionModel(Int_t dettype);

    virtual void SetDigitClassName(Int_t i,Char_t *digit) 
      {fDigClassName[i]=digit;}
    virtual void SetClusterClassName(Int_t i,Char_t *cluster)
      {fClusterClassName[i]=cluster;}
    virtual void SetRecPointClassName(Int_t i,Char_t *recpoint) 
      {fRecPointClassName[i]=recpoint;}
    
    Char_t* GetDigitClassName(Int_t i) const {return fDigClassName[i];}
    Char_t* GetClusterClassName(Int_t i) const {return fClusterClassName[i];}
    Char_t* GetRecPointClassName(Int_t i) const {return fRecPointClassName[i];}
    
    TObjArray* GetDigits() const {return fDigits;} 
    Int_t *Ndtype() {return fNdtype;}
    TClonesArray *DigitsAddress(Int_t id) const {return ((TClonesArray*)(*fDigits)[id]);}
    virtual void SelectVertexer(TString sel=" "){fSelectedVertexer = sel;}
    //
    virtual void MakeTreeC();
    virtual void GetTreeC(Int_t event);
    virtual void AddCluster(Int_t branch, AliITSRawCluster *c);
    virtual void ResetClusters(); 
    virtual void ResetClusters(Int_t branch);
    virtual void MakeBranchC();
    TBranch* MakeBranchInTree(TTree *tree, const char* name, const char *classname, void* address,Int_t size, Int_t splitlevel, const char */*file*/);

    TObjArray    *Ctype()  {return fCtype;}
    Int_t        *Nctype() {return fNctype;}
    TClonesArray *ClustersAddress(Int_t id) const {return ((TClonesArray*)(*fCtype)[id]);}
    virtual void ResetDigits();
    virtual void ResetDigits(Int_t branch);

    void MakeBranchR(const char *file, Option_t *opt=" ");
    void SetTreeAddressR(TTree *treeR);
    void AddRecPoint(const AliITSRecPoint &p);
    void ResetRecPoints(){if(fRecPoints) fRecPoints->Clear();fNRecPoints = 0;};
    // Return pointer to rec points 
    TClonesArray  *RecPoints()   {return fRecPoints;}
    void AddClusterV2(const AliITSclusterV2 &cl);
    void ResetClustersV2(){if(fClustersV2) fClustersV2->Clear();fNClustersV2=0;} 
    Int_t GetNClustersV2()const {return fNClustersV2;}

   TClonesArray *ClustersV2() {return fClustersV2;}

    void MakeBranchRF(const char *file){MakeBranchR(file,"Fast");}
    //    void HitsToFastRecPoints(Int_t evNumber,Int_t bgrev,Int_t size,
    //             Option_t *add, Option_t *det, const char *filename);
    void Digits2Reco(){
        DigitsToRecPoints(fLoader->GetRunLoader()->GetEventNumber(),0,"All");}
    void DigitsToRecPoints(Int_t evNumber,Int_t lastEntry,Option_t *det,Bool_t v2=kFALSE);
    void DigitsToRecPoints(AliRawReader* rawReader);

  private:

    static const Int_t fgkNdettypes; // number of det. types
    AliITSgeom   *fGeom;          // ITS geometry
    TObjArray    *fReconstruction;//! [NDet]
    TObjArray    *fSegmentation;  //! [NDet]
    TObjArray    *fCalibration;   //! [NMod]
    TObjArray    *fPreProcess;    //! [] e.g. Find Calibration values
    TObjArray    *fPostProcess;   //! [] e.g. find primary vertex
    TObjArray    *fDigits;        //! [NMod][NDigits]
    Int_t        *fNdtype;        //! detector types  
    Char_t*       fClusterClassName[3]; //! String with Cluster class name
    Char_t*       fDigClassName[3];     //! String with digit class name.
    Char_t*       fRecPointClassName[3];//! String with RecPoint class name

    TObjArray    *fCtype;      //! List of clusters
    Int_t        *fNctype;     //[fNDetTypes] Num. of clust. per type of det.

    TClonesArray *fRecPoints;  //! List of reconstructed points
    Int_t         fNRecPoints; // Number of rec points

    TClonesArray *fClustersV2; //!List of reconstructed clusters v2
    Int_t         fNClustersV2;    //Number of clusters v2
    TString fSelectedVertexer; // Vertexer selected in CreateVertexer
    AliITSLoader* fLoader;     // ITS loader

    ClassDef(AliITSDetTypeRec,1) // ITS Reconstruction structure
};

#endif
