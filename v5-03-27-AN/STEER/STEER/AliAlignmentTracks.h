#ifndef ALIALIGNMENTTRACKS_H
#define ALIALIGNMENTTRACKS_H

//*************************************************************************
// AliAlignmentTracks: main steering class which deals with the alignment *
// procedures based on reconstructed tracks.                              *
// More comments will come with the development of the interfaces and     *
// functionalities of the class.                                          *
//*************************************************************************

#include <TObject.h>

#include "AliAlignObj.h"

class TChain;
class AliTrackPointArray;
class AliAlignObj;
class AliTrackFitter;
class AliTrackResiduals;

class AliAlignmentTracks : public TObject {

 public:

  AliAlignmentTracks();
  AliAlignmentTracks(TChain *esdchain);
  AliAlignmentTracks(const char *esdfilename, const char *esdtreename = "esdTree");
  virtual ~AliAlignmentTracks();

  void AddESD(TChain *esdchain);
  void AddESD(const char *esdfilename, const char *esdtreename = "esdTree");

  void SetPointsFilename(const char *pointsfilename = "AliTrackPoints.root") { fPointsFilename = pointsfilename; }

  void ProcessESD(TSelector *selector);
  void ProcessESD(Bool_t onlyITS=kFALSE,Int_t minITSpts=0,
		  Bool_t cuts=kTRUE,
		  Float_t minAngleWrtITSModulePlanes=0.,
		  Float_t minMom=0.3,Float_t maxMom=1.e9,
		  Float_t minAbsSinPhi=0.,Float_t maxAbsSinPhi=1.,
		  Float_t minSinTheta=0.,Float_t maxSinTheta=1.);
  void ProcessESDCosmics(Bool_t onlyITS=kFALSE,Int_t minITSpts=0,
			 Float_t maxMatchingAngle=0.17, // 10 deg
			 Bool_t cuts=kTRUE,
			 Float_t minAngleWrtITSModulePlanes=0.,
			 Float_t minMom=0.3,Float_t maxMom=1.e9,
			 Float_t minAbsSinPhi=0.,Float_t maxAbsSinPhi=1.,
			 Float_t minSinTheta=0.,Float_t maxSinTheta=1.);

  void BuildIndex();

  Bool_t ReadAlignObjs(const char *alignObjFileName = "AlignObjs.root", const char* arrayName = "Alignment");

  void SetTrackFitter(AliTrackFitter *fitter) { fTrackFitter = fitter; }
  void SetMinimizer(AliTrackResiduals *minimizer) { fMinimizer = minimizer; }

  Bool_t AlignDetector(AliGeomManager::ELayerID firstLayer,
		       AliGeomManager::ELayerID lastLayer,
		       AliGeomManager::ELayerID layerRangeMin = AliGeomManager::kFirstLayer,
		       AliGeomManager::ELayerID layerRangeMax = AliGeomManager::kLastLayer,Int_t iterations = 1);
  Bool_t AlignLayer(AliGeomManager::ELayerID layer,
		    AliGeomManager::ELayerID layerRangeMin = AliGeomManager::kFirstLayer,
		    AliGeomManager::ELayerID layerRangeMax = AliGeomManager::kLastLayer,
		    Int_t iterations = 1);
  Bool_t AlignVolume(UShort_t volId, UShort_t volIdFit,
		     Int_t iterations);
  Bool_t AlignVolumes(const TArrayI *volids, const TArrayI *volidsfit = 0x0,
		      AliGeomManager::ELayerID layerRangeMin = AliGeomManager::kFirstLayer,
		      AliGeomManager::ELayerID layerRangeMax = AliGeomManager::kLastLayer,
		      Int_t iterations = 1);

  AliAlignObj* GetAlignObj(UShort_t volid) const {
    Int_t iModule;
    AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer(volid,iModule);
    return fAlignObjs[iLayer-AliGeomManager::kFirstLayer][iModule];
  }
  void    SetUpdate(Bool_t update){fDoUpdate = update;}
  void SetCovIsUsed(Bool_t covisused){fCovIsUsed=covisused;}
  Bool_t  GetUpdate() const { return fDoUpdate;}
  void WriteRealignObjArray(TString outfilename,AliGeomManager::ELayerID layerRangeMin,AliGeomManager::ELayerID layerRangeMax);
  Int_t GetLastIndex(Int_t iLayer,Int_t iModule) const { return fLastIndex[iLayer][iModule]; }  

  Bool_t Misalign(const char *misalignObjFileName, const char* arrayName);

 protected:

  void InitIndex();
  void ResetIndex();
  void DeleteIndex();

  void InitAlignObjs();
  void ResetAlignObjs();
  void DeleteAlignObjs();

  Int_t LoadPoints(const TArrayI *volids, AliTrackPointArray** &points,Int_t &pointsdim);
  void  UnloadPoints(Int_t n, AliTrackPointArray **points);

  AliTrackFitter *CreateFitter();
  AliTrackResiduals *CreateMinimizer();

  TChain           *fESDChain;       //! Chain with ESDs
  TString           fPointsFilename; //  Name of the file containing the track point arrays
  TFile            *fPointsFile;     //  File containing the track point arrays
  TTree            *fPointsTree;     //  Tree with the track point arrays
  Int_t           **fLastIndex;      //! Last filled index in volume arrays
  TArrayI        ***fArrayIndex;     //! Volume arrays which contains the tree index
  Bool_t            fIsIndexBuilt;   //  Is points tree index built
  AliAlignObj    ***fAlignObjs;      //  Array with alignment objects
  AliAlignObj    ***fMisalignObjs;   //  Array with alignment objects used to introduce misalignment of the space-points
  AliTrackFitter   *fTrackFitter;    //  Pointer to the track fitter
  AliTrackResiduals*fMinimizer;      //  Pointer to track residuals minimizer
  Bool_t            fDoUpdate;       //  Indicator - update Alignment object after minimization
  Bool_t            fCovIsUsed;      //  Indicator - use AlignObjs' Cov matrices

 private:
  AliAlignmentTracks(const AliAlignmentTracks & alignment);
  AliAlignmentTracks& operator= (const AliAlignmentTracks& alignment);

  ClassDef(AliAlignmentTracks,2)

};

#endif
