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
  AliAlignmentTracks(const AliAlignmentTracks & alignment);
  AliAlignmentTracks& operator= (const AliAlignmentTracks& alignment);
  virtual ~AliAlignmentTracks();

  void AddESD(TChain *esdchain);
  void AddESD(const char *esdfilename, const char *esdtreename = "esdTree");

  void SetPointsFilename(const char *pointsfilename = "AliTrackPoints.root") { fPointsFilename = pointsfilename; }

  void ProcessESD(TSelector *selector);
  void ProcessESD();

  void BuildIndex();

  Bool_t ReadAlignObjs(const char *alignObjFileName = "AlignObjs.root", const char* arrayName = "Alignment");

  void SetTrackFitter(AliTrackFitter *fitter) { fTrackFitter = fitter; }
  void SetMinimizer(AliTrackResiduals *minimizer) { fMinimizer = minimizer; }

  void AlignDetector(AliAlignObj::ELayerID firstLayer,
		     AliAlignObj::ELayerID lastLayer,
		     AliAlignObj::ELayerID layerRangeMin = AliAlignObj::kFirstLayer,
		     AliAlignObj::ELayerID layerRangeMax = AliAlignObj::kLastLayer,Int_t iterations = 1);
  void AlignLayer(AliAlignObj::ELayerID layer,
		  AliAlignObj::ELayerID layerRangeMin = AliAlignObj::kFirstLayer,
		  AliAlignObj::ELayerID layerRangeMax = AliAlignObj::kLastLayer,
		  Int_t iterations = 1);
  void AlignVolume(UShort_t volId, UShort_t volIdFit,
		   Int_t iterations);
  void AlignVolumes(const TArrayI *volids, const TArrayI *volidsfit = 0x0,
		   AliAlignObj::ELayerID layerRangeMin = AliAlignObj::kFirstLayer,
		   AliAlignObj::ELayerID layerRangeMax = AliAlignObj::kLastLayer,
		   Int_t iterations = 1);

  AliAlignObj* GetAlignObj(UShort_t volid) const {
    Int_t iModule;
    AliAlignObj::ELayerID iLayer = AliAlignObj::VolUIDToLayer(volid,iModule);
    return fAlignObjs[iLayer-AliAlignObj::kFirstLayer][iModule];
  }

 protected:

  void InitIndex();
  void ResetIndex();
  void DeleteIndex();

  void InitAlignObjs();
  void ResetAlignObjs();
  void DeleteAlignObjs();

  Int_t LoadPoints(const TArrayI *volids, AliTrackPointArray** &points);
  void  UnloadPoints(Int_t n, AliTrackPointArray **points);

  AliTrackFitter *CreateFitter();
  AliTrackResiduals *CreateMinimizer();

  Bool_t Misalign(const char *misalignObjFileName, const char* arrayName);

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

  ClassDef(AliAlignmentTracks,1)

};

#endif
