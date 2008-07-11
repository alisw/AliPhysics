#ifndef ALIITSREALIGNTRACKS_H
#define ALIITSREALIGNTRACKS_H

#include <TArray.h>
#include <TFile.h>
#include <TArray.h>
#include "AliGeomManager.h"
#include "AliAlignmentTracks.h"
#include "AliAlignObjParams.h"

class AliITSRealignTracks: public AliAlignmentTracks {
 public:
  
  AliITSRealignTracks():
    AliAlignmentTracks(),
    fSurveyObjs(0),
    fgeomfilename(),
    fmintracks(0),
    fCovIsUsed(kFALSE),
    fUpdateCov(kFALSE){}
  
  AliITSRealignTracks(TString minimizer,Int_t fit=0,Bool_t covUsed=kFALSE,TString fileintro="AliTrackPoints.root",TString geometryfile="geometry.root",TString misalignmentFile="",TString startingfile="");
  AliITSRealignTracks(const AliITSRealignTracks &realignTracks);
  AliITSRealignTracks& operator=(const AliITSRealignTracks& obj);
  ~AliITSRealignTracks();
  
  void InitAlignObjs();
  Bool_t InitSurveyObjs(Bool_t infinite=kFALSE,Double_t factor=1.,Bool_t fromfile=kFALSE,TString filename="",TString arrayName="");
  void ResetAlignObjs();
  void DeleteSurveyObjs();
  Bool_t SelectFitter(Int_t fit,Int_t minTrackPoint=2);
  Bool_t SelectMinimizer(TString minimizer,Int_t minpoints=1,const Bool_t *coord=0x0);
  void SetMinNtracks(Int_t minNtracks){fmintracks=minNtracks;}
  void SetCovUpdate(Bool_t covupdate){fUpdateCov=covupdate;}
  void SetGeomFilename(TString geomfilename){fgeomfilename=geomfilename;}
  //  Int_t LoadPoints(const TArrayI *volids, AliTrackPointArray** &points);
  Bool_t ReadAlignObjs(const char *alignObjFileName = "AlignObjs.root", const char* arrayName = "Alignment");
  void RealignITSVolIndependent(Int_t iter1,Int_t iterations,Int_t minNtracks,Int_t layer=0,Int_t minTrackPoint=6);
  void RealignITStracks(TString minimizer,Int_t fit,Int_t iter1,Int_t iterations,Int_t minNtracks,Int_t layer,Int_t minTrackPoint,Bool_t covUsed,TString misalignmentFile,TString startingfile,Int_t doGlobal);
  Bool_t AlignVolumesITS(const TArrayI *volids, const TArrayI *volidsfit,AliGeomManager::ELayerID layerRangeMin,AliGeomManager::ELayerID layerRangeMax,Int_t iterations);
  Bool_t FirstAlignmentSPD(Int_t minNtracks,Int_t iterations,TArrayI *volidsSet=0x0);
  Bool_t FirstAlignmentLayers(Bool_t *layers,Int_t minNtracks,Int_t iterations,TArrayI *volidsSet=0x0);
  Bool_t SPDmodulesAlignToSSD(Int_t minNtracks,Int_t iterations);
  Bool_t AlignSPDBarrel(Int_t iterations);
  Bool_t AlignSPDHalfBarrel(Int_t method,Int_t iterations);
  Bool_t AlignLayer(Int_t layer,Int_t iterations);
  Bool_t AlignLayersToLayers(Int_t *layer,Int_t iterations);
  Bool_t AlignLayerToSPDHalfBarrel(Int_t layer,Int_t updown,Int_t iterations);
  Bool_t AlignLayerToSector(Int_t layer,Int_t sector,Int_t iterations);
  Bool_t AlignSPDSectorToOuterLayers(Int_t sector,Int_t iterations);
  Bool_t AlignSPDSectorWithSectors(Int_t sector,Int_t iterations);
  Bool_t AlignSPDSectorsWithSectors(Int_t *sectorIN,Int_t *sectorFit,Int_t iterations);
  Bool_t AlignSPDHalfBarrelToHalfBarrel(Int_t updown,Int_t iterations); 
  Bool_t AlignSPDHalfBarrelToSectorRef(Int_t sector,Int_t iterations);
  Bool_t AlignSPD1SectorRef(Int_t sector,Int_t iterations);
  //masera  void AlignGlobalToSectRef(Int_t sector,Int_t minNtracks=100);
  TArrayI* GetLayersVolUID(Int_t *layer);
  AliAlignObjParams* MediateAlignObj(TArrayI *volIDs,Int_t lastVolid);
  TArrayI* GetSPDSectorsVolids(Int_t *sectors); 
  TArrayI* SelectLayerInVolids(const TArrayI *volidsIN,AliGeomManager::ELayerID layer);
  TArrayI* JoinVolArrays(const TArrayI *vol1,const TArrayI *vol2);
  TArrayI* IntersectVolArray(const TArrayI *vol1,const TArrayI *vol2);
  TArrayI* ExcludeVolidsFromVolidsArray(const TArrayI *volidsToExclude,const TArrayI *volStart);
  TArrayI* GetLayerVolumes(Int_t *layer);
  

 private:
  
  AliAlignObj    ***fSurveyObjs;   // Array with survey measurments 
  TString          fgeomfilename; // Geometry file name
  Double_t           fmintracks;   // minimum number of tracks to realign a set of volumes
  Bool_t             fCovIsUsed;   // indicates wheter AlignObj's cov. matrix is used in loading the points 
  Bool_t             fUpdateCov;   // Update of Covariance for AlignObjs

  ClassDef(AliITSRealignTracks,2)
    
    };
    
#endif
