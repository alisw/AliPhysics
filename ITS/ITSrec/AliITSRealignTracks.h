#ifndef ALIITSREALIGNTRACKS_H
#define ALIITSREALIGNTRACKS_H


//Class to perform the realignment if the Inner Tracking System 
//with an iterative approach based on track to cluster residuals
// minimization. More details in .cxx file
//
//Class by: A. Rossi, andrea,rossi@ts.infn.it


#include "AliGeomManager.h"
#include "AliAlignmentTracks.h"

class TArray;
class TGraph;
class TCanvas;
class TArray;
class TFile;
class AliAlignObjParams;


/* $Id$ */


class AliITSRealignTracks: public AliAlignmentTracks {
 public:
  
  AliITSRealignTracks():
    AliAlignmentTracks(),
    fSurveyObjs(0),
    fgeomfilename(),
    fmintracks(0),
    fUpdateCov(kFALSE),
    fVarySigmaY(kFALSE),
    fCorrModules(0),
    fLimitCorr(0.),
    fsigmaY(),
    fDraw(kFALSE),  
    fAlignDrawObjs(0), 
    fCanvPar(0), 
    fCanvGr(0), 
    fgrIterMeanX(0), 
    fgrIterRMSX(0),  
    fgrIterMeanY(0), 
    fgrIterRMSY(0),  
    fgrIterMeanZ(0), 
    fgrIterRMSZ(0),  
    fgrIterMeanPsi(0), 
    fgrIterRMSPsi(0),  
    fgrIterMeanTheta(0), 
    fgrIterRMSTheta(0),  
    fgrIterMeanPhi(0), 
    fgrIterRMSPhi(0)  
    {SetCovIsUsed(kFALSE);}
  
  AliITSRealignTracks(TString minimizer,Int_t fit=0,Bool_t covUsed=kFALSE,TString fileintro="AliTrackPoints.root",TString geometryfile="geometry.root",TString misalignmentFile="",TString startingfile="");
  AliITSRealignTracks(const AliITSRealignTracks &realignTracks);
  AliITSRealignTracks& operator=(const AliITSRealignTracks& obj);
  ~AliITSRealignTracks();
  
  void InitAlignObjs();
  Bool_t InitSurveyObjs(Bool_t infinite=kFALSE,Double_t factor=1.,TString filename="",TString arrayName="");
  void ResetAlignObjs(Bool_t all,TArrayI *volids=0x0);
  void ResetCorrModules();
  void DeleteSurveyObjs();
  void SetLimitCorr(Double_t limit=0.1){fLimitCorr=limit;}
  Int_t CheckWithSurvey(Double_t factor=2.,const TArrayI *volids=0x0);
  Bool_t SelectFitter(Int_t fit,Int_t minTrackPoint=2);
  Bool_t SelectMinimizer(TString minimizer,Int_t minpoints=1,const Bool_t *coord=0x0);
  void SetMinNtracks(Int_t minNtracks){fmintracks=minNtracks;}
  void SetCovUpdate(Bool_t covupdate){fUpdateCov=covupdate;}
  void SetVarySigmaY(Bool_t varysigmay,Double_t sigmaYfixed=1.);
  void SetGeomFilename(TString geomfilename){fgeomfilename=geomfilename;}
  //  Int_t LoadPoints(const TArrayI *volids, AliTrackPointArray** &points);
  Bool_t ReadAlignObjs(const char *alignObjFileName = "AlignObjs.root", const char* arrayName = "Alignment");
  void RealignITSVolIndependent(Int_t iter1,Int_t iterations,Int_t minNtracks,Int_t layer=0,Int_t minTrackPoint=6);
  void RealignITStracks(TString minimizer,Int_t fit,Int_t iter1,Int_t iterations,Int_t minNtracks,Int_t layer,Int_t minTrackPoint,Bool_t covUsed,TString misalignmentFile,TString startingfile,Int_t doGlobal);
  Bool_t AlignVolumesITS(const TArrayI *volids, const TArrayI *volidsfit,AliGeomManager::ELayerID layerRangeMin,AliGeomManager::ELayerID layerRangeMax,Int_t iterations);
  Bool_t FirstAlignmentSPD(Int_t minNtracks,Int_t iterations,Bool_t fitall=kTRUE,const TArrayI *volidsSet=0x0);
  Bool_t FirstAlignmentLayers(const Bool_t *layers,Int_t minNtracks,Int_t iterations,Bool_t fitall=kTRUE,const TArrayI *volidsSet=0x0);
  Bool_t SPDmodulesAlignToSSD(Int_t minNtracks,Int_t iterations);
  Bool_t AlignSPDBarrel(Int_t iterations);
  Bool_t AlignSPDHalfBarrel(Int_t method,Int_t iterations);
  Bool_t AlignLayer(Int_t layer,Int_t iterations);
  Bool_t AlignLayersToLayers(const Int_t *layer,Int_t iterations);
  Bool_t AlignLayerToSPDHalfBarrel(Int_t layer,Int_t updown,Int_t iterations);
  Bool_t AlignLayerToSector(Int_t layer,Int_t sector,Int_t iterations);
  Bool_t AlignSPDSectorToOuterLayers(Int_t sector,Int_t iterations);
  Bool_t AlignSPDSectorWithSectors(Int_t sector,Int_t iterations);
  Bool_t AlignSPDSectorsWithSectors(const Int_t *sectorIN,const Int_t *sectorFit,Int_t iterations);
  Bool_t AlignSPDStaves(const Int_t *staves,const Int_t *sectorsIN,const Int_t *sectorsFit,Int_t iterations);
  Bool_t AlignSPDHalfBarrelToHalfBarrel(Int_t updown,Int_t iterations); 
  Bool_t AlignSPDHalfBarrelToSectorRef(Int_t sector,Int_t iterations);
  Bool_t AlignSPD1SectorRef(Int_t sector,Int_t iterations);
  //masera  void AlignGlobalToSectRef(Int_t sector,Int_t minNtracks=100);
  TArrayI* GetLayersVolUID(const Int_t *layer);
  AliAlignObjParams* MediateAlignObj(const TArrayI *volIDs,Int_t lastVolid);
  TArrayI* GetSPDSectorsVolids(const Int_t *sectors); 
  TArrayI* GetSPDStavesVolids(const Int_t *sectors,const Int_t* staves);
  TArrayI* SelectLayerInVolids(const TArrayI *volidsIN,AliGeomManager::ELayerID layer);
  TArrayI* JoinVolArrays(const TArrayI *vol1,const TArrayI *vol2);
  TArrayI* IntersectVolArray(const TArrayI *vol1,const TArrayI *vol2);
  TArrayI* ExcludeVolidsFromVolidsArray(const TArrayI *volidsToExclude,const TArrayI *volStart);
  TArrayI* GetLayerVolumes(const Int_t *layer);
  TArrayI* GetAlignedVolumes(char *filename);
  /*  void AlignGlobalToSectRef(Int_t sector,Int_t minNtracks=100);
      AliAlignObjParams* MediateAlignObjs(AliAlignObj **alObjs,Int_t nObjs,const Bool_t *coords=0x0,TArrayI *volidArray=0x0,Bool_t local=kFALSE,const char* geometryfile=0x0);
      Bool_t MediateSectorsVolumes(char *filename,Bool_t local=kFALSE,char *geometryfile=0x0,Bool_t *coord=0x0);
  */
  void InitDrawHists();
  void SetDraw(Bool_t draw,Bool_t refresh);
  void UpdateDraw(TArrayI *volids,Int_t iter,Int_t color);
  void DeleteDrawHists();
  void WriteHists(const char *outfile);

 private:
  
  AliAlignObj    ***fSurveyObjs;   // Array with survey measurments 
  TString          fgeomfilename; // Geometry file name
  Int_t           fmintracks;   // minimum number of tracks to realign a set of volumes
  Bool_t             fUpdateCov;   // Update of Covariance for AlignObjs
  Bool_t            fVarySigmaY;   // If kTRUE the "sigmaY" parameter is changed accordingly to alignObj error
  Double_t          **fCorrModules;  //  Used to reduce correlations between modules
  Double_t           fLimitCorr;    // Maximum number of tracks shared between modules
  Double_t            fsigmaY;    // sigmaY parameter
  Bool_t              fDraw;     // flag to activate draw objects
  AliAlignObj  ***fAlignDrawObjs; //Array with reference objects for histograms
  TCanvas              *fCanvPar; //Canvas with iterations distributions
  TCanvas               *fCanvGr; //Canvas with iterations results
  TGraph           *fgrIterMeanX; // graph of Delta X
  TGraph           *fgrIterRMSX;  // graph of RMS X
  TGraph           *fgrIterMeanY; // graph of Delta Y
  TGraph           *fgrIterRMSY;  // graph of RMS Y
  TGraph           *fgrIterMeanZ; // graph of DeltaZ
  TGraph           *fgrIterRMSZ;  // TGraphs for displaying results during iterations
  TGraph           *fgrIterMeanPsi; // graphs during iterations
  TGraph           *fgrIterRMSPsi;  // graphs during iterations
  TGraph           *fgrIterMeanTheta; // graphs during iterations
  TGraph           *fgrIterRMSTheta;  // graphs during iterations
  TGraph           *fgrIterMeanPhi; // graphs during iterations
  TGraph           *fgrIterRMSPhi;  // graphs during iterations

  ClassDef(AliITSRealignTracks,3)
    
    };
    
#endif
