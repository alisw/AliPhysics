// @(#) $Id$

#ifndef ALILEVEL3_H
#define ALILEVEL3_H

#ifndef no_root
#include <TObject.h>
#include <TFile.h>
#endif

#ifdef use_newio
class AliRunLoader;
#endif

#include "AliHLTDigitData.h"
#include "AliHLTRootTypes.h"

class AliHLTSpacePointData;
class AliHLTDigitRowData;
class AliHLTTrackSegmentData;
class AliHLTDigitData;
class AliHLTConfMapper;
class AliHLTVertex;
class AliHLTVertexFinder;
class AliHLTTrackMerger;
class AliHLTGlobalMerger;
#ifndef no_root
class TDirectory;
#endif
class AliHLTClustFinderNew;
class AliHLTMerger;
class AliHLTInterMerger;
class AliHLTFileHandler;
class AliHLTMemHandler;
class AliHLTBenchmark;

#ifdef no_root
class AliLevel3 {
#else
class AliLevel3 : public TObject {
#endif

 private:
  UInt_t fNTrackData; //count data
  AliHLTTrackSegmentData* fTrackData; //!
  AliHLTConfMapper *fTracker; //!
  AliHLTVertex *fVertex; //! 
  AliHLTVertexFinder *fVertexFinder; //!
  AliHLTTrackMerger *fTrackMerger; //!
  AliHLTGlobalMerger *fGlobalMerger; //!
  AliHLTInterMerger *fInterMerger; //!
  AliHLTClustFinderNew *fClusterFinder; //! 
  AliHLTMemHandler *fFileHandler; //!
  AliHLTBenchmark *fBenchmark;//!

  Int_t fEvent;    //event number
  Int_t fNPatch;   //number of patches
  Int_t fRow[6][2];//rows
  Float_t fEta[2]; //eta
  
  Char_t *fInputFile;//!
#ifdef use_newio
  AliRunLoader *fRunLoader; //runloader
#endif
  Char_t fPath[256]; //path to aliroot
  Char_t fWriteOutPath[256]; //path to store
  
  Bool_t fDoRoi; //do region of interest
  Bool_t fFindVertex; //find vertex
  Bool_t fDoNonVertex;//do non vertex pass
  Bool_t fPileUp; //do pileup
  Bool_t fNoCF; //dont do cluster finder
  
  Bool_t fUseBinary; //use binary input
  Bool_t fWriteOut; //write tracks
  
  static Bool_t fgDoVertexFit; //do vertex fix

  Bool_t fClusterDeconv; //do cluster deconv
  Float_t fXYClusterError; //Cluster error
  Float_t fZClusterError; //Cluster error

  void WriteSpacePoints(UInt_t npoints,AliHLTSpacePointData *points,
                        Int_t slice,Int_t patch) const;
  Int_t WriteTracks(char *filename,AliHLTMerger *merger,char opt='o') const;  
  void WriteResults();
  void FitGlobalTracks();
  void SetPath(char *p){sprintf(fPath,"%s",p);}

 public:
  AliLevel3 ();
  AliLevel3(Char_t *infile);
#ifdef use_newio
  AliLevel3(AliRunLoader *rl);
#endif
  virtual ~AliLevel3();
  enum EFileType {kBinary, kBinary8, kRoot, kRaw, kDate, kRunLoader};
  void Init(Char_t *path,EFileType filetype=kBinary,Int_t npatches=6);
  void SetMergerParameters(Double_t maxy=1.2,Double_t maxz=1.6,Double_t maxkappa=0.003,
			   Double_t maxpsi=0.02,Double_t maxtgl=0.03);
  void SetTrackerParam(Int_t phi_segments=50,Int_t eta_segments=100,
		       Int_t trackletlength=3,Int_t tracklength=5,
		       Int_t rowscopetracklet=2,Int_t rowscopetrack=3,
		       Double_t min_pt_fit=0,Double_t maxangle=1.31,
		       Double_t goodDist=5,Double_t hitChi2Cut=10,
		       Double_t goodHitChi2=20,Double_t trackChi2Cut=50,
		       Int_t maxdist=50,Double_t maxphi=0.1,Double_t maxeta=0.1,
                       Bool_t vertexconstraint=kTRUE);
  void SetClusterFinderParam(Float_t fXYError=0.2,Float_t fZError=0.3,Bool_t deconv=kTRUE);

  void ProcessEvent(Int_t first,Int_t last,Int_t event=0);
  void ProcessSlice(Int_t slice);

  void DoMc(char* file="point_mc.dat");
  void DoNonVertexTracking() {fDoNonVertex=kTRUE;}
  void FindVertex() {fFindVertex=kTRUE;}
  void DoBench(char* name="benchmark");
  void DoPileup() {fPileUp = kTRUE;}
  void NoCF() {fNoCF=kTRUE;}
  void DoRoi(Float_t e0=0.4,Float_t e1=0.5){fEta[0]=e0;fEta[1]=e1;fDoRoi=kTRUE;}
  void WriteFiles(Char_t *path="./"){fWriteOut = kTRUE; sprintf(fWriteOutPath,"%s",path);}
  
  static void SetVertexFit(Bool_t f)   {fgDoVertexFit=f;}
  static Bool_t DoVertexFit()          {return fgDoVertexFit;}

  ClassDef(AliLevel3,1) //Interface class for Level3-tracking
};

#endif





