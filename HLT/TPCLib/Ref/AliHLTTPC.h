// @(#) $Id$

#ifndef ALILEVEL3_H
#define ALILEVEL3_H

#ifndef no_root
#include <TObject.h>
#include <TFile.h>
#endif

#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCRootTypes.h"

class AliHLTTPCSpacePointData;
class AliHLTTPCDigitRowData;
class AliHLTTPCTrackSegmentData;
class AliHLTTPCDigitData;
class AliHLTTPCConfMapper;
class AliHLTTPCVertex;
class AliHLTTPCVertexFinder;
class AliHLTTPCTrackMerger;
class AliHLTTPCGlobalMerger;
#ifndef no_root
class TDirectory;
#endif
class AliHLTTPCClustFinderNew;
class AliHLTTPCMerger;
class AliHLTTPCInterMerger;
class AliHLTTPCFileHandler;
class AliHLTTPCMemHandler;
class AliHLTTPCBenchmark;

#ifdef no_root
class AliHLTTPC {
#else
class AliHLTTPC : public TObject {
#endif

 private:
  UInt_t fNTrackData;
  AliHLTTPCTrackSegmentData* fTrackData; //!
  AliHLTTPCConfMapper *fTracker; //!
  AliHLTTPCVertex *fVertex; //! 
  AliHLTTPCVertexFinder *fVertexFinder; //!
  AliHLTTPCTrackMerger *fTrackMerger; //!
  AliHLTTPCGlobalMerger *fGlobalMerger; //!
  AliHLTTPCInterMerger *fInterMerger; //!
  AliHLTTPCClustFinderNew *fClusterFinder; //! 
  AliHLTTPCMemHandler *fFileHandler; //!
  AliHLTTPCBenchmark *fBenchmark;//!

  Int_t fEvent;
  Int_t fNPatch;
  Int_t fRow[6][2];
  Float_t fEta[2];
  
  Char_t *fInputFile;//!

  Char_t fPath[256];
  Char_t fWriteOutPath[256];
  
  Bool_t fDoRoi;
  Bool_t fFindVertex;
  Bool_t fDoNonVertex;
  Bool_t fPileUp;
  Bool_t fNoCF;
  
  Bool_t fUseBinary;
  Bool_t fWriteOut;
  
  //Define whether track parameters should be given at first point on track (default)
  //If not, the parameters will be given at the vertex.
  static Bool_t fSetTracks2FirstPoint; 

  Bool_t fClusterDeconv;
  Float_t fXYClusterError;
  Float_t fZClusterError;


  void WriteSpacePoints(UInt_t npoints,AliHLTTPCSpacePointData *points,
                                          Int_t slice,Int_t patch);
  Int_t WriteTracks(char *filename,AliHLTTPCMerger *merger,char opt='o');  
  void WriteResults();
  void FitGlobalTracks();
  void SetPath(char *p){sprintf(fPath,"%s",p);}

 public:
  AliHLTTPC ();
  AliHLTTPC(Char_t *infile);
  virtual ~AliHLTTPC();
  enum EFileType {kBinary, kBinary8, kRoot, kRaw, kDate};
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
  
  static void SetTracks2FirstPoint()   {fSetTracks2FirstPoint = kTRUE;}
  static void SetTracks2Vertex()       {fSetTracks2FirstPoint = kFALSE;}
  static Bool_t IsTracksAtFirstPoint() {return fSetTracks2FirstPoint;}
  
  ClassDef(AliHLTTPC,1) //Interface class for Level3-tracking
};

#endif





