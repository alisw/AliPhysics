#ifndef ALILEVEL3_H
#define ALILEVEL3_H

#include <TObject.h>
#include <TBuffer.h>
#include <TFile.h>

#include "AliL3DigitData.h"

class AliL3SpacePointData;
class AliL3DigitRowData;
class AliL3TrackSegmentData;
class AliL3DigitData;
class AliL3Transform;
class TClonesArray;
class AliTPCParam;
class AliL3ConfMapper;
class AliL3Vertex;
class AliL3VertexFinder;
class AliL3TrackMerger;
class AliL3GlobalMerger;
class TDirectory;
class AliL3Transform;
class AliL3ClustFinder;
class AliL3Merger;
class AliL3InterMerger;
class AliL3FileHandler;
class AliL3Benchmark;

class AliLevel3 : public TObject {

 private:
  UInt_t fNTrackData;
  AliL3TrackSegmentData* fTrackData; //!
  AliL3ConfMapper *fTracker; //!
  AliL3Vertex *fVertex; //! 
  AliL3VertexFinder *fVertexFinder; //!
  AliL3TrackMerger *fTrackMerger; //!
  AliL3GlobalMerger *fGlobalMerger; //!
  AliL3InterMerger *fInterMerger; //!
  AliL3ClustFinder *fClusterFinder; //! 
  AliL3FileHandler *fFileHandler; //!
  AliL3Benchmark *fBenchmark;//!
  Int_t fNPatch;
  Char_t fPath[256];
  AliL3Transform *fTransformer; //!
  TDirectory *savedir;
  TFile *fInputFile;
  TFile *fOutputFile;
  Bool_t fFindVertex;
  void Init();
  void WriteSpacePoints(UInt_t npoints,AliL3SpacePointData *points,
                                             Int_t slice,Int_t patch);
  void WriteResults();
  Int_t WriteTracks(char *filename,AliL3Merger *merger,char opt='o');  
  Float_t fEta[2];
  Bool_t fDoRoi;
  Bool_t fUseBinary;
  Bool_t fWriteOut;
  void SetPath(char *p){sprintf(fPath,"%s",p);}
 public:
  AliLevel3 ();
  AliLevel3(Char_t *infile,Char_t *outfile);
  AliLevel3(TFile *in, TFile *out);
  virtual ~AliLevel3();
  
  void SetTrackerParam(Int_t phi_segments=50,Int_t eta_segments=100,
		       Int_t trackletlength=3,Int_t tracklength=5,
		       Int_t rowscopetracklet=2,Int_t rowscopetrack=3,
		       Double_t min_pt_fit=0,Double_t maxangle=1.31,
		       Double_t goodDist=5,Double_t hitChi2Cut=10,
		       Double_t goodHitChi2=20,Double_t trackChi2Cut=50,
		       Int_t maxdist=50);


  void ProcessEvent(Int_t first,Int_t last);
  void ProcessSlice(Int_t slice);


  void UseBinaryInput(char *path){SetPath(path);fUseBinary=kTRUE;}
  void DoMc(char* file="point_mc.dat");
  void DoBench(char* name="benchmark");
  void DoRoi(Float_t e0=0.4,Float_t e1=0.5){fEta[0]=e0;fEta[1]=e1;fDoRoi=kTRUE;}
  void WriteFiles(){fWriteOut = kTRUE;}
  ClassDef(AliLevel3,1) //Interface class for Level3-tracking
};

#endif
