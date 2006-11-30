// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>, Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#ifndef no_root
#include <TFile.h>
#include <TDirectory.h>
#include <TClonesArray.h>
#include <TStopwatch.h>
#include <Riostream.h>
#endif

#ifdef use_newio
#include <AliRunLoader.h>
#endif

#include "AliHLTLogging.h"
#include "AliLevel3.h"
#include "AliHLTConfMapper.h"
#include "AliHLTVertex.h"
#include "AliHLTVertexFinder.h"
#include "AliHLTTrackMerger.h"
#include "AliHLTGlobalMerger.h"
#include "AliHLTInterMerger.h"
#include "AliHLTConfMapPoint.h"
#include "AliHLTConfMapTrack.h"
#include "AliHLTTransform.h"
#include "AliHLTClustFinderNew.h"
#include "AliHLTDigitData.h"
#include "AliHLTTrackArray.h"
#include "AliHLTMemHandler.h"
#include "AliHLTFitter.h"
#ifdef use_aliroot
#include "AliHLTFileHandler.h"
#endif
#include "AliHLTBenchmark.h"
#include "AliHLTDigitData.h"
#include "AliHLTTrackSegmentData.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTVertexData.h"
#include "AliHLTDDLDataFileHandler.h"

/** \class AliLevel3
<pre>
//_____________________________________________________________
//
//  AliLevel3
//
//  Interface class for Level3 tracker.
//  For example how to use, see exa/runtracker.C (root)
//  or programs/runtracker.cxx (standalone program).
//Begin_Html 
//<img src="tpcsectorsnb.gif">
//End_Html
</pre>
*/

ClassImp(AliLevel3)

Bool_t AliLevel3::fgDoVertexFit = kTRUE;//Include the vertex in the final track fit

AliLevel3::AliLevel3()
{
  //Default constructor. Should also be used when input is from binary files.
  //In that case the path to where the binary files are located has to be
  //passed to the AliLevel::Init function.
  
  fVertexFinder=0;
  fVertex=0;
  fTracker=0;
  fTrackMerger=0;
  fInterMerger=0;
  fFileHandler=0;
  fGlobalMerger=0;
  fInputFile=0;
#ifdef use_newio
  fRunLoader=0;
#endif
}

AliLevel3::AliLevel3(Char_t *infile)
{
  //Constructor to use for when input is anything else but binary files,
  //meaning rootfiles or raw files.
  
  fVertexFinder=0;
  fVertex=0;
  fTracker=0;
  fTrackMerger=0;
  fInterMerger=0;
  fFileHandler=0;
  fGlobalMerger=0;
  fInputFile = infile;
#ifdef use_newio
  fRunLoader=0;
#endif
}

#ifdef use_newio
AliLevel3::AliLevel3(AliRunLoader *rl)
{
  //Constructor to use when input is aliroot runloader
  fVertexFinder=0;
  fVertex=0;
  fTracker=0;
  fTrackMerger=0;
  fInterMerger=0;
  fFileHandler=0;
  fGlobalMerger=0;
  fInputFile=0;
  fRunLoader = rl;
}
#endif

void AliLevel3::Init(Char_t *path,EFileType filetype,Int_t npatches)
{
  //Init the whole standard tracker chain
#ifndef use_newio
  if (filetype==kRunLoader){
    LOG(AliHLTLog::kError,"AliLevel3::Init","Files")
	<<"You have not supplied the input rootfile; if you want "
	<<"to run with RunLoader use -Duse_newio for compiling!"<<ENDLOG;
  }
#endif

  if((filetype!=kBinary) && (filetype!=kDate) 
       && (filetype!=kRunLoader)&& !fInputFile)
    {
      LOG(AliHLTLog::kError,"AliLevel3::Init","Files")
	<<"You have not supplied the input rootfile; use the appropriate ctor!"<<ENDLOG;
      return;
    }
#if use_newio
  if((filetype==kRunLoader) && !fRunLoader)
    {
      LOG(AliHLTLog::kError,"AliLevel3::Init","Files")
	<<"You have not supplied the input runloader; use the appropriate ctor!"<<ENDLOG;
      return;
    }
#endif
  
  fWriteOut = kFALSE;
  fPileUp = kFALSE;
  fNoCF=kFALSE;
  fUseBinary = (filetype==kBinary);
  SetPath(path);
  
  fDoRoi = kFALSE;
  fDoNonVertex = kFALSE;
  fFindVertex = kFALSE;
  SetClusterFinderParam();

  fEta[0] = 0.;
  fEta[1] = 1.1;

  fEvent=0;
#ifdef use_aliroot /*just to be sure*/
  AliHLTFileHandler::CleanStaticIndex();
#endif

  switch(npatches){
  case 0:
    fNPatch = 1;
    fRow[0][0] = AliHLTTransform::GetFirstRow(3);
    fRow[0][1] = AliHLTTransform::GetLastRow(5);
    break;
  case 1:
    fNPatch = 1;        //number of patches change row in process
    fRow[0][0] = 0;
    fRow[0][1] = AliHLTTransform::GetLastRow(-1);
    break;
  case 2:
    fNPatch = 2;        //number of patches change row in process
    fRow[0][0] = 0;     // first row
    fRow[0][1] = AliHLTTransform::GetLastRow(1);
    fRow[1][0] = AliHLTTransform::GetFirstRow(2);
    fRow[1][1] = AliHLTTransform::GetLastRow(5);
    break;
  default: 
    fNPatch = 6;        
    fRow[0][0] = AliHLTTransform::GetFirstRow(0);
    fRow[0][1] = AliHLTTransform::GetLastRow(0);
    fRow[1][0] = AliHLTTransform::GetFirstRow(1);
    fRow[1][1] = AliHLTTransform::GetLastRow(1);
    fRow[2][0] = AliHLTTransform::GetFirstRow(2);
    fRow[2][1] = AliHLTTransform::GetLastRow(2);
    fRow[3][0] = AliHLTTransform::GetFirstRow(3);
    fRow[3][1] = AliHLTTransform::GetLastRow(3);
    fRow[4][0] = AliHLTTransform::GetFirstRow(4);
    fRow[4][1] = AliHLTTransform::GetLastRow(4);
    fRow[5][0] = AliHLTTransform::GetFirstRow(5);
    fRow[5][1] = AliHLTTransform::GetLastRow(5);
  }

  fVertexFinder = new AliHLTVertexFinder();
  fVertex = new AliHLTVertex();
  fTracker = new AliHLTConfMapper();
  fTrackMerger = new AliHLTTrackMerger(fNPatch);
  fInterMerger = new AliHLTInterMerger();
  fGlobalMerger = new AliHLTGlobalMerger();
  SetMergerParameters();//Set default merger parameters
#ifdef use_aliroot
  if(filetype==kRoot){
    fFileHandler = new AliHLTFileHandler(kTRUE); //static version
    fFileHandler->SetAliInput(fInputFile);
  }else if(filetype==kRaw){
    fFileHandler = new AliHLTDDLDataFileHandler();
    fFileHandler->SetReaderInput(fInputFile);
  }else if(filetype==kDate){
    fFileHandler = new AliHLTDDLDataFileHandler();
    fFileHandler->SetReaderInput(fInputFile,-1);
  }
#if use_newio
  else if(filetype==kRunLoader){
    fFileHandler = new AliHLTFileHandler(kTRUE); //static version
    fFileHandler->SetAliInput(fRunLoader);
  }
#endif
  else{
    fFileHandler = new AliHLTMemHandler();
  }
#else
  if(filetype==kRaw){
    fFileHandler = new AliHLTDDLDataFileHandler();
    fFileHandler->SetReaderInput(fInputFile);
  }else{
    fFileHandler = new AliHLTMemHandler();
  }
#endif
  fBenchmark = new AliHLTBenchmark();
}

void AliLevel3::DoBench(char* name)
{ 
  //dobench
  fBenchmark->Analyze(name);
  delete fBenchmark;
  fBenchmark = new AliHLTBenchmark();
}

void AliLevel3::DoMc(char* file)
{ 
  //domc
#ifdef use_aliroot
  if(!fFileHandler->IsDigit(fEvent))
    fFileHandler->SetMCOutput(file);
#endif
}

AliLevel3::~AliLevel3()
{
  //Destructor
  if(fVertexFinder) delete fVertexFinder;
  if(fVertex) delete fVertex;
  if(fTracker) delete fTracker;
  if(fTrackMerger) delete fTrackMerger;
  if(fInterMerger) delete fInterMerger;
  if(fFileHandler) delete fFileHandler;
  if(fGlobalMerger) delete fGlobalMerger;
}

void AliLevel3::SetClusterFinderParam(Float_t fXYError, Float_t fZError, Bool_t deconv)
{ 
  //set cluster finder parameter
  fXYClusterError=fXYError;
  fZClusterError=fZError;
  fClusterDeconv=deconv;
}

void AliLevel3::SetTrackerParam(Int_t phi_segments, Int_t eta_segments,
				Int_t trackletlength, Int_t tracklength,
				Int_t rowscopetracklet, Int_t rowscopetrack,
				Double_t min_pt_fit, Double_t maxangle,
				Double_t goodDist, Double_t hitChi2Cut,
				Double_t goodHitChi2, Double_t trackChi2Cut,
				Int_t maxdist,Double_t maxphi,Double_t maxeta,
                               Bool_t vertexconstraint)
{
  //Set parameters input to the tracker
  //If no arguments are given, default parameters will be used
  
  fTracker->SetNSegments(phi_segments,eta_segments);
  fTracker->SetMaxDca(min_pt_fit);
  fTracker->SetTrackCuts(hitChi2Cut,goodHitChi2,trackChi2Cut,maxdist,vertexconstraint);
  fTracker->SetTrackletCuts(maxangle,goodDist,vertexconstraint);
  if(vertexconstraint)
    fTracker->MainVertexSettings(trackletlength,tracklength,rowscopetracklet,rowscopetrack,maxphi,maxeta);
  else
    fTracker->NonVertexSettings(trackletlength,tracklength,rowscopetracklet,rowscopetrack);
  fTracker->InitVolumes();
}

void AliLevel3::SetMergerParameters(Double_t maxy,Double_t maxz,Double_t maxkappa,Double_t maxpsi,Double_t maxtgl)
{
  //set global merger parameter
  fGlobalMerger->SetParameter(maxy,maxz,maxkappa,maxpsi,maxtgl);
}

void AliLevel3::ProcessEvent(Int_t first,Int_t last,Int_t event)
{
  //Do tracking on all slices in region [first,last]
  //Slices numbering in TPC goes from 0-35, which means that one slice
  //corresponds to inner+outer sector.E.g. slice 2 corresponds to
  //inner=2 + outer=38.

  fGlobalMerger->Setup(first,last);
#ifdef use_aliroot
  if(fEvent!=event) AliHLTFileHandler::CleanStaticIndex();
#endif
  fEvent=event;
  for(Int_t i=first; i<=last; i++){
    ProcessSlice(i);
    fGlobalMerger->SetVertex(fVertex);
    fGlobalMerger->InitSlice(i);
    fGlobalMerger->FillTracks(fNTrackData,fTrackData);
    fFileHandler->Free();   //free the memory
    fNTrackData=0;
    fTrackData=0;
  }
  fBenchmark->Start("Global track merger");
  //fGlobalMerger->AddAllTracks();
  fGlobalMerger->Merge();
  //fGlobalMerger->SlowMerge(fWriteOutPath);
  fBenchmark->Stop("Global track merger");
  
  FitGlobalTracks();
  
  if(fWriteOut) WriteResults(); 
  fFileHandler->FreeDigitsTree();
}

void AliLevel3::ProcessSlice(Int_t slice)
{
  //process slice
  char name[256];
  Bool_t UseCF = kFALSE;
#ifdef use_aliroot
  UseCF = fFileHandler->IsDigit(fEvent);
#endif
  if(fUseBinary)
    UseCF = kTRUE;   //In case you are not using aliroot
  if(fNoCF == kTRUE) //In case you don't want to run with cluster finder
    UseCF = kFALSE;

  const Int_t kmaxpoints=120000;
  const Int_t kpointsize = kmaxpoints * sizeof(AliHLTSpacePointData);
  AliHLTMemHandler *memory = new AliHLTMemHandler();

  fTrackMerger->Reset();
  fTrackMerger->SetRows(fRow[0]);
  
  for(Int_t patch=fNPatch-1;patch>=0;patch--){
    fFileHandler->Init(slice,patch,&fRow[patch][0]);
    UInt_t npoints=0;
    AliHLTSpacePointData *points =0;
    UInt_t ndigits=0;
    AliHLTDigitRowData *digits =0;
    if(UseCF){
      if(fUseBinary){
        if(!fDoRoi){ 
          if(1){     //Binary to Memory
	    fFileHandler->Free();
            if(fNPatch == 1)
	      sprintf(name,"%sdigits_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	    else
	      sprintf(name,"%sdigits_%d_%d_%d.raw",fPath,fEvent,slice,patch);
	    if(!fFileHandler->SetBinaryInput(name)) return;
	    if(fPileUp)
	      { //Read binary files which are not RLE
		digits = (AliHLTDigitRowData*)fFileHandler->Allocate();
		fFileHandler->Binary2Memory(ndigits,digits); 
	      }
	    else //Read RLE binary files
	      digits= (AliHLTDigitRowData *)fFileHandler->CompBinary2Memory(ndigits);

	    fFileHandler->CloseBinaryInput(); 
          }

          if(0){     //Binary to Memory with Benchmark 
	    fFileHandler->Free();
            if(fNPatch == 1)
	      sprintf(name,"%sdigits_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	    else
	      sprintf(name,"%sdigits_%d_%d_%d.raw",fPath,fEvent,slice,patch);
            if(!memory->SetBinaryInput(name)) return;
            UInt_t compsize=memory->GetFileSize();
            UInt_t *comp=(UInt_t *)memory->Allocate(compsize);
            memory->CompBinary2CompMemory(ndigits,comp);
            memory->CloseBinaryInput();
            UInt_t datasize=memory->GetMemorySize(ndigits,comp);
            digits=(AliHLTDigitRowData *)fFileHandler->Allocate(datasize);
            fBenchmark->Start("Unpacker"); 
            fFileHandler->CompMemory2Memory(ndigits,digits,comp); 
            fBenchmark->Stop("Unpacker");
            memory->Free();
          }
  
          if(0){     //Binary to Memory with Random
            fFileHandler->Free();
            fFileHandler->ResetRandom();
            fFileHandler->SetRandomCluster(100);
            fFileHandler->SetNGenerate(100);
            if(fNPatch == 1)
	      sprintf(name,"%sdigits_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	    else
	      sprintf(name,"%sdigits_%d_%d_%d.raw",fPath,fEvent,slice,patch);
            if(!memory->SetBinaryInput(name)) return;
            UInt_t compsize=memory->GetFileSize();
            UInt_t *comp=(UInt_t *)memory->Allocate(compsize);
            memory->CompBinary2CompMemory(ndigits,comp);
            memory->CloseBinaryInput();
            UInt_t dsize=memory->GetMemorySize(ndigits,comp);
            UInt_t rsize=fFileHandler->GetRandomSize();       
            digits=(AliHLTDigitRowData*)fFileHandler->Allocate(dsize+rsize);
            fBenchmark->Start("Unpacker");
            fFileHandler->CompMemory2Memory(ndigits,digits,comp); 
            fBenchmark->Stop("Unpacker");
            memory->Free();
          }
        }

        else{     //Binary to Memory with Roi
          fFileHandler->Free();
          Int_t sli[2]={0,0};
          fFileHandler->SetROI(fEta,sli);
          if(fNPatch==1)
	    sprintf(name,"%sdigits_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	  else
	    sprintf(name,"%sdigits_%d_%d_%d.raw",fPath,fEvent,slice,patch);
          if(!memory->SetBinaryInput(name)) return;
          UInt_t compsize=memory->GetFileSize();
          UInt_t *comp=(UInt_t *)memory->Allocate(compsize);
          memory->CompBinary2CompMemory(ndigits,comp);
          memory->CloseBinaryInput();
          UInt_t datasize=memory->GetMemorySize(ndigits,comp);
          digits=(AliHLTDigitRowData *)fFileHandler->Allocate(datasize);
          fBenchmark->Start("Unpacker"); 
          datasize = fFileHandler->CompMemory2Memory(ndigits,digits,comp); 
          fBenchmark->Stop("Unpacker"); 
          memory->Free();
        }
      }//end UseBinary
      else{
#ifdef use_aliroot
        fBenchmark->Start("Dummy Unpacker");
        if(fNPatch==1)
	  sprintf(name,"digits_%d_%d_%d.raw",fEvent,slice,-1);
	else
	  sprintf(name,"digits_%d_%d_%d.raw",fEvent,slice,patch);
        fBenchmark->Stop("Dummy Unpacker");

        if(0){   //Ali to Binary
          fFileHandler->SetBinaryOutput(name);
          fFileHandler->AliDigits2CompBinary();
          fFileHandler->CloseBinaryOutput();
        }
  
        if(1){   //Ali to Memory
          digits=(AliHLTDigitRowData *)fFileHandler->AliAltroDigits2Memory(ndigits,fEvent);
          if(0){ //Memory to Binary
            fFileHandler->SetBinaryOutput(name);
            fFileHandler->Memory2CompBinary(ndigits,digits);
            fFileHandler->CloseBinaryOutput();
          }
        }
#endif
      }//end else UseBinary

      points = (AliHLTSpacePointData *) memory->Allocate(kpointsize);
      fClusterFinder = new AliHLTClustFinderNew();
      fClusterFinder->InitSlice(slice,patch,fRow[patch][0],fRow[patch][1],kmaxpoints);
      fClusterFinder->SetDeconv(fClusterDeconv);
      fClusterFinder->SetXYError(fXYClusterError);
      fClusterFinder->SetZError(fZClusterError);
      if((fXYClusterError>0)&&(fZClusterError>0))
	fClusterFinder->SetCalcErr(kFALSE);
      fClusterFinder->SetOutputArray(points);
      fBenchmark->Start("Cluster finder");
      fClusterFinder->Read(ndigits,digits);
      fClusterFinder->ProcessDigits();
      fBenchmark->Stop("Cluster finder");
      npoints = fClusterFinder->GetNumberOfClusters();
      delete fClusterFinder;
      fClusterFinder = 0;
      fFileHandler->Free();
      LOG(AliHLTLog::kInformational,"AliLevel3::ProcessSlice","Cluster Finder")
        <<AliHLTLog::kDec<<"Found "<<npoints<<" Points"<<ENDLOG;
    }//end UseCF
    else{// if not use Clusterfinder
      if(fUseBinary){//Binary to Memory
        memory->Free();
        if(fNPatch==1)
	  sprintf(name,"%s/points_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	else
	  sprintf(name,"%s/points_%d_%d_%d.raw",fPath,fEvent,slice,patch);
        if(!memory->SetBinaryInput(name)) return;
        points = (AliHLTSpacePointData *) memory->Allocate();
        memory->Binary2Memory(npoints,points);
        memory->CloseBinaryInput();
        LOG(AliHLTLog::kInformational,"AliLevel3::ProcessSlice","Read Cluster")
        <<AliHLTLog::kDec<<"Found "<<npoints<<" Points in File"<<ENDLOG;
      }
#ifdef use_aliroot
      else{
        points = fFileHandler->AliPoints2Memory(npoints);
      }
#endif
      fBenchmark->Start("Dummy Unpacker");
      fBenchmark->Stop("Dummy Unpacker");
      fBenchmark->Start("Dummy CF");
      fBenchmark->Stop("Dummy CF");
    }
   
    if(patch == fNPatch-1){
      // Vertex
      if(fFindVertex){
      // Vertex Finder
      
        fBenchmark->Start("Vertex Finder Read");
        fVertexFinder->Reset();
        fVertexFinder->Read(npoints,points);
        fBenchmark->Stop("Vertex Finder Read"); 
        fBenchmark->Start("Vertex Finder");
        fVertexFinder->Analyze();
        AliHLTVertexData vertex[1];
        fVertexFinder->Write(vertex);
        fVertex->Read(vertex);
        fBenchmark->Stop("Vertex Finder"); 
      }
      else{
        //use 0,0,0 for vertex
        fVertex->SetZero();
      }
      fTrackMerger->SetVertex(fVertex);
    }

    fTracker->InitSector(slice,fRow[patch],fEta);
    fTracker->SetVertex(fVertex);
    fBenchmark->Start("Tracker setup"); 
    fTracker->ReadHits(npoints,points);
    fTracker->MainVertexTracking_a();
    fBenchmark->Stop("Tracker setup");
    fBenchmark->Start("Track follower");
    fTracker->MainVertexTracking_b();
    fBenchmark->Stop("Track follower");
    if(fDoNonVertex)
      fTracker->NonVertexTracking();//Do a second pass for nonvertex tracks
    fBenchmark->Start("Sector track fitting");
    fTracker->FillTracks();
    fBenchmark->Stop("Sector track fitting");

    if(fWriteOut) 
       WriteSpacePoints(npoints, points, slice, patch); //do after Tracking

    //free memory
    if(UseCF)
      memory->Free();
    else
      fFileHandler->Free();

    UInt_t ntracks0 =0;
    AliHLTTrackSegmentData *trackdata0  = 
         (AliHLTTrackSegmentData *) memory->Allocate(fTracker->GetTracks());
    memory->TrackArray2Memory(ntracks0,trackdata0,fTracker->GetTracks());

    //write tracks
    if(fWriteOut){
      if(fNPatch==1)
	sprintf(name,"%s/tracks_tr_%d_%d_%d.raw",fWriteOutPath,fEvent,slice,-1);
      else
	sprintf(name,"%s/tracks_tr_%d_%d_%d.raw",fWriteOutPath,fEvent,slice,patch);
      memory->SetBinaryOutput(name);
      memory->Memory2Binary(ntracks0,trackdata0);
      memory->CloseBinaryOutput();
    }

    fInterMerger->Reset();
    fInterMerger->Init(fRow[patch],patch);

    fInterMerger->FillTracks(ntracks0,trackdata0);
    
    //fBenchmark->Start("Inter Merger");
    // fInterMerger->Merge();
    //    fInterMerger->SlowMerge();
    
    //fBenchmark->Stop("Inter Merger");
    /*
    //write inter merged tracks
    if(fWriteOut){
      sprintf(name,"%stracks_im_%d_%d.raw",fWriteOutPath,slice,patch);
      WriteTracks(name,fInterMerger,'i'); //write output of intermerger
      }
    */
    memory->Free();
    
    UInt_t ntracks1 =0;
    AliHLTTrackSegmentData *trackdata1 =
      (AliHLTTrackSegmentData *) memory->Allocate(fInterMerger->GetInTracks(0));
    memory->TrackArray2Memory(ntracks1,trackdata1,fInterMerger->GetInTracks(0));

    fTrackMerger->InitSector(slice,patch);
    fTrackMerger->FillTracks(ntracks1,trackdata1);

    memory->Free();
  }
  //fBenchmark->Start("Patch Merger");
  //fTrackMerger->SlowMerge();
  fTrackMerger->AddAllTracks();
  //fTrackMerger->Merge();
  //fBenchmark->Stop("Patch Merger");
  /*
  //write merged tracks
  if(fWriteOut){
    sprintf(name,"%stracks_tm_%d.raw",fWriteOutPath,slice);
    WriteTracks(name,fTrackMerger,'o'); //write output of trackmerger
  }
  */
  fTrackData = (AliHLTTrackSegmentData *) 
                         fFileHandler->Allocate(fTrackMerger->GetOutTracks());

  fFileHandler->TrackArray2Memory(fNTrackData,fTrackData,
                                                fTrackMerger->GetOutTracks());

  delete memory;
}

void AliLevel3::FitGlobalTracks()
{ 
  //fit global tracks
  AliHLTFitter *fitter = new AliHLTFitter(fVertex,AliLevel3::DoVertexFit());
  if(fNPatch==1)
    fitter->LoadClusters(fWriteOutPath,fEvent,kTRUE);
  else
    fitter->LoadClusters(fWriteOutPath,fEvent,kFALSE);
  
  fBenchmark->Start("Global track fitter");
  AliHLTTrackArray *tracks = fGlobalMerger->GetOutTracks();
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliHLTTrack *tr = tracks->GetCheckedTrack(i);
      if(!tr) continue;
      fitter->FitHelix(tr);
      tr->UpdateToFirstPoint();
    }
  fBenchmark->Stop("Global track fitter");
  delete fitter;
}

void AliLevel3::WriteSpacePoints(UInt_t npoints,AliHLTSpacePointData *points,
				 Int_t slice,Int_t patch) const
{ 
  //write space points
  char name[256];
  if(fNPatch==1)
    sprintf(name,"%s/points_%d_%d_%d.raw",fWriteOutPath,fEvent,slice,-1);
  else
    sprintf(name,"%s/points_%d_%d_%d.raw",fWriteOutPath,fEvent,slice,patch);
  AliHLTMemHandler * memory = new AliHLTMemHandler();
  memory->SetBinaryOutput(name);
  memory->Transform(npoints,points,slice);
  memory->Memory2Binary(npoints,points);
  memory->CloseBinaryOutput();
  delete  memory;
}

Int_t AliLevel3::WriteTracks(char *filename,AliHLTMerger *merger,char opt) const
{ 
  //write tracks
  AliHLTMemHandler *memory = new AliHLTMemHandler();
  memory->SetBinaryOutput(filename);
  if(opt=='a'||opt=='i'){  //add intracks
    for(Int_t i=0;i<merger->GetNIn();i++){
      AliHLTTrackArray *tr=merger->GetInTracks(i);
      memory->TrackArray2Binary(tr);
    }
  }

  if(opt=='o'||opt=='a'){
    AliHLTTrackArray *tr=merger->GetOutTracks();
    memory->TrackArray2Binary(tr);
  }

  memory->CloseBinaryOutput();
  delete memory;
  return 1;
}

void AliLevel3::WriteResults()
{
  //Write the resulting tracks to outputfile
  Char_t fname[256];
  sprintf(fname,"%s/tracks_%d.raw",fWriteOutPath,fEvent);
  WriteTracks(fname,fGlobalMerger,'a');
  //WriteTracks("tracks.raw",fGlobalMerger,'a');
  sprintf(fname,"%s/tracks_gl_%d.raw",fWriteOutPath,fEvent);
  WriteTracks(fname,fGlobalMerger,'o');
  //WriteTracks("tracks_gl.raw",fGlobalMerger,'o');
}
