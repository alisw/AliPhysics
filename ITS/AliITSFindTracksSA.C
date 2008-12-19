#if !defined(__CINT__) || defined(__MAKECINT__) 
  #include <Riostream.h>
  #include "AliITSgeom.h"
  #include "AliITS.h"
  #include "AliITStrackerSA.h"
  #include "AliITSVertexerFast.h"
  #include "AliRun.h"
  #include "AliRunLoader.h"
//  #include "AliTPCLoader.h"
  #include "AliITSLoader.h"
  #include "TStopwatch.h"
  #include "AliMagF.h"
#endif

Int_t AliITSFindTracksSA(Int_t evin=0,Int_t nevents=1,char *opt="onlyITS+6/6",const Char_t *clusterFileName="clusters.root", const Char_t *tracksFileName="ITS.TracksSA.root") {

  //This macro finds tracks in the ITS Stand Alone and writes them in
  //the file ITS.TracksSA.root as tracks of class AliITStracksV2.

  //This macro needs both AliITSRecPoint (to find the vertex) and 
  //AliITSclusterV2 reconstructed points (track finding). Clusters V2
  //must be saved in a file with a different name from that of RecPoint.

  //Options: write onlyITS to track only with the ITS
  //         without the option onlyITS combined tracking TPC+ITS 
  //         and ITS stand-alone will be performed
  //
  //         write 6/6 to accept only tracks with 6 clusters
  //         write 5/6 to accept tracks with 5 clusters good over 6

   
   if (gAlice) {
      delete AliRunLoader::GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }
 
   AliRunLoader* rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"AliITSFindTracksSA.C : Can not open session RL=NULL"<< endl;
      return 3;
   }
   
   Int_t retval = rl->LoadgAlice();
   if (retval) {
      cerr<<"AliITSFindTracksSA.C : LoadgAlice returned error"<<endl;
      delete rl;
      return 3;
   }
   
   retval = rl->LoadHeader();
   if (retval) {
      cerr<<"AliITSFindTracksSA.C : LoadHeader returned error"<<endl;
      delete rl;
      return 3;
   }
   gAlice=rl->GetAliRun();
       
   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliITSFindTracksSA.C : Can not get ITS loader"<<endl;
      return 4;
   }
   
   AliITS *iTS = (AliITS*)gAlice->GetDetector("ITS");
   if (!iTS) {
      cerr<<"AliITSFindTracksSA.C : Can not find the ITS detector !"<<endl;
      return 6;
   }
   AliITSgeom *geom = iTS->GetITSgeom();
   AliKalmanTrack::SetConvConst(1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField());
   
   TString choice(opt);
   Bool_t onlyITS=choice.Contains("onlyITS");

   TStopwatch timer;
  
   for(Int_t iev=evin;iev<nevents;iev++){
     rl->GetEvent(iev);
     itsl->LoadRecPoints();

     //AliITSVertexerPPZ* vertexer = new AliITSVertexerPPZ("vertici.root");
     Double_t smear[3]={0.0150,0.0150,0.0150};
     AliITSVertexerFast* vertexer = new AliITSVertexerFast(smear);
     TTree* cltree = itsl->TreeR();
     AliITSVertex* vert = vertexer->FindVertexForCurrentEvent(cltree);
     AliITStrackerSA tracker(geom,vert);  
     tracker.SetEventNumber(iev);
     
     itsl->UnloadRecPoints();
     itsl->SetRecPointsFileName(clusterFileName);
     itsl->LoadRecPoints();
    
     if(onlyITS){
       itsl->SetTracksFileName(tracksFileName);
       itsl->LoadTracks("recreate");
              
       TTree* treec = (TTree*)itsl->TreeR();
       TTree *itsTree=itsl->TreeT();
       if (!itsTree) {
	 itsl->MakeTree("T");
	 itsTree=itsl->TreeT();
       }
            
       tracker.FindTracks(treec,itsTree,iev,opt);
       itsl->WriteTracks("OVERWRITE");
     } 
     if(!onlyITS){
       itsl->LoadTracks("read");
       TTree *treev2=(TTree*)itsl->TreeT();
       TTree* treec = (TTree*)itsl->TreeR();
       tracker.UseFoundTracksV2(iev,treev2,treec);
       itsl->UnloadTracks();
       itsl->SetTracksFileName(tracksFileName);
       itsl->LoadTracks("recreate");
       TTree *itsTree=itsl->TreeT();
       if (!itsTree) {
	 itsl->MakeTree("T");
	 itsTree=itsl->TreeT();
       }
       tracker.FindTracks(treec,itsTree,iev,opt);
       itsl->WriteTracks("OVERWRITE");
       
    }
     
   }
   timer.Stop(); timer.Print();   
   delete geom; 

  
   return 0;
}





