//////////////////////////////////////////////////////////////////////////
//  Alice ITS first detector alignment program.                         //
//                                                                      //
// version: 0.0.0 Draft.                                                //
// Date: April 18 1999                                                  //
// By: Bjorn S. Nilsen                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

// Data structure to hold averaged clusts.
struct ClustAl_sl{
    Int_t lay,lad,det;
    Float_t xg,yg,zg,xl,yl,zl;
};
struct ClustAl_tl{
    Int_t    track,nclust;  // track number and number of data points.
    ClustAl_sl *clust;        // data points to fit.
    Float_t  a,b,c,d,qual;  // fit parameters and fit quality.
    Float_t  px,py,pz,p,pt;
    // x=a+bz and y=c+dz;
};

//
void AlignITSmacro (const char *filename="galice_ITS_B0.root",
		    const char *sfile="Align_ITS_B0",
		    Float_t sigma1=0.0,Float_t sigma2=0.0,Float_t sigma3=0.0,
		    Int_t evNumber=0) {
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//
//     Root > .L AlignITS.C    //this loads the macro in memory
//     Root > AlignITS();  //by default process first event   
//     Root > AlignITS("galice.root"); //process file galice.root
//     Root > AlignITS("galice.root",3); // process file galice.root
//                                          the first 4 eventst
//       or
//     Root > .x AlignITS.C;  //by default process first event   
//     Root > .x AlignITS.C("galice.root"); //process file galice.root
//     Root > .x AlignITS.C("galice.root",3); // process file galice.root
//                                               the first 4 eventst
/////////////////////////////////////////////////////////////////////////
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } // end if gClassTable...
      
   // Connect the Root Galice file containing Geometry, Kine and Clusts
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if(!file) file = new TFile(filename);
   printf("reading from file %s\n",filename);

   // Get AliRun object from file or create it if not on file
   if(!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if( gAlice) printf("AliRun object found on file\n");
      if(!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   } /* end if gAlice */

   for(Int_t evnt=0;evnt<=evNumber;evnt++){

      // define some variables for later use.
      Int_t      nparticles = gAlice->GetEvent(evnt);
      printf("nparticles %d\n",nparticles);
      if (nparticles <= 0) continue; /* get next event */

      // Get pointers to Alice detectors and Clusts containers
      AliITS *ITS    = (AliITS*)gAlice->GetDetector("ITS");
      if(!ITS) return;          /* error no ITS data exit */
      TTree  *TH     = gAlice->TreeH();
      Int_t  ntracks = TH->GetEntries();

      // Array (stucture) of clusts for the first and second layer
      // this should be replaced with either clusters or digits
      // when they are proporly defined.
      ClustAl_tl *trk = new ClustAl_tl[ntracks];
      Int_t      ntrk;
      Float_t    v0[3] = {0.0,0.0,0.0};

      printf("ntracks=%d\n",ntracks);

      HitsTo(trk,ntrk,ntracks,TH,ITS,sigma1,sigma2,sigma3,0.0,0.0,0.0);

      printf("Filled data structures ntrk=%d fitting lines next\n",ntrk);
//      return;
      // setup to save all created histograms.
      char Hfilename[80];
      sprintf(Hfilename,"%s.root",sfile);
      TFile *Hfile = (TFile*)gROOT->GetListOfFiles()->FindObject(Hfilename);
      if(Hfile) Hfile->Close();
      Hfile = new TFile(Hfilename,"RECREATE","Histograms from AlignITS.C");
      printf("Histograms saved to file %s\n",Hfilename);

      FitAllTracks(trk,ntrk,v0,sfile,Hfile); 
// fit all tracks and do a track quality hist.

      printf("Fitted tracks, finding vertex next\n");

      FitVertexAll(trk,ntrk,sfile,Hfile); 
// Find all 2 track vertecies and hist. them.

      printf("Event %d done\n",evnt);

      for(Int_t i=0;i<ntrk;i++) delete trk[i].clust;
      delete[] trk;
   } /* end for evnt */
   return;
}
