//////////////////////////////////////////////////////////////////////////
//  Alice ITS first detector alignment program.                         //
//                                                                      //
// version: 0.0.0 Draft.                                                //
// Date: April 18 1999                                                  //
// By: Bjorn S. Nilsen                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <fstream.h>

// Data structure to hold averaged clusts.
struct ClustAl_sl{
    Int_t lay,lad,det;
    Float_t xg,yg,zg,xl,yl,zl;
};
struct ClustAl_tl{
    Int_t    track,nclust;  // track number and number of data points.
    ClustAl_sl *clust;      // data points to fit.
    Float_t  a,b,c,d,a0,b0,c0,d0,qual;  // fit parameters and fit quality.
    Float_t  px,py,pz,p,pt;
    // x=a+bz and y=c+dz;
    // x=a0+b0*y and z=c0+d0*y in coordinate system of clust[0].lay,lad,det
};

//
void AlignITSmacro3(const char *Rfilename="galice_ITS_B0.root",
		    const char *sfile="Align_ITS_B0",
		    Float_t sigma1=0.0,Float_t sigma2=0.0,Float_t sigma3=0.0,
		    Int_t evNumber=0) {
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro illustrating how to 
//   read the output of GALICE and fill some histograms.
//
//     Root > .L AlignITSmacro2.C    // this loads the macro in memory
//     Root > AlignITSmacro2();      // by default process first event   
//     Root > AlignITSmacro2("galice.root");   // process file galice.root
//     Root > AlignITSmacro2("galice.root",3); // process file galice.root
//                                                the first 4 events
//       or
//     Root > .x AlignITSmacro2.C;  //by default process first event
//     Root > .x AlignITSmacro2.C("galice.root");   // process file galice.root
//     Root > .x AlignITSmacro2.C("galice.root",3); // process file galice.root
//                                                     the first 4 events
/////////////////////////////////////////////////////////////////////////
//
    gROOT->Reset();  // Reset root to it's default state
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } // end if gClassTable...

   // Connect the Root Galice file containing Geometry, Kine and Clusts
   TFile *Rfile = (TFile*)gROOT->GetListOfFiles()->FindObject(Rfilename);
   if(!Rfile) Rfile = new TFile(Rfilename);
   printf("reading from file %s\n",Rfilename);

   // Get AliRun object from file or create it if not on file
   if(!gAlice) {
      gAlice = (AliRun*)Rfile->Get("gAlice");
      if( gAlice) printf("AliRun object found on file\n");
      if(!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   } /* end if gAlice */

   Float_t    v0[3] = {0.0,0.0,0.0};
   Float_t    tran[3] = {0.0,0.0,0.0},rot[3] = {0.0,0.0,0.0};
   Float_t    trans[15] ={0.0E-0,1.0E-4,4.0E-4,7.0E-4,1.0E-3,
			  2.0E-3,4.0E-3,6.0E-3,8.0E-3,1.0E-2,
			  2.0E-2,3.0E-2,5.0E-2,7.5E-2,1.0E-1}; // cm
   Float_t    rots[15] ={0.0E-0,1.0E-4,4.0E-4,7.0E-4,1.0E-3,
			 2.0E-3,4.0E-3,6.0E-3,8.0E-3,1.0E-2,
			 2.0E-2,3.0E-2,5.0E-2,7.5E-2,1.0E-1}; // rad
   Int_t      nparticles,Ntrkp,ntrk,Itimes,Isigmas;
   AliITS     *ITS = 0;
   TTree      *TH = 0;
   AliITSgeom *gm,gm2;
   char       Hfilename[80];
//   char       Gfilename[80];
   TFile      *Hfile = 0;
   ClustAl_tl *trk = 0;
   FILE       *fp;

   for(Int_t evnt=0;evnt<=evNumber;evnt++){
      // define some variables for later use.
      nparticles = gAlice->GetEvent(evnt);
      printf("nparticles %d\n",nparticles);
      if (nparticles <= 0) continue; /* get next event */

      // Get pointers to Alice detectors and Clusts containers
      ITS   = (AliITS*)gAlice->GetDetector("ITS");
      if(!ITS) return;          /* error no ITS data exit */
      TH    = gAlice->TreeH();
      Ntrkp = TH->GetEntries();
      gm    = ITS->GetITSgeom();

      // Array (stucture) of clusts for the first and second layer
      // this should be replaced with either clusters or digits
      // when they are proporly defined.
      trk = new ClustAl_tl[Ntrkp];

      printf("Ntrkp=%d\n",Ntrkp);

      HitsToClustAl(trk,ntrk,Ntrkp,TH,ITS,1.0);
      printf("Filled data structures ntrk=%d fitting lines next\n",ntrk);
//
      for(Itimes=0;Itimes<1;Itimes++){
	  if(Itimes==0){
	      fp = fopen();
            fprintf(fp,"Rr,Rs,Rrx1,Rerx1,Rsrx1,Rserx1,Rrz1,Rerz1,Rsrz1,Rserz1,"
		             "Rrx2,Rerx2,Rsrx2,Rserx2,Rrz2,Rerz2,Rsrz2,Rserz2,"
		             "Rrx3,Rerx3,Rsrx3,Rserx3,Rrz3,Rerz3,Rsrz3,Rserz3,"
	                     "Rrx4,Rerx4,Rsrx4,Rserx4,Rrz4,Rerz4,Rsrz4,Rserz4,"
	                     "Rrx5,Rerx5,Rsrx5,Rserx5,Rrz5,Rerz5,Rsrz5,Rserz5,"
                            "Rrx6,Rerx6,Rsrx6,Rserx6,Rrz6,Rerz6,Rsrz6,Rserz6");
	  }else if(Itimes==1){
	      fp = fopen();
            fprintf(fp,"Fr,Fs,Frx1,Ferx1,Fsrx1,Fserx1,Frz1,Ferz1,Fsrz1,Fserz1,"
		             "Frx2,Ferx2,Fsrx2,Fserx2,Frz2,Ferz2,Fsrz2,Fserz2,"
		             "Frx3,Ferx3,Fsrx3,Fserx3,Frz3,Ferz3,Fsrz3,Fserz3,"
	                     "Frx4,Ferx4,Fsrx4,Fserx4,Frz4,Ferz4,Fsrz4,Fserz4,"
	                     "Frx5,Ferx5,Fsrx5,Fserx5,Frz5,Ferz5,Fsrz5,Fserz5,"
                            "Frx6,Ferx6,Fsrx6,Fserx6,Frz6,Ferz6,Fsrz6,Fserz6");
	  }else if(Itimes==2){
	      fp = fopen();
            fprintf(fp,"Zr,Zs,Zrx1,Zerx1,Zsrx1,Zserx1,Zrz1,Zerz1,Zsrz1,Zserz1,"
		             "Zrx2,Zerx2,Zsrx2,Zserx2,Zrz2,Zerz2,Zsrz2,Zserz2,"
		             "Zrx3,Zerx3,Zsrx3,Zserx3,Zrz3,Zerz3,Zsrz3,Zserz3,"
	                     "Zrx4,Zerx4,Zsrx4,Zserx4,Zrz4,Zerz4,Zsrz4,Zserz4,"
	                     "Zrx5,Zerx5,Zsrx5,Zserx5,Zrz5,Zerz5,Zsrz5,Zserz5,"
                            "Zrx6,Zerx6,Zsrx6,Zserx6,Zrz6,Zerz6,Zsrz6,Zserz6");
	  }else if(Itimes==3){
	      fp = fopen();
            fprintf(fp,"Ar,As,Arx1,Aerx1,Asrx1,Aserx1,Arz1,Aerz1,Asrz1,Aserz1,"
		             "Arx2,Aerx2,Asrx2,Aserx2,Arz2,Aerz2,Asrz2,Aserz2,"
		             "Arx3,Aerx3,Asrx3,Aserx3,Arz3,Aerz3,Asrz3,Aserz3,"
	                     "Arx4,Aerx4,Asrx4,Aserx4,Arz4,Aerz4,Asrz4,Aserz4,"
	                     "Arx5,Aerx5,Asrx5,Aserx5,Arz5,Aerz5,Asrz5,Aserz5,"
                            "Arx6,Aerx6,Asrx6,Aserx6,Arz6,Aerz6,Asrz6,Aserz6");
	  }else if(Itimes==4){
	      fp = fopen();
            fprintf(fp,"Br,Bs,Brx1,Berx1,Bsrx1,Bserx1,Brz1,Berz1,Bsrz1,Bserz1,"
		             "Brx2,Berx2,Bsrx2,Bserx2,Brz2,Berz2,Bsrz2,Bserz2,"
		             "Brx3,Berx3,Bsrx3,Bserx3,Brz3,Berz3,Bsrz3,Bserz3,"
	                     "Brx4,Berx4,Bsrx4,Bserx4,Brz4,Berz4,Bsrz4,Bserz4,"
	                     "Brx5,Berx5,Bsrx5,Bserx5,Brz5,Berz5,Bsrz5,Bserz5,"
                            "Brx6,Berx6,Bsrx6,Bserx6,Brz6,Berz6,Bsrz6,Bserz6");
	  }else if(Itimes==5){
	      fp = fopen();
            fprintf(fp,"Cr,Cs,Crx1,Cerx1,Csrx1,Cserx1,Crz1,Cerz1,Csrz1,Cserz1,"
		             "Crx2,Cerx2,Csrx2,Cserx2,Crz2,Cerz2,Csrz2,Cserz2,"
		             "Crx3,Cerx3,Csrx3,Cserx3,Crz3,Cerz3,Csrz3,Cserz3,"
	                     "Crx4,Cerx4,Csrx4,Cserx4,Crz4,Cerz4,Csrz4,Cserz4,"
	                     "Crx5,Cerx5,Csrx5,Cserx5,Crz5,Cerz5,Csrz5,Cserz5,"
                            "Crx6,Cerx6,Csrx6,Cserx6,Crz6,Cerz6,Csrz6,Cserz6")
	  } // end if
	  for(Isigmas=0;Isigmas<11;Isigmas++){
//
//	  tran[0] = sigma1;
//	  tran[1] = sigma2;
//	  tran[2] = sigma3;
	  if(Itimes==0){ tran[0] = trans[Isigmas];
	  }else tran[0] = 0.0;
	  if(Itimes==1){ tran[1] = trans[Isigmas];
	  }else tran[1] = 0.0;
	  if(Itimes==2){ tran[2] = trans[Isigmas];
	  }else tran[2] = 0.0;
	  if(Itimes==3){ rot[0] = rots[Isigmas];
	  }else rot[0] = 0.0;
	  if(Itimes==4){ rot[1] = rots[Isigmas];
	  }else rot[1] = 0.0;
	  if(Itimes==5){ rot[2] = rots[Isigmas];
	  }else rot[2] = 0.0;
	  printf("tran= %e %e %e (cm), rot=%e %e %e (rad)\n",
		 tran[0],tran[1],tran[2],rot[0],rot[1],rot[2]);
//
	  gm2 = *gm;
	  gm2.RandomCylindericalChange(tran,rot);

	  FillGlobalPositions(trk,ntrk,&gm2);
	  // setup to save all created histograms.
	  sprintf(Hfilename,"%s_%04.0fr%04.0fp%04.0fz%04.0fx%04.0fy%04.0fz.root",
		  sfile,
		  10000.*tran[0],10000.*tran[1],10000.*tran[2],
		  10000.* rot[0],10000.* rot[1],10000.* rot[2]);
	  Hfile = (TFile*)gROOT->GetListOfFiles()->FindObject(Hfilename);
	  if(Hfile) Hfile->Close();
	  Hfile = new TFile(Hfilename,"RECREATE","Histograms from AlignITS.C");
	  printf("Histograms saved to file %s\n",Hfilename);
          //
	  PlotGeomChanges(gm,&gm2,Hfile);
          //
	  // fit all tracks and do a track quality hist.
	  FitAllTracks(trk,ntrk,v0,&gm2,sfile,Hfile);
//	  FitVertexAll(trk,ntrk,sfile,Hfile); 
	  // Find all 2 track vertecies and hist. them.
	  Hfile->Close();
          //
      } // end for Isigmas
	  fp = fclose();
      } // end for Itimes
//
      printf("Event %d done\n",evnt);
//
      deleteClustAl(trk,ntrk); // subrotine to delet memory allocated
                               // inside HitsToclustAl.
      delete[] trk;            // now delet memory allocated above.
   } // end for evnt
   Rfile->Close();
   return;
}
