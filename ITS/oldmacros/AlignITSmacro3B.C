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
void AlignITSmacro3B(const char *Rfilename="galice_ITS_B0.root",
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
   Float_t    Rdta[6];
   Float_t    Fdta[24],Fdr[12],Fdta0[24],Fdr0[12];
   Int_t      Ndta[12];
   FILE       *fp1,*fp2;

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
      for(Itimes=4;Itimes<5;Itimes++){
	  if(Itimes==0){
	      fp1 = fopen("RvariationsR.csv","w");
	      fp2 = fopen("RvariationsN.csv","w");
            fprintf(fp1,"Rr,Rs,Rrx1,Rerx1,Rsrx1,Rserx1,Rrz1,Rerz1,Rsrz1,Rserz1,"
		    "Rrd1,Rerd1,Rsrd1,Resrd1,"
		             "Rrx2,Rerx2,Rsrx2,Rserx2,Rrz2,Rerz2,Rsrz2,Rserz2,"
		    "Rrd2,Rerd2,Rsrd2,Resrd2,"
		             "Rrx3,Rerx3,Rsrx3,Rserx3,Rrz3,Rerz3,Rsrz3,Rserz3,"
		    "Rrd3,Rerd3,Rsrd3,Resrd3,"
	                     "Rrx4,Rerx4,Rsrx4,Rserx4,Rrz4,Rerz4,Rsrz4,Rserz4,"
		    "Rrd4,Rerd4,Rsrd4,Resrd4,"
	                     "Rrx5,Rerx5,Rsrx5,Rserx5,Rrz5,Rerz5,Rsrz5,Rserz5,"
		    "Rrd5,Rerd5,Rsrd5,Resrd5,"
                            "Rrx6,Rerx6,Rsrx6,Rserx6,Rrz6,Rerz6,Rsrz6,Rserz6,"
		    "Rrd6,Rerd6,Rsrd6,Resrd6");
            fprintf(fp2,"Rr,Rs,RN1,RN2,RN3,RN4,RN5,RN6,RN7,RN8,RN9,RN10,"
		    "RN50,Rcut");
	  }else if(Itimes==1){
	      fp1 = fopen("RPhivariationsR.csv","w");
	      fp2 = fopen("RPhivariationsN.csv","w");
            fprintf(fp1,"Fr,Fs,Frx1,Ferx1,Fsrx1,Fserx1,Frz1,Ferz1,Fsrz1,Fserz1,"
		    "Frd1,Ferd1,Fsrd1,Fesrd1,"
		             "Frx2,Ferx2,Fsrx2,Fserx2,Frz2,Ferz2,Fsrz2,Fserz2,"
		    "Frd2,Ferd2,Fsrd2,Fesrd2,"
		             "Frx3,Ferx3,Fsrx3,Fserx3,Frz3,Ferz3,Fsrz3,Fserz3,"
		    "Frd3,Ferd3,Fsrd3,Fesrd3,"
	                     "Frx4,Ferx4,Fsrx4,Fserx4,Frz4,Ferz4,Fsrz4,Fserz4,"
		    "Frd4,Ferd4,Fsrd4,Fesrd4,"
	                     "Frx5,Ferx5,Fsrx5,Fserx5,Frz5,Ferz5,Fsrz5,Fserz5,"
		    "Frd5,Ferd5,Fsrd5,Fesrd5,"
                            "Frx6,Ferx6,Fsrx6,Fserx6,Frz6,Ferz6,Fsrz6,Fserz6,"
		    "Frd6,Ferd6,Fsrd6,Fesrd6");
            fprintf(fp2,"Fr,Fs,FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,"
		    "FN50,Fcut");
	  }else if(Itimes==2){
	      fp1 = fopen("ZvariationsR.csv","w");
	      fp2 = fopen("ZvariationsN.csv","w");
            fprintf(fp1,"Zr,Zs,Zrx1,Zerx1,Zsrx1,Zserx1,Zrz1,Zerz1,Zsrz1,Zserz1,"
		    "Zrd1,Zerd1,Zsrd1,Zesrd1,"
		             "Zrx2,Zerx2,Zsrx2,Zserx2,Zrz2,Zerz2,Zsrz2,Zserz2,"
		    "Zrd2,Zerd2,Zsrd2,Zesrd2,"
		             "Zrx3,Zerx3,Zsrx3,Zserx3,Zrz3,Zerz3,Zsrz3,Zserz3,"
		    "Zrd3,Zerd3,Zsrd3,Zesrd3,"
	                     "Zrx4,Zerx4,Zsrx4,Zserx4,Zrz4,Zerz4,Zsrz4,Zserz4,"
		    "Zrd4,Zerd4,Zsrd4,Zesrd4,"
	                     "Zrx5,Zerx5,Zsrx5,Zserx5,Zrz5,Zerz5,Zsrz5,Zserz5,"
		    "Zrd5,Zerd5,Zsrd5,Zesrd5,"
                            "Zrx6,Zerx6,Zsrx6,Zserx6,Zrz6,Zerz6,Zsrz6,Zserz6,"
		    "Zrd6,Zerd6,Zsrd6,Zesrd6");
            fprintf(fp2,"Zr,Zs,ZN1,ZN2,ZN3,ZN4,ZN5,ZN6,ZN7,ZN8,ZN9,ZN10,"
		    "ZN50,Zcut");
	  }else if(Itimes==3){
	      fp1 = fopen("AvariationsR.csv","w");
	      fp2 = fopen("AvariationsN.csv","w");
            fprintf(fp1,"Ar,As,Arx1,Aerx1,Asrx1,Aserx1,Arz1,Aerz1,Asrz1,Aserz1,"
		    "Ard1,Aerd1,Asrd1,Aesrd1,"
		             "Arx2,Aerx2,Asrx2,Aserx2,Arz2,Aerz2,Asrz2,Aserz2,"
		    "Ard2,Aerd2,Asrd2,Aesrd2,"
		             "Arx3,Aerx3,Asrx3,Aserx3,Arz3,Aerz3,Asrz3,Aserz3,"
		    "Ard3,Aerd3,Asrd3,Aesrd3,"
	                     "Arx4,Aerx4,Asrx4,Aserx4,Arz4,Aerz4,Asrz4,Aserz4,"
		    "Ard4,Aerd4,Asrd4,Aesrd4,"
	                     "Arx5,Aerx5,Asrx5,Aserx5,Arz5,Aerz5,Asrz5,Aserz5,"
		    "Ard5,Aerd5,Asrd5,Aesrd5,"
                            "Arx6,Aerx6,Asrx6,Aserx6,Arz6,Aerz6,Asrz6,Aserz6,"
		    "Ard6,Aerd6,Asrd6,Aesrd6");
            fprintf(fp2,"Ar,As,AN1,AN2,AN3,AN4,AN5,AN6,AN7,AN8,AN9,AN10,"
		    "AN50,Acut");
	  }else if(Itimes==4){
	      fp1 = fopen("BvariationsR.csv","w");
	      fp2 = fopen("BvariationsN.csv","w");
            fprintf(fp1,"Br,Bs,Brx1,Berx1,Bsrx1,Bserx1,Brz1,Berz1,Bsrz1,Bserz1,"
		    "Brd1,Berd1,Bsrd1,Besrd1,"
		             "Brx2,Berx2,Bsrx2,Bserx2,Brz2,Berz2,Bsrz2,Bserz2,"
		    "Brd2,Berd2,Bsrd2,Besrd2,"
		             "Brx3,Berx3,Bsrx3,Bserx3,Brz3,Berz3,Bsrz3,Bserz3,"
		    "Brd3,Berd3,Bsrd3,Besrd3,"
	                     "Brx4,Berx4,Bsrx4,Bserx4,Brz4,Berz4,Bsrz4,Bserz4,"
		    "Brd4,Berd4,Bsrd4,Besrd4,"
	                     "Brx5,Berx5,Bsrx5,Bserx5,Brz5,Berz5,Bsrz5,Bserz5,"
		    "Brd5,Berd5,Bsrd5,Besrd5,"
                            "Brx6,Berx6,Bsrx6,Bserx6,Brz6,Berz6,Bsrz6,Bserz6,"
		    "Brd6,Berd6,Bsrd6,Besrd6");
            fprintf(fp2,"Br,Bs,BN1,BN2,BN3,BN4,BN5,BN6,BN7,BN8,BN9,BN10,"
		    "BN50,Bcut");
	  }else if(Itimes==5){
	      fp1 = fopen("CvariationsR.csv","w");
	      fp2 = fopen("CvariationsN.csv","w");
            fprintf(fp1,"Cr,Cs,Crx1,Cerx1,Csrx1,Cserx1,Crz1,Cerz1,Csrz1,Cserz1,"
		    "Crd1,Cerd1,Csrd1,Cesrd1,"
		             "Crx2,Cerx2,Csrx2,Cserx2,Crz2,Cerz2,Csrz2,Cserz2,"
		    "Crd2,Cerd2,Csrd2,Cesrd2,"
		             "Crx3,Cerx3,Csrx3,Cserx3,Crz3,Cerz3,Csrz3,Cserz3,"
		    "Crd3,Cerd3,Csrd3,Cesrd3,"
	                     "Crx4,Cerx4,Csrx4,Cserx4,Crz4,Cerz4,Csrz4,Cserz4,"
		    "Crd4,Cerd4,Csrd4,Cesrd4,"
	                     "Crx5,Cerx5,Csrx5,Cserx5,Crz5,Cerz5,Csrz5,Cserz5,"
		    "Crd5,Cerd5,Csrd5,Cesrd5,"
		    "Crx6,Cerx6,Csrx6,Cserx6,Crz6,Cerz6,Csrz6,Cserz6,"
		    "Crd6,Cerd6,Csrd6,Cesrd6");
            fprintf(fp2,"Cr,Cs,CN1,CN2,CN3,CN4,CN5,CN6,CN7,CN8,CN9,CN10,"
		    "CN50,Ccut");
	  } // end if
	  fprintf(fp1,"\n");
	  fprintf(fp2,"\n");
	  for(Isigmas=0;Isigmas<15;Isigmas++){
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
	  PlotGeomChanges(gm,&gm2,Hfile,Rdta);
          //
	  // fit all tracks and do a track quality hist.
	  FitAllTracks(trk,ntrk,v0,&gm2,sfile,Hfile,Fdta,Ndta); 
//
	  for(Int_t l=0;l<6;l++){
	      Fdr[2*l+0] = sqrt(Fdta[4*l+0]*Fdta[4*l+0] +
				Fdta[4*l+2]*Fdta[4*l+2]);
	      Fdr[2*l+1] = sqrt(Fdta[4*l+0]*Fdta[4*l+0]*
				Fdta[4*l+1]*Fdta[4*l+1] +
				Fdta[4*l+2]*Fdta[4*l+2]*
				Fdta[4*l+3]*Fdta[4*l+3])/Fdr[2*l+0];
	  } // end for l
	  if(Isigmas==0){
	      for(Int_t fp=0;fp<24;fp++){
		  Fdta0[fp] = Fdta[fp];
		  if(fp<12) Fdr0[fp]  = Fdr[fp];
	      } // end for fp
	  } // end if Itimes==0&&Isigmas==0
	  if(Itimes==0) fprintf(fp1,"%e,%e,",tran[0],Rdta[0]);
	  else if(Itimes==1) fprintf(fp1,"%e,%e,",tran[1],Rdta[1]);
	  else if(Itimes==2) fprintf(fp1,"%e,%e,",tran[2],Rdta[2]);
	  else if(Itimes==3) fprintf(fp1,"%e,%e,",rot[0],Rdta[3]);
	  else if(Itimes==4) fprintf(fp1,"%e,%e,",rot[1],Rdta[4]);
	  else if(Itimes==5) fprintf(fp1,"%e,%e,",rot[2],Rdta[5]);
	  if(Itimes==0) fprintf(fp2,"%e,%e,",tran[0],Rdta[0]);
	  else if(Itimes==1) fprintf(fp2,"%e,%e,",tran[1],Rdta[1]);
	  else if(Itimes==2) fprintf(fp2,"%e,%e,",tran[2],Rdta[2]);
	  else if(Itimes==3) fprintf(fp2,"%e,%e,",rot[0],Rdta[3]);
	  else if(Itimes==4) fprintf(fp2,"%e,%e,",rot[1],Rdta[4]);
	  else if(Itimes==5) fprintf(fp2,"%e,%e,",rot[2],Rdta[5]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[0],Fdta[1],Fdta[0]/Fdta0[0],Fdta[1]/Fdta0[0]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[2],Fdta[3],Fdta[2]/Fdta0[2],Fdta[3]/Fdta0[2]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdr[0],Fdr[1],Fdr[0]/Fdr0[0],Fdr[1]/Fdr0[0]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[4],Fdta[5],Fdta[4]/Fdta0[4],Fdta[5]/Fdta0[4]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[6],Fdta[7],Fdta[6]/Fdta0[6],Fdta[7]/Fdta0[6]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdr[2],Fdr[3],Fdr[2]/Fdr0[2],Fdr[3]/Fdr0[2]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[8],Fdta[9],Fdta[8]/Fdta0[8],Fdta[9]/Fdta0[8]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[10],Fdta[11],Fdta[10]/Fdta0[10],Fdta[11]/Fdta0[10]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdr[4],Fdr[5],Fdr[4]/Fdr0[4],Fdr[5]/Fdr0[4]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[12],Fdta[12],Fdta[12]/Fdta0[12],Fdta[13]/Fdta0[12]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[14],Fdta[13],Fdta[14]/Fdta0[14],Fdta[15]/Fdta0[14]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdr[6],Fdr[7],Fdr[6]/Fdr0[6],Fdr[7]/Fdr0[6]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[16],Fdta[15],Fdta[16]/Fdta0[16],Fdta[17]/Fdta0[16]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[18],Fdta[17],Fdta[18]/Fdta0[18],Fdta[19]/Fdta0[18]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdr[8],Fdr[9],Fdr[8]/Fdr0[8],Fdr[9]/Fdr0[8]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[20],Fdta[19],Fdta[20]/Fdta0[20],Fdta[21]/Fdta0[20]);
	  fprintf(fp1,"%e,%e,%e,%e,",Fdta[22],Fdta[21],Fdta[22]/Fdta0[22],Fdta[23]/Fdta0[22]);
	  fprintf(fp1,"%e,%e,%e,%e\n",Fdr[10],Fdr[11],Fdr[10]/Fdr0[10],Fdr[11]/Fdr0[10]);
	  fprintf(fp2,"%d,%d,%d,%d,",Ndta[0],Ndta[1],Ndta[2],Ndta[3]);
	  fprintf(fp2,"%d,%d,%d,%d,",Ndta[4],Ndta[5],Ndta[6],Ndta[7]);
	  fprintf(fp2,"%d,%d,%d,%d\n",Ndta[8],Ndta[9],Ndta[10],Ndta[11]);
//	  FitVertexAll(trk,ntrk,sfile,Hfile); 
	  // Find all 2 track vertecies and hist. them.
	  Hfile->Close();
          //
      } // end for Isigmas
	  fclose(fp1);
	  fclose(fp2);
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
