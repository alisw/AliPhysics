#include "PMD_dig_reco.h"
#include <iostream.h> // define cout stream
#include <fstream.h>  // defines ofstream class
#include <stdlib.h>   // defines exit() functions
void PMD_dig_reco (Int_t evNumber=1) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//   
//
//  Created  by : Baba Potukuchi.V.K.S. ( JU / CERN ) 17/10/2000.  
//  Modified by : Baba Potukuchi.V.K.S. ( JU / CERN )  6/02/2001.
//  e-mail : potukuchi.baba@cern.ch 
//
/////////////////////////////////////////////////////////////////////////
// =========== Variables for Reconstruction==============================
  FILE *fp4;
  PMD_dig_reco* PPP = new PMD_dig_reco;
// ======================================================================
  Int_t padmult,cpvmult,padmult_gam,tot_trk,tot_trk_cpv;
  Int_t padmult_cut, cpvmult_cut;
  Int_t trk_Mult[100000];
  Int_t trk_Mult_cpv[100000];
  Float_t xPos,yPos,zPos,Edep,Vol1,Vol2,Vol2s,Vol3,Vol4,Vol5;
  Float_t  xpad,ypad;
  Int_t iVol4,iVol5 , ixpad,iypad;
  Int_t nbytes = 0;
  //  Float_t px,py,pz;
  Float_t px1,py1,pz1;
  Float_t mass;
  Float_t edep_tot;
  Int_t id_part;
  Int_t gam_in;
  Int_t jjj;
  Int_t pip_in,pim_in;
  TObjArray  *points;
  long double r,eta,px,py,pz;
  Float_t phi,Eparticle,theta1,theta;
  FILE *fp1;
  FILE *fp2;
  
  fp1=fopen("kine_vertex.out","w");
  fp2=fopen("kine_pmdhit.out","w");
  fp4=fopen("det1.out","w");

// Dynamical link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }
      
// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file
//   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
//   }
   // Just for test
//Loop over the requested number of events
  for (Int_t iEvt = 0; iEvt < evNumber; iEvt++) {

  printf("evt no %d \n",iEvt); 
   // End Just for test



// Import the Kine and Hits Trees for the event evNumber in the file
   Int_t nparticles = gAlice->GetEvent(iEvt);

    printf ("N of particles %d\n",nparticles);

    if (nparticles <= 0) return;

        Int_t nevt = iEvt + 1;

//         ----********---- vertex track search BEGIN

//Loop over the requested number of events
// Import the Kine and Hits Trees for the event evNumber in the file

	TObjArray *Particles = gAlice->Particles();
        Float_t phi;
   Float_t pl,cutmax,etamin,etamax,pmom,smin,smax;
	Float_t *pxyz;
	TObjArray  *points;
	TParticle *particle;
	AliPoints  *pm;
	//     for (Int_t it = 0; it < 100000; it++)trk_Mult[it]=0.; 

     for (Int_t in_track = 0; in_track < nparticles; in_track++) {
     //     for (Int_t in_track = 0; in_track < 10; in_track++) {
         particle = (TParticle*)gAlice->Particle(in_track);
         if (!particle) continue;
         id_part  = particle->GetPdgCode();
	 //      mass = TDatabasePDG::Instance()->GetParticle(id_part)->Mass();
	 //          printf("_____________--------%d %f\n",id_part,mass);

         px  = particle->Px();
         py  = particle->Py();
         pz  = particle->Pz();
         if(px !=0 || py!=0 || pz!=0) {
         if( (px!=0 || py!=0) && pz < 0 ){
//	   if(id_part == 22)gam_in++;
	 //	   printf("%d\n",gam_in);         
          r     = TMath::Sqrt(px*px + py*py);
	 theta = TMath::ATan2(r,TMath::Abs(pz));
	 phi = TMath::ATan2(py,px);
	 phi = (phi*180.)/3.14;
//         if(theta < -1.5707961)printf( " THETA=%f\n",theta) ;
         if(theta > -1.5707961)eta = -TMath::Log(TMath::Tan(0.5*theta));
         theta1=theta*180./3.14159;    

//         if(TMath::Abs(eta) > 2.3 && TMath::Abs(eta) < 3.5){
	 if(id_part == 22)gam_in++;
         if(id_part == 211)pip_in++;
	 if(id_part == -211)pim_in++;
     fprintf(fp1," %d %d %f %f %f %f %f %f \n",nevt,id_part,px,py,pz,eta,theta1,phi);
//            }
//	   printf("%d %d %d\n",gam_in,pip_in,pim_in);         
//          printf("%f %f %f %f %d\n",ptg,pxg,pyg,pzg,id_part);
     }
    }
  }
        
// 1999 file-- Get pointers to Alice detectors and Hits containers
     //------********----- vertex track ENDs

   Float_t x,y,z,mass,e;
   Int_t nbytes = 0;
   Int_t j,hit,ipart;
   Int_t npmd;
   Int_t mbytes = 0;
   Int_t nhits,gam_in;
   Int_t ipmd,id_part;
   Int_t ypad1;
   Int_t mpmd;
   const Int_t Npad2 = 72;

   Int_t dVol;
   //   Int_t iEvt;

   Int_t sector,plane;
   TParticle *particle;
   AliPMDhit  *pmdHit;

// Get pointers to Alice detectors and Hits containers
   AliDetector *PMD  = gAlice->GetDetector("PMD");
//   TClonesArray *fParticles = gAlice->Particles();
   if (PMD) TClonesArray *PMDhits   = PMD->Hits();
   if (PMD) TObjArray *points    = PMD->Points();

   TTree *TH = gAlice->TreeH();
   //   TTree *TK = gAlice->TreeK();

   Int_t ntracks    = TH->GetEntries();
   printf("No.of Tracksin the TreeH %d \n", ntracks);

// ===========   Initializing the array ============================
 for (Int_t iDet = 1; iDet <= 2 ; ++iDet) {
      printf(" Initializing for Detector %d \n",iDet) ;
   for (Int_t iBox = 1; iBox <= 30 ; ++iBox) {
     for (Int_t iCol = 1; iCol <= 72 ; ++iCol) {
        for (Int_t iRow = 1; iRow <= 72 ; ++iRow) {
       if(iDet==1)(*PPP).module1 [iBox] [iCol] [iRow] = 0.0 ;
       if(iDet==2)(*PPP).module2 [iBox] [iCol] [iRow] = 0.0 ;
        }
     }
   }
 }
//  Creating histograms .
   TH1F *edp = new TH1F("edp","Energy Deposited in each cell" ,100,100.,100000.);
   TH1F *Modnum = new TH1F("Modnum","Module number" ,100,1.,100.);
   TH2F *xycells = new TH2F("xycells","XY-of Cells in an Unit Modlule " ,100,0.,80.,100,0.,80.);
   TH2F *xyhits = new TH2F("xyhits","Hits of the cells " ,100,-150.,150.,100,-150.,150.);

//  ================================================================
    
// =======================================================
// Start loop on tracks in the hits containers
   for (Int_t track=0; track<ntracks;track++) {
        gAlice->ResetHits();

     //changes introduced ******

    nbytes += TH->GetEvent(track);
    //   mbytes += TK->GetEvent(track);
    //    printf("%d ******** %d \n" ,track,nbytes);
     // ======>Histogram PMD

     if (PMD) {
       npmd = PMDhits->GetEntriesFast();
       //printf("npmd=%d,  for track=%d %d \n" ,npmd,track,nbytes);
         for (ipmd = 0; ipmd < npmd; ipmd++) {
	 pmdHit = (AliPMDhit*) PMDhits->UncheckedAt(ipmd);
	 ipart = pmdHit->GetTrack();
      
      //      printf("ID part and hits on PMD --------%d %d \n" ,ipart,track);      
//  get kinematics of the particles
//	 particle = (TParticle*)fParticles->UncheckedAt(ipart);
          particle = gAlice->Particle(ipart);
//	 if (!particle) continue; 
	 //
	 id_part  = particle->GetPdgCode();
	 mass = TDatabasePDG::Instance()->GetParticle(id_part)->Mass();
	 //	 printf("ID part and hit on PMD %d %d \n" ,id_part,ipmd);      
	 px  = particle->Px();
	 py  = particle->Py();
	 pz  = particle->Pz();
	 Eparticle  = particle->Energy();
	 
	 r     = TMath::Sqrt(px*px + py*py);
	 theta = TMath::ATan2(r,TMath::Abs(pz));
	 //phi = TMath::ATan2(TMath::Abs(py),TMath::Abs(px));
	 phi = TMath::ATan2(py,px);
	 phi = (phi*180.)/3.14159;
	 eta = -TMath::Log(TMath::Tan(0.5*theta));
         if (pz < 0) eta = -eta;
	 theta=theta*180./3.14159;
	 
//         if(TMath::Abs(eta) > 2.3 && TMath::Abs(eta) < 3.5){
	 //
	 xPos = pmdHit->X();
	 yPos = pmdHit->Y();
	 zPos = pmdHit->Z();
	 Vol1 = pmdHit->fVolume[0];
	 Vol2 = pmdHit->fVolume[1];
	 Vol3 = pmdHit->fVolume[2];
	 Vol4 = pmdHit->fVolume[3];
	 Vol5 = pmdHit->fVolume[4];
         	 
	 Vol2s = Vol2-1;
	 ypad = Vol2s/Npad2 + 1;
	 ypad1 = ypad-1;
	 xpad = Vol2-ypad1*Npad2;
	 Edep = pmdHit->fEnergy;
         ixpad = xpad+0.2 ;
         iypad = ypad+0.2 ;
         iVol4 = Vol4+0.2 ;
         iVol5 = Vol5+0.2 ;
//         module [Vol4] [Vol5] [xpad] [ypad] = Edep ;
	 if(Edep>0)printf(" VOL4 %d %d %d %d %f \n",iVol4,iVol5,ixpad,iypad,Edep) ; 
         edp ->Fill(Edep) ;
         Modnum ->Fill(iVol5) ;
         xycells->Fill(ixpad,iypad);
         xyhits->Fill(xPos,yPos);
        
       if(iVol4==1)(*PPP).module1 [iVol5] [ixpad] [iypad] = Edep  ;
       if(iVol4==2)(*PPP).module2 [iVol5] [ixpad] [iypad] = Edep  ;
          Float_t XXP = (*PPP).module2 [iVol5] [ixpad] [iypad] ;
	 if(Edep>0)printf(" VOL4 %d %d %d %d %f \n",Vol4,Vol5,xpad,ypad,XXP) ; 
	 // Write to the kine.out file
	 fprintf(fp2," %d %d %d %f %f %f %d %d %d %d %f %f %f %f %f \n",
		nevt,ipart,id_part,xPos,yPos,zPos,xpad,ypad,Vol5,Vol4,Edep,px,py,pz,Eparticle);

//       Closing Eta loop of the incident particles.  
//         }
	 
       }
     }
   }
  }   

  pmdreco_(PPP);

// TCanvas *c1 = new TCanvas("c1","PMD Simulated digits(Edep) ",400,10,600,700);
// c1->Divide(1,2);
// c1->cd(1);
// gPad->SetFillColor(33);
// edp->SetFillColor(42);
// edp->Draw();
// c1->cd(2);
// gPad->SetFillColor(33);
// Modnum->SetFillColor(41);
// Modnum->Draw();
// c1->Print("baba.ps");
// TCanvas *c2 = new TCanvas("c2","PMD Simulated digits(Edep) ");
// gPad->SetFillColor(33);
// xycells->SetFillColor(42);
// xycells->Draw();
// c2->Print("baba1.ps");
// TCanvas *c3 = new TCanvas("c3","PMD Simulated digits(Edep) ");
// gPad->SetFillColor(33);
// xyhits->SetFillColor(42);
// xyhits->Draw();
// c3->Print("baba2.ps");
}

// ===========  RECONSTRUCTION ======================================
void pmdreco_(PMD_dig_reco* PPP) {
    Float_t  sqrth;
    sqrth = sqrt(3.) / 2.;
    
    printf(" I am in MAIN calling RanMAR & Intpts_ \n") ;
    rmarin_(PPP);
    intpts_(PPP);
//  Generate COORDINATES .
        for (Int_t i__ = 1; i__ <= 72; ++i__) {
            for (Int_t j = 1; j <= 72; ++j) {
                (*PPP).coord [1] [i__] [j] = i__ + j / 2.;
                (*PPP).coord [2] [i__] [j] = j * sqrth;
            }
        }
//-------------------------
  for (Int_t iDet = 1; iDet <= 2 ; ++iDet) {
       printf("idet= %d \n",iDet) ;
    for (Int_t iBox = 1; iBox <= 30 ; ++iBox) {
      for (Int_t iCol = 1; iCol <= 72 ; ++iCol) {
         for (Int_t iRow = 1; iRow <= 72 ; ++iRow) {
          Float_t PVAL = (*PPP).module1 [iBox] [iCol] [iRow];
//       printf("Det1 = %f %f %f %f %f \n",iDet, iBox, iCol, iRow,PVAL ) ;
          if(iDet==1)(*PPP).det1 [iCol] [iRow] = (*PPP).module1 [iBox] [iCol] [iRow] ;
          if(iDet==2)(*PPP).det2 [iCol] [iRow] = (*PPP).module2 [iBox] [iCol] [iRow] ;
       if(PVAL>0.0)printf("Det1 = %d %d %d %d %f \n",iDet, iBox, iCol, iRow,PVAL ) ;
//       fprintf(fp4," Det1 = %d %d %d %d %f \n",iDet, iBox, iCol, iRow,VAL ) ;

         }
      }
// Ending Reading pads in the Box .
    }
// Ending Reading pads for 1-CPV,2-PMD.
 }
/*     analyse detector 1 */
/*     CPV clustering. Ordering data */
       printf(" I am in MAIN Doing Clustering for CPV.  \n") ;
       printf (" I am in MAIN calling ==order== for CPV(det=1) \n");
       order_(1,PPP);
    printf(" I am in MAIN  after order for CPV(det=1) \n") ;
    printf(" I am in MAIN calling ==crclust== \n") ;
        crclust_(1, PPP) ;
    printf(" I am in MAIN after ==crclust==  & now call prints \n") ;
        prints_(1,PPP) ;
       printf(" I am in MAIN Doing Clustering for PMD.  \n") ;
printf (" I am in MAIN calling ==order== for PMD(det=2) \n");
        order_(2,PPP);
    printf(" I am in MAIN  after order for PMD(det=2) \n") ;
    printf(" I am in MAIN calling ==crclust== \n") ;
        crclust_(2, PPP) ;
    printf(" I am in MAIN after ==crclust==  & now call prints \n") ;
        prints_(2,PPP) ;
	(*PPP).incr = 0;
	for (Int_t i__ = 1; i__ <= 72; ++i__) {
	    for (j = 1; j <= 72; ++j) {
		if ((*PPP).infocl[1] [i__] [j]  == 1) {
		    ++(*PPP).incr;
		    (*PPP).inford[1] [(*PPP).incr] = (*PPP).infocl[2] [i__] [j] ;
		    (*PPP).inford[2] [(*PPP).incr] = i__;
		    (*PPP).inford[3] [(*PPP).incr] = j;
		}
	    }
	}
	for (i__ = 1; i__ <= 72; ++i__) {
	    for (j = 1; j <= 72; ++j) {
		if ((*PPP).infocl[1] [i__] [j]  == 2) {
		    ++(*PPP).incr;
		    (*PPP).inford[1] [(*PPP).incr] = (*PPP).infocl[2] [i__] [j];
		    (*PPP).inford[2] [(*PPP).incr] = i__;
		    (*PPP).inford[3] [(*PPP).incr] = j;
		}
	    }
	}
	arrange_(PPP);
/*     better clustering ( Gaussian fit ) */
	refclust_(PPP);
    return ;
} /* LOOP pmdreco */
// 
 /* >>>>>>>>>>>>  */void rmarin_(PMD_dig_reco* PPP)
{
    Int_t i__, j, k, l, m;
    Float_t s, t;
/*         Initializing routine for RANMAR, must be called before */
/*         generating any pseudorandom numbers with RANMAR. The input */
/*         values should be in the ranges:  0<=IJ<=31328 */
/*                                          0<=KL<=30081 */
/* To get the standard values in Marsaglia's paper, IJ=1802, KL=9373 */
    Int_t ij = 8775;
    Int_t kl = 5926;
    if (ij <= 0 || ij >= 31328) {
        goto L900;
    }
    if (kl <= 0 || kl >= 30081) {
        goto L900;
    }
 printf("I am in RanMAR \n");
    i__ = ij / 177 % 177 + 2;
    j = ij % 177 + 2;
    k = kl / 169 % 178 + 1;
    l = kl % 169;
    printf(" ranmar initialized: %d %d %d %d %d %d \n",ij,kl,i__,j,k,l);
    for (Int_t ii = 1; ii <= 97; ++ii) {
        s = 0.;
        t = .5;
        for (Int_t jj = 1; jj <= 24; ++jj) {
            m = (((i__ * j) % 179) * k) % 179;
            i__ = j;
            j = k;
            k = m ;
            l = l * 53 + 1 % 169;
            if (l * m % 64 >= 32) {
                s += t;
            }
            t *= .5;
        }
        (*PPP).u[ii - 1] = s;
    }
    (*PPP).c__ = .021602869033813477;
    (*PPP).cd  = .45623308420181274;
    (*PPP).cm  = .99999982118606567;
     printf(" Exiting from rmarin %f %f %f \n",(*PPP).c__,(*PPP).cd,(*PPP).cm) ;
     return ;
L900:
    printf( " +++++++ stop in rmarin: wrong params! ");
} /* rmarin_ */

/*>>>>>>>>>>>>>>>> */void order_(Int_t jdet,PMD_dig_reco* PPP)
{
    Int_t nmx,i__1, i__2, i__3;

    Float_t  adum;
    Int_t idum, i__, j, k, l;
    Float_t  d1[5185];
    Float_t  iord1[5185];
    Int_t itest, i1, i2, id ;

    printf("I am in order for DET=%d \n",jdet);
    nmx = 5184;
    for (i1 = 1; i1 <= 72; ++i1) {
	for (i2 = 1; i2 <= 72; ++i2) {
	    i__ = (i2 - 1) * 72 + i1;
            d1[i__] = 0.0 ;
	   iord1 [i__] = 0;
	   if(jdet==1) d1[i__] = (*PPP).det1 [i1] [i2] ;
	   if(jdet==2) d1[i__] = (*PPP).det2 [i1] [i2] ;
	   iord1 [i__] = i__;
	}
    }
    for (j = 2; j <= nmx; ++j) {
	itest = 0;
	adum = d1[j];
	idum = iord1 [j];
	i__2 = j - 1;
	for (k = 1; k <= i__2; ++k) {
	    if (adum > d1[k ] && itest == 0) {
		itest = 1 ; 
		i__3 = k;
		for (l = j - 1; l >= i__3; --l) {
		    d1[l+1] = d1[l];
		    iord1 [l+1] = iord1 [l];
		}
		d1[k] = adum;
		iord1 [k] = idum;
	    }
	}
    }
    i__1 = nmx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	id = iord1 [i__];
	i2 = id / 72 + 1;
	i1 = id - (i2 - 1) * 72;
	(*PPP).iorder[1] [i__] = i1;
	(*PPP).iorder[2] [i__] = i2;
     
     printf("I am in order for DET=%d %d %d %d %d %d \n",jdet,i__,i1,i2,(*PPP).iorder[1] [i__],(*PPP).iorder[2] [i__]);
    }
} /* order_ */

/* crclust_ */ void crclust_(Int_t ndet, PMD_dig_reco* PPP)
 { 
           Float_t d__ [73] [73] ;

    Int_t neibx[7] =  { 0,1,0,-1,-1,0,1 };
    Int_t neiby[7] =  { 0,0,1,1,0,-1,-1 };
    
    Int_t i__1;
    Float_t sqrth ;

    Int_t icld, idum, i__, j, k;
    Int_t nb, id1, id2, jd1, jd2, icle, icn, itr, nmx ;

   printf(" neibx %d %d %d %d %d %d \n",neibx[1],neibx[2],neibx[3],neibx[4],neibx[5],neibx[6]) ;
   printf(" neiby %d %d %d %d %d %d \n",neiby[1],neiby[2],neiby[3],neiby[4],neiby[5],neiby[6]) ;
  
	for (j = 1; j <= 72; ++j) {
	    for (k = 1; k <= 72; ++k) {
             d__ [j] [k] = 0.0 ;
             if(ndet==1) d__ [j] [k] = (*PPP).det1 [j] [k] ;
             if(ndet==2) d__ [j] [k] = (*PPP).det2 [j] [k] ;
             
            }
        }
    nmx = 5148;
    sqrth=TMath::Sqrt(3.)/2. ;
/*     initialise infocl */
    for (i__ = 1; i__ <= 2; ++i__) {
	for (j = 1; j <= 72; ++j) {
	    for (k = 1; k <= 72; ++k) {
		(*PPP).infocl[i__] [j] [k] = 0;
	    }
	}
    }
/*     compute number of pads with nonzero count ( nmx1 ) */
    i__1 = nmx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	id1 = (*PPP).iorder [1] [i__];
	id2 = (*PPP).iorder [2] [i__];
	if (d__[id1] [id2] > 0.0000000000001) {
	    (*PPP).nmx1 = i__;
          printf(" nmx1 -> %d %d %f %d \n",id1,id2, d__ [id1] [id2],(*PPP).nmx1 ) ;
	}
    }
    printf(" nmx1= %d \n",(*PPP).nmx1) ;
/*     crude clustering. start with pad with largest count ( cluster 1 ) 
*/
    icle = 1;
    idum = 0;
    id1 = (*PPP).iorder [1] [1];
    id2 = (*PPP).iorder [2] [1];
		(*PPP).infocl[1] [id1] [id2] = 1;
		(*PPP).infocl[2] [id1] [id2] = icle;
          printf(" infocl  %d %d %d %d %d \n",id1,id2,icle,(*PPP).infocl[1] [id1] [id2],(*PPP).infocl[2] [id1] [id2] ) ;
/*     check if neighbouring count in nonzero. If yes, include in the */
/*     cluster */
/*     infocl(1,i1,i2) =1 for primary ( peak ) =2 for secondary ( */
/*     neighbour ) */
/*     infocl(2,i1,i2) stores cluster no. */
    for (nb = 1; nb <= 6; ++nb) {
	jd1 = id1 + neibx[nb];
	jd2 = id2 + neiby[nb];
          printf(" infocl  %d %d %d %d %d %f %d %d\n",nb,jd1,id1,jd2,id2,d__[jd1] [jd2],icle,(*PPP).infocl[1] [jd1] [jd2] ) ;
	if (jd1 > 0 && jd1 <= 72 && jd2 > 0 && jd2 <= 72) {
	    if (d__[jd1] [jd2] > 0.0) {
		(*PPP).infocl[1] [jd1] [jd2] = 2;
		(*PPP).infocl[2] [jd1] [jd2] = icle;
	    }
	}
          printf(" infocl  %d %d %d %d %d %f %d %d\n",nb,jd1,id1,jd2,id2,d__[jd1] [jd2],icle,(*PPP).infocl[1] [jd1] [jd2] ) ;
    }
/*     start checking with next largest count and so on */
    i__1 = (*PPP).nmx1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	itr = 0;
	id1 = (*PPP).iorder[1] [i__];
	id2 = (*PPP).iorder[2] [i__];
       printf("iorder %d %d %d %d %d %d %d\n",i__,id1,id2,(*PPP).infocl [1] [id1] [id2],jd1,jd2,(*PPP).infocl [1] [jd1] [jd2] ) ; 
/*     nonzero infocl(1,id1,id2) means the pad is included in earlier 
*/
/*     cluster ( itr = 1 ) Otherwise itr = 0. */
	if ((*PPP).infocl [1] [id1] [id2] > 0) {
	    itr = 1;
	}
	if (itr ==  0) {
	    for (nb = 1; nb <= 6; ++nb) {
		jd1 = id1 + neibx[nb ];
		jd2 = id2 + neiby[nb ];
/*     check if a neighbouring pad belongs to previous cluster
. If yes, */
/*     icn = 1. Otherwise icn = 0. */
		icn = 0;
		if (jd1 > 0 && jd1 <= 72 && jd2 > 0 && jd2 <= 72) {
                    if ((*PPP).infocl [1] [jd1] [jd2] != 0) {
			icn = 1;
			icld = (*PPP).infocl [2] [jd1] [jd2] ;
		    }
		}
	    }
/*     itr = 0 and icn = 0 means new cluster. Store the cluster in
fo. in */
/*     infocl. */
	    if (TMath::Abs(icn) <= 0.1) {
		++icle;
		(*PPP).infocl[1] [id1] [id2] = 1;
		(*PPP).infocl[2] [id1] [id2] = icle;
              printf("ID1 icn %d %d %d %d %d\n",icn,icle,id1,id2,(*PPP).infocl[2] [id1] [id2]);
		for (nb = 1; nb <= 6; ++nb) {
/*     neighbours of the peak put with the peak. */
		    jd1 = id1 + neibx[nb ];
		    jd2 = id2 + neiby[nb ];
		    if (jd1 > 0 && jd1 <= 72 && jd2 > 0 && jd2 <= 72) {
	    if (d__[jd1] [jd2] != 0  && (*PPP).infocl [1] [jd1] [jd2] == 0) {
		(*PPP).infocl[1] [jd1] [jd2] = 2;
		(*PPP).infocl[2] [jd1] [jd2] = icle;
              printf("JD1 icn %d %d %d %d %d\n",icn,icle,jd1,jd2,(*PPP).infocl[2] [jd1] [jd2]);
			}
		    }
		}
	   } else {
/*     itr = 0 and icn = 1 means the pad belongs to cluster ic
ld */
		(*PPP).infocl[1] [id1] [id2] = 2;
		(*PPP).infocl[2] [id1] [id2] = icld;
		for (nb = 1; nb <= 6; ++nb) {
/*     and also the neighbours. */
		    jd1 = id1 + neibx[nb ];
		    jd2 = id2 + neiby[nb ];
		    if (jd1 > 0 && jd1 <= 72 && jd2 > 0 && jd2 <= 72) {
			if (d__[jd1] [jd2]  != 0. && (*PPP).infocl[1] [jd1] [jd2] == 0) {
			(*PPP).infocl[1] [jd1] [jd2] = 2;
			(*PPP).infocl[2] [jd1] [jd2] = icld;
			}
		    }
		}
         } else {
/*     itr = 0, so not a new cluster */
	    icld = (*PPP).infocl[2] [id1] [id2] ;
	    for (nb = 1; nb <= 6; ++nb) {
/*     neighbours put with the previous cluster. */
		jd1 = id1 + neibx[nb ];
		jd2 = id2 + neiby[nb ];
		if (jd1 > 0 && jd1 <= 72 && jd2 > 0 && jd2 <= 72) {
		if (d__[jd1] [jd2]  != 0. && (*PPP).infocl[1] [jd1] [jd2] == 0) {
			(*PPP).infocl[1] [jd1] [jd2] = 2;
			(*PPP).infocl[2] [jd1] [jd2] = icld;
		    }
		}
	    }
	}
    }
}
    printf(" ndet=%d,nmx1= %d , icle= %d icld= %d \n",ndet,(*PPP).nmx1,icle,icld) ;

/* L2: */
/* L3: */
    return ;
} /* crclust_ */

/* Subroutine */ void intpts_(PMD_dig_reco* PPP)
{
    Int_t i__1, i__2, i__3;

    Int_t ixpt;
    Float_t x, y, dx, dy;
    Int_t ix, iy;
    Float_t yl;
    Int_t idx, idy, ipt;
    Float_t wtx, wty;

    ixpt = 2;
    idx = -1;
    ipt = 0;
    dx = .25;
    i__1 = ixpt;
    for (ix = -ixpt; ix <= i__1; ++ix) {
	x = ix * dx;
	idx = -idx;
        i__2 = ix / ixpt ;
	wtx = (3. - idx -  (TMath::Abs(i__2) ))/ 3. * dx;
	yl = 1. / TMath::Sqrt(3.) * (1. - TMath::Abs(x));
	dy = yl * 2. / 4.;
	idy = -1;
	i__2 = ixpt;
	for (iy = -ixpt; iy <= i__2; ++iy) {
	    idy = -idy;
	    y = iy * dy;
            i__3 = iy / ixpt ;
	    wty = (3. - idy - (TMath::Abs(i__3))) / 3. * dy;
	    ++ipt;
	    (*PPP).pts[1] [ipt] = x;
	    (*PPP).pts[2] [ipt] = y;
	    (*PPP).wts[ipt] = wtx * wty;
	}
    }
     printf(" Exiting from intpts_ %f %f  \n", ipt,(*PPP).wts[ipt] );
    return ;
} /* intpts_ */

/* prints_  */ void prints_(Int_t ndet, PMD_dig_reco* PPP)
{
    Int_t i__, j;

    if (ndet == 1) {
       FILE *fp1;
       FILE *fp2;
       fp1=fopen("CPV_prime.out","w");
       fp2=fopen("CPV_secondary.out","w");
    }
    if (ndet == 2) {
       FILE *fp1;
       FILE *fp2;
       fp1=fopen("PMD_prime.out","w");
       fp2=fopen("PMD_secondary.out","w");
    }
    for (i__ = 1; i__ <= 72; ++i__) {
	for (j = 1; j <= 72; ++j) {
         if( (*PPP).infocl [1] [i__] [j] == 1) {
          fprintf(fp1," %f %f %d %d %d \n",(*PPP).coord [1] [i__] [j],(*PPP).coord [2] [i__] [j],i__,j,(*PPP).infocl [2] [i__] [j] ) ;
         }
         if( (*PPP).infocl [1] [i__] [j] == 2) {
          fprintf(fp2," %f %f %d %d %d \n",(*PPP).coord [1] [i__] [j],(*PPP).coord [2] [i__] [j],i__,j,(*PPP).infocl [2] [i__] [j] ) ;
         }
        }
    }
    return ;
} /* prints_ */

/* Subroutine */void  arrange_(PMD_dig_reco* PPP)
{
    Int_t i__1, i__2, i__3;

    Int_t ihld1, ihld2, ihld3, i__, j, k, itest;

    i__1 = (*PPP).incr;
    for (i__ = 2; i__ <= i__1; ++i__) {
	itest = 0;
	ihld1 = (*PPP).inford [1] [i__];
	ihld2 = (*PPP).inford [2] [i__];
	ihld3 = (*PPP).inford [3] [i__];
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    if (itest == 0 && (*PPP).inford [1] [j] > ihld1) {
		itest = 1;
		i__3 = j;
		for (k = i__ -1; k >= i__3; --k) {
                 (*PPP).inford [1] [k+1] = (*PPP).inford [1] [k] ;
                 (*PPP).inford [2] [k+1] = (*PPP).inford [2] [k] ;
                 (*PPP).inford [3] [k+1] = (*PPP).inford [3] [k] ;
		}
                 (*PPP).inford [1] [j] = ihld1 ;
                 (*PPP).inford [2] [j] = ihld2 ;
                 (*PPP).inford [3] [j] = ihld3 ;
	    }
	}
    }
     printf(" Exiting from arrange \n") ;
    return ;
} /* arrange_ */

/* refclust_ */void  refclust_(PMD_dig_reco* PPP)
{

    Int_t i__1, i__2;
    Float_t d__ [73] [73] ;
    Int_t npad[201], i__, j, k;
    Int_t id, id1, id2, iclb;

/*     using information from crude clustering, a refined algorithm for */
/*     clustering tried. At the moment using Gaussian fit with a number */
/*     of Gaussian functions. */
        for (j = 1; j <= 72; ++j) {
            for (k = 1; k <= 72; ++k) {
             d__ [j] [k] = (*PPP).det2 [j] [k] ;
            }
        }

    (*PPP).clno = 0;
/*     The No.of pads in each cluster = 200 */
    for (i__ = 1; i__ <= 200; ++i__) {
	npad[i__] = 0;
    }
    for (i__ = 1; i__ <= 5000; ++i__) {
	(*PPP).ncl[i__] = 0;
    }
/*     number of pads in each cluster determined. */
    i__1 = (*PPP).incr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	id = (*PPP).inford[1] [i__];
	++(*PPP).ncl[id];
        Int_t nclid = (*PPP).ncl[id] ;
    printf("id,ncl= %d %d\n",id,(*PPP).ncl[id]) ;
    }
    iclb = 0;
/*     number of clusters determined. */
    for (i__ = 1; i__ <= 5000; ++i__) {
	id = (*PPP).ncl[i__];
	if (id != 0 ) {
	    ++iclb;
	    ++npad[id];
	}
    }
    printf("iclb= %d \n",iclb) ;
    id = 0;
/*     define coordinates of each pad in a cluster and call gausspara */
/*     to determine Gaussian parameters. ( should it be eliminated ? ) */
    i__1 = iclb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	  i__2 = (*PPP).ncl[i__];
          printf("(*PPP).ncl(i)= %d %d\n",i__,(*PPP).ncl[i__]) ;
	for (j = 1; j <= i__2; ++j) {
	    ++id;
	    id1 = (*PPP).inford[2] [id];
	    id2 = (*PPP).inford[3] [id];
	    (*PPP).x[j] = (*PPP).coord[1] [id1] [id2];
	    (*PPP).y[j] = (*PPP).coord[2] [id1] [id2];
	    (*PPP).z__[j] = d__[id1] [id2];
	}
    printf(" I am in refclust calling ==gausspara_== \n") ;
	gausspara_(i__2,PPP );
    printf(" I am in refclust after ==gausspara_== \n") ;
    }
          printf("clno= %d \n",(*PPP).clno);
    return ;
} /* refclust_ */

/* gausspara_ */ void gausspara_(Int_t n,PMD_dig_reco* PPP)
{
    Int_t  i__1, i__2, i__3;
    Float_t d__1, d__2;

    Int_t   ihld, iord[501], idum, i__, j, k;
    Int_t   itest;
    Int_t   ig;
    Float_t pi, xi, yi, xj, yj, rrr;
    Float_t sum;

/*     single Gaussian of the form exp(-((x-x_0)^2+(y-y_0)^2)/r_0^2) */

    printf(" I am in ==gausspara_== %d  \n" , n) ;
    pi = 3.141593;
    Int_t nn = n ;
    if (nn == 1) {
	++(*PPP).clno;
/*     single pad cluster. Assumed fit to single Gaussian centered at the */
/*     center of the pad with Gaussian parameter r0 about 0.5 since */
/*     neighbouring pads don't fire. */
	(*PPP).x0[1] = (*PPP).x[1];
	(*PPP).y0[1] = (*PPP).y[1];
	(*PPP).r0[1] = 0.5;
/* Computing 2nd power */
	d__1 = (*PPP).r0[1];
	(*PPP).h0[1] = (*PPP).z__[1] / pi / (d__1 * d__1);
       printf("single pad ,clust no. = %d %d %f %f \n",n, (*PPP).clno, (*PPP).x0[1],(*PPP).y0[1] ) ;
    } else if (nn == 2) {
	++((*PPP).clno);
/*     two pad cluster. Assumed fit to single Gaussian centered at the */
/*     weighted average of the centers of two pads with r_0 about 0.5. */
/*     Could be a pair of clusters but no way of differentiating the two. */
/*         print *,'two pads   ','clust no.  ', clno */
       printf("two pad clusters  = %d \n",(*PPP).clno ) ;              
//       printf("z__(1),(2)  = %f %f \n",(*PPP).z__[1],(*PPP).z__[2] ) ; 
       Int_t btest = ((*PPP).z__[1])+((*PPP).z__[2]) ;
        if ( btest != 0 ) {
	(*PPP).x0[1] = ((*PPP).z__[1] * (*PPP).x[1] + (*PPP).z__[2] * (*PPP).x[2]) / ((*PPP).z__[1] + (*PPP).z__[2]);
	(*PPP).y0[1] = ((*PPP).z__[1] * (*PPP).y[1] + (*PPP).z__[2] * (*PPP).y[2]) / ((*PPP).z__[1] + (*PPP).z__[2]);
	(*PPP).r0[1] = 0.5;
        }
/* Computing 2nd power */
	d__1 = (*PPP).r0[1] = 0.5;
	(*PPP).h0[1] = ((*PPP).z__[1] + (*PPP).z__[2]) / pi / (d__1 * d__1);
//       printf("Two pad clusters  = %d %d %f %f \n",n, (*PPP).clno, (*PPP).x0[1],(*PPP).y0[1] ) ;
    } else {
/*     guessing possible no. of Gaussians required. */
/*     assuming that if pads are separated by more than one unit, they */
/*     are possible candidates for Gaussian centers. */
/*     do this from larger values of energy ( z here ). */
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iord[i__] = i__;
	}
	i__1 = nn;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    itest = 0;
	    ihld = iord[i__];
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		if (itest == 0 && (*PPP).z__[iord[j]] < (*PPP).z__[ihld]) {
		    itest = 1;
		    i__3 = j;
		    for (k = i__ - 1; k >= i__3; --k) {
			iord[k+1] = iord[k];
		    }
		    iord[j] = ihld;
		}
	    }
	}
/*     Gaussian centers determined. */
	ig = 1;
	(*PPP).icl[ig] = iord[ig];
	i__1 = nn;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    xi = (*PPP).x[iord[i__]];
	    yi = (*PPP).y[iord[i__]];
	    idum = 0;
	    i__2 = ig;
	    for (j = 1; j <= i__2; ++j) {
		xj = (*PPP).x[(*PPP).icl[j]];
		yj = (*PPP).y[(*PPP).icl[j]];
/* Computing 2nd power */
		d__1 = xi - xj;
/* Computing 2nd power */
		d__2 = yi - yj;
		rrr = TMath::Sqrt(d__1 * d__1 + d__2 * d__2);
		if (1.01 < rrr) {
		    ++idum;
		}
	    }
	    if (idum == ig) {
		++ig;
		(*PPP).icl[ig] = iord[i__];
	    }
	}
	(*PPP).clno += ig;
       printf("No. of pads = %d , No. of Gaussians = %f \n ",ig,(*PPP).clno ) ;
	sum = 0.;
	i__1 = ig;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    (*PPP).x0[i__] = (*PPP).x[(*PPP).icl[i__]];
	    (*PPP).y0[i__] = (*PPP).y[(*PPP).icl[i__]];
	    (*PPP).r0[i__] = .5;
	    sum += (*PPP).z__[(*PPP).icl[i__]];
	}
	i__1 = ig;
	for (i__ = 1; i__ <= i__1; ++i__) {
        if ( sum != 0 ) {
	    (*PPP).h0[i__] = (*PPP).z__[(*PPP).icl[i__]] / sum;
        }
	}
    printf(" I am in ==gausspara_== Calling gausfit_ \n") ;
	gausfit_(ig, nn, PPP );
    printf(" I am in ==gausspara_ after gausfit_ \n") ;
    }
    return ;
} /* gausspara_ */

/*  gausfit_ */ void gausfit_(Int_t ng, Int_t nn, PMD_dig_reco* PPP )
{
    Int_t  i__1;
    Float_t d__1;

    Int_t   iran;
    Float_t aint;
    Int_t   i__;
    Float_t s;
    Float_t  pi, ss, st;

    pi = 3.141593;
    printf(" ng=%d , nn=%d \n",ng,nn);
    s = 0.;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s +=(*PPP).z__[i__];
    }
/*      print *,s */
     printf(" s=%d  \n",s);
    i__1 = ng;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
//       printf("i__, X0,y0,r0 = %d %f %f %f \n",i__, (*PPP).x0[i__],(*PPP).y0[i__],(*PPP).r0[i__]) ;
	d__1 = (*PPP).r0[i__];
	(*PPP).h0[i__] = (*PPP).h0[i__] * s / pi / (d__1 * d__1);
    }
    st = 0.;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
        Int_t pbab = 1 ;
	aintg_(ng, i__, pbab, aint, PPP);
/* Computing 2nd power */
	d__1 = aint - (*PPP).z__[i__];
	if( s > 0. ) st += d__1 * d__1 / s;
    }
//   printf(" st=%d ",st);
    for (iran = 1; iran <= 5000; ++iran) {
        Float_t rand;
	i__1 = ng;
	for (i__ = 1; i__ <= i__1; ++i__) {
            ranmar_(rand ,PPP) ;
	    (*PPP).xx[i__] = (*PPP).x0[i__] + (rand - .5) * 1.;
            ranmar_(rand ,PPP) ;
	    (*PPP).yy[i__] = (*PPP).y0[i__] + (rand - .5) * 1.;
            ranmar_(rand ,PPP) ;
	    (*PPP).rr[i__] = (*PPP).r0[i__] * (((rand -.5) * .2 )+ 1.);
/* Computing 2nd power */
	    d__1 = (*PPP).r0[i__] / (*PPP).rr[i__];
	    (*PPP).hh[i__] = (*PPP).h0[i__] * (d__1 * d__1);
	}
	ss = 0.;
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
            Int_t pbab = 2 ;
	    aintg_(ng, i__, pbab, aint, PPP);
/* Computing 2nd power */
	    d__1 = aint - (*PPP).z__[i__];
	    if(s>0.0)ss += d__1 * d__1 / s;
	}
	if (ss < st) {
	    i__1 = ng;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		(*PPP).x0[i__] = (*PPP).xx[i__];
		(*PPP).y0[i__] = (*PPP).yy[i__];
		(*PPP).r0[i__] = (*PPP).rr[i__];
		(*PPP).h0[i__] = (*PPP).hh[i__];
	    }
	    st = ss;
    printf(" st=%d \n",st);
	}
    }
/*      print 1,((*PPP).x0(i),(*PPP).(*PPP).y0(i),(*PPP).r0(i),i=1,ng) */
/*      print *,st */
            i__1 = ng;
        for (i__ = 1; i__ <= i__1; ++i__) {
//             printf("x0(i),y0(i),r0(i) =%f %f %f \n",(*PPP).x0[i__1],(*PPP).y0[i__1],(*PPP).r0[i__1] ) ;
        }
    printf(" I am ending gausfit_ \n") ;
    return ;
} /* gausfit_ */

/* aintg_ */void aintg_(Int_t ng, Int_t n, Int_t bab, Float_t aint, PMD_dig_reco* PPP )
{
    Int_t i__1;
    Float_t  d__1, d__2, d__3;

    Int_t i__, j;
    Float_t xp, yp, fun;

    aint = 0.;
    for (i__ = 1; i__ <= 25; ++i__) {
	xp = (*PPP).pts[1] [i__] + (*PPP).x[n] ;
	yp = (*PPP).pts[2] [i__] + (*PPP).y[n] ;
// printf("x[n],y[n],pts[1] [i__],[2] [i__], xp,yp =%f %f %f %f %f %f\n",(*PPP).x[n],(*PPP).y[n],(*PPP).pts[1] [i__],(*PPP).pts[2] [i__],xp,yp) ;
	fun = 0.;
	i__1 = ng;
	for (j = 1; j <= i__1; ++j) {
//printf("J,D__1,D__2,D__3= %d %f %f %f %f\n",j,(*PPP).xx[j],xp,(*PPP).yy[j],(*PPP).rr[j]) ;
/* Computing 2nd power */
	    d__1 = xp - (*PPP).xx[j];
/* Computing 2nd power */
	    d__2 = yp - (*PPP).yy[j];
/* Computing 2nd power */
	    if( bab == 1) d__3 = (*PPP).r0[j];
	    if( bab == 2) d__3 = (*PPP).rr[j];
// printf("J,D__1,D__2,D__3= %d %f %f %f \n",j,d__1,d__2,d__3) ;
	    fun += TMath::exp(-(d__1 * d__1 + d__2 * d__2) / (d__3 * d__3)) * (*PPP).hh[j];
	}
	aint += fun * (*PPP).wts[i__];
//  printf("FUN,aint= %f %f \n",fun,aint) ;
    }
    return ;
} /* aintg_ */
/* ranmar() >*/void ranmar_(Float_t random , PMD_dig_reco* PPP)
{
    /* Initialized data */

    Int_t i__ = 97;
    Int_t j   = 33;
    Float_t    uni;

/* ===========================================C==========================C
 */
/* Universal random number generator proposed by Marsaglia and Zaman */
/* in report FSU-SCRI-87-50 */
    uni = (*PPP).u[i__ - 1] - (*PPP).u[j - 1];
    if (uni < 0.) {
        uni += 1.;
    }
    (*PPP).u[i__ ] = uni;
    --i__;
    if (i__ == 0) { i__ = 97; }
    --j;
    if (j   == 0) { j   = 97; }
    (*PPP).c__ -= (*PPP).cd;
    if ((*PPP).c__ < 0.) { (*PPP).c__ += (*PPP).cm; }
    uni -= (*PPP).c__;
    if (uni <= 0.) { uni += 1.; }
    random = uni;
} /* ranmar_ */

