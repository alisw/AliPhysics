/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 *                                                                        *
 * /////////////////////////////////////////////////////////////////////  *
 *                                                                        *
 * This class performs a fast fit of helices going through the <=6        *
 * points of the ITS, with the goal of studying tracking and              *
 * vertexing performances.                                                *
 * Generated kinematics is used to take into account different weights    *
 * associated to points in different layers (with different multiple      *
 * scattering-originated errors).                                         *
 *                                                                        *
 *   Based on the work by A. Strandlie, R. Fruhwirth                      *
 *                                                                        *
 *   First implementation by N. Bustreo, R. Turrisi - July 2000           *
 *                                                                        *
 *   Further modifications by A. Dainese, R. Turrisi                      *
 *                                                                        *
 *   Contact: Rosario Turrisi, rosario.turrisi@pd.infn.it                 *
 *                                                                        *
 * **************************************************************************/
//
//
//       Modified November, 7th 2001 by Rosario Turrisi 
//       (rosario.turrisi@pd.infn.it)
//
//       FitHelix returns different values. 0=ok, >0 =problem
//       void FitLinear -> Int_t FitLinear to give feedback of errors to FitHelix
//
//
//       Modified July, 30th 2001 by Rosario Turrisi 
//       (rosario.turrisi@pd.infn.it)
//       
//        Fit for z now in (z,s) plane.
//        Returns parameters in order to write the helix equation
//        and find the right phase/initial point.
//
//     "PROPER WEIGHTS":  (1+R^2)^2/(\sigma_x^2 + \sigma_y^2 + \sigma_MS^2)
//
#include <Riostream.h>
#include "AliITSRiemannFit.h"
#include "AliRun.h"
#include "TClonesArray.h"
#include "stdio.h"
#include "stdlib.h"
#include "Riostream.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TParticle.h"
#include "TFile.h"
#include "AliITSRecPoint.h"
#include "AliITSgeom.h"
#include "AliITSmodule.h"

ClassImp(AliITSRiemannFit)


AliITSRiemannFit::AliITSRiemannFit() {
  ///////////////////////////////////////////////////////////
  // Default constructor.
  // Set everything to zero.
  ////////////////////////////////////////////////////////////

  fSizeEvent     = 0;
  fPoints        = 0;
  fPrimaryTracks = 0;
  fPointRecs     = 0;
  //
  //  test erase
//    fspdi          = 0;
//    fspdo          = 0;
  for(Int_t i=0;i<6;i++)fPLay[i] = 0;
  
}
//----------------------------------------------------------------------

AliITSRiemannFit::~AliITSRiemannFit() {
  ///////////////////////////////////////////////////////////
  // Default destructor.
  // if arrays exist delete them. Then set everything to zero.
  ////////////////////////////////////////////////////////////
   if(fPointRecs!=0){
      for(Int_t i=0;i<fSizeEvent;i++) delete[] fPointRecs[i];
      delete[] fPointRecs;
   } // end if fPointRecs!=0
  fSizeEvent     = 0;  
  fPointRecs     = 0;
  fPoints        = 0;
  fPrimaryTracks = 0;
  //
  // test erase
//    fspdi          = 0;
//    fspdo          = 0;
  for(Int_t i=0;i<6;i++)fPLay[i] = 0;
  return;
}
//----------------------------------------------------------------------

AliITSRiemannFit::AliITSRiemannFit(Int_t size,Int_t ntracks) {
  ///////////////////////////////////////////////////////////
  // Constructor.
  // Set fSizeEvent to size and fPrimaryTracks to ntracks.
  // Others to zero.
  ////////////////////////////////////////////////////////////

  fSizeEvent     = size;
  fPoints        = 0;
  fPrimaryTracks = ntracks;
  //
  // test erase
//    fspdi          = 0;
//    fspdo          = 0;
  Point_tl *first = new Point_tl[fSizeEvent];
  Point_tl **PointRecs = new Point_tl*[fSizeEvent];
  for(Int_t i=0;i<6;i++)fPLay[i] = 0;
  for(Int_t j=0;j<fSizeEvent;j++)  // create an array of struct
    PointRecs[j] = &(first[j]);   
}
// ---------------------------------------------------------------------

void FillPoints(Point_tl **Points,Int_t &index,Float_t *xpoint,
		Float_t *error,
		TLorentzVector PE,TLorentzVector OT,Int_t *id,
		Int_t track,const Char_t *name,Int_t code,
		Float_t phiorigin){
  ///////////////////////////////////////////////////////////////////////
  // Fill the structure Point_tl with the proper data
  //
  //////////////////////////////////////////////////////////////////////
  Float_t PI2 = 2.0*TMath::Pi();
  Float_t phi,r,x,y,z;
  Int_t i;
  i = index;
  x = xpoint[0];
  y = xpoint[1];
  z = xpoint[2];
  r = sqrt(x*x+y*y);
  phi = TMath::ATan2(y,x);
  if(phi<0.0) phi += PI2;
  Points[i]->phi = phi;
  Points[i]->eta  = -0.5*tan(0.5*TMath::ATan2(r,z));
  Points[i]->fx = x;
  Points[i]->fy = y;
  Points[i]->fz = z;
  Points[i]->fdx = error[0];
  Points[i]->fdy = error[1];
  Points[i]->fdz = error[2];
  Points[i]->fr = r;
  Points[i]->track  = track;
  Points[i]->lay    = id[0],
  Points[i]->lad    = id[1];
  Points[i]->det    = id[2];
  Points[i]->fMomentum = PE;
  Points[i]->fOrigin   = OT;
  Points[i]->fPt   = sqrt(PE.X()*PE.X()+PE.Y()*PE.Y());
  Points[i]->fCode = code;
  Points[i]->fName = name;
  Points[i]->vertexPhi = phiorigin;
  index++;
  return;
  
}
// -----------------------------------------------------------------------

void AliITSRiemannFit::InitPoints(Int_t ntracks,AliITS *ITS,
				    TTree *TR,Int_t nparticles){
  //////////////////////////////////////////////////////////////////////
  // Fill the class member fPointRecs with the reconstructed points
  // Set All other members to the real values
  //
  /////////////////////////////////////////////////////////////////////
  printf("\n ************* Starting Init Points *************\n");
  TParticle *part;
  AliITSgeom *gm = (AliITSgeom*)ITS->GetITSgeom();
  //get pointer to modules array
  TObjArray *ITSmodules = ITS->GetModules();
  Int_t nmodules=ITSmodules->GetEntriesFast();
  printf("nmodules = %d \n",nmodules);
  // Get the points from points file
  AliITSmodule *itsModule;
  Int_t mod,irec;
  Stat_t nent;
  AliITSRecPoint  *recp;
  nent=TR->GetEntries();
  TClonesArray *ITSrec  = ITS->RecPoints();

  Int_t TotRP=0;
  for (mod=0; mod<nmodules; mod++) {
    itsModule=(AliITSmodule*)ITSmodules->At(mod);
    ITS->ResetRecPoints();
    TR->GetEvent(mod);
    Int_t nrecp = ITSrec->GetEntries();
    if(!nrecp) continue;
    TotRP += nrecp;
  }

  Int_t iMAX = TotRP;
  fPrimaryTracks = ntracks;
  fParticles     = nparticles;
  Point_tl *global = new Point_tl[iMAX];
  fPointRecs = new Point_tl*[iMAX];
  //
  // test erase
//    Point_tl *first = new Point_tl[iMAX];
//    Point_tl *second = new Point_tl[iMAX];
//    fspdi = new Point_tl*[iMAX];
//    fspdo = new Point_tl*[iMAX];
  for(Int_t j=0;j<iMAX;j++) {
    fPointRecs[j] = &(global[j]);
    //
    // test erase
//      fspdi[j]      = &(first[j]);
//      fspdo[j]      = &(second[j]);
  }
  
  Int_t ieta=0,ieta2=0;
  Int_t i,id[4],idold[4];
  Int_t track=0;//         // track  of hit
  Float_t xpoint[3],error_plus[3],error_minus[3],global_error[3];       // position and error of the point
  TLorentzVector OT,PE;
  Float_t locals[3],locals_error[3],locals_plus[3],locals_minus[3]; // local position and local errors
  Float_t Phi;
  Int_t code;
  const char *name;
  Int_t layer,ladder,detector;
  Float_t xcluster,zcluster;
  Int_t num=0,nspdi=0,nspdo=0,nsddi=0,nsddo=0,nssdi=0,nssdo=0;
 
  for (mod=0; mod<nmodules; mod++) {
    itsModule=(AliITSmodule*)ITSmodules->At(mod);
    ITS->ResetRecPoints();
    TR->GetEvent(mod);
    Int_t nrecp = ITSrec->GetEntries();
    if (!nrecp) continue;
    itsModule->GetID(layer,ladder,detector);

    for (irec=0;irec<nrecp;irec++) {
      recp   = (AliITSRecPoint*)ITSrec->UncheckedAt(irec);
      track=recp->fTracks[0];
      if(track <0 ) continue;
      xcluster=recp->GetX();     // x on cluster
      zcluster=recp->GetZ();     // z on cluster
      part   = (TParticle*) gAlice->Particle(track);    
      part->ProductionVertex(OT);  // set the vertex 
      part->Momentum(PE);          // set the vertex momentum
      name      = part->GetName();
      code      = part->GetPdgCode();
      Phi       = part->Phi();
      id[0]=layer;
      id[1]=ladder;
      id[2]=detector;
      id[3]=irec;
      locals[0]=xcluster;     // x on cluster
      locals[1]=0.0;          // y on cluster
      locals[2]=zcluster;     // z on cluster
      locals_error[0]=sqrt(recp->GetSigmaX2());
      locals_error[1]=0.0;
      locals_error[2]=sqrt(recp->GetSigmaZ2());
      locals_plus[0]=xcluster+sqrt(recp->GetSigmaX2());       // x on cluster
      if(layer==1||layer==2) locals_plus[1]=0.0150/2;         // y on cluster
      else if(layer==3||layer==4) locals_plus[1]=0.0280/2;    // y on cluster
      else if(layer==5||layer==6) locals_plus[1]=0.0300/2;    // y on cluster
      locals_plus[2]=zcluster+sqrt(recp->GetSigmaZ2());       // z on cluster
      locals_minus[0]=xcluster-sqrt(recp->GetSigmaX2());      // x on cluster
      locals_minus[1]=0.0;                                    // y on cluster
      locals_minus[2]=zcluster-sqrt(recp->GetSigmaZ2());      // z on cluster

      gm->LtoG(layer,ladder,detector,locals,xpoint);
      gm->LtoG(layer,ladder,detector,locals_plus,error_plus);
      gm->LtoG(layer,ladder,detector,locals_minus,error_minus);
      global_error[0]=0.5*TMath::Abs(error_plus[0]-error_minus[0]);
      global_error[1]=0.5*TMath::Abs(error_plus[1]-error_minus[1]);
      global_error[2]=0.5*TMath::Abs(error_plus[2]-error_minus[2]);
      if(track<ntracks) {
	if(TMath::Abs(part->Eta())<=1.0) ieta++;
	if(TMath::Abs(part->Eta())<=0.5) ieta2++;
      }
      if(!(id[0]==idold[0]&&id[1]==idold[1]&&
	   id[2]==idold[2]&&id[3]==idold[3])) {
	FillPoints(fPointRecs,num,xpoint,global_error,PE,OT,id,track,name,code,Phi);
	//
	// test erase	
	switch (idold[0]) {
	case 1:
	  nspdi++;
	  break;
	case 2:
	  nspdo++;
	  break;
	case 3:
  	  nsddi++;
	  break;
	case 4:
  	  nsddo++;
	  break;
	case 5:
  	  nssdi++;
	  break;
	case 6:
  	  nssdo++;
	  break;
	}
//  	if(idold[0]==1){
//  	  FillPoints(fspdi,nspdi,xpoint,global_error,PE,OT,id,track,name,code,Phi);  
//  	}
//  	if(idold[0]==2){
	 
//  	  FillPoints(fspdo,nspdo,xpoint,global_error,PE,OT,id,track,name,code,Phi);
//  	}
//  	if(idold[0]==3){
//  	  nsddi++;
//  	}
//  	if(idold[0]==4){
//  	  nsddo++;
//  	}
//  	if(idold[0]==5){
//  	  nssdi++;
//  	}
//  	if(idold[0]==6){
//  	  nssdo++;
//  	}
	for(i=0;i<4;i++) idold[i] = id[i];
	for(i=0;i<3;i++) xpoint[i]    = 0.0;
      } // end if id != idold
    } // end for irec
  }// end for mod

  fPoints = num;
  fSizeEvent = num;
  fPLay[0] = nspdi ;
  fPLay[1] = nspdo ;
  fPLay[2] = nsddi ;
  fPLay[3] = nsddo ;
  fPLay[4] = nssdi ;
  fPLay[5] = nssdo ;
  printf("%d primary tracks in eta=+-1\n",ieta);
  printf("%d primary tracks#2 in eta=+-0.5\n",ieta2);
  printf("\nInitPoints :\n\nPoints on Layer1 : %d on Layer2 : %d\n",nspdi,nspdo);
  printf("Points on Layer3 : %d on Layer4 : %d\n",nsddi,nsddo);
  printf("Points on Layer5 : %d on Layer6 : %d\n",nssdi,nssdo);
  printf("Points on all Layers: %d\n",num);
  printf("\n ************* Init Points Finished *************\n");
  return;
}
// ------------------------------------------------------------------------
///////////////////////////////////////////////////////////
// Functions for sorting the fPointRecs array
///////////////////////////////////////////////////////////
Bool_t SortZ(const Point_tl *s1,const Point_tl *s2){
  // Z sorting function for qsort.
   Float_t a;

   a = s1->fz - s2->fz;
   if(a<0.0) return kTRUE;
   if(a>0.0) return kFALSE;
   return kFALSE;
}
Bool_t SortTrack(const Point_tl *s1,const Point_tl *s2){
  // track sorting function for qsort.
   Float_t a;

   a = s1->track - s2->track;
   if(a<0.0) return kTRUE;
   if(a>0.0) return kFALSE;
   return kFALSE;
}
void hpsortTrack(Point_tl **ra,Int_t n){
   Int_t i,ir,j,l;
   Point_tl *rra;

   if(n<2) return;

   l  = ((n-1) >> 1) +1; // divide 2 + 1
   ir = n-1;
   for(;;){
     if(l>0){
        rra = ra[--l];  // decrement first
     }else{
        rra    = ra[ir];
        ra[ir] = ra[0];
        if(--ir == 0){  // decrement first
           ra[0] = rra;
           break;
        } // if --ra == 0 
     } // end l>0 
     i = l;
     j = l+1;
     while(j<=ir){
        if( j<ir && SortTrack(ra[j],ra[j+1]) ) j++;
        if( SortTrack(rra,ra[j]) ){
           ra[i] = ra[j];
           i = j;
           j <<= 1; // time 2.
        }else{
           break;
        } // end if func() 
     } // end while 
     ra[i] = rra;
   } // end for ever 
}
void hpsortZ(Point_tl **ra,Int_t n){
   Int_t i,ir,j,l;
   Point_tl *rra;

   if(n<2) return;

   l  = ((n-1) >> 1) +1; // devide 2 + 1
   ir = n-1;
   for(;;){
     if(l>0){
        rra = ra[--l];  // decrament first
     }else{
        rra    = ra[ir];
        ra[ir] = ra[0];
        if(--ir == 0){  // decrament first
           ra[0] = rra;
           break;
        } // if --ra == 0 
     } // end l>0 
     i = l;
     j = l+1;
     while(j<=ir){
        if( j<ir && SortZ(ra[j],ra[j+1]) ) j++;
        if( SortZ(rra,ra[j]) ){
           ra[i] = ra[j];
           i = j;
           j <<= 1; // time 2.
        }else{
           break;
        } // end if func() 
     } // end while 
     ra[i] = rra;
   } // end for ever 
}
//-----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////
//      Sorting functions
///////////////////////////////////////////////////////////////////
Int_t Partition(Int_t array[],Int_t left,Int_t right){
  Int_t val = array[left];
  Int_t lm = left - 1;
  Int_t rm = right + 1;
  for(;;) {
    do 
      rm--;
    while
      (array[rm]>val);
    do 
      lm++;
    while
      (array[lm]<val);
    if(lm<rm){
      Int_t tempr = array[rm];
      array[rm]=array[lm];
      array[lm]=tempr;
    }
    else
      return rm;
  }

  return 1;
}

///////////////////////////////////////////////////////////////////////

void AliITSRiemannFit::WritePoints(void) {
  /////////////////////////////////////////////////////////////////////
  // write the data in a file (temporary ascii)
  /////////////////////////////////////////////////////////////////////
  FILE *ascii= fopen("AsciiPoints.dat","w");
  for(Int_t i=0;i<fPoints;i++) {
    fprintf(ascii,"%d\t%d\t%f\t%f\t%f\n",fPointRecs[i]->lay,
	    fPointRecs[i]->track,fPointRecs[i]->fx,fPointRecs[i]->fy,
	    fPointRecs[i]->fz);
  }
  fclose(ascii);
  return;
}
//-----------------------------------------------------------------------

void AliITSRiemannFit::ReadPoints(void) {
  //////////////////////////////////////////////////////////////////////
  // read the filled array
  /////////////////////////////////////////////////////////////////////
  hpsortTrack(fPointRecs,fPoints);
  for(Int_t i=0;i<fPoints;i++) 
    printf("%d\t%d\t%d\t%f\t%f\t%f\t(%.0f,%.0f,%.0f)\t%.3f\t%s\n",
	   i,fPointRecs[i]->lay,fPointRecs[i]->track,fPointRecs[i]->fx,
	   fPointRecs[i]->fy,fPointRecs[i]->fz,fPointRecs[i]->fOrigin.X(),
	   fPointRecs[i]->fOrigin.Y(),fPointRecs[i]->fOrigin.Z(),
	   fPointRecs[i]->fPt,fPointRecs[i]->fName);
  return; 
}
//-----------------------------------------------------------------------

Int_t AliITSRiemannFit::SolveCubic(Double_t a,Double_t b,Double_t c,
				    Double_t &x1,Double_t &x2,Double_t &x3){
   //////////////////////////////////////////////
   ///  Solve cubic equation:
   ///  x^3 + a*x^2 +b*x + c
   ///  
   ///  returns  x1 , x2 , x3
   ////////////////////////////////////////
   
   Double_t Q = ((a*a - 3*b)/9);
   Double_t R = ((2*a*a*a - 9*a*b +27*c)/54);
   Double_t theta;
   Double_t F = -2*sqrt(Q);
   Double_t g = a/3;
   Double_t PI2 = TMath::Pi()*2;

  if( R*R>Q*Q*Q ) {
    cout<<"\nTrack "<<"Determinant :\n\t\t No Real Solutions !!!\n"<<endl;
    x1 = 9999999;
    x2 = 9999999;
    x3 = 9999999;
    return 0;
  }

  theta = TMath::ACos( R/sqrt(Q*Q*Q));

  x1 = (F*TMath::Cos(theta/3))-g;
  x2 = (F*TMath::Cos((theta+PI2)/3))-g;
  x3 = (F*TMath::Cos((theta-PI2)/3))-g;
  
  return 1;
}
//-----------------------------------------------------------------

void RiemannTransf(Int_t npoints,TVector3 **From,TVector3 **To) {
  ///////////////////////////////////////////////////////////////////////
  //   This function apllies the transformation in the Riemann sphere
  //   for xy plane
  ///////////////////////////////////////////////////////////////////////
  Float_t *R = new Float_t[npoints];
  Float_t *Theta = new Float_t[npoints];
  Float_t PI2 = 2*TMath::Pi();
  Float_t x=0,y=0,z=0;
  
  for(Int_t i=0;i<npoints;i++) {
    R[i]     = sqrt(From[i]->X()*From[i]->X()+From[i]->Y()*From[i]->Y());
    Theta[i] = TMath::ATan2(From[i]->Y(),From[i]->X());
    if(Theta[i]<0) Theta[i]+=PI2;
    x     = R[i]*cos(Theta[i])/(1+R[i]*R[i]);
    y     = R[i]*sin(Theta[i])/(1+R[i]*R[i]);
    z     = R[i]*R[i]/(1+R[i]*R[i]);
    To[i]->SetXYZ(x,y,z);
  }
  delete[] R;
  delete[] Theta;
  return;
}


//---------------------------------------------------------------------

Int_t FitLinear(Int_t npoints, TVector3 **input, TVector3 **errors, Double_t omega,
		Double_t &thu0, Double_t &thv0, Double_t &phi, TVector2 &zData, TVector3 &zError, 
		Double_t &CorrLin){
  ///////////////////////////////////////////////////////////////////////
  //   Fit the points in the (z,s) plane - helix 3rd equation
  //
  /////////////////////////////////////////////////////////////////////// 
  Int_t direction=0;
  //PH  Double_t z[npoints],x[npoints],y[npoints],s[npoints];
  //PH  Double_t ez[npoints],ex[npoints],ey[npoints],es[npoints];
  Double_t * z = new Double_t[npoints];
  Double_t * x = new Double_t[npoints];
  Double_t * y = new Double_t[npoints];
  Double_t * s = new Double_t[npoints];
  Double_t * ez = new Double_t[npoints];
  Double_t * ex = new Double_t[npoints];
  Double_t * ey = new Double_t[npoints];
  Double_t * es = new Double_t[npoints];
  Double_t z0=0.0,vpar=0.0,ez0=0.0,evpar=0.0, chisquare;

  //  Double_t chi=TMath::Pi()/2.0+phi;
  Double_t chi=-TMath::Pi()-phi;
  Double_t angold=0.0, tpang=0.0;
  for(Int_t k = 0; k<npoints; k++) {
    x[k] = 10.0*input[k]->X();   ex[k] = 10.0*errors[k]->X();
    y[k] = 10.0*input[k]->Y();   ey[k] = 10.0*errors[k]->Y();
    z[k] = 10.0*input[k]->Z();   ez[k] = 10.0*errors[k]->Z();
    if(TMath::Abs(x[k]-thu0)<1.0e-5) {  // should never happen, nor give troubles...
      chisquare=9999.99; 
      cerr<<"limit for  x-x_0 "<<x[k]<<" "<<thu0<<endl; 
      delete [] z;
      delete [] x;
      delete [] y;
      delete [] s;
      delete [] ez;
      delete [] ex;
      delete [] ey;
      delete [] es;
      return 12;
    }
    Double_t ang1=TMath::ATan2((y[k]-thv0),(x[k]-thu0));
    if( (x[k]-thu0)<0 ) {
      if (ang1*angold<0) {
	tpang=ang1-TMath::Sign(TMath::Pi()*2.0,ang1);
	ang1=tpang;
      }
    }
    angold=ang1;
    if (k>0) direction+=(z[k]>z[k-1] ? 1 : -1);
    s[k] = (ang1+chi)/omega;
    es[k]=TMath::Sqrt(ey[k]*ey[k]+ex[k]*ex[k]/TMath::Power((x[k]-thu0),4))*TMath::Abs(s[k]);
  }
  if ( TMath::Abs(direction) != (npoints-1) ) {return 11;} 

  TGraphErrors *fitHist = new TGraphErrors(npoints,s,z,es,ez);
  fitHist->Fit("pol1","Q");
  z0  = fitHist->GetFunction("pol1")->GetParameter(0);
  vpar  = fitHist->GetFunction("pol1")->GetParameter(1);
  ez0 = fitHist->GetFunction("pol1")->GetParError(0);
  evpar = fitHist->GetFunction("pol1")->GetParError(1);
  chisquare = fitHist->GetFunction("pol1")->GetChisquare();
  zData.Set(z0,vpar);
  zError.SetXYZ(ez0,evpar,chisquare);

  Double_t Sigmas=0.; 
  Double_t Sigmaz=0.;
  Double_t Avs=0.;
  Double_t Avz=0.;
  Double_t Avsz=0.;

  for(Int_t j = 0; j < npoints; j++) {
    Avs  += s[j];
    Avz  += z[j]; 
    Avsz += s[j]*z[j];
  }
    Avs  /= (Double_t)npoints;
    Avz  /= (Double_t)npoints; 
    Avsz /= (Double_t)npoints;

  for(Int_t l = 0; l < npoints; l++) {
    Sigmas += (s[l]-Avs)*(s[l]-Avs);
    Sigmaz += (z[l]-Avz)*(z[l]-Avz);
  }
  Sigmas /=(Double_t)npoints;
  Sigmaz /=(Double_t)npoints;

  Sigmas = sqrt(Sigmas);
  Sigmaz = sqrt(Sigmaz);

  CorrLin = (Avsz-Avs*Avz)/(Sigmas*Sigmaz);

  delete [] z;
  delete [] x;
  delete [] y;
  delete [] s;
  delete [] ez;
  delete [] ex;
  delete [] ey;
  delete [] es;
  
  return 0;
}

//-------------------------------------------------------------------
Int_t AliITSRiemannFit::FitHelix(Int_t tracknumber,Double_t Px,Double_t Py,Double_t Pz,Double_t& fd0,
				   Double_t& fphi,Double_t& u0, Double_t& v0, Double_t& rho,Double_t& omega, Double_t& z0,
				   Double_t& vpar,Double_t& chisql, Double_t& fCorrLin,Double_t& fFit,
				   Int_t first,Int_t second,Int_t third,Int_t fourth,Int_t fifth,Int_t sixth) {
  ///////////////////////////////////////////////////////////////////////
  //  This function finds the helix paramenters 
  //  d0  = impact parameter
  //  rho = radius of circle
  //  phi = atan(y0/x0)
  //  for the xy plane
  //  starting from the momentum and the outcome of
  //  the fit on the Riemann sphere (i.e. u0,v0,rho)
  //
  //   MIND !!!!   Here we assume both angular velocities be 1.0  (yes, one-dot-zero !)
  //
  //
  ///////////////////////////////////////////////////////////////////////
  //
  //  All this stuff relies on this hypothesis !!!
  //
//   FILE *pout=fopen("chisql.dat","a");
  Int_t ierr = 0, ierrl=0;
  omega = 1.0e-2;

  Int_t bitlay[6]={1,1,1,1,1,1}; 
  bitlay[0]*=first; bitlay[1]*=second; bitlay[2]*=third; bitlay[3]*=fourth; bitlay[4]*=fifth; bitlay[5]*=sixth;
  fd0 = -9999;   // No phisycs value
  u0 = -9999.9999; // parameters of helix - strange value...
  v0 = -9999.9999; // parameters of helix - strange value...
  rho = -9999.9999; // parameters of helix -unphysical strange value...
  Int_t Layer = 0;
  const Char_t* name = 0;
  Int_t i=0,k=0;
  Int_t iMAX = 50;
  Int_t N = 0;
  Int_t npl[6]={0,0,0,0,0,0};
  Double_t P = sqrt(Px*Px+Py*Py+Pz*Pz);
  Double_t Pt = sqrt(Px*Px+Py*Py);
  TVector3 zError;
  TVector2 zData;
  Double_t CorrLin;
  TVector3 *ori        = new TVector3[iMAX];
  TVector3 **original  = new TVector3*[iMAX];
  TVector3 *rie       = new TVector3[iMAX];
  TVector3 **Riemann  = new TVector3*[iMAX];
  TVector3 *err        = new TVector3[iMAX];
  TVector3 **errors    = new TVector3*[iMAX];
  TVector3 *linerr        = new TVector3[iMAX];
  TVector3 **linerrors    = new TVector3*[iMAX];
  //PH  Double_t Weight[iMAX];
  Double_t * Weight = new Double_t[iMAX];

  for(i=0;i<iMAX;i++){
    original[i]   = &(ori[i]);
    Riemann[i]   = &(rie[i]);
    errors[i]     = &(err[i]);
    linerrors[i]  = &(linerr[i]);
  }
  for(k =0;k<iMAX;k++) original[k]->SetXYZ(9999,9999,9999);
  Double_t A11,A12,A13,A21,A22,A23,A31,A32,A33;
  A11=0;A12=0;A13=0;A21=0;A22=0;A23=0;A31=0;A32=0;A33=0;
  Double_t xbar = 0;
  Double_t ybar = 0;
  Double_t zbar = 0;
  Double_t a,b,c,d;       // cubic parameters 
  Double_t roots[3]= {0.0,0.0,0.0};  // cubic solutions
  Double_t value = 0.0;   // minimum eigenvalue
  Double_t x1,x2,x3;      // eigenvector component
  Double_t n1,n2,n3,nr= 0;// unit eigenvector 
  Double_t Radiusdm[7] = {0.3,0.4,0.7,1.49,2.38,3.91,4.36}; // beam pipe and layers radii [dm]
  Double_t sigma_MS = 0;
  TVector3 Vec,VecNor;

// Select RecPoints belonging to the track
  for(k =0;k<fPoints;k++){
    if(fPointRecs[k]->track==tracknumber) {
      name = fPointRecs[k]->fName;
      Pt = fPointRecs[k]->fPt;
      Layer = fPointRecs[k]->lay;
      Int_t ilay = Layer-1;
      if(npl[ilay]!=0) continue;
      if(bitlay[ilay] == 1) {
	original[N]->SetXYZ(0.1*fPointRecs[k]->fx,0.1*fPointRecs[k]->fy,0.1*fPointRecs[k]->fz);
	errors[N]->SetXYZ(0.1*fPointRecs[k]->fdx,0.1*fPointRecs[k]->fdy,0.1*fPointRecs[k]->fdz);
        sigma_MS = (Radiusdm[Layer]-Radiusdm[0])*0.000724/P;// beam pipe contribution
	for(Int_t j=1;j<Layer;j++) {
	  sigma_MS += (Radiusdm[Layer]-Radiusdm[j])*0.00136/P;
	}
	Weight[N] = ( 1 + original[N]->Perp2() )*( 1+ original[N]->Perp2() )/
	  ( errors[N]->Perp2() + sigma_MS*sigma_MS );
	linerrors[N]->SetXYZ(errors[N]->X(),errors[N]->Y(),sqrt(errors[N]->Z()*errors[N]->Z()+sigma_MS*sigma_MS));
	N++;
	npl[ilay]++;
      }  // end if on layer 
    }    //end if track==tracknumber
  }      //end for k
  //
  //    6 points, no more, no less
  //
  if(original[5]->X() == 9999 || original[6]->X() != 9999)  
    {
      delete [] Weight;
      return 1;   // not enough points
    }
  
  //
  //
  //
  //  FIT ON THE RIEMANN SPHERE FOR (x,y) PLANE
  //  
  
  RiemannTransf(N,original,Riemann);  

  Double_t Sum_Weights = 0.0; // sum of weights factor

  for(Int_t j=0;j<N;j++){ // mean values for x[i],y[i],z[i]
    xbar+=Weight[j]*Riemann[j]->X();
    ybar+=Weight[j]*Riemann[j]->Y();
    zbar+=Weight[j]*Riemann[j]->Z();
    Sum_Weights+=Weight[j];
  }
  
  xbar /= Sum_Weights;
  ybar /= Sum_Weights;
  zbar /= Sum_Weights;
  
  for(Int_t j=0;j<N;j++) {  // Calculate the matrix elements
    A11 += Weight[j]*(Riemann[j]->X() - xbar)*(Riemann[j]->X() - xbar);
    A12 += Weight[j]*(Riemann[j]->X() - xbar)*(Riemann[j]->Y() - ybar);
    A22 += Weight[j]*(Riemann[j]->Y() - ybar)*(Riemann[j]->Y() - ybar);
    A23 += Weight[j]*(Riemann[j]->Y() - ybar)*(Riemann[j]->Z() - zbar);
    A13 += Weight[j]*(Riemann[j]->X() - xbar)*(Riemann[j]->Z() - zbar);
    A33 += Weight[j]*(Riemann[j]->Z() - zbar)*(Riemann[j]->Z() - zbar);
  }
  
  A11 /= N;
  A12 /= N;
  A22 /= N;
  A23 /= N;
  A13 /= N;
  A33 /= N;
  A21 = A12;
  A32 = A23;
  A31 = A13;
  
// **************  Determinant  parameters ********************
//   n.b. simplifications done keeping in mind symmetry of A
//
  a = 1;
  b = (-A11-A33-A22);
  c = (A11*(A22+A33)+A33*A22-A12*A21-A13*A31-A23*A32);
  d = (A31*A22*A13+(A12*A21-A11*A22)*A33-2.0*A23*A13*A12+A11*A23*A32);

// **************  Find the 3 eigenvalues *************************
  Int_t Check_Cubic = SolveCubic(b,c,d,roots[0],roots[1],roots[2]);
 
  if(Check_Cubic !=1 ){
    printf("Track %d Has no real solution continuing ...\n",tracknumber);
    delete [] Weight;
    return 2;
  }
  
// **************** Find the lowest eigenvalue *****************
  if(roots[0]<=roots[1] && roots[0]<=roots[2]) value = roots[0];
  if(roots[1]<=roots[0] && roots[1]<=roots[2]) value = roots[1];
  if(roots[2]<=roots[0] && roots[2]<=roots[1]) value = roots[2];

  // ************ Eigenvector relative to value **************
  x3 = 1;
  x2 = (A33*A21-A23*A31-value*A21)/(A22*A31-A32*A21-value*A31);
  x1 = (value-A33-A32*x2)/A31;
  Vec.SetXYZ(x1,x2,x3);
  VecNor = Vec.Unit();
  n1 = VecNor.X();
  n2 = VecNor.Y();
  n3 = VecNor.Z();
  nr = -n1*xbar-n2*ybar-n3*zbar; 

  u0  = -0.5*n1/(nr+n3);
  v0  = -0.5*n2/(nr+n3);
  rho = sqrt((n1*n1 + n2*n2 -4*nr*(nr+n3))/(4*(nr+n3)*(nr+n3)));

  fFit = 0.0;
  fFit += 10.*TMath::Abs(sqrt((original[0]->X()-u0)*(original[0]->X()-u0)+(original[0]->Y()-v0)*(original[0]->Y()-v0))-rho);
  fFit += 10.*TMath::Abs(sqrt((original[1]->X()-u0)*(original[1]->X()-u0)+(original[1]->Y()-v0)*(original[1]->Y()-v0))-rho);
  fFit += 10.*TMath::Abs(sqrt((original[2]->X()-u0)*(original[2]->X()-u0)+(original[2]->Y()-v0)*(original[2]->Y()-v0))-rho);
  fFit += 10.*TMath::Abs(sqrt((original[3]->X()-u0)*(original[3]->X()-u0)+(original[3]->Y()-v0)*(original[3]->Y()-v0))-rho);
  fFit += 10.*TMath::Abs(sqrt((original[4]->X()-u0)*(original[4]->X()-u0)+(original[4]->Y()-v0)*(original[4]->Y()-v0))-rho);
  fFit += 10.*TMath::Abs(sqrt((original[5]->X()-u0)*(original[5]->X()-u0)+(original[5]->Y()-v0)*(original[5]->Y()-v0))-rho);

  fd0      =   100000.*(TMath::Sqrt(u0*u0+v0*v0)-rho); // transverse impact parameter in microns 
  fphi     =   TMath::ATan2(v0,u0);

//**************************************************************************
//  LINEAR FIT IN (z,s) PLANE:    z = zData.X() + zData.Y()*s
//  strictly linear (no approximation)
//**************************************************************************

       ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //                                                                                                        //
     //      REMEMBER, HERE STILL LENGHTS IN DM'S FOR ___INPUT___ BUT zDATA PARAMETERS ARE RETURNED IN CM'S    //
    //       rho, u0, v0 parameters converted right now to cm's... it's a mess, I'll take care, sometimes...  //
   //                                                                                                        //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  rho    *=  10.0;
  u0     *=  10.0;
  v0     *=  10.0;
  ierrl=FitLinear(N,original,linerrors,omega,u0,v0,fphi,zData,zError,CorrLin);
  chisql=zError.Z();
//   fprintf(pout,"%f \n",chisql); 
  z0=zData.X();
  vpar=zData.Y();
  fCorrLin = CorrLin;
  ierr = (ierrl > ierr ? ierrl : ierr);
//   fclose(pout);
  delete [] Weight;
  return ierr;
}

//-----------------------------------------------------------------------------
void AliITSRiemannFit::Streamer(TBuffer &lRb){
////////////////////////////////////////////////////////////////////////
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. 
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSRiemannFit.
  Int_t i,j,n;
  Int_t ii,jj;
  n=20;

  if (lRb.IsReading()) {
    Version_t lRv = lRb.ReadVersion(); if (lRv) { }
    TObject::Streamer(lRb);
    lRb >> fSizeEvent;
    lRb >> fPrimaryTracks;
    lRb >> fPoints;
    for(i=0;i<6;i++) lRb >> fPLay[i];    
    if(fPointRecs!=0){
      for(i=0;i<fSizeEvent;i++) delete[] fPointRecs[i];
      delete[] fPointRecs;
    } // end if fPointRecs!=0
    fPointRecs = new Point_tl*[fSizeEvent];
    for(i=0;i<fSizeEvent;i++){
      fPointRecs[i] = new Point_tl[n];
      for(j=0;j<n;j++){
	lRb >> fPointRecs[i][j].lay;
	lRb >> fPointRecs[i][j].lad;
	lRb >> fPointRecs[i][j].det;
	lRb >> fPointRecs[i][j].track;
	lRb >> fPointRecs[i][j].fx;
	lRb >> fPointRecs[i][j].fy;
	lRb >> fPointRecs[i][j].fz;
	lRb >> fPointRecs[i][j].fr;
	lRb >> fPointRecs[i][j].fdE;
	lRb >> fPointRecs[i][j].fdx;
	lRb >> fPointRecs[i][j].fdy;
	lRb >> fPointRecs[i][j].fdz;
	lRb >> fPointRecs[i][j].fPt;
	lRb >> (Char_t*)fPointRecs[i][j].fName;
	for (ii=0;ii<4;ii++)              
	  lRb << fPointRecs[i][j].fOrigin[ii];
	for (jj=0;jj<4;jj++)              
	  lRb << fPointRecs[i][j].fMomentum[jj];
	lRb >> fPointRecs[i][j].fCode;
	lRb >> fPointRecs[i][j].phi;
	lRb >> fPointRecs[i][j].eta;
	lRb >> fPointRecs[i][j].vertexPhi;
      } //end for j
    } //end for i
//      if(fspdi!=0){
//        for(i=0;i<fSizeEvent/6;i++) delete[] fspdi[i];
//        delete[] fspdi;
//      } // end if fspdi!=0
//      fspdi = new Point_tl*[fSizeEvent/6];
//      for(i=0;i<fSizeEvent/6;i++){
//        fspdi[i] = new Point_tl[n];
//        for(j=0;j<n;j++){
//  	lRb >> fspdi[i][j].lay;
//  	lRb >> fspdi[i][j].lad;
//  	lRb >> fspdi[i][j].det;
//  	lRb >> fspdi[i][j].track;
//  	lRb >> fspdi[i][j].fx;
//  	lRb >> fspdi[i][j].fy;
//  	lRb >> fspdi[i][j].fz;
//  	lRb >> fspdi[i][j].fr;
//  	lRb >> fspdi[i][j].fdE;
//  	lRb >> fspdi[i][j].fdx;
//  	lRb >> fspdi[i][j].fdy;
//  	lRb >> fspdi[i][j].fdz;
//  	lRb >> fspdi[i][j].fPt;
//  	for (ii=0;ii<4;ii++)              
//  	  lRb << fspdi[i][j].fOrigin[ii];
//  	for (jj=0;jj<4;jj++)              
//  	  lRb << fspdi[i][j].fMomentum[jj];
//  	lRb >> fspdi[i][j].fCode;
//  	lRb >> (Char_t*)fspdi[i][j].fName;
//  	lRb >> fspdi[i][j].phi;
//  	lRb >> fspdi[i][j].eta;
//  	lRb >> fspdi[i][j].vertexPhi;
//        } //end for j
//      } //end for i
//      if(fspdo!=0){
//        for(i=0;i<fSizeEvent/6;i++) delete[] fspdo[i];
//        delete[] fspdo;
//      } // end if fspdo!=0
//      fspdo = new Point_tl*[fSizeEvent/6];
//      for(i=0;i<fSizeEvent/6;i++){
//        fspdo[i] = new Point_tl[n];
//        for(j=0;j<n;j++){
//  	lRb >> fspdo[i][j].lay;
//  	lRb >> fspdo[i][j].lad;
//  	lRb >> fspdo[i][j].det;
//  	lRb >> fspdo[i][j].track;
//  	lRb >> fspdo[i][j].fx;
//  	lRb >> fspdo[i][j].fy;
//  	lRb >> fspdo[i][j].fz;
//  	lRb >> fspdo[i][j].fr;
//  	lRb >> fspdo[i][j].fdE;
//  	lRb >> fspdo[i][j].fdx;
//  	lRb >> fspdo[i][j].fdy;
//  	lRb >> fspdo[i][j].fdz;
//  	lRb >> fspdo[i][j].fPt;
//  	for (ii=0;ii<4;ii++)              
//  	  lRb << fspdo[i][j].fOrigin[ii];
//  	for (jj=0;jj<4;jj++)              
//  	  lRb << fspdo[i][j].fMomentum[jj];
//  	lRb >> fspdo[i][j].fCode;
//  	lRb >> (Char_t*)fspdo[i][j].fName;
//  	lRb >> fspdo[i][j].phi;
//  	lRb >> fspdo[i][j].eta;
//  	lRb >> fspdo[i][j].vertexPhi;
//        } //end for j
//      } //end for i
  } else {
    lRb.WriteVersion(AliITSRiemannFit::IsA());
    TObject::Streamer(lRb);
    lRb << fSizeEvent;
    lRb << fPrimaryTracks;
    lRb << fPoints;
    for(i=0;i<6;i++) lRb >> fPLay[i];
    for(i=0;i<fSizeEvent;i++) for(j=0;j<n;j++){
      lRb << fPointRecs[i][j].lay;
      lRb << fPointRecs[i][j].lad;
      lRb << fPointRecs[i][j].det;
      lRb << fPointRecs[i][j].track;
      lRb << fPointRecs[i][j].fx;
      lRb << fPointRecs[i][j].fy;
      lRb << fPointRecs[i][j].fz;
      lRb << fPointRecs[i][j].fr;
      lRb << fPointRecs[i][j].fdE;
      lRb << fPointRecs[i][j].fdx;
      lRb << fPointRecs[i][j].fdy;
      lRb << fPointRecs[i][j].fdz;
      lRb << fPointRecs[i][j].fPt;
      for (ii=0;ii<4;ii++)              
	lRb << fPointRecs[i][j].fOrigin[ii];
      for (jj=0;jj<4;jj++)              
	lRb << fPointRecs[i][j].fMomentum[jj];
      lRb << fPointRecs[i][j].fCode;
      lRb << fPointRecs[i][j].fName;
      lRb << fPointRecs[i][j].phi;
      lRb << fPointRecs[i][j].eta;
      lRb << fPointRecs[i][j].vertexPhi;
    }
//      for(i=0;i<fSizeEvent/6;i++) for(j=0;j<n;j++){
//        lRb << fspdi[i][j].lay;
//        lRb << fspdi[i][j].lad;
//        lRb << fspdi[i][j].det;
//        lRb << fspdi[i][j].track;
//        lRb << fspdi[i][j].fx;
//        lRb << fspdi[i][j].fy;
//        lRb << fspdi[i][j].fz;
//        lRb << fspdi[i][j].fr;
//        lRb << fspdi[i][j].fdE;
//        lRb << fspdi[i][j].fdx;
//        lRb << fspdi[i][j].fdy;
//        lRb << fspdi[i][j].fdz;
//        lRb << fspdi[i][j].fPt;
//        for (ii=0;ii<4;ii++)              
//  	lRb << fspdi[i][j].fOrigin[ii];
//        for (jj=0;jj<4;jj++)              
//  	lRb << fspdi[i][j].fMomentum[jj];
//        lRb << fspdi[i][j].fCode;
//        lRb << fspdi[i][j].fName;
//        lRb << fspdi[i][j].phi;
//        lRb << fspdi[i][j].eta;
//        lRb << fspdi[i][j].vertexPhi;
//      }
//      for(i=0;i<fSizeEvent/6;i++) for(j=0;j<n;j++){
//        lRb << fspdo[i][j].lay;
//        lRb << fspdo[i][j].lad;
//        lRb << fspdo[i][j].det;
//        lRb << fspdo[i][j].track;
//        lRb << fspdo[i][j].fx;
//        lRb << fspdo[i][j].fy;
//        lRb << fspdo[i][j].fz;
//        lRb << fspdo[i][j].fr;
//        lRb << fspdo[i][j].fdE;
//        lRb << fspdo[i][j].fdx;
//        lRb << fspdo[i][j].fdy;
//        lRb << fspdo[i][j].fdz;
//        lRb << fspdo[i][j].fPt;
//        for (ii=0;ii<4;ii++)              
//  	lRb << fspdo[i][j].fOrigin[ii];
//        for (jj=0;jj<4;jj++)              
//  	lRb << fspdo[i][j].fMomentum[jj];
//        lRb << fspdo[i][j].fCode;
//        lRb << fspdo[i][j].fName;
//        lRb << fspdo[i][j].phi;
//        lRb << fspdo[i][j].eta;
//        lRb << fspdo[i][j].vertexPhi;
//    }
  } // end if reading
}
