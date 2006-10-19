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
 **************************************************************************/

// 
//                                                                        *
// This class performs a fast fit of helices going through the <=6        *
// points of the ITS, with the goal of studying tracking and              *
// vertexing performances.                                                *
// Generated kinematics is used to take into account different weights    *
// associated to points in different layers (with different multiple      *
// scattering-originated errors).                                         *
//                                                                        *
//   Based on the work by A. Strandlie, R. Fruhwirth                      *
//                                                                        *
//   First implementation by N. Bustreo, R. Turrisi - July 2000           *
//                                                                        *
//   Further modifications by A. Dainese, R. Turrisi                      *
//                                                                        *
//   Contact: Rosario Turrisi, rosario.turrisi@pd.infn.it                 *
//                                                                        *
// ************************************************************************
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



#include "AliITSRiemannFit.h"
#include "AliRun.h"
#include "TClonesArray.h"
#include "stdio.h"
#include "stdlib.h"
#include "Riostream.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TParticle.h"
#include "TTree.h"
#include "TVector3.h"
#include "AliITSRecPoint.h"
#include "AliITSgeom.h"
#include "AliMC.h"
#include "AliITSDetTypeRec.h"
#include "AliLog.h"

ClassImp(AliITSRiemannFit)


AliITSRiemannFit::AliITSRiemannFit():
fSizeEvent(0),
fPrimaryTracks(0),
fPoints(0),
fParticles(0),
fPointRecs(0) {
  ///////////////////////////////////////////////////////////
  // Default constructor.
  // Set everything to zero.
  ////////////////////////////////////////////////////////////


  for(Int_t i=0;i<6;i++)fPLay[i] = 0;
  
}

//______________________________________________________________________
AliITSRiemannFit::AliITSRiemannFit(const AliITSRiemannFit &rf) : TObject(rf),
fSizeEvent(rf.fSizeEvent),
fPrimaryTracks(rf.fPrimaryTracks),
fPoints(rf.fPoints),
fParticles(rf.fParticles),
fPointRecs(rf.fPointRecs) {
  // Copy constructor

}

//______________________________________________________________________
AliITSRiemannFit& AliITSRiemannFit::operator=(const AliITSRiemannFit& rf){
  // Assignment operator
  this->~AliITSRiemannFit();
  new(this) AliITSRiemannFit(rf);
  return *this;

}

//______________________________________________________________________
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

AliITSRiemannFit::AliITSRiemannFit(Int_t size,Int_t ntracks):
fSizeEvent(size),
fPrimaryTracks(ntracks),
fPoints(0),
fParticles(0),
fPointRecs(0)  {
  ///////////////////////////////////////////////////////////
  // Constructor.
  // Set fSizeEvent to size and fPrimaryTracks to ntracks.
  // Others to zero.
  ////////////////////////////////////////////////////////////


  AliPointtl *first = new AliPointtl[fSizeEvent];
  AliPointtl **pointRecs = new AliPointtl*[fSizeEvent];
  for(Int_t i=0;i<6;i++)fPLay[i] = 0;
  for(Int_t j=0;j<fSizeEvent;j++)  // create an array of struct
    pointRecs[j] = &(first[j]);   
}

// ---------------------------------------------------------------------
AliITSRiemannFit::AliPointtl::AliPointtl():
fLay(0),
fLad(0),
fDet(0),
fTrack(0),
fx(0),
fy(0),
fz(0),
fr(0),
fdE(0),
fdx(0),
fdy(0),
fdz(0),
fOrigin(0),
fMomentum(0),
fCode(0),
fName(0),
fPt(0),
fPhi(0),
fEta(0),fVertexPhi(0){
  // default constructor
  SetLay();
  SetLad();
  SetDet();
  SetTrack();
  SetX();
  SetY();
  SetZ();
  SetR();
  SetdE();
  SetdX();
  SetdY();
  SetdZ();
  SetOrigin();
  SetMomentum();
  SetCode();
  SetName();
  SetPt();
  SetPhi();
  SetEta();
  SetVertexPhi();
}

AliITSRiemannFit::AliPointtl::AliPointtl(const AliPointtl& ap):
fLay(ap.fLay),
fLad(ap.fLad),
fDet(ap.fDet),
fTrack(ap.fTrack),
fx(ap.fx),
fy(ap.fy),
fz(ap.fz),
fr(ap.fr),
fdE(ap.fdE),
fdx(ap.fdx),
fdy(ap.fdy),
fdz(ap.fdz),
fOrigin(ap.fOrigin),
fMomentum(ap.fMomentum),
fCode(ap.fCode),
fName(ap.fName),
fPt(ap.fPt),
fPhi(ap.fPhi),
fEta(ap.fEta),
fVertexPhi(ap.fVertexPhi){
  //copy constructor
}




// ---------------------------------------------------------------------

void FillPoints(AliITSRiemannFit::AliPointtl **Points,Int_t &index,Float_t *xpoint,
		Float_t *error,
		TLorentzVector pE,TLorentzVector oT,Int_t *id,
		Int_t track, Char_t *name,Int_t code,
		Float_t phiorigin){
  ///////////////////////////////////////////////////////////////////////
  // Fill the structure AliPointtl with the proper data
  //
  //////////////////////////////////////////////////////////////////////
  Float_t pPI2 = 2.0*TMath::Pi();
  Float_t phi,r,x,y,z;
  Int_t i;
  i = index;
  x = xpoint[0];
  y = xpoint[1];
  z = xpoint[2];
  r = sqrt(x*x+y*y);
  phi = TMath::ATan2(y,x);
  if(phi<0.0) phi += pPI2;
  Points[i]->SetPhi(phi);
  Points[i]->SetEta(-0.5*tan(0.5*TMath::ATan2(r,z)));
  Points[i]->SetX(x);
  Points[i]->SetY(y);
  Points[i]->SetZ(z);
  Points[i]->SetdX(error[0]);
  Points[i]->SetdY(error[1]);
  Points[i]->SetdZ(error[2]);
  Points[i]->SetR(r);
  Points[i]->SetTrack(track);
  Points[i]->SetLay(id[0]);
  Points[i]->SetLad(id[1]);
  Points[i]->SetDet(id[2]);
  Points[i]->SetMomentum(&pE);
  Points[i]->SetOrigin(&oT);
  Points[i]->SetPt(sqrt(pE.X()*pE.X()+pE.Y()*pE.Y()));
  Points[i]->SetCode(code);
  Points[i]->SetName(name);
  Points[i]->SetVertexPhi(phiorigin);
  index++;
  return;
  
}
// -----------------------------------------------------------------------

void AliITSRiemannFit::InitPoints(Int_t ntracks,TTree *TR,Int_t nparticles){
  //////////////////////////////////////////////////////////////////////
  // Fill the class member fPointRecs with the reconstructed points
  // Set All other members to the real values
  //
  /////////////////////////////////////////////////////////////////////
  printf("\n ************* Starting Init Points *************\n");
  TParticle *part;

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  rl->CdGAFile();
  AliITSLoader* loader = static_cast<AliITSLoader*>(rl->GetLoader("ITSLoader"));
  if (!loader) {
    Error("InitPoints", "ITS loader not found");
    return;
  }
  AliITSgeom* gm = loader->GetITSgeom();

  //get pointer to modules array
  Int_t nmodules = gm->GetIndexMax();
  // Get the points from points file
  Int_t mod,irec;
  Stat_t nent;
  AliITSRecPoint  *recp;
  nent=TR->GetEntries();
  AliITSDetTypeRec detTypeRec;
  TClonesArray *iTSrec  = detTypeRec.RecPoints();
  Int_t totRP=0;
  for (mod=0; mod<nmodules; mod++) {
    detTypeRec.ResetRecPoints();
    TR->GetEvent(mod);
    Int_t nrecp = iTSrec->GetEntries();
    if(!nrecp) continue;
    totRP += nrecp;
  }

  Int_t iMAX = totRP;
  fPrimaryTracks = ntracks;
  fParticles     = nparticles;
  AliITSRiemannFit::AliPointtl *global = new AliPointtl[iMAX];
  fPointRecs = new AliITSRiemannFit::AliPointtl*[iMAX];
  //
  for(Int_t j=0;j<iMAX;j++) {
    fPointRecs[j] = &(global[j]);
  }
  
  Int_t ieta=0,ieta2=0;
  Int_t i,id[4],idold[4];
  Int_t track=0;//         // track  of hit
  Float_t xpoint[3],errorPlus[3],errorMinus[3],globalError[3];       // position and error of the point
  TLorentzVector oT,pE;
  Float_t locals[3],localserror[3],localsplus[3],localsminus[3]; // local position and local errors
  Float_t pPhi;
  Int_t code;
  const char *name;
  Int_t layer,ladder,detector;
  Float_t xcluster,zcluster;
  Int_t num=0,nspdi=0,nspdo=0,nsddi=0,nsddo=0,nssdi=0,nssdo=0;
 
  for (mod=0; mod<nmodules; mod++) {
    //itsModule=(AliITSmodule*)iTSmodules->At(mod);
    //ITS->ResetRecPoints();
    detTypeRec.ResetRecPoints();
    TR->GetEvent(mod);
    Int_t nrecp = iTSrec->GetEntries();
    if (!nrecp) continue;
    //itsModule->GetID(layer,ladder,detector);
    gm->GetModuleId(mod,layer,ladder,detector);

    for (irec=0;irec<nrecp;irec++) {
      recp   = (AliITSRecPoint*)iTSrec->UncheckedAt(irec);
      track=recp->GetLabel(0);
      if(track <0 ) continue;
      xcluster=recp->GetDetLocalX();     // x on cluster
      zcluster=recp->GetDetLocalZ();     // z on cluster
      part   = (TParticle*) gAlice->GetMCApp()->Particle(track);    
      part->ProductionVertex(oT);  // set the vertex 
      part->Momentum(pE);          // set the vertex momentum
      name      = part->GetName();
      Char_t nam2[50];
      sprintf(nam2,"%s",name); 
      code      = part->GetPdgCode();
      pPhi       = part->Phi();
      id[0]=layer;
      id[1]=ladder;
      id[2]=detector;
      id[3]=irec;
      locals[0]=xcluster;     // x on cluster
      locals[1]=0.0;          // y on cluster
      locals[2]=zcluster;     // z on cluster
      localserror[0]=sqrt(recp->GetSigmaDetLocX2());
      localserror[1]=0.0;
      localserror[2]=sqrt(recp->GetSigmaZ2());
      localsplus[0]=xcluster+sqrt(recp->GetSigmaDetLocX2());       // x on cluster
      if(layer==1||layer==2) localsplus[1]=0.0150/2;         // y on cluster
      else if(layer==3||layer==4) localsplus[1]=0.0280/2;    // y on cluster
      else if(layer==5||layer==6) localsplus[1]=0.0300/2;    // y on cluster
      localsplus[2]=zcluster+sqrt(recp->GetSigmaZ2());       // z on cluster
      localsminus[0]=xcluster-sqrt(recp->GetSigmaDetLocX2());      // x on cluster
      localsminus[1]=0.0;                                    // y on cluster
      localsminus[2]=zcluster-sqrt(recp->GetSigmaZ2());      // z on cluster

      gm->LtoG(layer,ladder,detector,locals,xpoint);
      gm->LtoG(layer,ladder,detector,localsplus,errorPlus);
      gm->LtoG(layer,ladder,detector,localsminus,errorMinus);
      globalError[0]=0.5*TMath::Abs(errorPlus[0]-errorMinus[0]);
      globalError[1]=0.5*TMath::Abs(errorPlus[1]-errorMinus[1]);
      globalError[2]=0.5*TMath::Abs(errorPlus[2]-errorMinus[2]);
      if(track<ntracks) {
	if(TMath::Abs(part->Eta())<=1.0) ieta++;
	if(TMath::Abs(part->Eta())<=0.5) ieta2++;
      }
      if(!(id[0]==idold[0]&&id[1]==idold[1]&&
	   id[2]==idold[2]&&id[3]==idold[3])) {
	FillPoints(fPointRecs,num,xpoint,globalError,pE,oT,id,track,nam2,code,pPhi);
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
//  	  FillPoints(fspdi,nspdi,xpoint,globalError,pE,oT,id,track,name,code,pPhi);  
//  	}
//  	if(idold[0]==2){
	 
//  	  FillPoints(fspdo,nspdo,xpoint,globalError,pE,oT,id,track,name,code,pPhi);
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

  delete rl;
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
Bool_t SortZ(const AliITSRiemannFit::AliPointtl *s1,const AliITSRiemannFit::AliPointtl *s2){
  // Z sorting function for qsort.
   Float_t a;

   a = s1->GetZ() - s2->GetZ();
   if(a<0.0) return kTRUE;
   if(a>0.0) return kFALSE;
   return kFALSE;
}
Bool_t SortTrack(const AliITSRiemannFit::AliPointtl *s1,const AliITSRiemannFit::AliPointtl *s2){
  // track sorting function for qsort.
   Float_t a;

   a = s1->GetTrack() - s2->GetTrack();
   if(a<0.0) return kTRUE;
   if(a>0.0) return kFALSE;
   return kFALSE;
}
void hpsortTrack(AliITSRiemannFit::AliPointtl **ra,Int_t n){
   Int_t i,ir,j,l;
   AliITSRiemannFit::AliPointtl *rra;

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
void hpsortZ(AliITSRiemannFit::AliPointtl **ra,Int_t n){
   Int_t i,ir,j,l;
   AliITSRiemannFit::AliPointtl *rra;

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
    fprintf(ascii,"%d\t%d\t%f\t%f\t%f\n",fPointRecs[i]->GetLay(),
	    fPointRecs[i]->GetTrack(),fPointRecs[i]->GetX(),
            fPointRecs[i]->GetY(),fPointRecs[i]->GetZ());
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
	   i,fPointRecs[i]->GetLay(),fPointRecs[i]->GetTrack(),
           fPointRecs[i]->GetX(),fPointRecs[i]->GetY(),
           fPointRecs[i]->GetZ(),fPointRecs[i]->GetOrigin()->X(),
	   fPointRecs[i]->GetOrigin()->Y(),fPointRecs[i]->GetOrigin()->Z(),
	   fPointRecs[i]->GetPt(),fPointRecs[i]->GetName());
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
   
   Double_t qQ = ((a*a - 3*b)/9);
   Double_t rR = ((2*a*a*a - 9*a*b +27*c)/54);
   Double_t theta;
   Double_t aF = -2*sqrt(qQ);
   Double_t g = a/3;
   Double_t pPI2 = TMath::Pi()*2;

  if( rR*rR>qQ*qQ*qQ ) {
    cout<<"\nTrack "<<"Determinant :\n\t\t No Real Solutions !!!\n"<<endl;
    x1 = 9999999;
    x2 = 9999999;
    x3 = 9999999;
    return 0;
  }

  theta = TMath::ACos(rR/sqrt(qQ*qQ*qQ));

  x1 = (aF*TMath::Cos(theta/3))-g;
  x2 = (aF*TMath::Cos((theta+pPI2)/3))-g;
  x3 = (aF*TMath::Cos((theta-pPI2)/3))-g;
  
  return 1;
}
//-----------------------------------------------------------------

void RiemannTransf(Int_t npoints,TVector3 **From,TVector3 **To) {
  ///////////////////////////////////////////////////////////////////////
  //   This function apllies the transformation in the Riemann sphere
  //   for xy plane
  ///////////////////////////////////////////////////////////////////////
  Float_t *rR = new Float_t[npoints];
  Float_t *theta = new Float_t[npoints];
  Float_t pPI2 = 2*TMath::Pi();
  Float_t x=0,y=0,z=0;
  
  for(Int_t i=0;i<npoints;i++) {
    rR[i]     = sqrt(From[i]->X()*From[i]->X()+From[i]->Y()*From[i]->Y());
    theta[i] = TMath::ATan2(From[i]->Y(),From[i]->X());
    if(theta[i]<0) theta[i]+=pPI2;
    x     = rR[i]*cos(theta[i])/(1+rR[i]*rR[i]);
    y     = rR[i]*sin(theta[i])/(1+rR[i]*rR[i]);
    z     = rR[i]*rR[i]/(1+rR[i]*rR[i]);
    To[i]->SetXYZ(x,y,z);
  }
  delete[] rR;
  delete[] theta;
  return;
}


//---------------------------------------------------------------------

Int_t FitLinear(Int_t npoints, TVector3 **input, TVector3 **errors, Double_t omega,
		Double_t &thu0, Double_t &thv0, Double_t &phi, TVector2 &zData, TVector3 &zError, 
		Double_t &corrLin){
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
  fitHist->Fit("pol1","qQ");
  z0  = fitHist->GetFunction("pol1")->GetParameter(0);
  vpar  = fitHist->GetFunction("pol1")->GetParameter(1);
  ez0 = fitHist->GetFunction("pol1")->GetParError(0);
  evpar = fitHist->GetFunction("pol1")->GetParError(1);
  chisquare = fitHist->GetFunction("pol1")->GetChisquare();
  zData.Set(z0,vpar);
  zError.SetXYZ(ez0,evpar,chisquare);

  Double_t sigmas=0.; 
  Double_t sigmaz=0.;
  Double_t avs=0.;
  Double_t avz=0.;
  Double_t avsz=0.;

  for(Int_t j = 0; j < npoints; j++) {
    avs  += s[j];
    avz  += z[j]; 
    avsz += s[j]*z[j];
  }
    avs  /= (Double_t)npoints;
    avz  /= (Double_t)npoints; 
    avsz /= (Double_t)npoints;

  for(Int_t l = 0; l < npoints; l++) {
    sigmas += (s[l]-avs)*(s[l]-avs);
    sigmaz += (z[l]-avz)*(z[l]-avz);
  }
  sigmas /=(Double_t)npoints;
  sigmaz /=(Double_t)npoints;

  sigmas = sqrt(sigmas);
  sigmaz = sqrt(sigmaz);

  corrLin = (avsz-avs*avz)/(sigmas*sigmaz);

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
  Int_t pLayer = 0;
  const Char_t* name = 0;
  Int_t i=0,k=0;
  Int_t iMAX = 50;
  Int_t nN = 0;
  Int_t npl[6]={0,0,0,0,0,0};
  Double_t pP = sqrt(Px*Px+Py*Py+Pz*Pz);
  Double_t pt = sqrt(Px*Px+Py*Py);
  TVector3 zError;
  TVector2 zData;
  Double_t corrLin;
  TVector3 *ori        = new TVector3[iMAX];
  TVector3 **original  = new TVector3*[iMAX];
  TVector3 *rie       = new TVector3[iMAX];
  TVector3 **riemann  = new TVector3*[iMAX];
  TVector3 *err        = new TVector3[iMAX];
  TVector3 **errors    = new TVector3*[iMAX];
  TVector3 *linerr        = new TVector3[iMAX];
  TVector3 **linerrors    = new TVector3*[iMAX];
  //PH  Double_t weight[iMAX];
  Double_t * weight = new Double_t[iMAX];

  for(i=0;i<iMAX;i++){
    original[i]   = &(ori[i]);
    riemann[i]   = &(rie[i]);
    errors[i]     = &(err[i]);
    linerrors[i]  = &(linerr[i]);
  }
  for(k =0;k<iMAX;k++) original[k]->SetXYZ(9999,9999,9999);
  Double_t a11,a12,a13,a21,a22,a23,a31,a32,a33;
  a11=0;a12=0;a13=0;a21=0;a22=0;a23=0;a31=0;a32=0;a33=0;
  Double_t xbar = 0;
  Double_t ybar = 0;
  Double_t zbar = 0;
  Double_t a,b,c,d;       // cubic parameters 
  Double_t roots[3]= {0.0,0.0,0.0};  // cubic solutions
  Double_t value = 0.0;   // minimum eigenvalue
  Double_t x1,x2,x3;      // eigenvector component
  Double_t n1,n2,n3,nr= 0;// unit eigenvector 
  Double_t radiusdm[7] = {0.3,0.4,0.7,1.49,2.38,3.91,4.36}; // beam pipe and layers radii [dm]
  Double_t sigmaMS = 0;
  TVector3 vVec,vVecNor;

// Select RecPoints belonging to the track
  for(k =0;k<fPoints;k++){
    if(fPointRecs[k]->GetTrack()==tracknumber) {
      name = fPointRecs[k]->GetName();
      pt = fPointRecs[k]->GetPt();
      pLayer = fPointRecs[k]->GetLay();
      Int_t ilay = pLayer-1;
      if(npl[ilay]!=0) continue;
      if(bitlay[ilay] == 1) {
	original[nN]->SetXYZ(0.1*fPointRecs[k]->GetX(),0.1*fPointRecs[k]->GetY(),0.1*fPointRecs[k]->GetZ());
	errors[nN]->SetXYZ(0.1*fPointRecs[k]->GetdX(),0.1*fPointRecs[k]->GetdY(),0.1*fPointRecs[k]->GetdZ());
        sigmaMS = (radiusdm[pLayer]-radiusdm[0])*0.000724/pP;// beam pipe contribution
	for(Int_t j=1;j<pLayer;j++) {
	  sigmaMS += (radiusdm[pLayer]-radiusdm[j])*0.00136/pP;
	}
	weight[nN] = ( 1 + original[nN]->Perp2() )*( 1+ original[nN]->Perp2() )/
	  ( errors[nN]->Perp2() + sigmaMS*sigmaMS );
	linerrors[nN]->SetXYZ(errors[nN]->X(),errors[nN]->Y(),sqrt(errors[nN]->Z()*errors[nN]->Z()+sigmaMS*sigmaMS));
	nN++;
	npl[ilay]++;
      }  // end if on layer 
    }    //end if track==tracknumber
  }      //end for k
  //
  //    6 points, no more, no less
  //
  if(original[5]->X() == 9999 || original[6]->X() != 9999)  
    {
      delete [] weight;
      return 1;   // not enough points
    }
  
  //
  //
  //
  //  FIT ON THE RIEMANN SPHERE FOR (x,y) PLANE
  //  
  
  RiemannTransf(nN,original,riemann);  

  Double_t sumWeights = 0.0; // sum of weights factor

  for(Int_t j=0;j<nN;j++){ // mean values for x[i],y[i],z[i]
    xbar+=weight[j]*riemann[j]->X();
    ybar+=weight[j]*riemann[j]->Y();
    zbar+=weight[j]*riemann[j]->Z();
    sumWeights+=weight[j];
  }
  
  xbar /= sumWeights;
  ybar /= sumWeights;
  zbar /= sumWeights;
  
  for(Int_t j=0;j<nN;j++) {  // Calculate the matrix elements
    a11 += weight[j]*(riemann[j]->X() - xbar)*(riemann[j]->X() - xbar);
    a12 += weight[j]*(riemann[j]->X() - xbar)*(riemann[j]->Y() - ybar);
    a22 += weight[j]*(riemann[j]->Y() - ybar)*(riemann[j]->Y() - ybar);
    a23 += weight[j]*(riemann[j]->Y() - ybar)*(riemann[j]->Z() - zbar);
    a13 += weight[j]*(riemann[j]->X() - xbar)*(riemann[j]->Z() - zbar);
    a33 += weight[j]*(riemann[j]->Z() - zbar)*(riemann[j]->Z() - zbar);
  }
  
  a11 /= nN;
  a12 /= nN;
  a22 /= nN;
  a23 /= nN;
  a13 /= nN;
  a33 /= nN;
  a21 = a12;
  a32 = a23;
  a31 = a13;
  
// **************  Determinant  parameters ********************
//   n.b. simplifications done keeping in mind symmetry of A
//
  a = 1;
  b = (-a11-a33-a22);
  c = (a11*(a22+a33)+a33*a22-a12*a21-a13*a31-a23*a32);
  d = (a31*a22*a13+(a12*a21-a11*a22)*a33-2.0*a23*a13*a12+a11*a23*a32);

// **************  Find the 3 eigenvalues *************************
  Int_t checkCubic = SolveCubic(b,c,d,roots[0],roots[1],roots[2]);
 
  if(checkCubic !=1 ){
    printf("Track %d Has no real solution continuing ...\n",tracknumber);
    delete [] weight;
    return 2;
  }
  
// **************** Find the lowest eigenvalue *****************
  if(roots[0]<=roots[1] && roots[0]<=roots[2]) value = roots[0];
  if(roots[1]<=roots[0] && roots[1]<=roots[2]) value = roots[1];
  if(roots[2]<=roots[0] && roots[2]<=roots[1]) value = roots[2];

  // ************ Eigenvector relative to value **************
  x3 = 1;
  x2 = (a33*a21-a23*a31-value*a21)/(a22*a31-a32*a21-value*a31);
  x1 = (value-a33-a32*x2)/a31;
  vVec.SetXYZ(x1,x2,x3);
  vVecNor = vVec.Unit();
  n1 = vVecNor.X();
  n2 = vVecNor.Y();
  n3 = vVecNor.Z();
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
  ierrl=FitLinear(nN,original,linerrors,omega,u0,v0,fphi,zData,zError,corrLin);
  chisql=zError.Z();
//   fprintf(pout,"%f \n",chisql); 
  z0=zData.X();
  vpar=zData.Y();
  fCorrLin = corrLin;
  ierr = (ierrl > ierr ? ierrl : ierr);
//   fclose(pout);
  delete [] weight;
  return ierr;
}
Int_t AliITSRiemannFit::FitHelix(Int_t NPoints, TVector3** fPointRecs,TVector3** fPointRecErrors,Float_t& f1, Float_t& f2, Float_t& f3) {

  ///////////////////////////////////////////////////////////////////////
  //  This function finds the helix parameters 
  //  d0  = impact parameter
  //  rho = radius of circle
  //  phi = atan(y0/x0)
  //  for the xy plane
  //  starting from the momentum and the outcome of
  //  the fit on the Riemann sphere (i.e. u0,v0,rho)
  //
  //   MIND !!!!   Here we assume both angular velocities be 1.0e-2  (yes, 0.01 !)
  //
  //
  //   Also linear fit in (z,s) is performed, so it's 3-D !
  //   z0 and vpar are calculated (intercept and z-component of velocity, but 
  //   in units... you guess.
  //
  //
  //   Values calculated in addition:
  //
  //              - transverse impact parameter        fd0
  //              - sum of residuals in (x,y) plane    fFit
  //              - chisquare of linear fit            chisql
  //              - correlation coefficient            fCorrLin
  //
  //
  //
  //
  //
  ///////////////////////////////////////////////////////////////////////
  //
  //  All this stuff relies on this hypothesis !!!
  //
  Int_t ierr = 0, ierrl=0;
  const Double_t kOmega = 1.0e-2;

  


  Double_t  fd0 = -9999;      //  fake values
  Double_t  u0 = -9999.9999;  //  for eventual 
  Double_t  v0 = -9999.9999;  //  debugging
  Double_t  rho = -9999.9999; //
  Double_t fphi, fFit, chisql, z0, vpar, fCorrLin;

  //
  //  This info is no more there... to be re-considered... maybe
  //
  //   Double_t pP = sqrt(Px*Px+Py*Py+Pz*Pz);
  //   Double_t pt = sqrt(Px*Px+Py*Py);
  
  TVector3 zError;
  TVector2 zData;
  Double_t corrLin;
  TVector3 *ori        = new TVector3[NPoints];
  TVector3 **original  = new TVector3*[NPoints];
  TVector3 *rie        = new TVector3[NPoints];
  TVector3 **riemann   = new TVector3*[NPoints];
  TVector3 *err        = new TVector3[NPoints];
  TVector3 **errors    = new TVector3*[NPoints];
  TVector3 *linerr     = new TVector3[NPoints];
  TVector3 **linerrors = new TVector3*[NPoints];
  Double_t * weight = new Double_t[NPoints];
  
  for(Int_t i=0; i<NPoints; i++){

    original[i]   = &(ori[i]);
    riemann[i]   = &(rie[i]);
    errors[i]     = &(err[i]);
    linerrors[i]  = &(linerr[i]);

    original[i]->SetXYZ(9999,9999,9999);
  }

  //
  //  Riemann fit parameters
  //
  Double_t a11,a12,a13,a21,a22,a23,a31,a32,a33;
  a11=0;a12=0;a13=0;a21=0;a22=0;a23=0;a31=0;a32=0;a33=0;
  Double_t xbar = 0;
  Double_t ybar = 0;
  Double_t zbar = 0;
  //
  Double_t a,b,c,d;                  // cubic parameters 
  Double_t roots[3]= {0.0,0.0,0.0};  // cubic solutions
  Double_t value = 0.0;              // minimum eigenvalue
  Double_t x1,x2,x3;                 // eigenvector component
  Double_t n1,n2,n3,nr= 0;           // unit eigenvector 
  TVector3 vVec,vVecNor;
  
  for (Int_t ip=0; ip<NPoints; ip++) {
original[ip]->SetXYZ(0.1*fPointRecs[ip]->X(),0.1*fPointRecs[ip]->Y(),0.1*fPointRecs[ip]->Z());
   
errors[ip]->SetXYZ(0.1*fPointRecErrors[ip]->X(),0.1*fPointRecErrors[ip]->Y(),0.1*fPointRecErrors[ip]->Z());
    weight[ip] = (1+original[ip]->Perp2())*(1+original[ip]->Perp2())/(errors[ip]->Perp2());
    linerrors[ip]->SetXYZ(errors[ip]->X(),errors[ip]->Y(),errors[ip]->Z());
  }


  //
  //
  //  FIT ON THE RIEMANN SPHERE FOR (x,y) PLANE
  //  
  
  RiemannTransf(NPoints,original,riemann);  

  Double_t sumWeights = 0.0; // sum of weights factor

  for(Int_t j=0;j<NPoints;j++){ // mean values for x[i],y[i],z[i]
    xbar+=weight[j]*riemann[j]->X();
    ybar+=weight[j]*riemann[j]->Y();
    zbar+=weight[j]*riemann[j]->Z();
    sumWeights+=weight[j];
  }
  
  xbar /= sumWeights;
  ybar /= sumWeights;
  zbar /= sumWeights;
  
  for(Int_t j=0;j<NPoints;j++) {  // Calculate the matrix elements
    a11 += weight[j]*(riemann[j]->X() - xbar)*(riemann[j]->X() - xbar);
    a12 += weight[j]*(riemann[j]->X() - xbar)*(riemann[j]->Y() - ybar);
    a22 += weight[j]*(riemann[j]->Y() - ybar)*(riemann[j]->Y() - ybar);
    a23 += weight[j]*(riemann[j]->Y() - ybar)*(riemann[j]->Z() - zbar);
    a13 += weight[j]*(riemann[j]->X() - xbar)*(riemann[j]->Z() - zbar);
    a33 += weight[j]*(riemann[j]->Z() - zbar)*(riemann[j]->Z() - zbar);
  }
  //
  // this doesn't seem to work...
  //
  //   a11 /= sumWeights;
  //   a12 /= sumWeights;
  //   a22 /= sumWeights;
  //   a23 /= sumWeights;
  //   a13 /= sumWeights;
  //   a33 /= sumWeights;

  a11 /= NPoints;
  a12 /= NPoints;
  a22 /= NPoints;
  a23 /= NPoints;
  a13 /= NPoints;
  a33 /= NPoints;
  a21 = a12;
  a32 = a23;
  a31 = a13;
  
// **************  Determinant  parameters ********************
//   n.b. simplifications done keeping in mind symmetry of A
//
  a = 1;
  b = (-a11-a33-a22);
  c = (a11*(a22+a33)+a33*a22-a12*a21-a13*a31-a23*a32);
  d = (a31*a22*a13+(a12*a21-a11*a22)*a33-2.0*a23*a13*a12+a11*a23*a32);

// **************  Find the 3 eigenvalues *************************
  Int_t checkCubic = SolveCubic(b,c,d,roots[0],roots[1],roots[2]);
 
  if(checkCubic !=1 ){
    printf("No real solution. Check data.\n");
    delete [] weight;
    return 999;
  }
  
// **************** Find the lowest eigenvalue *****************
  if(roots[0]<=roots[1] && roots[0]<=roots[2]) value = roots[0];
  if(roots[1]<=roots[0] && roots[1]<=roots[2]) value = roots[1];
  if(roots[2]<=roots[0] && roots[2]<=roots[1]) value = roots[2];

  // ************ Eigenvector relative to value **************
  x3 = 1;
  x2 = (a33*a21-a23*a31-value*a21)/(a22*a31-a32*a21-value*a31);
  x1 = (value-a33-a32*x2)/a31;
  vVec.SetXYZ(x1,x2,x3);
  vVecNor = vVec.Unit();
  n1 = vVecNor.X();
  n2 = vVecNor.Y();
  n3 = vVecNor.Z();
  nr = -n1*xbar-n2*ybar-n3*zbar; 

  u0  = -0.5*n1/(nr+n3);
  v0  = -0.5*n2/(nr+n3);
  rho = sqrt((n1*n1 + n2*n2 -4*nr*(nr+n3))/(4*(nr+n3)*(nr+n3)));
  

  fFit = 0.0;
  for (Int_t i=0; i<NPoints; i++) {
  fFit += 10.*TMath::Abs(sqrt((original[i]->X()-u0)*(original[i]->X()-u0)+(original[i]->Y()-v0)*(original[i]->Y()-v0))-rho);
  }
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
  
  ierrl=LinearFit(NPoints,original,linerrors,kOmega,u0,v0,fphi,zData,zError,corrLin);
  if(ierrl==33) return 0;
  chisql=zError.Z();
//   fprintf(pout,"%f \n",chisql); 
  z0=zData.X();
  vpar=zData.Y();
  fCorrLin = corrLin;
  ierr = (ierrl > ierr ? ierrl : ierr);
//   fclose(pout);
  delete [] weight;
  
  f1=fphi;
  f2=vpar/(kOmega*TMath::Abs(rho));
  f3=1/rho;
  delete[] ori;
  delete[] rie;
  delete[] err;
  delete[] linerr;
  delete[] original;
  delete[] riemann;
  delete[] errors;
  delete[] linerrors;
  
  return 1;
  

}

//____________________________________________________________

Int_t AliITSRiemannFit::LinearFit(Int_t npoints, TVector3 **input, 
                                 TVector3 **errors, Double_t omega,
                                 Double_t &thu0, Double_t &thv0, Double_t &phi,TVector2 &zData, TVector3 &zError, 
                                 Double_t &corrLin){
  ///////////////////////////////////////////////////////////////////////
  //   Fit the points in the (z,s) plane - helix 3rd equation
  //
  /////////////////////////////////////////////////////////////////////// 
  //By R.Turrisi

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
  
  //  if(s[0]>-636 && s[0]<-625) return 33;

  TGraph* fitHist = new TGraph(npoints,s,z);
  TF1* f1 = new TF1("f1",Fitfunction,-100,100,2);

  f1->SetParameter(0,1);
  f1->SetParameter(1,1);
  f1->SetLineColor(2);
  fitHist->Fit(f1,"qQ"); 
  
  z0   = f1->GetParameter(0);
  vpar = f1->GetParameter(1);
  ez0  = f1->GetParError(0);
  evpar= f1->GetParError(1);
  chisquare=f1->GetChisquare();
  zData.Set(z0,vpar);
  zError.SetXYZ(ez0,evpar,chisquare);
  
  Double_t sigmas=0.; 
  Double_t sigmaz=0.;
  Double_t avs=0.;
  Double_t avz=0.;
  Double_t avsz=0.;

  for(Int_t j = 0; j < npoints; j++) {
    avs  += s[j];
    avz  += z[j]; 
    avsz += s[j]*z[j];
  }
    avs  /= (Double_t)npoints;
    avz  /= (Double_t)npoints; 
    avsz /= (Double_t)npoints;

  for(Int_t l = 0; l < npoints; l++) {
    sigmas += (s[l]-avs)*(s[l]-avs);
    sigmaz += (z[l]-avz)*(z[l]-avz);
  }
  sigmas /=(Double_t)npoints;
  sigmaz /=(Double_t)npoints;

  sigmas = sqrt(sigmas);
  sigmaz = sqrt(sigmaz);
  
  corrLin = (avsz-avs*avz)/(sigmas*sigmaz);
  


  delete [] z;
  delete [] x;
  delete [] y;
  delete [] s;
  delete [] ez;
  delete [] ex;
  delete [] ey;
  delete [] es;
  delete f1; delete fitHist;
  return 0;
}


//_______________________________________________________

Double_t AliITSRiemannFit::Fitfunction(Double_t *x, Double_t* par){
  // function used for fit
  return par[0]+(*x)*par[1];

}

