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
 **************************************************************************/
 
/* $Id$ */

// The purpose of this class is to permorm the ITS tracking. The 
// constructor has the task to inizialize some private members. The method 
// DoTracking is written to be called by a macro. It gets the event number,
// the minimum and maximum order number of TPC tracks that are to be tracked
// trough the ITS, and the file where the recpoints are registered. The 
// method Recursivetracking is a recursive function that performs the 
// tracking trough the ITS The method Intersection found the layer, ladder 
// and detector whre the intersection take place and caluclate the 
// cohordinates of this intersection. It returns an integer that is 0 if the 
// intersection has been found successfully. The two mwthods Kalmanfilter 
// and kalmanfiltervert operate the kalmanfilter without and with the vertex
// imposition respectively. The authors thank Mariana Bondila to have help 
// them to resolve some problems.  July-2000

#include <Riostream.h>
#include <TMath.h>
#include <TBranch.h>
#include <TVector.h>
#include <TFile.h>
#include <TRandom.h>
#include <TTree.h>

#include "TParticle.h"
#include "AliRun.h"
#include "AliITSDetTypeRec.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliKalmanTrack.h" 
#include "AliMagF.h"
#include "AliITSTrackV1.h"
#include "AliITSIOTrack.h"
#include "AliITSRad.h"   
#include "AliTPCtracker.h"
#include "AliITSTrackerV1.h"
#include "AliITSPid.h"
#include "AliMC.h"

ClassImp(AliITSTrackerV1)
 //______________________________________________________________________
AliITSTrackerV1::AliITSTrackerV1() {
  //Default constructor
  //fITS = 0;
  fresult = 0;
  fPtref=0.;
  fChi2max=0.;
  //fepsphi=0.;
  //fepsz=0.;
  frecPoints = 0;
  fvettid = 0;
  fflagvert=0;
  frl = 0;
  Int_t ia;
  for(ia=0; ia<6; ia++) {
  fNlad[ia]=0;
  fNdet[ia]=0;
  fAvrad[ia]=0.;
  fDetx[ia]=0.;
  fDetz[ia]=0.; 
  } // end for ia  
  fzmin = 0;
  fzmax = 0;
  fphimin = 0;
  fphimax = 0;
  fphidet = 0;
  fNRecPoints=0;
  fRecCylR=0;
  fRecCylPhi=0;
  fRecCylZ=0;
  fFieldFactor=0;
  fDetTypeRec = 0;
}
//______________________________________________________________________
AliITSTrackerV1::AliITSTrackerV1(Int_t evnumber, Bool_t flag) {
    //Origin   A. Badala' and G.S. Pappalardo:  
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
    // Class constructor. It does some initializations.

  //PH Initialisation taken from the default constructor
  //fITS      = IITTSS;
  fresult = 0;
  fPtref    = 0.;
  fChi2max  =0.; 
  frecPoints = 0;	  	 
  fvettid = 0;
  fflagvert = flag;	 
  frl = 0;
  fzmin = 0;
  fzmax = 0;
  fphimin = 0;
  fphimax = 0;
  fphidet = 0;
  
  Int_t imax = 200,jmax = 450;
  frl       = new AliITSRad(imax,jmax);
  
  //////////  gets information on geometry /////////////////////////////
  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  AliITSLoader* loader = (AliITSLoader*)runLoader->GetLoader("ITSLoader");
  if (!loader) {
    Error("AliITSTrackerV1", "ITS loader not found");
    return;
  }
  fDetTypeRec = new AliITSDetTypeRec(loader);
  AliITSgeom* g1 = loader->GetITSgeom();

  Int_t ll=1, dd=1;
  TVector det(9);
  
  Int_t ia;
  for(ia=0; ia<6; ia++) {
    fNlad[ia]=g1->GetNladders(ia+1);
    fNdet[ia]=g1->GetNdetectors(ia+1);
    //cout<<fNlad[i]<<" "<<fNdet[i]<<"\n"; 
  } // end for ia
  
    //cout<<" mean radius = ";
  Int_t ib;
  for(ib=0; ib<6; ib++) {  
    g1->GetCenterThetaPhi(ib+1,ll,dd,det);
    Double_t r1=TMath::Sqrt(det(0)*det(0)+det(1)*det(1));
    g1->GetCenterThetaPhi(ib+1,ll,dd+1,det);
    Double_t r2=TMath::Sqrt(det(0)*det(0)+det(1)*det(1));
    fAvrad[ib]=(r1+r2)/2.;
    //cout<<fAvrad[ib]<<" ";
  } // end for ib
    //cout<<"\n"; getchar();
  
  fDetx[0] = ((AliITSgeomSPD*)(g1->GetShape(1, ll, dd)))->GetDx();
  fDetz[0] = ((AliITSgeomSPD*)(g1->GetShape(1, ll, dd)))->GetDz();
  
  fDetx[1] = ((AliITSgeomSPD*)(g1->GetShape(2, ll, dd)))->GetDx();
  fDetz[1] = ((AliITSgeomSPD*)(g1->GetShape(2, ll, dd)))->GetDz();

    fDetx[2] = ((AliITSgeomSDD*)(g1->GetShape(3, ll, dd)))->GetDx();
    fDetz[2] = ((AliITSgeomSDD*)(g1->GetShape(3, ll, dd)))->GetDz();

    fDetx[3] = ((AliITSgeomSDD*)(g1->GetShape(4, ll, dd)))->GetDx();
    fDetz[3] = ((AliITSgeomSDD*)(g1->GetShape(4, ll, dd)))->GetDz();

    fDetx[4] = ((AliITSgeomSSD*)(g1->GetShape(5, ll, dd)))->GetDx();
    fDetz[4] = ((AliITSgeomSSD*)(g1->GetShape(5, ll, dd)))->GetDz();

    fDetx[5] = ((AliITSgeomSSD*)(g1->GetShape(6, ll, dd)))->GetDx();
    fDetz[5] = ((AliITSgeomSSD*)(g1->GetShape(6, ll, dd)))->GetDz();
    //cout<<"    Detx     Detz\n";
    //for(Int_t la=0; la<6; la++) cout<<"    "<<fDetx[la]<<"     "<<
    //                                 fDetz[la]<<endl;
    //getchar();
	 
    // allocate memory and define matrices fzmin, fzmax, fphimin and fphimax //
    Double_t epsz=1.2;
    Double_t epszdrift=0.05;

    fzmin = new Double_t*[6]; fzmax = new Double_t*[6];
    Int_t im1, im2, im2max;
    for(im1=0; im1<6; im1++) {
	im2max=fNdet[im1];
	fzmin[im1] = new Double_t[im2max]; fzmax[im1] = new Double_t[im2max];
    } // end for im1

    for(im1=0; im1<6; im1++) {
	im2max=fNdet[im1];
	for(im2=0; im2<im2max; im2++) {
	    g1->GetCenterThetaPhi(im1+1,1,im2+1,det);
	    if(im2!=0) fzmin[im1][im2]=det(2)-fDetz[im1];
	    else   
		fzmin[im1][im2]=det(2)-(fDetz[im1])*epsz;
	    if(im2!=(im2max-1)) fzmax[im1][im2]=det(2)+fDetz[im1];
	    else
		fzmax[im1][im2]=det(2)+fDetz[im1]*epsz;
	    if(im1==2 || im1==3) {
		fzmin[im1][im2]-=epszdrift;
		fzmax[im1][im2]+=epszdrift;
	    } // end if im1==2 || im1==3
	} // end for im2
    } // end for im1

    fphimin = new Double_t*[6]; fphimax = new Double_t*[6];
    for(im1=0;im1<6;im1++) {
	im2max=fNlad[im1];
	fphimin[im1] = new Double_t[im2max]; 
	fphimax[im1] = new Double_t[im2max];
    } // end for im1
  
    fphidet = new Double_t*[6];
    for(im1=0; im1<6; im1++) {
	im2max=fNlad[im1];
	fphidet[im1] = new Double_t[im2max];
    } // end for im1

	Double_t global[3],local[3];
    Double_t pigre=TMath::Pi();
    Double_t xmin,ymin,xmax,ymax;

    for(im1=0; im1<6; im1++) {
	im2max=fNlad[im1];
	for(im2=0; im2<im2max; im2++) {
	    Int_t idet=2;
	    g1->GetCenterThetaPhi(im1+1,im2+1,idet,det);
	    fphidet[im1][im2] = TMath::ATan2(Double_t(det(1)),
					     Double_t(det(0))); 
	    if(fphidet[im1][im2]<0.) fphidet[im1][im2]+=2.*pigre;  
	    local[1]=local[2]=0.;  
	    local[0]= -(fDetx[im1]);    
	    if(im1==0) local[0]= (fDetx[im1]); //to take into account 
	                                       // different reference system
	    g1->LtoG(im1+1,im2+1,idet,local,global);
	    xmax=global[0]; ymax=global[1];
	    local[0]= (fDetx[im1]);   
	    if(im1==0) local[0]= -(fDetx[im1]);//take into account different 
	                                       // reference system
	    g1->LtoG(im1+1,im2+1,idet,local,global);
	    xmin=global[0]; ymin=global[1];
	    fphimin[im1][im2]= TMath::ATan2(ymin,xmin);
	    if(fphimin[im1][im2]<0.) fphimin[im1][im2]+=2.*pigre; 
	    fphimax[im1][im2]= TMath::ATan2(ymax,xmax);
	    if(fphimax[im1][im2]<0.) fphimax[im1][im2]+=2.*pigre;
	} // end for im2
    } // end for im1
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// allocate memory and define vector fNRecPoints and matrices fRecCylR, fRecCylPhi, fRecCylZ /////////////
	gAlice->GetEvent(evnumber);
  Int_t numOfModules = g1->GetIndexMax();
  fRecCylR = new Double_t *[numOfModules];
  fRecCylPhi = new Double_t *[numOfModules]; 
  fRecCylZ = new Double_t *[numOfModules];
  AliITSRecPoint *recp;
  fNRecPoints = new Int_t[numOfModules];
   
		 for(Int_t module=0; module<numOfModules; module++) {				
		  fDetTypeRec->ResetRecPoints();		     
		  gAlice->TreeR()->GetEvent(module);		  
		  frecPoints=fDetTypeRec->RecPoints();
		  Int_t nRecPoints=fNRecPoints[module]=frecPoints->GetEntries();
		  fRecCylR[module] = new Double_t[nRecPoints];
		  fRecCylPhi[module] = new Double_t[nRecPoints];
		  fRecCylZ[module] = new  Double_t[nRecPoints];		  
		  Int_t ind;
		  for(ind=0; ind<fNRecPoints[module]; ind++) {	  
		    recp=(AliITSRecPoint*)frecPoints->UncheckedAt(ind);			 			 
			// Float_t global[3], local[3];
			 Double_t global[3], local[3];
	       local[0]=recp->GetX();
	       local[1]=0.;
	       local[2]= recp->GetZ();			 		
			 g1->LtoG(module,local,global);
			                                                                 
			 Double_t r = TMath::Sqrt(global[0]*global[0]+global[1]*global[1]);                     // r hit
			 Double_t phi = TMath::ATan2(global[1],global[0]); if(phi<0.) phi+=2.*TMath::Pi();      // phi hit			
          	 Double_t z = global[2];                                                                // z hit
			 																								             
			 fRecCylR[module][ind]=r;
			 fRecCylPhi[module][ind]=phi;
			 fRecCylZ[module][ind]=z;			 
		  } 		
		}	 
	 //}  
  //}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 

    ////////// gets magnetic field factor //////////////////////////////
    AliMagF * fieldPointer = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
    fFieldFactor =  fieldPointer ? fieldPointer->SolenoidField()/10/.2 : 0;
   // cout<< " field factor = "<<fFieldFactor<<"\n"; getchar();
    delete rl;
    
}
//______________________________________________________________________
AliITSTrackerV1::AliITSTrackerV1(const AliITSTrackerV1 &cobj) : TObject(cobj) {
    // Origin  A. Badala' and G.S. Pappalardo:
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
    // copy constructor
	 
  //	 *fITS = *cobj.fITS;
    *fresult = *cobj.fresult;
    fPtref = cobj.fPtref;
    fChi2max = cobj.fChi2max;    
    **fvettid = **cobj.fvettid;
    fflagvert = cobj.fflagvert;
    Int_t imax=200,jmax=450;
    frl = new AliITSRad(imax,jmax);	 
    *frl = *cobj.frl;
    fFieldFactor = cobj.fFieldFactor;
    Int_t i,im1,im2,im2max;
    for(i=0; i<6; i++) {
	fNlad[i] = cobj.fNlad[i];
	fNdet[i] = cobj.fNdet[i]; 
	fAvrad[i] = cobj.fAvrad[i];
	fDetx[i] = cobj.fDetx[i];
	fDetz[i] = cobj.fDetz[i];
    } // end or i
    fzmin = new Double_t*[6]; fzmax = new Double_t*[6];
    for(im1=0; im1<6; im1++) {
	im2max=fNdet[im1];
	fzmin[im1] = new Double_t[im2max];
	fzmax[im1] = new Double_t[im2max];
    } // end for im1
    fphimin = new Double_t*[6]; fphimax = new Double_t*[6];
    for(im1=0;im1<6;im1++) {
	im2max=fNlad[im1];
	fphimin[im1] = new Double_t[im2max];
	fphimax[im1] = new Double_t[im2max];
    } // end for im1
  
    fphidet = new Double_t*[6];
    for(im1=0; im1<6; im1++) {
	im2max=fNlad[im1];
	fphidet[im1] = new Double_t[im2max];
    } // end for im1
    for(im1=0; im1<6; im1++) {
	im2max=fNdet[im1];
	for(im2=0; im2<im2max; im2++) {
	    fzmin[im1][im2]=cobj.fzmin[im1][im2];
	    fzmax[im1][im2]=cobj.fzmax[im1][im2];
	} // end for im2
    } // end for im1
    for(im1=0; im1<6; im1++) {
	im2max=fNlad[im1];
	for(im2=0; im2<im2max; im2++) {
	    fphimin[im1][im2]=cobj.fphimin[im1][im2];
	    fphimax[im1][im2]=cobj.fphimax[im1][im2];
	    fphidet[im1][im2]=cobj.fphidet[im1][im2];  
	} // end for im2
    } // end for im2


    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    rl->CdGAFile();
    AliITSgeom* g1 = (AliITSgeom*)gDirectory->Get("AliITSgeom");

   Int_t numOfModules = g1->GetIndexMax();
	/*
  fRecCylR = new Float_t *[numOfModules];
  fRecCylPhi = new Float_t *[numOfModules]; 
  fRecCylZ = new Float_t *[numOfModules];
  */
  fRecCylR = new Double_t *[numOfModules];
  fRecCylPhi = new Double_t *[numOfModules]; 
  fRecCylZ = new Double_t *[numOfModules];  
  fNRecPoints = new Int_t[numOfModules];	
		for(Int_t module=0; module<numOfModules; module++) {		
		  Int_t nRecPoints=fNRecPoints[module]=cobj.fNRecPoints[module];
		  /*
		  fRecCylR[module] = new Float_t[nRecPoints];
		  fRecCylPhi[module] = new Float_t[nRecPoints];
		  fRecCylZ[module] = new Float_t[nRecPoints];
		  */
		  fRecCylR[module] = new Double_t[nRecPoints];
		  fRecCylPhi[module] = new Double_t[nRecPoints];
		  fRecCylZ[module] = new Double_t[nRecPoints];		  
		  Int_t ind;	
		  for(ind=0; ind<nRecPoints; ind++) {	    
			 fRecCylR[module][ind]=cobj.fRecCylR[module][ind];
			 fRecCylPhi[module][ind]=cobj.fRecCylPhi[module][ind];
			 fRecCylZ[module][ind]=cobj.fRecCylZ[module][ind];			 
		  } 		
		}	 
		delete rl;
 
}
/*
//______________________________________________________________________
void AliITSTrackerV1::DelMatrix(Int_t numOfModules) {
  // cleanup of some data members
  for(Int_t mod=0; mod<numOfModules; mod++) {
    delete fRecCylR[mod];
	 delete fRecCylPhi[mod];
	 delete fRecCylZ[mod];
  }
    delete fRecCylR;
	 delete fRecCylPhi;
	 delete fRecCylZ;
}
*/
//______________________________________________________________________
AliITSTrackerV1::~AliITSTrackerV1(){
    // Origin  A. Badala' and G.S. Pappalardo:
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
    // class destructor	 
  delete frl;
  delete fNRecPoints;
  for(Int_t i=0; i<6; i++) {
    delete fzmin[i];
	 delete fzmax[i];
	 delete fphimin[i];
	 delete fphimax[i];
	 delete fphidet[i];
  }

  delete fzmin;
  delete fzmax;
  delete fphimin;
  delete fphimax;
  delete fphidet;
	 
}
//______________________________________________________________________
AliITSTrackerV1 &AliITSTrackerV1::operator=(const AliITSTrackerV1 &obj) {
    // Origin  A. Badala' and G.S. Pappalardo:
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
    // assignement operator

  //	 *fITS = *obj.fITS;
    *fresult = *obj.fresult;
    fPtref = obj.fPtref;
    fChi2max = obj.fChi2max;      
    **fvettid = **obj.fvettid;
    fflagvert = obj.fflagvert;
    Int_t imax=200,jmax=450;
    frl = new AliITSRad(imax,jmax);	 
    *frl = *obj.frl;
    fFieldFactor = obj.fFieldFactor;
    Int_t i;
    for(i=0; i<6; i++) {
	fNlad[i] = obj.fNlad[i];
	fNdet[i] = obj.fNdet[i]; 
	fAvrad[i] = obj.fAvrad[i];
	fDetx[i] = obj.fDetx[i];
	fDetz[i] = obj.fDetz[i];
    } // end for i
    fzmin = new Double_t*[6];
    fzmax = new Double_t*[6];
    Int_t im1, im2, im2max;
    for(im1=0; im1<6; im1++) {
	im2max=fNdet[im1];
	fzmin[im1] = new Double_t[im2max]; fzmax[im1] = new Double_t[im2max];
    } // end for im1
    fphimin = new Double_t*[6]; fphimax = new Double_t*[6];
    for(im1=0;im1<6;im1++) {
	im2max=fNlad[im1];
	fphimin[im1] = new Double_t[im2max];
	fphimax[im1] = new Double_t[im2max];
    } // end for im1

    fphidet = new Double_t*[6];
    for(im1=0; im1<6; im1++) {
	im2max=fNlad[im1];
	fphidet[im1] = new Double_t[im2max];
    } // end for im1
    for(im1=0; im1<6; im1++) {
	im2max=fNdet[im1];
	for(im2=0; im2<im2max; im2++) {
	    fzmin[im1][im2]=obj.fzmin[im1][im2];
	    fzmax[im1][im2]=obj.fzmax[im1][im2];
	} // end for im2
    } // end for im1
    for(im1=0; im1<6; im1++) {
	im2max=fNlad[im1];
	for(im2=0; im2<im2max; im2++) {
	    fphimin[im1][im2]=obj.fphimin[im1][im2];
	    fphimax[im1][im2]=obj.fphimax[im1][im2];
	    fphidet[im1][im2]=obj.fphidet[im1][im2];  
	} // end for im2
    } // end for im1

    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    rl->CdGAFile();
    AliITSgeom* g1 = (AliITSgeom*)gDirectory->Get("AliITSgeom");
    //	AliITSgeom *g1 = fITS->GetITSgeom();  
   Int_t numOfModules = g1->GetIndexMax();
  fRecCylR = new Double_t *[numOfModules];
  fRecCylPhi = new Double_t *[numOfModules]; 
  fRecCylZ = new Double_t *[numOfModules];  
  fNRecPoints = new Int_t[numOfModules];  
	  for(Int_t module=0; module<numOfModules; module++) {		  
		  Int_t nRecPoints=fNRecPoints[module]=obj.fNRecPoints[module];
		  fRecCylR[module] = new Double_t[nRecPoints];
		  fRecCylPhi[module] = new Double_t[nRecPoints];
		  fRecCylZ[module] = new Double_t[nRecPoints];		  
		  Int_t ind;
		  for(ind=0; ind<nRecPoints; ind++) {	  
			 fRecCylR[module][ind]=obj.fRecCylR[module][ind];
			 fRecCylPhi[module][ind]=obj.fRecCylPhi[module][ind];
			 fRecCylZ[module][ind]=obj.fRecCylZ[module][ind];			 
		  } 		
		}	 
         
	 
	  delete rl;
    return *this;
}
//______________________________________________________________________
void AliITSTrackerV1::DoTracking(Int_t evNumber,Int_t minTr,Int_t maxTr,
				 TFile *file, Bool_t realmass) {
    // Origin   A. Badala' and G.S. Pappalardo:
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
    // The method needs the event number, the minimum and maximum order
    // number of TPC tracks that 
    // are to be tracked trough the ITS, and the file where the recpoints
    // are registered.
    // The method can be called by a macro. It preforms the tracking for
    // all good TPC tracks

    printf("begin DoTracking - file %p\n",(void*)file);

    gAlice->GetEvent(evNumber);  //modificato per gestire hbt
 
    TFile *cf=TFile::Open("AliTPCclusters.root");  
    AliTPCParam *digp= (AliTPCParam*)cf->Get("75x40_100x60_150x60");
    if (!digp) { cerr<<"TPC parameters have not been found !\n"; getchar();}

    cf->cd();
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
    TString foldname(itsl->GetEventFolder()->GetName());
    //TString foldname(fITS->GetLoader()->GetEventFolder()->GetName());
    delete rl;   
    printf("This method is not converted to the NewIO !\n");  //I.B.
    return; //I.B.
    AliTPCtracker *tracker = new AliTPCtracker(digp);            //I.B.
    //PH    tracker->SetEventNumber(evNumber);                           //I.B.

    // Load clusters
    printf("This method is not converted to the NewIO !\n");  //I.B.
    return; //I.B.
    tracker->LoadClusters(0);       //I.B.

    // Load tracks
    TFile *tf=TFile::Open("AliTPCtracksSorted.root");  
    if (!tf->IsOpen()) {
	cerr<<"Can't open AliTPCtracksSorted.root !\n";
	return ;
    } // end if
    TObjArray tracks(200000);
    char tname[100];                 
    sprintf(tname,"TreeT_TPC_%d",evNumber);

    TTree *tracktree=(TTree*)tf->Get(tname);   
    if (!tracktree) {cerr<<"Can't get a tree with TPC tracks !\n";}   
    TBranch *tbranch=tracktree->GetBranch("tracks");
    Int_t nentr=(Int_t)tracktree->GetEntries();
    Int_t kk;

    AliITSRecPoint *recp;    // oggi
	AliTPCtrack *ioTrackTPC=0;    
    for (kk=0; kk<nentr; kk++) {
	ioTrackTPC=new AliTPCtrack; 
	tbranch->SetAddress(&ioTrackTPC);
	tracktree->GetEvent(kk);    
	tracker->CookLabel(ioTrackTPC,0.1);       
	tracks.AddLast(ioTrackTPC);         
    } // end for kk
    delete tracker;      
    tf->Close();

    Int_t nt = tracks.GetEntriesFast();
    cerr<<"Number of found tracks "<<nt<<endl;

    TVector vec(5);
    TTree *tr=gAlice->TreeR();
    Int_t nent=(Int_t)tr->GetEntries();  
    frecPoints = fDetTypeRec->RecPoints();

    Int_t numbpoints;
    Int_t totalpoints=0;
    Int_t *np = new Int_t[nent];
    fvettid = new Int_t* [nent];
    Int_t mod;

    for (mod=0; mod<nent; mod++) {
	fvettid[mod]=0;
	fDetTypeRec->ResetRecPoints();  
	gAlice->TreeR()->GetEvent(mod); 
	numbpoints = frecPoints->GetEntries();
	totalpoints+=numbpoints;
	np[mod] = numbpoints;
	//cout<<" mod = "<<mod<<"   numbpoints = "<<numbpoints<<"\n";getchar();
	fvettid[mod] = new Int_t[numbpoints];
	Int_t ii;
	for (ii=0;ii<numbpoints; ii++) *(fvettid[mod]+ii)=0;
    } // end for mod

    AliTPCtrack *track=0;

    if(minTr < 0) {minTr = 0; maxTr = nt-1;}   

    TVector vgeant(3);

    TTree tracktree1("TreeT","Tree with ITS tracks");
    AliITSIOTrack *ioTrack=0;
	AliITSPid *pid=new AliITSPid(1000);  // oggi

    tracktree1.Branch("ITStracks","AliITSIOTrack",&ioTrack,32000,0);
  
    TDatabasePDG * db = new TDatabasePDG;   
  
    Int_t j;       
    for (j=minTr; j<=maxTr; j++) {     
	track=(AliTPCtrack*)tracks.UncheckedAt(j);
	if (!track) continue;
	
	///   mass definition ////////////////////////
    Double_t mass=0.13956995;
	 Int_t pcode=211;  // a pion by default
	 
	 if(realmass) {
    	if(TMath::Abs(pcode)<20443) mass=db->GetParticle(pcode)->Mass();
	 }
	 else {
		 mass = track->GetMass();
//		 cout << "Mass = " << mass << endl;
	 }
	 
	
	
	 //   new propagation to the end of TPC
    Double_t xk=80.;
   // track->PropagateTo(xk,0.,0.); //Ne if it's still there   //attenzione funziona solo se modifica in TPC
	// Double_t xk=77.415;	 
	 track->PropagateTo(xk, 28.94, 1.204e-3);
    xk-=0.005;
    track->PropagateTo(xk, 44.77,1.71); //Tedlar	 
    xk-=0.02;
    track->PropagateTo(xk, 44.86, 1.45);   //Kevlar
    xk-=2.0;
    track->PropagateTo(xk, 41.28, 0.029);//Nomex
    xk-=0.02;
    track->PropagateTo(xk, 44.86, 1.45);   //Kevlar
    xk-=0.005;
    track->PropagateTo(xk, 44.77, 1.71); //Tedlar

    xk=61.;
   // track->PropagateTo(xk,0.,0.); //C02
	 track->PropagateTo(xk,36.2,1.98e-3); //C02	   //attenzione funziona solo se modifica in TPC

    xk -=0.005;
    track->PropagateTo(xk, 24.01, 2.7);    //Al    
    xk -=0.005;
    track->PropagateTo(xk, 44.77, 1.71);  //Tedlar
    xk -=0.02;
    track->PropagateTo(xk, 44.86, 1.45);    //Kevlar
    xk -=0.5;
    track->PropagateTo(xk, 41.28, 0.029);  //Nomex    
    xk -=0.02;
    track->PropagateTo(xk, 44.86, 1.45);    //Kevlar
    xk -=0.005;
    track->PropagateTo(xk, 44.77, 1.71);  //Tedlar
    xk -=0.005;
    track->PropagateTo(xk, 24.01, 2.7);    //Al 
    
	////////////////////////////////////////////////////////////////////////////////////////////////////////  	
	//AliITSTrackV1 trackITS(*track);
	AliITSTrackV1 trackITS(*track, fFieldFactor);
	//cout<<" fFieldFactor = "<<fFieldFactor<<"\n";      
	trackITS.PutMass(mass);	  //new to add mass to track
	if(fresult){ delete fresult; fresult=0;}   	 
	fresult = new AliITSTrackV1(trackITS);	 

	AliITSTrackV1 primaryTrack(trackITS);
	vgeant=(*fresult).GetVertex();

	// Definition of dv and zv for vertex constraint	
	Double_t sigmaDv=0.0050;  Double_t sigmaZv=0.010;	
	//Double_t sigmaDv=0.0015;  Double_t sigmaZv=0.0015;
	Double_t uniform= gRandom->Uniform();
	Double_t signdv;
	if(uniform<=0.5) signdv=-1.;
	else
	    signdv=1.;

	Double_t vr=TMath::Sqrt(vgeant(0)*vgeant(0)+ vgeant(1)*vgeant(1));
	Double_t dv=gRandom->Gaus(signdv*vr,(Float_t)sigmaDv); 
	Double_t zv=gRandom->Gaus(vgeant(2),(Float_t)sigmaZv);		
	//cout<<" Dv e Zv = "<<dv<<" "<<zv<<"\n";
	trackITS.SetDv(dv);
	trackITS.SetZv(zv);
	trackITS.SetsigmaDv(sigmaDv);
	trackITS.SetsigmaZv(sigmaZv); 
	(*fresult).SetDv(dv);
	(*fresult).SetZv(zv);
	(*fresult).SetsigmaDv(sigmaDv);
	(*fresult).SetsigmaZv(sigmaZv);
	primaryTrack.SetDv(dv);
	primaryTrack.SetZv(zv);
	primaryTrack.SetsigmaDv(sigmaDv);
	primaryTrack.SetsigmaZv(sigmaZv);
	primaryTrack.PrimaryTrack();
	TVector  d2=primaryTrack.Getd2();
	TVector  tgl2=primaryTrack.Gettgl2();
	TVector  dtgl=primaryTrack.Getdtgl();
	trackITS.Setd2(d2); trackITS.Settgl2(tgl2);
	trackITS.Setdtgl(dtgl); 
	(*fresult).Setd2(d2); (*fresult).Settgl2(tgl2);
	(*fresult).Setdtgl(dtgl);
	/*
	  trackITS.SetVertex(vertex); trackITS.SetErrorVertex(ervertex);
	  (*result).SetVertex(vertex);   (*result).SetErrorVertex(ervertex);   
	*/
	TList *list= new TList();

	list->AddLast(&trackITS);
  
	fPtref=TMath::Abs( (trackITS).GetPt() );
	//cout<<" fPtref = " <<fPtref<<"\n";
	if(fPtref>1.0) fChi2max=40.;         
	if(fPtref<=1.0) fChi2max=20.;
	if(fPtref<0.4 ) fChi2max=100.;
	if(fPtref<0.2 ) fChi2max=40.;		  
	// if(fPtref<0.4 ) fChi2max=30.;				 
	// if(fPtref<0.2 ) fChi2max=20.;
	//if(fPtref<0.2 ) fChi2max=10.;
	//if(fPtref<0.1 ) fChi2max=5.;
	//cout << "\n Pt = " << fPtref <<"\n";  //stampa
	RecursiveTracking(list);   
	list->Delete();
	delete list;

	Int_t itot=-1;
	TVector vecTotLabRef(18);
	Int_t lay, k;
	for(lay=5; lay>=0; lay--) {
	    TVector vecLabRef(3); 
	    vecLabRef=(*fresult).GetLabTrack(lay);
	    Float_t clustZ=(*fresult).GetZclusterTrack( lay);   
	    for(k=0; k<3; k++){  
		Int_t lpp=(Int_t)vecLabRef(k);
		if(lpp>=0) {
		    TParticle *p=(TParticle*) gAlice->GetMCApp()->Particle(lpp);
		    Int_t pcode=p->GetPdgCode();
		    if(pcode==11) vecLabRef(k)=p->GetFirstMother();
		} // end if
		itot++; vecTotLabRef(itot)=vecLabRef(k);
		if(vecLabRef(k)==0. && clustZ == -1.) vecTotLabRef(itot) =-3.;
	    } // end for k
	} // end for lay
	Long_t labref;
	Int_t freq;  
	(*fresult).Search(vecTotLabRef, labref, freq);

	//if(freq < 6) labref=-labref;        // cinque - sei
	if(freq < 5) labref=-labref;        // cinque - sei	
	(*fresult).SetLabel(labref);

	// cout<<" progressive track number = "<<j<<"\r";
	// cout<<j<<"\r";
	Int_t numOfCluster=(*fresult).GetNumClust();  
	//cout<<" progressive track number = "<<j<<"\n";    // stampa
	Long_t labITS=(*fresult).GetLabel();
	//cout << " ITS track label = " << labITS << "\n"; 	// stampa	    
	Int_t lab=track->GetLabel();		    
	//cout << " TPC track label = " << lab <<"\n";      // stampa
	//propagation to vertex

	Double_t rbeam=3.;
	if((*fresult).DoNotCross(rbeam)) continue;  //no intersection with beampipe	
	(*fresult).Propagation(rbeam);
	Double_t c00,c10,c11,c20,c21,c22,c30,c31,c32,c33,c40,c41,c42,c43,c44;
	(*fresult).GetCElements(c00,
				c10,c11,
				c20,c21,c22,
				c30,c31,c32,c33,
				c40,c41,c42,c43,c44);
	 	 
	Double_t pt=TMath::Abs((*fresult).GetPt());
	Double_t dr=(*fresult).GetD();
	Double_t z=(*fresult).GetZ();
	Double_t tgl=(*fresult).GetTgl();
	Double_t c=(*fresult).GetC();
	Double_t cy=c/2.;
	Double_t dz=z-(tgl/cy)*TMath::ASin((*fresult).Arga(rbeam));
	dz-=vgeant(2);
	// cout<<" dr e dz alla fine = "<<dr<<" "<<dz<<"\n"; getchar();
	Double_t phi=(*fresult).Getphi();
	Double_t phivertex = phi - TMath::ASin((*fresult).ArgA(rbeam));
	Double_t duepi=2.*TMath::Pi();	 
	if(phivertex>duepi) phivertex-=duepi;
	if(phivertex<0.) phivertex+=duepi;
	/////////////////////////////////////////////////////////////
	Int_t idmodule,idpoint;
	if(numOfCluster >=5)  {            // cinque - sei
	    //if(numOfCluster ==6)  {            // cinque - sei
	    AliITSIOTrack outTrack;
	    ioTrack=&outTrack;
	    ioTrack->SetStatePhi(phi);
	    ioTrack->SetStateZ(z);
	    ioTrack->SetStateD(dr);
	    ioTrack->SetStateTgl(tgl);
	    ioTrack->SetStateC(c);
	    Double_t radius=(*fresult).Getrtrack();
	    ioTrack->SetRadius(radius);
	    Int_t charge;
	    if(c>0.) charge=-1;  else charge=1;
	    ioTrack->SetCharge(charge);
		Double_t trackmass=(*fresult).GetMass();         // oggi
		ioTrack->SetMass(trackmass);                     // oggi
	    ioTrack->SetCovMatrix(c00,
				  c10,c11,
				  c20,c21,c22,
				  c30,c31,c32,c33,
				  c40,c41,c42,c43,c44);
	    Double_t px=pt*TMath::Cos(phivertex);
	    Double_t py=pt*TMath::Sin(phivertex);
	    Double_t pz=pt*tgl;
	    Double_t xtrack=dr*TMath::Sin(phivertex);
	    Double_t ytrack=dr*TMath::Cos(phivertex);
	    Double_t ztrack=dz+vgeant(2);
	    ioTrack->SetPx(px);
	    ioTrack->SetPy(py);
	    ioTrack->SetPz(pz);
	    ioTrack->SetX(xtrack);
	    ioTrack->SetY(ytrack);
	    ioTrack->SetZ(ztrack);
	    ioTrack->SetLabel(labITS);
	    ioTrack->SetTPCLabel(lab);
		 ioTrack->SetDz(dz);
	    Int_t il;		
	    /*
		for(il=0;il<6; il++){
		ioTrack->SetIdPoint(il,(*fresult).GetIdPoint(il));
		ioTrack->SetIdModule(il,(*fresult).GetIdModule(il));
	    } // end for il
	    */
	    //tracktree1.Fill();
		Float_t q[4]={-1.,-1.,-1.,-1.};
	    Float_t  globaldedx=0.;	   
	    for (il=0;il<6;il++) {
		idpoint=(*fresult).GetIdPoint(il);
		idmodule=(*fresult).GetIdModule(il);
		if(idmodule>0.) *(fvettid[idmodule]+idpoint)=1;
		                                               
		ioTrack->SetIdPoint(il,idpoint);
		ioTrack->SetIdModule(il,idmodule);
		////  for q definition
  	        if(il>1){
	          if(idmodule>0.){			
	            fDetTypeRec->ResetRecPoints();
	            gAlice->TreeR()->GetEvent(idmodule);
	            recp=(AliITSRecPoint*)frecPoints->UncheckedAt(idpoint);
	            q[il-2]=recp->GetQ()*(*fresult).Getfcor(il-2);
	          }
	        }				
	    } // end for il	 
		q[0]/=280.; q[1]/=280.;
	    q[2]/=38.; q[3]/=38.;

 // cout<<" q prima = "<<q[0]<<" "<<q[1]<<" "<<q[2]<<" "<<q[3]<<"\n"; getchar(); 
     
     Int_t swap;
  do{
    swap=0;   
    for (il=0; il<3; il++) {
      if (q[il]<=q[il+1]) continue;
      Float_t tmp=q[il];
      q[il]=q[il+1]; q[il+1]=tmp;
      swap++;
    }
  } while(swap); 

 
  // cout<<" q dopo = "<<q[0]<<" "<<q[1]<<" "<<q[2]<<" "<<q[3]<<"\n"; getchar();
      
    if(q[0]<0.) {
      q[0]=q[1];
      q[1]=q[2];
      q[2]=q[3];
      q[3]=-1.; 	    
    } 
 
  // cout<<" q  dopo if = "<<q[0]<<" "<<q[1]<<" "<<q[2]<<" "<<q[3]<<"\n"; getchar(); 
     
    globaldedx=(q[0]+q[1])/2.;
    
   // if(q[3]> 0.) globaldedx=(q[0]+q[1]+q[2]+q[3])/4.;
   //      else    globaldedx=(q[0]+q[1]+q[2])/3.; 
   
    ioTrack->SetdEdx(globaldedx);
    ioTrack->SetPid(pid->GetPcode(ioTrack));
	    
	    tracktree1.Fill();	       
	} // end if on numOfCluster
	//gObjectTable->Print();    // stampa memoria     
    }  //  end for (int j=minTr; j<=maxTr; j++)
    delete db;    	 
    static Bool_t first=kTRUE;
    static TFile *tfile;
    if(first) {
	tfile=new TFile("itstracks.root","RECREATE");
	//cout<<"I have opened itstracks.root file "<<endl;
    } // end if	    
    first=kFALSE;
    tfile->cd();
    tfile->ls();
    char hname[30];
    sprintf(hname,"TreeT%d",evNumber);
	cout << "Number of saved ITS tracks " << tracktree1.GetEntries() << endl;
    tracktree1.Write(hname);
  
    TTree *fAli=gAlice->TreeK();
    TFile *fileAli=0;
    if (fAli) fileAli =fAli->GetCurrentFile();
    fileAli->cd();
    ////////////////////////////////////////////////////////////////////

    printf("delete vectors\n");
    if(np) delete [] np;
    if(fvettid) delete [] fvettid;
    if(fresult) {delete fresult; fresult=0;}
}
//______________________________________________________________________
void AliITSTrackerV1::RecursiveTracking(TList *trackITSlist) {
    //   This function perform the recursive tracking in ITS detectors
    // reference is a pointer to the final best track
    // Origin  A. Badala' and G.S. Pappalardo:
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
    // The authors thank Mariana Bondila to have help them to resolve some
    // problems.  July-2000

    //Rlayer[0]=4.; Rlayer[1]=7.;  Rlayer[2]=14.9;
    // Rlayer[3]=23.8;  Rlayer[4]=39.1;  Rlayer[5]=43.6; //vecchio

    ////////////////////// 
    Float_t sigmaphil[6], sigmazl[6];
    sigmaphil[0]=1.44e-6/(fAvrad[0]*fAvrad[0]);
    sigmaphil[1]=1.44e-6/(fAvrad[1]*fAvrad[1]);
    sigmaphil[2]=1.444e-5/(fAvrad[2]*fAvrad[2]);
    sigmaphil[3]=1.444e-5/(fAvrad[3]*fAvrad[3]);
    sigmaphil[4]=4e-6/(fAvrad[4]*fAvrad[4]);
    sigmaphil[5]=4e-6/(fAvrad[5]*fAvrad[5]);
    sigmazl[0]=1e-2;
    sigmazl[1]=1e-2;
    sigmazl[2]=7.84e-4;
    sigmazl[3]=7.84e-4;
    sigmazl[4]=0.6889;
    sigmazl[5]=0.6889;	
    ///////////////////////////////////////////////////////////
    Int_t index; 
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    rl->CdGAFile();
    AliITSgeom* g1 = (AliITSgeom*)gDirectory->Get("AliITSgeom");

    AliITSRecPoint *recp;        
    for(index =0; index<trackITSlist->GetSize(); index++) {
	AliITSTrackV1 *trackITS = (AliITSTrackV1 *) trackITSlist->At(index);
	if((*trackITS).GetLayer()==7) fresult->SetChi2(10.223e140);
	// cout <<" Layer inizio = "<<(*trackITS).GetLayer()<<"\n";
	//  cout<<"fvtrack =" <<"\n";
	//  cout << (*trackITS)(0) << " "<<(*trackITS)(1)<<" "
	//       <<(*trackITS)(2)<<" "<<(*trackITS)(3)<<" "
	//       <<(*trackITS)(4)<<"\n";
	//  cout<< " rtrack = "<<(*trackITS).Getrtrack()<<"\n";
	//  cout<< " Pt = "<<(*trackITS).GetPt()<<"\n";
	//  getchar();    
	Double_t chi2Now, chi2Ref;
	Float_t numClustRef = fresult->GetNumClust();		 
	if((*trackITS).GetLayer()==1 ) {
	    chi2Now = trackITS->GetChi2();
	    Float_t numClustNow = trackITS->GetNumClust();
	    if(trackITS->GetNumClust()) 
		chi2Now /= (Double_t)trackITS->GetNumClust();
	    chi2Ref = fresult->GetChi2();
	    if(fresult->GetNumClust()) 
		chi2Ref /= (Double_t)fresult->GetNumClust();
	    //cout<<" chi2Now and chi2Ref = "<<chi2Now<<" "<<chi2Ref<<"\n";
	    if( numClustNow > numClustRef ) {*fresult = *trackITS;} 
	    if((numClustNow == numClustRef )&& 
	       (chi2Now < chi2Ref))  {
		*fresult = *trackITS;
	    } // end if
	    continue;	
	} // end if

	if(trackITS->Getfnoclust()>=2)  continue;      
	Float_t numClustNow = trackITS->GetNumClust();
	if(numClustNow) { 
	    chi2Now = trackITS->GetChi2();

	    if(numClustNow<numClustRef && chi2Now>fresult->GetChi2()) continue;
	    //cout<<" chi2Now =  "<<chi2Now<<"\n";   
	      
	    chi2Now/=numClustNow;
	    if(fPtref > 1.0 && chi2Now > 30.) continue; 
	    if((fPtref >= 0.6 && fPtref<=1.0) && chi2Now > 40.) continue;
	    // if((fPtref <= 0.6 && fPtref>0.2)&& chi2Now > 40.) continue;
	    // if(fPtref <= 0.2 && chi2Now > 8.) continue;
	    if((fPtref <= 0.6 && fPtref>0.2)&& chi2Now > 30.) continue;
	    if(fPtref <= 0.2 && chi2Now > 7.) continue;     
	    /////////////////////////////
	} // end if

	Int_t layerInit = (*trackITS).GetLayer();
	Int_t layernew = layerInit - 2;// -1 for new layer, -1 for matrix index
	TList listoftrack;
	Int_t ladp, ladm, detp,detm,ladinters,detinters; 	
	Int_t layerfin=layerInit-1;
	// cout<<"Prima di intersection \n";
	Int_t  outinters=Intersection(*trackITS,layerfin,ladinters,detinters);
	// cout<<" outinters = "<<outinters<<"\n";
	//  cout<<" Layer ladder detector intersection ="
	//      <<layerfin<<" "<<ladinters<<" "<<detinters<<"\n";
	//  cout << " phiinters zinters = "<<(*trackITS)(0) 
	//       << " "<<(*trackITS)(1)<<"\n"; getchar();
	if(outinters==-1) continue;
	Int_t flaghit=0;
	(*trackITS).SetLayer(layerfin);  // oggi
	(*trackITS).Setfcor();           // oggi
	if(outinters==0){
	    TVector toucLad(9), toucDet(9);	 
	    Int_t lycur=layerfin;
	    ladp=ladinters+1;
	    ladm=ladinters-1;
	    if(ladm <= 0) ladm=fNlad[layerfin-1];    
	    if(ladp > fNlad[layerfin-1]) ladp=1;  
	    detp=detinters+1;
	    detm=detinters-1;
	    Int_t idetot=1;
	    /*
	      toucLad(0)=ladinters; toucLad(1)=ladm; toucLad(2)=ladp;
	      toucLad(3)=ladinters; toucLad(4)=ladm; toucLad(5)=ladp;
	      toucLad(6)=ladinters; toucLad(7)=ladm; toucLad(8)=ladp;
	      toucDet(0)=detinters; toucDet(1)=detinters; toucDet(2)=detinters;
	      if(detm > 0 && detp <= fNdet[layerfin-1]) {     
	           idetot=9;
	           toucDet(3)=detm; toucDet(4)=detm; toucDet(5)=detm;	   
	           toucDet(6)=detp; toucDet(7)=detp; toucDet(8)=detp;
	      } // end if
	      if(detm > 0 && detp > fNdet[layerfin-1]) {   
	           idetot=6;
		   toucDet(3)=detm; toucDet(4)=detm; toucDet(5)=detm;
	     } // end if
	     if(detm <= 0 && detp <= fNdet[layerfin-1]) {   
	          idetot=6;
	          toucDet(3)=detp; toucDet(4)=detp; toucDet(5)=detp;
	     } // end if
	    */
 	    Float_t epsphi=5.0, epsz=5.0;                  
	    if(fPtref<0.2) {epsphi=3.; epsz=3.;}     
	    // new definition of idetot e toucLad e toucDet to be
	    // transformed in a method
	    // these values could be modified
	    Float_t pigre=TMath::Pi();
	    Float_t rangephi=5., rangez=5.;
	    if(layerfin==1 || layerfin ==2){
		rangephi=40.*epsphi*TMath::Sqrt(sigmaphil[layerfin-1]+
						(*trackITS).GetSigmaphi());
		rangez = 40.*epsz*TMath::Sqrt(sigmazl[layerfin-1]+
					      (*trackITS).GetSigmaZ());
	    } // end if
	    if(layerfin==3 || layerfin ==4){
		//rangephi=30.*fepsphi*TMath::Sqrt(sigmaphil[layerfin-1]+
		//                                 (*trackITS).GetSigmaphi());
		//rangez = 40.*fepsz*TMath::Sqrt(sigmazl[layerfin-1]+
		//                               (*trackITS).GetSigmaZ());
		rangephi=40.*epsphi*TMath::Sqrt(sigmaphil[layerfin-1]+
						(*trackITS).GetSigmaphi());
		rangez = 50.*epsz*TMath::Sqrt(sigmazl[layerfin-1]+
					      (*trackITS).GetSigmaZ());
	    } // end if
	    if(layerfin==5 || layerfin ==6){
		rangephi=20.*epsphi*TMath::Sqrt(sigmaphil[layerfin-1]+
						(*trackITS).GetSigmaphi());
		rangez =5.*epsz*TMath::Sqrt(sigmazl[layerfin-1]+
					    (*trackITS).GetSigmaZ());
	    } // end if
	    Float_t phinters, zinters;
	    phinters=(*trackITS).Getphi();
	    zinters=(*trackITS).GetZ();
	    Float_t distz = 0.0;
	    Float_t phicm, phicp, distphim, distphip;
	    phicm=phinters;
		 if(phinters>fphimax[layerfin-1][ladm-1]) phicm=phinters-2*pigre;  //corretto il 20-11-2001
		 distphim=TMath::Abs(phicm-fphimax[layerfin-1][ladm-1]);  //corretto il 20-11-2001
	    phicp=phinters;
		 //cout<<" fNlad[layerfin-1] e ladp = "<<fNlad[layerfin-1]<<" "<<ladp<<endl;
		 if(phinters>fphimin[layerfin-1][ladp-1]) phicp=phinters-2.*pigre;   //corretto il 20-11-2001
		 distphip=TMath::Abs(phicp-fphimin[layerfin-1][ladp-1]);      //corretto il 20-11-2001
	    Int_t flagzmin=0;
	    Int_t flagzmax=0;
	    idetot=1;
	    toucLad(0)=ladinters; toucDet(0)=detinters;
	    if(detm>0) distz=TMath::Abs(zinters-fzmax[layerfin-1][detm-1]);
	    if(detm>0 && rangez>=distz){
		flagzmin=1; 
		idetot++; toucLad(idetot-1)=ladinters; toucDet(idetot-1)=detm;
		if(rangephi>=distphim){
		    idetot++; 
		    toucLad(idetot-1)=ladm;
		    toucDet(idetot-1)=detinters;
		    idetot++;
		    toucLad(idetot-1)=ladm;
		    toucDet(idetot-1)=detm;
		} // end if
		if(rangephi>=distphip){
		    idetot++;
		    toucLad(idetot-1)=ladp;
		    toucDet(idetot-1)=detinters;
		    idetot++;
		    toucLad(idetot-1)=ladp;
		    toucDet(idetot-1)=detm;
		} // end if
	    }  //end detm>0....		 
	    if(detp<=fNdet[layerfin-1]) 
		distz=TMath::Abs(zinters-fzmin[layerfin-1][detp-1]);
	    if(detp<=fNdet[layerfin-1] && rangez>=distz){
		flagzmax=1;
		idetot++; toucLad(idetot-1)=ladinters; toucDet(idetot-1)=detp;
		if(rangephi>=distphim){
		    idetot++; toucLad(idetot-1)=ladm; toucDet(idetot-1)=detp;
		    if(flagzmin == 0) {
			idetot++;
			toucLad(idetot-1)=ladm;
			toucDet(idetot-1)=detinters;
		    } // end if
		} // end if
		if(rangephi>=distphip){
		    idetot++;
		    toucLad(idetot-1)=ladp;
		    toucDet(idetot-1)=detp;
		    if(flagzmin == 0) {
			idetot++;
			toucLad(idetot-1)=ladp; 
			toucDet(idetot-1)=detinters;
		    } // end if
		} // end if
	    }  //end detm<fNdet[.......

	    if(flagzmin == 0 && flagzmax==0){
		if(rangephi>=distphim){
		    idetot++; 
		    toucLad(idetot-1)=ladm;
		    toucDet(idetot-1)=detinters;
		} // end if	
		if(rangephi>=distphip){
		    idetot++; 
		    toucLad(idetot-1)=ladp;
		    toucDet(idetot-1)=detinters;
		} // end if	  
	    } // end if
	    ////////////////////////////////////////////////////////////
	    Int_t iriv;
	    for (iriv=0; iriv<idetot; iriv++) {  //for on detectors
		///////////////////////////////////////////////////////
		/*** Rec points sorted by module *****/
		/**************************************/
		Int_t indexmod;       
		indexmod = g1->GetModuleIndex(lycur,(Int_t)toucLad(iriv),
					      (Int_t)toucDet(iriv)); 
		fDetTypeRec->ResetRecPoints();   
		gAlice->TreeR()->GetEvent(indexmod); 
		Int_t npoints=frecPoints->GetEntries();
       
		Int_t indnew;
		for(indnew=0; indnew<npoints; indnew++){
		    if (*(fvettid[indexmod]+indnew)==0)
			recp =(AliITSRecPoint*)frecPoints->UncheckedAt(indnew);
		    else
			continue;
		    TVector cluster(3),vecclust(9);
		    //vecclust(6)=vecclust(7)=vecclust(8)=-1.;
		    Double_t sigma[2];
  // now vecclust is with cylindrical cohordinates
	       vecclust(0)=(Float_t)fRecCylR[indexmod][indnew];     
	       vecclust(1)=(Float_t)fRecCylPhi[indexmod][indnew];
	       vecclust(2)=(Float_t)fRecCylZ[indexmod][indnew];	    			 
		    vecclust(3) = (Float_t)recp->fTracks[0]; 
		    vecclust(4) = (Float_t)indnew;			 
		    vecclust(5) = (Float_t)indexmod;    
		    vecclust(6) = (Float_t)recp->fTracks[0];
		    vecclust(7) = (Float_t)recp->fTracks[1];
		    vecclust(8) = (Float_t)recp->fTracks[2];
		    sigma[0] = (Double_t)  recp->GetSigmaX2();	   
		    sigma[1] = (Double_t) recp->GetSigmaZ2();
			 
			 cluster(0)=fRecCylR[indexmod][indnew];
          cluster(1)=fRecCylPhi[indexmod][indnew];
			 cluster(2)=fRecCylZ[indexmod][indnew];
			 
		    // cout<<" layer = "<<play<<"\n";
		    // cout<<" cluster prima = "<<vecclust(0)<<" "
		    //     <<vecclust(1)<<" "
		    //     <<vecclust(2)<<"\n"; getchar();    

		    Float_t sigmatotphi, sigmatotz;  
		    // Float_t epsphi=5.0, epsz=5.0;                 
		    //if(fPtref<0.2) {epsphi=3.; epsz=3.;}
		    Double_t rTrack=(*trackITS).Getrtrack();
		    Double_t sigmaphi=sigma[0]/(rTrack*rTrack);
		    sigmatotphi=epsphi*TMath::Sqrt(sigmaphi + 
						   (*trackITS).GetSigmaphi());
		    sigmatotz=epsz*TMath::Sqrt(sigma[1] + 
					       (*trackITS).GetSigmaZ());
		    //cout<<"cluster e sigmatotphi e track = "<<cluster(0)
		    //    <<" "<<cluster(1)<<" "<<sigmatotphi<<" "
		    //    <<vecclust(3)<<"\n";
		    //if(vecclust(3)==481) getchar();
		    if(cluster(1)<6. && (*trackITS).Getphi()>6.) 
			cluster(1)=cluster(1)+(2.*TMath::Pi());
		    if(cluster(1)>6. && (*trackITS).Getphi()<6.) 
			cluster(1)=cluster(1)-(2.*TMath::Pi());
		    if(TMath::Abs(cluster(1)-(*trackITS).Getphi())>sigmatotphi)
			continue;
		    // cout<<" supero sigmaphi \n";      
		    AliITSTrackV1 *newTrack = new AliITSTrackV1((*trackITS));
		    //(*newTrack).SetLayer((*trackITS).GetLayer()-1);
		    if (TMath::Abs(rTrack-cluster(0))/rTrack>1e-6) 
			(*newTrack).Correct(Double_t(cluster(0)));	
		    //cout<<" cluster(2) e(*newTrack).GetZ()="<<cluster(2)<<" "
		    //    << (*newTrack).GetZ()<<"\n";
		    if(TMath::Abs(cluster(2)-(*newTrack).GetZ()) > sigmatotz){ 
			delete newTrack;
			continue;
		    } // end if
		    Double_t sigmanew[2];
		    sigmanew[0]= sigmaphi;
		    sigmanew[1]=sigma[1];
		    Double_t m[2];
		    m[0]=cluster(1);
		    m[1]=cluster(2);
		    //  Double_t chi2pred=newTrack->GetPredChi2(m,sigmanew);   
		    // cout<<" chi2pred = "<<chi2pred<<"\n";
		    // if(chi2pred>fChi2max) continue; //aggiunto il 30-7-2001
		    if(iriv == 0) flaghit=1;
		    (*newTrack).AddMS();  // add the multiple scattering 
		                             //matrix to the covariance matrix 
		    (*newTrack).AddEL(1.,0);

		    if(fflagvert){
			KalmanFilterVert(newTrack,cluster,sigmanew);
			//KalmanFilterVert(newTrack,cluster,sigmanew,chi2pred);
		    }else{
			KalmanFilter(newTrack,cluster,sigmanew);
		    } // end if
		    (*newTrack).PutCluster(layernew, vecclust);
		    newTrack->AddClustInTrack();
		    listoftrack.AddLast(newTrack);
		}   // end for indnew
	    }  // end of for on detectors (iriv)
	}//end if(outinters==0)

	if(flaghit==0 || outinters==-2) {
	    AliITSTrackV1 *newTrack = new AliITSTrackV1(*trackITS);
	    (*newTrack).Setfnoclust();        	 
	    //(*newTrack).SetLayer((*trackITS).GetLayer()-1); 
	    (*newTrack).AddMS();  // add the multiple scattering matrix
	                             // to the covariance matrix  
	    (*newTrack).AddEL(1.,0);
	    listoftrack.AddLast(newTrack);
	} // end if

	//gObjectTable->Print();   // stampa memoria
	 
	RecursiveTracking(&listoftrack);          
	listoftrack.Delete();
    } // end of for on tracks (index)

    //gObjectTable->Print();   // stampa memoria
     delete rl;
    
}
//______________________________________________________________________
Int_t AliITSTrackerV1::Intersection(AliITSTrackV1 &track,Int_t layer,
				    Int_t &ladder,Int_t &detector) { 
    // Origin  A. Badala' and G.S. Pappalardo
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
    // Found the intersection and the detector 

    Double_t rk=fAvrad[layer-1];
    if(track.DoNotCross(rk)){ /*cout<< " Do not cross \n";*/ return -1;} 
    track.Propagation(rk);
    Double_t zinters=track.GetZ();
    Double_t phinters=track.Getphi();
    //cout<<"zinters = "<<zinters<<"  phinters = "<<phinters<<"\n";

    TVector det(9);
    TVector listDet(2);
    TVector distZCenter(2);

    Int_t iz=0; 
    Int_t iD;
    for(iD = 1; iD<= fNdet[layer-1]; iD++) {
	if(zinters > fzmin[layer-1][iD-1] && zinters <= fzmax[layer-1][iD-1]) {
	    if(iz>1) {
		cout<< " Errore su iz in Intersection \n";
		getchar();
	    }else {
		listDet(iz)= iD; distZCenter(iz)=TMath::Abs(zinters-det(2));
		iz++;
	    } // end if
	} // end if
    } // end for iD

    if(iz==0) {/* cout<< " No detector along Z \n";*/ return -2;}
    detector=Int_t (listDet(0));
    if(iz>1 && (distZCenter(0)>distZCenter(1)))   detector=Int_t (listDet(1));

    TVector listLad(2);
    TVector distPhiCenter(2);
    Int_t ip=0;
    Double_t pigre=TMath::Pi();
    Int_t iLd;   
    for(iLd = 1; iLd<= fNlad[layer-1]; iLd++) {
	Double_t phimin=fphimin[layer-1][iLd-1];
	Double_t phimax=fphimax[layer-1][iLd-1];
	Double_t phidet=fphidet[layer-1][iLd-1];
	Double_t phiconfr=phinters;
	if(phimin>phimax) {  
	    //if(phimin <5.5) {cout<<" Error in Intersection for phi \n";
	    // getchar();}
	    phimin-=(2.*pigre);
	    if(phinters>(1.5*pigre)) phiconfr=phinters-(2.*pigre); 
	    if(phidet>(1.5*pigre)) phidet-=(2.*pigre);
	} // end if
	if(phiconfr>phimin && phiconfr<= phimax) {
	    if(ip>1) {
		cout<< " Errore su ip in Intersection \n"; getchar();
	    }else  {
		listLad(ip)= iLd;
		distPhiCenter(ip)=TMath::Abs(phiconfr-phidet); ip++;
	    } // end if
	} // end if
    } // end for iLd
    if(ip==0) { cout<< " No detector along phi \n"; getchar();}
    ladder=Int_t (listLad(0));
    if(ip>1 && (distPhiCenter(0)>distPhiCenter(1)))  ladder=Int_t (listLad(1));
    return 0;
}
//______________________________________________________________________
void AliITSTrackerV1::KalmanFilter(AliITSTrackV1 *newTrack,TVector &cluster,
				   Double_t sigma[2]){ 
    //Origin  A. Badala' and G.S. Pappalardo:
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
    // Kalman filter without vertex constraint
    ////// Evaluation of the measurement vector ////////////////////////
    Double_t m[2];
    Double_t rk,phik,zk;
    rk=cluster(0);   phik=cluster(1);  zk=cluster(2);
    m[0]=phik;    m[1]=zk;
    //////////////////////// Evaluation of the error matrix V  /////////
    Double_t v00=sigma[0];
    Double_t v11=sigma[1];
    ////////////////////////////////////////////////////////////////////  
    Double_t cin00,cin10,cin20,cin30,cin40,cin11,cin21,cin31,cin41,cin22,
	     cin32,cin42,cin33,cin43,cin44;

    newTrack->GetCElements(cin00,
			   cin10,cin11,
			   cin20,cin21,cin22,
			   cin30,cin31,cin32,cin33,
			   cin40,cin41,cin42,cin43,cin44); //get C matrix
    Double_t rold00=cin00+v00;
    Double_t rold10=cin10;
    Double_t rold11=cin11+v11;
    ////////////////////// R matrix inversion  /////////////////////////
    Double_t det=rold00*rold11-rold10*rold10;
    Double_t r00=rold11/det;
    Double_t r10=-rold10/det;
    Double_t r11=rold00/det;
    ////////////////////////////////////////////////////////////////////
    Double_t k00=cin00*r00+cin10*r10;
    Double_t k01=cin00*r10+cin10*r11;
    Double_t k10=cin10*r00+cin11*r10;  
    Double_t k11=cin10*r10+cin11*r11;
    Double_t k20=cin20*r00+cin21*r10;  
    Double_t k21=cin20*r10+cin21*r11;  
    Double_t k30=cin30*r00+cin31*r10;  
    Double_t k31=cin30*r10+cin31*r11;  
    Double_t k40=cin40*r00+cin41*r10;
    Double_t k41=cin40*r10+cin41*r11;
    Double_t x0,x1,x2,x3,x4;
    newTrack->GetXElements(x0,x1,x2,x3,x4);     // get the state vector
    Double_t savex0=x0, savex1=x1;
    x0+=k00*(m[0]-savex0)+k01*(m[1]-savex1);
    x1+=k10*(m[0]-savex0)+k11*(m[1]-savex1);
    x2+=k20*(m[0]-savex0)+k21*(m[1]-savex1);
    x3+=k30*(m[0]-savex0)+k31*(m[1]-savex1);
    x4+=k40*(m[0]-savex0)+k41*(m[1]-savex1);
    Double_t c00,c10,c20,c30,c40,c11,c21,c31,c41,c22,c32,c42,c33,c43,c44;
    c00=cin00-k00*cin00-k01*cin10;
    c10=cin10-k00*cin10-k01*cin11;
    c20=cin20-k00*cin20-k01*cin21;
    c30=cin30-k00*cin30-k01*cin31;
    c40=cin40-k00*cin40-k01*cin41;
    c11=cin11-k10*cin10-k11*cin11;
    c21=cin21-k10*cin20-k11*cin21;
    c31=cin31-k10*cin30-k11*cin31;
    c41=cin41-k10*cin40-k11*cin41;
    c22=cin22-k20*cin20-k21*cin21;
    c32=cin32-k20*cin30-k21*cin31;
    c42=cin42-k20*cin40-k21*cin41;
    c33=cin33-k30*cin30-k31*cin31;
    c43=cin43-k30*cin40-k31*cin41;
    c44=cin44-k40*cin40-k41*cin41;
    newTrack->PutXElements(x0,x1,x2,x3,x4);  // put the new state vector
    newTrack->PutCElements(c00,
			   c10,c11,
			   c20,c21,c22,
			   c30,c31,c32,c33,
			   c40,c41,c42,c43,c44); // put in track the
                                                 // new cov matrix 
    Double_t vmcold00=v00-c00;
    Double_t vmcold10=-c10;
    Double_t vmcold11=v11-c11;
    ////////////////////// Matrix vmc inversion  ///////////////////////
    det=vmcold00*vmcold11-vmcold10*vmcold10;
    Double_t vmc00=vmcold11/det;
    Double_t vmc10=-vmcold10/det;
    Double_t vmc11=vmcold00/det;
    ////////////////////////////////////////////////////////////////////
  Double_t chi2=(m[0]-x0)*( vmc00*(m[0]-x0) + 2.*vmc10*(m[1]-x1) ) +
                (m[1]-x1)*vmc11*(m[1]-x1);
  newTrack->SetChi2(newTrack->GetChi2()+chi2);
}
//----------------------------------------------------------------------
//void AliITSTrackerV1::KalmanFilterVert(AliITSTrackV1 *newTrack,
//                                       TVector &cluster,Double_t sigma[2]){
void AliITSTrackerV1::KalmanFilterVert(AliITSTrackV1 *newTrack,
				       TVector &cluster,Double_t sigma[2]
                                       /*, Double_t chi2pred*/){
    //Origin  A. Badala' and G.S. Pappalardo:
    // e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
    // Kalman filter with vertex constraint
    ///////////////////// Evaluation of the measurement vector m ////////
    Double_t m[4];
    Double_t rk,phik,zk;
    rk=cluster(0);   phik=cluster(1);  zk=cluster(2);
    m[0]=phik;    m[1]=zk;
    Double_t cc=(*newTrack).GetC();
    Double_t zv=(*newTrack).GetZv(); 
    Double_t dv=(*newTrack).GetDv();
    Double_t cy=cc/2.;
    Double_t tgl= (zk-zv)*cy/TMath::ASin(cy*rk);
    m[2]=dv;    m[3]=tgl;
    /////////////////////// Evaluation of the error matrix V  //////////
    Int_t layer=newTrack->GetLayer();
    Double_t v00=sigma[0];
    Double_t v11=sigma[1];
    Double_t v31=sigma[1]/rk;
    Double_t sigmaDv=newTrack->GetsigmaDv();
    Double_t v22=sigmaDv*sigmaDv  + newTrack->Getd2(layer-1);
    Double_t v32=newTrack->Getdtgl(layer-1);
    Double_t sigmaZv=newTrack->GetsigmaZv();  
    Double_t v33=(sigma[1]+sigmaZv*sigmaZv)/(rk*rk)+newTrack->Gettgl2(layer-1);
    //////////////////////////////////////////////////////////////////
    Double_t cin00,cin10,cin11,cin20,cin21,cin22,
	     cin30,cin31,cin32,cin33,cin40,cin41,cin42,cin43,cin44;
    newTrack->GetCElements(cin00,
			   cin10,cin11,
			   cin20,cin21,cin22,
			   cin30,cin31,cin32,cin33,
			   cin40,cin41,cin42,cin43,cin44); //get C matrix
    Double_t r[4][4];
    r[0][0]=cin00+v00;
    r[1][0]=cin10;
    r[2][0]=cin20;
    r[3][0]=cin30;
    r[1][1]=cin11+v11;
    r[2][1]=cin21;
    r[3][1]=cin31+sigma[1]/rk;
    r[2][2]=cin22+sigmaDv*sigmaDv+newTrack->Getd2(layer-1);
    r[3][2]=cin32+newTrack->Getdtgl(layer-1);
    r[3][3]=cin33+(sigma[1]+sigmaZv*sigmaZv)/(rk*rk)+
	                            newTrack->Gettgl2(layer-1);
    r[0][1]=r[1][0]; r[0][2]=r[2][0]; r[0][3]=r[3][0]; 
    r[1][2]=r[2][1]; r[1][3]=r[3][1]; r[2][3]=r[3][2];
    /////////////////////  Matrix R inversion //////////////////////////
    const Int_t kn=4;
    Double_t big, hold;
    Double_t d=1.;
    Int_t ll[kn],mm[kn];
    Int_t i,j,k;

    for(k=0; k<kn; k++) {
	ll[k]=k;
	mm[k]=k;
	big=r[k][k];
	for(j=k; j<kn ; j++) {
	    for (i=j; i<kn; i++) {
		if(TMath::Abs(big) < TMath::Abs(r[i][j]) ) { 
		    big=r[i][j]; 
		    ll[k]=i; 
		    mm[k]=j; 
		} // end if
	    } // end for i
	} // end for j
	//
	j= ll[k];
	if(j > k) {
	    for(i=0; i<kn; i++) {
		hold=-r[k][i]; 
		r[k][i]=r[j][i]; 
		r[j][i]=hold;
	    } // end for i
	}// end if
	//
	i=mm[k];
	if(i > k ) { 
	    for(j=0; j<kn; j++) {
		hold=-r[j][k];
		r[j][k]=r[j][i]; 
		r[j][i]=hold; 
	    } // end for j
	} // end if
	//
	if(!big) {
	    d=0.;
	    cout << "Singular matrix\n"; 
	} // end if
	for(i=0; i<kn; i++) {
	    if(i == k) { continue; }    
	    r[i][k]=r[i][k]/(-big);
	} // end for i
	//
	for(i=0; i<kn; i++) {
	    hold=r[i][k];
	    for(j=0; j<kn; j++) {
		if(i == k || j == k) continue;
		r[i][j]=hold*r[k][j]+r[i][j];
	    } // end for j
	} // end of ri
	// 
	for(j=0; j<kn; j++) {
	    if(j == k) continue;
	    r[k][j]=r[k][j]/big;
	} // end for j
	//
	d=d*big;
	//
	r[k][k]=1./big;        
    } // end for k
    //
    for(k=kn-1; k>=0; k--) {
	i=ll[k];
	if(i > k) {
	    for (j=0; j<kn; j++) {
		hold=r[j][k];
		r[j][k]=-r[j][i]; 
		r[j][i]=hold;
	    } // end for j
	} // end if i
	j=mm[k];
	if(j > k) {
	    for (i=0; i<kn; i++) {
		hold=r[k][i]; 
		r[k][i]=-r[j][i]; 
		r[j][i]=hold;
	    } // end for i
	} // end if
    } // end for k
    ////////////////////////////////////////////////////////////////////
    Double_t k00=cin00*r[0][0]+cin10*r[1][0]+cin20*r[2][0]+cin30*r[3][0];
    Double_t k01=cin00*r[1][0]+cin10*r[1][1]+cin20*r[2][1]+cin30*r[3][1];
    Double_t k02=cin00*r[2][0]+cin10*r[2][1]+cin20*r[2][2]+cin30*r[3][2];
    Double_t k03=cin00*r[3][0]+cin10*r[3][1]+cin20*r[3][2]+cin30*r[3][3];
    Double_t k10=cin10*r[0][0]+cin11*r[1][0]+cin21*r[2][0]+cin31*r[3][0];  
    Double_t k11=cin10*r[1][0]+cin11*r[1][1]+cin21*r[2][1]+cin31*r[3][1];
    Double_t k12=cin10*r[2][0]+cin11*r[2][1]+cin21*r[2][2]+cin31*r[3][2];
    Double_t k13=cin10*r[3][0]+cin11*r[3][1]+cin21*r[3][2]+cin31*r[3][3];
    Double_t k20=cin20*r[0][0]+cin21*r[1][0]+cin22*r[2][0]+cin32*r[3][0];  
    Double_t k21=cin20*r[1][0]+cin21*r[1][1]+cin22*r[2][1]+cin32*r[3][1];  
    Double_t k22=cin20*r[2][0]+cin21*r[2][1]+cin22*r[2][2]+cin32*r[3][2];
    Double_t k23=cin20*r[3][0]+cin21*r[3][1]+cin22*r[3][2]+cin32*r[3][3];
    Double_t k30=cin30*r[0][0]+cin31*r[1][0]+cin32*r[2][0]+cin33*r[3][0];  
    Double_t k31=cin30*r[1][0]+cin31*r[1][1]+cin32*r[2][1]+cin33*r[3][1];  
    Double_t k32=cin30*r[2][0]+cin31*r[2][1]+cin32*r[2][2]+cin33*r[3][2];  
    Double_t k33=cin30*r[3][0]+cin31*r[3][1]+cin32*r[3][2]+cin33*r[3][3];
    Double_t k40=cin40*r[0][0]+cin41*r[1][0]+cin42*r[2][0]+cin43*r[3][0];
    Double_t k41=cin40*r[1][0]+cin41*r[1][1]+cin42*r[2][1]+cin43*r[3][1];
    Double_t k42=cin40*r[2][0]+cin41*r[2][1]+cin42*r[2][2]+cin43*r[3][2];  
    Double_t k43=cin40*r[3][0]+cin41*r[3][1]+cin42*r[3][2]+cin43*r[3][3];

    Double_t x0,x1,x2,x3,x4;
    newTrack->GetXElements(x0,x1,x2,x3,x4);     // get the state vector
    Double_t savex0=x0, savex1=x1, savex2=x2, savex3=x3;
    x0+=k00*(m[0]-savex0)+k01*(m[1]-savex1)+k02*(m[2]-savex2)+
	k03*(m[3]-savex3);
    x1+=k10*(m[0]-savex0)+k11*(m[1]-savex1)+k12*(m[2]-savex2)+
	k13*(m[3]-savex3);
    x2+=k20*(m[0]-savex0)+k21*(m[1]-savex1)+k22*(m[2]-savex2)+
	k23*(m[3]-savex3);
    x3+=k30*(m[0]-savex0)+k31*(m[1]-savex1)+k32*(m[2]-savex2)+
	k33*(m[3]-savex3);
    x4+=k40*(m[0]-savex0)+k41*(m[1]-savex1)+k42*(m[2]-savex2)+
	k43*(m[3]-savex3);
    Double_t c00,c10,c20,c30,c40,c11,c21,c31,c41,c22,c32,c42,c33,c43,c44;
    c00=cin00-k00*cin00-k01*cin10-k02*cin20-k03*cin30;
    c10=cin10-k00*cin10-k01*cin11-k02*cin21-k03*cin31;
    c20=cin20-k00*cin20-k01*cin21-k02*cin22-k03*cin32;
    c30=cin30-k00*cin30-k01*cin31-k02*cin32-k03*cin33;
    c40=cin40-k00*cin40-k01*cin41-k02*cin42-k03*cin43;
    c11=cin11-k10*cin10-k11*cin11-k12*cin21-k13*cin31;
    c21=cin21-k10*cin20-k11*cin21-k12*cin22-k13*cin32;
    c31=cin31-k10*cin30-k11*cin31-k12*cin32-k13*cin33;
    c41=cin41-k10*cin40-k11*cin41-k12*cin42-k13*cin43;
    c22=cin22-k20*cin20-k21*cin21-k22*cin22-k23*cin32;
    c32=cin32-k20*cin30-k21*cin31-k22*cin32-k23*cin33;
    c42=cin42-k20*cin40-k21*cin41-k22*cin42-k23*cin43;
    c33=cin33-k30*cin30-k31*cin31-k32*cin32-k33*cin33;
    c43=cin43-k30*cin40-k31*cin41-k32*cin42-k33*cin43;
    c44=cin44-k40*cin40-k41*cin41-k42*cin42-k43*cin43;

    newTrack->PutXElements(x0,x1,x2,x3,x4); // put the new state vector
    newTrack->PutCElements(c00,
			   c10,c11,
			   c20,c21,c22,
			   c30,c31,c32,c33,
			   c40,c41,c42,c43,c44); // put in track the
                                                 // new cov matrix
    Double_t vmc[4][4];
    vmc[0][0]=v00-c00; vmc[1][0]=-c10; vmc[2][0]=-c20; vmc[3][0]=-c30;
    vmc[1][1]=v11-c11; vmc[2][1]=-c21; vmc[3][1]=v31-c31;
    vmc[2][2]=v22-c22; vmc[3][2]=v32-c32;
    vmc[3][3]=v33-c33;
    vmc[0][1]=vmc[1][0]; vmc[0][2]=vmc[2][0]; vmc[0][3]=vmc[3][0];
    vmc[1][2]=vmc[2][1]; vmc[1][3]=vmc[3][1];
    vmc[2][3]=vmc[3][2];
    /////////////////////// vmc matrix inversion /////////////////////// 
    d=1.;
    for(k=0; k<kn; k++) {
	ll[k]=k;
	mm[k]=k;
	big=vmc[k][k];
	for(j=k; j<kn ; j++) {
	    for (i=j; i<kn; i++) {
		if(TMath::Abs(big) < TMath::Abs(vmc[i][j]) ) {
		    big=vmc[i][j]; 
		    ll[k]=i;
		    mm[k]=j;
		} // end if
	    } // end for i
	} // end for j
	//
	j= ll[k];
	if(j > k) {
	    for(i=0; i<kn; i++) {
		hold=-vmc[k][i]; 
		vmc[k][i]=vmc[j][i];
		vmc[j][i]=hold;
	    } // end for i
	} // end if
	//
	i=mm[k];
	if(i > k ) { 
	    for(j=0; j<kn; j++) { 
		hold=-vmc[j][k]; 
		vmc[j][k]=vmc[j][i]; 
		vmc[j][i]=hold; 
	    } // end for j
	} // end if
	//
	if(!big) {
	    d=0.;
	    cout << "Singular matrix\n"; 
	} // end if
	for(i=0; i<kn; i++) {
	    if(i == k) continue;    
	    vmc[i][k]=vmc[i][k]/(-big);
	}    // end for i
	//
	for(i=0; i<kn; i++) {
	    hold=vmc[i][k];
	    for(j=0; j<kn; j++) {
		if(i == k || j == k) continue;
		vmc[i][j]=hold*vmc[k][j]+vmc[i][j];
	    } // end for j
	} // end for i
	//  
	for(j=0; j<kn; j++) {
	    if(j == k) continue;
	    vmc[k][j]=vmc[k][j]/big;
	} // end for j
	//
	d=d*big;
	//
	vmc[k][k]=1./big;        
    }  // end for k
    //  
    for(k=kn-1; k>=0; k--) {
	i=ll[k];
	if(i > k) {
	    for (j=0; j<kn; j++) {
		hold=vmc[j][k]; 
		vmc[j][k]=-vmc[j][i];
		vmc[j][i]=hold;
	    } // end for j
	} // end if i>k
	j=mm[k];
	if(j > k) {
	    for (i=0; i<kn; i++) {
		hold=vmc[k][i]; 
		vmc[k][i]=-vmc[j][i]; 
		vmc[j][i]=hold;
	    } // end for i
	} // end if j>k
    } // end for k
    ////////////////////////////////////////////////////////////////////
    Double_t chi2=(m[0]-x0)*( vmc[0][0]*(m[0]-x0) + 2.*vmc[1][0]*(m[1]-x1) + 
			      2.*vmc[2][0]*(m[2]-x2)+ 2.*vmc[3][0]*(m[3]-x3) )+
	(m[1]-x1)* ( vmc[1][1]*(m[1]-x1) + 2.*vmc[2][1]*(m[2]-x2)+ 
		     2.*vmc[3][1]*(m[3]-x3) ) +
	(m[2]-x2)* ( vmc[2][2]*(m[2]-x2)+ 2.*vmc[3][2]*(m[3]-x3) ) +
	(m[3]-x3)*vmc[3][3]*(m[3]-x3);
    //cout<<" chi2 kalman = "<<chi2<<"\n";  getchar(); 
    newTrack->SetChi2(newTrack->GetChi2()+chi2);   
    //   newTrack->SetChi2(newTrack->GetChi2()+chi2pred);
}

