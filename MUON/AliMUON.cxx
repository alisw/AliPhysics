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

/*
$Log$
Revision 1.16  2000/04/26 10:17:31  fca
Changes in Lego for G4 compatibility

Revision 1.15  2000/01/19 17:16:56  fca
Introducing a list of lists of hits -- more hits allowed for detector now

Revision 1.14  1999/11/03 13:17:07  fca
Have ProdProcess return const char*

Revision 1.13  1999/10/26 06:04:48  fca
Introduce TLorentzVector in AliMC::GetSecondary. Thanks to I.Hrivnacova

Revision 1.12  1999/10/07 21:08:10  fca
Corrections by G.Chabratova

Revision 1.11  1999/10/05 17:15:45  fca
Minor syntax for the Alpha OSF

Revision 1.10  1999/10/01 09:24:40  fca
Protect against no current file in FinishEvent

Revision 1.9  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

////////////////////////////////////////////////
//  Manager and hits classes for set:MUON     //
////////////////////////////////////////////////

#include <TH1.h>
#include <TH2.h>
#include <TTUBE.h>
#include <TBRIK.h>
#include <TRotMatrix.h>
#include <TNode.h> 
#include <TTree.h> 
#include <TRandom.h> 
#include <TObject.h>
#include <TVector.h>
#include <TObjArray.h>
#include <TMinuit.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TDirectory.h>
#include <TObjectTable.h>
#include <AliPDG.h>

#include "AliMUON.h"
#include "TTUBE.h"
#include "AliMUONClusterFinder.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h" 

#ifndef WIN32 
# define reco_init       reco_init_
# define cutpxz          cutpxz_
# define sigmacut        sigmacut_
# define xpreci          xpreci_
# define ypreci          ypreci_
# define reconstmuon     reconstmuon_
# define trackf_read_geant     trackf_read_geant_
# define trackf_read_spoint     trackf_read_spoint_
# define chfill          chfill_
# define chfill2         chfill2_
# define chf1            chf1_
# define chfnt           chfnt_
# define hist_create     hist_create_
# define hist_closed     hist_closed_
# define rndm            rndm_
# define fcn             fcn_
# define trackf_fit      trackf_fit_
# define prec_fit        prec_fit_
# define fcnfit          fcnfit_
# define reco_term       reco_term_
#else 
# define reco_init       RECO_INIT
# define cutpxz          CUTPXZ
# define sigmacut        SIGMACUT
# define xpreci          XPRECI
# define ypreci          YPRECI
# define reconstmuon     RECONSTMUON
# define trackf_read_geant     TRACKF_READ_GEANT
# define trackf_read_spoint     TRACKF_READ_SPOINT
# define chfill          CHFILL
# define chfill2         CHFILL2
# define chf1            CHF1
# define chfnt           CHFNT
# define hist_create     HIST_CREATE
# define hist_closed     HIST_CLOSED
# define rndm            RNDM
# define fcn             FCN
# define trackf_fit      TRACKF_FIT
# define prec_fit        PREC_FIT
# define fcnfit          FCNFIT
# define reco_term       RECO_TERM
#endif 

extern "C"
{
void type_of_call reco_init(Double_t &, Double_t &, Double_t &);
void type_of_call reco_term();
void type_of_call cutpxz(Double_t &);
void type_of_call sigmacut(Double_t &);
void type_of_call xpreci(Double_t &);
void type_of_call ypreci(Double_t &);
void type_of_call reconstmuon(Int_t &, Int_t &, Int_t &, Int_t &, Int_t &);
void type_of_call trackf_read_geant(Int_t *, Double_t *, Double_t *, Double_t *, Int_t *, Int_t *, Double_t *, Double_t *, Double_t *, Double_t *,Int_t &, Double_t *, Double_t *, Double_t *, Int_t &, Int_t &, Double_t *, Double_t *, Double_t *, Double_t *); 
void type_of_call trackf_read_spoint(Int_t *, Double_t *, Double_t *, Double_t *, Int_t *, Int_t *, Double_t *, Double_t *, Double_t *, Double_t *,Int_t &, Double_t *, Double_t *, Double_t *, Int_t &, Int_t &, Double_t *, Double_t *, Double_t *, Double_t *); 
void type_of_call chfill(Int_t &, Float_t &, Float_t &, Float_t &);
void type_of_call chfill2(Int_t &, Float_t &, Float_t &, Float_t &);
void type_of_call chf1(Int_t &, Float_t &, Float_t &);
void type_of_call chfnt(Int_t &, Int_t &, Int_t *, Int_t *, Float_t *, Float_t *, Float_t *, Float_t *, Float_t *, Float_t *, Float_t *, Float_t *);
void type_of_call hist_create();
void type_of_call hist_closed();
void type_of_call fcnf(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);
void type_of_call fcn(Int_t &, Double_t *, Double_t &, Double_t *, Int_t &, Int_t &);
void type_of_call trackf_fit(Int_t &, Double_t *, Double_t *, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &);
void type_of_call prec_fit(Double_t &, Double_t &, Double_t &, Double_t &, Double_t&, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &);
void type_of_call fcnfitf(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);
void type_of_call fcnfit(Int_t &, Double_t *, Double_t &, Double_t *, Int_t &, Int_t &);
Float_t type_of_call rndm() {return gRandom->Rndm();}
}

void fcnfwrap(Int_t &i1, Double_t *d1, Double_t &d2,
	Double_t *d3, Int_t i2)
{
   fcnf(i1,d1,d2,d3,i2);
}

void fcnfitfwrap(Int_t &i1, Double_t *d1, Double_t &d2,
	Double_t *d3, Int_t i2)
{
   fcnfitf(i1,d1,d2,d3,i2);
}


// Static variables for the pad-hit iterator routines
static Int_t sMaxIterPad=0;
static Int_t sCurIterPad=0;
static TTree *TrH1;
static TTree *TK1;
static TClonesArray *fHits2;        //Listof hits for one track only
static TClonesArray *fClusters2;    //List of clusters for one track only
static TClonesArray *fParticles2;   //List of particles in the Kine tree
ClassImp(AliMUON)
//___________________________________________
AliMUON::AliMUON()
{
   fIshunt     = 0;
   fHits       = 0;
   fClusters   = 0;
   fNclusters  = 0;
   fDchambers  = 0;
   fNdch       = 0;
   fRawClusters= 0;
   fNrawch     = 0;
   fCathCorrel= 0;
   fNcorch     = 0;
   fTreeC = 0;

   // modifs perso
   fSPxzCut    = 0;
   fSSigmaCut  = 0;
   fSXPrec     = 0; 
   fSYPrec     = 0;
}
 
//___________________________________________
AliMUON::AliMUON(const char *name, const char *title)
       : AliDetector(name,title)
{
//Begin_Html
/*
<img src="gif/alimuon.gif">
*/
//End_Html
 
   fHits     = new TClonesArray("AliMUONhit",1000);
   gAlice->AddHitList(fHits);
   fClusters = new TClonesArray("AliMUONcluster",10000);
   fNclusters  =  0;
   fIshunt     =  0;

   fNdch      = new Int_t[10];

   fDchambers = new TObjArray(10);

   Int_t i;
   
   for (i=0; i<10 ;i++) {
       (*fDchambers)[i] = new TClonesArray("AliMUONdigit",10000); 
       fNdch[i]=0;
   }

   fNrawch      = new Int_t[10];

   fRawClusters = new TObjArray(10);

   for (i=0; i<10 ;i++) {
       (*fRawClusters)[i] = new TClonesArray("AliMUONRawCluster",10000); 
       fNrawch[i]=0;
   }

   fNcorch      = new Int_t[10];
   fCathCorrel = new TObjArray(10);
   for (i=0; i<10 ;i++) {
       (*fCathCorrel)[i] = new TClonesArray("AliMUONcorrelation",1000); 
       fNcorch[i]=0;
   }

   fTreeC = 0;

//   
// Transport angular cut
   fAccCut=0;
   fAccMin=2;
   fAccMax=9;

   // modifs perso
   fSPxzCut   = 3.0;
   fSSigmaCut = 1.0;
   fSXPrec    = 0.01; 
   fSYPrec    = 0.144;

   SetMarkerColor(kRed);
}
 
//___________________________________________
AliMUON::~AliMUON()
{

    printf("Calling AliMUON destructor !!!\n");
    
  Int_t i;
  fIshunt  = 0;
  delete fHits;
  delete fClusters;
  delete fTreeC;

  for (i=0;i<10;i++) {
      delete (*fDchambers)[i];
      fNdch[i]=0;
  }
  delete fDchambers;

  for (i=0;i<10;i++) {
      delete (*fRawClusters)[i];
      fNrawch[i]=0;
  }
  delete fRawClusters;

  for (i=0;i<10;i++) {
      delete (*fCathCorrel)[i];
      fNcorch[i]=0;
  }
  delete fCathCorrel;
}
 
//___________________________________________
void AliMUON::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONhit(fIshunt,track,vol,hits);
}
//___________________________________________
void AliMUON::AddCluster(Int_t *clhits)
{
   TClonesArray &lclusters = *fClusters;
   new(lclusters[fNclusters++]) AliMUONcluster(clhits);
}
//_____________________________________________________________________________
void AliMUON::AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Add a MUON digit to the list
    //

    TClonesArray &ldigits = *((TClonesArray*)(*fDchambers)[id]);
    new(ldigits[fNdch[id]++]) AliMUONdigit(tracks,charges,digits);
}

//_____________________________________________________________________________
void AliMUON::AddRawCluster(Int_t id, const AliMUONRawCluster& c)
{
    //
    // Add a MUON digit to the list
    //

    TClonesArray &lrawcl = *((TClonesArray*)(*fRawClusters)[id]);
    new(lrawcl[fNrawch[id]++]) AliMUONRawCluster(c);
}
//_____________________________________________________________________________
void AliMUON::AddCathCorrel(Int_t id, Int_t *idx, Float_t *x, Float_t *y)
{
    //
    // Add a MUON digit to the list
    //

    TClonesArray &lcorrel = *((TClonesArray*)(*fCathCorrel)[id]);
    new(lcorrel[fNcorch[id]++]) AliMUONcorrelation(idx,x,y);
}

//___________________________________________
void AliMUON::BuildGeometry()
{
    TNode *Node, *NodeF, *Top;
    const int kColorMUON = kBlue;
    //
    Top=gAlice->GetGeometry()->GetNode("alice");
// MUON
//
//  z-Positions of Chambers
    const Float_t cz[5]={511., 686., 971., 1245., 1445.};
//
//  inner diameter
    const Float_t dmi[5]={ 35.,  47.,  67.,   86.,  100.};
//
//  outer diameter
    const Float_t dma[5]={183., 245., 346.,  520.,  520.};

    TRotMatrix* rot000 = new TRotMatrix("Rot000"," ", 90,  0, 90, 90, 0, 0);
    TRotMatrix* rot090 = new TRotMatrix("Rot090"," ", 90, 90, 90,180, 0, 0);
    TRotMatrix* rot180 = new TRotMatrix("Rot180"," ", 90,180, 90,270, 0, 0);
    TRotMatrix* rot270 = new TRotMatrix("Rot270"," ", 90,270, 90,  0, 0, 0);
    

    float rmin, rmax, dx, dy, dz, dr, zpos;
    float dzc=4.;
    char NameChamber[9], NameSense[9], NameFrame[9], NameNode[8];
    for (Int_t i=0; i<5; i++) {
	for (Int_t j=0; j<2; j++) {
	    Int_t id=2*i+j+1;
	    if (j==0) {
		zpos=cz[i]-dzc;
	    } else {
		zpos=cz[i]+dzc;
	    }
	    
	    
	    sprintf(NameChamber,"C_MUON%d",id);
	    sprintf(NameSense,"S_MUON%d",id);
	    sprintf(NameFrame,"F_MUON%d",id);	
	    rmin = dmi[i]/2.-3;
	    rmax = dma[i]/2.+3;
	    new TTUBE(NameChamber,"Mother","void",rmin,rmax,0.25,1.);
	    rmin = dmi[i]/2.;
	    rmax = dma[i]/2.;
	    new TTUBE(NameSense,"Sens. region","void",rmin,rmax,0.25, 1.);
	    dx=(rmax-rmin)/2;
	    dy=3.;
	    dz=0.25;
	    TBRIK* FMUON = new TBRIK(NameFrame,"Frame","void",dx,dy,dz);
	    Top->cd();
	    sprintf(NameNode,"MUON%d",100+id);
	    Node = new TNode(NameNode,"ChamberNode",NameChamber,0,0,zpos,"");
	    Node->SetLineColor(kColorMUON);
	    fNodes->Add(Node);
	    Node->cd();
	    sprintf(NameNode,"MUON%d",200+id);
	    Node = new TNode(NameNode,"Sens. Region Node",NameSense,0,0,0,"");
	    Node->SetLineColor(kColorMUON);
	    fNodes->Add(Node);
	    Node->cd();
	    dr=dx+rmin;
	    sprintf(NameNode,"MUON%d",300+id);
	    NodeF = new TNode(NameNode,"Frame0",FMUON,dr, 0, 0,rot000,"");
	    NodeF->SetLineColor(kColorMUON);
	    fNodes->Add(NodeF);
	    Node->cd();
	    sprintf(NameNode,"MUON%d",400+id);
	    NodeF = new TNode(NameNode,"Frame1",FMUON,0 ,dr,0,rot090,"");
	    NodeF->SetLineColor(kColorMUON);
	    fNodes->Add(NodeF);
	    Node->cd();
	    sprintf(NameNode,"MUON%d",500+id);
	    NodeF = new TNode(NameNode,"Frame2",FMUON,-dr,0,0,rot180,"");
	    NodeF->SetLineColor(kColorMUON);
	    fNodes->Add(NodeF);
	    Node  ->cd();
	    sprintf(NameNode,"MUON%d",600+id);   
	    NodeF = new TNode(NameNode,"Frame3",FMUON,0,-dr,0,rot270,"");
	    NodeF->SetLineColor(kColorMUON);
	    fNodes->Add(NodeF);
	}
    }
}


//___________________________________________
Int_t AliMUON::DistancetoPrimitive(Int_t , Int_t )
{
   return 9999;
}

//___________________________________________
void AliMUON::MakeBranch(Option_t* option)
{
  // Create Tree branches for the MUON.
  
  const Int_t buffersize = 4000;
  char branchname[30];
  sprintf(branchname,"%sCluster",GetName());

  AliDetector::MakeBranch(option);

  if (fClusters   && gAlice->TreeH()) {
    gAlice->TreeH()->Branch(branchname,&fClusters, buffersize);
    printf("Making Branch %s for clusters\n",branchname);
  }

// one branch for digits per chamber
  Int_t i;
  
  for (i=0; i<10 ;i++) {
      sprintf(branchname,"%sDigits%d",GetName(),i+1);
      
      if (fDchambers   && gAlice->TreeD()) {
	  gAlice->TreeD()->Branch(branchname,&((*fDchambers)[i]), buffersize);
	  printf("Making Branch %s for digits in chamber %d\n",branchname,i+1);
      }	
  }

  //printf("Make Branch - TreeR address %p\n",gAlice->TreeR());

// one branch for raw clusters per chamber
  for (i=0; i<10 ;i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
      
      if (fRawClusters   && gAlice->TreeR()) {
	 gAlice->TreeR()->Branch(branchname,&((*fRawClusters)[i]), buffersize);
	 printf("Making Branch %s for raw clusters in chamber %d\n",branchname,i+1);
      }	
  }

}

//___________________________________________
void AliMUON::SetTreeAddress()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[30];
  AliDetector::SetTreeAddress();

  TBranch *branch;
  TTree *treeH = gAlice->TreeH();
  TTree *treeD = gAlice->TreeD();
  TTree *treeR = gAlice->TreeR();

  if (treeH) {
    if (fClusters) {
      branch = treeH->GetBranch("MUONCluster");
      if (branch) branch->SetAddress(&fClusters);
    }
  }

  if (treeD) {
      for (int i=0; i<10; i++) {
	  sprintf(branchname,"%sDigits%d",GetName(),i+1);
	  if (fDchambers) {
	      branch = treeD->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fDchambers)[i]));
	  }
      }
  }

  // printf("SetTreeAddress --- treeR address  %p \n",treeR);

  if (treeR) {
      for (int i=0; i<10; i++) {
	  sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	  if (fRawClusters) {
	      branch = treeR->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	  }
      }
  }

}
//___________________________________________
void AliMUON::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
  AliDetector::ResetHits();
  fNclusters = 0;
  if (fClusters) fClusters->Clear();
}

//____________________________________________
void AliMUON::ResetDigits()
{
    //
    // Reset number of digits and the digits array for this detector
    //
    for ( int i=0;i<10;i++ ) {
	if ((*fDchambers)[i])    ((TClonesArray*)(*fDchambers)[i])->Clear();
	if (fNdch)  fNdch[i]=0;
    }
}
//____________________________________________
void AliMUON::ResetRawClusters()
{
    //
    // Reset number of raw clusters and the raw clust array for this detector
    //
    for ( int i=0;i<10;i++ ) {
	if ((*fRawClusters)[i])    ((TClonesArray*)(*fRawClusters)[i])->Clear();
	if (fNrawch)  fNrawch[i]=0;
    }
}
//____________________________________________
void AliMUON::ResetCorrelation()
{
    //
    // Reset number of correl clusters and the correl clust array for 
    // this detector
    //
    for ( int i=0;i<10;i++ ) {
	if ((*fCathCorrel)[i])   ((TClonesArray*)(*fCathCorrel)[i])->Clear();
	if (fNcorch)  fNcorch[i]=0;
    }
}

//___________________________________________

void AliMUON::SetPADSIZ(Int_t id, Int_t isec, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])  ->SetPADSIZ(isec,p1,p2);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetPADSIZ(isec,p1,p2);
}

//___________________________________________
void AliMUON::SetChargeSlope(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])->SetChargeSlope(p1);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetChargeSlope(p1);
}

//___________________________________________
void AliMUON::SetChargeSpread(Int_t id, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])->SetChargeSpread(p1,p2);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetChargeSpread(p1,p2);
}

//___________________________________________
void AliMUON::SetSigmaIntegration(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])->SetSigmaIntegration(p1);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetSigmaIntegration(p1);
}

//___________________________________________
void AliMUON::SetMaxAdc(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])->SetMaxAdc(p1);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetMaxAdc(p1);
}

//___________________________________________
void AliMUON::SetMaxStepGas(Float_t p1)
{
     fMaxStepGas=p1;
}

//___________________________________________
void AliMUON::SetMaxStepAlu(Float_t p1)
{
    fMaxStepAlu=p1;
}

//___________________________________________
void AliMUON::SetMaxDestepGas(Float_t p1)
{
    fMaxDestepGas=p1;
}

//___________________________________________
void AliMUON::SetMaxDestepAlu(Float_t p1)
{
    fMaxDestepAlu=p1;
}
//___________________________________________
void AliMUON::SetMuonAcc(Bool_t acc, Float_t angmin, Float_t angmax)
{
   fAccCut=acc;
   fAccMin=angmin;
   fAccMax=angmax;
}
//___________________________________________
void   AliMUON::SetSegmentationModel(Int_t id, Int_t isec, AliMUONsegmentation *segmentation)
{
    ((AliMUONchamber*) (*fChambers)[id])->SegmentationModel(isec, segmentation);

}
//___________________________________________
void   AliMUON::SetResponseModel(Int_t id, AliMUONresponse *response)
{
    ((AliMUONchamber*) (*fChambers)[id])->ResponseModel(response);
}

void   AliMUON::SetReconstructionModel(Int_t id, AliMUONClusterFinder *reconst)
{
    ((AliMUONchamber*) (*fChambers)[id])->ReconstructionModel(reconst);
}

void   AliMUON::SetNsec(Int_t id, Int_t nsec)
{
    ((AliMUONchamber*) (*fChambers)[id])->SetNsec(nsec);
}


//___________________________________________

void AliMUON::StepManager()
{
    printf("Dummy version of muon step -- it should never happen!!\n");
    /*
    const Float_t kRaddeg = 180/TMath::Pi();
    Int_t nsec, ipart;
    TLorentzVector x, p;
    Float_t pt, th0, th2;
    char *proc;
    if(fAccCut) {
	if((nsec=gMC->NSecondaries())>0) {
	    proc=gMC->ProdProcess();
	    if((gMC->TrackPid()==443 || gMC->TrackPid()==553) && !strcmp(proc,"DCAY")) {
		//
		// Check angular acceptance
		// --- and have muons from resonance decays in the wanted window --- 
		if(nsec != 2) {
		    printf(" AliMUON::StepManager: Strange resonance Decay into %d particles\n",nsec);
		    gMC->StopEvent();
		} else {
		    gMC->GetSecondary(0,ipart,x,p);
		    pt  = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
		    th0 = TMath::ATan2(pt,p[2])*kRaddeg;
	 	    gMC->GetSecondary(1,ipart,x,p);
		    pt  = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
		    th2 = TMath::ATan2(pt,p[2])*kRaddeg;
		    if(!(fAccMin < th0 && th0 < fAccMax) ||
		       !(fAccMin < th2 && th2 < fAccMax)) 
			gMC->StopEvent();
		}
	    }
	}
    }
    */
}

void AliMUON::MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol)
{
//
//  Calls the charge disintegration method of the current chamber and adds
//  the simulated cluster to the root treee 
//
    Int_t clhits[7];
    Float_t newclust[6][500];
    Int_t nnew;
    
    
//
//  Integrated pulse height on chamber

    
    clhits[0]=fNhits+1;
//
//
    ((AliMUONchamber*) (*fChambers)[idvol])->DisIntegration(eloss, xhit, yhit, nnew, newclust);
//    printf("\n Add new clusters %d %f \n", nnew, eloss*1.e9);
    Int_t ic=0;
    
//
//  Add new clusters
    for (Int_t i=0; i<nnew; i++) {
	if (Int_t(newclust[3][i]) > 0) {
	    ic++;
// Cathode plane
	    clhits[1] = Int_t(newclust[5][i]);
//  Cluster Charge
	    clhits[2] = Int_t(newclust[0][i]);
//  Pad: ix
	    clhits[3] = Int_t(newclust[1][i]);
//  Pad: iy 
	    clhits[4] = Int_t(newclust[2][i]);
//  Pad: charge
	    clhits[5] = Int_t(newclust[3][i]);
//  Pad: chamber sector
	    clhits[6] = Int_t(newclust[4][i]);
	    
	    AddCluster(clhits);
	}
    }
//    printf("\n %d new clusters added", ic);
}

void AliMUON::Digitise(Int_t nev,Int_t bgr_ev,Option_t *option, Option_t *,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
  
    static Bool_t first=kTRUE;
//    static TTree *TrH1;
    static TFile *File;
    char *Add = strstr(option,"Add");
    //char *listoftracks = strstr(opt,"listoftracks");

    AliMUONchamber*  iChamber;
    AliMUONsegmentation*  segmentation;

    
    Int_t trk[50];
    Int_t chtrk[50];  
    TObjArray *list=new TObjArray;
    static TClonesArray *p_adr=0;
    if(!p_adr) p_adr=new TClonesArray("TVector",1000);
    Int_t digits[5]; 

    AliMUON *MUON  = (AliMUON *) gAlice->GetModule("MUON");
    AliMUONHitMap * HitMap[10];
    for (Int_t i=0; i<10; i++) {HitMap[i]=0;}
    if (Add ) {
	if(first) {
	    fFileName=filename;
	    cout<<"filename"<<fFileName<<endl;
	    File=new TFile(fFileName);
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
	    fHits2     = new TClonesArray("AliMUONhit",1000  );
	    fClusters2 = new TClonesArray("AliMUONcluster",10000);
	}	    
	first=kFALSE;
	File->cd();
	//File->ls();
	// Get Hits Tree header from file
	if(fHits2) fHits2->Clear();
	if(fClusters2) fClusters2->Clear();
	if(TrH1) delete TrH1;
	TrH1=0;
	
	char treeName[20];
	sprintf(treeName,"TreeH%d",bgr_ev);
	TrH1 = (TTree*)gDirectory->Get(treeName);
        //printf("TrH1 %p of treename %s for event %d \n",TrH1,treeName,bgr_ev);
	
	if (!TrH1) {
	    printf("ERROR: cannot find Hits Tree for event:%d\n",bgr_ev);
	}
	// Set branch addresses
	TBranch *branch;
	char branchname[20];
	sprintf(branchname,"%s",GetName());
	if (TrH1 && fHits2) {
	    branch = TrH1->GetBranch(branchname);
	    if (branch) branch->SetAddress(&fHits2);
	}
	if (TrH1 && fClusters2) {
	    branch = TrH1->GetBranch("MUONCluster");
	    if (branch) branch->SetAddress(&fClusters2);
	}
// test
	//Int_t ntracks1 =(Int_t)TrH1->GetEntries();
	//printf("background - ntracks1 - %d\n",ntracks1);
    }
    //
    // loop over cathodes
    //
    AliMUONHitMap* hm;
    Int_t countadr=0;
    for (int icat=0; icat<2; icat++) { 
	Int_t counter=0;
	for (Int_t i =0; i<10; i++) {
	    iChamber=(AliMUONchamber*) (*fChambers)[i];
	    if (iChamber->Nsec()==1 && icat==1) {
		continue;
	    } else {
		segmentation=iChamber->GetSegmentationModel(icat+1);
	    }
	    HitMap[i] = new AliMUONHitMapA1(segmentation, list);
	}
	//printf("Start loop over tracks \n");     
//
//   Loop over tracks
//

	TTree *TH = gAlice->TreeH();
	Int_t ntracks =(Int_t) TH->GetEntries();
        //printf("signal - ntracks %d\n",ntracks);
	Int_t nmuon[10]={0,0,0,0,0,0,0,0,0,0};
	Float_t xhit[10][2];
	Float_t yhit[10][2];
	
	for (Int_t track=0; track<ntracks; track++) {
	    gAlice->ResetHits();
	    TH->GetEvent(track);
	    
//
//   Loop over hits
	    for(AliMUONhit* mHit=(AliMUONhit*)MUON->FirstHit(-1); 
		mHit;
		mHit=(AliMUONhit*)MUON->NextHit()) 
	    {
		Int_t   nch   = mHit->fChamber-1;  // chamber number
		if (nch >9) continue;
		iChamber = &(MUON->Chamber(nch));
		Int_t rmin = (Int_t)iChamber->RInner();
		Int_t rmax = (Int_t)iChamber->ROuter();
                // new 17.07.99
		if (Add) {

		  if (mHit->fParticle == kMuonPlus || mHit->fParticle == kMuonMinus) {
		    xhit[nch][nmuon[nch]]=mHit->fX;
		    yhit[nch][nmuon[nch]]=mHit->fY;
		    nmuon[nch]++;
                    if (nmuon[nch] >2) printf("nmuon %d\n",nmuon[nch]);
		    
		  }
		}



		
//
// Loop over pad hits
		for (AliMUONcluster* mPad=
			 (AliMUONcluster*)MUON->FirstPad(mHit,fClusters);
		     mPad;
		     mPad=(AliMUONcluster*)MUON->NextPad(fClusters))
		{
		    Int_t cathode  = mPad->fCathode;    // cathode number
		    Int_t ipx      = mPad->fPadX;       // pad number on X
		    Int_t ipy      = mPad->fPadY;       // pad number on Y
		    Int_t iqpad    = Int_t(mPad->fQpad*kScale);// charge per pad
//		    Int_t iqpad    = mPad->fQpad;       // charge per pad
//
//
		    
		    if (cathode != (icat+1)) continue;
		    // fill the info array
		    Float_t thex, they;
		    segmentation=iChamber->GetSegmentationModel(cathode);
		    segmentation->GetPadCxy(ipx,ipy,thex,they);
		    Float_t rpad=TMath::Sqrt(thex*thex+they*they);
		    if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;

		    new((*p_adr)[countadr++]) TVector(2);
		    TVector &trinfo=*((TVector*) (*p_adr)[countadr-1]);
		    trinfo(0)=(Float_t)track;
		    trinfo(1)=(Float_t)iqpad;

		    digits[0]=ipx;
		    digits[1]=ipy;
		    digits[2]=iqpad;
		    digits[3]=iqpad;
		    if (mHit->fParticle == kMuonPlus || mHit->fParticle == kMuonMinus) {
		    digits[4]=mPad->fHitNumber;
		    } else digits[4]=-1;

		    AliMUONlist* pdigit;
		    // build the list of fired pads and update the info
		    if (!HitMap[nch]->TestHit(ipx, ipy)) {

			list->AddAtAndExpand(
			    new AliMUONlist(nch,digits),counter);
			
			HitMap[nch]->SetHit(ipx, ipy, counter);
			counter++;
			pdigit=(AliMUONlist*)list->At(list->GetLast());
			// list of tracks
			TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			trlist->Add(&trinfo);
		    } else {
			pdigit=(AliMUONlist*) HitMap[nch]->GetHit(ipx, ipy);
			// update charge
			(*pdigit).fSignal+=iqpad;
			(*pdigit).fPhysics+=iqpad;			
			// update list of tracks
			TObjArray* trlist=(TObjArray*)pdigit->TrackList();
			Int_t last_entry=trlist->GetLast();
			TVector *ptrk_p=(TVector*)trlist->At(last_entry);
			TVector &ptrk=*ptrk_p;
			Int_t last_track=Int_t(ptrk(0));
			Int_t last_charge=Int_t(ptrk(1));
			if (last_track==track) {
			    last_charge+=iqpad;
			    trlist->RemoveAt(last_entry);
			    trinfo(0)=last_track;
			    trinfo(1)=last_charge;
			    trlist->AddAt(&trinfo,last_entry);
			} else {
			    trlist->Add(&trinfo);
			}
			// check the track list
			Int_t nptracks=trlist->GetEntriesFast();
			if (nptracks > 2) {
			    for (Int_t tr=0;tr<nptracks;tr++) {
				TVector *pptrk_p=(TVector*)trlist->At(tr);
				TVector &pptrk=*pptrk_p;
				trk[tr]=Int_t(pptrk(0));
				chtrk[tr]=Int_t(pptrk(1));
			    }
			} // end if nptracks
		    } //  end if pdigit
		} //end loop over clusters
	    } // hit loop
	} // track loop
	
	//Int_t nentr1=list->GetEntriesFast();
	//printf(" \n counter, nentr1 %d %d\n",counter,nentr1);

	// open the file with background
       
	if (Add ) {
	    ntracks =(Int_t)TrH1->GetEntries();
	    //printf("background - icat,ntracks1  %d %d\n",icat,ntracks);
	    //printf("background - Start loop over tracks \n");     
//
//   Loop over tracks
//
	    for (Int_t track=0; track<ntracks; track++) {

		if (fHits2)       fHits2->Clear();
		if (fClusters2)   fClusters2->Clear();

		TrH1->GetEvent(track);
//
//   Loop over hits
		AliMUONhit* mHit;
		for(int i=0;i<fHits2->GetEntriesFast();++i) 
	{	
		    mHit=(AliMUONhit*) (*fHits2)[i];
		    Int_t   nch   = mHit->fChamber-1;  // chamber number
		    if (nch >9) continue;
		    iChamber = &(MUON->Chamber(nch));
		    Int_t rmin = (Int_t)iChamber->RInner();
		    Int_t rmax = (Int_t)iChamber->ROuter();
                    Float_t xbgr=mHit->fX;
		    Float_t ybgr=mHit->fY;
		    Bool_t cond=kFALSE;
		    
		    for (Int_t imuon =0; imuon < nmuon[nch]; imuon++) {
			Float_t dist= (xbgr-xhit[nch][imuon])*(xbgr-xhit[nch][imuon])
			    +(ybgr-yhit[nch][imuon])*(ybgr-yhit[nch][imuon]);
			if (dist<100) cond=kTRUE;
		    }
		    if (!cond) continue;
		    
//
// Loop over pad hits
		    for (AliMUONcluster* mPad=
			     (AliMUONcluster*)MUON->FirstPad(mHit,fClusters2);
			 mPad;
			 mPad=(AliMUONcluster*)MUON->NextPad(fClusters2))
		    {

			Int_t cathode  = mPad->fCathode;    // cathode number
			Int_t ipx      = mPad->fPadX;       // pad number on X
			Int_t ipy      = mPad->fPadY;       // pad number on Y
			Int_t iqpad    = Int_t(mPad->fQpad*kScale);// charge per pad
//			Int_t iqpad    = mPad->fQpad;       // charge per pad

			if (cathode != (icat+1)) continue;
			//if (!HitMap[nch]->CheckBoundary()) continue;
			// fill the info array
			Float_t thex, they;
			segmentation=iChamber->GetSegmentationModel(cathode);
			segmentation->GetPadCxy(ipx,ipy,thex,they);
			Float_t rpad=TMath::Sqrt(thex*thex+they*they);
			if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;

			    new((*p_adr)[countadr++]) TVector(2);
			    TVector &trinfo=*((TVector*) (*p_adr)[countadr-1]);
			    trinfo(0)=-1;  // tag background
			    trinfo(1)=-1;

			digits[0]=ipx;
			digits[1]=ipy;
			digits[2]=iqpad;
			digits[3]=0;
			digits[4]=-1;

			AliMUONlist* pdigit;
			// build the list of fired pads and update the info
			if (!HitMap[nch]->TestHit(ipx, ipy)) {
			    list->AddAtAndExpand(new AliMUONlist(nch,digits),counter);
			
			    HitMap[nch]->SetHit(ipx, ipy, counter);
			    counter++;
			    
			    pdigit=(AliMUONlist*)list->At(list->GetLast());
			    // list of tracks
				TObjArray *trlist=(TObjArray*)pdigit->
				                   TrackList();
				trlist->Add(&trinfo);
			} else {
			    pdigit=(AliMUONlist*) HitMap[nch]->GetHit(ipx, ipy);
			    // update charge
			    (*pdigit).fSignal+=iqpad;

			    // update list of tracks
				TObjArray* trlist=(TObjArray*)pdigit->
				                   TrackList();
				Int_t last_entry=trlist->GetLast();
				TVector *ptrk_p=(TVector*)trlist->
				                 At(last_entry);
				TVector &ptrk=*ptrk_p;
				Int_t last_track=Int_t(ptrk(0));
				if (last_track==-1) {
				    continue;
				} else {
				    trlist->Add(&trinfo);
				}
				// check the track list
				Int_t nptracks=trlist->GetEntriesFast();
				if (nptracks > 0) {
				    for (Int_t tr=0;tr<nptracks;tr++) {
					TVector *pptrk_p=(TVector*)trlist->At(tr);
					TVector &pptrk=*pptrk_p;
					trk[tr]=Int_t(pptrk(0));
					chtrk[tr]=Int_t(pptrk(1));
				    }
				} // end if nptracks
			} //  end if pdigit
		    } //end loop over clusters
		} // hit loop
	    } // track loop
	    //Int_t nentr2=list->GetEntriesFast();
	    //printf(" \n counter2, nentr2 %d %d \n",counter,nentr2);
	    TTree *fAli=gAlice->TreeK();
            TFile *file;
	    
	    if (fAli) file =fAli->GetCurrentFile();
	    file->cd();
	} // if Add	

	Int_t tracks[10];
	Int_t charges[10];
	//cout<<"start filling digits \n "<<endl;
	//	const Float_t zero_supm =    6.;
	Int_t nentries=list->GetEntriesFast();
	//printf(" \n \n nentries %d \n",nentries);
	// start filling the digits
	
	for (Int_t nent=0;nent<nentries;nent++) {
	    AliMUONlist *address=(AliMUONlist*)list->At(nent);
	    if (address==0) continue; 
	    Int_t ich=address->fChamber;
	    Int_t q=address->fSignal; 
	    iChamber=(AliMUONchamber*) (*fChambers)[ich];
	    AliMUONresponse * response=iChamber->GetResponseModel();
	    Int_t adcmax= (Int_t) response->MaxAdc();
	    // add white noise and do zero-suppression and signal truncation
	    Float_t MeanNoise = gRandom->Gaus(1, 0.2);
	    Float_t Noise     = gRandom->Gaus(0, MeanNoise);
	    q+=(Int_t)Noise; 
	    if (address->fPhysics !=0 ) address->fPhysics+=(Int_t)Noise; 
	    if ( q <= zero_supm ) continue;
	    if ( q > adcmax)  q=adcmax;
	    digits[0]=address->fPadX;
	    digits[1]=address->fPadY;
	    digits[2]=q;
	    digits[3]=address->fPhysics;
	    digits[4]=address->fHit;
            //printf("fSignal, fPhysics fTrack %d %d %d \n",digits[2],digits[3],digits[4]);
	    
	    TObjArray* trlist=(TObjArray*)address->TrackList();
	    Int_t nptracks=trlist->GetEntriesFast();
	    //printf("nptracks, trlist   %d  %p\n",nptracks,trlist);

		// this was changed to accomodate the real number of tracks
		if (nptracks > 10) {
		    cout<<"Attention - nptracks > 10 "<<nptracks<<endl;
		    nptracks=10;
		}
		if (nptracks > 2) {
		    printf("Attention - nptracks > 2  %d \n",nptracks);
		    printf("cat,ich,ix,iy,q %d %d %d %d %d \n",icat,ich,digits[0],digits[1],q);
		}
		for (Int_t tr=0;tr<nptracks;tr++) {
		    TVector *pp_p=(TVector*)trlist->At(tr);
		    if(!pp_p ) printf("pp_p - %p\n",pp_p);
		    TVector &pp  =*pp_p;
		    tracks[tr]=Int_t(pp(0));
		    charges[tr]=Int_t(pp(1));
                //printf("tracks, charges - %d %d\n",tracks[tr],charges[tr]);
		}      //end loop over list of tracks for one pad
            // Sort list of tracks according to charge
		if (nptracks > 1) {
		    SortTracks(tracks,charges,nptracks);
		}
		if (nptracks < 10 ) {
		    for (Int_t i=nptracks; i<10; i++) {
			tracks[i]=0;
			charges[i]=0;
		    }
		}

	    // fill digits
	    MUON->AddDigits(ich,tracks,charges,digits);
	}
	//cout<<"I'm out of the loops for digitisation"<<endl;
	gAlice->TreeD()->Fill();
	TTree *TD=gAlice->TreeD();

	Stat_t ndig=TD->GetEntries();
	cout<<"number of digits  "<<ndig<<endl;
	TClonesArray *fDch;
	for (int k=0;k<10;k++) {
	    fDch= MUON->DigitsAddress(k);
	    int ndig=fDch->GetEntriesFast();
	    printf (" i, ndig %d %d \n",k,ndig);
	}

	MUON->ResetDigits();
	list->Delete();
	for(Int_t ii=0;ii<10;++ii) {
	    if (HitMap[ii]) {
		hm=HitMap[ii];
		delete hm;
		HitMap[ii]=0;
	    }
	}
	
    } //end loop over cathodes

       char hname[30];
       sprintf(hname,"TreeD%d",nev);
       gAlice->TreeD()->Write(hname);
       // reset tree
       gAlice->TreeD()->Reset();
       delete list;
       //Int_t nadr=p_adr->GetEntriesFast();
       // printf(" \n \n nadr %d \n",nadr);

       p_adr->Clear();
       // gObjectTable->Print();
       
}

void AliMUON::SortTracks(Int_t *tracks,Int_t *charges,Int_t ntr)
{
  //
  // Sort the list of tracks contributing to a given digit
  // Only the 3 most significant tracks are acctually sorted
  //
  
  //
  //  Loop over signals, only 3 times
  //
  
  Int_t qmax;
  Int_t jmax;
  Int_t idx[3] = {-2,-2,-2};
  Int_t jch[3] = {-2,-2,-2};
  Int_t jtr[3] = {-2,-2,-2};
  Int_t i,j,imax;
  
  if (ntr<3) imax=ntr;
  else imax=3;
  for(i=0;i<imax;i++){
    qmax=0;
    jmax=0;
    
    for(j=0;j<ntr;j++){
      
      if((i == 1 && j == idx[i-1]) 
	 ||(i == 2 && (j == idx[i-1] || j == idx[i-2]))) continue;
      
      if(charges[j] > qmax) {
	qmax = charges[j];
	jmax=j;
      }       
    } 
    
    if(qmax > 0) {
      idx[i]=jmax;
      jch[i]=charges[jmax]; 
      jtr[i]=tracks[jmax]; 
    }
    
  } 
  
  for(i=0;i<3;i++){
    if (jtr[i] == -2) {
         charges[i]=0;
         tracks[i]=0;
    } else {
         charges[i]=jch[i];
         tracks[i]=jtr[i];
    }
  }

}

void AliMUON::FindClusters(Int_t nev,Int_t last_entry)
{

//
// Loop on chambers and on cathode planes
//
  for (Int_t icat=0;icat<2;icat++) {
	    gAlice->ResetDigits();
 	    gAlice->TreeD()->GetEvent(last_entry+icat); // spurious +1 ...
	    if (nev < 10) printf("last_entry , icat - %d %d \n",last_entry,icat);
	    //gAlice->TreeD()->GetEvent(icat+1); // spurious +1 ...

      for (Int_t ich=0;ich<10;ich++) {
	  AliMUONchamber* iChamber=(AliMUONchamber*) (*fChambers)[ich];
	  TClonesArray *MUONdigits  = this->DigitsAddress(ich);
	  if (MUONdigits == 0) continue;
          //
	  // Get ready the current chamber stuff
	  //
	  AliMUONresponse* response = iChamber->GetResponseModel();
	  AliMUONsegmentation*  seg = iChamber->GetSegmentationModel(icat+1);
	  AliMUONClusterFinder* rec = iChamber->GetReconstructionModel();
	  //printf("icat, ich, seg - %d %d %p\n",icat,ich,seg);
	  if (seg) {	  
	      rec->SetSegmentation(seg);
	      rec->SetResponse(response);
	      rec->SetDigits(MUONdigits);
	      rec->SetChamber(ich);
	      if (nev==0) rec->CalibrateCOG(); 
	      rec->FindRawClusters();
	  }  
          //printf("Finish FindRawClusters for cathode %d in chamber %d\n",icat,ich);
	  
          TClonesArray *fRch;
	  fRch=RawClustAddress(ich);
	  fRch->Sort();
          // it seems to work 
         

      } // for ich
      // fill the tree
      TTree *TR=gAlice->TreeR();

      gAlice->TreeR()->Fill();

      Stat_t nent=TR->GetEntries();
      cout<<"number of entries  "<<nent<<endl;
      TClonesArray *fRch;
      for (int i=0;i<10;i++) {
	  fRch=RawClustAddress(i);
	  int nraw=fRch->GetEntriesFast();
	  printf (" i, nraw %d %d \n",i,nraw);
      }
      ResetRawClusters();

  } // for icat

  char hname[30];
  sprintf(hname,"TreeR%d",nev);
  gAlice->TreeR()->Write(hname);
  gAlice->TreeR()->Reset();

  //gObjectTable->Print();

}
 
//______________________________________________________________________________
//_____________________________________________________________________________ 
void AliMUON::CathodeCorrelation(Int_t nev)
{

// Correlates the clusters on the two cathode planes and build a list of
// other possible combinations (potential ghosts) - for the moment use the
// criteria of minimum distance between the CoGs of the two correlated
// clusters


//
// Loop on chambers and on clusters on the cathode plane with the highest
// number of clusters

    static Bool_t first=kTRUE;

     AliMUONRawCluster  *mRaw1;
     AliMUONRawCluster  *mRaw2;
     AliMUONchamber     *iChamber;
     AliMUONsegmentation *seg;
     TArrayF x1, y1, x2, y2, q1, q2;
     x1.Set(5000);
     x2.Set(5000);     
     y1.Set(5000);
     y2.Set(5000);     
     q1.Set(5000);
     q2.Set(5000);     
     
// Get pointers to Alice detectors and Digits containers
     TTree *TR = gAlice->TreeR();
     Int_t nent=(Int_t)TR->GetEntries();
     if (nev < 10) printf("Found %d entries in the tree (must be one per cathode per event! + 1empty)\n",nent);
  

     Int_t idx[4]; 
     Float_t xc2[4],yc2[4];
     Float_t xrec2, yrec2;
     Float_t xd0, xdif, ydif;
     Float_t ysrch,xd,xmax,ymax;
     Int_t ilow, iup, iraw1, i;
     //
     Float_t xarray[50];
     Float_t xdarray[50];
     Float_t yarray[50];
     Float_t qarray[50];
     Int_t idx2[50];

     // Int_t nraw[2], entry,cathode;

     for (i=0;i<50;i++) {
         xdarray[i]=1100.;
         xarray[i]=0.;
         yarray[i]=0.;
         qarray[i]=0.;
         idx2[i]=-1;
     }
     for (i=0;i<4;i++) {
          idx[i]=-1;
          xc2[i]=0.;
          yc2[i]=0.;
     }

     // access to the Raw Clusters tree
     for (Int_t ich=0;ich<10;ich++) {
	 iChamber = &(Chamber(ich));
	 TClonesArray *MUONrawclust  = RawClustAddress(ich);
	 ResetRawClusters();
	 TR->GetEvent(nent-2);
	 //TR->GetEvent(1);
	 Int_t nrawcl1 = MUONrawclust->GetEntries();
	 // printf("Found %d raw clusters for cathode 1 in chamber %d \n"
	 //      ,nrawcl1,ich+1);
         if (!nrawcl1) continue;

	 seg = iChamber->GetSegmentationModel(1);
         // loop over raw clusters of first cathode
	 for (iraw1=0; iraw1<nrawcl1; iraw1++) {
	         mRaw1= (AliMUONRawCluster*)MUONrawclust->UncheckedAt(iraw1);
		 x1[iraw1]=mRaw1->fX;
		 y1[iraw1]=mRaw1->fY;
		 q1[iraw1]=(Float_t)mRaw1->fQ; //maybe better fPeakSignal
	 } // rawclusters cathode 1
//
         // Get information from 2nd cathode
	 ResetRawClusters();
	 TR->GetEvent(nent-1);
	 //TR->GetEvent(2);
	 Int_t nrawcl2 = MUONrawclust->GetEntries();
	 if (!nrawcl2) {
	     for (iraw1=0; iraw1<nrawcl1; iraw1++) {
		 idx[3]=iraw1;
		 xc2[3]=x1[iraw1];
		 yc2[3]=y1[iraw1];
                 //printf("nrawcl2 is zero - idx[0] %d\n",idx[0]);
		 
		 AddCathCorrel(ich,idx,xc2,yc2);
		 // reset
		 idx[3]=-1;
		 xc2[3]=0.;
		 yc2[3]=0.;
		 
	     } // store information from cathode 1 only 
	 } else {
	   //  printf("Found %d raw clusters for cathode 2 in chamber %d \n",
	   // nrawcl2, ich+1);

	     for (Int_t iraw2=0; iraw2<nrawcl2; iraw2++) {
	         mRaw2= (AliMUONRawCluster*)MUONrawclust->UncheckedAt(iraw2);
		 x2[iraw2]=mRaw2->fX;
		 y2[iraw2]=mRaw2->fY;	
		 q2[iraw2]=(Float_t)mRaw2->fQ;	
	     } // rawclusters cathode 2
//
// Initalisation finished
	     for (iraw1=0; iraw1<nrawcl1; iraw1++) {
	     // find the sector
                 Int_t ix,iy;
                 seg->GetPadIxy(x1[iraw1],y1[iraw1],ix,iy);   
                 Int_t isec=seg->Sector(ix,iy);
		 // range to look for ghosts ?!
                 if (ich < 5) {
		     ymax = seg->Dpy(isec)*7/2;
		     xmax = seg->Dpx(isec)*7/2;
                 } else {
		     ymax = seg->Dpy(isec)*13/2;
		     xmax = seg->Dpx(isec)*3/2;
		 }
		 ysrch=ymax+y1[iraw1];
		 
		 ilow = AliMUONRawCluster::
		     BinarySearch(ysrch-2*ymax,y2,0,nrawcl2+1);
		 iup=   AliMUONRawCluster::
		     BinarySearch(ysrch,y2,ilow,nrawcl2+1);
		 if (ilow<0 || iup <0 || iup>nrawcl2) continue;
		 Int_t counter=0;
		 for (Int_t iraw2=ilow; iraw2<=iup; iraw2++) {
		     xrec2=x2[iraw2];
		     yrec2=y2[iraw2];	
		     xdif=x1[iraw1]-xrec2;
		     ydif=y1[iraw1]-yrec2;
		     xd=TMath::Sqrt(xdif*xdif+ydif*ydif);
		     if (iraw2==ilow) { 
			 if (ilow==iup) 
			     xd0=TMath::
			     Sqrt(2*xmax*2*xmax+2*ymax*2*ymax);
			 else xd0=101.; 
		     } 
                     Float_t qdif=TMath::Abs(q1[iraw1]-q2[iraw2])/q1[iraw1];
		     
		     if (x1[iraw1]*xrec2 > 0) {
			 if (xd <= xd0 )  {
//			     printf("q1, q2 qdif % f %f %f \n",q1[iraw1],q2[iraw2],qdif);
//			     printf("x1, x2 y1 y2 % f %f %f %f \n",x1[iraw1],xrec2,y1[iraw1],yrec2);
			   //if (qdif <0.3) { //check this number
				 
				 xd0=xd;
				 idx2[counter]=iraw2;
				 xdarray[counter]=xd;
				 xarray[counter]=xdif;
				 yarray[counter]=ydif;
				 qarray[counter]=qdif;
				 counter++;
			   // }
			     
			 }
		     } // check for same quadrant
                 } // loop over 2nd cathode range 
		 
		 
                 if (counter >=2) {
		     AliMUONRawCluster::
			 SortMin(idx2,xdarray,xarray,yarray,qarray,counter);
		     if (xdarray[0]<seg->Dpx(isec) && xdarray[1]<seg->Dpx(isec)) {
			 if (qarray[0]>qarray[1]){
			     Int_t swap=idx2[0];
			     idx2[0]=idx2[1];
			     idx2[1]=swap;
			 }
		     }
		 }
                 int imax;
                 if (counter <3) imax=counter;
                 else imax=3;

                 for (int i=0;i<imax;i++) {
		     if (idx2[i] >= 0 && idx2[i] < nrawcl2) {
			 if (xarray[i] > xmax || yarray[i] > 2*ymax) 
			     continue;
			 idx[i]=idx2[i];
			 xc2[i]=x2[idx2[i]];
			 yc2[i]=y2[idx2[i]];
		     }
		 }
                 // add info about the cluster on the 'starting' cathode

                 idx[3]=iraw1;
                 xc2[3]=x1[iraw1];
                 yc2[3]=y1[iraw1];
		 //if (idx[0] <0)  printf("iraw1 imax idx2[0] idx[0] %d %d %d %d\n",iraw1,imax,idx2[0],idx[0]);
                 AddCathCorrel(ich,idx,xc2,yc2);
		 // reset
                 for (Int_t ii=0;ii<counter;ii++) {
		     xdarray[ii]=1100.;
		     xarray[ii]=0.;
		     yarray[ii]=0.;
		     qarray[ii]=0.;
		     idx2[ii]=-1;
		 }
                 for (Int_t iii=0;iii<3;iii++) {
		     idx[iii]=-1;
		     xc2[iii]=0.;
		     yc2[iii]=0.;
	         }
	     } // iraw1
	 }
	 x1.Reset();
	 x2.Reset();     
	 y1.Reset();
	 y2.Reset();     
	 q1.Reset();
	 q2.Reset();     
     } //ich
// 
     if (first) {
         MakeTreeC("C");
         first=kFALSE;
     }
     TTree *TC=TreeC();
     TC->Fill();
     //Int_t nentries=(Int_t)TC->GetEntries();
    //cout<<"number entries in tree of correlated clusters  "<<nentries<<endl;
     TClonesArray *fCch;
     static Int_t countev=0;
     Int_t countch=0;

     for (Int_t ii=0;ii<10;ii++) {
	   fCch= CathCorrelAddress(ii);
	   Int_t ncor=fCch->GetEntriesFast();
	   printf (" ii, ncor %d %d \n",ii,ncor);
           if (ncor>=2) countch++;
     }

     // write
     char hname[30];
     sprintf(hname,"TreeC%d",nev);
     TC->Write(hname);
     // reset tree
     ResetCorrelation();
     TC->Reset();

     if (countch==10) countev++;
     printf("countev - %d\n",countev);
    
//     gObjectTable->Print();
     
     
}


//_____________________________________________________________________________

void AliMUON::MakeTreeC(Option_t *option)
{
     char *C = strstr(option,"C");
     if (C && !fTreeC) fTreeC = new TTree("TC","CathodeCorrelation");

//  Create a branch for correlation 

     const Int_t buffersize = 4000;
     char branchname[30];

// one branch for correlation per chamber
     for (int i=0; i<10 ;i++) {
         sprintf(branchname,"%sCorrelation%d",GetName(),i+1);
      
         if (fCathCorrel   && fTreeC) {
	    TreeC()->Branch(branchname,&((*fCathCorrel)[i]), buffersize);
	    printf("Making Branch %s for correlation in chamber %d\n",branchname,i+1);
         }	
     }
}

//_____________________________________________________________________________
void AliMUON::GetTreeC(Int_t event)
{

    // set the branch address
    char treeName[20];
    char branchname[30];

    ResetCorrelation();
    if (fTreeC) {
	  delete fTreeC;
    }

    sprintf(treeName,"TreeC%d",event);
    fTreeC = (TTree*)gDirectory->Get(treeName);


    TBranch *branch;
    if (fTreeC) {
	for (int i=0; i<10; i++) {
	    sprintf(branchname,"%sCorrelation%d",GetName(),i+1);
	    if (fCathCorrel) {
		branch = fTreeC->GetBranch(branchname);
		if (branch) branch->SetAddress(&((*fCathCorrel)[i]));
	    }
	}
    } else {
	printf("ERROR: cannot find CathodeCorrelation Tree for event:%d\n",event);
    }

    // gObjectTable->Print();

}


void AliMUON::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliMUON.
      AliMUONchamber       *iChamber;
      AliMUONsegmentation  *segmentation;
      AliMUONresponse      *response;
      TClonesArray         *digitsaddress;
      TClonesArray         *rawcladdress;
      TClonesArray         *corcladdress;
      //      TObjArray            *clustaddress;
      
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDetector::Streamer(R__b);
      R__b >> fNclusters;
      R__b >> fClusters; // diff
      R__b >> fDchambers;
      R__b >> fRawClusters;
      R__b >> fCathCorrel;
      R__b.ReadArray(fNdch);
      R__b.ReadArray(fNrawch);
      R__b.ReadArray(fNcorch);
      //
      R__b >> fAccCut;
      R__b >> fAccMin;
      R__b >> fAccMax; 
      //   
      // modifs perso  
      R__b >> fSPxzCut;  
      R__b >> fSSigmaCut;
      R__b >> fSXPrec; 
      R__b >> fSYPrec;
      //
      R__b >> fChambers;
// Stream chamber related information
      for (Int_t i =0; i<10; i++) {
	  iChamber=(AliMUONchamber*) (*fChambers)[i];
	  iChamber->Streamer(R__b);
	  if (iChamber->Nsec()==1) {
	      segmentation=iChamber->GetSegmentationModel(1);
	      segmentation->Streamer(R__b);
	  } else {
	      segmentation=iChamber->GetSegmentationModel(1);
	      segmentation->Streamer(R__b);
	      segmentation=iChamber->GetSegmentationModel(2);
	      segmentation->Streamer(R__b);
	  }
          response=iChamber->GetResponseModel();
	  response->Streamer(R__b);	  
	  digitsaddress=(TClonesArray*) (*fDchambers)[i];
	  digitsaddress->Streamer(R__b);
	  rawcladdress=(TClonesArray*) (*fRawClusters)[i];
	  rawcladdress->Streamer(R__b);
	  corcladdress=(TClonesArray*) (*fCathCorrel)[i];
	  corcladdress->Streamer(R__b);
      }
      
   } else {
      R__b.WriteVersion(AliMUON::IsA());
      AliDetector::Streamer(R__b);
      R__b << fNclusters;
      R__b << fClusters; // diff
      R__b << fDchambers;
      R__b << fRawClusters;
      R__b << fCathCorrel;
      R__b.WriteArray(fNdch, 10);
      R__b.WriteArray(fNrawch, 10);
      R__b.WriteArray(fNcorch, 10);
      //
      R__b << fAccCut;
      R__b << fAccMin;
      R__b << fAccMax; 
      //   
      // modifs perso  
      R__b << fSPxzCut;  
      R__b << fSSigmaCut;
      R__b << fSXPrec; 
      R__b << fSYPrec;
      //
      R__b << fChambers;
//  Stream chamber related information
      for (Int_t i =0; i<10; i++) {
	  iChamber=(AliMUONchamber*) (*fChambers)[i];
	  iChamber->Streamer(R__b);
	  if (iChamber->Nsec()==1) {
	      segmentation=iChamber->GetSegmentationModel(1);
	      segmentation->Streamer(R__b);
	  } else {
	      segmentation=iChamber->GetSegmentationModel(1);
	      segmentation->Streamer(R__b);
	      segmentation=iChamber->GetSegmentationModel(2);
	      segmentation->Streamer(R__b);
	  }
          response=iChamber->GetResponseModel();
	  response->Streamer(R__b);
	  digitsaddress=(TClonesArray*) (*fDchambers)[i];
	  digitsaddress->Streamer(R__b);
	  rawcladdress=(TClonesArray*) (*fRawClusters)[i];
	  rawcladdress->Streamer(R__b);
	  corcladdress=(TClonesArray*) (*fCathCorrel)[i];
	  corcladdress->Streamer(R__b);
      }
   }
}
AliMUONcluster* AliMUON::FirstPad(AliMUONhit*  hit, TClonesArray *clusters) 
{
//
    // Initialise the pad iterator
    // Return the address of the first padhit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->fPHlast > 0) {
	sMaxIterPad=hit->fPHlast;
	sCurIterPad=hit->fPHfirst;
	return (AliMUONcluster*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}

AliMUONcluster* AliMUON::NextPad(TClonesArray *clusters) 
{
    sCurIterPad++;
    if (sCurIterPad <= sMaxIterPad) {
	return (AliMUONcluster*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}

//////////////////////////// modifs perso ///////////////

static TTree *ntuple_global;
static TFile *hfile_global;

// variables of the tracking ntuple 
struct { 
  Int_t ievr;           // number of event 
  Int_t ntrackr;        // number of tracks per event
  Int_t istatr[500];    // 1 = good muon, 2 = ghost, 0 = something else
  Int_t isignr[500];    // sign of the track
  Float_t pxr[500];     // x momentum of the reconstructed track
  Float_t pyr[500];     // y momentum of the reconstructed track
  Float_t pzr[500];     // z momentum of the reconstructed track
  Float_t zvr[500];     // z vertex 
  Float_t chi2r[500];   // chi2 of the fit of the track with the field map
  Float_t pxv[500];     // x momentum at vertex
  Float_t pyv[500];     // y momentum at vertex
  Float_t pzv[500];     // z momentum at vertex
} ntuple_st;

AliMUONRawCluster *AliMUON::RawCluster(Int_t ichamber, Int_t icathod, Int_t icluster)
{
    TClonesArray *MUONrawclust  = RawClustAddress(ichamber);
    ResetRawClusters();
    TTree *TR = gAlice->TreeR();
    Int_t nent=(Int_t)TR->GetEntries();
    TR->GetEvent(nent-2+icathod-1);
    //TR->GetEvent(icathod);
    //Int_t nrawcl = (Int_t)MUONrawclust->GetEntriesFast();

    AliMUONRawCluster * mRaw = (AliMUONRawCluster*)MUONrawclust->UncheckedAt(icluster);
    //printf("RawCluster _ nent nrawcl icluster mRaw %d %d %d%p\n",nent,nrawcl,icluster,mRaw);
    
    return  mRaw;
}

void AliMUON::Reconst(Int_t &ifit, Int_t &idebug, Int_t bgd_ev, Int_t &nev, Int_t &idres, Int_t &ireadgeant, Option_t *option,Text_t *filename)
{
  //
  // open kine and hits tree of background file for reconstruction of geant hits 
  // call tracking fortran program
  static Bool_t first=kTRUE;
  static TFile *File;
  char *Add = strstr(option,"Add");
  
  if (Add ) { // only in case of background with geant hits 
    if(first) {
      fFileName=filename;
      cout<<"filename  "<<fFileName<<endl;
      File=new TFile(fFileName);
      cout<<"I have opened "<<fFileName<<" file "<<endl;
      fHits2     = new TClonesArray("AliMUONhit",1000  );
      fParticles2 = new TClonesArray("TParticle",1000);
      first=kFALSE;
    }
    File->cd();
    if(fHits2) fHits2->Clear();
    if(fParticles2) fParticles2->Clear();
    if(TrH1) delete TrH1;
    TrH1=0;
    if(TK1) delete TK1;
    TK1=0;
    // Get Hits Tree header from file
    char treeName[20];
    sprintf(treeName,"TreeH%d",bgd_ev);
    TrH1 = (TTree*)gDirectory->Get(treeName);
    if (!TrH1) {
      printf("ERROR: cannot find Hits Tree for event:%d\n",bgd_ev);
    }
    // set branch addresses
    TBranch *branch;
    char branchname[30];
    sprintf(branchname,"%s",GetName());
    if (TrH1 && fHits2) {
      branch = TrH1->GetBranch(branchname);
      if (branch) branch->SetAddress(&fHits2);
    }
    TrH1->GetEntries();
    // get the Kine tree
    sprintf(treeName,"TreeK%d",bgd_ev);
    TK1 = (TTree*)gDirectory->Get(treeName);
    if (!TK1) {
      printf("ERROR: cannot find Kine Tree for event:%d\n",bgd_ev);
    }
    // set branch addresses
    if (TK1) 
      TK1->SetBranchAddress("Particles", &fParticles2);
    TK1->GetEvent(0);
    
    // get back to the first file
    TTree *TK = gAlice->TreeK();
    TFile *file1 = 0;
    if (TK) file1 = TK->GetCurrentFile();
    file1->cd();
    
  } // end if Add
  
  // call tracking fortran program
  reconstmuon(ifit,idebug,nev,idres,ireadgeant);
}


void AliMUON::InitTracking(Double_t &seff, Double_t &sb0, Double_t &sbl3)
{
  //
  // introduce in fortran program somme parameters and cuts for tracking 
  // create output file "reconst.root" (histos + ntuple)
  cutpxz(fSPxzCut);          // Pxz cut (GeV/c) to begin the track finding
  sigmacut(fSSigmaCut);      // Number of sigmas delimiting the searching areas
  xpreci(fSXPrec);           // Chamber precision in X (cm) 
  ypreci(fSYPrec);           // Chamber precision in Y (cm)
  reco_init(seff,sb0,sbl3);
}

void AliMUON::FinishEvent()
{
    TTree *TK = gAlice->TreeK();
    if (TK) {
      TFile *file1 = TK->GetCurrentFile();
      if(file1) file1->cd();
    }
}

void AliMUON::CloseTracking()
{
  //
  // write histos and ntuple to "reconst.root" file
    reco_term();
}

void chfill(Int_t &id, Float_t &x, Float_t &, Float_t &)
{
  //
  // fill histo like hfill in fortran
    char name[5];
    sprintf(name,"h%d",id);
    TH1F *h1 = (TH1F*) gDirectory->Get(name);
    h1->Fill(x);
}

void chfill2(Int_t &id, Float_t &x, Float_t &y, Float_t &w)
{
  //
  // fill histo like hfill2 in fortran
    char name[5];
    sprintf(name,"h%d",id);
    TH2F *h2 = (TH2F*) gDirectory->Get(name);
    h2->Fill(x,y,w);
}

void chf1(Int_t &id, Float_t &x, Float_t &w)
{
  //
  // fill histo like hf1 in fortran
    char name[5];
    sprintf(name,"h%d",id);
    TH1F *h1 = (TH1F*) gDirectory->Get(name);
    h1->Fill(x,w);
}

void hist_create()
{
  //
  // Create an output file ("reconst.root")
  // Create some histograms and an ntuple

    hfile_global = new TFile("reconst.root","RECREATE","Ntuple - reconstruction");

   ntuple_global = new TTree("ntuple","Reconst ntuple");
   ntuple_global->Branch("ievr",&ntuple_st.ievr,"ievr/I");
   ntuple_global->Branch("ntrackr",&ntuple_st.ntrackr,"ntrackr/I");
   ntuple_global->Branch("istatr",&ntuple_st.istatr[0],"istatr[500]/I");
   ntuple_global->Branch("isignr",&ntuple_st.isignr[0],"isignr[500]/I");
   ntuple_global->Branch("pxr",&ntuple_st.pxr[0],"pxr[500]/F");
   ntuple_global->Branch("pyr",&ntuple_st.pyr[0],"pyr[500]/F");
   ntuple_global->Branch("pzr",&ntuple_st.pzr[0],"pzr[500]/F");
   ntuple_global->Branch("zvr",&ntuple_st.zvr[0],"zvr[500]/F");
   ntuple_global->Branch("chi2r",&ntuple_st.chi2r[0],"chi2r[500]/F");
   ntuple_global->Branch("pxv",&ntuple_st.pxv[0],"pxv[500]/F");
   ntuple_global->Branch("pyv",&ntuple_st.pyv[0],"pyv[500]/F");
   ntuple_global->Branch("pzv",&ntuple_st.pzv[0],"pzv[500]/F");

   // test aliroot

  new TH1F("h100","particule id du hit geant",20,0.,20.);
  new TH1F("h101","position en x du hit geant",100,-200.,200.);
  new TH1F("h102","position en y du hit geant",100,-200.,200.);
  new TH1F("h103","chambre de tracking concernee",15,0.,14.);
  new TH1F("h104","moment ptot du hit geant",50,0.,100.);
  new TH1F("h105","px au vertex",50,0.,20.);
  new TH1F("h106","py au vertex",50,0.,20.);
  new TH1F("h107","pz au vertex",50,0.,20.);
  new TH1F("h108","position zv",50,-15.,15.);
  new TH1F("h109","position en x du hit reconstruit",100,-300.,300.);
  new TH1F("h110","position en y du hit reconstruit",100,-300.,300.);
  new TH1F("h111","delta x ",100,-0.4,0.4);
  new TH1F("h112","delta y ",100,-0.4,0.4);

  char hname[30];
  char hname1[30];
  for (int i=0;i<10;i++) {
    sprintf(hname,"deltax%d",i);
    sprintf(hname1,"h12%d",i);
    new TH1F(hname1,hname ,100,-0.4,0.4);
    sprintf(hname,"deltay%d",i);
    sprintf(hname1,"h13%d",i);
    new TH1F(hname1,hname ,100,-0.4,0.4);
  }
  new TH2F("h2000","VAR X st. 5",30,3.0,183.0,100,0.,25.);
  new TH2F("h2001","VAR Y st. 5",30,3.0,183.0,100,0.,25.);

  new TH2F("h2500","P vs X HHIT",30,3.0,183.0,200,0.,200.);
  new TH2F("h2501","P vs X HHIT**2",30,3.0,183.0,200,0.,5000.);
  new TH2F("h2502","P vs X EPH2 st. 5",30,3.0,183.0,100,0.,0.000005);
  new TH2F("h2503","P vs X EAL2 st. 5",30,3.0,183.0,100,0.,0.01);
  //new TH2F("h2504","P vs X EXM2 st. 5",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2504","P vs X EXM2 st. 5",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2505","P vs X EYM2 st. 5",30,3.0,183.0,100,0.,30.);

  new TH2F("h2507","P vs X EPH st. 5",30,3.0,183.0,100,0.,0.003);
  new TH2F("h2508","P vs X EAL st. 5",30,3.0,183.0,100,0.,0.3);
  //new TH2F("h2509","P vs X EXM st. 5",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2509","P vs X EXM st. 5",30,3.0,183.0,100,0.,0.4);
  new TH2F("h2510","P vs X EYM st. 5",30,3.0,183.0,100,0.,30.);

  new TH2F("h2511","P vs X EPH cut st. 5",30,3.0,183.0,100,0.,0.01);
  new TH2F("h2512","P vs X EAL cut st. 5",30,3.0,183.0,100,0.,0.3);
  //new TH2F("h2513","P vs X EXM cut st. 5",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2513","P vs X EXM cut st. 5",30,3.0,183.0,100,0.,0.4);
  new TH2F("h2514","P vs X EYM cut st. 5",30,3.0,183.0,100,0.,30.);
  // 4
  new TH2F("h2400","P vs X HHIT",30,3.0,183.0,200,0.,200.);
  new TH2F("h2401","P vs X HHIT**2",30,3.0,183.0,200,0.,5000.);
  new TH2F("h2402","P vs X EPH2 st. 4",30,3.0,183.0,100,0.,0.000005);
  new TH2F("h2403","P vs X EAL2 st. 4",30,3.0,183.0,100,0.,0.05);
  //new TH2F("h2404","P vs X EXM2 st. 4",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2404","P vs X EXM2 st. 4",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2405","P vs X EYM2 st. 4",30,3.0,183.0,100,0.,30.);

  new TH2F("h2407","P vs X EPH st. 4",30,3.0,183.0,100,0.,0.003);
  new TH2F("h2408","P vs X EAL st. 4",30,3.0,183.0,100,0.,0.3);
  //new TH2F("h2409","P vs X EXM st. 4",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2409","P vs X EXM st. 4",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2410","P vs X EYM st. 4",30,3.0,183.0,100,0.,30.);

  new TH2F("h2411","P vs X EPH cut st. 4",30,3.0,183.0,100,0.,0.01);
  new TH2F("h2412","P vs X EAL cut st. 4",30,3.0,183.0,100,0.,0.3);
  //new TH2F("h2413","P vs X EXM cut st. 4",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2413","P vs X EXM cut st. 4",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2414","P vs X EYM cut st. 4",30,3.0,183.0,100,0.,30.);
  // 3
  new TH1F("h2301","P2",30,3.0,183.0);
  new TH2F("h2302","P2 vs X EPH2 st. 3",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2303","P2 vs X EAL2 st. 3",30,3.0,183.0,100,0.,0.0005);
  //new TH2F("h2304","P2 vs X EXM2 st. 3",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2304","P2 vs X EXM2 st. 3",30,3.0,183.0,100,0.,2.);
  new TH2F("h2305","P2 vs X EYM2 st. 3",30,3.0,183.0,100,0.,3.);

  new TH2F("h2307","P vs X EPH2 st. 3",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2308","P vs X EAL2 st. 3",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2309","P vs X EXM2 st. 3",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2309","P vs X EXM2 st. 3",30,3.0,183.0,100,0.,2.);
  new TH2F("h2310","P vs X EYM2 st. 3",30,3.0,183.0,100,0.,3.);

  new TH2F("h2311","P vs X EPH cut st. 3",30,3.0,183.0,100,0.,0.06);
  new TH2F("h2312","P vs X EAL cut st. 3",30,3.0,183.0,100,0.,0.05);
  //new TH2F("h2313","P vs X EXM cut st. 3",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2313","P vs X EXM cut st. 3",30,3.0,183.0,100,0.,6.);
  new TH2F("h2314","P vs X EYM cut st. 3",30,3.0,183.0,100,0.,7.);

  new TH2F("h2315","P2 vs X EPH cut st. 3",30,3.0,183.0,100,0.,0.06);
  new TH2F("h2316","P2 vs X EAL cut st. 3",30,3.0,183.0,100,0.,0.05);
  //new TH2F("h2317","P2 vs X EXM cut st. 3",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2317","P2 vs X EXM cut st. 3",30,3.0,183.0,100,0.,6.);
  new TH2F("h2318","P2 vs X EYM cut st. 3",30,3.0,183.0,100,0.,7.);
  
  // 2
  new TH1F("h2201","P2",30,3.0,183.0);
  new TH2F("h2202","P2 vs X EPH2 st. 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2203","P2 vs X EAL2 st. 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2204","P2 vs X EXM2 st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2204","P2 vs X EXM2 st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h2205","P2 vs X EYM2 st. 2",30,3.0,183.0,100,0.,5.);

  new TH2F("h2207","P vs X EPH2 st. 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2208","P vs X EAL2 st. 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2209","P vs X EXM2 st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2209","P vs X EXM2 st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h2210","P vs X EYM2 st. 2",30,3.0,183.0,100,0.,5.);

  new TH2F("h2211","P vs X EPH cut st. 2",30,3.0,183.0,100,0.,0.05);
  new TH2F("h2212","P vs X EAL cut st. 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2213","P vs X EXM cut st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2213","P vs X EXM cut st. 2",30,3.0,183.0,100,0.,11.);
  new TH2F("h2214","P vs X EYM cut st. 2",30,3.0,183.0,100,0.,10.);

  new TH2F("h2215","P2 vs X EPH cut st. 2",30,3.0,183.0,100,0.,0.05);
  new TH2F("h2216","P2 vs X EAL cut st. 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2217","P2 vs X EXM cut st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2217","P2 vs X EXM cut st. 2",30,3.0,183.0,100,0.,11.);
  new TH2F("h2218","P2 vs X EYM cut st. 2",30,3.0,183.0,100,0.,10.);

  // 1
  new TH2F("h2102","P2 vs X EPH2 st. 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2103","P2 vs X EAL2 st. 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2104","P2 vs X EXM2 st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2104","P2 vs X EXM2 st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h2105","P2 vs X EYM2 st. 2",30,3.0,183.0,100,0.,7.);

  new TH2F("h2107","P vs X EPH2 st. 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2108","P vs X EAL2 st. 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2109","P vs X EXM2 st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2109","P vs X EXM2 st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h2110","P vs X EYM2 st. 2",30,3.0,183.0,100,0.,7.);

  new TH2F("h2111","P vs X EPH cut st. 2",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2112","P vs X EAL cut st. 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2113","P vs X EXM cut st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2113","P vs X EXM cut st. 2",30,3.0,183.0,100,0.,11.);
  new TH2F("h2114","P vs X EYM cut st. 2",30,3.0,183.0,100,0.,11.);

  new TH2F("h2115","P2 vs X EPH cut st. 2",30,3.0,183.0,100,0.,0.1);
  new TH2F("h2116","P2 vs X EAL cut st. 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2117","P2 vs X EXM cut st. 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2117","P2 vs X EXM cut st. 2",30,3.0,183.0,100,0.,11.);
  new TH2F("h2118","P2 vs X EYM cut st. 2",30,3.0,183.0,100,0.,11.);

  // 2,3,4,5
  new TH1F("h2701","P2 fit 2",30,3.0,183.0);
  new TH2F("h2702","P2 vs X EPH2 st. 1 fit 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2703","P2 vs X EAL2 st. 1 fit 2",30,3.0,183.0,100,0.,0.005);
  // new TH2F("h2704","P2 vs X EXM2 st. 1 fit 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2704","P2 vs X EXM2 st. 1 fit 2",30,3.0,183.0,100,0.,2.);
  new TH2F("h2705","P2 vs X EYM2 st. 1 fit 2",30,3.0,183.0,100,0.,3.);

  new TH2F("h2707","P vs X EPH2 st. 1 fit 2",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2708","P vs X EAL2 st. 1 fit 2",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2709","P vs X EXM2 st. 1 fit 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2709","P vs X EXM2 st. 1 fit 2",30,3.0,183.0,100,0.,2.);
  new TH2F("h2710","P vs X EYM2 st. 1 fit 2",30,3.0,183.0,100,0.,3.);

  new TH2F("h2711","P vs X EPH cut st. 1 fit 2",30,3.0,183.0,100,0.,0.07);
  new TH2F("h2712","P vs X EAL cut st. 1 fit 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2713","P vs X EXM cut st. 1 fit 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2713","P vs X EXM cut st. 1 fit 2",30,3.0,183.0,100,0.,6.);
  new TH2F("h2714","P vs X EYM cut st. 1 fit 2",30,3.0,183.0,100,0.,7.);

  new TH2F("h2715","P2 vs X EPH cut st. 1 fit 2",30,3.0,183.0,100,0.,0.07);
  new TH2F("h2716","P2 vs X EAL cut st. 1 fit 2",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2717","P2 vs X EXM cut st. 1 fit 2",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2717","P2 vs X EXM cut st. 1 fit 2",30,3.0,183.0,100,0.,6.);
  new TH2F("h2718","P2 vs X EYM cut st. 1 fit 2",30,3.0,183.0,100,0.,7.);

  // 1,3,4,5
  new TH1F("h2801","P2 fit 1",30,3.0,183.0);
  new TH2F("h2802","P2 vs X EPH2 st. 2 fit 1",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2803","P2 vs X EAL2 st. 2 fit 1",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2804","P2 vs X EXM2 st. 2 fit 1",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2804","P2 vs X EXM2 st. 2 fit 1",30,3.0,183.0,100,0.,2.);
  new TH2F("h2805","P2 vs X EYM2 st. 2 fit 1",30,3.0,183.0,100,0.,3.);

  new TH2F("h2807","P vs X EPH2 st. 2 fit 1",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h2808","P vs X EAL2 st. 2 fit 1",30,3.0,183.0,100,0.,0.005);
  //new TH2F("h2809","P vs X EXM2 st. 2 fit 1",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2809","P vs X EXM2 st. 2 fit 1",30,3.0,183.0,100,0.,2.);
  new TH2F("h2810","P vs X EYM2 st. 2 fit 1",30,3.0,183.0,100,0.,3.);

  new TH2F("h2811","P vs X EPH cut st. 2 fit 1",30,3.0,183.0,100,0.,0.05);
  new TH2F("h2812","P vs X EAL cut st. 2 fit 1",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2813","P vs X EXM cut st. 2 fit 1",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2813","P vs X EXM cut st. 2 fit 1",30,3.0,183.0,100,0.,5.);
  new TH2F("h2814","P vs X EYM cut st. 2 fit 1",30,3.0,183.0,100,0.,7.);

  new TH2F("h2815","P2 vs X EPH cut st. 2 fit 1",30,3.0,183.0,100,0.,0.05);
  new TH2F("h2816","P2 vs X EAL cut st. 2 fit 1",30,3.0,183.0,100,0.,0.2);
  //new TH2F("h2817","P2 vs X EXM cut st. 2 fit 1",30,3.0,183.0,100,0.,1.5);
  new TH2F("h2817","P2 vs X EXM cut st. 2 fit 1",30,3.0,183.0,100,0.,5.);
  new TH2F("h2818","P2 vs X EYM cut st. 2 fit 1",30,3.0,183.0,100,0.,7.);
  // fin de test

  new TH1F("h500","Acceptance en H st. 4",500,0.,500.);
  new TH1F("h600","Acceptance en H st. 5",500,0.,500.);
  new TH1F("h700","X vertex track found",200,-10.,10.);
  new TH1F("h701","Y vertex track found",200,-10.,10.);
  new TH1F("h800","Rap. muon gen.",100,0.,5.);
  new TH1F("h801","Rap. muon gen. recons.",100,0.,5.);
  new TH1F("h802","Rap. muon gen. ghost ",100,0.,5.);
  new TH1F("h900","Pt muon gen.",100,0.,20.);
  new TH1F("h901","Pt muon gen. recons.",100,0.,20.);
  new TH1F("h902","Pt muon gen. ghost",100,0.,20.);
  new TH1F("h910","phi muon gen.",100,-10.,10.);
  new TH1F("h911","phi muon gen. recons.",100,-10.,10.);
  new TH1F("h912","phi muon gen. ghost",100,-10.,10.);
  new TH2F("h1001","Y VS X hit st. 1",300,-300.,300.,300,-300.,300.);
  new TH2F("h1002","Y VS X hit st. 2",300,-300.,300.,300,-300.,300.);
  new TH2F("h1003","Y VS X hit st. 3",300,-300.,300.,300,-300.,300.);
  new TH2F("h1004","Y VS X hit st. 4",300,-300.,300.,300,-300.,300.);
  new TH2F("h1005","Y VS X hit st. 5",300,-300.,300.,300,-300.,300.);
  //  Histos variance dans 4      
  new TH2F("h11","VAR X st. 4",30,3.0,183.0,100,0.,2.);
  new TH2F("h12","VAR Y st. 4",30,3.0,183.0,100,0.,600.);
  new TH2F("h13","VAR PHI st. 4",30,3.0,183.0,100,0.,0.0001);
  new TH2F("h14","VAR ALM st. 4",30,3.0,183.0,100,0.,0.05);
  new TH1F("h15","P",30,3.0,183.0);
  new TH1F("h411","VAR X st. 4",100,-1.42,1.42);
  new TH1F("h412","VAR Y st. 4",100,-25.,25.);
  new TH1F("h413","VAR PHI st. 4",100,-0.01,0.01);
  new TH1F("h414","VAR ALM st. 4",100,-0.23,0.23);
  // histo2
  new TH2F("h211","histo2-VAR X st. 4",30,3.0,183.0,100,0.,2.);
  new TH2F("h212","histo2-VAR Y st. 4",30,3.0,183.0,100,0.,600.);
  new TH1F("h213","histo2-VAR X st. 4",100,-1.42,1.42);
  new TH1F("h214","histo2-VAR Y st. 4",100,-25.,25.);
  new TH1F("h215","histo2-P",30,3.0,183.0);

  //  Histos variance dans 2      
  new TH2F("h21","VAR X st. 2",30,3.0,183.0,100,0.,3.);
  new TH2F("h22","VAR Y st. 2",30,3.0,183.0,100,0.,7.);
  new TH2F("h23","VAR PHI st. 2",30,3.0,183.0,100,0.,0.006);
  new TH2F("h24","VAR ALM st. 2",30,3.0,183.0,100,0.,0.005);
  new TH1F("h25","P",30,3.0,183.0);
  new TH1F("h421","VAR X st. 2",100,-1.72,1.72);
  new TH1F("h422","VAR Y st. 2",100,-2.7,2.7);
  new TH1F("h423","VAR PHI st. 2",100,-0.08,0.08);
  new TH1F("h424","VAR ALM st. 2",100,-0.072,0.072);
  // histo2
  new TH2F("h221","histo2-VAR X st. 2",30,3.0,183.0,100,0.,3.);
  new TH2F("h222","histo2-VAR Y st. 2",30,3.0,183.0,100,0.,7.);
  new TH1F("h223","histo2-VAR X st. 2",100,-1.72,1.72);
  new TH1F("h224","histo2-VAR Y st. 2",100,-2.7,2.7);
  new TH1F("h225","histo2-P",30,3.0,183.0);

  //  Histos variance dans 1      
  new TH2F("h31","VAR X st. 1",30,3.0,183.0,100,0.,2.);
  new TH2F("h32","VAR Y st. 1",30,3.0,183.0,100,0.,0.5);
  new TH2F("h33","VAR PHI st. 1",30,3.0,183.0,100,0.,0.006);
  new TH2F("h34","VAR ALM st. 1",30,3.0,183.0,100,0.,0.005);
  new TH1F("h35","P",30,3.0,183.0);
  new TH1F("h431","VAR X st. 1",100,-1.42,1.42);
  new TH1F("h432","VAR Y st. 1",100,-0.72,0.72);
  new TH1F("h433","VAR PHI st. 1",100,-0.08,0.08);
  new TH1F("h434","VAR ALM st. 1",100,-0.072,0.072);
  //  Histos variance dans 1      
  new TH2F("h41","VAR X st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,4.);
  new TH2F("h42","VAR Y st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,20.);
  new TH2F("h43","VAR PHI st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,0.005);
  new TH2F("h44","VAR ALM st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,0.005);
  new TH1F("h45","P",30,3.0,183.0);
  new TH1F("h441","VAR X st. 1 fit 5,4,3,2,V",100,-2.,2.);
  new TH1F("h442","VAR Y st. 1 fit 5,4,3,2,V",100,-4.5,4.5);
  new TH1F("h443","VAR PHI st. 1 fit 5,4,3,2,V",100,-0.072,0.072);
  new TH1F("h444","VAR ALM st. 1 fit 5,4,3,2,V",100,-0.072,0.072);
  // histo2
  new TH2F("h241","histo2-VAR X st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,4.);
  new TH2F("h242","histo2-VAR Y st. 1 fit 5,4,3,2,V",30,3.0,183.0,100,0.,20.);
  new TH1F("h243","histo2-VAR X st. 1 fit 5,4,3,2,V",100,-2.,2.);
  new TH1F("h244","histo2-VAR Y st. 1 fit 5,4,3,2,V",100,-4.5,4.5);
  new TH1F("h245","histo2-P",30,3.0,183.0);

  //  Histos variance dans 2      
  new TH2F("h51","VAR X st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,0.5);
  new TH2F("h52","VAR Y st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,2.);
  new TH2F("h53","VAR PHI st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,0.005);
  new TH2F("h54","VAR ALM st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,0.01);
  new TH1F("h55","P",30,3.0,183.0);
  new TH1F("h451","VAR X st. 2 fit 5,4,3,1,V",100,-0.72,0.72);
  new TH1F("h452","VAR Y st. 2 fit 5,4,3,1,V",100,-1.42,1.42);
  new TH1F("h453","VAR PHI st. 2 fit 5,4,3,1,V",100,-0.072,0.072);
  new TH1F("h454","VAR ALM st. 2 fit 5,4,3,1,V",100,-0.1,0.1);
  new TH1F("h999","PTOT",30,3.0,183.0);
  // histo2
  new TH2F("h251","histo2-VAR X st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,0.5);
  new TH2F("h252","histo2-VAR Y st. 2 fit 5,4,3,1,V",30,3.0,183.0,100,0.,2.);
  new TH1F("h253","histo2-VAR X st. 2 fit 5,4,3,1,V",100,-0.72,0.72);
  new TH1F("h254","histo2-VAR Y st. 2 fit 5,4,3,1,V",100,-1.42,1.42);
  new TH1F("h255","histo2-P",30,3.0,183.0);
  //  Histos variance dans 3      
  new TH2F("h61","VAR X st. 3 fit 4,5,V",30,3.0,183.0,100,0.,5.);
  new TH2F("h62","VAR Y st. 3 fit 4,5,V",30,3.0,183.0,100,0.,2.);
  new TH2F("h63","VAR PHI st. 3 fit 4,5,V",30,3.0,183.0,100,0.,0.0006);
  new TH2F("h64","VAR ALM st. 3 fit 4,5,V",30,3.0,183.0,100,0.,0.0006);
  new TH1F("h65","P",30,3.0,183.0);
  new TH1F("h461","VAR X st. 3 fit 4,5,V",100,-2.25,2.25);
  new TH1F("h462","VAR Y st. 3 fit 4,5,V",100,-1.42,1.42);
  new TH1F("h463","VAR PHI st. 3 fit 4,5,V",100,-0.024,0.024);
  new TH1F("h464","VAR ALM st. 3 fit 4,5,V",100,-0.024,0.024);
  // histo2
  new TH2F("h261","histo2-VAR X st. 3 fit 4,5,V",30,3.0,183.0,100,0.,5.);
  new TH2F("h262","histo2-VAR Y st. 3 fit 4,5,V",30,3.0,183.0,100,0.,2.);
  new TH1F("h263","histo2-VAR X st. 3 fit 4,5,V",100,-2.25,2.25);
  new TH1F("h264","histo2-VAR Y st. 3 fit 4,5,V",100,-1.42,1.42);
  new TH1F("h265","Phisto2-",30,3.0,183.0);
  // Histos dx,dy distribution between chambers inside stations
  new TH1F("h71","DX in st. ID-70",100,-5.,5.);
  new TH1F("h81","DY in st. ID-80",100,-5.,5.);
  new TH1F("h72","DX in st. ID-70",100,-5.,5.);
  new TH1F("h82","DY in st. ID-80",100,-5.,5.);
  new TH1F("h73","DX in st. ID-70",100,-5.,5.);
  new TH1F("h83","DY in st. ID-80",100,-5.,5.);
  new TH1F("h74","DX in st. ID-70",100,-5.,5.);
  new TH1F("h84","DY in st. ID-80",100,-5.,5.);
  new TH1F("h75","DX in st. ID-70",100,-5.,5.);
  new TH1F("h85","DY in st. ID-80",100,-5.,5.);
}

void chfnt(Int_t &ievr, Int_t &ntrackr, Int_t *istatr, Int_t *isignr, Float_t *pxr, Float_t *pyr, Float_t *pzr, Float_t *zvr, Float_t *chi2r,  Float_t *pxv, Float_t *pyv, Float_t *pzv)
{
  //
  // fill the ntuple 
    ntuple_st.ievr = ievr;
    ntuple_st.ntrackr = ntrackr;
    for (Int_t i=0; i<500; i++) {
	ntuple_st.istatr[i] = istatr[i];
	ntuple_st.isignr[i] = isignr[i]; 
	ntuple_st.pxr[i]    = pxr[i]; 
	ntuple_st.pyr[i]    = pyr[i];
	ntuple_st.pzr[i]    = pzr[i];
	ntuple_st.zvr[i]    = zvr[i];
	ntuple_st.chi2r[i]  = chi2r[i];
	ntuple_st.pxv[i]    = pxv[i]; 
	ntuple_st.pyv[i]    = pyv[i];
	ntuple_st.pzv[i]    = pzv[i];
    }
    ntuple_global->Fill();   
}

void hist_closed()
{
  //
  // write histos and ntuple to "reconst.root" file
  hfile_global->Write();
}

void trackf_read_geant(Int_t *itypg, Double_t *xtrg, Double_t *ytrg, Double_t *ptotg, Int_t *idg, Int_t *izch, Double_t *pvert1g, Double_t *pvert2g, Double_t *pvert3g, Double_t *zvertg, Int_t &nhittot1, Double_t *cx, Double_t *cy, Double_t *cz, Int_t &ievr,Int_t &nev,Double_t *xgeant, Double_t *ygeant,Double_t *clsize1, Double_t *clsize2) 
{
  //
  // introduce aliroot variables in fortran common 
  // tracking study from geant hits 
  //

  AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
  
  //  TTree *TK = gAlice->TreeK();
  TTree *TH = gAlice->TreeH();
  Int_t ntracks = (Int_t)TH->GetEntries();
  cout<<"ntrack="<<ntracks<<endl;

  Int_t maxidg = 0;
  Int_t nres=0;
  
//
//  Loop over tracks
//

  for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      TH->GetEvent(track);
      
      if (MUON)  {
//
//  Loop over hits
//
	  for(AliMUONhit* mHit=(AliMUONhit*)MUON->FirstHit(-1); 
	      mHit;
	      mHit=(AliMUONhit*)MUON->NextHit()) 
	  {
	      if (maxidg<=20000) {
		
		if (mHit->fChamber > 10) continue;
		TClonesArray *fPartArray = gAlice->Particles();
		TParticle *Part;
		Int_t ftrack = mHit->fTrack;
		Int_t id = ((TParticle*) fPartArray->UncheckedAt(ftrack))->GetPdgCode();

		if (id==kMuonPlus||id==kMuonMinus) {
		    
		    // inversion de x et y car le champ est inverse dans le programme tracking
		    xtrg[maxidg]   = 0;       
		    ytrg[maxidg]   = 0;       
		    xgeant[maxidg]   = mHit->fY;             // x-pos of hit
		    ygeant[maxidg]   = mHit->fX;             // y-pos of hit
		    clsize1[maxidg]   = 0;     // cluster size on 1-st cathode
		    clsize2[maxidg]   = 0;     // cluster size on 2-nd cathode
		    cx[maxidg]     = mHit->fCyHit;            // Px/P of hit
		    cy[maxidg]     = mHit->fCxHit;            // Py/P of hit
		    cz[maxidg]     = mHit->fCzHit;            // Pz/P of hit
		    izch[maxidg]   = mHit->fChamber; 
		    /*      
		    Int_t pdgtype  = Int_t(mHit->fParticle); // particle number
		    itypg[maxidg]  = gMC->IdFromPDG(pdgtype);

		    */
                    if (id==kMuonPlus) itypg[maxidg]  = 5;
		    else  itypg[maxidg]  = 6;

		    ptotg[maxidg]  = mHit->fPTot;          // P of hit 
		    
		    Part = (TParticle*) fPartArray->UncheckedAt(ftrack);
		    Float_t thet = Part->Theta();
		    thet = thet*180./3.1416;
		    
		    Int_t iparent = Part->GetFirstMother();
		    if (iparent >= 0) {
			Int_t ip;
			while(1) {
			    ip=((TParticle*) fPartArray->UncheckedAt(iparent))->GetFirstMother();
			    if (ip < 0) {
				break;
			    } else {
				iparent = ip;
			    }
			}
		    }
		    //printf("iparent - %d\n",iparent);
		    Int_t id1  = ftrack; // numero de la particule generee au vertex
		    Int_t idum = track+1;
		    Int_t id2 = ((TParticle*) fPartArray->UncheckedAt(iparent))->GetPdgCode();

		    if (id2==443) id2=114;
		    else id2=116;
                   
                    if (id2==116) {
		      nres++;
		    }
		    //printf("id2 %d\n",id2);
		    idg[maxidg] = 30000*id1+10000*idum+id2;
		    
		    pvert1g[maxidg] = Part->Py();      // Px vertex
		    pvert2g[maxidg] = Part->Px();      // Py vertex  
		    pvert3g[maxidg] = Part->Pz();      // Pz vertex
		    zvertg[maxidg]  = Part->Vz();      // z vertex 
		    maxidg ++;

		}
	      }
	  } // hit loop
      } // if MUON
  } // track loop first file

  if (TrH1 && fHits2 ) { // if background file
    ntracks =(Int_t)TrH1->GetEntries();
    printf("Trackf_read - 2-nd file - ntracks %d\n",ntracks);

    //  Loop over tracks
    for (Int_t track=0; track<ntracks; track++) {
      
      if (fHits2) fHits2->Clear();
      TrH1->GetEvent(track);

      //  Loop over hits
      for (int i=0;i<fHits2->GetEntriesFast();i++) 
	{
	  AliMUONhit *mHit=(AliMUONhit*) (*fHits2)[i];
         
          if (mHit->fChamber > 10) continue;

	  if (maxidg<=20000) {
	    
	    // inversion de x et y car le champ est inverse dans le programme tracking !!!!
	    xtrg[maxidg]   = 0;                    // only for reconstructed point
	    ytrg[maxidg]   = 0;                    // only for reconstructed point
	    xgeant[maxidg]   = mHit->fY;           // x-pos of hit
	    ygeant[maxidg]   = mHit->fX;           // y-pos of hit
	    clsize1[maxidg]   = 0;           // cluster size on 1-st cathode
	    clsize2[maxidg]   = 0;           // cluster size on 2-nd cathode
	    cx[maxidg]     = mHit->fCyHit;         // Px/P of hit
	    cy[maxidg]     = mHit->fCxHit;         // Py/P of hit
	    cz[maxidg]     = mHit->fCzHit;         // Pz/P of hit
	    izch[maxidg]   = mHit->fChamber;       // chamber number
	    ptotg[maxidg]  = mHit->fPTot;          // P of hit 
	    
	    Int_t ftrack = mHit->fTrack;
	    Int_t id1  = ftrack;                   // track number 
	    Int_t idum = track+1;
	    
	    TClonesArray *fPartArray = fParticles2;
	    TParticle *Part;
	    Part = (TParticle*) fPartArray->UncheckedAt(ftrack);
	    Int_t id = ((TParticle*) fPartArray->UncheckedAt(ftrack))->GetPdgCode();
	    if (id==kMuonPlus||id==kMuonMinus) {
	        if (id==kMuonPlus) itypg[maxidg]  = 5;
	        else  itypg[maxidg]  = 6;
	    } else itypg[maxidg]=0;
	    
	    Int_t id2=0; // set parent to 0 for background !!
	    idg[maxidg] = 30000*id1+10000*idum+id2;
	    
	    pvert1g[maxidg] = Part->Py();      // Px vertex
	    pvert2g[maxidg] = Part->Px();      // Py vertex  
	    pvert3g[maxidg] = Part->Pz();      // Pz vertex
	    zvertg[maxidg]  = Part->Vz();      // z vertex 

	    maxidg ++;

	  } // check limits (maxidg)
	} // hit loop 
    } // track loop
  } // if TrH1

  ievr = nev;
  nhittot1 = maxidg ;
  cout<<"nhittot1="<<nhittot1<<endl;

  static Int_t nbres=0;
  if (nres>=19) nbres++;
  printf("nres, nbres %d %d \n",nres,nbres);

  hfile_global->cd();      

}



void trackf_read_spoint(Int_t *itypg, Double_t *xtrg, Double_t *ytrg, Double_t *ptotg, Int_t *idg, Int_t *izch, Double_t *pvert1g, Double_t *pvert2g, Double_t *pvert3g, Double_t *zvertg, Int_t &nhittot1, Double_t *cx, Double_t *cy, Double_t *cz, Int_t &ievr,Int_t &nev,Double_t *xgeant, Double_t *ygeant,Double_t *clsize1, Double_t *clsize2) 

{
  //
  // introduce aliroot variables in fortran common 
  // tracking study from reconstructed points 
  //
  AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");

  cout<<"numero de l'evenement "<<nev<<endl;
  
  MUON->GetTreeC(nev);
  TTree *TC=MUON->TreeC();
  TC->GetEntries();

  Int_t maxidg = 0;
  Int_t nres=0;
  Int_t nncor=0;
  static Int_t nuncor=0;
  static Int_t nbadcor=0;
  AliMUONRawCluster * mRaw;
  AliMUONRawCluster * mRaw1;
  TTree *TH = gAlice->TreeH();

  Int_t ihit;
  Int_t mult1, mult2;
  if (MUON) {
      for (Int_t ich=0;ich<10;ich++) {
	  TClonesArray *MUONcorrel  = MUON->CathCorrelAddress(ich);
	  MUON->ResetCorrelation();
	  TC->GetEvent();
	  Int_t ncor = (Int_t)MUONcorrel->GetEntries();
	  if (ncor>=2) nncor++;
	  if (!ncor) continue;

	  //  Loop over correlated clusters
	  for (Int_t icor=0;icor<ncor;icor++) {
	      AliMUONcorrelation * mCor = (AliMUONcorrelation*)MUONcorrel->UncheckedAt(icor);

              Int_t flag=0;   // = 1 if no information in the second cathode
	      Int_t index = mCor->fCorrelIndex[0]; // for the second cathode
              if (index >= 0) {
		  Int_t index1 = mCor->fCorrelIndex[3]; // for the 1-st cathode
		  mRaw1 = MUON->RawCluster(ich,1,index1);
                  mult1=mRaw1->fMultiplicity;
		  mRaw = MUON->RawCluster(ich,2,index);
		  mult2=mRaw->fMultiplicity;
              } else {
		  index = mCor->fCorrelIndex[3];
		  mRaw = MUON->RawCluster(ich,1,index);
		  mult1=mRaw->fMultiplicity;
		  mult2=0;
                  flag=1;
		  nuncor++;
	      }
	      if (!mRaw) continue;

	      Int_t ftrack1 = mRaw->fTracks[1]; // qui doit etre le meme pour 
	                                        // la cathode 1 et 2
              ihit= mRaw->fTracks[0];
	      //printf("icor, ftrack1 ihit %d %d %d\n",icor,ftrack1,ihit);

              if (mRaw->fClusterType == 0 ) {

		  if (maxidg<=20000) {
		      if (flag == 0) {
			  xtrg[maxidg]   = (Double_t) mCor->fY[3];
			  ytrg[maxidg]   = (Double_t) mCor->fX[0];
			  Int_t index1 = mCor->fCorrelIndex[3];
			  mRaw1 = MUON->RawCluster(ich,1,index1);
			  if (mRaw1->fClusterType==1 || mRaw1->fClusterType==2) {
			    Float_t xclust=mCor->fX[3];
			    Float_t yclust=mCor->fY[3];
			    AliMUONchamber *iChamber=&(MUON->Chamber(ich));
			    AliMUONsegmentation *seg = iChamber->GetSegmentationModel(1);
			    Int_t ix,iy;
			    seg->GetPadIxy(xclust,yclust,ix,iy);   
			    Int_t isec=seg->Sector(ix,iy);
			    printf("nev, CORRELATION with pure background in chamber sector %d  %d  %d !!!!!!!!!!!!!!!!!!!!!\n",nev,ich+1,isec);
			    nbadcor++;
			    
			  } // end if cluster type on cathode 1
		      }else {
			  xtrg[maxidg]   = (Double_t) mCor->fY[3];
			  ytrg[maxidg]   = (Double_t) mCor->fX[3];
		      } // if iflag
		      izch[maxidg]   = ich+1;
		      xgeant[maxidg] = 0;
		      ygeant[maxidg] = 0;
		      clsize1[maxidg] = mult1;
		      clsize2[maxidg] = mult2;

		      cx[maxidg]     = 0; // Px/P of hit
		      cy[maxidg]     = 0; // Py/P of hit
		      cz[maxidg]     = 0; // Pz/P of hit
		      itypg[maxidg]  = 0; // particle number
		      ptotg[maxidg]  = 0; // P of hit
		      idg[maxidg]    = 0;
		      pvert1g[maxidg] = 0; // Px vertex
		      pvert2g[maxidg] = 0; // Py vertex  
		      pvert3g[maxidg] = 0; // Pz vertex
		      zvertg[maxidg]  = 0; // z vertex	   
		      maxidg++;
		      
		  }// fin maxidg
		  
	      } else if (mRaw->fClusterType ==1 && ftrack1 < 0) // background + resonance
		{
                  nres++;
		  // get indexmap and loop over digits to find the signal
		  Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
		  gAlice->ResetDigits();
		  if (flag==0) {
		    //gAlice->TreeD()->GetEvent(2); // cathode 2
		    gAlice->TreeD()->GetEvent(nent-1); // cathode 2
		  } else {
		    //gAlice->TreeD()->GetEvent(1); // cathode 1
		    gAlice->TreeD()->GetEvent(nent-2); // cathode 1
		  }

		   TClonesArray *MUONdigits  = MUON->DigitsAddress(ich);
                   Int_t mul=mRaw->fMultiplicity;
                   Int_t trsign;
                   for (int i=0;i<mul;i++) {
		     Int_t idx=mRaw->fIndexMap[i];
                     AliMUONdigit *dig= (AliMUONdigit*)MUONdigits->UncheckedAt(idx);
		     trsign=dig->fTracks[0];
                     ihit=dig->fHit-1;
                     if (trsign > 0 && ihit >= 0) break;

		   } // loop over indexmap

		   //printf("trsign, ihit %d %d\n",trsign,ihit);
		   //printf("signal+background : trsign  %d\n",trsign);
                  
		   if (trsign < 0 || ihit < 0) { // no signal muon  was found
		     
		     if (maxidg<=20000) {
		       if (flag == 0) {
			 xtrg[maxidg]   = (Double_t) mCor->fY[3];
			 ytrg[maxidg]   = (Double_t) mCor->fX[0];
		       }else {
			 xtrg[maxidg]   = (Double_t) mCor->fY[3];
			 ytrg[maxidg]   = (Double_t) mCor->fX[3];
		       }
		       
		       izch[maxidg]   = ich+1;

		      // initialisation of informations which 
		      // can't be reached for background
		       
		       xgeant[maxidg] = 0; // only for resonances
		       ygeant[maxidg] = 0; // only for resonances
		       clsize1[maxidg] = mult1;
		       clsize2[maxidg] = mult2;

		       cx[maxidg]     = 0; // Px/P of hit
		       cy[maxidg]     = 0; // Py/P of hit
		       cz[maxidg]     = 0; // Pz/P of hit
		       itypg[maxidg]  = 0; // particle number
		       ptotg[maxidg]  = 0; // P of hit
		       idg[maxidg]    = 0;
		       pvert1g[maxidg] = 0; // Px vertex
		       pvert2g[maxidg] = 0; // Py vertex  
		       pvert3g[maxidg] = 0; // Pz vertex
		       zvertg[maxidg]  = 0;		   
		       maxidg++;
		       
		     }// fin maxidg
		   } else { // signal muon - retrieve info
		     //printf("inside trsign, ihit %d %d\n",trsign,ihit);
		     if (maxidg<=20000) {
		       if (flag == 0) {
			 xtrg[maxidg]   = (Double_t) mCor->fY[3];
			 ytrg[maxidg]   = (Double_t) mCor->fX[0];
		       }else {
			 xtrg[maxidg]   = (Double_t) mCor->fY[3];
			 ytrg[maxidg]   = (Double_t) mCor->fX[3];
		       }
		       izch[maxidg]   = ich+1;
		       clsize1[maxidg] = mult1;
		       clsize2[maxidg] = mult2;

		      // initialise and set to the correct values 
		      // if signal muons 
		       
		       xgeant[maxidg] = 0; // only for resonances
		       ygeant[maxidg] = 0; // only for resonances
		       
		       cx[maxidg]     = 0; // Px/P of hit
		       cy[maxidg]     = 0; // Py/P of hit
		       cz[maxidg]     = 0; // Pz/P of hit
		       itypg[maxidg]  = 0; // particle number
		       ptotg[maxidg]  = 0; // P of hit
		       idg[maxidg]    = 0;
		       pvert1g[maxidg] = 0; // Px vertex
		       pvert2g[maxidg] = 0; // Py vertex  
		       pvert3g[maxidg] = 0; // Pz vertex
		       zvertg[maxidg]  = 0;	
		       // try to retrieve info about signal muons	   
		       gAlice->ResetHits();
		       TH->GetEvent(trsign);

		       TClonesArray *MUONhits  = MUON->Hits();
		       AliMUONhit *mHit= (AliMUONhit*)MUONhits->
                                                        UncheckedAt(ihit);
			   TClonesArray *fPartArray = gAlice->Particles();
			   TParticle *Part;
			   Int_t nch=mHit->fChamber-1;
			   //printf("sig+bgr ich, nch %d %d \n",ich,nch);
			   if (nch==ich) {
			     Int_t ftrack = mHit->fTrack;
			     Int_t id = ((TParticle*) fPartArray->
                                        UncheckedAt(ftrack))->GetPdgCode();
			     if (id==kMuonPlus||id==kMuonMinus) {
				 xgeant[maxidg] = (Double_t) mHit->fY;
				 ygeant[maxidg] = (Double_t) mHit->fX;
				 cx[maxidg]     = (Double_t) mHit->fCyHit; 
				 cy[maxidg]     = (Double_t) mHit->fCxHit; 
				 cz[maxidg]     = (Double_t) mHit->fCzHit; 

				 if (id==kMuonPlus) {
				   itypg[maxidg]  = 5;
				 } else if (id==kMuonMinus) {
				   itypg[maxidg]  = 6;
				 } else itypg[maxidg]  = 0;
			     
				 ptotg[maxidg]  = (Double_t) mHit->fPTot;  
				 Part = (TParticle*) fPartArray->
				                     UncheckedAt(ftrack);
				 Int_t iparent = Part->GetFirstMother();
				 Int_t id2;
				 id2 = ((TParticle*) fPartArray->
					UncheckedAt(ftrack))->GetPdgCode();
			     
				 if (iparent >= 0) {
				   Int_t ip;
				   while(1) {
				     ip=((TParticle*) fPartArray->
				       UncheckedAt(iparent))->GetFirstMother();
				     if (ip < 0) {
				       id2 = ((TParticle*) fPartArray->
					   UncheckedAt(iparent))->GetPdgCode();
				       break;
				     } else {
				       iparent = ip;
				       id2 = ((TParticle*) fPartArray->
					   UncheckedAt(iparent))->GetPdgCode();
				     } // ip<0
				   } // while
				 }// iparent
				 Int_t id1  = ftrack; 
				 Int_t idum = trsign+1;
			     
				 if (id2==443 || id2==553) {
				   nres++;
				   if (id2==443) id2=114;
				   else id2=116;
				 }
			     
				 idg[maxidg] = 30000*id1+10000*idum+id2;
				 pvert1g[maxidg] = (Double_t) Part->Py(); 
				 pvert2g[maxidg] = (Double_t) Part->Px();   
				 pvert3g[maxidg] = (Double_t) Part->Pz(); 
				 zvertg[maxidg]  = (Double_t) Part->Vz();  
			     } //if muon			     
			   } //if nch
		     maxidg++;
		     } // check limits
		   } // sign+bgr, highest bgr
	      } 
              //pure resonance or mixed cluster with the highest 
	      //contribution coming from resonance
	      if (mRaw->fClusterType >= 1 && ftrack1>=0)  
		{		  	      
		  if (maxidg<=20000) {
		    if (flag == 0) {
		      xtrg[maxidg]   = (Double_t) mCor->fY[3];
		      ytrg[maxidg]   = (Double_t) mCor->fX[0];
		    }else {
		      xtrg[maxidg]   = (Double_t) mCor->fY[3];
		      ytrg[maxidg]   = (Double_t) mCor->fX[3];
		    }
		    clsize1[maxidg] = mult1;
		    clsize2[maxidg] = mult2;
		    izch[maxidg]   = ich+1;

		    Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
		    gAlice->ResetDigits();
		    if (flag==0) {
		      //gAlice->TreeD()->GetEvent(2); // cathode 2
		      gAlice->TreeD()->GetEvent(nent-1); // cathode 2
		    } else {
		      //gAlice->TreeD()->GetEvent(1);        // cathode 1
		      gAlice->TreeD()->GetEvent(nent-2); // cathode 1
		    }

		    TClonesArray *MUONdigits  = MUON->DigitsAddress(ich);
		    Int_t mul=mRaw->fMultiplicity;
		    for (int i=0;i<mul;i++) {
		      Int_t idx=mRaw->fIndexMap[i];
		      AliMUONdigit *dig= (AliMUONdigit*)MUONdigits->UncheckedAt(idx);
		      ihit=dig->fHit-1;
		      if (ihit >= 0) break;

		   } // loop over indexmap
		    //printf("fClusterType, ihit %d %d \n",mRaw->fClusterType,ihit);
		    if (ihit < 0) {
		       xgeant[maxidg] = 0; // only for resonances
		       ygeant[maxidg] = 0; // only for resonances
		       
		       cx[maxidg]     = 0; // Px/P of hit
		       cy[maxidg]     = 0; // Py/P of hit
		       cz[maxidg]     = 0; // Pz/P of hit
		       itypg[maxidg]  = 0; // particle number
		       ptotg[maxidg]  = 0; // P of hit
		       idg[maxidg]    = 0;
		       pvert1g[maxidg] = 0; // Px vertex
		       pvert2g[maxidg] = 0; // Py vertex  
		       pvert3g[maxidg] = 0; // Pz vertex
		       zvertg[maxidg]  = 0;	
		    } else {
		    gAlice->ResetHits();
		    TH->GetEvent(ftrack1);
		    TClonesArray *MUONhits  = MUON->Hits();
		    AliMUONhit *mHit= (AliMUONhit*)MUONhits->
                                                       UncheckedAt(ihit);
			   TClonesArray *fPartArray = gAlice->Particles();
			   TParticle *Part;
			   Int_t nch=mHit->fChamber-1;
			   //printf("signal ich, nch %d %d \n",ich,nch);
			   if (nch==ich) {
			     Int_t ftrack = mHit->fTrack;
			     Int_t id = ((TParticle*) fPartArray->
                                        UncheckedAt(ftrack))->GetPdgCode();
			     //printf("id %d \n",id);
			     if (id==kMuonPlus||id==kMuonMinus) {
				 xgeant[maxidg] = (Double_t) mHit->fY;
				 ygeant[maxidg] = (Double_t) mHit->fX;
				 cx[maxidg]     = (Double_t) mHit->fCyHit; 
				 cy[maxidg]     = (Double_t) mHit->fCxHit; 
     				 cz[maxidg]     = (Double_t) mHit->fCzHit; 

				 if (id==kMuonPlus) {
				   itypg[maxidg]  = 5;
				 } else if (id==kMuonMinus) {
				   itypg[maxidg]  = 6;
				 } else itypg[maxidg]  = 0;
			     
				 ptotg[maxidg]  = (Double_t) mHit->fPTot;  
				 Part = (TParticle*) fPartArray->
				                     UncheckedAt(ftrack);
				 Int_t iparent = Part->GetFirstMother();
				 Int_t id2;
				 id2 = ((TParticle*) fPartArray->
					UncheckedAt(ftrack))->GetPdgCode();
			     
				 if (iparent >= 0) {
				   Int_t ip;
				   while(1) {
				     ip=((TParticle*) fPartArray->
				       UncheckedAt(iparent))->GetFirstMother();
				     if (ip < 0) {
				       id2 = ((TParticle*) fPartArray->
					   UncheckedAt(iparent))->GetPdgCode();
				       break;
				     } else {
				       iparent = ip;
				       id2 = ((TParticle*) fPartArray->
					   UncheckedAt(iparent))->GetPdgCode();
				     } // ip<0
				   } // while
				 }// iparent
				 Int_t id1  = ftrack; 
				 Int_t idum = ftrack1+1;
			     
				 if (id2==443 || id2==553) {
				   nres++;
				   if (id2==443) id2=114;
				   else id2=116;
				 }
				 // printf("id2 %d\n",id2);
				 idg[maxidg] = 30000*id1+10000*idum+id2;
				 pvert1g[maxidg] = (Double_t) Part->Py(); 
				 pvert2g[maxidg] = (Double_t) Part->Px();   
				 pvert3g[maxidg] = (Double_t) Part->Pz(); 
				 zvertg[maxidg]  = (Double_t) Part->Vz();  
			     } //if muon			     
			   } //if nch
		    } // ihit
		    maxidg++;
		  } // check limits
		} // if cluster type
	  } // icor loop
      } // ich loop
  }// if MUON


  ievr = nev;
  cout<<"evenement "<<ievr<<endl;
  nhittot1 = maxidg ;
  cout<<"nhittot1="<<nhittot1<<endl;

  static Int_t nbres=0;
  static Int_t nbcor=0; 
  if (nres>=19) nbres++;
  printf("nres ,nncor - %d %d\n",nres,nncor);
  printf("nbres - %d\n",nbres);
  if (nncor>=20) nbcor++;
  printf("nbcor - %d\n",nbcor);
  printf("nuncor - %d\n",nuncor);
  printf("nbadcor - %d\n",nbadcor);
  
  TC->Reset();

  hfile_global->cd();
  
}

void trackf_fit(Int_t &ivertex, Double_t *pest, Double_t *pstep, Double_t &pxzinv, Double_t &tphi, Double_t &talam, Double_t &xvert, Double_t &yvert)
{
  //
  //  Fit a track candidate with the following input parameters: 
  //  INPUT :  IVERTEX  : vertex flag, if IVERTEX=1 (XVERT,YVERT) are free paramaters
  //                                   if IVERTEX=1 (XVERT,YVERT)=(0.,0.) 
  //           PEST(5)  : starting value of parameters (minuit)
  //           PSTEP(5) : step size for parameters (minuit)
  //  OUTPUT : PXZINV,TPHI,TALAM,XVERT,YVERT : fitted value of the parameters

  static Double_t arglist[10];
  static Double_t c[5] = {0.4, 0.45, 0.45, 90., 90.};
  static Double_t b1, b2, epxz, efi, exs, exvert, eyvert;
  TString chname;
  Int_t ierflg = 0;
  
  TMinuit *gMinuit = new TMinuit(5);
  gMinuit->mninit(5,10,7);
  gMinuit->SetFCN(fcnfwrap);  // constant m.f.

  arglist[0] = -1;
  
  gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);
  //      gMinuit->mnseti('track fitting');
  
  gMinuit->mnparm(0, "invmom",  pest[0], pstep[0], -c[0], c[0], ierflg);
  gMinuit->mnparm(1, "azimuth", pest[1], pstep[1], -c[1], c[1], ierflg);
  gMinuit->mnparm(2, "deep",    pest[2], pstep[2], -c[2], c[2], ierflg);
  if (ivertex==1) {
    gMinuit->mnparm(3, "x ", pest[3], pstep[3], -c[3], c[3], ierflg);
    gMinuit->mnparm(4, "y ", pest[4], pstep[4], -c[4], c[4], ierflg);  
  }   
  
  gMinuit->mnexcm("SET NOGR", arglist, 0, ierflg);
  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
  gMinuit->mnexcm("EXIT" , arglist, 0, ierflg);
  
  gMinuit->mnpout(0, chname, pxzinv, epxz, b1, b2, ierflg);
  gMinuit->mnpout(1, chname, tphi, efi, b1, b2, ierflg);
  gMinuit->mnpout(2, chname, talam, exs, b1, b2, ierflg);
  if (ivertex==1) {
    gMinuit->mnpout(3, chname, xvert, exvert, b1, b2, ierflg);
    gMinuit->mnpout(4, chname, yvert, eyvert, b1, b2, ierflg);
  }   
  
  delete gMinuit;

}
	   
void fcnf(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *pest, Int_t iflag)
{
  //
  // function called by trackf_fit
      Int_t futil = 0;
      fcn(npar,grad,fval,pest,iflag,futil);
}

void prec_fit(Double_t &pxzinv, Double_t &fis, Double_t &alams, Double_t &xvert, Double_t &yvert, Double_t &pxzinvf, Double_t &fif, Double_t &alf, Double_t &xvertf, Double_t &yvertf, Double_t &epxzinv, Double_t &efi, Double_t &exs, Double_t &exvert, Double_t &eyvert)
{
  // 
  // minuit fits for tracking finding 
                                                                            
      static Double_t arglist[10];
      static Double_t c1[5] = {0.001, 0.001, 0.001, 1., 1.};
      static Double_t c2[5] = {0.5, 0.5, 0.5, 120., 120.};
      static Double_t emat[9];
      static Double_t b1, b2;
      Double_t fmin, fedm, errdef; 
      Int_t npari, nparx, istat;

      TString chname;
      Int_t ierflg = 0;
      
      TMinuit *gMinuit = new TMinuit(5);
      gMinuit->mninit(5,10,7);
      gMinuit->SetFCN(fcnfitfwrap);

      arglist[0] = -1.;
      gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);
      
      //      gMinuit->mnseti('track fitting');

      gMinuit->mnparm(0,"invmom",   pxzinv, c1[0], -c2[0], c2[0], ierflg); // 0.003, 0.5
      gMinuit->mnparm(1,"azimuth ", fis,    c1[1], -c2[1], c2[1], ierflg);
      gMinuit->mnparm(2,"deep    ", alams,  c1[2], -c2[2], c2[2], ierflg);
      gMinuit->mnparm(3,"xvert",    xvert,  c1[3], -c2[3], c2[3], ierflg);
      gMinuit->mnparm(4,"yvert",    yvert,  c1[4], -c2[4], c2[4], ierflg);

      gMinuit->mnexcm("SET NOGR", arglist, 0, ierflg);
      arglist[0] = 2.;
      gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
      gMinuit->mnexcm("EXIT", arglist, 0, ierflg);
 
      gMinuit->mnpout(0, chname, pxzinvf, epxzinv, b1, b2, ierflg);
      gMinuit->mnpout(1, chname, fif, efi, b1, b2, ierflg);
      gMinuit->mnpout(2, chname, alf, exs, b1, b2, ierflg);
      gMinuit->mnpout(3, chname, xvertf, exvert, b1, b2, ierflg);
      gMinuit->mnpout(4, chname, yvertf, eyvert, b1, b2, ierflg);
  
      gMinuit->mnemat(emat, 3);
      gMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);

     delete gMinuit;
}

void fcnfitf(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *xval, Int_t iflag)
{
  //
  // function called by prec_fit 
      Int_t futil = 0;
      fcnfit(npar,grad,fval,xval,iflag,futil);
}

///////////////////// fin modifs perso //////////////////////

ClassImp(AliMUONcluster)
 
//___________________________________________
AliMUONcluster::AliMUONcluster(Int_t *clhits)
{
   fHitNumber=clhits[0];
   fCathode=clhits[1];
   fQ=clhits[2];
   fPadX=clhits[3];
   fPadY=clhits[4];
   fQpad=clhits[5];
   fRSec=clhits[6];
}
ClassImp(AliMUONdigit)
//_____________________________________________________________________________
AliMUONdigit::AliMUONdigit(Int_t *digits)
{
  //
  // Creates a MUON digit object to be updated
  //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fSignal      = digits[2];
    fPhysics     = digits[3];
    fHit       = digits[4];

}
//_____________________________________________________________________________
AliMUONdigit::AliMUONdigit(Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Creates a MUON digit object
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fSignal      = digits[2];
    fPhysics     = digits[3];
    fHit       = digits[4];
    for(Int_t i=0; i<10; i++) {
	fTcharges[i]  = charges[i];
	fTracks[i]    = tracks[i];
    }
}

AliMUONdigit::~AliMUONdigit()
{
    
}

ClassImp(AliMUONlist)
 
//____________________________________________________________________________
    AliMUONlist::AliMUONlist(Int_t ich, Int_t *digits): 
	AliMUONdigit(digits)
{
    //
    // Creates a MUON digit list object
    //

    fChamber     = ich;
    fTrackList   = new TObjArray;
 
}

ClassImp(AliMUONhit)
 
//___________________________________________
    AliMUONhit::AliMUONhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
	AliHit(shunt, track)
{
    fChamber=vol[0];
    fParticle=hits[0];
    fX=hits[1];
    fY=hits[2];
    fZ=hits[3];
    fTheta=hits[4];
    fPhi=hits[5];
    fTlength=hits[6];
    fEloss=hits[7];
    fPHfirst=(Int_t) hits[8];
    fPHlast=(Int_t) hits[9];

    // modifs perso
    fPTot=hits[10];
    fCxHit=hits[11];
    fCyHit=hits[12];
    fCzHit=hits[13];
}
ClassImp(AliMUONcorrelation)
//___________________________________________
//_____________________________________________________________________________
AliMUONcorrelation::AliMUONcorrelation(Int_t *idx, Float_t *x, Float_t *y)
{
    //
    // Creates a MUON correlation object
    //
    for(Int_t i=0; i<4; i++) {
	fCorrelIndex[i]  = idx[i];
	fX[i]    = x[i];
	fY[i]    = y[i];
    }
}
ClassImp(AliMUONRawCluster)
Int_t AliMUONRawCluster::Compare(TObject *obj)
{
  /*
         AliMUONRawCluster *raw=(AliMUONRawCluster *)obj;
	 Float_t r=GetRadius();
         Float_t ro=raw->GetRadius();
         if (r>ro) return 1;
         else if (r<ro) return -1;
         else return 0;
  */
         AliMUONRawCluster *raw=(AliMUONRawCluster *)obj;
	 Float_t y=fY;
         Float_t yo=raw->fY;
         if (y>yo) return 1;
         else if (y<yo) return -1;
         else return 0;

}

Int_t AliMUONRawCluster::
BinarySearch(Float_t y, TArrayF coord, Int_t from, Int_t upto)
{
   // Find object using a binary search. Array must first have been sorted.
   // Search can be limited by setting upto to desired index.

   Int_t low=from, high=upto-1, half;
   while(high-low>1) {
        half=(high+low)/2;
        if(y>coord[half]) low=half;
        else high=half;
   }
   return low;
}

void AliMUONRawCluster::SortMin(Int_t *idx,Float_t *xdarray,Float_t *xarray,Float_t *yarray,Float_t *qarray, Int_t ntr)
{
  //
  // Get the 3 closest points(cog) one can find on the second cathode 
  // starting from a given cog on first cathode
  //
  
  //
  //  Loop over deltax, only 3 times
  //
  
    Float_t xmin;
    Int_t jmin;
    Int_t id[3] = {-2,-2,-2};
    Float_t jx[3] = {0.,0.,0.};
    Float_t jy[3] = {0.,0.,0.};
    Float_t jq[3] = {0.,0.,0.};
    Int_t jid[3] = {-2,-2,-2};
    Int_t i,j,imax;
  
    if (ntr<3) imax=ntr;
    else imax=3;
    for(i=0;i<imax;i++){
        xmin=1001.;
        jmin=0;
    
        for(j=0;j<ntr;j++){
            if ((i == 1 && j == id[i-1]) 
	          ||(i == 2 && (j == id[i-1] || j == id[i-2]))) continue;
           if (TMath::Abs(xdarray[j]) < xmin) {
	      xmin = TMath::Abs(xdarray[j]);
	      jmin=j;
           }       
        } // j
        if (xmin != 1001.) {    
           id[i]=jmin;
           jx[i]=xarray[jmin]; 
           jy[i]=yarray[jmin]; 
           jq[i]=qarray[jmin]; 
           jid[i]=idx[jmin];
        } 
    
    }  // i
  
    for (i=0;i<3;i++){
        if (jid[i] == -2) {
            xarray[i]=1001.;
            yarray[i]=1001.;
            qarray[i]=1001.;
            idx[i]=-1;
        } else {
            xarray[i]=jx[i];
            yarray[i]=jy[i];
            qarray[i]=jq[i];
            idx[i]=jid[i];
        }
    }

}


Int_t AliMUONRawCluster::PhysicsContribution()
{
  Int_t iPhys=0;
  Int_t iBg=0;
  Int_t iMixed=0;
  for (Int_t i=0; i<fMultiplicity; i++) {
    if (fPhysicsMap[i]==2) iPhys++;
    if (fPhysicsMap[i]==1) iMixed++;
    if (fPhysicsMap[i]==0) iBg++;
  }
  if (iMixed==0 && iBg==0) {
    return 2;
  } else if ((iPhys != 0 && iBg !=0) || iMixed != 0) {
    return 1;
  } else {
    return 0;
  }
}

   
ClassImp(AliMUONreccluster) 
ClassImp(AliMUONsegmentation)
ClassImp(AliMUONresponse)	



















