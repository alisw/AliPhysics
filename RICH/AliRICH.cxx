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
  Revision 1.18  2000/06/12 15:15:46  jbarbosa
  Cleaned up version.

  Revision 1.17  2000/06/09 14:58:37  jbarbosa
  New digitisation per particle type

  Revision 1.16  2000/04/19 12:55:43  morsch
  Newly structured and updated version (JB, AM)

*/


////////////////////////////////////////////////
//  Manager and hits classes for set:RICH     //
////////////////////////////////////////////////

#include <TBRIK.h>
#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 
#include <TObject.h>
#include <TVector.h>
#include <TObjArray.h>
#include <TArrayF.h>
#include <TFile.h>
#include <TParticle.h>
#include <iostream.h>

#include "AliRICH.h"
#include "AliRICHSegmentation.h"
#include "AliRICHHit.h"
#include "AliRICHCerenkov.h"
#include "AliRICHPadHit.h"
#include "AliRICHDigit.h"
#include "AliRICHTransientDigit.h"
#include "AliRICHRawCluster.h"
#include "AliRICHRecHit.h"
#include "AliRICHHitMapA1.h"
#include "AliRICHClusterFinder.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliPoints.h"
#include "AliCallf77.h" 


// Static variables for the pad-hit iterator routines
static Int_t sMaxIterPad=0;
static Int_t sCurIterPad=0;
static TClonesArray *fClusters2;
static TClonesArray *fHits2;
static TTree *TrH1;
 
ClassImp(AliRICH)
    
//___________________________________________
AliRICH::AliRICH()
{
// Default constructor for RICH manager class

    fIshunt     = 0;
    fHits       = 0;
    fPadHits    = 0;
    fNPadHits   = 0;
    fNcerenkovs = 0;
    fDchambers  = 0;
    fCerenkovs  = 0;
    fNdch       = 0;
}

//___________________________________________
AliRICH::AliRICH(const char *name, const char *title)
    : AliDetector(name,title)
{
//Begin_Html
/*
  <img src="gif/alirich.gif">
*/
//End_Html
    
    fHits       = new TClonesArray("AliRICHHit",1000  );
    gAlice->AddHitList(fHits);
    fPadHits    = new TClonesArray("AliRICHPadHit",100000);
    fCerenkovs  = new TClonesArray("AliRICHCerenkov",1000);
    gAlice->AddHitList(fCerenkovs);
    //gAlice->AddHitList(fHits);
    fNPadHits   = 0;
    fNcerenkovs = 0;
    fIshunt     = 0;
    
    fNdch      = new Int_t[kNCH];
    
    fDchambers = new TObjArray(kNCH);

    fRecHits = new TObjArray(kNCH);
    
    Int_t i;
   
    for (i=0; i<kNCH ;i++) {
	(*fDchambers)[i] = new TClonesArray("AliRICHDigit",10000); 
	fNdch[i]=0;
    }

    fNrawch      = new Int_t[kNCH];
    
    fRawClusters = new TObjArray(kNCH);
    //printf("Created fRwClusters with adress:%p",fRawClusters);

    for (i=0; i<kNCH ;i++) {
      (*fRawClusters)[i] = new TClonesArray("AliRICHRawCluster",10000); 
      fNrawch[i]=0;
    }

    fNrechits      = new Int_t[kNCH];
    
    for (i=0; i<kNCH ;i++) {
	(*fRecHits)[i] = new TClonesArray("AliRICHRecHit",1000); 
    }
    //printf("Created fRecHits with adress:%p",fRecHits);

        
    SetMarkerColor(kRed);
}

AliRICH::AliRICH(const AliRICH& RICH)
{
// Copy Constructor
}


//___________________________________________
AliRICH::~AliRICH()
{

// Destructor of RICH manager class

    fIshunt  = 0;
    delete fHits;
    delete fPadHits;
    delete fCerenkovs;
}

//___________________________________________
void AliRICH::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{

//  
// Adds a hit to the Hits list
//

    TClonesArray &lhits = *fHits;
    new(lhits[fNhits++]) AliRICHHit(fIshunt,track,vol,hits);
}
//_____________________________________________________________________________
void AliRICH::AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs)
{

//
// Adds a RICH cerenkov hit to the Cerenkov Hits list
//

    TClonesArray &lcerenkovs = *fCerenkovs;
    new(lcerenkovs[fNcerenkovs++]) AliRICHCerenkov(fIshunt,track,vol,cerenkovs);
    //printf ("Done for Cerenkov %d\n\n\n\n",fNcerenkovs);
}
//___________________________________________
void AliRICH::AddPadHit(Int_t *clhits)
{

//
// Add a RICH pad hit to the list
//

    TClonesArray &lPadHits = *fPadHits;
    new(lPadHits[fNPadHits++]) AliRICHPadHit(clhits);
} 
//_____________________________________________________________________________
void AliRICH::AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{

  //
  // Add a RICH digit to the list
  //

    TClonesArray &ldigits = *((TClonesArray*)(*fDchambers)[id]);
    new(ldigits[fNdch[id]++]) AliRICHDigit(tracks,charges,digits);
}

//_____________________________________________________________________________
void AliRICH::AddRawCluster(Int_t id, const AliRICHRawCluster& c)
{
    //
    // Add a RICH digit to the list
    //

    TClonesArray &lrawcl = *((TClonesArray*)(*fRawClusters)[id]);
    new(lrawcl[fNrawch[id]++]) AliRICHRawCluster(c);
}

//_____________________________________________________________________________
void AliRICH::AddRecHit(Int_t id, Float_t *rechit, Float_t *photons, Int_t *padsx, Int_t* padsy)
{
  
  //
  // Add a RICH reconstructed hit to the list
  //

    TClonesArray &lrec = *((TClonesArray*)(*fRecHits)[id]);
    new(lrec[fNrechits[id]++]) AliRICHRecHit(id,rechit,photons,padsx,padsy);
}

//___________________________________________
void AliRICH::BuildGeometry()
    
{
  
  //
  // Builds a TNode geometry for event display
  //
    TNode *node, *top;
    
    const int kColorRICH = kGreen;
    //
    top=gAlice->GetGeometry()->GetNode("alice");
    
    
    new TBRIK("S_RICH","S_RICH","void",71.09999,11.5,73.15);
    
    top->cd();
    Float_t pos1[3]={0,471.8999,165.2599};
    //Chamber(0).SetChamberTransform(pos1[0],pos1[1],pos1[2],
    new TRotMatrix("rot993","rot993",90,0,70.69,90,19.30999,-90);
    node = new TNode("RICH1","RICH1","S_RICH",pos1[0],pos1[1],pos1[2],"rot993");
    

    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    
    Float_t pos2[3]={171,470,0};
    //Chamber(1).SetChamberTransform(pos2[0],pos2[1],pos2[2],
    new TRotMatrix("rot994","rot994",90,-20,90,70,0,0);
    node = new TNode("RICH2","RICH2","S_RICH",pos2[0],pos2[1],pos2[2],"rot994");
    
    
    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos3[3]={0,500,0};
    //Chamber(2).SetChamberTransform(pos3[0],pos3[1],pos3[2],
    new TRotMatrix("rot995","rot995",90,0,90,90,0,0);
    node = new TNode("RICH3","RICH3","S_RICH",pos3[0],pos3[1],pos3[2],"rot995");
    

    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos4[3]={-171,470,0};
    //Chamber(3).SetChamberTransform(pos4[0],pos4[1],pos4[2], 
    new TRotMatrix("rot996","rot996",90,20,90,110,0,0);  
    node = new TNode("RICH4","RICH4","S_RICH",pos4[0],pos4[1],pos4[2],"rot996");
    

    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos5[3]={161.3999,443.3999,-165.3};
    //Chamber(4).SetChamberTransform(pos5[0],pos5[1],pos5[2],
    new TRotMatrix("rot997","rot997",90,340,108.1999,70,18.2,70);
    node = new TNode("RICH5","RICH5","S_RICH",pos5[0],pos5[1],pos5[2],"rot997");
    
    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos6[3]={0., 471.9, -165.3,};
    //Chamber(5).SetChamberTransform(pos6[0],pos6[1],pos6[2],
    new TRotMatrix("rot998","rot998",90,0,109.3099,90,19.30999,90);
    node = new TNode("RICH6","RICH6","S_RICH",pos6[0],pos6[1],pos6[2],"rot998");
    
    
    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos7[3]={-161.399,443.3999,-165.3};
    //Chamber(6).SetChamberTransform(pos7[0],pos7[1],pos7[2],
    new TRotMatrix("rot999","rot999",90,20,108.1999,110,18.2,110);
    node = new TNode("RICH7","RICH7","S_RICH",pos7[0],pos7[1],pos7[2],"rot999");
    node->SetLineColor(kColorRICH);
    fNodes->Add(node); 
    
}

//___________________________________________
Int_t AliRICH::DistancetoPrimitive(Int_t , Int_t )
{

// Default value

    return 9999;
}

//___________________________________________
void AliRICH::MakeBranch(Option_t* option)
{
  // Create Tree branches for the RICH.
    
    const Int_t kBufferSize = 4000;
    char branchname[20];
    
    
    AliDetector::MakeBranch(option);
    sprintf(branchname,"%sCerenkov",GetName());
    if (fCerenkovs   && gAlice->TreeH()) {
	gAlice->TreeH()->Branch(branchname,&fCerenkovs, kBufferSize);
	printf("Making Branch %s for Cerenkov Hits\n",branchname);
    }
    
    sprintf(branchname,"%sPadHits",GetName());
    if (fPadHits   && gAlice->TreeH()) {
	gAlice->TreeH()->Branch(branchname,&fPadHits, kBufferSize);
	printf("Making Branch %s for PadHits\n",branchname);
    }
    
// one branch for digits per chamber
    Int_t i;
    
    for (i=0; i<kNCH ;i++) {
	sprintf(branchname,"%sDigits%d",GetName(),i+1);
	
	if (fDchambers   && gAlice->TreeD()) {
	    gAlice->TreeD()->Branch(branchname,&((*fDchambers)[i]), kBufferSize);
	    printf("Making Branch %s for digits in chamber %d\n",branchname,i+1);
	}	
    }

// one branch for raw clusters per chamber
  for (i=0; i<kNCH ;i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
      
      if (fRawClusters   && gAlice->TreeR()) {
	 gAlice->TreeR()->Branch(branchname,&((*fRawClusters)[i]), kBufferSize);
	 printf("Making Branch %s for raw clusters in chamber %d\n",branchname,i+1);
      }	
  }

  // one branch for rec hits per chamber
  for (i=0; i<kNCH ;i++) {
    sprintf(branchname,"%sRecHits%d",GetName(),i+1);
    
    if (fRecHits   && gAlice->TreeR()) {
      gAlice->TreeR()->Branch(branchname,&((*fRecHits)[i]), kBufferSize);
      printf("Making Branch %s for rec. hits in chamber %d\n",branchname,i+1);
    }	
  }
}

//___________________________________________
void AliRICH::SetTreeAddress()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[20];
  Int_t i;

    AliDetector::SetTreeAddress();
    
    TBranch *branch;
    TTree *treeH = gAlice->TreeH();
    TTree *treeD = gAlice->TreeD();
    TTree *treeR = gAlice->TreeR();
    
    if (treeH) {
	if (fPadHits) {
	    branch = treeH->GetBranch("RICHPadHits");
	    if (branch) branch->SetAddress(&fPadHits);
	}
	if (fCerenkovs) {
	    branch = treeH->GetBranch("RICHCerenkov");
	    if (branch) branch->SetAddress(&fCerenkovs);
	}
    }
    
    if (treeD) {
	for (int i=0; i<kNCH; i++) {
	    sprintf(branchname,"%sDigits%d",GetName(),i+1);
	    if (fDchambers) {
		branch = treeD->GetBranch(branchname);
		if (branch) branch->SetAddress(&((*fDchambers)[i]));
	    }
	}
    }
  if (treeR) {
      for (i=0; i<kNCH; i++) {
	  sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	  if (fRawClusters) {
	      branch = treeR->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	  }
      }
      
      for (i=0; i<kNCH; i++) {
	sprintf(branchname,"%sRecHits%d",GetName(),i+1);
	if (fRecHits) {
	  branch = treeR->GetBranch(branchname);
	  if (branch) branch->SetAddress(&((*fRecHits)[i]));
	  }
      }
      
  }
}
//___________________________________________
void AliRICH::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
    AliDetector::ResetHits();
    fNPadHits   = 0;
    fNcerenkovs = 0;
    if (fPadHits)  fPadHits->Clear();
    if (fCerenkovs) fCerenkovs->Clear();
}


//____________________________________________
void AliRICH::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  //
    for ( int i=0;i<kNCH;i++ ) {
	if ((*fDchambers)[i])   (*fDchambers)[i]->Clear();
	if (fNdch)  fNdch[i]=0;
    }
}

//____________________________________________
void AliRICH::ResetRawClusters()
{
  //
  // Reset number of raw clusters and the raw clust array for this detector
  //
    for ( int i=0;i<kNCH;i++ ) {
	if ((*fRawClusters)[i])    ((TClonesArray*)(*fRawClusters)[i])->Clear();
	if (fNrawch)  fNrawch[i]=0;
    }
}

//____________________________________________
void AliRICH::ResetRecHits()
{
  //
  // Reset number of raw clusters and the raw clust array for this detector
  //
  
  for ( int i=0;i<kNCH;i++ ) {
	if ((*fRecHits)[i])    ((TClonesArray*)(*fRecHits)[i])->Clear();
	if (fNrechits)  fNrechits[i]=0;
    }
}

//___________________________________________
void   AliRICH::SetGeometryModel(Int_t id, AliRICHGeometry *geometry)
{

//
// Setter for the RICH geometry model
//


    ((AliRICHChamber*) (*fChambers)[id])->GeometryModel(geometry);
}

//___________________________________________
void   AliRICH::SetSegmentationModel(Int_t id, AliRICHSegmentation *segmentation)
{

//
// Setter for the RICH segmentation model
//

    ((AliRICHChamber*) (*fChambers)[id])->SegmentationModel(segmentation);
}

//___________________________________________
void   AliRICH::SetResponseModel(Int_t id, AliRICHResponse *response)
{

//
// Setter for the RICH response model
//

    ((AliRICHChamber*) (*fChambers)[id])->ResponseModel(response);
}

void   AliRICH::SetReconstructionModel(Int_t id, AliRICHClusterFinder *reconst)
{

//
// Setter for the RICH reconstruction model (clusters)
//

    ((AliRICHChamber*) (*fChambers)[id])->ReconstructionModel(reconst);
}

void   AliRICH::SetNsec(Int_t id, Int_t nsec)
{

//
// Sets the number of padplanes
//

    ((AliRICHChamber*) (*fChambers)[id])->SetNsec(nsec);
}


//___________________________________________

void AliRICH::StepManager()
{

// Dummy step manager (should never be called)

}

void AliRICH::FindClusters(Int_t nev,Int_t lastEntry)
{

//
// Loop on chambers and on cathode planes
//
    for (Int_t icat=1;icat<2;icat++) {
	gAlice->ResetDigits();
	gAlice->TreeD()->GetEvent(1); // spurious +1 ...
	for (Int_t ich=0;ich<kNCH;ich++) {
	  AliRICHChamber* iChamber=(AliRICHChamber*) (*fChambers)[ich];
	  TClonesArray *pRICHdigits  = this->DigitsAddress(ich);
	  if (pRICHdigits == 0)	      
	      continue;
	  //
	  // Get ready the current chamber stuff
	  //
	  AliRICHResponse* response = iChamber->GetResponseModel();
	  AliRICHSegmentation*  seg = iChamber->GetSegmentationModel();
	  AliRICHClusterFinder* rec = iChamber->GetReconstructionModel();
	  if (seg) {	  
	      rec->SetSegmentation(seg);
	      rec->SetResponse(response);
	      rec->SetDigits(pRICHdigits);
	      rec->SetChamber(ich);
	      if (nev==0) rec->CalibrateCOG(); 
	      rec->FindRawClusters();
	  }  
	  TClonesArray *fRch;
	  fRch=RawClustAddress(ich);
	  fRch->Sort();
	} // for ich

	gAlice->TreeR()->Fill();
	TClonesArray *fRch;
	for (int i=0;i<kNCH;i++) {
	    fRch=RawClustAddress(i);
	    int nraw=fRch->GetEntriesFast();
	    printf ("Chamber %d, raw clusters %d\n",i,nraw);
	}
	
	ResetRawClusters();
	
    } // for icat
    
    char hname[30];
    sprintf(hname,"TreeR%d",nev);
    gAlice->TreeR()->Write(hname,kOverwrite,0);
    gAlice->TreeR()->Reset();
    
    //gObjectTable->Print();
}


//______________________________________________________________________________
void AliRICH::Streamer(TBuffer &R__b)
{
    // Stream an object of class AliRICH.
    AliRICHChamber       *iChamber;
    AliRICHSegmentation  *segmentation;
    AliRICHResponse      *response;
    TClonesArray         *digitsaddress;
    TClonesArray         *rawcladdress;
    TClonesArray         *rechitaddress;
      
    if (R__b.IsReading()) {
	Version_t R__v = R__b.ReadVersion(); if (R__v) { }
	AliDetector::Streamer(R__b);
	R__b >> fNPadHits;
	R__b >> fPadHits;   // diff
	R__b >> fNcerenkovs;
	R__b >> fCerenkovs; // diff
	R__b >> fDchambers;
	R__b >> fRawClusters;
	R__b >> fRecHits;  //diff
	R__b >> fDebugLevel;  //diff
	R__b.ReadArray(fNdch);
	R__b.ReadArray(fNrawch);
	R__b.ReadArray(fNrechits);
//
	R__b >> fChambers;
// Stream chamber related information
	for (Int_t i =0; i<kNCH; i++) {
	    iChamber=(AliRICHChamber*) (*fChambers)[i];
	    iChamber->Streamer(R__b);
	    segmentation=iChamber->GetSegmentationModel();
	    segmentation->Streamer(R__b);
	    response=iChamber->GetResponseModel();
	    response->Streamer(R__b);	  
	    rawcladdress=(TClonesArray*) (*fRawClusters)[i];
	    rawcladdress->Streamer(R__b);
	    rechitaddress=(TClonesArray*) (*fRecHits)[i];
	    rechitaddress->Streamer(R__b);
	    digitsaddress=(TClonesArray*) (*fDchambers)[i];
	    digitsaddress->Streamer(R__b);
	}
      
    } else {
	R__b.WriteVersion(AliRICH::IsA());
	AliDetector::Streamer(R__b);
	R__b << fNPadHits;
	R__b << fPadHits; // diff
	R__b << fNcerenkovs;
	R__b << fCerenkovs; // diff
	R__b << fDchambers;
	R__b << fRawClusters;
	R__b << fRecHits; //diff
	R__b << fDebugLevel; //diff
	R__b.WriteArray(fNdch, kNCH);
	R__b.WriteArray(fNrawch, kNCH);
	R__b.WriteArray(fNrechits, kNCH);
	//
	R__b << fChambers;
//  Stream chamber related information
	for (Int_t i =0; i<kNCH; i++) {
	    iChamber=(AliRICHChamber*) (*fChambers)[i];
	    iChamber->Streamer(R__b);
	    segmentation=iChamber->GetSegmentationModel();
	    segmentation->Streamer(R__b);
	    response=iChamber->GetResponseModel();
	    response->Streamer(R__b);
	    rawcladdress=(TClonesArray*) (*fRawClusters)[i];
	    rawcladdress->Streamer(R__b);
	    rechitaddress=(TClonesArray*) (*fRecHits)[i];
	    rechitaddress->Streamer(R__b);
	    digitsaddress=(TClonesArray*) (*fDchambers)[i];
	    digitsaddress->Streamer(R__b);
	}
    }
}
AliRICHPadHit* AliRICH::FirstPad(AliRICHHit*  hit,TClonesArray *clusters ) 
{
//
    // Initialise the pad iterator
    // Return the address of the first padhit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->fPHlast > 0) {
	sMaxIterPad=Int_t(hit->fPHlast);
	sCurIterPad=Int_t(hit->fPHfirst);
	return (AliRICHPadHit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
    
}

AliRICHPadHit* AliRICH::NextPad(TClonesArray *clusters) 
{

  // Iterates over pads
  
    sCurIterPad++;
    if (sCurIterPad <= sMaxIterPad) {
	return (AliRICHPadHit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}


void AliRICH::Digitise(Int_t nev, Int_t flag, Option_t *option,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !

    static Bool_t first=kTRUE;
    static TFile *pFile;
    char *addBackground = strstr(option,"Add");

    FILE* points; //these will be the digits...

    points=fopen("points.dat","w");

    AliRICHChamber*       iChamber;
    AliRICHSegmentation*  segmentation;

    Int_t digitse=0;
    Int_t trk[50];
    Int_t chtrk[50];  
    TObjArray *list=new TObjArray;
    static TClonesArray *pAddress=0;
    if(!pAddress) pAddress=new TClonesArray("TVector",1000);
    Int_t digits[5]; 
    
    AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
    AliRICHHitMap* pHitMap[10];
    Int_t i;
    for (i=0; i<10; i++) {pHitMap[i]=0;}
    if (addBackground ) {
	if(first) {
	    fFileName=filename;
	    cout<<"filename"<<fFileName<<endl;
	    pFile=new TFile(fFileName);
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
	    fHits2     = new TClonesArray("AliRICHHit",1000  );
	    fClusters2 = new TClonesArray("AliRICHPadHit",10000);
	    first=kFALSE;
	}
	pFile->cd();
	// Get Hits Tree header from file
	if(fHits2) fHits2->Clear();
	if(fClusters2) fClusters2->Clear();
	if(TrH1) delete TrH1;
	TrH1=0;

	char treeName[20];
	sprintf(treeName,"TreeH%d",nev);
	TrH1 = (TTree*)gDirectory->Get(treeName);
	if (!TrH1) {
	    printf("ERROR: cannot find Hits Tree for event:%d\n",nev);
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
	    branch = TrH1->GetBranch("RICHCluster");
	    if (branch) branch->SetAddress(&fClusters2);
	}
    }
    
    AliRICHHitMap* hm;
    Int_t countadr=0;
    Int_t counter=0;
    for (i =0; i<kNCH; i++) {
      iChamber=(AliRICHChamber*) (*fChambers)[i];
      segmentation=iChamber->GetSegmentationModel(1);
      pHitMap[i] = new AliRICHHitMapA1(segmentation, list);
    }
    //
    //   Loop over tracks
    //
    
    TTree *treeH = gAlice->TreeH();
    Int_t ntracks =(Int_t) treeH->GetEntries();
    for (Int_t track=0; track<ntracks; track++) {
      gAlice->ResetHits();
      treeH->GetEvent(track);
      //
      //   Loop over hits
      for(AliRICHHit* mHit=(AliRICHHit*)pRICH->FirstHit(-1); 
	  mHit;
	  mHit=(AliRICHHit*)pRICH->NextHit()) 
	{
	  
	  digitse=0;
	  
	  Int_t   nch   = mHit->fChamber-1;  // chamber number
	  if (nch >kNCH) continue;
	  iChamber = &(pRICH->Chamber(nch));
	  
	  TParticle *current = (TParticle*)(*gAlice->Particles())[track];
	  
	  Int_t particle = current->GetPdgCode();
	  
	  //printf("Flag:%d\n",flag);
	  //printf("Track:%d\n",track);
	  //printf("Particle:%d\n",particle);
	  
	  if (flag == 0)
	    digitse=1;
	  
	  if (flag == 1) 
	    if(TMath::Abs(particle) == 211 || TMath::Abs(particle) == 111)
	      digitse=1;
	  
	  if (flag == 2)
	    if(TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
	       || TMath::Abs(particle)==311)
	      digitse=1;
	  
	  if (flag == 3 && TMath::Abs(particle)==2212)
	    digitse=1;
	  
	  if (flag == 4 && TMath::Abs(particle)==13)
	    digitse=1;
	  
	  if (flag == 5 && TMath::Abs(particle)==11)
	    digitse=1;
	  
	  if (flag == 6 && TMath::Abs(particle)==2112)
	    digitse=1;
	  
	  
	  //printf ("Particle: %d, Flag: %d, Digitse: %d\n",particle,flag,digitse); 
	  
	  
	  if (digitse)
	    {
	      
	      //
	      // Loop over pad hits
	      for (AliRICHPadHit* mPad=
		     (AliRICHPadHit*)pRICH->FirstPad(mHit,fPadHits);
		   mPad;
		   mPad=(AliRICHPadHit*)pRICH->NextPad(fPadHits))
		{
		  Int_t cathode  = mPad->fCathode;    // cathode number
		  Int_t ipx      = mPad->fPadX;       // pad number on X
		  Int_t ipy      = mPad->fPadY;       // pad number on Y
		  Int_t iqpad    = mPad->fQpad;       // charge per pad
		  //
		  //
		  //printf("X:%d, Y:%d, Q:%d\n",ipx,ipy,iqpad);
		  
		  Float_t thex, they;
		  segmentation=iChamber->GetSegmentationModel(cathode);
		  segmentation->GetPadCxy(ipx,ipy,thex,they);
		  new((*pAddress)[countadr++]) TVector(2);
		  TVector &trinfo=*((TVector*) (*pAddress)[countadr-1]);
		  trinfo(0)=(Float_t)track;
		  trinfo(1)=(Float_t)iqpad;
		  
		  digits[0]=ipx;
		  digits[1]=ipy;
		  digits[2]=iqpad;
		  
		  AliRICHTransientDigit* pdigit;
		  // build the list of fired pads and update the info
		  if (!pHitMap[nch]->TestHit(ipx, ipy)) {
		    list->AddAtAndExpand(new AliRICHTransientDigit(nch,digits),counter);
		    pHitMap[nch]->SetHit(ipx, ipy, counter);
		    counter++;
		    pdigit=(AliRICHTransientDigit*)list->At(list->GetLast());
		    // list of tracks
		    TObjArray *trlist=(TObjArray*)pdigit->TrackList();
		    trlist->Add(&trinfo);
		  } else {
		    pdigit=(AliRICHTransientDigit*) pHitMap[nch]->GetHit(ipx, ipy);
		    // update charge
		    (*pdigit).fSignal+=iqpad;
		    // update list of tracks
		    TObjArray* trlist=(TObjArray*)pdigit->TrackList();
		    Int_t lastEntry=trlist->GetLast();
		    TVector *ptrkP=(TVector*)trlist->At(lastEntry);
		    TVector &ptrk=*ptrkP;
		    Int_t lastTrack=Int_t(ptrk(0));
		    Int_t lastCharge=Int_t(ptrk(1));
		    if (lastTrack==track) {
		      lastCharge+=iqpad;
		      trlist->RemoveAt(lastEntry);
		      trinfo(0)=lastTrack;
		      trinfo(1)=lastCharge;
		      trlist->AddAt(&trinfo,lastEntry);
		    } else {
		      trlist->Add(&trinfo);
		    }
		    // check the track list
		    Int_t nptracks=trlist->GetEntriesFast();
		    if (nptracks > 2) {
		      printf("Attention - tracks:  %d (>2)\n",nptracks);
		      //printf("cat,nch,ix,iy %d %d %d %d  \n",icat+1,nch,ipx,ipy);
		      for (Int_t tr=0;tr<nptracks;tr++) {
			TVector *pptrkP=(TVector*)trlist->At(tr);
			TVector &pptrk=*pptrkP;
			trk[tr]=Int_t(pptrk(0));
			chtrk[tr]=Int_t(pptrk(1));
		      }
		    } // end if nptracks
		  } //  end if pdigit
		} //end loop over clusters
	    }// track type condition
	} // hit loop
    } // track loop
    
    // open the file with background
    
    if (addBackground ) {
      ntracks =(Int_t)TrH1->GetEntries();
      //printf("background - icat,ntracks1  %d %d\n",icat,ntracks);
      //printf("background - Start loop over tracks \n");     
      //
      //   Loop over tracks
      //
      for (Int_t trak=0; trak<ntracks; trak++) {
	if (fHits2)       fHits2->Clear();
	if (fClusters2)   fClusters2->Clear();
	TrH1->GetEvent(trak);
	//
	//   Loop over hits
	AliRICHHit* mHit;
	for(int j=0;j<fHits2->GetEntriesFast();++j) 
	  {
	    mHit=(AliRICHHit*) (*fHits2)[j];
	    Int_t   nch   = mHit->fChamber-1;  // chamber number
	    if (nch >6) continue;
	    iChamber = &(pRICH->Chamber(nch));
	    Int_t rmin = (Int_t)iChamber->RInner();
	    Int_t rmax = (Int_t)iChamber->ROuter();
	    //
	    // Loop over pad hits
	    for (AliRICHPadHit* mPad=
		   (AliRICHPadHit*)pRICH->FirstPad(mHit,fClusters2);
		 mPad;
		 mPad=(AliRICHPadHit*)pRICH->NextPad(fClusters2))
	      {
		Int_t cathode  = mPad->fCathode;    // cathode number
		Int_t ipx      = mPad->fPadX;       // pad number on X
		Int_t ipy      = mPad->fPadY;       // pad number on Y
		Int_t iqpad    = mPad->fQpad;       // charge per pad
		
		Float_t thex, they;
		segmentation=iChamber->GetSegmentationModel(cathode);
		segmentation->GetPadCxy(ipx,ipy,thex,they);
		Float_t rpad=TMath::Sqrt(thex*thex+they*they);
		if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;
		new((*pAddress)[countadr++]) TVector(2);
		TVector &trinfo=*((TVector*) (*pAddress)[countadr-1]);
		trinfo(0)=-1;  // tag background
		trinfo(1)=-1;
		digits[0]=ipx;
		digits[1]=ipy;
		digits[2]=iqpad;
		if (trak <4 && nch==0)
		  printf("bgr - pHitMap[nch]->TestHit(ipx, ipy),trak %d %d\n",
			 pHitMap[nch]->TestHit(ipx, ipy),trak);
		AliRICHTransientDigit* pdigit;
		// build the list of fired pads and update the info
		if (!pHitMap[nch]->TestHit(ipx, ipy)) {
		  list->AddAtAndExpand(new AliRICHTransientDigit(nch,digits),counter);
		  
		  pHitMap[nch]->SetHit(ipx, ipy, counter);
		  counter++;
		  printf("bgr new elem in list - counter %d\n",counter);
		  
		  pdigit=(AliRICHTransientDigit*)list->At(list->GetLast());
		  // list of tracks
		  TObjArray *trlist=(TObjArray*)pdigit->TrackList();
		  trlist->Add(&trinfo);
		} else {
		  pdigit=(AliRICHTransientDigit*) pHitMap[nch]->GetHit(ipx, ipy);
		  // update charge
		  (*pdigit).fSignal+=iqpad;
		  // update list of tracks
		  TObjArray* trlist=(TObjArray*)pdigit->TrackList();
		  Int_t lastEntry=trlist->GetLast();
		  TVector *ptrkP=(TVector*)trlist->At(lastEntry);
		  TVector &ptrk=*ptrkP;
		  Int_t lastTrack=Int_t(ptrk(0));
		  if (lastTrack==-1) {
		    continue;
		  } else {
		    trlist->Add(&trinfo);
		  }
		  // check the track list
		  Int_t nptracks=trlist->GetEntriesFast();
		  if (nptracks > 0) {
		    for (Int_t tr=0;tr<nptracks;tr++) {
		      TVector *pptrkP=(TVector*)trlist->At(tr);
		      TVector &pptrk=*pptrkP;
		      trk[tr]=Int_t(pptrk(0));
		      chtrk[tr]=Int_t(pptrk(1));
		    }
		  } // end if nptracks
		} //  end if pdigit
	      } //end loop over clusters
	  } // hit loop
      } // track loop
	    TTree *fAli=gAlice->TreeK();
	    if (fAli) pFile =fAli->GetCurrentFile();
	    pFile->cd();
    } // if Add	
    
    Int_t tracks[10];
    Int_t charges[10];
    //cout<<"Start filling digits \n "<<endl;
    Int_t nentries=list->GetEntriesFast();
    //printf(" \n \n nentries %d \n",nentries);
    
    // start filling the digits
    
    for (Int_t nent=0;nent<nentries;nent++) {
      AliRICHTransientDigit *address=(AliRICHTransientDigit*)list->At(nent);
      if (address==0) continue; 
      
      Int_t ich=address->fChamber;
      Int_t q=address->fSignal; 
      iChamber=(AliRICHChamber*) (*fChambers)[ich];
      AliRICHResponse * response=iChamber->GetResponseModel();
      Int_t adcmax= (Int_t) response->MaxAdc();
      
      
      // add white noise and do zero-suppression and signal truncation (new electronics,old electronics gaus 1.2,0.2)
      //printf("Treshold: %d\n",iChamber->fTresh->GetHitIndex(address->fPadX,address->fPadY));
      Int_t pedestal = iChamber->fTresh->GetHitIndex(address->fPadX,address->fPadY);

      //printf("Pedestal:%d\n",pedestal);
      //Int_t pedestal=0;
      Float_t treshold = (pedestal + 4*1.7);
      
      Float_t meanNoise = gRandom->Gaus(1.7, 0.25);
      Float_t noise     = gRandom->Gaus(0, meanNoise);
      q+=(Int_t)(noise + pedestal);
      //q+=(Int_t)(noise);
      //          magic number to be parametrised !!! 
      if ( q <= treshold) 
	{
	  q = q - pedestal;
	  continue;
	}
      q = q - pedestal;
      if ( q >= adcmax) q=adcmax;
      digits[0]=address->fPadX;
      digits[1]=address->fPadY;
      digits[2]=q;
      
      TObjArray* trlist=(TObjArray*)address->TrackList();
      Int_t nptracks=trlist->GetEntriesFast();
      
      // this was changed to accomodate the real number of tracks
      if (nptracks > 10) {
	cout<<"Attention - tracks > 10 "<<nptracks<<endl;
	nptracks=10;
      }
      if (nptracks > 2) {
	printf("Attention - tracks > 2  %d \n",nptracks);
	//printf("cat,ich,ix,iy,q %d %d %d %d %d \n",
	//icat,ich,digits[0],digits[1],q);
      }
      for (Int_t tr=0;tr<nptracks;tr++) {
	TVector *ppP=(TVector*)trlist->At(tr);
	TVector &pp  =*ppP;
	tracks[tr]=Int_t(pp(0));
	charges[tr]=Int_t(pp(1));
      }      //end loop over list of tracks for one pad
      if (nptracks < 10 ) {
	for (Int_t t=nptracks; t<10; t++) {
	  tracks[t]=0;
	  charges[t]=0;
	}
      }
      //write file
      if (ich==2)
	fprintf(points,"%4d,      %4d,      %4d\n",digits[0],digits[1],digits[2]);
      
      // fill digits
      pRICH->AddDigits(ich,tracks,charges,digits);
    }	
    gAlice->TreeD()->Fill();
    
    list->Delete();
    for(Int_t ii=0;ii<kNCH;++ii) {
      if (pHitMap[ii]) {
	hm=pHitMap[ii];
	delete hm;
	pHitMap[ii]=0;
      }
    }
    
    //TTree *TD=gAlice->TreeD();
    //Stat_t ndig=TD->GetEntries();
    //cout<<"number of digits  "<<ndig<<endl;
    TClonesArray *fDch;
    for (int k=0;k<kNCH;k++) {
      fDch= pRICH->DigitsAddress(k);
      int ndigit=fDch->GetEntriesFast();
      printf ("Chamber %d digits %d \n",k,ndigit);
    }
    pRICH->ResetDigits();
    char hname[30];
    sprintf(hname,"TreeD%d",nev);
    gAlice->TreeD()->Write(hname,kOverwrite,0);
    
    // reset tree
    //    gAlice->TreeD()->Reset();
    delete list;
    pAddress->Clear();
    // gObjectTable->Print();
}

AliRICH& AliRICH::operator=(const AliRICH& rhs)
{
// Assignment operator
    return *this;
    
}

Int_t AliRICH::MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res)
{
//
//  Calls the charge disintegration method of the current chamber and adds
//  the simulated cluster to the root treee 
//
    Int_t clhits[kNCH];
    Float_t newclust[6][500];
    Int_t nnew;
    
//
//  Integrated pulse height on chamber
    
    clhits[0]=fNhits+1;
    
    ((AliRICHChamber*) (*fChambers)[idvol])->DisIntegration(eloss, xhit, yhit, nnew, newclust, res);
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
	    
	    AddPadHit(clhits);
	}
    }
return nnew;
}

