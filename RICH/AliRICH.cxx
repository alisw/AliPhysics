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

#include "AliRICH.h"
#include "AliRICHHitMap.h"
#include "AliRICHClusterFinder.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliPoints.h"
#include "iostream.h"
#include "AliCallf77.h" 
#include "TParticle.h"

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
    
    fNdch      = new Int_t[7];
    
    fDchambers = new TObjArray(7);

    fRecHits = new TObjArray(7);
    
    Int_t i;
   
    for (i=0; i<7 ;i++) {
	(*fDchambers)[i] = new TClonesArray("AliRICHDigit",10000); 
	fNdch[i]=0;
    }

    fNrawch      = new Int_t[7];
    
    fRawClusters = new TObjArray(7);
    //printf("Created fRwClusters with adress:%p",fRawClusters);

    for (i=0; i<7 ;i++) {
      (*fRawClusters)[i] = new TClonesArray("AliRICHRawCluster",10000); 
      fNrawch[i]=0;
    }

    fNrechits      = new Int_t[7];
    
    for (i=0; i<7 ;i++) {
	(*fRecHits)[i] = new TClonesArray("AliRICHRecHit",1000); 
    }
    //printf("Created fRecHits with adress:%p",fRecHits);

        
    SetMarkerColor(kRed);
}

//___________________________________________
AliRICH::~AliRICH()
{
    fIshunt  = 0;
    delete fHits;
    delete fPadHits;
    delete fCerenkovs;
}

//___________________________________________
void AliRICH::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
    TClonesArray &lhits = *fHits;
    new(lhits[fNhits++]) AliRICHHit(fIshunt,track,vol,hits);
}
//_____________________________________________________________________________
void AliRICH::AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs)
{
    TClonesArray &lcerenkovs = *fCerenkovs;
    new(lcerenkovs[fNcerenkovs++]) AliRICHCerenkov(fIshunt,track,vol,cerenkovs);
    //printf ("Done for Cerenkov %d\n\n\n\n",fNcerenkovs);
}
//___________________________________________
void AliRICH::AddPadHit(Int_t *clhits)
{
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
void AliRICH::AddRecHit(Int_t id, Float_t *rechit)
{
    //
    // Add a RICH reconstructed hit to the list
    //

    TClonesArray &lrec = *((TClonesArray*)(*fRecHits)[id]);
    new(lrec[fNrechits[id]++]) AliRICHRecHit(id,rechit);
}

//___________________________________________
void AliRICH::BuildGeometry()
    
{
    //
    // Builds a TNode geometry for event display
    //
    TNode *Node, *Top;
    
    const int kColorRICH = kGreen;
    //
    Top=gAlice->GetGeometry()->GetNode("alice");
    
    
    new TBRIK("S_RICH","S_RICH","void",71.09999,11.5,73.15);
    
    Top->cd();
    Float_t pos1[3]={0,471.8999,165.2599};
    //Chamber(0).SetChamberTransform(pos1[0],pos1[1],pos1[2],
    new TRotMatrix("rot993","rot993",90,0,70.69,90,19.30999,-90);
    Node = new TNode("RICH1","RICH1","S_RICH",pos1[0],pos1[1],pos1[2],"rot993");
    

    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    
    Float_t pos2[3]={171,470,0};
    //Chamber(1).SetChamberTransform(pos2[0],pos2[1],pos2[2],
    new TRotMatrix("rot994","rot994",90,-20,90,70,0,0);
    Node = new TNode("RICH2","RICH2","S_RICH",pos2[0],pos2[1],pos2[2],"rot994");
    
    
    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos3[3]={0,500,0};
    //Chamber(2).SetChamberTransform(pos3[0],pos3[1],pos3[2],
    new TRotMatrix("rot995","rot995",90,0,90,90,0,0);
    Node = new TNode("RICH3","RICH3","S_RICH",pos3[0],pos3[1],pos3[2],"rot995");
    

    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos4[3]={-171,470,0};
    //Chamber(3).SetChamberTransform(pos4[0],pos4[1],pos4[2], 
    new TRotMatrix("rot996","rot996",90,20,90,110,0,0);  
    Node = new TNode("RICH4","RICH4","S_RICH",pos4[0],pos4[1],pos4[2],"rot996");
    

    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos5[3]={161.3999,443.3999,-165.3};
    //Chamber(4).SetChamberTransform(pos5[0],pos5[1],pos5[2],
    new TRotMatrix("rot997","rot997",90,340,108.1999,70,18.2,70);
    Node = new TNode("RICH5","RICH5","S_RICH",pos5[0],pos5[1],pos5[2],"rot997");
    
    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos6[3]={0., 471.9, -165.3,};
    //Chamber(5).SetChamberTransform(pos6[0],pos6[1],pos6[2],
    new TRotMatrix("rot998","rot998",90,0,109.3099,90,19.30999,90);
    Node = new TNode("RICH6","RICH6","S_RICH",pos6[0],pos6[1],pos6[2],"rot998");
    
    
    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos7[3]={-161.399,443.3999,-165.3};
    //Chamber(6).SetChamberTransform(pos7[0],pos7[1],pos7[2],
    new TRotMatrix("rot999","rot999",90,20,108.1999,110,18.2,110);
    Node = new TNode("RICH7","RICH7","S_RICH",pos7[0],pos7[1],pos7[2],"rot999");
    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node); 
    
}

//___________________________________________
Int_t AliRICH::DistancetoPrimitive(Int_t , Int_t )
{
    return 9999;
}

//___________________________________________
void AliRICH::MakeBranch(Option_t* option)
{
    // Create Tree branches for the RICH.
    
    const Int_t buffersize = 4000;
    char branchname[20];
    
    
    AliDetector::MakeBranch(option);
    sprintf(branchname,"%sCerenkov",GetName());
    if (fCerenkovs   && gAlice->TreeH()) {
	gAlice->TreeH()->Branch(branchname,&fCerenkovs, buffersize);
	printf("Making Branch %s for Cerenkov Hits\n",branchname);
    }
    
    sprintf(branchname,"%sPadHits",GetName());
    if (fPadHits   && gAlice->TreeH()) {
	gAlice->TreeH()->Branch(branchname,&fPadHits, buffersize);
	printf("Making Branch %s for PadHits\n",branchname);
    }
    
// one branch for digits per chamber
    Int_t i;
    
    for (i=0; i<7 ;i++) {
	sprintf(branchname,"%sDigits%d",GetName(),i+1);
	
	if (fDchambers   && gAlice->TreeD()) {
	    gAlice->TreeD()->Branch(branchname,&((*fDchambers)[i]), buffersize);
	    printf("Making Branch %s for digits in chamber %d\n",branchname,i+1);
	}	
    }

// one branch for raw clusters per chamber
  for (i=0; i<7 ;i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
      
      if (fRawClusters   && gAlice->TreeR()) {
	 gAlice->TreeR()->Branch(branchname,&((*fRawClusters)[i]), buffersize);
	 printf("Making Branch %s for raw clusters in chamber %d\n",branchname,i+1);
      }	
  }

  // one branch for rec hits per chamber
  for (i=0; i<7 ;i++) {
    sprintf(branchname,"%sRecHits%d",GetName(),i+1);
    
    if (fRecHits   && gAlice->TreeR()) {
      gAlice->TreeR()->Branch(branchname,&((*fRecHits)[i]), buffersize);
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
	for (int i=0; i<7; i++) {
	    sprintf(branchname,"%sDigits%d",GetName(),i+1);
	    if (fDchambers) {
		branch = treeD->GetBranch(branchname);
		if (branch) branch->SetAddress(&((*fDchambers)[i]));
	    }
	}
    }
  if (treeR) {
      for (i=0; i<7; i++) {
	  sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	  if (fRawClusters) {
	      branch = treeR->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	  }
      }
      
      for (i=0; i<7; i++) {
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
    for ( int i=0;i<7;i++ ) {
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
    for ( int i=0;i<7;i++ ) {
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
  
  for ( int i=0;i<7;i++ ) {
	if ((*fRecHits)[i])    ((TClonesArray*)(*fRecHits)[i])->Clear();
	if (fNrechits)  fNrechits[i]=0;
    }
}

//___________________________________________
void   AliRICH::SetGeometryModel(Int_t id, AliRICHGeometry *geometry)
{
    ((AliRICHChamber*) (*fChambers)[id])->GeometryModel(geometry);
}

//___________________________________________
void   AliRICH::SetSegmentationModel(Int_t id, AliRICHSegmentation *segmentation)
{
    ((AliRICHChamber*) (*fChambers)[id])->SegmentationModel(segmentation);
}

//___________________________________________
void   AliRICH::SetResponseModel(Int_t id, AliRICHResponse *response)
{
    ((AliRICHChamber*) (*fChambers)[id])->ResponseModel(response);
}

void   AliRICH::SetReconstructionModel(Int_t id, AliRICHClusterFinder *reconst)
{
    ((AliRICHChamber*) (*fChambers)[id])->ReconstructionModel(reconst);
}

void   AliRICH::SetNsec(Int_t id, Int_t nsec)
{
    ((AliRICHChamber*) (*fChambers)[id])->SetNsec(nsec);
}


//___________________________________________

void AliRICH::StepManager()
{
}

void AliRICH::FindClusters(Int_t nev,Int_t last_entry)
{

//
// Loop on chambers and on cathode planes
//
    for (Int_t icat=1;icat<2;icat++) {
	gAlice->ResetDigits();
	gAlice->TreeD()->GetEvent(1); // spurious +1 ...
	for (Int_t ich=0;ich<7;ich++) {
	  AliRICHChamber* iChamber=(AliRICHChamber*) (*fChambers)[ich];
	  TClonesArray *RICHdigits  = this->DigitsAddress(ich);
	  if (RICHdigits == 0)	      
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
	      rec->SetDigits(RICHdigits);
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
	for (int i=0;i<7;i++) {
	    fRch=RawClustAddress(i);
	    int nraw=fRch->GetEntriesFast();
	    printf ("Chamber %d, raw clusters %d\n",i,nraw);
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
	R__b.ReadArray(fNdch);
	R__b.ReadArray(fNrawch);
	R__b.ReadArray(fNrechits);
//
	R__b >> fChambers;
// Stream chamber related information
	for (Int_t i =0; i<7; i++) {
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
	R__b.WriteArray(fNdch, 7);
	R__b.WriteArray(fNrawch, 7);
	R__b.WriteArray(fNrechits, 7);
	//
	R__b << fChambers;
//  Stream chamber related information
	for (Int_t i =0; i<7; i++) {
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
    static TFile *File;
    char *Add = strstr(option,"Add");

    FILE* points; //these will be the digits...

    points=fopen("points.dat","w");

    AliRICHChamber*       iChamber;
    AliRICHSegmentation*  segmentation;

    Int_t digitse=0;
    Int_t trk[50];
    Int_t chtrk[50];  
    TObjArray *list=new TObjArray;
    static TClonesArray *p_adr=0;
    if(!p_adr) p_adr=new TClonesArray("TVector",1000);
    Int_t digits[5]; 
    
    AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH");
    AliRICHHitMap* HitMap[10];
    Int_t i;
    for (i=0; i<10; i++) {HitMap[i]=0;}
    if (Add ) {
	if(first) {
	    fFileName=filename;
	    cout<<"filename"<<fFileName<<endl;
	    File=new TFile(fFileName);
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
	    fHits2     = new TClonesArray("AliRICHHit",1000  );
	    fClusters2 = new TClonesArray("AliRICHPadHit",10000);
	    first=kFALSE;
	}
	File->cd();
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
    //
    // loop over cathodes
    //
    AliRICHHitMap* hm;
    Int_t countadr=0;
    for (int icat=0; icat<1; icat++) { 
	Int_t counter=0;
	for (i =0; i<7; i++) {
	    iChamber=(AliRICHChamber*) (*fChambers)[i];
	    if (iChamber->Nsec()==1 && icat==1) {
		continue;
	    } else {
		segmentation=iChamber->GetSegmentationModel(icat+1);
	    }
	    HitMap[i] = new AliRICHHitMapA1(segmentation, list);
	}
//
//   Loop over tracks
//

	TTree *TH = gAlice->TreeH();
	Int_t ntracks =(Int_t) TH->GetEntries();
	for (Int_t track=0; track<ntracks; track++) {
	    gAlice->ResetHits();
	    TH->GetEvent(track);
//
//   Loop over hits
	    for(AliRICHHit* mHit=(AliRICHHit*)RICH->FirstHit(-1); 
		mHit;
		mHit=(AliRICHHit*)RICH->NextHit()) 
	    {
	      
	      digitse=0;
	      
	      Int_t   nch   = mHit->fChamber-1;  // chamber number
	      if (nch >7) continue;
	      iChamber = &(RICH->Chamber(nch));
	      
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
			 (AliRICHPadHit*)RICH->FirstPad(mHit,fPadHits);
		       mPad;
		       mPad=(AliRICHPadHit*)RICH->NextPad(fPadHits))
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
		      new((*p_adr)[countadr++]) TVector(2);
		      TVector &trinfo=*((TVector*) (*p_adr)[countadr-1]);
		      trinfo(0)=(Float_t)track;
		      trinfo(1)=(Float_t)iqpad;
		      
		      digits[0]=ipx;
		      digits[1]=ipy;
		      digits[2]=iqpad;
		      
		      AliRICHTransientDigit* pdigit;
		      // build the list of fired pads and update the info
		      if (!HitMap[nch]->TestHit(ipx, ipy)) {
			list->AddAtAndExpand(new AliRICHTransientDigit(nch,digits),counter);
			HitMap[nch]->SetHit(ipx, ipy, counter);
			counter++;
			pdigit=(AliRICHTransientDigit*)list->At(list->GetLast());
			// list of tracks
			TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			trlist->Add(&trinfo);
		      } else {
			pdigit=(AliRICHTransientDigit*) HitMap[nch]->GetHit(ipx, ipy);
			// update charge
			(*pdigit).fSignal+=iqpad;
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
			    printf("Attention - tracks:  %d (>2)\n",nptracks);
			    //printf("cat,nch,ix,iy %d %d %d %d  \n",icat+1,nch,ipx,ipy);
			    for (Int_t tr=0;tr<nptracks;tr++) {
				TVector *pptrk_p=(TVector*)trlist->At(tr);
				TVector &pptrk=*pptrk_p;
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
	
	if (Add ) {
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
		    iChamber = &(RICH->Chamber(nch));
		    Int_t rmin = (Int_t)iChamber->RInner();
		    Int_t rmax = (Int_t)iChamber->ROuter();
//
// Loop over pad hits
		    for (AliRICHPadHit* mPad=
			     (AliRICHPadHit*)RICH->FirstPad(mHit,fClusters2);
			 mPad;
			 mPad=(AliRICHPadHit*)RICH->NextPad(fClusters2))
		    {
			Int_t cathode  = mPad->fCathode;    // cathode number
			Int_t ipx      = mPad->fPadX;       // pad number on X
			Int_t ipy      = mPad->fPadY;       // pad number on Y
			Int_t iqpad    = mPad->fQpad;       // charge per pad
			if (trak==3 && nch==0 && icat==0) printf("bgr - trak,iqpad,ipx,ipy %d %d %d %d\n",trak,iqpad,ipx,ipy);
//
//
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
			if (trak <4 && icat==0 && nch==0)
			    printf("bgr - HitMap[nch]->TestHit(ipx, ipy),trak %d %d\n",
				   HitMap[nch]->TestHit(ipx, ipy),trak);
			AliRICHTransientDigit* pdigit;
			// build the list of fired pads and update the info
			if (!HitMap[nch]->TestHit(ipx, ipy)) {
			    list->AddAtAndExpand(new AliRICHTransientDigit(nch,digits),counter);
			    
			    HitMap[nch]->SetHit(ipx, ipy, counter);
			    counter++;
			    printf("bgr new elem in list - counter %d\n",counter);
			    
			    pdigit=(AliRICHTransientDigit*)list->At(list->GetLast());
			    // list of tracks
			    TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			    trlist->Add(&trinfo);
			} else {
			    pdigit=(AliRICHTransientDigit*) HitMap[nch]->GetHit(ipx, ipy);
			    // update charge
			    (*pdigit).fSignal+=iqpad;
			    // update list of tracks
			    TObjArray* trlist=(TObjArray*)pdigit->TrackList();
			    Int_t last_entry=trlist->GetLast();
			    TVector *ptrk_p=(TVector*)trlist->At(last_entry);
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
	    TTree *fAli=gAlice->TreeK();
	    if (fAli) File =fAli->GetCurrentFile();
	    File->cd();
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
	    Float_t MeanNoise = gRandom->Gaus(1.7, 0.25);
	    Float_t Noise     = gRandom->Gaus(0, MeanNoise);
	    q+=(Int_t)Noise;
//          magic number to be parametrised !!! 
	    if ( q <= 6.8) continue;
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
		TVector *pp_p=(TVector*)trlist->At(tr);
		TVector &pp  =*pp_p;
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
	    RICH->AddDigits(ich,tracks,charges,digits);
	}
	gAlice->TreeD()->Fill();

	list->Delete();
	for(Int_t ii=0;ii<7;++ii) {
	    if (HitMap[ii]) {
		hm=HitMap[ii];
		delete hm;
		HitMap[ii]=0;
	    }
	}
	
	//TTree *TD=gAlice->TreeD();
	//Stat_t ndig=TD->GetEntries();
	//cout<<"number of digits  "<<ndig<<endl;
	TClonesArray *fDch;
	for (int k=0;k<7;k++) {
	    fDch= RICH->DigitsAddress(k);
	    int ndigit=fDch->GetEntriesFast();
	    printf ("Chamber %d digits %d \n",k,ndigit);
	}
	RICH->ResetDigits();
    } //end loop over cathodes
    char hname[30];
    sprintf(hname,"TreeD%d",nev);
    gAlice->TreeD()->Write(hname);

// reset tree
//    gAlice->TreeD()->Reset();
    delete list;
    p_adr->Clear();
// gObjectTable->Print();
}



