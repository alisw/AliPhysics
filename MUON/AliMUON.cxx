////////////////////////////////////////////////
//  Manager and hits classes for set:MUON     //
////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 
#include <TObject.h>
#include <TVector.h>
#include <TObjArray.h>

#include "AliMUON.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h" 

// Static variables for the pad-hit iterator routines
static Int_t sMaxIterPad=0;
static Int_t sCurIterPad=0;
 
ClassImp(AliMUON)
 
//___________________________________________
AliMUON::AliMUON()
{
   fIshunt     = 0;
   fHits       = 0;
   fClusters   = 0;
   fNclusters  = 0;
   fDchambers  = 0;
   fRecClusters= 0;
   fNdch       = 0;
}
 
//___________________________________________
AliMUON::AliMUON(const char *name, const char *title)
       : AliDetector(name,title)
{
//Begin_Html
/*
<img src="picts/alimuon.gif">
*/
//End_Html
 
   fHits     = new TClonesArray("AliMUONhit",1000  );
   fClusters = new TClonesArray("AliMUONcluster",10000);
   fNclusters  =  0;
   fIshunt     =  0;

   fNdch      = new Int_t[11];

   fDchambers = new TObjArray(11);

   Int_t i;
   
   for (i=0; i<11 ;i++) {
       (*fDchambers)[i] = new TClonesArray("AliMUONdigit",10000); 
       fNdch[i]=0;
   }

   fRecClusters=new TObjArray(20);
   for (i=0; i<20;i++)
     (*fRecClusters)[i] = new TObjArray(1000);

//   
// Transport angular cut
   fAccCut=0;
   fAccMin=2;
   fAccMax=9;

   SetMarkerColor(kRed);
}
 
//___________________________________________
AliMUON::~AliMUON()
{
  fIshunt  = 0;
  delete fHits;
  delete fClusters;
//  for (Int_t i=0;i<10;i++) {
//      delete (*fDchambers)[i];
//      fNdch[i]=0;
//  }
//  delete fDchambers;
  for (Int_t i=0;i<20;i++) 
      delete (*fRecClusters)[i];
  delete fRecClusters;

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
void AliMUON::AddRecCluster(Int_t iCh, Int_t iCat, AliMUONRecCluster* Cluster)
{
    //
    // Add a MUON reconstructed cluster to the list
    //
    TObjArray* ClusterList = RecClusters(iCh,iCat);
    ClusterList->Add(Cluster);
}

//___________________________________________
void AliMUON::BuildGeometry()
{
  TNode *Node, *Top;
  const int kColorMUON = kBlue;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");

  // MUON
  const float cz[10] = { 511, 519, 686, 694, 971, 979, 1245, 1253, 1445, 1453};
  const float acc_min = TMath::Tan( 2*.0174532925199432958);
  const float acc_max = TMath::Tan(9*.0174532925199432958);
  float rmin, rmax;
  
  // Chamber 1
  rmin = (cz[0]+0.25)*acc_min;
  rmax = cz[0]*acc_max;
  new TTUBE("S_MUON1","MUON chamber 1","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON1","MUON chamber 1","S_MUON1",0,0,cz[0],"");
  Node->SetLineColor(kColorMUON);
  fNodes->Add(Node);
  
  // Chamber 2
  rmin = (cz[1]+0.25)*acc_min;
  rmax = cz[1]*acc_max;
  new TTUBE("S_MUON2","MUON chamber 2","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON2","MUON chamber 2","S_MUON2",0,0,cz[1],"");
  fNodes->Add(Node);
  Node->SetLineColor(kColorMUON);
  
  // Chamber 3
  rmin = (cz[2]+0.25)*acc_min;
  rmax =  cz[2]*acc_max;
  new TTUBE("S_MUON3","MUON chamber 3","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON3","MUON chamber 3","S_MUON3",0,0,cz[2],"");
  Node->SetLineColor(kColorMUON);
  fNodes->Add(Node);
  
  // Chamber 4
  rmin = (cz[3]+0.25)*acc_min;
  rmax = cz[3]*acc_max;
  new TTUBE("S_MUON4","MUON chamber 4","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON4","MUON chamber 4","S_MUON4",0,0,cz[3],"");
  Node->SetLineColor(kColorMUON);
  fNodes->Add(Node);

  // Chamber 5
  rmin = 30;
  rmax = cz[4]*acc_max;
  new TTUBE("S_MUON5","MUON chamber 5","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON5","MUON chamber 5","S_MUON5",0,0,cz[4],"");
  Node->SetLineColor(kColorMUON);
  fNodes->Add(Node);
  
  // Chamber 6
  rmin = 30;
  rmax = cz[5]*acc_max;
  new TTUBE("S_MUON6","MUON chamber 6","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON6","MUON chamber 6","S_MUON6",0,0,cz[5],"");
  fNodes->Add(Node);
  Node->SetLineColor(kColorMUON);
  
  // Chamber 7
  rmin = 30;
  rmax = cz[6]*acc_max;
  new TTUBE("S_MUON7","MUON chamber 7","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON7","MUON chamber 7","S_MUON7",0,0,cz[6],"");
  Node->SetLineColor(kColorMUON);
  fNodes->Add(Node);

  // Chamber 8
  rmin = 30;
  rmax = cz[7]*acc_max;
  new TTUBE("S_MUON8","MUON chamber 8","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON8","MUON chamber 8","S_MUON8",0,0,cz[7],"");
  Node->SetLineColor(kColorMUON);
  fNodes->Add(Node);

  // Chamber 9
  rmin = 30;
  rmax = cz[8]*acc_max;
  new TTUBE("S_MUON9","MUON chamber 9","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON9","MUON chamber 9","S_MUON9",0,0,cz[8],"");
  Node->SetLineColor(kColorMUON);
  fNodes->Add(Node);

  // Chamber 10
  rmin = 30;
  rmax = cz[9]*acc_max;
  new TTUBE("S_MUON10","MUON chamber 10","void",rmin,rmax,0.25);
  Top->cd();
  Node = new TNode("MUON10","MUON chamber 10","S_MUON10",0,0,cz[9],"");
  Node->SetLineColor(kColorMUON);
  fNodes->Add(Node);

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
  char branchname[20];
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
// one branch for rec clusters
  for (i=0; i<20; i++) {
      sprintf(branchname,"%sRecClus%d",GetName(),i+1);
      if (fRecClusters   && gAlice->TreeD()) {
	  gAlice->TreeR()
	      ->Branch(branchname,"TObjArray", 
		       &((*fRecClusters)[i]), buffersize,0);
	  printf("Making Branch %s for clusters in chamber %d\n",
		 branchname,i+1);
      }
  }
}

//___________________________________________
void AliMUON::SetTreeAddress()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[20];
  AliDetector::SetTreeAddress();

  TBranch *branch;
  TTree *treeH = gAlice->TreeH();
  TTree *treeD = gAlice->TreeD();

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
	if ((*fDchambers)[i])   (*fDchambers)[i]->Clear();
	if (fNdch)  fNdch[i]=0;
    }
}
//____________________________________________
void AliMUON::ResetRecClusters()
{
    //
    // Reset the rec clusters
    //
    for ( int i=0;i<20;i++ ) {
	if ((*fRecClusters)[i])   (*fRecClusters)[i]->Clear();
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
void AliMUON::SetMUCHSP(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])->SetMUCHSP(p1);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetMUCHSP(p1);
}

//___________________________________________
void AliMUON::SetMUSIGM(Int_t id, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])->SetMUSIGM(p1,p2);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetMUSIGM(p1,p2);
}

//___________________________________________
void AliMUON::SetRSIGM(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])->SetRSIGM(p1);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetRSIGM(p1);
}

//___________________________________________
void AliMUON::SetMAXADC(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONchamber*) (*fChambers)[i])->SetMAXADC(p1);
    ((AliMUONchamber*) (*fChambers)[i+1])->SetMAXADC(p1);
}

//___________________________________________
void AliMUON::SetSMAXAR(Float_t p1)
{
     fMaxStepGas=p1;
}

//___________________________________________
void AliMUON::SetSMAXAL(Float_t p1)
{
    fMaxStepAlu=p1;
}

//___________________________________________
void AliMUON::SetDMAXAR(Float_t p1)
{
    fMaxDestepGas=p1;
}

//___________________________________________
void AliMUON::SetDMAXAL(Float_t p1)
{
    fMaxDestepAlu=p1;
}
//___________________________________________
void AliMUON::SetMUONACC(Bool_t acc, Float_t angmin, Float_t angmax)
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

void   AliMUON::SetNsec(Int_t id, Int_t nsec)
{
    ((AliMUONchamber*) (*fChambers)[id])->SetNsec(nsec);
}


//___________________________________________

void AliMUON::StepManager()
{
    printf("Dummy version of muon step -- it should never happen!!\n");
    const Float_t kRaddeg = 180/TMath::Pi();
    AliMC* pMC = AliMC::GetMC();
    Int_t nsec, ipart;
    Float_t x[4], p[4];
    Float_t pt, th0, th1;
    char proc[5];
    if(fAccCut) {
	if((nsec=pMC->NSecondaries())>0) {
	    pMC->ProdProcess(proc);
	    if((pMC->TrackPid()==113 || pMC->TrackPid()==114) && !strcmp(proc,"DCAY")) {
		//
		// Check angular acceptance
		//* --- and have muons from resonance decays in the wanted window --- 
		if(nsec != 2) {
		    printf(" AliMUON::StepManager: Strange resonance Decay into %d particles\n",nsec);
		    pMC->StopEvent();
		} else {
		    pMC->GetSecondary(0,ipart,x,p);
		    pt  = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
		    th0 = TMath::ATan2(pt,p[2])*kRaddeg;
		    pMC->GetSecondary(1,ipart,x,p);
		    pt  = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
		    th1 = TMath::ATan2(pt,p[2])*kRaddeg;
		    if(!(fAccMin < th0 && th0 < fAccMax) ||
		       !(fAccMin < th1 && th1 < fAccMax)) 
			pMC->StopEvent();
		}
	    }
	}
    }
}
void AliMUON::ReconstructClusters()
{
//
// Initialize the necessary correspondance table
//
    static const Int_t kMaxNpadx = 600;
    static const Int_t kMaxNpady = 600;
    Int_t elem[kMaxNpadx*2][kMaxNpady*2];
//
// Loop on chambers and on cathode planes
//
    for (Int_t ich=0;ich<10;ich++)
	for (Int_t icat=0;icat<2;icat++) {
	    //
	    // Get ready the current chamber stuff
	    //
	    AliMUONchamber*  iChamber= &(this->Chamber(ich));
	    AliMUONsegmentation*  segmentation = 
		iChamber->GetSegmentationModel(icat+1);
	    if (!segmentation) 
		continue;
	    TClonesArray *MUONdigits  = this->DigitsAddress(ich);
	    if (MUONdigits == 0) 
		continue;
	    cout << "Npx " << segmentation->Npx() 
		 << " Npy " << segmentation->Npy() << endl;
	    //      
	    // Ready the digits
	    //  
	    gAlice->ResetDigits();
	    gAlice->TreeD()->GetEvent(icat+1); // spurious +1 ...
	    Int_t ndigits = MUONdigits->GetEntriesFast();
	    if (ndigits == 0) 
		continue;
	    printf("Found %d digits for cathode %d in chamber %d \n",
		   ndigits,icat,ich+1);
	    AliMUONdigit  *mdig;
	    AliMUONRecCluster *Cluster;
	    //
	    // Build the correspondance table
	    //
	    memset(elem,0,sizeof(Int_t)*kMaxNpadx*kMaxNpady*4);
	    Int_t digit;
	    for (digit=0; digit<ndigits; digit++) 
	    {
		mdig    = (AliMUONdigit*)MUONdigits->UncheckedAt(digit);
		elem[kMaxNpadx+mdig->fPadX][kMaxNpady+mdig->fPadY] = digit+1;
		// because default is 0
	    }
	    //
	    // Declare some useful variables
	    //
	    Int_t Xlist[10];
	    Int_t Ylist[10];
	    Int_t Nlist;
	    Int_t nclust=0;
	    //
	    // loop over digits
	    //
	    for (digit=0;digit<ndigits;digit++) {
		mdig    = (AliMUONdigit*)MUONdigits->UncheckedAt(digit);
		//
		// if digit still available, start clustering
		//
		if (elem[kMaxNpadx+mdig->fPadX][kMaxNpady+mdig->fPadY]) {
		    Cluster = new AliMUONRecCluster(digit, ich, icat);
		    elem[kMaxNpadx+mdig->fPadX][kMaxNpady+mdig->fPadY]=0;
		    //
		    // loop over the current list of digits 
                    // and look for neighbours
		    //
		    for(Int_t clusDigit=Cluster->FirstDigitIndex();
			clusDigit!=Cluster->InvalidDigitIndex();
			clusDigit=Cluster->NextDigitIndex()) {
			AliMUONdigit* pDigit=(AliMUONdigit*)MUONdigits
			    ->UncheckedAt(clusDigit);
			segmentation->Neighbours(pDigit->fPadX,pDigit->fPadY, 
						 &Nlist, Xlist, Ylist);
			for (Int_t Ilist=0;Ilist<Nlist;Ilist++) {
			    if (elem[kMaxNpadx+Xlist[Ilist]][kMaxNpady
							    +Ylist[Ilist]]) {
				//
				// Add the digit at the end and reset the table
				//
				Cluster->AddDigit(elem[kMaxNpadx+Xlist[Ilist]][kMaxNpady+Ylist[Ilist]]-1);
				elem[kMaxNpadx+Xlist[Ilist]][kMaxNpady
							    +Ylist[Ilist]] =0;
			    } // if elem
			} // for Ilist
		    } // for pDigit
		    //
		    // Store the cluster (good time to do Cluster polishing)
		    //
		    segmentation->FitXY(Cluster,MUONdigits);
		    nclust ++;
		    AddRecCluster(ich,icat,Cluster);
		}
	    }
	    printf("===> %d Clusters\n",nclust); 
	} // for icat
}

 
//______________________________________________________________________________
void AliMUON::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliMUON.
      AliMUONchamber       *iChamber;
      AliMUONsegmentation  *segmentation;
      AliMUONresponse      *response;
      TClonesArray         *digitsaddress;
      
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDetector::Streamer(R__b);
      R__b >> fNclusters;
      R__b >> fClusters; // diff
      R__b >> fDchambers;
      R__b.ReadArray(fNdch);
      //
      R__b >> fAccCut;
      R__b >> fAccMin;
      R__b >> fAccMax; 
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
      }
      
   } else {
      R__b.WriteVersion(AliMUON::IsA());
      AliDetector::Streamer(R__b);
      R__b << fNclusters;
      R__b << fClusters; // diff
      R__b << fDchambers;
      R__b.WriteArray(fNdch, 10);
      //
      R__b << fAccCut;
      R__b << fAccMin;
      R__b << fAccMax; 
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
      }
   }
}
AliMUONcluster* AliMUON::FirstPad(AliMUONhit*  hit) 
{
//
    // Initialise the pad iterator
    // Return the address of the first padhit for hit
    TClonesArray *theClusters = Clusters();
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->fPHlast > 0) {
	sMaxIterPad=hit->fPHlast;
	sCurIterPad=hit->fPHfirst;
	return (AliMUONcluster*) fClusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}

AliMUONcluster* AliMUON::NextPad() 
{
    sCurIterPad++;
    if (sCurIterPad <= sMaxIterPad) {
	return (AliMUONcluster*) fClusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}

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
    for(Int_t i=0; i<10; i++) {
	fTcharges[i]  = charges[i];
	fTracks[i]    = tracks[i];
    }
}

ClassImp(AliMUONlist)
 
//____________________________________________________________________________
    AliMUONlist::AliMUONlist(Int_t rpad, Int_t *digits): 
	AliMUONdigit(digits)
{
    //
    // Creates a MUON digit list object
    //

    fRpad = rpad;
    fTrackList   = new TObjArray;
 
}
//_____________________________________________________________________________


ClassImp(AliMUONhit)
 
//___________________________________________
    AliMUONhit::AliMUONhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
	AliHit(shunt, track)
{
    fChamber=vol[0];
    fParticle=(Int_t) hits[0];
    fX=hits[1];
    fY=hits[2];
    fZ=hits[3];
    fTheta=hits[4];
    fPhi=hits[5];
    fTlength=hits[6];
    fEloss=hits[7];
    fPHfirst=(Int_t) hits[8];
    fPHlast=(Int_t) hits[9];
}
ClassImp(AliMUONreccluster) 

ClassImp(AliMUONRecCluster)

//_____________________________________________________________________
AliMUONRecCluster::AliMUONRecCluster()
{
    fDigits=0;
    fNdigit=-1;
}

AliMUONRecCluster::
AliMUONRecCluster(Int_t FirstDigit,Int_t Ichamber, Int_t Icathod)
{
    fX = 0.;
    fY = 0.;
    fDigits = new TArrayI(10);
    fNdigit=0;
    AddDigit(FirstDigit);
    fChamber=Ichamber;
    fCathod=Icathod;
}

void AliMUONRecCluster::AddDigit(Int_t Digit)
{
    if (fNdigit==fDigits->GetSize()) {
	//enlarge the list by hand!
	Int_t *array= new Int_t[fNdigit*2];
	for(Int_t i=0;i<fNdigit;i++)
	    array[i] = fDigits->At(i);
	fDigits->Adopt(fNdigit*2,array);
    }
    fDigits->AddAt(Digit,fNdigit);
    fNdigit++;
}


AliMUONRecCluster::~AliMUONRecCluster()
{
    if (fDigits)
	delete fDigits;
}

Int_t AliMUONRecCluster::FirstDigitIndex()
{
    fCurrentDigit=0;
    return fDigits->At(fCurrentDigit);
}

Int_t AliMUONRecCluster::NextDigitIndex()
{
    fCurrentDigit++;
    if (fCurrentDigit<fNdigit)
	return fDigits->At(fCurrentDigit);
    else 
	return InvalidDigitIndex();
}

Int_t AliMUONRecCluster::NDigits()
{
    return fNdigit;
}

void AliMUONRecCluster::Finish()
{
    // In order to reconstruct coordinates, one has to
    // get back to the digits which is not trivial here,
    // because we don't know where digits are stored!
    // Center Of Gravity, or other method should be
    // a property of AliMUON class!
}











