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

#include "AliRICH.h"
#include "AliRICHHitMap.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h" 

// Static variables for the pad-hit iterator routines
static Int_t sMaxIterPad=0;
static Int_t sCurIterPad=0;
static TClonesArray *fClusters2;
static TClonesArray *fHits2;
 
ClassImp(AliRICH)
    
//___________________________________________
AliRICH::AliRICH()
{
    fIshunt     = 0;
    fHits       = 0;
    fClusters   = 0;
    fNclusters  = 0;
    fDchambers  = 0;
    fRecClusters= 0;
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
    
    fHits       = new TClonesArray("AliRICHhit",1000  );
    fClusters   = new TClonesArray("AliRICHcluster",10000);
    fCerenkovs  = new TClonesArray("AliRICHCerenkov",20000);
    fNclusters  =  0;
    fIshunt     =  0;
    
    fNdch      = new Int_t[7];
    
    fDchambers = new TObjArray(7);
    
    Int_t i;
   
    for (i=0; i<7 ;i++) {
	(*fDchambers)[i] = new TClonesArray("AliRICHdigit",10000); 
	fNdch[i]=0;
    }
    
    fRecClusters=new TObjArray(7);
    for (i=0; i<7;i++)
	(*fRecClusters)[i] = new TObjArray(1000);
    
//   
// Transport angular cut
    fAccCut=0;
    fAccMin=2;
    fAccMax=9;
    
    SetMarkerColor(kRed);
}

//___________________________________________
AliRICH::~AliRICH()
{
    fIshunt  = 0;
    delete fHits;
    delete fClusters;
    delete fCerenkovs;
    for (Int_t i=0;i<7;i++) 
	delete (*fRecClusters)[i];
    delete fRecClusters;
    
}

//___________________________________________
void AliRICH::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
    TClonesArray &lhits = *fHits;
    new(lhits[fNhits++]) AliRICHhit(fIshunt,track,vol,hits);
}
//_____________________________________________________________________________
void AliRICH::AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs)
{
    TClonesArray &lcerenkovs = *fCerenkovs;
    new(lcerenkovs[fNcerenkovs++]) AliRICHCerenkov(fIshunt,track,vol,cerenkovs);
}
//___________________________________________
void AliRICH::AddCluster(Int_t *clhits)
{
    TClonesArray &lclusters = *fClusters;
    new(lclusters[fNclusters++]) AliRICHcluster(clhits);
}
//_____________________________________________________________________________
void AliRICH::AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Add a RICH digit to the list
    //

    TClonesArray &ldigits = *((TClonesArray*)(*fDchambers)[id]);
    new(ldigits[fNdch[id]++]) AliRICHdigit(tracks,charges,digits);
}


//_____________________________________________________________________________
void AliRICH::AddRecCluster(Int_t iCh, Int_t iCat, AliRICHRecCluster* Cluster)
{
    //
    // Add a RICH reconstructed cluster to the list
    //
    TObjArray* ClusterList = RecClusters(iCh,iCat);
    ClusterList->Add(Cluster);
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
    Chamber(0).SetChamberTransform(pos1[0],pos1[1],pos1[2],new TRotMatrix("rot993","rot993",90,0,70.69,90,19.30999,-90));
    Node = new TNode("RICH1","RICH1","S_RICH",pos1[0],pos1[1],pos1[2],"rot993");
    

    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    
    Float_t pos2[3]={171,470,0};
    Chamber(1).SetChamberTransform(pos2[0],pos2[1],pos2[2],new TRotMatrix("rot994","rot994",90,-20,90,70,0,0));
    Node = new TNode("RICH2","RICH2","S_RICH",pos2[0],pos2[1],pos2[2],"rot994");
    
    
    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos3[3]={0,500,0};
    Chamber(2).SetChamberTransform(pos3[0],pos3[1],pos3[2],new TRotMatrix("rot995","rot995",90,0,90,90,0,0));
    Node = new TNode("RICH3","RICH3","S_RICH",pos3[0],pos3[1],pos3[2],"rot995");
    

    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos4[3]={-171,470,0};
    Chamber(3).SetChamberTransform(pos4[0],pos4[1],pos4[2], new TRotMatrix("rot996","rot996",90,20,90,110,0,0));  
    Node = new TNode("RICH4","RICH4","S_RICH",pos4[0],pos4[1],pos4[2],"rot996");
    

    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos5[3]={161.3999,443.3999,-165.3};
    Chamber(4).SetChamberTransform(pos5[0],pos5[1],pos5[2],new TRotMatrix("rot997","rot997",90,340,108.1999,70,18.2,70));
    Node = new TNode("RICH5","RICH5","S_RICH",pos5[0],pos5[1],pos5[2],"rot997");
    
    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos6[3]={0., 471.9, -165.3,};
    Chamber(5).SetChamberTransform(pos6[0],pos6[1],pos6[2],new TRotMatrix("rot998","rot998",90,0,109.3099,90,19.30999,90));
    Node = new TNode("RICH6","RICH6","S_RICH",pos6[0],pos6[1],pos6[2],"rot998");
    
    
    Node->SetLineColor(kColorRICH);
    fNodes->Add(Node);
    Top->cd();
    Float_t pos7[3]={-161.399,443.3999,-165.3};
    Chamber(6).SetChamberTransform(pos7[0],pos7[1],pos7[2],new TRotMatrix("rot999","rot999",90,20,108.1999,110,18.2,110));
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
    
    sprintf(branchname,"%sCluster",GetName());
    if (fClusters   && gAlice->TreeH()) {
	gAlice->TreeH()->Branch(branchname,&fClusters, buffersize);
	printf("Making Branch %s for clusters\n",branchname);
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
// one branch for rec clusters
    for (i=0; i<7; i++) {
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
void AliRICH::SetTreeAddress()
{
    // Set branch address for the Hits and Digits Tree.
    char branchname[20];
    AliDetector::SetTreeAddress();
    
    TBranch *branch;
    TTree *treeH = gAlice->TreeH();
    TTree *treeD = gAlice->TreeD();
    
    if (treeH) {
	if (fClusters) {
	    branch = treeH->GetBranch("RICHCluster");
	    if (branch) branch->SetAddress(&fClusters);
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
}
//___________________________________________
void AliRICH::ResetHits()
{
    // Reset number of clusters and the cluster array for this detector
    AliDetector::ResetHits();
    fNclusters = 0;
    if (fClusters) fClusters->Clear();
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
void AliRICH::ResetRecClusters()
{
    //
    // Reset the rec clusters
    //
    for ( int i=0;i<7;i++ ) {
	if ((*fRecClusters)[i])   (*fRecClusters)[i]->Clear();
    }
}
//___________________________________________

void AliRICH::SetPADSIZ(Int_t id, Int_t isec, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliRICHchamber*) (*fChambers)[i])  ->SetPADSIZ(isec,p1,p2);
    ((AliRICHchamber*) (*fChambers)[i+1])->SetPADSIZ(isec,p1,p2);
}

//___________________________________________
void AliRICH::SetMUCHSP(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliRICHchamber*) (*fChambers)[i])->SetMUCHSP(p1);
    ((AliRICHchamber*) (*fChambers)[i+1])->SetMUCHSP(p1);
}

//___________________________________________
void AliRICH::SetMUSIGM(Int_t id, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliRICHchamber*) (*fChambers)[i])->SetMUSIGM(p1,p2);
    ((AliRICHchamber*) (*fChambers)[i+1])->SetMUSIGM(p1,p2);
}

//___________________________________________
void AliRICH::SetRSIGM(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliRICHchamber*) (*fChambers)[i])->SetRSIGM(p1);
    ((AliRICHchamber*) (*fChambers)[i+1])->SetRSIGM(p1);
}

//___________________________________________
void AliRICH::SetMAXADC(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliRICHchamber*) (*fChambers)[i])->SetMAXADC(p1);
    ((AliRICHchamber*) (*fChambers)[i+1])->SetMAXADC(p1);
}

//___________________________________________
void AliRICH::SetSMAXAR(Float_t p1)
{
    fMaxStepGas=p1;
}

//___________________________________________
void AliRICH::SetSMAXAL(Float_t p1)
{
    fMaxStepAlu=p1;
}

//___________________________________________
void AliRICH::SetDMAXAR(Float_t p1)
{
    fMaxDestepGas=p1;
}

//___________________________________________
void AliRICH::SetDMAXAL(Float_t p1)
{
    fMaxDestepAlu=p1;
}
//___________________________________________
void AliRICH::SetRICHACC(Bool_t acc, Float_t angmin, Float_t angmax)
{
    fAccCut=acc;
    fAccMin=angmin;
    fAccMax=angmax;
}
//___________________________________________
void   AliRICH::SetSegmentationModel(Int_t id, Int_t isec, AliRICHsegmentation *segmentation)
{
    ((AliRICHchamber*) (*fChambers)[id])->SegmentationModel(isec, segmentation);

}
//___________________________________________
void   AliRICH::SetResponseModel(Int_t id, Response_t res,  AliRICHresponse *response)
{
    ((AliRICHchamber*) (*fChambers)[id])->ResponseModel(res, response);
}

void   AliRICH::SetNsec(Int_t id, Int_t nsec)
{
    ((AliRICHchamber*) (*fChambers)[id])->SetNsec(nsec);
}


//___________________________________________

void AliRICH::StepManager()
{
    printf("Dummy version of RICH step -- it should never happen!!\n");
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
		
		// Check angular acceptance
		//* --- and have muons from resonance decays in the wanted window --- 
		if(nsec != 2) {
		    printf(" AliRICH::StepManager: Strange resonance Decay into %d particles\n",nsec);
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
void AliRICH::ReconstructClusters()
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
    for (Int_t ich=0;ich<7;ich++)
	for (Int_t icat=0;icat<1;icat++) {
	    //
	    // Get ready the current chamber stuff
	    //
	    
	    printf ("Olarilole");
	    AliRICHchamber*  iChamber= &(this->Chamber(ich));
	    AliRICHsegmentation*  segmentation = 
		iChamber->GetSegmentationModel(icat+1);
	    if (!segmentation) 
		continue;
	    TClonesArray *RICHdigits  = this->DigitsAddress(ich);
	    if (RICHdigits == 0) 
		continue;
	    cout << "Npx " << segmentation->Npx() 
		 << " Npy " << segmentation->Npy() << endl;
	    //      
	    // Ready the digits
	    //  
	    gAlice->ResetDigits();
	    gAlice->TreeD()->GetEvent(icat+1); // spurious +1 ...
	    Int_t ndigits = RICHdigits->GetEntriesFast();
	    if (ndigits == 0) 
		continue;
	    printf("Found %d digits for cathode %d in chamber %d \n",
		   ndigits,icat,ich+1);
	    AliRICHdigit  *mdig;
	    AliRICHRecCluster *Cluster;
	    //
	    // Build the correspondance table
	    //
	    memset(elem,0,sizeof(Int_t)*kMaxNpadx*kMaxNpady*4);
	    Int_t digit;
	    for (digit=0; digit<ndigits; digit++) 
	    {
		mdig    = (AliRICHdigit*)RICHdigits->UncheckedAt(digit);
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
		mdig    = (AliRICHdigit*)RICHdigits->UncheckedAt(digit);
		//
		// if digit still available, start clustering
		//
		if (elem[kMaxNpadx+mdig->fPadX][kMaxNpady+mdig->fPadY]) {
		    Cluster = new AliRICHRecCluster(digit, ich, icat);
		    elem[kMaxNpadx+mdig->fPadX][kMaxNpady+mdig->fPadY]=0;
		    //
		    // loop over the current list of digits 
                    // and look for neighbours
		    //
		    for(Int_t clusDigit=Cluster->FirstDigitIndex();
			clusDigit!=Cluster->InvalidDigitIndex();
			clusDigit=Cluster->NextDigitIndex()) {
			AliRICHdigit* pDigit=(AliRICHdigit*)RICHdigits
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
		    segmentation->FitXY(Cluster,RICHdigits);
		    nclust ++;
		    AddRecCluster(ich,icat,Cluster);
		}
	    }
	    printf("===> %d Clusters\n",nclust); 
	} // for icat
}


//______________________________________________________________________________
void AliRICH::Streamer(TBuffer &R__b)
{
    // Stream an object of class AliRICH.
    AliRICHchamber       *iChamber;
    AliRICHsegmentation  *segmentation;
    AliRICHresponse      *response;
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
	for (Int_t i =0; i<7; i++) {
	    iChamber=(AliRICHchamber*) (*fChambers)[i];
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
	    response=iChamber->GetResponseModel(mip);
	    response->Streamer(R__b);	  
	    response=iChamber->GetResponseModel(cerenkov);
	    response->Streamer(R__b);	  
	    
	    digitsaddress=(TClonesArray*) (*fDchambers)[i];
	    digitsaddress->Streamer(R__b);
	}
      
    } else {
	R__b.WriteVersion(AliRICH::IsA());
	AliDetector::Streamer(R__b);
	R__b << fNclusters;
	R__b << fClusters; // diff
	R__b << fDchambers;
	R__b.WriteArray(fNdch, 7);
	//
	R__b << fAccCut;
	R__b << fAccMin;
	R__b << fAccMax; 
	//   
	R__b << fChambers;
//  Stream chamber related information
	for (Int_t i =0; i<7; i++) {
	    iChamber=(AliRICHchamber*) (*fChambers)[i];
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
	    response=iChamber->GetResponseModel(mip);
	    response->Streamer(R__b);
	    response=iChamber->GetResponseModel(cerenkov);
	    response->Streamer(R__b);
	    
	    digitsaddress=(TClonesArray*) (*fDchambers)[i];
	    digitsaddress->Streamer(R__b);
	}
    }
}
AliRICHcluster* AliRICH::FirstPad(AliRICHhit*  hit,TClonesArray *clusters ) 
{
//
    // Initialise the pad iterator
    // Return the address of the first padhit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->fPHlast > 0) {
	sMaxIterPad=hit->fPHlast;
	sCurIterPad=hit->fPHfirst;
	return (AliRICHcluster*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
    
}

AliRICHcluster* AliRICH::NextPad(TClonesArray *clusters) 
{
    sCurIterPad++;
    if (sCurIterPad <= sMaxIterPad) {
	return (AliRICHcluster*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}

ClassImp(AliRICHcluster)
    
//___________________________________________
AliRICHcluster::AliRICHcluster(Int_t *clhits)
{
    fHitNumber=clhits[0];
    fCathode=clhits[1];
    fQ=clhits[2];
    fPadX=clhits[3];
    fPadY=clhits[4];
    fQpad=clhits[5];
    fRSec=clhits[6];
}
ClassImp(AliRICHdigit)
//_____________________________________________________________________________
AliRICHdigit::AliRICHdigit(Int_t *digits)
{
    //
    // Creates a RICH digit object to be updated
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fSignal      = digits[2];
    
}
//_____________________________________________________________________________
AliRICHdigit::AliRICHdigit(Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Creates a RICH digit object
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fSignal      = digits[2];
    for(Int_t i=0; i<10; i++) {
	fTcharges[i]  = charges[i];
	fTracks[i]    = tracks[i];
    }
}

ClassImp(AliRICHlist)
    
//____________________________________________________________________________
AliRICHlist::AliRICHlist(Int_t ich, Int_t *digits): 
    AliRICHdigit(digits)
{
    //
    // Creates a RICH digit list object
    //
    
    fChamber = ich;
    fTrackList   = new TObjArray;
    
}
//_____________________________________________________________________________


ClassImp(AliRICHhit)
    
//___________________________________________
AliRICHhit::AliRICHhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
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
}
ClassImp(AliRICHreccluster)

ClassImp(AliRICHCerenkov)
//___________________________________________
AliRICHCerenkov::AliRICHCerenkov(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
    AliHit(shunt, track)
{
    fChamber=vol[0];
    fX=hits[1];
    fY=hits[2];
    fZ=hits[3];
    fTheta=hits[4];
    fPhi=hits[5];
    fTlength=hits[6];
    fPHfirst=(Int_t) hits[8];
    fPHlast=(Int_t) hits[9];
}

ClassImp(AliRICHRecCluster)
    
//_____________________________________________________________________
AliRICHRecCluster::AliRICHRecCluster()
{
    fDigits=0;
    fNdigit=-1;
}

AliRICHRecCluster::AliRICHRecCluster(Int_t FirstDigit,Int_t Ichamber, Int_t Icathod)
{
    fX = 0.;
    fY = 0.;
    fDigits = new TArrayI(10);
    fNdigit=0;
    AddDigit(FirstDigit);
    fChamber=Ichamber;
    fCathod=Icathod;
}

void AliRICHRecCluster::AddDigit(Int_t Digit)
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


AliRICHRecCluster::~AliRICHRecCluster()
{
    if (fDigits)
	delete fDigits;
}

Int_t AliRICHRecCluster::FirstDigitIndex()
{
    fCurrentDigit=0;
    return fDigits->At(fCurrentDigit);
}

Int_t AliRICHRecCluster::NextDigitIndex()
{
    fCurrentDigit++;
    if (fCurrentDigit<fNdigit)
	return fDigits->At(fCurrentDigit);
    else 
	return InvalidDigitIndex();
}

Int_t AliRICHRecCluster::NDigits()
{
    return fNdigit;
}

void AliRICHRecCluster::Finish()
{
    // In order to reconstruct coordinates, one has to
    // get back to the digits which is not trivial here,
    // because we don't know where digits are stored!
    // Center Of Gravity, or other method should be
    // a property of AliRICH class!
}



void AliRICH::Digitise(Int_t nev,Option_t *option,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
    
    static Bool_t first=true;
    static TTree *TH1;
    static TFile *File;
    char *Add = strstr(option,"Add");

    AliRICHchamber*  iChamber;
    AliRICHsegmentation*  segmentation;

    
    Int_t trk[50];
    Int_t chtrk[50];  
    TObjArray *list=new TObjArray;
    Int_t digits[3]; 
    
    AliRICH *RICH  = (AliRICH *) gAlice->GetDetector("RICH");
    AliRICHHitMap* HitMap[10];
    for (Int_t i=0; i<10; i++) {HitMap[i]=0;}
    if (Add ) {
	if(first) {
	    fFileName=filename;
	    cout<<"filename"<<fFileName<<endl;
	    File=new TFile(fFileName);
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
	    fHits2     = new TClonesArray("AliRICHhit",1000  );
	    fClusters2 = new TClonesArray("AliRICHcluster",10000);
	    first=false;
	}
	File->cd();
	File->ls();
	// Get Hits Tree header from file
	if(fHits2) fHits2->Clear();
	if(fClusters2) fClusters2->Clear();
	if(TH1) delete TH1;
	TH1=0;
	//
	char treeName[20];
	sprintf(treeName,"TreeH%d",nev);
	TH1 = (TTree*)gDirectory->Get(treeName);
	if (!TH1) {
	    printf("ERROR: cannot find Hits Tree for event:%d\n",nev);
	}
	// Set branch addresses
	TBranch *branch;
	char branchname[20];
	sprintf(branchname,"%s",GetName());
	if (TH1 && fHits2) {
	    branch = TH1->GetBranch(branchname);
	    if (branch) branch->SetAddress(&fHits2);
	}
	if (TH1 && fClusters2) {
	    branch = TH1->GetBranch("RICHCluster");
	    if (branch) branch->SetAddress(&fClusters2);
	}
    }
    //
    // loop over cathodes
    //
    AliRICHHitMap* hm;
    
    for (int icat=0; icat<1; icat++) { 
	for (Int_t i=0; i<7; i++) {
	    if (HitMap[i]) {
		hm=HitMap[i];
		delete hm;
		HitMap[i]=0;
	    }
	}
	Int_t counter=0;
	for (Int_t i =0; i<7; i++) {
	    iChamber=(AliRICHchamber*) (*fChambers)[i];
	    if (iChamber->Nsec()==1 && icat==1) {
		continue;
	    } else {
		segmentation=iChamber->GetSegmentationModel(icat+1);
	    }
	    HitMap[i] = new AliRICHHitMapA1(segmentation, list);
	}
	printf("Start loop over tracks \n");     
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
	    for(AliRICHhit* mHit=(AliRICHhit*)RICH->FirstHit(-1); 
		mHit;
		mHit=(AliRICHhit*)RICH->NextHit()) 
	    {
		Int_t   nch   = mHit->fChamber-1;  // chamber number
		if (nch >7) continue;
		iChamber = &(RICH->Chamber(nch));
				
//
// Loop over pad hits
		for (AliRICHcluster* mPad=
			 (AliRICHcluster*)RICH->FirstPad(mHit,fClusters);
		     mPad;
		     mPad=(AliRICHcluster*)RICH->NextPad(fClusters))
		{
		    Int_t cathode  = mPad->fCathode;    // cathode number
		    Int_t ipx      = mPad->fPadX;       // pad number on X
		    Int_t ipy      = mPad->fPadY;       // pad number on Y
		    Int_t iqpad    = mPad->fQpad;       // charge per pad
//
//
		    
		    if (cathode != (icat+1)) continue;
		    // fill the info array
		    Float_t thex, they;
		    segmentation=iChamber->GetSegmentationModel(cathode);
		    segmentation->GetPadCxy(ipx,ipy,thex,they);
		    TVector *trinfo_p= new TVector(2);
		    TVector &trinfo = *trinfo_p;
		    trinfo(0)=(Float_t)track;
		    trinfo(1)=(Float_t)iqpad;
		    
		    digits[0]=ipx;
		    digits[1]=ipy;
		    digits[2]=iqpad;
		    
		    AliRICHlist* pdigit;
		    // build the list of fired pads and update the info
		    if (!HitMap[nch]->TestHit(ipx, ipy)) {
			list->AddAtAndExpand(
			    new AliRICHlist(nch,digits),counter);
			HitMap[nch]->SetHit(ipx, ipy, counter);
			counter++;
			pdigit=(AliRICHlist*)list->At(list->GetLast());
			// list of tracks
			TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			trlist->Add(&trinfo);
		    } else {
			pdigit=(AliRICHlist*) HitMap[nch]->GetHit(ipx, ipy);
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
			    printf("Attention - nptracks > 2  %d \n",nptracks);
			    printf("cat,nch,ix,iy %d %d %d %d  \n",icat+1,nch,ipx,ipy);
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
	
	Int_t nentr1=list->GetEntriesFast();
	printf(" \n counter, nentr1 %d %d\n",counter,nentr1);
	
	// open the file with background
	
	if (Add ) {
	    ntracks =(Int_t)TH1->GetEntries();
	    printf("background - icat,ntracks1  %d %d\n",icat,ntracks);
	    printf("background - Start loop over tracks \n");     
//
//   Loop over tracks
//
	    for (Int_t track=0; track<ntracks; track++) {

		if (fHits2)       fHits2->Clear();
		if (fClusters2)   fClusters2->Clear();
		
		TH1->GetEvent(track);
//
//   Loop over hits
		AliRICHhit* mHit;
		for(int i=0;i<fHits2->GetEntriesFast();++i) 
		{
		    mHit=(AliRICHhit*) (*fHits2)[i];
		    Int_t   nch   = mHit->fChamber-1;  // chamber number
		    if (nch >9) continue;
		    iChamber = &(RICH->Chamber(nch));
		    Int_t rmin = (Int_t)iChamber->RInner();
		    Int_t rmax = (Int_t)iChamber->ROuter();
//
// Loop over pad hits
		    for (AliRICHcluster* mPad=
			     (AliRICHcluster*)RICH->FirstPad(mHit,fClusters2);
			 mPad;
			 mPad=(AliRICHcluster*)RICH->NextPad(fClusters2))
		    {
			Int_t cathode  = mPad->fCathode;    // cathode number
			Int_t ipx      = mPad->fPadX;       // pad number on X
			Int_t ipy      = mPad->fPadY;       // pad number on Y
			Int_t iqpad    = mPad->fQpad;       // charge per pad
			if (track==3 && nch==0 && icat==0) printf("bgr - track,iqpad,ipx,ipy %d %d %d %d\n",track,iqpad,ipx,ipy);
//
//
			if (cathode != (icat+1)) continue;
			// fill the info array
			Float_t thex, they;
			segmentation=iChamber->GetSegmentationModel(cathode);
			segmentation->GetPadCxy(ipx,ipy,thex,they);
			Float_t rpad=TMath::Sqrt(thex*thex+they*they);
			if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;
			
			TVector *trinfo_p;
			trinfo_p = new TVector(2);
			TVector &trinfo = *trinfo_p;
			trinfo(0)=-1;  // tag background
			trinfo(1)=-1;
			
			digits[0]=ipx;
			digits[1]=ipy;
			digits[2]=iqpad;

			
			if (track <4 && icat==0 && nch==0)
			    printf("bgr - HitMap[nch]->TestHit(ipx, ipy),track %d %d\n",
				   HitMap[nch]->TestHit(ipx, ipy),track);
			AliRICHlist* pdigit;
			// build the list of fired pads and update the info
			if (!HitMap[nch]->TestHit(ipx, ipy)) {
			    list->AddAtAndExpand(new AliRICHlist(nch,digits),counter);
			    
			    HitMap[nch]->SetHit(ipx, ipy, counter);
			    counter++;
			    printf("bgr new elem in list - counter %d\n",counter);
			    
			    pdigit=(AliRICHlist*)list->At(list->GetLast());
			    // list of tracks
			    TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			    trlist->Add(&trinfo);
			} else {
			    pdigit=(AliRICHlist*) HitMap[nch]->GetHit(ipx, ipy);
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
	    Int_t nentr2=list->GetEntriesFast();
	    printf(" \n counter2, nentr2 %d %d \n",counter,nentr2);
	    TTree *fAli=gAlice->TreeK();
	    if (fAli) File =fAli->GetCurrentFile();
	    File->cd();
	} // if Add	
	
	Int_t tracks[10];
	Int_t charges[10];
	cout<<"start filling digits \n "<<endl;
	Int_t nentries=list->GetEntriesFast();
	printf(" \n \n nentries %d \n",nentries);

	// start filling the digits
	
	for (Int_t nent=0;nent<nentries;nent++) {
	    AliRICHlist *address=(AliRICHlist*)list->At(nent);
	    if (address==0) continue; 
	    Int_t ich=address->fChamber;
	    Int_t q=address->fSignal; 
	    iChamber=(AliRICHchamber*) (*fChambers)[ich];
	    // add white noise and do zero-suppression and signal truncation
	    Float_t MeanNoise = gRandom->Gaus(1, 0.2);
	    Float_t ZeroSupp=5*MeanNoise;
	    Float_t Noise     = gRandom->Gaus(0, MeanNoise);
	    q+=(Int_t)Noise; 
	    if ( q <= ZeroSupp) continue;
	    digits[0]=address->fPadX;
	    digits[1]=address->fPadY;
	    digits[2]=q;
	    
	    TObjArray* trlist=(TObjArray*)address->TrackList();
	    Int_t nptracks=trlist->GetEntriesFast();
	    
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
		TVector &pp  =*pp_p;
		tracks[tr]=Int_t(pp(0));
		charges[tr]=Int_t(pp(1));
	    }      //end loop over list of tracks for one pad
	    if (nptracks < 10 ) {
		for (Int_t i=nptracks; i<10; i++) {
		    tracks[i]=0;
		    charges[i]=0;
		}
	    }
	    // fill digits
	    RICH->AddDigits(ich,tracks,charges,digits);
	    
	    delete address;
	}
	cout<<"I'm out of the loops for digitisation"<<endl;
	gAlice->TreeD()->Fill();
	TTree *TD=gAlice->TreeD();
	Stat_t ndig=TD->GetEntries();
	cout<<"number of digits  "<<ndig<<endl;
	TClonesArray *fDch;
	for (int i=0;i<7;i++) {
	    fDch= RICH->DigitsAddress(i);
	    int ndig=fDch->GetEntriesFast();
	    printf (" i, ndig %d %d \n",i,ndig);
	}
	RICH->ResetDigits();
	
	list->Clear();
	
    } //end loop over cathodes
    char hname[30];
    sprintf(hname,"TreeD%d",nev);
    gAlice->TreeD()->Write(hname);
}









