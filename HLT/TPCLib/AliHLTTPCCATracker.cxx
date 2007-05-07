// @(#) $Id$
/** \class AliHLTTPCCATracker
<pre>
//_____________________________________________________________
// AliHLTTPCCATracker
// !
// !
// !
// !
// !
//
// Author: Ivan Kisel
// Copyright &copy ALICE HLT Group
</pre>
*/

#define WRITETRACKS 1

#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCCATracker.h"
#include "AliHLTTPCTrackSegmentData.h"

#include "TH1.h"
#include "TH2.h"

#include "TFile.h"
#include "TMath.h"

#if __GNUC__ >= 3
using namespace std;
#endif

#include "TROOT.h"


#include "TCanvas.h"
#include "TMarker.h"
#include "TLine.h"
#include "TStyle.h"
#include "TApplication.h"

#include "AliHLTTPCConfMapPoint.h"
#include "AliHLTTPCTrackArray.h"
#include <iostream>

ClassImp(AliHLTTPCCATracker)

// ----------------------------------------------------------------------------------
AliHLTTPCCATracker::AliHLTTPCCATracker()
{
  //Default constructor
  fVertex = NULL;
  fTrack = NULL;
  fHit = NULL;
  fBench = (Bool_t)true;
  fVertexConstraint = (Bool_t)true;
  fParamSet[0]=0;
  fParamSet[1]=0;
  //JMT 2006/11/13
  fOutputPtr= NULL;
  fOutputNTracks=0;
  fOutputSize=0;
  fDRAW = 0;
  fAsk = 0;
}


// ----------------------------------------------------------------------------------
AliHLTTPCCATracker::AliHLTTPCCATracker(const AliHLTTPCCATracker &x )
{
  //Copy constructor
  fVertex = x.fVertex;
  fTrack = x.fTrack;
  fHit = x.fHit;
  fBench = x.fBench;
  fVertexConstraint = x.fVertexConstraint;
  fParamSet[0]=x.fParamSet[0];
  fParamSet[1]=x.fParamSet[1];
  //JMT 2006/11/13
  fOutputPtr= x.fOutputPtr;
  fOutputNTracks=x.fOutputNTracks;
  fOutputSize=x.fOutputSize;
  fDRAW = x.fDRAW;
  fAsk = x.fAsk;
}

// ----------------------------------------------------------------------------------
AliHLTTPCCATracker & AliHLTTPCCATracker::operator=( const AliHLTTPCCATracker &x )
{
  //Copy operator
  fVertex = x.fVertex;
  fTrack = x.fTrack;
  fHit = x.fHit;
  fBench = x.fBench;
  fVertexConstraint = x.fVertexConstraint;
  fParamSet[0]=x.fParamSet[0];
  fParamSet[1]=x.fParamSet[1];
  //JMT 2006/11/13
  fOutputPtr= x.fOutputPtr;
  fOutputNTracks=x.fOutputNTracks;
  fOutputSize=x.fOutputSize;
  fDRAW = x.fDRAW;
  fAsk = x.fAsk;
  return *this;
}

// ----------------------------------------------------------------------------------
AliHLTTPCCATracker::~AliHLTTPCCATracker()
{
  // Destructor.
  if(fHit) {
    delete [] fHit;
  }
  if(fTrack) {
    delete fTrack;
  }
  //JMT 2006/11/13
  if(fOutputPtr){
    fOutputPtr = NULL;
    delete fOutputPtr;
  }

}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CACreateHistos()
{
    // IK histogramming
    static bool firstCall = kTRUE;
    if (firstCall){
      firstCall = kFALSE;

      fHistEvent = 0;

      TDirectory *curdir = gDirectory;
      TDirectory *histodir = gROOT->mkdir("histodir");
      histodir->cd();
      
      fHistNClusters  = new TH1F("h_NClusters", "Number of clusters in event", 100, 0.0, 1000.0);

      fHistClustersXY[0] = new TH2F("h_ClustersXY_0", "Clusters in XY plane Ev 0", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[1] = new TH2F("h_ClustersXY_1", "Clusters in XY plane Ev 1", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[2] = new TH2F("h_ClustersXY_2", "Clusters in XY plane Ev 2", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[3] = new TH2F("h_ClustersXY_3", "Clusters in XY plane Ev 3", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[4] = new TH2F("h_ClustersXY_4", "Clusters in XY plane Ev 4", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[5] = new TH2F("h_ClustersXY_5", "Clusters in XY plane Ev 5", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[6] = new TH2F("h_ClustersXY_6", "Clusters in XY plane Ev 6", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[7] = new TH2F("h_ClustersXY_7", "Clusters in XY plane Ev 7", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[8] = new TH2F("h_ClustersXY_8", "Clusters in XY plane Ev 8", 100, 50.0, 300.0, 100, -100.0, 100.0);
      fHistClustersXY[9] = new TH2F("h_ClustersXY_9", "Clusters in XY plane Ev 9", 100, 50.0, 300.0, 100, -100.0, 100.0);

      fListHisto = gDirectory->GetList();
      curdir->cd();

      if( fDRAW ){

      fMyapp = new TApplication("myapp", 0, 0);
      
      gStyle->SetCanvasBorderMode(0);
      gStyle->SetCanvasBorderSize(1);
      gStyle->SetCanvasColor(0);


      fYX = new TCanvas ("YX", "YX (Pad-Row) window", -10, 0, 680, 400);
      fYX->Range(-50.0, 80.0, 50.0, 250.0);
      fYX->Clear();
      fYX->Draw();

      fZX = new TCanvas ("ZX", "ZX window", -700, 0, 680, 400);
      fZX->Range(-250.0, 80.0, 250.0, 250.0);
      fZX->Clear();
      fZX->Draw();

#ifdef XXX
      char symbol;
      if (fAsk){
        do{
          std::cin.get(symbol);
          if (symbol == 'r')
            fAsk = 0;
        } while (symbol != '\n');
      }
#endif //XXX

} //fDRAW

    } 
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAWriteHistos()
{
    // IK Open the output file and write histogramms

    TFile *outfile = new TFile("ConfMapper_histo.root", "RECREATE");
    TIter hiter(fListHisto);
    while (TObject *obj = hiter()) obj->Write();
    outfile->Close();
}

// ----------------------------------------------------------------------------------
Bool_t AliHLTTPCCATracker::ReadHits(UInt_t count, AliHLTTPCSpacePointData* hits )
{

  //read hits
  Int_t nhit=(Int_t)count; 

  vector<AliHLTTPCCAHit> fVecHits; 
  fVecHits.clear();

  for (Int_t i=0;i<nhit;i++)
    {	
      fHit[i+fClustersUnused].Reset();
      fHit[i+fClustersUnused].ReadHits(&(hits[i]));
    }

  fClustersUnused += nhit;

  return true;
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAInitialize()
{
  //!
  fVecHits.clear();
  fVecPatchTracks.clear();
  fVecSliceTracks.clear();

  return;
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAReadPatchHits(Int_t patch, UInt_t count, AliHLTTPCSpacePointData* hits )
{
  //read hits
  Int_t nhit=(Int_t)count; 

  fPatchFirstHitInd = fVecHits.size();
  fPatchLastHitInd  = fPatchFirstHitInd + nhit - 1;

  for (Int_t i=0;i<nhit;i++)
    {	
      AliHLTTPCSpacePointData* pSP = &(hits[i]);

      AliHLTTPCCAHit hit;
      hit.fX = pSP->fX;
      hit.fY = pSP->fY;
      hit.fZ = pSP->fZ;

      hit.fErrx = 1.0;//pSP->fSigmaY2;
      if (patch <= 3)
	hit.fErry = 10.0;//1.0;//pSP->fSigmaY2;
      else
	hit.fErry = 10.0;//pSP->fSigmaY2;
      hit.fErrz = 1.0;//pSP->fSigmaZ2;

      hit.fIndex = i;
      hit.fCounter = 0;

      fVecHits.push_back(hit);
    }

  if( fDRAW ){

  TMarker *mhit = new TMarker(0.0, 0.0, 6);
  mhit->SetMarkerColor(kBlue);

  for (Int_t i = fPatchFirstHitInd; i <= fPatchLastHitInd; i++)
    {
      fYX->cd();
      mhit->DrawMarker(fVecHits[i].fY, fVecHits[i].fX);
      
      fZX->cd();
      mhit->DrawMarker(fVecHits[i].fZ, fVecHits[i].fX);
    }

  delete mhit;

  fYX->cd(); fYX->Update(); fYX->Draw(); 
  fZX->cd(); fZX->Update(); fZX->Draw();     

  } //fDRAW

  return;
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAFindPatchTracks(Int_t patch)
{
//!

//if(fDRAW){

TLine *line = new TLine();
line->SetLineColor(kBlack);

//} //fDRAW

  Double_t dist01, dist12, dist02, mindist, chi2, minchi2;
  AliHLTTPCCAHit *ph0, *ph1, *ph2, *phbest, *next;

  Double_t x, y, z, ty, tz, fCovy, fCovty, fCovyty, fCovz, fCovtz, fCovztz, ye, ze, dx, dy, dz, sy2, sz2; 

  for (Int_t i1 = fPatchFirstHitInd; i1 < fPatchLastHitInd; i1++) // <nhit-1
    {	
      ph1  = &(fVecHits[i1]);
      next = &(fVecHits[i1+1]);

      mindist = 20.0;
      minchi2 = 1.0;
      phbest = NULL;

      if (ph1->fCounter != 0){ // the previous hit exists
	x = ph1->fX;
	y = ph1->fY;
	z = ph1->fZ;

	dx = (ph1->fX-ph0->fX); // assume different x !!!
	if (TMath::Abs(dx) < 0.0001) dx = 0.0001;
	ty = (ph1->fY-ph0->fY)/dx;
	tz = (ph1->fZ-ph0->fZ)/dx;

	fCovy   = ph1->fErry*ph1->fErry;
	fCovty  = (ph0->fErry*ph0->fErry + ph1->fErry*ph1->fErry)/(dx*dx);
	fCovyty = (ph1->fErry*ph1->fErry)/dx;

	fCovz   = ph1->fErrz*ph1->fErrz;
	fCovtz  = (ph0->fErrz*ph0->fErrz + ph1->fErrz*ph1->fErrz)/(dx*dx);
	fCovztz = (ph1->fErrz*ph1->fErrz)/dx;

      }

      for (Int_t i2 = i1+1; i2 <= fPatchLastHitInd; i2++)
	{	
	  ph2 = &(fVecHits[i2]);
	  if (ph2->fX < ph1->fX) continue; // no backward direction!

	  dist12 = TMath::Sqrt((ph1->fX-ph2->fX)*(ph1->fX-ph2->fX) + (ph1->fY-ph2->fY)*(ph1->fY-ph2->fY) + (ph1->fZ-ph2->fZ)*(ph1->fZ-ph2->fZ));
	  if (dist12 < mindist){
	    if (ph1->fCounter == 0){
	      mindist = dist12;
	      phbest  = ph2;
	    }
	    else{
	      // calculate chi2 distance between the hit and the line
	      dx = (ph2->fX-x); // assume different x !!!
	      if (TMath::Abs(dx) < 0.0001) dx = 0.0001;
	      
	      ye = y + ty*dx;
	      ze = z + tz*dx;

	      sy2 = fCovy + 2.0*fCovyty*dx + fCovty*dx*dx;
	      sz2 = fCovz + 2.0*fCovztz*dx + fCovtz*dx*dx;

	      dy = ph2->fY - ye;
	      dz = ph2->fZ - ze;
	      
	      chi2 = (dy*dy)/(sy2+ph2->fErry*ph2->fErry) + (dz*dz)/(sz2+ph2->fErrz*ph2->fErrz);

	      // check the straight line model
	      if (chi2 < minchi2){
		mindist = dist12;
		phbest  = ph2;	
	      }

	    }
	  }
	}

      if (phbest != NULL){// closest hit found 

	// exchange jbest with the next hit
	Double_t xtmp = next->fX;
	Double_t ytmp = next->fY;
	Double_t ztmp = next->fZ;

	Double_t errxtmp = next->fErrx;
	Double_t errytmp = next->fErry;
	Double_t errztmp = next->fErrz;

	Int_t indextmp = next->fIndex;
	Int_t countertmp = next->fCounter;

	next->fX = phbest->fX;
	next->fY = phbest->fY;
	next->fZ = phbest->fZ;

	next->fErrx = phbest->fErrx;
	next->fErry = phbest->fErry;
	next->fErrz = phbest->fErrz;

	next->fIndex = phbest->fIndex;
	next->fCounter = phbest->fCounter;

	phbest->fX = xtmp;
	phbest->fY = ytmp;
	phbest->fZ = ztmp;

	phbest->fErrx = errxtmp;
	phbest->fErry = errytmp;
	phbest->fErrz = errztmp;

	phbest->fIndex = indextmp;
	phbest->fCounter = countertmp;

	if (ph1->fCounter == 0)
	  ph1->fCounter = 1;
	next->fCounter = ph1->fCounter+1;

	ph0 = ph1;
	dist01 = mindist;

if(fDRAW){

	fYX->cd();
	line->DrawLine(ph1->fY, ph1->fX, next->fY, next->fX);
	
	fZX->cd();
	line->DrawLine(ph1->fZ, ph1->fX, next->fZ, next->fX);

}//fDRAW

      }

    }

  if(fDRAW){

  delete line;

  fYX->cd(); fYX->Update(); fYX->Draw(); 
  fZX->cd(); fZX->Update(); fZX->Draw();     

}//fDRAW
  
  AliHLTTPCCATrack tr;
  tr.fNhits = 0; 
  tr.fChi2  = 0.0; 
  tr.fNdf   = -4; 
  tr.fVecIHits.clear();

  fPatchFirstTrackInd = fVecPatchTracks.size(); 
  fPatchLastTrackInd  = fPatchFirstTrackInd - 1;

  AliHLTTPCCAHit *myphf; 

  // collect tracks
  for (Int_t ihr = fPatchLastHitInd; ihr >= fPatchFirstHitInd; ihr--){

    AliHLTTPCCAHit* phit = &(fVecHits[ihr]);

    if ((phit->fCounter < 5)&&(tr.fNhits == 0))
      continue; // do not collect short tracks

    if ((phit->fCounter != 0)&&(tr.fNhits == 0))// store the first hit
	myphf = phit;

    if (phit->fCounter != 0){//add hit to the track
      tr.fNhits++;
      tr.fVecIHits.push_back(ihr);

      if (tr.fNhits == 2){ // calculate initial track parameters
	AliHLTTPCCAHit* phit0 = &(fVecHits[tr.fVecIHits[0]]);
	AliHLTTPCCAHit* phit1 = &(fVecHits[tr.fVecIHits[1]]);

	tr.fX = phit1->fX;
	tr.fY = phit1->fY;
	tr.fZ = phit1->fZ;

	Double_t dx = (phit1->fX-phit0->fX); // assume different x !!!
	if (TMath::Abs(dx) < 0.0001) dx = 0.0001;
	tr.fTy = (phit1->fY-phit0->fY)/dx;
	tr.fTz = (phit1->fZ-phit0->fZ)/dx;

	tr.fCovy   = phit1->fErry*phit1->fErry;
	tr.fCovty  = (phit0->fErry*phit0->fErry + phit1->fErry*phit1->fErry)/(dx*dx);
	tr.fCovyty = (phit1->fErry*phit1->fErry)/dx;

	tr.fCovz   = phit1->fErrz*phit1->fErrz;
	tr.fCovtz  = (phit0->fErrz*phit0->fErrz + phit1->fErrz*phit1->fErrz)/(dx*dx);
	tr.fCovztz = (phit1->fErrz*phit1->fErrz)/dx;
      }

      if (tr.fNhits > 2){ 
	// propagate the track
	Double_t dx = (phit->fX-tr.fX); // assume different x !!!
	if (TMath::Abs(dx) < 0.0001) dx = 0.0001;
	
	Double_t ye = tr.fY + tr.fTy*dx;
	Double_t ze = tr.fZ + tr.fTz*dx;
	
	Double_t fCovy   = tr.fCovy + 2.0*tr.fCovyty*dx + tr.fCovty*dx*dx;
	Double_t fCovyty = tr.fCovyty + tr.fCovty*dx;
	Double_t fCovty  = tr.fCovty;

	Double_t fCovz   = tr.fCovz + 2.0*tr.fCovztz*dx + tr.fCovtz*dx*dx;
	Double_t fCovztz = tr.fCovztz + tr.fCovtz*dx;
	Double_t fCovtz  = tr.fCovtz;

	
	Double_t dy = phit->fY - ye;
	Double_t dz = phit->fZ - ze;
	
	// add the measurement

	Double_t w = 1.0/(phit->fErry*phit->fErry + fCovy);

	tr.fY  = ye + fCovy*dy*w;
	tr.fTy = tr.fTy + fCovyty*dy*w;

	tr.fCovy   = fCovy - fCovy*fCovy*w;
	tr.fCovyty = fCovyty - fCovy*fCovyty*w;
	tr.fCovty  = fCovty -fCovyty*fCovyty*w;

	tr.fChi2 += dy*dy*w;
	tr.fNdf  += 1;

	w = 1.0/(phit->fErrz*phit->fErrz + fCovz);

	tr.fZ  = ze + fCovz*dz*w;
	tr.fTz = tr.fTz + fCovztz*dz*w;

	tr.fCovz   = fCovz - fCovz*fCovz*w;
	tr.fCovztz = fCovztz - fCovz*fCovztz*w;
	tr.fCovtz  = fCovtz -fCovztz*fCovztz*w;

	tr.fX += dx;

	tr.fChi2 += dz*dz*w;
	tr.fNdf  += 1;
      }
    }
    if (phit->fCounter == 1){// store the track

      tr.fGood  =  1;
      tr.fUsed  =  0;
      tr.fNext  = -1;
      tr.fPatch = patch;

      Double_t trdist = TMath::Sqrt((myphf->fX-phit->fX)*(myphf->fX-phit->fX)+
			     (myphf->fY-phit->fY)*(myphf->fY-phit->fY)+
			     (myphf->fZ-phit->fZ)*(myphf->fZ-phit->fZ));

      if ((trdist > 1.0)&&(tr.fChi2/tr.fNdf < 10.0)){
	fVecPatchTracks.push_back(tr);
	fPatchLastTrackInd++;
      }

      // reset the track
      tr.fNhits = 0;
      tr.fChi2  = 0.0;
      tr.fNdf   = -4;
      tr.fVecIHits.clear();
    }
  }
  
  // sort tracks according (currently) to number of hits
  //sort(fVecPatchTracks.begin(), fVecPatchTracks.end(), compareNhits); //don't need to sort here


 

 if(fDRAW){

  TLine *trline = new TLine();
  trline->SetLineColor(kGreen);

  for (Int_t itr = fPatchFirstTrackInd; itr <= fPatchLastTrackInd; itr++){

    AliHLTTPCCATrack *ptr = &(fVecPatchTracks[itr]);    

    AliHLTTPCCAHit* phitf = &(fVecHits[ptr->fVecIHits.back()]);
    AliHLTTPCCAHit* phitl = &(fVecHits[ptr->fVecIHits.front()]);

    fYX->cd();
    Double_t yf = ptr->fY + ptr->fTy*(phitf->fX-ptr->fX);
    Double_t yl = ptr->fY + ptr->fTy*(phitl->fX-ptr->fX);
    trline->DrawLine(yf, phitf->fX, yl, phitl->fX);

    fZX->cd();
    Double_t zf = ptr->fZ + ptr->fTz*(phitf->fX-ptr->fX);
    Double_t zl = ptr->fZ + ptr->fTz*(phitl->fX-ptr->fX);
    trline->DrawLine(zf, phitf->fX, zl, phitl->fX);

  }

  delete trline;

  fYX->cd(); fYX->Update(); fYX->Draw();
  fZX->cd(); fZX->Update(); fZX->Draw();     

}//fDRAW

  return;
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAFindSliceTracks()
{
  //!

  Int_t ntracks = fVecPatchTracks.size();

  // merge patch tracks into slice tracks
  // start with the longest tracks in the first patches

  // sort tracks according (currently) to patch and within patch - number of hits
  sort(fVecPatchTracks.begin(), fVecPatchTracks.end(), CompareCATracks);

  for (Int_t itr1 = 0; itr1 < ntracks-1; itr1++){

    AliHLTTPCCATrack *ptr1 = &(fVecPatchTracks[itr1]);    

    Double_t maxdisty = 10.0;
    Double_t maxdistz =  5.0;

    Double_t mindist = 10.0;

    Int_t next = -1;
#ifdef XXX
    for (Int_t itr2 = itr1+1; itr2 < ntracks; itr2++){

      AliHLTTPCCATrack *ptr2 = &(fVecPatchTracks[itr2]);    

      if (ptr2->patch == ptr1->patch) continue; // the same patch - no prolongation
      if (ptr2->patch >  ptr1->patch+1) break;  // sorted tracks  - no prolongation skippeng over patches

      if (ptr2->fUsed == 1) continue; // already matched

      AliHLTTPCCAHit* phitf2 = &(fVecHits[ptr2->fVecIHits.back()]);
      AliHLTTPCCAHit* phitl2 = &(fVecHits[ptr2->fVecIHits.front()]);
	
      // extrapolate tr1 to the both ends of tr2

      Double_t dyf = ptr1->fY + ptr1->fty*(phitf2->fX-ptr1->fX) - phitf2->fY; 
      Double_t dyl = ptr1->fY + ptr1->fTy*(phitl2->fX-ptr1->fX) - phitl2->fY; 

      Double_t dzf = ptr1->fZ + ptr1->fTz*(phitf2->fX-ptr1->fX) - phitf2->fZ; 
      Double_t dzl = ptr1->fZ + ptr1->fTz*(phitl2->fX-ptr1->fX) - phitl2->fZ;

      // roughly similar tracks ?
      if ((TMath::Abs(dyf) > maxdisty)||(TMath::Abs(dyl) > maxdisty)||(TMath::Abs(dzf) > maxdistz)||(TMath::Abs(dzl) > maxdistz)) continue; 

      // roughly parallel tracks ?
      Double_t dist = TMath::Sqrt((dyf-dyl)*(dyf-dyl)+(dzf-dzl)*(dzf-dzl));
      if (dist > mindist) continue; 

      mindist = dist;
      next    = itr2;

    }
#endif
    
    // found track prolongations?
    if (next != -1){
      AliHLTTPCCATrack *ptr2 = &(fVecPatchTracks[next]);    

      ptr1->fNext = next;
      ptr2->fUsed = 1;


    }

  }

  AliHLTTPCCATrack tr;
  tr.fNhits = 0; 
  tr.fChi2  = 0.0; 
  tr.fNdf   = -4; 
  tr.fVecIHits.clear();


  //collect tracks
  for (Int_t itr = 0; itr < ntracks; itr++){

    AliHLTTPCCATrack *ptr = &(fVecPatchTracks[itr]);    
    if (ptr->fUsed) continue; // start with a track not used in prolongations

    Int_t first = 1;
    do{
      if (first == 1){
	first = 0;
      }
      else{
	ptr = &(fVecPatchTracks[ptr->fNext]);    
      }
      
      for (Int_t ih = 0; ih < ptr->fNhits; ih++){
	tr.fVecIHits.push_back(ptr->fVecIHits[ih]);
	tr.fNhits++;
      }
    }while(ptr->fNext != -1);

    //sort hits according to increasing x
    sort(tr.fVecIHits.begin(), tr.fVecIHits.end(), CompareCAHitsX);

    fVecSliceTracks.push_back(tr);

    tr.fNhits = 0; 
    tr.fChi2  = 0.0; 
    tr.fNdf   = -4; 
    tr.fVecIHits.clear();
    
  }

  //if(fDRAW){

  TLine *trline = new TLine();
  trline->SetLineColor(kRed);

  //} //fDRAW

#if WRITETRACKS
  UInt_t size = 0;
  UInt_t nTracks = 0;
#endif


  for (vector<AliHLTTPCCATrack>::iterator trIt = fVecSliceTracks.begin(); trIt != fVecSliceTracks.end(); trIt++){
    
    AliHLTTPCCAHit* phit0 = &(fVecHits[trIt->fVecIHits[trIt->fNhits-1]]);
    AliHLTTPCCAHit* phit1 = &(fVecHits[trIt->fVecIHits[trIt->fNhits-2]]);
    
    trIt->fX = phit1->fX;
    trIt->fY = phit1->fY;
    trIt->fZ = phit1->fZ;
    
    Double_t dx = (phit1->fX-phit0->fX); // assume different x !!!
    if (TMath::Abs(dx) < 0.0001) dx = 0.0001;
    trIt->fTy = (phit1->fY-phit0->fY)/dx;
    trIt->fTz = (phit1->fZ-phit0->fZ)/dx;
    
    trIt->fCovy   = phit1->fErry*phit1->fErry;
    trIt->fCovty  = (phit0->fErry*phit0->fErry + phit1->fErry*phit1->fErry)/(dx*dx);
    trIt->fCovyty = (phit1->fErry*phit1->fErry)/dx;
    
    trIt->fCovz   = phit1->fErrz*phit1->fErrz;
    trIt->fCovtz  = (phit0->fErrz*phit0->fErrz + phit1->fErrz*phit1->fErrz)/(dx*dx);
    trIt->fCovztz = (phit1->fErrz*phit1->fErrz)/dx;

    for (Int_t ih = trIt->fNhits-3; ih >= 0; ih--){

      AliHLTTPCCAHit* phit = &(fVecHits[trIt->fVecIHits[ih]]);

      // propagate the track
      Double_t dx = (phit->fX-trIt->fX); // assume different x !!!
      if (TMath::Abs(dx) < 0.0001) dx = 0.0001;
      
      Double_t ye = trIt->fY + trIt->fTy*dx;
      Double_t ze = trIt->fZ + trIt->fTz*dx;
      
      Double_t fCovy   = trIt->fCovy + 2.0*trIt->fCovyty*dx + trIt->fCovty*dx*dx;
      Double_t fCovyty = trIt->fCovyty + trIt->fCovty*dx;
      Double_t fCovty  = trIt->fCovty;
      
      Double_t fCovz   = trIt->fCovz + 2.0*trIt->fCovztz*dx + trIt->fCovtz*dx*dx;
      Double_t fCovztz = trIt->fCovztz + trIt->fCovtz*dx;
      Double_t fCovtz  = trIt->fCovtz;
      
      Double_t dy = phit->fY - ye;
      Double_t dz = phit->fZ - ze;
      
      // add the measurement
      
      Double_t w = 1.0/(phit->fErry*phit->fErry + fCovy);
      
      trIt->fY  = ye + fCovy*dy*w;
      trIt->fTy = trIt->fTy + fCovyty*dy*w;
      
      trIt->fCovy   = fCovy - fCovy*fCovy*w;
      trIt->fCovyty = fCovyty - fCovy*fCovyty*w;
      trIt->fCovty  = fCovty -fCovyty*fCovyty*w;
      
      trIt->fChi2 += dy*dy*w;
      trIt->fNdf  += 1;
      
      w = 1.0/(phit->fErrz*phit->fErrz + fCovz);
      
      trIt->fZ  = ze + fCovz*dz*w;
      trIt->fTz = trIt->fTz + fCovztz*dz*w;
      
      trIt->fCovz   = fCovz - fCovz*fCovz*w;
      trIt->fCovztz = fCovztz - fCovz*fCovztz*w;
      trIt->fCovtz  = fCovtz -fCovztz*fCovztz*w;
      
      trIt->fChi2 += dz*dz*w;
      trIt->fNdf  += 1;

      trIt->fX    += dx;
    }

    trIt->fGood  =  1;
    trIt->fUsed  =  0;
    trIt->fNext  = -1;
    trIt->fPatch = -1;
    
    //if (trIt->chi2/trIt->ndf < 1.0)
    //trIt->good  =  0;


    // JMT 2006/11/13 Write Tracks to container
#if WRITETRACKS
    AliHLTTPCCAHit* firstHit = &(fVecHits[trIt->fVecIHits[0]]);
    AliHLTTPCCAHit* lastHit = &(fVecHits[trIt->fVecIHits[trIt->fNhits-1]]);
    
    Float_t xFirst = (Float_t) firstHit->fX;
    Float_t yFirst = (Float_t) trIt->fY + trIt->fTy*(xFirst-trIt->fX);
    Float_t zFirst = (Float_t) trIt->fZ + trIt->fTz*(xFirst-trIt->fX);

    Float_t xLast = (Float_t) lastHit->fX;
    Float_t yLast = (Float_t) trIt->fY + trIt->fTy*(xLast-trIt->fX);
    Float_t zLast = (Float_t) trIt->fZ + trIt->fTz*(xLast-trIt->fX);

    fOutputPtr->fX = xFirst;
    fOutputPtr->fY = yFirst;
    fOutputPtr->fZ = zFirst;
    fOutputPtr->fPt = -9876.0;
    fOutputPtr->fPterr = 1;
    fOutputPtr->fLastX = xLast;
    fOutputPtr->fLastY = yLast;
    fOutputPtr->fLastZ = zLast;    
    fOutputPtr->fPsi = atan(trIt->fTy);
    fOutputPtr->fTgl = trIt->fTz;
    fOutputPtr->fPsierr = 1;
    fOutputPtr->fTglerr = 1;
    fOutputPtr->fCharge = 1;
    fOutputPtr->fNPoints = 0;

    Byte_t *tmpP = (Byte_t *)fOutputPtr;

    tmpP += sizeof(AliHLTTPCTrackSegmentData); //+fOutputPtr->fNPoints*sizeof(UInt_t);
    size += sizeof(AliHLTTPCTrackSegmentData); //+fOutputPtr->fNPoints*sizeof(UInt_t);
    fOutputPtr = (AliHLTTPCTrackSegmentData*)tmpP;
    nTracks++;
#endif

if(fDRAW){
    AliHLTTPCCAHit* phitf = &(fVecHits[trIt->fVecIHits[0]]);
    AliHLTTPCCAHit* phitl = &(fVecHits[trIt->fVecIHits[trIt->fNhits-1]]);
    
    Double_t xf = phitf->fX;
    Double_t yf = trIt->fY + trIt->fTy*(xf-trIt->fX);
    Double_t zf = trIt->fZ + trIt->fTz*(xf-trIt->fX);

    Double_t xl = phitl->fX;
    Double_t yl = trIt->fY + trIt->fTy*(xl-trIt->fX);
    Double_t zl = trIt->fZ + trIt->fTz*(xl-trIt->fX);

    fYX->cd();
    trline->DrawLine(yf, xf, yl, xl);
    
    fZX->cd();
    trline->DrawLine(zf, xf, zl, xl);

} //fDRAW
      
  }

#if WRITETRACKS
  fOutputNTracks = nTracks;
  fOutputSize = size;
  cout << "NTRACKS=" << nTracks << endl;
  cout << "SIZEoF=" <<  sizeof(AliHLTTPCTrackSegmentData) << endl;
  cout << "SIZE=" << fOutputSize << endl;
#endif

if(fDRAW){
    
  delete trline;

  fYX->cd(); fYX->Update(); fYX->Draw();
  fZX->cd(); fZX->Update(); fZX->Draw();     

 }//fDRAW
  
  return;
}
