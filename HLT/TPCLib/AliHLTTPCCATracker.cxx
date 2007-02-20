/** \class AliHLTTPCCATracker
<pre>
//_____________________________________________________________
// AliHLTTPCCATracker
//
//
// Author: Ivan Kisel
// Copyright &copy ALICE HLT Group
</pre>
*/

#define WRITETRACKS 1

#include <sys/time.h>


#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCLogging.h" 
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCConfMapTrack.h"
#include "AliHLTTPCConfMapPoint.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCCATracker.h"
#include "AliHLTTPCTrackSegmentData.h"

#include "TApplication.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TEllipse.h"
#include "TText.h"

#include <iostream>

#if __GNUC__ >= 3
using namespace std;
#endif

static TList *listHisto;
static TH1F *h_NClusters;
static TH2F *h_ClustersXY[10];

static Int_t h_Event;

//#define DRAW

#ifdef DRAW

static TApplication *myapp;
static TCanvas *YX, *ZX;

static bool ask = false;

#endif //DRAW

ClassImp(AliHLTTPCCATracker)

// ----------------------------------------------------------------------------------
AliHLTTPCCATracker::AliHLTTPCCATracker()
{
  //Default constructor
  fVertex = NULL;
  fTrack = NULL;
  fHit = NULL;
  fVolume = NULL;
  fRow = NULL;
  fBench = (Bool_t)true;
  fVertexConstraint = (Bool_t)true;
  fParamSet[0]=0;
  fParamSet[1]=0;
  //JMT 2006/11/13
  fOutputPtr= NULL;
  fOutputNTracks=0;
  fOutputSize=0;
}

// ----------------------------------------------------------------------------------
AliHLTTPCCATracker::~AliHLTTPCCATracker()
{
  // Destructor.
  if(fVolume) {
    delete [] fVolume;
  }
  if(fRow) {
    delete [] fRow;
  }
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
    static bool first_call = true;
    if (first_call){
      first_call = false;

      h_Event = 0;

      TDirectory *curdir = gDirectory;
      TDirectory *histodir = gROOT->mkdir("histodir");
      histodir->cd();
      
      h_NClusters  = new TH1F("h_NClusters", "Number of clusters in event", 100, 0.0, 1000.0);

      h_ClustersXY[0] = new TH2F("h_ClustersXY_0", "Clusters in XY plane Ev 0", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[1] = new TH2F("h_ClustersXY_1", "Clusters in XY plane Ev 1", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[2] = new TH2F("h_ClustersXY_2", "Clusters in XY plane Ev 2", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[3] = new TH2F("h_ClustersXY_3", "Clusters in XY plane Ev 3", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[4] = new TH2F("h_ClustersXY_4", "Clusters in XY plane Ev 4", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[5] = new TH2F("h_ClustersXY_5", "Clusters in XY plane Ev 5", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[6] = new TH2F("h_ClustersXY_6", "Clusters in XY plane Ev 6", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[7] = new TH2F("h_ClustersXY_7", "Clusters in XY plane Ev 7", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[8] = new TH2F("h_ClustersXY_8", "Clusters in XY plane Ev 8", 100, 50.0, 300.0, 100, -100.0, 100.0);
      h_ClustersXY[9] = new TH2F("h_ClustersXY_9", "Clusters in XY plane Ev 9", 100, 50.0, 300.0, 100, -100.0, 100.0);

      listHisto = gDirectory->GetList();
      curdir->cd();

#ifdef DRAW

      myapp = new TApplication("myapp", 0, 0);
      
      gStyle->SetCanvasBorderMode(0);
      gStyle->SetCanvasBorderSize(1);
      gStyle->SetCanvasColor(0);


      YX = new TCanvas ("YX", "YX (Pad-Row) window", -10, 0, 680, 400);
      YX->Range(-50.0, 80.0, 50.0, 250.0);
      YX->Clear();
      YX->Draw();

      ZX = new TCanvas ("ZX", "ZX window", -700, 0, 680, 400);
      ZX->Range(-250.0, 80.0, 250.0, 250.0);
      ZX->Clear();
      ZX->Draw();

#ifdef XXX
      char symbol;
      if (ask){
        do{
          std::cin.get(symbol);
          if (symbol == 'r')
            ask = false;
        } while (symbol != '\n');
      }
#endif //XXX

#endif //DRAW

    } 
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAWriteHistos()
{
    // IK Open the output file and write histogramms

    TFile *outfile = new TFile("ConfMapper_histo.root", "RECREATE");
    TIter hiter(listHisto);
    while (TObject *obj = hiter()) obj->Write();
    outfile->Close();
}

// ----------------------------------------------------------------------------------
Bool_t AliHLTTPCCATracker::ReadHits(UInt_t count, AliHLTTPCSpacePointData* hits )
{

  //read hits
  Int_t nhit=(Int_t)count; 

  vector<CAHit> vec_hits; 
  vec_hits.clear();

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
  vec_hits.clear();
  vec_patch_tracks.clear();
  vec_slice_tracks.clear();

  return;
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAReadPatchHits(Int_t patch, UInt_t count, AliHLTTPCSpacePointData* hits )
{
  //read hits
  Int_t nhit=(Int_t)count; 

  patch_first_hit_ind = vec_hits.size();
  patch_last_hit_ind  = patch_first_hit_ind + nhit - 1;

  for (Int_t i=0;i<nhit;i++)
    {	
      AliHLTTPCSpacePointData* pSP = &(hits[i]);

      CAHit hit;
      hit.x = pSP->fX;
      hit.y = pSP->fY;
      hit.z = pSP->fZ;

      hit.errx = 1.0;//pSP->fSigmaY2;
      if (patch <= 3)
	hit.erry = 10.0;//1.0;//pSP->fSigmaY2;
      else
	hit.erry = 10.0;//pSP->fSigmaY2;
      hit.errz = 1.0;//pSP->fSigmaZ2;

      hit.index = i;
      hit.counter = 0;

      vec_hits.push_back(hit);
    }

#ifdef DRAW

  TMarker *mhit = new TMarker(0.0, 0.0, 6);
  mhit->SetMarkerColor(kBlue);

  for (Int_t i = patch_first_hit_ind; i <= patch_last_hit_ind; i++)
    {
      YX->cd();
      mhit->DrawMarker(vec_hits[i].y, vec_hits[i].x);
      
      ZX->cd();
      mhit->DrawMarker(vec_hits[i].z, vec_hits[i].x);
    }

  delete mhit;

  YX->cd(); YX->Update(); YX->Draw(); 
  ZX->cd(); ZX->Update(); ZX->Draw();     

#endif //DRAW

  return;
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAFindPatchTracks(Int_t patch)
{

#ifdef DRAW

  TLine *line = new TLine();
  line->SetLineColor(kBlack);

#endif //DRAW

  Double_t dist01, dist12, dist02, mindist, chi2, minchi2;
  CAHit *ph0, *ph1, *ph2, *phbest, *next;

  Double_t x, y, z, ty, tz, cov_y, cov_ty, cov_yty, cov_z, cov_tz, cov_ztz, ye, ze, dx, dy, dz, sy2, sz2; 

  for (Int_t i1 = patch_first_hit_ind; i1 < patch_last_hit_ind; i1++) // <nhit-1
    {	
      ph1  = &(vec_hits[i1]);
      next = &(vec_hits[i1+1]);

      mindist = 20.0;
      minchi2 = 1.0;
      phbest = NULL;

      if (ph1->counter != 0){ // the previous hit exists
	x = ph1->x;
	y = ph1->y;
	z = ph1->z;

	dx = (ph1->x-ph0->x); // assume different x !!!
	if (abs(dx) < 0.0001) dx = 0.0001;
	ty = (ph1->y-ph0->y)/dx;
	tz = (ph1->z-ph0->z)/dx;

	cov_y   = ph1->erry*ph1->erry;
	cov_ty  = (ph0->erry*ph0->erry + ph1->erry*ph1->erry)/(dx*dx);
	cov_yty = (ph1->erry*ph1->erry)/dx;

	cov_z   = ph1->errz*ph1->errz;
	cov_tz  = (ph0->errz*ph0->errz + ph1->errz*ph1->errz)/(dx*dx);
	cov_ztz = (ph1->errz*ph1->errz)/dx;

      }

      for (Int_t i2 = i1+1; i2 <= patch_last_hit_ind; i2++)
	{	
	  ph2 = &(vec_hits[i2]);
	  if (ph2->x < ph1->x) continue; // no backward direction!

	  dist12 = sqrt((ph1->x-ph2->x)*(ph1->x-ph2->x) + (ph1->y-ph2->y)*(ph1->y-ph2->y) + (ph1->z-ph2->z)*(ph1->z-ph2->z));
	  if (dist12 < mindist){
	    if (ph1->counter == 0){
	      mindist = dist12;
	      phbest  = ph2;
	    }
	    else{
	      // calculate chi2 distance between the hit and the line
	      dx = (ph2->x-x); // assume different x !!!
	      if (abs(dx) < 0.0001) dx = 0.0001;
	      
	      ye = y + ty*dx;
	      ze = z + tz*dx;

	      sy2 = cov_y + 2.0*cov_yty*dx + cov_ty*dx*dx;
	      sz2 = cov_z + 2.0*cov_ztz*dx + cov_tz*dx*dx;

	      dy = ph2->y - ye;
	      dz = ph2->z - ze;
	      
	      chi2 = (dy*dy)/(sy2+ph2->erry*ph2->erry) + (dz*dz)/(sz2+ph2->errz*ph2->errz);

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
	Double_t x_tmp = next->x;
	Double_t y_tmp = next->y;
	Double_t z_tmp = next->z;

	Double_t errx_tmp = next->errx;
	Double_t erry_tmp = next->erry;
	Double_t errz_tmp = next->errz;

	Int_t index_tmp = next->index;
	Int_t counter_tmp = next->counter;

	next->x = phbest->x;
	next->y = phbest->y;
	next->z = phbest->z;

	next->errx = phbest->errx;
	next->erry = phbest->erry;
	next->errz = phbest->errz;

	next->index = phbest->index;
	next->counter = phbest->counter;

	phbest->x = x_tmp;
	phbest->y = y_tmp;
	phbest->z = z_tmp;

	phbest->errx = errx_tmp;
	phbest->erry = erry_tmp;
	phbest->errz = errz_tmp;

	phbest->index = index_tmp;
	phbest->counter = counter_tmp;

	if (ph1->counter == 0)
	  ph1->counter = 1;
	next->counter = ph1->counter+1;

	ph0 = ph1;
	dist01 = mindist;

#ifdef DRAW

	YX->cd();
	line->DrawLine(ph1->y, ph1->x, next->y, next->x);
	
	ZX->cd();
	line->DrawLine(ph1->z, ph1->x, next->z, next->x);

#endif //DRAW

      }

    }

#ifdef DRAW

  delete line;

  YX->cd(); YX->Update(); YX->Draw(); 
  ZX->cd(); ZX->Update(); ZX->Draw();     

#endif //DRAW
  
  CATrack tr;
  tr.nhits = 0; 
  tr.chi2  = 0.0; 
  tr.ndf   = -4; 
  tr.vec_ihits.clear();

  patch_first_track_ind = vec_patch_tracks.size(); 
  patch_last_track_ind  = patch_first_track_ind - 1;

  CAHit *my_phf; 

  // collect tracks
  for (Int_t ihr = patch_last_hit_ind; ihr >= patch_first_hit_ind; ihr--){

    CAHit* phit = &(vec_hits[ihr]);

    if ((phit->counter < 5)&&(tr.nhits == 0))
      continue; // do not collect short tracks

    if ((phit->counter != 0)&&(tr.nhits == 0))// store the first hit
	my_phf = phit;

    if (phit->counter != 0){//add hit to the track
      tr.nhits++;
      tr.vec_ihits.push_back(ihr);

      if (tr.nhits == 2){ // calculate initial track parameters
	CAHit* phit0 = &(vec_hits[tr.vec_ihits[0]]);
	CAHit* phit1 = &(vec_hits[tr.vec_ihits[1]]);

	tr.x = phit1->x;
	tr.y = phit1->y;
	tr.z = phit1->z;

	Double_t dx = (phit1->x-phit0->x); // assume different x !!!
	if (abs(dx) < 0.0001) dx = 0.0001;
	tr.ty = (phit1->y-phit0->y)/dx;
	tr.tz = (phit1->z-phit0->z)/dx;

	tr.cov_y   = phit1->erry*phit1->erry;
	tr.cov_ty  = (phit0->erry*phit0->erry + phit1->erry*phit1->erry)/(dx*dx);
	tr.cov_yty = (phit1->erry*phit1->erry)/dx;

	tr.cov_z   = phit1->errz*phit1->errz;
	tr.cov_tz  = (phit0->errz*phit0->errz + phit1->errz*phit1->errz)/(dx*dx);
	tr.cov_ztz = (phit1->errz*phit1->errz)/dx;
      }

      if (tr.nhits > 2){ 
	// propagate the track
	Double_t dx = (phit->x-tr.x); // assume different x !!!
	if (abs(dx) < 0.0001) dx = 0.0001;
	
	Double_t ye = tr.y + tr.ty*dx;
	Double_t ze = tr.z + tr.tz*dx;
	
	Double_t cov_y   = tr.cov_y + 2.0*tr.cov_yty*dx + tr.cov_ty*dx*dx;
	Double_t cov_yty = tr.cov_yty + tr.cov_ty*dx;
	Double_t cov_ty  = tr.cov_ty;

	Double_t cov_z   = tr.cov_z + 2.0*tr.cov_ztz*dx + tr.cov_tz*dx*dx;
	Double_t cov_ztz = tr.cov_ztz + tr.cov_tz*dx;
	Double_t cov_tz  = tr.cov_tz;

	
	Double_t dy = phit->y - ye;
	Double_t dz = phit->z - ze;
	
	// add the measurement

	Double_t w = 1.0/(phit->erry*phit->erry + cov_y);

	tr.y  = ye + cov_y*dy*w;
	tr.ty = tr.ty + cov_yty*dy*w;

	tr.cov_y   = cov_y - cov_y*cov_y*w;
	tr.cov_yty = cov_yty - cov_y*cov_yty*w;
	tr.cov_ty  = cov_ty -cov_yty*cov_yty*w;

	tr.chi2 += dy*dy*w;
	tr.ndf  += 1;

	w = 1.0/(phit->errz*phit->errz + cov_z);

	tr.z  = ze + cov_z*dz*w;
	tr.tz = tr.tz + cov_ztz*dz*w;

	tr.cov_z   = cov_z - cov_z*cov_z*w;
	tr.cov_ztz = cov_ztz - cov_z*cov_ztz*w;
	tr.cov_tz  = cov_tz -cov_ztz*cov_ztz*w;

	tr.x += dx;

	tr.chi2 += dz*dz*w;
	tr.ndf  += 1;
      }
    }
    if (phit->counter == 1){// store the track

      tr.good  =  1;
      tr.used  =  0;
      tr.next  = -1;
      tr.patch = patch;

      Double_t trdist = sqrt((my_phf->x-phit->x)*(my_phf->x-phit->x)+
			     (my_phf->y-phit->y)*(my_phf->y-phit->y)+
			     (my_phf->z-phit->z)*(my_phf->z-phit->z));

      if ((trdist > 1.0)&&(tr.chi2/tr.ndf < 10.0)){
	vec_patch_tracks.push_back(tr);
	patch_last_track_ind++;
      }

      // reset the track
      tr.nhits = 0;
      tr.chi2  = 0.0;
      tr.ndf   = -4;
      tr.vec_ihits.clear();
    }
  }
  
  // sort tracks according (currently) to number of hits
  //sort(vec_patch_tracks.begin(), vec_patch_tracks.end(), compareNhits); //don't need to sort here


#ifdef DRAW

  TLine *trline = new TLine();
  trline->SetLineColor(kGreen);

  for (Int_t itr = patch_first_track_ind; itr <= patch_last_track_ind; itr++){

    CATrack *ptr = &(vec_patch_tracks[itr]);    

    CAHit* phitf = &(vec_hits[ptr->vec_ihits.back()]);
    CAHit* phitl = &(vec_hits[ptr->vec_ihits.front()]);

    YX->cd();
    Double_t yf = ptr->y + ptr->ty*(phitf->x-ptr->x);
    Double_t yl = ptr->y + ptr->ty*(phitl->x-ptr->x);
    trline->DrawLine(yf, phitf->x, yl, phitl->x);

    ZX->cd();
    Double_t zf = ptr->z + ptr->tz*(phitf->x-ptr->x);
    Double_t zl = ptr->z + ptr->tz*(phitl->x-ptr->x);
    trline->DrawLine(zf, phitf->x, zl, phitl->x);

  }

  delete trline;

  YX->cd(); YX->Update(); YX->Draw();
  ZX->cd(); ZX->Update(); ZX->Draw();     

#endif //DRAW

  return;
}

// ----------------------------------------------------------------------------------
void AliHLTTPCCATracker::CAFindSliceTracks()
{

  Int_t ntracks = vec_patch_tracks.size();

  // merge patch tracks into slice tracks
  // start with the longest tracks in the first patches

  // sort tracks according (currently) to patch and within patch - number of hits
  sort(vec_patch_tracks.begin(), vec_patch_tracks.end(), compareCATracks);

  for (Int_t itr1 = 0; itr1 < ntracks-1; itr1++){

    CATrack *ptr1 = &(vec_patch_tracks[itr1]);    

    Double_t maxdisty = 10.0;
    Double_t maxdistz =  5.0;

    Double_t mindist = 10.0;

    Int_t next = -1;
#ifdef XXX
    for (Int_t itr2 = itr1+1; itr2 < ntracks; itr2++){

      CATrack *ptr2 = &(vec_patch_tracks[itr2]);    

      if (ptr2->patch == ptr1->patch) continue; // the same patch - no prolongation
      if (ptr2->patch >  ptr1->patch+1) break;  // sorted tracks  - no prolongation skippeng over patches

      if (ptr2->used == 1) continue; // already matched

      CAHit* phitf2 = &(vec_hits[ptr2->vec_ihits.back()]);
      CAHit* phitl2 = &(vec_hits[ptr2->vec_ihits.front()]);
	
      // extrapolate tr1 to the both ends of tr2

      Double_t dyf = ptr1->y + ptr1->ty*(phitf2->x-ptr1->x) - phitf2->y; 
      Double_t dyl = ptr1->y + ptr1->ty*(phitl2->x-ptr1->x) - phitl2->y; 

      Double_t dzf = ptr1->z + ptr1->tz*(phitf2->x-ptr1->x) - phitf2->z; 
      Double_t dzl = ptr1->z + ptr1->tz*(phitl2->x-ptr1->x) - phitl2->z;

      // roughly similar tracks ?
      if ((abs(dyf) > maxdisty)||(abs(dyl) > maxdisty)||(abs(dzf) > maxdistz)||(abs(dzl) > maxdistz)) continue; 

      // roughly parallel tracks ?
      Double_t dist = sqrt((dyf-dyl)*(dyf-dyl)+(dzf-dzl)*(dzf-dzl));
      if (dist > mindist) continue; 

      mindist = dist;
      next    = itr2;

    }
#endif
    
    // found track prolongations?
    if (next != -1){
      CATrack *ptr2 = &(vec_patch_tracks[next]);    

      ptr1->next = next;
      ptr2->used = 1;


    }

  }

  CATrack tr;
  tr.nhits = 0; 
  tr.chi2  = 0.0; 
  tr.ndf   = -4; 
  tr.vec_ihits.clear();


  //collect tracks
  for (Int_t itr = 0; itr < ntracks; itr++){

    CATrack *ptr = &(vec_patch_tracks[itr]);    
    if (ptr->used) continue; // start with a track not used in prolongations

    Int_t first = 1;
    do{
      if (first == 1){
	first = 0;
      }
      else{
	ptr = &(vec_patch_tracks[ptr->next]);    
      }
      
      for (Int_t ih = 0; ih < ptr->nhits; ih++){
	tr.vec_ihits.push_back(ptr->vec_ihits[ih]);
	tr.nhits++;
      }
    }while(ptr->next != -1);

    //sort hits according to increasing x
    sort(tr.vec_ihits.begin(), tr.vec_ihits.end(), compareCAHitsX);

    vec_slice_tracks.push_back(tr);

    tr.nhits = 0; 
    tr.chi2  = 0.0; 
    tr.ndf   = -4; 
    tr.vec_ihits.clear();
    
  }

#ifdef DRAW

  TLine *trline = new TLine();
  trline->SetLineColor(kRed);

#endif //DRAW

#if WRITETRACKS
  UInt_t size = 0;
  UInt_t nTracks = 0;
#endif


  for (vector<CATrack>::iterator trIt = vec_slice_tracks.begin(); trIt != vec_slice_tracks.end(); trIt++){
    
    CAHit* phit0 = &(vec_hits[trIt->vec_ihits[trIt->nhits-1]]);
    CAHit* phit1 = &(vec_hits[trIt->vec_ihits[trIt->nhits-2]]);
    
    trIt->x = phit1->x;
    trIt->y = phit1->y;
    trIt->z = phit1->z;
    
    Double_t dx = (phit1->x-phit0->x); // assume different x !!!
    if (abs(dx) < 0.0001) dx = 0.0001;
    trIt->ty = (phit1->y-phit0->y)/dx;
    trIt->tz = (phit1->z-phit0->z)/dx;
    
    trIt->cov_y   = phit1->erry*phit1->erry;
    trIt->cov_ty  = (phit0->erry*phit0->erry + phit1->erry*phit1->erry)/(dx*dx);
    trIt->cov_yty = (phit1->erry*phit1->erry)/dx;
    
    trIt->cov_z   = phit1->errz*phit1->errz;
    trIt->cov_tz  = (phit0->errz*phit0->errz + phit1->errz*phit1->errz)/(dx*dx);
    trIt->cov_ztz = (phit1->errz*phit1->errz)/dx;

    for (Int_t ih = trIt->nhits-3; ih >= 0; ih--){

      CAHit* phit = &(vec_hits[trIt->vec_ihits[ih]]);

      // propagate the track
      Double_t dx = (phit->x-trIt->x); // assume different x !!!
      if (abs(dx) < 0.0001) dx = 0.0001;
      
      Double_t ye = trIt->y + trIt->ty*dx;
      Double_t ze = trIt->z + trIt->tz*dx;
      
      Double_t cov_y   = trIt->cov_y + 2.0*trIt->cov_yty*dx + trIt->cov_ty*dx*dx;
      Double_t cov_yty = trIt->cov_yty + trIt->cov_ty*dx;
      Double_t cov_ty  = trIt->cov_ty;
      
      Double_t cov_z   = trIt->cov_z + 2.0*trIt->cov_ztz*dx + trIt->cov_tz*dx*dx;
      Double_t cov_ztz = trIt->cov_ztz + trIt->cov_tz*dx;
      Double_t cov_tz  = trIt->cov_tz;
      
      Double_t dy = phit->y - ye;
      Double_t dz = phit->z - ze;
      
      // add the measurement
      
      Double_t w = 1.0/(phit->erry*phit->erry + cov_y);
      
      trIt->y  = ye + cov_y*dy*w;
      trIt->ty = trIt->ty + cov_yty*dy*w;
      
      trIt->cov_y   = cov_y - cov_y*cov_y*w;
      trIt->cov_yty = cov_yty - cov_y*cov_yty*w;
      trIt->cov_ty  = cov_ty -cov_yty*cov_yty*w;
      
      trIt->chi2 += dy*dy*w;
      trIt->ndf  += 1;
      
      w = 1.0/(phit->errz*phit->errz + cov_z);
      
      trIt->z  = ze + cov_z*dz*w;
      trIt->tz = trIt->tz + cov_ztz*dz*w;
      
      trIt->cov_z   = cov_z - cov_z*cov_z*w;
      trIt->cov_ztz = cov_ztz - cov_z*cov_ztz*w;
      trIt->cov_tz  = cov_tz -cov_ztz*cov_ztz*w;
      
      trIt->chi2 += dz*dz*w;
      trIt->ndf  += 1;

      trIt->x    += dx;
    }

    trIt->good  =  1;
    trIt->used  =  0;
    trIt->next  = -1;
    trIt->patch = -1;
    
    //if (trIt->chi2/trIt->ndf < 1.0)
    //trIt->good  =  0;


    // JMT 2006/11/13 Write Tracks to container
#if WRITETRACKS
    CAHit* firstHit = &(vec_hits[trIt->vec_ihits[0]]);
    CAHit* lastHit = &(vec_hits[trIt->vec_ihits[trIt->nhits-1]]);
    
    Float_t xFirst = (Float_t) firstHit->x;
    Float_t yFirst = (Float_t) trIt->y + trIt->ty*(xFirst-trIt->x);
    float_t zFirst = (Float_t) trIt->z + trIt->tz*(xFirst-trIt->x);

    Float_t xLast = (Float_t) lastHit->x;
    Float_t yLast = (Float_t) trIt->y + trIt->ty*(xLast-trIt->x);
    Float_t zLast = (Float_t) trIt->z + trIt->tz*(xLast-trIt->x);

    fOutputPtr->fX = xFirst;
    fOutputPtr->fY = yFirst;
    fOutputPtr->fZ = zFirst;
    fOutputPtr->fPt = -9876.0;
    fOutputPtr->fPterr = 1;
    fOutputPtr->fLastX = xLast;
    fOutputPtr->fLastY = yLast;
    fOutputPtr->fLastZ = zLast;    
    fOutputPtr->fPsi = atan(trIt->ty);
    fOutputPtr->fTgl = trIt->tz;
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

#ifdef DRAW
    CAHit* phitf = &(vec_hits[trIt->vec_ihits[0]]);
    CAHit* phitl = &(vec_hits[trIt->vec_ihits[trIt->nhits-1]]);
    
    Double_t xf = phitf->x;
    Double_t yf = trIt->y + trIt->ty*(xf-trIt->x);
    Double_t zf = trIt->z + trIt->tz*(xf-trIt->x);

    Double_t xl = phitl->x;
    Double_t yl = trIt->y + trIt->ty*(xl-trIt->x);
    Double_t zl = trIt->z + trIt->tz*(xl-trIt->x);

    YX->cd();
    trline->DrawLine(yf, xf, yl, xl);
    
    ZX->cd();
    trline->DrawLine(zf, xf, zl, xl);

#endif //DRAW
      
  }

#if WRITETRACKS
  fOutputNTracks = nTracks;
  fOutputSize = size;
  cout << "NTRACKS=" << nTracks << endl;
  cout << "SIZEoF=" <<  sizeof(AliHLTTPCTrackSegmentData) << endl;
  cout << "SIZE=" << fOutputSize << endl;
#endif

#ifdef DRAW
    
  delete trline;

  YX->cd(); YX->Update(); YX->Draw();
  ZX->cd(); ZX->Update(); ZX->Draw();     

#endif //DRAW
  
  return;
}
