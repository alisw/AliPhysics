#define PointMaker_cxx
// The class definition in esdTree.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//

// To do: rewrite as AliAnalysisTaskSE.(like initiated by Alex at some point)

#include <TROOT.h>
#include <TStyle.h>
#include <TLinearFitter.h>
#include "AliAlignObjParams.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
#include "PointMaker.h"
//=============================================================================
PointMaker::PointMaker(char *outfil) :
  TSelector(),
  fChain(0),
  fESD(0),
  //  fESDfriend(0),
  fFile(0),
  fTree(0),
  fArray(0),
  fNevents(0),
  fNtracks(0),
  fNAcceptedTracks(0),
  fOutfil(outfil)
{
  // Constructor. Initialization of pointers
}
//=============================================================================
PointMaker::~PointMaker() {
  // Remove all pointers

  //  delete fESD;

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  delete fFile;
}
//=============================================================================
void PointMaker::Begin(TTree *)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  // create output file and tree  - TEMPORARY -until root trees are merged in the memory
  fFile = TFile::Open(fOutfil.Data(), "RECREATE",0);
  fTree = new TTree("spTree", "Tree with track space point arrays");
  fTree->Branch("SP","AliTrackPointArray", &fArray);
  printf("Begin called\n");
}
//=============================================================================
void PointMaker::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  Init(tree);

  AliAlignObjParams alobj;  // initialize align obj.  
  TString option = GetOption();

  // histograms to monitor cut efficiency and module population
  fCuttra = new TH1D("cuttra","cuttra",20,0.5,20.5); 
  fCutpoi = new TH1D("cutpoi","cutpoi",20,0.5,20.5);
  fModpop = new TH2D("modpop","modpop",90,-0.5,89.5,30,-0.5,29.5);
  fModpop->SetXTitle("module nr");
  fModpop->SetYTitle("layer nr");
  printf("SlaveBegin called\n");
}
//=============================================================================
void PointMaker::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  // Set branch addresses
  if (tree == 0) return;
  fChain = tree;
  
  /*
  fChain->SetBranchAddress("ESD",&fESD);
  fChain->SetBranchStatus("ESDfriend*",1);
  fChain->SetBranchAddress("ESDfriend.",&fESDfriend);
  */
  
  fESD = new AliESDEvent();
  fChain->SetBranchStatus("ESDfriend*",1);
  fESD->ReadFromTree(fChain);
  fESDfriend = (AliESDfriend*)fESD->FindListObject("AliESDfriend");
  if(!fESDfriend) fChain->SetBranchAddress("ESDfriend.",&fESDfriend); 
  printf("Init called\n");
}
//=============================================================================
Bool_t PointMaker::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.

  printf("Notify called\n");
  return kTRUE;
}
//=============================================================================
Bool_t PointMaker::IsIdenticalWithOneOf(AliTrackPoint *p, AliTrackPointArray *parray, int nmax) {

  // Is the point p identical with one of the points on the list parray?
  // This is a fix for aliroot 4-16-Rev-01 (and before) writing some 
  // spurious uninitialized points. 

  for (int i=0; i<parray->GetNPoints() && i<nmax; i++) {
    AliTrackPoint pa;
    parray->GetPoint(pa,i);
    //printf("comparing %7.3f with %7.3f\n",p->GetY(),pa.GetY());
    if (p->GetResidual(pa,0)<1e-8) return kTRUE;
    //printf("different\n");
  }
  return kFALSE;
}
//=============================================================================
Bool_t PointMaker::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fChain is the pointer to the TChain being processed,
  //  use fChain->GetTree()->GetEntry(entry).

  fChain->GetTree()->GetEntry(entry);
  fESD->SetESDfriend(fESDfriend);

  TLinearFitter fitter(2, "pol1");
  TLinearFitter fitterz(2, "pol1");

  // track cuts

  int tpc = 0; // require tpc
  int ptu = 0; // require certain pt's (magnetic field and tpc presumably on)

  const Float_t kMaxDelta       = 1;
  const Float_t kMinNcl         = 60;
  const Float_t kMinPtLow       = 0.2;
  const Float_t kMinNclLow      = 100;
  const Float_t kMinPt0         = 2;
  UInt_t status = AliESDtrack::kTRDrefit; 
  if (tpc) status |= AliESDtrack::kTPCrefit; 

  const Float_t kMinRadius2  = 2*2;
  const Float_t kMaxRadius2  = 400*400;
  const Float_t kDeadSpace   = 4;
  const Float_t kTan = TMath::Tan(10*TMath::DegToRad());
  Int_t ntracks = fESD->GetNumberOfTracks();
  const AliTrackPointArray *array = 0;
  AliTrackPointArray *tmpArray = 0;
  // trdarray contains all trd points in this event, for duplication detection
  AliTrackPointArray *trdarray = new AliTrackPointArray(1000);
  int ntrdarray = 0;

  int pr = (!gROOT->IsBatch());
  for (Int_t itrack=0; itrack < ntracks; itrack++) {
    if (pr) printf("\revent %d   track %d",fNevents,itrack);
    fCuttra->Fill(1);
    AliESDtrack * track = fESD->GetTrack(itrack);
    fNtracks++;
    if (!track) continue;
    fCuttra->Fill(2);
    //if ((track->GetStatus() & status) == 0) continue;
    fCuttra->Fill(3);
    if (track->GetKinkIndex(0)!=0) continue;
    fCuttra->Fill(4);
    if (tpc) if (track->GetTPCNcls()<kMinNcl) continue;
    fCuttra->Fill(5);
    if (ptu) if (track->GetP() < kMinPtLow) continue;
    fCuttra->Fill(6);
    if (ptu) if (track->GetP() < kMinPt0 && track->GetTPCNcls()<kMinNclLow) continue;
    fCuttra->Fill(7);
    //
    // select points
    //
    array = track->GetTrackPointArray();    
    if (!array) continue;
    Int_t npoints = array->GetNPoints();
    if (tmpArray) delete tmpArray;
    tmpArray = new AliTrackPointArray(npoints);
    Int_t current = 0;
    int ntpc = 0; // number of good TPC points
    int ntrd = 0; // number of good TRD points
    for (Int_t ipoint=0; ipoint<npoints; ipoint++){
      fCutpoi->Fill(1);
      AliTrackPoint p;
      array->GetPoint(p, ipoint);

      // filter buggy points

      UShort_t volid = array->GetVolumeID()[ipoint];
      Int_t iModule;
      AliGeomManager::ELayerID layer = AliGeomManager::VolUIDToLayer(volid,iModule);
      if ((layer < AliGeomManager::kFirstLayer) || (layer >= AliGeomManager::kLastLayer)) continue;
      fCutpoi->Fill(2);
      if ((iModule >= AliGeomManager::LayerSize(layer)) || (iModule < 0)) continue;
      fCutpoi->Fill(3);
      Float_t r2 = p.GetX()*p.GetX()+p.GetY()*p.GetY();
      if ( r2<kMinRadius2 || r2 > kMaxRadius2 ) continue;
      fCutpoi->Fill(4);

      // TPC clusters handling

      if (layer>=AliGeomManager::kTPC1 && layer<=AliGeomManager::kTPC2){
	// invalid covariance 
	if (p.GetCov()[0]<0 || p.GetCov()[3]<0 ||  p.GetCov()[5]<0) continue;
	fCutpoi->Fill(7);

	// remove edge clusters

	AliTrackPoint& plocal = p.MasterToLocal();
	Double_t ylocal  = plocal.GetY();
	Double_t zlocal  = plocal.GetZ();
	Double_t xlocal  = plocal.GetX();
	Float_t edgey = TMath::Abs(plocal.GetX()*kTan);
	Int_t nclose=0;
	fitter.ClearPoints();
	fitterz.ClearPoints();
	for (Int_t jpoint=ipoint-7; jpoint<=ipoint+7; jpoint++){
	  if (jpoint<0 || jpoint>=npoints) continue;
	  if (ipoint==jpoint) continue;
	  UShort_t volidL = array->GetVolumeID()[jpoint];
	  if (volidL!=volid) continue;
	  AliTrackPoint pc;	
	  array->GetPoint(pc, jpoint);
	  AliTrackPoint &pcl=  pc.MasterToLocal();		  
	  Double_t dx = pcl.GetX()-xlocal;
	  fitter.AddPoint(&dx,pcl.GetY(),1);	  
	  fitterz.AddPoint(&dx,pcl.GetZ(),1);	  
	  nclose++;
	}
	if (nclose<6) continue;
	fCutpoi->Fill(8);
	fitter.Eval();
	fitterz.Eval();
	Double_t fity =fitter.GetParameter(0); 
	Double_t fitz =fitterz.GetParameter(0); 
	if (TMath::Abs(ylocal-fity)>kMaxDelta) continue;  //outlier
	fCutpoi->Fill(9);
	if (TMath::Abs(zlocal-fitz)>kMaxDelta) continue;  //outlier
	fCutpoi->Fill(10);
	if (TMath::Abs(fity)>edgey-kDeadSpace) continue;   // remove edge clusters
	fCutpoi->Fill(11);
	ntpc++;
      }
      fCutpoi->Fill(12);

      // TRD track points handling

      if (layer>=AliGeomManager::kTRD1 && layer<=AliGeomManager::kTRD6){
	fCutpoi->Fill(14);
	if (IsIdenticalWithOneOf(&p,trdarray,ntrdarray)) continue; // bug fix
	trdarray->AddPoint(ntrdarray++,&p);
	fCutpoi->Fill(15);
	fModpop->Fill(iModule,layer);
	ntrd++;
      }
      //printf("event %4d  track %4d  volid %4d  layer %4d  module %4d  %10.3f %10.3f %10.3f\n",
      //     (int)entry,itrack,volid,layer,iModule,p.GetX(),p.GetY(),atan2(p.GetY(),p.GetX())*180/3.1416);
      tmpArray->AddPoint(current,&p);
      current++;
      fCutpoi->Fill(16);
    }
    fCuttra->Fill(8);
    if (ntpc < 100) continue;
    fCuttra->Fill(9);
    if (ntrd < 4) continue;
    fCuttra->Fill(10);
    if (fArray) delete fArray;
    fArray = new AliTrackPointArray(current);
    for (Int_t ipoint=0; ipoint<current; ipoint++){
      AliTrackPoint p;
      tmpArray->GetPoint(p, ipoint);
      fArray->AddPoint(ipoint,&p);
    }   
    fNAcceptedTracks++;
    fTree->Fill();
  }
  delete trdarray;
  fNevents++;
  return kTRUE;
}
//=============================================================================
void PointMaker::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  // Add the histograms to the output on each slave server
  AliInfo(Form("\nNumber of tracks:\tprocessed\t%d\taccepted\t%d", fNtracks, fNAcceptedTracks));
}
//=============================================================================
void PointMaker::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  //  fTree = dynamic_cast<TTree*>(fOutput->FindObject("spTree"));
  
  //  TFile* file = TFile::Open("AliTrackPoints.root", "RECREATE");
  fFile->cd();
  fTree->Write();
  fCuttra->Write();
  fCutpoi->Write();
  fModpop->Write();
  fModpop->Draw();
}
//=============================================================================
