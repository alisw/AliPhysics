#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDv0KineCuts.h"

//#include "AliMCEventHandler.h"
//#include "AliMCEvent.h"

//#include "AliKalmanTrack.h"

#include "AliTRDpadPlane.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"

#include "AliTRDdigitsTask.h"

#include <iostream>
#include <iomanip>
using namespace std;

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliTRDdigitsTask)

//________________________________________________________________________
AliTRDdigitsTask::AliTRDdigitsTask(const char *name) 
: AliAnalysisTaskSE(name),
  fESD(0), fOutputList(0),
  fV0cuts(0), fV0tags(0), fV0electrons(0), fV0pions(0),
  fDigitsInputFileName("TRD.Digits.root"), fDigitsOutputFileName(""),
  fDigitsInputFile(0), fDigitsOutputFile(0),
  fDigMan(0),fGeo(0), fEventNoInFile(-2), fDigitsLoadedFlag(kFALSE)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

  // create the digits manager
  fDigMan = new AliTRDdigitsManager;
  fDigMan->CreateArrays();
  
  // create a TRD geometry, needed for matching digits to tracks
  fGeo = new AliTRDgeometry;
  if (! fGeo) {
    AliFatal("cannot create geometry ");
  }

}

//_______________________________________________________________________
TFile* AliTRDdigitsTask::OpenDigitsFile(TString inputfile,
					TString digfile,
					TString opt)
{
  // we should check if we are reading ESDs or AODs - for now, only
  // ESDs are supported

  if (digfile == "") {
    return NULL;
  }

  // construct the name of the digits file from the input file
  inputfile.ReplaceAll("AliESDs.root", digfile);

  // open the file
  AliInfo( "opening digits file " + inputfile
           + " with option \"" + opt + "\"");
  //TFile* dfile = new TFile(inputfile, opt);
  TFile* dfile = TFile::Open(inputfile, opt);
  if (!dfile) {
    AliWarning("digits file '" + inputfile + "' cannot be opened");
  }

  return dfile;
}


//_______________________________________________________________________
Bool_t AliTRDdigitsTask::UserNotify()
{
  delete fDigitsInputFile;
  delete fDigitsOutputFile;

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if ( ! esdH ) return kFALSE;
  if ( ! esdH->GetTree() ) return kFALSE;
  if ( ! esdH->GetTree()->GetCurrentFile() ) return kFALSE;

  TString fname = esdH->GetTree()->GetCurrentFile()->GetName();

  fDigitsInputFile  = OpenDigitsFile(fname,fDigitsInputFileName,"");
  fDigitsOutputFile = OpenDigitsFile(fname,fDigitsOutputFileName,"RECREATE");

  fEventNoInFile = -1;

  return kTRUE;
}


//________________________________________________________________________
void AliTRDdigitsTask::UserCreateOutputObjects()
{

  // create list for output (QA) stuff
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  // At this point, additional histograms can be created, maybe via a
  // virtual function.

  
  PostData(1, fOutputList);
}



//________________________________________________________________________
Bool_t AliTRDdigitsTask::NextEvent(Bool_t preload) 
{
  fEventNoInFile++;
  fDigitsLoadedFlag = kFALSE;

  if (preload) {
    return ReadDigits();
  } else {
    return kTRUE;
  }
}

//________________________________________________________________________
void AliTRDdigitsTask::UserExec(Option_t *) 
{

  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
  // IMPORTANT: call NextEvent() for book-keeping
  // -----------------------------------------------------------------
  NextEvent();
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------


  if ( ! ReadDigits() ) return;
  
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
  // The following is a generic example how to access digits and match
  // them to tracks found in the ESD


  
  // -----------------------------------------------------------------
  // prepare event data structures
  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
}

//________________________________________________________________________
void AliTRDdigitsTask::AnalyseEvent()
{

  
  printf("There are %d tracks in this event\n", fESD->GetNumberOfTracks());

  if (fESD->GetNumberOfTracks() == 0) {
    // skip empty event
    return;
  }

  // make digits available
  ReadDigits();
   
  
  
  // -----------------------------------------------------------------
  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    
    // ---------------------------------------------------------------
    // gather track information

    // we always want the ESD track
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }

    // the friend and TRD tracks are only necessary if we want to
    // match via tracklets 
    //AliESDfriendTrack* friendtrack = fESDfriend->GetTrack(iTracks);
    // I do not handle missing friend tracks - I keep checking the
    // pointer when I access it

    //AliTRDtrackV1* trdtrack = FindTRDtrackV1(friendtrack);
    AliTRDtrackV1* trdtrack = NULL;
    //if (!trdtrack) {
    //  // this happens often, because not all tracks reach the TRD
    //  printf("NOTICE: Could not receive TRD track %d\n", iTracks);
    //  continue;
    //}

    // skip boring tracks
    if (trdtrack && trdtrack->GetNumberOfTracklets() == 0) continue;
    if (track->Pt() < 1.5) continue;

    // are there tracks without outer params?
    if ( ! track->GetOuterParam() ) {
      AliWarning(Form("Track %d has no OuterParam", iTracks));
      continue;
    }
      
    // print some info about the track
    cout << " ====== TRACK " << iTracks
	 << "   pT = " << track->Pt() << " GeV";
    if (trdtrack) {
      cout << ", " << trdtrack->GetNumberOfTracklets()
	   << " tracklets";
    }
    cout << " ======" << endl;

    // look for tracklets in all 5 layers
    for (int ly=0;ly<6;ly++) {
      Int_t det=-1;
      Int_t row=-1;
      Int_t col=-1;

      Int_t det1=-1,row1=-1,col1=-1;
      Float_t x=-999.,y=-999.,z=-999.;
      if (FindDigitsTrkl(trdtrack, ly, &det,&row,&col, &x,&y,&z) ) {

	det = det1;
	row = row1;
	col = col1;
	
	cout << "    tracklet: "
	     << det1 << ":" << row1 << ":" << col1 << "   "
	     << x << " / "<< y << " / "<< z
	     << endl;
		
      }


      Int_t det2,row2,col2;
      if ( FindDigits(track->GetOuterParam(),
		      fESD->GetMagneticField(), ly,
		      &det2,&row2,&col2) ) {

	if (det>=0 && det!=det2) {
	  AliWarning("DET mismatch between tracklet and extrapolation: " 
		     + TString(Form("%d != %d", det, det2)));

	  if (row>=0 && row!=row2) {
	    AliWarning("ROW mismatch between tracklet and extrapolation: " 
		       + TString(Form("%d != %d", row, row2)));
	  }
	}
	
	det = det2;
	row = row2;
	col = col2;
	
	cout << "    outparam: "
	     << det2 << ":" << row2 << ":" << col2 << "   "
	     << track->GetOuterParam()->GetX() << " / "
	     << track->GetOuterParam()->GetY() << " / "
	     << track->GetOuterParam()->GetZ() 
	     << endl;
      }
      

      if (det>=0) {
	cout << "Found tracklet at "
	     << det << ":" << row << ":" << col << endl;
	
	int np = 5;
	if ( col-np < 0 || col+np >= 144 )
	  continue;
	
	for (int c = col-np; c<=col+np;c++) {
	  cout << "  " << setw(3) << c << " ";
	  for (int t=0; t<fDigMan->GetDigits(det)->GetNtime(); t++) {
	    cout << setw(4) << fDigMan->GetDigitAmp(row,c,t,det);
	  }
	  cout << endl;
	}
      }

      cout << endl;

    }

  }

 
  for (int det=0; det<540; det++) {

    if (!fDigMan->GetDigits(det)) {
      AliWarning(Form("No digits found for detector %d", det));
      continue;
    }
    
    AliTRDpadPlane* padplane = fGeo->GetPadPlane(det);
    if (!padplane) {
      AliError(Form("AliTRDpadPlane for detector %d not found",det));
      continue;
    }
    
//    for (int row=0; row < padplane->GetNrows(); row++) {
//      for (int col=0; col < padplane->GetNcols(); col++) {
//	for (int tb=0; tb < fDigMan->GetDigits(det)->GetNtime(); tb++) {
//	  fhTrdAdc->Fill(fDigMan->GetDigitAmp(row,col,tb,det));
//	}
//      }
//    }

  }
   
//  PostData(1, fOutputList);
}      


  
//________________________________________________________________________
AliTRDtrackV1* AliTRDdigitsTask::FindTRDtrackV1(AliESDfriendTrack* friendtrack) 
{
  if (!friendtrack) {
    //AliWarning("ERROR: Could not receive friend track");
    return NULL;
  }

  // find AliTRDtrackV1
  TObject* fCalibObject = 0;
  AliTRDtrackV1* trdTrack = 0;
  // find TRD track
  int icalib=0;
  while ((fCalibObject = (TObject*)(friendtrack->GetCalibObject(icalib++)))){
    if(strcmp(fCalibObject->IsA()->GetName(), "AliTRDtrackV1") != 0)
      continue;
    trdTrack = (AliTRDtrackV1 *)fCalibObject;
  }

  return trdTrack;
}


//________________________________________________________________________
Int_t AliTRDdigitsTask::FindDigits(const AliExternalTrackParam* paramp,
				   Float_t bfield, Int_t layer,
				   Int_t* det, Int_t* row, Int_t* col)
{

  if ( ! paramp ) return kFALSE;
  
  // create a copy of the track params that we can propagate around
  AliExternalTrackParam par = *paramp;

  // find the sector
  Int_t sector = int(floor(par.GetAlpha()/TMath::Pi()*9.));
  if (sector<0) sector += 18;

  if ( ! par.PropagateTo(fGeo->GetTime0(layer), bfield) ) {
    AliWarning("Failed to propagate track param");
    return kFALSE;
  }
      
  
  Int_t st  = fGeo->GetStack(par.GetZ(),layer);
  if (st<0) {
    return kFALSE;
  }
  
  *det = 30*sector + 6*st + layer;
  AliTRDpadPlane* padplane = fGeo->GetPadPlane(layer,st);
  *row = padplane->GetPadRowNumber(par.GetZ());
  *col = padplane->GetPadColNumber(par.GetY());

  if (*row == -1) return kFALSE;
  if (*col == -1) return kFALSE;
  
  return kTRUE;
}


//________________________________________________________________________
Int_t AliTRDdigitsTask::FindDigitsTrkl(const AliTRDtrackV1* trdTrack,
				       Int_t layer,
				       Int_t* det, Int_t* row, Int_t* col,
				       Float_t* x, Float_t* y, Float_t* z)
{

  if ( ! trdTrack ) return kFALSE;
  
  // loop over tracklets
  for(Int_t itr = 0; itr < 6; ++itr) {
    
    AliTRDseedV1* tracklet = 0;
    
    if(!(tracklet = trdTrack->GetTracklet(itr)))
      continue;
    if(!tracklet->IsOK())
      continue;
    
    if ( tracklet->GetDetector()%6 == layer ) {

      AliTRDpadPlane *padplane = fGeo->GetPadPlane(tracklet->GetDetector());

      *det = tracklet->GetDetector();
      *row = padplane->GetPadRowNumber(tracklet->GetZ());
      *col = padplane->GetPadColNumber(tracklet->GetY());

      if (x!=NULL) *x = tracklet->GetX();
      if (y!=NULL) *y = tracklet->GetY();
      if (z!=NULL) *z = tracklet->GetZ();

      return kTRUE;
    }
  }

  return kFALSE; // no tracklet found
  
}  
    
//    // AliTrackPoint tp;
//    // for (int ip=0; ip<array->GetNPoints(); ip++) {
//    //   array->GetPoint(tp,ip);
//    //   cout << "  point " << ip << ": "
//    // 	   << tp.GetVolumeID() << "  -  "
//    // 	   << tp.GetX() << "/" << tp.GetY() << "/" << tp.GetZ()
//    // 	   << "   r = " << TMath::Hypot(tp.GetX(),tp.GetY())
//    // 	   << endl;
//    // }
//    
//
//    
//
//
//    
//    fHistPt->Fill(track->Pt());
//  } //track loop 
//
//  delete geo;
//


//______________________________________________________________________________
void AliTRDdigitsTask::FillV0PIDlist(){

  // no need to run if the V0 cuts are not set
  if (!fV0cuts) {
    cout << "FillV0PIDlist: skip event - no V0 cuts" << endl;
    return;
  }
  
  // basic sanity check...
  if (!fESD) {
    cout << "FillV0PIDlist: skip event - no ESD event" << endl;
    return;
  }

  const Int_t numTracks = fESD->GetNumberOfTracks();

  if (numTracks < 0) {
    AliFatal("negative number of tracks?!?");
  }
  
  // Ensure there is sufficient memory for the V0 tags
  if (!fV0tags) {
    //cout << "FillV0PIDlist: create fV0tags" << endl;
    fV0tags = new Int_t[numTracks];
  } else if ( sizeof(fV0tags)/sizeof(Int_t) < numTracks ) {
    //cout << "FillV0PIDlist: re-create fV0tags" << endl;
    delete fV0tags;
    fV0tags = new Int_t[numTracks];
  } else {
    // there is more than enough space to store the tags, no need to
    // do anything.
  }

  // Reset the V0 tags and reference particle arrays
  for (Int_t i = 0; i < numTracks; i++) {
    fV0tags[i] = 0;
  }

  fV0electrons->Clear();
  fV0pions->Clear();

  
  // V0 selection
  // loop over V0 particles
  fV0cuts->SetEvent(fESD);
  for(Int_t iv0=0; iv0<fESD->GetNumberOfV0s();iv0++){
    
    AliESDv0 *v0 = (AliESDv0 *) fESD->GetV0(iv0);

    if(!v0) continue;
    if(v0->GetOnFlyStatus()) continue;

    // Get the particle selection
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0, pdgP, pdgN;
    foundV0 = fV0cuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if(!foundV0) continue;

    Int_t iTrackP = v0->GetPindex();  // positive track inded
    Int_t iTrackN = v0->GetNindex();  // negative track

    if (fV0tags[iTrackP]) {
      printf("Warning: particle %d tagged more than once\n", iTrackP);
    }

    if (fV0tags[iTrackN]) {
      printf("Warning: particle %d tagged more than once\n", iTrackN);
    }

    
    // fill the Object arrays
    // positive particles
    if( pdgP == -11){
      //DigitsDictionary(iTrackP, iv0, pdgP);
      fV0electrons->Add(fESD->GetTrack(iTrackP));
      fV0tags[iTrackP] = pdgP;
    }
    else if( pdgP == 211){
      //DigitsDictionary(iTrackP, iv0, pdgP);
      fV0pions->Add(fESD->GetTrack(iTrackP));
      fV0tags[iTrackP] = pdgP;
    }


    // negative particles
    if( pdgN == 11){
      //DigitsDictionary(iTrackN, -iv0, pdgN);
      fV0electrons->Add(fESD->GetTrack(iTrackN));
      fV0tags[iTrackN] = pdgN;
    }
    else if( pdgN == -211){
      //DigitsDictionary(iTrackN, -iv0, pdgN);
      fV0pions->Add(fESD->GetTrack(iTrackN));
      fV0tags[iTrackN] = pdgN;
    }

  }

}



//________________________________________________________________________
void AliTRDdigitsTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

//  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
//  if (!fOutputList) {
//    printf("ERROR: Output list not available\n");
//    return;
//  }
//  
//  fHistPt = dynamic_cast<TH1F*> (fOutputList->At(0));
//  if (!fHistPt) {
//    printf("ERROR: fHistPt not available\n");
//    return;
//  }
//   
  //TCanvas *c1 = new TCanvas("AliTRDdigitsTask","Pt",10,10,510,510);
  //c1->cd(1)->SetLogy();
  //fHistPt->DrawCopy("E");
}


//________________________________________________________________________
Bool_t AliTRDdigitsTask::ReadDigits()
{

  // don't do anything if the digits have already been loaded
  if (fDigitsLoadedFlag) return kTRUE;
  
  if (!fDigMan) {
    AliError("no digits manager");
    return kFALSE;
  }
  
  // reset digit arrays
  for (Int_t det=0; det<540; det++) {
    fDigMan->ClearArrays(det);
    fDigMan->ClearIndexes(det);
  }


  if (!fDigitsInputFile) {
    AliError("digits file not available");
    return kFALSE;
  }

  
  // read digits from file
  TTree* tr = (TTree*)fDigitsInputFile->Get(Form("Event%d/TreeD",
                                                 fEventNoInFile));

  if (!tr) {
    AliWarning(Form("digits tree for event %d not found", fEventNoInFile));
    return kFALSE;
  }

  fDigMan->ReadDigits(tr);
  delete tr;

  // expand digits for use in this task
  for (Int_t det=0; det<540; det++) {
    if (fDigMan->GetDigits(det)) {
      fDigMan->GetDigits(det)->Expand();
    }
  }

  fDigitsLoadedFlag = kTRUE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliTRDdigitsTask::WriteDigits()
{
  // check for output file
  if (!fDigitsOutputFile) {
    AliError("digits output file not available");
    return kFALSE;
  }
  
  // compress digits for storage
  for (Int_t det=0; det<540; det++) {
    fDigMan->GetDigits(det)->Expand();
  }

  // create directory to store digits tree
  TDirectory* evdir =
    fDigitsOutputFile->mkdir(Form("Event%d", fEventNoInFile),
                             Form("Event%d", fEventNoInFile));

  evdir->Write();
  evdir->cd();

  // save digits tree
  TTree* tr = new TTree("TreeD", "TreeD");
  fDigMan->MakeBranch(tr);
  fDigMan->WriteDigits();
  delete tr;

  return kTRUE;
}

