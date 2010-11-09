#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"

#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliLog.h"
#include "AliESDTrdTrack.h"

#include "AliTRDtrackletMCM.h"
#include "AliVParticle.h"
#include "AliMCParticle.h" 
#include "AliTrackReference.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTRDtrackletWord.h"

#include "AliTRDonlineTrackletQA.h"

ClassImp(AliTRDonlineTrackletQA)

AliTRDonlineTrackletQA::AliTRDonlineTrackletQA(const char *name) :
  AliAnalysisTask(name, ""),
  fESD(0x0),
  fInputHandler(0x0),
  fInputEvent(0x0),
  fOutputAOD(0x0),
  fMCEvent(0x0),
  fTrackletsRaw(0x0),
  fTrackletsSim(0x0), 
  fOutputList(0x0),
  fHistYpos(0x0),
  fHistYres(0x0),
  fHistYresDy(0x0),
  fHistdY(0x0),
  fHistdYres(0x0),
  fHistYresESD(0x0),
  fHistdYresESD(0x0),
  fHistCanddY(0x0),
  fHistFounddY(0x0),
  fHistTrklPerRef(0x0),
  fHistdYdYref(0x0),
  fHistYposRaw(0x0),
  fHistdYRaw(0x0),
  fHistYdYRaw(0x0), 
  fHistZrow(0x0),
  fHistZrowRaw(0x0),
  fHistPid(0x0),
  fHistPidRaw(0x0),
  fHistYdiff(0x0),
  fHistdYdiff(0x0),
  fHistdYdYraw(0x0),
  fTreeTracklets(0x0),
  fY(0.),
  fDY(0.),
  fYdiff(0.),
  fDYdiff(0.),
  fQ0(0),
  fQ1(0),
  fNHits(0),
  fMinPt(1.),
  fGeo(new AliTRDgeometry),
  fNevent(0),
  fTrackletTree(0x0),
  fTrackletTreeRaw(0x0)
{
  // ctor

  DefineInput(0, TChain::Class());
  DefineInput(1, TTree::Class());

  DefineOutput(0, TTree::Class()); 
  DefineOutput(1, TList::Class());
}

AliTRDonlineTrackletQA::~AliTRDonlineTrackletQA()
{
  // dtor

  delete fGeo;
}

void AliTRDonlineTrackletQA::ConnectInputData(const Option_t * /* option */)
{
  // connect input data

  fInputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (fInputHandler)
    fInputEvent = fInputHandler->GetEvent();

  AliMCEventHandler *mcH = (AliMCEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
  if (mcH)
    fMCEvent = mcH->MCEvent();
}

void AliTRDonlineTrackletQA::CreateOutputObjects()
{
  // create output objects

  OpenFile(1); 
  
  fOutputList = new TList();

  fHistYpos        = new TH1F("ypos", 
			      "Tracklet (sim) y-position;y (cm);count", 
			      8192/32, -4096*160e-4, 4095*160e-4);
  fHistYposRaw     = new TH1F("ypos-raw", 
			      "Tracklet (raw) y-position;y (cm);count", 
			      130, -65, 65);
  fHistYres        = new TH1F("yres", 
			      "Tracklet (sim) #Deltay;y_{tracklet}-y_{MC} (cm);count", 
			      8192/32, -4096/32*160e-4, 4095/32*160e-4);
  fHistYresDy      = new TH2F("yresdy", 
			      "Tracklet (sim) #Deltay;y_{tracklet}-y_{MC} (cm);deflection (bin)", 
			      8192/32, -4096/32*160e-4, 4095/32*160e-4,
			      128, -64.5, 63.5);
  fHistYresESD     = new TH1F("yresesd", 
			      "Tracklet #Deltay;y (cm);count", 
			      100, -10, 10);
  fHistYdiff       = new TH1F("ydiff", 
			      "Tracklet #Deltay (sim - raw);y_{sim}-y_{raw} (cm);count", 
			      8192/32, -4096/32*160e-4, 4095/32*160e-4);
  for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
    fHistYlocal[iLayer] = new TH2F(Form("ylocal_%i", iLayer), 
				   Form("Tracklet local y, layer %i;y_{MC} (pad width);y_{trkl} (pad width)", iLayer),
				   100, -1, 1, 100, -1, 1);
  }

  fHistdY          = new TH1F("dy", 
			      "deflection (sim);dy (140 #mum)", 
			      128, -64.5, 63.5);
  fHistdYRaw       = new TH1F("dy-raw", 
			      "deflection (sim);dy (140 #mum)", 
			      128, -64.5, 63.5);
  fHistdYres       = new TH1F("dyres", 
			      "deflection residual;dy (cm)", 
			      128, -1., 1.);
  fHistdYresESD    = new TH1F("dyresesd", 
			      "deflection residual;dy (cm)", 
			      128, -1., 1.);
  fHistCanddY      = new TH1F("dycand", 
			      "deflection;dy (140 #mum)", 
			      128, -64.5, 63.5);
  fHistFounddY     = new TH1F("dyfound", 
			      "deflection;dy (140 #mum)", 
			      128, -64.5, 63.5);
  fHistdYdiff      = new TH1F("dydiff", 
			      "deflection #Deltady;dy_{sim}-dy_{raw} (140 #mum)", 
			      128, -64.5, 63.5);
  fHistdYdYraw     = new TH2F("dydyraw", 
			      "deflection from sim. vs raw;dy_{sim} (140 #mum);dy_{raw} (140 #mum)",
			      128, -64.5, 63.5, 128, -64.5, 63.5);

  fHistTrklPerRef  = new TH1F("trklperref", 
			      "No. of tracklets per track reference;no. of tracklets",
			      10, -0.5, 9.5);

  fHistdYdYref     = new TH2F("dydyref", 
			      "deflection vs. deflection from track reference;dy_{ref} (140 #mum);dy (140 #mum)",
			      128, -64.5, 63.5, 128, -64.5, 63.5);

  fHistZrow        = new TH1F("zrow", 
			      "z-position;pad row", 
			      16, -0.5, 15.5);
  fHistZrowRaw     = new TH1F("zrow-raw", 
			      "z-position;pad row", 
			      16, -0.5, 15.5);

  fHistPid         = new TH1F("pid", 
			      "pid", 
			      256, -0.5, 255.5);
  fHistPidRaw      = new TH1F("pid-raw", 
			      "pid", 
			      256, -0.5, 255.5);

  fHistYdYRaw      = new TH2F("ydyraw", 
			      "y vs dy (raw tracklets);y (cm);dy (140 #mum)", 
			      8192/32, -4096*160e-4, 4095*160e-4, 
			      128, -64.5, 63.5);

  fTreeTracklets   = new TTree("trkl", "trkl");
  fTreeTracklets->Branch("y", &fY);
  fTreeTracklets->Branch("dy", &fDY);
  fTreeTracklets->Branch("ydiff", &fYdiff);
  fTreeTracklets->Branch("dydiff", &fDYdiff);
  fTreeTracklets->Branch("q0", &fQ0);
  fTreeTracklets->Branch("q1", &fQ1);
  fTreeTracklets->Branch("nhits", &fNHits);

  fOutputList->Add(fHistYpos);
  fOutputList->Add(fHistdY);
  fOutputList->Add(fHistZrow);
  fOutputList->Add(fHistPid);

  fOutputList->Add(fHistYres);
  fOutputList->Add(fHistYresDy);
  fOutputList->Add(fHistCanddY);
  fOutputList->Add(fHistFounddY);
  fOutputList->Add(fHistTrklPerRef);
  fOutputList->Add(fHistdYres);
  fOutputList->Add(fHistYresESD);
  fOutputList->Add(fHistdYresESD);
  fOutputList->Add(fHistdYdYref);

  for (Int_t iLayer = 0; iLayer < 6; iLayer++) 
    fOutputList->Add(fHistYlocal[iLayer]);

  fOutputList->Add(fHistYposRaw);
  fOutputList->Add(fHistdYRaw);
  fOutputList->Add(fHistZrowRaw);
  fOutputList->Add(fHistPidRaw);
  fOutputList->Add(fHistYdYRaw);

  fOutputList->Add(fHistYdiff);
  fOutputList->Add(fHistdYdiff);
  fOutputList->Add(fHistdYdYraw);

  fOutputList->Add(fTreeTracklets);
}

void AliTRDonlineTrackletQA::Exec(const Option_t * /* option */)
{
  // execute this for each event

  // connect input data
  fTrackletTree = (TTree*) GetInputData(1);
  if (!fTrackletTree)
    return;

  fTrackletTree->SetBranchAddress("tracklets_sim", &fTrackletsSim);
  fTrackletTree->SetBranchAddress("tracklets_raw", &fTrackletsRaw);
  fTrackletTree->GetEntry(fTrackletTree->GetEntriesFast()-1);

  // prepare raw tracklets for comparison
  Int_t detRaw;
  Int_t robRaw; 
  Int_t mcmRaw;
  Int_t yRaw; 
  Int_t dyRaw;
  TTree trklRaw("raw tracklets", "raw tracklets");
  trklRaw.Branch("det", &detRaw);
  trklRaw.Branch("rob", &robRaw);
  trklRaw.Branch("mcm", &mcmRaw);
  trklRaw.Branch("y", &yRaw);
  trklRaw.Branch("dy", &dyRaw);
  trklRaw.SetDirectory(0x0);
  // prepare simulated tracklets for comparison
  Int_t detSim;
  Int_t robSim; 
  Int_t mcmSim;
  Int_t ySim; 
  Int_t dySim;
  TTree trklSim("sim tracklets", "sim tracklets");
  trklSim.Branch("det", &detSim);
  trklSim.Branch("rob", &robSim);
  trklSim.Branch("mcm", &mcmSim);
  trklSim.Branch("y", &ySim);
  trklSim.Branch("dy", &dySim);
  trklSim.SetDirectory(0x0);

  // ----- simulated tracklets -----
  AliTRDtrackletMCM *trkl = 0x0;
  if (fTrackletsSim) {
    for (Int_t iTracklet = 0; iTracklet < fTrackletsSim->GetEntries(); iTracklet++) {
      trkl = (AliTRDtrackletMCM*) (*fTrackletsSim)[iTracklet];
//      Int_t label = trkl->GetLabel();
//      if (label > -1 && label < maxTracks) 
//	mcTrackToTrackletMCM[label].idx[mcTrackToTrackletMCM[label].size < 10 ? mcTrackToTrackletMCM[label].size++ : 0] = iTracklet;
      fHistYpos->Fill(trkl->GetY());
      fHistdY->Fill(trkl->GetdY());
      fHistZrow->Fill(trkl->GetZbin());
      fHistPid->Fill(trkl->GetPID());

      detSim = trkl->GetDetector();
      robSim = trkl->GetROB();
      mcmSim = trkl->GetMCM();
      ySim   = trkl->GetYbin();
      dySim  = trkl->GetdY();
      trklSim.Fill();

      PlotMC(trkl);
      PlotESD(trkl);
    }
  }

  // ----- raw tracklets -----
  if (fTrackletsRaw) {
    for (Int_t iTracklet = 0; iTracklet < fTrackletsRaw->GetEntries(); iTracklet++) {
      AliTRDtrackletWord *trklWord = (AliTRDtrackletWord*) (*fTrackletsRaw)[iTracklet];
      // remove unwanted chambers
      if (trklWord->GetDetector() == 57 ||
	  trklWord->GetDetector() == 47 ||
	  trklWord->GetDetector() == 32)
	continue;

      fHistYposRaw->Fill(trklWord->GetY());
      fHistdYRaw->Fill(trklWord->GetdY());
      fHistZrowRaw->Fill(trklWord->GetZbin());
      fHistPidRaw->Fill(trklWord->GetPID());
      fHistYdYRaw->Fill(trklWord->GetY(), trklWord->GetdY());

      detRaw = trklWord->GetDetector();
      robRaw = trklWord->GetROB();
      mcmRaw = trklWord->GetMCM();
      yRaw   = trklWord->GetYbin();
      dyRaw  = trklWord->GetdY();
      trklRaw.Fill();
    }
  }

  // ----- tracklet comparison raw to simulated -----
  trklRaw.BuildIndex("(det+rob*540)*16+mcm", "det");
  trklSim.BuildIndex("(det+rob*540)*16+mcm", "det");
  trklRaw.AddFriend(&trklSim, (const char*) "sim");
  gDirectory->Add(fHistYdiff, kFALSE);
  Int_t ncomp = trklRaw.Draw("(sim.y-y)*160e-4>>+ydiff", "", "goff");
  printf("----- Compared %i tracklets -----\n", ncomp);
  gDirectory->Remove(fHistYdiff);
  gDirectory->Add(fHistdYdiff, kFALSE);
  trklRaw.Draw("(sim.dy-dy)>>+dydiff", "", "goff");
  gDirectory->Remove(fHistdYdiff);
  gDirectory->Add(fHistdYdYraw, kFALSE);
  trklRaw.Draw("dy:sim.dy>>+dydyraw", "", "goff");
  //  trklRaw.Scan("det:rob:mcm:y:dy:sim.dy", "sim.dy < 30 && dy > 30");
  gDirectory->Remove(fHistdYdYraw);

  // ----- MC tracks and track references -----
  // determine tracklet efficiency
  if (fMCEvent) {
    Int_t nTracksMC = fMCEvent->GetNumberOfTracks();
    AliInfo(Form("no. of MC tracks: %i", nTracksMC));
    for (Int_t iTrack = 0; iTrack < nTracksMC; iTrack++) {
      // we want primaries
      if (!fMCEvent->IsPhysicalPrimary(iTrack))
	continue;

      AliMCParticle *mcpart = (AliMCParticle*) fMCEvent->GetTrack(iTrack);

      // don't look at tracks with too low pt
      if (TMath::Abs(mcpart->Pt()) < fMinPt)
	continue;

      // look for two track references in a chamber
      Int_t nTrackRefs = mcpart->GetNumberOfTrackReferences();
      AliTrackReference *tr[2] = { 0x0 };
      Int_t nRef = 0;

      for (Int_t iTrackRef = 0; iTrackRef < nTrackRefs; iTrackRef++) {
	AliTrackReference *trackRef = mcpart->GetTrackReference(iTrackRef);
	if (trackRef->DetectorId() != AliTrackReference::kTRD)
	  continue;
	if (trackRef->Pt() < fMinPt)
	  continue;
	Int_t label = trackRef->Label();
	if (label < 0)
	  continue;

	// first track reference, remember it
	if (nRef == 0) {
	  tr[nRef] = trackRef;
	  nRef++;
	  continue;
	}
	else {
	  // next one is too far away, remember it but forget the previous one
	  if (TMath::Abs(trackRef->LocalX() - tr[0]->LocalX()) > 5.) {
	    tr[0] = trackRef;
	    nRef = 1;
	    continue;
	  }
	  // too close to previous track reference
	  // we don't want it
	  else if (TMath::Abs(trackRef->LocalX() - tr[0]->LocalX()) < .5) {
	    continue;
	  }
	  // then it must be ok, so we take it
	  else {
	    tr[1] = trackRef;
	    nRef++;
	  }
	}

	// calculation deflection from track references
	Float_t deflLength = 3. * (tr[1]->LocalY() - tr[0]->LocalY()) / (tr[1]->LocalX() - tr[0]->LocalX());
	// if it is too large we reject it
	if (deflLength < 1.) {
	  fHistCanddY->Fill(deflLength/140e-4);
	}
	else {
	  nRef = 0;
	  continue;
	}
	
	// now search for tracklets belonging to this track reference
	Int_t nTrackletsPerRef = 0;
	Int_t defl = 0.;
	for (Int_t iTracklet = 0; iTracklet < fTrackletsSim->GetEntries(); iTracklet++) {
	  trkl = (AliTRDtrackletMCM*) (*fTrackletsSim)[iTracklet];
	  // they must have the same label
	  if (!trkl->HasLabel(label)) 
	    continue;
	  // and be close enough in radial position
	  if (TMath::Abs(trackRef->LocalX() - trkl->GetX()) > 5.) 
	    continue;
	  // if they are close in position we accept it
	  if ((TMath::Abs(trackRef->LocalY() - trkl->GetY()) < 5.) &&
	      (TMath::Abs(trackRef->Z() - trkl->GetZ()) < 5.)) {
	    defl = trkl->GetdY();
	    nTrackletsPerRef++;
	  }
	}
	fHistTrklPerRef->Fill(nTrackletsPerRef);
	if (nTrackletsPerRef == 0) {
	  AliInfo(Form("Track ref without assigned tracklet: x=%4.2f, y=%4.2f, z=%4.2f, pt=%4.2f (%i)",
		       trackRef->X(), trackRef->Y(), trackRef->Z(), trackRef->Pt(), trackRef->Label()));
	}
	if (nTrackletsPerRef == 1) {
	  fHistdYdYref->Fill(deflLength/140.e-4, defl);
	  fHistFounddY->Fill(deflLength/140.e-4);
	}
	nRef = 0;
      }
    }
  }

  // ----- ESD tracks -----
  if (fESD) {
    for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
      AliESDtrack *esdTrack = fESD->GetTrack(iTrack);
      AliDebug(1, Form("ESD track pt: %7.2f", esdTrack->Pt()));
    }
  }

  PostData(1, fOutputList);  
}

void AliTRDonlineTrackletQA::LocalInit()
{

}

void AliTRDonlineTrackletQA::Terminate(const Option_t * /* option */)
{
  // upon terminate

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));

  if (!fOutputList) {
    AliError("No output objects found!");
    return;
  }

//  fHistYpos = dynamic_cast<TH1F*> (fOutputList->FindObject("ypos"));
//  if (fHistYpos) {
//    TCanvas *c = new TCanvas;
//    c->cd();
//    gPad->SetGridx();
//    gPad->SetGridy();
//    fHistYpos->DrawCopy();
//  }
}

void AliTRDonlineTrackletQA::PlotMC(AliTRDtrackletMCM *trkl)
{
  // compare the given tracklet to the MC information,
  // i.e. track references
  
  Int_t label = trkl->GetLabel();

  // check for valid label
  if (label < 0 ) {
    AliWarning("MC tracklet has no label");
    return;
  }
  if (label >= fMCEvent->GetNumberOfTracks()) {
    AliWarning("MC tracklet has invalid label");
    return;
  }

  // some cuts on the tracklet
//  if (!fMCEvent->IsPhysicalPrimary(label))
//    return;
//  if (TMath::Abs(trkl->GetdYdX()) > 0.05)
//    return;
//  if (trkl->GetDetector() % 6 != 5)
//    return;
  
  // get MC particle for this tracklet
  AliMCParticle *track = (AliMCParticle*) fMCEvent->GetTrack(label);

  // don't look at tracks with too low pt
  if (TMath::Abs(track->Pt() < fMinPt))
    return;

  // serach track references corresponding to the current tracklet
  AliTrackReference *tr[2] = { 0x0 };
  Int_t nTrackRefs = 0;

  for (Int_t iTrackRef = 0; iTrackRef < track->GetNumberOfTrackReferences(); iTrackRef++) {
    AliTrackReference *trackRef = track->GetTrackReference(iTrackRef);
    if (trackRef->DetectorId() != AliTrackReference::kTRD)
      continue;
    if (trackRef->Pt() < fMinPt)
      continue;
    
    if (TMath::Abs(trkl->GetX() - trackRef->LocalX()) > 5.)
      continue;

    tr[nTrackRefs++] = trackRef; 

    if (nTrackRefs == 2)
      break;
  }

  // if there were exactly 2 track references 
  // (otherwise something is strange and we want to look at clean cases)
  // compare tracklet to them
  if (nTrackRefs == 2) {
    // sanity check to be in the same sector
    if ( TMath::Abs((tr[0]->Alpha()*180./TMath::Pi()-10.)/20. - (trkl->GetDetector()/30)) > .1) {
      AliError("Track reference in different sector");
    }
    // require minimal distance in X and maximum deflection in Y
    // for the track references
    else if ((tr[1]->LocalX() - tr[0]->LocalX()) > 0.1 && TMath::Abs(tr[1]->LocalY() - tr[0]->LocalY()) < 1.) { 
      // calculate slope from track references
      // and check whether it's in the allowed range
      Float_t slope = 3. * (tr[1]->LocalY() - tr[0]->LocalY()) / (tr[1]->LocalX() - tr[0]->LocalX());
      if (TMath::Abs(slope) < 64*140e-4) {
	AliDebug(1,Form("x1: %f, x0: %f, y1: %f, y0:%f", 
			tr[1]->LocalX(), tr[0]->LocalX(), tr[1]->LocalY(), tr[0]->LocalY() ));
	// calculate y-position scaled to radial position of the tracklet
	// and consider the tilting angle of the pads
	// since the tracklets are affected by it
	Float_t yMC     = (tr[1]->LocalY() + (-0.5+trkl->GetX() - tr[1]->LocalX()) * (tr[1]->LocalY() - tr[0]->LocalY()) / (tr[1]->LocalX() - tr[0]->LocalX()));
	Float_t yMCtilt = yMC + (TMath::Tan(TMath::Power(-1, (trkl->GetDetector() % 6))*2.*TMath::Pi()/180.) * (tr[1]->Z() - trkl->GetZ()));
	if (TMath::Abs(trkl->GetY() - yMCtilt) > 10.) {
	  AliError(Form("Deviation too large for tracklet: 0x%08x in det. %i at x = %f, y = %f, z = %f, alpha = %f", 
			trkl->GetTrackletWord(), trkl->GetDetector(), trkl->GetX(), trkl->GetY(), trkl->GetZ(), tr[0]->Alpha()));
	}
	fHistYres->Fill(trkl->GetY() - yMCtilt);
	fHistYresDy->Fill(trkl->GetY() - yMCtilt, trkl->GetdY());
	// what about tilt correction here ???
	fHistdYres->Fill(3. * trkl->GetdYdX() - 
			 3. * (tr[1]->LocalY() - tr[0]->LocalY()) / (tr[1]->LocalX() - tr[0]->LocalX()));
	// plot position deviation in pad-coordinates
	// to study the influence of the position LUT
	Float_t padWidth = fGeo->GetPadPlane(trkl->GetDetector())->GetWidthIPad();
	Float_t yMClocal = yMCtilt/padWidth - floor(yMCtilt/padWidth) - padWidth/2.; 
	Int_t layer = trkl->GetDetector() % 6;
	fHistYlocal[layer]->Fill(yMClocal, 
				 trkl->GetY()/padWidth - floor(trkl->GetY()/padWidth) - padWidth/2. - yMClocal);
	// and fill everything to the tree
	fQ0 = trkl->GetQ0();
	fQ1 = trkl->GetQ1();
	fNHits = trkl->GetNHits();
	fYdiff = trkl->GetY() - yMCtilt;
	fDYdiff = 3. * trkl->GetdYdX() -
	  3. * (tr[1]->LocalY() - tr[0]->LocalY()) / (tr[1]->LocalX() - tr[0]->LocalX());
	fY = trkl->GetY();
	fDY = trkl->GetdYdX();
	fTreeTracklets->Fill();
	// output tracklets with large deviation
	if (TMath::Abs(fYdiff) > 0.5) {
	  printf("tracklet: y=%4.2f, dy=%4.2f, ydiff=%4.2f, dydiff=%4.2f, q0=%5d, q1=%5d, nhits=%2d, label=%i\n",
		 trkl->GetY(), trkl->GetdYdX(), fYdiff, fDYdiff, fQ0, fQ1, fNHits, label);
	}
      }
    }
  }
}


void AliTRDonlineTrackletQA::PlotESD(AliTRDtrackletMCM *trkl)
{
  // plot comparison to ESD

  Float_t xTrkl = trkl->GetX();
  Float_t yTrkl = trkl->GetY();
  Float_t zTrkl = trkl->GetZ();

  Float_t alpha = (trkl->GetDetector() / 30) * 20. + 10.;
  alpha *= TMath::Pi() / 180.;

  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*> (fInputEvent);
  if (!esdEvent)
    return;

  Float_t mag = ((AliMagF*) TGeoGlobalMagField::Instance()->GetField())->SolenoidField();

  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *track = esdEvent->GetTrack(iTrack);

    if (!track->GetOuterParam())
      continue;

    AliExternalTrackParam *param = new AliExternalTrackParam(*(track->GetOuterParam()));

    AliDebug(10, Form("track %i at x = %f, y = %f", 
		      iTrack, param->GetX(), param->GetY()));
    param->Propagate(alpha, xTrkl, mag);
    AliDebug(10, Form("after propagating track %i at x = %f, y = %f", 
		      iTrack, param->GetX(), param->GetY()));

    if ((TMath::Abs(xTrkl - param->GetX()) < 10.) &&
	(TMath::Abs(yTrkl - param->GetY()) < 5.) &&
	(TMath::Abs(zTrkl - param->GetZ()) < 10.)) {
      AliInfo(Form("match of tracklet-track: %i <-> %i", 
		   trkl->GetLabel(), track->GetLabel())); 
      AliDebug(5, Form("tracklet position: det: %3i  x = %f, y = %f, z = %f, alpha = %f",
		       trkl->GetDetector(), trkl->GetX(), trkl->GetY(), trkl->GetZ(), alpha));
      AliDebug(5, Form("after propagating track %i at x = %f, y = %f, z = %f", 
		      iTrack, param->GetX(), param->GetY(), param->GetZ()));

      fHistYresESD->Fill(yTrkl - param->GetY());
    }

    delete param;
  }

}

void AliTRDonlineTrackletQA::PlotESD(AliTRDtrackletWord *trkl) 
{
  // plot comparison to ESD

  Float_t xTrkl = trkl->GetX();
  Float_t yTrkl = trkl->GetY();
  Float_t zTrkl = trkl->GetZ();

  Float_t alpha = (trkl->GetDetector() / 30) * 20. + 10.;
  alpha *= TMath::Pi() / 180.;

  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*> (fInputEvent);
  if (!esdEvent)
    return;

  Float_t mag = ((AliMagF*) TGeoGlobalMagField::Instance()->GetField())->SolenoidField();

  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *track = esdEvent->GetTrack(iTrack);

    if (!track->GetOuterParam())
      continue;

    AliExternalTrackParam *param = new AliExternalTrackParam(*(track->GetOuterParam()));

    AliDebug(10, Form("track %i at x = %f, y = %f", 
		      iTrack, param->GetX(), param->GetY()));
    param->Propagate(alpha, xTrkl, mag);
    AliDebug(10, Form("after propagating track %i at x = %f, y = %f", 
		      iTrack, param->GetX(), param->GetY()));

    if ((TMath::Abs(xTrkl - param->GetX()) < 10.) &&
	(TMath::Abs(yTrkl - param->GetY()) < 5.) &&
	(TMath::Abs(zTrkl - param->GetZ()) < 10.)) {
      AliDebug(5, Form("tracklet position: det: %3i  x = %f, y = %f, z = %f, alpha = %f",
		       trkl->GetDetector(), trkl->GetX(), trkl->GetY(), trkl->GetZ(), alpha));
      AliDebug(5, Form("after propagating track %i at x = %f, y = %f, z = %f", 
		      iTrack, param->GetX(), param->GetY(), param->GetZ()));

      fHistYresESD->Fill(yTrkl - param->GetY());
    }

    delete param;
  }

}


Int_t AliTRDonlineTrackletQA::GetTrackletsForMC(Int_t /* label */, Int_t /*idx*/ [])
{
  // get tracklets for MC label
  // not implemented

  return 0;
}
