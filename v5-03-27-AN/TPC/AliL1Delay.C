
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// This is a macro which measures the time delay of L1 trigger in TPC      //
//                                                                         //
// It reads the ESD tracks and fills histograms with the distance          //
// between the primary vertex and the Z position                           //
// of the point of closest approach between the track and the beam axis.   //
// The processing of forward and backward tracks is done separately.       //
// The corresponding histograms are fitted with gaussian and the           //
// shifts of the Z positions are extracted.                                //
// The L1 time delay is taken as an average of the forward and backward    //
// (with opposite sign) shifts divided by the drift velocity.              //
// We assume that the ESD tracks are obtained by the TPC tracking stand    //
// alone.                                                                  //
// The propagation to the vertex can be done either with or without        //
// TGeoManager.                                                            //
//                                                                         //
// Cvetan.Cheshkov@cern.ch                                                 //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TMath.h>
  #include <TError.h>
  #include <Riostream.h>
  #include <TH1F.h>
  #include <TF1.h>
  #include <TTree.h>
  #include <TCanvas.h>
  #include <TFile.h>
  #include <TStyle.h>
  #include <TGeoManager.h>
  #include "AliRunLoader.h"
  #include "AliRun.h"
  #include "AliESD.h"
  #include "AliTracker.h"
  #include "AliTPCParam.h"
  #include "AliMagF.h"
  #include "AliITStrackV2.h"
#endif

Int_t CorrectForDeadZoneMaterial(AliITStrackV2 *t);
Bool_t PropagateToVertex(AliESDtrack *track, Float_t *vtxXYZ);
Bool_t PropagateToVertexG(AliESDtrack *track, Float_t *vtxXYZ);

void AliL1Delay(const Char_t *esdfilename = "./AliESDs.root", const Char_t *galicefilename = "./galice.root", Bool_t withTGeo = kFALSE, const Char_t *geomfilename = "./geometry.root")
{
  const Int_t kMinNClusters = 0;
  const Double_t kPtMin = 0.5;
  const Double_t kMaxImpact = 0.3;
  const Double_t kDeltaZRes = 0.3;
  const Double_t kZRange = 20.;

  // Book histograms
  Int_t NBins = (Int_t)(2.*kZRange/(0.5*kDeltaZRes));
  TH1F *hDeltaZForward = new TH1F("hDeltaZForward","Longitudinal distance to the vertex (forward tracks)",NBins,-kZRange,kZRange); 
  TH1F *hDeltaZBackward = new TH1F("hDeltaZBackward","Longitudinal distance to the vertex (backward tracks)",NBins,-kZRange,kZRange);
  TH1F *hTrImpact = new TH1F("hTrImpact","Transverse impact parameter (all tracks)",100,0,30.*kMaxImpact); 

  // Open RunLoader
  AliRunLoader *rl = AliRunLoader::Open(galicefilename);
  if (!rl) {
    ::Error("AliL1Delay.C","Can't start session !");
    return;
  }
  rl->LoadgAlice();
  gAlice = rl->GetAliRun();

  // Set magnetic field
  AliMagF* field = TGeoGlobalMagField::Instance()->GetField();
  AliExternalTrackParam::SetNonuniformFieldTracking();
  const Float_t sfield = field->SolenoidField();

  // GetGeo manager
  if(withTGeo)
    {
      TFile *geoFile = TFile::Open(geomfilename);
      if (!geoFile || !geoFile->IsOpen()) {
	::Error("AliL1Delay.C","Can't open %s !",geomfilename);
	return;
      }
      gGeoManager = (TGeoManager *)geoFile->Get("Geo");
      if(!gGeoManager) {
	::Error("AliL1Delay.C","Can't find TGeoManager !");
	return;
      }
    }

  // Get X positions of the pad-rows
  TDirectory* saveDir = gDirectory;
  rl->CdGAFile();
  AliTPCParam* param = (AliTPCParam*) gDirectory->Get("75x40_100x60_150x60");
  if (!param) {
    ::Error("AliL1Delay.C","No TPC parameters found !");
    return;
  }
  Int_t nRowLow = param->GetNRowLow();
  Int_t nRowUp = param->GetNRowUp();
  Float_t xRow[159];
  for (Int_t i=0;i<nRowLow;i++){
    xRow[i] = param->GetPadRowRadiiLow(i);
  }
  for (Int_t i=0;i<nRowUp;i++){
    xRow[i+nRowLow] = param->GetPadRowRadiiUp(i);
  }
  saveDir->cd();

  // Open file with ESDs
  TFile *esdFile = TFile::Open(esdfilename);
  if (!esdFile || !esdFile->IsOpen()) {
    ::Error("AliL1Delay.C","Can't open %s !",esdfilename);
    return;
  }
  AliESD* event = new AliESD;
  TTree* esdTree = (TTree*) esdFile->Get("esdTree");
  if (!esdTree) {
    ::Error("AliL1Delay.C","No ESD tree found !");
    return;
  }
  esdTree->SetBranchAddress("ESD", &event);

  // Init PID
  AliPID pid;

  Int_t nEvent = esdTree->GetEntries();

  // Loop over events
  for(Int_t iEvent = 0; iEvent<nEvent; iEvent++)
    {
      esdTree->GetEvent(iEvent);
      const AliESDVertex *esdVtx = event->GetVertex();
      Double_t zVtx = esdVtx->GetZv();
      if(zVtx == 0.) continue; // Vertex is not found. Skip the event.
      Int_t nTrk=event->GetNumberOfTracks();

      ::Info("AliL1Delay.C","Event %d : Z vertex = %f  |  Number of tracks = %d",iEvent,zVtx,nTrk);

      for(Int_t iTrk=0; iTrk<nTrk; iTrk++)
	{
	  AliESDtrack *track = event->GetTrack(iTrk);
	  if(!track) continue;

	  // Track selection 
	  if((track->GetStatus()&AliESDtrack::kTPCrefit)==0) continue;
	  if(track->GetTPCclusters((Int_t)0x0)<(kMinNClusters*3)) continue;
	  if(track->GetP()<kPtMin) continue;

	  Int_t firstRow = (Int_t)(track->GetTPCPoints(0)+0.5);
	  Int_t lastRow = (Int_t)(track->GetTPCPoints(2)+0.5);
	  Double_t firstXYZ[3];
	  Double_t lastXYZ[3];
	  track->GetXYZAt(xRow[firstRow], sfield, firstXYZ);
	  track->GetXYZAt(xRow[lastRow], sfield, lastXYZ);

	  // Select if the track is forward or backward one
	  Bool_t IsForward;
	  if(firstXYZ[2] > 0 && lastXYZ[2] > 0)
	    IsForward = kTRUE;
	  else if(firstXYZ[2] < 0 && lastXYZ[2] < 0)
	    IsForward = kFALSE;
	  else
	    continue;

	  // Propagate the track to the vertex
	  Float_t vtxXYZ[3];
	  if(withTGeo) {
	    if(PropagateToVertexG(track,vtxXYZ) == 0) continue;
	  }
	  else {
	    if(PropagateToVertex(track,vtxXYZ) == 0) continue;
	  }

	  // Get the position at point of closest approach to the vertex
	  Float_t impact = TMath::Sqrt(vtxXYZ[0]*vtxXYZ[0]+vtxXYZ[1]*vtxXYZ[1]);
	  hTrImpact->Fill(impact);
	  // Select primary tracks
	  if(impact > kMaxImpact) continue;

	  if (IsForward)
	    hDeltaZForward->Fill(vtxXYZ[2]-zVtx);
	  else
	    hDeltaZBackward->Fill(vtxXYZ[2]-zVtx);
	}
    }

  //  delete event;
  //  esdFile->Close();

  // Check track propagation
  TF1 *fexpo = new TF1("fexpo","expo");
  hTrImpact->Fit("fexpo","E");
  Double_t slope = fexpo->GetParameter(1);
  if(slope > 3.*kMaxImpact)
    ::Warning("AliL1Delay.C","The transverse impact parameter distribution is too broad: %f +- %f",slope,fexpo->GetParError(1));
  delete fexpo;

  // Fit histograms and extract L1 time delay
  Double_t shifts[2],shifterrs2[2];
  {
    Double_t params[3];
    params[0] = hDeltaZForward->GetMaximumBin();
    params[1] = hDeltaZForward->GetBinCenter(hDeltaZForward->GetMaximumBin());
    params[2] = kDeltaZRes;
    TF1 *fgaus = new TF1("fgaus","gaus");
    fgaus->SetParameters(params);
    hDeltaZForward->Fit("fgaus","E");
    shifts[0] = fgaus->GetParameter(1);
    shifterrs2[0] = fgaus->GetParError(1)*fgaus->GetParError(1);
    delete fgaus;
  }
  {
    Double_t params[3];
    params[0] = hDeltaZBackward->GetMaximumBin();
    params[1] = hDeltaZBackward->GetBinCenter(hDeltaZBackward->GetMaximumBin());
    params[2] = kDeltaZRes;
    TF1 *fgaus = new TF1("fgaus","gaus");
    fgaus->SetParameters(params);
    hDeltaZBackward->Fit("fgaus","E");
    shifts[1] = -fgaus->GetParameter(1);
    shifterrs2[1] = fgaus->GetParError(1)*fgaus->GetParError(1);
    delete fgaus;
  }

  // Check the consistency of the two time delays
  if(TMath::Abs(shifts[0]-shifts[1])>3.*TMath::Sqrt(shifterrs2[0]+shifterrs2[1]))
    ::Warning("AliL1Delay.C","The extracted delays are more than 3 sigma away from each other: %f %f",shifts[0],shifts[1]);

  // Calculated the overall time delay
  Double_t shifterr = 1./TMath::Sqrt(1./shifterrs2[0]+1./shifterrs2[1]);
  Double_t shift = (shifts[0]/shifterrs2[0]+shifts[1]/shifterrs2[1])*shifterr*shifterr;

  ::Info("AliL1Delay.C","");
  ::Info("AliL1Delay.C","=====================================================");
  ::Info("AliL1Delay.C","The measured L1 delay is (%f +- %f) cm",shift,shifterr);
  ::Info("AliL1Delay.C","The measured L1 delay is (%f +- %f) ns",shift*1.e9/param->GetDriftV(),shifterr*1.e9/param->GetDriftV());
  ::Info("AliL1Delay.C","=====================================================");

  delete event;
  esdFile->Close();

  gStyle->SetOptFit(1);

  TCanvas *c1 = new TCanvas("c1","",0,0,700,850);

  TPad *p1 = new TPad("p1","",0,0,1,0.5); p1->Draw();
  p1->cd(); p1->SetFillColor(42); p1->SetFrameFillColor(10);
  p1->SetLogy(1);
  hTrImpact->SetFillColor(4);  hTrImpact->SetXTitle("(cm)");
  hTrImpact->SetStats(kTRUE); hTrImpact->Draw(); c1->cd();

  TPad *p2 = new TPad("p2","",0,0.5,0.5,1); p2->Draw();
  p2->cd(); p2->SetFillColor(42); p2->SetFrameFillColor(10); 
  hDeltaZForward->SetFillColor(4);  hDeltaZForward->SetXTitle("(cm)");
  hDeltaZForward->SetStats(kTRUE); hDeltaZForward->Draw(); c1->cd();

  TPad *p3 = new TPad("p3","",0.5,0.5,1,1); p3->Draw();
  p3->cd(); p3->SetFillColor(42); p3->SetFrameFillColor(10); 
  hDeltaZBackward->SetFillColor(4);  hDeltaZBackward->SetXTitle("(cm)");
  hDeltaZBackward->SetStats(kTRUE); hDeltaZBackward->Draw(); c1->cd();

  TFile f("AliL1Delay.root","RECREATE");
  c1->Write();
  f.Close();

  delete rl;
}

Int_t CorrectForDeadZoneMaterial(AliITStrackV2 *t) {
  //--------------------------------------------------------------------
  // Correction for the material between the TPC and the ITS
  // (should it belong to the TPC code ?)
  //--------------------------------------------------------------------
  Double_t riw=80., diw=0.0053, x0iw=30; // TPC inner wall ? 
  Double_t rcd=61., dcd=0.0053, x0cd=30; // TPC "central drum" ?
  Double_t yr=12.8, dr=0.03; // rods ?
  Double_t zm=0.2, dm=0.40;  // membrane
  //Double_t rr=52., dr=0.19, x0r=24., yyr=7.77; //rails
  Double_t rs=50., ds=0.001; // something belonging to the ITS (screen ?)

  if (t->GetX() > riw) {
     if (!t->PropagateTo(riw,diw,x0iw)) return 1;
     if (TMath::Abs(t->GetY())>yr) t->CorrectForMaterial(dr);
     if (TMath::Abs(t->GetZ())<zm) t->CorrectForMaterial(dm);
     if (!t->PropagateTo(rcd,dcd,x0cd)) return 1;
     //Double_t x,y,z; t->GetGlobalXYZat(rr,x,y,z);
     //if (TMath::Abs(y)<yyr) t->PropagateTo(rr,dr,x0r); 
     if (!t->PropagateTo(rs,ds)) return 1;
  } else if (t->GetX() < rs) {
     if (!t->PropagateTo(rs,-ds)) return 1;
     //Double_t x,y,z; t->GetGlobalXYZat(rr,x,y,z);
     //if (TMath::Abs(y)<yyr) t->PropagateTo(rr,-dr,x0r); 
     if (!t->PropagateTo(rcd,-dcd,x0cd)) return 1;
     if (!t->PropagateTo(riw+0.001,-diw,x0iw)) return 1;
  } else {
  ::Error("CorrectForDeadZoneMaterial","track is already in the dead zone !");
    return 1;
  }
  
  return 0;
}

Bool_t PropagateToVertex(AliESDtrack *track, Float_t *vtxXYZ)
{
  // Make an ITS track and propagate it to the vertex
  AliITStrackV2 itstrack(*track);
  if (CorrectForDeadZoneMaterial(&itstrack) != 0) return 0;
  if (!itstrack.PropagateTo(3.,0.0028,65.19)) return 0;
  if (!itstrack.PropagateToVertex()) return 0;

  itstrack.GetXYZ(vtxXYZ);

  return 1;
}

Bool_t PropagateToVertexG(AliESDtrack *track, Float_t *vtxXYZ)
{
  // Make an ITS track and propagate it to the vertex using TGeoManager
  AliITStrackV2 itstrack(*track);
  AliExternalTrackParam etrack(itstrack);
  Double_t r = 3.;
  Double_t xyz0[3],xyz1[3],y,z;
  Double_t param[7];
  Double_t step = 2.;
  for (Double_t x = etrack.X()-step; x > r; x -= step) {
    etrack.GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);
    Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
    etrack.RotateTo(alpha);
    etrack.GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);      
    if (!etrack.GetProlongationAt(x,y,z)) return 0;
    xyz1[0] = x*TMath::Cos(alpha)-y*TMath::Sin(alpha);
    xyz1[1] = x*TMath::Sin(alpha)+y*TMath::Cos(alpha);
    xyz1[2] = z;
    AliTracker::MeanMaterialBudget(xyz0,xyz1,param);
    if (!etrack.PropagateTo(x,param[1],param[0])) return 0;
  }
  etrack.PropagateTo(r,1,0);
  if (!etrack.PropagateToDCA(0,0,param[1],0)) return 0;

  etrack.GetXYZ(vtxXYZ);

  return 1;
}
