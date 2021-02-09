/**************************************************************************
 * Contributors are not mentioned at all.                                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
//-----------------------------------------------------------------
//                 AliAnalysisTaskLNNv0Bkg class
//-----------------------------------------------------------------


class TTree;
class TParticle;
class TVector3;

#include "AliAnalysisManager.h"
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

class AliESDVertex;
class AliESDv0;

#include <iostream>
#include "AliAnalysisTaskSE.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TCutG.h"
#include "TF1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskLNNv0Bkg.h"
#include "AliCentrality.h"
#include "AliTRDPIDResponse.h"
#include "TString.h"
#include <TDatime.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
//#include <AliVTrack.h>

ClassImp (AliAnalysisTaskLNNv0Bkg)
 //________________________________________________________________________
 AliAnalysisTaskLNNv0Bkg::AliAnalysisTaskLNNv0Bkg ():
  AliAnalysisTaskSE (),
  fMaxPtPion(0.55),
  fMinPtTriton(0.75),
  fMC(kFALSE),
  fYear(2015),
  fBkgType(1),
  fEventCuts(),
  fListHist (0),
  fHistEventMultiplicity (0),
  fHistTrackMultiplicity (0),
  fhBB (0),
  fhTOF (0),
  fhMassTOF (0),
  fhBBPions (0),
  fhBBH3 (0),
  fhBBH3TofSel (0), fhTestNsigma (0), fTPCclusPID(0), fhTestQ(0), fNt (0), fPIDResponse(0)
{
 // Dummy Constructor
 //Printf("                   ****************** Default ctor \n");
}

//________________________________________________________________________
AliAnalysisTaskLNNv0Bkg::AliAnalysisTaskLNNv0Bkg (const char *name):
 AliAnalysisTaskSE (name),
 fMaxPtPion(0.55),
 fMinPtTriton(0.75),
 fMC(kFALSE),
 fYear(2015),
 fBkgType(1),
 fEventCuts(),
 fListHist (0),
 fHistEventMultiplicity (0),
 fHistTrackMultiplicity (0),
 fhBB (0),
 fhTOF (0),
 fhMassTOF (0),
 fhBBPions (0),
 fhBBH3 (0),
 fhBBH3TofSel (0),
 fhTestNsigma (0),
 fTPCclusPID(0),
 fhTestQ(0),
 fNt (0), fPIDResponse(0)
{

 // Define input and output slots here
 // Input slot #0 works with a TChain
 //DefineInput(0, TChain::Class());
 // Output slot #0 writes into a TList container ()

 DefineInput (0, TChain::Class ());

 DefineOutput (1, TList::Class ());
 //DefineOutput(2, TTree::Class());
 //DefineOutput(3, TTree::Class());
 //Printf("                   ****************** Right ctor \n");
}

//_______________________________________________________
AliAnalysisTaskLNNv0Bkg::~AliAnalysisTaskLNNv0Bkg ()
{
 // Destructor
 if (fListHist)
 {
  delete fListHist;
  fListHist = 0;
 }
}

//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskLNNv0Bkg::UserCreateOutputObjects ()
{

 fListHist = new TList ();
 fListHist->SetOwner ();	// IMPORTANT!
 if(fYear==2015) fEventCuts.SetupLHC15o();
 //if(fYear==2011) fEventCuts.SetupRun1PbPb();

 if (!fHistEventMultiplicity)
 {
  fHistEventMultiplicity =
   new TH1F ("fHistEventMultiplicity", "Nb of Events", 13, -0.5, 12.5);
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (1, "All Events");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (2, "Events w/PV");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (3, "Events w/|Vz|<10cm");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (4, "Central Events");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (5, "SemiCentral Events");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (6, "MB Events");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (7, "Central Events  w/|Vz|<10cm");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (8, "SemiCentral Events  w/|Vz|<10cm");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (9, "MB Events w/|Vz|<10cm");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (10,"Any Events");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (11,"Any Events w/|Vz|<10cm");
  fHistEventMultiplicity->GetXaxis ()->SetBinLabel (12,"other");

  fListHist->Add (fHistEventMultiplicity);
 }

 if (!fHistTrackMultiplicity)
 {
  fHistTrackMultiplicity =
   new TH2F ("fHistTrackMultiplicity", "Nb of Tracks", 2500, 0, 25000,
     210, -1, 104);
  fHistTrackMultiplicity->GetXaxis ()->SetTitle ("Number of tracks");
  fHistTrackMultiplicity->GetYaxis ()->SetTitle ("Percentile");
  fListHist->Add (fHistTrackMultiplicity);
 }

 Double_t pMax = 8;
 Double_t binWidth = 0.1;
 Int_t nBinBB = (Int_t) (2 * pMax / binWidth);

 if (!fhBB)
 {
  fhBB = new TH2F ("fhBB", "BetheBlochTPC", nBinBB, -pMax, pMax, 400, 0, 800);
  fhBB->GetXaxis ()->SetTitle ("p/z (GeV/#it{c})");
  fhBB->GetYaxis ()->SetTitle ("TPC Signal");
  fListHist->Add (fhBB);
 }

 if (!fhTOF)
 {
  fhTOF = new TH2F ("fhTOF", "Scatter Plot TOF", nBinBB, -pMax, pMax, 100, 0,1.2);
  fhTOF->GetXaxis ()->SetTitle ("p/z (GeV/#it{c})");
  fhTOF->GetYaxis ()->SetTitle ("#beta");
  fListHist->Add (fhTOF);
 }

 if (!fhMassTOF)
 {
  fhMassTOF =
   new TH2F ("fhMassTOF", "Particle Mass - TOF", nBinBB, 0, pMax, 400, 0,5);
  fhMassTOF->GetYaxis ()->SetTitle ("Mass (GeV/#it{c}^{2})");
  fhMassTOF->GetXaxis ()->SetTitle ("P (GeV/#it{c})");
  fListHist->Add (fhMassTOF);
 }

 if (!fhBBPions)
 {
  fhBBPions = new TH2F ("fhBBPions", "Bethe-Bloch TPC Pions", nBinBB, -pMax, pMax, 400, 0, 800);
  fhBBPions->GetXaxis ()->SetTitle ("p/z (GeV/#it{c})");
  fhBBPions->GetYaxis ()->SetTitle ("TPC Signal");
  fListHist->Add (fhBBPions);
 }

 if (!fhBBH3)
 {
  fhBBH3 = new TH2F ("fhBBH3", "Bethe-Bloch TPC ^3H", nBinBB, -pMax, pMax, 400, 0, 800);
  fhBBH3->GetXaxis ()->SetTitle ("p/z (GeV/#it{c})");
  fhBBH3->GetYaxis ()->SetTitle ("TPC Signal");
  fListHist->Add (fhBBH3);
 }
 if (!fhBBH3TofSel)
 {
  fhBBH3TofSel =
   new TH2F ("fhBBH3TofSel","Bethe-Bloch TPC #^{3}H after TOF 3#sigma cut", nBinBB,-pMax, pMax, 400, 0, 800);
  fhBBH3TofSel->GetXaxis ()->SetTitle ("p/z (GeV/#it{c})");
  fhBBH3TofSel->GetYaxis ()->SetTitle ("TPC Signal");
  fListHist->Add (fhBBH3TofSel);
 }

 if (!fhTestNsigma)
 {
  fhTestNsigma = new TH2F ("hNsigmaTri", "n #sigma distribution", 300, 0, 15, 100, -10, 10);
  fListHist->Add (fhTestNsigma);
 }
 if(!fTPCclusPID){
  fTPCclusPID = new TH2F("fTPCclusPID","triton track : clusters vs PID clusters",200,0,200,200,0,200);
  fListHist->Add (fTPCclusPID);

 }

 if(!fhTestQ){
  fhTestQ = new TH2F("htestQ","candidate charge as from TPC",16,-2,2,16,-2,2);
  fhTestQ->SetXTitle("#pi Id charge");
  fhTestQ->SetYTitle("3H Id charge");
  fListHist->Add(fhTestQ);
 }


 if(fMC){
          fNt = new TNtupleD ("nt", "V0 ntuple","piPx:piPy:piPz:triPx:triPy:triPz:nSpi:nStri:triTOFmass:piTPCsig:triTPCsig:v0P:ptArm:alphaArm:triDcaXY:triDCAZ:v0DcaD:decayL:decayLxy:v0Dca:CosP:nSElTof:sign:dcaPi:is3Hele:nSPiFromPiTof:nSPrTof:nSPiTof:nITSclus:triTRDPIDsig:triTRDsig:pdgPi:pdgTri:pdgMumPi:pdgMumTri");
 } else { fNt = new TNtupleD ("nt", "V0 ntuple","piPx:piPy:piPz:triPx:triPy:triPz:nSpi:nStri:triTOFmass:piTPCsig:triTPCsig:v0P:ptArm:alphaArm:triDcaXY:triDCAZ:v0DcaD:decayL:decayLxy:v0Dca:CosP:nSElTof:sign:dcaPi:is3Hele:nSPiFromPiTof:nSPrTof:nSPiTof:nITSclus:triTRDPIDsig:triTRDsig");}
 fListHist->Add (fNt);



 PostData (1, fListHist);
}				// end UserCreateOutputObjects


//====================== USER EXEC ========================

 void
AliAnalysisTaskLNNv0Bkg::UserExec (Option_t *)
{
 //------------------------------------------

 // Main loop

 AliVEvent *event = InputEvent ();
 if (!event)
 {
  Printf ("ERROR: Could not retrieve event");
  return;
 }


 if(fMC) Info ("AliAnalysisTaskLNNv0Bkg for MC", "Starting UserExec");
 else Info ("AliAnalysisTaskLNNv0Bkg for Data", "Starting UserExec");

 // Load ESD event
 AliESDEvent *lESDevent = dynamic_cast < AliESDEvent * >(event);
 if (!lESDevent)
 {
  AliError ("Cannot get the ESD event");
  return;
 }


 // Load MC if required
 AliStack *stack = 0;
 Int_t nbMcTracks=0;
 if(fMC){
  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent) {
   Printf("ERROR: Could not retrieve MC event");
   return;
  }
  stack = mcEvent->Stack();
  nbMcTracks = stack->GetNtrack();
  if(!stack) {
   Printf( "Stack not available, Exiting... \n");
   return;
  }
 }


 fHistEventMultiplicity->Fill (0);

 Double_t lMagneticField = lESDevent->GetMagneticField ();
 Int_t TrackNumber = -1;

 AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager ();
 AliInputEventHandler *inputHandler = (AliInputEventHandler *) (man->GetInputEventHandler ());
 TrackNumber = lESDevent->GetNumberOfTracks ();
 if (TrackNumber < 3)  return;
 if(!fEventCuts.AcceptEvent(lESDevent)) return;
 Float_t centrality=fEventCuts.GetCentrality(0); // 0 = V0M

 fHistTrackMultiplicity->Fill(TrackNumber,centrality);

 Int_t selCentrality = SelectCentrality(centrality);
 fHistEventMultiplicity->Fill (selCentrality);


 //*******//
 //  PID  //
 //*******//
 fPIDResponse = inputHandler->GetPIDResponse ();
 if(!fPIDResponse) {
  AliError("!!! Please check the PID task, it is not loaded!!\n");
  return;
 }
 fPIDResponse->SetCachePID (kTRUE);
 //===========================================

 // Primary vertex cut
 const AliESDVertex *vtx = lESDevent->GetPrimaryVertexTracks ();
 if (vtx->GetNContributors () < 1)
 {

  // SPD vertex cut
  vtx = lESDevent->GetPrimaryVertexSPD ();

  if (vtx->GetNContributors () < 1)
  {
   Info ("AliAnalysisTaskLNNv0Bkg","No good vertex, skip event");
   return;		// NO GOOD VERTEX, SKIP EVENT 
  }
 }

 fHistEventMultiplicity->Fill (1);	// analyzed events with PV


 if (TMath::Abs (vtx->GetZ()) > 10)  return;

 selCentrality = SelectCentrality(centrality,kTRUE);
 fHistEventMultiplicity->Fill (selCentrality);

 fHistEventMultiplicity->Fill (2);

 // **** list the tracks for further V0 finding *****// 
 //Int_t piMax = 1000, triMax = 0;
 Double_t Rmin=0.2;    //min radius of the fiducial volume
 Double_t Rmax=200.;   //max radius of the fiducial volume
 Double_t Dmin=0.05;  //min imp parameter for each daughter


 TArrayI pionId (0), tritonId (0);
 Int_t idPi = 0, idTri = 0;

 for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
  AliESDtrack *esdtrack=lESDevent->GetTrack(j);
  if(!esdtrack) {
   AliError(Form("ERROR: Could not retrieve esdtrack %d",j));
   continue;
  }
  // ************** Track cuts ****************
  if (!PassTrackCuts(esdtrack)) continue;
  //******************//
  // V0 daughter dca must be > 0.05 cm and < 200 cm See AliESDV0vertexer::Tracks2V0vertices
  Double_t d=esdtrack->GetD(vtx->GetX(),vtx->GetY(),lMagneticField);
  if (TMath::Abs(d)<Dmin) continue;
  if (TMath::Abs(d)>Rmax) continue;

  if (IsPionCandidate (esdtrack))
  {
   idPi++;
   //if (idPi > piMax) pionId.Set (idPi);
   pionId.Set (idPi);
   pionId.AddAt (j, idPi - 1);
  }

  if (IsTritonCandidate (esdtrack))
  {
   idTri++;
   //if (idTri > triMax) pionId.Set (idTri);
   tritonId.Set(idTri);
   tritonId.AddAt (j, idTri - 1);
  }

 }


 if(idTri==0) {
  printf("no triton identified, skipping event \n");
  return;
 }
 if(idTri>0) printf("\n\n **** identified tritons %i and pions %i *** \n",idTri,idPi); 
 // *************** Loop over tracks to recompute V0 by rotating pions  ***************
 printf("size Array pion %i and triton %i \n",pionId.GetSize(),tritonId.GetSize());

 TArrayI trdInfoId(idTri); // here the TRD tracklet information corresponding to the triton AliESDtrack is stored

 Int_t nTrdTracks = lESDevent->GetNumberOfTrdTracks();  // one has to loop over the online tracks

 for (Int_t iTrack = 0; iTrack < nTrdTracks; iTrack++) {
  AliVTrdTrack *trdTrack = lESDevent->GetTrdTrack(iTrack);
  if (!trdTrack) {
   AliError(Form("Failed to get track %i", iTrack));
   continue;
  }
  AliVTrack *matchedesdTrack = trdTrack->GetTrackMatch();  // checks if one has a matching offline track
  if(!matchedesdTrack) continue; 
  Int_t trdId = matchedesdTrack->GetID();
  for(Int_t iLoop =0; iLoop < idTri; iLoop++){
   AliESDtrack *trTmp = lESDevent->GetTrack(tritonId.At(iLoop));
   trdInfoId.AddAt(-999,iLoop);
   Int_t esdId = trTmp->GetID();
   if(trdId==esdId){
    trdInfoId.AddAt(trdTrack->GetPID(),iLoop); 
   }

  }
 }

 Int_t nv0 =0;
 Int_t nv0Tot =0;
 for (Int_t itri = 0; itri < idTri; itri++)
 {
  for (Int_t ipi = 0; ipi < idPi; ipi++)
  {

   nv0Tot++;
   Double_t xn, xp;
   Double_t dca = 0.;
   if(pionId.At(ipi)==tritonId.At(itri)) continue;
   AliESDtrack *ppTrack = lESDevent->GetTrack (pionId.At (ipi));
   AliESDtrack *ttTrack = lESDevent->GetTrack (tritonId.At(itri));
   if (TMath::Abs(ppTrack->GetD(vtx->GetX(),vtx->GetY(),lMagneticField))<0.05)
    if (TMath::Abs(ttTrack->GetD(vtx->GetX(),vtx->GetY(),lMagneticField))<0.05) continue;
   Double_t TRDPIDsignal = 10.*(trdInfoId.At(itri)/10.);
   AliExternalTrackParam trackInPi (*ppTrack);
   AliExternalTrackParam trackIn3H (*ttTrack);
   if(fBkgType==0) {
    if( ppTrack->GetSign ()*ttTrack->GetSign() > 0) continue;
   } else if(fBkgType==1) {
    if( ppTrack->GetSign ()*ttTrack->GetSign() < 0) continue;
   } else {
    if( ppTrack->GetSign ()*ttTrack->GetSign() > 0) continue;
    trackInPi.Set(ppTrack->GetX(), ppTrack->GetAlpha() + TMath::Pi(), ppTrack->GetParameter(),ppTrack->GetCovariance());
   }


   dca = trackInPi.GetDCA(ttTrack, lMagneticField, xn, xp);     //!dca (Neg to Pos)
   if (dca > 1) continue; // dca between v0 daughters < 1 cm!
   if ((xn+xp) > 2*Rmax) continue; // see v0 vertexer
   if ((xn+xp) < 2*Rmin) continue;


   AliExternalTrackParam *nt;
   AliExternalTrackParam *pt;
   if(fBkgType==1){
    nt=ppTrack;
    pt=ttTrack;
   } else if(fBkgType==2 || fBkgType==0){
    if(ppTrack->GetSign()>0){
     pt=ppTrack;
     nt=ttTrack;
    } else {
     nt=ppTrack;
     pt=ttTrack;
    }
   }
   Bool_t corrected=kFALSE;
   if ((nt->GetX() > 3.) && (xn < 3.)) {
    //correct for the beam pipe material
    corrected=kTRUE;
   }
   if ((pt->GetX() > 3.) && (xp < 3.)) {
    //correct for the beam pipe material
    corrected=kTRUE;
   }
   if (corrected) {
    dca=nt->GetDCA(pt,lMagneticField,xn,xp);
    if (dca > 1.) continue;
    if ((xn+xp) > 2*Rmax) continue;
    if ((xn+xp) < 2*Rmin) continue;
   }


   trackInPi.PropagateToDCA(&trackIn3H, lMagneticField);
   //  trackIn3H.PropagateTo(xp, lMagneticField);

   AliESDv0 v0tmp (trackInPi, pionId.At(ipi), trackIn3H,tritonId.At (itri));
   if (v0tmp.GetChi2V0 () > 30) continue;
   Double_t x=v0tmp.Xv(), y=v0tmp.Yv();
   Double_t r2=x*x + y*y;
   if (r2 < Rmin*Rmin) continue;
   if (r2 > Rmax*Rmax) continue;


   v0tmp.SetDcaV0Daughters(dca);
   Float_t CosPointingAngle=v0tmp.GetV0CosineOfPointingAngle(vtx->GetX(),vtx->GetY(),vtx->GetZ()); //PointingAngle
   v0tmp.SetV0CosineOfPointingAngle(CosPointingAngle);

   TVector3 globVtx (vtx->GetX(), vtx->GetY(), vtx->GetZ());
   TVector3 v0Vtx (v0tmp.Xv (), v0tmp.Yv (), v0tmp.Zv ());
   TVector3 decayL = globVtx - v0Vtx;
   Double_t path = decayL.Mag ();

   if (!Passv0Cuts(&v0tmp, path)) continue;

   nv0++;

   Double_t momPi[3], momTri[3];
   if (ppTrack->GetSign () > 0)
   {
    v0tmp.GetPPxPyPz (momPi[0], momPi[1], momPi[2]);
    v0tmp.GetNPxPyPz (momTri[0], momTri[1], momTri[2]);
   }
   else
   {
    v0tmp.GetNPxPyPz (momPi[0], momPi[1], momPi[2]);
    v0tmp.GetPPxPyPz (momTri[0], momTri[1], momTri[2]);
   }

   Double_t decayLXY = TMath::Sqrt(decayL.Mag2() - decayL.Z()*decayL.Z());
   Double_t ptArm = v0tmp.PtArmV0 ();
   Double_t alphaArm = v0tmp.AlphaV0();

   Float_t dcaTri[2] = { -100, -100 };
   Double_t Vtxpos[3] ={globVtx.X(), globVtx.Y(), globVtx.Z()};
   ttTrack->GetDZ(Vtxpos[0], Vtxpos[1], Vtxpos[2],lESDevent->GetMagneticField (), dcaTri);
   Float_t dcaPi= trackInPi.GetD(Vtxpos[0],Vtxpos[1],lESDevent->GetMagneticField());
   Double_t nSPi = fPIDResponse->NumberOfSigmasTPC (ppTrack,(AliPID::EParticleType) 2);
   Double_t nSTri = fPIDResponse->NumberOfSigmasTPC (ttTrack,(AliPID::EParticleType) 6);

   AliESDVertex vtxV0 = v0tmp.GetVertex();

   Double_t sigmaVtxV0[2] = { trackInPi.GetD(v0tmp.Xv(), v0tmp.Yv(),lMagneticField), ttTrack->GetD (v0tmp.Xv(),v0tmp.Yv(), lMagneticField) };
   Double_t err = TMath::Sqrt (sigmaVtxV0[0] * sigmaVtxV0[0] + sigmaVtxV0[1] * sigmaVtxV0[1]);


   Float_t isHele=0;
   Int_t nTrackletsPID=0;

   Float_t eleEff[3] = {0.85,0.9,0.95};
   AliPIDResponse::EDetPidStatus pidStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD,ttTrack);
   if(pidStatus!=AliPIDResponse::kDetPidOk) isHele = -1;
   else {
    for(Int_t iEff=0; iEff<3; iEff++)  {
     Bool_t isEle =  fPIDResponse->IdentifiedAsElectronTRD(ttTrack,nTrackletsPID,eleEff[iEff],centrality,AliTRDPIDResponse::kLQ2D);
     if(isEle && nTrackletsPID>3) isHele += (iEff+1)*TMath::Power(10,iEff);
    }
   }

   Double_t nSPiFromPiTof = fPIDResponse->NumberOfSigmasTOF (ppTrack,(AliPID::EParticleType) 2);
   Double_t nSPrTof = fPIDResponse->NumberOfSigmasTOF (ttTrack,(AliPID::EParticleType) 4); //check if 3H is identified as proton in TOF PID
   Double_t nSPiTof = fPIDResponse->NumberOfSigmasTOF (ttTrack,(AliPID::EParticleType) 2); // check if 3H is identified as pion TOF PID
   Double_t nSElTof = fPIDResponse->NumberOfSigmasTOF (ttTrack,(AliPID::EParticleType) 0); // check if 3H is identified as pion TOF PID


   if(fMC){
    Double_t pdgPi = -1 ;
    Double_t pdgTri = -1 ;
    Double_t pdgMumPi = -1 ;
    Double_t pdgMumTri = -1 ;
    TParticle *part;
    TParticle *partMum;
    if(TMath::Abs(ppTrack->GetLabel())<nbMcTracks){
     part = stack->Particle(TMath::Abs(ppTrack->GetLabel()));
     pdgPi = part->GetPdgCode();
     if(part->GetMother(0)>-1){
      partMum = stack->Particle(part->GetMother(0));
      pdgMumPi = partMum->GetPdgCode();
     }
    }

    if(TMath::Abs(ttTrack->GetLabel())<nbMcTracks){
     part = stack->Particle(TMath::Abs(ttTrack->GetLabel()));
     pdgTri = part->GetPdgCode();
     if(part->GetMother(0)>-1){
      partMum = stack->Particle(part->GetMother(0));
      pdgMumTri = partMum->GetPdgCode();
     }
    }

    Double_t ntupleMC[35] =
    {momPi[0], momPi[1], momPi[2], momTri[0], momTri[1], momTri[2], // 6 
     nSPi, nSTri, GetTOFmass (ttTrack),ppTrack->GetTPCsignal (), ttTrack->GetTPCsignal (),v0tmp.P (), //12
     ptArm, alphaArm, dcaTri[0], dcaTri[1], v0tmp.GetDcaV0Daughters(), decayL.Mag(), decayLXY, // 19
     v0tmp.GetD(vtx->GetX(), vtx->GetY(), vtx->GetZ()), CosPointingAngle,nSElTof,ppTrack->GetSign () + ttTrack->GetSign ()*10, //23
     dcaPi, isHele, nSPiFromPiTof, nSPrTof, nSPiTof, ppTrack->GetNumberOfITSClusters()+100.*ttTrack->GetNumberOfITSClusters(),//29
     TRDPIDsignal,ttTrack->GetTRDsignal(),//31
     pdgPi,pdgTri,pdgMumPi,pdgMumTri};//35

    printf(" ***** MC ntupleFilling with : TRD info %i or %f **** \n",trdInfoId.At(itri),TRDPIDsignal);
    fNt->Fill (ntupleMC);

   } else {

    Double_t ntuple[31] =
    { momPi[0], momPi[1], momPi[2], momTri[0], momTri[1], momTri[2], // 6 
     nSPi, nSTri, GetTOFmass (ttTrack),ppTrack->GetTPCsignal (), ttTrack->GetTPCsignal (),v0tmp.P (), //12
     ptArm, alphaArm, dcaTri[0], dcaTri[1], v0tmp.GetDcaV0Daughters(), decayL.Mag(), decayLXY, // 19
     v0tmp.GetD(vtx->GetX(), vtx->GetY(), vtx->GetZ()), CosPointingAngle,nSElTof,ppTrack->GetSign () + ttTrack->GetSign ()*10, //22
     dcaPi, isHele, nSPiFromPiTof, nSPrTof, nSPiTof, ppTrack->GetNumberOfITSClusters()+100.*ttTrack->GetNumberOfITSClusters(),TRDPIDsignal,ttTrack->GetTRDsignal()};
   fNt->Fill (ntuple);
   }
   printf(" ***** ntupleFilling with : TRD info %i or %f **** \n",trdInfoId.At(itri),TRDPIDsignal);

  }

 }
 printf("total number of V0 combinations  %i -> selected V0s %i \n",nv0Tot,nv0);

 PostData (1, fListHist);

}				//end userexec


//________________________________________________________________________

 void
AliAnalysisTaskLNNv0Bkg::Terminate (Option_t *)
{
 // Draw result to the screen
 // Called once at the end of the query
}

//________________________________________________________________________
Bool_t AliAnalysisTaskLNNv0Bkg::PassTrackCuts (AliESDtrack * tr)
{
 if(tr->GetTPCNcls() < 60 ) return kFALSE;
 if (!(tr->GetStatus () & AliESDtrack::kTPCrefit)) return kFALSE;
 if (Chi2perNDF (tr) > 5) return kFALSE;
 if (tr->GetKinkIndex (0) != 0) return kFALSE;
 if (TMath::Abs (tr->Eta ()) > 0.9) return kFALSE;
 if (tr->P () < 0.15) return kFALSE;
 if (tr->P () > 10) return kFALSE;
 return true;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskLNNv0Bkg::Passv0Cuts (AliESDv0 * v0, Double_t decayLength)
{
 if(v0->GetOnFlyStatus()==kTRUE) return false;
 if (v0->P () < 0.7) return kFALSE;
 if (v0->P () > 10) return kFALSE;
 if (v0->GetDcaV0Daughters () > 0.8) return kFALSE;
 if (v0->GetV0CosineOfPointingAngle () < 0.9995) return kFALSE;
 if (decayLength < 0.5) return false; // loose cut to reduce bkg of V0s coming from primary vertex( weak decay and the coice of 1 cm comes from MC study )

 return true;
}

//________________________________________________________________
Bool_t AliAnalysisTaskLNNv0Bkg::IsPionCandidate (AliESDtrack * tr)
{

 if(tr->Pt()> fMaxPtPion) return kFALSE;
 AliPIDResponse::EDetPidStatus statusTPC;
 Double_t nSigmaTPC = -999;
 statusTPC = fPIDResponse->NumberOfSigmas (AliPIDResponse::kTPC, tr, AliPID::kPion,nSigmaTPC);
 Bool_t z = kFALSE;
 if (statusTPC == AliPIDResponse::kDetPidOk && TMath::Abs (nSigmaTPC) <= 3.0) z = kTRUE;
 return z;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskLNNv0Bkg::IsTritonCandidate(AliESDtrack * tr)
{
 if(tr->Pt()<fMinPtTriton) return kFALSE;
 if(tr->P()>5) return kFALSE;
 AliPIDResponse::EDetPidStatus statusTPC;
 Double_t nSigmaTPC = -999;
 statusTPC = fPIDResponse->NumberOfSigmas (AliPIDResponse::kTPC, tr, AliPID::kTriton, nSigmaTPC);
 Bool_t z = kFALSE;
 if (statusTPC == AliPIDResponse::kDetPidOk && TMath::Abs (nSigmaTPC) <= 3.0) z = kTRUE;
 if(z) {
  fTPCclusPID->Fill(tr->GetTPCNcls(),tr->GetTPCsignalN());
  if(tr->GetTPCsignalN()<40) z = kFALSE;
 }
 return z;
}
//________________________________________________________________________
Double_t AliAnalysisTaskLNNv0Bkg::GetTOFmass(AliESDtrack * tr)
{
 // get TOF mass
 Double_t m = -999;		//mass
 Double_t b = -999;		//beta
 Double_t trackLeng = 0;
 if ((tr->GetStatus () & AliESDtrack::kTOFout) == AliESDtrack::kTOFout)
 {
  trackLeng = tr->GetIntegratedLength ();
  //Double_t p = tr->P ();
  Double_t p = tr->GetTPCmomentum();
  Double_t speedOfLight = TMath::C () * 1E2 * 1E-12;	// cm/ps
  Double_t timeTOF = tr->GetTOFsignal () - fPIDResponse->GetTOFResponse ().GetStartTime (p);	// ps
  b = trackLeng / (timeTOF * speedOfLight);
  if (b >= 1 || b == 0) return m;
  m = p * TMath::Sqrt (1 / (b * b) - 1);
 }
 return m;
}
//________________________________________________________________________
Double_t AliAnalysisTaskLNNv0Bkg::Chi2perNDF (AliESDtrack * track)
{
 // Calculate chi2 per ndf for track

 Int_t nClustersTPC = track->GetTPCNcls ();
 if (nClustersTPC > 5)  return (track->GetTPCchi2 () / Float_t (nClustersTPC - 5));
 else return (-1.);
}
//________________________________________________________________________
Int_t AliAnalysisTaskLNNv0Bkg::SelectCentrality(Float_t perc, Bool_t isPrimVtx){

 //  enum {kCentral=3, kSemiCentral=4, kPeripheral=5, kOther=9, kCentralPV=6,kSemiCentralPV=7, kPeripheralPV=8, kOtherPV=10 };
 Int_t result=-1;
 if(!isPrimVtx){
  if(perc<=10) result=3;
  else if(perc>10 && perc <=40) result=4;
  else if(perc>40 && perc <=90) result=5;
  else result=9;
 } else {
  if(perc<=10) result=6;
  else if(perc>10 && perc <=40) result=7;
  else if(perc>40 && perc <=90) result=8;
  else result=10;
 }

 return result;
}
