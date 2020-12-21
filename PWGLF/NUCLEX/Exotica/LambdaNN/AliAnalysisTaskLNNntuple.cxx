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
//                 AliAnalysisTaskLNNntuple class
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
#include "AliAnalysisTaskLNNntuple.h"
#include "AliCentrality.h"
#include "AliTRDPIDResponse.h"
#include "TString.h"
#include <TDatime.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
//#include <AliVTrack.h>

ClassImp (AliAnalysisTaskLNNntuple)
//________________________________________________________________________
AliAnalysisTaskLNNntuple::AliAnalysisTaskLNNntuple ():
AliAnalysisTaskSE (),
  fMaxPtPion(0.55),
  fMinPtTriton(0.75),
  fYear(2015),
  fMC(kTRUE),
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
AliAnalysisTaskLNNntuple::AliAnalysisTaskLNNntuple (const char *name, Bool_t mc):
  AliAnalysisTaskSE (name),
  fMaxPtPion(0.55),
  fMinPtTriton(0.75),
  fYear(2015),
  fMC(mc),
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
AliAnalysisTaskLNNntuple::~AliAnalysisTaskLNNntuple ()
{
  // Destructor
  if (fListHist)
    {
      delete fListHist;
      fListHist = 0;
    }
}

//==================DEFINITION OF OUTPUT OBJECTS==============================

void AliAnalysisTaskLNNntuple::UserCreateOutputObjects ()
{

  fListHist = new TList ();
  fListHist->SetOwner ();	// IMPORTANT!
  if(fYear==2015) fEventCuts.SetupLHC15o();
  if(fYear==2011) fEventCuts.SetupRun1PbPb();

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


  if (!fNt)
    {
      if(!fMC){

	fNt = new TNtupleD ("nt", "V0 ntuple","piPx:piPy:piPz:triPx:triPy:triPz:nSpi:nStri:triTOFmass:piTPCsig:triTPCsig:v0P:ptArm:alphaArm:triDcaXY:triDcaZ:v0DcaD:decayL:decayLxy:v0Dca:CosP:v0VtxErrSum:sign:dcaPi:is3Hele:nSPiFromPiTof:nSPrTof:nSPiTof:nITSclus:triTRDsig");
	fListHist->Add (fNt);
      } else {
	fNt = new TNtupleD ("nt", "V0 ntuple","piPx:piPy:piPz:triPx:triPy:triPz:nSpi:nStri:triTOFmass:piTPCsig:triTPCsig:v0P:ptArm:alphaArm:triDcaXY:triDCAZ:v0DcaD:decayL:decayLxy:v0Dca:CosP:v0VtxErrSum:sign:dcaPi:is3Hele:nSPiFromPiTof:nSPrTof:nSPiTof:nITSclus:piPdgCode:triPdgCode:piMumPdgCode:triMumPdgCode:triTRDsig");
	fListHist->Add (fNt);

      }
    }


  PostData (1, fListHist);
}				// end UserCreateOutputObjects


//====================== USER EXEC ========================

void
AliAnalysisTaskLNNntuple::UserExec (Option_t *)
{
  //------------------------------------------

  // Main loop

  AliVEvent *event = InputEvent ();
  if (!event)
    {
      Printf ("ERROR: Could not retrieve event");
      return;
    }


  if(fMC) Info ("AliAnalysisTaskLNNntuple for MC", "Starting UserExec");
  else Info ("AliAnalysisTaskLNNntuple for Data", "Starting UserExec");

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
    nbMcTracks = stack->GetNtrack();;
    if(!stack) { 
      Printf( "Stack not available, Exiting... \n"); 
      return; 
    }
  }


  // Load ESD event
  AliESDEvent *lESDevent = dynamic_cast < AliESDEvent * >(event);
  if (!lESDevent)
    {
      AliError ("Cannot get the ESD event");
      return;
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
	  Info ("AliAnalysisTaskLNNntuple","No good vertex, skip event");
	  return;		// NO GOOD VERTEX, SKIP EVENT 
	}
    }

  fHistEventMultiplicity->Fill (1);	// analyzed events with PV


  if (TMath::Abs (vtx->GetX()) > 10)  return;

  selCentrality = SelectCentrality(centrality,kTRUE);
  fHistEventMultiplicity->Fill (selCentrality);

  fHistEventMultiplicity->Fill (2);

  // track quality monitor plots
  for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
    AliESDtrack *esdtrack=lESDevent->GetTrack(j);
    if(!esdtrack) {
      AliError(Form("ERROR: Could not retrieve esdtrack %d",j));
      continue;
    }
    // ************** Track cuts ****************
    if (!PassTrackCuts(esdtrack)) continue;

    fhBB->Fill(esdtrack->P()/esdtrack->GetSign(),esdtrack->GetTPCsignal());

    fhTestNsigma->Fill(esdtrack->P(),fPIDResponse->NumberOfSigmasTPC(esdtrack,(AliPID::EParticleType)6));
    Double_t tofmass = -1;
    tofmass = GetTOFmass(esdtrack);
    if(tofmass>0)  {
      fhMassTOF->Fill(esdtrack->P(),tofmass);
    }
  }


  // *************** Loop over V0s ***************
  for(Int_t iv=0; iv<lESDevent->GetNumberOfV0s(); iv++){
    AliESDv0 * v0s = lESDevent->GetV0(iv);
    if(!v0s) continue;
    TVector3 globVtx(vtx->GetX(),vtx->GetY(),vtx->GetZ());
    TVector3 v0Vtx(v0s->Xv(),v0s->Yv(),v0s->Zv());
    TVector3 decayL = globVtx-v0Vtx;
    Double_t path = decayL.Mag();
    if(!Passv0Cuts(v0s,path)) continue;
    AliESDtrack* trackPos = lESDevent->GetTrack(v0s->GetPindex());
    AliESDtrack* trackNeg = lESDevent->GetTrack(v0s->GetNindex());
    if(!PassTrackCuts(trackPos)) continue;
    if(!PassTrackCuts(trackNeg)) continue;

    AliESDtrack *pion =0x0;
    AliESDtrack *triton=0x0;

    Double_t momPi[3], momTri[3];
    if(IsPionCandidate(trackPos) && IsTritonCandidate(trackNeg)) {
      pion = trackPos;
      triton = trackNeg;
      v0s->GetPPxPyPz(momPi[0],momPi[1],momPi[2]);
      v0s->GetNPxPyPz(momTri[0],momTri[1],momTri[2]);
    } else if(IsPionCandidate(trackNeg)&&IsTritonCandidate(trackPos)){
      pion = trackNeg; triton = trackPos;
      v0s->GetNPxPyPz(momPi[0],momPi[1],momPi[2]);
      v0s->GetPPxPyPz(momTri[0],momTri[1],momTri[2]);
    }
    else continue;
    // monitor dE/dx of tracks in selected V0s
    fhBBPions->Fill(pion->P()*pion->GetSign(),pion->GetTPCsignal());
    fhBBH3->Fill(triton->P()*triton->GetSign(),triton->GetTPCsignal());
    Double_t tofMass = GetTOFmass(triton);
    if(tofMass>2.5 && tofMass <3.5)  fhBBH3TofSel->Fill(triton->P()*triton->GetSign(),triton->GetTPCsignal()); 
    Double_t nSPi = fPIDResponse->NumberOfSigmasTPC (pion,(AliPID::EParticleType) 2);
    Double_t nSTri = fPIDResponse->NumberOfSigmasTPC (triton,(AliPID::EParticleType) 6);

    Double_t decayLXY = TMath::Sqrt(decayL.Mag2 () - decayL.Z () * decayL.Z ());

    Double_t ptArm = v0s->PtArmV0 ();
    Double_t alphaArm = v0s->AlphaV0 ();

    Float_t dcaTri[2] = { -100, -100 };
    Double_t Vtxpos[3] ={ globVtx.X (), globVtx.Y (), globVtx.Z () };
    triton->GetDZ (Vtxpos[0], Vtxpos[1], Vtxpos[2], lESDevent->GetMagneticField (), dcaTri);
    //dca to primary vertex  
    Float_t dcaPi= pion->GetD(Vtxpos[0],Vtxpos[1],lESDevent->GetMagneticField());
    //Float_t dcaTriTot = triton->GetD(Vtxpos[0],Vtxpos[1],lESDevent->GetMagneticField()); -> same as dcaTri[1]
    Double_t nSPiFromPiTof = fPIDResponse->NumberOfSigmasTOF (pion,(AliPID::EParticleType) 2);
    Double_t nSPrTof = fPIDResponse->NumberOfSigmasTOF (triton,(AliPID::EParticleType) 4); //check if 3H is identified as proton in TOF PID
    Double_t nSPiTof = fPIDResponse->NumberOfSigmasTOF (triton,(AliPID::EParticleType) 2); // check if 3H is identified as pion TOF PID

    AliESDVertex vtxV0 = v0s->GetVertex ();

    Double_t sigmaVtxV0[2] = { pion->GetD(vtxV0.GetX (), vtxV0.GetY (),lMagneticField), triton->GetD (vtxV0.GetX (), vtxV0.GetY (), lMagneticField) };
    Double_t err = TMath::Sqrt (sigmaVtxV0[0] * sigmaVtxV0[0] + sigmaVtxV0[1] * sigmaVtxV0[1]);

    fhTestQ->Fill(pion->Charge(),triton->Charge());




    Float_t isHele=0;
    Int_t nTrackletsPID=0;
  
    Float_t eleEff[3] = {0.85,0.9,0.95}; 
    AliPIDResponse::EDetPidStatus pidStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD,triton);
    if(pidStatus!=AliPIDResponse::kDetPidOk) isHele = -1;  
    else {
      for(Int_t iEff=0; iEff<3; iEff++)  {
	Bool_t isEle =  fPIDResponse->IdentifiedAsElectronTRD(triton,nTrackletsPID,eleEff[iEff],centrality,AliTRDPIDResponse::kLQ2D);
	if(isEle && nTrackletsPID>3) isHele += (iEff+1)*TMath::Power(10,iEff);
      }
    }

   
    if(fMC){
      Double_t pdgPion=-1, pdgTriton=-1, pdgPionMum=-1, pdgTritonMum=-1;
      if(TMath::Abs(pion->GetLabel()<nbMcTracks)) {
	TParticle *piMC = stack->Particle(TMath::Abs(pion->GetLabel()));
	pdgPion = piMC->GetPdgCode();
	Int_t mum = piMC->GetFirstMother();
	if(mum>0 && mum<nbMcTracks){
	  pdgPionMum = stack->Particle(mum)->GetPdgCode();
	}
      }
      if(TMath::Abs(triton->GetLabel()<nbMcTracks)) {
	TParticle *tritonMC = stack->Particle(TMath::Abs(triton->GetLabel()));
	pdgTriton = tritonMC->GetPdgCode();
	Int_t mum = tritonMC->GetFirstMother();
	if(mum>0 && mum<nbMcTracks){
	  pdgTritonMum = stack->Particle(mum)->GetPdgCode();
	}
      }
      Double_t ntuple[34] =
	{ momPi[0], momPi[1], momPi[2], momTri[0], momTri[1],momTri[2], nSPi, nSTri,
	  tofMass,pion->GetTPCsignal (), triton->GetTPCsignal (),
	  v0s->P(), ptArm, alphaArm, dcaTri[0], dcaTri[1],v0s->GetDcaV0Daughters (), decayL.Mag(), decayLXY,
	  v0s->GetD (vtx->GetX (), vtx->GetY (), vtx->GetZ ()),v0s->GetV0CosineOfPointingAngle (), err, pion->GetSign () + triton->GetSign ()*10,
	  dcaPi,isHele,nSPiFromPiTof,nSPrTof,nSPiTof,pion->GetNumberOfITSClusters()+100.*triton->GetNumberOfITSClusters(),
	  pdgPion,pdgTriton,pdgPionMum,pdgTritonMum,triton->GetTRDsignal()};

      fNt->Fill (ntuple);
    } else {
      Double_t ntuple[30] =
	{ momPi[0], momPi[1], momPi[2], momTri[0], momTri[1],momTri[2], nSPi, nSTri,
	  tofMass,pion->GetTPCsignal (), triton->GetTPCsignal (),
	  v0s->P(), ptArm, alphaArm, dcaTri[0], dcaTri[1], v0s->GetDcaV0Daughters (), decayL.Mag(), decayLXY,
	  v0s->GetD (vtx->GetX (), vtx->GetY (), vtx->GetZ ()),v0s->GetV0CosineOfPointingAngle (), err, pion->GetSign () + triton->GetSign ()*10,
	  dcaPi,isHele,nSPiFromPiTof,nSPrTof,nSPiTof,pion->GetNumberOfITSClusters()+100.*triton->GetNumberOfITSClusters(),triton->GetTRDsignal()};
      fNt->Fill (ntuple);
    }
  } // loop over V0s


  PostData (1, fListHist);

}				//end userexec


//________________________________________________________________________

void
AliAnalysisTaskLNNntuple::Terminate (Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
}

//________________________________________________________________________
Bool_t AliAnalysisTaskLNNntuple::PassTrackCuts (AliESDtrack * tr)
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
Bool_t AliAnalysisTaskLNNntuple::Passv0Cuts (AliESDv0 * v0, Double_t decayLength)
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
Bool_t AliAnalysisTaskLNNntuple::IsPionCandidate (AliESDtrack * tr)
{

  if(tr->Pt()> fMaxPtPion) return kFALSE;
  AliPIDResponse::EDetPidStatus statusTPC;
  Double_t nSigmaTPC = -999;
  statusTPC = fPIDResponse->NumberOfSigmas (AliPIDResponse::kTPC, tr, AliPID::kPion,nSigmaTPC);
  Bool_t z = kFALSE;
  if (statusTPC == AliPIDResponse::kDetPidOk && TMath::Abs (nSigmaTPC) <= 3.5) z = kTRUE;
  return z;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskLNNntuple::IsTritonCandidate(AliESDtrack * tr)
{
  if(tr->Pt()<fMinPtTriton) return kFALSE;
  AliPIDResponse::EDetPidStatus statusTPC;
  Double_t nSigmaTPC = -999;
  statusTPC = fPIDResponse->NumberOfSigmas (AliPIDResponse::kTPC, tr, AliPID::kTriton, nSigmaTPC);
  Bool_t z = kFALSE;
  if (statusTPC == AliPIDResponse::kDetPidOk && TMath::Abs (nSigmaTPC) <= 3.5) z = kTRUE;
  if(z) {
   fTPCclusPID->Fill(tr->GetTPCNcls(),tr->GetTPCsignalN());
    if(tr->GetTPCsignalN()<40) z = kFALSE;
   }
  return z;
}
//________________________________________________________________________
Double_t AliAnalysisTaskLNNntuple::GetTOFmass(AliESDtrack * tr)
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
Double_t AliAnalysisTaskLNNntuple::Chi2perNDF (AliESDtrack * track)
{
  // Calculate chi2 per ndf for track

  Int_t nClustersTPC = track->GetTPCNcls ();
  if (nClustersTPC > 5)  return (track->GetTPCchi2 () / Float_t (nClustersTPC - 5));
  else return (-1.);
}
//________________________________________________________________________
Int_t AliAnalysisTaskLNNntuple::SelectCentrality(Float_t perc, Bool_t isPrimVtx){

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
