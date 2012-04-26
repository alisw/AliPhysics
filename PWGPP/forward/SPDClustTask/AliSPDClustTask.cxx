/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
// Class AliSPDClustTask                                            //
// Analysis task to produce data and MC histos needed for tracklets      //
// dNdEta extraction in multiple bins in one go                          //
// Author:  ruben.shahoyan@cern.ch                                       //
///////////////////////////////////////////////////////////////////////////
/*
  Important parameters to set:
  1) make sure to initialize correct geometry in UserCreateOutputObjects
  2) The cut on signal selection variable (delta, dphi ...) should be decided beforehand
...
*/

#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TList.h"
#include "TObjArray.h"
#include "TGeoGlobalMagField.h"

#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliESDEvent.h"  
#include "AliESDInputHandler.h"
#include "AliESDInputHandlerRP.h"
#include "AliCDBPath.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliMagF.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliRunLoader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliCentrality.h"
#include "../ITS/AliITSRecPoint.h"
#include "../ITS/AliITSgeomTGeo.h"
#include "../ITS/AliITSMultReconstructor.h" 

#include "AliLog.h"

#include "AliPhysicsSelection.h"
#include "AliSPDClustTask.h"
#include "AliITSMultReconstructor.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliESDtrackCuts.h"

ClassImp(AliSPDClustTask)


//________________________________________________________________________
/*//Default constructor
AliSPDClustTask::AliSPDClustTask(const char *name)
  : AliAnalysisTaskSE(name),
*/  
//________________________________________________________________________
AliSPDClustTask::AliSPDClustTask(const char *name) 
  : AliAnalysisTaskSE(name), 
//
  fOutput(0), 
  fUseMC(kFALSE),
  fHistos(0),
//
  fEtaMin(-3.0),
  fEtaMax(3.0),
  fZVertexMin(-20),
  fZVertexMax( 20),
//
  fScaleDTBySin2T(kFALSE),
  fCutOnDThetaX(kFALSE),
  fNStdDev(1.),
  fDPhiWindow(0.08),
  fDThetaWindow(0.025),
  fDPhiShift(0.0045),
  fPhiOverlapCut(0.005),
  fZetaOverlap(0.05),
  fRemoveOverlaps(kFALSE),
//
  fDPhiSCut(0.06),
  fNStdCut(1.),
//
  fMultReco(0),
  fRPTree(0),
  fStack(0),
  fMCEvent(0),
  //
  fDontMerge(kFALSE),
  fhPtPionIn(0),
  fhPtKaonIn(0),
  fhPtProtonIn(0)
{
  // Constructor

  DefineOutput(1, TList::Class());
  //
  SetScaleDThetaBySin2T();
  SetNStdDev();
  SetPhiWindow();
  SetThetaWindow();
  SetPhiShift();
  SetPhiOverlapCut();
  SetZetaOverlapCut();
  SetRemoveOverlaps();
  //
}

//________________________________________________________________________
AliSPDClustTask::~AliSPDClustTask()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {  //RRR
    printf("Deleteing output\n");
    delete fOutput;
    fOutput = 0;
  }
  //
  delete fMultReco;
  delete fHistos;
  if (fhPtPionIn) delete fhPtPionIn;
  if (fhPtKaonIn) delete fhPtKaonIn;
  if (fhPtProtonIn) delete fhPtProtonIn;
  //
}

//________________________________________________________________________
void AliSPDClustTask::UserCreateOutputObjects() 
{
  //
  //  OpenFile(1);  fDontMerge = kTRUE;
  //  

  fOutput = new TList();
  fOutput->SetOwner(); 
  //
  AliCDBManager *man = AliCDBManager::Instance();
  if (fUseMC) {
    Bool_t newGeom = kTRUE;
    man->SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual");
    if (newGeom) {
      // new geom
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",130844,8);
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",130844,6,-1)) AliFatal("Failed to misalign geometry");
    }
    else {
      // old geom
      AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",130845,7);
      AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
      if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",130845,5,-1)) AliFatal("Failed to misalign geometry");
    }
  }
  else {
    man->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB"); //man->SetRun(137366);
    AliCDBEntry*  obj = man->Get("GRP/Geometry/Data",137366);
    AliGeomManager::SetGeometry((TGeoManager*) obj->GetObject());
    if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",137366,-1,-1)) AliFatal("Failed to misalign geometry");
  }
  //
  // Create histograms
  fHistos = BookHistos();
  for(Int_t i=0; i<fHistos->GetEntriesFast(); i++) {
    if(fHistos->At(i)){
      fOutput->AddLast(fHistos->At(i));
    }
  }

  PostData(1, fOutput);
  //
}

//________________________________________________________________________
void AliSPDClustTask::UserExec(Option_t *) 
{
  // Main loop
  //
  //
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  fRPTree = 0;
  AliESDInputHandlerRP *handRP = (AliESDInputHandlerRP*)anMan->GetInputEventHandler(); 
  if (!handRP) { AliError("No RP handler"); return; }
  //
  fRPTree = handRP->GetTreeR("ITS");
  if (!fRPTree) { AliError(" Invalid ITS cluster tree !\n"); return; }
  //
  AliESDEvent *esd  = handRP->GetEvent();
  if (!esd) { AliError("No AliESDEvent"); return; }
  //
  AliMCEventHandler* eventHandler = 0;
  fMCEvent = 0;
  fStack = 0;
  //
  if (fUseMC) {
    eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
    if (!eventHandler) { AliError("Could not retrieve MC event handler"); return; }
    fMCEvent = eventHandler->MCEvent();
    if (!fMCEvent) { AliError("Could not retrieve MC event"); return; }
    fStack = fMCEvent->Stack();
    if (!fStack) { AliError("Stack not available"); return; }
  }
  //
  const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
  if (vtxESD->GetNContributors()<1) return;
  //
  // ------------------RS: Put here additional event selection if needed ------->>
  // at the moment, I am just selecting on vertex
  TString vtxTyp = vtxESD->GetTitle();
  if (vtxTyp.Contains("vertexer: Z")) {
    if (vtxESD->GetDispersion()>0.04) return;
    if (vtxESD->GetZRes()>0.25) return;
  }
  //
  if (vtxESD->GetZ()<fZVertexMin || vtxESD->GetZ()>fZVertexMax) return;
  //
  // pile-up rejection
  if (esd->IsPileupFromSPD(3,0.8)) return;
  //
  // ------------------RS: Put here additional event selection if needed -------<<
  //
  fVtx[0] = vtxESD->GetX();
  fVtx[1] = vtxESD->GetY();
  fVtx[2] = vtxESD->GetZ();
  //
  // do we need to initialize the field?
  AliMagF* field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!field && !esd->InitMagneticField()) {printf("Failed to initialize the B field\n");return;}
  //
  AliStack *stack = NULL;
  if (fUseMC) {
    if (!fMCEvent) { AliError("Could not retrieve MC event"); return; }
    stack = fMCEvent->Stack();
    if (!stack) { AliError("ERROR: Could not get stack"); return; }
    for (Int_t iPart = 0; iPart < stack->GetNtrack(); ++iPart) {
      if (!stack->IsPhysicalPrimary(iPart)) continue;
      TParticle* p = stack->Particle(iPart);
      if (TMath::Abs(p->Y()) > 0.5) continue;
      if (TMath::Abs(p->GetPdgCode()) != 211) continue;
      TH1D *hPt = (TH1D*)fHistos->At(kHPt);
      hPt->Fill(p->Pt());
    }
  }
  //
  InitMultReco();
  fMultReco->Reconstruct(fRPTree, fVtx);
  //  AliMultiplicity* mlt = fMultReco->GetMultiplicity();
  FillHistos(stack);
  //
  delete fMultReco; 
  fMultReco = 0;
  //
}      

//________________________________________________________________________
void AliSPDClustTask::Terminate(Option_t *) 
{
  Printf("Terminating...");
}

//_________________________________________________________________________
void AliSPDClustTask::InitMultReco()
{
  // create mult reconstructor
  if (fMultReco) delete fMultReco;
  fMultReco = new AliITSMultReconstructor();
  fMultReco->SetCreateClustersCopy(kTRUE);
  fMultReco->SetScaleDThetaBySin2T(fScaleDTBySin2T);
  fMultReco->SetNStdDev(fNStdDev);
  fMultReco->SetPhiWindow( fDPhiWindow );
  fMultReco->SetThetaWindow( fDThetaWindow );
  fMultReco->SetPhiShift( fDPhiShift );
  fMultReco->SetRemoveClustersFromOverlaps(fRemoveOverlaps);
  fMultReco->SetPhiOverlapCut(fPhiOverlapCut);
  fMultReco->SetZetaOverlapCut(fZetaOverlap);
  fMultReco->SetHistOn(kFALSE); 
  //
}

//_________________________________________________________________________
TObjArray* AliSPDClustTask::BookHistos()
{
  // book set of histos
  //
  TObjArray* histos = new TObjArray();
  histos->SetOwner(kFALSE);
  // this array contains histos stored at slots:
  // 0-99    : tracklet related histos
  // 100-199 : SPD1 clusters not used by tracklets
  // 200-299 : SPD2 clusters not used by tracklets
  // 300-399 : SPD1 clusters used by tracklets
  // 400-499 : SPD2 clusters used by tracklets
  //
  int nEtaBins = int((fEtaMax-fEtaMin)/0.1);
  if (nEtaBins<1) nEtaBins = 1;
  //
  // just an example: cluster type vs eta
  for (int iused=0;iused<2;iused++) {
    for (int spd=0;spd<2;spd++) {
      TH2F* h = new TH2F(Form("clType_SPD%d_%s",spd,iused ? "used":"free"),
			 Form("clType SPD%d %s",spd,iused ? "used":"free"),
			 nEtaBins, fEtaMin,fEtaMax,
			 15, -0.5, 14.5);
      histos->AddAtAndExpand(h, kHClusters+iused*200+spd*100 + kClTypevsEta);
      TH1F* hZ = new TH1F(Form("clZ_SPD%d_%s",spd,iused ? "used":"free"),
			  Form("clZ SPD%d %s",spd,iused ? "used":"free"),
			  200, -15, 15);
      histos->AddAtAndExpand(hZ, kHClusters+iused*200+spd*100 + kClZ);
      TH1F* hEta = new TH1F(Form("clEta_SPD%d_%s",spd,iused ? "used":"free"),
			    Form("clEta SPD%d %s",spd,iused ? "used":"free"),
			    nEtaBins, fEtaMin,fEtaMax);
      histos->AddAtAndExpand(hEta, kHClusters+iused*200+spd*100 + kClEta);

      hZ = new TH1F(Form("clZ_pions_weighted_SPD%d_%s",spd,iused ? "used":"free"),
		    Form("clZ pions weighted SPD%d %s",spd,iused ? "used":"free"),
		    200, -15, 15);
      histos->AddAtAndExpand(hZ, kHClusters+iused*200+spd*100 + kClZPionsW);
      hEta = new TH1F(Form("clEta_pions_weighted_SPD%d_%s",spd,iused ? "used":"free"),
		      Form("clEta pions weighted SPD%d %s",spd,iused ? "used":"free"),
		      nEtaBins, fEtaMin,fEtaMax);
      histos->AddAtAndExpand(hEta, kHClusters+iused*200+spd*100 + kClEtaPionsW);

      hZ = new TH1F(Form("clZ_pions_SPD%d_%s",spd,iused ? "used":"free"),
		    Form("clZ pions SPD%d %s",spd,iused ? "used":"free"),
		    200, -15, 15);
      histos->AddAtAndExpand(hZ, kHClusters+iused*200+spd*100 + kClZPions);
      hEta = new TH1F(Form("clEta_pions_SPD%d_%s",spd,iused ? "used":"free"),
		      Form("clEta pions SPD%d %s",spd,iused ? "used":"free"),
		      nEtaBins, fEtaMin,fEtaMax);
      histos->AddAtAndExpand(hEta, kHClusters+iused*200+spd*100 + kClEtaPions);
    }
  }
  if (fUseMC) {
    TH1D *hPt = new TH1D("hPt","Pions Pt spectra (MC)",fhPtPionIn->GetXaxis()->GetNbins(),fhPtPionIn->GetXaxis()->GetXbins()->GetArray());
    histos->AddAtAndExpand(hPt, kHPt);
  }
  //
  return histos;
}

//_________________________________________________________________________
void AliSPDClustTask::FillHistos(AliStack *stack)
{
  // fill info on clusters, separately on associated to tracklets and singles
  const int kUsedBit = BIT(15);
  //
  // 0) --------- just in case: reset kUsedBit of clusters --->>>
  for (int ilr=0;ilr<2;ilr++) {
    for (int icl=fMultReco->GetNClustersLayer(ilr);icl--;) {  // loop on clusters of layer
      AliITSRecPoint *clus = fMultReco->GetRecPoint(ilr,icl);
      if (!clus) {
	AliError(Form("Failed to extract cluster %d of layer %d",icl,ilr));
	return;
      }
      clus->ResetBit(kUsedBit);
    } // loop on clusters of layer
  } // loop on layers
  //
  // 1) --------- mark used clusters by kUsedBit ------------>>>
  int ntr = fMultReco->GetNTracklets();
  for (int itr=ntr;itr--;) {
    Float_t *trc = fMultReco->GetTracklet(itr);
    AliITSRecPoint *cl0 = fMultReco->GetRecPoint(0,(int)trc[AliITSMultReconstructor::kClID1]); // cluster on SPD1
    AliITSRecPoint *cl1 = fMultReco->GetRecPoint(1,(int)trc[AliITSMultReconstructor::kClID2]); // cluster on SPD1
    //
    if (!cl0 || !cl1) {
      AliError(Form("Failed to extract clusters associated with tracklet %d",itr));
      continue;
    }
    cl0->SetBit(kUsedBit);
    cl1->SetBit(kUsedBit);
  }
  //
  // 2) --------- fill needed info for used and unuses clusters
  TH2F* h2d = 0;
  TH1F* h1d = 0;
  //
  for (int ilr=0;ilr<2;ilr++) { // loop on layers
    int ncl = fMultReco->GetNClustersLayer(ilr);
    for (int icl=ncl;icl--;) {  // loop on clusters of layer
      AliITSRecPoint *clus = fMultReco->GetRecPoint(ilr,icl); // extract the cluster to access clType
      Float_t* clInfo = fMultReco->GetClusterOfLayer(ilr,icl); // we already have in the MultReco cluster params extracted and
      // stored in the float array, see AliITSMultReconstructor
      //
      float phi = clInfo[AliITSMultReconstructor::kClPh];
      float eta = -TMath::Log(TMath::Tan(clInfo[AliITSMultReconstructor::kClTh]/2));
      float z   = clInfo[AliITSMultReconstructor::kClZ];
      //
      int used = clus->TestBit(kUsedBit) ? 1:0;
      //
      // Annalisa, here you should fill the info on used/unused clusters
      if (phi > 0 && phi < 1.6) {
	h2d = (TH2F*)fHistos->At(kHClusters+used*200+ilr*100 + kClTypevsEta);
	h2d->Fill(eta, clus->GetClusterType());
	h1d = (TH1F*)fHistos->At(kHClusters+used*200+ilr*100 + kClZ);
	h1d->Fill(z);
	h1d = (TH1F*)fHistos->At(kHClusters+used*200+ilr*100 + kClEta);
	h1d->Fill(eta);
	Double_t weight = 1.;
	if (stack) {
	  for(Int_t iLabel = 0; iLabel < 3; ++iLabel) {
	    Int_t clLabel = clus->GetLabel(iLabel);
	    if (clLabel <= 0) continue;
	    Int_t moLabel = FindMotherParticle(stack, clLabel);
	    if (moLabel <= 0) continue;
	    TParticle* p = stack->Particle(moLabel);
	    if ((TMath::Abs(p->GetPdgCode()) != 211) &&
		(TMath::Abs(p->GetPdgCode()) != 321) &&
		(TMath::Abs(p->GetPdgCode()) != 2212)) continue;
	    weight = PtWeight(p->Pt(),p->GetPdgCode());
	    h1d = (TH1F*)fHistos->At(kHClusters+used*200+ilr*100 + kClZPionsW);
	    h1d->Fill(z,weight);
	    h1d = (TH1F*)fHistos->At(kHClusters+used*200+ilr*100 + kClEtaPionsW);
	    h1d->Fill(eta,weight);
	    h1d = (TH1F*)fHistos->At(kHClusters+used*200+ilr*100 + kClZPions);
	    h1d->Fill(z);
	    h1d = (TH1F*)fHistos->At(kHClusters+used*200+ilr*100 + kClEtaPions);
	    h1d->Fill(eta);
	    break;
	  }
	}
      }
    } // loop on clusters of layer
  } // loop on layers
  //
}

//________________________________________________________________________
void AliSPDClustTask::SetInput(const char *filename)
{
  // Read pions Pt spectra
  // from a local input file

  AliInfo(Form("Reading input histograms from %s",filename));
  if (fUseMC) {
    TFile *fPt = TFile::Open("FitOutput.root");
    if (!fPt) {	AliError("Failed to open input file"); return; }
    fhPtPionIn = (TH1F*)fPt->Get("pionPtWeight")->Clone("fhPtPionIn");
    fhPtPionIn->SetDirectory(0);
    fhPtKaonIn = (TH1F*)fPt->Get("kaonPtWeight")->Clone("fhPtKaonIn");
    fhPtKaonIn->SetDirectory(0);
    fhPtProtonIn = (TH1F*)fPt->Get("protPtWeight")->Clone("fhPtProtonIn");
    fhPtProtonIn->SetDirectory(0);
    fPt->Close();
  }
}

//________________________________________________________________________
Int_t AliSPDClustTask::FindMotherParticle(AliStack* stack, Int_t i)
{
  if(stack->IsPhysicalPrimary(i)) return i;
  TParticle* p = stack->Particle(i);
  Int_t imo =  p->GetFirstMother();

  //  printf("imo = %d\n",imo);

  if(imo<=0)
    {
      AliInfo("particle is not a PhysPrim and has no mother");
      return imo;
    }
  return FindMotherParticle(stack, imo);

}

//________________________________________________________________________
Double_t AliSPDClustTask::PtWeight(Double_t pt, Int_t pdgCode)
{
  switch (TMath::Abs(pdgCode)) {
  case 211:
    return fhPtPionIn->GetBinContent(fhPtPionIn->FindBin(pt));
  case 321:
    return fhPtKaonIn->GetBinContent(fhPtKaonIn->FindBin(pt));
  case 2212:
    return fhPtProtonIn->GetBinContent(fhPtProtonIn->FindBin(pt));
  default:
    return 1.;
  }
}
