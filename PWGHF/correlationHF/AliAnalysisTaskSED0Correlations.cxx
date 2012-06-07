/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for D0 candidates (2Prongs)
// and hadrons correlations
//
// Authors: Chiara Bianchin, chiara.bianchin@pd.infn.it (invariant mass)
// Fabio Colamaria, fabio.colamaria@ba.infn.it (correlations)
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSED0Correlations.h"
#include "AliNormalizationCounter.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSED0Correlations)


//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations():
AliAnalysisTaskSE(),
  fNPtBinsCorr(0), 
  fBinLimsCorr(),
  fPtThreshLow(),
  fPtThreshUp(), 
  fEvents(0),
  fAlreadyFilled(kFALSE),
  fOutputMass(0),
  fOutputCorr(0),
  fOutputStudy(0),
  fNentries(0), 
  fCutsD0(0),
  fCutsTracks(0),
  fReadMC(0),
  fCounter(0),
  fNPtBins(1),
  fFillOnlyD0D0bar(0),
  fIsSelectedCandidate(0),
  fSys(0),
  fIsRejectSDDClusters(0),
  fFillGlobal(kTRUE)
{
  // Default constructor

}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations(const char *name,AliRDHFCutsD0toKpi* cutsD0, AliHFAssociatedTrackCuts* cutsTrk):
  AliAnalysisTaskSE(name),
  fNPtBinsCorr(0),  
  fBinLimsCorr(),
  fPtThreshLow(),
  fPtThreshUp(),
  fEvents(0),
  fAlreadyFilled(kFALSE),
  fOutputMass(0),
  fOutputCorr(0),
  fOutputStudy(0),
  fNentries(0),
  fCutsD0(0),
  fCutsTracks(cutsTrk),
  fReadMC(0),
  fCounter(0),
  fNPtBins(1),
  fFillOnlyD0D0bar(0),
  fIsSelectedCandidate(0),
  fSys(0),
  fIsRejectSDDClusters(0),
  fFillGlobal(kTRUE)
{
  // Default constructor

  fNPtBins=cutsD0->GetNPtBins();
    
  fCutsD0=cutsD0;

  // Output slot #1 writes into a TList container (mass with cuts)
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TH1F container (number of events)
  DefineOutput(2,TH1F::Class());  //My private output
  // Output slot #3 writes into a AliRDHFD0toKpi container (cuts)
  DefineOutput(3,AliRDHFCutsD0toKpi::Class());  //My private output
  // Output slot #4 writes Normalization Counter 
  DefineOutput(4,AliNormalizationCounter::Class());
  // Output slot #5 writes into a TList container (correl output)
  DefineOutput(5,TList::Class());  //My private output
  // Output slot #6 writes into a TList container (correl advanced)
  DefineOutput(6,TList::Class());  //My private output
  // Output slot #7 writes into a AliHFAssociatedTrackCuts container (cuts)
  DefineOutput(7,AliHFAssociatedTrackCuts::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::AliAnalysisTaskSED0Correlations(const AliAnalysisTaskSED0Correlations &source):
  AliAnalysisTaskSE(source),
  fNPtBinsCorr(source.fNPtBinsCorr), 
  fBinLimsCorr(source.fBinLimsCorr),
  fPtThreshLow(source.fPtThreshLow),
  fPtThreshUp(source.fPtThreshUp),
  fEvents(source.fEvents),
  fAlreadyFilled(source.fAlreadyFilled),
  fOutputMass(source.fOutputMass),
  fOutputCorr(source.fOutputCorr),
  fOutputStudy(source.fOutputStudy),
  fNentries(source.fNentries), 
  fCutsD0(source.fCutsD0),
  fCutsTracks(source.fCutsTracks),
  fReadMC(source.fReadMC),
  fCounter(source.fCounter),
  fNPtBins(source.fNPtBins),
  fFillOnlyD0D0bar(source.fFillOnlyD0D0bar),
  fIsSelectedCandidate(source.fIsSelectedCandidate),
  fSys(source.fSys),
  fIsRejectSDDClusters(source.fIsRejectSDDClusters),
  fFillGlobal(source.fFillGlobal)
{
  // Copy constructor
}

//________________________________________________________________________
AliAnalysisTaskSED0Correlations::~AliAnalysisTaskSED0Correlations()
{
  if (fOutputMass) {
    delete fOutputMass;
    fOutputMass = 0;
  }
  if (fOutputCorr) {
    delete fOutputCorr;
    fOutputCorr = 0;
  }
  if (fOutputStudy) {
    delete fOutputStudy;
    fOutputStudy = 0;
  }
  if (fCutsD0) {
    delete fCutsD0;
    fCutsD0 = 0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
  if(fCounter){
    delete fCounter;
    fCounter=0;
  }
}  

//______________________________________________________________________________
AliAnalysisTaskSED0Correlations& AliAnalysisTaskSED0Correlations::operator=(const AliAnalysisTaskSED0Correlations& orig)
{
// Assignment
  if (&orig == this) return *this; //if address is the same (same object), returns itself

  AliAnalysisTaskSE::operator=(orig); //Uses the AliAnalysisTaskSE operator to assign the inherited part of the class
  fNPtBinsCorr = orig.fNPtBinsCorr; 
  fBinLimsCorr = orig.fBinLimsCorr;
  fPtThreshLow = orig.fPtThreshLow;
  fPtThreshUp = orig.fPtThreshUp;
  fEvents = orig.fEvents;
  fAlreadyFilled = orig.fAlreadyFilled;
  fOutputMass = orig.fOutputMass;
  fOutputCorr = orig.fOutputCorr;
  fOutputStudy = orig.fOutputStudy;
  fNentries = orig.fNentries; 
  fCutsD0 = orig.fCutsD0;
  fCutsTracks = orig.fCutsTracks;
  fReadMC = orig.fReadMC;
  fCounter = orig.fCounter;
  fNPtBins = orig.fNPtBins;
  fFillOnlyD0D0bar = orig.fFillOnlyD0D0bar;
  fIsSelectedCandidate = orig.fIsSelectedCandidate;
  fSys = orig.fSys;
  fIsRejectSDDClusters = orig.fIsRejectSDDClusters;
  fFillGlobal = orig.fFillGlobal;

  return *this; //returns pointer of the class
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Correlations::Init() \n");
  
  //Copy of cuts objects
  AliRDHFCutsD0toKpi* copyfCutsD0 = new AliRDHFCutsD0toKpi(*fCutsD0);
  const char* nameoutput=GetOutputSlot(3)->GetContainer()->GetName();
  copyfCutsD0->SetName(nameoutput);

  // Post the data
  PostData(3,copyfCutsD0);
  PostData(7,fCutsTracks);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::UserCreateOutputObjects()
{

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Correlations::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("listMass");

  fOutputCorr = new TList();
  fOutputCorr->SetOwner();
  fOutputCorr->SetName("correlationslist");

  fOutputStudy = new TList();
  fOutputStudy->SetOwner();
  fOutputStudy->SetName("MCstudyplots");

  TString nameMass=" ",nameSgn=" ", nameBkg=" ", nameRfl=" ";

  for(Int_t i=0;i<fCutsD0->GetNPtBins();i++){

    nameMass="histMass_";
    nameMass+=i;
    nameSgn="histSgn_";
    nameSgn+=i;
    nameBkg="histBkg_";
    nameBkg+=i;
    nameRfl="histRfl_";
    nameRfl+=i;

    //histograms of invariant mass distributions

    //MC signal
    if(fReadMC){
      TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",120,1.5648,2.1648);
      tmpSt->Sumw2();

      //Reflection: histo filled with D0Mass which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
      TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries",120,1.5648,2.1648);
      TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",120,1.5648,2.1648);
      tmpBt->Sumw2();
      tmpRt->Sumw2();
      fOutputMass->Add(tmpSt);
      fOutputMass->Add(tmpRt);
      fOutputMass->Add(tmpBt);
    }

    //mass
    TH1F* tmpMt = new TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",120,1.5648,2.1648);
    tmpMt->Sumw2();
    fOutputMass->Add(tmpMt);
  }

  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();

  fNentries=new TH1F(nameoutput, "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of D0 selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 18,-0.5,17.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  fReadMC ? fNentries->GetXaxis()->SetBinLabel(3,"nD0Selected") : fNentries->GetXaxis()->SetBinLabel(3,"Dstar<-D0");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
  if(fSys==0) fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
  if(fReadMC && fSys==0){
    fNentries->GetXaxis()->SetBinLabel(12,"K");
    fNentries->GetXaxis()->SetBinLabel(13,"Lambda");
  }
  fNentries->GetXaxis()->SetBinLabel(14,"Pile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(15,"N. of 0SMH");
  if(fSys==1) fNentries->GetXaxis()->SetBinLabel(16,"Nev in centr");
  if(fIsRejectSDDClusters) fNentries->GetXaxis()->SetBinLabel(17,"SDD-Cls Rej");
  fNentries->GetXaxis()->SetBinLabel(18,"Phys.Sel.Rej");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(4)->GetContainer()->GetName()));
  fCounter->Init();

  CreateCorrelationsObjs(); //creates histos for correlations analysis

  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fNentries);
  PostData(4,fCounter);
  PostData(5,fOutputCorr);
  PostData(6,fOutputStudy);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  //cout<<"I'm in UserExec"<<endl;


  //cuts order
  //     printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
  //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
  //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
  //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
  //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
  //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
  //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
  //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
  //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
  

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  fEvents++;

  TString bname="D0toKpi";

  TClonesArray *inputArray=0;

  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
    }
  } else if(aod) {
    inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  }

  if(!inputArray || !aod) {
    printf("AliAnalysisTaskSED0Correlations::UserExec: input branch not found!\n");
    return;
  }

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSED0Correlations::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSED0Correlations::UserExec: MC header branch not found!\n");
      return;
    }
  }
  
  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fCutsD0,fReadMC); 

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);

  if(!fCutsD0->IsEventSelected(aod)) {
    if(fCutsD0->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
    if(fSys==1 && (fCutsD0->GetWhyRejection()==2 || fCutsD0->GetWhyRejection()==3)) fNentries->Fill(15);
    if(fCutsD0->GetWhyRejection()==7) fNentries->Fill(17);
    return;
  }

  // Check the Nb of SDD clusters
  if (fIsRejectSDDClusters) { 
    Bool_t skipEvent = kFALSE;
    Int_t ntracks = 0;
    if (aod) ntracks = aod->GetNTracks();
    for(Int_t itrack=0; itrack<ntracks; itrack++) { // loop on tacks
      //    ... get the track
      AliAODTrack * track = aod->GetTrack(itrack);
      if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
	skipEvent=kTRUE;
	fNentries->Fill(16);
	break;
      }
    }
    if (skipEvent) return;
  }
  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  Bool_t isGoodVtx=kFALSE;

  //vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fNentries->Fill(3);
  }

  //Reset flag for tracks distributions fill
  fAlreadyFilled=kFALSE;

  // loop over candidates
  Int_t nInD0toKpi = inputArray->GetEntriesFast();
  if(fDebug>2) printf("Number of D0->Kpi: %d\n",nInD0toKpi);

  if(fFillGlobal) { //loop on V0 and tracks for each event, to fill Pt distr. and InvMass distr.

    TClonesArray *v0array = (TClonesArray*)aod->GetList()->FindObject("v0s");
    Int_t pdgCodes[2] = {211,211};
    Int_t idArrayV0[v0array->GetEntriesFast()][2];
    for(int iV0=0; iV0<v0array->GetEntriesFast(); iV0++) {
      for(int j=0; j<2; j++) {idArrayV0[iV0][j]=-2;}
      AliAODv0 *v0 = (AliAODv0*)v0array->UncheckedAt(iV0);
      if(SelectV0(v0,vtx1,2,idArrayV0)) { //option 2 = for mass inv plots only
        if(fReadMC && (v0->MatchToMC(311,mcArray,2,pdgCodes)<0)) continue;
        ((TH2F*)fOutputStudy->FindObject("hK0MassInv"))->Fill(v0->MassK0Short(),v0->Pt()); //invariant mass plot
        ((TH1F*)fOutputStudy->FindObject("hist_Pt_K_AllEv"))->Fill(v0->Pt()); //pT distribution (in all events), K0 case
      }
    }

    for(Int_t itrack=0; itrack<aod->GetNTracks(); itrack++) { // loop on tacks
      AliAODTrack * track = aod->GetTrack(itrack);
      //rejection of tracks
      if(track->GetID() < 0) continue; //discard negative ID tracks
      if(track->Pt() < fPtThreshLow.at(0) || track->Pt() > fPtThreshUp.at(0)) continue; //discard tracks outside pt range for hadrons/K
      if(!fCutsTracks->IsHadronSelected(track,vtx1,aod->GetMagneticField())) continue;
      //pT distribution (in all events), charg and hadr cases
      ((TH1F*)fOutputStudy->FindObject("hist_Pt_Charg_AllEv"))->Fill(track->Pt()); 
      if(fCutsTracks->CheckKaonCompatibility(track,kFALSE,0,2)) ((TH1F*)fOutputStudy->FindObject("hist_Pt_Kcharg_AllEv"))->Fill(track->Pt());
    }
 
  } //end of loops for global plot fill

  Int_t nSelectedloose=0,nSelectedtight=0;  
  for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iD0toKpi);
 
    if(d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
	fNentries->Fill(2);
	continue; //skip the D0 from Dstar
      }
    
    if (fCutsD0->IsInFiducialAcceptance(d->Pt(),d->Y(421)) ) {
      nSelectedloose++;
      nSelectedtight++;      
      if(fSys==0){
	if(fCutsD0->IsSelected(d,AliRDHFCuts::kTracks,aod))fNentries->Fill(6);       
      }
      Int_t ptbin=fCutsD0->PtBin(d->Pt());
      if(ptbin==-1) {fNentries->Fill(4); continue;} //out of bounds

      fIsSelectedCandidate=fCutsD0->IsSelected(d,AliRDHFCuts::kAll,aod); //selected
      if(!fIsSelectedCandidate) continue;

      if(!fReadMC) CalculateCorrelations(aod,d); //correlations on real data
      else { //correlations on MC -> association of selected D0 to MCinfo with MCtruth
        Int_t pdgDgD0toKpi[2]={321,211};
	Int_t labD0 = d->MatchToMC(421,mcArray,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not
        if (labD0>-1) CalculateCorrelations(aod,d,labD0,mcArray);
      }
      
      FillMassHists(d,mcArray,fCutsD0,fOutputMass);
    }

  } //end for prongs
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);  
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);  

  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fNentries);
  PostData(4,fCounter);
  PostData(5,fOutputCorr);
  PostData(6,fOutputStudy);

  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0Correlations::FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi* cuts, TList *listout){
  //
  // function used in UserExec to fill mass histograms:
  //
  if (!fIsSelectedCandidate) return;

  if(fDebug>2)  cout<<"Candidate selected"<<endl;

  Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
  Int_t ptbin = cuts->PtBin(part->Pt());
  
  TString fillthis="";
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t labD0=-1;
  if (fReadMC) labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

  //count candidates selected by cuts
  fNentries->Fill(1);
  //count true D0 selected by cuts
  if (fReadMC && labD0>=0) fNentries->Fill(2);

  if ((fIsSelectedCandidate==1 || fIsSelectedCandidate==3) && fFillOnlyD0D0bar<2) { //D0

    if(fReadMC){
      if(labD0>=0) {
	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	if (pdgD0==421){ //D0
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	} else{ //it was a D0bar
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	}
      } else {//background
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      }
    }else{
      fillthis="histMass_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
    }
     
  }
  if (fIsSelectedCandidate>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)) { //D0bar

    if(fReadMC){
      if(labD0>=0) {
	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();

	if (pdgD0==-421){ //D0bar
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	} else{
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	}
      } else {//background or LS
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
      }
    }else{
      fillthis="histMass_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(invmassD0bar);
    }
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Correlations: Terminate() \n");

  fOutputMass = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputMass) {     
    printf("ERROR: fOutputMass not available\n");
    return;
  }

  fNentries = dynamic_cast<TH1F*>(GetOutputData(2));
  
  if(!fNentries){
    printf("ERROR: fNEntries not available\n");
    return;
  }

  fCutsD0 = dynamic_cast<AliRDHFCutsD0toKpi*>(GetOutputData(3));
  if(!fCutsD0){
    printf("ERROR: fCuts not available\n");
    return;
  }

  fCounter = dynamic_cast<AliNormalizationCounter*>(GetOutputData(4));    
  if (!fCounter) {
    printf("ERROR: fCounter not available\n");
    return;
  }
  fOutputCorr = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutputCorr) {     
    printf("ERROR: fOutputCorr not available\n");
    return;
  }
  fOutputStudy = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputStudy) {     
    printf("ERROR: fOutputStudy not available\n");
    return;
  }
  fCutsTracks = dynamic_cast<AliHFAssociatedTrackCuts*>(GetOutputData(7));
  if(!fCutsTracks){
    printf("ERROR: fCutsTracks not available\n");
    return;
  }

  return;
}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSED0Correlations::CheckD0Origin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
  printf(" AliAnalysisTaskSED0Correlations::CheckD0Origin() \n");
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;

  while (mother > 0){
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if(isQuarkFound) {
    if(isFromB) return 5;
    else return 4;
  }
  else return 1;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::CreateCorrelationsObjs() {
//

  TString namePlot = "";

  //These for limits in THnSparse (one per bin, same limits). 
  //Vars: DeltaPhi, InvMass, PtTrack
  Int_t nBinsPhi[4] = {32,150,30,3};
  Double_t binMinPhi[4] = {-1.6,1.6,0.,0.};  //is the minimum for all the bins
  Double_t binMaxPhi[4] = {4.8,2.2,3.0,3.};  //is the maximum for all the bins

  for(Int_t i=0;i<fNPtBinsCorr;i++){

    //THnSparse plots: correlations for various invariant mass (MC and data)
    namePlot="hPhi_K_Bin";
    namePlot+=i;

    THnSparseI *hPhiK = new THnSparseI(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
    hPhiK->Sumw2();
    fOutputCorr->Add(hPhiK);

    namePlot="hPhi_Kcharg_Bin";
    namePlot+=i;

    THnSparseI *hPhiH = new THnSparseI(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
    hPhiH->Sumw2();
    fOutputCorr->Add(hPhiH);

    namePlot="hPhi_Charg_Bin";
    namePlot+=i;

    THnSparseI *hPhiC = new THnSparseI(namePlot.Data(), "Azimuthal correlation; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
    hPhiC->Sumw2();
    fOutputCorr->Add(hPhiC);

    //histos for c/b origin for D0 (MC only)
    if (fReadMC) {

      //generic origin for tracks
      namePlot="hPhi_K_From_c_Bin";
      namePlot+=i;

      THnSparseI *hPhiK_c = new THnSparseI(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiK_c->Sumw2();
      fOutputCorr->Add(hPhiK_c);

      namePlot="hPhi_Kcharg_From_c_Bin";
      namePlot+=i;

      THnSparseI *hPhiH_c = new THnSparseI(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiH_c->Sumw2();
      fOutputCorr->Add(hPhiH_c);

      namePlot="hPhi_Charg_From_c_Bin";
      namePlot+=i;

      THnSparseI *hPhiC_c = new THnSparseI(namePlot.Data(), "Azimuthal correlation - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiC_c->Sumw2();
      fOutputCorr->Add(hPhiC_c);
  
      namePlot="hPhi_K_From_b_Bin";
      namePlot+=i;

      THnSparseI *hPhiK_b = new THnSparseI(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiK_b->Sumw2();
      fOutputCorr->Add(hPhiK_b);

      namePlot="hPhi_Kcharg_From_b_Bin";
      namePlot+=i;

      THnSparseI *hPhiH_b = new THnSparseI(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiH_b->Sumw2();
      fOutputCorr->Add(hPhiH_b);

      namePlot="hPhi_Charg_From_b_Bin";
      namePlot+=i;

      THnSparseI *hPhiC_b = new THnSparseI(namePlot.Data(), "Azimuthal correlation - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiC_b->Sumw2();
      fOutputCorr->Add(hPhiC_b);

      //HF-only tracks (c for c->D0, b for b->D0)
      namePlot="hPhi_K_HF_From_c_Bin";
      namePlot+=i;

      THnSparseI *hPhiK_HF_c = new THnSparseI(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiK_HF_c->Sumw2();
      fOutputCorr->Add(hPhiK_HF_c);

      namePlot="hPhi_Kcharg_HF_From_c_Bin";
      namePlot+=i;

      THnSparseI *hPhiH_HF_c = new THnSparseI(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiH_HF_c->Sumw2();
      fOutputCorr->Add(hPhiH_HF_c);

      namePlot="hPhi_Charg_HF_From_c_Bin";
      namePlot+=i;
      THnSparseI *hPhiC_HF_c = new THnSparseI(namePlot.Data(), "Azimuthal correlation HF - c origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiC_HF_c->Sumw2();
      fOutputCorr->Add(hPhiC_HF_c);

      namePlot="hPhi_K_HF_From_b_Bin";
      namePlot+=i;

      THnSparseI *hPhiK_HF_b = new THnSparseI(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiK_HF_b->Sumw2();
      fOutputCorr->Add(hPhiK_HF_b);

      namePlot="hPhi_Kcharg_HF_From_b_Bin";
      namePlot+=i;

      THnSparseI *hPhiH_HF_b = new THnSparseI(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiH_HF_b->Sumw2();
      fOutputCorr->Add(hPhiH_HF_b);

      namePlot="hPhi_Charg_HF_From_b_Bin";
      namePlot+=i;

      THnSparseI *hPhiC_HF_b = new THnSparseI(namePlot.Data(), "Azimuthal correlation HF - b origin; #Delta#phi; Inv. Mass (GeV/c^{2}); p_{t} (GeV/c)",4,nBinsPhi,binMinPhi,binMaxPhi);
      hPhiC_HF_b->Sumw2();
      fOutputCorr->Add(hPhiC_HF_b);
    }

    //leading hadron correlations
    namePlot="hPhi_Lead_Bin";
    namePlot+=i;

    TH2F *hCorrLead = new TH2F(namePlot.Data(), "Leading particle correlation; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
    hCorrLead->Sumw2();
    fOutputCorr->Add(hCorrLead);

    if (fReadMC) {
      namePlot="hPhi_Lead_From_c_Bin";
      namePlot+=i;

      TH2F *hCorrLead_c = new TH2F(namePlot.Data(), "Leading particle correlation - c origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      hCorrLead_c->Sumw2();
      fOutputCorr->Add(hCorrLead_c);

      namePlot="hPhi_Lead_From_b_Bin";
      namePlot+=i;

      TH2F *hCorrLead_b = new TH2F(namePlot.Data(), "Leading particle correlation - b origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      hCorrLead_b->Sumw2();
      fOutputCorr->Add(hCorrLead_b);

      namePlot="hPhi_Lead_HF_From_c_Bin";
      namePlot+=i;

      TH2F *hCorrLead_HF_c = new TH2F(namePlot.Data(), "Leading particle correlation HF - c origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      hCorrLead_HF_c->Sumw2();
      fOutputCorr->Add(hCorrLead_HF_c);

      namePlot="hPhi_Lead_HF_From_b_Bin";
      namePlot+=i;

      TH2F *hCorrLead_HF_b = new TH2F(namePlot.Data(), "Leading particle correlation HF - b origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      hCorrLead_HF_b->Sumw2();
      fOutputCorr->Add(hCorrLead_HF_b);
    }
    
    //pT weighted correlations
    namePlot="hPhi_Weig_Bin";
    namePlot+=i;

    TH2F *hCorrWeig = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted); #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
    fOutputCorr->Add(hCorrWeig);

    if (fReadMC) {
      namePlot="hPhi_Weig_From_c_Bin";
      namePlot+=i;

      TH2F *hCorrWeig_c = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted) - c origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeig_c);

      namePlot="hPhi_Weig_From_b_Bin";
      namePlot+=i;

      TH2F *hCorrWeig_b = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted) - b origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeig_b);

      namePlot="hPhi_Weig_HF_From_c_Bin";
      namePlot+=i;

      TH2F *hCorrWeig_HF_c = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted) HF - c origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeig_HF_c);

      namePlot="hPhi_Weig_HF_From_b_Bin";
      namePlot+=i;

      TH2F *hCorrWeig_HF_b = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted) HF - b origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeig_HF_b);
    }

    //pT weighted correlations (squared weights)
    namePlot="hPhi_WeigSqr_Bin";
    namePlot+=i;

    TH2F *hCorrWeigSqr = new TH2F(namePlot.Data(), "Charged particle correlation (pT Weighted - squared); #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
    fOutputCorr->Add(hCorrWeigSqr);

    if (fReadMC) {
      namePlot="hPhi_WeigSqr_From_c_Bin";
      namePlot+=i;

      TH2F *hCorrWeigSqr_c = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - squared) - c origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeigSqr_c);

      namePlot="hPhi_WeigSqr_From_b_Bin";
      namePlot+=i;

      TH2F *hCorrWeigSqr_b = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - squared) - b origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeigSqr_b);

      namePlot="hPhi_WeigSqr_HF_From_c_Bin";
      namePlot+=i;

      TH2F *hCorrWeigSqr_HF_c = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - squared) HF - c origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeigSqr_HF_c);

      namePlot="hPhi_WeigSqr_HF_From_b_Bin";
      namePlot+=i;

      TH2F *hCorrWeigSqr_HF_b = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - squared) HF - b origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeigSqr_HF_b);
    }

    //pT weighted correlations (inverse of pT distribution weights)
    namePlot="hPhi_WeigDist_Bin";
    namePlot+=i;

    TH2F *hCorrWeigDist = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - inv.dis.); #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
    fOutputCorr->Add(hCorrWeigDist);

    if (fReadMC) {
      namePlot="hPhi_WeigDist_From_c_Bin";
      namePlot+=i;

      TH2F *hCorrWeigDist_c = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - inv.dis.) - c origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeigDist_c);

      namePlot="hPhi_WeigDist_From_b_Bin";
      namePlot+=i;

      TH2F *hCorrWeigDist_b = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - inv.dis.) - b origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeigDist_b);

      namePlot="hPhi_WeigDist_HF_From_c_Bin";
      namePlot+=i;

      TH2F *hCorrWeigDist_HF_c = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - inv.dis.) HF - c origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeigDist_HF_c);

      namePlot="hPhi_WeigDist_HF_From_b_Bin";
      namePlot+=i;

      TH2F *hCorrWeigDist_HF_b = new TH2F(namePlot.Data(), "Charged particle correlation (pT weighted - inv.dis.) HF - b origin; #Delta#phi",32,-1.6,4.8,300,1.6,2.2);
      fOutputCorr->Add(hCorrWeigDist_HF_b);
    }

    //counters
    namePlot = "hist_Count_Charg_Bin"; namePlot+=i;
    TH1F *hCountC = new TH1F(namePlot.Data(), "Charged track counter; # Tracks",100,0.,100.);
    hCountC->SetMinimum(0);
    fOutputStudy->Add(hCountC);

    namePlot = "hist_Count_Kcharg_Bin"; namePlot+=i;
    TH1F *hCountH = new TH1F(namePlot.Data(), "Hadrons counter; # Tracks",100,0.,100.);
    hCountH->SetMinimum(0);
    fOutputStudy->Add(hCountH);

    namePlot = "hist_Count_K_Bin"; namePlot+=i;
    TH1F *hCountK = new TH1F(namePlot.Data(), "Kaons counter; # Tracks",10,0.,10.);
    hCountK->SetMinimum(0);
    fOutputStudy->Add(hCountK);

    //pT distribution histos
    namePlot = "hist_Pt_Charg_Bin"; namePlot+=i;
    TH1F *hPtC = new TH1F(namePlot.Data(), "Charged track pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
    hPtC->SetMinimum(0);
    fOutputStudy->Add(hPtC);

    namePlot = "hist_Pt_Kcharg_Bin"; namePlot+=i;
    TH1F *hPtH = new TH1F(namePlot.Data(), "Hadrons pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
    hPtH->SetMinimum(0);
    fOutputStudy->Add(hPtH);

    namePlot = "hist_Pt_K_Bin"; namePlot+=i;
    TH1F *hPtK = new TH1F(namePlot.Data(), "Kaons pT (in D0 evs); p_{T} (GeV/c)",240,0.,12.);
    hPtK->SetMinimum(0);
    fOutputStudy->Add(hPtK);

    //D* feeddown pions rejection histos
    namePlot = "hDstarPions_Bin"; namePlot+=i;
    TH1F *hDstarPions = new TH1F(namePlot.Data(), "Tracks rejected for D* inv.mass cut; # Tracks",2,0.,2.);
    hDstarPions->GetXaxis()->SetBinLabel(1,"Not rejected");
    hDstarPions->GetXaxis()->SetBinLabel(2,"Rejected");
    hDstarPions->SetMinimum(0);
    fOutputStudy->Add(hDstarPions);
  }

  //out of bin loop
  TH1F *hRejTracks = new TH1F("hRejTracks", "Tracks accepted/rejected; # Tracks",2,0.,2.);
  hRejTracks->SetMinimum(0);
  hRejTracks->GetXaxis()->SetBinLabel(1,"Accepted");
  hRejTracks->GetXaxis()->SetBinLabel(2,"Rejected");
  fOutputStudy->Add(hRejTracks);

  if (fFillGlobal) { //all-events plots
    //pt distributions
    TH1F *hPtCAll = new TH1F("hist_Pt_Charg_AllEv", "Charged track pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtCAll->SetMinimum(0);
    fOutputStudy->Add(hPtCAll);

    TH1F *hPtHAll = new TH1F("hist_Pt_Kcharg_AllEv", "Hadrons pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtHAll->SetMinimum(0);
    fOutputStudy->Add(hPtHAll);

    TH1F *hPtKAll = new TH1F("hist_Pt_K_AllEv", "Kaons  pT (All); p_{T} (GeV/c)",240,0.,12.);
    hPtKAll->SetMinimum(0);
    fOutputStudy->Add(hPtKAll);

    //K0 Invariant Mass plots
    TH2F *hK0MassInv = new TH2F("hK0MassInv", "K0 invariant mass; Invariant mass (MeV/c^{2}); pT (GeV/c)",200,0.4,0.6,100,0.,10.);
    hK0MassInv->SetMinimum(0);
    fOutputStudy->Add(hK0MassInv);
  }

  //for MC analysis only
  if (fReadMC) {

    for(Int_t i=0;i<fNPtBinsCorr;i++){

      //displacement histos
      namePlot="histDispl_K_Bin"; namePlot+=i;
      TH1F *hDisplK = new TH1F(namePlot.Data(), "Kaons Displacement; DCA",150,0.,0.15);
      hDisplK->SetMinimum(0);
      fOutputStudy->Add(hDisplK);
  
      namePlot="histDispl_K_HF_Bin";  namePlot+=i;
      TH1F *hDisplK_HF = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplK_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF);

      namePlot="histDispl_Kcharg_Bin"; namePlot+=i;
      TH1F *hDisplHadr = new TH1F(namePlot.Data(), "Hadrons Displacement; DCA",150,0.,0.15);
      hDisplHadr->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr);
  
      namePlot="histDispl_Kcharg_HF_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplHadr_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF);

      namePlot="histDispl_Charg_Bin"; namePlot+=i;
      TH1F *hDisplCharg = new TH1F(namePlot.Data(), "Charged tracks Displacement; DCA",150,0.,0.15);
      hDisplCharg->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg);
  
      namePlot="histDispl_Charg_HF_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only); DCA",150,0.,0.15);
      hDisplCharg_HF->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF);

      namePlot="histDispl_K_From_c_Bin"; namePlot+=i;
      TH1F *hDisplK_c = new TH1F(namePlot.Data(), "Kaons Displacement - c origin; DCA",150,0.,0.15);
      hDisplK_c->SetMinimum(0);
      fOutputStudy->Add(hDisplK_c);
  
      namePlot="histDispl_K_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplK_HF_c = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplK_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF_c);

      namePlot="histDispl_Kcharg_From_c_Bin"; namePlot+=i;
      TH1F *hDisplHadr_c = new TH1F(namePlot.Data(), "Hadrons Displacement - c origin; DCA",150,0.,0.15);
      hDisplHadr_c->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_c);
  
      namePlot="histDispl_Kcharg_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF_c = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplHadr_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF_c);

      namePlot="histDispl_Charg_From_c_Bin"; namePlot+=i;
      TH1F *hDisplCharg_c = new TH1F(namePlot.Data(), "Charged tracks Displacement - c origin; DCA",150,0.,0.15);
      hDisplCharg_c->Sumw2();
      hDisplCharg_c->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_c);
  
      namePlot="histDispl_Charg_HF_From_c_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF_c = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only) - c origin; DCA",150,0.,0.15);
      hDisplCharg_HF_c->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF_c);

      namePlot="histDispl_K_From_b_Bin"; namePlot+=i;
      TH1F *hDisplK_b = new TH1F(namePlot.Data(), "Kaons Displacement - b origin; DCA",150,0.,0.15);
      hDisplK_b->SetMinimum(0);
      fOutputStudy->Add(hDisplK_b);
  
      namePlot="histDispl_K_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplK_HF_b = new TH1F(namePlot.Data(), "Kaons Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplK_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplK_HF_b);

      namePlot="histDispl_Kcharg_From_b_Bin"; namePlot+=i;
      TH1F *hDisplHadr_b = new TH1F(namePlot.Data(), "Hadrons Displacement - b origin; DCA",150,0.,0.15);
      hDisplHadr_b->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_b);

      namePlot="histDispl_Kcharg_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplHadr_HF_b = new TH1F(namePlot.Data(), "Hadrons Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplHadr_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplHadr_HF_b);

      namePlot="histDispl_Charg_From_b_Bin"; namePlot+=i;
      TH1F *hDisplCharg_b = new TH1F(namePlot.Data(), "Charged tracks Displacement - b origin; DCA",150,0.,0.15);
      hDisplCharg_b->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_b);
  
      namePlot="histDispl_Charg_HF_From_b_Bin";  namePlot+=i;
      TH1F *hDisplCharg_HF_b = new TH1F(namePlot.Data(), "Charged tracks Displacement (from HF decay only) - b origin; DCA",150,0.,0.15);
      hDisplCharg_HF_b->SetMinimum(0);
      fOutputStudy->Add(hDisplCharg_HF_b);

      //origin of tracks histos
      namePlot="histOrig_Charg_Bin";  namePlot+=i;
      TH1F *hOrigin_Charm = new TH1F(namePlot.Data(), "Origin of charged tracks",9,0.,9.);
      hOrigin_Charm->SetMinimum(0);
      hOrigin_Charm->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_Charm->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(4,"B->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(5,"B->X-># (X!=D)");
      hOrigin_Charm->GetXaxis()->SetBinLabel(6,"B->D->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(7,"B->D->X->#");
      hOrigin_Charm->GetXaxis()->SetBinLabel(8,"c hadr.");
      hOrigin_Charm->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_Charm);

      namePlot="histOrig_Kcharg_Bin";  namePlot+=i;
      TH1F *hOrigin_Kcharg = new TH1F(namePlot.Data(), "Origin of hadrons",9,0.,9.);
      hOrigin_Kcharg->SetMinimum(0);
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(4,"B->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(5,"B->X-># (X!=D)");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(6,"B->D->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(7,"B->D->X->#");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(8,"c hadr.");
      hOrigin_Kcharg->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_Kcharg);

      namePlot="histOrig_K_Bin";  namePlot+=i;
      TH1F *hOrigin_K = new TH1F(namePlot.Data(), "Origin of kaons",9,0.,9.);
      hOrigin_K->SetMinimum(0);
      hOrigin_K->GetXaxis()->SetBinLabel(1,"Not HF");
      hOrigin_K->GetXaxis()->SetBinLabel(2,"D->#");
      hOrigin_K->GetXaxis()->SetBinLabel(3,"D->X->#");
      hOrigin_K->GetXaxis()->SetBinLabel(4,"B->#");
      hOrigin_K->GetXaxis()->SetBinLabel(5,"B->X-># (X!=D)");
      hOrigin_K->GetXaxis()->SetBinLabel(6,"B->D->#");
      hOrigin_K->GetXaxis()->SetBinLabel(7,"B->D->X->#");
      hOrigin_K->GetXaxis()->SetBinLabel(8,"c hadr.");
      hOrigin_K->GetXaxis()->SetBinLabel(9,"b hadr.");
      fOutputStudy->Add(hOrigin_K);

      //origin of D0 histos
      namePlot="histOrig_D0_Bin";  namePlot+=i;
      TH1F *hOrigin_D0 = new TH1F(namePlot.Data(), "Origin of D0",2,0.,2.);
      hOrigin_D0->SetMinimum(0);
      hOrigin_D0->GetXaxis()->SetBinLabel(1,"From c");
      hOrigin_D0->GetXaxis()->SetBinLabel(2,"From b");
      fOutputStudy->Add(hOrigin_D0);
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::CalculateCorrelations(AliAODEvent* aod, AliAODRecoDecayHF2Prong* d, Int_t labD0, TClonesArray* mcArray) {
//
// Method for correlations D0-hadrons study
//

  Double_t mD0, mD0bar, deltaphi, d0, d0err;
  d->InvMassD0(mD0,mD0bar);
  Int_t ptbin = PtBinCorr(d->Pt());
  if(ptbin < 0) return;
  AliAODVertex *vtx = (AliAODVertex*)aod->GetPrimaryVertex();

  Int_t N_Charg = 0, N_Kcharg = 0, N_Kaons = 0;

  if(fReadMC == 0) {
    Int_t idDaughs[2] = {((AliVTrack*)d->GetDaughter(0))->GetID(),((AliVTrack*)d->GetDaughter(1))->GetID()}; //IDs of daughters to be skipped
    Double_t highPt = 0; Double_t lead[2] = {0,0};  //infos for leading particle (pt,deltaphi)

    for(Int_t itrack=0; itrack<aod->GetNTracks(); itrack++) { // loop on tracks
      AliAODTrack * track = aod->GetTrack(itrack);
     
      if(!TrackSelectedInLoop(track,d,aod,ptbin,idDaughs)) continue; //rejection of track (daughter of D0/quality cuts not passed/soft pion/negative ID

      GetImpParameter(track,aod,d0,d0err); //evaluates d0 and sigma_d0

      ((TH1F*)fOutputStudy->FindObject("hRejTracks"))->Fill(0); //track accepted by quality cuts
      deltaphi = d->Phi()-track->Phi();
      if (deltaphi < -1.571) deltaphi+=6.283;
      if (deltaphi > 4.712) deltaphi-=6.283;
      Double_t pttrack = track->Pt();

      if (pttrack > highPt) {highPt = pttrack; lead[0] = pttrack; lead[1] = deltaphi;}  //for leading particle

      //Lines needed to include overflow into THnSparse projections!
      Double_t ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->GetAxis(2)->GetXmax(); //all plots have same axes...
      Double_t displLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->GetAxis(3)->GetXmax();
      if(pttrack > ptLim_Sparse) pttrack = ptLim_Sparse-0.01;
      if(d0/d0err > displLim_Sparse) d0 = (displLim_Sparse-0.001)*d0err;

      //variables for filling histos
      Double_t fillSpPhiD0[4] = {deltaphi,mD0,pttrack,d0/d0err};
      Double_t fillSpPhiD0bar[4] = {deltaphi,mD0bar,pttrack,d0/d0err};
   
      //generic charged tracks (NO PID selection)
      if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) { //D0
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->Fill(fillSpPhiD0);
      }
      if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) { //D0bar
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->Fill(fillSpPhiD0bar);
      }
      if(!fAlreadyFilled) ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_Charg_Bin%d",ptbin)))->Fill(track->Pt());
      N_Charg++;
        //pT_weighted track correlations fill
        if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) { //D0
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(deltaphi,mD0,pttrack);
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr_Bin%d",ptbin)))->Fill(deltaphi,mD0,pow(pttrack,2.));
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist_Bin%d",ptbin)))->Fill(deltaphi,mD0,PtWeig(ptbin,pttrack));
        }
        if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) { //D0bar
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(deltaphi,mD0bar,pttrack);
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr_Bin%d",ptbin)))->Fill(deltaphi,mD0bar,pow(pttrack,2.));
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist_Bin%d",ptbin)))->Fill(deltaphi,mD0bar,PtWeig(ptbin,pttrack));
        }       

      //hadron identification
      if(fCutsTracks->CheckKaonCompatibility(track,kFALSE,0,2)) {
        if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) { //D0
           ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg_Bin%d",ptbin)))->Fill(fillSpPhiD0);
        }
        if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) { //D0bar
           ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg_Bin%d",ptbin)))->Fill(fillSpPhiD0bar);
        }
        if(!fAlreadyFilled) ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_Kcharg_Bin%d",ptbin)))->Fill(track->Pt());
        N_Kcharg++;
      } //end hadron case

    } //end tracks loop

    //kaon identification
    TClonesArray *v0array = (TClonesArray*)aod->GetList()->FindObject("v0s");
    Int_t idArrayV0[v0array->GetEntriesFast()][2]; //array for skipping K0 double-counting
    for(int iV0=0; iV0<v0array->GetEntriesFast(); iV0++) {
      for(int j=0; j<2; j++) {idArrayV0[iV0][j]=-2;}

      AliAODv0 *v0 = (AliAODv0*)v0array->UncheckedAt(iV0);
      if(v0->Pt() < fPtThreshLow.at(ptbin) || v0->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
      if(SelectV0(v0,vtx,1,idArrayV0)) continue; //option 1 = for correlations

      Double_t ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_Bin%d",ptbin)))->GetAxis(2)->GetXmax(); //all plots have same axes...
      Double_t ptV0=v0->Pt(), deltaphiV0=d->Phi()-v0->Phi();
      if (ptV0 > ptLim_Sparse) ptV0 = ptLim_Sparse-0.01;
      deltaphiV0 = d->Phi()-v0->Phi();
      if (deltaphiV0 < -1.571) deltaphiV0+=6.283;
      if (deltaphiV0 > 4.712) deltaphiV0-=6.283;

      Double_t fillSpPhiD0K0[4] = {deltaphiV0,mD0,ptV0,0.};
      Double_t fillSpPhiD0barK0[4] = {deltaphiV0,mD0bar,ptV0,0.};

      if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) { //D0
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_Bin%d",ptbin)))->Fill(fillSpPhiD0K0);
      }
      if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) { //D0bar
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_Bin%d",ptbin)))->Fill(fillSpPhiD0barK0);
      }
      if(!fAlreadyFilled) ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_K_Bin%d",ptbin)))->Fill(v0->Pt());
      N_Kaons++;
    } // end kaon case

    //Leading track correlations fill
    if(fIsSelectedCandidate == 1 || fIsSelectedCandidate == 3) { //D0
     ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(lead[1],mD0);
    }
    if(fIsSelectedCandidate == 2 || fIsSelectedCandidate == 3) { //D0bar
     ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(lead[1],mD0bar);
    }

  } //end fReadMC = 0

  if(fReadMC == 1) {
    Int_t idDaughs[2] = {((AliVTrack*)d->GetDaughter(0))->GetID(),((AliVTrack*)d->GetDaughter(1))->GetID()}; //IDs of daughters to be skipped
    Double_t highPt = 0; Double_t lead[3] = {0,0,0};  //infos for leading particle (pt,deltaphi,origtrack)

    //Origin of D0
    TString orig=""; Int_t origD0=CheckD0Origin(mcArray,(AliAODMCParticle*)mcArray->At(labD0));
    switch (CheckD0Origin(mcArray,(AliAODMCParticle*)mcArray->At(labD0)))
    {
      case 4:
        orig = "_From_c";
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(0.);
        break;
      case 5:
        orig = "_From_b";
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_D0_Bin%d",ptbin)))->Fill(1.);
        break;
      default:
        return;
    }

    for(Int_t itrack=0; itrack<aod->GetNTracks(); itrack++) { // loop on tracks
      AliAODTrack * track = aod->GetTrack(itrack);

      if(!TrackSelectedInLoop(track,d,aod,ptbin,idDaughs)) continue; //rejection of track (daughter of D0/quality cuts not passed/soft pion/negative ID

      Int_t label = track->GetLabel();
      if (label<0) continue; //discard negative label tracks
      Int_t PDGtrack = ((AliAODMCParticle*)mcArray->At(label))->GetPdgCode();

      GetImpParameter(track,aod,d0,d0err); //evaluates d0 and sigma_d0

      //Infos on track (origin, phi, eta)
      Int_t origTr = CheckTrackOrigin(mcArray,(AliAODMCParticle*)mcArray->At(label));
      deltaphi = d->Phi()-track->Phi();
      if (deltaphi < -1.571) deltaphi+=6.283;
      if (deltaphi > 4.712) deltaphi-=6.283;
      Double_t pttrack = track->Pt();

      if (pttrack > highPt) {highPt = pttrack; lead[0] = pttrack; lead[1] = deltaphi; lead[2] = origTr;}  //for leading particle

      //Lines needed to include overflow into THnSparse projections!
      Double_t d0orig = d0;
      Double_t ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->GetAxis(2)->GetXmax(); //all plots have same axes...
      Double_t displLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->GetAxis(3)->GetXmax();
      if(pttrack > ptLim_Sparse) pttrack = ptLim_Sparse-0.01;
      if(d0/d0err > displLim_Sparse) d0 = (displLim_Sparse-0.001)*d0err;

      //variables for filling histos
      Double_t fillSpPhiD0[4] = {deltaphi,mD0,pttrack,d0/d0err};
      Double_t fillSpPhiD0bar[4] = {deltaphi,mD0bar,pttrack,d0/d0err};

      //generic charged tracks case
      if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==421) { //D0 (from MCTruth)
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->Fill(fillSpPhiD0);
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0);
         if(origD0==4&&origTr>=1&&origTr<=2) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0);//c tr
         if(origD0==5&&origTr>=3&&origTr<=6) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0);//b tr
      }
      if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==-421) { //D0bar
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_Bin%d",ptbin)))->Fill(fillSpPhiD0bar);
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0bar);
         if(origD0==4&&origTr>=1&&origTr<=2) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0bar);
         if(origD0==5&&origTr>=3&&origTr<=6) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Charg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0bar); 
      }
      if(!fAlreadyFilled) {
	((TH1F*)fOutputStudy->FindObject(Form("histDispl_Charg_Bin%d",ptbin)))->Fill(d0orig); //Fills displacement histos
        if (origTr>=1&&origTr<=6) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_Charg_HF_Bin%d",ptbin)))->Fill(d0orig);
        if (origTr>=1&&origTr<=6) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_Charg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(d0orig);
	((TH1F*)fOutputStudy->FindObject(Form("histDispl_Charg%s_Bin%d",orig.Data(),ptbin)))->Fill(d0orig); //Fills displacement histos
        ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_Charg_Bin%d",ptbin)))->Fill(track->Pt());
        ((TH1F*)fOutputStudy->FindObject(Form("histOrig_Charg_Bin%d",ptbin)))->Fill(origTr);
      }
      N_Charg++;
        //pT_weighted track correlations fill
        if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==421) { //D0
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(deltaphi,mD0,pttrack);
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,pttrack);
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr_Bin%d",ptbin)))->Fill(deltaphi,mD0,pow(pttrack,2.));
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,pow(pttrack,2.));
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist_Bin%d",ptbin)))->Fill(deltaphi,mD0,PtWeig(ptbin,pttrack));
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,PtWeig(ptbin,pttrack));
           if(origD0==4&&origTr>=1&&origTr<=2) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,pttrack);//c tr
           if(origD0==5&&origTr>=3&&origTr<=6) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,pttrack);//b tr
           if(origD0==4&&origTr>=1&&origTr<=2) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,pow(pttrack,2.));//c tr
           if(origD0==5&&origTr>=3&&origTr<=6) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,pow(pttrack,2.));//b tr
           if(origD0==4&&origTr>=1&&origTr<=2) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,PtWeig(ptbin,pttrack));//c tr
           if(origD0==5&&origTr>=3&&origTr<=6) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0,PtWeig(ptbin,pttrack));//b tr
        }
        if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==-421) { //D0bar
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig_Bin%d",ptbin)))->Fill(deltaphi,mD0bar,pttrack);
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,pttrack);
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr_Bin%d",ptbin)))->Fill(deltaphi,mD0bar,pow(pttrack,2.));
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,pow(pttrack,2.));
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist_Bin%d",ptbin)))->Fill(deltaphi,mD0bar,PtWeig(ptbin,pttrack));
           ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,PtWeig(ptbin,pttrack));
           if(origD0==4&&origTr>=1&&origTr<=2) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,pttrack);//c tr
           if(origD0==5&&origTr>=3&&origTr<=6) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Weig_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,pttrack);//b tr
           if(origD0==4&&origTr>=1&&origTr<=2) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,pow(pttrack,2.));//c tr
           if(origD0==5&&origTr>=3&&origTr<=6) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigSqr_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,pow(pttrack,2.));//b tr
           if(origD0==4&&origTr>=1&&origTr<=2) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,PtWeig(ptbin,pttrack));//c tr
           if(origD0==5&&origTr>=3&&origTr<=6) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_WeigDist_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(deltaphi,mD0bar,PtWeig(ptbin,pttrack));//b tr
        } 

      //hadron identification (K/pi/p from MCTruth)
      if(TMath::Abs(PDGtrack) == 321) {  
        if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==421) { //D0 (from MCTruth)
           ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg_Bin%d",ptbin)))->Fill(fillSpPhiD0);
           ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0);
           if(origD0==4&&origTr>=1&&origTr<=2) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0);
           if(origD0==5&&origTr>=3&&origTr<=6) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0);
        }
        if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==-421) { //D0bar
           ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg_Bin%d",ptbin)))->Fill(fillSpPhiD0bar);
           ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0bar);
           if(origD0==4&&origTr>=1&&origTr<=2) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0bar);
           if(origD0==5&&origTr>=3&&origTr<=6) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_Kcharg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0bar); 
        }
        if(!fAlreadyFilled) {
	  ((TH1F*)fOutputStudy->FindObject(Form("histDispl_Kcharg_Bin%d",ptbin)))->Fill(d0orig); //Fills displacement histos
  	  if (origTr>=1&&origTr<=6) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_Kcharg_HF_Bin%d",ptbin)))->Fill(d0orig);
          if (origTr>=1&&origTr<=6) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_Kcharg_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(d0orig);
	  ((TH1F*)fOutputStudy->FindObject(Form("histDispl_Kcharg%s_Bin%d",orig.Data(),ptbin)))->Fill(d0orig); //Fills displacement histos
          ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_Kcharg_Bin%d",ptbin)))->Fill(track->Pt());
          ((TH1F*)fOutputStudy->FindObject(Form("histOrig_Kcharg_Bin%d",ptbin)))->Fill(origTr);
	}
        N_Kcharg++;

      } //end hadron case

    } //end tracks loop

    //Kaon identification
    TClonesArray *v0array = (TClonesArray*)aod->GetList()->FindObject("v0s");
    Int_t idArrayV0[v0array->GetEntriesFast()][2]; //array for skipping K0 double-counting
    for(int iV0=0; iV0<v0array->GetEntriesFast();iV0++) {
      for(int j=0; j<2; j++) {idArrayV0[iV0][j]=-2;}
      AliAODv0 *v0 = (AliAODv0*)v0array->UncheckedAt(iV0);

      if(v0->Pt() < fPtThreshLow.at(ptbin) || v0->Pt() > fPtThreshUp.at(ptbin)) continue; //discard tracks outside pt range for hadrons/K
      if(SelectV0(v0,vtx,1,idArrayV0)) continue; //option 1 = for correlations
      Int_t pdgCodes[2] = {211,211};
      Int_t labV0 = v0->MatchToMC(311,mcArray,2,pdgCodes);
      if(labV0<=0) continue;

      Double_t ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_Bin%d",ptbin)))->GetAxis(2)->GetXmax(); //all plots have same axes...
      Double_t ptV0=v0->Pt(), deltaphiV0=d->Phi()-v0->Phi();
      deltaphiV0 = d->Phi()-v0->Phi();
      if (deltaphiV0 < -1.571) deltaphiV0+=6.283;
      if (deltaphiV0 > 4.712) deltaphiV0-=6.283;
      if (ptV0 > ptLim_Sparse) ptV0 = ptLim_Sparse-0.01;

      Double_t fillSpPhiD0K0[4] = {deltaphiV0,mD0,ptV0,0.};
      Double_t fillSpPhiD0barK0[4] = {deltaphiV0,mD0bar,ptV0,0.};

      Int_t origV0 = CheckTrackOrigin(mcArray,(AliAODMCParticle*)mcArray->At(labV0));

      if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==421) { //D0 (from MCTruth)
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_Bin%d",ptbin)))->Fill(fillSpPhiD0K0);
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0K0);
         if(origD0==4&&origV0>=1&&origV0<=2) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0K0);  
         if(origD0==5&&origV0>=3&&origV0<=6) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0K0);  
      }
      if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==-421) { //D0bar
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_Bin%d",ptbin)))->Fill(fillSpPhiD0barK0);
         ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0barK0);
         if(origD0==4&&origV0>=1&&origV0<=2) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0barK0);  
         if(origD0==5&&origV0>=3&&origV0<=6) ((THnSparseI*)fOutputCorr->FindObject(Form("hPhi_K_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(fillSpPhiD0barK0); 
      }
      if(!fAlreadyFilled) {
         ((TH1F*)fOutputStudy->FindObject(Form("histDispl_K_Bin%d",ptbin)))->Fill(0.); //Fills displacement histos
	 if (origV0>=1&&origV0<=6) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_K_HF_Bin%d",ptbin)))->Fill(0.);
         if (origV0>=1&&origV0<=6) ((TH1F*)fOutputStudy->FindObject(Form("histDispl_K_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(0.);
	 ((TH1F*)fOutputStudy->FindObject(Form("histDispl_K%s_Bin%d",orig.Data(),ptbin)))->Fill(0.); //Fills displacement histos
         ((TH1F*)fOutputStudy->FindObject(Form("hist_Pt_K_Bin%d",ptbin)))->Fill(v0->Pt());
         ((TH1F*)fOutputStudy->FindObject(Form("histOrig_K_Bin%d",ptbin)))->Fill(origV0);
      }
      N_Kaons++;
    } // end kaon case

    //leading track correlations fill
    if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==421) { //D0
      ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(lead[1],mD0); //c and b D0
      ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(lead[1],mD0); //c or b D0
      if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=2) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(lead[1],mD0);  
      if(origD0==5&&(int)lead[2]>=3&&(int)lead[2]<=6) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(lead[1],mD0);  
    }
    if(((AliAODMCParticle*)mcArray->At(labD0))->GetPdgCode()==-421) { //D0bar
      ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead_Bin%d",ptbin)))->Fill(lead[1],mD0bar);
      ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead%s_Bin%d",orig.Data(),ptbin)))->Fill(lead[1],mD0bar); //c or b D0
      if(origD0==4&&(int)lead[2]>=1&&(int)lead[2]<=2) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(lead[1],mD0bar);  
      if(origD0==5&&(int)lead[2]>=3&&(int)lead[2]<=6) ((TH2F*)fOutputCorr->FindObject(Form("hPhi_Lead_HF%s_Bin%d",orig.Data(),ptbin)))->Fill(lead[1],mD0bar); 
    }

  } //end fReadMC = 1

  if (!fAlreadyFilled) {
    ((TH1F*)fOutputStudy->FindObject(Form("hist_Count_Charg_Bin%d",ptbin)))->Fill(N_Charg);
    ((TH1F*)fOutputStudy->FindObject(Form("hist_Count_Kcharg_Bin%d",ptbin)))->Fill(N_Kcharg);
    ((TH1F*)fOutputStudy->FindObject(Form("hist_Count_K_Bin%d",ptbin)))->Fill(N_Kaons);
  }

  fAlreadyFilled=kTRUE; //distribution plots for tracks filled  

}

//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSED0Correlations::CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checks on particle (#) origin:
  // 0) Not HF
  // 1) D->#
  // 2) D->X->#
  // 3) B->#
  // 4) B->X-># (X!=D)
  // 5) B->D->#
  // 6) B->D->X->#
  // 7) c hadronization
  // 8) b hadronization
  //
  printf(" AliAnalysisTaskSED0Correlations::CheckD0Origin() \n");
	
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isDdaugh=kFALSE;
  Bool_t isDchaindaugh=kFALSE;
  Bool_t isBdaugh=kFALSE;
  Bool_t isBchaindaugh=kFALSE;
  Bool_t isQuarkFound=kFALSE;

  while (mother > 0){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isBchaindaugh=kTRUE;
        if(istep==1) isBdaugh=kTRUE;
      }
      if ((abspdgGranma > 400 && abspdgGranma < 500) || (abspdgGranma > 4000 && abspdgGranma < 5000)){
	isDchaindaugh=kTRUE;
        if(istep==1) isDdaugh=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) {isQuarkFound=kTRUE; if(abspdgGranma==5) isFromB = kTRUE;}
      mother = mcGranma->GetMother();
    }else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  //decides what to return based on the flag status
  if(isQuarkFound) {
    if(!isFromB) {  //charm
      if(isDdaugh) return 1; //charm immediate
      else if(isDchaindaugh) return 2; //charm chain
      else return 7; //charm hadronization
    }
    else { //beauty
      if(isBdaugh) return 3; //b immediate
      else if(isBchaindaugh) { //b chain
        if(isDchaindaugh) {
          if(isDdaugh) return 5; //d immediate
          return 6; //d chain
          }
        else return 4; //b, not d
      }
      else return 8; //b hadronization
    }
  }
  else return 0; //no HF quark
}


//________________________________________________________________________
Int_t AliAnalysisTaskSED0Correlations::PtBinCorr(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin=-1;
  if(pt<fBinLimsCorr.at(0)) return ptbin; //out of bounds
  
  Int_t i = 0;
  while(pt>fBinLimsCorr.at(i)) {ptbin=i; i++;}
  
  return ptbin;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSED0Correlations::PtWeig(Int_t ptbin, Double_t x) const {
  //
  //gives the weight for Weighted Corrs (inverse of pT distribution, from 3 D0 pt bins)
  //
  if(x>11.5) x=11.5; if(x<0.3) x=0.3; //bounds for fit of distributions!
  if(ptbin >= 3 && ptbin <= 5)  return 1/(0.22958*(1.507e+03-5.035e+02*x+5.681e+01*x*x-2.186e+00*x*x*x+exp(1.336e+01-2.146e+00*x+1.206e-01*x*x)));
  if(ptbin >= 6 && ptbin <= 8)  return 1/(1.71828*(1.984e+02-6.113e+01*x+6.444e+00*x*x-2.342e-01*x*x*x+exp(1.023e+01-2.059e+00*x+1.194e-01*x*x)));
  if(ptbin >= 9 && ptbin <= 10) return 1/(1.25905*(1.990e+02-6.306e+01*x+6.995e+00*x*x-2.681e-01*x*x*x+exp(9.125e+00-2.053e+00*x+1.276e-01*x*x)));

  return 0; //for other bins plot is disabled!
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::GetImpParameter(AliAODTrack *track, AliAODEvent* aod, Double_t &d0, Double_t &d0err) const {
  //
  //calculates d0 and error on d0 for the track
  //      
  Double_t  d0z0[2],covd0z0[3];
  AliESDtrack esdTrack(track);  // AliESDTrack conversion of AOD track
  esdTrack.PropagateToDCA(aod->GetPrimaryVertex(),aod->GetMagneticField(),10000.,d0z0,covd0z0); 
  d0 = TMath::Abs(d0z0[0]); // impact parameter
  d0err = TMath::Sqrt(covd0z0[0]); // resolution of impact parameter
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSED0Correlations::TrackSelectedInLoop(AliAODTrack* track, AliAODRecoDecayHF2Prong *d, AliAODEvent *aod, Int_t ptbin, Int_t idDaughs[]) const {
  //
  //rejection of tracks in loop for correlations
  //      
  Bool_t output = kTRUE;
  AliAODVertex *vtx = (AliAODVertex*)aod->GetPrimaryVertex();
  Double_t bz = aod->GetMagneticField();
  if(track->GetID() == idDaughs[0] || track->GetID() == idDaughs[1] || track->GetID() < 0) output = kFALSE; //discards daughters of candidate
  if(track->Pt() < fPtThreshLow.at(ptbin) || track->Pt() > fPtThreshUp.at(ptbin)) output = kFALSE; //discard tracks outside pt range for hadrons/K
  if(!fCutsTracks->IsHadronSelected(track,vtx,bz)) { //track discarded by quality cuts
    ((TH1F*)fOutputStudy->FindObject("hRejTracks"))->Fill(1);
    output = kFALSE;
  } else ((TH1F*)fOutputStudy->FindObject("hRejTracks"))->Fill(0); //track accepted by quality cuts
  if(!fCutsTracks->InvMassDstarRejection(d,track,fIsSelectedCandidate)) {
    ((TH1F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(1);
    output = kFALSE; 
  } else ((TH1F*)fOutputStudy->FindObject(Form("hDstarPions_Bin%d",ptbin)))->Fill(0); //rejects possible pions from D* using inv.mass

  return output;
}

//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskSED0Correlations::SelectV0(AliAODv0* v0, AliAODVertex *vtx, Int_t opt, Int_t idArrayV0[][2]) const
{
  //
  // Selection for K0 hypotheses
  // options: 1 = selects mass invariant about 3 sigma inside the peak + threshold of 0.3 GeV
  // 	      2 = no previous selections

  if(!fCutsTracks->IsKZeroSelected(v0,vtx)) return kFALSE;

  AliAODTrack *v0Daug1 = (AliAODTrack*)v0->GetDaughter(0);
  AliAODTrack *v0Daug2 = (AliAODTrack*)v0->GetDaughter(1);

  if(opt==1) { //additional cuts for correlations (V0 has to be closer than 3 sigma from K0 mass)
    if(v0->Pt() < 3. && TMath::Abs(v0->MassK0Short()-0.4976) > 3*0.0040) return kFALSE;
    if(v0->Pt() > 3. && v0->Pt() < 6. && TMath::Abs(v0->MassK0Short()-0.4976) > 3*0.0052) return kFALSE;
    if(v0->Pt() > 6. && TMath::Abs(v0->MassK0Short()-0.4976) > 3*0.0075) return kFALSE;
  }

  //This part removes double counting for swapped tracks!
  Int_t i = 0;  //while loop (until the last-written entry pair of ID!
  while(idArrayV0[i][0]!=-2 && idArrayV0[i][1]!=-2) {
    if((v0Daug1->GetID()==idArrayV0[i][0] && v0Daug2->GetID()==idArrayV0[i][1])||
       (v0Daug1->GetID()==idArrayV0[i][1] && v0Daug2->GetID()==idArrayV0[i][0])) return kFALSE;
    i++;
  }
  idArrayV0[i][0]=v0Daug1->GetID();
  idArrayV0[i][1]=v0Daug2->GetID();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Correlations::PrintBinsAndLimits() {

  cout << "--------------------------\n";
  cout << "PtBins = " << fNPtBinsCorr << "\n";
  cout << "PtBin limits--------------\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << "Bin "<<i+1<<" = "<<fBinLimsCorr.at(i)<<" to "<<fBinLimsCorr.at(i+1)<<"\n";
  }
  cout << "\n--------------------------\n";
  cout << "PtBin tresh. tracks low---\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fPtThreshLow.at(i) << "; ";
  }
  cout << "PtBin tresh. tracks up----\n";
  for (int i=0; i<fNPtBinsCorr; i++) {
    cout << fPtThreshUp.at(i) << "; ";
  }
  cout << "\n--------------------------\n";
  cout << "MC Truth = "<<fReadMC<<"\n";
  cout << "--------------------------\n";
}

