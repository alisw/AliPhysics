//------------------------------------------------------------------------------
// Implementation of AliPerformanceEff class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is 
// a data member of AliPerformanceEff.
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();
  TFile f("Output.root");
  AliPerformanceEff * compObj = (AliPerformanceEff*)coutput->FindObject("AliPerformanceEff");

  // Analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderEff" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_Eff.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include <TAxis.h>
#include <TH1D.h>

// 
#include "AliESDtrack.h"
#include "AliRecInfoCuts.h" 
#include "AliMCInfoCuts.h" 
#include "AliLog.h" 
#include "AliESDVertex.h" 
#include "AliExternalTrackParam.h" 
#include "AliTracker.h" 
#include "AliESDEvent.h" 
#include "AliMCEvent.h" 
#include "AliMCParticle.h" 
#include "AliHeader.h" 
#include "AliGenEventHeader.h" 
#include "AliStack.h" 
#include "AliPerformanceEff.h" 

using namespace std;

ClassImp(AliPerformanceEff)

//_____________________________________________________________________________
AliPerformanceEff::AliPerformanceEff():
  AliPerformanceObject("AliPerformanceEff"),

  // histograms
  fEffHisto(0),
  fEffSecHisto(0),

  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),

  // histogram folder 
  fAnalysisFolder(0)
{
  // default consttructor	
  Init();
}

//_____________________________________________________________________________
AliPerformanceEff::AliPerformanceEff(Char_t* name="AliPerformanceEff",Char_t*title="AliPerformanceEff",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE):
  AliPerformanceObject(name,title),

  // histograms
  fEffHisto(0),
  fEffSecHisto(0),

  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),

  // histogram folder 
  fAnalysisFolder(0)
{
  // named constructor
  //
  SetAnalysisMode(analysisMode);
  SetHptGenerator(hptGenerator);

  Init();
}


//_____________________________________________________________________________
AliPerformanceEff::~AliPerformanceEff()
{
// destructor

  if(fEffHisto)  delete  fEffHisto; fEffHisto=0;
  if(fEffSecHisto)  delete  fEffSecHisto; fEffSecHisto=0;
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliPerformanceEff::Init()
{
  // Init histograms
  //
  // set pt bins
  Int_t nPtBins = 50;
  Double_t ptMin = 1.e-2, ptMax = 10.;

  Double_t *binsPt = 0;
  if (IsHptGenerator())  { 
    nPtBins = 100; ptMax = 100.;
    binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);
  } else {
    binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);
  }

  /*
  Int_t nPtBins = 31;
  Double_t binsPt[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
  Double_t ptMin = 0., ptMax = 10.;

  if(IsHptGenerator() == kTRUE) {
    nPtBins = 100;
    ptMin = 0.; ptMax = 100.;
  }
  */

  //mceta:mcphi:mcpt:pid:recStatus:findable
  Int_t binsEffHisto[6]={30,90,nPtBins,5,2,2};
  Double_t minEffHisto[6]={-1.5,0.,ptMin,0.,0.,0.};
  Double_t maxEffHisto[6]={ 1.5,2.*TMath::Pi(), ptMax,5.,2.,2.};

  fEffHisto = new THnSparseF("fEffHisto","mceta:mcphi:mcpt:pid:recStatus:findable",6,binsEffHisto,minEffHisto,maxEffHisto);
  fEffHisto->SetBinEdges(2,binsPt);

  fEffHisto->GetAxis(0)->SetTitle("#eta_{mc}");
  fEffHisto->GetAxis(1)->SetTitle("#phi_{mc} (rad)");
  fEffHisto->GetAxis(2)->SetTitle("p_{Tmc} (GeV/c)");
  fEffHisto->GetAxis(3)->SetTitle("pid");
  fEffHisto->GetAxis(4)->SetTitle("recStatus");
  fEffHisto->GetAxis(5)->SetTitle("findable");
  fEffHisto->Sumw2();

  //mceta:mcphi:mcpt:pid:recStatus:findable:mcR:mother_phi:mother_eta
  Int_t binsEffSecHisto[9]={30,60,nPtBins,5,2,2,100,60,30};
  Double_t minEffSecHisto[9]={-1.5,0.,ptMin,0.,0.,0.,0.,0.,-1.5};
  Double_t maxEffSecHisto[9]={ 1.5,2.*TMath::Pi(), ptMax,5.,2.,2.,200,2.*TMath::Pi(),1.5};

  fEffSecHisto = new THnSparseF("fEffSecHisto","mceta:mcphi:mcpt:pid:recStatus:findable:mcR:mother_phi:mother_eta",9,binsEffSecHisto,minEffSecHisto,maxEffSecHisto);
  fEffSecHisto->SetBinEdges(2,binsPt);

  fEffSecHisto->GetAxis(0)->SetTitle("#eta_{mc}");
  fEffSecHisto->GetAxis(1)->SetTitle("#phi_{mc} (rad)");
  fEffSecHisto->GetAxis(2)->SetTitle("p_{Tmc} (GeV/c)");
  fEffSecHisto->GetAxis(3)->SetTitle("pid");
  fEffSecHisto->GetAxis(4)->SetTitle("recStatus");
  fEffSecHisto->GetAxis(5)->SetTitle("findable");
  fEffSecHisto->GetAxis(6)->SetTitle("mcR (cm)");
  fEffSecHisto->GetAxis(7)->SetTitle("mother_phi (rad)");
  fEffSecHisto->GetAxis(8)->SetTitle("mother_eta");
  fEffSecHisto->Sumw2();

  // init cuts
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderEff","Analysis Efficiency Folder");
}

//_____________________________________________________________________________
void AliPerformanceEff::ProcessTPC(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent)
{
  // Fill TPC only efficiency comparison information 
  Int_t *labelsRec =  new Int_t[esdEvent->GetNumberOfTracks()];
  if(!labelsRec) 
     AliDebug(AliLog::kError, "Cannot create labelsRec");

  Int_t *labelsAllRec =  new Int_t[esdEvent->GetNumberOfTracks()];
  if(!labelsAllRec) 
     AliDebug(AliLog::kError, "Cannot create labelsAllRec");

  // loop over rec. tracks
  AliESDtrack *track=0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    track = esdEvent->GetTrack(iTrack);
    if(!track) continue;
    if(track->Charge()==0) continue;
    Int_t label = TMath::Abs(track->GetLabel()); 
    labelsAllRec[iTrack]=label;

    // TPC only
    if(IsRecTPC(track) != 0) 
      labelsRec[iTrack]=label;
  }

  // 
  // MC histograms for efficiency studies
  //
  if(!mcEvent) return;
 
  AliStack *stack = mcEvent->Stack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    return;
  }

  //Int_t nPart  = stack->GetNtrack();
  Int_t nPart  = stack->GetNprimary();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;
    if (particle->GetPDG()->Charge() == 0.0) continue;
      
    // physical primary
    Bool_t prim = stack->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    Bool_t findable = kFALSE;
    for(Int_t iRec=0; iRec<esdEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check findable
      if(iMc == labelsAllRec[iRec]) 
      {
        findable = IsFindable(mcEvent,iMc);
	break;
      }
    }  

    Bool_t recStatus = kFALSE;
    for(Int_t iRec=0; iRec<esdEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check reconstructed
      if(iMc == labelsRec[iRec]) 
      {
        recStatus = kTRUE;
        break;
      }
    }

    // Only 5 charged particle species (e,mu,pi,K,p)
    if (fCutsMC->IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) continue; 

    // transform Pdg to Pid
    Int_t pid = TransformToPID(particle);

    Float_t mceta =  particle->Eta();
    Float_t mcphi =  particle->Phi();
    if(mcphi<0) mcphi += 2.*TMath::Pi();
    Float_t mcpt = particle->Pt();


    // Fill histograms
    Double_t vEffHisto[6] = {mceta, mcphi, mcpt, pid, recStatus, findable}; 
    fEffHisto->Fill(vEffHisto);
  }

  if(labelsRec) delete [] labelsRec; labelsRec = 0;
  if(labelsAllRec) delete [] labelsAllRec; labelsAllRec = 0;
}

//_____________________________________________________________________________
void AliPerformanceEff::ProcessTPCSec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent)
{
  // Fill TPC only efficiency comparison information for secondaries
  Int_t *labelsRec =  new Int_t[esdEvent->GetNumberOfTracks()];
  if(!labelsRec) 
     AliDebug(AliLog::kError, "Cannot create labelsRec");

  Int_t *labelsAllRec =  new Int_t[esdEvent->GetNumberOfTracks()];
  if(!labelsAllRec) 
     AliDebug(AliLog::kError, "Cannot create labelsAllRec");

  // loop over rec. tracks
  AliESDtrack *track=0;
  Int_t multAll=0, multRec=0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    track = esdEvent->GetTrack(iTrack);
    if(!track) continue;
    if(track->Charge()==0) continue;
    Int_t label = TMath::Abs(track->GetLabel()); 
    labelsAllRec[multAll]=label;
    multAll++;

    // TPC only
    if(IsRecTPC(track) != 0) {
      labelsRec[multRec]=label;
      multRec++;
    }
  }

  // 
  // MC histograms for efficiency studies
  //
  if(!mcEvent) return;
 
  AliStack *stack = mcEvent->Stack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    return;
  }

  Int_t nPart  = stack->GetNtrack();
  //Int_t nPart  = stack->GetNprimary();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;
    if (particle->GetPDG()->Charge() == 0.0) continue;
      
    // physical primary
    Bool_t prim = stack->IsPhysicalPrimary(iMc);

    // only secondaries which can be reconstructed at TPC
    if(prim) continue;

    //Float_t radius = TMath::Sqrt(particle->Vx()*particle->Vx()+particle->Vy()*particle->Vy()+particle->Vz()*particle->Vz());
    //if(radius > fCutsMC->GetMaxR()) continue;

    // only secondary electrons from gamma conversion
    //if( TMath::Abs(particle->GetPdgCode())!=fCutsMC->GetEM() ||   particle->GetUniqueID() != 5) continue;

    Bool_t findable = kFALSE;
    for(Int_t iRec=0; iRec<multAll; ++iRec) 
    {
      // check findable
      if(iMc == labelsAllRec[iRec]) 
      {
        findable = IsFindable(mcEvent,iMc);
	break;
      }
    }  

    Bool_t recStatus = kFALSE;
    for(Int_t iRec=0; iRec<multRec; ++iRec) 
    {
      // check reconstructed
      if(iMc == labelsRec[iRec]) 
      {
        recStatus = kTRUE;
        break;
      }
    }

    // Only 5 charged particle species (e,mu,pi,K,p)
    if (fCutsMC->IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) continue; 

    // transform Pdg to Pid
    Int_t pid = TransformToPID(particle);

    Float_t mceta =  particle->Eta();
    Float_t mcphi =  particle->Phi();
    if(mcphi<0) mcphi += 2.*TMath::Pi();
    Float_t mcpt = particle->Pt();
    Float_t mcR = particle->R();

    // get info about mother
    Int_t motherLabel = particle->GetMother(0);
    if(motherLabel < 0) continue;
    TParticle *mother = stack->Particle(motherLabel);
    if(!mother) continue; 

    Float_t mother_eta = mother->Eta();
    Float_t mother_phi = mother->Phi();
    if(mother_phi<0) mother_phi += 2.*TMath::Pi();

    // Fill histograms
    Double_t vEffSecHisto[9] = { mceta, mcphi, mcpt, pid, recStatus, findable, mcR, mother_phi, mother_eta }; 
    fEffSecHisto->Fill(vEffSecHisto);
  }

  if(labelsRec) delete [] labelsRec; labelsRec = 0;
  if(labelsAllRec) delete [] labelsAllRec; labelsAllRec = 0;
}




//_____________________________________________________________________________
void AliPerformanceEff::ProcessTPCITS(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent)
{
  // Fill efficiency comparison information
  Int_t *labelsRec =  new Int_t[esdEvent->GetNumberOfTracks()];
  if(!labelsRec) 
     AliDebug(AliLog::kError, "Cannot create labelsRec");

  Int_t *labelsAllRec =  new Int_t[esdEvent->GetNumberOfTracks()];
  if(!labelsAllRec) 
     AliDebug(AliLog::kError, "Cannot create labelsAllRec");

  // loop over rec. tracks
  AliESDtrack *track=0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    track = esdEvent->GetTrack(iTrack);
    if(!track) continue;
    if(track->Charge()==0) continue;
    Int_t label = TMath::Abs(track->GetLabel()); 
    labelsAllRec[iTrack]=label;

    // iTPC+ITS
    if(IsRecTPCITS(track) != 0) 
      labelsRec[iTrack]=label;
  }

  // 
  // MC histograms for efficiency studies
  //
  if(!mcEvent) return;
 
  AliStack *stack = mcEvent->Stack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    return;
  }

  //Int_t nPart  = stack->GetNtrack();
  Int_t nPart  = stack->GetNprimary();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;
    if (particle->GetPDG()->Charge() == 0.0) continue;
      
    // physical primary
    Bool_t prim = stack->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    Bool_t findable = kFALSE;
    for(Int_t iRec=0; iRec<esdEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check findable
      if(iMc == labelsAllRec[iRec]) 
      {
        findable = IsFindable(mcEvent,iMc);
	break;
      }
    }  

    Bool_t recStatus = kFALSE;
    for(Int_t iRec=0; iRec<esdEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check reconstructed
      if(iMc == labelsRec[iRec]) 
      {
        recStatus = kTRUE;
        break;
      }
    }

    // Only 5 charged particle species (e,mu,pi,K,p)
    if (fCutsMC->IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) continue; 

    // transform Pdg to Pid
    Int_t pid = TransformToPID(particle);

    Float_t mceta =  particle->Eta();
    Float_t mcphi =  particle->Phi();
    if(mcphi<0) mcphi += 2.*TMath::Pi();
    Float_t mcpt = particle->Pt();

    // Fill histograms
    Double_t vEffHisto[6] = { mceta, mcphi, mcpt, pid, recStatus, findable}; 
    fEffHisto->Fill(vEffHisto);
  }

  if(labelsRec) delete [] labelsRec; labelsRec = 0;
  if(labelsAllRec) delete [] labelsAllRec; labelsAllRec = 0;
}

//_____________________________________________________________________________
void AliPerformanceEff::ProcessConstrained(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent)
{
  // Process comparison information 
  Int_t *labelsRec =  new Int_t[esdEvent->GetNumberOfTracks()];
  if(!labelsRec) 
     AliDebug(AliLog::kError, "Cannot create labelsRec");

  Int_t *labelsAllRec =  new Int_t[esdEvent->GetNumberOfTracks()];
  if(!labelsAllRec) 
     AliDebug(AliLog::kError, "Cannot create labelsAllRec");

  // loop over rec. tracks
  AliESDtrack *track=0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    track = esdEvent->GetTrack(iTrack);
    if(!track) continue;
    if(track->Charge()==0) continue;
    Int_t label = TMath::Abs(track->GetLabel()); 
    labelsAllRec[iTrack]=label;

    // Constrained
    if(IsRecConstrained(track) != 0) 
      labelsRec[iTrack]=label;

  }

  // 
  // MC histograms for efficiency studies
  //
  if(!mcEvent) return;
 
  AliStack *stack = mcEvent->Stack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    return;
  }

  //Int_t nPart  = stack->GetNtrack();
  Int_t nPart  = stack->GetNprimary();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;
    if (particle->GetPDG()->Charge() == 0.0) continue;
      
    // physical primary
    Bool_t prim = stack->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    Bool_t findable = kFALSE;
    for(Int_t iRec=0; iRec<esdEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check findable
      if(iMc == labelsAllRec[iRec]) 
      {
        findable = IsFindable(mcEvent,iMc);
	break;
      }
    }  

    Bool_t recStatus = kFALSE;
    for(Int_t iRec=0; iRec<esdEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check reconstructed
      if(iMc == labelsRec[iRec]) 
      {
        recStatus = kTRUE;
        break;
      }
    }

    // Only 5 charged particle species (e,mu,pi,K,p)
    if (fCutsMC->IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) continue; 

    // transform Pdg to Pid
    Int_t pid = TransformToPID(particle);

    Float_t mceta =  particle->Eta();
    Float_t mcphi =  particle->Phi();
    if(mcphi<0) mcphi += 2.*TMath::Pi();
    Float_t mcpt = particle->Pt();

    // Fill histograms
    Double_t vEffHisto[6] = { mceta, mcphi, mcpt, pid, recStatus, findable}; 
    fEffHisto->Fill(vEffHisto);
  }

  if(labelsRec) delete [] labelsRec; labelsRec = 0;
  if(labelsAllRec) delete [] labelsAllRec; labelsAllRec = 0;
}

//_____________________________________________________________________________
void AliPerformanceEff::Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, const Bool_t bUseMC, const Bool_t bUseESDfriend)
{
  // Process comparison information 
  //
  if(!esdEvent) 
  {
    Error("Exec","esdEvent not available");
    return;
  }
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  
  if(bUseMC)
  {
    if(!mcEvent) {
      Error("Exec","mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      Error("Exec","Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      Error("Exec","Stack not available");
      return;
    }
    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      Error("Exec","Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);
  } 
  else {
    Error("Exec","MC information required!");
    return;
  } 

  // use ESD friends
  if(bUseESDfriend) {
    if(!esdFriend) {
      Error("Exec","esdFriend not available");
      return;
    }
  }

  //
  //  Process events
  //
  if(GetAnalysisMode() == 0) ProcessTPC(mcEvent,esdEvent);
  else if(GetAnalysisMode() == 1) ProcessTPCITS(mcEvent,esdEvent);
  else if(GetAnalysisMode() == 2) ProcessConstrained(mcEvent,esdEvent);
  else if(GetAnalysisMode() == 5) ProcessTPCSec(mcEvent,esdEvent);
  else {
    printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
    return;
  }
}

//_____________________________________________________________________________
Int_t AliPerformanceEff::TransformToPID(TParticle *particle) 
{
// transform Pdg to Pid
// Pdg convension is different for hadrons and leptons 
// (e.g. K+/K- = 321/-321; e+/e- = -11/11 ) 

  Int_t pid = -1;
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC->GetEM() ) pid = 0; 
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC->GetMuM() ) pid = 1; 
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC->GetPiP() ) pid = 2; 
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC->GetKP() ) pid = 3; 
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC->GetProt() ) pid = 4; 

return pid;
}

//_____________________________________________________________________________
Bool_t AliPerformanceEff::IsFindable(AliMCEvent *mcEvent, Int_t label) 
{
if(!mcEvent) return kFALSE;

  AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(label);
  if(!mcParticle) return kFALSE;

  Int_t counter; 
  Float_t tpcTrackLength = mcParticle->GetTPCTrackLength(AliTracker::GetBz(),0.05,counter,3.0); 
  //printf("tpcTrackLength %f \n", tpcTrackLength);

return (tpcTrackLength>fCutsMC->GetMinTrackLength());    
}

//_____________________________________________________________________________
Bool_t AliPerformanceEff::IsRecTPC(AliESDtrack *esdTrack) 
{
if(!esdTrack) return kFALSE;

  const AliExternalTrackParam *track = esdTrack->GetTPCInnerParam();
  if(!track) return kFALSE;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParametersTPC(dca,cov);

  Bool_t recStatus = kFALSE;
  if(esdTrack->GetTPCNcls()>fCutsRC->GetMinNClustersTPC()) recStatus = kTRUE; 

  /*
  if( TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && 
      TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ())
  {
    recStatus = kTRUE;
  }
  */

return recStatus;
}

//_____________________________________________________________________________
Bool_t AliPerformanceEff::IsRecTPCITS(AliESDtrack *esdTrack) 
{
if(!esdTrack) return kFALSE;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParameters(dca,cov);

  Bool_t recStatus = kFALSE;

  if ((esdTrack->GetStatus()&AliESDtrack::kTPCrefit)==0) return kFALSE; // TPC refit
  if (esdTrack->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return kFALSE; // min. nb. TPC clusters
  //if ((esdTrack->GetStatus()&AliESDtrack::kITSrefit)==0) return kFALSE; // ITS refit
  //Int_t clusterITS[200];
  //if(esdTrack->GetITSclusters(clusterITS)<fCutsRC->GetMinNClustersITS()) return kFALSE;  // min. nb. ITS clusters

  recStatus = kTRUE;
  /*
  if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && 
     TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ())
  {
    recStatus = kTRUE;
  }
  */

return recStatus;
}

//_____________________________________________________________________________
Bool_t AliPerformanceEff::IsRecConstrained(AliESDtrack *esdTrack) 
{
  if(!esdTrack) return kFALSE;

  const AliExternalTrackParam * track = esdTrack->GetConstrainedParam();
  if(!track) return kFALSE;

  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  esdTrack->GetImpactParameters(dca,cov);
  //Int_t label = TMath::Abs(esdTrack->GetLabel()); 

  Bool_t recStatus = kFALSE;

  if ((esdTrack->GetStatus()&AliESDtrack::kTPCrefit)==0) return kFALSE; // TPC refit
  if (esdTrack->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return kFALSE; // min. nb. TPC clusters
  Int_t clusterITS[200];
  if(esdTrack->GetITSclusters(clusterITS)<fCutsRC->GetMinNClustersITS()) return kFALSE;  // min. nb. ITS clusters

  recStatus = kTRUE;

  /*
  if(TMath::Abs(dca[0])<fCutsRC->GetMaxDCAToVertexXY() && 
     TMath::Abs(dca[1])<fCutsRC->GetMaxDCAToVertexZ())
  {
    recStatus = kTRUE;
  }
  */

return recStatus;
}

//_____________________________________________________________________________
Long64_t AliPerformanceEff::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  // collection of generated histograms

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliPerformanceEff* entry = dynamic_cast<AliPerformanceEff*>(obj);
    if (entry == 0) continue; 
  
     fEffHisto->Add(entry->fEffHisto);
     fEffSecHisto->Add(entry->fEffSecHisto);
  count++;
  }

return count;
}
 
//_____________________________________________________________________________
void AliPerformanceEff::Analyse() 
{
  // Analyse comparison information and store output histograms
  // in the folder "folderEff" 
  //
  TH1::AddDirectory(kFALSE);
  TObjArray *aFolderObj = new TObjArray;
  char title[256];

  //
  // efficiency vs pt
  //

  if(GetAnalysisMode() != 5) {

  fEffHisto->GetAxis(0)->SetRangeUser(-0.9,0.9); // eta range
  fEffHisto->GetAxis(2)->SetRangeUser(0.1,10.); // pt range

  // rec efficiency vs pt
  TH1D *ptAll = fEffHisto->Projection(2);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.);  // reconstructed 
  TH1D *ptRec = fEffHisto->Projection(2);
  TH1D *ptRec_c = (TH1D*)ptRec->Clone();
  ptRec_c->Divide(ptRec,ptAll,1,1,"B");
  ptRec_c->SetName("ptRecEff");

  ptRec_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(2)->GetTitle());
  ptRec_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency",fEffHisto->GetAxis(2)->GetTitle());
  ptRec_c->SetTitle(title);

  ptRec_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRec_c);

  // rec efficiency vs pid vs pt

  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions

  TH1D *ptAllPi = fEffHisto->Projection(2);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecPi = fEffHisto->Projection(2);
  TH1D *ptRecPi_c = (TH1D*)ptRecPi->Clone();
  ptRecPi_c->Divide(ptRecPi,ptAllPi,1,1,"B");
  ptRecPi_c->SetName("ptRecEffPi");

  ptRecPi_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(2)->GetTitle());
  ptRecPi_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (pions)",fEffHisto->GetAxis(2)->GetTitle());
  ptRecPi_c->SetTitle(title);

  ptRecPi_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecPi_c);

  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(3.,3.); // kaons
  TH1D *ptAllK = fEffHisto->Projection(2);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecK = fEffHisto->Projection(2);

  TH1D *ptRecK_c = (TH1D*)ptRecK->Clone();
  ptRecK_c->Divide(ptRecK,ptAllK,1,1,"B");
  ptRecK_c->SetName("ptRecEffK");

  ptRecK_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(2)->GetTitle());
  ptRecK_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (kaons)",fEffHisto->GetAxis(2)->GetTitle());
  ptRecK_c->SetTitle(title);


  ptRecK_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecK_c);

  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(4.,4.); // protons
  TH1D *ptAllP = fEffHisto->Projection(2);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecP = fEffHisto->Projection(2);
  TH1D *ptRecP_c = (TH1D*)ptRecP->Clone();
  ptRecP_c->Divide(ptRecP,ptAllP,1,1,"B");
  ptRecP_c->SetName("ptRecEffP");

  ptRecP_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(2)->GetTitle());
  ptRecP_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (protons)",fEffHisto->GetAxis(2)->GetTitle());
  ptRecP_c->SetTitle(title);

  ptRecP_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecP_c);

  // findable efficiency vs pt

  fEffHisto->GetAxis(3)->SetRangeUser(0.,4.); 
  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // findable
  TH1D *ptAllF = fEffHisto->Projection(2);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.);
  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.);

  TH1D *ptRecF = fEffHisto->Projection(2); // rec findable
  TH1D *ptRecF_c = (TH1D*)ptRecF->Clone();
  ptRecF_c->Divide(ptRecF,ptAllF,1,1,"B");
  ptRecF_c->SetName("ptRecEffF");

  ptRecF_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(2)->GetTitle());
  ptRecF_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (findable)",fEffHisto->GetAxis(2)->GetTitle());
  ptRecF_c->SetTitle(title);

  ptRecF_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecF_c);

  //
  // efficiency vs eta
  //

  fEffHisto->GetAxis(0)->SetRangeUser(-1.5,1.5); // eta range
  fEffHisto->GetAxis(2)->SetRangeUser(0.2,10.); // pt range
  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.);   // all
  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.);   // all

  TH1D *etaAll = fEffHisto->Projection(0);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.);  // reconstructed 
  TH1D *etaRec = fEffHisto->Projection(0);
  TH1D *etaRec_c = (TH1D*)etaRec->Clone();
  etaRec_c->Divide(etaRec,etaAll,1,1,"B");
  etaRec_c->SetName("etaRecEff");

  etaRec_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(0)->GetTitle());
  etaRec_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency",fEffHisto->GetAxis(0)->GetTitle());
  etaRec_c->SetTitle(title);

  aFolderObj->Add(etaRec_c);

  // rec efficiency vs pid vs eta
  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions

  TH1D *etaAllPi = fEffHisto->Projection(0);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecPi = fEffHisto->Projection(0);
  TH1D *etaRecPi_c = (TH1D*)etaRecPi->Clone();
  etaRecPi_c->Divide(etaRecPi,etaAllPi,1,1,"B");
  etaRecPi_c->SetName("etaRecEffPi");

  etaRecPi_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(0)->GetTitle());
  etaRecPi_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (pions)",fEffHisto->GetAxis(0)->GetTitle());
  etaRecPi_c->SetTitle(title);

  aFolderObj->Add(etaRecPi_c);

  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(3.,3.); // kaons
  TH1D *etaAllK = fEffHisto->Projection(0);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecK = fEffHisto->Projection(0);

  TH1D *etaRecK_c = (TH1D*)etaRecK->Clone();
  etaRecK_c->Divide(etaRecK,etaAllK,1,1,"B");
  etaRecK_c->SetName("etaRecEffK");

  etaRecK_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(0)->GetTitle());
  etaRecK_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (kaons)",fEffHisto->GetAxis(0)->GetTitle());
  etaRecK_c->SetTitle(title);


  aFolderObj->Add(etaRecK_c);

  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(4.,4.); // protons
  TH1D *etaAllP = fEffHisto->Projection(0);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecP = fEffHisto->Projection(0);
  TH1D *etaRecP_c = (TH1D*)etaRecP->Clone();
  etaRecP_c->Divide(etaRecP,etaAllP,1,1,"B");
  etaRecP_c->SetName("etaRecEffP");

  etaRecP_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(0)->GetTitle());
  etaRecP_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (protons)",fEffHisto->GetAxis(0)->GetTitle());
  etaRecP_c->SetTitle(title);

  aFolderObj->Add(etaRecP_c);

  // findable efficiency vs eta

  fEffHisto->GetAxis(3)->SetRangeUser(0.,4.); 
  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // findable
  TH1D *etaAllF = fEffHisto->Projection(0);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.);
  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.);

  TH1D *etaRecF = fEffHisto->Projection(0); // rec findable
  TH1D *etaRecF_c = (TH1D*)etaRecF->Clone();
  etaRecF_c->Divide(etaRecF,etaAllF,1,1,"B");
  etaRecF_c->SetName("etaRecEffF");

  etaRecF_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(0)->GetTitle());
  etaRecF_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (findable)",fEffHisto->GetAxis(0)->GetTitle());
  etaRecF_c->SetTitle(title);

  aFolderObj->Add(etaRecF_c);

  //
  // efficiency vs phi
  //

  fEffHisto->GetAxis(0)->SetRangeUser(-0.9,0.9); // eta range
  fEffHisto->GetAxis(2)->SetRangeUser(0.2,10.); // pt range
  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.);   // all
  fEffHisto->GetAxis(5)->SetRangeUser(0.,1.);   // all

  TH1D *phiAll = fEffHisto->Projection(1);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.);  // reconstructed 
  TH1D *phiRec = fEffHisto->Projection(1);
  TH1D *phiRec_c = (TH1D*)phiRec->Clone();
  phiRec_c->Divide(phiRec,phiAll,1,1,"B");
  phiRec_c->SetName("phiRecEff");

  phiRec_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(1)->GetTitle());
  phiRec_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency",fEffHisto->GetAxis(1)->GetTitle());
  phiRec_c->SetTitle(title);

  aFolderObj->Add(phiRec_c);

  // rec efficiency vs pid vs phi
  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions

  TH1D *phiAllPi = fEffHisto->Projection(1);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *phiRecPi = fEffHisto->Projection(1);
  TH1D *phiRecPi_c = (TH1D*)phiRecPi->Clone();
  phiRecPi_c->Divide(phiRecPi,phiAllPi,1,1,"B");
  phiRecPi_c->SetName("phiRecEffPi");

  phiRecPi_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(1)->GetTitle());
  phiRecPi_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (pions)",fEffHisto->GetAxis(1)->GetTitle());
  phiRecPi_c->SetTitle(title);

  aFolderObj->Add(phiRecPi_c);

  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(3.,3.); // kaons
  TH1D *phiAllK = fEffHisto->Projection(1);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *phiRecK = fEffHisto->Projection(1);

  TH1D *phiRecK_c = (TH1D*)phiRecK->Clone();
  phiRecK_c->Divide(phiRecK,phiAllK,1,1,"B");
  phiRecK_c->SetName("phiRecEffK");

  phiRecK_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(1)->GetTitle());
  phiRecK_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (kaons)",fEffHisto->GetAxis(1)->GetTitle());
  phiRecK_c->SetTitle(title);


  aFolderObj->Add(phiRecK_c);

  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(3)->SetRangeUser(4.,4.); // protons
  TH1D *phiAllP = fEffHisto->Projection(1);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *phiRecP = fEffHisto->Projection(1);
  TH1D *phiRecP_c = (TH1D*)phiRecP->Clone();
  phiRecP_c->Divide(phiRecP,phiAllP,1,1,"B");
  phiRecP_c->SetName("phiRecEffP");

  phiRecP_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(1)->GetTitle());
  phiRecP_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (protons)",fEffHisto->GetAxis(1)->GetTitle());
  phiRecP_c->SetTitle(title);

  aFolderObj->Add(phiRecP_c);

  // findable efficiency vs phi

  fEffHisto->GetAxis(3)->SetRangeUser(0.,4.); 
  fEffHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.); // findable
  TH1D *phiAllF = fEffHisto->Projection(1);

  fEffHisto->GetAxis(4)->SetRangeUser(1.,1.);
  fEffHisto->GetAxis(5)->SetRangeUser(1.,1.);

  TH1D *phiRecF = fEffHisto->Projection(1); // rec findable
  TH1D *phiRecF_c = (TH1D*)phiRecF->Clone();
  phiRecF_c->Divide(phiRecF,phiAllF,1,1,"B");
  phiRecF_c->SetName("phiRecEffF");

  phiRecF_c->GetXaxis()->SetTitle(fEffHisto->GetAxis(1)->GetTitle());
  phiRecF_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (findable)",fEffHisto->GetAxis(1)->GetTitle());
  phiRecF_c->SetTitle(title);

  aFolderObj->Add(phiRecF_c);
  }
  else {
  // 
  Float_t minEta=-1.5, maxEta=1.5;
  Float_t minR=0.0, maxR=150.0;
  Float_t minPt=0.15, maxPt=100.0;

  // mother eta range
  fEffSecHisto->GetAxis(8)->SetRangeUser(minEta,maxEta);

  // particle creation radius range 
  fEffSecHisto->GetAxis(6)->SetRangeUser(minR,maxR);

  //
  fEffSecHisto->GetAxis(0)->SetRangeUser(minEta,maxEta);
  fEffSecHisto->GetAxis(2)->SetRangeUser(minPt,maxPt);

  // rec efficiency vs pt
  TH1D *ptAll = fEffSecHisto->Projection(2);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.);  // reconstructed 
  TH1D *ptRec = fEffSecHisto->Projection(2);
  TH1D *ptRec_c = (TH1D*)ptRec->Clone();
  ptRec_c->Divide(ptRec,ptAll,1,1,"B");
  ptRec_c->SetName("ptRecEff");

  ptRec_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(2)->GetTitle());
  ptRec_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency",fEffSecHisto->GetAxis(2)->GetTitle());
  ptRec_c->SetTitle(title);

  ptRec_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRec_c);

  // rec efficiency vs pid vs pt
  
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(0.,0.); // electrons

  TH1D *ptAllEle = fEffSecHisto->Projection(2);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecEle = fEffSecHisto->Projection(2);
  TH1D *ptRecEle_c = (TH1D*)ptRecEle->Clone();
  ptRecEle_c->Divide(ptRecEle,ptAllEle,1,1,"B");
  ptRecEle_c->SetName("ptRecEffEle");

  ptRecEle_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecEle_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (electrons)",fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecEle_c->SetTitle(title);

  ptRecEle_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecEle_c);

  //

  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions

  TH1D *ptAllPi = fEffSecHisto->Projection(2);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecPi = fEffSecHisto->Projection(2);
  TH1D *ptRecPi_c = (TH1D*)ptRecPi->Clone();
  ptRecPi_c->Divide(ptRecPi,ptAllPi,1,1,"B");
  ptRecPi_c->SetName("ptRecEffPi");

  ptRecPi_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecPi_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (pions)",fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecPi_c->SetTitle(title);

  ptRecPi_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecPi_c);

  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(3.,3.); // kaons
  TH1D *ptAllK = fEffSecHisto->Projection(2);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecK = fEffSecHisto->Projection(2);

  TH1D *ptRecK_c = (TH1D*)ptRecK->Clone();
  ptRecK_c->Divide(ptRecK,ptAllK,1,1,"B");
  ptRecK_c->SetName("ptRecEffK");

  ptRecK_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecK_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (kaons)",fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecK_c->SetTitle(title);


  ptRecK_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecK_c);

  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(4.,4.); // protons
  TH1D *ptAllP = fEffSecHisto->Projection(2);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *ptRecP = fEffSecHisto->Projection(2);
  TH1D *ptRecP_c = (TH1D*)ptRecP->Clone();
  ptRecP_c->Divide(ptRecP,ptAllP,1,1,"B");
  ptRecP_c->SetName("ptRecEffP");

  ptRecP_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecP_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (protons)",fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecP_c->SetTitle(title);

  ptRecP_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecP_c);

  // findable efficiency vs pt

  fEffSecHisto->GetAxis(3)->SetRangeUser(0.,4.); 
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(5)->SetRangeUser(1.,1.); // findable
  TH1D *ptAllF = fEffSecHisto->Projection(2);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.);
  fEffSecHisto->GetAxis(5)->SetRangeUser(1.,1.);

  TH1D *ptRecF = fEffSecHisto->Projection(2); // rec findable
  TH1D *ptRecF_c = (TH1D*)ptRecF->Clone();
  ptRecF_c->Divide(ptRecF,ptAllF,1,1,"B");
  ptRecF_c->SetName("ptRecEffF");

  ptRecF_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecF_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (findable)",fEffSecHisto->GetAxis(2)->GetTitle());
  ptRecF_c->SetTitle(title);

  ptRecF_c->SetBit(TH1::kLogX);
  aFolderObj->Add(ptRecF_c);

  //
  // efficiency vs eta
  //
  fEffSecHisto->GetAxis(2)->SetRangeUser(minPt,maxPt);
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.);   // all
  fEffSecHisto->GetAxis(5)->SetRangeUser(0.,1.);   // all

  TH1D *etaAll = fEffSecHisto->Projection(0);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.);  // reconstructed 
  TH1D *etaRec = fEffSecHisto->Projection(0);
  TH1D *etaRec_c = (TH1D*)etaRec->Clone();
  etaRec_c->Divide(etaRec,etaAll,1,1,"B");
  etaRec_c->SetName("etaRecEff");

  etaRec_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(0)->GetTitle());
  etaRec_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency",fEffSecHisto->GetAxis(0)->GetTitle());
  etaRec_c->SetTitle(title);

  aFolderObj->Add(etaRec_c);

  // rec efficiency vs pid vs eta
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(0.,0.); // electrons

  TH1D *etaAllEle = fEffSecHisto->Projection(0);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecEle = fEffSecHisto->Projection(0);
  TH1D *etaRecEle_c = (TH1D*)etaRecEle->Clone();
  etaRecEle_c->Divide(etaRecEle,etaAllEle,1,1,"B");
  etaRecEle_c->SetName("etaRecEffEle");

  etaRecEle_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecEle_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (electrons)",fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecEle_c->SetTitle(title);

  aFolderObj->Add(etaRecEle_c);

  //
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions

  TH1D *etaAllPi = fEffSecHisto->Projection(0);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecPi = fEffSecHisto->Projection(0);
  TH1D *etaRecPi_c = (TH1D*)etaRecPi->Clone();
  etaRecPi_c->Divide(etaRecPi,etaAllPi,1,1,"B");
  etaRecPi_c->SetName("etaRecEffPi");

  etaRecPi_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecPi_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (pions)",fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecPi_c->SetTitle(title);

  aFolderObj->Add(etaRecPi_c);

  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(3.,3.); // kaons
  TH1D *etaAllK = fEffSecHisto->Projection(0);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecK = fEffSecHisto->Projection(0);

  TH1D *etaRecK_c = (TH1D*)etaRecK->Clone();
  etaRecK_c->Divide(etaRecK,etaAllK,1,1,"B");
  etaRecK_c->SetName("etaRecEffK");

  etaRecK_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecK_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (kaons)",fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecK_c->SetTitle(title);


  aFolderObj->Add(etaRecK_c);

  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(4.,4.); // protons
  TH1D *etaAllP = fEffSecHisto->Projection(0);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *etaRecP = fEffSecHisto->Projection(0);
  TH1D *etaRecP_c = (TH1D*)etaRecP->Clone();
  etaRecP_c->Divide(etaRecP,etaAllP,1,1,"B");
  etaRecP_c->SetName("etaRecEffP");

  etaRecP_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecP_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (protons)",fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecP_c->SetTitle(title);

  aFolderObj->Add(etaRecP_c);

  // findable efficiency vs eta

  fEffSecHisto->GetAxis(3)->SetRangeUser(0.,4.); 
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(5)->SetRangeUser(1.,1.); // findable
  TH1D *etaAllF = fEffSecHisto->Projection(0);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.);
  fEffSecHisto->GetAxis(5)->SetRangeUser(1.,1.);

  TH1D *etaRecF = fEffSecHisto->Projection(0); // rec findable
  TH1D *etaRecF_c = (TH1D*)etaRecF->Clone();
  etaRecF_c->Divide(etaRecF,etaAllF,1,1,"B");
  etaRecF_c->SetName("etaRecEffF");

  etaRecF_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecF_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (findable)",fEffSecHisto->GetAxis(0)->GetTitle());
  etaRecF_c->SetTitle(title);

  aFolderObj->Add(etaRecF_c);

  //
  // efficiency vs phi
  //

  fEffSecHisto->GetAxis(0)->SetRangeUser(minEta,maxEta);
  fEffSecHisto->GetAxis(2)->SetRangeUser(minPt,maxPt);

  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.);   // all
  fEffSecHisto->GetAxis(5)->SetRangeUser(0.,1.);   // all

  TH1D *phiAll = fEffSecHisto->Projection(1);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.);  // reconstructed 
  TH1D *phiRec = fEffSecHisto->Projection(1);
  TH1D *phiRec_c = (TH1D*)phiRec->Clone();
  phiRec_c->Divide(phiRec,phiAll,1,1,"B");
  phiRec_c->SetName("phiRecEff");

  phiRec_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(1)->GetTitle());
  phiRec_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency",fEffSecHisto->GetAxis(1)->GetTitle());
  phiRec_c->SetTitle(title);

  aFolderObj->Add(phiRec_c);

  // rec efficiency vs pid vs phi
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions

  TH1D *phiAllEle = fEffSecHisto->Projection(1);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *phiRecEle = fEffSecHisto->Projection(1);
  TH1D *phiRecEle_c = (TH1D*)phiRecEle->Clone();
  phiRecEle_c->Divide(phiRecEle,phiAllEle,1,1,"B");
  phiRecEle_c->SetName("phiRecEffEle");

  phiRecEle_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecEle_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (electrons)",fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecEle_c->SetTitle(title);

  aFolderObj->Add(phiRecEle_c);

  //
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(2.,2.); // pions

  TH1D *phiAllPi = fEffSecHisto->Projection(1);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *phiRecPi = fEffSecHisto->Projection(1);
  TH1D *phiRecPi_c = (TH1D*)phiRecPi->Clone();
  phiRecPi_c->Divide(phiRecPi,phiAllPi,1,1,"B");
  phiRecPi_c->SetName("phiRecEffPi");

  phiRecPi_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecPi_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (pions)",fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecPi_c->SetTitle(title);

  aFolderObj->Add(phiRecPi_c);

  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(3.,3.); // kaons
  TH1D *phiAllK = fEffSecHisto->Projection(1);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *phiRecK = fEffSecHisto->Projection(1);

  TH1D *phiRecK_c = (TH1D*)phiRecK->Clone();
  phiRecK_c->Divide(phiRecK,phiAllK,1,1,"B");
  phiRecK_c->SetName("phiRecEffK");

  phiRecK_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecK_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (kaons)",fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecK_c->SetTitle(title);


  aFolderObj->Add(phiRecK_c);

  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(3)->SetRangeUser(4.,4.); // protons
  TH1D *phiAllP = fEffSecHisto->Projection(1);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.); // reconstructed
  TH1D *phiRecP = fEffSecHisto->Projection(1);
  TH1D *phiRecP_c = (TH1D*)phiRecP->Clone();
  phiRecP_c->Divide(phiRecP,phiAllP,1,1,"B");
  phiRecP_c->SetName("phiRecEffP");

  phiRecP_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecP_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (protons)",fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecP_c->SetTitle(title);

  aFolderObj->Add(phiRecP_c);

  // findable efficiency vs phi

  fEffSecHisto->GetAxis(3)->SetRangeUser(0.,4.); 
  fEffSecHisto->GetAxis(4)->SetRangeUser(0.,1.); 
  fEffSecHisto->GetAxis(5)->SetRangeUser(1.,1.); // findable
  TH1D *phiAllF = fEffSecHisto->Projection(1);

  fEffSecHisto->GetAxis(4)->SetRangeUser(1.,1.);
  fEffSecHisto->GetAxis(5)->SetRangeUser(1.,1.);

  TH1D *phiRecF = fEffSecHisto->Projection(1); // rec findable
  TH1D *phiRecF_c = (TH1D*)phiRecF->Clone();
  phiRecF_c->Divide(phiRecF,phiAllF,1,1,"B");
  phiRecF_c->SetName("phiRecEffF");

  phiRecF_c->GetXaxis()->SetTitle(fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecF_c->GetYaxis()->SetTitle("efficiency");
  sprintf(title,"%s vs %s","rec. efficiency (findable)",fEffSecHisto->GetAxis(1)->GetTitle());
  phiRecF_c->SetTitle(title);

  aFolderObj->Add(phiRecF_c);
  }

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliPerformanceEff::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceEff * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
     // get name and title from old folder
     name = folder->GetName();  
     title = folder->GetTitle();  

	 // delete old one
     delete folder;

	 // create new one
     newFolder = CreateFolder(name.Data(),title.Data());
     newFolder->SetOwner();

	 // add objects to folder
     while(i < size) {
	   newFolder->Add(array->At(i));
	   i++;
	 }
  }

return newFolder;
}


//_____________________________________________________________________________
TFolder* AliPerformanceEff::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
