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
// Changes by J.Salzwedel 10/30/2014
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/Macros/LoadMyLibs.C");
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
#include <TFile.h>
#include <TTree.h>
#include "THnSparse.h"

// 
#include "AliExternalTrackParam.h"
#include "AliLog.h" 
#include "AliTracker.h" 
#include "AliVEvent.h"
#include "AliMCEvent.h" 
#include "AliMCParticle.h" 
#include "AliHeader.h" 
#include "AliGenEventHeader.h" 
#include "AliPerformanceEff.h" 

using namespace std;

#define TPC_PRIMARY_MIN_PT 0.05

ClassImp(AliPerformanceEff)

bool AliPerformanceEff::fReadNClsTree = false;

//_____________________________________________________________________________
AliPerformanceEff::AliPerformanceEff(TRootIOCtor* b):
  AliPerformanceObject(b),
  // histograms
  fEffHisto(NULL),
  fEffSecHisto(NULL),
  fTrackPtNCls(NULL),
  fTrackNClsFound(NULL),
  // histogram folder 
  fAnalysisFolder(NULL),
  fNEvent(0),
  fNClsTreeFile(NULL),
  fNClsTree(NULL),
  fNClsVec(NULL)
{
  // io constructor
}

//_____________________________________________________________________________
AliPerformanceEff::AliPerformanceEff(const Char_t* name, const Char_t* title, Int_t analysisMode, Bool_t hptGenerator):  AliPerformanceObject(name,title),

  // histograms
  fEffHisto(NULL),
  fEffSecHisto(NULL),
  fTrackPtNCls(NULL),
  fTrackNClsFound(NULL),

  // histogram folder 
  fAnalysisFolder(NULL),
  fNEvent(0),
  fNClsTreeFile(NULL),
  fNClsTree(NULL),
  fNClsVec(NULL)
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
  if(fTrackPtNCls) delete fTrackPtNCls; fTrackPtNCls=0;
  if(fTrackNClsFound) delete fTrackNClsFound; fTrackNClsFound=0;
}

//_____________________________________________________________________________
void AliPerformanceEff::Init()
{
  // Init histograms
  //

  // set pt bins
  Int_t nPtBins = 50;
  Double_t ptMin = 0.015, ptMax = 20.;
  if (GetAnalysisMode() == 0) ptMin = TPC_PRIMARY_MIN_PT;

  Double_t *binsPt = 0;

  if (IsHptGenerator())  { 
        ptMax = 100.;
  } 
  binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);

  /*
  Int_t nPtBins = 31;
  Double_t binsPt[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
  Double_t ptMin = 0., ptMax = 10.;

  if(IsHptGenerator() == kTRUE) {
    nPtBins = 100;
    ptMin = 0.; ptMax = 100.;
  }
  */

  //mceta:mcphi:mcpt:pid:recStatus:findable:charge:clones:fakes
  Int_t binsEffHisto[9]={30,144,nPtBins,5,2,2,3,fgkMaxClones+1,fgkMaxFakes+1};
  Double_t minEffHisto[9]={-1.5,0.,ptMin,0.,0.,0.,-1.5,0,0};
  Double_t maxEffHisto[9]={ 1.5,2.*TMath::Pi(), ptMax,5.,2.,2.,1.5,fgkMaxClones+1,fgkMaxFakes+1};

  fEffHisto = new THnSparseF("fEffHisto","mceta:mcphi:mcpt:pid:recStatus:findable:charge:nclones:nfakes",9,binsEffHisto,minEffHisto,maxEffHisto);
  fEffHisto->SetBinEdges(2,binsPt);

  fEffHisto->GetAxis(0)->SetTitle("#eta_{mc}");
  fEffHisto->GetAxis(1)->SetTitle("#phi_{mc} (rad)");
  fEffHisto->GetAxis(2)->SetTitle("p_{Tmc} (GeV/c)");
  fEffHisto->GetAxis(3)->SetTitle("pid");
  fEffHisto->GetAxis(4)->SetTitle("recStatus");
  fEffHisto->GetAxis(5)->SetTitle("findable");
  fEffHisto->GetAxis(6)->SetTitle("charge");
  fEffHisto->GetAxis(7)->SetTitle("nClones");
  fEffHisto->GetAxis(8)->SetTitle("nFakes");
  fEffHisto->Sumw2();

  //mceta:mcphi:mcpt:pid:recStatus:findable:mcR:mother_phi:mother_eta:charge:clones:fakes
  Int_t binsEffSecHisto[12]={30,60,nPtBins,5,2,2,100,60,30,3,fgkMaxClones+1,fgkMaxFakes+1};
  Double_t minEffSecHisto[12]={-1.5,0.,ptMin,0.,0.,0.,0.,0.,-1.5,-1.5,0,0};
  Double_t maxEffSecHisto[12]={ 1.5,2.*TMath::Pi(), ptMax,5.,2.,2.,200,2.*TMath::Pi(),1.5,1.5,fgkMaxClones+1,fgkMaxFakes+1};

  fEffSecHisto = new THnSparseF("fEffSecHisto","mceta:mcphi:mcpt:pid:recStatus:findable:mcR:mother_phi:mother_eta:charge:nclones:nfakes",12,binsEffSecHisto,minEffSecHisto,maxEffSecHisto);
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
  fEffSecHisto->GetAxis(9)->SetTitle("charge");
  fEffSecHisto->GetAxis(10)->SetTitle("nClones");
  fEffSecHisto->GetAxis(11)->SetTitle("nFakes");
  fEffSecHisto->Sumw2();
  
  // init folder
  fAnalysisFolder = CreateFolder("folderEff","Analysis Efficiency Folder");
  
  if (fReadNClsTree)
  {
    fNClsTreeFile = new TFile("nclTPCperMCtrack.root");
    fNClsTree = (TTree*) fNClsTreeFile->Get("nclPerMCtrack");
    if (fNClsTree == NULL)
    {
      AliDebug(AliLog::kFatal, "Stack not available for event");
    }
    fNClsTree->SetBranchAddress("nclPerTrack", &fNClsVec);
    fNEvent = -1;

    fTrackPtNCls = new TH1D("fTrackPtNCls", "total number of cluster vs track Pt", nPtBins, binsPt);
    fTrackNClsFound = new TH2D("fTrackNClsFound", "found number of clusters vs mc number of clusters", 80, 0., 320, 40, 0., 160);
  }
}

//_____________________________________________________________________________
void AliPerformanceEff::ProcessTPC(AliMCEvent* const mcEvent, AliVEvent *const vEvent)
{
  // Fill TPC only efficiency comparison information 
  if(!vEvent) return;
  if(!mcEvent) return;

  Int_t *labelsRec =  NULL;
  labelsRec =  new Int_t[vEvent->GetNumberOfTracks()];
  if(!labelsRec) 
  {
     Printf("Cannot create labelsRec");
     return;
  }
  for(Int_t i=0;i<vEvent->GetNumberOfTracks();i++) { labelsRec[i] = 0; }

  Int_t *labelsAllRec =  NULL;
  labelsAllRec =  new Int_t[vEvent->GetNumberOfTracks()];
  if(!labelsAllRec) { 
    delete  [] labelsRec;
    Printf("Cannot create labelsAllRec");
    return;
  }
  for(Int_t i=0;i<vEvent->GetNumberOfTracks();i++) { labelsAllRec[i] = 0; }

  // loop over rec. tracks
  AliVTrack *vTrack=0;
  for (Int_t iTrack = 0; iTrack < vEvent->GetNumberOfTracks(); iTrack++) 
  { 
    vTrack = dynamic_cast<AliVTrack*>(vEvent->GetTrack(iTrack));
    if(!vTrack) continue;
    if(vTrack->Charge()==0) continue;

    // if not fUseKinkDaughters don't use tracks with kink index > 0
    if(!fUseKinkDaughters && vTrack->GetKinkIndex(0) > 0) continue;

    //Int_t label = TMath::Abs(vTrack->GetLabel()); 
    Int_t label = vTrack->GetTPCLabel(); //Use TPC-only label for TPC-only efficiency analysis 
    labelsAllRec[iTrack]=label;
    
    // TPC only
    if(IsRecTPC(vTrack) != 0) { 
      labelsRec[iTrack]=label;
    }
    
  }
  
  // 
  // MC histograms for efficiency studies
  //
  Int_t nPart  = mcEvent->GetNumberOfTracks();
  Int_t nPrim  = mcEvent->GetNumberOfPrimaries();
  if (fReadNClsTree)
  {
    fNClsTree->GetEntry(fNEvent);
    if (nPart != fNClsVec->size())
    {
      AliDebug(AliLog::kFatal, "nMCTracks mismatch");
    }
  }
  //Int_t nPart  = stack->GetNprimary();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) {
    if (iMc == 0) continue;		//Cannot distinguish between track or fake track
    TParticle* particle = ((AliMCParticle*)mcEvent->GetTrack(iMc))->Particle();
    if (!particle) continue;
    if (!particle->GetPDG()) continue; 
    if (particle->GetPDG()->Charge() == 0.0) continue;

    Float_t mcpt = particle->Pt();
    if (fReadNClsTree) fTrackPtNCls->Fill(mcpt, (*fNClsVec)[iMc]);

    // physical primary
    Bool_t prim = mcEvent->IsPhysicalPrimary(iMc);
    if(!prim) continue;
    
    // --- check for double filling in stack
    // use only particles with no daughters in the list of primaries
    Int_t nDaughters = 0;// particle->GetNDaughters();
    
    //RS we just need to know if there are prim. daughters
    if (particle->GetFirstDaughter() != -1 && particle->GetFirstDaughter() < nPrim) break;
    // --- check for double filling in stack
    
    /*Bool_t findable = kFALSE;
    for(Int_t iRec=0; iRec<esdEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check findable
      if(iMc == labelsAllRec[iRec]) 
      {
        findable = IsFindable(mcEvent,iMc);
        break;
      }
    }*/
    if (!HasTPCReference(mcEvent, iMc)) continue; //Filter particles that do not touch the TPC
    if (fReadNClsTree && (*fNClsVec)[iMc] == 0) continue;
    Bool_t findable = IsFindable(mcEvent,iMc);

    Bool_t recStatus = kFALSE;
    Int_t nClones = 0, nFakes = 0;
    Int_t nClsRec = 0;
    for(Int_t iRec=0; iRec<vEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check reconstructed
      if(iMc == labelsRec[iRec]) {
        if (recStatus && nClones < fgkMaxClones) nClones++;
        recStatus = kTRUE;
        if (fReadNClsTree) nClsRec = ((AliVTrack*) vEvent->GetTrack(iRec))->GetTPCNcls();
      }
      //In order to relate the fake track to track parameters, we assign it to the best matching ESD track.
      if (labelsRec[iRec] < 0 && -labelsRec[iRec] == iMc && nFakes < fgkMaxFakes) nFakes++;
    }

    // Only 5 charged particle species (e,mu,pi,K,p)
    if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) continue; 


    
    // transform Pdg to Pid
    Int_t pid = TransformToPID(particle);
    
    Float_t mceta =  particle->Eta();
    Float_t mcphi =  particle->Phi();
    if(mcphi<0) mcphi += 2.*TMath::Pi();
    Float_t charge = 0.;
    if (particle->GetPDG()->Charge() < 0)  charge = -1.;    
    else if (particle->GetPDG()->Charge() > 0)  charge = 1.;
    if (fReadNClsTree && recStatus) fTrackNClsFound->Fill((*fNClsVec)[iMc], nClsRec);

    // Fill histograms
    Double_t vEffHisto[9] = {mceta, mcphi, mcpt, static_cast<Double_t>(pid), static_cast<Double_t>(recStatus), static_cast<Double_t>(findable), static_cast<Double_t>(charge), static_cast<Double_t>(nClones), static_cast<Double_t>(nFakes)}; 
    fEffHisto->Fill(vEffHisto);
  }
  if(labelsRec) delete [] labelsRec; labelsRec = 0;
  if(labelsAllRec) delete [] labelsAllRec; labelsAllRec = 0;
}

//_____________________________________________________________________________
void AliPerformanceEff::ProcessTPCSec(AliMCEvent* const mcEvent, AliVEvent *const vEvent)
{
  // Fill TPC only efficiency comparison information for secondaries

  if(!vEvent) return;
  if(!mcEvent) return;

  AliStack *stack = mcEvent->Stack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available for event");
    return;
  }

  Int_t *labelsRecSec =  NULL;
  labelsRecSec =  new Int_t[vEvent->GetNumberOfTracks()];
  if(!labelsRecSec) 
  {
     Printf("Cannot create labelsRecSec");
     return;
  }
  for(Int_t i=0;i<vEvent->GetNumberOfTracks();i++) { labelsRecSec[i] = 0; }

  Int_t *labelsAllRecSec =  NULL;
  labelsAllRecSec =  new Int_t[vEvent->GetNumberOfTracks()];
  if(!labelsAllRecSec) { 
     delete [] labelsRecSec;
     Printf("Cannot create labelsAllRecSec");
     return;
  }
  for(Int_t i=0;i<vEvent->GetNumberOfTracks();i++) { labelsAllRecSec[i] = 0; }

  // loop over rec. tracks
  AliVTrack *vTrack=0;
  Int_t multAll=0, multRec=0;
  for (Int_t iTrack = 0; iTrack < vEvent->GetNumberOfTracks(); iTrack++) 
  {
    vTrack = dynamic_cast<AliVTrack*>(vEvent->GetTrack(iTrack));
    if(!vTrack) continue;

    if(vTrack->Charge()==0) continue;

    // if not fUseKinkDaughters don't use tracks with kink index > 0
    if(!fUseKinkDaughters && vTrack->GetKinkIndex(0) > 0) continue;

    //Int_t label = TMath::Abs(track->GetLabel());
    Int_t label = vTrack->GetTPCLabel(); //Use TPC-only label for TPC-only efficiency analysis
    labelsAllRecSec[multAll]=label;
    multAll++;

    // TPC only
    if(IsRecTPC(vTrack) != 0) {
      labelsRecSec[multRec]=label;
      multRec++;
    }
  }

  // 
  // MC histograms for efficiency studies
  //
  if(mcEvent)  { 
 
    Int_t nPart  = mcEvent->GetNumberOfTracks();
    //Int_t nPart  = stack->GetNprimary();
    if (fReadNClsTree)
    {
      fNClsTree->GetEntry(fNEvent);
      if (nPart != fNClsVec->size())
      {
        AliDebug(AliLog::kFatal, "nMCTracks mismatch");
      }
    }
    for (Int_t iMc = 0; iMc < nPart; ++iMc) 
      {
	if (iMc == 0) continue;		//Cannot distinguish between track or fake track
	TParticle* particle = ((AliMCParticle*)mcEvent->GetTrack(iMc))->Particle();
	if (!particle) continue;
	if (!particle->GetPDG()) continue; 
	if (particle->GetPDG()->Charge() == 0.0) continue;
	
	// physical primary
	Bool_t prim = mcEvent->IsPhysicalPrimary(iMc);
	
	// only secondaries which can be reconstructed at TPC
	if(prim) continue;
	
	//Float_t radius = TMath::Sqrt(particle->Vx()*particle->Vx()+particle->Vy()*particle->Vy()+particle->Vz()*particle->Vz());
	//if(radius > fCutsMC->GetMaxR()) continue;
	
	// only secondary electrons from gamma conversion
	//if( TMath::Abs(particle->GetPdgCode())!=fCutsMC->GetEM() ||   particle->GetUniqueID() != 5) continue;
	
	/*Bool_t findable = kFALSE;
	  for(Int_t iRec=0; iRec<multAll; ++iRec) 
	  {
	  // check findable
	  if(iMc == labelsAllRecSec[iRec]) 
	  {
	  findable = IsFindable(mcEvent,iMc);
	  break;
	  }
	  }*/
	if (!HasTPCReference(mcEvent, iMc)) continue; //Filter particles that do not touch the TPC
	double xyz[3];
	((AliMCParticle*) mcEvent->GetTrack(iMc))->XvYvZv(xyz);
	if ((xyz[2] >= 250 && particle->Pz() > 0) || (xyz[2] <= -250 && particle->Pz() < 0)) continue; //Filter particles created at high Z in the TPC readout electronics and flying away
	if (sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]) - fabs(666. / particle->Pt()) > 250) continue; //will never touch the TPC
	if (fReadNClsTree && (*fNClsVec)[iMc] == 0) continue;
	
	Bool_t findable = IsFindable(mcEvent,iMc);
	
	Bool_t recStatus = kFALSE;
	Int_t nClones = 0, nFakes = 0;
	Int_t nClsRec = 0;
	for(Int_t iRec=0; iRec<multRec; ++iRec) 
	  {
	    // check reconstructed
	    if(iMc == labelsRecSec[iRec])
	      {
		if (recStatus && nClones < fgkMaxClones) nClones++;
		recStatus = kTRUE;
		if (fReadNClsTree) nClsRec = ((AliVTrack*) vEvent->GetTrack(iRec))->GetTPCNcls();
	      }
	    //In order to relate the fake track to track parameters, we assign it to the best matching ESD track.
	    if (labelsRecSec[iRec] < 0 && -labelsRecSec[iRec] == iMc && nFakes < fgkMaxFakes) nFakes++;
	  }
	
	// Only 5 charged particle species (e,mu,pi,K,p)
	if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) continue; 
	
	// transform Pdg to Pid
	Int_t pid = TransformToPID(particle);
	
	Float_t mceta =  particle->Eta();
	Float_t mcphi =  particle->Phi();
	if(mcphi<0) mcphi += 2.*TMath::Pi();
	Float_t mcpt = particle->Pt();
	Float_t mcR = particle->R();
	if (fReadNClsTree && recStatus) fTrackNClsFound->Fill((*fNClsVec)[iMc], nClsRec);
	
	// get info about mother
	Int_t motherLabel = particle->GetMother(0);
	if(motherLabel < 0) continue;
	TParticle *mother = ((AliMCParticle*)mcEvent->GetTrack(motherLabel))->Particle(); 
	if(!mother) continue; 

	Float_t mother_eta = mother->Eta();
	Float_t mother_phi = mother->Phi();
	if(mother_phi<0) mother_phi += 2.*TMath::Pi();
	
	Float_t charge = 0.;
	if (particle->GetPDG()->Charge() < 0)  charge = -1.;    
	else if (particle->GetPDG()->Charge() > 0)  charge = 1.;
	
	// Fill histograms
	Double_t vEffSecHisto[12] = { mceta, mcphi, mcpt, static_cast<Double_t>(pid), static_cast<Double_t>(recStatus), static_cast<Double_t>(findable), mcR, mother_phi, mother_eta, static_cast<Double_t>(charge), static_cast<Double_t>(nClones), static_cast<Double_t>(nFakes) }; 
	fEffSecHisto->Fill(vEffSecHisto);
      }
  }
  
  if(labelsRecSec) delete [] labelsRecSec; labelsRecSec = 0;
  if(labelsAllRecSec) delete [] labelsAllRecSec; labelsAllRecSec = 0;
}
 



//_____________________________________________________________________________
void AliPerformanceEff::ProcessTPCITS(AliMCEvent* const mcEvent, AliVEvent *const vEvent)
{
  // Fill efficiency comparison information

  if(!vEvent) return;
  if(!mcEvent) return;

  Int_t *labelsRecTPCITS =  NULL;
  labelsRecTPCITS =  new Int_t[vEvent->GetNumberOfTracks()];
  if(!labelsRecTPCITS) 
  {
     Printf("Cannot create labelsRecTPCITS");
     return;
  }
  for(Int_t i=0;i<vEvent->GetNumberOfTracks();i++) { labelsRecTPCITS[i] = 0; }

  Int_t *labelsAllRecTPCITS =  NULL;
  labelsAllRecTPCITS =  new Int_t[vEvent->GetNumberOfTracks()];
  if(!labelsAllRecTPCITS) { 
     delete [] labelsRecTPCITS;
     Printf("Cannot create labelsAllRecTPCITS");
     return;
  }
  for(Int_t i=0;i<vEvent->GetNumberOfTracks();i++) { labelsAllRecTPCITS[i] = 0; }

  // loop over rec. tracks
  AliVTrack *vTrack=0;
  for (Int_t iTrack = 0; iTrack < vEvent->GetNumberOfTracks(); iTrack++) 
  { 
    vTrack = dynamic_cast<AliVTrack*>(vEvent->GetTrack(iTrack));
    if(!vTrack) continue;
    if(vTrack->Charge()==0) continue;

    // if not fUseKinkDaughters don't use tracks with kink index > 0
    if(!fUseKinkDaughters && vTrack->GetKinkIndex(0) > 0) continue;

    //Int_t label = TMath::Abs(track->GetLabel()); 
    Int_t label = vTrack->GetLabel();  //Use global label for combined efficiency analysis
    labelsAllRecTPCITS[iTrack]=label;

    // iTPC+ITS
    if(IsRecTPCITS(vTrack) != 0) 
      labelsRecTPCITS[iTrack]=label;
  }

  // 
  // MC histograms for efficiency studies
  //
  //Int_t nPart  = stack->GetNtrack();
  Int_t nPart  = mcEvent->GetNumberOfPrimaries();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
    if (iMc == 0) continue;		//Cannot distinguish between track or fake track
    TParticle* particle =  ((AliMCParticle*)mcEvent->GetTrack(iMc))->Particle();
    if (!particle) continue;
    if (!particle->GetPDG()) continue; 
    if (particle->GetPDG()->Charge() == 0.0) continue;
      
    // physical primary
    Bool_t prim = mcEvent->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    Bool_t findable = IsFindable(mcEvent,iMc);

    Bool_t recStatus = kFALSE;
    Int_t nClones = 0, nFakes = 0;
    for(Int_t iRec=0; iRec<vEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check reconstructed
      if(iMc == labelsRecTPCITS[iRec]) 
      {
        if (recStatus && nClones < fgkMaxClones) nClones++;
	recStatus = kTRUE;
      }
      //In order to relate the fake track to track parameters, we assign it to the best matching reconstructed track.
      if (labelsRecTPCITS[iRec] < 0 && -labelsRecTPCITS[iRec] == iMc && nFakes < fgkMaxFakes) nFakes++;
    }

    // Only 5 charged particle species (e,mu,pi,K,p)
    if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) continue; 

    // transform Pdg to Pid
    Int_t pid = TransformToPID(particle);

    Float_t mceta =  particle->Eta();
    Float_t mcphi =  particle->Phi();
    if(mcphi<0) mcphi += 2.*TMath::Pi();
    Float_t mcpt = particle->Pt();

    Float_t charge = 0.;
    if (particle->GetPDG()->Charge() < 0)  charge = -1.;    
    else if (particle->GetPDG()->Charge() > 0)  charge = 1.;
    
    // Fill histograms
    Double_t vEffHisto[9] = { mceta, mcphi, mcpt, static_cast<Double_t>(pid), static_cast<Double_t>(recStatus), static_cast<Double_t>(findable), static_cast<Double_t>(charge), static_cast<Double_t>(nClones), static_cast<Double_t>(nFakes)}; 
    fEffHisto->Fill(vEffHisto);
  }

  if(labelsRecTPCITS) delete [] labelsRecTPCITS; labelsRecTPCITS = 0;
  if(labelsAllRecTPCITS) delete [] labelsAllRecTPCITS; labelsAllRecTPCITS = 0;
}

//_____________________________________________________________________________
void AliPerformanceEff::ProcessConstrained(AliMCEvent* const mcEvent, AliVEvent *const vEvent)
{
  // Process MC vs reconstructed comparison information for TPC-only tracks
  // constrained to the primary vertex
  if(!vEvent) return;
  if(!mcEvent) return;

  Int_t *labelsRecConstrained =  NULL;
  labelsRecConstrained =  new Int_t[vEvent->GetNumberOfTracks()];
  if(!labelsRecConstrained) 
  {
     Printf("Cannot create labelsRecConstrained");
     return;
  }
  for(Int_t i=0;i<vEvent->GetNumberOfTracks();i++) { labelsRecConstrained[i] = 0; }

  Int_t *labelsAllRecConstrained =  NULL;
  labelsAllRecConstrained =  new Int_t[vEvent->GetNumberOfTracks()];
  if(!labelsAllRecConstrained) { 
     delete [] labelsRecConstrained;
     Printf("Cannot create labelsAllRecConstrained");
     return;
  }
  for(Int_t i=0;i<vEvent->GetNumberOfTracks();i++) { labelsAllRecConstrained[i] = 0; }

  // loop over rec. tracks
  AliVTrack *vTrack=0;
  for (Int_t iTrack = 0; iTrack < vEvent->GetNumberOfTracks(); iTrack++) 
  { 
    vTrack = dynamic_cast<AliVTrack*>(vEvent->GetTrack(iTrack));
    if(!vTrack) continue;

    // if not fUseKinkDaughters don't use tracks with kink index > 0
    if(!fUseKinkDaughters && vTrack->GetKinkIndex(0) > 0) continue;

    //Int_t label = TMath::Abs(track->GetLabel()); 
    Int_t label = vTrack->GetLabel();
    labelsAllRecConstrained[iTrack]=label;

    // Constrained
    if(IsRecConstrained(vTrack) != 0) 
      labelsRecConstrained[iTrack]=label;
  }

  // 
  // MC histograms for efficiency studies
  //
 

  //Int_t nPart  = stack->GetNtrack();
  Int_t nPart  = mcEvent->GetNumberOfPrimaries();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
    if (iMc == 0) continue;		//Cannot distinguish between track or fake track
    TParticle* particle = ((AliMCParticle*)mcEvent->GetTrack(iMc))->Particle();
    if (!particle) continue;
    if (!particle->GetPDG()) continue; 
    if (particle->GetPDG()->Charge() == 0.0) continue;
      
    // physical primary
    Bool_t prim = mcEvent->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    Bool_t findable = IsFindable(mcEvent,iMc);

    Bool_t recStatus = kFALSE;
    Int_t nClones = 0, nFakes = 0;
    for(Int_t iRec=0; iRec<vEvent->GetNumberOfTracks(); ++iRec) 
    {
      // check reconstructed
      if(iMc == labelsRecConstrained[iRec]) 
      {
        if (recStatus && nClones < fgkMaxClones) nClones++;
	recStatus = kTRUE;
      }
      //In order to relate the fake track to track parameters, we assign it to the best matching reconstructed track.
      if (labelsRecConstrained[iRec] < 0 && -labelsRecConstrained[iRec] == iMc && nFakes < fgkMaxFakes) nFakes++;
    }

    // Only 5 charged particle species (e,mu,pi,K,p)
    if (fCutsMC.IsPdgParticle(TMath::Abs(particle->GetPdgCode())) == kFALSE) continue; 

    // transform Pdg to Pid
    Int_t pid = TransformToPID(particle);
    Float_t mceta =  particle->Eta();
    Float_t mcphi =  particle->Phi();
    if(mcphi<0) mcphi += 2.*TMath::Pi();
    Float_t mcpt = particle->Pt();
    Float_t charge = 0.;
    if (particle->GetPDG()->Charge() < 0)  charge = -1.;    
    else if (particle->GetPDG()->Charge() > 0)  charge = 1.;

    // Fill histograms
    Double_t vEffHisto[9] = { mceta, mcphi, mcpt, static_cast<Double_t>(pid), static_cast<Double_t>(recStatus), static_cast<Double_t>(findable), static_cast<Double_t>(charge), static_cast<Double_t>(nClones), static_cast<Double_t>(nFakes) }; 
    fEffHisto->Fill(vEffHisto);
  }

  if(labelsRecConstrained) delete [] labelsRecConstrained; labelsRecConstrained = 0;
  if(labelsAllRecConstrained) delete [] labelsAllRecConstrained; labelsAllRecConstrained = 0;
}

//_____________________________________________________________________________
void AliPerformanceEff::Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vFriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend)
{
  // Process comparison information 
  
  fNEvent++;
  if(!vEvent) 
  {
    Error("Exec","vEvent not available");
    return;
  }
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
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

  // use V friends
  if(bUseVfriend) {
    if(!vFriendEvent) {
      Error("Exec","esdFriend not available");
      return;
    }
  }

  //
  //  Process events
  //
  if(GetAnalysisMode() == 0) ProcessTPC(mcEvent,vEvent);
  else if(GetAnalysisMode() == 1) ProcessTPCITS(mcEvent,vEvent);
  else if(GetAnalysisMode() == 2) ProcessConstrained(mcEvent,vEvent);
  else if(GetAnalysisMode() == 5) ProcessTPCSec(mcEvent,vEvent);
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
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC.GetEM() ) pid = 0; 
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC.GetMuM() ) pid = 1; 
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC.GetPiP() ) pid = 2; 
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC.GetKP() ) pid = 3; 
  if( TMath::Abs(particle->GetPdgCode())==fCutsMC.GetProt() ) pid = 4; 

return pid;
}

//_____________________________________________________________________________
Bool_t AliPerformanceEff::IsFindable(const AliMCEvent *mcEvent, Int_t label) 
{
  //
  // Determine if track is findable
  //
  if (fReadNClsTree)
  {
      return((*fNClsVec)[label] > 70);
  }
  if(!mcEvent) return kFALSE;

  AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(label);
  if(!mcParticle) return kFALSE;

  Int_t counter;
  Float_t tpcTrackLength = mcParticle->GetTPCTrackLength(AliTracker::GetBz(),TPC_PRIMARY_MIN_PT,counter,3.0); 
  //printf("tpcTrackLength %f \n", tpcTrackLength);
  return (tpcTrackLength>fCutsMC.GetMinTrackLength());
}

Bool_t AliPerformanceEff::HasTPCReference(const AliMCEvent *mcEvent, Int_t label) 
{
  //
  // Determine if track is findable
  //
  if(!mcEvent) return kFALSE;

  AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(label);
  if(!mcParticle) return kFALSE;

  for (int i = 0;i < mcParticle->GetNumberOfTrackReferences();i++)
  {
    if (mcParticle->GetTrackReference(i)->DetectorId() == AliTrackReference::kTPC) return(true);
  }
  return(false);
}

//_____________________________________________________________________________
Bool_t AliPerformanceEff::IsRecTPC(AliVTrack *vTrack) 
{
//
// Check whether track is reconstructed in TPC
//
  if(!vTrack) return kFALSE;

  Bool_t recStatus = kFALSE;
  const AliExternalTrackParam *track = vTrack->GetTPCInnerParam();
  if(!track) return recStatus;
  
  if(vTrack->GetTPCNcls()>fCutsRC.GetMinNClustersTPC()) recStatus = kTRUE;

  return recStatus;
}

//_____________________________________________________________________________
Bool_t AliPerformanceEff::IsRecTPCITS(AliVTrack *vTrack) 
{
//
// Check whether track is reconstructed in TPCITS
//
  if(!vTrack) return kFALSE;

  Bool_t recStatus = kFALSE;

  if ((vTrack->GetStatus()&AliVTrack::kTPCrefit)==0) return recStatus; // TPC refit
  if ((vTrack->GetStatus()&AliVTrack::kITSrefit)==0) return recStatus; // ITS refit
  if (vTrack && vTrack->GetITSclusters(0)<fCutsRC.GetMinNClustersITS()) return recStatus;  // min. nb. ITS clusters

  recStatus = kTRUE;
  return recStatus;
}

//_____________________________________________________________________________
Bool_t AliPerformanceEff::IsRecConstrained(AliVTrack *vTrack) 
{
//
// Check whether track is reconstructed in IsRecConstrained.
// Here we are looking for TPC-only tracks constrained to the primary
// vertex.
//

  if(!vTrack) return kFALSE;
  Bool_t recStatus = kFALSE;

  if ((vTrack->GetStatus()&AliVTrack::kTPCrefit)==0) return kFALSE; // TPC refit
  if (vTrack->GetTPCNcls()<fCutsRC.GetMinNClustersTPC()) return kFALSE; // min. nb. TPC clusters
  const AliExternalTrackParam * track = vTrack->GetConstrainedParam();
  if(!track) return recStatus;
  if(vTrack->GetITSclusters(0)<fCutsRC.GetMinNClustersITS()) return kFALSE;  // min. nb. ITS clusters
  recStatus = kTRUE;
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
    if (fTrackPtNCls && entry->fTrackPtNCls) fTrackPtNCls->Add(entry->fTrackPtNCls);
    if (fTrackNClsFound && entry->fTrackNClsFound) fTrackNClsFound->Add(entry->fTrackNClsFound);
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
  if(!aFolderObj) return;
  //  char title[256];

  //
  // efficiency vs pt
  //

  if(GetAnalysisMode() != 5) {

  fEffHisto->GetAxis(0)->SetRangeUser(-0.9,0.89); // eta range
  fEffHisto->GetAxis(2)->SetRangeUser(TPC_PRIMARY_MIN_PT,19.99);   // pt range  , primaries below 0.1 do not exist in TPC

  // rec efficiency vs pt
  fEffHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(2, "ptRecEff", "rec. efficiency", 0));
  aFolderObj->Add(AddHistoEff(2, "ptClone", "clone rate", 1));
  aFolderObj->Add(AddHistoEff(2, "ptFake", "fake rate", 2));

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(2, "ptRecEffNeg", "rec. efficiency neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(2, "ptRecEffPos", "rec. efficiency pos.", 0));

  // rec efficiency vs pid vs pt
  fEffHisto->GetAxis(3)->SetRange(3,3);    // pions
  
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(2, "ptRecEffPi", "rec. efficiency (pions)", 0));
 
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(2, "ptRecEffPiNeg", "rec. efficiency (pions) neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(2, "ptRecEffPiPos", "rec. efficiency (pions) pos.", 0));


  fEffHisto->GetAxis(3)->SetRange(4,4);    // kaons
  
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(2, "ptRecEffK", "rec. efficiency (kaons)", 0));
  
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(2, "ptRecEffKNeg", "rec. efficiency (kaons) neg.", 0));
 
  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(2, "ptRecEffKPos", "rec. efficiency (kaons) pos.", 0));
 
  
  fEffHisto->GetAxis(3)->SetRange(5,5);    // protons

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(2, "ptRecEffP", "rec. efficiency (protons)", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(2, "ptRecEffPNeg", "rec. efficiency (protons) neg.", 0));
 
  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(2, "ptRecEffPPos", "rec. efficiency (protons) pos.", 0));

  // findable efficiency vs pt
  fEffHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset
  fEffHisto->GetAxis(5)->SetRange(2,2); // findable

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(2, "ptRecEffF", "rec. efficiency (findable)", 0));
  aFolderObj->Add(AddHistoEff(2, "ptCloneF", "clone rate (findable)", 1));
  aFolderObj->Add(AddHistoEff(2, "ptFakeF", "fake rate (findable)", 2));


  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(2, "ptRecEffFNeg", "rec. efficiency (findable) neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(2, "ptRecEffFPos", "rec. efficiency (findable) pos.", 0));

  //
  // efficiency vs eta
  //

  fEffHisto->GetAxis(0)->SetRangeUser(-1.5,1.49); // eta range
  fEffHisto->GetAxis(2)->SetRangeUser(0.1,19.99); // pt range
  fEffHisto->GetAxis(5)->SetRange(1,0); //reset

  // rec efficiency vs eta
  fEffHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(0, "etaRecEff", "rec. efficiency", 0));
  aFolderObj->Add(AddHistoEff(0, "etaClone", "clone rate", 1));
  aFolderObj->Add(AddHistoEff(0, "etaFake", "fake rate", 2));


  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(0, "etaRecEffNeg", "rec. efficiency neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(0, "etaRecEffPos", "rec. efficiency pos.", 0));

  // rec efficiency vs pid vs eta
  fEffHisto->GetAxis(3)->SetRange(3,3);    // pions
  
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(0, "etaRecEffPi", "rec. efficiency (pions)", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(0, "etaRecEffPiNeg", "rec. efficiency (pions) neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(0, "etaRecEffPiPos", "rec. efficiency (pions) pos.", 0));


  fEffHisto->GetAxis(3)->SetRange(4,4);    // kaons
  
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(0, "etaRecEffK", "rec. efficiency (kaons)", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(0, "etaRecEffKNeg", "rec. efficiency (kaons) neg.", 0));
  
  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(0, "etaRecEffKPos", "rec. efficiency (kaons) pos.", 0));
  
  
  fEffHisto->GetAxis(3)->SetRange(5,5);    // protons

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(0, "etaRecEffP", "rec. efficiency (protons)", 0));
 
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(0, "etaRecEffPNeg", "rec. efficiency (protons) neg.", 0));
  
  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(0, "etaRecEffPPos", "rec. efficiency (protons) pos.", 0));

  // findable efficiency vs eta
  fEffHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset
  fEffHisto->GetAxis(5)->SetRange(2,2); // findable

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(0, "etaRecEffF", "rec. efficiency (findable)", 0));
  aFolderObj->Add(AddHistoEff(0, "etaCloneF", "clone rate (findable)", 1));
  aFolderObj->Add(AddHistoEff(0, "etaFakeF", "fake rate (findable)", 2));


  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(0, "etaRecEffFNeg", "rec. efficiency (findable) neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(0, "etaRecEffFPos", "rec. efficiency (findable) pos.", 0));

  //
  // efficiency vs phi
  //

  fEffHisto->GetAxis(0)->SetRangeUser(-0.9,0.89); // eta range
  fEffHisto->GetAxis(2)->SetRangeUser(0.1,19.99); // pt range
  fEffHisto->GetAxis(5)->SetRange(1,0); //reset

  // rec efficiency vs phi
  fEffHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(1, "phiRecEff", "rec. efficiency", 0));
  aFolderObj->Add(AddHistoEff(1, "phiClone", "clone rate", 1));
  aFolderObj->Add(AddHistoEff(1, "phiFake", "fake rate", 2));

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(1, "phiRecEffNeg", "rec. efficiency neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(1, "phiRecEffPos", "rec. efficiency pos.", 0));

  // rec efficiency vs pid vs phi
  fEffHisto->GetAxis(3)->SetRange(3,3);    // pions
  
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(1, "phiRecEffPi", "rec. efficiency (pions)", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(1, "phiRecEffPiNeg", "rec. efficiency (pions) neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(1, "phiRecEffPiPos", "rec. efficiency (pions) pos.", 0));


  fEffHisto->GetAxis(3)->SetRange(4,4);    // kaons
  
  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(1, "phiRecEffK", "rec. efficiency (kaons)", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(1, "phiRecEffKNeg", "rec. efficiency (kaons) neg.", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(1, "phiRecEffKPos", "rec. efficiency (kaons) pos.", 0));
 
  
  fEffHisto->GetAxis(3)->SetRange(5,5);    // protons

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(1, "phiRecEffP", "rec. efficiency (protons)", 0));

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(1, "phiRecEffPNeg", "rec. efficiency (protons) neg.", 0));
 
  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(1, "phiRecEffPPos", "rec. efficiency (protons) pos.", 0));

  // findable efficiency vs phi
  fEffHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset
  fEffHisto->GetAxis(5)->SetRange(2,2); // findable

  fEffHisto->GetAxis(6)->SetRangeUser(-2.,2.);   // charge all
  aFolderObj->Add(AddHistoEff(1, "phiRecEffF", "rec. efficiency (findable)", 0));
  aFolderObj->Add(AddHistoEff(1, "phiCloneF", "clone rate (findable)", 1));
  aFolderObj->Add(AddHistoEff(1, "phiFakeF", "fake rate (findable)", 2));


  fEffHisto->GetAxis(6)->SetRangeUser(-2.,0.);   // charge negativ
  aFolderObj->Add(AddHistoEff(1, "phiRecEffFNeg", "rec. efficiency (findable) neg.", 0));
 
  fEffHisto->GetAxis(6)->SetRangeUser(0.,2.);    // charge positiv
  aFolderObj->Add(AddHistoEff(1, "phiRecEffFPos", "rec. efficiency (findable) pos.", 0));
  }
  else {
  // 
  Float_t minEta=-1.5, maxEta=1.49;
  Float_t minR=0.0, maxR=150.0;
  Float_t minPt=15e-3, maxPt=19.99;

  // mother eta range
  fEffSecHisto->GetAxis(8)->SetRangeUser(minEta,maxEta);

  // particle creation radius range 
  fEffSecHisto->GetAxis(6)->SetRangeUser(minR,maxR);

  // particle eta and pt cut
  fEffSecHisto->GetAxis(0)->SetRangeUser(-0.9,0.89);
  fEffSecHisto->GetAxis(2)->SetRangeUser(minPt,maxPt);

  // rec efficiency vs pt

  aFolderObj->Add(AddHistoEff(2, "ptRecEff", "rec. efficiency", 0, 1));
  aFolderObj->Add(AddHistoEff(2, "ptClone", "clone rate", 1, 1));
  aFolderObj->Add(AddHistoEff(2, "ptFake", "fake rate", 2, 1));

  // rec efficiency vs pid vs pt
  
  fEffSecHisto->GetAxis(3)->SetRange(1,1); // electrons
  aFolderObj->Add(AddHistoEff(2, "ptRecEffEle", "rec. efficiency (electrons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(3,3); // pions
  aFolderObj->Add(AddHistoEff(2, "ptRecEffPi", "rec. efficiency (pions)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(4,4); // kaons
  aFolderObj->Add(AddHistoEff(2, "ptRecEffK", "rec. efficiency (kaons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(5,5); // protons
  aFolderObj->Add(AddHistoEff(2, "ptRecEffP", "rec. efficiency (protons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset
  
  // findable efficiency vs pt

  fEffSecHisto->GetAxis(5)->SetRange(2,2); // findable
  aFolderObj->Add(AddHistoEff(2, "ptRecEffF", "rec. efficiency (findable)", 0, 1));
  aFolderObj->Add(AddHistoEff(2, "ptCloneF", "clone rate (findable)", 1, 1));
  aFolderObj->Add(AddHistoEff(2, "ptFakeF", "fake rate (findable)", 2, 1));

  //
  // efficiency vs eta
  //
  fEffSecHisto->GetAxis(0)->SetRangeUser(minEta,maxEta);
  fEffSecHisto->GetAxis(2)->SetRangeUser(0.1,maxPt);
  fEffSecHisto->GetAxis(5)->SetRange(1,0); //reset

  aFolderObj->Add(AddHistoEff(0, "etaRecEff", "rec. efficiency", 0, 1));
  aFolderObj->Add(AddHistoEff(0, "etaClone", "clone rate", 1, 1));
  aFolderObj->Add(AddHistoEff(0, "etaFake", "fake rate", 2, 1));

  // rec efficiency vs pid vs eta
  fEffSecHisto->GetAxis(3)->SetRange(1,1); // electrons
  aFolderObj->Add(AddHistoEff(0, "etaRecEffEle", "rec. efficiency (electrons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(3,3); // pions
  aFolderObj->Add(AddHistoEff(0, "etaRecEffPi", "rec. efficiency (pions)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(4,4); // kaons
  aFolderObj->Add(AddHistoEff(0, "etaRecEffK", "rec. efficiency (kaons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(5,5); // protons
  aFolderObj->Add(AddHistoEff(0, "etaRecEffP", "rec. efficiency (protons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset

  // findable efficiency vs eta

  fEffSecHisto->GetAxis(5)->SetRange(2,2); // findable
  aFolderObj->Add(AddHistoEff(0, "etaRecEffF", "rec. efficiency (findable)", 0, 1));
  aFolderObj->Add(AddHistoEff(0, "etaCloneF", "clone rate (findable)", 1, 1));
  aFolderObj->Add(AddHistoEff(0, "etaFakeF", "fake rate (findable)", 2, 1));

  //
  // efficiency vs phi
  //

  fEffSecHisto->GetAxis(0)->SetRangeUser(-0.9,0.89);
  fEffSecHisto->GetAxis(2)->SetRangeUser(0.1,maxPt);
  fEffSecHisto->GetAxis(5)->SetRange(1,0); //reset

  aFolderObj->Add(AddHistoEff(1, "phiRecEff", "rec. efficiency", 0, 1));
  aFolderObj->Add(AddHistoEff(1, "phiClone", "clone rate", 1, 1));
  aFolderObj->Add(AddHistoEff(1, "phiFake", "fake rate", 2, 1));

  // rec efficiency vs pid vs phi
  fEffSecHisto->GetAxis(3)->SetRange(1,1);
  aFolderObj->Add(AddHistoEff(1, "phiRecEffEle", "rec. efficiency (electrons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(3,3); // pions
  aFolderObj->Add(AddHistoEff(1, "phiRecEffPi", "rec. efficiency (pions)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(4,4); // kaons
  aFolderObj->Add(AddHistoEff(1, "phiRecEffK", "rec. efficiency (kaons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRange(5,5); // protons
  aFolderObj->Add(AddHistoEff(1, "phiRecEffP", "rec. efficiency (protons)", 0, 1));

  fEffSecHisto->GetAxis(3)->SetRangeUser(0.,5.); //reset 

  // findable efficiency vs phi

  fEffSecHisto->GetAxis(5)->SetRange(2,2); // findable
  aFolderObj->Add(AddHistoEff(1, "phiRecEffF", "rec. efficiency (findable)", 0, 1));
  aFolderObj->Add(AddHistoEff(1, "phiCloneF", "clone rate (findable)", 1, 1));
  aFolderObj->Add(AddHistoEff(1, "phiFakeF", "fake rate (findable)", 2, 1));
  }

  THnSparseF* h = (GetAnalysisMode() != 5) ? fEffHisto : fEffSecHisto;
  for (Int_t i = 0;i < h->GetNdimensions();i++)
  {
	  h->GetAxis(i)->SetRange(1,0);				//Reset Range
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

TH1D* AliPerformanceEff::WeightedProjection(THnSparseF* src, Int_t axis, Int_t nWeights, Int_t* weightCoords)
{
    THnSparseF* tmp = (THnSparseF*) src->Clone();
    Int_t i;
    for (i = 0;i < tmp->GetNbins();i++)
    {
        Int_t coords[12];
        tmp->GetBinContent(i, coords);
        Int_t weight = 0, j;
        for (j = 0;j < nWeights;j++)
        {
			//The coordinates range from 1 to maxClones / maxFakes + 1, so we have to subtract one
            weight += coords[weightCoords[j]] - 1;
        }
        tmp->SetBinContent(i, weight);
    }
    
    TH1D* ret = tmp->Projection(axis);
    delete tmp;
    return(ret);
}

//_____________________________________________________________________________
TH1D* AliPerformanceEff::AddHistoEff(Int_t axis, const Char_t *name, const Char_t* vsTitle,
	const Int_t type, const Int_t secondary) {
  // Create and add rec efficiency vs pt, eta, phi
  
  char title[256];

  TH1D *recc = NULL;

  THnSparseF* EffHisto = secondary ? fEffSecHisto : fEffHisto;

  Int_t axis_clone = secondary ? 10 : 7;
  Int_t axis_fake = secondary ? 11 : 8;
  Int_t axis_all[3] = {4, axis_clone, axis_fake};


  if (type == 0) // Efficiency
  {
      EffHisto->GetAxis(4)->SetRange(1.,2.);    // all
      TH1D *all = EffHisto->Projection(axis);

      EffHisto->GetAxis(4)->SetRange(2.,2.);    // reconstructed 
      TH1D *rec = EffHisto->Projection(axis);
      recc = (TH1D*)rec->Clone();

      if(recc) 
      {
        recc->Divide(rec,all,1,1,"B");
        recc->GetYaxis()->SetTitle("efficiency");
      }
  }
  else if (type == 1) // Clone Rate
  {
      EffHisto->GetAxis(4)->SetRange(2.,2.);    // reconstructed
      TH1D *all = WeightedProjection(EffHisto, axis, 3, axis_all);

      EffHisto->GetAxis(4)->SetRange(2.,2.);    // reconstructed
      TH1D *clone = WeightedProjection(EffHisto, axis, 1, &axis_clone);
      recc = (TH1D*) clone->Clone();

      if(recc) 
      {
        recc->Divide(clone,all,1,1,"B");
        recc->GetYaxis()->SetTitle("clone rate");
      }
  }
  else if (type == 2) // Fake Rate
  {
      EffHisto->GetAxis(4)->SetRange(2.,2.);    // reconstructed
      TH1D *all = WeightedProjection(EffHisto, axis, 3, axis_all);

      EffHisto->GetAxis(4)->SetRange(2.,2.);    // reconstructed 
      TH1D *fake = WeightedProjection(EffHisto, axis, 1, &axis_fake);
      recc = (TH1D*) fake->Clone();

      if(recc) 
      {
        recc->Divide(fake,all,1,1,"B");
        recc->GetYaxis()->SetTitle("fake rate");
      }
  }

  EffHisto->GetAxis(4)->SetRange(1,0);				//Reset Range

  if (recc) { // Coverity fix
    recc->SetName(name);
    
    recc->GetXaxis()->SetTitle(fEffHisto->GetAxis(axis)->GetTitle());

    snprintf(title,256,"%s vs %s",vsTitle, fEffHisto->GetAxis(axis)->GetTitle());  
    recc->SetTitle(title);

    if (axis == 2 ) recc->SetBit(TH1::kLogX);
  }
  return recc;
}
