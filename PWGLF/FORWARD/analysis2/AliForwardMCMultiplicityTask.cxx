// 
// Calculate the multiplicity in the forward regions event-by-event 
// 
// Inputs: 
//   - AliESDEvent 
//   - Kinematics
//   - Track references
//
// Outputs: 
//   - AliAODForwardMult 
// 
// Histograms 
//   
// Corrections used 
// 
#include "AliForwardMCMultiplicityTask.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliForwardCorrectionManager.h"
#include "AliAnalysisManager.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TROOT.h>

//====================================================================
AliForwardMCMultiplicityTask::AliForwardMCMultiplicityTask()
  : AliForwardMultiplicityBase(),
    fHData(0),
    fESDFMD(),
    fHistos(),
    fAODFMD(),
    fAODEP(),
    fMCESDFMD(),
    fMCHistos(),
    fMCAODFMD(),
    fRingSums(),
    fMCRingSums(),
    fPrimary(0),
    fEventInspector(),
    fSharingFilter(),
    fDensityCalculator(),
    fCorrections(),
    fHistCollector(),
    fEventPlaneFinder(),
    fList(0),
    fListVertexBins(0)	    
	
{
  // 
  // Constructor
  //
}

//____________________________________________________________________
AliForwardMCMultiplicityTask::AliForwardMCMultiplicityTask(const char* name)
  : AliForwardMultiplicityBase(name), 
    fHData(0),
    fESDFMD(),
    fHistos(),
    fAODFMD(kFALSE),
    fAODEP(kFALSE),
    fMCESDFMD(),
    fMCHistos(),
    fMCAODFMD(kTRUE),
    fRingSums(),
    fMCRingSums(),
    fPrimary(0),
    fEventInspector("event"),
    fSharingFilter("sharing"), 
    fDensityCalculator("density"),
    fCorrections("corrections"),
    fHistCollector("collector"),
    fEventPlaneFinder("eventplane"),
    fList(0),
    fListVertexBins(0)	  	
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________________
AliForwardMCMultiplicityTask::AliForwardMCMultiplicityTask(const AliForwardMCMultiplicityTask& o)
  : AliForwardMultiplicityBase(o),
    fHData(o.fHData),
    fESDFMD(o.fESDFMD),
    fHistos(o.fHistos),
    fAODFMD(o.fAODFMD),
    fAODEP(o.fAODEP),
    fMCESDFMD(o.fMCESDFMD),
    fMCHistos(o.fMCHistos),
    fMCAODFMD(o.fMCAODFMD),
    fRingSums(o.fRingSums),
    fMCRingSums(o.fMCRingSums),
    fPrimary(o.fPrimary),
    fEventInspector(o.fEventInspector),
    fSharingFilter(o.fSharingFilter),
    fDensityCalculator(o.fDensityCalculator),
    fCorrections(o.fCorrections),
    fHistCollector(o.fHistCollector),
    fEventPlaneFinder(o.fEventPlaneFinder),
    fList(o.fList),	
    fListVertexBins(o.fListVertexBins)	  	 
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________________
AliForwardMCMultiplicityTask&
AliForwardMCMultiplicityTask::operator=(const AliForwardMCMultiplicityTask& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
  if (&o == this) return *this;
  AliForwardMultiplicityBase::operator=(o);

  fHData             = o.fHData;
  fEventInspector    = o.fEventInspector;
  fSharingFilter     = o.fSharingFilter;
  fDensityCalculator = o.fDensityCalculator;
  fCorrections       = o.fCorrections;
  fHistCollector     = o.fHistCollector;
  fEventPlaneFinder  = o.fEventPlaneFinder;
  fHistos            = o.fHistos;
  fAODFMD            = o.fAODFMD;
  fAODEP             = o.fAODEP;
  fMCHistos          = o.fMCHistos;
  fMCAODFMD          = o.fMCAODFMD;
  fRingSums          = o.fRingSums;
  fMCRingSums        = o.fMCRingSums;
  fPrimary           = o.fPrimary;
  fList              = o.fList;
  fListVertexBins    = o.fListVertexBins; 	
  return *this;
}

//____________________________________________________________________
void
AliForwardMCMultiplicityTask::SetDebug(Int_t dbg)
{
  // 
  // Set debug level 
  // 
  // Parameters:
  //    dbg debug level
  //
  fEventInspector.SetDebug(dbg);
  fSharingFilter.SetDebug(dbg);
  fDensityCalculator.SetDebug(dbg);
  fCorrections.SetDebug(dbg);
  fHistCollector.SetDebug(dbg);
  fEventPlaneFinder.SetDebug(dbg);
}
//____________________________________________________________________
void
AliForwardMCMultiplicityTask::SetOnlyPrimary(Bool_t use)
{
  fSharingFilter.GetTrackDensity().SetUseOnlyPrimary(use);
  fCorrections.SetSecondaryForMC(!use);
}

//____________________________________________________________________
Bool_t
AliForwardMCMultiplicityTask::SetupForData()
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  const TAxis* pe = 0;
  const TAxis* pv = 0;

  if (!ReadCorrections(pe,pv,true)) return false;

  fHistos.Init(*pe);
  fAODFMD.Init(*pe);
  fAODEP.Init(*pe);
  fMCHistos.Init(*pe);
  fMCAODFMD.Init(*pe);
  fRingSums.Init(*pe);
  fMCRingSums.Init(*pe);

  fHData = static_cast<TH2D*>(fAODFMD.GetHistogram().Clone("d2Ndetadphi"));
  fHData->SetStats(0);
  fHData->SetDirectory(0);

  fList->Add(fHData);

  TList* rings = new TList;
  rings->SetName("ringSums");
  rings->SetOwner();
  fList->Add(rings);

  rings->Add(fRingSums.Get(1, 'I'));
  rings->Add(fRingSums.Get(2, 'I'));
  rings->Add(fRingSums.Get(2, 'O'));
  rings->Add(fRingSums.Get(3, 'I'));
  rings->Add(fRingSums.Get(3, 'O'));
  fRingSums.Get(1, 'I')->SetMarkerColor(AliForwardUtil::RingColor(1, 'I'));
  fRingSums.Get(2, 'I')->SetMarkerColor(AliForwardUtil::RingColor(2, 'I'));
  fRingSums.Get(2, 'O')->SetMarkerColor(AliForwardUtil::RingColor(2, 'O'));
  fRingSums.Get(3, 'I')->SetMarkerColor(AliForwardUtil::RingColor(3, 'I'));
  fRingSums.Get(3, 'O')->SetMarkerColor(AliForwardUtil::RingColor(3, 'O'));

  TList* mcRings = new TList;
  mcRings->SetName("mcRingSums");
  mcRings->SetOwner();
  fList->Add(mcRings);

  mcRings->Add(fMCRingSums.Get(1, 'I'));
  mcRings->Add(fMCRingSums.Get(2, 'I'));
  mcRings->Add(fMCRingSums.Get(2, 'O'));
  mcRings->Add(fMCRingSums.Get(3, 'I'));
  mcRings->Add(fMCRingSums.Get(3, 'O'));
  fMCRingSums.Get(1, 'I')->SetMarkerColor(AliForwardUtil::RingColor(1, 'I'));
  fMCRingSums.Get(2, 'I')->SetMarkerColor(AliForwardUtil::RingColor(2, 'I'));
  fMCRingSums.Get(2, 'O')->SetMarkerColor(AliForwardUtil::RingColor(2, 'O'));
  fMCRingSums.Get(3, 'I')->SetMarkerColor(AliForwardUtil::RingColor(3, 'I'));
  fMCRingSums.Get(3, 'O')->SetMarkerColor(AliForwardUtil::RingColor(3, 'O'));


  for(int i=1;i<=pv->GetNbins();i++)	
  {
	TString nametmp=Form("vtxbin%03d",i);
	//TList* lbin= new TList();
	//lbin->SetName(nametmp.Data());
	//lbin->SetOwner();
	//fListVertexBins->Add(lbin);
	AliForwardUtil::Histos* bin=new AliForwardUtil::Histos();
	bin->Init(*pe);
	bin->Get(1, 'I')->SetName(Form("%s%s",bin->Get(1, 'I')->GetName(),nametmp.Data()));
	bin->Get(2, 'I')->SetName(Form("%s%s",bin->Get(2, 'I')->GetName(),nametmp.Data()));
	bin->Get(2, 'O')->SetName(Form("%s%s",bin->Get(2, 'O')->GetName(),nametmp.Data())); 
	bin->Get(3, 'I')->SetName(Form("%s%s",bin->Get(3, 'I')->GetName(),nametmp.Data()));
	bin->Get(3, 'O')->SetName(Form("%s%s",bin->Get(3, 'O')->GetName(),nametmp.Data()));
	fList->Add(bin->Get(1, 'I'));
	fList->Add(bin->Get(2, 'I'));
	fList->Add(bin->Get(2, 'O'));
	fList->Add(bin->Get(3, 'I'));
	fList->Add(bin->Get(3, 'O'));
	fListVertexBins->Add(bin);

}
  fEventInspector.SetupForData(*pv);
  fSharingFilter.SetupForData(*pe);
  fDensityCalculator.SetupForData(*pe);
  fCorrections.SetupForData(*pe);
  fHistCollector.SetupForData(*pv,*pe);
  fEventPlaneFinder.SetupForData(*pe);

  this->Print();

  return true;
}

//____________________________________________________________________
void
AliForwardMCMultiplicityTask::UserCreateOutputObjects()
{
  // 
  // Create output objects 
  // 
  //
  fList = new TList;
  fList->SetOwner();

  fListVertexBins=new TList();
  fListVertexBins->SetOwner();		
  //fList->Add(fListVertexBins);
	
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (ah) 
  {  
    //AliFatal("No AOD output handler set in analysis manager");

  	TObject* obj = &fAODFMD;
 	 ah->AddBranch("AliAODForwardMult", &obj);

  	TObject* mcobj = &fMCAODFMD;
  	ah->AddBranch("AliAODForwardMult", &mcobj);

  	TObject* epobj = &fAODEP;
  	ah->AddBranch("AliAODForwardEP", &epobj);
  }
  fPrimary = new TH2D("primary", "MC Primaries", 
		      200, -4, 6, 20, 0, 2*TMath::Pi());
  fPrimary->SetXTitle("#eta");
  fPrimary->SetYTitle("#varphi [radians]");
  fPrimary->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  fPrimary->Sumw2();
  fPrimary->SetStats(0);
  fPrimary->SetDirectory(0);
  if(ah)	
 	 ah->AddBranch("TH2D", &fPrimary);

  fEventInspector.CreateOutputObjects(fList);
  fSharingFilter.CreateOutputObjects(fList);
  fDensityCalculator.CreateOutputObjects(fList);
  fCorrections.CreateOutputObjects(fList);
  fHistCollector.CreateOutputObjects(fList);
  fEventPlaneFinder.CreateOutputObjects(fList);


 	

  PostData(1, fList);
}
//____________________________________________________________________
void
AliForwardMCMultiplicityTask::UserExec(Option_t*)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  

  // Read production details 
  if (fFirstEvent) 
    fEventInspector.ReadProductionDetails(MCEvent());
    
  // Get the input data 
  AliESDEvent* esd     = GetESDEvent();
  AliMCEvent*  mcEvent = MCEvent();
  if (!esd || !mcEvent) return;

  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  fAODFMD.Clear();
  fAODEP.Clear();
  fMCHistos.Clear();
  fMCESDFMD.Clear();
  fMCAODFMD.Clear();
  fPrimary->Reset();

  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip(1024, 1024, 0);
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(esd, triggers, lowFlux, 
					       ivz, ip, cent, nClusters);
  UShort_t ivzMC    = 0;
  Double_t vzMC     = 0;
  Double_t phiR     = 0;
  Double_t b        = 0;
  Double_t cMC      = 0;
  Int_t    npart    = 0;
  Int_t    nbin     = 0;
  // UInt_t   foundMC  = 
  fEventInspector.ProcessMC(mcEvent, triggers, ivzMC, vzMC, b, cMC,
			    npart, nbin, phiR);
  fEventInspector.CompareResults(ip.Z(), vzMC, cent, cMC, b, npart, nbin);
  
  //Store all events
  MarkEventForStore();
  
  Bool_t isAccepted = true;
  if (found & AliFMDEventInspector::kNoEvent)    isAccepted = false; // return;
  if (found & AliFMDEventInspector::kNoTriggers) isAccepted = false; // return;
  //MarkEventForStore();
  // Always set the B trigger - each MC event _is_ a collision 
  triggers |= AliAODForwardMult::kB;
  // Set trigger bits, and mark this event for storage 
  fAODFMD.SetTriggerBits(triggers);
  fAODFMD.SetSNN(fEventInspector.GetEnergy());
  fAODFMD.SetSystem(fEventInspector.GetCollisionSystem());
  fAODFMD.SetCentrality(cent);
  fAODFMD.SetNClusters(nClusters);

  fMCAODFMD.SetTriggerBits(triggers);
  fMCAODFMD.SetSNN(fEventInspector.GetEnergy());
  fMCAODFMD.SetSystem(fEventInspector.GetCollisionSystem());
  fMCAODFMD.SetCentrality(cent);
  fMCAODFMD.SetNClusters(nClusters);
  
  //All events should be stored - HHD
  //if (isAccepted) MarkEventForStore();

  // Disable this check on SPD - will bias data 
  // if (found & AliFMDEventInspector::kNoSPD)  isAccepted = false; // return;
  if (found & AliFMDEventInspector::kNoFMD)     isAccepted = false; // return;
  if (found & AliFMDEventInspector::kNoVertex)  isAccepted = false; // return;

  if (isAccepted) {
    fAODFMD.SetIpZ(ip.Z());
    fMCAODFMD.SetIpZ(ip.Z());
  }
  if (found & AliFMDEventInspector::kBadVertex) isAccepted = false; // return;

  // We we do not want to use low flux specific code, we disable it here. 
  if (!fEnableLowFlux) lowFlux = false;
  
  

  // Get FMD data 
  AliESDFMD*  esdFMD  = esd->GetFMDData();

  // Apply the sharing filter (or hit merging or clustering if you like)
  if (isAccepted && !fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD,ip.Z())) { 
    AliWarning("Sharing filter failed!");
    return;
  }
  if (!fSharingFilter.FilterMC(*esdFMD, *mcEvent, ip.Z(),fMCESDFMD,fPrimary)) { 
    AliWarning("MC Sharing filter failed!");
    return;
  }
  if (!isAccepted) return; // Exit on MC event w/o trigger, vertex, data
  // HHD if (!isAccepted) return; // Exit on MC event w/o trigger, vertex, data
  
  //MarkEventForStore();
  fSharingFilter.CompareResults(fESDFMD, fMCESDFMD);

  // Calculate the inclusive charged particle density 
  if (!fDensityCalculator.Calculate(fESDFMD, fHistos, lowFlux, cent, ip)) { 
    AliWarning("Density calculator failed!");
    return;
  }
  if (!fDensityCalculator.CalculateMC(fMCESDFMD, fMCHistos)) { 
    AliWarning("MC Density calculator failed!");
    return;
  }
  fDensityCalculator.CompareResults(fHistos, fMCHistos);
  
  if (fEventInspector.GetCollisionSystem() == AliFMDEventInspector::kPbPb) {
    if (!fEventPlaneFinder.FindEventplane(esd, fAODEP, &(fAODFMD.GetHistogram()) , &fHistos))
      AliWarning("Eventplane finder failed!");
  }

  // Do the secondary and other corrections. 
  if (!fCorrections.Correct(fHistos, ivz)) { 
    AliWarning("Corrections failed");
    return;
  }
  if (!fCorrections.CorrectMC(fMCHistos, ivz)) { 
    AliWarning("MC Corrections failed");
    return;
  }
  fCorrections.CompareResults(fHistos, fMCHistos);
    
  if (!fHistCollector.Collect(fHistos, fRingSums, 
			      ivz, fAODFMD.GetHistogram(),0x0,-1,fListVertexBins)) {
    AliWarning("Histogram collector failed");
    return;
  }
  if (!fHistCollector.Collect(fMCHistos, fMCRingSums, 
			      ivz, fMCAODFMD.GetHistogram())) {
    AliWarning("MC Histogram collector failed");
    return;
  }

  if (fAODFMD.IsTriggerBits(AliAODForwardMult::kInel))
    fHData->Add(&(fAODFMD.GetHistogram()));

  PostData(1, fList);
}

//____________________________________________________________________
void
AliForwardMCMultiplicityTask::Terminate(Option_t*)
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  TList* list = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError(Form("No output list defined (%p)", GetOutputData(1)));
    if (GetOutputData(1)) GetOutputData(1)->Print();
    return;
  }

  // Output list 
  TList* output = new TList;
  output->SetOwner();
  output->SetName(Form("%sResults", GetName()));

  Double_t nTr = 0, nTrVtx = 0, nAcc = 0;
  MakeSimpledNdeta(list, output, nTr, nTrVtx, nAcc);
  MakeRingdNdeta(list, "ringSums", output, "ringResults");
  MakeRingdNdeta(list, "mcRingSums", output, "mcRingResults", 24);

  fSharingFilter.Terminate(list,output,Int_t(nTr));
  fDensityCalculator.Terminate(list,output,Int_t(nTrVtx));
  fCorrections.Terminate(list,output,Int_t(nTrVtx));

  PostData(2, output);
}


//
// EOF
//
