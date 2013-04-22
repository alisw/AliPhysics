// 
// Calculate the multiplicity in the forward regions event-by-event 
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODForwardMult 
// 
// Histograms 
//   
// Corrections used 
//
#include "AliForwardMultiplicityTask.h"
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
#include <TH3D.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TROOT.h>


//====================================================================
AliForwardMultiplicityTask::AliForwardMultiplicityTask()
  : AliForwardMultiplicityBase(),
    fHData(0),
    fESDFMD(),
    fHistos(),
    fAODFMD(),
    fAODEP(),
    fRingSums(),
    fEventInspector(),
    fSharingFilter(),
    fDensityCalculator(),
    fCorrections(),
    fHistCollector(),
    fEventPlaneFinder(),
    fFMD1icent(0),
    fFMD2icent(0),
    fFMD2ocent(0),
    fFMD3icent(0),
    fFMD3ocent(0),
    fList(0),	
    fListVertexBins(0)

{
  // 
  // Constructor
  //
  DGUARD(fDebug, 3,"Default CTOR of AliForwardMultiplicityTask");
}

//____________________________________________________________________
AliForwardMultiplicityTask::AliForwardMultiplicityTask(const char* name)
  : AliForwardMultiplicityBase(name),
    fHData(0),
    fESDFMD(),
    fHistos(),
    fAODFMD(false),
    fAODEP(false),
    fRingSums(),
    fEventInspector("event"),
    fSharingFilter("sharing"), 
    fDensityCalculator("density"),
    fCorrections("corrections"),
    fHistCollector("collector"),
    fEventPlaneFinder("eventplane"),
    fFMD1icent(0),
    fFMD2icent(0),
    fFMD2ocent(0),
    fFMD3icent(0),
    fFMD3ocent(0),
    fList(0),
    fListVertexBins(0)	


{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DGUARD(fDebug, 3,"named CTOR of AliForwardMultiplicityTask: %s", name);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________________
AliForwardMultiplicityTask::AliForwardMultiplicityTask(const AliForwardMultiplicityTask& o)
  : AliForwardMultiplicityBase(o),
    fHData(o.fHData),
    fESDFMD(o.fESDFMD),
    fHistos(o.fHistos),
    fAODFMD(o.fAODFMD),
    fAODEP(o.fAODEP),
    fRingSums(o.fRingSums),
    fEventInspector(o.fEventInspector),
    fSharingFilter(o.fSharingFilter),
    fDensityCalculator(o.fDensityCalculator),
    fCorrections(o.fCorrections),
    fHistCollector(o.fHistCollector),
    fEventPlaneFinder(o.fEventPlaneFinder),
     fFMD1icent(o.fFMD1icent),
    fFMD2icent(o.fFMD2icent),
    fFMD2ocent(o.fFMD2ocent),
    fFMD3icent(o.fFMD3icent),
    fFMD3ocent(o.fFMD3ocent),
    fList(o.fList),
    fListVertexBins(o.fListVertexBins)	

{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DGUARD(fDebug, 3,"Copy CTOR of AliForwardMultiplicityTask");
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________________
AliForwardMultiplicityTask&
AliForwardMultiplicityTask::operator=(const AliForwardMultiplicityTask& o)
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
  DGUARD(fDebug,3,"Assignment to AliForwardMultiplicityTask");
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
  fRingSums          = o.fRingSums;
  fFMD1icent	     = o.fFMD1icent;
  fFMD2icent	     = o.fFMD2icent;
  fFMD2ocent	     = o.fFMD2ocent;
  fFMD3icent	     = o.fFMD3icent;
  fFMD3ocent	     = o.fFMD3ocent;
  fList              = o.fList;
  fListVertexBins    =o.fListVertexBins;
  
  return *this;
}

//____________________________________________________________________
void
AliForwardMultiplicityTask::SetDebug(Int_t dbg)
{
  // 
  // Set debug level 
  // 
  // Parameters:
  //    dbg Debug level
  //
  fEventInspector.SetDebug(dbg);
  fSharingFilter.SetDebug(dbg);
  fDensityCalculator.SetDebug(dbg);
  fCorrections.SetDebug(dbg);
  fHistCollector.SetDebug(dbg);
  fEventPlaneFinder.SetDebug(dbg);
}

//____________________________________________________________________
Bool_t
AliForwardMultiplicityTask::SetupForData()
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  DGUARD(fDebug,1,"Initialize sub-algorithms");
  const TAxis* pe = 0;
  const TAxis* pv = 0;

  if (!ReadCorrections(pe,pv)) return false;

  fHistos.Init(*pe);
  fAODFMD.Init(*pe);
  fAODEP.Init(*pe);
  fRingSums.Init(*pe);

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
  
  fFMD1icent=new TH3D("FMD1Ietavcent","FMD1ietavcent;#eta;cent",pe->GetNbins(),pe->GetXmin(),pe->GetXmax(),101,-0.5,100.5,1,0,1);
  fFMD2icent=new TH3D("FMD2Ietavcent","FMD2ietavcent;#eta;cent",pe->GetNbins(),pe->GetXmin(),pe->GetXmax(),101,-0.5,100.5,1,0,1);
  fFMD2ocent=new TH3D("FMD2Oetavcent","FMD2oetavcent;#eta;cent",pe->GetNbins(),pe->GetXmin(),pe->GetXmax(),101,-0.5,100.5,1,0,1);
  fFMD3icent=new TH3D("FMD3Ietavcent","FMD3ietavcent;#eta;cent",pe->GetNbins(),pe->GetXmin(),pe->GetXmax(),101,-0.5,100.5,1,0,1);
  fFMD3ocent=new TH3D("FMD3Oetavcent","FMD3oetavcent;#eta;cent",pe->GetNbins(),pe->GetXmin(),pe->GetXmax(),101,-0.5,100.5,1,0,1);
  fList->Add(fFMD1icent);
  fList->Add(fFMD2icent);
  fList->Add(fFMD2ocent);
  fList->Add(fFMD3icent);
  fList->Add(fFMD3ocent);

	

  this->Print();
  return true;
}

//____________________________________________________________________
void
AliForwardMultiplicityTask::UserCreateOutputObjects()
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create user ouput");
  fList = new TList;
  fList->SetOwner();
  fListVertexBins=new TList();
  fListVertexBins->SetOwner();	
  
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  //if (!ah) AliFatal("No AOD output handler set in analysis manager");
    if (ah)  
   {
	TObject* obj = &fAODFMD;
	ah->AddBranch("AliAODForwardMult", &obj);
  	TObject* epobj = &fAODEP;
 	ah->AddBranch("AliAODForwardEP", &epobj);

   }
    

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
AliForwardMultiplicityTask::UserExec(Option_t*)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  

  DGUARD(fDebug,1,"Process the input event");
  // static Int_t cnt = 0;
  // cnt++;
  // Get the input data 
  AliESDEvent* esd = GetESDEvent();
  if (!esd) return;

  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  fAODFMD.Clear();
  fAODEP.Clear();
  
  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(esd, triggers, lowFlux, 
					       ivz, ip, cent, nClusters);
  
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;

  // Set trigger bits, and mark this event for storage 
  fAODFMD.SetTriggerBits(triggers);
  fAODFMD.SetSNN(fEventInspector.GetEnergy());
  fAODFMD.SetSystem(fEventInspector.GetCollisionSystem());
  fAODFMD.SetCentrality(cent);
  fAODFMD.SetNClusters(nClusters);
  MarkEventForStore();
 
  if (found & AliFMDEventInspector::kNoSPD)      return;
  if (found & AliFMDEventInspector::kNoFMD)      return;
  if (found & AliFMDEventInspector::kNoVertex)   return;
  
  if (triggers & AliAODForwardMult::kPileUp) return;
  
  fAODFMD.SetIpZ(ip.Z());

  if (found & AliFMDEventInspector::kBadVertex) return;

  // We we do not want to use low flux specific code, we disable it here. 
  if (!fEnableLowFlux) lowFlux = false;

  // Get FMD data 
  AliESDFMD* esdFMD = esd->GetFMDData();
  //  // Apply the sharing filter (or hit merging or clustering if you like)
  if (!fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD, ip.Z())) { 
    AliWarning("Sharing filter failed!");
    return;
  }
  
  // Calculate the inclusive charged particle density 
  if (!fDensityCalculator.Calculate(fESDFMD, fHistos, lowFlux, cent, ip)) { 
    // if (!fDensityCalculator.Calculate(*esdFMD, fHistos, ivz, lowFlux)) { 
    AliWarning("Density calculator failed!");
    return;
  }

  if (fEventInspector.GetCollisionSystem() == AliFMDEventInspector::kPbPb) {
    if (!fEventPlaneFinder.FindEventplane(esd, fAODEP, 
					  &(fAODFMD.GetHistogram()), &fHistos))
      AliWarning("Eventplane finder failed!");
  }
  
  // Do the secondary and other corrections. 
  if (!fCorrections.Correct(fHistos, ivz)) { 
    AliWarning("Corrections failed");
    return;
  }

  if (!fHistCollector.Collect(fHistos, fRingSums, 
			      ivz, fAODFMD.GetHistogram(),fList,fAODFMD.GetCentrality(),fListVertexBins)) {
    AliWarning("Histogram collector failed");
    return;
  }

  if (fAODFMD.IsTriggerBits(AliAODForwardMult::kInel))
    fHData->Add(&(fAODFMD.GetHistogram()));

  PostData(1, fList);
}

//____________________________________________________________________
void
AliForwardMultiplicityTask::FinishTaskOutput()
{
  if (!fList) 
    Warning("FinishTaskOutput", "No list defined");
  else {
    if (fDebug) 
      fList->ls();
  }
  AliAnalysisTaskSE::FinishTaskOutput();
}

//____________________________________________________________________
void
AliForwardMultiplicityTask::Terminate(Option_t*)
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  DGUARD(fDebug,1,"Processing the merged results");

  TList* list = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError(Form("No output list defined (%p)", GetOutputData(1)));
    if (GetOutputData(1)) GetOutputData(1)->Print();
    return;
  }

  TList* output = new TList;
  output->SetName(Form("%sResults", GetName()));
  output->SetOwner();

  Double_t nTr = 0, nTrVtx = 0, nAcc = 0;
  MakeSimpledNdeta(list, output, nTr, nTrVtx, nAcc);
  MakeRingdNdeta(list, "ringSums", output, "ringResults");

  fSharingFilter.Terminate(list,output,Int_t(nTr));
  fDensityCalculator.Terminate(list,output,Int_t(nTrVtx));
  fCorrections.Terminate(list,output,Int_t(nTrVtx));

  PostData(2, output);
}

//
// EOF
//
