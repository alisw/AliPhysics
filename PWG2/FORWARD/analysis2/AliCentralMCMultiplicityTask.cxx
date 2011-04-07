//====================================================================
// 
// Base class for classes that calculate the multiplicity in the
// central region event-by-event
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODCentralMult 
// 
// Histograms 
//   
// Corrections used 
#include "AliCentralMCMultiplicityTask.h"
#include "AliForwardCorrectionManager.h"
#include "AliForwardUtil.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliFMDEventInspector.h"
#include <AliMCEvent.h>
#include <AliTrackReference.h>
#include <AliStack.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TError.h>
#include <iostream>
#include <iomanip>

//====================================================================
AliCentralMCMultiplicityTask::AliCentralMCMultiplicityTask(const char* name) 
  : AliCentralMultiplicityTask(name),
    fAODMCCentral(kTRUE),
    fMinR(3.5), 
    fMaxR(4.5),
    fMinZ(-14.1), 
    fMaxZ(+14.1),
    fRZ(0), 
    fXYZ(0),
    fNRefs(0)
{
  // 
  // Constructor 
  //   
}
//____________________________________________________________________
AliCentralMCMultiplicityTask::AliCentralMCMultiplicityTask() 
  : AliCentralMultiplicityTask(),
    fAODMCCentral(kTRUE),
    fMinR(3.5), 
    fMaxR(4.5),
    fMinZ(-14.1), 
    fMaxZ(+14.1),
    fRZ(0), 
    fXYZ(0),
    fNRefs(0)
{
  // 
  // Constructor 
  // 
}
//____________________________________________________________________
AliCentralMCMultiplicityTask::AliCentralMCMultiplicityTask(const AliCentralMCMultiplicityTask& o)
  : AliCentralMultiplicityTask(o),
    fAODMCCentral(o.fAODMCCentral),
    fMinR(o.fMinR), 
    fMaxR(o.fMaxR),
    fMinZ(o.fMinZ), 
    fMaxZ(o.fMaxZ),
    fRZ(o.fRZ), 
    fXYZ(o.fXYZ),
    fNRefs(o.fNRefs)
{
  //
  // Copy constructor 
  // 
}
//____________________________________________________________________
AliCentralMCMultiplicityTask&
AliCentralMCMultiplicityTask::operator=(const AliCentralMCMultiplicityTask& o)
{
  // 
  // Assignment operator 
  //
  AliCentralMultiplicityTask::operator=(o);
  fAODMCCentral     = o.fAODMCCentral;
  fMinR             = o.fMinR;
  fMaxR             = o.fMaxR;
  fMinZ             = o.fMinZ;
  fMaxZ             = o.fMaxZ;
  fRZ               = o.fRZ;
  fXYZ              = o.fXYZ;
  fNRefs            = o.fNRefs;
  return *this;
}
//____________________________________________________________________
void AliCentralMCMultiplicityTask::UserCreateOutputObjects() 
{
  // 
  // Create output objects 
  // 
  //
  AliCentralMultiplicityTask::UserCreateOutputObjects();

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) AliFatal("No AOD output handler set in analysis manager");
  
  
  TObject* obj = &fAODMCCentral;
  ah->AddBranch("AliAODCentralMult", &obj);

  fRZ = new TH2D("rz", "(r,z) of used track references", 
		 50, 0, 5, 30, -15, 15);
  fRZ->SetXTitle("z [cm]");
  fRZ->SetYTitle("r [cm]");
  fRZ->SetDirectory(0);
  fList->Add(fRZ);

  fXYZ = new TH3D("xyz", "(x,y,z) of used track references", 
		  100, -5., 5., 100, -5., 5., 30, -15., 15.);
  fXYZ->SetXTitle("x [cm]");
  fXYZ->SetYTitle("y [cm]");
  fXYZ->SetZTitle("z [cm]");
  fXYZ->SetDirectory(0);
  fList->Add(fXYZ);

  fNRefs = new TH1D("nrefs", "Number of references used per track", 
		    11, -.5, 10.5); 
  fNRefs->SetXTitle("# of references per track");
  fNRefs->SetDirectory(0);
  fList->Add(fNRefs);
}
//____________________________________________________________________
void AliCentralMCMultiplicityTask::UserExec(Option_t* option) 
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  fAODMCCentral.Clear("");
  
  // Call base class 
  AliCentralMultiplicityTask::UserExec(option);

  // check if we need this event 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah)  
    AliFatal("No AOD output handler set in analysis manager");

  // if base class did not want this event, then neither to we 
  if (!ah->GetFillAOD() || fIvz <= 0) return;
  
  AliMCEvent*  mcEvent = MCEvent();
  TH2D&        hist    = fAODMCCentral.GetHistogram();

  ProcessMC(hist, mcEvent);
  CorrectData(hist, fIvz);

}
//____________________________________________________________________
void AliCentralMCMultiplicityTask::ProcessMC(TH2D& hist, 
					     const AliMCEvent* event) const
{
  Double_t vz = GetManager().GetSecMap()->GetVertexAxis().GetBinCenter(fIvz);

  AliStack* stack = const_cast<AliMCEvent*>(event)->Stack();
  Int_t nTracks   = stack->GetNtrack();//event.GetNumberOfTracks();
  // Int_t nPrim     = stack->GetNprimary();//event.GetNumberOfPrimary();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { 
    AliMCParticle* particle = 
      static_cast<AliMCParticle*>(event->GetTrack(iTr));
    
    // Check the returned particle 
    if (!particle) continue;
    
    // Check if this charged and a primary 
    Bool_t isCharged = particle->Charge() != 0;
    if (!isCharged) continue;

    // Pseudo rapidity and azimuthal angle 
    // Double_t eta = particle->Eta();
    // Double_t phi = particle->Phi();

    Int_t    nTrRef  = particle->GetNumberOfTrackReferences();
    Int_t    nRef    = 0;
    for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) { 
      AliTrackReference* ref = particle->GetTrackReference(iTrRef);
      
      // Check existence 
      if (!ref) continue;

      // Check that we hit an FMD element 
      if (ref->DetectorId() != AliTrackReference::kITS) 
	continue;
      
      // Get radius and z where the track reference was made 
      Double_t r = ref->R();
      Double_t x = ref->X();
      Double_t y = ref->Y();
      Double_t z = ref->Z();
      if (r > fMaxR || r < fMinR) continue;
      if (z > fMaxZ || z < fMinZ) continue;

      fRZ->Fill(z,r);
      fXYZ->Fill(x, y, z);
      // Only fill first reference 
      if (nRef == 0) { 
	Double_t zr = z-vz;
	Double_t th = TMath::ATan2(r,zr);
	if (th < 0) th += 2*TMath::Pi();
	Double_t et = -TMath::Log(TMath::Tan(th/2));
	Double_t ph = TMath::ATan2(y,x);
	if (ph < 0) ph += 2*TMath::Pi();
	hist.Fill(et,ph);
      }
      nRef++;
    }
    fNRefs->Fill(nRef);
  }
}

//____________________________________________________________________
void AliCentralMCMultiplicityTask::Terminate(Option_t* option) 
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  AliCentralMultiplicityTask::Terminate(option);
}
//____________________________________________________________________
void
AliCentralMCMultiplicityTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  AliCentralMultiplicityTask::Print(option);
}
//
// EOF
//
