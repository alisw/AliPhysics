#include "TChain.h"

#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliTwoParticlePIDCorrKine.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliVEvent.h"
#include "AliVParticle.h"


ClassImp(AliTwoParticlePIDCorrKine)

//____________________________________________| Constructor
AliTwoParticlePIDCorrKine::AliTwoParticlePIDCorrKine():
fEvent(0x0),
fMcHandler(0x0),
  fHistEventsProcessed(0x0),
  fOutputList(0),
  fHistZvtx(0),
  fHistPt(0),
  fHistImpact(0x0),
fZvtxLim(10),
fCentralityFrom("Impact"),
fCentralityEstimator("V0M")

{
  //Default constructor 
}

//____________________________________________| Specific Constructor
AliTwoParticlePIDCorrKine::AliTwoParticlePIDCorrKine(const Char_t* name) :
  AliAnalysisTaskSE(name),
fEvent(0x0),
fMcHandler(0x0),
  fHistEventsProcessed(0x0),
  fOutputList(0),
  fHistZvtx(0),
  fHistPt(0),
  fHistImpact(0x0),
fZvtxLim(10),
fCentralityFrom("Impact"),
fCentralityEstimator("V0M")
{
  // Constructor. Initialization of Inputs and Outputs
  Info("AliTwoParticlePIDCorrKine","Calling Constructor");
  // Output slot #1 writes into a TList container (nevents histogram)
    DefineInput(0, TChain::Class());
    DefineOutput(1,TList::Class()); // Basic output slot (more needed)
}

//____________________________________________| Destructor
AliTwoParticlePIDCorrKine::~AliTwoParticlePIDCorrKine()
{
  // Destructor
  Info("~AliTwoParticlePIDCorrKine","Calling Destructor");
  if (fHistEventsProcessed) delete fHistEventsProcessed;
  if (fOutputList) delete fOutputList;
  if (fHistZvtx) delete fHistZvtx;
  if (fHistPt) delete fHistPt;

}

//___________________________________________________________________________
void AliTwoParticlePIDCorrKine::UserCreateOutputObjects()
{

  //	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  //      AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    
    Info("CreateOutputObjects","CreateOutputObjects of task %s", GetName());
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    //    fOutputList->SetName("KineTest");


    fHistEventsProcessed = new TH1F("fHistNEvents","fHistEventsProcessed",3,-0.5,2.5) ;
    fHistEventsProcessed->GetXaxis()->SetBinLabel(1,"All events");
    fHistEventsProcessed->GetXaxis()->SetBinLabel(2,"Event within |Ztx| 10cm");
    fHistEventsProcessed->GetXaxis()->SetBinLabel(3,"Good Reconstructed events");
    
    
    fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
    fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPt->SetMarkerStyle(kFullCircle);
    
    fHistZvtx = new TH1F("fHistZvtx", "Zvtx distribution", 40, -20, 20);
    fHistZvtx->GetXaxis()->SetTitle("ZVtx (cm)");
    fHistZvtx->GetYaxis()->SetTitle("Nch");
    fHistZvtx->SetMarkerStyle(kFullCircle);

    fHistImpact = new TH1F ("fHistImpact","ImpactParameter_Dist", 100,0.0,20.0);
    
    fOutputList->Add(fHistEventsProcessed);
    fOutputList->Add(fHistZvtx);
    fOutputList->Add(fHistPt);
    fOutputList->Add(fHistImpact);

    PostData(1, fOutputList);

    
    return;
}



//____________________________________________| User Exec
void AliTwoParticlePIDCorrKine::UserExec(Option_t *)
{
  
 // fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); //get handler
//if(!fMCHandler) return;

 fEvent = dynamic_cast<AliVEvent*>(MCEvent());
  if(!fEvent)
    {
      //AliError("event not available");
      return;
    }


    fHistEventsProcessed->Fill(0.0); // # of Event after passing MC cuts
    
    


       
    
    const AliVVertex *vtxMC = fEvent->GetPrimaryVertex();
    Float_t zVtx = vtxMC->GetZ();
    if(TMath::Abs(zVtx) < fZvtxLim)fHistZvtx->Fill(zVtx);
    else return;

fHistEventsProcessed->Fill(1.0);

  Double_t gImpactParameter = 0.;
  //Double_t  gMultiplicity = 0.;
  //Double_t gReactionPlane=0.0;
  //Double_t gCentrality=0.0;

 AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(fEvent);
  if(gMCEvent){
AliCollisionGeometry* headerH = dynamic_cast<AliCollisionGeometry*>(gMCEvent->GenEventHeader());
    if(headerH){
      gImpactParameter = headerH->ImpactParameter();
      //gMultiplicity =    GenMultiplicity(event); //calculate the multiplicity depending on the choice of centrality estimator                                           
      // gReactionPlane=     headerH->ReactionPlaneAngle();
     // gCentrality = (CentralityFrom =="Multiplicity") ? gImpactParameter :gMultiplicity;
    }
}
  
	fHistImpact->Fill(gImpactParameter);
  
     // # of Event after passing MC cuts

    //Printf("MC particles: %d", fEvent->GetNumberOfTracks());
    
    for (Int_t iTracks = 0; iTracks < fEvent->GetNumberOfTracks(); iTracks++) {
        AliVParticle* track = fEvent->GetTrack(iTracks);
       if (!track) {
      //AliError(Form("Could not receive particle %d", iTracks));
      continue;
    }

        fHistPt->Fill(track->Pt());
    }
    
    
    fHistEventsProcessed->Fill(2.0); // # of Event after passing MC cuts

  PostData(1, fOutputList);
  return;
}


//___________________________________________________________________________
void AliTwoParticlePIDCorrKine::Terminate(Option_t*)
{
  
  Info("Terminate","Start and end of Method");
  AliAnalysisTaskSE::Terminate();
    
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
        Printf("ERROR: Output list not available");
        return;
  }

}
