
//
// Calculate flow in the forward and central regions using the Q cumulants method.
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root or forward_flow.root
//
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TList.h>
#include <THn.h>

#include "AliLog.h"
#include "AliForwardFlowRun2Task.h"
#include "AliForwardQCumulantRun2.h"
#include "AliForwardGenericFramework.h"

#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"

#include "AliForwardUtil.h"

#include "AliVVZERO.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"

#include "AliAnalysisFilter.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

using namespace std;
ClassImp(AliForwardFlowRun2Task)
#if 0
; // For emacs 
#endif

//_____________________________________________________________________
AliForwardFlowRun2Task::AliForwardFlowRun2Task() : AliAnalysisTaskSE(),
  fAOD(0),           // input event
  fOutputList(0),    // output list
  fStdQCList(0), 
  fGFList(0),
  fEventList(0),
  fRandom(0),
  fSettings(),
  fEventCuts(),
  fMultTOFLowCut(),
  fMultTOFHighCut(),
  fMultCentLowCut()
  {
  //
  //  Default constructor
  //
  }

  //_____________________________________________________________________
  AliForwardFlowRun2Task::AliForwardFlowRun2Task(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),           // input event
  fOutputList(0),    // output list
  fStdQCList(0), 
  fGFList(0),
  fEventList(0),
  fRandom(0),
  fSettings(),
  fEventCuts(),
  fMultTOFLowCut(),
  fMultTOFHighCut(),
  fMultCentLowCut()
  {
  // 
  //  Constructor
  //
  //  Parameters:
  //   name: Name of task
  //
    DefineOutput(1, TList::Class());
  }

//_____________________________________________________________________
  void AliForwardFlowRun2Task::UserCreateOutputObjects()
  {
  //
  //  Create output objects
  //




    fOutputList = new TList();          // the final output list

    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested
    
    



//..defining Alex cuts
    /*
    fMultTOFLowCut = new TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFHighCut = new TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultCentLowCut = new TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);

    fMultTOFLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    fMultTOFHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    fMultCentLowCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
  */
  //........

        //..adding QA plots from AliEventCuts class
  //    fEventCuts.fMinVtz = -5.f;
  //fEventCuts.fMaxVtz = 5.f;
  fEventCuts.AddQAplotsToList(fOutputList);






    TRandom r = TRandom();              // random integer to use for creation of samples (used for error bars). 
                                        // Needs to be created here, otherwise it will draw the same random number.

    //int nC2Bins = fSettings.k3pWeight;
    //int nC4Bins = fSettings.kSinphi1phi2phi3p;


    fStdQCList = new TList(); 
    fGFList = new TList();
    fStdQCList->SetName("StdQC"); 
    fGFList->SetName("GF");

    ///TList* fRefList   = new TList();
    //fDiffList  = new TList();
    fEventList = new TList();

    fEventList->Add(new TH1D("Centrality","Centrality",10,0,100));
    fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fEventList->Add(new TH2F("hOutliers","Maximum #sigma from mean N_{ch} pr. bin", 
     20, 0., 100., 500, 0., 5.)); //((fFlags & kMC) ? 15. : 5. // Sigma <M> histogram 
    fEventList->Add(new TH1D("FMDHits","FMDHits",100,0,10));



    fEventList->SetName("EventInfo");
    //fRefList->SetName("Reference");
    //fDiffList->SetName("Differential");

    fStdQCList->Add(new TList());
    fStdQCList->Add(new TList());
    static_cast<TList*>(fStdQCList->At(0))->SetName("Reference");
    static_cast<TList*>(fStdQCList->At(1))->SetName("Differential");  

    fGFList->Add(new TList());
    fGFList->Add(new TList());   
    fGFList->Add(new TList());   
    static_cast<TList*>(fGFList->At(0))->SetName("Reference");
    static_cast<TList*>(fGFList->At(1))->SetName("Differential"); 
    static_cast<TList*>(fGFList->At(2))->SetName("AutoCorrection"); 

    static_cast<TList*>(fGFList->At(2))->Add(new TH1F("fQcorrfactor", "fQcorrfactor", 1, -6.0, 6.0)); //(eta, n)
    static_cast<TList*>(fGFList->At(2))->Add(new TH1F("fpcorrfactor", "fpcorrfactor", fSettings.fNDiffEtaBins, -6.0, 6.0)); //(eta, n)

  //fQcorrfactor = new THnD("fQcorrfactor", "fQcorrfactor", corrdimensions, rcorrbins, cxmin, cxmax);
  //fpcorrfactor = new THnD("fpcorrfactor", "fpcorrfactor", corrdimensions, dcorrbins, cxmin, cxmax);
  //fqcorrfactor = new THnD("fqcorrfactor", "fqcorrfactor", corrdimensions, dcorrbins, cxmin, cxmax);



    fOutputList->Add(fStdQCList);
    fOutputList->Add(fGFList);

    fOutputList->Add(fEventList);



    // do analysis to a maximum of v_5
    Int_t fMaxMoment = 5;
    // create a THn for each harmonic
    for (Int_t n = 2; n <= fMaxMoment; n++) {

      Int_t dimensions = 5;

      Int_t dbins[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, fSettings.fCentBins, fSettings.kSinphi1phi2phi3p+1} ;
      Int_t rbins[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNRefEtaBins, fSettings.fCentBins, fSettings.kSinphi1phi2phi3p+1} ;
      Double_t xmin[5] = {0,fSettings.fZVtxAcceptanceLowEdge, -6.0, 0, 0};
      Double_t xmax[5] = {10,fSettings.fZVtxAcceptanceUpEdge, 6, 100, fSettings.kSinphi1phi2phi3p+1};

      // The THn has dimensions [random samples, vertex position, eta, centrality, kind of variable to store]
      //static_cast<TList*>(fStdQCList->At(0))->Add(new THnD(Form("cumuRef_v%d", n), Form("cumuRef_v%d", n), dimensions, rbins, xmin, xmax));
      //static_cast<TList*>(fStdQCList->At(1))->Add(new THnD(Form("cumuDiff_v%d", n),Form("cumuDiff_%d", n), dimensions, dbins, xmin, xmax));
      static_cast<TList*>(fGFList->At(0))->Add(new THnD(Form("cumuRef_v%d", n), Form("cumuRef_v%d", n), dimensions, rbins, xmin, xmax));
      static_cast<TList*>(fGFList->At(1))->Add(new THnD(Form("cumuDiff_v%d", n),Form("cumuDiff_%d", n), dimensions, dbins, xmin, xmax));
      static_cast<TList*>(fStdQCList->At(0))->Add(new THnD(Form("cumuRef_v%d", n), Form("cumuRef_v%d", n), dimensions, rbins, xmin, xmax));
      static_cast<TList*>(fStdQCList->At(1))->Add(new THnD(Form("cumuDiff_v%d", n),Form("cumuDiff_%d", n), dimensions, dbins, xmin, xmax));

      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(4)->SetName("identifier");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(4)->SetName("identifier");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(4)->SetName("identifier");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(4)->SetName("identifier");


    }
    PostData(1, fOutputList);
  }


//_____________________________________________________________________
void AliForwardFlowRun2Task::UserExec(Option_t */*option*/)
{
  //
  //  Analyses the event with use of the helper class AliForwardQCumulantRun2
  //
  //  Parameters:
  //   option: Not used
  //

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
    
  
//..select trigger
  /*UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isTrigselected = false;
  isTrigselected = fSelectMask&AliVEvent::kINT7;
  if(isTrigselected == false) return;
*/
  //..check if I have AOD
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD){
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }
  
  //..AliEventCuts selection
  if(!fEventCuts.AcceptEvent(fInputEvent)) {
    PostData(1, fOutputList);
    return;
  }  


  AliAODVertex* aodVtx = fAOD->GetPrimaryVertex();
/*
if (fabs( aodVtx->GetZ()) > 5) {
      PostData(1, fOutputList);
    return;
}
*/
  //..get variables for additional event selection cuts (from Alex)
  float v0Centr = 0;
  //float cl1Centr = 0;
  //float cl0Centr = 0;

  //AliMultSelection *MultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  //v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
  

  //cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
  //cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
  //std::cout << v0Centr << std::endl;
/*
  int multTPC = 0;
  int multTPC32 = 0;
  int multTOF = 0;
  int multESD = 0;
  int multTrk = 0;
  const int nTracks = fAOD->GetNumberOfTracks();
  for(int it=0; it<nTracks; it++)
  {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(it);
    if(!track)
    {
      delete track;
      continue;
    }

    if(track->TestFilterBit(128)) multTPC++;

    if(track->TestFilterBit(32))  {
      multTPC32++;
      if ( TMath::Abs(track->GetTOFsignalDz()) <= 10. && track->GetTOFsignal() >= 12000. && track->GetTOFsignal() <= 25000.) multTOF++;
      //if((TMath::Abs(track->Eta())) < fEtaCut && (track->GetTPCNcls() >= fTPCclusters) && (track->Pt() >= fMinPt) && (track->Pt() < fMaxPt)) multTrk++;
    }
  }
  multESD = ((AliAODHeader*)fInputEvent->GetHeader())->GetNumberOfESDTracks();
  float multESDTPCdif = multESD - multTPC*3.38;

  //..do the Alex event selection
  if(multESDTPCdif > 15000) return;
  if(float(multTOF) < fMultTOFLowCut->Eval(float(multTPC32))) return;
  if(float(multTOF) > fMultTOFHighCut->Eval(float(multTPC32))) return;
  //if(((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return;
  if(float(multTrk) < fMultCentLowCut->Eval(v0Centr)) return;
  */
  //..after Alex event selection

  /*fCentralityDis->Fill(v0Centr);

  fCentralityV0MCL1->Fill(v0Centr, cl1Centr);
  fCentralityV0MCL0->Fill(v0Centr, cl0Centr);
  fCentralityCL0CL1->Fill(cl0Centr, cl1Centr);

  fMult128vsCentr->Fill(v0Centr, multTPC);
  fMultvsCentr->Fill(v0Centr, multTrk);
  fMultTPCvsTOF->Fill(multTPC32, multTOF);
  fMultTPCvsESD->Fill(multTPC, multESD);
*/
//..the rest of the analysis



  /*
  UInt_t maskIsSelected =
    dynamic_cast<AliInputEventHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())
    ->IsEventSelected();

  Bool_t isINT7selected = maskIsSelected& AliVEvent::kINT7;

  if (!isINT7selected) return;
*/
  // Get detector objects
  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));

  //AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClusters")); // only exists if created by user from ESDs
  //TH2D spddNdedp = TH2D("spddNdedp","spddNdedp",400,-1.5,1.5,400,0,2*TMath::Pi()); // Histogram to contain the central tracks
TH2D spddNdedp = TH2D("spddNdedp","spddNdedp",400,-1.5,1.5,400,0,2*TMath::Pi()); // Histogram to contain the central tracks


/*
  AliAODTracklets* aodTracklets = fAOD->GetTracklets();

  for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
    spddNdedp.Fill(aodTracklets->GetEta(i),aodTracklets->GetPhi(i), 1);
  }
 */




  Int_t  iTracks(fAOD->GetNumberOfTracks());
  for(Int_t i(0); i < iTracks; i++) {

    // loop  over  all  the  tracks
    AliAODTrack* track = static_cast<AliAODTrack *>(fAOD->GetTrack(i));
    if (track->TestFilterBit(kHybrid)){
      if (track->Pt() >= 0.2 && track->Pt() <= 5){
        spddNdedp.Fill(track->Eta(),track->Phi(), 1);
      }
    }
  }


  //const AliAODTracklets* spdmult = fAOD->GetMultiplicity();

  TH2D& forwarddNdedp = aodfmult->GetHistogram(); // also known as dNdetadphi

  //Double_t cent = (Double_t)aodfmult->GetCentrality();
 // AliAODVertex* aodVtx = fAOD->GetPrimaryVertex();

  if(!fAOD) return;              
  //AliMultSelection* MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  Float_t lPerc = v0Centr; 
  //std::cout << lPerc << std::endl;
  //nonsense 
 /* if ( MultSelection ) {
    lPerc = MultSelection->GetMultiplicityPercentile("V0M");
    //Quality check
    Int_t lEvSelCode = MultSelection->GetEvSelCode();
    if( lEvSelCode > 0 ) lPerc = lEvSelCode; // also if lEvSelCode > 200, the event is probably useless
    //disregard!
  } 
  else
  {
    //If this happens, re-check if AliMultSelectionTask ran before your task!
    AliInfo("Didn't find MultSelection!"); 
  }*/
    
    UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);

      static_cast<TH1D*>(fEventList->FindObject("Centrality"))->Fill(lPerc);
      static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(aodVtx->GetZ());

      //AliForwardQCumulantRun2 calculator = AliForwardQCumulantRun2();
      AliForwardGenericFramework calculator = AliForwardGenericFramework();
      calculator.fSettings = fSettings;

      //double vertexpos = aodVtx->GetZ(); //*spddNdedp

      calculator.CumulantsAccumulate(spddNdedp,fOutputList, lPerc, aodVtx->GetZ(),"central");
      //calculator.fCumuDiff.Reset();
      //calculator.fpvector->Reset();
      if (calculator.useEvent) calculator.CumulantsAccumulate(forwarddNdedp, fOutputList, lPerc, aodVtx->GetZ(),"forward");
      if (calculator.useEvent) calculator.saveEvent(fOutputList, lPerc, aodVtx->GetZ(),  randomInt);
      calculator.reset();

  PostData(1, fOutputList); 

  return;
}




//_____________________________________________________________________
void AliForwardFlowRun2Task::Terminate(Option_t */*option*/)
{
  return;
}


//_____________________________________________________________________
//
//
// EOF