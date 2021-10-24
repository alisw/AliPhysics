/*
 * AliAnalysisTaskGeorgiosNTuple.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */
#include "AliAnalysisTaskGeorgiosNTuple.h"
#include "AliNanoAODTrack.h"
#include "TRandom3.h"


ClassImp(AliAnalysisTaskGeorgiosNTuple)
AliAnalysisTaskGeorgiosNTuple::AliAnalysisTaskGeorgiosNTuple()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fisMC(false),  
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fv0(nullptr),
      fLambda(nullptr),
      fLambdaList(nullptr),
      fAntiLambda(nullptr),
      fAntiLambdaList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fXiBGR(nullptr),
      fXiBGRList(nullptr),
      fAntiXiBGR(nullptr),
      fAntiXiBGRList(nullptr),
#ifdef MONTECARLO
      fLambdaMCList(nullptr),
      fAntiLambdaMCList(nullptr),
      fXiMCList(nullptr),
      fAntiXiMCList(nullptr),
      fXiBGRMCList(nullptr),
      fAntiXiBGRMCList(nullptr),
#endif
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fGeorgiosTree(0),
      fGTI(nullptr) {
}
//former version had in the constructor bool XiiTreeFlag
AliAnalysisTaskGeorgiosNTuple::AliAnalysisTaskGeorgiosNTuple(const char* name, bool isMC)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fisMC(isMC),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fv0(nullptr),
      fLambda(nullptr),
      fLambdaList(nullptr),
      fAntiLambda(nullptr),
      fAntiLambdaList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fXiBGR(nullptr),
      fXiBGRList(nullptr),
      fAntiXiBGR(nullptr),
      fAntiXiBGRList(nullptr),
#ifdef MONTECARLO
      fLambdaMCList(nullptr),
      fAntiLambdaMCList(nullptr),
      fXiMCList(nullptr),
      fAntiXiMCList(nullptr),
      fXiBGRMCList(nullptr),
      fAntiXiBGRMCList(nullptr),
#endif
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fGeorgiosTree(0),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Lambda Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiLambda Cuts
  DefineOutput(4, TList::Class());  //Output for the XiBGR Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiXiBGR Cuts
  DefineOutput(6, TList::Class());  //Output for the Xi Cuts
  DefineOutput(7, TList::Class());  //Output for the AntiXi Cuts
  DefineOutput(8, TList::Class());  //Output for the Results
  DefineOutput(9, TList::Class());  //Output for the Results QA
  DefineOutput(10, TTree::Class());  // XiTree (former OmegaTree)

#ifdef MONTECARLO
   DefineOutput(11, TList::Class());  //Output for the v0Cuts MC
   DefineOutput(12, TList::Class());  //Output for the Antiv0Cuts MC
   DefineOutput(13, TList::Class());  //Output for the XiBGRCuts MC
   DefineOutput(14, TList::Class());  //Output for the AntiXiBGRCuts MC
   DefineOutput(15, TList::Class());  //Output for the XiCuts MC
   DefineOutput(16, TList::Class());  //Output for the AntiXiCuts MC
#endif
}

AliAnalysisTaskGeorgiosNTuple::~AliAnalysisTaskGeorgiosNTuple() {
  if (fEvent) {
    delete fEvent;
  }
  if (fEventCuts) {
    delete fEventCuts;
  }
//Lambda
  if (fv0) {
    delete fv0;
  }
  if (fLambda) {
    delete fLambda;
  }
  if (fAntiLambda) {
    delete fAntiLambda;
  }
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
  if (fGeorgiosTree) {
    delete fGeorgiosTree;
  }
//Xi
  if (fCascade) {
    delete fCascade;
  }
  if (fXi) {
    delete fXi;
  }
  if (fAntiXi) {
    delete fAntiXi;
  }
  if (fXiBGR) {
    delete fXiBGR;
  }
  if (fAntiXiBGR) {
    delete fAntiXiBGR;
  }

}

void AliAnalysisTaskGeorgiosNTuple::UserCreateOutputObjects() {           
  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }
//Lambda
  if (!fLambda) {
    AliError("No Lambda cuts \n");
  } else {
    fLambda->Init();
  }
  if (!fAntiLambda) {
    AliError("No AntiLambda cuts \n");
  } else {
    fAntiLambda->Init();
  }

//Xi
  if (!fXiBGR) {
    AliError("No XiBGR cuts \n");
  } else {
    fXiBGR->Init();
  }
  if (!fAntiXiBGR) {
    AliError("No AntiXiBGR cuts \n");
  } else {
    fAntiXiBGR->Init();
  }
  if (!fXi) {
    AliError("No Xi cuts \n");
  } else {
    fXi->Init();
  }
  if (!fAntiXi) {
    AliError("No AntiXi cuts \n");
  } else {
    fAntiXi->Init();
  }


  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 2,
                                                fConfig->GetMinimalBookingME());
  }
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), false);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());


//Lambda
  fv0 = new AliFemtoDreamv0();
  fv0->SetUseMCInfo(fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());
  fv0->SetPDGCode(3122);
  fv0->SetPDGDaughterPos(2212);
  fv0->GetPosDaughter()->SetUseMCInfo(fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());
  fv0->SetPDGDaughterNeg(-211);
  fv0->GetNegDaughter()->SetUseMCInfo(fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());

//Xi
  fCascade = new AliFemtoDreamCascade();
  fCascade->SetUseMCInfo(fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGCode(3312);
  fCascade->SetPDGDaugPos(2212);            //Proton
  fCascade->GetPosDaug()->SetUseMCInfo(fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGDaugNeg(-211);             //pi^-
  fCascade->GetNegDaug()->SetUseMCInfo(fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGDaugBach(-211);            //EDIT pi^-
  fCascade->GetBach()->SetUseMCInfo(fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->Setv0PDGCode(3122);             //Lambda


  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }


//Lambda
if (!fLambda->GetMinimalBooking()) {
  fLambdaList = fLambda->GetQAHists();
} else {
  fLambdaList = new TList();
  fLambdaList->SetName("LambdaCuts");
  fLambdaList->SetOwner();
}
if (!fAntiLambda->GetMinimalBooking()) {
  fAntiLambdaList = fAntiLambda->GetQAHists();
} else {
  fAntiLambdaList = new TList();
  fAntiLambdaList->SetName("AntiLambdaCuts");
  fAntiLambdaList->SetOwner();
}

//Xi
  if (!fXiBGR->GetMinimalBooking()) {
    fXiBGRList = fXiBGR->GetQAHists();
  } else {
    fXiBGRList = new TList();
    fXiBGRList->SetName("XiBGRCuts");
    fXiBGRList->SetOwner();
  }
  if (!fAntiXiBGR->GetMinimalBooking()) {
    fAntiXiBGRList = fAntiXiBGR->GetQAHists();
  } else {
    fAntiXiBGRList = new TList();
    fAntiXiBGRList->SetName("AntiXiBGRCuts");
    fAntiXiBGRList->SetOwner();
  }

  if (!fXi->GetMinimalBooking()) {
    fXiList = fXi->GetQAHists();
  } else {
    fXiList = new TList();
    fXiList->SetName("XiCuts");
    fXiList->SetOwner();
  }
  if (!fAntiXi->GetMinimalBooking()) {
    fAntiXiList = fAntiXi->GetQAHists();
  } else {
    fAntiXiList = new TList();
    fAntiXiList->SetName("AntiXiCuts");
    fAntiXiList->SetOwner();
  }


  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");
  if (!fConfig->GetMinimalBookingME()) {
    fResults = fPartColl->GetHistList();
    fResultsQA->Add(fPartColl->GetQAList());
    fResultsQA->Add(fPairCleaner->GetHistList());
  } else {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }

  //book tree:
  ////////////
  fGeorgiosTree = new TTree("oTTree","a simple TTree");
  fGeorgiosTree->SetDirectory(0); // This is to force a memory-resident Tree, and avoid errors. // ????? is this necessary? does it create memory problems?
  fGeorgiosTree->Branch("RunNumber",&fTRunNumber,"fTRunNumber/I");
  fGeorgiosTree->Branch("Vz",&fTVz,"fTVz/F");
  fGeorgiosTree->Branch("Mult",&fTMult,"fTMult/I");


//Lambda
fGeorgiosTree->Branch("nv0",&fTnv0,"fTnv0/I");
fGeorgiosTree->Branch("v0Charge",&fTv0Charge,"fTv0Charge[fTnv0]/S");
fGeorgiosTree->Branch("v0DCA",&fTv0DCA,"fTv0DCA[fTnv0]/F");
fGeorgiosTree->Branch("v0DaughtersDCA",&fTv0DaughtersDCA,"fTv0DaughtersDCA[fTnv0]/F");
//fGeorgiosTree->Branch("v0LambdaMass",&fTv0LambdaMass,"fTv0LambdaMass[fTnv0]/F");
fGeorgiosTree->Branch("v0LambdaVr",&fTv0LambdaVr,"fTv0LambdaVr[fTnv0]/F");
fGeorgiosTree->Branch("v0LambdaPA",&fTv0LambdaPA,"fTv0LambdaPA[fTnv0]/F");
fGeorgiosTree->Branch("Trackv0Px",&fTTrackv0Px,"fTTrackv0Px[fTnv0][2]/F");
fGeorgiosTree->Branch("Trackv0Py",&fTTrackv0Py,"fTTrackv0Py[fTnv0][2]/F");
fGeorgiosTree->Branch("Trackv0Pz",&fTTrackv0Pz,"fTTrackv0Pz[fTnv0][2]/F");
//fGeorgiosTree->Branch("Trackv0Eta",&fTTrackv0Eta,"fTTrackv0Eta[fTnv0][2]/F");
fGeorgiosTree->Branch("Trackv0Charge",&fTTrackv0Charge,"fTTrackv0Charge[fTnv0][2]/S");
fGeorgiosTree->Branch("Trackv0DCA",&fTTrackv0DCA,"fTTrackv0DCA[fTnv0][2]/F");
fGeorgiosTree->Branch("Trackv0TPCsigma",&fTTrackv0TPCsigma,"fTTrackv0TPCsigma[fTnv0][2]/F");
//fGeorgiosTree->Branch("Trackv0TOFsigma",&fTTrackv0TOFsigma,"fTTrackv0TOFsigma[fTnv0][2]/F");
fGeorgiosTree->Branch("Trackv0Ncl",&fTTrackv0Ncl,"fTTrackv0Ncl[fTnv0][2]/I");
/*
 * reduce disk space for background
fGeorgiosTree->Branch("Trackv0CrR",&fTTrackv0CrR,"fTTrackv0CrR[fTnv0][2]/F");
fGeorgiosTree->Branch("Trackv0CrF",&fTTrackv0CrF,"fTTrackv0CrF[fTnv0][2]/F");
fGeorgiosTree->Branch("Trackv0FilterBit",&fTTrackv0FilterBit,"fTTrackv0FilterBit[fTnv0][2]/O");
fGeorgiosTree->Branch("Trackv0Phi",&fTTrackv0Phi,"fTTrackv0Phi[fTnv0][2]/F");
fGeorgiosTree->Branch("Trackv0ITStime",&fTTrackv0ITStime,"fTTrackv0ITStime[fTnv0][2]/O");
fGeorgiosTree->Branch("Trackv0TOFtime",&fTTrackv0TOFtime,"fTTrackv0TOFtime[fTnv0][2]/O");
*/
fGeorgiosTree->Branch("Trackv0ID",&fTTrackv0ID,"fTTrackv0ID[fTnv0][2]/I");
	
#ifdef MONTECARLO
//Lambda MC
  fGeorgiosTree->Branch("v0LambdaPxMC",&fTv0LambdaPxMC,"fTv0LambdaPxMC[fTnv0]/F");
  fGeorgiosTree->Branch("v0LambdaPyMC",&fTv0LambdaPyMC,"fTv0LambdaPyMC[fTnv0]/F");
  fGeorgiosTree->Branch("v0LambdaPzMC",&fTv0LambdaPzMC,"fTv0LambdaPzMC[fTnv0]/F");
  fGeorgiosTree->Branch("v0LambdaPDG",&fTv0LambdaPDG,"fTv0LambdaPDG[fTnv0]/I");
  fGeorgiosTree->Branch("v0LambdaMotherPDG",&fTv0LambdaMotherPDG,"fTv0LambdaMotherPDG[fTnv0]/I");
  fGeorgiosTree->Branch("v0LambdaMotherWeakPDG",&fTv0LambdaMotherWeakPDG,"fTv0LambdaMotherWeakPDG[fTnv0]/I");
  fGeorgiosTree->Branch("v0LambdaOrigin",&fTv0LambdaOrigin,"fTv0LambdaOrigin[fTnv0]/I");
  fGeorgiosTree->Branch("Trackv0PxMC",&fTTrackv0PxMC,"fTTrackv0PxMC[fTnv0][2]/F");
  fGeorgiosTree->Branch("Trackv0PyMC",&fTTrackv0PyMC,"fTTrackv0PyMC[fTnv0][2]/F");
  fGeorgiosTree->Branch("Trackv0PzMC",&fTTrackv0PzMC,"fTTrackv0PzMC[fTnv0][2]/F");
  fGeorgiosTree->Branch("Trackv0PDG",&fTTrackv0PDG,"fTTrackv0PDG[fTnv0][2]/I");
  fGeorgiosTree->Branch("Trackv0MotherPDG",&fTTrackv0MotherPDG,"fTTrackv0MotherPDG[fTnv0][2]/I");
  fGeorgiosTree->Branch("Trackv0MotherWeakPDG",&fTTrackv0MotherWeakPDG,"fTTrackv0MotherWeakPDG[fTnv0][2]/I");
  fGeorgiosTree->Branch("Trackv0Origin",&fTTrackv0Origin,"fTTrackv0Origin[fTnv0][2]/I");
#endif

 //Xis:
 fGeorgiosTree->Branch("nCascade",&fTnCascade,"fTnCascade/I"); //
 fGeorgiosTree->Branch("CascadeCharge",&fTCascadeCharge,"fTCascadeCharge[fTnCascade]/S");
 fGeorgiosTree->Branch("CascadeDCA",&fTCascadeDCA,"fTCascadeDCA[fTnCascade]/F");
 fGeorgiosTree->Branch("CascadeDaughtersDCA",&fTCascadeDaughtersDCA,"fTCascadeDaughtersDCA[fTnCascade]/F");
// fGeorgiosTree->Branch("CascadeXiMass",&fTCascadeXiMass,"fTCascadeXiMass[fTnCascade]/F");
// fGeorgiosTree->Branch("CascadeOmegaMass",&fTCascadeOmegaMass,"fTCascadeOmegaMass[fTnCascade]/F");
 fGeorgiosTree->Branch("CascadeVr",&fTCascadeVr,"fTCascadeVr[fTnCascade]/F");
 fGeorgiosTree->Branch("CascadePA",&fTCascadePA,"fTCascadePA[fTnCascade]/F");
 fGeorgiosTree->Branch("LambdaDCA",&fTLambdaDCA,"fTLambdaDCA[fTnCascade]/F");
 fGeorgiosTree->Branch("LambdaDaughtersDCA",&fTLambdaDaughtersDCA,"fTLambdaDaughtersDCA[fTnCascade]/F");
// fGeorgiosTree->Branch("LambdaMass",&fTLambdaMass,"fTLambdaMass[fTnCascade]/F");
 fGeorgiosTree->Branch("LambdaVr",&fTLambdaVr,"fTLambdaVr[fTnCascade]/F");
 fGeorgiosTree->Branch("LambdaPA",&fTLambdaPA,"fTLambdaPA[fTnCascade]/F");
 fGeorgiosTree->Branch("TrackPx",&fTTrackPx,"fTTrackPx[fTnCascade][3]/F");
 fGeorgiosTree->Branch("TrackPy",&fTTrackPy,"fTTrackPy[fTnCascade][3]/F");
 fGeorgiosTree->Branch("TrackPz",&fTTrackPz,"fTTrackPz[fTnCascade][3]/F");
// fGeorgiosTree->Branch("TrackEta",&fTTrackEta,"fTTrackEta[fTnCascade][3]/F");
 fGeorgiosTree->Branch("TrackCharge",&fTTrackCharge,"fTTrackCharge[fTnCascade][3]/S");
 fGeorgiosTree->Branch("TrackDCA",&fTTrackDCA,"fTTrackDCA[fTnCascade][3]/F");
 fGeorgiosTree->Branch("TrackTPCsigma",&fTTrackTPCsigma,"fTTrackTPCsigma[fTnCascade][3]/F");
// fGeorgiosTree->Branch("TrackTOFsigma",&fTTrackTOFsigma,"fTTrackTOFsigma[fTnCascade][3]/F");
 fGeorgiosTree->Branch("TrackNcl",&fTTrackNcl,"fTTrackNcl[fTnCascade][3]/I");
/*
 * reducing disk space for background
 fGeorgiosTree->Branch("TrackCrR",&fTTrackCrR,"fTTrackCrR[fTnCascade][3]/F");
 fGeorgiosTree->Branch("TrackCrF",&fTTrackCrF,"fTTrackCrF[fTnCascade][3]/F");
 fGeorgiosTree->Branch("TrackFilterBit",&fTTrackFilterBit,"fTTrackFilterBit[fTnCascade][3]/O");
*/
// fGeorgiosTree->Branch("TrackPhi",&fTTrackPhi,"fTTrackPhi[fTnCascade][3]/F");
 fGeorgiosTree->Branch("TrackITStime",&fTTrackITStime,"fTTrackITStime[fTnCascade][3]/O");
 fGeorgiosTree->Branch("TrackTOFtime",&fTTrackTOFtime,"fTTrackTOFtime[fTnCascade][3]/O");
 fGeorgiosTree->Branch("TrackID",&fTTrackID,"fTTrackID[fTnCascade][3]/I");

#ifdef MONTECARLO
 //XisMC:
  fGeorgiosTree->Branch("CascadePxMC",&fTCascadePxMC,"fTCascadePxMC[fTnCascade]/F");
  fGeorgiosTree->Branch("CascadePyMC",&fTCascadePyMC,"fTCascadePyMC[fTnCascade]/F");
  fGeorgiosTree->Branch("CascadePzMC",&fTCascadePzMC,"fTCascadePzMC[fTnCascade]/F");
  fGeorgiosTree->Branch("CascadePDG",&fTCascadePDG,"fTCascadePDG[fTnCascade]/I");
  fGeorgiosTree->Branch("CascadeMotherPDG",&fTCascadeMotherPDG,"fTCascadeMotherPDG[fTnCascade]/I");
  //fGeorgiosTree->Branch("CascadeMotherWeakPDG",&fTCascadeMotherWeakPDG,"fTCascadeMotherWeakPDG[fTnCascade]/I");
  fGeorgiosTree->Branch("CascadeOrigin",&fTCascadeOrigin,"fTCascadeOrigin[fTnCascade]/I");
  fGeorgiosTree->Branch("TrackPxMC",&fTTrackPxMC,"fTTrackPxMC[fTnCascade][3]/F");
  fGeorgiosTree->Branch("TrackPyMC",&fTTrackPyMC,"fTTrackPyMC[fTnCascade][3]/F");
  fGeorgiosTree->Branch("TrackPzMC",&fTTrackPzMC,"fTTrackPzMC[fTnCascade][3]/F");
  fGeorgiosTree->Branch("TrackPDG",&fTTrackPDG,"fTTrackPDG[fTnCascade][3]/I");
  fGeorgiosTree->Branch("TrackMotherPDG",&fTTrackMotherPDG,"fTTrackMotherPDG[fTnCascade][3]/I");
  fGeorgiosTree->Branch("TrackMotherWeakPDG",&fTTrackMotherWeakPDG,"fTTrackMotherWeakPDG[fTnCascade][3]/I");
  fGeorgiosTree->Branch("TrackOrigin",&fTTrackOrigin,"fTTrackOrigin[fTnCascade][3]/I");
#endif

 PostData(1, fEvtList);
 PostData(2, fLambdaList);
 PostData(3, fAntiLambdaList);
 PostData(4, fXiBGRList);
 PostData(5, fAntiXiBGRList);      //
 PostData(6, fXiList);
 PostData(7, fAntiXiList);
 PostData(8, fResults);
 PostData(9, fResultsQA);
 PostData(10, fGeorgiosTree);


#ifdef MONTECARLO 
//LambdaMC
if (fLambda->GetIsMonteCarlo()) {
 if (!fLambda->GetMinimalBooking()) {
   fLambdaMCList = fLambda->GetMCQAHists();
 } else {
   fLambdaMCList = new TList();
   fLambdaMCList->SetName("MCLambdaCuts");
   fLambdaMCList->SetOwner();
 }
 PostData(11, fLambdaMCList);
}
if (fAntiLambda->GetIsMonteCarlo()) {
 if (!fAntiLambda->GetMinimalBooking()) {
   fAntiLambdaMCList = fAntiLambda->GetMCQAHists();
 } else {
   fAntiLambdaMCList = new TList();
   fAntiLambdaMCList->SetName("MCAntiLambdaCuts");
   fAntiLambdaMCList->SetOwner();
 }
 PostData(12, fAntiLambdaMCList);

}

//Xi
if (fXiBGR->GetIsMonteCarlo()) {
 if (!fXiBGR->GetMinimalBooking()) {
   fXiBGRMCList = fXiBGR->GetMCQAHists();
 } else {
   fXiBGRMCList = new TList();
   fXiBGRMCList->SetName("MCXiBGRCuts");
   fXiBGRMCList->SetOwner();
 }
 PostData(13, fXiBGRMCList);
}

if (fAntiXiBGR->GetIsMonteCarlo()) {
 if (!fAntiXiBGR->GetMinimalBooking()) {
   fAntiXiBGRMCList = fAntiXiBGR->GetMCQAHists();
 } else {
   fAntiXiBGRMCList = new TList();
   fAntiXiBGRMCList->SetName("MCAntiXiBGRCuts");
   fAntiXiBGRMCList->SetOwner();
 }
 PostData(14, fAntiXiBGRMCList);
}

if (fXi->GetIsMonteCarlo()) {
 if (!fXi->GetMinimalBooking()) {
   fXiMCList = fXi->GetMCQAHists();
 } else {
   fXiMCList = new TList();
   fXiMCList->SetName("MCXiCuts");
   fXiMCList->SetOwner();
 }
 PostData(15, fXiMCList);
}

if (fAntiXi->GetIsMonteCarlo()) {
 if (!fAntiXi->GetMinimalBooking()) {
   fAntiXiMCList = fAntiXi->GetMCQAHists();
 } else {
   fAntiXiMCList = new TList();
   fAntiXiMCList->SetName("MCAntiXiCuts");
   fAntiXiMCList->SetOwner();
 }
 PostData(16, fAntiXiMCList);
}
#endif

}


void AliAnalysisTaskGeorgiosNTuple::UserExec(Option_t *option) {
//  AliVEvent *fInputEvent = InputEvent();
  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }

  fTRunNumber = 0.; //For NanoAOD filtering trains <100, no info

  Double_t PrimVtx[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(PrimVtx);
  fTVz = PrimVtx[2];
  fTMult = fEvent->GetMultiplicity();

  for(int ii=0; ii<MAXv0; ii++){
   fTv0Charge[ii]=-10;
   fTv0DCA[ii]=-100000.;
   fTv0DaughtersDCA[ii]=-100000.;
//   fTv0LambdaMass[ii]=-100000.;
   fTv0LambdaVr[ii]=-100000.;
   fTv0LambdaPA[ii]=-100000.;
 
   fTv0Vx[ii]=-100000.;
   fTv0Vy[ii]=-100000.;
   fTv0Vz[ii]=-100000.;
  for (int jj=0; jj<2; jj++){
    fTTrackv0Px[ii][jj]=-100000.;
    fTTrackv0Py[ii][jj]=-100000.;
    fTTrackv0Pz[ii][jj]=-100000.;
    fTTrackv0TPCmom[ii][jj]=-100000.;
//    fTTrackv0Eta[ii][jj]=-100000.;
    fTTrackv0Charge[ii][jj]=-10;
    fTTrackv0DCA[ii][jj]=-100000.;
    fTTrackv0TPCsigma[ii][jj]=-100000.;
    //fTTrackv0TOFsigma[ii][jj]=-100000.;
    fTTrackv0Ncl[ii][jj]=-100000;
/*
 * reduce disk space for background
    fTTrackv0CrR[ii][jj]=-100000.;
    fTTrackv0CrF[ii][jj]=-100000.;
    fTTrackv0Shared[ii][jj]=-100000;
    fTTrackv0TPCchi2[ii][jj]=-100000.;
    fTTrackv0TPConly[ii][jj]=kFALSE;
    fTTrackv0ITScomplementary[ii][jj]=kFALSE;
    fTTrackv0ITSpure[ii][jj]=kFALSE;
    fTTrackv0GLOBAL[ii][jj]=kFALSE;
    fTTrackv0FilterBit[ii][jj]=0;
    fTTrackv0Phi[ii][jj]=-100000.;
    fTTrackv0ITStime[ii][jj]=kFALSE;
    fTTrackv0TOFtime[ii][jj]=kFALSE;
*/
    fTTrackv0ID[ii][jj]=-100000;

#ifdef MONTECARLO
     fTTrackv0PxMC[ii][jj]=-100000.;
     fTTrackv0PyMC[ii][jj]=-100000.;
     fTTrackv0PzMC[ii][jj]=-100000.;
     fTTrackv0PDG[ii][jj]=0;
     fTTrackv0MotherPDG[ii][jj]=0;
     fTTrackv0MotherWeakPDG[ii][jj]=0;
     fTTrackv0Origin[ii][jj]=10;
#endif 
  }
  
#ifdef MONTECARLO
    fTv0LambdaPxMC[ii]=-100000.;
    fTv0LambdaPyMC[ii]=-100000.;
    fTv0LambdaPzMC[ii]=-100000.;
    fTv0LambdaPDG[ii]=0;
    fTv0LambdaMotherPDG[ii]=0;
    fTv0LambdaMotherWeakPDG[ii]=0;
    fTv0LambdaOrigin[ii]=10;
#endif 
  }
  fTnv0=0;




  // Pre-track loop to resetglobaltrackreference and storeglobaltrackreference
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }



  //Lambdas
  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;
  AliAODEvent* aodEvt = dynamic_cast<AliAODEvent*>(fInputEvent);

  fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  for (int iv0 = 0;iv0< static_cast<TClonesArray *>(aodEvt->GetV0s())->GetEntriesFast();++iv0) {  
    Bool_t IsLambda = kFALSE;
    Bool_t IsAntiLambda = kFALSE;

    AliAODv0* v0aod = aodEvt->GetV0(iv0);
    fv0->Setv0(fInputEvent, v0aod);

    if (fLambda->isSelected(fv0)) {
      Lambdas.push_back(*fv0);
      IsLambda = kTRUE;
    }

    if(IsLambda){Fillv0(fv0, 1);}  

    if (fAntiLambda->isSelected(fv0)) {
      AntiLambdas.push_back(*fv0);
      IsAntiLambda = kTRUE;
    }

    if(IsAntiLambda){Fillv0(fv0, -1);}
  }




  //init tree Xi (Cascades)
  for(int ii=0;ii<MAXCASCADES;ii++){
   
   fTCascadeP[ii]=-100000.;
   fTCascadePx[ii]=-100000.;
   fTCascadePy[ii]=-100000.;
   fTCascadePz[ii]=-100000.;
   fTCascademT[ii]=-100000.;
   fTCascadeCharge[ii]=-10;
   fTCascadeDCA[ii]=-100000.;
   fTCascadeDaughtersDCA[ii]=-100000.;
//   fTCascadeXiMass[ii]=-100000.;
//   fTCascadeOmegaMass[ii]=-100000.;
   fTCascadeVx[ii]=-100000.;
   fTCascadeVy[ii]=-100000.;
   fTCascadeVz[ii]=-100000.;
   fTCascadeVr[ii]=-100000.;
   fTCascadePA[ii]=-100000.;
   //Lambda (form Cascade) 
   fTLambdaDCA[ii]=-100000.;
   fTLambdaDaughtersDCA[ii]=-100000.;
//   fTLambdaMass[ii]=-100000.;
//   fTLambdaK0Mass[ii]=-100000.;
   fTLambdaVx[ii]=-100000.;
   fTLambdaVy[ii]=-100000.;
   fTLambdaVz[ii]=-100000.;
   fTLambdaVr[ii]=-100000.;
   fTLambdaPA[ii]=-100000.;
   //Tracks (of cascade): 0 proton, 1 pion, 2 bachelor
   for(int jj=0;jj<3;jj++){
    fTTrackPx[ii][jj]=-100000.;
    fTTrackPy[ii][jj]=-100000.;
    fTTrackPz[ii][jj]=-100000.;
    fTTrackTPCmom[ii][jj]=-100000.;
//    fTTrackEta[ii][jj]=-100000.;
    fTTrackCharge[ii][jj]=-10;
    fTTrackDCA[ii][jj]=-100000.;
    fTTrackTPCsigma[ii][jj]=-100000.;
//    fTTrackTOFsigma[ii][jj]=-100000.;
    fTTrackNcl[ii][jj]=-100000;
/*
 * disk space for background
    fTTrackCrR[ii][jj]=-100000.;
    fTTrackCrF[ii][jj]=-100000.;
    fTTrackShared[ii][jj]=-100000;
    fTTrackTPCchi2[ii][jj]=-100000.;
    fTTrackTPConly[ii][jj]=kFALSE;
    fTTrackITScomplementary[ii][jj]=kFALSE;
    fTTrackITSpure[ii][jj]=kFALSE;
    fTTrackGLOBAL[ii][jj]=kFALSE;
    fTTrackFilterBit[ii][jj]=0;
*/
//    fTTrackPhi[ii][jj]=-100000.;
    fTTrackITStime[ii][jj]=kFALSE;
    fTTrackTOFtime[ii][jj]=kFALSE;
    fTTrackID[ii][jj]=-100000;
   
#ifdef MONTECARLO
     fTTrackPxMC[ii][jj]=-100000.;
     fTTrackPyMC[ii][jj]=-100000.;
     fTTrackPzMC[ii][jj]=-100000.;
     fTTrackPDG[ii][jj]=0;
     fTTrackMotherPDG[ii][jj]=0;
     fTTrackMotherWeakPDG[ii][jj]=0;
     fTTrackOrigin[ii][jj]=10;
#endif
   }

#ifdef MONTECARLO  
    fTCascadePxMC[ii]=-100000.;
    fTCascadePyMC[ii]=-100000.;
    fTCascadePzMC[ii]=-100000.;
    fTCascadePDG[ii]=0;
    fTCascadeMotherPDG[ii]=0;
    //fTCascadeMotherWeakPDG[ii]=0;
    fTCascadeOrigin[ii]=10;
#endif
  }
  fTnCascade=0;


  //Xi (bkg) and omegas (cascade loop)
  std::vector<AliFemtoDreamBasePart> Xis;       
  std::vector<AliFemtoDreamBasePart> AntiXis;
  std::vector<AliFemtoDreamBasePart> XisBGR;
  std::vector<AliFemtoDreamBasePart> AntiXisBGR;
  for (int iCasc = 0;iCasc< static_cast<TClonesArray *>(aodEvt->GetCascades())->GetEntriesFast();++iCasc) {
    Bool_t IsXiBGR = kFALSE;
    Bool_t IsAntiXiBGR = kFALSE;
    Bool_t IsXi = kFALSE;
    Bool_t IsAntiXi = kFALSE;

    AliAODcascade* casc = aodEvt->GetCascade(iCasc);
    fCascade->SetCascade(fInputEvent, casc);           

    if (fXiBGR->isSelected(fCascade)) {
      Xis.push_back(*fCascade);
      IsXiBGR = kTRUE;
    }
    if (fAntiXiBGR->isSelected(fCascade)) {
      AntiXisBGR.push_back(*fCascade);
      IsAntiXiBGR = kTRUE;
    }

    if (fXi->isSelected(fCascade)) {
      Xis.push_back(*fCascade);
      IsXi = kTRUE;
    }
    if (fAntiXi->isSelected(fCascade)) {
      AntiXis.push_back(*fCascade);
      IsAntiXi = kTRUE;
    }

    if(IsXiBGR||IsAntiXiBGR||IsXi||IsAntiXi) FillCascade(fCascade);   //

  }
  //fill Tree
  if(fTnv0>0&&fTnCascade>0) fGeorgiosTree->Fill(); //Select here events with at least on Xi and one Lambda
  //if(fTnv0>1) fGeorgiosTree->Fill();

  //pair cleaner
  fPairCleaner->ResetArray();

//  fPairCleaner->CleanTrackAndDecay(&Protons, &Lambdas, 0);   // to be set
//  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiLambdas, 1); // to be set

  fPairCleaner->CleanDecay(&Lambdas, 0);
  fPairCleaner->CleanDecay(&AntiLambdas, 1);
 
  fPairCleaner->StoreParticle(Lambdas);
  fPairCleaner->StoreParticle(AntiLambdas);


  fPairCleaner->CleanTrackAndDecay(&Lambdas, &XisBGR, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiLambdas, &AntiXisBGR, 1);
  fPairCleaner->CleanTrackAndDecay(&Lambdas, &Xis, 0); //this is apparently working
  fPairCleaner->CleanTrackAndDecay(&AntiLambdas, &AntiXis, 1); //this is apparently working

  fPairCleaner->CleanDecay(&Xis, 0);
  fPairCleaner->CleanDecay(&AntiXis, 1);
  fPairCleaner->CleanDecay(&XisBGR, 0);
  fPairCleaner->CleanDecay(&AntiXisBGR, 1);
  fPairCleaner->StoreParticle(Xis);
  fPairCleaner->StoreParticle(AntiXis);
  fPairCleaner->StoreParticle(XisBGR);
  fPairCleaner->StoreParticle(AntiXisBGR);


  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
  PostData(1, fEvtList);
  PostData(2, fLambdaList);
  PostData(3, fAntiLambdaList);
  PostData(4, fXiBGRList);
  PostData(5, fAntiXiBGRList);      //
  PostData(6, fXiList);
  PostData(7, fAntiXiList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fGeorgiosTree);

#ifdef MONTECARLO  
  if(fLambda->GetIsMonteCarlo()) {
   PostData(11, fLambdaMCList);
  }
  if(fAntiLambda->GetIsMonteCarlo()) {
   PostData(12, fAntiLambdaMCList);
  }
  if(fXiBGR->GetIsMonteCarlo()) {
   PostData(13, fXiBGRMCList);
  }
  if(fAntiXiBGR->GetIsMonteCarlo()) {
   PostData(14, fAntiXiBGRMCList);
  }
  if(fXi->GetIsMonteCarlo()) {
   PostData(15, fXiMCList);
  }
  if(fAntiXi->GetIsMonteCarlo()) {
   PostData(16, fAntiXiMCList);
  }
#endif
}

//____________________________________________________________________________________________________
void AliAnalysisTaskGeorgiosNTuple::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskGeorgiosNTuple::StoreGlobalTrackReference(AliVTrack *track) {
  // see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
	    return;
  }

  if (fGTI[trackID]) {
    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}






//________________________________________________________________________________
Bool_t AliAnalysisTaskGeorgiosNTuple::Fillv0(AliFemtoDreamv0 *Thev0, int Thev0Charge) {
 Bool_t Filled = kFALSE;

 fTv0DCA[fTnv0] = Thev0->GetDCAPrimVtx();
 fTv0Charge[fTnv0] = Thev0->GetCharge().at(0);
// fTv0LambdaMass[fTnv0] = Thev0->Getv0Mass();
 fTv0LambdaVr[fTnv0] = Thev0->GetTransverseRadius();
 fTv0LambdaPA[fTnv0] = Thev0->GetCPA();
 fTv0DaughtersDCA[fTnv0] = Thev0-> GetDaugDCA();

#ifdef MONTECARLO 
  TVector3 momMotherMC;
  momMotherMC = Thev0->GetMCMomentum();
  fTv0LambdaPxMC[fTnv0] = momMotherMC.X();
  fTv0LambdaPyMC[fTnv0] = momMotherMC.Y();
  fTv0LambdaPzMC[fTnv0] = momMotherMC.Z();
  fTv0LambdaPDG[fTnv0] = Thev0->GetMCPDGCode();
  fTv0LambdaMotherPDG[fTnv0] = Thev0->GetMotherPDG();
  fTv0LambdaMotherWeakPDG[fTnv0] = Thev0->GetMotherWeak();
  fTv0LambdaOrigin[fTnv0] =  Thev0->GetParticleOrigin();
#endif

 AliFemtoDreamTrack* TheTrackv0 = Thev0->GetPosDaughter();

 for(int jj=0;jj<2;jj++){
  if(jj==0){ //select the proton
   if(Thev0Charge==1) {TheTrackv0 = Thev0->GetPosDaughter();fTTrackv0DCA[fTnv0][jj]= Thev0->GetDCADaugPosVtx();}
   if(Thev0Charge==-1) {TheTrackv0 = Thev0->GetNegDaughter();fTTrackv0DCA[fTnv0][jj]= Thev0->GetDCADaugNegVtx();}
   fTTrackv0TPCsigma[fTnv0][jj]=(TheTrackv0->GetnSigmaTPC((int) (AliPID::kProton)));
   //fTTrackv0TOFsigma[fTnv0][jj]=(TheTrackv0->GetnSigmaTOF((int) (AliPID::kProton)));

  }
  else if(jj==1) { //select the pion
   if(Thev0Charge==1) {TheTrackv0 = Thev0->GetNegDaughter();fTTrackv0DCA[fTnv0][jj]= Thev0->GetDCADaugNegVtx();}
   if(Thev0Charge==-1) {TheTrackv0 = Thev0->GetPosDaughter();fTTrackv0DCA[fTnv0][jj]= Thev0->GetDCADaugPosVtx();}
   fTTrackv0TPCsigma[fTnv0][jj]=(TheTrackv0->GetnSigmaTPC((int) (AliPID::kPion)));
   //fTTrackv0TOFsigma[fTnv0][jj]=(TheTrackv0->GetnSigmaTOF((int) (AliPID::kPion)));

   }
  
#ifdef MONTECARLO
    TVector3 momMC;
    momMC = TheTrackv0->GetMCMomentum();
    fTTrackv0PxMC[fTnv0][jj] = momMC.X();
    fTTrackv0PyMC[fTnv0][jj] = momMC.Y();
    fTTrackv0PzMC[fTnv0][jj] = momMC.Z();
    fTTrackv0PDG[fTnv0][jj] = TheTrackv0->GetMCPDGCode();
    fTTrackv0MotherPDG[fTnv0][jj] = TheTrackv0->GetMotherPDG();
    fTTrackv0MotherWeakPDG[fTnv0][jj] = TheTrackv0->GetMotherWeak();
    fTTrackv0Origin[fTnv0][jj] = TheTrackv0->GetParticleOrigin();
#endif

   TVector3 mom;
   mom = TheTrackv0->GetMomentum();
   fTTrackv0Px[fTnv0][jj] = mom.X();
   fTTrackv0Py[fTnv0][jj] = mom.Y();
   fTTrackv0Pz[fTnv0][jj] = mom.Z();
//   fTTrackv0Eta[fTnv0][jj] = TheTrackv0->GetEta().at(0);
   fTTrackv0Charge[fTnv0][jj] = TheTrackv0->GetCharge().at(0);

   fTTrackv0Ncl[fTnv0][jj] = TheTrackv0->GetNClsTPC();
   //fTTrackv0CrF[fTnv0][jj] = TheTrackv0->GetRatioCr();
   //fTTrackv0CrR[fTnv0][jj] = TheTrackv0->GetTPCCrossedRows();
   //fTTrackv0FilterBit[fTnv0][jj] = TheTrackv0->GetFilterMap();

//   fTTrackv0Phi[fTnv0][jj] = (TheTrackv0->GetPhiAtRaidius().at(0)).at(0);//phi for r=85.cm ???
//   fTTrackv0ITStime[fTnv0][jj] = TheTrackv0->GetHasITSHit();
//   fTTrackv0TOFtime[fTnv0][jj] = TheTrackv0->GetTOFTimingReuqirement();
   fTTrackv0ID[fTnv0][jj] = TheTrackv0->GetIDTracks().at(0);
  }

  fTnv0++;
  Filled = kTRUE;
  return Filled;
 }


//________________________________________________________________________________
Bool_t AliAnalysisTaskGeorgiosNTuple::FillCascade(AliFemtoDreamCascade *TheCasc) {
 Bool_t Filled = kFALSE;

 fTCascadeCharge[fTnCascade] = TheCasc->GetCharge().at(0);
 fTCascadeDCA[fTnCascade] = TheCasc->GetDCAXiPrimVtx();
 fTCascadeDaughtersDCA[fTnCascade] = TheCasc->GetXiDCADaug();
// fTCascadeXiMass[fTnCascade] = TheCasc->GetXiMass();
// fTCascadeOmegaMass[fTnCascade] = TheCasc->GetOmegaMass();
 fTCascadeVr[fTnCascade] = TheCasc->GetXiTransverseRadius();
 fTCascadePA[fTnCascade] = TheCasc->GetCPA();

 fTLambdaDCA[fTnCascade] = TheCasc->Getv0DCAPrimVtx();
 fTLambdaDaughtersDCA[fTnCascade] = TheCasc->Getv0DCADaug();
// fTLambdaMass[fTnCascade] = TheCasc->Getv0Mass();
 fTLambdaVr[fTnCascade] = TheCasc->Getv0TransverseRadius();
 fTLambdaPA[fTnCascade] = TheCasc->Getv0CPA();

#ifdef MONTECARLO
  TVector3 momMotherMC;
  momMotherMC = TheCasc->GetMCMomentum();
  fTCascadePxMC[fTnCascade] = momMotherMC.X();
  fTCascadePyMC[fTnCascade] = momMotherMC.Y();
  fTCascadePzMC[fTnCascade] = momMotherMC.Z();
  fTCascadePDG[fTnCascade] = TheCasc->GetMCPDGCode();
  fTCascadeMotherPDG[fTnCascade] = TheCasc->GetMotherPDG();
  //fTCascadeMotherPDG[fTnCascade] = TheCasc->GetMotherWeak();
  fTCascadeOrigin[fTnCascade] =  TheCasc->GetParticleOrigin();
#endif

 AliFemtoDreamTrack* TheTrack = TheCasc->GetPosDaug();
 for(int jj=0;jj<3;jj++){
  if(jj==0){ //select the proton
   if(fTCascadeCharge[fTnCascade]==-1) {TheTrack = TheCasc->GetPosDaug();fTTrackDCA[fTnCascade][jj]= TheCasc->Getv0PosToPrimVtx();}
   if(fTCascadeCharge[fTnCascade]==1) {TheTrack = TheCasc->GetNegDaug();fTTrackDCA[fTnCascade][jj]= TheCasc->Getv0NegToPrimVtx();}
   fTTrackTPCsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kProton)));
/*
   fTTrackTOFsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kProton)));
*/
  }else if(jj==1) { //select the pion
   if(fTCascadeCharge[fTnCascade]==-1) {TheTrack = TheCasc->GetNegDaug();fTTrackDCA[fTnCascade][jj]= TheCasc->Getv0NegToPrimVtx();}
   if(fTCascadeCharge[fTnCascade]==1) {TheTrack = TheCasc->GetPosDaug();fTTrackDCA[fTnCascade][jj]= TheCasc->Getv0PosToPrimVtx();}
   fTTrackTPCsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kPion)));
/*  
   fTTrackTOFsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kPion)));
*/   
  }else if(jj==2) { //select the bachelor
   TheTrack = TheCasc->GetBach();
   fTTrackDCA[fTnCascade][jj]= TheCasc->BachDCAPrimVtx();
   fTTrackTPCsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTPC((int) (AliPID::kKaon)));
/*  
   fTTrackTOFsigma[fTnCascade][jj]=(TheTrack->GetnSigmaTOF((int) (AliPID::kKaon)));
*/   
  }

#ifdef MONTECARLO
    TVector3 momMC;
    momMC = TheTrack->GetMCMomentum();
    fTTrackPxMC[fTnCascade][jj] = momMC.X();
    fTTrackPyMC[fTnCascade][jj] = momMC.Y();
    fTTrackPzMC[fTnCascade][jj] = momMC.Z();
    fTTrackPDG[fTnCascade][jj] = TheTrack->GetMCPDGCode();
    fTTrackMotherPDG[fTnCascade][jj] = TheTrack->GetMotherPDG();
    fTTrackMotherWeakPDG[fTnCascade][jj] = TheTrack->GetMotherWeak();
    fTTrackOrigin[fTnCascade][jj] = TheTrack->GetParticleOrigin();
#endif

  TVector3 mom;
  mom = TheTrack->GetMomentum();
  fTTrackPx[fTnCascade][jj] = mom.X();
  fTTrackPy[fTnCascade][jj] = mom.Y();
  fTTrackPz[fTnCascade][jj] = mom.Z();
//  fTTrackEta[fTnCascade][jj] = TheTrack->GetEta().at(0);
  fTTrackCharge[fTnCascade][jj] = TheTrack->GetCharge().at(0);
  fTTrackNcl[fTnCascade][jj] = TheTrack->GetNClsTPC();
/*
 *
  fTTrackCrF[fTnCascade][jj] = TheTrack->GetRatioCr();
  fTTrackCrR[fTnCascade][jj] = TheTrack->GetTPCCrossedRows();
  fTTrackFilterBit[fTnCascade][jj] = TheTrack->GetFilterMap();
*/  
//  fTTrackPhi[fTnCascade][jj] = (TheTrack->GetPhiAtRaidius().at(0)).at(0);//phi for r=85.cm ???
  fTTrackITStime[fTnCascade][jj] = TheTrack->GetHasITSHit();
  fTTrackTOFtime[fTnCascade][jj] = TheTrack->GetTOFTimingReuqirement();
  fTTrackID[fTnCascade][jj] = TheTrack->GetIDTracks().at(0);
 }

 fTnCascade++;
 Filled = kTRUE;
 return Filled;
}


