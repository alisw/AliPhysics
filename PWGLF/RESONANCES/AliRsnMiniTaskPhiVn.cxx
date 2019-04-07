//
// Analysis task for 'mini' sub-package
// Contains all definitions needed for running an analysis:
// -- global event cut
// -- a list of track cuts (any number)
// -- definitions of output histograms
// -- values to be computed.
// Each one must be defined using the "CREATE" methods, which
// add directly a new element in the task collections, and don't
// need an external object to be passed to the task itself.
//

#include <Riostream.h>

#include <TH1.h>
#include <TList.h>
#include <TTree.h>
#include <TStopwatch.h>
#include "TRandom.h"
#include "TComplex.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "AliLog.h"
#include "AliEventplane.h"
#include "AliMultiplicity.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliRsnCutSet.h"
#include "AliRsnMiniPair.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniParticle.h"

#include "AliRsnMiniTaskPhiVn.h"
#include "AliMultSelection.h"


ClassImp(AliRsnMiniTaskPhiVn)

//__________________________________________________________________________________________________
AliRsnMiniTaskPhiVn::AliRsnMiniTaskPhiVn() :
AliAnalysisTaskSE(),
  fUseMC(kFALSE),
  fEvNum(0),
  fTriggerMask(0),
  fUseCentrality(kFALSE),
  fCentralityType("QUALITY"),
  fUseAOD049CentralityPatch(kFALSE),
  fUseCentralityPatchPbPb2011(0),
  fContinuousMix(kTRUE),
  fNMix(0),
  fMaxDiffMult(10),
  fMaxDiffVz(1.0),
  fMaxDiffAngle(1E20),
  fOutput(0x0),//new
  fHistograms("AliRsnMiniOutput", 0),
  fValues("AliRsnMiniValue", 0),
  fHEventStat(0x0),
  fHAEventsVsMulti(0x0),
  fHAEventsVsTracklets(0x0),
  fHQVectorPosReTest(0x0),
  fHQVectorPosImTest(0x0),
  fHQVectorNegReTest(0x0),
  fHQVectorNegImTest(0x0),
//fHPVectorTest(0x0),
  fHCentrality(0x0),
  fHMass1(0x0),
  fHMass2(0x0),
  fHMassPhi(0x0),
  fHTestPtMassCentrality(0x0),
  fHAEventVz(0x0),
  fHAEventMultiCent(0x0),
  fHAEventPlane(0x0),
  fSel1(0),
  fSel2(0),
  fEventCuts(0x0),
  fTrackCuts(0),
  fRsnEvent(),
  fEvBuffer(0x0),
  fTriggerAna(0x0),
  fESDtrackCuts(0x0),
  fMiniEvent(0x0),
  fMotherMass(0.0),
  fBigOutput(kFALSE),
  fMixPrintRefresh(-1),
  fMaxNDaughters(-1),
  fCheckP(kFALSE),
  fCheckFeedDown(kFALSE),   
  fOriginDselection(kFALSE),
  fKeepDfromB(kFALSE),
  fKeepDfromBOnly(kFALSE),
  fRejectIfNoQuark(kFALSE),
  fMotherAcceptanceCutMinPt(0.0),
  fMotherAcceptanceCutMaxEta(0.9),
  fKeepMotherInAcceptance(kFALSE),
  fnHarmToProcess(10)
{
//
// Dummy constructor ALWAYS needed for I/O.
//
   
for (int i = 0; i < fnHarmToProcess; ++i)//new
  {
fHQQVectorDenominator[i] = 0;
fHQQVectorCentralityNorm[i] = 0;
fHQQVectorCentrality[i] = 0;
fHPVectorPosTest[i] = 0;
fHPVectorNegTest[i] = 0;
      
for(int j = 0; j < 10; j++){
fHPposQnegCharged[i][j] = 0;
fHPposQnegChargedWeight[i][j] = 0;
fHPnegQposCharged[i][j] = 0;
fHPnegQposChargedWeight[i][j] = 0;

fHPposQneg[i][j] = 0;
fHPnegQpos[i][j] = 0;
fHPposQnegWeight[i][j] = 0;
fHPnegQposWeight[i][j] = 0;
}   
} 

fCharge[0] = fCharge[1] = 0;
fDaughter[0] = fDaughter[1] = AliRsnDaughter::kUnknown;
fUseStoredMass[0] = fUseStoredMass[1] = kFALSE;

     
}

  //__________________________________________________________________________________________________
  AliRsnMiniTaskPhiVn::AliRsnMiniTaskPhiVn(const char *name, Bool_t useMC) :
    AliAnalysisTaskSE(name),
    fUseMC(useMC),
    fEvNum(0),
    fTriggerMask(AliVEvent::kMB),
    fUseCentrality(kFALSE),
    fCentralityType("QUALITY"),
    fUseAOD049CentralityPatch(kFALSE),
    fUseCentralityPatchPbPb2011(0),
    fContinuousMix(kTRUE),
    fNMix(0),
    fMaxDiffMult(10),
    fMaxDiffVz(1.0),
    fMaxDiffAngle(1E20),
    fOutput(0x0),//new
    fHistograms("AliRsnMiniOutput", 0),
    fValues("AliRsnMiniValue", 0),
    fHEventStat(0x0),
    fHAEventsVsMulti(0x0),
    fHAEventsVsTracklets(0x0),
    fHQVectorPosReTest(0x0),
    fHQVectorPosImTest(0x0),
    fHQVectorNegReTest(0x0),
    fHQVectorNegImTest(0x0),
    //fHPVectorTest(0x0),
    fHCentrality(0x0),
    fHMass1(0x0),
    fHMass2(0x0),
    fHMassPhi(0x0),
    fHTestPtMassCentrality(0x0),
    fHAEventVz(0x0),
    fHAEventMultiCent(0x0),
    fHAEventPlane(0x0),
    fSel1(0),
    fSel2(0),   
    fEventCuts(0x0),
    fTrackCuts(0),
    fRsnEvent(),
    fEvBuffer(0x0),
    fTriggerAna(0x0),
    fESDtrackCuts(0x0),
    fMiniEvent(0x0),
    fMotherMass(0.0),
    fBigOutput(kFALSE),
    fMixPrintRefresh(-1),
    fMaxNDaughters(-1),
    fCheckP(kFALSE),
    fCheckFeedDown(kFALSE),   
    fOriginDselection(kFALSE),
    fKeepDfromB(kFALSE),
    fKeepDfromBOnly(kFALSE),
    fRejectIfNoQuark(kFALSE),
    fMotherAcceptanceCutMinPt(0.0),
    fMotherAcceptanceCutMaxEta(0.9),
    fKeepMotherInAcceptance(kFALSE),
    fnHarmToProcess(10) 
  {
//
// Default constructor.
// Define input and output slots here (never in the dummy constructor)
// Input slot #0 works with a TChain - it is connected to the default input container
// Output slot #1 writes into a TH1 container
//

    DefineOutput(1, TList::Class());
   
    for (int i = 0; i < fnHarmToProcess; ++i)//new
      {
	fHQQVectorDenominator[i] = 0;
	fHQQVectorCentralityNorm[i] = 0;
	fHQQVectorCentrality[i] = 0;
	fHPVectorPosTest[i] = 0;
	fHPVectorNegTest[i] = 0;
      
	for(int j = 0; j < 10; j++){
	  fHPposQnegCharged[i][j] = 0;
	  fHPposQnegChargedWeight[i][j] = 0;
	  fHPnegQposCharged[i][j] = 0;
	  fHPnegQposChargedWeight[i][j] = 0;
         
	  fHPposQneg[i][j] = 0;
	  fHPnegQpos[i][j] = 0;
	  fHPposQnegWeight[i][j] = 0;
	  fHPnegQposWeight[i][j] = 0;
	}
      } 

    fCharge[0] = fCharge[1] = 0;
    fDaughter[0] = fDaughter[1] = AliRsnDaughter::kUnknown;
    fUseStoredMass[0] = fUseStoredMass[1] = kFALSE;    
  }

//__________________________________________________________________________________________________
AliRsnMiniTaskPhiVn::AliRsnMiniTaskPhiVn(const AliRsnMiniTaskPhiVn &copy) :
  AliAnalysisTaskSE(copy),
  fUseMC(copy.fUseMC),
  fEvNum(0),
  fTriggerMask(copy.fTriggerMask),
  fUseCentrality(copy.fUseCentrality),
  fCentralityType(copy.fCentralityType),
  fUseAOD049CentralityPatch(copy.fUseAOD049CentralityPatch),
  fUseCentralityPatchPbPb2011(copy.fUseCentralityPatchPbPb2011),
  fContinuousMix(copy.fContinuousMix),
  fNMix(copy.fNMix),
  fMaxDiffMult(copy.fMaxDiffMult),
  fMaxDiffVz(copy.fMaxDiffVz),
  fMaxDiffAngle(copy.fMaxDiffAngle),
  fOutput(0x0),
  fHistograms(copy.fHistograms),
  fValues(copy.fValues),
  fHEventStat(0x0),
  fHAEventsVsMulti(0x0),
  fHAEventsVsTracklets(0x0),
  fHQVectorPosReTest(0x0),
  fHQVectorPosImTest(0x0),
  fHQVectorNegReTest(0x0),
  fHQVectorNegImTest(0x0),
  //fHPVectorTest(0x0),
  fHCentrality(0x0),
  fHMass1(0x0),
  fHMass2(0x0),
  fHMassPhi(0x0),
  fHTestPtMassCentrality(0x0),
  fHAEventVz(0x0),
  fHAEventMultiCent(0x0),
  fHAEventPlane(0x0),
  fSel1(0),
  fSel2(0),
  fEventCuts(copy.fEventCuts),
  fTrackCuts(copy.fTrackCuts),
  fRsnEvent(),
  fEvBuffer(0x0),
  fTriggerAna(copy.fTriggerAna),
  fESDtrackCuts(copy.fESDtrackCuts),
  fMiniEvent(0x0),
  fMotherMass(copy.fMotherMass),
  fBigOutput(copy.fBigOutput),
  fMixPrintRefresh(copy.fMixPrintRefresh),
  fMaxNDaughters(copy.fMaxNDaughters),
  fCheckP(copy.fCheckP),
  fCheckFeedDown(copy.fCheckFeedDown),   
  fOriginDselection(copy.fOriginDselection),
  fKeepDfromB(copy.fOriginDselection),
  fKeepDfromBOnly(copy.fKeepDfromBOnly),
  fRejectIfNoQuark(copy.fRejectIfNoQuark),
  fMotherAcceptanceCutMinPt(copy.fMotherAcceptanceCutMinPt),
  fMotherAcceptanceCutMaxEta(copy.fMotherAcceptanceCutMaxEta),
  fKeepMotherInAcceptance(copy.fKeepMotherInAcceptance),
  fnHarmToProcess(copy.fnHarmToProcess) 
{
  //
  // Copy constructor.
  // Implemented as requested by C++ standards.
  // Can be used in PROOF and by plugins.
  //
   
  for (int i = 0; i < fnHarmToProcess; ++i)//new
    {
      fHQQVectorDenominator[i] = 0;
      fHQQVectorCentralityNorm[i] = 0;
      fHQQVectorCentrality[i] = 0;
      fHPVectorPosTest[i] = 0;
      fHPVectorNegTest[i] = 0;
      
      for(int j = 0; j < 10; j++){
	fHPposQnegCharged[i][j] = 0;
	fHPposQnegChargedWeight[i][j] = 0;
	fHPnegQposCharged[i][j] = 0;
	fHPnegQposChargedWeight[i][j] = 0;
         
	fHPposQneg[i][j] = 0;
	fHPnegQpos[i][j] = 0;
	fHPposQnegWeight[i][j] = 0;
	fHPnegQposWeight[i][j] = 0;
      }
    } 

  fCharge[0] = fCharge[1] = 0;
  fDaughter[0] = fDaughter[1] = AliRsnDaughter::kUnknown; 
  fUseStoredMass[0] = fUseStoredMass[1] = kFALSE;     
}

//__________________________________________________________________________________________________
AliRsnMiniTaskPhiVn &AliRsnMiniTaskPhiVn::operator=(const AliRsnMiniTaskPhiVn &copy)
{
  //
  // Assignment operator.
  // Implemented as requested by C++ standards.
  // Can be used in PROOF and by plugins.
  //

  AliAnalysisTaskSE::operator=(copy);
  if (this == &copy)
    return *this;
  fUseMC = copy.fUseMC;
  fEvNum = copy.fEvNum;
  fTriggerMask = copy.fTriggerMask;
  fUseCentrality = copy.fUseCentrality;
  fCentralityType = copy.fCentralityType;
  fUseAOD049CentralityPatch = copy.fUseAOD049CentralityPatch;
  fUseCentralityPatchPbPb2011 = copy.fUseCentralityPatchPbPb2011;
  fContinuousMix = copy.fContinuousMix;
  fNMix = copy.fNMix;
  fMaxDiffMult = copy.fMaxDiffMult;
  fMaxDiffVz = copy.fMaxDiffVz;
  fMaxDiffAngle = copy.fMaxDiffAngle;
  fHistograms = copy.fHistograms;
  fValues = copy.fValues;
  fHEventStat = copy.fHEventStat;
  fHAEventsVsMulti = copy.fHAEventsVsMulti;
  fHAEventsVsTracklets = copy.fHAEventsVsTracklets;
  fHQVectorPosReTest = copy.fHQVectorPosReTest;
  fHQVectorPosImTest = copy.fHQVectorPosImTest;
  fHQVectorNegReTest = copy.fHQVectorNegReTest;
  fHQVectorNegImTest = copy.fHQVectorNegImTest;
  //fHPVectorTest = copy.fHPVectorTest;   
  fHCentrality = copy.fHCentrality;
  fHMass1 = copy.fHMass1;
  fHMass2 = copy.fHMass2;
  fHMassPhi = copy.fHMassPhi;
  fHTestPtMassCentrality = copy.fHTestPtMassCentrality;
  fHAEventVz = copy.fHAEventVz;
  fHAEventMultiCent = copy.fHAEventMultiCent;
  fHAEventPlane = copy.fHAEventPlane;
  fSel1.Set(0);
  fSel2.Set(0);
  fEventCuts = copy.fEventCuts;
  fTrackCuts = copy.fTrackCuts;
  fTriggerAna = copy.fTriggerAna;
  fESDtrackCuts = copy.fESDtrackCuts;
  fMotherMass = copy.fMotherMass;
  fBigOutput = copy.fBigOutput;
  fMixPrintRefresh = copy.fMixPrintRefresh;
  fMaxNDaughters = copy.fMaxNDaughters;
  fCheckP = copy.fCheckP;
  fCheckFeedDown = copy.fCheckFeedDown;
  fOriginDselection = copy.fOriginDselection;
  fKeepDfromB = copy.fOriginDselection;
  fKeepDfromBOnly = copy.fKeepDfromBOnly;
  fRejectIfNoQuark = copy.fRejectIfNoQuark;
  fMotherAcceptanceCutMinPt = copy.fMotherAcceptanceCutMinPt;
  fMotherAcceptanceCutMaxEta = copy.fMotherAcceptanceCutMaxEta;
  fKeepMotherInAcceptance = copy.fKeepMotherInAcceptance;
  fnHarmToProcess = copy.fnHarmToProcess;
  return (*this);

  Int_t i;
  for (i = 0; i < 2; i++) {
    fDaughter[i] = copy.fDaughter[i];
    fCharge[i] = copy.fCharge[i];
    fUseStoredMass[i] = copy.fUseStoredMass[i];      
  }   
}

//__________________________________________________________________________________________________
AliRsnMiniTaskPhiVn::~AliRsnMiniTaskPhiVn()
{
  //
  // Destructor.
  // Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  //


  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    delete fEvBuffer;
  }
}

//__________________________________________________________________________________________________
Int_t AliRsnMiniTaskPhiVn::AddTrackCuts(AliRsnCutSet *cuts)
{
  //
  // Add a new cut set for a new criterion for track selection.
  // A user can add as many as he wants, and each one corresponds
  // to one of the available bits in the AliRsnMiniParticle mask.
  // The only check is the following: if a cut set with the same name
  // as the argument is there, this is not added.
  // Return value is the array position of this set.
  //

  TObject *obj = fTrackCuts.FindObject(cuts->GetName());

  if (obj) {
    AliInfo(Form("A cut set named '%s' already exists", cuts->GetName()));
    return fTrackCuts.IndexOf(obj);
  } else {
    fTrackCuts.AddLast(cuts);
    return fTrackCuts.IndexOf(cuts);
  }
}

//__________________________________________________________________________________________________
void AliRsnMiniTaskPhiVn::UserCreateOutputObjects()
{
  //
  // Initialization of outputs.
  // This is called once per worker node.
  //

  // reset counter
  fEvNum = -1;

  // message
  AliInfo(Form("Selected event characterization: %s (%s)", (fUseCentrality ? "centrality" : "multiplicity"), fCentralityType.Data()));

  // initialize trigger analysis
  if (fTriggerAna) delete fTriggerAna;
  fTriggerAna = new AliTriggerAnalysis;

  // initialize ESD quality cuts
  if (fESDtrackCuts) delete fESDtrackCuts;
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

  // create list and set it as owner of its content (MANDATORY)
  if (fBigOutput) OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

  //new
  fHQVectorPosReTest = new TH1F("QVector Pos Re","QVector Pos Re",1000, -1.0, 1.0);
  fOutput->Add(fHQVectorPosReTest);

  fHQVectorPosImTest = new TH1F("QVector Pos Im","QVector Pos Im",1000, -1.0, 1.0);
  fOutput->Add(fHQVectorPosImTest);

  fHQVectorNegReTest = new TH1F("QVector Neg Re","QVector Neg Re",1000, -1.0, 1.0);
  fOutput->Add(fHQVectorNegReTest);

  fHQVectorNegImTest = new TH1F("QVector Neg Im","QVector Neg Im",1000, -1.0, 1.0);
  fOutput->Add(fHQVectorNegImTest);

  //fHPVectorTest = new TH1F("PVector Test", "PVector Test", 1000, -1000., 1000.);
  //fOutput->Add(fHPVectorTest);

  fHCentrality = new TH1F("Centrality", "Centrality", 100, 0., 100.);
  fOutput->Add(fHCentrality);

  fHMass1 = new TH1F("Mass1", "Mass1", 1000, 0., 10.);
  fOutput->Add(fHMass1);

  fHMass2 = new TH1F("Mass2", "Mass2", 1000, 0., 10.);
  fOutput->Add(fHMass2);

  fHMassPhi = new TH1F("MassPhi", "MassPhi", 1000, 0., 10.);
  fOutput->Add(fHMassPhi);

  fHTestPtMassCentrality = new TH3F("PtMassCentrality", "PtMassCentrality", 100, 0., 9.9, 1000, 0.0, 1.998, 10, 0., 100.);
  fOutput->Add(fHTestPtMassCentrality);

  //initialize flow mode histos -->> Segmentation Violation wegen fnHarmToProcess, mit [10] scheint es zu funktionieren. Was da los?

  Int_t dCbins[11] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

  for (Int_t harm=0;harm<1;harm++) {


    fHQQVectorDenominator[harm] = new TH1F(Form("QxQ* of v%d",harm+2),Form("QxQ* of v%d",harm+2),1000, -1.0, 1.0);
    fOutput->Add(fHQQVectorDenominator[harm]);

    fHQQVectorCentralityNorm[harm] = new TProfile(Form("QxQ* Weighted vs Centrality of v%d",harm+2), Form("QxQ* Weighted vs Centrality of v%d",harm+2), 10, 0., 100., -1., 1.);
    fOutput->Add(fHQQVectorCentralityNorm[harm]);

    fHQQVectorCentrality[harm] = new TProfile(Form("QxQ* vs Centrality of v%d",harm+2), Form("QxQ* vs Centrality of v%d",harm+2), 10, 0., 100., -1., 1.);
    fOutput->Add(fHQQVectorCentrality[harm]);  

    fHPVectorPosTest[harm] = new TH1F(Form("PPos Vector Charged of v%d",harm+2), Form("PPos Vector Charged of v%d",harm+2), 1000, -1, 1);
    fOutput->Add(fHPVectorPosTest[harm]);

    fHPVectorNegTest[harm] =  new TH1F(Form("PNeg Vector Charged of v%d",harm+2), Form("PNeg Vector Charged of v%d",harm+2), 1000, -1, 1);
    fOutput->Add(fHPVectorNegTest[harm]);         

    for(Int_t c=0;c<10;c++){

      fHPposQnegCharged[harm][c] = new TProfile(Form("Ppos Qneg Charged of v%d in %d to %d",harm+2, dCbins[c], dCbins[c+1]), Form("Ppos Qneg Charged of v%d",harm+2), 100, 0., 9.9, -1, 1);
      fOutput->Add(fHPposQnegCharged[harm][c]);

      fHPposQnegChargedWeight[harm][c] = new TProfile(Form("Ppos Qneg Charged Weighted of v%d in %d to %d",harm+2, dCbins[c], dCbins[c+1]), Form("Ppos Qneg Charged Weighted of v%d",harm+2), 100, 0., 9.9, -1, 1);
      fOutput->Add(fHPposQnegChargedWeight[harm][c]);

      fHPnegQposCharged[harm][c] = new TProfile(Form("Pneg Qpos Charged of v%d in %d to %d",harm+2, dCbins[c], dCbins[c+1]), Form("Pneg Qpos Charged of v%d",harm+2), 100, 0., 9.9, -1, 1);
      fOutput->Add(fHPnegQposCharged[harm][c]);

      fHPnegQposChargedWeight[harm][c] =  new TProfile(Form("Pneg Qpos Charged Weighted of v%d in %d to %d",harm+2, dCbins[c], dCbins[c+1]), Form("Pneg Qpos Charged Weighted of v%d",harm+2), 100, 0., 9.9, -1, 1);
      fOutput->Add(fHPnegQposChargedWeight[harm][c]);

      fHPposQneg[harm][c] = new TProfile2D(Form("Ppos Qneg of v%d in %d to %d",harm+2, dCbins[c], dCbins[c+1]), Form("Ppos Qneg of v%d",harm+2), 100, 0., 9.9, 1000, 0.0, 1.998, -1, 1);
      fOutput->Add(fHPposQneg[harm][c]);

      fHPnegQpos[harm][c] = new TProfile2D(Form("Pneg Qpos of v%d in %d to %d",harm+2, dCbins[c], dCbins[c+1]), Form("Pneg Qpos of v%d",harm+2), 100, 0., 9.9, 1000, 0.0, 1.998, -1, 1);
      fOutput->Add(fHPnegQpos[harm][c]);

      fHPposQnegWeight[harm][c] = new TProfile2D(Form("Ppos Qneg Weighted of v%d in %d to %d",harm+2, dCbins[c], dCbins[c+1]), Form("Ppos Qneg Weighted of v%d",harm+2), 100, 0., 9.9, 1000, 0.0, 1.998, -1, 1);
      fOutput->Add(fHPposQnegWeight[harm][c]);

      fHPnegQposWeight[harm][c] = new TProfile2D(Form("Pneg Qpos Weighted of v%d in %d to %d",harm+2, dCbins[c], dCbins[c+1]), Form("Pneg Qpos Weighted of v%d",harm+2), 100, 0., 9.9, 1000, 0.0, 1.998, -1, 1);
      fOutput->Add(fHPnegQposWeight[harm][c]);
    }
  }
 
  //end new

  // initialize event statistics counter
  fHEventStat = new TH1F("hEventStat", "Event statistics", 8, 0.0, 8.0);
  fHEventStat->GetXaxis()->SetBinLabel(1, "CINT1B");
  fHEventStat->GetXaxis()->SetBinLabel(2, "V0AND");
  fHEventStat->GetXaxis()->SetBinLabel(3, "Candle");
  fHEventStat->GetXaxis()->SetBinLabel(4, "Accepted");
  fHEventStat->GetXaxis()->SetBinLabel(5, "Not Accepted - Total");
  fHEventStat->GetXaxis()->SetBinLabel(6, "Not Accepted - No Track Vertex");
  fHEventStat->GetXaxis()->SetBinLabel(7, "Not Accepted - Not Enough Contributors");
  fHEventStat->GetXaxis()->SetBinLabel(8, "Not Accepted - No Vertex inside |z| < 10 cm");
   
  fOutput->Add(fHEventStat);

  fHAEventsVsMulti = new TH1F("hAEventsVsMulti", "Accepted events vs Centrality", 100, 0, 100.0);
  fOutput->Add(fHAEventsVsMulti);
   
  fHAEventsVsTracklets = new TH1F("hAEventsVsTracklets", "Accepted events vs Tracklet Number",1000, 0, 1000.0);
  fOutput->Add(fHAEventsVsTracklets);

  if(fHAEventVz) fOutput->Add(fHAEventVz);
  if(fHAEventMultiCent) fOutput->Add(fHAEventMultiCent);
  if(fHAEventPlane) fOutput->Add(fHAEventPlane);

  TIter next(&fTrackCuts);
  AliRsnCutSet *cs;
  while ((cs = (AliRsnCutSet *) next())) {
    cs->Init(fOutput);
  }

  // create temporary tree for filtered events
  if (fMiniEvent) SafeDelete(fMiniEvent);
  fEvBuffer = new TTree("EventBuffer", "Temporary buffer for mini events");
  fMiniEvent = new AliRsnMiniEvent();
  fEvBuffer->Branch("events", "AliRsnMiniEvent", &fMiniEvent);

  // create one histogram per each stored definition (event histograms)
  Int_t i, ndef = fHistograms.GetEntries();
  AliRsnMiniOutput *def = 0x0;
  for (i = 0; i < ndef; i++) {
    def = (AliRsnMiniOutput *)fHistograms[i];
    if (!def) continue;
    if (!def->Init(GetName(), fOutput)) {
      AliError(Form("Def '%s': failed initialization", def->GetName()));
      continue;
    }
  }

  // post data for ALL output slots >0 here, to get at least an empty histogram
  PostData(1, fOutput);
}
//_______________________________________________________________________________
void AliRsnMiniTaskPhiVn::UserExec(Option_t *)
{
  //
  // Computation loop.
  // In this case, it checks if the event is acceptable, and eventually
  // creates the corresponding mini-event and stores it in the buffer.
  // The real histogram filling is done at the end, in "FinishTaskOutput".
  //

  // increment event counter
  fEvNum++;

  // check current event
  Char_t check = CheckCurrentEvent();
  if (!check) return;


  // setup PID response
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler *)man->GetInputEventHandler();
  fRsnEvent.SetPIDResponse(inputHandler->GetPIDResponse());

  // fill a mini-event from current
  // and skip this event if no tracks were accepted
  FillMiniEvent(check);

  //new: loop over all tracks per event apply filterbit and calculate QVector   
  Double_t dEtaGap = 0.4;
  //Double_t dCentralityBins[11] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
  Int_t iFlowNumHarmonicsMax = 3;
  Double_t iHarm = 2.;                                                                 // iHarm = harmonic (n = 2,3,4,5,...)
  const Int_t iDiffHarm = iFlowNumHarmonicsMax - int(iHarm);
   
  TComplex FlowVecQneg[iDiffHarm];
  TComplex FlowVecQpos[iDiffHarm];
  TComplex FlowVecQProduct[iDiffHarm];

  Double_t    MQPos[iDiffHarm];
  Double_t    MQNeg[iDiffHarm];

  Double_t FlowVecQnegRe[iDiffHarm];
  Double_t FlowVecQnegIm[iDiffHarm];
  Double_t FlowVecQposRe[iDiffHarm];
  Double_t FlowVecQposIm[iDiffHarm];

  for (Int_t i=0; i<iDiffHarm; i++){
    FlowVecQneg[i] = TComplex(0,0,kFALSE);
    FlowVecQpos[i] = TComplex(0,0,kFALSE);
    FlowVecQProduct[i] = TComplex(0,0,kFALSE);
    MQPos[i] = 0.0;
    MQNeg[i] = 0.0;
    FlowVecQnegRe[i] = 0.0;
    FlowVecQnegIm[i] = 0.0;
    FlowVecQposRe[i] = 0.0;
    FlowVecQposIm[i] = 0.0;
  }
 
  if(!fInputEvent) return;                                                                // if the pointer to the event is empty (getting it failed) skip this event

  Int_t iTracks(fInputEvent->GetNumberOfTracks());                                        // see how many tracks there are in the event

  Float_t lPercentile = 300; 
  AliMultSelection *MultSelection = 0x0; 
  MultSelection = (AliMultSelection * ) fInputEvent->FindListObject("MultSelection");
  if( !MultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  }else{
    //This is the Centrality
    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }

  if(lPercentile > 100.) return;

  //calculating index of centrality
  Int_t icentralityIndex = CentralityIndex(lPercentile, 100.0, 0.0, 10.);


  for(int j = 0; j < iDiffHarm; j++) {

    for(Int_t i(0); i < iTracks; i++) {                                           // loop over all these tracks
         
      AliAODTrack* track = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));         // get a track (type AliAODTrack) from the event
         
      if(!track || !track->TestFilterBit(5)) continue;                           // 5 for all charged particles
            
      if (track->Pt() < 0.2 || track->Pt() > 5.0)
	{
	  continue;
	}
      Double_t weights = 1.;
      if(track->Eta() > dEtaGap / 2){   
	//RFP in positive eta acceptance
	Double_t dQcosPos = weights*TMath::Cos(iHarm * track->Phi());
	Double_t dQsinPos = weights*TMath::Sin(iHarm * track->Phi());
               
	Double_t FlowQPosRe = TComplex(dQcosPos,dQsinPos,kFALSE).Re();
	FlowVecQposRe[j] += FlowQPosRe;

	Double_t FlowQPosIm = TComplex(dQcosPos,dQsinPos,kFALSE).Im();
	FlowVecQposIm[j] += FlowQPosIm;

	Double_t FlowQPos = TComplex(dQcosPos,dQsinPos,kFALSE);
	FlowVecQpos[j] += FlowQPos;         
 
	MQPos[j] += TComplex(1.,0.,kFALSE).Re();
      }
         
      if (track->Eta() < -dEtaGap / 2) {
	//RFP in negative eta acceptance
	Double_t dQcosNeg = weights*TMath::Cos(iHarm * track->Phi());
	Double_t dQsinNeg = weights*TMath::Sin(iHarm * track->Phi());

	Double_t FlowQNegRe = TComplex(dQcosNeg,dQsinNeg,kFALSE).Re();
	FlowVecQnegRe[j] += FlowQNegRe;

	Double_t FlowQNegIm = TComplex(dQcosNeg,dQsinNeg,kFALSE).Im();
	FlowVecQnegIm[j] += FlowQNegIm;

	Double_t FlowQNeg = TComplex(dQcosNeg,dQsinNeg,kFALSE);
	FlowVecQneg[j] += FlowQNeg;
               
	MQNeg[j] += TComplex(1.,0.,kFALSE).Re();
      }
    }

    Double_t weightingQ = MQNeg[j]*MQPos[j];

    FlowVecQProduct[j] = FlowVecQpos[j]*TComplex::Conjugate(FlowVecQneg[j]);

    //Filling
    if (weightingQ != 0)
      {
	fHQVectorPosReTest->Fill(FlowVecQposRe[j]/weightingQ);
	fHQVectorPosImTest->Fill(FlowVecQposIm[j]/weightingQ);

	fHQVectorNegReTest->Fill(FlowVecQnegRe[j]/weightingQ);
	fHQVectorNegImTest->Fill(FlowVecQnegIm[j]/weightingQ);
      }

    fHCentrality->Fill(lPercentile);

    fHQQVectorDenominator[j]->Fill((FlowVecQProduct[j].Re())/weightingQ);
      
    if (weightingQ != 0)
      {
	fHQQVectorCentralityNorm[j]->Fill(lPercentile,(FlowVecQProduct[j].Re())/weightingQ, weightingQ);
	fHQQVectorCentrality[j]->Fill(lPercentile,(FlowVecQProduct[j].Re())/weightingQ);  
      }

      
    //PVector
    //AliRsnMiniEvent *event1 = fMiniEvent; 
    //AliRsnMiniEvent *event2 = fMiniEvent; 
    Double_t weightsP = 1.;
    Double_t weightsN = 1.;
    const Int_t iPtBins = 100;
    const Int_t iInvMass = 1000;
    Double_t dPtBins[iPtBins];

    Double_t MPPos[iDiffHarm][iPtBins][iInvMass];
    Double_t MPNeg[iDiffHarm][iPtBins][iInvMass];
    Double_t MPPosCharged[iDiffHarm][iPtBins];
    Double_t MPNegCharged[iDiffHarm][iPtBins];
    TComplex FlowPVecPosCharged[iDiffHarm][iPtBins];
    TComplex FlowPVecNegCharged[iDiffHarm][iPtBins];
    TComplex PVecPosCharged[iDiffHarm][iPtBins];
    TComplex PVecNegCharged[iDiffHarm][iPtBins];

      
    for (Int_t i=0; i<iDiffHarm; i++){
      for (Int_t j=0; j<iPtBins; j++){
	for (Int_t k=0; k<iInvMass; k++){
	  MPPos[i][j][k] = 0.0;
	  MPNeg[i][j][k] = 0.0;

	  MPPosCharged[i][j] = 0.0;
	  MPNegCharged[i][j] = 0.0;

	  FlowPVecPosCharged[i][j] = TComplex(0,0,kFALSE);
	  FlowPVecNegCharged[i][j] = TComplex(0,0,kFALSE);
	  PVecPosCharged[i][j] = TComplex(0,0,kFALSE);
	  PVecNegCharged[i][j] = TComplex(0,0,kFALSE); 
	}
      }
    }
 


    for (int i = 0; i < 100; ++i)
      {dPtBins[i] = double(i)/10;} //0.0 - 9.9

    Double_t dInvMassBins[iInvMass];
   
    for(int b = 0; b < iInvMass; ++b)
      {dInvMassBins[b] = 2 * double(b)/1000;} // 0.0 - 1.998   



    for(Int_t tr(0); tr < iTracks; tr++) {                                           // loop over all these tracks
    
      AliAODTrack* tracktr = static_cast<AliAODTrack*>(fInputEvent->GetTrack(tr));    // get a track (type AliAODTrack) from the event
         
      if(!tracktr || !tracktr->TestFilterBit(5)) continue;                              // 5 for all charged particles

      Double_t dPtCharged = tracktr->Pt();
      Double_t dEtaCharged = tracktr->Eta();
      Double_t dPhiCharged = tracktr->Phi();
      //Double_t dEnergyCharged = tracktr->E();
      //Double_t dPCharged = tracktr->P();


      if (dPtCharged > 9.9 || dPtCharged < 0.0)
	{
	  continue;
	}


      Int_t ptIndexCharged = PtIndex(dPtCharged,9.9, 0.0, 100.);


      if(dEtaCharged > dEtaGap / 2){   
	//RFP in positive eta acceptance
	Double_t dPCosPosCharged = TMath::Cos(iHarm * dPhiCharged); 
	Double_t dPSinPosCharged = TMath::Sin(iHarm * dPhiCharged);
	FlowPVecPosCharged[j][ptIndexCharged] += TComplex(weightsP*dPCosPosCharged,weightsP*dPSinPosCharged,kFALSE);
	MPPosCharged[j][ptIndexCharged] += TComplex(1.,0.,kFALSE).Re();
	PVecPosCharged[j][ptIndexCharged] += TComplex(weightsP*TMath::Cos(iHarm * dPhiCharged),weightsP*TMath::Sin(iHarm * dPhiCharged),kFALSE);
      }
         
      if (dEtaCharged < -dEtaGap / 2) {
	//RFP in negative eta acceptance   
	Double_t dPCosNegCharged = TMath::Cos(iHarm * dPhiCharged); 
	Double_t dPSinNegCharged = TMath::Sin(iHarm * dPhiCharged);
	FlowPVecNegCharged[j][ptIndexCharged] += TComplex(weightsN*dPCosNegCharged,weightsN*dPSinNegCharged,kFALSE);
	MPNegCharged[j][ptIndexCharged] += TComplex(1.,0.,kFALSE).Re(); 
	PVecNegCharged[j][ptIndexCharged] += TComplex(weightsP*TMath::Cos(iHarm * dPhiCharged),weightsP*TMath::Sin(iHarm * dPhiCharged),kFALSE);
      }
    }

    TComplex tPosNegCharged = 999;
    TComplex tNegPosCharged = 999;

    for (int pCharged = 0; pCharged < 100; ++pCharged)
      { 
	if (PVecPosCharged[j][pCharged].Re() != 0.)
	  {
            fHPVectorPosTest[j]->Fill(PVecPosCharged[j][pCharged].Re()/MPPosCharged[j][pCharged]);
	  }
            
	if (PVecNegCharged[j][pCharged].Re() != 0.)
	  {
            fHPVectorNegTest[j]->Fill(PVecNegCharged[j][pCharged].Re()/MPNegCharged[j][pCharged]);
	  }

	tPosNegCharged = FlowPVecPosCharged[j][pCharged]*TComplex::Conjugate(FlowVecQneg[j]);
	tNegPosCharged = FlowPVecNegCharged[j][pCharged]*TComplex::Conjugate(FlowVecQpos[j]);
        
	Double_t weightingPposQnegCharged = MPPosCharged[j][pCharged]*MQNeg[j];
	Double_t weightingPnegQposCharged = MPNegCharged[j][pCharged]*MQPos[j];

	if (weightingPposQnegCharged == 0. || weightingPnegQposCharged == 0.)
	  {
            continue;
	  }

	if(TMath::Abs(tPosNegCharged.Re()/weightingPposQnegCharged) < 1)
	  {   
            fHPposQnegCharged[j][icentralityIndex]->Fill(dPtBins[pCharged],tPosNegCharged.Re()/weightingPposQnegCharged);
            fHPposQnegChargedWeight[j][icentralityIndex]->Fill(dPtBins[pCharged], tPosNegCharged.Re()/weightingPposQnegCharged, weightingPposQnegCharged);
	  }

	if(TMath::Abs(tNegPosCharged.Re()/weightingPnegQposCharged) < 1)
	  { 
            fHPnegQposCharged[j][icentralityIndex]->Fill(dPtBins[pCharged], tNegPosCharged.Re()/weightingPnegQposCharged);
            fHPnegQposChargedWeight[j][icentralityIndex]->Fill(dPtBins[pCharged], tNegPosCharged.Re()/weightingPnegQposCharged, weightingPnegQposCharged);
	  }
      }


    //TComplex PVector[100];
    TComplex FlowPVecPos[iDiffHarm][iPtBins][iInvMass];
    TComplex FlowPVecNeg[iDiffHarm][iPtBins][iInvMass];

    for (Int_t i=0; i<iDiffHarm; i++){
      for (Int_t j=0; j<iPtBins; j++){
	for (Int_t k=0; k<iInvMass; k++){
	  FlowPVecPos[i][j][k] = TComplex(0,0,kFALSE);
	  FlowPVecNeg[i][j][k] = TComplex(0,0,kFALSE);
	}
      }
    }
 


    AliRsnMiniPair *dummyPair = new AliRsnMiniPair();

    // loop variables
    Int_t i1, i2, start;
    AliRsnMiniParticle *p1, *p2;
    Double_t mass1, mass2;

    // it is necessary to know if criteria for the two daughters are the same
    // and if the two events are the same or not (mixing)
    Bool_t sameCriteria = ((fCharge[0] == fCharge[1]) && (fDaughter[0] == fDaughter[1]));
    TString selList1  = "";
    TString selList2  = "";

    TObject *obj = fTrackCuts.FindObject("cutNEW");
    Int_t fCutID = fTrackCuts.IndexOf(obj);
    
    if (fMiniEvent->IsEmpty())	return;
    
    Int_t   n1 = fMiniEvent->CountParticles(fSel1, fCharge[0], fCutID);
    //Int_t   n1 = event1->CountParticles(fSel1, fCharge[0], fCutID);
    Int_t   n2 = fMiniEvent->CountParticles(fSel2, fCharge[1], fCutID);
    //Int_t   n2 = event2->CountParticles(fSel2, fCharge[1], fCutID);
    for (i1 = 0; i1 < n1; i1++) selList1.Append(Form("%d ", fSel1[i1]));
    for (i2 = 0; i2 < n2; i2++) selList2.Append(Form("%d ", fSel2[i2]));
    AliDebugClass(1, Form("[%10s] Part #1: [%s] -- evID %6d -- charge = %c -- cut ID = %d --> %4d tracks (%s)", GetName(), (fMiniEvent == fMiniEvent ? "def" : "mix"), fMiniEvent->ID(), fCharge[0], fCutID, n1, selList1.Data()));
    AliDebugClass(1, Form("[%10s] Part #2: [%s] -- evID %6d -- charge = %c -- cut ID = %d --> %4d tracks (%s)", GetName(), (fMiniEvent == fMiniEvent ? "def" : "mix"), fMiniEvent->ID(), fCharge[1], fCutID, n2, selList2.Data()));
    if (!n1 || !n2) {
      AliDebugClass(1, "No pairs to mix");
      return;
    }

    // external loop
    for (i1 = 0; i1 < n1; i1++) {
      p1 = fMiniEvent->GetParticle(fSel1[i1]);
      //p1 = event1->GetParticle(fSel1[i1]);
      // define starting point for inner loop
      // if daughter selection criteria (charge, cuts) are the same
      // and the two events coincide, internal loop must start from
      // the first track *after* current one;
      // otherwise it starts from the beginning
      start = (sameCriteria ? i1 + 1 : 0);
      AliDebugClass(2, Form("Start point = %d", start));
      // internal loop

      for (i2 = start; i2 < n2; i2++) {
	p2 = fMiniEvent->GetParticle(fSel2[i2]);
	//p2 = event2->GetParticle(fSel2[i2]);
	// avoid to mix a particle with itself
	if (p1->Index() == p2->Index()) {
	  AliDebugClass(2, "Skipping same index");
	  continue;
	}
	
	// sum momenta  
	mass1 = p1->StoredMass(kFALSE);
	if(!fUseStoredMass[0] || mass1 < 0.0) mass1 = GetMass(0);
	mass2 = p2->StoredMass(kFALSE);
	if(!fUseStoredMass[1] || mass2 < 0.0) mass2 = GetMass(1);
	
	fHMass1->Fill(mass1);
	fHMass2->Fill(mass2);
	dummyPair->Fill(p1, p2, mass1, mass2, fMotherMass);
	
	Double_t pairPt = dummyPair->Pt(kFALSE);
	Double_t pairPhi = dummyPair->PhiV(kFALSE);
	Double_t pairEta = dummyPair->Eta(kFALSE);
	Double_t pairInvMass = dummyPair->InvMass(kFALSE);

	fHMassPhi->Fill(pairInvMass);

	if (pairPt > 9.9 || pairPt < 0.0) continue;
	if (pairInvMass > 1.5 || pairInvMass < 0.0) continue;
	
	Int_t ptIndex = PtIndex(pairPt, 9.9, 0.0, 100.);
	Int_t invmassindex = InvMassIndex(pairInvMass, 1.998, 0.0, 1000.);
	
	fHTestPtMassCentrality->Fill(pairPt, pairInvMass, lPercentile);
            
	if (pairEta > dEtaGap/2) {
	  Double_t dPCosPos = TMath::Cos(iHarm * pairPhi); 
	  Double_t dPSinPos = TMath::Sin(iHarm * pairPhi);
	  FlowPVecPos[j][ptIndex][invmassindex] += TComplex(weightsP*dPCosPos,weightsP*dPSinPos,kFALSE);
	  MPPos[j][ptIndex][invmassindex] += TComplex(1.,0.,kFALSE).Re();
	}
            
	if (pairEta < -dEtaGap/2) {
	  Double_t dPCosNeg = TMath::Cos(iHarm * pairPhi); 
	  Double_t dPSinNeg = TMath::Sin(iHarm * pairPhi);
	  FlowPVecNeg[j][ptIndex][invmassindex] += TComplex(weightsN*dPCosNeg,weightsN*dPSinNeg,kFALSE);
	  MPNeg[j][ptIndex][invmassindex] += TComplex(1.,0.,kFALSE).Re();
	}      
      }
    }  
      
    TComplex tPosNeg = 999;
    TComplex tNegPos = 999;

    for (int p = 0; p < 100; ++p)
      {
	for (int m = 0; m < 1000; ++m)
	  {
            TComplex tPosNeg = FlowPVecPos[j][p][m]*TComplex::Conjugate(FlowVecQneg[j]);
            TComplex tNegPos = FlowPVecNeg[j][p][m]*TComplex::Conjugate(FlowVecQpos[j]);

            Double_t weightingPposQneg = MPPos[j][p][m]*MQNeg[j];
            Double_t weightingPnegQpos = MPNeg[j][p][m]*MQPos[j];

            if (weightingPposQneg == 0. || weightingPnegQpos == 0.)
	      {
		continue;
	      }
            
            if (TMath::Abs(tPosNeg.Re()/weightingPposQneg)  < 1)
	      {
		fHPposQneg[j][icentralityIndex]->Fill(dPtBins[p], dInvMassBins[m], tPosNeg.Re()/weightingPposQneg);
		fHPposQnegWeight[j][icentralityIndex]->Fill(dPtBins[p], dInvMassBins[m], tPosNeg.Re()/weightingPposQneg, weightingPposQneg);
	      }

            if (TMath::Abs(tNegPos.Re()/weightingPnegQpos)  < 1)
	      {
		fHPnegQpos[j][icentralityIndex]->Fill(dPtBins[p], dInvMassBins[m], tNegPos.Re()/weightingPnegQpos);
		fHPnegQposWeight[j][icentralityIndex]->Fill(dPtBins[p], dInvMassBins[m], tNegPos.Re()/weightingPnegQpos, weightingPnegQpos);
	      }
	  }
      }
  }  

  //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fill MC based histograms on mothers,
  // which do need the original event
  // if (fUseMC) {
  //   if (fRsnEvent.IsESD() && fMCEvent)
  //     FillTrueMotherESD(fMiniEvent);
  //   else if (fRsnEvent.IsAOD() && fRsnEvent.GetAODList())
  //     FillTrueMotherAOD(fMiniEvent);
  // }

  // if the event is not empty, store it
  if (fMiniEvent->IsEmpty()) {
    AliDebugClass(2, Form("Rejecting empty event #%d", fEvNum));
  } else {
    Int_t id = fEvBuffer->GetEntries();
    AliDebugClass(2, Form("Adding event #%d with ID = %d", fEvNum, id));
    fMiniEvent->ID() = id;
    fEvBuffer->Fill();
  }
  
  // post data for computed stuff
  PostData(1, fOutput);
}

//__________________________________________________________________________________________________


//__________________________________________________________________________________________________
void AliRsnMiniTaskPhiVn::FinishTaskOutput()
{
  //
  // This function is called at the end of the loop on available events,
  // and then the buffer will be full with all the corresponding mini-events,
  // each one containing all tracks selected by each of the available track cuts.
  // Here a loop is done on each of these events, and both single-event and mixing are computed
  //

  // security code: reassign the buffer to the mini-event cursor
  fEvBuffer->SetBranchAddress("events", &fMiniEvent);
  TStopwatch timer;
  // prepare variables
  Int_t ievt, nEvents = (Int_t)fEvBuffer->GetEntries();
  Int_t idef, nDefs   = fHistograms.GetEntries();
  Int_t imix, iloop, ifill;
  AliRsnMiniOutput *def = 0x0;
  AliRsnMiniOutput::EComputation compType;

  Int_t printNum = fMixPrintRefresh;
  if (printNum < 0) {
    if (nEvents>1e5) printNum=nEvents/100;
    else if (nEvents>1e4) printNum=nEvents/10;
    else printNum = 0;
  }

  // loop on events, and for each one fill all outputs
  // using the appropriate procedure depending on its type
  // only mother-related histograms are filled in UserExec,
  // since they require direct access to MC event
  timer.Start();
  for (ievt = 0; ievt < nEvents; ievt++) {
    // get next entry
    fEvBuffer->GetEntry(ievt);
    if (printNum&&(ievt%printNum==0)) {
      AliInfo(Form("[%s] Std.Event %d/%d",GetName(), ievt,nEvents));
      timer.Stop(); timer.Print(); fflush(stdout); timer.Start(kFALSE);
    }
    // fill
    for (idef = 0; idef < nDefs; idef++) {
      def = (AliRsnMiniOutput *)fHistograms[idef];
      if (!def) continue;
      compType = def->GetComputation();
      // execute computation in the appropriate way
      switch (compType) {
      case AliRsnMiniOutput::kEventOnly:
	//AliDebugClass(1, Form("Event %d, def '%s': event-value histogram filling", ievt, def->GetName()));
	ifill = 1;
	def->FillEvent(fMiniEvent, &fValues);
	break;
      case AliRsnMiniOutput::kTruePair:
	//AliDebugClass(1, Form("Event %d, def '%s': true-pair histogram filling", ievt, def->GetName()));
	ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
	break;
      case AliRsnMiniOutput::kTrackPair:
	//AliDebugClass(1, Form("Event %d, def '%s': pair-value histogram filling", ievt, def->GetName()));
	ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
	break;
      case AliRsnMiniOutput::kTrackPairRotated1:
	//AliDebugClass(1, Form("Event %d, def '%s': rotated (1) background histogram filling", ievt, def->GetName()));
	ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
	break;
      case AliRsnMiniOutput::kTrackPairRotated2:
	//AliDebugClass(1, Form("Event %d, def '%s': rotated (2) background histogram filling", ievt, def->GetName()));
	ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
	break;
      default:
	// other kinds are processed elsewhere
	ifill = 0;
	AliDebugClass(2, Form("Computation = %d", (Int_t)compType));
      }
      // message
      AliDebugClass(1, Form("Event %6d: def = '%15s' -- fills = %5d", ievt, def->GetName(), ifill));
    }
  }

  // if no mixing is required, stop here and post the output
  if (fNMix < 1) {
    AliDebugClass(2, "Stopping here, since no mixing is required");
    PostData(1, fOutput);
    return;
  }

  // initialize mixing counter
  Int_t    nmatched[nEvents];
  TString *smatched = new TString[nEvents];
  for (ievt = 0; ievt < nEvents; ievt++) {
    smatched[ievt] = "|";
    nmatched[ievt] = 0;
  }


  AliInfo(Form("[%s] Std.Event %d/%d",GetName(), nEvents,nEvents));
  timer.Stop(); timer.Print(); timer.Start(); fflush(stdout);

  // search for good matchings
  for (ievt = 0; ievt < nEvents; ievt++) {
    if (printNum&&(ievt%printNum==0)) {
      AliInfo(Form("[%s] EventMixing searching %d/%d",GetName(),ievt,nEvents));
      timer.Stop(); timer.Print(); timer.Start(kFALSE); fflush(stdout);
    }
    if (nmatched[ievt] >= fNMix) continue;
    fEvBuffer->GetEntry(ievt);
    AliRsnMiniEvent evMain(*fMiniEvent);
    for (iloop = 1; iloop < nEvents; iloop++) {
      imix = ievt + iloop;
      if (imix >= nEvents) imix -= nEvents;
      if (imix == ievt) continue;
      // text next entry
      fEvBuffer->GetEntry(imix);
      // skip if events are not matched
      if (!EventsMatch(&evMain, fMiniEvent)) continue;
      // check that the array of good matches for mixed does not already contain main event
      if (smatched[imix].Contains(Form("|%d|", ievt))) continue;
      // check that the found good events has not enough matches already
      if (nmatched[imix] >= fNMix) continue;
      // add new mixing candidate
      smatched[ievt].Append(Form("%d|", imix));
      nmatched[ievt]++;
      nmatched[imix]++;
      if (nmatched[ievt] >= fNMix) break;
    }
    AliDebugClass(1, Form("Matches for event %5d = %d [%s] (missing are declared above)", evMain.ID(), nmatched[ievt], smatched[ievt].Data()));
  }

  AliInfo(Form("[%s] EventMixing searching %d/%d",GetName(),nEvents,nEvents));
  timer.Stop(); timer.Print(); fflush(stdout); timer.Start();

  // perform mixing
  TObjArray *list = 0x0;
  TObjString *os = 0x0;
  for (ievt = 0; ievt < nEvents; ievt++) {
    if (printNum&&(ievt%printNum==0)) {
      AliInfo(Form("[%s] EventMixing %d/%d",GetName(),ievt,nEvents));
      timer.Stop(); timer.Print(); timer.Start(kFALSE); fflush(stdout);
    }
    ifill = 0;
    fEvBuffer->GetEntry(ievt);
    AliRsnMiniEvent evMain(*fMiniEvent);
    list = smatched[ievt].Tokenize("|");
    TObjArrayIter next(list);
    while ( (os = (TObjString *)next()) ) {
      imix = os->GetString().Atoi();
      fEvBuffer->GetEntry(imix);
      for (idef = 0; idef < nDefs; idef++) {
	def = (AliRsnMiniOutput *)fHistograms[idef];
	if (!def) continue;
	if (!def->IsTrackPairMix()) continue;
	ifill += def->FillPair(&evMain, fMiniEvent, &fValues, kTRUE);
	if (!def->IsSymmetric()) {
	  AliDebugClass(2, "Reflecting non symmetric pair");
	  ifill += def->FillPair(fMiniEvent, &evMain, &fValues, kFALSE);
	}
      }
    }
    delete list;
  }

  delete [] smatched;

  AliInfo(Form("[%s] EventMixing %d/%d",GetName(),nEvents,nEvents));
  timer.Stop(); timer.Print(); fflush(stdout);

  /*
    OLD
    ifill = 0;
    for (iloop = 1; iloop < nEvents; iloop++) {
    imix = ievt + iloop;
    // restart from beginning if reached last event
    if (imix >= nEvents) imix -= nEvents;
    // avoid to mix an event with itself
    if (imix == ievt) continue;
    // skip all events already mixed enough times
    if (fNMixed[ievt] >= fNMix) break;
    if (fNMixed[imix] >= fNMix) continue;
    fEvBuffer->GetEntry(imix);
    // skip if events are not matched
    if (TMath::Abs(evMain.Vz()    - fMiniEvent->Vz()   ) > fMaxDiffVz   ) continue;
    if (TMath::Abs(evMain.Mult()  - fMiniEvent->Mult() ) > fMaxDiffMult ) continue;
    if (TMath::Abs(evMain.Angle() - fMiniEvent->Angle()) > fMaxDiffAngle) continue;
    // found a match: increment counter for both events
    AliDebugClass(1, Form("Event %d, def '%s': event mixing (%d with %d)", ievt, def->GetName(), ievt, imix));
    fNMixed[ievt]++;
    fNMixed[imix]++;
    // process mixing
    ifill += def->FillPair(&evMain, fMiniEvent, &fValues);
    // stop if mixed enough times
    if (fNMixed[ievt] >= fNMix) break;
    }
    break;
    // print number of mixings done with each event
    for (ievt = 0; ievt < nEvents; ievt++) {
    AliDebugClass(2, Form("Event %6d: mixed %2d times", ievt, fNMixed[ievt]));
    }
  */

  // post computed data
  PostData(1, fOutput);
}

//__________________________________________________________________________________________________
void AliRsnMiniTaskPhiVn::Terminate(Option_t *)
{
  //
  // Draw result to screen, or perform fitting, normalizations
  // Called once at the end of the query
  //

  fOutput = dynamic_cast<TList *>(GetOutputData(1));
  if (!fOutput) {
    AliError("Could not retrieve TList fOutput");
    return;
  }
}

//__________________________________________________________________________________________________
Char_t AliRsnMiniTaskPhiVn::CheckCurrentEvent()
{
  //
  // This method checks if current event is OK for analysis.
  // In case it is, the pointers of the local AliRsnEvent data member
  // will point to it, in order to allow cut checking, otherwise the
  // function exits with a failure message.
  // ---
  // ESD events must pass the physics selection, AOD are supposed to do.
  // ---
  // While checking the event, a histogram is filled to count the number
  // of CINT1B, V0AND and CANDLE events, which are needed for normalization
  // ---
  // Return values can be:
  //    -- 'E' if the event is accepted and is ESD
  //    -- 'A' if the event is accepted and is AOD
  //    --  0  if the event is not accepted
  //

  // string to sum messages
  TString msg("");

  // check input type
  // exit points are provided in all cases an event is bad
  // if this block is passed, an event can be rejected only
  // if it does not pass the set of event cuts defined in the task
  Char_t output = 0;
  Bool_t isSelected;
  if (fInputEvent->InheritsFrom(AliESDEvent::Class())) {
    // type ESD
    output = 'E';
    // ESD specific check: Physics Selection
    // --> if this is failed, the event is rejected
    isSelected = (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);

    if (!isSelected) {
      AliDebugClass(2, "Event does not pass physics selections");
      fRsnEvent.SetRef(0x0);
      fRsnEvent.SetRefMC(0x0);
      return 0;
    }
    // set reference to input
    fRsnEvent.SetRef(fInputEvent);
    // add MC if requested and available
    if (fUseMC) {
      if (fMCEvent)
	fRsnEvent.SetRefMC(fMCEvent);
      else {
	AliWarning("MC event requested but not available");
	fRsnEvent.SetRefMC(0x0);
      }
    }
  } else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
    // type AOD
    output = 'A';
    // set reference to input
    fRsnEvent.SetRef(fInputEvent);
    // add MC if requested and available (it is in the same object)
    if (fUseMC) {
      fRsnEvent.SetRefMC(fInputEvent);
      if (!fRsnEvent.GetAODList()) {
	AliWarning("MC event requested but not available");
	fRsnEvent.SetRefMC(0x0);
      }
    }
  } else {
    AliError(Form("Bad input event class: %s", fInputEvent->ClassName()));
    // reset pointers in local AliRsnEvent object
    fRsnEvent.SetRef(0x0);
    fRsnEvent.SetRefMC(0x0);
    return 0;
  }

  // fill counter of accepted events
  fHEventStat->Fill(0.1);

  // check if it is V0AND
  // --> uses a cast to AliESDEvent even if the input is an AliAODEvent
  Bool_t v0A = fTriggerAna->IsOfflineTriggerFired((AliESDEvent *)fInputEvent, AliTriggerAnalysis::kV0A);
  Bool_t v0C = fTriggerAna->IsOfflineTriggerFired((AliESDEvent *)fInputEvent, AliTriggerAnalysis::kV0C);
  if (v0A && v0C) {
    msg += " -- VOAND = YES";
    fHEventStat->Fill(1.1);
  } else {
    msg += " -- VOAND = NO ";
  }

  // check candle
  // --> requires at least one good quality track with Pt > 0.5 and |eta| <= 0.8
  Int_t iTrack, ntracksLoop = fInputEvent->GetNumberOfTracks();
  Bool_t candle = kFALSE;
  for (iTrack = 0; iTrack < ntracksLoop; iTrack++) {
    AliVTrack   *track = (AliVTrack *)fInputEvent->GetTrack(iTrack);
    AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
    AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
    if (track->Pt() < 0.5) continue;
    if(TMath::Abs(track->Eta()) > 0.8) continue;
    if (esdt) if (!fESDtrackCuts->AcceptTrack(esdt)) continue;
    if (aodt) if (!aodt->TestFilterBit(5)) continue;
    candle = kTRUE;
    break;
  }
  if (candle) {
    msg += " -- CANDLE = YES";
    fHEventStat->Fill(2.1);
  } else {
    msg += " -- CANDLE = NO ";
  }

  // if event cuts are defined, they are checked here
  // final decision on the event depends on this
  isSelected = kTRUE;
  if (fEventCuts) {
    if (!fEventCuts->IsSelected(&fRsnEvent)) {
      msg += " -- Local cuts = REJECTED";
      isSelected = kFALSE;
    } else {
      msg += " -- Local cuts = ACCEPTED";
      isSelected = kTRUE;
    }
  } else {
    msg += " -- Local cuts = NONE";
    isSelected = kTRUE;
  }

  // if the above exit point is not taken, the event is accepted
  AliDebugClass(2, Form("Stats: %s", msg.Data()));
  if (isSelected) {
    fHEventStat->Fill(3.1);
    //Double_t multi = ComputeCentrality((output == 'E'));
    Double_t multi = ComputeMultiplicity(kFALSE, "ALIMULTSELECTION_V0M");
    Double_t tracklets = ComputeTracklets();
    fHAEventsVsMulti->Fill(multi);
    fHAEventsVsTracklets->Fill(tracklets);
    if(fHAEventVz) fHAEventVz->Fill(multi,fInputEvent->GetPrimaryVertex()->GetZ());
    if(fHAEventMultiCent) fHAEventMultiCent->Fill(multi,ComputeMultiplicity(output == 'E',fHAEventMultiCent->GetYaxis()->GetTitle()));
    if(fHAEventPlane) fHAEventPlane->Fill(multi,ComputeAngle());
    return output;
  } else {
    fHEventStat->Fill(4.1);
    const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
    if(!vertex) fHEventStat->Fill(5.1);
    else{
      TString title=vertex->GetTitle();
      if( (title.Contains("Z")) || (title.Contains("3D")) ) fHEventStat->Fill(5.1);
      if(vertex->GetNContributors()<1.) fHEventStat->Fill(6.1);
      if(TMath::Abs(vertex->GetZ())>10.) fHEventStat->Fill(7.1);
    }
    return 0;
  }
}

//__________________________________________________________________________________________________
void AliRsnMiniTaskPhiVn::FillMiniEvent(Char_t evType)
{
  //
  // Refresh cursor mini-event data member to fill with current event.
  // Returns the total number of tracks selected.
  //

  // assign event-related values
  //if (fMiniEvent) delete fMiniEvent;
  //fMiniEvent = new AliRsnMiniEvent;
  fMiniEvent->Clear();
  fMiniEvent->Mult()  = ComputeMultiplicity(kFALSE, "ALIMULTSELECTION_V0M");  
  fMiniEvent->Tracklets() = ComputeTracklets();
  fMiniEvent->SetRef(fRsnEvent.GetRef());
  fMiniEvent->SetRefMC(fRsnEvent.GetRefMC());
  fMiniEvent->Vz() = fInputEvent->GetPrimaryVertex()->GetZ();
  fMiniEvent->Angle() = ComputeAngle();
  AliDebugClass(2, Form("Event %d: type = %c -- vz = %f -- mult = %f -- angle = %f", fEvNum, evType, fMiniEvent->Vz(), fMiniEvent->Mult(), fMiniEvent->Angle()));
  
  // loop on daughters and assign track-related values
  Int_t ic, ncuts = fTrackCuts.GetEntries();
  Int_t ip, npart = fRsnEvent.GetAbsoluteSum();
  Int_t npos = 0, nneg = 0, nneu = 0;
  AliRsnDaughter cursor;
  AliRsnMiniParticle miniParticle;
  for (ip = 0; ip < npart; ip++) {
    // point cursor to next particle
    fRsnEvent.SetDaughter(cursor, ip);
    // copy momentum and MC info if present
    miniParticle.CopyDaughter(&cursor);
    miniParticle.Index() = ip;
    // switch on the bits corresponding to passed cuts
    for (ic = 0; ic < ncuts; ic++) {
      AliRsnCutSet *cuts = (AliRsnCutSet *)fTrackCuts[ic];
      if (cuts->IsSelected(&cursor)) miniParticle.SetCutBit(ic);
    }
    // if a track passes at least one track cut, it is added to the pool
    if (miniParticle.CutBits()) {
      fMiniEvent->AddParticle(miniParticle);
      if (miniParticle.Charge() == '+') npos++;
      else if (miniParticle.Charge() == '-') nneg++;
      else nneu++;
    }
  }

  // get number of accepted tracks
  AliDebugClass(1, Form("Event %6d: total = %5d, accepted = %4d (pos %4d, neg %4d, neu %4d)", fEvNum, npart, (Int_t)fMiniEvent->Particles().GetEntriesFast(), npos, nneg, nneu));
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniTaskPhiVn::ComputeAngle()
{
  //
  // Get the plane angle
  //

  AliEventplane *plane = 0x0;

  if (fInputEvent->InheritsFrom(AliESDEvent::Class()))
    plane = fInputEvent->GetEventplane();
  else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
    AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
    plane = ((AliVAODHeader*)aodEvent->GetHeader())->GetEventplaneP();
  }

  if (plane)
    return plane->GetEventplane("Q");
  else {
    AliWarning("No event plane defined");
    return 1E20;
  }
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniTaskPhiVn::ComputeCentrality(Bool_t isESD)
{
  //
  // Computes event centrality/multiplicity according to the criterion defined
  // by two elements: (1) choice between multiplicity and centrality and
  // (2) the string defining what criterion must be used for specific computation.
  //


  if (fUseCentrality) {
      
    if (fCentralityType.CompareTo("ALIMULTSELECTION_V0M"))
      {
	Float_t lPercentile = 300; 
	AliMultSelection *MultSelection = 0x0; 
	MultSelection = (AliMultSelection * ) fInputEvent->FindListObject("MultSelection");
      
	if( !MultSelection) {
	  //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
	  AliWarning("AliMultSelection object not found!");
	}else{
	  //This is the Centrality
	  lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
	}
	return lPercentile;
      }
      

    if ((!fUseMC) && (fUseCentralityPatchPbPb2011)) {
      return ApplyCentralityPatchPbPb2011();//
    }
    if ((!fUseMC) && (!isESD) && (fUseAOD049CentralityPatch)) {
      return ApplyCentralityPatchAOD049();
    } else {
      AliCentrality *centrality = fInputEvent->GetCentrality();
      if (!centrality) {
	AliError("Cannot compute centrality!");
	return -1.0;
      }
      return centrality->GetCentralityPercentile(fCentralityType.Data());
    }
  } else {
    if (!fCentralityType.CompareTo("TRACKS"))
      return fInputEvent->GetNumberOfTracks();
    else if (!fCentralityType.CompareTo("QUALITY"))
      if (isESD)
	return AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent *)fInputEvent, kTRUE);
      else {
	Double_t count = 0.;
	Int_t iTrack, ntracksLoop = fInputEvent->GetNumberOfTracks();
	for (iTrack = 0; iTrack < ntracksLoop; iTrack++) {
	  AliVTrack   *track = (AliVTrack *)fInputEvent->GetTrack(iTrack);
	  AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
	  if (!aodt) continue;
	  if (!aodt->TestFilterBit(5)) continue;
	  count++;
	}
	return count;
      }
    else if (!fCentralityType.CompareTo("TRACKLETS")) {
      if (isESD) {
	const AliMultiplicity *mult = ((AliESDEvent *)fInputEvent)->GetMultiplicity();
	Float_t nClusters[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
	for(Int_t ilay = 0; ilay < 6; ilay++) nClusters[ilay] = (Float_t)mult->GetNumberOfITSClusters(ilay);
	return AliESDUtils::GetCorrSPD2(nClusters[1], fInputEvent->GetPrimaryVertex()->GetZ());
      } else {
	AliWarning("Cannot compute multiplicity with SPD tracklets from AOD");
	return 1E20;
      }
    } else {
      AliError(Form("String '%s' does not define a possible multiplicity/centrality computation", fCentralityType.Data()));
      return -1.0;
    }
  }
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniTaskPhiVn::ComputeMultiplicity(Bool_t isESD,TString type)
{
  //
  // Computes event multiplicity according to the string defining
  // what criterion must be used for specific computation.
  //

  type.ToUpper();


  //if (type.CompareTo("ALIMULTSELECTION_V0M"))
  if (type.CompareTo("AliMultSelection_V0M"))
    {
      Float_t lPercentile = 300; 
      AliMultSelection *MultSelection = 0x0; 
      MultSelection = (AliMultSelection * ) fInputEvent->FindListObject("MultSelection");
      
      if( !MultSelection) {
	//If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
	AliWarning("AliMultSelection object not found!");
      }else{
	//This is the Centrality
	lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
      }
      return lPercentile;
    }

  if (!type.CompareTo("TRACKS"))
    return fInputEvent->GetNumberOfTracks();
  else if (!type.CompareTo("QUALITY"))
    if (isESD)
      return AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent *)fInputEvent, kTRUE);
    else {
      Double_t count = 0.;
      Int_t iTrack, ntracksLoop = fInputEvent->GetNumberOfTracks();
      for (iTrack = 0; iTrack < ntracksLoop; iTrack++) {
	AliVTrack   *track = (AliVTrack *)fInputEvent->GetTrack(iTrack);
	AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
	if (!aodt) continue;
	if (!aodt->TestFilterBit(5)) continue;
	count++;
      }
      return count;
    }
  else if (!type.CompareTo("TRACKLETS")) {
    if (isESD) {
      const AliMultiplicity *mult = ((AliESDEvent *)fInputEvent)->GetMultiplicity();
      Float_t nClusters[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
      for(Int_t ilay = 0; ilay < 6; ilay++) nClusters[ilay] = (Float_t)mult->GetNumberOfITSClusters(ilay);
      return AliESDUtils::GetCorrSPD2(nClusters[1], fInputEvent->GetPrimaryVertex()->GetZ());
    } else {
      AliWarning("Cannot compute multiplicity with SPD tracklets from AOD");
      return 1E20;
    }
  } else {
    AliError(Form("String '%s' does not define a possible multiplicity/centrality computation", type.Data()));
    return -1.0;
  }
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniTaskPhiVn::ComputeTracklets()
{
  //
  // Get number of tracklets
  //

  Double_t count = 100;

  if (fInputEvent->InheritsFrom(AliESDEvent::Class())){
    AliESDEvent *esdEvent = (AliESDEvent *)fInputEvent;
    const AliMultiplicity *spdmult = esdEvent->GetMultiplicity();
    count = 1.0*spdmult->GetNumberOfTracklets();
  }
  else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
    AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
    AliAODTracklets *spdmult = aodEvent->GetTracklets();
    count = 1.0*spdmult->GetNumberOfTracklets();
  }

  return count;
}

//__________________________________________________________________________________________________
void AliRsnMiniTaskPhiVn::FillTrueMotherESD(AliRsnMiniEvent *miniEvent)
{
  //
  // Fills the histograms with true mother (ESD version)
  //

  Bool_t okMatch;
  Int_t id, ndef = fHistograms.GetEntries();
  Int_t ip, label1, label2, npart = fMCEvent->GetNumberOfTracks();
  static AliRsnMiniPair miniPair;
  AliMCParticle *daughter1, *daughter2;
  TLorentzVector p1, p2;
  AliRsnMiniOutput *def = 0x0;

  for (id = 0; id < ndef; id++) {
    def = (AliRsnMiniOutput *)fHistograms[id];
    if (!def) continue;
    if (!def->IsMother() && !def->IsMotherInAcc()) continue;
    for (ip = 0; ip < npart; ip++) {
      AliMCParticle *part = (AliMCParticle *)fMCEvent->GetTrack(ip);
      //get mother pdg code
      if (part->Particle()->GetPdgCode() != def->GetMotherPDG()) continue;
      // check that daughters match expected species
      if (part->Particle()->GetNDaughters() < 2) continue; 
      if (fMaxNDaughters > 0 && part->Particle()->GetNDaughters() > fMaxNDaughters) continue;
      //label1 = part->Particle()->GetDaughter(0); // Before change in accessing MC infor in AliRoot v5-09-46
      //label2 = part->Particle()->GetDaughter(1); // Before change in accessing MC infor in AliRoot v5-09-46
      label1 = part->GetDaughterLabel(0);
      label2 = part->GetDaughterLabel(1);
      daughter1 = (AliMCParticle *)fMCEvent->GetTrack(label1);
      daughter2 = (AliMCParticle *)fMCEvent->GetTrack(label2);
      okMatch = kFALSE;
      if (TMath::Abs(daughter1->Particle()->GetPdgCode()) == def->GetPDG(0) && TMath::Abs(daughter2->Particle()->GetPdgCode()) == def->GetPDG(1)) {
	okMatch = kTRUE;
	p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(0));
	p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(1));
      } else if (TMath::Abs(daughter1->Particle()->GetPdgCode()) == def->GetPDG(1) && TMath::Abs(daughter2->Particle()->GetPdgCode()) == def->GetPDG(0)) {
	okMatch = kTRUE;
	p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(1));
	p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(0));
      }
      if (!okMatch) continue;
      if(fCheckP && (TMath::Abs(part->Px()-(daughter1->Px()+daughter2->Px()))/(TMath::Abs(part->Px())+1.e-13)) > 0.00001 &&   
	 (TMath::Abs(part->Py()-(daughter1->Py()+daughter2->Py()))/(TMath::Abs(part->Py())+1.e-13)) > 0.00001 &&
	 (TMath::Abs(part->Pz()-(daughter1->Pz()+daughter2->Pz()))/(TMath::Abs(part->Pz())+1.e-13)) > 0.00001 ) continue;
      if(fCheckFeedDown){
	Int_t pdgGranma = 0;
	Int_t mother = 0;
	mother = part->GetMother();
	Int_t istep = 0;
	Int_t abspdgGranma =0;
	Bool_t isFromB=kFALSE;
	Bool_t isQuarkFound=kFALSE;
	while (mother >=0 ){
	  istep++;
	  AliDebug(2,Form("mother at step %d = %d", istep, mother));
	  AliMCParticle* mcGranma = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(mother));
	  if (mcGranma){
            pdgGranma = mcGranma->PdgCode();
            AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
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
	if(fRejectIfNoQuark && !isQuarkFound) pdgGranma = -99999;
	if(isFromB){
	  if (!fKeepDfromB) pdgGranma = -9999; //skip particle if come from a B meson.
	}
	else{ 
	  if (fKeepDfromBOnly) pdgGranma = -999;
        
	  if (pdgGranma == -99999){
	    AliDebug(2,"This particle does not have a quark in his genealogy\n");
	    continue;
	  }
	  if (pdgGranma == -9999){
	    AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only the prompt charm particles\n");   
	    continue;
	  }  
      
	  if (pdgGranma == -999){
	    AliDebug(2,"This particle come from a prompt charm particles but according to the settings of the task, we want only the ones coming from B\n");  
	    continue;
	  }  

	}
      }
      // assign momenta to computation object
      miniPair.Sum(0) = miniPair.Sum(1) = (p1 + p2);
      miniPair.FillRef(def->GetMotherMass());
      // do computations
      def->FillMother(&miniPair, miniEvent, &fValues);
      if(fKeepMotherInAcceptance){
	if(daughter1->Pt()<fMotherAcceptanceCutMinPt || daughter2->Pt()<fMotherAcceptanceCutMinPt || TMath::Abs(daughter1->Eta())>fMotherAcceptanceCutMaxEta ||  TMath::Abs(daughter2->Eta())>fMotherAcceptanceCutMaxEta) continue;
	def->FillMotherInAcceptance(&miniPair, miniEvent, &fValues);
      }  
    
    
    }
  }
}

//__________________________________________________________________________________________________
void AliRsnMiniTaskPhiVn::FillTrueMotherAOD(AliRsnMiniEvent *miniEvent)
{
  //
  // Fills the histograms with true mother (AOD version)
  //

  Bool_t okMatch;
  TClonesArray *list = fRsnEvent.GetAODList();
  Int_t id, ndef = fHistograms.GetEntries();
  Int_t ip, label1, label2, npart = list->GetEntries();
  static AliRsnMiniPair miniPair;
  AliAODMCParticle *daughter1, *daughter2;
  TLorentzVector p1, p2;
  AliRsnMiniOutput *def = 0x0;

  for (id = 0; id < ndef; id++) {
    def = (AliRsnMiniOutput *)fHistograms[id];
    if (!def) continue;
    if (!def->IsMother() && !def->IsMotherInAcc()) continue;
    for (ip = 0; ip < npart; ip++) {
      AliAODMCParticle *part = (AliAODMCParticle *)list->At(ip);
      if (part->GetPdgCode() != def->GetMotherPDG()) continue;
      // check that daughters match expected species
      if (part->GetNDaughters() < 2) continue;
      if (fMaxNDaughters > 0 && part->GetNDaughters() > fMaxNDaughters) continue;
      label1 = part->GetDaughterLabel(0);
      label2 = part->GetDaughterLabel(1);
      daughter1 = (AliAODMCParticle *)list->At(label1);
      daughter2 = (AliAODMCParticle *)list->At(label2);
      okMatch = kFALSE;
      if (TMath::Abs(daughter1->GetPdgCode()) == def->GetPDG(0) && TMath::Abs(daughter2->GetPdgCode()) == def->GetPDG(1)) {
	okMatch = kTRUE;
	p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(0));
	p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(1));
      } else if (TMath::Abs(daughter1->GetPdgCode()) == def->GetPDG(1) && TMath::Abs(daughter2->GetPdgCode()) == def->GetPDG(0)) {
	okMatch = kTRUE;
	p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(1));
	p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(0));
      }
      if (!okMatch) continue;
      if(fCheckP && (TMath::Abs(part->Px()-(daughter1->Px()+daughter2->Px()))/(TMath::Abs(part->Px())+1.e-13)) > 0.00001 &&   
	 (TMath::Abs(part->Py()-(daughter1->Py()+daughter2->Py()))/(TMath::Abs(part->Py())+1.e-13)) > 0.00001 &&
	 (TMath::Abs(part->Pz()-(daughter1->Pz()+daughter2->Pz()))/(TMath::Abs(part->Pz())+1.e-13)) > 0.00001 ) continue;
      if(fCheckFeedDown){
	Int_t pdgGranma = 0;
	Int_t mother = 0;
	mother = part->GetMother();
	Int_t istep = 0;
	Int_t abspdgGranma =0;
	Bool_t isFromB=kFALSE;
	Bool_t isQuarkFound=kFALSE;
	while (mother >=0 ){
	  istep++;
	  AliDebug(2,Form("mother at step %d = %d", istep, mother));
	  AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(list->At(mother));
	  if (mcGranma){
            pdgGranma = mcGranma->GetPdgCode();
            AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
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
	if(fRejectIfNoQuark && !isQuarkFound) pdgGranma = -99999;
	if(isFromB){
	  if (!fKeepDfromB) pdgGranma = -9999; //skip particle if come from a B meson.
	}
	else{ 
	  if (fKeepDfromBOnly) pdgGranma = -999;
        }
        
	if (pdgGranma == -99999){
	  AliDebug(2,"This particle does not have a quark in his genealogy\n");
	  continue;
	}
	if (pdgGranma == -9999){
	  AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only the prompt charm particles\n");   
	  continue;
	}  
      
	if (pdgGranma == -999){
	  AliDebug(2,"This particle come from a prompt charm particles but according to the settings of the task, we want only the ones coming from B\n");  
	  continue;
	}  
      } 
      // assign momenta to computation object
      miniPair.Sum(0) = miniPair.Sum(1) = (p1 + p2);
      miniPair.FillRef(def->GetMotherMass());
      // do computations
      def->FillMother(&miniPair, miniEvent, &fValues);
      if(fKeepMotherInAcceptance){
	if(daughter1->Pt()<fMotherAcceptanceCutMinPt || daughter2->Pt()<fMotherAcceptanceCutMinPt || TMath::Abs(daughter1->Eta())>fMotherAcceptanceCutMaxEta ||  TMath::Abs(daughter2->Eta())>fMotherAcceptanceCutMaxEta) continue;
	def->FillMotherInAcceptance(&miniPair, miniEvent, &fValues);
      }
    }
  }
}
//___________________________________________________________
void AliRsnMiniTaskPhiVn::SetDselection(UShort_t originDselection)
{
  // setting the way the D0 will be selected
  // 0 --> only from c quarks
  // 1 --> only from b quarks
  // 2 --> from both c quarks and b quarks
      
  fOriginDselection = originDselection;
   
  if (fOriginDselection == 0) {
    fKeepDfromB = kFALSE;
    fKeepDfromBOnly = kFALSE;
  }
   
  if (fOriginDselection == 1) {
    fKeepDfromB = kTRUE;
    fKeepDfromBOnly = kTRUE;
  }
   
  if (fOriginDselection == 2) {
    fKeepDfromB = kTRUE;
    fKeepDfromBOnly = kFALSE;
  }
   
  return;  
}
//__________________________________________________________________________________________________
Bool_t AliRsnMiniTaskPhiVn::EventsMatch(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2)
{
  //
  // Check if two events are compatible.
  // If the mixing is continuous, this is true if differences in vz, mult and angle are smaller than
  // the specified values.
  // If the mixing is binned, this is true if the events are in the same bin.
  //

  if (!event1 || !event2) return kFALSE;
  Int_t ivz1, ivz2, imult1, imult2, iangle1, iangle2;
  Double_t dv, dm, da;

  if (fContinuousMix) {
    dv = TMath::Abs(event1->Vz()    - event2->Vz()   );
    dm = TMath::Abs(event1->Mult()  - event2->Mult() );
    da = TMath::Abs(event1->Angle() - event2->Angle());
    if (dv > fMaxDiffVz) {
      //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Vz = %f", event1->ID(), event2->ID(), dv));
      return kFALSE;
    }
    if (dm > fMaxDiffMult ) {
      //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Mult = %f", event1->ID(), event2->ID(), dm));
      return kFALSE;
    }
    if (da > fMaxDiffAngle) {
      //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Angle = %f", event1->ID(), event2->ID(), da));
      return kFALSE;
    }
    return kTRUE;
  } else {
    ivz1 = (Int_t)(event1->Vz() / fMaxDiffVz);
    ivz2 = (Int_t)(event2->Vz() / fMaxDiffVz);
    imult1 = (Int_t)(event1->Mult() / fMaxDiffMult);
    imult2 = (Int_t)(event2->Mult() / fMaxDiffMult);
    iangle1 = (Int_t)(event1->Angle() / fMaxDiffAngle);
    iangle2 = (Int_t)(event2->Angle() / fMaxDiffAngle);
    if (ivz1 != ivz2) return kFALSE;
    if (imult1 != imult2) return kFALSE;
    if (iangle1 != iangle2) return kFALSE;
    return kTRUE;
  }
}

//---------------------------------------------------------------------
Double_t AliRsnMiniTaskPhiVn::ApplyCentralityPatchPbPb2011(){
  //This part rejects randomly events such that the centrality gets flat for LHC11h Pb-Pb data
  //for 0-5% and 10-20% centrality bin
  
  if (fCentralityType!="V0M") {
    AliWarning("Wrong value (not centrality from V0).");
    return -999.0;
  }
  
  AliCentrality *centrality = fInputEvent->GetCentrality();
  if (!centrality) {
    AliWarning("Cannot get centrality from AOD event.");
    return -999.0;
  }
  
  Double_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));               
  Double_t rnd_hc = -1., testf = 0.0, ff = 0, N1 = -1., N2 = -1.;

  if(fUseCentralityPatchPbPb2011==510){
    N1 = 1.9404e+06;
    N2 = 1.56435e+06; //N2 is the reference 
    ff = 5.04167e+06 - 1.49885e+06*cent + 2.35998e+05*cent*cent -1.22873e+04*cent*cent*cent;
  } else {
    if(fUseCentralityPatchPbPb2011==1020){
      N2 = 2.0e+05; //N2 is the reference
      N1 = 3.7e+05;
      ff = -1.73979e+06 - 3.05316e+06*cent + 1.05517e+06*cent*cent - 133205*cent*cent*cent + 8187.45*cent*cent*cent*cent - 247.875*cent*cent*cent*cent*cent + 2.9676*cent*cent*cent*cent*cent*cent;
    } else {
      AliError(Form("Patch for the requested centrality (%i) is not available", fUseCentralityPatchPbPb2011));
      return -999.0;
    }
  }
  testf = ( N2 + (N1-ff) ) / N1;
  rnd_hc = gRandom->Rndm();

  //AliDebugClass(1, Form("Flat Centrality %d", fUseCentralityPatchPbPb2011));

  if (rnd_hc < 0 || rnd_hc > 1 ) 
    {
      AliWarning("Wrong Random number generated");
      return -999.0;
    }
  
  if (rnd_hc < testf)
    return cent;
  else
    return -999.0;
}
//---------------------------------------------------------------------
Double_t AliRsnMiniTaskPhiVn::ApplyCentralityPatchAOD049()
{
  //
  //Apply centrality patch for AOD049 outliers
  //
  if (fInputEvent->InheritsFrom(AliESDEvent::Class())) {
    AliWarning("Requested patch for AOD049 for ESD. ");
    return -999.0;
  }

  if (fCentralityType!="V0M") {
    AliWarning("Requested patch forAOD049 for wrong value (not centrality from V0).");
    return -999.0;
  }

  AliCentrality *centrality = fInputEvent->GetCentrality();
  if (!centrality) {
    AliWarning("Cannot get centrality from AOD event.");
    return -999.0;
  }

  Float_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));
  /*
    Bool_t isSelRun = kFALSE;
    Int_t selRun[5] = {138364, 138826, 138828, 138836, 138871};
    if(cent<0){
    Int_t quality = centrality->GetQuality();
    if(quality<=1){
    cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
    } else {
    Int_t runnum=aodEvent->GetRunNumber();
    for(Int_t ir=0;ir<5;ir++){
    if(runnum==selRun[ir]){
    isSelRun=kTRUE;
    break;
    }
    }
    if((quality==8||quality==9)&&isSelRun) cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
    }
    }
  */

  if(cent>=0.0) {
    Float_t v0 = 0.0;
    AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
    AliAODVZERO *aodV0 = (AliAODVZERO *) aodEvent->GetVZEROData();
    v0+=aodV0->GetMTotV0A();
    v0+=aodV0->GetMTotV0C();
    if ( (cent==0) && (v0<19500) ) {
      AliDebug(3, Form("Filtering issue in centrality -> cent = %5.2f",cent));
      return -999.0;
    }
    Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
    Float_t val = 1.30552 +  0.147931 * v0;

    Float_t tklSigma[101] = {176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86,
			     120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654,
			     92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334,
			     68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224,
			     51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255,
			     37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398,
			     26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235,
			     19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504,
			     12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544,
			     13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544
    };

    if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)cent] )  {
      AliDebug(3, Form("Outlier event in centrality -> cent = %5.2f",cent));
      return -999.0;
    }
  } else {
    //force it to be -999. whatever the negative value was
    cent = -999.;
  }
  return cent;
}

//----------------------------------------------------------------------------------
void AliRsnMiniTaskPhiVn::SetEventQAHist(TString type,TH2F *histo)
{
  if(!histo) {
    AliWarning(Form("event QA histogram pointer not defined for slot %s",type.Data()));
    return;
  }

  type.ToLower();

  if(!type.CompareTo("vz")) fHAEventVz = histo;
  else if(!type.CompareTo("multicent")) {
    TString mtype(histo->GetYaxis()->GetTitle());
    mtype.ToUpper();
    if(mtype.CompareTo("QUALITY") && mtype.CompareTo("TRACKS") && mtype.CompareTo("TRACKLETS")) {
      AliWarning(Form("multiplicity vs. centrality histogram y-axis %s unknown, setting to TRACKS",mtype.Data()));
      histo->GetYaxis()->SetTitle("TRACKS");
    }
    fHAEventMultiCent = histo;
  }
  else if(!type.CompareTo("eventplane")) fHAEventPlane = histo;
  else AliWarning(Form("event QA histogram slot %s undefined",type.Data()));

  return;
}

//----------------------------------------------------------------------------------
Int_t AliRsnMiniTaskPhiVn::CreateValue(AliRsnMiniValue::EType type, Bool_t useMC)
{
  //
  // Create a new value in the task,
  // and returns its ID, which is needed for setting up histograms.
  // If that value was already initialized, returns its ID and does not recreate it.
  //

  Int_t valID = ValueID(type, useMC);
  if (valID >= 0 && valID < fValues.GetEntries()) {
    AliInfo(Form("Value '%s' is already created in slot #%d", AliRsnMiniValue::ValueName(type, useMC), valID));
  } else {
    valID = fValues.GetEntries();
    AliInfo(Form("Creating value '%s' in slot #%d", AliRsnMiniValue::ValueName(type, useMC), valID));
    new (fValues[valID]) AliRsnMiniValue(type, useMC);
  }

  return valID;
}

//----------------------------------------------------------------------------------
Int_t AliRsnMiniTaskPhiVn::ValueID(AliRsnMiniValue::EType type, Bool_t useMC)
{
  //
  // Searches if a value computation is initialized
  //

  const char *name = AliRsnMiniValue::ValueName(type, useMC);
  TObject *obj = fValues.FindObject(name);
  if (obj)
    return fValues.IndexOf(obj);
  else
    return -1;
}
//__________________________________________________________________________________
Int_t AliRsnMiniTaskPhiVn::PtIndex(Double_t pt, Double_t ptmax, Double_t ptmin, Double_t totalBinsPt){
  Double_t widthpt = (ptmax - ptmin)/(totalBinsPt);
  Int_t iIndex = int(pt/widthpt);
  return iIndex;
}
//__________________________________________________________________________________
Int_t AliRsnMiniTaskPhiVn::InvMassIndex(Double_t invmass, Double_t invmassmax, Double_t invmassmin, Double_t totalBinsInvMass){
  Double_t widthinvmass = (invmassmax - invmassmin)/(totalBinsInvMass);
  Int_t iIndex = int(invmass/widthinvmass);
  return iIndex;
}
//__________________________________________________________________________________
Int_t AliRsnMiniTaskPhiVn::CentralityIndex(Double_t centrality, Double_t centmax, Double_t centmin, Double_t totalBinsCentrality){
  Double_t widthcent = (centmax - centmin)/(totalBinsCentrality);
  Int_t iIndex = int(centrality/widthcent);
  return iIndex;
}
//__________________________________________________________________________________
AliRsnMiniOutput *AliRsnMiniTaskPhiVn::CreateOutput(const char *name, AliRsnMiniOutput::EOutputType type, AliRsnMiniOutput::EComputation src)
{
  //
  // Create a new histogram definition in the task,
  // which is then returned to the user for its configuration
  //

  Int_t n = fHistograms.GetEntries();
  AliRsnMiniOutput *newDef = new (fHistograms[n]) AliRsnMiniOutput(name, type, src);

  return newDef;
}

//----------------------------------------------------------------------------------
AliRsnMiniOutput *AliRsnMiniTaskPhiVn::CreateOutput(const char *name, const char *outType, const char *compType)
{
  //
  // Create a new histogram definition in the task,
  // which is then returned to the user for its configuration
  //

  Int_t n = fHistograms.GetEntries();
  AliRsnMiniOutput *newDef = new (fHistograms[n]) AliRsnMiniOutput(name, outType, compType);

  return newDef;
}
