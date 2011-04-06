//
// Implementation file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#include "Riostream.h"
#include <iomanip>

#include "TH1.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include "AliLog.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include "AliITSPIDResponse.h"
#include "AliCFContainer.h"
#include "AliVEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliMultiInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliRsnEvent.h"
#include "AliRsnTarget.h"
#include "AliRsnDaughter.h"

#include "AliRsnAnalysisPhi7TeV.h"

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::AliRsnAnalysisPhi7TeV(const char *name, Bool_t isMC) :
   AliAnalysisTaskSE(name),
   fIsMC(isMC),
   fAddTPC(kTRUE),
   fAddITS(kFALSE),
   fCheckITS(kTRUE),
   fCheckTPC(kTRUE),
   fCheckTOF(kTRUE),
   fStep(1000),
   fPhiMass(1.019455),
   fKaonMass(0.49677),
   fMaxVz(10.0),
   fITSNSigma(3.0),
   fTPCMomentumThreshold(0.350),
   fTOFNSigma(3.0),
   fESDtrackCutsTPC(0),
   fESDtrackCutsITS(0),
   fESDpid(0),
   fOutList(0x0),
   fCFunlike(0),
   fCFlikePP(0),
   fCFlikeMM(0),
   fCFtrues(0),
   fCFkaons(0),
   fHEvents(0x0)
{
//
// Constructor
//

   fTPCNSigma[0] = 5.0;
   fTPCNSigma[1] = 3.0;
   
   fVertexX[0] = fVertexX[1] = 0;
   fVertexY[0] = fVertexZ[1] = 0;
   fVertexZ[0] = fVertexZ[1] = 0;
   
   DefineOutput(1, TList::Class());
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::AliRsnAnalysisPhi7TeV(const AliRsnAnalysisPhi7TeV& copy) :
   AliAnalysisTaskSE(copy),
   fIsMC(copy.fIsMC),
   fAddTPC(copy.fAddTPC),
   fAddITS(copy.fAddITS),
   fCheckITS(copy.fCheckITS),
   fCheckTPC(copy.fCheckTPC),
   fCheckTOF(copy.fCheckTOF),
   fStep(copy.fStep),
   fPhiMass(copy.fPhiMass),
   fKaonMass(copy.fKaonMass),
   fMaxVz(copy.fMaxVz),
   fITSNSigma(copy.fITSNSigma),
   fTPCMomentumThreshold(copy.fTPCMomentumThreshold),
   fTOFNSigma(copy.fTOFNSigma),
   fESDtrackCutsTPC(copy.fESDtrackCutsTPC),
   fESDtrackCutsITS(copy.fESDtrackCutsITS),
   fESDpid(copy.fESDpid),
   fOutList(0x0),
   fCFunlike(0),
   fCFlikePP(0),
   fCFlikeMM(0),
   fCFtrues(0),
   fCFkaons(0),
   fHEvents(0x0)
{
//
// Copy constructor
//

   fTPCNSigma[0] = copy.fTPCNSigma[0];
   fTPCNSigma[1] = copy.fTPCNSigma[1];
   
   fVertexX[0] = fVertexX[1] = 0;
   fVertexY[0] = fVertexZ[1] = 0;
   fVertexZ[0] = fVertexZ[1] = 0;
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV& AliRsnAnalysisPhi7TeV::operator=(const AliRsnAnalysisPhi7TeV& copy)
{
//
// Assignment operator
//

   fIsMC = copy.fIsMC;
   fAddTPC = copy.fAddTPC;
   fAddITS = copy.fAddITS;
   fCheckITS = copy.fCheckITS;
   fCheckTPC = copy.fCheckTPC;
   fCheckTOF = copy.fCheckTOF;
   fStep = copy.fStep;
   fPhiMass = copy.fPhiMass;
   fKaonMass = copy.fKaonMass;
   fMaxVz = copy.fMaxVz;
   fITSNSigma = copy.fITSNSigma;
   fTPCNSigma[0] = copy.fTPCNSigma[0];
   fTPCNSigma[1] = copy.fTPCNSigma[1];
   fTPCMomentumThreshold = copy.fTPCMomentumThreshold;
   fTOFNSigma = copy.fTOFNSigma;
   fESDtrackCutsTPC = copy.fESDtrackCutsTPC;
   fESDtrackCutsITS = copy.fESDtrackCutsITS;
   fESDpid = copy.fESDpid;

   return (*this);
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::~AliRsnAnalysisPhi7TeV()
{
//
// Destructor
//
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::UserCreateOutputObjects()
{
//
// Create the output data container
//

   // create output list
   OpenFile(1);
   fOutList = new TList;
   
   // event related histograms
   fHEvents    = new TH1I("hEvents", "Event details", kVertexTypes, 0, kVertexTypes);
   fVertexX[0] = new TH1F("hVertexTRKX", "X position of primary vertex (tracks)", 200,  -2,  2);
   fVertexY[0] = new TH1F("hVertexTRKY", "Y position of primary vertex (tracks)", 200,  -2,  2);
   fVertexZ[0] = new TH1F("hVertexTRKZ", "Z position of primary vertex (tracks)", 400, -20, 20);
   fVertexX[1] = new TH1F("hVertexSPDX", "X position of primary vertex (SPD)"   , 200,  -2,  2);
   fVertexY[1] = new TH1F("hVertexSPDY", "Y position of primary vertex (SPD)"   , 200,  -2,  2);
   fVertexZ[1] = new TH1F("hVertexSPDZ", "Z position of primary vertex (SPD)"   , 400, -20, 20);
   
   fHEvents->GetXaxis()->SetBinLabel(1, "Good vertex with tracks");
   fHEvents->GetXaxis()->SetBinLabel(2, "Good vertex with SPD");
   fHEvents->GetXaxis()->SetBinLabel(3, "Far vertex with tracks");
   fHEvents->GetXaxis()->SetBinLabel(4, "Far vertex with SPD");
   fHEvents->GetXaxis()->SetBinLabel(5, "No good vertex");
   fHEvents->GetXaxis()->SetBinLabel(6, "Empty event");

   // define axes:
   // pairs [0] = inv. mass
   //       [1] = inv. mass resolution
   //       [2] = pt
   //       [3] = rec. multiplicity
   //       [4] = MC multiplicity (if available)
   //
   // kaons [0] = pt
   //       [1] = pTPC
   //       [2] = DCAr
   //       [3] = DCAz
   //       [4] = signalITS
   //       [5] = signalTPC
   //       [6] = betaTOF
   //       [7] = rec. multiplicity
   //       [8] = MC multiplicity (if available)
   Double_t minIM     =  0.9, maxIM    =   1.4;
   Double_t minIMRes  = -5.0, maxIMRes =   5.0;
   Double_t minMom    =  0.0, maxMom   =   5.0;
   Double_t minDCAr   =  0.0, maxDCAr  =   1.0;
   Double_t minDCAz   =  0.0, maxDCAz  =   5.0;
   Double_t minITS    =  0.0, maxITS   = 800.0; // raw signal
   Double_t minTPC    =  0.0, maxTPC   = 800.0; // raw signal
   Double_t minTOF    =  0.0, maxTOF   =   1.0; // beta = length / TOF time / c
   Double_t mult[]    =  {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 500.};
   Int_t    nIM       =  500;
   Int_t    nIMRes    =  100;
   Int_t    nMom      =   50;
   Int_t    nDCAr     =  100;
   Int_t    nDCAz     =  100;
   Int_t    nITS      =  800;
   Int_t    nTPC      =  800;
   Int_t    nTOF      =  100;
   Int_t    nMult     =  sizeof(mult) / sizeof(mult[0]) - 1;
   Int_t    nBins1[9] =  {nMom, nMom, nDCAr, nDCAz, nITS, nTPC, nTOF, nMult, nMult};
   Int_t    nBins2[5] =  {nIM , nIMRes, nMom, nMult, nMult};
      
   // containers for analysis pairs have 2 steps (quality, quality+PID)
   // and those for efficiency have 3 (add full MC)
   // REMINDER: to use variable bins, call container->SetBinLimits(iaxis, Double_t *array);
   fCFunlike = new AliCFContainer("CFUnlike", "", 2, 5, nBins2);
   fCFlikePP = new AliCFContainer("CFLikePP", "", 2, 5, nBins2);
   fCFlikeMM = new AliCFContainer("CFLikeMM", "", 2, 5, nBins2);
   fCFtrues  = new AliCFContainer("CFtrues" , "", 3, 5, nBins2);
   fCFkaons  = new AliCFContainer("CFkaons" , "", 3, 9, nBins1);
   
   // create a unique temp array to all pair containers
   AliCFContainer *CFpair[4];
   CFpair[0] = fCFunlike;
   CFpair[1] = fCFlikePP;
   CFpair[2] = fCFlikeMM;
   CFpair[3] = fCFtrues;
   
   // define axes:
   // pairs [0] = inv. mass
   //       [1] = inv. mass resolution
   //       [2] = pt
   //       [3] = rec. multiplicity
   //       [4] = MC multiplicity (if available)
   //
   // kaons [0] = pt
   //       [1] = pTPC
   //       [2] = DCAr
   //       [3] = DCAz
   //       [4] = signalITS
   //       [5] = signalTPC
   //       [6] = betaTOF
   //       [7] = rec. multiplicity
   //       [8] = MC multiplicity (if available)
   Int_t i;
   for (i = 0; i < 4; i++) {
      CFpair[i]->SetBinLimits(0, minIM   , maxIM   );
      CFpair[i]->SetBinLimits(1, minIMRes, maxIMRes);
      CFpair[i]->SetBinLimits(2, minMom  , maxMom  );
      CFpair[i]->SetBinLimits(3, mult);
      CFpair[i]->SetBinLimits(4, mult);
   }
   fCFkaons->SetBinLimits(0, minMom , maxMom );
   fCFkaons->SetBinLimits(1, minMom , maxMom );
   fCFkaons->SetBinLimits(2, minDCAr, maxDCAr);
   fCFkaons->SetBinLimits(3, minDCAz, maxDCAz);
   fCFkaons->SetBinLimits(4, minITS , maxITS );
   fCFkaons->SetBinLimits(5, minTPC , maxTPC );
   fCFkaons->SetBinLimits(6, minTOF , maxTOF );
   fCFkaons->SetBinLimits(7, mult);
   fCFkaons->SetBinLimits(8, mult);

   // add all outputs to the list
   fOutList->Add(fHEvents);
   fOutList->Add(fVertexX[0]);
   fOutList->Add(fVertexY[0]);
   fOutList->Add(fVertexZ[0]);
   fOutList->Add(fVertexX[1]);
   fOutList->Add(fVertexY[1]);
   fOutList->Add(fVertexZ[1]);
   fOutList->Add(fCFunlike);
   fOutList->Add(fCFlikePP);
   fOutList->Add(fCFlikeMM);
   fOutList->Add(fCFtrues);
   fOutList->Add(fCFkaons);
   
   // update histogram container
   PostData(1, fOutList);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::UserExec(Option_t *)
{
//
// Main execution function.
// Fills the fHEvents data member with the following legenda:
// 0 -- event OK, prim vertex with tracks
// 1 -- event OK, prim vertex with SPD
// 2 -- event OK but vz large
// 3 -- event bad
//

   static Int_t evNum = 0;
   evNum++;
   if (evNum % fStep == 0) AliInfo(Form("Processing event #%7d", evNum));

   // retrieve ESD event and related stack (if available)
   AliESDEvent *esd = static_cast<AliESDEvent*>(fInputEvent);
   AliMCEvent  *mc  = static_cast<AliMCEvent*>(fMCEvent);
   
   // if the AliESDpid object is not assigned, do this now
   // find an AliESDpid from input handler
   if (!fESDpid) {
      AliAnalysisManager *mgr   = AliAnalysisManager::GetAnalysisManager();
      AliVEventHandler   *evh   = mgr->GetInputEventHandler();
      AliESDInputHandler *esdin = 0x0;
      if (evh->IsA() == AliMultiInputEventHandler::Class()) {
         AliMultiInputEventHandler *multh = static_cast<AliMultiInputEventHandler*>(evh);
         AliInputEventHandler      *input = multh->GetFirstInputEventHandler();
         if (input->IsA() == AliESDInputHandler::Class()) {
            esdin = static_cast<AliESDInputHandler*>(input);
         }
      }
      else if (evh->IsA() == AliESDInputHandler::Class()) {
         esdin = static_cast<AliESDInputHandler*>(evh);
      }
      if (!esdin) {
         AliFatal("Required ESD input handler");
         return;
      }
      fESDpid = esdin->GetESDpid();
   }

   // check the event
   EVertexType eval = EventEval(esd);
   fHEvents->Fill((Int_t)eval);

   // if the event is good for analysis, process it
   if (eval == kGoodTracksPrimaryVertex || eval == kGoodSPDPrimaryVertex) {
      ProcessESD(esd, mc);
      if (mc) ProcessMC (esd, mc);
   }

   // update histogram container
   PostData(1, fOutList);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::Terminate(Option_t *)
{
//
// Terminate
//
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::EVertexType AliRsnAnalysisPhi7TeV::EventEval(AliESDEvent *esd)
{
//
// Checks if the event is good for analysis.
// Returns one of the flag values defined in the header
//

   // reject empty events
   Int_t ntracks = esd->GetNumberOfTracks();
   if (!ntracks) return kEmptyEvent;

   // get the best primary vertex:
   // first try the one with tracks
   const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
   const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
   Double_t            vzTrk = 1000000.0;
   Double_t            vzSPD = 1000000.0;
   Int_t               ncTrk = -1;
   Int_t               ncSPD = -1;
   if (vTrk) ncTrk = (Int_t)vTrk->GetNContributors();
   if (vSPD) ncSPD = (Int_t)vSPD->GetNContributors();
   if (vTrk) vzTrk = TMath::Abs(vTrk->GetZv());
   if (vSPD) vzSPD = TMath::Abs(vSPD->GetZv());
   if (vTrk && ncTrk > 0) {
      // fill the histograms
      fVertexX[0]->Fill(vTrk->GetXv());
      fVertexY[0]->Fill(vTrk->GetYv());
      fVertexZ[0]->Fill(vTrk->GetZv());

      // check VZ position
      if (vzTrk <= fMaxVz)
         return kGoodTracksPrimaryVertex;
      else
         return kFarTracksPrimaryVertex;
   } else if (vSPD && ncSPD > 0) {
      // fill the histograms
      fVertexX[1]->Fill(vSPD->GetXv());
      fVertexY[1]->Fill(vSPD->GetYv());
      fVertexZ[1]->Fill(vSPD->GetZv());

      // check VZ position
      if (vzSPD <= fMaxVz)
         return kGoodSPDPrimaryVertex;
      else
         return kFarSPDPrimaryVertex;
   } else
      return kNoGoodPrimaryVertex;
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::ProcessESD
(AliESDEvent *esd, AliMCEvent *mc)
{
//
// This function works with the ESD object
// and fills all containers
//

   // some utility arrays, declared static
   // to avoid creting/destroying at each event
   static TArrayC quality; // will be 1 or 0 depending if cut is passed or not
   static TArrayC pid;     // will be 1 or 0 depending if cut is passed or not

   // get stack, if possible
   AliStack *stack = 0x0;
   if (mc) stack = mc->Stack();

   // get number of tracks, which is used for setting arrays
   Int_t itrack, ntracks = esd->GetNumberOfTracks();
   quality.Set(ntracks);
   pid    .Set(ntracks);
   
   // loop on tracks and sets the flags
   for (itrack = 0; itrack < ntracks; itrack++) {
      
      // get track
      AliESDtrack *track = esd->GetTrack(itrack);
            
      // check cuts
      quality[itrack] = (Char_t)OkQuality(track);
      pid    [itrack] = (Char_t)OkPID(track);
   }
   
   // now, we can fill the histograms
   Int_t i1, i2, label;
   Double_t var[5];
   TLorentzVector p1[2], p2[2], sum[2];
   AliESDtrack *tr1, *tr2;
   TParticle   *pr1, *pr2;
   ETrackType type1, type2;
   
   // multiplicity is defined for whole event
   var[3] = AliESDtrackCuts::GetReferenceMultiplicity(esd, kTRUE);
   var[4] = (mc ? mc->GetNumberOfTracks() : 0.0);
   
   for (i1 = 0; i1 < ntracks; i1++) {
      
      // check only quality
      if (!quality[i1]) continue;
      
      // skip tracks which are not of a good type
      tr1 = esd->GetTrack(i1);
      type1 = TrackType(tr1);
      if (type1 >= kTrackTypes) continue;
      
      // skip neutral
      tr1 = esd->GetTrack(i1);
      if (tr1->Charge() == 0) continue;
      
      // try to retrieve the MC
      pr1 = 0x0;
      label = TMath::Abs(tr1->GetLabel());
      if (stack) {
         if (label >= 0 && label < stack->GetNtrack()) {
            pr1 = stack->Particle(label);
         }
      }
      
      // initialize 4-momenta
      p1[0].SetXYZM(tr1->Px(), tr1->Py(), tr1->Pz(), fKaonMass);
      if (pr1) p1[1].SetXYZM(pr1->Px(), pr1->Py(), pr1->Pz(), fKaonMass);
      
      for (i2 = i1+1; i2 < ntracks; i2++) {
         
         // check only quality
         if (!quality[i2]) continue;
         
         // skip tracks which are not of a good type
         tr2 = esd->GetTrack(i2);
         type2 = TrackType(tr2);
         if (type2 >= kTrackTypes) continue;
         
         // skip neutrals
         tr2 = esd->GetTrack(i2);
         if (tr2->Charge() == 0) continue;
         
         // try to retrieve MC
         pr2 = 0x0;
         label = TMath::Abs(tr2->GetLabel());
         if (stack) {
            if (label >= 0 && label < stack->GetNtrack()) {
               pr2 = stack->Particle(label);
            }
         }
         
         // initialize 4-momenta
         p2[0].SetXYZM(tr2->Px(), tr2->Py(), tr2->Pz(), fKaonMass);
         if (pr2) p1[1].SetXYZM(pr2->Px(), pr2->Py(), pr2->Pz(), fKaonMass);
         
         // sum 4-momenta
         for (Int_t i = 0; i < 2; i++) sum[i] = p1[i] + p2[i];
         
         // compute values
         // [0] = inv. mass
         // [1] = inv. mass resolution
         // [2] = pt
         // [3] = rec. multiplicity
         // [4] = MC multiplicity (if available)
         var[0] = sum[0].M();
         if (mc) var[1] = (sum[1].M() - sum[0].M()) / sum[1].M();
         var[2] = sum[0].Perp();
         
         // fill containers
         AliCFContainer *cf = 0x0;
         if (tr1->Charge() != tr2->Charge())
            cf = fCFunlike;
         else if (tr1->Charge() > 0)
            cf = fCFlikePP;
         else
            cf = fCFlikeMM;
            
         // fill steps: 0 = quality onl1, 1 = quality + PID
         cf->Fill(var, 0);
         if (pid[i1] && pid[i2]) cf->Fill(var, 1);
      }
   }
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::ProcessMC(AliESDEvent *esd, AliMCEvent *mc)
{
//
// Searches true pairs and fill efficiency for kaons
//

   // get stack, if possible
   AliStack *stack = 0x0;
   if (mc) stack = mc->Stack(); else return;
   
   // count tracks
   Int_t ipart, itrack, npart = stack->GetNprimary();
   AliESDtrack *matched[2] = {0x0, 0x0};
   Double_t var1[10], var2[5];
   Float_t  b[2], bcov[3];
   
   var1[7] = var2[3] = AliESDtrackCuts::GetReferenceMultiplicity(esd, kTRUE);
   var1[8] = var2[4] = npart;

   // loop on kaons
   for (ipart = 0; ipart < npart; ipart++) {
      TParticle *part = stack->Particle(ipart);
      if (!stack->IsPhysicalPrimary(ipart)) continue;
      if (TMath::Abs(part->GetPdgCode()) != 321) continue;
      
      // search best reconstructed track
      // which matches that particle ID
      // and checks it against cuts
      matched[0] = 0x0;
      itrack = BestMatchedTrack(esd, ipart);
      if (itrack >= 0) matched[0] = esd->GetTrack(itrack);
      
      // compute all variables fo single-tracks
      // [0] = pt
      // [1] = pTPC
      // [2] = DCAr
      // [3] = DCAz
      // [4] = signalITS
      // [5] = signalTPC
      // [6] = betaTOF
      // [7] = rec. multiplicity
      // [8] = MC multiplicity (if available)
      var1[0] = part->Pt();
      var1[1] = -1.0;
      for (Int_t i = 0; i < 6; i++) var1[i] = -1.0;
      if (matched[0]) {
         matched[0]->GetImpactParameters(b, bcov);
         if (IsTPC(matched[0])) var1[1] = matched[0]->GetInnerParam()->P();
         var1[2] = TMath::Abs((Double_t)b[0]);
         var1[3] = TMath::Abs((Double_t)b[1]);
         var1[4] = matched[0]->GetITSsignal();
         var1[5] = matched[0]->GetTPCsignal();
         if (MatchTOF(matched[0])) {
            var1[6]  = matched[0]->GetIntegratedLength();
            var1[6] /= matched[0]->GetTOFsignal();
            var1[6] /= 0.299792458E8;
         }
      }
      
      // fill container
      fCFkaons->Fill(var1, 0);
      
      // replace MC momentum with rec
      if (matched[0]) {
         var1[0] = matched[0]->Pt();
         if (OkQuality(matched[0])) {
            fCFkaons->Fill(var1, 1);
            if (OkPID(matched[0])) fCFkaons->Fill(var1, 2);
         }
      }
   }
   
   // loop on phi's
   Int_t id1, id2;
   TParticle *d1, *d2;
   TLorentzVector p1[2], p2[2], sum[2];
   for (ipart = 0; ipart < npart; ipart++) {
      TParticle *part = stack->Particle(ipart);
      
      // check that particle is a phi
      // and that has decayed into kaons
      if (TMath::Abs(part->GetPdgCode()) != 333) continue;
      id1 = part->GetFirstDaughter();
      id2 = part->GetLastDaughter();
      d1 = d2 = 0x0;
      if (id1 >= 0 && id1 < npart) d1 = stack->Particle(id1);
      if (id2 >= 0 && id2 < npart) d2 = stack->Particle(id2);
      if (!d1 || !d2) continue;
      if (TMath::Abs(d1->GetPdgCode()) != 321) continue;
      if (TMath::Abs(d2->GetPdgCode()) != 321) continue;
      
      // get matched tracks
      matched[0] = matched[1] = 0x0;
      itrack = BestMatchedTrack(esd, id1);
      if (itrack >= 0) matched[0] = esd->GetTrack(itrack);
      itrack = BestMatchedTrack(esd, id2);
      if (itrack >= 0) matched[1] = esd->GetTrack(itrack);
      
      // compute values
      // [0] = inv. mass
      // [1] = inv. mass resolution
      // [2] = pt
      // [3] = rec. multiplicity
      // [4] = MC multiplicity (if available)
      p1[0].SetXYZM(d1->Px(), d1->Py(), d1->Pz(), fKaonMass);
      p2[0].SetXYZM(d2->Px(), d2->Py(), d2->Pz(), fKaonMass);
      p1[1].SetXYZM(0.0, 0.0, 0.0, 0.0);
      p2[1].SetXYZM(0.0, 0.0, 0.0, 0.0);
      if (matched[0]) p1[1].SetXYZM(matched[0]->Px(), matched[0]->Py(), matched[0]->Pz(), fKaonMass);
      if (matched[1]) p2[1].SetXYZM(matched[1]->Px(), matched[1]->Py(), matched[1]->Pz(), fKaonMass);
      for (Int_t i = 0; i < 2; i++) sum[i] = p1[i] + p2[i];
      var2[0] = sum[0].M();
      var2[1] = (sum[0].M() - sum[1].M()) / sum[0].M();
      var2[2] = sum[0].Perp();
      
      // fill container
      fCFtrues->Fill(var2, 0);
      
      // replace MC momenta with rec
      if (matched[0] && matched[1]) {
         var2[0] = sum[1].M();
         var2[2] = sum[1].Perp();
         if (OkQuality(matched[0]) && OkQuality(matched[1])) {
            fCFtrues->Fill(var2, 1);
            if (OkPID(matched[0]) && OkPID(matched[1])) fCFtrues->Fill(var2, 2);
         }
      }
   }
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhi7TeV::OkPIDITS(AliESDtrack *track, AliPID::EParticleType pid)
{
//
// Check ITS particle identification with 3sigma cut
//

   // reject not ITS standalone tracks
   if (!IsITS(track)) return kFALSE;

   // count PID layers and reject if they are too few
   Int_t   k, nITSpidLayers = 0;
   UChar_t itsCluMap = track->GetITSClusterMap();
   for (k = 2; k < 6; k++) if (itsCluMap & (1 << k)) ++nITSpidLayers;
   if (nITSpidLayers < 3) {
      AliDebug(AliLog::kDebug + 2, "Rejecting track with too few ITS pid layers");
      return kFALSE;
   }

   // check the track type (ITS+TPC or ITS standalone)
   // and reject it if it is of none of the allowed types
   Bool_t isSA = kFALSE;
   if (IsTPC(track)) isSA = kFALSE;
   else if (IsITS(track)) isSA = kTRUE;
   else {
      AliWarning("Track is neither ITS+TPC nor ITS standalone");
      return kFALSE;
   }

   // create the PID response object and compute nsigma
   AliITSPIDResponse &itsrsp = fESDpid->GetITSResponse();
   Double_t mom    = track->P();
   Double_t nSigma = itsrsp.GetNumberOfSigmas(mom, track->GetITSsignal(), pid, nITSpidLayers, isSA);

   // evaluate the cut
   Bool_t ok = (TMath::Abs(nSigma) <= fITSNSigma);

   // debug message
   AliDebug(AliLog::kDebug + 2, Form("ITS nsigma = %f -- max = %f -- cut %s", nSigma, fITSNSigma, (ok ? "passed" : "failed")));

   // outcome
   return ok;
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhi7TeV::OkPIDTPC(AliESDtrack *track, AliPID::EParticleType pid)
{
//
// Check TPC particle identification with {3|5}sigmacut,
// depending on the track total momentum.
//

   // reject not TPC tracks
   if (!IsTPC(track)) return kFALSE;
   if (!track->GetInnerParam()) return kFALSE;

   // setup TPC PID response
   AliTPCPIDResponse &tpcrsp = fESDpid->GetTPCResponse();

   // get momentum and number of sigmas and choose the reference band
   Double_t mom       = track->GetInnerParam()->P();
   Double_t nSigma    = tpcrsp.GetNumberOfSigmas(mom, track->GetTPCsignal(), track->GetTPCsignalN(), pid);
   Double_t maxNSigma = fTPCNSigma[0];
   if (mom > fTPCMomentumThreshold) maxNSigma = fTPCNSigma[1];

   // evaluate the cut
   Bool_t ok = (TMath::Abs(nSigma) <= maxNSigma);

   // debug message
   AliDebug(AliLog::kDebug + 2, Form("TPC nsigma = %f -- max = %f -- cut %s", nSigma, maxNSigma, (ok ? "passed" : "failed")));

   // outcome
   return ok;
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhi7TeV::OkPIDTOF(AliESDtrack *track, AliPID::EParticleType pid)
{
//
// Check TOF particle identification if matched there.
//

   // reject not TOF-matched tracks
   if (!MatchTOF(track)) return kFALSE;

   // setup TOF PID response
   AliTOFPIDResponse &tofrsp = fESDpid->GetTOFResponse();

   // get info for computation
   Double_t momentum = track->P();
   Double_t time     = track->GetTOFsignal();
   Double_t timeint[AliPID::kSPECIES];
   tofrsp.GetStartTime(momentum);
   track->GetIntegratedTimes(timeint);

   // check the cut
   Double_t timeDiff = time - timeint[(Int_t)pid];
   Double_t sigmaRef = tofrsp.GetExpectedSigma(momentum, timeint[(Int_t)pid], AliPID::ParticleMass(pid));
   Double_t nSigma   = timeDiff / sigmaRef;

   // evaluate the cut
   Bool_t ok = (TMath::Abs(nSigma) <= fTOFNSigma);

   // debug message
   AliDebug(AliLog::kDebug + 2, Form("TOF nsigma = %f -- nsigma = %f -- cut %s", nSigma, fTOFNSigma, (ok ? "passed" : "failed")));

   // outcome
   return ok;
}

//______________________________________________________________________________
Int_t AliRsnAnalysisPhi7TeV::MatchedTrack(AliESDEvent *esd, Int_t label, Int_t &npassed, Int_t start)
{
//
// Starting from 'start' searches a track in ESD which has label 'label',
// and if it is found, records how many cuts it passes:
// 0 = none
// 1 = only quality
// 2 = quality and PID
// Returns the ID of that track (-1 if not found)
//

   Int_t it, ntracks = esd->GetNumberOfTracks();
   AliESDtrack *track = 0x0;
   
   for (it = start; it < ntracks; it++) {
      track = esd->GetTrack(it);
      if (TMath::Abs(track->GetLabel() != label)) continue;
      
      npassed = 0;
      if (OkQuality(track)) npassed++;
      if (OkPID    (track)) npassed++;
      return it;
   }
   
   return -1;
}

//______________________________________________________________________________
Int_t AliRsnAnalysisPhi7TeV::BestMatchedTrack(AliESDEvent *esd, Int_t label)
{
//
// Searches the best track which matches the specified label
// where 'best' means the one which passes most cuts
//

   // first attempt:
   // if it fails, there are no matched tracks
   Int_t ncuts = 0;
   Int_t itrack = MatchedTrack(esd, label, ncuts);
   /*if (itrack < 0) return -1;
   
   // if it succeeds, use it as a start for a loop
   Int_t bestCuts = ncuts;
   Int_t start    = itrack + 1;
   for (;;) {
      it = MatchedTrack(esd, label, ncuts, start);
      if (it < 0) break;
      
      if (ncuts > bestCuts) {
         bestCuts = ncuts;
         itrack = it;
         start = it + 1;
      }
   }*/
   
   return itrack;
}
