// Ntuplizer of event data into a TTree.

#include <TH1D.h>
#include <TTree.h>
#include <TObjArray.h>

#include <AliLog.h>
#include <AliVEvent.h>
#include <AliEventCuts.h>
#include <AliMultSelection.h>
#include <AliMultSelectionTask.h>
#include <AliVEventHandler.h>


// bins of fhStats
#define ANALYZED  1
#define ERRORS    2
#define WRITTEN   3
#define BADCENT   4
#define INCOMPDAQ 5

#include "AliAnalysisTaskNtuplizer.h"
#include <AliMultSelection.h>
#include <AliVVertex.h>
#include <AliVZDC.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <TSystem.h>

#include <AliAODv0.h>
#include <AliVVertex.h>
#include <AliAODEvent.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliPIDResponse.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliVVZERO.h>

float rawV0[64];


// maximum number of V0s to keep per event
#define NN 65535

Int_t nv;             // number of V0s (up to NN entries)
float x      [NN];    // X coordinate of decay point
float y      [NN];    // Y coordinate of decay point
float z      [NN];    // Z coordinate of decay point
float pt     [NN];    // transverse momentum
float eta    [NN];    // direction in pseudorapidity
float phi    [NN];    // direction in the phi angle
float dca    [NN];    // DCA of V0 to the primary vertex
float angle  [NN];    // angle between decay direction and momentum direction
float path   [NN];    // distance travelled before decay
float anglepn[NN];    // angle between the daughters
float dcapn  [NN];    // DCA between the daughters
Short_t v    [NN];    // particle type: -1/anti-Lambda0, 0/kaon, 1/Lambda0
float m      [NN];    // invariant mass of two daughters
float ptp    [NN];    // positive track: transverse momentum
float etap   [NN];    // positive track: direction in pseudorapidity
float phip   [NN];    // positive track: direction in the phi angle
float dcap   [NN];    // positive track: DCA to primary vertex
float nrowp  [NN];    // positive track: number of TPC hits
float nsgmp  [NN];    // positive track: nsigmas energy loss in TPC
float ptn    [NN];    // negative track: transverse momentum
float etan   [NN];    // negative track: direction in pseudorapidity
float phin   [NN];    // negative track: direction in the phi angle
float dcan   [NN];    // negative track: DCA to primary vertex
float nrown  [NN];    // negative track: number of TPC hits
float nsgmn  [NN];    // negative track: nsigmas energy loss in TPC

AliPIDResponse* PIDResponse = NULL;


// flow from TPC-only tracks (charge positive-negative, eta left-right)
struct flowTPC_t {
   float xpl[7], ypl[7];
   float xpr[7], ypr[7];
   float xnl[7], ynl[7];
   float xnr[7], ynr[7];
} __attribute__((packed)) flowTPC;



struct rawZDC_t {
   float nA[5];
   float nC[5];
} __attribute__((packed)) rawZDC;



Int_t    run;            // run number
float    magfld;         // magnetic field
UInt_t   orbit;          // linear (non-wrapped) orbit number
UShort_t bx;             // bunch crossing number
Int_t    numSPDClus1;    // number of fired SPD clusters, layer 1
Int_t    numSPDClus2;    // number of fired SPD clusters, layer 2
Int_t    numSPDTrklets;  // number of reconstructed SPD-only tracks

float centV0;      // centrality from V0 multiplicity
float centCL0;     // centrality from SPD first layer
float centCL1;     // centrality from SPD second layer

// best and SPD primary vertices
struct vtx_t {
   float vx, vy, vz;    // X, Y and Z coordinates
   Int_t vncontr;       // number of contributors
   float vchindf;       // chi^2/NDF
   float vcovzz;        // ZZ component of covariance matrix
   Int_t vtype;         // 0/unknown, 1/trksConstr, 2/zConstr, 3/trksNoConstr, 4/noTitle
} __attribute__((packed)) vtx, vtxSPD;


ClassImp(AliAnalysisTaskNtuplizer)

//_____________________________________________________________________________________________
AliAnalysisTaskNtuplizer::AliAnalysisTaskNtuplizer() :
   AliAnalysisTaskSE(),
   fhStats(NULL),
   fHistos(NULL),
   fTree(NULL)
{
   // ROOT I/O constructor.
}

//_____________________________________________________________________________________________
AliAnalysisTaskNtuplizer::AliAnalysisTaskNtuplizer(const char* name) :
   AliAnalysisTaskSE(name),
   fhStats(NULL),
   fHistos(NULL),
   fTree(NULL)
{
   DefineOutput(1, TObjArray::Class());
   DefineOutput(2, TTree::Class());
}

//_____________________________________________________________________________________________
AliAnalysisTaskNtuplizer::~AliAnalysisTaskNtuplizer()
{
   if (fHistos) delete fHistos;
   if (fTree) delete fTree;
}
//_________________________________________________________________________________________
void init_event(TTree* fTree)
{
   fTree->Branch("run",    &run);
   fTree->Branch("magfld", &magfld);
   fTree->Branch("orbit",  &orbit);
   fTree->Branch("bx",     &bx);
   fTree->Branch("numSPDClus1", &numSPDClus1);
   fTree->Branch("numSPDClus2", &numSPDClus2);
   fTree->Branch("numSPDTrklets", &numSPDTrklets);
}
//_________________________________________________________________________________________
void fill_event(AliVEvent* vEvent, AliAnalysisTaskNtuplizer* task)
{
   UInt_t period = vEvent->GetPeriodNumber();
   UInt_t orbit24 = vEvent->GetOrbitNumber(); // wrapped down to 24 bits

   if (period > 255) { // 8 bits
      task->WarnIncErr("invalid vEvent->GetPeriodNumber()");
      period = 255;
      orbit24 = (1<<24) - 1;
   }
   if (orbit24 >= (1<<24)) { // 24 bits
      task->WarnIncErr("invalid vEvent->GetOrbitNumber()");
      period = 255;
      orbit24 = (1<<24) - 1;
   }

   run    = vEvent->GetRunNumber();
   magfld = (float) vEvent->GetMagneticField();
   orbit  = period * (1<<24) + orbit24;
   bx     = vEvent->GetBunchCrossNumber();

   numSPDClus1 = vEvent->GetNumberOfITSClusters(0);
   numSPDClus2 = vEvent->GetNumberOfITSClusters(1);
   numSPDTrklets = vEvent->GetMultiplicity()->GetNumberOfTracklets();
}
//_____________________________________________________________________________________________
void init_cent(TTree* fTree)
{
   fTree->Branch("centV0",  &centV0);
   fTree->Branch("centCL0", &centCL0);
   fTree->Branch("centCL1", &centCL1);
}

//_____________________________________________________________________________________________
bool fill_cent(AliVEvent* vEvent, AliAnalysisTaskNtuplizer* task)
{
   AliMultSelection* cent = dynamic_cast<AliMultSelection*>(vEvent->FindListObject("MultSelection"));
   if (!cent) {
      task->WarnIncErr("cent is NULL");
      return true;
   }

   centV0  = cent->GetMultiplicityPercentile("V0M");
   centCL0 = cent->GetMultiplicityPercentile("CL0");
   centCL1 = cent->GetMultiplicityPercentile("CL1");

   return false;
}
//_____________________________________________________________________________________________
void init_vtx(TTree* fTree, bool fillSPD)
{
   fTree->Branch("vtx", &vtx, "vx/F:vy/F:vz/F:vncontr/I:vchindf/F:vcovzz/F:vtype/I");

   if (fillSPD)
      fTree->Branch("vtxSPD", &vtxSPD, "vx2/F:vy2/F:vz2/F:vncontr2/I:vchindf2/F:vcovzz2/F:vtype2/I");
}

//_____________________________________________________________________________________________
void FillPrimaryVertex(vtx_t& vtx, const AliVVertex* pv, AliAnalysisTaskNtuplizer* task)
{
   // Fills vtx with pv info.

   // primary vertex
   vtx.vx = (float) pv->GetX();
   vtx.vy = (float) pv->GetY();
   vtx.vz = (float) pv->GetZ();

   vtx.vncontr = pv->GetNContributors();
   vtx.vchindf = (float) pv->GetChi2perNDF();

   double covmtx[6];
   pv->GetCovarianceMatrix(covmtx);
   vtx.vcovzz = (float) covmtx[5];

   if (strcmp(pv->GetTitle(), "VertexerTracksWithConstraint") == 0)
      vtx.vtype = 1;
   else if (strcmp(pv->GetTitle(), "vertexer: Z") == 0)
      vtx.vtype = 2;
   else if (strcmp(pv->GetTitle(), "VertexerTracksNoConstraint") == 0)
      vtx.vtype = 3;
   else if (strcmp(pv->GetTitle(), "") == 0)
      vtx.vtype = 4;
   else {
      task->WarnIncErr(Form("Unknown pv title \"%s\"", pv->GetTitle()));
      vtx.vtype = 0;
   }
}

//_____________________________________________________________________________________________
bool fill_vtx(AliVEvent* vEvent, AliAnalysisTaskNtuplizer* task)
{
   const AliVVertex* pv    = vEvent->GetPrimaryVertex();
   const AliVVertex* pvSPD = vEvent->GetPrimaryVertexSPD();

   if (!pv || !pvSPD)
      return true;

   FillPrimaryVertex(vtx,    pv,    task);
   FillPrimaryVertex(vtxSPD, pvSPD, task);

   return false;
}
//_____________________________________________________________________________________________
void init_zdc(TTree* fTree)
{
   fTree->Branch("rawZDC", &rawZDC, "znA[5]/F:znC[5]/F");
}

//_____________________________________________________________________________________________
void fill_zdc(AliVEvent* vEvent, AliAnalysisTaskNtuplizer* task)
{
   AliVZDC* zdc = vEvent->GetZDCData();
   if (!zdc) {
      task->WarnIncErr("zdc is NULL, rawZDC[] filled with zeroes");
      memset(&rawZDC, 0, sizeof(rawZDC));
      return;
   }

   for (int i = 0; i < 4; i++) {
      rawZDC.nA[i] = (float) zdc->GetZNATowerEnergy()[i + 1];
      rawZDC.nC[i] = (float) zdc->GetZNCTowerEnergy()[i + 1];
   }
   // 5th channel = total deposited energies
   rawZDC.nA[4] = (float) zdc->GetZNATowerEnergy()[0];
   rawZDC.nC[4] = (float) zdc->GetZNCTowerEnergy()[0];
}
//_____________________________________________________________________________________________
void init_tpcflow(TTree* fTree)
{
   fTree->Branch("flowTPC", &flowTPC,
                 "txpl[7]:typl[7]:txpr[7]:typr[7]:txnl[7]:tynl[7]:txnr[7]:tynr[7]");
}

//_____________________________________________________________________________________________
void fill_tpcflow(AliVEvent* vEvent, AliAnalysisTaskNtuplizer* task)
{
   // Fills flow info from TPC-only tracks.

   memset(&flowTPC, 0, sizeof(flowTPC));

   AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(vEvent);
   if (!aodEvent) {
      task->WarnIncErr("aodEvent is NULL, TPC flow filled with zeroes");
      return;
   }

   // harmonics loop
   for (int n = 0; n < 7; n++) {
      for (Int_t i = 0; i < aodEvent->GetNumberOfTracks(); i++) {
         AliAODTrack* trk = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
         if (!trk) {
            task->WarnIncErr("trk is NULL, skipped");
            continue;
         }

//          bool skip = false;
//
//          if (trk->GetKinkIndex(0) > 0)
//             skip = true;
//          if (!(trk->GetStatus() & AliESDtrack::kTPCrefit))
//             skip = true;
//
//          float TPCCrossedRows = trk->GetTPCClusterInfo(2, 1);
//          UShort_t TPCFindableCls = trk->GetTPCNclsF();
//
//          if (TPCCrossedRows < 70 || TPCFindableCls < 1)
//             skip = true;
//          if (TPCCrossedRows < 0.8 * TPCFindableCls)
//             skip = true;
//
//          if (skip)
//             continue;

         // TPC-only primary vertex constrained tracks
         if (!trk->TestFilterBit(128)) continue;

         if (trk->Pt() < 0.2 || trk->Pt() > 3) continue;
         if (fabs(trk->Eta()) > 0.8) continue;

         double c = (n == 0 ? 1         : cos(n * trk->Phi()));
         double s = (n == 0 ? trk->Pt() : sin(n * trk->Phi()));

         if (trk->Eta() < -0.1 && trk->Charge() < 0) {
            flowTPC.xnl[n] += c;
            flowTPC.ynl[n] += s;
         } else if (trk->Eta() < -0.1 && trk->Charge() > 0) {
            flowTPC.xpl[n] += c;
            flowTPC.ypl[n] += s;
         } else if (trk->Eta() > 0.1 && trk->Charge() < 0) {
            flowTPC.xnr[n] += c;
            flowTPC.ynr[n] += s;
         } else if (trk->Eta() > 0.1 && trk->Charge() > 0) {
            flowTPC.xpr[n] += c;
            flowTPC.ypr[n] += s;
         }
      }
   }
}
//_____________________________________________________________________________________________
void init_lambdas(TTree* fTree)
{
   // initialize PIDResponse
   AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Printf("FATAL: mgr is NULL");
      gSystem->Exit(1);
   }

   AliInputEventHandler* inh = dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
   if (!inh) {
      Printf("FATAL: inh is NULL");
      gSystem->Exit(1);
   }

   PIDResponse = inh->GetPIDResponse();
   if (!PIDResponse) {
      Printf("FATAL: PIDResponse is NULL");
      gSystem->Exit(1);
   }

   fTree->Branch("nv",      &nv);
   //fTree->Branch("x",       x,       "x[nv]/F");
   //fTree->Branch("y",       y,       "y[nv]/F");
   //fTree->Branch("z",       z,       "z[nv]/F");
   fTree->Branch("pt",      pt,      "pt[nv]/F");
   fTree->Branch("eta",     eta,     "eta[nv]/F");
   fTree->Branch("phi",     phi,     "phi[nv]/F");
   fTree->Branch("dca",     dca,     "dca[nv]/F");
   fTree->Branch("angle",   angle,   "angle[nv]/F");
   fTree->Branch("path",    path,    "path[nv]/F");

   //fTree->Branch("anglepn", anglepn, "anglepn[nv]/F");
   fTree->Branch("dcapn",   dcapn,   "dcapn[nv]/F");
   fTree->Branch("v",       v,       "v[nv]/S");
   fTree->Branch("m",       m,       "m[nv]/F");

   fTree->Branch("ptp",     ptp,     "ptp[nv]/F");
   fTree->Branch("etap",    etap,    "etap[nv]/F");
   fTree->Branch("phip",    phip,    "phip[nv]/F");
   fTree->Branch("dcap",    dcap,    "dcap[nv]/F");
   //fTree->Branch("nrowp",   nrowp,   "nrowp[nv]/F");
   fTree->Branch("nsgmp",   nsgmp,   "nsgmp[nv]/F");

   fTree->Branch("ptn",     ptn,     "ptn[nv]/F");
   fTree->Branch("etan",    etan,    "etan[nv]/F");
   fTree->Branch("phin",    phin,    "phin[nv]/F");
   fTree->Branch("dcan",    dcan,    "dcan[nv]/F");
   //fTree->Branch("nrown",   nrown,   "nrown[nv]/F");
   fTree->Branch("nsgmn",   nsgmn,   "nsgmn[nv]/F");
}

//_____________________________________________________________________________________________
void fill_lambdas(AliVEvent* vEvent, AliAnalysisTaskNtuplizer* task)
{
   // Fills tree branches with secondary neutral vertices.

   // cleanup from previous event
   nv = 0;

   AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(vEvent);
   if (!aodEvent) {
      task->WarnIncErr("aodEvent is NULL, lambdas not filled");
      return;
   }

   const AliVVertex* pv = vEvent->GetPrimaryVertex();
   if (!pv) {
      task->WarnIncErr("pv is NULL, lambdas not filled");
      return;
   }

   // V0 loop
   for (Int_t i = 0; i < vEvent->GetNumberOfV0s(); i++) {
      if (nv >= NN) {
         task->WarnIncErr("too many V0s, some V0s skipped");
         break;
      }

      AliAODv0* aodV0 = aodEvent->GetV0(i);
      if (!aodV0) {
         task->WarnIncErr("aodV0 is NULL, V0 skipped");
         continue;
      }

       // on-the-fly V0s
      if (!aodV0->GetOnFlyStatus())
         continue;

      if (aodV0->GetNProngs() != 2) {
         task->WarnIncErr("aodV0->GetNProngs() != 2, V0 skipped");
         continue;
      }
      

      // NOTE: trkN must correspond to aodV0->*Neg*() calls below
      AliAODTrack* trkN = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
      AliAODTrack* trkP = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
      if (!trkN || !trkP) {
         task->WarnIncErr("trkN or trkP is NULL, V0 skipped");
         continue;
      }
      
         // skip kinks
      if (trkN->GetKinkIndex(0) > 0 || trkP->GetKinkIndex(0) > 0)
         continue;

      // require TPC refit on tracks
      if (!(trkP->GetStatus() & AliESDtrack::kTPCrefit) ||
          !(trkN->GetStatus() & AliESDtrack::kTPCrefit))
         continue;


      int chargeN = trkN->Charge();
      int chargeP = trkP->Charge();
      if (abs(chargeN) != 1 || abs(chargeP) != 1) {
         task->WarnIncErr("Invalid charge of a daughter, V0 skipped");
         continue;
      }

      if (chargeN != -chargeP)
         continue;

     
      TVector3 vectN, vectP;
      float dca_, dcapn_, dcap_, dcan_;

      dca_   = (float) aodV0->DcaV0ToPrimVertex();
      dcapn_ = (float) aodV0->DcaV0Daughters();

      if (chargeN < 0) {
         dcap_ = (float) aodV0->DcaPosToPrimVertex();
         dcan_ = (float) aodV0->DcaNegToPrimVertex();
         vectP.SetXYZ(aodV0->MomPosX(), aodV0->MomPosY(), aodV0->MomPosZ());
         vectN.SetXYZ(aodV0->MomNegX(), aodV0->MomNegY(), aodV0->MomNegZ());
      } else { // swap daughters
         dcap_ = (float) aodV0->DcaNegToPrimVertex();
         dcan_ = (float) aodV0->DcaPosToPrimVertex();
         vectP.SetXYZ(aodV0->MomNegX(), aodV0->MomNegY(), aodV0->MomNegZ());
         vectN.SetXYZ(aodV0->MomPosX(), aodV0->MomPosY(), aodV0->MomPosZ());

         AliAODTrack* trk3 = trkN;
         trkN = trkP;
         trkP = trk3;
      }

     

          // cuts to reduce output file size
      if (dca_ > 2) continue; // in cm //****************************************check the cut value
      if (dcap_ < 0.05 || dcan_ < 0.05) continue; // in cm //****************************************check the cut value
      if (fabs(trkN->Eta()) > 0.8 || fabs(trkP->Eta()) > 0.8) continue;//****************************************check the cut value

      float pTPCCrossedRows = trkP->GetTPCClusterInfo(2, 1);
      float nTPCCrossedRows = trkN->GetTPCClusterInfo(2, 1);
      UShort_t pTPCFindableCls = trkP->GetTPCNclsF();
      UShort_t nTPCFindableCls = trkN->GetTPCNclsF();

      // standard track quality cuts
      if (pTPCCrossedRows < 70 || nTPCCrossedRows < 70) continue;
      if (pTPCFindableCls < 1 || nTPCFindableCls < 1) continue;
      if (pTPCCrossedRows < 0.8 * pTPCFindableCls) continue;
      if (nTPCCrossedRows < 0.8 * nTPCFindableCls) continue;

      // bonus cut
      if (vectP.Pt() > 20 || vectN.Pt() > 20) continue;

      // evaluate invariant masses
      double mass[3]; // 0=anti-Lambda0, 1=kaon, 2=Lambda0
      double mass0[3] = {1.116, 0.498, 1.116}; // nominal mass values
      // double width[3] = {0.017, 0.04, 0.017};  // accepted region widths
      double width[3] = {0.013, 0.033, 0.013};  // accepted region widths

      TLorentzVector L1, L2;

      L1.SetVectM(vectN, 0.93827);
      L2.SetVectM(vectP, 0.13957);
      mass[0] = (L1 + L2).M();

      L1.SetVectM(vectN, 0.13957);
      L2.SetVectM(vectP, 0.13957);
      mass[1] = (L1 + L2).M();

      L1.SetVectM(vectN, 0.13957);
      L2.SetVectM(vectP, 0.93827);
      mass[2] = (L1 + L2).M();

      // mass hypothesis loop
      for (int k = 0; k < 3; k++) {
         if (fabs(mass[k] - mass0[k]) > width[k])
            continue;

         x[nv] = (float) aodV0->Xv();
         y[nv] = (float) aodV0->Yv();
         z[nv] = (float) aodV0->Zv();

         TVector3 vect(aodV0->Px(), aodV0->Py(), aodV0->Pz());
         pt[nv]  = (float) vect.Pt();
         eta[nv] = (float) vect.Eta();
         phi[nv] = (float) vect.Phi();

         TVector3 vectPV(aodV0->Xv() - pv->GetX(),
                         aodV0->Yv() - pv->GetY(),
                         aodV0->Zv() - pv->GetZ());
         angle[nv] = vectPV.Angle(vect);
         path[nv] = vectPV.Mag();
	 anglepn[nv] = vectP.Angle(vectN);

         // more cuts
         if (angle[nv] > 0.07) continue;//****************************************check the cut value
         if (path[nv] > 180) continue;//******************************************check the cut value

         dca[nv] = dca_;
         dcapn[nv] = dcapn_;
	 if (dcapn[nv]>0.5) continue;//****************************************check the cut value

         v[nv] = k - 1;
         m[nv] = mass[k];

         ptp[nv]  = (float) vectP.Pt();
         etap[nv] = (float) vectP.Eta();
         phip[nv] = (float) vectP.Phi();
         dcap[nv] = dcap_;
	 nrowp[nv] = pTPCCrossedRows;

         ptn[nv]  = (float) vectN.Pt();
         etan[nv] = (float) vectN.Eta();
         phin[nv] = (float) vectN.Phi();
         dcan[nv] = dcan_;
         nrown[nv] = nTPCCrossedRows;

         nsgmp[nv] = 0;
         nsgmn[nv] = 0;

         // particle identification
         AliPID::EParticleType tP = (k == 2 ? AliPID::kProton : AliPID::kPion);
         AliPID::EParticleType tN = (k == 0 ? AliPID::kProton : AliPID::kPion);

         double val;

         // skip candidates without TCP dE/dx
         if (PIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, trkP, tP, val) ==
            AliPIDResponse::kDetPidOk) {
            nsgmp[nv] = (float) val;
         } else
            continue;

         if (PIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, trkN, tN, val) ==
            AliPIDResponse::kDetPidOk) {
            nsgmn[nv] = (float) val;
         } else
            continue;

         // dE/dx cuts to reduce output file size
         if (fabs(nsgmp[nv]) > 3 || fabs(nsgmn[nv]) > 3)
            continue;

         nv++;
      }

   } // V0s
}
//_____________________________________________________________________________________________
void init_v0(TTree* fTree)
{
   fTree->Branch("rawV0", rawV0, "rawV0[64]/F");
}

//_____________________________________________________________________________________________
void fill_v0(AliVEvent* vEvent, AliAnalysisTaskNtuplizer* task)
{
   AliVVZERO* vzero = vEvent->GetVZEROData();
   if (!vzero) {
      task->WarnIncErr("vzero is NULL, rawV0[] filled with zeroes");
      memset(rawV0, 0, sizeof(rawV0));
      return;
   }

   for (int c = 0; c < 64; c++)
      rawV0[c] = vzero->GetMultiplicity(c);
}
//_____________________________________________________________________________________________
void AliAnalysisTaskNtuplizer::UserCreateOutputObjects()
{
   // Creates output histogram(s) and output TTree, initializes necessary persistent objects.

   fHistos = new TObjArray;
   fHistos->SetOwner(true);

   fhStats = new TH1D("hStats", "Event statistics", 10, 1, 11);
   fhStats->GetXaxis()->SetBinLabel(ANALYZED,  "Analyzed");
   fhStats->GetXaxis()->SetBinLabel(ERRORS,    "Errors");
   fhStats->GetXaxis()->SetBinLabel(BADCENT,   "Bad centrality");
   fhStats->GetXaxis()->SetBinLabel(INCOMPDAQ, "Incomplete DAQ");
   fhStats->GetXaxis()->SetBinLabel(WRITTEN,   "Written");
   fHistos->Add(fhStats);

   if (!OpenFile(2)) AliFatal("OpenFile(2) failed");
   fTree = new TTree("events", "TTree with filtered event data");

   // associate tree branches with variables, initialize other
   init_event(fTree);
   init_cent(fTree);
   init_vtx(fTree, true);  // false=do not fill SPD vertex
   init_zdc(fTree);
   init_v0(fTree);
   // init_fmd(fTree, fHistos);
//    init_ad0(fTree);
//   init_chtracks(fTree);
   init_tpcflow(fTree);
   init_lambdas(fTree);
//    init_chkaons(fTree);

   PostData(1, fHistos);
   PostData(2, fTree);
}

//_____________________________________________________________________________________________
void AliAnalysisTaskNtuplizer::UserExec(Option_t*)
{
   // Main analysis function.

   fhStats->Fill(ANALYZED);

   AliVEvent* Event = InputEvent();
   if (!Event) {
     WarnIncErr("vEvent is NULL, event skipped");
     PostData(1, fHistos);
     PostData(2, fTree);
     return;
   }


   vEvent = dynamic_cast<AliAODEvent*>(InputEvent());
   if (!vEvent)
     {
       return;
     }


   if (!fEventCut_d.AcceptEvent(vEvent))//all the event cuts designed for LHC15o in AliVEventCuts class                          
     {
       PostData(1, fHistos);
       PostData(2, fTree);
       return;
     }



 if (vEvent->IsIncompleteDAQ()) {
      fhStats->Fill(INCOMPDAQ);
      PostData(1, fHistos);
      PostData(2, fTree);
      return;
   }


 
   if (fill_cent(vEvent, this)) {
      fhStats->Fill(BADCENT);
      PostData(1, fHistos);
      PostData(2, fTree);
      return;
   }

   if (fill_vtx(vEvent, this)) {
      WarnIncErr("PV not filled, event skipped");
      PostData(1, fHistos);
      PostData(2, fTree);
      return;
   }

   fill_event(vEvent, this);
   fill_zdc(vEvent, this);
   fill_v0(vEvent, this);
   // fill_fmd(vEvent, this, 1);  // 1 = harmonics number
   // fill_ad0(vEvent, this);
//    fill_chtracks(vEvent, this);
    fill_tpcflow(vEvent, this);
    fill_lambdas(vEvent, this);
//    fill_chkaons(vEvent, this);

   fhStats->Fill(WRITTEN);
   fTree->Fill();

   PostData(1, fHistos);
   PostData(2, fTree);
}

//_____________________________________________________________________________________________
void AliAnalysisTaskNtuplizer::WarnIncErr(const char* msg)
{
   // Prints out message and increases error counter in fhStats.

   AliErrorF("%s", msg);
   fhStats->Fill(ERRORS);
}
