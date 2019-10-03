//-----------------------------------------------------------------
// AliAnalysisTaskV0ForRAA class
// This task is for analysing Lambda and K0s pt spectra in PbPb and
// pp as well as with MC. The flag for pp and MC  must be set
// accordingly, default is PbPb data.
// It works with ESD files only.
//-----------------------------------------------------------------

#ifndef ALIANALYSISTASKV0FORRAA_H
#define ALIANALYSISTASKV0FORRAA_H


class TH1F;
class TH2F;
//class TH3F;

class Tlist;

class AliESDv0;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpid;
class AliESDEvent;
class AliMCEvent;
class AliPIDResponse;
class AliStack;
//#include "THn.h"

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif



class AliAnalysisTaskV0ForRAA : public AliAnalysisTaskSE {
 public:
  
  AliAnalysisTaskV0ForRAA(); 
  AliAnalysisTaskV0ForRAA(const char *name);
  virtual ~AliAnalysisTaskV0ForRAA();
  
 
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);



  //-- MC truth/reco --//
  void SetMCMode(Bool_t mcmode)                               {fMCMode            = mcmode;  if(fMCMode) Printf("AliAnalysisTaskV0ForRAA::running mc mode: histos of MC reco");}
  void SetMCTruthMode(Bool_t mcmode)                          {fMCTruthMode       = mcmode;  if(fMCTruthMode) Printf("AliAnalysisTaskV0ForRAA::running mc mode: histos of MC truth");}
  void SelectInjected(Bool_t injected)                        {fSelectInjected    = injected;if(fSelectInjected) Printf("AliAnalysisTaskV0ForRAA::only injected MC particles");}
  void SelectMBMotherMC(Bool_t mbmother)                      {fSelectMBMotherMC  = mbmother;if(mbmother)  Printf("AliAnalysisTaskV0ForRAA::only MB mother MC for sec lambdas selected");}
  void SelectOnlyPosLabelMC(Bool_t poslabel)                  {fCheckNegLabelReco = poslabel;if(poslabel) Printf("AliAnalysisTaskV0ForRAA::Select only MC truth and reco with pos label reco");}
  void SelectOnlyFoundRecoV0MC(Bool_t found)                  {fOnlyFoundRecoV0   = found;   if(found) Printf("AliAnalysisTaskV0ForRAA::Select only MC truth with found reco V0");}


  //-- Centrality  --//
  // use centrality - if yes, which one
  void  SetUseCentrality(Int_t cent)                          {fUseCentrality      = cent; Printf("AliAnalysisTaskV0ForRAA::centrality selected for detector %i (0=off, 1=VZERO, 2=SPD)",cent);}
  // set range
  void  SetUseCentralityRange(Int_t range)                    {fUseCentralityRange = range;if(fUseCentrality) Printf("AliAnalysisTaskV0::centrality range %i",fUseCentralityRange);}
  // centrality bin to be used
  void  SetUseCentralityBin(Int_t bin)                        {fUseCentralityBin   = bin; if(fUseCentrality) Printf("AliAnalysisTaskV0ForRAA::centrality selected for bin %i",fUseCentralityBin); }


  //-- event cuts --//
  void SetPrimVertexZCut(Double_t vtxcut,Bool_t status)       {fVertexZCut = vtxcut;fVtxStatus = status; Printf("AliAnalysisTaskV0ForRAA::SetPrimVertexZCut %3.2f",vtxcut);}
  void SetAnapp(Bool_t anapp)                                 {fAnapp = anapp ;if(fAnapp) Printf("AliAnalysisTaskV0ForRAA::analysing pp!!!");}
  void SetRejectPileUpSPD(Bool_t rejectPU = kFALSE)           {fRejectPileUpSPD = rejectPU;if(fRejectPileUpSPD) Printf("AliAnalysisTaskV0ForRAA::reject pileup events from SDP in pp");}
  void SelectWithSDD(Bool_t sdd)                              {fSelSDD =sdd; if(sdd) Printf("AliAnalysisTaskV0ForRAA:: only events with SDD selected!");}
  void SelectWithNoSDD(Bool_t sdd)                            {fSelNoSDD =sdd; if(sdd) Printf("AliAnalysisTaskV0ForRAA:: only events with NO SDD selected!");}

  //-- track cuts --//
  void SetESDTrackCuts(Int_t ncr, Double_t chi2=4, Bool_t tpcrefit=kTRUE,Bool_t itsrefit=kFALSE)       {fNcr=ncr;fChi2cls=chi2,fTPCrefit=tpcrefit;fITSrefit=itsrefit;Printf("AliAnalysisTaskV0ForRAA::AliESDtrackCuts for V0s set ncr %i, chi2 %1.2f, TPC refit %i, ITS refit %i",ncr,chi2,tpcrefit,itsrefit);}
  void SetESDTrackCutsCharged(Int_t ncr, Double_t chi2=4, Bool_t tpcrefit=kTRUE,Bool_t itsrefit=kFALSE) {fNcrCh=ncr;fChi2clsCh=chi2,fTPCrefitCh=tpcrefit;fITSrefitCh=itsrefit;Printf("AliAnalysisTaskV0ForRAA::AliESDtrackCuts for charged particles setncr %i, chi2 %1.2f, TPC refit %i, ITS refit %i",ncr,chi2,tpcrefit,itsrefit);}
  void SetESDTrackCutsLowPt(Int_t ncr, Double_t chi2=4, Bool_t tpcrefit=kTRUE,Bool_t itsrefit=kFALSE)  {fNcrLpt=ncr;fChi2clsLpt=chi2,fTPCrefitLpt=tpcrefit;fITSrefitLpt = itsrefit;Printf("AliAnalysisTaskV0ForRAA::AliESDtrackCuts for low pt particles set ncr %i, chi2 %1.2f, TPC refit %i, ITS refit %i",ncr,chi2,tpcrefit,itsrefit);}

  void SetUseOnthefly(Bool_t useonthefly)                     {fOntheFly = useonthefly; if(!fOntheFly) Printf("AliAnalysisTaskV0ForRAA::offline V0s");}
  void SetUsePID(Bool_t usepid,Double_t nsigma=100.0,Double_t pcut=100.0,Bool_t pidpion=kFALSE,Double_t nsigma2=100.0) {fUsePID = usepid;fNSigma = nsigma;fPPIDcut = pcut; fUsePIDPion = pidpion;fNSigma2 = nsigma2; if(fUsePID) Printf("AliAnalysisTaskV0ForRAA::proton PID! of %4.2f for p: %4.2f, also pion? %i nsig2=%4.2f",fNSigma,pcut,pidpion,fNSigma2);}
  void SetCutMoreNclsThanRows(Bool_t cut)                     {fMoreNclsThanRows=cut; if(cut) Printf("AliAnalysisTaskV0ForRAA::cut on more ncls than crossed rows");}  
  void SetCutMoreNclsThanFindable(Bool_t cut)                 {fMoreNclsThanFindable=cut; if(cut) Printf("AliAnalysisTaskV0ForRAA::cut on more ncls than ncls findable");}
  void SetCutMoreNclsThanFindableMax(Bool_t cut)              {fMoreNclsThanFindableMax = cut; if(cut) Printf("AliAnalysisTaskV0ForRAA::cut on more ncls than ncls findable max");}

  void SetRatioFoundOverFindable(Double_t cut)                {fRatioFoundOverFindable = cut; Printf("AliAnalysisTaskV0ForRAA::cut on found over finable clusters %f",cut);}
  void SetRatioMaxCRowsOverFindable(Double_t cut)             {fRatioMaxCRowsOverFindable = cut;  Printf("AliAnalysisTaskV0ForRAA::cut on max crossed rows over finable clusters %f",cut);}

  void SetLowPtTPCCutAliESDTrackCut(Double_t pt)              {fPtTPCCut=pt;Printf("AliAnalysisTaskV0ForRAA::SetLowPtTPCCutAliESDTrackCut pt=%2.2f",pt);} 
   
  void SetMaxChi2PerITSCluster(Double_t chi2)                 {fChi2PerClusterITS = chi2; Printf("AliAnalysisTaskV0ForRAA::max chi2 per ITS cluster %3.2f",chi2);}
  void SetRapidityCutMother(Bool_t cut,Double_t val=5.0)      {fRapCutV0 = cut; fRap = val; if(cut) Printf("AliAnalysisTaskV0ForRAA::cut on truth secondary rapidity %2.2f",val);}
  void SetMinPt(Double_t minPt=0.0)                           {fMinPt = minPt; if(minPt>0.0) Printf("AliAnalysisTaskV0ForRAA::cut on min pt %2.2f",minPt);}
  void SetPtShift(const Double_t shiftVal) {
    //user defined shift in charge/pt
    if(shiftVal) { fShift=kTRUE; fDeltaInvP = shiftVal; Printf("AliAnalysisTaskV0ForRAA::WARNING!!!!!!!!!!!!!! pt shift introduced %2.8f!",fDeltaInvP);}
  }
  
  void SetDCAV0ToVertexK0(Double_t dcaTovertex)               {fDCAToVertexK0  = dcaTovertex; Printf("AliAnalysisTaskV0ForRAA::dca to vertex K0s %2.3f",dcaTovertex);}
  void SetDCAV0ToVertexL(Double_t dcaTovertex)                {fDCAToVertexL   = dcaTovertex; Printf("AliAnalysisTaskV0ForRAA::dca to vertex L/AL %2.3f",dcaTovertex);}
  void SetDCADaughtersL(Double_t dcaDaughters)                {fDCADaughtersL  = dcaDaughters; Printf("AliAnalysisTaskV0:ForRAA:dca daughters L %2.3f",dcaDaughters);}
  void SetDCADaughtersAL(Double_t dcaDaughters)               {fDCADaughtersAL = dcaDaughters; Printf("AliAnalysisTaskV0ForRAA::dca daughters AL %2.3f",dcaDaughters);}
  void SetDCADaughtersK0(Double_t dcaDaughters)               {fDCADaughtersK0 = dcaDaughters; Printf("AliAnalysisTaskV0ForRAA::dca daughters K0s %2.3f",dcaDaughters);}
  void SetDCADaughtersLargeToVertex(Double_t dcaDaughtersVtx) {fDCADaughtersToVtxLarge = dcaDaughtersVtx; Printf("AliAnalysisTaskV0ForRAA::dca daughters to vertex large %2.3f",dcaDaughtersVtx);}
  void SetDCADaughtersSmallToVertex(Double_t dcaDaughtersVtx) {fDCADaughtersToVtxSmall = dcaDaughtersVtx; Printf("AliAnalysisTaskV0ForRAA::dca daughters to vertex small %2.3f",dcaDaughtersVtx);}
  void SetDecayRadiusXYMinMax(Double_t decMin,Double_t decMax,Double_t pt=100000.0){fDecayRadXYMin  = decMin;fDecayRadXYMax = decMax;fPtDecRadMin =pt;Printf("AliAnalysisTaskV0ForRAA::min xy decay radius %2.3f max %2.3f for max pt %2.2f",decMin,decMax,pt);}
  void SetCosOfPointingAngleL(Double_t pointAng,Double_t ptMaxCut=100.0) {fCosPointAngL = pointAng;fCPAPtCutL = ptMaxCut;Printf("AliAnalysisTaskV0ForRAA::SetCosOfPointingAngleL %1.5f and pt max %2.2f",pointAng,ptMaxCut);} 
  void SetCosOfPointingAngleK(Double_t pointAng,Double_t ptMaxCut=100.0) {fCosPointAngK = pointAng;fCPAPtCutK0 = ptMaxCut; Printf("AliAnalysisTaskV0ForRAA::SetCosOfPointingAngleK  %1.5f and pt max %2.2f",pointAng,ptMaxCut);}
  void SetOpeningAngleCut(Double_t opang, Double_t maxpt)     {fOpengAngleDaughters=opang; fOpAngPtCut = maxpt,Printf("AliAnalysisTaskV0::cut on opening angle %1.3f up to pt= %2.2f",opang,maxpt);}

  void SetMaxDecayLength(Double_t decLength)                  {fDecayLengthMax = decLength; Printf("AliAnalysisTaskV0ForRAA::SetMaxDecayLength %2.3f",decLength);}
  void SetMinDecayLength(Double_t decLength)                  {fDecayLengthMin = decLength; Printf("AliAnalysisTaskV0ForRAA::SetMinDecayLength %2.3f",decLength);}
  void SetDCAXK0(Double_t dcaXK)                              {fDCAXK = dcaXK; Printf("AliAnalysisTaskV0ForRAA::SetDCAXK0 %2.3f",dcaXK);}
  void SetDCAYK0(Double_t dcaYK)                              {fDCAYK = dcaYK; Printf("AliAnalysisTaskV0ForRAA::SetDCAYK0 %2.3f",dcaYK);}
  void SetDCAXLambda(Double_t dcaXL)                          {fDCAXL = dcaXL; Printf("AliAnalysisTaskV0ForRAA::SetDCAXLambda %2.3f",dcaXL);}
  void SetDCAYLambda(Double_t dcaYL)                          {fDCAXL = dcaYL; Printf("AliAnalysisTaskV0ForRAA::SetDCAYLambda %2.3f",dcaYL);}
  void SetDCAZ(Double_t dcaZ)                                 {fDCAZ  = dcaZ;  Printf("AliAnalysisTaskV0ForRAA::SetDCAZ %2.3f",dcaZ);}
  void SetChi2CutKf(Double_t chi2)                            {fChiCutKf = chi2; Printf("AliAnalysisTaskV0ForRAA::SetChi2CutKf %3.2f",chi2);}
  void SetArmenterosCutAlpha(Double_t alfaMin)                {fAlfaCut = alfaMin;Printf("AliAnalysisTaskV0ForRAA::SetArmenterosCut a=%1.3f",alfaMin);}
  void SetArmenterosCutQt(Double_t ptmin,Double_t ptmax,Bool_t k0s,Bool_t la,Double_t slope=0.2,Double_t qtLinear=0.0){fQtCutPt = ptmax;fQtCutPtLow = ptmin, fArmQtSlope = slope,fArmCutK0 = k0s;fArmCutL = la;fQtCut = qtLinear;Printf("AliAnalysisTaskV0ForRAA::SetArmenterosCut ptmin = %3.2f ptmax = %3.2f. slope: %1.2f.  Is K0s? %i La? %i, qt linear: %3.2f",ptmin,ptmax,slope,k0s,la,qtLinear);}
  void SetMinMassDiffLK0s(Double_t diffK,Double_t diffL)      {fExcludeLambdaFromK0s = diffK;fExcludeK0sFromLambda = diffL; Printf("AliAnalysisTaskV0ForRAA::SetMaxMassDifferenceLK0s for K0s %1.3f  K0s for L %1.3f",diffK,diffL);}
  void SetMinMassDiffPhoton(Double_t diffK,Double_t diffL)      {fExcludePhotonsFromK0s = diffK;fExcludePhotonsFromLambda = diffL; Printf("AliAnalysisTaskV0ForRAA::SetMaxMassDifferencePhoton for K0s %1.3f  K0s for L %1.3f",diffK,diffL);}

  void SetCtauCut(Double_t ctK0s, Double_t ctL,Double_t ptK0=100.0,Double_t ptL=100.0) {fCtauK0s = ctK0s*2.6842; fCtauL = ctL*7.89;fCtauPtCutK0 = ptK0; fCtauPtCutL = ptL;
    Printf("AliAnalysisTaskV0ForRAA::SetCtauCut ctK=%2.2f, ctL = %2.2f for ptK= %5.2f ptL=%5.2f",ctK0s,ctL,ptK0,ptL);}
  void SetDoEtaOfMCDaughtersCut(Bool_t doCut,Double_t eta=5.0){fEtaCutMCDaughters = doCut; fEtaCutMCDaughtersVal=eta; Printf("AliAnalysisTaskV0ForRAA::eta cut on V0 (MC truth ? %i) daughters %1.3f !",doCut,eta);}
   void SetEtaSignCut(Double_t etasign)                        {fEtaSignCut = etasign;Printf("AliAnalysisTaskV0ForRAA::eta cut sign on  daughters %2.2f !",etasign);}
  void SetLowHighMassCut(Double_t lowK=0.25,Double_t highK=0.75,Double_t lowL=1.05,Double_t highL=1.25){fK0sLowMassCut = lowK; fK0sHighMassCut = highK; fLLowMassCut = lowL; fLHighMassCut = highL; Printf("AliAnalysisTaskV0ForRAA::SetLowHighMassCut K0s: low = %1.3f  high = %1.3f  Lambda: low = %1.3f  high = %1.3f",lowK,highK,lowL,highL);}
  void SetMinMaxNCLSITS(Int_t minP,Int_t maxP,Int_t minN,Int_t maxN,Bool_t switchCase=kFALSE,Double_t radmin=0.0000,Double_t radmax=10000.0){fMinNCLSITSPos = minP; fMaxNCLSITSPos = maxP;fMinNCLSITSNeg = minN; fMaxNCLSITSNeg = maxN;fSwitchCaseITSCls = switchCase;fDecRadCutITSMin=radmin;fDecRadCutITSMax=radmax;Printf("AliAnalysisTaskV0ForRAA::SetMinMaxNCLSITS for V0 daughters minPos %i, maxPos %i, minNeg %i, maxNeg %i switch case %i for 2D decay rad. min: %3.2f  max: %3.2f",minP,maxP,minN,maxN,switchCase,radmin,radmax);}
  
  void SetTPCTrackCutsMI(Bool_t tlength=kFALSE, Bool_t crows=kFALSE, Bool_t ncls=kFALSE,Double_t lf1=1.0,Double_t lf2=0.85){fCutMITrackLength = tlength; fCutMICrossedR=crows;  fCutMITPCncls=ncls; fCutMITrackLengthLengthF=lf1;fCutMICrossedRLengthF=lf2;Printf("AliAnalysisTaskV0ForRAA::SetTPCTrackCutsMI track length %i  crossed rows %i  ncls %i factor length %1.2f factor ncr %1.2f",fCutMITrackLength, fCutMICrossedR,fCutMITPCncls,lf1,lf2);}

  void SetFillDetHistoAL(Bool_t fillAL = kFALSE)              {fSetFillDetAL = fillAL; if(fillAL) Printf("AliAnalysisTaskV0ForRAA::SetFillDetHistoAL fill detetctor histos with AL instead L");}
  void SetFillPt(Bool_t fillpt = kFALSE)                      {fSetPtDepHist = fillpt; if(fillpt) Printf("AliAnalysisTaskV0ForRAA::SetFillPt fill pt instead of mass");}
  void SetMinDistTPCInner(Double_t dist = 1000000.0)          {fDistanceTPCInner = dist; Printf("AliAnalysisTaskV0ForRAA::SetMinDistTPCInner set dist min to %2.2f",dist); }

  void SetStopRecoLoop(Bool_t stop)                           {fStopLoop = stop;  Printf("AliAnalysisTaskV0ForRAA::SetStopRecoLoop %i",stop);}

  void SetDistDauForCheck(Double_t dist)                      {fDistDauForCheck = dist;Printf("AliAnalysisTaskV0ForRAA::SetDistDauForCheck %2.2f",dist);}

  void SetUseXiOmega(Bool_t xiM,Bool_t xi0,Bool_t omega){fUseXiM = xiM; fUseXi0 = xi0; fUseOmega =omega;Printf("AliAnalysisTaskV0ForRAA::SetUseXiOmega xi minus/plus %i xi0 %i omega %i",xiM,xi0,omega);}
  void SetDoRapCutXi(Bool_t rap ){fCutRapXi = rap;Printf("AliAnalysisTaskV0ForRAA::SetDoRapCutXi  %i <0.5",rap); }
 private:
   
  //----------------------------functions --------------------------------------------//

  void   Process();                  // process event
  void   V0RecoLoop(Int_t id0,Int_t id1,Int_t isSecd,Int_t what,Double_t ptV0MC,Int_t pdgMother,Double_t ptXiMother,Double_t decaylengthMCV0); // loop over reconstructed V0 (data or MC)
  void   V0MCTruthLoop();            // loop over MC truth V0s
  Int_t  CalculateCentralityBin();   // get the centrality bin from multiplicity
  Bool_t GetMCTruthPartner(AliESDtrack *pos,AliESDtrack *neg,Int_t id0,Int_t id1);// find MC truth partner for reconstructed track
  Bool_t CheckMultipleV0Candidates(Int_t part1,Int_t part2,Int_t iV0MI,Int_t trackID[][2]);//check if V0 was already found
  Int_t  FindPDGCode(AliStack *stackRec,AliESDtrack *trackPos,AliESDtrack *trackNeg);
  void   CheckDistanceOfDaughters(Int_t iV0MI,Float_t V0ID[],Double_t magF,Int_t particle);
  //----------------------------- objects ----------------------------------------------//

  //event
  AliESDEvent     *fESD;                //ESD event object
  AliMCEvent      *fMCev;               //MC event object

  
  //PID and track cuts
  AliPIDResponse  *fESDpid;             //pid object
  AliESDtrackCuts *fESDTrackCuts;       //esd track cuts for daughters
  AliESDtrackCuts *fESDTrackCutsCharged;//esd track cuts for all charged particles
  AliESDtrackCuts *fESDTrackCutsLowPt;  //esd track cuts for daughters at low pt

  TList           *fOutputContainer;    // output data container
   
  //----------------------------histograms --------------------------------------------//
  /*
    THnF *fTHnFK0s; 
    THnF *fTHnFL; 
    THnF *fTHnFAL; 
   
    THnF *fTHnFK0sDauEta; 
    THnF *fTHnFLDauEta; 
    THnF *fTHnFALDauEta; 
    THnF *fTHnFK0sDauPhi; 
    THnF *fTHnFLDauPhi; 
    THnF *fTHnFALDauPhi; 
  */
  //-------------------event histos -------------------//
  TH1F   *fHistITSLayerHits;                        // pp 2.76 TeV analysis: check hist on div. ITS layer
  TH1F   *fHistOneHitWithSDD;                       // pp 2.76 TeV analysis: check hist on at least one ITS layer
  TH1F   *fHistNEvents;                             // count number of events for each event cut
  TH2F   *fHistPrimVtxZESDVSNContributors;          // count contributors to ESD vertex
  TH2F   *fHistPrimVtxZESDTPCVSNContributors;       // count contributors to TPC vertex
  TH2F   *fHistPrimVtxZESDSPDVSNContributors;       // count contributors to SPD vertex

  TH1F   *fHistPrimVtxZESD;                         // primary ESD vertex position z after cuts and processing
  TH1F   *fHistPrimVtxZESDTPC;                      // primary TPC vertex position z after cuts and processing
  TH1F   *fHistPrimVtxZESDSPD;                      // primary SPD vertex position z after cuts and processing

  TH1F   *fHistESDVertexZ;                          // primary TPC vertex position z before cuts
     
  TH1F   *fHistMuliplicity;                         // number of particles from centrality selection
  TH1F   *fHistMuliplicityRaw;                      // number of particles from centrality selection before processing
  TH1F   *fHistCentBinRaw;                          // events per centralitybin before centrality selection
  TH1F   *fHistCentBin;                             // events per centralitybin
  TH1F   *fHistMultiplicityPrimary;                 // number of charged particles
  TH1F   *fHistNPrim;                               // number of contributors to the prim vertex

  //------------------------ single V0 histos --------------------------//
  //  TH3F   *fHistPiPiPhiPosVsPtPosVsMass;//xxx
  // TH3F   *fHistPiPPhiPosVsPtPosVsMass;//xxx
  //TH3F   *fHistPiAPPhiPosVsPtPosVsMass;//xxx
  TH2F   *fHistPiPiK0sVsLambdaMass;                     // K0s mass vs Lamba mass for all pt for K0s
  TH2F   *fHistPiPiK0sVsALambdaMass;                    // K0s mass vs ALamba mass for all pt for K0s
  TH2F   *fHistPiPK0sVsLambdaMass;                      // K0s mass vs Lamba mass for all pt for Lambda
  TH2F   *fHistPiAPK0sVsALambdaMass;                    // K0s mass vs ALamba mass for all pt for ALambda
  TH2F   *fHistPiPALambdaVsLambdaMass;                  // ALambda mass vs Lambda for Lambda
  TH2F   *fHistPiAPLambdaVsALambdaMass;                 // Lambda mass vs ALambda for ALambda

  //----------------------- K0 ----------------------------------------//
  TH1F   *fHistPiPiMass;                                // pi+pi- InvMass spectrum
  TH2F   *fHistPiPiMassVSPt;                            // pi+pi- InvMass spectrum vs pt
  TH2F   *fHistPiPiMassVSPtMCTruth;                     // pi+pi- InvMass spectrum vs pt MC truth 
  TH2F   *fHistPiPiMassVSPtPosMCTruth;                  // pi+ daughter InvMass spectrum vs pt MC truth
  TH2F   *fHistPiPiMassVSPtNegMCTruth;                  // pi- daughter InvMass spectrum vs pt MC truth
  TH2F   *fHistPiPiMassVSY;                             // pi+pi- InvMass spectrum vs rapidity
  TH2F   *fHistPiPiPtVSY;                               // pi+pi- pt vs rapidity

  // TH2F   *fHistPiPiMassVSAlpha;                        // pi+pi- InvMass spectrum vs armenteros alpha
  TH2F   *fHistPiPiRadiusXY;                            // pi+pi- opening angle vs mass
  TH2F   *fHistPiPiCosPointAng;                         // pi+pi- cosine of pointing angle vs pt or dca to vertex
  TH2F   *fHistPiPiDCADaughterPosToPrimVtxVSMass;       // dca of pos. K0s daughter to prim vtx vs mass
  TH2F   *fHistPiPiDecayLengthVsPt;                     // pi+pi- decay lenght vs pt
  TH2F   *fHistPiPiDecayLengthVsMass;                   // pi+pi- decay lenght vs pt
  TH2F   *fHistPiPiDecayLengthVsCtau;                   // pi+pi- decay lenght vs pt

  /*
  TH2F   *fHistPiPiDistDaughtersPos[3];                 // dist of pos daughter tracks for K0s cand. vs mass for 3 pt bins
  TH2F   *fHistPiPiDistDaughtersNeg[3];                 // dist of neg daughter tracks for K0s cand. vs mass for 3 pt bins
  TH2F   *fHistPiPiDCADaughtersPos[3];                  // weighted dca of pos daughter tracks for K0s cand. vs mass for 3 pt bins
  TH2F   *fHistPiPiDCADaughtersNeg[3];                  // weighted dca of neg daughter tracks for K0s cand. vs mass for 3 pt bins
  TH2F   *fHistPiPiRadAtDCA5cmDaughtersPos[3];          // radius in TPC at 5cm track distancd of pos daughter tracks for K0s cand. vs mass for 3 pt bins
  TH2F   *fHistPiPiRadAtDCA5cmDaughtersNeg[3];          // radius in TPC at 5cm track distancd of neg daughter tracks for K0s cand. vs mass for 3 pt bins
  */
  //TH2F   *fHistPiPiMassVSPtK0L;                       // K0L InvMass vs pt distribution
  TH2F   *fHistPiPiDCADaughters;                        // pi+pi- dca between daughters
  // TH2F   *fHistPiPiPtDaughters;                         // pi+pi- daughters pt pos vs pt neg 
  TH2F   *fHistPiPiDCAVSMass;                           // pi+pi- dca xy to prim vtx vs mass
  TH2F   *fHistPiPiDCAZVSMass;                           // pi+pi- dca z to prim vtx vs mass
  TH2F   *fHistPiPiDCAZPos;                             // dca z component of pos K0s daughter
  TH2F   *fHistPiPiDCAZNeg;                             // dca z component of neg K0s daughter
  TH2F   *fHistPiPiTrackLengthPosVsMass;                // track length of pos K0s daughter in TPC
  TH2F   *fHistPiPiTrackLengthNegVsMass;                // track length of neg K0s daughter in TPC  
  TH1F   *fHistPiPiMonitorCuts;                         // pi+pi- cut monitor
  TH1F   *fHistPiPiMonitorMCCuts;                       // pi+pi- cut monitor mc
  TH2F   *fHistPiPiDecayLengthResolution;               // pi+pi- decay length resolution: mcreco vs mctruth
  //detectors
  TH2F   *fHistNclsITSPosK0;                            // number of clusters from ITS of positive K0s daughters
  TH2F   *fHistNclsITSNegK0;                            // number of clusters from ITS of negative K0s daughters
  TH2F   *fHistNclsTPCPosK0;                            // number of clusters from TPC of positive K0s daughters
  TH2F   *fHistNclsTPCNegK0;                            // number of clusters from TPC of negative K0s daughters
  TH2F   *fHistChi2PerNclsITSPosK0;                     // chi^2 per number of clusters ITS of positive K0s daughters
  TH2F   *fHistChi2PerNclsITSNegK0;                     // chi^2 per number of clusters ITS of negative K0s daughters  
  TH2F   *fHistNCRowsTPCPosK0;                          // no of crossed rows for K0s pos daughter
  TH2F   *fHistNCRowsTPCNegK0;                          // no of crossed rows for K0s neg daughter
  TH2F   *fHistRatioFoundOverFinableTPCK0Pos;           // ratio of ncls findable over found TPC K0s daughters
  TH2F   *fHistRatioFoundOverFinableTPCK0Neg;           // ratio of ncls findable over found TPC K0s daughters
  // TH2F   *fHistPiPiDistDaughtersTPCEntrVsMass;
  //---Lambda Antilambda -----//
  //  TH2F* fHistPiPDistDaughtersTPCEntrVsMass;
  // TH2F* fHistPiAPDistDaughtersTPCEntrVsMass;
  //------------------------- MC only histos ---------------------------------------------------//
  TH2F   *fHistPrimVtxZESDVSNContributorsMC;        // count contributors to ESD vertex MC
  TH2F   *fHistPrimVtxZESDTPCVSNContributorsMC;     // count contributors to TPC vertex MC
  TH2F   *fHistPrimVtxZESDSPDVSNContributorsMC;     // count contributors to SPD vertex MC
  TH1F   *fHistMCVertexZ;                           // primary MC vertex position z 
  TH1F   *fHistPiPiPDGCode;                         // PDG code of K0 mothers
  TH1F   *fHistPiPPDGCode;                          // PDG code of Lambda mothers
  TH1F   *fHistPiAPPDGCode;                         // PDG code of Lambda mothers
  /*
  //-- BG of K0s
  TH2F   *fHistPiPiGA;
  TH2F   *fHistPiPiKch;
  TH2F   *fHistPiPiPhi;
  TH2F   *fHistPiPiL;
  TH2F   *fHistPiPiPi0;
  TH2F   *fHistPiPiPich;
  TH2F   *fHistPiPiRoh;
  TH2F   *fHistPiPiOmega;
  TH2F   *fHistPiPiKStar;
  TH2F   *fHistPiPiNoMother;
  TH2F *fHistPiPiK0s;
  TH2F *fHistPiPiK0L;
  TH2F *fHistPiPiN;
  TH2F *fHistPiPiSigma;
  TH2F *fHistPiPiXi;
  TH2F *fHistPiPiDelta;
  TH2F *fHistPiPiB;
  TH2F *fHistPiPiD;
  TH2F *fHistPiPiEta;
  //-- BG of Lambda
  TH2F   *fHistPiPGA;
  TH2F   *fHistPiPKch;
  TH2F   *fHistPiPK0s;
  TH2F   *fHistPiPPi0;
  TH2F   *fHistPiPPich;
  TH2F   *fHistPiPKStar;
  TH2F   *fHistPiPN;
  TH2F   *fHistPiPNoMother;
  TH2F   *fHistPiPL;
  */
  //others for (A)Lambda
  TH2F   *fHistPiPCosPointAngXiVsPt;                // cosine of pointing angle of xis vs pt
  TH2F   *fHistPiAPCosPointAngXiVsPt;               // cosine of pointing angle of xis vs pt
  TH2F   *fHistPiPMassVSPtSecXiMCTruth;
  TH2F   *fHistPiPMassVSPtSecOmegaMCTruth;
  TH2F   *fHistPiAPMassVSPtSecXiMCTruth;
  TH2F   *fHistPiAPMassVSPtSecOmegaMCTruth;

  //--------------------------------- histos with secondaries' histo------------------------------//
  TH2F   *fHistV0RadiusZ[2];                        // V0 decay radius z
  TH2F   *fHistV0RadiusZVSPt[2];                    // V0 decay radius z vs pt
  TH2F   *fHistV0RadiusXY[2];                       // V0 decay radius x vs y
  TH2F   *fHistV0RadiusXYVSY[2];                    // V0 decay radius xy vs rapidity
   
  TH2F   *fHistArmenteros[2];                       // armenteros

  //------------------------------------- Lambda -------------------------------------------------//
  TH1F   *fHistPiPMass[2];                          // p+pi- InvMass spectrum
  TH2F   *fHistPiPMassVSPt[2];                      // p+pi- InvMass spectrum vs pt
  TH2F   *fHistPiPMassVSPtMCTruth[2];               // p+pi- InvMass spectrum vs pt MC truth
  TH2F   *fHistPiPMassVSPtPosMCTruth[2];            // p+ daughter InvMass spectrum vs pt MC truth
  TH2F   *fHistPiPMassVSPtNegMCTruth[2];            // pi- daughter InvMass spectrum vs pt MC truth
  TH2F   *fHistPiPMassVSY[2];                       // p+pi- InvMass spectrum vs rapidity
  TH2F   *fHistPiPPtVSY[2];                         // p+pi- pt vs rapidity
  TH2F   *fHistPiPRadiusXY[2];                      // p+pi- opening angle vs mass
  TH2F   *fHistPiPCosPointAng[2];                   // p+pi- cosine of pointing angle vs pt  or dca to vertex
  TH2F   *fHistPiPDCADaughterPosToPrimVtxVSMass[2]; // dca of pos. Lambda daughter to prim vtx vs mass
  TH2F   *fHistPiPDCADaughterNegToPrimVtxVSMass[2]; // dca of neg. Lambda daughter to prim vtx vs mass
  TH2F   *fHistPiPDecayLengthVsPt[2];               // p+pi- decay lenght vs pt
  TH2F   *fHistPiPDecayLengthVsMass[2];             // p+pi- decay lenght vs pt
  TH2F   *fHistPiPDecayLengthVsCtau[2];             // p+pi- decay lenght vs pt
  /*
  TH2F   *fHistPiPDistDaughtersPos[3];              // dist of pos daughter tracks for Lambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiPDistDaughtersNeg[3];              // dist of neg daughter tracks for Lambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiPDCADaughtersPos[3];               // weighted dca of pos daughter tracks for Lambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiPDCADaughtersNeg[3];               // weighted dca of neg daughter tracks for Lambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiPRadAtDCA5cmDaughtersPos[3];       // radius in TPC at 5cm track distancd of pos daughter tracks for Lambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiPRadAtDCA5cmDaughtersNeg[3];       // radius in TPC at 5cm track distancd of neg daughter tracks for Lambda cand. vs mass for 3 pt bins
  */
  TH2F   *fHistPiPDCADaughters[2];                  // p+pi- dca between daughters
  //TH2F   *fHistPiPPtDaughters[2];                   // p+pi- daughters pt pos vs pt neg 
  TH2F   *fHistPiPDCAVSMass[2];                     // p+pi- dca xy to prim vtx vs mass
  TH2F   *fHistPiPDCAZVSMass[2];                     // p+pi- dca z to prim vtx vs mass
  TH1F   *fHistPiPMonitorCuts[2];                   // p+pi- cut monitor
  TH1F   *fHistPiPMonitorMCCuts[2];                 // p+pi- cut monitor mc
  TH2F   *fHistPiPMassVSPtSecSigma[2];              // InvMass distribution vs pt of secondary lambdas from sigma truth(0) reco(1)
  TH2F   *fHistPiPMassVSPtSecXi[2];                 // InvMass distribution vs pt of secondary lambdas from xi MC truth(0) reco(1)
  TH2F   *fHistPiPMassVSPtSecOmega[2];              // InvMass distribution vs pt of secondary lambdas from omega MC truth(0) reco(1)
  TH2F   *fHistPiPMassVSYSecXi[2];                  // InvMass distribution vs rapidity of secondary lambdas from xi MC truth(0) reco(1)
  TH2F   *fHistPiPXi0PtVSLambdaPt[2] ;              // pt of xi0 vs pt lambda truth(0) reco(1)
  TH2F   *fHistPiPXiMinusPtVSLambdaPt[2];           // pt of ximinus vs pt lambda truth(0) reco(1)
  TH2F   *fHistPiPOmegaPtVSLambdaPt[2];             // pt of omega plus vs pt alambda truth(0) reco(1)
  TH2F   *fHistPiPDecayLengthResolution[2];         // Lambda decay length resolution MCreco vs MC truth
  TH2F   *fHistPiPDCAZPos[2];                       // dca z component of pos Lambda daughter
  TH2F   *fHistPiPDCAZNeg[2];                       // dca z component of neg Lambda daughter
  TH2F   *fHistPiPTrackLengthPosVsMass[2];          // track length of pos Lambda daughter in TPC
  TH2F   *fHistPiPTrackLengthNegVsMass[2];          // track length of neg Lambda daughter in TPC

  //---------------------------------------- Antilambda --------------------------------------------------------------//
  TH1F   *fHistPiAPMass[2];                         // pi+p- InvMass spectrum
  TH2F   *fHistPiAPMassVSPt[2];                     // pi+p- InvMass spectrum vs pt
  TH2F   *fHistPiAPMassVSPtMCTruth[2];              // pi+p- InvMass spectrum vs pt MC Truth
  TH2F   *fHistPiAPMassVSPtPosMCTruth[2];           // pi+ daughter InvMass spectrum vs pt MC Truth
  TH2F   *fHistPiAPMassVSPtNegMCTruth[2];           // p- daughter InvMass spectrum vs pt MC Truth
  TH2F   *fHistPiAPMassVSY[2];                      // pi+p- InvMass spectrum vs rapidity
  TH2F   *fHistPiAPPtVSY[2];                        // pi+p- pt vs rapidity
  TH2F   *fHistPiAPRadiusXY[2];                     // pi+p- opening angle vs mass
  TH2F   *fHistPiAPCosPointAng[2];                  // pi+p- cosine of pointing angle vs pt  or dca to vertex
  TH2F   *fHistPiAPDCADaughterPosToPrimVtxVSMass[2];// dca of pos ALambda daughter to prim vtx vs mass
  TH2F   *fHistPiAPDCADaughterNegToPrimVtxVSMass[2];// dca of neg ALambda daughter to prim vtx vs mass
  TH2F   *fHistPiAPDecayLengthVsPt[2];              // pi+p- decay lenght vs pt
  TH2F   *fHistPiAPDecayLengthVsMass[2];            // pi+p- decay lenght vs pt
  TH2F   *fHistPiAPDecayLengthVsCtau[2];            // pi+p- decay lenght vs pt
 
  /*
  TH2F   *fHistPiAPDistDaughtersPos[3];              // dist of pos daughter tracks for ALambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiAPDistDaughtersNeg[3];              // dist of neg daughter tracks for ALambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiAPDCADaughtersPos[3];               // weighted dca of pos daughter tracks for ALambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiAPDCADaughtersNeg[3];               // weighted dca of neg daughter tracks for aLambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiAPRadAtDCA5cmDaughtersPos[3];       // radius in TPC at 5cm track distancd of pos daughter tracks for ALambda cand. vs mass for 3 pt bins
  TH2F   *fHistPiAPRadAtDCA5cmDaughtersNeg[3];       // radius in TPC at 5cm track distancd of neg daughter tracks for ALambda cand. vs mass for 3 pt bins
  */
  TH2F   *fHistPiAPDCADaughters[2];                 // pi+p- dca between daughters
  // TH2F   *fHistPiAPPtDaughters[2];                  // pi+p- daughters pt pos vs pt neg 
  TH2F   *fHistPiAPDCAVSMass[2];                    // pi+p- dca xy prim vtx vs mass
  TH2F   *fHistPiAPDCAZVSMass[2];                    // pi+p- dca z  prim vtx vs mass
  TH1F   *fHistPiAPMonitorCuts[2];                  // pi+p- cut monitor
  TH1F   *fHistPiAPMonitorMCCuts[2];                // pi+p- cut monitor mc
  TH2F   *fHistPiAPMassVSPtSecSigma[2];             // InvMass distribution vs pt of secondary alambdas from sigma truth(0) reco(1)
  TH2F   *fHistPiAPMassVSPtSecXi[2];                // InvMass distribution vs pt of secondary alambdas from xi MC truth(0) reco(1)
  TH2F   *fHistPiAPMassVSPtSecOmega[2];             // InvMass distribution vs pt of secondary alambdas from omega MC truth(0) reco(1)
  TH2F   *fHistPiAPMassVSYSecXi[2];                 // InvMass distribution vs rapidity of secondary alambdas from xi MC truth(0) reco(1)
  TH2F   *fHistPiAPXi0PtVSLambdaPt[2] ;             // pt of xi0 vs pt alambda truth(0) reco(1)
  TH2F   *fHistPiAPXiMinusPtVSLambdaPt[2];          // pt of ximinus vs pt alambda truth(0) reco(1)
  TH2F   *fHistPiAPOmegaPtVSLambdaPt[2];            // pt of omega plus vs pt alambda truth(0) reco(1)
  TH2F   *fHistPiAPDecayLengthResolution[2];        // ALambda decay length resolution MCreco vs MC truth
  //  TH2F   *fHistPiAPDCAZPos[2];                      // dca z component of pos ALambda daughter
  //TH2F   *fHistPiAPDCAZNeg[2];                      // dca z component of neg ALambda daughter
  TH2F   *fHistPiAPTrackLengthPosVsMass[2];         // track length of pos ALambda daughter in TPC
  TH2F   *fHistPiAPTrackLengthNegVsMass[2];         // track length of neg ALambda daughter in TPC


  //-------------------------------------------------------- others --------------------------------------------------//
  //dEdx
  TH2F   *fHistDedxSecProt[2];                      // dedx from proton cadidates vs pt
  TH2F   *fHistDedxSecAProt[2];                     // dedx from antiproton candidates vs pt
  TH2F   *fHistDedxSecPiMinus[2];                   // dedx from pi minus candidates vs pt
  TH2F   *fHistDedxSecPiPlus[2];                    // dedx from pi plus candidates vs pt
  TH2F   *fHistDedxProt[2];                         // dedx from proton cadidates vs pt before pidcut
  TH2F   *fHistDedxAProt[2];                        // dedx from antiproton candidates vs pt before pidcut
  TH2F   *fHistDedxPiMinus[2];                      // dedx from pi minus candidates vs pt before pidcut
  TH2F   *fHistDedxPiPlus[2];                       // dedx from pi plus candidates vs pt before pidcut
   
  //clusters
  TH2F   *fHistNclsITS[2];                          // number of clusters ITS pos vs neg daughters
  TH2F   *fHistNclsTPC[2];                          // number of clusters TPC  neg daughters vs number of crossed rows
  TH2F   *fHistNclsITSPosL[2];                      // number of clusters from ITS of positive lambda daughters
  TH2F   *fHistNclsITSNegL[2];                      // number of clusters from ITS of negative lambda daughters
  TH2F   *fHistNclsTPCPosL[2];                      // number of clusters from TPC of positive lambda daughters
  TH2F   *fHistNclsTPCNegL[2];                      // number of clusters from TPC of negative lambda daughters
  TH2F   *fHistChi2PerNclsITSPosL[2];               // chi^2 per number of clusters ITS of positive lambda daughters
  TH2F   *fHistChi2PerNclsITSNegL[2];               // chi^2 per number of clusters ITS of negative lambda daughters
  TH2F   *fHistNCRowsTPCPosL[2];                    // number of crossed rows for Lambda pos daughter
  TH2F   *fHistNCRowsTPCNegL[2];                    // number of crossed rows for Lambda neg daughter
  TH2F   *fHistRatioFoundOverFinableTPCLPos[2];     // ratio of ncls findable over found TPC L daughters
  TH2F   *fHistRatioFoundOverFinableTPCLNeg[2];     // ratio of ncls findable over found TPC L daughters
  TH2F   *fHistPiPiEtaDMC[2];                       // eta of daughters vs pt K0s MC truth raw(0) after cuts(1)
  TH2F   *fHistPiPEtaDMC[2];                        // eta of daughters vs pt lambda MC truth raw(0) after cuts(1)
  TH2F   *fHistPiPiEtaDReco[2];                     // eta of daughters ESD track vs eta AliESDv0 or vs pt K0s raw(0) after cuts(1)
  TH2F   *fHistPiPEtaDReco[2];                      // eta of daughters ESD track vs eta AliESDv0 or vs  pt (a)lambda raw(0) after cuts(1)

  /*
  //user shift
  TH1F   *fHistUserPtShift;//monitor user defined charge/pt shift
  */


   
  //---------------------------------- Variables--------------------------------------------//

  //--cut options --//
  //MC only
  Bool_t    fMCMode;                   // run over MC general yes/no
  Bool_t    fMCTruthMode;              // MC truth selection yes/no
  Bool_t    fSelectInjected;           // for MC with injected signals, select only injected
  Bool_t    fSelectMBMotherMC;         // for MC with injected signals, select only MB MC mother for sec. Lambdas
  Bool_t    fCheckNegLabelReco;        // reject MC truth and reco for neg labels in reco
  Bool_t    fOnlyFoundRecoV0;          // reject MC truth if reco V0 not found

  // Calculate centrality
  Int_t     fUseCentrality;            // use centrality (0=off(default),1=VZERO,2=SPD)
  Int_t     fUseCentralityBin;         // centrality bin to be used 
  Int_t     fUseCentralityRange;       // use centrality (0=off(default),1=VZERO,2=SPD) 

  //pp analysis
  Bool_t    fAnapp;                    // flag for pp analysis
  Bool_t    fRejectPileUpSPD;          // reject pileup events from SPD 
  Bool_t    fSelSDD;                   // select pp events with SDD (for pp 2.76TeV LHC11a)
  Bool_t    fSelNoSDD;                 // select pp events with no SDD (for pp 2.76TeV LHC11a)
  //onthefly
  Bool_t    fOntheFly;                 // true if onfly finder shall be used

  //vertex
  Double_t  fVertexZCut;               // z vertex cut value
  Bool_t    fVtxStatus;                // vertex cut on/off

  //esdtrackcuts
  Int_t     fNcr;                      // esd track cuts: number of crossed rows TPC for V0 daughters
  Double_t  fChi2cls;                  // esd track cuts: chi2 per cluster TPC for V0 daughters
  Bool_t    fTPCrefit;                 // esd track cuts: tpc refit  for V0 daughters
  Bool_t    fITSrefit;                 // esd track cuts: its refit  for V0 daughters for offline finder

  Int_t     fNcrCh;                    // esd track cuts: number of crossed rows TPC for charged
  Double_t  fChi2clsCh;                // esd track cuts: chi2 per cluster TPC for charged
  Bool_t    fTPCrefitCh;               // esd track cuts: tpc refit for charged
  Bool_t    fITSrefitCh;               // esd track cuts: its refit  for V0 daughters for offline finder

  Int_t     fNcrLpt;                   // esd track cuts: number of crossed rows TPC for low pt
  Double_t  fChi2clsLpt;               // esd track cuts: chi2 per cluster TPC for low pt
  Bool_t    fTPCrefitLpt;              // esd track cuts: tpc refit for low pt
  Bool_t    fITSrefitLpt;              // esd track cuts: its refit  for V0 daughters for offline finder

  //PID
  Bool_t    fUsePID;                   // use proton pid yes/no
  Bool_t    fUsePIDPion;               // use pion pid yes/no
  Double_t  fNSigma;                   // set nsigma value
  Double_t  fNSigma2;                  // set nsigma 2 value
  Double_t  fPPIDcut;                  // set max momentum for pid cut usage 
  Double_t  fPtTPCCut;                 // low pt limit cut for TPC cluster cuts from AliESDtrackCuts
  Bool_t    fMoreNclsThanRows;         // cut on ncls>ncrossed rows yes/no
  Bool_t    fMoreNclsThanFindable;     // cut on ncls>nfindable cls yes/no
  Bool_t    fMoreNclsThanFindableMax;  // cut on ncls>nfindable max cls yes/no
  Double_t  fRatioFoundOverFindable;   // cut on found over findable clusters TPC
  Double_t  fRatioMaxCRowsOverFindable;// cut on crossed rows over finable max
  Double_t  fChi2PerClusterITS;        // cut on chi2 per ITS cluster
  Double_t  fDistanceTPCInner;         // cut on distance of daughters at TPC entrance
  Int_t     fMinNCLSITSPos;            // min ncls ITS of pos daugter cut
  Int_t     fMinNCLSITSNeg;            // min ncls ITS of neg daugter cut
  Int_t     fMaxNCLSITSPos;            // max ncls ITS of pos daugter cut
  Int_t     fMaxNCLSITSNeg;            // max ncls ITS of neg daugter cut
  Bool_t    fSwitchCaseITSCls;         // apply pos and neg ITS cls cluster cut with 
                                       // or for both daughters for at least one of the daughters shall have ...
  Bool_t    fCutMITrackLength;         // cut on geom track length in TPC as Marian Ivanov sugg.
  Bool_t    fCutMICrossedR;            // cut on crossed rows in TPC as Marian Ivanov sugg.
  Bool_t    fCutMITPCncls;             // cut on ncls in TPC as Marian Ivanov sugg.
  Double_t  fCutMITrackLengthLengthF;  // cut on track length in TPC as Marian Ivanov sugg. length factor
  Double_t  fCutMICrossedRLengthF;     // cut on crossed rows in TPC as Marian Ivanov sugg. length factor

  //rapidity
  Bool_t    fRapCutV0;                 // use rapidity cut for V0 yes/no
  Double_t  fRap;                      // user defined value for rapidity cut

  //eta and pt
  Double_t  fEtaCutMCDaughters;        // eta cut for MC daughters on/off
  Double_t  fEtaCutMCDaughtersVal;     // eta cut value for MC daughters
  Double_t  fEtaSignCut;            // eta cutsign daughters

  //for FD and sec. variables calc
  Bool_t    fUseXi0;                   // take xi0 into account
  Bool_t    fUseXiM;                   // take xi minus/plus into account
  Bool_t    fUseOmega;                 // take omega into account
  Bool_t    fCutRapXi;                 // do rapidity cut for Xi
  Double_t  fMinPt;                    // pt min cut value 

  //armenteros
  Double_t  fAlfaCut;                  // set alpha armenteros cut value
  Double_t  fQtCut;                    // set ptmax for qt armenteros cut 
  Double_t  fQtCutPt;                  // set ptmax for  qt armenteros cut
  Double_t  fQtCutPtLow;               // set ptmin for  qt armenteros cut
  Bool_t    fArmCutK0;                 // set armenteros cut on/off for K0s
  Bool_t    fArmCutL;                  // set armenteros cut on/off for Lambda
  Double_t  fArmQtSlope;               // slope for armenteros K0s cut: qt = alpha*slope
  //others
  Double_t  fExcludeLambdaFromK0s;     // exlude Lambda mass from K0s throuh mass difference below this value
  Double_t  fExcludeK0sFromLambda;     // exlude K0s mass from Lambda throuh mass difference below this value
  Double_t  fExcludePhotonsFromK0s;    // exlude photons from K0s throuh mass difference below this value
  Double_t  fExcludePhotonsFromLambda; // exlude photons from K0s throuh mass difference below this value
  Double_t  fDCAToVertexK0;            // dca of V0 to vertex cut value K0s
  Double_t  fDCAToVertexL;             // dca of V0 to vertex cut value L/AL
  Double_t  fDCAXK;                    // dca in x of K0s to vertex cut value
  Double_t  fDCAYK;                    // dca in y of K0s to vertex cut value
  Double_t  fDCAXL;                    // dca in x of Lambda to vertex cut value
  Double_t  fDCAYL;                    // dca in y of Lambda to vertex cut value
  Double_t  fDCAZ;                     // dca in z of V0 to vertex cut value
  
  Double_t  fDCADaughtersL;            // dca between Lambda daughters cut value
  Double_t  fDCADaughtersAL;           // dca between ALambda daughters cut value
  Double_t  fDCADaughtersK0;           // dca between K0s daughters cut value
  
  Double_t  fDCADaughtersToVtxLarge;   // dca large between V0 daughters and vertex cut value
  Double_t  fDCADaughtersToVtxSmall;   // dca small between V0 daughters and vertex cut value
  
  Double_t  fDecayRadXYMin;            // minmal decay radius in x-y cut value
  Double_t  fDecayRadXYMax;            // maximal decay radius in x-y cut value
  Double_t  fPtDecRadMin;              // pt cut for max pt of radius cut usage
  Double_t  fCosPointAngL;             // cosine of pointing angle cut value for Lambda and ALambda
  Double_t  fCosPointAngK;             // cosine of pointing angle cut value for K0s
  Double_t  fCPAPtCutK0;               // pt max for cosine of pointing angle cut K0s
  Double_t  fCPAPtCutL;                // pt max for cosine of pointing angle cut Lambda
  Double_t  fOpengAngleDaughters;      // cut on opening angle between V0 daughters
  Double_t  fOpAngPtCut;               // max pt for using the  opening angle between V0 daughters cut
    
  Double_t  fDecayLengthMax;           // maximal decay length in x-y-z cut value
  Double_t  fDecayLengthMin;           // minimal decay length in x-y-z cut value

  Double_t  fDecRadCutITSMin;          // radius min for ITS cluster cut
  Double_t  fDecRadCutITSMax;          // radius max for ITS cluster cut

  //ctau
  Double_t  fCtauK0s;                  // multiple of ctau cut value for K0s
  Double_t  fCtauL;                    // multiple of ctau cut value for Lambda
  Double_t  fCtauPtCutK0;              // pt max for ctau cut usage for K0s
  Double_t  fCtauPtCutL;               // pt max for ctau cut usage for Lambda

  //KF particle chi cut
  Double_t  fChiCutKf;            //cut value of chi2 of AliKFParticle
    // Bool_t    fChiCutKf;                 //cut value of chi2 of AliKFParticle

  Double_t  fK0sLowMassCut;            //lower cut on K0s mass
  Double_t  fK0sHighMassCut;           //higher cut on K0s mass

  Double_t  fLLowMassCut;              //lower cut on Lambda mass
  Double_t  fLHighMassCut;             //higher cut on lambda mass


  Bool_t   fSetFillDetAL;              // fill det histo with AL instead of Lambda
  Bool_t   fSetPtDepHist;              // fill pt instead of mass

  Bool_t   fStopLoop;                  // set stop reco loop to reject multiple times found V0s

  Double_t fDistDauForCheck;            // upper limit for distance of daughters check

  
  // option for user defined charge/pt shift
  Bool_t     fShift;// shift yes/no
  Double_t   fDeltaInvP;//define shift value
  
 

  AliAnalysisTaskV0ForRAA(const AliAnalysisTaskV0ForRAA&);
  AliAnalysisTaskV0ForRAA&operator=(const AliAnalysisTaskV0ForRAA&);
   
  ClassDef(AliAnalysisTaskV0ForRAA, 1); 
};
#endif
