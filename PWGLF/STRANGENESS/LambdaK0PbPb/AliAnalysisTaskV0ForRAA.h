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

class Tlist;

class AliESDv0;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpid;
class AliESDEvent;
class AliMCEvent;
class AliPIDResponse;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif


class AliAnalysisTaskV0ForRAA : public AliAnalysisTaskSE {
 public:
   
  AliAnalysisTaskV0ForRAA(const char *name="AliAnalysisTaskV0ForRAA");
  virtual ~AliAnalysisTaskV0ForRAA();
  

  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);



  //-- MC truth/reco --//
  void SetMCMode(Bool_t mcmode)                               {fMCMode = mcmode; if(fMCMode) Printf("AliAnalysisTaskV0ForRAA::running mc mode: histos of MC reco");}
  void SetMCTruthMode(Bool_t mcmode)                          {fMCTruthMode = mcmode; if(fMCTruthMode) Printf("AliAnalysisTaskV0ForRAA::running mc mode: histos of MC truth");}
  void SelectInjected(Bool_t injected)                        {fSelectInjected = injected; if(fSelectInjected) Printf("AliAnalysisTaskV0ForRAA::only injected MC particles");}
  void SelectMBMotherMC(Bool_t mbmother)                      {fSelectMBMotherMC = mbmother;if(mbmother)  Printf("AliAnalysisTaskV0ForRAA::only MB mother MC for sec lambdas selected");}
  void SelectOnlyPosLabelMC(Bool_t poslabel)                  {fCheckNegLabelReco = poslabel;if(poslabel) Printf("AliAnalysisTaskV0ForRAA::Select only MC truth and reco with pos label reco");}

  void SelectOnlyFoundRecoV0MC(Bool_t found)                  {fOnlyFoundRecoV0 = found;if(found) Printf("AliAnalysisTaskV0ForRAA::Select only MC truth with found reco V0");}


  //-- Centrality  --//
  // use centrality - if yes, which one
  void  SetUseCentrality(Int_t cent)                          {fUseCentrality = cent; Printf("AliAnalysisTaskV0ForRAA::centrality selected for detector %i (0=off, 1=VZERO, 2=SPD)",cent);}
  // set range
  void  SetUseCentralityRange(Int_t range)                    {fUseCentralityRange = range;if(fUseCentrality) Printf("AliAnalysisTaskV0::centrality range %i",fUseCentralityRange);}
  // centrality bin to be used
  void  SetUseCentralityBin(Int_t bin)                        {fUseCentralityBin = bin; if(fUseCentrality) Printf("AliAnalysisTaskV0ForRAA::centrality selected for bin %i",fUseCentralityBin); }


  //-- event cuts --//
  void SetPrimVertexZCut(Double_t vtxcut,Bool_t status)       {fVertexZCut = vtxcut;fVtxStatus = status; Printf("AliAnalysisTaskV0ForRAA::SetPrimVertexZCut %3.2f",vtxcut);}
  void SetAnapp(Bool_t anapp)                                 {fAnapp = anapp ;if(fAnapp) Printf("AliAnalysisTaskV0ForRAA::analysing pp!!!");}
  void SelectWithSDD(Bool_t sdd)                              {fSelSDD =sdd; if(sdd) Printf("AliAnalysisTaskV0ForRAA:: only events with SDD selected!");}
  void SelectWithNoSDD(Bool_t sdd)                              {fSelNoSDD =sdd; if(sdd) Printf("AliAnalysisTaskV0ForRAA:: only events with NO SDD selected!");}

  //-- track cuts --//
  void SetESDTrackCuts(AliESDtrackCuts *esdcuts =NULL)        {fESDTrackCuts = esdcuts;Printf("AliAnalysisTaskV0ForRAA::AliESDtrackCuts for V0s set");}
  void SetESDTrackCutsCharged(AliESDtrackCuts *esdcuts=NULL)  {fESDTrackCutsCharged = esdcuts;Printf("AliAnalysisTaskV0ForRAA::AliESDtrackCuts for charged particles set");}
  void SetESDTrackCutsLowPt(AliESDtrackCuts *esdcuts=NULL)    {fESDTrackCutsLowPt = esdcuts;Printf("AliAnalysisTaskV0ForRAA::AliESDtrackCuts for low pt particles set");}
  void SetUseOnthefly(Bool_t useonthefly)                     {fOntheFly = useonthefly; if(!fOntheFly) Printf("AliAnalysisTaskV0ForRAA::offline V0s");}
  void SetUsePID(Bool_t usepid,Double_t nsigma=100.0,Double_t pcut=100.0) {fUsePID = usepid;fNSigma = nsigma;fPPIDcut = pcut; if(fUsePID) Printf("AliAnalysisTaskV0ForRAA::PID! of %4.2f for p: %4.2f",fNSigma,pcut);}
  void SetCutMoreNclsThanRows(Bool_t cut)                     {fMoreNclsThanRows=cut; if(cut) Printf("AliAnalysisTaskV0ForRAA::cut on more ncls than crossed rows");}  
  void SetCutMoreNclsThanFindable(Bool_t cut)                 {fMoreNclsThanFindable=cut; if(cut) Printf("AliAnalysisTaskV0ForRAA::cut on more ncls than ncls findable");}
  void SetCutMoreNclsThanFindableMax(Bool_t cut)              {fMoreNclsThanFindableMax = cut; if(cut) Printf("AliAnalysisTaskV0ForRAA::cut on more ncls than ncls findable max");}

  void SetRatioFoundOverFindable(Double_t cut)                 {fRatioFoundOverFindable = cut; Printf("AliAnalysisTaskV0ForRAA::cut on found over finable clusters %f",cut);}
  void SetRatioMaxCRowsOverFindable(Double_t cut)              {fRatioMaxCRowsOverFindable = cut;  Printf("AliAnalysisTaskV0ForRAA::cut on max crossed rows over finable clusters %f",cut);}

  void SetLowPtTPCCutAliESDTrackCut(Double_t pt)              {fPtTPCCut=pt;Printf("AliAnalysisTaskV0ForRAA::SetLowPtTPCCutAliESDTrackCut pt=%2.2f",pt);} 
   
  void SetMaxChi2PerITSCluster(Double_t chi2)                 {fChi2PerClusterITS = chi2; Printf("AliAnalysisTaskV0ForRAA::max chi2 per ITS cluster %3.2f",chi2);}
  void SetRapidityCutMother(Bool_t cut,Double_t val=5.0)      {fRapCutV0 = cut; fRap = val; if(cut) Printf("AliAnalysisTaskV0ForRAA::cut on mother rapidity %2.2f",val);}
  void SetMinPt(Double_t minPt=0.0)                           {fMinPt = minPt; if(minPt>0.0) Printf("AliAnalysisTaskV0ForRAA::cut on min pt %2.2f",minPt);}
  /*  void SetPtShift(const Double_t shiftVal) {
  //user defined shift in charge/pt
  if(shiftVal) { fShift=kTRUE; fDeltaInvP = shiftVal; Printf("AliAnalysisTaskV0::WARNING!!!!!!!!!!!!!! pt shift introduced!");}
  }
  */
  void SetDCAV0ToVertexK0(Double_t dcaTovertex)               {fDCAToVertexK0 = dcaTovertex; Printf("AliAnalysisTaskV0ForRAA::dca to vertex K0s %2.3f",dcaTovertex);}
  void SetDCAV0ToVertexL(Double_t dcaTovertex)                {fDCAToVertexL = dcaTovertex; Printf("AliAnalysisTaskV0ForRAA::dca to vertex L/AL %2.3f",dcaTovertex);}
  void SetDCADaughtersL(Double_t dcaDaughters)                {fDCADaughtersL = dcaDaughters; Printf("AliAnalysisTaskV0:ForRAA:dca daughters L %2.3f",dcaDaughters);}
  void SetDCADaughtersAL(Double_t dcaDaughters)               {fDCADaughtersAL = dcaDaughters; Printf("AliAnalysisTaskV0ForRAA::dca daughters AL %2.3f",dcaDaughters);}
  void SetDCADaughtersK0(Double_t dcaDaughters)               {fDCADaughtersK0 = dcaDaughters; Printf("AliAnalysisTaskV0ForRAA::dca daughters K0s %2.3f",dcaDaughters);}
  void SetDCADaughtersLargeToVertex(Double_t dcaDaughtersVtx) {fDCADaughtersToVtxLarge = dcaDaughtersVtx; Printf("AliAnalysisTaskV0ForRAA::dca daughters to vertex large %2.3f",dcaDaughtersVtx);}
  void SetDCADaughtersSmallToVertex(Double_t dcaDaughtersVtx) {fDCADaughtersToVtxSmall = dcaDaughtersVtx; Printf("AliAnalysisTaskV0ForRAA::dca daughters to vertex small %2.3f",dcaDaughtersVtx);}
  void SetDecayRadiusXYMinMax(Double_t decMin,Double_t decMax){fDecayRadXYMin = decMin;fDecayRadXYMax = decMax; Printf("AliAnalysisTaskV0ForRAA::min xy decay radius %2.3f max %2.3f",decMin,decMax);}
  void SetCosOfPointingAngleL(Double_t pointAng,Double_t ptMaxCut=100.0) {fCosPointAngL=pointAng;fCPAPtCutL = ptMaxCut;Printf("AliAnalysisTaskV0ForRAA::SetCosOfPointingAngleL %1.5f and pt max %2.2f",pointAng,ptMaxCut);} 
  void SetCosOfPointingAngleK(Double_t pointAng,Double_t ptMaxCut=100.0) {fCosPointAngK=pointAng;fCPAPtCutK0 = ptMaxCut; Printf("AliAnalysisTaskV0ForRAA::SetCosOfPointingAngleK  %1.5f and pt max %2.2f",pointAng,ptMaxCut);}
  void SetOpeningAngleCut(Double_t opang, Double_t maxpt)     {fOpengAngleDaughters=opang; fOpAngPtCut=maxpt,Printf("AliAnalysisTaskV0::cut on opening angle %1.3f up to pt= %2.2f",opang,maxpt);}

  void SetMaxDecayLength(Double_t decLength)                  {fDecayLengthMax = decLength; Printf("AliAnalysisTaskV0ForRAA::SetMaxDecayLength %2.3f",decLength);}
  void SetMinDecayLength(Double_t decLength)                  {fDecayLengthMin = decLength; Printf("AliAnalysisTaskV0ForRAA::SetMinDecayLength %2.3f",decLength);}
  void SetDCAXK0(Double_t dcaXK)                              {fDCAXK = dcaXK; Printf("AliAnalysisTaskV0ForRAA::SetDCAXK0 %2.3f",dcaXK);}
  void SetDCAYK0(Double_t dcaYK)                              {fDCAYK = dcaYK; Printf("AliAnalysisTaskV0ForRAA::SetDCAYK0 %2.3f",dcaYK);}
  void SetDCAXLambda(Double_t dcaXL)                          {fDCAXL = dcaXL; Printf("AliAnalysisTaskV0ForRAA::SetDCAXLambda %2.3f",dcaXL);}
  void SetDCAYLambda(Double_t dcaYL)                          {fDCAXL = dcaYL; Printf("AliAnalysisTaskV0ForRAA::SetDCAYLambda %2.3f",dcaYL);}
  void SetDCAZ(Double_t dcaZ)                                 {fDCAZ = dcaZ; Printf("AliAnalysisTaskV0ForRAA::SetDCAZ %2.3f",dcaZ);}
  void SetChi2CutKf(Bool_t chi2){ fChiCutKf = chi2; Printf("AliAnalysisTaskV0ForRAA::SetChi2CutKf %i",chi2);}
  //Double_t chi2)                            {fChiCutKf = chi2; Printf("AliAnalysisTaskV0ForRAA::SetChi2CutKf %3.2f",chi2);}
  void SetArmenterosCutAlpha(Double_t alfaMin)                {fAlfaCut=alfaMin;Printf("AliAnalysisTaskV0ForRAA::SetArmenterosCut a=%1.3f",alfaMin);}
  void SetArmenterosCutQt(Double_t ptmin,Double_t ptmax,Bool_t k0s,Bool_t la,Double_t slope=0.2){fQtCut = ptmax;fQtCutPtLow=ptmin, fArmQtSlope=slope,fArmCutK0=k0s;fArmCutL=la;Printf("AliAnalysisTaskV0ForRAA::SetArmenterosCut ptmin = %3.2f ptmax = %3.2f. slope: %1.2f.  Is K0s? %i La? %i",ptmin,ptmax,slope,k0s,la);}
  void SetMinMassDiffLK0s(Double_t diffK,Double_t diffL)                {fExcludeLambdaFromK0s = diffK;fExcludeK0sFromLambda = diffL; Printf("AliAnalysisTaskV0ForRAA::SetMaxMassDifferenceL for K0s %1.3f  K0s for L %1.3f",diffK,diffL);}

  void SetCtauCut(Double_t ctK0s, Double_t ctL,Double_t ptK0=100.0,Double_t ptL=100.0) {fCtauK0s = ctK0s*2.6842; fCtauL = ctL*7.89;fCtauPtCutK0=ptK0; fCtauPtCutL=ptL;
    Printf("AliAnalysisTaskV0ForRAA::SetCtauCut ctK=%2.2f, ctL = %2.2f for ptK= %5.2f ptL=%5.2f",ctK0s,ctL,ptK0,ptL);}
  void SetDoEtaOfMCDaughtersCut(Bool_t doCut,Double_t eta=5.0){fEtaCutMCDaughters =doCut; fEtaCutMCDaughtersVal=eta; Printf("AliAnalysisTaskV0ForRAA::eta cut on V0 (MC truth ? %i) daughters %1.3f !",doCut,eta);}
  //  void SetEtaSignCut(Double_t etasign)                        {fEtaSignCut = etasign;Printf("AliAnalysisTaskV0ForRAA::eta cut sign on  daughters %2.2f !",etasign);}

  
 private:
   
  //----------------------------functions --------------------------------------------//

  void   Process();                                                                                                   // process event
  void   V0RecoLoop(Int_t id0,Int_t id1,Int_t isSecd,Int_t what,Double_t ptV0MC,Int_t pdgMother,Double_t ptXiMother,Double_t decaylengthMCV0); // loop over reconstructed V0 (data or MC)
  void   V0MCTruthLoop();                                                                                             // loop over MC truth V0s
  Int_t  CalculateCentralityBin();                                                                                    // get the centrality bin from multiplicity
  Bool_t GetMCTruthPartner(AliESDtrack *pos,AliESDtrack *neg,Int_t id0,Int_t id1);                                    // find MC truth partner for reconstructed track


   
  //----------------------------- objects ----------------------------------------------//

  //event
  AliESDEvent     *fESD;                 // ESD event object
  AliMCEvent      *fMCev;                // MC event object

  
  //PID and track cuts
  AliPIDResponse  *fESDpid;              // pid object
  AliESDtrackCuts *fESDTrackCuts;        // esd track cuts for daughters
  AliESDtrackCuts *fESDTrackCutsCharged; // esd track cuts for all charged particles
  AliESDtrackCuts *fESDTrackCutsLowPt;   // esd track cuts for daughters at low pt

  TList           *fOutputContainer;     // output data container
   


  //----------------------------histograms --------------------------------------------//

  //---------------------------event histos --------------------------//
  TH1F   *fHistITSLayerHits;                        // pp 2.76 TeV analysis: check hist on div. ITS layer
  TH1F   *fHistOneHitWithSDD;                       // pp 2.76 TeV analysis: check hist on at least one ITS layer
  TH1F   *fHistNEvents;                             // count number of events for each event cut
  TH2F   *fHistPrimVtxZESDVSNContributors;          // count contributors to ESD vertex
  TH2F   *fHistPrimVtxZESDTPCVSNContributors;       // count contributors to TPC vertex
  TH2F   *fHistPrimVtxZESDSPDVSNContributors;       // count contributors to SPD vertex

  TH2F   *fHistPrimVtxZESDVSNContributorsMC;        // count contributors to ESD vertex MC
  TH2F   *fHistPrimVtxZESDTPCVSNContributorsMC;     // count contributors to TPC vertex MC
  TH2F   *fHistPrimVtxZESDSPDVSNContributorsMC;     // count contributors to SPD vertex MC

  TH1F   *fHistPrimVtxZESD;                         // primary ESD vertex position z after cuts and processing
  TH1F   *fHistPrimVtxZESDTPC;                      // primary TPC vertex position z after cuts and processing
  TH1F   *fHistPrimVtxZESDSPD;                      // primary SPD vertex position z after cuts and processing

  TH1F   *fHistESDVertexZ;                          // primary TPC vertex position z before cuts
  TH1F   *fHistMCVertexZ;                           // primary MC vertex position z 
   
  TH1F   *fHistMuliplicity;                         // number of particles from centrality selection
  TH1F   *fHistMuliplicityRaw;                      // number of particles from centrality selection before processing
  TH1F   *fHistCentBinRaw;                          // events per centralitybin before centrality selection
  TH1F   *fHistCentBin;                             // events per centralitybin
  TH1F   *fHistMultiplicityPrimary;                 // number of charged particles
   
  TH1F   *fHistNPrim;                               // number of contributors to the prim vertex


  //------------------------ single V0 histos MC case--------------------------//
  //K0s
  TH1F   *fHistPiPiPDGCode;                         // PDG code of K0 mothers
 
  //Lambda
  TH1F   *fHistPiPPDGCode;                          // PDG code of Lambda mothers 
  TH2F   *fHistPiPCosPointAngXiVsPt;                // p+pi- cosine of pointing angle of xis vs pt
  TH2F   *fHistPiPMassVSPtSecXiMCTruth;             // p+pi- InvMass spectrum vs Xi (-,0) pt MC truth
  TH2F   *fHistPiPMassVSPtSecOmegaMCTruth;          // p+pi- InvMass spectrum vs Omega (-)  pt MC truth
  
  //AntiLambda
  TH1F   *fHistPiAPPDGCode;                         // PDG code of AntiLambda mothers
  TH2F   *fHistPiAPCosPointAngXiVsPt;               // p-pi+ cosine of pointing angle of xis vs pt
  TH2F   *fHistPiAPMassVSPtSecXiMCTruth;            // p-pi+ InvMass spectrum vs Xi (+,anti 0) pt MC truth
  TH2F   *fHistPiAPMassVSPtSecOmegaMCTruth;         // p-pi+ InvMass spectrum vs Omega (+)  pt MC truth


  //----------------------------- V0 histos --------------------------------------//
  TH2F   *fHistV0RadiusZ[2];                        // V0 decay radius z filled for K0s and Lambda candidates
  TH2F   *fHistV0RadiusZVSPt[2];                    // V0 decay radius z vs pt filled for K0s and Lambda candidates
  TH2F   *fHistV0RadiusXY[2];                       // V0 decay radius x vs y filled for K0s and Lambda candidates
  TH2F   *fHistV0RadiusXYVSY[2];                    // V0 decay radius xy vs rapidity filled for K0s and Lambda candidates
   
  TH2F   *fHistArmenteros[2];                       // armenteros podolanski filled for K0s and Lambda candidates
     
  //-- K0 --//
  TH1F   *fHistPiPiMass[1];                         // pi+pi- InvMass spectrum
  TH2F   *fHistPiPiPtVSY[1];                        // pi+pi- InvMass spectrum vs rapidity
  TH2F   *fHistPiPiMassVSPt[1];                     // pi+pi- InvMass spectrum vs pt
  TH2F   *fHistPiPiMassVSPtMCTruth[1];              // pi+pi- InvMass spectrum vs pt MC truth
  // TH2F   *fHistPiPiMassVSAlpha[1];               // pi+pi- InvMass spectrum vs armenteros alpha
  TH2F   *fHistPiPiRadiusXY[1];                     // pi+pi- opening angle vs mass
  TH2F   *fHistPiPiCosPointAng[1];                  // pi+pi- cosine of pointing angle vs pt or dca to vertex
  TH2F   *fHistPiPiDCADaughterPosToPrimVtxVSMass[1];// dca of pos. K0s daughter to prim vtx vs mass
  TH2F   *fHistPiPiDecayLengthVsPt[1];              // pi+pi- decay lenght vs pt
  TH2F   *fHistPiPiDecayLengthVsMass[1];            // pi+pi- decay lenght vs pt
  TH2F   *fHistPiPiDecayLengthVsCtau[1];            // pi+pi- decay lenght vs pt
  TH2F   *fHistPiPiDCADaughters[1];                 // pi+pi- dca between daughters
  //  TH2F   *fHistPiPiPtDaughters[1];              // pi+pi- daughters pt pos vs pt neg 
  TH2F   *fHistPiPiDCAVSMass[1];                    // pi+pi- dca to prim vtx vs mass
  TH1F   *fHistPiPiMonitorCuts[1];                  // pi+pi- cut monitor
  TH1F   *fHistPiPiMonitorMCCuts[1];                // pi+pi- cut monitor mc
  TH2F   *fHistPiPiDecayLengthResolution[1];        // decay length mc reco vs mc truth K0s

  //-- Lambda --//
  TH1F   *fHistPiPMass[2];                          // p+pi- InvMass spectrum
  TH2F   *fHistPiPPtVSY[2];                         // p+pi- InvMass spectrum vs rapidity
  TH2F   *fHistPiPMassVSPt[2];                      // p+pi- InvMass spectrum vs pt
  TH2F   *fHistPiPMassVSPtMCTruth[2];               // p+pi- InvMass spectrum vs pt MC truth
  TH2F   *fHistPiPRadiusXY[2];                      // p+pi- opening angle vs mass
  TH2F   *fHistPiPCosPointAng[2];                   // p+pi- cosine of pointing angle vs pt  or dca to vertex
  TH2F   *fHistPiPDCADaughterPosToPrimVtxVSMass[2]; // dca of pos. Lambda daughter to prim vtx vs mass
  TH2F   *fHistPiPDecayLengthVsPt[2];               // p+pi- decay lenght vs pt
  TH2F   *fHistPiPDecayLengthVsMass[2];             // p+pi- decay lenght vs pt
  TH2F   *fHistPiPDecayLengthVsCtau[2];             // p+pi- decay lenght vs pt
  TH2F   *fHistPiPDCADaughters[2];                  // p+pi- dca between daughters
  //  TH2F   *fHistPiPPtDaughters[2];               // p+pi- daughters pt pos vs pt neg 
  TH2F   *fHistPiPDCAVSMass[2];                     // p+pi- dca to prim vtx vs mass
  TH1F   *fHistPiPMonitorCuts[2];                   // p+pi- cut monitor
  TH1F   *fHistPiPMonitorMCCuts[2];                 // p+pi- cut monitor mc
  TH2F   *fHistPiPMassVSPtSecSigma[2];              // InvMass distribution vs pt of secondary lambdas from sigma truth(0) reco(1)
  TH2F   *fHistPiPMassVSPtSecXi[2];                 // InvMass distribution vs pt of secondary lambdas from xi MC truth(0) reco(1)
  TH2F   *fHistPiPMassVSPtSecOmega[2];              // InvMass distribution vs pt of secondary lambdas from omega MC truth(0) reco(1)
  TH2F   *fHistPiPMassVSYSecXi[2];                  // InvMass distribution vs rapidity of secondary lambdas from xi MC truth(0) reco(1)
  TH2F   *fHistPiPXi0PtVSLambdaPt[2] ;              // pt of xi0 vs pt lambda truth(0) reco(1)
  TH2F   *fHistPiPXiMinusPtVSLambdaPt[2];           // pt of ximinus vs pt lambda truth(0) reco(1)
  TH2F   *fHistPiPOmegaPtVSLambdaPt[2];             // pt of omega plus vs pt alambda truth(0) reco(1)
  TH2F   *fHistPiPDecayLengthResolution[2];         // decay length mc reco vs mc truth lambda

  //-- AntiLambda --//
  TH1F   *fHistPiAPMass[2];                         // pi+p- InvMass spectrum
  TH2F   *fHistPiAPPtVSY[2];                        // pi+p- InvMass spectrum vs rapidity
  TH2F   *fHistPiAPMassVSPt[2];                     // pi+p- InvMass spectrum vs pt
  TH2F   *fHistPiAPMassVSPtMCTruth[2];              // pi+p- InvMass spectrum vs pt MC Truth
  TH2F   *fHistPiAPRadiusXY[2];                     // pi+p- opening angle vs mass
  TH2F   *fHistPiAPCosPointAng[2];                  // pi+p- cosine of pointing angle vs pt  or dca to vertex
  TH2F   *fHistPiAPDCADaughterPosToPrimVtxVSMass[2];// dca of pos. Lambda daughter to prim vtx vs mass
  TH2F   *fHistPiAPDecayLengthVsPt[2];              // pi+p- decay lenght vs pt
  TH2F   *fHistPiAPDecayLengthVsMass[2];            // pi+p- decay lenght vs pt
  TH2F   *fHistPiAPDecayLengthVsCtau[2];            // pi+p- decay lenght vs pt
  TH2F   *fHistPiAPDCADaughters[2];                 // pi+p- dca between daughters
  //  TH2F   *fHistPiAPPtDaughters[2];              // pi+p- daughters pt pos vs pt neg 
  TH2F   *fHistPiAPDCAVSMass[2];                    // pi+p- dca to prim vtx vs mass
  TH1F   *fHistPiAPMonitorCuts[2];                  // pi+p- cut monitor
  TH1F   *fHistPiAPMonitorMCCuts[2];                // pi+p- cut monitor mc
  TH2F   *fHistPiAPMassVSPtSecSigma[2];             // InvMass distribution vs pt of secondary alambdas from sigma truth(0) reco(1)
  TH2F   *fHistPiAPMassVSPtSecXi[2];                // InvMass distribution vs pt of secondary alambdas from xi MC truth(0) reco(1)
  TH2F   *fHistPiAPMassVSPtSecOmega[2];             // InvMass distribution vs pt of secondary alambdas from omega MC truth(0) reco(1)
  TH2F   *fHistPiAPMassVSYSecXi[2];                 // InvMass distribution vs rapidity of secondary alambdas from xi MC truth(0) reco(1)
  TH2F   *fHistPiAPXi0PtVSLambdaPt[2] ;             // pt of xi0 vs pt alambda truth(0) reco(1)
  TH2F   *fHistPiAPXiMinusPtVSLambdaPt[2];          // pt of ximinus vs pt alambda truth(0) reco(1)
  TH2F   *fHistPiAPOmegaPtVSLambdaPt[2];            // pt of omega plus vs pt alambda truth(0) reco(1)
  TH2F   *fHistPiAPDecayLengthResolution[2];        // decay length mc reco vs mc truth antilambda

  //-- others --//
  //dEdx
  TH2F   *fHistDedxSecProt[2];                      // dedx from proton cadidates vs pt
  TH2F   *fHistDedxSecAProt[2];                     // dedx from antiproton candidates vs pt
  TH2F   *fHistDedxSecPiMinus[2];                   // dedx from pi minus candidates vs pt
  TH2F   *fHistDedxSecPiPlus[2];                    // dedx from pi plus candidates vs pt

  TH2F   *fHistDedxProt[2];                         // dedx from proton cadidates vs pt before pidcut
  TH2F   *fHistDedxAProt[2];                        // dedx from antiproton candidates vs pt before pidcut
  TH2F   *fHistDedxPiMinus[2];                      // dedx from pi minus candidates vs pt before pidcut
  TH2F   *fHistDedxPiPlus[2];                       // dedx from pi plus candidates vs pt before pidcut
   
  //----------clusters and TPC var.------------//
  //K0s
  TH1F   *fHistNclsITSPosK0[1];                     // number of clusters from ITS of positive K0s daughters
  TH1F   *fHistNclsITSNegK0[1];                     // number of clusters from ITS of negative K0s daughters
  TH1F   *fHistNclsTPCPosK0[1];                     // number of clusters from TPC of positive K0s daughters
  TH1F   *fHistNclsTPCNegK0[1];                     // number of clusters from TPC of negative K0s daughters
  TH1F   *fHistChi2PerNclsITSPosK0[1];              // chi^2 per number of clusters ITS of positive K0s daughters
  TH1F   *fHistChi2PerNclsITSNegK0[1];              // chi^2 per number of clusters ITS of negative K0s daughters
  //Lambda
  TH1F   *fHistNclsITSPosL[2];                      // number of clusters from ITS of positive lambda daughters
  TH1F   *fHistNclsITSNegL[2];                      // number of clusters from ITS of negative lambda daughters
  TH1F   *fHistNclsTPCPosL[2];                      // number of clusters from TPC of positive lambda daughters
  TH1F   *fHistNclsTPCNegL[2];                      // number of clusters from TPC of negative lambda daughters
  TH1F   *fHistChi2PerNclsITSPosL[2];               // chi^2 per number of clusters ITS of positive lambda daughters
  TH1F   *fHistChi2PerNclsITSNegL[2];               // chi^2 per number of clusters ITS of negative lambda daughters
  //general
  TH2F   *fHistNclsITSPos[2];                       // number of clusters from ITS of positive daughters vs pt dautghter
  TH2F   *fHistNclsITSNeg[2];                       // number of clusters from ITS of negative daughters vs pt dautghter
  TH2F   *fHistNclsTPCPos[2];                       // number of clusters from TPC of positive daughters vs number of finabale clutsters
  TH2F   *fHistNclsTPCNeg[2];                       // number of clusters from TPC of negative daughters vs number of finabale clutsters
  TH2F   *fHistChi2PerNclsITSPos[2];                // chi^2 per number of clusters ITS of positive daughters vs pt of daughter
  TH2F   *fHistChi2PerNclsITSNeg[2];                // chi^2 per number of clusters ITS of negative daughters  vs pt of daughter
  TH2F   *fHistNclsITS[2];                          // number of clusters ITS pos vs neg daughters
  TH2F   *fHistNclsTPC[2];                          // number of clusters TPC  neg daughters vs number of crossed rows
  TH2F   *fHistNCRowsTPCPos[2];                     // number of crossed rows TPC pos. vs pt of daughter
  TH2F   *fHistNCRowsTPCNeg[2];                     // number of crossed rows TPC neg. vs pt of daughter
  TH2F   *fHistRatioFoundOverFinableTPCK0[2];       // ratio of ncls findable over found TPC K0s daughters
  TH2F   *fHistRatioFoundOverFinableTPCL[2];        // ratio of ncls findable over found TPC L daughters
  //eta all 
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
  Bool_t    fSelSDD;                   // select pp events with SDD (for pp 2.76TeV LHC11a)
  Bool_t    fSelNoSDD;                 // select pp events with no SDD (for pp 2.76TeV LHC11a)

  //onthefly
  Bool_t    fOntheFly;                 // true if onfly finder shall be used

  //vertex
  Double_t  fVertexZCut;               // z vertex cut value
  Bool_t    fVtxStatus;                // vertex cut on/off

  //PID
  Bool_t    fUsePID;                   // use pid yes/no
  Double_t  fNSigma;                   // set nsigma value
  Double_t  fPPIDcut;                  // set max momentum for pid cut usage 
  Double_t  fPtTPCCut;                 // low pt limit cut for TPC cluster cuts from AliESDtrackCuts
  Bool_t    fMoreNclsThanRows;         // cut on ncls>ncrossed rows yes/no
  Bool_t    fMoreNclsThanFindable;     // cut on ncls>nfindable cls yes/no
  Bool_t    fMoreNclsThanFindableMax;  // cut on ncls>nfindable max cls yes/no
  Double_t  fRatioFoundOverFindable;   // cut on found over findable clusters TPC
  Double_t  fRatioMaxCRowsOverFindable;// cut on crossed rows over finable max
  Double_t  fChi2PerClusterITS;        // cut on chi2 per ITS cluster
   
  //rapidity
  Bool_t    fRapCutV0;                 // use rapidity cut for V0 yes/no
  Double_t  fRap;                      // user defined value for rapidity cut

  //eta and pt
  Double_t  fEtaCutMCDaughters;        // eta cut for MC daughters on/off
  Double_t  fEtaCutMCDaughtersVal;     // eta cut value for MC daughters
  // Double_t  fEtaSignCut;            // eta cutsign daughters
  Double_t  fMinPt;                    // pt min cut value 

  //armenteros
  Double_t  fAlfaCut;                  // set alpha armenteros cut value
  Double_t  fQtCut;                    // set ptmax for qt armenteros cut 
  Double_t  fQtCutPtLow;               // set ptmin for  qt armenteros cut
  Bool_t    fArmCutK0;                 // set armenteros cut on/off for K0s
  Bool_t    fArmCutL;                  // set armenteros cut on/off for Lambda
  Double_t  fArmQtSlope;               // slope for armenteros K0s cut: qt = alpha*slope
  //others
  Double_t  fExcludeLambdaFromK0s;     // exlude Lambda mass from K0s throuh mass difference below this value
  Double_t  fExcludeK0sFromLambda;     // exlude K0s mass from Lambda throuh mass difference below this value
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
  
  Double_t  fCosPointAngL;             // cosine of pointing angle cut value for Lambda and ALambda
  Double_t  fCosPointAngK;             // cosine of pointing angle cut value for K0s
  Double_t  fCPAPtCutK0;               // pt max for cosine of pointing angle cut K0s
  Double_t  fCPAPtCutL;                // pt max for cosine of pointing angle cut Lambda
  Double_t  fOpengAngleDaughters;      // cut on opening angle between V0 daughters
  Double_t  fOpAngPtCut;               // max pt for using the  opening angle between V0 daughters cut
    
  Double_t  fDecayLengthMax;           // maximal decay length in x-y-z cut value
  Double_t  fDecayLengthMin;           // minimal decay length in x-y-z cut value

  //ctau
  Double_t  fCtauK0s;                  // multiple of ctau cut value for K0s
  Double_t  fCtauL;                    // multiple of ctau cut value for Lambda
  Double_t  fCtauPtCutK0;              // pt max for ctau cut usage for K0s
  Double_t  fCtauPtCutL;               // pt max for ctau cut usage for Lambda

  //KF particle chi cut
  //   Double_t  fChiCutKf;            // cut value of chi2 of AliKFParticle
  Bool_t  fChiCutKf;                   // cut value of chi2 of AliKFParticle

  
  /*
  // option for user defined charge/pt shift
  Bool_t     fShift;// shift yes/no
  Double_t   fDeltaInvP;//define shift value
  */
 

  AliAnalysisTaskV0ForRAA(const AliAnalysisTaskV0ForRAA&);
  AliAnalysisTaskV0ForRAA&operator=(const AliAnalysisTaskV0ForRAA&);
   
  ClassDef(AliAnalysisTaskV0ForRAA, 0); 
};
#endif
