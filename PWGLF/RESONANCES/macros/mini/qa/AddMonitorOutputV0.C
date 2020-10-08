/***************************************************************************
              anders.knospe@cern.ch - last modified on 15/05/2019

 Macro to add monitoring histograms for V0s to the rsn task

Options ("opt" argument):
- dim1 --> use TH1 only (no p or pt dependence)
***************************************************************************/

#ifndef ALIRSNADDMONITOROUTPUTV0_C
#define ALIRSNADDMONITOROUTPUTV0_C

#if !defined (__CINT__) || defined (__CLING__)
#include "AliRsnValueEvent.h"
#include "AliRsnValueDaughter.h"
#include "AliRsnListOutput.h"
#include "AliRsnLoopDaughter.h"
#endif

void AddMonitorOutputV0(Bool_t useMCMon, TObjArray *mon=0, TString v0name="K0S", TString opt="", AliRsnLoopDaughter *lm=0)
{
  //V0 type
  v0name.ToLower();
  Bool_t isK0 = v0name.Contains("k0");
  Bool_t isLambda = v0name.Contains("lambda");
  Bool_t isAntiLambda = (isLambda && (v0name.Contains("anti") || v0name.Contains("bar")));
  if(isAntiLambda) isLambda=false;
  if(!isK0 && !isLambda && !isAntiLambda){
    cerr<<"Error in AddMonitorOutputV0: unknown V0 type "<<v0name<<endl;
    isK0=true;
  }

  //Set options
  Bool_t useTH1 = opt.Contains("dim1");
  Bool_t doKine = !opt.Contains("nokine");//use "nokine" to skip kinematic histograms

  //Multiplicity/centrality
  AliRsnValueEvent* multi = new AliRsnValueEvent("multi",AliRsnValueEvent::kMult);
  multi->SetBins(0.0, 100.0, 1);
  //Momentum
  AliRsnValueDaughter* axisMomP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kP);
  axisMomP->SetBins(0.0, 10.0, 0.02);
  //pT
  AliRsnValueDaughter* axisMomPt = new AliRsnValueDaughter("pt", AliRsnValueDaughter::kPt);
  axisMomPt->SetBins(0.0,10.0,0.02);
  //Eta
  AliRsnValueDaughter* axisMomEta = new AliRsnValueDaughter("eta", AliRsnValueDaughter::kEta);
  axisMomEta->SetBins(-1.0, 1.0, 0.1);
  //Phi
  AliRsnValueDaughter* axisPhi = new AliRsnValueDaughter("phi", AliRsnValueDaughter::kPhi);
  axisPhi->SetBins(0.0, 360.0, 1.0);

  //V0 Momentum
  AliRsnValueDaughter* axisV0MomP = new AliRsnValueDaughter("v0p", AliRsnValueDaughter::kV0P);
  axisV0MomP->SetBins(0.0, 10.0, 0.02);
  //V0 pT
  AliRsnValueDaughter* axisV0MomPt = new AliRsnValueDaughter("v0pt", AliRsnValueDaughter::kV0Pt);
  axisV0MomPt->SetBins(0.0,10.0,0.02);

  //Mass
  AliRsnValueDaughter* axisMass = new AliRsnValueDaughter("Mass",AliRsnValueDaughter::kV0Mass);
  if(isK0) axisMass->SetBins(0.4,0.6,0.001);
  else axisMass->SetBins(1.08,1.16,0.001);//Lambda
  //Mass (K0S hypothesis)
  AliRsnValueDaughter* axisK0SMass = new AliRsnValueDaughter("K0SMass",AliRsnValueDaughter::kK0SMass);
  if(isK0) axisK0SMass->SetBins(0.4,0.6,0.001);
  else axisK0SMass->SetBins(0.25,1.5,0.005);//misidentified
  //Mass (Lambda hypothesis)
  AliRsnValueDaughter* axisLambdaMass = new AliRsnValueDaughter("LambdaMass",AliRsnValueDaughter::kLambdaMass);
  if(isLambda) axisLambdaMass->SetBins(1.08,1.16,0.001);
  else axisLambdaMass->SetBins(1.,2.,0.005);//misidentified
  //Mass (anti-Lambda hypothesis)
  AliRsnValueDaughter* axisAntiLambdaMass = new AliRsnValueDaughter("AntiLambdaMass",AliRsnValueDaughter::kAntiLambdaMass);
  if(isAntiLambda) axisAntiLambdaMass->SetBins(1.08,1.16,0.001);
  else axisAntiLambdaMass->SetBins(1.,2.,0.005);//misidentified
  //DCA
  AliRsnValueDaughter* axisDCA = new AliRsnValueDaughter("DCAtoVertex",AliRsnValueDaughter::kV0DCA);
  axisDCA->SetBins(0.0,4.0,0.001);
  //Radius
  AliRsnValueDaughter* axisRadius = new AliRsnValueDaughter("Radius",AliRsnValueDaughter::kV0Radius);
  axisRadius->SetBins(0.0,200.0,0.1);
  //Lifetime
  AliRsnValueDaughter* axisLifetime = new AliRsnValueDaughter("Lifetime",AliRsnValueDaughter::kV0Lifetime);
  axisLifetime->SetBins(0.0,200.0,0.1);
  //Daughter DCA
  AliRsnValueDaughter* axisDDCA = new AliRsnValueDaughter("DaughtersDCA",AliRsnValueDaughter::kDaughterDCA);
  axisDDCA->SetBins(0.0,2.0,0.001);
  //Cosine of Pointing Angle
  AliRsnValueDaughter* axisCPA = new AliRsnValueDaughter("CosPointAng",AliRsnValueDaughter::kCosPointAng);
  axisCPA->SetBins(0.96,1.0,0.0002);

  //Negative Daughter pT
  AliRsnValueDaughter* axisNpt = new AliRsnValueDaughter("NegDaughterPt",AliRsnValueDaughter::kV0NPt);
  axisNpt->SetBins(0.0,10.0,0.05);
  //Positive Daughter pT
  AliRsnValueDaughter* axisPpt = new AliRsnValueDaughter("PosDaughterPt",AliRsnValueDaughter::kV0PPt);
  axisPpt->SetBins(0.0,10.0,0.05);
  //DCA of Daughters to Primary Vertex
  AliRsnValueDaughter* axisDCA2T = new AliRsnValueDaughter("DaughtersDCAPrimVtx",AliRsnValueDaughter::kV0DCAXY);
  axisDCA2T->SetBins(-10.0,10.0,0.01);
  //pi- TPC PID
  AliRsnValueDaughter* axisPimPID = new AliRsnValueDaughter("V0PimTPC",AliRsnValueDaughter::kLambdaPionPIDCut);
  axisPimPID->SetBins(0.0,10.0,0.05);
  //pi+ TPC PID
  AliRsnValueDaughter* axisPipPID = new AliRsnValueDaughter("V0PipTPC",AliRsnValueDaughter::kAntiLambdaAntiPionPIDCut);
  axisPipPID->SetBins(0.0,10.0,0.05);
  //Proton TPC PID
  AliRsnValueDaughter* axisProtonPID = new AliRsnValueDaughter("V0ProtonPID",AliRsnValueDaughter::kLambdaProtonPIDCut);
  axisProtonPID->SetBins(0.0,10.0,0.05);
  //Antiproton TPC PID
  AliRsnValueDaughter* axisAntiprotonPID = new AliRsnValueDaughter("V0AntiprotonPID",AliRsnValueDaughter::kAntiLambdaAntiProtonPIDCut);
  axisAntiprotonPID->SetBins(0.0,10.0,0.05);

  /****************************************************************/
  /***************         KINEMATICS          ********************/
  /****************************************************************/
  if(doKine) {
    // output: 1D histogram of momentum
    AliRsnListOutput* outMonitorP = new AliRsnListOutput("P", AliRsnListOutput::kHistoDefault);
    outMonitorP->AddValue(axisMomP);
    if (mon) mon->Add(outMonitorP);
    if (lm) lm->AddOutput(outMonitorP);

    // output: 1D histogram of pt
    AliRsnListOutput* outMonitorPt = new AliRsnListOutput("Pt", AliRsnListOutput::kHistoDefault);
    outMonitorPt->AddValue(axisMomPt);
    if (mon) mon->Add(outMonitorPt);
    if (lm) lm->AddOutput(outMonitorPt);

    // output: 1D histogram of pseudorapidity
    AliRsnListOutput* outMonitorEta = new AliRsnListOutput("Eta", AliRsnListOutput::kHistoDefault);
    outMonitorEta->AddValue(axisMomEta);
    if (mon) mon->Add(outMonitorEta);
    if (lm) lm->AddOutput(outMonitorEta);

    // output: 1D histogram of phi at vertex
    AliRsnListOutput* outMonitorPhi = new AliRsnListOutput("Phi", AliRsnListOutput::kHistoDefault);
    outMonitorPhi->AddValue(axisPhi);
    if (mon) mon->Add(outMonitorPhi);
    if (lm) lm->AddOutput(outMonitorPhi);

    // output: 2D histogram of phi vs pt
    AliRsnListOutput* outMonitorPhiVsPt = new AliRsnListOutput("PhiVsPt", AliRsnListOutput::kHistoDefault);
    outMonitorPhiVsPt->AddValue(axisMomPt);
    outMonitorPhiVsPt->AddValue(axisPhi);
    if (mon) mon->Add(outMonitorPhiVsPt);
    if (lm) lm->AddOutput(outMonitorPhiVsPt);

    // output: 1D histogram of V0 momentum
    AliRsnListOutput* outMonitorV0P = new AliRsnListOutput("V0P", AliRsnListOutput::kHistoDefault);
    outMonitorV0P->AddValue(axisV0MomP);
    if (mon) mon->Add(outMonitorV0P);
    if (lm) lm->AddOutput(outMonitorV0P);

    // output: 1D histogram of V0 pt
    AliRsnListOutput* outMonitorV0Pt = new AliRsnListOutput("V0Pt", AliRsnListOutput::kHistoDefault);
    outMonitorV0Pt->AddValue(axisV0MomPt);
    if (mon) mon->Add(outMonitorV0Pt);
    if (lm) lm->AddOutput(outMonitorV0Pt);
  }

  /****************************************************************/
  /***************    TOPOLOGICAL VARIABLES    ********************/
  /****************************************************************/
  // output: 2D histogram of mass vs pt
  AliRsnListOutput* outMonitorMass = new AliRsnListOutput("V0MassVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorMass->AddValue(axisMomPt);
  outMonitorMass->AddValue(axisMass);
  if (mon) mon->Add(outMonitorMass);
  if (lm) lm->AddOutput(outMonitorMass);

  // output: 2D histogram of mass (K0S hypothesis) vs pt
  AliRsnListOutput* outMonitorK0SMass = new AliRsnListOutput("K0SMassVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorK0SMass->AddValue(axisMomPt);
  outMonitorK0SMass->AddValue(axisK0SMass);
  if (mon) mon->Add(outMonitorK0SMass);
  if (lm) lm->AddOutput(outMonitorK0SMass);

  // output: 2D histogram of mass (Lambda hypothesis) vs pt
  AliRsnListOutput* outMonitorLambdaMass = new AliRsnListOutput("LambdaMassVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorLambdaMass->AddValue(axisMomPt);
  outMonitorLambdaMass->AddValue(axisLambdaMass);
  if (mon) mon->Add(outMonitorLambdaMass);
  if (lm) lm->AddOutput(outMonitorLambdaMass);

  // output: 2D histogram of mass (anti-Lambda hypothesis) vs pt
  AliRsnListOutput* outMonitorAntiLambdaMass = new AliRsnListOutput("AntiLambdaMassVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorAntiLambdaMass->AddValue(axisMomPt);
  outMonitorAntiLambdaMass->AddValue(axisAntiLambdaMass);
  if (mon) mon->Add(outMonitorAntiLambdaMass);
  if (lm) lm->AddOutput(outMonitorAntiLambdaMass);

  // output: 2D histogram of V0 DCA vs pt
  AliRsnListOutput* outMonitorDCA = new AliRsnListOutput("V0DCAtoVertexVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorDCA->AddValue(axisMomPt);
  outMonitorDCA->AddValue(axisDCA);
  if (mon) mon->Add(outMonitorDCA);
  if (lm) lm->AddOutput(outMonitorDCA);

  // output: 2D histogram of radius vs pt
  AliRsnListOutput* outMonitorRadius = new AliRsnListOutput("V0RadiusVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorRadius->AddValue(axisMomPt);
  outMonitorRadius->AddValue(axisRadius);
  if (mon) mon->Add(outMonitorRadius);
  if (lm) lm->AddOutput(outMonitorRadius);

  // output: 2D histogram of lifetime vs p
  AliRsnListOutput* outMonitorLifetime = new AliRsnListOutput("V0LifetimeVsP", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorLifetime->AddValue(axisMomP);
  outMonitorLifetime->AddValue(axisLifetime);
  if (mon) mon->Add(outMonitorLifetime);
  if (lm) lm->AddOutput(outMonitorLifetime);

  // output: 2D histogram of Daughter DCA (to each other) vs pt
  AliRsnListOutput* outMonitorDDCA = new AliRsnListOutput("V0DaughtersDCAVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorDDCA->AddValue(axisMomPt);
  outMonitorDDCA->AddValue(axisDDCA);
  if (mon) mon->Add(outMonitorDDCA);
  if (lm) lm->AddOutput(outMonitorDDCA);

  // output: 2D histogram of Cosine of Pointing Angle vs p
  AliRsnListOutput* outMonitorCPA = new AliRsnListOutput("V0CosPointAngVsP", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorCPA->AddValue(axisMomP);
  outMonitorCPA->AddValue(axisCPA);
  if (mon) mon->Add(outMonitorCPA);
  if (lm) lm->AddOutput(outMonitorCPA);

  /****************************************************************/
  /***************          DAUGHTERS          ********************/
  /****************************************************************/
  // output: 1D histogram of negative daughter pt
  AliRsnListOutput* outMonitorNPt = new AliRsnListOutput("NegDaughterPt", AliRsnListOutput::kHistoDefault);
  outMonitorNPt->AddValue(axisNpt);
  if (mon) mon->Add(outMonitorNPt);
  if (lm) lm->AddOutput(outMonitorNPt);

  // output: 1D histogram of negative daughter pt
  AliRsnListOutput* outMonitorPPt = new AliRsnListOutput("PosDaughterPt", AliRsnListOutput::kHistoDefault);
  outMonitorPPt->AddValue(axisPpt);
  if (mon) mon->Add(outMonitorPPt);
  if (lm) lm->AddOutput(outMonitorPPt);

  // output: 1D histogram of DCA of daughters to primary vertex
  AliRsnListOutput* outMonitorDCA2T = new AliRsnListOutput("V0DaughtersDCAPrimVtx", AliRsnListOutput::kHistoDefault);
  outMonitorDCA2T->AddValue(axisDCA2T);
  if (mon) mon->Add(outMonitorDCA2T);
  if (lm) lm->AddOutput(outMonitorDCA2T);

  if(isK0 || isLambda){
    // output: 2D histogram of pi- TPC PID vs pt
    AliRsnListOutput* outMonitorPimPID = new AliRsnListOutput("V0PimPIDVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorPimPID->AddValue(axisNpt);
    outMonitorPimPID->AddValue(axisPimPID);
    if (mon) mon->Add(outMonitorPimPID);
    if (lm) lm->AddOutput(outMonitorPimPID);
  }

  if(isK0 || isAntiLambda){
    // output: 2D histogram of pi+ TPC PID vs pt
    AliRsnListOutput* outMonitorPipPID = new AliRsnListOutput("V0PipPIDVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorPipPID->AddValue(axisPpt);
    outMonitorPipPID->AddValue(axisPipPID);
    if (mon) mon->Add(outMonitorPipPID);
    if (lm) lm->AddOutput(outMonitorPipPID);
  }

  if(isLambda){
    // output: 2D histogram of proton TPC PID vs pt
    AliRsnListOutput* outMonitorProtonPID = new AliRsnListOutput("V0ProtonPIDVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorProtonPID->AddValue(axisPpt);
    outMonitorProtonPID->AddValue(axisProtonPID);
    if (mon) mon->Add(outMonitorProtonPID);
    if (lm) lm->AddOutput(outMonitorProtonPID);
  }

  if(isAntiLambda){
    // output: 2D histogram of antiproton TPC PID vs pt
    AliRsnListOutput* outMonitorAntiprotonPID = new AliRsnListOutput("V0AntiprotonPIDVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorAntiprotonPID->AddValue(axisNpt);
    outMonitorAntiprotonPID->AddValue(axisAntiprotonPID);
    if (mon) mon->Add(outMonitorAntiprotonPID);
    if (lm) lm->AddOutput(outMonitorAntiprotonPID);
  }

  /****************************************************************/
  /***************       MONITOR MC            ********************/
  /****************************************************************/
  if (useMCMon) {
    //Momentum
    AliRsnValueDaughter* axisMomPMC = new AliRsnValueDaughter("pMC", AliRsnValueDaughter::kP);
    axisMomPMC->SetUseMCInfo(kTRUE);
    axisMomPMC->SetBins(0.0,10.0,0.01);
    //Pt
    AliRsnValueDaughter* axisMomPtMC = new AliRsnValueDaughter("ptMC", AliRsnValueDaughter::kPt);
    axisMomPtMC->SetUseMCInfo(kTRUE);
    axisMomPtMC->SetBins(0.0,10.0,0.01);
    //Eta
    AliRsnValueDaughter* axisEtaMC = new AliRsnValueDaughter("etaMC", AliRsnValueDaughter::kEta);
    axisEtaMC->SetUseMCInfo(kTRUE);
    axisEtaMC->SetBins(-1.0,1.0,0.1);

    // output: 2D histo for kine acceptance
    AliRsnListOutput* outMonitorPtVsEtaMC = new AliRsnListOutput("V0Pt_VsEtaMC", AliRsnListOutput::kHistoDefault);
    outMonitorPtVsEtaMC->AddValue(axisMomPtMC);
    outMonitorPtVsEtaMC->AddValue(axisEtaMC);
    if (mon) mon->Add(outMonitorPtVsEtaMC);
    if (lm) lm->AddOutput(outMonitorPtVsEtaMC);

    // output: 1D histo pt from MC
    AliRsnListOutput* outMonitorPtMC = new AliRsnListOutput("V0PtMC", AliRsnListOutput::kHistoDefault);
    outMonitorPtMC->AddValue(axisMomPtMC);
    if (mon) mon->Add(outMonitorPtMC);
    if (lm) lm->AddOutput(outMonitorPtMC);
 
    // output: 1D histo eta from MC
    AliRsnListOutput* outMonitorEtaMC = new AliRsnListOutput("V0EtaMC", AliRsnListOutput::kHistoDefault);
    outMonitorEtaMC->AddValue(axisEtaMC);
    if (mon) mon->Add(outMonitorEtaMC);
    if (lm) lm->AddOutput(outMonitorEtaMC);
  }

  return;
}

#endif
