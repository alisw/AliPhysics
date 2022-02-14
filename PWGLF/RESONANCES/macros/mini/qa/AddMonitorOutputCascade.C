/***************************************************************************
              anders.knospe@cern.ch - last modified on 16/05/2019

 Macro to add monitoring histograms for V0s to the rsn task

Options ("opt" argument):
- dim1 --> use TH1 only (no p or pt dependence)
***************************************************************************/

#ifndef ALIRSNADDMONITOROUTPUTCASCADE_C
#define ALIRSNADDMONITOROUTPUTCASCADE_C

#if !defined (__CINT__) || defined (__CLING__)
#include "AliRsnValueEvent.h"
#include "AliRsnValueDaughter.h"
#include "AliRsnListOutput.h"
#include "AliRsnLoopDaughter.h"
#endif

void AddMonitorOutputCascade(Bool_t useMCMon, TObjArray *mon=0, TString caname="Xi", TString opt="", AliRsnLoopDaughter *lm=0)
{
  //Cascade type
  caname.ToLower();
  Bool_t isXi = caname.Contains("xi");
  Bool_t isOmega = caname.Contains("omega");
  //Short_t charge = (caname.Contains("plus") || caname.Contains("anti") || caname.Contains("bar")) ? 1 : -1;
  if(!isXi && !isOmega){
    cerr<<"Error in AddMonitorOutputCascade: unknown cascade type "<<caname<<endl;
    isXi=true;
  }

  //Set options
  Bool_t useTH1 = opt.Contains("dim1");

  //Multiplicity/centrality
  AliRsnValueEvent* multi = new AliRsnValueEvent("multi",AliRsnValueEvent::kMult);
  multi->SetBins(0.0, 100.0, 1);
  //Momentum
  AliRsnValueDaughter* axisMomP = new AliRsnValueDaughter("p", AliRsnValueDaughter::kCascadeP);
  axisMomP->SetBins(0.0, 10.0, 0.02);
  //pT
  AliRsnValueDaughter* axisMomPt = new AliRsnValueDaughter("pt", AliRsnValueDaughter::kCascadePt);
  axisMomPt->SetBins(0.0,10.0,0.02);
  //Eta (for AODs: this gives the V0 eta)
  AliRsnValueDaughter* axisMomEta = new AliRsnValueDaughter("eta", AliRsnValueDaughter::kEta);
  axisMomEta->SetBins(-2.0, 2.0, 0.1);
  //Phi (for AODs: this gives the V0 phi)
  AliRsnValueDaughter* axisPhi = new AliRsnValueDaughter("phi", AliRsnValueDaughter::kPhi);
  axisPhi->SetBins(0.0, 360.0, 1.0);

  //Mass (Xi hypothesis)
  AliRsnValueDaughter* axisXiMass = new AliRsnValueDaughter("XiMass",AliRsnValueDaughter::kXiMass);
  if(isXi) axisXiMass->SetBins(1.25,1.4,0.001);
  else axisXiMass->SetBins(1.25,2.,0.005);//misidentified
  //Mass (Omega hypothesis)
  AliRsnValueDaughter* axisOmegaMass = new AliRsnValueDaughter("OmegaMass",AliRsnValueDaughter::kOmegaMass);
  if(isOmega) axisOmegaMass->SetBins(1.6,1.9,0.001);
  else axisOmegaMass->SetBins(1.6,2.5,0.005);//misidentified
  //DCA
  AliRsnValueDaughter* axisDCA = new AliRsnValueDaughter("DCAtoVertex",AliRsnValueDaughter::kCascadeDCA);
  axisDCA->SetBins(-0.4,0.4,0.001);
  //Radius
  AliRsnValueDaughter* axisRadius = new AliRsnValueDaughter("Radius",AliRsnValueDaughter::kCascadeRadius);
  axisRadius->SetBins(0.0,200,0.2);
  //Daughter DCA
  AliRsnValueDaughter* axisDDCA = new AliRsnValueDaughter("CascadeDaughtersDCA",AliRsnValueDaughter::kCascadeDaughterDCA);
  axisDDCA->SetBins(-2.0,2.0,0.001);
  //Cosine of Pointing Angle
  AliRsnValueDaughter* axisCPA = new AliRsnValueDaughter("CascadeCosPointAng",AliRsnValueDaughter::kCascadeCosPointAng);
  axisCPA->SetBins(0.96,1.,0.001);
  //Cosine of V0 Pointing Angle to Cascade Vertex
  AliRsnValueDaughter* axisV0CPA = new AliRsnValueDaughter("CascadeV0CosPointAng",AliRsnValueDaughter::kCascadeV0CosPointAng);
  axisV0CPA->SetBins(0.96,1.,0.001);
  //V0 Lifetime
  AliRsnValueDaughter* axisV0Lifetime = new AliRsnValueDaughter("CascadeV0Lifetime",AliRsnValueDaughter::kCascadeV0Lifetime);
  axisV0Lifetime->SetBins(0.0,200,0.1);

  //V0 pT
  AliRsnValueDaughter* axisV0pt = new AliRsnValueDaughter("CascadeV0Pt",AliRsnValueDaughter::kCascadeV0Pt);
  axisV0pt->SetBins(0.,10.0,0.05);
  //Bachelor Track pT
  AliRsnValueDaughter* axisBpt = new AliRsnValueDaughter("BachelorPt",AliRsnValueDaughter::kBachelorPt);
  axisBpt->SetBins(0.,10.0,0.05);
  //pion TPC PID
  AliRsnValueDaughter* axisPiPID = new AliRsnValueDaughter("CascadeBachelorPionTPC",AliRsnValueDaughter::kBachelorPionTPCnsigma);
  axisPiPID->SetBins(-10.,10.,0.05);
  //kaon TPC PID
  AliRsnValueDaughter* axisKPID = new AliRsnValueDaughter("CascadeBachelorKaonTPC",AliRsnValueDaughter::kBachelorKaonTPCnsigma);
  axisKPID->SetBins(-10.,10.,0.05);

  /****************************************************************/
  /***************         KINEMATICS          ********************/
  /****************************************************************/
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

  /****************************************************************/
  /***************    TOPOLOGICAL VARIABLES    ********************/
  /****************************************************************/
  // output: 2D histogram of mass (Xi hypothesis) vs pt
  AliRsnListOutput* outMonitorXiMass = new AliRsnListOutput("XiMassVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorXiMass->AddValue(axisMomPt);
  outMonitorXiMass->AddValue(axisXiMass);
  if (mon) mon->Add(outMonitorXiMass);
  if (lm) lm->AddOutput(outMonitorXiMass);

  // output: 2D histogram of mass (Omega hypothesis) vs pt
  AliRsnListOutput* outMonitorOmegaMass = new AliRsnListOutput("OmegaMassVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorOmegaMass->AddValue(axisMomPt);
  outMonitorOmegaMass->AddValue(axisOmegaMass);
  if (mon) mon->Add(outMonitorOmegaMass);
  if (lm) lm->AddOutput(outMonitorOmegaMass);

  // output: 2D histogram of Cascade DCA vs pt
  AliRsnListOutput* outMonitorDCA = new AliRsnListOutput("CascadeDCAtoVertexVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorDCA->AddValue(axisMomPt);
  outMonitorDCA->AddValue(axisDCA);
  if (mon) mon->Add(outMonitorDCA);
  if (lm) lm->AddOutput(outMonitorDCA);

  // output: 2D histogram of radius vs pt
  AliRsnListOutput* outMonitorRadius = new AliRsnListOutput("CascadeRadiusVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorRadius->AddValue(axisMomPt);
  outMonitorRadius->AddValue(axisRadius);
  if (mon) mon->Add(outMonitorRadius);
  if (lm) lm->AddOutput(outMonitorRadius);

  // output: 2D histogram of Daughter DCA (to each other) vs pt
  AliRsnListOutput* outMonitorDDCA = new AliRsnListOutput("CascadeDaughtersDCAVsPt", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorDDCA->AddValue(axisMomPt);
  outMonitorDDCA->AddValue(axisDDCA);
  if (mon) mon->Add(outMonitorDDCA);
  if (lm) lm->AddOutput(outMonitorDDCA);

  // output: 2D histogram of Cosine of Pointing Angle vs p
  AliRsnListOutput* outMonitorCPA = new AliRsnListOutput("CascadeCosPointAngVsP", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorCPA->AddValue(axisMomP);
  outMonitorCPA->AddValue(axisCPA);
  if (mon) mon->Add(outMonitorCPA);
  if (lm) lm->AddOutput(outMonitorCPA);

  // output: 2D histogram of Cosine of V0 Pointing Angle to Cascade Vertex vs p
  AliRsnListOutput* outMonitorV0CPA = new AliRsnListOutput("CascadeV0CosPointAngVsP", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorV0CPA->AddValue(axisMomP);
  outMonitorV0CPA->AddValue(axisV0CPA);
  if (mon) mon->Add(outMonitorV0CPA);
  if (lm) lm->AddOutput(outMonitorV0CPA);

  // output: 2D histogram of V0 lifetime vs p
  AliRsnListOutput* outMonitorV0Lifetime = new AliRsnListOutput("CascadeV0LifetimeVsP", AliRsnListOutput::kHistoDefault);
  if (!useTH1) outMonitorV0Lifetime->AddValue(axisMomP);
  outMonitorV0Lifetime->AddValue(axisV0Lifetime);
  if (mon) mon->Add(outMonitorV0Lifetime);
  if (lm) lm->AddOutput(outMonitorV0Lifetime);

  /****************************************************************/
  /***************          DAUGHTERS          ********************/
  /****************************************************************/
  // output: 1D histogram of pt of V0 from cascade
  AliRsnListOutput* outMonitorV0Pt = new AliRsnListOutput("CascadeV0Pt", AliRsnListOutput::kHistoDefault);
  outMonitorV0Pt->AddValue(axisV0pt);
  if (mon) mon->Add(outMonitorV0Pt);
  if (lm) lm->AddOutput(outMonitorV0Pt);

  // output: 1D histogram of bachelor track pt
  AliRsnListOutput* outMonitorBPt = new AliRsnListOutput("BachelorPt", AliRsnListOutput::kHistoDefault);
  outMonitorBPt->AddValue(axisBpt);
  if (mon) mon->Add(outMonitorBPt);
  if (lm) lm->AddOutput(outMonitorBPt);

  if(isXi){
    // output: 2D histogram of pion TPC PID vs pt
    AliRsnListOutput* outMonitorPionPID = new AliRsnListOutput("CascadeBachelorPionTPCVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorPionPID->AddValue(axisBpt);
    outMonitorPionPID->AddValue(axisPiPID);
    if (mon) mon->Add(outMonitorPionPID);
    if (lm) lm->AddOutput(outMonitorPionPID);
  }

  if(isOmega){
    // output: 2D histogram of kaon TPC PID vs pt
    AliRsnListOutput* outMonitorKaonPID = new AliRsnListOutput("CascadeBachelorPionTPCVsPt", AliRsnListOutput::kHistoDefault);
    if (!useTH1) outMonitorKaonPID->AddValue(axisBpt);
    outMonitorKaonPID->AddValue(axisKPID);
    if (mon) mon->Add(outMonitorKaonPID);
    if (lm) lm->AddOutput(outMonitorKaonPID);
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
    AliRsnListOutput* outMonitorPtVsEtaMC = new AliRsnListOutput("CascadePt_VsEtaMC", AliRsnListOutput::kHistoDefault);
    outMonitorPtVsEtaMC->AddValue(axisMomPtMC);
    outMonitorPtVsEtaMC->AddValue(axisEtaMC);
    if (mon) mon->Add(outMonitorPtVsEtaMC);
    if (lm) lm->AddOutput(outMonitorPtVsEtaMC);

    // output: 1D histo pt from MC
    AliRsnListOutput* outMonitorPtMC = new AliRsnListOutput("CascadePtMC", AliRsnListOutput::kHistoDefault);
    outMonitorPtMC->AddValue(axisMomPtMC);
    if (mon) mon->Add(outMonitorPtMC);
    if (lm) lm->AddOutput(outMonitorPtMC);
 
    // output: 1D histo eta from MC
    AliRsnListOutput* outMonitorEtaMC = new AliRsnListOutput("CascadeEtaMC", AliRsnListOutput::kHistoDefault);
    outMonitorEtaMC->AddValue(axisEtaMC);
    if (mon) mon->Add(outMonitorEtaMC);
    if (lm) lm->AddOutput(outMonitorEtaMC);
  }
    
  return;
}

#endif
