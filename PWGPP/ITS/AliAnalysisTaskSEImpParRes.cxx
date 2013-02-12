/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the study of the impact parameter resolution
//
// Authors:A.Dainese,    andrea.dainese@pd.infn.it
//     and Xianbao Yuan, yuanxb@iopp.ccnu.edu.cn; xianbao.yuan@pd.infn.it
/////////////////////////////////////////////////////////

#include <TList.h>
#include <TH1F.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "AliGeomManager.h"
#include "AliMultiplicity.h"
#include "AliTriggerClass.h"
#include "AliTriggerCluster.h"
#include "AliTriggerConfiguration.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTrackPointArray.h"
#include "AliMCEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"   
#include "AliVertexerTracks.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliPID.h"
#include "AliAnalysisTaskSEImpParRes.h"

ClassImp(AliAnalysisTaskSEImpParRes)

//________________________________________________________________________
AliAnalysisTaskSEImpParRes::AliAnalysisTaskSEImpParRes():
  AliAnalysisTaskSE(),
  fIsAOD(kFALSE),
  fReadMC(kFALSE),
  fSelectedPdg(-1),
  fUseDiamond(kFALSE),
  fSkipTrack(kTRUE),
  fMinMult(0),
  fMaxMult(1000000),
  fCheckSDDIsIn(0),
  fTriggerClass(""),
  fTrigConfig(0),
  fOCDBPath(""),
  fESDtrackCuts(0),
  fOutputitspureSARec(0),
  fOutputitspureSASkip(0), 
  fOutputallPointRec(0),
  fOutputallPointSkip(0),
  fOutputpartPointRec(0),
  fOutputpartPointSkip(0),
  fOutputonepointSPDRec(0),
  fOutputonepointSPDSkip(0),
  fOutputpostvTracRec(0),
  fOutputpostvTracSkip(0),
  fOutputnegtvTracRec(0),
  fOutputnegtvTracSkip(0),
  fOutputpullAllpointRec(0),
  fOutputpullAllpointSkip(0),
  fOutputOnlyRefitRec(0),
  fOutputOnlyRefitSkip(0),
  fOutputSinThetaRec(0),
  fOutputSinThetaSkip(0),
  fOutputallPointTrue(0),
  fOutputpostvTracTrue(0),
  fOutputnegtvTracTrue(0),
  fOutputpullAllpointTrue(0),
  fOutputphiAllpointSkip(0),
  fOutputphiPostvtracSkip(0),
  fOutputphiNegtvtracSkip(0),
  fOutputparticlePID(0),
  fOutputWithTrackCuts(0),
  fOutputPt(0),
  fNentries(0),
  fEstimVtx(0)
{
  //
  // Default constructor
  //
}

//________________________________________________________________________
AliAnalysisTaskSEImpParRes::AliAnalysisTaskSEImpParRes(const char *name):
  AliAnalysisTaskSE(name),
  fIsAOD(kFALSE),
  fReadMC(kFALSE),
  fSelectedPdg(-1),
  fUseDiamond(kFALSE),
  fSkipTrack(kTRUE),
  fMinMult(0),
  fMaxMult(1000000),
  fCheckSDDIsIn(0),
  fTriggerClass(""),
  fTrigConfig(0),
  fOCDBPath(""),
  fESDtrackCuts(0),
  fOutputitspureSARec(0),
  fOutputitspureSASkip(0), 
  fOutputallPointRec(0),
  fOutputallPointSkip(0),
  fOutputpartPointRec(0),
  fOutputpartPointSkip(0),
  fOutputonepointSPDRec(0),
  fOutputonepointSPDSkip(0),
  fOutputpostvTracRec(0),
  fOutputpostvTracSkip(0),
  fOutputnegtvTracRec(0),
  fOutputnegtvTracSkip(0),
  fOutputpullAllpointRec(0),
  fOutputpullAllpointSkip(0),
  fOutputOnlyRefitRec(0),
  fOutputOnlyRefitSkip(0),
  fOutputSinThetaRec(0),
  fOutputSinThetaSkip(0),
  fOutputallPointTrue(0),
  fOutputpostvTracTrue(0),
  fOutputnegtvTracTrue(0),
  fOutputpullAllpointTrue(0),
  fOutputphiAllpointSkip(0),
  fOutputphiPostvtracSkip(0),
  fOutputphiNegtvtracSkip(0),
  fOutputparticlePID(0),
  fOutputWithTrackCuts(0),
  fOutputPt(0),
  fNentries(0),
  fEstimVtx(0)
{
  //
  // Default constructor
  //

  DefineOutput(1, TList::Class());  //My private output
  DefineOutput(2, TList::Class());  //My private output
  DefineOutput(3, TList::Class());  //My private output
  DefineOutput(4, TList::Class());  //My private output
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());  //My private output
  DefineOutput(7, TList::Class());
  DefineOutput(8, TList::Class());  //My private output
  DefineOutput(9, TList::Class());  //My private output
  DefineOutput(10, TList::Class());  //My private output
  DefineOutput(11, TList::Class());  //My private output
  DefineOutput(12, TList::Class());
  DefineOutput(13, TList::Class());  //My private output
  DefineOutput(14, TList::Class());
  DefineOutput(15, TList::Class());  //My private output
  DefineOutput(16, TList::Class());
  DefineOutput(17, TList::Class());  //My private output
  DefineOutput(18, TList::Class());
  DefineOutput(19, TList::Class());  //My private output
  DefineOutput(20, TList::Class());  //My private output
  DefineOutput(21, TList::Class());
  DefineOutput(22, TList::Class());  //My private output
  DefineOutput(23, TList::Class());
  DefineOutput(24, TList::Class());  //My private output
  DefineOutput(25, TList::Class());
  DefineOutput(26, TList::Class());  //My private output
  DefineOutput(27, TList::Class());
  DefineOutput(28, TH1F::Class());  //My private output
  DefineOutput(29, TH1F::Class());
  DefineOutput(30, TList::Class());
  
}

//________________________________________________________________________
AliAnalysisTaskSEImpParRes::~AliAnalysisTaskSEImpParRes()
{
  //
  // default distructor  
  // 
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return; // RS
  //
  if (fESDtrackCuts) {    delete fESDtrackCuts;  fESDtrackCuts = 0;  }
  if (fOutputitspureSARec)                      { delete fOutputitspureSARec; fOutputitspureSARec=0x0;}
  if (fOutputitspureSASkip)                   { delete fOutputitspureSASkip; fOutputitspureSASkip=0x0;}
  if (fOutputallPointRec)                        { delete fOutputallPointRec; fOutputallPointRec=0x0; }
  if (fOutputallPointSkip)                     { delete fOutputallPointSkip; fOutputallPointSkip=0x0; }
  if (fOutputpartPointRec)                     { delete fOutputpartPointRec; fOutputpartPointRec=0x0; }
  if (fOutputpartPointSkip)                  { delete fOutputpartPointSkip; fOutputpartPointSkip=0x0; }
  if (fOutputonepointSPDRec)                 { delete fOutputonepointSPDRec;fOutputonepointSPDRec=0x0;}
  if (fOutputonepointSPDSkip)              { delete fOutputonepointSPDSkip;fOutputonepointSPDSkip=0x0;}
  if (fOutputpostvTracRec)                      { delete fOutputpostvTracRec; fOutputpostvTracRec=0x0;}
  if (fOutputpostvTracSkip)                  {  delete fOutputpostvTracSkip; fOutputpostvTracSkip=0x0;}
  if (fOutputnegtvTracRec)                      { delete fOutputnegtvTracRec; fOutputnegtvTracRec=0x0;}
  if (fOutputnegtvTracSkip)                   { delete fOutputnegtvTracSkip; fOutputnegtvTracSkip=0x0;}
  if (fOutputpullAllpointRec)              {delete fOutputpullAllpointRec; fOutputpullAllpointRec=0x0;}
  if (fOutputpullAllpointSkip)           {delete fOutputpullAllpointSkip; fOutputpullAllpointSkip=0x0;}
  if (fOutputOnlyRefitRec)                       {delete fOutputOnlyRefitRec; fOutputOnlyRefitRec=0x0;}
  if (fOutputOnlyRefitSkip)                    {delete fOutputOnlyRefitSkip; fOutputOnlyRefitSkip=0x0;}
  if (fOutputSinThetaRec)                          {delete fOutputSinThetaRec; fOutputSinThetaRec=0x0;}  
  if (fOutputSinThetaSkip)                       {delete fOutputSinThetaSkip; fOutputSinThetaSkip=0x0;}
  if (fOutputallPointTrue)                       {delete fOutputallPointTrue; fOutputallPointTrue=0x0;}
  if (fOutputpostvTracTrue)                     {delete fOutputpostvTracTrue;fOutputpostvTracTrue=0x0;}
  if (fOutputnegtvTracTrue)                     {delete fOutputnegtvTracTrue;fOutputnegtvTracTrue=0x0;}
  if (fOutputpullAllpointTrue)            {delete fOutputpullAllpointTrue;fOutputpullAllpointTrue=0x0;}
  if (fOutputphiAllpointSkip)               {delete fOutputphiAllpointSkip;fOutputphiAllpointSkip=0x0;}
  if (fOutputphiPostvtracSkip)            {delete fOutputphiPostvtracSkip;fOutputphiPostvtracSkip=0x0;}
  if (fOutputphiNegtvtracSkip)            {delete fOutputphiNegtvtracSkip;fOutputphiNegtvtracSkip=0x0;}
  if (fOutputparticlePID)                           {delete fOutputparticlePID;fOutputparticlePID=0x0;}
  if (fOutputWithTrackCuts) {  delete fOutputWithTrackCuts; fOutputWithTrackCuts=0;}
  if (fOutputPt)                                                      {delete fOutputPt;fOutputPt=0x0;}
  if (fNentries)                                           { delete fNentries;     fNentries    =0x0; }
  if (fEstimVtx)                                           { delete fEstimVtx;     fEstimVtx    =0x0; }

} 
//______________________________________________________________________________________________________
void AliAnalysisTaskSEImpParRes::UserCreateOutputObjects()
{
  // 
  // Create the output container
  //
  
  if(fDebug>1) printf("AnalysisTaskSEImpParRes::UserCreateOutputObjects() \n");
  
  // Several histograms are more conveniently managed in a TList
  if (!fOutputitspureSARec) {
    fOutputitspureSARec = new TList();
    fOutputitspureSARec->SetOwner();
    fOutputitspureSARec->SetName("ITSpureSARec");
  }

  if (!fOutputitspureSASkip) {
    fOutputitspureSASkip = new TList();
    fOutputitspureSASkip->SetOwner();
    fOutputitspureSASkip->SetName("ITSpureSASkip");
  }

  if (!fOutputallPointRec) {
    fOutputallPointRec = new TList();
    fOutputallPointRec->SetOwner();
    fOutputallPointRec->SetName("allpointRec");
  }

  if (!fOutputallPointSkip) {
    fOutputallPointSkip = new TList();
    fOutputallPointSkip->SetOwner();
    fOutputallPointSkip->SetName("allpointSkip");
  }

  if (!fOutputpartPointRec) {
    fOutputpartPointRec = new TList();
    fOutputpartPointRec->SetOwner();
    fOutputpartPointRec->SetName("partpointRec");
  }

  if (!fOutputpartPointSkip) {
    fOutputpartPointSkip = new TList();
    fOutputpartPointSkip->SetOwner();
    fOutputpartPointSkip->SetName("partpointSkip");
  }

  if (!fOutputonepointSPDRec) {
    fOutputonepointSPDRec = new TList();
    fOutputonepointSPDRec->SetOwner();
    fOutputonepointSPDRec->SetName("onepointSPDRec");
  }

  if (!fOutputonepointSPDSkip) {
    fOutputonepointSPDSkip = new TList();
    fOutputonepointSPDSkip->SetOwner();
    fOutputonepointSPDSkip->SetName("onepointSPDSkip");
  }

  if (!fOutputpostvTracRec) {
    fOutputpostvTracRec = new TList();
    fOutputpostvTracRec->SetOwner();
    fOutputpostvTracRec->SetName("postvtracRec");
  }

  if (!fOutputpostvTracSkip) {
    fOutputpostvTracSkip = new TList();
    fOutputpostvTracSkip->SetOwner();
    fOutputpostvTracSkip->SetName("postvtracSkip");
  }
 
  if (!fOutputnegtvTracRec) {
    fOutputnegtvTracRec = new TList();
    fOutputnegtvTracRec->SetOwner();
    fOutputnegtvTracRec->SetName("negtvtracRe");
  }

  if (!fOutputnegtvTracSkip) {
    fOutputnegtvTracSkip = new TList();
    fOutputnegtvTracSkip->SetOwner();
    fOutputnegtvTracSkip->SetName("negtvtracSkip");
  }
  
  if (!fOutputpullAllpointSkip) {
    fOutputpullAllpointSkip = new TList();
    fOutputpullAllpointSkip->SetOwner();
    fOutputpullAllpointSkip->SetName("pullAllpointSkip");
  }
  
  if (!fOutputpullAllpointRec) {
    fOutputpullAllpointRec = new TList();
    fOutputpullAllpointRec->SetOwner();
    fOutputpullAllpointRec->SetName("pullAllpointRec");
  }
  
  if (!fOutputOnlyRefitRec) {
    fOutputOnlyRefitRec = new TList();
    fOutputOnlyRefitRec->SetOwner();
    fOutputOnlyRefitRec->SetName("onlyRefitRec");
  }
  
  if (!fOutputOnlyRefitSkip) {
    fOutputOnlyRefitSkip = new TList();
    fOutputOnlyRefitSkip->SetOwner();
    fOutputOnlyRefitSkip->SetName("onlyRefitRec");
  }
  
  if (!fOutputallPointTrue) {
    fOutputallPointTrue = new TList();
    fOutputallPointTrue->SetOwner();
    fOutputallPointTrue->SetName("allpointTrue");
  }
  
  if (!fOutputpostvTracTrue) {
    fOutputpostvTracTrue = new TList();
    fOutputpostvTracTrue->SetOwner();
    fOutputpostvTracTrue->SetName("postvtracTrue");
  }
  
  if (!fOutputnegtvTracTrue) {
    fOutputnegtvTracTrue = new TList();
    fOutputnegtvTracTrue->SetOwner();
    fOutputnegtvTracTrue->SetName("negtvtracTrue");
  }
  
  if (!fOutputpullAllpointTrue) {
    fOutputpullAllpointTrue = new TList();
    fOutputpullAllpointTrue->SetOwner();
    fOutputpullAllpointTrue->SetName("pullAllpointTrue");
  } 
  
  
  if (!fOutputparticlePID) {
    fOutputparticlePID = new TList();
    fOutputparticlePID->SetOwner();
    fOutputparticlePID->SetName("particlePID");
  }
 
  if (!fOutputPt) {
    fOutputPt = new TList();
    fOutputPt->SetOwner();
    fOutputPt->SetName("Pt");
  }

  if (!fOutputWithTrackCuts) {
    fOutputWithTrackCuts = new TList();
    fOutputWithTrackCuts->SetOwner();
    fOutputWithTrackCuts->SetName("OutputWithESDTrackCuts");
  }

  const Int_t nhist=26;
  const TString d0ITSpureSArphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0ITSpureSAzTitle = "d_{0} Distribution_z; d_{0} [#mum]; Entries";
  const TString d0allpointrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0allpointzTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";
  const TString d0partpointrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0partpointzTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";
  const TString d0onepointSPDrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0onepointSPDzTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0postvtracrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0postvtraczTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";
  const TString d0negtvtracrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0negtvtraczTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";
  const TString d0pullAllpointrphiTitle = "d_{0} Pull Distribution_rphi; d_{0} pull; Entries";
  const TString d0pullAllpointzTitle  = "d_{0} Pull Distribution_z; d_{0} pull; Entries";
  const TString d0onlyRefitrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0onlyRefitzTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";
  const TString d0ptTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0rphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0zTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";
  const TString d0rphiParticlPID = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0zPrtilePID  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";

  TString  named0itspureSArphiRec,named0itspureSAzRec,named0allpointrphiRec, named0allpointzRec,named0partpointrphiRec, named0partpointzRec,named0onepointSPDrphiRec, named0onepointSPDzRec,named0postvtracrphiRec, named0postvtraczRec,named0negtvtracrphiRec, named0negtvtraczRec,named0pt,named0pullAllpointrphiRec,named0pullAllpointzRec,named0onlyRefitrphiRec,named0onlyRefitzRec,named0pionPIDrphiRec, named0pionPIDzRec,named0kaonPIDrphiRec, named0kaonPIDzRec,named0protonPIDrphiRec, named0protonPIDzRec;
 
  TH1F *d0ITSpureSArphiRec=0,*d0ITSpureSAzRec=0,*d0AllpointrphiRec=0, *d0AllpointzRec=0,*d0PartpointrphiRec=0, *d0PartpointzRec=0,
    *d0OnepointSPDrphiRec=0,*d0OnepointSPDzRec=0,*d0PostvtracrphiRec=0, *d0PostvtraczRec=0,*d0NegtvtracrphiRec=0, *d0NegtvtraczRec=0,*d0Pt=0,*d0PullAllpointrphiRec=0,*d0PullAllpointzRec=0,*d0OnlyRefitrphiRec=0,*d0OnlyRefitzRec=0,*d0PionPIDrphiRec=0,*d0PionPIDzRec=0,*d0KaonPIDrphiRec=0,*d0KaonPIDzRec=0,*d0ProtonPIDrphiRec=0,*d0ProtonPIDzRec=0;

  TString  named0itspureSArphiSkip,named0itspureSAzSkip,named0allpointrphiSkip, named0allpointzSkip,named0partpointrphiSkip, named0partpointzSkip,named0onepointSPDrphiSkip, named0onepointSPDzSkip,named0postvtracrphiSkip, named0postvtraczSkip,named0negtvtracrphiSkip, named0negtvtraczSkip,named0ptSkip,named0pullAllpointrphiSkip,named0pullAllpointzSkip,named0onlyRefitrphiSkip,named0onlyRefitzSkip,named0allpointrphiTrue, named0allpointzTrue,named0postvtracrphiTrue, named0postvtraczTrue,named0negtvtracrphiTrue, named0negtvtraczTrue,named0pullAllpointrphiTrue,named0pullAllpointzTrue,named0pionPIDrphiSkip, named0pionPIDzSkip,named0kaonPIDrphiSkip, named0kaonPIDzSkip,named0protonPIDrphiSkip, named0protonPIDzSkip;

  TH1F *d0ITSpureSArphiSkip=0,*d0ITSpureSAzSkip=0,*d0AllpointrphiSkip=0, *d0AllpointzSkip=0,*d0PartpointrphiSkip=0, *d0PartpointzSkip=0,*d0OnepointSPDrphiSkip=0,*d0OnepointSPDzSkip=0,*d0PostvtracrphiSkip=0, *d0PostvtraczSkip=0,*d0NegtvtracrphiSkip=0,*d0NegtvtraczSkip=0,*d0PullAllpointrphiSkip=0,*d0PullAllpointzSkip=0,*d0OnlyRefitrphiSkip=0,*d0OnlyRefitzSkip=0,*d0AllpointrphiTrue=0, *d0AllpointzTrue=0,*d0PostvtracrphiTrue=0, *d0PostvtraczTrue=0,*d0NegtvtracrphiTrue=0,*d0NegtvtraczTrue=0,*d0PullAllpointrphiTrue,*d0PullAllpointzTrue,*d0PionPIDrphiSkip=0,*d0PionPIDzSkip=0,*d0KaonPIDrphiSkip=0,*d0KaonPIDzSkip=0,*d0ProtonPIDrphiSkip=0,*d0ProtonPIDzSkip=0;

  TString named0DistrESDTCrphiRec, named0DistrESDTCrphiSkip, named0DistrESDTCrphiTrue, named0DistrESDTCzRec, named0DistrESDTCzSkip, named0DistrESDTCzTrue, named0PullESDTCrphiRec, named0PullESDTCrphiSkip, named0PullESDTCrphiTrue, named0PullESDTCzRec, named0PullESDTCzSkip, named0PullESDTCzTrue, named0ptESDTC;

  TH1F *d0DistrESDTCrphiRec=0, *d0DistrESDTCrphiSkip=0, *d0DistrESDTCrphiTrue=0, *d0DistrESDTCzRec=0, *d0DistrESDTCzSkip=0, *d0DistrESDTCzTrue=0, *d0PullESDTCrphiRec=0, *d0PullESDTCrphiSkip=0, *d0PullESDTCrphiTrue=0, *d0PullESDTCzRec=0, *d0PullESDTCzSkip=0, *d0PullESDTCzTrue=0, *d0PtESDTC;

  for(Int_t i=1; i<=nhist; i++) {
   
    named0itspureSArphiRec = "d0itspureSArphiRec_";
    named0itspureSArphiRec += i;
    named0itspureSAzRec = "d0itspureSAzRec_";
    named0itspureSAzRec += i;
    d0ITSpureSArphiRec = new TH1F(named0itspureSArphiRec.Data(), d0ITSpureSArphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0ITSpureSArphiRec->Sumw2();
    d0ITSpureSArphiRec->SetMinimum(0);  
    fOutputitspureSARec->Add(d0ITSpureSArphiRec);
    d0ITSpureSAzRec = new TH1F(named0itspureSAzRec.Data(), d0ITSpureSAzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0ITSpureSAzRec->Sumw2();
    d0ITSpureSAzRec->SetMinimum(0);  
    fOutputitspureSARec->Add(d0ITSpureSAzRec);

    named0itspureSArphiSkip = "d0itspureSArphiSkip_";
    named0itspureSArphiSkip += i;
    named0itspureSAzSkip = "d0itspureSAzSkip_";
    named0itspureSAzSkip += i;
    d0ITSpureSArphiSkip = new TH1F(named0itspureSArphiSkip.Data(), d0ITSpureSArphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));    d0ITSpureSArphiSkip->Sumw2();
    d0ITSpureSArphiSkip->SetMinimum(0);  
    fOutputitspureSASkip->Add(d0ITSpureSArphiSkip);
    d0ITSpureSAzSkip = new TH1F(named0itspureSAzSkip.Data(), d0ITSpureSAzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0ITSpureSAzSkip->Sumw2();
    d0ITSpureSAzSkip->SetMinimum(0);  
    fOutputitspureSASkip->Add(d0ITSpureSAzSkip);

    named0allpointrphiRec = "d0allpointrphiRec_";
    named0allpointrphiRec += i;
    named0allpointzRec = "d0allpointzRec_";
    named0allpointzRec += i;
    d0AllpointrphiRec = new TH1F(named0allpointrphiRec.Data(), d0allpointrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0AllpointrphiRec->Sumw2();
    d0AllpointrphiRec->SetMinimum(0);  
    fOutputallPointRec->Add(d0AllpointrphiRec);
    d0AllpointzRec= new TH1F(named0allpointzRec.Data(), d0allpointzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0AllpointzRec->Sumw2();
    d0AllpointzRec->SetMinimum(0);  
    fOutputallPointRec->Add(d0AllpointzRec);

    named0allpointrphiSkip = "d0allpointrphiSkip_";
    named0allpointrphiSkip += i;
    named0allpointzSkip = "d0allpointzSkip_";
    named0allpointzSkip += i;
    d0AllpointrphiSkip = new TH1F(named0allpointrphiSkip.Data(), d0allpointrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0AllpointrphiSkip->Sumw2();
    d0AllpointrphiSkip->SetMinimum(0);  
    fOutputallPointSkip->Add(d0AllpointrphiSkip);
    d0AllpointzSkip = new TH1F(named0allpointzSkip.Data(), d0allpointzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0AllpointzSkip->Sumw2();
    d0AllpointzSkip->SetMinimum(0);  
    fOutputallPointSkip->Add(d0AllpointzSkip);

    named0partpointrphiRec = "d0partpointrphiRec_";
    named0partpointrphiRec += i;
    named0partpointzRec = "d0partpointzRec_";
    named0partpointzRec += i;
    d0PartpointrphiRec = new TH1F(named0partpointrphiRec.Data(), d0partpointrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PartpointrphiRec->Sumw2();
    d0PartpointrphiRec->SetMinimum(0);  
    fOutputpartPointRec->Add(d0PartpointrphiRec);
    d0PartpointzRec = new TH1F(named0partpointzRec.Data(), d0partpointzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PartpointzRec->Sumw2();
    d0PartpointzRec->SetMinimum(0);  
    fOutputpartPointRec->Add(d0PartpointzRec);

    named0partpointrphiSkip = "d0partpointrphiSkip_";
    named0partpointrphiSkip += i;
    named0partpointzSkip = "d0partpointzSkip_";
    named0partpointzSkip += i;
    d0PartpointrphiSkip = new TH1F(named0partpointrphiSkip.Data(), d0partpointrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PartpointrphiSkip->Sumw2();
    d0PartpointrphiSkip->SetMinimum(0);  
    fOutputpartPointSkip->Add(d0PartpointrphiSkip);
    d0PartpointzSkip = new TH1F(named0partpointzSkip.Data(), d0partpointzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PartpointzSkip->Sumw2();
    d0PartpointzSkip->SetMinimum(0);  
    fOutputpartPointSkip->Add(d0PartpointzSkip);

    named0onepointSPDrphiRec = "d0onepointSPDrphiRec_";
    named0onepointSPDrphiRec += i;
    named0onepointSPDzRec = "d0onepointSPDzRec_";
    named0onepointSPDzRec += i;
    d0OnepointSPDrphiRec = new TH1F(named0onepointSPDrphiRec.Data(), d0onepointSPDrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0OnepointSPDrphiRec->Sumw2();
    d0OnepointSPDrphiRec->SetMinimum(0);  
    fOutputonepointSPDRec->Add(d0OnepointSPDrphiRec);
    d0OnepointSPDzRec = new TH1F(named0onepointSPDzRec.Data(), d0onepointSPDzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0OnepointSPDzRec->Sumw2();
    d0OnepointSPDzRec->SetMinimum(0);  
    fOutputonepointSPDRec->Add(d0OnepointSPDzRec);

    named0onepointSPDrphiSkip = "d0onepointSPDrphiSkip_";
    named0onepointSPDrphiSkip += i;
    named0onepointSPDzSkip = "d0onepointSPDzSkip_";
    named0onepointSPDzSkip += i;
    d0OnepointSPDrphiSkip = new TH1F(named0onepointSPDrphiSkip.Data(), d0onepointSPDrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0OnepointSPDrphiSkip->Sumw2();
    d0OnepointSPDrphiSkip->SetMinimum(0);  
    fOutputonepointSPDSkip->Add(d0OnepointSPDrphiSkip);
    d0OnepointSPDzSkip = new TH1F(named0onepointSPDzSkip.Data(), d0onepointSPDzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0OnepointSPDzSkip->Sumw2();
    d0OnepointSPDzSkip->SetMinimum(0);  
    fOutputonepointSPDSkip->Add(d0OnepointSPDzSkip);

    named0postvtracrphiRec = "d0postvtracrphiRec_";
    named0postvtracrphiRec += i;
    named0postvtraczRec = "d0postvtraczRec_";
    named0postvtraczRec += i;
    d0PostvtracrphiRec = new TH1F(named0postvtracrphiRec.Data(), d0postvtracrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PostvtracrphiRec->Sumw2();
    d0PostvtracrphiRec->SetMinimum(0);  
    fOutputpostvTracRec->Add(d0PostvtracrphiRec);
    d0PostvtraczRec = new TH1F(named0postvtraczRec.Data(), d0postvtraczTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PostvtraczRec->Sumw2();
    d0PostvtraczRec->SetMinimum(0);  
    fOutputpostvTracRec->Add(d0PostvtraczRec);

    named0postvtracrphiSkip = "d0postvtracrphiSkip_";
    named0postvtracrphiSkip += i;
    named0postvtraczSkip = "d0postvtraczSkip_";
    named0postvtraczSkip += i;
    d0PostvtracrphiSkip = new TH1F(named0postvtracrphiSkip.Data(), d0postvtracrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PostvtracrphiSkip->Sumw2();
    d0PostvtracrphiSkip->SetMinimum(0);  
    fOutputpostvTracSkip->Add(d0PostvtracrphiSkip);
    d0PostvtraczSkip = new TH1F(named0postvtraczSkip.Data(), d0postvtraczTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PostvtraczSkip->Sumw2();
    d0PostvtraczSkip->SetMinimum(0);  
    fOutputpostvTracSkip->Add(d0PostvtraczSkip);

    named0negtvtracrphiRec = "d0negtvtracrphiRec_";
    named0negtvtracrphiRec += i;
    named0negtvtraczRec = "d0negtvtraczRec_";
    named0negtvtraczRec += i;
    d0NegtvtracrphiRec = new TH1F(named0negtvtracrphiRec.Data(), d0negtvtracrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0NegtvtracrphiRec->Sumw2();
    d0NegtvtracrphiRec->SetMinimum(0);  
    fOutputnegtvTracRec->Add(d0NegtvtracrphiRec);
    d0NegtvtraczRec = new TH1F(named0negtvtraczRec.Data(), d0negtvtraczTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0NegtvtraczRec->Sumw2();
    d0NegtvtraczRec->SetMinimum(0);  
    fOutputnegtvTracRec->Add(d0NegtvtraczRec);

    named0negtvtracrphiSkip = "d0negtvtracrphiSkip_";
    named0negtvtracrphiSkip += i;
    named0negtvtraczSkip = "d0negtvtraczSkip_";
    named0negtvtraczSkip += i;
    d0NegtvtracrphiSkip = new TH1F(named0negtvtracrphiSkip.Data(), d0negtvtracrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0NegtvtracrphiSkip->Sumw2();
    d0NegtvtracrphiSkip->SetMinimum(0);  
    fOutputnegtvTracSkip->Add(d0NegtvtracrphiSkip);
    d0NegtvtraczSkip = new TH1F(named0negtvtraczSkip.Data(), d0negtvtraczTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0NegtvtraczSkip->Sumw2();
    d0NegtvtraczSkip->SetMinimum(0);  
    fOutputnegtvTracSkip->Add(d0NegtvtraczSkip);

    named0pullAllpointrphiSkip = "d0pullAllpointrphiSkip_";
    named0pullAllpointrphiSkip +=i;
    named0pullAllpointzSkip = "d0pullAllpointzSkip_";
    named0pullAllpointzSkip +=i;
    d0PullAllpointrphiSkip = new TH1F(named0pullAllpointrphiSkip.Data(),d0pullAllpointrphiTitle.Data(),400,-10.,10.);
    d0PullAllpointrphiSkip->Sumw2();
    d0PullAllpointrphiSkip->SetMinimum(0);
    fOutputpullAllpointSkip->Add(d0PullAllpointrphiSkip);
    d0PullAllpointzSkip = new TH1F(named0pullAllpointzSkip.Data(),d0pullAllpointzTitle.Data(),400,-10.,10.);
    d0PullAllpointzSkip->Sumw2();
    d0PullAllpointzSkip->SetMinimum(0);
    fOutputpullAllpointSkip->Add(d0PullAllpointzSkip);

    named0pullAllpointrphiRec = "d0pullAllpointrphiRec_";
    named0pullAllpointrphiRec +=i;
    named0pullAllpointzRec = "d0pullAllpointzRec_";
    named0pullAllpointzRec +=i;
    d0PullAllpointrphiRec = new TH1F(named0pullAllpointrphiRec.Data(),d0pullAllpointrphiTitle.Data(),400,-10.,10.);
    d0PullAllpointrphiRec->Sumw2();
    d0PullAllpointrphiRec->SetMinimum(0);
    fOutputpullAllpointRec->Add(d0PullAllpointrphiRec);
    d0PullAllpointzRec = new TH1F(named0pullAllpointzRec.Data(),d0pullAllpointzTitle.Data(),400,-10.,10.);
    d0PullAllpointzRec->Sumw2();
    d0PullAllpointzRec->SetMinimum(0);
    fOutputpullAllpointRec->Add(d0PullAllpointzRec);

    named0onlyRefitrphiRec = "d0onlyrefitrphiRec_";
    named0onlyRefitrphiRec +=i;
    named0onlyRefitzRec = "d0onlyrefitzRec_";
    named0onlyRefitzRec +=i;
    d0OnlyRefitrphiRec = new TH1F(named0onlyRefitrphiRec.Data(),d0onlyRefitrphiTitle.Data(),400,-Getd0HistRange(i),Getd0HistRange(i));
    d0OnlyRefitrphiRec->Sumw2();
    d0OnlyRefitrphiRec->SetMinimum(0);
    fOutputOnlyRefitRec->Add(d0OnlyRefitrphiRec);
    d0OnlyRefitzRec = new TH1F(named0onlyRefitzRec.Data(),d0onlyRefitzTitle.Data(),400,-Getd0HistRange(i),Getd0HistRange(i));
    d0OnlyRefitzRec->Sumw2();
    d0OnlyRefitzRec->SetMinimum(0);
    fOutputOnlyRefitRec->Add(d0OnlyRefitzRec);

    named0onlyRefitrphiSkip = "d0onlyrefitrphiSkip_";
    named0onlyRefitrphiSkip +=i;
    named0onlyRefitzSkip = "d0onlyrefitzSkip_";
    named0onlyRefitzSkip +=i;
    d0OnlyRefitrphiSkip = new TH1F(named0onlyRefitrphiSkip.Data(),d0onlyRefitrphiTitle.Data(),400,-Getd0HistRange(i),Getd0HistRange(i));
    d0OnlyRefitrphiSkip->Sumw2();
    d0OnlyRefitrphiSkip->SetMinimum(0);
    fOutputOnlyRefitSkip->Add(d0OnlyRefitrphiSkip);
    d0OnlyRefitzSkip = new TH1F(named0onlyRefitzSkip.Data(),d0onlyRefitzTitle.Data(),400,-Getd0HistRange(i),Getd0HistRange(i));
    d0OnlyRefitzSkip->Sumw2();
    d0OnlyRefitzSkip->SetMinimum(0);
    fOutputOnlyRefitSkip->Add(d0OnlyRefitzSkip);

    named0allpointrphiTrue = "d0allpointrphiTrue_";
    named0allpointrphiTrue += i;
    named0allpointzTrue = "d0allpointzTrue_";
    named0allpointzTrue += i;
    d0AllpointrphiTrue = new TH1F(named0allpointrphiTrue.Data(), d0allpointrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0AllpointrphiTrue->Sumw2();
    d0AllpointrphiTrue->SetMinimum(0);  
    fOutputallPointTrue->Add(d0AllpointrphiTrue);
    d0AllpointzTrue = new TH1F(named0allpointzTrue.Data(), d0allpointzTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0AllpointzTrue->Sumw2();
    d0AllpointzTrue->SetMinimum(0);  
    fOutputallPointTrue->Add(d0AllpointzTrue);

    named0postvtracrphiTrue = "d0postvtracrphiTrue_";
    named0postvtracrphiTrue += i;
    named0postvtraczTrue = "d0postvtraczTrue_";
    named0postvtraczTrue += i;
    d0PostvtracrphiTrue = new TH1F(named0postvtracrphiTrue.Data(), d0postvtracrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PostvtracrphiTrue->Sumw2();
    d0PostvtracrphiTrue->SetMinimum(0);  
    fOutputpostvTracTrue->Add(d0PostvtracrphiTrue);
    d0PostvtraczTrue = new TH1F(named0postvtraczTrue.Data(), d0postvtraczTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PostvtraczTrue->Sumw2();
    d0PostvtraczTrue->SetMinimum(0);  
    fOutputpostvTracTrue->Add(d0PostvtraczTrue);

    named0negtvtracrphiTrue = "d0negtvtracrphiTrue_";
    named0negtvtracrphiTrue += i;
    named0negtvtraczTrue = "d0negtvtraczTrue_";
    named0negtvtraczTrue += i;
    d0NegtvtracrphiTrue = new TH1F(named0negtvtracrphiTrue.Data(), d0negtvtracrphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0NegtvtracrphiTrue->Sumw2();
    d0NegtvtracrphiTrue->SetMinimum(0);  
    fOutputnegtvTracTrue->Add(d0NegtvtracrphiTrue);
    d0NegtvtraczTrue = new TH1F(named0negtvtraczTrue.Data(), d0negtvtraczTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0NegtvtraczTrue->Sumw2();
    d0NegtvtraczTrue->SetMinimum(0);  
    fOutputnegtvTracTrue->Add(d0NegtvtraczTrue);

    named0pullAllpointrphiTrue = "d0pullAllpointrphiTrue_";
    named0pullAllpointrphiTrue +=i;
    named0pullAllpointzTrue = "d0pullAllpointzTrue_";
    named0pullAllpointzTrue +=i;
    d0PullAllpointrphiTrue = new TH1F(named0pullAllpointrphiTrue.Data(),d0pullAllpointrphiTitle.Data(),400,-10.,10.);
    d0PullAllpointrphiTrue->Sumw2();
    d0PullAllpointrphiTrue->SetMinimum(0);
    fOutputpullAllpointTrue->Add(d0PullAllpointrphiTrue);
    d0PullAllpointzTrue = new TH1F(named0pullAllpointzTrue.Data(),d0pullAllpointzTitle.Data(),400,-10.,10.);
    d0PullAllpointzTrue->Sumw2();
    d0PullAllpointzTrue->SetMinimum(0);
    fOutputpullAllpointTrue->Add(d0PullAllpointzTrue);


    named0pionPIDrphiRec = "d0pionPIDrphiRec_";
    named0pionPIDrphiRec += i;
    named0pionPIDzRec = "d0pionPIDzRec_";
    named0pionPIDzRec += i;
    d0PionPIDrphiRec = new TH1F(named0pionPIDrphiRec.Data(), d0rphiParticlPID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PionPIDrphiRec->Sumw2();
    d0PionPIDrphiRec->SetMinimum(0);  
    fOutputparticlePID->Add(d0PionPIDrphiRec);
    d0PionPIDzRec = new TH1F(named0pionPIDzRec.Data(), d0zPrtilePID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PionPIDzRec->Sumw2();
    d0PionPIDzRec->SetMinimum(0);  
    fOutputparticlePID->Add(d0PionPIDzRec);

    named0pionPIDrphiSkip = "d0pionPIDrphiSkip_";
    named0pionPIDrphiSkip += i;
    named0pionPIDzSkip = "d0pionPIDzSkip_";
    named0pionPIDzSkip += i;
    d0PionPIDrphiSkip = new TH1F(named0pionPIDrphiSkip.Data(), d0rphiParticlPID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PionPIDrphiSkip->Sumw2();
    d0PionPIDrphiSkip->SetMinimum(0);  
    fOutputparticlePID->Add(d0PionPIDrphiSkip);
    d0PionPIDzSkip = new TH1F(named0pionPIDzSkip.Data(), d0zPrtilePID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0PionPIDzSkip->Sumw2();
    d0PionPIDzSkip->SetMinimum(0);  
    fOutputparticlePID->Add(d0PionPIDzSkip);

    named0kaonPIDrphiRec = "d0kaonPIDrphiRec_";
    named0kaonPIDrphiRec += i;
    named0kaonPIDzRec = "d0kaonPIDzRec_";
    named0kaonPIDzRec += i;
    d0KaonPIDrphiRec = new TH1F(named0kaonPIDrphiRec.Data(), d0rphiParticlPID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0KaonPIDrphiRec->Sumw2();
    d0KaonPIDrphiRec->SetMinimum(0);  
    fOutputparticlePID->Add(d0KaonPIDrphiRec);
    d0KaonPIDzRec = new TH1F(named0kaonPIDzRec.Data(), d0zPrtilePID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0KaonPIDzRec->Sumw2();
    d0KaonPIDzRec->SetMinimum(0);  
    fOutputparticlePID->Add(d0KaonPIDzRec);

    named0kaonPIDrphiSkip = "d0kaonPIDrphiSkip_";
    named0kaonPIDrphiSkip += i;
    named0kaonPIDzSkip = "d0kaonPIDzSkip_";
    named0kaonPIDzSkip += i;
    d0KaonPIDrphiSkip = new TH1F(named0kaonPIDrphiSkip.Data(), d0rphiParticlPID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0KaonPIDrphiSkip->Sumw2();
    d0KaonPIDrphiSkip->SetMinimum(0);  
    fOutputparticlePID->Add(d0KaonPIDrphiSkip);
    d0KaonPIDzSkip = new TH1F(named0kaonPIDzSkip.Data(), d0zPrtilePID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0KaonPIDzSkip->Sumw2();
    d0KaonPIDzSkip->SetMinimum(0);  
    fOutputparticlePID->Add(d0KaonPIDzSkip);

    named0protonPIDrphiRec = "d0protonPIDrphiRec_";
    named0protonPIDrphiRec += i;
    named0protonPIDzRec = "d0protonPIDzRec_";
    named0protonPIDzRec += i;
    d0ProtonPIDrphiRec = new TH1F(named0protonPIDrphiRec.Data(), d0rphiParticlPID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0ProtonPIDrphiRec->Sumw2();
    d0ProtonPIDrphiRec->SetMinimum(0);  
    fOutputparticlePID->Add(d0ProtonPIDrphiRec);
    d0ProtonPIDzRec = new TH1F(named0protonPIDzRec.Data(), d0zPrtilePID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0ProtonPIDzRec->Sumw2();
    d0ProtonPIDzRec->SetMinimum(0);  
    fOutputparticlePID->Add(d0ProtonPIDzRec);

    named0protonPIDrphiSkip = "d0protonPIDrphiSkip_";
    named0protonPIDrphiSkip += i;
    named0protonPIDzSkip = "d0protonPIDzSkip_";
    named0protonPIDzSkip += i;
    d0ProtonPIDrphiSkip = new TH1F(named0protonPIDrphiSkip.Data(), d0rphiParticlPID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0ProtonPIDrphiSkip->Sumw2();
    d0ProtonPIDrphiSkip->SetMinimum(0);  
    fOutputparticlePID->Add(d0ProtonPIDrphiSkip);
    d0ProtonPIDzSkip = new TH1F(named0protonPIDzSkip.Data(), d0zPrtilePID.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0ProtonPIDzSkip->Sumw2();
    d0ProtonPIDzSkip->SetMinimum(0);  
    fOutputparticlePID->Add(d0ProtonPIDzSkip);

    named0pt = "d0pt_";
    named0pt += i;
    d0Pt = new TH1F(named0pt.Data(), d0ptTitle.Data(), 100, 0, 35.);
    d0Pt->Sumw2();
    d0Pt->SetMinimum(0);  
    fOutputPt->Add(d0Pt);

    //foutputwithtrackcuts
    named0DistrESDTCrphiRec = "d0DistrESDTCrphiRec_";
    named0DistrESDTCrphiRec +=i;
    named0DistrESDTCzRec = "d0DistrESDTCzRec_";
    named0DistrESDTCzRec +=i;
    d0DistrESDTCrphiRec = new TH1F(named0DistrESDTCrphiRec.Data(),d0allpointrphiTitle.Data(),400, -Getd0HistRange(i), Getd0HistRange(i));
    d0DistrESDTCrphiRec->Sumw2();
    d0DistrESDTCrphiRec->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0DistrESDTCrphiRec);
    d0DistrESDTCzRec = new TH1F(named0DistrESDTCzRec.Data(),d0allpointzTitle.Data(),400, -Getd0HistRange(i), Getd0HistRange(i));
    d0DistrESDTCzRec->Sumw2();
    d0DistrESDTCzRec->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0DistrESDTCzRec);

    named0DistrESDTCrphiSkip = "d0DistrESDTCrphiSkip_";
    named0DistrESDTCrphiSkip +=i;
    named0DistrESDTCzSkip = "d0DistrESDTCzSkip_";
    named0DistrESDTCzSkip +=i;
    d0DistrESDTCrphiSkip = new TH1F(named0DistrESDTCrphiSkip.Data(),d0allpointrphiTitle.Data(),400, -Getd0HistRange(i), Getd0HistRange(i));
    d0DistrESDTCrphiSkip->Sumw2();
    d0DistrESDTCrphiSkip->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0DistrESDTCrphiSkip);
    d0DistrESDTCzSkip = new TH1F(named0DistrESDTCzSkip.Data(),d0allpointzTitle.Data(),400, -Getd0HistRange(i), Getd0HistRange(i));
    d0DistrESDTCzSkip->Sumw2();
    d0DistrESDTCzSkip->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0DistrESDTCzSkip);

    named0DistrESDTCrphiTrue = "d0DistrESDTCrphiTrue_";
    named0DistrESDTCrphiTrue +=i;
    named0DistrESDTCzTrue = "d0DistrESDTCzTrue_";
    named0DistrESDTCzTrue +=i;
    d0DistrESDTCrphiTrue = new TH1F(named0DistrESDTCrphiTrue.Data(),d0allpointrphiTitle.Data(),400, -Getd0HistRange(i), Getd0HistRange(i));
    d0DistrESDTCrphiTrue->Sumw2();
    d0DistrESDTCrphiTrue->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0DistrESDTCrphiTrue);
    d0DistrESDTCzTrue = new TH1F(named0DistrESDTCzTrue.Data(),d0allpointzTitle.Data(),400, -Getd0HistRange(i), Getd0HistRange(i));
    d0DistrESDTCzTrue->Sumw2();
    d0DistrESDTCzTrue->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0DistrESDTCzTrue);

    named0PullESDTCrphiRec = "d0PullESDTCrphiRec_";
    named0PullESDTCrphiRec +=i;
    named0PullESDTCzRec = "d0PullESDTCzRec_";
    named0PullESDTCzRec +=i;
    d0PullESDTCrphiRec = new TH1F(named0PullESDTCrphiRec.Data(),d0pullAllpointrphiTitle.Data(),400,-10.,10.);
    d0PullESDTCrphiRec->Sumw2();
    d0PullESDTCrphiRec->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0PullESDTCrphiRec);
    d0PullESDTCzRec = new TH1F(named0PullESDTCzRec.Data(),d0pullAllpointzTitle.Data(),400,-10.,10.);
    d0PullESDTCzRec->Sumw2();
    d0PullESDTCzRec->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0PullESDTCzRec);

    named0PullESDTCrphiSkip = "d0PullESDTCrphiSkip_";
    named0PullESDTCrphiSkip +=i;
    named0PullESDTCzSkip = "d0PullESDTCzSkip_";
    named0PullESDTCzSkip +=i;
    d0PullESDTCrphiSkip = new TH1F(named0PullESDTCrphiSkip.Data(),d0pullAllpointrphiTitle.Data(),400,-10.,10.);
    d0PullESDTCrphiSkip->Sumw2();
    d0PullESDTCrphiSkip->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0PullESDTCrphiSkip);
    d0PullESDTCzSkip = new TH1F(named0PullESDTCzSkip.Data(),d0pullAllpointzTitle.Data(),400,-10.,10.);
    d0PullESDTCzSkip->Sumw2();
    d0PullESDTCzSkip->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0PullESDTCzSkip);

    named0PullESDTCrphiTrue = "d0PullESDTCrphiTrue_";
    named0PullESDTCrphiTrue +=i;
    named0PullESDTCzTrue = "d0PullESDTCzTrue_";
    named0PullESDTCzTrue +=i;
    d0PullESDTCrphiTrue = new TH1F(named0PullESDTCrphiTrue.Data(),d0pullAllpointrphiTitle.Data(),400,-10.,10.);
    d0PullESDTCrphiTrue->Sumw2();
    d0PullESDTCrphiTrue->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0PullESDTCrphiTrue);
    d0PullESDTCzTrue = new TH1F(named0PullESDTCzTrue.Data(),d0pullAllpointzTitle.Data(),400,-10.,10.);
    d0PullESDTCzTrue->Sumw2();
    d0PullESDTCzTrue->SetMinimum(0);
    fOutputWithTrackCuts->Add(d0PullESDTCzTrue);

    named0ptESDTC = "d0ptESDTC_";
    named0ptESDTC += i;
    d0PtESDTC = new TH1F(named0ptESDTC.Data(), d0ptTitle.Data(), 100, 0, 35.);
    d0PtESDTC->Sumw2();
    d0PtESDTC->SetMinimum(0);  
    fOutputWithTrackCuts->Add(d0PtESDTC);


  }

  
  if (!fOutputSinThetaRec){
    fOutputSinThetaRec = new TList();
    fOutputSinThetaRec->SetOwner();
    fOutputSinThetaRec->SetName("thetaRec");
  }
  
  if (!fOutputSinThetaSkip){
    fOutputSinThetaSkip = new TList();
    fOutputSinThetaSkip->SetOwner();
    fOutputSinThetaSkip->SetName("thetaSkip");
  }
  
  if (!fOutputphiAllpointSkip) {
    fOutputphiAllpointSkip = new TList();
    fOutputphiAllpointSkip->SetOwner();
    fOutputphiAllpointSkip->SetName("phiallpointSkip");
  }
  
  if (!fOutputphiPostvtracSkip) {
    fOutputphiPostvtracSkip = new TList();
    fOutputphiPostvtracSkip->SetOwner();
    fOutputphiPostvtracSkip->SetName("postvtracSkip");
  }
  
  if (!fOutputphiNegtvtracSkip) {
    fOutputphiNegtvtracSkip = new TList();
    fOutputphiNegtvtracSkip->SetOwner();
    fOutputphiNegtvtracSkip->SetName("negtvtracSkip");
  }
  
  
  const TString d0sinThetarphiTitle ="d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0sinThetazTitle ="d_{0} Distribution_z; d_{0} [#mum]; Entries";  
  const TString d0phiAllpointrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0phiAllpointzTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";    
  const TString  d0phiPostvtracrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0phiPostvtraczTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";  
  const TString  d0phiNegtvtracrphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0phiNegtvtraczTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";   
  TString named0sinThetaonerphiRec,named0sinThetaonezRec,named0sinThetatworphiRec,named0sinThetatwozRec,named0sinThetathreerphiRec,named0sinThetathreezRec,named0sinThetafourrphiRec,named0sinThetafourzRec,named0thetaForwardrphiRec,named0thetaForwardzRec,named0thetaBackwardrphiRec,named0thetaBackwardzRec;
  
  TH1F *d0SinThetaonerphiRec,*d0SinThetaonezRec,*d0SinThetatworphiRec,*d0SinThetatwozRec,*d0SinThetathreerphiRec,*d0SinThetathreezRec,*d0SinThetafourrphiRec,*d0SinThetafourzRec,*d0ThetaforwardrphiRec,*d0ThetaforwardzRec,*d0ThetabackwardrphiRec,*d0ThetabackwardzRec;
  
  TString  named0sinThetaonerphiSkip,named0sinThetaonezSkip,named0sinThetatworphiSkip,named0sinThetatwozSkip,named0sinThetathreerphiSkip,named0sinThetathreezSkip,named0sinThetafourrphiSkip,named0sinThetafourzSkip,named0phiAllpointrphiSkip, named0phiAllpointzSkip,named0phiPostvtracrphiSkip, named0phiPostvtraczSkip,named0phiNegtvtracrphiSkip,named0phiNegtvtraczSkip,named0thetaForwardrphiSkip,named0thetaForwardzSkip,named0thetaBackwardrphiSkip,named0thetaBackwardzSkip;
  
  TH1F*d0SinThetaonerphiSkip,*d0SinThetaonezSkip,*d0SinThetatworphiSkip,*d0SinThetatwozSkip,*d0SinThetathreerphiSkip,*d0SinThetathreezSkip,*d0SinThetafourrphiSkip,*d0SinThetafourzSkip, *d0PhiAllpointrphiSkip,*d0PhiAllpointzSkip,*d0PhiPostvtracrphiSkip,*d0PhiPostvtraczSkip,*d0PhiNegtvtracrphiSkip,*d0PhiNegtvtraczSkip,*d0ThetaforwardrphiSkip,*d0ThetaforwardzSkip,*d0ThetabackwardrphiSkip,*d0ThetabackwardzSkip;
  
  const Int_t nhistm=10;
  for(Int_t i=0; i<=nhistm; i++) {
    named0sinThetaonerphiRec = "d0sinthetaonerphiRec_";
    named0sinThetaonerphiRec += i;
    named0sinThetaonezRec ="d0sinthetaonezRec_";
    named0sinThetaonezRec += i;
    d0SinThetaonerphiRec = new TH1F(named0sinThetaonerphiRec.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0SinThetaonerphiRec->Sumw2();
    d0SinThetaonerphiRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0SinThetaonerphiRec);
    d0SinThetaonezRec = new TH1F(named0sinThetaonezRec.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0SinThetaonezRec->Sumw2();
    d0SinThetaonezRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0SinThetaonezRec);
    
    named0sinThetatworphiRec = "d0sinthetatworphiRec_";
    named0sinThetatworphiRec += i;
    named0sinThetatwozRec ="d0sinthetatwozRec_";
    named0sinThetatwozRec += i;
    d0SinThetatworphiRec = new TH1F(named0sinThetatworphiRec.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0SinThetatworphiRec->Sumw2();
    d0SinThetatworphiRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0SinThetatworphiRec);
    d0SinThetatwozRec = new TH1F(named0sinThetatwozRec.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0SinThetatwozRec->Sumw2();
    d0SinThetatwozRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0SinThetatwozRec);
    
    named0sinThetathreerphiRec = "d0sinthetathreerphiRec_";
    named0sinThetathreerphiRec += i;
    named0sinThetathreezRec ="d0sinthetathreezRec_";
    named0sinThetathreezRec += i;
    
    d0SinThetathreerphiRec = new TH1F(named0sinThetathreerphiRec.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0SinThetathreerphiRec->Sumw2();
    d0SinThetathreerphiRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0SinThetathreerphiRec);
    d0SinThetathreezRec = new TH1F(named0sinThetathreezRec.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0SinThetathreezRec->Sumw2();
    d0SinThetathreezRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0SinThetathreezRec);
    
    named0sinThetafourrphiRec = "d0sinthetafourrphiRec_";
    named0sinThetafourrphiRec += i;
    named0sinThetafourzRec ="d0sinthetafourzRec_";
    named0sinThetafourzRec += i;
    d0SinThetafourrphiRec = new TH1F(named0sinThetafourrphiRec.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0SinThetafourrphiRec->Sumw2();
    d0SinThetafourrphiRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0SinThetafourrphiRec);
    d0SinThetafourzRec = new TH1F(named0sinThetafourzRec.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0SinThetafourzRec->Sumw2();
    d0SinThetafourzRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0SinThetafourzRec);

    named0thetaForwardrphiRec = "d0thetaforwardrphiRec_";
    named0thetaForwardrphiRec += i;
    named0thetaForwardzRec ="d0thetaforwardzRec_";
    named0thetaForwardzRec += i;
    d0ThetaforwardrphiRec = new TH1F(named0thetaForwardrphiRec.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0ThetaforwardrphiRec->Sumw2();
    d0ThetaforwardrphiRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0ThetaforwardrphiRec);
    d0ThetaforwardzRec = new TH1F(named0thetaForwardzRec.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0ThetaforwardzRec->Sumw2();
    d0ThetaforwardzRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0ThetaforwardzRec);

    named0thetaBackwardrphiRec = "d0thetabackwardrphiRec_";
    named0thetaBackwardrphiRec += i;
    named0thetaBackwardzRec ="d0thetabackwardzRec_";
    named0thetaBackwardzRec += i;
    d0ThetabackwardrphiRec = new TH1F(named0thetaBackwardrphiRec.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0ThetabackwardrphiRec->Sumw2();
    d0ThetabackwardrphiRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0ThetabackwardrphiRec);
    d0ThetabackwardzRec = new TH1F(named0thetaBackwardzRec.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0ThetabackwardzRec->Sumw2();
    d0ThetabackwardzRec->SetMinimum(0);
    fOutputSinThetaRec->Add(d0ThetabackwardzRec);
    
    named0sinThetaonerphiSkip = "d0sinthetaonerphiSkip_";
    named0sinThetaonerphiSkip += i;
    named0sinThetaonezSkip ="d0sinthetaonezSkip_";
    named0sinThetaonezSkip += i;
    d0SinThetaonerphiSkip = new TH1F(named0sinThetaonerphiSkip.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0SinThetaonerphiSkip->Sumw2();
    d0SinThetaonerphiSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0SinThetaonerphiSkip);
    d0SinThetaonezSkip = new TH1F(named0sinThetaonezSkip.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0SinThetaonezSkip->Sumw2();
    d0SinThetaonezSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0SinThetaonezSkip);
    
    named0sinThetatworphiSkip = "d0sinthetatworphiSkip_";
    named0sinThetatworphiSkip += i;
    named0sinThetatwozSkip ="d0sinthetatwozSkip_";
    named0sinThetatwozSkip += i;
    d0SinThetatworphiSkip = new TH1F(named0sinThetatworphiSkip.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0SinThetatworphiSkip->Sumw2();
    d0SinThetatworphiSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0SinThetatworphiSkip);
    d0SinThetatwozSkip = new TH1F(named0sinThetatwozSkip.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0SinThetatwozSkip->Sumw2();
    d0SinThetatwozSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0SinThetatwozSkip);
    
    named0sinThetathreerphiSkip = "d0sinthetathreerphiSkip_";
    named0sinThetathreerphiSkip += i;
    named0sinThetathreezSkip ="d0sinthetathreezSkip_";
    named0sinThetathreezSkip += i;
    
    d0SinThetathreerphiSkip = new TH1F(named0sinThetathreerphiSkip.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0SinThetathreerphiSkip->Sumw2();
    d0SinThetathreerphiSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0SinThetathreerphiSkip);
    d0SinThetathreezSkip = new TH1F(named0sinThetathreezSkip.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0SinThetathreezSkip->Sumw2();
    d0SinThetathreezSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0SinThetathreezSkip);
    
    named0sinThetafourrphiSkip = "d0sinthetafourrphiSkip_";
    named0sinThetafourrphiSkip += i;
    named0sinThetafourzSkip ="d0sinthetafourzSkip_";
    named0sinThetafourzSkip += i;
    d0SinThetafourrphiSkip = new TH1F(named0sinThetafourrphiSkip.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0SinThetafourrphiSkip->Sumw2();
    d0SinThetafourrphiSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0SinThetafourrphiSkip);
    d0SinThetafourzSkip = new TH1F(named0sinThetafourzSkip.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0SinThetafourzSkip->Sumw2();
    d0SinThetafourzSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0SinThetafourzSkip);

    named0thetaForwardrphiSkip = "d0thetaforwardrphiSkip_";
    named0thetaForwardrphiSkip += i;
    named0thetaForwardzSkip ="d0thetaforwardzSkip_";
    named0thetaForwardzSkip += i;
    d0ThetaforwardrphiSkip = new TH1F(named0thetaForwardrphiSkip.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0ThetaforwardrphiSkip->Sumw2();
    d0ThetaforwardrphiSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0ThetaforwardrphiSkip);
    d0ThetaforwardzSkip = new TH1F(named0thetaForwardzSkip.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0ThetaforwardzSkip->Sumw2();
    d0ThetaforwardzSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0ThetaforwardzSkip);

    named0thetaBackwardrphiSkip = "d0thetabackwardrphiSkip_";
    named0thetaBackwardrphiSkip += i;
    named0thetaBackwardzSkip ="d0thetabackwardzSkip_";
    named0thetaBackwardzSkip += i;
    d0ThetabackwardrphiSkip = new TH1F(named0thetaBackwardrphiSkip.Data(),d0sinThetarphiTitle.Data(),400,-2000,2000);
    d0ThetabackwardrphiSkip->Sumw2();
    d0ThetabackwardrphiSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0ThetabackwardrphiSkip);
    d0ThetabackwardzSkip = new TH1F(named0thetaBackwardzSkip.Data(),d0sinThetazTitle.Data(),400,-2000,2000);
    d0ThetabackwardzSkip->Sumw2();
    d0ThetabackwardzSkip->SetMinimum(0);
    fOutputSinThetaSkip->Add(d0ThetabackwardzSkip);

  }

  const Int_t nhistphi=20;
  for(Int_t i=0; i<=nhistphi; i++) {

    named0phiAllpointrphiSkip = "d0phiallpointrphiSkip_";
    named0phiAllpointrphiSkip += i;
    named0phiAllpointzSkip ="d0phiallpointzSkip_";
    named0phiAllpointzSkip += i;
    d0PhiAllpointrphiSkip = new TH1F(named0phiAllpointrphiSkip.Data(),d0phiAllpointrphiTitle.Data(),400,-2000,2000);
    d0PhiAllpointrphiSkip->Sumw2();
    d0PhiAllpointrphiSkip->SetMinimum(0);
    fOutputphiAllpointSkip->Add(d0PhiAllpointrphiSkip);
    d0PhiAllpointzSkip = new TH1F(named0phiAllpointzSkip.Data(),d0phiAllpointzTitle.Data(),400,-2000,2000);
    d0PhiAllpointzSkip->Sumw2();
    d0PhiAllpointzSkip->SetMinimum(0);
    fOutputphiAllpointSkip->Add(d0PhiAllpointzSkip);


    named0phiPostvtracrphiSkip = "d0phipostvtracrphiSkip_";
    named0phiPostvtracrphiSkip += i;
    named0phiPostvtraczSkip ="d0phipostvtraczSkip_";
    named0phiPostvtraczSkip += i;
    d0PhiPostvtracrphiSkip = new TH1F(named0phiPostvtracrphiSkip.Data(),d0phiPostvtracrphiTitle.Data(),400,-2000,2000);
    d0PhiPostvtracrphiSkip->Sumw2();
    d0PhiPostvtracrphiSkip->SetMinimum(0);
    fOutputphiPostvtracSkip->Add(d0PhiPostvtracrphiSkip);
    d0PhiPostvtraczSkip = new TH1F(named0phiPostvtraczSkip.Data(),d0phiPostvtraczTitle.Data(),400,-2000,2000);
    d0PhiPostvtraczSkip->Sumw2();
    d0PhiPostvtraczSkip->SetMinimum(0);
    fOutputphiPostvtracSkip->Add(d0PhiPostvtraczSkip);


    named0phiNegtvtracrphiSkip = "d0phinegtvtracrphiSkip_";
    named0phiNegtvtracrphiSkip += i;
    named0phiNegtvtraczSkip ="d0phinegtvtraczSkip_";
    named0phiNegtvtraczSkip += i;
    d0PhiNegtvtracrphiSkip = new TH1F(named0phiNegtvtracrphiSkip.Data(),d0phiNegtvtracrphiTitle.Data(),400,-2000,2000);
    d0PhiNegtvtracrphiSkip->Sumw2();
    d0PhiNegtvtracrphiSkip->SetMinimum(0);
    fOutputphiNegtvtracSkip->Add(d0PhiNegtvtracrphiSkip);
    d0PhiNegtvtraczSkip = new TH1F(named0phiNegtvtraczSkip.Data(),d0phiNegtvtraczTitle.Data(),400,-2000,2000);
    d0PhiNegtvtraczSkip->Sumw2();
    d0PhiNegtvtraczSkip->SetMinimum(0);
    fOutputphiNegtvtracSkip->Add(d0PhiNegtvtraczSkip);
  }

  if(!fNentries) fNentries = new TH1F("hNentries", "number of entries", 26, 0., 40.);
  if(!fEstimVtx) fEstimVtx = new TH1F("vtxRes","Resolution of vertex",1000,-5000.,5000);
  PostData(1, fOutputitspureSARec);
  PostData(2, fOutputitspureSASkip);
  PostData(3, fOutputallPointRec);
  PostData(4, fOutputallPointSkip);
  PostData(5, fOutputpartPointRec);
  PostData(6, fOutputpartPointSkip);
  PostData(7, fOutputonepointSPDRec);
  PostData(8, fOutputonepointSPDSkip);
  PostData(9, fOutputpostvTracRec);
  PostData(10, fOutputpostvTracSkip);
  PostData(11, fOutputnegtvTracRec);
  PostData(12, fOutputnegtvTracSkip);
  PostData(13, fOutputpullAllpointRec);
  PostData(14, fOutputpullAllpointSkip);
  PostData(15, fOutputOnlyRefitRec);
  PostData(16, fOutputOnlyRefitSkip);
  PostData(17, fOutputSinThetaRec);
  PostData(18, fOutputSinThetaSkip);
  PostData(19, fOutputallPointTrue);
  PostData(20, fOutputpostvTracTrue);
  PostData(21, fOutputnegtvTracTrue);
  PostData(22, fOutputpullAllpointTrue);
  PostData(23, fOutputphiAllpointSkip);
  PostData(24, fOutputphiPostvtracSkip);
  PostData(25, fOutputphiNegtvtracSkip);
  PostData(26, fOutputparticlePID);
  PostData(27, fOutputPt);
  PostData(28, fNentries);
  PostData(29, fEstimVtx);  
  PostData(30, fOutputWithTrackCuts);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEImpParRes::UserExec(Option_t */*option*/)
{
  //
  // Track selection and filling of d0 histograms
  //
  AliVEvent* event = dynamic_cast<AliVEvent*>(InputEvent());
  if (!event) {
    AliError("event not found. Nothing done!");
    return;
  }

  // only events in the requested multiplicity range
  TString firedTriggerClasses="";
  Int_t runNumber=0;
  if(fIsAOD){
    Int_t nclsITS = 0;
    runNumber=((AliAODEvent*)event)->GetRunNumber();
    nclsITS = ((AliAODEvent*)event)->GetHeader()->GetNumberOfITSClusters(1);
    if(nclsITS<fMinMult || nclsITS>fMaxMult) return;
    firedTriggerClasses=((AliAODEvent*)event)->GetFiredTriggerClasses();
    if(!firedTriggerClasses.Contains(fTriggerClass.Data())) return;
  }
  else{
    runNumber=((AliESDEvent*)event)->GetRunNumber();
    if(!IsSelectedCentrality(((AliESDEvent*)event))) return;
    firedTriggerClasses=((AliESDEvent*)event)->GetFiredTriggerClasses();
    if(!firedTriggerClasses.Contains(fTriggerClass.Data())) return;
  }

 

  Bool_t sddIsIn=kTRUE;
  if(fCheckSDDIsIn) {

    if(!fTrigConfig) {    
      AliCDBManager* man = AliCDBManager::Instance();
      if(fOCDBPath.Contains("OCDB")) { // when running in the QAtrain this is not called (OCBD is already set)
        man->SetDefaultStorage(fOCDBPath.Data());
        man->SetRun(runNumber);
      }
      if(!man) {      
        AliFatal("CDB not set but needed by AliAnalysisTaskITSTrackingCheck");
        return;    
      }  
      AliCDBEntry* eT=(AliCDBEntry*)man->Get("GRP/CTP/Config");    
      if(eT) {      
        fTrigConfig=(AliTriggerConfiguration*)eT->GetObject();    
      }    
      if(!eT || !fTrigConfig) {      
        AliError("Cannot retrieve CDB entry for GRP/CTP/Config");      
        return;    
      }
    }

    if(fIsAOD){
      const TObjArray& classesArray=fTrigConfig->GetClasses();
      ULong64_t trigMask=((AliAODEvent*)event)->GetTriggerMask();
      Int_t nclasses = classesArray.GetEntriesFast();
      for(Int_t iclass=0; iclass < nclasses; iclass++ ) 
	{
	  AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
	  ULong64_t classMask=trclass->GetMask();
	  if(trigMask & classMask)
	    {
	      TString detList=trclass->GetCluster()->GetDetectorsInCluster();
	      if(detList.Contains("ITSSDD")) sddIsIn = kTRUE;
	      else sddIsIn = kFALSE;
	    }
	}
      //sddIsIn = kFALSE;
    }
    else sddIsIn=((AliESDEvent*)event)->IsDetectorInTriggerCluster("ITSSDD",fTrigConfig);
    if(fCheckSDDIsIn==1 && !sddIsIn) return;
    if(fCheckSDDIsIn==-1 && sddIsIn) return;
  }

  fNentries->Fill(1);


  Int_t nTrks = event->GetNumberOfTracks();
  Bool_t highMult=(nTrks>500 ? kTRUE : kFALSE);

  Double_t vtxTrue[3];
  AliStack *stack=0;
  TClonesArray *mcArray=0;
  AliESDVertex *vtxESDTrue=0;
  AliVVertex *vtxVSkip=0;
  AliVVertex *vtxVRec=0;  
  AliVVertex* primaryVtx=0;


  // event primary vertex
  AliVertexerTracks vertexer0(event->GetMagneticField());
  vertexer0.SetITSMode();
  vertexer0.SetMinClusters(3);
  if(highMult) vertexer0.SetITSMode(0.1,0.1,0.5,5,1,3.,100.,1000.,3.,30.,1,1);
  if(fUseDiamond){
    // diamond constraint
    Float_t diamondcovxy[3];
    event->GetDiamondCovXY(diamondcovxy);
    Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
    Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.};
    AliESDVertex diamond(pos,cov,1.,1);
    vertexer0.SetVtxStart(&diamond);
  }
  vtxVRec=(AliVVertex*)vertexer0.FindPrimaryVertex(event);
  if(!vtxVRec) return;
  if(vtxVRec->GetNContributors()<1){
    delete vtxVRec; vtxVRec=NULL;
    return;
  }
  
  if (fReadMC) {
    if (fIsAOD){
      mcArray = dynamic_cast<TClonesArray*>(((AliAODEvent*)event)->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!mcArray){
	AliError("Clould not find Monte-Carlo in AOD");
	return;
      }
      AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(((AliAODEvent*)(event))->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
      if (!mcHeader) {
	AliError("Could not find MC Header in AOD");
	return;
      }
      
      Double_t mcVertex[3]={9999.,9999.,9999.};
      mcHeader->GetVertex(mcVertex);
      vtxTrue[0]=mcVertex[0];vtxTrue[1]=mcVertex[1];vtxTrue[2]=mcVertex[2];
      Double_t sigmaTrue[3]={0., 0., 0.,};
      vtxESDTrue = new AliESDVertex(vtxTrue,sigmaTrue);
    }//end if isAOD
    else{
      AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!eventHandler) {
	Printf("ERROR: Could not retrieve MC event handler");
	return;
      }
      
      AliMCEvent* mcEvent = eventHandler->MCEvent();
      if (!mcEvent) {
	Printf("ERROR: Could not retrieve MC event");
	return;
      }
      
      stack = mcEvent->Stack();
      if (!stack) {
	AliDebug(AliLog::kError, "Stack not available");
	return;
      }
      
      //load MC header for ESD;//see $ALICE_ROOT/PWGPP/global/AliAnalysisTaskSEVertexESD.cxx
      AliHeader *mcHeader = eventHandler->MCEvent()->Header();
      if (!mcHeader) {
	AliDebug(AliLog::kError, "Header not available");
	return;
      }

      AliGenEventHeader* genHeader = mcHeader->GenEventHeader();  
      TArrayF mcVertex(3);
      mcVertex[0]=9999.; mcVertex[1]=9999.; mcVertex[2]=9999.;
      genHeader->PrimaryVertex(mcVertex);
      vtxTrue[0]=mcVertex[0];vtxTrue[1]=mcVertex[1];vtxTrue[2]=mcVertex[2];
      Double_t sigmaTrue[3]={0., 0., 0.,};
      //mcHeader->GetVertex(vtxTrue);//note the vtxTrue is void here,so must need the next line.
      //AliESDVertex *vtxESDTrue = new AliESDVertex(vtxTrue,sigmaTrue);
      vtxESDTrue = new AliESDVertex(vtxTrue,sigmaTrue);
      
    }//end else (!isAOD)
  }
  
  Double_t beampiperadius=3.;
  AliVTrack *vtrack = 0;
  Int_t pdgCode=0;
  Int_t trkLabel;
  TParticle  *part =0;
  AliAODMCParticle *AODpart=0;
  Int_t npointsITS=0,npointsSPD=0;
  Int_t skipped[2];
  Double_t dzRec[2], covdzRec[3], dzRecSkip[2], covdzRecSkip[3],dzTrue[2], covdzTrue[3];
  Double_t pt;
  Int_t bin;
  Int_t nClsTotTPC=0;
  Bool_t haskITSrefit=kFALSE;
  Bool_t haskTPCrefit=kFALSE;
  Int_t charge=0;
  Double_t phi=0.;
  Double_t theta=0.;
  Double_t eta=0.;

  
  for (Int_t it=0; it<nTrks; it++){ //start loop over tracks
    vtrack = (AliVTrack*)event->GetTrack(it);
    if(!vtrack) continue;

    eta = vtrack->Eta();
    if(eta<-0.8 || eta>0.8) continue;
    
    npointsITS=0; npointsSPD=0;
    if(fIsAOD){
      haskITSrefit=(((AliAODTrack*)vtrack)->GetStatus()&AliESDtrack::kITSrefit);
      haskTPCrefit=(((AliAODTrack*)vtrack)->GetStatus()&AliESDtrack::kTPCrefit);
      nClsTotTPC=((AliAODTrack*)vtrack)->GetTPCNcls();
      if(!haskITSrefit) continue;
      for(Int_t ilayer=0; ilayer<6; ilayer++){
	if (ilayer<2 && ((AliAODTrack*)vtrack)->HasPointOnITSLayer(ilayer)) npointsSPD++;
	if (((AliAODTrack*)vtrack)->HasPointOnITSLayer(ilayer)) npointsITS++;  
      }
    }
    else {
      haskITSrefit=(((AliESDtrack*)vtrack)->GetStatus()&AliESDtrack::kITSrefit);
      haskTPCrefit=(((AliESDtrack*)vtrack)->GetStatus()&AliESDtrack::kTPCrefit);
      nClsTotTPC=((AliESDtrack*)vtrack)->GetTPCNcls();
      if(!haskITSrefit) continue;
      for (Int_t ilayer=0; ilayer<6; ilayer++){ 
	if (ilayer<2 && ((AliESDtrack*)vtrack)->HasPointOnITSLayer(ilayer)) npointsSPD++;
	if (((AliESDtrack*)vtrack)->HasPointOnITSLayer(ilayer)) npointsITS++;  
      }
    }
    charge=vtrack->Charge();
    phi=vtrack->Phi();
    theta=vtrack->Theta();

    //MC
    if (fReadMC){
      trkLabel = vtrack->GetLabel();
      if(trkLabel<0) continue;
      if(fIsAOD && mcArray){
	AODpart = (AliAODMCParticle*)mcArray->At(trkLabel);
	if(!AODpart) printf("NOPART\n");
	pdgCode = TMath::Abs(AODpart->GetPdgCode());	
      }
      if(!fIsAOD && stack) {
	part = (TParticle*)stack->Particle(trkLabel);
	pdgCode = TMath::Abs(part->GetPdgCode());
      }
      //pdgCode = TMath::Abs(part->GetPdgCode());
      //printf("pdgCode===%d\n", pdgCode);
      if(fSelectedPdg>0 && pdgCode!=fSelectedPdg) continue;
    }
    
      
      //Get specific primary vertex--Reconstructed primary vertex do not include the track considering.
    AliVertexerTracks vertexer(event->GetMagneticField());
    vertexer.SetITSMode();
    vertexer.SetMinClusters(3);
    if(fUseDiamond){
      Float_t diamondcovxy[3];
      event->GetDiamondCovXY(diamondcovxy);
      Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
      Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.};
      AliESDVertex diamond(pos,cov,1.,1);
      vertexer.SetVtxStart(&diamond);
    }
    skipped[0] = (Int_t)vtrack->GetID();
    vertexer.SetSkipTracks(1,skipped);      
    // create vertex with new!
    if(!highMult && fSkipTrack) {
      vtxVSkip = (AliVVertex*)vertexer.FindPrimaryVertex(event);
      if(!vtxVSkip) continue;
      if(vtxVSkip->GetNContributors()<1) {
	delete vtxVSkip; vtxVSkip=NULL;
	continue;
      }
    } // else {
      // vtxVSkip = new AliVVertex(); produce error!!!
      // } 
 
    pt = vtrack->Pt();
    bin = PtBin(pt);
 
    if(bin==-1) {
      delete vtxVSkip; vtxVSkip=NULL;
      continue;
    }

    // Select primary particle if MC event (for ESD event), Rprod < 1 micron
    if(fReadMC){
      if(fIsAOD){
	if((AODpart->Xv()-vtxTrue[0])*(AODpart->Xv()-vtxTrue[0])+
	   (AODpart->Yv()-vtxTrue[1])*(AODpart->Yv()-vtxTrue[1])
	   > 0.0001*0.0001) {
	  delete vtxVSkip; vtxVSkip=NULL;
	  continue;
	}
      }
      else{
	if((part->Vx()-vtxTrue[0])*(part->Vx()-vtxTrue[0])+
	   (part->Vy()-vtxTrue[1])*(part->Vy()-vtxTrue[1])
	   > 0.0001*0.0001) {
	  delete vtxVSkip; vtxVSkip=NULL;
	  continue;
	}
      }
    }
    
    
    // compute impact patameters
    // wrt event vertex
    vtrack->PropagateToDCA(vtxVRec, event->GetMagneticField(), beampiperadius, dzRec, covdzRec);
    // wrt event vertex without this track
    if(!highMult && fSkipTrack) {
      vtrack->PropagateToDCA(vtxVSkip, event->GetMagneticField(), beampiperadius, dzRecSkip, covdzRecSkip);
    } else if(!fSkipTrack) {
      dzRecSkip[0]=dzRec[0]; 
      dzRecSkip[1]=dzRec[1];
      covdzRecSkip[0]=covdzRec[0];
      covdzRecSkip[1]=covdzRec[1];
      covdzRecSkip[2]=covdzRec[2];
    } else {
      dzRecSkip[0]=0; 
      dzRecSkip[1]=0;
      covdzRecSkip[0]=0;
      covdzRecSkip[1]=0;
      covdzRecSkip[2]=0;
    }
    //delete vtxVSkip; vtxVSkip=NULL; // not needed anymore
 
    if(fReadMC) vtrack->PropagateToDCA(vtxESDTrue, event->GetMagneticField(), beampiperadius, dzTrue, covdzTrue);
    if(covdzRec[0]<1.e-13 || covdzRec[2]<1.e-13 || covdzRecSkip[0]<1.e-13 || covdzRecSkip[2]<1.e-13) continue;
    if(fReadMC && (covdzTrue[0]<1.e-13 || covdzTrue[2]<1.e-13)) continue;


    // Bayesian PID only for ESD
    if(!fIsAOD && (npointsITS==6 || (npointsITS==4 && !sddIsIn))){
      Double_t prob[AliPID::kSPECIES];
      ((AliESDtrack*)(vtrack))->GetESDpid(prob);
      Double_t priors[5] = {0.01, 0.01, 0.85, 0.10, 0.05};
      
      
      AliPID pid;
      pid.SetPriors(priors);
      pid.SetProbabilities(prob);
      
      // identify particle as the most probable    
      Double_t pelectron = pid.GetProbability(AliPID::kElectron);
      Double_t pmuon = pid.GetProbability(AliPID::kMuon);
      Double_t ppion = pid.GetProbability(AliPID::kPion);
      Double_t pkaon = pid.GetProbability(AliPID::kKaon);
      Double_t pproton = pid.GetProbability(AliPID::kProton);  
      
      if (ppion > pelectron &&
	  ppion > pmuon  &&
	  ppion > pkaon &&
	  ppion > pproton ) {
	//esdPid =-kPDGelectron;
	char *named0PionPIDrphiRec = Form("d0pionPIDrphiRec_%d", bin);
	char *named0PionPIDzRec = Form("d0pionPIDzRec_%d", bin);
	char *named0PionPIDrphiSkip = Form("d0pionPIDrphiSkip_%d", bin);
	char *named0PionPIDzSkip = Form("d0pionPIDzSkip_%d", bin);
	((TH1F*)(fOutputparticlePID->FindObject(named0PionPIDrphiRec)))->Fill(10000.*dzRec[0]);
	((TH1F*)(fOutputparticlePID->FindObject(named0PionPIDzRec)))->Fill(10000.*dzRec[1]);
	((TH1F*)(fOutputparticlePID->FindObject(named0PionPIDrphiSkip)))->Fill(10000.*dzRecSkip[0]);
	((TH1F*)(fOutputparticlePID->FindObject(named0PionPIDzSkip)))->Fill(10000.*dzRecSkip[1]);
      }
      
      
      if (pkaon > pelectron &&
	  pkaon > pmuon  &&
	  pkaon > ppion &&
	  pkaon > pproton ) {
	//esdPid =-kPDGelectron;
	char *named0KaonPIDrphiRec = Form("d0kaonPIDrphiRec_%d", bin);
	char *named0KaonPIDzRec = Form("d0kaonPIDzRec_%d", bin);
	char *named0KaonPIDrphiSkip = Form("d0kaonPIDrphiSkip_%d", bin);
	char *named0KaonPIDzSkip = Form("d0kaonPIDzSkip_%d", bin);
	((TH1F*)(fOutputparticlePID->FindObject(named0KaonPIDrphiRec)))->Fill(10000.*dzRec[0]);
	((TH1F*)(fOutputparticlePID->FindObject(named0KaonPIDzRec)))->Fill(10000.*dzRec[1]);
	((TH1F*)(fOutputparticlePID->FindObject(named0KaonPIDrphiSkip)))->Fill(10000.*dzRecSkip[0]);
	((TH1F*)(fOutputparticlePID->FindObject(named0KaonPIDzSkip)))->Fill(10000.*dzRecSkip[1]);
      }
      
      
      if (pproton > pelectron &&
	  pproton >pmuon  &&
	  pproton > ppion &&
	  pproton > pkaon ) {
	//esdPid =-kPDGelectron;
	//if(p<0.5 && fReadMC){fEstimVtx->Fill(pdgCode);}
	char *named0ProtonPIDrphiRec = Form("d0protonPIDrphiRec_%d", bin);
	char *named0ProtonPIDzRec = Form("d0protonPIDzRec_%d", bin);
	char *named0ProtonPIDrphiSkip = Form("d0protonPIDrphiSkip_%d", bin);
	char *named0ProtonPIDzSkip = Form("d0protonPIDzSkip_%d", bin);
	((TH1F*)(fOutputparticlePID->FindObject(named0ProtonPIDrphiRec)))->Fill(10000.*dzRec[0]);
	((TH1F*)(fOutputparticlePID->FindObject(named0ProtonPIDzRec)))->Fill(10000.*dzRec[1]);
	((TH1F*)(fOutputparticlePID->FindObject(named0ProtonPIDrphiSkip)))->Fill(10000.*dzRecSkip[0]);
	((TH1F*)(fOutputparticlePID->FindObject(named0ProtonPIDzSkip)))->Fill(10000.*dzRecSkip[1]);
      }
    }

    // ESD TRACK CUTS
    if(fReadMC) primaryVtx=vtxESDTrue;
    else if(fSkipTrack) primaryVtx=vtxVSkip;
    else primaryVtx=vtxVRec;

    if(IsTrackSelected(vtrack,primaryVtx,fESDtrackCuts)){


      char *named0PtESDTC = Form("d0ptESDTC_%d",bin);
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0PtESDTC)))->Fill(pt);

       	char *named0DistrESDTCrphiRec = Form("d0DistrESDTCrphiRec_%d", bin);
	char *named0DistrESDTCrphiSkip = Form("d0DistrESDTCrphiSkip_%d", bin);
	char *named0DistrESDTCrphiTrue = Form("d0DistrESDTCrphiTrue_%d", bin);
	char *named0DistrESDTCzRec = Form("d0DistrESDTCzRec_%d", bin);
	char *named0DistrESDTCzSkip = Form("d0DistrESDTCzSkip_%d", bin);
	char *named0DistrESDTCzTrue = Form("d0DistrESDTCzTrue_%d", bin);
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0DistrESDTCrphiRec)))->Fill(10000.*dzRec[0]);
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0DistrESDTCzRec)))->Fill(10000.*dzRec[1]);
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0DistrESDTCrphiSkip)))->Fill(10000.*dzRecSkip[0]);
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0DistrESDTCzSkip)))->Fill(10000.*dzRecSkip[1]);
	
 
	if(fReadMC) {
	  ((TH1F*)(fOutputWithTrackCuts->FindObject(named0DistrESDTCrphiTrue)))->Fill(10000.*dzTrue[0]);
	  ((TH1F*)(fOutputWithTrackCuts->FindObject(named0DistrESDTCzTrue)))->Fill(10000.*dzTrue[1]);
	}
	
	// pulls
	char *named0PullESDTCrphiRec = Form("d0PullESDTCrphiRec_%d", bin);
	char *named0PullESDTCrphiSkip = Form("d0PullESDTCrphiSkip_%d", bin);
	char *named0PullESDTCrphiTrue = Form("d0PullESDTCrphiTrue_%d", bin);
	char *named0PullESDTCzRec = Form("d0PullESDTCzRec_%d", bin);
	char *named0PullESDTCzSkip = Form("d0PullESDTCzSkip_%d", bin);
	char *named0PullESDTCzTrue = Form("d0PullESDTCzTrue_%d", bin);
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0PullESDTCrphiRec)))->Fill(dzRec[0]/TMath::Sqrt(covdzRec[0]));    
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0PullESDTCzRec)))->Fill(dzRec[1]/TMath::Sqrt(covdzRec[2]));
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0PullESDTCrphiSkip)))->Fill(dzRecSkip[0]/TMath::Sqrt(covdzRecSkip[0]));
	((TH1F*)(fOutputWithTrackCuts->FindObject(named0PullESDTCzSkip)))->Fill(dzRecSkip[1]/TMath::Sqrt(covdzRecSkip[2]));
	if(fReadMC) {
	  ((TH1F*)(fOutputWithTrackCuts->FindObject(named0PullESDTCrphiTrue)))->Fill(dzTrue[0]/TMath::Sqrt(covdzTrue[0]));
	  ((TH1F*)(fOutputWithTrackCuts->FindObject(named0PullESDTCzTrue)))->Fill(dzTrue[1]/TMath::Sqrt(covdzTrue[2]));
	}
	
	}

    

    // ITS standalone
    if (nClsTotTPC==0 && haskITSrefit && npointsSPD>0 && npointsITS>=4) {
      char *named0ITSpureSArphiRec = Form("d0itspureSArphiRec_%d", bin);
      char *named0ITSpureSArphiSkip = Form("d0itspureSArphiSkip_%d", bin);
      char *named0ITSpureSAzRec = Form("d0itspureSAzRec_%d", bin);
      char *named0ITSpureSAzSkip = Form("d0itspureSAzSkip_%d", bin);
      ((TH1F*)(fOutputitspureSARec->FindObject(named0ITSpureSArphiRec)))->Fill(10000.*dzRec[0]);
      ((TH1F*)(fOutputitspureSARec->FindObject(named0ITSpureSAzRec)))->Fill(10000.*dzRec[1]);
      ((TH1F*)(fOutputitspureSASkip->FindObject(named0ITSpureSArphiSkip)))->Fill(10000.*dzRecSkip[0]);
      ((TH1F*)(fOutputitspureSASkip->FindObject(named0ITSpureSAzSkip)))->Fill(10000.*dzRecSkip[1]);
    }


   
    // ask for TPC refit
    if (!haskTPCrefit || nClsTotTPC<70) continue;
   
    // only ITS and TPC refit
    char *named0OnlyrefitrphiRec = Form("d0onlyrefitrphiRec_%d", bin);
    char *named0OnlyrefitrphiSkip = Form("d0onlyrefitrphiSkip_%d", bin);
    char *named0OnlyrefitzRec = Form("d0onlyrefitzRec_%d", bin);
    char *named0OnlyrefitzSkip = Form("d0onlyrefitzSkip_%d", bin);
    ((TH1F*)(fOutputOnlyRefitRec->FindObject(named0OnlyrefitrphiRec)))->Fill(10000.*dzRec[0]);
    ((TH1F*)(fOutputOnlyRefitRec->FindObject(named0OnlyrefitzRec)))->Fill(10000.*dzRec[1]);
    ((TH1F*)(fOutputOnlyRefitSkip->FindObject(named0OnlyrefitrphiSkip)))->Fill(10000.*dzRecSkip[0]);
    ((TH1F*)(fOutputOnlyRefitSkip->FindObject(named0OnlyrefitzSkip)))->Fill(10000.*dzRecSkip[1]);  
   

    if(npointsITS>=4 && npointsSPD>0) {
      char *named0PartpointrphiRec = Form("d0partpointrphiRec_%d", bin);
      char *named0PartpointrphiSkip = Form("d0partpointrphiSkip_%d", bin);
      char *named0PartpointzRec = Form("d0partpointzRec_%d", bin);
      char *named0PartpointzSkip = Form("d0partpointzSkip_%d", bin);
      ((TH1F*)(fOutputpartPointRec->FindObject(named0PartpointrphiRec)))->Fill(10000.*dzRec[0]);
      ((TH1F*)(fOutputpartPointRec->FindObject(named0PartpointzRec)))->Fill(10000.*dzRec[1]);
      ((TH1F*)(fOutputpartPointSkip->FindObject(named0PartpointrphiSkip)))->Fill(10000.*dzRecSkip[0]);
      ((TH1F*)(fOutputpartPointSkip->FindObject(named0PartpointzSkip)))->Fill(10000.*dzRecSkip[1]);
    }

    if(npointsSPD>0) {
      char *named0OnepointSPDrphiRec = Form("d0onepointSPDrphiRec_%d", bin);
      char *named0OnepointSPDrphiSkip = Form("d0onepointSPDrphiSkip_%d", bin);
      char *named0OnepointSPDzRec = Form("d0onepointSPDzRec_%d", bin);
      char *named0OnepointSPDzSkip = Form("d0onepointSPDzSkip_%d", bin);
      ((TH1F*)(fOutputonepointSPDRec->FindObject(named0OnepointSPDrphiRec)))->Fill(10000.*dzRec[0]);
      ((TH1F*)(fOutputonepointSPDRec->FindObject(named0OnepointSPDzRec)))->Fill(10000.*dzRec[1]);
      ((TH1F*)(fOutputonepointSPDSkip->FindObject(named0OnepointSPDrphiSkip)))->Fill(10000.*dzRecSkip[0]);
      ((TH1F*)(fOutputonepointSPDSkip->FindObject(named0OnepointSPDzSkip)))->Fill(10000.*dzRecSkip[1]);
    }

   
    if(npointsITS==6 || (npointsITS==4 && !sddIsIn)) {
      //pt
      char *named0Pt = Form("d0pt_%d",bin);
      ((TH1F*)(fOutputPt->FindObject(named0Pt)))->Fill(pt);


      // allpoint
      char *named0AllpointrphiRec = Form("d0allpointrphiRec_%d", bin);
      char *named0AllpointrphiSkip = Form("d0allpointrphiSkip_%d", bin);
      char *named0AllpointrphiTrue = Form("d0allpointrphiTrue_%d", bin);
      char *named0AllpointzRec = Form("d0allpointzRec_%d", bin);
      char *named0AllpointzSkip = Form("d0allpointzSkip_%d", bin);
      char *named0AllpointzTrue = Form("d0allpointzTrue_%d", bin);
      ((TH1F*)(fOutputallPointRec->FindObject(named0AllpointrphiRec)))->Fill(10000.*dzRec[0]);
      ((TH1F*)(fOutputallPointRec->FindObject(named0AllpointzRec)))->Fill(10000.*dzRec[1]);
      ((TH1F*)(fOutputallPointSkip->FindObject(named0AllpointrphiSkip)))->Fill(10000.*dzRecSkip[0]);
      ((TH1F*)(fOutputallPointSkip->FindObject(named0AllpointzSkip)))->Fill(10000.*dzRecSkip[1]);
      if(fReadMC) {
        ((TH1F*)(fOutputallPointTrue->FindObject(named0AllpointrphiTrue)))->Fill(10000.*dzTrue[0]);
        ((TH1F*)(fOutputallPointTrue->FindObject(named0AllpointzTrue)))->Fill(10000.*dzTrue[1]);
      }
     
      // pulls
      char *named0PullAllpointrphiRec = Form("d0pullAllpointrphiRec_%d", bin);
      char *named0PullAllpointrphiSkip = Form("d0pullAllpointrphiSkip_%d", bin);
      char *named0PullAllpointrphiTrue = Form("d0pullAllpointrphiTrue_%d", bin);
      char *named0PullAllpointzRec = Form("d0pullAllpointzRec_%d", bin);
      char *named0PullAllpointzSkip = Form("d0pullAllpointzSkip_%d", bin);
      char *named0PullAllpointzTrue = Form("d0pullAllpointzTrue_%d", bin);
      ((TH1F*)(fOutputpullAllpointRec->FindObject(named0PullAllpointrphiRec)))->Fill(dzRec[0]/TMath::Sqrt(covdzRec[0]));    
      ((TH1F*)(fOutputpullAllpointRec->FindObject(named0PullAllpointzRec)))->Fill(dzRec[1]/TMath::Sqrt(covdzRec[2]));
      ((TH1F*)(fOutputpullAllpointSkip->FindObject(named0PullAllpointrphiSkip)))->Fill(dzRecSkip[0]/TMath::Sqrt(covdzRecSkip[0]));
      ((TH1F*)(fOutputpullAllpointSkip->FindObject(named0PullAllpointzSkip)))->Fill(dzRecSkip[1]/TMath::Sqrt(covdzRecSkip[2]));
      if(fReadMC) {
        ((TH1F*)(fOutputpullAllpointTrue->FindObject(named0PullAllpointrphiTrue)))->Fill(dzTrue[0]/TMath::Sqrt(covdzTrue[0]));
        ((TH1F*)(fOutputpullAllpointTrue->FindObject(named0PullAllpointzTrue)))->Fill(dzTrue[1]/TMath::Sqrt(covdzTrue[2]));
      }
      //postive and negative track
      //Int_t charge=esdtrack->Charge();
      if(charge==1) {
        char *named0PostvtracrphiRec = Form("d0postvtracrphiRec_%d", bin);
        char *named0PostvtracrphiSkip = Form("d0postvtracrphiSkip_%d", bin);
        char *named0PostvtracrphiTrue = Form("d0postvtracrphiTrue_%d", bin);
        char *named0PostvtraczRec = Form("d0postvtraczRec_%d", bin);
        char *named0PostvtraczSkip = Form("d0postvtraczSkip_%d", bin);
        char *named0PostvtraczTrue = Form("d0postvtraczTrue_%d", bin);
        ((TH1F*)(fOutputpostvTracRec->FindObject(named0PostvtracrphiRec)))->Fill(10000.*dzRec[0]);
        ((TH1F*)(fOutputpostvTracRec->FindObject(named0PostvtraczRec)))->Fill(10000.*dzRec[1]);
        ((TH1F*)(fOutputpostvTracSkip->FindObject(named0PostvtracrphiSkip)))->Fill(10000.*dzRecSkip[0]);
        ((TH1F*)(fOutputpostvTracSkip->FindObject(named0PostvtraczSkip)))->Fill(10000.*dzRecSkip[1]);
        if(fReadMC) {
          ((TH1F*)(fOutputpostvTracTrue->FindObject(named0PostvtracrphiTrue)))->Fill(10000.*dzTrue[0]);
          ((TH1F*)(fOutputpostvTracTrue->FindObject(named0PostvtraczTrue)))->Fill(10000.*dzTrue[1]);
        }
      }
     
      if(charge==-1) {
        char *named0NegtvtracrphiRec = Form("d0negtvtracrphiRec_%d", bin);
        char *named0NegtvtracrphiSkip = Form("d0negtvtracrphiSkip_%d", bin);
        char *named0NegtvtracrphiTrue = Form("d0negtvtracrphiTrue_%d", bin);
        char *named0NegtvtraczRec = Form("d0negtvtraczRec_%d", bin);
        char *named0NegtvtraczSkip = Form("d0negtvtraczSkip_%d", bin);
        char *named0NegtvtraczTrue = Form("d0negtvtraczTrue_%d", bin);
        ((TH1F*)(fOutputnegtvTracRec->FindObject(named0NegtvtracrphiRec)))->Fill(10000.*dzRec[0]);
        ((TH1F*)(fOutputnegtvTracRec->FindObject(named0NegtvtraczRec)))->Fill(10000.*dzRec[1]);
        ((TH1F*)(fOutputnegtvTracSkip->FindObject(named0NegtvtracrphiSkip)))->Fill(10000.*dzRecSkip[0]);
        ((TH1F*)(fOutputnegtvTracSkip->FindObject(named0NegtvtraczSkip)))->Fill(10000.*dzRecSkip[1]);
        if(fReadMC) {
          ((TH1F*)(fOutputnegtvTracTrue->FindObject(named0NegtvtracrphiTrue)))->Fill(10000.*dzTrue[0]);
          ((TH1F*)(fOutputnegtvTracTrue->FindObject(named0NegtvtraczTrue)))->Fill(10000.*dzTrue[1]);   
        }    
      }
     
      // SinTheta
      //Double_t theta=esdtrack->Theta();
      Double_t Sintheta=TMath::Sin(theta);
      Double_t pi=TMath::Pi();
      Double_t halfpi=0.5*pi;
      Int_t thetabin = SinThetaBin(Sintheta);
      if(thetabin<0) continue;
      if(bin==4 && theta<halfpi){
        char *named0ThetaforwardrphiRec = Form("d0thetaforwardrphiRec_%d", thetabin);
        char *named0ThetaforwardzRec = Form("d0thetaforwardzRec_%d", thetabin);
        char *named0ThetaforwardrphiSkip = Form("d0thetaforwardrphiSkip_%d", thetabin);
        char *named0ThetaforwardzSkip = Form("d0thetaforwardzSkip_%d", thetabin);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0ThetaforwardrphiRec)))->Fill(10000*dzRec[0]);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0ThetaforwardzRec)))->Fill(10000*dzRec[1]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0ThetaforwardrphiSkip)))->Fill(10000*dzRecSkip[0]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0ThetaforwardzSkip)))->Fill(10000*dzRecSkip[1]);
      }

      if(bin==4 && theta>halfpi){
        char *named0ThetabackwardrphiRec = Form("d0thetabackwardrphiRec_%d", thetabin);
        char *named0ThetabackwardzRec = Form("d0thetabackwardzRec_%d", thetabin);
        char *named0ThetabackwardrphiSkip = Form("d0thetabackwardrphiSkip_%d", thetabin);
        char *named0ThetabackwardzSkip = Form("d0thetabackwardzSkip_%d", thetabin);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0ThetabackwardrphiRec)))->Fill(10000*dzRec[0]);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0ThetabackwardzRec)))->Fill(10000*dzRec[1]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0ThetabackwardrphiSkip)))->Fill(10000*dzRecSkip[0]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0ThetabackwardzSkip)))->Fill(10000*dzRecSkip[1]);
      }
     
      if(bin==1) {
        char *named0SinthetaonerphiRec = Form("d0sinthetaonerphiRec_%d", thetabin);
        char *named0SinthetaonezRec = Form("d0sinthetaonezRec_%d", thetabin);
        char *named0SinthetaonerphiSkip = Form("d0sinthetaonerphiSkip_%d", thetabin);
        char *named0SinthetaonezSkip = Form("d0sinthetaonezSkip_%d", thetabin);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0SinthetaonerphiRec)))->Fill(10000*dzRec[0]);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0SinthetaonezRec)))->Fill(10000*dzRec[1]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0SinthetaonerphiSkip)))->Fill(10000*dzRecSkip[0]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0SinthetaonezSkip)))->Fill(10000*dzRecSkip[1]);
      }
     
      if(bin==5) {
        char *named0SinthetatworphiRec = Form("d0sinthetatworphiRec_%d", thetabin);
        char *named0SinthetatwozRec = Form("d0sinthetatwozRec_%d", thetabin);
        char *named0SinthetatworphiSkip = Form("d0sinthetatworphiSkip_%d", thetabin);
        char *named0SinthetatwozSkip = Form("d0sinthetatwozSkip_%d", thetabin);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0SinthetatworphiRec)))->Fill(10000*dzRec[0]);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0SinthetatwozRec)))->Fill(10000*dzRec[1]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0SinthetatworphiSkip)))->Fill(10000*dzRecSkip[0]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0SinthetatwozSkip)))->Fill(10000*dzRecSkip[1]);
      }
     
      if(bin==10) {
        char *named0SinthetathreerphiRec = Form("d0sinthetathreerphiRec_%d", thetabin);
        char *named0SinthetathreezRec = Form("d0sinthetathreezRec_%d", thetabin);
        char *named0SinthetathreerphiSkip = Form("d0sinthetathreerphiSkip_%d", thetabin);
        char *named0SinthetathreezSkip = Form("d0sinthetathreezSkip_%d", thetabin);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0SinthetathreerphiRec)))->Fill(10000*dzRec[0]);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0SinthetathreezRec)))->Fill(10000*dzRec[1]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0SinthetathreerphiSkip)))->Fill(10000*dzRecSkip[0]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0SinthetathreezSkip)))->Fill(10000*dzRecSkip[1]);
      }
     
      if(bin==15) {
        char *named0SinthetafourrphiRec = Form("d0sinthetafourrphiRec_%d", thetabin);
        char *named0SinthetafourzRec = Form("d0sinthetafourzRec_%d", thetabin);
        char *named0SinthetafourrphiSkip = Form("d0sinthetafourrphiSkip_%d", thetabin);
        char *named0SinthetafourzSkip = Form("d0sinthetafourzSkip_%d", thetabin);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0SinthetafourrphiRec)))->Fill(10000*dzRec[0]);
        ((TH1F*)(fOutputSinThetaRec->FindObject(named0SinthetafourzRec)))->Fill(10000*dzRec[1]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0SinthetafourrphiSkip)))->Fill(10000*dzRecSkip[0]);
        ((TH1F*)(fOutputSinThetaSkip->FindObject(named0SinthetafourzSkip)))->Fill(10000*dzRecSkip[1]);
      }
     
      //Phi
      //Double_t phi=esdtrack->Phi();
      //Double_t pi=TMath::Pi();
      Int_t phibin=PhiBin(phi);
      if(phibin<0) continue;
      if(pt>0.34 && pt<0.5) {
        char *named0PhiallpointrphiSkip =Form("d0phiallpointrphiSkip_%d",phibin);
        char *named0PhiallpointzSkip = Form("d0phiallpointzSkip_%d",phibin);
        char *named0PhipostvtracrphiSkip =Form("d0phipostvtracrphiSkip_%d",phibin);
        char *named0PhipostvtraczSkip = Form("d0phipostvtraczSkip_%d",phibin);
        char *named0PhinegtvtracrphiSkip =Form("d0phinegtvtracrphiSkip_%d",phibin);
        char *named0PhinegtvtraczSkip = Form("d0phinegtvtraczSkip_%d",phibin);
        ((TH1F*)(fOutputphiAllpointSkip->FindObject(named0PhiallpointrphiSkip)))->Fill(10000*dzRecSkip[0]);
        ((TH1F*)(fOutputphiAllpointSkip->FindObject(named0PhiallpointzSkip)))->Fill(10000*dzRecSkip[1]);
        if(charge==+1) {
          ((TH1F*)(fOutputphiPostvtracSkip->FindObject(named0PhipostvtracrphiSkip)))->Fill(10000*dzRecSkip[0]);
          ((TH1F*)(fOutputphiPostvtracSkip->FindObject(named0PhipostvtraczSkip)))->Fill(10000*dzRecSkip[1]);
        }
        if(charge==-1) {
          ((TH1F*)(fOutputphiNegtvtracSkip->FindObject(named0PhinegtvtracrphiSkip)))->Fill(10000*dzRecSkip[0]);
          ((TH1F*)(fOutputphiNegtvtracSkip->FindObject(named0PhinegtvtraczSkip)))->Fill(10000*dzRecSkip[1]);
        }
      }
      
    }

     
  }//end loop over tracks

  //delete esdTrackCuts; esdTrackCuts=NULL;
  //delete primaryVtx; primaryVtx=NULL;
  delete vtxVSkip; vtxVSkip=NULL;
  delete vtxVRec;  vtxVRec=NULL;
  delete vtxESDTrue; vtxESDTrue=NULL;
  PostData(1, fOutputitspureSARec);
  PostData(2, fOutputitspureSASkip);
  PostData(3, fOutputallPointRec);
  PostData(4, fOutputallPointSkip);
  PostData(5, fOutputpartPointRec);
  PostData(6, fOutputpartPointSkip);
  PostData(7, fOutputonepointSPDRec);
  PostData(8, fOutputonepointSPDSkip);
  PostData(9, fOutputpostvTracRec);
  PostData(10, fOutputpostvTracSkip);
  PostData(11, fOutputnegtvTracRec);
  PostData(12, fOutputnegtvTracSkip);
  PostData(13, fOutputpullAllpointRec);
  PostData(14, fOutputpullAllpointSkip);
  PostData(15, fOutputOnlyRefitRec);
  PostData(16, fOutputOnlyRefitSkip);
  PostData(17, fOutputSinThetaRec);
  PostData(18, fOutputSinThetaSkip);
  PostData(19, fOutputallPointTrue);
  PostData(20, fOutputpostvTracTrue);
  PostData(21, fOutputnegtvTracTrue);
  PostData(22, fOutputpullAllpointTrue);
  PostData(23, fOutputphiAllpointSkip);
  PostData(24, fOutputphiPostvtracSkip);
  PostData(25, fOutputphiNegtvtracSkip);
  PostData(26, fOutputparticlePID);
  PostData(27, fOutputPt);
  PostData(28, fNentries);
  PostData(29, fEstimVtx); 
  PostData(30, fOutputWithTrackCuts);

  return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSEImpParRes::PtBin(Double_t pt) const {
  //
  // return the number of the pt bin
  //

  if (pt>0.22 && pt<0.23) return 1;
  if (pt>0.26 && pt<0.27) return 2; 
  if (pt>0.345 && pt<0.355) return 3;
  if (pt>0.45 && pt<0.46) return 4; 
  if (pt>0.55 && pt<0.56) return 5;
  if (pt>0.65 && pt<0.66) return 6; 
  if (pt>0.75 && pt<0.76) return 7; 
  if (pt>0.85 && pt<0.865) return 8; 
  if (pt>1.05 && pt<1.07) return 9;
  if (pt>1.25 && pt<1.30) return 10; 
  if (pt>1.4 && pt<1.55) return 11; 
  if (pt>1.6 && pt<1.8) return 12; 
  if (pt>1.8 && pt<2.0) return 13; 
  if (pt>2.1 && pt<2.3) return 14; 
  if (pt>2.34 && pt<2.64) return 15; 
  if (pt>2.65 && pt<3.0) return 16; 
  if (pt>3.1 && pt<4.) return 17; 
  if (pt>4.1 && pt<5.2) return 18; 
  if (pt>5.3 && pt<6.8)  return 19; 
  if (pt>7.0 && pt<8.8) return 20; 
  if (pt>9. && pt<11.) return 21;
  if (pt>11.1 && pt<14.) return 22; 
  if (pt>14.1  && pt<17.)  return 23;
  if (pt>17.2  && pt<21.8) return 24; 
  if (pt>22.1  && pt<29.)  return 25;
  if (pt>29.05  && pt<35.)  return 26;
  /*
    if (pt>0.22 && pt<0.23) return 1 ;
    if (pt>0.26 && pt<0.27) return 2 ; 
    if (pt>0.35 && pt<0.36) return 3 ;
    if (pt>0.45 && pt<0.46) return 4 ; 
    if (pt>0.55 && pt<0.56) return 5 ;
    if (pt>0.65 && pt<0.66) return 6 ; 
    if (pt>0.75 && pt<0.76) return 7 ; 
    if (pt>0.85 && pt<0.86) return 8 ; 
    if (pt>1.05 && pt<1.06) return 9 ;
    if (pt>1.25 && pt<1.27) return 10; 
    if (pt>1.45 && pt<1.47) return 11; 
    if (pt>1.65 && pt<1.67) return 12; 
    if (pt>1.85 && pt<1.87) return 13; 
    if (pt>2.15 && pt<2.17) return 14; 
    if (pt>2.45 && pt<2.48) return 15; 
    if (pt>2.65 && pt<2.67) return 16; 
    if (pt>2.85 && pt<2.87) return 17; 
    if (pt>3.25 && pt<3.27) return 18; 
    if (pt>3.75 && pt<3.8)  return 19; 
    if (pt>4.15 && pt<4.20) return 20; 
    if (pt>4.95 && pt<5.15) return 21;
    if (pt>5.35 && pt<5.55) return 22; 
    if (pt>6.0  && pt<6.8)  return 23;
    if (pt>8.5  && pt<10.5) return 24; 
    if (pt>12.  && pt<19.)  return 25;
    if (pt>21.  && pt<32.)  return 26;
  */  
  return -1; 
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEImpParRes::Getd0HistRange(Int_t i) const {
  //
  // Return the range of the d0 histograms for each pt bin
  //
  if (i==1) return  2500.; 
  if (i==2) return  1800.;
  if (i==3) return  1750.;
  if (i==4) return  1200.;
  if (i==5) return  1000.;
  if (i==6) return  900.;
  if (i==7) return  850.;
  if (i==8) return  700.;
  if (i==9) return  650.;
  if (i==10) return 600.;
  if (i==11) return 550.;
  if (i==12) return 500.;
  if (i==13) return 450.;
  if (i==14) return 400.;
  if (i==15) return 390.;
  if (i==16) return 380.;
  if (i==17) return 380.;
  if (i==18) return 350.;
  if (i==19) return 320.;
  if (i==20) return 300.;
  if (i==21) return 290.;
  if (i==22) return 270.;
  if (i==23) return 250.;
  if (i==24) return 270.;
  if (i==25) return 279.;
  if (i==26) return 270.;

  return 2000.;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSEImpParRes::SinThetaBin(Double_t sintheta) const {  
  //
  // Return the number of the sinTheta bin
  // 
  if(sintheta>0.7 && sintheta<0.73) return 1;
  if(sintheta>0.73 && sintheta<0.76) return 2; 
  if(sintheta>0.76 && sintheta<0.79) return 3; 
  if(sintheta>0.79 && sintheta<0.82) return 4;
  if(sintheta>0.82 && sintheta<0.85) return 5;       
  if(sintheta>0.85 && sintheta<0.88) return 6; 
  if(sintheta>0.88 && sintheta<0.91) return 7;
  if(sintheta>0.91 && sintheta<0.94) return 8;       
  if(sintheta>0.94 && sintheta<0.97) return 9; 
  if(sintheta>0.97 && sintheta<1.0) return 10;  
  return -1;
}

//___________________________________________________________________________
Int_t AliAnalysisTaskSEImpParRes::PhiBin(Double_t phi) const { 
  Double_t pi=TMath::Pi();
  if(phi>2.*pi || phi<0.) return -1;
  if(phi<0.1*pi) return 1;
  if(phi<0.2*pi) return 2;
  if(phi<0.3*pi) return 3;
  if(phi<0.4*pi) return 4;
  if(phi<0.5*pi) return 5;
  if(phi<0.6*pi) return 6;
  if(phi<0.7*pi) return 7;
  if(phi<0.8*pi) return 8;
  if(phi<0.9*pi) return 9;
  if(phi<1.0*pi) return 10;
  if(phi<1.1*pi) return 11;
  if(phi<1.2*pi) return 12;
  if(phi<1.3*pi) return 13;
  if(phi<1.4*pi) return 14;
  if(phi<1.5*pi) return 15;
  if(phi<1.6*pi) return 16;
  if(phi<1.7*pi) return 17;
  if(phi<1.8*pi) return 18;
  if(phi<1.9*pi) return 19;
  if(phi<2.0*pi) return 20;
  return -1;
}
//___________________________________________________________________________
void AliAnalysisTaskSEImpParRes::Terminate(Option_t */*option*/) {
  //
  // Terminate analysis
  //

  if (fDebug>1) printf("AnalysisTaskSEImpParRes: Terminate() \n");

  return;
}
//__________________________________________________________________________
Int_t AliAnalysisTaskSEImpParRes::ClusterTypeOnITSLayer(AliESDtrack *track,
							Int_t layer) const {
  //
  // Returns cluster type on ITS layer. Returns -1 if no cluster on this layer
  //
  Int_t ctype=-1;

  if(layer<0 || layer>5) return ctype;
  if(!track->HasPointOnITSLayer(layer)) return ctype;
  
  const AliTrackPointArray *array = track->GetTrackPointArray();
  if(!array) {
    //    printf("No tracks points avaialble: check ESDfriends\n");
    return ctype;
  }
  AliTrackPoint point;
  Int_t ipt,volId,modId,layerId;
  for(ipt=0; ipt<array->GetNPoints(); ipt++) {
    array->GetPoint(point,ipt);
    volId = point.GetVolumeID();
    if(volId<=0) continue;
    layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if(layerId==layer+1 && !point.IsExtra()) {
      ctype = point.GetClusterType();
      break;
    }
  }
  return ctype;
}
//---------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEImpParRes::IsSelectedCentrality(AliESDEvent *esd) const
{
  //
  // check if events is in the required multiplicity range
  //

  const AliMultiplicity *alimult = esd->GetMultiplicity();
  Int_t ntrklets=1;
  Int_t nclsSPDouter=0;
  if(alimult) {
    ntrklets = alimult->GetNumberOfTracklets();
    nclsSPDouter = alimult->GetNumberOfITSClusters(1);
  }

  if(nclsSPDouter<fMinMult || nclsSPDouter>fMaxMult) return kFALSE;


  return kTRUE;
}

//----------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEImpParRes::IsTrackSelected(AliVTrack *track, AliVVertex *primary, AliESDtrackCuts *cuts) const{

  if(!cuts) return kTRUE;
  Bool_t retval = kTRUE;
  if(fIsAOD) {
    AliESDtrack esdTrack(track);
    esdTrack.SetTPCClusterMap(((AliAODTrack*)track)->GetTPCClusterMap());
    esdTrack.SetTPCSharedMap(((AliAODTrack*)track)->GetTPCSharedMap());
    esdTrack.SetTPCPointsF(((AliAODTrack*)track)->GetTPCNclsF());
    esdTrack.RelateToVertex((AliESDVertex*)primary,0.,3.);
    if(!cuts->IsSelected(&esdTrack)) retval = kFALSE;
  }
  else {
    AliESDtrack *esdTrack = (AliESDtrack*)track;
    if(!cuts->IsSelected(esdTrack)) retval = kFALSE;
  }
  return retval;
}
