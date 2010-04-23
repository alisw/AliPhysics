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

#include "AliAnalysisManager.h"
#include "AliGeomManager.h"
#include "AliTrackPointArray.h"
#include "AliMCEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"   
#include "AliVertexerTracks.h"
#include "AliAnalysisTaskSEImpParRes.h"

ClassImp(AliAnalysisTaskSEImpParRes)

//________________________________________________________________________
AliAnalysisTaskSEImpParRes::AliAnalysisTaskSEImpParRes():
AliAnalysisTaskSE(),
fReadMC(kFALSE),
fSelectedPdg(-1),
fUseDiamond(kFALSE),
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
fOutputclusterTypeSPD01Skip(0),
fOutputclusterTypeSPD02Skip(0),
fOutputclusterTypeSPD03Skip(0),
fOutputclusterTypeSPD11Skip(0),
fOutputclusterTypeSPD12Skip(0),
fOutputclusterTypeSPD13Skip(0),
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
fReadMC(kFALSE),
fSelectedPdg(-1),
fUseDiamond(kFALSE),
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
fOutputclusterTypeSPD01Skip(0),
fOutputclusterTypeSPD02Skip(0),
fOutputclusterTypeSPD03Skip(0),
fOutputclusterTypeSPD11Skip(0),
fOutputclusterTypeSPD12Skip(0),
fOutputclusterTypeSPD13Skip(0),
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
  DefineOutput(28, TList::Class());  //My private output
  DefineOutput(29, TList::Class());
  DefineOutput(30, TList::Class());  //My private output
  DefineOutput(31, TList::Class());
  DefineOutput(32, TList::Class());  //My private output
  DefineOutput(33, TH1F::Class());
  DefineOutput(34, TH1F::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEImpParRes::~AliAnalysisTaskSEImpParRes()
{
  //
  // default distructor  
  // 
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
  if (fOutputphiPostvtracSkip)             {delete fOutputphiPostvtracSkip;fOutputphiPostvtracSkip=0x0;}
  if (fOutputphiNegtvtracSkip)            {delete fOutputphiNegtvtracSkip;fOutputphiNegtvtracSkip=0x0;}
  if (fOutputclusterTypeSPD01Skip){delete fOutputclusterTypeSPD01Skip;fOutputclusterTypeSPD01Skip=0x0;}
  if (fOutputclusterTypeSPD02Skip){delete fOutputclusterTypeSPD02Skip;fOutputclusterTypeSPD02Skip=0x0;}
  if (fOutputclusterTypeSPD03Skip){delete fOutputclusterTypeSPD03Skip;fOutputclusterTypeSPD03Skip=0x0;}
  if (fOutputclusterTypeSPD11Skip){delete fOutputclusterTypeSPD11Skip;fOutputclusterTypeSPD11Skip=0x0;}
  if (fOutputclusterTypeSPD12Skip){delete fOutputclusterTypeSPD12Skip;fOutputclusterTypeSPD12Skip=0x0;}
  if (fOutputclusterTypeSPD13Skip){delete fOutputclusterTypeSPD13Skip;fOutputclusterTypeSPD13Skip=0x0;}
  if (fOutputPt)                                                      {delete fOutputPt;fOutputPt=0x0;}
  if (fNentries)                                           { delete fNentries;     fNentries    =0x0; }
  if (fEstimVtx)                                           { delete fEstimVtx;     fEstimVtx    =0x0; }

} 

//________________________________________________________________________
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
  
  if (!fOutputclusterTypeSPD01Skip) {
    fOutputclusterTypeSPD01Skip = new TList();
    fOutputclusterTypeSPD01Skip->SetOwner();
    fOutputclusterTypeSPD01Skip->SetName("clustertypeSPD01Skip");
  }
  
  
  if (!fOutputclusterTypeSPD02Skip) {
    fOutputclusterTypeSPD02Skip = new TList();
    fOutputclusterTypeSPD02Skip->SetOwner();
    fOutputclusterTypeSPD02Skip->SetName("clustertypeSPD02Skip");
  }
  
  if (!fOutputclusterTypeSPD03Skip) {
    fOutputclusterTypeSPD03Skip = new TList();
    fOutputclusterTypeSPD03Skip->SetOwner();
    fOutputclusterTypeSPD03Skip->SetName("clustertypeSPD03Skip");
  }
  
  if (!fOutputclusterTypeSPD11Skip) {
    fOutputclusterTypeSPD11Skip = new TList();
    fOutputclusterTypeSPD11Skip->SetOwner();
    fOutputclusterTypeSPD11Skip->SetName("clustertypeSPD11Skip");
  }
  
  if (!fOutputclusterTypeSPD12Skip) {
    fOutputclusterTypeSPD12Skip = new TList();
    fOutputclusterTypeSPD12Skip->SetOwner();
    fOutputclusterTypeSPD12Skip->SetName("clustertypeSPD12Skip");
  }
  
  if (!fOutputclusterTypeSPD13Skip) {
    fOutputclusterTypeSPD13Skip = new TList();
    fOutputclusterTypeSPD13Skip->SetOwner();
    fOutputclusterTypeSPD13Skip->SetName("clustertypeSPD13Skip");
  }

  if (!fOutputPt) {
    fOutputPt = new TList();
    fOutputPt->SetOwner();
    fOutputPt->SetName("Pt");
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
  const TString d0clusterTypeSPD01rphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD01zTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries"; 
  const TString d0clusterTypeSPD02rphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD02zTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD03rphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD03zTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD11rphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD11zTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD12rphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD12zTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD13rphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0clusterTypeSPD13zTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0rphiTitle = "d_{0} Distribution_rphi; d_{0} [#mum]; Entries";
  const TString d0zTitle  = "d_{0} Distribution_z; d_{0} [#mum]; Entries";

  TString  named0itspureSArphiRec,named0itspureSAzRec,named0allpointrphiRec, named0allpointzRec,named0partpointrphiRec, named0partpointzRec,named0onepointSPDrphiRec, named0onepointSPDzRec,named0postvtracrphiRec, named0postvtraczRec,named0negtvtracrphiRec, named0negtvtraczRec,named0pt,named0pullAllpointrphiRec,named0pullAllpointzRec,named0onlyRefitrphiRec,named0onlyRefitzRec;
 
  TH1F *d0ITSpureSArphiRec=0,*d0ITSpureSAzRec=0,*d0AllpointrphiRec=0, *d0AllpointzRec=0,*d0PartpointrphiRec=0, *d0PartpointzRec=0,
    *d0OnepointSPDrphiRec=0,*d0OnepointSPDzRec=0,*d0PostvtracrphiRec=0, *d0PostvtraczRec=0,*d0NegtvtracrphiRec=0, *d0NegtvtraczRec=0,*d0Pt=0,*d0PullAllpointrphiRec=0,*d0PullAllpointzRec=0,*d0OnlyRefitrphiRec=0,*d0OnlyRefitzRec=0;

  TString  named0itspureSArphiSkip,named0itspureSAzSkip,named0allpointrphiSkip, named0allpointzSkip,named0partpointrphiSkip, named0partpointzSkip,named0onepointSPDrphiSkip, named0onepointSPDzSkip,named0postvtracrphiSkip, named0postvtraczSkip,named0negtvtracrphiSkip, named0negtvtraczSkip,named0ptSkip,named0pullAllpointrphiSkip,named0pullAllpointzSkip,named0onlyRefitrphiSkip,named0onlyRefitzSkip,named0allpointrphiTrue, named0allpointzTrue,named0postvtracrphiTrue, named0postvtraczTrue,named0negtvtracrphiTrue, named0negtvtraczTrue,named0pullAllpointrphiTrue,named0pullAllpointzTrue,named0clusterTypeSPD01rphiSkip,named0clusterTypeSPD01zSkip,named0clusterTypeSPD02rphiSkip,named0clusterTypeSPD02zSkip,named0clusterTypeSPD03rphiSkip,named0clusterTypeSPD03zSkip,named0clusterTypeSPD11rphiSkip,named0clusterTypeSPD11zSkip,named0clusterTypeSPD12rphiSkip,named0clusterTypeSPD12zSkip,named0clusterTypeSPD13rphiSkip,named0clusterTypeSPD13zSkip;

  TH1F *d0ITSpureSArphiSkip=0,*d0ITSpureSAzSkip=0,*d0AllpointrphiSkip=0, *d0AllpointzSkip=0,*d0PartpointrphiSkip=0, *d0PartpointzSkip=0,*d0OnepointSPDrphiSkip=0,*d0OnepointSPDzSkip=0,*d0PostvtracrphiSkip=0, *d0PostvtraczSkip=0,*d0NegtvtracrphiSkip=0,*d0NegtvtraczSkip=0,*d0PullAllpointrphiSkip=0,*d0PullAllpointzSkip=0,*d0OnlyRefitrphiSkip=0,*d0OnlyRefitzSkip=0,*d0AllpointrphiTrue=0, *d0AllpointzTrue=0,*d0PostvtracrphiTrue=0, *d0PostvtraczTrue=0,*d0NegtvtracrphiTrue=0,*d0NegtvtraczTrue=0,*d0PullAllpointrphiTrue,*d0PullAllpointzTrue,*d0ClustertypeSPD01rphiSkip=0,*d0ClustertypeSPD01zSkip=0,*d0ClustertypeSPD02rphiSkip=0,*d0ClustertypeSPD02zSkip=0,*d0ClustertypeSPD03rphiSkip=0,*d0ClustertypeSPD03zSkip=0,*d0ClustertypeSPD11rphiSkip=0,*d0ClustertypeSPD11zSkip=0,*d0ClustertypeSPD12rphiSkip=0,*d0ClustertypeSPD12zSkip=0,*d0ClustertypeSPD13rphiSkip=0,*d0ClustertypeSPD13zSkip=0;

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
    d0ITSpureSArphiSkip = new TH1F(named0itspureSArphiSkip.Data(), d0ITSpureSArphiTitle.Data(), 400, -Getd0HistRange(i), Getd0HistRange(i));
    d0ITSpureSArphiSkip->Sumw2();
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

    named0clusterTypeSPD01rphiSkip = "d0clustertypeSPD01rphiSkip_";
    named0clusterTypeSPD01rphiSkip +=i;
    named0clusterTypeSPD01zSkip = "d0clustertypeSPD01zSkip_";
    named0clusterTypeSPD01zSkip +=i;
    d0ClustertypeSPD01rphiSkip = new TH1F(named0clusterTypeSPD01rphiSkip.Data(),d0clusterTypeSPD01rphiTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD01rphiSkip->Sumw2();
    d0ClustertypeSPD01rphiSkip->SetMinimum(0);
    fOutputclusterTypeSPD01Skip->Add(d0ClustertypeSPD01rphiSkip);
    d0ClustertypeSPD01zSkip = new TH1F(named0clusterTypeSPD01zSkip.Data(),d0clusterTypeSPD01zTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD01zSkip->Sumw2();
    d0ClustertypeSPD01zSkip->SetMinimum(0);
    fOutputclusterTypeSPD01Skip->Add(d0ClustertypeSPD01zSkip);

    named0clusterTypeSPD02rphiSkip = "d0clustertypeSPD02rphiSkip_";
    named0clusterTypeSPD02rphiSkip +=i;
    named0clusterTypeSPD02zSkip = "d0clustertypeSPD02zSkip_";
    named0clusterTypeSPD02zSkip +=i;
    d0ClustertypeSPD02rphiSkip = new TH1F(named0clusterTypeSPD02rphiSkip.Data(),d0clusterTypeSPD02rphiTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD02rphiSkip->Sumw2();
    d0ClustertypeSPD02rphiSkip->SetMinimum(0);
    fOutputclusterTypeSPD02Skip->Add(d0ClustertypeSPD02rphiSkip);
    d0ClustertypeSPD02zSkip = new TH1F(named0clusterTypeSPD02zSkip.Data(),d0clusterTypeSPD02zTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD02zSkip->Sumw2();
    d0ClustertypeSPD02zSkip->SetMinimum(0);
    fOutputclusterTypeSPD02Skip->Add(d0ClustertypeSPD02zSkip);

    named0clusterTypeSPD03rphiSkip = "d0clustertypeSPD03rphiSkip_";
    named0clusterTypeSPD03rphiSkip +=i;
    named0clusterTypeSPD03zSkip = "d0clustertypeSPD03zSkip_";
    named0clusterTypeSPD03zSkip +=i;
    d0ClustertypeSPD03rphiSkip = new TH1F(named0clusterTypeSPD03rphiSkip.Data(),d0clusterTypeSPD03rphiTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD03rphiSkip->Sumw2();
    d0ClustertypeSPD03rphiSkip->SetMinimum(0);
    fOutputclusterTypeSPD03Skip->Add(d0ClustertypeSPD03rphiSkip);
    d0ClustertypeSPD03zSkip = new TH1F(named0clusterTypeSPD03zSkip.Data(),d0clusterTypeSPD03zTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD03zSkip->Sumw2();
    d0ClustertypeSPD03zSkip->SetMinimum(0);
    fOutputclusterTypeSPD03Skip->Add(d0ClustertypeSPD03zSkip);

    named0clusterTypeSPD11rphiSkip = "d0clustertypeSPD11rphiSkip_";
    named0clusterTypeSPD11rphiSkip +=i;
    named0clusterTypeSPD11zSkip = "d0clustertypeSPD11zSkip_";
    named0clusterTypeSPD11zSkip +=i;
    d0ClustertypeSPD11rphiSkip = new TH1F(named0clusterTypeSPD11rphiSkip.Data(),d0clusterTypeSPD11rphiTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD11rphiSkip->Sumw2();
    d0ClustertypeSPD11rphiSkip->SetMinimum(0);
    fOutputclusterTypeSPD11Skip->Add(d0ClustertypeSPD11rphiSkip);
    d0ClustertypeSPD11zSkip = new TH1F(named0clusterTypeSPD11zSkip.Data(),d0clusterTypeSPD11zTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD11zSkip->Sumw2();
    d0ClustertypeSPD11zSkip->SetMinimum(0);
    fOutputclusterTypeSPD11Skip->Add(d0ClustertypeSPD11zSkip);

    named0clusterTypeSPD12rphiSkip = "d0clustertypeSPD12rphiSkip_";
    named0clusterTypeSPD12rphiSkip +=i;
    named0clusterTypeSPD12zSkip = "d0clustertypeSPD12zSkip_";
    named0clusterTypeSPD12zSkip +=i;
    d0ClustertypeSPD12rphiSkip = new TH1F(named0clusterTypeSPD12rphiSkip.Data(),d0clusterTypeSPD12rphiTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD12rphiSkip->Sumw2();
    d0ClustertypeSPD12rphiSkip->SetMinimum(0);
    fOutputclusterTypeSPD12Skip->Add(d0ClustertypeSPD12rphiSkip);
    d0ClustertypeSPD12zSkip = new TH1F(named0clusterTypeSPD12zSkip.Data(),d0clusterTypeSPD12zTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD12zSkip->Sumw2();
    d0ClustertypeSPD12zSkip->SetMinimum(0);
    fOutputclusterTypeSPD12Skip->Add(d0ClustertypeSPD12zSkip);

    named0clusterTypeSPD13rphiSkip = "d0clustertypeSPD13rphiSkip_";
    named0clusterTypeSPD13rphiSkip +=i;
    named0clusterTypeSPD13zSkip = "d0clustertypeSPD13zSkip_";
    named0clusterTypeSPD13zSkip +=i;
    d0ClustertypeSPD13rphiSkip = new TH1F(named0clusterTypeSPD13rphiSkip.Data(),d0clusterTypeSPD13rphiTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD13rphiSkip->Sumw2();
    d0ClustertypeSPD13rphiSkip->SetMinimum(0);
    fOutputclusterTypeSPD13Skip->Add(d0ClustertypeSPD13rphiSkip);
    d0ClustertypeSPD13zSkip = new TH1F(named0clusterTypeSPD13zSkip.Data(),d0clusterTypeSPD13zTitle.Data(),400,-10.,10.);
    d0ClustertypeSPD13zSkip->Sumw2();
    d0ClustertypeSPD13zSkip->SetMinimum(0);
    fOutputclusterTypeSPD13Skip->Add(d0ClustertypeSPD13zSkip);

    named0pt = "d0pt_";
    named0pt += i;
    d0Pt = new TH1F(named0pt.Data(), d0ptTitle.Data(), 100, 0, 35.);
    d0Pt->Sumw2();
    d0Pt->SetMinimum(0);  
    fOutputPt->Add(d0Pt);
  }

  const Int_t nhistm=10;
  for(Int_t i=0; i<=nhistm; i++) {
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
    TString named0sinThetaonerphiRec,named0sinThetaonezRec,named0sinThetatworphiRec,named0sinThetatwozRec,named0sinThetathreerphiRec,named0sinThetathreezRec,named0sinThetafourrphiRec,named0sinThetafourzRec;
    
    TH1F *d0SinThetaonerphiRec,*d0SinThetaonezRec,*d0SinThetatworphiRec,*d0SinThetatwozRec,*d0SinThetathreerphiRec,*d0SinThetathreezRec,*d0SinThetafourrphiRec,*d0SinThetafourzRec;
    
    TString  named0sinThetaonerphiSkip,named0sinThetaonezSkip,named0sinThetatworphiSkip,named0sinThetatwozSkip,named0sinThetathreerphiSkip,named0sinThetathreezSkip,named0sinThetafourrphiSkip,named0sinThetafourzSkip,named0phiAllpointrphiSkip, named0phiAllpointzSkip,named0phiPostvtracrphiSkip, named0phiPostvtraczSkip,named0phiNegtvtracrphiSkip,named0phiNegtvtraczSkip;
    
    TH1F*d0SinThetaonerphiSkip,*d0SinThetaonezSkip,*d0SinThetatworphiSkip,*d0SinThetatwozSkip,*d0SinThetathreerphiSkip,*d0SinThetathreezSkip,*d0SinThetafourrphiSkip,*d0SinThetafourzSkip, *d0PhiAllpointrphiSkip,*d0PhiAllpointzSkip,*d0PhiPostvtracrphiSkip,*d0PhiPostvtraczSkip,*d0PhiNegtvtracrphiSkip,*d0PhiNegtvtraczSkip;

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
  if(!fEstimVtx) fEstimVtx = new TH1F("vtxRes","Resolution of vertex",400,-1000.,1000);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEImpParRes::UserExec(Option_t */*option*/)
{
  //
  // Track selection and filling of d0 histograms
  //
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliError("ESD event not found. Nothing done!");
    return;
  }

  Float_t diamondcovxy[3];
  esd->GetDiamondCovXY(diamondcovxy);
  Double_t pos[3]={esd->GetDiamondX(),esd->GetDiamondY(),0.};
  //diamondcovxy[0]=0.05*0.05;
  //diamondcovxy[2]=0.05*0.05;
  Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.};
  AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);

  AliVertexerTracks *vertexer0 = new AliVertexerTracks(esd->GetMagneticField());
  vertexer0->SetITSMode();
  vertexer0->SetMinClusters(4);  
  if(fUseDiamond) vertexer0->SetVtxStart(diamond);
  AliESDVertex *vtxESDRec = (AliESDVertex*)vertexer0->FindPrimaryVertex(esd);
  delete vertexer0; vertexer0=0;
  if(vtxESDRec->GetNContributors()<1) return;

  AliStack *stack=0;
  AliESDVertex *vtxESDTrue=0;
   if (fReadMC) { 
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

    //load MC header for ESD;//see $ALICE_ROOT/PWG1/global/AliAnalysisTaskSEVertexESD.cxx
    AliHeader *mcHeader = eventHandler->MCEvent()->Header();
    if (!mcHeader) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }

    AliGenEventHeader* genHeader = mcHeader->GenEventHeader();
    Double_t vtxTrue[3];
    TArrayF mcVertex(3);
    mcVertex[0]=9999.; mcVertex[1]=9999.; mcVertex[2]=9999.;
    genHeader->PrimaryVertex(mcVertex);
    vtxTrue[0]=mcVertex[0];vtxTrue[1]=mcVertex[1];vtxTrue[2]=mcVertex[2];
    Double_t sigmaTrue[3]={0., 0., 0.,};
    //mcHeader->GetVertex(vtxTrue);//note the vtxTrue is void here,so must need the next line.
    //AliESDVertex *vtxESDTrue = new AliESDVertex(vtxTrue,sigmaTrue);
    vtxESDTrue = new AliESDVertex(vtxTrue,sigmaTrue);
    }

  Int_t nTrks = esd->GetNumberOfTracks();
  for (Int_t it=0; it<nTrks; it++) {  // loop over tracks
    // read track
    AliESDtrack *esdtrack = esd->GetTrack(it);

    // ask for ITS refit
    if (!(esdtrack->GetStatus()&AliESDtrack::kITSrefit)) {
      continue;
    }
    
    Int_t pdgCode=0; 
    if(fReadMC && stack) { 
      Int_t trkLabel = esdtrack->GetLabel(); 
      //Int_t trkLabel = TMath::Abs(esdtrack->GetLabel()); 
      if(trkLabel<0) continue;
      TParticle  *part = (TParticle*)stack->Particle(trkLabel); 
      pdgCode = TMath::Abs(part->GetPdgCode()); 
      //printf("pdgCode===%d\n", pdgCode);
      if(fSelectedPdg>0 && pdgCode!=fSelectedPdg) continue;
    }


    Int_t npointsITS=0,npointsSPD=0;
    for (Int_t i=0;i<6;i++){ 
      if (i<2 && esdtrack->HasPointOnITSLayer(i)) npointsSPD +=1;
      if (esdtrack->HasPointOnITSLayer(i)) npointsITS +=1;  
    }

    //Get specific primary vertex--Reconstructed primary vertex do not include the track considering.
    AliVertexerTracks *vertexer = new AliVertexerTracks(esd->GetMagneticField());
    vertexer->SetITSMode();
    vertexer->SetMinClusters(4);
    if(fUseDiamond) vertexer->SetVtxStart(diamond);
    Int_t skipped[2];
    skipped[0] = (Int_t)esdtrack->GetID();
    vertexer->SetSkipTracks(1,skipped);      
    AliESDVertex *vtxESDSkip = (AliESDVertex*)vertexer->FindPrimaryVertex(esd); 
    delete vertexer; vertexer=NULL;
    if(vtxESDSkip->GetNContributors()<1) continue;
     

    Double_t dzRec[2], covdzRec[3], dzRecSkip[2], covdzRecSkip[3],dzTrue[2], covdzTrue[3];
    Double_t beampiperadius=3.;
    Double_t pt= esdtrack->Pt();
    Int_t bin = PtBin(pt);

    if(bin==-1) {
      delete vtxESDSkip; vtxESDSkip = 0x0;
      continue;
    }

    esdtrack->PropagateToDCA(vtxESDRec, esd->GetMagneticField(), beampiperadius, dzRec, covdzRec);
    esdtrack->PropagateToDCA(vtxESDSkip, esd->GetMagneticField(), beampiperadius, dzRecSkip, covdzRecSkip);
    if(fReadMC) esdtrack->PropagateToDCA(vtxESDTrue, esd->GetMagneticField(), beampiperadius, dzTrue, covdzTrue);
    if(covdzRec[0]<0 || covdzRec[2]<0 || covdzRecSkip[0]<0 || covdzRecSkip[2]<0) {
      delete vtxESDSkip; vtxESDSkip = 0x0;
      continue;
    }
      
    //printf("Pt: %f GeV/c; Impact parameter: rphi %f cm  z %f cm\n", pt, dzRec[0], dzRec[1]);

    /*
    // RUBEN'S METHOD, NOT FULLY IMPLEMENTED YET    
    Double_t      fIPCenIni[3], xyzDCA[3];
    for (int i=3;i--;) fIPCenIni[i] = 0.;
    //Int_t nTracks = nTrks - 1;
    esdtrack->GetXYZ(xyzDCA);
    //double pTrack = esdtrack->GetP();
    double phiTrack = esdtrack->Phi();
    double cs = TMath::Cos(phiTrack);
    double sn = TMath::Sin(phiTrack);
    double trDCA = (xyzDCA[0]-fIPCenIni[0])         *sn - (xyzDCA[1]-fIPCenIni[1])         *cs;  // track signed DCA to origin
    double vtDCA = (vtxESDSkip->GetXv()-fIPCenIni[0])*sn - (vtxESDSkip->GetYv()-fIPCenIni[1])*cs;  // vertex signed DCA to origin
    
    // update the estimator values
    //double estIP  = rvD*rtD;
    double estVtx = vtDCA*( vtDCA- trDCA);
    //double estTrc = rtD*(rtD - rvD);
    //
    //if (nTracks >= fMinTracksForIP) fEstimIP->Fill(phiTrack, estIP);
    fEstimVtx->Fill(10000*estVtx);
    //if (pTrack<1e-6) pTrack = GetTrackMinP()+1e6;
    //fEstimTrc->Fill(1./pTrack,estTrc);
    */


    //fill the histgram with all particle 
    //-------------------------------------------1----------------------------------------------
    
    /*
    TString named0AllrphiRec=" ",named0AllzRec =" ";//named0AllrphiTrue=" ",named0AllzTrue =" ";
      //named0AllrphiSkip=" ",named0AllzSkip =" ";
      
      named0AllrphiRec ="d0allrphiRec_";
      named0AllrphiRec += bin;
      named0AllzRec = "d0allzRec_";
      named0AllzRec += bin;
      ((TH1F*)(fOutputallRec->FindObject(named0AllrphiRec)))->Fill(10000*dzRec[0]);
      ((TH1F*)(fOutputallRec->FindObject(named0AllzRec)))->Fill(10000*dzRec[1]);
    */

    //-------------------------------------------2----------------------------------------------
    //TSting named0AllrphiRec = "d0allrphiRec_" + bin;
    //TSting named0AllzRec = "d0allzRec_" + bin;
    //((TH1F*)(fOutputallRec->FindObject(named0AllrphiRec.Data())))->Fill(10000.*dzRec[0]);
    //((TH1F*)(fOutputallRec->FindObject(named0AllzRec.Data())))->Fill(10000.*dzRec[1]);

    //-------------------------------------------3------------------------------------------------


    // ITS standalone
    if (!(esdtrack->GetStatus()&AliESDtrack::kTPCrefit) && (esdtrack->GetStatus()&AliESDtrack::kITSpureSA)
	&& npointsSPD>0 && npointsITS>=4) {
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
    if (!(esdtrack->GetStatus()&AliESDtrack::kTPCrefit) || esdtrack->GetNcls(1)<70) {
      delete vtxESDSkip; vtxESDSkip = 0x0;
      continue;
    }

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

    // with 6 ITS points (including different selection)
    if(npointsITS==6) {

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
      ((TH1F*)(fOutputallPointTrue->FindObject(named0AllpointrphiTrue)))->Fill(10000.*dzTrue[0]);
      ((TH1F*)(fOutputallPointTrue->FindObject(named0AllpointzTrue)))->Fill(10000.*dzTrue[1]);
      
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
      ((TH1F*)(fOutputpullAllpointTrue->FindObject(named0PullAllpointrphiTrue)))->Fill(dzTrue[0]/TMath::Sqrt(covdzTrue[0]));
      ((TH1F*)(fOutputpullAllpointTrue->FindObject(named0PullAllpointzTrue)))->Fill(dzTrue[1]/TMath::Sqrt(covdzTrue[2]));
      
      //postive and negative track
      Short_t charge=esdtrack->Charge();
      if(charge==1 ) {
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
	((TH1F*)(fOutputpostvTracTrue->FindObject(named0PostvtracrphiTrue)))->Fill(10000.*dzTrue[0]);
	((TH1F*)(fOutputpostvTracTrue->FindObject(named0PostvtraczTrue)))->Fill(10000.*dzTrue[1]);
      }
      
      if(charge==-1 ) {
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
	((TH1F*)(fOutputnegtvTracTrue->FindObject(named0NegtvtracrphiTrue)))->Fill(10000.*dzTrue[0]);
	((TH1F*)(fOutputnegtvTracTrue->FindObject(named0NegtvtraczTrue)))->Fill(10000.*dzTrue[1]);	
      }    
      
      // SinTheta 
      Double_t theta=esdtrack->Theta(); 
      Double_t Sintheta=TMath::Sin(theta);
      Int_t thetabin = SinThetaBin(Sintheta);
      if(thetabin<0) {delete vtxESDSkip;continue;}
      
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
      Double_t phi=esdtrack->Phi(); 
      //Double_t pi=TMath::Pi();
      Int_t phibin=PhiBin(phi);
      if(phibin<0) {delete vtxESDSkip;continue;}
      if(pt>0.4 && pt<0.6){
	char *named0PhiallpointrphiSkip =Form("d0phiallpointrphiSkip_%d",phibin);
	char *named0PhiallpointzSkip = Form("d0phiallpointzSkip_%d",phibin);
	char *named0PhipostvtracrphiSkip =Form("d0phipostvtracrphiSkip_%d",phibin);
	char *named0PhipostvtraczSkip = Form("d0phipostvtraczSkip_%d",phibin);
	char *named0PhinegtvtracrphiSkip =Form("d0phinegtvtracrphiSkip_%d",phibin);
	char *named0PhinegtvtraczSkip = Form("d0phinegtvtraczSkip_%d",phibin);
	((TH1F*)(fOutputphiAllpointSkip->FindObject(named0PhiallpointrphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputphiAllpointSkip->FindObject(named0PhiallpointzSkip)))->Fill(10000*dzRecSkip[1]);
	((TH1F*)(fOutputphiPostvtracSkip->FindObject(named0PhipostvtracrphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputphiPostvtracSkip->FindObject(named0PhipostvtraczSkip)))->Fill(10000*dzRecSkip[1]);
	((TH1F*)(fOutputphiNegtvtracSkip->FindObject(named0PhinegtvtracrphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputphiNegtvtracSkip->FindObject(named0PhinegtvtraczSkip)))->Fill(10000*dzRecSkip[1]);
      }
      
      //cluster type
      Int_t ctypeSPD0 = ClusterTypeOnITSLayer(esdtrack,0);
      if(ctypeSPD0==-1) {delete vtxESDSkip;continue;}
      Int_t ctypeSPD1 = ClusterTypeOnITSLayer(esdtrack,1);
      if(ctypeSPD1==-1) {delete vtxESDSkip;continue;}

      if(ctypeSPD0==1) {
	char *named0ClustertypeSPD01rphiSkip = Form("d0clustertypeSPD01rphiSkip_%d",phibin);
	char *named0ClustertypeSPD01zSkip = Form("d0clustertypeSPD01rphiSkip_%d",phibin);
	((TH1F*)(fOutputclusterTypeSPD01Skip->FindObject(named0ClustertypeSPD01rphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputclusterTypeSPD01Skip->FindObject(named0ClustertypeSPD01zSkip)))->Fill(10000*dzRecSkip[1]);
      }
      
      if(ctypeSPD0==2) {
	char *named0ClustertypeSPD02rphiSkip = Form("d0clustertypeSPD02rphiSkip_%d",phibin);
	char *named0ClustertypeSPD02zSkip = Form("d0clustertypeSPD02rphiSkip_%d",phibin);
	((TH1F*)(fOutputclusterTypeSPD02Skip->FindObject(named0ClustertypeSPD02rphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputclusterTypeSPD02Skip->FindObject(named0ClustertypeSPD02zSkip)))->Fill(10000*dzRecSkip[1]);
      }
      
      if(ctypeSPD0==3) {
	char *named0ClustertypeSPD03rphiSkip = Form("d0clustertypeSPD03rphiSkip_%d",phibin);
	char *named0ClustertypeSPD03zSkip = Form("d0clustertypeSPD03zSkip_%d",phibin);
	((TH1F*)(fOutputclusterTypeSPD03Skip->FindObject(named0ClustertypeSPD03rphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputclusterTypeSPD03Skip->FindObject(named0ClustertypeSPD03zSkip)))->Fill(10000*dzRecSkip[1]);
      }
      
      if(ctypeSPD1==1){
	char *named0ClustertypeSPD11rphiSkip = Form("d0clustertypeSPD11rphiSkip_%d",phibin);
	char *named0ClustertypeSPD11zSkip = Form("d0clustertypeSPD11rphiSkip_%d",phibin);
	((TH1F*)(fOutputclusterTypeSPD11Skip->FindObject(named0ClustertypeSPD11rphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputclusterTypeSPD11Skip->FindObject(named0ClustertypeSPD11zSkip)))->Fill(10000*dzRecSkip[1]);
      } 
      if(ctypeSPD1==2){
	char *named0ClustertypeSPD12rphiSkip = Form("d0clustertypeSPD12rphiSkip_%d",phibin);
	char *named0ClustertypeSPD12zSkip = Form("d0clustertypeSPD12rphiSkip_%d",phibin);
	((TH1F*)(fOutputclusterTypeSPD12Skip->FindObject(named0ClustertypeSPD12rphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputclusterTypeSPD12Skip->FindObject(named0ClustertypeSPD12zSkip)))->Fill(10000*dzRecSkip[1]);
      } 
      
      if(ctypeSPD1==3){
	char *named0ClustertypeSPD13rphiSkip = Form("d0clustertypeSPD13rphiSkip_%d",phibin);
	char *named0ClustertypeSPD13zSkip = Form("d0clustertypeSPD13rphiSkip_%d",phibin);
	((TH1F*)(fOutputclusterTypeSPD13Skip->FindObject(named0ClustertypeSPD13rphiSkip)))->Fill(10000*dzRecSkip[0]);
	((TH1F*)(fOutputclusterTypeSPD13Skip->FindObject(named0ClustertypeSPD13zSkip)))->Fill(10000*dzRecSkip[1]);
      }  
    }
    
    delete vtxESDSkip; vtxESDSkip = 0x0;
  }  // end loop of tracks
  
  delete diamond; diamond=NULL;
  delete vtxESDRec; vtxESDRec = 0x0;
  delete vtxESDTrue;vtxESDTrue = 0x0;
  fNentries->Fill(1);
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
  PostData(26, fOutputclusterTypeSPD01Skip);
  PostData(27, fOutputclusterTypeSPD02Skip);
  PostData(28, fOutputclusterTypeSPD03Skip);
  PostData(29, fOutputclusterTypeSPD11Skip);
  PostData(30, fOutputclusterTypeSPD12Skip);
  PostData(31, fOutputclusterTypeSPD13Skip);
  PostData(32, fOutputPt);
  PostData(33, fNentries);
  PostData(34, fEstimVtx);
  
  return;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSEImpParRes::PtBin(Double_t pt) const {
  //
  // return the number of the pt bin
  //

  if (pt>0.22 && pt<0.23) return 1;
  if (pt>0.26 && pt<0.27) return 2; 
  if (pt>0.35 && pt<0.36) return 3;
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
  if(phi<0.2*pi) return 1;
  if(phi<0.4*pi) return 2;
  if(phi<0.6*pi) return 3;
  if(phi<0.8*pi) return 4;
  if(phi<1.0*pi) return 5;
  if(phi<1.2*pi) return 6;
  if(phi<1.4*pi) return 7;
  if(phi<1.6*pi) return 8;
  if(phi<1.8*pi) return 9;
  if(phi<2.*pi) return 10;
  return -1;
}

//___________________________________________________________________________
void AliAnalysisTaskSEImpParRes::Terminate(Option_t */*option*/) {
  //
  // Terminate analysis
  //

  if (fDebug>1) printf("AnalysisTaskSEImpParRes: Terminate() \n");

  fOutputitspureSARec = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutputitspureSARec) printf("ERROR: fOutputitspureSARec not available\n");

  fOutputitspureSASkip = dynamic_cast<TList*>(GetOutputData(2));
  if (!fOutputitspureSASkip) printf("ERROR: fOutputitspureSASkip not available\n");

  fOutputallPointRec = dynamic_cast<TList*>(GetOutputData(3));
  if (!fOutputallPointRec) printf("ERROR: fOutputallPointRec not available\n");

  fOutputallPointSkip = dynamic_cast<TList*>(GetOutputData(4));
  if (!fOutputallPointSkip) printf("ERROR: fOutputallPointSkip not available\n");

  fOutputpartPointRec = dynamic_cast<TList*>(GetOutputData(5));
  if (!fOutputpartPointRec) printf("ERROR: fOutputpartPointRec not available\n");

  fOutputpartPointSkip = dynamic_cast<TList*>(GetOutputData(6));
  if (!fOutputpartPointSkip) printf("ERROR: fOutputpartPointSkip not available\n");

  fOutputonepointSPDRec = dynamic_cast<TList*>(GetOutputData(7));
  if (!fOutputonepointSPDRec) printf("ERROR: fOutputonepointSPDSRec not available\n");

  fOutputonepointSPDSkip = dynamic_cast<TList*>(GetOutputData(8));
  if (!fOutputonepointSPDSkip) printf("ERROR: fOutputonepointSPDSkip not available\n");

  fOutputpostvTracRec = dynamic_cast<TList*>(GetOutputData(9));
  if (!fOutputpostvTracRec) printf("ERROR: fOutputpostvTracRec not available\n");

  fOutputpostvTracSkip = dynamic_cast<TList*>(GetOutputData(10));
  if (!fOutputpostvTracSkip) printf("ERROR: fOutputpostvTracSkip not available\n");

  fOutputnegtvTracRec = dynamic_cast<TList*>(GetOutputData(11));
  if (!fOutputnegtvTracRec) printf("ERROR: fOutputnegtvTracRec not available\n");

  fOutputnegtvTracSkip = dynamic_cast<TList*>(GetOutputData(12));
  if (!fOutputnegtvTracSkip) printf("ERROR: fOutputnegtvTracSkip not available\n");


  fOutputpullAllpointRec = dynamic_cast<TList*>(GetOutputData(13));
  if (!fOutputpullAllpointRec) printf("ERROR: fOutputpullAllpointRec not available\n");

  fOutputpullAllpointSkip = dynamic_cast<TList*>(GetOutputData(14));
  if (!fOutputpullAllpointSkip) printf("ERROR: fOutputpullAllpointSkip not available\n");

  fOutputOnlyRefitRec = dynamic_cast<TList*>(GetOutputData(15));
  if (!fOutputOnlyRefitRec) printf("ERROR: fOutputOnlyRefitRec not available\n");

  fOutputOnlyRefitSkip = dynamic_cast<TList*>(GetOutputData(16));
  if (!fOutputOnlyRefitSkip) printf("ERROR: fOutputOnlyRefitSkip not available\n");
 
  fOutputSinThetaRec = dynamic_cast<TList*>(GetOutputData(17));
  if (!fOutputSinThetaRec) printf("ERROR:fOutputSinThetaRec not available\n");

  fOutputSinThetaSkip = dynamic_cast<TList*>(GetOutputData(18));
  if (!fOutputSinThetaSkip) printf("ERROR:fOutputSinThetaSkip not available\n");
  
  fOutputallPointTrue = dynamic_cast<TList*>(GetOutputData(19));
  if (!fOutputallPointTrue) printf("ERROR:fOutputallPointTrue not available\n");
  
  fOutputpostvTracTrue = dynamic_cast<TList*>(GetOutputData(20));
  if (!fOutputpostvTracTrue) printf("ERROR:fOutputpostvTracTrue not available\n");
  
  fOutputnegtvTracTrue = dynamic_cast<TList*>(GetOutputData(21));
  if (!fOutputnegtvTracTrue) printf("ERROR:fOutputnegtvTracTrue not available\n");
  
  fOutputpullAllpointTrue = dynamic_cast<TList*>(GetOutputData(22));
  if (!fOutputpullAllpointTrue) printf("ERROR:fOutputpullAllpointTrue not available\n");
  
  fOutputphiAllpointSkip = dynamic_cast<TList*>(GetOutputData(23));
  if (!fOutputphiAllpointSkip) printf("ERROR:fOutputphiAllpointSkip not available\n");
  
  fOutputphiPostvtracSkip = dynamic_cast<TList*>(GetOutputData(24));
  if (!fOutputphiPostvtracSkip) printf("ERROR:fOutputphiPostvtracSkip not available\n");
  
  fOutputphiNegtvtracSkip = dynamic_cast<TList*>(GetOutputData(25));
  if (!fOutputphiNegtvtracSkip) printf("ERROR:fOutputphiNegtvtracSkip not available\n");
  
  fOutputclusterTypeSPD01Skip = dynamic_cast<TList*>(GetOutputData(26));
  if (!fOutputclusterTypeSPD01Skip) printf("ERROR:fOutputclusterTypeSPD01Skip not available\n");
  
  fOutputclusterTypeSPD02Skip = dynamic_cast<TList*>(GetOutputData(27));
  if (!fOutputclusterTypeSPD02Skip) printf("ERROR:fOutputclusterTypeSPD02Skip not available\n");
  
  fOutputclusterTypeSPD03Skip = dynamic_cast<TList*>(GetOutputData(28));
  if (!fOutputclusterTypeSPD03Skip) printf("ERROR:fOutputclusterTypeSPD03Skip not available\n");
  
  fOutputclusterTypeSPD11Skip = dynamic_cast<TList*>(GetOutputData(29));
  if (!fOutputclusterTypeSPD11Skip) printf("ERROR:fOutputclusterTypeSPD11Skip not available\n");
  
  fOutputclusterTypeSPD12Skip = dynamic_cast<TList*>(GetOutputData(30));
  if (!fOutputclusterTypeSPD12Skip) printf("ERROR:fOutputclusterTypeSPD12Skip not available\n");
  
  fOutputclusterTypeSPD13Skip = dynamic_cast<TList*>(GetOutputData(31));
  if (!fOutputclusterTypeSPD13Skip) printf("ERROR:fOutputclusterTypeSPD13Skip not available\n");
  
  fOutputPt = dynamic_cast<TList*>(GetOutputData(32));
  if (!fOutputPt) printf("ERROR: fOutputPt not available\n");
  
  fNentries = dynamic_cast<TH1F*>(GetOutputData(33));
  if (!fNentries) printf("ERROR: fNentries not available\n");
  
  fEstimVtx = dynamic_cast<TH1F*>(GetOutputData(34));
  if (!fEstimVtx) printf("ERROR: fEstimVtx not available\n");
  return;
}


Int_t AliAnalysisTaskSEImpParRes::ClusterTypeOnITSLayer(AliESDtrack *track,Int_t layer) const {
//
// Returns cluster type on ITS layer. Returns -1 if no cluster on this layer
//
   Int_t ctype=-1;

  if(layer<0 || layer>5) return ctype;
  if(!track->HasPointOnITSLayer(layer)) return ctype;

  const AliTrackPointArray *array = track->GetTrackPointArray();
    if(!array) {
      printf("No tracks points avaialble: check ESDfriends\n");
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

