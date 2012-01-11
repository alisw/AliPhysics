
/**************************************************************************
 *  Authors : Iouri Belikov, Antonin Maire
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskGlobalQA class
// This task is for running the GlobalQA over already existing ESDs
//          Origin:  I.Belikov, Iouri.Belikov@cern.ch, June 2009
//-----------------------------------------------------------------

#include <TPDGCode.h>

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDInputHandler.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h" 


#include "AliAnalysisTaskGlobalQA.h"

//
// Run the GlobalQA analysis over already existing ESDs
// Origin: Iouri.Belikov@cern.ch
//

ClassImp(AliAnalysisTaskGlobalQA)

//________________________________________________________________________
AliAnalysisTaskGlobalQA::AliAnalysisTaskGlobalQA() : 
 AliAnalysisTaskSE("GlobalQA"), 
 fArrayQA(0)
{
  // Default Constructor

  DefineOutput(1, TObjArray::Class());
}

//________________________________________________________________________
void AliAnalysisTaskGlobalQA::UserCreateOutputObjects()
{
  // Create the histograms
  // Called once

  fArrayQA = new TObjArray(kLast);


  {// Event related QA
    const Char_t *name[]={
      "hGlobalPrimaryVertex"
    };
    const Char_t *title[]={
      "Z-distribution of the primary vertex"
    };
    Add2ESDsList(new TH1F(name[0],title[0],100,-20.,20.),kEvt0);
  }
 
  {// Cluster related QA
    const Char_t *name[]={
      "hGlobalFractionAssignedClustersITS",
      "hGlobalFractionAssignedClustersTPC",
      "hGlobalFractionAssignedClustersTRD",
      "hGlobalClustersPerITSModule"
     };
  const Char_t *title[]={
    "Fraction of the assigned clusters in ITS",
    "Fraction of the assigned clusters in TPC",
    "Fraction of the assigned clusters in TRD",
    "Number of clusters per an ITS module"
  };
  Add2ESDsList(new TH1F(name[0],title[0],100,0.,2.),kClr0);
  Add2ESDsList(new TH1F(name[1],title[1],100,0.,2.),kClr1);
  Add2ESDsList(new TH1F(name[2],title[2],100,0.,2.),kClr2);
  Add2ESDsList(new TH1F(name[3],title[3],2201,-0.5,2200.5),kClr3);
  }

  {// Track related QA
    const Char_t *name[]={
      "hGlobalTrackAzimuthe",                               // kTrk0
      "hGlobalTrackEta",                                    // kTrk1
      "hGlobalTPCTrackpT",                                  // kTrk2
      "hGlobalTPCITSMatchedpT",                             // kTrk3
      "hGlobalTPCTOFMatchedpT",                             // kTrk4
      "hGlobalTPCITSMatchingProbability",                   // kTrk5
      "hGlobalTPCTOFMatchingProbability",                   // kTrk6
      "hGlobalTPCsideAposDCA",                              // kTrk7
      "hGlobalTPCsideAnegDCA",                              // kTrk8
      "hGlobalTPCsideCposDCA",                              // kTrk9
      "hGlobalTPCsideCnegDCA"                               // kTrk10
  };
  const Char_t *title[]={
    "Track azimuthal distribution (rad)",                   // kTrk0
    "Track pseudo-rapidity distribution",                   // kTrk1
    "TPC: track momentum distribution (GeV)",               // kTrk2
    "TPC-ITS matched: track momentum distribution (GeV)",   // kTrk3
    "TPC-TOF matched: track momentum distribution (GeV)",   // kTrk4
    "TPC-ITS track-matching probability",                   // kTrk5
    "TPC-TOF track-matching probability",                   // kTrk6
    "TPC side A: DCA for the positive tracks (mm)",         // kTrk7
    "TPC side A: DCA for the negative tracks (mm)",         // kTrk8
    "TPC side C: DCA for the positive tracks (mm)",         // kTrk9
    "TPC side C: DCA for the negative tracks (mm)"          // kTrk10
  };
  Add2ESDsList(new TH1F(name[0],title[0],100, 0.,TMath::TwoPi()),kTrk0);
  Add2ESDsList(new TH1F(name[1],title[1],100,-2.00,2.00),kTrk1);
  Add2ESDsList(new TH1F(name[2],title[2],50,  0.20,5.00),kTrk2);
  Add2ESDsList(new TH1F(name[3],title[3],50,  0.20,5.00),kTrk3);
  Add2ESDsList(new TH1F(name[4],title[4],50,  0.20,5.00),kTrk4);
  Add2ESDsList(new TH1F(name[5],title[5],50,  0.20,5.00),kTrk5);
  Add2ESDsList(new TH1F(name[6],title[6],50,  0.20,5.00),kTrk6);
  Add2ESDsList(new TH1F(name[7],title[7],50, -25.0,25.0),kTrk7);
  Add2ESDsList(new TH1F(name[8],title[8],50, -25.0,25.0),kTrk8);
  Add2ESDsList(new TH1F(name[9],title[9],50, -25.0,25.0),kTrk9);
  Add2ESDsList(new TH1F(name[10],title[10],50, -25.0,25.0),kTrk10);
  }

  {// V0 related QA
    const Char_t *name[]={
      "hGlobalPromptK0sMass",
      "hGlobalOfflineK0sMass",
      "hGlobalPromptLambda0Lambda0BarMass",
      "hGlobalOfflineLambda0Lambda0BarMass"
    };
  const Char_t *title[]={
    "On-the-fly K0s mass (GeV)",
    "Offline K0s mass (GeV)",
    "On-the-fly Lambda0 + Lambda0Bar mass (GeV)",
    "Offline Lambda0 + Lambda0Bar mass (GeV)"
  };
  Add2ESDsList(new TH1F(name[0],title[0],50,  0.4477,0.5477),kK0on);
  Add2ESDsList(new TH1F(name[1],title[1],50,  0.4477,0.5477),kK0off);
  Add2ESDsList(new TH1F(name[2],title[2],50,  1.0657,1.1657),kL0on);
  Add2ESDsList(new TH1F(name[3],title[3],50,  1.0657,1.1657),kL0off);
  }

  {// PID related QA
  const Char_t *name[]={
    "hGlobalITSdEdx",
    "hGlobalTPCdEdx",
    "hGlobalTOFTrackingvsMeasured",
    "hGlobalTPCdEdxvsMomentum"
  };
  const Char_t *title[]={
    "ITS: dEdx (A.U.) for particles with momentum 0.4 - 0.5 (GeV)",
    "TPC: dEdx (A.U.) for particles with momentum 0.4 - 0.5 (GeV)",
    "TOF: tracking - measured (ps)",
    "TPC: dEdx (A.U.) vs momentum (GeV)"
  };
  Add2ESDsList(new TH1F(name[0],title[0],50,0.00,200.),kPid0);
  Add2ESDsList(new TH1F(name[1],title[1],50,0.00,100.),kPid1);
  Add2ESDsList(new TH1F(name[2],title[2],50,-3500.,3500.),kPid2);
  Add2ESDsList(new TH2F(name[3],title[3],1500,0.05,15.,700,0.,700.),kPid3);
  }

  {// Multiplicity related QA
    const Char_t *name[]={
      "hGlobalV0AvsITS",
      "hGlobalV0CvsITS"
    };
    const Char_t *title[]={
      "Multiplicity: V0A vs ITS",
      "Multiplicity: V0C vs ITS"
    };
    TH2F *h0=new TH2F(name[0],title[0],41,-0.5,40.5, 33,-0.5,32.5);
    Add2ESDsList(h0,kMlt0);
    TH2F *h1=new TH2F(name[1],title[1],41,-0.5,40.5, 33,-0.5,32.5);
    Add2ESDsList(h1,kMlt1);
  }

}

//________________________________________________________________________
void AliAnalysisTaskGlobalQA::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  const AliESDEvent *esd=(const AliESDEvent *)InputEvent();

  if (!esd) {
    Printf("ERROR: ESD is not available");
    return;
  }

  // Event related QA
  const AliESDVertex *vtx=esd->GetPrimaryVertex();
  if (!vtx->GetStatus()) return;

  Double_t xv=vtx->GetXv();
  Double_t yv=vtx->GetYv();
  Double_t zv=vtx->GetZv();
  GetESDsData(kEvt0)->Fill(zv);


  for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* track = esd->GetTrack(iTracks);
      if (!track) {
         Printf("ERROR: Could not receive track %d", iTracks);
         continue;
      }


    // Cluster related QA
    if (track->IsOn(AliESDtrack::kITSrefit)) {
      Int_t n=track->GetITSclusters(0);
      GetESDsData(kClr0)->Fill(Float_t(n)/6.); //6 is the number of ITS layers
    }

    for (Int_t i=0; i<6; i++) {
      Int_t idet, sts;
      Float_t xloc,zloc;
      if (!track->GetITSModuleIndexInfo(i,idet,sts,xloc,zloc)) continue;
      if (i>=2) idet+=240;
      if (i>=4) idet+=260;
      if ((sts==1)||(sts==2)||(sts==4)) GetESDsData(kClr3)->Fill(idet);
    }

    if (track->IsOn(AliESDtrack::kTPCrefit)) {
      Int_t n =track->GetTPCNcls();
      Int_t nf=track->GetTPCNclsF();      // number of crossed TPC pad rows
      if (nf>0) {
        Double_t val = n*1.0/nf; 
        GetESDsData(kClr1)->Fill(val); 
      }
    }

    if (track->IsOn(AliESDtrack::kTRDrefit)) {
      Int_t n=track->GetTRDclusters(0);
      GetESDsData(kClr2)->Fill(Float_t(n)/(6*24));//(6*24) is the number of TRD time bins
    }

    Double_t p=track->GetP();

    // Track related QA
    if (track->IsOn(AliESDtrack::kTPCrefit)) {
      Float_t dz[2]; 
      track->GetDZ(xv,yv,zv,esd->GetMagneticField(),dz); 
      if ((TMath::Abs(dz[0])<3.) && (TMath::Abs(dz[1])<3.)) { // beam pipe
        Double_t phi=track->Phi();
	GetESDsData(kTrk0)->Fill(phi);
	Double_t y=track->Eta();
	GetESDsData(kTrk1)->Fill(y);

        if (TMath::Abs(y)<0.9) {
	   GetESDsData(kTrk2)->Fill(p);
	   if (track->IsOn(AliESDtrack::kITSrefit)) GetESDsData(kTrk3)->Fill(p);
	  //if (track->IsOn(AliESDtrack::kTOFout)) GetESDsData(kTrk4)->Fill(p);
	   if (track->GetTOFsignal()>0) GetESDsData(kTrk4)->Fill(p);
	}
      }
    }
    const AliExternalTrackParam *tpcTrack=track->GetTPCInnerParam();
    const AliExternalTrackParam *innTrack=track->GetInnerParam();
    if (tpcTrack)
    if (innTrack) {
       Float_t dz[2];
       tpcTrack->GetDZ(xv,yv,zv,esd->GetMagneticField(),dz);
       dz[0]*=10.; // in mm
       if (innTrack->GetZ()  > 0)
       if (innTrack->GetTgl()> 0) { // TPC side A
          if (tpcTrack->GetSign() > 0) GetESDsData(kTrk7)->Fill(dz[0]);
          else                         GetESDsData(kTrk8)->Fill(dz[0]);
       }
       if (innTrack->GetZ()  < 0)
       if (innTrack->GetTgl()< 0) { // TPC side C
          if (tpcTrack->GetSign() > 0) GetESDsData(kTrk9)->Fill(dz[0]);
          else                         GetESDsData(kTrk10)->Fill(dz[0]);
       }
    }

    // PID related QA
    if ((p>0.4) && (p<0.5)) {
      if (track->IsOn(AliESDtrack::kITSpid)) {
	Double_t dedx=track->GetITSsignal();
        GetESDsData(kPid0)->Fill(dedx);
      }
      if (track->IsOn(AliESDtrack::kTPCpid)) {
	Double_t dedx=track->GetTPCsignal();
        GetESDsData(kPid1)->Fill(dedx);
      }
    }
    if (p>1.0) {
      if (track->IsOn(AliESDtrack::kITSrefit))
      if (track->IsOn(AliESDtrack::kTPCrefit))
      if (track->IsOn(AliESDtrack::kTOFout)) {
         Float_t dz[2];
         track->GetDZ(xv,yv,zv,esd->GetMagneticField(),dz);
         if (dz[1]<3.) {
            Double_t times[10];
            track->GetIntegratedTimes(times);
            Double_t tof=track->GetTOFsignal()/*-847055 -1771207*/;
            GetESDsData(kPid2)->Fill(times[2]-tof);
	 }
      }
    }
    const AliExternalTrackParam *par=track->GetInnerParam();
    if (par) {
      Double_t pp=par->GetP();
      Double_t dedx=track->GetTPCsignal();
      TH2F *h = dynamic_cast<TH2F*>(GetESDsData(kPid3));
      if (h) h->Fill(pp,dedx);
    }
  }

  // Multiplicity related QA
  AliESDVZERO     *mltV0 =esd->GetVZEROData();
  const AliMultiplicity *mltITS=esd->GetMultiplicity();
  if (mltV0)
    if (mltITS) {
       Short_t nv0a=mltV0->GetNbPMV0A();
       Short_t nv0c=mltV0->GetNbPMV0C();
       Int_t   nits=mltITS->GetNumberOfTracklets();
       TH2F *h0=dynamic_cast<TH2F*>(GetESDsData(kMlt0));
       if (h0) h0->Fill(nits,nv0a);
       TH2F *h1=dynamic_cast<TH2F*>(GetESDsData(kMlt1));
       if (h1) h1->Fill(nits,nv0c);
    }

  // V0 related QA
  Int_t nV0=esd->GetNumberOfV0s();
  for (Int_t i=0; i<nV0; i++) {
    Double_t mass;
    AliESDv0 v0(*esd->GetV0(i));

    Int_t nidx=TMath::Abs(v0.GetNindex());
    AliESDtrack *ntrack1=esd->GetTrack(nidx);
    if (!ntrack1->IsOn(AliESDtrack::kTPCrefit)) continue;

    Int_t pidx=TMath::Abs(v0.GetPindex());
    AliESDtrack *ptrack1=esd->GetTrack(pidx);
    if (!ptrack1->IsOn(AliESDtrack::kTPCrefit)) continue;

    v0.ChangeMassHypothesis(kK0Short);
    mass=v0.GetEffMass();
    if (v0.GetOnFlyStatus())
       GetESDsData(kK0on)->Fill(mass);
    else
       GetESDsData(kK0off)->Fill(mass);

    v0.ChangeMassHypothesis(kLambda0);
    mass=v0.GetEffMass();
    if (v0.GetOnFlyStatus())
       GetESDsData(kL0on)->Fill(mass);
    else
       GetESDsData(kL0off)->Fill(mass);

    v0.ChangeMassHypothesis(kLambda0Bar);
    mass=v0.GetEffMass();
    if (v0.GetOnFlyStatus())
       GetESDsData(kL0on)->Fill(mass);
    else
       GetESDsData(kL0off)->Fill(mass);
  }

  // Post output data.
  PostData(1, fArrayQA);
}      

//________________________________________________________________________
void AliAnalysisTaskGlobalQA::Terminate(Option_t *) 
{
  // Draw the results on the screen
  // Called once at the end of the query

  fArrayQA=(TObjArray*)GetOutputData(1);

  TH1 *tpc=GetESDsData(kTrk2); tpc->Sumw2();
  TH1 *its=GetESDsData(kTrk3); its->Sumw2();
  TH1 *tof=GetESDsData(kTrk4); tof->Sumw2();
  GetESDsData(kTrk5)->Divide(its,tpc,1,1.,"b");
  GetESDsData(kTrk6)->Divide(tof,tpc,1,1.,"b");

  TH1 *hTPCdEdxMIP = GetESDsData(kPid1);
  if (!hTPCdEdxMIP) {
    Printf("ERROR: hTPCdEdxMIP not available");
    return;
  }

  TH2 *hTPCdEdxVsP = dynamic_cast<TH2*>(GetESDsData(kPid3));
  if (!hTPCdEdxVsP) {
    Printf("ERROR: hTPCdEdxVsP not available");
    return;
  }

  TCanvas *c2=new TCanvas("c2","",320,32,530,590);

  TPad *p6=new TPad("p6","",0.,0.,1.,.5); p6->Draw(); p6->cd(); 
  p6->SetFillColor(42); p6->SetFrameFillColor(10);
  hTPCdEdxMIP->SetFillColor(2); hTPCdEdxMIP->SetFillStyle(3005);
  if (hTPCdEdxMIP->GetEntries()<333) 
      hTPCdEdxMIP->DrawCopy("E"); 
  else 
      hTPCdEdxMIP->Fit("gaus"); 
  c2->cd();

  TPad *p7=new TPad("p7","",0.,0.5,1.,1.); p7->Draw(); p7->cd(); p7->SetLogx();
  p7->SetFillColor(42); p7->SetFrameFillColor(10);
  hTPCdEdxVsP->DrawCopy();

}
