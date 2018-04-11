/*
 * AliFemtoDreamv0Hist.cxx
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamv0Hist.h"
#include "TMath.h"
ClassImp(AliFemtoDreamv0Hist)
AliFemtoDreamv0Hist::AliFemtoDreamv0Hist()
:fMinimalBooking(false)
,fHistList(0)
,fConfig(0)
,fCutCounter(0)
,fInvMassBefKaonRej(0)
,fInvMassKaon(0)
,fInvMassBefSelection(0)
,fInvMassPt(0)
,fCPAPtBins(0)
,fInvMassPerRunNumber()
{
  for (int i=0;i<2;++i) {
    fv0CutQA[i]=nullptr;
    fOnFly[i]=nullptr;
    fpTDist[i]=nullptr;
    fetaDist[i]=nullptr;
    fDecayVtxv0X[i]=nullptr;
    fDecayVtxv0Y[i]=nullptr;
    fDecayVtxv0Z[i]=nullptr;
    fTransRadius[i]=nullptr;
    fDCAPosDaugToPrimVtx[i]=nullptr;
    fDCANegDaugToPrimVtx[i]=nullptr;
    fDCADaugToVtx[i]=nullptr;
    fCPA[i]=nullptr;
    fInvMass[i]=nullptr;
  }
}

AliFemtoDreamv0Hist::AliFemtoDreamv0Hist(
    int MassNBins,float MassMin,float MassMax,bool CPAPlots,bool perRunnumber, int iRunMin=0, int iRunMax=0)
:fMinimalBooking(false)
{
  TString sName[2]={"before","after"};

  fHistList=new TList();
  fHistList->SetName("v0Cuts");
  fHistList->SetOwner();

  fConfig = new TProfile("TrackCutConfig", "Track Cut Config", 20, 0, 20);
  fConfig->SetStats(0);
  fConfig->GetXaxis()->SetBinLabel(1, "OnFly Status");
  fConfig->GetXaxis()->SetBinLabel(2, "Charge");
  fConfig->GetXaxis()->SetBinLabel(3, "p_{T, Min}");
  fConfig->GetXaxis()->SetBinLabel(4, "p_{T, Max}");
  fConfig->GetXaxis()->SetBinLabel(5, "K0 Rejection");
  fConfig->GetXaxis()->SetBinLabel(6, "Max Decay Vertex XYZ");
  fConfig->GetXaxis()->SetBinLabel(7, "Min Transverse Radius");
  fConfig->GetXaxis()->SetBinLabel(8, "Max Transverse Radius");
  fConfig->GetXaxis()->SetBinLabel(9, "Min DCA to PV of Daug");
  fConfig->GetXaxis()->SetBinLabel(10,"Max Distance Daug to Vertex");
  fConfig->GetXaxis()->SetBinLabel(11,"Inv Mass Cut");
  fConfig->GetXaxis()->SetBinLabel(12,"Min Cos Pointing Angle");

  fHistList->Add(fConfig);

  fCutCounter = new TH1F("CutCounter", "Cut Counter", 20, 0, 20);
  fCutCounter->GetXaxis()->SetBinLabel(1, "Input");
  fCutCounter->GetXaxis()->SetBinLabel(2, "Has Daug");
  fCutCounter->GetXaxis()->SetBinLabel(3, "Daugters pass Cuts");
  fCutCounter->GetXaxis()->SetBinLabel(4, "On Fly Status");
  fCutCounter->GetXaxis()->SetBinLabel(5, "Charge");
  fCutCounter->GetXaxis()->SetBinLabel(6, "p_{T}");
  fCutCounter->GetXaxis()->SetBinLabel(7, "v0 DCA PV");
  fCutCounter->GetXaxis()->SetBinLabel(8, "Transverse Radius");
  fCutCounter->GetXaxis()->SetBinLabel(9, "Daug PV");
  fCutCounter->GetXaxis()->SetBinLabel(10, "Daug Vtx");
  fCutCounter->GetXaxis()->SetBinLabel(11, "K0 Rejection");
  fCutCounter->GetXaxis()->SetBinLabel(12, "Inv Mass Cut");
  fCutCounter->GetXaxis()->SetBinLabel(13, "Cos Pointing Angle");
  fCutCounter->GetXaxis()->SetBinLabel(14, "D1&D2 right");
  fCutCounter->GetXaxis()->SetBinLabel(15, "D1&D2 pass cuts");
  fCutCounter->GetXaxis()->SetBinLabel(16, "D1&D2 wrong");
  fCutCounter->GetXaxis()->SetBinLabel(16, "D1&D2 pass cuts");

  fHistList->Add(fCutCounter);

  fInvMassBefKaonRej = new TH1F("InvMassBefK0Rej","InvMassBefK0Rej",
                                MassNBins,MassMin,MassMax);
  fInvMassBefKaonRej->Sumw2();
  fInvMassBefKaonRej->GetXaxis()->SetName("m_{Pair}");
  fHistList->Add(fInvMassBefKaonRej);

  fInvMassKaon = new TH1F("InvMassKaon","InvMassKaon",400,0.4,0.6);
  fInvMassKaon->Sumw2();
  fInvMassKaon->GetXaxis()->SetName("m_{Pair}");
  fHistList->Add(fInvMassKaon);

  fInvMassBefSelection = new TH1F("InvMasswithCuts","InvMasswithCuts",MassNBins,MassMin,MassMax);
  fInvMassBefSelection->Sumw2();
  fInvMassBefSelection->GetXaxis()->SetName("m_{Pair}");
  fHistList->Add(fInvMassBefSelection);

  fInvMassPt=new TH2F("InvMassPt","Invariant Mass in Pt Bins",
                      8,0.3,4.3,MassNBins,MassMin,MassMax);
  fInvMassPt->Sumw2();
  fInvMassPt->GetXaxis()->SetTitle("P_{T}");
  fInvMassPt->GetYaxis()->SetTitle("m_{Pair}");
  fHistList->Add(fInvMassPt);

  for(int i=0;i<2;++i){
    fv0CutQA[i] = new TList();
    fv0CutQA[i]->SetName(sName[i].Data());
    fv0CutQA[i]->SetOwner();
    fHistList->Add(fv0CutQA[i]);

    TString OnFlyName = Form("OnFly_%s",sName[i].Data());
    fOnFly[i] = new TH1F(OnFlyName.Data(),OnFlyName.Data(), 2, 0, 2.);
    fOnFly[i]->Sumw2();
    fOnFly[i]->GetXaxis()->SetBinLabel(1,"Online");
    fOnFly[i]->GetXaxis()->SetBinLabel(2,"Offline");
    fv0CutQA[i]->Add(fOnFly[i]);

    TString ptname = Form("pTDist_%s",sName[i].Data());
    fpTDist[i] = new TH1F(ptname.Data(),ptname.Data(),100,0,10.);
    fpTDist[i]->Sumw2();
    fpTDist[i]->GetXaxis()->SetTitle("P_{T}");
    fv0CutQA[i]->Add(fpTDist[i]);

    TString etaname = Form("EtaDist_%s",sName[i].Data());
    fetaDist[i] = new TH1F(etaname.Data(),etaname.Data(),200,-10.,10.);
    fetaDist[i]->Sumw2();
    fetaDist[i]->GetXaxis()->SetTitle("#eta");
    fv0CutQA[i]->Add(fetaDist[i]);

    TString phiname = Form("PhiDist_%s",sName[i].Data());
    fPhiDist[i] = new TH1F(phiname.Data(),phiname.Data(),100,0.,2*TMath::Pi());
    fPhiDist[i]->Sumw2();
    fPhiDist[i]->GetXaxis()->SetTitle("#phi");
    fv0CutQA[i]->Add(fPhiDist[i]);

    TString decayVtxXname = Form("DecayVtxXPV_%s",sName[i].Data());
    fDecayVtxv0X[i] = new TH1F(decayVtxXname.Data(),decayVtxXname.Data(),
                               400,0.,200);
    fDecayVtxv0X[i]->Sumw2();
    fDecayVtxv0X[i]->GetXaxis()->SetTitle("Decay Vtx To PV X");
    fv0CutQA[i]->Add(fDecayVtxv0X[i]);

    TString decayVtxYname = Form("DecayVtxYPV_%s",sName[i].Data());
    fDecayVtxv0Y[i] = new TH1F(decayVtxYname.Data(),decayVtxYname.Data(),
                               400,0.,200);
    fDecayVtxv0Y[i]->Sumw2();
    fDecayVtxv0Y[i]->GetXaxis()->SetTitle("Decay Vtx To PV Y");
    fv0CutQA[i]->Add(fDecayVtxv0Y[i]);

    TString decayVtxZname = Form("DecayVtxZPV_%s",sName[i].Data());
    fDecayVtxv0Z[i] = new TH1F(decayVtxZname.Data(),
                               decayVtxZname.Data(),400,0.,200);
    fDecayVtxv0Z[i]->Sumw2();
    fDecayVtxv0Z[i]->GetXaxis()->SetTitle("Decay Vtx To PV Z");
    fv0CutQA[i]->Add(fDecayVtxv0Z[i]);

    TString transverseRadname = Form("TransverseRadius_%s",sName[i].Data());
    fTransRadius[i] = new TH1F(transverseRadname.Data(),
                               transverseRadname.Data(),750,0,150);
    fTransRadius[i]->Sumw2();
    fTransRadius[i]->GetXaxis()->SetTitle("Transverse Radius");
    fv0CutQA[i]->Add(fTransRadius[i]);

    TString DCADauPVPname = Form("DCADauPToPV_%s",sName[i].Data());
    fDCAPosDaugToPrimVtx[i] = new TH1F(DCADauPVPname.Data(),
                                       DCADauPVPname.Data(),500,0,100);
    fDCAPosDaugToPrimVtx[i]->Sumw2();
    fDCAPosDaugToPrimVtx[i]->GetXaxis()->SetTitle("DCADaugther P to PV");
    fv0CutQA[i]->Add(fDCAPosDaugToPrimVtx[i]);

    TString DCADauPVNname = Form("DCADauNToPV_%s", sName[i].Data());
    fDCANegDaugToPrimVtx[i] = new TH1F(DCADauPVNname.Data(),
                                       DCADauPVNname.Data(),500,0,100);
    fDCANegDaugToPrimVtx[i]->Sumw2();
    fDCANegDaugToPrimVtx[i]->GetXaxis()->SetTitle("DCADaugther N to PV");
    fv0CutQA[i]->Add(fDCANegDaugToPrimVtx[i]);

    TString DCADaugVtxname = Form("DCADauToVtx_%s", sName[i].Data());
    fDCADaugToVtx[i] = new TH1F(DCADaugVtxname.Data(),DCADaugVtxname.Data(),
                                100,0,10);
    fDCADaugToVtx[i]->Sumw2();
    fDCADaugToVtx[i]->GetXaxis()->SetTitle("DCA Daug to Vtx");
    fv0CutQA[i]->Add(fDCADaugToVtx[i]);

    TString cosPointName = Form("PointingAngle_%s", sName[i].Data());
    fCPA[i] = new TH1F(cosPointName.Data(),cosPointName.Data(),500,0.8,1.001);
    fCPA[i]->Sumw2();
    fCPA[i]->GetXaxis()->SetTitle("Cos Pointing Angle");
    fv0CutQA[i]->Add(fCPA[i]);

    TString invMassName = Form("InvariantMass_%s", sName[i].Data());
    fInvMass[i] = new TH1F(invMassName.Data(),invMassName.Data(),
                           MassNBins,MassMin,MassMax);
    fInvMass[i]->Sumw2();
    fInvMass[i]->GetXaxis()->SetTitle("m_{Pair}");
    fv0CutQA[i]->Add(fInvMass[i]);
  }

  if (CPAPlots) {
      fCPAPtBins=new TH2F("CPAPtBinsTot","CPAPtBinsTot",
                          8,0.3,4.3,1000,0.90,1.);
      fCPAPtBins->Sumw2();
      fCPAPtBins->GetXaxis()->SetTitle("P_{T}");
      fCPAPtBins->GetYaxis()->SetTitle("CPA");
      fHistList->Add(fCPAPtBins);
  } else {
      fCPAPtBins=nullptr;
  }
  if (perRunnumber) {
    int nBins=iRunMax-iRunMin;
    TString InvMassRunNumbName="InvMassPerRunnumber";
    fInvMassPerRunNumber= new TH2F(InvMassRunNumbName.Data(),InvMassRunNumbName.Data(),
                                   nBins,iRunMin,iRunMax,MassNBins,MassMin,MassMax);
    fHistList->Add(fInvMassPerRunNumber);
  } else {
    fInvMassPerRunNumber=nullptr;
  }
}

AliFemtoDreamv0Hist::AliFemtoDreamv0Hist(TString MinimalBooking,int MassNBins,float MassMin,float MassMax)
:fMinimalBooking(true)
,fConfig(0)
,fCutCounter(0)
,fInvMassBefKaonRej(0)
,fInvMassKaon(0)
,fInvMassBefSelection(0)
,fInvMassPt(0)
,fCPAPtBins(0)
,fInvMassPerRunNumber()
{
  for (int i=0;i<2;++i) {
    fv0CutQA[i]=nullptr;
    fOnFly[i]=nullptr;
    fpTDist[i]=nullptr;
    fetaDist[i]=nullptr;
    fDecayVtxv0X[i]=nullptr;
    fDecayVtxv0Y[i]=nullptr;
    fDecayVtxv0Z[i]=nullptr;
    fTransRadius[i]=nullptr;
    fDCAPosDaugToPrimVtx[i]=nullptr;
    fDCANegDaugToPrimVtx[i]=nullptr;
    fDCADaugToVtx[i]=nullptr;
    fCPA[i]=nullptr;
    fInvMass[i]=nullptr;
  }
  fHistList=new TList();
  fHistList->SetOwner();
  fHistList->SetName(MinimalBooking.Data());

  fInvMassPt=new TH2F("InvMassPt","Invariant Mass in Pt Bins",
                      8,0.3,4.3,MassNBins,MassMin,MassMax);
  fInvMassPt->Sumw2();
  fInvMassPt->GetXaxis()->SetTitle("P_{T}");
  fInvMassPt->GetYaxis()->SetTitle("m_{Pair}");
  fHistList->Add(fInvMassPt);
}

AliFemtoDreamv0Hist::~AliFemtoDreamv0Hist() {
  delete fHistList;
}

