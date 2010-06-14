/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
// AliAnalysisTaskSE for D0 candidates invariant mass histogram
// and comparison with the MC truth and cut variables distributions.
//
// Authors: A.Dainese, andrea.dainese@lnl.infn.it
// Chiara Bianchin, chiara.bianchin@pd.infn.it (invariant mass)
// Carmelo Di Giglio, carmelo.digiglio@ba.infn.it (like sign)
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSED0Mass.h"


ClassImp(AliAnalysisTaskSED0Mass)


//________________________________________________________________________
AliAnalysisTaskSED0Mass::AliAnalysisTaskSED0Mass():
AliAnalysisTaskSE(),
fOutputMass(0),
fDistr(0),
fNentries(0), 
fChecks(0),
fCuts(0),
fArray(0),
fReadMC(0),
fCutOnDistr(0),
fNPtBins(1),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.)


{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::AliAnalysisTaskSED0Mass(const char *name,AliRDHFCutsD0toKpi* cuts):
AliAnalysisTaskSE(name),
fOutputMass(0), 
fDistr(0),
fNentries(0),
fChecks(0),
fCuts(0),
fArray(0),
fReadMC(0),
fCutOnDistr(0),
fNPtBins(1),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.)


{
  // Default constructor

  fNPtBins=cuts->GetNPtBins();
  fTotPosPairs=new Int_t[fNPtBins];
  fTotNegPairs=new Int_t[fNPtBins];
  for(Int_t i=0;i<fNPtBins;i++) {fTotPosPairs[i]=0; fTotNegPairs[i]=0;}
    
  fCuts=cuts;

  // Output slot #1 writes into a TList container (mass with cuts)
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TList container (distributions)
  DefineOutput(2,TList::Class());  //My private output
  // Output slot #3 writes into a TH1F container (number of events)
  DefineOutput(3,TH1F::Class());  //My private output
  // Output slot #4 writes into a TList container (quality check)
  DefineOutput(4,TList::Class());  //My private output
  // Output slot #5 writes into a TList container (cuts)
  DefineOutput(5,AliRDHFCutsD0toKpi::Class());  //My private output

}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::~AliAnalysisTaskSED0Mass()
{
   if (fOutputMass) {
    delete fOutputMass;
    fOutputMass = 0;
  }
  if (fDistr) {
    delete fDistr;
    fDistr = 0;
  }
  if (fChecks) {
    delete fChecks;
    fChecks = 0;
  }
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
 
 
}  

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Mass::Init() \n");

  
  AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*fCuts);
  const char* nameoutput=GetOutputSlot(5)->GetContainer()->GetName();
  copyfCuts->SetName(nameoutput);
  // Post the data
  PostData(5,copyfCuts);


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserCreateOutputObjects()
{

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("listMass");

  fDistr = new TList();
  fDistr->SetOwner();
  fDistr->SetName("distributionslist");

  fChecks = new TList();
  fChecks->SetOwner();
  fChecks->SetName("checkHistograms");

  TString nameMass=" ",nameSgn27=" ",nameSgn=" ", nameBkg=" ", nameRfl=" ",nameMassNocutsS =" ",nameMassNocutsB =" ", namedistr=" ";

  for(Int_t i=0;i<fCuts->GetNPtBins();i++){

    nameMass="histMass_";
    nameMass+=i;
    nameSgn27="histSgn27_";
    nameSgn27+=i;
    nameSgn="histSgn_";
    nameSgn+=i;
    nameBkg="histBkg_";
    nameBkg+=i;
    nameRfl="histRfl_";
    nameRfl+=i;
    nameMassNocutsS="hMassS_";
    nameMassNocutsS+=i;
    nameMassNocutsB="hMassB_";
    nameMassNocutsB+=i;

    //histograms of cut variable distributions

    //  pT
    namedistr="hptpiS_";
    namedistr+=i;
    TH1F *hptpiS = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptKS_";
    namedistr+=i;
    TH1F *hptKS = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptB_";
    namedistr+=i;
    TH1F *hptB = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

    //  pT no mass cut
    namedistr="hptpiSnoMcut_";
    namedistr+=i;
    TH1F *hptpiSnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptKSnoMcut_";
    namedistr+=i;
    TH1F *hptKSnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptB1prongnoMcut_";
    namedistr+=i;
    TH1F *hptB1pnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptB2prongsnoMcut_";
    namedistr+=i;
    TH1F *hptB2pnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

    
    //  dca
    namedistr="hdcaS_";
    namedistr+=i;
    TH1F *hdcaS = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);
    namedistr="hdcaB_";
    namedistr+=i;
    TH1F *hdcaB = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);

    //  costhetastar
    namedistr="hcosthetastarS_";
    namedistr+=i;
    TH1F *hcosthetastarS = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);
    namedistr="hcosthetastarB_";
    namedistr+=i;
    TH1F *hcosthetastarB = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);

    // impact parameter
    namedistr="hd0piS_";
    namedistr+=i;
    TH1F *hd0piS = new TH1F(namedistr.Data(), "Impact parameter distribution (pions);d0(#pi) [cm]",200,-0.1,0.1);

    namedistr="hd0KS_";
    namedistr+=i;
    TH1F *hd0KS = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons);d0(K) [cm]",200,-0.1,0.1);

    namedistr="hd0B_";
    namedistr+=i;
    TH1F *hd0B = new TH1F(namedistr.Data(), "Impact parameter distribution (both);d0 [cm]",200,-0.1,0.1);

    namedistr="hd0p0B_";
    namedistr+=i;
    TH1F *hd0p0B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong +);d0 [cm]",200,-0.1,0.1);

    namedistr="hd0p1B_";
    namedistr+=i;
    TH1F *hd0p1B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong -);d0 [cm]",200,-0.1,0.1);

    namedistr="hd0vpiS_";
    namedistr+=i;
    TH1F *hd0vpiS = new TH1F(namedistr.Data(), "Impact parameter distribution (pions)(vtx w/o these tracks);d0(#pi) [cm]",200,-0.1,0.1);

    namedistr="hd0vKS_";
    namedistr+=i;
    TH1F *hd0vKS = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons) (vtx w/o these tracks);d0(K) [cm]",200,-0.1,0.1);
    namedistr="hd0vB_";
    namedistr+=i;
    TH1F *hd0vB = new TH1F(namedistr.Data(), "Impact parameter distribution (vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

    namedistr="hd0vp0B_";
    namedistr+=i;
    TH1F *hd0vp0B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong + ** vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

    namedistr="hd0vp1B_";
    namedistr+=i;
    TH1F *hd0vp1B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong - ** vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

    namedistr="hd0d0S_";
    namedistr+=i;
    TH1F *hd0d0S = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);
    namedistr="hd0d0B_";
    namedistr+=i;
    TH1F *hd0d0B = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

     namedistr="hd0d0vS_";
    namedistr+=i;
    TH1F *hd0d0vS = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution (vtx w/o these tracks);d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);
    namedistr="hd0d0vB_";
    namedistr+=i;
    TH1F *hd0d0vB = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution (vtx w/o these tracks);d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

    //decay lenght
    namedistr="hdeclS_";
    namedistr+=i;
    TH1F *hdeclengthS = new TH1F(namedistr.Data(), "Decay Length distribution;Decay Length [cm]",200,0,0.6);

    namedistr="hdeclB_";
    namedistr+=i;
    TH1F *hdeclengthB = new TH1F(namedistr.Data(), "Decay Length distribution;Decay Length [cm]",200,0,0.6);

    namedistr="hnormdeclS_";
    namedistr+=i;
    TH1F *hnormdeclengthS = new TH1F(namedistr.Data(), "Normalized Decay Length distribution;Decay Length/Err ",200,0,10.);

    namedistr="hnormdeclB_";
    namedistr+=i;
    TH1F *hnormdeclengthB = new TH1F(namedistr.Data(), "Normalized Decay Length distribution;Decay Length/Err ",200,0,10.);

    namedistr="hdeclvS_";
    namedistr+=i;
    TH1F *hdeclengthvS = new TH1F(namedistr.Data(), "Decay Length distribution (vtx w/o tracks);Decay Length [cm]",200,0,0.6);

    namedistr="hdeclvB_";
    namedistr+=i;
    TH1F *hdeclengthvB = new TH1F(namedistr.Data(), "Decay Length distribution (vtx w/o tracks);Decay Length [cm]",200,0,0.6);

    namedistr="hnormdeclvS_";
    namedistr+=i;
    TH1F *hnormdeclengthvS = new TH1F(namedistr.Data(), "Normalized Decay Length distribution (vtx w/o tracks);Decay Length/Err ",200,0,10.);

    namedistr="hnormdeclvB_";
    namedistr+=i;
    TH1F *hnormdeclengthvB = new TH1F(namedistr.Data(), "Normalized Decay Length distribution (vtx w/o tracks);Decay Length/Err ",200,0,10.);

   //  costhetapoint
    namedistr="hcosthetapointS_";
    namedistr+=i;
    TH1F *hcosthetapointS = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);
    namedistr="hcosthetapointB_";
    namedistr+=i;
    TH1F *hcosthetapointB = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);

   //  costhetapoint vs d0 or d0d0
    namedistr="hcosthpointd0S_";
    namedistr+=i;
    TH2F *hcosthpointd0S= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0};cos#theta_{Point};d_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

    namedistr="hcosthpointd0B_";
    namedistr+=i;
    TH2F *hcosthpointd0B= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0};cos#theta_{Point};d_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

    namedistr="hcosthpointd0d0S_";
    namedistr+=i;
    TH2F *hcosthpointd0d0S= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);
    namedistr="hcosthpointd0d0B_";
    namedistr+=i;
    TH2F *hcosthpointd0d0B= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

    fDistr->Add(hptpiS);
    fDistr->Add(hptKS);
    fDistr->Add(hptB);

    fDistr->Add(hptpiSnoMcut);
    fDistr->Add(hptKSnoMcut);
    fDistr->Add(hptB1pnoMcut);
    fDistr->Add(hptB2pnoMcut);

    fDistr->Add(hdcaS);
    fDistr->Add(hdcaB);

    fDistr->Add(hd0piS);
    fDistr->Add(hd0KS);
    fDistr->Add(hd0B);
    fDistr->Add(hd0p0B);
    fDistr->Add(hd0p1B);

    fDistr->Add(hd0vpiS);
    fDistr->Add(hd0vKS);
    fDistr->Add(hd0vB);
    fDistr->Add(hd0vp0B);
    fDistr->Add(hd0vp1B);

    fDistr->Add(hd0d0S);
    fDistr->Add(hd0d0B);

    fDistr->Add(hd0d0vS);
    fDistr->Add(hd0d0vB);

    fDistr->Add(hcosthetastarS);
    fDistr->Add(hcosthetastarB);

    fDistr->Add(hcosthetapointS);
    fDistr->Add(hcosthetapointB);

    fDistr->Add(hdeclengthS);
    fDistr->Add(hdeclengthB);

    fDistr->Add(hnormdeclengthS);
    fDistr->Add(hnormdeclengthB);

    fDistr->Add(hdeclengthvS);
    fDistr->Add(hdeclengthvB);

    fDistr->Add(hnormdeclengthvS);
    fDistr->Add(hnormdeclengthvB);

    fDistr->Add(hcosthpointd0S);
    fDistr->Add(hcosthpointd0B);

    fDistr->Add(hcosthpointd0d0S);
    fDistr->Add(hcosthpointd0d0B);


    //histograms of invariant mass distributions

    TH1F* tmpMt = new TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpMl=(TH1F*)tmpMt->Clone();
    tmpMt->Sumw2();
    tmpMl->Sumw2();

    //to compare with AliAnalysisTaskCharmFraction
    TH1F* tmpS27t = new TH1F(nameSgn27.Data(),"D^{0} invariant mass in M(D^{0}) +/- 27 MeV - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpS27l=(TH1F*)tmpS27t->Clone();
    tmpS27t->Sumw2();
    tmpS27l->Sumw2();
 
    //distribution w/o cuts
    //   TH1F* tmpMS = new TH1F(nameMassNocutsS.Data(),"D^{0} invariant mass; M [GeV]; Entries",300,0.7,3.);
    TH1F* tmpMS = new TH1F(nameMassNocutsS.Data(),"D^{0} invariant mass; M [GeV]; Entries",300,1.56484,2.16484); //range (MD0-300MeV, mD0 + 300MeV)
    TH1F *tmpMB=(TH1F*)tmpMS->Clone();
    tmpMB->SetName(nameMassNocutsB.Data());
    tmpMS->Sumw2();
    tmpMB->Sumw2();

    //MC signal and background
    TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpSl=(TH1F*)tmpSt->Clone();
    tmpSt->Sumw2();
    tmpSl->Sumw2();

    TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpBl=(TH1F*)tmpBt->Clone();
    tmpBt->Sumw2();
    tmpBl->Sumw2();

    //Reflection: histo filled with D0Mass which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
    TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpRl=(TH1F*)tmpRt->Clone();
    tmpRt->Sumw2();
    tmpRl->Sumw2();
  //  printf("Created histograms %s\t%s\t%s\t%s\n",tmpM->GetName(),tmpS->GetName(),tmpB->GetName(),tmpR->GetName());

    fOutputMass->Add(tmpMt);
    fOutputMass->Add(tmpSt);
    fOutputMass->Add(tmpS27t);
    fOutputMass->Add(tmpBt);
    fOutputMass->Add(tmpRt);

    fDistr->Add(tmpMS);
    fDistr->Add(tmpMB);


  }

  //histograms for vertex checking and TOF checking
  TString checkname="hptGoodTr";

  TH1F* hptGoodTr=new TH1F(checkname.Data(),"Pt distribution of 'good' tracks;p_{t}[GeV];Number",200,0.,8.);
  hptGoodTr->SetTitleOffset(1.3,"Y");
  checkname="hdistrGoodTr";

  TH1F* hdistrGoodTr=new TH1F(checkname.Data(),"Distribution of number of good tracks per event;no.good-tracks/ev;Entries",31,0,31);
  hdistrGoodTr->SetTitleOffset(1.3,"Y");

  checkname="hTOFsig";
  TH1F* hTOFsig=new TH1F(checkname.Data(),"Distribution of TOF signal;TOF time [ps];Entries", 100, -2.e3,40.e3);

  checkname="hTOFtimeKaonHyptime";
  TH2F* hTOFtimeKaonHyptime=new TH2F(checkname.Data(),"TOFtime - timeHypothesisForKaon;p_{t}[GeV/c];TOFtime - timeHypothesisForKaon [ps]",200,0.,4.,1000,-20000.,20000.);
 
 checkname="hTOFtime";
  TH1F* hTOFtime=new TH1F(checkname.Data(),"Distribution of TOF time Kaon;TOF time(Kaon) [ps];Entries", 1000, 0.,50000.);
  

  fChecks->Add(hptGoodTr);
  fChecks->Add(hdistrGoodTr);
  fChecks->Add(hTOFsig);
  fChecks->Add(hTOFtimeKaonHyptime);
  fChecks->Add(hTOFtime);

  const char* nameoutput=GetOutputSlot(3)->GetContainer()->GetName();

  fNentries=new TH1F(nameoutput, "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of D0 selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 8,0.,8.);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  fNentries->GetXaxis()->SetBinLabel(3,"nD0Selected");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtx");
  fNentries->GetXaxis()->SetBinLabel(5,"nEventsGoodVtx+>2tracks");
  fNentries->GetXaxis()->SetBinLabel(6,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(7,"no daughter");
  fNentries->GetXaxis()->SetBinLabel(8,"nCandSel(Tr)");

  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

 
  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fDistr);
  PostData(3,fNentries);
  PostData(4,fChecks);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  //cout<<"I'm in UserExec"<<endl;


    //cuts order
    //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
  

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  TString bname;
  if(fArray==0){ //D0 candidates
    // load D0->Kpi candidates
    //cout<<"D0 candidates"<<endl;
    bname="D0toKpi";
  } else { //LikeSign candidates
    // load like sign candidates
    //cout<<"LS candidates"<<endl;
    bname="LikeSign2Prong";
  }

  TClonesArray *inputArray=0;
 
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
    }
  } else {
    inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  }


  if(!inputArray) {
    printf("AliAnalysisTaskSED0Mass::UserExec: input branch not found!\n");
    return;
  }
  
  
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSED0Mass::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSED0Mass::UserExec: MC header branch not found!\n");
      return;
    }
  }
  
  //printf("VERTEX Z %f %f\n",vtx1->GetZ(),mcHeader->GetVtxZ());
  
  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
    

  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  Bool_t isGoodVtx=kFALSE;

   //vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fNentries->Fill(3);
  }

  //cout<<"Start checks"<<endl;
  Int_t ntracks=0,isGoodTrack=0;

  if(aod) ntracks=aod->GetNTracks();

  //cout<<"ntracks = "<<ntracks<<endl;
 
  //loop on tracks in the event
  for (Int_t k=0;k<ntracks;k++){
    AliAODTrack* track=aod->GetTrack(k);

    //check TOF

    if(!(track->GetStatus()&AliESDtrack::kTPCrefit &&
	 track->GetStatus()&AliESDtrack::kITSrefit && 
	 track->GetTPCNcls() >=70 &&
	 track->GetStatus()&AliESDtrack::kTOFpid && 
	 track->GetStatus()&AliESDtrack::kTOFout &&
	 track->GetStatus()&AliESDtrack::kTIME)) continue;
    AliAODPid *pid = track->GetDetPid();
    if(!pid)  {if (fDebug>1)cout<<"No AliAODPid found"<<endl; continue;}

    Double_t times[5];
    pid->GetIntegratedTimes(times);

    ((TH1F*)fChecks->FindObject("hTOFtime"))->Fill(times[3]);
    ((TH1F*)fChecks->FindObject("hTOFtimeKaonHyptime"))->Fill(track->P(),pid->GetTOFsignal()-times[3]); //3 is kaon

    ((TH1F*)fChecks->FindObject("hTOFsig"))->Fill(pid->GetTOFsignal());
    if (pid->GetTOFsignal()< 0) ((TH1F*)fChecks->FindObject("hTOFsig"))->Fill(-1);


    //check clusters of the tracks
    Int_t nclsTot=0,nclsSPD=0;
    
    for(Int_t l=0;l<6;l++) {
      if(TESTBIT(track->GetITSClusterMap(),l)) {
	nclsTot++; if(l<2) nclsSPD++;
      }
    }
    
    if (track->Pt()>0.3 &&
	track->GetStatus()&AliESDtrack::kTPCrefit &&
	track->GetStatus()&AliESDtrack::kITSrefit &&
	nclsTot>3 &&
	nclsSPD>0) {//fill hist good tracks
     
      ((TH1F*)fChecks->FindObject("hptGoodTr"))->Fill(track->Pt());
      isGoodTrack++;
    }
    //cout<<"isGoodTrack = "<<isGoodTrack<<endl;
    ((TH1F*)fChecks->FindObject("hdistrGoodTr"))->Fill(isGoodTrack);
  }
  //number of events with good vertex and at least 2 good tracks
  if (isGoodTrack>=2 && isGoodVtx) fNentries->Fill(4);


  // loop over candidates
  Int_t nInD0toKpi = inputArray->GetEntriesFast();
  if(fDebug>2) printf("Number of D0->Kpi: %d\n",nInD0toKpi);

  for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
    //Int_t nPosPairs=0, nNegPairs=0;
    //cout<<"inside the loop"<<endl;
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iD0toKpi);

    //check daughters
    if(!(d->GetDaughter(0) || d->GetDaughter(1))) {
      AliDebug(1,"at least one daughter not found!");
      fNentries->Fill(6);
      return;
    }

    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
  
    
    //check reco daughter in acceptance
    /*
    Double_t eta0=d->EtaProng(0);
    Double_t eta1=d->EtaProng(1);
    */
    if ( fCuts->IsInFiducialAcceptance(d->Pt(),d->Y(421)) ) {
      //if( TMath::Abs(eta0)<0.9 && TMath::Abs(eta1)<0.9 ){
       //apply cuts on tracks
      Int_t isSelected = fCuts->IsSelected(d,AliRDHFCuts::kTracks);
      if(((AliAODTrack*)d->GetDaughter(0))->GetTPCNcls() < 70 || ((AliAODTrack*)d->GetDaughter(1))->GetTPCNcls() < 70) isSelected=kFALSE;
      if (!isSelected) return;
      fNentries->Fill(7);       
      if(fDebug>2) cout<<"tracks selected"<<endl;

      FillVarHists(aod,d,mcArray,fCuts,fDistr);
      FillMassHists(d,mcArray,fCuts,fOutputMass);
    }
  
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  } //end for prongs


  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fDistr);
  PostData(3,fNentries);
  PostData(4,fChecks);

  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0Mass::FillVarHists(AliAODEvent* aod,AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout){
  //
  // function used in UserExec to fill variable histograms:
  //

  //add distr here
  UInt_t pdgs[2];
    
  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  pdgs[0]=211;
  pdgs[1]=321;
  Double_t minvD0 = part->InvMassD0();
  pdgs[1]=211;
  pdgs[0]=321;
  Double_t minvD0bar = part->InvMassD0bar();
  //cout<<"inside mass cut"<<endl;

  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t lab=-9999;
  if(fReadMC) lab=part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
  //Double_t pt = d->Pt(); //mother pt
  Bool_t isSelected=kFALSE;


  if(fCutOnDistr){
    isSelected = cuts->IsSelected(part,AliRDHFCuts::kCandidate);
    if (!isSelected){
      //cout<<"Not Selected"<<endl;
      return;
    }
  }

  Double_t invmasscut=0.03;

  TString fillthispi="",fillthisK="",fillthis="";

  Int_t ptbin=cuts->PtBin(part->Pt());
  if(ptbin==-1) {fNentries->Fill(5); return;} //out of bounds

  //recalculate vertex
  AliAODVertex *vtxProp=0x0;
  vtxProp=GetPrimaryVtxSkipped(aod,part);
  Double_t dz1[2],dz2[2],covar1[3],covar2[3];//,d0xd0proper,errd0xd0proper;
  dz1[0]=-99; dz2[0]=-99;
  if(vtxProp) {
    part->SetOwnPrimaryVtx(vtxProp);
    //Bool_t unsetvtx=kTRUE;
    //Calculate d0 for daughter tracks with Vtx Skipped
    AliESDtrack *esdtrack1=new AliESDtrack((AliVTrack*)part->GetDaughter(0));
    AliESDtrack *esdtrack2=new AliESDtrack((AliVTrack*)part->GetDaughter(1));
    esdtrack1->PropagateToDCA(vtxProp,aod->GetMagneticField(),1.,dz1,covar1);
    esdtrack2->PropagateToDCA(vtxProp,aod->GetMagneticField(),1.,dz2,covar2);
    delete esdtrack1;
    delete esdtrack2;
  }

  Double_t d0[2]={dz1[0],dz2[0]};

  if(!fCutOnDistr || (fCutOnDistr && isSelected)){ //if no cuts or cuts passed 
    //printf("\nif no cuts or cuts passed\n");
    if(lab>=0 && fReadMC){ //signal

      //check pdg of the prongs
      AliAODTrack *prong0=(AliAODTrack*)part->GetDaughter(0);
      AliAODTrack *prong1=(AliAODTrack*)part->GetDaughter(1);
      if(!prong0 || !prong1) {
	return;
      }
      Int_t labprong[2];
      labprong[0]=prong0->GetLabel();
      labprong[1]=prong1->GetLabel();
      AliAODMCParticle *mcprong=0;
      Int_t pdgProng[2]={0,0};
      
      for (Int_t iprong=0;iprong<2;iprong++){
	if(labprong[iprong]>=0)  mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
	pdgProng[iprong]=mcprong->GetPdgCode();
      }

      //no mass cut ditributions: ptbis
	
      fillthispi="hptpiSnoMcut_";
      fillthispi+=ptbin;

      fillthisK="hptKSnoMcut_";
      fillthisK+=ptbin;

      if (TMath::Abs(pdgProng[0]) == 211 && TMath::Abs(pdgProng[1]) == 321){
	((TH1F*)listout->FindObject(fillthispi))->Fill(part->PtProng(0));
	((TH1F*)listout->FindObject(fillthisK))->Fill(part->PtProng(1));
      }else {
	if (TMath::Abs(pdgProng[0]) == 321 && TMath::Abs(pdgProng[1]) == 211){
	  ((TH1F*)listout->FindObject(fillthisK))->Fill(part->PtProng(0));
	  ((TH1F*)listout->FindObject(fillthispi))->Fill(part->PtProng(1));
	}
      }
      
      //no mass cut ditributions: mass
      fillthis="hMassS_";
      fillthis+=ptbin;
      
      if (((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421){//D0
	((TH1F*)listout->FindObject(fillthis))->Fill(minvD0);
      }
      else { //D0bar
	((TH1F*)listout->FindObject(fillthis))->Fill(minvD0bar);
      }

      //apply cut on invariant mass on the pair
      if(TMath::Abs(minvD0-mPDG)<invmasscut || TMath::Abs(minvD0bar-mPDG)<invmasscut){
	
	if(fArray==1) cout<<"LS signal: ERROR"<<endl;
	for (Int_t iprong=0; iprong<2; iprong++){
	  AliAODTrack *prong=(AliAODTrack*)part->GetDaughter(iprong);
	  labprong[iprong]=prong->GetLabel();
	  
	  //cout<<"prong name = "<<prong->GetName()<<" label = "<<prong->GetLabel()<<endl;
	  if(labprong[iprong]>=0)  mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
	  Int_t pdgprong=mcprong->GetPdgCode();
	  
	  if(TMath::Abs(pdgprong)==211) {
	    //cout<<"pi"<<endl;
	    
	    fillthispi="hptpiS_";
	    fillthispi+=ptbin;
	    ((TH1F*)listout->FindObject(fillthispi))->Fill(part->PtProng(iprong));
	    fillthispi="hd0piS_";
	    fillthispi+=ptbin;
	    ((TH1F*)listout->FindObject(fillthispi))->Fill(part->Getd0Prong(iprong));
	    if(d0[iprong] != -99) {

	      fillthispi="hd0vpiS_";
	      fillthispi+=ptbin;
	      ((TH1F*)listout->FindObject(fillthispi))->Fill(d0[iprong]);
	    }

	  }
	  
	  if(TMath::Abs(pdgprong)==321) {
	    //cout<<"kappa"<<endl;
	    
	    fillthisK="hptKS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(part->PtProng(iprong));
	    fillthisK="hd0KS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(part->Getd0Prong(iprong));
	    if (d0[iprong] != -99){
	      fillthisK="hd0vKS_";
	      fillthisK+=ptbin;
	      ((TH1F*)listout->FindObject(fillthisK))->Fill(d0[iprong]);
	    }
	  }

	  fillthis="hcosthpointd0S_";
	  fillthis+=ptbin;	  
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle(),part->Getd0Prong(iprong));

	} //end loop on prongs

	fillthis="hdcaS_";
	fillthis+=ptbin;	  
	((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDCA());

	fillthis="hcosthetapointS_";
	fillthis+=ptbin;	  
	((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle());

	fillthis="hcosthpointd0d0S_";
	fillthis+=ptbin;	  
	((TH2F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle(),part->Prodd0d0());
	  
	fillthis="hcosthetastarS_";
	fillthis+=ptbin;
	if (((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421)((TH1F*)listout->FindObject(fillthis))->Fill(part->CosThetaStarD0());
	else ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosThetaStarD0bar());

	fillthis="hd0d0S_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Prodd0d0());

	if(d0[0] != -99 && d0[1] != -99 ){
	  fillthis="hd0d0vS_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]*d0[1]);
	}

	fillthis="hdeclS_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->DecayLength());

	fillthis="hnormdeclS_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->NormalizedDecayLength());

	if(vtxProp) {
	  fillthis="hdeclvS_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->AliAODRecoDecay::DecayLength(vtxProp));

	  fillthis="hnormdeclvS_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->AliAODRecoDecay::NormalizedDecayLength(vtxProp));
	}
     } //end mass cut
    
    } else{ //Background or LS
      //cout<<"is background"<<endl;
     
      //no mass cut distributions: mass, ptbis
      fillthis="hMassB_";
      fillthis+=ptbin;
      
      if (!fCutOnDistr || (fCutOnDistr && (isSelected==1 || isSelected==3))) ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0);
      if (!fCutOnDistr || (fCutOnDistr && isSelected>1)) ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0bar);

      fillthis="hptB1prongnoMcut_";
      fillthis+=ptbin;
      
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(0));
      
      fillthis="hptB2prongsnoMcut_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(0));
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(1));

      //apply cut on invariant mass on the pair
      if(TMath::Abs(minvD0-mPDG)<invmasscut || TMath::Abs(minvD0bar-mPDG)<invmasscut){


	AliAODTrack *prong=(AliAODTrack*)part->GetDaughter(0);
	if(!prong) {
	  if(fDebug>2) cout<<"No daughter found";
	  return;
	}
	else{
	  if(prong->Charge()==1) {fTotPosPairs[ptbin]++;} else {fTotNegPairs[ptbin]++;}
	}
	
	//normalise pt distr to half afterwards
	fillthis="hptB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(0));
	((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(1));

	fillthis="hd0p0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(0));
	fillthis="hd0p1B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(1));
	fillthis="hd0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(0));
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(1));

	fillthis="hd0vp0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]);
	fillthis="hd0vp1B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0[1]);

	fillthis="hd0vB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]);
	((TH1F*)listout->FindObject(fillthis))->Fill(d0[1]);
  

	fillthis="hcosthpointd0B_";
	fillthis+=ptbin;	  
	((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle(),part->Getd0Prong(0));
	((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle(),part->Getd0Prong(1));


	fillthis="hdcaB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDCA());

	fillthis="hcosthetastarB_";
	fillthis+=ptbin;
	if (!fCutOnDistr || (fCutOnDistr && (isSelected==1 || isSelected==3)))((TH1F*)listout->FindObject(fillthis))->Fill(part->CosThetaStarD0());
	if (!fCutOnDistr || (fCutOnDistr && isSelected>1))((TH1F*)listout->FindObject(fillthis))->Fill(part->CosThetaStarD0bar());	

	fillthis="hd0d0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Prodd0d0());

	if(d0[0] != -99 && d0[1] != -99 ){
	  fillthis="hd0d0vB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]*d0[1]);
	}

	fillthis="hcosthetapointB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle());

	fillthis="hcosthpointd0d0B_";
	fillthis+=ptbin;
	((TH2F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle(),part->Prodd0d0());

	fillthis="hdeclB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->DecayLength());

	fillthis="hnormdeclB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->NormalizedDecayLength());

	if(vtxProp) {
	  fillthis="hdeclvB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->AliAODRecoDecay::DecayLength(vtxProp));

	  fillthis="hnormdeclvB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->AliAODRecoDecay::NormalizedDecayLength(vtxProp));
	}
     }//mass cut	
    }//else (background)
  }
  return;
}
//____________________________________________________________________________

void AliAnalysisTaskSED0Mass::FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi* cuts, TList *listout){
  //
  // function used in UserExec to fill mass histograms:
  //


  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  Int_t isSelected=cuts->IsSelected(part,AliRDHFCuts::kCandidate); //selected
  //cout<<"is selected = "<<isSelected<<endl;

  //cout<<"check cuts = "<<endl;
  //cuts->PrintAll();
  if (!isSelected){
    //cout<<"Not Selected"<<endl;
    return;
  }

  if(fDebug>2)  cout<<"Candidate selected"<<endl;

  Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
  //printf("SELECTED\n");
  Int_t ptbin=cuts->PtBin(part->Pt());

  AliAODTrack *prong=(AliAODTrack*)part->GetDaughter(0);
  if(!prong) {
    AliDebug(2,"No daughter found");
    return;
  }
  else{
    if(prong->Charge()==1) {fTotPosPairs[ptbin]++;} else {fTotNegPairs[ptbin]++;}
  }

  TString fillthis="";
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t labD0=-1;
  if (fReadMC) labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

  //count candidates selected by cuts
  fNentries->Fill(1);
  //count true D0 selected by cuts
  if (fReadMC && labD0>=0) fNentries->Fill(2);
  //PostData(3,fNentries);

  if (isSelected==1 || isSelected==3) { //D0
    fillthis="histMass_";
    fillthis+=ptbin;
    //cout<<"Filling "<<fillthis<<endl;

    //printf("Fill mass with D0");
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
 
    if(labD0>=0) {
      if(fArray==1) cout<<"LS signal ERROR"<<endl;

      AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
      Int_t pdgD0 = partD0->GetPdgCode();
      //cout<<"pdg = "<<pdgD0<<endl;
      if (pdgD0==421){ //D0
	//cout<<"Fill S with D0"<<endl;
	fillthis="histSgn_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	if(TMath::Abs(invmassD0 - mPDG) < 0.027){
	  fillthis="histSgn27_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	}
      } else{ //it was a D0bar
	fillthis="histRfl_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      }
    } else {//background
      fillthis="histBkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
    }
      
  }
  if (isSelected>1) { //D0bar
    fillthis="histMass_";
    fillthis+=ptbin;
    //printf("Fill mass with D0bar");
    ((TH1F*)listout->FindObject(fillthis))->Fill(invmassD0bar);
 
    if(labD0>=0) {
      if(fArray==1) cout<<"LS signal ERROR"<<endl;
      AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
      Int_t pdgD0 = partD0->GetPdgCode();
      //cout<<" pdg = "<<pdgD0<<endl;
      if (pdgD0==-421){ //D0bar
	fillthis="histSgn_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	if (TMath::Abs(invmassD0bar - mPDG) < 0.027){
	  fillthis="histSgn27_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	}

	  
      } else{
	fillthis="histRfl_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
      }
    } else {//background or LS
      fillthis="histBkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
    }
  }

  return;
}

//__________________________________________________________________________
AliAODVertex* AliAnalysisTaskSED0Mass::GetPrimaryVtxSkipped(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *d){
  //Calculate the primary vertex w/o the daughter tracks of the candidate
  
  AliESDVertex *vertexESD=0x0;
  AliAODVertex *vertexAOD=0x0;
  AliVertexerTracks *vertexer = new AliVertexerTracks(aodev->GetMagneticField());
  
  Int_t skipped[2];
  Int_t nTrksToSkip=2;
  AliAODTrack *dgTrack = (AliAODTrack*)d->GetDaughter(0);
  if(!dgTrack){
    AliDebug(2,"no daughter found!");
    return 0x0;
  }
  skipped[0]=dgTrack->GetID();
  dgTrack = (AliAODTrack*)d->GetDaughter(1);
  if(!dgTrack){
    AliDebug(2,"no daughter found!");
    return 0x0;
  }
  skipped[1]=dgTrack->GetID();

 
  //
  vertexer->SetSkipTracks(nTrksToSkip,skipped);
  vertexer->SetMinClusters(4);  
  vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(aodev); 
  if(!vertexESD) return vertexAOD;
  if(vertexESD->GetNContributors()<=0) { 
    AliDebug(2,"vertexing failed"); 
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  
  delete vertexer; vertexer=NULL;
  
  
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;
  
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);
  return vertexAOD;
  
}


//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass: Terminate() \n");


  fOutputMass = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputMass) {     
    printf("ERROR: fOutputMass not available\n");
    return;
  }
  fDistr = dynamic_cast<TList*> (GetOutputData(2));
  if (!fDistr) {
    printf("ERROR: fDistr not available\n");
    return;
  }
 
  fNentries = dynamic_cast<TH1F*>(GetOutputData(3));
  
  if(!fNentries){
    printf("ERROR: fNEntries not available\n");
    return;
  }

  fChecks = dynamic_cast<TList*> (GetOutputData(4));
  if (!fChecks) {
    printf("ERROR: fChecks not available\n");
    return;
  }
  
 
  for(Int_t ipt=0;ipt<5;ipt++){ //change 5 in GetNPtBins when sure it is written and check


    if(fArray==1){ 
      fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs[ipt]*fTotNegPairs[ipt]); 


      if(fLsNormalization>1e-6) {
	
	TString massName="histMass_";
	massName+=ipt;
	((TH1F*)fOutputMass->FindObject(massName))->Scale((1/fLsNormalization)*((TH1F*)fOutputMass->FindObject(massName))->GetEntries());

      }
    

      fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs[4]*fTotNegPairs[4]);

      if(fLsNormalization>1e-6) {

	TString nameDistr="hptB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hdcaB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthetastarB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hd0B_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hd0d0B_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthetapointB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthpointd0d0B_";
	nameDistr+=ipt;
	((TH2F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH2F*)fDistr->FindObject(nameDistr))->GetEntries());

      }
    }
  }
  TString cvname,cstname;

  if (fArray==0){
    cvname="D0invmass";
    cstname="cstat0";
  } else {
    cvname="LSinvmass";
    cstname="cstat1";
  }

  TCanvas *cMass=new TCanvas(cvname,cvname);
  cMass->cd();
  ((TH1F*)fOutputMass->FindObject("histMass_3"))->Draw();

  TCanvas* cStat=new TCanvas(cstname,Form("Stat%s",fArray ? "LS" : "D0"));
  cStat->cd();
  cStat->SetGridy();
  fNentries->Draw("htext0");

  // TCanvas *ccheck=new TCanvas(Form("cc%d",fArray),Form("cc%d",fArray));
  // ccheck->cd();

  return;
}

