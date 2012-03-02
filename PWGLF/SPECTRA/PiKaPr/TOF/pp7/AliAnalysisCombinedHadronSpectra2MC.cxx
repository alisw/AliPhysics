/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Analysis for identified charged hadron spectra. TPC/TOF               //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////




#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisFilter.h"
#include "AliCFContainer.h"
#include "AliCFFrame.h"
#include "AliCFGridSparse.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisFilter.h"
#include "AliCFContainer.h"
#include "AliCFFrame.h"
#include "AliCFGridSparse.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TSystem.h"
#include "AliMultiplicity.h"


#include "AliLog.h"
#include "AliESDUtils.h"

#include "AliAnalysisCombinedHadronSpectra2MC.h"


ClassImp(AliAnalysisCombinedHadronSpectra2MC)

//________________________________________________________________________
AliAnalysisCombinedHadronSpectra2MC::AliAnalysisCombinedHadronSpectra2MC() 
: AliAnalysisTaskSE("TaskChargedHadron"), multiplicity(0),vert(0), fESD(0), fListHist(0), fESDtrackCuts(0),fESDpid(0),
  fMCtrue(0),
  fAlephParameters(),
  TOFCheck(0),  
  calibrateESD(kTRUE), correctTExp(kTRUE), useT0TOF(kTRUE), timeResolution(96.), tuneTOFMC(kFALSE),
  fTreeTrack(0),
  fTreeEv(0),
  fLoadOCDB(kTRUE),
  frunOld(0),
  frun(0),
  tofCalib(0),
  t0maker(0),
  fT0TOF0(-999), fT0TOF1(-999), fT0TOF2(-999), fT0TOF3(-999), fT0TOF4(-999), fT0TOF5(-999), fT0TOF6(-999), fT0TOF7(-999), XPrimVertex(-999), YPrimVertex(-999), ZPrimVertex(-999), NContrPrimVertex(-999), rapidityMC(0), fDCAXY(-999), fDCAZ(-999), fcut(-999), fTOFout(-999), ftrdout(-999), ftime(-999), ftpcclust(-999), flength(-999), fsign(-999), ftimetof(-999), ftofchan(-999), feta(-999), fphi(-999), fmomtrasv(-999),sigmapi(-999), sigmaka(-999), sigmapr(-999), fTot(-999), fmom(-999), fexptimepi(-999), fexptimeka(-999), fexptimepr(-999), ftofz(-999),ftofx(-999), t0track(-999), TPCSignal(-999), TPCSigmaPI(-999), TPCSigmaKA(-999), TPCSigmaPR(-999), fmatch(-999), fPhiout(-999), fXout(-999), fYout(-999), fZout(-999), fTimeZeroType(AliESDpid::kTOF_T0), fMCtracks(-999), fMCPrimaries(-999), spdCorr(-999), treeMCP(-999), treeMCPt(-999), treeMCEta(-999), treeMCPhi(-999), treeMCPdg(-999), treeMCPBis(-999), treeMCPtBis(-999), treeMCEtaBis(-999), treeMCPhiBis(-999), treeMCPdgBis(-999), t0trackSigma(-999), fptMC(-999), fphiMC(-999), fetaMC(-999), fPdgcode(-999),
pad (0x0), resx(0x0), resz(0x0), tofres(0x0), tofresTOF(0x0), tofresgood (0x0), hNumMatch(0x0), hNumMatchPos(0x0),  hNumMatchNeg(0x0),  hDenMatch(0x0),  
  hNumMatchPip(0x0), hNumMatchPim(0x0),hNumMatchKap(0x0),  hNumMatchKam(0x0),hNumMatchPrp(0x0),  hNumMatchPrm(0x0), hDenMatchPip(0x0),  hDenMatchPim(0x0), hDenMatchKap(0x0),  hDenMatchKam(0x0), hDenMatchPrp(0x0), hDenMatchPrm(0x0),

hDenMatchPos(0x0),  hDenMatchNeg(0x0),  hNumMatchEta(0x0),  hNumMatchPosEta(0x0),  hNumMatchNegEta(0x0),  hDenMatchEta(0x0),  hDenMatchPosEta(0x0),   hDenMatchNegEta(0x0),  hNumMatchphiOut(0x0),  hNumMatchPosphiOut(0x0),  hNumMatchNegphiOut(0x0),  hDenMatchphiOut(0x0),  hDenMatchPosphiOut(0x0),  hDenMatchNegphiOut(0x0),  hNumMatchEtaPtMa(0x0),  hNumMatchPosEtaPtMa(0x0),  hNumMatchNegEtaPtMa(0x0),  hDenMatchEtaPtMa(0x0),  hDenMatchPosEtaPtMa(0x0),  hDenMatchNegEtaPtMa(0x0),  hNumMatchphiOutPtMa(0x0),  hNumMatchPosphiOutPtMa(0x0),  hNumMatchNegphiOutPtMa(0x0),  hDenMatchphiOutPtMa(0x0),  hDenMatchPosphiOutPtMa(0x0),  hDenMatchNegphiOutPtMa(0x0), hNumMatchTRDOut(0x0), hNumMatchPosTRDOut(0x0), hNumMatchNegTRDOut(0x0), hDenMatchTRDOut(0x0), hDenMatchPosTRDOut(0x0), hDenMatchNegTRDOut(0x0), hNumMatchNoTRDOut(0x0), hNumMatchPosNoTRDOut(0x0), hNumMatchNegNoTRDOut(0x0), hDenMatchNoTRDOut(0x0), hDenMatchPosNoTRDOut(0x0), hDenMatchNegNoTRDOut(0x0), hNumMatchTPCpip(0x0), hNumMatchTPCkap(0x0), hNumMatchTPCprp(0x0), hDenMatchTPCpip(0x0), hDenMatchTPCkap(0x0), hDenMatchTPCprp(0x0), hNumMatchTPCpim(0x0), hNumMatchTPCkam(0x0), hNumMatchTPCprm(0x0), hDenMatchTPCpim(0x0), hDenMatchTPCkam(0x0), hDenMatchTPCprm(0x0), hNumEv(0x0) 
{
  // default Constructor
  for(Int_t i=0; i<5;i++){r1[i]=-999;}
  
  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrk[mult][part]=(0);
      hDenMatchMultTrk[mult][part]=(0);
      hDenTrkMultTrk[mult][part]=(0);
    }
  }
  
  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrkTRDOut[mult][part]=(0);
      hDenMatchMultTrkTRDOut[mult][part]=(0);
      hDenTrkMultTrkTRDOut[mult][part]=(0);
    }
  }

  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrkNoTRDOut[mult][part]=(0);
      hDenMatchMultTrkNoTRDOut[mult][part]=(0);
      hDenTrkMultTrkNoTRDOut[mult][part]=(0);
    }
  }

  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultSPD[mult][part]=(0);
      hDenMatchMultSPD[mult][part]=(0);
      hDenTrkMultSPD[mult][part]=(0);
    }
  }
  
  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<2;part++){
      hNumMatchMultTrkInc[mult][part]=(0);
      hDenMatchMultTrkInc[mult][part]=(0);
      hNumMatchMultSPDInc[mult][part]=(0);
      hDenMatchMultSPDInc[mult][part]=(0);
    }
  }
  
  for(Int_t part=0; part<6;part++){
    hDenTrkVertMultTrk[part]=(0);
  }
  for(Int_t part=0; part<6;part++){
    hDenTrkTriggerMultTrk[part]=(0);
  }
  
}


//________________________________________________________________________
AliAnalysisCombinedHadronSpectra2MC::AliAnalysisCombinedHadronSpectra2MC(const char *name)  : AliAnalysisTaskSE(name), multiplicity(0),vert(0), fESD(0), fListHist(0), fESDtrackCuts(0),fESDpid(0),
  fMCtrue(0),
  fAlephParameters(),
  TOFCheck(0),  
  calibrateESD(kTRUE), correctTExp(kTRUE), useT0TOF(kTRUE), timeResolution(96.), tuneTOFMC(kFALSE),
  fTreeTrack(0),
  fTreeEv(0),
  fLoadOCDB(kTRUE),
  frunOld(0),
  frun(0),
  tofCalib(0),
  t0maker(0),
											      fT0TOF0(-999), fT0TOF1(-999), fT0TOF2(-999), fT0TOF3(-999), fT0TOF4(-999), fT0TOF5(-999), fT0TOF6(-999), fT0TOF7(-999), XPrimVertex(-999), YPrimVertex(-999), ZPrimVertex(-999), NContrPrimVertex(-999), rapidityMC(0), fDCAXY(-999), fDCAZ(-999), fcut(-999), fTOFout(-999), ftrdout(-999), ftime(-999), ftpcclust(-999), flength(-999), fsign(-999), ftimetof(-999), ftofchan(-999), feta(-999), fphi(-999), fmomtrasv(-999),sigmapi(-999), sigmaka(-999), sigmapr(-999), fTot(-999), fmom(-999), fexptimepi(-999), fexptimeka(-999), fexptimepr(-999), ftofz(-999),ftofx(-999), t0track(-999), TPCSignal(-999), TPCSigmaPI(-999), TPCSigmaKA(-999), TPCSigmaPR(-999), fmatch(-999), fPhiout(-999), fXout(-999), fYout(-999), fZout(-999), fTimeZeroType(AliESDpid::kTOF_T0), fMCtracks(-999), fMCPrimaries(-999), spdCorr(-999), treeMCP(-999), treeMCPt(-999), treeMCEta(-999), treeMCPhi(-999), treeMCPdg(-999), treeMCPBis(-999), treeMCPtBis(-999), treeMCEtaBis(-999), treeMCPhiBis(-999), treeMCPdgBis(-999), t0trackSigma(-999), fptMC(-999), fphiMC(-999), fetaMC(-999), fPdgcode(-999),
pad (0x0), resx(0x0), resz(0x0), tofres(0x0), tofresTOF(0x0), tofresgood (0x0), hNumMatch(0x0), hNumMatchPos(0x0),  hNumMatchNeg(0x0),  hDenMatch(0x0), 
 hNumMatchPip(0x0), hNumMatchPim(0x0),hNumMatchKap(0x0),  hNumMatchKam(0x0),hNumMatchPrp(0x0),  hNumMatchPrm(0x0), hDenMatchPip(0x0),  hDenMatchPim(0x0), hDenMatchKap(0x0),  hDenMatchKam(0x0), hDenMatchPrp(0x0), hDenMatchPrm(0x0),
 hDenMatchPos(0x0),  hDenMatchNeg(0x0),  hNumMatchEta(0x0),  hNumMatchPosEta(0x0),  hNumMatchNegEta(0x0),  hDenMatchEta(0x0),  hDenMatchPosEta(0x0),   hDenMatchNegEta(0x0),  hNumMatchphiOut(0x0),  hNumMatchPosphiOut(0x0),  hNumMatchNegphiOut(0x0),  hDenMatchphiOut(0x0),  hDenMatchPosphiOut(0x0),  hDenMatchNegphiOut(0x0),  hNumMatchEtaPtMa(0x0),  hNumMatchPosEtaPtMa(0x0),  hNumMatchNegEtaPtMa(0x0),  hDenMatchEtaPtMa(0x0),  hDenMatchPosEtaPtMa(0x0),  hDenMatchNegEtaPtMa(0x0),  hNumMatchphiOutPtMa(0x0),  hNumMatchPosphiOutPtMa(0x0),  hNumMatchNegphiOutPtMa(0x0),  hDenMatchphiOutPtMa(0x0),  hDenMatchPosphiOutPtMa(0x0),  hDenMatchNegphiOutPtMa(0x0), hNumMatchTRDOut(0x0), hNumMatchPosTRDOut(0x0), hNumMatchNegTRDOut(0x0), hDenMatchTRDOut(0x0), hDenMatchPosTRDOut(0x0), hDenMatchNegTRDOut(0x0), hNumMatchNoTRDOut(0x0), hNumMatchPosNoTRDOut(0x0), hNumMatchNegNoTRDOut(0x0), hDenMatchNoTRDOut(0x0), hDenMatchPosNoTRDOut(0x0), hDenMatchNegNoTRDOut(0x0), hNumMatchTPCpip(0x0), hNumMatchTPCkap(0x0), hNumMatchTPCprp(0x0), hDenMatchTPCpip(0x0), hDenMatchTPCkap(0x0), hDenMatchTPCprp(0x0), hNumMatchTPCpim(0x0), hNumMatchTPCkam(0x0), hNumMatchTPCprm(0x0), hDenMatchTPCpim(0x0), hDenMatchTPCkam(0x0), hDenMatchTPCprm(0x0),  hNumEv(0x0)
{
  
  //
  // standard constructur which should be used
  //
  Printf("*** CONSTRUCTOR CALLED ****");


  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrk[mult][part]=(0);
      hDenMatchMultTrk[mult][part]=(0);
      hDenTrkMultTrk[mult][part]=(0);
    }
  }

  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrkTRDOut[mult][part]=(0);
      hDenMatchMultTrkTRDOut[mult][part]=(0);
      hDenTrkMultTrkTRDOut[mult][part]=(0);
    }
  }

  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrkNoTRDOut[mult][part]=(0);
      hDenMatchMultTrkNoTRDOut[mult][part]=(0);
      hDenTrkMultTrkNoTRDOut[mult][part]=(0);
    }
  }

  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultSPD[mult][part]=(0);
      hDenMatchMultSPD[mult][part]=(0);
      hDenTrkMultSPD[mult][part]=(0);
    }
  }
  
  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<2;part++){
      hNumMatchMultTrkInc[mult][part]=(0);
      hDenMatchMultTrkInc[mult][part]=(0);
      hNumMatchMultSPDInc[mult][part]=(0);
      hDenMatchMultSPDInc[mult][part]=(0);
      }
  }

  for(Int_t part=0; part<6;part++){
    hDenTrkVertMultTrk[part]=(0);
  }
  for(Int_t part=0; part<6;part++){
    hDenTrkTriggerMultTrk[part]=(0);
  }
  

  // 
  /* real */
  //   for(Int_t i=0; i<5;i++){r1[i]=-999;}
  //   fAlephParameters[0] = 0.0283086;
  //   fAlephParameters[1] = 2.63394e+01;
  //   fAlephParameters[2] = 5.04114e-11;
  //   fAlephParameters[3] = 2.12543e+00;
  //   fAlephParameters[4] = 4.88663e+00;
  //
  // initialize PID object
  //
  
  //   tofCalib = new AliTOFcalib();
  
  //   fESDpid = new AliESDpid();
  
  //   t0maker = new AliTOFT0maker(fESDpid, tofCalib); 
  
  //   if(AliPID::ParticleMass(0) == 0) new AliPID(); 
  
  //
  // create track cuts
  //
  AliESDtrackCuts* ESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCuts = (AliESDtrackCuts*) ESDtrackCuts->GetStandardITSTPCTrackCuts2010(kTRUE);
  
  
  DefineOutput(1, TList::Class());

}


//________________________________________________________________________
Int_t AliAnalysisCombinedHadronSpectra2MC::Mult()
{
  AliESDtrackCuts* esdTrackCutsMult = new AliESDtrackCuts;
  
  // TPC
  esdTrackCutsMult->SetMinNClustersTPC(70);
  esdTrackCutsMult->SetMaxChi2PerClusterTPC(4);
  esdTrackCutsMult->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsMult->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCutsMult->SetRequireITSRefit(kTRUE);
  esdTrackCutsMult->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					     AliESDtrackCuts::kAny);
  
  esdTrackCutsMult->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  
  esdTrackCutsMult->SetMaxDCAToVertexZ(2);
  
  esdTrackCutsMult->SetDCAToVertex2D(kFALSE);
  esdTrackCutsMult->SetRequireSigmaToVertex(kFALSE);
  
  esdTrackCutsMult->SetEtaRange(-0.8,+0.8);
  esdTrackCutsMult->SetPtRange(0.15, 1e10);
  
  // Here is what we finally want:
  Int_t multipli = esdTrackCutsMult->CountAcceptedTracks(fESD);
  delete esdTrackCutsMult;
  return multipli;
  
 }


//________________________________________________________________________
void AliAnalysisCombinedHadronSpectra2MC::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  
  //OpenFile(2);

  fTreeTrack = new TTree("TreeTrack","Track Properties");

  fTreeTrack->Branch("XPrimVertex",&XPrimVertex,"XPrimVertex/D");  
  fTreeTrack->Branch("YPrimVertex",&YPrimVertex,"YPrimVertex/D");  
  fTreeTrack->Branch("ZPrimVertex",&ZPrimVertex,"ZPrimVertex/D"); 
  
  fTreeTrack->Branch("NContrPrimVertex",&NContrPrimVertex,"NContrPrimVertex/I");  
  fTreeTrack->Branch("multiplicity",&multiplicity,"multiplicity/I");  
  fTreeTrack->Branch("T0TOF0",&fT0TOF0,"T0TOF0/D");
  fTreeTrack->Branch("T0TOF1",&fT0TOF1,"T0TOF1/D");
  fTreeTrack->Branch("T0TOF7",&fT0TOF7,"T0TOF7/D");
  fTreeTrack->Branch("DCAXY",&fDCAXY,"fDCAXY/F");  
  fTreeTrack->Branch("DCAZ",&fDCAZ,"fDCAZ/F");  
  fTreeTrack->Branch("cut",&fcut,"cut/I");
  fTreeTrack->Branch("TOFout",&fTOFout,"TOFout/I");
  fTreeTrack->Branch("ftrdout",&ftrdout,"ftrdout/I");
  fTreeTrack->Branch("ftime",&ftime,"ftime/I");
  fTreeTrack->Branch("fmatch",&fmatch,"fmatch/I");
  fTreeTrack->Branch("fmom",&fmom,"fmom/D");
  fTreeTrack->Branch("tpcclust",&ftpcclust,"tpcclust/I");
  fTreeTrack->Branch("length",&flength,"length/D");
  fTreeTrack->Branch("sign",&fsign,"sign/I");
  fTreeTrack->Branch("timetof",&ftimetof,"timetof/D");
  fTreeTrack->Branch("exptimepi",&fexptimepi,"exptimepi/D");
  fTreeTrack->Branch("exptimeka",&fexptimeka,"exptimeka/D");
  fTreeTrack->Branch("exptimepr",&fexptimepr,"exptimepr/D");
  fTreeTrack->Branch("tofz",&ftofz,"tofz/D"); 
  fTreeTrack->Branch("tofx",&ftofx,"tofx/D");
  fTreeTrack->Branch("t0track",&t0track,"t0track/F");
  fTreeTrack->Branch("tofchan",&ftofchan,"tofchan/I");
  fTreeTrack->Branch("eta",&feta,"eta/D");
  fTreeTrack->Branch("phi",&fphi,"phi/D");
  fTreeTrack->Branch("TOFtot",&fTot,"TOFtot/F");
  fTreeTrack->Branch("momtrasv",&fmomtrasv,"momtrasv/D");
  fTreeTrack->Branch("sigmapi",&sigmapi,"sigmapi/D");
  fTreeTrack->Branch("sigmaka",&sigmaka,"sigmaka/D");
  fTreeTrack->Branch("sigmapr",&sigmapr,"sigmapr/D");
  fTreeTrack->Branch("TPCsignal",&TPCSignal,"TPCsignal/D");
  fTreeTrack->Branch("TPCSigmaPI",&TPCSigmaPI,"TPCSigmaPI/F");
  fTreeTrack->Branch("TPCSigmaKA",&TPCSigmaKA,"TPCSigmaKA/F");
  fTreeTrack->Branch("TPCSigmaPR",&TPCSigmaPR,"TPCSigmaPR/F");
  fTreeTrack->Branch("r10",&r1[0],"r10/D");
  fTreeTrack->Branch("r11",&r1[1],"r11/D");
  fTreeTrack->Branch("r12",&r1[2],"r12/D");
  fTreeTrack->Branch("r13",&r1[3],"r13/D");
  fTreeTrack->Branch("r14",&r1[4],"r14/D");
  
  fTreeTrack->Branch("fPhiout",&fPhiout,"fPhiout/D");
  fTreeTrack->Branch("fXout",&fXout,"fXout/D");
  fTreeTrack->Branch("fYout",&fYout,"fYout/D");
  fTreeTrack->Branch("fZout",&fZout,"fZout/D");

  //PostData(2, fTreeTrack);

  //OpenFile(1);

  fTreeEv = new TTree("TreeEv","Event Properties");
  fTreeEv->Branch("multiplicity",&multiplicity,"multiplicity/I");  
  fTreeEv->Branch("vert",&vert,"vert/I");  

  //PostData(1, fTreeEv);

  OpenFile(1);

  if (!TOFCheck) TOFCheck = new TList();
  TOFCheck->SetOwner();

  Double_t fBinLim0[47];
  for(Int_t i=0;i<=46;i++){
    // pt [0.20,5] GeV/c - not uniform
    if(i<=16) fBinLim0[i] = 0.2 + i*0.05;
    if((i>16)&&(i<=36))fBinLim0[i] = 1.0 + (i-16)*0.1;
    if(i>36) fBinLim0[i] = 3.0 + (i-36)*0.2;
  }


  pad=new TH2D("pad", "pad", 400, -10, 10, 400, -10, 10);
  pad->GetXaxis()->SetTitle("#DeltaX_{pad} (cm)");
  pad->GetYaxis()->SetTitle("#DeltaZ_{pad} (cm)");
  TOFCheck->AddLast(pad);
  resx=new TH1D("resx", "resx", 150, -10, 10);
  resx->GetXaxis()->SetTitle("#DeltaX_{pad} (cm)");
  TOFCheck->AddLast(resx);
  resz=new TH1D("resz", "resz", 150, -10, 10);
  resz->GetXaxis()->SetTitle("#DeltaZ_{pad} (cm)");
  TOFCheck->AddLast(resz);
  tofres=new TH1D("tofres", "tofres", 100, -500, 500);
  tofres->GetXaxis()->SetTitle("t_{TOF}-t_{0}-t_{exp #pi} (ps)");
  TOFCheck->AddLast(tofres);
  tofresTOF=new TH1D("tofresTOF", "tofresTOF", 100, -500, 500);
  tofresTOF->GetXaxis()->SetTitle("t_{TOF}-t_{0}-t_{exp #pi} (ps)");
  TOFCheck->AddLast(tofresTOF);
  tofresgood=new TH1D("tofresgood", "tofresgood", 100, -500, 500);
  tofresgood->Sumw2();
  tofresgood->GetXaxis()->SetTitle("t_{TOF}-t_{0}-t_{exp #pi} (ps)");
  TOFCheck->AddLast(tofresgood);

  //PT

  hNumMatch=new TH1F("hNumMatch","",46,fBinLim0);
  hNumMatch->Sumw2();
  TOFCheck->AddLast(hNumMatch);
  hNumMatchPos=new TH1F("hNumMatchPos","",46,fBinLim0);
  hNumMatchPos->Sumw2();
  TOFCheck->AddLast(hNumMatchPos);
  hNumMatchNeg=new TH1F("hNumMatchNeg","",46,fBinLim0);
  hNumMatchNeg->Sumw2();
  TOFCheck->AddLast(hNumMatchNeg);
  hDenMatch=new TH1F("hDenMatch","",46,fBinLim0);
  hDenMatch->Sumw2();
  TOFCheck->AddLast(hDenMatch);
  hDenMatchPos=new TH1F("hDenMatchPos","",46,fBinLim0);
  hDenMatchPos->Sumw2();
  TOFCheck->AddLast(hDenMatchPos);
  hDenMatchNeg=new TH1F("hDenMatchNeg","",46,fBinLim0);
  hDenMatchNeg->Sumw2();
  TOFCheck->AddLast(hDenMatchNeg);

  //PT per SPECIE
  
  hNumMatchPip=new TH1F("hNumMatchPip","",46,fBinLim0);
  hNumMatchPip->Sumw2();
  TOFCheck->AddLast(hNumMatchPip);
  hNumMatchPim=new TH1F("hNumMatchPim","",46,fBinLim0);
  hNumMatchPim->Sumw2();
  TOFCheck->AddLast(hNumMatchPim);
  hNumMatchKap=new TH1F("hNumMatchKap","",46,fBinLim0);
  hNumMatchKap->Sumw2();
  TOFCheck->AddLast(hNumMatchKap);
  hNumMatchKam=new TH1F("hNumMatchKam","",46,fBinLim0);
  hNumMatchKam->Sumw2();
  TOFCheck->AddLast(hNumMatchKam);
  hNumMatchPrp=new TH1F("hNumMatchPrp","",46,fBinLim0);
  hNumMatchPrp->Sumw2();
  TOFCheck->AddLast(hNumMatchPrp);
  hNumMatchPrm=new TH1F("hNumMatchPrm","",46,fBinLim0);
  hNumMatchPrm->Sumw2();
  TOFCheck->AddLast(hNumMatchPrm);

  hDenMatchPip=new TH1F("hDenMatchPip","",46,fBinLim0);
  hDenMatchPip->Sumw2();
  TOFCheck->AddLast(hDenMatchPip);
  hDenMatchPim=new TH1F("hDenMatchPim","",46,fBinLim0);
  hDenMatchPim->Sumw2();
  TOFCheck->AddLast(hDenMatchPim);
  hDenMatchKap=new TH1F("hDenMatchKap","",46,fBinLim0);
  hDenMatchKap->Sumw2();
  TOFCheck->AddLast(hDenMatchKap);
  hDenMatchKam=new TH1F("hDenMatchKam","",46,fBinLim0);
  hDenMatchKam->Sumw2();
  TOFCheck->AddLast(hDenMatchKam);
  hDenMatchPrp=new TH1F("hDenMatchPrp","",46,fBinLim0);
  hDenMatchPrp->Sumw2();
  TOFCheck->AddLast(hDenMatchPrp);
  hDenMatchPrm=new TH1F("hDenMatchPrm","",46,fBinLim0);
  hDenMatchPrm->Sumw2();
  TOFCheck->AddLast(hDenMatchPrm);


  //ETA

  
  hNumMatchEta=new TH1F("hNumMatchEta","",36,-0.9,0.9);
  hNumMatchEta->Sumw2();
  TOFCheck->AddLast(hNumMatchEta);
  hNumMatchPosEta=new TH1F("hNumMatchPosEta","",36,-0.9,0.9);
  hNumMatchPosEta->Sumw2();
  TOFCheck->AddLast(hNumMatchPosEta);
  hNumMatchNegEta=new TH1F("hNumMatchNegEta","",36,-0.9,0.9);
  hNumMatchNegEta->Sumw2();
  TOFCheck->AddLast(hNumMatchNegEta);
  hDenMatchEta=new TH1F("hDenMatchEta","",36,-0.9,0.9);
  hDenMatchEta->Sumw2();
  TOFCheck->AddLast(hDenMatchEta);
  hDenMatchPosEta=new TH1F("hDenMatchPosEta","",36,-0.9,0.9);
  hDenMatchPosEta->Sumw2();
  TOFCheck->AddLast(hDenMatchPosEta);
  hDenMatchNegEta=new TH1F("hDenMatchNegEta","",36,-0.9,0.9);
  hDenMatchNegEta->Sumw2();
  TOFCheck->AddLast(hDenMatchNegEta);

  //PHI OUT

  hNumMatchphiOut=new TH1F("hNumMatchphiOut","",18,0,360);
  hNumMatchphiOut->Sumw2();
  TOFCheck->AddLast(hNumMatchphiOut);
  hNumMatchPosphiOut=new TH1F("hNumMatchPosphiOut","",18,0,360);
  hNumMatchPosphiOut->Sumw2();
  TOFCheck->AddLast(hNumMatchPosphiOut);
  hNumMatchNegphiOut=new TH1F("hNumMatchNegphiOut","",18,0,360);
  hNumMatchNegphiOut->Sumw2();
  TOFCheck->AddLast(hNumMatchNegphiOut);
  hDenMatchphiOut=new TH1F("hDenMatchphiOut","",18,0,360);
  hDenMatchphiOut->Sumw2();
  TOFCheck->AddLast(hDenMatchphiOut);
  hDenMatchPosphiOut=new TH1F("hDenMatchPosphiOut","",18,0,360);
  hDenMatchPosphiOut->Sumw2();
  TOFCheck->AddLast(hDenMatchPosphiOut);
  hDenMatchNegphiOut=new TH1F("hDenMatchNegphiOut","",18,0,360);
  hDenMatchNegphiOut->Sumw2();
  TOFCheck->AddLast(hDenMatchNegphiOut);

  //ETA pt>0.5

  hNumMatchEtaPtMa=new TH1F("hNumMatchEtaPtMa","",36,-0.9,0.9);
  hNumMatchEtaPtMa->Sumw2();
  TOFCheck->AddLast(hNumMatchEtaPtMa);
  hNumMatchPosEtaPtMa=new TH1F("hNumMatchPosEtaPtMa","",36,-0.9,0.9);
  hNumMatchPosEtaPtMa->Sumw2();
  TOFCheck->AddLast(hNumMatchPosEtaPtMa);
  hNumMatchNegEtaPtMa=new TH1F("hNumMatchNegEtaPtMa","",36,-0.9,0.9);
  hNumMatchNegEtaPtMa->Sumw2();
  TOFCheck->AddLast(hNumMatchNegEtaPtMa);
  hDenMatchEtaPtMa=new TH1F("hDenMatchEtaPtMa","",36,-0.9,0.9);
  hDenMatchEtaPtMa->Sumw2();
  TOFCheck->AddLast(hDenMatchEtaPtMa);
  hDenMatchPosEtaPtMa=new TH1F("hDenMatchPosEtaPtMa","",36,-0.9,0.9);
  hDenMatchPosEtaPtMa->Sumw2();
  TOFCheck->AddLast(hDenMatchPosEtaPtMa);
  hDenMatchNegEtaPtMa=new TH1F("hDenMatchNegEtaPtMa","",36,-0.9,0.9);
  hDenMatchNegEtaPtMa->Sumw2();
  TOFCheck->AddLast(hDenMatchNegEtaPtMa);

  //PHI OUT pt>0.5

  hNumMatchphiOutPtMa=new TH1F("hNumMatchphiOutPtMa","",18,0,360);
  hNumMatchphiOutPtMa->Sumw2();
  TOFCheck->AddLast(hNumMatchphiOutPtMa);
  hNumMatchPosphiOutPtMa=new TH1F("hNumMatchPosphiOutPtMa","",18,0,360);
  hNumMatchPosphiOutPtMa->Sumw2();
  TOFCheck->AddLast(hNumMatchPosphiOutPtMa);
  hNumMatchNegphiOutPtMa=new TH1F("hNumMatchNegphiOutPtMa","",18,0,360);
  hNumMatchNegphiOutPtMa->Sumw2();
  TOFCheck->AddLast(hNumMatchNegphiOutPtMa);
  hDenMatchphiOutPtMa=new TH1F("hDenMatchphiOutPtMa","",18,0,360);
  hDenMatchphiOutPtMa->Sumw2();
  TOFCheck->AddLast(hDenMatchphiOutPtMa);
  hDenMatchPosphiOutPtMa=new TH1F("hDenMatchPosphiOutPtMa","",18,0,360);
  hDenMatchPosphiOutPtMa->Sumw2();
  TOFCheck->AddLast(hDenMatchPosphiOutPtMa);
  hDenMatchNegphiOutPtMa=new TH1F("hDenMatchNegphiOutPtMa","",18,0,360);
  hDenMatchNegphiOutPtMa->Sumw2();
  TOFCheck->AddLast(hDenMatchNegphiOutPtMa);
 
  //TRD OUT

  hNumMatchTRDOut=new TH1F("hNumMatchTRDOut","",46,fBinLim0);
  hNumMatchTRDOut->Sumw2();
  TOFCheck->AddLast(hNumMatchTRDOut);
  hNumMatchPosTRDOut=new TH1F("hNumMatchPosTRDOut","",46,fBinLim0);
  hNumMatchPosTRDOut->Sumw2();
  TOFCheck->AddLast(hNumMatchPosTRDOut);
  hNumMatchNegTRDOut=new TH1F("hNumMatchNegTRDOut","",46,fBinLim0);
  hNumMatchNegTRDOut->Sumw2();
  TOFCheck->AddLast(hNumMatchNegTRDOut);
  hDenMatchTRDOut=new TH1F("hDenMatchTRDOut","",46,fBinLim0);
  hDenMatchTRDOut->Sumw2();
  TOFCheck->AddLast(hDenMatchTRDOut);
  hDenMatchPosTRDOut=new TH1F("hDenMatchPosTRDOut","",46,fBinLim0);
  hDenMatchPosTRDOut->Sumw2();
  TOFCheck->AddLast(hDenMatchPosTRDOut);
  hDenMatchNegTRDOut=new TH1F("hDenMatchNegTRDOut","",46,fBinLim0);
  hDenMatchNegTRDOut->Sumw2();
  TOFCheck->AddLast(hDenMatchNegTRDOut);

  hNumMatchNoTRDOut=new TH1F("hNumMatchNoTRDOut","",46,fBinLim0);
  hNumMatchNoTRDOut->Sumw2();
  TOFCheck->AddLast(hNumMatchNoTRDOut);
  hNumMatchPosNoTRDOut=new TH1F("hNumMatchPosNoTRDOut","",46,fBinLim0);
  hNumMatchPosNoTRDOut->Sumw2();
  TOFCheck->AddLast(hNumMatchPosNoTRDOut);
  hNumMatchNegNoTRDOut=new TH1F("hNumMatchNegNoTRDOut","",46,fBinLim0);
  hNumMatchNegNoTRDOut->Sumw2();
  TOFCheck->AddLast(hNumMatchNegNoTRDOut);
  hDenMatchNoTRDOut=new TH1F("hDenMatchNoTRDOut","",46,fBinLim0);
  hDenMatchNoTRDOut->Sumw2();
  TOFCheck->AddLast(hDenMatchNoTRDOut);
  hDenMatchPosNoTRDOut=new TH1F("hDenMatchPosNoTRDOut","",46,fBinLim0);
  hDenMatchPosNoTRDOut->Sumw2();
  TOFCheck->AddLast(hDenMatchPosNoTRDOut);
  hDenMatchNegNoTRDOut=new TH1F("hDenMatchNegNoTRDOut","",46,fBinLim0);
  hDenMatchNegNoTRDOut->Sumw2();
  TOFCheck->AddLast(hDenMatchNegNoTRDOut);

  //confronto con TPC

  hNumMatchTPCpip=new TH1F("hNumMatchTPCpip","",46,fBinLim0);
  hNumMatchTPCpip->Sumw2();
  TOFCheck->AddLast(hNumMatchTPCpip);
  hNumMatchTPCkap=new TH1F("hNumMatchTPCkap","",46,fBinLim0);
  hNumMatchTPCkap->Sumw2();
  TOFCheck->AddLast(hNumMatchTPCkap);
  hNumMatchTPCprp=new TH1F("hNumMatchTPCprp","",46,fBinLim0);
  hNumMatchTPCprp->Sumw2();
  TOFCheck->AddLast(hNumMatchTPCprp);
  
  hDenMatchTPCpip=new TH1F("hDenMatchTPCpip","",46,fBinLim0);
  hDenMatchTPCpip->Sumw2();
  TOFCheck->AddLast(hDenMatchTPCpip);
  hDenMatchTPCkap=new TH1F("hDenMatchTPCkap","",46,fBinLim0);
  hDenMatchTPCkap->Sumw2();
  TOFCheck->AddLast(hDenMatchTPCkap);
  hDenMatchTPCprp=new TH1F("hDenMatchTPCprp","",46,fBinLim0);
  hDenMatchTPCprp->Sumw2();
  TOFCheck->AddLast(hDenMatchTPCprp);

  hNumMatchTPCpim=new TH1F("hNumMatchTPCpim","",46,fBinLim0);
  hNumMatchTPCpim->Sumw2();
  TOFCheck->AddLast(hNumMatchTPCpim);
  hNumMatchTPCkam=new TH1F("hNumMatchTPCkam","",46,fBinLim0);
  hNumMatchTPCkam->Sumw2();
  TOFCheck->AddLast(hNumMatchTPCkam);
  hNumMatchTPCprm=new TH1F("hNumMatchTPCprm","",46,fBinLim0);
  hNumMatchTPCprm->Sumw2();
  TOFCheck->AddLast(hNumMatchTPCprm);
  
  hDenMatchTPCpim=new TH1F("hDenMatchTPCpim","",46,fBinLim0);
  hDenMatchTPCpim->Sumw2();
  TOFCheck->AddLast(hDenMatchTPCpim);
  hDenMatchTPCkam=new TH1F("hDenMatchTPCkam","",46,fBinLim0);
  hDenMatchTPCkam->Sumw2();
  TOFCheck->AddLast(hDenMatchTPCkam);
  hDenMatchTPCprm=new TH1F("hDenMatchTPCprm","",46,fBinLim0);
  hDenMatchTPCprm->Sumw2();
  TOFCheck->AddLast(hDenMatchTPCprm);


  //selected events
  hNumEv=new TH1F("NumEv","NumEv",4,1,5);
  TOFCheck->AddLast(hNumEv);


  
  //efficiency
  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<2;part++){
      hNumMatchMultTrkInc[mult][part]=new TH1F(Form("hNumMatch_Inc_MultTrk%i_Charge%i",mult,part),"",46,fBinLim0);
      hNumMatchMultTrkInc[mult][part]->Sumw2();
      TOFCheck->AddLast(hNumMatchMultTrkInc[mult][part]);
      hDenMatchMultTrkInc[mult][part]=new TH1F(Form("hDenMatch_Inc_MultTrk%i_Charge%i",mult,part),"",46,fBinLim0);
      hDenMatchMultTrkInc[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenMatchMultTrkInc[mult][part]);
    }
  }
  

  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrk[mult][part]=new TH1F(Form("hNumMatch_MultTrk%i_Part%i",mult,part),"",46,fBinLim0);
      hNumMatchMultTrk[mult][part]->Sumw2();
      TOFCheck->AddLast(hNumMatchMultTrk[mult][part]);
      hDenMatchMultTrk[mult][part]=new TH1F(Form("hDenMatch_MultTrk%i_Part%i",mult,part),"",46,fBinLim0);
      hDenMatchMultTrk[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenMatchMultTrk[mult][part]);
      hDenTrkMultTrk[mult][part]=new TH1F(Form("hDenTrk_MultTrk%i_Part%i",mult,part),"",46,fBinLim0);
      hDenTrkMultTrk[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenTrkMultTrk[mult][part]);
    }
  }
  

   for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrkTRDOut[mult][part]=new TH1F(Form("hNumMatch_MultTrk%i_Part%iTRDOut",mult,part),"",46,fBinLim0);
      hNumMatchMultTrkTRDOut[mult][part]->Sumw2();
      TOFCheck->AddLast(hNumMatchMultTrkTRDOut[mult][part]);
      hDenMatchMultTrkTRDOut[mult][part]=new TH1F(Form("hDenMatch_MultTrk%i_Part%iTRDOut",mult,part),"",46,fBinLim0);
      hDenMatchMultTrkTRDOut[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenMatchMultTrkTRDOut[mult][part]);
      hDenTrkMultTrkTRDOut[mult][part]=new TH1F(Form("hDenTrk_MultTrk%i_Part%iTRDOut",mult,part),"",46,fBinLim0);
      hDenTrkMultTrkTRDOut[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenTrkMultTrkTRDOut[mult][part]);
    }
  }

     for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultTrkNoTRDOut[mult][part]=new TH1F(Form("hNumMatch_MultTrk%i_Part%iNoTRDOut",mult,part),"",46,fBinLim0);
      hNumMatchMultTrkNoTRDOut[mult][part]->Sumw2();
      TOFCheck->AddLast(hNumMatchMultTrkNoTRDOut[mult][part]);
      hDenMatchMultTrkNoTRDOut[mult][part]=new TH1F(Form("hDenMatch_MultTrk%i_Part%iNoTRDOut",mult,part),"",46,fBinLim0);
      hDenMatchMultTrkNoTRDOut[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenMatchMultTrkNoTRDOut[mult][part]);
      hDenTrkMultTrkNoTRDOut[mult][part]=new TH1F(Form("hDenTrk_MultTrk%i_Part%iNoTRDOut",mult,part),"",46,fBinLim0);
      hDenTrkMultTrkNoTRDOut[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenTrkMultTrkNoTRDOut[mult][part]);
    }
  }

  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<2;part++){
      hNumMatchMultSPDInc[mult][part]=new TH1F(Form("hNumMatch_Inc_MultSPD%i_Charge%i",mult,part),"",46,fBinLim0);
      hNumMatchMultSPDInc[mult][part]->Sumw2();
      TOFCheck->AddLast(hNumMatchMultSPDInc[mult][part]);
      hDenMatchMultSPDInc[mult][part]=new TH1F(Form("hDenMatch_Inc_MultSPD%i_Charge%i",mult,part),"",46,fBinLim0);
      hDenMatchMultSPDInc[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenMatchMultSPDInc[mult][part]);
    }
  }

  for(Int_t mult=0; mult<7;mult++){
    for(Int_t part=0; part<6;part++){
      hNumMatchMultSPD[mult][part]=new TH1F(Form("hNumMatch_MultSPD%i_Part%i",mult,part),"",46,fBinLim0);
      hNumMatchMultSPD[mult][part]->Sumw2();
      TOFCheck->AddLast(hNumMatchMultSPD[mult][part]);
      hDenMatchMultSPD[mult][part]=new TH1F(Form("hDenMatch_MultSPD%i_Part%i",mult,part),"",46,fBinLim0);
      hDenMatchMultSPD[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenMatchMultSPD[mult][part]);
      hDenTrkMultSPD[mult][part]=new TH1F(Form("hDenTrk_MultSPD%i_Part%i",mult,part),"",46,fBinLim0);
      hDenTrkMultSPD[mult][part]->Sumw2();
      TOFCheck->AddLast(hDenTrkMultSPD[mult][part]);
    }
  }

  for(Int_t part=0; part<6;part++){
    hDenTrkVertMultTrk[part]=new TH1F(Form("hDenTrkVert_Part%i",part),"",46,fBinLim0);
    hDenTrkVertMultTrk[part]->Sumw2();
    TOFCheck->AddLast(hDenTrkVertMultTrk[part]);

  }
  
  for(Int_t part=0; part<6;part++){
    hDenTrkTriggerMultTrk[part]=new TH1F(Form("hDenTrkTrigger_Part%i",part),"",46,fBinLim0);
    hDenTrkTriggerMultTrk[part]->Sumw2();
    TOFCheck->AddLast(hDenTrkTriggerMultTrk[part]);
  }
  
  PostData(1, TOFCheck);

}

//________________________________________________________________________
void AliAnalysisCombinedHadronSpectra2MC::UserExec(Option_t *) 
{
  //
  // main event loop
  //
 
  
  multiplicity=-999, vert=-999, XPrimVertex=-999, YPrimVertex=-999, ZPrimVertex=-999, NContrPrimVertex=-999, fMCtracks=-999, fMCPrimaries=-999, spdCorr=-1.0;
  fT0TOF0=-999, fT0TOF1=-999, fT0TOF2=-999, fT0TOF3=-999, fT0TOF4=-999, fT0TOF5=-999, fT0TOF6=-999, fT0TOF7=-999, frun=-999;
  
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  if (!fESDtrackCuts) {
    Printf("ERROR: fESDtrackCuts not available");
    return;
  }
  
  hNumEv->Fill(1);
  
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    Printf("ERROR: Could not retrieve MC event handler");
    return;
  }
  
  AliESDInputHandler* esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (esdH)
    fESDpid = esdH->GetESDpid();
  
  AliMCEvent* mcEvent = 0x0;
  AliStack* stack = 0x0;
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  
  stack = mcEvent->Stack();
  if (!stack) return;

  //
  // check if event is selected by physics selection class
  //

  //trigger efficiency correction
  //loop on primary MC tracks
  for(Int_t i = 0; i < stack->GetNtrack(); i++) {
    treeMCPBis=-999; treeMCPtBis=-999; treeMCEtaBis=-999; treeMCPhiBis=-999; treeMCPdgBis=-999;    
    if (!stack->IsPhysicalPrimary(i)) continue;
    TParticle * trackMC = stack->Particle(i);
    //Double_t rapidityMC=-999;
    
    treeMCPBis=trackMC->P();
    treeMCPtBis=trackMC->Pt();
    treeMCEtaBis=trackMC->Eta();
    //if(TMath::Abs(treeMCEtaBis)>=0.9){continue;}
    treeMCPhiBis=trackMC->Phi()* 180 / TMath::Pi();
    treeMCPdgBis = trackMC->GetPdgCode();
    if(TMath::Abs(trackMC->Y())>=0.5) {continue;}
    if((TMath::Abs(treeMCPdgBis)!=211)&&(TMath::Abs(treeMCPdgBis)!=321)&&(TMath::Abs(treeMCPdgBis)!=2212)){continue;}
  
    Int_t PartTypeMC=-5;
    if((TMath::Abs(treeMCPdgBis)==211)&&(treeMCPdgBis>0)){PartTypeMC=0;}
    if((TMath::Abs(treeMCPdgBis)==211)&&(treeMCPdgBis<0)){PartTypeMC=1;}
    if((TMath::Abs(treeMCPdgBis)==321)&&(treeMCPdgBis>0)){PartTypeMC=2;}
    if((TMath::Abs(treeMCPdgBis)==321)&&(treeMCPdgBis<0)){PartTypeMC=3;}
    if((TMath::Abs(treeMCPdgBis)==2212)&&(treeMCPdgBis>0)){PartTypeMC=4;}
    if((TMath::Abs(treeMCPdgBis)==2212)&&(treeMCPdgBis<0)){PartTypeMC=5;}
  

    //fill Histo
    hDenTrkTriggerMultTrk[PartTypeMC]->Fill(treeMCPtBis);    

  }  



  Bool_t isSelected = kFALSE;
  isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()& AliVEvent::kMB);

  if (!isSelected) {
    return;
  }
  
  hNumEv->Fill(2);
  
  
  
  //vertex efficiency correction+senza taglio in eta
  //loop on primary MC tracks
  for(Int_t i = 0; i < stack->GetNtrack(); i++) {
    treeMCPBis=-999; treeMCPtBis=-999; treeMCEtaBis=-999; treeMCPhiBis=-999; treeMCPdgBis=-999;    
    if (!stack->IsPhysicalPrimary(i)) continue;
    TParticle * trackMC = stack->Particle(i);
    //Double_t rapidityMC=-999;
    
    treeMCPBis=trackMC->P();
    treeMCPtBis=trackMC->Pt();
    treeMCEtaBis=trackMC->Eta();
    //if(TMath::Abs(treeMCEtaBis)>=0.9){continue;}
    treeMCPhiBis=trackMC->Phi()* 180 / TMath::Pi();
    treeMCPdgBis = trackMC->GetPdgCode();
    if(TMath::Abs(trackMC->Y())>=0.5) {continue;}
    if((TMath::Abs(treeMCPdgBis)!=211)&&(TMath::Abs(treeMCPdgBis)!=321)&&(TMath::Abs(treeMCPdgBis)!=2212)){continue;}
  
    Int_t PartTypeMC=-5;
    if((TMath::Abs(treeMCPdgBis)==211)&&(treeMCPdgBis>0)){PartTypeMC=0;}
    if((TMath::Abs(treeMCPdgBis)==211)&&(treeMCPdgBis<0)){PartTypeMC=1;}
    if((TMath::Abs(treeMCPdgBis)==321)&&(treeMCPdgBis>0)){PartTypeMC=2;}
    if((TMath::Abs(treeMCPdgBis)==321)&&(treeMCPdgBis<0)){PartTypeMC=3;}
    if((TMath::Abs(treeMCPdgBis)==2212)&&(treeMCPdgBis>0)){PartTypeMC=4;}
    if((TMath::Abs(treeMCPdgBis)==2212)&&(treeMCPdgBis<0)){PartTypeMC=5;}
  

    //fill Histo
    hDenTrkVertMultTrk[PartTypeMC]->Fill(treeMCPtBis);    

  }     
       
 


  //
  // monitor vertex position
  //
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks(); //! Primary vertex estimated using ESD tracks
  if(vertex->GetNContributors()<1) { // # of tracklets/tracks used for the estimate
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD(); //! Primary vertex estimated by the SPD
    if(vertex->GetNContributors()<1) vertex = 0x0;
  }  

  if (!vertex) {
    vert=0;
  } else {vert=1;}

  if (!vertex) {
    return;
  }

  hNumEv->Fill(3);

  
  if (vertex) {ZPrimVertex=vertex->GetZv();}
  if(TMath::Abs(ZPrimVertex)>10){return;}

  hNumEv->Fill(4);
  multiplicity=Mult();

  fTreeEv->Fill();


  if (vertex) {XPrimVertex=vertex->GetXv(); YPrimVertex=vertex->GetYv(); ZPrimVertex=vertex->GetZv(); NContrPrimVertex=vertex->GetNContributors(); 
    vert=1;}  

  //multiplicity as defined by Marek
  const AliMultiplicity *mult = fESD->GetMultiplicity();
  Float_t nClusters[6]={0.0,0.0,0.0,0.0,0.0,0.0};
  for(Int_t ilay=0; ilay<6; ilay++)
    {
      nClusters[ilay] = (Float_t)mult->GetNumberOfITSClusters(ilay);
    }
  //cambio
  spdCorr = AliESDUtils::GetCorrSPD2(nClusters[1],vertex->GetZ());
  //spdCorr=50;
  // end cambio
  
  Int_t imult=-5;
  if((multiplicity>=0)&&(multiplicity<=5)){imult=1;}
  if((multiplicity>=6)&&(multiplicity<=9)){imult=2;}
  if((multiplicity>=10)&&(multiplicity<=14)){imult=3;}
  if((multiplicity>=15)&&(multiplicity<=22)){imult=4;}
  if((multiplicity>=23)&&(multiplicity<=32)){imult=5;}
  if(multiplicity>=33){imult=6;}
  
  Int_t imultSPD=-5;
  if((spdCorr>=0)&&(spdCorr<=16)){imultSPD=1;}
  if((spdCorr>=17)&&(spdCorr<=30)){imultSPD=2;}
  if((spdCorr>=31)&&(spdCorr<=45)){imultSPD=3;}
  if((spdCorr>=46)&&(spdCorr<=68)){imultSPD=4;}
  if((spdCorr>=69)&&(spdCorr<=97)){imultSPD=5;}
  if(spdCorr>=98){imult=6;}

  fMCtracks=mcEvent->GetNumberOfTracks();
  fMCPrimaries=mcEvent->GetNumberOfPrimaries();
  
  
  //cout<<"trk  "<<fMCtracks<<endl;
  //cout<<"trk primaries  "<<fMCPrimaries<<endl;
  
  //  //TOF settings done in the TOF tender
  
  //   frun = fESD->GetRunNumber();
  //   if(frun==frunOld){fLoadOCDB=kFALSE;}else {fLoadOCDB=kTRUE;}
  //   //if (tuneTOFMC) calibrateESD = kFALSE;
  //   Double_t *T0TOF;
  //   if(fLoadOCDB){
  //     AliCDBManager *cdb = AliCDBManager::Instance();
  //     cdb->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
  //     //cdb->SetDefaultStorage("raw://");
  //     cdb->SetRun(frun);}
  
  //   /* init TOF calibration */
  //   if (correctTExp)
  //     tofCalib->SetCorrectTExp(kTRUE);
  
  //   tofCalib->Init(frun);
  
  //   /* init TOF T0-maker */
  //   t0maker->SetTimeResolution(timeResolution);
  
  //   /* calibrate ESD */
  //   if (calibrateESD)
  //     tofCalib->CalibrateESD(fESD);
  
  //   T0TOF=t0maker->ComputeT0TOF(fESD);// calcola il t0 solo col tof in 10 bin di pt e lo setta in AliTOFPidResponse 
  
  //   //scrive i valori precedenti nel TOFHeader
  //   t0maker->WriteInESD(fESD); 
  
  //   //setta T0_TOF in AliTOFPidResponse ovvero i valori settati in AliTOFHeader. Se non ci sono setta T0spread
  //   fTimeZeroType=AliESDpid::kTOF_T0;
  
  //   fESDpid->SetTOFResponse(fESD,(AliESDpid::EStartTimeType_t)fTimeZeroType);
  
  
  //   fESDpid->MakePID(fESD,kFALSE); //calcola la sigma e le gi in più definisce la flag kTOFmismatch
  

  
  

//   fT0TOF0=T0TOF[0];
//   fT0TOF1=T0TOF[1];
//   fT0TOF2=T0TOF[2];
//   fT0TOF3=T0TOF[3];
//   fT0TOF4=T0TOF[4];
//   fT0TOF5=T0TOF[5];
//   fT0TOF6=T0TOF[6];
//   fT0TOF7=T0TOF[7];
 
  //
  // track loop
  //
  
  
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    
    AliESDtrack *track =fESD->GetTrack(i); 
    if (!track){continue;}
    // start TOF analysis
    
    rapidityMC=-999, fDCAXY=-999, fDCAZ=-999, fcut=-999, fTOFout=-999, ftrdout=-999, ftime=-999, ftpcclust=-999, flength=-999, fsign=-999, ftimetof=-999, ftofchan=-999, feta=-999, fphi=-999, fmomtrasv=-999,sigmapi=-999, sigmaka=-999, sigmapr=-999, fTot=-999, fmom=-999, fexptimepi=-999, fexptimeka=-999, fexptimepr=-999, ftofz=-999,ftofx=-999, TPCSignal=-999, TPCSigmaPI=-999, TPCSigmaKA=-999, TPCSigmaPR=-999, r1[0]=-999,r1[1]=-999,r1[2]=-999,r1[3]=-999,r1[4]=-999, fmatch=-999, fXout=-999, fYout=-999, fZout=-999, fPhiout=-999;
    

    if (!fESDtrackCuts->AcceptTrack(track)) {continue;}
    if (!(track->GetStatus()&AliESDtrack::kTOFout)==0) {fTOFout=1;}else {fTOFout=0;}
    if (!(track->GetStatus()&AliESDtrack::kTIME)==0) {ftime=1;}else {ftime=0;}  
    if (!(track->GetStatus()&AliESDtrack::kTRDout)==0) {ftrdout=1;}else {ftrdout=0;}  
    if (!(track->GetStatus()&AliESDtrack::kTOFmismatch)==0) {fmatch=0;}else {fmatch=1;}  
    track->GetImpactParameters(fDCAXY, fDCAZ);

    track->GetTOFpid(r1);
    fmomtrasv=track->Pt(); 
    feta=track->Eta();// return pseudorapidity return -TMath::Log(TMath::Tan(0.5 * Theta())); 
    if(!(TMath::Abs(feta)<0.9))continue;
    
    fmom=track->GetP();  // This function returns the track momentum
    Double_t *trackT0;
    // trackT0 = t0maker->GetT0p(fmom);// [0]=to -- [1] = sigma T0
//     fT0meas=trackT0[0]; 
//     fT0sigma=trackT0[1];
    ftpcclust=track->GetNcls(1);
    flength=track->GetIntegratedLength();
    fsign=track->GetSign();
    ftimetof=track->GetTOFsignal();
    Double_t inttime[5]; 
    track->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
    fexptimepi=inttime[2];
    fexptimeka=inttime[3];
    fexptimepr=inttime[4];

    fTot = track->GetTOFsignalToT();
    ftofz=track->GetTOFsignalDz(); // local z  of track's impact on the TOF pad 
    ftofx=track->GetTOFsignalDx(); // local x  of track's impact on the TOF pad 
    ftofchan=track->GetTOFCalChannel(); // Channel Index of the TOF Signal
    
    fphi=track->Phi()* 180 / TMath::Pi();// Returns the azimuthal angle of momentum  0 <= phi < 2*pi
    Double_t sigma[5];
    t0track = fESDpid->GetTOFResponse().GetStartTime(fmom); // T0best time
    t0trackSigma = fESDpid->GetTOFResponse().GetStartTimeRes(fmom); // T0best time

    for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++){
      sigma[ipart] = fESDpid->GetTOFResponse().GetExpectedSigma(fmom, inttime[ipart], AliPID::ParticleMass(ipart));}
    sigmapi=sigma[2];
    sigmaka=sigma[3];
    sigmapr=sigma[4];
    TPCSignal = track->GetTPCsignal();
    TPCSigmaPI=fESDpid->NumberOfSigmasTPC(track,AliPID::kPion);
    TPCSigmaKA=fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon);
    TPCSigmaPR=fESDpid->NumberOfSigmasTPC(track,AliPID::kProton);

    AliExternalTrackParam *exttrack=(AliExternalTrackParam *)track->GetOuterParam();
    if(exttrack){
    fPhiout=exttrack->Phi();
    fXout=exttrack->GetX();
    fYout=exttrack->GetY();
    fZout=exttrack->GetZ();}

    //start to take MC data
    Int_t flab=0;
    flab=track->GetLabel(); /*The Get*Label() getters return the label of the associated MC particle. The absolute value of this label is the index of the particle within the MC stack. If the label is negative, this track was assigned a certain number of clusters that did not in fact belong to this track. */
    Int_t abslab=TMath::Abs(flab);
    AliMCParticle *MCpart =(AliMCParticle* ) mcEvent->GetTrack(abslab);
    if(MCpart){
      TParticle *part = MCpart->Particle();
      if(part){
        fptMC=part->Pt();
	fphiMC=part->Phi()* 180 / TMath::Pi(); //angolo tra 0 e 2pi
	fetaMC=part->Eta();
	fPdgcode = part->GetPdgCode();
      }
    }

    //inizio plot match eff
    //numeratore
    
    if((TMath::Abs(TPCSigmaPI)<3)&&(TMath::Abs(TPCSigmaKA)>3)&&(TMath::Abs(TPCSigmaPR)>3)&&(fsign>0)){hDenMatchTPCpip->Fill(fmomtrasv);}
    if((TMath::Abs(TPCSigmaPI)>3)&&(TMath::Abs(TPCSigmaKA)<3)&&(TMath::Abs(TPCSigmaPR)>3)&&(fsign>0)){hDenMatchTPCkap->Fill(fmomtrasv);}
    if((TMath::Abs(TPCSigmaPI)>3)&&(TMath::Abs(TPCSigmaKA)>3)&&(TMath::Abs(TPCSigmaPR)<3)&&(fsign>0)){hDenMatchTPCprp->Fill(fmomtrasv);}
    if((TMath::Abs(TPCSigmaPI)<3)&&(TMath::Abs(TPCSigmaKA)>3)&&(TMath::Abs(TPCSigmaPR)>3)&&(fsign<0)){hDenMatchTPCpim->Fill(fmomtrasv);}
    if((TMath::Abs(TPCSigmaPI)>3)&&(TMath::Abs(TPCSigmaKA)<3)&&(TMath::Abs(TPCSigmaPR)>3)&&(fsign<0)){hDenMatchTPCkam->Fill(fmomtrasv);}
    if((TMath::Abs(TPCSigmaPI)>3)&&(TMath::Abs(TPCSigmaKA)>3)&&(TMath::Abs(TPCSigmaPR)<3)&&(fsign<0)){hDenMatchTPCprm->Fill(fmomtrasv);}

    hDenMatch->Fill(fmomtrasv);
    if(fsign>0) hDenMatchPos->Fill(fmomtrasv);
    if(fsign<0) hDenMatchNeg->Fill(fmomtrasv);
    hDenMatchEta->Fill(feta);
    if(fsign>0) hDenMatchPosEta->Fill(feta);
    if(fsign<0) hDenMatchNegEta->Fill(feta);

    if((TMath::Abs(fPdgcode)==211)&&(fsign>0)){hDenMatchPip->Fill(fmomtrasv);}
    if((TMath::Abs(fPdgcode)==211)&&(fsign<0)){hDenMatchPim->Fill(fmomtrasv);}
    if((TMath::Abs(fPdgcode)==321)&&(fsign>0)){hDenMatchKap->Fill(fmomtrasv);}
    if((TMath::Abs(fPdgcode)==321)&&(fsign<0)){hDenMatchKam->Fill(fmomtrasv);}
    if((TMath::Abs(fPdgcode)==2212)&&(fsign>0)){hDenMatchPrp->Fill(fmomtrasv);}
    if((TMath::Abs(fPdgcode)==2212)&&(fsign<0)){hDenMatchPrm->Fill(fmomtrasv);}

    if(fmomtrasv>0.5){
      hDenMatchEtaPtMa->Fill(feta);
      if(fsign>0) hDenMatchPosEtaPtMa->Fill(feta);
      if(fsign<0) hDenMatchNegEtaPtMa->Fill(feta);
      hDenMatchphiOutPtMa->Fill(fPhiout*180/TMath::Pi());
      if(fsign>0) hDenMatchPosphiOutPtMa->Fill(fPhiout*180/TMath::Pi());
      if(fsign<0) hDenMatchNegphiOutPtMa->Fill(fPhiout*180/TMath::Pi());
    }
    //Int_t NSM=tofchan/8736;
    hDenMatchphiOut->Fill(fPhiout*180/TMath::Pi());
    if(fsign>0) hDenMatchPosphiOut->Fill(fPhiout*180/TMath::Pi());
    if(fsign<0) hDenMatchNegphiOut->Fill(fPhiout*180/TMath::Pi());
    
    if(ftrdout==1){
      hDenMatchTRDOut->Fill(fmomtrasv);
      if(fsign>0) hDenMatchPosTRDOut->Fill(fmomtrasv);
      if(fsign<0) hDenMatchNegTRDOut->Fill(fmomtrasv);
    }

    if(ftrdout==0){
      hDenMatchNoTRDOut->Fill(fmomtrasv);
      if(fsign>0) hDenMatchPosNoTRDOut->Fill(fmomtrasv);
      if(fsign<0) hDenMatchNegNoTRDOut->Fill(fmomtrasv);
    }

    //denominatore
    if((fTOFout==1)&&(ftime==1)&&(flength>350)&&(ftimetof>10000)&&(fexptimepi>10000)&&(fexptimeka>10000)&&(fexptimepr>10000)&&(ftimetof<80000)){
      
      if((TMath::Abs(TPCSigmaPI)<3)&&(TMath::Abs(TPCSigmaKA)>3)&&(TMath::Abs(TPCSigmaPR)>3)&&(fsign>0)){hNumMatchTPCpip->Fill(fmomtrasv);}
      if((TMath::Abs(TPCSigmaPI)>3)&&(TMath::Abs(TPCSigmaKA)<3)&&(TMath::Abs(TPCSigmaPR)>3)&&(fsign>0)){hNumMatchTPCkap->Fill(fmomtrasv);}
      if((TMath::Abs(TPCSigmaPI)>3)&&(TMath::Abs(TPCSigmaKA)>3)&&(TMath::Abs(TPCSigmaPR)<3)&&(fsign>0)){hNumMatchTPCprp->Fill(fmomtrasv);}
      if((TMath::Abs(TPCSigmaPI)<3)&&(TMath::Abs(TPCSigmaKA)>3)&&(TMath::Abs(TPCSigmaPR)>3)&&(fsign<0)){hNumMatchTPCpim->Fill(fmomtrasv);}
      if((TMath::Abs(TPCSigmaPI)>3)&&(TMath::Abs(TPCSigmaKA)<3)&&(TMath::Abs(TPCSigmaPR)>3)&&(fsign<0)){hNumMatchTPCkam->Fill(fmomtrasv);}
      if((TMath::Abs(TPCSigmaPI)>3)&&(TMath::Abs(TPCSigmaKA)>3)&&(TMath::Abs(TPCSigmaPR)<3)&&(fsign<0)){hNumMatchTPCprm->Fill(fmomtrasv);}
      
      hNumMatch->Fill(fmomtrasv);
      if(fsign>0) hNumMatchPos->Fill(fmomtrasv);
      if(fsign<0) hNumMatchNeg->Fill(fmomtrasv);
      hNumMatchEta->Fill(feta);
      if(fsign>0) hNumMatchPosEta->Fill(feta);
      if(fsign<0) hNumMatchNegEta->Fill(feta);
      
      if((TMath::Abs(fPdgcode)==211)&&(fsign>0)){hNumMatchPip->Fill(fmomtrasv);}
      if((TMath::Abs(fPdgcode)==211)&&(fsign<0)){hNumMatchPim->Fill(fmomtrasv);}
      if((TMath::Abs(fPdgcode)==321)&&(fsign>0)){hNumMatchKap->Fill(fmomtrasv);}
      if((TMath::Abs(fPdgcode)==321)&&(fsign<0)){hNumMatchKam->Fill(fmomtrasv);}
      if((TMath::Abs(fPdgcode)==2212)&&(fsign>0)){hNumMatchPrp->Fill(fmomtrasv);}
      if((TMath::Abs(fPdgcode)==2212)&&(fsign<0)){hNumMatchPrm->Fill(fmomtrasv);}
      
      if(fmomtrasv>0.5){
	hNumMatchEtaPtMa->Fill(feta);
	if(fsign>0) hNumMatchPosEtaPtMa->Fill(feta);
	if(fsign<0) hNumMatchNegEtaPtMa->Fill(feta);
	hNumMatchphiOutPtMa->Fill(fPhiout*180/TMath::Pi());
	if(fsign>0) hNumMatchPosphiOutPtMa->Fill(fPhiout*180/TMath::Pi());
	if(fsign<0) hNumMatchNegphiOutPtMa->Fill(fPhiout*180/TMath::Pi());
      }
      hNumMatchphiOut->Fill(fPhiout*180/TMath::Pi());
      if(fsign>0) hNumMatchPosphiOut->Fill(fPhiout*180/TMath::Pi());
      if(fsign<0) hNumMatchNegphiOut->Fill(fPhiout*180/TMath::Pi());
      
      if(ftrdout==1){
	hNumMatchTRDOut->Fill(fmomtrasv);
	if(fsign>0) hNumMatchPosTRDOut->Fill(fmomtrasv);
	if(fsign<0) hNumMatchNegTRDOut->Fill(fmomtrasv);
      }
      
      if(ftrdout==0){
	hNumMatchNoTRDOut->Fill(fmomtrasv);
	if(fsign>0) hNumMatchPosNoTRDOut->Fill(fmomtrasv);
	if(fsign<0) hNumMatchNegNoTRDOut->Fill(fmomtrasv);
      }
    }
    
    //fine plot test matching efficiency

    pad->Fill(ftofx,ftofz);
    resx->Fill(ftofx);
    resz->Fill(ftofz);
    Float_t deltat;
    if((fmom>0.9)&&(fmom<1.1)){
      deltat=ftimetof-t0track-fexptimepi;
      tofres->Fill(deltat);
      if(t0track!=0){tofresTOF->Fill(deltat);}
      if((TMath::Abs(ftofx)<1.25)&&(TMath::Abs(ftofz)<1.75)){
	tofresgood->Fill(deltat);
      }
    }
    
    
    //inizio plot efficienze matching + tracking spettri
    Int_t ip=0;
    if(TMath::Abs(fPdgcode)==211){ip=2;}
    if(TMath::Abs(fPdgcode)==321){ip=3;}
    if(TMath::Abs(fPdgcode)==2212){ip=4;}
    if((ip!=2)&&(ip!=3)&&(ip!=4)){continue;}

    Float_t mass;
    mass=AliPID::ParticleMass(ip);//GeV
    Double_t momlung;
    momlung=TMath::Sqrt(fmom*fmom-fmomtrasv*fmomtrasv);
    Double_t transvmass;
    transvmass=TMath::Sqrt(mass*mass+fmomtrasv*fmomtrasv);
    rapidityMC=TMath::ASinH(momlung/transvmass);
    if(rapidityMC>=0.5){continue;}


    Int_t PartType=-5;
    if((TMath::Abs(fPdgcode)==211)&&(fsign>0)){PartType=0;}
    if((TMath::Abs(fPdgcode)==211)&&(fsign<0)){PartType=1;}
    if((TMath::Abs(fPdgcode)==321)&&(fsign>0)){PartType=2;}
    if((TMath::Abs(fPdgcode)==321)&&(fsign<0)){PartType=3;}
    if((TMath::Abs(fPdgcode)==2212)&&(fsign>0)){PartType=4;}
    if((TMath::Abs(fPdgcode)==2212)&&(fsign<0)){PartType=5;}

    Int_t Pos=-5;
    if((PartType==0)||(PartType==2)||(PartType==4)){Pos=0;}
    if((PartType==1)||(PartType==3)||(PartType==5)){Pos=1;}
     
    if((fTOFout==1)&&(ftime==1))
      {
	
	
	hNumMatchMultTrk[imult][PartType]->Fill(fmomtrasv);
        hNumMatchMultSPD[imultSPD][PartType]->Fill(fmomtrasv); 
	hNumMatchMultTrk[0][PartType]->Fill(fmomtrasv);
        hNumMatchMultSPD[0][PartType]->Fill(fmomtrasv); 
	
	hNumMatchMultTrkInc[imult][Pos]->Fill(fmomtrasv);
        hNumMatchMultSPDInc[imultSPD][Pos]->Fill(fmomtrasv); 
	hNumMatchMultTrkInc[0][Pos]->Fill(fmomtrasv);
        hNumMatchMultSPDInc[0][Pos]->Fill(fmomtrasv);
	
	if(ftrdout==1) {hNumMatchMultTrkTRDOut[0][PartType]->Fill(fmomtrasv);}
	if(ftrdout==0) {hNumMatchMultTrkNoTRDOut[0][PartType]->Fill(fmomtrasv);}
		
      }
    
    hDenMatchMultTrk[imult][PartType]->Fill(fmomtrasv);
    hDenMatchMultSPD[imultSPD][PartType]->Fill(fmomtrasv); 
    hDenMatchMultTrk[0][PartType]->Fill(fmomtrasv);
    hDenMatchMultSPD[0][PartType]->Fill(fmomtrasv); 
    
    hDenMatchMultTrkInc[imult][Pos]->Fill(fmomtrasv);
    hDenMatchMultSPDInc[imultSPD][Pos]->Fill(fmomtrasv); 
    hDenMatchMultTrkInc[0][Pos]->Fill(fmomtrasv);
    hDenMatchMultSPDInc[0][Pos]->Fill(fmomtrasv);

  
   
    

    fTreeTrack->Fill(); 
    
    /*
      kTOFout = TOF matching
      kTIME = good integrated time
      100000 > track->GetTOFsignal() > 12000 = TOF time reasanble range
      tracklength > 365 = should be greater than the TOF radius (370 cm)
    */
    
   
    
  } // end of track loop
  

   //secondo modo per trovare le primarie
  //loop on primary MC tracks
  for(Int_t i = 0; i < stack->GetNtrack(); i++) {
    treeMCPBis=-999; treeMCPtBis=-999; treeMCEtaBis=-999; treeMCPhiBis=-999; treeMCPdgBis=-999;    
    if (!stack->IsPhysicalPrimary(i)) continue;
    TParticle * trackMC = stack->Particle(i);
    //Double_t rapidityMC=-999;
    
    treeMCPBis=trackMC->P();
    treeMCPtBis=trackMC->Pt();
    treeMCEtaBis=trackMC->Eta();
    if(TMath::Abs(treeMCEtaBis)>=0.9){continue;}
    treeMCPhiBis=trackMC->Phi()* 180 / TMath::Pi();
    treeMCPdgBis = trackMC->GetPdgCode();
    if(TMath::Abs(trackMC->Y())>=0.5) {continue;}
    if((TMath::Abs(treeMCPdgBis)!=211)&&(TMath::Abs(treeMCPdgBis)!=321)&&(TMath::Abs(treeMCPdgBis)!=2212)){continue;}
  
    Int_t PartTypeMC=-5;
    if((TMath::Abs(treeMCPdgBis)==211)&&(treeMCPdgBis>0)){PartTypeMC=0;}
    if((TMath::Abs(treeMCPdgBis)==211)&&(treeMCPdgBis<0)){PartTypeMC=1;}
    if((TMath::Abs(treeMCPdgBis)==321)&&(treeMCPdgBis>0)){PartTypeMC=2;}
    if((TMath::Abs(treeMCPdgBis)==321)&&(treeMCPdgBis<0)){PartTypeMC=3;}
    if((TMath::Abs(treeMCPdgBis)==2212)&&(treeMCPdgBis>0)){PartTypeMC=4;}
    if((TMath::Abs(treeMCPdgBis)==2212)&&(treeMCPdgBis<0)){PartTypeMC=5;}
  

    //fill Histo
    hDenTrkMultTrk[imult][PartTypeMC]->Fill(treeMCPtBis);
    hDenTrkMultSPD[imultSPD][PartTypeMC]->Fill(treeMCPtBis);
    hDenTrkMultTrk[0][PartTypeMC]->Fill(treeMCPtBis);
    hDenTrkMultSPD[0][PartTypeMC]->Fill(treeMCPtBis); 
    

  }     
       
 


  frunOld=frun;
  

}      


//________________________________________________________________________
void AliAnalysisCombinedHadronSpectra2MC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf("*** CONSTRUCTOR CALLED ****");
   TOFCheck = dynamic_cast<TList*> (GetOutputData(1));
  if (!TOFCheck) {
    Printf("ERROR: TOFCheck not available");
    return;   
  } 
  
  system("touch ok.job");
}

