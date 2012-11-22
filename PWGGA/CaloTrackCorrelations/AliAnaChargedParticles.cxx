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

//_________________________________________________________________________
//
// Class for track selection and identification (not done now)
// Tracks from the CTS are kept in the AOD.
// Few histograms produced.
//
//-- Author: Gustavo Conesa (INFN-LNF)
//_________________________________________________________________________


// --- ROOT system ---
#include "TParticle.h"
#include "TH2F.h"

//---- AliRoot system ----
#include "AliAnaChargedParticles.h"
#include "AliCaloTrackReader.h"
#include "AliAODPWG4Particle.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"

ClassImp(AliAnaChargedParticles)
  
//__________________________________________________
  AliAnaChargedParticles::AliAnaChargedParticles() :
    AliAnaCaloTrackCorrBaseClass(),
    fPdg(0),           fFillPileUpHistograms(0),
    fhNtracks(0),      fhPt(0),
    fhPhiNeg(0),       fhEtaNeg(0), 
    fhPhiPos(0),       fhEtaPos(0), 
    fhEtaPhiPos(0),    fhEtaPhiNeg(0),
    //MC
    fhPtPion(0),       fhPhiPion(0),         fhEtaPion(0),
    fhPtProton(0),     fhPhiProton(0),       fhEtaProton(0),
    fhPtElectron(0),   fhPhiElectron(0),     fhEtaElectron(0),
    fhPtKaon(0),       fhPhiKaon(0),         fhEtaKaon(0),
    fhPtUnknown(0),    fhPhiUnknown(0),      fhEtaUnknown(0),
    fhTOFSignal(0),    fhTOFSignalPtCut(0),  fhTOFSignalBCOK(0),
    fhPtTOFSignal(0),  fhPtTOFStatus0(0),    fhEtaPhiTOFStatus0(0),
    fhEtaPhiTOFBC0(0), fhEtaPhiTOFBCPlus(0), fhEtaPhiTOFBCMinus(0),
    fhEtaPhiTOFBC0PileUpSPD(0),
    fhEtaPhiTOFBCPlusPileUpSPD(0),
    fhEtaPhiTOFBCMinusPileUpSPD(0)
{
  //Default Ctor

  for(Int_t i = 0; i < 7; i++)
  {
    fhPtPileUp         [i] = 0;
    fhPtTOFSignalPileUp[i] = 0;
  }
  
  for(Int_t i = 0; i < 3; i++)
  {
    fhPtDCA               [i] = 0 ;
    //fhPtDCAVtxOutBC0      [i] = 0 ;
    fhPtDCAPileUp         [i] = 0 ;
    //fhPtDCAVtxOutBC0PileUp[i] = 0 ;
    
    fhPtDCATOFBC0         [i] = 0 ;
    fhPtDCAPileUpTOFBC0   [i] = 0 ;

    fhPtDCANoTOFHit               [i] = 0 ;
    //fhPtDCAVtxOutBC0NoTOFHit      [i] = 0 ;
    fhPtDCAPileUpNoTOFHit         [i] = 0 ;
    //fhPtDCAVtxOutBC0PileUpNoTOFHit[i] = 0 ;
  }
  
  //Initialize parameters
  InitParameters();

}

//_______________________________________________________
TList *  AliAnaChargedParticles::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in fOutputContainer
  
  
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ExampleHistos") ; 
  
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins(); Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins(); Int_t netabins = GetHistogramRanges()->GetHistoEtaBins();
  Float_t ptmax  = GetHistogramRanges()->GetHistoPtMax();  Float_t phimax = GetHistogramRanges()->GetHistoPhiMax();  Float_t etamax = GetHistogramRanges()->GetHistoEtaMax();
  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();  Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();  Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();	

  fhNtracks  = new TH1F ("hNtracks","# of tracks", 1000,0,1000); 
  fhNtracks->SetXTitle("# of tracks");
  outputContainer->Add(fhNtracks);
    
  fhPt  = new TH1F ("hPt","p_T distribution", nptbins,ptmin,ptmax); 
  fhPt->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPt);
  
  fhPhiNeg  = new TH2F ("hPhiNegative","#phi of negative charges distribution",
                        nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
  fhPhiNeg->SetYTitle("#phi (rad)");
  fhPhiNeg->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPhiNeg);
  
  fhEtaNeg  = new TH2F ("hEtaNegative","#eta of negative charges distribution",
                        nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  fhEtaNeg->SetYTitle("#eta ");
  fhEtaNeg->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhEtaNeg);
  
  fhPhiPos  = new TH2F ("hPhiPositive","#phi of positive charges distribution",
                        nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
  fhPhiPos->SetYTitle("#phi (rad)");
  fhPhiPos->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPhiPos);
  
  fhEtaPos  = new TH2F ("hEtaPositive","#eta of positive charges distribution",
                        nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  fhEtaPos->SetYTitle("#eta ");
  fhEtaPos->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhEtaPos);
  
  fhEtaPhiPos  = new TH2F ("hEtaPhiPositive","pt/eta/phi of positive charge",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiPos->SetXTitle("#eta ");
  fhEtaPhiPos->SetYTitle("#phi (rad)");  
  outputContainer->Add(fhEtaPhiPos);
  
  fhEtaPhiNeg  = new TH2F ("hEtaPhiNegative","eta vs phi of negative charge",netabins,etamin,etamax, nphibins,phimin,phimax); 
  fhEtaPhiNeg->SetXTitle("#eta ");
  fhEtaPhiNeg->SetYTitle("#phi (rad)");  
  outputContainer->Add(fhEtaPhiNeg);
  
  fhProductionVertexBC      = new TH1F("hProductionVertexBC", "tracks production vertex bunch crossing ", 18 , -9 , 9 ) ;
  fhProductionVertexBC->SetYTitle("# tracks");
  fhProductionVertexBC->SetXTitle("Bunch crossing");
  outputContainer->Add(fhProductionVertexBC);

  Int_t ntofbins = 1000;
  Int_t mintof = -500;
  Int_t maxtof =  500;
  
  fhTOFSignal  = new TH1F ("hTOFSignal","TOF signal", ntofbins,mintof,maxtof);
  fhTOFSignal->SetXTitle("TOF signal (ns)");
  outputContainer->Add(fhTOFSignal);

  fhTOFSignalBCOK  = new TH1F ("hTOFSignalBCOK","TOF signal", ntofbins,mintof,maxtof);
  fhTOFSignalBCOK->SetXTitle("TOF signal (ns)");
  outputContainer->Add(fhTOFSignalBCOK);
  
  fhTOFSignalPtCut  = new TH1F ("hTOFSignalPtCut","TOF signal", ntofbins,mintof,maxtof);
  fhTOFSignalPtCut->SetXTitle("TOF signal (ns)");
  outputContainer->Add(fhTOFSignalPtCut);

  fhPtTOFSignal  = new TH2F ("hPtTOFSignal","TOF signal", nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
  fhPtTOFSignal->SetYTitle("TOF signal (ns)");
  fhPtTOFSignal->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtTOFSignal);

  if(fFillPileUpHistograms)
  {    
    TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;
    
    for(Int_t i = 0 ; i < 7 ; i++)
    {
      fhPtPileUp[i]  = new TH1F(Form("hPtPileUp%s",pileUpName[i].Data()),
                                Form("Track p_{T} distribution, %s Pile-Up event",pileUpName[i].Data()),
                                nptbins,ptmin,ptmax);
      fhPtPileUp[i]->SetXTitle("p_{T} (GeV/c)");
      outputContainer->Add(fhPtPileUp[i]);
            
      fhPtTOFSignalPileUp[i]  = new TH2F(Form("hPtTOFSignalPileUp%s",pileUpName[i].Data()),
                                         Form("Track TOF vs p_{T} distribution, %s Pile-Up event",pileUpName[i].Data()),
                                         nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
      fhPtTOFSignalPileUp[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtTOFSignalPileUp[i]->SetXTitle("TOF signal (ns)");
      outputContainer->Add(fhPtTOFSignalPileUp[i]);
    }
 
    fhEtaPhiTOFBC0  = new TH2F ("hEtaPhiTOFBC0","eta-phi for tracks with hit on TOF, and tof corresponding to BC=0",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiTOFBC0->SetXTitle("#eta ");
    fhEtaPhiTOFBC0->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhiTOFBC0);
    
    fhEtaPhiTOFBCPlus  = new TH2F ("hEtaPhiTOFBCPlus","eta-phi for tracks with hit on TOF, and tof corresponding to BC>0",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiTOFBCPlus->SetXTitle("#eta ");
    fhEtaPhiTOFBCPlus->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhiTOFBCPlus);
    
    fhEtaPhiTOFBCMinus  = new TH2F ("hEtaPhiTOFBCMinus","eta-phi for tracks with hit on TOF, and tof corresponding to BC<0",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiTOFBCMinus->SetXTitle("#eta ");
    fhEtaPhiTOFBCMinus->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhiTOFBCMinus);

    fhEtaPhiTOFBC0PileUpSPD  = new TH2F ("hEtaPhiTOFBC0PileUpSPD","eta-phi for tracks with hit on TOF, and tof corresponding to BC=0, SPD pile-up",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiTOFBC0PileUpSPD->SetXTitle("#eta ");
    fhEtaPhiTOFBC0PileUpSPD->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhiTOFBC0PileUpSPD);
    
    fhEtaPhiTOFBCPlusPileUpSPD  = new TH2F ("hEtaPhiTOFBCPlusPileUpSPD","eta-phi for tracks with hit on TOF, and tof corresponding to BC>0, SPD pile-up",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiTOFBCPlusPileUpSPD->SetXTitle("#eta ");
    fhEtaPhiTOFBCPlusPileUpSPD->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhiTOFBCPlusPileUpSPD);
    
    fhEtaPhiTOFBCMinusPileUpSPD  = new TH2F ("hEtaPhiTOFBCMinusPileUpSPD","eta-phi for tracks with hit on TOF, and tof corresponding to BC<0, SPD pile-up",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiTOFBCMinusPileUpSPD->SetXTitle("#eta ");
    fhEtaPhiTOFBCMinusPileUpSPD->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhiTOFBCMinusPileUpSPD);
    
  }
  
  fhPtTOFStatus0  = new TH1F ("hPtTOFStatus0","p_T distribution of tracks not hitting TOF", nptbins,ptmin,ptmax);
  fhPtTOFStatus0->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtTOFStatus0);
  
  
  fhEtaPhiTOFStatus0  = new TH2F ("hEtaPhiTOFStatus0","eta-phi for tracks without hit on TOF",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiTOFStatus0->SetXTitle("#eta ");
  fhEtaPhiTOFStatus0->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiTOFStatus0);

  TString dcaName[] = {"xy","z","Cons"} ;
  Int_t ndcabins = 800;
  Int_t mindca = -4;
  Int_t maxdca =  4;
  
  for(Int_t i = 0 ; i < 3 ; i++)
  {
    
    fhPtDCA[i]  = new TH2F(Form("hPtDCA%s",dcaName[i].Data()),
                           Form("Track DCA%s vs p_{T} distribution",dcaName[i].Data()),
                           nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCA[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCA[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCA[i]);
    
    fhPtDCATOFBC0[i]  = new TH2F(Form("hPtDCA%sTOFBC0",dcaName[i].Data()),
                           Form("Track DCA%s vs p_{T} distribution",dcaName[i].Data()),
                           nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCATOFBC0[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCATOFBC0[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCATOFBC0[i]);

    fhPtDCANoTOFHit[i]  = new TH2F(Form("hPtDCA%sNoTOFHit",dcaName[i].Data()),
                           Form("Track (no TOF hit) DCA%s vs p_{T} distribution",dcaName[i].Data()),
                           nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCANoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCANoTOFHit[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCANoTOFHit[i]);

//    fhPtDCAVtxOutBC0[i]  = new TH2F(Form("hPtDCA%sVtxOutBC0",dcaName[i].Data()),
//                           Form("Track DCA%s vs p_{T} distribution, vertex with BC!=0",dcaName[i].Data()),
//                           nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
//    fhPtDCAVtxOutBC0[i]->SetXTitle("p_{T} (GeV/c)");
//    fhPtDCAVtxOutBC0[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
//    outputContainer->Add(fhPtDCAVtxOutBC0[i]);
//    
//    fhPtDCAVtxOutBC0NoTOFHit[i]  = new TH2F(Form("hPtDCA%sVtxOutBC0NoTOFHit",dcaName[i].Data()),
//                                   Form("Track (no TOF hit) DCA%s vs p_{T} distribution, vertex with BC!=0",dcaName[i].Data()),
//                                   nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
//    fhPtDCAVtxOutBC0NoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
//    fhPtDCAVtxOutBC0NoTOFHit[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
//    outputContainer->Add(fhPtDCAVtxOutBC0NoTOFHit[i]);
    
    if(fFillPileUpHistograms)
    {
      fhPtDCAPileUp[i]  = new TH2F(Form("hPtDCA%sPileUp",dcaName[i].Data()),
                             Form("Track DCA%s vs p_{T} distribution, SPD Pile-Up",dcaName[i].Data()),
                             nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAPileUp[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAPileUp[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAPileUp[i]);
      
      fhPtDCAPileUpTOFBC0[i]  = new TH2F(Form("hPtDCA%sPileUpTOFBC0",dcaName[i].Data()),
                                   Form("Track DCA%s vs p_{T} distribution",dcaName[i].Data()),
                                   nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAPileUpTOFBC0[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAPileUpTOFBC0[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAPileUpTOFBC0[i]);
      
      fhPtDCAPileUpNoTOFHit[i]  = new TH2F(Form("hPtDCA%sPileUpNoTOFHit",dcaName[i].Data()),
                                     Form("Track (no TOF hit) DCA%s vs p_{T} distribution, SPD Pile-Up, vertex with BC!=0",dcaName[i].Data()),
                                     nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAPileUpNoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAPileUpNoTOFHit[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAPileUpNoTOFHit[i]);
      
//      fhPtDCAVtxOutBC0PileUp[i]  = new TH2F(Form("hPtDCA%sPileUpVtxOutBC0",dcaName[i].Data()),
//                                   Form("Track DCA%s vs p_{T} distribution, SPD Pile-Up",dcaName[i].Data()),
//                                   nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
//      fhPtDCAVtxOutBC0PileUp[i]->SetXTitle("p_{T} (GeV/c)");
//      fhPtDCAVtxOutBC0PileUp[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
//      outputContainer->Add(fhPtDCAVtxOutBC0PileUp[i]);
//      
//      fhPtDCAVtxOutBC0PileUpNoTOFHit[i]  = new TH2F(Form("hPtDCA%sVtxOutBC0PileUpNoTOFHit",dcaName[i].Data()),
//                                           Form("Track (no TOF hit) DCA%s vs p_{T} distribution, SPD Pile-Up, vertex with BC!=0",dcaName[i].Data()),
//                                           nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
//      fhPtDCAVtxOutBC0PileUpNoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
//      fhPtDCAVtxOutBC0PileUpNoTOFHit[i]->SetXTitle(Form("DCA_{%s}",dcaName[i].Data()));
//      outputContainer->Add(fhPtDCAVtxOutBC0PileUpNoTOFHit[i]);
    }
  }

  
  
  if(IsDataMC()){
    
    fhPtPion  = new TH1F ("hPtMCPion","p_T distribution from #pi", nptbins,ptmin,ptmax); 
    fhPtPion->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtPion);
    
    fhPhiPion  = new TH2F ("hPhiMCPion","#phi distribution from #pi",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiPion->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiPion);
    
    fhEtaPion  = new TH2F ("hEtaMCPion","#eta distribution from #pi",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaPion->SetXTitle("#eta ");
    outputContainer->Add(fhEtaPion);
    
    fhPtProton  = new TH1F ("hPtMCProton","p_T distribution from proton", nptbins,ptmin,ptmax); 
    fhPtProton->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtProton);
    
    fhPhiProton  = new TH2F ("hPhiMCProton","#phi distribution from proton",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiProton->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiProton);
    
    fhEtaProton  = new TH2F ("hEtaMCProton","#eta distribution from proton",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaProton->SetXTitle("#eta ");
    outputContainer->Add(fhEtaProton);
    
    fhPtKaon  = new TH1F ("hPtMCKaon","p_T distribution from kaon", nptbins,ptmin,ptmax); 
    fhPtKaon->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtKaon);
    
    fhPhiKaon  = new TH2F ("hPhiMCKaon","#phi distribution from kaon",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiKaon->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiKaon);
    
    fhEtaKaon  = new TH2F ("hEtaMCKaon","#eta distribution from kaon",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaKaon->SetXTitle("#eta ");
    outputContainer->Add(fhEtaKaon);
    
    fhPtElectron  = new TH1F ("hPtMCElectron","p_T distribution from electron", nptbins,ptmin,ptmax); 
    fhPtElectron->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtElectron);
    
    fhPhiElectron  = new TH2F ("hPhiMCElectron","#phi distribution from electron",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiElectron->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiElectron);
    
    fhEtaElectron  = new TH2F ("hEtaMCElectron","#eta distribution from electron",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaElectron->SetXTitle("#eta ");
    outputContainer->Add(fhEtaElectron);
    
    fhPtUnknown  = new TH1F ("hPtMCUnknown","p_T distribution from unknown", nptbins,ptmin,ptmax); 
    fhPtUnknown->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtUnknown);
    
    fhPhiUnknown  = new TH2F ("hPhiMCUnknown","#phi distribution from unknown",nptbins,ptmin,ptmax, nphibins,phimin,phimax); 
    fhPhiUnknown->SetXTitle("#phi (rad)");
    outputContainer->Add(fhPhiUnknown);
    
    fhEtaUnknown  = new TH2F ("hEtaMCUnknown","#eta distribution from unknown",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
    fhEtaUnknown->SetXTitle("#eta ");
    outputContainer->Add(fhEtaUnknown);
    
  }
  
  return outputContainer;

}

//___________________________________________
void AliAnaChargedParticles::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  SetOutputAODClassName("AliAODPWG4Particle");
  SetOutputAODName("PWG4Particle");

  AddToHistogramsName("AnaCharged_");

  fPdg = -1; //Select all tracks 
  
}

//____________________________________________________________
void AliAnaChargedParticles::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");	
	
  printf("Min Pt = %3.2f\n", GetMinPt());
  printf("Max Pt = %3.2f\n", GetMaxPt());
  printf("Select clusters with pdg %d \n",fPdg);
  
} 

//_________________________________
void AliAnaChargedParticles::Init()
{  
  //Init
  //Do some checks
  if(!GetReader()->IsCTSSwitchedOn()){
    printf("AliAnaChargedParticles::Init() - STOP!: You want to use CTS tracks in analysis but not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
}

//_________________________________________________
void  AliAnaChargedParticles::MakeAnalysisFillAOD() 
{
  //Do analysis and fill aods
  if(!GetCTSTracks() || GetCTSTracks()->GetEntriesFast() == 0) return ;
  
  Int_t ntracks = GetCTSTracks()->GetEntriesFast();
  Double_t vert[3] = {0,0,0}; //vertex ;
  
  //Some prints
  if(GetDebug() > 0)
    printf("AliAnaChargedParticles::MakeAnalysisFillAOD() - In CTS aod entries %d\n", ntracks);
  
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (GetReader()->GetInputEvent());
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (GetReader()->GetInputEvent());
  
  Double_t bz = GetReader()->GetInputEvent()->GetMagneticField();

  //Fill AODParticle with CTS aods
  TVector3 p3;
  Int_t evtIndex = 0;
  for(Int_t i = 0; i < ntracks; i++){
    
    AliVTrack * track =  (AliVTrack*) (GetCTSTracks()->At(i));
    
     //Fill AODParticle after some selection
    Double_t mom[3] = {track->Px(),track->Py(),track->Pz()};
    p3.SetXYZ(mom[0],mom[1],mom[2]);
    
    //TOF
    ULong_t status = track->GetStatus();
    Bool_t okTOF = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
    Double32_t tof = track->GetTOFsignal()*1e-3;
    //if( tof < 0) printf("TOF Signal %e, status %d, pt %f\n", tof,status,status2,p3.Pt());
    
    Double_t dcaCons = -999;
    Int_t vtxBC = -999;
    if     (esdevent) vtxBC = esdevent->GetPrimaryVertex()->GetBC();
    else if(aodevent) vtxBC = aodevent->GetPrimaryVertex()->GetBC();
    
    if(vtxBC!=AliVTrack::kTOFBCNA)printf("BC primary %d",vtxBC);
    
    AliAODTrack * aodTrack = dynamic_cast<AliAODTrack*>(track);
    if(aodTrack)
    {
      dcaCons = aodTrack->DCA();
      vtxBC   = aodTrack->GetProdVertex()->GetBC();
    }

    if(vtxBC!=AliVTrack::kTOFBCNA) printf(" - production %d\n",vtxBC);

    fhProductionVertexBC->Fill(vtxBC);
    
    Double_t dca[2]   = {1e6,1e6};
    Double_t covar[3] = {1e6,1e6,1e6};
    track->PropagateToDCA(GetReader()->GetInputEvent()->GetPrimaryVertex(),bz,100.,dca,covar);
    
    if(dcaCons == -999)
    {
      fhPtDCA[0]->Fill(p3.Pt(), dca[0]);
      fhPtDCA[1]->Fill(p3.Pt(), dca[1]);
    }
    else fhPtDCA[2]->Fill(p3.Pt(), dcaCons);
    
//    if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
//    {
//      if(dcaCons == -999)
//      {
//        fhPtDCAVtxOutBC0[0]->Fill(p3.Pt(), dca[0]);
//        fhPtDCAVtxOutBC0[1]->Fill(p3.Pt(), dca[1]);
//      }
//      else
//        fhPtDCAVtxOutBC0[2]->Fill(p3.Pt(), dcaCons);
//    }
    
    if(fFillPileUpHistograms && GetReader()->IsPileUpFromSPD())
    {
      if(dcaCons == -999)
      {
        fhPtDCAPileUp[0]->Fill(p3.Pt(), dca[0]);
        fhPtDCAPileUp[1]->Fill(p3.Pt(), dca[1]);
      }
      else
        fhPtDCAPileUp[2]->Fill(p3.Pt(), dcaCons);
      
//      if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
//      {
//        if(dcaCons == -999)
//        {
//          fhPtDCAVtxOutBC0PileUp[0]->Fill(p3.Pt(), dca[0]);
//          fhPtDCAVtxOutBC0PileUp[1]->Fill(p3.Pt(), dca[1]);
//        }
//        else fhPtDCAVtxOutBC0PileUp[2]->Fill(p3.Pt(), dcaCons);
//      }
    }

    if(!okTOF)
    {
      if(dcaCons == -999)
      {
        fhPtDCANoTOFHit[0]->Fill(p3.Pt(), dca[0]);
        fhPtDCANoTOFHit[1]->Fill(p3.Pt(), dca[1]);
      }
      else
        fhPtDCANoTOFHit[2]->Fill(p3.Pt(), dcaCons);
      
//      if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
//      {
//        if(dcaCons == -999)
//        {
//          fhPtDCAVtxOutBC0NoTOFHit[0]->Fill(p3.Pt(), dca[0]);
//          fhPtDCAVtxOutBC0NoTOFHit[1]->Fill(p3.Pt(), dca[1]);
//        }
//        else
//          fhPtDCAVtxOutBC0NoTOFHit[2]->Fill(p3.Pt(), dcaCons);
//      }
      
      if(fFillPileUpHistograms && GetReader()->IsPileUpFromSPD())
      {
        if(dcaCons == -999)
        {
          fhPtDCAPileUpNoTOFHit[0]->Fill(p3.Pt(), dca[0]);
          fhPtDCAPileUpNoTOFHit[1]->Fill(p3.Pt(), dca[1]);
        }
        else
          fhPtDCAPileUpNoTOFHit[2]->Fill(p3.Pt(), dcaCons);
        
//        if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
//        {
//          if(dcaCons == -999)
//          {
//            fhPtDCAVtxOutBC0PileUpNoTOFHit[0]->Fill(p3.Pt(), dca[0]);
//            fhPtDCAVtxOutBC0PileUpNoTOFHit[1]->Fill(p3.Pt(), dca[1]);
//          }
//          else
//            fhPtDCAVtxOutBC0PileUpNoTOFHit[2]->Fill(p3.Pt(), dcaCons);
//        }
      }
    }
    
    //printf("track pT %2.2f, DCA Cons %f, DCA1 %f, DCA2 %f, TOFBC %d, oktof %d, tof %f\n",
    //      p3.Pt(),dcaCons,dca[0],dca[1],track->GetTOFBunchCrossing(bz),okTOF, tof);
    
    Int_t trackBC = track->GetTOFBunchCrossing(bz);
    
    if(okTOF && trackBC==0)
    {
      fhTOFSignalBCOK->Fill(tof);
      
      if(dcaCons == -999)
      {
        fhPtDCATOFBC0[0]->Fill(p3.Pt(), dca[0]);
        fhPtDCATOFBC0[1]->Fill(p3.Pt(), dca[1]);
      }
      else
        fhPtDCATOFBC0[2]->Fill(p3.Pt(), dcaCons);
      
      if(fFillPileUpHistograms && GetReader()->IsPileUpFromSPD())
      {
        if(dcaCons == -999)
        {
          fhPtDCAPileUpTOFBC0[0]->Fill(p3.Pt(), dca[0]);
          fhPtDCAPileUpTOFBC0[1]->Fill(p3.Pt(), dca[1]);
        }
        else
          fhPtDCAPileUpTOFBC0[2]->Fill(p3.Pt(), dcaCons);
      }
    }
    
    if(okTOF && fFillPileUpHistograms)
    {
      fhTOFSignal  ->Fill(tof);
      fhPtTOFSignal->Fill(p3.Pt(), tof);
            
      if(GetReader()->IsPileUpFromSPD())               fhPtTOFSignalPileUp[0]->Fill(p3.Pt(), tof); 
      if(GetReader()->IsPileUpFromEMCal())             fhPtTOFSignalPileUp[1]->Fill(p3.Pt(), tof);
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtTOFSignalPileUp[2]->Fill(p3.Pt(), tof);
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtTOFSignalPileUp[3]->Fill(p3.Pt(), tof);
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtTOFSignalPileUp[4]->Fill(p3.Pt(), tof);
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtTOFSignalPileUp[5]->Fill(p3.Pt(), tof);
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtTOFSignalPileUp[6]->Fill(p3.Pt(), tof);
      
      if      (trackBC ==0)  { fhEtaPhiTOFBC0    ->Fill(track->Eta(),track->Phi()); if(GetReader()->IsPileUpFromSPD()) fhEtaPhiTOFBC0PileUpSPD    ->Fill(track->Eta(),track->Phi()); }
      else if (trackBC < 0)  { fhEtaPhiTOFBCPlus ->Fill(track->Eta(),track->Phi()); if(GetReader()->IsPileUpFromSPD()) fhEtaPhiTOFBCPlusPileUpSPD ->Fill(track->Eta(),track->Phi()); }
      else if (trackBC > 0)  { fhEtaPhiTOFBCMinus->Fill(track->Eta(),track->Phi()); if(GetReader()->IsPileUpFromSPD()) fhEtaPhiTOFBCMinusPileUpSPD->Fill(track->Eta(),track->Phi()); }

    }
        
    Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,"CTS") ;
    
    if(GetDebug() > 1) 
      printf("AliAnaChargedParticles::MakeAnalysisFillAOD() - Track pt %2.2f, phi %2.2f, eta %2.2f in fiducial cut %d\n",p3.Pt(), p3.Phi(), p3.Eta(),in);
    
    //Acceptance selection
    if(IsFiducialCutOn() && ! in ) continue ;
    
    // Momentum selection
    if(track->Pt() < GetMinPt() || track->Pt() > GetMaxPt()) continue;
    
    if(okTOF) fhTOFSignalPtCut->Fill(tof); 
    else
    {
      fhPtTOFStatus0    ->Fill(track->Pt());
      fhEtaPhiTOFStatus0->Fill(track->Eta(),track->Phi());
    }
    
    //Keep only particles identified with fPdg
    //Selection not done for the moment
    //Should be done here.
    
    // Mixed event
    if (GetMixedEvent())
    {
      evtIndex = GetMixedEvent()->EventIndex(track->GetID()) ;
    }
    
    GetVertex(vert,evtIndex); 
    if(TMath::Abs(vert[2])> GetZvertexCut()) return; 
        
    AliAODPWG4Particle tr = AliAODPWG4Particle(mom[0],mom[1],mom[2],0);
    tr.SetDetector("CTS");
    tr.SetLabel(track->GetLabel());
    tr.SetTrackLabel(track->GetID(),-1);
    tr.SetChargedBit(track->Charge()>0);
		
    AddAODParticle(tr);
    
  }//loop
  
  if(GetDebug() > 0) 	
    printf("AliAnaChargedParticles::MakeAnalysisFillAOD() - Final aod branch entries %d\n", GetOutputAODBranch()->GetEntriesFast());   
} 

//__________________________________________________________________
void  AliAnaChargedParticles::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  
  //Loop on stored AODParticles
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  
  fhNtracks->Fill(GetReader()->GetTrackMultiplicity()) ;
  
  if(GetDebug() > 0) 
    printf("AliAnaChargedParticles::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* tr =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
        
    fhPt->Fill(tr->Pt());
    
    if(tr->GetChargedBit()){
      fhPhiPos   ->Fill(tr->Pt(), tr->Phi());
      fhEtaPos   ->Fill(tr->Pt(), tr->Eta());
      fhEtaPhiPos->Fill(tr->Eta(),tr->Phi());
    }
    else{
      fhPhiNeg   ->Fill(tr->Pt(), tr->Phi());
      fhEtaNeg   ->Fill(tr->Pt(), tr->Eta());
      fhEtaPhiNeg->Fill(tr->Eta(),tr->Phi());
    }
    
    if(fFillPileUpHistograms)
    {
      if(GetReader()->IsPileUpFromSPD())               {fhPtPileUp[0]->Fill(tr->Pt());}
      if(GetReader()->IsPileUpFromEMCal())             {fhPtPileUp[1]->Fill(tr->Pt());}
      if(GetReader()->IsPileUpFromSPDOrEMCal())        {fhPtPileUp[2]->Fill(tr->Pt());}
      if(GetReader()->IsPileUpFromSPDAndEMCal())       {fhPtPileUp[3]->Fill(tr->Pt());}
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    {fhPtPileUp[4]->Fill(tr->Pt());}
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    {fhPtPileUp[5]->Fill(tr->Pt());}
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) {fhPtPileUp[6]->Fill(tr->Pt());}
    }

    
    if(IsDataMC())
    {
      //Play with the MC stack if available		
      Int_t mompdg = -1;
      Int_t label  = tr->GetLabel();
      if(label >= 0)
      {
        if( GetReader()->ReadStack() && label < GetMCStack()->GetNtrack())
        {
          TParticle * mom = GetMCStack()->Particle(label);
          mompdg =TMath::Abs(mom->GetPdgCode());
        }
        else if(GetReader()->ReadAODMCParticles())
        {
          AliAODMCParticle * aodmom = 0;
          //Get the list of MC particles
          aodmom = (AliAODMCParticle*) (GetReader()->GetAODMCParticles(tr->GetInputFileIndex()))->At(label);
          mompdg =TMath::Abs(aodmom->GetPdgCode());
        }
      }
      
      if(mompdg==211)
      {
        fhPtPion ->Fill(tr->Pt());
        fhPhiPion->Fill(tr->Pt(), tr->Phi());
        fhEtaPion->Fill(tr->Pt(), tr->Eta());
      }
      else if(mompdg==2212)
      {
        fhPtProton ->Fill(tr->Pt());
        fhPhiProton->Fill(tr->Pt(), tr->Phi());
        fhEtaProton->Fill(tr->Pt(), tr->Eta());
      }
      else if(mompdg==321)
      {
        fhPtKaon ->Fill(tr->Pt());
        fhPhiKaon->Fill(tr->Pt(), tr->Phi());
        fhEtaKaon->Fill(tr->Pt(), tr->Eta());
      }
      else if(mompdg==11)
      {
        fhPtElectron ->Fill(tr->Pt());
        fhPhiElectron->Fill(tr->Pt(), tr->Phi());
        fhEtaElectron->Fill(tr->Pt(), tr->Eta());
      }
      else {
        //printf("unknown pdg %d\n",mompdg);
        fhPtUnknown ->Fill(tr->Pt());
        fhPhiUnknown->Fill(tr->Pt(), tr->Phi());
        fhEtaUnknown->Fill(tr->Pt(), tr->Eta());
      }
      
    }//Work with stack also
    
  }// aod branch loop
  
}
