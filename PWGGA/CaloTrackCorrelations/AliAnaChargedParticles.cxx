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
    fFillPileUpHistograms(0),
    fFillVertexBC0Histograms(0),
    //Histograms
    fhNtracks(0),      fhPt(0),            fhPtNoCut(0),
    fhPtCutDCA(0),     fhPtCutDCABCOK(0),
    fhPhiNeg(0),       fhEtaNeg(0),
    fhPhiPos(0),       fhEtaPos(0), 
    fhEtaPhiPos(0),    fhEtaPhiNeg(0),
    fhPtVtxOutBC0(0),  fhEtaPhiVtxOutBC0(0),
    fhPtVtxInBC0(0),   fhEtaPhiVtxInBC0(0),
    fhPtSPDRefit(0),         fhPtNoSPDRefit(0),         fhPtNoSPDNoRefit(0),
    fhEtaPhiSPDRefitPt02(0), fhEtaPhiNoSPDRefitPt02(0), fhEtaPhiNoSPDNoRefitPt02(0),
    fhEtaPhiSPDRefitPt3(0),  fhEtaPhiNoSPDRefitPt3(0),  fhEtaPhiNoSPDNoRefitPt3(0),
    //MC
    fhPtPion(0),       fhPhiPion(0),         fhEtaPion(0),
    fhPtProton(0),     fhPhiProton(0),       fhEtaProton(0),
    fhPtElectron(0),   fhPhiElectron(0),     fhEtaElectron(0),
    fhPtKaon(0),       fhPhiKaon(0),         fhEtaKaon(0),
    fhPtUnknown(0),    fhPhiUnknown(0),      fhEtaUnknown(0),
    //TOF
    fhTOFSignal(0),    fhTOFSignalPtCut(0),  fhTOFSignalBCOK(0),
    fhPtTOFSignal(0),  fhPtTOFSignalDCACut(0),
    fhPtTOFSignalVtxOutBC0(0), fhPtTOFSignalVtxInBC0(0),
    fhPtTOFStatus0(0), fhEtaPhiTOFStatus0(0),
    fhEtaPhiTOFBC0(0), fhEtaPhiTOFBCPlus(0), fhEtaPhiTOFBCMinus(0),
    fhEtaPhiTOFBC0PileUpSPD(0),
    fhEtaPhiTOFBCPlusPileUpSPD(0),
    fhEtaPhiTOFBCMinusPileUpSPD(0),
    fhProductionVertexBC(0)
{
  //Default Ctor

  for(Int_t i = 0; i < 7; i++)
  {
    fhPtPileUp         [i] = 0;
    fhPtTOFSignalPileUp[i] = 0;
    fhPtTOFSignalVtxOutBC0PileUp[i] = 0;
    fhPtTOFSignalVtxInBC0PileUp [i] = 0;
    fhProductionVertexBCPileUp  [i] = 0;
  }
  
  for(Int_t i = 0; i < 3; i++)
  {
    fhPtDCA               [i] = 0 ;
    
    fhPtDCASPDRefit       [i] = 0 ;
    fhPtDCANoSPDRefit     [i] = 0 ;
    fhPtDCANoSPDNoRefit   [i] = 0 ;

    fhPtDCAPileUp         [i] = 0 ;
    fhPtDCATOFBC0         [i] = 0 ;
    fhPtDCATOFBCOut       [i] = 0 ;
    fhPtDCAPileUpTOFBC0   [i] = 0 ;
    fhPtDCANoTOFHit       [i] = 0 ;
    fhPtDCAPileUpNoTOFHit [i] = 0 ;

    fhPtDCAVtxOutBC0      [i] = 0 ;
    fhPtDCAVtxInBC0       [i] = 0 ;
    fhPtDCAVtxOutBC0PileUp[i] = 0 ;
    fhPtDCAVtxInBC0PileUp [i] = 0 ;
    fhPtDCAVtxOutBC0NoTOFHit[i] = 0 ;
    fhPtDCAVtxInBC0NoTOFHit [i] = 0 ;
    fhPtDCAVtxOutBC0PileUpNoTOFHit[i] = 0 ;
    fhPtDCAVtxInBC0PileUpNoTOFHit [i] = 0 ;
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
  
  fhPtNoCut  = new TH1F ("hPtNoCut","p_T distribution, raw tracks", nptbins,ptmin,ptmax);
  fhPtNoCut->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtNoCut);
  
  fhPtCutDCA  = new TH1F ("hPtCutDCA","p_T distribution, cut DCA", nptbins,ptmin,ptmax);
  fhPtCutDCA->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtCutDCA);
  
  fhPtCutDCABCOK  = new TH1F ("hPtCutDCABCOK","p_T distribution, DCA cut, track BC=0 or -100", nptbins,ptmin,ptmax);
  fhPtCutDCABCOK->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtCutDCABCOK);

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
  
  if(fFillVertexBC0Histograms)
  {
    fhPtVtxOutBC0  = new TH1F ("hPtVtxOutBC0","p_T distribution, vertex in BC=0", nptbins,ptmin,ptmax);
    fhPtVtxOutBC0->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtVtxOutBC0);
    
    fhEtaPhiVtxOutBC0  = new TH2F ("hEtaPhiVtxOutBC0","eta vs phi of all charges with vertex in BC=0",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiVtxOutBC0->SetXTitle("#eta ");
    fhEtaPhiVtxOutBC0->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhiVtxOutBC0);
    
    fhPtVtxInBC0  = new TH1F ("hPtVtxInBC0","p_T distribution, vertex in BC=0", nptbins,ptmin,ptmax);
    fhPtVtxInBC0->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtVtxInBC0);
    
    fhEtaPhiVtxInBC0  = new TH2F ("hEtaPhiVtxInBC0","eta vs phi of all charges with vertex in BC=0",netabins,etamin,etamax, nphibins,phimin,phimax);
    fhEtaPhiVtxInBC0->SetXTitle("#eta ");
    fhEtaPhiVtxInBC0->SetYTitle("#phi (rad)");
    outputContainer->Add(fhEtaPhiVtxInBC0);
  }
  
  fhPtSPDRefit  = new TH1F ("hPtSPDRefit","p_T distribution of tracks with SPD and ITS refit", nptbins,ptmin,ptmax);
  fhPtSPDRefit->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtSPDRefit);

  fhEtaPhiSPDRefitPt02  = new TH2F ("hEtaPhiSPDRefitPt02","eta vs phi of tracks with SPD and ITS refit, p_{T} < 2 GeV/c",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiSPDRefitPt02->SetXTitle("#eta ");
  fhEtaPhiSPDRefitPt02->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiSPDRefitPt02);
  
  fhEtaPhiSPDRefitPt3  = new TH2F ("hEtaPhiSPDRefitPt3","eta vs phi of tracks with SPD and ITS refit, p_{T} > 3 GeV/c",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiSPDRefitPt3->SetXTitle("#eta ");
  fhEtaPhiSPDRefitPt3->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiSPDRefitPt3);

  fhPtNoSPDRefit  = new TH1F ("hPtNoSPDRefit","p_T distribution of constrained tracks no SPD and with ITSRefit", nptbins,ptmin,ptmax);
  fhPtNoSPDRefit->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtNoSPDRefit);
  
  fhEtaPhiNoSPDRefitPt02  = new TH2F ("hEtaPhiNoSPDRefitPt02","eta vs phi of constrained tracks no SPD and with ITSRefit, p_{T} < 2 GeV/c",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiNoSPDRefitPt02->SetXTitle("#eta ");
  fhEtaPhiNoSPDRefitPt02->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiNoSPDRefitPt02);
  
  fhEtaPhiNoSPDRefitPt3  = new TH2F ("hEtaPhiNoSPDRefitPt3","eta vs phi of of constrained tracks no SPD and with ITSRefit, p_{T} > 3 GeV/c",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiNoSPDRefitPt3->SetXTitle("#eta ");
  fhEtaPhiNoSPDRefitPt3->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiNoSPDRefitPt3);
  
  fhPtNoSPDNoRefit  = new TH1F ("hPtNoSPDNoRefit","p_T distribution of constrained tracks with no SPD requierement and without ITSRefit", nptbins,ptmin,ptmax);
  fhPtNoSPDNoRefit->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtNoSPDNoRefit);
  
  fhEtaPhiNoSPDNoRefitPt02  = new TH2F ("hEtaPhiNoSPDNoRefitPt02","eta vs phi of constrained tracks with no SPD requierement and without ITSRefit, p_{T} < 2 GeV/c",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiNoSPDNoRefitPt02->SetXTitle("#eta ");
  fhEtaPhiNoSPDNoRefitPt02->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiNoSPDNoRefitPt02);
  
  fhEtaPhiNoSPDNoRefitPt3  = new TH2F ("hEtaPhiNoSPDNoRefitPt3","eta vs phi of constrained tracks with no SPD requierement and without ITSRefit, p_{T} > 3 GeV/c",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiNoSPDNoRefitPt3->SetXTitle("#eta ");
  fhEtaPhiNoSPDNoRefitPt3->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiNoSPDNoRefitPt3);

  if(fFillVertexBC0Histograms)
  {
    fhProductionVertexBC      = new TH1F("hProductionVertexBC", "tracks production vertex bunch crossing ", 41 , -20 , 20 ) ;
    fhProductionVertexBC->SetYTitle("# tracks");
    fhProductionVertexBC->SetXTitle("Bunch crossing");
    outputContainer->Add(fhProductionVertexBC);
  }
  
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
  
  fhPtTOFSignalDCACut  = new TH2F ("hPtTOFSignalDCACut","TOF signal after DCA cut", nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
  fhPtTOFSignalDCACut->SetYTitle("TOF signal (ns)");
  fhPtTOFSignalDCACut->SetXTitle("p_{T} (GeV/c)");
  outputContainer->Add(fhPtTOFSignalDCACut);

  if(fFillVertexBC0Histograms)
  {
    fhPtTOFSignalVtxOutBC0  = new TH2F ("hPtTOFSignalVtxOutBC0","TOF signal, vtx BC!=0", nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
    fhPtTOFSignalVtxOutBC0->SetYTitle("TOF signal (ns)");
    fhPtTOFSignalVtxOutBC0->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtTOFSignalVtxOutBC0);
    
    fhPtTOFSignalVtxInBC0  = new TH2F ("hPtTOFSignalVtxInBC0","TOF signal, vtx BC=0", nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
    fhPtTOFSignalVtxInBC0->SetYTitle("TOF signal (ns)");
    fhPtTOFSignalVtxInBC0->SetXTitle("p_{T} (GeV/c)");
    outputContainer->Add(fhPtTOFSignalVtxInBC0);
  }
  
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
      
      if(fFillVertexBC0Histograms)
      {
        fhPtTOFSignalVtxOutBC0PileUp[i]  = new TH2F(Form("hPtTOFSignalVtxOutBC0PileUp%s",pileUpName[i].Data()),
                                           Form("Track TOF vs p_{T} distribution, %s Pile-Up event, vtx BC!=0",pileUpName[i].Data()),
                                           nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
        fhPtTOFSignalVtxOutBC0PileUp[i]->SetXTitle("p_{T} (GeV/c)");
        fhPtTOFSignalVtxOutBC0PileUp[i]->SetXTitle("TOF signal (ns)");
        outputContainer->Add(fhPtTOFSignalVtxOutBC0PileUp[i]);

        fhPtTOFSignalVtxInBC0PileUp[i]  = new TH2F(Form("hPtTOFSignalVtxInBC0PileUp%s",pileUpName[i].Data()),
                                                    Form("Track TOF vs p_{T} distribution, %s Pile-Up event, vtx BC=0",pileUpName[i].Data()),
                                                    nptbins,ptmin,ptmax,ntofbins,mintof,maxtof);
        fhPtTOFSignalVtxInBC0PileUp[i]->SetXTitle("p_{T} (GeV/c)");
        fhPtTOFSignalVtxInBC0PileUp[i]->SetXTitle("TOF signal (ns)");
        outputContainer->Add(fhPtTOFSignalVtxInBC0PileUp[i]);
      }
      
      if(fFillVertexBC0Histograms)
      {
        fhProductionVertexBCPileUp[i]      = new TH1F(Form("hProductionVertexBCPileUp%s",pileUpName[i].Data()),
                                                Form("tracks production vertex bunch crossing, %s Pile-Up event",pileUpName[i].Data()),
                                                41 , -20 , 20 ) ;
        fhProductionVertexBCPileUp[i]->SetYTitle("# tracks");
        fhProductionVertexBCPileUp[i]->SetXTitle("Bunch crossing");
        outputContainer->Add(fhProductionVertexBCPileUp[i]);
      }
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
    fhPtDCA[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCA[i]);
    
    fhPtDCASPDRefit[i]  = new TH2F(Form("hPtDCA%sSPDRefit",dcaName[i].Data()),
                                        Form("Track DCA%s vs p_{T} distribution of tracks with SPD and ITS refit",dcaName[i].Data()),
                                        nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCASPDRefit[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCASPDRefit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCASPDRefit[i]);

    fhPtDCANoSPDRefit[i]  = new TH2F(Form("hPtDCA%sNoSPDRefit",dcaName[i].Data()),
                                 Form("Track DCA%s vs p_{T} distributionof constrained tracks no SPD and with ITSRefit",dcaName[i].Data()),
                                 nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCANoSPDRefit[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCANoSPDRefit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCANoSPDRefit[i]);

    fhPtDCANoSPDNoRefit[i]  = new TH2F(Form("hPtDCA%sNoSPDNoRefit",dcaName[i].Data()),
                           Form("Track DCA%s vs p_{T} distribution, constrained tracks with no SPD requierement and without ITSRefit",dcaName[i].Data()),
                           nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCANoSPDNoRefit[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCANoSPDNoRefit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCANoSPDNoRefit[i]);
    
    fhPtDCATOFBC0[i]  = new TH2F(Form("hPtDCA%sTOFBC0",dcaName[i].Data()),
                           Form("Track DCA%s vs p_{T} distribution, BC=0",dcaName[i].Data()),
                           nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCATOFBC0[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCATOFBC0[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCATOFBC0[i]);

    fhPtDCATOFBCOut[i]  = new TH2F(Form("hPtDCA%sTOFBCOut",dcaName[i].Data()),
                                 Form("Track DCA%s vs p_{T} distribution, BC!=0",dcaName[i].Data()),
                                 nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCATOFBCOut[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCATOFBCOut[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCATOFBCOut[i]);
    
    fhPtDCANoTOFHit[i]  = new TH2F(Form("hPtDCA%sNoTOFHit",dcaName[i].Data()),
                           Form("Track (no TOF hit) DCA%s vs p_{T} distribution",dcaName[i].Data()),
                           nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
    fhPtDCANoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
    fhPtDCANoTOFHit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
    outputContainer->Add(fhPtDCANoTOFHit[i]);

    if(fFillVertexBC0Histograms)
    {
      fhPtDCAVtxOutBC0[i]  = new TH2F(Form("hPtDCA%sVtxOutBC0",dcaName[i].Data()),
                                      Form("Track DCA%s vs p_{T} distribution, vertex with BC!=0",dcaName[i].Data()),
                                      nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAVtxOutBC0[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAVtxOutBC0[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAVtxOutBC0[i]);
      
      fhPtDCAVtxOutBC0NoTOFHit[i]  = new TH2F(Form("hPtDCA%sVtxOutBC0NoTOFHit",dcaName[i].Data()),
                                              Form("Track (no TOF hit) DCA%s vs p_{T} distribution, vertex with BC!=0",dcaName[i].Data()),
                                              nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAVtxOutBC0NoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAVtxOutBC0NoTOFHit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAVtxOutBC0NoTOFHit[i]);
      
      fhPtDCAVtxInBC0[i]  = new TH2F(Form("hPtDCA%sVtxInBC0",dcaName[i].Data()),
                                      Form("Track DCA%s vs p_{T} distribution, vertex with BC==0",dcaName[i].Data()),
                                      nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAVtxInBC0[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAVtxInBC0[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAVtxInBC0[i]);
      
      fhPtDCAVtxInBC0NoTOFHit[i]  = new TH2F(Form("hPtDCA%sVtxInBC0NoTOFHit",dcaName[i].Data()),
                                              Form("Track (no TOF hit) DCA%s vs p_{T} distribution, vertex with BC==0",dcaName[i].Data()),
                                              nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAVtxInBC0NoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAVtxInBC0NoTOFHit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAVtxInBC0NoTOFHit[i]);
    }
    
    if(fFillPileUpHistograms)
    {
      fhPtDCAPileUp[i]  = new TH2F(Form("hPtDCA%sPileUp",dcaName[i].Data()),
                             Form("Track DCA%s vs p_{T} distribution, SPD Pile-Up",dcaName[i].Data()),
                             nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAPileUp[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAPileUp[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAPileUp[i]);
      
      fhPtDCAPileUpTOFBC0[i]  = new TH2F(Form("hPtDCA%sPileUpTOFBC0",dcaName[i].Data()),
                                   Form("Track DCA%s vs p_{T} distribution",dcaName[i].Data()),
                                   nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAPileUpTOFBC0[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAPileUpTOFBC0[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAPileUpTOFBC0[i]);
      
      fhPtDCAPileUpNoTOFHit[i]  = new TH2F(Form("hPtDCA%sPileUpNoTOFHit",dcaName[i].Data()),
                                     Form("Track (no TOF hit) DCA%s vs p_{T} distribution, SPD Pile-Up, vertex with BC!=0",dcaName[i].Data()),
                                     nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
      fhPtDCAPileUpNoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
      fhPtDCAPileUpNoTOFHit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
      outputContainer->Add(fhPtDCAPileUpNoTOFHit[i]);
      
      if(fFillVertexBC0Histograms)
      {
        fhPtDCAVtxOutBC0PileUp[i]  = new TH2F(Form("hPtDCA%sPileUpVtxOutBC0",dcaName[i].Data()),
                                              Form("Track DCA%s vs p_{T} distribution, SPD Pile-Up, vertex with BC!=0",dcaName[i].Data()),
                                              nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
        fhPtDCAVtxOutBC0PileUp[i]->SetXTitle("p_{T} (GeV/c)");
        fhPtDCAVtxOutBC0PileUp[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
        outputContainer->Add(fhPtDCAVtxOutBC0PileUp[i]);
        
        fhPtDCAVtxOutBC0PileUpNoTOFHit[i]  = new TH2F(Form("hPtDCA%sVtxOutBC0PileUpNoTOFHit",dcaName[i].Data()),
                                                      Form("Track (no TOF hit) DCA%s vs p_{T} distribution, SPD Pile-Up, vertex with BC!=0",dcaName[i].Data()),
                                                      nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
        fhPtDCAVtxOutBC0PileUpNoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
        fhPtDCAVtxOutBC0PileUpNoTOFHit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
        outputContainer->Add(fhPtDCAVtxOutBC0PileUpNoTOFHit[i]);
        
        fhPtDCAVtxInBC0PileUp[i]  = new TH2F(Form("hPtDCA%sPileUpVtxInBC0",dcaName[i].Data()),
                                              Form("Track DCA%s vs p_{T} distribution, SPD Pile-Up,vertex with BC==0",dcaName[i].Data()),
                                              nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
        fhPtDCAVtxInBC0PileUp[i]->SetXTitle("p_{T} (GeV/c)");
        fhPtDCAVtxInBC0PileUp[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
        outputContainer->Add(fhPtDCAVtxInBC0PileUp[i]);
        
        fhPtDCAVtxInBC0PileUpNoTOFHit[i]  = new TH2F(Form("hPtDCA%sVtxInBC0PileUpNoTOFHit",dcaName[i].Data()),
                                                      Form("Track (no TOF hit) DCA%s vs p_{T} distribution, SPD Pile-Up, vertex with BC==0",dcaName[i].Data()),
                                                      nptbins,ptmin,ptmax,ndcabins,mindca,maxdca);
        fhPtDCAVtxInBC0PileUpNoTOFHit[i]->SetXTitle("p_{T} (GeV/c)");
        fhPtDCAVtxInBC0PileUpNoTOFHit[i]->SetYTitle(Form("DCA_{%s}",dcaName[i].Data()));
        outputContainer->Add(fhPtDCAVtxInBC0PileUpNoTOFHit[i]);
        
      }
    }
  }

  if(IsDataMC())
  {
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
  
  AliVEvent * event = GetReader()->GetInputEvent();

  Int_t vtxBC = GetReader()->GetVertexBC();
  if(!GetReader()->IsDCACutOn()) vtxBC = GetReader()->GetVertexBC(event->GetPrimaryVertex());

  if(fFillVertexBC0Histograms)
  {
    fhProductionVertexBC->Fill(vtxBC);
    if(fFillPileUpHistograms)
    {
      if(GetReader()->IsPileUpFromSPD())               fhProductionVertexBCPileUp[0]->Fill(vtxBC);
      if(GetReader()->IsPileUpFromEMCal())             fhProductionVertexBCPileUp[1]->Fill(vtxBC);
      if(GetReader()->IsPileUpFromSPDOrEMCal())        fhProductionVertexBCPileUp[2]->Fill(vtxBC);
      if(GetReader()->IsPileUpFromSPDAndEMCal())       fhProductionVertexBCPileUp[3]->Fill(vtxBC);
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhProductionVertexBCPileUp[4]->Fill(vtxBC);
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhProductionVertexBCPileUp[5]->Fill(vtxBC);
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhProductionVertexBCPileUp[6]->Fill(vtxBC);
    }
  }
  
  //printf("AliAnaChargedParticles::MakeAnalysisFillAOD() - primary vertex BC %d\n",vtxBC);
  
  Double_t bz = event->GetMagneticField();

  //Fill AODParticle with CTS aods
  Float_t pt  = 0;
  Float_t phi = 0;
  Float_t eta = 0;
  Int_t evtIndex = 0;
  for(Int_t i = 0; i < ntracks; i++)
  {
    AliVTrack * track =  (AliVTrack*) (GetCTSTracks()->At(i));
        
    pt  = track->Pt();
    eta = track->Eta();
    phi = track->Phi();

    fhPtNoCut->Fill(pt);
    
    AliAODTrack * aodTrack = dynamic_cast<AliAODTrack*>(track);
    AliESDtrack * esdTrack = dynamic_cast<AliESDtrack*>(track);
    
    //TOF
    ULong_t status = track->GetStatus();
    Bool_t okTOF = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
    Double32_t tof = track->GetTOFsignal()*1e-3;    
    
    //DCA
    Double_t dcaCons = -999;
    if(aodTrack)
    {
      dcaCons = aodTrack->DCA();
      //vtxBC   = aodTrack->GetProdVertex()->GetBC();
    }
    
    Double_t dca[2]   = {1e6,1e6};
    Double_t covar[3] = {1e6,1e6,1e6};
    track->PropagateToDCA(GetReader()->GetInputEvent()->GetPrimaryVertex(),bz,100.,dca,covar);
    
    Float_t trackDCA = dca[0];
    
    if(dcaCons == -999)
    {
      fhPtDCA[0]->Fill(pt, dca[0]);
      fhPtDCA[1]->Fill(pt, dca[1]);
    }
    else
    {
      trackDCA = dcaCons;
      fhPtDCA[2]->Fill(pt, dcaCons);
    }
    
    if(GetReader()->AcceptDCA(pt,trackDCA)) fhPtCutDCA->Fill(pt);
    
    if(fFillVertexBC0Histograms)
    {
      if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
      {        
        fhPtVtxOutBC0->Fill(pt);
        fhEtaPhiVtxOutBC0->Fill(eta,phi);
        
        if(dcaCons == -999)
        {
          fhPtDCAVtxOutBC0[0]->Fill(pt, dca[0]);
          fhPtDCAVtxOutBC0[1]->Fill(pt, dca[1]);
        }
        else
          fhPtDCAVtxOutBC0[2]->Fill(pt, dcaCons);
      }
      else
      {
        fhPtVtxInBC0->Fill(pt);
        fhEtaPhiVtxInBC0->Fill(eta,phi);
        if(GetReader()->AcceptDCA(pt,trackDCA)) fhPtCutDCABCOK->Fill(pt);
        
        if(dcaCons == -999)
        {
          fhPtDCAVtxInBC0[0]->Fill(pt, dca[0]);
          fhPtDCAVtxInBC0[1]->Fill(pt, dca[1]);
        }
        else
          fhPtDCAVtxInBC0[2]->Fill(pt, dcaCons);
        
      }
    }
    
    if(fFillPileUpHistograms && GetReader()->IsPileUpFromSPD())
    {
      if(dcaCons == -999)
      {
        fhPtDCAPileUp[0]->Fill(pt, dca[0]);
        fhPtDCAPileUp[1]->Fill(pt, dca[1]);
      }
      else
        fhPtDCAPileUp[2]->Fill(pt, dcaCons);

      if(fFillVertexBC0Histograms)
      {
        if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
        {
          if(dcaCons == -999)
          {
            fhPtDCAVtxOutBC0PileUp[0]->Fill(pt, dca[0]);
            fhPtDCAVtxOutBC0PileUp[1]->Fill(pt, dca[1]);
          }
          else fhPtDCAVtxOutBC0PileUp[2]->Fill(pt, dcaCons);
        }
        else
        {
          if(dcaCons == -999)
          {
            fhPtDCAVtxInBC0PileUp[0]->Fill(pt, dca[0]);
            fhPtDCAVtxInBC0PileUp[1]->Fill(pt, dca[1]);
          }
          else fhPtDCAVtxInBC0PileUp[2]->Fill(pt, dcaCons);
        }
      }
    }

    if(!okTOF)
    {
      if(dcaCons == -999)
      {
        fhPtDCANoTOFHit[0]->Fill(pt, dca[0]);
        fhPtDCANoTOFHit[1]->Fill(pt, dca[1]);
      }
      else
        fhPtDCANoTOFHit[2]->Fill(pt, dcaCons);
      
      if(fFillVertexBC0Histograms)
      {
        if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
        {
          if(dcaCons == -999)
          {
            fhPtDCAVtxOutBC0NoTOFHit[0]->Fill(pt, dca[0]);
            fhPtDCAVtxOutBC0NoTOFHit[1]->Fill(pt, dca[1]);
          }
          else
            fhPtDCAVtxOutBC0NoTOFHit[2]->Fill(pt, dcaCons);
        }
        else
        {
          if(dcaCons == -999)
          {
            fhPtDCAVtxInBC0NoTOFHit[0]->Fill(pt, dca[0]);
            fhPtDCAVtxInBC0NoTOFHit[1]->Fill(pt, dca[1]);
          }
          else
            fhPtDCAVtxInBC0NoTOFHit[2]->Fill(pt, dcaCons);
          
        }
      }
      
      if(fFillPileUpHistograms && GetReader()->IsPileUpFromSPD())
      {
        if(dcaCons == -999)
        {
          fhPtDCAPileUpNoTOFHit[0]->Fill(pt, dca[0]);
          fhPtDCAPileUpNoTOFHit[1]->Fill(pt, dca[1]);
        }
        else
          fhPtDCAPileUpNoTOFHit[2]->Fill(pt, dcaCons);
        
        if(fFillVertexBC0Histograms)
        {
          if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
          {
            if(dcaCons == -999)
            {
              fhPtDCAVtxOutBC0PileUpNoTOFHit[0]->Fill(pt, dca[0]);
              fhPtDCAVtxOutBC0PileUpNoTOFHit[1]->Fill(pt, dca[1]);
            }
            else
              fhPtDCAVtxOutBC0PileUpNoTOFHit[2]->Fill(pt, dcaCons);
          }
          else
          {
            if(dcaCons == -999)
            {
              fhPtDCAVtxInBC0PileUpNoTOFHit[0]->Fill(pt, dca[0]);
              fhPtDCAVtxInBC0PileUpNoTOFHit[1]->Fill(pt, dca[1]);
            }
            else
              fhPtDCAVtxInBC0PileUpNoTOFHit[2]->Fill(pt, dcaCons);
            
          }
        }
      }
    }
    
    //printf("track pT %2.2f, DCA Cons %f, DCA1 %f, DCA2 %f, TOFBC %d, oktof %d, tof %f\n",
    //      pt,dcaCons,dca[0],dca[1],track->GetTOFBunchCrossing(bz),okTOF, tof);
    
    Int_t trackBC = track->GetTOFBunchCrossing(bz);
    
//    if( vtxBC == 0 && trackBC !=0 && trackBC!=AliVTrack::kTOFBCNA)
//      printf("TOF Signal %e, BC %d, pt %f, dca_xy %f, dca_z %f, dca_tpc %f \n", tof,trackBC, pt,dca[0],dca[1],dcaCons);
    
    
    if(okTOF)
    {
      fhTOFSignal  ->Fill(tof);
      fhPtTOFSignal->Fill(pt, tof);
      if(GetReader()->AcceptDCA(pt,trackDCA)) fhPtTOFSignalDCACut->Fill(pt, tof);
        
      if(fFillVertexBC0Histograms)
      {
        if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
          fhPtTOFSignalVtxOutBC0->Fill(pt, tof);
        else
          fhPtTOFSignalVtxInBC0->Fill(pt, tof);
      }
      
      if(trackBC==0)
      {
        fhTOFSignalBCOK->Fill(tof);
        
        if(dcaCons == -999)
        {
          fhPtDCATOFBC0[0]->Fill(pt, dca[0]);
          fhPtDCATOFBC0[1]->Fill(pt, dca[1]);
        }
        else
          fhPtDCATOFBC0[2]->Fill(pt, dcaCons);
        
        if(fFillPileUpHistograms && GetReader()->IsPileUpFromSPD())
        {
          if(dcaCons == -999)
          {
            fhPtDCAPileUpTOFBC0[0]->Fill(pt, dca[0]);
            fhPtDCAPileUpTOFBC0[1]->Fill(pt, dca[1]);
          }
          else
            fhPtDCAPileUpTOFBC0[2]->Fill(pt, dcaCons);
        }
      }
      else if(trackBC!=AliVTrack::kTOFBCNA)
      {
        if(dcaCons == -999)
        {
          fhPtDCATOFBCOut[0]->Fill(pt, dca[0]);
          fhPtDCATOFBCOut[1]->Fill(pt, dca[1]);
        }
        else
          fhPtDCATOFBCOut[2]->Fill(pt, dcaCons);
      
      }

      if(fFillPileUpHistograms)
      {
        if(GetReader()->IsPileUpFromSPD())               fhPtTOFSignalPileUp[0]->Fill(pt, tof);
        if(GetReader()->IsPileUpFromEMCal())             fhPtTOFSignalPileUp[1]->Fill(pt, tof);
        if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtTOFSignalPileUp[2]->Fill(pt, tof);
        if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtTOFSignalPileUp[3]->Fill(pt, tof);
        if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtTOFSignalPileUp[4]->Fill(pt, tof);
        if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtTOFSignalPileUp[5]->Fill(pt, tof);
        if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtTOFSignalPileUp[6]->Fill(pt, tof);
        
        if      (trackBC ==0)  { fhEtaPhiTOFBC0    ->Fill(eta,phi); if(GetReader()->IsPileUpFromSPD()) fhEtaPhiTOFBC0PileUpSPD    ->Fill(eta,phi); }
        else if (trackBC < 0)  { fhEtaPhiTOFBCPlus ->Fill(eta,phi); if(GetReader()->IsPileUpFromSPD()) fhEtaPhiTOFBCPlusPileUpSPD ->Fill(eta,phi); }
        else if (trackBC > 0)  { fhEtaPhiTOFBCMinus->Fill(eta,phi); if(GetReader()->IsPileUpFromSPD()) fhEtaPhiTOFBCMinusPileUpSPD->Fill(eta,phi); }
        
        if(fFillVertexBC0Histograms)
        {
          if(TMath::Abs(vtxBC) > 0 && vtxBC!=AliVTrack::kTOFBCNA)
          {            
            if(GetReader()->IsPileUpFromSPD())               fhPtTOFSignalVtxOutBC0PileUp[0]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromEMCal())             fhPtTOFSignalVtxOutBC0PileUp[1]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtTOFSignalVtxOutBC0PileUp[2]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtTOFSignalVtxOutBC0PileUp[3]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtTOFSignalVtxOutBC0PileUp[4]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtTOFSignalVtxOutBC0PileUp[5]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtTOFSignalVtxOutBC0PileUp[6]->Fill(pt, tof);
          }
          else
          {            
            if(GetReader()->IsPileUpFromSPD())               fhPtTOFSignalVtxInBC0PileUp[0]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromEMCal())             fhPtTOFSignalVtxInBC0PileUp[1]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromSPDOrEMCal())        fhPtTOFSignalVtxInBC0PileUp[2]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromSPDAndEMCal())       fhPtTOFSignalVtxInBC0PileUp[3]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhPtTOFSignalVtxInBC0PileUp[4]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhPtTOFSignalVtxInBC0PileUp[5]->Fill(pt, tof);
            if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhPtTOFSignalVtxInBC0PileUp[6]->Fill(pt, tof);
          }
        }
      }
    }
    
    //Fill AODParticle after some selection
    Double_t mom[3] = {track->Px(),track->Py(),track->Pz()};
    
    Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,"CTS") ;
    
    if(GetDebug() > 1) 
      printf("AliAnaChargedParticles::MakeAnalysisFillAOD() - Track pt %2.2f, eta %2.2f, phi %2.2f in fiducial cut %d\n",pt,eta,phi,in);
    
    //Acceptance selection
    if(IsFiducialCutOn() && ! in ) continue ;
    
    // Momentum selection
    if(pt < GetMinPt() || pt > GetMaxPt()) continue;
    
    if(okTOF) fhTOFSignalPtCut->Fill(tof); 
    else
    {
      fhPtTOFStatus0    ->Fill(pt);
      fhEtaPhiTOFStatus0->Fill(eta,phi);
    }
    
    Bool_t bITSRefit    = (status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit;
    Bool_t bConstrained = kFALSE;
    if     (aodTrack) bConstrained = aodTrack->IsGlobalConstrained();
    else if(esdTrack) bConstrained = (!esdTrack->HasPointOnITSLayer(0) && !esdTrack->HasPointOnITSLayer(1));
    //printf("Track %d, pt %2.2f, eta %2.2f, phi %2.2f, SPDRefit %d, refit %d, dcaCons %2.2f\n",
    //       i, pt, eta, phi, bConstrained, bITSRefit, dcaCons);
    
    if(bConstrained)
    {      
      if(bITSRefit)
      {
        fhPtNoSPDRefit->Fill(pt);
        if(pt < 2)fhEtaPhiNoSPDRefitPt02->Fill(eta,phi);
        if(pt > 3)fhEtaPhiNoSPDRefitPt3 ->Fill(eta,phi);
        
        if(dcaCons == -999)
        {
          fhPtDCANoSPDRefit[0]->Fill(pt, dca[0]);
          fhPtDCANoSPDRefit[1]->Fill(pt, dca[1]);
        }
        else
          fhPtDCANoSPDRefit[2]->Fill(pt, dcaCons);
        
      }
      else
      {
        fhPtNoSPDNoRefit->Fill(pt);
        if(pt < 2)fhEtaPhiNoSPDNoRefitPt02->Fill(eta,phi);
        if(pt > 3)fhEtaPhiNoSPDNoRefitPt3 ->Fill(eta,phi);
        if(dcaCons == -999)
        {
          fhPtDCANoSPDNoRefit[0]->Fill(pt, dca[0]);
          fhPtDCANoSPDNoRefit[1]->Fill(pt, dca[1]);
        }
        else
          fhPtDCANoSPDNoRefit[2]->Fill(pt, dcaCons);

      }
    }
    else
    {
      fhPtSPDRefit->Fill(pt);
      if(pt < 2)fhEtaPhiSPDRefitPt02->Fill(eta,phi);
      if(pt > 3)fhEtaPhiSPDRefitPt3 ->Fill(eta,phi);
      if(dcaCons == -999)
      {
        fhPtDCASPDRefit[0]->Fill(pt, dca[0]);
        fhPtDCASPDRefit[1]->Fill(pt, dca[1]);
      }
      else
        fhPtDCASPDRefit[2]->Fill(pt, dcaCons);
    }
        
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
  
  Float_t pt  = 0;
  Float_t phi = 0;
  Float_t eta = 0;
  
  for(Int_t iaod = 0; iaod < naod ; iaod++)
  {
    AliAODPWG4Particle* track =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    
    pt  = track->Pt();
    eta = track->Eta();
    phi = track->Phi();
    
    fhPt->Fill(pt);
    
    if(track->GetChargedBit())
    {
      fhPhiPos   ->Fill(pt, phi);
      fhEtaPos   ->Fill(pt, eta);
      fhEtaPhiPos->Fill(eta,phi);
    }
    else
    {
      fhPhiNeg   ->Fill(pt, phi);
      fhEtaNeg   ->Fill(pt, eta);
      fhEtaPhiNeg->Fill(eta,phi);
    }
    
    if(fFillPileUpHistograms)
    {
      if(GetReader()->IsPileUpFromSPD())               {fhPtPileUp[0]->Fill(pt);}
      if(GetReader()->IsPileUpFromEMCal())             {fhPtPileUp[1]->Fill(pt);}
      if(GetReader()->IsPileUpFromSPDOrEMCal())        {fhPtPileUp[2]->Fill(pt);}
      if(GetReader()->IsPileUpFromSPDAndEMCal())       {fhPtPileUp[3]->Fill(pt);}
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    {fhPtPileUp[4]->Fill(pt);}
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    {fhPtPileUp[5]->Fill(pt);}
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) {fhPtPileUp[6]->Fill(pt);}
    }

    
    if(IsDataMC())
    {
      //Play with the MC stack if available		
      Int_t mompdg = -1;
      Int_t label  = track->GetLabel();
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
          aodmom = (AliAODMCParticle*) (GetReader()->GetAODMCParticles())->At(label);
          mompdg =TMath::Abs(aodmom->GetPdgCode());
        }
      }
      
      if(mompdg==211)
      {
        fhPtPion ->Fill(pt);
        fhPhiPion->Fill(pt, phi);
        fhEtaPion->Fill(pt, eta);
      }
      else if(mompdg==2212)
      {
        fhPtProton ->Fill(pt);
        fhPhiProton->Fill(pt, phi);
        fhEtaProton->Fill(pt, eta);
      }
      else if(mompdg==321)
      {
        fhPtKaon ->Fill(pt);
        fhPhiKaon->Fill(pt, phi);
        fhEtaKaon->Fill(pt, eta);
      }
      else if(mompdg==11)
      {
        fhPtElectron ->Fill(pt);
        fhPhiElectron->Fill(pt, phi);
        fhEtaElectron->Fill(pt, eta);
      }
      else
      {
        fhPtUnknown ->Fill(pt);
        fhPhiUnknown->Fill(pt, phi);
        fhEtaUnknown->Fill(pt, eta);
      }
      
    }//Work with stack also
    
  }// aod branch loop
  
}
