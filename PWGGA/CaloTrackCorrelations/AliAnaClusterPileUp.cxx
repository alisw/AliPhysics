/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//_________________________________________________________________________
//
// Class for the study of Pile-up effect on
// Calorimeter clusters.
// Open time cuts in reader.
//
//-- Author: Gustavo Conesa (CNRS-LPSC-Grenoble)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TH2F.h>
#include <TClonesArray.h>
#include <TObjString.h>

// --- Analysis system ---
#include "AliAnaClusterPileUp.h"
#include "AliCaloTrackReader.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"

ClassImp(AliAnaClusterPileUp)

//___________________________________________
AliAnaClusterPileUp::AliAnaClusterPileUp() :
AliAnaCaloTrackCorrBaseClass(),
fCalorimeter(""),                     fNCellsCut(0),
// Histograms
fhTimePtNoCut(0),                     fhTimePtSPD(0),
fhTimeNPileUpVertSPD(0),              fhTimeNPileUpVertTrack(0),
fhTimeNPileUpVertContributors(0),
fhTimePileUpMainVertexZDistance(0),   fhTimePileUpMainVertexZDiamond(0),
fhClusterMultSPDPileUp(),             fhClusterMultNoPileUp(),
fhEtaPhiBC0(0),  fhEtaPhiBCPlus(0),   fhEtaPhiBCMinus(0),
fhEtaPhiBC0PileUpSPD(0),
fhEtaPhiBCPlusPileUpSPD(0),           fhEtaPhiBCMinusPileUpSPD(0),
fhPtNPileUpSPDVtx(0),                 fhPtNPileUpTrkVtx(0),
fhPtNPileUpSPDVtxTimeCut(0),          fhPtNPileUpTrkVtxTimeCut(0),
fhPtNPileUpSPDVtxTimeCut2(0),         fhPtNPileUpTrkVtxTimeCut2(0)
{
  //default ctor

  for(Int_t i = 0; i < 7; i++)
  {
    fhPtPileUp       [i] = 0;
    fhPtNeutralPileUp[i] = 0;
    
    fhLambda0PileUp       [i] = 0;
    fhLambda0NeutralPileUp[i] = 0;
    
    fhClusterEFracLongTimePileUp  [i] = 0;
    
    fhClusterCellTimePileUp       [i] = 0;
    fhClusterTimeDiffPileUp       [i] = 0;
    fhClusterTimeDiffNeutralPileUp[i] = 0;
    
  }
  
  for(Int_t i = 0; i < 4; i++)
  {
    fhClusterMultSPDPileUp[i] = 0;
    fhClusterMultNoPileUp [i] = 0;
  }
  
  //Initialize parameters
  InitParameters();
  
}

//___________________________________________
TObjString *  AliAnaClusterPileUp::GetAnalysisCuts()
{
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaClusterPileUp ---\n") ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  
  //Get parameters set in base class.
  //parList += GetBaseParametersList() ;
  
  return new TObjString(parList) ;
}

//________________________________________________________________________
TList *  AliAnaClusterPileUp::GetCreateOutputObjects()
{
  // Create histograms to be saved in output file and
  // store them in outputContainer
  TList * outputContainer = new TList() ;
  outputContainer->SetName("PhotonHistos") ;
	
  Int_t nptbins  = GetHistogramRanges()->GetHistoPtBins();  Float_t ptimecluster  = GetHistogramRanges()->GetHistoPtMax();  Float_t ptmin  = GetHistogramRanges()->GetHistoPtMin();
  Int_t nphibins = GetHistogramRanges()->GetHistoPhiBins(); Float_t phimax = GetHistogramRanges()->GetHistoPhiMax(); Float_t phimin = GetHistogramRanges()->GetHistoPhiMin();
  Int_t netabins = GetHistogramRanges()->GetHistoEtaBins(); Float_t etamax = GetHistogramRanges()->GetHistoEtaMax(); Float_t etamin = GetHistogramRanges()->GetHistoEtaMin();
  Int_t ssbins   = GetHistogramRanges()->GetHistoShowerShapeBins();  Float_t ssmax   = GetHistogramRanges()->GetHistoShowerShapeMax();  Float_t ssmin   = GetHistogramRanges()->GetHistoShowerShapeMin();
  Int_t ntimebins= GetHistogramRanges()->GetHistoTimeBins();         Float_t timemax = GetHistogramRanges()->GetHistoTimeMax();         Float_t timemin = GetHistogramRanges()->GetHistoTimeMin();
  
  
  fhTimePtNoCut  = new TH2F ("hTimePt_NoCut","time of cluster vs pT of clusters, no event selection", nptbins,ptmin,ptimecluster, ntimebins,timemin,timemax);
  fhTimePtNoCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhTimePtNoCut->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimePtNoCut);
  
  fhTimePtSPD  = new TH2F ("hTimePt_SPD","time of cluster vs pT of clusters, SPD Pile-up events", nptbins,ptmin,ptimecluster, ntimebins,timemin,timemax);
  fhTimePtSPD->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fhTimePtSPD->SetYTitle("#it{time} (ns)");
  outputContainer->Add(fhTimePtSPD);
  
  fhTimeNPileUpVertSPD  = new TH2F ("hTime_NPileUpVertSPD","time of cluster vs N pile-up SPD vertex", ntimebins,timemin,timemax,20,0,20);
  fhTimeNPileUpVertSPD->SetYTitle("# vertex ");
  fhTimeNPileUpVertSPD->SetXTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeNPileUpVertSPD);
  
  fhTimeNPileUpVertTrack  = new TH2F ("hTime_NPileUpVertTracks","time of cluster vs N pile-up Tracks vertex", ntimebins,timemin,timemax, 20,0,20 );
  fhTimeNPileUpVertTrack->SetYTitle("# vertex ");
  fhTimeNPileUpVertTrack->SetXTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeNPileUpVertTrack);
  
  fhTimeNPileUpVertContributors  = new TH2F ("hTime_NPileUpVertContributors","time of cluster vs N constributors to pile-up SPD vertex", ntimebins,timemin,timemax,50,0,50);
  fhTimeNPileUpVertContributors->SetYTitle("# vertex ");
  fhTimeNPileUpVertContributors->SetXTitle("#it{time} (ns)");
  outputContainer->Add(fhTimeNPileUpVertContributors);
  
  fhTimePileUpMainVertexZDistance  = new TH2F ("hTime_PileUpMainVertexZDistance","time of cluster vs distance in Z pile-up SPD vertex - main SPD vertex",ntimebins,timemin,timemax,100,0,50);
  fhTimePileUpMainVertexZDistance->SetYTitle("distance Z (cm) ");
  fhTimePileUpMainVertexZDistance->SetXTitle("#it{time} (ns)");
  outputContainer->Add(fhTimePileUpMainVertexZDistance);
  
  fhTimePileUpMainVertexZDiamond  = new TH2F ("hTime_PileUpMainVertexZDiamond","time of cluster vs distance in Z pile-up SPD vertex - z diamond",ntimebins,timemin,timemax,100,0,50);
  fhTimePileUpMainVertexZDiamond->SetYTitle("diamond distance Z (cm) ");
  fhTimePileUpMainVertexZDiamond->SetXTitle("#it{time} (ns)");
  outputContainer->Add(fhTimePileUpMainVertexZDiamond);

  fhEtaPhiBC0  = new TH2F ("hEtaPhiBC0","eta-phi for clusters tof corresponding to BC=0",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiBC0->SetXTitle("#eta ");
  fhEtaPhiBC0->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiBC0);
  
  fhEtaPhiBCPlus  = new TH2F ("hEtaPhiBCPlus","eta-phi for clusters tof corresponding to BC>0",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiBCPlus->SetXTitle("#eta ");
  fhEtaPhiBCPlus->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiBCPlus);
  
  fhEtaPhiBCMinus  = new TH2F ("hEtaPhiBCMinus","eta-phi for clusters tof corresponding to BC<0",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiBCMinus->SetXTitle("#eta ");
  fhEtaPhiBCMinus->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiBCMinus);
  
  fhEtaPhiBC0PileUpSPD  = new TH2F ("hEtaPhiBC0PileUpSPD","eta-phi for clusters tof corresponding to BC=0, SPD pile-up",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiBC0PileUpSPD->SetXTitle("#eta ");
  fhEtaPhiBC0PileUpSPD->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiBC0PileUpSPD);
  
  fhEtaPhiBCPlusPileUpSPD  = new TH2F ("hEtaPhiBCPlusPileUpSPD","eta-phi for clusters tof corresponding to BC>0, SPD pile-up",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiBCPlusPileUpSPD->SetXTitle("#eta ");
  fhEtaPhiBCPlusPileUpSPD->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiBCPlusPileUpSPD);
  
  fhEtaPhiBCMinusPileUpSPD  = new TH2F ("hEtaPhiBCMinusPileUpSPD","eta-phi for clusters tof corresponding to BC<0, SPD pile-up",netabins,etamin,etamax, nphibins,phimin,phimax);
  fhEtaPhiBCMinusPileUpSPD->SetXTitle("#eta ");
  fhEtaPhiBCMinusPileUpSPD->SetYTitle("#phi (rad)");
  outputContainer->Add(fhEtaPhiBCMinusPileUpSPD);
  
  TString title[] = {"no |t diff| cut","|t diff|<20 ns","|t diff|>20 ns","|t diff|>40 ns"};
  TString name [] = {"TDiffNoCut","TDiffSmaller20ns","TDiffLarger20ns","TDiffLarger40ns"};
  for(Int_t i = 0; i < 4; i++)
  {
    fhClusterMultSPDPileUp[i] = new TH2F(Form("fhClusterMultSPDPileUp_%s", name[i].Data()),
                                         Form("Number of clusters per pile up event with #it{E} > 0.5 and %s respect cluster max vs cluster max E ",title[i].Data()),
                                         nptbins,ptmin,ptimecluster,100,0,100);
    fhClusterMultSPDPileUp[i]->SetYTitle("n clusters ");
    fhClusterMultSPDPileUp[i]->SetXTitle("#it{E}_{cluster max} (GeV)");
    outputContainer->Add(fhClusterMultSPDPileUp[i]) ;
    
    fhClusterMultNoPileUp[i] = new TH2F(Form("fhClusterMultNoPileUp_%s", name[i].Data()),
                                        Form("Number of clusters per non pile up event with #it{E} > 0.5 and %s respect cluster max vs cluster max E ",title[i].Data()),
                                        nptbins,ptmin,ptimecluster,100,0,100);
    fhClusterMultNoPileUp[i]->SetYTitle("n clusters ");
    fhClusterMultNoPileUp[i]->SetXTitle("#it{E}_{cluster max} (GeV)");
    outputContainer->Add(fhClusterMultNoPileUp[i]) ;
  }
  
  fhPtNPileUpSPDVtx  = new TH2F ("hPt_NPileUpVertSPD","pT of cluster vs N pile-up SPD vertex",
                                 nptbins,ptmin,ptimecluster,20,0,20);
  fhPtNPileUpSPDVtx->SetYTitle("# vertex ");
  fhPtNPileUpSPDVtx->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtNPileUpSPDVtx);
  
  fhPtNPileUpTrkVtx  = new TH2F ("hPt_NPileUpVertTracks","pT of cluster vs N pile-up Tracks vertex",
                                 nptbins,ptmin,ptimecluster, 20,0,20 );
  fhPtNPileUpTrkVtx->SetYTitle("# vertex ");
  fhPtNPileUpTrkVtx->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtNPileUpTrkVtx);
  
  fhPtNPileUpSPDVtxTimeCut  = new TH2F ("hPt_NPileUpVertSPD_TimeCut","pT of cluster vs N pile-up SPD vertex, |tof| < 25 ns",
                                        nptbins,ptmin,ptimecluster,20,0,20);
  fhPtNPileUpSPDVtxTimeCut->SetYTitle("# vertex ");
  fhPtNPileUpSPDVtxTimeCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtNPileUpSPDVtxTimeCut);
  
  fhPtNPileUpTrkVtxTimeCut  = new TH2F ("hPt_NPileUpVertTracks_TimeCut","pT of cluster vs N pile-up Tracks vertex, |tof| < 25 ns",
                                        nptbins,ptmin,ptimecluster, 20,0,20 );
  fhPtNPileUpTrkVtxTimeCut->SetYTitle("# vertex ");
  fhPtNPileUpTrkVtxTimeCut->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtNPileUpTrkVtxTimeCut);
  
  fhPtNPileUpSPDVtxTimeCut2  = new TH2F ("hPt_NPileUpVertSPD_TimeCut2","pT of cluster vs N pile-up SPD vertex, -25 < tof < 75 ns",
                                         nptbins,ptmin,ptimecluster,20,0,20);
  fhPtNPileUpSPDVtxTimeCut2->SetYTitle("# vertex ");
  fhPtNPileUpSPDVtxTimeCut2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtNPileUpSPDVtxTimeCut2);
  
  fhPtNPileUpTrkVtxTimeCut2  = new TH2F ("hPt_NPileUpVertTracks_TimeCut2","pT of cluster vs N pile-up Tracks vertex, -25 < tof < 75 ns",
                                         nptbins,ptmin,ptimecluster, 20,0,20 );
  fhPtNPileUpTrkVtxTimeCut2->SetYTitle("# vertex ");
  fhPtNPileUpTrkVtxTimeCut2->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  outputContainer->Add(fhPtNPileUpTrkVtxTimeCut2);
  
  
  TString pileUpName[] = {"SPD","EMCAL","SPDOrEMCAL","SPDAndEMCAL","SPDAndNotEMCAL","EMCALAndNotSPD","NotSPDAndNotEMCAL"} ;
  
  for(Int_t i = 0 ; i < 7 ; i++)
  {
    fhPtPileUp[i]  = new TH1F(Form("hPtPileUp%s",pileUpName[i].Data()),
                              Form("Cluster  #it{p}_{T} distribution, %s Pile-Up event",pileUpName[i].Data()), nptbins,ptmin,ptimecluster);
    fhPtPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtPileUp[i]);
    
    fhPtNeutralPileUp[i]  = new TH1F(Form("hPtNeutralPileUp%s",pileUpName[i].Data()),
                                     Form("Neutral clusters #it{p}_{T} distribution, %s Pile-Up event",pileUpName[i].Data()), nptbins,ptmin,ptimecluster);
    fhPtNeutralPileUp[i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    outputContainer->Add(fhPtNeutralPileUp[i]);
    
    fhClusterEFracLongTimePileUp[i]  = new TH2F(Form("hClusterEFracLongTimePileUp%s",pileUpName[i].Data()),
                                                Form("Cluster E vs fraction of cluster energy from large T cells, %s Pile-Up event",pileUpName[i].Data()),
                                                nptbins,ptmin,ptimecluster,200,0,1);
    fhClusterEFracLongTimePileUp[i]->SetXTitle("#it{E} (GeV)");
    fhClusterEFracLongTimePileUp[i]->SetYTitle("E(large time) / E");
    outputContainer->Add(fhClusterEFracLongTimePileUp[i]);
    
    fhClusterCellTimePileUp[i]  = new TH2F(Form("hClusterCellTimePileUp%s",pileUpName[i].Data()),
                                           Form("Cluster E vs cell time in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                           nptbins,ptmin,ptimecluster,ntimebins,timemin,timemax);
    fhClusterCellTimePileUp[i]->SetXTitle("#it{E} (GeV)");
    fhClusterCellTimePileUp[i]->SetYTitle("t_{cell} (ns)");
    outputContainer->Add(fhClusterCellTimePileUp[i]);
    
    fhClusterTimeDiffPileUp[i]  = new TH2F(Form("hClusterTimeDiffPileUp%s",pileUpName[i].Data()),
                                           Form("Cluster E vs t_{max}-t_{cell} in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                           nptbins,ptmin,ptimecluster,400,-200,200);
    fhClusterTimeDiffPileUp[i]->SetXTitle("#it{E} (GeV)");
    fhClusterTimeDiffPileUp[i]->SetYTitle("t_{max}-t_{cell} (ns)");
    outputContainer->Add(fhClusterTimeDiffPileUp[i]);
    
    fhClusterTimeDiffNeutralPileUp[i]  = new TH2F(Form("hClusterTimeDiffNeutralPileUp%s",pileUpName[i].Data()),
                                                  Form("Neutral clusters E vs t_{max}-t_{cell} in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                                  nptbins,ptmin,ptimecluster,400,-200,200);
    fhClusterTimeDiffNeutralPileUp[i]->SetXTitle("#it{E} (GeV)");
    fhClusterTimeDiffNeutralPileUp[i]->SetYTitle("t_{max}-t_{cell} (ns)");
    outputContainer->Add(fhClusterTimeDiffNeutralPileUp[i]);
    
    fhLambda0PileUp[i]  = new TH2F(Form("hLambda0PileUp%s",pileUpName[i].Data()),
                                   Form("Cluster E vs #lambda^{2}_{0} in cluster, %s Pile-Up event",pileUpName[i].Data()),
                                   nptbins,ptmin,ptimecluster,ssbins,ssmin,ssmax);
    fhLambda0PileUp[i]->SetXTitle("#it{E} (GeV)");
    fhLambda0PileUp[i]->SetYTitle("#lambda^{2}_{0}");
    outputContainer->Add(fhLambda0PileUp[i]);
    
    fhLambda0NeutralPileUp[i]  = new TH2F(Form("hLambda0NeutralPileUp%s",pileUpName[i].Data()),
                                          Form("Neutral clusters E vs #lambda^{2}_{0}in cluster, %s Pile-Up event",pileUpName[i].Data()), nptbins,ptmin,ptimecluster,ssbins,ssmin,ssmax);
    fhLambda0NeutralPileUp[i]->SetXTitle("#it{E} (GeV)");
    fhLambda0NeutralPileUp[i]->SetYTitle("#lambda^{2}_{0}");
    outputContainer->Add(fhLambda0NeutralPileUp[i]);
    
  }
  
  return outputContainer ;
  
}

//_______________________
void AliAnaClusterPileUp::Init()
{
  //Init
  
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn())
    AliFatal("You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!");
  
  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn())
    AliFatal("You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!");
  
  if(GetReader()->GetDataType() == AliCaloTrackReader::kMC)
    AliFatal("You want to use MC data in analysis but this is not possible in pile-up!!");
  
}

//____________________________________________________________________________
void AliAnaClusterPileUp::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  AddToHistogramsName("AnaClusterPileUp_");
  
  fCalorimeter = "EMCAL" ;
 	fNCellsCut   = 2;
}

//__________________________________________________________________
void  AliAnaClusterPileUp::MakeAnalysisFillHistograms()
{
  // Do cluster analysis
  // Remember to open time cuts in reader
  
  //Select the calorimeter
  TObjArray * pl = 0x0;
  AliVCaloCells* cells    = 0;
  if      (fCalorimeter == "PHOS" )
  {
    pl    = GetPHOSClusters();
    cells = GetPHOSCells();
  }
  else if (fCalorimeter == "EMCAL")
  {
    pl    = GetEMCALClusters();
    cells = GetEMCALCells();
  }
  
  if(!pl)
  {
    Info("MakeAnalysisFillAOD","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }
  
  AliVEvent  * event = GetReader()->GetInputEvent();
  AliESDEvent* esdEv = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodEv = dynamic_cast<AliAODEvent*> (event);
	
  //-------------------
  // N pile up vertices
  Int_t nVtxSPD = -1;
  Int_t nVtxTrk = -1;
	
  if      (esdEv)
  {
		nVtxSPD = esdEv->GetNumberOfPileupVerticesSPD();
		nVtxTrk = esdEv->GetNumberOfPileupVerticesTracks();
  }//ESD
  else if (aodEv)
  {
		nVtxSPD = aodEv->GetNumberOfPileupVerticesSPD();
		nVtxTrk = aodEv->GetNumberOfPileupVerticesTracks();
  }//AOD

  //-------------------
  // Loop on clusters
  Int_t nCaloClusters = pl->GetEntriesFast();

  if(GetDebug() > 0) printf("AliAnaClusterPileUp::MakeAnalysisFillAOD() - input %s cluster entries %d\n", fCalorimeter.Data(), nCaloClusters);
  //Init variables
  TLorentzVector mom;
  Int_t   idMax = 0;
  Float_t ptMax = 0;
  Float_t  tMax = 0;
  
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++)
  {
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));
    //printf("calo %d, %f\n",icalo,calo->E());
    
    if(!calo)  continue; // it should not happen, but just in case
    
    calo->GetMomentum(mom,GetVertex(0)) ;
  
    Float_t  ecluster  = mom.E();
    Float_t ptcluster  = mom.Pt();
    Float_t l0cluster  = calo->GetM02();
    Float_t etacluster = mom.Eta();
    Float_t phicluster = mom.Phi();
    if(phicluster < 0) phicluster+=TMath::TwoPi();
    Float_t tofcluster   = calo->GetTOF()*1.e9;
    
    Bool_t matched = IsTrackMatched(calo,GetReader()->GetInputEvent());

    //.......................................
    //If too small or big energy, skip it
    if(ecluster < GetMinEnergy() || ecluster > GetMaxEnergy() ) continue ;

    //.......................................
    if(calo->GetNCells() <= fNCellsCut && GetReader()->GetDataType() != AliCaloTrackReader::kMC) continue;
    
     //.......................................
    //Check acceptance selection
    if(IsFiducialCutOn())
    {
      Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter) ;
      if(! in ) continue;
    }

    // Select highest pt cluster passing the cuts
    if(ptcluster > ptMax && tofcluster < 30)
    {
      ptMax = ptcluster;
			tMax  = tofcluster;
      idMax = icalo;
    }
    
    //-------------------------------------
    // Cluster timing for different pile-up
    
    fhTimePtNoCut->Fill(ptcluster,tofcluster);
    if(GetReader()->IsPileUpFromSPD()) fhTimePtSPD->Fill(ptcluster,tofcluster);
    
    //----------------------------------------
    // correlate cluster and number of vertices
    
    fhPtNPileUpSPDVtx->Fill(ptcluster,nVtxSPD);
		fhPtNPileUpTrkVtx->Fill(ptcluster,nVtxTrk);
    
		if(TMath::Abs(tofcluster) < 30)
		{
			fhPtNPileUpSPDVtxTimeCut->Fill(ptcluster,nVtxSPD);
			fhPtNPileUpTrkVtxTimeCut->Fill(ptcluster,nVtxTrk);
		}
    
    if(tofcluster < 75 && tofcluster > -30)
    {
      fhPtNPileUpSPDVtxTimeCut2->Fill(ptcluster,nVtxSPD);
      fhPtNPileUpTrkVtxTimeCut2->Fill(ptcluster,nVtxTrk);
    }
    
    // Loop on the vertices arrays, correlate with timing
    // only for sufficiently large cluster energy
    if(ecluster > 8)
    {
      fhTimeNPileUpVertSPD  ->Fill(tofcluster,nVtxSPD);
      fhTimeNPileUpVertTrack->Fill(tofcluster,nVtxTrk);
      
      Int_t ncont = -1;
      Float_t z1 = -1, z2 = -1;
      Float_t diamZ = -1;
      for(Int_t iVert=0; iVert<nVtxSPD;iVert++)
      {
        if      (esdEv)
        {
          const AliESDVertex* pv=esdEv->GetPileupVertexSPD(iVert);
          ncont=pv->GetNContributors();
          z1 = esdEv->GetPrimaryVertexSPD()->GetZ();
          z2 = pv->GetZ();
          diamZ = esdEv->GetDiamondZ();
        }//ESD
        else if (aodEv)
        {
          AliAODVertex *pv=aodEv->GetVertex(iVert);
          if(pv->GetType()!=AliAODVertex::kPileupSPD) continue;
          ncont=pv->GetNContributors();
          z1=aodEv->GetPrimaryVertexSPD()->GetZ();
          z2=pv->GetZ();
          diamZ = aodEv->GetDiamondZ();
        }// AOD
        
        Double_t distZ  = TMath::Abs(z2-z1);
        diamZ  = TMath::Abs(z2-diamZ);
        
        fhTimeNPileUpVertContributors  ->Fill(tofcluster,ncont);
        fhTimePileUpMainVertexZDistance->Fill(tofcluster,distZ);
        fhTimePileUpMainVertexZDiamond ->Fill(tofcluster,diamZ);
        
      }// vertex loop
    }
    
    //------------------------------------
    // Eta-Phi cluster position depending on timing
    // Continue only for BC0
    if      (tofcluster > 28)
    {
      fhEtaPhiBCPlus ->Fill(etacluster,phicluster);
      if(GetReader()->IsPileUpFromSPD()) fhEtaPhiBCPlusPileUpSPD ->Fill(etacluster,phicluster);
      continue;
    }
    else if (tofcluster <-28)
    {
      fhEtaPhiBCMinus->Fill(etacluster,phicluster);
      if(GetReader()->IsPileUpFromSPD()) fhEtaPhiBCMinusPileUpSPD->Fill(etacluster,phicluster);
      continue ;
    }
    
    //--------------------------------------
    // Fill histograms for clusters in BC=0
    
    fhEtaPhiBC0->Fill(etacluster,phicluster); if(GetReader()->IsPileUpFromSPD()) fhEtaPhiBC0PileUpSPD    ->Fill(etacluster,phicluster);
    
    if(GetReader()->IsPileUpFromSPD())               {fhPtPileUp[0]->Fill(ptcluster); fhLambda0PileUp[0]->Fill(ptcluster,l0cluster); }
    if(GetReader()->IsPileUpFromEMCal())             {fhPtPileUp[1]->Fill(ptcluster); fhLambda0PileUp[1]->Fill(ptcluster,l0cluster); }
    if(GetReader()->IsPileUpFromSPDOrEMCal())        {fhPtPileUp[2]->Fill(ptcluster); fhLambda0PileUp[2]->Fill(ptcluster,l0cluster); }
    if(GetReader()->IsPileUpFromSPDAndEMCal())       {fhPtPileUp[3]->Fill(ptcluster); fhLambda0PileUp[3]->Fill(ptcluster,l0cluster); }
    if(GetReader()->IsPileUpFromSPDAndNotEMCal())    {fhPtPileUp[4]->Fill(ptcluster); fhLambda0PileUp[4]->Fill(ptcluster,l0cluster); }
    if(GetReader()->IsPileUpFromEMCalAndNotSPD())    {fhPtPileUp[5]->Fill(ptcluster); fhLambda0PileUp[5]->Fill(ptcluster,l0cluster); }
    if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) {fhPtPileUp[6]->Fill(ptcluster); fhLambda0PileUp[6]->Fill(ptcluster,l0cluster); }
    
    if(!matched)
    {
      if(GetReader()->IsPileUpFromSPD())               {fhPtNeutralPileUp[0]->Fill(ptcluster); fhLambda0NeutralPileUp[0]->Fill(ptcluster,l0cluster); }
      if(GetReader()->IsPileUpFromEMCal())             {fhPtNeutralPileUp[1]->Fill(ptcluster); fhLambda0NeutralPileUp[1]->Fill(ptcluster,l0cluster); }
      if(GetReader()->IsPileUpFromSPDOrEMCal())        {fhPtNeutralPileUp[2]->Fill(ptcluster); fhLambda0NeutralPileUp[2]->Fill(ptcluster,l0cluster); }
      if(GetReader()->IsPileUpFromSPDAndEMCal())       {fhPtNeutralPileUp[3]->Fill(ptcluster); fhLambda0NeutralPileUp[3]->Fill(ptcluster,l0cluster); }
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())    {fhPtNeutralPileUp[4]->Fill(ptcluster); fhLambda0NeutralPileUp[4]->Fill(ptcluster,l0cluster); }
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())    {fhPtNeutralPileUp[5]->Fill(ptcluster); fhLambda0NeutralPileUp[5]->Fill(ptcluster,l0cluster); }
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) {fhPtNeutralPileUp[6]->Fill(ptcluster); fhLambda0NeutralPileUp[6]->Fill(ptcluster,l0cluster); }
    }

    //----------------------------------------------------------------------------
    // Loop on cells inside cluster, max cell must be over 100 MeV and time in BC=0
    // Get the fraction of the cluster energy that carries the cell with highest
    // energy and its absId
    
    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, calo,maxCellFraction);
    
    Float_t clusterLongTimePt = 0;
    Float_t clusterOKTimePt   = 0;

    if(cells->GetCellAmplitude(absIdMax) < 0.1) continue ;
    
    for (Int_t ipos = 0; ipos < calo->GetNCells(); ipos++)
    {
      Int_t absId  = calo->GetCellsAbsId()[ipos];
      
      if( absId == absIdMax ) continue ;
      
      Double_t time  = cells->GetCellTime(absId);
      Float_t  amp   = cells->GetCellAmplitude(absId);
      Int_t    bc    = GetReader()->GetInputEvent()->GetBunchCrossNumber();
      GetCaloUtils()->GetEMCALRecoUtils()->AcceptCalibrateCell(absId,bc,amp,time,cells);
      time*=1e9;
      
      Float_t diff = (tofcluster-time);
      
      if(GetReader()->IsInTimeWindow(time,amp)) clusterOKTimePt   += amp;
      else                                      clusterLongTimePt += amp;
      
      if( cells->GetCellAmplitude(absIdMax) < 0.1 ) continue ;
      
      if(GetReader()->IsPileUpFromSPD())
      {
        fhClusterCellTimePileUp[0]->Fill(ptcluster, time);
        fhClusterTimeDiffPileUp[0]->Fill(ptcluster, diff);
        if(!matched) fhClusterTimeDiffNeutralPileUp[0]->Fill(ptcluster, diff);
      }
      
      if(GetReader()->IsPileUpFromEMCal())
      {
        fhClusterCellTimePileUp[1]->Fill(ptcluster, time);
        fhClusterTimeDiffPileUp[1]->Fill(ptcluster, diff);
        if(!matched) fhClusterTimeDiffNeutralPileUp[1]->Fill(ptcluster, diff);
      }
      
      if(GetReader()->IsPileUpFromSPDOrEMCal())
      {
        fhClusterCellTimePileUp[2]->Fill(ptcluster, time);
        fhClusterTimeDiffPileUp[2]->Fill(ptcluster, diff);
        if(!matched) fhClusterTimeDiffNeutralPileUp[2]->Fill(ptcluster, diff);
      }
      
      if(GetReader()->IsPileUpFromSPDAndEMCal())
      {
        fhClusterCellTimePileUp[3]->Fill(ptcluster, time);
        fhClusterTimeDiffPileUp[3]->Fill(ptcluster, diff);
        if(!matched) fhClusterTimeDiffNeutralPileUp[3]->Fill(ptcluster, diff);
      }
      
      if(GetReader()->IsPileUpFromSPDAndNotEMCal())
      {
        fhClusterCellTimePileUp[4]->Fill(ptcluster, time);
        fhClusterTimeDiffPileUp[4]->Fill(ptcluster, diff);
        if(!matched) fhClusterTimeDiffNeutralPileUp[4]->Fill(ptcluster, diff);
      }
      
      if(GetReader()->IsPileUpFromEMCalAndNotSPD())
      {
        fhClusterCellTimePileUp[5]->Fill(ptcluster, time);
        fhClusterTimeDiffPileUp[5]->Fill(ptcluster, diff);
        if(!matched) fhClusterTimeDiffNeutralPileUp[5]->Fill(ptcluster, diff);
      }
      
      if(GetReader()->IsPileUpFromNotSPDAndNotEMCal())
      {
        fhClusterCellTimePileUp[6]->Fill(ptcluster, time);
        fhClusterTimeDiffPileUp[6]->Fill(ptcluster, diff);
        if(!matched) fhClusterTimeDiffNeutralPileUp[6]->Fill(ptcluster, diff);
      }
    }//loop
    
    Float_t frac = 0;
    if(clusterLongTimePt+clusterOKTimePt > 0.001)
      frac = clusterLongTimePt/(clusterLongTimePt+clusterOKTimePt);
    //printf("E long %f, E OK %f, Fraction large time %f, E %f\n",clusterLongTimePt,clusterOKTimePt,frac,ptcluster);
    
    if(GetReader()->IsPileUpFromSPD())               fhClusterEFracLongTimePileUp[0]->Fill(ptcluster,frac);
    if(GetReader()->IsPileUpFromEMCal())             fhClusterEFracLongTimePileUp[1]->Fill(ptcluster,frac);
    if(GetReader()->IsPileUpFromSPDOrEMCal())        fhClusterEFracLongTimePileUp[2]->Fill(ptcluster,frac);
    if(GetReader()->IsPileUpFromSPDAndEMCal())       fhClusterEFracLongTimePileUp[3]->Fill(ptcluster,frac);
    if(GetReader()->IsPileUpFromSPDAndNotEMCal())    fhClusterEFracLongTimePileUp[4]->Fill(ptcluster,frac);
    if(GetReader()->IsPileUpFromEMCalAndNotSPD())    fhClusterEFracLongTimePileUp[5]->Fill(ptcluster,frac);
    if(GetReader()->IsPileUpFromNotSPDAndNotEMCal()) fhClusterEFracLongTimePileUp[6]->Fill(ptcluster,frac);
  }//loop
  
  //----------------------------------------------------------------------------------------------
  // Loop again on clusters to compare this max cluster t and the rest of the clusters, if E > 0.3
  Int_t n20  = 0;
  Int_t n40  = 0;
  Int_t n    = 0;
  Int_t nOK  = 0;
  
  for(Int_t icalo = 0; icalo < nCaloClusters; icalo++)
  {
	  AliVCluster * calo =  (AliVCluster*) (pl->At(icalo));
    
    if(!calo || calo->E() < 0.3 || icalo == idMax) continue;
    
    Float_t tdiff = TMath::Abs(tMax-calo->GetTOF()*1e9);
    n++;
    if(tdiff < 25) nOK++;
    else
    {
      n20++;
      if(tdiff > 40 ) n40++;
    }
  }
  
  // Check pile-up and fill histograms depending on the different cluster multiplicities
  if(GetReader()->IsPileUpFromSPD())
  {
    fhClusterMultSPDPileUp[0]->Fill(ptMax,n  );
    fhClusterMultSPDPileUp[1]->Fill(ptMax,nOK);
    fhClusterMultSPDPileUp[2]->Fill(ptMax,n20);
    fhClusterMultSPDPileUp[3]->Fill(ptMax,n40);
  }
  else
  {
    fhClusterMultNoPileUp[0]->Fill(ptMax,n  );
    fhClusterMultNoPileUp[1]->Fill(ptMax,nOK);
    fhClusterMultNoPileUp[2]->Fill(ptMax,n20);
    fhClusterMultNoPileUp[3]->Fill(ptMax,n40);
  }

  
  if(GetDebug() > 1) printf("AliAnaClusterPileUp::MakeAnalysisFillHistograms()  End fill histograms\n");
  
}



//__________________________________________________________________
void AliAnaClusterPileUp::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
  printf("    \n") ;
	
}
