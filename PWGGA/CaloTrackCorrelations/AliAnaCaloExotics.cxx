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

// --- ROOT system ---
#include <TObjArray.h>
#include <TDatabasePDG.h>
#include <TH3F.h>
#include <TObjString.h>

//---- AliRoot system ----
#include "AliAnaCaloExotics.h"
#include "AliCaloTrackReader.h"
#include "AliVCaloCells.h"
#include "AliFiducialCut.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliMCAnalysisUtils.h"
#include "TCustomBinning.h"

// --- Detectors --- 
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliAnaCaloExotics) ;
/// \endcond

//__________________________________________
/// Default Constructor. Initialize parameters.
/// Init histogram arrays to 0.
//__________________________________________
AliAnaCaloExotics::AliAnaCaloExotics() :
AliAnaCaloTrackCorrBaseClass(),  


fCellAmpMin(),                         fEMinForExo(0),
fFill1CellHisto(0),                    fConstantTimeShift(0),                 
fClusterMomentum(),         

// Histograms
fhExoticityEClus(0),                   fhExoticityEMaxCell(0),
fhExoticityEClusTrackMatch(0),         fhExoticity1Cell(0),

fhNCellsPerCluster(0),                 fhNCellsPerClusterTrackMatch(0),                             
fhNCellsPerClusterExo(0),              fhNCellsPerClusterExoTrackMatch(0),  
fhNCellsPerClusterM02(0),

fhEtaPhiExo(0),                        fhEtaPhiExoHighE(0),                    fhEtaPhi1Cell(0),
fhTimeEnergyExo(0),                    fhTimeEnergy1Cell(0),                
fhTimeDiffClusCellExo(0),              fhTimeDiffClusCellExoWide(0),           fhTimeDiffAmpClusCellExo(0),   
fhTimeEnergyM02(0),                    fhTimeDiffClusCellM02(0),               fhTimeDiffClusCellM02Wide(0),
fhM02EnergyNCell(0),                   fhM02EnergyExo(0),                     
fhM02EnergyExoZoomIn(0),               fhM20EnergyExoM02MinCut(0),                                                 

// Track matching vs exoticity
fhTrackMatchedDEtaNegExo(0),           fhTrackMatchedDPhiNegExo(0),            fhTrackMatchedDEtaDPhiNegExo(0),
fhTrackMatchedDEtaPosExo(0),           fhTrackMatchedDPhiPosExo(0),            fhTrackMatchedDEtaDPhiPosExo(0),
fhEOverPExo(0),

// Track matching of 1 cell clusters
fhTrackMatchedDEtaNeg1Cell(0),         fhTrackMatchedDPhiNeg1Cell(0),          fhTrackMatchedDEtaDPhiNeg1Cell(0),
fhTrackMatchedDEtaPos1Cell(0),         fhTrackMatchedDPhiPos1Cell(0),          fhTrackMatchedDEtaDPhiPos1Cell(0),
fhEOverP1Cell(0),

// Cells
fhCellExoAmp(0),                        fhCellExoAmpLowGain(0),                       
fhCellExoAmpTime(0),                    fhCellExoGrid(0),                      fhCellExoGridHighE(0)

{        
  AddToHistogramsName("AnaCaloExotic_");
  
  fCellAmpMin = 0.5;
  fEMinForExo = 6.0;
}

//
//____________________________________________________________
/// Fill histograms related to cells only.
/// \param cells: cells info list container
//____________________________________________________________
void AliAnaCaloExotics::CellHistograms(AliVCaloCells *cells)
{  
  Int_t ncells = cells->GetNumberOfCells();
  if( ncells <= 0 ) return;
  
  AliDebug(1,Form("%s cell entries %d", GetCalorimeterString().Data(), ncells));
  
  Int_t    icol   = -1, icolAbs = -1;
  Int_t    irow   = -1, irowAbs = -1;
  Int_t    iRCU   = -1;
  Float_t  amp    = 0.;
  Double_t time   = 0.;
  Int_t    id     = -1;
  Bool_t   highG  = kFALSE;
  Float_t  exoticity = -1000;
  
  Int_t    bc     = GetReader()->GetInputEvent()->GetBunchCrossNumber();
  
  for (Int_t iCell = 0; iCell < cells->GetNumberOfCells(); iCell++)
  {
    if ( cells->GetCellNumber(iCell) < 0 ||  cells->GetAmplitude(iCell) < fCellAmpMin ) continue; // CPV 
    
    AliDebug(2,Form("Cell : amp %f, absId %d", cells->GetAmplitude(iCell), cells->GetCellNumber(iCell)));
    
    Int_t nModule = GetModuleNumberCellIndexesAbsCaloMap(cells->GetCellNumber(iCell),GetCalorimeter(), 
                                                         icol   , irow, iRCU,
                                                         icolAbs, irowAbs    );
    
    AliDebug(2,Form("\t module %d, column %d (%d), row %d (%d)", nModule,icolAbs,icol,irowAbs,irow));
    
    amp     = cells->GetAmplitude(iCell);
    
    time    = cells->GetTime(iCell);
    time *= 1.0e9;
    time-=fConstantTimeShift;
    
    id      = cells->GetCellNumber(iCell);
    
    highG   = cells->GetCellHighGain(id);
    if ( IsDataMC() ) highG = kTRUE; // MC does not distinguish High and Low, put them all in high
    
    exoticity =  1-GetCaloUtils()->GetECross(id,cells,bc)/amp;
    
    // Fill histograms
    
    fhCellExoAmp    ->Fill(amp       , exoticity, GetEventWeight());
    fhCellExoAmpTime->Fill(amp , time, exoticity, GetEventWeight());
    
    if ( !highG )
      fhCellExoAmpLowGain->Fill(amp  , exoticity, GetEventWeight());
    
    if ( amp > fEMinForExo )
      fhCellExoGrid->Fill(icolAbs, irowAbs, exoticity, GetEventWeight());
    
    if ( amp > 100 && fEMinForExo < 100 )
      fhCellExoGridHighE->Fill(icolAbs, irowAbs, exoticity, GetEventWeight());

  } // Cell loop
  
}


//____________________________________________________________________________
/// Fill clusters related histograms, execute here the loop of clusters
/// apply basic selection cuts (track matching, goondess, exoticity, timing)
/// and the call to the different methods
/// filling different type of histograms:
/// * Basic cluster histograms for good or bad clusters
/// * Exotic cluster histograms
/// * Cells in cluster
/// * Invariant mass
/// * Matched clusters histograms
/// \param caloClusters: full list of clusters
/// \param cells: full list of cells
//____________________________________________________________________________
void AliAnaCaloExotics::ClusterHistograms(const TObjArray *caloClusters,
                                                AliVCaloCells* cells)
{
  Int_t  nCaloClusters = caloClusters->GetEntriesFast() ;
  Int_t  bc            = GetReader()->GetInputEvent()->GetBunchCrossNumber();
  
  // Get vertex, not needed.
  Double_t v[3] = {0,0,0};
  //GetReader()->GetVertex(v);
  
  AliDebug(1,Form("In %s there are %d clusters", GetCalorimeterString().Data(), nCaloClusters));
  
  // Loop over CaloClusters
  for(Int_t iclus = 0; iclus < nCaloClusters; iclus++)
  {
    AliDebug(1,Form("Cluster: %d/%d, data %d",iclus+1,nCaloClusters,GetReader()->GetDataType()));
    
    AliVCluster* clus =  (AliVCluster*) caloClusters->At(iclus);
    
    // Get cluster kinematics
    clus->GetMomentum(fClusterMomentum,v);
    
    Float_t en  = fClusterMomentum.E();
    Float_t eta = fClusterMomentum.Eta();
    Float_t phi = GetPhi(fClusterMomentum.Phi());
    
    // Check only certain regions
    Bool_t in = kTRUE;
    if ( IsFiducialCutOn() ) 
      in =  GetFiducialCut()->IsInFiducialCut(fClusterMomentum.Eta(),fClusterMomentum.Phi(),GetCalorimeter()) ;
    
    if ( !in )
    {
      AliDebug(1,Form("Remove cluster with phi %2.2f and eta %2.2f",phi*TMath::RadToDeg(),eta));
      continue;
    }
    
    Int_t    nCaloCellsPerCluster = clus->GetNCells();
    //Int_t  nLabel = clus->GetNLabels();
    //Int_t *labels = clus->GetLabels();
    
    // Cluster mathed with track?
    Bool_t matched = GetCaloPID()->IsTrackMatched(clus,GetCaloUtils(), GetReader()->GetInputEvent());
     
    // Get the fraction of the cluster energy that carries the cell with highest energy and its absId
    Float_t maxCellFraction = 0.;
    Int_t absIdMax = GetCaloUtils()->GetMaxEnergyCell(cells, clus,maxCellFraction);
    
    // Get time of max cell/cluster
    Double_t tmax  = cells->GetCellTime(absIdMax);
    tmax*=1.e9;
    tmax-=fConstantTimeShift;
    
    Float_t ampMax = cells->GetCellAmplitude(absIdMax);
    
    Float_t exoticity = 1-GetCaloUtils()->GetECross(absIdMax,cells,bc)/ampMax;
    
    Float_t m02  = clus->GetM02();
    Float_t m20  = clus->GetM20();
    
    AliDebug(1,Form("cluster: E %2.3f, F+ %2.3f, eta %2.3f, phi %2.3f, ncells %d,"
                    "match %d; cell max: id %d, en %2.3f, time %2.3f, m02 %2.2f",
                    en,exoticity,eta,phi*TMath::RadToDeg(),nCaloCellsPerCluster,
                    matched, absIdMax,ampMax,tmax,m02));  

    //
    // Fill histograms related to single cluster 
    //
    
    // Exoticity
    //
    if ( nCaloCellsPerCluster > 1 )
    {
      fhExoticityEClus   ->Fill(en    , exoticity, GetEventWeight());
      fhExoticityEMaxCell->Fill(ampMax, exoticity, GetEventWeight());
      
      if ( matched )
        fhExoticityEClusTrackMatch->Fill(en, exoticity, GetEventWeight());
    }
    else if ( fFill1CellHisto )
      fhExoticity1Cell->Fill(en, exoticity, GetEventWeight());
    
    // N cells per cluster
    //
    fhNCellsPerCluster               ->Fill(en, nCaloCellsPerCluster           , GetEventWeight());
    fhNCellsPerClusterExo            ->Fill(en, nCaloCellsPerCluster, exoticity, GetEventWeight());
    fhNCellsPerClusterM02            ->Fill(en, nCaloCellsPerCluster, m02      , GetEventWeight());
    
    if ( matched )
    {
      fhNCellsPerClusterTrackMatch   ->Fill(en, nCaloCellsPerCluster           , GetEventWeight());
      fhNCellsPerClusterExoTrackMatch->Fill(en, nCaloCellsPerCluster, exoticity, GetEventWeight());
    }
    
    // Acceptance
    //
    if ( en > fEMinForExo )
    {
      if ( nCaloCellsPerCluster > 1 )
      {
        fhEtaPhiExo  ->Fill(eta, phi, exoticity, GetEventWeight());
        if ( fEMinForExo < 100 && en > 100 )
          fhEtaPhiExoHighE->Fill(eta, phi, exoticity, GetEventWeight());
      }
      else if ( fFill1CellHisto )
        fhEtaPhi1Cell->Fill(eta, phi,            GetEventWeight());
    }
    
    // Timing
    //
    if ( nCaloCellsPerCluster > 1 )
    {
      fhTimeEnergyExo  ->Fill(en, tmax, exoticity, GetEventWeight());
      fhTimeEnergyM02  ->Fill(en, tmax, m02      , GetEventWeight());
    }
    else if ( fFill1CellHisto )
      fhTimeEnergy1Cell->Fill(en, tmax,            GetEventWeight());

    for (Int_t ipos = 0; ipos < nCaloCellsPerCluster; ipos++) 
    {
      Int_t   absId = clus->GetCellsAbsId()[ipos];  
      Float_t amp   = cells->GetCellAmplitude(absId);
      
      if ( absId == absIdMax || amp < 0.1 ) continue;
      
      Double_t time  = cells->GetCellTime(absId);
      time*=1.e9;
      time-=fConstantTimeShift;
      
      Float_t tdiff = tmax-time;
      
      fhTimeDiffClusCellExo->Fill(en , tdiff, exoticity, GetEventWeight());
      fhTimeDiffClusCellM02->Fill(en , tdiff, m02      , GetEventWeight());  
      
      fhTimeDiffClusCellExoWide->Fill(en , tdiff, exoticity, GetEventWeight());
      fhTimeDiffClusCellM02Wide->Fill(en , tdiff, m02      , GetEventWeight());
      
      if ( en > fEMinForExo )
        fhTimeDiffAmpClusCellExo->Fill(amp, tdiff, exoticity, GetEventWeight());
      
    } // Fill cell-cluster histogram loop
    
    // Shower shape
    if ( nCaloCellsPerCluster > 1 )
    {
      fhM02EnergyNCell->Fill(en, m02, nCaloCellsPerCluster, GetEventWeight());
      
      fhM02EnergyExo->Fill(en, m02, exoticity, GetEventWeight());
      
      fhM02EnergyExoZoomIn->Fill(en, m02, exoticity, GetEventWeight());
      
      if ( m02 > 0.1 )
        fhM20EnergyExoM02MinCut->Fill(en, m20, exoticity, GetEventWeight());
    }
    
    // Track matching
    //
    if ( matched )
    {
      AliVTrack *track = GetCaloUtils()->GetMatchedTrack(clus, GetReader()->GetInputEvent());
      
      if(!track) continue;
      
      Double_t tmom   = track->P();
      Double_t eOverP = en/tmom;

      Bool_t positive = kFALSE;
      if(track) positive = (track->Charge()>0);
      
      // Residuals
      Float_t deta  = clus->GetTrackDz();
      Float_t dphi  = clus->GetTrackDx();
      
      if ( nCaloCellsPerCluster > 1 )
      {
        fhEOverPExo->Fill(en, eOverP, exoticity, GetEventWeight());
        
        if ( positive )
        {
          fhTrackMatchedDEtaPosExo->Fill(en, deta, exoticity, GetEventWeight());
          fhTrackMatchedDPhiPosExo->Fill(en, dphi, exoticity, GetEventWeight());
          
          if ( en > fEMinForExo ) 
          {
            fhTrackMatchedDEtaDPhiPosExo->Fill(deta, dphi, exoticity, GetEventWeight());
          }
        }
        else 
        {
          fhTrackMatchedDEtaNegExo->Fill(en, deta, exoticity, GetEventWeight());
          fhTrackMatchedDPhiNegExo->Fill(en, dphi, exoticity, GetEventWeight());
          
          if ( en > fEMinForExo ) 
          {
            fhTrackMatchedDEtaDPhiNegExo->Fill(deta, dphi, exoticity, GetEventWeight());
          }
        }
      
      } // more than 1 cell
      else if ( fFill1CellHisto )
      {
        fhEOverP1Cell->Fill(en, eOverP, GetEventWeight());
        
        if ( positive )
        {
          fhTrackMatchedDEtaPos1Cell->Fill(en, deta, GetEventWeight());
          fhTrackMatchedDPhiPos1Cell->Fill(en, dphi, GetEventWeight());
          
          if ( en > fEMinForExo ) 
          {
            fhTrackMatchedDEtaDPhiPos1Cell->Fill(deta, dphi, GetEventWeight());
          }
        }
        else 
        {
          fhTrackMatchedDEtaNeg1Cell->Fill(en, deta, GetEventWeight());
          fhTrackMatchedDPhiNeg1Cell->Fill(en, dphi, GetEventWeight());
          
          if(en > fEMinForExo) 
          {
            fhTrackMatchedDEtaDPhiNeg1Cell->Fill(deta, dphi, GetEventWeight());
          }
        }
      } // 1 cell clusters
      
    } // matched
    
  } // Cluster loop
  
}


//_________________________________________________
/// Save parameters used for analysis in a string.
//_________________________________________________
TObjString * AliAnaCaloExotics::GetAnalysisCuts()
{  	
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaCaloExotics ---:") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Calorimeter: %s;",GetCalorimeterString().Data()) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"Cell Amplitude > %2.1f GeV;",fCellAmpMin) ;
  parList+=onePar ;

  snprintf(onePar,buffersize,"E min Exo > %2.1f GeV;",fEMinForExo) ;
  parList+=onePar ;
 
  snprintf(onePar,buffersize,"fill 1 cell histo: %d;",fFill1CellHisto) ;
  parList+=onePar ;
  
  //Get parameters set in base class.
  //parList += GetBaseParametersList() ;
  
  //Get parameters set in FiducialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
	
  return new TObjString(parList) ;
}

//___________________________________________________
/// Create histograms to be saved in output file and
/// store them in the output container.
//___________________________________________________
TList * AliAnaCaloExotics::GetCreateOutputObjects()
{
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ExoticHistos") ; 
  
  // Init the number of modules, set in the class AliCalorimeterUtils
  //
  InitCaloParameters(); // See AliCaloTrackCorrBaseClass
    
  // Histogram binning and ranges
  // 
  Int_t nptbins       = GetHistogramRanges()->GetHistoPtBins(); 	        
  Float_t ptmax       = GetHistogramRanges()->GetHistoPtMax();           
  Float_t ptmin       = GetHistogramRanges()->GetHistoPtMin();
  
  Int_t nphibins      = GetHistogramRanges()->GetHistoPhiBins();     	    
  Float_t phimax      = GetHistogramRanges()->GetHistoPhiMax();          
  Float_t phimin      = GetHistogramRanges()->GetHistoPhiMin();
  
  Int_t netabins      = GetHistogramRanges()->GetHistoEtaBins();          
  Float_t etamax      = GetHistogramRanges()->GetHistoEtaMax();          
  Float_t etamin      = GetHistogramRanges()->GetHistoEtaMin();	
  
  Int_t ntimebins     = GetHistogramRanges()->GetHistoTimeBins();         
  Float_t timemax     = GetHistogramRanges()->GetHistoTimeMax();         
  Float_t timemin     = GetHistogramRanges()->GetHistoTimeMin();       
 
  Int_t nceclbins     = GetHistogramRanges()->GetHistoNClusterCellBins(); 
  Int_t   nceclmax    = GetHistogramRanges()->GetHistoNClusterCellMax(); 
  Int_t   nceclmin    = GetHistogramRanges()->GetHistoNClusterCellMin(); 
  
  Int_t ssbins        = GetHistogramRanges()->GetHistoShowerShapeBins();  
  Float_t ssmax       = GetHistogramRanges()->GetHistoShowerShapeMax();  
  Float_t ssmin       = GetHistogramRanges()->GetHistoShowerShapeMin();
  
  Int_t tdbins        = GetHistogramRanges()->GetHistoDiffTimeBins() ;    
  Float_t tdmax       = GetHistogramRanges()->GetHistoDiffTimeMax();     
  Float_t tdmin       = GetHistogramRanges()->GetHistoDiffTimeMin();
  
  // TM residuals
  Int_t   nresetabins = GetHistogramRanges()->GetHistoTrackResidualEtaBins();
  Float_t resetamax   = GetHistogramRanges()->GetHistoTrackResidualEtaMax();
  Float_t resetamin   = GetHistogramRanges()->GetHistoTrackResidualEtaMin();
  Int_t   nresphibins = GetHistogramRanges()->GetHistoTrackResidualPhiBins();
  Float_t resphimax   = GetHistogramRanges()->GetHistoTrackResidualPhiMax();
  Float_t resphimin   = GetHistogramRanges()->GetHistoTrackResidualPhiMin();
  
  Int_t nPoverEbins   = GetHistogramRanges()->GetHistoEOverPBins();       
  Float_t eOverPmax   = GetHistogramRanges()->GetHistoEOverPMax();       
  Float_t eOverPmin   = GetHistogramRanges()->GetHistoEOverPMin();
  
  // Exoticity
  Int_t nexobins  = 200; Float_t exomin  = -1 ; Float_t exomax  = 1;
  // For TH3, reduced binning
  Int_t nexobinsS = 50 ; Float_t exominS = 0.5; Float_t exomaxS = 1;
  
  // Cell column-row histograms, see base class for data members setting
  //fNMaxColsFull+2,-1.5,fNMaxColsFull+0.5, fNMaxRowsFull+2,-1.5,fNMaxRowsFull+0.5
  Int_t   ncolcell    = fNMaxColsFull+2;
  Float_t colcellmin  = -1.5;
  Float_t colcellmax  = fNMaxColsFull+0.5;
  
  Int_t   nrowcell    = fNMaxRowsFullMax-fNMaxRowsFullMin+2;
  Float_t rowcellmin  = fNMaxRowsFullMin-1.5;
  Float_t rowcellmax  = fNMaxRowsFullMax+0.5;
  
  //
  // Init histograms
  //
  
  // Cluster Exoticity 2D
  //
  fhExoticityEClus = new TH2F 
  ("hExoticityEClus","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoticityEClus->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoticityEClus->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEClus);    
  
  fhExoticityEClusTrackMatch = new TH2F 
  ("hExoticityEClusTrackMatch","cell #it{F}_{+} vs #it{E}_{cluster}, #it{n}_{cluster}^{cell} > 1, track matched",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoticityEClusTrackMatch->SetXTitle("#it{E}_{cluster} (GeV) ");
  fhExoticityEClusTrackMatch->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEClusTrackMatch);    
  
  fhExoticityEMaxCell = new TH2F 
  ("hExoticityEMaxCell","cell #it{F}_{+} vs #it{E}_{cell}^{max}, #it{n}_{cluster}^{cell} > 1",
   nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
  fhExoticityEMaxCell->SetXTitle("#it{E}_{cell}^{max} (GeV) ");
  fhExoticityEMaxCell->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhExoticityEMaxCell);    
  
  if ( fFill1CellHisto )
  {
    fhExoticity1Cell = new TH2F 
    ("hExoticity1Cell","cell #it{F}_{+} vs #it{E}, #it{n}_{cluster}^{cell} = 1",
     nptbins,ptmin,ptmax, nexobins,exomin,exomax); 
    fhExoticity1Cell->SetXTitle("#it{E} (GeV) ");
    fhExoticity1Cell->SetYTitle("#it{F}_{+}");
    outputContainer->Add(fhExoticity1Cell);    
  }
  
  // N cells per cluster
  fhNCellsPerCluster  = new TH2F 
  ("hNCellsPerCluster","# cells per cluster vs #it{E}_{cluster}",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
  fhNCellsPerCluster->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerCluster->SetYTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsPerCluster);

  fhNCellsPerClusterTrackMatch  = new TH2F 
  ("hNCellsPerClusterTrackMatch","# cells per cluster vs #it{E}_{cluster}, track-matched",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax); 
  fhNCellsPerClusterTrackMatch->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterTrackMatch->SetYTitle("#it{n}_{cells}");
  outputContainer->Add(fhNCellsPerClusterTrackMatch);

  fhNCellsPerClusterExo  = new TH3F 
  ("hNCellsPerClusterExo","# cells per cluster vs #it{E}_{cluster} vs exoticity",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax,nexobinsS,exominS,exomaxS); 
  fhNCellsPerClusterExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterExo->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhNCellsPerClusterExo);
  
  fhNCellsPerClusterExoTrackMatch  = new TH3F 
  ("hNCellsPerClusterExoTrackMatch","# cells per cluster vs #it{E}_{cluster}, track-matched",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax, nexobinsS,exominS,exomaxS); 
  fhNCellsPerClusterExoTrackMatch->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterExoTrackMatch->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterExoTrackMatch->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhNCellsPerClusterExoTrackMatch);

  fhNCellsPerClusterM02  = new TH3F 
  ("hNCellsPerClusterM02","# cells per cluster vs #it{E}_{cluster} vs  #sigma^{2}_{long}",
   nptbins,ptmin,ptmax, nceclbins,nceclmin,nceclmax,100,0,0.5); 
  fhNCellsPerClusterM02->SetXTitle("#it{E}_{cluster} (GeV)");
  fhNCellsPerClusterM02->SetYTitle("#it{n}_{cells}");
  fhNCellsPerClusterM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhNCellsPerClusterM02);
  
  // Cluster acceptance
  fhEtaPhiExo  = new TH3F 
  ("hEtaPhiExo",
   Form("#eta vs #varphi vs #it{F}_{+}, #it{E}_{cluster}> %2.1f, #it{n}_{cells}>1",fEMinForExo),
   netabins,etamin,etamax,nphibins,phimin,phimax,nexobinsS,exominS,exomaxS); 
  fhEtaPhiExo->SetXTitle("#eta ");
  fhEtaPhiExo->SetYTitle("#varphi (rad)");
  fhEtaPhiExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhEtaPhiExo);

  if ( fEMinForExo < 100 )
  {
    fhEtaPhiExoHighE  = new TH3F 
    ("hEtaPhiExoHighE",
     "#eta vs #varphi vs #it{F}_{+}, #it{E}_{cluster}> 100, #it{n}_{cells}>1",
     netabins,etamin,etamax,nphibins,phimin,phimax,nexobinsS,exominS,exomaxS); 
    fhEtaPhiExoHighE->SetXTitle("#eta ");
    fhEtaPhiExoHighE->SetYTitle("#varphi (rad)");
    fhEtaPhiExoHighE->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhEtaPhiExoHighE);
  }
  
  if ( fFill1CellHisto )
  {
    fhEtaPhi1Cell  = new TH2F 
    ("hEtaPhi1Cell",
     Form("#eta vs #varphi, #it{E}_{cluster}> %2.1f, #it{n}_{cells}=1",fEMinForExo),
     netabins,etamin,etamax,nphibins,phimin,phimax); 
    fhEtaPhi1Cell->SetXTitle("#eta ");
    fhEtaPhi1Cell->SetYTitle("#varphi (rad)");
    outputContainer->Add(fhEtaPhi1Cell);
  }
  
  // Timing and energy
  fhTimeEnergyExo  = new TH3F 
  ("hTimeEnergyExo","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{F}_{+}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, ntimebins,timemin,timemax, nexobinsS,exominS,exomaxS); 
  fhTimeEnergyExo->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyExo->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeEnergyExo);
  
  if ( fFill1CellHisto )
  {
    fhTimeEnergy1Cell  = new TH2F 
    ("hTimeEnergy1Cell","#it{E}_{cluster} vs #it{t}_{cluster} vs #it{F}_{+}, #it{n}_{cells}=1",
     nptbins,ptmin,ptmax, ntimebins,timemin,timemax); 
    fhTimeEnergy1Cell->SetXTitle("#it{E} (GeV)");
    fhTimeEnergy1Cell->SetYTitle("#it{t}_{cluster} (ns)");
    outputContainer->Add(fhTimeEnergy1Cell);
  }
  
  fhTimeDiffClusCellExo  = new TH3F 
  ("hTimeDiffClusCellExo","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
  fhTimeDiffClusCellExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffClusCellExo);
 
  fhTimeDiffClusCellExoWide  = new TH3F 
  ("hTimeDiffClusCellExoWide","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, tdbins,-1000,1000, nexobinsS,exominS,exomaxS); 
  fhTimeDiffClusCellExoWide->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellExoWide->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellExoWide->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffClusCellExoWide);
  
  fhTimeDiffAmpClusCellExo  = new TH3F 
  ("hTimeDiffAmpClusCellExo",
   Form("#it{E}_{cell i} vs #it{t}_{cell max}-#it{t}_{cell i} vs #it{F}_{+}, #it{n}_{cells}>1"),
   nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, nexobinsS,exominS,exomaxS); 
  fhTimeDiffAmpClusCellExo->SetXTitle("#it{E}_{cell i} (GeV)");
  fhTimeDiffAmpClusCellExo->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffAmpClusCellExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhTimeDiffAmpClusCellExo);
  
  fhTimeEnergyM02  = new TH3F 
  ("hTimeEnergyM02","#it{E}_{cluster} vs #it{t}_{cluster} vs #sigma^{2}_{long}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, ntimebins,timemin,timemax, 100,0,0.5); 
  fhTimeEnergyM02->SetXTitle("#it{E} (GeV)");
  fhTimeEnergyM02->SetYTitle("#it{t}_{cluster} (ns)");
  fhTimeEnergyM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhTimeEnergyM02);
  
  fhTimeDiffClusCellM02  = new TH3F 
  ("hTimeDiffClusCellM02","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #sigma^{2}_{long}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, tdbins,tdmin,tdmax, 100,0,0.5); 
  fhTimeDiffClusCellM02->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellM02->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellM02->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhTimeDiffClusCellM02);
 
  fhTimeDiffClusCellM02Wide = new TH3F 
  ("hTimeDiffClusCellM02Wide","#it{E}_{cluster} vs #it{t}_{cell max}-#it{t}_{cell i} vs #sigma^{2}_{long}, #it{n}_{cells}>1",
   nptbins,ptmin,ptmax, tdbins,-1000,1000, 100,0,0.5); 
  fhTimeDiffClusCellM02Wide->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTimeDiffClusCellM02Wide->SetYTitle("#Delta #it{t}_{cell max-i} (ns)");
  fhTimeDiffClusCellM02Wide->SetZTitle("#sigma^{2}_{long}");
  outputContainer->Add(fhTimeDiffClusCellM02Wide);
  
  // Shower shape
  //
  fhM02EnergyNCell  = new TH3F 
  ("hM02EnergyNCell","shower shape, #sigma^{2}_{long} vs #it{E}_{cluster} vs #it{n}_{cells}",
   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax, nceclbins,nceclmin,nceclmax); 
  fhM02EnergyNCell->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM02EnergyNCell->SetYTitle("#sigma^{2}_{long}");
  fhM02EnergyNCell->SetZTitle("#it{n}_{cells}");
  outputContainer->Add(fhM02EnergyNCell); 
  
  fhM02EnergyExo  = new TH3F 
  ("hM02EnergyExo","shower shape, #sigma^{2}_{long} vs #it{E}_{cluster} vs #it{F}_{+}",
   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax, nexobinsS,exominS,exomaxS); 
  fhM02EnergyExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM02EnergyExo->SetYTitle("#sigma^{2}_{long}");
  fhM02EnergyExo->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhM02EnergyExo); 

  fhM02EnergyExoZoomIn  = new TH3F 
  ("hM02EnergyExoZoomIn","shower shape, #sigma^{2}_{long} vs #it{E}_{cluster} vs #it{F}_{+}",
   nptbins,ptmin,ptmax,100,0,0.5,200,0.8,1.); 
  fhM02EnergyExoZoomIn->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM02EnergyExoZoomIn->SetYTitle("#sigma^{2}_{long}");
  fhM02EnergyExoZoomIn->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhM02EnergyExoZoomIn); 
  
  fhM20EnergyExoM02MinCut  = new TH3F 
  ("hM20EnergyExoM02MinCut","shower shape, #sigma^{2}_{short} vs #it{E}_{cluster} vs #it{F}_{+}, #sigma^{2}_{long} > 0.1",
   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax/2, nexobinsS,exominS,exomaxS); 
  fhM20EnergyExoM02MinCut->SetXTitle("#it{E}_{cluster} (GeV)");
  fhM20EnergyExoM02MinCut->SetYTitle("#sigma^{2}_{short}");
  fhM20EnergyExoM02MinCut->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhM20EnergyExoM02MinCut);  
  
  // Track matching
  
  fhTrackMatchedDEtaNegExo  = new TH3F
  ("hTrackMatchedDEtaNegExo","d#eta of cluster-negative track vs cluster energy vs #it{F}_{+}",
   nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax, nexobinsS,exominS,exomaxS);
  fhTrackMatchedDEtaNegExo->SetYTitle("d#eta");
  fhTrackMatchedDEtaNegExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTrackMatchedDEtaNegExo->SetZTitle("#it{F}_{+}");

  fhTrackMatchedDPhiNegExo  = new TH3F
  ("hTrackMatchedDPhiNegExo","d#varphi of cluster-negative track vs cluster energy vs #it{F}_{+}",
   nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
  fhTrackMatchedDPhiNegExo->SetYTitle("d#varphi (rad)");
  fhTrackMatchedDPhiNegExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTrackMatchedDPhiNegExo->SetZTitle("#it{F}_{+}");

  fhTrackMatchedDEtaDPhiNegExo  = new TH3F
  ("hTrackMatchedDEtaDPhiNegExo",
   Form("d#eta vs d#varphi of cluster- negative track vs #it{F}_{+}, E > %2.1f GeV",fEMinForExo),
   nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
  fhTrackMatchedDEtaDPhiNegExo->SetYTitle("d#varphi (rad)");
  fhTrackMatchedDEtaDPhiNegExo->SetXTitle("d#eta");
  fhTrackMatchedDEtaDPhiNegExo->SetZTitle("#it{F}_{+}");

  fhTrackMatchedDEtaPosExo  = new TH3F
  ("hTrackMatchedDEtaPosExo","d#eta of cluster-positive track vs cluster energy vs #it{F}_{+}",
   nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax, nexobinsS,exominS,exomaxS);
  fhTrackMatchedDEtaPosExo->SetYTitle("d#eta");
  fhTrackMatchedDEtaPosExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTrackMatchedDEtaPosExo->SetZTitle("#it{F}_{+}");

  fhTrackMatchedDPhiPosExo  = new TH3F
  ("hTrackMatchedDPhiPosExo","d#varphi of cluster-positive track vs cluster energy vs #it{F}_{+}",
   nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
  fhTrackMatchedDPhiPosExo->SetYTitle("d#varphi (rad)");
  fhTrackMatchedDPhiPosExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhTrackMatchedDPhiNegExo->SetZTitle("#it{F}_{+}");

  fhTrackMatchedDEtaDPhiPosExo  = new TH3F
  ("hTrackMatchedDEtaDPhiPosExo",
   Form("d#eta vs d#varphi of cluster-positive track vs #it{F}_{+}, E > %2.1f GeV",fEMinForExo),
   nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax, nexobinsS,exominS,exomaxS);
  fhTrackMatchedDEtaDPhiPosExo->SetYTitle("d#varphi (rad)");
  fhTrackMatchedDEtaDPhiPosExo->SetXTitle("d#eta");
  fhTrackMatchedDEtaDPhiNegExo->SetZTitle("#it{F}_{+}");

  fhEOverPExo = new TH3F
  ("hEOverPExo",
   "Track matches #it{E}/#it{p} vs #it{F}_{+}",
   nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax, nexobinsS,exominS,exomaxS);
  fhEOverPExo->SetYTitle("#it{E}/#it{p}");
  fhEOverPExo->SetXTitle("#it{E}_{cluster} (GeV)");
  fhEOverPExo->SetZTitle("#it{F}_{+}");

  outputContainer->Add(fhTrackMatchedDEtaNegExo) ;
  outputContainer->Add(fhTrackMatchedDPhiNegExo) ;
  outputContainer->Add(fhTrackMatchedDEtaPosExo) ;
  outputContainer->Add(fhTrackMatchedDPhiPosExo) ;
  outputContainer->Add(fhTrackMatchedDEtaDPhiNegExo) ;
  outputContainer->Add(fhTrackMatchedDEtaDPhiPosExo) ;
  outputContainer->Add(fhEOverPExo);
  
  if ( fFill1CellHisto )
  {
    fhTrackMatchedDEtaNeg1Cell  = new TH2F
    ("hTrackMatchedDEtaNeg1Cell",
     "d#eta of cluster-negative track vs cluster energy, #it{n}_{cell}=1",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEtaNeg1Cell->SetYTitle("d#eta");
    fhTrackMatchedDEtaNeg1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
    
    fhTrackMatchedDPhiNeg1Cell  = new TH2F
    ("hTrackMatchedDPhiNeg1Cell",
     "d#varphi of cluster-negative track vs cluster energy, #it{n}_{cell}=1",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhiNeg1Cell->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDPhiNeg1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhiNeg1Cell  = new TH2F
    ("hTrackMatchedDEtaDPhiNeg1Cell",
     Form("d#eta vs d#varphi of cluster-negative track, E > %2.2f, #it{n}_{cell}=1",fEMinForExo),
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhiNeg1Cell->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDEtaDPhiNeg1Cell->SetXTitle("d#eta");
    
    fhTrackMatchedDEtaPos1Cell  = new TH2F
    ("hTrackMatchedDEtaPos1Cell",
     "d#eta of cluster-positive track vs cluster energy, #it{n}_{cell}=1",
     nptbins,ptmin,ptmax,nresetabins,resetamin,resetamax);
    fhTrackMatchedDEtaPos1Cell->SetYTitle("d#eta");
    fhTrackMatchedDEtaPos1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
    
    fhTrackMatchedDPhiPos1Cell  = new TH2F
    ("hTrackMatchedDPhiPos1Cell",
     "d#varphi of cluster-positive track vs cluster energy, #it{n}_{cell}=1",
     nptbins,ptmin,ptmax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDPhiPos1Cell->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDPhiPos1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
    
    fhTrackMatchedDEtaDPhiPos1Cell  = new TH2F
    ("hTrackMatchedDEtaDPhiPos1Cell",
     Form("d#eta vs d#varphi of cluster-positive track, E > %2.2f, #it{n}_{cell}=1",fEMinForExo),
     nresetabins,resetamin,resetamax,nresphibins,resphimin,resphimax);
    fhTrackMatchedDEtaDPhiPos1Cell->SetYTitle("d#varphi (rad)");
    fhTrackMatchedDEtaDPhiPos1Cell->SetXTitle("d#eta");
    
    fhEOverP1Cell = new TH2F
    ("hEOverP1Cell",
     "Track matches #it{E}/#it{p}, #it{n}_{cell}=1",
     nptbins,ptmin,ptmax, nPoverEbins,eOverPmin,eOverPmax);
    fhEOverP1Cell->SetYTitle("#it{E}/#it{p}");
    fhEOverP1Cell->SetXTitle("#it{E}_{cluster} (GeV)");
    
    outputContainer->Add(fhTrackMatchedDEtaNeg1Cell) ;
    outputContainer->Add(fhTrackMatchedDPhiNeg1Cell) ;
    outputContainer->Add(fhTrackMatchedDEtaPos1Cell) ;
    outputContainer->Add(fhTrackMatchedDPhiPos1Cell) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiNeg1Cell) ;
    outputContainer->Add(fhTrackMatchedDEtaDPhiPos1Cell) ;
    outputContainer->Add(fhEOverP1Cell);
  }
  
  // Calorimeter cells
  //
  fhCellExoAmp     = new TH2F 
  ("hCellExoAmp","cell #it{F}_{+} vs #it{E}_{cell}",
   nptbins,ptmin,ptmax/2, nexobins,exomin,exomax); 
  fhCellExoAmp->SetXTitle("#it{E}_{cell} (GeV) ");
  fhCellExoAmp->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhCellExoAmp);    

  fhCellExoAmpLowGain     = new TH2F 
  ("hCellExoAmpLowGain","Cell #it{F}_{+} vs #it{E}_{cell}, low gain",
   nptbins,ptmin,ptmax/2, nexobins,exomin,exomax); 
  fhCellExoAmpLowGain->SetXTitle("#it{E}_{cell} (GeV) ");
  fhCellExoAmpLowGain->SetYTitle("#it{F}_{+}");
  outputContainer->Add(fhCellExoAmpLowGain);    
  
  fhCellExoAmpTime = new TH3F 
  ("hCellExoAmpTime","Cell #it{F}_{+} vs #it{E}_{cell} vs time",
   nptbins,ptmin,ptmax/2, ntimebins,timemin,timemax, nexobinsS,exominS,exomaxS); 
  fhCellExoAmpTime->SetXTitle("#it{E}_{cell} (GeV) ");
  fhCellExoAmpTime->SetYTitle("#it{t}_{cell} (ns)");
  fhCellExoAmpTime->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellExoAmpTime);    
  
  fhCellExoGrid    = new TH3F 
  ("hCellExoGrid",
   Form("Cell hits row-column vs #it{F}_{+} for #it{E}_{cell} > %2.1f",fEMinForExo), 
   ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, nexobinsS,exominS,exomaxS); 
  fhCellExoGrid->SetYTitle("row (phi direction)");
  fhCellExoGrid->SetXTitle("column (eta direction)");
  fhCellExoGrid->SetZTitle("#it{F}_{+}");
  outputContainer->Add(fhCellExoGrid);
    
  if ( fEMinForExo < 100 )
  {
    fhCellExoGridHighE    = new TH3F 
    ("hCellExoGridHighE",
     "Cell hits row-column vs #it{F}_{+} for #it{E}_{cell} > 100", 
    ncolcell,colcellmin,colcellmax,nrowcell,rowcellmin,rowcellmax, nexobinsS,exominS,exomaxS); 
    fhCellExoGridHighE->SetYTitle("row (phi direction)");
    fhCellExoGridHighE->SetXTitle("column (eta direction)");
    fhCellExoGridHighE->SetZTitle("#it{F}_{+}");
    outputContainer->Add(fhCellExoGridHighE);
  }
  
  //  for(Int_t i = 0; i < outputContainer->GetEntries() ; i++)
  //    printf("i=%d, name= %s\n",i,outputContainer->At(i)->GetName());
  
  return outputContainer;
}


//______________________________
/// Check if the calorimeter setting is ok, if not abort.
//______________________________
void AliAnaCaloExotics::Init()
{
  if(GetCalorimeter() != kPHOS && GetCalorimeter() !=kEMCAL)
    AliFatal(Form("Wrong calorimeter name <%s>", GetCalorimeterString().Data()));
  
  if(GetReader()->GetDataType()== AliCaloTrackReader::kMC)
    AliFatal("Analysis of reconstructed data, MC reader not aplicable");
}


//_________________________________________________________
/// Print some relevant parameters set for the analysis.
//_________________________________________________________
void AliAnaCaloExotics::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaCaloTrackCorrBaseClass::Print(" ");
  
  printf("Select Calorimeter %s \n",GetCalorimeterString().Data());
  printf("Min Amplitude : %2.1f GeV/c\n", fCellAmpMin) ;
  printf("Min Energy for exotic : %2.1f GeV/c\n", fEMinForExo) ;
  printf("Fill 1 cell histo : %d GeV/c\n", fFill1CellHisto) ;
}


//_____________________________________________________
/// Main task method, call all the methods filling QA histograms.
//_____________________________________________________
void  AliAnaCaloExotics::MakeAnalysisFillHistograms()
{
  AliDebug(1,"Start");

  // Get List with CaloClusters , calo Cells, init min amplitude
  TObjArray     * caloClusters = NULL;
  AliVCaloCells * cells        = 0x0;
  
  if      (GetCalorimeter() == kPHOS)
  {
    caloClusters = GetPHOSClusters();
    cells        = GetPHOSCells();
  }
  else if (GetCalorimeter() == kEMCAL)
  {
    caloClusters = GetEMCALClusters();
    cells        = GetEMCALCells();
  }
  
  if( !caloClusters || !cells )
  {
    AliWarning(Form("AliAnaCaloExotics::MakeAnalysisFillHistograms() - No CaloClusters or CaloCells available"));
    return; 
  }

  if(caloClusters->GetEntriesFast() == 0) return ;
  
  //printf("Exotic Task: N cells %d, N clusters %d \n",cells->GetNumberOfCells(),caloClusters->GetEntriesFast());
  
  // Clusters
  ClusterHistograms(caloClusters,cells);
  
  // Cells  
  CellHistograms(cells);
  
  AliDebug(1,"End");
}

