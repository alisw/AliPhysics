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
// Split clusters with some criteria and calculate invariant mass
// to identify them as pi0 or conversion
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble)  
//_________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TList.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TH3F.h>
//#include "Riostream.h"

// --- Analysis system --- 
#include "AliAnaInsideClusterInvariantMass.h" 
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"
#include "AliStack.h"
#include "AliFiducialCut.h"
#include "TParticle.h"
#include "AliVCluster.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliEMCALGeoParams.h"

ClassImp(AliAnaInsideClusterInvariantMass)
  
//__________________________________________________________________
AliAnaInsideClusterInvariantMass::AliAnaInsideClusterInvariantMass() : 
  AliAnaPartCorrBaseClass(),  
  fCalorimeter(""),  
  fM02Cut(0),
  fMinNCells(0)
{
  //default ctor
  
  // Init array of histograms
  for(Int_t i = 0; i < 7; i++){
    fhMassNLocMax1[i]  = 0;
    fhMassNLocMax2[i]  = 0;
    fhMassNLocMaxN[i]  = 0;
    fhNLocMax[i]       = 0;
    fhNLocMaxM02Cut[i] = 0;
    fhM02NLocMax1[i]   = 0;
    fhM02NLocMax2[i]   = 0;
    fhM02NLocMaxN[i]   = 0;
    fhNCellNLocMax1[i] = 0;
    fhNCellNLocMax2[i] = 0;
    fhNCellNLocMaxN[i] = 0;
    fhM02Pi0[i]        = 0;
    fhM02Eta[i]        = 0;
    fhM02Con[i]        = 0;
    fhInvMassAllCells[i]=0;
  }
   
  InitParameters();
  
}

//_______________________________________________________________
TObjString *  AliAnaInsideClusterInvariantMass::GetAnalysisCuts()
{	
	//Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliAnaInsideClusterInvariantMass ---\n") ;
  parList+=onePar ;	
  
  snprintf(onePar,buffersize,"Calorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fM02Cut =%f \n"   ,fM02Cut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fMinNCells =%f \n",fMinNCells) ;
  parList+=onePar ;  
  
  return new TObjString(parList) ;
  
}

//____________________________________________________________________________________________________
Bool_t AliAnaInsideClusterInvariantMass::AreNeighbours( const Int_t absId1, const Int_t absId2 ) const
{
  // Tells if (true) or not (false) two digits are neighbours
  // A neighbour is defined as being two digits which share a corner
	
  Bool_t areNeighbours = kFALSE ;
  Int_t nSupMod =0, nModule =0, nIphi =0, nIeta =0;
  Int_t nSupMod1=0, nModule1=0, nIphi1=0, nIeta1=0;
  Int_t relid1[2],  relid2[2] ;
  Int_t rowdiff=0,  coldiff=0;
  
  areNeighbours = kFALSE ;
  
  GetEMCALGeometry()->GetCellIndex(absId1, nSupMod,nModule,nIphi,nIeta);
  GetEMCALGeometry()->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, relid1[0],relid1[1]);
  
  GetEMCALGeometry()->GetCellIndex(absId2, nSupMod1,nModule1,nIphi1,nIeta1);
  GetEMCALGeometry()->GetCellPhiEtaIndexInSModule(nSupMod1,nModule1,nIphi1,nIeta1, relid2[0],relid2[1]);
  
  // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
  // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
  if(nSupMod1!=nSupMod){ 
    if(nSupMod1%2) relid1[1]+=AliEMCALGeoParams::fgkEMCALCols;
    else           relid2[1]+=AliEMCALGeoParams::fgkEMCALCols;
  }
	
  rowdiff = TMath::Abs( relid1[0] - relid2[0] ) ;  
  coldiff = TMath::Abs( relid1[1] - relid2[1] ) ;  
  
  if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0)) 
    areNeighbours = kTRUE ;
  
  return areNeighbours;
}

//_____________________________________________________________________________________
TLorentzVector AliAnaInsideClusterInvariantMass::GetCellMomentum(const Int_t absId,
                                                                 AliVCaloCells * cells)
{

  // Cell momentum calculation for invariant mass
  
  Double_t cellpos[] = {0, 0, 0};
  GetEMCALGeometry()->GetGlobal(absId, cellpos);
  
  if(GetVertex(0)){//calculate direction from vertex
    cellpos[0]-=GetVertex(0)[0];
    cellpos[1]-=GetVertex(0)[1];
    cellpos[2]-=GetVertex(0)[2];  
  }
  
  Double_t r = TMath::Sqrt(cellpos[0]*cellpos[0]+cellpos[1]*cellpos[1]+cellpos[2]*cellpos[2] ) ; 
  
  Float_t en = cells->GetCellAmplitude(absId);
  RecalibrateCellAmplitude(en,absId);  
  TLorentzVector cellMom ;   
  cellMom.SetPxPyPzE( en*cellpos[0]/r,  en*cellpos[1]/r, en*cellpos[2]/r,  en) ;   

  return cellMom;
  
}

//________________________________________________________________
TList * AliAnaInsideClusterInvariantMass::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("InsideClusterHistos") ; 
  
  Int_t nptbins  = GetHistoPtBins();           Float_t ptmax  = GetHistoPtMax();           Float_t ptmin  = GetHistoPtMin();
  Int_t ssbins   = GetHistoShowerShapeBins();  Float_t ssmax  = GetHistoShowerShapeMax();  Float_t ssmin  = GetHistoShowerShapeMin();
  Int_t mbins    = GetHistoMassBins();         Float_t mmax   = GetHistoMassMax();         Float_t mmin   = GetHistoMassMin();
  Int_t ncbins   = GetHistoNClusterCellBins(); Int_t   ncmax  = GetHistoNClusterCellMax(); Int_t   ncmin  = GetHistoNClusterCellMin(); 

  TString ptype[] ={"","#gamma","#gamma->e^{#pm}","#pi^{0}","#eta","e^{#pm}", "hadron"}; 
  TString pname[] ={"","Photon","Conversion",     "Pi0",    "Eta", "Electron","Hadron"};
  
  Int_t n = 1;
  
  if(IsDataMC()) n = 7;
  
  for(Int_t i = 0; i < n; i++){ 

    fhInvMassAllCells[i]  = new TH2F(Form("hInvMassAllCells%s",pname[i].Data()),
                                  Form("Invariant mass of all cells %s",ptype[i].Data()),
                                  nptbins,ptmin,ptmax,mbins,mmin,mmax); 
    fhInvMassAllCells[i]->SetYTitle("M (MeV/c^2)");
    fhInvMassAllCells[i]->SetXTitle("E (GeV)");
    outputContainer->Add(fhInvMassAllCells[i]) ;   
    
    
    fhMassNLocMax1[i]  = new TH2F(Form("hMassNLocMax1%s",pname[i].Data()),
                                  Form("Invariant mass of 2 highest energy cells %s",ptype[i].Data()),
                                  nptbins,ptmin,ptmax,mbins,mmin,mmax); 
    fhMassNLocMax1[i]->SetYTitle("M (MeV/c^2)");
    fhMassNLocMax1[i]->SetXTitle("E (GeV)");
    outputContainer->Add(fhMassNLocMax1[i]) ;   
    
    fhMassNLocMax2[i]  = new TH2F(Form("hMassNLocMax2%s",pname[i].Data()),
                                  Form("Invariant mass of 2 local maxima cells %s",ptype[i].Data()),
                                  nptbins,ptmin,ptmax,mbins,mmin,mmax); 
    fhMassNLocMax2[i]->SetYTitle("M (MeV/c^2)");
    fhMassNLocMax2[i]->SetXTitle("E (GeV)");
    outputContainer->Add(fhMassNLocMax2[i]) ;   

    fhMassNLocMaxN[i]  = new TH2F(Form("hMassNLocMaxN%s",pname[i].Data()),
                                  Form("Invariant mass of 2 local maxima cells %s",ptype[i].Data()),
                                  nptbins,ptmin,ptmax,mbins,mmin,mmax); 
    fhMassNLocMaxN[i]->SetYTitle("M (MeV/c^2)");
    fhMassNLocMax2[i]->SetXTitle("E (GeV)");
    outputContainer->Add(fhMassNLocMaxN[i]) ;   
    
    
    fhNLocMax[i]     = new TH2F(Form("hNLocMax%s",pname[i].Data()),
                             Form("Number of local maxima in cluster %s",ptype[i].Data()),
                             nptbins,ptmin,ptmax,5,0,5); 
    fhNLocMax[i]   ->SetYTitle("N maxima");
    fhNLocMax[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhNLocMax[i]) ; 

    fhNLocMaxM02Cut[i] = new TH2F(Form("hNLocMaxM02Cut%s",pname[i].Data()),
                              Form("Number of local maxima in cluster %s for M02 > %2.2f",ptype[i].Data(),fM02Cut),
                              nptbins,ptmin,ptmax,5,0,5); 
    fhNLocMaxM02Cut[i]->SetYTitle("N maxima");
    fhNLocMaxM02Cut[i]->SetXTitle("E (GeV)");
    outputContainer->Add(fhNLocMaxM02Cut[i]) ; 
    
    
    fhM02NLocMax1[i]     = new TH2F(Form("hM02NLocMax1%s",pname[i].Data()),
                                     Form("#lambda_{0}^{2} vs E for N max  = 1 %s",ptype[i].Data()),
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhM02NLocMax1[i]   ->SetYTitle("#lambda_{0}^{2}");
    fhM02NLocMax1[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhM02NLocMax1[i]) ; 
    
    fhM02NLocMax2[i]     = new TH2F(Form("hM02NLocMax2%s",pname[i].Data()),
                                     Form("#lambda_{0}^{2} vs E for N max  = 2 %s",ptype[i].Data()),
                                     nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhM02NLocMax2[i]   ->SetYTitle("#lambda_{0}^{2}");
    fhM02NLocMax2[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhM02NLocMax2[i]) ; 
    
    
    fhM02NLocMaxN[i]    = new TH2F(Form("hM02NLocMaxN%s",pname[i].Data()),
                                   Form("#lambda_{0}^{2} vs E for N max  > 2 %s",ptype[i].Data()),
                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhM02NLocMaxN[i]   ->SetYTitle("#lambda_{0}^{2}");
    fhM02NLocMaxN[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhM02NLocMaxN[i]) ; 

    
    fhNCellNLocMax1[i]  = new TH2F(Form("hNCellNLocMax1%s",pname[i].Data()),
                                   Form("#lambda_{0}^{2} vs E for N max  = 1 %s",ptype[i].Data()),
                                   nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
    fhNCellNLocMax1[i] ->SetYTitle("N cells");
    fhNCellNLocMax1[i] ->SetXTitle("E (GeV)");
    outputContainer->Add(fhNCellNLocMax1[i]) ; 
    
    fhNCellNLocMax2[i]     = new TH2F(Form("hNCellNLocMax2%s",pname[i].Data()),
                                    Form("#lambda_{0}^{2} vs E for N max  = 2 %s",ptype[i].Data()),
                                    nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
    fhNCellNLocMax2[i]   ->SetYTitle("N cells");
    fhNCellNLocMax2[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhNCellNLocMax2[i]) ; 
    
    
    fhNCellNLocMaxN[i]     = new TH2F(Form("hNCellNLocMaxN%s",pname[i].Data()),
                                    Form("#lambda_{0}^{2} vs E for N max  > 2 %s",ptype[i].Data()),
                                    nptbins,ptmin,ptmax,ncbins,ncmin,ncmax); 
    fhNCellNLocMaxN[i]   ->SetYTitle("N cells");
    fhNCellNLocMaxN[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhNCellNLocMaxN[i]) ;
    
    
    fhM02Pi0[i]     = new TH2F(Form("hM02Pi0%s",pname[i].Data()),
                                    Form("#lambda_{0}^{2} vs E ffor mass [0.1-0.2] MeV/c2 %s",ptype[i].Data()),
                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhM02Pi0[i]   ->SetYTitle("#lambda_{0}^{2}");
    fhM02Pi0[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhM02Pi0[i]) ; 
    
    fhM02Eta[i]     = new TH2F(Form("hM02Eta%s",pname[i].Data()),
                                    Form("#lambda_{0}^{2} vs E for mass [0.4-0.7] MeV/c2, %s",ptype[i].Data()),
                                    nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhM02Eta[i]   ->SetYTitle("#lambda_{0}^{2}");
    fhM02Eta[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhM02Eta[i]) ; 
    
    
    fhM02Con[i]    = new TH2F(Form("hM02Con%s",pname[i].Data()),
                                   Form("#lambda_{0}^{2} vs E for mass < 0.5 MeV/c2, %s",ptype[i].Data()),
                                   nptbins,ptmin,ptmax,ssbins,ssmin,ssmax); 
    fhM02Con[i]   ->SetYTitle("#lambda_{0}^{2}");
    fhM02Con[i]   ->SetXTitle("E (GeV)");
    outputContainer->Add(fhM02Con[i]) ; 
    
  }
  
  return outputContainer ;
  
}

//________________________________________________________________________________________________________
Int_t AliAnaInsideClusterInvariantMass::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
                                                               Int_t *absIdList,     Float_t *maxEList) 
{
  // Find local maxima in cluster
    
  Float_t locMaxCut = 0; // not used for the moment
  
  Int_t iDigitN = 0 ;
  Int_t iDigit  = 0 ;
  Int_t absId1 = -1 ;
  Int_t absId2 = -1 ;
  const Int_t nCells = cluster->GetNCells();
  //printf("cluster :");
  for(iDigit = 0; iDigit < nCells ; iDigit++){
    absIdList[iDigit] = cluster->GetCellsAbsId()[iDigit]  ; 
    
    //Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
    //RecalibrateCellAmplitude(en,absIdList[iDigit]);  
    //Int_t icol = -1, irow = -1, iRCU = -1;
    //Int_t sm = GetCaloUtils()->GetModuleNumberCellIndexes(absIdList[iDigit], fCalorimeter, icol, irow, iRCU) ;

    //printf("\t cell %d, id %d, sm %d, col %d, row %d, e %f\n", iDigit, absIdList[iDigit], sm, icol, irow, en );
  }
  
  for(iDigit = 0 ; iDigit < nCells; iDigit++) {   
    if(absIdList[iDigit]>=0) {
      
      absId1 = absIdList[iDigit] ;
      Float_t en1 = cells->GetCellAmplitude(absId1);
      RecalibrateCellAmplitude(en1,absId1);  
      
      for(iDigitN = 0; iDigitN < nCells; iDigitN++) {	
        absId2 = absIdList[iDigitN] ;     

        Float_t en2 = cells->GetCellAmplitude(absId2);
        RecalibrateCellAmplitude(en2,absId2);
        
        if ( AreNeighbours(absId1, absId2) ) {
          
          if (en1 > en2 ) {    
            absIdList[iDigitN] = -1 ;
            // but may be digit too is not local max ?
            if(en1 < en2 + locMaxCut) 
              absIdList[iDigit] = -1 ;
          }
          else {
            absIdList[iDigit] = -1 ;
            // but may be digitN too is not local max ?
            if(en1 > en2 - locMaxCut) 
              absIdList[iDigitN] = -1 ; 
          } 
        } // if Areneighbours
      } // while digitN
    } // slot not empty
  } // while digit
  
  iDigitN = 0 ;
  for(iDigit = 0; iDigit < nCells; iDigit++) { 
    if(absIdList[iDigit]>=0 ){
      absIdList[iDigitN] = absIdList[iDigit] ;
      Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
      RecalibrateCellAmplitude(en,absIdList[iDigit]);  
      if(en < 0.1) continue; // Maxima only with seed energy at least
      maxEList[iDigitN] = en ;
      //printf("Local max %d, id %d, en %f\n", iDigit,absIdList[iDigitN],en);
      iDigitN++ ; 
    }
  }
  
  return iDigitN ;
  
}


//___________________________________________
void AliAnaInsideClusterInvariantMass::Init()
{ 
  //Init
  //Do some checks
  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn() && NewOutputAOD()){
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn() && NewOutputAOD()){
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
    abort();
  }
  
  if( GetReader()->GetDataType() == AliCaloTrackReader::kMC ){
    printf("AliAnaInsideClusterInvariantMass::Init() - !!STOP: You want to use pure MC data!!\n");
    abort();
    
  }
  
}

//_____________________________________________________
void AliAnaInsideClusterInvariantMass::InitParameters()
{
  //Initialize the parameters of the analysis.  
  AddToHistogramsName("AnaPi0InsideClusterInvariantMass_");
  
  fCalorimeter = "EMCAL" ;
  fM02Cut      = 0.26 ;
  fMinNCells   = 4 ;
}


//__________________________________________________________________
void  AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() 
{
  //Search for pi0 in fCalorimeter with shower shape analysis 
  
  TObjArray * pl       = 0x0; 
  AliVCaloCells* cells = 0x0;

  //Select the Calorimeter of the photon
  if(fCalorimeter == "PHOS"){
    pl    = GetPHOSClusters();
    cells = GetPHOSCells();
  }
  else if (fCalorimeter == "EMCAL"){
    pl    = GetEMCALClusters();
    cells = GetEMCALCells();
  }
  
  if(!pl || !cells) {
    Info("MakeAnalysisFillHistograms","TObjArray with %s clusters is NULL!\n",fCalorimeter.Data());
    return;
  }  
  
	if(fCalorimeter == "PHOS") return; // Not implemented for PHOS yet

  for(Int_t icluster = 0; icluster < pl->GetEntriesFast(); icluster++){
    AliVCluster * cluster = (AliVCluster*) (pl->At(icluster));	

    // Study clusters with large shape parameter
    Float_t en = cluster->E();
    Float_t l0 = cluster->GetM02();
    Int_t   nc = cluster->GetNCells();
    
    //If too small or big E or low number of cells, skip it
    if( ( en < GetMinEnergy() || en > GetMaxEnergy() ) && nc < fMinNCells) continue ; 
  
    Int_t    absId1    = -1; Int_t absId2 = -1;
    Int_t   *absIdList = new Int_t  [nc]; 
    Float_t *maxEList  = new Float_t[nc]; 
    Int_t    nMax      = GetNumberOfLocalMaxima(cluster, cells, absIdList, maxEList) ;
    
    if (nMax <= 0) {
      printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - No local maximum found!");
      delete absIdList ;
      delete maxEList  ;
      return;
    }
    
    fhNLocMax[0]->Fill(en,nMax);
    
    if     ( nMax == 1  ) { fhM02NLocMax1[0]->Fill(en,l0) ; fhNCellNLocMax1[0]->Fill(en,nc) ; }
    else if( nMax == 2  ) { fhM02NLocMax2[0]->Fill(en,l0) ; fhNCellNLocMax2[0]->Fill(en,nc) ; }
    else if( nMax >= 3  ) { fhM02NLocMaxN[0]->Fill(en,l0) ; fhNCellNLocMaxN[0]->Fill(en,nc) ; }
    else printf("N max smaller than 1 -> %d \n",nMax);

    // Play with the MC stack if available
    // Check origin of the candidates
    Int_t mcindex = -1;
    if(IsDataMC()){
      
      Int_t tag	= GetMCAnalysisUtils()->CheckOrigin(cluster->GetLabels(),cluster->GetNLabels(), GetReader(), 0);
            
      if      ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0)  )      mcindex = mcPi0;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEta)  )      mcindex = mcEta;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && 
               !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = mcPhoton;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton) && 
                GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) mcindex = mcConversion;
      else if ( GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))   mcindex = mcElectron;
      else                                                                                mcindex = mcHadron;
      
/*      printf("AliAnaInsideClusterInvariantMass::FillAnalysisMakeHistograms() - tag %d, photon %d, prompt %d, frag %d, isr %d, pi0 decay %d, eta decay %d, other decay %d \n conv %d, pi0 %d, hadron %d, electron %d, unk %d, muon %d,pion %d, proton %d, neutron %d, kaon %d, antiproton %d, antineutron %d, bad %d\n",tag,
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPhoton),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPrompt),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCFragmentation),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCISR),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay),

             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOther),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCUnknown),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCMuon), 
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPion),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCProton), 
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCKaon), 
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiProton), 
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCAntiNeutron),
             GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCBadLabel)
);
*/      
      
      fhNLocMax[mcindex]->Fill(en,nMax);
      
      if     (nMax == 1 ) { fhM02NLocMax1[mcindex]->Fill(en,l0) ; fhNCellNLocMax1[mcindex]->Fill(en,nc) ; }
      else if(nMax == 2 ) { fhM02NLocMax2[mcindex]->Fill(en,l0) ; fhNCellNLocMax2[mcindex]->Fill(en,nc) ; }
      else if(nMax >= 3 ) { fhM02NLocMaxN[mcindex]->Fill(en,l0) ; fhNCellNLocMaxN[mcindex]->Fill(en,nc) ; }
      
    }  
    
    if( l0 < fM02Cut) continue;    
        
    // Get the 2 max indeces and do inv mass
    
    absId1 = absIdList[0];
    TLorentzVector cellMomi = GetCellMomentum(absId1, cells);

    if     ( nMax == 2 ) absId2 = absIdList[1];
    else if( nMax == 1 ){
      //Find second highest energy cell
      Int_t enmax = 0 ;
      for(Int_t iDigit = 0 ; iDigit < cluster->GetNCells() ; iDigit++){
        Int_t absId = cluster->GetCellsAbsId()[iDigit];
        if( absId == absId1 ) continue ; 
        Float_t endig = cells->GetCellAmplitude(absId);
        RecalibrateCellAmplitude(endig,absId); 
        if(endig > enmax) {
          enmax  = endig ;
          absId2 = absId ;
        }
        
        //Get mass of all cell in cluster combinations
        
        
        Float_t enj = cells->GetCellAmplitude(absId);
        RecalibrateCellAmplitude(enj,absId); 
        
        if(enj<0.3) continue;
        
        TLorentzVector cellMomj = GetCellMomentum(absId, cells);
        
        fhInvMassAllCells[0]->Fill(en,(cellMomj+cellMomi).M());
        
        if(IsDataMC()) fhInvMassAllCells[mcindex]->Fill(en,(cellMomj+cellMomi).M());
        
      }// cell loop
      
    }
    else { // loop on maxima, find 2 highest
      
      Int_t enmax = 0 ;
      for(Int_t iDigit = 0 ; iDigit < nMax ; iDigit++){
        Float_t endig = maxEList[iDigit];
        if(endig > enmax) {
          endig  = enmax ;
          absId2 = absIdList[iDigit];
        }
      }// maxima loop
      
      // First max is not highest, check if there is other higher
      if(maxEList[0] < enmax){
      
        for(Int_t iDigit = 0 ; iDigit < nMax ; iDigit++){
          if(absId2 == absIdList[iDigit]) continue;
          Float_t endig = maxEList[iDigit];
          if(endig > enmax) {
            endig  = enmax ;
            absId1 = absIdList[iDigit];
          }
        }// maxima loop
        
      }
      
    }
    
    fhNLocMaxM02Cut[0]->Fill(en,nMax);
    
    TLorentzVector cellMom1 = GetCellMomentum(absId1, cells);
    TLorentzVector cellMom2 = GetCellMomentum(absId2, cells);
    
    Float_t mass = (cellMom1+cellMom2).M();
    
    if     (nMax==1) fhMassNLocMax1[0]->Fill(en,mass);    
    else if(nMax==2) fhMassNLocMax2[0]->Fill(en,mass);
    else if(nMax >2) fhMassNLocMaxN[0]->Fill(en,mass);    

    if     (mass < 0.1)               fhM02Con[0]->Fill(en,l0);
    else if(mass < 0.2)               fhM02Pi0[0]->Fill(en,l0);
    else if(mass < 0.6 && mass > 0.4) fhM02Eta[0]->Fill(en,l0);

    if(IsDataMC()){
            
      fhNLocMaxM02Cut[mcindex]->Fill(en,nMax);

      if     (nMax==1) fhMassNLocMax1[mcindex]->Fill(en,mass);    
      else if(nMax==2) fhMassNLocMax2[mcindex]->Fill(en,mass);
      else if(nMax >2) fhMassNLocMaxN[mcindex]->Fill(en,mass);
      
      if     (mass < 0.1)               fhM02Con[mcindex]->Fill(en,l0);
      else if(mass < 0.2)               fhM02Pi0[mcindex]->Fill(en,l0);
      else if(mass < 0.6 && mass > 0.4) fhM02Eta[mcindex]->Fill(en,l0);
      
      
    }//Work with MC truth first
  
    delete absIdList ;
    delete maxEList  ;

  }//loop
  
  if(GetDebug() > 1) printf("AliAnaInsideClusterInvariantMass::MakeAnalysisFillHistograms() - END \n");  
  
}

//______________________________________________________________________
void AliAnaInsideClusterInvariantMass::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print("");
  printf("Calorimeter     =     %s\n", fCalorimeter.Data()) ;
  printf("lambda 0 sqared >  %2.1f\n", fM02Cut);
  printf("    \n") ;
  
} 

//____________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::RecalibrateCellAmplitude(Float_t & amp, const Int_t id)
{
  //Recaculate cell energy if recalibration factor
  
  Int_t icol     = -1; Int_t irow     = -1; Int_t iRCU     = -1;
  Int_t nModule  = GetModuleNumberCellIndexes(id,fCalorimeter, icol, irow, iRCU);
  
  if (GetCaloUtils()->IsRecalibrationOn()) {
    if(fCalorimeter == "PHOS") {
      amp *= GetCaloUtils()->GetPHOSChannelRecalibrationFactor(nModule,icol,irow);
    }
    else		                   {
      amp *= GetCaloUtils()->GetEMCALChannelRecalibrationFactor(nModule,icol,irow);
    }
  }
}

//____________________________________________________________________________________________
void AliAnaInsideClusterInvariantMass::SplitEnergy(const Int_t absId1, const Int_t absId2, AliVCaloCells* cells,
                                                   Float_t & e1, Float_t & e2 )
{
  
  // Split energy of cluster between the 2 local maxima.
  
  Int_t icol1 = -1, irow1 = -1, iRCU1 = -1;
  Int_t sm1 = GetCaloUtils()->GetModuleNumberCellIndexes(absId1, fCalorimeter, icol1, irow1, iRCU1) ;
  Int_t icol2 = -1, irow2 = -1, iRCU2 = -1;
  Int_t sm2 = GetCaloUtils()->GetModuleNumberCellIndexes(absId2, fCalorimeter, icol2, irow2, iRCU2) ;
  
  if(sm1!=sm2) {
    if(sm1%2) icol1+=AliEMCALGeoParams::fgkEMCALCols;
    else      icol2+=AliEMCALGeoParams::fgkEMCALCols;
  }
  
/// continue here  
  
}


