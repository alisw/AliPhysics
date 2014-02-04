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
// Class containing methods for the isolation cut. 
// An AOD candidate (AliAODPWG4ParticleCorrelation type)
// is passed. Look in a cone around the candidate and study
// the hadronic activity inside to decide if the candidate is isolated
//
//
//*-- Author: Gustavo Conesa (LNF-INFN) 

//-Yaxian Mao (add the possibility for different IC method with different pt range, 01/10/2010)
//-Yaxian Mao (check the candidate particle is the leading particle or not at the same hemishere)

//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TLorentzVector.h>
#include <TObjArray.h>

// --- AliRoot system --- 
#include "AliIsolationCut.h" 
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliCalorimeterUtils.h"
#include "AliAODTrack.h"
#include "AliVCluster.h"
#include "AliCaloTrackReader.h"
#include "AliMixedEvent.h"
#include "AliCaloPID.h"

ClassImp(AliIsolationCut)
  
//____________________________________
AliIsolationCut::AliIsolationCut() : 
TObject(),
fConeSize(0.),
fPtThreshold(0.),
fPtThresholdMax(10000.), 
fSumPtThreshold(0.), 
fPtFraction(0.), 
fICMethod(0),
fPartInCone(0),
fDebug(-1),
fFracIsThresh(1)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();
  
}

//_________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateUEBandClusterNormalization(AliCaloTrackReader * /*reader*/, Float_t   etaC, Float_t /*phiC*/,
                                                           Float_t   phiUEptsumCluster,     Float_t   etaUEptsumCluster,
                                                           Float_t & phiUEptsumClusterNorm, Float_t & etaUEptsumClusterNorm,
                                                           Float_t & excessFracEta,         Float_t & excessFracPhi         ) const
{
  // Normalize cluster background band
  
  Float_t coneA     = fConeSize*fConeSize*TMath::Pi(); // A = pi R^2, isolation cone area
  
  //Careful here if EMCal limits changed .. 2010 (4 SM) to 2011-12 (10 SM), for the moment consider 100 deg in phi
  Float_t emcEtaSize = 0.7*2; // TO FIX
  Float_t emcPhiSize = TMath::DegToRad()*100.; // TO FIX
  
  /* //Catherine code
   if(((((2*fConeSize*emcPhiSize)-coneA))*phiBandBadCellsCoeff)!=0)phiUEptsumClusterNorm = phiUEptsumCluster*(coneA*coneBadCellsCoeff / (((2*fConeSize*emcPhiSize)-coneA))*phiBandBadCellsCoeff); // pi * R^2 / (2 R * 2 100 deg) -  trigger cone
   if(((((2*(fConeSize-excess)*emcPhiSize)-(coneA-excessFracEta))*etaBandBadCellsCoeff))!=0)phiUEptsumClusterNorm = phiUEptsumCluster*(coneA *coneBadCellsCoeff/ (((2*(fConeSize-excess)*emcPhiSize)-(coneA/excessFracEta))*etaBandBadCellsCoeff));
   if(((2*(fConeSize-excess)*emcEtaSize)-(coneA-excessFracPhi))*phiBandBadCellsCoeff!=0) etaUEptsumClusterNorm = etaUEptsumCluster*(coneA*coneBadCellsCoeff / (((2*(fConeSize-excess)*emcEtaSize)-(coneA/excessFracPhi))*phiBandBadCellsCoeff));
   */
  
  if((2*fConeSize*emcPhiSize-coneA)!=0) phiUEptsumClusterNorm = phiUEptsumCluster*(coneA / (((2*fConeSize*emcPhiSize)-coneA))); // pi * R^2 / (2 R * 2 100 deg) -  trigger cone
  if((2*fConeSize*emcEtaSize-coneA)!=0) etaUEptsumClusterNorm = etaUEptsumCluster*(coneA / (((2*fConeSize*emcEtaSize)-coneA))); // pi * R^2 / (2 R * 2*0.7)  -  trigger cone
  
  //out of eta acceptance
  excessFracEta = 1;
  excessFracPhi = 1;

  if(TMath::Abs(etaC)+fConeSize > emcEtaSize/2.)
  {
    Float_t excess = TMath::Abs(etaC) + fConeSize - emcEtaSize/2.;
    excessFracEta  = CalculateExcessAreaFraction(excess);
    
    if ( excessFracEta != 0) coneA /=  excessFracEta;
    
    //UE band is also out of acceptance, need to estimate corrected area
    if(((2*fConeSize-excess)*emcPhiSize-coneA) != 0 ) phiUEptsumClusterNorm = phiUEptsumCluster*(coneA / ((((2*fConeSize-excess)*emcPhiSize)-coneA)));
    if(( 2*fConeSize        *emcEtaSize-coneA) != 0 ) etaUEptsumClusterNorm = etaUEptsumCluster*(coneA / ((( 2*fConeSize        *emcEtaSize)-coneA)));
  }
  
}

//________________________________________________________________________________________________________________________________
void AliIsolationCut::CalculateUEBandTrackNormalization  (AliCaloTrackReader * reader,    Float_t   etaC, Float_t /*phiC*/,
                                                          Float_t   phiUEptsumTrack,      Float_t   etaUEptsumTrack,
                                                          Float_t & phiUEptsumTrackNorm,  Float_t & etaUEptsumTrackNorm,
                                                          Float_t & excessFracEta,        Float_t & excessFracPhi          ) const
{
  // Normalize track background band
  
  Float_t coneA     = fConeSize*fConeSize*TMath::Pi(); // A = pi R^2, isolation cone area
  
  // Get the cut used for the TPC tracks in the reader, +-0.8, +-0.9 ...
  // Only valid in simple fidutial cut case and if the cut is applied, careful!
  Float_t tpcEtaSize = reader->GetFiducialCut()->GetCTSFidCutMaxEtaArray()->At(0) -
  reader->GetFiducialCut()->GetCTSFidCutMinEtaArray()->At(0) ;
  Float_t tpcPhiSize = TMath::TwoPi();
  
  /*//Catherine code
   //phiUEptsumTrackNorm = phiUEptsumTrack*(coneA*coneBadCellsCoeff / (((2*fConeSize*tpcPhiSize)-coneA))*phiBandBadCellsCoeff); // pi * R^2 / (2 R * 2 pi) -  trigger cone
   //etaUEptsumTrackNorm = etaUEptsumTrack*(coneA*coneBadCellsCoeff / (((2*fConeSize*tpcEtaSize)-coneA))*etaBandBadCellsCoeff); // pi * R^2 / (2 R * 1.6)  -  trigger cone
   if((2*fConeSize*tpcPhiSize-coneA)!=0)phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / (((2*fConeSize*tpcPhiSize)-coneA))); // pi * R^2 / (2 R * 2 pi) -  trigger cone
   if((2*fConeSize*tpcEtaSize-coneA)!=0)etaUEptsumTrackNorm = etaUEptsumTrack*(coneA / (((2*fConeSize*tpcEtaSize)-coneA))); // pi * R^2 / (2 R * 1.6)  -  trigger cone
   if((2*(fConeSize-excess)*tpcPhiSize)-(coneA-excessFracEta)!=0)phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / (((2*(fConeSize-excess)*tpcPhiSize)-(coneA/excessFracEta))));
   */ //end Catherine code
  
  //correct out of eta acceptance
  excessFracEta = 1;
  excessFracPhi = 1;

  if((2*fConeSize*tpcPhiSize-coneA)!=0) phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / (((2*fConeSize*tpcPhiSize)-coneA))); // pi * R^2 / (2 R * 2 pi) -  trigger cone
  if((2*fConeSize*tpcEtaSize-coneA)!=0) etaUEptsumTrackNorm = etaUEptsumTrack*(coneA / (((2*fConeSize*tpcEtaSize)-coneA))); // pi * R^2 / (2 R * 1.6)  -  trigger cone
  
  if(TMath::Abs(etaC)+fConeSize > tpcEtaSize/2.)
  {
    Float_t excess = TMath::Abs(etaC) + fConeSize - tpcEtaSize/2.;
    excessFracEta  = CalculateExcessAreaFraction(excess);
    if (excessFracEta != 0) coneA /=  excessFracEta;
    
    //UE band is also out of acceptance, need to estimate corrected area
    if(((2*fConeSize-excess)*tpcPhiSize - coneA) !=0 ) phiUEptsumTrackNorm = phiUEptsumTrack*(coneA / ((((2*fConeSize-excess)*tpcPhiSize)-coneA)));
    if(( 2*fConeSize        *tpcEtaSize - coneA) !=0 ) etaUEptsumTrackNorm = etaUEptsumTrack*(coneA / ((( 2*fConeSize        *tpcEtaSize)-coneA)));
  }
  
}

//________________________________________________________________________
Float_t AliIsolationCut::CalculateExcessAreaFraction(Float_t excess) const
{
  // Area of a circunference segment segment 1/2 R^2 (angle-sin(angle)), angle = 2*ACos((R-excess)/R)
  
  
  Float_t angle   = 2*TMath::ACos( (fConeSize-excess) / fConeSize );
  
  Float_t coneA   = fConeSize*fConeSize*TMath::Pi(); // A = pi R^2, isolation cone area
  
  Float_t excessA = fConeSize*fConeSize / 2 * (angle-TMath::Sin(angle));
  
  if(coneA > excessA) return coneA / (coneA-excessA);
  else
  {
    printf("AliIsolationCut::CalculateExcessAreaFraction() - Please Check : Excess Track %2.3f, coneA %2.2f,  excessA %2.2f, angle %2.2f,factor %2.2f\n",
           excess,coneA, excessA, angle*TMath::RadToDeg(), coneA / (coneA-excessA));
    return  1;
  }
}

//_______________________________________________________________________________________
Float_t AliIsolationCut::GetCellDensity(AliAODPWG4ParticleCorrelation * pCandidate,
                                        AliCaloTrackReader * reader) const
{
  // Get good cell density (number of active cells over all cells in cone)
  
  Double_t coneCells    = 0.; //number of cells in cone with radius fConeSize
  Double_t coneCellsBad = 0.; //number of bad cells in cone with radius fConeSize
  Double_t cellDensity  = 1.;

  Float_t phiC  = pCandidate->Phi() ;
  if(phiC<0) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;
  
  if(pCandidate->GetDetector()=="EMCAL")
  {
    AliEMCALGeometry* eGeom = AliEMCALGeometry::GetInstance();
    AliCalorimeterUtils *cu = reader->GetCaloUtils();
    
    Int_t absId = -999;
    if (eGeom->GetAbsCellIdFromEtaPhi(etaC,phiC,absId))
    {
      //Get absolute (col,row) of candidate
      Int_t iEta=-1, iPhi=-1, iRCU = -1;      
      Int_t nSupMod = cu->GetModuleNumberCellIndexes(absId, pCandidate->GetDetector(), iEta, iPhi, iRCU);
      
      Int_t colC = iEta;
      if (nSupMod % 2) colC =  AliEMCALGeoParams::fgkEMCALCols + iEta ;
      Int_t rowC = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(nSupMod/2);
      
      Int_t sqrSize = int(fConeSize/0.0143) ; // Size of cell in radians
      //loop on cells in a square of side fConeSize to check cells in cone    
      for(Int_t icol = colC-sqrSize; icol < colC+sqrSize;icol++)
      {
        for(Int_t irow = rowC-sqrSize; irow < rowC+sqrSize; irow++)
        {
          if (Radius(colC, rowC, icol, irow) < sqrSize)
          {
            coneCells += 1.;
            
            Int_t cellSM  = -999;
            Int_t cellEta = -999;
            Int_t cellPhi = -999;
            if(icol > AliEMCALGeoParams::fgkEMCALCols-1) 
            {
              cellSM = 0+int(irow/AliEMCALGeoParams::fgkEMCALRows)*2;
              cellEta = icol-AliEMCALGeoParams::fgkEMCALCols;
              cellPhi = irow-AliEMCALGeoParams::fgkEMCALRows*int(cellSM/2);
            }
            if(icol < AliEMCALGeoParams::fgkEMCALCols) 
            {
              cellSM = 1+int(irow/AliEMCALGeoParams::fgkEMCALRows)*2;
              cellEta = icol;
              cellPhi = irow-AliEMCALGeoParams::fgkEMCALRows*int(cellSM/2);
            }
            
            //Count as bad "cells" out of EMCAL acceptance
            if(icol < 0 || icol > AliEMCALGeoParams::fgkEMCALCols*2 || 
               irow < 0 || irow > AliEMCALGeoParams::fgkEMCALRows*16./3) //5*nRows+1/3*nRows
            {
              coneCellsBad += 1.;
            }
            //Count as bad "cells" marked as bad in the DataBase
            else if (cu->GetEMCALChannelStatus(cellSM,cellEta,cellPhi)==1) 
            {
              coneCellsBad += 1. ;
            }
          }
        }
      }//end of cells loop
    }
    
    else if(fDebug>0) printf("cluster with bad (eta,phi) in EMCal for energy density calculation\n");
    
    if (coneCells > 0.) 
    {
      cellDensity = (coneCells-coneCellsBad)/coneCells;
      //printf("Energy density = %f\n", cellDensity);
    }
  }

  return cellDensity;
  
}

//__________________________________________________________________________________
void AliIsolationCut::GetCoeffNormBadCell(AliAODPWG4ParticleCorrelation * pCandidate,
                                          AliCaloTrackReader * reader,
                                          Float_t &  coneBadCellsCoeff,
                                          Float_t &  etaBandBadCellsCoeff,
                                          Float_t & phiBandBadCellsCoeff)
{
  // Get good cell density (number of active cells over all cells in cone)
  
  Double_t coneCells    = 0.; //number of cells in cone with radius fConeSize
  Double_t phiBandCells = 0.; //number of cells in band phi
  Double_t etaBandCells = 0.; //number of cells in band eta
  
  Float_t phiC  = pCandidate->Phi() ;
  if(phiC<0) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;
  
  if(pCandidate->GetDetector()=="EMCAL")
  {
    AliEMCALGeometry* eGeom = AliEMCALGeometry::GetInstance();
    AliCalorimeterUtils *cu = reader->GetCaloUtils();
    
    Int_t absId = -999;
    if (eGeom->GetAbsCellIdFromEtaPhi(etaC,phiC,absId))
    {
      //Get absolute (col,row) of candidate
      Int_t iEta=-1, iPhi=-1, iRCU = -1;
      Int_t nSupMod = cu->GetModuleNumberCellIndexes(absId, pCandidate->GetDetector(),
                                                     iEta, iPhi, iRCU);
      
      Int_t colC = iEta;
      if (nSupMod % 2) colC =  AliEMCALGeoParams::fgkEMCALCols + iEta ;
      Int_t rowC = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(nSupMod/2);
      
      Int_t sqrSize = int(fConeSize/0.0143) ; // Size of cell in radians
      for(Int_t icol = 0; icol < 2*AliEMCALGeoParams::fgkEMCALCols-1;icol++)
      {
        for(Int_t irow = 0; irow < 5*AliEMCALGeoParams::fgkEMCALRows -1; irow++)
        {
          //loop on cells in a square of side fConeSize to check cells in cone
          if     ( Radius(colC, rowC, icol, irow) < sqrSize ) { coneCells    += 1.; }
          else if( icol>colC-sqrSize  &&  icol<colC+sqrSize ) { phiBandCells += 1 ; }
          else if( irow>rowC-sqrSize  &&  irow<rowC+sqrSize ) { etaBandCells += 1 ; }
          
          Int_t cellSM  = -999;
          Int_t cellEta = -999;
          Int_t cellPhi = -999;
          if(icol > AliEMCALGeoParams::fgkEMCALCols-1)
          {
            cellSM = 0+int(irow/AliEMCALGeoParams::fgkEMCALRows)*2;
            cellEta = icol-AliEMCALGeoParams::fgkEMCALCols;
            cellPhi = irow-AliEMCALGeoParams::fgkEMCALRows*int(cellSM/2);
          }
          if(icol < AliEMCALGeoParams::fgkEMCALCols)
          {
            cellSM = 1+int(irow/AliEMCALGeoParams::fgkEMCALRows)*2;
            cellEta = icol;
            cellPhi = irow-AliEMCALGeoParams::fgkEMCALRows*int(cellSM/2);
          }
          
          if( (icol < 0 || icol > AliEMCALGeoParams::fgkEMCALCols*2-1 ||
               irow < 0 || irow > AliEMCALGeoParams::fgkEMCALRows*5 - 1) //5*nRows+1/3*nRows //Count as bad "cells" out of EMCAL acceptance
             || (cu->GetEMCALChannelStatus(cellSM,cellEta,cellPhi)==1))  //Count as bad "cells" marked as bad in the DataBase
          {
            if     ( Radius(colC, rowC, icol, irow) < sqrSize ) coneBadCellsCoeff    += 1.;
            else if( icol>colC-sqrSize  &&  icol<colC+sqrSize ) phiBandBadCellsCoeff += 1 ;
	          else if( irow>rowC-sqrSize  &&  irow<rowC+sqrSize ) etaBandBadCellsCoeff += 1 ;
          }
        }
      }//end of cells loop
    }
    
    else if(fDebug > 0) printf("cluster with bad (eta,phi) in EMCal for energy density coeff calculation\n");
    
    if (coneCells > 0.)
    {
      //   printf("Energy density coneBadCellsCoeff= %.2f coneCells%.2f\n", coneBadCellsCoeff,coneCells);
      coneBadCellsCoeff = (coneCells-coneBadCellsCoeff)/coneCells;
      //  printf("coneBadCellsCoeff= %.2f\n", coneBadCellsCoeff);
    }
    if (phiBandCells > 0.)
    {
      // printf("Energy density phiBandBadCellsCoeff = %.2f phiBandCells%.2f\n", phiBandBadCellsCoeff,phiBandCells);
      phiBandBadCellsCoeff = (phiBandCells-phiBandBadCellsCoeff)/phiBandCells;
      // printf("phiBandBadCellsCoeff = %.2f\n", phiBandBadCellsCoeff);
    }
    if (etaBandCells > 0.)
    {
      //printf("Energy density etaBandBadCellsCoeff = %.2f etaBandCells%.2f\n", etaBandBadCellsCoeff,etaBandCells);
      etaBandBadCellsCoeff = (etaBandCells-etaBandBadCellsCoeff)/etaBandCells;
      // printf("etaBandBadCellsCoeff = %.2f\n",etaBandBadCellsCoeff);
    }
    
  }
  
}

//____________________________________________
TString AliIsolationCut::GetICParametersList()
{
  //Put data member values in string to keep in output container
  
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  
  snprintf(onePar,buffersize,"--- AliIsolationCut ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"fConeSize: (isolation cone size) %1.2f\n",fConeSize) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPtThreshold =%1.2f (isolation pt threshold) \n",fPtThreshold) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPtFraction=%1.2f (isolation pt threshold fraction ) \n",fPtFraction) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fICMethod=%d (isolation cut case) \n",fICMethod) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fPartInCone=%d \n",fPartInCone) ;
  parList+=onePar ;
 snprintf(onePar,buffersize,"fFracIsThresh=%i \n",fFracIsThresh) ;
  parList+=onePar ;
 
  return parList; 
}

//____________________________________
void AliIsolationCut::InitParameters()
{
  //Initialize the parameters of the analysis.
  
  fConeSize       = 0.4 ; 
  fPtThreshold    = 1.  ;
  fPtThresholdMax = 10000.  ;
  fSumPtThreshold = 0.5 ; 
  fPtFraction     = 0.1 ; 
  fPartInCone     = kOnlyCharged;
  fICMethod       = kSumPtFracIC; // 0 pt threshol method, 1 cone pt sum method
  fFracIsThresh   = 1; 
}

//________________________________________________________________________________
void  AliIsolationCut::MakeIsolationCut(TObjArray * plCTS, 
                                        TObjArray * plNe, 
                                        AliCaloTrackReader * reader, 
                                        AliCaloPID * pid,
                                        Bool_t bFillAOD,
                                        AliAODPWG4ParticleCorrelation  *pCandidate, 
                                        TString aodArrayRefName,
                                        Int_t   & n, 
                                        Int_t   & nfrac, 
                                        Float_t & coneptsum,  
                                        Bool_t  & isolated) const
{
  //Search in cone around a candidate particle if it is isolated 
  Float_t ptC   = pCandidate->Pt() ;
  Float_t phiC  = pCandidate->Phi() ;
  if(phiC<0) phiC+=TMath::TwoPi();
  Float_t etaC  = pCandidate->Eta() ;
  
  Float_t pt    = -100. ;
  Float_t eta   = -100. ;
  Float_t phi   = -100. ;
  Float_t rad   = -100. ;
  
  Float_t coneptsumCluster = 0;
  Float_t coneptsumTrack   = 0;
  
  Float_t  etaBandPtSumTrack   = 0;
  Float_t  phiBandPtSumTrack   = 0;
  Float_t  etaBandPtSumCluster = 0;
  Float_t  phiBandPtSumCluster = 0;
  
  n         = 0 ;
  nfrac     = 0 ;
  isolated  = kFALSE;
  
  if(fDebug>0) 
  {
    printf("AliIsolationCut::MakeIsolationCut() - Cadidate pT %2.2f, eta %2.2f, phi %2.2f, cone %1.2f, thres %2.2f, Fill AOD? %d",
           pCandidate->Pt(), pCandidate->Eta(), pCandidate->Phi()*TMath::RadToDeg(), fConeSize,fPtThreshold,bFillAOD);
    if(plCTS) printf(", nTracks %d"  ,plCTS->GetEntriesFast());
    if(plNe)  printf(", nClusters %d",plNe ->GetEntriesFast());
    
    printf("\n");
  }
  
  //Initialize the array with refrences
  TObjArray * refclusters  = 0x0;
  TObjArray * reftracks    = 0x0;
  Int_t       ntrackrefs   = 0;
  Int_t       nclusterrefs = 0;
  
  //Check charged particles in cone.
  if(plCTS && 
     (fPartInCone==kOnlyCharged || fPartInCone==kNeutralAndCharged))
  {
    TVector3 p3;
    for(Int_t ipr = 0;ipr < plCTS->GetEntries() ; ipr ++ )
    {
      AliVTrack* track = dynamic_cast<AliVTrack*>(plCTS->At(ipr)) ; 
      
      if(track)
      {
        //Do not count the candidate (pion, conversion photon) or the daughters of the candidate
        if(track->GetID() == pCandidate->GetTrackLabel(0) || track->GetID() == pCandidate->GetTrackLabel(1) || 
           track->GetID() == pCandidate->GetTrackLabel(2) || track->GetID() == pCandidate->GetTrackLabel(3)   ) continue ;
        
        p3.SetXYZ(track->Px(),track->Py(),track->Pz());
        pt  = p3.Pt();
        eta = p3.Eta();
        phi = p3.Phi() ;
      }
      else
      {// Mixed event stored in AliAODPWG4Particles
        AliAODPWG4Particle * trackmix = dynamic_cast<AliAODPWG4Particle*>(plCTS->At(ipr)) ; 
        if(!trackmix)
        {
          printf("AliIsolationCut::MakeIsolationCut() - Wrong track data type, continue\n");
          continue;
        }
        
        pt  = trackmix->Pt();
        eta = trackmix->Eta();
        phi = trackmix->Phi() ;
      }
      
      if( phi < 0 ) phi+=TMath::TwoPi();
      
      rad = Radius(etaC, phiC, eta, phi);
      
      // ** For the background out of cone **
      
      if(rad > fConeSize)
      {
        if(eta > (etaC-fConeSize) && eta < (etaC+fConeSize)) phiBandPtSumTrack += pt;
        if(phi > (phiC-fConeSize) && phi < (phiC+fConeSize)) etaBandPtSumTrack += pt;
      }
      
      // ** For the isolated particle **
      
      // Only loop the particle at the same side of candidate
      if(TMath::Abs(phi-phiC) > TMath::PiOver2()) continue ;
      
      // If at the same side has particle larger than candidate, 
      // then candidate can not be the leading, skip such events
      if(pt > ptC)
      {
        n         = -1;
        nfrac     = -1;
        coneptsumTrack = -1;
        isolated  = kFALSE;
        
        pCandidate->SetLeadingParticle(kFALSE);
        
        if(bFillAOD && reftracks) 
        {
          reftracks->Clear(); 
          delete reftracks;
        }
        
        return ;
      }
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      
      if( fDebug > 0 ) 
        printf("\t track %d, pT %2.2f, eta %1.2f, phi %2.2f, R candidate %2.2f", ipr,pt,eta,phi,rad);
      
      if(rad < fConeSize)
      {
        if(fDebug > 0)  printf(" -  inside candidate cone");
        
        if(bFillAOD)
        {
          ntrackrefs++;
          if(ntrackrefs == 1)
          {
            reftracks = new TObjArray(0);
            //reftracks->SetName(Form("Tracks%s",aodArrayRefName.Data()));
            TString tempo(aodArrayRefName)  ; 
            tempo += "Tracks" ; 
            reftracks->SetName(tempo);
            reftracks->SetOwner(kFALSE);
          }
          reftracks->Add(track);
        }
        
        coneptsumTrack+=pt;
        if(pt > fPtThreshold && pt < fPtThresholdMax)  n++;
        if(pt > fPtFraction*ptC ) nfrac++;  
        
      } // Inside cone
      
      if(fDebug>0)  printf("\n");
      
    }// charged particle loop
    
    
  }//Tracks
  
  
  //Check neutral particles in cone.  
  if(plNe &&
     (fPartInCone==kOnlyNeutral || fPartInCone==kNeutralAndCharged))
  {
    TLorentzVector mom ;
    
    for(Int_t ipr = 0;ipr < plNe->GetEntries() ; ipr ++ )
    {
      AliVCluster * calo = dynamic_cast<AliVCluster *>(plNe->At(ipr)) ;
      
      if(calo)
      {
        //Get the index where the cluster comes, to retrieve the corresponding vertex
        Int_t evtIndex = 0 ; 
        if (reader->GetMixedEvent()) 
          evtIndex=reader->GetMixedEvent()->EventIndexForCaloCluster(calo->GetID()) ; 
        
        
        //Do not count the candidate (photon or pi0) or the daughters of the candidate
        if(calo->GetID() == pCandidate->GetCaloLabel(0) || 
           calo->GetID() == pCandidate->GetCaloLabel(1)   ) continue ;      
        
        //Skip matched clusters with tracks in case of neutral+charged analysis
        if( fPartInCone == kNeutralAndCharged && 
           pid->IsTrackMatched(calo,reader->GetCaloUtils(),reader->GetInputEvent()) ) continue ;
        
        //Assume that come from vertex in straight line
        calo->GetMomentum(mom,reader->GetVertex(evtIndex)) ;
        
        pt  = mom.Pt()  ;
        eta = mom.Eta() ;
        phi = mom.Phi() ;
      }
      else 
      {// Mixed event stored in AliAODPWG4Particles
        AliAODPWG4Particle * calomix = dynamic_cast<AliAODPWG4Particle*>(plNe->At(ipr)) ; 
        if(!calomix)
        {
          printf("AliIsolationCut::MakeIsolationCut() - Wrong calo data type, continue\n");
          continue;
        }
        
        pt  = calomix->Pt();
        eta = calomix->Eta();
        phi = calomix->Phi() ;
      }
      
      if( phi < 0 ) phi+=TMath::TwoPi();
      
      rad = Radius(etaC, phiC, eta, phi);
      
      // ** For the background out of cone **
      
      if(rad > fConeSize)
      {
        if(eta > (etaC-fConeSize) && eta < (etaC+fConeSize)) phiBandPtSumCluster += pt;
        if(phi > (phiC-fConeSize) && phi < (phiC+fConeSize)) etaBandPtSumCluster += pt;
      }
      
      // ** For the isolated particle **
      
      // Only loop the particle at the same side of candidate
      if(TMath::Abs(phi-phiC)>TMath::PiOver2()) continue ;
      
      // If at the same side has particle larger than candidate, 
      // then candidate can not be the leading, skip such events
      if(pt > ptC)
      {
        n         = -1;
        nfrac     = -1;
        coneptsumCluster = -1;
        isolated  = kFALSE;
        
        pCandidate->SetLeadingParticle(kFALSE);
        
        if(bFillAOD)
        {
          if(reftracks)
          {  
            reftracks  ->Clear();
            delete reftracks;
          }
          
          if(refclusters)
          {
            refclusters->Clear(); 
            delete refclusters;
          }
        }
        return ;
      }
      
      //Check if there is any particle inside cone with pt larger than  fPtThreshold
      
      if(fDebug > 0 ) 
        printf("\t cluster %d, pT %2.2f, eta %1.2f, phi %2.2f, R candidate %2.2f", ipr,pt,eta,phi,rad);
      
      if(rad < fConeSize)
      {
        if(fDebug > 0 )  printf(" - inside candidate cone");
        
        if(bFillAOD) 
        {
          nclusterrefs++;
          if(nclusterrefs==1)
          {
            refclusters = new TObjArray(0);
            //refclusters->SetName(Form("Clusters%s",aodArrayRefName.Data()));
            TString tempo(aodArrayRefName)  ; 
            tempo += "Clusters" ; 
            refclusters->SetName(tempo);
            refclusters->SetOwner(kFALSE);
          }
          refclusters->Add(calo);
        }
        
        coneptsumCluster+=pt;
        if(pt > fPtThreshold && pt < fPtThresholdMax)  n++;
        //if fPtFraction*ptC<fPtThreshold then consider the fPtThreshold directly
        if(fFracIsThresh)
        {
          if( fPtFraction*ptC<fPtThreshold)
          {
            if(pt>fPtThreshold)    nfrac++ ;
          }
          else 
          {
            if(pt>fPtFraction*ptC) nfrac++; 
          }
        }
        else
        {
          if(pt>fPtFraction*ptC) nfrac++;   
        }
        
      }//in cone
      
      if(fDebug>0)  printf("\n");
      
    }// neutral particle loop
    
  }//neutrals
  
  
  //Add reference arrays to AOD when filling AODs only
  if(bFillAOD) 
  {
    if(refclusters)	pCandidate->AddObjArray(refclusters);
    if(reftracks)	  pCandidate->AddObjArray(reftracks);
  }
  
  coneptsum = coneptsumCluster+coneptsumTrack;
  
  //Check isolation, depending on selected isolation criteria
  if( fICMethod == kPtThresIC)
  {
    if(n==0) isolated = kTRUE ;
  }
  else if( fICMethod == kSumPtIC)
  {
    if(coneptsum < fSumPtThreshold)
      isolated  =  kTRUE ;
  }
  else if( fICMethod == kPtFracIC)
  {
    if(nfrac==0) isolated = kTRUE ;
  }
  else if( fICMethod == kSumPtFracIC)
  {
    //when the fPtFraction*ptC < fSumPtThreshold then consider the later case
 // printf("photon analysis IsDataMC() ?%i\n",IsDataMC());
    if(fFracIsThresh )
    {
      if( fPtFraction*ptC < fSumPtThreshold  && coneptsum < fSumPtThreshold) isolated  =  kTRUE ;
      if( fPtFraction*ptC > fSumPtThreshold  && coneptsum < fPtFraction*ptC) isolated  =  kTRUE ;
    }
    else 
    {
      if(coneptsum < fPtFraction*ptC) isolated  =  kTRUE ;
    }
  }
  else if( fICMethod == kSumDensityIC)
  {    
    // Get good cell density (number of active cells over all cells in cone)
    // and correct energy in cone
    
    Float_t cellDensity = GetCellDensity(pCandidate,reader);
    
    if(coneptsum < fSumPtThreshold*cellDensity)
      isolated = kTRUE;
  }
  else if( fICMethod == kSumBkgSubIC)
  {
    Double_t coneptsumBkg = 0.;
    Float_t  etaBandPtSumTrackNorm   = 0;
    Float_t  phiBandPtSumTrackNorm   = 0;
    Float_t  etaBandPtSumClusterNorm = 0;
    Float_t  phiBandPtSumClusterNorm = 0;
    
    Float_t  excessFracEtaTrack   = 1;
    Float_t  excessFracPhiTrack   = 1;
    Float_t  excessFracEtaCluster = 1;
    Float_t  excessFracPhiCluster = 1;
    
    // Normalize background to cone area
    if     (fPartInCone != kOnlyCharged       )
      CalculateUEBandClusterNormalization(reader, etaC, phiC,
                                          phiBandPtSumCluster    , etaBandPtSumCluster,
                                          phiBandPtSumClusterNorm, etaBandPtSumClusterNorm,
                                          excessFracEtaCluster   , excessFracPhiCluster    );
    if     (fPartInCone != kOnlyNeutral       )
      CalculateUEBandTrackNormalization(reader, etaC, phiC,
                                          phiBandPtSumTrack    , etaBandPtSumTrack  ,
                                          phiBandPtSumTrackNorm, etaBandPtSumTrackNorm,
                                          excessFracEtaTrack   , excessFracPhiTrack    );
    
    if     (fPartInCone == kOnlyCharged       ) coneptsumBkg = etaBandPtSumTrackNorm;
    else if(fPartInCone == kOnlyNeutral       ) coneptsumBkg = etaBandPtSumClusterNorm;
    else if(fPartInCone == kNeutralAndCharged ) coneptsumBkg = etaBandPtSumClusterNorm + etaBandPtSumTrackNorm;
    
    //coneptsumCluster*=(coneBadCellsCoeff*excessFracEtaCluster*excessFracPhiCluster) ; // apply this correction earlier???
    // line commented out in last modif!!!
    
    coneptsum = coneptsumCluster+coneptsumTrack;
    
    coneptsum -= coneptsumBkg;
    if(coneptsum < fSumPtThreshold)
      isolated  =  kTRUE ;
  } 
  
}

//_____________________________________________________
void AliIsolationCut::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print %s %s **** \n", GetName(), GetTitle() ) ;
  
  printf("IC method          =     %d\n",    fICMethod   ) ; 
  printf("Cone Size          =     %1.2f\n", fConeSize   ) ; 
  printf("pT threshold       =     %2.1f\n", fPtThreshold) ;
  printf("pT fraction        =     %3.1f\n", fPtFraction ) ;
  printf("particle type in cone =  %d\n",    fPartInCone ) ;
  printf("using fraction for high pt leading instead of frac ? %i\n",fFracIsThresh);
  printf("    \n") ;
  
} 

//___________________________________________________________________________
Float_t AliIsolationCut::Radius(Float_t etaC, Float_t phiC,
                                Float_t eta , Float_t phi) const
{
  // Calculate the distance to trigger from any particle

  Float_t dEta = etaC-eta;
  Float_t dPhi = phiC-phi;
  
  if(TMath::Abs(dPhi) >= TMath::Pi()) 
    dPhi = TMath::TwoPi()-TMath::Abs(dPhi);
  
  return TMath::Sqrt( dEta*dEta + dPhi*dPhi );
  
}



