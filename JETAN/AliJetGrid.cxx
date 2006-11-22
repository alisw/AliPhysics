/**************************************************************************
 * Copyright(c) 2001-2002, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliJetGrid.cxx,v 1.0 07/09/2006 */

//=========================================================================
//  *** 07 September 2006
//  This class allows to fill a really general grid in (eta,phi) 
//  Two types of grid can be setted : 
//     if fGrid==0 -> the complete acceptance of a rectangular detector acceptance is filled
//     if fGrid==1 -> the previous grid - an other rectangular detector acceptance is filled 
//  A (eta,phi(rad)) position can be extracted from an index value without looping over all entries
//  An index position can be extracted from a (eta,phi(rad)) position without looping over all entries
//
//  How to run it :
//  > AliJetGrid *grid = new AliJetGrid((NbinPhi-1),(NbinEta-1),PhiMin,PhiMax,EtaMin,EtaMax);
//  > grid->SetGridType(...); // 0 or 1 for the moment
//  > grid->SetMatrixIndexes();
//  > grid->SetIndexIJ();
//
//  Author : magali.estienne@ires.in2p3.fr
//=========================================================================

// Standard headers 
#include <Riostream.h>
// Root headers
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TArrayI.h>
// AliRoot headers
#include "AliJetGrid.h"

ClassImp(AliJetGrid)

//__________________________________________________________
AliJetGrid::AliJetGrid() {

  // Default constructor

  fNphi = 0;
  fNeta = 0;
  fPhi = 0;
  fEta = 0;
  fIndex = 0x0;
  fPhiMin = 0.; // acceptance removed inside
  fPhiMax = 0.; // acceptance removed inside
  fEtaMin = 0.; // acceptance removed inside
  fEtaMax = 0.; // acceptance removed inside
  fMinPhi = 0.; // total acceptance
  fMaxPhi = 0.; // total acceptance
  fMinEta = 0.; // total acceptance
  fMaxEta = 0.; // total acceptance
  fEtaBinInTPCAcc = 0;
  fPhiBinInTPCAcc = 0;
  fEtaBinInEMCalAcc = 0;
  fPhiBinInEMCalAcc = 0;
  fNbinPhi = 0;

  fDebug = 1;
}

//__________________________________________________________
AliJetGrid::AliJetGrid(Int_t nphi,Int_t neta,Double_t phiMin,Double_t phiMax,Double_t etaMin,Double_t etaMax) {

  // Standard constructor

  fNphi = nphi;
  fNeta = neta;
  fPhiMin = 0.; // rad - acceptance removed inside
  fPhiMax = 0.; // rad - acceptance removed inside
  fEtaMin = 0.; //  acceptance removed inside
  fEtaMax = 0.; //  acceptance removed inside
  fNbinPhi = 0;

  fMaxPhi = phiMax; // total acceptance
  fMinPhi = phiMin; // total acceptance
  fMaxEta = etaMax; // total acceptance
  fMinEta = etaMin; // total acceptance

  fPhi  = new TArrayD(fNphi+1);
  fEta  = new TArrayD(fNeta+1);
  fIndexI = new TArrayI((fNeta+1)*(fNphi+1)+1);
  fIndexJ = new TArrayI((fNeta+1)*(fNphi+1)+1);

  for(Int_t i=0; i<fNphi+1; i++) (*fPhi)[i] = (phiMax-phiMin)/fNphi*i+phiMin;
  for(Int_t i=0; i<fNeta+1; i++) (*fEta)[i] = (etaMax-etaMin)/fNeta*i+etaMin;

  if(fDebug > 3){
    for(Int_t i=0; i<(fNphi+1); i++)  cout << (*fPhi)[i] << endl;
    for(Int_t i=0; i<(fNeta+1); i++)  cout << (*fEta)[i] << endl;
  }

  fIndex = new TMatrixD(fNphi+1,fNeta+1);

}

//__________________________________________________________
AliJetGrid::AliJetGrid(const AliJetGrid& grid):TNamed(grid) {

  // Copy constructor

  fNphi = grid.fNphi;
  fNeta = grid.fNeta;
  fPhiMin = grid.fPhiMin;
  fPhiMax = grid.fPhiMax;
  fEtaMin = grid.fEtaMin;
  fEtaMax = grid.fEtaMax;
  fEtaBinInTPCAcc = grid.fEtaBinInTPCAcc;
  fPhiBinInTPCAcc = grid.fPhiBinInTPCAcc;
  fEtaBinInEMCalAcc = grid.fEtaBinInEMCalAcc;
  fPhiBinInEMCalAcc = grid.fPhiBinInEMCalAcc;
  fNbinPhi = grid.fNbinPhi;
  fMaxPhi = grid.fMaxPhi;
  fMinPhi = grid.fMinPhi;
  fMaxEta = grid.fMaxEta;
  fMinEta = grid.fMinEta;

  fPhi = new TArrayD(fNphi+1);
  for(Int_t i=0; i<fNphi+1; i++) (*fPhi)[i] = grid.fPhi->At(i);
  fEta = new TArrayD(fNeta+1);
  for(Int_t i=0; i<fNeta+1; i++) (*fEta)[i] = grid.fEta->At(i);

  fIndex = new TMatrixD(fNphi+1,fNeta+1);
  for(Int_t i=0; i<fNphi+1; i++) {
    for(Int_t j=0; j<fNeta+1; j++) (*fIndex)(i,j)=(*grid.fIndex)(i,j);
  }
}

//__________________________________________________________
AliJetGrid::~AliJetGrid() {

  // Destructor
  delete fPhi;
  delete fEta;
  delete fIndexI;
  delete fIndexJ;
  delete fIndex;
}

//__________________________________________________________
void AliJetGrid::InitParams(Double_t phiMinCal,Double_t phiMaxCal,Double_t etaMinCal,Double_t etaMaxCal) 
{ // To set initial parameters

  fPhiMin = phiMinCal; // rad
  fPhiMax = phiMaxCal; // rad
  fEtaMin = etaMinCal;
  fEtaMax = etaMaxCal;
  fNbinPhi = static_cast<int>(fPhiMin/(2*TMath::Pi()/fNphi));

  // Define some binning numbers
  if(fGrid==0){
    for(Int_t i=0; i<fNphi+1; i++) fPhiBinInTPCAcc++;
    fPhiBinInEMCalAcc = 0;
    
    for(Int_t i=0; i<fNeta+1; i++) fEtaBinInTPCAcc++;
    fEtaBinInEMCalAcc = 0;
  }

  if(fGrid==1){
    for(Int_t i=0; i<fNphi+1; i++) {
      fPhiBinInTPCAcc++;
      if(fPhi->At(i) >= fPhiMin &&
	 fPhi->At(i) <= fPhiMax)
	fPhiBinInEMCalAcc++;
    }
    for(Int_t i=0; i<fNeta+1; i++) {
      fEtaBinInTPCAcc++;
      if((fEta->At(i) >= fEtaMin) &&
	 (fEta->At(i) <= fEtaMax))
	fEtaBinInEMCalAcc++;
    }
  }

}

//__________________________________________________________
TArrayD* AliJetGrid::GetArrayEta() 
{ // Returns an array with the eta points

  return fEta;

}

//__________________________________________________________
TArrayD* AliJetGrid::GetArrayPhi() 
{ // Returns an array with the phi points

  return fPhi;

}

//__________________________________________________________
TMatrixD* AliJetGrid::GetIndexObject()
{ // Returns a pointer to the matrix

  return fIndex;

}

//__________________________________________________________
void AliJetGrid::GetAccParam(Int_t &nphi, Int_t &neta, Float_t &minphi, Float_t &maxphi, 
				Float_t &mineta, Float_t &maxeta)
{ // Returns EMCAL acceptance initially setted

  nphi = fNphi;
  neta = fNeta;
  minphi = fPhiMin;
  maxphi = fPhiMax;
  mineta = fEtaMin;
  maxeta = fEtaMax;
}

//__________________________________________________________
void AliJetGrid::GetBinParam(Int_t &phibintpc, Int_t &etabintpc, 
				Int_t &phibinemc, Int_t &etabinemc, Int_t &nbinphi)
{ // Returns number of bins in TPC and EMCAL geometry

  etabintpc = fEtaBinInTPCAcc;
  phibintpc = fPhiBinInTPCAcc;
  etabinemc = fEtaBinInEMCalAcc;
  phibinemc = fPhiBinInEMCalAcc;
  nbinphi = fNbinPhi;
}

//__________________________________________________________
Int_t AliJetGrid::GetIndexFromEtaPhi(Double_t phi,Double_t eta) const 
{ // Tells the index value of a corresponding (eta,phi) real position
  // Loop over all entries -> takes time. 
  // Used one time at the begining to fill the grids

  /*   this is how bins are numbered
   
       in all TPC        or      in TPC - EMCAL
                                ... ... ... ..
                                ---------------
          ...  ... .            10 |   |   | 11
          ---+---+---           ---------------
    ^      6 | 7 | 8     or      8 |   |   | 9 
    |     ---+---+---           --------------- 
           3 | 4 | 5             4 | 5 | 6 | 7
   phi    ---+---+---           ---------------
           0 | 1 | 2             0 | 1 | 2 | 3  

       
             eta ->
  */

  Int_t etaBin=0,phiBin=0,absID=0;
  Int_t etaBin2=0,etaBin3=0;

  // Fill all the grid in eta/phi (all TPC acceptance)
  //-----------------------------------------------------
  if(fGrid==0){ 
    if(eta <= fEta->At(0)) {
      etaBin = 0;
    } else if(eta >= fEta->At(fNeta)) {
      etaBin = fNeta;
    } else {
      for(Int_t i=0; i<fNeta+1; i++) {
	if(eta < fEta->At(i)) {
	  etaBin = i-1;
	  break;
	} 
      }
    }
    if(phi <= fPhi->At(0)) {
      phiBin = 0;
    } else if(phi >= fPhi->At(fNphi)) {
      phiBin = fNphi;
    } else {
      for(Int_t i=0; i<fNphi+1; i++) {
	if(phi < fPhi->At(i)) {
	  phiBin = i-1;
	  break;
	} 
      }
    }
    
    // Calculate absolute id
    absID = phiBin*(fNeta+1) + etaBin;

  }

  // Fill the grid but do not count id in EMCal acceptance
  //------------------------------------------------------
  if((eta >= fEtaMin && eta <= fEtaMax) &&
     (phi >= fPhiMin && phi <= fPhiMax)){
    etaBin = etaBin2 = etaBin3 = -1;
    phiBin = -1;
  }  

  if(fGrid == 1){ 
    if(phi<fPhiMin) {
      if(eta <= fEta->At(0)) {
	etaBin = 0;
      } else if(eta >= fEta->At(fNeta)) {
	etaBin = fNeta;
      } else {
	for(Int_t i=0; i<fNeta+1; i++) {
	  if(eta < fEta->At(i)) {
	    etaBin = i-1;
	    break;
	  } 
	}
      }
    }
    else {
      if((phi>=fPhiMin) && (phi<=fPhiMax)) {
	if(eta <= fEta->At(0)) {
	  etaBin2 = 0;
	} else if(eta >= fEta->At(fNeta)) {
	  etaBin2 = fNeta-fEtaBinInEMCalAcc;
	} else {
	  for(Int_t i=0; i<fNeta+1; i++) {
	    if(eta < fEta->At(i) && ((fEta->At(0)<eta && eta<fEtaMin) || 
				     (fEtaMax<eta && eta<fEta->At(fNeta)))) {
	      //		 cout << "i : " << i << endl;
	      if(eta<fEtaMin)
		etaBin2 = i-1;
	      else etaBin2 = i-1-fEtaBinInEMCalAcc;
	      break;
	    } 
	  }
	}
      }
      else {
	if(eta <= fEta->At(0)) {
	  etaBin3 = 0;
	} else if(eta >= fEta->At(fNeta)) {
	  etaBin3 = fNeta;
	} else {
	  for(Int_t i=0; i<fNeta+1; i++) {
	    if(eta < fEta->At(i)) {
	      etaBin3 = i-1;
	      break;
	    } 
	  }
	}
      }
    }
    
    if(phi <= fPhi->At(0)) {
      phiBin = 0;
    } else if(phi >= fPhi->At(fNphi)) {
      phiBin = fNphi;
    } else {
      for(Int_t i=0; i<fNphi+1; i++) {
	if(phi < fPhi->At(i)) {
	  phiBin = i-1;
	  break;
	} 
      }
    }

    // Calculate absID
    //-------------------------------------------------------
    if(phi<fPhiMin)	 
      absID = phiBin*(fNeta+1) + etaBin;
    if(phi>=fPhiMin && phi<=fPhiMax) {
      if(eta >= fEtaMin && eta <= fEtaMax) absID = -1;
      else{
	absID = (fNbinPhi+1)*(fNeta+1)
	  + (phiBin-fNbinPhi-1)*((fEtaBinInTPCAcc-fEtaBinInEMCalAcc))
	  + etaBin2;
      }
    }
    if(phi>fPhiMax)
      absID = (fNbinPhi+1)*(fNeta+1)
	+ fPhiBinInEMCalAcc*((fEtaBinInTPCAcc-fEtaBinInEMCalAcc)) 
	+ (phiBin-(fNbinPhi+1+fPhiBinInEMCalAcc))*(fNeta+1)
	+ etaBin3;
    
  } // END OPTION==1    
  
  return absID;
  
}

//__________________________________________________________
void AliJetGrid::GetEtaPhiFromIndex(Int_t index, Float_t &eta, Float_t &phi)
{ // Get (eta,phi) position for a given index BUT loop over all entries (takes time)

  for(Int_t j=0; j<fNphi+1; j++) {
    for(Int_t i=0; i<fNeta+1; i++) {

      // TPC grid only 
      //-------------------------------------
      if(fGrid==0) {	
	if(j*(fNeta+1)+i == index) {
	  eta = fEta->At(i); 
	  phi = fPhi->At(j);
	}
      }

      // TPC-EMCAL grid
      //-------------------------------------
      Int_t ii = 0;
      if(i==0) ii = 0;
      if(i>0 && i<(fEtaBinInTPCAcc-fEtaBinInEMCalAcc)/2) ii = i; 
      if(i>=(fEtaBinInTPCAcc+fEtaBinInEMCalAcc)/2 && i<fNeta+1) ii = i-fEtaBinInEMCalAcc;

      if(fGrid==1) {
	if(j<(fNbinPhi+1) && j*(fNeta+1)+i == index) {
	  eta = fEta->At(i);
	  phi = fPhi->At(j);
	}  

	if((j>=(fNbinPhi+1) && j<(fNbinPhi+1+fPhiBinInEMCalAcc)) && 
	   ((fNbinPhi+1)*(fNeta+1) + (j-fNbinPhi-1)*(fEtaBinInTPCAcc-fEtaBinInEMCalAcc) + ii)== index ) {
	  if(ii==0) {Int_t ind = 0; eta = fEta->At(ind);}
	  else eta = fEta->At(i);
	  phi = fPhi->At(j);
	}

	if(j>=(fNbinPhi+1+fPhiBinInEMCalAcc) && 
	   ((fNbinPhi+1)*(fNeta+1)+fPhiBinInEMCalAcc*((fEtaBinInTPCAcc-fEtaBinInEMCalAcc))
	    +(j-(fNbinPhi+1+fPhiBinInEMCalAcc))*(fNeta+1)+i == index)) {
	  eta = fEta->At(i);
	  phi = fPhi->At(j);
	}
      }
    }
  }
}

//__________________________________________________________
Int_t AliJetGrid::GetIndex(Double_t phi, Double_t eta) 
{ // Get index value for a (eta,phi) position - Direct value

  Int_t ieta = GetIndexJFromEta(eta);
  Int_t iphi = GetIndexIFromPhi(phi);

  Int_t index = GetMatrixIndex(iphi,ieta);

  if(fDebug>10){
    cout << "(phi,eta) : " << phi << ", " << eta << endl;
    cout << "index : " << index << endl;
  }
  return index;

}

//__________________________________________________________
Int_t AliJetGrid::GetIndexJFromEta(Double_t eta)
{ // Get eta id 
  // Eta discretized

  Int_t idEta =0;
  Double_t temp = (eta+fMaxEta)/(fMaxEta-fMinEta)*fNeta;
  if(fDebug>20)
    {
      cout << "eta : " << eta << endl;
      cout << "fMaxEta : " << fMaxEta << endl;
      cout << "fMinEta : " << fMinEta << endl;
      cout << "fNeta : " << fNeta << endl;
      cout << "temp eta before cast : " << temp << endl;
    }
  idEta = static_cast<Int_t>(temp+0.5);
  if(fDebug>20) cout << "temp eta after cast : " << idEta << endl;
  return idEta;
}
//__________________________________________________________
Int_t AliJetGrid::GetIndexIFromPhi(Double_t phi)
{ // Get phi id
  // Phi discretized

  Int_t idPhi = 0;
  Double_t temp = 0.;
  if(fMinPhi==0) temp = phi/(fMaxPhi-fMinPhi)*fNphi;
  else temp = (phi-fMinPhi)/(fMaxPhi-fMinPhi)*fNphi;

  if(fDebug>20)
    {
      cout << "phi : " << phi << endl;
      cout << "fMaxPhi : " << fMaxPhi << endl;
      cout << "fMinPhi : " << fMinPhi << endl;
      cout << "fNphi : " << fNphi << endl;
      cout << "temp phi before cast : " << temp << endl;
    }
  idPhi = static_cast<Int_t>(temp+0.5);
  if(fDebug>20) cout << "temp phi after cast : " << idPhi << endl;
  return idPhi;


}

//__________________________________________________________
void AliJetGrid::SetMatrixIndex(Int_t i,Double_t par) 
{ // Allows to set parameters using only one index (if fGrid==0) !!
  // Not used !

  Int_t iphi = (Int_t)i/fNeta;
  Int_t ieta = i-iphi*fNeta;
  SetMatrixIndex(iphi,ieta,par);

  return;
}

//__________________________________________________________
void AliJetGrid::SetMatrixIndexes() 
{ // Fill the final matrix object with the corresponding index in eta/phi

  for(Int_t i=0; i<fNphi+1; i++){
    for(Int_t j=0; j<fNeta+1; j++){
      (*fIndex)(i,j) = GetIndexFromEtaPhi(fPhi->At(i),fEta->At(j))+1;
      if(fDebug>2){
	cout << "(*fIndex)(" << i << "," << j << ") : " << (*fIndex)(i,j) << 
	  ", phi : " << fPhi->At(i) << ", eta : " << fEta->At(j) << endl;
      }
    }
  }
  printf("##########################################################\n");
  printf("TMatrix object filled !\n");  
  printf("Size of the object : phi x eta = (fNphi+1) x (fNeta+1) = %d\n",(fNphi+1)*(fNeta+1));
  printf("##########################################################\n");
}

//__________________________________________________________
void AliJetGrid::SetIndexIJ()
{ // 

  for(Int_t i=0; i<fNphi+1; i++){
    for(Int_t j=0; j<fNeta+1; j++){
      Int_t id = static_cast<Int_t>((*fIndex)(i,j));

      if(id!=-1)
	{
	  (*fIndexI)[id] = i;
	  (*fIndexJ)[id] = j;
	}
    }
  }

  printf("##########################################################\n");
  printf("     In SetIndexIJ - Grid indexes setted !\n");
  printf("##########################################################\n");
}

//__________________________________________________________
void AliJetGrid::GetIJFromIndex(Int_t index, Int_t i, Int_t j)
{ // returns i position id of eta and j position id of phi for a given grid index
  i = (*fIndexI)[index];
  j = (*fIndexJ)[index];
}

//__________________________________________________________
void AliJetGrid::GetEtaPhiFromIndex2(Int_t index, Float_t &phi, Float_t &eta)
{ // returns eta, phi values for a given grid index

  phi = fPhi->At((*fIndexI)[index]);
  eta = fEta->At((*fIndexJ)[index]);
}

//__________________________________________________________
Int_t AliJetGrid::GetNEntries()
{ // Returns the number of entries of the grid

  if(fDebug>20){
    cout << "fMaxPhi : " << fMaxPhi << endl;
    cout << "fMaxEta : " << fMaxEta << endl;
  }

  Int_t indexNum = GetIndex(fMaxPhi,fMaxEta);
  if(fDebug>20) cout << "indexNum : " << indexNum << endl;
  return indexNum;

}

//__________________________________________________________
Int_t AliJetGrid::GetNEntries2()
{ // Returns the number of entries of the grid

  Int_t indexNum = GetIndex(fMaxPhi-1.,fMaxEta-0.5);
    if(fDebug>20) cout << "indexNum : " << indexNum << endl;
  return indexNum;

}






