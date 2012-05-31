/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
////////////////////////////////////////////////////////////////////////////
//            Implementation of the ITS clusterer V2 class                //
//                                                                        //
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch            //
//          Unfolding switch from AliITSRecoParam: D. Elia, INFN Bari     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include "AliITSCalibrationSPD.h"
#include "AliHLTITSClusterFinderSPD.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliITSDetTypeRec.h"
#include "AliITSReconstructor.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSPD.h"
#include <TClonesArray.h>
#include "AliITSdigitSPD.h"
#include "AliITSFOSignalsSPD.h"

ClassImp(AliHLTITSClusterFinderSPD)

//__________________________________________________________________________
AliHLTITSClusterFinderSPD::AliHLTITSClusterFinderSPD(AliITSDetTypeRec* dettyp)
  :
  TObject(),
  fRecoParam(0),
  fDetTypeRec( dettyp ),
  fNModules(AliITSgeomTGeo::GetNModules()),
  fLastSPD1(AliITSgeomTGeo::GetModuleIndex(2,1,1)-1),
  fNySPD(256),
  fNzSPD(160),
  fNzBins( fNzSPD + 2 ),
  fNyBins( fNySPD + 2),
  fMaxBin( fNzBins * fNyBins ),
  fYpitchSPD(0.0050),
  fZ1pitchSPD(0.0425),
  fZ2pitchSPD(0.0625),
  fHwSPD(0.64),
  fHlSPD(3.48),
  fNSignals(0),
  fSignal2Bin(0),
  fBin2Signal(0)
{

  //Default constructor

  fRecoParam = (AliITSRecoParam*) AliITSReconstructor::GetRecoParam();
  if(!fRecoParam){
    fRecoParam = AliITSRecoParam::GetHighFluxParam();
    AliWarning("Using default AliITSRecoParam class");
  }

 //
 // Initialisation of ITS geometry
 //
  Int_t mmax=AliITSgeomTGeo::GetNModules();
  for (Int_t m=0; m<mmax; m++) {
    Int_t lay,lad,det; AliITSgeomTGeo::GetModuleId(m,lay,lad,det);
    fNdet[m] = (lad-1)*AliITSgeomTGeo::GetNDetectors(lay) + (det-1);
    fNlayer[m] = lay-1;
  }

  fYSPD[0]=0.5*fYpitchSPD;
  for (Int_t m=1; m<fNySPD; m++) fYSPD[m]=fYSPD[m-1]+fYpitchSPD; 
  fZSPD[0]=fZ1pitchSPD;
  for (Int_t m=1; m<fNzSPD; m++) {
    Double_t dz=fZ1pitchSPD;
    if (m==31 || m==32 || m==63  || m==64  || m==95 || m==96 || 
        m==127 || m==128) dz=fZ2pitchSPD; 
    fZSPD[m]=fZSPD[m-1]+dz;
  }
  for (Int_t m=0; m<fNzSPD; m++) {
    Double_t dz=0.5*fZ1pitchSPD;
    if (m==31 || m==32 || m==63  || m==64  || m==95 || m==96 || 
        m==127 || m==128) dz=0.5*fZ2pitchSPD; 
    fZSPD[m]-=dz;
  }

  fSignal2Bin = new UShort_t[fMaxBin];
  fBin2Signal = new UShort_t [fMaxBin];  
  for( int i=0; i<fMaxBin; i++ ) fBin2Signal[i] = 0;

}

AliHLTITSClusterFinderSPD::AliHLTITSClusterFinderSPD( const AliHLTITSClusterFinderSPD &)
  :
  TObject(),
  fRecoParam(0),
  fDetTypeRec(0),
  fNModules(0),
  fLastSPD1(0),
  fNySPD(256),
  fNzSPD(160),
  fNzBins( fNzSPD + 2 ),
  fNyBins( fNySPD + 2),
  fMaxBin( fNzBins * fNyBins ),
  fYpitchSPD(0.0050),
  fZ1pitchSPD(0.0425),
  fZ2pitchSPD(0.0625),
  fHwSPD(0.64),
  fHlSPD(3.48),
  fNSignals(0),
  fSignal2Bin(0),
  fBin2Signal(0)
{
  // dummy
  for (Int_t m=0; m<2200; m++) {
    fNdet[m] = 0;
    fNlayer[m] = 0;
  }

  for (Int_t m=0; m<fNySPD; m++) fYSPD[m]=0;
  for (Int_t m=0; m<fNzSPD; m++) fZSPD[m]=0;  
}

AliHLTITSClusterFinderSPD &AliHLTITSClusterFinderSPD::operator=( const AliHLTITSClusterFinderSPD &)
{
  // dummy
  return *this;
}

void AliHLTITSClusterFinderSPD::RawdataToClusters(AliRawReader* rawReader, std::vector<AliITSRecPoint> & clusters){
  //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------
  rawReader->Reset();
  AliITSRawStreamSPD inputSPD(rawReader);
  FindClustersSPD(&inputSPD, clusters);
}

//__________________________________________________________________________
void AliHLTITSClusterFinderSPD::FindClustersSPD(AliITSRawStreamSPD* input, 
						  std::vector<AliITSRecPoint> & clusters) 
{
  //------------------------------------------------------------
  // SPD cluster finder for raw data (this method is called once per event)
  //------------------------------------------------------------
  
  fNSignals = 0;

  // read raw data input stream
  while (kTRUE) {
    Bool_t next = input->Next();
    if (!next || input->IsNewModule()) {
      Int_t iModule = input->GetPrevModuleID();
      // when all data from a module was read, search for clusters
      if ( fNSignals > 0) { 
	ClustersSPD( clusters, iModule );	
	for( int i=0; i<fNSignals; i++ ){
	  fBin2Signal[fSignal2Bin[i]] = 0;
	}
	fNSignals = 0;
      }
      if (!next) break;      
    }
    // fill the current digit into the bins array
    UShort_t  index = (UShort_t ) ( (input->GetCoord2()+1) * fNzBins + (input->GetCoord1()+1) );
    fSignal2Bin[fNSignals] = index;
    fBin2Signal[index] = fNSignals+1;
    fNSignals++;  
  }  
  //cout<<clusters.size()<<endl;
}

void AliHLTITSClusterFinderSPD::FindCluster(Int_t k,Int_t &n,Int_t *idx) {
  //------------------------------------------------------------
  // returns an array of indices of digits belonging to the cluster
  // (needed when the segmentation is not regular) 
  //------------------------------------------------------------
  if (n<200) idx[n++]=fBin2Signal[k];
  fBin2Signal[k]=0;

  if (fBin2Signal[k-fNzBins] ) FindCluster(k-fNzBins,n,idx);
  if (fBin2Signal[k-1   ] ) FindCluster(k-1 ,n,idx);
  if (fBin2Signal[k+fNzBins] ) FindCluster(k+fNzBins,n,idx);
  if (fBin2Signal[k+1   ] ) FindCluster(k+1   ,n,idx);
}


//__________________________________________________________________________
Int_t AliHLTITSClusterFinderSPD::ClustersSPD( std::vector<AliITSRecPoint> & clusters, Int_t iModule ){

  //Cluster finder for SPD (from digits and from rawdata)

  const TGeoHMatrix *mT2L=AliITSgeomTGeo::GetTracking2LocalMatrix(iModule);

   if (fRecoParam->GetSPDRemoveNoisyFlag()) {
    // Loop on noisy pixels and reset them
    AliITSCalibrationSPD *cal =  
      (AliITSCalibrationSPD*) fDetTypeRec->GetCalibrationModel(iModule);
    for(Int_t ipix = 0; ipix<cal->GetNrBad(); ipix++){
      Int_t row, col;
      cal->GetBadPixel(ipix,row,col);
      Int_t index = (row+1) * fNzBins + (col+1);      
      fBin2Signal[index] = 0;
    }
  }
  
   if (fRecoParam->GetSPDRemoveDeadFlag()) {
     // Loop on dead pixels and reset them
     AliITSCalibrationSPD *cal =  
       (AliITSCalibrationSPD*) fDetTypeRec->GetSPDDeadModel(iModule); 
     if (cal->IsBad()) return 0; // if all ladder is dead, return to save time
     for(Int_t ipix = 0; ipix<cal->GetNrBad(); ipix++){
       Int_t row, col;
       cal->GetBadPixel(ipix,row,col);
       Int_t index = (row+1) * fNzBins + (col+1);
       fBin2Signal[index] = 0;
     }
   }
  
  Int_t nclu=0;

  for(Int_t is =0; is<fNSignals; is++ ){
    int iBin = fSignal2Bin[is];
    if(fBin2Signal[iBin]==0) continue;
    Int_t nBins = 0; 
    Int_t idxBins[200];
    Int_t idxSignals[200];
    FindCluster(iBin,nBins,idxSignals );
    if (nBins >= 199 ){
      //Error("ClustersSPD","SPD Too big cluster !\n"); 
      continue;
    }
    idxBins[0] = 0;
    for( int i=0; i<nBins; i++ ) idxBins[i] = fSignal2Bin[idxSignals[i]-1];
    
    Int_t ymin,ymax,zmin,zmax;
    ymin = (idxBins[0] / fNzBins) - 1;
    ymax = ymin;
    zmin = (idxBins[0] % fNzBins) - 1;
    zmax = zmin;
    
    for (Int_t idx = 0; idx < nBins; idx++) {
      Int_t iy;
      Int_t iz; 
      iy  = (idxBins[idx] / fNzBins) - 1;
      iz  = (idxBins[idx] % fNzBins) - 1;      
      if (ymin > iy) ymin = iy;
      if (ymax < iy) ymax = iy;
      if (zmin > iz) zmin = iz;
      if (zmax < iz) zmax = iz;
    }
    
    Int_t idy =0; //max 2 clusters
    if((iModule <= fLastSPD1) &&idy<3) idy=3;
    if((iModule > fLastSPD1) &&idy<4) idy=4;
    Int_t idz=3;

    // Switch the unfolding OFF/ON
    if(!fRecoParam->GetUseUnfoldingInClusterFinderSPD()) {
      idy=ymax-ymin+1;
      idz=zmax-zmin+1;
    }
    
    for(Int_t iiz=zmin; iiz<=zmax;iiz+=idz){
      for(Int_t iiy=ymin;iiy<=ymax;iiy+=idy){
	Int_t ndigits=0;
	Float_t y=0.,z=0.,q=0.;
	for(Int_t idx=0;idx<nBins;idx++){
	  Int_t iy;
	  Int_t iz; 
	  iy  = (idxBins[idx] / fNzBins)-1;
	  iz  = (idxBins[idx] % fNzBins)-1;
	  if(zmax-zmin>=idz || ymax-ymin>=idy){
	    if(TMath::Abs(iy-iiy)>0.75*idy) continue;
	    if(TMath::Abs(iz-iiz)>0.75*idz) continue;
	  }
	  ndigits++;
	  Float_t qBin=0.;
	  qBin = 1;
	  y+= qBin * fYSPD[iy];
	  z+= qBin * fZSPD[iz];
	  q+= qBin;	
	}// for idx
	if(ndigits==0) continue;
         
	y /= q;
	z /= q;
	y -= fHwSPD;
	z -= fHlSPD;

	Float_t hit[6]; //y,z,sigma(y)^2, sigma(z)^2, charge
        {
        Double_t loc[3]={y,0.,z},trk[3]={0.,0.,0.};
        mT2L->MasterToLocal(loc,trk);
        hit[0]=trk[1];
        hit[1]=trk[2];
	}
	hit[2] = fYpitchSPD*fYpitchSPD/12.;
	hit[3] = fZ1pitchSPD*fZ1pitchSPD/12.;
	hit[4] = 1.;
	hit[5] = 0;

	Int_t info[3] = {ymax-ymin+1,zmax-zmin+1,fNlayer[iModule]};
	Int_t label[4]={-2,-2,-2, fNdet[iModule] };
	AliITSRecPoint cl(label, hit,info);
	cl.SetType(nBins);
	clusters.push_back(cl);       
	nclu++;
      }// for iiy
    }// for iiz
  }//end for iBin
  return nclu;  
}
