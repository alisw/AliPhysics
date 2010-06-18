
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTCALOCLUSTERDATASTRUCT_H
#define ALIHLTCALOCLUSTERDATASTRUCT_H

/**
 * Calo cluster struct for  HLT
 *
 * @file   AliHLTCaloClusterDataStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Calo cluster struct for HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliPID.h"
#include "TArrayI.h"
#include "Rtypes.h"

/**
 * @struct AliHLTCaloClusterHeaderStruct
 * Calorimeter cluster header describing the number of 
 * clusters in the following block
 *
 * @ingroup alihlt_phos
 */


enum ESDClu_t { 	kUndef = -2,
	kPHOSCluster,
	kEMCALPseudoCluster,
	kEMCALClusterv1
};	

struct AliHLTCaloClusterHeaderStruct
{
  Short_t fNClusters;
  Short_t fNDigits;
};

struct AliHLTCaloCellDataStruct
{
   Short_t fCellsAbsId;
   Float_t fCellsAmpFraction;
};


/**
 * @struct AliHLTCaloClusterDataStruct
 * Calorimeter cluster data struct for  HLT
 * Similar to the AliESDCaloCluster class
 * @ingroup alihlt_phos
 */

struct AliHLTCaloClusterDataStruct
{
   
   /** Set the ID */
  void SetID(Int_t id) {fID = id;}      //COMMENT
  
  /** Get the ID */
  Int_t GetID() const {return fID;}  //COMMENT
  
  /** Get the cluster type */
  void SetClusterType(Int_t type) { fClusterType = type; }  //COMMENT
  
  /** Set the cluster type */
  Char_t GetClusterType() const {return fClusterType; } //COMMENT

   /** Is it an EMCAL cluster? */
  Bool_t IsEMCAL() const {return (fClusterType == kEMCALClusterv1);}  //COMMENT
  
  /** Is it a PHOS cluster */
  Bool_t IsPHOS() const {return (fClusterType == kPHOSCluster);} //COMMENT

  /** Set the global postion */
  void SetPosition(const Float_t *pos) {
    fGlobalPos[0] = pos[0]; fGlobalPos[1] = pos[1]; fGlobalPos[2] = pos[2];
  }
  
  /** Get the global position */
  void GetPosition(Float_t *pos) const {
    pos[0] = fGlobalPos[0]; pos[1] = fGlobalPos[1]; pos[2] = fGlobalPos[2];
  }
  /** Set the energy */
  void SetE(Float_t ene) { fEnergy = ene;} //COMMENT
  
  /** Get the energy */
  Double_t E() const   { return fEnergy;} //COMMENT

   /** Set the cluster dispersion */
  void SetClusterDisp(Float_t disp)  { fDispersion = disp; } //COMMENT
  
  /** Get the cluster dispersion */
  Double_t GetClusterDisp() const     { return fDispersion; } //COMMENT

  /** Set the cluster chi2 */
  void SetClusterChi2(Float_t chi2)  { fChi2 = chi2; } //COMMENT
  
  /** Get the cluster chi2 */
  Double_t GetClusterChi2() const     { return fChi2; } //COMMENT

  /** Set the PID data */
  void SetPid(const Float_t *p) 
  {
   // Sets the probability of each particle type
  // Copied from AliESDCaloCluster
  // This function copies "n" PID weights from "scr" to "dest"
  // and normalizes their sum to 1 thus producing conditional
  // probabilities.
  // The negative weights are set to 0.
  // In case all the weights are non-positive they are replaced by
  // uniform probabilities

  Int_t n = AliPID::kSPECIESN;

  Float_t uniform = 1./(Float_t)n;

  Float_t sum = 0;
  for (Int_t i=0; i<n; i++)
    if (p[i]>=0) {
      sum+=p[i];
      fPID[i] = p[i];
    }
    else {
      fPID[i] = 0;
    }

  if(sum>0)
    for (Int_t i=0; i<n; i++) fPID[i] /= sum;
  else
    for (Int_t i=0; i<n; i++) fPID[i] = uniform;

}

  /** Get the PID */
  Float_t *GetPid() {return fPID;} //COMMENT

   /** Set the M20 */
  void SetM20(Float_t m20)                { fM20 = m20; } //COMMENT
  
  /** Get the M20 */
  Double_t GetM20() const                  { return fM20; } //COMMENT

  /** Set the the M02 */
  void SetM02(Float_t m02)                { fM02 = m02; } //COMMENT
  
  /** Get the M02  */
  Double_t GetM02() const                  { return fM02; } //COMMENT

  /** Set number of ex-maxima */
  void SetNExMax(UChar_t nExMax)         { fNExMax = nExMax; } //COMMENT
  
  /** Get the number of ex maxima */
  UChar_t GetNExMax() const              { return fNExMax; } //COMMENT

  /** Set the EMC CPV distance */
  void SetEmcCpvDistance(Float_t dEmcCpv) { fEmcCpvDistance = dEmcCpv; } //COMMENT
  
  /** Get the EMC CPV distance */
  Double_t GetEmcCpvDistance() const       { return fEmcCpvDistance; } //COMMENT
  
  /** Set the distance to track in x and z dimensions */
  void SetTrackDistance(Double_t dx, Double_t dz){fTrackDx=dx; fTrackDz=dz;}
  
  /** Get the distance to track in x */
  Double_t GetTrackDx(void)const {return fTrackDx;}  //COMMENT
  
  /** Get the distance to track in z */
  Double_t GetTrackDz(void)const {return fTrackDz;} //COMMENT
  
  /** Set the distance to closest bad channel */
  void SetDistanceToBadChannel(Float_t dist) {fDistToBadChannel=dist;}
  
  /** Get the distance to closest bad channel */
  Double_t GetDistanceToBadChannel() const {return fDistToBadChannel;}

  /** Set the TOF */
  void SetTOF(Double_t tof) { fTOF = tof; } //COMMENT
  
  /** Geth the TOF */
  Double_t GetTOF() const { return fTOF; } //COMMENT
  
  /** Add an array of tracks */
   void AddTracksMatched(TArrayI & array)  
   { 
      fNTracksMatched = array.GetSize();
       for(Int_t t = 0; (t < fNTracksMatched) && (t < 10); t++) //TODO: remove hard coded 10
       {
 	 fTracksMatched[t] = array[t];
       }
   }
   
  /** 
  * Get the array of the matched tracks 
  */	
   Int_t * GetTracksMatched()  
   {
      return fTracksMatched;
   }
   
   /** Get the best match */
   Int_t GetTrackMatched() const   
  {
    if( fTracksMatched[0] >0)  return  fTracksMatched[0]; 
    else return -1;
  } //Most likely the track associated to the cluster

/** Get the number of tracks matched */
  Int_t GetNTracksMatched() const 
  { 
    for(int i = 0; i < 10; i++) 
      {
        if (fTracksMatched[i] < 0)
         {
     	    return i;
         }
    }
    return 10;
  }

  Int_t GetNCells() const
  {
    return ( Int_t ) (fNCells);
  }

  /** Number of cells in the cluster */
  UInt_t fNCells;                                //COMMENT

  /** Global position */
  Float_t fGlobalPos[3];                      //COMMENT

  /** The total energy of the cell */
  Float_t fEnergy;                            //COMMENT

  /** The time of flight */
  Float_t fTOF;                               //COMMENT

  /** Dispersion */
  Float_t fDispersion;                        //COMMENT
  
  /** Chi2 */
  Float_t fChi2; //COMMENT

  /** Quality of cluster fit */
  Float_t fFitQuality;                        //COMMENT

  /** Second moment along the main eigen axis */
  Float_t fM20;                               //COMMENT

  /** Second moment along the second eigen axis */ 
  Float_t fM02;                               //COMMENT

  /** Distance to closest CPV rec point */
  Float_t fEmcCpvDistance;                    //COMMENT

  /** Distance to nearest bad channel */
  Float_t fDistToBadChannel;                  //COMMENT

  /** Distance to closest track in x direction */
   Float_t fTrackDx; 					//COMMENT

  /** Distance to closest track in z direction */
  Float_t fTrackDz; //COMMENT

  /** PID */
  Float_t fPID[AliPID::kSPECIESN];            //COMMENT

  /** Unique ID of the cluster*/
  Int_t fID;                                     //COMMENT

  /** Number of (Ex) Maxima */
  UChar_t fNExMax;                               //COMMENT 

  /** Flag for differtent cluster type/versions */
  Char_t fClusterType;                           //COMMENT

  /** Distance to nearest bad channel */
  Float_t fDistanceToBadChannel;              //COMMENT

  /** Number of matched tracks */
  Int_t fNTracksMatched; //COMMENT Obsolete?

  /** the matced tracks */
  Int_t fTracksMatched[10];           //COMMENT TODO: remove hardcoded 10

  /** The module */
  Int_t fModule; //COMMENT

  /** Struct containing cell ID and amplitude fraction for the cells */
  AliHLTCaloCellDataStruct fCaloCells;   //COMMENT
  
};


#endif
