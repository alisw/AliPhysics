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

/* $Id$*/

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton  
// EMCAL consists of layers of scintillator and lead
// with scintillator fiber arranged as "shish-kebab" skewers 
// Places the the Barrel Geometry of The EMCAL at Midrapidity
// between 80 and 180(or 190) degrees of Phi and
// -0.7 to 0.7 in eta 
//
//     EMCAL geometry tree:
//     EMCAL -> superModule -> module -> tower(cell)
//     Indexes
//     absId -> nSupMod     -> nModule -> (nIphi,nIeta)
//
//   Name choices: 
//   EMCAL_PDC06 (geometry used for PDC06 simulations, kept for backward compatibility)
//      = equivalent to SHISH_77_TRD1_2X2_FINAL_110DEG in old notation
//   EMCAL_COMPLETE (geometry for expected complete detector)
//      = equivalent to SHISH_77_TRD1_2X2_FINAL_110DEG scTh=0.176 pbTh=0.144
//          in old notation
//   EMCAL_WSUC (Wayne State test stand)
//      = no definite equivalent in old notation, was only used by
//          Aleksei, but kept for testing purposes
//
//   etc.
//
//
//
//*-- Author: Sahal Yacoob (LBL / UCT)
//     and  : Yves Schutz (SUBATECH)
//     and  : Jennifer Klay (LBL)
//     and  : Aleksei Pavlinov (WSU) 
//

//--- Root header files ---
#include <TVector2.h>
#include <TVector3.h>
//-- ALICE Headers.
#include "AliLog.h"

// // --- EMCAL headers
#include "AliEMCALGeometry.h"
#include "AliEMCALShishKebabTrd1Module.h"
#include "AliEMCALRecPoint.h"
//#include "AliEMCALHistoUtilities.h"

ClassImp(AliEMCALGeometry)

// these initialisations are needed for a singleton
AliEMCALGeometry  *AliEMCALGeometry::fgGeom      = 0;
const Char_t*      AliEMCALGeometry::fgkDefaultGeometryName = "EMCAL_COMPLETE";
//
// Usage: 
//        You can create the AliEMCALGeometry object independently from anything.
//        You have to use just the correct name of geometry. If name is empty string the
//        default name of geometry will be used.
//         
//  AliEMCALGeometry* g = AliEMCALGeometry::GetInstance(name,title); // first time
//  ..
//  g = AliEMCALGeometry::GetInstance();                             // after first time
//
//  MC:   If you work with MC data you have to get geometry the next way: 
//  ==                                      =============================
//  AliRunLoader    *rl   = AliRunLoader::Instance();
//  AliEMCALGeometry *geom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
//  TGeoManager::Import("geometry.root");

AliEMCALGeometry::AliEMCALGeometry() 
  : AliEMCALGeoUtils()
{ 
  // default ctor only for internal usage (singleton)
  // must be kept public for root persistency purposes, 
  // but should never be called by the outside world    

  AliDebug(2, "AliEMCALGeometry : default ctor ");
}
//______________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry(const Text_t* name, const Text_t* title) 
  : AliEMCALGeoUtils(name, title)
{
  // ctor only for internal usage (singleton)
  AliDebug(2, Form("AliEMCALGeometry(%s,%s) ", name,title));

}
//______________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry(const AliEMCALGeometry& geom)
  : AliEMCALGeoUtils(geom)
{
  //copy ctor
}

//______________________________________________________________________
AliEMCALGeometry::~AliEMCALGeometry(void){
    // dtor
}


//______________________________________________________________________
AliEMCALGeometry *  AliEMCALGeometry::GetInstance(){ 
  // Returns the pointer of the unique instance
  
  AliEMCALGeometry * rv = static_cast<AliEMCALGeometry *>( fgGeom );
  return rv; 
}

//______________________________________________________________________
AliEMCALGeometry* AliEMCALGeometry::GetInstance(const Text_t* name,
						const Text_t* title){
    // Returns the pointer of the unique instance

    AliEMCALGeometry * rv = 0; 
    if ( fgGeom == 0 ) {
      if ( strcmp(name,"") == 0 ) { // get default geometry
	 fgGeom = new AliEMCALGeometry(fgkDefaultGeometryName, title);
      } else {
	 fgGeom = new AliEMCALGeometry(name, title);
      }  // end if strcmp(name,"")
      if ( AliEMCALEMCGeometry::fgInit ) rv = (AliEMCALGeometry * ) fgGeom;
      else {
	 rv = 0; 
	 delete fgGeom; 
	 fgGeom = 0; 
      } // end if fgInit
    }else{
	if ( strcmp(fgGeom->GetName(), name) != 0) {
	  printf("\ncurrent geometry is %s : ", fgGeom->GetName());
	  printf(" you cannot call %s ",name);  
	}else{
	  rv = (AliEMCALGeometry *) fgGeom; 
	} // end 
    }  // end if fgGeom
    return rv; 
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, Double_t distEff, Double_t &xr, Double_t &yr, Double_t &zr) const
{
  // Jul 30, 2007 - taking into account position of shower max
  // Look to see what the relative
  // position inside a given cell is
  // for a recpoint.
  // In:
  // absId   - cell is as in Geant,     0<= absId   < fNCells;
  // e       - cluster energy
  // OUT:
  // xr,yr,zr - x,y,z coordinates of cell with absId inside SM 

  // Shift index taking into account the difference between standard SM 
  // and SM of half size in phi direction
  const  Int_t kphiIndexShift = fCentersOfCellsPhiDir.GetSize()/4; // Nov 22, 2006; was 6 for cas 2X2
  static Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
  static Int_t iphim, ietam;
  static AliEMCALShishKebabTrd1Module *mod = 0;
  static TVector2 v;
  if(!CheckAbsCellId(absId)) return kFALSE;
  
  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  GetModulePhiEtaIndexInSModule(nSupMod, nModule, iphim, ietam);
  GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta); 
  
  //Get eta position. Careful with ALICE conventions (increase index decrease eta)	
  if(nSupMod%2 == 0) {             
    ietam = (fCentersOfCellsEtaDir.GetSize()/2-1)-ietam;// 47-ietam, revert the ordering on A side in order to keep convention.
    if(nIeta == 0) nIeta = 1;
    else	   nIeta = 0;
  }
  mod = GetShishKebabModule(ietam);
  mod ->GetPositionAtCenterCellLine(nIeta, distEff, v); 
  xr = v.Y() - fParSM[0];
  zr = v.X() - fParSM[2];
  
  //Get phi position. Careful with ALICE conventions (increase index increase phi)
  Int_t iphi2 = iphi;
  if(nSupMod<10) { 
    if(nSupMod%2 != 0) 
      iphi2 = (fCentersOfCellsPhiDir.GetSize()-1)-iphi;// 23-iphi, revert the ordering on C side in order to keep convention.
    yr = fCentersOfCellsPhiDir.At(iphi2);
    
  } else {
    if(nSupMod%2 != 0) 
      iphi2 = (fCentersOfCellsPhiDir.GetSize()/2-1)-iphi;// 11-iphi, revert the ordering on C side in order to keep convention.
    yr = fCentersOfCellsPhiDir.At(iphi2 + kphiIndexShift);
  }
  
  AliDebug(1,Form("absId %i nSupMod %i iphi %i ieta %i xr %f yr %f zr %f ",absId,nSupMod,iphi,ieta,xr,yr,zr));
  
  return kTRUE;
}

//Not in use, comment for the moment
//________________________________________________________________________________________________
//Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, Int_t maxAbsId, Double_t distEff, Double_t &xr, Double_t &yr, Double_t &zr) const
//{
//  // Jul 31, 2007 - taking into account position of shower max and apply coor2.
//  // Look to see what the relative
//  // position inside a given cell is
//  // for a recpoint.
//  // In:
//  // absId     - cell is as in Geant,     0<= absId   < fNCells;
//  // maxAbsId  - abs id of cell with highest energy
//  // e         - cluster energy
//  // OUT:
//  // xr,yr,zr - x,y,z coordinates of cell with absId inside SM 
//  
//  // Shift index taking into account the difference between standard SM 
//  // and SM of half size in phi direction
//  const  Int_t kphiIndexShift = fCentersOfCellsPhiDir.GetSize()/4; // Nov 22, 2006; was 6 for cas 2X2
//  static Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
//  static Int_t iphim, ietam;
//  static AliEMCALShishKebabTrd1Module *mod = 0;
//  static TVector2 v;
//
//  static Int_t nSupModM, nModuleM, nIphiM, nIetaM, iphiM, ietaM;
//  static Int_t iphimM, ietamM, maxAbsIdCopy=-1;
//  static AliEMCALShishKebabTrd1Module *modM = 0;
//  static Double_t distCorr;
//
//  if(!CheckAbsCellId(absId)) return kFALSE;
//  
//  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
//  GetModulePhiEtaIndexInSModule(nSupMod, nModule, iphim, ietam);
//  GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta); 
//  
//  //Get eta position. Careful with ALICE conventions (increase index decrease eta)	
//  if(nSupMod%2 == 0) {             
//    ietam = (fCentersOfCellsEtaDir.GetSize()/2-1)-ietam;// 23-ietam, revert the ordering on A side in order to keep convention.
//    if(nIeta == 0) nIeta = 1;
//    else		   nIeta = 0;
//  }
//  
//  mod = GetShishKebabModule(ietam);
//  
//  if(absId != maxAbsId) {
//    distCorr = 0.;
//    if(maxAbsIdCopy != maxAbsId) {
//      GetCellIndex(maxAbsId, nSupModM, nModuleM, nIphiM, nIetaM);
//      GetModulePhiEtaIndexInSModule(nSupModM, nModuleM, iphimM, ietamM);
//      GetCellPhiEtaIndexInSModule(nSupModM,nModuleM,nIphiM,nIetaM, iphiM, ietaM); 
//      //Careful with ALICE conventions (increase index decrease eta)	
//      if(nSupModM%2 == 0) {             
//	ietamM = (fCentersOfCellsEtaDir.GetSize()/2-1)-ietamM;// 47-ietam, revert the ordering on A side in order to keep convention.
//      }
//      
//      modM = GetShishKebabModule(ietamM); // do I need this ?
//      maxAbsIdCopy = maxAbsId;
//    }
//    
//    if(ietamM !=0) {
//      distCorr = fEMCGeometry->GetEtaModuleSize()*(ietam-ietamM)/TMath::Tan(modM->GetTheta()); // Stay here
//      //printf(" distCorr %f | dist %f | ietam %i -> etamM %i\n", distCorr, dist, ietam, ietamM);  
//    }
//    // distEff += distCorr;
//  }
//  // Bad resolution in this case, strong bias vs phi
//  // distEff = 0.0; 
//  mod->GetPositionAtCenterCellLine(nIeta, distEff, v); // Stay here
//  xr = v.Y() - fParSM[0];
//  zr = v.X() - fParSM[2];
//  
//  //Get phi position. Careful with ALICE conventions (increase index increase phi)
//  Int_t iphi2 = iphi;
//  if(nSupMod<10) { 
//    if(nSupMod%2 != 0) 
//      iphi2 = (fCentersOfCellsPhiDir.GetSize()-1)-iphi;// 23-iphi, revert the ordering on C side in order to keep convention.
//    yr = fCentersOfCellsPhiDir.At(iphi2);
//    
//  } else {
//    if(nSupMod%2 != 0) 
//      iphi2 = (fCentersOfCellsPhiDir.GetSize()/2-1)-iphi;// 11-iphi, revert the ordering on C side in order to keep convention.
//    yr = fCentersOfCellsPhiDir.At(iphi2 + kphiIndexShift);
//  }
//  AliDebug(1,Form("absId %i nSupMod %i iphi %i ieta %i xr %f yr %f zr %f ",absId,nSupMod,iphi,ieta,xr,yr,zr));
//  
//  return kTRUE;
//}
//

//
// == Shish-kebab cases ==
//

//
////_________________________________________________________________________________
//void AliEMCALGeometry::GetGlobalEMCAL(const AliEMCALRecPoint *rp, TVector3 &vglob) const
//{
//  // Figure out the global numbering
//  // of a given supermodule from the
//  // local numbering for RecPoints
//
//  static TVector3 vloc;
//  static Int_t nSupMod, nModule, nIphi, nIeta;
//
//  const AliEMCALRecPoint *rpTmp = rp;
//  const AliEMCALRecPoint *rpEmc = rpTmp;
//
//  GetCellIndex(rpEmc->GetAbsId(0), nSupMod, nModule, nIphi, nIeta);
//  rpTmp->GetLocalPosition(vloc);
//  GetGlobal(vloc, vglob, nSupMod);
//}

