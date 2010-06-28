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

/* $Id: AliEMCALGeoUtils.cxx 25590 2008-05-06 07:09:11Z prsnko $ */

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton  
//                  
// -- Author: Magali Estienne (magali.estienne@subatech.in2p3.fr)

//
// Usage: 
//        You can create the AliEMCALGeoUtils object independently from anything.
//        You have to use just the correct name of geometry. If name is empty string the
//        default name of geometry will be used.
//         
//  AliEMCALGeoUtils* geom = new AliEMCALGeoUtils("EMCAL_COMPLETE","EMCAL");
//  TGeoManager::Import("geometry.root");
//
//  MC:   If you work with MC data you have to get geometry the next way: 
//  ==                                      =============================
// !!!!!!!!! This part has to be modified
//  AliRunLoader    *rl   = AliRunLoader::GetRunLoader();
//  AliEMCALEMCGeometry *geom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
//  TGeoManager::Import("geometry.root");


// --- ROOT system ---

#include <TParticle.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TList.h>
#include <TBrowser.h>

// --- Standard library ---
//#include <Riostream.h>

// --- AliRoot header files ---
#include "AliEMCALGeoUtils.h"
#include "AliEMCALShishKebabTrd1Module.h"

ClassImp(AliEMCALGeoUtils)

//____________________________________________________________________________
AliEMCALGeoUtils::AliEMCALGeoUtils():
  fEMCGeometry(0x0),fGeoName(0),
  fKey110DEG(0),fNCellsInSupMod(0),fNETAdiv(0),fNPHIdiv(0),
  fNCellsInModule(0),fPhiBoundariesOfSM(0x0),fPhiCentersOfSM(0x0),
  fPhiCentersOfCells(0x0),fCentersOfCellsEtaDir(0x0),
  fCentersOfCellsPhiDir(0x0),fEtaCentersOfCells(0x0),
  fNCells(0),fNPhi(0),fCentersOfCellsXDir(0x0),fArm1EtaMin(0),
  fArm1EtaMax(0),fArm1PhiMin(0),fArm1PhiMax(0),fEtaMaxOfTRD1(0),
  fShishKebabTrd1Modules(0),fParSM(0x0),fPhiModuleSize(0.),
  fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fNZ(0),
  fIPDistance(0.),fLongModuleSize(0.),fShellThickness(0.),
  fZLength(0.),fSampling(0.)
{
  // default ctor 
  // must be kept public for root persistency purposes, but should never be called by the outside world
  fEnvelop[0] = 0.;
  fEnvelop[1] = 0.;
  fEnvelop[2] = 0.;

}  

//____________________________________________________________________________
AliEMCALGeoUtils::AliEMCALGeoUtils(const AliEMCALGeoUtils & geo)
  : TNamed(geo),
    fEMCGeometry(geo.fEMCGeometry),fGeoName(geo.fGeoName),
    fKey110DEG(geo.fKey110DEG),fNCellsInSupMod(geo.fNCellsInSupMod),fNETAdiv(geo.fNETAdiv),fNPHIdiv(geo.fNPHIdiv),
    fNCellsInModule(geo.fNCellsInModule),fPhiBoundariesOfSM(geo.fPhiBoundariesOfSM),fPhiCentersOfSM(geo.fPhiCentersOfSM),
    fPhiCentersOfCells(geo.fPhiCentersOfCells),fCentersOfCellsEtaDir(geo.fCentersOfCellsEtaDir),
    fCentersOfCellsPhiDir(geo.fCentersOfCellsPhiDir),fEtaCentersOfCells(geo.fEtaCentersOfCells),
    fNCells(geo.fNCells),fNPhi(geo.fNPhi),fCentersOfCellsXDir(geo.fCentersOfCellsXDir),fArm1EtaMin(geo.fArm1EtaMin),
    fArm1EtaMax(geo.fArm1EtaMax),fArm1PhiMin(geo.fArm1PhiMin),fArm1PhiMax(geo.fArm1PhiMax),fEtaMaxOfTRD1(geo.fEtaMaxOfTRD1),
    fShishKebabTrd1Modules(geo.fShishKebabTrd1Modules),fParSM(geo.fParSM),fPhiModuleSize(geo.fPhiModuleSize),
    fEtaModuleSize(geo.fEtaModuleSize),fPhiTileSize(geo.fPhiTileSize),fEtaTileSize(geo.fEtaTileSize),fNZ(geo.fNZ),
    fIPDistance(geo.fIPDistance),fLongModuleSize(geo.fLongModuleSize),fShellThickness(geo.fShellThickness),
    fZLength(geo.fZLength),fSampling(geo.fSampling)
{
  fEnvelop[0] = geo.fEnvelop[0];
  fEnvelop[1] = geo.fEnvelop[1];
  fEnvelop[2] = geo.fEnvelop[2];
}

//____________________________________________________________________________
AliEMCALGeoUtils::AliEMCALGeoUtils(const Text_t* name, const Text_t* title) 
  : TNamed(name, title),
    fEMCGeometry(0x0),fGeoName(0),
    fKey110DEG(0),fNCellsInSupMod(0),fNETAdiv(0),fNPHIdiv(0),
    fNCellsInModule(0),fPhiBoundariesOfSM(0x0),fPhiCentersOfSM(0x0),
    fPhiCentersOfCells(0x0),fCentersOfCellsEtaDir(0x0),
    fCentersOfCellsPhiDir(0x0),fEtaCentersOfCells(0x0),
    fNCells(0),fNPhi(0),fCentersOfCellsXDir(0x0),fArm1EtaMin(0),
    fArm1EtaMax(0),fArm1PhiMin(0),fArm1PhiMax(0),fEtaMaxOfTRD1(0),
    fShishKebabTrd1Modules(0),fParSM(0x0),fPhiModuleSize(0.),
    fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fNZ(0),
    fIPDistance(0.),fLongModuleSize(0.),fShellThickness(0.),
    fZLength(0.),fSampling(0.)
{ 

  // ctor only for normal usage 

  fEMCGeometry = new AliEMCALEMCGeometry(name,title);

  fGeoName = fEMCGeometry->GetGeoName();
  fKey110DEG = fEMCGeometry->GetKey110DEG();
  fNCellsInSupMod = fEMCGeometry->GetNCellsInSupMod();
  fNETAdiv = fEMCGeometry->GetNETAdiv();
  fNPHIdiv = fEMCGeometry->GetNPHIdiv();
  fNCellsInModule = fNPHIdiv*fNETAdiv;
  static int i;
  Int_t nSMod = fEMCGeometry->GetNumberOfSuperModules();
  fPhiBoundariesOfSM.Set(nSMod);
  fPhiCentersOfSM.Set(nSMod/2);
  for(Int_t sm=0; sm<nSMod; sm++) {
    i = sm/2;
    fEMCGeometry->GetPhiBoundariesOfSM(sm,fPhiBoundariesOfSM[2*i],fPhiBoundariesOfSM[2*i+1]);
  }

  Double_t phiMin =  0.;
  Double_t phiMax =  0.;
  for(Int_t sm=0; sm<nSMod; sm++) {
    fEMCGeometry->GetPhiBoundariesOfSM(sm,phiMin,phiMax);
	i=sm/2;
    fPhiCentersOfSM[i] = fEMCGeometry->GetPhiCenterOfSM(sm);
  }
  fNCells = fEMCGeometry->GetNCells();
  fNPhi = fEMCGeometry->GetNPhi();
  fEnvelop[0] = fEMCGeometry->GetEnvelop(0);
  fEnvelop[1] = fEMCGeometry->GetEnvelop(1);
  fEnvelop[2] = fEMCGeometry->GetEnvelop(2);
  fArm1EtaMin = fEMCGeometry->GetArm1EtaMin();
  fArm1EtaMax = fEMCGeometry->GetArm1EtaMax();
  fArm1PhiMin = fEMCGeometry->GetArm1PhiMin();
  fArm1PhiMax = fEMCGeometry->GetArm1PhiMax();
  fShellThickness = fEMCGeometry->GetShellThickness();
  fZLength = fEMCGeometry->GetZLength();
  fSampling = fEMCGeometry->GetSampling();
  fParSM = fEMCGeometry->GetSuperModulesPars();
  fEtaModuleSize = fEMCGeometry->GetEtaModuleSize();
  fPhiModuleSize = fEMCGeometry->GetPhiModuleSize();
  fEtaTileSize = fEMCGeometry->GetEtaTileSize();
  fPhiTileSize = fEMCGeometry->GetPhiTileSize();
  fNZ = fEMCGeometry->GetNZ();
  fIPDistance = fEMCGeometry->GetIPDistance();
  fLongModuleSize = fEMCGeometry->GetLongModuleSize();

  CreateListOfTrd1Modules();

  for(Int_t smod=0; smod < fEMCGeometry->GetNumberOfSuperModules(); smod++)
		fkSModuleMatrix[smod]=0 ;	
	
  if (AliDebugLevel()>=2) {
    fEMCGeometry->Print();
    PrintGeometry();
  }

}

//____________________________________________________________________________
AliEMCALGeoUtils & AliEMCALGeoUtils::operator = (const AliEMCALGeoUtils  & /*rvalue*/) { 

  Fatal("assignment operator", "not implemented") ; 
    return *this ;
}

//____________________________________________________________________________
AliEMCALGeoUtils::~AliEMCALGeoUtils(void)
{
  // dtor
  if(fEMCGeometry){
    delete fEMCGeometry; fEMCGeometry = 0 ;
  }
}


//________________________________________________________________________________________________
void AliEMCALGeoUtils::Browse(TBrowser* b)
{
  //Browse the modules
  if(fShishKebabTrd1Modules) b->Add(fShishKebabTrd1Modules);
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::IsFolder() const
{
  //Check if fShishKebabTrd1Modules is in folder
  if(fShishKebabTrd1Modules) return kTRUE;
  else                       return kFALSE;
}

//________________________________________________________________________________________________
void AliEMCALGeoUtils::GetGlobal(const Double_t *loc, Double_t *glob, int ind) const
{
  // Figure out the global numbering
  // of a given supermodule from the
  // local numbering and the transformation
  // matrix stored by the geometry manager (allows for misaligned
  // geometry)
	
	const TGeoHMatrix* m = GetMatrixForSuperModule(ind);
    if(m) {
      m->LocalToMaster(loc, glob);
    } else {
      AliFatal("Geo matrixes are not loaded \n") ;
    }
}

//________________________________________________________________________________________________
void AliEMCALGeoUtils::GetGlobal(const TVector3 &vloc, TVector3 &vglob, int ind) const
{
  //Figure out the global numbering
  //of a given supermodule from the
  //local numbering given a 3-vector location

  static Double_t tglob[3], tloc[3];
  vloc.GetXYZ(tloc);
  GetGlobal(tloc, tglob, ind);
  vglob.SetXYZ(tglob[0], tglob[1], tglob[2]);
}

//________________________________________________________________________________________________
void AliEMCALGeoUtils::GetGlobal(Int_t absId , double glob[3]) const
{
  // Alice numbering scheme - Jun 03, 2006
  static Int_t nSupMod, nModule, nIphi, nIeta;
  static double loc[3];

  glob[0]=glob[1]=glob[2]=0.0; // bad case
  if(RelPosCellInSModule(absId, loc)) {
    GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);

	  const TGeoHMatrix* m = GetMatrixForSuperModule(nSupMod);
	  if(m) {
      m->LocalToMaster(loc, glob);
    } else {
      AliFatal("Geo matrixes are not loaded \n") ;
    }
  }
}

//___________________________________________________________________
void AliEMCALGeoUtils::GetGlobal(Int_t absId , TVector3 &vglob) const
{
  // Alice numbering scheme - Jun 03, 2006
  static Double_t glob[3];

  GetGlobal(absId, glob);
  vglob.SetXYZ(glob[0], glob[1], glob[2]);

}


//______________________________________________________________________
void AliEMCALGeoUtils::PrintCellIndexes(Int_t absId, int pri, const char *tit) const
{
  // Service methods
  Int_t nSupMod, nModule, nIphi, nIeta;
  Int_t iphi, ieta;
  TVector3 vg;

  GetCellIndex(absId,  nSupMod, nModule, nIphi, nIeta);
  printf(" %s | absId : %i -> nSupMod %i nModule %i nIphi %i nIeta %i \n", tit, absId,  nSupMod, nModule, nIphi, nIeta);
  if(pri>0) {
    GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
    printf(" local SM index : iphi %i : ieta %i \n", iphi,ieta);
    GetGlobal(absId, vg);
    printf(" vglob : mag %7.2f : perp %7.2f : z %7.2f : eta %6.4f : phi %6.4f(%6.2f) \n", 
	   vg.Mag(), vg.Perp(), vg.Z(), vg.Eta(), vg.Phi(), vg.Phi()*TMath::RadToDeg());
  }
}

//________________________________________________________________________________________________
void AliEMCALGeoUtils::EtaPhiFromIndex(Int_t absId,Double_t &eta,Double_t &phi) const
{
  // Nov 16, 2006- float to double
  // version for TRD1 only
  static TVector3 vglob;
  GetGlobal(absId, vglob);
  eta = vglob.Eta();
  phi = vglob.Phi();
}

//________________________________________________________________________________________________
void AliEMCALGeoUtils::EtaPhiFromIndex(Int_t absId,Float_t &eta,Float_t &phi) const
{
  // Nov 16,2006 - should be discard in future
  static TVector3 vglob;
  GetGlobal(absId, vglob);
  eta = float(vglob.Eta());
  phi = float(vglob.Phi());
}

//
// == Shish-kebab cases ==
//
//________________________________________________________________________________________________
Int_t AliEMCALGeoUtils::GetAbsCellId(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta) const
{ 
  // 27-aug-04; 
  // corr. 21-sep-04; 
  //       13-oct-05; 110 degree case
  // May 31, 2006; ALICE numbering scheme:
  // 0 <= nSupMod < fNumberOfSuperModules
  // 0 <= nModule  < fNPHI * fNZ ( fNPHI * fNZ/2 for fKey110DEG=1)
  // 0 <= nIphi   < fNPHIdiv
  // 0 <= nIeta   < fNETAdiv
  // 0 <= absid   < fNCells
  static Int_t id=0; // have to change from 0 to fNCells-1
  if(fKey110DEG == 1 && nSupMod >= 10) { // 110 degree case; last two supermodules
    id  = fNCellsInSupMod*10 + (fNCellsInSupMod/2)*(nSupMod-10);
  } else {
    id  = fNCellsInSupMod*nSupMod;
  }
  id += fNCellsInModule *nModule;
  id += fNPHIdiv *nIphi;
  id += nIeta;
  if(id<0 || id >= fNCells) {
//     printf(" wrong numerations !!\n");
//     printf("    id      %6i(will be force to -1)\n", id);
//     printf("    fNCells %6i\n", fNCells);
//     printf("    nSupMod %6i\n", nSupMod);
//     printf("    nModule  %6i\n", nModule);
//     printf("    nIphi   %6i\n", nIphi);
//     printf("    nIeta   %6i\n", nIeta);
    id = -TMath::Abs(id); // if negative something wrong
  }
  return id;
}

//________________________________________________________________________________________________
void  AliEMCALGeoUtils::GetModuleIndexesFromCellIndexesInSModule(Int_t nSupMod, Int_t iphi, Int_t ieta, 
			Int_t &iphim, Int_t &ietam, Int_t &nModule) const
{
  // Transition from cell indexes (ieta,iphi) to module indexes (ietam,iphim, nModule)
  static Int_t nphi;
  nphi  = GetNumberOfModuleInPhiDirection(nSupMod);  

  ietam  = ieta/fNETAdiv;
  iphim  = iphi/fNPHIdiv;
  nModule = ietam * nphi + iphim; 
}

//________________________________________________________________________________________________
Int_t  AliEMCALGeoUtils::GetAbsCellIdFromCellIndexes(Int_t nSupMod, Int_t iphi, Int_t ieta) const
{
  // Transition from super module number(nSupMod) and cell indexes (ieta,iphi) to absId
  static Int_t ietam, iphim, nModule;
  static Int_t nIeta, nIphi; // cell indexes in module

  GetModuleIndexesFromCellIndexesInSModule(nSupMod, iphi, ieta, ietam, iphim, nModule);

  nIeta = ieta%fNETAdiv;
  nIeta = fNETAdiv - 1 - nIeta;
  nIphi = iphi%fNPHIdiv;

  return GetAbsCellId(nSupMod, nModule, nIphi, nIeta);
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::SuperModuleNumberFromEtaPhi(Double_t eta, Double_t phi, Int_t &nSupMod) const
{ 
  // Return false if phi belongs a phi cracks between SM
 
  static Int_t i;

  if(TMath::Abs(eta) > fEtaMaxOfTRD1) return kFALSE;

  phi = TVector2::Phi_0_2pi(phi); // move phi to (0,2pi) boundaries
  for(i=0; i<6; i++) {
	
	//Check if it is not the complete geometry
	if (i >= fEMCGeometry->GetNumberOfSuperModules()/2) return kFALSE;

    if(phi>=fPhiBoundariesOfSM[2*i] && phi<=fPhiBoundariesOfSM[2*i+1]) {
      nSupMod = 2*i;
      if(eta < 0.0) nSupMod++;
      AliDebug(1,Form("eta %f phi %f(%5.2f) : nSupMod %i : #bound %i", eta,phi,phi*TMath::RadToDeg(), nSupMod,i));
      return kTRUE;
    }
  }
  return kFALSE;
}


//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::GetAbsCellIdFromEtaPhi(Double_t eta, Double_t phi, Int_t &absId) const
{
  // Nov 17,2006
  // stay here - phi problem as usual 
  static Int_t nSupMod, i, ieta, iphi, etaShift, nphi;
  static Double_t absEta=0.0, d=0.0, dmin=0.0, phiLoc;
  absId = nSupMod = - 1;
  if(SuperModuleNumberFromEtaPhi(eta, phi, nSupMod)) {
    // phi index first
    phi    = TVector2::Phi_0_2pi(phi);
    phiLoc = phi - fPhiCentersOfSM[nSupMod/2];
    nphi   = fPhiCentersOfCells.GetSize();
    if(nSupMod>=10) {
      phiLoc = phi - 190.*TMath::DegToRad();
      nphi  /= 2;
    }

    dmin   = TMath::Abs(fPhiCentersOfCells[0]-phiLoc);
    iphi   = 0;
    for(i=1; i<nphi; i++) {
      d = TMath::Abs(fPhiCentersOfCells[i] - phiLoc);
      if(d < dmin) {
        dmin = d;
        iphi = i;
      }
      //      printf(" i %i : d %f : dmin %f : fPhiCentersOfCells[i] %f \n", i, d, dmin, fPhiCentersOfCells[i]);
    }
    // odd SM are turned with respect of even SM - reverse indexes
    AliDebug(2,Form(" iphi %i : dmin %f (phi %f, phiLoc %f ) ", iphi, dmin, phi, phiLoc));
    // eta index
    absEta   = TMath::Abs(eta);
    etaShift = iphi*fCentersOfCellsEtaDir.GetSize();
    dmin     = TMath::Abs(fEtaCentersOfCells[etaShift]-absEta);
    ieta     = 0;
    for(i=1; i<fCentersOfCellsEtaDir.GetSize(); i++) {
      d = TMath::Abs(fEtaCentersOfCells[i+etaShift] - absEta);
      if(d < dmin) {
        dmin = d;
        ieta = i;
      }
    }
    AliDebug(2,Form(" ieta %i : dmin %f (eta=%f) : nSupMod %i ", ieta, dmin, eta, nSupMod));

    if(eta<0) iphi = (nphi-1) - iphi;
	  
	//patch for mapping following alice convention  
	if(nSupMod%2 == 0)		  
		  ieta = (fCentersOfCellsEtaDir.GetSize()-1)-ieta;// 47-ieta, revert the ordering on A side in order to keep convention.
	else {
		if(nSupMod<10) 
				iphi = (fCentersOfCellsPhiDir.GetSize()-1)  -iphi;// 23-iphi, revert the ordering on C side in order to keep convention.
		else 
				iphi = (fCentersOfCellsPhiDir.GetSize()/2-1)-iphi;// 11-iphi, revert the ordering on C side in order to keep convention.
	}
  
    absId = GetAbsCellIdFromCellIndexes(nSupMod, iphi, ieta);

    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t  AliEMCALGeoUtils::CheckAbsCellId(Int_t absId) const
{ 
  // May 31, 2006; only trd1 now
  if(absId<0 || absId >= fNCells) return kFALSE;
  else                            return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::GetCellIndex(Int_t absId,Int_t &nSupMod,Int_t &nModule,Int_t &nIphi,Int_t &nIeta) const
{ 
  // 21-sep-04; 19-oct-05;
  // May 31, 2006; ALICE numbering scheme:
  // 
  // In:
  // absId   - cell is as in Geant,     0<= absId   < fNCells;
  // Out:
  // nSupMod - super module(SM) number, 0<= nSupMod < fNumberOfSuperModules;
  // nModule  - module number in SM,     0<= nModule  < fNCellsInSupMod/fNCellsInSupMod or(/2) for tow last SM (10th and 11th);
  // nIphi   - cell number in phi driection inside module; 0<= nIphi < fNPHIdiv; 
  // nIeta   - cell number in eta driection inside module; 0<= nIeta < fNETAdiv; 
  // 
  static Int_t tmp=0, sm10=0;
  if(!CheckAbsCellId(absId)) return kFALSE;

  sm10 = fNCellsInSupMod*10;
  if(fKey110DEG == 1 && absId >= sm10) { // 110 degree case; last two supermodules  
    nSupMod = (absId-sm10) / (fNCellsInSupMod/2) + 10;
    tmp     = (absId-sm10) % (fNCellsInSupMod/2);
  } else {
    nSupMod = absId / fNCellsInSupMod;
    tmp     = absId % fNCellsInSupMod;
  }

  nModule  = tmp / fNCellsInModule;
  tmp     = tmp % fNCellsInModule;
  nIphi   = tmp / fNPHIdiv;
  nIeta   = tmp % fNPHIdiv;

  return kTRUE;
}

//________________________________________________________________________________________________
Int_t  AliEMCALGeoUtils::GetSuperModuleNumber(Int_t absId)  const
{
  // Return the number of the  supermodule given the absolute
  // ALICE numbering id

  static Int_t nSupMod, nModule, nIphi, nIeta;
  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  return nSupMod;
} 

//________________________________________________________________________________________________
void AliEMCALGeoUtils::GetModulePhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule,  int &iphim, int &ietam) const
{ 
  // added nSupMod; - 19-oct-05 !
  // Alice numbering scheme        - Jun 01,2006 
  // ietam, iphi - indexes of module in two dimensional grid of SM
  // ietam - have to change from 0 to fNZ-1
  // iphim - have to change from 0 to nphi-1 (fNPhi-1 or fNPhi/2-1)
  static Int_t nphi;

  if(fKey110DEG == 1 && nSupMod>=10) nphi = fNPhi/2;
  else                               nphi = fNPhi;

  ietam = nModule/nphi;
  iphim = nModule%nphi;
}

//________________________________________________________________________________________________
void AliEMCALGeoUtils::GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta, 
int &iphi, int &ieta) const
{ 
  // 
  // Added nSupMod; Nov 25, 05
  // Alice numbering scheme  - Jun 01,2006 
  // IN:
  // nSupMod - super module(SM) number, 0<= nSupMod < fNumberOfSuperModules;
  // nModule  - module number in SM,     0<= nModule  < fNCellsInSupMod/fNCellsInSupMod or(/2) for tow last SM (10th and 11th);
  // nIphi   - cell number in phi driection inside module; 0<= nIphi < fNPHIdiv; 
  // nIeta   - cell number in eta driection inside module; 0<= nIeta < fNETAdiv; 
  // 
 // OUT:
  // ieta, iphi - indexes of cell(tower) in two dimensional grid of SM
  // ieta - have to change from 0 to (fNZ*fNETAdiv-1)
  // iphi - have to change from 0 to (fNPhi*fNPHIdiv-1 or fNPhi*fNPHIdiv/2-1)
  //
  static Int_t iphim, ietam;

  GetModulePhiEtaIndexInSModule(nSupMod,nModule, iphim, ietam); 
  //  ieta  = ietam*fNETAdiv + (1-nIeta); // x(module) = -z(SM) 
  ieta  = ietam*fNETAdiv + (fNETAdiv - 1 - nIeta); // x(module) = -z(SM) 
  iphi  = iphim*fNPHIdiv + nIphi;     // y(module) =  y(SM) 

  if(iphi<0 || ieta<0)
  AliDebug(1,Form(" nSupMod %i nModule %i nIphi %i nIeta %i => ieta %i iphi %i\n", 
  nSupMod, nModule, nIphi, nIeta, ieta, iphi));
}


// Methods for AliEMCALRecPoint - Feb 19, 2006
//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::RelPosCellInSModule(Int_t absId, Double_t &xr, Double_t &yr, Double_t &zr) const
{
  // Look to see what the relative
  // position inside a given cell is
  // for a recpoint.
  // Alice numbering scheme - Jun 08, 2006
  // In:
  // absId   - cell is as in Geant,     0<= absId   < fNCells;
  // OUT:
  // xr,yr,zr - x,y,z coordinates of cell with absId inside SM 

  // Shift index taking into account the difference between standard SM 
  // and SM of half size in phi direction
  const Int_t kphiIndexShift = fCentersOfCellsPhiDir.GetSize()/4; // Nov 22, 2006; was 6 for cas 2X2
  static Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
  if(!CheckAbsCellId(absId)) return kFALSE;

  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta); 
	
  //Get eta position. Careful with ALICE conventions (increase index decrease eta)	
  Int_t ieta2 = ieta;
  if(nSupMod%2 == 0)		  
	  ieta2 = (fCentersOfCellsEtaDir.GetSize()-1)-ieta;// 47-ieta, revert the ordering on A side in order to keep convention.
  zr = fCentersOfCellsEtaDir.At(ieta2); 
  xr = fCentersOfCellsXDir.At(ieta2);

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

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::RelPosCellInSModule(Int_t absId, Double_t loc[3]) const
{
  // Look to see what the relative
  // position inside a given cell is
  // for a recpoint.	// Alice numbering scheme - Jun 03, 2006
  loc[0] = loc[1] = loc[2]=0.0;
  if(RelPosCellInSModule(absId, loc[0],loc[1],loc[2])) {
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::RelPosCellInSModule(Int_t absId, TVector3 &vloc) const
{
  // Look to see what the relative
  // position inside a given cell is
  // for a recpoint.  
  // Alice numbering scheme - Jun 03, 2006
  static Double_t loc[3];
  if(RelPosCellInSModule(absId,loc)) {
    vloc.SetXYZ(loc[0], loc[1], loc[2]);
    return kTRUE;
  } else {
    vloc.SetXYZ(0,0,0);
    return kFALSE;
  }
}

//________________________________________________________________________________________________
void AliEMCALGeoUtils::CreateListOfTrd1Modules()
{
  // Generate the list of Trd1 modules
  // which will make up the EMCAL
  // geometry

  AliDebug(2,Form(" AliEMCALGeometry::CreateListOfTrd1Modules() started "));

  AliEMCALShishKebabTrd1Module *mod=0, *mTmp=0; // current module
  if(fShishKebabTrd1Modules == 0) {
    fShishKebabTrd1Modules = new TList;
    fShishKebabTrd1Modules->SetName("ListOfTRD1");
    for(int iz=0; iz< fEMCGeometry->GetNZ(); iz++) {
      if(iz==0) {
	//        mod  = new AliEMCALShishKebabTrd1Module(TMath::Pi()/2.,this);
        mod  = new AliEMCALShishKebabTrd1Module(TMath::Pi()/2.,fEMCGeometry);
      } else {
        mTmp  = new AliEMCALShishKebabTrd1Module(*mod);
        mod   = mTmp;
      }
      fShishKebabTrd1Modules->Add(mod);
    }
  } else {
    AliDebug(2,Form(" Already exits : "));
  }
  mod = (AliEMCALShishKebabTrd1Module*)fShishKebabTrd1Modules->At(fShishKebabTrd1Modules->GetSize()-1);
  fEtaMaxOfTRD1 = mod->GetMaxEtaOfModule(0);

  AliDebug(2,Form(" fShishKebabTrd1Modules has %i modules : max eta %5.4f \n",
                  fShishKebabTrd1Modules->GetSize(),fEtaMaxOfTRD1));
  // Feb 20,2006;
  // Jun 01, 2006 - ALICE numbering scheme
  // define grid for cells in eta(z) and x directions in local coordinates system of SM
  // Works just for 2x2 case only -- ?? start here
  //
  //
  // Define grid for cells in phi(y) direction in local coordinates system of SM
  // as for 2X2 as for 3X3 - Nov 8,2006
  //
  AliDebug(2,Form(" Cells grid in phi directions : size %i\n", fCentersOfCellsPhiDir.GetSize()));
  Int_t ind=0; // this is phi index
  Int_t ieta=0, nModule=0, iphiTemp;
  Double_t xr=0., zr=0., theta=0., phi=0., eta=0., r=0., x=0.,y=0.;
  TVector3 vglob;
  Double_t ytCenterModule=0.0, ytCenterCell=0.0;

  fCentersOfCellsPhiDir.Set(fNPhi*fNPHIdiv);
  fPhiCentersOfCells.Set(fNPhi*fNPHIdiv);

  Double_t r0 = fIPDistance + fLongModuleSize/2.;
  for(Int_t it=0; it<fNPhi; it++) { // cycle on modules
    ytCenterModule = -fParSM[1] + fPhiModuleSize*(2*it+1)/2;  // center of module
    for(Int_t ic=0; ic<fNPHIdiv; ic++) { // cycle on cells in module
      if(fNPHIdiv==2) {
        ytCenterCell = ytCenterModule + fPhiTileSize *(2*ic-1)/2.;
      } else if(fNPHIdiv==3){
        ytCenterCell = ytCenterModule + fPhiTileSize *(ic-1);
      } else if(fNPHIdiv==1){
        ytCenterCell = ytCenterModule;
      }
      fCentersOfCellsPhiDir.AddAt(ytCenterCell,ind);
      // Define grid on phi direction
      // Grid is not the same for different eta bin;
      // Effect is small but is still here
      phi = TMath::ATan2(ytCenterCell, r0);
      fPhiCentersOfCells.AddAt(phi, ind);

      AliDebug(2,Form(" ind %2.2i : y %8.3f ", ind, fCentersOfCellsPhiDir.At(ind)));
      ind++;
    }
  }

  fCentersOfCellsEtaDir.Set(fNZ *fNETAdiv);
  fCentersOfCellsXDir.Set(fNZ *fNETAdiv);
  fEtaCentersOfCells.Set(fNZ *fNETAdiv * fNPhi*fNPHIdiv);
  AliDebug(2,Form(" Cells grid in eta directions : size %i\n", fCentersOfCellsEtaDir.GetSize()));
  for(Int_t it=0; it<fNZ; it++) {
    AliEMCALShishKebabTrd1Module *trd1 = GetShishKebabModule(it);
    nModule = fNPhi*it;
    for(Int_t ic=0; ic<fNETAdiv; ic++) {
      if(fNPHIdiv==2) {
        trd1->GetCenterOfCellInLocalCoordinateofSM(ic, xr, zr);      // case of 2X2
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta);
      } if(fNPHIdiv==3) {
        trd1->GetCenterOfCellInLocalCoordinateofSM3X3(ic, xr, zr);  // case of 3X3
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta);
      } if(fNPHIdiv==1) {
        trd1->GetCenterOfCellInLocalCoordinateofSM1X1(xr, zr);      // case of 1X1
        GetCellPhiEtaIndexInSModule(0, nModule, 0, ic, iphiTemp, ieta);
      }
      fCentersOfCellsXDir.AddAt(float(xr) - fParSM[0],ieta);
      fCentersOfCellsEtaDir.AddAt(float(zr) - fParSM[2],ieta);
      // Define grid on eta direction for each bin in phi
      for(int iphi=0; iphi<fCentersOfCellsPhiDir.GetSize(); iphi++) {
        x = xr + trd1->GetRadius();
        y = fCentersOfCellsPhiDir[iphi];
        r = TMath::Sqrt(x*x + y*y + zr*zr);
        theta = TMath::ACos(zr/r);
        eta   = AliEMCALShishKebabTrd1Module::ThetaToEta(theta);
        //        ind   = ieta*fCentersOfCellsPhiDir.GetSize() + iphi;
        ind   = iphi*fCentersOfCellsEtaDir.GetSize() + ieta;
        fEtaCentersOfCells.AddAt(eta, ind);
      }
      //printf(" ieta %i : xr + trd1->GetRadius() %f : zr %f : eta %f \n", ieta, xr + trd1->GetRadius(), zr, eta);
    }
  }
  for(Int_t i=0; i<fCentersOfCellsEtaDir.GetSize(); i++) {
    AliDebug(2,Form(" ind %2.2i : z %8.3f : x %8.3f", i+1,
                    fCentersOfCellsEtaDir.At(i),fCentersOfCellsXDir.At(i)));
  }

}


//________________________________________________________________________________________________
AliEMCALShishKebabTrd1Module* AliEMCALGeoUtils::GetShishKebabModule(Int_t neta) const
{
  //This method was too long to be
  //included in the header file - the
  //rule checker complained about it's
  //length, so we move it here.  It returns the
  //shishkebabmodule at a given eta index point.

  static AliEMCALShishKebabTrd1Module* trd1=0;
  if(fShishKebabTrd1Modules && neta>=0 && neta<fShishKebabTrd1Modules->GetSize()) {
    trd1 = (AliEMCALShishKebabTrd1Module*)fShishKebabTrd1Modules->At(neta);
  } else trd1 = 0;
  return trd1;
}

//___________________________________________________________________
void AliEMCALGeoUtils::PrintGeometry()
{
  //Print information from geometry
  fEMCGeometry->PrintGeometry();

  printf(" fShishKebabTrd1Modules has %i modules : max eta %5.4f \n", 
	 fShishKebabTrd1Modules->GetSize(),fEtaMaxOfTRD1);
  
  printf("\n Cells grid in eta directions : size %i\n", fCentersOfCellsEtaDir.GetSize());
  for(Int_t i=0; i<fCentersOfCellsEtaDir.GetSize(); i++) {
    printf(" ind %2.2i : z %8.3f : x %8.3f \n", i, 
	   fCentersOfCellsEtaDir.At(i),fCentersOfCellsXDir.At(i));
    int ind=0; // Nov 21,2006
    for(Int_t iphi=0; iphi<fCentersOfCellsPhiDir.GetSize(); iphi++) {
      ind = iphi*fCentersOfCellsEtaDir.GetSize() + i;
      printf("%6.4f ", fEtaCentersOfCells[ind]);
      if((iphi+1)%12 == 0) printf("\n");
    }
    printf("\n");
    
  }

  printf("\n Cells grid in phi directions : size %i\n", fCentersOfCellsPhiDir.GetSize());
  for(Int_t i=0; i<fCentersOfCellsPhiDir.GetSize(); i++) {
    double phi=fPhiCentersOfCells.At(i);
    printf(" ind %2.2i : y %8.3f : phi %7.5f(%6.2f) \n", i, fCentersOfCellsPhiDir.At(i), 
	   phi, phi*TMath::RadToDeg());
  }

}

//____________________________________________________________________________
Bool_t  AliEMCALGeoUtils::Impact(const TParticle * particle) const 
{
  // Tells if a particle enters EMCAL
  Bool_t in=kFALSE;
  Int_t absID=0;
  TVector3 vtx(particle->Vx(),particle->Vy(),particle->Vz());
  TVector3 vimpact(0,0,0);
  ImpactOnEmcal(vtx,particle->Theta(),particle->Phi(),absID,vimpact);
  if(absID >=0) 
    in=kTRUE;
  return in;
}
//____________________________________________________________________________
void AliEMCALGeoUtils::ImpactOnEmcal(TVector3 vtx, Double_t theta, Double_t phi, 
				     Int_t & absId, TVector3 & vimpact) const
{
  // calculates the impact coordinates on EMCAL (centre of a tower/not on EMCAL surface) 
  // of a neutral particle  
  // emitted in the vertex vtx[3] with direction theta and phi in the ALICE global coordinate system

  TVector3 p(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta)) ;

  vimpact.SetXYZ(0,0,0);
  absId=-1;
  if(phi==0 || theta==0) return;

   TVector3 direction;
   Double_t factor = (fIPDistance-vtx[1])/p[1];
  direction = vtx + factor*p;

  //from particle direction -> tower hitted
  GetAbsCellIdFromEtaPhi(direction.Eta(),direction.Phi(),absId);
  
  //tower absID hitted -> tower/module plane (evaluated at the center of the tower)
  Int_t nSupMod, nModule, nIphi, nIeta;
  Double_t loc[3],loc2[3],loc3[3];
  Double_t glob[3]={},glob2[3]={},glob3[3]={};
  
  if(!RelPosCellInSModule(absId,loc)) return;
  
  //loc is cell center of tower
  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);

  //look at 2 neighbours-s cell using nIphi={0,1} and nIeta={0,1}
  Int_t nIphi2,nIeta2,absId2,absId3;
  if(nIeta==0) nIeta2=1;
  else nIeta2=0;
  absId2=GetAbsCellId(nSupMod,nModule,nIphi,nIeta2);  
  if(nIphi==0) nIphi2=1;
  else nIphi2=0;
  absId3=GetAbsCellId(nSupMod,nModule,nIphi2,nIeta);

  //2nd point on emcal cell plane
  if(!RelPosCellInSModule(absId2,loc2)) return;
    
  //3rd point on emcal cell plane
  if(!RelPosCellInSModule(absId3,loc3)) return;
    
  // Get Matrix
  const TGeoHMatrix* m = GetMatrixForSuperModule(nSupMod);
  if(m) {
    m->LocalToMaster(loc, glob);
    m->LocalToMaster(loc2, glob2);
    m->LocalToMaster(loc3, glob3);
  } else {
    AliFatal("Geo matrixes are not loaded \n") ;
  }

  //Equation of Plane from glob,glob2,glob3 (Ax+By+Cz+D=0)
  Double_t a = glob[1]*(glob2[2]-glob3[2]) + glob2[1]*(glob3[2]-glob[2]) + glob3[1]*(glob[2]-glob2[2]);
  Double_t b = glob[2]*(glob2[0]-glob3[0]) + glob2[2]*(glob3[0]-glob[0]) + glob3[2]*(glob[0]-glob2[0]);
  Double_t c = glob[0]*(glob2[1]-glob3[1]) + glob2[0]*(glob3[1]-glob[1]) + glob3[0]*(glob[1]-glob2[1]);
  Double_t d = glob[0]*(glob2[1]*glob3[2]-glob3[1]*glob2[2]) + glob2[0]*(glob3[1]*glob[2]-glob[1]*glob3[2]) + glob3[0]*(glob[1]*glob2[2]-glob2[1]*glob[2]);
  d=-d;
  
  //shift equation of plane from tower/module center to surface along vector (A,B,C) normal to tower/module plane
  Double_t dist = fLongModuleSize/2.;
  Double_t norm = TMath::Sqrt(a*a+b*b+c*c);
  Double_t glob4[3]={};
  TVector3 dir(a,b,c);
  TVector3 point(glob[0],glob[1],glob[2]); 
  if(point.Dot(dir)<0) dist*=-1;
  glob4[0]=glob[0]-dist*a/norm;
  glob4[1]=glob[1]-dist*b/norm;
  glob4[2]=glob[2]-dist*c/norm;
  d = glob4[0]*a +  glob4[1]*b +  glob4[2]*c ;
  d = -d;

  //Line determination (2 points for equation of line : vtx and direction)
  //impact between line (particle) and plane (module/tower plane)
  Double_t den = a*(vtx(0)-direction(0)) + b*(vtx(1)-direction(1)) + c*(vtx(2)-direction(2));
  if(den==0){
    printf("ImpactOnEmcal() No solution :\n");
    return;
  }
  
  Double_t length = a*vtx(0)+b*vtx(1)+c*vtx(2)+d;
  length /=den;
  
  vimpact.SetXYZ(vtx(0)+length*(direction(0)-vtx(0)),vtx(1)+length*(direction(1)-vtx(1)),vtx(2)+length*(direction(2)-vtx(2)));
  
  //shift vimpact from tower/module surface to center along vector (A,B,C) normal to tower/module plane
  vimpact.SetXYZ(vimpact(0)+dist*a/norm,vimpact(1)+dist*b/norm,vimpact(2)+dist*c/norm);
  
  return;
}

//_____________________________________________________________________________
Bool_t AliEMCALGeoUtils::IsInEMCAL(Double_t x, Double_t y, Double_t z) const {
  // Checks whether point is inside the EMCal volume, used in AliEMCALv*.cxx
  //
  // Code uses cylindrical approximation made of inner radius (for speed)
  //
  // Points behind EMCAl, i.e. R > outer radius, but eta, phi in acceptance 
  // are considered to inside

  Double_t r=sqrt(x*x+y*y);

  if ( r > fEnvelop[0] ) {
     Double_t theta;
     theta  =    TMath::ATan2(r,z);
     Double_t eta;
     if(theta == 0) 
       eta = 9999;
     else 
       eta    =   -TMath::Log(TMath::Tan(theta/2.));
     if (eta < fArm1EtaMin || eta > fArm1EtaMax)
       return 0;
 
     Double_t phi = TMath::ATan2(y,x) * 180./TMath::Pi();
     if (phi < 0) phi += 360;  // phi should go from 0 to 360 in this case
     if (phi > fArm1PhiMin && phi < fArm1PhiMax)
       return 1;
  }
  return 0;
}

//________________________________________________________________________________________________
Int_t AliEMCALGeoUtils::GetAbsTRUNumberFromNumberInSm(const Int_t row, const Int_t col, const Int_t sm) const
{ // Nov 6, 2007
  // Get TRU absolute number from column, row and Super Module number
  Int_t itru = row + col*fEMCGeometry->GetNModulesInTRUPhi() + sm*fEMCGeometry->GetNTRU();
  // printf("  GetAbsTRUNumberFromNumberInSm : row %2i col %2i sm %2i -> itru %2i\n", row, col, sm, itru); 
  return itru;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::GetAbsFastORIndexFromTRU(const Int_t iTRU, const Int_t iADC, Int_t& id) const
{
	//Trigger mapping method, get  FastOr Index from TRU

    if (iTRU > 31 || iTRU < 0 || iADC > 95 || iADC < 0) 
	{
		AliError("TRU out of range!");
		return kFALSE;
	}
				 
	id = iADC + iTRU * 96;
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::GetTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iADC) const
{

	//Trigger mapping method, get TRU number from FastOr Index

	if (id > 3071 || id < 0)
	{
		AliError("Id out of range!");
		return kFALSE;
	}
	
	iTRU = id / 96;
	iADC = id % 96;
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::GetPositionInTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi) const
{
	//Trigger mapping method, get position in TRU from FasOr Index
	
	Int_t iADC;
	
	Bool_t isOK = GetTRUFromAbsFastORIndex(id, iTRU, iADC);
	
	if (!isOK) return kFALSE;
	
	Int_t x = iADC / 4;
	Int_t y = iADC % 4;
	
	if ( int( iTRU / 3 ) % 2 ) // C side 
	{
		iEta = 23 - x;
		iPhi =      y;
	}
	else                       // A side
	{
		iEta =      x;
		iPhi =  3 - y;
	}
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::GetPositionInSMFromAbsFastORIndex(const Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi) const
{
	//Trigger mapping method, get position in Super Module from FasOr Index

	Int_t iTRU;
	Bool_t isOK = GetPositionInTRUFromAbsFastORIndex(id, iTRU, iEta, iPhi);
	
	if (!isOK) return kFALSE;
	
	iSM  = iTRU / 3;
	
	if ( int( iTRU / 3 ) % 2 ) // C side
	{
		iPhi = iPhi + 4 * ( 2 - ( iTRU % 3 ) );
	}
	else                       // A side
	{
		iPhi = iPhi + 4 * (       iTRU % 3   );
	}
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeoUtils::GetAbsFastORIndexFromPositionInTRU(const Int_t iTRU, const Int_t iEta, const Int_t iPhi, Int_t& id) const
{
	//Trigger mapping method, get Index if FastOr from Position in TRU

	if (iTRU < 0 || iTRU > 31 || iEta < 0 || iEta > 23 || iPhi < 0 || iPhi > 3) return kFALSE;
	
	if ( int( iTRU / 3 ) % 2 ) // C side
	{
		id =      iPhi  + 4 * ( 23 - iEta ) + iTRU * 96;
	}
	else 
	{
		id = (3 - iPhi) + 4 *        iEta + iTRU * 96;
	}
	
	return kTRUE;
}

//____________________________________________________________________________
const TGeoHMatrix * AliEMCALGeoUtils::GetMatrixForSuperModule(Int_t smod) const {

	//Provides shift-rotation matrix for EMCAL
	
	if(smod < 0 || smod > fEMCGeometry->GetNumberOfSuperModules()) 
		AliFatal(Form("Wrong supermodule index -> %d",smod));
		
	//If GeoManager exists, take matrixes from it
	
	//
	//    if(fKey110DEG && ind>=10) {
	//    }
	//
	//    if(!gGeoManager->cd(volpath.Data()))
	//      AliFatal(Form("AliEMCALGeometry::GeoManager cannot find path %s!",volpath.Data()));
	//
	//    TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
	
	if(gGeoManager){
		char path[255] ;
		sprintf(path,"/ALIC_1/XEN1_1/SMOD_%d",smod+1) ;
		//TString volpath = "ALIC_1/XEN1_1/SMOD_";
	    //volpath += smod+1;

		if(fKey110DEG && smod >= 10){
			  sprintf(path,"/ALIC_1/XEN1_1/SM10_%d",smod-10+1) ;
			//volpath = "ALIC_1/XEN1_1/SM10_";
			//volpath += smod-10+1;
		}
		if (!gGeoManager->cd(path)){
			AliFatal(Form("Geo manager can not find path %s!\n",path));
		}
		return gGeoManager->GetCurrentMatrix();
	}

	if(fkSModuleMatrix[smod]){
		return fkSModuleMatrix[smod] ;
	}
	else{
		AliInfo("Stop:");
		printf("\t Can not find EMCAL misalignment matrixes\n") ;
		printf("\t Either import TGeoManager from geometry.root or \n");
		printf("\t read stored matrixes from AliESD Header:  \n") ;   
		printf("\t AliEMCALGeoUtils::SetMisalMatrixes(header->GetEMCALMisalMatrix()) \n") ;
		abort() ;
	}
	return 0 ;
}

//______________________________________________________________________
void AliEMCALGeoUtils::GetModulePhiEtaIndexInSModuleFromTRUIndex(Int_t itru, Int_t iphitru, Int_t ietatru, Int_t &iphiSM, Int_t &ietaSM) const 
{
  
  // This method transforms the (eta,phi) index of module in a 
  // TRU matrix into Super Module (eta,phi) index.
  
  // Calculate in which row and column where the TRU are 
  // ordered in the SM

  Int_t col = itru/fEMCGeometry->GetNTRUPhi() ; // indexes of TRU in SM
  Int_t row = itru - col*fEMCGeometry->GetNTRUPhi();
   
  iphiSM = fEMCGeometry->GetNModulesInTRUPhi()*row + iphitru  ;
  ietaSM = fEMCGeometry->GetNModulesInTRUEta()*col + ietatru  ; 
  //printf(" GetModulePhiEtaIndexInSModuleFromTRUIndex : itru %2i iphitru %2i ietatru %2i iphiSM %2i ietaSM %2i \n", 
  // itru, iphitru, ietatru, iphiSM, ietaSM);
}
