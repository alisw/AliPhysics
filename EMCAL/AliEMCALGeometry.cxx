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
//   EMCAL_FIRSTYEARV1 - geometry for December 2009 to December 2010 run period; 
//                fixed bug for positions of modules inside SM
//                (first module has tilt 0.75 degree);
//                the sizes updated with last information from production
//                drawing (end of October 2010). 
//      
//   EMCAL_COMPLETEV1: Same fixes as FIRSTYEAR and 10 SM instead of 10+2 half SM
//
//   EMCAL_WSUC (Wayne State test stand)
//      = no definite equivalent in old notation, was only used by
//          Aleksei, but kept for testing purposes
//
//   etc.

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
//
//*-- Author: Sahal Yacoob (LBL / UCT)
//     and  : Yves Schutz (SUBATECH)
//     and  : Jennifer Klay (LBL)
//     and  : Alexei Pavlinov (WSU) 
//
//  Implementation for analysis usage, before AliEMCALGeometry now (06/2011) merged again
//  in AliEMCALGeometry
//                  
// -- Author: Magali Estienne (magali.estienne@subatech.in2p3.fr)
//
//
// Usage: 
//        You can create the AliEMCALGeometry object independently from anything.
//        You have to use just the correct name of geometry. If name is empty string the
//        default name of geometry will be used.
//         
//  AliEMCALGeometry* geom = new AliEMCALGeometry("EMCAL_COMPLETEV1","EMCAL");
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
#include <TGeoBBox.h>
#include <TList.h>
#include <TBrowser.h>

// --- Standard library ---
//#include <Riostream.h>

// --- AliRoot header files ---
#include "AliEMCALGeometry.h"
#include "AliEMCALShishKebabTrd1Module.h"

ClassImp(AliEMCALGeometry)

// these initialisations are needed for a singleton
AliEMCALGeometry  *AliEMCALGeometry::fgGeom      = 0;
const Char_t*      AliEMCALGeometry::fgkDefaultGeometryName = "EMCAL_COMPLETEV1";

//____________________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry():
  fEMCGeometry(0x0),fGeoName(0),
  fKey110DEG(0),fNCellsInSupMod(0),fNETAdiv(0),fNPHIdiv(0),
  fNCellsInModule(0),fPhiBoundariesOfSM(0x0),fPhiCentersOfSM(0x0),
  fPhiCentersOfCells(0x0),fCentersOfCellsEtaDir(0x0),
  fCentersOfCellsPhiDir(0x0),fEtaCentersOfCells(0x0),
  fNCells(0),fNPhi(0),fCentersOfCellsXDir(0x0),fArm1EtaMin(0),
  fArm1EtaMax(0),fArm1PhiMin(0),fArm1PhiMax(0),fEtaMaxOfTRD1(0),
  fShishKebabTrd1Modules(0),fPhiModuleSize(0.),
  fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fNZ(0),
  fIPDistance(0.),fLongModuleSize(0.),fShellThickness(0.),
  fZLength(0.),fSampling(0.),fUseExternalMatrices(kFALSE)
{
  // default ctor 
  // must be kept public for root persistency purposes, but should never be called by the outside world
  fEnvelop[0] = 0.;
  fEnvelop[1] = 0.;
  fEnvelop[2] = 0.;
  fParSM[0]   = 0.;
  fParSM[1]   = 0.;
  fParSM[2]   = 0.;
  for (Int_t i=0;i<AliEMCALGeoParams::fgkEMCALModules;i++)
    fkSModuleMatrix[i]=0 ;

  for (Int_t i = 0; i < 48; i++)
	for (Int_t j = 0; j < 64; j++) fFastOR2DMap[i][j] = -1;
}  

//____________________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry(const AliEMCALGeometry & geo)
  : TNamed(geo),
    fEMCGeometry(geo.fEMCGeometry),fGeoName(geo.fGeoName),
    fKey110DEG(geo.fKey110DEG),fNCellsInSupMod(geo.fNCellsInSupMod),fNETAdiv(geo.fNETAdiv),fNPHIdiv(geo.fNPHIdiv),
    fNCellsInModule(geo.fNCellsInModule),fPhiBoundariesOfSM(geo.fPhiBoundariesOfSM),fPhiCentersOfSM(geo.fPhiCentersOfSM),
    fPhiCentersOfCells(geo.fPhiCentersOfCells),fCentersOfCellsEtaDir(geo.fCentersOfCellsEtaDir),
    fCentersOfCellsPhiDir(geo.fCentersOfCellsPhiDir),fEtaCentersOfCells(geo.fEtaCentersOfCells),
    fNCells(geo.fNCells),fNPhi(geo.fNPhi),fCentersOfCellsXDir(geo.fCentersOfCellsXDir),fArm1EtaMin(geo.fArm1EtaMin),
    fArm1EtaMax(geo.fArm1EtaMax),fArm1PhiMin(geo.fArm1PhiMin),fArm1PhiMax(geo.fArm1PhiMax),fEtaMaxOfTRD1(geo.fEtaMaxOfTRD1),
    fShishKebabTrd1Modules(geo.fShishKebabTrd1Modules),fPhiModuleSize(geo.fPhiModuleSize),
    fEtaModuleSize(geo.fEtaModuleSize),fPhiTileSize(geo.fPhiTileSize),fEtaTileSize(geo.fEtaTileSize),fNZ(geo.fNZ),
    fIPDistance(geo.fIPDistance),fLongModuleSize(geo.fLongModuleSize),fShellThickness(geo.fShellThickness),
    fZLength(geo.fZLength),fSampling(geo.fSampling),fUseExternalMatrices(geo.fUseExternalMatrices)
{
  fEnvelop[0] = geo.fEnvelop[0];
  fEnvelop[1] = geo.fEnvelop[1];
  fEnvelop[2] = geo.fEnvelop[2];
  fParSM[0]   = geo.fParSM[0];
  fParSM[1]   = geo.fParSM[1];
  fParSM[2]   = geo.fParSM[2];
  for (Int_t i=0;i<AliEMCALGeoParams::fgkEMCALModules;i++)
    fkSModuleMatrix[i]=0 ;
  
  for (Int_t i = 0; i < 48; i++)
	for (Int_t j = 0; j < 64; j++) fFastOR2DMap[i][j] = geo.fFastOR2DMap[i][j];
}

//____________________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry(const Text_t* name, const Text_t* title) 
  : TNamed(name, title),
    fEMCGeometry(0x0),fGeoName(0),
    fKey110DEG(0),fNCellsInSupMod(0),fNETAdiv(0),fNPHIdiv(0),
    fNCellsInModule(0),fPhiBoundariesOfSM(0x0),fPhiCentersOfSM(0x0),
    fPhiCentersOfCells(0x0),fCentersOfCellsEtaDir(0x0),
    fCentersOfCellsPhiDir(0x0),fEtaCentersOfCells(0x0),
    fNCells(0),fNPhi(0),fCentersOfCellsXDir(0x0),fArm1EtaMin(0),
    fArm1EtaMax(0),fArm1PhiMin(0),fArm1PhiMax(0),fEtaMaxOfTRD1(0),
    fShishKebabTrd1Modules(0),fPhiModuleSize(0.),
    fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fNZ(0),
    fIPDistance(0.),fLongModuleSize(0.),fShellThickness(0.),
    fZLength(0.),fSampling(0.), fUseExternalMatrices(kFALSE)
{ 

  // ctor only for normal usage 

  fEMCGeometry = new AliEMCALEMCGeometry(name,title);

  fGeoName = fEMCGeometry->GetGeoName();
  fKey110DEG = fEMCGeometry->GetKey110DEG();
  fNCellsInSupMod = fEMCGeometry->GetNCellsInSupMod();
  fNETAdiv = fEMCGeometry->GetNETAdiv();
  fNPHIdiv = fEMCGeometry->GetNPHIdiv();
  fNCellsInModule = fNPHIdiv*fNETAdiv;
  static int i=0;
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
  fParSM[0]   = fEMCGeometry->GetSuperModulesPar(0);
  fParSM[1]   = fEMCGeometry->GetSuperModulesPar(1);
  fParSM[2]   = fEMCGeometry->GetSuperModulesPar(2);
  fArm1EtaMin = fEMCGeometry->GetArm1EtaMin();
  fArm1EtaMax = fEMCGeometry->GetArm1EtaMax();
  fArm1PhiMin = fEMCGeometry->GetArm1PhiMin();
  fArm1PhiMax = fEMCGeometry->GetArm1PhiMax();
  fShellThickness = fEMCGeometry->GetShellThickness();
  fZLength    = fEMCGeometry->GetZLength();
  fSampling   = fEMCGeometry->GetSampling();
  fEtaModuleSize = fEMCGeometry->GetEtaModuleSize();
  fPhiModuleSize = fEMCGeometry->GetPhiModuleSize();
  fEtaTileSize = fEMCGeometry->GetEtaTileSize();
  fPhiTileSize = fEMCGeometry->GetPhiTileSize();
  fNZ          = fEMCGeometry->GetNZ();
  fIPDistance  = fEMCGeometry->GetIPDistance();
  fLongModuleSize = fEMCGeometry->GetLongModuleSize();

  CreateListOfTrd1Modules();

  for(Int_t smod=0; smod < AliEMCALGeoParams::fgkEMCALModules; smod++)
    fkSModuleMatrix[smod]=0 ;	
	
  if (AliDebugLevel()>=2) {
    fEMCGeometry->Print();
    PrintGeometryGeoUtils();
  }

  for (Int_t ix = 0; ix < 48; ix++)
	for (Int_t jx = 0; jx < 64; jx++) fFastOR2DMap[ix][jx] = -1;

  BuildFastOR2DMap();
}

//____________________________________________________________________________
AliEMCALGeometry & AliEMCALGeometry::operator = (const AliEMCALGeometry  & /*rvalue*/) { 
  //assing operator
  Fatal("assignment operator", "not implemented") ; 
  return *this ;
}

//____________________________________________________________________________
AliEMCALGeometry::~AliEMCALGeometry(void)
{
  // dtor
  if (this==fgGeom) {
    AliError("Do not call delete on me");
    return;
  }
  if (fEMCGeometry){ 
    for(Int_t smod = 0 ; smod < fEMCGeometry->GetNumberOfSuperModules(); smod++){
      if(fkSModuleMatrix[smod])
        delete fkSModuleMatrix[smod] ;
        fkSModuleMatrix[smod]=0 ;
    }
    delete fEMCGeometry; fEMCGeometry = 0 ;
  }
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
void AliEMCALGeometry::Browse(TBrowser* b)
{
  //Browse the modules
  if(fShishKebabTrd1Modules) b->Add(fShishKebabTrd1Modules);
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::IsFolder() const
{
  //Check if fShishKebabTrd1Modules is in folder
  if(fShishKebabTrd1Modules) return kTRUE;
  else                       return kFALSE;
}

//________________________________________________________________________________________________
void AliEMCALGeometry::GetGlobal(const Double_t *loc, Double_t *glob, int ind) const
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
void AliEMCALGeometry::GetGlobal(const TVector3 &vloc, TVector3 &vglob, int ind) const
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
void AliEMCALGeometry::GetGlobal(Int_t absId , double glob[3]) const
{
  // Alice numbering scheme - Jun 03, 2006
  static Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1;
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
void AliEMCALGeometry::GetGlobal(Int_t absId , TVector3 &vglob) const
{
  // Alice numbering scheme - Jun 03, 2006
  static Double_t glob[3];

  GetGlobal(absId, glob);
  vglob.SetXYZ(glob[0], glob[1], glob[2]);

}


//______________________________________________________________________
void AliEMCALGeometry::PrintCellIndexes(Int_t absId, int pri, const char *tit) const
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

void AliEMCALGeometry::PrintLocalTrd1(Int_t pri) const
{
  // For comparing with numbers from drawing
  for(Int_t i=0; i<GetShishKebabTrd1Modules()->GetSize(); i++){
    printf(" %s | ", GetShishKebabModule(i)->GetName());
    if(i==0 && pri<1) GetShishKebabModule(i)->PrintShish(1);
    else     GetShishKebabModule(i)->PrintShish(pri);
  }
}

//________________________________________________________________________________________________
void AliEMCALGeometry::EtaPhiFromIndex(Int_t absId,Double_t &eta,Double_t &phi) const
{
  // Nov 16, 2006- float to double
  // version for TRD1 only
  static TVector3 vglob;
  GetGlobal(absId, vglob);
  eta = vglob.Eta();
  phi = vglob.Phi();
}

//________________________________________________________________________________________________
void AliEMCALGeometry::EtaPhiFromIndex(Int_t absId,Float_t &eta,Float_t &phi) const
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
Int_t AliEMCALGeometry::GetAbsCellId(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta) const
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
void  AliEMCALGeometry::GetModuleIndexesFromCellIndexesInSModule(Int_t nSupMod, Int_t iphi, Int_t ieta, 
			Int_t &iphim, Int_t &ietam, Int_t &nModule) const
{
  // Transition from cell indexes (ieta,iphi) to module indexes (ietam,iphim, nModule)
  static Int_t nphi=-1;
  nphi  = GetNumberOfModuleInPhiDirection(nSupMod);  

  ietam  = ieta/fNETAdiv;
  iphim  = iphi/fNPHIdiv;
  nModule = ietam * nphi + iphim; 
}

//________________________________________________________________________________________________
Int_t  AliEMCALGeometry::GetAbsCellIdFromCellIndexes(Int_t nSupMod, Int_t iphi, Int_t ieta) const
{
  // Transition from super module number(nSupMod) and cell indexes (ieta,iphi) to absId
  static Int_t ietam=-1, iphim=-1, nModule=-1;
  static Int_t nIeta=-1, nIphi=-1; // cell indexes in module

  GetModuleIndexesFromCellIndexesInSModule(nSupMod, iphi, ieta, ietam, iphim, nModule);

  nIeta = ieta%fNETAdiv;
  nIeta = fNETAdiv - 1 - nIeta;
  nIphi = iphi%fNPHIdiv;

  return GetAbsCellId(nSupMod, nModule, nIphi, nIeta);
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::SuperModuleNumberFromEtaPhi(Double_t eta, Double_t phi, Int_t &nSupMod) const
{ 
  // Return false if phi belongs a phi cracks between SM
 
  static Int_t i=0;

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
Bool_t AliEMCALGeometry::GetAbsCellIdFromEtaPhi(Double_t eta, Double_t phi, Int_t &absId) const
{
  // Nov 17,2006
  // stay here - phi problem as usual 
  static Int_t nSupMod=-1, i=0, ieta=-1, iphi=-1, etaShift=0, nphi=-1;
  static Double_t absEta=0.0, d=0.0, dmin=0.0, phiLoc=0;
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
Bool_t  AliEMCALGeometry::CheckAbsCellId(Int_t absId) const
{ 
  // May 31, 2006; only trd1 now
  if(absId<0 || absId >= fNCells) return kFALSE;
  else                            return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetCellIndex(Int_t absId,Int_t &nSupMod,Int_t &nModule,Int_t &nIphi,Int_t &nIeta) const
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
Int_t  AliEMCALGeometry::GetSuperModuleNumber(Int_t absId)  const
{
  // Return the number of the  supermodule given the absolute
  // ALICE numbering id

  static Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1;
  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  return nSupMod;
} 

//________________________________________________________________________________________________
void AliEMCALGeometry::GetModulePhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule,  int &iphim, int &ietam) const
{ 
  // added nSupMod; - 19-oct-05 !
  // Alice numbering scheme        - Jun 01,2006 
  // ietam, iphi - indexes of module in two dimensional grid of SM
  // ietam - have to change from 0 to fNZ-1
  // iphim - have to change from 0 to nphi-1 (fNPhi-1 or fNPhi/2-1)
  static Int_t nphi=-1;

  if(fKey110DEG == 1 && nSupMod>=10) nphi = fNPhi/2;
  else                               nphi = fNPhi;

  ietam = nModule/nphi;
  iphim = nModule%nphi;
}

//________________________________________________________________________________________________
void AliEMCALGeometry::GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta, 
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
  static Int_t iphim=-1, ietam=-1;

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
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, Double_t &xr, Double_t &yr, Double_t &zr) const
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
  static Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1;
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
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, Double_t loc[3]) const
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
Bool_t AliEMCALGeometry::RelPosCellInSModule(Int_t absId, TVector3 &vloc) const
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
  static Int_t nSupMod=0, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1;
  static Int_t iphim=-1, ietam=-1;
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


//________________________________________________________________________________________________
void AliEMCALGeometry::CreateListOfTrd1Modules()
{
  // Generate the list of Trd1 modules
  // which will make up the EMCAL
  // geometry
  // key: look to the AliEMCALShishKebabTrd1Module::

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
AliEMCALShishKebabTrd1Module* AliEMCALGeometry::GetShishKebabModule(Int_t neta) const
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
void AliEMCALGeometry::PrintGeometryGeoUtils()
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
Bool_t  AliEMCALGeometry::Impact(const TParticle * particle) const 
{
  // Tells if a particle enters EMCAL
  Bool_t in=kFALSE;
  Int_t absID=0;
  TVector3 vtx(particle->Vx(),particle->Vy(),particle->Vz());
  TVector3 vimpact(0,0,0);
  ImpactOnEmcal(vtx,particle->Theta(),particle->Phi(),absID,vimpact);
  if(absID>=0) 
    in=kTRUE;
  return in;
}
//____________________________________________________________________________
void AliEMCALGeometry::ImpactOnEmcal(TVector3 vtx, Double_t theta, Double_t phi, 
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
  Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1;
  Double_t loc[3],loc2[3],loc3[3];
  Double_t glob[3]={},glob2[3]={},glob3[3]={};
  
  if(!RelPosCellInSModule(absId,loc)) return;
  
  //loc is cell center of tower
  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);

  //look at 2 neighbours-s cell using nIphi={0,1} and nIeta={0,1}
  Int_t nIphi2=-1,nIeta2=-1,absId2=-1,absId3=-1;
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
Bool_t AliEMCALGeometry::IsInEMCAL(Double_t x, Double_t y, Double_t z) const {
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
Int_t AliEMCALGeometry::GetAbsTRUNumberFromNumberInSm(const Int_t row, const Int_t col, const Int_t sm) const
{ // Nov 6, 2007
  // Get TRU absolute number from column, row and Super Module number
  Int_t itru = row + col*fEMCGeometry->GetNModulesInTRUPhi() + sm*fEMCGeometry->GetNTRU();
  // printf("  GetAbsTRUNumberFromNumberInSm : row %2i col %2i sm %2i -> itru %2i\n", row, col, sm, itru); 
  return itru;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetAbsFastORIndexFromTRU(const Int_t iTRU, const Int_t iADC, Int_t& id) const
{
	//Trigger mapping method, get  FastOr Index from TRU

  if (iTRU > 31 || iTRU < 0 || iADC > 95 || iADC < 0) 
	{
		AliError("TRU out of range!");
		return kFALSE;
	}
	
	id  = ( iTRU % 2 ) ? iADC%4 + 4 * (23 - int(iADC/4)) : (3 - iADC%4) + 4 * int(iADC/4);

	id += iTRU * 96;
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iADC) const
{

	//Trigger mapping method, get TRU number from FastOr Index

	if (id > 3071 || id < 0)
	{
		AliError("Id out of range!");
		return kFALSE;
	}
	
	iTRU = id / 96;
	
	iADC = id % 96;
	
	iADC = ( iTRU % 2 ) ? iADC%4 + 4 * (23 - int(iADC/4)) : (3 - iADC%4) + 4 * int(iADC/4);
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetPositionInTRUFromAbsFastORIndex(const Int_t id, Int_t& iTRU, Int_t& iEta, Int_t& iPhi) const
{
	//Trigger mapping method, get position in TRU from FasOr Index
	
	Int_t iADC=-1;	
	if (!GetTRUFromAbsFastORIndex(id, iTRU, iADC)) return kFALSE;
	
	Int_t x = iADC / 4;
	Int_t y = iADC % 4;
	
	if ( iTRU % 2 ) // C side 
	{
		iEta = 23 - x;
		iPhi =      y;
	}
	else            // A side
	{
		iEta =      x;
		iPhi =  3 - y;
	}
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetPositionInSMFromAbsFastORIndex(const Int_t id, Int_t& iSM, Int_t& iEta, Int_t& iPhi) const
{
	//Trigger mapping method, get position in Super Module from FasOr Index

	Int_t iTRU=-1;
		
	if (!GetPositionInTRUFromAbsFastORIndex(id, iTRU, iEta, iPhi)) return kFALSE;
	
	if (iTRU % 2) // C side
	{
		iSM  = 2 * ( int( int(iTRU / 2) / 3 ) ) + 1;
	}
	else            // A side
	{
		iSM  = 2 * ( int( int(iTRU / 2) / 3 ) );
	}

	iPhi += 4 * int((iTRU % 6) / 2);
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetPositionInEMCALFromAbsFastORIndex(const Int_t id, Int_t& iEta, Int_t& iPhi) const
{
  //Trigger mapping method, get position in EMCAL from FastOR index

	Int_t iSM=-1;
	
	if (GetPositionInSMFromAbsFastORIndex(id, iSM, iEta, iPhi))
	{
		if (iSM % 2) iEta += 24; 
		
		iPhi += 12 * int(iSM / 2);
		
		return kTRUE;
	}
	
	return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetAbsFastORIndexFromPositionInTRU(const Int_t iTRU, const Int_t iEta, const Int_t iPhi, Int_t& id) const
{
	//Trigger mapping method, get Index if FastOr from Position in TRU

	if (iTRU < 0 || iTRU > 31 || iEta < 0 || iEta > 23 || iPhi < 0 || iPhi > 3) 
	{
		AliError("Out of range!");	
		return kFALSE;
	}
	
	id =  iPhi  + 4 * iEta + iTRU * 96;
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetAbsFastORIndexFromPositionInSM(const Int_t  iSM, const Int_t iEta, const Int_t iPhi, Int_t& id) const
{
  //Trigger mapping method, from position in SM Index get FastOR index 

	if (iSM < 0 || iSM > 11 || iEta < 0 || iEta > 23 || iPhi < 0 || iPhi > 11) 
	{
		AliError("Out of range!");
		return kFALSE;
	}
	
	Int_t x = iEta;
	Int_t y = iPhi % 4;	
	
	Int_t iOff = (iSM % 2) ? 1 : 0;
	Int_t iTRU = 2 * int(iPhi / 4) + 6 * int(iSM / 2) + iOff;

	if (GetAbsFastORIndexFromPositionInTRU(iTRU, x, y, id))
	{
		return kTRUE;
	}
	
	return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetAbsFastORIndexFromPositionInEMCAL(const Int_t iEta, const Int_t iPhi, Int_t& id) const
{
  //Trigger mapping method, from position in EMCAL Index get FastOR index 

	if (iEta < 0 || iEta > 47 || iPhi < 0 || iPhi > 63 ) 
	{
		AliError("Out of range!");
		return kFALSE;
	}
	
	if (fFastOR2DMap[iEta][iPhi] == -1) 
	{
		AliError("Invalid index!");
		return kFALSE;
	}
	
	id = fFastOR2DMap[iEta][iPhi];
	
	return kTRUE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetFastORIndexFromCellIndex(const Int_t id, Int_t& idx) const
{
  //Trigger mapping method, from cell index get FastOR index 

	Int_t iSupMod, nModule, nIphi, nIeta, iphim, ietam;
	
	Bool_t isOK = GetCellIndex( id, iSupMod, nModule, nIphi, nIeta );
	
	GetModulePhiEtaIndexInSModule( iSupMod, nModule, iphim, ietam );
	
	if (isOK && GetAbsFastORIndexFromPositionInSM(iSupMod, ietam, iphim, idx))
	{
		return kTRUE;
	}
	
	return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetCellIndexFromFastORIndex(const Int_t id, Int_t idx[4]) const
{
  //Trigger mapping method, from FASTOR index get cell index 

  Int_t iSM=-1, iEta=-1, iPhi=-1;
	if (GetPositionInSMFromAbsFastORIndex(id, iSM, iEta, iPhi))
	{
		Int_t ix = 2 * iEta;
		Int_t iy = 2 * iPhi;
		
		for (Int_t i=0; i<2; i++)
		{
			for (Int_t j=0; j<2; j++)
			{
				idx[2*i+j] = GetAbsCellIdFromCellIndexes(iSM, iy + i, ix + j);
			}
		}
		
		return kTRUE;
	}
	
	return kFALSE;
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetTRUIndexFromSTUIndex(const Int_t id, Int_t& idx) const
{
  //Trigger mapping method, from STU index get TRU index 

	if (id > 31 || id < 0) 
	{
		AliError(Form("TRU index out of range: %d",id));
		return kFALSE;
	}
	
	idx = (id > 15) ? 2 * (31 - id) : 2 * (15 - id) + 1;
	
	return kTRUE;
}

//________________________________________________________________________________________________
Int_t AliEMCALGeometry::GetTRUIndexFromSTUIndex(const Int_t id) const
{
  //Trigger mapping method, from STU index get TRU index 

	if (id > 31 || id < 0) 
	{
		AliError(Form("TRU index out of range: %d",id));
	}
	
	Int_t idx = (id > 15) ? 2 * (31 - id) : 2 * (15 - id) + 1;
	
	return idx;
}

//________________________________________________________________________________________________
void AliEMCALGeometry::BuildFastOR2DMap()
{
	// Needed by STU
	for (Int_t i = 0; i < 32; i++)
	{
		for (Int_t j = 0; j < 24; j++)
		{
			for (Int_t k = 0; k < 4; k++)
			{
				Int_t id;
				if (GetAbsFastORIndexFromPositionInTRU(i, j, k, id))
				{
					Int_t x = j, y = k + 4 * int(i / 2);
				
					if (i % 2) x += 24;
				
					fFastOR2DMap[x][y] = id;
				}
			}			
		}
	}
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetFastORIndexFromL0Index(const Int_t iTRU, const Int_t id, Int_t idx[], const Int_t size) const
{
  //Trigger mapping method, from L0 index get FastOR index 
	if (size <= 0 ||size > 4)
	{
		AliError("Size not supported!");
		return kFALSE;
	}
		
	Int_t motif[4] = {0, 1, 4, 5};
	
	switch (size)
	{
		case 1: // Cosmic trigger
			if (!GetAbsFastORIndexFromTRU(iTRU, id, idx[1])) return kFALSE;
			break;
		case 4: // 4 x 4
			for (Int_t k = 0; k < 4; k++)
			{
				Int_t iADC = motif[k] + 4 * int(id / 3) + (id % 3);
				
				if (!GetAbsFastORIndexFromTRU(iTRU, iADC, idx[k])) return kFALSE;
			}
			break;
		default:
			break;
	}
	
	return kTRUE;
}

//____________________________________________________________________________
const TGeoHMatrix * AliEMCALGeometry::GetMatrixForSuperModule(Int_t smod) const {

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
  
  //Use matrices set externally
	if(!gGeoManager || (gGeoManager && fUseExternalMatrices)){
    if(fkSModuleMatrix[smod]){
      return fkSModuleMatrix[smod] ;
    }
    else{
      AliInfo("Stop:");
      printf("\t Can not find EMCAL misalignment matrixes\n") ;
      printf("\t Either import TGeoManager from geometry.root or \n");
      printf("\t read stored matrixes from AliESD Header:  \n") ;   
      printf("\t AliEMCALGeometry::SetMisalMatrixes(header->GetEMCALMisalMatrix()) \n") ;
      abort() ;
    }  
  }//external matrices
  
	if(gGeoManager){
    const Int_t buffersize = 255;
		char path[buffersize] ;
		snprintf(path,buffersize,"/ALIC_1/XEN1_1/SMOD_%d",smod+1) ;
		//TString volpath = "ALIC_1/XEN1_1/SMOD_";
	    //volpath += smod+1;

		if(fKey110DEG && smod >= 10){
			  snprintf(path,buffersize,"/ALIC_1/XEN1_1/SM10_%d",smod-10+1) ;
			//volpath = "ALIC_1/XEN1_1/SM10_";
			//volpath += smod-10+1;
		}
		if (!gGeoManager->cd(path)){
			AliFatal(Form("Geo manager can not find path %s!\n",path));
		}
		return gGeoManager->GetCurrentMatrix();
	}

	return 0 ;
}

//______________________________________________________________________
void AliEMCALGeometry::GetModulePhiEtaIndexInSModuleFromTRUIndex(Int_t itru, Int_t iphitru, Int_t ietatru, Int_t &iphiSM, Int_t &ietaSM) const 
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

//__________________________________________________________________________________________________________________
void AliEMCALGeometry::RecalculateTowerPosition(Float_t drow, Float_t dcol, const Int_t sm, const Float_t depth,
                                                const Float_t misaligTransShifts[15], const Float_t misaligRotShifts[15], Float_t global[3]) const
{ //Transform clusters cell position into global with alternative method, taking into account the depth calculation.
  //Input are: the tower indeces, 
  //           supermodule, 
  //           particle type (photon 0, electron 1, hadron 2 )
  //           misalignment shifts to global position in case of need.
  // Federico.Ronchetti@cern.ch
  
    
  // To use in a print later
  Float_t droworg = drow;
  Float_t dcolorg = dcol;
  
  if(gGeoManager){
    //Recover some stuff

    const Int_t nSMod = fEMCGeometry->GetNumberOfSuperModules();
 
    gGeoManager->cd("ALIC_1/XEN1_1");
    TGeoNode        *geoXEn1 = gGeoManager->GetCurrentNode();
    TGeoNodeMatrix  *geoSM[nSMod];        
    TGeoVolume      *geoSMVol[nSMod];     
    TGeoShape       *geoSMShape[nSMod];    
    TGeoBBox        *geoBox[nSMod];        
    TGeoMatrix      *geoSMMatrix[nSMod];       
    
    for(int iSM = 0; iSM < nSMod; iSM++) {  
      geoSM[iSM]       = dynamic_cast<TGeoNodeMatrix *>(geoXEn1->GetDaughter(iSM));
      geoSMVol[iSM]    = geoSM[iSM]->GetVolume(); 
      geoSMShape[iSM]  = geoSMVol[iSM]->GetShape();
      geoBox[iSM]      = dynamic_cast<TGeoBBox *>(geoSMShape[iSM]);
      geoSMMatrix[iSM] = geoSM[iSM]->GetMatrix();
    }
    
    if(sm % 2 == 0) {
      dcol = 47. - dcol;
      drow = 23. - drow;
    }
    
    Int_t istrip = 0;
    Float_t z0   = 0;
    Float_t zb   = 0;
    Float_t z_is = 0;
    
    Float_t x,y,z; // return variables in terry's RF
    
    //***********************************************************
    //Do not like this: too many hardcoded values, is it not already stored somewhere else?
    //                : need more comments in the code 
    //***********************************************************
    
    Float_t dz = 6.0;   // base cell width in eta
    Float_t dx = 6.004; // base cell width in phi
    
    
    //Float_t L = 26.04; // active tower length for hadron (lead+scint+paper)
    // we use the geant numbers 13.87*2=27.74
    Float_t teta1 = 0.;
      
    //Do some basic checks
    if (dcol >= 47.5 || dcol<-0.5) {
      AliError(Form("Bad tower coordinate dcol=%f, where dcol >= 47.5 || dcol<-0.5; org: %f", dcol, dcolorg));
      return;
    }
    if (drow >= 23.5 || drow<-0.5) {
      AliError(Form("Bad tower coordinate drow=%f, where drow >= 23.5 || drow<-0.5; org: %f", drow, droworg));
      return;
    }
    if (sm >= nSMod || sm < 0) {
      AliError(Form("Bad SM number sm=%d, where sm >= %d || sm < 0", nSMod, sm));
      return;
    }    
    
    istrip = int ((dcol+0.5)/2);
    
    // tapering angle
    teta1 = TMath::DegToRad() * istrip * 1.5;
    
    // calculation of module corner along z 
    // as a function of strip
    
    for (int is=0; is<= istrip; is++) {
      
      teta1 = TMath::DegToRad() * (is*1.5 + 0.75);
      if(is==0)
        z_is = z_is + 2*dz*TMath::Cos(teta1);
      else
        z_is = z_is + 2*dz*TMath::Cos(teta1) + 2*dz*TMath::Sin(teta1)*TMath::Tan(teta1-0.75*TMath::DegToRad());
      
    }
    
    z0 = dz*(dcol-2*istrip+0.5);
    zb = (2*dz-z0-depth*TMath::Tan(teta1));
    
    z = z_is - zb*TMath::Cos(teta1);
    y = depth/TMath::Cos(teta1) + zb*TMath::Sin(teta1);
    
    x = (drow + 0.5)*dx;
    
    // moving the origin from terry's RF
    // to the GEANT one
    
    double xx =  y - geoBox[sm]->GetDX();
    double yy = -x + geoBox[sm]->GetDY(); 
    double zz =  z - geoBox[sm]->GetDZ(); 
    const double localIn[3] = {xx, yy, zz};
    double dglobal[3];
    //geoSMMatrix[sm]->Print();
    //printf("TFF Local    (row = %d, col = %d, x = %3.2f,  y = %3.2f, z = %3.2f)\n", iroworg, icolorg, localIn[0], localIn[1], localIn[2]);
    geoSMMatrix[sm]->LocalToMaster(localIn, dglobal);
    //printf("TFF Global   (row = %2.0f, col = %2.0f, x = %3.2f,  y = %3.2f, z = %3.2f)\n", drow, dcol, dglobal[0], dglobal[1], dglobal[2]);
    
    //apply global shifts
    if(sm == 2 || sm == 3) {//sector 1
      global[0] = dglobal[0] + misaligTransShifts[3] + misaligRotShifts[3]*TMath::Sin(TMath::DegToRad()*20) ; 
      global[1] = dglobal[1] + misaligTransShifts[4] + misaligRotShifts[4]*TMath::Cos(TMath::DegToRad()*20) ; 
      global[2] = dglobal[2] + misaligTransShifts[5];
    }
    else if(sm == 0 || sm == 1){//sector 0
      global[0] = dglobal[0] + misaligTransShifts[0]; 
      global[1] = dglobal[1] + misaligTransShifts[1]; 
      global[2] = dglobal[2] + misaligTransShifts[2];
    }
    else {
      AliInfo("Careful, correction not implemented yet!");
      global[0] = dglobal[0] ;
      global[1] = dglobal[1] ;
      global[2] = dglobal[2] ;
    }
    
    
  }
  else{
    AliFatal("Geometry boxes information, check that geometry.root is loaded\n");
  }
  
}
