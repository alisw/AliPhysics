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
#include <TParticle.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TList.h>
#include <TBrowser.h>

// --- Standard library ---
//#include <Riostream.h>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALShishKebabTrd1Module.h"
#include "AliEMCALTriggerMappingV1.h"
#include "AliEMCALTriggerMappingV2.h"

ClassImp(AliEMCALGeometry)

// these initialisations are needed for a singleton
AliEMCALGeometry  *AliEMCALGeometry::fgGeom      = 0;
const Char_t*      AliEMCALGeometry::fgkDefaultGeometryName = "EMCAL_COMPLETE12SMV1_DCAL_8SM";

//____________________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry():
  fEMCGeometry(0x0),fTriggerMapping(0x0),fGeoName(0),fEMCSMSystem(0x0),
  fKey110DEG(0),fnSupModInDCAL(0),fNCellsInSupMod(0),fNETAdiv(0),fNPHIdiv(0),
  fNCellsInModule(0),fPhiBoundariesOfSM(0x0),fPhiCentersOfSM(0x0),
  fPhiCentersOfSMSec(0x0),fPhiCentersOfCells(0x0),fCentersOfCellsEtaDir(0x0),
  fCentersOfCellsPhiDir(0x0),fEtaCentersOfCells(0x0),
  fNCells(0),fNPhi(0),fCentersOfCellsXDir(0x0),fArm1EtaMin(0),
  fArm1EtaMax(0),fArm1PhiMin(0),fArm1PhiMax(0),fEtaMaxOfTRD1(0),
  fDCALPhiMin(0),fDCALPhiMax(0),fEMCALPhiMax(0),fDCALStandardPhiMax(0),
  fDCALInnerExtandedEta(0),fShishKebabTrd1Modules(0),fPhiModuleSize(0.),
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
}  

//____________________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry(const AliEMCALGeometry & geo)
  : TNamed(geo),
    fEMCGeometry(geo.fEMCGeometry),fTriggerMapping(geo.fTriggerMapping),fGeoName(geo.fGeoName),fEMCSMSystem(geo.fEMCSMSystem),
    fKey110DEG(geo.fKey110DEG),fnSupModInDCAL(geo.fnSupModInDCAL),fNCellsInSupMod(geo.fNCellsInSupMod),fNETAdiv(geo.fNETAdiv),fNPHIdiv(geo.fNPHIdiv),
    fNCellsInModule(geo.fNCellsInModule),fPhiBoundariesOfSM(geo.fPhiBoundariesOfSM),fPhiCentersOfSM(geo.fPhiCentersOfSM),
    fPhiCentersOfSMSec(geo.fPhiCentersOfSMSec),fPhiCentersOfCells(geo.fPhiCentersOfCells),fCentersOfCellsEtaDir(geo.fCentersOfCellsEtaDir),
    fCentersOfCellsPhiDir(geo.fCentersOfCellsPhiDir),fEtaCentersOfCells(geo.fEtaCentersOfCells),
    fNCells(geo.fNCells),fNPhi(geo.fNPhi),fCentersOfCellsXDir(geo.fCentersOfCellsXDir),fArm1EtaMin(geo.fArm1EtaMin),
    fArm1EtaMax(geo.fArm1EtaMax),fArm1PhiMin(geo.fArm1PhiMin),fArm1PhiMax(geo.fArm1PhiMax),fEtaMaxOfTRD1(geo.fEtaMaxOfTRD1),
    fDCALPhiMin(geo.fDCALPhiMin),fDCALPhiMax(geo.fDCALPhiMax),fEMCALPhiMax(geo.fEMCALPhiMax),fDCALStandardPhiMax(geo.fDCALStandardPhiMax),
    fDCALInnerExtandedEta(geo.fDCALInnerExtandedEta),fShishKebabTrd1Modules(geo.fShishKebabTrd1Modules),fPhiModuleSize(geo.fPhiModuleSize),
    fEtaModuleSize(geo.fEtaModuleSize),fPhiTileSize(geo.fPhiTileSize),fEtaTileSize(geo.fEtaTileSize),fNZ(geo.fNZ),
    fIPDistance(geo.fIPDistance),fLongModuleSize(geo.fLongModuleSize),fShellThickness(geo.fShellThickness),
    fZLength(geo.fZLength),fSampling(geo.fSampling),fUseExternalMatrices(geo.fUseExternalMatrices)
{
  // Copy constarctor
  fEnvelop[0] = geo.fEnvelop[0];
  fEnvelop[1] = geo.fEnvelop[1];
  fEnvelop[2] = geo.fEnvelop[2];
  fParSM[0]   = geo.fParSM[0];
  fParSM[1]   = geo.fParSM[1];
  fParSM[2]   = geo.fParSM[2];
  for (Int_t i=0;i<AliEMCALGeoParams::fgkEMCALModules;i++)
    fkSModuleMatrix[i]=0 ;
}

//____________________________________________________________________________
AliEMCALGeometry::AliEMCALGeometry(const Text_t* name,   const Text_t* title,
                                   const Text_t* mcname, const Text_t* mctitle) 
  : TNamed(name, title),
    fEMCGeometry(0x0),fTriggerMapping(0x0),fGeoName(0),fEMCSMSystem(0x0),
    fKey110DEG(0),fnSupModInDCAL(0),fNCellsInSupMod(0),fNETAdiv(0),fNPHIdiv(0),
    fNCellsInModule(0),fPhiBoundariesOfSM(0x0),fPhiCentersOfSM(0x0),
    fPhiCentersOfSMSec(0x0),fPhiCentersOfCells(0x0),fCentersOfCellsEtaDir(0x0),
    fCentersOfCellsPhiDir(0x0),fEtaCentersOfCells(0x0),
    fNCells(0),fNPhi(0),fCentersOfCellsXDir(0x0),fArm1EtaMin(0),
    fArm1EtaMax(0),fArm1PhiMin(0),fArm1PhiMax(0),fEtaMaxOfTRD1(0),
    fDCALPhiMin(0),fDCALPhiMax(0),fEMCALPhiMax(0),fDCALStandardPhiMax(0),
    fDCALInnerExtandedEta(0),fShishKebabTrd1Modules(0),fPhiModuleSize(0.),
    fEtaModuleSize(0.),fPhiTileSize(0.),fEtaTileSize(0.),fNZ(0),
    fIPDistance(0.),fLongModuleSize(0.),fShellThickness(0.),
    fZLength(0.),fSampling(0.), fUseExternalMatrices(kFALSE)
{ 
  // ctor only for normal usage 
  
  fEMCGeometry = new AliEMCALEMCGeometry(name,title,mcname,mctitle);
  fGeoName = fEMCGeometry->GetGeoName();
  fEMCSMSystem = fEMCGeometry->GetEMCSystem();
  fKey110DEG = fEMCGeometry->GetKey110DEG();
  fnSupModInDCAL = fEMCGeometry->GetnSupModInDCAL();
  fNCellsInSupMod = fEMCGeometry->GetNCellsInSupMod();
  fNETAdiv = fEMCGeometry->GetNETAdiv();
  fNPHIdiv = fEMCGeometry->GetNPHIdiv();
  fNCellsInModule = fNPHIdiv*fNETAdiv;
  static int i=0;
  Int_t nSMod = fEMCGeometry->GetNumberOfSuperModules();
  fPhiBoundariesOfSM.Set(nSMod);
  fPhiCentersOfSM.Set(nSMod/2);
  fPhiCentersOfSMSec.Set(nSMod/2);
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
    fPhiCentersOfSMSec[i] = fEMCGeometry->GetPhiCenterOfSMSec(sm);
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
  fDCALPhiMin = fEMCGeometry->GetDCALPhiMin();
  fDCALPhiMax = fEMCGeometry->GetDCALPhiMax();
  fEMCALPhiMax = fEMCGeometry->GetEMCALPhiMax();
  fDCALStandardPhiMax = fEMCGeometry->GetDCALStandardPhiMax();
  fDCALInnerExtandedEta = fEMCGeometry->GetDCALInnerExtandedEta();
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
  
  AliLog::Message(AliLog::kInfo, Form("Name <<%s>>",name),
                  MODULENAME(), "AliEMCALGeometry", FUNCTIONNAME(), __FILE__, __LINE__);  
  
  if ((fEMCGeometry->GetGeoName()).Contains("DCAL")) {
    fTriggerMapping = new AliEMCALTriggerMappingV2(52, this);
    AliLog::Message(AliLog::kInfo, "EMCAL Trigger Mapping Version V2 Enabled",
                    MODULENAME(), "AliEMCALGeometry", FUNCTIONNAME(), __FILE__, __LINE__);

  } else { 
    fTriggerMapping = new AliEMCALTriggerMappingV1(32, this);
    AliLog::Message(AliLog::kInfo, "EMCAL Trigger Mapping Version V1 Enabled",
                    MODULENAME(), "AliEMCALGeometry", FUNCTIONNAME(), __FILE__, __LINE__);
  }
}

//____________________________________________________________________________
AliEMCALGeometry & AliEMCALGeometry::operator = (const AliEMCALGeometry  & /*rvalue*/) 
{ 
  //assing operator
  Fatal("assignment operator", "not implemented") ; 
  return *this ;
}

//____________________________________________________________________________
AliEMCALGeometry::~AliEMCALGeometry(void)
{
  // dtor
  if (this==fgGeom)
  {
    AliError("Do not call delete on me");
    return;
  }
  
  if (fEMCGeometry)
  {
    for(Int_t smod = 0 ; smod < fEMCGeometry->GetNumberOfSuperModules(); smod++)
    {
      if(fkSModuleMatrix[smod])
        delete fkSModuleMatrix[smod] ;
      
      fkSModuleMatrix[smod]=0 ;
    }
    
    delete fEMCGeometry; // fEMCGeometry = 0 ;
  }
  
  if (fTriggerMapping) delete fTriggerMapping;
}

//______________________________________________________________________
AliEMCALGeometry *  AliEMCALGeometry::GetInstance()
{ 
  // Returns the pointer of the unique instance
  
  AliEMCALGeometry * rv = static_cast<AliEMCALGeometry *>( fgGeom );
  return rv; 
}

//
/// \return the pointer of the unique instance
///
//______________________________________________________________________
AliEMCALGeometry* AliEMCALGeometry::GetInstance(const Text_t* name,   const Text_t* title,
                                                const Text_t* mcname, const Text_t* mctitle )
{
  AliEMCALGeometry * rv = 0; 
  
  if ( fgGeom == 0 ) 
  {
    if ( strcmp(name,"") == 0 ) 
    { // get default geometry
      fgGeom = new AliEMCALGeometry(fgkDefaultGeometryName, title, mcname, mctitle);
    } 
    else 
    {
      fgGeom = new AliEMCALGeometry(name, title,mcname,mctitle);
    }  // end if strcmp(name,"")
    
    if ( AliEMCALEMCGeometry::fgInit ) 
    {
      rv = (AliEMCALGeometry * ) fgGeom;
    }
    else 
    {
      rv = 0; 
      delete fgGeom; 
      fgGeom = 0; 
    } // end if fgInit
  }
  else
  {
    if ( strcmp(fgGeom->GetName(), name) != 0)
    {
      printf("\n current geometry is %s : ", fgGeom->GetName());
      printf(" you should not call %s \n",name);
    } // end 
    
    rv = (AliEMCALGeometry *) fgGeom; 
  }  // end if fgGeom
  
  return rv; 
}

///
/// Instanciate geometry depending on the run number
///
/// \return the pointer of the unique instance
///
//___________________________________________________________________________
AliEMCALGeometry* AliEMCALGeometry::GetInstanceFromRunNumber(Int_t runNumber,       TString geoName,
                                                             const Text_t* mcname,  const Text_t* mctitle )
{
  //printf("AliEMCALGeometry::GetInstanceFromRunNumber() - run %d, geoName <<%s>> \n",runNumber,geoName.Data());
  
  static bool showInfo = !(getenv("HLT_ONLINE_MODE") && strcmp(getenv("HLT_ONLINE_MODE"), "on") == 0);
  
  if      ( runNumber >= 104064 && runNumber < 140000  ) 
  {
    // 2009-2010 runs
    // First year geometry, 4 SM.
    
    if (showInfo)
    {
      if(!geoName.Contains("FIRSTYEARV1") && geoName!="")
      { 
        printf("AliEMCALGeometry::GetInstanceFromRunNumber() *** ATTENTION *** \n");
        printf("\t Specified geometry name <<%s>> for run %d is not considered! \n",geoName.Data(),runNumber);
        printf("\t In use <<EMCAL_FIRSTYEARV1>>, check run number and year\n");
      }
      else 
      {
        printf("AliEMCALGeometry::GetInstanceFromRunNumber() - Initialized geometry with name <<EMCAL_FIRSTYEARV1>>\n");
      }
    }
    
    return AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1","EMCAL",mcname,mctitle) ;
  }
  else if ( runNumber >= 140000 && runNumber <= 170593 )
  {
    // Almost complete EMCAL geometry, 10 SM. Year 2011 configuration
    
    if (showInfo)
    {
      if(!geoName.Contains("COMPLETEV1") && geoName!="")
      {
        printf("AliEMCALGeometry::GetInstanceFromRunNumber() *** ATTENTION *** \n");
        printf("\t Specified geometry name <<%s>> for run %d is not considered! \n",geoName.Data(),runNumber);
        printf("\t In use <<EMCAL_COMPLETEV1>>, check run number and year\n");
      }
      else 
      {
        printf("AliEMCALGeometry::GetInstanceFromRunNumber() - Initialized geometry with name <<EMCAL_COMPLETEV1>>\n");
      }
    }
    return AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1","EMCAL",mcname,mctitle) ;
  }
  else if ( runNumber >  176000 && runNumber <= 197692 )
  {
    // Complete EMCAL geometry, 12 SM. Year 2012 and on
    // The last 2 SM were not active, anyway they were there.
    
    if (showInfo)
    {
      if(!geoName.Contains("COMPLETE12SMV1") && geoName!="")
      {
        printf("AliEMCALGeometry::GetInstanceFromRunNumber() *** ATTENTION *** \n");
        printf("\t Specified geometry name <<%s>> for run %d is not considered! \n",geoName.Data(),runNumber);
        printf("\t In use <<EMCAL_COMPLETE12SMV1>>, check run number and year\n");
      }
      else 
      {
        printf("AliEMCALGeometry::GetInstanceFromRunNumber() - Initialized geometry with name <<EMCAL_COMPLETE12SMV1>>\n");
      }
    }
    return AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1","EMCAL",mcname,mctitle) ;
  }
  else // Run 2
  {
    // EMCAL + DCAL geometry, 20 SM. Year 2015 and on
    
    if (showInfo)
    {
      if(!geoName.Contains("DCAL_8SM") && geoName!="")
      {
        printf("AliEMCALGeometry::GetInstanceFromRunNumber() *** ATTENTION *** \n");
        printf("\t Specified geometry name <<%s>> for run %d is not considered! \n",geoName.Data(),runNumber);
        printf("\t In use <<EMCAL_COMPLETE12SMV1_DCAL_8SM>>, check run number and year\n");
      }
      else 
      {
        printf("AliEMCALGeometry::GetInstanceFromRunNumber() - Initialized geometry with name <<EMCAL_COMPLETE12SMV1_DCAL_8SM>>\n");
      }
    }
    return AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM","EMCAL",mcname,mctitle) ;
  }  
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
  Int_t id=0; // have to change from 0 to fNCells-1
  for( int i = 0 ; i < nSupMod; i++) {
    if(      GetSMType(i) == kEMCAL_Standard) id += fNCellsInSupMod;
    else if( GetSMType(i) == kEMCAL_Half)     id += fNCellsInSupMod/2;
    else if( GetSMType(i) == kEMCAL_3rd)      id += fNCellsInSupMod/3;
    else if( GetSMType(i) == kDCAL_Standard)  id += 2*fNCellsInSupMod/3;
    else if( GetSMType(i) == kDCAL_Ext)       id += fNCellsInSupMod/3;
    else {
      AliError(Form("Uknown SuperModule Type !!"));
    }
  }
  
  id += fNCellsInModule *nModule;
  id += fNPHIdiv *nIphi;
  id += nIeta;
  if( !CheckAbsCellId(id) ) {
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
  
  // Check if the indeces correspond to existing SM or tower indeces
  if(iphi    < 0 || iphi    >= AliEMCALGeoParams::fgkEMCALRows || 
     ieta    < 0 || ieta    >= AliEMCALGeoParams::fgkEMCALCols ||
     nSupMod < 0 || nSupMod >= GetNumberOfSuperModules()         )
  {
    AliDebug(1,Form("Wrong cell indexes : SM %d, column (eta) %d, row (phi) %d", nSupMod,ieta,iphi));
    return -1 ;
  }
  
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
  Int_t nphism = fEMCGeometry->GetNumberOfSuperModules()/2;
  for(i=0; i<nphism; i++) {
    if(phi>=fPhiBoundariesOfSM[2*i] && phi<=fPhiBoundariesOfSM[2*i+1]) {
      nSupMod = 2*i;
      if(eta < 0.0) nSupMod++;
      if( GetSMType(nSupMod) == kDCAL_Standard) {// Gap between DCAL
        if(TMath::Abs(eta) < GetNEta()/3*(GetEMCGeometry()->GetTrd1Angle())*TMath::DegToRad()) return kFALSE;
      }
      AliDebug(1,Form("eta %f phi %f(%5.2f) : nSupMod %i : #bound %i", eta,phi,phi*TMath::RadToDeg(), nSupMod,i));
      return kTRUE;
    }
  }
  return kFALSE;
}

///
/// Online mapping and numbering is the same for EMCal and DCal SMs but:
///  - DCal odd SM (13,15,17) has online cols: 16-47; offline cols 0-31.
///  - Even DCal SMs have the same numbering online and offline 0-31.
///  - DCal 1/3 SM (18,19), online rows 16-23; offline rows 0-7
///
/// Here shift the online cols or rows depending on the
/// super-module number to match the offline mapping.
///
/// \param sm: super module number of the channel/cell
/// \param iphi: row/phi cell index, modified for DCal
/// \param ieta: column/eta index, modified for DCal
///
//________________________________________________________________________________________________
void AliEMCALGeometry::ShiftOnlineToOfflineCellIndexes(Int_t sm, Int_t & iphi, Int_t & ieta) const 
{
  if ( sm == 13 || sm == 15 || sm == 17 )
  {
    // DCal odd SMs
    ieta -= 16; // Same cabling mapping as for EMCal, not considered offline.
  }
  else if ( sm == 18 || sm == 19 )
  {
    // DCal 1/3 SMs
    iphi -= 16; // Needed due to cabling mistake.
  }
}

///
/// Here shift the DCal online cols or rows depending on the
/// super-module number to match the online mapping. 
///
/// Reverse procedure to the one in the method above
/// ShiftOnlineToOfflineCellIndexes().
///
/// \param sm: super module number of the channel/cell
/// \param iphi: row/phi cell index, modified for DCal
/// \param ieta: column/eta index, modified for DCal
///
//________________________________________________________________________________________________
void AliEMCALGeometry::ShiftOfflineToOnlineCellIndexes(Int_t sm, Int_t & iphi, Int_t & ieta) const 
{
  if ( sm == 13 || sm == 15 || sm == 17 )
  {
    // DCal odd SMs
    ieta += 16; // Same cabling mapping as for EMCal, not considered offline.
  }
  else if ( sm == 18 || sm == 19 )
  {
    // DCal 1/3 SMs
    iphi += 16; // Needed due to cabling mistake.
  }
}

//________________________________________________________________________________________________
Bool_t AliEMCALGeometry::GetAbsCellIdFromEtaPhi(Double_t eta, Double_t phi, Int_t &absId) const
{

  // Nov 17,2006
  // stay here - phi problem as usual 
  static Int_t nSupMod=-1, i=0, ieta=-1, iphi=-1, etaShift=0, neta=-1, nphi=-1;
  static Double_t absEta=0.0, d=0.0, dmin=0.0, phiLoc=0;
  absId = nSupMod = - 1;
  if(SuperModuleNumberFromEtaPhi(eta, phi, nSupMod)) {
    // phi index first
    phi    = TVector2::Phi_0_2pi(phi);
    phiLoc = phi - fPhiCentersOfSMSec[nSupMod/2];
    nphi   = fPhiCentersOfCells.GetSize();
    if (     GetSMType(nSupMod) == kEMCAL_Half ) nphi  /= 2;
    else if( GetSMType(nSupMod) == kEMCAL_3rd )  nphi  /= 3;
    else if( GetSMType(nSupMod) == kDCAL_Ext )   nphi  /= 3;
    
    dmin   = TMath::Abs(fPhiCentersOfCells[0]-phiLoc);
    iphi   = 0;
    for(i=1; i<nphi; i++) {
      d = TMath::Abs(fPhiCentersOfCells[i] - phiLoc);
      if(d < dmin) {
        dmin = d;
        iphi = i;
      }
      //printf(" i %i : d %f : dmin %f : fPhiCentersOfCells[i] %f \n", i, d, dmin, fPhiCentersOfCells[i]);
    }
    // odd SM are turned with respect of even SM - reverse indexes
    AliDebug(2,Form(" iphi %i : dmin %f (phi %f, phiLoc %f ) ", iphi, dmin, phi, phiLoc));

    // eta index
    absEta   = TMath::Abs(eta);
    neta     = fCentersOfCellsEtaDir.GetSize();
    etaShift = iphi*neta;
    ieta     = 0;
    if( GetSMType(nSupMod) == kDCAL_Standard) ieta += 16; //jump 16 cells for DCSM
    dmin     = TMath::Abs(fEtaCentersOfCells[etaShift + ieta]-absEta);
    for(i= ieta+1 ; i<neta; i++) {
      d = TMath::Abs(fEtaCentersOfCells[i+etaShift] - absEta);
      if(d < dmin) {
        dmin = d;
        ieta = i;
      }
    }
    if( GetSMType(nSupMod) == kDCAL_Standard) ieta -= 16; //jump 16 cells for DCSM

    AliDebug(2,Form(" ieta %i : dmin %f (eta=%f) : nSupMod %i ", ieta, dmin, eta, nSupMod));
     
   //patch for mapping following alice convention  
    if(nSupMod%2 == 0) {// 47 + 16 -ieta for DCSM, 47 - ieta for others, revert the ordering on A side in order to keep convention.
      ieta = (neta -1)-ieta;
      if( GetSMType(nSupMod) == kDCAL_Standard) ieta -= 16; //recover cells for DCSM
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
  if(!CheckAbsCellId(absId)) return kFALSE;

  static Int_t tmp = absId;
  Int_t test = absId;
 
  for(nSupMod = -1; test >= 0; ) {
    nSupMod++;
    tmp = test;
    if(      GetSMType(nSupMod) == kEMCAL_Standard) test -= fNCellsInSupMod;
    else if( GetSMType(nSupMod) == kEMCAL_Half)     test -= fNCellsInSupMod/2;
    else if( GetSMType(nSupMod) == kEMCAL_3rd)      test -= fNCellsInSupMod/3;
    else if( GetSMType(nSupMod) == kDCAL_Standard)  test -= 2*fNCellsInSupMod/3;
    else if( GetSMType(nSupMod) == kDCAL_Ext)       test -= fNCellsInSupMod/3;
    else {
      AliError(Form("Uknown SuperModule Type !!"));
      return kFALSE;
    }
  }
  nModule = tmp / fNCellsInModule;
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
  if(      GetSMType(nSupMod) == kEMCAL_Half )  nphi = fNPhi/2; // halfSM
  else if( GetSMType(nSupMod) == kEMCAL_3rd  )  nphi = fNPhi/3; // 1/3 SM
  else if( GetSMType(nSupMod) == kDCAL_Ext   )  nphi = fNPhi/3; // 1/3 SM
  else                                          nphi = fNPhi;   // full SM
  
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
  // and SM of half (or one third) size in phi direction
 
  const Int_t kNphiIndex = fCentersOfCellsPhiDir.GetSize(); 
  Double_t  zshift = 0.5*GetDCALInnerEdge();
      
  static Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1;
  if(!CheckAbsCellId(absId)) return kFALSE;

  GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
  GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta); 
	
  //Get eta position. Careful with ALICE conventions (increase index decrease eta)	
  Int_t ieta2 = ieta;
  if(nSupMod%2 == 0) {
    ieta2 = (fCentersOfCellsEtaDir.GetSize()-1)-ieta;// 47-ieta, revert the ordering on A side in order to keep convention.
  }
  if( GetSMType(nSupMod) == kDCAL_Standard && nSupMod%2 ) ieta2 += 16; // DCAL revert the ordering on C side ...
  zr = fCentersOfCellsEtaDir.At(ieta2); 
  if( GetSMType(nSupMod) == kDCAL_Standard ) zr -= zshift; // DCAL shift (SMALLER SM)
  xr = fCentersOfCellsXDir.At(ieta2);

  //Get phi position. Careful with ALICE conventions (increase index increase phi)
  Int_t iphi2 = iphi;
  if( GetSMType(nSupMod) == kDCAL_Ext ) {
  if(nSupMod%2 != 0)   iphi2 = (kNphiIndex/3 -1)-iphi;  // 7-iphi [1/3SM], revert the ordering on C side in order to keep convention.
    yr = fCentersOfCellsPhiDir.At(iphi2 + kNphiIndex/3);
  } else if( GetSMType(nSupMod) == kEMCAL_Half ){
    if(nSupMod%2 != 0)   iphi2 = (kNphiIndex/2 -1)-iphi;  //11-iphi [1/2SM], revert the ordering on C side in order to keep convention.
    yr = fCentersOfCellsPhiDir.At(iphi2 + kNphiIndex/4);
  } else if( GetSMType(nSupMod) == kEMCAL_3rd ){
    if(nSupMod%2 != 0)   iphi2 = (kNphiIndex/3 -1)-iphi;  // 7-iphi [1/3SM], revert the ordering on C side in order to keep convention.
    yr = fCentersOfCellsPhiDir.At(iphi2 + kNphiIndex/3);
  } else {
    if(nSupMod%2 != 0)   iphi2 = (kNphiIndex   -1)-iphi;// 23-iphi, revert the ordering on C side in order to keep conventi
    yr = fCentersOfCellsPhiDir.At(iphi2);
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
  // and SM of half (or one third) size in phi direction
  
  const Int_t kNphiIndex = fCentersOfCellsPhiDir.GetSize();
  Double_t  zshift = 0.5*GetDCALInnerEdge();
  Int_t kDCalshift = 8;//wangml DCal cut first 8 modules(16 cells)
   
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
    ietam = (fCentersOfCellsEtaDir.GetSize()/2-1)-ietam;// 24-ietam, revert the ordering on A side in order to keep convention.
    if(nIeta == 0) nIeta = 1;
    else	   nIeta = 0;
  }
  if( GetSMType(nSupMod) == kDCAL_Standard && nSupMod%2) ietam += kDCalshift; // DCAL revert the ordering on C side ....
  mod = GetShishKebabModule(ietam);
  mod ->GetPositionAtCenterCellLine(nIeta, distEff, v); 
  xr = v.Y() - fParSM[0];
  zr = v.X() - fParSM[2];
  if( GetSMType(nSupMod) == kDCAL_Standard ) zr -= zshift; // DCAL shift (SMALLER SM)
 
  //Get phi position. Careful with ALICE conventions (increase index increase phi)
  Int_t iphi2 = iphi;
  if( GetSMType(nSupMod) == kDCAL_Ext ) {
     if(nSupMod%2 != 0)  iphi2 = (kNphiIndex/3 -1)-iphi;  // 7-iphi [1/3SM], revert the ordering on C side in order to keep convention.
     yr = fCentersOfCellsPhiDir.At(iphi2 + kNphiIndex/3);
   } else if( GetSMType(nSupMod) == kEMCAL_Half ){
     if(nSupMod%2 != 0)  iphi2 = (kNphiIndex/2 -1)-iphi;  //11-iphi [1/2SM], revert the ordering on C side in order to keep convention.
     yr = fCentersOfCellsPhiDir.At(iphi2 + kNphiIndex/2);
   } else if( GetSMType(nSupMod) == kEMCAL_3rd ){
     if(nSupMod%2 != 0)  iphi2 = (kNphiIndex/3 -1)-iphi;  // 7-iphi [1/3SM], revert the ordering on C side in order to keep convention.
     yr = fCentersOfCellsPhiDir.At(iphi2 + kNphiIndex/3);
   } else {
     if(nSupMod%2 != 0)  iphi2 = (kNphiIndex   -1)-iphi;// 23-iphi, revert the ordering on C side in order to keep convention.
     yr = fCentersOfCellsPhiDir.At(iphi2);
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
Bool_t AliEMCALGeometry::IsInEMCAL(Double_t x, Double_t y, Double_t z) const 
{
  // Checks whether point is inside the EMCal volume 
  if( IsInEMCALOrDCAL(x,y,z) == 1 ) return kTRUE;
  else return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliEMCALGeometry::IsInDCAL(Double_t x, Double_t y, Double_t z) const 
{
  // Checks whether point is inside the DCal volume
  if( IsInEMCALOrDCAL(x,y,z) == 2 ) return kTRUE;
  else return kFALSE;
}

//_____________________________________________________________________________
Int_t AliEMCALGeometry::IsInEMCALOrDCAL(Double_t x, Double_t y, Double_t z) const 
{
  // Checks whether point is inside the EMCal volume (included DCal), used in AliEMCALv*.cxx
  //
  // Code uses cylindrical approximation made of inner radius (for speed)
  //
  // Points behind EMCAl/DCal, i.e. R > outer radius, but eta, phi in acceptance 
  // are considered to inside

  Double_t r=sqrt(x*x+y*y);

  if ( r <= fEnvelop[0] ) return 0;
  else {
    Double_t theta = TMath::ATan2(r,z);
    Double_t eta;
    if(theta == 0)  eta = 9999;
    else            eta = -TMath::Log(TMath::Tan(theta/2.));
    if (eta < fArm1EtaMin || eta > fArm1EtaMax) return 0;

    Double_t phi = TMath::ATan2(y,x) * 180./TMath::Pi();
    if (phi < 0) phi += 360;  // phi should go from 0 to 360 in this case

    if (      phi >= fArm1PhiMin         && phi <= fEMCALPhiMax ) return 1;
    else if ( phi >= fDCALPhiMin         && phi <= fDCALStandardPhiMax && TMath::Abs(eta) > fDCALInnerExtandedEta ) return 2;
    else if ( phi > fDCALStandardPhiMax  && phi <= fDCALPhiMax  ) return 2;
    else return 0;
  } 
}

///
/// Provides shift-rotation matrix for EMCAL from externally set matrix or 
/// from TGeoManager
/// \return alignment matrix for a super module number
/// \param smod: super module number
///
//____________________________________________________________________________
const TGeoHMatrix * AliEMCALGeometry::GetMatrixForSuperModule(Int_t smod) const 
{	
  if(smod < 0 || smod > fEMCGeometry->GetNumberOfSuperModules()) 
    AliFatal(Form("Wrong supermodule index -> %d",smod));
		
  // Use matrices set externally
  if(!gGeoManager || (gGeoManager && fUseExternalMatrices))
  {
    if(fkSModuleMatrix[smod])
    {
      return fkSModuleMatrix[smod] ;
    }
    else
    {
      AliInfo("Stop:");
      printf("\t Can not find EMCAL misalignment matrixes\n") ;
      printf("\t Either import TGeoManager from geometry.root or \n");
      printf("\t read stored matrixes from AliESD Header:  \n") ;   
      printf("\t AliEMCALGeometry::SetMisalMatrixes(header->GetEMCALMisalMatrix()) \n") ;
      AliFatal("") ;
    }  
  }//external matrices
  
  // If gGeoManager exists, take matrix from it
  if(gGeoManager) return GetMatrixForSuperModuleFromGeoManager(smod);
  
  return 0 ;
}


///
/// Provides shift-rotation matrix for EMCAL from fkSModuleMatrix[smod]
/// Unsafe method, not to be used in reconstruction, just check there is 
/// something in the array of matrices without crashing, for EVE checks.
///
/// \return alignment matrix for a super module number
/// \param smod: super module number
///
//______________________________________________________________________________________
const TGeoHMatrix * AliEMCALGeometry::GetMatrixForSuperModuleFromArray(Int_t smod) const 
{	
  if(smod < 0 || smod > fEMCGeometry->GetNumberOfSuperModules()) 
    AliFatal(Form("Wrong supermodule index -> %d",smod));
  
  return fkSModuleMatrix[smod] ;
}


///
/// Provides shift-rotation matrix for EMCAL from the TGeoManager.
/// \return alignment matrix for a super module number
/// \param smod: super module number
///
//____________________________________________________________________________
const TGeoHMatrix * AliEMCALGeometry::GetMatrixForSuperModuleFromGeoManager(Int_t smod) const 
{  
  const Int_t buffersize = 255;
  char  path[buffersize] ;
  Int_t tmpType = -1;
  Int_t smOrder = 0;
  
  //Get the order for SM
  for( Int_t i = 0; i < smod+1; i++)
  {
    if(GetSMType(i) == tmpType) 
    {
      smOrder++;
    } 
    else 
    {
      tmpType = GetSMType(i);
      smOrder = 1;
    }
  } 
  
  Int_t   smType = GetSMType(smod);
  TString smName = "";
  
  if      ( smType == kEMCAL_Standard ) smName = "SMOD";
  else if ( smType == kEMCAL_Half )     smName = "SM10";
  else if ( smType == kEMCAL_3rd )      smName = "SM3rd";
  else if ( smType == kDCAL_Standard )  smName = "DCSM";
  else if ( smType == kDCAL_Ext )       smName = "DCEXT";
  else AliError("Unkown SM Type!!");
  
  snprintf(path,buffersize,"/ALIC_1/XEN1_1/%s_%d", smName.Data(), smOrder) ;
  
  if (!gGeoManager->cd(path))
    AliFatal(Form("Geo manager can not find path %s!\n",path));
  
  return gGeoManager->GetCurrentMatrix();
}


//__________________________________________________________________________________________________________________
void AliEMCALGeometry::RecalculateTowerPosition(Float_t drow, Float_t dcol, const Int_t sm, const Float_t depth,
                                                const Float_t misaligTransShifts[15], const Float_t misaligRotShifts[15], Float_t global[3]) const
{ 
  //Transform clusters cell position into global with alternative method, taking into account the depth calculation.
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
    Float_t zIs = 0;
    
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
        zIs = zIs + 2*dz*TMath::Cos(teta1);
      else
        zIs = zIs + 2*dz*TMath::Cos(teta1) + 2*dz*TMath::Sin(teta1)*TMath::Tan(teta1-0.75*TMath::DegToRad());
      
    }
    
    z0 = dz*(dcol-2*istrip+0.5);
    zb = (2*dz-z0-depth*TMath::Tan(teta1));
    
    z = zIs - zb*TMath::Cos(teta1);
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

//__________________________________________________________________________________________________________________
void AliEMCALGeometry::SetMisalMatrix(const TGeoHMatrix * m, Int_t smod) 
{
  // Method to set shift-rotational matrixes from ESDHeader
  // Move from header due to coding violations : Dec 2,2011 by PAI
  fUseExternalMatrices = kTRUE;

  if (smod >= 0 && smod < fEMCGeometry->GetNumberOfSuperModules()){
    if(!fkSModuleMatrix[smod]) fkSModuleMatrix[smod] = new TGeoHMatrix(*m) ; //Set only if not set yet
  } else AliFatal(Form("Wrong supermodule index -> %d",smod));
}

//__________________________________________________________________________________________________________________
Bool_t AliEMCALGeometry::IsDCALSM(Int_t iSupMod) const
{
  if( fEMCSMSystem[iSupMod] == kDCAL_Standard || fEMCSMSystem[iSupMod] == kDCAL_Ext ) return kTRUE;
  return kFALSE;
}

//__________________________________________________________________________________________________________________
Bool_t AliEMCALGeometry::IsDCALExtSM(Int_t iSupMod) const
{
  if( fEMCSMSystem[iSupMod] == kDCAL_Ext ) return kTRUE;
  return kFALSE;
}
