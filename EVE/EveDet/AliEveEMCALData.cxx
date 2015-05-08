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


#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TRefArray.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "AliRunLoader.h"
#include "AliESDVertex.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"

#include "AliEMCAL.h"
#include "AliEMCALLoader.h"
#include "AliEMCALHit.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"

#include "AliEveEMCALData.h"
#include "AliEveEMCALSModuleData.h"

class Riostream;
class TObject;
class TEveUtil;
class TEvePointSet;

class AliRun;
class AliESDEvent;
class AliLog;

class AliEMCAL;
class AliEMCALGeometry;

class AliEveEMCALSModule;

/// \cond CLASSIMP
ClassImp(AliEveEMCALData) ;
/// \endcond

///
/// Default constructor.
///
//______________________________________________________________________________
AliEveEMCALData::AliEveEMCALData():
  TObject(),
  TEveRefCnt(),
  fGeom(0x0),
  fNode(0x0),
//  fHMatrix(0),
  fESD(0x0),
  fRunLoader(0x0),
  fNsm(20),
  fSM(20),
//  fPoint(0),
  fClusterMom()
{
}

///
/// Constructor.
///
//______________________________________________________________________________
AliEveEMCALData::AliEveEMCALData(AliRunLoader* rl, TGeoNode* node): //, TGeoHMatrix* m):
  TObject(),
  TEveRefCnt(),
  fGeom(0x0),
  fNode(node),
//  fHMatrix(m),
  fESD(0x0),
  fRunLoader(rl),
  fNsm(20),
  fSM(20),
//  fPoint(0),
  fClusterMom()
{  
  InitEMCALGeom();
  
  CreateAllSModules();
}

///
/// Destructor.
///
//______________________________________________________________________________
AliEveEMCALData::~AliEveEMCALData()
{
  DeleteSuperModules();
  //  delete fPoint;
}

///
/// Copy constructor.
///
//______________________________________________________________________________
AliEveEMCALData::AliEveEMCALData(const AliEveEMCALData &edata) :
  TObject(edata),
  TEveRefCnt(edata),
  fGeom(edata.fGeom),
  fNode(edata.fNode),
//  fHMatrix(edata.fHMatrix),
  fESD(edata.fESD),
  fRunLoader(edata.fRunLoader),
  fNsm    (edata.fNsm),
  fSM     (edata.fSM),
//  fPoint(edata.fPoint),
  fClusterMom(edata.fClusterMom)
{
  CreateAllSModules();
}

///
/// Assignment operator.
///
//______________________________________________________________________________
AliEveEMCALData& AliEveEMCALData::operator=(const AliEveEMCALData &edata)
{
  if (this != &edata) 
  {
  }

  return *this;
}

///
/// Set data members for EMCAL geometry.
///
//______________________________________________________________________________
void AliEveEMCALData::InitEMCALGeom()
{
  if(fRunLoader && fRunLoader->GetAliRun())
  {
    AliEMCAL * emcal = (AliEMCAL*) fRunLoader->GetAliRun()->GetDetector("EMCAL");
    
    if ( emcal ) 
      fGeom  = (AliEMCALGeometry*) emcal->GetGeometry();
  }
  
  if(!fGeom)
  {
    // Use default geometry, Run2. In case of running Run1, explicitely create
    // the instance for the corresponding Run1 geometry in the macro executing the
    // display or before.
    
    AliInfo("Set EMCAL geometry to default");
    
    fGeom = AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  }
  
  if(!fGeom) AliFatal("EMCAL geometry pointer is NULL");
  
  // Get the number of super modules from the recovered nodes
  // in case OCDB file used does not contain same SM as the requested geometry.
  //fNsm = fGeom->GetNumberOfSuperModules();
  fNsm = fNode->GetNdaughters();
  
  if(fNsm != fGeom->GetNumberOfSuperModules()) 
    AliWarning(Form("Number of nodes (%d) is different to the number of expected super modules (%d)",
                    fNsm,fGeom->GetNumberOfSuperModules()));
}

///
/// Get geometrical information from hit/digit/cluster absolute id.
///
//______________________________________________________________________________
void AliEveEMCALData::GetGeomInfo(Int_t id, Int_t &iSupMod, Double_t& x, Double_t& y, Double_t& z)
{
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;

  // Geometry methods
  fGeom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);

  // Gives SuperModule and Tower numbers
  fGeom->RelPosCellInSModule(id, x, y, z);
}

///
/// Create all fNsm super modules.
///
//______________________________________________________________________________
void  AliEveEMCALData::CreateAllSModules()
{
  for(Int_t sm = 0; sm < fNsm; sm++)
  {
    if ( fSM[sm] == 0 ) 
      fSM[sm] = new AliEveEMCALSModuleData(sm, fGeom, fNode);//,fHMatrix);
  }
}

///
/// Delete all super module data.
///
//______________________________________________________________________________
void AliEveEMCALData::DeleteSuperModules()
{
  for (Int_t sm = 0; sm < 20; sm++)
  {
    fSM[sm] = 0;
    delete fSM[sm];
  }  
}

///
/// Get hit information from AliRunLoader.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadHits()
{
  /*
   // These are global coordinates !
   char form[1000];
   const char *selection = "";
   const char *varexp = "fX:fY:fZ";
   sprintf(form,"EMCAL Hits '%s'", selection);
   fPoint = new TEvePointSet(form);
   
   TEvePointSelector ps(t, fPoint, varexp, selection);
   ps.Select();
   
   if (fPoint->Size() == 0) {
   Warning("emcal_hits", Form("No hits match '%s'", selection));
   delete fPoint;
   //    return 0;
   }
   */
  
  fRunLoader->LoadHits("EMCAL");
  
  TTree * t = fRunLoader->GetTreeH("EMCAL",kFALSE);
 
  TObjArray *harr=NULL;
  
  TBranch *hbranch=t->GetBranch("EMCAL");
  
  hbranch->SetAddress(&harr);
  
  if(hbranch->GetEvent(0)) 
  {
    for(Int_t ih = 0; ih < harr->GetEntriesFast(); ih++) 
    {
      AliEMCALHit* hit =(AliEMCALHit*)harr->UncheckedAt(ih);
      
      if( !hit ) continue ; 
      
      Int_t    id  = hit->GetId();
      Double_t amp = hit->GetEnergy();

      AliDebug(1,Form("Hit info %d, energy %2.3f",id, amp));
            
      // Get global coordinates
      Double_t x = hit->X();
      Double_t y = hit->Y();
      Double_t z = hit->Z();
      
      Int_t iSupMod = 0;
      // These are local coordinates
      Double_t xl = 0.; Double_t yl = 0.; Double_t zl = 0.;

      // Get SM Id
      GetGeomInfo(id,iSupMod,xl,yl,zl);
      
      if(iSupMod >= fNsm ) continue ;
      
      fSM[iSupMod]->RegisterHit(id,iSupMod,amp,x,y,z);
//      if     ( iSupMod < 10 ) fSMfull [iSupMod]   ->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 12 ) fSMhalf [iSupMod-10]->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 18 ) fSMfullD[iSupMod-12]->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 20 ) fSMhalfD[iSupMod-18]->RegisterHit(id,iSupMod,amp,x,y,z);
    }
  }
  
  fRunLoader->UnloadHits("EMCAL");
  
  AliDebug(1,"Hits loop done");  
}

///
/// Get hit information from AliEMCALLoader.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadHitsFromEMCALLoader()
{      
  AliEMCALLoader * emcl = dynamic_cast<AliEMCALLoader*> (fRunLoader->GetDetectorLoader("EMCAL"));

  fRunLoader->LoadHits("EMCAL");
  
  AliEMCALHit* hit;
  
  // Fill array of hits                                                                        
  TClonesArray *hits = 0;//(TClonesArray*)emcl->Hits();
 
  TTree *treeH = emcl->TreeH();	
  
  if (!treeH) return ; 
  
  Int_t nTrack = treeH->GetEntries();  // TreeH has array of hits for every primary
  
  TBranch * branchH = treeH->GetBranch("EMCAL");
  //if(fHits)fHits->Clear();
  
  branchH->SetAddress(&hits);
  
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) 
  {
    branchH->GetEntry(iTrack);
    
    // Get hits from the list                                                                    
    for(Int_t ihit = 0; ihit< hits->GetEntries();ihit++)
    {
      hit = static_cast<AliEMCALHit *>(hits->At(ihit)) ;
      
      if(!hit) continue ;
      
      Int_t    id  = hit->GetId();
      Double_t amp = hit->GetEnergy();

      AliDebug(1,Form("Hit info %d, energy %2.3f",id, amp));
      
      // These are local coordinates
      Double_t xl = 0.; Double_t yl = 0.; Double_t zl = 0.;
     
      // Get global coordinates
      Double_t x = hit->X();
      Double_t y = hit->Y();
      Double_t z = hit->Z();
      
      Int_t iSupMod = 0;
      
      // Get SM Id
      GetGeomInfo(id,iSupMod,xl,yl,zl);
      
      if(iSupMod >= fNsm ) continue ;
      
      fSM[iSupMod]->RegisterHit(id,iSupMod,amp,x,y,z);
    }// hit loop

    hits->Clear();
  }// track loop
  
  fRunLoader->UnloadHits("EMCAL");
  
  AliDebug(1,"Hits loop done");  
}

///
/// Get digit information from AliRunLoader.
///
//______________________________________________________________________________
void  AliEveEMCALData::LoadDigits()
{
  fRunLoader->LoadDigits("EMCAL");
  
  TTree * dt = fRunLoader->GetTreeD("EMCAL", kFALSE);

  TClonesArray *digits = 0;
  dt->SetBranchAddress("EMCAL", &digits);
  dt->GetEntry(0);
  
  Int_t nEnt = digits->GetEntriesFast();

  AliDebug(1,Form("Number of digits %d",nEnt));
  
  AliEMCALDigit * dig;
  
  Double_t ampFlo   = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
  
  for (Int_t idig = 0; idig < nEnt; idig++)
  {
    dig = static_cast<AliEMCALDigit *>(digits->At(idig));
    
    if(!dig)
    {
      AliWarning(Form("Digit %d in array is null",idig));
      continue; 
    }
    
    id      = dig->GetId() ;       // cell (digit) label
    ampFlo  = dig->GetAmplitude(); // amplitude in cell (digit) ADC
                                   // acess OCDB and get calibration factor to GeV?	
    
    // Do not add too low ADC values (3 times pedestal)
    if(ampFlo <= 3 ) continue ;
    
    GetGeomInfo(id,iSupMod,x,y,z);
    
    if(iSupMod >= fNsm ) continue ;
    
    fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
  } // end loop digits
  
  fRunLoader->UnloadDigits("EMCAL");
  
  AliDebug(1,"Digits loop done");  
}

///
/// Get digit information from AliEMCALLoader.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadDigitsFromEMCALLoader()
{  
  AliEMCALLoader * emcl = dynamic_cast<AliEMCALLoader*> (fRunLoader->GetDetectorLoader("EMCAL"));

  fRunLoader->LoadDigits("EMCAL");
  
  AliEMCALDigit* dig;

  // Fill array of digits                                                                        
  TClonesArray *digits = (TClonesArray*)emcl->Digits();
  
  // Get digits from the list  
  Double_t ampFlo   = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
  
  for(Int_t idig = 0; idig< digits->GetEntries();idig++)
  {
    dig = static_cast<AliEMCALDigit *>(digits->At(idig)) ;
    
    if(!dig)
    {
      AliWarning(Form("Digit %d in array is null",idig));
      continue; 
    }
    
    id      = dig->GetId() ;       // cell (digit) label
    ampFlo  = dig->GetAmplitude(); // amplitude in cell (digit) ADC
                                   // acess OCDB and get calibration factor to GeV?	
    
    // Do not add too low ADC values (3 times pedestal)
    if(ampFlo <= 3 ) continue ;
    
    AliDebug(1,Form("Digit info %d, energy %2.3f",id, ampFlo));
    
    GetGeomInfo(id,iSupMod,x,y,z);
    
    if(iSupMod >= fNsm ) continue ;
    
    fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
  } // end loop on digits
 
  fRunLoader->UnloadDigits("EMCAL");
  
  AliDebug(1,"Digits loop done");    
}

///
/// Get digit information from ESDs.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadDigitsFromESD()
{
  AliESDCaloCells &cells= *(fESD->GetEMCALCells());
  
  Int_t ncell = cells.GetNumberOfCells() ;  
  
  AliDebug(1,Form("Number of ESD CaloCells %d",ncell));
  
  Int_t iSupMod =  0 ;
  
  Double_t x, y, z;
  
  // Extract digit information from the ESDs
  for (Int_t icell=  0; icell <  ncell; icell++) 
  {
    Int_t id         = cells.GetCellNumber(icell);
    Double_t ampFlo  = cells.GetAmplitude (icell); // GeV
    
    AliDebug(1,Form("CaloCell %d, ID %d, energy %2.3f",icell,id,ampFlo));
    
    GetGeomInfo(id,iSupMod,x,y,z);
    
    if(iSupMod >= fNsm ) continue ;
    
    fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
  } // end loop cells
  
  AliDebug(1,"CaloCells loop done");  
}

///
/// Get rec point information from AliRunLoader.
/// To be improved!!!, size and shape of cluster to be implemented.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadRecPoints()
{
  fRunLoader->LoadRecPoints("EMCAL");

  TTree * t = fRunLoader->GetTreeR("EMCAL",kFALSE);
  
  // From TTreeR
  TObjArray * carr = NULL;
  TBranch * cbranch=t->GetBranch("EMCALECARP");
  cbranch->SetAddress(&carr);
  
  if(!cbranch->GetEvent(0)) return;
  
  for(Int_t ic = 0; ic < carr->GetEntriesFast(); ic++) 
  {
    AliEMCALRecPoint * rp =(AliEMCALRecPoint*)carr->UncheckedAt(ic);
    
    if(!rp) continue ;
    
    Int_t    iSupMod = rp->GetSuperModuleNumber();
    Double_t amp     = (Double_t)rp->GetEnergy();
    
    if(iSupMod >= fNsm ) continue ;
    
    AliDebug(1,Form("RecPoint info, Id %d, energy %2.3f",rp->GetAbsId(0), amp));

    TVector3 lpos;
    rp->GetLocalPosition(lpos);
        
    fSM[iSupMod]->RegisterCluster(iSupMod,amp,lpos[0],lpos[1],lpos[2]);
  }

  fRunLoader->UnloadRecPoints("EMCAL");
  
  AliDebug(1,"RecPoints loop done");    
}

///
/// Get rec point information from AliEMCALLoader.
/// To be improved!!!, size and shape of cluster to be implemented.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadRecPointsFromEMCALLoader()
{  
  AliEMCALLoader * emcl = dynamic_cast<AliEMCALLoader*> (fRunLoader->GetDetectorLoader("EMCAL"));
  
  fRunLoader->LoadRecPoints("EMCAL");
  
  AliEMCALRecPoint* rp = 0;
  
  // Fill array of clusters                                                                        
  TClonesArray *clusters = (TClonesArray*)emcl->RecPoints();
  
  // Get clusters from the list                                                                    
  for(Int_t iclu = 0; iclu < clusters->GetEntries(); iclu++)
  {
    rp = static_cast<AliEMCALRecPoint *>(clusters->At(iclu)) ;
        
    if(!rp) continue ;
    
    Int_t    iSupMod = rp->GetSuperModuleNumber();
    Double_t amp     = (Double_t)rp->GetEnergy();
    
    if(iSupMod >= fNsm ) continue ;
    
    AliDebug(1,Form("RecPoint info, Id %d, energy %2.3f",rp->GetAbsId(0), amp));
    
    TVector3 lpos;
    rp->GetLocalPosition(lpos);
    
    fSM[iSupMod]->RegisterCluster(iSupMod,amp,lpos[0],lpos[1],lpos[2]);
  }
  
  fRunLoader->UnloadRecPoints("EMCAL");
  
  AliDebug(1,"RecPoints loop done");    
}

///
/// Get cluster information from ESDs.
/// To be improved!!!, size and shape of cluster to be implemented.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadRecPointsFromESD()
{
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
  //  Int_t iSM =  0 ;
  //  Int_t iT  =  0 ;
  //  Int_t iIp =  0 ;
  //  Int_t iIe =  0 ;
  //  Double_t xd, yd, zd;
  //  Float_t pos[3] ; 
  
  // Get reconstructed vertex position
  AliESDVertex* primVertex =(AliESDVertex*) fESD->GetVertex();
  Double_t vertexPosition[3] ; 
  primVertex->GetXYZ(vertexPosition) ; 
  
  // Get the CaloClusters
  // select EMCAL clusters only 
  TRefArray * caloClusters  = new TRefArray();
  fESD->GetEMCALClusters(caloClusters);
  
  Int_t nclus = caloClusters->GetEntries();

  AliDebug(1,Form("Number of clusters %d",nclus));
  
  for (Int_t iclus =  0; iclus <  nclus; iclus++) 
  {
    AliESDCaloCluster *clus = (AliESDCaloCluster *) caloClusters->At(iclus) ; 
    //Get the cluster info
    
    Double_t energy = clus->E() ;  
    Double_t disp   = clus->GetDispersion() ;
    
    //clus->GetPosition(pos) ; // Global position
    
    clus->GetMomentum(fClusterMom,vertexPosition);

    Double_t eta = fClusterMom.Eta();
    Double_t phi = ( (fClusterMom.Phi()) < 0) ? (fClusterMom.Phi()) + 2. * TMath::Pi() : (fClusterMom.Phi());
    
    Int_t mult = clus->GetNCells() ;
    
    AliDebug(2,Form("In cluster %d, ncells %d, energy %2.2f, disp %2.2f, eta %2.2f, phi %2.2f",
                    iclus,mult,energy,disp,eta,phi));
    
    Int_t clusId = clus->GetCellsAbsId()[0];
    
    GetGeomInfo(clusId,iSupMod,x,y,z);
    
    if(iSupMod >= fNsm ) continue ;
    
//    //******** Not used yet but will come  ********
//    // AliESDCaloCells &cells= *(fESD->GetEMCALCells());
//    Int_t     digMult = clus->GetNCells() ;
//    UShort_t *digID   = clus->GetCellsAbsId() ;
//    for(Int_t i=0; i<digMult; i++)
//    {
//      fGeom->RelPosCellInSModule(digID[i], xd, yd, zd);
//      fGeom->GetCellIndex(digID[i],iSM,iT,iIp,iIe);
//    } // end digit loop
//      //*********************************************
  
    fSM[iSupMod]->RegisterCluster(iSupMod,energy,x,y,z);
  } // end cluster loop
}

///
/// Return super module data.
///
//______________________________________________________________________________
AliEveEMCALSModuleData* AliEveEMCALData::GetSModuleData(Int_t sm)
{
  if (sm < 0 || sm > fNsm) 
  {
    AliWarning(Form("The number of super modules must be lower or equal to %d, requested %d",fNsm,sm));
    return 0;
  }
  
  return fSM[sm];
}

