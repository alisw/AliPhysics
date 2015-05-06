#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TRefArray.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "AliRunLoader.h"
#include "AliEMCAL.h"
#include "AliEMCALLoader.h"
#include "AliESDVertex.h"
#include "AliEMCALHit.h"
#include "AliEMCALDigit.h"

#include "AliEMCALRecPoint.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"

#include "AliEveEMCALData.h"
#include "AliEveEMCALSModuleData.h"

class Riostream;
class TObject;
class TEveUtil;
class TEvePointSet;
class AliRun;
class AliESDEvent;
class AliEMCAL;
class AliEMCALGeometry;
class AliEveEMCALSModule;
class AliLog;

/// \cond CLASSIMP
ClassImp(AliEveEMCALData) ;
/// \endcond

///
/// Default constructor
///
//______________________________________________________________________________
AliEveEMCALData::AliEveEMCALData():
  TObject(),
  TEveRefCnt(),
  fEmcal(0x0),
  fGeom(0x0),
  fNode(0x0),
  fHMatrix(0),
  fTree(0x0),
  fESD(0x0),
  fNsm(20),
  fNsmfull(10),
  fNsmhalf(2),  
  fNsmfullD(6),
  fNsmhalfD(2),
  fSM(20),
  fSMfull(10),
  fSMhalf(2),
  fSMfullD(6),
  fSMhalfD(2),
  fRunLoader(0),
  fDebug(0),
  fPoint(0)
{
  CreateAllSModules();
}

///
/// Constructor
///
//______________________________________________________________________________
AliEveEMCALData::AliEveEMCALData(AliRunLoader* rl, TGeoNode* node, TGeoHMatrix* m):
  TObject(),
  TEveRefCnt(),
  fEmcal(0x0),
  fGeom(0x0),
  fNode(node),
  fHMatrix(m),
  fTree(0x0),
  fESD(0x0),
  fNsm(20),
  fNsmfull(10),
  fNsmhalf(2),  
  fNsmfullD(6),
  fNsmhalfD(2),
  fSM(20),
  fSMfull(10),
  fSMhalf(2),
  fSMfullD(6),
  fSMhalfD(2),
  fRunLoader(rl),
  fDebug(0),
  fPoint(0)
{
  InitEMCALGeom(rl);
  CreateAllSModules();
}

///
/// Destructor
///
//______________________________________________________________________________
AliEveEMCALData::~AliEveEMCALData()
{
  DeleteSuperModules();
  delete fTree;
  // delete fEmcal; // deleted by run-loader
  delete fGeom;
  delete fNode;
  delete fHMatrix;
  delete fPoint;
}

///
/// Copy constructor
///
//______________________________________________________________________________
AliEveEMCALData::AliEveEMCALData(const AliEveEMCALData &edata) :
  TObject(edata),
  TEveRefCnt(edata),
  fEmcal(edata.fEmcal),
  fGeom(edata.fGeom),
  fNode(edata.fNode),
  fHMatrix(edata.fHMatrix),
  fTree(edata.fTree),
  fESD(edata.fESD),
  fNsm     (edata.fNsm     ),
  fNsmfull (edata.fNsmfull ),
  fNsmhalf (edata.fNsmhalf ),
  fNsmfullD(edata.fNsmfullD),
  fNsmhalfD(edata.fNsmhalfD),
  fSM     (edata.fSM),
  fSMfull (edata.fSMfull ),
  fSMhalf (edata.fSMhalf ),  
  fSMfullD(edata.fSMfullD),
  fSMhalfD(edata.fSMhalfD),
  fRunLoader(edata.fRunLoader),
  fDebug(edata.fDebug),
  fPoint(edata.fPoint)
{
  InitEMCALGeom(edata.fRunLoader);
  CreateAllSModules();
}

///
/// Assignment operator
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
/// Set digit-tree to be used for digit retrieval. 
/// Data is loaded on demand.
///
//______________________________________________________________________________
void AliEveEMCALData::SetTree(TTree* const tree)
{
  fTree = tree;
}

///
/// Set esd
///
//______________________________________________________________________________
void AliEveEMCALData::SetESD(AliESDEvent* const esd)
{
  fESD = esd;
}

///
/// Set node
///
//______________________________________________________________________________
void AliEveEMCALData::SetNode(TGeoNode* const node)
{
  fNode = node;
}

///
/// Set data members for EMCAL geometry
///
//______________________________________________________________________________
void AliEveEMCALData::InitEMCALGeom(AliRunLoader* const rl)
{
  fEmcal = (AliEMCAL*) rl->GetAliRun()->GetDetector("EMCAL");
  fGeom  = (AliEMCALGeometry*) fEmcal->GetGeometry();
  
  if(!fGeom) AliFatal("EMCAL geometry pointer is NULL");
  
  // Get the number of super modules from geometry
  fNsm = fGeom->GetNumberOfSuperModules();
}

///
/// Get geometrical information from hit/digit/cluster absolute id
///
//______________________________________________________________________________
void AliEveEMCALData::GetGeomInfo(Int_t id, Int_t &iSupMod, Double_t& x, Double_t& y, Double_t& z)
{
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;

  //Geometry methods
  fGeom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);

  //Gives SuperModule and Tower numbers
  fGeom->RelPosCellInSModule(id, x, y, z);
}

///
/// Create all fNsm super modules
///
//______________________________________________________________________________
void  AliEveEMCALData::CreateAllSModules()
{
  for(Int_t sm = 0; sm < fNsm; sm++)
    CreateSModule(sm);
}

///
/// Create super-module-data for SM if it does not exist already.
///
//______________________________________________________________________________
void  AliEveEMCALData::CreateSModule(Int_t sm)
{
  if     (fSM[sm] == 0)                    fSM     [sm]    = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
  if     (fSMfull [sm]    == 0 && sm < 10) fSMfull [sm]    = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
  else if(fSMhalf [sm-10] == 0 && sm < 12) fSMhalf [sm-10] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
  else if(fSMfullD[sm-12] == 0 && sm < 18) fSMfullD[sm-12] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
  else if(fSMhalfD[sm-18] == 0 && sm < 20) fSMhalfD[sm-18] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
}

///
/// Drop data of all existing sectors.
///
//______________________________________________________________________________
void AliEveEMCALData::DropAllSModules()
{
  for (Int_t sm = 0; sm < fNsm; sm++) 
  {
    if (fSM[sm] != 0)
      fSM[sm]->DropData();
  }
}

///
/// Delete all super module data
///
//______________________________________________________________________________
void AliEveEMCALData::DeleteSuperModules()
{
  for (Int_t sm = 0; sm < fNsm; sm++)
  {
    fSM[sm] = 0;
    delete fSM[sm];
  }
  
  for(Int_t smf = 0; smf < fNsmfull; smf++) 
  {
    fSMfull[smf] = 0;
    delete fSMfull[smf];
  }
  
  for(Int_t smh = 0; smh < fNsmhalf; smh++)
  {
    fSMhalf[smh] = 0;
    delete fSMhalf[smh];
  }

  for(Int_t smd = 0; smd < fNsmfullD; smd++) 
  {
    fSMfullD[smd] = 0;
    delete fSMfullD[smd];
  }
  
  for(Int_t smh = 0; smh < fNsmhalfD; smh++)
  {
    fSMhalfD[smh] = 0;
    delete fSMhalfD[smh];
  }
}

///
/// Get hit information from RunLoader
///
//______________________________________________________________________________
void AliEveEMCALData::LoadHits(TTree* const t)
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
      
      fSM[iSupMod]->RegisterHit(id,iSupMod,amp,x,y,z);
//      if     ( iSupMod < 10 ) fSMfull [iSupMod]   ->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 12 ) fSMhalf [iSupMod-10]->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 18 ) fSMfullD[iSupMod-12]->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 20 ) fSMhalfD[iSupMod-18]->RegisterHit(id,iSupMod,amp,x,y,z);
    }
  }
}

///
/// Get hit information from EMCAL Loader
///
//______________________________________________________________________________
void AliEveEMCALData::LoadHitsFromEMCALLoader(AliEMCALLoader* const emcl)
{  
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
      
      fSM[iSupMod]->RegisterHit(id,iSupMod,amp,x,y,z);
      
//      if     ( iSupMod < 10 ) fSMfull [iSupMod]   ->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 12 ) fSMhalf [iSupMod-10]->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 18 ) fSMfullD[iSupMod-12]->RegisterHit(id,iSupMod,amp,x,y,z);
//      else if( iSupMod < 20 ) fSMhalfD[iSupMod-18]->RegisterHit(id,iSupMod,amp,x,y,z);
    }//hit loop

    hits->Clear();
  }// track loop
}

///
/// Get digit information from RunLoader
///
//______________________________________________________________________________
void AliEveEMCALData::LoadDigits(TTree *t)
{
  TClonesArray *digits = 0;
  t->SetBranchAddress("EMCAL", &digits);
  t->GetEntry(0);
  
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
    
    GetGeomInfo(id,iSupMod,x,y,z);
    
    fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
//	  fSM[iSupMod]->SaveDigit(dig);
//    if     ( iSupMod < 10 ) fSMfull [iSupMod]   ->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
//    else if( iSupMod < 12 ) fSMhalf [iSupMod-10]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
//    else if( iSupMod < 18 ) fSMfullD[iSupMod-12]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
//    else if( iSupMod < 20 ) fSMhalfD[iSupMod-18]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);  
  } // end loop digits

  AliDebug(1,("Digits loop done"));
}

///
/// Get digit information from EMCAL Loader
///
//______________________________________________________________________________
void AliEveEMCALData::LoadDigitsFromEMCALLoader(AliEMCALLoader* const emcl)
{
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
    
    AliDebug(1,Form("Digit info %d, energy %2.3f",id, ampFlo));
    
    GetGeomInfo(id,iSupMod,x,y,z);
    
    fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
    
  } // end loop on digits
}

///
/// Get digit information from esd
///
//______________________________________________________________________________
void AliEveEMCALData::LoadDigitsFromESD()
{
  AliESDCaloCells &cells= *(fESD->GetEMCALCells());
  Int_t ncell = cells.GetNumberOfCells() ;  
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
  
  // Extract digit information from the ESDs
  for (Int_t icell=  0; icell <  ncell; icell++) 
  {
    Int_t id         = cells.GetCellNumber(icell);
    Double_t ampFlo  = cells.GetAmplitude (icell); // GeV
    
    GetGeomInfo(id,iSupMod,x,y,z);
    
    fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
    
    if     ( iSupMod < 10 ) fSMfull [iSupMod]   ->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
    else if( iSupMod < 12 ) fSMhalf [iSupMod-10]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
    else if( iSupMod < 18 ) fSMfullD[iSupMod-12]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
    else if( iSupMod < 20 ) fSMhalfD[iSupMod-18]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
  } // end loop cells
}

///
/// Get rec point information from RunLoader.
/// To be improved, size and shape of cluster to be implemented.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadRecPoints(TTree* const t)
{
  //*************************************************
  // To be improved !!!!!
  // Size and shape of cluster to be implemented
  // 
  //*************************************************
  
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

    AliDebug(1,Form("RecPoint info, Id %d, energy %2.3f",rp->GetAbsId(0), amp));

    TVector3 lpos;
    rp->GetLocalPosition(lpos);
    
    fSM[iSupMod]->RegisterCluster(iSupMod,amp,lpos[0],lpos[1],lpos[2]);
  }
}

///
/// Get rec point information from EMCalLoader.
/// To be improved, size and shape of cluster to be implemented.
///
//______________________________________________________________________________
void AliEveEMCALData::LoadRecPointsFromEMCALLoader(AliEMCALLoader* const emcl)
{
  //*************************************************
  // To be improved !!!!!
  // Size and shape of cluster to be implemented
  // 
  //*************************************************
  
  // From EMCALLoader
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
    
    AliDebug(1,Form("RecPoint info, Id %d, energy %2.3f",rp->GetAbsId(0), amp));
    
    TVector3 lpos;
    rp->GetLocalPosition(lpos);
    
    fSM[iSupMod]->RegisterCluster(iSupMod,amp,lpos[0],lpos[1],lpos[2]);
  }
}

///
/// Get cluster information from esd
///
//______________________________________________________________________________
void AliEveEMCALData::LoadRecPointsFromESD()
{
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
  Int_t iSM =  0 ;
  Int_t iT  =  0 ;
  Int_t iIp =  0 ;
  Int_t iIe =  0 ;
  Double_t xd, yd, zd;
  Float_t pos[3] ; 
  
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
    
    clus->GetPosition(pos) ; // Global position
    TVector3 vpos(pos[0],pos[1],pos[2]) ;
    
    TLorentzVector p4 ;
    clus->GetMomentum(p4,vertexPosition);
    
    TVector3 p3;
    p3.SetXYZ(p4[0],p4[1],p4[2]);
    
    Double_t eta = p3.Eta();
    Double_t phi = ( (p3.Phi()) < 0) ? (p3.Phi()) + 2. * TMath::Pi() : (p3.Phi());
    
    Int_t mult = clus->GetNCells() ;
    
    AliDebug(2,Form("In cluster %d, ncells %d, energy %2.2f, disp %2.2f, eta %2.2f, phi %2.2f",
                    iclus,mult,energy,disp,eta,phi));
    
    Int_t clusId = 0;
    fGeom->GetAbsCellIdFromEtaPhi(eta,phi,clusId);
    
    AliDebug(2,Form("Cluster AbsId %d, x %2.2f, y %2.2f, z %2.2f",
                    clusId,pos[0],pos[1],pos[2]));
    
    GetGeomInfo(clusId,iSupMod,x,y,z);
    
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
  
    fSM[iSupMod]->RegisterCluster(iSM,energy,x,y,z);
    
  } // end cluster loop
}

///
/// Return super module data
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

///
/// Get raw information
///
/// To be implemented !
//
//______________________________________________________________________________
void AliEveEMCALData::LoadRaw() const
{
}

