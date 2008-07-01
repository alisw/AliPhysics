//*********************************************************************
// - AliEVE implementation -
// Fill containers for visualisation of EMCAL data structures
//    - read and store MC Hits
//    - read and store digits from esds or runloader
//    - read and store clusters from esds or runloader 
//
// Author: Magali Estienne (magali.estienne@cern.ch)
// June 30 2008
//*********************************************************************

#include <TObject.h>
#include <TEveUtil.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TRefArray.h>
#include <TClonesArray.h>
#include <TEvePointSet.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliEMCALLoader.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliEMCAL.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALHit.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"

#include "AliEveEMCALData.h"
#include "AliEveEMCALSModule.h"
#include "AliEveEMCALSModuleData.h"

ClassImp(AliEveEMCALData)

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
  fNsm(12),
  fNsmfull(10),
  fNsmhalf(2),
  fSM(12),
  fSMfull(10),
  fSMhalf(2),
  fRunLoader(0),
  fDebug(0),
  fPoint(0)
{
  
  //
  // Constructor
  //
  //   for(Int_t i=0; i<12; i++)
  //     fSM[i] = 0x0;
  //   for(Int_t i=0; i<10; i++)
  //     fSMfull[i] = 0x0;
  //   fSMhalf[0] = fSMhalf[1] = 0x0;
  //  InitEMCALGeom();
  CreateAllSModules();


}

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
  fNsm(12),
  fNsmfull(10),
  fNsmhalf(2),
  fSM(12),
  fSMfull(10),
  fSMhalf(2),
  fRunLoader(rl),
  fDebug(0),
  fPoint(0)
{
  
  //
  // Constructor
  //
//   for(Int_t i=0; i<12; i++)
//     fSM[i] = 0x0;
//   for(Int_t i=0; i<10; i++)
//     fSMfull[i] = 0x0;
//   fSMhalf[0] = fSMhalf[1] = 0x0;
  InitEMCALGeom(rl);
  CreateAllSModules();


}

//______________________________________________________________________________
AliEveEMCALData::~AliEveEMCALData()
{
  //
  // Destructor
  //

  DeleteSuperModules();
  delete fTree;
  delete fEmcal;
  delete fGeom;
  delete fNode;
  delete fHMatrix;
  delete fPoint;
}

//______________________________________________________________________________
void AliEveEMCALData::Reset()
{

}

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
  fNsm(edata.fNsm),
  fNsmfull(edata.fNsmfull),
  fNsmhalf(edata.fNsmhalf),
  fSM(edata.fSM),
  fSMfull(edata.fSMfull),
  fSMhalf(edata.fSMhalf),
  fRunLoader(edata.fRunLoader),
  fDebug(edata.fDebug),
  fPoint(edata.fPoint)
{
  //
  // Copy constructor
  //
  InitEMCALGeom(edata.fRunLoader);
  CreateAllSModules();
}

//______________________________________________________________________________
AliEveEMCALData& AliEveEMCALData::operator=(const AliEveEMCALData &edata)
{
  //
  // Assignment operator
  //

  if (this != &edata) {

  }

  return *this;

}

//______________________________________________________________________________
void AliEveEMCALData::SetTree(TTree* tree)
{

  // Set digit-tree to be used for digit retrieval. Data is loaded on
  // demand.

  fTree = tree;

}

//______________________________________________________________________________
void AliEveEMCALData::SetESD(AliESDEvent* esd)
{
  fESD = esd;
}

//______________________________________________________________________________
void AliEveEMCALData::SetNode(TGeoNode* node)
{
  fNode = node;
}

//______________________________________________________________________________
void AliEveEMCALData::InitEMCALGeom(AliRunLoader* rl)
{
//   // Handle an individual point selection from GL.

//   fParent->SpawnEditor();
  fEmcal = (AliEMCAL*) rl->GetAliRun()->GetDetector("EMCAL");
  fGeom  = (AliEMCALGeometry*) fEmcal->GetGeometry();

}

//______________________________________________________________________________
void AliEveEMCALData::GetGeomInfo(Int_t id, Int_t &iSupMod, Double_t& x, Double_t& y, Double_t& z)
{
  // Get geometrical information from hit/digit/cluster absolute id

  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;
  Int_t iphi    =  0 ;
  Int_t ieta    =  0 ;

  //Geometry methods
  fGeom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
  //Gives SuperModule and Tower numbers
  fGeom->RelPosCellInSModule(id, x, y, z);

}

//______________________________________________________________________________
void  AliEveEMCALData::CreateAllSModules()
{

  // Create all fNsm super modules
  for(Int_t sm = 0; sm < fNsm; sm++)
    CreateSModule(sm);

}

//______________________________________________________________________________
void  AliEveEMCALData::CreateSModule(Int_t sm)
{
  // Create super-module-data for sm if it does not exist already.
  if(fSM[sm] == 0) fSM[sm] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
  if(fSMfull[sm] == 0 && sm < 10) fSMfull[sm] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
  if(fSMhalf[sm-10] == 0 && sm > 10) fSMhalf[sm-10] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
}

//______________________________________________________________________________
void AliEveEMCALData::DropAllSModules()
{
  // Drop data of all existing sectors.

  for (Int_t sm = 0; sm < fNsm; sm++) {
    if (fSM[sm] != 0)
      fSM[sm]->DropData();
  }
}

//______________________________________________________________________________
void AliEveEMCALData::DeleteSuperModules()
{
  //
  // delete all super module data
  //

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
  
}

//______________________________________________________________________________
void AliEveEMCALData::LoadHits(TTree* t)
{

  // With TTreeH

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


//   TObjArray *harr=NULL;
//   TBranch *hbranch=t->GetBranch("EMCAL");
//   hbranch->SetAddress(&harr);
  
//   if(hbranch->GetEvent(0)) {
//     for(Int_t ih = 0; ih < harr->GetEntriesFast(); ih++) {
//       AliEMCALHit* hit =(AliEMCALHit*)harr->UncheckedAt(ih);
//       if(hit != 0){
// 	cout << "Hit info " << hit->GetId() << " " << hit->GetEnergy() << endl;
// 	fSM[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
//       }
//     }
//   }
//   //********************************
//   // To be completed !!!
//   //********************************

}

//______________________________________________________________________________
void AliEveEMCALData::LoadHitsFromEMCALLoader(AliEMCALLoader* emcl)
{

  // Need to know the local position in each sm

  // With EMCALLoader
  AliEMCALHit* hit;
  
  //Fill array of hits                                                                        
  TClonesArray *hits = (TClonesArray*)emcl->Hits();

  //Get hits from the list                                                                    
  for(Int_t ihit = 0; ihit< hits->GetEntries();ihit++){
    //cout<<">> idig "<<idig<<endl;                                                             
    hit = static_cast<AliEMCALHit *>(hits->At(ihit)) ;
    
    if(hit != 0){
      cout << "Hit info " << hit->GetId() << " " << hit->GetEnergy() << endl;
    }
  }
  
  //********************************
  // To be completed !!!
  //********************************



}

//______________________________________________________________________________
void AliEveEMCALData::LoadDigits(TTree *t)
{
  //
  // load digits from the TreeD
  //

  TClonesArray *digits = 0;
  t->SetBranchAddress("EMCAL", &digits);
  t->GetEntry(0);

  Int_t nEnt = digits->GetEntriesFast();
  cout << "nEnt: " << nEnt << endl;
  AliEMCALDigit * dig;

  Float_t amp   = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Double_t x, y, z;

  for (Int_t idig = 0; idig < nEnt; idig++)
    {
      dig = static_cast<AliEMCALDigit *>(digits->At(idig));
      
      if(dig != 0) {
	id   = dig->GetId() ; //cell (digit) label
	amp  = dig->GetAmp(); //amplitude in cell (digit)
	
	GetGeomInfo(id,iSupMod,x,y,z);

	fSM[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
//	fSM[iSupMod]->SaveDigit(dig);
// 	if(iSupMod<fNsmfull) fSMfull[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
// 	if(iSupMod>fNsmfull) fSMhalf[iSupMod-10]->RegisterDigit(id,iSupMod,amp,x,y,z);
      }
      else {
	cout << "Digit object empty" << endl;
	return;
      }
    } // end loop digits
  cout << "after loop on digits !" << endl;
}

//______________________________________________________________________________
void AliEveEMCALData::LoadDigitsFromEMCALLoader(AliEMCALLoader* emcl)
{

  //
  // load digits from EMCAL Loader
  //

  AliEMCALDigit* dig;
  
  //Fill array of digits                                                                        
  TClonesArray *digits = (TClonesArray*)emcl->Digits();
  
  //Get digits from the list  
  
  Float_t amp   = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
                                               
   for(Int_t idig = 0; idig< digits->GetEntries();idig++){

     dig = static_cast<AliEMCALDigit *>(digits->At(idig)) ;

     if(dig != 0){
       if(fDebug>1) cout << "Digit info " << dig->GetId() << " " << dig->GetAmp() << endl;
       id   = dig->GetId() ; //cell (digit) label
       amp  = dig->GetAmp(); //amplitude in cell (digit)
       
       GetGeomInfo(id,iSupMod,x,y,z);
       
       fSM[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
     }
      else {
	cout << "Digit object empty" << endl;
	return;
      }
   } // end loop on digits
   
}

//______________________________________________________________________________
void AliEveEMCALData::LoadDigitsFromESD()
{
  //
  // Get digit information from esd
  //
  
  AliESDCaloCells &cells= *(fESD->GetEMCALCells());
  Int_t ncell = cells.GetNumberOfCells() ;  
  Int_t type = cells.GetType();
  Float_t pos[3] ; 
  Int_t iSupMod =  0 ;
  Double_t x, y, z;

  // Extract digit information from the ESDs
  for (Int_t icell=  0; icell <  ncell; icell++) 
    {
      Int_t id   = cells.GetCellNumber(icell);
      Float_t amp  = cells.GetAmplitude(icell);

      GetGeomInfo(id,iSupMod,x,y,z);

      fSM[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
      if(iSupMod<fNsmfull) fSMfull[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
      if(iSupMod>fNsmfull) fSMhalf[iSupMod-10]->RegisterDigit(id,iSupMod,amp,x,y,z);

    } // end loop cells
}

//______________________________________________________________________________
void AliEveEMCALData::LoadRecPoints(TTree* t)
{
  //*************************************************
  // To be improved !!!!!
  // Size and shape of cluster to be implemented
  // 
  //*************************************************
        
  // From TTreeR
  TObjArray *carr=NULL;
  TBranch *cbranch=t->GetBranch("EMCALECARP");
  cbranch->SetAddress(&carr);
  
  if(cbranch->GetEvent(0)) {
    for(Int_t ic = 0; ic < carr->GetEntriesFast(); ic++) {
      AliEMCALRecPoint* rp =(AliEMCALRecPoint*)carr->UncheckedAt(ic);
      if(rp){
	cout << "RecPoint info " << rp->GetAbsId() << " " << rp->GetEnergy() << endl;
	Int_t iSupMod = rp->GetSuperModuleNumber();
	Float_t amp = rp->GetEnergy();
	TVector3 lpos;
	rp->GetLocalPosition(lpos);

	fSM[iSupMod]->RegisterCluster(iSupMod,amp,lpos[0],lpos[1],lpos[2]);
      }
    }
  }
  
}

//______________________________________________________________________________
void AliEveEMCALData::LoadRecPointsFromEMCALLoader(AliEMCALLoader* emcl)
{
  //*************************************************
  // To be improved !!!!!
  // Size and shape of cluster to be implemented
  // 
  //*************************************************

  // From EMCALLoader
  AliEMCALRecPoint* rp;
  
  //Fill array of clusters                                                                        
  TClonesArray *clusters = (TClonesArray*)emcl->RecPoints();
  
  //Get clusters from the list                                                                    
  for(Int_t iclu = 0; iclu< clusters->GetEntries();iclu++){

    rp = static_cast<AliEMCALRecPoint *>(clusters->At(iclu)) ;
    
    if(rp){
       cout << "RecPoint info " << rp->GetAbsId() << " " << rp->GetEnergy() << endl;
       Int_t iSupMod = rp->GetSuperModuleNumber();
       Float_t amp = rp->GetEnergy();
       TVector3 lpos;
       rp->GetLocalPosition(lpos);
       
       fSM[iSupMod]->RegisterCluster(iSupMod,amp,lpos[0],lpos[1],lpos[2]);
    }
  }
  
}

//______________________________________________________________________________
void AliEveEMCALData::LoadRecPointsFromESD()
{
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
  Int_t iSM =  0 ;
  Int_t iT  =  0 ;
  Int_t iIp   =  0 ;
  Int_t iIe   =  0 ;
  Int_t ip    =  0 ;
  Int_t ie    =  0 ;
  Double_t xd, yd, zd;
  Float_t pos[3] ; 
  
  // Get reconstructed vertex position
  AliESDVertex* primVertex =(AliESDVertex*) fESD->GetVertex();
  Double_t vertex_position[3] ; 
  primVertex->GetXYZ(vertex_position) ; 

  //Get the CaloClusters
  //select EMCAL clusters only 
  TRefArray * caloClusters  = new TRefArray();
  fESD->GetEMCALClusters(caloClusters);
  Int_t nclus = caloClusters->GetEntries();
  cout << "nclus: " << nclus << endl; 
  
  if(!caloClusters) return;

  for (Int_t iclus =  0; iclus <  nclus; iclus++) 
    {
      AliESDCaloCluster *clus = (AliESDCaloCluster *) caloClusters->At(iclus) ; 
      //Get the cluster info
      
      Float_t energy = clus->E() ;
      Int_t   eneInt = (Int_t) energy*500+0.5;
      Float_t disp   = clus->GetClusterDisp() ;
      Int_t   iprim  = clus->GetLabel();
      
      clus->GetPosition(pos) ; // Global position
      TVector3 vpos(pos[0],pos[1],pos[2]) ;
      TLorentzVector p4 ;
      TVector3 p3;
      clus->GetMomentum(p4,vertex_position);
      p3.SetXYZ(p4[0],p4[1],p4[2]);
      Float_t eta = p3.Eta();
      Float_t phi = ( (p3.Phi()) < 0) ? (p3.Phi()) + 2. * TMath::Pi() : (p3.Phi());
      
      Int_t mult = clus->GetNCells() ;
      if(fDebug>1) {
	cout << "In cluster: " << iclus << ", ncells: " << mult << ", energy: " << 
	  eneInt << ", disp: " << disp << endl;
	cout << "Cluster " << iclus << ", eta: " << eta << ", phi: " << phi << endl;
      }

      Int_t clusId = 0;
      fGeom->GetAbsCellIdFromEtaPhi(eta,phi,clusId);
      if(fDebug>1) {
	cout << "Abs Cluster Id: " << clusId << ", xc: " << pos[0] << 
	  ", yc: " << pos[1] << ", zc: " << pos[2] << endl;
      }

      GetGeomInfo(clusId,iSupMod,x,y,z);
      
      //******** Not used yet but will come  ********
      Int_t     digMult = clus->GetNCells() ;
      UShort_t *digID   = clus->GetCellsAbsId() ;
      for(Int_t i=0; i<digMult; i++){
	fGeom->RelPosCellInSModule(digID[i], xd, yd, zd);
	//Geometry methods
	fGeom->GetCellIndex(digID[i],iSM,iT,iIp,iIe);
	//Gives SuperModule and Tower numbers

      } // end digit loop
      //*********************************************

      fSM[iSupMod]->RegisterCluster(iSM,energy,x,y,z);

    } // end cluster loop
}

//______________________________________________________________________________
AliEveEMCALSModuleData* AliEveEMCALData::GetSModuleData(Int_t sm)
{
  //
  // return super module data
  //

  if (sm < 0 || sm > fNsm) 
    {
      printf("The number of super modules must be lower or equal to %d",fNsm);
      return 0;
    }
  
  return fSM[sm];
}

//______________________________________________________________________________
void AliEveEMCALData::LoadRaw()
{

}

