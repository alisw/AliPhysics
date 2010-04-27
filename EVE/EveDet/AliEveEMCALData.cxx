//
// Fill containers for visualisation of EMCAL data structures
//    - read and store MC Hits    - read and store digits from esds or runloader
//    - read and store clusters from esds or runloader 
//
// Author: Magali Estienne (magali.estienne@cern.ch)
// June 30 2008
//

//#include <Riostream.h>
//#include <vector>

#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TRefArray.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "AliRunLoader.h"
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
  // delete fEmcal; // deleted by run-loader
  delete fGeom;
  delete fNode;
  delete fHMatrix;
  delete fPoint;
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
void AliEveEMCALData::SetTree(TTree* const tree)
{
  //
  // Set digit-tree to be used for digit retrieval. 
  // Data is loaded on demand.
  // 

  fTree = tree;

}

//______________________________________________________________________________
void AliEveEMCALData::SetESD(AliESDEvent* const esd)
{
  //
  // Set esd
  //

  fESD = esd;
}

//______________________________________________________________________________
void AliEveEMCALData::SetNode(TGeoNode* const node)
{
  //
  // Set node
  //

  fNode = node;
}

//______________________________________________________________________________
void AliEveEMCALData::InitEMCALGeom(AliRunLoader* const rl)
{
  //
  // Set data members for EMCAL geometry
  //

  fEmcal = (AliEMCAL*) rl->GetAliRun()->GetDetector("EMCAL");
  fGeom  = (AliEMCALGeometry*) fEmcal->GetGeometry();

}

//______________________________________________________________________________
void AliEveEMCALData::GetGeomInfo(Int_t id, Int_t &iSupMod, Double_t& x, Double_t& y, Double_t& z)
{
  //
  // Get geometrical information from hit/digit/cluster absolute id
  //

  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;

  //Geometry methods
  fGeom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
  //Gives SuperModule and Tower numbers
  fGeom->RelPosCellInSModule(id, x, y, z);

}

//______________________________________________________________________________
void  AliEveEMCALData::CreateAllSModules()
{
  //
  // Create all fNsm super modules
  //

  for(Int_t sm = 0; sm < fNsm; sm++)
    CreateSModule(sm);

}

//______________________________________________________________________________
void  AliEveEMCALData::CreateSModule(Int_t sm)
{
  //
  // Create super-module-data for SM if it does not exist already.
  //

  if(fSM[sm] == 0) fSM[sm] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
  if(fSMfull[sm] == 0 && sm < 10) fSMfull[sm] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
  if(fSMhalf[sm-10] == 0 && sm > 10) fSMhalf[sm-10] = new AliEveEMCALSModuleData(sm,fGeom,fNode,fHMatrix);
}

//______________________________________________________________________________
void AliEveEMCALData::DropAllSModules()
{
  //
  // Drop data of all existing sectors.
  //

  for (Int_t sm = 0; sm < fNsm; sm++) {
    if (fSM[sm] != 0)
      fSM[sm]->DropData();
  }
}

//______________________________________________________________________________
void AliEveEMCALData::DeleteSuperModules()
{
  //
  // Delete all super module data
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
void AliEveEMCALData::LoadHits(TTree* const t)
{
  //
  // Get hit information from RunLoader
  //

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
  
  if(hbranch->GetEvent(0)) {
    for(Int_t ih = 0; ih < harr->GetEntriesFast(); ih++) {
      AliEMCALHit* hit =(AliEMCALHit*)harr->UncheckedAt(ih);
      if(hit != 0){
	if(fDebug>1) cout << "Hit info " << hit->GetId() << " " << hit->GetEnergy() << endl;
	Int_t id = hit->GetId();
	// These are local coordinates
	Double_t xl = 0.; Double_t yl = 0.; Double_t zl = 0.;
	// Get global coordinates
	Double_t x = hit->X();
	Double_t y = hit->Y();
	Double_t z = hit->Z();
	Double_t amp = hit->GetEnergy();
	Int_t iSupMod = 0;
	// Get SM Id
	GetGeomInfo(id,iSupMod,xl,yl,zl);
	fSM[iSupMod]->RegisterHit(id,iSupMod,amp,x,y,z);
      }
    }
  }
}

//______________________________________________________________________________
void AliEveEMCALData::LoadHitsFromEMCALLoader(AliEMCALLoader* const emcl)
{
  //
  // Get hit information from EMCAL Loader
  //

  AliEMCALHit* hit;
  
  //Fill array of hits                                                                        
  TClonesArray *hits = (TClonesArray*)emcl->Hits();

  //Get hits from the list                                                                    
  for(Int_t ihit = 0; ihit< hits->GetEntries();ihit++){

    hit = static_cast<AliEMCALHit *>(hits->At(ihit)) ;
    
    if(hit != 0){
      if(fDebug>1) cout << "Hit info " << hit->GetId() << " " << hit->GetEnergy() << endl;

      Int_t id = hit->GetId();
      // These are local coordinates
      Double_t xl = 0.; Double_t yl = 0.; Double_t zl = 0.;
      // Get global coordinates
      Double_t x = hit->X();
      Double_t y = hit->Y();
      Double_t z = hit->Z();
      Double_t amp = hit->GetEnergy();
      Int_t iSupMod = 0;
      // Get SM Id
      GetGeomInfo(id,iSupMod,xl,yl,zl);
      fSM[iSupMod]->RegisterHit(id,iSupMod,amp,x,y,z);
    }
  }
  
}

//______________________________________________________________________________
void AliEveEMCALData::LoadDigits(TTree *t)
{
  //
  // Get digit information from RunLoader
  //

  TClonesArray *digits = 0;
  t->SetBranchAddress("EMCAL", &digits);
  t->GetEntry(0);

  Int_t nEnt = digits->GetEntriesFast();
  cout << "nEnt: " << nEnt << endl;
  AliEMCALDigit * dig;

  //  Double_t amp   = -1 ;
  Double_t ampFlo   = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Double_t x, y, z;

  for (Int_t idig = 0; idig < nEnt; idig++)
    {
      dig = static_cast<AliEMCALDigit *>(digits->At(idig));
      
      if(dig != 0) {
	id   = dig->GetId() ; //cell (digit) label
	// adc
	ampFlo  = dig->GetAmplitude(); //amplitude in cell (digit)
	// GeV
	// 	amp = ampFlo*0.0153; // To be modified with correct OCDB conversion	

	GetGeomInfo(id,iSupMod,x,y,z);

// 	// GeV
// 	fSM[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
// //	fSM[iSupMod]->SaveDigit(dig);
// // 	if(iSupMod<fNsmfull) fSMfull[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
// // 	if(iSupMod>fNsmfull) fSMhalf[iSupMod-10]->RegisterDigit(id,iSupMod,amp,x,y,z);
	fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
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
void AliEveEMCALData::LoadDigitsFromEMCALLoader(AliEMCALLoader* const emcl)
{

  //
  // Get digit information from EMCAL Loader
  //

  AliEMCALDigit* dig;
  
  //Fill array of digits                                                                        
  TClonesArray *digits = (TClonesArray*)emcl->Digits();
  
  //Get digits from the list  
  
  //  Double_t amp   = -1 ;
  Double_t ampFlo   = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Double_t x, y, z;
                                               
   for(Int_t idig = 0; idig< digits->GetEntries();idig++){

     dig = static_cast<AliEMCALDigit *>(digits->At(idig)) ;

     if(dig != 0){
       if(fDebug>1) cout << "Digit info " << dig->GetId() << " " << dig->GetAmp() << endl;
       id   = dig->GetId() ; //cell (digit) label
       // adc
       ampFlo  = dig->GetAmplitude(); //amplitude in cell (digit)
       // GeV
       //       amp = ampFlo*0.0153.; // To be modified with correct OCDB conversion

       GetGeomInfo(id,iSupMod,x,y,z);
       
       //       // GeV
       //       fSM[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
       // adc
       fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
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
  Int_t iSupMod =  0 ;
  Double_t x, y, z;

  // Extract digit information from the ESDs
  for (Int_t icell=  0; icell <  ncell; icell++) 
    {
      Int_t id   = cells.GetCellNumber(icell);
      // adc
      Double_t ampFlo  = cells.GetAmplitude(icell);
      // GeV
      //      Double_t amp = ampFlo*0.0153; // To be modified with correct OCDB conversion

      GetGeomInfo(id,iSupMod,x,y,z);

//       // GeV
//       fSM[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
//       if(iSupMod<fNsmfull) fSMfull[iSupMod]->RegisterDigit(id,iSupMod,amp,x,y,z);
//       if(iSupMod>fNsmfull) fSMhalf[iSupMod-10]->RegisterDigit(id,iSupMod,amp,x,y,z);
      // adc
      fSM[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
      if(iSupMod<fNsmfull) fSMfull[iSupMod]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);
      if(iSupMod>fNsmfull) fSMhalf[iSupMod-10]->RegisterDigit(id,iSupMod,ampFlo,x,y,z);

    } // end loop cells
}

//______________________________________________________________________________
void AliEveEMCALData::LoadRecPoints(TTree* const t)
{
  //
  // Get rec point information from RunLoader
  //

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
	if(fDebug>1) cout << "RecPoint info " << rp->GetAbsId() << " " << rp->GetEnergy() << endl;
	Int_t iSupMod = rp->GetSuperModuleNumber();
	// GeV
	Double_t amp = (Double_t)rp->GetEnergy();
	// adc
	Double_t ampFlo = amp/0.0153; // To be modified with correct OCDB conversion
	TVector3 lpos;
	rp->GetLocalPosition(lpos);

// 	// GeV
// 	fSM[iSupMod]->RegisterCluster(iSupMod,amp,lpos[0],lpos[1],lpos[2]);
        // adc
	fSM[iSupMod]->RegisterCluster(iSupMod,ampFlo,lpos[0],lpos[1],lpos[2]);
      }
    }
  }
  
}

//______________________________________________________________________________
void AliEveEMCALData::LoadRecPointsFromEMCALLoader(AliEMCALLoader* const emcl)
{
  //
  // Get rec point information from EMCAL Loader
  //

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
       if(fDebug>1) cout << "RecPoint info " << rp->GetAbsId() << " " << rp->GetEnergy() << endl;
       Int_t iSupMod = rp->GetSuperModuleNumber();
       Double_t amp = (Double_t)rp->GetEnergy();
       Double_t ampFlo = amp/0.0153; // To be modified with correct OCDB conversion
       TVector3 lpos;
       rp->GetLocalPosition(lpos);
       
//        // GeV
//        fSM[iSupMod]->RegisterCluster(iSupMod,amp,lpos[0],lpos[1],lpos[2]);
       // adc
       fSM[iSupMod]->RegisterCluster(iSupMod,ampFlo,lpos[0],lpos[1],lpos[2]);
    }
  }
  
}

//______________________________________________________________________________
void AliEveEMCALData::LoadRecPointsFromESD()
{
  //
  // Get cluster information from esd
  //

 Int_t iSupMod =  0 ;
  Double_t x, y, z;
  Int_t iSM =  0 ;
  Int_t iT  =  0 ;
  Int_t iIp   =  0 ;
  Int_t iIe   =  0 ;
  Double_t xd, yd, zd;
  Float_t pos[3] ; 
  
  // Get reconstructed vertex position
  AliESDVertex* primVertex =(AliESDVertex*) fESD->GetVertex();
  Double_t vertexPosition[3] ; 
  primVertex->GetXYZ(vertexPosition) ; 

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
      
      Double_t energy = clus->E() ;
      // adc
      //      Int_t   eneInt = (Int_t)energy*500+0.5;
      Double_t eneInt = energy/0.0153; // To be modified with correct OCDB conversion
      Double_t disp   = clus->GetClusterDisp() ;
      
      clus->GetPosition(pos) ; // Global position
      TVector3 vpos(pos[0],pos[1],pos[2]) ;
      TLorentzVector p4 ;
      TVector3 p3;
      clus->GetMomentum(p4,vertexPosition);
      p3.SetXYZ(p4[0],p4[1],p4[2]);
      Double_t eta = p3.Eta();
      Double_t phi = ( (p3.Phi()) < 0) ? (p3.Phi()) + 2. * TMath::Pi() : (p3.Phi());
      
      Int_t mult = clus->GetNCells() ;
      if(fDebug>2) {
	cout << "In cluster: " << iclus << ", ncells: " << mult << ", energy : " << energy << 
	  ", disp: " << disp << endl;
	cout << "Cluster " << iclus << ", eta: " << eta << ", phi: " << phi << endl;
      }

      Int_t clusId = 0;
      fGeom->GetAbsCellIdFromEtaPhi(eta,phi,clusId);
      if(fDebug>2) {
	cout << "Abs Cluster Id: " << clusId << ", xc: " << pos[0] << 
	  ", yc: " << pos[1] << ", zc: " << pos[2] << endl;
      }

      GetGeomInfo(clusId,iSupMod,x,y,z);
      
      //******** Not used yet but will come  ********
      // AliESDCaloCells &cells= *(fESD->GetEMCALCells());
      Int_t     digMult = clus->GetNCells() ;
      UShort_t *digID   = clus->GetCellsAbsId() ;
      for(Int_t i=0; i<digMult; i++){
	// Float_t  digitAmp     = cells.GetCellAmplitude(digID[i]) ;
	fGeom->RelPosCellInSModule(digID[i], xd, yd, zd);
	//Geometry methods
	fGeom->GetCellIndex(digID[i],iSM,iT,iIp,iIe);
	//Gives SuperModule and Tower numbers

      } // end digit loop
      //*********************************************
      //      // GeV
      //      fSM[iSupMod]->RegisterCluster(iSM,energy,x,y,z);
      // adc
      fSM[iSupMod]->RegisterCluster(iSM,eneInt,x,y,z);

    } // end cluster loop
}

//______________________________________________________________________________
AliEveEMCALSModuleData* AliEveEMCALData::GetSModuleData(Int_t sm)
{
  //
  // Return super module data
  //

  if (sm < 0 || sm > fNsm) 
    {
      printf("The number of super modules must be lower or equal to %d",fNsm);
      return 0;
    }
  
  return fSM[sm];
}

//______________________________________________________________________________
void AliEveEMCALData::LoadRaw() const
{
  //
  // Get raw information
  //

  // To be implemented !
}

