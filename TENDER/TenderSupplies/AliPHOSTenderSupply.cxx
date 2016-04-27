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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  PHOS tender, recalibrate PHOS clusters                                   //
//  and do track matching                                                    //
//  Author : Dmitri Peressounko (RRC KI)                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom.h"

#include <AliLog.h>
#include <AliVEvent.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliTender.h>
#include <AliCDBManager.h>
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"

#include "AliVCluster.h"
#include "AliPHOSTenderSupply.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSAodCluster.h"
#include "AliOADBContainer.h"
#include "AliAODCaloCells.h"
#include "AliESDCaloCells.h"

ClassImp(AliPHOSTenderSupply)

AliPHOSTenderSupply::AliPHOSTenderSupply() :
  AliTenderSupply()
  ,fOCDBpass("local://OCDB")
  ,fNonlinearityVersion("Default")
  ,fPHOSGeo(0x0)
  ,fRecoPass(-1)  //to be defined
  ,fUsePrivateBadMap(0)
  ,fUsePrivateCalib(0)
  ,fAddNoiseMC(0)
  ,fNoiseMC(0.001)
  ,fPHOSCalibData(0x0)
  ,fTask(0x0)
  ,fIsMC(kFALSE)
  ,fMCProduction("")  
{
	//
	// default ctor
	//
   for(Int_t i=0;i<10;i++)fNonlinearityParams[i]=0. ;
   for(Int_t mod=0;mod<6;mod++)fPHOSBadMap[mod]=0x0 ;
}

//_____________________________________________________
AliPHOSTenderSupply::AliPHOSTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender)
  ,fOCDBpass("alien:///alice/cern.ch/user/p/prsnko/PHOSrecalibrations/")
  ,fNonlinearityVersion("Default")
  ,fPHOSGeo(0x0)
  ,fRecoPass(-1)  //to be defined
  ,fRunNumber(-1) //to be defined
  ,fUsePrivateBadMap(0)
  ,fUsePrivateCalib(0)
  ,fAddNoiseMC(0)
  ,fNoiseMC(0.001)
  ,fPHOSCalibData(0x0)
  ,fTask(0x0)
  ,fIsMC(kFALSE)
  ,fMCProduction("")  
{
	//
	// named ctor
	//
   for(Int_t i=0;i<10;i++)fNonlinearityParams[i]=0. ;
   for(Int_t mod=0;mod<6;mod++)fPHOSBadMap[mod]=0x0 ;
}

//_____________________________________________________
AliPHOSTenderSupply::~AliPHOSTenderSupply()
{
  //Destructor
  if(fPHOSCalibData)
    delete fPHOSCalibData;
  fPHOSCalibData=0x0 ;
}

//_____________________________________________________
void AliPHOSTenderSupply::InitTender()
{   
  //
  // Initialise PHOS tender
  //
  AliESDEvent *esd = 0x0 ; 
  AliAODEvent *aod = 0x0 ;
  if(fTender){
    esd = fTender->GetEvent();
  }
  else{
    if(fTask){
      esd = dynamic_cast<AliESDEvent*>(fTask->InputEvent()) ;
      aod = dynamic_cast<AliAODEvent*>(fTask->InputEvent()) ;
    }
  }     
  
  if(fTender)
    fRunNumber = fTender->GetRun();
  else{
    if(!fTask){
      AliError("Neither Tender not Taks was not set") ;
      return ;
    }
    if(aod)
      fRunNumber = aod->GetRunNumber() ;
    else{
      if(esd)
        fRunNumber = esd->GetRunNumber() ;
      else{
        AliError("Taks does not contain neither ESD nor AOD") ;
        return ;
      }
    }   
  }

  //In MC always reco pass 1
  if(fIsMC)
    fRecoPass=1 ;
  
  if(fRecoPass<0){ //not defined yet
    // read if from filename.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    TTree * t = mgr->GetTree();
    if(t){  
      TFile * f = t->GetCurrentFile() ;
      if(f){  
        TString fname(f->GetName());
        if(fname.Contains("pass1"))
	   fRecoPass=1;
        else 
	  if(fname.Contains("pass2"))
	   fRecoPass=2;
          else 
	    if(fname.Contains("pass3")) 
  	      fRecoPass=3;
            else 
	      if(fname.Contains("pass4")) 
  	        fRecoPass=4;
      }
    }
    if(fRecoPass<0){
      AliError("Can not find pass number from file name, set it manually");
    }
  }
   
  //Init geometry 
  if(!fPHOSGeo){
    if(fRunNumber<209122) //Run1
      fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
    else
      fPHOSGeo =  AliPHOSGeometry::GetInstance("Run2") ;      
    AliOADBContainer geomContainer("phosGeo");
    if(fIsMC){ //use excatly the same geometry as in simulation, stored in esd
      if(esd){
        for(Int_t mod=0; mod<6; mod++) {
           const TGeoHMatrix * m = esd->GetPHOSMatrix(mod) ;
           if(m){
             fPHOSGeo->SetMisalMatrix(m,mod) ;
             printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
             m->Print() ;
	   }
	}
      } 
      if(aod){ //To be fixed
        geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCGeometry.root","PHOSMCRotationMatrixes");
        TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
        for(Int_t mod=0; mod<6; mod++) {
          if(!matrixes->At(mod)) continue;
          fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
          printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
          ((TGeoHMatrix*)matrixes->At(mod))->Print() ;
        } 
      }
    }
    else{ //Use best approaximation to real geometry
      geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
      TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
      for(Int_t mod=0; mod<6; mod++) {
        if(!matrixes->At(mod)) continue;
        fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
        printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
        ((TGeoHMatrix*)matrixes->At(mod))->Print() ;
      } 
    }
  }
  
  //Init Bad channels map
  if(!fUsePrivateBadMap){
   AliOADBContainer badmapContainer(Form("phosBadMap"));
    badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root","phosBadMap");
    TObjArray *maps = (TObjArray*)badmapContainer.GetObject(fRunNumber,"phosBadMap");
    if(!maps){
      AliError(Form("Can not read Bad map for run %d. \n You may choose to use your map with ForceUsingBadMap()\n",fRunNumber)) ;    
    }
    else{
      AliInfo(Form("Setting PHOS bad map with name %s \n",maps->GetName())) ;
      for(Int_t mod=0; mod<6;mod++){
        if(fPHOSBadMap[mod]) 
          delete fPHOSBadMap[mod] ;
        TH2I * h = (TH2I*)maps->At(mod) ;      
	if(h)
          fPHOSBadMap[mod]=new TH2I(*h) ;
      }
    }    
  } 

  if(!fUsePrivateCalib){
    if(fIsMC){ //re/de-calibration for MC productions
      //Init recalibration
      AliOADBContainer calibContainer("phosRecalibration");
      calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCCalibrations.root","phosRecalibration");

      AliInfo(Form("Reading PHOS MC recalibration object for production %s, run=%d", fMCProduction.Data(),fRunNumber)) ;      
      TObjArray *recalib = (TObjArray*)calibContainer.GetObject(fRunNumber,"PHOSRecalibration",fMCProduction.Data());
      if(!recalib){
        AliFatal(Form("Can not read calibrations for run %d and name >%s<\n. You may choose your specific calibration with ForceUsingCalibration()\n",fRunNumber,fMCProduction.Data())) ;
      }
      else{
	//Now try to find object with proper name
	for(Int_t i=0; i<recalib->GetEntriesFast(); i++){
	  AliPHOSCalibData * tmp = (AliPHOSCalibData*)recalib->At(i) ;
	  if(fMCProduction.CompareTo(tmp->GetName())==0){
            fPHOSCalibData = tmp ;
	    break ;
	  }
	}
        if(!fPHOSCalibData) {
          AliFatal(Form("Can not find calibration for run %d, and name %s \n",fRunNumber, fMCProduction.Data())) ;
        }
      }
      
    }
    else{ //real data
      //Init recalibration
      //Check the pass1-pass2-pass3 reconstruction
      AliOADBContainer calibContainer("phosRecalibration");
      calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");
      TObjArray *recalib = (TObjArray*)calibContainer.GetObject(fRunNumber,"PHOSRecalibration");
      if(!recalib){
        AliFatal(Form("Can not read calibrations for run %d\n. You may choose your specific calibration with ForceUsingCalibration()\n",fRunNumber)) ;
      }
      else{
        fPHOSCalibData = (AliPHOSCalibData*)recalib->At(fRecoPass-1) ;
        if(!fPHOSCalibData) {
          AliFatal(Form("Can not find calibration for run %d, pass %d \n",fRunNumber, fRecoPass)) ;
        }
      }
    }
  }
  
}

//_____________________________________________________
void AliPHOSTenderSupply::ProcessEvent()
{
  //Choose PHOS clusters and recalibrate them
  //that it recalculate energy, position and distance 
  //to closest track extrapolation	

  AliESDEvent *esd = 0x0 ; 
  AliAODEvent *aod = 0x0 ;
  if(fTender){
    esd = fTender->GetEvent();
    if(!esd)
      return ;
  }
  else{
    if(!fTask){
      return ;
    }
    esd = dynamic_cast<AliESDEvent*>(fTask->InputEvent()) ;
    aod = dynamic_cast<AliAODEvent*>(fTask->InputEvent()) ;
    if(!esd && !aod)
      return ;
  }     
    
  if(!fPHOSCalibData 
    || (fTender && fTender->RunChanged())){ //In case of Task init called automatically
    InitTender();
    
  }

  TVector3 vertex ;
  if(esd){
    const AliESDVertex *esdVertex = esd->GetPrimaryVertex();
    vertex.SetXYZ(esdVertex->GetX(),esdVertex->GetY(),esdVertex->GetZ());
  }
  else{//AOD
    const AliAODVertex *aodVertex = aod->GetPrimaryVertex();
    vertex.SetXYZ(aodVertex->GetX(),aodVertex->GetY(),aodVertex->GetZ());
  }
  if(vertex.Mag()>99.) //vertex not defined?
    vertex.SetXYZ(0.,0.,0.) ;


  //For re-calibration
  const Double_t logWeight=4.5 ;  
  
  Int_t phosCluSelection =  AliVCluster::kPHOSNeutral ;  
//   if(fRunNumber>=209122) //Run2
//      phosCluSelection
//      
//   (fClusterType == kPHOSNeutral || fClusterType == kPHOSCharged)     

  if(esd){ //To avoid multiple if in loops we made 
           //almost identical pecies of code. Please apply changes to both!!!
    Int_t multClust=esd->GetNumberOfCaloClusters();
    AliESDCaloCells * cells = esd->GetPHOSCells() ;
    
    //Make copy of phos Energy and add noise
    if(fIsMC && fAddNoiseMC){
      Short_t ncell= cells->GetNumberOfCells() ;
      Short_t cellNumber;
      Double_t amplitude=0., time=0., efrac=0.;
      Int_t mclabel;
      for(Short_t pos=0; pos<ncell; pos++){
	 cells->GetCell(pos, cellNumber, amplitude, time, mclabel, efrac) ;
	 amplitude=TMath::Max(0.,amplitude+gRandom->Gaus(0,fNoiseMC)) ;
	 Bool_t isHG=cells->GetHighGain(pos) ;
         cells->SetCell(pos, cellNumber, amplitude, time,  mclabel,  efrac, isHG);
      }      
    }
 
    for (Int_t i=0; i<multClust; i++) {
      AliESDCaloCluster *clu = esd->GetCaloCluster(i);    
      if (clu->GetType()!=phosCluSelection) continue;
      
      //remove clusters from bad modules without re-calibration
      Int_t relid[4] ;
      fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(0), relid) ; //shold be at lease one digit
      if(!fPHOSBadMap[relid[0]]){
        clu->SetE(0) ;
	continue ;
      }
     
      //Apply re-Calibreation
      AliPHOSEsdCluster cluPHOS(*clu);
      cluPHOS.Recalibrate(fPHOSCalibData,cells); // modify the cell energies
      cluPHOS.EvalAll(logWeight,vertex);         // recalculate the cluster parameters
      cluPHOS.SetE(CorrectNonlinearity(cluPHOS.E()));// Users's nonlinearity

      Float_t  position[3];
      cluPHOS.GetPosition(position);
      clu->SetPosition(position);                       //rec.point position in MARS      
      TVector3 global(position) ;
      Int_t relId[4] ;
      fPHOSGeo->GlobalPos2RelId(global,relId) ;
      Int_t mod  = relId[0] ;
      Int_t cellX = relId[2];
      Int_t cellZ = relId[3] ;
      if ( !IsGoodChannel(mod,cellX,cellZ) ) {
        clu->SetE(0.) ;
        continue ;
      }  
            
      Double_t ecore=CoreEnergy(&cluPHOS) ; 
      ecore=CorrectNonlinearity(ecore) ;
      
      clu->SetE(cluPHOS.E());                      //total particle energy
      clu->SetCoreEnergy(ecore);                            //core particle energy
           
      //Eval FullDispersion
      clu->SetDispersion(TestFullLambda(clu->E(),cluPHOS.GetM20(),cluPHOS.GetM02())) ;
      //Eval CoreDispersion
      Double_t m02=0.,m20=0.;
      EvalLambdas(&cluPHOS,m02, m20);   
      clu->SetChi2(TestCoreLambda(clu->E(),m20,m02));                     //not yet implemented     
      clu->SetM02(m02) ;               //second moment M2x
      clu->SetM20(m20) ;               //second moment M2z
      
      //correct distance to track      
      Double_t dx=0.,dz=0. ;
      fPHOSGeo->GlobalPos2RelId(global,relId) ;
      TVector3 locPos;
      fPHOSGeo->Global2Local(locPos,global,mod) ;

      Double_t pttrack=0.;
      Int_t charge=0;
      FindTrackMatching(mod,&locPos,dx,dz,pttrack,charge) ;
      Double_t r=TestCPV(dx, dz, pttrack,charge) ;
      clu->SetTrackDistance(dx,dz); 
     
      clu->SetEmcCpvDistance(r);    
      
      Double_t tof=EvalTOF(&cluPHOS,cells); 
//      if(TMath::Abs(tof-clu->GetTOF())>100.e-9) //something wrong in cell TOF!
//	tof=clu->GetTOF() ;
      clu->SetTOF(tof);       
      Double_t minDist=clu->GetDistanceToBadChannel() ;//Already calculated
      DistanceToBadChannel(mod,&locPos,minDist);
      clu->SetDistanceToBadChannel(minDist) ;

      Double_t ecross = EvalEcross(&cluPHOS);  
      clu->SetMCEnergyFraction(ecross) ;
    }
    
  }
  else{//AOD
    Int_t multClust=aod->GetNumberOfCaloClusters();
    AliAODCaloCells * cells = aod->GetPHOSCells() ;
    //Add noise
    if(fIsMC && fAddNoiseMC){
      Short_t ncell= cells->GetNumberOfCells() ;
      Short_t cellNumber;
      Double_t amplitude=0., time=0., efrac=0.;
      Int_t mclabel;
      for(Short_t pos=0; pos<ncell; pos++){
	 cells->GetCell(pos, cellNumber, amplitude, time, mclabel, efrac) ;
	 amplitude=TMath::Max(0.,amplitude+gRandom->Gaus(0,fNoiseMC)) ;
	 Bool_t isHG=cells->GetHighGain(pos) ;
         cells->SetCell(pos, cellNumber, amplitude, time,  mclabel,  efrac, isHG);
      }      
    }
  
    for (Int_t i=0; i<multClust; i++) {
      AliAODCaloCluster *clu = aod->GetCaloCluster(i);    
//      if ( !clu->IsPHOS()) continue;
      if (clu->GetType()!=phosCluSelection) continue;
      
    
      //Apply re-Calibreation
      AliPHOSAodCluster cluPHOS(*clu);
      cluPHOS.Recalibrate(fPHOSCalibData,cells); // modify the cell energies
      cluPHOS.EvalAll(logWeight,vertex);         // recalculate the cluster parameters
      cluPHOS.SetE(CorrectNonlinearity(cluPHOS.E()));// Users's nonlinearity

      Float_t  position[3];
      cluPHOS.GetPosition(position);
      clu->SetPosition(position);                       //rec.point position in MARS
      TVector3 global(position) ;
      Int_t relId[4] ;
      fPHOSGeo->GlobalPos2RelId(global,relId) ;
      Int_t mod  = relId[0] ;
      Int_t cellX = relId[2];
      Int_t cellZ = relId[3] ;
      if ( !IsGoodChannel(mod,cellX,cellZ) ) {
        clu->SetE(0.) ;
        continue ;
      }  
      TVector3 locPosOld; //Use it to re-calculate distance to track
      fPHOSGeo->Global2Local(locPosOld,global,mod) ;
      
      Double_t ecore=CoreEnergy(&cluPHOS) ; 
      ecore=CorrectNonlinearity(ecore) ;
     
      clu->SetE(cluPHOS.E());                           //total particle energy
      clu->SetCoreEnergy(ecore);                  //core particle energy

      //Eval FullDispersion
      clu->SetDispersion(TestFullLambda(clu->E(),cluPHOS.GetM20(),cluPHOS.GetM02())) ;
      //Eval CoreDispersion
      Double_t m02=0.,m20=0.;
      EvalLambdas(&cluPHOS,m02, m20);   
      clu->SetChi2(TestCoreLambda(clu->E(),m20,m02));                     //not yet implemented
      clu->SetM02(m02) ;               //second moment M2x
      clu->SetM20(m20) ;               //second moment M2z
      
      //correct distance to track
      Double_t dx=clu->GetTrackDx() ;
      Double_t dz=clu->GetTrackDz() ;
      TVector3 locPos;
      fPHOSGeo->Global2Local(locPos,global,mod) ;
      if(dx!=-999.){ //there is matched track
        dx+=locPos.X()-locPosOld.X() ;
        dz+=locPos.Z()-locPosOld.Z() ;      
        clu->SetTrackDistance(dx,dz);
      }
      Double_t r = 999. ; //Big distance
      int nTracksMatched = clu->GetNTracksMatched();
      if(nTracksMatched > 0) {
        AliVTrack* track = dynamic_cast<AliVTrack*> (clu->GetTrackMatched(0));
        if ( track ) {
          Double_t pttrack = track->Pt();
          Short_t charge = track->Charge();
          r=TestCPV(dx, dz, pttrack,charge) ;
	}
      }
      clu->SetEmcCpvDistance(r); //Distance in sigmas


      Double_t tof=EvalTOF(&cluPHOS,cells); 
//      if(TMath::Abs(tof-clu->GetTOF())>100.e-9) //something wrong in cell TOF!
//	tof=clu->GetTOF() ;
      clu->SetTOF(tof);       
      Double_t minDist=clu->GetDistanceToBadChannel() ;//Already calculated
      DistanceToBadChannel(mod,&locPos,minDist);
      clu->SetDistanceToBadChannel(minDist) ;

      Double_t ecross = EvalEcross(&cluPHOS);  
      clu->SetMCEnergyFraction(ecross) ;      
    }
    
  }

}
//___________________________________________________________________________________________________
void AliPHOSTenderSupply::FindTrackMatching(Int_t mod,TVector3 *locpos,
					    Double_t &dx, Double_t &dz,
					    Double_t &pt,Int_t &charge){
  //Find track with closest extrapolation to cluster
  AliESDEvent *esd = 0x0 ;
  if(fTender)
    esd= fTender->GetEvent();
  else{ 
    esd= dynamic_cast<AliESDEvent*>(fTask->InputEvent());
  }
  
  if(!esd){
    AliError("ESD is not found") ;
    return ;
  }
  Double_t  magF = esd->GetMagneticField();
 
  Double_t magSign = 1.0;
  if(magF<0)magSign = -1.0;
  
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliError("Margnetic filed was not initialized, use default") ;
    AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  Int_t nt=0;
  nt = esd->GetNumberOfTracks();
  //Calculate actual distance to PHOS module
  TVector3 globaPos ;
  fPHOSGeo->Local2Global(mod, 0.,0., globaPos) ;
  const Double_t rPHOS = globaPos.Pt() ; //Distance to center of  PHOS module
  const Double_t kYmax = 72.+10. ; //Size of the module (with some reserve) in phi direction
  const Double_t kZmax = 64.+10. ; //Size of the module (with some reserve) in z direction
  const Double_t kAlpha0=330./180.*TMath::Pi() ; //First PHOS module angular direction
  const Double_t kAlpha= 20./180.*TMath::Pi() ; //PHOS module angular size
  Double_t minDistance = 1.e6;


  Double_t gposTrack[3] ; 

  Double_t bz = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->SolenoidField();
  bz = TMath::Sign(0.5*kAlmost0Field,bz) + bz;

  Double_t b[3]; 

    for (Int_t i=0; i<nt; i++) {
      AliESDtrack *esdTrack=esd->GetTrack(i);

      // Skip the tracks having "wrong" status (has to be checked/tuned)
      ULong_t status = esdTrack->GetStatus();
      if ((status & AliESDtrack::kTPCout)   == 0) continue;
     
      //Continue extrapolation from TPC outer surface
      const AliExternalTrackParam *outerParam=esdTrack->GetOuterParam();
      if (!outerParam) continue;
      AliExternalTrackParam t(*outerParam);
     
      t.GetBxByBz(b) ;
      //Direction to the current PHOS module
      Double_t phiMod=kAlpha0-kAlpha*mod ;
      if(!t.Rotate(phiMod))
        continue ;
    
      Double_t y;                       // Some tracks do not reach the PHOS
      if (!t.GetYAt(rPHOS,bz,y)) continue; //    because of the bending
      
      Double_t z; 
      if(!t.GetZAt(rPHOS,bz,z))
        continue ;
      if (TMath::Abs(z) > kZmax) 
        continue; // Some tracks miss the PHOS in Z
      if(TMath::Abs(y) < kYmax){
        t.PropagateToBxByBz(rPHOS,b);        // Propagate to the matching module
        //t.CorrectForMaterial(...); // Correct for the TOF material, if needed
        t.GetXYZ(gposTrack) ;
        TVector3 globalPositionTr(gposTrack) ;
        TVector3 localPositionTr ;
        fPHOSGeo->Global2Local(localPositionTr,globalPositionTr,mod) ;
        Double_t ddx = locpos->X()-localPositionTr.X();
        Double_t ddz = locpos->Z()-localPositionTr.Z();
        Double_t d2 = ddx*ddx + ddz*ddz;
        if(d2 < minDistance) {
	  dx = ddx ;
  	  dz = ddz ;
	  minDistance=d2 ;
	  pt=esdTrack->Pt() ;
	  charge=esdTrack->Charge() ;
        }
      }
    }//Scanned all tracks
 
}
//____________________________________________________________
Float_t AliPHOSTenderSupply::CorrectNonlinearity(Float_t en){

  //For backward compatibility, if no RecoParameters found
  if(fNonlinearityVersion=="Default"){
    return 0.0241+1.0504*en+0.000249*en*en ;
  }
  if(fNonlinearityVersion=="MC"){ //Default + some correction
    return (0.0241+1.0504*en+0.000249*en*en)*fNonlinearityParams[0]*(1+fNonlinearityParams[1]/(1.+en*en/fNonlinearityParams[2]/fNonlinearityParams[2])) ;
  }

  if(fNonlinearityVersion=="NoCorrection"){
    return en ;
  }
  if(fNonlinearityVersion=="Gustavo2005"){
    return fNonlinearityParams[0]+fNonlinearityParams[1]*en + fNonlinearityParams[2]*en*en ;
  }
  if(fNonlinearityVersion=="Henrik2010"){
    return en*(fNonlinearityParams[0]+fNonlinearityParams[1]*TMath::Exp(-en*fNonlinearityParams[2]))*(1.+fNonlinearityParams[3]*TMath::Exp(-en*fNonlinearityParams[4]))*(1.+fNonlinearityParams[6]/(en*en+fNonlinearityParams[5])) ;
  }

  return en ;
}
//_____________________________________________________________________________
Double_t AliPHOSTenderSupply::TestCoreLambda(Double_t pt,Double_t l1,Double_t l2){
//Parameterization for core dispersion   
//For R=4.5
  Double_t   l1Mean  = 1.150200 + 0.097886/(1.+1.486645*pt+0.000038*pt*pt) ;
  Double_t   l2Mean = 1.574706 + 0.997966*exp(-0.895075*pt)-0.010666*pt ;
  Double_t   l1Sigma = 0.100255 + 0.337177*exp(-0.517684*pt)+0.001170*pt ;
  Double_t   l2Sigma = 0.232580 + 0.573401*exp(-0.735903*pt)-0.002325*pt ;
  Double_t   c = -0.110983 -0.017353/(1.-1.836995*pt+0.934517*pt*pt) ;

/*
  //Parameterizatino for full dispersion
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
*/
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return R2 ;
  
}
//_____________________________________________________________________________
Double_t AliPHOSTenderSupply::TestFullLambda(Double_t pt,Double_t l1,Double_t l2){
//Parameterization for full dispersion   
  //Parameterizatino for full dispersion
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;

  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return R2 ;
  
}
//____________________________________________________________________________
Double_t AliPHOSTenderSupply::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC10h period
  //_true if neutral_
  
  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
              6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
  
  Double_t mf = 0.; //Positive for ++ and negative for --
  if(fTender){
    AliESDEvent *esd = fTender->GetEvent();
    mf = esd->GetMagneticField();
  }
  else{ 
    if(fTask){
      AliESDEvent *esd= dynamic_cast<AliESDEvent*>(fTask->InputEvent());
      if(esd)
         mf = esd->GetMagneticField();
      else{
        AliAODEvent *aod= dynamic_cast<AliAODEvent*>(fTask->InputEvent());
	if(aod)
          mf = aod->GetMagneticField();
      }
    }else{
       AliError("Neither Tender nor Task defined") ;    
    }
  }
  
  if(mf<0.){ //field --
    meanZ = -0.468318 ;
    if(charge>0)
      meanX=TMath::Min(7.3, 3.89994*1.20679*1.20679/(pt*pt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01)) ;
    else
      meanX=-TMath::Min(7.7,3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-pt*2.57070e+01)) ;
  }
  else{ //Field ++
    meanZ= -0.468318;
    if(charge>0)
      meanX=-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01)) ;
    else
      meanX= TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01)) ;     
  }

  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz) ;
}

//________________________________________________________________________
Bool_t AliPHOSTenderSupply::IsGoodChannel(Int_t mod, Int_t ix, Int_t iz)
{
  //Check if this channel belogs to the good ones
  
  if(!fPHOSBadMap[mod]){
     AliError(Form("No Bad map for PHOS module %d",mod)) ;
     return kFALSE ;
  }
  if(fPHOSBadMap[mod]->GetBinContent(ix,iz)>0)
    return kFALSE ;
  else
    return kTRUE ;
}
//________________________________________________________________________
void AliPHOSTenderSupply::ForceUsingBadMap(const char * filename){
  //Read TH2I histograms with bad maps from local or alien file 
  TFile * fbm = TFile::Open(filename) ;
  if(!fbm || !fbm->IsOpen()){
    AliError(Form("Can not open BadMaps file %s",filename)) ;
    return ;
  }
  gROOT->cd() ;
  char key[55] ;
  for(Int_t mod=1;mod<6; mod++){
    snprintf(key,55,"PHOS_BadMap_mod%d",mod) ;
    TH2I * h = (TH2I*)fbm->Get(key) ;
    if(h)
       fPHOSBadMap[mod] = new TH2I(*h) ;
  }
  fbm->Close() ;
  fUsePrivateBadMap=kTRUE ;
}
//________________________________________________________________________
void AliPHOSTenderSupply::ForceUsingCalibration(const char * filename){
  //Read PHOS recalibration parameters from the file.
  //We assume that file contains single entry: AliPHOSCalibData
  TFile * fc = TFile::Open(filename) ;
  if(!fc || !fc->IsOpen()){
    AliFatal(Form("Can not open Calibration file %s",filename)) ;
    return ;
  }
  fPHOSCalibData = (AliPHOSCalibData*)fc->Get("PHOSCalibration") ;
  fc->Close() ;
  fUsePrivateCalib=kTRUE; 
}
//________________________________________________________________________
void AliPHOSTenderSupply::CorrectPHOSMisalignment(TVector3 &global,Int_t mod){
   //Correct for PHOS modules misalignment 
  
    //correct misalignment
    const Float_t shiftX[6]={0.,-2.3,-2.11,-1.53,0.,0.} ;
    const Float_t shiftZ[6]={0.,-0.4, 0.52, 0.8,0.,0.} ;
    TVector3 localPos ;
    fPHOSGeo->Global2Local(localPos,global,mod) ;
    fPHOSGeo->Local2Global(mod,localPos.X()+shiftX[mod],localPos.Z()+shiftZ[mod],global);  
}
//________________________________________________________________________
void AliPHOSTenderSupply::EvalLambdas(AliVCluster * clu, Double_t &m02, Double_t &m20){ 
  //calculate dispecrsion of the cluster in the circle with radius distanceCut around the maximum
    
  const Double_t rCut=4.5 ;
  
  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  Double_t xc[100]={0} ;
  Double_t zc[100]={0} ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  Int_t mulDigit=TMath::Min(100,clu->GetNCells()) ;
  const Double_t logWeight=4.5 ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    if (clu->E()>0 && elist[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( elist[iDigit] / clu->E() ) ) ;
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
     
  wtot = 0. ;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
  Double_t xCut = 0. ;
  Double_t zCut = 0. ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    if (clu->E()>0 && elist[iDigit]>0.) {
        Double_t w = TMath::Max( 0., logWeight + TMath::Log( elist[iDigit] / clu->E() ) ) ;
        Double_t xi= xc[iDigit] ;
        Double_t zi= zc[iDigit] ;
	if((xi-x)*(xi-x)+(zi-z)*(zi-z) < rCut*rCut){
          xCut += w * xi ;
          zCut += w * zi ; 
          dxx  += w * xi * xi ;
          dzz  += w * zi * zi ;
          dxz  += w * xi * zi ; 
          wtot += w ;
	}
    }
    
  }
  if (wtot>0) {
    xCut/= wtot ;
    zCut/= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dxx -= xCut * xCut ;
    dzz -= zCut * zCut ;
    dxz -= xCut * zCut ;

    m02 =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    m20 =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
  }
  else {
    m20=m02=0.;
  }

}
//____________________________________________________________________________
Double_t  AliPHOSTenderSupply::CoreEnergy(AliVCluster * clu){  
  //calculate energy of the cluster in the circle with radius distanceCut around the maximum
  
  //Can not use already calculated coordinates?
  //They have incidence correction...
  const Double_t distanceCut =3.5 ;
  const Double_t logWeight=4.5 ;
  
  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  Double_t xc[100]={0} ;
  Double_t zc[100]={0} ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  Int_t mulDigit=TMath::Min(100,clu->GetNCells()) ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    if (clu->E()>0 && elist[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( elist[iDigit] / clu->E() ) ) ;
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
  Double_t coreE=0. ;
  for(Int_t iDigit=0; iDigit < mulDigit; iDigit++) {
    Double_t distance = TMath::Sqrt((xc[iDigit]-x)*(xc[iDigit]-x)+(zc[iDigit]-z)*(zc[iDigit]-z)) ;
    if(distance < distanceCut)
      coreE += elist[iDigit] ;
  }
  //Apply non-linearity correction
  return coreE ;
}
//____________________________________________________________________________
Double_t AliPHOSTenderSupply::EvalEcross(AliVCluster * clu){  
  //Calculate propoerion of the cluster energy in cross around the 
  //cell with maximal energy deposition. Can be used to reject exotic clusters 
  
  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
  Int_t mulDigit=clu->GetNCells() ;
  // Calculates the center of gravity in the local PHOS-module coordinates
  //Find cell with max E
  Double_t eMax=0.;
  Int_t iMax=0;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    if(elist[iDigit]>eMax){
      eMax=elist[iDigit] ;
      iMax=iDigit ;
    }
  }
  //Calculate e in cross
  Double_t eCross=0 ;
  Int_t relidMax[4] ;
  fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iMax), relidMax) ;
  Int_t relid[4] ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {    
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
    if(TMath::Abs(relid[2]-relidMax[2])+TMath::Abs(relid[3]-relidMax[3])==1)
      eCross+= elist[iDigit] ; 
  }
  if(eMax>0)
    return 1.-eCross/eMax ;
  else
    return 0 ;
}


//________________________________________________________________________
Double_t AliPHOSTenderSupply::EvalTOF(AliVCluster * clu,AliVCaloCells * cells){ 
  //Evaluate TOF of the cluster after re-calibration
  //TOF here is weighted average of digits
  // -within 50ns from the most energetic cell
  // -not too soft.
    
  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
  Int_t mulDigit=clu->GetNCells() ;

  Float_t tMax= 0.; //Time at the maximum
  Float_t eMax=0. ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t absId=clu->GetCellAbsId(iDigit) ;
    Bool_t isHG=cells->GetCellHighGain(absId) ;
    if( elist[iDigit]>eMax){
      tMax=CalibrateTOF(cells->GetCellTime(absId),absId,isHG) ;
      eMax=elist[iDigit] ;
    }
  }
  
   //Try to improve accuracy 
  //Do not account time of soft cells:
  //  const Double_t part=0.5 ;
  Double_t eMin=TMath::Min(0.5,0.2*eMax) ;
  Float_t wtot = 0.;
  Double_t t = 0. ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t absId=clu->GetCellAbsId(iDigit) ;
    Bool_t isHG=cells->GetCellHighGain(absId) ;
      
    Double_t ti=CalibrateTOF(cells->GetCellTime(absId),absId,isHG) ;
    if(TMath::Abs(ti-tMax)>50.e-9) //remove soft cells with wrong time
      continue ;
    
    //Remove too soft cells
    if(elist[iDigit]<eMin)
      continue ;
    
    if(elist[iDigit]>0){ 
      //weight = 1./sigma^2
      //Sigma is parameterization of TOF resolution 16.05.2013
      Double_t wi2=0.;
      if(isHG)
	wi2=1./(2.4 + 3.9/elist[iDigit]) ;
      else
	wi2=1./(2.4 + 3.9/(0.1*elist[iDigit])) ; //E of LG digit is 1/16 of correcponding HG  
      t+=ti*wi2 ;
      wtot+=wi2 ;
    }
  }
  if(wtot>0){
    t=t/wtot ;
  }
  else{
   t=tMax ; 
  }  
  
  return t ;
     
} 
//________________________________________________________________________
Double_t AliPHOSTenderSupply::CalibrateTOF(Double_t tof, Int_t absId, Bool_t isHG){
  //Apply time re-calibration separately for HG and LG channels
  //By default (if not filled) shifts are zero.  
    
  Int_t relId[4];
  fPHOSGeo->AbsToRelNumbering(absId,relId) ;
  Int_t   module = relId[0];
  Int_t   column = relId[3];
  Int_t   row    = relId[2];
  if(isHG)
    tof-=fPHOSCalibData->GetTimeShiftEmc(module, column, row);
  else{
    tof-=fPHOSCalibData->GetLGTimeShiftEmc(module, column, row);
  }
  return tof ;
  
}
//________________________________________________________________________
void AliPHOSTenderSupply::DistanceToBadChannel(Int_t mod, TVector3 * locPos, Double_t &minDist){
  //Check if distance to bad channel was reduced
  Int_t range = minDist/2.2 +1 ; //Distance at which bad channels should be serached
  
  Int_t relid[4]={0,0,0,0} ;
  fPHOSGeo->RelPosToRelId(mod, locPos->X(), locPos->Z(), relid) ; 
  Int_t xmin=TMath::Max(1,relid[2]-range) ;
  Int_t xmax=TMath::Min(64,relid[2]+range) ;
  Int_t zmin=TMath::Max(1,relid[3]-range) ;
  Int_t zmax=TMath::Min(56,relid[3]+range) ;
  
  Float_t x=0.,z=0.;
  for(Int_t ix=xmin;ix<=xmax;ix++){
    for(Int_t iz=zmin;iz<=zmax;iz++){
      if(fPHOSBadMap[mod] && fPHOSBadMap[mod]->GetBinContent(ix,iz)>0){ //Bad channel
        Int_t relidBC[4]={mod,0,ix,iz} ;
        fPHOSGeo->RelPosInModule(relidBC,x,z); 
        Double_t dist = TMath::Sqrt((x-locPos->X())*(x-locPos->X()) + (z-locPos->Z())*(z-locPos->Z()));
        if(dist<minDist) minDist = dist;
      }
    }  
  }
  
}


