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

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliTender.h>
#include <AliCDBManager.h>
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"

#include "AliESDCaloCluster.h"
#include "AliPHOSTenderSupply.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliOADBContainer.h"

ClassImp(AliPHOSTenderSupply)

AliPHOSTenderSupply::AliPHOSTenderSupply() :
  AliTenderSupply()
  ,fOCDBpass("local://OCDB")
  ,fNonlinearityVersion("Default")
  ,fPHOSGeo(0x0)
  ,fRunNumber(0)
  ,fRecoPass(-1)  //to be defined
  ,fUsePrivateBadMap(0)
  ,fUsePrivateCalib(0)
  ,fPHOSCalibData(0x0)
{
	//
	// default ctor
	//
   for(Int_t i=0;i<10;i++)fNonlinearityParams[i]=0. ;
   for(Int_t mod=0;mod<5;mod++)fPHOSBadMap[mod]=0x0 ;
}

//_____________________________________________________
AliPHOSTenderSupply::AliPHOSTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender)
  ,fOCDBpass("alien:///alice/cern.ch/user/p/prsnko/PHOSrecalibrations/")
  ,fNonlinearityVersion("Default")
  ,fPHOSGeo(0x0)
  ,fRunNumber(0)
  ,fRecoPass(-1)  //to be defined
  ,fUsePrivateBadMap(0)
  ,fUsePrivateCalib(0)
  ,fPHOSCalibData(0x0)
{
	//
	// named ctor
	//
   for(Int_t i=0;i<10;i++)fNonlinearityParams[i]=0. ;
   for(Int_t mod=0;mod<5;mod++)fPHOSBadMap[mod]=0x0 ;
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
  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
	      
  fRunNumber=event->GetRunNumber();	      
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
      }
    }
    if(fRecoPass<0){
      AliError("Can not find pass number from file name, set it manually");
    }
  }
   
  //Init geometry 
  if(!fPHOSGeo){
    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
    fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
    for(Int_t mod=0; mod<5; mod++) {
      if(!matrixes->At(mod)) continue;
      fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
      printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
    }
  }
  
  //Init Bad channels map
  if(!fUsePrivateBadMap){
   AliOADBContainer badmapContainer(Form("phosBadMap"));
    badmapContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSBadMaps.root","phosBadMap");
    TObjArray *maps = (TObjArray*)badmapContainer.GetObject(fRunNumber,"phosBadMap");
    if(!maps){
      AliError(Form("Can not read Bad map for run %d. \n You may choose to use your map with ForceUsingBadMap()\n",fRunNumber)) ;    
    }
    else{
      AliInfo(Form("Setting PHOS bad map with name %s \n",maps->GetName())) ;
      for(Int_t mod=0; mod<5;mod++){
        if(fPHOSBadMap[mod]) 
          delete fPHOSBadMap[mod] ;
        TH2I * h = (TH2I*)maps->At(mod) ;      
	if(h)
          fPHOSBadMap[mod]=new TH2I(*h) ;
      }
    }    
  } 

  if(!fUsePrivateCalib){
    //Init recalibration
    //Check the pass1-pass2-pass3 reconstruction
    AliOADBContainer calibContainer("phosRecalibration");
    calibContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");
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

//_____________________________________________________
void AliPHOSTenderSupply::ProcessEvent()
{
  //Choose PHOS clusters and recalibrate them
  //that it recalculate energy, position and distance 
  //to closest track extrapolation	

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
	      
  fRunNumber=event->GetRunNumber();	      

  if(!fPHOSCalibData || fTender->RunChanged()){
    InitTender();
    
  }


  const AliESDVertex *esdVertex = event->GetPrimaryVertex();
  AliESDCaloCells * cells = event->GetPHOSCells() ;
  TVector3 vertex(esdVertex->GetX(),esdVertex->GetY(),esdVertex->GetZ());
  if(vertex.Mag()>99.) //vertex not defined?
    vertex.SetXYZ(0.,0.,0.) ;

  //For re-calibration
  const Double_t logWeight=4.5 ;  

  Int_t multClust = event->GetNumberOfCaloClusters();
  for (Int_t i=0; i<multClust; i++) {
    AliESDCaloCluster *clu = event->GetCaloCluster(i);
    if ( !clu->IsPHOS()) continue;
    
    Float_t  position[3];
    clu->GetPosition(position);
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
   //Apply re-Calibreation
    AliPHOSEsdCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cells); // modify the cell energies
    cluPHOS1.EvalAll(logWeight,vertex);         // recalculate the cluster parameters
    cluPHOS1.SetE(CorrectNonlinearity(cluPHOS1.E()));// Users's nonlinearity
    
    Float_t  xyz[3];
    cluPHOS1.GetPosition(xyz);
    clu->SetPosition(xyz);                       //rec.point position in MARS
    clu->SetE(cluPHOS1.E());                           //total or core particle energy
    clu->SetDispersion(cluPHOS1.GetDispersion());  //cluster dispersion
    //    ec->SetPID(rp->GetPID()) ;            //array of particle identification
    clu->SetM02(cluPHOS1.GetM02()) ;               //second moment M2x
    clu->SetM20(cluPHOS1.GetM20()) ;               //second moment M2z
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
    clu->SetChi2(TestLambda(clu->E(),clu->GetM20(),clu->GetM02()));                     //not yet implemented
    clu->SetTOF(cluPHOS1.GetTOF());       

  }


}
//___________________________________________________________________________________________________
void AliPHOSTenderSupply::FindTrackMatching(Int_t mod,TVector3 *locpos,
					    Double_t &dx, Double_t &dz,
					    Double_t &pt,Int_t &charge){
  //Find track with closest extrapolation to cluster
  
  AliESDEvent *event= fTender->GetEvent();
  Double_t  magF = event->GetMagneticField();
  Double_t magSign = 1.0;
  if(magF<0)magSign = -1.0;
  
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  Int_t nt=event->GetNumberOfTracks();
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
    AliESDtrack *esdTrack=event->GetTrack(i);

    // Skip the tracks having "wrong" status (has to be checked/tuned)
    ULong_t status = esdTrack->GetStatus();
    if ((status & AliESDtrack::kTPCout)   == 0) continue;
    //     if ((status & AliESDtrack::kTRDout)   == 0) continue;
    //     if ((status & AliESDtrack::kTRDrefit) == 1) continue;
    
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
  } //Scanned all tracks
  
}
//____________________________________________________________
Float_t AliPHOSTenderSupply::CorrectNonlinearity(Float_t en){

  //For backward compatibility, if no RecoParameters found
  if(fNonlinearityVersion=="Default"){
    return 0.0241+1.0504*en+0.000249*en*en ;
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
Double_t AliPHOSTenderSupply::TestLambda(Double_t pt,Double_t l1,Double_t l2){
  
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<2.5*2.5) ;
  
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
  AliESDEvent *event=fTender->GetEvent();
  Double_t mf = event->GetMagneticField(); //Positive for ++ and negative for --

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
  
  if(mod>5 || mod<1){
    AliError(Form("No bad map for PHOS module %d ",mod)) ;
    return kTRUE ;
  }
  if(!fPHOSBadMap[mod]){
     AliError(Form("No Bad map for PHOS module %d",mod)) ;
     return kTRUE ;
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
    AliError(Form("Can not open BadMaps file %s \n",filename)) ;
    return ;
  }
  gROOT->cd() ;
  char key[55] ;
  for(Int_t mod=1;mod<4; mod++){
    sprintf(key,"PHOS_BadMap_mod%d",mod) ;
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
    AliFatal(Form("Can not open Calibration file %s \n",filename)) ;
    return ;
  }
  fPHOSCalibData = (AliPHOSCalibData*)fc->Get("PHOSCalibration") ;
  fc->Close() ;
  fUsePrivateCalib=kTRUE; 
}