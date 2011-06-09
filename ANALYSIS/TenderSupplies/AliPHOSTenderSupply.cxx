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

AliPHOSTenderSupply::AliPHOSTenderSupply() :
  AliTenderSupply()
  ,fOCDBpass("local://OCDB")
  ,fNonlinearityVersion("Default")
  ,fPHOSGeo(0x0)
  ,fPHOSCalibData(0x0)
{
	//
	// default ctor
	//
}

//_____________________________________________________
AliPHOSTenderSupply::AliPHOSTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender)
  ,fOCDBpass("alien:///alice/cern.ch/user/p/prsnko/PHOSrecalibrations/")
  ,fNonlinearityVersion("Default")
  ,fPHOSGeo(0x0)
  ,fPHOSCalibData(0x0)
{
	//
	// named ctor
	//
}

//_____________________________________________________
AliPHOSTenderSupply::~AliPHOSTenderSupply()
{
  //Destructor
  if(fPHOSCalibData)delete fPHOSCalibData;
}

//_____________________________________________________
void AliPHOSTenderSupply::Init()
{
  //
  // Initialise PHOS tender
  //
    

  
}

//_____________________________________________________
void AliPHOSTenderSupply::ProcessEvent()
{
  //Choose PHOS clusters and recalibrate them
  //that it recalculate energy, position and distance 
  //to closest track extrapolation	

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;

  // Init goemetry
  if(!fPHOSGeo){
    fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
    for(Int_t mod=0; mod<5; mod++) {
      if(!event->GetPHOSMatrix(mod)) continue;
      fPHOSGeo->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
    }
  }


  if(!fPHOSCalibData || fTender->RunChanged()){
    AliCDBManager * man = AliCDBManager::Instance();
    man->SetRun(event->GetRunNumber()) ;
    //    man->SetDefaultStorage("local://OCDB");
    man->SetSpecificStorage("PHOS/Calib/EmcGainPedestals",fOCDBpass);
    if(fPHOSCalibData) delete fPHOSCalibData; 
    fPHOSCalibData = new AliPHOSCalibData();
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
    Double_t r=0.,dx=0.,dz=0. ;
    TVector3 locPos;
    TVector3 globPos(xyz) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(globPos,relId) ;
    Int_t mod  = relId[0] ;
    fPHOSGeo->Global2Local(locPos,globPos,mod) ;

    FindTrackMatching(mod,&locPos,r,dx,dz) ;
    clu->SetEmcCpvDistance(r);    
    clu->SetTrackDistance(dx,dz); 
    //    clu->SetChi2(-1);                     //not yet implemented
    clu->SetTOF(cluPHOS1.GetTOF());       

  }


}
//___________________________________________________________________________________________________
void AliPHOSTenderSupply::FindTrackMatching(Int_t mod,TVector3 *locpos,Double_t &r,Double_t &dx, Double_t &dz){
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
      }
    }
  } //Scanned all tracks
  r=TMath::Sqrt(minDistance) ;
  
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




