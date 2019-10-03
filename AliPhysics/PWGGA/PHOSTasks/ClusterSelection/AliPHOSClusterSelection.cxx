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

/* $Id$ */

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultiInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDCaloCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAODCluster.h"
#include "AliVCaloCells.h"
#include "AliOADBContainer.h"
#include "AliPHOSGeometry.h"
#include "AliAnalysisManager.h"



#include "AliPHOSClusterSelection.h"

AliPHOSClusterSelection::AliPHOSClusterSelection()
  : fMinChargedParticleTrackDistance(-1.), 
    fNotUnfolded(false),
    fMaxDispR2(-1.),
    fIsCore(false),
    fMaxTOF(-1.)
{
  // Defaults to the most lenient selection allowable
	return;
}

AliPHOSClusterSelection::~AliPHOSClusterSelection()
{
}

Bool_t AliPHOSClusterSelection::IsSelected(AliVCluster* cluster) const
{
  return IsSelectedCPV(cluster)
    && IsSelectedUnfolded(cluster)
    && IsSelectedDisp(cluster)
    && IsSelectedTOF(cluster);
}

Bool_t AliPHOSClusterSelection::IsSelectedCPV(AliVCluster* cluster) const
{
  if( 0 > fMinChargedParticleTrackDistance )//No selection on CPV
    return true; 

  else{
    Bool_t cpvBit = true; //Changed to false if there is a matching track in the requested radius.

    Double_t dx=cluster->GetTrackDx();//Position for Track matching
    Double_t dz=cluster->GetTrackDz();


    AliESDEvent* EventESD = dynamic_cast<AliESDEvent*> (AliPHOSClusterSelection::GetCurrentEvent());//dynamic cast to test for ESD or AOD
    AliAODEvent* EventAOD = dynamic_cast<AliAODEvent*> (AliPHOSClusterSelection::GetCurrentEvent());
    Double_t mf = 0.; //
    if(EventAOD) mf = EventAOD->GetMagneticField(); //Positive for ++ and negative for -- 
    else if(EventESD) mf = EventESD->GetMagneticField(); //Positive for ++ and negative for --


    if(EventESD){//if ESD
      AliESDCaloCluster * ESDcluster = static_cast<AliESDCaloCluster*> (cluster);//Know its ESD so static cast is fine?
      TArrayI * itracks = ESDcluster-> GetTracksMatched() ;
      if(itracks->GetSize()>0){
    	Int_t iTr = itracks->At(0);
    	if(iTr>=0 && iTr<EventESD->GetNumberOfTracks()){
    	  AliVParticle* track = EventESD->GetTrack(iTr);
    	  Double_t pt = track->Pt() ;
    	  Short_t charge = track->Charge() ;
    	  Double_t r = AliPHOSClusterSelection::TestCPV(dx, dz, pt, charge, mf);
    	  cpvBit=(r>fMinChargedParticleTrackDistance) ;
    	}
      }
    }


    if(EventAOD){//if AOD
      AliAODCluster * AODcluster = static_cast<AliAODCluster*> (cluster);
      int nTracksMatched = AODcluster->GetNTracksMatched();
      if(nTracksMatched > 0) {
	AliVTrack* track = dynamic_cast<AliVTrack*> (cluster->GetTrackMatched(0));
	if ( track ) {
	  Double_t pt = track->Pt();
	  Short_t charge = track->Charge();
	  Double_t r = AliPHOSClusterSelection::TestCPV(dx, dz, pt, charge, mf);
	  cpvBit=(r>fMinChargedParticleTrackDistance) ;
	}
      }
    }
    return cpvBit;
  }
}


Bool_t AliPHOSClusterSelection::IsSelectedUnfolded(AliVCluster* cluster) const
{
  if(!fNotUnfolded)
    return true;
  else{
    Bool_t NotUnfolded = cluster->GetNExMax()<2;//True if it was not unfolded
    return !NotUnfolded;
  }
}

Bool_t AliPHOSClusterSelection::IsSelectedDisp(AliVCluster* cluster) const
{
  if(0 > fMaxDispR2)
    return true;
  else{
    Double_t m02 = 0.,m20 = 0.;
    if(!fIsCore){//No core calculation
      m02 = cluster->GetM02();
      m20 = cluster->GetM20();
      return AliPHOSClusterSelection::TestLambda(cluster->E(),m20,m02) ;
    }
    else{//DispCore
      AliVCaloCells* cells = static_cast<AliVCaloCells*> (AliPHOSClusterSelection::GetCurrentEvent()->GetPHOSCells());//Need the cells
      AliPHOSClusterSelection::EvalCoreLambdas(cluster, cells, m02, m20);
      return AliPHOSClusterSelection::TestLambda(cluster->E(),m20,m02);
    }
  }

}

Bool_t AliPHOSClusterSelection::IsSelectedTOF(AliVCluster* cluster) const
{
  if(0 > fMaxTOF)
    return true;
  else{ 
    // Time of Flight (TOF)
    Double_t tof = cluster->GetTOF();//Time of Flight for the cluster
    return  TMath::Abs(tof) < fMaxTOF;//True if the cut is passed
  }
}

AliPHOSClusterSelection* AliPHOSClusterSelection::SetMinChargedParticleTrackDistance(Float_t distance)
{
  // 'distance' set the minimal allowable distance between the cluster
  // and the nearest extrapolated track.
  // If 'distance' is negative, then all clusters are sellected, the selection
  // being "not applied" or "disabled".
  
  fMinChargedParticleTrackDistance = distance;
  return this;
}

AliPHOSClusterSelection* AliPHOSClusterSelection::SetNotUnfolded(Bool_t notUnfolded)
{
  //if notUnfolded true, it rejects Unfolded Clusters
  fNotUnfolded = notUnfolded;
  return this;
}

AliPHOSClusterSelection* AliPHOSClusterSelection::SetMaxDispR2(Float_t maxR2)
{
  // 'maxR2' sets the maximum allowed dispersion.
  // If 'maxR2' is negative, then all clusters are selected, the selection
  // being "not applied" or "disabled".
  
  fMaxDispR2 = maxR2;
  return this;
}

AliPHOSClusterSelection* AliPHOSClusterSelection::SetIsCore(Bool_t isCore)
{
  // 'isCore' sets wether core version of Disp is used. -1 gives no core.
  
  fIsCore = isCore;
  return this;
}

AliPHOSClusterSelection* AliPHOSClusterSelection::SetMaxTOF(Float_t maxTOF)
{
  // 'maxTOF' sets the maximum allowed time of flight for the cluster.
  // If 'maxTOF' is negative, all clusters are selected and the selection is "disabled".
  
  fMaxTOF = maxTOF;
  return this;
}

TString AliPHOSClusterSelection::ToString() const
{
  // returns a string an quasi-unique string for whatever selection
  // parameters the instance contains. The uniqueness of the string
  // is limited by the precision given in the formatting of the string.
  // Take care that the precision is sufficient for your needs.

  TString string ="v1";
  string+=Form("_%.2f", fMinChargedParticleTrackDistance);
  string+=Form("_%i", fNotUnfolded);
  string+=Form("_%.2f", fMaxDispR2);
  string+=Form("_%i", fIsCore);
  string+=Form("_%.2f", fMaxTOF);
  return string;
}


Float_t AliPHOSClusterSelection::GetMinChargedParticleTrackDistance(const TString& string)
{
  TObjArray* array = string.Tokenize("_");
  TObjString* objString = (TObjString*) array->At(kMinChargedParticleTrackDistance);
  Float_t flt = objString->String().Atof();
  delete array;
  return flt;
}

Bool_t AliPHOSClusterSelection::GetUnfolded(const TString& string)
{
  TObjArray* array = string.Tokenize("_");
  TObjString* objString = (TObjString*) array->At(kNotUnfolded);
  Float_t flt = objString->String().Atoi();
  delete array;
  return flt;
}

Float_t AliPHOSClusterSelection::GetMaxDispR2(const TString& string)
{
  TObjArray* array = string.Tokenize("_");
  TObjString* objString = (TObjString*) array->At(kMaxDispR2);
  Float_t flt = objString->String().Atof();
  delete array;
  return flt;
}

Bool_t AliPHOSClusterSelection::GetIsCore(const TString& string)
{
  TObjArray* array = string.Tokenize("_");
  TObjString* objString = (TObjString*) array->At(kIsCore);
  Float_t flt = objString->String().Atoi();
  delete array;
  return flt;
}

Float_t AliPHOSClusterSelection::GetMaxTOF(const TString& string)
{
  TObjArray* array = string.Tokenize("_");
  TObjString* objString = (TObjString*) array->At(kMaxTOF);
  Float_t flt = objString->String().Atof();
  delete array;
  return flt;
}


Double_t AliPHOSClusterSelection::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge, Double_t mf) const {
  //Parameterization of LHC10h period
  //_true if neutral_
  //Copied from Pi0Flow task

  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
			 6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
 
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

//____________________________________________________________________________
void  AliPHOSClusterSelection::EvalCoreLambdas(AliVCluster  * clu, AliVCaloCells * cells,Double_t &m02, Double_t &m20) const
{ 
  //calculate dispecrsion of the cluster in the circle with radius distanceCut around the maximum
  //Copied from pi0flowtask
  AliVEvent* vevent = AliPHOSClusterSelection::GetCurrentEvent();
  Int_t runNumber = vevent->GetRunNumber();
  const Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
  // Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  const Int_t mulDigit=clu->GetNCells() ;
  Double_t xc[mulDigit] ;
  Double_t zc[mulDigit] ;
  Double_t wi[mulDigit] ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  const Double_t logWeight=4.5 ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi=0. ;
    Float_t zi=0. ;
    Int_t absId = clu->GetCellAbsId(iDigit) ;

    AliOADBContainer geomContainer("phosGeo");//Initialize Geometry
    geomContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(runNumber,"PHOSRotationMatrixes");
    AliPHOSGeometry * fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
    for(Int_t mod=0; mod<5; mod++) {
      if(!matrixes->At(mod)) {
	AliInfo(Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
	continue;
      }


      fPHOSGeo->AbsToRelNumbering(absId, relid) ;
      fPHOSGeo->RelPosInModule(relid, xi, zi);
      xc[iDigit]=xi ;
      zc[iDigit]=zi ;
      Double_t ei = elist[iDigit]*cells->GetCellAmplitude(absId) ;
      wi[iDigit]=0. ;
      if (clu->E()>0 && ei>0) {
	wi[iDigit] = TMath::Max( 0., logWeight + TMath::Log( ei / clu->E() ) ) ;
	Double_t w=wi[iDigit];
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
    for(Int_t iDigit1=0; iDigit1<mulDigit; iDigit1++) {
      Double_t w=wi[iDigit1];
      if (w>0.) {
        Double_t xi1= xc[iDigit1] ;
        Double_t zi1= zc[iDigit1] ;
	if((xi1-x)*(xi1-x)+(zi1-z)*(zi1-z) < 4.5*4.5){
          xCut += w * xi1 ;
          zCut += w * zi1 ; 
          dxx  += w * xi1 * xi1 ;
          dzz  += w * zi1 * zi1 ;
          dxz  += w * xi1 * zi1 ; 
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
}

//_____________________________________________________________________________
Bool_t AliPHOSClusterSelection::TestLambda(Double_t pt,Double_t l1,Double_t l2) const 
{
  //Evaluates if lambdas correspond to photon cluster
  //Tuned using pp data 
  //copied from Pi0FlowTask
  Double_t l2Mean, l1Mean, l2Sigma, l1Sigma, c, R2;
  if(! fIsCore ){
    l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
    l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
    l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
    l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
    c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
    R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma +
      0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
      0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  }
  else{
    //For core radius R=4.5
    l1Mean  = 1.150200 + 0.097886/(1.+1.486645*pt+0.000038*pt*pt) ;
    l2Mean = 1.574706 + 0.997966*exp(-0.895075*pt)-0.010666*pt ;
    l1Sigma = 0.100255 + 0.337177*exp(-0.517684*pt)+0.001170*pt ;
    l2Sigma = 0.232580 + 0.573401*exp(-0.735903*pt)-0.002325*pt ;
    c = -0.110983 -0.017353/(1.-1.836995*pt+0.934517*pt*pt) ;
    R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
      0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
      0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  }
  return (R2<fMaxDispR2*fMaxDispR2);//TODO: Check if this should be 2. order
}


AliVEvent* AliPHOSClusterSelection::GetCurrentEvent() const
{
  // Hackish way of getting the current event.
  // Its probably not appropriate to call this function outside of
  // AliAnalysisTaskSE::UserExec
  
  AliAnalysisManager* analysisManager = dynamic_cast<AliAnalysisManager*>(AliAnalysisManager::GetAnalysisManager());
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(analysisManager->GetInputEventHandler());
  AliMultiInputEventHandler *multiInputHandler = dynamic_cast<AliMultiInputEventHandler *>(inputHandler);
  if (multiInputHandler)
    inputHandler = dynamic_cast<AliInputEventHandler *>(multiInputHandler->GetFirstInputEventHandler());
  
  AliVEvent* inputEvent = dynamic_cast<AliVEvent*>(inputHandler->GetEvent());
  if( ! inputEvent ) 
    AliError("Was not able to retrieve event!");
  
  return inputEvent;
}
