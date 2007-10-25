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
//_________________________________________________________________________
//  Reconstructed Points for the EMCAL
//  A RecPoint is a cluster of digits
//*-- Author: Yves Schutz (SUBATECH)
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)
//*-- Author: Heather Gray (LBL) merged AliEMCALRecPoint and AliEMCALTowerRecPoint 02/04

// --- ROOT system ---
#include "TPad.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"
#include "TGeoPhysicalNode.h"

// --- Standard library ---
#include <Riostream.h>

// --- AliRoot header files ---
//#include "AliGenerator.h"
class AliGenerator;
class AliEMCAL;
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALHit.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"

ClassImp(AliEMCALRecPoint)

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint()
  : AliCluster(), fGeomPtr(0),
    fAmp(0), fIndexInList(-1), //to be set when the point is already stored
    fLocPos(0,0,0), fLocPosM(0),
    fMaxDigit(100), fMulDigit(0), fMaxTrack(200),
    fMulTrack(0), fDigitsList(0), fTracksList(0),
    fClusterType(-1), fCoreEnergy(0), fDispersion(0),
    fEnergyList(0), fTimeList(0), fAbsIdList(0),
    fTime(0.), fCoreRadius(10),  //HG check this
    fDETracksList(0), fMulParent(0), fMaxParent(0),
    fParentsList(0), fDEParentsList(0), fSuperModuleNumber(0),
    fDigitIndMax(-1)
{
  // ctor
  fGeomPtr = AliEMCALGeometry::GetInstance();
  
  fLambda[0] = 0;
  fLambda[1] = 0;

}

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint(const char *) 
  : AliCluster(), fGeomPtr(0),
    fAmp(0), fIndexInList(-1), //to be set when the point is already stored
    fLocPos(0,0,0), fLocPosM(new TMatrixF(3,3)),
    fMaxDigit(100), fMulDigit(0), fMaxTrack(1000), fMulTrack(0),
    fDigitsList(new Int_t[fMaxDigit]), fTracksList(new Int_t[fMaxTrack]),
    fClusterType(-1), fCoreEnergy(0), fDispersion(0),
    fEnergyList(new Float_t[fMaxDigit]), fTimeList(new Float_t[fMaxDigit]), 
    fAbsIdList(new Int_t[fMaxDigit]), fTime(-1.), fCoreRadius(10),
    fDETracksList(new Float_t[fMaxTrack]), fMulParent(0), fMaxParent(1000),
    fParentsList(new Int_t[fMaxParent]), fDEParentsList(new Float_t[fMaxParent]),
    fSuperModuleNumber(0), fDigitIndMax(-1)
{
  // ctor
  for (Int_t i = 0; i < fMaxTrack; i++)
    fDETracksList[i] = 0;
  for (Int_t i = 0; i < fMaxParent; i++) {
    fParentsList[i] = -1;
    fDEParentsList[i] = 0;
  }

  fGeomPtr = AliEMCALGeometry::GetInstance();
  fLambda[0] = 0;
  fLambda[1] = 0;
}

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint(const AliEMCALRecPoint & rp) 
  : AliCluster(rp), fGeomPtr(rp.fGeomPtr),
    fAmp(rp.fAmp), fIndexInList(rp.fIndexInList),
    fLocPos(rp.fLocPos), fLocPosM(rp.fLocPosM),
    fMaxDigit(rp.fMaxDigit), fMulDigit(rp.fMulDigit),
    fMaxTrack(rp.fMaxTrack), fMulTrack(rp.fMaxTrack),
    fDigitsList(new Int_t[rp.fMaxDigit]), fTracksList(new Int_t[rp.fMaxTrack]),
    fClusterType(rp.fClusterType), fCoreEnergy(rp.fCoreEnergy), 
    fDispersion(rp.fDispersion),
    fEnergyList(new Float_t[rp.fMaxDigit]), fTimeList(new Float_t[rp.fMaxDigit]), 
    fAbsIdList(new Int_t[rp.fMaxDigit]), fTime(rp.fTime), fCoreRadius(rp.fCoreRadius),
    fDETracksList(new Float_t[rp.fMaxTrack]), fMulParent(rp.fMulParent), 
    fMaxParent(rp.fMaxParent), fParentsList(new Int_t[rp.fMaxParent]), 
    fDEParentsList(new Float_t[rp.fMaxParent]),
    fSuperModuleNumber(rp.fSuperModuleNumber), fDigitIndMax(rp.fDigitIndMax)
{
  //copy ctor
  fLambda[0] = rp.fLambda[0];
  fLambda[1] = rp.fLambda[1];

  for(Int_t i = 0; i < rp.fMulDigit; i++) {
    fEnergyList[i] = rp.fEnergyList[i];
    fTimeList[i] = rp.fTimeList[i];
    fAbsIdList[i] = rp.fAbsIdList[i];
  }

  for(Int_t i = 0; i < rp.fMulTrack; i++) fDETracksList[i] = rp.fDETracksList[i];

  for(Int_t i = 0; i < rp.fMulParent; i++) {
    fParentsList[i] = rp.fParentsList[i];
    fDEParentsList[i] = rp.fDEParentsList[i];
  }

}
//____________________________________________________________________________
AliEMCALRecPoint::~AliEMCALRecPoint()
{
  // dtor
  if ( fEnergyList )
    delete[] fEnergyList ; 
  if ( fTimeList )
    delete[] fTimeList ; 
  if ( fAbsIdList )
    delete[] fAbsIdList ; 
   if ( fDETracksList)
    delete[] fDETracksList;
   if ( fParentsList)
    delete[] fParentsList;
   if ( fDEParentsList)
    delete[] fDEParentsList;

   delete fLocPosM ;
   delete [] fDigitsList ;
   delete [] fTracksList ;
}

//____________________________________________________________________________
AliEMCALRecPoint& AliEMCALRecPoint::operator= (const AliEMCALRecPoint &rp)
{
  if(&rp == this) return *this;

  fGeomPtr = rp.fGeomPtr;
  fAmp = rp.fAmp;
  fIndexInList = rp.fIndexInList;
  fLocPos = rp.fLocPos;
  fLocPosM = rp.fLocPosM;
  fMaxDigit = rp.fMaxDigit;
  fMulDigit = rp.fMulDigit;
  fMaxTrack = rp.fMaxTrack;
  fMulTrack = rp.fMaxTrack;
  for(Int_t i = 0; i<fMaxDigit; i++) fDigitsList[i] = rp.fDigitsList[i];
  for(Int_t i = 0; i<fMaxTrack; i++) fTracksList[i] = rp.fTracksList[i];
  fClusterType = rp.fClusterType;
  fCoreEnergy = rp.fCoreEnergy; 
  fDispersion = rp.fDispersion;
  for(Int_t i = 0; i<fMaxDigit; i++) {
    fEnergyList[i] = rp.fEnergyList[i];
    fTimeList[i] = rp.fTimeList[i]; 
    fAbsIdList[i] = rp.fAbsIdList[i];
  }
  fTime = rp.fTime;
  fCoreRadius = rp.fCoreRadius;
  for(Int_t i = 0; i < fMaxTrack; i++) fDETracksList[i] = rp.fDETracksList[i];
  fMulParent = rp.fMulParent;
  fMaxParent = rp.fMaxParent;
  for(Int_t i = 0; i < fMaxParent; i++) {
    fParentsList[i] = rp.fParentsList[i]; 
    fDEParentsList[i] = rp.fDEParentsList[i];
  }
  fSuperModuleNumber = rp.fSuperModuleNumber;
  fDigitIndMax = rp.fDigitIndMax;

  fLambda[0] = rp.fLambda[0];
  fLambda[1] = rp.fLambda[1];

  return *this;

}

//____________________________________________________________________________
void AliEMCALRecPoint::AddDigit(AliEMCALDigit & digit, Float_t Energy)
{
  // Adds a digit to the RecPoint
  // and accumulates the total amplitude and the multiplicity 
  
  if(fEnergyList == 0)
    fEnergyList =  new Float_t[fMaxDigit]; 
  if(fTimeList == 0)
    fTimeList =  new Float_t[fMaxDigit]; 
  if(fAbsIdList == 0) {
    fAbsIdList =  new Int_t[fMaxDigit];
  }

  if ( fMulDigit >= fMaxDigit ) { // increase the size of the lists 
    fMaxDigit*=2 ; 
    Int_t * tempo = new Int_t[fMaxDigit]; 
    Float_t * tempoE =  new Float_t[fMaxDigit];
    Float_t * tempoT =  new Float_t[fMaxDigit];
    Int_t * tempoId = new Int_t[fMaxDigit]; 

    Int_t index ;     
    for ( index = 0 ; index < fMulDigit ; index++ ){
      tempo[index]   = fDigitsList[index] ;
      tempoE[index]  = fEnergyList[index] ; 
      tempoT[index]  = fTimeList[index] ; 
      tempoId[index] = fAbsIdList[index] ; 
    }
    
    delete [] fDigitsList ;
    delete [] fEnergyList ;
    delete [] fTimeList ;
    delete [] fAbsIdList ;

    fDigitsList = tempo;
    fEnergyList = tempoE; 
    fTimeList = tempoT;
    fAbsIdList = tempoId;
  } // if
  
  fDigitsList[fMulDigit]   = digit.GetIndexInList()  ; 
  fEnergyList[fMulDigit]   = Energy ;
  fTimeList[fMulDigit]     = digit.GetTimeR() ;
  fAbsIdList[fMulDigit]    = digit.GetId();
  fMulDigit++ ; 
  fAmp += Energy ; 

  //JLK 10-Oct-2007 this hasn't been filled before because it was in
  //the wrong place in previous versions.
  //Now we evaluate it only if the supermodulenumber for this recpoint
  //has not yet been set (or is the 0th one)
  if(fSuperModuleNumber == 0)
    fSuperModuleNumber = fGeomPtr->GetSuperModuleNumber(digit.GetId());

}
//____________________________________________________________________________
Bool_t AliEMCALRecPoint::AreNeighbours(AliEMCALDigit * digit1, AliEMCALDigit * digit2 ) const
{
  // Tells if (true) or not (false) two digits are neighbours
  // A neighbour is defined as being two digits which share a corner
  
  static Bool_t areNeighbours = kFALSE ;
  static Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  static int nSupMod1=0, nModule1=0, nIphi1=0, nIeta1=0;
  static Int_t relid1[2] , relid2[2] ; // ieta, iphi
  static Int_t rowdiff=0, coldiff=0;

  areNeighbours = kFALSE ;

  fGeomPtr->GetCellIndex(digit1->GetId(), nSupMod,nModule,nIphi,nIeta);
  fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, relid1[0],relid1[1]);

  fGeomPtr->GetCellIndex(digit2->GetId(), nSupMod1,nModule1,nIphi1,nIeta1);
  fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod1,nModule1,nIphi1,nIeta1, relid2[0],relid2[1]);
  
  rowdiff = TMath::Abs( relid1[0] - relid2[0] ) ;  
  coldiff = TMath::Abs( relid1[1] - relid2[1] ) ;  

  if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0)) 
  areNeighbours = kTRUE ;
  
  return areNeighbours;
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::Compare(const TObject * obj) const
{
  // Compares two RecPoints according to their position in the EMCAL modules

  Float_t delta = 1 ; //Width of "Sorting row". If you change this 
                      //value (what is senseless) change as well delta in
                      //AliEMCALTrackSegmentMakerv* and other RecPoints...
  Int_t rv ; 

  AliEMCALRecPoint * clu = (AliEMCALRecPoint *)obj ; 

  TVector3 locpos1; 
  GetLocalPosition(locpos1);
  TVector3 locpos2;  
  clu->GetLocalPosition(locpos2);  

  Int_t rowdif = (Int_t)(TMath::Ceil(locpos1.X()/delta)-TMath::Ceil(locpos2.X()/delta)) ;
  if (rowdif> 0) 
    rv = 1 ;
  else if(rowdif < 0) 
    rv = -1 ;
  else if(locpos1.Y()>locpos2.Y()) 
    rv = -1 ;
  else 
    rv = 1 ; 

  return rv ; 
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::DistancetoPrimitive(Int_t px, Int_t py)
{
  // Compute distance from point px,py to  a AliEMCALRecPoint considered as a Tmarker
  // Compute the closest distance of approach from point px,py to this marker.
  // The distance is computed in pixels units.
  // HG Still need to update -> Not sure what this should achieve

  TVector3 pos(0.,0.,0.) ;
  GetLocalPosition(pos) ;
  Float_t x =  pos.X() ;
  Float_t y =  pos.Y() ;
  const Int_t kMaxDiff = 10;
  Int_t pxm  = gPad->XtoAbsPixel(x);
  Int_t pym  = gPad->YtoAbsPixel(y);
  Int_t dist = (px-pxm)*(px-pxm) + (py-pym)*(py-pym);
  
  if (dist > kMaxDiff) return 9999;
  return dist;
}

//___________________________________________________________________________
 void AliEMCALRecPoint::Draw(Option_t *option)
 {
   // Draw this AliEMCALRecPoint with its current attributes
   
   AppendPad(option);
 }

//______________________________________________________________________________
void AliEMCALRecPoint::ExecuteEvent(Int_t /*event*/, Int_t, Int_t)
{
  // Execute action corresponding to one event
  // This member function is called when a AliEMCALRecPoint is clicked with the locator
  //
  // If Left button is clicked on AliEMCALRecPoint, the digits are switched on    
  // and switched off when the mouse button is released.

  //  static Int_t pxold, pyold;

  /*  static TGraph *  digitgraph = 0 ;
  static TPaveText* clustertext = 0 ;
  
  if (!gPad->IsEditable()) return;
  
  switch (event) {
    
    
  case kButton1Down:{
    AliEMCALDigit * digit ;

    Int_t iDigit;
    Int_t relid[2] ;
  
    const Int_t kMulDigit=AliEMCALRecPoint::GetDigitsMultiplicity() ;
    Float_t * xi = new Float_t [kMulDigit] ; 
    Float_t * zi = new Float_t [kMulDigit] ;
    
    for(iDigit = 0; iDigit < kMulDigit; iDigit++) {
      Fatal("AliEMCALRecPoint::ExecuteEvent", " -> Something wrong with the code"); 
      digit = 0 ; //dynamic_cast<AliEMCALDigit *>((fDigitsList)[iDigit]);
      fGeomPtr->AbsToRelNumbering(digit->GetId(), relid) ;
      fGeomPtr->PosInAlice(relid, xi[iDigit], zi[iDigit]) ;
    }
    
    if (!digitgraph) {
      digitgraph = new TGraph(fMulDigit,xi,zi);
      digitgraph-> SetMarkerStyle(5) ; 
      digitgraph-> SetMarkerSize(1.) ;
      digitgraph-> SetMarkerColor(1) ;
      digitgraph-> Draw("P") ;
    }
    if (!clustertext) {
      
      TVector3 pos(0.,0.,0.) ;
      GetLocalPosition(pos) ;
      clustertext = new TPaveText(pos.X()-10,pos.Z()+10,pos.X()+50,pos.Z()+35,"") ;
      Text_t  line1[40] ;
      Text_t  line2[40] ;
      sprintf(line1,"Energy=%1.2f GeV",GetEnergy()) ;
      sprintf(line2,"%d Digits",GetDigitsMultiplicity()) ;
      clustertext ->AddText(line1) ;
      clustertext ->AddText(line2) ;
      clustertext ->Draw("");
    }
    gPad->Update() ; 
    Print("") ;
    delete[] xi ; 
    delete[] zi ; 
   }
  
break;
  
  case kButton1Up:
    if (digitgraph) {
      delete digitgraph  ;
      digitgraph = 0 ;
    }
    if (clustertext) {
      delete clustertext ;
      clustertext = 0 ;
    }
    
    break;
    
    }*/
}
//____________________________________________________________________________
void AliEMCALRecPoint::EvalAll(Float_t logWeight,TClonesArray * digits) 
{
  // Evaluates all shower parameters
  EvalLocalPosition(logWeight, digits) ;
  EvalElipsAxis(logWeight, digits) ;
  EvalDispersion(logWeight, digits) ;
  //EvalCoreEnergy(logWeight, digits);
  EvalTime(digits) ;
  EvalPrimaries(digits) ;
  EvalParents(digits);

  //Called last because it sets the global position of the cluster?
  EvalLocal2TrackingCSTransform();

}

//____________________________________________________________________________
void  AliEMCALRecPoint::EvalDispersion(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the dispersion of the shower at the origin of the RecPoint
  // in cell units - Nov 16,2006

  Double_t d = 0., wtot = 0., w = 0.;
  Int_t iDigit=0, nstat=0;
  AliEMCALDigit * digit ;
  
  // Calculates the dispersion in cell units 
  Double_t etai, phii, etaMean=0.0, phiMean=0.0; 
  int nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  int iphi=0, ieta=0;
  // Calculate mean values
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;

    if (fAmp>0 && fEnergyList[iDigit]>0) {
      fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nModule,nIphi,nIeta);
      fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;

      if(w>0.0) {
        phiMean += phii*w;
        etaMean += etai*w;
        wtot    += w;
      }
    }
  }
  if (wtot>0) {
    phiMean /= wtot ;
    etaMean /= wtot ;
  } else AliError(Form("Wrong weight %f\n", wtot));

  // Calculate dispersion
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;

    if (fAmp>0 && fEnergyList[iDigit]>0) {
      fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nModule,nIphi,nIeta);
      fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;

      if(w>0.0) {
        nstat++;
        d += w*((etai-etaMean)*(etai-etaMean)+(phii-phiMean)*(phii-phiMean));
      }
    }
  }
  
  if ( wtot > 0 && nstat>1) d /= wtot ;
  else                      d = 0. ; 

  fDispersion = TMath::Sqrt(d) ;
}

//____________________________________________________________________________
void AliEMCALRecPoint::EvalLocalPosition(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the center of gravity in the local EMCAL-module coordinates 
  //  Info("Print", " logWeight %f : cluster energy %f ", logWeight, fAmp); // for testing
  
  static Double_t dist;

  AliEMCALDigit * digit;
  Int_t i=0, nstat=0, idMax=-1;
  Double_t clXYZ[3]={0.,0.,0.}, clRmsXYZ[3]={0.,0.,0.}, xyzi[3], wtot=0., w=0.;

  //printf(" dist : %f  e : %f \n", dist, fAmp);
  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;
    if(iDigit==0) {
      idMax = digit->GetId(); // is it correct
      dist  = TmaxInCm(Double_t(fAmp));
    }
    fGeomPtr->RelPosCellInSModule(digit->GetId(), idMax, dist, xyzi[0], xyzi[1], xyzi[2]);
    //printf(" Id %i : dist %f Local x,y,z %f %f %f \n", digit->GetId(), dist, xyzi[0], xyzi[1], xyzi[2]);

    //fGeomPtr->RelPosCellInSModule(digit->GetId(), xyzi[0], xyzi[1], xyzi[2]);
    //printf(" Id %i : dist %f Local x,y,z %f %f %f \n", digit->GetId(), 0.0, xyzi[0], xyzi[1], xyzi[2]);
    //    if(fAmp>102.) assert(0);

    if(logWeight > 0.0)  w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ));
    else  w = fEnergyList[iDigit]; // just energy

    if(w>0.0) {
      wtot += w ;
      nstat++;
      for(i=0; i<3; i++ ) {
        clXYZ[i]    += (w*xyzi[i]);
        clRmsXYZ[i] += (w*xyzi[i]*xyzi[i]);
      }
    }
  }
  //  cout << " wtot " << wtot << endl;
  if ( wtot > 0 ) { 
    //    xRMS   = TMath::Sqrt(x2m - xMean*xMean);
    for(i=0; i<3; i++ ) {
      clXYZ[i] /= wtot;
      if(nstat>1) {
        clRmsXYZ[i] /= (wtot*wtot);  
        clRmsXYZ[i]  =  clRmsXYZ[i] - clXYZ[i]*clXYZ[i];
	if(clRmsXYZ[i] > 0.0) {
          clRmsXYZ[i] =  TMath::Sqrt(clRmsXYZ[i]);
        } else      clRmsXYZ[i] = 0;
      } else        clRmsXYZ[i] = 0;
    }    
  } else {
    for(i=0; i<3; i++ ) {
      clXYZ[i] = clRmsXYZ[i] = -1.;
    }
  }
  // clRmsXYZ[i] ??
  fLocPos.SetX(clXYZ[0]);
  fLocPos.SetY(clXYZ[1]);
  fLocPos.SetZ(clXYZ[2]);
    
//  if (gDebug==2)
//    printf("EvalLocalPosition: eta,phi,r = %f,%f,%f", fLocPos.X(), fLocPos.Y(), fLocPos.Z()) ; 
  fLocPosM = 0 ; // covariance matrix
}

//____________________________________________________________________________
void AliEMCALRecPoint::EvalLocalPositionFit(Double_t deff, Double_t logWeight, 
Double_t phiSlope, TClonesArray * digits)
{
  // Aug 14-16, 2007 - for fit 
  // Aug 31 - should be static ??
  static Double_t dist, ycorr;
  static AliEMCALDigit *digit;

  Int_t i=0, nstat=0, idMax=-1;
  Double_t clXYZ[3]={0.,0.,0.}, clRmsXYZ[3]={0.,0.,0.}, xyzi[3], wtot=0., w=0.; 

  for(Int_t iDigit=0; iDigit<digits->GetEntries(); iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;
    if(iDigit==0) {
      idMax = digit->GetId(); // is it correct
      dist  = TmaxInCm(Double_t(fAmp));
    }

    dist = deff;
    fGeomPtr->RelPosCellInSModule(digit->GetId(), idMax, dist, xyzi[0], xyzi[1], xyzi[2]);

    if(logWeight > 0.0)  w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ));
    else                 w = fEnergyList[iDigit]; // just energy

    if(w>0.0) {
      wtot += w ;
      nstat++;
      for(i=0; i<3; i++ ) {
        clXYZ[i]    += (w*xyzi[i]);
        clRmsXYZ[i] += (w*xyzi[i]*xyzi[i]);
      }
    }
  }
  //  cout << " wtot " << wtot << endl;
  if ( wtot > 0 ) { 
    //    xRMS   = TMath::Sqrt(x2m - xMean*xMean);
    for(i=0; i<3; i++ ) {
      clXYZ[i] /= wtot;
      if(nstat>1) {
        clRmsXYZ[i] /= (wtot*wtot);  
        clRmsXYZ[i]  =  clRmsXYZ[i] - clXYZ[i]*clXYZ[i];
	if(clRmsXYZ[i] > 0.0) {
          clRmsXYZ[i] =  TMath::Sqrt(clRmsXYZ[i]);
        } else      clRmsXYZ[i] = 0;
      } else        clRmsXYZ[i] = 0;
    }    
  } else {
    for(i=0; i<3; i++ ) {
      clXYZ[i] = clRmsXYZ[i] = -1.;
    }
  }
  // clRmsXYZ[i] ??
  if(phiSlope != 0.0 && logWeight > 0.0 && wtot) { 
    // Correction in phi direction (y - coords here); Aug 16;
    // May be put to global level or seperate method
    ycorr = clXYZ[1] * (1. + phiSlope);
    //printf(" y %f : ycorr %f : slope %f \n", clXYZ[1], ycorr, phiSlope); 
    clXYZ[1] = ycorr;
  }
  fLocPos.SetX(clXYZ[0]);
  fLocPos.SetY(clXYZ[1]);
  fLocPos.SetZ(clXYZ[2]);
    
//  if (gDebug==2)
//    printf("EvalLocalPosition: eta,phi,r = %f,%f,%f", fLocPos.X(), fLocPos.Y(), fLocPos.Z()) ; 
  fLocPosM = 0 ; // covariance matrix
}

//_____________________________________________________________________________
Bool_t AliEMCALRecPoint::EvalLocalPosition2(TClonesArray * digits, TArrayD &ed)
{
  // Evaluated local position of rec.point using digits 
  // and parametrisation of w0 and deff
  //printf(" <I> AliEMCALRecPoint::EvalLocalPosition2() \n"); 
  return AliEMCALRecPoint::EvalLocalPositionFromDigits(digits, ed, fLocPos);
}

//_____________________________________________________________________________
Bool_t AliEMCALRecPoint::EvalLocalPositionFromDigits(TClonesArray *digits, TArrayD &ed, TVector3 &locPos)
{
  // Used when digits should be recalibrated
  static Double_t deff, w0, esum;
  static Int_t iDigit;
  //  static AliEMCALDigit *digit;

  if(ed.GetSize() && (digits->GetEntries()!=ed.GetSize())) return kFALSE;

  // Calculate sum energy of digits
  esum = 0.0;
  for(iDigit=0; iDigit<ed.GetSize(); iDigit++) esum += ed[iDigit];

  GetDeffW0(esum, deff, w0);
  
  return EvalLocalPositionFromDigits(esum, deff, w0, digits, ed, locPos); 
}

//_____________________________________________________________________________
Bool_t AliEMCALRecPoint::EvalLocalPositionFromDigits(const Double_t esum, const Double_t deff, const Double_t w0, TClonesArray *digits, TArrayD &ed, TVector3 &locPos)
{
  static AliEMCALDigit *digit;

  Int_t i=0, nstat=0, idMax=-1;
  Double_t clXYZ[3]={0.,0.,0.}, xyzi[3], wtot=0., w=0.; 

  // Get pointer to EMCAL geometry
  // (can't use fGeomPtr in static method)
  AliEMCALGeometry* geo = AliEMCALGeometry::GetInstance(); 

  for(Int_t iDigit=0; iDigit<digits->GetEntries(); iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(iDigit));

    geo->RelPosCellInSModule(digit->GetId(), idMax, deff, xyzi[0], xyzi[1], xyzi[2]);

    if(w0 > 0.0)  w = TMath::Max( 0., w0 + TMath::Log(ed[iDigit] / esum));
    else          w = ed[iDigit]; // just energy

    if(w>0.0) {
      wtot += w ;
      nstat++;
      for(i=0; i<3; i++ ) {
        clXYZ[i] += (w*xyzi[i]);
      }
    }
  }
  //  cout << " wtot " << wtot << endl;
  if (wtot > 0) { 
    for(i=0; i<3; i++ ) {
      clXYZ[i] /= wtot;
    }
    locPos.SetX(clXYZ[0]);
    locPos.SetY(clXYZ[1]);
    locPos.SetZ(clXYZ[2]);
    return kTRUE;
  } else {
    return kFALSE;
  }

}

//_____________________________________________________________________________
void AliEMCALRecPoint::GetDeffW0(const Double_t esum , Double_t &deff,  Double_t &w0)
{
  //
  // Aug 31, 2001 
  // Applied for simulation data with threshold 3 adc
  // Calculate efective distance (deff) and weigh parameter (w0) 
  // for coordinate calculation; 0.5 GeV < esum <100 GeV.
  // Look to:  http://rhic.physics.wayne.edu/~pavlinov/ALICE/SHISHKEBAB/RES/CALIB/GEOMCORR/deffandW0VaEgamma_2.gif
  //
  static Double_t e=0.0;
  const  Double_t dp0=9.25147, dp1=1.16700; // Hard coded now
  const  Double_t wp0=4.83713, wp1=-2.77970e-01, wp2 = 4.41116;

  // No extrapolation here
  e = esum<0.5?0.5:esum;
  e = e>100.?100.:e;

  deff = dp0 + dp1*TMath::Log(e);
  w0   = wp0 / (1. + TMath::Exp(wp1*(e+wp2)));
  //printf("<I> AliEMCALRecPoint::GetDeffW0 esum %5.2f : deff %5.2f : w0 %5.2f \n", esum, deff, w0); 
}

//______________________________________________________________________________
void AliEMCALRecPoint::EvalCoreEnergy(Float_t logWeight, TClonesArray * digits)
{
  // This function calculates energy in the core, 
  // i.e. within a radius rad = fCoreEnergy around the center. Beyond this radius
  // in accordance with shower profile the energy deposition 
  // should be less than 2%
  // Unfinished - Nov 15,2006
  // Distance is calculate in (phi,eta) units

  AliEMCALDigit * digit ;

  Int_t iDigit;

  if (!fLocPos.Mag()) {
    EvalLocalPosition(logWeight, digits);
  }
  
  Double_t phiPoint = fLocPos.Phi(), etaPoint = fLocPos.Eta();
  Double_t eta, phi, distance;
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) ( digits->At(fDigitsList[iDigit]) ) ;
    
    eta = phi = 0.0;
    fGeomPtr->EtaPhiFromIndex(digit->GetId(),eta, phi) ;
    phi = phi * TMath::DegToRad();
  
    distance = TMath::Sqrt((eta-etaPoint)*(eta-etaPoint)+(phi-phiPoint)*(phi-phiPoint));
    if(distance < fCoreRadius)
      fCoreEnergy += fEnergyList[iDigit] ;
  }
  
}
//____________________________________________________________________________
void  AliEMCALRecPoint::EvalElipsAxis(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the axis of the shower ellipsoid in eta and phi
  // in cell units

  static TString gn(fGeomPtr->GetName());

  Double_t wtot = 0.;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;

  AliEMCALDigit * digit = 0;

  Double_t etai , phii, w; 
  int nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  int iphi=0, ieta=0;
  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
    etai = phii = 0.; 
    if(gn.Contains("SHISH")) { 
    // Nov 15,2006 - use cell numbers as coordinates
    // Copied for shish-kebab geometry, ieta,iphi is cast as double as eta,phi
    // We can use the eta,phi(or coordinates) of cell
      nSupMod = nModule = nIphi = nIeta = iphi = ieta = 0;

      fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nModule,nIphi,nIeta);
      fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
    } else { // 
      fGeomPtr->EtaPhiFromIndex(digit->GetId(), etai, phii);
      phii = phii * TMath::DegToRad();
    }

    w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    // fAmp summed amplitude of digits, i.e. energy of recpoint 
    // Gives smaller value of lambda than log weight  
    // w = fEnergyList[iDigit] / fAmp; // Nov 16, 2006 - try just energy

    dxx  += w * etai * etai ;
    x    += w * etai ;
    dzz  += w * phii * phii ;
    z    += w * phii ; 

    dxz  += w * etai * phii ; 

    wtot += w ;
  }

  if ( wtot > 0 ) { 
    dxx /= wtot ;
    x   /= wtot ;
    dxx -= x * x ;
    dzz /= wtot ;
    z   /= wtot ;
    dzz -= z * z ;
    dxz /= wtot ;
    dxz -= x * z ;

    fLambda[0] =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    if(fLambda[0] > 0)
      fLambda[0] = TMath::Sqrt(fLambda[0]) ;
    else
      fLambda[0] = 0;
    
    fLambda[1] =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;

    if(fLambda[1] > 0) //To avoid exception if numerical errors lead to negative lambda.
      fLambda[1] = TMath::Sqrt(fLambda[1]) ;
    else
      fLambda[1]= 0. ;
  } else { 
    fLambda[0]= 0. ;
    fLambda[1]= 0. ;
  }

  //  printf("Evalaxis: lambdas  = %f,%f", fLambda[0],fLambda[1]) ; 

}

//______________________________________________________________________________
void  AliEMCALRecPoint::EvalPrimaries(TClonesArray * digits)
{
  // Constructs the list of primary particles (tracks) which 
  // have contributed to this RecPoint and calculate deposited energy 
  // for each track
  
  AliEMCALDigit * digit ;
  Int_t * primArray = new Int_t[fMaxTrack] ;
  Float_t * dEPrimArray = new Float_t[fMaxTrack] ;

  Int_t index ;  
  for ( index = 0 ; index < GetDigitsMultiplicity() ; index++ ) { // all digits
    digit = dynamic_cast<AliEMCALDigit *>(digits->At( fDigitsList[index] )) ; 
    Int_t nprimaries = digit->GetNprimary() ;
    if ( nprimaries == 0 ) continue ;
    Int_t jndex ;
    for ( jndex = 0 ; jndex < nprimaries ; jndex++ ) { // all primaries in digit
      if ( fMulTrack > fMaxTrack ) {
	fMulTrack = fMaxTrack ;
	Error("EvalPrimaries", "increase fMaxTrack ")  ;
	break ;
      }
      Int_t newPrimary = digit->GetPrimary(jndex+1);
      Float_t dEPrimary = digit->GetDEPrimary(jndex+1);
      Int_t kndex ;
      Bool_t already = kFALSE ;
      for ( kndex = 0 ; kndex < fMulTrack ; kndex++ ) { //check if not already stored
	if ( newPrimary == primArray[kndex] ){
	  already = kTRUE ;
	  dEPrimArray[kndex] += dEPrimary; 
	  break ;
	}
      } // end of check
      if ( !already && (fMulTrack < fMaxTrack)) { // store it
	primArray[fMulTrack] = newPrimary ; 
	dEPrimArray[fMulTrack] = dEPrimary ; 
	fMulTrack++ ;
      } // store it
    } // all primaries in digit
  } // all digits

  Int_t *sortIdx = new Int_t[fMulTrack];
  TMath::Sort(fMulTrack,dEPrimArray,sortIdx); 
  for(index = 0; index < fMulTrack; index++) {
    fTracksList[index] = primArray[sortIdx[index]] ;    
    fDETracksList[index] = dEPrimArray[sortIdx[index]] ;
  }
  delete [] sortIdx;
  delete [] primArray ;
  delete [] dEPrimArray ;

}

//______________________________________________________________________________
void  AliEMCALRecPoint::EvalParents(TClonesArray * digits)
{
  // Constructs the list of parent particles (tracks) which have contributed to this RecPoint
 
  AliEMCALDigit * digit ;
  Int_t * parentArray = new Int_t[fMaxTrack] ;
  Float_t * dEParentArray = new Float_t[fMaxTrack] ;

  Int_t index ;  
  for ( index = 0 ; index < GetDigitsMultiplicity() ; index++ ) { // all digits
    if (fDigitsList[index] >= digits->GetEntries() || fDigitsList[index] < 0)
       AliError(Form("Trying to get invalid digit %d (idx in WriteRecPoint %d)",fDigitsList[index],index));
    digit = dynamic_cast<AliEMCALDigit *>(digits->At( fDigitsList[index] )) ; 
    Int_t nparents = digit->GetNiparent() ;
    if ( nparents == 0 ) continue ;

    Int_t jndex ;
    for ( jndex = 0 ; jndex < nparents ; jndex++ ) { // all primaries in digit
      if ( fMulParent > fMaxParent ) {
	fMulTrack = - 1 ;
	Error("EvalParents", "increase fMaxParent")  ;
	break ;
      }
      Int_t newParent = digit->GetIparent(jndex+1) ;
      Float_t newdEParent = digit->GetDEParent(jndex+1) ;
      Int_t kndex ;
      Bool_t already = kFALSE ;
      for ( kndex = 0 ; kndex < fMulParent ; kndex++ ) { //check if not already stored
	if ( newParent == parentArray[kndex] ){
	  dEParentArray[kndex] += newdEParent;
	  already = kTRUE ;
	  break ;
	}
      } // end of check
      if ( !already && (fMulParent < fMaxParent)) { // store it
	parentArray[fMulParent] = newParent ; 
	dEParentArray[fMulParent] = newdEParent ; 
	fMulParent++ ;
      } // store it
    } // all parents in digit
  } // all digits

  if (fMulParent>0) {
    Int_t *sortIdx = new Int_t[fMulParent];
    TMath::Sort(fMulParent,dEParentArray,sortIdx); 
    for(index = 0; index < fMulParent; index++) {
      fParentsList[index] = parentArray[sortIdx[index]] ;      
      fDEParentsList[index] = dEParentArray[sortIdx[index]] ;
    }
    delete [] sortIdx;
  }
 
  delete [] parentArray;
  delete [] dEParentArray;
}

//____________________________________________________________________________
void AliEMCALRecPoint::GetLocalPosition(TVector3 & lpos) const
{
  // returns the position of the cluster in the local reference system
  // of the sub-detector
  
  lpos = fLocPos;
}

//____________________________________________________________________________
void AliEMCALRecPoint::GetGlobalPosition(TVector3 & gpos) const
{
  // returns the position of the cluster in the global reference system of ALICE
  // These are now the Cartesian X, Y and Z
  //  cout<<" geom "<<geom<<endl;
  fGeomPtr->GetGlobal(fLocPos, gpos, fSuperModuleNumber);

}

//____________________________________________________________________________
void AliEMCALRecPoint::GetGlobalPosition(TVector3 & gpos, TMatrixF & gmat) const
{
  // returns the position of the cluster in the global reference system of ALICE
  // These are now the Cartesian X, Y and Z
  //  cout<<" geom "<<geom<<endl;

  //To be implemented
  fGeomPtr->GetGlobalEMCAL(this, gpos, gmat);

}

//_____________________________________________________________________________
void AliEMCALRecPoint::EvalLocal2TrackingCSTransform()
{
  //Evaluates local to "tracking" c.s. transformation (B.P.).
  //All evaluations should be completed before calling for this
  //function.                           
  //See ALICE PPR Chapter 5 p.18 for "tracking" c.s. definition,
  //or just ask Jouri Belikov. :) 

  SetVolumeId(AliGeomManager::LayerToVolUID(AliGeomManager::kEMCAL,GetSuperModuleNumber()));

  const TGeoHMatrix* tr2loc = GetTracking2LocalMatrix();
  if(!tr2loc) AliFatal(Form("No Tracking2LocalMatrix found."));

  Double_t lxyz[3] = {fLocPos.X(),fLocPos.Y(),fLocPos.Z()};
  Double_t txyz[3] = {0,0,0};

  tr2loc->MasterToLocal(lxyz,txyz);
  SetX(txyz[0]); SetY(txyz[1]); SetZ(txyz[2]);

  if(AliLog::GetGlobalDebugLevel()>0) {
    TVector3 gpos; TMatrixF gmat;
    GetGlobalPosition(gpos,gmat);
    Float_t gxyz[3];
    GetGlobalXYZ(gxyz);
    AliInfo(Form("lCS-->(%.3f,%.3f,%.3f), tCS-->(%.3f,%.3f,%.3f), gCS-->(%.3f,%.3f,%.3f),  gCScalc-\
->(%.3f,%.3f,%.3f), supermodule %d",
                 fLocPos.X(),fLocPos.Y(),fLocPos.Z(),
                 GetX(),GetY(),GetZ(),
                 gpos.X(),gpos.Y(),gpos.Z(),
                 gxyz[0],gxyz[1],gxyz[2],GetSuperModuleNumber()));
  }

}

//____________________________________________________________________________
Float_t AliEMCALRecPoint::GetMaximalEnergy(void) const
{
  // Finds the maximum energy in the cluster
  
  Float_t menergy = 0. ;

  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
 
    if(fEnergyList[iDigit] > menergy) 
      menergy = fEnergyList[iDigit] ;
  }
  return menergy ;
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::GetMultiplicityAtLevel(Float_t H) const
{
  // Calculates the multiplicity of digits with energy larger than H*energy 
  
  Int_t multipl   = 0 ;
  Int_t iDigit ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {

    if(fEnergyList[iDigit] > H * fAmp) 
      multipl++ ;
  }
  return multipl ;
}

//____________________________________________________________________________
Int_t  AliEMCALRecPoint::GetNumberOfLocalMax(AliEMCALDigit **  maxAt, Float_t * maxAtEnergy,
					       Float_t locMaxCut,TClonesArray * digits) const
{ 
  // Calculates the number of local maxima in the cluster using fLocalMaxCut as the minimum
  // energy difference between two local maxima

  AliEMCALDigit * digit ;
  AliEMCALDigit * digitN ;
  
  Int_t iDigitN ;
  Int_t iDigit ;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++)
    maxAt[iDigit] = (AliEMCALDigit*) digits->At(fDigitsList[iDigit])  ;
  
  for(iDigit = 0 ; iDigit < fMulDigit; iDigit++) {   
    if(maxAt[iDigit]) {
      digit = maxAt[iDigit] ;
          
      for(iDigitN = 0; iDigitN < fMulDigit; iDigitN++) {	
	digitN = (AliEMCALDigit *) digits->At(fDigitsList[iDigitN]) ; 
	
	if ( AreNeighbours(digit, digitN) ) {
	  if (fEnergyList[iDigit] > fEnergyList[iDigitN] ) {    
	    maxAt[iDigitN] = 0 ;
	    // but may be digit too is not local max ?
	    if(fEnergyList[iDigit] < fEnergyList[iDigitN] + locMaxCut) 
	      maxAt[iDigit] = 0 ;
	  }
	  else {
	    maxAt[iDigit] = 0 ;
	    // but may be digitN too is not local max ?
	    if(fEnergyList[iDigit] > fEnergyList[iDigitN] - locMaxCut) 
	      maxAt[iDigitN] = 0 ; 
	  } 
	} // if Areneighbours
      } // while digitN
    } // slot not empty
  } // while digit
  
  iDigitN = 0 ;
  for(iDigit = 0; iDigit < fMulDigit; iDigit++) { 
    if(maxAt[iDigit] ){
      maxAt[iDigitN] = maxAt[iDigit] ;
      maxAtEnergy[iDigitN] = fEnergyList[iDigit] ;
      iDigitN++ ; 
    }
  }
  return iDigitN ;
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::GetPrimaryIndex() const  
{
  // Get the primary track index in TreeK which deposits the most energy 
  // in Digits which forms RecPoint. 

  if (fMulTrack)
    return fTracksList[0];
  return -12345;
}

//____________________________________________________________________________
void AliEMCALRecPoint::EvalTime(TClonesArray * digits){
  // time is set to the time of the digit with the maximum energy

  Float_t maxE = 0;
  Int_t maxAt = 0;
  for(Int_t idig=0; idig < fMulDigit; idig++){
    if(fEnergyList[idig] > maxE){
      maxE = fEnergyList[idig] ;
      maxAt = idig;
    }
  }
  fTime = ((AliEMCALDigit*) digits->At(fDigitsList[maxAt]))->GetTime() ;
  
}

//______________________________________________________________________________
void AliEMCALRecPoint::Paint(Option_t *)
{
  // Paint this ALiRecPoint as a TMarker  with its current attributes
  
  TVector3 pos(0.,0.,0.)  ;
  GetLocalPosition(pos)   ;
  Coord_t x = pos.X()     ;
  Coord_t y = pos.Z()     ;
  Color_t markercolor = 1 ;
  Size_t  markersize = 1. ;
  Style_t markerstyle = 5 ;
  
  if (!gPad->IsBatch()) {
    gVirtualX->SetMarkerColor(markercolor) ;
    gVirtualX->SetMarkerSize (markersize)  ;
    gVirtualX->SetMarkerStyle(markerstyle) ;
  }
  gPad->SetAttMarkerPS(markercolor,markerstyle,markersize) ;
  gPad->PaintPolyMarker(1,&x,&y,"") ;
}

//_____________________________________________________________________
Double_t AliEMCALRecPoint::TmaxInCm(const Double_t e , const Int_t key)
{ 
  // e energy in GeV)
  // key  =  0(gamma, default)
  //     !=  0(electron)
  static Double_t ca = 4.82;  // shower max parameter - first guess; ca=TMath::Log(1000./8.07)
  static Double_t X0 = 1.23;  // radiation lenght (cm)
  static Double_t tmax = 0.;   // position of electromagnetic shower max in cm

  tmax = 0.0;
  if(e>0.1) {
    tmax = TMath::Log(e) + ca;
    if      (key==0) tmax += 0.5; 
    else             tmax -= 0.5;
    tmax *= X0; // convert to cm
  }
  return tmax;
}

//______________________________________________________________________________
Float_t AliEMCALRecPoint::EtaToTheta(Float_t arg) const
{
  //Converts Theta (Radians) to Eta(Radians)
  return (2.*TMath::ATan(TMath::Exp(-arg)));
}

//______________________________________________________________________________
Float_t AliEMCALRecPoint::ThetaToEta(Float_t arg) const
{
  //Converts Eta (Radians) to Theta(Radians)
  return (-1 * TMath::Log(TMath::Tan(0.5 * arg)));
}

//____________________________________________________________________________
void AliEMCALRecPoint::Print(Option_t *opt) const
{
  // Print the list of digits belonging to the cluster
  if(strlen(opt)==0) return;
  TString message ; 
  message  = "AliEMCALRecPoint:\n" ;
  message +=  " digits # = " ; 
  Info("Print", message.Data()) ; 

  Int_t iDigit;
  for(iDigit=0; iDigit<fMulDigit; iDigit++)
    printf(" %d ", fDigitsList[iDigit] ) ;  
  printf("\n");

  Info("Print", " Energies = ") ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf(" %f ", fEnergyList[iDigit] ) ;
  printf("\n");

  Info("Print", "\n Abs Ids = ") ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf(" %i ", fAbsIdList[iDigit] ) ;
  printf("\n");

  Info("Print", " Primaries  ") ;
  for(iDigit = 0;iDigit < fMulTrack; iDigit++)
    printf(" %d ", fTracksList[iDigit]) ;

  printf("\n Local x %6.2f y %7.2f z %7.1f \n", fLocPos[0], fLocPos[1], fLocPos[2]);

  message  = "       ClusterType     = %d" ;
  message += "       Multiplicity    = %d" ;
  message += "       Cluster Energy  = %f" ; 
  message += "       Core energy     = %f" ; 
  message += "       Core radius     = %f" ; 
  message += "       Number of primaries %d" ; 
  message += "       Stored at position %d" ; 
  Info("Print", message.Data(), fClusterType, fMulDigit, fAmp, fCoreEnergy, fCoreRadius, fMulTrack, GetIndexInList() ) ;  
}

//___________________________________________________________
Double_t  AliEMCALRecPoint::GetPointEnergy() const
{
  static double e;
  e=0.0;
  for(int ic=0; ic<GetMultiplicity(); ic++) e += double(fEnergyList[ic]);
  return e;
}
