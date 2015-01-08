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
 
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TVector3.h>
#include <TMath.h>
#include "AliGeomManager.h"
#include "AliITSsegmentationUpgrade.h"
#include "AliLog.h"
#include <TFile.h>
//////////////////////////////////////////////////////
// Segmentation class for                           //
// ITS upgrade                                      //
//                                                  //
//////////////////////////////////////////////////////

ClassImp(AliITSsegmentationUpgrade)

//_____________________________________________________________________________

  AliITSsegmentationUpgrade::AliITSsegmentationUpgrade(): 
TObject(),
							  fCellSizeX(0),
							  fCellSizeZ(0),
 							  fNSectors(0),
							  fMinRadius(0),
							  fMaxRadius(0),
							  fHalfLength(0)
{ 

  AliDebug(2,"Default constructor is called");
  // Default constructor

  if(!gGeoManager) AliGeomManager::LoadGeometry("geometry.root");
  TFile *f=TFile::Open("Segmentation.root");
  TArrayD *x=0;
  TArrayD *z=0; 
  TArrayD *nS=0;

  if(!f){
    AliError("Segmentation not available");
    return;
  }else {
    x=(TArrayD*)f->Get("CellSizeX");  
    z=(TArrayD*)f->Get("CellSizeZ");
    nS=(TArrayD*)f->Get("nSectors");
  }
  f->Close(); 
    
  Int_t i=0;
  while(gGeoManager->GetVolume(Form("LayerSilicon%i",i))){
    TGeoVolume *vol = gGeoManager->GetVolume(Form("LayerSilicon%i",i));
    if(!vol) {
      AliInfo(Form("the volume %s has not been found... exiting!",Form("LayerSilicon%i",i)));
      return;
    }

    TGeoTube *shape = (TGeoTube*)vol->GetShape();
    if(!shape) {
      AliInfo(Form("the volume %s has not shape defined... exiting!",vol->GetName()));
      return;
    }   
   // AliInfo(Form(" volume index %i ", i));
    // setting the geometry parameters (needed for trasformations Global-Local)
    fMaxRadius.Set(i+1);   fMaxRadius.AddAt(shape->GetRmax(),i);
    fMinRadius.Set(i+1);   fMinRadius.AddAt(shape->GetRmin(),i);
    fHalfLength.Set(i+1); fHalfLength.AddAt(shape->GetDz(),i);
    fCellSizeX.Set(i+1);
    fCellSizeZ.Set(i+1);
    fCellSizeX.AddAt(x->At(i),i);
    fCellSizeZ.AddAt(z->At(i),i); 
    i++;  
  }
 fNSectors = (Int_t)(nS->At(0));

}
//_______________________________________________________________
AliITSsegmentationUpgrade::AliITSsegmentationUpgrade(TArrayD radii, TArrayD widths, TArrayD Length): 
TObject(),
												     fCellSizeX(0),
												     fCellSizeZ(0),
												     fNSectors(110),
												     fMinRadius(0),
												     fMaxRadius(0),
												     fHalfLength(0)
{
  for(Int_t i=0; i<radii.GetSize();i++){
    fMaxRadius.Set(i+1);
    fMaxRadius.AddAt(radii.At(i)+widths.At(i),i);
    fMinRadius.Set(i+1);
    fMinRadius.AddAt(radii.At(i),i);
    fHalfLength.Set(i+1);
    fHalfLength.AddAt(Length.At(i),i);
  }
}
//_____________________________________________________________________________
void AliITSsegmentationUpgrade::SetSegmentation(Int_t ilayer, Double_t xsize, Double_t zsize){
  
  // x/z size in microns
  if(fCellSizeX.GetSize()<ilayer+1) fCellSizeX.Set(ilayer+1);
  if(fCellSizeZ.GetSize()<ilayer+1) fCellSizeZ.Set(ilayer+1); 
  AliDebug(10,Form("xsize %f zsize %f ilayer %i", xsize,zsize,ilayer));
  fCellSizeX.AddAt(xsize,ilayer); 
  fCellSizeZ.AddAt(zsize,ilayer);  
  AliDebug(10,Form("fCellsizeX %f fCellSizeZ %f", fCellSizeX.At(ilayer), fCellSizeZ.At(ilayer)));       
}
//_____________________________________________________________________________
void AliITSsegmentationUpgrade::SetFullSegmentation(TArrayD xsize, TArrayD zsize){
 
  if(xsize.GetSize()!=zsize.GetSize())AliDebug(1,"Be Carefull Array Size Differ!!");
    
  if(xsize.GetSize()!=fMaxRadius.GetSize())AliDebug(10,Form("Be Carefull Segmentation Array (%i) and Geometry Array Differ (%i)!!",xsize.GetSize(),fMaxRadius.GetSize()));
   
  for(Int_t ilayer=0; ilayer<fMaxRadius.GetSize(); ilayer++)SetSegmentation(ilayer, xsize.At(ilayer),  zsize.At(ilayer));
}
//_________________________________________________________________________________
void AliITSsegmentationUpgrade::GetSegmentation(Int_t ilayer, Double_t &xsize, Double_t &zsize ) const{
  
  xsize = fCellSizeX.At(ilayer);
  zsize =fCellSizeZ.At(ilayer);
}
//_____________________________________________________________________________
Bool_t  AliITSsegmentationUpgrade::GlobalToDet(Int_t ilayer, Double_t x,Double_t y,Double_t z,Double_t &xl,Double_t &zl) const{
  

  if(TMath::Abs(z)>fHalfLength.At(ilayer)) return kFALSE;
  
  zl = z;
  Double_t xyproj= TMath::Sqrt(x*x+y*y);
  if(xyproj==0 || xyproj-0.001 > fMaxRadius.At(ilayer))return kFALSE;  
  Double_t alpha= TMath::ATan2(y,x);
  if(x>0 && y<0)alpha=alpha+2.*(TMath::Pi());
  if(x<0 && y<0)alpha=alpha+2.*(TMath::Pi());
  xl = xyproj*alpha;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliITSsegmentationUpgrade::DetToGlobal(Int_t ilayer, Double_t xl,Double_t zl,Double_t &x,Double_t &y, Double_t &z) const {
  
  if(TMath::Abs(z)>fHalfLength.At(ilayer)) return kFALSE;
  z = zl;
  Double_t alpha = xl/fMaxRadius.At(ilayer);

  
  x = fMaxRadius.At(ilayer)*TMath::Cos(alpha);
  y = fMaxRadius.At(ilayer)*TMath::Sin(alpha);

  return kTRUE;
}
//_______________________________________________________________________________
void AliITSsegmentationUpgrade::GetNpad(Int_t ilayer, Int_t &nx, Int_t &nz){
      
  if(fCellSizeX.At(ilayer)==0)AliDebug(1,"Attention! Check your X Segmentation!!");
  if(fCellSizeZ.At(ilayer)==0)AliDebug(1,"Attention! Check your Z Segmentation!!");
     
  nx=(Int_t)(2.*TMath::Pi()*fMaxRadius.At(ilayer)/fCellSizeX.At(ilayer));
  nz=(Int_t)(2.*fHalfLength.At(ilayer)/fCellSizeZ.At(ilayer));     

}
//________________________________________________________________________________
Int_t AliITSsegmentationUpgrade::GetNLayers(){

  if(!gGeoManager) TGeoManager::Import("geometry.root");

  Int_t nlay=0;

  for(Int_t ivol=0; ivol<gGeoManager->GetListOfVolumes()->GetEntries();ivol++){
    TString volname = gGeoManager->GetListOfVolumes()->At(ivol)->GetName();
    if(!volname.Contains("Silicon")) continue;

    nlay++;

  }

  return nlay;




}
//__________________________________________________________________________________
Bool_t AliITSsegmentationUpgrade::DetToGlobal(Int_t ilayer, Int_t module, Double_t xl,Double_t zl,Double_t &x,Double_t &y, Double_t &z) const {
if(module<0){ DetToGlobal(ilayer,xl,zl, x, y, z);
return kTRUE;
}

 if(TMath::Abs(z)>fHalfLength.At(ilayer)) return kFALSE;
  z = zl;
  Double_t alpha = (xl+module*(fMaxRadius.At(ilayer)*((TMath::Pi()*2)/fNSectors)))/fMaxRadius.At(ilayer);


  x = fMaxRadius.At(ilayer)*TMath::Cos(alpha);
  y = fMaxRadius.At(ilayer)*TMath::Sin(alpha);

return kTRUE;
}
//_______________________________________________________________________________________________
Bool_t AliITSsegmentationUpgrade::GlobalToDet(Int_t ilayer, Double_t x,Double_t y,Double_t z,Double_t &xl,Double_t &zl,Int_t &module) const {
 if(module<0){
  GlobalToDet(ilayer,x,y,z, xl, zl);
  return kTRUE;
 }

  if(TMath::Abs(z)>fHalfLength.At(ilayer)) return kFALSE;

  zl = z;
 
 Double_t xyproj= TMath::Sqrt(x*x+y*y);
  if(xyproj==0 || xyproj-0.001 > fMaxRadius.At(ilayer)){
return kFALSE;
  }
  Double_t alpha= TMath::ATan2(y,x);
  if(x>0 && y<0)alpha=alpha+2.*(TMath::Pi());
  if(x<0 && y<0)alpha=alpha+2.*(TMath::Pi());
  
  module = (Int_t)(alpha*(fNSectors/(TMath::Pi()*2)));
  Double_t alphaSect= module*((TMath::Pi()*2)/fNSectors);
  xl = xyproj*(alpha-alphaSect);//with respect the module
  //AliInfo(Form(" globtodet con modulo x,y,z = %f, %f, %f,  alpha  %f  module %i alpha del modulo  %f  xlocal %f xlocal rispettp al modulo %f ", x,y,z,alpha,module, alphaSect, xyproj*alpha, xl));
  return kTRUE;
}
//____________________________________________________________________________________________________
Bool_t AliITSsegmentationUpgrade::DetToPixID(Double_t xl, Double_t zl,Int_t layer, Int_t &nx, Int_t &nz) const {

   if(xl>fMaxRadius.At(layer)*(TMath::Pi()*2/fNSectors))return kFALSE;
   if(TMath::Abs(zl)>fHalfLength.At(layer)) return kFALSE;

    nx=(Int_t)(xl/fCellSizeX.At(layer));
    //printf(" half %f \n ", fHalfLength.At(layer));
    nz=(Int_t)((zl+fHalfLength.At(layer))/fCellSizeZ.At(layer));
  return kTRUE; 
}
//_________________________________________________________________________________________________________
Bool_t AliITSsegmentationUpgrade::DetToTrack2(Int_t layer,Int_t module, Double_t xl,Double_t zl, Double_t &ytr, Double_t &ztr)const {
     //sistema di tracciamento la x=0, y è la rphi rispetto al centro del modulo
     //sistema locale è la y=0 x è rphi rispetto all'inizio del modulo 

     //voglio andare dalle locali alle tracciamento
      Double_t x=0.;
      Double_t y=0.;
      Double_t z=0.;
//     printf( " det to track xl %f  zl %f \n ", xl,  zl); 
      ztr = zl;
      Bool_t ch=kFALSE;
      ch=DetToGlobal(layer,module,xl,zl,x,y,z);      
      Double_t radius =fMaxRadius.At(layer);
      Double_t phiGrad;
      phiGrad=TMath::ATan2(y,x);//cluster phi angle (rad)
      if (phiGrad<0) phiGrad+=TMath::TwoPi();//from 0 to 360
      else if (phiGrad>=TMath::TwoPi()) phiGrad-=TMath::TwoPi();//

      Double_t alpha= (module*(((TMath::Pi()*2)/fNSectors))+(((TMath::Pi()*2)/fNSectors)/2));//angle at the center of the ladder (rad)

      Double_t phiTrk= (phiGrad-alpha);// cluster angle in the rotated system (rad)

      ytr=radius*phiTrk; // tracking system coordinate: r*phi
 //     AliInfo(Form( " det : x %f |  track y %f,  zeta %f ", xl, ytr,ztr));
 return kTRUE;
}
//__________________________________________
Bool_t AliITSsegmentationUpgrade::DetToTrack(Int_t layer,Int_t module, Double_t xl,Double_t zl, Double_t &ytr, Double_t &ztr)const {
     //sistema di tracciamento la x=0, y è la rphi rispetto al centro del modulo
     //sistema locale è la y=0 x è rphi rispetto all'inizio del modulo

     //voglio andare dalle locali alle tracciamento
      ztr = zl;
      Double_t radius =fMaxRadius.At(layer);
      Double_t alphaloc = module*(TMath::Pi()*2)/fNSectors; //alpha di riferimento per il modulo locale in rad     
      Double_t xlrif = alphaloc*radius;//xlocale di riferimento per questo modulo

      Double_t xlvera = xlrif + xl; // x locale vera rispetto a tutto il cilindro, per sapere a che angolo sto effettivamente

      Double_t alphaGlob = xlvera/radius; //in rad 

//da alpha vera posso ottenereinfo di tracciamento 

      Double_t alpha_track= (module*(((TMath::Pi()*2)/fNSectors))+(((TMath::Pi()*2)/fNSectors)/2));//angle at the center of the ladder (rad)

      Double_t phiTrk= (alphaGlob-alpha_track);// cluster angle in the rotated system (rad)

      ytr=radius*phiTrk; // tracking system coordinate: r*phi
//      AliInfo(Form( " det : x %f |  track y %f,  zeta %f ", xl, ytr,ztr));
 return kTRUE;
}
//_______________________
Int_t AliITSsegmentationUpgrade::GetModule(Double_t x,Double_t y)const{
  

  Double_t alpha= TMath::ATan2(y,x);
  if(x>0 && y<0)alpha=alpha+2.*(TMath::Pi());
  if(x<0 && y<0)alpha=alpha+2.*(TMath::Pi());

return GetModule(alpha);
}
//_______________________
Int_t AliITSsegmentationUpgrade::GetModule(Double_t phi) const {
// 
// sector corresponding to that phi [0,2Pi]
//

Int_t module = -1;
if(phi<0 || phi > TMath::TwoPi()) {
  AliError(Form("Phi is not within the acceptable range [0,2Pi] %f ",phi));
  }
else {
  module = (Int_t) (phi*((Double_t)fNSectors/(TMath::TwoPi())));
 } 

return module;

}
//_______________________
Double_t AliITSsegmentationUpgrade::GetAlpha(Int_t module)const{
//
// angle of the sector center
//
Double_t angle = TMath::DegToRad()*(360./fNSectors);
return module*angle+0.5*angle;

}



