#ifndef AliRICHHelix_h
#define AliRICHHelix_h

#include <TObject.h>              //base class
#include <TVector3.h>             //used extensively
#include "AliRICHParam.h"         //RichIntersect() 

class AliRICHHelix: public TObject
{
public:
  AliRICHHelix():TObject(),
    fX0(TVector3(0,0,0)),
    fP0(TVector3(0,0,0)),
    fX(TVector3(0,0,0)),
    fP(TVector3(0,0,0)),
    fLen(0),
    fQ(0),
    fBz(0),
    fPosRad(TVector2(0,0)),
    fPosPc(TVector2(0,0)),
    fPloc(TVector3(0,0,0))            {}
  AliRICHHelix(Double_t p,Double_t theta,Double_t phi,Double_t bz=0.2);//p [GeV], theta,phi [deg], Bz [Tesla];
  AliRICHHelix(const TVector3 &x0,const TVector3 &p0,Int_t q=1,Double_t b=0.2):TObject(),
    fX0(x0),
    fP0(p0),
    fX(x0),
    fP(p0),
    fLen(0),
    fQ(q),
    fBz(b),
    fPosRad(TVector2(0,0)),
    fPosPc(TVector2(0,0)),
    fPloc(TVector3(0,0,0))            {}
  virtual ~AliRICHHelix()                                                                                                                   {}        
           
         void     Draw          (const Option_t *opt=""              );                              //from TObject, draw helix
         void     Print         (const Option_t *opt=""              )const;                         //from TObject, print status
//private part++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  inline Bool_t   Intersection  (TVector3 planePnt,TVector3 planeNorm);                              //intersection with plane given by point and normal vector
  inline void     Propagate     (Double_t len                        );                              //propogate helix by given length along it
  inline Int_t    RichIntersect (AliRICHParam *pParam                );                              //search intersection with any RICH chamber 
         TVector2 PosRad        (                                    )const{return fPosRad;}         //intersection with radiator LORS
         TVector2 PosPc         (                                    )const{return fPosPc;}          //returns position of intersection with PC (local system)
         TVector3 Ploc          (                                    )const{return fPloc;}           //returns momentum at the position of intersection with radiator   
//         Double_t Length        ()                const{return fLen;}             //returns length of the track from initial point to RICH
         TVector3 X             (                                    )const{return fX;}
protected:
  TVector3 fX0;            //helix position in point of definition, [cm]    in MARS
  TVector3 fP0;            //helix momentum in point of definition, [GeV/c] in MARS
  TVector3 fX;             //helix position in point of interest,  [cm]    in MARS
  TVector3 fP;             //helix momentum in point of interest,  [GeV/c] in MARS
  Double_t fLen;           //helix length from point of definition to point of interest, [cm]    
  Int_t    fQ;             //sign of track charge (value not provided by current ESD)  
  Double_t fBz;            //magnetic field along z, [Tesla] 
  TVector2 fPosRad;        //helix intersection with radiator entrance, LORS [cm]
  TVector2 fPosPc;         //helix intersection with PC, LORS [cm]
  TVector3 fPloc;          //helix momentum, LORS [GeV/c]
  ClassDef(AliRICHHelix,0) //General helix
};//class AliRICHHelix
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHHelix::Propagate(Double_t len)
{
// Propogates the helix from inintial point by a given distance along helix. Assumes uniform magnetic field along z direction.  
// Arguments: len - distance to propagate by, [cm]
//   Returns: none    
  if(fBz==0){//no magnetic field->straight line
    fX=fX0+fP0.Unit()*len;   
  }else{
    const Double_t c = 0.00299792458;//this speed of light value provides that coordinates are in cm momentum in GeV/c
    Double_t a = -c*fBz*fQ;
     Double_t rho = a/fP0.Mag();
    fX.SetX( fX0.X()+fP0.X()*TMath::Sin(rho*len)/a-fP0.Y()*(1-TMath::Cos(rho*len))/a  );
    fX.SetY( fX0.Y()+fP0.Y()*TMath::Sin(rho*len)/a+fP0.X()*(1-TMath::Cos(rho*len))/a  ); 
    fX.SetZ( fX0.Z()+fP0.Z()*len/fP0.Mag()                                            );
    fP.SetX( fP0.X()*TMath::Cos(rho*len)-fP0.Y()*TMath::Sin(rho*len)                  );
    fP.SetY( fP0.Y()*TMath::Cos(rho*len)+fP0.X()*TMath::Sin(rho*len)                  );
    fP.SetZ( fP0.Z()                                                                  );
    fLen=len;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliRICHHelix::Intersection(TVector3 planePoint,TVector3 planeNorm)
{
// Finds point of intersection (if exists) of the helix with the plane. Stores result in fX and fP.   
// Arguments: planePoint,planeNorm - the plane defined by any plane's point and vector, normal to the plane
//   Returns: kTrue if helix intersects the plane, kFALSE otherwise.
  
  Double_t s=(planePoint-fX0)*planeNorm,dist=99999,distPrev=dist;//estimates initial distance to plane

  while(TMath::Abs(dist)>0.00001){
    Propagate(s);                        //calculates helix at the distance s from x0 ALONG the helix
    dist=(fX-planePoint)*planeNorm;      //distance between current helix position and plane
    if(TMath::Abs(dist) >= TMath::Abs(distPrev)) { return kFALSE;}
    distPrev=dist;
    s-=dist;
  }
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHHelix::RichIntersect(AliRICHParam *pParam)
{
// Searchs for intersection of this helix with all RICH chambers, returns chamber number or 0 if no intersection
// On exit fPosRad contain position of intersection in radiator LORS (cm)    
//         fPosPc  contains the same for photocathode  
  for(Int_t iCh=1;iCh<=AliRICHParam::kNch;iCh++){//chambers loop
    TVector3 norm =pParam->Lors2MarsVec(iCh,TVector3(0,0,1));
    if(Intersection(pParam->Lors2Mars(iCh,0,0,AliRICHParam::kRad),norm)){//there is intersection with radiator plane
      fPosRad=pParam->Mars2Lors(iCh,fX,AliRICHParam::kRad);//position on radiator plane
      if(pParam->IsAccepted(fPosRad)){//intersection within radiator (even if in dead zone)
        
        if(Intersection(pParam->Lors2Mars(iCh,0,0,AliRICHParam::kPc),norm)){//there is intersection with photocathode
          fPosPc=pParam->Mars2Lors(iCh,fX,AliRICHParam::kPc);//position on radiator plane
          if(pParam->IsAccepted(fPosPc)){//intersection within pc (even if in dead zone)
            
            fPloc=pParam->Mars2LorsVec(iCh,fP);//trasform p to local system
            return iCh;
          }//if inside PC
        }//if intersects PC
        
      }//if inside radiator
    }//if for radiator       
  }//chambers loop
  return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
