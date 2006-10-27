#ifndef AliRICHHelix_h
#define AliRICHHelix_h

#include <TObject.h>              //base class
#include <TVector3.h>             //used extensively

class AliRICHHelix: public TObject
{
public:
  AliRICHHelix(                                                            ):TObject(),fX0(TVector3(0,0,0)),fP0(TVector3(0,0,0)),fQ(0),fBz(0 )  {}
  AliRICHHelix(const TVector3 &x0,const TVector3 &p0,Int_t q=1,Double_t b=2):TObject(),fX0(x0             ),fP0(p0             ),fQ(q),fBz(b )  {}
  AliRICHHelix(Double_t p,Double_t theta,Double_t phi,Double_t bz=2        ):TObject(),fX0(TVector3(0,0,0)),fP0(TVector3(0,0,0)),fQ(0),fBz(bz)  
                   {fP0.SetMagThetaPhi(p,theta*TMath::DegToRad(),phi*TMath::DegToRad());} //p [GeV], theta,phi [deg], Bz [Tesla];
  virtual ~AliRICHHelix()                                                                                                                   {}        
           
         void     Draw          (const Option_t *opt=""              );                              //from TObject, draw helix
         void     Print    (const Option_t *opt=""              )const;                         //from TObject, print status
//private part++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  inline Bool_t   Intersect(            TVector3 &pnt,TVector3 &norm );          //intersection with plane given by point and normal vector
  inline void     Propagate(Float_t len,TVector3 &x,  TVector3 &p    );          //propogate helix by given length along it            
protected:
  TVector3 fX0;            //helix position in point of definition, [cm]    in MARS
  TVector3 fP0;            //helix momentum in point of definition, [GeV/c] in MARS
  Int_t    fQ;             //sign of track charge  
  Float_t  fBz;            //magnetic field along z, [kGaus] 
  ClassDef(AliRICHHelix,0) //General helix
};//class AliRICHHelix
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHHelix::Propagate(Float_t len,TVector3 &x,TVector3 &p)
{
// Propogates the helix from inintial point by a given distance along helix. Assumes uniform magnetic field along z direction.  
// Arguments: len - distance to propagate by, [cm]
//   Returns: none    
  if(fBz==0){//no magnetic field->straight line
    x=fX0+fP0.Unit()*len;   
    p=fP0;
  }else{
    const Float_t c = 0.000299792458;//this speed of light value provides that coordinates are in [cm] momentum in [GeV/c] field in [kGaus]
    Float_t a = -c*fBz*fQ;
    Float_t rho = a/fP0.Mag();
    x.SetX( fX0.X()+fP0.X()*TMath::Sin(rho*len)/a-fP0.Y()*(1-TMath::Cos(rho*len))/a  );
    x.SetY( fX0.Y()+fP0.Y()*TMath::Sin(rho*len)/a+fP0.X()*(1-TMath::Cos(rho*len))/a  ); 
    x.SetZ( fX0.Z()+fP0.Z()*len/fP0.Mag()                                            );
    x.SetX( fP0.X()*TMath::Cos(rho*len)-fP0.Y()*TMath::Sin(rho*len)                  );
    p.SetY( fP0.Y()*TMath::Cos(rho*len)+fP0.X()*TMath::Sin(rho*len)                  );
    p.SetZ( fP0.Z()                                                                  );
  }
}//Propagate()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliRICHHelix::Intersect(TVector3 &pnt,TVector3 &norm)
{
// Finds point of intersection (if exists) of the helix with the plane. Stores result in fX and fP.   
// Arguments: pnt  - arbitrary point of the plane, [cm] in MARS
//            norm - vector, normal to the plane, [cm] in MARS
//   Returns:      - kTrue if helix intersects the plane, kFALSE otherwise.
//                 - pnt contains the point of intersection, [cm] in MARS  
  TVector3 x,p;                                                    //current helix position and momentum
  Double_t s=(pnt-fX0)*norm,dist=99999,distPrev=dist;              //estimates initial distance to plane
  while(TMath::Abs(dist)>0.00001){                                 //loop while the distance is less then precision   
    Propagate(s,x,p);                                              //calculates helix at the distance s from x0 ALONG the helix
    dist=(x-pnt)*norm;                                             //distance between current helix position and plane
    if(TMath::Abs(dist) >= TMath::Abs(distPrev)) { return kFALSE;} //if distance increases then no intersection 
    distPrev=dist;
    s-=dist;
  }  
  norm=p;
  pnt=x;                                    
  return kTRUE;
}//Intersect()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
