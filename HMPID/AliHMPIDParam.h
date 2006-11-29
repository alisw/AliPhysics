#ifndef AliHMPIDParam_h
#define AliHMPIDParam_h

#include <TNamed.h>        //base class
#include <TGeoManager.h>   //Instance()
#include <TVector3.h>      //Lors2Mars() Mars2Lors()
 
static const int kCerenkov=50000050;  //??? go to something more general like TPDGCode
static const int kFeedback=50000051;  //??? go to something more general like TPDGCode

// Class providing all the needed parametrised information
// to construct the geometry, to define segmentation and to provide response model
// In future will also provide all the staff needed for alignment and calibration

class AliHMPIDParam :public TNamed  
{
public:
//ctor&dtor    
  virtual        ~AliHMPIDParam()                                    {for(Int_t i=0;i<7;i++) delete fM[i]; delete fgInstance; fgInstance=0;}
         void     Print(Option_t *opt="") const;                                         //print current parametrization
  static inline AliHMPIDParam* Instance();                                //pointer to AliHMPIDParam singleton
  
                Double_t   MeanIdxRad              () {return 1.29204;}//???????????
                Double_t   MeanIdxWin              () {return 1.57819;}//???????????
  static        Int_t      Stack(Int_t evt=-1,Int_t tid=-1);              //Print stack info for event and tid
  static        Int_t      StackCount(Int_t pid,Int_t evt);               //Counts stack particles of given sort in given event  
//trasformation methodes
  void     Lors2Mars   (Int_t c,Float_t x,Float_t y,Double_t *m,Int_t pl=kPc)const{Double_t z=0; switch(pl){case kPc:z=8.0;break; case kAnod:z=7.806;break; case kRad:z=-1.25; break;}   Double_t l[3]={x-fX,y-fY,z};  fM[c]->LocalToMaster(l,m); }    
  TVector3 Lors2Mars   (Int_t c,Float_t x,Float_t y,            Int_t pl=kPc)const{Double_t m[3];Lors2Mars(c,x,y,m,pl); return TVector3(m);    }//MRS->LRS  
  void     Mars2Lors   (Int_t c,Double_t *m,Float_t &x,Float_t &y           )const{Double_t l[3];fM[c]->MasterToLocal(m,l);x=l[0]+fX;y=l[1]+fY;}//MRS->LRS
  void     Mars2LorsVec(Int_t c,Double_t *m,Float_t &th,Float_t &ph         )const{Double_t l[3]; fM[c]->MasterToLocalVect(m,l); Float_t pt=TMath::Sqrt(l[0]*l[0]+l[1]*l[1]); th=TMath::ATan(l[3]/pt); ph=TMath::ATan(l[0]/pt);}    
  TVector3 Norm        (Int_t c                                             )const{Double_t n[3]; Norm(c,n); return TVector3(n);               }//norm 
  void     Norm        (Int_t c,Double_t *n                                 )const{Double_t l[3]={0,0,1};fM[c]->LocalToMasterVect(l,n);        }//norm

  enum EPlaneId {kPc,kRad,kAnod};            //3 planes in chamber 
protected:
         AliHMPIDParam();             //default ctor is protected to enforce it to be singleton
  static AliHMPIDParam *fgInstance;   //static pointer  to instance of AliHMPIDParam singleton
  TGeoHMatrix *fM[7];                //poiners to matrices defining HMPID chambers rotations-translations
  Float_t fX;                        //x shift of LORS with respect to rotated MARS 
  Float_t fY;                        //y shift of LORS with respect to rotated MARS   
  ClassDef(AliHMPIDParam,0)           //HMPID main parameters class
};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDParam* AliHMPIDParam::Instance()
{
// Return pointer to the AliHMPIDParam singleton. 
// Arguments: none
//   Returns: pointer to the instance of AliHMPIDParam or 0 if no geometry       
  if(!fgInstance)
    if(gGeoManager) new AliHMPIDParam; 
    else            Printf("AliHMPIDParam> Error:: No geometry defined!");
  return fgInstance;  
}//Instance()    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif
