#ifndef AliRICH_h
#define AliRICH_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDetector.h>  //inheritance
#include <TClonesArray.h> 
#include <TObjArray.h>
#include <TVector.h>
#include <TVector3.h>

#include "AliRICHParam.h"
#include "AliRICHCluster.h"
#include "AliRICHHit.h"

//__________________AliRICH_________________________________________________________________________
class AliESD;

class AliRICH : public AliDetector 
{
public:
//ctor & dtor    
            AliRICH();                                            
            AliRICH(const char *name, const char *title);
            AliRICH(const AliRICH& RICH):AliDetector(RICH) {;}  //copy ctor 
  virtual  ~AliRICH();                                            
          
  AliRICH&  operator=(const AliRICH&)                 {return *this;}
//framework part  
  virtual Int_t         IsVersion()                           const =0;                                  //interface from         
  virtual void          StepManager()                               =0;                                  //interface from AliMC
  virtual void          SetTreeAddress();                                                                //interface from AliLoader
  virtual void          MakeBranch(Option_t *opt=" ");                                                   //interface from AliLoader
  virtual void          CreateMaterials();                                                               //interface from AliMC
  virtual void          CreateGeometry();                                                                //interface from AliMC
  virtual void          BuildGeometry();                                                                 //interface 
  virtual void          Print(Option_t *option="")               const;       //prints current RICH status
//private part  
          void          GeomPadPanelFrame()const;                                                        //defines PPF geometry
          void          GeomAmpGap()       const;                                                        //defines gap geometry + anod wires
          void          GeomRadiators()    const;                                                        //defines radiators geometry
          void          GeomSandBox()      const;                                                        //defines sandbox geometry
          void          GeomRadioSrc()     const;                                                        //defines radio source geometry
          void          GeomAerogel()      const;                                                        //defines aerogel geometry
  static  Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);       //deals with Fresnel absorption
  
         AliRICHHit*   Hit           (Int_t tid                                              )const;                                                 //first hit of given TID 
  inline void          HitAdd        (Int_t c,Int_t tid,TVector3 in,TVector3 out,Double_t e=0);                                                      //add new hit
  inline void          HitsCreate    (                                                       );                                                      //create hits container
         void          HitsPrint     (Int_t iEvent=0                                         )const;                                                 //prints hits
         void          HitsQA        (Double_t cut=0,Double_t cutele=0,Double_t cutR=999);
            
         TClonesArray* SDigits       (                                                       )const{return fSdigits;}                                //pointer to sdigits list 
  inline void          SDigitAdd     (Int_t c,TVector pad,Double_t q,Int_t pid,Int_t tid     )     ;                                                 //add new sdigit
  inline void          SDigitsCreate (                                                       )     ;                                                 //create sdigits container
         void          SDigitsReset  (                                                       )     {fNsdigits=0;  if(fSdigits)  fSdigits ->Clear();} //clean a list of sdigits                                
         void          SDigitsPrint  (Int_t iEvent=0                                         )const;                                                 //prints sdigits 
         
    using AliDetector::Digits;  
         TClonesArray* Digits        (Int_t iC                                               )const{return fDigs ? (TClonesArray *)fDigs->At(iC-1):0;}
  inline void          DigitAdd      (Int_t c,TVector pad,int q,int cfm,int *tid             )     ;                                                 //add new digit
  inline void          DigitAdd      (AliRICHDigit &dif                                      )     ;                                                 //add new digit
  inline void          DigitsCreate  (                                                       )     ;                                                 //create digits
         void          DigitsReset   (                                                       )     {if(fDigs)for(int i=0;i<kNchambers;i++){fDigs->At(i)->Clear();fNdigs[i]=0;}} //virtual
         void          DigitsPrint   (Int_t iEvent=0                                         )const;                                                 //prints digits
          
         TClonesArray* Clusters      (Int_t iC                                               )const{if(fClus)  return (TClonesArray *)fClus->At(iC-1);else return 0;}
  inline void          ClusterAdd    (AliRICHCluster &cl                                     )     ;                                                  //add new cluster                        
  inline void          ClustersCreate(                                                       )     ;                                                  //create clusters container
         void          ClustersReset (                                                       )     {if(fClus)for(int i=0;i<kNchambers;i++){fClus ->At(i)->Clear();fNclus[i]=0;}}
         void          ClustersPrint (Int_t iEvent=0                                         )const;                        //prints a list of clusters for a given event

         void          OccupancyPrint(Int_t iEvent=-1                                        )const;
         void          SummaryOfEvent(Int_t iEvent=0                                         )const;
         
  AliRICHChamber* C(Int_t iC)           const{return fParam->C(iC);}   //provides pointer to a given chamber
  AliRICHParam*   P()                   const{return fParam;}          //provides pointer to a RICH params
  AliRICH*        R()                        {return this;}             //provides pointer to RICH main object
  TVector         Counters()            const{return fCounters;}        //provides a set of counters
  void            ControlPlots();                                       //creates ~/RCP.root with a set of QA plots
  void            Display()const; //Display event
  void            DisplayEvent(Int_t,Int_t)const; //Display event
  static Int_t     Nparticles(Int_t iPid,Int_t iEventN,AliRunLoader *pRunLoader); //counts total number of particle with iPid
  void            PrintTracks  (Int_t iEvent=0);                        //prints a list of tracks for a given event
  void            ReadESD(Int_t iEventN, Int_t iChamber)const;
  void            DrawRing(TVector3 entrance,TVector3 vectorTrack,Double_t thetaCer)const;

protected:  
  enum                  EMedia {kAir=1,kRoha,kSiO2,kC6F14,kCH4,kCsI,kGridCu,kOpSiO2,kGap,kAl,kGlass,kCu,kW,kSteel,kPerpex,kSr90,kMylar,kGel,kReflector};
  enum                  ECounters {kStepManager=0,kCerProdTot,kCerProdRad,kCerKillTot,kCerKillRad,kCerKillRef,kEleProdTot};
  AliRICHParam         *fParam;                    //main RICH parametrization     
                                                   //fHits and fDigits belong to AliDetector
  TClonesArray         *fSdigits;                  //! list of sdigits  
  Int_t                 fNsdigits;                 //! current number of sdigits
  
  TObjArray            *fDigs;                     //! each chamber holds it's one lists of digits
  Int_t                 fNdigs[7];                 //! array of current numbers of digits
  
  TObjArray            *fClus;                     //! each chamber holds it's one lists of clusters 
  Int_t                 fNclus[7];                 //! array of current numbers of raw clusters
  
  TVector               fCounters;                 //Particle history counters, explanation in StepManager() 
  
  ClassDef(AliRICH,9)                              //Main RICH class 
};//class AliRICH  

//__________________________________________________________________________________________________
void AliRICH::HitsCreate()
{
  if(fHits) return;
  AliDebug(1,"creating hits container.");
  fHits=new TClonesArray("AliRICHHit",10000);   fNhits=0;
}
//__________________________________________________________________________________________________
void AliRICH::HitAdd(Int_t c,Int_t tid,TVector3 i3,TVector3 o3,Double_t eloss)
{
//add new RICH hit to the list of hits  
  TClonesArray &tmp=*fHits;
  new(tmp[fNhits++])AliRICHHit(c,tid,i3,o3,eloss);
}//AddHit()
//__________________________________________________________________________________________________
void AliRICH::SDigitsCreate()
{
  if(fSdigits) return;
  AliDebug(1,"creating sdigits container.");
  fSdigits=new TClonesArray("AliRICHDigit",10000); fNsdigits=0;
}
//__________________________________________________________________________________________________
void AliRICH::SDigitAdd(Int_t c,TVector pad,Double_t q,Int_t pid,Int_t tid) 
{ 
  Int_t cfm;  
  switch(pid){
    case 50000050: cfm=1000000;break;//cerenkov
    case 50000051: cfm=1000;   break;//feedback
    default:       cfm=1;      break;//mip
  }   
  TClonesArray &tmp=*fSdigits; new(tmp[fNsdigits++])AliRICHDigit(c,pad,q,cfm,tid,-1,-1);
}
//__________________________________________________________________________________________________
void AliRICH::DigitsCreate()
{
  if(fDigs) return;
  AliDebug(1,"creating digits containers.");
  fDigs = new TObjArray(kNchambers);  
  for(Int_t i=0;i<kNchambers;i++) {fDigs->AddAt(new TClonesArray("AliRICHDigit",10000), i); fNdigs[i]=0;}
}
//__________________________________________________________________________________________________
void AliRICH::DigitAdd(AliRICHDigit &dig)
{
//special for digit formed from raw  
  TClonesArray &tmp=*((TClonesArray*)fDigs->At(dig.Chamber()-1));
  new(tmp[fNdigs[dig.Chamber()-1]++])AliRICHDigit(dig);
}    
//__________________________________________________________________________________________________
void AliRICH::DigitAdd(int c,TVector pad,int q,int cfm,int *tid)
{
  TClonesArray &tmp=*((TClonesArray*)fDigs->At(c-1));
  new(tmp[fNdigs[c-1]++])AliRICHDigit(c,pad,q,cfm,tid[0],tid[1],tid[2]);
}    
//__________________________________________________________________________________________________
void AliRICH::ClustersCreate()
{
  if(fClus) return;
  AliDebug(1,"creating clusters containers.");
  fClus = new TObjArray(kNchambers);  
  for(Int_t i=0;i<kNchambers;i++) {fClus->AddAt(new TClonesArray("AliRICHCluster",10000), i); fNclus[i]=0;}
}
//__________________________________________________________________________________________________
void AliRICH::ClusterAdd(AliRICHCluster &cl)                     
{
  Int_t c=cl.C()-1;TClonesArray &tmp=*((TClonesArray*)fClus->At(c));
  new(tmp[fNclus[c]++])AliRICHCluster(cl);
}
//__________________________________________________________________________________________________
#endif//#ifndef AliRICH_h
