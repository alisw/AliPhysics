#ifndef AliRICH_h
#define AliRICH_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDetector.h>  //base class
#include <TClonesArray.h>  
#include <TObjArray.h>
#include <TVector.h>
#include <TVector3.h>

#include "AliRICHCluster.h"
#include "AliRICHHit.h"

//__________________AliRICH_________________________________________________________________________
class AliESD;

class AliRICH : public AliDetector //TObject-TNamed-AliModule-AliDetector-AliRICH
{
public:
//ctor & dtor    
            AliRICH();                                            
            AliRICH(const char *name, const char *title);
            AliRICH(const AliRICH& RICH                ):AliDetector(RICH),fSdig(0),fSdigCnt(0),fDig(0),fClu(0),fCounters(0),fNcham(0) {;}  //copy ctor not implemented
  virtual  ~AliRICH();                                            
          
  AliRICH&  operator=(const AliRICH&)                 {return *this;}
//framework part  
          void  BuildGeometry   (                   );            //from AliModule invoked from AliMC::InitGeometry() to build geometry for event display
  virtual void  CreateMaterials (                   )=0;          //from AliModule invoked from AliMC::ConstructGeometry() to define detector materials
  virtual void  CreateGeometry  (                   )=0;          //from AliModule invoked from AliMC::ConstructGeometry() to build detector for simulation
  virtual Int_t IsVersion       (                   )const=0;     //from AliModule not used        
  virtual void  Init            (                   )=0;          //from AliModule invoked from AliMC::InitGeometry() after CreateGeometry() to do VolID initialization
          void  MakeBranch      (Option_t *opt=""   );            //from AliModule invokde from AliRun::Tree2Tree() to make requested RICH branch
//  virtual void  Print           (const Option_t *opt)const=0;     //from TObject
  virtual void  StepManager     (                   )=0;          //from AliModule invoked from AliMC
          void  SetTreeAddress  (                   );            //from AliModule invoked from AliRun::GetEvent(), AliLoader::SetTAddrInDet()
//private part  
  void    HitAdd    (Int_t c,Int_t tid,Int_t pid,TVector3 in,TVector3 out,Double_t e=0        ){new((*fHits)[fNhits++])AliRICHHit(c,tid,in,out,e,pid);          } 
  void    HitAdd    (Int_t tid,Double_t e,Int_t pad,Double_t x,Double_t y,Double_t z,Int_t pid){new((*fHits)[fNhits++])AliRICHHit(tid,e,pad,x,y,z,pid);         } 
  void    HitCreate (                                                                         ){if(fHits)return; fHits=new TClonesArray("AliRICHHit"); fNhits=0;}
  void    HitPrint  (Int_t iEvent=0                                                           )const;                                                 
  void    HitQA     (Double_t cut=0,Double_t cutele=0,Double_t cutR=999                       );
            
  TClonesArray* SDigs         (                                                       )const{return fSdig;}                                  //pointer to sdigits list 
  inline void          SDigAdd       (Int_t c,TVector pad,Double_t q,Int_t pid,Int_t tid     )     ;                                                 //add new sdigit
  void   SDigCreate (                                                       )   {if(fSdig) return;  fSdig=new TClonesArray("AliRICHDigit"); fSdigCnt=0;}
  void   SDigReset  (                                                       )   {fSdigCnt=0;  if(fSdig)  fSdig ->Clear();}     //clean a list of sdigits                                
  void   SDigPrint  (Int_t iEvent=0                                         )const;                                                 //prints sdigits 
         
         TClonesArray* Digs          (Int_t iC                                               )const{return fDig ? (TClonesArray *)fDig->At(iC-1):0;}
  inline void          DigAdd        (Int_t c,TVector pad,int q,int cfm,int *tid             )     ;                                                 //add new digit
  inline void          DigAdd        (AliRICHDigit &dif                                      )     ;                                                 //add new digit
  inline void          DigCreate     (                                                       )     ;                                                 //create digits
         void          DigReset      (                                                       )     {if(fDig)for(int i=0;i<fNcham;i++){fDig->At(i)->Clear();fDigCnt[i]=0;}} 
         void          DigPrint      (Int_t iEvent=0                                         )const;                                                 //prints digits
          
         TClonesArray* Clus          (Int_t iC                                               )const{return fClu ? (TClonesArray *)fClu->At(iC-1):0;}
  inline void          CluCreate     (                                                       )     ;                                                  //create clusters container
         void          CluReset      (                                                       )     {if(fClu)for(int i=0;i<fNcham;i++){fClu->At(i)->Clear();fCluCnt[i]=0;}}
         void          CluPrint      (Int_t iEvent=0                                         )const;                        //prints a list of clusters for a given event
         
         void     Display      (                                                        )const;                    //Display event
         void     DisplayEvent (Int_t,Int_t                                             )const;                    //Display event
         void     DrawRing     (TVector3 entrance,TVector3 vectorTrack,Double_t thetaCer)const;
         void     OccupancyPrint(Int_t iEvent=-1                                        )const;                    //print chambers occupancy 
         void     ReadESD      (Int_t iEventN, Int_t iChamber                           )const;                    
         void     SummaryOfEvent(Int_t iEvent=0                                         )const;

protected:  
  TClonesArray         *fSdig;                    //! list of sdigits  
  Int_t                 fSdigCnt;                 //! current number of sdigits
  
  TObjArray            *fDig;                     //! each chamber holds it's one list of digits
  Int_t                 fDigCnt[7];               //! array of current numbers of digits
  
  TObjArray            *fClu;                     //! each chamber holds it's one list of clusters 
  Int_t                 fCluCnt[7];               //! array of current numbers of raw clusters
  
  TVector               fCounters;                // Particle history counters, explanation in StepManager() 
  Int_t                 fNcham;                   // Number of RICH chambers during simulation
  
  ClassDef(AliRICH,11)                            //Main RICH class 
};//class AliRICH  

//__________________________________________________________________________________________________
void AliRICH::SDigAdd(Int_t c,TVector pad,Double_t q,Int_t pid,Int_t tid) 
{ 
  Int_t cfm;  
  switch(pid){
    case 50000050: cfm=1000000;break;//cerenkov
    case 50000051: cfm=1000;   break;//feedback
    default:       cfm=1;      break;//mip
  }   
  new((*fSdig)[fSdigCnt++])AliRICHDigit(c,pad,q,cfm,tid,-1,-1);
}
//__________________________________________________________________________________________________
void AliRICH::DigCreate()
{
  if(fDig) return;
  fDig = new TObjArray(fNcham);  
  for(Int_t i=0;i<fNcham;i++) {fDig->AddAt(new TClonesArray("AliRICHDigit"), i); fDigCnt[i]=0;}
}
//__________________________________________________________________________________________________
void AliRICH::DigAdd(AliRICHDigit &dig)
{
//special for digit formed from raw  
  TClonesArray &tmp=*((TClonesArray*)fDig->At(dig.C()-1));
  new(tmp[fDigCnt[dig.C()-1]++])AliRICHDigit(dig);
}    
//__________________________________________________________________________________________________
void AliRICH::DigAdd(int c,TVector pad,int q,int cfm,int *tid)
{
  TClonesArray &tmp=*((TClonesArray*)fDig->At(c-1));
  new(tmp[fDigCnt[c-1]++])AliRICHDigit(c,pad,q,cfm,tid[0],tid[1],tid[2]);
}    
//__________________________________________________________________________________________________
void AliRICH::CluCreate()
{
  if(fClu) return;
  fClu = new TObjArray(fNcham);  
  for(Int_t i=0;i<fNcham;i++) {fClu->AddAt(new TClonesArray("AliRICHCluster"), i); fCluCnt[i]=0;}
}
#endif//#ifndef AliRICH_h
