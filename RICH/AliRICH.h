#ifndef AliRICH_h
#define AliRICH_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDetector.h>  //base class
#include <TClonesArray.h> //XxxCreate() 
#include <TObjArray.h>    //fClu field



class AliRICH : public AliDetector //TObject-TNamed-AliModule-AliDetector-AliRICH
{
public:
//ctor & dtor    
            AliRICH(const char *nm,const char *ttl);                                              //named ctor
            AliRICH(                              ):AliDetector(    ),fSdi(0),fDig(0),fClu(0) {}  //default ctor          
  virtual  ~AliRICH();                                            
//framework part  
          void  BuildGeometry   (                ) {}          //from AliModule invoked from AliMC::InitGeometry() to build geometry for old event display
  virtual void  CreateMaterials (                )=0;          //from AliModule invoked from AliMC::ConstructGeometry() to define detector materials
  virtual void  CreateGeometry  (                )=0;          //from AliModule invoked from AliMC::ConstructGeometry() to build detector for simulation
  virtual Int_t IsVersion       (                )const=0;     //from AliModule not used        
  virtual void  Init            (                )=0;          //from AliModule invoked from AliMC::InitGeometry() after CreateGeometry() to do VolID initialization
          void  MakeBranch      (Option_t *opt="");            //from AliModule invokde from AliRun::Tree2Tree() to make requested RICH branch
          void  SetTreeAddress  (                );            //from AliModule invoked from AliRun::GetEvent(), AliLoader::SetTAddrInDet()
  virtual void  StepManager     (                )=0;          //from AliModule invoked from AliMC
//private part +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  void          HitCreate(         )     {if(fHits)return; fHits=new TClonesArray("AliRICHHit"); fNhits=0;      }//create hits list
  void          HitPrint (Int_t evt)const;                                                                       //print hits list
              
  TClonesArray* SdiLst   (         )const{return fSdi;                                                          }//get sdigits list 
  void          SdiCreate(         )     {if(fSdi)return; fSdi=new TClonesArray("AliRICHDigit");                }//create sdigits list
  void          SdiReset (         )     {if(fSdi)  fSdi ->Clear();                                             }//clean sdigits list
  void          SdiPrint (Int_t evt)const;                                                                       //print sdigits 
         
  TObjArray*    DigLst   (         )const{return fDig;                                                          }//get digits list for all chambers
  TClonesArray* DigLst   (Int_t c  )const{return fDig ? (TClonesArray *)fDig->At(c):0;                          }//get digits list for chamber
  void          DigCreate(         )     {fDig=new TObjArray(7);for(Int_t i=0;i<7;i++)fDig->AddAt(new TClonesArray("AliRICHDigit"),i);}//create digits list
  void          DigReset (         )     {if(fDig)for(int i=0;i<7;i++)fDig->At(i)->Clear();                     }//clean digits list 
  void          DigPrint (Int_t evt)const;                                                                       //print digits
          
  TClonesArray* CluLst   (Int_t c  )const{return fClu ? (TClonesArray *)fClu->At(c):0;                          }//get clusters list for chamber
  inline void   CluCreate(         )     {fClu=new TObjArray(7); for(Int_t i=0;i<7;i++)fClu->AddAt(new TClonesArray("AliRICHCluster"),i);}//create clusters list
         void   CluReset (         )     {if(fClu)for(int i=0;i<7;i++)fClu->At(i)->Clear();                     }//clean clusters list
         void   CluPrint (Int_t evt)const;                                                                       //print clusters list
         
  void          OccupancyPrint(Int_t evt=-1);                    //print chambers occupancy 
  void          SummaryOfEvent(Int_t evt=0)const;

protected:  
  TClonesArray         *fSdi;                     //! list of sdigits  
  TObjArray            *fDig;                     //! each chamber holds it's one list of digits
  TObjArray            *fClu;                     //! each chamber holds it's one list of clusters 
  
 private:
  AliRICH(const AliRICH &rich           );
  AliRICH&  operator=(const AliRICH&);

  ClassDef(AliRICH,11)                            //Main RICH class 
};//class AliRICH  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
