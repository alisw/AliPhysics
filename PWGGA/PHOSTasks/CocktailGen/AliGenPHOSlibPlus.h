#include "AliGenPHOSlib.h"
#include "TF1.h"

class AliGenPHOSlibPlus: public AliGenPHOSlib
{
  
public:
  AliGenPHOSlibPlus() 
  {
    fPart=111;
    fFun=NULL;
  }
  
  AliGenPHOSlibPlus(Int_t part, TF1* ptSpect)
  {
    fPart = part;
    fFun = ptSpect;
  }
  
  GenFunc   GetPt(Int_t , const char* ) const {return PtPion ;}
  GenFunc   GetY (Int_t , const char* ) const {return YPion ;}
  GenFuncIp GetIp(Int_t , const char* ) const {return IpPion ;}
  
  static Double_t PtPion(const Double_t *px, const Double_t *)
  {
    // printf("PtPion: return %f of %s\n",fFun->Eval((*px)),fFun->GetName());
    return fFun->Eval((*px));
  }
  
  static Int_t    IpPion(TRandom* ){return fPart ;}
  static Double_t YPion( const Double_t *, const Double_t *){return 1.;}
  
public:
  
  static Int_t fPart ;
  static TF1* fFun; // pT-spectrum
  
  ClassDef(AliGenPHOSlibPlus,1)
  
} ;

Int_t AliGenPHOSlibPlus::fPart = 111 ;
TF1*  AliGenPHOSlibPlus::fFun = new TF1();
