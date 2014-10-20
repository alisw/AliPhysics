#ifndef ALIVMULTIPLICITY_H
#define ALIVMULTIPLICITY_H

#include <TNamed.h>
#include <TMath.h>


//////////////////////////////////////////////////////////
//                                                      //
//     Base virtual class for multiplicity information  //
//                                                      //
//////////////////////////////////////////////////////////

class AliVMultiplicity : public TNamed {

 public:
  //
  enum {kMultTrackRefs  =BIT(14), // in new format (old is default for bwd.comp.) multiple cluster->track references are allowed
	kScaleDThtbySin2=BIT(15), // scale Dtheta by 1/sin^2(theta). Default is DON'T scale, for bwd.comp.
	kSPD2Sng        =BIT(16)  // are SPD2 singles stored?
  };   
  //
  AliVMultiplicity() {}
 AliVMultiplicity(const char* name, const char* title) : TNamed(name,title) {}
  AliVMultiplicity(const AliVMultiplicity& m) : TNamed(m) {}
  AliVMultiplicity& operator=(const AliVMultiplicity& m) {if (this!=&m) TNamed::operator=(m); return *this;}
  virtual ~AliVMultiplicity() {}
  //
  // methods to access tracklet information
  Bool_t  GetMultTrackRefs()                          const {return TestBit(kMultTrackRefs);}
  Bool_t  GetScaleDThetaBySin2T()                     const {return TestBit(kScaleDThtbySin2);}
  void    SetMultTrackRefs(Bool_t v)                        {SetBit(kMultTrackRefs,v);}
  void    SetScaleDThetaBySin2T(Bool_t v)                   {SetBit(kScaleDThtbySin2,v);}
  //
  virtual void Clear(Option_t* opt="");
  //
  virtual Int_t    GetNumberOfTracklets()             const = 0;
  virtual Double_t GetTheta(Int_t i)                  const = 0;
  virtual Double_t GetPhi(Int_t i)                    const = 0;
  virtual Double_t GetDeltaPhi(Int_t i)               const = 0;
  virtual Int_t    GetLabel(Int_t i, Int_t layer)     const = 0;
  virtual void     SetLabel(Int_t i, Int_t layer, Int_t label) = 0;
  Double_t         GetEta(Int_t i)                    const 
  { 
    if(i>=0 && i<GetNumberOfTracklets()) return -TMath::Log(TMath::Tan(GetTheta(i)/2.));
    Error("GetEta","Invalid track number %d",i); return -9999.;
  }
  //
  // array getters
  virtual Double_t* GetTheta()                        const = 0;
  virtual Double_t* GetPhi()                          const = 0;
  virtual Double_t* GetDeltPhi()                      const = 0;
  virtual Int_t*    GetLabels()                       const = 0;
  virtual Int_t*    GetLabels2()                      const = 0;
  //
  virtual void Print(Option_t *opt="")                const = 0;
  //
  ClassDef(AliVMultiplicity,1);
};


#endif
