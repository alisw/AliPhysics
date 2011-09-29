#ifndef ALICONVERSIONPHOTONBASE_H
#define ALICONVERSIONPHOTONBASE_H

#include "TMath.h"
#include "TParticle.h"
#include "AliStack.h"
#include "AliLog.h"
#include "TObject.h"

using namespace std;

class AliConversionPhotonBase {

 public: 

  //Constructors
  AliConversionPhotonBase();    

  //Copy Constructor
  AliConversionPhotonBase(const AliConversionPhotonBase & g);           
  //assignment operator
  AliConversionPhotonBase & operator = (const AliConversionPhotonBase & g);

  //Destructor
  virtual ~AliConversionPhotonBase();

  ///Set the tag for decay meson
  void SetTag( Bool_t tagged ) { fTagged = tagged; }
  Bool_t IsTagged(){return fTagged;}

  //Get the Chi2 of particle
  void SetChi2perNDF(Float_t chi2) {fChi2perNDF = chi2;}
  Float_t GetChi2perNDF() const {return fChi2perNDF;}

  ///Track labels
  void SetLabelPositive(Int_t label){fLabel[0] = label;}
  void SetLabelNegative(Int_t label){fLabel[1] = label;}
  void SetTrackLabels(Int_t label1, Int_t label2){fLabel[0] = label1; fLabel[1] = label2;}
  Int_t GetTrackLabelPositive() const{return fLabel[0];}
  Int_t GetTrackLabelNegative() const {return fLabel[1];}
  Int_t GetTrackLabel(Int_t i) const {return fLabel[i];}

  // MC Label

  void SetMCLabel(Int_t* label){fMCLabel[0]=label[0];fMCLabel[1]=label[1];}
  void SetMCLabelPositive(Int_t label){fMCLabel[0]=label;}
  void SetMCLabelNegative(Int_t label){fMCLabel[1]=label;}
  Int_t GetMCLabel(Int_t i) const{return fMCLabel[i];}
  Int_t GetMCLabelPositive() const{return fMCLabel[0];}
  Int_t GetMCLabelNegative() const{return fMCLabel[1];}
  Int_t GetMCParticleLabel(AliStack *fMCStack);

  // GetMCParticle

  TParticle *GetMCParticle(AliStack *fMCStack);
  TParticle *GetPositiveMCDaughter(AliStack *fMCStack){return GetMCDaughter(fMCStack,0);};
  TParticle *GetNegativeMCDaughter(AliStack *fMCStack){return GetMCDaughter(fMCStack,1);};
  TParticle *GetMCDaughter(AliStack *fMCStack,Int_t label);

  // V0Index
  Int_t GetV0Index(){return fV0Index;}
  void SetV0Index(Int_t index){fV0Index=index;}

  // Conversion Point
   void SetConversionPoint(Double_t convpoint[3]){fConversionPoint[0]=convpoint[0];fConversionPoint[1]=convpoint[1];fConversionPoint[2]=convpoint[2];}
  void GetConversionPoint(Double_t convpoint[3]){convpoint[0]=fConversionPoint[0];convpoint[1]=fConversionPoint[1];convpoint[2]=fConversionPoint[2];}
  Double_t GetConversionRadius(){return TMath::Sqrt(fConversionPoint[0]*fConversionPoint[0]+fConversionPoint[1]*fConversionPoint[1]);}
  Double_t GetConversionX(){return fConversionPoint[0];}
  Double_t GetConversionY(){return fConversionPoint[1];}
  Double_t GetConversionZ(){return fConversionPoint[2];}

  // Armenteros Qt Alpha
  void GetArmenterosQtAlpha(Double_t qtalpha[2]){qtalpha[0]=fArmenteros[0];qtalpha[1]=fArmenteros[1];}
  Double_t GetArmenterosQt(){return fArmenteros[0];}
  Double_t GetArmenterosAlpha(){return fArmenteros[1];}

  // dEdx
  void SetNSigmadEdx(Float_t positive[5],Float_t negative[5]){for(Int_t jj=0;jj<5;jj++){fNSigmadEdxPositive[jj]=positive[jj];fNSigmadEdxNegative[jj]=negative[jj];}}
  void GetNSigmadEdx(Float_t positive[5],Float_t negative[5]){for(Int_t jj=0;jj<5;jj++){positive[jj]=fNSigmadEdxPositive[jj];negative[jj]=fNSigmadEdxNegative[jj];}}
  Float_t GetNSigmadEdxPositive(Int_t species){return fNSigmadEdxPositive[species];}
  Float_t GetNSigmadEdxNegative(Int_t species){return fNSigmadEdxNegative[species];}
  Float_t GetNSigmadEdx(Int_t label,Int_t species){if(label==0){return fNSigmadEdxPositive[species];}if(label==1){return fNSigmadEdxNegative[species];}printf("label not defined");;return 1000;}

  // virtual functions to be implemented in KF/AOD classes

  virtual Double_t GetPhotonMass() const = 0;
  virtual Double_t GetPhotonPt()const = 0;
  virtual Double_t GetPhotonP() const = 0;
  virtual Double_t GetPhotonEta() const = 0;


 protected:

  Int_t fLabel[2]; // Electron/Positron Track Label
  Int_t fV0Index; // Index of the V0
  Int_t fMCLabel[2]; // Electron/Positron MC Label
  Float_t fChi2perNDF; // Chi2perNDF
  Double_t fArmenteros[2]; // Armenteros Paramters
  Float_t fNSigmadEdxPositive[5]; // N Sigma to lines
  Float_t fNSigmadEdxNegative[5]; // N Sigma to lines
  Double_t fConversionPoint[3]; // Conversion Point
  Bool_t fTagged; // Is it tagged as decay pion (only for gammas)

  ClassDef(AliConversionPhotonBase,1)
};


#endif



