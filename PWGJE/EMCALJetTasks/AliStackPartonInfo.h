#ifndef ALISTACKPARTONINFO_H
#define ALISTACKPARTONINFO_H

#include <TMath.h>
#include <TNamed.h>

class AliStackPartonInfo : public TNamed{

 public:
  AliStackPartonInfo();
  AliStackPartonInfo(const char* name); 

  void SetPartonFlag6(Int_t flag6) {fPartonFlag6 = flag6;}
  void SetPartonPt6(Float_t pt6) {fPartonPt6 = pt6;}
  void SetPartonEta6(Float_t eta6) {fPartonEta6 = eta6;}
  void SetPartonPhi6(Float_t phi6) {fPartonPhi6 = phi6;}

  void SetPartonFlag7(Int_t flag7) {fPartonFlag7 = flag7;}
  void SetPartonPt7(Float_t pt7) {fPartonPt7 = pt7;}
  void SetPartonEta7(Float_t eta7) {fPartonEta7 = eta7;}
  void SetPartonPhi7(Float_t phi7) {fPartonPhi7 = phi7;}

  
  Int_t GetPartonFlag6() {return fPartonFlag6;}
  Float_t GetPartonPt6() {return fPartonPt6;}
  Float_t GetPartonEta6() {return fPartonEta6;}
  Float_t GetPartonPhi6() {return fPartonPhi6;}

  Int_t GetPartonFlag7() {return fPartonFlag7;}
  Float_t GetPartonPt7() {return fPartonPt7;}
  Float_t GetPartonEta7() {return fPartonEta7;}
  Float_t GetPartonPhi7() {return fPartonPhi7;}

 private: 
  Int_t fPartonFlag6; //! parton flag 
  Float_t fPartonPt6; //! pT parton 
  Float_t fPartonEta6; //!eta parton 
  Float_t fPartonPhi6; //! phi parton

  Int_t fPartonFlag7; //! parton flag 
  Float_t fPartonPt7; //! pT parton 
  Float_t fPartonEta7; //!eta parton 
  Float_t fPartonPhi7; //! phi parton
    
  AliStackPartonInfo(const AliStackPartonInfo&);
  AliStackPartonInfo& operator=(const AliStackPartonInfo&);
  
  ClassDef(AliStackPartonInfo, 1);

};
#endif
