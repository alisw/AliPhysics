#ifndef AliHeader_H
#define AliHeader_H

#include "TObject.h"
 
class AliHeader : public TObject {
protected:
  Int_t         fRun;         //Run number
  Int_t         fNvertex;     //Number of vertices
  Int_t         fNprimary;    //Number of primary tracks
  Int_t         fNtrack;      //Number of tracks
  Int_t         fEvent;       //Event number

public:
  AliHeader();
  AliHeader(Int_t run, Int_t event);
  ~AliHeader() {;}

  virtual void Reset(Int_t run, Int_t event);

  virtual  void  SetRun(Int_t run) {fRun = run;}
  virtual  Int_t GetRun() const {return fRun;}
  
  virtual  void  SetNprimary(Int_t nprimary) {fNprimary = nprimary;}
  virtual  Int_t GetNprimary() const {return fNprimary;}
  
  virtual  void  SetNvertex(Int_t vertex) {fNvertex = vertex;}
  virtual  Int_t GetNvertex() const {return fNvertex;}
  
  virtual  void  SetNtrack(Int_t ntrack) {fNtrack = ntrack;}
  virtual  Int_t GetNtrack() const {return fNtrack;}
  
  virtual  void  SetEvent(Int_t event) {fEvent = event;}
  virtual  Int_t GetEvent() const {return fEvent;}

  virtual void Dump();
  
  ClassDef(AliHeader,1) //Alice event header
    
};

#endif
