#ifndef AliModule_H
#define AliModule_H
#include <TNamed.h>
#include <TClonesArray.h>
#include <TBrowser.h>
#include <TAttLine.h>
#include <TAttMarker.h>
#include <TArrayI.h>
#include <AliHit.h>

class AliModule : public TNamed , public TAttLine, public TAttMarker {

  // Data members
protected:      
  
  TString       fEuclidMaterial;  //Name of the Euclid file for materials (if any)
  TString       fEuclidGeometry;  //Name of the Euclid file for geometry (if any)
  
  TArrayI      *fIdtmed;      //List of tracking medium numbers
  TArrayI      *fIdmate;      //List of material numbers
  Int_t         fLoMedium;   //Minimum tracking medium ID for this Module
  Int_t         fHiMedium;   //Maximum tracking medium ID for this Module

  Bool_t        fActive;      //Detector activity flag
  TList        *fHistograms;  //List of histograms
  TList        *fNodes;       //List of geometry nodes

public:

  // Creators - distructors
  AliModule(const char* name, const char *title);
  AliModule();
  virtual ~AliModule();

  // Inline functions
  virtual  int           GetNdigits() {return 0;}
  virtual  int           GetNhits()   {return 0;}
  virtual  TArrayI      *GetIdtmed()   {return fIdtmed;}
  virtual  TList        *Histograms() {return fHistograms;}
  virtual  TList        *Nodes()  {return fNodes;}
  virtual  TClonesArray *Digits() {return 0;}
  virtual  TClonesArray *Hits()   {return 0;}
  virtual  TObjArray    *Points() {return 0;}
  virtual  Int_t         GetIshunt() {return 0;}
  virtual  void          SetIshunt(Int_t) {}
  virtual  Bool_t        IsActive() {return fActive;}
  virtual  Bool_t        IsFolder() {return kTRUE;}
  virtual  Int_t&        LoMedium() {return fLoMedium;}
  virtual  Int_t&        HiMedium() {return fHiMedium;}

  // Module composition
  virtual void  AliMaterial(Int_t, const char*, Float_t, Float_t, Float_t, Float_t,
			    Float_t, Float_t* buf=0, Int_t nwbuf=0) const;
  virtual void  AliMixture(Int_t, const char*, Float_t*, Float_t*, Float_t, Int_t, Float_t*) const;
  virtual void  AliMedium(Int_t, const char*, Int_t, Int_t, Int_t, Float_t, Float_t, 
		   Float_t, Float_t, Float_t, Float_t, Float_t* ubuf=0, Int_t nbuf=0) const;
  virtual void  AliMatrix(Int_t&, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) const;
  
  // Virtual methods
  virtual void  BuildGeometry() {};
  virtual Int_t IsVersion() const =0;

  // Other methods
  virtual void        AddDigit(Int_t*, Int_t*){}
  virtual void        AddHit(Int_t, Int_t*, Float_t *) {}
  virtual void        Browse(TBrowser *) {}
  virtual void        CreateGeometry() {}
  virtual void        CreateMaterials() {}
  virtual void        Disable();
  virtual Int_t       DistancetoPrimitive(Int_t px, Int_t py);
  virtual void        Enable();
  virtual void        PreTrack(){}
  virtual void        PostTrack(){}
  virtual void        FinishEvent() {}
  virtual void        FinishRun() {}
  virtual void        Hits2Digits() {}
  virtual void        Init() {}
  virtual void        LoadPoints(Int_t ) {}
  virtual void        MakeBranch(Option_t *) {}
  virtual void        Paint(Option_t *) {}
  virtual void        ResetDigits() {}
  virtual void        ResetHits() {}
  virtual void        ResetPoints() {}
  virtual void        SetTreeAddress() {}
  virtual void        SetTimeGate(Float_t) {}
  virtual Float_t     GetTimeGate() {return 1.e10;}
  virtual void        StepManager() {}
  virtual AliHit*     FirstHit(Int_t) {return 0;}
  virtual AliHit*     NextHit() {return 0;}
  virtual void        SetBufferSize(Int_t) {}  
  virtual void        SetEuclidFile(char*,char*geometry=0);
 
  ClassDef(AliModule,1)  //Base class for ALICE Modules
};
#endif
