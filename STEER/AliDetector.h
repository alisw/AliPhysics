#ifndef AliDetector_H
#define AliDetector_H
#include <TNamed.h>
#include <TClonesArray.h>
#include <TBrowser.h>
#include <TAttLine.h>
#include <TAttMarker.h>
#include <TArrayI.h>
#include <Gtypes.h>
#include <AliHit.h>

class AliDetector : public TNamed , public TAttLine, public TAttMarker {

  // Data members
protected:      
  
  TString       fEuclidMaterial;  //Name of the Euclid file for materials (if any)
  TString       fEuclidGeometry;  //Name of the Euclid file for geometry (if any)
  
  TArrayI      *fIdtmed;      //List of tracking medium numbers
  TArrayI      *fIdmate;      //List of material numbers
  Int_t         fLoMedium;   //Minimum tracking medium ID for this detector
  Int_t         fHiMedium;   //Maximum tracking medium ID for this detector

  Float_t       fTimeGate;    //Time gate in seconds

  Bool_t        fActive;      //Detector activity flag
  Int_t         fIshunt;      //1 if the hit is attached to the primary
  Int_t         fNhits;       //Number of hits
  Int_t         fNdigits;     //Number of digits
  Int_t         fBufferSize;  //buffer size for Tree detector branches
  TList        *fHistograms;  //List of histograms
  TList        *fNodes;       //List of geometry nodes
  TClonesArray *fHits;        //List of hits for one track only
  TClonesArray *fDigits;      //List of digits for this detector
  TObjArray    *fPoints;      //Array of points for each track (all tracks in memory)

public:

  // Creators - distructors
  AliDetector(const char* name, const char *title);
  AliDetector();
  virtual ~AliDetector();

  // Inline functions
  inline virtual int   GetNdigits() {return fNdigits;}
  inline virtual int   GetNhits()   {return fNhits;}
  inline TList        *Histograms() {return fHistograms;}
  inline TList        *Nodes()  {return fNodes;}
  inline TClonesArray *Digits() {return fDigits;}
  inline TClonesArray *Hits()   {return fHits;}
  inline TObjArray    *Points() {return fPoints;}
  inline Int_t         GetIshunt() {return fIshunt;}
  inline void          SetIshunt(Int_t ishunt) {fIshunt=ishunt;}
  inline Bool_t        IsActive() {return fActive;}
  inline Bool_t        IsFolder() {return kTRUE;}
  inline Int_t&        LoMedium() {return fLoMedium;}
  inline Int_t&        HiMedium() {return fHiMedium;}

  // Detector composition
  virtual void  AliMaterial(Int_t, const char*, Float_t, Float_t, Float_t, Float_t,
			    Float_t, Float_t* buf=0, Int_t nwbuf=0) const;
  virtual void  AliMixture(Int_t, const char*, Float_t*, Float_t*, Float_t, Int_t, Float_t*) const;
  virtual void  AliMedium(Int_t, const char*, Int_t, Int_t, Int_t, Float_t, Float_t, 
		   Float_t, Float_t, Float_t, Float_t, Float_t* ubuf=0, Int_t nbuf=0) const;
  virtual void  AliMatrix(Int_t&, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) const;
  
  // Virtual methods
  virtual void  BuildGeometry()=0;
  virtual Int_t IsVersion() const =0;

  // Other methods
  virtual void        AddDigit(Int_t*, Int_t*){}
  virtual void        AddHit(Int_t, Int_t*, Float_t *) {}
  virtual void        Browse(TBrowser *b);
  virtual void        CreateGeometry() {}
  virtual void        CreateMaterials() {}
  virtual void        Disable();
  virtual Int_t       DistancetoPrimitive(Int_t px, Int_t py);
  virtual void        Enable();
  virtual inline void PreTrack(){}
  virtual inline void PostTrack(){}
  virtual inline void FinishEvent() {}
  virtual void        FinishRun();
  virtual void        Hits2Digits() {}
  virtual void        Init() {}
  virtual void        LoadPoints(Int_t track);
  virtual void        MakeBranch(Option_t *opt=" ");
  virtual void        Paint(Option_t *) {}
  virtual void        ResetDigits();
  virtual void        ResetHits();
  virtual void        ResetPoints();
  virtual void        SetTreeAddress();
  virtual void        SetTimeGate(Float_t gate) {fTimeGate=gate;}
  virtual Float_t     GetTimeGate() {return fTimeGate;}
  virtual void        StepManager();
  virtual void        DrawDetector() {}
  virtual AliHit*     FirstHit(Int_t);
  virtual AliHit*     NextHit();
  virtual void        SetBufferSize(Int_t bufsize=8000) {fBufferSize = bufsize;}  
  virtual void        SetEuclidFile(char*,char*geometry=0);
 
  ClassDef(AliDetector,1)  //Base class for ALICE detectors
};
#endif
