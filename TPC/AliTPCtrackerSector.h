#ifndef ALITPCTRACKERSECTOR_H
#define ALITPCTRACKERSECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: AliTPCtrackerSector.h 25837 2008-05-16 16:39:00Z marian $ */

//-------------------------------------------------------
//   TPC tracker - helper classes for cluster storing
//                 and navigation                   
//   
//
//   Origin: 
//-------------------------------------------------------


//class TFile;
class AliTPCParam;
class TTreeSRedirector;



class AliTPCtrackerRow : public TObject{
public:
  AliTPCtrackerRow();
  ~AliTPCtrackerRow();
  void InsertCluster(const AliTPCclusterMI *c, UInt_t index);
  void ResetClusters();
  operator int() const {return fN;}
  Int_t GetN() const {return fN;}
  const AliTPCclusterMI* operator[](Int_t i) const {return fClusters[i];}
  UInt_t GetIndex(Int_t i) const {return fIndex[i];}
  Int_t Find(Double_t z) const; 
  AliTPCclusterMI *  FindNearest(Double_t y, Double_t z, Double_t roady, Double_t roadz) const;
  AliTPCclusterMI *  FindNearest2(Double_t y, Double_t z, Double_t roady, Double_t roadz, UInt_t & index) const;
  
  void SetX(Double_t x) {fX=x;}
  Double_t GetX() const {return fX;}
  Float_t GetDeadZone() const {return fDeadZone;}
  void SetDeadZone(Float_t d) {fDeadZone=d;}
  Int_t GetN1() const {return fN1;}
  void SetN1(Int_t n) {fN1=n;}
  Int_t GetN2() const {return fN2;}
  void SetN2(Int_t n) {fN2=n;}
  TClonesArray* GetClusters1() const {return fClusters1;}
  TClonesArray* GetClusters2() const {return fClusters2;}
  void SetCluster1(Int_t i, const AliTPCclusterMI &cl);
  void SetCluster2(Int_t i, const AliTPCclusterMI &cl);
  AliTPCclusterMI* GetCluster1(Int_t i) const {return (AliTPCclusterMI*) fClusters1->At(i);}
  AliTPCclusterMI* GetCluster2(Int_t i) const {return (AliTPCclusterMI*) fClusters2->At(i);}
  Short_t GetFastCluster(Int_t i) const {return fFastCluster[i];}
  void SetFastCluster(Int_t i, Short_t cl);
  Int_t IncrementN1() { return ++fN1;}
  Int_t IncrementN2() { return ++fN2;}
  
private:  
  AliTPCtrackerRow & operator=(const AliTPCtrackerRow & );
  AliTPCtrackerRow(const AliTPCtrackerRow& /*r*/);           //dummy copy constructor
  Float_t fDeadZone;  // the width of the dead zone
  TClonesArray *fClusters1; //array with clusters 1
  Int_t fN1;  //number of clusters on left side
  TClonesArray *fClusters2; //array with clusters 2
  Int_t fN2; // number of clusters on right side of the TPC
  Short_t fFastCluster[510];   //index of the nearest cluster at given position
  Int_t fN;                                          //number of clusters 
  const AliTPCclusterMI *fClusters[kMaxClusterPerRow]; //pointers to clusters
  // indexes for cluster at given position z  
  // AliTPCclusterMI *fClustersArray;                     // 
  UInt_t fIndex[kMaxClusterPerRow];                  //indeces of clusters
  Double_t fX;                                 //X-coordinate of this row  
  ClassDef(AliTPCtrackerRow,0)
};


//**************** Internal tracker class ********************** 
class AliTPCtrackerSector: public TObject {
 public:
  AliTPCtrackerSector():
    fN(0),
    fRow(0),
    fAlpha(0.),
    fAlphaShift(0.),
    fPadPitchWidth(0.),
    fPadPitchLength(0.),
    f1PadPitchLength(0.),
    f2PadPitchLength(0.) {}
    ~AliTPCtrackerSector() { delete[] fRow; }
    AliTPCtrackerRow& operator[](Int_t i) const { return *(fRow+i); }
    Int_t GetNRows() const { return fN; }
    void Setup(const AliTPCParam *par, Int_t flag);
    Double_t GetX(Int_t l) const {return fRow[l].GetX();}
    Double_t GetMaxY(Int_t l) const {
      return GetX(l)*TMath::Tan(0.5*GetAlpha());
    } 
    Double_t GetAlpha() const {return fAlpha;}
    Double_t GetAlphaShift() const {return fAlphaShift;}     
    //Int_t GetFirst(){return fFirstRow;}
    Int_t GetRowNumber(Double_t  x) const;
    Double_t GetPadPitchWidth()  const {return fPadPitchWidth;}
    Double_t GetPadPitchLength() const {return fPadPitchLength;}
    Double_t GetPadPitchLength(Float_t x) const {return (x<200) ? fPadPitchLength:f2PadPitchLength ;}
    
    void InsertCluster(AliTPCclusterMI *cl, Int_t size, const AliTPCParam *par);

 private:
    AliTPCtrackerSector & operator=(const AliTPCtrackerSector & );
    AliTPCtrackerSector(const AliTPCtrackerSector &/*s*/);           //dummy copy contructor 
    Int_t fN;                           //number of pad rows 
    //Int_t fFirstRow;                  //offset
    AliTPCtrackerRow *fRow;             //array of pad rows
    Double_t fAlpha;                    //opening angle
    Double_t fAlphaShift;               //shift angle;
    Double_t fPadPitchWidth;            //pad pitch width
    Double_t fPadPitchLength;           //pad pitch length
    Double_t f1PadPitchLength;          //pad pitch length
    Double_t f2PadPitchLength;          //pad pitch length  
    ClassDef(AliTPCtrackerSector,1)
};



#endif


