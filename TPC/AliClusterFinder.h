#ifndef TCLUSTERFINDER_H
#define TCLUSTERFINDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// include files and class forward declarations
#include "TObject.h"
#include "TH2.h"
class TClonesArray;
class AliCluster;
class AliCell;
class AliArrayI;
class AliDetectorParam;
class TMinuit;
class AliH2F;
class  AliClusterFinder : public TObject {
 
public:      
  AliClusterFinder(); 
  // constructor which create cluster finder  object  
  ~AliClusterFinder();
  void GetHisto(TH2F * his2);
  //reset object to include histograms values
  TClonesArray * FindPeaks1( TClonesArray *arr=0);  
  TClonesArray * FindPeaks2( TClonesArray *arr=0);  
  TClonesArray * FindPeaks3( TClonesArray *arr=0);  

  void FindMaxima();
  // if at point is local maximum return cell with maximum information  (for testing only
  Int_t & GetNType(){return fNType;}  //return type of neighborow for max determinatio 
  void SetThreshold(Float_t threshold) { fThreshold = threshold;}
  void SetNoise(Float_t noise) {fNoiseTh =noise;}
  void SetDirSigmaFac(Float_t fac) {fDirSigmaFac = fac;}
  void SetDirAmpFac(Float_t fac) {fDirAmpFac = fac;}

  void SetDetectorParam(AliDetectorParam*param) {fDetectorParam = param;}
  //set Detector parameters -necesssary to estimate cluster size
  void SetDetectorIndex(Int_t *index) {fDetectorIndex = index;}
  //set index of described detector
  Bool_t  SetSigma2(Int_t i, Int_t j, Float_t & sigmax2, Float_t & sigmay2);
  //set sigmax and sigma y accordig i and j position of cell 
  void SetMulSigma(Float_t sigma) {fMulSigma2= sigma*sigma;}
  AliArrayI * GetStack(){return fStack;}
  Int_t GetStackIndex(){return fStackIndex;}
  void SetBFit(Bool_t fit) {fBFit = fit;}
  TMinuit * GetMinuit() {return fMinuit;}
  AliH2F *  Draw( const char *option=0,Float_t x1=-1, Float_t x2=-1, Float_t y1=-1, Float_t y2=-1); 
           //draw digits
  void DrawCluster(Int_t color=5, Int_t size=5, Int_t style=4);
  AliH2F *  DrawBorders( const char *option=0,  AliH2F *his=0, Int_t type =0, Float_t x1=-1, Float_t x2=-1, Float_t y1=-1, Float_t y2=-1); 
  //draw digits
public:
  Bool_t   IsMaximum(Int_t i, Int_t  j);
  Bool_t   IsVirtualMaximum(Float_t x, Float_t  y);
  
  void ResetSignal();   //reset signals to 0
  void ResetStatus();   //reset status of signals to not used

  AliCell  * GetCell(Int_t i, Int_t j);   
  //return reference to the cell with index i,j 
  void SetBlockIndex(Int_t *index); //calculate which indexes we must check for border
public:
  void AddToStack(Int_t i, Int_t j, Int_t signal);	  
  //add given cell to the stack of particles
  void GetClusterStatistic(AliDigitCluster & cluster);
  //go through the cluster and calculate statistic
  void GetClusterFit(AliDigitCluster & cluster);
  Bool_t CheckIfDirBorder(Float_t x, Float_t y, Int_t i,Int_t j);
  //check if given cell is border
  void   Adjacent(Int_t i, Int_t j);  
  //recursion procedure 
  Float_t ItoX(Float_t i) {return (fX1)+(i+0.5)*(fX2-fX1)/fDimX;}
  Float_t JtoY(Float_t j) {return (fY1)+(j+0.5)*(fY2-fY1)/fDimY;}

  inline Bool_t IsChecked(Int_t index, Int_t i, Int_t j);
  inline Bool_t IsBorder(Int_t index, Int_t i, Int_t j);
  inline Bool_t IsThBorder(Int_t index, Int_t i, Int_t j);
  inline Bool_t IsDirBorder(Int_t index, Int_t i, Int_t j);
  inline Bool_t IsMaximum(Int_t index, Int_t i, Int_t j);
  inline void  SetChecked(Int_t index, Int_t i, Int_t j);
  inline void  SetBorder(Int_t index, Int_t i, Int_t j);
  inline void  SetThBorder(Int_t index, Int_t i, Int_t j);
  inline void  SetDirBorder(Int_t index, Int_t i, Int_t j);
  inline void  SetMaximum(Int_t index, Int_t i, Int_t j);

  
  
  inline Int_t  GetSignal(Int_t i, Int_t j); 
  Float_t  GetVirtualSignal(Float_t ri, Float_t rj); 
  //create new virtual cell and interpolate signal at  position ri,rj

  void Transform(AliDigitCluster *c);
private:
  Int_t fNType;  //type of neighborow for maximum determination
  Float_t fCurrentMaxX; //!current cluster maximum X index
  Float_t fCurrentMaxY; //!current cluster maximum X index
  Float_t fCurrentSigmaX2; //!current sigmax2 according detector (updated by function ...)
  Float_t fCurrentSigmaY2;//!
  Float_t fCurrentMaxAmp;//!current cluster maximum amplitude
  Bool_t fBDistType;  //
  Bool_t fBFit;  //
  Float_t  fMulSigma2; //
  Float_t fDirSigmaFac;  //!for direction border calculation
  Float_t fDirAmpFac;    //!for direction border calculation

  Float_t fThreshold; //treshold;
  Float_t  fNoiseTh;  //noise threshoshol to accept maximum    
  AliCell   * fDigits;  //field with all cell digits   
 
 
  Int_t      fIndex;        //!index of current  cluster
  AliArrayI  * fStack;      //!stack with digits index
  Int_t      fStackIndex;   //!stack index
  TMinuit    *fMinuit;      //!minuit object
  AliDetectorParam * fDetectorParam; //pointer to detector param  - finder is not owner
  Int_t *    fDetectorIndex; // detector index -  
  //original frame 
  Float_t   fX1;
  Float_t   fY1;
  Float_t   fX2;
  Float_t   fY2;
  Int_t      fDimX;
  Int_t      fDimY;
  Int_t      fOver;     

  //
  TClonesArray * fClustersArray; //array with current clusters
  Bool_t     rOK;       
  //signalize that all fields were initialised 
  ClassDef(AliClusterFinder,2)
};  



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//objec AliCell 

const Int_t krCheck     = 1;
const Int_t krBorder    = 2;
const Int_t krThBorder  = 4;
const Int_t krDirBorder = 8;
const Int_t krMaximum = 16;
const Int_t krIndexNull =0x1F;

class AliCell{
public :
  AliCell(Int_t signal =0, Int_t status = 0){fSignal =signal;fStatus = status;}  
  //at the begining set  
  void SetSignal(Int_t signal){fSignal = signal;}
  void SetStatus(Int_t status){fStatus = status;}  
  void SetChecked(Int_t index) {fStatus &=krIndexNull; fStatus+=(index<<5); fStatus|=krCheck;}
  void SetChecked() {fStatus|=krCheck;}
  void SetBorder(Int_t index) {fStatus &=krIndexNull; fStatus |= krBorder;fStatus+=(index<<5);}
  void SetThBorder(Int_t index) {fStatus &=krIndexNull;fStatus|=krBorder|krThBorder;fStatus+=(index<<5);}
  void SetDirBorder(Int_t index) {fStatus &=krIndexNull;fStatus|=krBorder|krDirBorder;fStatus+=(index<<5);}
  void SetMaximum(Int_t index) {fStatus &=krIndexNull;fStatus|=krMaximum;fStatus+=(index<<5);}


  void SetUnChecked(){if (fStatus&krCheck) fStatus-=krCheck;}
  void SetUnBorder(){if (fStatus&krBorder) fStatus-=krBorder;}
  void SetThUnBorder(){SetUnBorder();if (fStatus&krBorder) fStatus-=krThBorder+krBorder;}
  void SetDirUnBorder(){SetUnBorder();if (fStatus&krBorder) fStatus-=krDirBorder+krBorder;}
  
  Bool_t IsChecked() {return fStatus&&krBorder;}
  Bool_t IsChecked(Int_t index) {return ((fStatus>>5)==index);}
  
  Bool_t  IsBorder() {return ((fStatus&krBorder)!=0);}
  Bool_t  IsBorder(Int_t index) {return ( ((fStatus&krBorder)!=0) && ((fStatus>>5)==index));}  
  Bool_t  IsDirBorder() {return ((fStatus&krDirBorder)!=0);}        
  Bool_t  IsDirBorder(Int_t index) {return ( ((fStatus&krDirBorder)!=0) && ((fStatus>>5)==index));}
  Bool_t  IsThBorder() {return ((fStatus&krThBorder)!=0);}        
  Bool_t  IsThBorder(Int_t index) {return ( ((fStatus&krThBorder)!=0) && ((fStatus>>5)==index));}

  Bool_t  IsMaximum() {return ((fStatus&krMaximum)!=0);}
  Bool_t  IsMaximum(Int_t index) {return ( ((fStatus&krMaximum)!=0) && ((fStatus>>5)==index));}  

  void Reset() { fStatus = 0; fSignal =0;} 
  Int_t GetSignal() {return fSignal;}
  Int_t GetStatus() {return fStatus;}

  
private:
  Int_t fSignal;
  Int_t fStatus;
};



Int_t AliClusterFinder::GetSignal(Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  Int_t res;
  if (c==0) res = -1;
  else  res = c->GetSignal();
  return res;
}



Bool_t AliClusterFinder::IsBorder(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  Bool_t res;
  if (c==0) res = kFALSE;
  else {
    if (index == 0) res = c->IsBorder();
    else res = c->IsBorder(index);
  }
  return res;
}



Bool_t AliClusterFinder::IsThBorder(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  Bool_t res;
  if (c==0) res = kFALSE;
  else {
    if (index == 0) res = c->IsThBorder();
    else res = c->IsThBorder(index);
  }
  return res;
}

Bool_t AliClusterFinder::IsDirBorder(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  Bool_t res;
  if (c==0) res = kFALSE;
  else {
    if (index == 0) res = c->IsDirBorder();
    else res = c->IsDirBorder(index);
  }
  return res;
}

Bool_t AliClusterFinder::IsChecked(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  Bool_t res;
  if (c==0) res = kTRUE;
  else {
    if (index == 0) res = c->IsChecked();
    else res = c->IsChecked(index);
  }
  return res;
}

Bool_t AliClusterFinder::IsMaximum(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  Bool_t res;
  if (c==0) res = kTRUE;
  else {
    if (index == 0) res = c->IsMaximum();
    else res = c->IsMaximum(index);
  }
  return res;
}


void AliClusterFinder::SetChecked(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  if (c!=0) {
    if (index>0) c->SetChecked(index);
    else c->SetChecked();
  }
}
 
void AliClusterFinder::SetBorder(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  if (c!=0) {
    if (index>0) c->SetBorder(index);
    //    else c->SetBorder();
  }
}


void AliClusterFinder::SetThBorder(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  if (c!=0) {
    if (index>0) c->SetThBorder(index);
    else c->SetThBorder(0);
  }
}


void AliClusterFinder::SetDirBorder(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  if (c!=0) {
    if (index>0) c->SetDirBorder(index);
    else c->SetDirBorder(0);
  }
}

 
void AliClusterFinder::SetMaximum(Int_t index, Int_t i, Int_t j)
{
  AliCell *c = GetCell(i,j);
  if (c!=0) {
    if (index>0) c->SetMaximum(index);
    else c->SetMaximum(0);
  }
}

#endif /* TCLUSTERFINDER_H */
