#ifndef ALITPCCLUSTERERMI_H
#define ALITPCCLUSTERERMI_H

//-------------------------------------------------------
//                       TPC clusterer
//
//   Origin: Marian Ivanov  
//-------------------------------------------------------
#include <Rtypes.h>
#include <TObject.h>
#define kMAXCLUSTER 2500

class TFile;
class AliTPCParam;
class AliTPCclusterMI;
class AliTPCClustersRow;
class AliSimDigits;
class TTree;

class AliTPCclustererMI : public TObject{
public:
  AliTPCclustererMI();
  virtual void Digits2Clusters(const AliTPCParam *par, Int_t eventn=1);
  virtual void SetInput(TTree * tree);  // set input tree with digits    
  virtual void SetOutput(TTree * tree); //set output tree with 
private:
  Bool_t IsMaximum(Int_t k, Int_t max, const Int_t *bins); 
  void MakeCluster2(Int_t k,Int_t max,Int_t *bins,UInt_t m,
   AliTPCclusterMI &c);  
  void MakeCluster(Int_t k,Int_t max,Int_t *bins,UInt_t m,
   AliTPCclusterMI &c); 
  Float_t  GetSigmaY2(Int_t iz);
  Float_t  GetSigmaZ2(Int_t iz);
  Float_t  FitMax(Float_t vmatrix[5][5], Float_t y, Float_t z, Float_t sigmay, Float_t sigmaz);
  void AddCluster(AliTPCclusterMI &c);  // add the cluster to the array
  void UnfoldCluster(Int_t * matrix[7], Float_t recmatrix[5][5], 
		     Float_t & meani, Float_t & meanj, Float_t & sum, Float_t &overlap );



  Int_t * fBins;       //!digits array
  Int_t * fResBins;    //!digits array with res. after 1 finder
  Int_t fLoop;         //loop - cf in 2 loops
  Int_t fMaxBin;
  Int_t fMaxTime;
  Int_t fMaxPad;
  Int_t fSector;      //!current sector
  Float_t fSign;      //!current sign 
  Float_t fRx;        // current radius
  Float_t fPadWidth;  // the width of the pad
  Float_t fPadLength;  // the width of the pad
  Float_t fZWidth;     //the z bin width

  TTree * fInput;   //!input  tree with digits - object not owner
  TTree * fOutput;   //!output tree with digits - object not owner
  AliTPCClustersRow * fRowCl;  //! current cluster row
  AliSimDigits * fRowDig;      //! current digits row
  const AliTPCParam * fParam;        //! tpc parameters
  Int_t fNcluster;             // number of clusters - for given row
  ClassDef(AliTPCclustererMI,1)  // Time Projection Chamber digits
};

inline Bool_t AliTPCclustererMI::IsMaximum(Int_t q,Int_t max,const Int_t *bins){
  //is this a local maximum ?
  if (bins[-max] >= q) return kFALSE;
  if (bins[-1  ] >= q) return kFALSE; 
  if (bins[+max] > q) return kFALSE; 
  if (bins[+1  ] > q) return kFALSE; 
  if (bins[-max-1] >= q) return kFALSE;
  if (bins[+max-1] >= q) return kFALSE; 
  if (bins[+max+1] > q) return kFALSE; 
  if (bins[-max+1] >= q) return kFALSE;
  return kTRUE; 
}



//-----------------------------------------------------------------

#endif


