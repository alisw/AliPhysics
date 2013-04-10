#ifndef ALIITSUCLUSTERIZER_H
#define ALIITSUCLUSTERIZER_H

#include <TClonesArray.h>
class TTree;
class TObjAray;
class AliITSUSegmentationPix;
class AliITSdigit;
class AliCluster;
class AliITSUClusterPix;
class AliITSURecoParam;

class AliITSUClusterizer : public TObject 
{
  //
 public:
  enum {kDigitChunkSize=1024, kMaxLabels=20,kMaxLabInCluster=3};
  //
  AliITSUClusterizer(Int_t nrowInit=0);
  virtual ~AliITSUClusterizer();
  void SetRawData(Bool_t v=kTRUE)                      {fRawData = v;}
  void Clusterize();
  void SetSegmentation(const AliITSUSegmentationPix *segm);
  void SetRecoParam(const AliITSURecoParam* param)     {fRecoParam = param;}
  void SetLayerID(Int_t id)                            {fLayerID = id;}
  void SetVolID(Int_t id)                              {fVolID = id;}
  void SetNRow(Int_t nrow);
  void SetAllowDiagonalClusterization(Bool_t v)        {fAllowDiagonalClusterization = v;}
  void PrepareLorentzAngleCorrection(Double_t bz);
  //
  // interface methods
  void MakeRecPointBranch(TTree */*treeR*/)            {};
  void SetRecPointTreeAddress(TTree */*treeR*/)        {};
  void DigitsToRecPoints(const TObjArray */*digList*/) {};
  
  void SetDigits(const TClonesArray *digits)       {fInputDigits=digits;fInputDigitsReadIndex=0;}
  void SetClusters(TClonesArray *clusters)         {fOutputClusters=clusters;}
  //
  // labeling methods
  void AddLabel(int label);
  void CheckLabels();
  //
 protected: // transient data types
  struct AliITSUClusterizerClusterDigit {
    AliITSUClusterizerClusterDigit *fNext;
    AliITSdigit *fDigit;
  };
  
  struct AliITSUClusterizerClusterCand;
  struct AliITSUClusterizerClusterPart {
    AliITSUClusterizerClusterPart *fNextInRow;
    AliITSUClusterizerClusterPart *fPrevInCluster;
    AliITSUClusterizerClusterPart *fNextInCluster;
    AliITSUClusterizerClusterCand *fParent;
    Int_t fUBegin;
    Int_t fUEnd;
  };
  
  struct AliITSUClusterizerClusterCand {
    AliITSUClusterizerClusterCand  *fNext; // only used for memory management
    AliITSUClusterizerClusterPart  *fFirstPart;
    AliITSUClusterizerClusterDigit *fFirstDigit;
    AliITSUClusterizerClusterDigit *fLastDigit ;
  };

 protected:
  //
  // allocation and deallocation
  AliITSUClusterizerClusterDigit* AllocDigitFreelist();
  AliITSUClusterizerClusterCand*  AllocCand();
  void                            DeallocCand(AliITSUClusterizerClusterCand *cand);
  AliITSUClusterizerClusterPart*  AllocPart();
  void                            DeallocPart(AliITSUClusterizerClusterPart *part) {DeallocParts(part,part); }
  void                            DeallocParts(AliITSUClusterizerClusterPart *first,AliITSUClusterizerClusterPart *last);
  AliITSUClusterizerClusterDigit* AllocDigit();
  void                            DeallocDigit(AliITSUClusterizerClusterDigit *digit) {DeallocDigits(digit,digit);}
  void                            DeallocDigits(AliITSUClusterizerClusterDigit *first,AliITSUClusterizerClusterDigit *last);

  // input "iterator"
  AliITSUClusterizerClusterDigit* NextDigit();
  // output "iterator"
  AliCluster*                     NextCluster() {return (AliCluster*)fOutputClusters->New(fOutputClusters->GetEntriesFast());}
  
  // modifiers
  void SetAllowDiagonalClusterization();

  void AttachDigitToCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterDigit *digit);
  void AttachPartToCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterPart *part);
  void DetachPartFromCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterPart *part);
  void MergeCands(AliITSUClusterizerClusterCand *a,AliITSUClusterizerClusterCand *b);
  void Transform(AliITSUClusterPix *cluster, AliITSUClusterizerClusterCand *cand);
  void CloseCand(AliITSUClusterizerClusterCand *cand);
  void ClosePart(AliITSUClusterizerClusterPart *part);

  void CloseRemainingParts(AliITSUClusterizerClusterPart *part);
  //
 protected:
  //
  Int_t fVolID;                             // Volume id (module index)
  Bool_t fAllowDiagonalClusterization;      // allow clusters with pixels having common corners only
  const AliITSUSegmentationPix* fSegm;      // Segmentation or local coord calc.
  const AliITSURecoParam*       fRecoParam; // reco params
  //
  // Digit Input
  const TClonesArray *fInputDigits;         // supplied digits
  Int_t         fInputDigitsReadIndex;      // digits counter
  Int_t         fLayerID;                   // current layer id
  //
  Int_t         fCurrLabels[kMaxLabels];    // labels collected for current cluster
  Int_t         fNLabels;                   // number of collected labels
  Bool_t        fRawData;                   // is raw data processed?
  //
  Double_t      fLorAngCorrection;          // Lorentz Angle correction for current layer
  // Cluster Output
  TClonesArray *fOutputClusters;            // external container to store clusters
  //
  // temporary variables
  AliITSUClusterizerClusterDigit *fDigitFreelist    ; //! pool of local digits
  AliITSUClusterizerClusterPart  *fPartFreelist     ; //! pool of unfinished clusters
  AliITSUClusterizerClusterCand  *fCandFreelist     ; //! pool of clusters
  AliITSUClusterizerClusterDigit *fDigitFreelistBptrFirst; //! pointer in the pool
  AliITSUClusterizerClusterDigit *fDigitFreelistBptrLast ; //! pointer in the pool
  AliITSUClusterizerClusterPart  *fPartFreelistBptr ; //! pointer in the pool
  AliITSUClusterizerClusterCand  *fCandFreelistBptr ; //!pointer in the pool
  //
 private:
  AliITSUClusterizer(const AliITSUClusterizer&); //Not implemented
  AliITSUClusterizer& operator=(const AliITSUClusterizer&); //Not implemented
  //
  ClassDef(AliITSUClusterizer,0)
};


//_______________________________________________________________________________
inline void AliITSUClusterizer::DeallocCand(AliITSUClusterizerClusterCand *cand)
{
  // free candidate
  cand->fNext=fCandFreelist;
  fCandFreelist=cand;
}

//_______________________________________________________________________________
inline void AliITSUClusterizer::DeallocParts(AliITSUClusterizerClusterPart *first,AliITSUClusterizerClusterPart *last)
{
  // free cluster part
  last->fNextInRow=fPartFreelist;
  fPartFreelist=first;
}

//_______________________________________________________________________________
inline AliITSUClusterizer::AliITSUClusterizerClusterDigit* AliITSUClusterizer::AllocDigit()
{
  // allocate digits
  if (!fDigitFreelist) fDigitFreelist = AllocDigitFreelist();
  AliITSUClusterizerClusterDigit *digit = fDigitFreelist;
  fDigitFreelist = fDigitFreelist->fNext;
  return digit;
}

//_______________________________________________________________________________
inline AliITSUClusterizer::AliITSUClusterizerClusterPart* AliITSUClusterizer::AllocPart()
{
  // allocate cluster part
  AliITSUClusterizerClusterPart *part=fPartFreelist;
  fPartFreelist=fPartFreelist->fNextInRow;
  return part;
}

//_______________________________________________________________________________
inline AliITSUClusterizer::AliITSUClusterizerClusterCand* AliITSUClusterizer::AllocCand()
{
  // allocate cluster 
  AliITSUClusterizerClusterCand *cand=fCandFreelist;
  fCandFreelist=fCandFreelist->fNext;
  return cand;
}

//_______________________________________________________________________________
inline void AliITSUClusterizer::DeallocDigits(AliITSUClusterizerClusterDigit *first, AliITSUClusterizerClusterDigit *last) 
{
  // free digit
  last->fNext = fDigitFreelist;
  fDigitFreelist = first;
}

//_______________________________________________________________________________
inline void AliITSUClusterizer::AttachDigitToCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterDigit *digit) 
{
  // attach digit
  digit->fNext = cand->fFirstDigit;
  cand->fFirstDigit = digit;
}

//_______________________________________________________________________________
inline void AliITSUClusterizer::DetachPartFromCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterPart *part) 
{
  // remove cluster part 
  if (part->fPrevInCluster)    part->fPrevInCluster->fNextInCluster=part->fNextInCluster;
  else                         cand->fFirstPart=part->fNextInCluster;
  if (part->fNextInCluster)    part->fNextInCluster->fPrevInCluster=part->fPrevInCluster;
}

//______________________________________________________________________________
inline void AliITSUClusterizer::AddLabel(int label)
{
  // add new label
  if (fNLabels==kMaxLabels) return;
  for (int i=fNLabels;i--;) if (fCurrLabels[i]==label) return;
  fCurrLabels[fNLabels++] = label;
}


#endif

