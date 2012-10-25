#ifndef ALIITSUCLUSTERIZER_H
#define ALIITSUCLUSTERIZER_H

#include <AliCluster.h>
#include <AliITSdigit.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <AliITSUSegmentationPix.h>

class AliITSUClusterizer : public TObject 
{
  //
 public:
  enum {kDigitChunkSize=1024};
  //
  AliITSUClusterizer(Int_t nrowInit=0);
  virtual ~AliITSUClusterizer();
  void Clusterize();
  void SetSegmentation(const AliITSUSegmentationPix *segm);
  void SetVolID(Int_t id)                              {fVolID = id;}
  void SetNRow(Int_t nrow);
  // interface methods
  void MakeRecPointBranch(TTree */*treeR*/)            {};
  void SetRecPointTreeAddress(TTree */*treeR*/)        {};
  void DigitsToRecPoints(const TObjArray */*digList*/) {};
  
  void SetDigits(TClonesArray *digits)             {fInputDigits=digits;fInputDigitsReadIndex=0;}
  void SetClusters(TClonesArray *clusters)         {fOutputClusters=clusters;}

 protected: // transient data types
  struct AliITSUClusterizerClusterDigit {
    AliITSUClusterizerClusterDigit *fNext;
    AliITSdigit *fDigit;
    Int_t fU;
    Int_t fV;
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
  AliCluster*                     NextCluster() {return new( (*fOutputClusters)[fOutputClusters->GetEntries()] ) AliCluster();}
  
  // modifiers
  void AttachDigitToCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterDigit *digit);
  void AttachPartToCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterPart *part);
  void DetachPartFromCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterPart *part);
  void MergeCands(AliITSUClusterizerClusterCand *a,AliITSUClusterizerClusterCand *b);

  void Transform(AliCluster *cluster,AliITSUClusterizerClusterCand *cand);
  void CloseCand(AliITSUClusterizerClusterCand *cand);
  void ClosePart(AliITSUClusterizerClusterPart *part);

  void CloseRemainingParts(AliITSUClusterizerClusterPart *part);
  //
 protected:
  //
  Int_t fVolID;                             // Volume id (module index)
  const AliITSUSegmentationPix* fSegm;      // Segmentation or local coord calc.
  //
  // Digit Input
  TClonesArray *fInputDigits;
  Int_t         fInputDigitsReadIndex;
  // Cluster Output
  TClonesArray *fOutputClusters;
    
  // temporary variables
  AliITSUClusterizerClusterDigit *fDigitFreelist    ; //! 
  AliITSUClusterizerClusterPart  *fPartFreelist     ; //!
  AliITSUClusterizerClusterCand  *fCandFreelist     ; //!
  AliITSUClusterizerClusterDigit *fDigitFreelistBptrFirst; //!
  AliITSUClusterizerClusterDigit *fDigitFreelistBptrLast ; //!
  AliITSUClusterizerClusterPart  *fPartFreelistBptr ; //!
  AliITSUClusterizerClusterCand  *fCandFreelistBptr ; //!
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


#endif

