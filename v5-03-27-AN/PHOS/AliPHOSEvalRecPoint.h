#ifndef AliPHOSEvalRecPoint_H
#define AliPHOSEvalRecPoint_H

// --- ROOT system ---
#include "TObjArray.h"


#include "AliPHOSCpvRecPoint.h"
class AliPHOSClusterizer;
class AliPHOSRecPoint;
class AliPHOSRecManager;
class TString;

class AliPHOSEvalRecPoint : public AliPHOSCpvRecPoint
{
public:

  AliPHOSEvalRecPoint();
  AliPHOSEvalRecPoint(Bool_t cpv,AliPHOSEvalRecPoint* parent);
  AliPHOSEvalRecPoint(Int_t cluster, Bool_t cpv);
  virtual ~AliPHOSEvalRecPoint() {}

  Bool_t TooClose(AliPHOSRecPoint* pt) const ;
  Bool_t NeedToSplit() const ;

  void MergeClosePoint();
  void SplitMergedShowers();
  Int_t UnfoldLocalMaxima();
  void EvaluatePosition();
  Bool_t KillWeakPoint();

  Bool_t IsEmc(void) const { return fIsEmc; }
  Bool_t IsCPV(void) const { return fIsCpv; }

  void SetLocalPosition(TVector3& pos) { fLocPos=pos; }
  void UpdateWorkingPool();
  void DeleteParent();

  Int_t HasChild(TObjArray& children);

  void MakeJob();

  AliPHOSClusterizer* GetClusterizer();
  AliPHOSRecManager* GetReconstructionManager() const;
  
  void PrintPoint();

  AliPHOSRecPoint* Parent();
  Float_t Chi2Dof() const;

  const TObject* GetWorkingPool();

  void AddToWorkingPool(TObject* obj);
  TObject* GetFromWorkingPool(Int_t i);
  Int_t InWorkingPool();
  void RemoveFromWorkingPool(TObject* obj);
  void PrintWorkingPool();

  enum RecPointType {kEmc,kCpv};

private:

  void Init();
  void InitTwoGam(Float_t* gamma1, Float_t* gamma2);
  void TwoGam(Float_t* gamma1, Float_t* gamma2);
  void UnfoldTwoMergedPoints(Float_t* gamma1, Float_t* gamma2);

private:

  Bool_t fIsEmc;    // true if EMC
  Bool_t fIsCpv;    // true if CPV
  Int_t fParent;    // parent
  Float_t fChi2Dof; // chi^2 per degree of freedom
  TString fEventFolderName; // Folder name

  ClassDef(AliPHOSEvalRecPoint,1)

};

#endif
