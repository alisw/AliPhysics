// -*- C++ -*-

#ifndef PWGPP_VDM_ALIVDMTREE_H
#define PWGPP_VDM_ALIVDMTREE_H

#include <map>
#include <cmath>
#include <cfloat>

#include <TObject.h>
#include <TTree.h>
#include <TBranch.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TString.h>

// -----------------------------------------------------------------------------
// this class holds data for a single scan and a single trigger combination
class AliVdMTree : public TObject {
public:
  class DefaultBranchData {
  public:
    DefaultBranchData()
      : fTime(2)
      , fSep(2)
      , fBC(0)
      , fCounter(30) {}

    DefaultBranchData& operator=(const DefaultBranchData& d) {
      for (Int_t i=0; i<2; ++i) {
        fTime[i] = d.fTime[i];
        fSep[i]  = d.fSep[i];
      }
      fBC      = d.fBC;
      fCounter = d.fCounter;
      return *this;
    }
    Double_t& StartTime() { return fTime(0); }
    Double_t& EndTime() { return fTime(1); }
    Short_t& BCID() { return fBC; }
    Double_t& Counter(Int_t idx=0) { return fCounter(idx); }

    Double_t StartTime() const { return fTime(0); }
    Double_t EndTime() const { return fTime(1); }
    Short_t BCID() const { return fBC; }
    Double_t Counter(Int_t idx=0) const { return fCounter(idx); }
    const TVectorD& Counters() const { return fCounter; }

    Double_t DeltaT() const { return EndTime() - StartTime(); }

    Double_t Sep(Int_t idx) const { return fSep(idx); }

    void Branch(TTree* t) {
      t->Branch("fTime",     fTime.GetMatrixArray(),    "start/D:stop");
      t->Branch("fSep",      fSep.GetMatrixArray(),     "X/D:Y");
      t->Branch("fBC",      &fBC);
      t->Branch("fCounter",  fCounter.GetMatrixArray(), Form("counter[%d]/D", fCounter.GetNoElements()));
    }
  protected:
  private:
    TVectorD    fTime;      // time interval [begin,end]
    TVectorD    fSep;       // (X,Y) [mm]
    Short_t     fBC;        // BCID
    TVectorD    fCounter;   // number of triggers for current and 9 past BCs
  } ;

  class ValErr {
  public:
    ValErr(Double_t val=0, Double_t err=0)
      : fArray{val,err} {}

    // ValErr(const ValErr& v)
    //   : ValErr(v.val(), v.err()) {}

    // ValErr& operator=(const ValErr& v) {
    //   fArray[0] = v.fArray[0];
    //   fArray[1] = v.fArray[1];
    //   return *this;
    // }
    Bool_t isInf() const {
      return (std::isinf(fArray[0]) || std::isinf(fArray[1]) ||
              std::isnan(fArray[0]) || std::isnan(fArray[1]));
    }

    Double_t& val() { return fArray[0]; }
    Double_t& err() { return fArray[1]; }

    Double_t val() const { return fArray[0]; }
    Double_t err() const { return fArray[1]; }

    ValErr& operator/=(const ValErr& v) {
      const Double_t r = val() / v.val();
      err() = r * TMath::Sqrt(relErr2() + v.relErr2());
      val() = r;
      return *this;
    }
    ValErr& operator*=(const ValErr& v) {
      const Double_t p = val() * v.val();
      err() = p * TMath::Sqrt(relErr2() + v.relErr2());
      val() = p;
      return *this;
    }
    Double_t* GetArray() { return fArray; }
    const Double_t* GetArray() const { return fArray; }
  protected:
    Double_t relErr() const { return val() ? err()/val() : 0.0; }
    Double_t relErr2() const { return relErr()*relErr(); }
  private:
    Double_t fArray[2];
  } ;

  typedef std::map<TString, ValErr> branchMap_t;

  AliVdMTree()
    : TObject()
    , fTree(nullptr)
    , fDefaultBranchData()
    , fBranchMap() {}

  virtual ~AliVdMTree() {
    if (fTree)
      fTree->ResetBranchAddresses();
     SafeDelete(fTree);
  }

  TTree* GetTTree() { return fTree; }

  const branchMap_t& GetBranchMap() const { return fBranchMap; }
  branchMap_t& GetBranchMap()  { return fBranchMap; }

  // construction of fTree and its default branches
  TTree* CreateDefaultBranches(TString treeName) {
    SafeDelete(fTree); // pedantic
    fTree = new TTree(treeName, treeName);
    fTree->SetDirectory(nullptr);
    fDefaultBranchData.Branch(fTree);
    return fTree;
  }
  // fills the default branches of fTree
  void FillDefaultBranches(const DefaultBranchData &defBranchData) {
    fDefaultBranchData = defBranchData;
    fTree->Fill();
  }

  // this method iterates over all entries in fTree and adds a new branch
  template<typename FUNC_ITER>
  void AddBranch(const char* branchName, const FUNC_ITER& funcIter) {
    const std::vector<std::string> bs{branchName};
    AddBranch(bs, funcIter);
  }
  template<typename FUNC_ITER>
  void AddBranch(const std::vector<std::string>& branchNames, const FUNC_ITER& funcIter) {
    std::vector<TBranch*> branches;
    for (const std::string& bName : branchNames) {
      TBranch *b = fTree->Branch(bName.c_str(), fBranchMap[bName].GetArray(), "val/D:err");
      branches.push_back(b);
    }
    for (Long64_t i=0, n=fTree->GetEntries(); i<n; ++i) {
      fTree->GetEntry(i);
      funcIter(fDefaultBranchData, fBranchMap);
      for (TBranch *b : branches)
        b->Fill();
    }
  }

  // Zip(this, t2, t3)
  template<typename FUNC_ITER>
  void Zip3(const FUNC_ITER& funcIter, AliVdMTree& t2, AliVdMTree& t3) {
    // assert(fTree->GetEntries() == t2.fTree->GetEntries());
    // assert(fTree->GetEntries() == t3.fTree->GetEntries());
    for (Long64_t i=0, n=fTree->GetEntries(); i<n; ++i) {
      fTree->GetEntry(i);
      t2.GetEntry(i);
      t3.GetEntry(i);
      funcIter(fDefaultBranchData, GetBranchMap(), t2.GetBranchMap(), t3.GetBranchMap());
    }
  }
  void GetEntry(Long64_t iEntry) {
    if (fTree)
      fTree->GetEntry(iEntry);
  }

protected:
private:
  TTree            *fTree;              //
  DefaultBranchData fDefaultBranchData; // default branch data
  branchMap_t       fBranchMap;         // holds named ValErr branch data
  ClassDef(AliVdMTree, 1);
} ;

#endif // PWGPP_VDM_ALIVDMTREE_H
