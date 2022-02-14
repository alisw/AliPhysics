///
/// \file AliFemtoUser/AliFemtoKtBinnedCorrFunc.h
///

#pragma once

#ifndef ALIFEMTOKTBINNEDCORRFUNC_H_
#define ALIFEMTOKTBINNEDCORRFUNC_H_

#include "AliFemtoCorrFctn.h"
#include <TString.h>
#include <vector>

class TH1D;


/// \class AliFemtoKtBinnedCorrFunc
/// \brief A "wrapper class" that wraps a correlation function
///        with code that checks the kt of the pair.
///
class AliFemtoKtBinnedCorrFunc : public AliFemtoCorrFctn {
public:

  /// Construct around a correlation function
  ///
  /// The supplied correlation function object will be cloned for
  /// each bin added
  AliFemtoKtBinnedCorrFunc(const TString& name, AliFemtoCorrFctn*);

  /// Destructor
  ///
  /// Deletes the prototype correlation function
  ///
  virtual ~AliFemtoKtBinnedCorrFunc();

  /// Add a k_{T} bin to the correlation function
  ///
  /// This clones the prototype correlation function and adds a new
  /// range.
  /// It returns the id number of the correlation function.
  ///
  /// Note: This sorts ranges based on the lower bound, so if there
  /// is a range contained within another range, the order inserted
  /// matters.
  UInt_t AddKtRange(float low, float high);

  /// Add continuous ranges from an array of floating points.
  ///
  /// This assumes a "continuous" range, where each element of the
  /// array is the beginning of the next range.
  /// The stop_at parameter is a sentinal value to tell the loop
  /// when to stop (default -1).
  ///
  /// eg.
  ///  cf->AddKtRanges({0, 0.1, 0.2, 0.5, -1.0})
  ///  // adds ranges: {0, 0.1}, {0.1, 0.2}, {0.2, 0.5}
  ///
  void AddKtRanges(const float data[], float stop_at = -1.0);

  /// Add low-high kt ranges in vector
  template <typename Iterator>
  void AddKtRanges(const Iterator &start, const Iterator &stop)
    {
      for (Iterator it = start; it != stop; ++it) {
        AddKtRange(it->first, it->second);
      }
    }

  /// Add low-high kt ranges in vector
  void AddKtRanges(const std::vector<std::pair<float, float>> &data)
  {
    AddKtRanges(data.begin(), data.end());
  };

  /// Add a named k_{T} bin to the correlation function
  ///
  /// This clones the prototype correlation function and adds a new
  /// range.
  /// It returns the id number of the kt range.
  ///
  /// The provided name parameter will be added to the end of the new
  /// output objects
  UInt_t AddKtRange(const TString& name, float low, float high);

  /// Selects the appropriate correlation function and forwards the
  /// pair to the corresponding method.
  virtual void AddRealPair(AliFemtoPair *pair);
  virtual void AddMixedPair(AliFemtoPair *pair);

  /// Given a pair, return the index of the matching correlation
  /// function in the buffer.
  ///
  /// If the pair's kT is outside of any of the stored ranges,
  /// return the NPos constant.
  UInt_t FindKtBin(const AliFemtoPair *pair);

  /// Return TList of output objects
  virtual TList* GetOutputList();
  virtual void AddOutputObjectsTo(TCollection &);

  virtual AliFemtoString Report()
    { return ""; }

  virtual void Finish()
    { }

  virtual AliFemtoCorrFctn* Clone() const
    { return new AliFemtoKtBinnedCorrFunc(*this); }

  /// Constant to return if no correlation is found
  static const UInt_t NPos = static_cast<UInt_t>(-1);

  virtual void EventBegin(const AliFemtoEvent *event);
  virtual void EventEnd(const AliFemtoEvent *event);

private:
  AliFemtoKtBinnedCorrFunc(const AliFemtoKtBinnedCorrFunc&);
  AliFemtoKtBinnedCorrFunc& operator=(const AliFemtoKtBinnedCorrFunc&);

protected:

  void AddPair(AliFemtoPair *, bool);

  /// Name of the output TObjArray
  TString fName;

  /// Prototype Correlation function which will be Clone()'d to create the
  /// objects in each bin.
  AliFemtoCorrFctn *fPrototypeCF;

  /// The vector of correlation functions
  std::vector<AliFemtoCorrFctn*> fCFBuffer;

  /// Internal vector of ranges.
  std::vector<std::pair<Float_t, Float_t> > fRanges;

  /// Histogram monitoring kT distribution
  TH1D* fKtMonitor;
};


#endif /* end of include guard: ALIFEMTOKTBINNEDCORRFUNC_H_ */
