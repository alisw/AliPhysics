#ifndef ALITPCCLUSTERINFO_H
#define ALITPCCLUSTERINFO_H

/// \class AliTPCclusterInfo
/// \brief TPC Cluster Class
///
/// Information for debugging puposes
///
/// \author Marian Ivanov

/* $Id$ */


#include "TObject.h"
//_____________________________________________________________________________
class AliTPCclusterInfo : public TObject {
public:
  AliTPCclusterInfo();
  AliTPCclusterInfo(Bool_t withGraph);
  AliTPCclusterInfo(Float_t *matrix, Int_t nbins, Float_t* graph);
  AliTPCclusterInfo(const  AliTPCclusterInfo & info);
  AliTPCclusterInfo& operator=(const AliTPCclusterInfo& info);
  virtual ~AliTPCclusterInfo();
  UChar_t GetNPads() const { return fNPads;}
  UChar_t GetNTimeBins() const { return fNTimeBins;}
  UChar_t GetNPads(Float_t threshold) const;
  UChar_t GetNTimeBins(Float_t threshold) const;
  Float_t* GetMatrix(){ return fMatrix;}
  void    SetGraph(Float_t * graph, Int_t nbins){ fGraph = graph; fNBins=nbins;}
  void    SetNPadsTimes(UChar_t npads, UChar_t ntimes){ fNPads = npads; fNTimeBins = ntimes;}
 protected:
  Float_t  fMatrix[25];               ///< matrix of amplitude arround center pad - time
  UChar_t  fNPads;                ///< number of pads      in cluster
  UChar_t  fNTimeBins;            ///< number of time bins in cluster
  Int_t  fNBins;                ///< number of bins in graph
  /// signal time dependence graph
  Float_t *fGraph;               //[fNBins]
private:
  /// \cond CLASSIMP
  ClassDef(AliTPCclusterInfo,1)  // Time Projection Chamber clusters Inofrmation
  /// \endcond
};




#endif


