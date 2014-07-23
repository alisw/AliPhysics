#ifndef ALIITSUTRACKERSAAUX_H
#define ALIITSUTRACKERSAAUX_H

#include <vector>
using std::vector;
#include <algorithm>
using std::sort;
#include "AliExternalTrackParam.h"
#include "AliITSUTrackCooked.h"
#include "AliITSUAux.h"

#include <TMath.h>

#include <Riostream.h>
using std::cout;
using std::endl;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template <typename T> 
struct Comparison { //Adapted from TMath ROOT code
  Comparison(vector<T> *d) : fData(d) {}

  bool operator()(int i1, int i2) {
    return fData->at(i1).CmpVar() < fData->at(i2).CmpVar();
  }
  vector<T> *fData;
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
struct Cell {
  
  Cell() : Level(1),
           Param(),
           Points(),
           N(0),
           Neighbours(),
           Label(-1){ 
      // Inlined standard constructor
  } ;

  Cell(int i1, int i2) : Level(1),
                         Param(),
                         Points(),
                         N(2),
                         Neighbours(),
                         Label(-1){ 
      // Inlined standard constructor
      Points[0] = i1;
      Points[1] = i2;
  } ;

  Cell(int i1, int i2, int i3) : Level(1),
                                 Param(),
                                 Points(),
                                 N(3),
                                 Neighbours(),
                                 Label(-1) {  
				  
      // Inlined standard constructor
      Points[0] = i1;
      Points[1] = i2;
      Points[2] = i3;
  } ;
  
  int Level;
  float Param[3];
  int Points[3];
  int N;
  vector<int> Neighbours;
  int Label;
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
struct SpacePoint {
  SpacePoint(const float *xyz, const float &alpha) : XYZ(),
                                                     Cov(),
                                                     Phi(),
                                                     Alpha(alpha),
                                                     Label(),
                                                     Used(false) 
  {
      XYZ[0] = xyz[0];
      XYZ[1] = xyz[1];
      XYZ[2] = xyz[2];
      Cov[0]=Cov[1]=Cov[2]=99999.f;
  };

  SpacePoint(const SpacePoint &sp) : XYZ(),
                                     Cov(),
                                     Phi(sp.Phi),
                                     Alpha(sp.Alpha),
                                     Label(),
                                     Used(sp.Used) 
  {
    for(int i=0; i<3; ++i) {
      XYZ[i]=sp.XYZ[i];
      Label[i]=sp.Label[i];
      Cov[i]=sp.Cov[i];
    }
  };

  float CmpVar() const { return Phi; }

  float XYZ[3];
  float Cov[3];
  float Phi;
  float Alpha;
  int Label[3];
  bool Used;
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
struct Layer {
  Layer() : Index(),
            N(0),
            Points() {};

  int operator()(const int &i) { return Index[i]; } // Just for compatibility with old code
  SpacePoint& operator[](const int &i) { return Points[Index[i]]; } 
  ~Layer() { Clear(); }

  void AddPoint(const float xyz[3], const float cov[3], const float &phi, const float &alpha) {
    Index.push_back(N++);
    Points.push_back(SpacePoint(xyz,alpha));
    Points.back().Cov[0] = cov[0];
    Points.back().Cov[1] = cov[1];
    Points.back().Cov[2] = cov[2];
    Points.back().Phi = phi;
  }

  void Sort() {
    Comparison<SpacePoint> cmp(&Points);
    sort(Index.begin(),Index.end(),cmp); 
  }

  void Clear() { 
    Index.clear();
    Points.clear();
    N=0;
  }

  vector<int> Index;
  int N;
  vector<SpacePoint> Points;
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
struct Road {

  Road() : Elements(), N(0), Label(-1) {
    ResetElements(); 
  }

  Road(int layer, int id) : Elements(), N(1), Label(-1) {
    ResetElements();
    Elements[layer] = id;
  }

  Road(const Road& copy) : Elements(), N(copy.N), Label(copy.Label) {
    for ( int i=0; i<5; ++i ) Elements[i] = copy.Elements[i];
  }

  int &operator[] (const int &i) {
    return Elements[i];
  }

  void ResetElements() {
    for ( int i=0; i<5; ++i ) 
      Elements[i] = -1;   
  }

  void AddElement(int i, int el) {
    ++N;
    Elements[i] = el;
  }
  
  int Elements[5];
  int N;
  int Label;
};

// class AliITSUTrackCA : public AliITSUTrackCooked {
// public:
//   AliITSUTrackCA() : AliITSUTrackCooked() {}
//   AliITSUTrackCA(const AliITSUTrackCA &t) : AliITSUTrackCooked((AliITSUTrackCooked)t) {}
//   AliITSUTrackCA(const AliESDtrack &t) : AliITSUTrackCooked(t) {}
//   AliITSUTrackCA &operator=(const AliITSUTrackCooked &t);
//   virtual ~AliITSUTrackCA();
//   void SetChi2(Double_t chi2) { fChi2=chi2; }
//   ClassDef(AliITSUTrackCA,1)
// };

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> 
struct Comparison <AliITSUTrackCooked> {
  Comparison(vector<AliITSUTrackCooked> *d) : fData(d) {}

  bool operator()(int i1, int i2) {
    return (fData->at(i1).GetChi2() < fData->at(i2).GetChi2());
  }
  vector<AliITSUTrackCooked> *fData;
};

// //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// class Candidate : public AliKalmanTrack {
//   public :

//   Candidate() : AliKalmanTrack(), 
//                 fChi2( 0. ), 
//                 fPoints(), 
//                 fNPoints(0), 
//                 fInnermostLayer(-1), 
//                 fOutermostLayer(-1) 
//   {
//     for ( unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i ) fPoints[i]=-1;
//   };

//   Candidate(const Candidate &copy) : AliExternalTrackParam(), 
//                                      fChi2(copy.fChi2), 
//                                      fPoints(),
//                                      fNPoints(copy.fNPoints), 
//                                      fInnermostLayer(copy.fInnermostLayer), 
//                                      fOutermostLayer(copy.fOutermostLayer)
//   {
//     for ( unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i ) fPoints[i]=copy.fPoints[i];
//   };

//   Candidate(int points[7]) : AliExternalTrackParam(), 
//                              fChi2( 0. ), 
//                              fPoints(),
//                              fNPoints( 0 ), 
//                              fInnermostLayer( -1 ), 
//                              fOutermostLayer( -1 )
//   { 
//     bool outer=false;
//     ResetPoints();
//     for ( int i=6; i>=0; --i ) {
//       if (points[i]!=-1) {
//         if (! outer ) { 
//           outer = true;
//           fOutermostLayer = points[i];
//         }
//         fInnermostLayer = points[i];
//         ++fNPoints;
//       }
//       fPoints[i<<0x1] = points[i];
//     }
//   }

//   int GetClusterIndex(const int &i) const {
//     return ((fInnermostLayer+i)<<28)+fPoints[(fInnermostLayer+i)<<0x1];
//   }

//   Double_t CmpVar() const { return fChi2; }

//   void ResetPoints() { for(unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i) fPoints[i]=-1; }
//   //
//   // Double_t fChi2;
//   // Double_t fFakeRatio;
//   Int_t fPoints[2*AliITSUAux::kMaxLayers];
//   // Int_t fNPoints;
//   Int_t fLabel;
//   Int_t fInnermostLayer;
//   Int_t fOutermostLayer;

// };

#endif
