#ifndef ALITESTAD_H
#define ALITESTAD_H

#include "AliDisplacedVertexSelectionAD.h"
#include "AliAnalysisTaskSE.h"
class TList;
class TH2;

class AliTestAD : public AliAnalysisTaskSE
{
public:
  AliTestAD()
    : AliAnalysisTaskSE(),
      fList(0),
      fAD(),
      fIPz0(0),
      fIPz1(0),
      fIPz2(0),
      fIPzSat0(0),
      fIPzSat1(0),
      fIPzSat2(0),
      fIPzSatP(0),
      fIsInit(false)
  {}
  AliTestAD(const char*); //  : AliAnalysisTaskSE("ad"), fList(0) {}
  AliTestAD(const AliTestAD& o)
    : AliAnalysisTaskSE(o),
      fList(0),
      fAD(o.fAD),
      fIPz0(0),
      fIPz1(0),
      fIPz2(0),
      fIPzSat0(0),
      fIPzSat1(0),
      fIPzSat2(0),
      fIPzSatP(0),
      fIsInit(false)
  {}
  ~AliTestAD() {}
  AliTestAD& operator=(const AliTestAD&) { return *this; }
  void UserCreateOutputObjects();
  void UserExec(Option_t* option);
  void Connect();
private:
  TList* fList; //!
  AliDisplacedVertexSelectionAD fAD;
  TH2*   fIPz0;    //!
  TH2*   fIPz1;    //!
  TH2*   fIPz2;    //!
  TH1*   fIPzSat0; //! 
  TH1*   fIPzSat1; //! 
  TH1*   fIPzSat2; //! 
  TH1*   fIPzSatP; //! 
  Bool_t fIsInit;
  ClassDef(AliTestAD,1);
};

#endif
//
// Local Variables:
//  mode: C++
// End:
//
