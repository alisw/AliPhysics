#ifndef ALIFMDSURVEYTOALIGNOBJS_H
#define ALIFMDSURVEYTOALIGNOBJS_H
#include <AliSurveyToAlignObjs.h>

// Forward decl
class TVector3;


class AliFMDSurveyToAlignObjs : public AliSurveyToAlignObjs
{
public:
  AliFMDSurveyToAlignObjs() : AliSurveyToAlignObjs() {}
  void Run();
  Bool_t CreateAlignObjs() { return kTRUE; }
protected:
  Bool_t DoFMD1();
  Bool_t GetFMD1Plane(Double_t* rot, Double_t* trans) const;
  Bool_t DoFMD2();
  Bool_t GetFMD2Plane(Double_t* rot, Double_t* trans) const;

  Double_t GetUnitFactor() const;
  Bool_t   GetPoint(const char* name, TVector3& p, TVector3& e) const;
  Bool_t   CalculatePlane(const     TVector3& a, 
			  const     TVector3& b,
			  const     TVector3& c, 
			  Double_t* trans,
			  Double_t* rot) const;
  static void PrintVector(const char* text, const Double_t* v);
  static void PrintVector(const char* text, const TVector3& v);
  static void PrintRotation(const char* text, const Double_t* rot);

  
  ClassDef(AliFMDSurveyToAlignObjs,0) // Convert FMD survey to alignments
};


#endif
//____________________________________________________________________
//
// Local Variables:
//  mode: C++
// End:
//

