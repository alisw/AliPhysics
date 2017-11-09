///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoModelWeightGeneratorLednicky : the most advanced weight       //
// generator available. Supports a large number of different pair types  //
// and interaction types. Can calculate pair weights coming from         //
// quantum statistics, coulomb interation and strong interaction ot any  //
// combination of the three, as applicable.                              //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELWEIGHTGENERATORLEDNICKY_H
#define ALIFEMTOMODELWEIGHTGENERATORLEDNICKY_H

#include "AliFemtoTypes.h"
#include "AliFemtoModelWeightGenerator.h"

class AliFemtoModelWeightGeneratorLednicky : public  AliFemtoModelWeightGenerator {
 public: 
// --- Constructor
  AliFemtoModelWeightGeneratorLednicky(); // call SetDefaultCalcPar
  AliFemtoModelWeightGeneratorLednicky(const AliFemtoModelWeightGeneratorLednicky &aWeight); // call SetDefaultCalcPar
// --- Destructor : nothing to explicitly delete
  AliFemtoModelWeightGeneratorLednicky& operator=(const AliFemtoModelWeightGeneratorLednicky& aWeight);
  ~AliFemtoModelWeightGeneratorLednicky();

  virtual Double_t GenerateWeight(AliFemtoPair *aPair);

  virtual void     SetPairType(Int_t aPairType);
  virtual void     SetPairTypeFromPair(AliFemtoPair *aPair);
  virtual Int_t    GetPairType() const;

  virtual Double_t GetKStar() const;
  virtual Double_t GetKStarOut() const;
  virtual Double_t GetKStarSide() const;
  virtual Double_t GetKStarLong() const;
  virtual Double_t GetRStar() const;
  virtual Double_t GetRStarOut() const;
  virtual Double_t GetRStarSide() const;
  virtual Double_t GetRStarLong() const;

  virtual AliFemtoModelWeightGenerator* Clone() const;

// --- Setting

// >>> Calculation mode
  void SetDefaultCalcPar(); // Default is CoulOn, QuantumOn, StrongOn, 3BodyOff, Square, T0ApproxOff
  void SetCoulOn();
  void SetCoulOff();

  void SetQuantumOn();
  void SetQuantumOff();
  void SetStrongOn();
  void SetStrongOff();
  void Set3BodyOn();
  void Set3BodyOff();
  void SetSphere(); // use Spherical wave approximation
  void SetSquare(); // use use Square potential (only for p-p and pi+Pi-) otherwise, use spherical wave approx
  void SetT0ApproxOff();//only with  Spherical wave Approximation - this is default mode
  void SetT0ApproxOn(); 
 
// Test Lambda parameters
  void PrintLambdas(){;}
  
  void SetNuclCharge(const double aNuclCharge); // for 3-body calculation
  void SetNuclMass(const double aNuclMass);

  void SetKpKmModelType(const int aModelType, const int aPhi_OffOn);  // K+K- model type,Phi off/on

  virtual AliFemtoString Report();

protected:
  // Fsi weight output
  double  fWei;  // normal weight
  double  fWein; // weight with nuclear influence
  double  fWeif; // weight
  double  fWeightDen; // weight for the denominator

  // Setting parameters
  int fItest;    // if set to 1 default parameters will be used

  //int mNs;
  int    fIch;        // switch coulomb interaction on/off
  int    fIqs;        // switch quantum statistics on/off
  int    fIsi;        // switch strong interaction on/off
  int    fI3c;        // switch 3rd body influence on/off
  double fNuclMass;   // mass of the third body
  double fNuclCharge; // charge of the third body

  bool   fSphereApp;       // use spherical approximation
  bool   fT0App;           // use square well approximation

  //Pair identification
  int       fLL;             // internal pair type code
  short     fNuclChargeSign; // sign of the 3rd body charge
  bool      fSwap;           // are particle in right order ? 
  int const fLLMax;          // number of supported pairs
  char**    fLLName;         // name of the system
  int *     fNumProcessPair; // number of process pairs of each type
  int       fNumbNonId;      // Number of unidentified pairs

  //K+K- model type
  int       fKpKmModel;      //ij (i=1..4, j=1..4; see AliFemtoFsiWeightLednicky.F)
  int       fPhi_OffOn;      //0->Phi Off,1->Phi On
  int       fNS_4;           //set NS is equal to 4

  // Interface to the fortran functions
  void FsiSetKpKmModelType();  //// initialize K+K- model type
  void FsiInit();
  void FsiSetLL();
  void FsiNucl();
  bool SetPid(const int aPid1,const int aPid2);

#ifdef __ROOT__
  ClassDef(AliFemtoModelWeightGeneratorLednicky,1)
#endif
};

#endif
