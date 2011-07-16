#ifndef ALISIGNALPROCESOR_H
#define ALISIGNALPROCESOR_H

#include <TObject.h>
class TF1;

class AliSignalProcesor : public TObject{
 public: 
  TF1 * GetAsymGauss();
  void SplineSmoother(const Double_t *ampin, Double_t *ampout, Int_t n) const;
  void TailCancelationALTRO(const Double_t *ampin, Double_t *ampout, Float_t K, Float_t L, 
			    Int_t n) const;
  void TailCancelationTRD(const Double_t *ampin, Double_t *ampout, Float_t r, Float_t c, 
			  Int_t n) const;
  void TailCancelationALTRO1(Double_t *ampin, Double_t *ampout, Float_t norm, Float_t lambda, 
			   Int_t n);

  void TailCancelationTRD1(Double_t *ampin, Double_t *ampout, Float_t norm, Float_t lambda, 
			Int_t n);

  void TailCancelationMI(const Double_t *ampin, Double_t *ampout, Float_t norm, Float_t lambda, 
		       Int_t n) const;

  void TailMaker(const Double_t *ampin, Double_t *ampout, Float_t lambda, 
	       Int_t n) const;

  void TailMakerSpline(const Double_t *ampin, Double_t *ampout, Float_t lambda, 
	       Int_t n) const;
  ClassDef(AliSignalProcesor,1)
};

#endif
