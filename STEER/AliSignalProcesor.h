#ifndef ALISIGNALPROCESOR_H
#define ALISIGNALPROCESOR_H

class AliSignalProcesor : public TObject{
 public: 
  TF1 * GetAsymGauss();
  void SplineSmoother(Double_t *ampin, Double_t *ampout, Int_t n);
  void TailCancelationALTRO(Double_t *ampin, Double_t *ampout, Float_t K, Float_t L, 
			    Int_t n);
  void TailCancelationTRD(Double_t *ampin, Double_t *ampout, Float_t r, Float_t c, 
			  Int_t n);
  void TailCancelationALTRO1(Double_t *ampin, Double_t *ampout, Float_t norm, Float_t lambda, 
			   Int_t n);

  void TailCancelationTRD1(Double_t *ampin, Double_t *ampout, Float_t norm, Float_t lambda, 
			Int_t n);

  void TailCancelationMI(Double_t *ampin, Double_t *ampout, Float_t norm, Float_t lambda, 
		       Int_t n);

  void TailMaker(Double_t *ampin, Double_t *ampout, Float_t lambda, 
	       Int_t n);

  void TailMakerSpline(Double_t *ampin, Double_t *ampout, Float_t lambda, 
	       Int_t n);
  ClassDef(AliSignalProcesor,1)
};

#endif
