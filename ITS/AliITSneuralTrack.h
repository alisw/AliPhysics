#ifndef ALIITSNEURALTRACK_H
#define ALIITSNEURALTRACK_H

class AliITSglobalRecPoint;

class AliITSneuralTrack : public TObject {
	
	friend class AliITSneuralTracker;

public:
	         AliITSneuralTrack();
	virtual ~AliITSneuralTrack();
	
	Double_t GetSqChi() {return fSqChi;}
	Double_t GetR(Double_t &Xc, Double_t &Yc)
		{Xc=fFitXC; Yc=fFitYC; return fFitRadius;}
		
	Int_t  CheckMe(Bool_t verbose);
	Int_t  EvaluateTrack(Bool_t verbose, Int_t min, Int_t* &good);
		
	void GetCoords(Double_t* &x, Double_t* &y, Double_t* &z);
	void CopyPoint(AliITSglobalRecPoint *p);
	void Print(Option_t *option, Int_t min);
	void Kinks(Int_t &pos, Int_t &neg, Int_t &decr, Int_t &incr);
			
	AliITSglobalRecPoint* GetPoint(Int_t i) { return fPoint[i]; }
		
private:
	
	Double_t FitXY(Double_t VX, Double_t VY); //!
		
	Double_t fFitXC;
	Double_t fFitYC;
	Double_t fFitRadius;
	Double_t fFitTanL;
	Double_t fSqChi;
	
	AliITSglobalRecPoint *fPoint[6]; 
	
	ClassDef(AliITSneuralTrack, 1)
};

#endif
