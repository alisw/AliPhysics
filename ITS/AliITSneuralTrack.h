#ifndef ALIITSNEURALTRACK_H
#define ALIITSNEURALTRACK_H

class AliITSglobalRecPoint;

class AliITSneuralTrack : public TObject {
	
	friend class AliITSneuralTracker;

public:
	         AliITSneuralTrack();
		 AliITSneuralTrack(const AliITSneuralTrack &src);
		 AliITSneuralTrack &operator=(const AliITSneuralTrack &src);
	virtual ~AliITSneuralTrack();
	
	Double_t GetSqChi() const {return fSqChi;}
	Double_t GetR(Double_t &Xc, Double_t &Yc) const
		{Xc=fFitXC; Yc=fFitYC; return fFitRadius;}
		
	Int_t  CheckMe(Bool_t verbose) const ;
	Int_t  EvaluateTrack(Bool_t verbose, Int_t min, Int_t* &good);
		
	void GetCoords(Double_t* &x, Double_t* &y, Double_t* &z) const;
	void CopyPoint(AliITSglobalRecPoint *p);
	void Print(Option_t *option, Int_t min);
	void Kinks(Int_t &pos, Int_t &neg, Int_t &decr, Int_t &incr);
			
	AliITSglobalRecPoint* GetPoint(Int_t i) { return fPoint[i]; }
		
private:
	
	Double_t FitXY(Double_t VX, Double_t VY); //!
		
	Double_t fFitXC;  // center x coord.
	Double_t fFitYC;  // center y coord.
	Double_t fFitRadius; // curv. radius
	Double_t fFitTanL; // tangent of L
	Double_t fSqChi; // chi squared
	
	AliITSglobalRecPoint *fPoint[6];  // points on 6 layers
	
	ClassDef(AliITSneuralTrack, 1)
};

#endif
