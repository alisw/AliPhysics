#ifndef ALIITSGLOBALRECPOINT_H
#define ALIITSGLOBALRECPOINT_H

class AliITSglobalRecPoint : public TObject {

	public:
			
		AliITSglobalRecPoint();
      AliITSglobalRecPoint(Double_t gx, Double_t gy, Double_t gz, Double_t gsx, Double_t gsy, Double_t gsz, Int_t l);
		
		virtual ~AliITSglobalRecPoint() { }
		
		// difference in azymuthal angle
		Double_t DPhi  (AliITSglobalRecPoint *p); 
		// difference in polar angle
		Double_t DTheta(AliITSglobalRecPoint *p); 
		// checks if the pt contains a given positive label
		Bool_t   HasID (Int_t ID)      
		{ if (ID < 0) return kFALSE; else return (fLabel[0]==ID || fLabel[1]==ID || fLabel[2]==ID); }
		// checks if the pt shares at least a positive label with another one
		Bool_t   SharesID(AliITSglobalRecPoint *p) 
		{ return (HasID(p->fLabel[0]) || HasID(p->fLabel[1]) || HasID(p->fLabel[2])); }
		
		// Parameters for sorting
		Bool_t  IsSortable() const { return kTRUE; }
		Int_t   Compare(const TObject *O) const;
		
	public:
			
		Double_t fGX, fGY, fGZ; // (x,y,z) in the global reference
			
		Double_t fR2;    // = sqrt(x^2 + y^2)
		Double_t fR3;    // = sqrt(x^2 + y^2 + z^2)
		Double_t fPhi;   // = atan(y/x)
		Double_t fTheta; // = acos(z/r3)
		
		Double_t fGSX; //
		Double_t fGSY; // sigmas of global coords
		Double_t fGSZ; //
		
		Int_t fModule;      // ITS module containing the point
		Int_t fPosInModule; // position in TClonesArray of TreeR
		Int_t fLayer;       // ITS layer containing the point
		Int_t fLabel[3];    // Generated tracks containing the point
		Int_t fKalmanLabel; // Kalman recognized track containing the point
		
		Int_t fSign;        //   0 if the point doesn't belong to a Kalman track
		                    // + 1 if the point belongs to a Kalman good track
		                    // - 1 if the point belongs to a Kalman fake track
		
		Int_t fUsed;        // a flag to avoid point recycling
		
		ClassDef(AliITSglobalRecPoint, 1) // AliITSglobalRecPoints class
};

#endif
