#ifndef ALIITSVERTEX_H
#define ALIITSVERTEX_H

class TTree;
class TFile;
class AliITSgeom;
class AliITSRecPoint;
class TH1F;
class TF1;
class TClonesArray;
class TObject;
class AliGenerator;

class AliITSVertex : public TObject  {
 
 public:
 
	AliITSVertex();
	~AliITSVertex();
	Double_t PhiFunc(Float_t p[]);

//	At present this class determines vertex position, resolution and signal 
//	to noise ratio only for the z coordinate. For x and y coordinates it 
//	gives the values x and y setted in Config.C with resolution and signal 
//      to noise ratio values = 0.
//      The cases of beam off-set and magnetic field = 0.4 T are included.

	Double_t GetZv() {return (Double_t)fPosition[2];}
	Double_t GetZRes() {return fResolution[2];}
	Double_t GetZSNR() {return fSNR[2];}
	Double_t GetYv() {return (Double_t)fPosition[1];}
	Double_t GetYRes() {return fResolution[1];}
	Double_t GetYSNR() {return fSNR[1];}
	Double_t GetXv() {return (Double_t)fPosition[0];}
	Double_t GetXRes() {return fResolution[0];}
	Double_t GetXSNR() {return fSNR[0];}
        
 private:
    
	Double_t *fPosition;
	Double_t *fResolution;
	Double_t *fSNR;

ClassDef(AliITSVertex,1) // Class for Vertex finder
};

#endif
