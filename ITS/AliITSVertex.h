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
        virtual ~AliITSVertex();
		  void Exec();
        Double_t PhiFunc(Float_t p[]);

//      This class determines 3D vertex position, resolution and signal 
//      to noise ratio, for arbitrary location (x,y,z) of the vertex.
//      Tests have been carried out with vertex locations up to radial 
//      distances of 10 mm in the transverse plane and up to 15 cm along z.
//      The procedure has been tested also in case of high magnetic fields
//      in ALICE, up to B = 0.5 T.

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
