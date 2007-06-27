#ifndef AliTPCCALIBTRACKS_H
#define AliTPCCALIBTRACKS_H


#include <TNamed.h>


class AliTPCClusterParam;
class TTreeSRedirector;
class AliTPCROC;
class AliTPCseed;
class AliESDtrack;
class TH3F;
class TH1F;
class TH1I;


class AliTPCcalibTracks : public TNamed {
public :

   // List of branches

   AliTPCcalibTracks();
   virtual ~AliTPCcalibTracks() {;}
   virtual void    ProofSlaveBegin(TList * output);
   virtual void    ProcessTrack(AliTPCseed * seed);
   //
   //
   //
   Float_t         TPCBetheBloch(Float_t bg);
   Bool_t          AcceptTrack(AliTPCseed * track); //
   void            FillHistoCluster(AliTPCseed * track);
   void            FillResolutionHistoLocal(AliTPCseed * track);
   void            AlignUpDown(AliTPCseed * track, AliESDtrack *esd);
   static  Int_t   GetBin(Float_t q,Int_t pad);
   static  Int_t   GetBin(Int_t  iq,Int_t pad);
   static  Float_t GetQ(Int_t bin);
   static  Float_t GetPad(Int_t bin){return bin%3;}

private:
   AliTPCClusterParam *fClusterParam; //pointer to cluster parameterization
   TTreeSRedirector   *fDebugStream;  //debug stream for 
   TList          *fOutput;            //output list 	
   //
   TObjArray *     fArrayAmpRow;//array with amplitudes versus row for given sector 
   TObjArray *     fArrayAmp;   //array with amplitude for sectors
   TObjArray *     fArrayQDY;   //q binned delta Y histograms
   TObjArray *     fArrayQDZ;   //q binned delta Z histograms 
   TObjArray *     fArrayQRMSY;   //q binned delta Y histograms
   TObjArray *     fArrayQRMSZ;   //q binned delta Z histograms 
   TH1F      *     fDeltaY;      // integrated delta y histo
   TH1F      *     fDeltaZ;      // integrated delta z histo
   TObjArray *     fResolY;      // array of resolution histograms Y
   TObjArray *     fResolZ;      // array of resolution histograms Z
   TObjArray *     fRMSY;        // array of RMS histograms Y
   TObjArray *     fRMSZ;        // array of RMS histograms Z
   //   
   TH1I *fHclus;             //!   
   AliTPCROC *fROC;          //!
   Int_t fNRows;             //!
   Int_t fNSect;             //!  
   Int_t fFileNo;            //!
       
   ClassDef(AliTPCcalibTracks,1);
};



#endif
