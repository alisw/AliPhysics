#include "AliPHOSCpvPedProducer.h"
#include "AliPHOSCpvParam.h"
#include "AliPHOSCpvRawStream.h"
#include <fstream>
#include <iostream>
#include <TTree.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TSystem.h>
#include <TTimeStamp.h>

#include "TFile.h"

using namespace std;

using std::ifstream;
using std::ofstream;
ClassImp(AliPHOSCpvPedProducer) ;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliPHOSCpvPedProducer::AliPHOSCpvPedProducer(Int_t sigcut):
  TObject(),
  fSigCut(sigcut),
  fTurbo(kTRUE),
  fhErrors(0),
  fRawStream(0)
{
  //
  //constructor
  //
  for(Int_t iDDL=0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL++) {//iDDL
    fPedMeanMap[iDDL]=0;
    fPedSigMap [iDDL]=0;
    f1DPedMean [iDDL]=0;
    f1DPedSigma[iDDL]=0;
    fPermanentBadMap[iDDL]=0x0;
    for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++)
      for(Int_t iY=1; iY<AliPHOSCpvParam::kPadPcY; iY++)
	fPadAdc[iDDL][iX][iY]=0;
  }//iDDL

  CreateErrHist();
}  //constructor
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliPHOSCpvPedProducer::~AliPHOSCpvPedProducer()
{
  //
  //destructor
  //
  for(Int_t iDDL=0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL++) {//iDDL
    delete fPedMeanMap[iDDL];
    delete fPedSigMap [iDDL];
    delete f1DPedMean [iDDL];
    delete f1DPedSigma[iDDL];
    delete fPermanentBadMap[iDDL];
    for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++)
      for(Int_t iY=1; iY<AliPHOSCpvParam::kPadPcY; iY++)
	delete fPadAdc[iDDL][iX][iY];
  }//iDDL

  //delete fhErrors;
}  //destructor
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvPedProducer::SetPermanentBadMap(TH2* badMap, int iDDL = 0)
{
  if(badMap!=0x0){
    if(iDDL>=0&&iDDL<2*AliPHOSCpvParam::kNDDL){
      fPermanentBadMap[iDDL] = (TH2I*)badMap->Clone();
    }
    else cout<<"DDL number "<<iDDL<<" is not valid"<<endl;
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvPedProducer::SetTurbo(Bool_t turbo)
{
  fTurbo = turbo;
  if(fRawStream) fRawStream->SetTurbo(fTurbo);
}
//--------------------------------------------------------------------------------------
Bool_t AliPHOSCpvPedProducer::LoadNewEvent(AliRawReader *& rawReader)
{
  if(fRawStream) delete fRawStream;
  fRawStream = new AliPHOSCpvRawStream(rawReader);
  if(fRawStream) {
    fRawStream->SetTurbo(fTurbo);
    return kTRUE;
  }
  fhErrors->Fill(0);
  return kFALSE;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvPedProducer::FillPedestal(Int_t abspad,Float_t q)
{
  //
  //Called from the CpvdaFillPedestal() and fills the pedestal values
  //Arguments: absolute pad number as from AliPHOSCpvParam and q-charge
  //
  if(q<0) {
   AliError("Negative charge is read!!!!!!");
   return kFALSE;
  }
  if(AliPHOSCpvParam::IsValidAbs(abspad) && q>0) {
    Int_t iDDL=AliPHOSCpvParam::A2DDL(abspad),
            iX=AliPHOSCpvParam::A2X(abspad),
            iY=AliPHOSCpvParam::A2Y(abspad);
    if(!fPadAdc [iDDL][iX][iY]) CreateDDLHistos(iDDL); 
    fPadAdc [iDDL][iX][iY] -> Fill(q);
    return kTRUE;
  }
  return kFALSE;
}//FillPedestal(int,float)
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvPedProducer::FillPedestal()
{
  //
  //Called from the Cpvda
  //
  while(fRawStream->Next()){
    for(Int_t iPad=0;iPad<fRawStream->GetNPads();iPad++) {
      Int_t charge = fRawStream->GetChargeArray()[iPad];
      Int_t aPad = fRawStream -> GetPadArray()[iPad];
      if(charge){
	if(!AliPHOSCpvParam::IsValidAbs(aPad)) continue;
	if(!FillPedestal(aPad, (Float_t)charge)) return kFALSE;
      }
    }
  }
  return kTRUE;
}//FillPedestal(TClonesArray*)
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliPHOSCpvPedProducer::CalcPedestal(Int_t iDDL)
{
  //
  //Calculate pedestal for each pad
  //Arguments: nDDL-DDL number, nEv - number of the read events
  //Retutns: kTRUE/kFALSE
  //
  //cout<<"Now we are going to calculate pedestals"<<endl;

  if(fPedMeanMap[iDDL]){
    for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++) {
      for(Int_t iY=0; iY<AliPHOSCpvParam::kPadPcY; iY++) {
	//cout<<"Ped["<<iX<<"]["<<iY<<"] = " << fPadAdc[iDDL][iX][iY]->GetMean()<<endl;
	if(fPermanentBadMap[iDDL]!=0x0){
	  if(fPermanentBadMap[iDDL]->GetBinContent(iX+1,iY+1)>0){//bad channel
	    fPedMeanMap[iDDL] -> Fill(iX, iY, fMaxThr);
	    fPedSigMap [iDDL] -> Fill(iX, iY, 0);
	    f1DPedMean [iDDL] -> Fill(fMaxThr);
	    f1DPedSigma[iDDL] -> Fill(0);
	  }
	  else{
	    fPedMeanMap[iDDL] -> Fill(iX, iY, fPadAdc[iDDL][iX][iY]->GetMean());
	    fPedSigMap [iDDL] -> Fill(iX, iY, fPadAdc[iDDL][iX][iY]->GetRMS ());
	    f1DPedMean [iDDL] -> Fill(fPadAdc[iDDL][iX][iY]->GetMean());
	    f1DPedSigma[iDDL] -> Fill(fPadAdc[iDDL][iX][iY]->GetRMS ());
	  }
	}
	else{
	fPedMeanMap[iDDL] -> Fill(iX, iY, fPadAdc[iDDL][iX][iY]->GetMean());
	fPedSigMap [iDDL] -> Fill(iX, iY, fPadAdc[iDDL][iX][iY]->GetRMS ());
	f1DPedMean [iDDL] -> Fill(fPadAdc[iDDL][iX][iY]->GetMean());
	f1DPedSigma[iDDL] -> Fill(fPadAdc[iDDL][iX][iY]->GetRMS ());
	}
      }

    }
    return kTRUE;
  }
  else return kFALSE;
}//CaclPedestal()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvPedProducer::WriteAllHistsToFile(const char * name) const
{
  TFile * rootF = TFile::Open(name,"RECREATE");
  printf("Root file created \n");
  //rootF->cd();
  for(Int_t iDDL=0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL++) {
    // for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++) {
    //   for(Int_t iY=0; iY<AliPHOSCpvParam::kPadPcY; iY++) {
    // 	//fPadAdc[iDDL][iX][iY]->Write();
    //   }
    // }
	//Printf("iDDL = %d\n", iDDL);
    if ( fPedMeanMap[iDDL])
      rootF->WriteObject(fPedMeanMap[iDDL], Form("fPedMeanMap%d",iDDL));
    if ( fPedSigMap[iDDL])
      rootF->WriteObject(fPedSigMap [iDDL], Form("fPedSigMap%d",iDDL));
    if ( f1DPedMean[iDDL])
      rootF->WriteObject(f1DPedMean [iDDL], Form("f1DPedMean%d",iDDL));
    if ( f1DPedSigma[iDDL])
      rootF->WriteObject(f1DPedSigma[iDDL], Form("f1DPedSig%d",iDDL));
    //printf("Write here something \n");
  }
  //if(fhErrors) fhErrors -> Write();

  for(Int_t iDDL=0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL++)
    for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++)
      for(Int_t iY=0; iY<AliPHOSCpvParam::kPadPcY; iY++);
  //fPadAdc[iDDL][iX][iY]->Write();

  rootF->Close();
} //WriteAllHistsToFile()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvPedProducer::WritePedFiles(Int_t iDDL) const
{
  //
  // Write pedestal files to load them to RCB card
  // One file per each column controler
  //

  // accordance to RCB format, pedestals must be written in blocks of 64 lines (not AliPHOSCpvParam::kNPadAdd!)
  Int_t block = 64;

  //cout<<"pedestal files now will be created!"<< endl;

  if(!fPedMeanMap[iDDL]) {
    Printf("No pedestals found for DDL %d !\n");
    return;
  }

  for(Int_t iCC=0; iCC<AliPHOSCpvParam::kNRows; iCC++) {
    FILE * pedFile;
    FILE* pedFileForRCB;
    pedFile = fopen(Form("thr%d_%02d.dat",iDDL,iCC),"w");
    if(!pedFile) {
      Printf("AliPHOSCpvPedProducer::WritePedFiles: Error, file thr%d_%02d.dat could not be open",iDDL,iCC);
    }
    // create and initialize arrays for ped and sigmas
    Int_t ped[AliPHOSCpvParam::kN3GAdd][block],
      sig[AliPHOSCpvParam::kN3GAdd][block];
    for(Int_t i3g=0; i3g<AliPHOSCpvParam::kN3GAdd; i3g++) {
      for(Int_t iPad=0; iPad<block; iPad++) {
	ped[i3g][iPad] = fMaxThr;
	sig[i3g][iPad] = 0;
      }
    }
    Int_t iXmin, iXmax;
    AliPHOSCpvParam::GetLimOfCConX(iCC,iXmin,iXmax);
    //cout<<iXmin<<iXmax<<endl;
    for(Int_t iY=0; iY<AliPHOSCpvParam::kPadPcY; iY++) {
      Int_t g3 = AliPHOSCpvParam::Y23G(iY);
      for(Int_t iX=iXmin; iX<=iXmax; iX++) {
	Int_t pad = AliPHOSCpvParam::XY2Pad(iX,iY);
	if(fPermanentBadMap[iDDL]!=0x0 ){
	  if(fPermanentBadMap[iDDL]->GetBinContent(iX+1,iY+1)>0){
	    ped[g3][pad] = fMaxThr;
	    sig[g3][pad] = 0;
	  }
	  else{
	    ped[g3][pad] = (Int_t) fPadAdc[iDDL][iX][iY]->GetMean();
	    sig[g3][pad] = (Int_t) fPadAdc[iDDL][iX][iY]->GetRMS ();
	  }
	}
	  else{
	    ped[g3][pad] = (Int_t) fPadAdc[iDDL][iX][iY]->GetMean();
	    sig[g3][pad] = (Int_t) fPadAdc[iDDL][iX][iY]->GetRMS ();
	  }
	//cout<< "ped is " << fPadAdc[iDDL][iX][iY]->GetMean()<<endl;
      }
    }

    // write arrays to file
    for(Int_t i3g=0; i3g<AliPHOSCpvParam::kN3GAdd; i3g++) {
      for(Int_t pad=0; pad<block; pad++) {
	// first 10 bit for pedestal
	// last 9 bit for fSigCut*sigma
	Int_t write,writeRCB;
	if ((ped[i3g][pad] + fSigCut*sig[i3g][pad]) > fMaxThr)
	  write = (fMaxThr<<9) + 0;
	else
	  write = ((ped[i3g][pad]+fSigCut * sig[i3g][pad])<<9) + fSigCut * sig[i3g][pad];

	fprintf(pedFile, "0x%05x\n", write);
      }
    }
    fclose(pedFile);
  } // iCC
} // WritePedFiles(iDDL)
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliPHOSCpvPedProducer::CreateErrHist()
{
  Int_t nErrors = AliPHOSCpvRawStream::GetNErrors();
  const char * errNames[nErrors];
  for(Int_t i=0; i<nErrors; i++) {
    errNames[i] = AliPHOSCpvRawStream::GetErrName(i);
  }
  fhErrors = new TH1I("errorTypes","Errors occured during processing",nErrors+1,0,nErrors+1);
  TAxis* x = fhErrors->GetXaxis();
  x->SetBinLabel(1, "Can't get event");
  for(Int_t i=0; i<nErrors; i++) {
    x->SetBinLabel(i+2,errNames[i]);
  }
}
//--------------------------------------------------------------------------------------
void AliPHOSCpvPedProducer::CreateDDLHistos(Int_t iDDL)
{
  // creating histograms
  fPedMeanMap[iDDL] = new TH2F(Form("hPedMeanMap%d",iDDL),Form("2D pedestal value map, DDL = %d",iDDL) ,AliPHOSCpvParam::kPadPcX,0.,AliPHOSCpvParam::kPadPcX,AliPHOSCpvParam::kPadPcY,0.,AliPHOSCpvParam::kPadPcY);
  fPedSigMap [iDDL] = new TH2F(Form("hPedSigMap%d" ,iDDL),Form("2D pedestal sigma map, DDL = %d",iDDL),AliPHOSCpvParam::kPadPcX,0.,AliPHOSCpvParam::kPadPcX,AliPHOSCpvParam::kPadPcY,0.,AliPHOSCpvParam::kPadPcY);
  f1DPedMean [iDDL] = new TH1F(Form("h1DPedMean%d" ,iDDL),Form("pedestal value distribution, DDL = %d",iDDL) ,5000,0,5000);
  f1DPedSigma[iDDL] = new TH1F(Form("h1DPedSigma%d",iDDL),Form("pedestal sigma distribution, DDL = %d",iDDL),1000 ,0,100 );

  // initialization of arrays
  int adr;
  for(Int_t iX=0; iX<AliPHOSCpvParam::kPadPcX; iX++) {
    for(Int_t iY=0; iY<AliPHOSCpvParam::kPadPcY; iY++) {
      adr = AliPHOSCpvParam::XY2A(iDDL,iX,iY);
      fPadAdc[iDDL][iX][iY] = new TH1F(Form("hPad%d_%d_%d",iDDL,iX,iY),Form("Amplitudes CC=%d, 3Gass=%d, pad=%d",AliPHOSCpvParam::A2CC(adr),AliPHOSCpvParam::A23G(adr),AliPHOSCpvParam::A2Pad(adr)),5000,0,5000);
    }//iY
  }//iX

}
