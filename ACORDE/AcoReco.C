//_____________________________________________________//
//                                                     //
//     This macro reads ACORDE Raw-Rootfied files      //
//     and draws Histograms of Acorde's Bitpattern     //    
//                                                     //
//     Author:                                         //
//                                                     //
//   Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> //
//                                                     //
//     Also comments to:                               //
//                                                     //
//   Mario Sitta <sitta@to.infn.it>                    //
//   Arturo Fernandez <afernan@mail.cern.ch>           //
//   Luciano Diaz <luciano.diaz@nucleares.unam.mx>     //
//   Eleazar Cuautle <ecuautle@nucleares.unam.mx>      //
//____________________________________________________ //


#include <TClonesArray.h>
#include <stdio.h>

//void AcoReco(char* fileName="08000020614001.20.root")
void AcoReco(char* fileName)
{
  TStopwatch timer;
  timer.Start();
  AliRawReader* rawReader = new AliRawReaderRoot(fileName); 
  AliRawReader* rCount = new AliRawReaderRoot(fileName);
  AliACORDERawStream* rawStream  = new AliACORDERawStream(rawReader);    
  TH1D *h1 = new TH1D("h1","ACORDE - Single Muon Hits",60,1,60);
  TH1D *h2 = new TH1D("h2","ACORDE - Hit Multiplicity",60,1,60);
  size_t contSingle=0;
  size_t contMulti=0;
  UInt_t dy[4];
  bool kroSingle[60],kroMulti[60];
  UInt_t tmpDy;
  Int_t dyma[60],DyM,nEvents=0;
  for(Int_t m=0;m<60;m++) {kroSingle[m]=0;kroMulti[m]=0;dyma[m]=0;}
  
 
  Int_t nEvents = rawStream->GetNEvents(fileName);
  
  printf("Numero de eventos: %d \n",nEvents);

for (Int_t i=1; i<=nEvents; i++) 
{

//       printf("No Event %d",i);	
       if (!rawReader->NextEvent()) break;
	rawStream->Reset();
      // printf("No Event %d",i);
        
	if (!rawStream->Next()) continue;
	dy[0]=rawStream->GetWord(0);
	dy[1]=rawStream->GetWord(1);
	dy[2]=rawStream->GetWord(2);
	dy[3]=rawStream->GetWord(3);
	tmpDy=dy[0];
	for(Int_t r=0;r<30;++r)
	{
		kroSingle[r] = tmpDy & 1;
		tmpDy>>=1;
	}
	tmpDy=dy[1];
	for(Int_t r=30;r<60;++r)
	{
		kroSingle[r] = tmpDy & 1;
		tmpDy>>=1;
	}
	tmpDy=dy[2];
	for(Int_t r=0;r<30;++r)
	{
		kroMulti[r] = tmpDy & 1;
		tmpDy>>=1;
	}
	tmpDy=dy[3];
	for(Int_t r=30;r<60;++r)
	{
		kroMulti[r] = tmpDy & 1;
		tmpDy>>=1;
	}
	contSingle=0;	
	for(Int_t r=0;r<60;++r) 
	{
		if(kroSingle[r]==1) 
		{
			dyma[r]=dyma[r]+1;
			h1->Fill(r+1);
			contSingle=contSingle+1;
		}
		
	}h2->Fill(contSingle);
}

TCanvas *acorde = new TCanvas("ACORDE","ACORDE-Hist-Real Data",1);
acorde->Divide(2,1);
acorde->cd(1);
h1->GetXaxis()->SetTitle("No. of Channel");
h1->GetYaxis()->SetTitle("No. of Hits");
h1->SetFillColor(kRed);
h1->Draw();
acorde->cd(2);
gPad->SetLogy();
h2->GetXaxis()->SetTitle("No. of Modules");
h2->GetYaxis()->SetTitle("Multiplicity");
h2->SetFillColor(kBlue);
h2->Draw();

delete rawReader;
delete rawStream;
timer.Stop();
timer.Print();
}
