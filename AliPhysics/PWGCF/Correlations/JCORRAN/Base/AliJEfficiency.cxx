#include "AliJEfficiency.h"
#include <TSystem.h>
#include <iostream>
#include <TGrid.h>

// AliJEfficiency
// ...
// ...
// ...
// ...
// TODO

using namespace std;

AliJEfficiency::AliJEfficiency():
  fMode(kAuto),
  fPeriod(-1),
  fTrackCut(),
  fRunTable(),
  fDataPath(""),
  fName(""),
  fPeriodStr(""),
  fMCPeriodStr(""),
  fRunNumber(0),
  fTag(""),
  fInputRootName(""),
  fInputRoot(NULL),
  fCentBin(0x0)
{
  for (int i=0; i<3; i++) fEffDir[i] = NULL;
}

AliJEfficiency::AliJEfficiency(const AliJEfficiency& obj) :
  fMode(obj.fMode),
  fPeriod(obj.fPeriod),
  fTrackCut(obj.fTrackCut),
  fRunTable(obj.fRunTable),
  fDataPath(obj.fDataPath),
  fName(obj.fName),
  fPeriodStr(obj.fPeriodStr),
  fMCPeriodStr(obj.fMCPeriodStr),
  fRunNumber(obj.fRunNumber),
  fTag(obj.fTag),
  fInputRootName(obj.fInputRootName),
  fInputRoot(obj.fInputRoot),
  fCentBin(obj.fCentBin)
{
  // copy constructor TODO: handling of pointer members
  JUNUSED(obj);
  for (int i=0; i<3; i++) fEffDir[i] = obj.fEffDir[i];
}

AliJEfficiency& AliJEfficiency::operator=(const AliJEfficiency& obj){
  // equal sign operator TODO: content
  JUNUSED(obj);
  return *this;
}

TString AliJEfficiency::GetEffName() {
  /*
     1. kNotUse : no Load, efficiency is 1 always
     2. has fInputRootName : Load that or crash
     3. has fName : Load fName [+runnumber] or crash
     4. has runnumber : Find Good MC period from AliJRunTable, or crash
     3. has period : Find Good MC period from AliJRunTable, or crash

*/
  if( fMode == kNotUse ) {
    cout<<"J_WARNING : Eff Mode is \"NOTUSE\". eff is 1 !!!"<<endl;
    return "";
  }

  //==== 1. fInputRootName
  //==== 2. Eff-fName-fPeriod-fMCPeriod-fRunNumber-fTag.root
  if( fInputRootName.Length() == 0 ){

    //==== SELECT Period : fPeriod, fPeriodStr
    // 1. Use fPeriodStr if it is
    // 2. Use fPeriod    if it is
    // 3. Use RunNumber
    if( fPeriodStr.Length() == 0 ){
      if( fRunNumber > 0  && fPeriod < 0 ) {
        fPeriod = fRunTable.GetRunNumberToPeriod( fRunNumber );
      }
      fPeriodStr = fRunTable.GetPeriodName( fPeriod );
    }
    //==== Select McPeriod
    //==== 1. Use fMCPeriodStr
    //==== 2. Use fPeriod
    if( fMCPeriodStr.Length() == 0 ){
      fMCPeriodStr =  fRunTable.GetMCPeriod( fPeriod );
    }
    //==== Select fRunNumber
    //==== MODE 1 : runnumber = 0;
    if( fMode == kPeriod ) fRunNumber = 0;
    else if( fRunNumber < 0 ){
      cout<< "J_ERROR : Runumber must be >=0. Input runnumber is "<<fRunNumber<<endl;
      gSystem->Exit(1);
    }

    fInputRootName = Form("Eff-%s-%s-%s-%d-%s.root", fName.Data(), fPeriodStr.Data(), fMCPeriodStr.Data(), int(fRunNumber), fTag.Data() );
  }

  return fInputRootName;
}

TString AliJEfficiency::GetEffFullName() {
  GetEffName();
  fInputRootName = fDataPath + "/" + fInputRootName;
  return fInputRootName;
}


bool AliJEfficiency::Load(){
  // Load Efficiency File based on fMode
  if( fMode == kNotUse ) {
    cout<<"J_WARNING : Eff Mode is \"NOTUSE\". eff is 1 !!!"<<endl;
    return true;
  }
  GetEffFullName();
  if (TString(fInputRootName).BeginsWith("alien:"))  TGrid::Connect("alien:");
  fInputRoot = TFile::Open( fInputRootName);
  //fInputRoot = new TFile( fInputRootName,"READ");
  if( !fInputRoot ) {
	  cout<< "J_ERROR : "<<fInputRootName <<" is not exist"<<endl;
	  gSystem->Exit(1);
      return 0;
  }

  //fEffDir[0] = (TDirectory*)fInputRoot->Get("EffRE");
  ///fEffDir[1] = (TDirectory*)fInputRoot->Get("EffMC");
  fEffDir[2] = (TDirectory*)fInputRoot->Get("Efficiency");
  //if( fEffDir[0] && fEffDir[1] && fEffDir[2] )
  if( !fEffDir[2] )
  {
	  cout<< "J_ERROR : Directory EFF is not exist"<<endl;
	  gSystem->Exit(1);
  }

  fCentBin = (TAxis*)fEffDir[2]->Get("CentralityBin");
  if( !fCentBin ){
	  cout<< "J_ERROR : No CentralityBin in directory"<<endl;
	  gSystem->Exit(1);
      return 0;
  }


  int nVtx = 1;
  int nCentBin = fCentBin->GetNbins();
  for( int ivtx=0;ivtx<nVtx;ivtx++ ){
	  for( int icent=0;icent<nCentBin;icent++ ){
		  for( int icut=0;icut<fTrackCut.GetNCut();icut++ ){
			  fCorrection[ivtx][icent][icut] 
				  = (TGraphErrors*) fEffDir[2]->Get(Form("gCor%02d%02d%02d", ivtx,icent,icut));
			  //cout<<"J_LOG : Eff graph - "<<Form("gCor%02d%02d%02d", ivtx,icent,icut)<<" - "<<g<<endl;
		  }
	  }
  }
  cout<<"J_LOG : Eff file is "<<fInputRootName<<endl;
  cout<<"J_LOG : Eff Cent Bins are ";
  for( int i=0;i<=nCentBin;i++ ){
	  cout<<fCentBin->GetXbins()->At(i)<<" ";
  }
  cout<<endl;
  return true;
}

double AliJEfficiency::GetCorrection( double pt, int icut , double cent ) const {
	// TODO : Function mode
	if( fMode == kNotUse ) return 1;
	int icent = fCentBin->FindBin( cent ) -1 ;
	if( icent < 0 || icent > fCentBin->GetNbins()-1 ) {
		cout<<"J_WARNING : Centrality "<<cent<<" is out of CentBinBorder"<<endl;
		return 1;
	}
	// TODO error check for icent;
	int ivtx = 0;
	if( ! fCorrection[ivtx][icent][icut] ) {
		cout<<"J_WARNING : No Eff Info "<<pt<<"\t"<<icut<<"\t"<<cent<<"\t"<<icent<<endl;
		return 1;
	}
	TGraphErrors * gr = fCorrection[ivtx][icent][icut];
	//=== TEMPERORY SETTING. IT will be removed soon.
	if( pt > 30 ) pt = 30; // Getting eff of 30GeV for lager pt
	double cor = gr->Eval(pt);
	if ( cor < 0.2 ) cor = 0.2;
	return cor;
}

void AliJEfficiency::Write(){
	// Write Efficiency information to root file 
	if( fMode == kNotUse ){
		cout<<"J_LOG : Efficiency mode is \"NotUse\", nothing will be Written" <<endl;
		return;
	}
	cout<<"J_LOG : Efficiency Write to "<<gDirectory->GetName()<<endl;
	TDirectory *cwd = gDirectory;
	TDirectory *eff = gDirectory->mkdir("Efficiency");
	eff->cd();
	int nVtx = 1;
	int nCentBin = fCentBin->GetNbins();
	cout<<nCentBin<<endl;
	for( int ivtx=0;ivtx<nVtx;ivtx++ ){
		for( int icent=0;icent<nCentBin;icent++ ){
			for( int icut=0;icut<fTrackCut.GetNCut();icut++ ){
				cout<<fCorrection[ivtx][icent][icut]<<endl;
				if( !fCorrection[ivtx][icent][icut]) continue;
				fCorrection[ivtx][icent][icut]->Write(Form("gCor%02d%02d%02d", ivtx,icent,icut));
			}
		}
	}
	fCentBin->Write("CentralityBin");
	cwd->cd();
}

