/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 **************************************************************************/


#include "AliLHCData.h"
#include "TMap.h"


ClassImp(AliLHCData)

const Char_t* AliLHCData::fgkDCSNames[] = {
  "lhcMon_LHCIntensityBeam%d.totalIntensity",
  "lhcInst_BeamIntensityPerBunchBeam%d.averageBeamIntensity",
  "lhcInst_BeamIntensityPerBunchBeam%d_Avgerage.averageBunchIntensities",
  "lhcMon_LHCLumAverageBRANB_4%c2.meanLuminosity",
  "lhcInst_BeamLuminosityPerBunchBRANB_4%c2_Average.bunchByBunchLuminosity",
  "lhcMon_LHCLumAverageBRANB_4%c2.meanCrossingAngle",
  "lhcMon_LHCCirculatingBunchConfig_Beam%d.value",
  "lhcMon_LHCFillNumber.payload",
  //
  "lhcMon_LHCBeamSizeBeam%d.planeSet%d",
  "lhcMon_LHCBeamSizeBeam%d.amplitudeSet%d",
  "lhcMon_LHCBeamSizeBeam%d.positionSet%d",
  "lhcMon_LHCBeamSizeBeam%d.sigmaSet%d"
};

const Char_t* AliLHCData::fgkDCSColNames[] = {
  "lhcMon_LHCCollimatorPos_TCTVB_4L2.%s",
  "lhcMon_LHCCollimatorPos_TCTVB_4R2.%s",  
  "lhcMon_LHCCollimatorPos_TCLIA_4R2.%s"
};

const Char_t* AliLHCData::fgkDCSColJaws[] = {
  "lvdt_gap_downstream","lvdt_gap_upstream","lvdt_left_downstream",
  "lvdt_left_upstream","lvdt_right_downstream","lvdt_right_upstream"};

//___________________________________________________________________
AliLHCData::AliLHCData(const TMap* dcsMap, double tmin, double tmax, int avPeriod)
  : fPeriod(avPeriod),fTMin(tmin),fTMax(tmax)  
{
  FillData(dcsMap,tmin,tmax);
}

//___________________________________________________________________
Bool_t AliLHCData::FillData(const TMap* dcsMap, double tmin, double tmax)
{
  // process DCS map and fill all fields. 
  // Accept only entries with timestamp between tmin and tmax
  const double ktReal = 1200000000.;
  const double kCollTolerance = 100e-4; // tolerance on collimator move (cm)
  char buff[100];
  TObjArray* arr;
  AliDCSArray* dcsVal;
  Double_t tPeriodEnd=0;
  Int_t dcsSize,nEntries,iEntry;
  //
  SetTMin(tmin);
  SetTMax(tmax);
  //
  // -------------------------- extract Fill Number
  arr = GetDCSEntry(dcsMap,fgkDCSNames[kRecFillNum],iEntry,tmin,tmax);
  if (arr && iEntry>=0) SetFillNumber( ((AliDCSArray*)arr->At(iEntry))->GetInt(0) );
  //
  // -------------------------- extract total intensities for beam 1 and 2 (DC BCT: slow average)
  for (int ibm=0;ibm<2;ibm++) {
    //
    sprintf(buff,fgkDCSNames[kRecTotInt],ibm+1);  
    if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,tmin,tmax)) || iEntry<0 ) continue;
    AliLHCDipValD* curVal;
    tPeriodEnd = 0;
    //
    nEntries = arr->GetEntriesFast();
    while (iEntry<nEntries) {
      dcsVal = (AliDCSArray*) arr->At(iEntry++);
      double tstamp = dcsVal->GetTimeStamp();
      if (tstamp>tmax) break;
      if (tstamp>tPeriodEnd) {
	curVal = new AliLHCDipValD(1,0.,0);  // start new period
	fIntTot[ibm].Add(curVal);
	// if tmin is provided, count periods from it, otherwise from the 1st timestamp
	if (tPeriodEnd<1) tPeriodEnd = ((tmin>ktReal) ? tmin : tstamp);
	tPeriodEnd += fPeriod;
      }
      *curVal += *dcsVal;
    }
    for (int i=fIntTot[ibm].GetEntries();i--;) ((AliLHCDipValD*)(fIntTot[ibm])[i])->Average();
  }
  //
  // -------------------------- extract total intensities for beam 1 and 2 (BCTFR: fast average)
  for (int ibm=0;ibm<2;ibm++) {
    //
    sprintf(buff,fgkDCSNames[kRecTotIntBunch],ibm+1);  
    if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,tmin,tmax)) || iEntry<0 ) continue;
    AliLHCDipValD* curVal;
    tPeriodEnd = 0;
    //
    nEntries = arr->GetEntriesFast();
    while (iEntry<nEntries) {
      dcsVal = (AliDCSArray*) arr->At(iEntry++);
      double tstamp = dcsVal->GetTimeStamp();
      if (tstamp>tmax) break;
      if (tstamp>tPeriodEnd) {
	curVal = new AliLHCDipValD(1,0.,0);  // start new period
	fIntTotBunch[ibm].Add(curVal);
	// if tmin is provided, count periods from it, otherwise from the 1st timestamp
	if (tPeriodEnd<1) tPeriodEnd = ((tmin>ktReal) ? tmin : tstamp);
	tPeriodEnd += fPeriod;
      }
      *curVal += *dcsVal;
    }
    for (int i=fIntTotBunch[ibm].GetEntries();i--;) ((AliLHCDipValD*)(fIntTotBunch[ibm])[i])->Average();
  }
  //  
  // -------------------------- extract total luminosities according L and R detectors
  for (int ilr=0;ilr<2;ilr++) {
    //
    sprintf(buff,fgkDCSNames[kRecTotLum],ilr ? 'L':'R');  
    if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,tmin,tmax)) || iEntry<0 ) continue;
    AliLHCDipValD* curVal = 0;
    tPeriodEnd = 0;
    //
    nEntries = arr->GetEntriesFast();
    while (iEntry<nEntries) {
      dcsVal = (AliDCSArray*) arr->At(iEntry++);
      double tstamp = dcsVal->GetTimeStamp();
      if (tstamp>tmax) break;
      if (tstamp>tPeriodEnd) {
	curVal = new AliLHCDipValD(1,0.,0);  // start new period
	fLuminTot[ilr].Add(curVal);
	// if tmin is provided, count periods from it, otherwise from the 1st timestamp
	if (tPeriodEnd<1) tPeriodEnd = ((tmin>ktReal) ? tmin : tstamp);
	tPeriodEnd += fPeriod;
      }
      *curVal += *dcsVal;
    }
    for (int i=fLuminTot[ilr].GetEntries();i--;) ((AliLHCDipValD*)(fLuminTot[ilr])[i])->Average();
  }
  //
  // -------------------------- extract mean crossing angles according to L and R detectors
  for (int ilr=0;ilr<2;ilr++) {
    //
    sprintf(buff,fgkDCSNames[kRecCrossAngle],ilr ? 'L':'R');  
    if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,tmin,tmax)) || iEntry<0 ) continue;
    AliLHCDipValD* curVal=0;
    tPeriodEnd = 0;
    //
    nEntries = arr->GetEntriesFast();
    while (iEntry<nEntries) {
      dcsVal = (AliDCSArray*) arr->At(iEntry++);
      double tstamp = dcsVal->GetTimeStamp();
      if (tstamp>tmax) break;
      if (tstamp>tPeriodEnd) {
	curVal = new AliLHCDipValD(1,0.,0);  // start new period
	fCrossAngle[ilr].Add(curVal);
	// if tmin is provided, count periods from it, otherwise from the 1st timestamp
	if (tPeriodEnd<1) tPeriodEnd = ((tmin>ktReal) ? tmin : tstamp);
	tPeriodEnd += fPeriod;
      }
      *curVal += *dcsVal;
    }
    for (int i=fCrossAngle[ilr].GetEntries();i--;) ((AliLHCDipValD*)(fCrossAngle[ilr])[i])->Average();
  }
  //
  //
  // ------------------------- extract bunch configuration for beam 1 and 2
  int nbunch[2];
  for (int ibm=0;ibm<2;ibm++) {
    //
    nbunch[ibm] = -1;
    sprintf(buff,fgkDCSNames[kRecBunchConf],ibm+1);  
    if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,tmin,tmax)) || iEntry<0 ) continue;
    dcsVal = (AliDCSArray*) arr->At(iEntry);
    nEntries = dcsVal->GetNEntries();     // count number of actual bunches
    nbunch[ibm] = 0;
    while(nbunch[ibm]<nEntries && dcsVal->GetInt(nbunch[ibm])) nbunch[ibm]++;
    if (!nbunch[ibm]) {
      AliWarning(Form("Beam%d bunches configuration record is present but empty",ibm+1));
      continue;
    }
    fBunchConfig[ibm].SetSize(nbunch[ibm]);
    fBunchConfig[ibm] += *dcsVal;
  }
  // -------------------------- extract intensities per bunch for beam 1 and 2
  for (int ibm=0;ibm<2;ibm++) {
    //
    sprintf(buff,fgkDCSNames[kRecBunchInt],ibm+1);  
    if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,tmin,tmax)) || iEntry<0 ) continue;
    AliLHCDipValD* curVal=0;
    tPeriodEnd = 0;
    //
    dcsVal = (AliDCSArray*)arr->At(iEntry);
    nEntries = dcsVal->GetNEntries();     // count number of actual bunches
    dcsSize = 0;
    while(dcsSize<nEntries && dcsVal->GetDouble(dcsSize)>0) dcsSize++;
    if (!dcsSize) {
      AliWarning(Form("Beam%d bunch intensities record is present but empty",ibm+1));
      continue;
    }
    if (nbunch[ibm]>0) { // bunch pattern was provided
      if (dcsSize>nbunch[ibm]) {
 	AliWarning(Form("Beam%d declares %d bunches but %d bunch intensities are non-0. Take first %d",
			ibm+1,nbunch[ibm],dcsSize,nbunch[ibm]));
	dcsSize = nbunch[ibm];
      }
      else if (dcsSize<nbunch[ibm]) {
 	AliWarning(Form("Beam%d declares %d bunches but %d bunch intensities are non-0. Skip",
			ibm+1,nbunch[ibm],dcsSize));
	continue;
      }
    }
    //
    nEntries = arr->GetEntriesFast();
    while (iEntry<nEntries) {
      dcsVal = (AliDCSArray*) arr->At(iEntry++);
      double tstamp = dcsVal->GetTimeStamp();
      if (tstamp>tmax) break;
      if (tstamp>tPeriodEnd) {
	curVal = new AliLHCDipValD(dcsSize,0.,0);  // start new period
	fIntBunch[ibm].Add(curVal);
	// if tmin is provided, count periods from it, otherwise from the 1st timestamp
	if (tPeriodEnd<1) tPeriodEnd = ((tmin>ktReal) ? tmin : tstamp);
	tPeriodEnd += fPeriod;
      }
      *curVal += *dcsVal;
    }
    for (int i=fIntBunch[ibm].GetEntries();i--;) ((AliLHCDipValD*)(fIntBunch[ibm])[i])->Average();
  }
  //
  // -------------------------- extract per bunch luminosities according L and R detectors
  for (int ilr=0;ilr<2;ilr++) {
    //
    sprintf(buff,fgkDCSNames[kRecBunchLum],ilr ? 'L':'R');  
    if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,tmin,tmax)) || iEntry<0 ) continue;
    AliLHCDipValD* curVal=0;
    tPeriodEnd = 0;
    //
    dcsVal = (AliDCSArray*) arr->At(iEntry);
    nEntries = dcsVal->GetNEntries();     // count number of actual bunches
    dcsSize = 0;
    while(dcsSize<nEntries && dcsVal->GetDouble(dcsSize)>0) dcsSize++;
    if (!dcsSize) {
      AliWarning(Form("Probe%c bunch luminosities record is present but empty",ilr ? 'R':'L'));
      continue;
    }
    //
    if (nbunch[ilr]>0) { // bunch pattern was provided
      if (dcsSize>nbunch[ilr]) {
 	AliWarning(Form("Beam%d declares %d bunches but %d bunch luminosities are non-0. Take first %d",
			ilr+1,nbunch[ilr],dcsSize,nbunch[ilr]));
	dcsSize = nbunch[ilr];
      }
      else if (dcsSize<nbunch[ilr]) {
 	AliWarning(Form("Beam%d declares %d bunches but %d bunch luminosities are non-0. Skip",
			ilr+1,nbunch[ilr],dcsSize));
	continue;
      }
    }
    //
    nEntries = arr->GetEntriesFast();
    while (iEntry<nEntries) {
      dcsVal = (AliDCSArray*) arr->At(iEntry++);
      double tstamp = dcsVal->GetTimeStamp();
      if (tstamp>tmax) break;
      if (tstamp>tPeriodEnd) {
	curVal = new AliLHCDipValD(dcsSize,0.,0);  // start new period
	fLuminBunch[ilr].Add(curVal);
	// if tmin is provided, count periods from it, otherwise from the 1st timestamp
	if (tPeriodEnd<1) tPeriodEnd = ((tmin>ktReal) ? tmin : tstamp);
	tPeriodEnd += fPeriod;
      }
      *curVal += *dcsVal;
    }
    for (int i=fLuminBunch[ilr].GetEntries();i--;) ((AliLHCDipValD*)(fLuminBunch[ilr])[i])->Average();
  }
  //
  // ------------------------- extract gaussian fit params for beam 1 and 2 profiles
  for (int ibm=0;ibm<2;ibm++) {
    for (int ixy=0;ixy<2;ixy++) {
      // determine which projection corresponds actually to given ixy
      sprintf(buff,fgkDCSNames[kRecPrfPrID],ibm+1,ixy+1); 
      if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,0/*tmin*/,tmax)) || iEntry<0 ) continue;
      dcsVal = (AliDCSArray*) arr->At(iEntry);
      int proj = dcsVal->GetInt(0)-1;          // beam projection
      //
      if (proj!=kX && proj!=kY) {
	AliError(Form("Unknown beam projection %d for %s",proj,buff));
	continue;
      }
      // Amp,Pos,Sig - each one have separate entry and time stamp (though come from the single fit)
      int entPar[3];
      TObjArray *arrPar[3];
      AliDCSArray *dcsPar[3];
      AliLHCDipValD* curVal=0;
      double tstamp = 0;
      int npars = 0;
      for (int ipar=0;ipar<3;ipar++) {
	sprintf(buff,fgkDCSNames[ipar+kRecPrfAmp],ibm+1,ixy+1);
	if ( !(arrPar[ipar]=GetDCSEntry(dcsMap,buff,entPar[ipar],tmin,tmax)) || entPar[ipar]<0 ) break;
	dcsPar[ipar] = (AliDCSArray*) arrPar[ipar]->At(entPar[ipar]);	
	if (dcsPar[ipar]->GetTimeStamp()>tstamp) tstamp = dcsPar[ipar]->GetTimeStamp(); // max time among 1st entries
	npars++;
      }
      if (npars<3) continue; // ignore incomplete data
      //
      tPeriodEnd = 0;
      // start recording from max timeStamp: 
      // the entries for different params must correspond to same timestamp
      while(1) { 
	//
	// read next timestamp for which all 3 params are present
	npars = 0; // align the first entries to read to same timestamp
	for (int ipar=0;ipar<3;ipar++) {
	  while(entPar[ipar]<arrPar[ipar]->GetEntriesFast()) {
	    dcsPar[ipar] = (AliDCSArray*) arrPar[ipar]->At(entPar[ipar]);
	    double df = dcsPar[ipar]->GetTimeStamp() - tstamp;
	    if (TMath::Abs(df)<0.5) { // same time stamp, ok
	      npars++; 
	      break;
	    }
	    if (df<0) entPar[ipar]++;  // check next entry
	    else {
	      tstamp = dcsPar[ipar]->GetTimeStamp();
	      ipar = -1; // reference tstamp was changed, check all arrays again
	      npars = 0;
	      break;
	    }
	  }
	} // 
	if (npars<3) break; // no more data
	for (int ipar=0;ipar<3;ipar++) entPar[ipar]++;
	//
	if (tstamp>tmax) break;
	if (tstamp>tPeriodEnd) {
	  curVal = new AliLHCDipValD(3,0.,0);  // start new period
	  fBeamPos[ibm][proj].Add(curVal);
	  // if tmin is provided, count periods from it, otherwise from the 1st timestamp
	  if (tPeriodEnd<1) tPeriodEnd = ((tmin>ktReal) ? tmin : tstamp);
	  tPeriodEnd += fPeriod;
	}
	int nsamp = curVal->GetNSamplesUsed()+1;
	curVal->SetTimeStamp( (tstamp + curVal->GetNSamplesUsed()*curVal->GetTimeStamp())/nsamp);
	curVal->SetNSamplesUsed(nsamp);
	for (int ipar=3;ipar--;) (*curVal)[ipar] += dcsPar[ipar]->GetDouble(0);
	//
      }
      //
      for (int i=fBeamPos[ibm][proj].GetEntriesFast();i--;) ((AliLHCDipValD*)(fBeamPos[ibm][proj])[i])->Average();
      //
    } // projection
  } // beam
  //
  // ------------------------- extract collimators data
  for (int icl=0;icl<kNCollimators;icl++) {
    for (int jaw=0;jaw<kNJaws;jaw++) {
      sprintf(buff,fgkDCSColNames[icl],fgkDCSColJaws[jaw]);        
      if ( !(arr=GetDCSEntry(dcsMap,buff,iEntry,tmin,tmax)) || iEntry<0 ) continue;
      dcsVal = (AliDCSArray*) arr->At(iEntry);	      
      AliLHCDipValD* curVal = new AliLHCDipValD(1,dcsVal->GetTimeStamp(),1);
      (*curVal)[0] = dcsVal->GetDouble(0)/10;  // gap in cm
      fCollimators[icl][jaw].Add(curVal);
      //
      // now track the changes above the threshold (100 um?)
      nEntries = arr->GetEntriesFast();
      while(++iEntry<nEntries) {
	dcsVal = (AliDCSArray*) arr->At(iEntry);
	if (dcsVal->GetTimeStamp() > tmax) break;  // out of time
	double val = dcsVal->GetDouble(0)/10;
	if ( TMath::Abs(val-curVal->GetValue(0))<kCollTolerance ) continue;  // no significant change
	// need to create a new entry
	curVal = new AliLHCDipValD(1,dcsVal->GetTimeStamp(),1);
	(*curVal)[0] = val;  // gap in cm
	fCollimators[icl][jaw].Add(curVal);	
      }
    } // jaws
  } // collimators
  //
  return kTRUE;
}

//___________________________________________________________________
TObjArray* AliLHCData::GetDCSEntry(const TMap* dcsMap,const char* key,int &entry,double tmin,double tmax) const
{
  // extract array from the DCS map and find the first entry within the time limits
  TObjArray* arr = (TObjArray*)dcsMap->GetValue(key);
  if (!arr || !arr->GetEntriesFast()) { 
    AliWarning(Form("No data for %s",key)); 
    return 0;
  }
  int ntot = arr->GetEntriesFast();
  for (entry=0;entry<ntot;entry++) {
    AliDCSArray* ent = (AliDCSArray*)arr->At(entry);
    if (ent->GetTimeStamp()>=tmin && ent->GetTimeStamp()<=tmax) break;
  }
  if (entry==ntot) {
    entry = -1;
    TString str;
    str += AliLHCDipValD::TimeAsString(tmin);
    str += " : ";
    str += AliLHCDipValD::TimeAsString(tmax);
    AliWarning(Form("All entries for %s are outside the requested range:\n%s",key,str.Data()));
  }
  return arr;
}

//___________________________________________________________________
void AliLHCData::Print(const Option_t* opt) const
{
  // print everything
  printf("LHC DIP Data | Fill Number#%d (Averaging time: %d sec.)\n",GetFillNumber(),fPeriod);
  printf("Validity period: %s : %s\n\n",
	 fTMin<1.23e9 ? "N/A": AliLHCDipValD::TimeAsString(fTMin),
	 fTMax>7.00e9 ? "N/A": AliLHCDipValD::TimeAsString(fTMax) );
  //
  int n=0;
  for (int ibm=0;ibm<2;ibm++) {
    printf("*** Bunch Configuration for Beam%d: %s\n",ibm,(n=fBunchConfig[ibm].GetSize()) ? "":"N/A");
    if (n) fBunchConfig[ibm].Print(opt);
  }
  printf("\n");
  //
  for (int ibm=0;ibm<2;ibm++) {
    printf("*** Average total intensity for Beam%d (DCBCT): %s\n",ibm,(n=fIntTot[ibm].GetEntriesFast()) ? "":"N/A");
    for (int i=0;i<n;i++) (fIntTot[ibm])[i]->Print(opt);
  }
  printf("\n");
  //
  for (int ibm=0;ibm<2;ibm++) {
    printf("*** Total intensity for Beam%d (BCTFR): %s\n",ibm,(n=fIntTotBunch[ibm].GetEntriesFast()) ? "":"N/A");
    for (int i=0;i<n;i++) (fIntTotBunch[ibm])[i]->Print(opt);
  }
  printf("\n");  //
  for (int ibm=0;ibm<2;ibm++) {
    printf("*** Bunch intensities for Beam%d: %s\n",ibm,(n=fIntBunch[ibm].GetEntriesFast()) ? "":"N/A");
    for (int i=0;i<n;i++) (fIntBunch[ibm])[i]->Print(opt);
  }
  printf("\n");
  //
  for (int ibm=0;ibm<2;ibm++) {
    printf("*** Total luminosity for probe%c: %s\n",ibm ? 'R':'L',(n=fLuminTot[ibm].GetEntriesFast()) ? "":"N/A");
    for (int i=0;i<n;i++) (fLuminTot[ibm])[i]->Print(opt);
  }
  printf("\n");
  //
  for (int ibm=0;ibm<2;ibm++) {
    printf("*** Bunch luminosities for probe%c: %s\n",ibm ? 'R':'L',(n=fLuminBunch[ibm].GetEntriesFast()) ? "":"N/A");
    for (int i=0;i<n;i++) (fLuminBunch[ibm])[i]->Print(opt);
  }
  printf("\n");
  //
  for (int ibm=0;ibm<2;ibm++) {
    printf("*** Crossing angle for probe%c: %s\n",ibm ? 'L':'R',(n=fCrossAngle[ibm].GetEntriesFast()) ? "":"N/A");
    for (int i=0;i<n;i++) (fCrossAngle[ibm])[i]->Print(opt);
  }
  printf("\n");
  //
  for (int ibm=0;ibm<2;ibm++) {
    for (int ixy=0;ixy<2;ixy++) {
      printf("*** Gaussian fit for Beam%d %c profile: %s\n",ibm,ixy ? 'Y':'X',(n=fBeamPos[ibm][ixy].GetEntriesFast()) ? "":"N/A");
      for (int i=0;i<n;i++) (fBeamPos[ibm][ixy])[i]->Print(opt);
    }
  }
  //
  for (int icl=0;icl<kNCollimators;icl++) {
    printf("\n");
    for (int ij=0;ij<kNJaws;ij++) {
      printf(fgkDCSColNames[icl],fgkDCSColJaws[ij]);
      printf(": %s\n",(n=fCollimators[icl][ij].GetEntriesFast()) ? "":"N/A");
      for (int i=0;i<n;i++) (fCollimators[icl][ij])[i]->Print(opt);
    }
  }
  //
}
