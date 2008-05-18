
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
 **************************************************************************/

/****************************************************************************
 * TPC Calibration Class for GOOFIE values. Drift velocity, gas composition *
 * and the gain.                                                            *
 ****************************************************************************/
 
#include "AliTPCGoofieValues.h"
#include <iostream>

/*****************************************************************************
*  The class AliTPCGoofieValues allows the access to GoofieValues. The only  *
*  construtor is load a data from ASCI file. The methods make Tgraphs and    *
*  TSplines of the time dependace of th values. One method allows save the   *
*  the graphs and spline togather with tree of allvalues into file.          *

Current example usage:

AliTPCGoofieValues *goofieVal = new AliTPCGoofieValues("Goofie_data_january_run_01_08.txt");
TFile f("goofieValues.root","recreate");
goofieVal->Write("goofie");
TBrowser b;
And now you can browse


*****************************************************************************/

ClassImp(AliTPCGoofieValues)

//____________________________________________________________________________
AliTPCGoofieValues::AliTPCGoofieValues():
  TNamed(),
  fLinesInFile(0),///< lines in ASCII file
  fStartTime(0),///< StartTime
  fEndTime(0),///< EndTime
  fTimeOfRun(0),///<TimeOfRun
  fTempGrad(0),///<TempGrad
  fAverageTempGraph(0), ///<graph of average temperatures 
  fTempGradGraph(0),///<graph of temperature gradients
  fPressGraph(0),///<graph of pressures
  fVdriftGraph(0),///<graph of drift velocities
  fVdriftcorGraph(0),///<graph of drift velocities corrected
  fGainFGraph(0),///<graph of near gain
  fGainNGraph(0),///<graph of far gain
  fCO2Graph(0),///<graph of CO2 content
  fN2Graph(0),///<graph of N2 content
  
  fAverageTempSpline(0),///< spline of average temperatures 
  fTempGradSpline(0),///<spline of  temperature gradients
  fPressSpline(0),///<spline of pressures
  fVdriftSpline(0),///<spline of drift velocities
  fVdriftcorSpline(0),///<spline of drift velocities corrected
  fGainFSpline(0),///<spline of near gain
  fGainNSpline(0),///<spline of far gain
  fCO2Spline(0),///<spline of CO2 content
  fN2Spline(0),
  fGoofieValues(0)///<spline of  N2 content
{
  //
  // Default constructor
  //
} 
  
  //____________________________________________________________________________
AliTPCGoofieValues::AliTPCGoofieValues(const char *fname):
  TNamed(),
  fLinesInFile(0),///< lines in ASCII file
  fStartTime(0),///< StartTime
  fEndTime(0),///< EndTime
  fTimeOfRun(0),///<TimeOfRun
  fTempGrad(0),///<TempGrad
  fAverageTempGraph(0), ///<graph of average temperatures 
  fTempGradGraph(0),///<graph of temperature gradients
  fPressGraph(0),///<graph of pressures
  fVdriftGraph(0),///<graph of drift velocities
  fVdriftcorGraph(0),///<graph of drift velocities corrected
  fGainFGraph(0),///<graph of near gain
  fGainNGraph(0),///<graph of far gain
  fCO2Graph(0),///<graph of CO2 content
  fN2Graph(0),///<graph of N2 content
  
  fAverageTempSpline(0),///< spline of average temperatures 
  fTempGradSpline(0),///<spline of  temperature gradients
  fPressSpline(0),///<spline of pressures
  fVdriftSpline(0),///<spline of drift velocities
  fVdriftcorSpline(0),///<spline of drift velocities corrected
  fGainFSpline(0),///<spline of near gain
  fGainNSpline(0),///<spline of far gain
  fCO2Spline(0),///<spline of CO2 content
  fN2Spline(0),
  fGoofieValues(0)///<spline of  N2 content
{
    /**
       Constructor take a values from ASCI file with raw values. <br>
       example: <b>fname = AliTPCGoofie_run_001.txt </b><br>
       Read the next values, in this order (each line of the file is a data point) <br>
       fGoofieTime fGoofieTempF fGoofieTempN fGoofiePress fGoofieVdrift fGoofieVdriftcor fGoofieAreaF fGoofieAreaN fGoofieCO2 fGoofieN2 

     */

    fLinesInFile= 0 ;
    fStartTime= 0 ;
    fEndTime= 0 ;
    fTimeOfRun= 0 ;
    fTempGrad= 0 ;
 
    fGoofieValues = new TTree("tree","goofie values"); 
    fLinesInFile  = fGoofieValues->ReadFile(fname,"fGoofieTime/D:fGoofieTempF/D:fGoofieTempN/D:fGoofiePress/D:fGoofieVdrift/D:fGoofieVdriftcor/D:fGoofieAreaF/D:fGoofieAreaN/D:fGoofieCO2/D:fGoofieN2/D"); 

    fAverageTempGraph =  new TGraph();
    fTempGradGraph    =  new TGraph();
    fPressGraph       =  new TGraph();
    fVdriftGraph      =  new TGraph();
    fVdriftcorGraph   =  new TGraph();
    fGainFGraph       =  new TGraph();
    fGainNGraph       =  new TGraph();
    fCO2Graph         =  new TGraph();
    fN2Graph          =  new TGraph();
 
    FillAllGraphs();
    //Splines are allocated via FillAllSplines
    FillAllSplines();
  }
  //________________________________________________________________________________
  AliTPCGoofieValues::~AliTPCGoofieValues()
  {
    /**
       AliTPCGoofieValues destructor
     */
 
  }
  //_________________________________________________________________________________
  Long64_t AliTPCGoofieValues::GetLinesInFile(){
    return fLinesInFile;
  }
  //_________________________________________________________________________________
  Double_t AliTPCGoofieValues::GetStartTime()
  {
    /**
      take a time in beging of run from Goofie in seconds 
     */
 
    return fStartTime;
  }
  //__________________________________________________________________________________
  Double_t AliTPCGoofieValues::GetEndTime()
  {
    /**
       take the time in the end of run from Goofie in seconds
     */
  
    return fEndTime;
  }
  
  //__________________________________________________________________________________
  Double_t AliTPCGoofieValues::GetTimeOfRun()
  {
    /**
       return time of run in seconds
    */
  
    Double_t time = GetEndTime(/*const char *fname*/) - GetStartTime(/*const char *fname*/);
    return time;
  }
  //__________________________________________________________________________________
  Double_t AliTPCGoofieValues::GetTempGrad(Double_t timeSec)
  {

    /**
       gradient of temperature in the chosen time in run
     */
  
    Double_t tempGrad;
    FillTempGradGraph();
    tempGrad = EvalTempGrad(timeSec);
    return tempGrad;
  }
  //__________________________________________________________________________________
  void AliTPCGoofieValues::FillAllGraphs(){
    /**
       Fill ALL the Graphs.<br>
       There are individual methods, to do the same.<br>
       This function is anyway called in the ctor. 
       Selection has to be implemented:
       //
       Example selection  :
       
       AliTPCGoofieValues *a = new AliTPCGoofieValues("Goofie_data_january_run_01_08.txt");
       TEventList list("listGood","listGood");
       a->GetTree()->Draw(">>listGood","abs(fGoofieAreaN-2600)<500&&abs(fGoofieAreaF-2900)<300&&abs(fGoofieAreaF/fGoofieAreaN-1.1)<0.2","");
       a->GetTree()->Draw("fGoofieVdrift:fGoofieTime");
       a->GetTree()->SetEventList(&list);
       a->GetTree()->SetMarkerColor(2);
       //
       a->GetTree()->Draw("fGoofieVdrift:fGoofieTime","","same*");
       TGraph gr(200,a.GetTree()->GetV2(),a.GetTree()->GetV1());
       //example residuals
       AliSplineFit fit;
       fit.SetGraph(&gr)
       fit->SetMinPoints(200);
       fit->InitKnots(&gr,15,0,0.2)
       fit.SplineFit(1)
       fit.MakeDiffHisto(&gr)->Draw();

       
     */
    Long64_t nevent = fGoofieValues->GetEntries();
    
    //temporal for reading the branches
    Double_t fGoofieTime   = 0; 
    Double_t fGoofieTempN  = 0; 
    Double_t fGoofieTempF  = 0; 
    Double_t fGoofiePress  = 0;   
    Double_t fGoofieAreaF  = 0; 
    Double_t fGoofieAreaN  = 0; 
    Double_t fGoofieVdriftcor  = 0; 
    Double_t fGoofieVdrift     = 0; 
    Double_t fGoofieCO2  = 0; 
    Double_t fGoofieN2   = 0; 
    fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);
    fGoofieValues->SetBranchAddress("fGoofieTempN",&fGoofieTempN); 
    fGoofieValues->SetBranchAddress("fGoofieTempF",&fGoofieTempF);  
    fGoofieValues->SetBranchAddress("fGoofiePress",&fGoofiePress);
    fGoofieValues->SetBranchAddress("fGoofieAreaF",&fGoofieAreaF);
    fGoofieValues->SetBranchAddress("fGoofieAreaN",&fGoofieAreaN);
    fGoofieValues->SetBranchAddress("fGoofieVdriftcor",&fGoofieVdriftcor);
    fGoofieValues->SetBranchAddress("fGoofieVdrift",&fGoofieVdrift);
    fGoofieValues->SetBranchAddress("fGoofieCO2",&fGoofieCO2);
    fGoofieValues->SetBranchAddress("fGoofieN2",&fGoofieN2); 
    fStartTime = 0; //fGoofieValues->GetBranch("fGoofieTime")->GetEntry(0);

    for (int i = 0; i < nevent; ++i){
      if (fGoofieValues->GetEvent(i)<0){
	// cout<< "you're done, man  !" << endl;
	continue;
      }
      if (i==0) fStartTime = fGoofieTime;
      else if (i == (nevent-1)) fEndTime = fGoofieTime;
      //// cout<<  " reading " << i << endl;
      fAverageTempGraph->SetPoint(i,(fGoofieTime - fStartTime),(fGoofieTempF+fGoofieTempN)/2);
      fTempGradGraph->SetPoint(i,(fGoofieTime - fStartTime),(fGoofieTempF-fGoofieTempN)/25); 
      fPressGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofiePress);
      fGainFGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieAreaF);
      fGainNGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieAreaN);
      fVdriftGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieVdrift);  
      fVdriftcorGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieVdriftcor);     
      fCO2Graph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieCO2);
      fN2Graph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieN2);
    }

  }

  void AliTPCGoofieValues::FillAllSplines(){ 
    /**
       Fill ALL the Splines.<br>
       There are individual methods, to do the same.<br>
       This function is anyway called in the ctor. 
     */
   
    fAverageTempSpline = new TSpline3("temperature",fAverageTempGraph);
    fTempGradSpline    = new TSpline3("temperature gradient",fTempGradGraph);
    fPressSpline       = new TSpline3("pressure",fPressGraph);
    fVdriftcorSpline   = new TSpline3("vdriftcor",fVdriftcorGraph);
    fVdriftSpline      = new TSpline3("vdrift",fVdriftGraph);  
    fGainFSpline       = new TSpline3("gainF",fGainFGraph);
    fGainNSpline       = new TSpline3("gainN",fGainNGraph);
    fCO2Spline         = new TSpline3("co2",fCO2Graph);
    fN2Spline          = new TSpline3("n2",fN2Graph);
    
  }

  //___________________________________________________________________________________
  void AliTPCGoofieValues::FillAverageTempGraph() 
  {
    /**
       graph of temperature vs time. <br>
       If the graph is already filled, it's only putting a title.
     */

    fAverageTempGraph->SetTitle("time dependence of temperature");    

    if (!fAverageTempGraph->GetN()==0){
      // cout<< " graph already filled !!! " << endl;
    }
    else {
      Double_t fGoofieTime   = 0; 
      Double_t fGoofieTempN  = 0; 
      Double_t fGoofieTempF  = 0; 
      fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);
      fGoofieValues->SetBranchAddress("fGoofieTempN",&fGoofieTempN); 
      fGoofieValues->SetBranchAddress("fGoofieTempF",&fGoofieTempF);  
      
      Long64_t nevent = fGoofieValues->GetEntries();
      
      for (int i = 0; i < nevent; ++i){
	if (fGoofieValues->GetEvent(i)<0){
	  // cout<< "you're done, man  !" << endl;
	  continue;
	} 
	if (i==0) fStartTime = fGoofieTime;
	else if (i == (nevent-1)) fEndTime = fGoofieTime;
	fAverageTempGraph->SetPoint(i,(fGoofieTime - fStartTime),(fGoofieTempF+fGoofieTempN)/2);
      }
    }
  }
  //__________________________________________________________________________________
  void AliTPCGoofieValues::FillTempGradGraph()
  {
    /**
       graph of temperature gradient [K/cm]<br>
       If the graph is already filled, it's only putting a title.
     */
  
    fTempGradGraph->SetTitle("time dpendance of Temperature`s gradient");
    
    if(! fTempGradGraph->GetN()==0){ 
      // cout<< " graph already filled !!! " << endl;
    }
    else{
      Double_t fGoofieTime   = 0; 
      Double_t fGoofieTempN  = 0; 
      Double_t fGoofieTempF  = 0; 
      fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);
      fGoofieValues->SetBranchAddress("fGoofieTempN",&fGoofieTempN); 
      fGoofieValues->SetBranchAddress("fGoofieTempF",&fGoofieTempF);  
      
      Long64_t nevent = fGoofieValues->GetEntries();
      for (int i = 0; i < nevent; ++i){
	if (fGoofieValues->GetEvent(i)<0){
	  // cout<< "you're done, man  !" << endl;
	  continue;
	} 
	if (i==0) fStartTime = fGoofieTime;
	else if (i == (nevent-1)) fEndTime = fGoofieTime;    
	fTempGradGraph->SetPoint(i,(fGoofieTime - fStartTime),(fGoofieTempF-fGoofieTempN)/25); 
      }
    }
  }

  //____________________________________________________________________________________
  void AliTPCGoofieValues::FillPressGraph()
  {
    /**
       Graph of pressure<br>
       If the graph is already filled, it's only putting a title. 
     */

    fPressGraph->SetTitle("time dpendance of Pressure");     
    if(!fPressGraph->GetN()==0){
      // cout<< " graph already filled !!! " << endl; 
    }
    else{
      Double_t fGoofieTime   = 0; 
      Double_t fGoofiePress  = 0;      
      fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);
      fGoofieValues->SetBranchAddress("fGoofiePress",&fGoofiePress);
     
      Long64_t nevent = fGoofieValues->GetEntries();
      for (int i = 0; i < nevent; ++i){
	if (fGoofieValues->GetEvent(i)<0){
	  // cout<< "you're done, man  !" << endl;
	  continue;
	} 
	if (i==0) fStartTime = fGoofieTime;
	else if (i == (nevent-1)) fEndTime = fGoofieTime;      
	fPressGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofiePress);
      }
    }
  }

  
//________________________________________________________________________________________________
  void AliTPCGoofieValues::FillGainFGraph()
  {

    /**
       return graph of Gain in the far point of Goofie`s detector<br>
       If the graph is already filled, it's only putting a title.
     */

   fGainFGraph->SetTitle("time dpendance of Gain in the far point in natural unit");

   if(!fGainFGraph->GetN()==0){
     // cout<< " graph already filled !!! " << endl; 
   }
   else{
     Double_t fGoofieTime   = 0;    
     Double_t fGoofieAreaF  = 0; 
     fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);
     fGoofieValues->SetBranchAddress("fGoofieAreaF",&fGoofieAreaF);
     Long64_t nevent = fGoofieValues->GetEntries();
     for (int i = 0; i < nevent; ++i){
       if (fGoofieValues->GetEvent(i)<0){
	 // cout<< "you're done, man  !" << endl;
	 continue;
       } 
       if (i==0) fStartTime = fGoofieTime;
       else if (i == (nevent-1)) fEndTime = fGoofieTime;          
       fGainFGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieAreaF);
     }
   }
  }
  
//_____________________________________________________________________________________
  void AliTPCGoofieValues::FillGainNGraph()
  {
    /**
       Graph of Gain in the near point of Goofie`s detector<br>
       If the graph is already filled, it's only putting a title.
     */

   fGainNGraph->SetTitle("time dependance of Gain in the far poin in natural unit");
   if(!fGainNGraph->GetN()==0){
     // cout<< " graph already filled !!! " << endl; 
   }
   else{

     Double_t fGoofieTime   = 0;    
     Double_t fGoofieAreaN  = 0; 
     fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);
     fGoofieValues->SetBranchAddress("fGoofieAreaN",&fGoofieAreaN);
     Long64_t nevent = fGoofieValues->GetEntries();
     for (int i = 0; i < nevent; ++i){
       if (fGoofieValues->GetEvent(i)<0){
	 // cout<< "you're done, man  !" << endl;
	 continue;
       } 
       if (i==0) fStartTime = fGoofieTime;
       else if (i == (nevent-1)) fEndTime = fGoofieTime;          
       fGainNGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieAreaN);
     }
   }
   
  }
//_______________________________________________________________________________________
void AliTPCGoofieValues::FillVdriftGraph()
  {
    /**
       Graph of raw Vdrift<br>
       If the graph is already filled, it's only putting a title.
    */
 
    fVdriftGraph->SetTitle("time dpendance of raw Vdrift");
    if(!fVdriftGraph->GetN()==0){
      // cout<< " graph already filled !!! " << endl; 
    }
    else{  
      Double_t fGoofieTime   = 0;      
      Double_t fGoofieVdrift = 0; 
      fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);     
      fGoofieValues->SetBranchAddress("fGoofieVdrift",&fGoofieVdrift);
      Long64_t nevent = fGoofieValues->GetEntries();
      for (int i = 0; i < nevent; ++i){
	if (fGoofieValues->GetEvent(i)<0){
	  // cout<< "you're done, man  !" << endl;
	  continue;
	} 
	if (i==0) fStartTime = fGoofieTime;
	else if (i == (nevent-1)) fEndTime = fGoofieTime;         
	fVdriftGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieVdrift);   
      }
    }
  }  
//_____________________________________________________________________________________________
void AliTPCGoofieValues::FillVdriftcorGraph()
  {
   /**
       Graph of Vdrift corrected by temperature and pressure<br>
       If the graph is already filled, it's only putting a title. 
   */
 
    fVdriftcorGraph->SetTitle("time dpendance of raw Vdriftcor");   
    if(!fVdriftcorGraph->GetN()==0){
      // cout<< " graph already filled !!! " << endl; 
    }
    else{   
      Double_t fGoofieTime   = 0;      
      Double_t fGoofieVdriftcor = 0; 
      fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);     
      fGoofieValues->SetBranchAddress("fGoofieVdriftcor",&fGoofieVdriftcor);
      Long64_t nevent = fGoofieValues->GetEntries();
      for (int i = 0; i < nevent; ++i){
	if (fGoofieValues->GetEvent(i)<0){
	  // cout<< "you're done, man  !" << endl;
	  continue;
	} 
	if (i==0) fStartTime = fGoofieTime;
	else if (i == (nevent-1)) fEndTime = fGoofieTime;         
	fVdriftcorGraph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieVdriftcor);   
      }
    }
  }
//____________________________________________________________________________________________
void AliTPCGoofieValues::FillCO2Graph()
  {
    /**
       Graph of CO2<br>
       If the graph is already filled, it's only putting a title. 
    */
 
    fCO2Graph->SetTitle("time dependance of CO2");
    if(!fCO2Graph->GetN()==0){
      // cout<< " graph already filled !!! " << endl; 
    }
    else{     
      Double_t fGoofieTime   = 0;    
      fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);        
      Double_t fGoofieCO2  = 0;     
      fGoofieValues->SetBranchAddress("fGoofieCO2",&fGoofieCO2);
      Long64_t nevent = fGoofieValues->GetEntries();
      for (int i = 0; i < nevent; ++i){
	if (fGoofieValues->GetEvent(i)<0){
	  // cout<< "you're done, man  !" << endl;
	  continue;
	} 
	if (i==0) fStartTime = fGoofieTime;
	else if (i == (nevent-1)) fEndTime = fGoofieTime;       
	fCO2Graph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieCO2);  
      }
    }
  }
//_____________________________________________________________________________________________
void AliTPCGoofieValues::FillN2Graph()
  {
    /**
       Graph of N2<br>
       If the graph is already filled, it's only putting a title.
     */

  
    fN2Graph->SetTitle("time dependance of N2");
    if(!fN2Graph->GetN()==0){
      // cout<< " graph already filled !!! " << endl; 
    }
    else{     
      Double_t fGoofieTime   = 0;    
      fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);        
      Double_t fGoofieN2  = 0;     
      fGoofieValues->SetBranchAddress("fGoofieN2",&fGoofieN2);
      Long64_t nevent = fGoofieValues->GetEntries();
      for (int i = 0; i < nevent; ++i){
	if (fGoofieValues->GetEvent(i)<0){
	  // cout<< "you're done, man  !" << endl;
	  continue;
	} 
	if (i==0) fStartTime = fGoofieTime;
	else if (i == (nevent-1)) fEndTime = fGoofieTime;       
	fN2Graph->SetPoint(i,(fGoofieTime - fStartTime),fGoofieN2);  
      }
    } 
  }
  
//_________________________________________________________________________________________
void AliTPCGoofieValues::FillAverageTempSpline()
{
  /**
   TSpline object of average temperature<br>
   If the spline is already filled, it's printing a message.
  */
  if(!fAverageTempSpline)
    fAverageTempSpline = new TSpline3("temperature",fAverageTempGraph);
  //else
    // cout<< " spline already filled !!! " << endl; 
 
}

//_________________________________________________________________________________________
void AliTPCGoofieValues::FillTempGradSpline()
{

  /**
     TSpline object of temperature gradient<br>
     If the spline is already filled, it's printing a message.
  */
  if(!fTempGradSpline)
    fTempGradSpline = new TSpline3("temperature gradient",fTempGradGraph);
  //else
    // cout<< " spline already filled !!! " << endl; 
}

//_________________________________________________________________________________________
void AliTPCGoofieValues::FillPressSpline()
{

  /**
     TSpline object of pressure<br>
     If the spline is already filled, it's printing a message.
  */
  if(!fPressSpline)
    fPressSpline = new TSpline3("pressure",fPressGraph);
  //else
    // cout<< " spline already filled !!! " << endl; 
  
}

//_________________________________________________________________________________________
void AliTPCGoofieValues::FillVdriftSpline()
{

  /**
     TSpline object of drift velocity<br>
     If the spline is already filled, it's printing a message.
  */
  if(!fVdriftSpline)
    fVdriftSpline = new TSpline3("vdrift",fVdriftGraph);
  //else{
    // cout<< " spline already filled !!! " << endl; 
  //}
}

//_________________________________________________________________________________________
void AliTPCGoofieValues::FillVdriftcorSpline()
{


  /**
     TSpline object of drift velocit corrected<br>
     If the spline is already filled, it's printing a message.
  */
  if(!fVdriftcorSpline)
    fVdriftcorSpline = new TSpline3("vdriftcor",fVdriftcorGraph);
  //else
    // cout<< " spline already filled !!! " << endl; 

}

//_________________________________________________________________________________________
void AliTPCGoofieValues::FillGainFSpline()
{

  /**
     return TSpline object of gain in far point<br>
     If the spline is already filled, it's printing a message.
  */
  if(!fGainFSpline)
    fGainFSpline = new TSpline3("gainF",fGainFGraph);
  //else
    // cout<< " spline already filled !!! " << endl; 

}

//_________________________________________________________________________________________
void AliTPCGoofieValues::FillGainNSpline()
{

  /**
     TSpline object of gain in near point<br>
     If the spline is already filled, it's printing a message.
  */
  if(!fGainNSpline)
    fGainNSpline = new TSpline3("gainN",fGainNGraph);
  //else
    // cout<< " spline already filled !!! " << endl; 
 
}

//_________________________________________________________________________________________
void AliTPCGoofieValues::FillCO2Spline()
{

  /**
     TSpline object of average temperature<br>
     If the spline is already filled, it's printing a message.
  */  
  if(!fCO2Spline)
    fCO2Spline = new TSpline3("co2",fCO2Graph);
  //else
    // cout<< " spline already filled !!! " << endl; 
}

//_________________________________________________________________________________________
void AliTPCGoofieValues::FillN2Spline()
{

  /**
     TSpline object of average temperature<br>
     If the spline is already filled, it's printing a message.
  */ 
  if(!fN2Spline)
    fN2Spline = new TSpline3("n2",fN2Graph);
  //else
    // cout<< " spline already filled !!! " << endl; 
}

//__________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalTempGrad(Double_t timeSec)
{
  if(!fTempGradSpline)
    fTempGradSpline = new TSpline3("temperature gradient",fTempGradGraph);
  Double_t a = fTempGradSpline->Eval(timeSec);
  return a;
}

//___________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalAverageTemp(Double_t timeSec)
{
  if(!fAverageTempSpline)
    fAverageTempSpline = new TSpline3("temperature",fAverageTempGraph);
  Double_t a = fAverageTempSpline->Eval(timeSec);
  return a;
}

//___________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalPress(Double_t timeSec)
{
  if(!fPressSpline)
    fPressSpline = new TSpline3("pressure",fPressGraph);
  Double_t a = fPressSpline->Eval(timeSec);
  return a;
}

//___________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalVdrift(Double_t timeSec)
{
  if(!fVdriftSpline)
    fVdriftSpline = new TSpline3("vdrift",fVdriftGraph);
  Double_t a = fVdriftSpline->Eval(timeSec);
  return a;
}

//___________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalVdriftcor(Double_t timeSec)
{
  if(!fVdriftcorSpline)
    fVdriftcorSpline = new TSpline3("vdriftcor",fVdriftcorGraph);
  Double_t a = fVdriftcorSpline->Eval(timeSec);
  return a;
}

//___________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalGainF(Double_t timeSec)
{
  if(!fGainFSpline)
    fGainFSpline = new TSpline3("gainF",fGainFGraph);
  Double_t a = fGainFSpline->Eval(timeSec);
  return a;
}

//___________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalGainN(Double_t timeSec)
{
  if(!fGainNSpline)
    fGainNSpline = new TSpline3("gainN",fGainNGraph);
  Double_t a  = fGainNSpline->Eval(timeSec);
  return a;
}

//___________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalCO2(Double_t timeSec)
{
  if(!fCO2Spline)
    fCO2Spline  = new TSpline3("co2",fCO2Graph);
  Double_t a = fCO2Spline->Eval(timeSec);
  return a;
}

//___________________________________________________________________________________________
Double_t AliTPCGoofieValues::EvalN2(Double_t timeSec)
{
  if(!fN2Spline)
    fN2Spline = new TSpline3("n2",fN2Graph);
  Double_t a = fN2Spline->Eval(timeSec);
  return a;
}

void AliTPCGoofieValues::PrintTree(){
  /** 
      Testing function: it prints all the information stored in the tree
   */

  Long64_t nevent = fGoofieValues->GetEntries();
  // cout<< " number of entries: " << nevent << endl;
  //temporal for reading the branches
  Double_t fGoofieTime   = 0; 
  Double_t fGoofieTempN  = 0; 
  Double_t fGoofieTempF  = 0; 
  Double_t fGoofiePress  = 0;   
  Double_t fGoofieAreaF  = 0; 
  Double_t fGoofieAreaN  = 0; 
  
  Double_t fGoofieVdriftcor  = 0; 
  Double_t fGoofieVdrift     = 0; 
  Double_t fGoofieCO2  = 0; 
  Double_t fGoofieN2   = 0; 

  fGoofieValues->SetBranchAddress("fGoofieTime",&fGoofieTime);
  fGoofieValues->SetBranchAddress("fGoofieTempN",&fGoofieTempN); 
  fGoofieValues->SetBranchAddress("fGoofieTempF",&fGoofieTempF);  
  fGoofieValues->SetBranchAddress("fGoofiePress",&fGoofiePress);
  
  fGoofieValues->SetBranchAddress("fGoofieAreaF",&fGoofieAreaF);
  fGoofieValues->SetBranchAddress("fGoofieAreaN",&fGoofieAreaN);
  
  fGoofieValues->SetBranchAddress("fGoofieVdriftcor",&fGoofieVdriftcor);
  fGoofieValues->SetBranchAddress("fGoofieVdrift",&fGoofieVdrift);
  
  fGoofieValues->SetBranchAddress("fGoofieCO2",&fGoofieCO2);
  fGoofieValues->SetBranchAddress("fGoofieN2",&fGoofieN2); 

  for (int j = 0; j < nevent; ++j){
    if (fGoofieValues->GetEvent(j)<0){
      // cout<< "you're done, man  !" << endl;
      continue;
    } 
//     // cout<< "  reading event: " << j << " " << fGoofieTime << " " 
// 	 << fGoofieTempN << " " << fGoofieTempF << " " 
// 	 << fGoofiePress << " " << fGoofiePress << " " 
// 	 << fGoofieAreaF << " " << fGoofieAreaN << " " 
// 	 << fGoofieVdriftcor << " " << fGoofieVdrift << " " 
// 	 <<  fGoofieCO2 << " " << fGoofieN2 << endl;
  }
}

  
  
  
  
  
  
  
