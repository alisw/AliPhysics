
/////////////////////////////////////////////////////////////////////////////////
//                                                                             
// AliTRDCalibra                                                               
//                                                                             
// This class is for the TRD calibration of the relative gain factor, the drift velocity,
// the time 0 and the pad response function.        
// It can be used for the calibration per chamber but also per group of pads and eventually per pad.
// The user has to choose with the functions SetNz and SetNrphi the precision of the calibration. 
//Begin_Html
//<br>
//<CENTER>
//<TABLE border=1>
//<TR><TD><center>Nz</center></TD><TD><center> 0 </center></TD><TD><center> 1 </center></TD><TD><center> 2 </center></TD><TD><center> 3 </center></TD><TD><center> 4 </center></TD></TR>
//<TR><TD><CENTER>group of row pads per detector</CENTER></TD><TD><CENTER>1</CENTER></TD><TD><CENTER>2</CENTER></TD><TD><CENTER>4</CENTER></TD><TD><CENTER>6(chamb2)<br> 8(others chambers)</CENTER></TD><TD><CENTER>12 (chamb2)<br> 16 (chamb0)</CENTER></TD></TR>
//<TR><TD><CENTER>row pads per group</CENTER></TD><TD><CENTER>12 (chamb2)<br> 16 (chamb0)</CENTER></TD><TD><CENTER>6 (chamb2)<br> 8 (chamb0)</CENTER></TD><TD><CENTER>3 (chamb2)<br> 4 (chamb0)</CENTER></TD><TD><CENTER>2</CENTER></TD><TD><CENTER>1</CENTER></TD></TR>
//<TR><TD><CENTER>~distance [cm]</CENTER></TD><TD><CENTER>106 (chamb2)<br> 130 (chamb0)</CENTER></TD><TD><CENTER>53 (chamb2)<br> 65 (chamb0)</CENTER></TD><TD><CENTER>26.5 (chamb2)<br> 32.5 (chamb0)</CENTER></TD><TD><CENTER>17 (chamb2)<br> 17 (chamb0)</CENTER></TD><TD><CENTER>9 (chamb2)<br> 9 (chamb0)</CENTER></TD></TR>
//<CAPTION>In the z direction</CAPTION>
//</TABLE>
//</CENTER>
//<CENTER>
//<br>
//<TABLE border=1>
//<TR><TD><center>Nrphi</center></TD><TD><center> 0 </center></TD><TD><center> 1 </center></TD><TD><center> 2 </center></TD><TD><center> 3 </center></TD><TD><center> 4 </center></TD><TD><center> 5 </center></TD><TD><center> 6 </center></TD></TR>
//<TR><TD><CENTER>group of col pads per detector</CENTER></TD><TD><CENTER>1</CENTER></TD><TD><CENTER>2</CENTER></TD><TD><CENTER>4</CENTER></TD><TD><CENTER>8</CENTER></TD><TD><CENTER>16</CENTER></TD><TD><center>36</center></TD><TD><center>144</center></TD></TR>
//<TR><TD><CENTER>col pads per group</CENTER></TD><TD><CENTER>144</CENTER></TD><TD><CENTER>72</CENTER></TD><TD><CENTER>36</CENTER></TD><TD><CENTER>18</CENTER></TD><TD><CENTER>9</CENTER></TD><TD><center>4</center></TD><TD><center>1</center></TD></TR>
//<TR><TD><CENTER>~distance [cm]</CENTER></TD><TD><CENTER>113.4</CENTER></TD><TD><CENTER>56.7</CENTER></TD><TD><CENTER>25.3</CENTER></TD><TD><CENTER>14.3</CENTER></TD><TD><CENTER>7.25</CENTER></TD><TD><center>3.2</center></TD><TD><center>0.8</center></TD></TR>
//<CAPTION>In the rphi direction</CAPTION>
//</TABLE>
//</CENTER>
//<br>
//End_Html 
//
// Fill histograms or vectors
//----------------------------
//   
// 2D Histograms (Histo2d) or vectors (Vector2d), then converted in Trees, will be filled
// from RAW DATA in a run or from reconstructed TRD tracks during the offline tracking 
// in the function "FollowBackProlongation" (AliTRDtracker)
// Per default the functions to fill are off.                                   
//
//Begin_Html
//Example of 2D histos for the relative gain (left) and the drift velocity (right) calibration of the sector 13 <br>
//<center>
//<img src="./gif/2dhisto.gif" width="600" height="350"><br>
//</center>
//End_Html    
//
// Fit the histograms to find the coefficients 
//---------------------------------------------
//
// These 2D histograms or vectors (first converted in 1D histos) will be fitted  
// if they have enough entries, otherwise the (default) value of the choosen database 
// will be put. For the relative gain calibration the resulted factors will be globally 
// normalized to the gain factors of the choosen database. It unables to precise     
// previous calibration procedure.
// The function SetDebug enables the user to see:                                     
// _fDebug = 0: nothing, only the values are written in the tree if wanted
// _fDebug = 1: a comparaison of the coefficients found and the default values 
//              in the choosen database.
//              fCoef , histogram of the coefs as function of the calibration group number
//              fDelta , histogram of the relative difference of the coef with the default
//                        value in the database as function of the calibration group number
//              fError , dirstribution of this relative difference
// _fDebug = 2: only the fit of the choosen calibration group fFitVoir (SetFitVoir)
// _fDebug = 3: The coefficients in the choosen detector fDet (SetDet) as function of the
//              pad row and col number
// _fDebug = 4; The coeffcicients in the choosen detector fDet (SetDet) like in the 3 but with
//              also the comparaison histograms of the 1 for this detector
//
//Begin_Html
//Example of fCoef for the relative gain calibration of the sector 13 <br>
//<center>
//<img src="./gif/coef.gif"  width="400" height="460">
//</center><br>
//Example of fDelta (right) and fError (left) for the relative gain calibration of the sector 13 <br>
//<center>
//<img src="./gif/detlaerror.gif"  width="550" height="380"><br>
//</center>
//End_Html 
//
// Author:
//   R. Bailhache (R.Bailhache@gsi.de)
//                            
//////////////////////////////////////////////////////////////////////////////////////


#include "AliTRDCalibra.h"

using namespace std;


ClassImp(AliTRDCalibra)

AliTRDCalibra* AliTRDCalibra::fgInstance = 0;
Bool_t AliTRDCalibra::fgTerminated = kFALSE;

//_____________singleton implementation_________________________________________________
AliTRDCalibra* AliTRDCalibra::Instance()
{
  //
  // Singleton implementation
  //

  if(fgTerminated != kFALSE)
    return 0;

  if (fgInstance == 0)
    fgInstance = new AliTRDCalibra();

  return fgInstance;

}

//______________________________________________________________________________________
void AliTRDCalibra::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class
  //

  fgTerminated = kTRUE;

  if (fgInstance != 0) {
    delete fgInstance;
    fgInstance = 0;
  }

}

//______________________________________________________________________________________
AliTRDCalibra::AliTRDCalibra()
  :TObject()
  ,fMItracking(kFALSE)
  ,fmcmtracking(kFALSE)
  ,fmcmcorrectangle(kFALSE)
  ,fCH2dOn(kFALSE)
  ,fPH2dOn(kFALSE)
  ,fPRF2dOn(kFALSE)
  ,fHisto2d(kFALSE)
  ,fVector2d(kFALSE)
  ,fRelativeScale(0)
  ,fCountRelativeScale(0)
  ,fRelativeScaleAuto(kFALSE)
  ,fThresholddigit(0)
  ,fThresholdClusterPRF1(0.0)
  ,fThresholdClusterPRF2(0.0)
  ,fCenterOfflineCluster(kFALSE)
  ,fTraMaxPad(kFALSE)
  ,fWriteNameCoef(0)
  ,fWriteName(0)
  ,fFitPHOn(kFALSE)
  ,fFitPHPeriode(0)
  ,fBeginFitCharge(0.0)
  ,fRangeFitPRF(0.0)
  ,fMeanChargeOn(kFALSE)
  ,fFitChargeBisOn(kFALSE)
  ,fT0Shift(0.0)
  ,fDebug(0)
  ,fFitVoir(0)
  ,fPRF(0)
  ,fGain(0)
  ,fT0(0)
  ,fVdrift(0)
  ,fVdriftDetector(0)
  ,fVdriftPad(0x0)
  ,fT0Detector(0)
  ,fT0Pad(0x0)
  ,fPRFDetector(0)
  ,fPRFPad(0x0)
  ,fcoefCH(0x0)
  ,fDetectorAliTRDtrack(kFALSE)
  ,fChamberAliTRDtrack(-1)
  ,fDetectorprevioustrack(-1)
  ,fGoodTrack(kTRUE)
  ,famptotal(0x0)
  ,fPHplace(0x0)
  ,fPHvalue(0x0)
  ,fNumberClusters(0)
  ,fprocent(0.0)
  ,fdifference(0)
  ,fNumbertrack(0)
  ,fDeltaPRF(0)
  ,fErrorPRF(0)
  ,fCoefPRFDB(0)
  ,fTimeMax(0)
  ,fSf(0.0)
  ,fScalefitfactor(0.0)
  ,fMinEntries(0)
  ,fEntriesCurrent(0)
  ,fl3P0(0.0)
  ,fl3P2(0.0)
  ,fg3P2(0.0)
  ,fVectorPH(0)
  ,fPlaPH(0)
  ,fNumberBinCharge(0)
  ,fVectorCH(0)
  ,fPlaCH(0)
  ,fVectorFitCH(0)
  ,fNumberBinPRF(0)
  ,fVectorPRF(0)
  ,fPlaPRF(0)
  ,fPH2d(0x0)
  ,fPRF2d(0x0)
  ,fCH2d(0x0)
  ,fRebin(0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 3; i++) {
    fNz[i]    = 0;
    fNrphi[i] = 0;
  }

  for (Int_t k = 0; k < 3; k++) {
    fNtotal[k]    = 0;
    fdetChamb2[k] = 0;
    fdetChamb0[k] = 0;
  }

  // Write
  for (Int_t i = 0; i < 3; i++) {
    fWriteCoef[i] = kFALSE;
    fWrite[i]     = kFALSE;
  }

  // Debug Mode
  for (Int_t k = 0; k < 3; k++) {
    fDet[k] = 0;
  }

  for (Int_t i = 0; i < 2; i++) {
    fPhd[i] = 0.0;
  }
  fPhd[3] = 0.0;

  // Init
  Init();
  
}

//______________________________________________________________________________________
AliTRDCalibra::AliTRDCalibra(const AliTRDCalibra &c)
  :TObject(c)
  ,fMItracking(kFALSE)
  ,fmcmtracking(kFALSE)
  ,fmcmcorrectangle(kFALSE)
  ,fCH2dOn(kFALSE)
  ,fPH2dOn(kFALSE)
  ,fPRF2dOn(kFALSE)
  ,fHisto2d(kFALSE)
  ,fVector2d(kFALSE)
  ,fRelativeScale(0)
  ,fCountRelativeScale(0)
  ,fRelativeScaleAuto(kFALSE)
  ,fThresholddigit(0)
  ,fThresholdClusterPRF1(0.0)
  ,fThresholdClusterPRF2(0.0)
  ,fCenterOfflineCluster(kFALSE)
  ,fTraMaxPad(kFALSE)
  ,fWriteNameCoef(0)
  ,fWriteName(0)
  ,fFitPHOn(kFALSE)
  ,fFitPHPeriode(0)
  ,fBeginFitCharge(0.0)
  ,fRangeFitPRF(0.0)
  ,fMeanChargeOn(kFALSE)
  ,fFitChargeBisOn(kFALSE)
  ,fT0Shift(0.0)
  ,fDebug(0)
  ,fFitVoir(0)
  ,fPRF(0)
  ,fGain(0)
  ,fT0(0)
  ,fVdrift(0)
  ,fVdriftDetector(0)
  ,fVdriftPad(0x0)
  ,fT0Detector(0)
  ,fT0Pad(0x0)
  ,fPRFDetector(0)
  ,fPRFPad(0x0)
  ,fcoefCH(0x0)
  ,fDetectorAliTRDtrack(kFALSE)
  ,fChamberAliTRDtrack(-1)
  ,fDetectorprevioustrack(-1)
  ,fGoodTrack(kTRUE)
  ,famptotal(0x0)
  ,fPHplace(0x0)
  ,fPHvalue(0x0)
  ,fNumberClusters(0)
  ,fprocent(0.0)
  ,fdifference(0)
  ,fNumbertrack(0)
  ,fDeltaPRF(0)
  ,fErrorPRF(0)
  ,fCoefPRFDB(0)
  ,fTimeMax(0)
  ,fSf(0.0)
  ,fScalefitfactor(0.0)
  ,fMinEntries(0)
  ,fEntriesCurrent(0)
  ,fl3P0(0.0)
  ,fl3P2(0.0)
  ,fg3P2(0.0)
  ,fVectorPH(0)
  ,fPlaPH(0)
  ,fNumberBinCharge(0)
  ,fVectorCH(0)
  ,fPlaCH(0)
  ,fVectorFitCH(0)
  ,fNumberBinPRF(0)
  ,fVectorPRF(0)
  ,fPlaPRF(0)
  ,fPH2d(0x0)
  ,fPRF2d(0x0)
  ,fCH2d(0x0)
  ,fRebin(0)
{
  //
  // Copy constructor
  //

}

//____________________________________________________________________________________
AliTRDCalibra::~AliTRDCalibra()
{
  //
  // AliTRDCalibra destructor
  //

  ClearHistos();
  ClearTree();
  ClearFile();

}

//_____________________________________________________________________________
void AliTRDCalibra::Destroy() 
{
  //
  // delete instance 
  //

  if (fgInstance) {
    //fgInstance->Delete();
    delete fgInstance;
    fgInstance = 0x0;
  }
}

//_____________________________________________________________________________
void AliTRDCalibra::ClearHistos() 
{
  //
  // delete the histos
  //

  if (fPH2dOn && fPH2d) {
    delete fPH2d;
    fPH2d = 0x0;
  }
  if (fCH2dOn && fCH2d) {
    delete fCH2d;
    fCH2d = 0x0;
  }
  if (fPRF2dOn && fPRF2d) {
    delete fPRF2d;
    fPRF2d = 0x0;
  }

}

//_____________________________________________________________________________
void AliTRDCalibra::ClearTree() 
{
  //
  // delete the trees
  //
  
  if (fPRF2dOn && fPRF && fWriteCoef[3]) {
    delete fPRF;
    fPRF = 0x0;
  }
  if (fCH2dOn && fGain && fWriteCoef[0]) {
    delete fGain;
    fGain = 0x0;
  }
  if (fPH2dOn && fT0 && fWriteCoef[2]) {
    delete fT0;
    fT0 = 0x0;
  }
  if (fPH2dOn && fVdrift && fWriteCoef[1]) {
    delete fVdrift;
    fVdrift = 0x0;
  }

}

//_____________________________________________________________________________
void AliTRDCalibra::Init() 
{
  //
  // Init some default values
  //

  //How to fill the 2D
  fThresholddigit = 5;
  fThresholdClusterPRF1 = 2.0;
  fThresholdClusterPRF2 = 20.0;
  

  //Store the Info
  fNumberBinCharge = 100;
  fNumberBinPRF = 20;
  
  //write
  fWriteName = "TRD.calibration.root";
  fWriteNameCoef = "TRD.coefficient.root";
  
  //Fit
  fFitPHPeriode = 1;
  fBeginFitCharge = 3.5;
  fRangeFitPRF = 0.5;
  fMinEntries = 800;
  fT0Shift = 0.143397;
  
  
  //Internal variables************************
  
  
  //Fill the 2D histos in the offline tracking
  fDetectorprevioustrack = -1;
  fChamberAliTRDtrack = -1;
  fGoodTrack = kTRUE;

  fprocent = 6.0;
  fdifference = 17;
  fNumberClusters = 18;
  fNumbertrack = 0;
  fNumberusedch[0] = 0;
  fNumberusedch[1] = 0;
  fNumberusedph[0] = 0;
  fNumberusedph[1] = 0;
  
  //For debugging in ClearHistos() 
  
  //variables in the loop
  for(Int_t k = 0; k < 4; k++){
    fChargeCoef[k] = 1.0;
    fVdriftCoef[k] = 1.5;
    fT0Coef[k] = -1.0;
  }
  for(Int_t i = 0; i < 2; i++){
    fPRFCoef[i] = -1.0;
  }
  

  //pad calibration
  for(Int_t i = 0; i < 3; i++){
    frowmin[i] = -1;
    frowmax[i] = -1;
    fcolmax[i] = -1;
    fcolmin[i] = -1;
    fNnz[i] = -1;
    fNnrphi[i] = -1;
    fNfragz[i] = -1;
    fNfragrphi[i] = -1;
    fXbins[i] = -1;
  }
  
  
  //local database to be changed
  fRebin = 1;
  
}

//_______________________Functions fit Online___________________________________________

//____________Functions fit Online CH2d___________________________________________________________

Bool_t AliTRDCalibra::FitCHOnline( TH2I * ch)
{
  //
  // Fit the 1D histos, projections of the 2D ch on the Xaxis, for each calibration group
  // normalized the resulted coefficients (to 1 normally) and write the results in a tree
  //

   
  //Number of Xbins(detectors or groups of pads)*******************************
  TAxis *xch = ch->GetXaxis();
  Int_t Nbins = xch->GetNbins();
  TAxis *yph = ch->GetYaxis();
  Int_t NYbins = yph->GetNbins();
  Double_t lowedge = xch->GetBinLowEdge(1);
  Double_t upedge = xch->GetBinUpEdge(xch->GetNbins());
  if(!InitFit(Nbins,lowedge,upedge,0)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;



  //DB Setting************************************************************************
  //DB and so one****************************************

  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
   
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo(Form("Could not get CommonParam Manager"));
    return kFALSE;
  }
 
  //Init fcountdet and fcount***************************************
  Initfcountdetandfcount(0);
  
  //Beginning of the loop betwwen dect1 and dect2**************************************************************
  for (Int_t idect = fdect1[0]; idect < fdect2[0]; idect++) {

    
    TH1I* projch = (TH1I*)ch->ProjectionY("projch",idect+1,idect+1,(Option_t *)"e");
    projch->SetDirectory(0);


    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,0);
    
    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 0);
    
    
    //Number of entries for this calibration group
    Double_t nentries = 0.0;
    for(Int_t k = 0; k < NYbins; k++){
      nentries += ch->GetBinContent(ch->GetBin(idect+1,k+1));
    }
      
    //Rebin and statistic stuff**********************************
    //Rebin
    if(fRebin > 1) projch = ReBin((TH1I *)projch);
    //This detector has not enough statistics or was off
    if(nentries < fMinEntries) {
      
      //Fill with the default infos***********************************
      NotEnoughStatistic(idect,0);
      //Memory!!!**************************************************
      if(fDebug != 2){
	delete projch;
      }
      continue;
    }
    
    //Statistic of the group fitted************************************
    AliInfo(Form("For the group number %d there are %f stats", idect, nentries));
    statisticmean += nentries;
    numberfit ++; 
    
    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitCH((TH1 *) projch, (Int_t) (idect-fdect1[0]));
    //Method fit bis*********************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist
    if(fFitChargeBisOn) { 
      FitBisCH((TH1 *) projch, (Int_t) (idect-fdect1[0]));
    }
    
    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) FillCoefChargeDB();
    //Fill Infos Fit********************************************
    FillInfosFit(idect,0);
    
   //Memory!!!**************************************************
    if(fDebug != 2){
      delete projch;
    }
    
    
  }//boucle object


  //Normierungcharge************************************************
  if(fDebug != 2) NormierungCharge();
  
  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) ErrorCH();
  
 
  //Plot***********************************************************
  //0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotCH();
   }
  if((fDebug == 4) || (fDebug == 3)){
    PlotCHDB();
   }

  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }
  
  //Write the things!**************************************************
  ConvertVectorFitCHTree();
  if(fWriteCoef[0]){
    WriteFitInfos(0);       
  }
  return kTRUE;
  
}
//____________Functions fit Online CH2d___________________________________________________________

Bool_t AliTRDCalibra::FitCHOnline()
{
  //
  // Reconstruct a 1D histo from the vectorCH for each calibration group, fit the histo, 
  // normalized the resulted coefficients (to 1 normally) and write the results in a tree
  //


   
  //Number of Xbins(detectors or groups of pads)*******************************
  if(!InitFit(0,0,0,0)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;

  //DB Setting************************************************************************
  //DB and so one****************************************
  
  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
 
  
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }

  //Init fcountdet and fcount*************************************************
  Initfcountdetandfcount(0);

  //Beginning of the loop betwwen dect1 and dect2**************************************************************
  for (Int_t idect = fdect1[0]; idect < fdect2[0]; idect++) {
    
    //Search if the group is in the VectorCH
    Int_t place = SearchInVector(idect,0);
       
    //Is in
    TH1F* projch = 0x0;
    TString name("CH");
    name += idect;
    if(place != -1){
      //Variable
      TCTInfo *fCHInfo = new TCTInfo();
      fCHInfo->fentries = new UShort_t[fNumberBinCharge];
      //retrieve
      fCHInfo = fVectorCH[place];
      projch = ConvertVectorCTHisto(fCHInfo,(const char *)name);
      projch->SetDirectory(0);
      delete fCHInfo;
    }

    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,0);
    
    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 0);
    
    //Number of entries
    Double_t nentries = 0.0;
    for(Int_t k = 0; k < fNumberBinCharge; k++){
      nentries += projch->GetBinContent(k+1);
    }
  

    //Rebin and statistic stuff**********************************
    //Rebin
    if((fRebin > 1) && (place != -1)) projch = ReBin((TH1F *) projch);
    //This detector has not enough statistics or was not found in VectorCH
    if((place == -1) || ((place != -1) && (nentries < fMinEntries))) {
      
      //Fill with the default infos**************************************
      NotEnoughStatistic(idect,0);

      //Memory!!!**************************************************
      if(fDebug != 2){
	delete projch;
      }     
      
      continue;
    }

    //Statistic of the histos fitted*********************************************
    AliInfo(Form("For the group number %d there are %f stats", idect, nentries));
    numberfit ++;
    statisticmean += nentries;
    
    
    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitCH((TH1 *) projch, (Int_t) (idect-fdect1[0]));
    
    //Method fit bis*********************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist
    if(fFitChargeBisOn) { 
      FitBisCH((TH1 *) projch, (Int_t) (idect-fdect1[0]));
    }
    
    
    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) FillCoefChargeDB();
    
    //Fill Infos Fit****************************************************
    FillInfosFit(idect,0); 
    
    //Memory!!!**************************************************
    if(fDebug != 2){
      delete projch;
    }
  
    
  }//boucle object


  //Normierungcharge************************************************
  if(fDebug != 2) NormierungCharge();
  
  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) ErrorCH();
  
  
  //Plot***********************************************************
  //0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotCH();
  }
  if((fDebug == 4) || (fDebug == 3)){
    PlotCHDB();
  }
  
  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }
  
  
  //Write the things!**************************************************
  ConvertVectorFitCHTree();
  if(fWriteCoef[0]){
    WriteFitInfos(0);      
  }
  return kTRUE;
  
}
//____________Functions fit Online CH2d___________________________________________________________

Bool_t AliTRDCalibra::FitCHOnline(TTree *tree)
{
  //
  // Look if the calibration group can be found in the tree, if yes take the histo, fit it, 
  // normalized the resulted coefficients (to 1 normally) and write the results in a tree
  //
   
  //Number of Xbins(detectors or groups of pads)*******************************
  if(!InitFit(0,0,0,0)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;
  
  
  //DB Setting************************************************************************
  //DB and so one****************************************
  
  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
  
  
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
  
  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }
  
  //Init fcountdet and fcount*************************************************
  Initfcountdetandfcount(0);
  TH1F* projch = 0x0;
  tree->SetBranchAddress("histo",&projch);
  std::vector<Int_t> vectorplace = ConvertTreeVector(tree);

  //Beginning of the loop betwwen dect1 and dect2**************************************************************
  for (Int_t idect = fdect1[0]; idect < fdect2[0]; idect++) {
    
    //Search if the group is in the VectorCH
    Int_t place = SearchInTreeVector(vectorplace,idect);
    
    //Is in
    if(place != -1){
      //Variable
      tree->GetEntry(place);
    }
    
    
    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,0);
    
    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 0);
    
    //Number of entries
    Double_t nentries = 0.0;
    for(Int_t k = 0; k < fNumberBinCharge; k++){
      nentries += projch->GetBinContent(k+1);
    }
    
    //Rebin and statistic stuff**********************************
    //Rebin
    if((fRebin > 1) && (place != -1)) projch = ReBin((TH1F *) projch);
    //This detector has not enough statistics or was not found in VectorCH
    if((place == -1) || ((place != -1) && (nentries < fMinEntries))) {
      
      //Fill with the default infos**************************************
      NotEnoughStatistic(idect,0);
      
      continue;
    }
    

    //Statistic of the group fitted***************************************
    AliInfo(Form("For the group number %d there are %f stats", idect, nentries));
    numberfit ++;
    statisticmean += nentries;
   

    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitCH((TH1 *) projch, (Int_t) (idect-fdect1[0]));

    //Method fit bis*********************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist
    if(fFitChargeBisOn) { 
      FitBisCH((TH1 *) projch, (Int_t) (idect-fdect1[0]));
    }
  

    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) FillCoefChargeDB();
    
    //Fill Infos Fit****************************************************
    FillInfosFit(idect,0); 
 
    
  }//boucle object


  //Normierungcharge************************************************
  if(fDebug != 2) NormierungCharge();

  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) ErrorCH();
 
 
  //Plot***********************************************************
  //0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotCH();
   }
  if((fDebug == 4) || (fDebug == 3)){
    PlotCHDB();
   }

  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }


  //Write the things!**************************************************
  ConvertVectorFitCHTree();
  if(fWriteCoef[0]){
    WriteFitInfos(0);      
  }
  return kTRUE;
  
}
//________________functions fit Online PH2d____________________________________________________________________

Bool_t AliTRDCalibra::FitPHOnline( TProfile2D *PH)
{
  //
  // Take the 1D profiles (average pulse height), projections of the 2D PH on the Xaxis, for each calibration group
  // Fit or use the slope of the average pulse height to reconstruct the drift velocity
  // write the results in a tree
  // A first calibration of T0 is also made  using the same method (slope method)
  //
 
  //Number of Xbins(detectors or groups of pads)********************************
  TAxis *xph = PH->GetXaxis();
  TAxis *yph = PH->GetYaxis();
  Int_t Nbins = xph->GetNbins();
  Int_t NYbins = yph->GetNbins();
  Double_t lowedge = xph->GetBinLowEdge(1);
  Double_t upedge = xph->GetBinUpEdge(xph->GetNbins());
  if(!InitFit(Nbins,lowedge,upedge,1)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;



  //DB Setting************************************************************************
  //DB and so one****************************************

  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
 
  
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }
  //Init fcountdet and fcount***************************************************
  Initfcountdetandfcount(1);
  
  
  //beginning of the loop***************************************************************
  
  for (Int_t idect = fdect1[1]; idect < fdect2[1]; idect++) {

      
    TH1D* projph = (TH1D *)PH->ProjectionY("projph",idect+1,idect+1,(Option_t *) "e");
    projph->SetDirectory(0); 

    //Number of entries for this calibration group
    Double_t nentries = 0;
    for(Int_t k = 0; k < NYbins; k++){
      nentries += PH->GetBinEntries(PH->GetBin(idect+1,k+1));
    }
  


    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,1);

    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 1);

    
    //Rebin and statistic stuff*********************************
    //This detector has not enough statistics or was off
    if(nentries  < fMinEntries) {
      
      //Fill with the default values**************************************
      NotEnoughStatistic(idect,1);     
      
      //Memory!!!**************************************************
      if(fDebug != 2){
	delete projph;
      }
           
      continue;
    }
    
    //Statistic of the histos fitted*********************************************
    AliInfo(Form("For the group number %d there are %f stats", idect, nentries));
    numberfit ++;
    statisticmean += nentries;
    
    
    //Calcul of "real" coef****************************************
    CalculVdriftCoefMean(fcountdet[1],(Int_t) (idect-fdect1[1]));
    CalculT0CoefMean(fcountdet[1],(Int_t) (idect-fdect1[1]));
    

    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitPente((TH1 *) projph, (Int_t) (idect-fdect1[1]));

    //Method fit bis*********************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist
    if(fFitPHOn) { 
      FitPH((TH1 *) projph, (Int_t) (idect-fdect1[1]));
    }
     
    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) {
      FillCoefVdriftDB();
      FillCoefT0DB();
    }

    //Fill the tree if end of a detector or only the pointer to the branch!!!*******************************************
    FillInfosFit(idect, 1);
    
    
    
    //Memory!!!**************************************************
    if(fDebug != 2){
      delete projph;
    }
  
    
  }//boucle object
  


  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) {
    ErrorPH();
    ErrorT0();
  }
 
 
  //Plot***********************************************************
  //0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotPH();
    PlotT0();
   }
  if((fDebug == 4) || (fDebug == 3)){
    PlotPHDB();
    PlotT0DB();
   }

  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }



  //Write the things!**************************************************
  if(fWriteCoef[1]){
    WriteFitInfos(1);
  }
  return kTRUE;
  
}
//________________functions fit Online PH2d____________________________________________________________________

Bool_t AliTRDCalibra::FitPHOnline()
{
  //
  // Reconstruct the average pulse height from the vectorPH for each calibration group
  // Fit or use the slope of the average pulse height to reconstruct the drift velocity
  // write the results in a tree
  // A first calibration of T0 is also made  using the same method (slope method)
  //

   
  //Number of Xbins(detectors or groups of pads)********************************
  if(!InitFit(0,0,0,1)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;



  //DB Setting************************************************************************
  //DB and so one****************************************

  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
 
  
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }
  //Init fcountdet and fcount********************************************************
  Initfcountdetandfcount(1);

    //beginning of the loop***************************************************************
  
  for (Int_t idect = fdect1[1]; idect < fdect2[1]; idect++) {


    //Search if the group is in the VectorCH
    Int_t place = SearchInVector(idect,1);
    
    //Is in
    TH1F* projph = 0x0;
    TString name("PH");
    name += idect;
    if(place != -1){
      //Variable
      TPInfo *fPHInfo = new TPInfo();
      fPHInfo->fsum = new Float_t[fTimeMax];
      fPHInfo->fsumsquare = new Float_t[fTimeMax];
      fPHInfo->fentries = new UShort_t[fTimeMax]; 
      //retrieve
      fPHInfo = fVectorPH[place];
      projph = CorrectTheError((TGraphErrors *) ConvertVectorPHisto(fPHInfo,(const char *)name));
      projph->SetDirectory(0);
      delete fPHInfo;
    }


    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,1);
    
    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 1);

    //Rebin and statistic stuff*********************************
    //This detector has not enough statistics or was off
    if((place == -1) || ((place != -1) && (fEntriesCurrent < fMinEntries))) {
      
      //Fill with the default values*******************************************
      NotEnoughStatistic(idect,1);

      //Memory!!!**************************************************
      if(fDebug != 2){
	delete projph;
      }
      continue;
    }

    //Statistic of the histos fitted*********************************************
    AliInfo(Form("For the group number %d there are %d stats", idect, fEntriesCurrent));
    numberfit ++;
    statisticmean += fEntriesCurrent;

    //Calcul of "real" coef****************************************
    CalculVdriftCoefMean(fcountdet[1],(Int_t) (idect-fdect1[1]));
    CalculT0CoefMean(fcountdet[1],(Int_t) (idect-fdect1[1]));
    

    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitPente((TH1 *) projph, (Int_t) (idect-fdect1[1]));

    //Method fit bis*********************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist
    if(fFitPHOn) { 
      FitPH((TH1 *) projph, (Int_t) (idect-fdect1[1]));
    }
  
    
    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) {
      FillCoefVdriftDB();
      FillCoefT0DB();
    }
    
    //Fill the tree if end of a detector or only the pointer to the branch!!!*******************************************
    FillInfosFit(idect,1);
    
    
    //Memory!!!**************************************************
    if(fDebug != 2){
      delete projph;
    }
    
    
  }//boucle object
  

  
  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) {
    ErrorPH();
    ErrorT0();
  }
 
  
  //Plot***********************************************************
  //0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotPH();
    PlotT0();
  }
  if((fDebug == 4) || (fDebug == 3)){
    PlotPHDB();
    PlotT0DB();
  }
  
  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }
  

  //Write the things!**************************************************
  if(fWriteCoef[1]){
    WriteFitInfos(1);
  }
  return kTRUE;
  
}

//________________functions fit Online PH2d____________________________________________________________________

Bool_t AliTRDCalibra::FitPHOnline(TTree *tree)
{
  //
  // Look if the calibration group can be found in the tree, if yes take the histo, fit it, 
  // and write the results in a tree
  // A first calibration of T0 is also made  using the same method (slope method)
  //


   
  //Number of Xbins(detectors or groups of pads)********************************
  if(!InitFit(0,0,0,1)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;
  


  //DB Setting************************************************************************
  //DB and so one****************************************

  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
  
  
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
  
  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }

  //Init fcountdet and fcount********************************************************
  Initfcountdetandfcount(1);
  TGraphErrors* projphtree = 0x0;
  tree->SetBranchAddress("histo",&projphtree);
  std::vector<Int_t> vectorplace = ConvertTreeVector(tree);
  
  //beginning of the loop***************************************************************
  
  for (Int_t idect = fdect1[1]; idect < fdect2[1]; idect++) {


    //Search if the group is in the VectorCH
    Int_t place = SearchInTreeVector(vectorplace,idect);
    
    TH1F *projph = 0x0;
    //Is in
    if(place != -1){
      //Variable
      tree->GetEntry(place);
      projph = CorrectTheError(projphtree);
    }
    

    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,1);

    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 1);

    //Rebin and statistic stuff*********************************
    //This detector has not enough statistics or was off
    if((place == -1) || ((place != -1) && (fEntriesCurrent < fMinEntries))) {
      
      //Fill with the default values*******************************************
      NotEnoughStatistic(idect,1);
      
      //Memory!!!**************************************************
      if(fDebug != 2){
	delete projph;
      }
      
      continue;
    }
    
    //Statistic of the histos fitted*********************************************
    AliInfo(Form("For the group number %d there are %d stats", idect, fEntriesCurrent));
    numberfit ++;
    statisticmean += fEntriesCurrent;


    //Calcul of "real" coef****************************************
    CalculVdriftCoefMean(fcountdet[1],(Int_t) (idect-fdect1[1]));
    CalculT0CoefMean(fcountdet[1],(Int_t) (idect-fdect1[1]));
    
    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitPente((TH1 *) projph, (Int_t) (idect-fdect1[1]));
    //Method fit bis*********************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist
    if(fFitPHOn) { 
      FitPH((TH1 *) projph, (Int_t) (idect-fdect1[1]));
    }
  

    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) {
      FillCoefVdriftDB();
      FillCoefT0DB();
    }

    //Fill the tree if end of a detector or only the pointer to the branch!!!*******************************************
    FillInfosFit(idect,1);
    

    //Memory!!!**************************************************
    if(fDebug != 2){
      delete projph;
    }

       
  }//boucle object
  


  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) {
    ErrorPH();
    ErrorT0();
  }
 
 
  //Plot***********************************************************
  //0 no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotPH();
    PlotT0();
  }
  if((fDebug == 4) || (fDebug == 3)){
    PlotPHDB();
    PlotT0DB();
  }

  
  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }
  


  //Write the things!**************************************************
  if(fWriteCoef[1]){
    WriteFitInfos(1);
  }
  return kTRUE;
  
}
//________________Functions fit Online PRF2d____________________________________________________________________

Bool_t AliTRDCalibra::FitPRFOnline( TProfile2D *Prf)
{
  //
  //Take the 1D profiles (pad response function), projections of the 2D PRF on the Xaxis, for each calibration group
  //Fit with a gaussian to reconstruct the sigma of the pad response function
  //write the results in a tree
  //

  
  //Number of Xbins(detectors or groups of pads)*********************************
  TAxis *xprf = Prf->GetXaxis();
  TAxis *yprf = Prf->GetYaxis();
  Int_t NYbins = yprf->GetNbins();
  Int_t Nbins = xprf->GetNbins();
  Double_t lowedge = xprf->GetBinLowEdge(1);
  Double_t upedge = xprf->GetBinUpEdge(xprf->GetNbins());
  if(!InitFit(Nbins,lowedge,upedge,2)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;
 
 

  //DB and so one****************************************
  
  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
  
  
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }
  //Init fcountdet and fcount*****************************************************
  Initfcountdetandfcount(2);
  
  //beginning of the loop***************************************************************
  
  for (Int_t idect = fdect1[2]; idect < fdect2[2]; idect++) {

   
    TH1D* projprf = (TH1D*)Prf->ProjectionY("projprf",idect+1,idect+1,(Option_t *) "e");
    projprf->SetDirectory(0);
    
    //Number of entries for this calibration group
    Double_t nentries = 0;
    for(Int_t k = 0; k < NYbins; k++){
      nentries += Prf->GetBinEntries(Prf->GetBin(idect+1,k+1));
    }
    
    
    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,2);
    
    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 2);
    
    
    //Rebin and statistic stuff*********************************
    //This detector has not enough statistics or was off
    if(nentries < fMinEntries) {
      
      //Fill with the default values******************************
      NotEnoughStatistic(idect,2);
      
      //Memory!
      if(fDebug != 2){
	delete projprf;
      }
      
      continue;
    }
    
    
    //Statistic of the histos fitted*********************************************
    AliInfo(Form("For the group number %d there are %f stats", idect, nentries));
    numberfit ++;
    statisticmean += nentries;

    
    //Calcul of "real" coef****************************************
    if((fDebug == 1) || (fDebug == 4)){
      CalculPRFCoefMean(fcountdet[2],(Int_t) (idect-fdect1[2]));
    }
    
    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitPRF((TH1 *) projprf, (Int_t) (idect-fdect1[2]));
    
    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) FillCoefPRFDB();

    //Fill the tree if end of a detector or only the pointer to the branch!!!*******************************************
    FillInfosFit(idect,2);
    
    //Memory!!!**************************************************
    if(fDebug != 2){
      delete projprf;
    }
    
  }//boucle object

  
  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) ErrorPRF();
 
 
  //Plot***********************************************************
  //no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotPRF();
   }
  if((fDebug == 4) || (fDebug == 3)){
    PlotPRFDB();
   }

  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }


  //Write the things!**************************************************
  if(fWriteCoef[2]){
    WriteFitInfos(2);
  }
  return kTRUE;
  
}
//________________Functions fit Online PRF2d____________________________________________________________________

Bool_t AliTRDCalibra::FitPRFOnline(TTree *tree)
{
  //
  //Look if the calibration group can be found in the tree, if yes take the histo, fit it, 
  //and write the results in a tree
  //
  
  //Number of Xbins(detectors or groups of pads)*********************************
  if(!InitFit(0,0,0,2)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;
  

  //DB and so one****************************************

  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
   
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }

  //Init fcountdet and fcount*********************************************************
  Initfcountdetandfcount(2);
  TGraphErrors* projprftree = 0x0;
  tree->SetBranchAddress("histo",&projprftree);
  std::vector<Int_t> vectorplace = ConvertTreeVector(tree);

  //beginning of the loop***************************************************************
  
  for (Int_t idect = fdect1[2]; idect < fdect2[2]; idect++) {

    //Search if the group is in the VectorCH
    Int_t place = SearchInTreeVector(vectorplace,idect);
    
    
    //Is in   
    TH1F *projprf = 0x0;
    if(place != -1){
      //Variable
      tree->GetEntry(place);
      projprf = CorrectTheError(projprftree);
    }

    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,2);

    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 2);


    //Rebin and statistic stuff*********************************
    //This detector has not enough statistics or was off
    if((place == -1) ||((place != -1) &&(fEntriesCurrent < fMinEntries))) {
      
      //Fill with the default values*******************************
      NotEnoughStatistic(idect,2);
      
      //Memory!!!**************************************************
      if(fDebug != 2){
	delete projprf;
      }
      
      continue;
    }

    //Statistic of the histos fitted*********************************************
    AliInfo(Form("For the group number %d there are %d stats", idect, fEntriesCurrent));
    numberfit ++;
    statisticmean += fEntriesCurrent;
        
    //Calcul of "real" coef****************************************
    if((fDebug == 1) || (fDebug == 4)){
      CalculPRFCoefMean(fcountdet[2],(Int_t) (idect-fdect1[2]));
    }
    
    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitPRF((TH1 *) projprf, (Int_t) (idect-fdect1[2]));
    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) FillCoefPRFDB();
    //Fill the tree if end of a detector or only the pointer to the branch!!!*******************************************
    FillInfosFit(idect,2);
    

    //Memory!!!**************************************************
    if(fDebug != 2){
      delete projprf;
    }                
    
  }//boucle object


  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) ErrorPRF();
  
  
  //Plot***********************************************************
  //no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotPRF();
  }
  if((fDebug == 4) || (fDebug == 3)){
    PlotPRFDB();
  }
  
  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }

  
  //Write the things!**************************************************
  if(fWriteCoef[2]){
    WriteFitInfos(2);
  }
  return kTRUE;
  
}
//________________Functions fit Online PRF2d____________________________________________________________________

Bool_t AliTRDCalibra::FitPRFOnline()
{
  //
  //Reconstruct the 1D histo (pad response function) from the vectorPRD for each calibration group
  //Fit with a gaussian to reconstruct the sigma of the pad response function
  //write the results in a tree
  //
  

  
  //Number of Xbins(detectors or groups of pads)*********************************
  if(!InitFit(0,0,0,2)) return kFALSE;
  Double_t statisticmean = 0.0;
  Int_t numberfit = 0;


  //DB and so one****************************************

  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
 
  
  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }
  //Init fcountdet and fcount*********************************************************
  Initfcountdetandfcount(2);

    //beginning of the loop***************************************************************
  
  for (Int_t idect = fdect1[2]; idect < fdect2[2]; idect++) {

    //Search if the group is in the VectorCH
    Int_t place = SearchInVector(idect,2);
    
    
    //Is in
    TH1F* projprf = 0x0;
    TString name("PRF");
    name += idect;
    if(place != -1){
      //Variable
      TPInfo *fPRFInfo = new TPInfo();
      fPRFInfo->fsum = new Float_t[fNumberBinPRF];
      fPRFInfo->fsumsquare = new Float_t[fNumberBinPRF];
      fPRFInfo->fentries = new UShort_t[fNumberBinPRF]; 
      //retrieve
      fPRFInfo = fVectorPRF[place];
      projprf = CorrectTheError((TGraphErrors *) ConvertVectorPHisto(fPRFInfo, (const char *)name));
      projprf->SetDirectory(0);
      delete fPRFInfo;
    }

    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
    Updatefcountdetandfcount(idect,2);

    //Reconstruction of the row and pad group: rowmin, row max..**********************
    Reconstructfitrowminrowmax(idect, 2);


    //Rebin and statistic stuff*********************************
    //This detector has not enough statistics or was off
    if((place == -1) ||((place != -1) &&(fEntriesCurrent < fMinEntries))) {

      //Fill with the default values*******************************
      NotEnoughStatistic(idect,2);

      //Memory************************
      if(fDebug != 2){
	delete projprf;
      }
     

      continue;
    }


    //Statistic of the histos fitted*********************************************
    AliInfo(Form("For the group number %d there are %d stats", idect, fEntriesCurrent));
    numberfit ++;
    statisticmean += fEntriesCurrent;


    //Calcul of "real" coef****************************************
    if((fDebug == 1) || (fDebug == 4)){
      CalculPRFCoefMean(fcountdet[2],(Int_t) (idect-fdect1[2]));
    }
    //Method Mean and fit*****************************************
    //idect is egal for fDebug = 0 and 2, only to fill the hist 
    FitPRF((TH1 *) projprf, (Int_t) (idect-fdect1[2]));
    //Visualise the detector for fDebug 3 or 4*******************************
    //here is the reconstruction of the pad and row group is used!
    if(fDebug >= 3) FillCoefPRFDB();
    //Fill the tree if end of a detector or only the pointer to the branch!!!*******************************************
    FillInfosFit(idect,2);
                    
    //Memory!!!**************************************************
    if(fDebug != 2){
      delete projprf;
    }
    
  }//boucle object



  //Error*******************************************************
  if((fDebug == 1) || (fDebug == 4)) ErrorPRF();
 
 
  //Plot***********************************************************
  //no plot, 1 and 4 error plot, 3 and 4 DB plot
  if((fDebug == 1) || (fDebug == 4)){
    PlotPRF();
   }
  if((fDebug == 4) || (fDebug == 3)){
    PlotPRFDB();
   }

  //Mean Statistic********************************************************
  if(numberfit > 0) {
    AliInfo(Form("There is a mean statistic of: %d", (Int_t)statisticmean/numberfit));
  }
  else {
    AliInfo("There is no fit!");
  }


  //Write the things!**************************************************
  if(fWriteCoef[2]){
    WriteFitInfos(2);
  }
  return kTRUE;
  
}


//____________Functions for initialising the AliTRDCalibra in the code_________________________________________________________________________-

Bool_t AliTRDCalibra::Init2Dhistos()
{
//
//For the offline tracking
//This function will be called in the function AliReconstruction::Run() 
//Init the calibration mode (Nz, Nrphi), the 2D histograms if fHisto2d = kTRUE, 
//



  //DB Setting ************************************************************************************** 
  //Get cal
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

 

  // Some parameters
  fTimeMax     = cal->GetNumberOfTimeBins();
  fSf =  cal->GetSamplingFrequency();
  if(fRelativeScaleAuto) fRelativeScale = 0;
  else fRelativeScale = 20;




  //Create the 2D histos corresponding to the pad groupCalibration mode***************************************************************
  if(fCH2dOn) {
    AliInfo(Form("The pad calibration mode for the relative gain calibration: Nz %d, and Nrphi %d", fNz[0], fNrphi[0]));
    
    //Calcul the number of Xbins
    fNtotal[0] = 0;
    ModePadCalibration(2,0);
    ModePadFragmentation(0,2,0,0);
    fdetChamb2[0] = fNfragz[0]*fNfragrphi[0];
    if(fDebug == 4) {
      AliInfo(Form("For the chamber 2: %d", fdetChamb2[0]));
    }
    fNtotal[0] += 6*18*fdetChamb2[0];
    ModePadCalibration(0,0);
    ModePadFragmentation(0,0,0,0);
    fdetChamb0[0] = fNfragz[0]*fNfragrphi[0];
    if(fDebug == 4) {
      AliInfo(Form("For the other chamber 0: %d", fdetChamb0[0]));
    }
    fNtotal[0] += 6*4*18*fdetChamb0[0];
    AliInfo(Form("Total number of Xbins: %d", fNtotal[0]));


    //Create the 2D histo
    if(fHisto2d) CreateCH2d(fNtotal[0]);

    //Variable
    famptotal = new Float_t[TMath::Max(fdetChamb2[0],fdetChamb0[0])];
    for(Int_t k = 0; k < TMath::Max(fdetChamb2[0],fdetChamb0[0]); k++){
      famptotal[k] = 0.0;
    } 


  }

  if(fPH2dOn) {
    AliInfo(Form("The pad calibration mode for the drift velocity calibration: Nz %d, and Nrphi %d", fNz[1], fNrphi[1]));
    
    //Calcul the number of Xbins
    fNtotal[1] = 0;
    ModePadCalibration(2,1);
    ModePadFragmentation(0,2,0,1);
    fdetChamb2[1] = fNfragz[1]*fNfragrphi[1];
    if(fDebug == 4) {
      AliInfo(Form("For the chamber 2: %d", fdetChamb2[1]));
    }
    fNtotal[1] += 6*18*fdetChamb2[1];
    ModePadCalibration(0,1);
    ModePadFragmentation(0,0,0,1);
    fdetChamb0[1] = fNfragz[1]*fNfragrphi[1];
    if(fDebug == 4) {
      AliInfo(Form("For the chamber 0: %d", fdetChamb0[1]));
    }
    fNtotal[1] += 6*4*18*fdetChamb0[1];
    AliInfo(Form("Total number of Xbins: %d", fNtotal[1]));


    //Create the 2D histo
    if(fHisto2d) CreatePH2d(fNtotal[1]);

   
    //Variable
    fPHplace = new Short_t[fTimeMax];
    for(Int_t k = 0; k < fTimeMax; k++){
      fPHplace[k] = -1;
    } 
    fPHvalue = new Float_t[fTimeMax];
    for(Int_t k = 0; k < fTimeMax; k++){
      fPHvalue[k] = -1.0;
    } 


   


  }

  if(fPRF2dOn) {
    AliInfo(Form("The pad calibration mode for the PRF calibration: Nz %d, and Nrphi %d", fNz[2], fNrphi[2]));
    
    //Calcul the number of Xbins
    fNtotal[2] = 0;
    ModePadCalibration(2,2);
    ModePadFragmentation(0,2,0,2);
    fdetChamb2[2] = fNfragz[2]*fNfragrphi[2];
    if(fDebug == 4) {
      AliInfo(Form("For the chamber 2: %d", fdetChamb2[2]));
    }
    fNtotal[2] += 6*18*fdetChamb2[2];
    ModePadCalibration(0,2);
    ModePadFragmentation(0,0,0,2);
    fdetChamb0[2] = fNfragz[2]*fNfragrphi[2];
    if(fDebug == 4) {
      AliInfo(Form("For the chamber 0: %d", fdetChamb0[2]));
    }
    fNtotal[2] += 6*4*18*fdetChamb0[2];
    AliInfo(Form("Total number of Xbins: %d", fNtotal[2]));


    //Create the 2D histo
    if(fHisto2d) CreatePRF2d(fNtotal[2]);
 
  
  }

  return kTRUE;

}


//__________Functions for filling the histos in the code_____________________________________________________________________________

//_______Offine tracking in the AliTRDtracker_________________________________________________
Bool_t AliTRDCalibra::Resettrack()
{
  //
  //For the offline tracking
  //This function will be called in the function AliTRDtracker::FollowBackPropagation() at the beginning 
  //Reset the parameter to know we have a new TRD track
  //
  
  fDetectorAliTRDtrack = kFALSE;
  return kTRUE;
}


//______________Offline tracking in the AliTRDtracker_________________________________________________________________________________

Bool_t AliTRDCalibra::UpdateHistograms(AliTRDcluster *cl, AliTRDtrack *t)
{
  //
  //For the offline tracking
  //This function will be called in the function AliTRDtracker::FollowBackPropagation()
  //in the loop over the clusters of TRD tracks 
  //Fill the 2D histos or the vectors with the info of the clusters at the end of a detectors
  //if the track is "good"
  //

  
  // Get the parameter object****************************************************
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam");
    return kFALSE;
  }

  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
 
  //Localisation of the detector**************************************************
  Int_t detector = cl->GetDetector();
  Int_t chamber = GetChamber(detector);
  Int_t plane = GetPlane(detector);
  

  //Fill the infos for the previous clusters if not the same detector anymore or if not the same track***********
  if(((detector != fDetectorprevioustrack) || (!fDetectorAliTRDtrack)) && (fDetectorprevioustrack != -1)) {
    
    fNumbertrack++;   
    
    //If the same track, then look if the previous detector is in the same plane, if yes: not a good track
    if(fDetectorAliTRDtrack && (GetPlane(detector) <= GetPlane(fDetectorprevioustrack))) {
      fGoodTrack = kFALSE;
    }

    //Fill only if the track doesn't touch a masked pad or doesn't appear in the middle (fGoodTrack)***************************
    if(fGoodTrack){
     

      //gain calibration************************************
      if(fCH2dOn){
	FillTheInfoOfTheTrackCH();
      }//if CH2dOn


      //PH calibration************************************
      if(fPH2dOn){
	FillTheInfoOfTheTrackPH();    
      }//if PH2dOn
    
    }//if a good track
    
    Resetfvariables();
   

  }//fill at the end the charge
  
  
  //Calcul the position of the detector*********************************************************
  if(detector != fDetectorprevioustrack) {
    LocalisationdetectorXbins(detector);
  }

  //Reset the good track for the PRF***********************************************************8
  Bool_t good = kTRUE;
  
  
  //Localisation of the cluster**************************************************************
  Double_t pos[3] = {0.0,0.0,0.0};
  pos[0] = cl->GetX();
  pos[1] = cl->GetY();
  pos[2] = cl->GetZ();
  Int_t time = cl->GetLocalTimeBin();
  
  //Reset the detector ***********************************************************************
  fDetectorprevioustrack = detector;
  fDetectorAliTRDtrack = kTRUE;
 
  
  //Position of the cluster***********************************************************************
  AliTRDpadPlane *Padplane = parCom->GetPadPlane(plane,chamber);
  Int_t row = Padplane->GetPadRowNumber(pos[2]);
  Double_t OffsetZ = Padplane->GetPadRowOffset(row,pos[2]);
  Double_t OffsetTilt = Padplane->GetTiltOffset(OffsetZ);
  Int_t col = Padplane->GetPadColNumber(pos[1]+OffsetTilt,OffsetZ);
 
  
  //See if we are not near a masked pad******************************************************************
  if(!IsPadOn(detector,col,row)){
    good = kFALSE;
    fGoodTrack = kFALSE;
  }

  if(col > 0){
    if(!IsPadOn(detector,col-1,row)){
      fGoodTrack = kFALSE;
      good = kFALSE;
    }
  }

  if(col < 143){
    if(!IsPadOn(detector,col+1,row)){
      fGoodTrack = kFALSE;
      good = kFALSE;
    }
  }

  //row of the cluster and position in the pad groups***********************************************
  Int_t posr[3] = {0,0,0};
  if((fCH2dOn) && (fNnz[0] != 0)) posr[0] = (Int_t) row/fNnz[0];
  if((fPH2dOn) && (fNnz[1] != 0)) posr[1] = (Int_t) row/fNnz[1];
  if((fPRF2dOn) && (fNnz[2] != 0)) posr[2] = (Int_t) row/fNnz[2];
      
      
  //col of the cluster and position in the pad groups************************************************
  Int_t posc[3] = {0,0,0};
  if((fCH2dOn) && (fNnrphi[0] != 0)) posc[0] = (Int_t) col/fNnrphi[0];
  if((fPH2dOn) && (fNnrphi[1] != 0)) posc[1] = (Int_t) col/fNnrphi[1];
  if((fPRF2dOn) && (fNnrphi[2] != 0)) posc[2] = (Int_t) col/fNnrphi[2];
 
      
  //Charge in the cluster***************************************************************************
  //For the moment take the abs
  Float_t Q = TMath::Abs(cl->GetQ());
  Short_t* signals = cl->GetSignals();
 


  //Correction due to the track angle*************************************************************************
  Float_t correction = 1.0;
  Float_t normalisation = 6.67;
  if((Q >0) && (t->GetNdedx() > 0)){
    correction = t->GetClusterdQdl((t->GetNdedx() - 1))/(Q*normalisation);
  }
 
  //Fill the famptotal with the charge*****************************************************************
  if(fCH2dOn) {
    if(!fTraMaxPad){ 
      famptotal[(Int_t) (posc[0]*fNfragz[0]+posr[0])] += Q*correction;
    }//fTraMaxPad
    else{
      famptotal[(Int_t) (posc[0]*fNfragz[0]+posr[0])] += ((Float_t) signals[3])*correction;
    }//no fTraMaxPad
  }//CH2dOn
  

  //Fill the fPHplace and value *********************************************************************
  if(fPH2dOn) {
    fPHplace[time] = posc[1]*fNfragz[1]+posr[1];
    if(!fTraMaxPad) {
      fPHvalue[time] = Q*correction;
    }//fTraMaxPad
    
    else fPHvalue[time] = ((Float_t) signals[3])*correction;
    
  }//PH2dOn



  //Fill direct the PRF*******************************************************************************  
  if((fPRF2dOn) && (good)){
    Float_t yminus = 0.0;
    Float_t xcenter = 0.0;
    Float_t ycenter = 0.0;
    Float_t ymax = 0.0;
    Bool_t echec = kFALSE;
    
    
    if((cl->From3pad()) && (!cl->IsUsed())){ 
         
      //center 3 balanced
      if((((Float_t)signals[3]) > fThresholdClusterPRF2) && (((Float_t)signals[2]) > fThresholdClusterPRF2) && (((Float_t)signals[4]) > fThresholdClusterPRF2) && (((Float_t)signals[1]) < fThresholdClusterPRF1) && (((Float_t)signals[5]) < fThresholdClusterPRF1) && (Q*correction > 130.0)){
	//col correspond to signals[3]
	if(fCenterOfflineCluster) xcenter = cl->GetCenter();
	else{
	  //security of the denomiateur is 0
	  if((((Float_t)(((Float_t)signals[3])*((Float_t)signals[3])))/((Float_t)(((Float_t)signals[2])*((Float_t)signals[4])))) != 1.0){
	    xcenter = 0.5*(TMath::Log((Float_t)(((Float_t)signals[4])/((Float_t)signals[2]))))/(TMath::Log(((Float_t)(((Float_t)signals[3])*((Float_t)signals[3])))/((Float_t)(((Float_t)signals[2])*((Float_t)signals[4])))));
	  }
	  else xcenter = -100.0;
	}
	if((xcenter > -0.5) && (xcenter < 0.5)) {
	  ycenter = (Float_t)(((Float_t)signals[3])/(((Float_t)signals[2])+((Float_t)signals[3])+(((Float_t)signals[4]))));
	  yminus = (Float_t)(((Float_t)signals[2])/(((Float_t)signals[2])+((Float_t)signals[3])+(((Float_t)signals[4]))));
	  ymax = (Float_t)(((Float_t)signals[4])/(((Float_t)signals[2])+((Float_t)signals[3])+(((Float_t)signals[4]))));
	  if((ycenter > 0.485) && (TMath::Abs(((Float_t)signals[2])+((Float_t)signals[3])+(((Float_t)signals[4])) - Q) < 8.0)) echec = kTRUE;
	}
      }
      
      //fill only if it is in the drift region!
      if((((Float_t)(((Float_t)time)/fSf)) > 0.3) && (echec)) {
	if(fHisto2d) {
	  fPRF2d->Fill((fXbins[2]+posc[2]*fNfragz[2]+posr[2]+0.5),xcenter,ycenter);
	  if(xcenter < 0.0) fPRF2d->Fill((fXbins[2]+posc[2]*fNfragz[2]+posr[2]+0.5),-(xcenter+1.0),yminus);
	  if(xcenter > 0.0) fPRF2d->Fill((fXbins[2]+posc[2]*fNfragz[2]+posr[2]+0.5),1.0-xcenter,ymax);
	}
	if(fVector2d) {
	  UpdateVectorPRF(fXbins[2]+posc[2]*fNfragz[2]+posr[2],xcenter,ycenter);
	  if(xcenter < 0.0) UpdateVectorPRF(fXbins[2]+posc[2]*fNfragz[2]+posr[2],-(xcenter+1.0),yminus);
	  if(xcenter > 0.0) UpdateVectorPRF(fXbins[2]+posc[2]*fNfragz[2]+posr[2],1.0-xcenter,ymax);
	}
	
      }//If in the drift region
    }//cluster isole
  }//PRF2dOn	
  
  return kTRUE;
  
}

//___________________________________Online trackling in AliTRDtrigger________________________________________________________________________________________________________________________________________

Bool_t AliTRDCalibra::UpdateHistogramcm(AliTRDmcmTracklet *fTrk)
{
  //
  //For the tracking
  //This function will be called in the function AliTRDtrigger::TestTracklet before applying the pt cut on the tracklets 
  //Fill the infos for the tracklets fTrkTest if the tracklets is "good"
  //
  
  //Localisation of the Xbins involved********************************************
  Int_t idect = fTrk->GetDetector();
  LocalisationdetectorXbins(idect);

  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
   
  //Reset************************************************************************8
  Resetfvariables();
  
    
  //row of the tracklet and position in the pad groups******************************
  Int_t row = fTrk->GetRow();
  Int_t posr[3] = {0,0,0};
  if((fCH2dOn) && (fNnz[0] != 0)) posr[0] = (Int_t) row/fNnz[0];
  if((fPH2dOn) && (fNnz[1] != 0)) posr[1] = (Int_t) row/fNnz[1];
  if((fPRF2dOn) && (fNnz[2] != 0)) posr[2] = (Int_t) row/fNnz[2];
     
  //Eventuelle correction due to track angle in z direction**********************************
  Float_t correction = 1.0;
  if(fmcmcorrectangle){
    Float_t z = fTrk->GetRowz();
    Float_t r = fTrk->GetTime0();
    correction = r/TMath::Sqrt((r*r+z*z));
  }
  
    
  //boucle sur les clusters*********************************************************
  //Condition on number of cluster: don't come from the middle of the detector******************************************************
  if(fTrk->GetNclusters() >= fNumberClusters){

    for (Int_t icl = 0; icl < fTrk->GetNclusters(); icl++) {
      Float_t amp[3] = {0.0,0.0,0.0};
      Int_t time   = fTrk->GetClusterTime(icl);
      Int_t col = fTrk->GetClusterCol(icl);
            
      amp[0] = fTrk->GetClusterADC(icl)[0]*correction;
      amp[1] = fTrk->GetClusterADC(icl)[1]*correction;
      amp[2] = fTrk->GetClusterADC(icl)[2]*correction;
           
      if (amp[0] < 0.0 || amp[1] < 0.0 || amp[2] < 0.0) continue;
      
      //col of cluster and position in the pad groups**********************************
      Int_t posc[3] = {0,0,0};
      if((fCH2dOn) && (fNnrphi[0] != 0)) posc[0] = (Int_t) col/fNnrphi[0];
      if((fPH2dOn) && (fNnrphi[1] != 0)) posc[1] = (Int_t) col/fNnrphi[1];
      if((fPRF2dOn) && (fNnrphi[2] != 0)) posc[2] = (Int_t) col/fNnrphi[2];
      
           
      //See if we are not near a masked pad**********************************************
      Bool_t good = kTRUE;
      if(!IsPadOn(idect,col,row)){
	
	fGoodTrack = kFALSE;
	good = kFALSE;
	
      }
      
      if(col > 0){
	if(!IsPadOn(idect,col-1,row)){
	  fGoodTrack = kFALSE;
	  good = kFALSE;
	}
      }
      
      if(col < 143){
	if(!IsPadOn(idect,col+1,row)){
	  fGoodTrack = kFALSE;
	  good = kFALSE;
	}
      }
           
      
      // Total spectrum********************************************************************************
      if(fPH2dOn)  fPHplace[time] = posc[1]*fNfragz[1]+posr[1];
      if(!fTraMaxPad){
	
	if(fCH2dOn) famptotal[(Int_t) (posc[0]*fNfragz[0]+posr[0])] += (Float_t) (amp[0]+amp[1]+amp[2]);

	if(fPH2dOn) fPHvalue[time] = (Float_t)(amp[0]+amp[1]+amp[2]);
      
      }
      
      else{
	
	if(fCH2dOn) famptotal[(Int_t) (posc[0]*fNfragz[0]+posr[0])] += (Float_t) amp[1];

	if(fPH2dOn) fPHvalue[time] = amp[1];
       
      }
            
      //Fill PRF direct**********************************************************************************8 
      if(fPRF2dOn && good){
	if((amp[0] > fThresholdClusterPRF2) && (amp[1] > fThresholdClusterPRF2) && (amp[2] > fThresholdClusterPRF2) && ((amp[0]+amp[1]+amp[2]) > 130.0)){
	//security of the denomiateur is 0
	  if((((Float_t)(((Float_t)amp[1])*((Float_t)amp[1])))/((Float_t)(((Float_t)amp[0])*((Float_t)amp[2])))) != 1.0){
	    Float_t xcenter = 0.5*(TMath::Log(amp[2]/amp[0]))/(TMath::Log((amp[1]*amp[1])/(amp[0]*amp[2])));
	    Float_t ycenter = amp[1]/(amp[0]+amp[1]+amp[2]);
	    if((xcenter > -0.5) && (xcenter < 0.5) && (ycenter > 0.485)) {
	      Float_t yminus = amp[0]/(amp[0]+amp[1]+amp[2]);
	      Float_t ymax = amp[2]/(amp[0]+amp[1]+amp[2]);
	      
	      //fill only if it is in the drift region!
	      if(((Float_t)time/fSf) > 0.3) {
		if(fHisto2d) 
		  {
		    fPRF2d->Fill((fXbins[2]+posc[2]*fNfragz[2]+posr[2]+0.5),xcenter,ycenter);
		    if(xcenter < 0.0) fPRF2d->Fill((fXbins[2]+posc[2]*fNfragz[2]+posr[2]+0.5),-(xcenter+1.0),yminus);
		    if(xcenter > 0.0) fPRF2d->Fill((fXbins[2]+posc[2]*fNfragz[2]+posr[2]+0.5),(1.0-xcenter),ymax);
		  }
		if(fVector2d) {
		  UpdateVectorPRF((fXbins[2]+posc[2]*fNfragz[2]+posr[2]),xcenter,ycenter);
		  if(xcenter < 0.0) UpdateVectorPRF(fXbins[2]+posc[2]*fNfragz[2]+posr[2],-(xcenter+1.0),yminus);
		  if(xcenter > 0.0) UpdateVectorPRF(fXbins[2]+posc[2]*fNfragz[2]+posr[2],(1.0-xcenter),ymax);
		}
	      } 
	    }
	  }
	}
      }
      
    }//boucle clusters
    
    
    //Fill the charge************************************************************************************
    if(fCH2dOn && fGoodTrack){
      FillTheInfoOfTheTrackCH();
    }//if CH2don

    //PH calibration************************************
      if(fPH2dOn && fGoodTrack){
	FillTheInfoOfTheTrackPH();	
      }//if PH2dOn
        
  }//condition on number of clusters
  return kTRUE;
  
}

//______________Functions for seeing if the pad is really okey__________________________________________________________________________________________

//___________________________________________________________________________
Bool_t AliTRDCalibra::IsPadOn( Int_t detector, Int_t col, Int_t row)
{
  //
  //look in the choosen database if the pad is On. If no the track will be "not good"
  //
  

  // Get the parameter object
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
  
  Int_t npads = 18;
  Int_t colmcm = (Int_t) col/npads;
  
  if(!cal->IsSuperModuleInstalled(GetSector(detector)) || !cal->IsChamberInstalled(detector) || cal->IsSuperModuleMasked(GetSector(detector)) || cal->IsChamberMasked(detector) || cal->IsMCMMasked(detector,colmcm,row) || cal->IsPadMasked(detector,col,row)){
    
    return kFALSE;
  }
  
  else return kTRUE;
  
}

//______________Functions for plotting the 2D__________________________________________________________________________________________

//___________________________________________________________________________
void AliTRDCalibra::Plot2d()
{
  //
  // Plot the 2D histos 
  //
 
  if(fCH2dOn) PlotCH2d();
  if(fPH2dOn) PlotPH2d();
  if(fPRF2dOn) PlotPRF2d();

}
//____________________Writing the 2D_________________________________________________________________________________________________

//________________________________________________________________________________
Bool_t AliTRDCalibra::Write2d()
{
  //
  //Write the 2D histograms or the vectors converted in trees in the file "TRD.calibration.root" 
  //
  
  TFile *fout = TFile::Open(fWriteName,"RECREATE");
  //Check if the file could be opened
  if (!fout || !fout->IsOpen()) {
    AliInfo("No File found!");
    return kFALSE;
  }
  AliInfo(Form(" Numbertrack: %d Numberusedch[0]: %d, Numberusedch[1]: %d Numberusedph[0]: %d, Numberusedph[1]: %d"
              ,fNumbertrack,fNumberusedch[0],fNumberusedch[1],fNumberusedph[0],fNumberusedph[1]));
  
  
  TStopwatch stopwatch;
  stopwatch.Start();
  AliInfo("Write2d");
  
  if((fCH2dOn ) && (fWrite[0])) {
    if(fHisto2d) fout->WriteTObject(fCH2d);
    if(fVector2d){
      TTree *treeCH2d = ConvertVectorCTTreeHisto(fVectorCH,fPlaCH,"treeCH2d");
      fout->WriteTObject(treeCH2d);
    }
  }
  if((fPH2dOn ) && (fWrite[1])) {
    if(fHisto2d) fout->WriteTObject(fPH2d);
    if(fVector2d){
      TTree *treePH2d = ConvertVectorPTreeHisto(fVectorPH,fPlaPH,"treePH2d");
      fout->WriteTObject(treePH2d);
    }
  }
  if((fPRF2dOn ) && (fWrite[2])) {
    if(fHisto2d) fout->WriteTObject(fPRF2d);
    if(fVector2d){
      TTree *treePRF2d = ConvertVectorPTreeHisto(fVectorPRF,fPlaPRF,"treePRF2d");
      fout->WriteTObject(treePRF2d);
    }
  }
  
  fout->Close();
  
  AliInfo(Form("Execution time Write2d: R:%.2fs C:%.2fs",
	       stopwatch.RealTime(),stopwatch.CpuTime()));
  
  
  return kTRUE;
  
}
//______________________________________________________________________________________________________
AliTRDCalDet* AliTRDCalibra::CreateDetObjectTree(TTree *tree, Int_t i)
{
  //
  // It Creates the AliTRDCalDet object from the tree of the coefficient for the calibration i (i != 2)
  // It takes the mean value of the coefficients per detector 
  // This object has to be written in the database
  //
  
  //Create the DetObject
  AliTRDCalDet *object = 0x0;
  if(i == 0) object = new AliTRDCalDet("ChamberGainFactor","GainFactor (detector value)");
  if(i == 1) object = new AliTRDCalDet("ChamberVdrift","TRD drift velocities (detector value)");
  else object = new AliTRDCalDet("ChamberT0","T0 (detector value)");
  
  //Read the Tree
  Int_t detector = -1;
  Float_t values[2304];
  tree->SetBranchAddress("detector",&detector);
  if(i == 0) tree->SetBranchAddress("gainPad",values);
  if(i == 1) tree->SetBranchAddress("vdrift",values);
  if(i == 3) tree->SetBranchAddress("t0",values);
  
  //For calculating the mean
  Float_t mean = 0.0;
  Int_t Nto = 0;
  Int_t numberofentries = tree->GetEntries();
  
  if(numberofentries != 540) AliInfo("The tree is not complete");
  
  for (Int_t det= 0; det < numberofentries; ++det){
    tree->GetEntry(det);
    if(GetChamber(detector) == 2) Nto = 1728;
    else Nto = 2304;
    mean = 0.0;
    for(Int_t k = 0; k < Nto; k++){
      mean += TMath::Abs(values[k])/Nto;  
    } 
    object->SetValue(detector,mean);
  }
  return object;
}
//______________________________________________________________________________________________________
TObject* AliTRDCalibra::CreatePadObjectTree(TTree *tree, Int_t i, AliTRDCalDet *detobject)
{
  //
  // It Creates the AliTRDCalPad object from the tree of the coefficient for the calibration i (i != 2)
  // You need first to create the object for the detectors, where the mean value is put.
  // This object has to be written in the database
  //
  
  
  //Create the DetObject
  AliTRDCalPad *object = 0x0;
  if(i == 0) object = new AliTRDCalPad("GainFactor","GainFactor (local variations)");
  if(i == 1) object = new AliTRDCalPad("LocalVdrift","TRD drift velocities (local variations)");
  else object = new AliTRDCalPad("LocalT0","T0 (local variations)");
  
  //Read the Tree
  Int_t detector = -1;
  Float_t values[2304];
  tree->SetBranchAddress("detector",&detector);
  if(i == 0) tree->SetBranchAddress("gainPad",values);
  if(i == 1) tree->SetBranchAddress("vdrift",values);
  if(i == 3) tree->SetBranchAddress("t0",values);
  
  //Variables
  Float_t mean = 0.0;
  Int_t numberofentries = tree->GetEntries();
  
  if(numberofentries != 540) AliInfo("The tree is not complete");
  
  for (Int_t det= 0; det < numberofentries; ++det){
    tree->GetEntry(det);
    AliTRDCalROC *calROC = object->GetCalROC(detector);
    mean = detobject->GetValue(detector);
    if(mean == 0) continue;
    Int_t rowMax = calROC->GetNrows();
    Int_t colMax = calROC->GetNcols();
    for(Int_t row = 0; row < rowMax; ++row){
      for(Int_t col = 0; col < colMax; ++col){
	calROC->SetValue(col,row,TMath::Abs(values[(Int_t) (col*rowMax+row)])/mean);
      }//col
    }//row
  }
  return object;
}
//______________________________________________________________________________________________________
TObject* AliTRDCalibra::CreatePadObjectTree(TTree *tree)
{
  //
  // It Creates the AliTRDCalPad object from the tree of the coefficient for the calibration PRF (i = 2)
  // This object has to be written in the database
  //
  
  
  //Create the DetObject
  AliTRDCalPad *object = new AliTRDCalPad("PRFWidth","PRFWidth");

  //Read the Tree
  Int_t detector = -1;
  Float_t values[2304];
  tree->SetBranchAddress("detector",&detector);
  tree->SetBranchAddress("width",values);
   
  //Variables
  Int_t numberofentries = tree->GetEntries();

  if(numberofentries != 540) AliInfo("The tree is not complete");

  for (Int_t det= 0; det < numberofentries; ++det){
    tree->GetEntry(det);
    AliTRDCalROC *calROC = object->GetCalROC(detector);
    Int_t rowMax = calROC->GetNrows();
    Int_t colMax = calROC->GetNcols();
    for(Int_t row = 0; row < rowMax; ++row){
      for(Int_t col = 0; col < colMax; ++col){
	calROC->SetValue(col,row,TMath::Abs(values[(Int_t) (col*rowMax+row)]));
      }//col
    }//row
  }
  return object;
}
//______________________________________________________________________________________________________
TTree* AliTRDCalibra::Sum2Trees(const char* filename1, const char* filename2, const char *variablecali)
{
  //
  // It returns the sum of two trees with the name variablecali in the files filenam1 and filename2
  // equivalent of merging two 2D histos
  // The name of the resulting tree is the same as the two input trees
  // variablecali can be treeCH2d, treePH2d or treePRF2d 
  //
  
  
  //Variables
  TChain *treeChain = new TChain(variablecali);
  std::vector<Int_t> vectorplace;
  std::vector<std::vector<Int_t> > where;
 
 
  
  //first tree*****************************
  //Take the tree
  TFile *file1 = (TFile *) TFile::Open(filename1,"READ");
  TTree *tree1 = (TTree *) file1->Get(variablecali);
  //Take the places
  vectorplace = ConvertTreeVector(tree1);
  //Say where it is in tree 1
  for(Int_t jui = 0; jui < (Int_t) vectorplace.size(); jui++){
    std::vector<Int_t> chainplace;
    chainplace.push_back((Int_t) jui);
    where.push_back(chainplace);
  }
  //Add to the chain
  treeChain->Add(filename1);
  delete file1;
  



  //second tree*****************************
  //Take the tree
  TFile *file2 = (TFile *) TFile::Open(filename2,"READ");
  TTree *tree2 = (TTree *) file2->Get(variablecali);
  //Take the places
  std::vector<Int_t> vector2 = ConvertTreeVector(tree2);
  Int_t j = treeChain->GetEntries();
  for(Int_t jui = 0; jui < (Int_t) vector2.size(); jui++){
    //Search if already found
    Int_t place = SearchInTreeVector(vectorplace,vector2[jui]);
    //Create a new element in the two std vectors
    if(place == -1){
      std::vector<Int_t> chainplace;
      chainplace.push_back((Int_t) (j+jui));
      vectorplace.push_back((Int_t) (vector2[jui]));
      where.push_back(chainplace);
	   
    }
    //Update the element at the place "place" in the std vector whereinthechain
    else {
      std::vector<Int_t> chainplace = where[place];
      chainplace.push_back((Int_t) (j+jui));
      std::vector<std::vector<Int_t> >::iterator it = where.begin()+place;
      where.erase(it);
      where.insert(it,chainplace);
    }
  }
  //Add to the Chain
  treeChain->Add(filename2);
  delete file2; 
  

  //Take care of the profile
  const char* pattern = "P";
  TTree *tree = 0x0;

  if(!strstr(variablecali,pattern)){
    //Ready to read the chain
    TH1F *his = 0x0;
    treeChain->SetBranchAddress("histo",&his);
    
    //Initialise the final tree
    Int_t group = -1;
    TH1F *histsum = 0x0;
   
    tree = new TTree(variablecali,variablecali);
    tree->Branch("groupnumber",&group,"groupnumber/I");
    tree->Branch("histo","TH1F",&histsum,32000,0);
    
    //Init histsum
    if(treeChain->GetEntries() < 1) return tree1; 
    
    for(Int_t h = 0; h < (Int_t) vectorplace.size(); h++){
      group = vectorplace[h];
      std::vector<Int_t> chainplace = where[h];
      treeChain->GetEntry(chainplace[0]);
      //Init for the first time
      if(h == 0)  {
	histsum = new TH1F("","",his->GetXaxis()->GetNbins(),his->GetXaxis()->GetBinLowEdge(1),his->GetXaxis()->GetBinUpEdge(his->GetXaxis()->GetNbins()));
	 histsum->Sumw2();
      }
      //Reset for each new group
      histsum->SetEntries(0.0);
      for(Int_t l = 0; l <= histsum->GetXaxis()->GetNbins(); l++){
	histsum->SetBinContent(l,0.0);
	histsum->SetBinError(l,0.0);
      }
      histsum->Add(his,1);
      if((Int_t) chainplace.size() > 1){
	for(Int_t s = 1; s < (Int_t)chainplace.size(); s++){
	  treeChain->GetEntry(chainplace[s]);
	  histsum->Add(his,1);
	}
      }
      tree->Fill();
    }
  }
  else{
    //Ready to read the chain
    TGraphErrors *his = 0x0;
    treeChain->SetBranchAddress("histo",&his);
    
    //Initialise the final tree
    Int_t group = -1;
    TGraphErrors *histsum = 0x0;
    Double_t *xref = 0x0;
  
    tree = new TTree(variablecali,variablecali);
    tree->Branch("groupnumber",&group,"groupnumber/I");
    tree->Branch("histo","TGraphErrors",&histsum,32000,0);
    
    //Init histsum
    if(treeChain->GetEntries() < 1) return tree1; 
    
    for(Int_t h = 0; h < (Int_t) vectorplace.size(); h++){
      group = vectorplace[h];
      std::vector<Int_t> chainplace = where[h];
      treeChain->GetEntry(chainplace[0]);
      //Init or reset for a new group
      Int_t nbins = his->GetN();
      Double_t *x;
      x = new Double_t[nbins];
      xref = his->GetX();
      Double_t *ex;
      ex = new Double_t[nbins];
      Double_t *y;
      y = new Double_t[nbins];
      Double_t *ey;
      ey = new Double_t[nbins];
      
      for(Int_t lo = 0; lo < nbins; lo++){
	x[lo] = xref[lo];
	ex[lo] = 0.0;
	y[lo]  = 0.0;
	ey[lo] = 0.0;
      }
      delete histsum;
      histsum = new TGraphErrors(nbins,x,y,ex,ey);
     
          
      //Add the first
      histsum = AddProfiles(his,histsum);
      if((Int_t) chainplace.size() > 1){
	for(Int_t s = 1; s < (Int_t)chainplace.size(); s++){
	  treeChain->GetEntry(chainplace[s]);
	  histsum = AddProfiles(his,histsum);
	}
      }
      tree->Fill();
    }
  }
  
  
  return tree;
}
//____________Function fill 2D for the moment out of the code_____________________________________________________________________________________________

//____________Function fill 2D all objects from digits_________________________________________________________________________

 Bool_t AliTRDCalibra::Create2DDiSimOnline(Int_t iev1, Int_t iev2)
{
  //
  // Only for simulations, after the simulation, create the 2D histos from the digits stored in the file "TRD.Digits.root" 
  // Only for CH and PH
  //

  
  const Int_t kNplan    = 6;
  const Int_t kNcham    = 5;


  //RunLoader and so on*****************************************************************
  if (gAlice) {
    delete gAlice->GetRunLoader();
    delete gAlice;
    gAlice=0;
  }

 
  AliRunLoader* rl = AliRunLoader::Open("TRD_test.root");
  if (!rl) {
    return kFALSE;
  }
  
  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  if (!gAlice) {
    return kFALSE;
  }

  // Import the Trees for the event nEvent in the file
  rl->LoadKinematics();
  rl->GetEvent(0);
  rl->LoadHeader();
  
  AliLoader* loader = rl->GetLoader("TRDLoader");
  if (!loader) {
    AliInfo("No TRDLLoader found!");
    return kFALSE;
  }

  // Get the pointer to the TRD detector 
  AliTRD *trd = (AliTRD *) gAlice->GetDetector("TRD");
  if (!trd) {
    AliInfo("No TRD detector found");
    return kFALSE;
  }
  // Get the pointer to the geometry object
  AliTRDgeometry *geo;
  if (trd) {
    geo = trd->GetGeometry();
  }
  else {
    AliInfo("No TRD geometry found");
    return kFALSE;
  }


  //DB Setting********************************************************************************
  //DB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
 
  
  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam");
    return kFALSE;
  }
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

 
  // Some parameters*****************************************************************
  fTimeMax     = cal->GetNumberOfTimeBins();
  fSf = (Float_t) cal->GetSamplingFrequency();


  //Create the 2D histos corresponding to the pad group calibration mode***********************************************************************
  if(fCH2dOn) {
    AliInfo(Form("We will fill the CH2d histo with the pad calibration mode: Nz %d, and Nrphi %d", fNz[0], fNrphi[0]));
    
    //Calcul the number of Xbins
    fNtotal[0] = 0;
    ModePadCalibration(2,0);
    ModePadFragmentation(0,2,0,0);
    fdetChamb2[0] = fNfragz[0]*fNfragrphi[0];
    fNtotal[0] += 6*18*fdetChamb2[0];
    ModePadCalibration(0,0);
    ModePadFragmentation(0,0,0,0);
    fdetChamb0[0] = fNfragz[0]*fNfragrphi[0];
    fNtotal[0] += 6*4*18*fdetChamb0[0];
    AliInfo(Form("Total number of Xbins: %d", fNtotal[0]));
    
    
    //Create the 2D histo
    if(fHisto2d) CreateCH2d(fNtotal[0]);
    
  }
  
  if(fPH2dOn) {
    AliInfo(Form("We will fill the PH2d histo with the pad calibration mode: Nz %d, and Nrphi %d", fNz[1], fNrphi[1]));
    
    //Calcul the number of Xbins
    fNtotal[1] = 0;
    ModePadCalibration(2,1);
    ModePadFragmentation(0,2,0,1);
    fdetChamb2[1] = fNfragz[1]*fNfragrphi[1];
    fNtotal[1] += 6*18*fdetChamb2[1];
    ModePadCalibration(0,1);
    ModePadFragmentation(0,0,0,1);
    fdetChamb0[1] = fNfragz[1]*fNfragrphi[1];
    fNtotal[1] += 6*4*18*fdetChamb0[1];
    AliInfo(Form("Total number of Xbins: %d", fNtotal[1]));


    //Create the 2D histo
    if(fHisto2d) CreatePH2d(fNtotal[1]);

  }

  loader->LoadDigits();
  AliInfo("LoadDigits ");
  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();
  //digitsManager->SetDebug(1);

  //Loop on event**************************************************************************************
  for (Int_t ievent = iev1; ievent < iev2; ievent++) {
    AliInfo(Form("Process event %d",ievent));
    rl->GetEvent(ievent);
    if (!loader->TreeD()){
      AliInfo("loader Loading Digits ... ");
      loader->LoadDigits();
    }
    digitsManager->ReadDigits(loader->TreeD());
    AliInfo("digitsManager Read Digits Done");
    // Read the digits from the file
    if (!(digitsManager->ReadDigits(loader->TreeD()))) {
      return kFALSE;
    }


    //Loop on detector************************************************************************
    for(Int_t iSect = 0; iSect < 18; iSect++){
      for(Int_t iPlane = 0; iPlane < kNplan; iPlane++){
	for(Int_t iChamb = 0; iChamb < kNcham; iChamb++){
	  
	  //A little geometry:
	  Int_t iDet   = geo->GetDetector(iPlane,iChamb,iSect);
	  Int_t rowMax = parCom->GetRowMax(iPlane,iChamb,iSect);
	  Int_t colMax = parCom->GetColMax(iPlane);

	  //first Xbins of the detector
	  if(fCH2dOn) CalculXBins(iDet,0);
	  if(fPH2dOn) CalculXBins(iDet,1);
	  //fragmentation of idect
	  for(Int_t i = 0; i < 2; i++){
	    ModePadCalibration(iChamb,i);
	    ModePadFragmentation(iPlane, iChamb, iSect,i);
	  }
	  //In the cas of charge
	  Float_t *amptotal;
	  amptotal = new Float_t[fNfragrphi[0]*fNfragz[0]];
	  if(fCH2dOn){
	    for(Int_t k = 0; k < fNfragrphi[0]*fNfragz[0]; k++){
	      amptotal[k] = 0.0;
	    }
	  }

	  
	  // Loop through the detector pixel*****************************************
	  for (Int_t time = 0; time < fTimeMax; time++) {
	    for (Int_t  col = 0;  col <  colMax;  col++) {
	      for (Int_t  row = 0;  row <  rowMax;  row++) {
		

		//amplitude and position  in pad group 
		AliTRDdigit *digit    = digitsManager->GetDigit(row,col,time,iDet);
		Int_t        amp      = digit->GetAmp();
		Int_t posr[2] = {0,0};
		Int_t posc[2] = {0,0};
		if((fCH2dOn) && (fNnz[0] != 0)) posr[0] = (Int_t) row/fNnz[0];
		if((fCH2dOn) && (fNnrphi[0] != 0)) posc[0] = (Int_t) col/fNnrphi[0];
		if((fPH2dOn) && (fNnz[1] != 0)) posr[1] = (Int_t) row/fNnz[1];
		if((fPH2dOn) && (fNnrphi[1] != 0)) posc[1] = (Int_t) col/fNnrphi[1];

	
		// Total spectrum
		if(fCH2dOn){
		  if(amp < fThresholddigit) amp = 0;
		  amptotal[(Int_t) (fXbins[0]+posc[0]*fNfragz[0]+posr[0])] += amp;
		}
		if(fPH2dOn ) {
		  if(fHisto2d) fPH2d->Fill((fXbins[1]+posc[1]*fNfragz[1]+posr[1])+0.5,(Float_t)time/fSf,amp);
		  if(fVector2d) UpdateVectorPH((fXbins[1]+posc[1]*fNfragz[1]+posr[1]), time, amp);

		}
		

		//memory stuff
		delete digit;
		
	      }//boucle row
	    }//boucle col
	  }//boucle time


	  if(fCH2dOn ){
	    for(Int_t k = 0; k < fNfragz[0]*fNfragrphi[0]; k++){
	      if(fHisto2d) fCH2d->Fill(fXbins[0]+k+0.5, amptotal[k]/fRelativeScale); 
	      if(fVector2d) UpdateVectorCH(fXbins[0]+k, amptotal[k]/fRelativeScale);
	    }
	  }

	  delete amptotal;
	  

	  
	}//boucle chamber
      }//boucle plane
    }//boucle sect
   
    
    loader->UnloadDigits();  
    
  }//boucle event
  
  if(fDebug == 1){

    if(fPH2dOn && fHisto2d) PlotPH2d();
    if(fCH2dOn && fHisto2d) PlotCH2d();
       
   
  }
  
  if(fWrite[0] || fWrite[1]) {

    TFile *fout = TFile::Open(fWriteName,"UPDATE");
    //Check if the file could be opened
    if (!fout || !fout->IsOpen()) {
      AliInfo("<No File found!");
      return kFALSE;
    }
  

    if(fCH2dOn && fHisto2d && fWrite[0]) fout->WriteTObject(fCH2d);
    if(fPH2dOn && fHisto2d && fWrite[1]) fout->WriteTObject(fPH2d);

    if(fVector2d && fCH2dOn && fWrite[0]){
      TTree *treeCH2d = ConvertVectorCTTreeHisto(fVectorCH,fPlaCH,"treeCH2d");
      fout->WriteTObject(treeCH2d);
    }

    if(fVector2d && fPH2dOn && fWrite[1]){
      TTree *treePH2d = ConvertVectorPTreeHisto(fVectorPH,fPlaPH,"treePH2d");
      fout->WriteTObject(treePH2d);
    }
  
  }
 
  return kTRUE;

}
//_______________Function fill 2D all objects from Raw Data______________________________________________________________________

Bool_t AliTRDCalibra::Create2DRaDaOnline(Int_t iev1, Int_t iev2)
{
  //
  // After having written the RAW DATA in the current directory, create the 2D histos from these RAW DATA  
  // Only for CH and PH
  //
  
  const Int_t kNplan    = 6;
  const Int_t kNcham    = 5;
  TString dirname(".");
  
  //DB Setting***************************************************************
  //CDB setting
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
  
  
  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam");
    return kFALSE;
  }
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

 
  // Some parameters
  fTimeMax     = cal->GetNumberOfTimeBins();
  fSf = (Float_t) cal->GetSamplingFrequency();


  //Create the 2D histo corresponding to the pad group calibration mode***************************************************************
  if(fCH2dOn) {
    AliInfo(Form("We will fill the CH2d histo with the pad calibration mode: Nz %d, and Nrphi %d", fNz[0], fNrphi[0]));
    
    //Calcul the number of Xbins
    fNtotal[0] = 0;
    ModePadCalibration(2,0);
    ModePadFragmentation(0,2,0,0);
    fdetChamb2[0] = fNfragz[0]*fNfragrphi[0];
    fNtotal[0] += 6*18*fdetChamb2[0];
    ModePadCalibration(0,0);
    ModePadFragmentation(0,0,0,0);
    fdetChamb0[0] = fNfragz[0]*fNfragrphi[0];
    fNtotal[0] += 6*4*18*fdetChamb0[0];
    AliInfo(Form("Total number of Xbins: %d", fNtotal[0]));
    
    
    //Create the 2D histo
    if(fHisto2d) CreateCH2d(fNtotal[0]);
    
  }
  
  if(fPH2dOn) {
    AliInfo(Form("We will fill the PH2d histo with the pad calibration mode: Nz %d, and Nrphi %d", fNz[1], fNrphi[1]));
    
    //Calcul the number of Xbins
    fNtotal[1] = 0;
    ModePadCalibration(2,1);
    ModePadFragmentation(0,2,0,1);
    fdetChamb2[1] = fNfragz[1]*fNfragrphi[1];
    fNtotal[1] += 6*18*fdetChamb2[1];
    ModePadCalibration(0,1);
    ModePadFragmentation(0,0,0,1);
    fdetChamb0[1] = fNfragz[1]*fNfragrphi[1];
    fNtotal[1] += 6*4*18*fdetChamb0[1];
    AliInfo(Form("Total number of Xbins: %d", fNtotal[1]));


    //Create the 2D histo
    if(fHisto2d) CreatePH2d(fNtotal[1]);

  }
  
   
  AliTRDrawData *rawdata = new AliTRDrawData();
  AliInfo("AliTRDrawData object created ");
  
  
  //Loop on events*************************************************************
  for (Int_t ievent = iev1; ievent < iev2; ievent++) {

   
    //AliRawReaderFile
    AliRawReaderFile *readerfile = new AliRawReaderFile(dirname,ievent);
    if(!readerfile) {
      AliInfo("No readerfile found!");
      return kFALSE;
    }
 

    AliTRDdigitsManager *digitsManager = rawdata->Raw2Digits((AliRawReader *) readerfile);
    if(!digitsManager) {
      AliInfo("No DigitsManager done!");
      return kFALSE;
    }
 
    //Loop on detectors**********************************************************   
    for(Int_t iSect = 0; iSect < 18; iSect++){
      for(Int_t iPlane = 0; iPlane < kNplan; iPlane++){
	for(Int_t iChamb = 0; iChamb < kNcham; iChamb++){
	  
	  //A little geometry:
	  Int_t iDet   = AliTRDgeometry::GetDetector(iPlane,iChamb,iSect);
	  Int_t rowMax = parCom->GetRowMax(iPlane,iChamb,iSect);
	  Int_t colMax = parCom->GetColMax(iPlane);

	  //first Xbins of the detector
	  if(fCH2dOn) CalculXBins(iDet,0);
	  if(fPH2dOn) CalculXBins(iDet,1);
	  //fragmentation of idect
	  for(Int_t i = 0; i < 2; i++){
	    ModePadCalibration(iChamb,i);
	    ModePadFragmentation(iPlane, iChamb, iSect,i);
	  }
	  //In the cas of charge
	  Float_t *amptotal;
	  amptotal = new Float_t[fNfragrphi[0]*fNfragz[0]];
	  if(fCH2dOn){
	    for(Int_t k = 0; k < fNfragrphi[0]*fNfragz[0]; k++){
	      amptotal[k] = 0.0;
	    }
	  }
	 
	   
	 
	  
	  // Loop through the detector pixel*****************************************
	  for (Int_t time = 0; time < fTimeMax; time++) {
	    for (Int_t  col = 0;  col <  colMax;  col++) {
	      for (Int_t  row = 0;  row <  rowMax;  row++) {

		//amplitude and position of the digit
		AliTRDdigit *digit    = digitsManager->GetDigit(row,col,time,iDet);
		Int_t        amp      = digit->GetAmp();
		Int_t posr[2] = {0,0};
		Int_t posc[2] = {0,0};
		if((fCH2dOn) && (fNnz[0] != 0)) posr[0] = (Int_t) row/fNnz[0];
		if((fCH2dOn) && (fNnrphi[0] != 0)) posc[0] = (Int_t) col/fNnrphi[0];
		if((fPH2dOn) && (fNnz[1] != 0)) posr[1] = (Int_t) row/fNnz[1];
		if((fPH2dOn) && (fNnrphi[1] != 0)) posc[1] = (Int_t) col/fNnrphi[1];
		
		
		// Total spectrum
		if(fCH2dOn){
		  if(amp < fThresholddigit) amp = 0;
		  amptotal[(Int_t) (fXbins[0]+posc[0]*fNfragz[0]+posr[0])] += amp;
		}
		if(fPH2dOn ) {
		  if(fHisto2d) fPH2d->Fill((fXbins[1]+posc[1]*fNfragz[1]+posr[1])+0.5,(Float_t)time/fSf,amp);
		  if(fVector2d) UpdateVectorPH(fXbins[1]+posc[1]*fNfragz[1]+posr[1], time, amp);
		}


		delete digit;
	      }//boucle row
	    }//boucle col
	  }//boucle time



	  if(fCH2dOn ){
	    for(Int_t k = 0; k < fNfragz[0]*fNfragrphi[0]; k++){
	      if(fHisto2d) fCH2d->Fill(fXbins[0]+k+0.5, amptotal[k]/fRelativeScale);
	      if(fVector2d) UpdateVectorCH(fXbins[0]+k, amptotal[k]/fRelativeScale); 
	    }
	  }

	  delete amptotal;
	  

	}//boucle chamber
      }//boucle plane
    }//boucle sect
  
    delete digitsManager;
    delete readerfile;
    
  }//boucle event
  
  if(fDebug == 1){

    if(fPH2dOn && fHisto2d) PlotPH2d();
    if(fCH2dOn && fHisto2d) PlotCH2d();   
    
  }

 
  if(fWrite[0] || fWrite[1]) {
    
    TFile *fout = TFile::Open(fWriteName,"UPDATE");
    //Check if the file could be opened
    if (!fout || !fout->IsOpen()) {
      AliInfo("<No File found!");
      return kFALSE;
    }
  

    if(fCH2dOn && fHisto2d && fWrite[0]) fout->WriteTObject(fCH2d);
    if(fPH2dOn && fHisto2d && fWrite[1]) fout->WriteTObject(fPH2d);

    if(fVector2d && fCH2dOn && fWrite[0]){
      TTree *treeCH2d = ConvertVectorCTTreeHisto(fVectorCH,fPlaCH,"treeCH2d");
      fout->WriteTObject(treeCH2d);
    }

    if(fVector2d && fPH2dOn && fWrite[1]){
      TTree *treePH2d = ConvertVectorPTreeHisto(fVectorPH,fPlaPH,"treePH2d");
      fout->WriteTObject(treePH2d);
    }
  
  }
  
  return kTRUE;

}
//_________________________________Pad Calibration Public____________________________________________________

//________________Define the number of pads per group for one detector and one calibration________________________________________________________________
void AliTRDCalibra::ModePadCalibration(Int_t iChamb, Int_t i)
{
  //
  // Definition of the calibration mode
  // from Nz and Nrphi, the number of row and col pads per calibration groups are setted
  //


  fNnz[i] = 0;
  fNnrphi[i] = 0;
  
  if((fNz[i] == 0) && (iChamb == 2)) {
    fNnz[i] = 12;
  }
  if((fNz[i] == 0) && (iChamb != 2)) {
    fNnz[i] = 16;
  }  
 
  if((fNz[i] == 1) && (iChamb == 2)) {
    fNnz[i] = 6;
  }
  if((fNz[i] == 1) && (iChamb != 2)) {
    fNnz[i] = 8;
  }

  if((fNz[i] == 2) && (iChamb == 2)) {
    fNnz[i] = 3;
  }
  if((fNz[i] == 2) && (iChamb != 2)) {
    fNnz[i] = 4;
  }
  if(fNz[i] == 3) {
    fNnz[i] = 2;
  }
  if(fNz[i] == 4) {
    fNnz[i] = 1;
  }
   
  if(fNrphi[i] == 0) {
    fNnrphi[i] = 144;
  }
  if(fNrphi[i] == 1) {
    fNnrphi[i] = 72;
  } 
  
  if(fNrphi[i] == 2) {
    fNnrphi[i] = 36;
  } 
  
  if(fNrphi[i] == 3) {
    fNnrphi[i] = 18;
  } 

  if(fNrphi[i] == 4) {
    fNnrphi[i] = 9;
  } 
  if(fNrphi[i] == 5) {
    fNnrphi[i] = 4;
  } 
  if(fNrphi[i] == 6) {
    fNnrphi[i] = 1;
  } 
}
//__________Define the number of pad groups in one detector for one calibration___________________________________________________________________________________________
Bool_t AliTRDCalibra::ModePadFragmentation(Int_t iPlane,Int_t iChamb, Int_t iSect, Int_t i)
{
  //
  // Definition of the calibration mode
  // from the number of row and col pads per calibration groups the number of calibration groups are setted
  //

  fNfragz[i] = 0;
  fNfragrphi[i] = 0;
  
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }

  //A little geometry:
  Int_t rowMax = parCom->GetRowMax(iPlane,iChamb,iSect);
  Int_t colMax = parCom->GetColMax(iPlane);
  
  
  //The fragmentation
  if (fNnz[i] != 0)  fNfragz[i] = (Int_t) rowMax/fNnz[i];
 

  if(fNnrphi[i] != 0)  fNfragrphi[i] = (Int_t) colMax/fNnrphi[i];
  

  return kTRUE;
}

//______________________________________________________________________________________________________
//______________Protected Functions____________________________________________________________________________________

//_______________________Create the 2D histo to be filled online__________________________________________________

//_________________________________________________________________________________
void AliTRDCalibra::CreatePRF2d(Int_t Nn)
{
  //
  // Create the 2D histos
  //

  fPRF2d = new TProfile2D("PRF2d","", Nn, 0, Nn, fNumberBinPRF, -1.0, 1.0);
  fPRF2d->SetXTitle("Det/pad groups");
  fPRF2d->SetYTitle("Position x/W [pad width units]");
  fPRF2d->SetZTitle("Q_{i}/Q_{total}");
  fPRF2d->SetStats(0);
  fPRF2d->Sumw2();

}

//_________________________________________________________________________________
void AliTRDCalibra::CreatePH2d(Int_t Nn)
{
  //
  // Create the 2D histos
  //

  fPH2d = new TProfile2D("PH2d","", Nn, 0, Nn, fTimeMax, -0.5/fSf,(Float_t)(fTimeMax-0.5)/fSf);
  fPH2d->SetXTitle("Det/pad groups");
  fPH2d->SetYTitle("time [#mus]");
  fPH2d->SetZTitle("<PH> [a.u.]");
  fPH2d->SetStats(0);
  fPH2d->Sumw2();

}

//__________________________________________________________________________________
void AliTRDCalibra::CreateCH2d(Int_t Nn)
{
  //
  // Create the 2D histos
  //

  fCH2d = new TH2I("CH2d", "", Nn, 0, Nn, fNumberBinCharge, 0, 300);
  fCH2d->SetXTitle("Det/pad groups");
  fCH2d->SetYTitle("charge deposit [a.u]");
  fCH2d->SetZTitle("counts");
  fCH2d->SetStats(0);
  fCH2d->Sumw2();

}
//_______Offine tracking in the AliTRDtracker_________________________________________________
void AliTRDCalibra::FillTheInfoOfTheTrackCH()
{
  //
  // For the offline tracking or mcm tracklets
  // This function will be called in the functions UpdateHistogram... 
  // to fill the info of a track for the relativ gain calibration
  //
	
  Int_t Nb = 0;//nombre de zones traversees
  Int_t fd = -1;//premiere zone non nulle
  
  
  //See if the track goes through different zones
  for(Int_t k = 0; k < fNfragz[0]*fNfragrphi[0]; k++){
    if(famptotal[k] > 0.0) {
      Nb++;
      if(Nb == 1) fd = k;
    }
  }
 
  //if automatic scale
  if((fCountRelativeScale < 100) && (fRelativeScaleAuto)) {
    //Take only the one zone track
    if(Nb == 1){
      fRelativeScale += famptotal[fd]*0.014*0.01;
      fCountRelativeScale ++;
    }
  }
	  
  
  //We fill the CH2d after having scale with the first 100
  if((fCountRelativeScale >= 100) && (fRelativeScaleAuto)){
    //case of track with only one zone
    if(Nb == 1){
      if(fHisto2d) fCH2d->Fill(fXbins[0]+fd+0.5, famptotal[fd]/fRelativeScale);
      if(fVector2d) UpdateVectorCH(fXbins[0]+fd,famptotal[fd]/fRelativeScale);
    }//case 1 zone
    //case of track with two zones
    if(Nb == 2){
      //two zones voisines sinon rien!
      if(famptotal[fd] > 0.0 && famptotal[fd+1] > 0.0){
	//one of the two very big
	if(famptotal[fd] > fprocent*famptotal[fd+1]) {
	  if(fHisto2d) fCH2d->Fill(fXbins[0]+fd+0.5, famptotal[fd]/fRelativeScale);
	  if(fVector2d) UpdateVectorCH(fXbins[0]+fd,famptotal[fd]/fRelativeScale);
	}
	if(famptotal[fd+1] > fprocent*famptotal[fd])  {
	  if(fHisto2d) fCH2d->Fill(fXbins[0]+fd+1.5, famptotal[fd+1]/fRelativeScale);
	  if(fVector2d) UpdateVectorCH(fXbins[0]+fd,famptotal[fd+1]/fRelativeScale);
	}
      }
    }//case 2 zones
  }
  
  
  //Fill with no automatic scale
  if(!fRelativeScaleAuto){
    //case of track with only one zone
    if(Nb == 1){
      fNumberusedch[0]++;
      if(fHisto2d) fCH2d->Fill(fXbins[0]+fd+0.5, famptotal[fd]/fRelativeScale);
      if(fVector2d) UpdateVectorCH(fXbins[0]+fd,famptotal[fd]/fRelativeScale);
    }//case 1 zone
    //case of track with two zones
    if(Nb == 2){
      //two zones voisines sinon rien!
      //case 1
      if(famptotal[fd] > 0.0 && famptotal[fd+1] > 0.0){
	//one of the two very big
	if(famptotal[fd] > fprocent*famptotal[fd+1])  {
	  if(fHisto2d) fCH2d->Fill(fXbins[0]+fd+0.5, famptotal[fd]/fRelativeScale);
	  if(fVector2d) UpdateVectorCH(fXbins[0]+fd,famptotal[fd]/fRelativeScale);
	  fNumberusedch[1]++;
	}
	if(famptotal[fd+1] > fprocent*famptotal[fd])  {
	  if(fHisto2d) fCH2d->Fill(fXbins[0]+fd+1.5, famptotal[fd+1]/fRelativeScale);
	  if(fVector2d) UpdateVectorCH(fXbins[0]+fd+1,famptotal[fd+1]/fRelativeScale);
	  fNumberusedch[1]++;
	}
      }
      //case 2
      if(fNfragz[0] > 1){
	if(famptotal[fd] > 0.0){
	  if((fd+fNfragz[0])< (fNfragz[0]*fNfragrphi[0])){
	    if(famptotal[fd+fNfragz[0]] > 0.0){
	      
	      //one of the two very big
	      if(famptotal[fd] > fprocent*famptotal[fd+fNfragz[0]])  {
		if(fHisto2d) fCH2d->Fill(fXbins[0]+fd+0.5, famptotal[fd]/fRelativeScale);
		if(fVector2d) UpdateVectorCH(fXbins[0]+fd,famptotal[fd]/fRelativeScale);
		fNumberusedch[1]++;
	      }
	      if(famptotal[fd+fNfragz[0]] > fprocent*famptotal[fd])  {
		if(fHisto2d) fCH2d->Fill(fXbins[0]+fd+fNfragz[0]+0.5, famptotal[fd+fNfragz[0]]/fRelativeScale);
		fNumberusedch[1]++;
		if(fVector2d) UpdateVectorCH(fXbins[0]+fd+fNfragz[0],famptotal[fd+fNfragz[0]]/fRelativeScale);
	      }
	    }
	  }
	}
      }
    }//case 2 zones
  }
}
//_______Offine tracking in the AliTRDtracker_________________________________________________
void AliTRDCalibra::Resetfvariables()
{
  //
  // Reset values of famptotal, fPHvalue and fPHplace for the updateHistogram... functions
  //

  //Reset the good track******************************************************************
  fGoodTrack = kTRUE;
  
  //Reset the famptotal where we put value*****************************************************************
  if(fCH2dOn){
    for(Int_t k = 0; k < fNfragz[0]*fNfragrphi[0]; k++){
      famptotal[k] = 0.0;
    }//for
  }//if
  
  //Reset the fPHvalue*****************************************************************
  if(fPH2dOn){
    for(Int_t k = 0; k < fTimeMax; k++){
      fPHvalue[k] = -1.0;
      fPHplace[k] = -1;
    }//for
  }//if
}
//_______Offine tracking in the AliTRDtracker_________________________________________________
void AliTRDCalibra::FillTheInfoOfTheTrackPH()
{
  //
  // For the offline tracking or mcm tracklets
  // This function will be called in the functions UpdateHistogram... 
  // to fill the info of a track for the drift velocity  calibration
  //
    
  Int_t Nb = 1;//nombre de zones traversees 1, 2 ou plus de 3
  Int_t fd1 = -1;//premiere zone non nulle
  Int_t fd2 = -1;//deuxieme zone non nulle
  Int_t k1 = -1;//debut de la premiere zone
  Int_t k2 = -1;//debut de la seconde zone
  
  
  //See if the track goes through different zones
  for(Int_t k = 0; k < fTimeMax; k++){
    if(fPHvalue[k] > 0.0) {
      if(fd1 == -1) {
	fd1 = fPHplace[k];
	k1 = k;	      
      }
      if(fPHplace[k] != fd1) {
	if(fd2 == -1) {
	  k2 = k;
	  fd2 = fPHplace[k];
	  Nb = 2;
	}
	if(fPHplace[k] != fd2) Nb = 3;
      }
    }
  }
  
  //Fill 
  //case of track with only one zone
  if(Nb == 1){
    fNumberusedph[0]++;
    for(Int_t i = 0; i < fTimeMax; i++){
      if(fPHvalue[i] > 0.0) {
	if(fHisto2d) fPH2d->Fill((fXbins[1]+fPHplace[i])+0.5,(Float_t)i/fSf,(Float_t) fPHvalue[i]);
	if(fDebug == 13) {
	  AliInfo(Form("WRITE Nb %d ,place final: %d, fPHplace[i]: %d, i: %d, fPHvalue[i]: %f", Nb, fXbins[1]+fPHplace[i], fPHplace[i], i, fPHvalue[i]));
	}
	if(fVector2d) UpdateVectorPH(fXbins[1]+fPHplace[i],i,fPHvalue[i]);
      }
    }
  }//case 1 zone
  //case of track with two zones
  if(Nb == 2){
    //two zones voisines sinon rien!
    //case 1
    if((fd1 == fd2+1) || (fd2 == fd1+1)){
      //one of the two fast all the think
      if(k2 > (k1+fdifference)){
	fNumberusedph[1]++;
	for(Int_t i = k1; i < k2; i++){
	  if(fPHvalue[i] > 0.0) {
	    if(fHisto2d) fPH2d->Fill((fXbins[1]+fPHplace[i])+0.5,(Float_t)i/fSf,(Float_t) fPHvalue[i]);
	    if(fVector2d) UpdateVectorPH(fXbins[1]+fPHplace[i],i,fPHvalue[i]);
	  }
	}
      }
      if((k2+fdifference) < fTimeMax){
	fNumberusedph[1]++;
	for(Int_t i = k2; i < fTimeMax; i++){
	  if(fPHvalue[i] > 0.0) {
	    if(fHisto2d) fPH2d->Fill((fXbins[1]+fPHplace[i])+0.5,(Float_t)i/fSf,(Float_t) fPHvalue[i]);
	    if(fVector2d) UpdateVectorPH(fXbins[1]+fPHplace[i],i,fPHvalue[i]);
	  }
	}
      }
    }
    //two zones voisines sinon rien!
    if(fNfragz[1] > 1){
      //case 2
      if((fd1+fNfragz[1]) < (fNfragz[1]*fNfragrphi[1])){
	if(fd2 == (fd1+fNfragz[1])){
	  //one of the two fast all the think
	  if(k2 > (k1+fdifference)){
	    fNumberusedph[1]++;
	    for(Int_t i = k1; i < k2; i++){
	      if(fPHvalue[i] > 0.0) {
		if(fHisto2d) fPH2d->Fill((fXbins[1]+fPHplace[i])+0.5,(Float_t)i/fSf,(Float_t) fPHvalue[i]);
		if(fVector2d) UpdateVectorPH(fXbins[1]+fPHplace[i],i,fPHvalue[i]);
	      }
	    }
	  }
	  if((k2+fdifference) < fTimeMax){
	    fNumberusedph[1]++;
	    for(Int_t i = k2; i < fTimeMax; i++){
	      if(fPHvalue[i] > 0.0) {
		if(fHisto2d) fPH2d->Fill((fXbins[1]+fPHplace[i])+0.5,(Float_t)i/fSf,(Float_t) fPHvalue[i]);
		if(fVector2d) UpdateVectorPH(fXbins[1]+fPHplace[i],i,fPHvalue[i]);
	      }
	    }
	  }
	}
      }
      //two zones voisines sinon rien!
      //case 3
      if((fd1-fNfragz[1]) >= 0){
	if(fd2 == (fd1-fNfragz[1])){
	  //one of the two fast all the think
	  if(k2 > (k1+fdifference)){
	    fNumberusedph[1]++;
	    for(Int_t i = k1; i < k2; i++){
	      if(fPHvalue[i] > 0.0) {
		if(fHisto2d) fPH2d->Fill((fXbins[1]+fPHplace[i])+0.5,(Float_t)i/fSf,(Float_t) fPHvalue[i]);
		if(fVector2d) UpdateVectorPH(fXbins[1]+fPHplace[i],i,fPHvalue[i]);
	      }
	    }
	  }
	  if((k2+fdifference) < fTimeMax){
	    fNumberusedph[1]++;
	    for(Int_t i = k2; i < fTimeMax; i++){
	      if(fPHvalue[i] > 0.0) {
		if(fHisto2d) fPH2d->Fill((fXbins[1]+fPHplace[i])+0.5,(Float_t)i/fSf,(Float_t) fPHvalue[i]);
		if(fVector2d) UpdateVectorPH(fXbins[1]+fPHplace[i],i,fPHvalue[i]);
	      }
	    }
	  }
	}
      }
    }
  }//case 2 zones
}

//_____________Set the pad calibration variables for the detector_____________________________________________________________________

Bool_t AliTRDCalibra::LocalisationdetectorXbins(Int_t detector)
{
  //
  //For the detector calcul the first Xbins and set the number of row and col pads per calibration groups, the number of calibration groups in the detector.
  //
  
  //first Xbins of the detector
  if(fCH2dOn) CalculXBins(detector,0);
  if(fPH2dOn) CalculXBins(detector,1);
  if(fPRF2dOn) CalculXBins(detector,2);
  //fragmentation of idect
  for(Int_t i = 0; i < 3; i++){
    ModePadCalibration((Int_t) GetChamber(detector),i);
    ModePadFragmentation((Int_t) GetPlane(detector), (Int_t) GetChamber(detector), (Int_t) GetSector(detector),i);
  }
  
  return kTRUE;
}

//____________________________Plot the 2D histos filled Online______________________________________________________________________

//___________________________________________________________________________________________
void AliTRDCalibra::PlotPH2d()
{
  //
  // Plot the 2D histo 
  //

  TCanvas *cph2d = new TCanvas("cph2d","",50,50,600,800);
  cph2d->cd();
  fPH2d->Draw("LEGO");

}

//______________________________________________________________________________________________
void AliTRDCalibra::PlotCH2d()
{
  //
  // Plot the 2D histos
  //

  TCanvas *cch2d = new TCanvas("cch2d","",50,50,600,800);
  cch2d->cd();
  fCH2d->Draw("LEGO");
  
}

//______________________________________________________________________________________________
void AliTRDCalibra::PlotPRF2d()
{
  //
  // Plot the 2D histos
  //

  TCanvas *cPRF2d = new TCanvas("cPRF2d","",50,50,600,800);
  cPRF2d->cd();
  fPRF2d->Draw("LEGO");
      
}

//_________Fit_______________________________________________________________________

//______________________Create histos if fDebug == 1 or fDebug >= 3_____________________________________________________________________

//___________________________________________________________________________________________
void AliTRDCalibra::CreateFitHistoPH(Int_t Nbins, Double_t Low, Double_t High)
{
  //
  // Create the histos for fDebug = 1 and fDebug = 4 (Fit functions)
  //

  //****Histograms to store the coef
  fCoefVdrift[0] = new TH1F("coefvdrift0","",Nbins,Low,High);
  fCoefVdrift[1] = new TH1F("coefvdrift1","",Nbins,Low,High);
  fCoefVdrift[2] = new TH1F("coefvdrift2","",Nbins,Low,High);
 

  //Histograms for Debug 
  fDeltaVdrift[0] = new TH1F("deltavdrift0","",Nbins,Low,High);
  fDeltaVdrift[1] = new TH1F("deltavdrift1","",Nbins,Low,High);
  fErrorVdrift[0] = new TH1I("errorvdrift0","",300,-0.5,0.5);
  fErrorVdrift[1] = new TH1I("errorvdrift1","",300,-0.5,0.5);
 



  fCoefVdrift[0]->SetXTitle("Det/pad groups");
  fCoefVdrift[0]->SetYTitle("Vdrift [cm/#mus]");
  fCoefVdrift[1]->SetXTitle("Det/pad groups");
  fCoefVdrift[1]->SetYTitle("Vdrift [cm/#mus]");
  fCoefVdrift[2]->SetXTitle("Det/pad groups");
  fCoefVdrift[2]->SetYTitle("Vdrift [cm/#mus]");


  fDeltaVdrift[0]->SetXTitle("Det/pad groups");
  fDeltaVdrift[0]->SetYTitle("#Deltav/v_{sim}");
  fDeltaVdrift[1]->SetXTitle("Det/pad groups");
  fDeltaVdrift[1]->SetYTitle("#Deltav/v_{sim}");


  fErrorVdrift[0]->SetXTitle("#Deltav/v_{sim}");
  fErrorVdrift[0]->SetYTitle("counts");
  fErrorVdrift[1]->SetXTitle("#Deltav/v_{sim}");
  fErrorVdrift[1]->SetYTitle("counts");


  fCoefVdrift[0]->SetStats(0);
  fCoefVdrift[1]->SetStats(0);
  fCoefVdrift[2]->SetStats(0);
  fDeltaVdrift[0]->SetStats(0);
  fDeltaVdrift[1]->SetStats(0);
  fErrorVdrift[0]->SetStats(0);
  fErrorVdrift[1]->SetStats(0);
 
  fCoefVdrift[0]->SetMarkerColor(6);
  fCoefVdrift[0]->SetMarkerStyle(26);
  fCoefVdrift[0]->SetLineColor(6);
  fCoefVdrift[1]->SetMarkerColor(2);
  fCoefVdrift[1]->SetMarkerStyle(24);
  fCoefVdrift[1]->SetLineColor(2);
  fCoefVdrift[2]->SetLineColor(4);
  
 
  fDeltaVdrift[1]->SetMarkerColor(2);
  fDeltaVdrift[1]->SetMarkerStyle(24);
  fDeltaVdrift[1]->SetLineColor(2);
  fDeltaVdrift[0]->SetMarkerColor(6);
  fDeltaVdrift[0]->SetMarkerStyle(26);
  fDeltaVdrift[0]->SetLineColor(6); 


  fErrorVdrift[1]->SetLineColor(2);
  fErrorVdrift[1]->SetLineStyle(2);
  fErrorVdrift[0]->SetLineColor(6);
  fErrorVdrift[0]->SetLineStyle(1);
  
  
}
//___________________________________________________________________________________________
void AliTRDCalibra::CreateFitHistoT0(Int_t Nbins, Double_t Low, Double_t High)
{
  //
  // Create the histos for fDebug = 1 and fDebug = 4 (Fit functions)
  //

  //****Histograms to store the coef
  fCoefT0[0] = new TH1F("coefT00","",Nbins,Low,High);
  fCoefT0[1] = new TH1F("coefT01","",Nbins,Low,High);
  fCoefT0[2] = new TH1F("coefT02","",Nbins,Low,High);
 

  //Histograms for Debug 
  fDeltaT0[0] = new TH1F("deltaT00","",Nbins,Low,High);
  fDeltaT0[1] = new TH1F("deltaT01","",Nbins,Low,High);
  fErrorT0[0] = new TH1I("errorT00","",300,-0.1,0.1);
  fErrorT0[1] = new TH1I("errorT01","",300,-0.1,0.1);
 



  fCoefT0[0]->SetXTitle("Det/pad groups");
  fCoefT0[0]->SetYTitle("t0 [#mus]");
  fCoefT0[1]->SetXTitle("Det/pad groups");
  fCoefT0[1]->SetYTitle("t0 [#mus]");
  fCoefT0[2]->SetXTitle("Det/pad groups");
  fCoefT0[2]->SetYTitle("t0 [#mus]");


  fDeltaT0[0]->SetXTitle("Det/pad groups");
  fDeltaT0[0]->SetYTitle("#Deltat0 [#mus]");
  fDeltaT0[1]->SetXTitle("Det/pad groups");
  fDeltaT0[1]->SetYTitle("#Deltat0 [#mus]");


  fErrorT0[0]->SetXTitle("#Deltat0 [#mus]");
  fErrorT0[0]->SetYTitle("counts");
  fErrorT0[1]->SetXTitle("#Deltat0 [#mus]");
  fErrorT0[1]->SetYTitle("counts");


  fCoefT0[0]->SetStats(0);
  fCoefT0[1]->SetStats(0);
  fCoefT0[2]->SetStats(0);
  fDeltaT0[0]->SetStats(0);
  fDeltaT0[1]->SetStats(0);
  fErrorT0[0]->SetStats(0);
  fErrorT0[1]->SetStats(0);
 
  fCoefT0[0]->SetMarkerColor(6);
  fCoefT0[0]->SetMarkerStyle(26);
  fCoefT0[0]->SetLineColor(6);
  fCoefT0[1]->SetMarkerColor(2);
  fCoefT0[1]->SetMarkerStyle(24);
  fCoefT0[1]->SetLineColor(2);
  fCoefT0[2]->SetLineColor(4);
  
 
  fDeltaT0[1]->SetMarkerColor(2);
  fDeltaT0[1]->SetMarkerStyle(24);
  fDeltaT0[1]->SetLineColor(2);
  fDeltaT0[0]->SetMarkerColor(6);
  fDeltaT0[0]->SetMarkerStyle(26);
  fDeltaT0[0]->SetLineColor(6); 


  fErrorT0[1]->SetLineColor(2);
  fErrorT0[1]->SetLineStyle(2);
  fErrorT0[0]->SetLineColor(6);
  fErrorT0[0]->SetLineStyle(1);
  
  
}

//___________________________________________________________________________________________
void AliTRDCalibra::CreateFitHistoCH(Int_t Nbins, Double_t Low, Double_t High)
{
  //
  // Create the histos for fDebug = 1 and fDebug = 4 (Fit functions)
  //

  //Histograms to store the coef
  fCoefCharge[0] = new TH1F("coefcharge0","",Nbins, Low, High);
  fCoefCharge[1] = new TH1F("coefcharge1","",Nbins, Low, High);
  fCoefCharge[2] = new TH1F("coefcharge2","",Nbins, Low, High);
  fCoefCharge[3] = new TH1F("coefcharge3","",Nbins, Low, High);
 

  //Histograms for Debug 
  fDeltaCharge[0] = new TH1F("deltacharge0","",Nbins, Low, High);
  fDeltaCharge[1] = new TH1F("deltacharge1","",Nbins, Low, High);
  fDeltaCharge[2] = new TH1F("deltacharge2","",Nbins, Low, High);
 
 
  fErrorCharge[0] = new TH1I("errorcharge0","",100,-0.5,0.5);
  fErrorCharge[1] = new TH1I("errorcharge1","",100,-0.5,0.5);
  fErrorCharge[2] = new TH1I("errorcharge2","",100,-0.5,0.5);


  fCoefCharge[0]->SetXTitle("Det/Pad groups");
  fCoefCharge[0]->SetYTitle("gain factor");
  fCoefCharge[1]->SetXTitle("Det/Pad groups");
  fCoefCharge[1]->SetYTitle("gain factor");
  fCoefCharge[2]->SetXTitle("Det/Pad groups");
  fCoefCharge[2]->SetYTitle("gain factor");
  fCoefCharge[3]->SetXTitle("Det/Pad groups");
  fCoefCharge[3]->SetYTitle("gain factor");

  fDeltaCharge[0]->SetXTitle("Det/Pad groups");
  fDeltaCharge[0]->SetYTitle("#Deltag/g_{sim}");
  fDeltaCharge[1]->SetXTitle("Det/Pad groups");
  fDeltaCharge[1]->SetYTitle("#Deltag/g_{sim}");
  fDeltaCharge[2]->SetXTitle("Det/Pad groups");
  fDeltaCharge[2]->SetYTitle("#Deltag/g_{sim}");
  fDeltaCharge[0]->SetAxisRange(-0.5,0.5,"Y");
  fDeltaCharge[1]->SetAxisRange(-0.5,0.5,"Y");
  fDeltaCharge[2]->SetAxisRange(-0.5,0.5,"Y");

  fErrorCharge[0]->SetXTitle("#Deltag/g_{sim}");
  fErrorCharge[0]->SetYTitle("counts"); 
  fErrorCharge[1]->SetXTitle("#Deltag/g_{sim}");
  fErrorCharge[1]->SetYTitle("counts"); 
  fErrorCharge[2]->SetXTitle("#Deltag/g_{sim}");
  fErrorCharge[2]->SetYTitle("counts"); 
 

  fDeltaCharge[1]->SetMarkerColor(2);
  fDeltaCharge[1]->SetMarkerStyle(24);
  fDeltaCharge[1]->SetLineColor(2);
  fErrorCharge[1]->SetLineColor(2);
  fErrorCharge[1]->SetLineStyle(2);
  fDeltaCharge[2]->SetMarkerColor(8);
  fDeltaCharge[2]->SetLineColor(8);
  fDeltaCharge[2]->SetMarkerStyle(9);
  fErrorCharge[2]->SetLineColor(8);
  fErrorCharge[2]->SetLineStyle(5);
  fDeltaCharge[0]->SetMarkerColor(6);
  fDeltaCharge[0]->SetLineColor(6);
  fDeltaCharge[0]->SetMarkerStyle(26);
  fErrorCharge[0]->SetLineColor(6);
  fErrorCharge[0]->SetLineStyle(1);

  fCoefCharge[3]->SetLineColor(4);
  fCoefCharge[1]->SetMarkerColor(2);
  fCoefCharge[1]->SetLineColor(2);
  fCoefCharge[1]->SetMarkerStyle(24);
  fCoefCharge[2]->SetMarkerColor(8);
  fCoefCharge[2]->SetLineColor(8);
  fCoefCharge[2]->SetMarkerStyle(9);
  fCoefCharge[0]->SetMarkerColor(6);
  fCoefCharge[0]->SetLineColor(6);
  fCoefCharge[0]->SetMarkerStyle(26);
 
  fErrorCharge[2]->SetLineWidth(3);

 
  fDeltaCharge[1]->SetStats(0);
  fDeltaCharge[2]->SetStats(0);
  fDeltaCharge[0]->SetStats(0);
  fErrorCharge[1]->SetStats(0);
  fErrorCharge[2]->SetStats(0);
  fErrorCharge[0]->SetStats(0);
  fCoefCharge[1]->SetStats(0);
  fCoefCharge[0]->SetStats(0);
  fCoefCharge[3]->SetStats(0); 
  fCoefCharge[2]->SetStats(0);
    

}
//___________________________________________________________________________________________
void AliTRDCalibra::CreateFitHistoPRF(Int_t Nbins, Double_t Low, Double_t High)
{
  //
  // Create the histos for fDebug = 1 and fDebug = 4 (Fit functions)
  //

  //Histograms to store the coef
  fCoefPRF[0] = new TH1F("coefPRF0","",Nbins, Low, High);
  fCoefPRF[1] = new TH1F("coefPRF1","",Nbins, Low, High);
 
  //Histograms for Debug 
  fDeltaPRF = new TH1F("deltaPRF","",Nbins, Low, High);

  fErrorPRF = new TH1I("errorPRF","",300,-0.5,0.5);
 
  fDeltaPRF->SetMarkerColor(6);
  fDeltaPRF->SetMarkerStyle(26);
  fDeltaPRF->SetLineColor(6);
  fErrorPRF->SetLineColor(6);
  fErrorPRF->SetLineStyle(2);

  fCoefPRF[1]->SetLineColor(4);
  fCoefPRF[0]->SetMarkerColor(6);
  fCoefPRF[0]->SetMarkerStyle(26);
  fCoefPRF[0]->SetLineColor(6);

  fCoefPRF[0]->SetXTitle("Det/Pad groups");
  fCoefPRF[0]->SetYTitle("#sigma_{PRF}"); 
  fCoefPRF[1]->SetXTitle("Det/Pad groups");
  fCoefPRF[1]->SetYTitle("#sigma_{PRF}"); 

  fDeltaPRF->SetXTitle("Det/Pad groups");
  fDeltaPRF->SetYTitle("#Delta#sigma/#sigma_{sim}");

  fErrorPRF->SetXTitle("#Delta#sigma/#sigma_{sim}");
  fErrorPRF->SetYTitle("counts");  

 
  fDeltaPRF->SetStats(0);
  fErrorPRF->SetStats(0);
  fCoefPRF[1]->SetStats(0);
  fCoefPRF[0]->SetStats(0);

}
//_________________________________________________________________________________________
void AliTRDCalibra::CreateFitHistoPRFDB(Int_t rowMax, Int_t colMax)
{
  //
  // Create the histos for fDebug = 3 and fDebug = 4 (Fit functions)
  //
  fCoefPRFDB = new TH2F("coefPRF","", rowMax, 0,rowMax, colMax, 0, colMax);

  fCoefPRFDB->SetStats(0);
  fCoefPRFDB->SetXTitle("row Number");
  fCoefPRFDB->SetYTitle("col Number");
  fCoefPRFDB->SetZTitle("PRF width [pad width units]");
 

  fCoefPRFDB->SetFillColor(6);
  fCoefPRFDB->SetLineColor(6);


}

//_________________________________________________________________________________________
void AliTRDCalibra::CreateFitHistoCHDB(Int_t rowMax, Int_t colMax){
  //
  // Create the histos for fDebug = 3 and fDebug = 4 (Fit functions)
  //

  fCoefChargeDB[0] = new TH2F("coefchargedb0","", rowMax, 0,rowMax, colMax, 0, colMax);
  fCoefChargeDB[1] = new TH2F("coefchargedb1","", rowMax, 0,rowMax, colMax, 0, colMax);
  fCoefChargeDB[2] = new TH2F("coefchargedb2","", rowMax, 0,rowMax, colMax, 0, colMax);

  fCoefChargeDB[0]->SetStats(0);
  fCoefChargeDB[1]->SetStats(0);
  fCoefChargeDB[2]->SetStats(0);
  fCoefChargeDB[0]->SetXTitle("row Number");
  fCoefChargeDB[0]->SetYTitle("col Number");
  fCoefChargeDB[1]->SetXTitle("row Number");
  fCoefChargeDB[1]->SetYTitle("col Number");
  fCoefChargeDB[2]->SetXTitle("row Number");
  fCoefChargeDB[2]->SetYTitle("col Number");
  fCoefChargeDB[0]->SetZTitle("f_{g} Fit method");
  fCoefChargeDB[1]->SetZTitle("f_{g} Mean method");
  fCoefChargeDB[2]->SetZTitle("f_{g} Fitbis method");
 


  fCoefChargeDB[0]->SetFillColor(6);
  fCoefChargeDB[0]->SetLineColor(6);
  fCoefChargeDB[0]->SetLineColor(6);
  fCoefChargeDB[1]->SetFillColor(2);
  fCoefChargeDB[1]->SetLineColor(2);
  fCoefChargeDB[1]->SetLineColor(2);
  fCoefChargeDB[2]->SetFillColor(8);
  fCoefChargeDB[2]->SetLineColor(8);
  fCoefChargeDB[2]->SetLineColor(8);

}

//_________________________________________________________________________________________
void AliTRDCalibra::CreateFitHistoPHDB(Int_t rowMax, Int_t colMax)
{
  //
  // Create the histos for fDebug = 3 and fDebug = 4 (Fit functions)
  //

  fCoefVdriftDB[0] = new TH2F("coefvdriftdb0","", rowMax, 0,rowMax, colMax, 0, colMax);
  fCoefVdriftDB[1] = new TH2F("coefvdriftdb1","", rowMax, 0,rowMax, colMax, 0, colMax);
 
  
  fCoefVdriftDB[0]->SetStats(0);
  fCoefVdriftDB[1]->SetStats(0);
  fCoefVdriftDB[0]->SetXTitle("row Number");
  fCoefVdriftDB[0]->SetYTitle("col Number");
  fCoefVdriftDB[1]->SetXTitle("row Number");
  fCoefVdriftDB[1]->SetYTitle("col Number");
  fCoefVdriftDB[0]->SetZTitle("v_{drift} Fit method");
  fCoefVdriftDB[1]->SetZTitle("v_{drift} slope method");
  
  fCoefVdriftDB[0]->SetFillColor(6);
  fCoefVdriftDB[0]->SetLineColor(6);
  fCoefVdriftDB[0]->SetLineColor(6);
  fCoefVdriftDB[1]->SetFillColor(2);
  fCoefVdriftDB[1]->SetLineColor(2);
  fCoefVdriftDB[1]->SetLineColor(2);
  
}
//_________________________________________________________________________________________
void AliTRDCalibra::CreateFitHistoT0DB(Int_t rowMax, Int_t colMax)
{
  //
  // Create the histos for fDebug = 3 and fDebug = 4 (Fit functions)
  //

  fCoefT0DB[0] = new TH2F("coefT0db0","", rowMax, 0,rowMax, colMax, 0, colMax);
  fCoefT0DB[1] = new TH2F("coefT0db1","", rowMax, 0,rowMax, colMax, 0, colMax);
 
  
  fCoefT0DB[0]->SetStats(0);
  fCoefT0DB[1]->SetStats(0);
  fCoefT0DB[0]->SetXTitle("row Number");
  fCoefT0DB[0]->SetYTitle("col Number");
  fCoefT0DB[1]->SetXTitle("row Number");
  fCoefT0DB[1]->SetYTitle("col Number");
  fCoefT0DB[0]->SetZTitle("t0 Fit method");
  fCoefT0DB[1]->SetZTitle("t0 slope method");
  
  fCoefT0DB[0]->SetFillColor(6);
  fCoefT0DB[0]->SetLineColor(6);
  fCoefT0DB[0]->SetLineColor(6);
  fCoefT0DB[1]->SetFillColor(2);
  fCoefT0DB[1]->SetLineColor(2);
  fCoefT0DB[1]->SetLineColor(2);
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibra::FillVectorFitCH(Int_t countdet)
{
  //
  // For the Fit functions fill the vector FitCH special for the gain calibration
  //
  TFitCHInfo *fFitCHInfo = new TFitCHInfo();
  Int_t Ntotal = 1;
  if(GetChamber(countdet) == 2) Ntotal = 1728;
  else Ntotal = 2304;
  fFitCHInfo->fcoef = new Float_t[Ntotal];
  for(Int_t i = 0; i < Ntotal; i++){
    fFitCHInfo->fcoef[i] = fcoefCH[i];
  }
  fFitCHInfo->fDetector = countdet;
  fVectorFitCH.push_back(fFitCHInfo);

  return kTRUE;
}
//____________Functions for initialising the AliTRDCalibra in the code_________________________________________________________________________

Bool_t AliTRDCalibra::InitFit(Int_t Nbins, Double_t lowedge, Double_t upedge, Int_t i)
{
  //
  //Init the calibration mode (Nz, Nrphi), the histograms for debugging the fit methods if fDebug > 0, 
  //
  
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.01);
  
  //DB Setting **************************************************************************************
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    AliInfo("Could not get CDB Manager");
    return kFALSE;
  }
  
  
  //Get cal
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }
  
  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }
  
  
  
  
  // Some parameters
  fTimeMax     = cal->GetNumberOfTimeBins();
  fSf =  cal->GetSamplingFrequency();
  
  
  //Mode groups of pads: the total number of bins!*******************************************************
  Int_t numberofbinsexpected = 0;
  ModePadCalibration(2,i);
  ModePadFragmentation(0,2,0,i);
  fdetChamb2[i] = fNfragz[i]*fNfragrphi[i];
  if(fDebug == 1) AliInfo(Form("For the chamber 2: %d", fdetChamb2[i]));
  numberofbinsexpected += 6*18*fdetChamb2[i];
  ModePadCalibration(0,i);
  ModePadFragmentation(0,0,0,i);
  fdetChamb0[i] = fNfragz[i]*fNfragrphi[i];
  if(fDebug == 1) AliInfo(Form("For the other chamber 0: %d", fdetChamb0[i]));
  numberofbinsexpected += 6*4*18*fdetChamb0[i];
  
  
  //Quick verification that we have the good pad calibration mode if 2D histos!
  if(Nbins != 0){
    if(numberofbinsexpected != Nbins){
      AliInfo("It doesn't correspond to the mode of pad group calibration!");
      return kFALSE;
    }
  }

  //Security for fDebug 3 and 4
  if((fDebug >= 3) && ((fDet[0] > 5) || (fDet[1] > 4) || (fDet[2] > 17))){
    AliInfo("This detector doesn't exit!");
    return kFALSE;
    
  }
  
  
  //Determine fDet1 and fDet2***************************************************************
  fdect1[i] = -1;
  fdect2[i] = -1;
  if(fDebug == 2) {
    fdect1[i] = fFitVoir;
    fdect2[i] = fdect1[i] +1;
  }
  if(fDebug <= 1){
    fdect1[i] = 0;
    fdect2[i] = numberofbinsexpected;
  }
  if(fDebug >= 3){
    CalculXBins(AliTRDgeometry::GetDetector(fDet[0],fDet[1],fDet[2]),i);
    fdect1[i] = fXbins[i];
    CalculXBins((AliTRDgeometry::GetDetector(fDet[0],fDet[1],fDet[2])+1),i);
    fdect2[i] = fXbins[i];
  }
  
  //Create the histos for debugging***************************************************************
  
  //CH
  if(i == 0) {
    
    //Init the VectorFitCH********************************************************
    fcoefCH = new Float_t[2304];
    for(Int_t k = 0; k < 2304; k++){
      fcoefCH[k] = 0.0;    
    }
    fScalefitfactor = 0.0;
    
    //Number of Xbins(detectors or groups of pads) if Vector2d*******************************
    //Quick verification that we are not out of range!
    if(Nbins == 0){
      if((Int_t) fVectorCH.size() > numberofbinsexpected){
	AliInfo("ch doesn't correspond to the mode of pad group calibration!");
	return kFALSE;
      }
      
      if((Int_t) fVectorCH.size() != (Int_t) fPlaCH.size()){
	AliInfo("VectorCH doesn't correspond to PlaCH!");
	return kFALSE;
      }
    }
          
    //Debugging____Create the histos********************************************
    
    //fDebug == 0 nothing
    
    //fDebug == 1 
    if(fDebug == 1){
      if(Nbins != 0){
	//create the histos replique de ch if histos2D
	CreateFitHistoCH(Nbins, lowedge, upedge);
      }
      else{
	//create the histos replique de ch vector2d
	CreateFitHistoCH(numberofbinsexpected, 0, numberofbinsexpected);
      }
    }
    
    
    //fDebug == 2 and fFitVoir no histo
    if(fDebug == 2){
      if(fFitVoir < numberofbinsexpected) {
	AliInfo(Form("We will see the fit of the object %d",fFitVoir));
      }
      else {
	AliInfo("fFitVoir is out of range of the histo!");
	return kFALSE;
      }
    }
     
    
    //fDebug == 3  or 4 and fDet
    if(fDebug >= 3){
      if((fNz[0] == 0) && (fNrphi[0] == 0)) {
	AliInfo("Do you really want to see one detector without pad groups?");
	return kFALSE;
      }
      else  {
	AliInfo(Form("You will see the detector: iPlane %d, iChamb %d, iSect %d",fDet[0], fDet[1], fDet[2]));
	



	//A little geometry:
	Int_t rowMax = parCom->GetRowMax(fDet[0],fDet[1],fDet[2]);
	Int_t colMax = parCom->GetColMax(fDet[0]);
	

	//Create the histos to visualise
	CreateFitHistoCHDB(rowMax, colMax);
	if(fDebug == 4) CreateFitHistoCH((Int_t) (fdect2[0]-fdect1[0]), fdect1[0], fdect2[0]);
	
      }
    }
  }
    
    
    //PH and T0
  if(i == 1) {
    
    //Number of Xbins(detectors or groups of pads) if vector2d*******************************
    //Quick verification that we are not out of range!
    if(Nbins == 0){
      if((Int_t) fVectorPH.size() > numberofbinsexpected){
	AliInfo("ch doesn't correspond to the mode of pad group calibration!");
	return kFALSE;
      }
      if((Int_t) fVectorPH.size() != (Int_t) fPlaPH.size()){
	AliInfo("VectorPH doesn't correspond to PlaPH!");
	return kFALSE;
      }
    }
        
    //init tree**********************************************
    
    InittreePH();
    InittreeT0();
        
    
    
    //Debugging____Create the histos********************************************
    
    //fDebug == 0 nothing
    
    //fDebug == 1 
    if(fDebug == 1){
      if(Nbins != 0){
	//create the histos replique de ch
	CreateFitHistoPH(Nbins,lowedge,upedge);
	CreateFitHistoT0(Nbins,lowedge,upedge);
      }
      else {
	//create the histos replique de ch if vector2d
	CreateFitHistoPH(numberofbinsexpected,0,numberofbinsexpected);
	CreateFitHistoT0(numberofbinsexpected,0,numberofbinsexpected);
      }
    }
    
    
    //fDebug == 2 and fFitVoir no histo
    if(fDebug == 2){
      if(fFitVoir < numberofbinsexpected) {
	AliInfo(Form("We will see the fit of the object %d",fFitVoir));
      }
      else {
	AliInfo("fFitVoir is out of range of the histo!");
	return kFALSE;
      }
    }
    
    
    
    //fDebug == 3  or 4 and fDet
    if(fDebug >= 3){
      if((fNz[1] == 0) && (fNrphi[1] == 0)) {
	AliInfo("Do you really want to see one detector without pad groups?");
	return kFALSE;
      }
      else  {
	AliInfo(Form("You will see the detector: iPlane %d, iChamb %d, iSect %d",fDet[0], fDet[1], fDet[2]));
	
	//A little geometry:
	Int_t rowMax = parCom->GetRowMax(fDet[0],fDet[1],fDet[2]);
	Int_t colMax = parCom->GetColMax(fDet[0]);
	
	//Create the histos to visualise
	CreateFitHistoPHDB(rowMax, colMax);
	CreateFitHistoT0DB(rowMax, colMax);
	if(fDebug == 4) {
	  CreateFitHistoPH((Int_t) (fdect2[1]-fdect1[1]), fdect1[1], fdect2[1]);
	  CreateFitHistoT0((Int_t) (fdect2[1]-fdect1[1]), fdect1[1], fdect2[1]);
	}
	
      }
    }
  }
  //PRF
  if(i == 2) {
    
    //Number of Xbins(detectors or groups of pads) if vector2d*******************************
    if(Nbins == 0){
      //Quick verification that we are not out of range!
      if((Int_t) fVectorPRF.size() > numberofbinsexpected){
	AliInfo("ch doesn't correspond to the mode of pad group calibration!");
	return kFALSE;
      }
      if((Int_t) fVectorPRF.size() != (Int_t) fPlaPRF.size()){
	AliInfo("VectorPRF doesn't correspond to PlaCH!");
	return kFALSE;
      }
    }
    
    
    //init tree**********************************************
    InittreePRF();
    
    
    
    //Debugging____Create the histos********************************************
    
    //fDebug == 0 nothing
    
    //fDebug == 1 
    if(fDebug == 1){
      if(Nbins != 0){
	//create the histos replique de ch
	CreateFitHistoPRF(Nbins,lowedge,upedge);
      }
      else {
	//create the histos replique de ch
	CreateFitHistoPRF(numberofbinsexpected,0,numberofbinsexpected);
      }
    }
    
    
    //fDebug == 2 and fFitVoir no histo
    if(fDebug == 2){
      if(fFitVoir < numberofbinsexpected) {
	AliInfo(Form("We will see the fit of the object %d",fFitVoir));
      }
      else {
	AliInfo("fFitVoir is out of range of the histo!");
	return kFALSE;
      }
    }
    
    
    
    //fDebug == 3  or 4 and fDet
    if(fDebug >= 3){
      if((fNz[2] == 0) && (fNrphi[2] == 0)) {
	AliInfo("Do you really want to see one detector without pad groups?");
	return kFALSE;
      }
      else  {
	AliInfo(Form("You will see the detector: iPlane %d, iChamb %d, iSect %d",fDet[0], fDet[1], fDet[2]));
	
	//A little geometry:
	Int_t rowMax = parCom->GetRowMax(fDet[0],fDet[1],fDet[2]);
	Int_t colMax = parCom->GetColMax(fDet[0]);
	
	//Create the histos to visualise
	CreateFitHistoPRFDB(rowMax, colMax);
	if(fDebug == 4) CreateFitHistoPRF((Int_t) (fdect2[2]-fdect1[2]), fdect1[2], fdect2[2]);
	
      }
    }
  }

  return kTRUE;
  
}

//____________Functions for initialising the AliTRDCalibra in the code_________________________________________________________________________

void AliTRDCalibra::Initfcountdetandfcount(Int_t i)
{
  //
  //Init the current detector where we are fcountdet and the next fcount for the functions Fit... 
  //
  
  //loop on the Xbins of ch!!**********************************************************************
  fcountdet[i] = -1;//current detector
  fcount[i] = 0;//to find the next detector
  
  //if fDebug >= 3*******************************************************************
  if(fDebug >= 3){
    
    //Set countdet to the detector
    fcountdet[i] = AliTRDgeometry::GetDetector(fDet[0],fDet[1], fDet[2]);
    
    
    //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi
    ModePadCalibration(fDet[1], i);
    ModePadFragmentation(fDet[0],fDet[1],fDet[2],i);
    
    
    //Set counter to write at the end of the detector
    fcount[i] = fdect1[i] + fNfragz[i]*fNfragrphi[i];
    
    
  }
}

//____________Functions for initialising the AliTRDCalibra in the code_________________________________________________________________________

void AliTRDCalibra::Updatefcountdetandfcount(Int_t idect, Int_t i)
{
  //
  //See if we are in a new detector and update the variables fNfragz and fNfragrphi if yes 
  //

  
  //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi***********
  //if fDebug == 1 or 0
  if((fDebug == 0) || (fDebug == 1)){
    if(fcount[i] == idect) {
      
      //on en est au detector
      fcountdet[i] += 1;
      
      //Determination of fNnz, fNnrphi, fNfragz and fNfragrphi
      ModePadCalibration((Int_t) GetChamber(fcountdet[i]), i);
      ModePadFragmentation((Int_t) GetPlane(fcountdet[i]),(Int_t) GetChamber(fcountdet[i]),(Int_t) GetSector(fcountdet[i]),i);
      
      
      //Set for the next detector
      fcount[i] += fNfragz[i]*fNfragrphi[i];
      
    }
  }
  
}
//____________Functions for initialising the AliTRDCalibra in the code_________________________________________________________________________

void AliTRDCalibra::Reconstructfitrowminrowmax(Int_t idect, Int_t i)
{
  //
  // Reconstruct the min pad row, max pad row, min pad col and max pad col of the calibration group for the Fit functions
  //
  
  if(fDebug < 2) {
    ReconstructionRowPadGroup((Int_t) (idect-(fcount[i]-(fNfragz[i]*fNfragrphi[i]))), i);
  }
  if(fDebug >= 3){
    ReconstructionRowPadGroup((Int_t) (idect - fdect1[i]), i);
  }
}


//____________Functions for initialising the AliTRDCalibra in the code_________________________________________________________________________

Bool_t AliTRDCalibra::NotEnoughStatistic(Int_t idect, Int_t i)
{
  //
  // For the case where there are not enough entries in the histograms of the calibration group, the value present in the choosen database will be put. A negativ sign enables to know that a fit was not possible.
  //
  
  
  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }
  
  //Get cal
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }


  if(fDebug != 2) {
    AliInfo(Form("The element %d in this detector %d has not enough statistic to be fitted",idect-(fcount[i]-(fNfragz[i]*fNfragrphi[i])),fcountdet[i]));
  }
  if(fDebug == 2) {
    AliInfo("The element has not enough statistic to be fitted");
  }

  if( (i == 0) && (fDebug != 2)){
    
    //Calcul the coef from the database choosen
    CalculChargeCoefMean(fcountdet[0],(Int_t) (idect-fdect1[0]), kFALSE);
    
    //Fill the coefCH[2304] with negative value to say: not fitted
    for(Int_t k = frowmin[0]; k < frowmax[0]; k++){
      for(Int_t j = fcolmin[0]; j < fcolmax[0]; j++){
	if(GetChamber(fcountdet[0]) == 2) fcoefCH[(Int_t)(j*12+k)]=-TMath::Abs(fChargeCoef[3]);
	if(GetChamber(fcountdet[0]) != 2) fcoefCH[(Int_t)(j*16+k)]=-TMath::Abs(fChargeCoef[3]);
      }
    }                
    //end of one detector
    if((idect == (fcount[0]-1))) {
      FillVectorFitCH((Int_t) fcountdet[0]);
      //Reset
      for(Int_t k = 0; k < 2304; k++){
	fcoefCH[k] = 0.0;    
      }
    }
  }
  
  if((i == 1) && (fDebug != 2)){

    CalculVdriftCoefMean(fcountdet[1],(Int_t) (idect-fdect1[1]));
    CalculT0CoefMean(fcountdet[1],(Int_t) (idect-fdect1[1]));
    
    
    //Put the default value 
    if((fDebug == 1) || (fDebug == 4)) {
      if(fFitPHOn) {
	fCoefVdrift[0]->SetBinContent(idect-fdect1[1]+1,fVdriftCoef[2]);
	fCoefT0[0]->SetBinContent(idect-fdect1[1]+1,fT0Coef[2]);
      }
      fCoefVdrift[1]->SetBinContent(idect-fdect1[1]+1,fVdriftCoef[2]);
      fCoefT0[1]->SetBinContent(idect-fdect1[1]+1,fT0Coef[2]);
    }
    
    //Put the default value
    if(fDebug >= 3){
      fVdriftCoef[0] = fVdriftCoef[2];
      fVdriftCoef[1] = fVdriftCoef[2];
      FillCoefVdriftDB();
      fT0Coef[0] = fT0Coef[2];
      fT0Coef[1] = fT0Coef[2];
      FillCoefT0DB();
    }
   

    //Fill the tree if end of a detector. The pointer to the branch stays with the default value 1.5!!!*******
    //PH*******************************
      //pointer to the branch
    for(Int_t k = frowmin[1]; k < frowmax[1]; k++){
      for(Int_t j = fcolmin[1]; j < fcolmax[1]; j++){
	if(GetChamber(fcountdet[1]) == 2) fVdriftPad[(Int_t)(j*12+k)]=-TMath::Abs(fVdriftCoef[2]);
	if(GetChamber(fcountdet[1]) != 2) fVdriftPad[(Int_t)(j*16+k)]=-TMath::Abs(fVdriftCoef[2]);
      }
    }     
    //end of one detector
    if((idect == (fcount[1]-1)) && (fDebug != 2)) FilltreeVdrift((Int_t) fcountdet[1]);
    
    
    //T0************************************
    //Fill the tree if end of a detector. The pointer to the branch stays with the default value 1.5!!!*******
    //pointer to the branch
    for(Int_t k = frowmin[1]; k < frowmax[1]; k++){
      for(Int_t j = fcolmin[1]; j < fcolmax[1]; j++){
	if(GetChamber(fcountdet[1]) == 2) fT0Pad[(Int_t)(j*12+k)]=-TMath::Abs(fT0Coef[2]);
	if(GetChamber(fcountdet[1]) != 2) fT0Pad[(Int_t)(j*16+k)]=-TMath::Abs(fT0Coef[2]);
      }
    }     
    //end of one detector
    if((idect == (fcount[1]-1)) && (fDebug != 2)) FilltreeT0((Int_t) fcountdet[1]);
  }
    
  
  if((i == 2) && (fDebug != 2)){

    CalculPRFCoefMean(fcountdet[2],(Int_t) (idect-fdect1[2]));
      
    if((fDebug == 1) || (fDebug == 4)){
      fCoefPRF[0]->SetBinContent(idect-fdect1[2]+1,fPRFCoef[1]);
    }
    if(fDebug >= 3){
      fPRFCoef[0] = fPRFCoef[1];
      FillCoefPRFDB();
    }

    //Fill the tree if end of a detector. The pointer to the branch stays with the default value 1.5!!!*******
    //pointer to the branch
    for(Int_t k = frowmin[2]; k < frowmax[2]; k++){
      for(Int_t j = fcolmin[2]; j < fcolmax[2]; j++){
	if((parCom->GetColMax(GetPlane(fcountdet[2])) != (j+1)) && (j != 0)){
	  if(GetChamber(fcountdet[2]) == 2) fPRFPad[(Int_t)(j*12+k)]=-fPRFCoef[1];
	  if(GetChamber(fcountdet[2]) != 2) fPRFPad[(Int_t)(j*16+k)]=-fPRFCoef[1];
	}
	else{
	  if(GetChamber(fcountdet[2]) == 2) fPRFPad[(Int_t)(j*12+k)]=-((Float_t)cal->GetPRFWidth(fcountdet[2],j,k));
	  if(GetChamber(fcountdet[2]) != 2) fPRFPad[(Int_t)(j*16+k)]=-((Float_t) cal->GetPRFWidth(fcountdet[2],j,k));
	}
      }
    }      
    //end of one detector
    if((idect == (fcount[2]-1)) && (fDebug != 2)) FilltreePRF((Int_t) fcountdet[2]);
  }
  
  return kTRUE;
  
}

//____________Functions for initialising the AliTRDCalibra in the code_________________________________________________________________________
Bool_t AliTRDCalibra::FillInfosFit(Int_t idect, Int_t i)
{
  //
  // Fill the coefficients found with the fits or other methods from the Fit functions
  //
  



  // Get the parameter object
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam Manager");
    return kFALSE;
  }

  //Get cal
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB");
    return kFALSE;
  }

  if((i == 0) && (fDebug != 2)){
    //Fill the coefCH[2304] with fChargeCoef[0] that would be negativ only if the fit failed totally*******************************
    for(Int_t k = frowmin[0]; k < frowmax[0]; k++){
      for(Int_t j = fcolmin[0]; j < fcolmax[0]; j++){
	if(GetChamber(fcountdet[0]) == 2) fcoefCH[(Int_t)(j*12+k)]= fChargeCoef[0];
	if(GetChamber(fcountdet[0]) != 2) fcoefCH[(Int_t)(j*16+k)]= fChargeCoef[0];
      }
    }                
    //end of one detector
    if((idect == (fcount[0]-1))) {
      FillVectorFitCH((Int_t) fcountdet[0]);
      //Reset
      for(Int_t k = 0; k < 2304; k++){
	fcoefCH[k] = 0.0;    
      }
    }
  }
  if((i == 1) && (fDebug != 2)){

    //PH*****************************************************************
    //pointer to the branch: fVdriftCoef[1] will ne negativ only if the fit failed totally 
    for(Int_t k = frowmin[1]; k < frowmax[1]; k++){
      for(Int_t j = fcolmin[1]; j < fcolmax[1]; j++){
	if(GetChamber(fcountdet[1]) == 2) fVdriftPad[(Int_t)(j*12+k)]=fVdriftCoef[1];
	if(GetChamber(fcountdet[1]) != 2) fVdriftPad[(Int_t)(j*16+k)]=fVdriftCoef[1];
      }
    }                
    //end of one detector
    if((idect == (fcount[1]-1)) && (fDebug != 2)) FilltreeVdrift((Int_t) fcountdet[1]);




    //T0*****************************************************************8
    //pointer to the branch: fT0Coef[1] will ne negativ only if the fit failed totally 
    for(Int_t k = frowmin[1]; k < frowmax[1]; k++){
      for(Int_t j = fcolmin[1]; j < fcolmax[1]; j++){
	if(GetChamber(fcountdet[1]) == 2) fT0Pad[(Int_t)(j*12+k)]=fT0Coef[1];
	if(GetChamber(fcountdet[1]) != 2) fT0Pad[(Int_t)(j*16+k)]=fT0Coef[1];
      }
    }                
    //end of one detector
    if((idect == (fcount[1]-1)) && (fDebug != 2)) FilltreeT0((Int_t) fcountdet[1]);
    
  }
  if((i == 2) && (fDebug != 2)){
    //pointer to the branch
    for(Int_t k = frowmin[2]; k < frowmax[2]; k++){
      for(Int_t j = fcolmin[2]; j < fcolmax[2]; j++){
	if((parCom->GetColMax(GetPlane(fcountdet[2])) != (j+1)) && (j != 0) && (parCom->GetRowMax(GetPlane(fcountdet[2]), GetChamber(fcountdet[2]),GetSector(fcountdet[2])) != (k+1)) && (k != 0)){
	    if(GetChamber(fcountdet[2]) == 2) fPRFPad[(Int_t)(j*12+k)]=fPRFCoef[0];
	    if(GetChamber(fcountdet[2]) != 2) fPRFPad[(Int_t)(j*16+k)]=fPRFCoef[0];
	    
	  }
	  else{
	    if(GetChamber(fcountdet[2]) == 2) fPRFPad[(Int_t)(j*12+k)]=(Float_t) cal->GetPRFWidth(fcountdet[2],j,k);
	    if(GetChamber(fcountdet[2]) != 2) fPRFPad[(Int_t)(j*16+k)]=(Float_t) cal->GetPRFWidth(fcountdet[2],j,k);
          }
	}
      }                
      //end of one detector
      if((idect == (fcount[2]-1)) && (fDebug != 2)) FilltreePRF((Int_t) fcountdet[2]);
  }
  return kTRUE;

}
//____________Functions for initialising the AliTRDCalibra in the code_________________________________________________________________________

Bool_t AliTRDCalibra::WriteFitInfos(Int_t i)
{
  //
  // In the case the user wants to write a file with a tree of the found coefficients for the calibration before putting them in the database
  //



  TFile *fout = TFile::Open(fWriteNameCoef,"UPDATE");
  //Check if the file could be opened
  if (!fout || !fout->IsOpen()) {
    AliInfo("No File found!");
    return kFALSE;
  }
  if((i == 0) && (fDebug != 2)){
    //The error stuff  
    if((fDebug == 1) || (fDebug == 4)) WriteCH(fout);
    //The DB stuff
    if((fDebug == 4) || (fDebug == 3)) WriteCHDB(fout);
    //The tree
    fout->WriteTObject(fGain, fGain->GetName(), (Option_t *) "writedelete");
  }
  if((i == 1) && (fDebug != 2)){
    //PH********************************************
    //The error stuff  
    if((fDebug == 1) || (fDebug == 4)) WritePH(fout);
    //The DB stuff
    if((fDebug == 4) || (fDebug == 3)) WritePHDB(fout);
    //The tree
    fout->WriteTObject(fVdrift, fVdrift->GetName(), (Option_t *) "writedelete");
    //T0********************************************
    //The error stuff  
    if((fDebug == 1) || (fDebug == 4)) WriteT0(fout);
    //The DB stuff
    if((fDebug == 4) || (fDebug == 3)) WriteT0DB(fout);
    //The tree
    fout->WriteTObject(fT0, fT0->GetName(), (Option_t *) "writedelete");
  }
  if((i == 2) && (fDebug != 2)){
    //The error stuff  
    if((fDebug == 1) || (fDebug == 4)) WritePRF(fout);
    //The DB stuff
    if((fDebug == 4) || (fDebug == 3)) WritePRFDB(fout);
    //The tree
    fout->WriteTObject(fPRF, fPRF->GetName(), (Option_t *) "writedelete");
  }
  fout->Close();
  return kTRUE;
}
//____Fill the Error histos in case of fDebug == 1________________________________________________________________________________

//__________________________________________________________________________________________
void AliTRDCalibra::ErrorPRF()
{
  //
  // Fill the error histos for fDebug = 1 and fDebug = 4 from the delta histos
  //

  for(Int_t k= 0; k < fDeltaPRF->GetNbinsX(); k++){
    if(fDeltaPRF->GetBinContent(k+1) != 0.0) fErrorPRF->Fill(fDeltaPRF->GetBinContent(k+1));
  }
 
}
//__________________________________________________________________________________________
void AliTRDCalibra::ErrorCH()
{
  //
  // Fill the error histos for fDebug = 1 and fDebug = 4 from the delta histos
  //
  for(Int_t k= 0; k < fDeltaCharge[0]->GetNbinsX(); k++){
    if(fDeltaCharge[0]->GetBinContent(k+1) != 0.0) fErrorCharge[0]->Fill(fDeltaCharge[0]->GetBinContent(k+1));
  }
  if(fMeanChargeOn){
    for(Int_t k= 0; k < fDeltaCharge[1]->GetNbinsX(); k++){
      if(fDeltaCharge[1]->GetBinContent(k+1) != 0.0) fErrorCharge[1]->Fill(fDeltaCharge[1]->GetBinContent(k+1));
    }
  }
  if(fFitChargeBisOn ){
    for(Int_t k= 0; k < fDeltaCharge[2]->GetNbinsX(); k++){
      if(fDeltaCharge[2]->GetBinContent(k+1) != 0.0) fErrorCharge[2]->Fill(fDeltaCharge[2]->GetBinContent(k+1));
    }
  }
 
}
//______________________________________________________________________________________________
void AliTRDCalibra::ErrorPH()
{
  //
  // Fill the error histos for fDebug = 1 and fDebug = 4 from the delta histos
  //
  if(fFitPHOn ){
    for(Int_t k= 0; k < fDeltaVdrift[0]->GetNbinsX(); k++){     
      if(fDeltaVdrift[0]->GetBinContent(k+1) != 0.0) fErrorVdrift[0]->Fill(fDeltaVdrift[0]->GetBinContent(k+1));
    }
  }
  for(Int_t k= 0; k < fDeltaVdrift[1]->GetNbinsX(); k++){
    if(fDeltaVdrift[1]->GetBinContent(k+1) != 0.0) fErrorVdrift[1]->Fill(fDeltaVdrift[1]->GetBinContent(k+1));
  }
}
//______________________________________________________________________________________________
void AliTRDCalibra::ErrorT0()
{
  //
  // Fill the error histos for fDebug = 1 and fDebug = 4 from the delta histos
  //
  if(fFitPHOn ){
    for(Int_t k= 0; k < fDeltaT0[0]->GetNbinsX(); k++){     
      if(fDeltaT0[0]->GetBinContent(k+1) != 0.0) fErrorT0[0]->Fill(fDeltaT0[0]->GetBinContent(k+1));
    }
  }
  for(Int_t k= 0; k < fDeltaT0[1]->GetNbinsX(); k++){
    if(fDeltaT0[1]->GetBinContent(k+1) != 0.0) fErrorT0[1]->Fill(fDeltaT0[1]->GetBinContent(k+1));
  }
}
//_____Fill Coef DB in case of visualisation of one detector________________________________________________________________________

//________________________________________________________________________________________
void AliTRDCalibra::FillCoefVdriftDB()
{
  //
  // Fill the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //
 
  for(Int_t row = frowmin[1]; row < frowmax[1]; row++){
    for(Int_t col = fcolmin[1]; col < fcolmax[1]; col++){
      fCoefVdriftDB[1]->SetBinContent(row+1,col+1,TMath::Abs(fVdriftCoef[1]));
      if(fFitPHOn ) fCoefVdriftDB[0]->SetBinContent(row+1,col+1,TMath::Abs(fVdriftCoef[0]));
    }
  }
}
//________________________________________________________________________________________
void AliTRDCalibra::FillCoefT0DB()
{
  //
  // Fill the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //
 
  for(Int_t row = frowmin[1]; row < frowmax[1]; row++){
    for(Int_t col = fcolmin[1]; col < fcolmax[1]; col++){
      fCoefT0DB[1]->SetBinContent(row+1,col+1,TMath::Abs(fT0Coef[1]));
      if(fFitPHOn ) fCoefT0DB[0]->SetBinContent(row+1,col+1,TMath::Abs(fT0Coef[0]));
    }
  }
}
//__________________________________________________________________________________________
void AliTRDCalibra::FillCoefChargeDB()
{
  //
  // Fill the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  for(Int_t row = frowmin[0]; row < frowmax[0]; row++){
    for(Int_t col = fcolmin[0]; col < fcolmax[0]; col++){
      if(fMeanChargeOn) fCoefChargeDB[1]->SetBinContent(row+1,col+1,TMath::Abs(fChargeCoef[1]));
      if(fFitChargeBisOn ) fCoefChargeDB[2]->SetBinContent(row+1,col+1,TMath::Abs(fChargeCoef[2]));
      fCoefChargeDB[0]->SetBinContent(row+1,col+1,TMath::Abs(fChargeCoef[0]));
    }
  }   
}
//__________________________________________________________________________________________
void AliTRDCalibra::FillCoefPRFDB()
{
  //
  // Fill the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //
  for(Int_t row = frowmin[2]; row < frowmax[2]; row++){
    for(Int_t col = fcolmin[2]; col < fcolmax[2]; col++){
      fCoefPRFDB->SetBinContent(row+1,col+1,fPRFCoef[0]);
     
    }
  }   
}

//_________Plot histos CoefPRF....__________________________________________________________________________________________


//___________________________________________________________________________________________
void AliTRDCalibra::PlotCH()
{
  //
  // Plot the histos for fDebug = 1 and fDebug = 4 for the errors
  //
  TLatex t;
  t.SetTextAlign(22);
  t.SetTextSize(0.1);
  TCanvas *cch1 = new TCanvas("cch1","",50,50,600,800);
  cch1->cd();
  TLegend *legch1 = new TLegend(0.4,0.6,0.89,0.89);
  legch1->AddEntry(fCoefCharge[3],"f_{g} simulated","l");
  if(fMeanChargeOn) legch1->AddEntry(fCoefCharge[1],"f_{g} mean","p");
  legch1->AddEntry(fCoefCharge[0],"f_{g} fit","p");
  if(fFitChargeBisOn ) legch1->AddEntry(fCoefCharge[2],"f_{g} fitbis","p");
 
  fCoefCharge[0]->Draw("E2");
  if(fMeanChargeOn)  fCoefCharge[1]->Draw("E2 same");
  if(fFitChargeBisOn ) fCoefCharge[2]->Draw("E2 same");
  fCoefCharge[3]->Draw("same");
  legch1->Draw("same");
  
  
  
  TCanvas *cch2 = new TCanvas("cch2","",50,50,600,800);
  cch2->Divide(2,1);
  cch2->cd(1);
  TLegend *legch2 = new TLegend(0.4,0.6,0.89,0.89);
  if(fMeanChargeOn) legch2->AddEntry(fErrorCharge[1],"f_{g} mean","l");
  legch2->AddEntry(fErrorCharge[0],"f_{g} fit","l");
  if(fFitChargeBisOn ) legch2->AddEntry(fErrorCharge[2],"f_{g} fitbis","l");
  fErrorCharge[0]->Draw();
  if(fMeanChargeOn) fErrorCharge[1]->Draw("same");
  if(fFitChargeBisOn ) fErrorCharge[2]->Draw("same");
  legch2->Draw("same");
  cch2->cd(2);
  TLegend *legch3 = new TLegend(0.4,0.6,0.89,0.89);
  if(fMeanChargeOn) legch3->AddEntry(fDeltaCharge[1],"mean","p");
  legch3->AddEntry(fDeltaCharge[0],"fit","p");
  if(fFitChargeBisOn ) legch3->AddEntry(fDeltaCharge[2],"fit","p");
  fDeltaCharge[0]->Draw("E2");
  if(fMeanChargeOn) fDeltaCharge[1]->Draw("E2 same");
  if(fFitChargeBisOn ) fDeltaCharge[2]->Draw("E2 same");
  legch3->Draw("same");
  
}
//__________________________________________________________________________
void AliTRDCalibra::PlotPH()
{
  //
  // Plot the histos for fDebug = 1 and fDebug = 4 for the errors
  //

  TCanvas *cph1 = new TCanvas("cph1","",50,50,600,800);
  cph1->cd();
  TLegend *legph1 = new TLegend(0.4,0.6,0.89,0.89);
  legph1->AddEntry(fCoefVdrift[2],"v_{real} simulated","l");
  legph1->AddEntry(fCoefVdrift[1],"v_{sm} slope 1 method","p");
 
  if(fFitPHOn ) legph1->AddEntry(fCoefVdrift[0],"v_{fit} fit","p");
  fCoefVdrift[1]->Draw("E2");
  fCoefVdrift[2]->Draw("same");
  if(fFitPHOn )fCoefVdrift[0]->Draw("E2 same");
  legph1->Draw("same");
  

  

  TCanvas *cph2 = new TCanvas("cph2","",50,50,600,800);
  cph2->Divide(2,1);
  cph2->cd(1);
  TLegend *legph2 = new TLegend(0.4,0.6,0.89,0.89);
  legph2->AddEntry(fErrorVdrift[1],"v_{sm} slope 1 method","l");
  if(fFitPHOn ) legph2->AddEntry(fErrorVdrift[0],"v_{fit} fit","l");
  fErrorVdrift[1]->Draw();
  if(fFitPHOn ) fErrorVdrift[0]->Draw("l,same");
  legph2->Draw("same");
  cph2->cd(2);
  TLegend *legph3 = new TLegend(0.4,0.6,0.89,0.89);
  legph3->AddEntry(fDeltaVdrift[1],"v_{sm} slope 1 method","p");
  if(fFitPHOn ) legph3->AddEntry(fDeltaVdrift[0],"v_{fit} fit","p");
  fDeltaVdrift[1]->Draw("E2");
  if(fFitPHOn ) fDeltaVdrift[0]->Draw("E2 same");
  legph3->Draw("same");
}
//__________________________________________________________________________
void AliTRDCalibra::PlotT0()
{
  //
  // Plot the histos for fDebug = 1 and fDebug = 4 for the errors
  //

  TCanvas *ct01 = new TCanvas("ct01","",50,50,600,800);
  ct01->cd();
  TLegend *legt01 = new TLegend(0.4,0.6,0.89,0.89);
  legt01->AddEntry(fCoefT0[2],"t0 simulated","l");
  legt01->AddEntry(fCoefT0[1],"t0 slope 1 method","p");
 
  if(fFitPHOn ) legt01->AddEntry(fCoefT0[0],"t0 fit","p");
  fCoefT0[1]->Draw("E2");
  fCoefT0[2]->Draw("same");
  if(fFitPHOn )fCoefT0[0]->Draw("E2 same");
  legt01->Draw("same");
  

  TCanvas *ct02 = new TCanvas("ct02","",50,50,600,800);
  ct02->Divide(2,1);
  ct02->cd(1);
  TLegend *legt02 = new TLegend(0.4,0.6,0.89,0.89);
  legt02->AddEntry(fErrorT0[1],"t0 slope 1 method","l");
  if(fFitPHOn ) legt02->AddEntry(fErrorT0[0],"t0 fit","l");
  fErrorT0[1]->Draw();
  if(fFitPHOn ) fErrorT0[0]->Draw("l,same");
  legt02->Draw("same");
  ct02->cd(2);
  TLegend *legt03 = new TLegend(0.4,0.6,0.89,0.89);
  legt03->AddEntry(fDeltaT0[1],"t0 slope 1 method","p");
  if(fFitPHOn ) legt03->AddEntry(fDeltaT0[0],"t0 fit","p");
  fDeltaT0[1]->Draw("E2");
  if(fFitPHOn ) fDeltaT0[0]->Draw("E2 same");
  legt03->Draw("same");
}

//___________________________________________________________________________
void AliTRDCalibra::PlotPRF()
{
  //
  // Plot the histos for fDebug = 1 and fDebug = 4 for the errors
  //

  TCanvas *cprf1 = new TCanvas("cprf1","",50,50,600,800);
  cprf1->cd();
  TLegend *legprf1 = new TLegend(0.4,0.6,0.89,0.89);
  legprf1->AddEntry(fCoefPRF[1],"#sigma_{real} simulated","l");
  legprf1->AddEntry(fCoefPRF[0],"#sigma_{fit} reconstructed","p");

  fCoefPRF[0]->Draw("E2");
  fCoefPRF[1]->Draw("same");
  legprf1->Draw("same");
  

  

  TCanvas *cprf2 = new TCanvas("cprf2","",50,50,600,800);
  cprf2->Divide(2,1);
  cprf2->cd(1);
  TLegend *legprf2 = new TLegend(0.4,0.6,0.89,0.89);
  legprf2->AddEntry(fErrorPRF,"#sigma_{fit} reconstructed","l");
  fErrorPRF->Draw("");
  legprf2->Draw("same");
  cprf2->cd(2);
  TLegend *legprf3 = new TLegend(0.4,0.6,0.89,0.89);
  legprf3->AddEntry(fDeltaPRF,"#sigma_{fit} reconstructed","p");
  fDeltaPRF->Draw("E2");
  legprf3->Draw("same");

}

//___________Plot histos DB ___________________________________________________________________________________________________
//___________________________________________________________________________________________
void AliTRDCalibra::PlotCHDB()
{
  //
  // Plot the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  TCanvas *cchdb = new TCanvas("cchdb","",50,50,600,800);
  if((fFitChargeBisOn) && (fMeanChargeOn)){
    cchdb->Divide(3,1);
    cchdb->cd(1);
    fCoefChargeDB[0]->Draw("LEGO");
    cchdb->cd(2);
    fCoefChargeDB[1]->Draw("LEGO");
    cchdb->cd(3);
    fCoefChargeDB[2]->Draw("LEGO");
  }
  if((!fFitChargeBisOn) && (fMeanChargeOn)){
    cchdb->Divide(2,1);
    cchdb->cd(1);
    fCoefChargeDB[0]->Draw("LEGO");
    cchdb->cd(2);
    fCoefChargeDB[1]->Draw("LEGO");
  }
  else{
    cchdb->cd();
    fCoefChargeDB[0]->Draw("LEGO");
  }

}
//_______________________________________________________________________________________

void AliTRDCalibra::PlotPHDB()
{
  //
  // Plot the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  TCanvas *cphdb = new TCanvas("cphdb","",50,50,600,800);
  if(fFitPHOn ){
    cphdb->Divide(2,1);
    cphdb->cd(1);
    fCoefVdriftDB[0]->Draw("LEGO");
    cphdb->cd(2);
    fCoefVdriftDB[1]->Draw("LEGO");
  }
  else{
    cphdb->cd();
    fCoefVdriftDB[1]->Draw("LEGO");
  }
}
//_______________________________________________________________________________________

void AliTRDCalibra::PlotT0DB()
{
  //
  // Plot the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  TCanvas *ct0db = new TCanvas("ct0db","",50,50,600,800);
  if(fFitPHOn ){
    ct0db->Divide(2,1);
    ct0db->cd(1);
    fCoefT0DB[0]->Draw("LEGO");
    ct0db->cd(2);
    fCoefT0DB[1]->Draw("LEGO");
  }
  else{
    ct0db->cd();
    fCoefT0DB[1]->Draw("LEGO");
  }
}
//_______________________________________________________________________________________

void AliTRDCalibra::PlotPRFDB()
{
  //
  // Plot the histos for fDebug = 3 and fDebug = 4 to visualise the detector
  //

  TCanvas *cprfdb = new TCanvas("cprfdb","",50,50,600,800);
  cprfdb->cd();
  fCoefPRFDB->Draw("LEGO");
}

//______________Write histos Coef________________________________________________________________________________________________________

//________________________________________________________________________________
void AliTRDCalibra::WriteCH(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 1 and fDebug = 4
  //

  fout->WriteTObject(fCoefCharge[0], fCoefCharge[0]->GetName(),(Option_t *)"OverWrite");
  if(fMeanChargeOn) fout->WriteTObject(fCoefCharge[1], fCoefCharge[1]->GetName(),(Option_t *)"OverWrite");
  if(fFitChargeBisOn ) fout->WriteTObject(fCoefCharge[2], fCoefCharge[2]->GetName(),(Option_t *)"OverWrite");
  
  fout->WriteTObject(fCoefCharge[3], fCoefCharge[3]->GetName(),(Option_t *)"OverWrite");
  
  fout->WriteTObject(fDeltaCharge[0], fDeltaCharge[0]->GetName(),(Option_t *)"OverWrite");
  if(fMeanChargeOn) fout->WriteTObject(fDeltaCharge[1], fDeltaCharge[1]->GetName(),(Option_t *)"OverWrite");
  if(fFitChargeBisOn ) fout->WriteTObject(fDeltaCharge[2], fDeltaCharge[2]->GetName(),(Option_t *)"OverWrite");
  
  
  fout->WriteTObject(fErrorCharge[0], fErrorCharge[0]->GetName(),(Option_t *)"OverWrite");
  if(fMeanChargeOn) fout->WriteTObject(fErrorCharge[1], fErrorCharge[1]->GetName(),(Option_t *)"OverWrite");
  if(fFitChargeBisOn ) fout->WriteTObject(fErrorCharge[2], fErrorCharge[2]->GetName(),(Option_t *)"OverWrite");
  

}
//________________________________________________________________________________
void AliTRDCalibra::WritePH(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 1 and fDebug = 4
  //

  if (fFitPHOn )fout->WriteTObject(fCoefVdrift[0], fCoefVdrift[0]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fCoefVdrift[1], fCoefVdrift[1]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fCoefVdrift[2], fCoefVdrift[2]->GetName(),(Option_t *)"OverWrite");

    
  if(fFitPHOn ) fout->WriteTObject(fDeltaVdrift[0], fDeltaVdrift[0]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fDeltaVdrift[1], fDeltaVdrift[1]->GetName(),(Option_t *)"OverWrite");
 
    
  if(fFitPHOn ) fout->WriteTObject(fErrorVdrift[0], fErrorVdrift[0]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fErrorVdrift[1], fErrorVdrift[1]->GetName(),(Option_t *)"OverWrite");
 
  
}
//________________________________________________________________________________
void AliTRDCalibra::WriteT0(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 1 and fDebug = 4
  //

  if (fFitPHOn )fout->WriteTObject(fCoefT0[0], fCoefT0[0]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fCoefT0[1], fCoefT0[1]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fCoefT0[2], fCoefT0[2]->GetName(),(Option_t *)"OverWrite");

    
  if(fFitPHOn ) fout->WriteTObject(fDeltaT0[0], fDeltaT0[0]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fDeltaT0[1], fDeltaT0[1]->GetName(),(Option_t *)"OverWrite");
 
    
  if(fFitPHOn ) fout->WriteTObject(fErrorT0[0], fErrorT0[0]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fErrorT0[1], fErrorT0[1]->GetName(),(Option_t *)"OverWrite");
 
  
}
//________________________________________________________________________________
void AliTRDCalibra::WritePRF(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 1 and fDebug = 4
  //
  
  fout->WriteTObject(fCoefPRF[0], fCoefPRF[0]->GetName(),(Option_t *)"OverWrite");
     
  fout->WriteTObject(fCoefPRF[1], fCoefPRF[1]->GetName(),(Option_t *)"OverWrite");

  fout->WriteTObject(fDeltaPRF, fDeltaPRF->GetName(),(Option_t *)"OverWrite");
        
  fout->WriteTObject(fErrorPRF, fErrorPRF->GetName(),(Option_t *)"OverWrite");
}


//_________________________Write DB Histos____________________________________________________________________________________


//________________________________________________________________________________
void AliTRDCalibra::WriteCHDB(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 3 and fDebug = 4
  //

  fout->WriteTObject(fCoefChargeDB[0], fCoefChargeDB[0]->GetName(),(Option_t *)"OverWrite");
  if(fMeanChargeOn) fout->WriteTObject(fCoefChargeDB[1], fCoefChargeDB[1]->GetName(),(Option_t *)"OverWrite");
  if(fFitChargeBisOn ) fout->WriteTObject(fCoefChargeDB[2], fCoefChargeDB[2]->GetName(),(Option_t *)"OverWrite");
  

}
//________________________________________________________________________________
void AliTRDCalibra::WritePHDB(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 3 and fDebug = 4
  //

  if (fFitPHOn )fout->WriteTObject(fCoefVdriftDB[0], fCoefVdriftDB[0]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fCoefVdriftDB[1], fCoefVdriftDB[1]->GetName(),(Option_t *)"OverWrite");



}
//________________________________________________________________________________
void AliTRDCalibra::WriteT0DB(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 3 and fDebug = 4
  //

  if (fFitPHOn )fout->WriteTObject(fCoefT0DB[0], fCoefT0DB[0]->GetName(),(Option_t *)"OverWrite");
  fout->WriteTObject(fCoefT0DB[1], fCoefT0DB[1]->GetName(),(Option_t *)"OverWrite");



}
//________________________________________________________________________________
void AliTRDCalibra::WritePRFDB(TFile *fout)
{
  //
  // If wanted, write the debug histos for fDebug = 3 and fDebug = 4
  //

  fout->WriteTObject(fCoefPRFDB, fCoefPRFDB->GetName(),(Option_t *)"OverWrite");

}


//_______________Calcul Coef Mean__________________________________________________________________________________________________


//__________________________________________________________________________________________
Bool_t AliTRDCalibra::CalculT0CoefMean(Int_t Dect, Int_t idect)
{
  //
  // For the detector Dect calcul the mean time 0 for the calibration group idect from the choosen database
  //

  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB  Manager");
    return kFALSE;
  }

  fT0Coef[2] = 0.0;

  if(fDebug != 2){

    for(Int_t row = frowmin[1]; row < frowmax[1]; row++){
      for(Int_t col = fcolmin[1]; col< fcolmax[1]; col++){
	//groups of pads
	if((fNz[1] > 0) && (fNrphi[1] > 0)) {
	  fT0Coef[2] += (Float_t) cal->GetT0(Dect,col,row);
	}
	//per detectors
	else fT0Coef[2] += (Float_t) cal->GetT0Average(Dect);
      }
    }
    fT0Coef[2] = fT0Coef[2]/((fcolmax[1]-fcolmin[1])*(frowmax[1]-frowmin[1]));
    if((fDebug == 1) || (fDebug == 4)) fCoefT0[2]->SetBinContent(idect+1,fT0Coef[2]);
  }
  return kTRUE;
  
}
//__________________________________________________________________________________________
Bool_t AliTRDCalibra::CalculChargeCoefMean(Int_t Dect, Int_t idect, Bool_t vrai)
{
  //
  // For the detector Dect calcul the mean gain factor for the calibration group idect from the choosen database
  //


  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB  Manager");
    return kFALSE;
  }
  AliTRDCommonParam     *parCom    = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam  Manager");
    return kFALSE;
  }

  fChargeCoef[3] = 0.0;
  if(fDebug != 2){

    for(Int_t row = frowmin[0]; row < frowmax[0]; row++){
      for(Int_t col = fcolmin[0]; col< fcolmax[0]; col++){
	//groups of pads
	if((fNz[0] > 0) || (fNrphi[0] > 0)) {
	  fChargeCoef[3] += (Float_t) cal->GetGainFactor(Dect,col,row);
	  if(vrai) fScalefitfactor += (Float_t) cal->GetGainFactor(Dect,col,row);
	}
	//per detectors
	else {
	  fChargeCoef[3] += (Float_t) cal->GetGainFactorAverage(Dect);
	  if(vrai) fScalefitfactor += ((Float_t) cal->GetGainFactorAverage(Dect));
	}
      }
    }
    fChargeCoef[3] = fChargeCoef[3]/((fcolmax[0]-fcolmin[0])*(frowmax[0]-frowmin[0]));
    if((fDebug == 1) || (fDebug == 4)) fCoefCharge[3]->SetBinContent(idect+1,fChargeCoef[3]);
  }
    return kTRUE;
    


}
//__________________________________________________________________________________________
Bool_t AliTRDCalibra::CalculPRFCoefMean(Int_t Dect, Int_t idect)
{
  //
  // For the detector Dect calcul the mean sigma of pad response function for the calibration group idect from the choosen database
  //


  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB  Manager");
    return kFALSE;
  }

  AliTRDCommonParam     *parCom    = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam  Manager");
    return kFALSE;
  }


  fPRFCoef[1] = 0.0;
  Int_t cot = 0;
  if(fDebug != 2){
    
    for(Int_t row = frowmin[2]; row < frowmax[2]; row++){
      for(Int_t col = fcolmin[2]; col< fcolmax[2]; col++){
	if((parCom->GetColMax(GetPlane(Dect)) != (col+1)) && (col != 0)){
	  cot++;
	  fPRFCoef[1] += (Float_t) cal->GetPRFWidth(Dect,col,row);
	}
      }
    }
    if(cot > 0){
      fPRFCoef[1] = fPRFCoef[1]/cot;
      if((fDebug == 1) ||(fDebug == 4)) fCoefPRF[1]->SetBinContent(idect+1,fPRFCoef[1]);
    }
    else{
      if((fDebug == 1) ||(fDebug == 4)) fCoefPRF[1]->SetBinContent(idect+1, cal->GetPRFWidth(Dect,fcolmin[2],frowmin[2]));
    }
  }
  return kTRUE;
  
}
//__________________________________________________________________________________________
Bool_t AliTRDCalibra::CalculVdriftCoefMean(Int_t Dect, Int_t idect)
{
  //
  // For the detector Dect calcul the mean drift velocity for the calibration group idect from the choosen database
  //


  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();
  if (!cal) {
    AliInfo("Could not get calibDB  Manager");
    return kFALSE;
  }


  fVdriftCoef[2] = 0.0;

  if(fDebug != 2){
    for(Int_t row = frowmin[1]; row < frowmax[1]; row++){
      for(Int_t col = fcolmin[1]; col< fcolmax[1]; col++){
	//groups of pads
	if((fNz[1] > 0) || (fNrphi[1] > 0)) {
	  fVdriftCoef[2] += (Float_t) cal->GetVdrift(Dect,col,row);
	}
	//per detectors
	else fVdriftCoef[2] += (Float_t) cal->GetVdriftAverage(Dect);
      }
    }
    fVdriftCoef[2] = fVdriftCoef[2]/((fcolmax[1]-fcolmin[1])*(frowmax[1]-frowmin[1]));
    if((fDebug == 1) || (fDebug == 4)) fCoefVdrift[2]->SetBinContent(idect+1,fVdriftCoef[2]);
  }
  return kTRUE;
  
}

//_________Pad group calibration mode____________________________________________________________________________________________________________
//________________________________________________________________________________________
void AliTRDCalibra::ReconstructionRowPadGroup(Int_t idect, Int_t i)
{
  //
  // For the calibration group idect in a detector calculate the first and last row pad and col pad.
  // The pads in the interval will have the same calibrated coefficients
  //


  Int_t posc = -1;
  Int_t posr = -1;
  frowmin[i] = -1;
  frowmax[i] = -1;
  fcolmin[i] = -1;
  fcolmax[i] = -1;
  
  if(fNfragz[i]!= 0) posc = (Int_t) idect/fNfragz[i];
  if(fNfragrphi[i] != 0) posr = (Int_t) idect%fNfragz[i];
  frowmin[i] = posr*fNnz[i];
  frowmax[i] = (posr+1)*fNnz[i];
  fcolmin[i] = posc*fNnrphi[i];
  fcolmax[i] = (posc+1)*fNnrphi[i];
}

//_____________________________________________________________________________
void AliTRDCalibra::CalculXBins(Int_t idect, Int_t i)
{
  //
  // For the detector idect calcul the first Xbins
  //


  fXbins[i] = 0;
  if(fDebug == 4) {
    AliInfo(Form("detector: %d", idect));
  }

  //In which sector?
  Int_t sector = GetSector(idect);
  fXbins[i] += sector*(6*fdetChamb2[i]+6*4*fdetChamb0[i]);
 
  //In which chamber?
  Int_t chamber = GetChamber(idect);
  Int_t kc = 0;
  while(kc < chamber){
    if(kc == 2) fXbins[i] += 6*fdetChamb2[i];
    else fXbins[i] += 6*fdetChamb0[i];
    kc ++;
  }
  
  //In wich plane?
  Int_t plane = GetPlane(idect);
  if(chamber == 2) fXbins[i] += plane*fdetChamb2[i];
  else fXbins[i] += plane*fdetChamb0[i];
 
}
//_____________________________________________________________________________
Int_t AliTRDCalibra::SearchInVector(Int_t group, Int_t i)
{
  //
  // Search if the calibration group "group" has already been initialised by a previous track in the vector
  //



  if(i == 0){
    for(Int_t k = 0; k < (Int_t) fPlaCH.size(); k++){
      if(fPlaCH[k] == group) return k;
    }
    return -1;
  }

  if(i == 1){
    for(Int_t k = 0; k < (Int_t) fPlaPH.size(); k++){
      if(fPlaPH[k] == group) return k;
    }
    return -1;
  }

  if(i == 2){
    for(Int_t k = 0; k < (Int_t) fPlaPRF.size(); k++){
      if(fPlaPRF[k] == group) return k;
    }
    return -1;
  }
  return -1;
}
//_____________________________________________________________________________
Int_t AliTRDCalibra::SearchInTreeVector(std::vector<Int_t> vectorplace, Int_t group)
{
  //
  // Search if the calibration group "group" is present in the tree
  //

  for(Int_t k = 0; k < (Int_t) vectorplace.size(); k++){
    if(vectorplace[k] == group) return k;
  }
  return -1;

}
//_____________________________________________________________________________
Int_t AliTRDCalibra::SearchBin(Float_t value, Int_t i)
{
  //
  // Search the bin
  //

  Int_t reponse = 0;
  Int_t fbinmin = 0;
  Int_t fbinmax = (Int_t) value;
  Int_t fNumberOfBin = -1;

  //charge
  if(i == 0){
    fbinmax = 300;
    fbinmin = 0;
    fNumberOfBin = fNumberBinCharge;
  }

  //PRF
  if(i == 2){
    fbinmax = 1;
    fbinmin = -1;
    fNumberOfBin = fNumberBinPRF;
  }

  //return -1 if out
  if((value >= fbinmax) || (value < fbinmin)) return -1;
  //sinon
  else{
    reponse = (Int_t)((fNumberOfBin*(value-fbinmin))/(fbinmax-fbinmin));
  }

  return reponse;
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibra::UpdateVectorCH(Int_t group, Float_t value)
{
  //
  // Fill the vector if a new calibration group "group" or update the values of the calibration group "group" if already here  
  //


  //Search bin
  Int_t bin = SearchBin(value,0);
  //out
  if((bin < 0) || (bin >= fNumberBinCharge)) return kFALSE; 
  //Search place
  Int_t place = SearchInVector(group,0);
  //new group
  if(place == -1){
    fPlaCH.push_back(group);
    //Variable
    TCTInfo *fCHInfo = new TCTInfo();
    fCHInfo->fentries = new UShort_t[fNumberBinCharge];
    //Initialise first
    for(Int_t k = 0; k < fNumberBinCharge; k++){
      fCHInfo->fentries[k] = 0;
    }
    //Add the value
    fCHInfo->fentries[bin]= 1;
  
    //Set in the vector
    fVectorCH.push_back(fCHInfo);
   
    
  }
  //group already exits
  else{
    //Variable
    TCTInfo *fCHInfo = new TCTInfo();
    fCHInfo->fentries = new UShort_t[fNumberBinCharge];
    //retrieve
    fCHInfo = fVectorCH[place];
    //add
    fCHInfo->fentries[bin]++;
    //update the vector
    std::vector<TCTInfo *>::iterator it = fVectorCH.begin()+place;
    fVectorCH.erase(it);
    fVectorCH.insert(it,fCHInfo);
   
   
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDCalibra::UpdateVectorPRF(Int_t group, Float_t x, Float_t y)
{
  //
  // Fill the vector if a new calibration group "group" or update the values of the calibration group "group" if already here  
  //

  //Search bin
  Int_t bin = SearchBin(x,2);
  //out
  if((bin < 0) || (bin >= fNumberBinPRF)) return kFALSE; 
  //Search place
  Int_t place = SearchInVector(group,2);
  //new group
  if(place == -1){
    fPlaPRF.push_back(group);
    TPInfo *fPRFInfo = new TPInfo();
    fPRFInfo->fsum = new Float_t[fNumberBinPRF];
    fPRFInfo->fsumsquare = new Float_t[fNumberBinPRF];
    fPRFInfo->fentries = new UShort_t[fNumberBinPRF];
    //Initialise first
    for(Int_t k = 0; k < fNumberBinPRF; k++){
      fPRFInfo->fsum[k] = 0.0;
      fPRFInfo->fsumsquare[k] = 0.0;
      fPRFInfo->fentries[k] = 0;
    }
    //Add the value
    fPRFInfo->fsum[bin]+= y;
    fPRFInfo->fsumsquare[bin] += y*y;
    fPRFInfo->fentries[bin]++;

    //Set in the vector
    fVectorPRF.push_back(fPRFInfo);
        
  }
  //group already exits
  else{
    TPInfo *fPRFInfo = new TPInfo();
    fPRFInfo->fsum = new Float_t[fNumberBinPRF];
    fPRFInfo->fsumsquare = new Float_t[fNumberBinPRF];
    fPRFInfo->fentries = new UShort_t[fNumberBinPRF];
    //retrieve
    fPRFInfo = fVectorPRF[place];
    //add
    Double_t calcul = (((Double_t)fPRFInfo->fentries[bin])*((Double_t)fPRFInfo->fsum[bin])+(Double_t)y)/(((Double_t)fPRFInfo->fentries[bin])+1);
    fPRFInfo->fsum[bin] = (Float_t) calcul;
    Double_t calculsquare = (((Double_t)fPRFInfo->fsumsquare[bin])*((Double_t)fPRFInfo->fentries[bin])+((Double_t)y)*((Double_t)y))/(((Double_t)fPRFInfo->fentries[bin])+1);
    fPRFInfo->fsumsquare[bin] = (Float_t)calculsquare;
    fPRFInfo->fentries[bin]++;
     
    //update the vector
    std::vector<TPInfo *>::iterator it = fVectorPRF.begin()+place;
    fVectorPRF.erase(it);
    fVectorPRF.insert(it,fPRFInfo);

  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDCalibra::UpdateVectorPH(Int_t group, Int_t time, Float_t value)
{
  //
  // Fill the vector if a new calibration group "group" or update the values of the calibration group "group" if already here  
  //

  //Search bin
  Int_t bin = time;
  //out
  if((bin < 0) || (bin >= fTimeMax)) return kFALSE; 
  //Search place
  Int_t place = SearchInVector(group,1);
  //new group
  if(place == -1){
    fPlaPH.push_back(group);
    TPInfo *fPHInfo = new TPInfo();
    fPHInfo->fsum = new Float_t[fTimeMax];
    fPHInfo->fsumsquare = new Float_t[fTimeMax];
    fPHInfo->fentries = new UShort_t[fTimeMax];
    //Initialise first
    for(Int_t k = 0; k < fTimeMax; k++){
      fPHInfo->fsum[k] = 0.0;
      fPHInfo->fsumsquare[k] = 0.0;
      fPHInfo->fentries[k] = 0;
    }
    //Add the value
    fPHInfo->fsum[bin]+= value;
    fPHInfo->fsumsquare[bin] += value*value;
    fPHInfo->fentries[bin]++;

    //Set in the vector
    fVectorPH.push_back(fPHInfo);
    
  }
  //group already exits
  else{
    TPInfo *fPHInfo = new TPInfo();
    fPHInfo->fsum = new Float_t[fTimeMax];
    fPHInfo->fsumsquare = new Float_t[fTimeMax];
    fPHInfo->fentries = new UShort_t[fTimeMax];
    //retrieve
    fPHInfo = fVectorPH[place];
    //add

    Double_t calcul = (((Double_t)fPHInfo->fentries[bin])*((Double_t)fPHInfo->fsum[bin])+(Double_t)value)/(((Double_t)fPHInfo->fentries[bin])+1);
    fPHInfo->fsum[bin] = (Float_t) calcul;
    Double_t calculsquare = ((((Double_t)fPHInfo->fsumsquare[bin])*((Double_t)fPHInfo->fentries[bin]))+(((Double_t)value)*((Double_t)value)))/(((Double_t)fPHInfo->fentries[bin])+1);
    fPHInfo->fsumsquare[bin] = (Float_t)calculsquare;
    fPHInfo->fentries[bin]++;
    //update the vector
    std::vector<TPInfo *>::iterator it = fVectorPH.begin()+place;
    fVectorPH.erase(it);
    fVectorPH.insert(it,fPHInfo);
  }

  return kTRUE;

}
  
//_____________________________________________________________________________
TGraphErrors *AliTRDCalibra::ConvertVectorPHisto(TPInfo *fPInfo, const char* name)
{
  //
  // Convert the PInfo in a 1D grapherror, name must contains "PRF" if PRF calibration and not "PRF" for Vdrift calibration
  //

  TGraphErrors *histo;
  const char *pattern1="PRF";

  //Axis
  Double_t *x;
  Double_t *y;
  Double_t *ex;
  Double_t *ey;
  Double_t step = 0.0;
  Double_t min = 0.0;

  //Ntimes
  Int_t Ntimes = 0;
  if(strstr(name,pattern1)) Ntimes = fNumberBinPRF;
  else Ntimes = fTimeMax;
  x = new Double_t[Ntimes];// xaxis
  y = new Double_t[Ntimes]; // mean
  ex = new Double_t[Ntimes]; // nentries
  ey = new Double_t[Ntimes]; // sum of square/nentries


  //init histo
  if(!strstr(name,pattern1)){
    step = 1/fSf;
    min = 0.0;
  }
  else {
    step = (1.0-(-1.0))/fNumberBinPRF;
    min = -1.0+step/2;
  }

  //Fill histo
  for(Int_t k = 0; k < Ntimes; k++){
    x[k] = min + k*step;
    y[k] = 0.0;
    ex[k] = 0.0;
    ey[k] = 0.0;
    //fill only if there is more than 0 something
    if(fPInfo->fentries[k] > 0){
      ex[k] = fPInfo->fentries[k];
      y[k] = fPInfo->fsum[k];
      ey[k] =  fPInfo->fsumsquare[k];
    }

  }

  //Define the TGraphErrors
  histo = new TGraphErrors(Ntimes, x, y, ex, ey);
  histo->SetTitle(name); 
  return histo;

 

}
 

//_____________________________________________________________________________
TH1F *AliTRDCalibra::ConvertVectorCTHisto(TCTInfo *fCTInfo, const char* name)
{
  //
  // Convert the CTInfo in a 1D histo
  //

  TH1F *histo;
  
  Int_t Ntimes = 0;
  Ntimes = fNumberBinCharge;
  

  //init histo
  histo = new TH1F(name,name,fNumberBinCharge,0,300);
  histo->Sumw2();
  //Fill histo
  for(Int_t k = 0; k < Ntimes; k++){
    histo->SetBinContent(k+1,fCTInfo->fentries[k]);
    histo->SetBinError(k+1,TMath::Sqrt(TMath::Abs(fCTInfo->fentries[k])));
  }
 


  return histo;

 

}
//_____________________________________________________________________________
TTree *AliTRDCalibra::ConvertVectorCTTreeHisto(std::vector<TCTInfo *> VectorCT, std::vector<Int_t> PlaCT, const char* name)
{
  //
  // Convert the vector in a tree with two branchs: the group number and the TH1F histo reconstructed from the vector
  //

  //Size of the things
  Int_t Ntotal = (Int_t) PlaCT.size();
  if(Ntotal == 0) {
    AliInfo("nothing to write!");
    TTree *TreeCT = new TTree(name,name);
    return TreeCT;
  }


  //Variable of the tree
  Int_t groupnumber = -1; //group calibration
  TH1F* histo = 0x0;

 
  

  //Init the tree
  TTree *TreeCT = new TTree(name,name);
  TreeCT->Branch("groupnumber",&groupnumber,"groupnumber/I");
  TreeCT->Branch("histo","TH1F",&histo,32000,0);
 
 
  
  //Fill
  Int_t k = 0;
  while(k < Ntotal){
    TString nome(name);
    groupnumber = PlaCT[0];
    nome += groupnumber;
    histo = ConvertVectorCTHisto(VectorCT[0],nome);
    
    TreeCT->Fill();
    std::vector<TCTInfo* >::iterator it = VectorCT.begin();
    VectorCT.erase(it);
    std::vector<Int_t>::iterator it2 = PlaCT.begin();
    PlaCT.erase(it2);
    k++;
 
  } 


  return TreeCT;

 

}

//_____________________________________________________________________________
TTree *AliTRDCalibra::ConvertVectorPTreeHisto(std::vector<TPInfo *> VectorP, std::vector<Int_t> PlaP, const char* name)
{
  //
  // Convert the vector in a tree with two branch: the group number and the TGraphErrors histo reconstructed from the vector.The name must contain "PRF" for PRF calibration and not "PRF" for Vdrift calibration
  //

  //Size of the things
  Int_t Ntotal = (Int_t) PlaP.size();
  if(Ntotal == 0) {
    AliInfo("nothing to write!");
    TTree *TreeP = new TTree(name,name);
    return TreeP;
  }


  //Variable of the tree
  Int_t groupnumber = -1; //group calibration
  TGraphErrors* histo = 0x0;

 
  

  //Init the tree
  TTree *TreeP = new TTree(name,name);
  TreeP->Branch("groupnumber",&groupnumber,"groupnumber/I");
  TreeP->Branch("histo","TGraphErrors",&histo,32000,0);
 
 
  
  //Fill
  Int_t k = 0;
  while(k < Ntotal){
    TString nome(name);
    groupnumber = PlaP[0];
    nome += groupnumber;
    histo = ConvertVectorPHisto(VectorP[0],nome);
    
    TreeP->Fill();
    std::vector<TPInfo*>::iterator it = VectorP.begin();
    VectorP.erase(it);
    std::vector<Int_t>::iterator it2 = PlaP.begin();
    PlaP.erase(it2);
    k++;
 
  } 


  return TreeP;

 

}

//_____________________________________________________________________________
std::vector<Int_t> AliTRDCalibra::ConvertTreeVector(TTree *tree)
{
  //
  // Convert the branch groupnumber of the tree taken from TRD.calibration.root in case of vector method in a std::vector to be faster
  //

  //initialise
  std::vector<Int_t> vectorplace;
  


  //Variable of the tree
  Int_t groupnumber = -1; //group calibration
   

  //Set the branch
  tree->SetBranchAddress("groupnumber",&groupnumber);

  
  //Fill
  Int_t Ntotal = tree->GetEntries();
  for(Int_t k = 0; k < Ntotal; k++){
    tree->GetEntry(k);
    vectorplace.push_back(groupnumber);
  } 


  return vectorplace;

}

//_____________________________________________________________________________
Bool_t AliTRDCalibra::MergeVectorCT(std::vector<TCTInfo *> VectorCT2, std::vector<Int_t> PlaCT2)
{
  //
  // Add the two vectors and place the result in the first
  //


  if(((Int_t) PlaCT2.size()) != ((Int_t) VectorCT2.size())){
    AliInfo("VectorCT2 doesn't correspond to PlaCT2!");
    return kFALSE;
  }

 
  //CH case
  for(Int_t k = 0; k < (Int_t) fPlaCH.size(); k++){
    
    //Look if PlaCT1[k] it is also in the second vector
    Int_t place = -1;
    for(Int_t j = 0; j < (Int_t) PlaCT2.size(); j++){
      if(PlaCT2[j] == fPlaCH[k]) {
	place = j;
	break;
      }
    }
    
    //If not in the second vector nothing to do
    
    
    
    //If in the second vector
    if(place != -1){
      
      
      TCTInfo *fCTInfo = new TCTInfo();
      fCTInfo->fentries = new UShort_t[fNumberBinCharge];
      
      
      for(Int_t nu = 0; nu < fNumberBinCharge; nu++){
	fCTInfo->fentries[nu] = fVectorCH[fPlaCH[k]]->fentries[nu]+VectorCT2[fPlaCH[k]]->fentries[nu];
      }
      
      //nothing to do on PlaCT1
      
      //update the vector 
      std::vector<TCTInfo *>::iterator it = fVectorCH.begin()+fPlaCH[k];
      
      fVectorCH.erase(it);
      fVectorCH.insert(it,fCTInfo);
      
    }
    
  }
    
    
  //And at the end the vector in CT2 but not in CH1
  for(Int_t k = 0; k < (Int_t) PlaCT2.size(); k++){
    
    //Look if PlaCT2[k] it is also in the second vector
    Int_t place = -1;
    for(Int_t j = 0; j < (Int_t) fPlaCH.size(); j++){
      if(fPlaCH[j] == PlaCT2[k]) {
	place = j;
	break;
      }
    }
    
    //If not in the first vector
    if(place == -1){
      
      TCTInfo *fCTInfo = new TCTInfo();
      fCTInfo->fentries = new UShort_t[fNumberBinCharge];
      
      fCTInfo = VectorCT2[PlaCT2[k]];
      
      //Add at the end 
      fPlaCH.push_back(PlaCT2[k]);
      fVectorCH.push_back(fCTInfo);
      
    }
    
  }

  
  return kTRUE;
  
  
}

//_____________________________________________________________________________
Bool_t AliTRDCalibra::MergeVectorP(std::vector<TPInfo *> VectorP2, std::vector<Int_t> PlaP2, Int_t i)
{
  //
  // Add the two vectors and place the result in the first
  //
  

   if(((Int_t) PlaP2.size()) != ((Int_t) VectorP2.size())){
     AliInfo("VectorP2 doesn't correspond to PlaP2!");
     return kFALSE;
   }

   //PH case
   if(i == 1){
     for(Int_t k = 0; k < (Int_t) fPlaPH.size(); k++){
       
       //Look if fPlaPH[k] it is also in the second vector
       Int_t place = -1;
       for(Int_t j = 0; j < (Int_t) PlaP2.size(); j++){
	 if(PlaP2[j] == fPlaPH[k]) {
	   place = j;
	   break;
	 }
       }
       
       //If not in the second vector nothing to do
          
       
       //If in the second vector
       if(place != -1){
	
	 TPInfo *fPInfo = new TPInfo();
	 fPInfo->fentries = new UShort_t[fTimeMax];
	 fPInfo->fsum = new Float_t[fTimeMax];
	 fPInfo->fsumsquare = new Float_t[fTimeMax];
      
     
	 for(Int_t nu = 0; nu < fTimeMax; nu++){
	   
	   fPInfo->fentries[nu] = fVectorPH[fPlaPH[k]]->fentries[nu]+VectorP2[fPlaPH[k]]->fentries[nu];
	   
	   Double_t calcul = ((((Double_t)fVectorPH[fPlaPH[k]]->fsum[nu])*((Double_t)fVectorPH[fPlaPH[k]]->fentries[nu]))+(((Double_t)VectorP2[fPlaPH[k]]->fsum[nu])*((Double_t)VectorP2[fPlaPH[k]]->fentries[nu])))/((Double_t)fPInfo->fentries[nu]);
	   
	   fPInfo->fsum[nu] = (Float_t) calcul;
	   
	   Double_t calculsquare = ((((Double_t)fVectorPH[fPlaPH[k]]->fsumsquare[nu])*((Double_t)fVectorPH[fPlaPH[k]]->fentries[nu]))+(((Double_t)VectorP2[fPlaPH[k]]->fsumsquare[nu])*((Double_t)VectorP2[fPlaPH[k]]->fentries[nu])))/((Double_t)fPInfo->fentries[nu]);
	   
	   
	   fPInfo->fsumsquare[nu] = calculsquare;
	 }
	 
	 //nothing to do on PlaCT1
	 
	 //update the vector VectorCT1
	 std::vector<TPInfo *>::iterator it = fVectorPH.begin()+fPlaPH[k];
	 fVectorPH.erase(it);
	 fVectorPH.insert(it,fPInfo);
	 
       }
       
     }
     
     
     //And at the end the vector in P2 but not in CH1
     for(Int_t k = 0; k < (Int_t) PlaP2.size(); k++){
       
       //Look if PlaCT2[k] it is also in the second vector
       Int_t place = -1;
       for(Int_t j = 0; j < (Int_t) fPlaPH.size(); j++){
	 if(fPlaPH[j] == PlaP2[k]) {
	   place = j;
	   break;
	 }
       }
       
       //If not in the first vector
       if(place == -1){
	 	 
	 TPInfo *fPInfo = new TPInfo();
	 fPInfo->fentries = new UShort_t[fTimeMax];
	 fPInfo->fsum = new Float_t[fTimeMax];
	 fPInfo->fsumsquare = new Float_t[fTimeMax];
	 
	 fPInfo = VectorP2[PlaP2[k]];
	 
	 //Add at the end of CH1
	 fPlaPH.push_back(PlaP2[k]);
	 fVectorPH.push_back(fPInfo);
	 
       }
     }
   }
   

   //PRF case
   if(i == 1){
     for(Int_t k = 0; k < (Int_t) fPlaPRF.size(); k++){

       //Look if fPlaPRF[k] it is also in the second vector
       Int_t place = -1;
       for(Int_t j = 0; j < (Int_t) PlaP2.size(); j++){
	 if(PlaP2[j] == fPlaPRF[k]) {
	   place = j;
	   break;
	 }
       }
       
       //If not in the second vector nothing to do
       
       
       
       //If in the second vector
       if(place != -1){
	
	 TPInfo *fPInfo = new TPInfo();
	 fPInfo->fentries = new UShort_t[fNumberBinPRF];
	 fPInfo->fsum = new Float_t[fNumberBinPRF];
	 fPInfo->fsumsquare = new Float_t[fNumberBinPRF];
      
     
	 for(Int_t nu = 0; nu < fNumberBinPRF; nu++){
	   
	   fPInfo->fentries[nu] = fVectorPRF[fPlaPRF[k]]->fentries[nu]+VectorP2[fPlaPRF[k]]->fentries[nu];
	   
	   Double_t calcul = ((((Double_t)fVectorPRF[fPlaPRF[k]]->fsum[nu])*((Double_t)fVectorPRF[fPlaPRF[k]]->fentries[nu]))+(((Double_t)VectorP2[fPlaPRF[k]]->fsum[nu])*((Double_t)VectorP2[fPlaPRF[k]]->fentries[nu])))/((Double_t)fPInfo->fentries[nu]);
	   
	   fPInfo->fsum[nu] = (Float_t) calcul;
	   
	   Double_t calculsquare = ((((Double_t)fVectorPRF[fPlaPRF[k]]->fsumsquare[nu])*((Double_t)fVectorPRF[fPlaPRF[k]]->fentries[nu]))+(((Double_t)VectorP2[fPlaPRF[k]]->fsumsquare[nu])*((Double_t)VectorP2[fPlaPRF[k]]->fentries[nu])))/((Double_t)fPInfo->fentries[nu]);
	   
	   
	   fPInfo->fsumsquare[nu] = calculsquare;
	 }
	 
	 //nothing to do on PlaCT1
	 
	 //update the vector VectorCT1
	 std::vector<TPInfo *>::iterator it = fVectorPRF.begin()+fPlaPRF[k];
	 fVectorPRF.erase(it);
	 fVectorPRF.insert(it,fPInfo);
	 
       }
       
     }
     
     
     //And at the end the vector in P2 but not in CH1
     for(Int_t k = 0; k < (Int_t) PlaP2.size(); k++){
       
       //Look if PlaCT2[k] it is also in the second vector
       Int_t place = -1;
       for(Int_t j = 0; j < (Int_t) fPlaPRF.size(); j++){
	 if(fPlaPRF[j] == PlaP2[k]) {
	   place = j;
	   break;
	 }
       }
       
       //If not in the first vector
       if(place == -1){
	 
	 TPInfo *fPInfo = new TPInfo();
	 fPInfo->fentries = new UShort_t[fNumberBinPRF];
	 fPInfo->fsum = new Float_t[fNumberBinPRF];
	 fPInfo->fsumsquare = new Float_t[fNumberBinPRF];
	 
	 fPInfo = VectorP2[PlaP2[k]];
	 
	 //Add at the end of CH1
	 fPlaPRF.push_back(PlaP2[k]);
	 fVectorPRF.push_back(fPInfo);
	 
       }
       
     }
   } 
   
   
   return kTRUE;
   
   
}
//_____________Fit Methods______________________________________________________________________________________________________________________

//______________________________________________________________________________________________________
void AliTRDCalibra::FitPente(TH1* projPH, Int_t idect)
{
  //
  // Slope methode for the drift velocity
  //
  
  //constants
  const Float_t kDrWidth = AliTRDgeometry::DrThick();
  Int_t binmax = 0;
  Int_t binmin = 0;
  fPhd[0] = 0.0;
  fPhd[1] = 0.0;
  fPhd[2] = 0.0;
  Int_t ju = 0;
  fVdriftCoef[1] = 0.0;
  fT0Coef[1] = 0.0;
  TLine *line = new TLine();
   
  

  //beginning of the signal
  TH1D *pentea = new TH1D("pentea", "pentea", projPH->GetNbinsX(),0,(Float_t)(fTimeMax)/fSf);
  for(Int_t k = 1; k <  projPH->GetNbinsX(); k++){
    pentea->SetBinContent(k,(Double_t)(projPH->GetBinContent(k+1)-projPH->GetBinContent(k)));
  }
  
  binmax = (Int_t)pentea->GetMaximumBin();
  if(binmax == 1) {
    binmax = 2;
    AliInfo("Put the binmax from 1 to 2 to enable the fit");
  }
  pentea->Fit("pol2","0MR","",TMath::Max(pentea->GetBinCenter(binmax-1),0.0),pentea->GetBinCenter(binmax+1));
  Float_t l3P1am = pentea->GetFunction("pol2")->GetParameter(1);
  Float_t l3P2am = pentea->GetFunction("pol2")->GetParameter(2);
  
  if(l3P2am != 0) fPhd[0] = -(l3P1am/(2*l3P2am));
  
  //amplification region
  binmax = 0;
  ju = 0;
  for(Int_t kbin = 1; kbin < projPH->GetNbinsX(); kbin ++){
    if(((projPH->GetBinContent(kbin+1)-projPH->GetBinContent(kbin)) <= 0.0) && (ju == 0)){
      binmax = kbin;
      ju = 1;
    }
  }
  if(binmax == 1) {
    binmax = 2;
    AliInfo("Put the binmax from 1 to 2 to enable the fit");
  }
  projPH->Fit("pol2","0MR","",TMath::Max(projPH->GetBinCenter(binmax-1),0.0),projPH->GetBinCenter(binmax+1));
  Float_t l3P1amf = projPH->GetFunction("pol2")->GetParameter(1);
  Float_t l3P2amf = projPH->GetFunction("pol2")->GetParameter(2);
  
  if(l3P2amf != 0) fPhd[1] = -(l3P1amf/(2*l3P2amf));
  

  //Drift region
  TH1D *pente = new TH1D("pente", "pente", projPH->GetNbinsX(),0,(Float_t)(fTimeMax)/fSf);
  for(Int_t k = binmax+4; k <  projPH->GetNbinsX(); k++){
    pente->SetBinContent(k,(Double_t)(projPH->GetBinContent(k+1)-projPH->GetBinContent(k)));
  }
  binmin = (Int_t)pente->GetMinimumBin();
  if(binmin == 1) {
    binmin = 2;
    AliInfo("Put the binmax from 1 to 2 to enable the fit");
  }
  pente->Fit("pol2","0MR","",TMath::Max(pente->GetBinCenter(binmin-1),0.0),TMath::Min(pente->GetBinCenter(binmin+2),(Double_t)(fTimeMax)/fSf));
  Float_t l3P1dr = pente->GetFunction("pol2")->GetParameter(1);
  Float_t l3P2dr = pente->GetFunction("pol2")->GetParameter(2);
  if(l3P2dr != 0) fPhd[2] = -(l3P1dr/(2*l3P2dr));


  if((fPhd[2] > fPhd[0]) && (fPhd[2] > fPhd[1])){
    fVdriftCoef[1] = (kDrWidth)/(fPhd[2]-fPhd[1]);
    if(fPhd[0] >= 0.0) fT0Coef[1] = fPhd[0]-fT0Shift;
    else fT0Coef[1] = -TMath::Abs(fT0Coef[2]);
  }
  else{
    fVdriftCoef[1] = -TMath::Abs(fVdriftCoef[2]);
    fT0Coef[1] = -TMath::Abs(fT0Coef[2]);
  }
   
  if((fDebug == 1) || (fDebug == 4)){

    fCoefVdrift[1]->SetBinContent(idect+1,TMath::Abs(fVdriftCoef[1]));
    fCoefT0[1]->SetBinContent(idect+1,TMath::Abs(fT0Coef[1]));
    if(fVdriftCoef[1] > 0.0){
      if(fVdriftCoef[2] != 0.0) fDeltaVdrift[1]->SetBinContent(idect+1,(fVdriftCoef[1]-fVdriftCoef[2])/fVdriftCoef[2]);
      fDeltaT0[1]->SetBinContent(idect+1,(fT0Coef[1]-fT0Coef[2]));
    }
  }
  

  if(fDebug == 2){
    TCanvas *cpentei = new TCanvas("cpentei","cpentei",50,50,600,800);
    cpentei->cd();
    projPH->Draw();
    line->SetLineColor(2);
    line->DrawLine(fPhd[0],0,fPhd[0],projPH->GetMaximum());
    line->DrawLine(fPhd[1],0,fPhd[1],projPH->GetMaximum());
    line->DrawLine(fPhd[2],0,fPhd[2],projPH->GetMaximum());
    
    AliInfo(Form("fPhd[0] (beginning of the signal): %f fPhd[1] (end of the amplification region): %f fPhd[2] (end of the drift region): %f fVriftCoef[1] (with only the drift region(default)): %f",(Float_t) fPhd[0], (Float_t) fPhd[1],(Float_t) fPhd[2], (Float_t) fVdriftCoef[1]));
    
  }
 if(fDebug != 2) delete pentea;
 if(fDebug != 2) delete pente;
  
}

//__________________________________________________________________________________________
void AliTRDCalibra::FitPH(TH1* projPH, Int_t idect)
{
  //
  // Fit methode for the drift velocity
  //
  
  //constants
  const Float_t kDrWidth = AliTRDgeometry::DrThick();   

  TF1 *fPH = new TF1("fPH", AliTRDCalibra::PH,-0.05,3.2,6);
  fPH->SetParameter(0,0.469);  // scaling
  fPH->SetParameter(1,0.18);  // start 
  fPH->SetParameter(2,0.0857325);  // AR
  fPH->SetParameter(3,1.89);  // DR
  fPH->SetParameter(4,0.08);  // QA/QD
  fPH->SetParameter(5,0.0);//baseline

  TLine *line = new TLine();

  fVdriftCoef[0] = 0.0;
  fT0Coef[0] = 0.0;
 
  if(idect%fFitPHPeriode == 0){
    AliInfo(Form("<AliTRDCalibra::FitPH> The detector %d will be fitted", idect));
    fPH->SetParameter(0,(projPH->Integral()-(projPH->GetBinContent(1)*projPH->GetNbinsX()))*0.00028);  // scaling
    fPH->SetParameter(1,fPhd[0]-0.1);  // start 
    fPH->SetParameter(2,fPhd[1]-fPhd[0]);  // AR
    fPH->SetParameter(3,fPhd[2]-fPhd[1]);  // DR
    fPH->SetParameter(4,0.225);  // QA/QD
    fPH->SetParameter(5,(Float_t) projPH->GetBinContent(1));
    
    if(fDebug != 2){
      projPH->Fit(fPH,"0M","",0.0,(Float_t)fTimeMax/fSf);
    }

    if(fDebug == 2){
      TCanvas *cpente = new TCanvas("cpente","cpente",50,50,600,800);
      cpente->cd();
      projPH->Fit(fPH,"M+","",0.0,(Float_t)fTimeMax/fSf);
      projPH->Draw("E0");
      line->SetLineColor(4);
      line->DrawLine(fPH->GetParameter(1),0,fPH->GetParameter(1),projPH->GetMaximum());
      line->DrawLine(fPH->GetParameter(1)+fPH->GetParameter(2),0,fPH->GetParameter(1)+fPH->GetParameter(2),projPH->GetMaximum());
      line->DrawLine(fPH->GetParameter(1)+fPH->GetParameter(2)+fPH->GetParameter(3),0,fPH->GetParameter(1)+fPH->GetParameter(2)+fPH->GetParameter(3),projPH->GetMaximum());
    }
    

    if(fPH->GetParameter(3) !=0){
      fVdriftCoef[0] = kDrWidth/(fPH->GetParameter(3));
      fT0Coef[0] = fPH->GetParameter(1);
    } 
    else {
      fVdriftCoef[0] = -TMath::Abs(fVdriftCoef[2]);
      fT0Coef[0] = -TMath::Abs(fT0Coef[2]);
    }
    if((fDebug == 1) || (fDebug == 4)) {
      fCoefVdrift[0]->SetBinContent(idect+1,TMath::Abs(fVdriftCoef[0]));
      fCoefT0[0]->SetBinContent(idect+1,TMath::Abs(fT0Coef[0]));
    if(fVdriftCoef[0] > 0.0){
	if(fVdriftCoef[2] != 0.0) fDeltaVdrift[0]->SetBinContent(idect+1,(fVdriftCoef[0]-fVdriftCoef[2])/fVdriftCoef[2]);
	fDeltaT0[0]->SetBinContent(idect+1,(fT0Coef[0]-fT0Coef[2]));
      }
    }
    if(fDebug == 2){
      AliInfo(Form("fVdriftCoef[0]: %f", (Float_t) fVdriftCoef[0]));
    }
    
  }

  else{
    //Put the default value 
    if((fDebug <= 1) || (fDebug = 4)) {
      fCoefVdrift[0]->SetBinContent(idect+1,fVdriftCoef[2]);
      fCoefT0[0]->SetBinContent(idect+1,fT0Coef[2]);
    }

  }

  if(fDebug != 2) delete fPH;
  
}
//__________________________________________________________________________________________
void AliTRDCalibra::FitPRF(TH1* projPRF, Int_t idect)
{
  //
  // Fit methode for the sigma of the pad response function
  //
  
  fPRFCoef[0] = 0.0;

 
  
  if(fDebug != 2){
    projPRF->Fit("gaus","0M","", -fRangeFitPRF, fRangeFitPRF);
  }
  
  if(fDebug == 2){
    TCanvas *cfit = new TCanvas("cfit","cfit",50,50,600,800);
    cfit->cd();
    projPRF->Fit("gaus","M+","", -fRangeFitPRF, fRangeFitPRF);
    projPRF->Draw();
    
  }
  
  
  fPRFCoef[0] = projPRF->GetFunction("gaus")->GetParameter(2);
 

  if((fDebug == 1) || (fDebug == 4)){
    fCoefPRF[0]->SetBinContent(idect+1,fPRFCoef[0]);
    if(fPRFCoef[1] != 0.0){
      fDeltaPRF->SetBinContent(idect+1,(fPRFCoef[0]-fPRFCoef[1])/fPRFCoef[1]);
    }
  }
  if(fDebug == 2){
    AliInfo(Form("fPRFCoef[0]: %f", (Float_t) fPRFCoef[0]));
  }
  
}
//________________________________________________________________________________________
void AliTRDCalibra::FitCH( TH1* projch, Int_t idect)
{
  //
  // Fit methode for the gain factor
  //
 
  fChargeCoef[0] = 0.0;
  fChargeCoef[1] = 0.0;
  TF1 *fLandauGaus = new TF1("fLandauGaus",funcLandauGaus,0,300,5);

  fChargeCoef[1] = projch->GetMean();
  projch->Fit("landau","0","",(Float_t) fChargeCoef[1]/fBeginFitCharge,projch->GetBinCenter(projch->GetNbinsX()));
  fl3P0 = projch->GetFunction("landau")->GetParameter(0);
  Double_t l3P1 = projch->GetFunction("landau")->GetParameter(1);
  fl3P2 = projch->GetFunction("landau")->GetParameter(2);
    
  projch->Fit("gaus","0","",(Float_t) fChargeCoef[1]/fBeginFitCharge,projch->GetBinCenter(projch->GetNbinsX()));
  Double_t g3P0 = projch->GetFunction("gaus")->GetParameter(0);
  fg3P2 = projch->GetFunction("gaus")->GetParameter(2);
   
        
  fLandauGaus->SetParameters(fl3P0,l3P1,fl3P2,g3P0,fg3P2);
  if((fDebug <= 1) || (fDebug >= 3)){
    projch->Fit("fLandauGaus","0","",(Float_t) fChargeCoef[1]/fBeginFitCharge,projch->GetBinCenter(projch->GetNbinsX()));
  }
  if(fDebug == 2){
    TCanvas *cp = new TCanvas("cp","cp",50,50,600,800);
    cp->cd();
    projch->Fit("fLandauGaus","+","",(Float_t) fChargeCoef[1]/fBeginFitCharge,projch->GetBinCenter(projch->GetNbinsX()));
    projch->Draw();
    fLandauGaus->Draw("same");
  }
  
  if(projch->GetFunction("fLandauGaus")->GetParameter(1) > 0){

    //Calcul of "real" coef****************************************
    CalculChargeCoefMean(fcountdet[0],(Int_t) idect, kTRUE);
    fChargeCoef[0] = projch->GetFunction("fLandauGaus")->GetParameter(1);}
  
  else {
    
    //Calcul of "real" coef****************************************
    CalculChargeCoefMean(fcountdet[0],(Int_t) idect, kFALSE);
    fChargeCoef[0] = -TMath::Abs(fChargeCoef[3]);
  }
  if(fDebug == 2){
    AliInfo(Form("fChargeCoef[0]: %f", (Float_t) fChargeCoef[0]));
    AliInfo(Form("fChargeCoef[1]: %f",(Float_t) fChargeCoef[1]));
  }
  
  if((fDebug == 1) || (fDebug == 4)) {
    if(fChargeCoef[0] > 0.0){
      fCoefCharge[0]->SetBinContent(idect+1,fChargeCoef[0]);
      fCoefCharge[1]->SetBinContent(idect+1,fChargeCoef[1]);
      fDeltaCharge[0]->SetBinContent(idect+1,fChargeCoef[0]);
      fDeltaCharge[1]->SetBinContent(idect+1,fChargeCoef[1]);
    }
  }
  fl3P0 = fLandauGaus->Integral(0.3*projch->GetMean(),3*projch->GetMean());
  fg3P2 = fLandauGaus->GetParameter(2);
  fl3P2 = fLandauGaus->GetParameter(4);
   
  if(fDebug != 2) delete fLandauGaus;
}
//_________________________________________________________________________________
void AliTRDCalibra::FitBisCH( TH1* projch, Int_t idect)
{
  //
  // Fit methode for the gain factor more time consuming
  //
  

 
  // Setting fit range and start values
  Double_t fr[2];
  //Double_t sv[4] = {l3P2,fChargeCoef[1],projch->Integral("width"),fg3P2};
  Double_t sv[4] = {fl3P2,fChargeCoef[1],fl3P0,fg3P2};
  Double_t  pllo[4] = {0.001,0.001,0.001,0.001};
  Double_t plhi[4] = {300.0,300.0,100000000.0,300.0};
  Double_t fp[4] = {1.0,1.0,1.0,1.0};
  Double_t  fpe[4] = {1.0,1.0,1.0,1.0};
  fr[0]=0.3*projch->GetMean();
  fr[1]=3.0*projch->GetMean();
  fChargeCoef[2] = 0.0;
    
       
  Double_t chisqr;
  Int_t    ndf;
  TF1 *fitsnr = langaufit(projch,&fr[0],&sv[0],&pllo[0],&plhi[0],&fp[0],&fpe[0],&chisqr,&ndf);
    
  Double_t projchPeak, projchFWHM;
  langaupro(fp,projchPeak,projchFWHM);
  
  if(fp[1] > 0){
    fChargeCoef[2] = fp[1];} 
  else fChargeCoef[2] = -TMath::Abs(fChargeCoef[3]);
  
  
  if(fDebug == 2){
    AliInfo(Form("fChargeCoef[2]: %f", (Float_t) fChargeCoef[2]));
    TCanvas *cpy = new TCanvas("cpy","cpy",50,50,600,800);
    cpy->cd();
    projch->Draw();
    fitsnr->Draw("same");
  }
 
  if((fDebug == 1) || (fDebug ==4)) {
    if(fChargeCoef[2] > 0.0){
      fCoefCharge[2]->SetBinContent(idect+1,fChargeCoef[2]);
      fDeltaCharge[2]->SetBinContent(idect+1,fChargeCoef[2]);
    }
  }
  
  if(fDebug != 2) delete fitsnr;
  
} 
//________________________________________________________________________________
void AliTRDCalibra::NormierungCharge()
{
  //
  // normalisation of the gain factor resulting for the fits
  //
  
  //Calcul of the mean of the fit
  Double_t sum = 0.0;
  Float_t scalefactor = 1.0;
  for(Int_t k = 0; k < (Int_t) fVectorFitCH.size(); k++){
    Int_t total = 0;
    if(GetChamber(fVectorFitCH[k]->fDetector) == 2) total = 1728;
    if(GetChamber(fVectorFitCH[k]->fDetector) != 2) total = 2304;
    for(Int_t j = 0; j < total; j++){
      if(fVectorFitCH[k]->fcoef[j] >= 0){
	sum += (Double_t)fVectorFitCH[k]->fcoef[j];
      }
    }
  }
 
  if(sum > 0){
    fScalefitfactor = fScalefitfactor/sum;
  }
  else fScalefitfactor = 1.0;
  //Scale the histo
  if((fDebug == 1) || (fDebug == 4)){
    if((fCoefCharge[0]->GetEntries() > 0.0) && (fCoefCharge[0]->GetSumOfWeights() > 0.0)){
      scalefactor = fCoefCharge[0]->GetEntries()/fCoefCharge[0]->GetSumOfWeights();
      fCoefCharge[0]->Scale(scalefactor);
      fDeltaCharge[0]->Scale(scalefactor);
    }
    if((fMeanChargeOn) && (fCoefCharge[1]->GetEntries() > 0.0) && (fCoefCharge[1]->GetSumOfWeights() > 0.0)){
      fCoefCharge[1]->Scale(fCoefCharge[1]->GetEntries()/fCoefCharge[1]->GetSumOfWeights());
    }
    if((fFitChargeBisOn) && (fCoefCharge[2]->GetEntries() > 0.0) && (fCoefCharge[2]->GetSumOfWeights() > 0.0)){
      fCoefCharge[2]->Scale(fCoefCharge[2]->GetEntries()/fCoefCharge[2]->GetSumOfWeights());
    }
    if((fMeanChargeOn) && (fDeltaCharge[1]->GetEntries() > 0.0) && (fDeltaCharge[1]->GetSumOfWeights() > 0.0)){
      fDeltaCharge[1]->Scale(fDeltaCharge[1]->GetEntries()/fDeltaCharge[1]->GetSumOfWeights());
    }
    if((fFitChargeBisOn) && (fDeltaCharge[2]->GetEntries() > 0.0) && (fDeltaCharge[2]->GetSumOfWeights() > 0.0)){
      fDeltaCharge[2]->Scale(fDeltaCharge[2]->GetEntries()/fDeltaCharge[2]->GetSumOfWeights());
    }
    
  }

  if((fDebug == 3) || (fDebug == 4)){
    fCoefChargeDB[0]->Scale(scalefactor);
    if((fMeanChargeOn)  && (fCoefChargeDB[1]->GetEntries() > 0.0) && (fCoefChargeDB[1]->GetSumOfWeights() > 0.0))fCoefChargeDB[1]->Scale(fCoefChargeDB[1]->GetEntries()/fCoefChargeDB[1]->GetSumOfWeights());
    if((fFitChargeBisOn) && (fCoefChargeDB[2]->GetEntries() > 0.0) && (fCoefChargeDB[2]->GetSumOfWeights() > 0.0))fCoefChargeDB[2]->Scale(fCoefChargeDB[2]->GetEntries()/fCoefChargeDB[2]->GetSumOfWeights());
  }
  
  
  
  if((fDebug == 1) || (fDebug == 4)){  
    fDeltaCharge[0]->Add(fCoefCharge[3],-1);
    fDeltaCharge[0]->Divide(fCoefCharge[3]);
    
    if(fMeanChargeOn) fDeltaCharge[1]->Add(fCoefCharge[3],-1);
    if(fMeanChargeOn) fDeltaCharge[1]->Divide(fCoefCharge[3]);
    
    if(fFitChargeBisOn) fDeltaCharge[2]->Add(fCoefCharge[3],-1);
    if(fFitChargeBisOn) fDeltaCharge[2]->Divide(fCoefCharge[3]);
  }
  
}

//_______________________________________________________________________________________________
TH1I *AliTRDCalibra::ReBin(TH1I *hist)
{
  //
  // Rebin of the 1D histo for the gain calibration if needed.
  // you have to choose fRebin, divider of fNumberBinCharge
  //

 TAxis *xhist = hist->GetXaxis();
 TH1I* rehist = new TH1I("projrebin","",(Int_t) xhist->GetNbins()/fRebin, xhist->GetBinLowEdge(1), xhist->GetBinUpEdge(xhist->GetNbins()));
 AliInfo(Form("fRebin: %d", fRebin));
 Int_t i = 1;
 for(Int_t k = 1; k <= (Int_t) xhist->GetNbins()/fRebin; k++){
   Double_t sum = 0.0;
   for(Int_t ji = i; ji < i+fRebin; ji++){
     sum += hist->GetBinContent(ji);
   }
   sum = sum/fRebin;
   rehist->SetBinContent(k,sum);
   i +=fRebin;
   
   
 }
 if(fDebug == 2){
   TCanvas *crebin = new TCanvas("crebin","",50,50,600,800);
   crebin->cd();
   rehist->Draw();
 }
 
 return rehist;

}

//_______________________________________________________________________________________________
TH1F *AliTRDCalibra::ReBin(TH1F *hist)
{
  //
  // Rebin of the 1D histo for the gain calibration if needed
  // you have to choose fRebin divider of fNumberBinCharge
  //

  TAxis *xhist = hist->GetXaxis();
  TH1F* rehist = new TH1F("projrebin","",(Int_t) xhist->GetNbins()/fRebin, xhist->GetBinLowEdge(1), xhist->GetBinUpEdge(xhist->GetNbins()));
  AliInfo(Form("fRebin: %d", fRebin));
  Int_t i = 1;
  for(Int_t k = 1; k <= (Int_t) xhist->GetNbins()/fRebin; k++){
    Double_t sum = 0.0;
    for(Int_t ji = i; ji < i+fRebin; ji++){
      sum += hist->GetBinContent(ji);
    }
    sum = sum/fRebin;
    rehist->SetBinContent(k,sum);
    i +=fRebin;
    
    
  }
  if(fDebug == 2){
    TCanvas *crebin = new TCanvas("crebin","",50,50,600,800);
    crebin->cd();
    rehist->Draw();
  }
  
  return rehist;
  
}
//_______________________________________________________________________________________________
TH1F *AliTRDCalibra::CorrectTheError(TGraphErrors *hist)
{
  //
  // In the case of the vectors method the trees contains TGraphErrors for PH and PRF
  // to be able to add them after
  // We convert it to a TH1F to be able to applied the same fit function method
  // After having called this function you can not add the statistics anymore
  //

  TH1F* rehist = 0x0;

  Int_t nbins = hist->GetN();
  Double_t *x = hist->GetX();
  Double_t *entries = hist->GetEX();
  Double_t *mean = hist->GetY();
  Double_t *square = hist->GetEY();
  fEntriesCurrent = 0;


  if(nbins < 2) return rehist; 
  Double_t step = x[1]-x[0]; 
  Double_t minvalue = x[0]-step/2;
  Double_t maxvalue = x[(nbins-1)]+step/2;


  rehist = new TH1F("projcorrecterror","",nbins,minvalue,maxvalue);

  for(Int_t k = 0; k < nbins; k++){
    rehist->SetBinContent(k+1,mean[k]);
    if(entries[k] > 0.0){
      fEntriesCurrent += (Int_t) entries[k];
      Double_t d = TMath::Abs(square[k]- (mean[k]*mean[k]));
      rehist->SetBinError(k+1,TMath::Sqrt(d/entries[k]));
    }
    else rehist->SetBinError(k+1,0.0);
   
  }
 return rehist;
 
}
//_______________________________________________________________________________________________
TGraphErrors *AliTRDCalibra::AddProfiles(TGraphErrors *hist1, TGraphErrors *hist2)
{
  //
  // In the case of the vectors method we use TGraphErrors for PH and PRF
  // to be able to add the them after
  // Here we add the TGraphErrors  
 //

  
  //First TGraphErrors
  Int_t nbins1 = hist1->GetN();
  Double_t *x1 = hist1->GetX();
  Double_t *ex1 = hist1->GetEX();
  Double_t *y1 = hist1->GetY();
  Double_t *ey1 = hist1->GetEY();

  TGraphErrors* rehist = new TGraphErrors(nbins1);

  

  //Second TGraphErrors
  Double_t *ex2 = hist2->GetEX();
  Double_t *y2 = hist2->GetY();
  Double_t *ey2 = hist2->GetEY();
   
  //Define the Variables for the new TGraphErrors
  Double_t x, ex, y, ey; // Xaxis
  
  
 for(Int_t k = 0; k < nbins1; k++){
   x = x1[k];
   Double_t nentries = 0.0;
   y = 0.0;
   ey = 0.0;
   ex = 0.0;
   if((ex2[k] == 0.0) && (ex1[k] == 0.0)) nentries = 0.0;
   if((ex2[k] == 0.0) && (ex1[k] > 0.0)) {
     nentries = ex1[k];
     y = y1[k];
     ey = ey1[k];
     ex = ex1[k];
   }
   if((ex2[k] > 0.0) && (ex1[k] == 0.0)) {
     nentries = ex2[k];
     y = y2[k];
     ey = ey2[k];
     ex = ex2[k];
   }
   if((ex2[k] > 0.0) && (ex1[k] > 0.0)){ 
     nentries = ex1[k]+ex2[k];
     y = (y1[k]*ex1[k]+y2[k]*ex2[k])/nentries;
     ey = (ey1[k]*ex1[k]+ey2[k]*ex2[k])/nentries;
     ex = nentries;
   }
   rehist->SetPoint(k,x,y);
   rehist->SetPointError(k,ex,ey);
    
 }


 return rehist;
 
}


//_________________Some basic geometry function______________________________________________________________________________________________________

//_____________________________________________________________________________
Int_t AliTRDCalibra::GetPlane(Int_t d) const
{
  //
  //Reconstruct the plane number from the detector number
  //

  return ((Int_t) (d % 6));

}

//_____________________________________________________________________________
Int_t AliTRDCalibra::GetChamber(Int_t d) const
{
  //
  //Reconstruct the chamber number from the detector number
  //
  Int_t fgkNplan = 6;

  return ((Int_t) (d % 30) / fgkNplan);

}

//_____________________________________________________________________________
Int_t AliTRDCalibra::GetSector(Int_t d) const
{
  //
  //Reconstruct the sector number from the detector number
  //
  Int_t fg = 30;

  return ((Int_t) (d / fg));

}

//____________________________Fill and Init tree Gain, PRF, Vdrift and T0__________________________________________________________________________________________________________________

//_____________________________________________________________________________
void AliTRDCalibra::InittreePRF()
{
  //
  // Init the tree where the coefficients from the fit methods can be stored
  //
 
  fPRFPad = new Float_t[2304];
  fPRF = new TTree("PRF","PRF");
  fPRF->Branch("detector",&fPRFDetector,"detector/I");
  fPRF->Branch("width",fPRFPad,"width[2304]/F");
  //Set to default value for the plane 0 supposed to be the first one
  for(Int_t k = 0; k < 2304; k++){
    fPRFPad[k] = 0.515;
  }
  fPRFDetector = -1;


 

}
//_____________________________________________________________________________
void AliTRDCalibra::FilltreePRF(Int_t countdet)
{
  //
  // Fill the tree with the sigma of the pad response function for the detector countdet
  //
  
  Int_t numberofgroup = 0; 
  fPRFDetector = countdet;
  fPRF->Fill();
  if(GetChamber((Int_t)(countdet+1)) == 2) numberofgroup = 1728;
  else numberofgroup = 2304;
  //Reset to default value for the next
  for(Int_t k = 0; k < numberofgroup; k++){
    if(GetPlane((Int_t) (countdet+1)) == 0) fPRFPad[k] = 0.515;
    if(GetPlane((Int_t) (countdet+1)) == 1) fPRFPad[k] = 0.502;
    if(GetPlane((Int_t) (countdet+1)) == 2) fPRFPad[k] = 0.491;
    if(GetPlane((Int_t) (countdet+1)) == 3) fPRFPad[k] = 0.481;
    if(GetPlane((Int_t) (countdet+1)) == 4) fPRFPad[k] = 0.471;
    if(GetPlane((Int_t) (countdet+1)) == 5) fPRFPad[k] = 0.463;
  }
  fPRFDetector = -1;

}

//_____________________________________________________________________________
void AliTRDCalibra::ConvertVectorFitCHTree()
{
  //
  // Convert the vector stuff to a tree of 1D histos if the user want to write it after the fill functions
  //
  
  Int_t detector = -1;
  Int_t numberofgroup = 1;
  Float_t gainPad[2304];
  fGain = new TTree("Gain","Gain");
  fGain->Branch("detector",&detector,"detector/I");
  fGain->Branch("gainPad",gainPad,"gainPad[2304]/F");
  Int_t loop = (Int_t) fVectorFitCH.size();
  for(Int_t k = 0; k < loop; k++){
    detector = fVectorFitCH[k]->fDetector;
    if(GetChamber((Int_t)fVectorFitCH[k]->fDetector) == 2) numberofgroup = 1728;
    else numberofgroup = 2304;
    for(Int_t i = 0; i < numberofgroup; i++){
      if(fVectorFitCH[k]->fcoef[i] >= 0) gainPad[i] = fVectorFitCH[k]->fcoef[i]*fScalefitfactor;
      else {
	gainPad[i] = (Float_t) fVectorFitCH[k]->fcoef[i];
      }
    }
    fGain->Fill();
  }
 
  

 

}
//_____________________________________________________________________________
void AliTRDCalibra::FilltreeVdrift(Int_t countdet)
{
  //
  // Fill the tree with the drift velocities for the detector countdet
  //
  Int_t numberofgroup = 0;
  fVdriftDetector = countdet;
  for(Int_t k = 0; k < 2304; k++){
  }
  fVdrift->Fill();
  if(GetChamber((Int_t)(countdet+1)) == 2) numberofgroup = 1728;
  else numberofgroup = 2304;
  //Reset to default value the gain coef
  for(Int_t k = 0; k < numberofgroup; k++){
    fVdriftPad[k] = -1.5;
  }
  fVdriftDetector = -1;

}

//_____________________________________________________________________________
void AliTRDCalibra::InittreePH()
{
  //
  // Init the tree where the coefficients from the fit methods can be stored
  //
  
  fVdriftPad = new Float_t[2304];
  fVdrift = new TTree("Vdrift","Vdrift");
  fVdrift->Branch("detector",&fVdriftDetector,"detector/I");
  fVdrift->Branch("vdrift",fVdriftPad,"vdrift[2304]/F");
  //Set to default value for the plane 0 supposed to be the first one
  for(Int_t k = 0; k < 2304; k++){
    fVdriftPad[k] = -1.5;
  }
  fVdriftDetector = -1;

}
//_____________________________________________________________________________
void AliTRDCalibra::FilltreeT0(Int_t countdet)
{
  //
  // Fill the tree with the t0 value for the detector countdet
  //
  Int_t numberofgroup = 0;
  fT0Detector = countdet;
  for(Int_t k = 0; k < 2304; k++){
  }
  fT0->Fill();
  if(GetChamber((Int_t)(countdet+1)) == 2) numberofgroup = 1728;
  else numberofgroup = 2304;
  //Reset to default value the gain coef
  for(Int_t k = 0; k < numberofgroup; k++){
    fT0Pad[k] = 0.0;
  }
  fT0Detector = -1;

}

//_____________________________________________________________________________
void AliTRDCalibra::InittreeT0()
{
  //
  // Init the tree where the coefficients from the fit methods can be stored
  //
  
  fT0Pad = new Float_t[2304];
  fT0 = new TTree("T0","T0");
  fT0->Branch("detector",&fT0Detector,"detector/I");
  fT0->Branch("t0",fT0Pad,"t0[2304]/F");
  //Set to default value for the plane 0 supposed to be the first one
  for(Int_t k = 0; k < 2304; k++){
    fT0Pad[k] = 0.0;
  }
  fT0Detector = -1;

}

//________________private Functions______________________________________________________________________________________________________________


//_________________________________________________________________________________

Double_t AliTRDCalibra::PH(Double_t* x, Double_t* par) 
{
  //
  // Function for the fit
  //

  //TF1 * fAsymmGauss = new TF1("fAsymmGauss",AsymmGauss,0,4,6);

  //PARAMETERS FOR FIT PH
  // PASAv.4
  //fAsymmGauss->SetParameter(0, 0.113755);
  //fAsymmGauss->SetParameter(1, 0.350706);
  //fAsymmGauss->SetParameter(2, 0.0604244);
  //fAsymmGauss->SetParameter(3, 7.65596);
  //fAsymmGauss->SetParameter(4, 1.00124);
  //fAsymmGauss->SetParameter(5, 0.870597);  //no tail cancelation

  Double_t xx = x[0];
  
  if (xx < par[1]) return par[5];

  Double_t dx = 0.005;
  Double_t xs = par[1];
  Double_t ss = 0.0;
  Double_t paras[2] = {0.0,0.0};
  while (xs < xx) {
    if (xs >= par[1] && xs < (par[1]+par[2])) {
      //fAsymmGauss->SetParameter(0,par[0]);
      //fAsymmGauss->SetParameter(1,xs);
      //ss += fAsymmGauss->Eval(xx);
      paras[0] = par[0];
      paras[1] = xs;
      ss += AsymmGauss(&xx,paras);
    }
    if (xs >= (par[1]+par[2]) && xs < (par[1]+par[2]+par[3])) {
      //fAsymmGauss->SetParameter(0,par[0]*par[4]);
      //fAsymmGauss->SetParameter(1,xs);
      //ss += fAsymmGauss->Eval(xx);
      paras[0] = par[0]*par[4];
      paras[1] = xs;
      ss += AsymmGauss(&xx,paras);
    }
    xs += dx;
  }
  
  return ss+par[5];

}

//_________________________________________________________________________________

Double_t AliTRDCalibra::AsymmGauss(Double_t* x, Double_t* par) 
{
  //
  // Function for the fit
  //

  //par[0] = normalization
  //par[1] = mean
  //par[2] = sigma
  //norm0  = 1
  //par[3] = lambda0
  //par[4] = norm1
  //par[5] = lambda1
  //
  
  Double_t par1save = par[1];    
  //Double_t par2save = par[2];
  Double_t par2save = 0.0604244;
  //Double_t par3save = par[3];
  Double_t par3save = 7.65596;
  //Double_t par5save = par[5];
  Double_t par5save = 0.870597;
  Double_t dx   = x[0]-par1save;
  //
  //
  Double_t  sigma2  = par2save*par2save;
  Double_t  sqrt2   = TMath::Sqrt(2.);
  Double_t  exp1    = par3save*TMath::Exp(-par3save*(dx-0.5*par3save*sigma2))
    *(1-TMath::Erf((par3save*sigma2-dx)/(sqrt2*par2save)));

  Double_t  exp2    = par5save*TMath::Exp(-par5save*(dx-0.5*par5save*sigma2))
    *(1-TMath::Erf((par5save*sigma2-dx)/(sqrt2*par2save)));


  //return par[0]*(exp1+par[4]*exp2);
  return par[0]*(exp1+1.00124*exp2);
}

//__________________________________________________________________________

Double_t AliTRDCalibra::funcLandauGaus(Double_t* x, Double_t* par)
{
  //
  //sum landau + gaus with identical mean
  //

  Double_t valLandau = par[0]*TMath::Landau(x[0],par[1],par[2]);
 //Double_t valGaus   = par[3]*TMath::Gaus(x[0],par[4],par[5]);
  Double_t valGaus   = par[3]*TMath::Gaus(x[0],par[1],par[4]);

  Double_t val = valLandau+valGaus;
  return val;
}
//________________________________________________________________________
Double_t AliTRDCalibra::langaufun(Double_t *x, Double_t *par) 
{
  //
  // Function for the fit
  //
  
  
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}
//______________________________________________________________________________________________
TF1 * AliTRDCalibra::langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  //
  // Function for the fit
  //
  
  
  Int_t i;
  Char_t FunName[100];
  
  AliInfo(Form(FunName,"Fitfcn_%s",his->GetName()));
  
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;
  
  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma");
  
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  
  his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
  
   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
     fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function
   
}
//________________________________________________________________________________________
Int_t AliTRDCalibra::langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
  //
  // Function for the fit
  //
  
  
  Double_t p,x,fy,fxr,fxl;
  Double_t step;
  Double_t l,lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;
  
  
  // Search for maximum
  
  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l    = -1.0;
  
  
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    
    lold = l;
    x = p + step;
    l = langaufun(&x,params);
    
    if (l < lold)
      step = -step/10;
    
    p += step;
  }
  
  if (i == MAXCALLS)
    return (-1);
  
  maxx = x;

  fy = l/2;
  
  
  // Search for right x location of fy
  
  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  
  
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    
    if (l > lold)
      step = -step/10;
 
    p += step;
  }
  
  if (i == MAXCALLS)
    return (-2);
  
  fxr = x;
  
  
  // Search for left x location of fy
  
  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;

    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    
    if (l > lold)
      step = -step/10;
    
    p += step;
  }
  
  if (i == MAXCALLS)
      return (-3);
  

  fxl = x;
  
  FWHM = fxr - fxl;
  return (0);
}


