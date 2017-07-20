#include <signal.h>
//#include <sys/time.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iosfwd>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>

#include <map>
#include <utility>
#include <iterator>


#include <TSystem.h> //To have gSystem known at compilation time
#include <TApplication.h>
#include <TChain.h> 
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h> 
#include <TCanvas.h>
#include <TFrame.h>
#include <TPostScript.h>
#include <TLine.h> 
#include <TGaxis.h> 
#include <TStyle.h> 
#include <TGraphErrors.h> 
#include <TMath.h>
#include <TMatrixF.h>
#include <TText.h>


//#include "/cebaf/faivre/recherche/utilities/defineMyPalette2011.C"


namespace std {} using namespace std;


const double kMaxHV = 395.; // default max voltage limit; could possibly be relaxed in future
const double kMinHV = 210.; // Hexa encoding limit
const double coefFactorWanted=0.0162;


// voir http://dsilverm.web.cern.ch/dsilverm/fee/addrP2.html
 char SMP2Name[][100]={"SMA0","SMC0","SMA1","SMC1","SMA2","SMC2","SMA3","SMC3","SMA4","SMC4","SMA5","SMC5","SMA9","SMC9","SMA10","SMC10","SMA11","SMC11","SMA12","SMC12"};
 char SMcalibName[][100]={"US2","US1","EU2","EU1","US3","US5","US4","EU3","US7","US6","US8C0","US8C2","DCN1","DJP1","DUS2","DUS1","DJP2","DUS3","CN1A","CN1C"};
 char SMnumber[][100]={"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"};

 enum detType {kEMCAL,kEMCALthird,kDCAL,kDCALthird};
 int detTypeType[]={kEMCAL,kEMCALthird,kDCAL,kDCALthird};
 char detTypeString[][100]={"EMCAL","EMCALthird","DCAL","DCALthird"};
 int SMdetType[]={kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCALthird,kEMCALthird,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCALthird,kDCALthird};
 const int kNbColEMCAL=48;
 const int kNbRowEMCAL=24;
 const int kNbSMEMCAL=10;
 const int kNbColEMCALthird=kNbColEMCAL;
 const int kNbRowEMCALthird=(int)(kNbRowEMCAL/3);
 const int kNbSMEMCALthird=2;
 const int kNbColDCAL=32;
 const int kNbRowDCAL=kNbRowEMCAL;
 const int kNbSMDCAL=6;
 const int kNbColDCALthird=kNbColEMCALthird;
 const int kNbRowDCALthird=kNbRowEMCALthird;
 const int kNbSMDCALthird=2;
 const int kNbSMtot=kNbSMEMCAL+kNbSMEMCALthird+kNbSMDCAL+kNbSMDCALthird;
 const int kTabNbCol[4]={kNbColEMCAL,kNbColEMCALthird,kNbColDCAL,kNbColDCALthird};
 const int kTabNbRow[4]={kNbRowEMCAL,kNbRowEMCALthird,kNbRowDCAL,kNbRowDCALthird};
 const int kTabNbSM[4]={kNbSMEMCAL,kNbSMEMCALthird,kNbSMDCAL,kNbSMDCALthird};
 const int kNbColMax=kNbColEMCAL;
 const int kNbRowMax=kNbRowEMCAL;
 const int kNbColOffsetDCAL=kNbColEMCAL-kNbColDCAL;

 int lastSM;




///
/// \file HVrecalculation_EMCALthirds.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief HV recalculation
///
/// HV recalculation for 1/3 SMs, see HVrecalculation.C for full
///
/// How to run :
/// Customize fchNameBadGainIn, then run :
///
/// root -b -q 'macros/HVrecalculation_EMCALthirds.C++("ident")'
///
///
///
/// Stored in the "badGain file" :
/// SM col row choice gain1 gain2 gain3 OCDBcoeff
/// "choice" is 1, 2 or 3 ;
///    gain1 = default value calculated by the program ;
///    gain2 = gain 1.000 (no HV change) ;
///    gain3 = a value to be set in the file by the user.
/// Choice 4 [obsolete] = same as choice 3 (picks up the 3rd gain value) but is used as a tag 
///                       to recognize values copy-pasted by the code from a pass to the next 
///                       one (as opposed to manual setting of the choice by hand in the file).
/// Choices 11, 12, 13 : same as 1, 2, 3, except that they've been copied from the "badGainIn" file.
/// Choice -3 = force take gain3 and OCDBcoeff.
/// Choice 000 = bad HV scan not detected by the code.
///
///
///
/// This code reads HV's for year N-1, OCDB coeffs for year N-1 from pi0 calibration, (p0,p1,p2) parameters from a 
/// recent HV scan, and calculates HV's for year N (to come) with associated OCDB coeffs (most often 0.0162).
///
/// Mapping convention :
/// - All text HVs (input and output)(NewBias.txt files) are in Grenoble cosmic = electronic mapping.
/// - All OCDB files (input and output) are in Alice offline mapping (identical for side A, flipped rows and cols for side C).
/// - The file with the calib coeffs is in Alice offline mapping (identical for side A, flipped rows and cols for side C).
/// - Inside the code, everything is in Grenoble cosmic = electronic mapping, i.e. the OCDB tower coordinates are flipped 
///   at reading and writing, and the calib coeffs are flipped at reading.
/// - Col/row plots ("hSpace") are drawn with Grenoble cosmic = electronic mapping for all SMs.
/// - (p0,p1,p2) parameters file : Grenoble cosmic = electronic mapping.
/// (Caution : flipped doesn't mean col <--> row, but iCol <--> (kNbCol-1)-iCol and iRow <--> (kNbRow-1)-iRow).
///
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///


///
/// Add comment
///
void CalcNewHVperSM(TString hvFileName,TString paramFileName,int smNb,char *fchNameBase,FILE *outFile,FILE *badGainFileOut,double coeffOCDB[kNbSMtot][kNbColMax][kNbRowMax],double coeffCalib[kNbSMtot][kNbColMax][kNbRowMax],double forcedGainChange[kNbSMtot][kNbColMax][kNbRowMax],double forcedOCDBcoeff[kNbSMtot][kNbColMax][kNbRowMax],int flag[kNbSMtot][kNbColMax][kNbRowMax],TPostScript *ps,TCanvas *c1)
{int nbTowersAbove395,nbTowersBelowMin,nbTowersLargeModif,nbTowersBadHVscan,nbTowersProblem;
 int icol,irow,expectSameE;
 double hv,p0,p1,p2,gainChange,newHV,tmpNewHV;
 double HVtab[kNbColMax][kNbRowMax]={{0.}};
 double p0tab[kNbColMax][kNbRowMax]={{0.}};
 double p1tab[kNbColMax][kNbRowMax]={{0.}};
 double p2tab[kNbColMax][kNbRowMax]={{0.}};
 double ampliIn,Eout1,Eout2,Ediff,initialGain,newGain;
 //CUSTOMIZE customize every year :
 double kMaxDeltaGainLook=2.0; //Will printf the tower above this.
 double kMaxDeltaGainAuthorize=3.0; //Will suggest that number to any tower beyond that.
 double kMaxDeltaHVgraph=50.;
 //CUSTOMIZE customize every pass :
 //Pass 0, 1, 2 :
 double kMaxDeltaHVcut=10.;
 
 double kMaxEmismatchCut=0.01; //0.001 is OK but with flag -3 it's too hard. Put 0.01 instead (that's in %).
 
 printf("----- Calculating HV's for SM %d (%s)",smNb,detTypeString[SMdetType[smNb]]);
 fprintf(outFile,"\n-------------- Calculating HV's for SM %d (%s)",smNb,detTypeString[SMdetType[smNb]]);
 if ((smNb%2) == 0) 
    {printf(" (side A)\n");
     fprintf(outFile," (side A) -----\n\n");
     }
   else
    {printf(" (side C)\n");
     fprintf(outFile," (side C) -----\n\n");
     }
 
 // Check file existence and load the table of HV's :
 FILE *inFile,*outHVFile,*outParamFile;
 inFile=fopen(hvFileName.Data(),"r");
 if(!inFile) {printf("File %s can not be found\n",hvFileName.Data());exit(-1);}
 while (fscanf(inFile,"%d %d %lf\n",&icol,&irow,&hv)>0)
    {if (smNb == 10) HVtab[icol][irow]=hv;
     else HVtab[icol][irow-16]=hv;
     }
 fclose(inFile);

 inFile=fopen(paramFileName.Data(),"r");
 if(!inFile) {printf("File %s can not be found\n",paramFileName.Data());exit(-1);}
 while(fscanf(inFile,"%d %d %lf %lf %lf\n",&icol,&irow,&p0,&p1,&p2)>0)
    {if (smNb == 10)
        {p0tab[icol][irow]=p0;
         p1tab[icol][irow]=p1;
         p2tab[icol][irow]=p2;
         }
     else
        {p0tab[icol][irow-16]=p0;
         p1tab[icol][irow-16]=p1;
         p2tab[icol][irow-16]=p2;
         }
     }
 fclose(inFile);

 // Create some histograms :
 TH1F *hEChange = new TH1F(Form("hEChange%d",smNb),Form("hEChange%d",smNb),80,-50.,50.);
 hEChange->SetXTitle("(E_{2}-E_{1})/E_{1} (%)");
 hEChange->SetYTitle("Counts");
 hEChange->SetStats(0);
 TH1F *hEChangeLogP = new TH1F(Form("hEChangeLogP%d",smNb),Form("hEChangeLogP%d",smNb),80,-3.,3.);
 hEChangeLogP->SetXTitle("log [(E_{2}-E_{1})/E_{1} (%)]");
 hEChangeLogP->SetYTitle("Counts");
 hEChangeLogP->SetStats(0);
 TH1F *hEChangeLogN = new TH1F(Form("hEChangeLogN%d",smNb),Form("hEChangeLogN%d",smNb),80,-3.,3.);
 hEChangeLogN->SetXTitle("log [(E_{2}-E_{1})/E_{1} (%)]");
 hEChangeLogN->SetYTitle("Counts");
 hEChangeLogN->SetStats(0);
 TH2F *h2HVdiffVsOCDB = new TH2F(Form("h2HVdiffVsOCDB%d",smNb),Form("h2HVdiffVsOCDB%d",smNb),80,0.,5.,80,-kMaxDeltaHVgraph,kMaxDeltaHVgraph);
 h2HVdiffVsOCDB->SetXTitle("OCDB factor / 0.0162");
 h2HVdiffVsOCDB->SetYTitle("NewHV - OldHV (V)");
 h2HVdiffVsOCDB->SetStats(0);
 h2HVdiffVsOCDB->SetContour(30);
 TH2F *h2HVdiffVsOCDBZm = new TH2F(Form("h2HVdiffVsOCDBZm%d",smNb),Form("h2HVdiffVsOCDBZm%d",smNb),80,0.5,2.,80,-30.,30.);
 h2HVdiffVsOCDBZm->SetXTitle("OCDB factor / 0.0162");
 h2HVdiffVsOCDBZm->SetYTitle("NewHV - OldHV (V)");
 h2HVdiffVsOCDBZm->SetStats(0);
 h2HVdiffVsOCDBZm->SetContour(30);
 TH2F *hSpaceHVbefore = new TH2F(Form("hSpaceHVbefore%d",smNb),Form("hSpaceHVbefore%d",smNb),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
 hSpaceHVbefore->SetXTitle("Column");
 hSpaceHVbefore->SetYTitle("Row");
 hSpaceHVbefore->SetStats(0);
 hSpaceHVbefore->SetContour(30);
 TH2F *hSpaceHVafter = new TH2F(Form("hSpaceHVafter%d",smNb),Form("hSpaceHVafter%d",smNb),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
 hSpaceHVafter->SetXTitle("Column");
 hSpaceHVafter->SetYTitle("Row");
 hSpaceHVafter->SetStats(0);
 hSpaceHVafter->SetContour(30);
 TH2F *hSpaceOCDBbefore = new TH2F(Form("hSpaceOCDBbefore%d",smNb),Form("hSpaceOCDBbefore%d",smNb),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
 hSpaceOCDBbefore->SetXTitle("Column");
 hSpaceOCDBbefore->SetYTitle("Row");
 hSpaceOCDBbefore->SetStats(0);
 hSpaceOCDBbefore->SetContour(30);
 TH2F *hSpaceOCDBafter = new TH2F(Form("hSpaceOCDBafter%d",smNb),Form("hSpaceOCDBafter%d",smNb),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
 hSpaceOCDBafter->SetXTitle("Column");
 hSpaceOCDBafter->SetYTitle("Row");
 hSpaceOCDBafter->SetStats(0);
 hSpaceOCDBafter->SetContour(30);
 TH2F *hSpaceAboveBefore = new TH2F(Form("hSpaceAboveBefore%d",smNb),Form("hSpaceAboveBefore%d",smNb),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
 hSpaceAboveBefore->SetXTitle("Column");
 hSpaceAboveBefore->SetYTitle("Row");
 hSpaceAboveBefore->SetStats(0);
 hSpaceAboveBefore->SetContour(30);
 TH2F *hSpaceAboveAfter = new TH2F(Form("hSpaceAboveAfter%d",smNb),Form("hSpaceAboveAfter%d",smNb),kNbColMax,-0.5,kNbColMax-0.5,kNbRowMax,-0.5,kNbRowMax-0.5);
 hSpaceAboveAfter->SetXTitle("Column");
 hSpaceAboveAfter->SetYTitle("Row");
 hSpaceAboveAfter->SetStats(0);
 hSpaceAboveAfter->SetContour(30);
 
 // Setup a few counters :
 ampliIn=1.; //Nb of "photoelectrons" (arbitrary unit)
 nbTowersAbove395=0;
 nbTowersBelowMin=0;
 nbTowersLargeModif=0;
 nbTowersBadHVscan=0;
 nbTowersProblem=0;
 
 // Loop over the towers to calculate the new HV's :
 for(icol=0;icol<kTabNbCol[SMdetType[smNb]];icol++)
    {for(irow=0;irow<kTabNbRow[SMdetType[smNb]];irow++)
 //for(icol=4;icol<5;icol++)
 //   {for(irow=3;irow<4;irow++)
        {hv=HVtab[icol][irow];
         p0=p0tab[icol][irow];
         p1=p1tab[icol][irow];
         p2=p2tab[icol][irow];
//printf("%d %d %d    HV %f   coeff %f %f\n",smNb,icol,irow,hv,coeffCalib[smNb][icol][irow],0.0162/coeffOCDB[smNb][icol][irow]);
         expectSameE=1;
         hSpaceHVbefore->Fill(icol,irow,hv);
         hSpaceOCDBbefore->Fill(icol,irow,coeffOCDB[smNb][icol][irow]/coefFactorWanted);
         if (hv == kMaxHV) hSpaceAboveBefore->Fill(icol,irow,coeffOCDB[smNb][icol][irow]/coefFactorWanted);
         
         gainChange=coeffOCDB[smNb][icol][irow]/coefFactorWanted;
         // Suggest discarding the towers with calib coeff = 1 (half of the time dead towers) :
         if ((coeffCalib[smNb][icol][irow] == 1.0) && (flag[smNb][icol][irow] == 5)) //Don't escape again if tower has already seen at previous pass !
            {fprintf(outFile,"######## Tower with calib coeff 1.0 (%d,%d,%d) ### do nothing.\n",smNb,icol,irow);
             fprintf(badGainFileOut,"%d %2d %2d XXX2 %f %4.2f %f %f\n",smNb,icol,irow,gainChange,1.,forcedGainChange[smNb][icol][irow],forcedOCDBcoeff[smNb][icol][irow]);
             nbTowersProblem++;
             continue;
             }
         
         // Check that the p0, p1, p2 parameters look OK (raw check) :
         if ((p1==0.) && (p2==0.)) // LED HV-scan not good.
            {//fprintf(outFile,"######## Bad HV scan for tower (%d,%d,%d), HV %7.3f  ### do nothing.\n",smNb,icol,irow,hv);
             //For 1/3-EMCAL-SMs : only this error appears with the HV scans. Calculate new HV with default parameters :
             fprintf(outFile,"######## Bad HV scan for tower (%d,%d,%d), HV %7.3f  ### use default parameters.\n",smNb,icol,irow,hv);
             fprintf(outFile,"######## CAUTION : EMCAL thirds 2015 only !!!!!!\n");
             p0=0.;
             p1=1.;
             p2=0.02;
             nbTowersBadHVscan++;
             nbTowersProblem++;
             //continue;
             }
         if ((p1<=0.) || (p2<=0.))
            {fprintf(outFile,"######## Bad HV scan parameters for tower (%d,%d,%d), HV %7.3f (params : %f %e %f) ### do nothing.\n",smNb,icol,irow,hv,p0,p1,p2);
             nbTowersBadHVscan++;
             nbTowersProblem++;
             continue;
             }
         if ( (p0 + p1 * exp(p2*(hv-kMaxDeltaHVcut-10.))) < 0. )
            {fprintf(outFile,"######## Bad HV scan value for tower (%d,%d,%d), HV %7.3f (params : %f %e %f) ### do nothing.\n",smNb,icol,irow,hv,p0,p1,p2);
             nbTowersBadHVscan++;
             nbTowersProblem++;
             continue;
             }
         
         // Check for a chosen, forced value.
         // If there isn't, check that the required gain change isn't too far from 1.0.
         if (flag[smNb][icol][irow] == 0) //Undecided tower (shouldn't be found...), copy as is.
            {fprintf(badGainFileOut,"%d %2d %2d XXX %f %4.2f %f %f\n",smNb,icol,irow,gainChange,1.,forcedGainChange[smNb][icol][irow],forcedOCDBcoeff[smNb][icol][irow]);
             nbTowersProblem++;
             continue;
             }
         if (flag[smNb][icol][irow] == -3) //HV change was too large, forced gain and OCDB value.
            {fprintf(badGainFileOut,"%d %2d %2d -3 %f %4.2f %f %f\n",smNb,icol,irow,gainChange,1.,forcedGainChange[smNb][icol][irow],forcedOCDBcoeff[smNb][icol][irow]);
             gainChange=forcedGainChange[smNb][icol][irow];
             initialGain = p0 + p1 * exp(p2*hv);
             newGain = initialGain * gainChange;
             Eout1=ampliIn*initialGain*coeffOCDB[smNb][icol][irow];
             //expectSameE=0; //If the "-3" results from a too large HV change (most of the cases), the energies must be equal ! Some outliers may be caught by the printf.
             if ((newGain-p0) > 0)
                {newHV = log ( (newGain - p0)/p1 ) / p2;
                 HVtab[icol][irow]=newHV;
                 coeffOCDB[smNb][icol][irow]=forcedOCDBcoeff[smNb][icol][irow];
                 Eout2=ampliIn*(p0 + p1 * exp(p2*newHV))*coeffOCDB[smNb][icol][irow];
                 h2HVdiffVsOCDB->Fill(gainChange,newHV-hv);
                 h2HVdiffVsOCDBZm->Fill(gainChange,newHV-hv);
                 Ediff=100.*(Eout2-Eout1)/Eout1;
                 if (expectSameE == 1)
                    {hEChange->Fill(Ediff);
                     if (Eout2 > Eout1) hEChangeLogP->Fill(TMath::Log10(TMath::Abs(Ediff)));
                     if (Eout2 < Eout1) hEChangeLogN->Fill(TMath::Log10(TMath::Abs(Ediff)));
                     if (TMath::Abs(Ediff) > kMaxEmismatchCut)
                        {printf("##### E mismatch tower (%d,%d,%d) : Eout1 = %f, Eout2 = %f     HV %7.3f -> %7.3f    =>  %7.3f w/ OCDB %f   flag %d\n",smNb,icol,irow,Eout1,Eout2,hv,newHV,HVtab[icol][irow],coeffOCDB[smNb][icol][irow],flag[smNb][icol][irow]);
                         fprintf(outFile,"E mismatch tower (%d,%d,%d) : Eout1 = %f, Eout2 = %f     HV %7.3f -> %7.3f    =>  %7.3f w/ OCDB %f   flag %d\n",smNb,icol,irow,Eout1,Eout2,hv,newHV,HVtab[icol][irow],coeffOCDB[smNb][icol][irow],flag[smNb][icol][irow]);
                         }
                     }
                 continue;
                 }
                else // LED HV-scan not good (can happen if changed HV scan between an execution and the next one) : tell it and move to the next tower.
                {fprintf(outFile,"######## Negative in log for forced tower (%d,%d,%d), HV %7.3f and gain %5.3f (params : %f %e %f) ### do nothing.\n",smNb,icol,irow,hv,gainChange,p0,p1,p2);
                 continue;
                 }
             }
         if (forcedGainChange[smNb][icol][irow] != 0.)
            {if ((flag[smNb][icol][irow] > 9) || (flag[smNb][icol][irow] == -4)) fprintf(badGainFileOut,"%d %2d %2d %d %f %4.2f %f %f\n",smNb,icol,irow,flag[smNb][icol][irow],gainChange,1.,forcedGainChange[smNb][icol][irow],forcedOCDBcoeff[smNb][icol][irow]); //Can't put these lines later because the older "gainChange" must be known. That causes this line to be issued a 2nd time (with different values) in badGainFileOut when deltaHV is too large. The latter line overwrites the previous anyway when the file is read back => that's OK even if the user doesn't cancel it manually. The line automatically disappears from the "bad gain" file at the next pass.
                else fprintf(badGainFileOut,"%d %2d %2d %d %f %4.2f %f %f\n",smNb,icol,irow,flag[smNb][icol][irow]+10,gainChange,1.,forcedGainChange[smNb][icol][irow],forcedOCDBcoeff[smNb][icol][irow]);
             gainChange=forcedGainChange[smNb][icol][irow];
             }
            else
            {if (gainChange < 1./kMaxDeltaGainAuthorize)
                {fprintf(outFile,"######## Gain change < factor 1/%4.2f for tower (%d,%d,%d), HV %7.3f  and gain %5.3f ### do nothing.\n",kMaxDeltaGainAuthorize,smNb,icol,irow,hv,gainChange);
                 fprintf(badGainFileOut,"%d %2d %2d XXX-4 %f %4.2f %f %f\n",smNb,icol,irow,gainChange,1.,1./kMaxDeltaGainAuthorize,coefFactorWanted*1./kMaxDeltaGainAuthorize); //Have to setup both the gain change and the OCDB factor, as *that* OCDB factor is still the desired change (not yet the factor that will be written to file).
                 nbTowersLargeModif++;
                 nbTowersProblem++;
                 continue;
                 }
             if (gainChange > kMaxDeltaGainAuthorize)
                {fprintf(outFile,"######## Gain change > factor %4.2f for tower (%d,%d,%d), HV %7.3f  and gain %5.3f ### do nothing.\n",kMaxDeltaGainAuthorize,smNb,icol,irow,hv,gainChange);
                 fprintf(badGainFileOut,"%d %2d %2d XXX-4 %f %4.2f %f %f\n",smNb,icol,irow,gainChange,1.,kMaxDeltaGainAuthorize,coefFactorWanted*kMaxDeltaGainAuthorize); //Have to setup both the gain change and the OCDB factor, as *that* OCDB factor is still the desired change (not yet the factor that will be written to file).
                 nbTowersLargeModif++;
                 nbTowersProblem++;
                 continue;
                 }
             if ((gainChange < 1./kMaxDeltaGainLook) || (gainChange > kMaxDeltaGainLook))
                {fprintf(outFile,"######## Gain change > factor %4.2f for tower (%d,%d,%d), HV %7.3f  and gain %5.3f ### do nothing.\n",kMaxDeltaGainLook,smNb,icol,irow,hv,gainChange);
                 fprintf(badGainFileOut,"%d %2d %2d XXX1 %f %4.2f %4.2f %4.2f\n",smNb,icol,irow,gainChange,1.,0.,0.);
                 nbTowersLargeModif++;
                 nbTowersProblem++;
                 continue;
                 }
             }
         if ((flag[smNb][icol][irow] == 2) || (flag[smNb][icol][irow] == 12) || (flag[smNb][icol][irow] == 3) || (flag[smNb][icol][irow] == 13)) expectSameE=0;
         
         // Calculate the new HV :
         newHV=hv;
         initialGain = p0 + p1 * exp(p2*hv);
         newGain = initialGain * gainChange;
         Eout1=ampliIn*initialGain*coeffOCDB[smNb][icol][irow];
         if ((newGain-p0) <= 0) // LED HV-scan not good.
            {fprintf(outFile,"######## Negative in log for tower (%d,%d,%d), HV %7.3f and gain %5.3f (params : %f %e %f) ### do nothing.\n",smNb,icol,irow,hv,gainChange,p0,p1,p2);
             nbTowersBadHVscan++;
             nbTowersProblem++;
             continue;
             }
             else
             newHV = log ( (newGain - p0)/p1 ) / p2;
         if (flag[smNb][icol][irow] == -4) coeffOCDB[smNb][icol][irow]=forcedOCDBcoeff[smNb][icol][irow]; //GainChange has been levelled off at kMaxDeltaGainAuthorize
         
         // Check the calculated newHV :
         if (newHV < 0) // conversion failed:  let's just keep the old custom value then
            {fprintf(outFile,"######## Neg HV returned for tower (%d,%d,%d), HV %7.3f and gain %5.3f (params : %f %e %f) ### keeping previous value.\n",smNb,icol,irow,hv,gainChange,p0,p1,p2);
             nbTowersProblem++;
             continue;
             }
         if (TMath::Abs(newHV-hv) > kMaxDeltaHVcut)
            {fprintf(outFile,"######## Too large HV change for tower (%d,%d,%d), HV %7.3f -> %7.3f with gain %5.3f (params : %f %e %f) ### do nothing, ",smNb,icol,irow,hv,newHV,gainChange,p0,p1,p2);
             fprintf(badGainFileOut,"%d %2d %2d XXX-3 %f %4.2f ",smNb,icol,irow,gainChange,1.);
             if (newHV > hv) tmpNewHV=TMath::Min(kMaxHV,hv+(kMaxDeltaHVcut-0.01)); //Use val-0.01 instead of val so that we're sure it doesn't get trapped again after correcting the gain.
                else tmpNewHV=TMath::Max(kMinHV,hv-(kMaxDeltaHVcut-0.01));
             newGain= p0 + p1 * exp(p2*tmpNewHV);
             gainChange=newGain/initialGain;
             fprintf(outFile,"suggest gain %f w/ OCDBcoeff %f\n",gainChange,coeffOCDB[smNb][icol][irow]/gainChange);
             fprintf(badGainFileOut,"%f %f\n",gainChange,coeffOCDB[smNb][icol][irow]/gainChange);
             nbTowersProblem++;
             continue;
             }
         if (newHV < kMinHV)
            {fprintf(outFile,"HV below min voltage for tower (%d,%d,%d), HV %7.3f -> %7.3f with gain %5.3f -- set to supply voltage %7.3f.\n",smNb,icol,irow,hv,newHV,gainChange,kMinHV);
             nbTowersBelowMin++;
             nbTowersProblem++;
             newGain= p0 + p1 * exp(p2*kMinHV);
             gainChange=newGain/initialGain;
             coeffOCDB[smNb][icol][irow]=coeffOCDB[smNb][icol][irow]/gainChange;
             HVtab[icol][irow]=kMinHV;
             Eout2=ampliIn*newGain*coeffOCDB[smNb][icol][irow];
             Ediff=100.*(Eout2-Eout1)/Eout1;
             if (expectSameE == 1) //Energies are not equal when we force the gain change to be different from the calculated OCDB coeff.
                {hEChange->Fill(Ediff);
                 if (Eout2 > Eout1) hEChangeLogP->Fill(TMath::Log10(TMath::Abs(Ediff)));
                 if (Eout2 < Eout1) hEChangeLogN->Fill(TMath::Log10(TMath::Abs(Ediff)));
                 if (TMath::Abs(Ediff) > kMaxEmismatchCut)
                    {printf("##### E mismatch tower (%d,%d,%d) : Eout1 = %f, Eout2 = %f     HV %7.3f -> %7.3f    =>  %7.3f w/ OCDB %f   flag %d\n",smNb,icol,irow,Eout1,Eout2,hv,newHV,HVtab[icol][irow],coeffOCDB[smNb][icol][irow],flag[smNb][icol][irow]);
                     fprintf(outFile,"E mismatch tower (%d,%d,%d) : Eout1 = %f, Eout2 = %f     HV %7.3f -> %7.3f    =>  %7.3f w/ OCDB %f   flag %d\n",smNb,icol,irow,Eout1,Eout2,hv,newHV,HVtab[icol][irow],coeffOCDB[smNb][icol][irow],flag[smNb][icol][irow]);
                     }
                 }
             continue;
             }
         if (newHV>kMaxHV) // we reached a too high voltage: let's keep the max then
            {fprintf(outFile,"HV above supply voltage for tower (%d,%d,%d), HV %7.3f -> %7.3f with gain %5.3f -- set to supply voltage %7.3f.\n",smNb,icol,irow,hv,newHV,gainChange,kMaxHV);
             nbTowersAbove395++;
             nbTowersProblem++;
             newGain= p0 + p1 * exp(p2*kMaxHV);
             gainChange=newGain/initialGain;
             coeffOCDB[smNb][icol][irow]=coeffOCDB[smNb][icol][irow]/gainChange;
             HVtab[icol][irow]=kMaxHV;
             Eout2=ampliIn*newGain*coeffOCDB[smNb][icol][irow];
             Ediff=100.*(Eout2-Eout1)/Eout1;
             if (expectSameE == 1) //Energies are not equal when we force the gain change to be different from the calculated OCDB coeff.
                {hEChange->Fill(Ediff);
                 if (Eout2 > Eout1) hEChangeLogP->Fill(TMath::Log10(TMath::Abs(Ediff)));
                 if (Eout2 < Eout1) hEChangeLogN->Fill(TMath::Log10(TMath::Abs(Ediff)));
                 if (TMath::Abs(Ediff) > kMaxEmismatchCut)
                    {printf("##### E mismatch tower (%d,%d,%d) : Eout1 = %f, Eout2 = %f     HV %7.3f -> %7.3f    =>  %7.3f w/ OCDB %f   flag %d\n",smNb,icol,irow,Eout1,Eout2,hv,newHV,HVtab[icol][irow],coeffOCDB[smNb][icol][irow],flag[smNb][icol][irow]);
                     fprintf(outFile,"E mismatch tower (%d,%d,%d) : Eout1 = %f, Eout2 = %f     HV %7.3f -> %7.3f    =>  %7.3f w/ OCDB %f   flag %d\n",smNb,icol,irow,Eout1,Eout2,hv,newHV,HVtab[icol][irow],coeffOCDB[smNb][icol][irow],flag[smNb][icol][irow]);
                     }
                 }
             continue;
             }
         //everything ok
         HVtab[icol][irow]=newHV;
         coeffOCDB[smNb][icol][irow]=coefFactorWanted;
         Eout2=ampliIn*(p0 + p1 * exp(p2*newHV))*coeffOCDB[smNb][icol][irow];
//if ((gainChange < 0.95) && (newHV-hv > -1)) printf("##### BAAAAD tower (%d,%d,%d) HV's %f %f, gainChange %f, p's %e %e %e\n",smNb,icol,irow,hv,newHV,gainChange,p0,p1,p2); /////////
         h2HVdiffVsOCDB->Fill(gainChange,newHV-hv);
         h2HVdiffVsOCDBZm->Fill(gainChange,newHV-hv);
         Ediff=100.*(Eout2-Eout1)/Eout1;
         if (expectSameE == 1) //Energies are not equal when we force the gain change to be different from the calculated OCDB coeff.
            {hEChange->Fill(Ediff);
             if (Eout2 > Eout1) hEChangeLogP->Fill(TMath::Log10(TMath::Abs(Ediff)));
             if (Eout2 < Eout1) hEChangeLogN->Fill(TMath::Log10(TMath::Abs(Ediff)));
             if (TMath::Abs(Ediff) > kMaxEmismatchCut)
                {printf("##### E mismatch tower (%d,%d,%d) : Eout1 = %f, Eout2 = %f     HV %7.3f -> %7.3f    =>  %7.3f w/ OCDB %f   flag %d\n",smNb,icol,irow,Eout1,Eout2,hv,newHV,HVtab[icol][irow],coeffOCDB[smNb][icol][irow],flag[smNb][icol][irow]);
                 fprintf(outFile,"E mismatch tower (%d,%d,%d) : Eout1 = %f, Eout2 = %f     HV %7.3f -> %7.3f    =>  %7.3f w/ OCDB %f   flag %d\n",smNb,icol,irow,Eout1,Eout2,hv,newHV,HVtab[icol][irow],coeffOCDB[smNb][icol][irow],flag[smNb][icol][irow]);
                 }
             }
         }
     }

 // Write to file the new HV values and the new OCDB parameters :
 outHVFile=fopen(Form("%s/%s/NewBias.txt",fchNameBase,SMP2Name[smNb]),"w");
 outParamFile=fopen(Form("%s/%s/NewParamFactor.txt",fchNameBase,SMP2Name[smNb]),"w");
 for(icol=0;icol<kTabNbCol[SMdetType[smNb]];icol++)
    {for(irow=0;irow<kTabNbRow[SMdetType[smNb]];irow++)
        {fprintf(outHVFile,"%d %d %f\n",icol,irow,HVtab[icol][irow]);
         //Flip OCDB params from Grenoble/electonic mapping to Alice-offline mapping. HVs are kept in Grenoble/electonic mapping.
         if ((smNb%2) == 0) fprintf(outParamFile,"%d %d %d %f\n",smNb,icol,irow,coeffOCDB[smNb][icol][irow]); //Side A.
         else fprintf(outParamFile,"%d %d %d %f\n",smNb,(kTabNbCol[SMdetType[smNb]]-1)-icol,(kTabNbRow[SMdetType[smNb]]-1)-irow,coeffOCDB[smNb][icol][irow]); //Side C.
         hSpaceHVafter->Fill(icol,irow,HVtab[icol][irow]);
         hSpaceOCDBafter->Fill(icol,irow,coeffOCDB[smNb][icol][irow]/coefFactorWanted);
         if (HVtab[icol][irow] == kMaxHV) hSpaceAboveAfter->Fill(icol,irow,coeffOCDB[smNb][icol][irow]/coefFactorWanted);
         }
     }
 fclose(outHVFile);
 fclose(outParamFile);
 
 fprintf(outFile,"\n\nNb towers above 395 V :       %d\n",nbTowersAbove395);
 fprintf(outFile,"Nb towers below min voltage : %d\n",nbTowersBelowMin);
 fprintf(outFile,"Nb towers w/ too large gain : %d <--\n",nbTowersLargeModif);
 fprintf(outFile,"Nb towers w/ bad HV scan :    %d <--\n",nbTowersBadHVscan);
 fprintf(outFile,"Nb other towers to look at :  %d <--\n",nbTowersProblem-nbTowersAbove395-nbTowersBelowMin-nbTowersLargeModif-nbTowersBadHVscan);
 fprintf(outFile,"Tot nb towers w/ problem :    %d\n",nbTowersProblem);
 fprintf(outFile,"\n");

 // Draw a few plots :
 TLine *lineH = new TLine(h2HVdiffVsOCDBZm->GetXaxis()->GetXmin(),0.,h2HVdiffVsOCDBZm->GetXaxis()->GetXmax(),0.);
 TLine *lineV = new TLine(1.,h2HVdiffVsOCDBZm->GetYaxis()->GetXmin(),1.,h2HVdiffVsOCDBZm->GetYaxis()->GetXmax());
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(1);
 hEChange->Draw();
 gPad->SetLogy();
 hEChange->Write();
 c1->cd(3);
 hEChangeLogP->Draw();
 gPad->SetLogy();
 hEChangeLogP->Write();
 c1->cd(5);
 hEChangeLogN->Draw();
 gPad->SetLogy();
 hEChangeLogN->Write();
 c1->cd(2);
 h2HVdiffVsOCDB->Draw("COLZ");
 h2HVdiffVsOCDB->Write();
 lineH->Draw();
 lineV->Draw();
 c1->cd(4);
 h2HVdiffVsOCDBZm->Draw("COLZ");
 h2HVdiffVsOCDBZm->Write();
 lineH->Draw();
 lineV->Draw();
 c1->Update(); 
 
 ps->NewPage();
 c1->Clear();
 c1->Divide(2,3);
 c1->cd(1);
 hSpaceHVbefore->SetMinimum(kMinHV);
 hSpaceHVbefore->SetMaximum(kMaxHV+5.);
 hSpaceHVbefore->Draw("COLZ");
 hSpaceHVbefore->Write();
 c1->cd(2);
 hSpaceHVafter->SetMinimum(kMinHV);
 hSpaceHVafter->SetMaximum(kMaxHV+5.);
 hSpaceHVafter->Draw("COLZ");
 hSpaceHVafter->Write();
 c1->cd(3);
 hSpaceOCDBbefore->Draw("COLZ");
 hSpaceOCDBbefore->Write();
 c1->cd(4);
 hSpaceOCDBafter->Draw("COLZ");
 hSpaceOCDBafter->Write();
 c1->cd(5);
 hSpaceAboveBefore->Draw("COLZ");
 hSpaceAboveBefore->Write();
 c1->cd(6);
 hSpaceAboveAfter->Draw("COLZ");
 hSpaceAboveAfter->Write();
 c1->Update();
 
 if (smNb < (lastSM-1)) //Otherwise, the plots on the last page don't show up. Not sure why.
    {delete hEChange;
     delete hEChangeLogP;
     delete hEChangeLogN;
     delete h2HVdiffVsOCDB;
     delete h2HVdiffVsOCDBZm;
     delete hSpaceHVbefore;
     delete hSpaceHVafter;
     delete hSpaceOCDBbefore;
     delete hSpaceOCDBafter;
     delete hSpaceAboveBefore;
     delete hSpaceAboveAfter;
     } 

 return;
 }





///
/// Main method
///
/// \param choiceArg indicates which SM are used :
///  10^0 = EMCAL, 10^1 = EMCAL thirds, 10^2 = DCAL, 10^3 = DCAL thirds,
///  0 = no, 1 = yes.
/// Therefore : old EMCAL SMs only -> choiceArg = 0b0001 ;
///             thirds only -> choiceArg = 0b1010 ;
///             DCAL only -> choiceArg = 0b0100.
/// Caution, do not write eg '0001', because integers with a leading zero in C are interpreted as octal values.
/// Instead, write '0b0001' to state it's binary.
///
void HVrecalculation_EMCALthirds(TString ident,int choiceArg=0b0010)
{int i,ism,icol,irow,choice,oldFlag;
 int tabChoiceCalos[4],tmpChoice;
 double paramFactor;
 double val1,val2,val3,ocdbCoeff;
 
  //defineMyPalette2011(30,5);
 
 for (i=4-1;i>=0;i--)
    {tmpChoice=choiceArg>>i;
     tabChoiceCalos[i]=tmpChoice;
     choiceArg-=tmpChoice*pow(2,i);
     }
 printf("\n\n---------------------------------\n| Running for :");
 for (i=0;i<4;i++)
    {if (tabChoiceCalos[i] == 1) printf(" %s",detTypeString[i]);
     }
 printf("\n");
 if (tabChoiceCalos[0] == 1) lastSM=kNbSMEMCAL;
 if (tabChoiceCalos[1] == 1) lastSM=kNbSMEMCAL+kNbSMEMCALthird;
 if (tabChoiceCalos[2] == 1) lastSM=kNbSMEMCAL+kNbSMEMCALthird+kNbSMDCAL;
 if (tabChoiceCalos[3] == 1) lastSM=kNbSMtot;
 
 // Setup filenames, including :
 //    - new directory created to store everything,
 //    - name of the file with the OCDB parameters,
 //    - name of the file (to be created) with the "bad [large] gain towers",
 //    - name of the file (**to provide**) with the wanted gain for towers with initially "bad [large] gain".
 char fchNameBase[250],fchNameRoot[250],fchNamePs[250],fchNameOut[250],fchNameBadGainIn[250],fchNameBadGainOut[250];
 FILE *inFile,*outFile,*badGainFileIn,*badGainFileOut,*OCDBcoeffsFile,*calibCoeffsFile;
 sprintf(fchNameBase,"output_HVrecalculation_EMCALthirds_%s",ident.Data());
 sprintf(fchNameRoot,"%s/%s.root",fchNameBase,fchNameBase);
 sprintf(fchNamePs,"%s/%s.ps",fchNameBase,fchNameBase);
 sprintf(fchNameOut,"%s/%s.out",fchNameBase,fchNameBase);
 sprintf(fchNameBadGainOut,"%s/%s_badGain.txt",fchNameBase,fchNameBase);
 gSystem->Exec(Form("mkdir %s",fchNameBase));
 for (ism=0;ism<kNbSMtot;ism++) gSystem->Exec(Form("mkdir %s/%s",fchNameBase,SMP2Name[ism]));
 outFile=fopen(fchNameOut,"w");

 inFile=fopen(fchNameBadGainOut,"r");
 if (inFile)
    {printf("File %s already exists, exiting.\n",fchNameBadGainOut);
     fclose(inFile); //Must be inside the if-loop, otherwise breaks (tries to close a file that it didn't open because it didn't exist).
     return;
     }

 TFile *rootFileOut = new TFile(fchNameRoot,"RECREATE");
 const int cWidth=500;
 const int cHeight=(int)(500*(29./21.));
 TCanvas *c1 = new TCanvas("c1","EMCal cosmics analysis",cWidth,cHeight);
 TPostScript *ps = new TPostScript(fchNamePs,111);
 
 //CUSTOMIZE customize every year (OCDB coeffs + calib coeffs) :
 const char fchNameOCDBcoeffsFile[] = "/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/OCDBparamsAfterCalib2015.txt";
 const char fchNameCalibCoeffsFile[] = "/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/output_calibPi0_coeffs_clean_finalFile4HVcalculation_setUntrustedToOne.txt"; //This is used only to check which towers are dead (calib coeff = 1).
 
 //CUSTOMIZE customize at every pass :
 //***** To be set by the user at every execution : *****//
 //For pass0 :
 //sprintf(fchNameBadGainIn,"000.txt");
 //For pass1 :
 //sprintf(fchNameBadGainIn,"/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/output_HVrecalculation_EMCALthirds_pass0/output_HVrecalculation_EMCALthirds_pass0_badGain_Custom.txt");
 //For pass2 :
 sprintf(fchNameBadGainIn,"/cebaf/cebaf/EMCAL/calibPi0_run2/recalculateHV_4_with2015data/output_HVrecalculation_EMCALthirds_pass1/output_HVrecalculation_EMCALthirds_pass1_badGain.txt"); //Nothing to customize, all the lines were OK.

 printf("\n***** Will be using 'bad gain' correction file %s. PLEASE CHECK.\n",fchNameBadGainIn);
 fprintf(outFile,"\nUse 'bad gain' correction file %s.\n",fchNameBadGainIn);
 printf("Number of lines with 'XXX' :\n");
 fprintf(outFile,"Number of lines with 'XXX' :\n");
 gSystem->Exec(Form("grep -c XXX %s",fchNameBadGainIn));
 gSystem->Exec(Form("grep -c XXX %s >> %s",fchNameBadGainIn,fchNameOut));
 printf("\n\n");
 fprintf(outFile,"\n\n");
 fprintf(outFile,"Tower numbers are given in Grenoble/electronic convention.\n\n");

 // Load the calibration parameters from the OCDB text file :
 double coeffOCDB[kNbSMtot][kNbColMax][kNbRowMax]={{{0.}}};
 double coeffCalib[kNbSMtot][kNbColMax][kNbRowMax]={{{0.}}};
 OCDBcoeffsFile = fopen(fchNameOCDBcoeffsFile,"r");
 if(!OCDBcoeffsFile) {printf("File %s can not be found\n",fchNameOCDBcoeffsFile);exit(-1);}
 while (fscanf(OCDBcoeffsFile,"%d %d %d %lf",&ism,&icol,&irow,&paramFactor)>0) //Flip from Alice-offline mapping to Grenoble/electronic mapping.
    {if ((ism%2) == 0) coeffOCDB[ism][icol][irow]=paramFactor; //Side A.
     else coeffOCDB[ism][(kTabNbCol[SMdetType[ism]]-1)-icol][(kTabNbRow[SMdetType[ism]]-1)-irow]=paramFactor; //Side C.
     }
 fclose(OCDBcoeffsFile);
 // Load the calib coefficients used to produce the OCDB file :
 calibCoeffsFile = fopen(fchNameCalibCoeffsFile,"r");
 if(!calibCoeffsFile) {printf("File %s can not be found\n",fchNameCalibCoeffsFile);exit(-1);}
 while (fscanf(calibCoeffsFile,"%d %d %d %d %lf",&oldFlag,&ism,&icol,&irow,&paramFactor)>0) //Flip from Alice-offline mapping to Grenoble/electronic mapping.
    {if ((ism%2) == 0) coeffCalib[ism][icol][irow]=paramFactor; //Side A.
     else coeffCalib[ism][(kTabNbCol[SMdetType[ism]]-1)-icol][(kTabNbRow[SMdetType[ism]]-1)-irow]=paramFactor; //Side C.
     }
 fclose(calibCoeffsFile);

 // Read bad gain correction file and store the desired gains :
 printf("\n\nRead 'bad gain' correction file.\n");
 double forcedGainChange[kNbSMtot][kNbColMax][kNbRowMax]={{{0.}}};
 double forcedOCDBcoeff[kNbSMtot][kNbColMax][kNbRowMax]={{{0.}}};
 int flag[kNbSMtot][kNbColMax][kNbRowMax]={{{0}}};
 //Need to initialize by hand because {{{5}}} doesn't put "5" everywhere ; it's still initialized to 0 !
 for (ism=0;ism<kNbSMtot;ism++)
    {for (icol=0;icol<kNbColMax;icol++)
        {for (irow=0;irow<kNbRowMax;irow++) flag[ism][icol][irow]=5;
         }
     }
 //In table "flag" :
 // -3 -> force gain *and* OCDB coeff,
 // 0 -> 'bad gain' tower but not customized yet,
 // 1 -> calculated OCDB coeff OK,
 // 2 -> don't change HV,
 // 3 -> apply hand-calculated HV change,
 // 4 [obsolete] -> 'bad gain' tower already looked at in a previous pass,
 // 11,12,13 -> 'bad gain' tower already looked at in a previous pass, with flag resp. 1, 2, 3.
 // 5 -> not a 'bad gain' tower.
 badGainFileIn=fopen(fchNameBadGainIn,"r");
 if (!badGainFileIn)
    {printf("##### 'Bad gain' correction file %s can not be found.\n##### CAUTION : is this the 1st pass ? #####\n\n",fchNameBadGainIn);
     fprintf(outFile,"##### 'Bad gain' correction file %s can not be found.\n##### CAUTION : is this the 1st pass ? #####\n\n",fchNameBadGainIn);
     }
    else
    {while (fscanf(badGainFileIn,"%d %d %d %d %lf %lf %lf %lf",&ism,&icol,&irow,&choice,&val1,&val2,&val3,&ocdbCoeff)>0) //Caution, do use %lf with doubles, not %f.
        {switch (TMath::Abs(choice))
            {case 1 : if (val3 != 0.0) printf("#### Check 'bad gain' correction file (%d,%d,%d) : choice %d while val3 %f (should be 0.)\n",ism,icol,irow,choice,val3);
             case 11 : forcedGainChange[ism][icol][irow]=val1;
                       if (val2 != 1.0) printf("#### Check 'bad gain' correction file (%d,%d,%d) : choice %d while val2 %f (should be 1.)\n",ism,icol,irow,choice,val2);
                       if (ocdbCoeff != 0.0) printf("#### Check 'bad gain' correction file (%d,%d,%d) : choice %d while OCDB coeff %f (should be 0.)\n",ism,icol,irow,choice,ocdbCoeff);
                       break;
             case 2 : 
             case 12 : forcedGainChange[ism][icol][irow]=val2;
                       break;
             case 0 : 
             case 3 : 
             case 4 : 
             case 13 : forcedGainChange[ism][icol][irow]=val3;
                       if ((choice == -3) && (ocdbCoeff == 0.)) printf("#### Check 'bad gain' correction file (%d,%d,%d) : choice %d while OCDB coeff %f.\n",ism,icol,irow,choice,ocdbCoeff);
                       if (val3 == 0.) printf("#### Check 'bad gain' correction file (%d,%d,%d) : choice %d while val3 %f.\n",ism,icol,irow,choice,val3);
                       break;
             default : printf("*** Unknown choice %d for (SM,col,row) (%d,%d,%d)\n",choice,ism,icol,irow);
                       fprintf(outFile,"*** Unknown choice %d for (SM,col,row) (%d,%d,%d)\n",choice,ism,icol,irow);
             }
         forcedOCDBcoeff[ism][icol][irow]=ocdbCoeff;
         flag[ism][icol][irow]=choice;
         }
     fclose(badGainFileIn);
     }
 printf("\n");
 
 // Launch subroutine to calculate new HV's and new OCDB params for every SM :
 badGainFileOut=fopen(fchNameBadGainOut,"w");
 //CUSTOMIZE customize every year (HV files + LED scans) :
 //Merging of HV scans LED 2013, 2012 and 2011 :
 for (ism=0;ism<kNbSMtot;ism++)
 //for (ism=0;ism<1;ism++)
    {if (tabChoiceCalos[SMdetType[ism]] == 0) continue;
     //CalcNewHVperSM(Form("/cebaf/faivre/recherche/calibPi0/biasSentByDavidMay2015/bias/%s/NewBias.txt",SMP2Name[ism]),Form("/cebaf/faivre/recherche/calibPi0/recalculateDCAL_HV_July2015/HVscanParametersGrenobleDCAL/HVscanParamsGrenobleDCAL/HvscanSM%s_%s_%s/parameters.txt",SMnumber[ism],SMcalibName[ism],SMP2Name[ism]),ism,fchNameBase,outFile,badGainFileOut,coeffOCDB,coeffCalib,forcedGainChange,forcedOCDBcoeff,flag,ps,c1); //This is WRONG.
     CalcNewHVperSM(Form("/cebaf/faivre/recherche/calibPi0/recalculateDCAL_HV_July2015/output_HVrecalculationDCAL2015July_mergedWithEMCalHV_lowerTemperature_createBias/%s/NewBias.txt",SMP2Name[ism]),Form("/cebaf/faivre/recherche/calibPi0/recalculateDCAL_HV_July2015/HVscanParametersGrenobleDCAL/HVscanParamsGrenobleDCAL/HvscanSM%s_%s_%s/parameters.txt",SMnumber[ism],SMcalibName[ism],SMP2Name[ism]),ism,fchNameBase,outFile,badGainFileOut,coeffOCDB,coeffCalib,forcedGainChange,forcedOCDBcoeff,flag,ps,c1);
     }

 fclose(badGainFileOut);
 
 // Close what has been opened :
 ps->Close();
 rootFileOut->Close();
 fclose(outFile);

 return;
 }








