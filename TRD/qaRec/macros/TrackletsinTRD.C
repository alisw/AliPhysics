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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors: Yvonne C. Pachmayer <pachmay@physi.uni-heidelberg.de>        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#define TrackletsinTRD_cxx
#include "TrackletsinTRD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void TrackletsinTRD::Loop()
{
    //   In a ROOT session, you can do:
    //      Root > .L TrackletsinTRD.C
    //      Root > TrackletsinTRD t
    //      Root > t.GetEntry(12); // Fill t data members with entry number 12
    //      Root > t.Show();       // Show values of entry 12
    //      Root > t.Show(16);     // Read and show values of entry 16
    //      Root > t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch

    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);


    if (fChain == 0) return;

    // open run loader and load gAlice, kinematics and header
    AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
    if (!runLoader)
    {
        printf("Error: readKine:\nGetting run loader from file \"%s\" failed", filename);
        return;
    }

    runLoader->LoadHeader();
    runLoader->LoadRecPoints("TRD");
    TObjArray *module = new TObjArray();

    runLoader->CdGAFile();

    Int_t nEvents = runLoader->GetNumberOfEvents();
    const Int_t nEventsarray=nEvents;
    AliTRDcluster *cls = 0;

    Double_t x_clus[540][90];
    Double_t y_clus[540][90];

    for(Int_t ev = 0; ev < nEvents - 1; ev++)
    {
        TTree *tree = runLoader->GetTreeR("TRD", 0);
        tree->SetBranchAddress("TRDcluster", &module);

        int N = tree->GetEntries(); // number of chamber, max 540
        // Check number of clusters
        for(Int_t ind = 0; ind < N; ind++)
        {
            tree->GetEntry(ind);
            Int_t m = module->GetEntries();

            for (Int_t j = 0; j < m; j++)    // loop over clusters of one chamber
            {
                if(j>=90) break;
                if (cls != 0) delete cls;
                cls = (AliTRDcluster*)module->At(j);
                x_clus[ind][j]=cls->GetX();
                y_clus[ind][j]=cls->GetY();

            }
        }



        // loop over debug file and analysis

        Float_t xbmin = 0., xbmax = 380.;
        /* // histogram for debugging purpose
        Int_t nxbins = 3000, nybins = 280;
        Float_t ybmin = -70., ybmax = 70.;
        hxy = new TH2F("hxy",";x;y",nxbins,xbmin,xbmax,nybins,ybmin,ybmax);
        */

        Long64_t nentries = fChain->GetEntriesFast();

        Long64_t nbytes = 0, nb = 0;
        Float_t xpos[6];
        Float_t ypos[6];
        Int_t counter_test;
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);   nbytes += nb;
            xpos[layer]=xtrack;
            ypos[layer]=ytrack;


            if(layer==5)
            {
                Float_t x[6]= {xpos[0],xpos[1],xpos[2],xpos[3],xpos[4],xpos[5]};
                Float_t y[6]= {ypos[0],ypos[1],ypos[2],ypos[3],ypos[4],ypos[5]};

                gxy = new TGraph(6,x,y);

                TF1 *f1 = new TF1("f1","pol1",xbmin,xbmax);
                gxy->Fit("f1","RQ");
                // resulting function: y=ax+b --> x*f1->GetParameter(1)+f1->GetParameter(0)

                // graph and fit plotting only for checking
                /*
                 c1 = new TCanvas("c1","hough Transform",0,0,800,800);
                 c1->cd(1);
                 hxy->Draw();
                 gxy->SetMarkerStyle(3);
                 gxy->Draw("p");
                 */

                Float_t px;
                Float_t py;

                if(eventcounter==ev)
                {
                    for(Int_t b=0;b<540;b++)
                    {
                        if(b==dettracklet)
                        {
                            Int_t counter_distcalc;
                            Float_t distance_sum=0;

                            for(Int_t c=0;c<90;c++)
                            {
                                px=0;
                                py=0;

                                px=x_clus[b][c];
                                py=y_clus[b][c];

                                if(px!=0 &&  py!=0)
                                {
                                    // Function to calculate the distance of a point to the above fitted function
                                    Double_t distance = 0;
                                    if((TMath::Sqrt(f1->GetParameter(1)*f1->GetParameter(1)+1))!=0) distance=TMath::Abs(f1->GetParameter(1)*px-py+f1->GetParameter(0))/(TMath::Sqrt(f1->GetParameter(1)*f1->GetParameter(1)+1));
                                    // if(distance<10) cout << eventcounter << " " << b << " " << c << " " <<   distance << endl;
                                    if(distance>0.6 && distance<=2)
                                    {
                                        counter_distcalc++;
                                        distance_sum+=distance;
                                    }
                                }

                            }
                            if(counter_distcalc!=0)
                            {
                                // these are tracks which have additional clusters close by: 0.6 < distance <= 2
                            }
                            else
                            {
                                // these are good tracks no additional clusters close by
                               cout <<" good tracks " <<  endl;
                            }
                            counter_distcalc=0;
                        }
                    }
                }

                if(gxy) delete gxy;
                if(f1) delete f1;
            } // end of 5 layer loop


           


        }  // end loop over debug file and analysis

        runLoader->GetNextEvent();
    }

}
 