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

/*
$Log$
*/

#include "AliCalorimeter.h"
 
ClassImp(AliCalorimeter) // Class implementation to enable ROOT I/O
 
AliCalorimeter::AliCalorimeter()
{
// Default constructor, all parameters set to 0
 fNrows=0;
 fNcolumns=0;
 fNsignals=0;
 fNclusters=0;
 fMatrix=0;
 fClusters=0;
 fModules=0;
 fHmodules=0;
 fHclusters=0;
 fNvetos=0;
 fVetos=0;
}
///////////////////////////////////////////////////////////////////////////
AliCalorimeter::~AliCalorimeter()
{
// Destructor to delete memory allocated for matrix and cluster array
 if (fMatrix)
 {
  for (Int_t i=0; i<fNrows; i++)
  {
   delete [] fMatrix[i];
  }
   delete [] fMatrix;
 }
 fMatrix=0;
 if (fModules) delete fModules;
 fModules=0;
 if (fClusters)
 {
  fClusters->Delete();
  delete fClusters;
  fClusters=0;
 }
 if (fHmodules) delete fHmodules;
 fHmodules=0;
 if (fHclusters) delete fHclusters;
 fHclusters=0;
 if (fVetos)
 {
  fVetos->Delete();
  delete fVetos;
  fVetos=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliCalorimeter::AliCalorimeter(Int_t nrow, Int_t ncol)
{
// Create a calorimeter module matrix
 fNrows=nrow;
 fNcolumns=ncol;
 fNsignals=0;
 fNclusters=0;
 fClusters=0;
 fMatrix=new AliCalmodule*[nrow];
 for (Int_t i=0; i<nrow; i++)
 {
  fMatrix[i]=new AliCalmodule[ncol];
 }
 // Mark the edge modules
 for (Int_t j=0; j<ncol; j++)
 {
  fMatrix[0][j].SetEdgeOn();
  fMatrix[nrow-1][j].SetEdgeOn();
 }
 for (Int_t k=0; k<nrow; k++)
 {
  fMatrix[k][0].SetEdgeOn();
  fMatrix[k][ncol-1].SetEdgeOn();
 }
 
 fModules=new TObjArray();  // Default size, expanded automatically
 
 fHmodules=0;
 fHclusters=0;

 fNvetos=0;
 fVetos=0;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNrows()
{
// Provide the number of rows for the calorimeter module matrix
 return fNrows;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNcolumns()
{
// Provide the number of columns for the calorimeter module matrix
 return fNcolumns;
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetSignal(Int_t row, Int_t col, Float_t sig)
{
// Set the signal for a certain calorimeter module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  if (fMatrix[row-1][col-1].GetSignal() <= 0.) // only count new modules
  {
   fNsignals++;
   fModules->Add(&(fMatrix[row-1][col-1]));
  }
  fMatrix[row-1][col-1].SetSignal(row,col,sig);
 }
 else
 {
  cout << " *AliCalorimeter::SetSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::AddSignal(Int_t row, Int_t col, Float_t sig)
{
// Add the signal to a certain calorimeter module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  if (fMatrix[row-1][col-1].GetSignal() <= 0.) // only count new modules
  {
   fNsignals++;
   fModules->Add(&(fMatrix[row-1][col-1]));
  }
  fMatrix[row-1][col-1].AddSignal(row,col,sig);
 }
 else
 {
  cout << " *AliCalorimeter::AddSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Reset(Int_t row, Int_t col)
{
// Reset the signal for a certain calorimeter module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  fMatrix[row-1][col-1].SetSignal(row,col,0);
  fNsignals--;
  fModules->Remove(&(fMatrix[row-1][col-1]));
 }
 else
 {
  cout << " *AliCalorimeter::Reset* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Reset()
{
// Reset the signals for the complete calorimeter
// Normally this is done to prepare for the data of the next event
// Note : Module gains, edge and dead flags remain unchanged
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 fNsignals=0;
 for (Int_t i=0; i<fNrows; i++)
 {
  for (Int_t j=0; j<fNcolumns; j++)
  {
   fMatrix[i][j].SetSignal(i+1,j+1,0);
  }
 }
 if (fModules) fModules->Clear();

 fNclusters=0;
 if (fClusters)
 {
  fClusters->Delete();
  delete fClusters;
  fClusters=0;
 }

 fNvetos=0;
 if (fVetos)
 {
  fVetos->Delete();
  delete fVetos;
  fVetos=0;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalorimeter::GetSignal(Int_t row, Int_t col)
{
// Provide the signal of a certain calorimeter module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  return fMatrix[row-1][col-1].GetSignal();
 }
 else
 {
  cout << " *AliCalorimeter::GetSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetEdgeOn(Int_t row, Int_t col)
{
// Indicate a certain calorimeter module as 'edge module'
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  fMatrix[row-1][col-1].SetEdgeOn();
 }
 else
 {
  cout << " *AliCalorimeter::SetEdgeOn* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetEdgeOff(Int_t row, Int_t col)
{
// Indicate a certain calorimeter module as 'non-edge module'
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  fMatrix[row-1][col-1].SetEdgeOff();
 }
 else
 {
  cout << " *AliCalorimeter::SetEdgeOff* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetDead(Int_t row, Int_t col)
{
// Indicate a certain calorimeter module as 'dead module'
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  fMatrix[row-1][col-1].SetDead();
 
  // Increase the 'edge value' of surrounding modules
  Int_t rlow=row-1;
  Int_t rup=row+1;
  Int_t clow=col-1;
  Int_t cup=col+1;
 
  if (rlow < 1) rlow=row;
  if (rup > fNrows) rup=fNrows;
  if (clow < 1) clow=col;
  if (cup > fNcolumns) cup=fNcolumns;
 
  for (Int_t i=rlow; i<=rup; i++)
  {
   for (Int_t j=clow; j<=cup; j++)
   {
    fMatrix[i-1][j-1].EdgeUp();
   }
  }
 
  // Correct the 'edge value' for the dead module itself
  fMatrix[row-1][col-1].EdgeDown();
 }
 else
 {
  cout << " *AliCalorimeter::SetDead* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetAlive(Int_t row, Int_t col)
{
// Indicate a certain calorimeter module as 'active module'
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  fMatrix[row-1][col-1].SetAlive();
 
  // Decrease the 'edge value' of surrounding modules
  Int_t rlow=row-1;
  Int_t rup=row+1;
  Int_t clow=col-1;
  Int_t cup=col+1;
 
  if (rlow < 1) rlow=row;
  if (rup > fNrows) rup=fNrows;
  if (clow < 1) clow=col;
  if (cup > fNcolumns) cup=fNcolumns;
 
  for (Int_t i=rlow; i<=rup; i++)
  {
   for (Int_t j=clow; j<=cup; j++)
   {
    fMatrix[i-1][j-1].EdgeDown();
   }
  }
 
  // Correct the 'edge value' for the dead module itself
  fMatrix[row-1][col-1].EdgeUp();
 }
 else
 {
  cout << " *AliCalorimeter::SetAlive* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetGain(Int_t row, Int_t col, Float_t gain)
{
// Set the gain value for a certain calorimeter module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  fMatrix[row-1][col-1].SetGain(gain);
 }
 else
 {
  cout << " *AliCalorimeter::SetGain* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetPosition(Int_t row,Int_t col,Float_t* vec,TString f)
{
// Set the position in user coordinates for a certain calorimeter module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  fMatrix[row-1][col-1].SetPosition(vec,f);
 }
 else
 {
  cout << " *AliCalorimeter::SetPosition* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetEdgeValue(Int_t row, Int_t col)
{
// Provide the value of the edge flag of a certain module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  return fMatrix[row-1][col-1].GetEdgeValue();
 }
 else
 {
  cout << " *AliCalorimeter::GetEdgeValue* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetDeadValue(Int_t row, Int_t col)
{
// Provide the value of the dead flag of a certain module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  return fMatrix[row-1][col-1].GetDeadValue();
 }
 else
 {
  cout << " *AliCalorimeter::GetDeadValue* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalorimeter::GetGain(Int_t row, Int_t col)
{
// Provide the gain value of a certain module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  return fMatrix[row-1][col-1].GetGain();
 }
 else
 {
  cout << " *AliCalorimeter::GetGain* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::GetPosition(Int_t row,Int_t col,Float_t* vec,TString f)
{
// Return the position in user coordinates for a certain calorimeter module
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  fMatrix[row-1][col-1].GetPosition(vec,f);
 }
 else
 {
  cout << " *AliCalorimeter::GetPosition* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalorimeter::GetClusteredSignal(Int_t row, Int_t col)
{
// Provide the module signal after clustering
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
 {
  return fMatrix[row-1][col-1].GetClusteredSignal();
 }
 else
 {
  cout << " *AliCalorimeter::GetClusteredSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNsignals()
{
// Provide the number of modules that contain a signal
// Note : The number of modules marked 'dead' but which had a signal
//        are included.
 return fNsignals;
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Group(Int_t n)
{
// Group the individual modules into clusters
// Module signals of n rings around the central module will be grouped
 
 if (fNsignals > 0) // Directly return if no modules fired
 {
  if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
  if (fNclusters > 0) Ungroup(); // Restore unclustered situation if needed
 
  // Order the modules with decreasing signal
  AliCalmodule* ordered=new AliCalmodule[fNsignals]; // temp. array for ordered modules
  Sortm(ordered);
 
  // Clustering of modules. Start with the highest signal.
  if (fClusters)
  {
   fClusters->Delete();
   delete fClusters;
   fClusters=0;
  }
  fClusters=new TObjArray();
  fNclusters=0;
  Int_t row=0;
  Int_t col=0;
  AliCalcluster* c=0;
  for (Int_t i=0; i<fNsignals; i++)
  {
   row=ordered[i].GetRow();    // row number of cluster center
   col=ordered[i].GetColumn(); // column number of cluster center
   if (row>0 && row<=fNrows && col>0 && col<=fNcolumns)
   {
    // only use modules not yet used in a cluster
    if (fMatrix[row-1][col-1].GetClusteredSignal() > 0.)
    {
     c=new AliCalcluster;
     c->Start(fMatrix[row-1][col-1]); // module to start the cluster
     if (c->GetNmodules() > 0)        // cluster started successfully (no edge)
     {
      fClusters->Add(c);
      fNclusters++;       // update cluster counter
      AddRing(row,col,n); // add signals of n rings around the center
     }
     else
     {
      if (c) delete c;
      c=0;
     }
    }
   }
  }
 
  // Delete the temp. array
  delete [] ordered;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Sortm(AliCalmodule* ordered)
{
// Order the modules with decreasing signal
 
 Int_t nord=0;
 for (Int_t i=0; i<fNrows; i++) // loop over all modules of the matrix
 {
  for (Int_t ii=0; ii<fNcolumns; ii++)
  {
   if (fMatrix[i][ii].GetSignal() <= 0.) continue; // only take modules with a signal
 
   if (nord == 0) // store the first module with a signal at the first ordered position
   {
    nord++;
    ordered[nord-1]=fMatrix[i][ii];
    continue;
   }
 
   for (Int_t j=0; j<=nord; j++) // put module in the right ordered position
   {
    if (j == nord) // module has smallest signal seen so far
    {
     nord++;
     ordered[j]=fMatrix[i][ii]; // add module at the end
     break; // go for next matrix module
    }
 
    if (fMatrix[i][ii].GetSignal() < ordered[j].GetSignal()) continue;
 
    nord++;
    for (Int_t k=nord-1; k>j; k--) {ordered[k]=ordered[k-1];} // create empty position
    ordered[j]=fMatrix[i][ii]; // put module at empty position
    break; // go for next matrix module
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::AddRing(Int_t row, Int_t col, Int_t n)
{
// Add module signals of 1 ring around (row,col) to current cluster
// n denotes the maximum number of rings around cluster center
// Note : This function is used recursively
 
 if (n >= 1) // Check if any rings left for recursive calls
 {
  Float_t signal=GetSignal(row,col); // signal of (row,col) module
 
  Int_t lrow=row-1; if (lrow < 1) lrow=1;                 // row lowerbound for ring
  Int_t urow=row+1; if (urow > fNrows) urow=fNrows;       // row upperbound for ring
  Int_t lcol=col-1; if (lcol < 1) lcol=1;                 // col lowerbound for ring
  Int_t ucol=col+1; if (ucol > fNcolumns) ucol=fNcolumns; // row upperbound for ring
 
  for (Int_t i=lrow; i<=urow; i++)
  {
   for (Int_t j=lcol; j<=ucol; j++)
   {
    // add module(i,j) to cluster if the signal <= signal(row,col)
    if (fMatrix[i-1][j-1].GetSignal() <= signal)
    {
     ((AliCalcluster*)fClusters->At(fNclusters-1))->Add(fMatrix[i-1][j-1]);
    }
    AddRing(i,j,n-1); // Go for ring of modules around this (i,j) one
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNclusters()
{
// Provide the number of clusters
 return fNclusters;
}
///////////////////////////////////////////////////////////////////////////
AliCalcluster* AliCalorimeter::GetCluster(Int_t j)
{
// Provide cluster number j
// Note : j=1 denotes the first cluster
 if ((j >= 1) && (j <= fNclusters))
 {
  return (AliCalcluster*)fClusters->At(j-1);
 }
 else
 {
  cout << " *AliCalorimeter::GetCluster* cluster number : " << j
       << " out of range." << endl;
  cout << " -- Cluster number 1 (if any) returned " << endl;
  return (AliCalcluster*)fClusters->At(0);
 }
}
///////////////////////////////////////////////////////////////////////////
AliCalmodule* AliCalorimeter::GetModule(Int_t j)
{
// Provide 'fired' module number j
// Note : j=1 denotes the first 'fired' module
 if ((j >= 1) && (j <= fNsignals))
 {
  return (AliCalmodule*)fModules->At(j-1);
 }
 else
 {
  cout << " *AliCalorimeter::GetModule* module number : " << j
       << " out of range." << endl;
  cout << " -- Fired module number 1 (if any) returned " << endl;
  return (AliCalmodule*)fModules->At(0);
 }
}
///////////////////////////////////////////////////////////////////////////
TH2F* AliCalorimeter::DrawModules()
{
// Provide a lego plot of the module signals
 
 if (fHmodules)
 {
  fHmodules->Reset();
 }
 else
 {
  fHmodules=new TH2F("fHmodules","Module signals",
            fNcolumns,0.5,float(fNcolumns)+0.5,fNrows,0.5,float(fNrows)+0.5);
 
  fHmodules->SetDirectory(0); // Suppress global character of histo pointer
 }
 
 AliCalmodule* m;
 Float_t row,col,signal;
 for (Int_t i=0; i<fNsignals; i++)
 {
  m=(AliCalmodule*)fModules->At(i);
  if (m)
  {
   row=float(m->GetRow());
   col=float(m->GetColumn());
   signal=m->GetSignal();
   if (signal>0.) fHmodules->Fill(col,row,signal);
  }
 }
 
 fHmodules->Draw("lego");
 return fHmodules;
}
///////////////////////////////////////////////////////////////////////////
TH2F* AliCalorimeter::DrawClusters()
{
// Provide a lego plot of the cluster signals
 
 if (fHclusters)
 {
  fHclusters->Reset();
 }
 else
 {
  fHclusters=new TH2F("fHclusters","Cluster signals",
            fNcolumns,0.5,float(fNcolumns)+0.5,fNrows,0.5,float(fNrows)+0.5);
 
  fHclusters->SetDirectory(0); // Suppress global character of histo pointer
 }
 
 AliCalcluster* c;
 Float_t row,col,signal;
 for (Int_t i=0; i<fNclusters; i++)
 {
  c=(AliCalcluster*)fClusters->At(i);
  if (c)
  {
   row=float(c->GetRow());
   col=float(c->GetColumn());
   signal=c->GetSignal();
   if (signal>0.) fHclusters->Fill(col,row,signal);
  }
 }
 
 fHclusters->Draw("lego");
 return fHclusters;
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::LoadMatrix()
{
// Load the Calorimeter module matrix data back from the TObjArray
 
 // Create the module matrix space
 if (fMatrix)
 {
  for (Int_t k=0; k<fNrows; k++)
  {
   delete [] fMatrix[k];
  }
   delete [] fMatrix;
 }
 fMatrix=new AliCalmodule*[fNrows];
 for (Int_t i=0; i<fNrows; i++)
 {
  fMatrix[i]=new AliCalmodule[fNcolumns];
 }
 
 // Copy the module data back into the matrix
 AliCalmodule* m;
 Int_t row;
 Int_t col;
 for (Int_t j=0; j<fNsignals; j++)
 {
  m=(AliCalmodule*)fModules->At(j);
  row=m->GetRow();
  col=m->GetColumn();
  fMatrix[row-1][col-1]=*m;
  fModules->AddAt(&(fMatrix[row-1][col-1]),j); // Store new pointer
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Ungroup()
{
// Set the module signals back to the non-clustered situation
 
 if (!fMatrix) LoadMatrix(); // Restore matrix data in case of reading input
 
 Float_t signal=0;
 for (Int_t i=0; i<fNrows; i++)
 {
  for (Int_t j=0; j<fNcolumns; j++)
  {
   signal=fMatrix[i][j].GetSignal();
   fMatrix[i][j].SetClusteredSignal(signal);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::AddVetoSignal(Float_t* r,TString f,Float_t s)
{
// Associate an (extrapolated) AliSignal at location r as veto to the cal.
// Note : The default signal value (s) is 0
 if (!fVetos)
 {
  fNvetos=0;
  fVetos=new TObjArray();
 } 

 fVetos->Add(new AliSignal);
 fNvetos++;

 ((AliSignal*)fVetos->At(fNvetos-1))->SetPosition(r,f);
 ((AliSignal*)fVetos->At(fNvetos-1))->SetSignal(s);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNvetos()
{
// Provide the number of veto signals associated to the calorimeter
 return fNvetos;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliCalorimeter::GetVetoSignal(Int_t i)
{
// Provide access to the i-th veto signal of this calorimeter
// Note : The first hit corresponds to i=1

 if (i>0 && i<=fNvetos)
 {
  return (AliSignal*)fVetos->At(i-1);
 }
 else
 {
  cout << " *AliCalorimeter::GetVetoSignal* Signal number " << i
       << " out of range." << endl;
  cout << " --- First signal (if any) returned." << endl;
  return (AliSignal*)fVetos->At(0);
 }
}
///////////////////////////////////////////////////////////////////////////
