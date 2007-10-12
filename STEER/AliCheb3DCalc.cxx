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

/* $Id$ */

// Author: ruben.shahoyan@cern.ch   09/09/2006
//
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliCheb3D produces the interpolation of the user 3D->NDimOut arbitrary     //
// function supplied in "void (*fcn)(float* inp,float* out)" format           //
// either in a separate macro file or as a function pointer.                  //
// Only coefficients needed to guarantee the requested precision are kept.    //
//                                                                            //
// The user-callable methods are:                                             //
// To create the interpolation use:                                           //
// AliCheb3D(const char* funName,  // name of the file with user function     //
//          or                                                                //
// AliCheb3D(void (*ptr)(float*,float*),// pointer on the  user function      //
//        Int_t     DimOut,     // dimensionality of the function's output    // 
//        Float_t  *bmin,       // lower 3D bounds of interpolation domain    // 
//        Float_t  *bmax,       // upper 3D bounds of interpolation domain    // 
//        Int_t    *npoints,    // number of points in each of 3 input        //
//                              // dimension, defining the interpolation grid //
//        Float_t   prec=1E-6); // requested max.absolute difference between  //
//                              // the interpolation and any point on grid    //
//                                                                            //
// To test obtained parameterization use the method                           //
// TH1* TestRMS(int idim,int npoints = 1000,TH1* histo=0);                    // 
// it will compare the user output of the user function and interpolation     //
// for idim-th output dimension and fill the difference in the supplied       //
// histogram. If no histogram is supplied, it will be created.                //
//                                                                            //
// To save the interpolation data:                                            //
// SaveData(const char* filename, Bool_t append )                             //
// write text file with data. If append is kTRUE and the output file already  //
// exists, data will be added in the end of the file.                         //
// Alternatively, SaveData(FILE* stream) will write the data to               //
// already existing stream.                                                   //
//                                                                            //
// To read back already stored interpolation use either the constructor       // 
// AliCheb3D(const char* inpFile);                                            //
// or the default constructor AliCheb3D() followed by                         //
// AliCheb3D::LoadData(const char* inpFile);                                  //
//                                                                            //
// To compute the interpolation use Eval(float* par,float *res) method, with  //
// par being 3D vector of arguments (inside the validity region) and res is   //
// the array of DimOut elements for the output.                               //
//                                                                            //
// If only one component (say, idim-th) of the output is needed, use faster   //
// Float_t Eval(Float_t *par,int idim) method.                                //
//                                                                            //
// void Print(option="") will print the name, the ranges of validity and      //
// the absolute precision of the parameterization. Option "l" will also print //
// the information about the number of coefficients for each output           //
// dimension.                                                                 //
//                                                                            //
// NOTE: during the evaluation no check is done for parameter vector being    //
// outside the interpolation region. If there is such a risk, use             //
// Bool_t IsInside(float *par) method. Chebyshev parameterization is not      //
// good for extrapolation!                                                    //
//                                                                            //
// For the properties of Chebyshev parameterization see:                      //
// H.Wind, CERN EP Internal Report, 81-12/Rev.                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliCheb3DCalc.h"

ClassImp(AliCheb3DCalc)


AliCheb3DCalc::AliCheb3DCalc():
    TNamed("", ""),
    fNCoefs(0),  
    fNRows(0),
    fNCols(0),
    fNElemBound2D(0),
    fNColsAtRow(0),
    fColAtRowBg(0),
    fCoefBound2D0(0),
    fCoefBound2D1(0),
    fCoefs(0),
    fTmpCf1(0),
    fTmpCf0(0)
{
    // Default constructor
    Init0();
}

AliCheb3DCalc::AliCheb3DCalc(FILE* stream):
    TNamed("", ""),
    fNCoefs(0),  
    fNRows(0),
    fNCols(0),
    fNElemBound2D(0),
    fNColsAtRow(0),
    fColAtRowBg(0),
    fCoefBound2D0(0),
    fCoefBound2D1(0),
    fCoefs(0),
    fTmpCf1(0),
    fTmpCf0(0)
{
    // Default constructor
    Init0();
    LoadData(stream);
}


AliCheb3DCalc::AliCheb3DCalc(const AliCheb3DCalc& src) :
    TNamed(src), 
    fNCoefs(src.fNCoefs), 
    fNRows(src.fNRows), 
    fNCols(src.fNCols),
    fNElemBound2D(src.fNElemBound2D), 
    fNColsAtRow(0), 
    fColAtRowBg(0), 
    fCoefBound2D0(0), 
    fCoefBound2D1(0), 
    fCoefs(0), 
    fTmpCf1(0), 
    fTmpCf0(0)
{
    // Copy constructor
    if (src.fNColsAtRow) {
	fNColsAtRow = new Int_t[fNRows]; 
	for (int i=fNRows;i--;) fNColsAtRow[i] = src.fNColsAtRow[i];
    }
    if (src.fColAtRowBg) {
	fColAtRowBg = new Int_t[fNRows]; 
	for (int i=fNRows;i--;) fColAtRowBg[i] = src.fColAtRowBg[i];
    }
    if (src.fCoefBound2D0) {
	fCoefBound2D0 = new Int_t[fNElemBound2D];
	for (int i=fNElemBound2D;i--;) fCoefBound2D0[i] = src.fCoefBound2D0[i];
    }
    if (src.fCoefBound2D1) {
	fCoefBound2D1 = new Int_t[fNElemBound2D];
	for (int i=fNElemBound2D;i--;) fCoefBound2D1[i] = src.fCoefBound2D1[i];
    }
    if (src.fCoefs) {
	fCoefs = new Float_t[fNCoefs];
	for (int i=fNCoefs;i--;) fCoefs[i] = src.fCoefs[i];
    }
    if (src.fTmpCf1) fTmpCf1 = new Float_t[fNCols];
    if (src.fTmpCf0) fTmpCf0 = new Float_t[fNRows];
}

AliCheb3DCalc& AliCheb3DCalc::operator=(const AliCheb3DCalc& rhs)
{
    // Assignment operator
    if (this != &rhs) {
	Clear();
	SetName(rhs.GetName());
	SetTitle(rhs.GetTitle());
	fNCoefs = rhs.fNCoefs;
	fNRows  = rhs.fNRows;
	fNCols  = rhs.fNCols;    
	if (rhs.fNColsAtRow) {
	    fNColsAtRow = new Int_t[fNRows]; 
	    for (int i=fNRows;i--;) fNColsAtRow[i] = rhs.fNColsAtRow[i];
	}
	if (rhs.fColAtRowBg) {
	    fColAtRowBg = new Int_t[fNRows]; 
	    for (int i=fNRows;i--;) fColAtRowBg[i] = rhs.fColAtRowBg[i];
	}
	if (rhs.fCoefBound2D0) {
	    fCoefBound2D0 = new Int_t[fNElemBound2D];
	    for (int i=fNElemBound2D;i--;) fCoefBound2D0[i] = rhs.fCoefBound2D0[i];
	}
	if (rhs.fCoefBound2D1) {
	    fCoefBound2D1 = new Int_t[fNElemBound2D];
	    for (int i=fNElemBound2D;i--;) fCoefBound2D1[i] = rhs.fCoefBound2D1[i];
	}
	if (rhs.fCoefs) {
	    fCoefs = new Float_t[fNCoefs];
	    for (int i=fNCoefs;i--;) fCoefs[i] = rhs.fCoefs[i];
	}
	if (rhs.fTmpCf1) fTmpCf1 = new Float_t[fNCols];
	if (rhs.fTmpCf0) fTmpCf0 = new Float_t[fNRows];    
    }
    return *this;
}

//__________________________________________________________________________________________
void AliCheb3DCalc::Clear(Option_t*)
{
  // delete all dynamycally allocated structures
  if (fTmpCf1)       { delete[] fTmpCf1;  fTmpCf1 = 0;}
  if (fTmpCf0)       { delete[] fTmpCf0;  fTmpCf0 = 0;}
  if (fCoefs)        { delete[] fCoefs;   fCoefs  = 0;}
  if (fCoefBound2D0) { delete[] fCoefBound2D0; fCoefBound2D0 = 0; }
  if (fCoefBound2D1) { delete[] fCoefBound2D1; fCoefBound2D1 = 0; }
  if (fNColsAtRow)   { delete[] fNColsAtRow;   fNColsAtRow = 0; }
  if (fColAtRowBg)   { delete[] fColAtRowBg;   fColAtRowBg = 0; }
  //
}

//__________________________________________________________________________________________
void AliCheb3DCalc::Init0()
{
  // reset everything to 0
  fNCoefs = fNRows = fNCols = fNElemBound2D = 0;
  fCoefs = 0;
  fCoefBound2D0 = fCoefBound2D1 = 0;
  fNColsAtRow = fColAtRowBg = 0;
  fTmpCf0 = fTmpCf1 = 0;
}

//__________________________________________________________________________________________
void AliCheb3DCalc::Print(Option_t* ) const
{
  // Print the Chebychev paramterization data
  printf("Chebyshev parameterization data %s for 3D->1 function.\n",GetName());
  int nmax3d = 0; 
  for (int i=fNElemBound2D;i--;) if (fCoefBound2D0[i]>nmax3d) nmax3d = fCoefBound2D0[i];
  printf("%d coefficients in %dx%dx%d matrix\n",fNCoefs,fNRows,fNCols,nmax3d);
  //
}

//__________________________________________________________________________________________
Float_t  AliCheb3DCalc::Eval(Float_t  *par) const
{
  // evaluate Chebyshev parameterization for 3D function.
  // VERY IMPORTANT: par must contain the function arguments ALREADY MAPPED to [-1:1] interval
  Float_t  &z = par[2];
  Float_t  &y = par[1];
  Float_t  &x = par[0];
  //
  int ncfRC;
  for (int id0=fNRows;id0--;) {
    int nCLoc = fNColsAtRow[id0];                   // number of significant coefs on this row
    int col0  = fColAtRowBg[id0];                   // beginning of local column in the 2D boundary matrix
    for (int id1=nCLoc;id1--;) {
      int id = id1+col0;
      fTmpCf1[id1] = (ncfRC=fCoefBound2D0[id]) ? ChebEval1D(z,fCoefs + fCoefBound2D1[id], ncfRC) : 0.0;
    }
    fTmpCf0[id0] = nCLoc>0 ? ChebEval1D(y,fTmpCf1,nCLoc):0.0;
  }
  return ChebEval1D(x,fTmpCf0,fNRows);
  //
}

//_______________________________________________
#ifdef _INC_CREATION_ALICHEB3D_
void AliCheb3DCalc::SaveData(const char* outfile,Bool_t append) const
{
  // writes coefficients data to output text file, optionallt appending on the end of existing file
  TString strf = outfile;
  gSystem->ExpandPathName(strf);
  FILE* stream = fopen(strf,append ? "a":"w");
  SaveData(stream);
  fclose(stream);
  //
}
#endif

//_______________________________________________
#ifdef _INC_CREATION_ALICHEB3D_
void AliCheb3DCalc::SaveData(FILE* stream) const
{
  // writes coefficients data to existing output stream
  // Note: fNCols, fNElemBound2D and fColAtRowBg is not stored, will be computed on fly during the loading of this file
  fprintf(stream,"#\nSTART %s\n",GetName());
  fprintf(stream,"# Number of rows\n%d\n",fNRows);
  //
  fprintf(stream,"# Number of columns per row\n");
  for (int i=0;i<fNRows;i++) fprintf(stream,"%d\n",fNColsAtRow[i]);
  //
  fprintf(stream,"# Number of Coefs in each significant block of third dimension\n");
  for (int i=0;i<fNElemBound2D;i++) fprintf(stream,"%d\n",fCoefBound2D0[i]);
  //
  fprintf(stream,"# Coefficients\n");
  for (int i=0;i<fNCoefs;i++) fprintf(stream,"%+.8e\n",fCoefs[i]);
  fprintf(stream,"END %s\n",GetName());
  //
}
#endif

//_______________________________________________
void AliCheb3DCalc::LoadData(FILE* stream)
{
  // Load coefs. from the stream 
  if (!stream) {Error("LoadData","No stream provided.\nStop"); exit(1);}
  TString buffs;
  Clear();
  ReadLine(buffs,stream);
  if (!buffs.BeginsWith("START")) {Error("LoadData","Expected: \"START <fit_name>\", found \"%s\"\nStop\n",buffs.Data());exit(1);}
  if (buffs.First(' ')>0) SetName(buffs.Data()+buffs.First(' ')+1);
  //
  ReadLine(buffs,stream); // NRows
  fNRows = buffs.Atoi(); 
  if (fNRows<1) {Error("LoadData","Expected: '<number_of_rows>', found \"%s\"\nStop\n",buffs.Data());exit(1);}
  //
  fNCols = 0;
  fNElemBound2D = 0;
  InitRows(fNRows);
  //
  for (int id0=0;id0<fNRows;id0++) {
    ReadLine(buffs,stream);               // n.cols at this row
    fNColsAtRow[id0] = buffs.Atoi();
    fColAtRowBg[id0] = fNElemBound2D;     // begining of this row in 2D boundary surface
    fNElemBound2D   += fNColsAtRow[id0];
    if (fNCols<fNColsAtRow[id0]) fNCols = fNColsAtRow[id0];
  }
  InitCols(fNCols);
  //  
  fNCoefs = 0; 
  InitElemBound2D(fNElemBound2D);
  //
  for (int i=0;i<fNElemBound2D;i++) {
    ReadLine(buffs,stream);               // n.coeffs at 3-d dimension for the given column/row
    fCoefBound2D0[i] = buffs.Atoi();
    fCoefBound2D1[i] = fNCoefs;
    fNCoefs += fCoefBound2D0[i];
  }
  //
  if (fNCoefs<=0) {Error("LoadData","Negtive (%d) number of Chebychef coeffs. is obtained.\nStop\n",fNCoefs);exit(1);}
  //
  InitCoefs(fNCoefs);
  for (int i=0;i<fNCoefs;i++) {
    ReadLine(buffs,stream);
    fCoefs[i] = buffs.Atof();
  }
  // check end_of_data record
  ReadLine(buffs,stream);
  if (!buffs.BeginsWith("END") || !buffs.Contains(GetName())) {
    Error("LoadData","Expected \"END %s\", found \"%s\".\nStop\n",GetName(),buffs.Data());
    exit(1);
  }
  //
}

//_______________________________________________
void AliCheb3DCalc::ReadLine(TString& str,FILE* stream) 
{
  // read single line from the stream, skipping empty and commented lines. EOF is not expected
  while (str.Gets(stream)) {
    str = str.Strip(TString::kBoth,' ');
    if (str.IsNull()||str.BeginsWith("#")) continue;
    return;
  }
  fprintf(stderr,"AliCheb3D::ReadLine: Failed to read from stream.\nStop");exit(1); // normally, should not reach here
}

//_______________________________________________
void AliCheb3DCalc::InitCols(int nc)
{
  // Set max.number of significant columns in the coefs matrix
  fNCols = nc;
  if (fTmpCf1) delete[] fTmpCf1;
  fTmpCf1 = new Float_t [fNCols];
}

//_______________________________________________
void AliCheb3DCalc::InitRows(int nr)
{
  // Set max.number of significant rows in the coefs matrix
  if (fNColsAtRow) delete[] fNColsAtRow;
  if (fColAtRowBg) delete[] fColAtRowBg;
  if (fTmpCf0)     delete[] fTmpCf0;
  fNRows = nr;
  fNColsAtRow = new Int_t[fNRows];
  fTmpCf0     = new Float_t [fNRows];
  fColAtRowBg = new Int_t[fNRows];
  for (int i=fNRows;i--;) fNColsAtRow[i] = fColAtRowBg[i] = 0;
}

//_______________________________________________
void AliCheb3DCalc::InitElemBound2D(int ne)
{
  // Set max number of significant coefs for given row/column of coefs 3D matrix
  if (fCoefBound2D0) delete[] fCoefBound2D0; 
  if (fCoefBound2D1) delete[] fCoefBound2D1; 
  fNElemBound2D = ne;
  fCoefBound2D0 = new Int_t[fNElemBound2D];
  fCoefBound2D1 = new Int_t[fNElemBound2D];
  for (int i=fNElemBound2D;i--;) fCoefBound2D0[i] = fCoefBound2D1[i] = 0;
}

//_______________________________________________
void AliCheb3DCalc::InitCoefs(int nc)
{
  // Set total number of significant coefs
  if (fCoefs) delete[] fCoefs; 
  fNCoefs = nc;
  fCoefs = new Float_t [fNCoefs];
  for (int i=fNCoefs;i--;) fCoefs[i] = 0.0;
}


