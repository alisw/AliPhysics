#include "AliTPCPoissonSolver.h"



/// \cond CLASSIMP
ClassImp(AliTPCPoissonSolver)
/// \endcond


// FIXME: the following values should come from the database
const Double_t AliTPCPoissonSolver::fgkTPCZ0    = 249.7;     ///< nominal gating grid position
const Double_t AliTPCPoissonSolver::fgkIFCRadius=  83.5;     ///< radius which renders the "18 rod manifold" best -> compare calc. of Jim Thomas
// compare gkIFCRadius=  83.05: Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
const Double_t AliTPCPoissonSolver::fgkOFCRadius= 254.5;     ///< Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
const Double_t AliTPCPoissonSolver::fgkZOffSet  =   0.2;     ///< Offset from CE: calculate all distortions closer to CE as if at this point
const Double_t AliTPCPoissonSolver::fgkCathodeV = -100000.0; ///< Cathode Voltage (volts)
const Double_t AliTPCPoissonSolver::fgkGG       =     -70.0; ///< Gating Grid voltage (volts)

const Double_t AliTPCPoissonSolver::fgkdvdE = 0.0024; ///< [cm/V] drift velocity dependency on the E field (from Magboltz for NeCO2N2 at standard environment)

const Double_t AliTPCPoissonSolver::fgkEM = -1.602176487e-19/9.10938215e-31; ///< charge/mass in [C/kg]
const Double_t AliTPCPoissonSolver::fgke0 = 8.854187817e-12;                 ///< vacuum permittivity [A·s/(V·m)]

Double_t AliTPCPoissonSolver::fgExactErr = 1e-4;
Double_t AliTPCPoissonSolver::fgConvErr = 1e-3;


AliTPCPoissonSolver::AliTPCPoissonSolver()
  : TNamed("poisson solver","solver"), fStrategy(kRelaxation)
{
	/**
	if (!fTPCParam) {
		fTPCParam=AliTPCcalibDB::Instance()->GetParameters();
		if (!fTPCParam)
		{ AliError("Could not get the TPC parameters"); return; }
	}	
	**/
	// default strategy
	fStrategy = kMultigrid;
	fExactPresent = kFALSE;
	fErrorConvergenceNorm2 = new TVectorD(fMgParameters.NMGCYCLE);
	fErrorConvergenceNormInf = new TVectorD(fMgParameters.NMGCYCLE);
	fError = new TVectorD(fMgParameters.NMGCYCLE);
	
	//fExactSolution == NULL;
}

AliTPCPoissonSolver::AliTPCPoissonSolver(const char *name,const char *title)
  : TNamed(name,title)
{
	fExactPresent = kFALSE;
	fErrorConvergenceNorm2 = new TVectorD(fMgParameters.NMGCYCLE);
	fErrorConvergenceNormInf = new TVectorD(fMgParameters.NMGCYCLE);
	fError = new TVectorD(fMgParameters.NMGCYCLE);
  /// constructor
}

AliTPCPoissonSolver::~AliTPCPoissonSolver() {
  /// virtual destructor
	delete fErrorConvergenceNorm2;
	delete fErrorConvergenceNormInf;
	delete fError;

}


/// Provides poisson solver in 2D
///
/// Based on the strategy (relaxation, multigrid or FFT)
///
/// \param arrayV TMatrixD& potential in matrix 
/// \param chargeDensity TMatrixD& charge density in matrix (side effect 
/// \param rows Int_t number of rows in the grid for discritization
/// \param columns Int_t number of columns in the grid for discritization
/// \param maxIterations Int_t maximum iteration for relaxation method
///
/// \return A fixed number that has nothing to do with what the function does
void AliTPCPoissonSolver::PoissonSolver2D(TMatrixD &arrayV, TMatrixD &chargeDensity, Int_t rows, Int_t columns, Int_t maxIterations) {

	switch (fStrategy) {
		case kMultigrid:
			PoissonMultigrid2D(arrayV,chargeDensity,rows,columns);
			break;
		case kFFT:
			break;
		default:
			PoissonRelaxation2D(arrayV,chargeDensity,rows,columns,maxIterations);
	}


}

/// Provides poisson solver in Cylindrical 3D (TPC geometry)
///
/// Strategy based on parameter settings (fStrategy and fMgParameters)provided
/// * Cascaded multi grid with S.O.R
/// * Geometric Multigrid
///		* Cycles: V, W, Full
///		* Relaxation: Jacobi, Weighted-Jacobi, Gauss-Seidel
///		* Grid transfer operators: Full, Half
/// * Spectral Methods (TODO)
///
/// \param arrayofArrayV TMatrixD** potential in 3D matrix 
/// \param arrayofChargeDensities TMatrixD** charge density in 3D matrix (side effect)
/// \param rows Int_t number of rows in the r direction of TPC
/// \param columns Int_t number of columns in z direction of TPC
/// \param phislices Int_t number of phislices in phi direction of T{C
/// \param maxIterations Int_t maximum iteration for relaxation method
/// \param symmetry Int_t symmetry or not
///
/// \pre Charge density distribution in **arrayofChargeDensities** is known and boundary values for **arrayofArrayV** are set
/// \post Numerical solution for potential distribution is calculated and stored in each rod at **arrayofArrayV**
void AliTPCPoissonSolver::PoissonSolver3D( 
	TMatrixD**arrayofArrayV, 
	TMatrixD**arrayofChargeDensities,
	Int_t rows, 
	Int_t columns,  
	Int_t phislices, 
	Int_t maxIterations,  
	Int_t symmetry
) 
{
	switch (fStrategy) {
		case kMultigrid:
			if (fMgParameters.isFull3D)
				PoissonMultigrid3D(arrayofArrayV,arrayofChargeDensities,rows,columns,phislices,symmetry);
			else
				PoissonMultigrid3D2D(arrayofArrayV,arrayofChargeDensities,rows,columns,phislices,symmetry);
			break;
		case kFFT:
			break;
		default:
			PoissonRelaxation3D(arrayofArrayV,arrayofChargeDensities,rows,columns,phislices,maxIterations,symmetry);
	}

}

/// Solve Poisson's Equation by Relaxation Technique in 2D (assuming cylindrical symmetry)
///
/// Solve Poissons equation in a cylindrical coordinate system. The arrayV matrix must be filled with the
/// boundary conditions on the first and last rows, and the first and last columns.  The remainder of the
/// array can be blank or contain a preliminary guess at the solution.  The Charge density matrix contains
/// the enclosed spacecharge density at each point. The charge density matrix can be full of zero's if
/// you wish to solve Laplaces equation however it should not contain random numbers or you will get
/// random numbers back as a solution.
/// Poisson's equation is solved by iteratively relaxing the matrix to the final solution.  In order to
/// speed up the convergence to the best solution, this algorithm does a binary expansion of the solution
/// space.  First it solves the problem on a very sparse grid by skipping rows and columns in the original
/// matrix.  Then it doubles the number of points and solves the problem again.  Then it doubles the
/// number of points and solves the problem again.  This happens several times until the maximum number
/// of points has been included in the array.
///
/// NOTE: In order for this algorithm to work, the number of rows and columns must be a power of 2 plus one.
/// So rows == 2**M + 1 and columns == 2**N + 1.  The number of rows and columns can be different.
///
/// Method for relaxation: S.O.R Weighted Jacobi
///
/// \param arrayV TMatrixD& potential in matrix 
/// \param chargeDensity TMatrixD& charge density in matrix (side effect 
/// \param rows Int_t number of rows in the grid for discritization
/// \param columns Int_t number of columns in the grid for discritization
/// \param maxIterations Int_t maximum iteration for relaxation method
///
/// \return A fixed number that has nothing to do with what the function does
///
///
/// Original code by Jim Thomas (STAR TPC Collaboration)
void AliTPCPoissonSolver::PoissonRelaxation2D(TMatrixD &arrayV, TMatrixD &chargeDensity, Int_t rows, Int_t columns, Int_t maxIterations) {
  const Float_t  gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (rows-1) ;  
  const Float_t  gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (columns-1) ;  
  const Float_t  ratio       =  gridSizeR*gridSizeR / (gridSizeZ*gridSizeZ) ;

  TMatrixD  arrayEr(rows,columns) ;
  TMatrixD  arrayEz(rows,columns) ;

  //Check that number of rows and columns is suitable for a binary expansion

  if ( !IsPowerOfTwo(rows-1) ) {
    AliError("PoissonRelaxation - Error in the number of rows. Must be 2**M - 1");
    return;
  }
  
  if ( !IsPowerOfTwo(columns-1) ) {
    AliError("PoissonRelaxation - Error in the number of columns. Must be 2**N - 1");
    return;
  }

  // Solve Poisson's equation in cylindrical coordinates by relaxation technique
  // Allow for different size grid spacing in R and Z directions
  // Use a binary expansion of the size of the matrix to speed up the solution of the problem

  Int_t iOne = (rows-1)/4 ;
  Int_t jOne = (columns-1)/4 ;
 
	// Coarse until loops
  Int_t loops = 1 + (int) ( 0.5 + TMath::Log2( (double) TMath::Max(iOne,jOne) ) ) ;
  
  
	// Loop while the matrix expands & the resolution increases.
  for ( Int_t count = 0 ; count < loops ; count++ ) {  

    Float_t tempGridSizeR = gridSizeR * iOne ;
    Float_t tempRatio     = ratio * iOne * iOne / ( jOne * jOne ) ;
    Float_t tempFourth    = 1.0 / (2.0 + 2.0*tempRatio) ;

    // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
    std::vector<float> coef1(rows) ;
    std::vector<float> coef2(rows) ;

    for ( Int_t i = iOne ; i < rows-1 ; i+=iOne ) {
      Float_t radius = AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
      coef1[i] = 1.0 + tempGridSizeR/(2*radius);
      coef2[i] = 1.0 - tempGridSizeR/(2*radius);
    }

    TMatrixD sumChargeDensity(rows,columns) ;


		// average charge at the coarse point
    for ( Int_t i = iOne ; i < rows-1 ; i += iOne ) {
      Float_t radius = AliTPCPoissonSolver::fgkIFCRadius + iOne*gridSizeR ;
      for ( Int_t j = jOne ; j < columns-1 ; j += jOne ) {	
				if ( iOne == 1 && jOne == 1 ) sumChargeDensity(i,j) = chargeDensity(i,j) ;
				else {
					// Add up all enclosed charge density contributions within 1/2 unit in all directions
					Float_t weight = 0.0 ;
					Float_t sum    = 0.0 ;
					sumChargeDensity(i,j) = 0.0 ;
					for ( Int_t ii = i-iOne/2 ; ii <= i+iOne/2 ; ii++ ) {
						for ( Int_t jj = j-jOne/2 ; jj <= j+jOne/2 ; jj++ ) {
							if ( ii == i-iOne/2 || ii == i+iOne/2 || jj == j-jOne/2 || jj == j+jOne/2 ) weight = 0.5 ;
							else
								weight = 1.0 ;					
							sumChargeDensity(i,j) += chargeDensity(ii,jj)*weight*radius ;
							sum += weight*radius ;
						}
					}
					sumChargeDensity(i,j) /= sum ;
				}
				sumChargeDensity(i,j) *= tempGridSizeR*tempGridSizeR; // just saving a step later on
			}
		}

		// Iterate on the current level
    for ( Int_t k = 1 ; k <= maxIterations; k++ ) {
      // Solve Poisson's Equation
      // Over-relaxation index, must be >= 1 but < 2.  Arrange for it to evolve from 2 => 1
      // as interations increase.
      Float_t overRelax   = 1.0 + TMath::Sqrt( TMath::Cos( (k*TMath::PiOver2())/maxIterations ) ) ;
      Float_t overRelaxM1 = overRelax - 1.0 ;
      Float_t overRelaxtempFourth, overRelaxcoef5 ;
      overRelaxtempFourth = overRelax * tempFourth ;
      overRelaxcoef5 = overRelaxM1 / overRelaxtempFourth ;

      for ( Int_t i = iOne ; i < rows-1 ; i += iOne ) {
				for ( Int_t j = jOne ; j < columns-1 ; j += jOne ) {
						// S.O.R
						// 
						arrayV(i,j) = (   coef2[i]       *   arrayV(i-iOne,j)
							+ tempRatio      * ( arrayV(i,j-jOne) + arrayV(i,j+jOne) )
							- overRelaxcoef5 *   arrayV(i,j)
							+ coef1[i]       *   arrayV(i+iOne,j)
							+ sumChargeDensity(i,j)
						) * overRelaxtempFourth;
				}
      }
	
			// if already at maxIterations 
			// TODO: stop when it converged
      if ( k == maxIterations ) {
	
				// After full solution is achieved, copy low resolution solution into higher res array
				// Interpolate solution
				for ( Int_t i = iOne ; i < rows-1 ; i += iOne ) {
					for ( Int_t j = jOne ; j < columns-1 ; j += jOne ) {
						if ( iOne > 1 ) {
							arrayV(i+iOne/2,j)                    =  ( arrayV(i+iOne,j) + arrayV(i,j)     ) / 2 ;
							if ( i == iOne )  arrayV(i-iOne/2,j) =  ( arrayV(0,j)       + arrayV(iOne,j) ) / 2 ;
						}
						if ( jOne > 1 ) {
							arrayV(i,j+jOne/2)                    =  ( arrayV(i,j+jOne) + arrayV(i,j) )     / 2 ;
							if ( j == jOne )  arrayV(i,j-jOne/2) =  ( arrayV(i,0)       + arrayV(i,jOne) ) / 2 ;
						}
						if ( iOne > 1 && jOne > 1 ) {
							arrayV(i+iOne/2,j+jOne/2) =  ( arrayV(i+iOne,j+jOne) + arrayV(i,j) ) / 2 ;
							if ( i == iOne ) arrayV(i-iOne/2,j-jOne/2) =   ( arrayV(0,j-jOne) + arrayV(iOne,j) ) / 2 ;
							if ( j == jOne ) arrayV(i-iOne/2,j-jOne/2) =   ( arrayV(i-iOne,0) + arrayV(i,jOne) ) / 2 ;
							// Note that this leaves a point at the upper left and lower right corners uninitialized.
							// -> Not a big deal.
						}
					}
				}
      }
    }

    iOne = iOne / 2 ; if ( iOne < 1 ) iOne = 1 ;
    jOne = jOne / 2 ; if ( jOne < 1 ) jOne = 1 ;
    sumChargeDensity.Clear();
  }
}

/// Solve Poisson's Equation by Multigrid Technique in 2D (assuming cylindrical symmetry)
///
/// NOTE: In order for this algorithmto work, the number of rows and columns must be a power of 2 plus one.
/// So rows == 2**M + 1 and columns == 2**N + 1.  The number of rows and columns can be different.
///
/// \param arrayV TMatrixD& potential in matrix 
/// \param chargeDensity TMatrixD& charge density in matrix (side effect 
/// \param rows Int_t number of rows in the grid for discritization
/// \param columns Int_t number of columns in the grid for discritization
/// \param maxIterations Int_t maximum iteration for relaxation method
///
/// \return A fixed number that has nothing to do with what the function does
void AliTPCPoissonSolver::PoissonMultigrid2D(TMatrixD &arrayV, TMatrixD &chargeDensity, Int_t rows, Int_t columns) {
	/// Geometry of TPC -- should be use AliTPCParams instead
  const Float_t  gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (rows-1) ;  
  const Float_t  gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (columns-1) ;  
  const Float_t  ratio       =  gridSizeR*gridSizeR / (gridSizeZ*gridSizeZ) ;
		
	
	Int_t ngridrows = 0; // number grid
  Int_t ngridcols = 0; // number grid
  Int_t nnrows;
  Int_t nncols;
	
  nnrows = rows;
  while (nnrows >>= 1) ngridrows++;
  
  nncols = columns;
  while (nncols >>= 1) ngridcols++;
    
  //Check that number of rows and columns is suitable for multi grid  
  if ( !IsPowerOfTwo(rows-1) ) {
    AliError("PoissonMultigrid - Error in the number of rows. Must be 2**M - 1");
    return;
  }
  if ( !IsPowerOfTwo(columns-1) ) {
    AliError("PoissonMultigrid - Error in the number of columns. Must be 2**N - 1");
    return;
  }
  
  
  Int_t loops = TMath::Max(ngridrows, ngridcols) ;      // Calculate the number of loops for the binary expansion
	
	AliInfo(Form("ngridrows=%d, ngridcols=%d, loops=%d, nmgcycles=%d",ngridrows,ngridcols,loops,fMgParameters.NMGCYCLE));
	
	
	Float_t h,h2,radius;
	Int_t iOne = 1; // in/dex
	Int_t jOne = 1; // index
	Int_t trows = rows,tcolumns = columns;
	Int_t count;
	Float_t tempRatio,tempFourth;
	
	// Vector for storing multi grid array
	
	
		
	std::vector<TMatrixD *> tvChargeFMG(loops);	
		
	std::vector<TMatrixD *> tvArrayV(loops);
	std::vector<TMatrixD *> tvCharge(loops);
	std::vector<TMatrixD *> tvResidu(loops);
	
	// Alocate memory for temprary grid
	for ( count = 1 ; count <= loops ; count++ ) {		
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;		
		// if one just address to arrayV
		tvResidu[count - 1] = new TMatrixD(trows,tcolumns);
		if (count == 1) {
			tvChargeFMG[count - 1] = &chargeDensity;
			tvArrayV[count - 1] = &arrayV;
			tvCharge[count - 1] = &chargeDensity;
		} else {
			tvArrayV[count - 1] = new TMatrixD(trows,tcolumns);
			tvCharge[count - 1] = new TMatrixD(trows,tcolumns);
			tvChargeFMG[count - 1] = new TMatrixD(trows,tcolumns);
			Restrict2D(*tvChargeFMG[count - 1],*tvChargeFMG[count - 2],trows,tcolumns);
			
		}
		iOne = 2*iOne;
		jOne = 2*jOne;				
	}
	
	/// full multi grid
	if (fMgParameters.cycleType == kFCycle) {
		
		AliInfo(Form("Do full cycle"));
		// FMG
		// 1) Relax on the coarsest grid
		iOne = iOne/2;
		jOne = jOne/2;
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		h = gridSizeR * count;
		h2 = h*h;
		tempRatio     = ratio * iOne * iOne / ( jOne * jOne ) ;
		tempFourth    = 1.0 / (2.0 + 2.0*tempRatio) ;
	
		std::vector<float> coef1(trows);
		std::vector<float> coef2(trows);
	
		for ( Int_t i = 1 ; i < trows-1 ; i++ ) {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h ;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);		
		}
	
		Relax2D(*tvArrayV[loops - 1],*tvChargeFMG[loops - 1],trows,tcolumns,h2,tempFourth,tempRatio,coef1,coef2);			
		
		
		// Do Vcycle from loops H to h
		for (  count = loops-2 ; count >=0; count-- ) {
		
			iOne = iOne/2;
			jOne = jOne/2;
		
			trows = iOne == 1 ? rows : rows/iOne + 1;
			tcolumns = jOne == 1 ? columns : columns/jOne + 1;
	
			//AliInfo(Form("Grid information: ione=%d, jone=%d, trows=%d, tcols=%d, count=%d\n", iOne,jOne,trows,tcolumns,count));
			Interp2D(*tvArrayV[count],*tvArrayV[count+1],trows, tcolumns);
			// Copy the relax charge to the tvCharge
			*tvCharge[count] = *tvChargeFMG[count]; //copy
			//tvCharge[count]->Print();
			// Do V cycle
		
			for (Int_t mgcycle=0;mgcycle < fMgParameters.NMGCYCLE  ;mgcycle++) {
				
				VCycle2D(rows,columns,count+1,loops,fMgParameters.NPRE, fMgParameters.NPOST, gridSizeR, ratio, tvArrayV, tvCharge, tvResidu);
				
				/// error computation
				//if (count == 0) {
				//	if (mgcycle < 100) {
				//	TMatrixD ExactMinStep = TMatrixD(rows,columns);		
				//	ExactMinStep = (*fExactSolution) - arrayV;
				//	(*fErrorProgress)(mgcycle) = ExactMinStep.E2Norm() / (rows*columns);
				//	}
				//}
				///
				
			}		
			
		}
	} else if (fMgParameters.cycleType == kVCycle) {
		// 2. VCycle
		AliInfo(Form("Do V cycle"));
		// VCycle(int gridfrom, int gridto, float gridsizer,float ratio, tvArrayV, tvCharge, tvResidu)
		Int_t gridfrom = 1;		
		
		Int_t gridto = loops;	
	
		// Do MGCycle Ntime
		for (Int_t mgcycle=0;mgcycle < fMgParameters.NMGCYCLE;mgcycle++) {
			VCycle2D(rows,columns,gridfrom,gridto,fMgParameters.NPRE, fMgParameters.NPOST,gridSizeR, ratio,tvArrayV,tvCharge,tvResidu);			
			
		}
	}	 else if (fMgParameters.cycleType == kWCycle) {
	
		// 3. W Cycle (TODO:)
		
		Int_t gridfrom = 1;
		
		//loops = loops >= 4 ? 4 : loops;
		
		Int_t gridto = loops;	
		//Int_t gamma = 1;
	
		// Do MGCycle Ntime
		for (Int_t mgcycle=0;mgcycle < fMgParameters.NMGCYCLE;mgcycle++) {
			WCycle2D(rows,columns,gridfrom,gridto,fMgParameters.gamma,fMgParameters.NPRE, fMgParameters.NPOST,gridSizeR, ratio,tvArrayV,tvCharge,tvResidu);			
			
		}
	}	
	
	
	
	// Dealocate memory
	for ( count = loops; count >= 1  ; count-- ) {		
		// if one just address to arrayV
		if (count > 1) {
			delete tvArrayV[count - 1];
			delete tvCharge[count - 1];
			delete tvChargeFMG[count - 1];
		}
		delete tvResidu[count - 1];
  }
}
			



/// 3D - Solve Poisson's Equation in 3D by Relaxation Technique
///
///    NOTE: In order for this algorith to work, the number of rows and columns must be a power of 2 plus one.
///    The number of rows and COLUMNS can be different.
///
///    ROWS       ==  2**M + 1
///    COLUMNS    ==  2**N + 1
///    PHISLICES  ==  Arbitrary but greater than 3
///
///    DeltaPhi in Radians
///
///    SYMMETRY = 0 if no phi symmetries, and no phi boundary conditions
/// = 1 if we have reflection symmetry at the boundaries (eg. sector symmetry or half sector symmetries).
///
/// \param arrayofArrayV TMatrixD** potential in 3D matrix 
/// \param arrayofChargeDensities TMatrixD** charge density in 3D matrix (side effect)
/// \param rows Int_t number of rows in the r direction of TPC
/// \param columns Int_t number of columns in z direction of TPC
/// \param phislices Int_t number of phislices in phi direction of T{C
/// \param maxIterations Int_t maximum iteration for relaxation method
/// \param symmetry Int_t symmetry or not
///
void AliTPCPoissonSolver::PoissonRelaxation3D( TMatrixD**arrayofArrayV, TMatrixD**arrayofChargeDensities,
		Int_t rows, Int_t columns,  Int_t phislices, Int_t maxIterations,  Int_t symmetry) {
			
  const Float_t  gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (rows-1) ;
  const Float_t  gridSizePhi =  TMath::TwoPi()/phislices; 
  const Float_t  gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (columns-1) ;
  const Float_t  ratioPhi    =  gridSizeR*gridSizeR / (gridSizePhi*gridSizePhi) ;
  const Float_t  ratioZ      =  gridSizeR*gridSizeR / (gridSizeZ*gridSizeZ) ;

	AliInfo(Form("in Poisson Solver 3D relaxation rows=%d, cols=%d, phislices=%d \n",rows,columns,phislices));
	
	
	
  //printf("(%f,%f,%f,%f,%f)\n",gridSizeR,gridSizePhi,gridSizeZ,ratioPhi,ratioZ);
  // Check that the number of rows and columns is suitable for a binary expansion
  if ( !IsPowerOfTwo((rows-1))    ) {
    AliError("Poisson3DRelaxation - Error in the number of rows. Must be 2**M - 1");
    return; }
  if ( !IsPowerOfTwo((columns-1)) ) {
    AliError("Poisson3DRelaxation - Error in the number of columns. Must be 2**N - 1");
    return; }
  if ( phislices <= 3   )  {
    AliError("Poisson3DRelaxation - Error in the number of phislices. Must be larger than 3");
    return; }
  if  ( phislices > 1000 ) {
    AliError("Poisson3D  phislices > 1000 is not allowed (nor wise) ");
    return; }

  // Solve Poisson's equation in cylindrical coordinates by relaxation technique
  // Allow for different size grid spacing in R and Z directions
  // Use a binary expansion of the matrix to speed up the solution of the problem

  Int_t loops, mplus, mminus, signplus, signminus;
  Int_t ione = (rows-1)/4 ;
  Int_t jone = (columns-1)/4 ;
  loops = TMath::Max(ione, jone) ;      // Calculate the number of loops for the binary expansion
  loops = 1 + (int) ( 0.5 + TMath::Log2((double)loops) ) ;  // Solve for N in 2**N

//	printf("loops = %d\n",loops);
  

  TMatrixD* arrayofSumChargeDensities[1000] ;    // Create temporary arrays to store low resolution charge arrays
 // printf("creating arrayofSumChargeDensity = %d\n",loops);
  

  // array of for storing convergence error (|V^(i) - V^{i-1}|)
  // temporary arrays to create the boundary conditions
  //TMatrixD *arrayofArrayVError[phislices];
  //TVectorD *vectorError = new TVectorD(phislices);
  //TVectorD *vectorErrorInf = new TVectorD(phislices);
  //TVectorD *vectorErrorNorm2Iter = new TVectorD(maxIterations);
  //TVectorD *vectorErrorInfIter = new TVectorD(maxIterations);
  
  //for ( Int_t k = 0 ; k < phislices ; k++ ) {
  //  arrayofArrayVError[k]     =   new TMatrixD(rows,columns) ;
 // }
  // end of array init for error
  
	//("allocating standard vector = %d\n",loops);
	std::vector<float> coef1(rows) ;  // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
	std::vector<float> coef2(rows) ;  // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
	std::vector<float> coef3(rows) ;  // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
	std::vector<float> coef4(rows) ;  // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
	std::vector<float> overRelaxcoef4(rows) ;  // Do this the standard C++ way to avoid gcc extensions
	std::vector<float> overRelaxcoef5(rows) ;  // Do this the standard C++ way to avoid gcc extensions
	
//	printf("finishing allocating standard vector = %d\n",loops);
  

  // do for each phi slices create the matrix 		
  for ( Int_t i = 0 ; i < phislices ; i++ ) 
  { arrayofSumChargeDensities[i] = new TMatrixD(rows,columns) ; } 
  
  
//	printf("finishing allocating standard vector each slice = %d\n",loops);
  


  AliSysInfo::AddStamp("3DInit", 10,0,0);

	///// Test of Convergence
	TMatrixD* prevArrayV[phislices] ;
	
	for (Int_t m=0;m<phislices;m++)
		prevArrayV[m] = new TMatrixD(rows,columns);
	/////


//	printf("finishing test prevArray = %d\n",loops);
  

	for ( Int_t count = 0 ; count < loops ; count++ ) {      // START the master loop and do the binary expansion
	AliSysInfo::AddStamp("3Diter", 20,count,0);
	
//	printf("count = %d\n", count);

	
	Float_t  tempgridSizeR   =  gridSizeR  * ione ;
	Float_t  tempratioPhi    =  ratioPhi * ione * ione ; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
	Float_t  tempratioZ      =  ratioZ   * ione * ione / ( jone * jone ) ;

//	printf("coef = %d\n", count);


	for ( Int_t i = ione ; i < rows-1 ; i+=ione )  {
		Float_t radius = AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
		coef1[i] = 1.0 + tempgridSizeR/(2*radius);
		coef2[i] = 1.0 - tempgridSizeR/(2*radius);
		coef3[i] = tempratioPhi/(radius*radius);
		coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
	}
	
//	printf("finish coef = %d\n", count);

	for ( Int_t m = 0 ; m < phislices ; m++ ) {
		TMatrixD &chargeDensity    = *arrayofChargeDensities[m] ;
		TMatrixD &sumChargeDensity = *arrayofSumChargeDensities[m] ;
		for ( Int_t i = ione ; i < rows-1 ; i += ione ) {
			Float_t radius = AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
			for ( Int_t j = jone ; j < columns-1 ; j += jone ) {
				if ( ione == 1 && jone == 1 ) sumChargeDensity(i,j) = chargeDensity(i,j) ;
				else {           // Add up all enclosed charge density contributions within 1/2 unit in all directions
					Float_t weight = 0.0 ;
					Float_t sum    = 0.0 ;
					sumChargeDensity(i,j) = 0.0 ;
					for ( Int_t ii = i-ione/2 ; ii <= i+ione/2 ; ii++ ) {
						for ( Int_t jj = j-jone/2 ; jj <= j+jone/2 ; jj++ ) {
							if ( ii == i-ione/2 || ii == i+ione/2 || jj == j-jone/2 || jj == j+jone/2 ) weight = 0.5 ;
							else
								weight = 1.0 ;
							sumChargeDensity(i,j) += chargeDensity(ii,jj)*weight*radius ;
							sum += weight*radius ;
						}
					}
					sumChargeDensity(i,j) /= sum ;
				}
				sumChargeDensity(i,j) *= tempgridSizeR*tempgridSizeR; // just saving a step later on
			}
		}
	}
	
	

	for ( Int_t k = 1 ; k <= maxIterations; k++ ) {
		
		
				/////	
				if (count==loops-1) {
					//// Test of Convergence
					for (Int_t m=0;m<phislices;m++)
						(*prevArrayV[m]) = (*arrayofArrayV[m]);			
					////
			
				}
	
		
		// over-relaxation index, >= 1 but < 2
		Float_t overRelax   = 1.0 + TMath::Sqrt( TMath::Cos( (k*TMath::PiOver2())/maxIterations ) ) ;
		Float_t overRelaxM1 = overRelax - 1.0 ;

		
		for ( Int_t i = ione ; i < rows-1 ; i+=ione ) {
			overRelaxcoef4[i] = overRelax * coef4[i] ;
			overRelaxcoef5[i] = overRelaxM1 / overRelaxcoef4[i] ;
		}

		for ( Int_t m = 0 ; m < phislices ; m++ ) {

			mplus  = m + 1;   signplus  = 1 ;
			mminus = m - 1 ;  signminus = 1 ;
			if (symmetry==1) {  // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
				if ( mplus  > phislices-1 ) mplus  = phislices - 2 ;
				if ( mminus < 0 )           mminus = 1 ;
			}
			else if (symmetry==-1) {   // Anti-symmetry in phi
				if ( mplus  > phislices-1 ) { mplus  = phislices - 2 ; signplus  = -1 ; }
				if ( mminus < 0 )           { mminus = 1 ;	         signminus = -1 ; }
			}
			else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
				if ( mplus  > phislices-1 ) mplus  = m + 1 - phislices ;
				if ( mminus < 0 )           mminus = m - 1 + phislices ;
			}
			
			TMatrixD& arrayV    =  *arrayofArrayV[m] ;
			TMatrixD& arrayVP   =  *arrayofArrayV[mplus] ;
			TMatrixD& arrayVM   =  *arrayofArrayV[mminus] ;
			TMatrixD& sumChargeDensity =  *arrayofSumChargeDensities[m] ;
			Double_t *arrayVfast = arrayV.GetMatrixArray();
			Double_t *arrayVPfast = arrayVP.GetMatrixArray();
			Double_t *arrayVMfast = arrayVM.GetMatrixArray();
			Double_t *sumChargeDensityFast=sumChargeDensity.GetMatrixArray();
			
			/// error declaration
		//	TMatrixD& arrayVError = *arrayofArrayVError[m];				
		//	Double_t *arrayVfastError = arrayVError.GetMatrixArray();
			///

			if (1)	{
			// slow implementation
				for ( Int_t i = ione ; i < rows-1 ; i+=ione )  {
					for ( Int_t j = jone ; j < columns-1 ; j+=jone ) {

						arrayV(i,j) = (   coef2[i]          *   arrayV(i-ione,j)
						+ tempratioZ        * ( arrayV(i,j-jone)  +  arrayV(i,j+jone) )
						- overRelaxcoef5[i] *   arrayV(i,j)
						+ coef1[i]          *   arrayV(i+ione,j)
						+ coef3[i]          * ( signplus*arrayVP(i,j)       +  signminus*arrayVM(i,j) )
						+ sumChargeDensity(i,j)
						) * overRelaxcoef4[i] ;
				// Note: over-relax the solution at each step.  This speeds up the convergance.
					}
				}
			}	else	{
				for ( Int_t i = ione ; i < rows-1 ; i+=ione )  {
					Double_t *arrayVfastI = &(arrayVfast[i*columns]);
					Double_t *arrayVPfastI = &(arrayVPfast[i*columns]);
					Double_t *arrayVMfastI = &(arrayVMfast[i*columns]);
					Double_t *sumChargeDensityFastI=&(sumChargeDensityFast[i*columns]);
					
					/// Error
	//				Double_t *arrayVfastIError = &(arrayVfastError[i*columns]);
					/// 
					
					for ( Int_t j = jone ; j < columns-1 ; j+=jone ) {
						Double_t /*resSlow*/resFast;
						
// 	      resSlow  = (   coef2[i]          *   arrayV(i-ione,j)
// 				+ tempratioZ        * ( arrayV(i,j-jone)  +  arrayV(i,j+jone) )
// 				- overRelaxcoef5[i] *   arrayV(i,j)
// 				+ coef1[i]          *   arrayV(i+ione,j)
// 				+ coef3[i]          * ( signplus*arrayVP(i,j)       +  signminus*arrayVM(i,j) )
// 				+ sumChargeDensity(i,j)
// 				) * overRelaxcoef4[i] ;
						resFast   = (   coef2[i]          *   arrayVfastI[j-columns*ione]
						+ tempratioZ        * ( arrayVfastI[j-jone]  +  arrayVfastI[j+jone] )
						- overRelaxcoef5[i] *   arrayVfastI[j]
						+ coef1[i]          * arrayVfastI[j+columns*ione]
						+ coef3[i]          * ( signplus* arrayVPfastI[j]      +  signminus*arrayVMfastI[j])
						+ sumChargeDensityFastI[j]
						) * overRelaxcoef4[i] ;
// 	      if (resSlow!=resFast){
// 		printf("problem\t%d\t%d\t%f\t%f\t%f\n",i,j,resFast,resSlow,resFast-resSlow);
// 	      }
//						arrayVfastIError[j] = resFast - arrayVfastI[j];
					
						arrayVfastI[j]=resFast;
					
					
	      // Note: over-relax the solution at each step.  This speeds up the convergance.
					} // end j
				} //end i
			} // end phi

			if ( k == maxIterations ) {   // After full solution is achieved, copy low resolution solution into higher res array
				for ( Int_t i = ione ; i < rows-1 ; i+=ione )  {
					for ( Int_t j = jone ; j < columns-1 ; j+=jone ) {

						if ( ione > 1 ) {
							arrayV(i+ione/2,j)                    =  ( arrayV(i+ione,j) + arrayV(i,j)     ) / 2 ;
							if ( i == ione )  arrayV(i-ione/2,j) =  ( arrayV(0,j)       + arrayV(ione,j) ) / 2 ;
						}
						if ( jone > 1 ) {
							arrayV(i,j+jone/2)                    =  ( arrayV(i,j+jone) + arrayV(i,j) )     / 2 ;
							if ( j == jone )  arrayV(i,j-jone/2) =  ( arrayV(i,0)       + arrayV(i,jone) ) / 2 ;
						}
						if ( ione > 1 && jone > 1 ) {
							arrayV(i+ione/2,j+jone/2) =  ( arrayV(i+ione,j+jone) + arrayV(i,j) ) / 2 ;
							if ( i == ione ) arrayV(i-ione/2,j-jone/2) =   ( arrayV(0,j-jone) + arrayV(ione,j) ) / 2 ;
							if ( j == jone ) arrayV(i-ione/2,j-jone/2) =   ( arrayV(i-ione,0) + arrayV(i,jone) ) / 2 ;
						// Note that this leaves a point at the upper left and lower right corners uninitialized. Not a big deal.
						}
					}
				}
			}
			
			// for each phi we calculate the norm of error
			//if (count == (loops-1)) {
		//		(*vectorError)(m) = arrayofArrayVError[m]->E2Norm();
		//		(*vectorErrorInf)(m) = arrayofArrayVError[m]->NormInf();
				// AliInfo(Form("count = %d, iter = %d\n",count,k));
				
		//	}
	//		if (count == (loops-1)) {
		//		AliInfo(Form("Error Norm2 = %f\n",vectorError->Norm2Sqr()));
			//	AliInfo(Form("Error Inf = %f\n",vectorErrorInf->NormInf()));
			//}
		}
		//if (count == (loops-1)) {
		//	(*vectorErrorNorm2Iter)(k-1) = vectorError->Norm2Sqr();
		//	(*vectorErrorInfIter)(k-1) = vectorErrorInf->NormInf();
		//}
		
	
					if (count == loops-1)
					{
						
						(*fErrorConvergenceNormInf)(k-1) = GetErrorConv(arrayofArrayV,prevArrayV,phislices);
						(*fError)(k-1) = GetErrorExact(arrayofArrayV,prevArrayV,phislices);
						
						// if error already achieved then stop mg iteration
						fIterations = k-1;
						if ((*fErrorConvergenceNormInf)(k-1)  <= fgConvErr) {
							AliInfo(Form("Exact Err: %f, Iteration : %d", (*fError)(k - 1), k - 1));							
							break;
						} 
						
						if (k == maxIterations) {
							AliInfo(Form("Exact Err: %f, Iteration : %d", (*fError)(k - 1), k - 1));							
						}
					
					}
	  }

    ione = ione / 2 ; if ( ione < 1 ) ione = 1 ;
    jone = jone / 2 ; if ( jone < 1 ) jone = 1 ;

  }
  
  //DrawError("Norm 2 - at last loop (grid size=513)",vectorErrorNorm2Iter);
  //DrawError("Norm Inf - at last loop (grid size=513)",vectorErrorInfIter);
  
  for ( Int_t k = 0 ; k < phislices ; k++ )
  {
      arrayofSumChargeDensities[k]->Delete() ;

	  //delete arrayofArrayVError[k];	  
  }
	
	for (Int_t m=0;m<phislices;m++)
		delete prevArrayV[m];
		//delete prevArrayV;
  //delete vectorError;
  //delete vectorErrorInf;
  //delete vectorErrorNorm2Iter;
  //delete vectorErrorInfIter;
			
	//fIterations = maxIterations;
}
			



			
/// 3D - Solve Poisson's Equation in 3D by Multigrid with constant phi slices
/// 
///    NOTE: In order for this algorith to work, the number of rows and columns must be a power of 2 plus one.
///    The number of rows and COLUMNS can be different.
///
///    ROWS       ==  2**M + 1
///    COLUMNS    ==  2**N + 1
///    PHISLICES  ==  Arbitrary but greater than 3
///
///		 Solving: \f$  \nabla^{2}V(r,\phi,z) = - f(r,\phi,z) \f$
///
/// Algorithm for MultiGrid Full Cycle (FMG)
/// - Relax on the coarsest grid
/// - Do from coarsest to finest
/// 	- Interpolate potential from coarse -> fine
///   - Do V-Cycle to the current coarse level to the coarsest
///   - Stop if converged
///
/// DeltaPhi in Radians
/// \param arrayofArrayV TMatrixD** potential in 3D matrix \f$ V(r,\phi,z) \f$
/// \param arrayofChargeDensities TMatrixD** charge density in 3D matrix (side effect) \f$ - f(r,\phi,z) \f$
/// \param rows Int_t number of rows in the r direction of TPC
/// \param columns Int_t number of columns in z direction of TPC
/// \param phislices Int_t number of phislices in phi direction of T{C
/// \param maxIterations Int_t maximum iteration for relaxation method (NOT USED)
/// \param symmetry Int_t symmetry (TODO for symmetry = 1)
//
///    SYMMETRY = 0 if no phi symmetries, and no phi boundary condition
///    = 1 if we have reflection symmetry at the boundaries (eg. sector symmetry or half sector symmetries).
///
void AliTPCPoissonSolver::PoissonMultigrid3D2D
( 
	TMatrixD**arrayofArrayV, 
	TMatrixD**arrayofChargeDensities,
	Int_t rows, 
	Int_t columns,  
	Int_t phislices,  
	Int_t symmetry
) 
{
	
	const Float_t  gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (rows-1); // h_{r}
  const Float_t  gridSizePhi =  TMath::TwoPi()/phislices;  // h_{phi}
  const Float_t  gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (columns-1) ; // h_{z}
  const Float_t  ratioPhi    =  gridSizeR*gridSizeR / (gridSizePhi*gridSizePhi) ;  // ratio_{phi} = gridsize_{r} / gridsize_{phi}
  const Float_t  ratioZ      =  gridSizeR*gridSizeR / (gridSizeZ*gridSizeZ) ; // ratio_{Z} = gridsize_{r} / gridsize_{z}
	
	// error tolerate
	//const Float_t  ERR = 1e-8;
	Double_t  convErr;

  AliInfo(Form("in Poisson Solver 3D multigrid semi coarsening rows=%d, cols=%d, phislices=%d \n",rows,columns,phislices));
  
  // Check that the number of rows and columns is suitable for a binary expansion
  if ( !IsPowerOfTwo((rows-1))    ) {
    AliError("Poisson3DMultigrid - Error in the number of rows. Must be 2**M + 1");
    return; }
  if ( !IsPowerOfTwo((columns-1)) ) {
    AliError("Poisson3DMultigrid - Error in the number of columns. Must be 2**N - 1");
    return; }
  if ( phislices <= 3   )  {
    AliError("Poisson3DMultigrid - Error in the number of phislices. Must be larger than 3");
    return; }
  if  ( phislices > 1000 ) {
    AliError("Poisson3D  phislices > 1000 is not allowed (nor wise) ");
    return; }

  // Solve Poisson's equation in cylindrical coordinates by multigrid technique
  // Allow for different size grid spacing in R and Z directions 
	
	Int_t ngridrows = 0; // number grid
  Int_t ngridcols = 0; // number grid
  Int_t nnrows;
  Int_t nncols;
		
  nnrows = rows;
  while (nnrows >>= 1) ngridrows++;  
  nncols = columns;
  while (nncols >>= 1) ngridcols++;
  
  
	// memory allocation for multigrids	
  //AliInfo(Form("nnrows=%d,nncols=%d,ngridrows=%d,ngridcols=%d\n",nnrows,nncols,ngridrows,ngridcols));
	
  
  Int_t loops = TMath::Max(ngridrows, ngridcols) ;      // Calculate the number of loops for the binary expansion
	
	//AliInfo(Form("Allocating memory=%d \n",loops));
	
	//AliInfo(Form("MAXLOOP=%d \n",fMgParameters.MAXLOO//P));
	
	
  loops = (loops > fMgParameters.MAXLOOP) ?  fMgParameters.MAXLOOP : loops;
	
	//AliInfo(Form("Allocating memory=%d \n",loops));
	
	// memory allocation for multigrids	
  //AliInfo(Form("Allocating memory=%d \n",loops));
	
	
	// for tracking the power of two
	Int_t count;
	
	// Vector for storing multi grid array
	Int_t iOne = 1; // index i in gridsize r (original)
	Int_t jOne = 1; // index j in gridsize z (original)
	Int_t trows = rows,tcolumns = columns;
	
	
	// vectors for storing multigrid charges, potential, residues, error
	std::vector<TMatrixD **> tvChargeFMG(loops);			 	// charge is restricted in full multigrid
	std::vector<TMatrixD **> tvArrayV(loops);						// potential <--> error
	std::vector<TMatrixD **> tvCharge(loops);						// charge <--> residu
	std::vector<TMatrixD **> tvResidu(loops);						// residu calculation
	std::vector<TMatrixD **> tvPrevArrayV(loops);				// error calculation
	
	// memory allocation for multigrids	
  //AliInfo(Form("Allocating memory=%d \n",loops));
	
	
	for ( count = 1 ; count <= loops ; count++ ) {		
		// trows,tcolumns in new grid
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;				
		
		// allocate memory for residu
		// memory allocation for multigrids	
		//AliInfo(Form("Allocating memory for residu=%d \n",phislices));
		//AliInfo(Form("Count=%d \n",count - 1));
	
		tvResidu[count - 1] = new TMatrixD*[phislices];
		tvPrevArrayV[count - 1] = new TMatrixD*[phislices];		
		for ( Int_t k = 0 ; k < phislices ; k++ ) {
			tvResidu[count - 1][k] = new TMatrixD(trows,tcolumns);			
			tvPrevArrayV[count - 1][k] = new TMatrixD(trows,tcolumns);					
		}
		
		// memory for the finest grid is from parameters
		if (count == 1) {
			tvChargeFMG[count - 1] = arrayofChargeDensities;
			tvArrayV[count - 1] = arrayofArrayV;
			tvCharge[count - 1] = arrayofChargeDensities;
		} else {
			// allocate for coarser grid
			tvChargeFMG[count - 1] = new TMatrixD*[phislices];
			tvArrayV[count - 1] = new TMatrixD*[phislices];
			tvCharge[count - 1] = new TMatrixD*[phislices];
			for ( Int_t k = 0 ; k < phislices ; k++ ) {
				tvArrayV[count - 1][k] = new TMatrixD(trows,tcolumns) ;
				tvCharge[count - 1][k] = new TMatrixD(trows,tcolumns) ;
				tvChargeFMG[count - 1][k] = new TMatrixD(trows,tcolumns) ;				
			}
			
			// restrict charge to coarser grids 
			Restrict3D(tvChargeFMG[count - 1],tvChargeFMG[count - 2], trows, tcolumns,  phislices, phislices);
			
			// pass boundary information to coarse level
			RestrictBoundary3D(tvArrayV[count - 1],tvArrayV[count - 2], trows, tcolumns,  phislices, phislices);

		}
		
		iOne = 2*iOne; // doubling
		jOne = 2*jOne; // doubling				
	}
	
	
	
//	AliInfo(Form("Finish allocation memeor\n"));
	
	Float_t h,h2,radius;
	Float_t tempratioPhi,tempratioZ;
	
	
//	AliInfo(Form("vectprs memeor: %d\n",rows));
	
	
	std::vector<float> coef1(rows);  // coef1(rows) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
	std::vector<float> coef2(rows);  // coef2(rows) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
	std::vector<float> coef3(rows);  // coef3(rows) for storing (1/r_{i}^2) from central differences in phi direction
	std::vector<float> coef4(rows);  // coef4(rows) for storing  1/2
	std::vector<float> icoef4(rows);  // inverse of coef4(rows)
	
				
	// Case full multi grid (FMG)
	if (fMgParameters.cycleType == kFCycle) {
				
		// 1) Relax on the coarsest grid
		iOne = iOne/2;
		jOne = jOne/2;
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
			
		h = gridSizeR *  iOne;
		h2 = h*h;
		
		tempratioPhi    =  ratioPhi * iOne * iOne ; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
		tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;

		for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
		}
	
		
		// relax on the coarsest level
		Relax3D(tvArrayV[loops - 1],tvChargeFMG[loops - 1],trows,tcolumns,phislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
		
		// 2) Do multigrid v-cycle from coarsest to finest
		for (  count = loops-2 ; count >=0; count-- ) {
			
			// move to finer grid
			iOne = iOne/2;
			jOne = jOne/2;
		
			trows = iOne == 1 ? rows : rows/iOne + 1;
			tcolumns = jOne == 1 ? columns : columns/jOne + 1;
	
			// 2) a) Interpolate potential for h -> 2h (coarse -> fine)
			Interp3D(tvArrayV[count],tvArrayV[count+1],trows, tcolumns,phislices,phislices);
	
	
			// 2) c) Copy the restricted charge to charge for calculation
			for (Int_t m=0;m<phislices;m++) {
				*tvCharge[count][m] = *tvChargeFMG[count][m]; //copy
			}
	
			// 2) c) Do V cycle fMgParameters.NMGCYCLE times at most
			for (Int_t mgcycle=0;mgcycle < fMgParameters.NMGCYCLE ;mgcycle++) {
				
				// Copy the potential to temp array for convergence calculation
				for (Int_t m=0;m<phislices;m++) {
					*tvPrevArrayV[count][m] = *tvArrayV[count][m]; //copy
				}	
					
			
				// 2) c) i) Call V cycle from grid count+1 (current fine level) to loops (coarsest)
				VCycle3D2D(rows,columns,phislices, symmetry,count+1,loops,fMgParameters.NPRE, fMgParameters.NPOST,gridSizeR, ratioZ,ratioPhi,tvArrayV, tvCharge, tvResidu,coef1,coef2,coef3,coef4,icoef4);
				
				
				convErr = GetErrorConv(tvArrayV[count],tvPrevArrayV[count],phislices);
				//AliInfo(Form("Conv Err: %1.8f, count=%d, MG Iteration : %d", convErr, count, mgcycle));
					
				//// error counting /////
				if (count == 0){					
					(*fErrorConvergenceNormInf)(mgcycle) = convErr;
					(*fError)(mgcycle) = GetErrorExact(arrayofArrayV,tvPrevArrayV[count],phislices);
				}
				/// if already converge just break move to finer grid
				if (convErr <= fgConvErr) {
					//AliInfo(Form("Conv Err: %f, MG Iteration : %d", convErr, mgcycle));
					fIterations = mgcycle+1;
					break;
				}	
			}
		}				
	}  // Case V multi grid (VMG)
	else if (fMgParameters.cycleType == kVCycle) {
		Int_t gridfrom = 1;
		Int_t gridto = loops;	
		
		// do v cycle fMgParameters.NMGCYCLE from the coarsest to finest
		for (Int_t mgcycle=0;mgcycle < fMgParameters.NMGCYCLE;mgcycle++) {			
			

			// copy to store previous potential
			for (Int_t m=0;m<phislices;m++) {
				*tvPrevArrayV[0][m] = *tvArrayV[0][m]; //copy
			}	
			
			// Do V Cycle for constant phislices
			VCycle3D2D(rows,columns,phislices, symmetry,gridfrom,gridto,fMgParameters.NPRE, fMgParameters.NPOST,gridSizeR, ratioZ,ratioPhi,tvArrayV, tvCharge, tvResidu,coef1,coef2,coef3,coef4,icoef4);
			
			// convergence error
			convErr = GetErrorConv(tvArrayV[0],tvPrevArrayV[0],phislices);
			//AliInfo(Form("Conv Err: %f, MG Iteration : %d", convErr, mgcycle));
				
			(*fErrorConvergenceNormInf)(mgcycle) = convErr;
			(*fError)(mgcycle) = GetErrorExact(arrayofArrayV,tvPrevArrayV[0],phislices);
					
			// if error already achieved then stop mg iteration
			if (convErr <= fgConvErr) {
				//AliInfo(Form("Exact Err: %f, MG Iteration : %d", (*fError)(mgcycle), mgcycle));
				//AliInfo(Form("Conv Err: %f, MG Iteration : %d", convErr, mgcycle));
				
				fIterations = mgcycle+1;
				break;
			}
		}
	}		
	
	
	// Dealocate memory
	for ( count = 1 ; count <= loops ; count++ ) {		
		delete[] tvResidu[count - 1];
		delete[] tvPrevArrayV[count - 1];
		
		if (count > 1) {
			delete[] tvChargeFMG[count - 1];
			delete[] tvArrayV[count - 1];
			delete[] tvCharge[count - 1];
		}	
	}
	
}			
		
		




		

/// 3D - Solve Poisson's Equation in 3D in all direction by Multigrid
///
///    NOTE: In order for this algorith to work, the number of rows and columns must be a power of 2 plus one.
///    The number of rows and COLUMNS can be different.
///
///    ROWS       ==  2**M + 1
///    COLUMNS    ==  2**N + 1
///    PHISLICES  ==  Arbitrary but greater than 3
///
///		 Solving: \f$  \nabla^{2}V(r,\phi,z) = - f(r,\phi,z) \f$
///
///  Algorithm for MultiGrid Full Cycle (FMG)
/// - Relax on the coarsest grid
/// - Do from coarsest to finest
/// 	- Interpolate potential from coarse -> fine
///   - Do V-Cycle to the current coarse level to the coarsest
///   - Stop if converged
///
///    DeltaPhi in Radians
/// \param arrayofArrayV TMatrixD** potential in 3D matrix 
/// \param arrayofChargeDensities TMatrixD** charge density in 3D matrix (side effect)
/// \param rows Int_t number of rows in the r direction of TPC
/// \param columns Int_t number of columns in z direction of TPC
/// \param phislices Int_t number of phislices in phi direction of T{C
/// \param maxIterations Int_t maximum iteration for relaxation method
/// \param symmetry Int_t symmetry or not
//
///    SYMMETRY = 0 if no phi symmetries, and no phi boundary condition
/// = 1 if we have reflection symmetry at the boundaries (eg. sector symmetry or half sector symmetries).
///
void AliTPCPoissonSolver::PoissonMultigrid3D( TMatrixD**arrayofArrayV, TMatrixD**arrayofChargeDensities,
		Int_t rows, Int_t columns,  Int_t phislices,  Int_t symmetry) {
	
	const Float_t  gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (rows-1); // h_{r}  
  const Float_t  gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (columns-1) ; // h_{z}
  const Float_t  ratioZ      =  gridSizeR*gridSizeR / (gridSizeZ*gridSizeZ) ; // ratio_{Z} = gridsize_{r} / gridsize_{z}
	
	Float_t  gridSizePhi =  TMath::TwoPi()/phislices;  // h_{phi}
	//Float_t  ratioPhi    =  gridSizeR*gridSizeR / (gridSizePhi*gridSizePhi) ;  // ratio_{phi} = gridsize_{r} / gridsize_{phi}
	
	Float_t h,h2,radius;
	Float_t tempratioPhi,tempratioZ;
	
	
	Float_t convErr; // Convergence error

  AliInfo(Form("in Poisson Solver 3D multigrid full coarsening  rows=%d, cols=%d, phislices=%d \n",rows,columns,phislices));
  
  // Check that the number of rows and columns is suitable for a binary expansion
  if ( !IsPowerOfTwo((rows-1))    ) {
    AliError("Poisson3DMultigrid - Error in the number of rows. Must be 2**M + 1");
    return; }
  if ( !IsPowerOfTwo((columns-1)) ) {
    AliError("Poisson3DMultigrid - Error in the number of columns. Must be 2**N - 1");
    return; }
  if ( phislices <= 3   )  {
    AliError("Poisson3DMultigrid - Error in the number of phislices. Must be larger than 3");
    return; }
  if  ( phislices > 1000 ) {
    AliError("Poisson3D  phislices > 1000 is not allowed (nor wise) ");
    return; }

  // Solve Poisson's equation in cylindrical coordinates by multigrid technique
  // Allow for different size grid spacing in R and Z directions 
	
	Int_t ngridrows = 0; // number grid
  Int_t ngridcols = 0; // number grid
	Int_t ngridphis = 0;
	
  Int_t nnrows;
  Int_t nncols;
	Int_t nnphis;
	
  nnrows = rows;
  while (nnrows >>= 1) ngridrows++;
  
  nncols = columns;
  while (nncols >>= 1) ngridcols++;
  
	nnphis = phislices;
	
	while (nnphis % 2 == 0) { 
		ngridphis++;
		nnphis /= 2;
	}
	
	AliInfo(Form("ngridrows=%d, ngridcols=%d, ngridphis=%d",ngridrows,ngridcols,ngridphis));
  Int_t loops = TMath::Max(ngridrows, ngridcols) ;      // Calculate the number of loops for the binary expansion
	loops = TMath::Max(loops,ngridphis);
	
	
	
	
	// Vector for storing multi grid array
	Int_t iOne = 1; // index i in gridsize r (original)
	Int_t jOne = 1; // index j in gridsize z (original)
	Int_t kOne = 1; // index k in gridsize phi
	Int_t trows = rows,tcolumns = columns, tphislices = phislices,otphislices;



	// 1)	Memory allocation for multigrids
	std::vector<TMatrixD **> tvChargeFMG(loops);			 	// charge is restricted in full multigrid
	std::vector<TMatrixD **> tvArrayV(loops);						// potential <--> error
	std::vector<TMatrixD **> tvCharge(loops);						// charge <--> residu
	std::vector<TMatrixD **> tvResidu(loops);						// residu calculation
	std::vector<TMatrixD **> tvPrevArrayV(loops);				// error calculation
	
	// these vectors for storing the coefs in smoother
	std::vector<float> coef1(rows);  // coef1(rows) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
	std::vector<float> coef2(rows);  // coef2(rows) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
	std::vector<float> coef3(rows);  // coef3(rows) for storing (1/r_{i}^2) from central differences in phi direction
	std::vector<float> coef4(rows);  // coef4(rows) for storing  1/2
	std::vector<float> icoef4(rows);  // inverse of coef4(rows)
	
	for (Int_t count = 1 ; count <= loops ; count++ ) {		
		
		// trows,tcolumns in new grid
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;				
		tphislices = kOne == 1 ? phislices : phislices/kOne;
		tphislices = tphislices < nnphis ? nnphis : tphislices;
		
		// allocate memory for residu
		tvResidu[count - 1] = new TMatrixD*[tphislices];		
		tvPrevArrayV[count - 1] = new TMatrixD*[tphislices];		
		for ( Int_t k = 0 ; k < tphislices ; k++ ) {
			tvResidu[count - 1][k] = new TMatrixD(trows,tcolumns);			
			tvPrevArrayV[count - 1][k] = new TMatrixD(trows,tcolumns);			
		}
		
		// memory for the finest grid is from parameters
		if (count == 1) {
			tvChargeFMG[count - 1] = arrayofChargeDensities;
			tvArrayV[count - 1] = arrayofArrayV;
			tvCharge[count - 1] = arrayofChargeDensities;
		} else {
			// allocate for coarser grid
			tvChargeFMG[count - 1] = new TMatrixD*[tphislices];
			tvArrayV[count - 1] = new TMatrixD*[tphislices];
			tvCharge[count - 1] = new TMatrixD*[tphislices];
			for ( Int_t k = 0 ; k < tphislices ; k++ ) {
				tvArrayV[count - 1][k] = new TMatrixD(trows,tcolumns) ;
				tvCharge[count - 1][k] = new TMatrixD(trows,tcolumns) ;
				tvChargeFMG[count - 1][k] = new TMatrixD(trows,tcolumns) ;				
			}
		}
		iOne = 2*iOne; // doubling
		jOne = 2*jOne; // doubling				
		kOne = 2*kOne;
	}

	
	// Case full multi grid (FMG)
	
	if (fMgParameters.cycleType == kFCycle)  {
		// Restrict the charge to coarser grid
		iOne = 2;
		jOne = 2;
		kOne = 2;
		otphislices = phislices;
		
		// 1) Restrict Charge and Boundary to coarser grid
		for (Int_t count = 2 ; count <= loops ; count++ ) {		
			// trows,tcolumns in new grid
			trows = iOne == 1 ? rows : rows/iOne + 1;
			tcolumns = jOne == 1 ? columns : columns/jOne + 1;				
			tphislices = kOne == 1 ? phislices : phislices/kOne;
			tphislices = tphislices < nnphis ? nnphis : tphislices;
			
			
			AliInfo(Form("Restrict3D, trows=%d, tcolumns=%d, newphislices=%d, oldphislices=%d\n",trows,tcolumns,tphislices,otphislices));
			Restrict3D(tvChargeFMG[count - 1],tvChargeFMG[count - 2], trows, tcolumns, tphislices, otphislices);			
			// copy boundary values of V
			RestrictBoundary3D(tvArrayV[count - 1],tvArrayV[count - 2], trows, tcolumns, tphislices, otphislices);			
			otphislices = tphislices;
			
			
			
			iOne = 2*iOne; // doubling
			jOne = 2*jOne; // doubling				
			kOne = 2*kOne;
		}
		
		// Relax on the coarsest grid
		// FMG
		// 2) Relax on the coarsest grid
		
		// move to the coarsest + 1
		iOne = iOne/2;
		jOne = jOne/2;
		kOne = kOne/2;
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		tphislices = kOne == 1 ? phislices : phislices/kOne;
		tphislices = tphislices < nnphis ? nnphis : tphislices;
		otphislices = tphislices;
			
			
		h = gridSizeR *  iOne;
		h2 = h*h;		
		gridSizePhi =  TMath::TwoPi()/tphislices;  // h_{phi}		
		tempratioPhi    =  h*h / (gridSizePhi*gridSizePhi) ;  // ratio_{phi} = gridsize_{r} / gridsize_{phi}		
		tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;

		for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
		}
	
		
	
		// 3) Relax on the coarsest grid
		Relax3D(tvArrayV[loops - 1],tvChargeFMG[loops - 1],trows,tcolumns,tphislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
		
		
		
		// 4) V Cycle from coarsest to finest
		for (Int_t  count = loops-2 ; count >=0; count-- ) {			
			// move to finer grid
			
			
			coef1.clear();
			coef2.clear();
			coef3.clear();
			coef4.clear();
			icoef4.clear();
		
			iOne = iOne/2;
			jOne = jOne/2;
			kOne = kOne/2;
		
			trows = iOne == 1 ? rows : rows/iOne + 1;
			tcolumns = jOne == 1 ? columns : columns/jOne + 1;	
			tphislices = kOne == 1 ? phislices : phislices/kOne;
			tphislices = tphislices < nnphis ? nnphis : tphislices;
			
			
		
			// 4) a) interpolate from 2h --> h grid
			//AliInfo(Form("Interp3D, trows=%d, tcolumns=%d, newphislices=%d, oldphislices=%d\n",trows,tcolumns,tphislices,otphislices));
			
			Interp3D(tvArrayV[count],tvArrayV[count+1],trows, tcolumns,tphislices,otphislices);
			
			// Copy the relax charge to the tvCharge 
			if (count > 0) {
				for (Int_t m=0;m<tphislices;m++) {
					*tvCharge[count][m] = *tvChargeFMG[count][m]; //copy
				}				
			}
			
		
			// do mgcycle in NcMGCycle times
			for (Int_t mgcycle=0;mgcycle < fMgParameters.NMGCYCLE ;mgcycle++) {
			
				// copy to store previous potential
				for (Int_t m=0;m<tphislices;m++) {
					*tvPrevArrayV[count][m] = *tvArrayV[count][m]; //copy
				}	
			
				VCycle3D(rows,columns,phislices,symmetry,count+1,loops,fMgParameters.NPRE, fMgParameters.NPOST,gridSizeR, ratioZ,tvArrayV, 
				  tvCharge, tvResidu,coef1,coef2,coef3,coef4,icoef4);
				
				
				/// converge error
				convErr = GetErrorConv(tvArrayV[count],tvPrevArrayV[count],tphislices);
				//// error counting /////
				if (count == 0){
					(*fErrorConvergenceNormInf)(mgcycle) = convErr;
					(*fError)(mgcycle) = GetErrorExact(arrayofArrayV,tvPrevArrayV[count],phislices);
				}
				/// if already converge just break move to finer grid
				if (convErr <= fgConvErr) {
					//AliInfo(Form("Conv Err: %f, MG Iteration : %d", convErr, mgcycle));
					fIterations = mgcycle+1;
					break;				
				}
			}		
			// keep old slice information
			otphislices = tphislices;
		}			
		
	} else if (fMgParameters.cycleType == kVCycle)  {
		// V-cycle
		Int_t gridfrom = 1;
		Int_t gridto = loops;	
		
		
	
		for (Int_t mgcycle=0;mgcycle < fMgParameters.NMGCYCLE;mgcycle++) {
			// copy to store previous potential
			for (Int_t m=0;m<phislices;m++) {
				*tvPrevArrayV[0][m] = *tvArrayV[0][m]; //copy
			}				
			
			// Do V Cycle from the coarsest to finest grid	
			
			VCycle3D(rows,columns,phislices, symmetry,gridfrom,gridto,fMgParameters.NPRE, fMgParameters.NPOST,gridSizeR, ratioZ,tvArrayV, tvCharge, tvResidu,
			coef1,coef2,coef3,coef4,icoef4);
		
					
			// convergence error
			convErr = GetErrorConv(tvArrayV[0],tvPrevArrayV[0],phislices);
			
			(*fErrorConvergenceNormInf)(mgcycle) = convErr;
			(*fError)(mgcycle) = GetErrorExact(arrayofArrayV,tvPrevArrayV[0],phislices);
					
			// if error already achieved then stop mg iteration
			if (convErr <= fgConvErr) {
				//AliInfo(Form("Exact Err: %f, MG Iteration : %d", (*fError)(mgcycle), mgcycle));
				fIterations = mgcycle+1;
				break;
			}
					
			
		}
	}


	// dealocate memory for multigrid
	for (Int_t count = 1 ; count <= loops ; count++ ) {		
		delete[] tvResidu[count - 1];
		delete[] tvPrevArrayV[count - 1];
		if (count > 1) {
			delete[] tvChargeFMG[count - 1];
			delete[] tvArrayV[count - 1];
			delete[] tvCharge[count - 1];
		}	
	}
}			
		

/// Helper function to check if the integer is equal to a power of two
/// \param i Int_t the number
/// \return 1 if it is a power of two, else 0
Int_t AliTPCPoissonSolver::IsPowerOfTwo(Int_t i) const {
  Int_t j = 0;
  while( i > 0 ) { j += (i&1) ; i = (i>>1) ; }
  if ( j == 1 ) return(1) ;  // True
  return(0) ;                // False
}


/// Relax3D
///
///    Relaxation operation for multigrid
///		 relaxation used 7 stencil in cylindrical coordinate
///
/// Using the following equations
/// \f$ U_{i,j,k} = (1 + \frac{1}{r_{i}h_{r}}}) U_{i+1,j,k}  + (1 - \frac{1}{r_{i}h_{r}}}) U_{i+1,j,k}  \f$
///
/// \param curArrayV TMatrixD** potential in 3D (matrices of matrix)
/// \param curCharge TMatrixD** charge in 3D
/// \param rows const Int_t number of rows in the r direction of TPC
/// \param columns const Int_t number of columns in z direction of TPC
/// \param phislices const Int_t number of phislices in phi direction of TPC
/// \param symmetry const Int_t is the cylinder has symmetry
/// \param h2 const Float_t \f$  h_{r}^{2} \f$
/// \param tempRatioZ const Float_t ration between grid size in z-direction and r-direction
/// \param coef1 std::vector<float> coef for \f$  V_{x+1,y,z} \f$
/// \param coef2 std::vector<float> coef for \f$  V_{x-1,y,z} \f$
/// \param coef3 std::vector<float> coef for z 
/// \param coef4 std::vector<float> coef for f(r,\phi,z)
///
void AliTPCPoissonSolver::Relax3D(TMatrixD **curArrayV,TMatrixD **curCharge,const Int_t trows,const Int_t tcolumns, 
		const Int_t phislices, const Int_t symmetry,const Float_t h2, 
		const Float_t tempRatioZ,std::vector<float> &coef1,std::vector<float> &coef2, std::vector<float> &coef3,std::vector<float> &coef4) {
			
	Int_t mplus,mminus, signplus,signminus;			
	TMatrixD * arrayV;
	TMatrixD * arrayVP;
	TMatrixD * arrayVM;
	TMatrixD * arrayCharge;		
	
	// Gauss-Seidel (Read Black}
	if (fMgParameters.relaxType == kGaussSeidel) {
		// for each slice
		Int_t isw,jsw,msw;	
		msw = 1;		
		for ( Int_t ipass=1;ipass<=2;ipass++,msw = 3 - msw) {	
			jsw = msw;
			for ( Int_t m = 0 ; m < phislices ; m++, jsw = 3 - jsw ) {				
				mplus  = m + 1;   signplus  = 1 ;
				mminus = m - 1 ;  signminus = 1 ;			
				// Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
				if (symmetry==1) {  
					if ( mplus  > phislices-1 ) mplus  = phislices - 2 ;
					if ( mminus < 0 )           mminus = 1 ;
				}		
				// Anti-symmetry in phi
				else if (symmetry==-1) {   
					if ( mplus  > phislices-1 ) { mplus  = phislices - 2 ; signplus  = -1 ; }
					if ( mminus < 0 )           { mminus = 1 ;	         signminus = -1 ; }
				}
				else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
					if ( mplus  > phislices-1 ) mplus  = m + 1 - phislices ;
					if ( mminus < 0 )           mminus = m - 1 + phislices ;
				}
				arrayV    =  curArrayV[m] ;
				arrayVP   =  curArrayV[mplus] ; // slice
				arrayVM   =  curArrayV[mminus] ; // slice
				arrayCharge =  curCharge[m] ;

				isw = jsw;
				for ( Int_t j = 1; j < tcolumns-1 ; j++, isw = 3-isw) {
					for ( Int_t i =  isw; i < trows-1 ; i += 2 ) {	
						//AliInfo(Form("Doing slice %d, z=%d, r=%d", m,j,i));			
						(*arrayV)(i,j) = ( coef2[i]          *   (*arrayV)(i-1,j)
							+ tempRatioZ        * ( (*arrayV)(i,j-1)  +  (*arrayV)(i,j+1) )						
							+ coef1[i]          *   (*arrayV)(i+1,j)
							+ coef3[i]          * ( signplus* (*arrayVP)(i,j)       +  signminus* (*arrayVM)(i,j) )
							+ (h2 * (*arrayCharge)(i,j))
							) * coef4[i] ;						
					} // end cols
				}	// end rows
			} // end phi
		} // end sweep		
	} 	
else if (fMgParameters.relaxType == kJacobi ) {
		// for each slice
		for ( Int_t m = 0 ; m < phislices ; m++ ) {

			mplus  = m + 1;   signplus  = 1 ;
			mminus = m - 1 ;  signminus = 1 ;
			
			// Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
			if (symmetry==1) {  
				if ( mplus  > phislices-1 ) mplus  = phislices - 2 ;
				if ( mminus < 0 )           mminus = 1 ;
			}		
			// Anti-symmetry in phi
			else if (symmetry==-1) {   
				if ( mplus  > phislices-1 ) { mplus  = phislices - 2 ; signplus  = -1 ; }
				if ( mminus < 0 )           { mminus = 1 ;	         signminus = -1 ; }
			}
			else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
				if ( mplus  > phislices-1 ) mplus  = m + 1 - phislices ;
				if ( mminus < 0 )           mminus = m - 1 + phislices ;
			}

			
			arrayV    =  curArrayV[m] ;
			arrayVP   =  curArrayV[mplus] ; // slice
			arrayVM   =  curArrayV[mminus] ; // slice
			arrayCharge =  curCharge[m] ;
		
			// Jacobian
			for ( Int_t j = 1; j < tcolumns-1 ; j++) {
				for ( Int_t i =  1; i < trows-1 ; i++) {
					(*arrayV)(i,j) = ( coef2[i]          *   (*arrayV)(i-1,j)
						+ tempRatioZ        * ( (*arrayV)(i,j-1)  +  (*arrayV)(i,j+1) )						
						+ coef1[i]          *   (*arrayV)(i+1,j)
						+ coef3[i]          * ( signplus* (*arrayVP)(i,j)       +  signminus* (*arrayVM)(i,j) )
						+ (h2 * (*arrayCharge)(i,j))
					) * coef4[i] ;
					
					
				} // end cols
			}	// end rows
			
			
		} // end phi
		
	}	 	else  {
		// Case weighted Jacobi
		// TODO
	}
	
}




/// Relax2D
///
///    Relaxation operation for multigrid
///		 relaxation used 5 stencil in cylindrical coordinate
///
/// Using the following equations
/// \f$ U_{i,j,k} = (1 + \frac{1}{r_{i}h_{r}}}) U_{i+1,j,k}  + (1 - \frac{1}{r_{i}h_{r}}}) U_{i+1,j,k}  \f$
///
/// \param curArrayV TMatrixD& potential in 3D (matrices of matrix)
/// \param curCharge TMatrixD& charge in 3D
/// \param rows const Int_t number of rows in the r direction of TPC
/// \param columns const Int_t number of columns in z direction of TPC
/// \param phislices const Int_t number of phislices in phi direction of TPC
/// \param symmetry const Int_t is the cylinder has symmetry
/// \param h2 const Float_t \f$  h_{r}^{2} \f$
/// \param tempFourth const Float_t coef for h
/// \param tempRatio const Float_t ratio between grid size in z-direction and r-direction
/// \param coef1 std::vector<float> coef for \f$  V_{x+1,y,z} \f$
/// \param coef2 std::vector<float> coef for \f$  V_{x-1,y,z} \f$
///
void AliTPCPoissonSolver::Relax2D(TMatrixD &curArrayV,TMatrixD &curCharge,const Int_t trows,const Int_t tcolumns,const Float_t h2, const Float_t tempFourth, const Float_t tempRatio,std::vector<float> &coef1,std::vector<float> &coef2) {
			
		// Gauss-Seidel
		if (fMgParameters.relaxType == kGaussSeidel ) {
			
			Int_t isw,jsw=1;
			for ( Int_t ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
				isw = jsw;		
				for ( Int_t j = 1; j < tcolumns-1 ; j++, isw = 3-isw) {
					for ( Int_t i =  isw; i < trows-1 ; i += 2 ) {
							curArrayV(i,j) =  tempFourth *(coef1[i] * curArrayV(i+1,j) +  coef2[i]  * curArrayV(i-1,j) 
									+ tempRatio * (curArrayV(i,j+1) + curArrayV(i,j-1)) + (h2 * curCharge(i,j)));						
					} // end cols
				}	// end rows
			} // end pass red-black	
		} else if (fMgParameters.relaxType == kJacobi )  {
				// Jacobian
			for ( Int_t j = 1; j < tcolumns-1 ; j++) {
				for ( Int_t i =  1; i < trows-1 ; i++) {
						curArrayV(i,j) =  tempFourth *(coef1[i] * curArrayV(i+1,j) +   coef2[i] * curArrayV(i-1,j) 
								+ tempRatio * (curArrayV(i,j+1) + curArrayV(i,j-1)) + (h2 * curCharge(i,j)));						
				} // end cols
			}	// end rows
		} else if (fMgParameters.relaxType == kWeightedJacobi){
				// Weighted Jacobi
				// TODO
		}
}




/// Residu3D
///
///    Compute residu from V(.) where V(.) is numerical potential and f(.).
///		 residu used 7 stencil in cylindrical coordinate
///
/// Using the following equations
/// \f$ U_{i,j,k} = (1 + \frac{1}{r_{i}h_{r}}}) U_{i+1,j,k}  + (1 - \frac{1}{r_{i}h_{r}}}) U_{i+1,j,k}  \f$
///
/// \param residu TMatrixD** residu in 3D (matrices of matrix)
/// \param curArrayV TMatrixD** potential in 3D (matrices of matrix)
/// \param curCharge TMatrixD** charge in 3D
/// \param rows const Int_t number of rows in the r direction of TPC
/// \param columns const Int_t number of columns in z direction of TPC
/// \param phislices const Int_t number of phislices in phi direction of TPC
/// \param symmetry const Int_t is the cylinder has symmetry
/// \param ih2 const Float_t \f$ 1/ h_{r}^{2} \f$
/// \param tempRatioZ const Float_t ration between grid size in z-direction and r-direction
/// \param coef1 std::vector<float> coef for \f$  V_{x+1,y,z} \f$
/// \param coef2 std::vector<float> coef for \f$  V_{x-1,y,z} \f$
/// \param coef3 std::vector<float> coef for z 
/// \param icoef4 std::vector<float> inverse coef for f(r,\phi,z)
///
void AliTPCPoissonSolver::Residu3D(TMatrixD **residu,TMatrixD **curArrayV,TMatrixD **curCharge,const Int_t trows,const Int_t tcolumns,const Int_t phislices, const Int_t symmetry, const Float_t ih2, 
		const Float_t tempRatioZ,std::vector<float> &coef1,std::vector<float> &coef2,std::vector<float> &coef3,std::vector<float> &icoef4) {
	Int_t mplus,mminus, signplus,signminus;		
	for ( Int_t m = 0 ; m < phislices ; m++ ) {

		mplus  = m + 1;   signplus  = 1 ;
		mminus = m - 1 ;  signminus = 1 ;
			
		// Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
		if (symmetry==1) {  
			if ( mplus  > phislices-1 ) mplus  = phislices - 2 ;
			if ( mminus < 0 )           mminus = 1 ;
		}		
		// Anti-symmetry in phi
		else if (symmetry==-1) {   
			if ( mplus  > phislices-1 ) { mplus  = phislices - 2 ; signplus  = -1 ; }
			if ( mminus < 0 )           { mminus = 1 ;	         signminus = -1 ; }
		}
		else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
			if ( mplus  > phislices-1 ) mplus  = m + 1 - phislices ;
			if ( mminus < 0 )           mminus = m - 1 + phislices ;
		}
		
		//TMatrixD& aresidu
		TMatrixD& arrayResidu = *residu[m];
		TMatrixD& arrayV    =  *curArrayV[m] ;
		TMatrixD& arrayVP   =  *curArrayV[mplus] ; // slice
		TMatrixD& arrayVM   =  *curArrayV[mminus] ; // slice
		TMatrixD& arrayCharge =  *curCharge[m] ;
		
		for ( Int_t j = 1; j < tcolumns-1 ; j++) {
			for ( Int_t i =  1; i < trows-1 ; i++) {		
				
				arrayResidu(i,j) = ih2 * (coef2[i] *   arrayV(i-1,j) + tempRatioZ * ( arrayV(i,j-1)  +  arrayV(i,j+1) )
						+ coef1[i] * arrayV(i+1,j) + coef3[i] * (signplus*arrayVP(i,j)  +  signminus*arrayVM(i,j))  - icoef4[i] * arrayV(i,j)) 
						+	arrayCharge(i,j);
			
			} // end cols
		}	// end rows
		
		//arrayResidu.Print();
	}
	
}

/// Residu2D
///
///    Compute residu from V(.) where V(.) is numerical potential and f(.).
///		 residu used 5 stencil in cylindrical coordinate
///
/// Using the following equations
/// \f$ U_{i,j,k} = (1 + \frac{1}{r_{i}h_{r}}}) U_{i+1,j,k}  + (1 - \frac{1}{r_{i}h_{r}}}) U_{i+1,j,k}  \f$
///
/// \param residu TMatrixD& potential in 2D 
/// \param curArrayV TMatrixD& potential in 2D 
/// \param curCharge TMatrixD& charge in 2D
/// \param rows const Int_t number of rows in the r direction of TPC
/// \param columns const Int_t number of columns in z direction of TPC
/// \param phislices const Int_t number of phislices in phi direction of TPC
/// \param symmetry const Int_t is the cylinder has symmetry
/// \param h2 const Float_t \f$  h_{r}^{2} \f$
/// \param tempFourth const Float_t coef for h
/// \param tempRatio const Float_t ratio between grid size in z-direction and r-direction
/// \param coef1 std::vector<float> coef for \f$  V_{x+1,y,z} \f$
/// \param coef2 std::vector<float> coef for \f$  V_{x-1,y,z} \f$
///
void AliTPCPoissonSolver::Residu2D(TMatrixD &residu,TMatrixD &curArrayV,TMatrixD &curCharge,const Int_t trows,const Int_t tcolumns,const Float_t ih2, const Float_t itempFourth, const Float_t tempRatio,std::vector<float> &coef1,std::vector<float> &coef2) {
		for ( Int_t i = 1; i < trows -1; i++) {
			for ( Int_t j =  1; j < tcolumns -1; j++ ) {
				residu(i,j) = ih2 * (coef1[i] * curArrayV(i+1,j) +   coef2[i] * curArrayV(i-1,j) 
									+ tempRatio * (curArrayV(i,j+1) + curArrayV(i,j-1)) - itempFourth * curArrayV(i,j)) + 
									curCharge(i,j);
		
			} // end cols
		}	// end rows
		
		//Boundary points.
		for (Int_t i=0;i<trows;i++) 
				residu(i,0) =residu(i,tcolumns-1) = 0.0;
				
		for (Int_t j=0;j<tcolumns;j++) 
				residu(0,j) =residu(trows-1,j) = 0.0;
		
}



/// Restrict2D
///
///    Grid transfer operator, restrict from fine -> coarse grid
///		 provide full-half weighting
///		
///		\f$ \[ \frac{1}{16}\left( \begin{array}{ccc}
///      1 & 2 & 1 
///      2 & 4 & 2
///      1 & 2 & 1 \end{array} \right) \]  \f$ 
///
/// \param curCharge TMatrixD& coarse grid (2h)
/// \param residu TMatrixD& fine grid  (h)
/// \param rows const Int_t number of rows in the r direction of TPC
/// \param columns const Int_t number of columns in z direction of TPC
///
void AliTPCPoissonSolver::Restrict2D(TMatrixD &curCharge,TMatrixD &residu,const Int_t trows,const Int_t tcolumns) {
		
		for (Int_t i = 1, ii=2; i < trows-1; i++,ii+=2) {
			for  ( Int_t j = 1, jj = 2; j < tcolumns-1; j++,jj+=2){
				if (fMgParameters.gtType ==  kHalf) {		
					// half
						curCharge(i,j) =  0.5 * residu(ii,jj) + 
															0.125 * (residu(ii+1,jj) + residu(ii-1,jj) +  residu(ii,jj+1) +  residu(ii,jj-1));
					
									
				} else 
				// full
						if (fMgParameters.gtType ==  kFull) {
						curCharge(i,j) =  0.25 * residu(ii,jj) + 
															 0.125 * (residu(ii+1,jj) + residu(ii-1,jj) +  residu(ii,jj+1) +  residu(ii,jj-1)) +
															 0.0625 * (residu(ii+1,jj+1) + residu(ii-1,jj+1) +  residu(ii+1,jj-1) +  residu(ii-1,jj-1));
				
				}
				
			} // end cols
		}	// end rows	
		
		// boundary
		// for boundary
		for ( Int_t j = 0, jj =0; j < tcolumns ; j++,jj+=2) {
			curCharge(0,j) = residu(0,jj);
			curCharge(trows - 1,j) = residu((trows - 1)*2,jj);
		}
			
			
		// for boundary
		for ( Int_t i = 0, ii =0; i < tcolumns ; i++,ii+=2) {
			curCharge(i,0) = residu(ii,0);
			curCharge(i,tcolumns-1) = residu(ii,(tcolumns - 1)*2);
		}
		
		
}


/// RestrictBoundary2D
///
///    Boundary transfer  restrict from fine -> coarse grid
///
/// \param curCharge TMatrixD& coarse grid (2h)
/// \param residu TMatrixD& fine grid  (h)
/// \param rows const Int_t number of rows in the r direction of TPC
/// \param columns const Int_t number of columns in z direction of TPC
///
void AliTPCPoissonSolver::RestrictBoundary2D(TMatrixD &curCharge,TMatrixD &residu,const Int_t trows,const Int_t tcolumns) {
		// for boundary
		for ( Int_t j = 0, jj =0; j < tcolumns ; j++,jj+=2) {
			curCharge(0,j) = residu(0,jj);
			curCharge(trows - 1,j) = residu((trows - 1)*2,jj);
		}
				
		// for boundary
		for ( Int_t i = 0, ii =0; i < tcolumns ; i++,ii+=2) {
			curCharge(i,0) = residu(ii,0);
			curCharge(i,tcolumns-1) = residu(ii,(tcolumns - 1)*2);
		}	
}


/// Restriction in 3D
///
/// Restriction is a map from fine grid (h) to coarse grid (2h)
///
/// In case of 3D
/// Full weighting:
/// \f[ (R u)_{i,j,k} = \frac{1}{2} u_{2i,2j,2k} + \frac{1}{4} S_{1} + \frac{1}{8} S_{2} + \frac{1}{16} S_{3}\f]
///
///
/// Restriction in all direction r-phi-z
/// restriction in phi only if oldphi == 2*newphi
/// \param curCharge TMatrixD** coarser grid 2h
/// \param residu TMatrixD ** fine grid h
/// \param trows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tcolumns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param newphislices Int_t number of phislices (in phi-direction) for coarser grid
/// \param oldphislices Int_t number of phislices (in phi-direction) for finer grid
///
void AliTPCPoissonSolver::Restrict3D(TMatrixD **curCharge,TMatrixD **residu,const Int_t trows,const Int_t tcolumns, const Int_t newphislices, const Int_t oldphislices) {
	
	Double_t s1,s2,s3;
	
	
	if(2*newphislices ==oldphislices) {	
		
		Int_t mplus,mminus;		
		Int_t mm=0;
				
		for ( Int_t m = 0 ; m < newphislices ; m++, mm+=2 ) {

			// assuming no symmetry
			mplus  = mm + 1;  
			mminus = mm - 1 ; 
			
			if ( mplus  > (oldphislices)-1 ) mplus  = mm + 1 - (oldphislices) ;
			if ( mminus < 0 )           mminus = mm - 1 + (oldphislices) ;
		
			TMatrixD& arrayResidu 	=  *residu[mm];
			TMatrixD& arrayResiduP  =  *residu[mplus] ;
			TMatrixD& arrayResiduM  =  *residu[mminus] ; // slice
			TMatrixD& arrayCharge 	=  *curCharge[m] ;
		
			
			for ( Int_t i =  1, ii=2; i < trows-1 ; i++,ii+=2) {		
				for ( Int_t j = 1, jj =2; j < tcolumns-1 ; j++,jj+=2) {
	
				
					// at the same plane
					s1 = arrayResidu(ii+1,jj) +  arrayResidu(ii-1,jj) +  arrayResidu(ii,jj+1) +  arrayResidu(ii,jj-1) + arrayResiduP(ii,jj) + arrayResiduM(ii,jj); 
					s2 = (arrayResidu(ii+1,jj+1)  +  arrayResidu(ii+1,jj-1) 	+ arrayResiduP(ii+1,jj) +  arrayResiduM(ii+1,jj)) +
							 (arrayResidu(ii-1,jj-1)  + arrayResidu(ii-1,jj+1)+ arrayResiduP(ii-1,jj) + arrayResiduM(ii-1,jj)) +  
								arrayResiduP(ii,jj-1) +  	arrayResiduM(ii,jj+1) + arrayResiduM(ii,jj-1) + arrayResiduP(ii,jj+1); 
								
					s3 = (arrayResiduP(ii+1,jj+1) + arrayResiduP(ii+1,jj-1)  + arrayResiduM(ii+1,jj+1) + arrayResiduM(ii+1,jj-1)) +
							 (arrayResiduM(ii-1,jj-1) + arrayResiduM(ii-1,jj+1) + arrayResiduP(ii-1,jj-1) + arrayResiduP(ii-1,jj+1));
								
					arrayCharge(i,j) = 0.125 * arrayResidu(ii,jj) + 0.0625 * s1 +  0.03125 * s2 + 0.015625 * s3;
				} // end cols
			}	// end rows		
			
			// for boundary
			for ( Int_t j = 0, jj =0; j < tcolumns ; j++,jj+=2) {
					arrayCharge(0,j) = arrayResidu(0,jj);
					arrayCharge(trows - 1,j) = arrayResidu((trows - 1)*2,jj);
			}
			
			
			// for boundary
			for ( Int_t i = 0, ii =0; i < tcolumns ; i++,ii+=2) {
					arrayCharge(i,0) = arrayResidu(ii,0);
					arrayCharge(i,tcolumns-1) = arrayResidu(ii,(tcolumns - 1)*2);
			}
		}// end phis
			
	} else {		
		for (int m=0;m<newphislices;m++) {
			Restrict2D(*curCharge[m],*residu[m],trows,tcolumns);
		}		
	}
	
}



/// Restrict Boundary in 3D 
///
/// Pass boundary information to coarse grid
///
/// \param curCharge TMatrixD** coarser grid 2h
/// \param residu TMatrixD ** fine grid h
/// \param trows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tcolumns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param newphislices Int_t number of phislices (in phi-direction) for coarser grid
/// \param oldphislices Int_t number of phislices (in phi-direction) for finer grid
///
void AliTPCPoissonSolver::RestrictBoundary3D(TMatrixD **curCharge,TMatrixD **residu,const Int_t trows,const Int_t tcolumns, const Int_t newphislices, const Int_t oldphislices) {
	
	// in case of full 3d and the phislices is also coarsening
	
	
	if(2*newphislices ==oldphislices) {	
		
		for ( Int_t m = 0,mm=0 ; m < newphislices ; m++, mm+=2 ) {
		
			
			
			TMatrixD& arrayResidu 	=  *residu[mm];
			TMatrixD& arrayCharge 	=  *curCharge[m] ;
			// for boundary
			for ( Int_t j = 0, jj =0; j < tcolumns ; j++,jj+=2) {
					arrayCharge(0,j) = arrayResidu(0,jj);
					arrayCharge(trows - 1,j) = arrayResidu((trows - 1)*2,jj);
			}
			
			
			// for boundary
			for ( Int_t i = 0, ii =0; i < tcolumns ; i++,ii+=2) {
					arrayCharge(i,0) = arrayResidu(ii,0);
					arrayCharge(i,tcolumns-1) = arrayResidu(ii,(tcolumns - 1)*2);
			}
		}// end phis
		
		// dont forget to copy boundary 
			
			
	} else {
		
		
		for (int m=0;m<newphislices;m++) {
			RestrictBoundary2D(*curCharge[m],*residu[m],trows,tcolumns);
		}		
	}
	
}



/// Prolongation with Addition for 2D
///
/// Interpolation with addition from coarse level (2h) -->  fine level (h)
///
/// Interpolation in all direction r-phi-z
/// \param curArrayV TMatrixD& fine grid h
/// \param curArrayVC TMatrixD& coarse grid 2h
/// \param trows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tcolumns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1a
///
void AliTPCPoissonSolver::AddInterp2D(TMatrixD &curArrayV,TMatrixD &curArrayVC,const Int_t trows,const Int_t tcolumns) {
		for (Int_t j=2;j < tcolumns-1;j +=2 ) {
				for (Int_t i=2;i < trows-1;i += 2) {
					curArrayV(i,j) = curArrayV(i,j)  +  curArrayVC(i/2,j/2); 
				}
			}
		
		for (Int_t j=1;j < tcolumns-1;j +=2 ) {
			for (Int_t i=2;i < trows-1;i += 2) {
				curArrayV(i,j) = curArrayV(i,j)  +  0.5 * (curArrayVC(i/2,j/2) +curArrayVC(i/2,j/2+1)); 
			}
		}

		for (Int_t j=2;j < tcolumns-1;j +=2 ) {
			for (Int_t i=1;i < trows-1;i += 2) {
				curArrayV(i,j) = curArrayV(i,j)  +  0.5 * (curArrayVC(i/2,j/2) +curArrayVC(i/2+1,j/2)); 
			}
		}

		// only if full 
		if (fMgParameters.gtType == kFull ) {				
			for (Int_t j=1;j < tcolumns-1;j +=2 ) {
				for (Int_t i=1;i < trows-1;i += 2) {
					curArrayV(i,j) = curArrayV(i,j)  +  0.25 * (curArrayVC(i/2,j/2) + curArrayVC(i/2,j/2+1)+ curArrayVC(i/2+1,j/2) + curArrayVC(i/2+1,j/2+1)); 
			}
		}
	}
}

/// Prolongation with Addition for 3D
///
/// Interpolation with addition from coarse level (2h) -->  fine level (h) 
///
/// Interpolation in all direction r-phi-z
/// Interpolation in phi only if oldphi == 2*newphi
/// \param curArrayV TMatrixD& fine grid h
/// \param curArrayVC TMatrixD& coarse grid 2h
/// \param trows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tcolumns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1a
/// \param newphislices Int_t number of phislices (in phi-direction) for coarser grid
/// \param oldphislices Int_t number of phislices (in phi-direction) for finer grid
///
void AliTPCPoissonSolver::AddInterp3D(TMatrixD **curArrayV,TMatrixD **curArrayVC,const Int_t trows,const Int_t tcolumns, const Int_t newphislices, const Int_t oldphislices) {
	// Do restrict 2 D for each slice
	
	//const Float_t  h   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / ((trows-1)/2); // h_{r}  
	//Float_t radius,ratio;
	//std::vector<float> coef1((trows-1) / 2 );  // coef1(rows) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
	//std::vector<float> coef2((trows-1) / 2);  // coef2(rows) for storing (1 + h_{r}/2r_{i}) from central differences in r direction
	
	
	if (newphislices == 2*oldphislices) {
		Int_t mplus,mmplus;		
		Int_t mm=0;
		
		for ( Int_t m = 0 ; m < newphislices ; m += 2 ) {

			// assuming no symmetry
			mm = m/2;
			mmplus  = mm + 1; 
			mplus = m + 1;
			
			// round
			if ( mmplus  > (oldphislices)-1 ) mmplus  = mm + 1 - (oldphislices) ;			
			if ( mplus  > (newphislices)-1 ) mplus  = m + 1 - (newphislices) ;
		
			
			TMatrixD& fineV		=  *curArrayV[m];
			TMatrixD& fineVP	=  *curArrayV[mplus] ;
			TMatrixD& coarseV =  *curArrayVC[mm] ;
			TMatrixD& coarseVP =  *curArrayVC[mmplus] ;	
			
			
			for (Int_t j=2;j < tcolumns-1;j +=2 ) {
				for (Int_t i=2;i < trows-1;i += 2) {
					fineV(i,j) +=  coarseV(i/2,j/2); 				
					// point on corner lines at phi direction
					fineVP(i,j) += 0.5 * (coarseV(i/2,j/2) + coarseVP(i/2,j/2));								
				}
			}
		
			for (Int_t j=1;j < tcolumns-1;j +=2 ) {
				for (Int_t i=2;i < trows-1;i += 2) {
					fineV(i,j) +=   0.5 * (coarseV(i/2,j/2) + coarseV(i/2,j/2+1)); 					
					// point on corner lines at phi direction
					fineVP(i,j) +=  0.25 *(coarseV(i/2,j/2) + coarseV(i/2,j/2+1) + coarseVP(i/2,j/2) + coarseVP(i/2,j/2+1));
					
				}
			}

			for (Int_t j=2;j < tcolumns-1;j +=2 ) {
				for (Int_t i=1;i < trows-1;i += 2) {
					fineV(i,j) +=  0.5 * (coarseV(i/2,j/2) + coarseV(i/2+1,j/2)); 
					
					// point on line at phi direction
					fineVP(i,j) +=  0.25 *((coarseV(i/2,j/2)  + coarseVP(i/2,j/2)) + (coarseVP(i/2+1,j/2) + coarseV(i/2+1,j/2)));
					
				}
			}

			for (Int_t j=1;j < tcolumns-1;j +=2 ) {
				for (Int_t i=1;i < trows-1;i += 2) {
					fineV(i,j) +=   0.25 * ((coarseV(i/2,j/2) + coarseV(i/2,j/2+1)) + (coarseV(i/2+1,j/2) + coarseV(i/2+1,j/2+1))); 
					
					// point at the center at phi direction
					fineVP(i,j) +=   0.125 * ((coarseV(i/2,j/2) + coarseV(i/2,j/2+1)  + coarseVP(i/2,j/2) + coarseVP(i/2,j/2+1)) + (coarseV(i/2+1,j/2) + coarseV(i/2+1,j/2+1)+ coarseVP(i/2+1,j/2) + coarseVP(i/2+1,j/2+1))); 
				}
			}		
		}
		
		
	} else {
		for (int m=0;m<newphislices;m++) {
			AddInterp2D(*curArrayV[m],*curArrayVC[m],trows,tcolumns);		
		}		
	}
	
}



/// Interpolation/Prolongation in 3D
///
/// Interpolation is a map from coarse grid (h) to fine grid (2h)
///
/// In case of 3D
/// Full weighting:
/// \f[ (R u)_{i,j,k} = \frac{1}{2} u_{2i,2j,2k} + \frac{1}{4} S_{1} + \frac{1}{8} S_{2} + \frac{1}{16} S_{3}\f]
///
///
/// Restriction in all direction r-phi-z
/// restriction in phi only if oldphi == 2*newphi
/// \param curArrayV TMatrixD** finer grid h
/// \param curArrayCV TMatrixD ** coarse grid 2h
/// \param trows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tcolumns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param newphislices Int_t number of phislices (in phi-direction) for coarser grid
/// \param oldphislices Int_t number of phislices (in phi-direction) for finer grid
///
void AliTPCPoissonSolver::Interp3D(TMatrixD **curArrayV,TMatrixD **curArrayVC,const Int_t trows,const Int_t tcolumns, const Int_t newphislices, const Int_t oldphislices) {
	
	// Do restrict 2 D for each slice
	if (newphislices == 2*oldphislices) {
		Int_t mplus,mmplus;		
		Int_t mm=0;
		
		for ( Int_t m = 0 ; m < newphislices ; m += 2 ) {

			// assuming no symmetry
			mm = m/2;
			mmplus  = mm + 1; 
			mplus = m + 1;
			
			// round
			if ( mmplus  > (oldphislices)-1 ) mmplus  = mm + 1 - (oldphislices) ;			
			if ( mplus  > (newphislices)-1 ) mplus  = m + 1 - (newphislices) ;
		
			
			TMatrixD& fineV		=  *curArrayV[m];
			TMatrixD& fineVP	=  *curArrayV[mplus] ;
			TMatrixD& coarseV =  *curArrayVC[mm] ;
			TMatrixD& coarseVP =  *curArrayVC[mmplus] ;
			
			
			for (Int_t j=2;j < tcolumns-1;j +=2 ) {
				for (Int_t i=2;i < trows-1;i += 2) {
					fineV(i,j) =  coarseV(i/2,j/2); 
					
					// point on corner lines at phi direction
					fineVP(i,j) = 0.5 * (coarseV(i/2,j/2) + coarseVP(i/2,j/2));				
					
				}
			}
		
			for (Int_t j=1;j < tcolumns-1;j +=2 ) {
				for (Int_t i=2;i < trows-1;i += 2) {
					fineV(i,j) =   0.5 * (coarseV(i/2,j/2) + coarseV(i/2,j/2+1)); 
					
					// point on corner lines at phi direction
					fineVP(i,j) =  0.25 *(coarseV(i/2,j/2) + coarseV(i/2,j/2+1) + coarseVP(i/2,j/2) + coarseVP(i/2,j/2+1));
					
				}
			}

			for (Int_t j=2;j < tcolumns-1;j +=2 ) {
				for (Int_t i=1;i < trows-1;i += 2) {
					fineV(i,j) =  0.5 * ( coarseV(i/2,j/2) + coarseV(i/2+1,j/2)); 
					
					// point on line at phi direction
					fineVP(i,j) =  0.25 *((coarseV(i/2,j/2)  + coarseVP(i/2,j/2)) + (coarseVP(i/2+1,j/2) + coarseV(i/2+1,j/2)));
					
				}
			}

			for (Int_t j=1;j < tcolumns-1;j +=2 ) {
				for (Int_t i=1;i < trows-1;i += 2) {
					fineV(i,j) =   0.25 * ( (coarseV(i/2,j/2) + coarseV(i/2,j/2+1)) + (coarseV(i/2+1,j/2) + coarseV(i/2+1,j/2+1))); 
					
					// point at the center at phi direction
					fineVP(i,j) =   0.125 * ( (coarseV(i/2,j/2) + coarseV(i/2,j/2+1)  + coarseVP(i/2,j/2) + coarseVP(i/2,j/2+1)) +  (coarseV(i/2+1,j/2) + coarseV(i/2+1,j/2+1)+ coarseVP(i/2+1,j/2) + coarseVP(i/2+1,j/2+1))); 
				}
			}		
		}
		
		
	} else {
		for (int m=0;m<newphislices;m++) {
			Interp2D(*curArrayV[m],*curArrayVC[m],trows,tcolumns);		
		}		
	}
}



/// Interpolation/Prolongation in 2D
///
/// Interpolation is a map from coarse grid (h) to fine grid (2h)
///
/// In case of 2D
/// Full weighting:
/// \f[ (R u)_{i,j,k} = \frac{1}{2} u_{2i,2j,2k} + \frac{1}{4} S_{1} + \frac{1}{8} S_{2} + \frac{1}{16} S_{3}\f]
///
///
/// Restriction in all direction r-phi-z
/// \param curArrayV TMatrixD** finer grid h
/// \param curArrayCV TMatrixD ** coarse grid 2h
/// \param trows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param tcolumns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
///
void AliTPCPoissonSolver::Interp2D(TMatrixD &curArrayV,TMatrixD &curArrayVC,const Int_t trows,const Int_t tcolumns) {
		 for (Int_t j=2;j < tcolumns-1;j +=2 ) {
				for (Int_t i=2;i < trows-1;i += 2) {
					curArrayV(i,j) =  curArrayVC(i/2,j/2); 
				}
			}
		
		for (Int_t j=1;j < tcolumns-1;j +=2 ) {
			for (Int_t i=2;i < trows-1;i += 2) {
				curArrayV(i,j) =   0.5 * (curArrayVC(i/2,j/2) +curArrayVC(i/2,j/2+1)); 
			}
		}

		for (Int_t j=2;j < tcolumns-1;j +=2 ) {
			for (Int_t i=1;i < trows-1;i += 2) {
				curArrayV(i,j) =  0.5 * (curArrayVC(i/2,j/2) +curArrayVC(i/2+1,j/2)); 
			}
		}

		// only if full 
		if (fMgParameters.gtType == kFull) {		
			for (Int_t j=1;j < tcolumns-1;j +=2 ) {
				for (Int_t i=1;i < trows-1;i += 2) {
					curArrayV(i,j) =   0.25 * (curArrayVC(i/2,j/2) + curArrayVC(i/2,j/2+1)+ curArrayVC(i/2+1,j/2) + curArrayVC(i/2+1,j/2+1)); 
				}
			}
		}
		
}



/// V-Cycle 2D
///
/// Implementation non-recursive V-cycle for 2D
/// 
///	Algorithms:
///
/// \param rows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param columns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param gridfrom const Int_t finest level of grid
/// \param gridto const Int_t coarsest level of grid
/// \param NPRE const Int_t number of smoothing before coarsening
/// \param NPOST const Int_t number of smoothing after coarsening
/// \param gridSizeR const Float_t grid size in r direction (OPTION,  recalculate)
/// \param ratio const Float_t ratio between square of grid r and grid z (OPTION,  recalculate)
/// \param tvArrayV vector<TMatrixD *> vector of V potential in different grids
/// \param tvCharge vector<TMatrixD *> vector of charge distribution in different grids
/// \param tvResidu vector<TMatrixD *> vector of residu calculation in different grids
///
void AliTPCPoissonSolver::VCycle2D(const Int_t rows, const Int_t columns, const Int_t gridfrom,const  Int_t gridto,
				const Int_t NPRE, const Int_t NPOST, const  Float_t gridSizeR,const  Float_t ratio,std::vector<TMatrixD *> & tvArrayV,
				std::vector<TMatrixD *> & tvCharge,std::vector<TMatrixD *> & tvResidu) {

	Float_t h,h2,ih2,tempRatio,tempFourth,itempFourth,radius;
	TMatrixD *curArrayV,*curArrayVC;
	TMatrixD *curCharge;
	TMatrixD *residu;
	Int_t iOne,jOne, trows, tcolumns, count;
	iOne = 1 << (gridfrom-1); 
	jOne = 1 << (gridfrom-1); 	
	
	
	curArrayV = NULL;
	curArrayVC = NULL;
	curCharge = NULL;
	residu = NULL;
	
	trows = iOne == 1 ? rows : rows/iOne + 1;
	tcolumns = jOne == 1 ? columns : columns/jOne + 1;
	
	// allocate vector for coefessions
	std::vector<float> coef1(rows);
	std::vector<float> coef2(columns);		
	
	// 1) Go to coarsest level
	for (  count = gridfrom ; count <= gridto-1 ; count++ ) {		
		h = gridSizeR *  iOne;
		h2 = h*h;
		ih2 = 1.0/h2;
		tempRatio     = ratio * iOne * iOne / ( jOne * jOne ) ;
		tempFourth    = 1.0 / (2.0 + 2.0*tempRatio) ;
		itempFourth 	= 1.0 / tempFourth;			
		for ( Int_t i = 1 ; i < trows-1 ; i++ ) {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h ;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);			
		}		
		curArrayV = tvArrayV[count - 1];
		curCharge = tvCharge[count - 1];		
		residu = tvResidu[count - 1];	
		
		// 1) Pre-Smoothing: Gauss-Seidel Relaxation or Jacobi		
		for (Int_t jpre=1;jpre <= NPRE;jpre++) {				
			Relax2D(*curArrayV,*curCharge,trows,tcolumns,h2,tempFourth,tempRatio,coef1,coef2);			
		} 
		
		// 2) Residu calculation 
		Residu2D(*residu,*curArrayV,*curCharge,trows,tcolumns,ih2,itempFourth,tempRatio,coef1,coef2);			
					
		iOne = 2*iOne;
		jOne = 2*jOne;
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		
		curCharge =  tvCharge[count];		
		curArrayV = tvArrayV[count];
		
		//3) Restriction 
		Restrict2D(*curCharge,*residu,trows,tcolumns);			
		
		//4) Zeroing coarser V
		curArrayV->Zero();		
	}
	
	
	// 5) coarsest grid
	h = gridSizeR *  iOne;
	h2 = h*h;
	tempRatio     = ratio * iOne * iOne / ( jOne * jOne ) ;
	tempFourth    = 1.0 / (2.0 + 2.0*tempRatio) ;
	
	for ( Int_t i = 1 ; i < trows-1 ; i++ ) {
		radius = AliTPCPoissonSolver::fgkIFCRadius + i*h ;
		coef1[i] = 1.0 + h/(2*radius);
		coef2[i] = 1.0 - h/(2*radius);		
	}
	
	Relax2D(*curArrayV,*curCharge,trows,tcolumns,h2,tempFourth,tempRatio,coef1,coef2);			
	
	// Go to finest grid
	for (  count = gridto-1 ; count >= gridfrom ; count-- ) {
		
		iOne = iOne/2;
		jOne = jOne/2;
		
		h = gridSizeR * iOne;
		h2 = h*h;
		ih2 = 1.0/h2;
		tempRatio     = ratio * iOne * iOne / ( jOne * jOne ) ;
		tempFourth    = 1.0 / (2.0 + 2.0*tempRatio) ;
		itempFourth 	= 1.0 / tempFourth;
	
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		curCharge =  tvCharge[count-1];				
		curArrayV = tvArrayV[count-1];
		curArrayVC = tvArrayV[count];
	
		
		
		// 6) Interpolation/Prolongation
		AddInterp2D(*curArrayV,*curArrayVC,trows,tcolumns);
			
						
		for ( Int_t i = 1 ; i < trows-1 ; i++ ) {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);			
		}
		
		// 7) Post-Smoothing: Gauss-Seidel Relaxation
		for (Int_t jpost=1;jpost <= NPOST;jpost++) {			
			Relax2D(*curArrayV,*curCharge,trows,tcolumns,h2,tempFourth,tempRatio,coef1,coef2);									
		} // end post smoothing	
		
		//// DEBUG ////
		//AliInfo(Form("Count %d", count));
		//AliInfo(Form("Exact Err: %f, MG Iteration : %d", (*fError)(mgcycle), mgcycle));
		//curArrayV->Print();
		//curCharge->Print();
	}	
}
	


/// W-Cycle 2D
///
/// Implementation non-recursive W-cycle for 2D
/// 
///	Algorithms:
///
/// \param rows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param columns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param gridfrom const Int_t finest level of grid
/// \param gridto const Int_t coarsest level of grid
/// \param gamma const Int_t number of iterations at coarsest level
/// \param NPRE const Int_t number of smoothing before coarsening
/// \param NPOST const Int_t number of smoothing after coarsening
/// \param gridSizeR const Float_t grid size in r direction (OPTION,  recalculate)
/// \param ratio const Float_t ratio between square of grid r and grid z (OPTION,  recalculate)
/// \param tvArrayV vector<TMatrixD *> vector of V potential in different grids
/// \param tvCharge vector<TMatrixD *> vector of charge distribution in different grids
/// \param tvResidu vector<TMatrixD *> vector of residu calculation in different grids
///
void AliTPCPoissonSolver::WCycle2D(const Int_t rows, const Int_t columns, const Int_t gridfrom,const  Int_t gridto, const int gamma,
				const Int_t NPRE, const Int_t NPOST, const  Float_t gridSizeR,const  Float_t ratio,std::vector<TMatrixD *> & tvArrayV,
				std::vector<TMatrixD *> & tvCharge,std::vector<TMatrixD *> & tvResidu) {

	Float_t h,h2,ih2,tempRatio,tempFourth,itempFourth,radius;
	TMatrixD *curArrayV,*curArrayVC;
	TMatrixD *curCharge;
	TMatrixD *residu;
	Int_t iOne,jOne, trows, tcolumns, count;
	iOne = 1 << (gridfrom-1); 
	jOne = 1 << (gridfrom-1); 	
	
	trows = iOne == 1 ? rows : rows/iOne + 1;
	tcolumns = jOne == 1 ? columns : columns/jOne + 1;
	
	// allocate vector for coefessions
	std::vector<float> coef1(rows);
	std::vector<float> coef2(columns);		
	
	// 1) Go to coarsest level
	for (  count = gridfrom ; count <= gridto-2 ; count++ ) {		
		h = gridSizeR *  iOne;
		h2 = h*h;
		ih2 = 1.0/h2;
		tempRatio     = ratio * iOne * iOne / ( jOne * jOne ) ;
		tempFourth    = 1.0 / (2.0 + 2.0*tempRatio) ;
		itempFourth 	= 1.0 / tempFourth;			
		for ( Int_t i = 1 ; i < trows-1 ; i++ ) {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h ;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);			
		}		
		curArrayV = tvArrayV[count - 1];
		curCharge = tvCharge[count - 1];		
		residu = tvResidu[count - 1];	
		
		// 1) Pre-Smoothing: Gauss-Seidel Relaxation or Jacobi		
		for (Int_t jpre=1;jpre <= NPRE;jpre++) {				
			Relax2D(*curArrayV,*curCharge,trows,tcolumns,h2,tempFourth,tempRatio,coef1,coef2);			
		} 
		
		// 2) Residu calculation 
		Residu2D(*residu,*curArrayV,*curCharge,trows,tcolumns,ih2,itempFourth,tempRatio,coef1,coef2);			
					
		iOne = 2*iOne;
		jOne = 2*jOne;
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		
		curCharge =  tvCharge[count];		
		curArrayV = tvArrayV[count];
		
		//3) Restriction 
		Restrict2D(*curCharge,*residu,trows,tcolumns);			
		
		//4) Zeroing coarser V
		curArrayV->Zero();		
	}
	
	// Do V cylce from: gridto-1 to gridto gamma times
	for (Int_t igamma=0;igamma < gamma;igamma++) {
		VCycle2D(rows,columns,gridto - 1,gridto,
				NPRE,NPOST,gridSizeR,ratio,tvArrayV,
				tvCharge,tvResidu);
	}
	
	
		
	
	
	
	// Go to finest grid
	for (  count = gridto-2 ; count >= gridfrom ; count-- ) {
		
		iOne = iOne/2;
		jOne = jOne/2;
		
		h = gridSizeR * iOne;
		h2 = h*h;
		ih2 = 1.0/h2;
		tempRatio     = ratio * iOne * iOne / ( jOne * jOne ) ;
		tempFourth    = 1.0 / (2.0 + 2.0*tempRatio) ;
		itempFourth 	= 1.0 / tempFourth;
	
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		curCharge =  tvCharge[count-1];				
		curArrayV = tvArrayV[count-1];
		curArrayVC = tvArrayV[count];
	
		
		// 6) Interpolation/Prolongation
		AddInterp2D(*curArrayV,*curArrayVC,trows,tcolumns);
			
						
		for ( Int_t i = 1 ; i < trows-1 ; i++ ) {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);			
		}
		
		// 7) Post-Smoothing: Gauss-Seidel Relaxation
		for (Int_t jpost=1;jpost <= NPOST;jpost++) {			
			Relax2D(*curArrayV,*curCharge,trows,tcolumns,h2,tempFourth,tempRatio,coef1,coef2);									
		} // end post smoothing	
	}	
}


/// VCycle 3D2D, V Cycle 3D in multigrid with constant phislices
/// fine-->coarsest-->fine, propagating the residu to correct initial guess of V
///
/// Algorithm:
///
///    NOTE: In order for this algorith to work, the number of rows and columns must be a power of 2 plus one.
///    The number of rows and COLUMNS can be different.
///
///    ROWS       ==  2**M + 1
///    COLUMNS    ==  2**N + 1
///    PHISLICES  ==  Arbitrary but greater than 3
///
///    DeltaPhi in Radians
/// \param rows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param columns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param gridfrom const Int_t finest level of grid
/// \param gridto const Int_t coarsest level of grid
/// \param NPRE const Int_t number of smoothing before coarsening
/// \param NPOST const Int_t number of smoothing after coarsening
/// \param gridSizeR const Float_t grid size in r direction (OPTION,  recalculate)
/// \param ratio const Float_t ratio between square of grid r and grid z (OPTION,  recalculate)
/// \param tvArrayV vector<TMatrixD *> vector of V potential in different grids
/// \param tvCharge vector<TMatrixD *> vector of charge distribution in different grids
/// \param tvResidu vector<TMatrixD *> vector of residu calculation in different grids
/// \param coef1 std::vector<float>& coefissien for relaxation (r direction)
/// \param coef2 std::vector<float>& coefissien for relaxation (r direction)
/// \param coef3 std::vector<float>& coefissien for relaxation (ratio r/z)
/// \param coef4 std::vector<float>& coefissien for relaxation (ratio for grid_r)
/// \param icoef4 std::vector<float>& coefissien for relaxation (inverse coef4)
///
void AliTPCPoissonSolver::VCycle3D2D(const Int_t rows, const Int_t columns, const Int_t phislices, const Int_t symmetry,const Int_t gridfrom,const  Int_t gridto,const Int_t NPRE, const Int_t NPOST,const Float_t gridSizeR, const Float_t ratioZ, const Float_t ratioPhi,
		std::vector<TMatrixD **> & tvArrayV,	std::vector<TMatrixD **> & tvCharge,std::vector<TMatrixD **> & tvResidu,std::vector<float> &coef1,std::vector<float> &coef2, std::vector<float> &coef3,std::vector<float> &coef4,std::vector<float> &icoef4) {
					
	Float_t h,h2,ih2,tempratioZ, tempratioPhi,radius;
	TMatrixD **curArrayV,**curArrayVC;
	TMatrixD **curCharge;
	TMatrixD **residu;
	Int_t iOne,jOne, trows, tcolumns, count;
	
	curArrayV = NULL;
	curArrayVC = NULL;
	curCharge = NULL;
	residu = NULL;
	
	
	iOne = 1 << (gridfrom-1); 
	jOne = 1 << (gridfrom-1); 	
						
	trows = iOne == 1 ? rows : rows/iOne + 1;
	tcolumns = jOne == 1 ? columns : columns/jOne + 1;

	
	
	for (  count = gridfrom ; count <= gridto-1 ; count++ ) {		
		h   =  gridSizeR  * iOne ;
		h2 = h*h;
		ih2 = 1.0/h2;
			
		tempratioPhi    =  ratioPhi * iOne * iOne ; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
		tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;


		for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
			icoef4[i] = 1.0 / coef4[i];
		}
		
		curArrayV = tvArrayV[count - 1];
		curCharge = tvCharge[count - 1];		
		residu = tvResidu[count - 1];
	
		
		//AliInfo("Before Pre-smoothing");
		//curArrayV->Print();
		
		// 1) Pre-Smoothing: Gauss-Seidel Relaxation or Jacobi		
		for (Int_t jpre=1;jpre <= NPRE;jpre++) {				
			Relax3D(curArrayV,curCharge,trows,tcolumns,phislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
		} // end pre smoothing
		
		// 2) Residu calculation 
		Residu3D(residu,curArrayV,curCharge,trows,tcolumns,phislices,symmetry,ih2,tempratioZ,coef1,coef2,coef3,icoef4);			
				
		iOne = 2*iOne;
		jOne = 2*jOne;
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		
		curCharge =  tvCharge[count];		
		curArrayV = tvArrayV[count];
		
		//3) Restriction 
		//Restrict2D(*curCharge,*residu,trows,tcolumns);			
		Restrict3D(curCharge, residu, trows, tcolumns,  phislices, phislices);
		
		//4) Zeroing coarser V
		for(Int_t m=0;m<phislices;m++) {
				curArrayV[m]->Zero();		
		}		
	}
	
	// coarsest grid
	h   =  gridSizeR  * iOne ;
	h2 = h*h;
			
	tempratioPhi    =  ratioPhi * iOne * iOne ; 
	tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;

	
	for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
	}
	
	// 3) Relax on the coarsest grid
	Relax3D(curArrayV,curCharge,trows,tcolumns,phislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
	
	// back to fine
	for (  count = gridto-1 ; count >= gridfrom ; count-- ) {
		iOne = iOne/2;
		jOne = jOne/2;
		
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		
		h   =  gridSizeR  * iOne ;
		h2 = h*h;
			
		tempratioPhi    =  ratioPhi * iOne * iOne ; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
		tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;
	
		curCharge =  tvCharge[count-1];				
		curArrayV = tvArrayV[count-1];
		curArrayVC = tvArrayV[count];
	
		// 4) Interpolation/Prolongation
		AddInterp3D(curArrayV,curArrayVC,trows,tcolumns,phislices,phislices);
			

		
		for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
		}		
	
		
		// 5) Post-Smoothing: Gauss-Seidel Relaxation
		for (Int_t jpost=1;jpost <= NPOST;jpost++) {			
			Relax3D(curArrayV,curCharge,trows,tcolumns,phislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
		} // end post smoothing	
	}			
}
	



/// VCycle 3D, V Cycle in multigrid, fine-->coarsest-->fine, propagating the residu to correct initial guess of V
///
///    NOTE: In order for this algorith to work, the number of rows and columns must be a power of 2 plus one.
///    The number of rows and COLUMNS can be different.
///
///    ROWS       ==  2**M + 1
///    COLUMNS    ==  2**N + 1
///    PHISLICES  ==  Arbitrary but greater than 3
///
///    DeltaPhi in Radians
///
/// \param rows Int_t number of rows in the r direction of TPC
/// \param columns Int_t number of columns in z direction of TPC
/// \param phislices Int_t number of phislices in phi direction of T{C
/// \param symmetry Int_t symmetry or not
/// \param gridfrom const Int_t finest level of grid
/// \param gridto const Int_t coarsest level of grid
/// \param NPRE const Int_t number of smoothing before coarsening
/// \param NPOST const Int_t number of smoothing after coarsening
/// \param gridSizeR const Float_t grid size in r direction (OPTION,  recalculate)
/// \param ratioz const Float_t ratio between square of grid r and grid z (OPTION,  recalculate)
/// \param tvArrayV vector<TMatrixD *> vector of V potential in different grids
/// \param tvCharge vector<TMatrixD *> vector of charge distribution in different grids
/// \param tvResidu vector<TMatrixD *> vector of residu calculation in different grids
/// \param coef1 std::vector<float>& coefissien for relaxation (r direction)
/// \param coef2 std::vector<float>& coefissien for relaxation (r direction)
/// \param coef3 std::vector<float>& coefissien for relaxation (ratio r/z)
/// \param coef4 std::vector<float>& coefissien for relaxation (ratio for grid_r)
/// \param icoef4 std::vector<float>& coefissien for relaxation (inverse coef4)
///
void AliTPCPoissonSolver::VCycle3D(const Int_t rows, const Int_t columns, const Int_t phislices, const Int_t symmetry,const Int_t gridfrom,const  Int_t gridto,
		const Int_t NPRE, const Int_t NPOST,const Float_t gridSizeR, const Float_t ratioZ,
		std::vector<TMatrixD **> & tvArrayV,	std::vector<TMatrixD **> & tvCharge,std::vector<TMatrixD **> & tvResidu,
		std::vector<float> &coef1,std::vector<float> &coef2, std::vector<float> &coef3,std::vector<float> &coef4,std::vector<float> &icoef4) {
					
	Float_t h,h2,ih2,tempratioZ, tempratioPhi,radius, tempgridSizePhi;
	TMatrixD **curArrayV,**curArrayVC;
	TMatrixD **curCharge;
	TMatrixD **residu;
	Int_t iOne,jOne, kOne, trows, tcolumns, tphislices, otphislices,count, nnphis;
	
	
	curArrayV = NULL;
	curArrayVC = NULL;
	curCharge = NULL;
	residu = NULL;
	
	iOne = 1 << (gridfrom-1); 
	jOne = 1 << (gridfrom-1); 	
	kOne = 1 << (gridfrom-1); 
	
	nnphis = phislices;
	
	while (nnphis % 2 == 0) { 
		nnphis /= 2;
	}
	
						
	trows = iOne == 1 ? rows : rows/iOne + 1;
	tcolumns = jOne == 1 ? columns : columns/jOne + 1;	
	tphislices = kOne == 1 ? phislices : phislices/kOne;
	tphislices = tphislices < nnphis ? nnphis : tphislices;
		

	//AliInfo(Form("Grid information: trows=%d, tcols=%d, tphislices=%d\n", trows,tcolumns,tphislices));
	

	for (  count = gridfrom ; count <= gridto-1 ; count++ ) {		
		otphislices = tphislices;
		
		
		h   =  gridSizeR  * iOne ;
		h2 = h*h;
		ih2 = 1.0/h2;
		tempgridSizePhi =  TMath::TwoPi()/tphislices;  // phi now is multigrid
		
		tempratioPhi    =  h*h / (tempgridSizePhi*tempgridSizePhi) ;  // ratio_{phi} = gridsize_{r} / gridsize_{phi}
			
		tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;
		
		for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
			icoef4[i] = 1.0 / coef4[i];
		}
		
		curArrayV = tvArrayV[count - 1];
		curCharge = tvCharge[count - 1];		
		residu = tvResidu[count - 1];
	
		
		//AliInfo("Before Pre-smoothing");
		//curArrayV->Print();
		
		// 1) Pre-Smoothing: Gauss-Seidel Relaxation or Jacobi		
		for (Int_t jpre=1;jpre <= NPRE;jpre++) {				
			
			//AliInfo(Form("Call Relax3D: trows=%d, tcols=%d, tphislices=%d\n", trows,tcolumns,tphislices));
			
			Relax3D(curArrayV,curCharge,trows,tcolumns,tphislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
		} // end pre smoothing
		
		// 2) Residu calculation 
		//AliInfo(Form("Call Residu3D: trows=%d, tcols=%d, tphislices=%d\n", trows,tcolumns,tphislices));
		
		Residu3D(residu,curArrayV,curCharge,trows,tcolumns,tphislices,symmetry,ih2,tempratioZ,coef1,coef2,coef3,icoef4);			
				
		iOne = 2*iOne;
		jOne = 2*jOne;
		kOne = 2*kOne;
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		tphislices = phislices/kOne;
		tphislices = tphislices < nnphis ? nnphis : tphislices;
	
		
		curCharge =  tvCharge[count];		
		curArrayV = tvArrayV[count];
		//AliInfo(Form("Call Restrict3D: trows=%d, tcols=%d, tphislices=%d, otphislices\n", trows,tcolumns,tphislices,otphislices));
			
		//3) Restriction 
		Restrict3D(curCharge, residu, trows, tcolumns,  tphislices, otphislices);
		
		//4) Zeroing coarser V
		for(Int_t m=0;m<tphislices;m++) {
				curArrayV[m]->Zero();		
		}		

	}
	
	
	// coarsest grid
	h   =  gridSizeR  * iOne ;
	h2 = h*h;
	tempgridSizePhi =  TMath::TwoPi()/tphislices;  // phi now is multigrid
		
	tempratioPhi    =  h*h / (tempgridSizePhi*tempgridSizePhi) ;  // ratio_{phi} = gridsize_{r} / gridsize_{phi}
	tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;

	
	for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
	}
	
	// 3) Relax on the coarsest grid
	//AliInfo(Form("Call Relax3D: trows=%d, tcols=%d, tphislices=%d, otphislices\n", trows,tcolumns,tphislices,otphislices));
		
	Relax3D(curArrayV,curCharge,trows,tcolumns,tphislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
	
	
	// back to fine
	for (  count = gridto-1 ; count >= gridfrom ; count-- ) {
		otphislices = tphislices;
		
		iOne = iOne/2;
		jOne = jOne/2;
		kOne = kOne/2;
		
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		tphislices = kOne == 1 ? phislices : phislices/kOne;
		tphislices = tphislices < nnphis ? nnphis : tphislices;
	
		
		h   =  gridSizeR  * iOne ;
		h2 = h*h;
		tempgridSizePhi =  TMath::TwoPi()/tphislices;  // phi now is multigrid
		
		tempratioPhi    =  h*h / (tempgridSizePhi*tempgridSizePhi) ;  // ratio_{phi} = gridsize_{r} / gridsize_{phi}
			
		tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;
	
		curCharge =  tvCharge[count-1];				
		curArrayV = tvArrayV[count-1];
		curArrayVC = tvArrayV[count];
	
		// 4) Interpolation/Prolongation
		//AliInfo(Form("Call AddInterp3D: trows=%d, tcols=%d, tphislices=%d, otphislices\n", trows,tcolumns,tphislices,otphislices));
		
	
		AddInterp3D(curArrayV,curArrayVC,trows,tcolumns,tphislices,otphislices);
			

		
		for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
		}		
	
		
		// 5) Post-Smoothing: Gauss-Seidel Relaxation
		for (Int_t jpost=1;jpost <= NPOST;jpost++) {			
			//AliInfo(Form("Call Relax 3D: trows=%d, tcols=%d, tphislices=%d, otphislices\n", trows,tcolumns,tphislices,otphislices));
		
			Relax3D(curArrayV,curCharge,trows,tcolumns,tphislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
		} // end post smoothin
	}
}


///
/// Set matrix exact solution for relative error calculation
///
/// \param exactSolution TMatrixD** pointer to exact solution (potential) in 3D
/// \param fPhiSlices const Int_t number of phi slices
///
void AliTPCPoissonSolver::SetExactSolution(TMatrixD ** exactSolution,const Int_t fPhiSlices) {
	Double_t maxAbs;
	fExactSolution = exactSolution;
	fExactPresent = kTRUE;
	fMaxExact = 0.0;
	for (Int_t m=0;m<fPhiSlices;m++) {
		maxAbs = TMath::Max(TMath::Abs((*fExactSolution[m]).Max()), TMath::Abs((*fExactSolution[m]).Min()));
		if (maxAbs > fMaxExact)
			fMaxExact = maxAbs;
	}
}

///
/// Relative error calculation: comparison with exact solution
///
/// \param curArrayV TMatrixD** current potential (numerical solution)
/// \param tempArrayV TMatrixD** temporary matrix for calculating error
/// \param rows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param columns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param phislice const Int_t phi slices
///
Double_t AliTPCPoissonSolver::GetErrorExact(TMatrixD **curArrayV,TMatrixD **tempArrayV,  const Int_t phislices) {
	Double_t error = 0.0;
	
	if (fExactPresent == kTRUE) {
		for (Int_t m=0;m<phislices;m++) {
			(*tempArrayV[m]) = (*fExactSolution[m]) - (*curArrayV[m]);						
			(*tempArrayV[m]) *= 1.0 / GetMaxExact();
			if (tempArrayV[m]->E2Norm() > error) error = tempArrayV[m]->E2Norm();					
			//printf("%f\n",tempArrayV[m]->E2Norm();
		}	
	}
	return error;
}


///
/// Relative error calculation: comparison with exact solution
///
/// \param curArrayV TMatrixD** current potential (numerical solution)
/// \param tempArrayV TMatrixD** temporary matrix for calculating error
/// \param rows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param columns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param phislice const Int_t phi slices
///
Double_t AliTPCPoissonSolver::GetErrorConv(TMatrixD **curArrayV,TMatrixD **prevArrayV,  const Int_t phislices) {
	Double_t error = 0.0;
	
	for (Int_t m=0;m<phislices;m++) {
		
		// absolute 
		(*prevArrayV[m]) = (*prevArrayV[m]) - (*curArrayV[m]);						
		
		
		
		if (  prevArrayV[m]->E2Norm() > error) error = prevArrayV[m]->E2Norm() ;							
	}
	return error;
}


///////////////////// interface for GPU ///////////////////

/// VCycle 3D2D, V Cycle 3D in multigrid with constant phislices
/// fine-->coarsest-->fine, propagating the residu to correct initial guess of V
///
/// Algorithm:
///
///    NOTE: In order for this algorith to work, the number of rows and columns must be a power of 2 plus one.
///    The number of rows and COLUMNS can be different.
///
///    ROWS       ==  2**M + 1
///    COLUMNS    ==  2**N + 1
///    PHISLICES  ==  Arbitrary but greater than 3
///
///    DeltaPhi in Radians
/// \param rows Int_t number of grid in rows (in r-direction) for coarser grid should be 2^N + 1, finer grid in 2^{N+1} + 1
/// \param columns Int_t number of grid in columns (in z-direction) for coarser grid should be  2^M + 1, finer grid in 2^{M+1} + 1
/// \param gridfrom const Int_t finest level of grid
/// \param gridto const Int_t coarsest level of grid
/// \param NPRE const Int_t number of smoothing before coarsening
/// \param NPOST const Int_t number of smoothing after coarsening
/// \param gridSizeR const Float_t grid size in r direction (OPTION,  recalculate)
/// \param ratio const Float_t ratio between square of grid r and grid z (OPTION,  recalculate)
/// \param tvArrayV vector<TMatrixD *> vector of V potential in different grids
/// \param tvCharge vector<TMatrixD *> vector of charge distribution in different grids
/// \param tvResidu vector<TMatrixD *> vector of residu calculation in different grids
/// \param coef1 std::vector<float>& coefissien for relaxation (r direction)
/// \param coef2 std::vector<float>& coefissien for relaxation (r direction)
/// \param coef3 std::vector<float>& coefissien for relaxation (ratio r/z)
/// \param coef4 std::vector<float>& coefissien for relaxation (ratio for grid_r)
/// \param icoef4 std::vector<float>& coefissien for relaxation (inverse coef4)
///
void AliTPCPoissonSolver::VCycle3D2DGPU
(
	const Int_t rows, 
	const Int_t columns, 
	const Int_t phislices, 
	const Int_t symmetry,
	const Int_t gridfrom,
	const  Int_t gridto,
	const Int_t NPRE, 
	const Int_t NPOST,
	const Float_t gridSizeR, 
	const Float_t ratioZ, 
	const Float_t ratioPhi,
	std::vector<TMatrixD **> & tvArrayV,	
	std::vector<TMatrixD **> & tvCharge,
	std::vector<TMatrixD **> & tvResidu,
	std::vector<float> &coef1,std::vector<float> &coef2, 
	std::vector<float> &coef3,std::vector<float> &coef4,
	std::vector<float> &icoef4
) 
{
					
	Float_t h,h2,ih2,tempratioZ, tempratioPhi,radius;
	TMatrixD **curArrayV,**curArrayVC;
	TMatrixD **curCharge;
	TMatrixD **residu;
	Int_t iOne,jOne, trows, tcolumns, count;
	
	curArrayV = NULL;
	curArrayVC = NULL;
	curCharge = NULL;
	residu = NULL;
	
	
	iOne = 1 << (gridfrom-1); 
	jOne = 1 << (gridfrom-1); 	
						
	trows = iOne == 1 ? rows : rows/iOne + 1;
	tcolumns = jOne == 1 ? columns : columns/jOne + 1;

	
	
	for (  count = gridfrom ; count <= gridto-1 ; count++ ) {		
		h   =  gridSizeR  * iOne ;
		h2 = h*h;
		ih2 = 1.0/h2;
			
		tempratioPhi    =  ratioPhi * iOne * iOne ; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
		tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;


		for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
			icoef4[i] = 1.0 / coef4[i];
		}
		
		curArrayV = tvArrayV[count - 1];
		curCharge = tvCharge[count - 1];		
		residu = tvResidu[count - 1];
	
		
		//AliInfo("Before Pre-smoothing");
		//curArrayV->Print();
		
		// 1) Pre-Smoothing: Gauss-Seidel Relaxation or Jacobi		
		for (Int_t jpre=1;jpre <= NPRE;jpre++) {				
			Relax3D(curArrayV,curCharge,trows,tcolumns,phislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
		} // end pre smoothing
		
		// 2) Residu calculation 
		Residu3D(residu,curArrayV,curCharge,trows,tcolumns,phislices,symmetry,ih2,tempratioZ,coef1,coef2,coef3,icoef4);			
				
		iOne = 2*iOne;
		jOne = 2*jOne;
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		
		curCharge =  tvCharge[count];		
		curArrayV = tvArrayV[count];
		
		//3) Restriction 
		//Restrict2D(*curCharge,*residu,trows,tcolumns);			
		Restrict3D(curCharge, residu, trows, tcolumns,  phislices, phislices);
		
		//4) Zeroing coarser V
		for(Int_t m=0;m<phislices;m++) {
				curArrayV[m]->Zero();		
		}		
	}
	
	// coarsest grid
	h   =  gridSizeR  * iOne ;
	h2 = h*h;
			
	tempratioPhi    =  ratioPhi * iOne * iOne ; 
	tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;

	
	for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
	}
	
	// 3) Relax on the coarsest grid
	Relax3D(curArrayV,curCharge,trows,tcolumns,phislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
	
	// back to fine
	for (  count = gridto-1 ; count >= gridfrom ; count-- ) {
		iOne = iOne/2;
		jOne = jOne/2;
		
		trows = iOne == 1 ? rows : rows/iOne + 1;
		tcolumns = jOne == 1 ? columns : columns/jOne + 1;
		
		h   =  gridSizeR  * iOne ;
		h2 = h*h;
			
		tempratioPhi    =  ratioPhi * iOne * iOne ; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
		tempratioZ      =  ratioZ   * iOne * iOne / ( jOne * jOne ) ;
	
		curCharge =  tvCharge[count-1];				
		curArrayV = tvArrayV[count-1];
		curArrayVC = tvArrayV[count];
	
		// 4) Interpolation/Prolongation
		AddInterp3D(curArrayV,curArrayVC,trows,tcolumns,phislices,phislices);
			

		
		for ( Int_t i = 1 ; i < trows-1 ; i++ )  {
			radius = AliTPCPoissonSolver::fgkIFCRadius + i*h;
			coef1[i] = 1.0 + h/(2*radius);
			coef2[i] = 1.0 - h/(2*radius);
			coef3[i] = tempratioPhi/(radius*radius);
			coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
		}		
	
		
		// 5) Post-Smoothing: Gauss-Seidel Relaxation
		for (Int_t jpost=1;jpost <= NPOST;jpost++) {			
			Relax3D(curArrayV,curCharge,trows,tcolumns,phislices,symmetry, h2, tempratioZ,coef1,coef2,coef3,coef4);			
		} // end post smoothing	
	}			
}

