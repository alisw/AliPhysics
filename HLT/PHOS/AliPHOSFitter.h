#ifndef ALIPHOSFITTER_H
#define ALIPHOSFITTER_H
/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



/**
 * The PhosFitter class is the class for extracting the basic signal parameters
 * "timing" and "energy" from the PHOS raw data. Physical data will for a given readout channel be
 * a sequense of ADC digitized 10 bit integer values, however for performance reasons all values used in
 * calculation is of type double.
 **/
class AliPHOSFitter
      {
      public:
	/**
	 * Default constructor
	 **/
	AliPHOSFitter();

	/**
	 * Main constructor
	 * @param dataPtr Data array for wich a subarray will be taken to perform the fit
	 * @param fs the sampling frequency in entities of MHz. Needed in order to calculate physical time
	 **/
	AliPHOSFitter(double *dataPtr, double fs);

	/**
	 * Destructor
	 **/
	~AliPHOSFitter();


	/**
	 * Attemps to level the basline to zero.
	 * The baseline will be calculated from the pretrigger samples and subtracted from the 
	 * data array. 
	 * If pretrigger samples are not present then the basline correction will be incorrect. 
	 * @param dataPtr array for wich to correct the basline 
	 * @param N the number of pretrigger samples used to calculate the baseline.
	 **/
	void BaselineCorrection(double *dataPtr, int N);

	/**
	 * Shifts the basline with the amount given by baselineValue
	 * If pretrigger samples are not present then the basline correction will be incorrect. 
	 * @param dataPtr array for wich to correct the basline 
	 * @param BaslineValue the basline value to subtract..
	 **/
	void BaselineCorrection(double *dataPtr, double baselineValue);

	/**
	 * Extraction of timing and energy using a least square fit assuming uncorrelated measurement errors.
	 * This is also called a "Chi square fit". The default method is the Levenberg Marquardt nonlinear fit. 
	 * If the errors are correlated (which is typically the case)
	 * the timing and energy will not be mimimum variance lower bound. Correlated errors might also give a
	 * systematic error. The parameters "start" and "length" defines the subarray of the data array set in the
	 * constructor that will be used in the fit. "start + length" cannot exeed the length of the data array.
	 * The baseline must be subtracted before performing the fit othervise the result can be biased.
	 * the initial guess parameters is found by the method "MakeInitialGuess".	
	 * @param start the start index of the subarray of the data array. 
	 * @param length the number of samples to use starting from index 
 	 * @param start the start index of the subarray of the data array. 
	 * @param length the number of samples to use starting from index 
	 **/
	void FitChiSquare(int start, int lenght);

	/** 
	 * Extraction of timing and energy using a least square fit assuming uncorrelated measurement errors.
	 * This is also called a "Chi square fit". The default method is the Levenberg Marquardt nonlinear fit. 
	 * If the errors are correlated (which is typically the case)
	 * the timing and energy will not be mimimum variance lower bound. Correlated errors might also give a
	 * systematic error. The parameters "start" and "length" defines the subarray of the data array set in the
	 * constructor that will be used in the fit. "start + length" cannot exeed the length of the data array.
	 * The baseline must be subtracted before performing the fit othervise the result can be biased.
	 * A good initial guess greatly enchanes performance. A poor initial guess might give a solution that is a local 
	 * minima instead of a global minima. If the startindex of the pulse is the first sample above a low (2-3 ADC levels)
	 * above the baseline, then the initial guess of t0 an be set to zero. A good initial guess for amplitude will
	 * for a gamma 2 function be for instance (Max sample valu - baseline)*exp(2). If a Gamma N function is used
	 * multiply with exp(N).
 	 * @param start the start index of the subarray of the data array. 
	 * @param length the number of samples to use starting from index 
	 * @param tGuess initial guess for  timing (in entities of samples)
	 * @param aGuess initial guess for energy in entities of (ADC channels)*exp(2) (for a gamma 2 fuction)
	 **/
	void FitChiSquare(int start, int lenght, double tGuess, double aGuess);

	/**
	 * Extraction of timing an energy using the
	 * K-level method.
	 * @param kLevel the K-level
	 * "start" and "length" defines the subarray of the data array (set in the constructor)
	 * @param start index of samples to use relative to the the data array
	 * @param length the number of samples to use (starting at "start")
	 **/
	void FitKLevel(double kLevel, int start, int lenght);
	
	/** 
	 * This method gives the statistically best possible unbiased estimate of t0 and amplitude provided that a good
	 * estimate of the covariance matrix exists.
	 * Extraction of timing and energy using a least square fit assuming correlated errors.
	 * The user must supply the autocovariance matrix. For N number of samples this matrix
	 * will be a NxN matrix. 
	 * The default method is the Levenberg Marquardt nonlinear fit. 
	 * The parameters "start" and "length" defines the subarray of the data array set in the
	 * constructor that will be used in the fit. "start + length" cannot exeed the length of the data array.
	 * The baseline must be subtracted before performing the fit othervise the result wil biased.
	 * If the startindex of the pulse is the first sample above a low (2-3 ADC levels)
	 * above the baseline, then the initial guess of t0 an be set to zero. A good initial guess for amplitude will
	 * for a gamma 2 function be for instance (Max sample valu - baseline)*exp(2). If a Gamma N function is used
	 * multiply with exp(N).
	 * The correlation matrix for the parameters is written to the matrix pointe to by pCovar 
	 * @param start the start index of the subarray of the data array. 
	 * @param length the number of samples to use starting from "start" 
	 * @param kmCovarPtrPtr the measurement correlation matrix 	
	 * @param pCovarPtrPtr the correlation matrix of the estimated parameters
	 * @param tGuess initial guess of t0 in entities of sampling intervals
	 * @param aGuess initial guess of the amplitude in entities of ADC levels*exp(2) 
	**/
	void FitLeastMeanSquare(int start, int lenght, const double **kmCovarPtrPtr, double **pCovarPtrPtr);

	/** 
	 * This method gives the statistically best possible unbiased estimate of t0 and amplitude provided that a good
	 * estimate of the covariance matrix exists.
	 * For N number of samples this covariance matrix must be an NxN matrix.
	 * The default method is the Levenberg Marquardt nonlinear fit. 
	 * The parameters "start" and "length" defines the subarray of the data array set in the
	 * constructor that will be used in the fit. "start + length" cannot exeed the length of the data array.
	 * The baseline must be subtracted before performing the fit othervise the result wil biased.
	 * A good initial guess will greatly enchance performance.
	 * a poor initial guess will slow down the algorithm and possibly lead to convergens to local minima
	 * of the cost function instead of the global minima. A very bad intial guess might give divergence of the solution.
	 * If the startindex of the pulse is the first sample above a low (2-3 ADC levels)
	 * above the baseline, then the initial guess of t0 can be set to zero. A good initial guess for amplitude will
	 * for a gamma 2 function be for instance (Max sample value - baseline)*exp(2). If a Gamma N function is used
	 * multiply with exp(N).
	 * @param start the start index of the subarray of the data array. 
	 * @param length the number of samples to use starting from index 
	 * @param mCovar the measurement correlation matrix 	
	 * @param pCovar the correlation matrix of the estimated parameters
	 * @param tGuess initial guess of t0 in entities of sampling intervals
	 * @param aGuess initial guess of the amplitude in entities of ADC levels*exp(2)
	**/
	void FitLeastMeanSquare(int start, int lenght, const double **kmCovarPtrPtr, double **pCovarPtrPtr, double tGuess, double aGuess);

	/**
	 * Extraction of timing and energy using the Peakfinde Algorithm.
	 * The. The parameters "start" and "length" defines a sub array  of the data array
	 * that will be used for the the fit. If start+length must not exeed the total length
	 * of the Data array. "start" must be chosen as close as possible to t0.
	 * The baseline must also be subtracted.
	 * The length of "tVector" and "aVector" mus be equal to length.
	 * "index + length" must not exeed the length of the data array set in the constructor.
	 * @param start the start index of the subarray of the data array. 
	 * @param length the number of samples to use starting from index 
	 * @param tVector the peakfinder vector for timing
	 * @param aVector the peakfinder vector for amplitude (energy)
	 **/
	void FitPeakFinder(int start, int lenght, double *tVector, double *aVector);

	/**
	 * This method finds the start index of the pulse (relative to the data array) by serching for
	 * three or more continious samples above trheshold.
	 * @param treshold treshold to use when searchin for the start of the pulse
	**/
	int FindStartIndex(double treshold);
	
	/**
	 * Gives the timing in entities of sample indexes
	 * Physical time is found by multiplying  with the sampling intervall (Ts).
	 **/	
	float GetTiming();

	/**
	 * Gives the time in entities of ADC channels (quantization levels).  
	 * Absolute enrgy is found by multiplying with offline calibration constants.
	**/
	float GetEnergy();

	//	/**
	//	 * Set data array. Overrrides the data array set in the constructor.
	//	 **/
	//	void SetData(int *data);

	/**
	 * Set data array. Overrides data data array set in the constructor.
	 **/
	void SetData(double *data);


      private:

	/**
	 * This function applies only to the Chi and Least mean square fit. An initial guess is made
	 * based on the average of the first 5 samples and the first value exeeding this value.
	 **/
	void MakeInitialGuess();

	/**
	 * This function applies only to the Chi and Least mean square fit. An initial guess is made
	 * based on the average of the first 5 samples and the first value exeeding threshold + this value.
	 * @param treshold The index of the first value above treshold is ntaken to be the first value.
	 **/
	void MakeInitialGuess(int treshold);
	
	double    *fFloatDataPtr;    /**<Float representation of data that should be fitted */
	double     fSampleFrequency; /**<The ADC sample frequency in MHz used under data taking */
	double     fTau;	     /**<The risetime in micro seconds*/		 
	double     fDTof;            /**<Time of flight in entities of sample intervals */
	double     fDAmpl;           /**<Amplitude in entities of ADC levels*/
	double     fDTofGuess;       /**<Initial guess for t0*/
	double     fDAmplGuess;      /**<Initial guess for amplitude*/
	double   **kfMCovarPtrPtr;   /**<Covariance matrix of the measurements*/
	double   **fPCovarPtrPtr;    /**<Covariance matrix of the estimated parameters*/
      };

#endif
