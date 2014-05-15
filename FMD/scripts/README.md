Content of this directory
=========================

*   `ApplyAlignment.C`
       Read in the geometry, and get alignment data from CDB, and  apply
       that to the geometry. 
*   `CheckAlign.C`
       Check if alignment is properly applied in simulations 
*   `CheckCalibData.C`
       Check if a file contains the right calibration data
*   `CheckOverlaps.C`
       Run overlap checker on a loaded geometry 
*   `CheckRaw.C`
       Script that contains a class to compare the raw data written to the
       digits it's created from.
*   `checkSizes.sh`
       Obsolete shell script
*   `Compare.C`
       Compare energy loss deposition from Fluka and GEANT3.21
*   `CompareMatDep.C`
       (obsolete) Compare the secondary correction from variations in
       material density 
*   `CompareTrackRefsAndHits.C`
       Compare number of track references and hits from simulation
*   `Compile.C`
       Compile other scripts using AcLIC
*   `Convert2Raw.C`
       Convert digits to raw data using AliFMDRawWriter 
*   `ConvertGeom.C`
       Write out all medium parameters in geometry 
*   `deadChannels.C`
       (obsolete) Extract the dead channel information from OCDB files 
*   `Digits2Raw.C`
       (obsolete) Convert digits to raw data using AliFMD
*   `DisplayDigits.C`
       Show digits using AliFMDDisplay 
*   `DisplayESD.C`
       Show ESD info using AliFMDDisplay 
*   `DisplayHits.C`
       Show MC Hit info using AliFMDDisplay 
*   `DisplayRecs.C`
       Show RecPoint info using AliFMDDisplay 
*   `Document.C`
       (obsolete) Use THtml to document the FMD code
*   `DoDrawCalib.C`
       Script that calls DrawCalib 
*   `DrawBothDigits.C`
       Draw Digits versus SDigits using AliFMDInput 
*   `DrawCalib.C`
       GUI to draw calibrations 
*   `DrawCalibRaw.C`
       Draw ELoss dist from calibrated raw data using AliFMDInput
*   `DrawDigits.C`
       Draw the ELoss from digits using AliFMDInput 
*   `DrawDigitsRecs.C`
       Draw the ELoss from digits versus multiplicity from rec points
       using AliFMDInput 
*   `DrawESD.C`
       Draw the ELoss dist from ESDs using AliFMDInput - somewhat obsolete
*   `DrawGeometry.C`
       Draw the geometry 
*   `DrawHits.C`
       Draw the ELoss dist versus beta=p/m from MC hits
*   `DrawHitsDigits.C`
       Draw the MC ELoss versus ADCs 
*   `DrawHitsRecs.C`
       Draw the MC ELoss versus "multiplicity" 
*   `DrawHitsSDigits.C`
       Draw the MC ELoss versus sum-able ADCs 
*   `DrawLego.C`
       Draw summary LEGO (absorbtion, radiation length, ...) plots 
*   `DrawSDigits.C`
       Draw the sum-able ADCs 
*   `DrawTrackRefs.C`
       Draw number of track references in each ring 
*   `DrawXsection.C`
       Draw the x-section for ELoss from GEANT3.21 X-section run 
*   `Dummy.C`
       Make dummy detector classes 
*   `DummyConfig.C`
       An MC configuration that does no stepping in anything but FMD
*   `dump10bit.C`
       Dump 10bit words from a DDL file 
*   `ESDDataSize.C`
       Estimate the size of the FMD ESD object 
*   `esdQA.C`
       (obsolete) Do an ESD QA run 
*   `FancyDigits.C`
       Fancy display of digits using AliFMDFancy
*   `FancyHits.C`
       Fancy display of MC hits using AliFMDFancy
*   `FindCommonModeNoise.C`
       Script to estimate the common mode noise
*   `FixOCDBEntries.C`
       Script to correct OCDB objects according to file name 
*   `FullMapping.C`
       Map from Hardware address to XYZ coordinates 
*   `GeoGeometry.C`
       Script I used for rapid prototyping of the FMD3 geometry - in
       particular the support cone 
*   `GetMedia.C`
       Script that contains a class to get the media where a certain track
       was created.   This is used for background studies. 
*   `GetXsection.C`
       Script to get the various cross sections, energy loss, and ranges
       of a particular particle type in a particular medium.  See also
       DrawXsecttion.C 
*   `Hits2Digits.C`
       (obsolete) Convert MC hits to digits using AliFMD
*   `Hits2SDigits.C`
       (obsolete) Convert MC hits to sum-able digits using AliFMD
*   `LoadDummy.C`
       Load the Dummy detector class 
*   `MakeAlignment.C`
       Make the default alignment 
*   `MakeCalibration.C`
       Make dummy calibrations 
*   `MakeFakeDigits.C`
       Make some fake digits 
*   `MakeFakeHits.C`
       Make some fake hits 
*   `MakeLego.C`
       Create LEGO plots from geometry 
*   `makelego.sh`
       Steer running MakeLego.C and DrawLego.C 
*   `MakeMap.C`
       Script to make a class derived from AliFMDMap. 
*   `MakeRecoParam.C`
       Make default reconstruction parameters 
*   `MakeResidualAlignment.C`
       Make fake residual alignment data.
*   `MakeSpectra.C`
       Compile ELoss spectra in output file 
*   `MakeXsection.C`
       Steer generation of X-section data (see GetXsection and DrawXsection)
*   `MediaTable.C`
       Script to make a file with a list of medium numbers for all active
       detectors.
*   `NodeGeometry.C`
       Script I used for rapid prototyping of the FMD3 geometry - in
       particular the support cone 
*   `PatternDigits.C`
       Display the pattern of digits 
*   `PatternESD.C`
       Display the pattern of ESD information 
*   `PatternHits.C`
       Display the pattern of MC hits
*   `PatternRaw.C`
       Display the pattern of Raw digits
*   `PatternRecs.C`
       Display the pattern of RecPoints
*   `PatternSDigits.C`
       Display the pattern of MC sum-able digits
*   `pdc06_config.C`
       MC config 
*   `pdc06_rec.C`
       MC reconstruction 
*   `pdc06_sim.C`
       MC simulation
*   `Poisson.C`
       Script that illustrates the poisson Nch estimation from ESD 
*   `PoissonHit.C`
       Check Poisson estimate of Nch versus hits 
*   `polar.C`
       A test 
*   `PrintAlignment.C`
       Print alignment data
*   `PrintCalibration.C`
       Print OCDB calibration data
*   `PrintSensorVertices.C`
       Print senser vertices from geometry 
*   `Raw2ESD.C`
       Convert raw data to ESD using the AliFMDReconstructor 
*   `RawTest.C`
       (obsolete) Small script to test consistency of writing and
       reading raw data. 
*   `ReadRaw.C`
       Read raw data using AliRawReader and AliFMDRawReader 
*   `Resolution.C`
       Plot resolution of sensors 
*   `RunAnaESD.C`
       (obsolete) Run analysis over ESD 
*   `RunAnaKine.C`
       (obsolete) Run analysis over MC data
*   `RunQATest.C`
       Run a test of QA framework
*   `RunSimpleChain.C`
       A script that will run the most simple chain possible.  Uses:
        - MakeFakeHits.C
        - Hits2Digits.C
        - Hits2SDigits.C
        - Digits2Raw.C
        - Raw2ESD.C 
*   `ShowCoordinates.C`
      A script to dump the physical coordinates as given by the
      geometry.
*   `ShowFMDITS.C`
      Draw FMD and ITS geometry 
*   `SpectraMonitor.C`
      Monitor data 
*   `TestAcc.C`
      Acceptance correction from strip length 
*   `TestAlignable.C`
      Test align-able volumes 
*   `TestAltroMapping.C`
      Check integrety of Hardware2Detector and Detector2Hardware
*   `TestESD.C`
      Test Read/Write of AliESDFMD object 
*   `TestESDCompat.C`
      Same as above 
*   `TestESDPhi.C`
      Test the geometry from ESD and from full geometry 
*   `TestFloatMap.C`
      Test I/O of floating point FMD map 
*   `TestGainDA.C`
      Run the gain DA
*   `TestHWMap.C`
      Test hardware address map by converting from detector coordinates
      to hardware address and then back again.
*   `TestIndex.C`
      Test of AliFMDIndex and AliFMDObjIndex 
*   `TestMapAccess.C`
      Test of ForAll access to AliFMDMap 
*   `TestMapAlgebra.C`
      Test algebra on map 
*   `TestMap.C`
      Test I/O of ALiFMDMap
*   `TestMapIO.C`
      Test I/O of ALiFMDMap
*   `TestPedestalDA.C`
      Run the pedestal DA
*   `TestPreprocessor.C`
      Test the SHUTTLE pre-processor using fake data
*   `TestRaw2SDigits.C`
      Test creation of sum-able digits using AliSimulation 
*   `TestRawIO.C`
      Test I/O of raw data using AliFMDRawWriter and AliFMDRawReader 
*   `TestRawReader.C`
      Check iterative usage of AliFMDRawReader 
*   `TestShaping.C`
      Create an example shaping curve 
*   `TestSurveyToAlignObjs.C`
      Convert survey data into alignment objects 
*   `TestZeroSuppress.C`
      Test the zero-suppression filter of AliFMDRawWriter 
*   `VA1Response.C`
      Script to try to fit the reponse function of the VA1 signals, based
      on a finite number of ALTRO samples. 
*   `VA1Train.C`
      Small script that shows a signal train from a VA1 pre-amp.
*   `Wafer.C`
      Small script that I used to make some intial testing of the wafer
      layout and geometry. 
*   `WriteMedArrays.C`
      Write medium map 
*   `XSection.C`
       Script to get the various cross sections, energy loss, and ranges
       of a particular particle type in a particular medium. 
*   `dqm`
    *   `fmd_online.C`
	      Event display using EVE 
    *   `PatternCalib.C`
	      Event display using AliFMDPattern 
    *   `PatternRaw.C`
	      Event display using AliFMDPattern     
    *   `SpectraRaw.C`
	      Event display using AliFMDPattern 
