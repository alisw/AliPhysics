#error This is not for compilation 
/** 
 * @page pwglf_fwd_mult_doc The basic multiplicity code 
 *
 * Module: @ref pwglf_forward_aod
 *
 * @section pwglf_fwd_mult_intro Introduction 
 *
 * The code in this section defines methods to generate a fully
 * corrected 2-dimensional histogram of the charged particle
 * multiplicity in the forward regions _per_ _event_.  This histogram
 * is then stored in an object of class AliAODForwardMult on a
 * seperate branch on the AOD tree. This object can then be used for
 * sub-sequent detail analysis, like @f$dN_{ch}/d\eta@f$,
 * @f$v_n{m}@f$, @f$P(N_{ch})@f$, and @f$Corr(n_f,n_b)@f$. 
 *
 * The input to this analysis is the ESD information from the FMD,
 * plus some corrections stored in AliROOT (should be migrated to
 * OADB).  
 * 
 * Also defined in this module, is code to produce a similar object
 * AliAODCentralMult, but this time the input information comes from
 * the SPD clusters on the first layer.  The two analysis are similar
 * in methodology. 
 * 
 * @section pwglf_fwd_mult_structure Structure 
 *
 * There are classes for containing data, classes that represent
 * analysis tasks, and classes that perform calculations, as well as
 * specialized classes for analysis of simulation (MC) output. 
 *
 * @subsection pwglf_fwd_mult_struct_data Data structures 
 *
 * The classes AliAODForwardMult and AliAODCentralMult each contain a
 * 2-dimensional @f$(\eta,\phi)@f$ map of the charged particle
 * multiplicity _per_ _event_.  The data is fully corrected for
 * acceptance and secondary particle production, though this can be
 * turned off if needed. 
 *
 * @subsection pwglf_fwd_mult_struct_tasks Tasks 
 * 
 * For the forward analysis, there are two tasks:
 * AliForwardMultiplicityTask and AliForwardMCMultiplicityTask - both
 * of which derive from the base class AliForwardMultiplicityBase.
 * AliForwardMultiplicityTask is intented for analysis of collision
 * data, while AliForwardMCMultiplicityTask is intented simulation
 * data, and the ::AddTaskForwardMult script function automatically
 * selects the appropriate one.  Both of these tasks uses worker
 * classes to do the actual computations. 
 *
 * For the analysis of the SPD cluster, there are also two tasks:
 * AliCentralMultiplicityTask and AliCentralMCMultiplicity task, where
 * the latter derives from the first.  Again, the ::AddTaskCentralMult
 * script function automatically choses the appropriate implementation
 * based on the set-up of the train. 
 *
 * The special task AliForwardMCCorrectionsTask is used to generate
 * the various simulation based corrections needed by
 * AliForwardMultiplicityBase and it's derivatives, and the task
 * AliFMDEnergyFitterTask generates corrections for both real and
 * simulated input used by AliForwardMultiplicityTask and
 * AliForwardMCMultiplicityTask respectively.
 *
 * The module also defines the utility task AliCopyHeaderTask, which
 * copies the ESD header information into the AOD header.  This was
 * defined to address a short coming in the overall ANALYSIS
 * framework.  This task can be added by the ::AddTaskCopyHeader
 * script function. 
 * 
 * @subsection pwglf_fwd_mult_workers Workers 
 *
 * <dl>
 *  <dt>AliFMDEventInspector</dt>
 *  <dd>A worker that inspects the global event properties, such as
 *  location of the interaction point, the centrality (if applicable),
 *  and the fired triggers (both on-line and off-line).  The result of
 *  the inspection is returned to the caller task and diagnostics
 *  histograms are filled.</dd>
 *  <dt>AliFMDSharingFilter</dt>
 *  <dd>A worker that merges signals in the FMD silicon strips to
 *  correct for the situation where a charged particles energy loss is
 *  shared over 2 or 3 strips as a consquence of it's incident angle.
 *  The worker returns a new AliESDFMD object with the corrected
 *  signals and fills diagnostic histograms.</dd>
 *  <dt>AliFMDDensityCalculator</dt>
 *  <dd>This worker calculates the number of charged particles per
 *  @f$(\eta,\varphi)@f$ bin using one of two methods: Energy loss
 *  deconvolution or Poissonian statistics.  The user selects which
 *  method to use at start-up but both methods are fired and compared
 *  in the diagnostics histograms.  The output is 5 histograms - one
 *  for each sub-ring - of the event charged particle multiplicity in
 *  @f$(\eta,\varphi)@f$ bins.</dd>
 *  <dt>AliFMDCorrector</dt>
 *  <dd>Applies various corrections, such as @f$\phi@f$ acceptance,
 *  correction for secondary particle production, merging efficiency,
 *  and so on.  The exact corrections applied is selected by the
 *  user. The output is 5 histograms of the event charged particle
 *  multiplicity in @f$(\eta,\varphi)@f$ bins corrected for the
 *  selected effects.</dd>
 *  <dt>AliFMDHistCollector</dt>
 *  <dd>Collects the 5 histograms of the event charged particle
 *  multiplicity in @f$(\eta,\varphi)@f$ bins into one `super'
 *  histogram of the event charged particle multiplicity, taking care
 *  of overlaps and differences in @f$\varphi@f$ acceptance of the 5
 *  sub-rings.</dd>
 *  <dt>AliFMDEventPlaneFinder</dt> <dd>Calculates the reaction plane
 *  of the event from the FMD data</dd>
 * </dl>
 *
 * @subsection pwglf_fwd_mult_mc_workers MC Workers 
 *
 * @subsection pwglf_fwd_mult_misc Misc.
 *
 */
//
// EOF
//
