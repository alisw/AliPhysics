#error This is not for compilation 
/** 
 * @page mult_doc The basic multiplicity code 
 *
 * Module: @ref pwglf_forward_aod
 *
 * @tableofcontents 
 *
 * @section mult_intro Introduction 
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
 * @section mult_structure Structure 
 *
 * There are classes for containing data, classes that represent
 * analysis tasks, and classes that perform calculations, as well as
 * specialized classes for analysis of simulation (MC) output. 
 *
 * @subsection mult_struct_data Data structures 
 *
 * The classes AliAODForwardMult and AliAODCentralMult each contain a
 * 2-dimensional @f$(\eta,\phi)@f$ map of the charged particle
 * multiplicity _per_ _event_.  The data is fully corrected for
 * acceptance and secondary particle production, though this can be
 * turned off if needed. 
 *
 * @subsection mult_struct_tasks Tasks 
 * 
 * For the forward analysis, there are two tasks:
 * AliForwardMultiplictyTask and AliForwardMCMultiplicityTask - both
 * of which derive from the base class AliForwardMultiplicityBase.
 * AliForwardMultiplictyTask is intented for analysis of collision
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
 * AliForwardMultiplicityBase and it's derivatives. 
 *
 * The module also defines the utility task AliCopyHeaderTask, which
 * copies the ESD header information into the AOD header.  This was
 * defined to address a short coming in the overall ANALYSIS
 * framework.  This task can be added by the ::AddTaskCopyHeader
 * script function. 
 * 
 * @subsection mult_workers Workers 
 *
 * @subsection mult_mc_workers MC Workers 
 *
 * @subsection mult_misc Misc.
 *
 */
//
// EOF
//
