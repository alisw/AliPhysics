#error This is not for compilation 
/** 
 * @page pwglf_fwd_density_doc Charge Particle Multiplicity Densities 
 *
 * Module: @ref pwglf_forward_dndeta
 *
 * @tableofcontents 
 *
 * @section pwglf_fwd_density_intro Introduction 
 *
 * This code uses the AOD produced by the @ref pwglf_fwd_mult_doc code to
 * produce results on @f$ 1/N_{ev} dN_{ch}/d\eta@f$ in pp, PbPb, and
 * pPb collisions.
 * 
 * @section pwglf_fwd_density_tasks Tasks 
 *
 * <dl>
 * <dt>AliBasedNdetaTask</dt>
 * <dd>Base class for all other @f$ 1/N_{ev} dN_{ch}/d\eta@f$ tasks.
 * Provides a number of services and common calculations. </dd>
 * <dt>AliCentraldNdetaTask</dt>
 * <dd>Calculates the @f$ 1/N_{ev} dN_{ch}/d\eta@f$ in the central
 * region from the AliAODCentralMult objects stored in the AOD input.
 * </dd>
 * <dt>AliForwarddNdetaTask</dt>
 * <dd>Calculates the @f$ 1/N_{ev} dN_{ch}/d\eta@f$ in the forward
 * regions from the AliAODForwardMult objects stored in the AOD input.
 * </dd>
 * <dt>AliMCTruthdNdetaTask</dt>
 * <dd>Calculates the @f$ 1/N_{ev} dN_{ch}/d\eta@f$ from the MC `truth'
 * TH2D objects stored in the AOD input.
 * </dl>
 *
 * @section pwglf_fwd_density_scripts Scripts 
 *
 * <dl>
 * <dt>DrawdNdeta.C</dt>
 * <dd>Script to draw the results</dd>
 * </dl>
 */
//
// EOF
//

