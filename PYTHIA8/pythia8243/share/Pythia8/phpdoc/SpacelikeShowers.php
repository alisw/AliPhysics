<html>
<head>
<title>Spacelike Showers</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='SpacelikeShowers.php'>
 
<h2>Spacelike Showers</h2> 
<ol id="toc">
  <li><a href="#section0">Main variables</a></li>
  <li><a href="#section1">Dipole showers</a></li>
  <li><a href="#section2">Weak showers</a></li>
  <li><a href="#section3">Further variables</a></li>
  <li><a href="#section4">Technical notes</a></li>
</ol>

 
The PYTHIA algorithm for spacelike initial-state showers is 
based on the article [<a href="Bibliography.php#refSjo05" target="page">Sjo05</a>], where a 
transverse-momentum-ordered backwards evolution scheme is introduced, 
with the extension to fully interleaved evolution covered in 
[<a href="Bibliography.php#refCor10a" target="page">Cor10a</a>]. 
This algorithm is a further development of the virtuality-ordered one 
presented in [<a href="Bibliography.php#refSjo85" target="page">Sjo85</a>], with matching to first-order matrix 
element for <i>Z^0</i>, <i>W^+-</i> and Higgs (in the 
<i>m_t &rarr; infinity</i> limit) production as introduced in 
[<a href="Bibliography.php#refMiu99" target="page">Miu99</a>]. 
 
<p/> 
The normal user is not expected to call <code>SpaceShower</code> 
directly, but only have it called from <code>Pythia</code>, 
via <code>PartonLevel</code>. Nonetheless, some of the parameters below, 
in particular <code>SpaceShower:alphaSvalue</code>, 
would be of interest for uncertainty estimates and tuning exercises. 
Note that 
PYTHIA also incorporates an 
<?php $filepath = $_GET["filepath"];
echo "<a href='Variations.php?filepath=".$filepath."' target='page'>";?>automated framework</a> for shower 
uncertainty variations. 
 
<a name="section0"></a> 
<h3>Main variables</h3> 
 
The maximum <i>pT</i> to be allowed in the shower evolution is 
related to the nature of the hard process itself. It involves a 
delicate balance between not double-counting and not leaving any 
gaps in the coverage. The best procedure may depend on information 
only the user has: how the events were generated and mixed (e.g. with 
Les Houches Accord external input), and how they are intended to be 
used. Therefore a few options are available, with a sensible default 
behaviour. 
 
<br/><br/><table><tr><td><strong>SpaceShower:pTmaxMatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Way in which the maximum shower evolution scale is set to match the 
scale of the hard process itself. 
<br/>
<input type="radio" name="1" value="0" checked="checked"><strong>0 </strong>: <b>(i)</b> if the final state of the hard process  (not counting subsequent resonance decays) contains at least one quark  (<ei>u, d, s, c ,b</ei>), gluon or photon then <ei>pT_max</ei>  is chosen to be the factorization scale for internal processes  and the <code>scale</code> value for Les Houches input;  <b>(ii)</b> if not, emissions are allowed to go all the way up to  the kinematical limit.  The reasoning is that in the former set of processes the ISR  emission of yet another quark, gluon or photon could lead to  double-counting, while no such danger exists in the latter case.  <br/>
<input type="radio" name="1" value="1"><strong>1 </strong>: always use the factorization scale for an internal  process and the <code>scale</code> value for Les Houches input,  i.e. the lower value. This should avoid double-counting, but  may leave out some emissions that ought to have been simulated.  (Also known as wimpy showers.)  <br/>
<input type="radio" name="1" value="2"><strong>2 </strong>: always allow emissions up to the kinematical limit.  This will simulate all possible event topologies, but may lead to  double-counting.  (Also known as power showers.)  <br/>
<br/><b>Note 1:</b> Some processes contain matrix-element matching 
to the first emission; this is the case notably for single 
<ei>gamma^*/Z^0, W^+-</ei> and <ei>H^0</ei> production. Then default 
and option 2 give the correct result, while option 1 should never 
be used. 
<br/><b>Note 2:</b> as enumerated in the text, these options take effect 
both for internal and external processes. Whether a particular option 
makes sense depends on the context. For instance, if events for the same 
basic process to different orders are to be matched, then option 1 would 
be a reasonable first guess. Note, however, that a program like the 
POWHEG BOX uses a <ei>pT</ei> definition for ISR and FSR that does not 
quite agree with the PYTHIA evolution scale, and thus there will be some 
amount of mismatch. In more sophisticated descriptions, therefore, 
option 2 could be combined with <code>UserHooks</code> vetoes on emissions 
that would lead to double-counting, using more flexible phase space 
boundaries. Further details are found in the 
<aloc href="MatchingAndMerging">Matching and Merging</aloc> description, 
with an example in <code>examples/main31</code>. 
Option 0, finally, may be most realistic when only Born-level 
processes are involved, possibly in combination with a nonzero 
<code>SpaceShower:pTdampMatch</code>. The rules used for avoiding 
double-counting are not foolproof, however. As an example, for the 
<ei>t</ei>-channel process <ei>gamma gamma &rarr; e^+ e^-</ei> its <ei>pT</ei> 
scale is the plausible upper shower limit, with only dampened emissions 
above it. But the initial state is not checked and, had only incoming 
quarks and gluons been taken into account, only the <ei>s</ei>-channel 
process <ei>q qbar &rarr; gamma^*/Z^0 &rarr; e^+ e^-</ei> would have 
been possible, where indeed the whole phase space should be populated. 
So this is erroneously used, giving too much emissions. 
<br/><b>Note 3:</b> These options only apply to the hard interaction. 
If a "second hard" process is present, the two are analyzed and 
set separately for the default 0 option, while both are affected 
the same way for non-default options 1 and 2. 
Emissions off subsequent multiparton interactions are always constrained 
to be below the factorization scale of each process itself. 
 
<br/><br/><table><tr><td><strong>SpaceShower:pTmaxFudge </td><td></td><td> <input type="text" name="2" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 2.0</code>)</td></tr></table>
In cases where the above <code>pTmaxMatch</code> rules would imply 
that <i>pT_max = pT_factorization</i>, <code>pTmaxFudge</code> 
introduces a multiplicative factor <i>f</i> such that instead 
<i>pT_max = f * pT_factorization</i>. Only applies to the hardest 
interaction in an event, and a "second hard" if there is such a one, 
cf. below. It is strongly suggested that <i>f = 1</i>, but variations 
around this default can be useful to test this assumption. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:pTmaxFudgeMPI </td><td></td><td> <input type="text" name="3" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 2.0</code>)</td></tr></table>
A multiplicative factor <i>f</i> such that 
<i>pT_max = f * pT_factorization</i>, as above, but here for the 
non-hardest interactions (when multiparton interactions are allowed). 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:pTdampMatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
These options only take effect when a process is allowed to radiate up 
to the kinematical limit by the above <code>pTmaxMatch</code> choice, 
and no matrix-element corrections are available. Then, in many processes, 
the fall-off in <ei>pT</ei> will be too slow by one factor of <ei>pT^2</ei>. 
That is, while showers have an approximate <ei>dpT^2/pT^2</ei> shape, often 
it should become more like <ei>dpT^2/pT^4</ei> at <ei>pT</ei> values above 
the scale of the hard process. Whether this actually is the case 
depends on the particular process studied, e.g. if <ei>t</ei>-channel 
gluon exchange is likely to dominate. If so, the options below could 
provide a reasonable high-<ei>pT</ei> behaviour without requiring 
higher-order calculations. 
<br/>
<input type="radio" name="4" value="0"><strong>0 </strong>: emissions go up to the kinematical limit,  with no special dampening.  <br/>
<input type="radio" name="4" value="1"><strong>1 </strong>: emissions go up to the kinematical limit,  but dampened by a factor <ei>k^2 Q^2_fac/(pT^2 + k^2 Q^2_fac)</ei>,  where <ei>Q_fac</ei> is the factorization scale and <ei>k</ei> is a  multiplicative fudge factor stored in <code>pTdampFudge</code> below.  <br/>
<input type="radio" name="4" value="2"><strong>2 </strong>: emissions go up to the kinematical limit,  but dampened by a factor <ei>k^2 Q^2_ren/(pT^2 + k^2 Q^2_ren)</ei>,  where <ei>Q_ren</ei> is the renormalization scale and <ei>k</ei> is a  multiplicative fudge factor stored in <code>pTdampFudge</code> below.  <br/>
<input type="radio" name="4" value="3" checked="checked"><strong>3 </strong>: as option 1, but in addition to the standard requirements  for dampening it is further necessary to have ar least two top or  beyond-the-Standard-Model coloured particles in the final state.  Examples include <ei>t tbar</ei> and <ei>squark gluino</ei> production.  <br/>
<input type="radio" name="4" value="4"><strong>4 </strong>: as option 2, but in addition to the standard requirements  for dampening it is further necessary to have ar least two top or  beyond-the-Standard-Model coloured particles in the final state.  Examples include <ei>t tbar</ei> and <ei>squark gluino</ei> production.  <br/>
<br/><b>Note:</b> These options only apply to the hard interaction. 
Specifically, a "second hard" interaction would not be affected. 
Emissions off subsequent multiparton interactions are always constrained 
to be below the factorization scale of the process itself. 
 
<br/><br/><table><tr><td><strong>SpaceShower:pTdampFudge </td><td></td><td> <input type="text" name="5" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 4.0</code>)</td></tr></table>
In cases 1 and 2 above, where a dampening is imposed at around the 
factorization or renormalization scale, respectively, this allows the 
<i>pT</i> scale of dampening of radiation by a half to be shifted 
by this factor relative to the default <i>Q_fac</i> or <i>Q_ren</i>. 
This number ought to be in the neighbourhood of unity, but variations 
away from this value could do better in some processes. 
   
 
<p/> 
The amount of QCD radiation in the shower is determined by 
<br/><br/><table><tr><td><strong>SpaceShower:alphaSvalue </td><td></td><td> <input type="text" name="6" value="0.1365" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1365</strong></code>; <code>minimum = 0.06</code>; <code>maximum = 0.25</code>)</td></tr></table>
The <i>alpha_strong</i> value at scale <code>M_Z^2</code>. 
   
 
<p/> 
The actual value is then regulated by the running to the scale 
<i>pT^2</i>, at which it is evaluated 
<br/><br/><table><tr><td><strong>SpaceShower:alphaSorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Order at which <ei>alpha_strong</ei> runs, 
<br/>
<input type="radio" name="7" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_strong</ei> is kept  fixed.<br/>
<input type="radio" name="7" value="1" checked="checked"><strong>1 </strong>: first order, which is the normal value.<br/>
<input type="radio" name="7" value="2"><strong>2 </strong>: second order. Since other parts of the code do  not go to second order there is no strong reason to use this option,  but there is also nothing wrong with it.<br/>
 
<p/> 
The CMW rescaling of <i>Lambda_QCD</i> (see the section on 
<?php $filepath = $_GET["filepath"];
echo "<a href='StandardModelParameters.php?filepath=".$filepath."' target='page'>";?>StandardModelParameters</a>) 
can be applied to the <i>alpha_strong</i> values used for spacelike showers. 
Note that tunes using this option need lower values of 
<i>alpha_strong(m_Z^2)</i> than tunes that do not. 
<br/><br/><strong>SpaceShower:alphaSuseCMW</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
<br/><code>option </code><strong> off</strong> : Do not apply the CMW rescaling.    
<br/><code>option </code><strong> on</strong> : Apply the CMW rescaling, increasing 
<i>Lambda_QCD</i> for spacelike showers by a factor roughly 1.6. 
   
   
 
<p/> 
QED radiation is regulated by the <i>alpha_electromagnetic</i> 
value at the <i>pT^2</i> scale of a branching. 
 
<br/><br/><table><tr><td><strong>SpaceShower:alphaEMorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = -1</code>; <code>maximum = 1</code>)</td></tr></table>
The running of <ei>alpha_em</ei>. 
<br/>
<input type="radio" name="9" value="1" checked="checked"><strong>1 </strong>: first-order running, constrained to agree with  <code>StandardModel:alphaEMmZ</code> at the <ei>Z^0</ei> mass.  <br/>
<input type="radio" name="9" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept  fixed at its value at vanishing momentum transfer.<br/>
<input type="radio" name="9" value="-1"><strong>-1 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept  fixed, but at <code>StandardModel:alphaEMmZ</code>, i.e. its value  at the <ei>Z^0</ei> mass.  <br/>
 
<p/> 
The natural scale for couplings and PDFs is <i>pT^2</i>. To explore 
uncertainties it is possibly to vary around this value, however, in 
analogy with what can be done for 
<?php $filepath = $_GET["filepath"];
echo "<a href='CouplingsAndScales.php?filepath=".$filepath."' target='page'>";?>hard processes</a>. (Note that 
there is also an <?php $filepath = $_GET["filepath"];
echo "<a href='Variations.php?filepath=".$filepath."' target='page'>";?>automated framework</a> for shower 
uncertainties.) 
 
<br/><br/><table><tr><td><strong>SpaceShower:renormMultFac </td><td></td><td> <input type="text" name="10" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.</code>)</td></tr></table>
The default <i>pT^2</i> renormalization scale is multiplied by 
this prefactor. For QCD this is equivalent to a change of 
<i>Lambda^2</i> in the opposite direction, i.e. to a change of 
<i>alpha_strong(M_Z^2)</i> (except that flavour thresholds 
remain at fixed scales). Below, when <i>pT^2 + pT_0^2</i> is used 
as scale, it is this whole expression that is multiplied by the prefactor. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:factorMultFac </td><td></td><td> <input type="text" name="11" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.</code>)</td></tr></table>
The default <i>pT^2</i> factorization scale is multiplied by 
this prefactor. 
   
 
<p/> 
There are two complementary ways of regularizing the small-<i>pT</i> 
divergence, a sharp cutoff and a smooth dampening. These can be 
combined as desired but it makes sense to coordinate with how the 
same issue is handled in multiparton interactions. 
 
<br/><br/><strong>SpaceShower:samePTasMPI</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Regularize the <i>pT &rarr; 0</i> divergence using the same sharp cutoff 
and smooth dampening parameters as used to describe multiparton interactions. 
That is, the <code>MultipartonInteractions:pT0Ref</code>, 
<code>MultipartonInteractions:ecmRef</code>, 
<code>MultipartonInteractions:ecmPow</code> and 
<code>MultipartonInteractions:pTmin</code> parameters are used to regularize 
all ISR QCD radiation, rather than the corresponding parameters below. 
This is a sensible physics ansatz, based on the assumption that colour 
screening effects influence both MPI and ISR in the same way. Photon 
radiation is regularized separately in either case. 
<br/><b>Note:</b> For photon-photon collisions these parameters are set 
as in <?php $filepath = $_GET["filepath"];
echo "<a href='Photoproduction.php?filepath=".$filepath."' target='page'>";?>Photoproduction</a>. 
<br/><b>Warning:</b> if a large <code>pT0</code> is picked for multiparton 
interactions, such that the integrated interaction cross section is 
below the nondiffractive inelastic one, this <code>pT0</code> will 
automatically be scaled down to cope. Information on such a rescaling 
does NOT propagate to <code>SpaceShower</code>, however. 
   
 
<p/> 
The actual <i>pT0</i> parameter used at a given CM energy scale, 
<i>ecmNow</i>, is obtained from a power law or a logarithmic 
parametrization. The former is default with hadron beams and 
the latter for photon-photon collisions. 
 
<br/><br/><table><tr><td><strong>SpaceShower:pT0parametrization  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
Choice of <ei>pT0</ei> parametrization. 
<br/>
<input type="radio" name="13" value="0" checked="checked"><strong>0 </strong>:  Power law dependence on <ei>ecmNow</ei>:<br/>  <ei>      pT0 = pT0(ecmNow) = pT0Ref * (ecmNow / ecmRef)^ecmPow  </ei>  <br/>
<input type="radio" name="13" value="1"><strong>1 </strong>:  Logarithmic dependence on <ei>ecmNow</ei>:  <eq>      pT0 = pT0(ecmNow) = pT0Ref + ecmPow * log (ecmNow / ecmRef)  </eq>  <br/>
where <ei>pT0Ref</ei>, <ei>ecmRef</ei> and <ei>ecmPow</ei> are the 
three parameters below. 
 
<br/><br/><table><tr><td><strong>SpaceShower:pT0Ref </td><td></td><td> <input type="text" name="14" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 10.0</code>)</td></tr></table>
Regularization of the divergence of the QCD emission probability for 
<i>pT &rarr; 0</i> is obtained by a factor <i>pT^2 / (pT0^2 + pT^2)</i>, 
and by using an <i>alpha_s(pT0^2 + pT^2)</i>. An energy dependence 
of the <i>pT0</i> choice is introduced by the next two parameters, 
so that <i>pT0Ref</i> is the <i>pT0</i> value for the reference 
cm energy, <i>pT0Ref = pT0(ecmRef)</i>. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:ecmRef </td><td></td><td> <input type="text" name="15" value="7000.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7000.0</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
The <i>ecmRef</i> reference energy scale introduced above. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:ecmPow </td><td></td><td> <input type="text" name="16" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 0.5</code>)</td></tr></table>
The <i>ecmPow</i> energy rescaling pace introduced above. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:pTmin </td><td></td><td> <input type="text" name="17" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.0</code>)</td></tr></table>
Lower cutoff in <i>pT</i>, below which no further ISR branchings 
are allowed. Normally the <i>pT0</i> above would be used to 
provide the main regularization of the branching rate for 
<i>pT &rarr; 0</i>, in which case <i>pTmin</i> is used  mainly for 
technical reasons. It is possible, however, to set <i>pT0Ref = 0</i> 
and use <i>pTmin</i> to provide a step-function regularization, 
or to combine them in intermediate approaches. Currently <i>pTmin</i> 
is taken to be energy-independent. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:pTminChgQ </td><td></td><td> <input type="text" name="18" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.01</code>)</td></tr></table>
Parton shower cut-off <i>pT</i> for photon coupling to a coloured 
particle. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:pTminChgL </td><td></td><td> <input type="text" name="19" value="0.0005" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0005</strong></code>; <code>minimum = 0.0001</code>)</td></tr></table>
Parton shower cut-off mass for pure QED branchings. 
Assumed smaller than (or equal to) <i>pTminChgQ</i>. 
   
 
<br/><br/><strong>SpaceShower:rapidityOrder</strong>  <input type="radio" name="20" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="20" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Force emissions, after the first,  to be ordered in rapidity, 
i.e. in terms of decreasing angles in a backwards-evolution sense. 
Could be used to probe sensitivity to unordered emissions. 
Only affects QCD emissions, and only the hard subcollision of an event. 
(For the case "soft QCD" processes the first MPI counts as the hard 
subcollision.) 
 
   
 
<br/><br/><strong>SpaceShower:rapidityOrderMPI</strong>  <input type="radio" name="21" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="21" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Same as the last switch, but this time only emissions in secondary 
scattering systems from MPIs are forced to be ordered in rapidity. 
Each MPI is ordered separately from the others. 
   
 
<a name="section1"></a> 
<h3>Dipole showers</h3> 
 
By default the recoil of an ISR emission is taken by the whole final 
state. The option below gives an alternative approach with local recoils, 
where only one final-state parton takes the recoil of an emission. 
See [<a href="Bibliography.php#refCab17" target="page">Cab17</a>] for further details on the philosophy and 
implementation. 
 
<p/> 
The existing initial-initial global recoil scheme is maintained for 
an emission off a colour line that stretches through the hard process, 
so it is the handling of initial-final dipole ends that is changed. 
Here the single recoiler is picked based on the colour flow of the 
hard process. Additionally the description unifies the emission of 
a gluon from the initial-final and final-initial dipole ends, and 
handles both as part of the ISR framework. Therefore the separation 
into ISR and FSR is not a meaningful classification, and either both 
should be simulated or none. 
 
<p/> 
Note that this option should not be combined with the global option 
for FSR, <code>TimeShower:globalRecoil</code>. Further some settings 
are neglected internally to ensure the same behaviour as obtained for 
<code>TimeShower:allowBeamRecoil = on</code>, 
<code>TimeShower:dampenBeamRecoil = off</code>, and 
<code>SpaceShower:phiIntAsym = off</code>. 
 
<p/> 
The dipole recoil option for the first time allows the simulation 
of Deeply Inelastic Scattering processes in PYTHIA 8, see the 
<code>main36.cc</code> example. Note that the simultaneous emission of 
photons off the lepton leg has not yet been implemented, so you need to 
set <code>PDF:lepton = off</code> and 
<code>TimeShower:QEDshowerByL = off</code>. You are further recommended 
to set <code>SpaceShower:pTmaxMatch = 2</code> to fill the whole phase 
space with parton showers. This is allowed since the shower and 
matrix-element behaviours match well over the whole phase space 
(at least for the first emission). 
 
<br/><br/><strong>SpaceShower:dipoleRecoil</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Option to switch on the dipole-recoil scheme as described above. 
   
 
<a name="section2"></a> 
<h3>Weak showers</h3> 
 
The emission of weak gauge bosons is an integrated part of the initial- 
and final-state radiation, see <?php $filepath = $_GET["filepath"];
echo "<a href='WeakShowers.php?filepath=".$filepath."' target='page'>";?>Weak Showers</a>. 
The following settings are those specifically related to the initial-state 
weak radiation, while common settings are found in the 
<?php $filepath = $_GET["filepath"];
echo "<a href='WeakShowers.php?filepath=".$filepath."' target='page'>";?>Weak Showers</a> description. 
 
<br/><br/><strong>SpaceShower:weakShower</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow a weak shower, yes or no. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:weakShowerMode  </td><td></td><td> <input type="text" name="24" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Determine which branchings are allowed. 
<br/><code>option </code><strong> 0</strong> :  both <i>W^+-</i> and <i>Z^0</i> branchings. 
   
<br/><code>option </code><strong> 1</strong> :  only <i>W^+-</i> branchings.    
<br/><code>option </code><strong> 2</strong> :  only <i>Z^0</i> branchings.    
   
 
<br/><br/><table><tr><td><strong> </td><td></td><td> <input type="text" name="25" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 2.0</code>)</td></tr></table>
Parton shower cut-off <i>pT</i> for weak branchings. 
   
 
<a name="section3"></a> 
<h3>Further variables</h3> 
 
These should normally not be touched. Their only function is for 
cross-checks. 
 
<p/> 
There are three flags you can use to switch on or off selected 
branchings in the shower: 
 
<br/><br/><strong>SpaceShower:QCDshower</strong>  <input type="radio" name="26" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="26" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow a QCD shower; on/off = true/false. 
   
 
<br/><br/><strong>SpaceShower:QEDshowerByQ</strong>  <input type="radio" name="27" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="27" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow quarks to radiate photons; on/off = true/false. 
   
 
<br/><br/><strong>SpaceShower:QEDshowerByL</strong>  <input type="radio" name="28" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="28" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow leptons to radiate photons; on/off = true/false. 
   
 
<p/> 
There are some further possibilities to modify the shower: 
 
<br/><br/><strong>SpaceShower:MEcorrections</strong>  <input type="radio" name="29" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="29" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use of matrix element corrections; on/off = true/false. 
   
 
<br/><br/><strong>SpaceShower:MEafterFirst</strong>  <input type="radio" name="30" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="30" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use of matrix element corrections also after the first emission, 
for dipole ends of the same system that did not yet radiate. 
Only has a meaning if <code>MEcorrections</code> above is 
switched on. 
   
 
<br/><br/><strong>SpaceShower:phiPolAsym</strong>  <input type="radio" name="31" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="31" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Azimuthal asymmetry induced by gluon polarization; on/off = true/false. 
   
 
<br/><br/><strong>SpaceShower:phiPolAsymHard</strong>  <input type="radio" name="32" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="32" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Extend the above azimuthal asymmetry (if on) also back to gluons produced 
in the hard process itself, where feasible; on/off = true/false. 
   
 
<br/><br/><strong>SpaceShower:phiIntAsym</strong>  <input type="radio" name="33" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="33" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Azimuthal asymmetry induced by interference; on/off = true/false. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:strengthIntAsym </td><td></td><td> <input type="text" name="34" value="0.7" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.7</strong></code>; <code>minimum = 0.</code>; <code>maximum = 0.9</code>)</td></tr></table>
Size of asymmetry induced by interference. Natural value of order 0.5; 
expression would blow up for a value of 1. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:nQuarkIn  </td><td></td><td> <input type="text" name="35" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed quark flavours in <i>g &rarr; q qbar</i> branchings, 
when kinematically allowed, and thereby also in incoming beams. 
Changing it to 4 would forbid <i>g &rarr; b bbar</i>, etc. 
   
 
<br/><br/><strong>SpaceShower:useFixedFacScale</strong>  <input type="radio" name="36" value="on"><strong>On</strong>
<input type="radio" name="36" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow the possibility to use a fixed factorization scale, set by 
the <code>parm</code> below. This option is unphysical and only 
intended for toy-model and debug studies. 
   
 
<br/><br/><table><tr><td><strong>SpaceShower:fixedFacScale </td><td></td><td> <input type="text" name="37" value="100." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100.</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
The fixed factorization scale, in GeV, that would be used in the 
evaluation of parton densities if the <code>flag</code> above is on. 
   
 
<a name="section4"></a> 
<h3>Technical notes</h3> 
 
Almost everything is equivalent to the algorithm in 
[<a href="Bibliography.php#refSjo05" target="page">Sjo05</a>,<a href="Bibliography.php#refCor10a" target="page">Cor10a</a>]. Minor changes are as follows. 
<ul> 
<li> 
It is now possible to have a second-order running <i>alpha_s</i>, 
in addition to fixed or first-order running. 
</li> 
<li> 
The description of heavy flavour production in the threshold region 
has been modified, so as to be more forgiving about mismatches 
between the <i>c/b</i>  masses used in Pythia relative to those 
used in a respective PDF parametrization. The basic idea is that, 
in the threshold region of a heavy quark <i>Q</i>, <i>Q = c/b</i>, 
the effect of subsequent <i>Q &rarr; Q g</i> branchings is negligible. 
If so, then 
<br/><i> 
   f_Q(x, pT2) = integral_mQ2^pT2  dpT'2/pT'2 * alpha_s(pT'2)/2pi 
      * integral P(z) g(x', pT'2) delta(x - z x') 
</i><br/> 
so use this to select the <i>pT2</i> of the <i>g &rarr; Q Qbar</i> 
branching. In the old formalism the same kind of behaviour should 
be obtained, but by a cancellation of a <i>1/f_Q</i> that diverges 
at the threshold and a Sudakov that vanishes. 
<br/> 
The strategy therefore is that, once <i>pT2 &lt; f * mQ2</i>, with 
<i>f</i> a parameter of the order of 2, a <i>pT2</i> is chosen 
like <i>dpT2/pT2</i> between <i>mQ2</i> and <i>f * mQ2</i>, a 
nd a <i>z</i> flat in the allowed range. Thereafter acceptance 
is based on the product of three factors, representing the running 
of <i>alpha_strong</i>, the splitting kernel (including the mass term) 
and the gluon density weight. At failure, a new <i>pT2</i> is chosen 
in the same  range, i.e. is not required to be lower since no Sudakov 
is involved. 
</li> 
<li> 
The QED algorithm now allows for hadron beams with non-zero photon 
content. The backwards-evolution of a photon in a hadron is identical 
to that of a gluon, with <i>CF &rarr; eq^2</i> and <i>CA &rarr; 0</i>. 
Note that this will only work in conjunction with parton distributions 
that explicitly include photons as part of the hadron structure, such 
as the NNPDF2.3 QCD+QED sets. The possibility of a fermion 
backwards-evolving to a photon has not yet been included, nor has 
photon backwards-evolution in lepton beams. 
</li> 
</ul> 
 
<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "0")
{
$data = "SpaceShower:pTmaxMatch = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "1.0")
{
$data = "SpaceShower:pTmaxFudge = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1.0")
{
$data = "SpaceShower:pTmaxFudgeMPI = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "3")
{
$data = "SpaceShower:pTdampMatch = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.0")
{
$data = "SpaceShower:pTdampFudge = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.1365")
{
$data = "SpaceShower:alphaSvalue = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1")
{
$data = "SpaceShower:alphaSorder = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "SpaceShower:alphaSuseCMW = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1")
{
$data = "SpaceShower:alphaEMorder = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "1.")
{
$data = "SpaceShower:renormMultFac = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1.")
{
$data = "SpaceShower:factorMultFac = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "SpaceShower:samePTasMPI = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0")
{
$data = "SpaceShower:pT0parametrization = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "2.0")
{
$data = "SpaceShower:pT0Ref = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "7000.0")
{
$data = "SpaceShower:ecmRef = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.0")
{
$data = "SpaceShower:ecmPow = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.2")
{
$data = "SpaceShower:pTmin = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0.5")
{
$data = "SpaceShower:pTminChgQ = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.0005")
{
$data = "SpaceShower:pTminChgL = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "on")
{
$data = "SpaceShower:rapidityOrder = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "on")
{
$data = "SpaceShower:rapidityOrderMPI = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "SpaceShower:dipoleRecoil = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "SpaceShower:weakShower = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0")
{
$data = "SpaceShower:weakShowerMode = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "1.0")
{
$data = " = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "on")
{
$data = "SpaceShower:QCDshower = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "on")
{
$data = "SpaceShower:QEDshowerByQ = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "on")
{
$data = "SpaceShower:QEDshowerByL = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "on")
{
$data = "SpaceShower:MEcorrections = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "on")
{
$data = "SpaceShower:MEafterFirst = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "on")
{
$data = "SpaceShower:phiPolAsym = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "on")
{
$data = "SpaceShower:phiPolAsymHard = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "on")
{
$data = "SpaceShower:phiIntAsym = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "0.7")
{
$data = "SpaceShower:strengthIntAsym = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "5")
{
$data = "SpaceShower:nQuarkIn = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "off")
{
$data = "SpaceShower:useFixedFacScale = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "100.")
{
$data = "SpaceShower:fixedFacScale = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
