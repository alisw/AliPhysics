<html>
<head>
<title>Timelike Showers</title>
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

<form method='post' action='TimelikeShowers.php'>
 
<h2>Timelike Showers</h2> 
 
The PYTHIA algorithm for timelike final-state showers is based on 
the article [<a href="Bibliography.php" target="page">Sjo05</a>], where a transverse-momentum-ordered 
evolution scheme is introduced, with the extension to fully interleaved 
evolution covered in [<a href="Bibliography.php" target="page">Cor10a</a>]. This algorithm is influenced by 
the previous mass-ordered algorithm in PYTHIA [<a href="Bibliography.php" target="page">Ben87</a>] and by 
the dipole-emission formulation in Ariadne [<a href="Bibliography.php" target="page">Gus86</a>]. From the 
mass-ordered algorithm it inherits a merging procedure for first-order 
gluon-emission matrix elements in essentially all two-body decays 
in the standard model and its minimal supersymmetric extension 
[<a href="Bibliography.php" target="page">Nor01</a>]. 
 
<p/> 
The normal user is not expected to call <code>TimeShower</code> directly, 
but only have it called from <code>Pythia</code>. Some of the parameters 
below, in particular <code>TimeShower:alphaSvalue</code>, would be of 
interest for a tuning exercise, however. 
 
<h3>Main variables</h3> 
 
Often the maximum scale of the FSR shower evolution is understood from the 
context. For instance, in a resonance decay half the resonance mass sets an 
absolute upper limit. For a hard process in a hadronic collision the choice 
is not as unique. Here the <?php $filepath = $_GET["filepath"];
echo "<a href='CouplingsAndScales.php?filepath=".$filepath."' target='page'>";?>factorization 
scale</a> has been chosen as the maximum evolution scale. This would be 
the <i>pT</i> for a <i>2 &rarr; 2</i> process, supplemented by mass terms 
for massive outgoing particles. For some special applications we do allow 
an alternative. 
 
<br/><br/><table><tr><td><strong>TimeShower:pTmaxMatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Way in which the maximum shower evolution scale is set to match the 
scale of the hard process itself. 
<br/>
<input type="radio" name="1" value="0"><strong>0 </strong>: <b>(i)</b> if the final state of the hard process  (not counting subsequent resonance decays) contains at least one quark  (<ei>u, d, s, c ,b</ei>), gluon or photon then <ei>pT_max</ei>  is chosen to be the factorization scale for internal processes  and the <code>scale</code> value for Les Houches input;  <b>(ii)</b> if not, emissions are allowed to go all the way up to  the kinematical limit (i.e. to half the dipole mass).  This option agrees with the corresponding one for  <aloc href="SpacelikeShowers">spacelike showers</aloc>. There the  reasoning is that in the former set of processes the ISR  emission of yet another quark, gluon or photon could lead to  double-counting, while no such danger exists in the latter case.  The argument is less compelling for timelike showers, but could  be a reasonable starting point.  <br/>
<input type="radio" name="1" value="1" checked="checked"><strong>1 </strong>: always use the factorization scale for an internal  process and the <code>scale</code> value for Les Houches input,  i.e. the lower value. This should avoid double-counting, but  may leave out some emissions that ought to have been simulated.  (Also known as wimpy showers.)  <br/>
<input type="radio" name="1" value="2"><strong>2 </strong>: always allow emissions up to the kinematical limit  (i.e. to half the dipole mass). This will simulate all possible event  topologies, but may lead to double-counting.  (Also known as power showers.)  <br/>
<br/><b>Note 1:</b> as enumerated in the text, these options take effect 
both for internal and external processes. Whether a particular option 
makes sense depends on the context. For instance, if events for the same 
basic process to different orders are to be matched, then option 1 would 
be a reasonable first guess. But in more sophisticated descriptions 
option 2 could be combined with <code>UserHooks</code> vetoes on 
emissions that would lead to double-counting, using more flexible 
phase space boundaries. Further details are found in the 
<aloc href="MatchingAndMerging">Matching and Merging</aloc> description, 
with an example in <code>examples/main31</code>. 
Option 0, finally, may be most realistic when only Born-level processes 
are involved, possibly in combination with a nonzero 
<code>TimeShower:pTdampMatch</code>. 
<br/><b>Note 2:</b> These options only apply to the hard interaction. 
If a "second hard" process is present, the two are analyzed and 
set separately for the default 0 option, while both are affected 
the same way for non-default options 1 and 2. 
Emissions off subsequent multiparton interactions are always constrained 
to be below the factorization scale of each process itself. The options 
also assume that you use interleaved evolution, so that FSR is in direct 
competition with ISR for the hardest emission. If you already 
generated a number of ISR partons at low <ei>pT</ei>, it would not 
make sense to have a later FSR shower up to the kinematical limit 
for all of them. 
<br/><b>Note 3:</b> Recall that resonance decays are not affected by 
this mode, but that showers there are always set to fill the full phase 
space, often with built-in matrix-element-matching that give a NLO 
accuracy. A modification of this behaviour would require you to work with 
<code>UserHooks</code>. However, for Les Houches input the optional 
<code><aloc href="BeamParameters">Beams:strictLHEFscale = on</aloc></code> 
setting restricts all emissions, also in resonance decays, to be below 
the input <code>scale</code> value. 
 
<br/><br/><table><tr><td><strong>TimeShower:pTmaxFudge </td><td></td><td> <input type="text" name="2" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 2.0</code>)</td></tr></table>
In cases where the above <code>pTmaxMatch</code> rules would imply 
that <i>pT_max = pT_factorization</i>, <code>pTmaxFudge</code> 
introduces a multiplicative factor <i>f</i> such that instead 
<i>pT_max = f * pT_factorization</i>. Only applies to the hardest 
interaction in an event, and a "second hard" if there is such a one, 
cf. below. It is strongly suggested that <i>f = 1</i>, but variations 
around this default can be useful to test this assumption. 
<br/><b>Note:</b>Scales for resonance decays are not affected, but can 
be set separately by <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>user hooks</a>. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:pTmaxFudgeMPI </td><td></td><td> <input type="text" name="3" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 2.0</code>)</td></tr></table>
A multiplicative factor <i>f</i> such that 
<i>pT_max = f * pT_factorization</i>, as above, but here for the 
non-hardest interactions (when multiparton interactions are allowed). 
   
 
<br/><br/><table><tr><td><strong>TimeShower:pTdampMatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 4</code>)</td></tr></table>
These options only take effect when a process is allowed to radiate up 
to the kinematical limit by the above <code>pTmaxMatch</code> choice, 
and no matrix-element corrections are available. Then, in many processes, 
the fall-off in <ei>pT</ei> will be too slow by one factor of <ei>pT^2</ei>. 
That is, while showers have an approximate <ei>dpT^2/pT^2</ei> shape, often 
it should become more like <ei>dpT^2/pT^4</ei> at <ei>pT</ei> values above 
the scale of the hard process. This argument is more obvious and relevant 
for ISR, where emissions could go the the kinematical limit, whereas they 
are constrained by the respective dipole mass for FSR. Nevertheless this 
matching option is offered for FSR to have a (semi-)symmetric description. 
Note that a dampening factor is applied to all dipoles in the final state 
of the hard process, which is somewhat different from the ISR implementation. 
<br/>
<input type="radio" name="4" value="0" checked="checked"><strong>0 </strong>: emissions go up to the kinematical limit,  with no special dampening.  <br/>
<input type="radio" name="4" value="1"><strong>1 </strong>: emissions go up to the kinematical limit,  but dampened by a factor <ei>k^2 Q^2_fac/(pT^2 + k^2 Q^2_fac)</ei>,  where <ei>Q_fac</ei> is the factorization scale and <ei>k</ei> is a  multiplicative fudge factor stored in <code>pTdampFudge</code> below.  <br/>
<input type="radio" name="4" value="2"><strong>2 </strong>: emissions go up to the kinematical limit,  but dampened by a factor <ei>k^2 Q^2_ren/(pT^2 + k^2 Q^2_ren)</ei>,  where <ei>Q_ren</ei> is the renormalization scale and <ei>k</ei> is a  multiplicative fudge factor stored in <code>pTdampFudge</code> below.  <br/>
<input type="radio" name="4" value="3"><strong>3 </strong>: as option 1, but in addition to the standard requirements for dampening it is further necessary to have ar least two top or  beyond-the-Standard-Model coloured particles in the final state.  Examples include <ei>t tbar</ei> and <ei>squark gluino</ei> production.   <br/>
<input type="radio" name="4" value="4"><strong>4 </strong>: as option 2, but in addition to the standard requirements for dampening it is further necessary to have ar least two top or  beyond-the-Standard-Model coloured particles in the final state.  Examples include <ei>t tbar</ei> and <ei>squark gluino</ei> production.  <br/>
<br/><b>Note:</b> These options only apply to the hard interaction. 
Specifically, a "second hard" interaction would not be affected. 
Emissions off subsequent multiparton interactions are always constrained 
to be below the factorization scale of the process itself. 
 
<br/><br/><table><tr><td><strong>TimeShower:pTdampFudge </td><td></td><td> <input type="text" name="5" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 4.0</code>)</td></tr></table>
In cases 1 and 2 above, where a dampening is imposed at around the 
factorization or renormalization scale, respectively, this allows the 
<i>pT</i> scale of dampening of radiation by a half to be shifted 
by this factor relative to the default <i>Q_fac</i> or <i>Q_ren</i>. 
This number ought to be in the neighbourhood of unity, but variations 
away from this value could do better in some processes. 
   
 
<p/> 
The amount of QCD radiation in the shower is determined by 
<br/><br/><table><tr><td><strong>TimeShower:alphaSvalue </td><td></td><td> <input type="text" name="6" value="0.1365" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1365</strong></code>; <code>minimum = 0.06</code>; <code>maximum = 0.25</code>)</td></tr></table>
The <i>alpha_strong</i> value at scale <i>M_Z^2</i>. The default 
value corresponds to a crude tuning to LEP data, to be improved. 
   
 
<p/> 
The actual value is then regulated by the running to the scale 
<i>pT^2</i>, at which the shower evaluates <i>alpha_strong</i>. 
 
<br/><br/><table><tr><td><strong>TimeShower:alphaSorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Order at which <ei>alpha_strong</ei> runs, 
<br/>
<input type="radio" name="7" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_strong</ei> is kept  fixed.<br/>
<input type="radio" name="7" value="1" checked="checked"><strong>1 </strong>: first order, which is the normal value.<br/>
<input type="radio" name="7" value="2"><strong>2 </strong>: second order. Since other parts of the code do  not go to second order there is no strong reason to use this option,  but there is also nothing wrong with it.<br/>
 
<p/> 
The CMW rescaling of <i>Lambda_QCD</i> (see the section on 
<?php $filepath = $_GET["filepath"];
echo "<a href='StandardModelParameters.php?filepath=".$filepath."' target='page'>";?>StandardModelParameters</a>) 
can be applied to the <i>alpha_strong</i> values used for 
timelike showers. Note that tunes using this option need lower values of 
<i>alpha_strong(m_Z^2)</i> than tunes that do not. 
<br/><br/><strong>TimeShower:alphaSuseCMW</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>false</strong></code>)<br/>
<br/><code>option </code><strong> false</strong> : Do not apply the CMW rescaling.    
<br/><code>option </code><strong> true</strong> : Apply the CMW rescaling, increasing 
 <i>Lambda_QCD</i> for timelike showers by a factor roughly 1.6. 
   
   
 
<p/> 
QED radiation is regulated by the <i>alpha_electromagnetic</i> 
value at the <i>pT^2</i> scale of a branching. 
 
<br/><br/><table><tr><td><strong>TimeShower:alphaEMorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = -1</code>; <code>maximum = 1</code>)</td></tr></table>
The running of <ei>alpha_em</ei>. 
<br/>
<input type="radio" name="9" value="1" checked="checked"><strong>1 </strong>: first-order running, constrained to agree with  <code>StandardModel:alphaEMmZ</code> at the <ei>Z^0</ei> mass.  <br/>
<input type="radio" name="9" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept  fixed at its value at vanishing momentum transfer.<br/>
<input type="radio" name="9" value="-1"><strong>-1 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept  fixed, but at <code>StandardModel:alphaEMmZ</code>, i.e. its value  at the <ei>Z^0</ei> mass.  <br/>
 
<p/> 
The natural scale for couplings, and PDFs for dipoles stretching out 
to the beam remnants, is <i>pT^2</i>. To explore uncertainties it 
is possibly to vary around this value, however, in analogy with what 
can be done for <?php $filepath = $_GET["filepath"];
echo "<a href='CouplingsAndScales.php?filepath=".$filepath."' target='page'>";?>hard processes</a>. 
 
<br/><br/><table><tr><td><strong>TimeShower:renormMultFac </td><td></td><td> <input type="text" name="10" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.</code>)</td></tr></table>
The default <i>pT^2</i> renormalization scale is multiplied by 
this prefactor. For QCD this is equivalent to a change of 
<i>Lambda^2</i> in the opposite direction, i.e. to a change of 
<i>alpha_strong(M_Z^2)</i> (except that flavour thresholds 
remain at fixed scales). 
   
 
<br/><br/><table><tr><td><strong>TimeShower:factorMultFac </td><td></td><td> <input type="text" name="11" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.</code>)</td></tr></table>
The default <i>pT^2</i> factorization scale is multiplied by 
this prefactor. 
   
 
<p/> 
The rate of radiation if divergent in the <i>pT &rarr; 0</i> limit. Here, 
however, perturbation theory is expected to break down. Therefore an 
effective <i>pT_min</i> cutoff parameter is introduced, below which 
no emissions are allowed. The cutoff may be different for QCD and QED 
radiation off quarks, and is mainly a technical parameter for QED 
radiation off leptons. 
 
<br/><br/><table><tr><td><strong>TimeShower:pTmin </td><td></td><td> <input type="text" name="12" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 2.0</code>)</td></tr></table>
Parton shower cut-off <i>pT</i> for QCD emissions. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:pTminChgQ </td><td></td><td> <input type="text" name="13" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 2.0</code>)</td></tr></table>
Parton shower cut-off <i>pT</i> for photon coupling to coloured particle. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:pTminChgL </td><td></td><td> <input type="text" name="14" value="1e-6" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1e-6</strong></code>; <code>minimum = 1e-10</code>; <code>maximum = 2.0</code>)</td></tr></table>
Parton shower cut-off <i>pT</i> for pure QED branchings. 
Assumed smaller than (or equal to) <code>pTminChgQ</code>. 
   
 
<p/> 
Shower branchings <i>gamma &rarr; f fbar</i>, where <i>f</i> is a 
quark or lepton, in part compete with the hard processes involving 
<i>gamma^*/Z^0</i> production. In order to avoid overlap it makes 
sense to correlate the maximum <i>gamma</i> mass allowed in showers 
with the minimum <i>gamma^*/Z^0</i> mass allowed in hard processes. 
In addition, the shower contribution only contains the pure 
<i>gamma^*</i> contribution, i.e. not the <i>Z^0</i> part, so 
the mass spectrum above 50 GeV or so would not be well described. 
 
<br/><br/><table><tr><td><strong>TimeShower:mMaxGamma </td><td></td><td> <input type="text" name="15" value="10.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10.0</strong></code>; <code>minimum = 0.001</code>; <code>maximum = 5000.0</code>)</td></tr></table>
Maximum invariant mass allowed for the created fermion pair in a 
<i>gamma &rarr; f fbar</i> branching in the shower. 
   
 
<h3>Interleaved evolution</h3> 
 
Multiparton interactions (MPI) and initial-state showers (ISR) are 
always interleaved, as follows. Starting from the hard interaction, 
the complete event is constructed by a set of steps. In each step 
the <i>pT</i> scale of the previous step is used as starting scale 
for a downwards evolution. The MPI and ISR components each make 
their respective Monte Carlo choices for the next lower <i>pT</i> 
value. The one with larger <i>pT</i> is allowed to carry out its 
proposed action, thereby modifying the conditions for the next steps. 
This is relevant since the two components compete for the energy 
contained in the beam remnants: both an interaction and an emission 
take away some of the energy, leaving less for the future. The end 
result is a combined chain of decreasing <i>pT</i> values, where 
ones associated with new interactions and ones with new emissions 
are interleaved. 
 
<p/> 
There is no corresponding requirement for final-state radiation (FSR) 
to be interleaved. Such an FSR emission does not compete directly for 
beam energy (but see below), and also can be viewed as occurring after 
the other two components in some kind of time sense. Interleaving is 
allowed, however, since it can be argued that a high-<i>pT</i> FSR 
occurs on shorter time scales than a low-<i>pT</i> MPI, say. 
Backwards evolution of ISR is also an example that physical time 
is not the only possible ordering principle, but that one can work 
with conditional probabilities: given the partonic picture at a 
specific <i>pT</i> resolution scale, what possibilities are open 
for a modified picture at a slightly lower <i>pT</i> scale, either 
by MPI, ISR or FSR? Complete interleaving of the three components also 
offers advantages if one aims at matching to higher-order matrix 
elements above some given scale. 
 
<br/><br/><strong>TimeShower:interleave</strong>  <input type="radio" name="16" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="16" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on, final-state emissions are interleaved in the same 
decreasing-<i>pT</i> chain as multiparton interactions and initial-state 
emissions. If off, final-state emissions are only addressed after the 
multiparton interactions and initial-state radiation have been considered. 
   
 
<p/> 
As an aside, it should be noted that such interleaving does not affect 
showering in resonance decays, such as a <i>Z^0</i>. These decays are 
only introduced after the production process has been considered in full, 
and the subsequent FSR is carried out inside the resonance, with 
preserved resonance mass. 
 
<p/> 
One aspect of FSR for a hard process in hadron collisions is that often 
colour dipoles are formed between a scattered parton and a beam remnant, 
or rather the hole left behind by an incoming partons. If such holes 
are allowed as dipole ends and take the recoil when the scattered parton 
undergoes a branching then this translates into the need to take some 
amount of remnant energy also in the case of FSR, i.e. the roles of 
ISR and FSR are not completely decoupled. The energy taken away is 
bookkept by increasing the <i>x</i> value assigned to the incoming 
scattering parton, and a reweighting factor 
<i>x_new f(x_new, pT^2) / x_old f(x_old, pT^2)</i> 
in the emission probability ensures that not unphysically large 
<i>x_new</i> values are reached. Usually such <i>x</i> changes are 
small, and they can be viewed as a higher-order effect beyond the 
accuracy of the leading-log initial-state showers. 
 
<p/> 
This choice is not unique, however. As an alternative, if nothing else 
useful for cross-checks, one could imagine that the FSR is completely 
decoupled from the ISR and beam remnants. 
 
<br/><br/><strong>TimeShower:allowBeamRecoil</strong>  <input type="radio" name="17" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="17" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on, the final-state shower is allowed to borrow energy from 
the beam remnants as described above, thereby changing the mass of the 
scattering subsystem. If off, the partons in the scattering subsystem 
are constrained to borrow energy from each other, such that the total 
four-momentum of the system is preserved. This flag has no effect 
on resonance decays, where the shower always preserves the resonance 
mass, cf. the comment above about showers for resonances never being 
interleaved. 
   
 
<br/><br/><strong>TimeShower:dampenBeamRecoil</strong>  <input type="radio" name="18" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="18" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When beam recoil is allowed there is still some ambiguity how far 
into the beam end of the dipole that emission should be allowed. 
It is dampened in the beam region, but probably not enough. 
When on an additional suppression factor 
<i>4 pT2_hard / (4 pT2_hard + m2)</i> is multiplied on to the 
emission probability. Here <i>pT_hard</i> is the transverse momentum 
of the radiating parton and <i>m</i> the off-shell mass it acquires 
by the branching, <i>m2 = pT2/(z(1-z))</i>. Note that 
<i>m2 = 4 pT2_hard</i> is the kinematical limit for a scattering 
at 90 degrees without beam recoil. 
   
 
<p/> 
When there is no interleaving, a number of MPIs may have been generated
before FSR is considered. In principle there could be colour correlations 
between the MPIs, such that a final-state colour of one MPI could be
matched by the corresponding final-state anticolour of another MPI.
These thereby would form a colour dipole, but one that does not come out
from a common vertex, and therefore presumably could not radiate in full.
Currently the standard procedure is to match colours between MPIs
only after FSR, so MPI systems would radiate independently, with 
recoils taken by the beam remnant, where necessary. This could change,
however, and the following switch would then regulate the choice of 
behaviour.
 
<br/><br/><strong>TimeShower:allowMPIdipole</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, and if interleaving is off, then dipoles are allowed to be
formed between matching final-state colour-anticolour pairs also
between two different MPIs. Else dipoles can normally only form
inside the same MPI, and the could-have-been dipoles between different
MPIs instead appear as dipoles stretched to the beam remnants. 
In either case a dipole can still form between two MPIs if a final-state 
colour cannot be matched inside the same MPI. This should normally
not happen, except if rescattering is allowed, whereby two or more 
MPIs get interconnected.
   
 
<h3>Global recoil</h3> 
 
The final-state algorithm is based on dipole-style recoils, where 
one single parton takes the full recoil of a branching. This is unlike 
the initial-state algorithm, where the complete already-existing 
final state shares the recoil of each new emission. As an alternative, 
also the final-state algorithm contains an option where the recoil 
is shared between all partons in the final state. Thus the radiation 
pattern is unrelated to colour correlations. This is especially 
convenient for some matching algorithms, like MC@NLO, where a full 
analytic knowledge of the shower radiation pattern is needed to avoid 
double-counting. (The <i>pT</i>-ordered shower is described in 
[<a href="Bibliography.php" target="page">Sjo05</a>], and the corrections for massive radiator and recoiler 
in [<a href="Bibliography.php" target="page">Nor01</a>].) 
 
<p/> 
Technically, the radiation pattern is most conveniently represented 
in the rest frame of the final state of the hard subprocess. Then, for 
each parton at a time, the rest of the final state can be viewed as 
a single effective parton. This "parton" has a fixed invariant mass 
during the emission process, and takes the recoil without any changed 
direction of motion. The momenta of the individual new recoilers are 
then obtained by a simple common boost of the original ones. 
 
<p/> 
This alternative approach will miss out on the colour coherence 
phenomena. Specifically, with the whole subcollision mass as "dipole" 
mass, the phase space for subsequent emissions is larger than for 
the normal dipole algorithm. The phase space difference grows as 
more and more gluons are created, and thus leads to a way too steep 
multiplication of soft gluons. Therefore the main application is 
for the first one or few emissions of the shower, where a potential 
overestimate of the emission rate is to be corrected for anyway, 
by matching to the relevant matrix elements. Thereafter, subsequent 
emissions should be handled as before, i.e. with dipoles spanned 
between nearby partons. Furthermore, only the first (hardest) 
subcollision is handled with global recoils, since subsequent MPI's 
would not be subject to matrix element corrections anyway. 
 
<p/> 
In order for the mid-shower switch from global to local recoils 
to work, colours are traced and bookkept just as for normal showers; 
it is only that this information is not used in those steps where 
a global recoil is requested. (Thus, e.g., a gluon is still bookkept 
as one colour and one anticolour dipole end, with half the charge 
each, but with global recoil those two ends radiate identically.) 
 
<br/><br/><strong>TimeShower:globalRecoil</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Alternative approach as above, where all final-state particles share 
the recoil of an emission. 
<br/>If off, then use the standard dipole-recoil approach. 
<br/>If on, use the alternative global recoil, but only for the first 
interaction, and only while the number of particles in the final state 
is at most <code>TimeShower:nMaxGlobalRecoil</code> before the 
branching. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:nMaxGlobalRecoil  </td><td></td><td> <input type="text" name="21" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>)</td></tr></table>
Represents the maximum number of particles in the final state for which 
the next final-state emission can be performed with the global recoil 
strategy. This number counts all particles, whether they are 
allowed to radiate or not, e.g. also <i>Z^0</i>. Also partons 
created by initial-state radiation emissions counts towards this sum, 
as part of the interleaved evolution. Without interleaved evolution 
this option would not make sense, since then a varying and large 
number of partons could already have been created by the initial-state 
radiation before the first final-state one, and then there is not 
likely to be any matrix elements available for matching. 
   
 
<p/> 
Two variations of the scheme outlined above are also available, 
(motivated by comparative studies within aMC@NLO). These studies indicate 
that global recoils should be used as sparsely as possible, in order to 
retain desirable features of the radiation pattern produced with the local 
recoil prescription. 
 
<br/><br/><table><tr><td><strong>TimeShower:globalRecoilMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Choice which splittings are produced with the global recoil approach. 
<br/>
<input type="radio" name="22" value="0" checked="checked"><strong>0 </strong>: Global recoil mode as outlined above, i.e. using global  recoils until the number of final state particles exceeds  <code>TimeShower:nMaxGlobalRecoil</code>.<br/>
<input type="radio" name="22" value="1"><strong>1 </strong>: Global recoil only for the first branching of  final state legs that have an ancestor in the hard process, and  if the maximal number of branchings generated according to the global  recoil scheme (see <code>TimeShower:nMaxGlobalBranch</code> below) has  not yet been reached.<br/>
<input type="radio" name="22" value="2"><strong>2 </strong>: Global recoil only if the first branching in  the whole evolution is a timelike splitting of a parton in an  event with Born-like kinematics (i.e.\ an S-event).  The impact of global recoils should be minimal in this case.  This option is only sensible for interleaved evolution.  <br/>
 
<br/><br/><table><tr><td><strong>TimeShower:nMaxGlobalBranch  </td><td></td><td> <input type="text" name="23" value="-1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1</strong></code>)</td></tr></table>
The maximum number of splittings in the final state for which 
the next final-state emission can be performed with the global recoil 
strategy. This number has to be set if <code>TimeShower:globalRecoilMode = 1 
</code> or <code>TimeShower:globalRecoilMode = 2</code> 
   
 
<br/><br/><table><tr><td><strong>TimeShower:nPartonsInBorn  </td><td></td><td> <input type="text" name="24" value="-1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1</strong></code>)</td></tr></table>
The number of partons for Born-like phase space points. This number needs 
to be set if a different treatment of S-events (with Born-like kinematics) 
and H-events (with real-emission kinematics) is desired. This number has 
to be set if <code>TimeShower:globalRecoilMode = 2</code>. 
   
 
<br/><br/><strong>TimeShower:limitPTmaxGlobal</strong>  <input type="radio" name="25" value="on"><strong>On</strong>
<input type="radio" name="25" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, limit the maximal pT produced in branchings in the global recoil scheme 
exactly as in the default (local) scheme. This means that the mass of the 
splitting dipole will set an upper bound for the pT of an emission. 
To be more explicit, this disallows emissions with pT larger than 
<i>min{&mu;<sub>start </sub><sup>2</sup>, m<sub>D</sub><sup>2</sup>/4}</i>, 
with <i>m<sub>D</sub><sup>2</sup> = 
(&radic;<span style="text-decoration:overline;color:green">&nbsp;(p<sub>r 
</sub>+p<sub>s</sub>)<sup>2</sup>&nbsp;</span>-m<sub>0,s</sub>)<sup>2</sup> 
- m<sub>0,r</sub><sup>2</sup> </i>, where 
the shower starting scale is <i>&mu;<sub>start</sub></i> (i.e. SCALUP when 
reading LHE files, and <code> Info.QFac()</code> otherwise), <i>r</i> the 
radiating parton, and <i>s</i> the recoiling particle that would have been 
used in the local recoil scheme. This option is only used if wimpy showers are 
enabled. 
   
 
<p/> 
The global-recoil machinery does not work well with rescattering in the 
MPI machinery, since then the recoiling system is not uniquely defined. 
<code>MultipartonInteractions:allowRescatter = off</code> by default, 
so this is not a main issue. If both options are switched on, 
rescattering will only be allowed to kick in after the global recoil 
has ceased to be active, i.e. once the <code>nMaxGlobalRecoil</code> 
limit has been exceeded. This should not be a major conflict, 
since rescattering is mainly of interest at later stages of the 
downwards <i>pT</i> evolution. 
 
<p/> 
Further, it is strongly recommended to set 
<code>TimeShower:MEcorrections = off</code> (not default!), i.e. not 
to correct the emission probability to the internal matrix elements. 
The internal ME options do not cover any cases relevant for a multibody 
recoiler anyway, so no guarantees are given what prescription would 
come to be used. Instead, without ME corrections,  a process-independent 
emission rate is obtained, and  <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>user hooks</a> 
can provide the desired process-specific rejection factors. 
 
<h3>Radiation off octet onium states</h3> 
 
In the current implementation, charmonium and bottomonium production 
can proceed either through colour singlet or colour octet mechanisms, 
both of them implemented in terms of <i>2 &rarr; 2</i> hard processes 
such as <i>g g &rarr; (onium) g</i>. 
In the former case the state does not radiate and the onium therefore 
is produced in isolation, up to normal underlying-event activity. In 
the latter case the situation is not so clear, but it is sensible to 
assume that a shower can evolve. (Assuming, of course, that the 
transverse momentum of the onium state is sufficiently high that 
radiation is of relevance.) 
 
<p/> 
There could be two parts to such a shower. Firstly a gluon (or even a 
quark, though less likely) produced in a hard <i>2 &rarr; 2</i> process 
can undergo showering into many gluons, whereof one branches into the 
heavy-quark pair. Secondly, once the pair has been produced, each quark 
can radiate further gluons. This latter kind of emission could easily 
break up a semibound quark pair, but might also create a new semibound 
state where before an unbound pair existed, and to some approximation 
these two effects should balance in the onium production rate. 
The showering "off an onium state" as implemented here therefore should 
not be viewed as an accurate description of the emission history 
step by step, but rather as an effective approach to ensure that the 
octet onium produced "in the hard process" is embedded in a realistic 
amount of jet activity. 
Of course both the isolated singlet and embedded octet are likely to 
be extremes, but hopefully the mix of the two will strike a reasonable 
balance. However, it is possible that some part of the octet production 
occurs in channels where it should not be accompanied by (hard) radiation. 
Therefore reducing the fraction of octet onium states allowed to radiate 
is a valid variation to explore uncertainties. 
 
<p/> 
If an octet onium state is chosen to radiate, the simulation of branchings 
is based on the assumption that the full radiation is provided by an 
incoherent sum of radiation off the quark and off the antiquark of the 
onium state. Thus the splitting kernel is taken to be the normal 
<i>q &rarr; q g</i> one, multiplied by a factor of two. Obviously this is 
a simplification of a more complex picture, averaging over factors pulling 
in different directions. Firstly, radiation off a gluon ought 
to be enhanced by a factor 9/4 relative to a quark rather than the 2 
now used, but this is a minor difference. Secondly, our use of the 
<i>q &rarr; q g</i> branching kernel is roughly equivalent to always 
following the harder gluon in a <i>g &rarr; g g</i> branching. This could 
give us a bias towards producing too hard onia. A soft gluon would have 
little phase space to branch into a heavy-quark pair however, so the 
bias may not be as big as it would seem at first glance. Thirdly, 
once the gluon has branched into a quark pair, each quark carries roughly 
only half of the onium energy. The maximum energy per emitted gluon should 
then be roughly half the onium energy rather than the full, as it is now. 
Thereby the energy of radiated gluons is exaggerated, i.e. onia become too 
soft. So the second and the third points tend to cancel each other. 
 
<p/> 
Finally, note that the lower cutoff scale of the shower evolution depends 
on the onium mass rather than on the quark mass, as it should be. Gluons 
below the octet-onium scale should only be part of the octet-to-singlet 
transition. 
 
<br/><br/><table><tr><td><strong>TimeShower:octetOniumFraction </td><td></td><td> <input type="text" name="26" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
Allow colour-octet charmonium and bottomonium states to radiate gluons. 
0 means that no octet-onium states radiate, 1 that all do, with possibility 
to interpolate between these two extremes. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:octetOniumColFac </td><td></td><td> <input type="text" name="27" value="2." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 4.</code>)</td></tr></table>
The colour factor used used in the splitting kernel for those octet onium 
states that are allowed to radiate, normalized to the <i>q &rarr; q g</i> 
splitting kernel. Thus the default corresponds to twice the radiation 
off a quark. The physically preferred range would be between 1 and 9/4. 
   
 
<h3>Weak showers</h3> 
 
The emission of weak gauge bosons is an integrated part of the initial- 
and final-state radiation, see <?php $filepath = $_GET["filepath"];
echo "<a href='WeakShowers.php?filepath=".$filepath."' target='page'>";?>Weak Showers</a>. 
The following settings are those specifically related to the final-state 
weak radiation, while common settings are found in the 
<?php $filepath = $_GET["filepath"];
echo "<a href='WeakShowers.php?filepath=".$filepath."' target='page'>";?>Weak Showers</a> description. 
 
<br/><br/><strong>TimeShower:weakShower</strong>  <input type="radio" name="28" value="on"><strong>On</strong>
<input type="radio" name="28" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow a weak shower, yes or no. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:weakShowerMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Determine which branchings are allowed. 
<br/>
<input type="radio" name="29" value="0" checked="checked"><strong>0 </strong>:  both <ei>W^+-</ei> and <ei>Z^0</ei> branchings.  <br/>
<input type="radio" name="29" value="1"><strong>1 </strong>:  only <ei>W^+-</ei> branchings. <br/>
<input type="radio" name="29" value="2"><strong>2 </strong>:  only <ei>Z^0</ei> branchings. <br/>
 
<br/><br/><table><tr><td><strong> </td><td></td><td> <input type="text" name="30" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 2.0</code>)</td></tr></table>
Parton shower cut-off <i>pT</i> for weak branchings. 
   
 
<h3>Further variables</h3> 
 
There are several possibilities you can use to switch on or off selected 
branching types in the shower, or in other respects simplify the shower. 
These should normally not be touched. Their main function is for 
cross-checks. 
 
<br/><br/><strong>TimeShower:QCDshower</strong>  <input type="radio" name="31" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="31" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow a QCD shower, i.e. branchings <i>q &rarr; q g</i>, 
<i>g &rarr; g g</i> and <i>g &rarr; q qbar</i>; on/off = true/false. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:nGluonToQuark  </td><td></td><td> <input type="text" name="32" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed quark flavours in <i>g &rarr; q qbar</i> branchings 
(phase space permitting). A change to 4 would exclude 
<i>g &rarr; b bbar</i>, etc. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:weightGluonToQuark  </td><td>  &nbsp;&nbsp;(<code>default = <strong>4</strong></code>; <code>minimum = 1</code>; <code>maximum = 8</code>)</td></tr></table>
Different options to assign kinematics distributions and weights 
for <ei>g &rarr; q qbar</ei> branchings, notably for charm and bottom 
quarks. These options also have the corresponding effect on 
<ei>gamma &rarr; f fbar</ei> branchings. The rationale for the options 
is described in <a href="../pdfdoc/g2qqbarsplit.pdf">this note</a>. 
<br/>Notation: <ei>r_q = m_q^2/m_qq^2</ei>, <ei>beta = sqrt(1 - 4r_q)</ei>, 
with <ei>m_q</ei> the quark mass and <ei>m_qq</ei> the <ei>q qbar</ei> pair 
invariant mass. The scale factor <ei>k</ei> is described below, 
<code>TimeShower:scaleGluonToQuark</code>. 
<br/>
<input type="radio" name="33" value="1"><strong>1 </strong>: same splitting kernel <ei>(1/2) (z^2 + (1-z)^2)</ei> for  massive as massless quarks, only with an extra <ei>beta</ei> phase  space factor.<br/>
<input type="radio" name="33" value="2"><strong>2 </strong>: a splitting kernel  <ei>(beta/2) (z^2 + (1-z)^2 + 8r_q z(1-z))</ei>.<br/>
<input type="radio" name="33" value="3"><strong>3 </strong>: a splitting kernel <ei>z^2 + (1-z)^2 + 8r_q z(1-z)</ei>,  normalized so that the <ei>z</ei>-integrated rate is  <ei>(beta/3) (1 + r/2)</ei>.<br/>
<input type="radio" name="33" value="4" checked="checked"><strong>4 </strong>: same as 3, but additionally a suppression factor  <ei>(1 - m_qq^2/m_dipole^2)^3</ei>, which reduces the rate of high-mass  <ei>q qbar</ei> pairs.<br/>
<input type="radio" name="33" value="5"><strong>5 </strong>: same as 1, but reweighted to an <ei>alpha_s(k m_qq^2)</ei>  rather than the normal <ei>alpha_s(pT^2)</ei>.<br/>
<input type="radio" name="33" value="6"><strong>6 </strong>: same as 2, but reweighted to an <ei>alpha_s(k m_qq^2)</ei>  rather than the normal <ei>alpha_s(pT^2)</ei>.<br/>
<input type="radio" name="33" value="7"><strong>7 </strong>: same as 3, but reweighted to an <ei>alpha_s(k m_qq^2)</ei>  rather than the normal <ei>alpha_s(pT^2)</ei>.<br/>
<input type="radio" name="33" value="8"><strong>8 </strong>: same as 4, but reweighted to an <ei>alpha_s(k m_qq^2)</ei>  rather than the normal <ei>alpha_s(pT^2)</ei>.<br/>
 
<br/><br/><table><tr><td><strong>TimeShower:scaleGluonToQuark </td><td></td><td> <input type="text" name="34" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 1.0</code>)</td></tr></table>
Extra scale parameter <i>k</i> for 
<code>TimeShower:weightGluonToQuark</code> options 5 - 8. Comes on top of 
<code>TimeShower:renormMultFac</code>, which affects <i>alpha_s(pT^2)</i> 
alike. 
   
 
<br/><br/><strong>TimeShower:QEDshowerByQ</strong>  <input type="radio" name="35" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="35" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow quarks to radiate photons, i.e. branchings <i>q &rarr; q gamma</i>; 
on/off = true/false. 
   
 
<br/><br/><strong>TimeShower:QEDshowerByL</strong>  <input type="radio" name="36" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="36" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow leptons to radiate photons, i.e. branchings <i>l &rarr; l gamma</i>; 
on/off = true/false. 
   
 
<br/><br/><strong>TimeShower:QEDshowerByGamma</strong>  <input type="radio" name="37" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="37" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow photons to branch into lepton or quark pairs, i.e. branchings 
<i>gamma &rarr; l+ l-</i> and <i>gamma &rarr; q qbar</i>; 
on/off = true/false. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:nGammaToQuark  </td><td></td><td> <input type="text" name="38" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed quark flavours in <i>gamma &rarr; q qbar</i> branchings 
(phase space permitting). A change to 4 would exclude 
<i>g &rarr; b bbar</i>, etc. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:nGammaToLepton  </td><td></td><td> <input type="text" name="39" value="3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Number of allowed lepton flavours in <i>gamma &rarr; l+ l-</i> branchings 
(phase space permitting). A change to 2 would exclude 
<i>gamma &rarr; tau+ tau-</i>, and a change to 1 also 
<i>gamma &rarr; mu+ mu-</i>. 
   
 
<br/><br/><strong>TimeShower:MEcorrections</strong>  <input type="radio" name="40" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="40" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use of matrix element corrections where available; on/off = true/false. 
   
 
<br/><br/><strong>TimeShower:MEafterFirst</strong>  <input type="radio" name="41" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="41" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use of matrix element corrections also after the first emission, 
for dipole ends of the same system that did not yet radiate. 
Only has a meaning if <code>MEcorrections</code> above is 
switched on. 
   
 
<br/><br/><strong>TimeShower:phiPolAsym</strong>  <input type="radio" name="42" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="42" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Azimuthal asymmetry induced by gluon polarization; on/off = true/false. 
   
 
<br/><br/><strong>TimeShower:recoilToColoured</strong>  <input type="radio" name="43" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="43" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
In the decays of coloured resonances, say <i>t &rarr; b W</i>, it is not 
possible to set up dipoles with matched colours. Originally the 
<i>b</i> radiator therefore has <i>W</i> as recoiler, and that 
choice is unique. Once a gluon has been radiated, however, it is 
possible either to have the unmatched colour (inherited by the gluon) 
still recoiling against the <i>W</i> (<code>off</code>), or else 
let it recoil against the <i>b</i> also for this dipole 
(<code>on</code>). Before version 8.160 the former was the only 
possibility, which could give unphysical radiation patterns. It is 
kept as an option to check backwards compatibility. The same issue 
exists for QED radiation, but obviously is less significant. Consider 
the example <i>W &rarr; e nu</i>, where originally the <i>nu</i> 
takes the recoil. In the old (<code>off</code>) scheme the <i>nu</i> 
would remain recoiler, while in the new (<code>on</code>) instead 
each newly emitted photon becomes the new recoiler. 
   
 
<br/><br/><strong>TimeShower:useFixedFacScale</strong>  <input type="radio" name="44" value="on"><strong>On</strong>
<input type="radio" name="44" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow the possibility to use a fixed factorization scale, set by 
the <code>parm</code> below. This option is unphysical and only 
intended for toy-model and debug studies. 
   
 
<br/><br/><table><tr><td><strong>TimeShower:fixedFacScale </td><td></td><td> <input type="text" name="45" value="100." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100.</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
The fixed factorization scale, in GeV, that would be used in the 
evaluation of parton densities if the <code>flag</code> above is on. 
   
 
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

if($_POST["1"] != "1")
{
$data = "TimeShower:pTmaxMatch = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "1.0")
{
$data = "TimeShower:pTmaxFudge = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1.0")
{
$data = "TimeShower:pTmaxFudgeMPI = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0")
{
$data = "TimeShower:pTdampMatch = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.0")
{
$data = "TimeShower:pTdampFudge = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.1365")
{
$data = "TimeShower:alphaSvalue = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1")
{
$data = "TimeShower:alphaSorder = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "false")
{
$data = "TimeShower:alphaSuseCMW = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1")
{
$data = "TimeShower:alphaEMorder = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "1.")
{
$data = "TimeShower:renormMultFac = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1.")
{
$data = "TimeShower:factorMultFac = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.5")
{
$data = "TimeShower:pTmin = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.5")
{
$data = "TimeShower:pTminChgQ = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "1e-6")
{
$data = "TimeShower:pTminChgL = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "10.0")
{
$data = "TimeShower:mMaxGamma = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "on")
{
$data = "TimeShower:interleave = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "on")
{
$data = "TimeShower:allowBeamRecoil = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "on")
{
$data = "TimeShower:dampenBeamRecoil = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "TimeShower:allowMPIdipole = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "TimeShower:globalRecoil = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "2")
{
$data = "TimeShower:nMaxGlobalRecoil = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0")
{
$data = "TimeShower:globalRecoilMode = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "-1")
{
$data = "TimeShower:nMaxGlobalBranch = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "-1")
{
$data = "TimeShower:nPartonsInBorn = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "off")
{
$data = "TimeShower:limitPTmaxGlobal = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "1.")
{
$data = "TimeShower:octetOniumFraction = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "2.")
{
$data = "TimeShower:octetOniumColFac = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off")
{
$data = "TimeShower:weakShower = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "0")
{
$data = "TimeShower:weakShowerMode = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "1.0")
{
$data = " = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "on")
{
$data = "TimeShower:QCDshower = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "5")
{
$data = "TimeShower:nGluonToQuark = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "4")
{
$data = "TimeShower:weightGluonToQuark = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "1.0")
{
$data = "TimeShower:scaleGluonToQuark = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "on")
{
$data = "TimeShower:QEDshowerByQ = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "on")
{
$data = "TimeShower:QEDshowerByL = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "on")
{
$data = "TimeShower:QEDshowerByGamma = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "5")
{
$data = "TimeShower:nGammaToQuark = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "3")
{
$data = "TimeShower:nGammaToLepton = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "on")
{
$data = "TimeShower:MEcorrections = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "on")
{
$data = "TimeShower:MEafterFirst = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "on")
{
$data = "TimeShower:phiPolAsym = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "on")
{
$data = "TimeShower:recoilToColoured = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
if($_POST["44"] != "off")
{
$data = "TimeShower:useFixedFacScale = ".$_POST["44"]."\n";
fwrite($handle,$data);
}
if($_POST["45"] != "100.")
{
$data = "TimeShower:fixedFacScale = ".$_POST["45"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
