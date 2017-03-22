<html>
<head>
<title>A Second Hard Process</title>
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

<form method='post' action='ASecondHardProcess.php'>
 
<h2>A Second Hard Process</h2> 
 
When you have selected a set of hard processes for hadron beams, the 
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>multiparton interactions</a> 
framework can add further interactions to build up a realistic 
underlying event. These further interactions can come from a wide 
variety of processes, and will occasionally be quite hard. They 
do represent a realistic random mix, however, which means one cannot 
predetermine what will happen. Occasionally there may be cases 
where one wants to specify also the second hard interaction rather 
precisely. The options on this page allow you to do precisely that. 
 
<br/><br/><strong>SecondHard:generate</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Generate two hard scatterings in a collision between hadron beams. 
The hardest process can be any combination of internal processes, 
available in the normal <?php $filepath = $_GET["filepath"];
echo "<a href='ProcessSelection.php?filepath=".$filepath."' target='page'>";?>process 
selection</a> machinery, or external input. Here you must further 
specify which set of processes to allow for the second hard one, see 
the following. 
   
 
<h3>Process Selection</h3> 
 
In principle the whole <?php $filepath = $_GET["filepath"];
echo "<a href='ProcessSelection.php?filepath=".$filepath."' target='page'>";?>process 
selection</a> allowed for the first process could be repeated 
for the second one. However, this would probably be overkill. 
Therefore here a more limited set of prepackaged process collections 
are made available, that can then be further combined at will. 
Since the description is almost completely symmetric between the 
first and the second process, you always have the possibility 
to pick one of the two processes according to the complete list 
of possibilities. 
 
<p/> 
Here comes the list of allowed sets of processes, to combine at will: 
 
<br/><br/><strong>SecondHard:TwoJets</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Standard QCD <i>2 &rarr; 2</i> processes involving gluons and 
<i>d, u, s, c, b</i> quarks. 
   
 
<br/><br/><strong>SecondHard:PhotonAndJet</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
A prompt photon recoiling against a quark or gluon jet. 
 
<br/><br/><strong>SecondHard:TwoPhotons</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Two prompt photons recoiling against each other. 
 
<br/><br/><strong>SecondHard:Charmonium</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Production of charmonium via colour singlet and colour octet channels. 
 
<br/><br/><strong>SecondHard:Bottomonium</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Production of bottomonium via colour singlet and colour octet channels. 
 
<br/><br/><strong>SecondHard:SingleGmZ</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar &rarr; gamma^*/Z^0</i>, with full interference 
between the <i>gamma^*</i> and <i>Z^0</i>. 
   
 
<br/><br/><strong>SecondHard:SingleW</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar' &rarr; W^+-</i>. 
   
 
<br/><br/><strong>SecondHard:GmZAndJet</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar &rarr; gamma^*/Z^0 g</i> and 
<i>q g &rarr; gamma^*/Z^0 q</i>. 
   
 
<br/><br/><strong>SecondHard:WAndJet</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar' &rarr; W^+- g</i> and 
<i>q g &rarr; W^+- q'</i>. 
   
 
<br/><br/><strong>SecondHard:TopPair</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Production of a top pair, either via QCD processes or via an 
intermediate <i>gamma^*/Z^0</i> resonance. 
   
 
<br/><br/><strong>SecondHard:SingleTop</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Production of a single top, either via a <i>t-</i> or 
an <i>s-</i>channel <i>W^+-</i> resonance. 
   
 
<p/> 
A further process collection comes with a warning flag: 
 
<br/><br/><strong>SecondHard:TwoBJets</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
The <i>q qbar &rarr; b bbar</i> and <i>g g &rarr; b bbar</i> processes. 
These are already included in the <code>TwoJets</code> sample above, 
so it would be double-counting to include both, but we assume there 
may be cases where the <i>b</i> subsample will be of special interest. 
This subsample does not include flavour-excitation or gluon-splitting 
contributions to the <i>b</i> rate, however, so, depending 
on the topology if interest, it may or may not be a good approximation. 
   
 
<h3>Cuts and scales</h3> 
 
The second hard process obeys exactly the same selection rules for 
<?php $filepath = $_GET["filepath"];
echo "<a href='PhaseSpaceCuts.php?filepath=".$filepath."' target='page'>";?>phase space cuts</a> and 
<?php $filepath = $_GET["filepath"];
echo "<a href='CouplingsAndScales.php?filepath=".$filepath."' target='page'>";?>couplings and scales</a> 
as the first one does. Specifically, a <i>pTmin</i> cut for 
<i>2 &rarr; 2</i> processes would apply to the first and the second hard 
process alike, and ballpark half of the time the second could be 
generated with a larger <i>pT</i> than the first. (Exact numbers 
depending on the relative shape of the two cross sections.) That is, 
first and second is only used as an administrative distinction between 
the two, not as a physics ordering one. 
 
<p/> 
Optionally it is possible to pick the mass and <i>pT</i> 
<?php $filepath = $_GET["filepath"];
echo "<a href='PhaseSpaceCuts.php?filepath=".$filepath."' target='page'>";?>phase space cuts</a> separately for 
the second hard interaction. The main application presumably would 
be to allow a second process that is softer than the first, but still 
hard. But one is also free to make the second process harder than the 
first, if desired. So long as the two <i>pT</i> (or mass) ranges 
overlap the ordering will not be the same in all events, however. 
 
<h3>Cross-section calculation</h3> 
 
As an introduction, a brief reminder of Poissonian statistics. 
Assume a stochastic process in time, for now not necessarily a 
high-energy physics one, where the probability for an event to occur 
at any given time is independent of what happens at other times. 
Then the probability for <i>n</i> events to occur in a finite 
time interval is 
<br/><i> 
P_n = &lt;n&gt;^n exp(-&lt;n&gt;) / n! 
</i><br/> 
where <i>&lt;n&gt;</i> is the average number of events. If this 
number is small we can approximate <i>exp(-&lt;n&gt;) = 1 </i>, 
so that <i>P_1 = &lt;n&gt;</i> and 
<i>P_2 = &lt;n&gt;^2 / 2 = P_1^2 / 2</i>. 
 
<p/> 
Now further assume that the events actually are of two different 
kinds <i>a</i> and <i>b</i>, occurring independently of each 
other, such that <i>&lt;n&gt; = &lt;n_a&gt; + &lt;n_b&gt;</i>. 
It then follows that the probability of having one event of type 
<i>a</i> (or <i>b</i>) and nothing else is 
<i>P_1a = &lt;n_a&gt;</i> (or <i>P_1b = &lt;n_b&gt;</i>). 
From 
<br/><i> 
P_2 = (&lt;n_a&gt; + &lt;n_b&gt)^2 / 2 = (P_1a + P_1b)^2 / 2 = 
(P_1a^2 + 2 P_1a P_1b + P_1b^2) / 2 
</i><br/> 
it is easy to read off that the probability to have exactly two 
events of kind <i>a</i> and none of <i>b</i> is 
<i>P_2aa = P_1a^2 / 2</i> whereas that of having one <i>a</i> 
and one <i>b</i> is <i>P_2ab = P_1a P_1b</i>. Note that the 
former, with two identical events, contains a factor <i>1/2</i> 
while the latter, with two different ones, does not. If viewed 
in a time-ordered sense, the difference is that the latter can be 
obtained two ways, either first an <i>a</i> and then a <i>b</i> 
or else first a <i>b</i> and then an <i>a</i>. 
 
<p/> 
To translate this language into cross-sections for high-energy 
events, we assume that interactions can occur at different <i>pT</i> 
values independently of each other inside inelastic nondiffractive 
(sometimes equated with "minbias") events. Then the above probabilities 
translate into 
<i>P_n = sigma_n / sigma_ND</i> where <i>sigma_ND</i> is the 
total nondiffractive cross section. Again we want to assume that 
<i>exp(-&lt;n&gt;)</i> is close to unity, i.e. that the total 
hard cross section above <i>pTmin</i> is much smaller than 
<i>sigma_ND</i>. The hard cross section is dominated by QCD 
jet production, and a reasonable precaution is to require a 
<i>pTmin</i> of at least 20 GeV at LHC energies. 
(For <i>2 &rarr; 1</i> processes such as 
<i>q qbar &rarr; gamma^*/Z^0 (&rarr; f fbar)</i> one can instead make a 
similar cut on mass.) Then the generic equation 
<i>P_2 = P_1^2 / 2</i> translates into 
<i>sigma_2/sigma_ND = (sigma_1 / sigma_ND)^2 / 2</i> or 
<i>sigma_2 = sigma_1^2 / (2 sigma_ND)</i>. 
 
<p/> 
Again different processes <i>a, b, c, ...</i> contribute, 
and by the same reasoning we obtain 
<i>sigma_2aa = sigma_1a^2 / (2 sigma_ND)</i>, 
<i>sigma_2ab = sigma_1a sigma_1b / sigma_ND</i>, 
and so on. 
 
<p/> 
There is one important correction to this picture: all collisions 
do no occur under equal conditions. Some are more central in impact 
parameter, others more peripheral. This leads to a further element of 
variability: central collisions are likely to have more activity 
than the average, peripheral less. Integrated over impact 
parameter standard cross sections are recovered, but correlations 
are affected by a "trigger bias" effect: if you select for events 
with a hard process you favour events at small impact parameter 
which have above-average activity, and therefore also increased 
chance for further interactions. (In PYTHIA this is the origin 
of the "pedestal effect", i.e. that events with a hard interaction 
have more underlying activity than the level found in minimum-bias 
events.) 
 
<p/> 
When you specify a matter overlap profile in the multiparton-interactions 
scenario, such an enhancement/depletion factor <i>f_impact</i> is 
chosen event-by-event and can be averaged during the course of the run. 
As an example, the double Gaussian form used in Tune A gives 
approximately <i>&lt;f_impact&gt; = 2.5</i>. In general, the more 
uneven the distribution the higher the <i>&lt;f_impact&gt;</i>. 
Also the <i>pT0</i> parameter value has an impact, even if it is 
less important over a realistic range of values, although it implies 
that <i>&lt;f_impact&gt;</i> is energy-dependent. The origin of this 
effect is as follows. A lower <i>pT0</i> implies more MPI activity 
at all impact parameters, so that the nondiffractive cross section 
<i>sigma_ND</i> increases, or equivalently the proton size. But if 
<i>sigma_ND</i> is fixed by data then the input radius of the matter 
overlap profile (not explicitly specified but implicitly adjusted at 
initialization) has to be shrunk  so that the output value can stay 
constant. This means that the proton matter is more closely packed and 
therefore <i>&lt;f_impact&gt;</i> goes up. 
 
<p/> 
The above equations therefore have to be modified to 
<i>sigma_2aa = &lt;f_impact&gt; sigma_1a^2 / (2 sigma_ND)</i>, 
<i>sigma_2ab = &lt;f_impact&gt; sigma_1a sigma_1b / sigma_ND</i>. 
Experimentalists often instead use the notation 
<i>sigma_2ab = sigma_1a sigma_1b / sigma_eff</i>, 
from which we see that PYTHIA "predicts" 
<i>sigma_eff = sigma_ND / &lt;f_impact&gt;</i>. 
When the generation of multiparton interactions is switched off it is 
not possible to calculate <i>&lt;f_impact&gt;</i> and therefore 
it is set to unity. 
 
<p/> 
When this recipe is to be applied to calculate 
actual cross sections, it is useful to distinguish three cases, 
depending on which set of processes are selected to study for 
the first and second interaction. 
 
<p/> 
(1) The processes <i>a</i> for the first interaction and 
<i>b</i> for the second one have no overlap at all. 
For instance, the first could be <code>TwoJets</code> and the 
second <code>TwoPhotons</code>. In that case, the two interactions 
can be selected independently, and cross sections tabulated 
for each separate subprocess in the two above classes. At the 
end of the run, the cross sections in <i>a</i> should be multiplied 
by <i>&lt;f_impact&gt; sigma_1b / sigma_ND</i> to bring them to 
the correct overall level, and those in <i>b</i> by 
<i>&lt;f_impact&gt; sigma_1a / sigma_ND</i>. 
 
<p/> 
(2) Exactly the same processes <i>a</i> are selected for the 
first and second interaction. In that case it works as above, 
with <i>a = b</i>, and it is only necessary to multiply by an 
additional factor <i>1/2</i>. A compensating factor of 2 
is automatically obtained for picking two different subprocesses, 
e.g. if <code>TwoJets</code> is selected for both interactions, 
then the combination of the two subprocesses <i>q qbar &rarr; g g</i> 
and <i>g g &rarr; g g</i> can trivially be obtained two ways. 
 
<p/> 
(3) The list of subprocesses partly but not completely overlap. 
For instance, the first process is allowed to contain <i>a</i> 
or <i>c</i> and the second <i>b</i> or <i>c</i>, where 
there is no overlap between <i>a</i> and <i>b</i>. Then, 
when an independent selection for the first and second interaction 
both pick one of the subprocesses in <i>c</i>, half of those 
events have to be thrown, and the stored cross section reduced 
accordingly. Considering the four possible combinations of first 
and second process, this gives a 
<br/><i> 
sigma'_1 = sigma_1a + sigma_1c * (sigma_2b + sigma_2c/2) / 
(sigma_2b + sigma_2c) 
</i><br/> 
with the factor <i>1/2</i> for the <i>sigma_1c sigma_2c</i> term. 
At the end of the day, this <i>sigma'_1</i> should be multiplied 
by the normalization factor 
<br/><i> 
f_1norm = &lt;f_impact&gt; (sigma_2b + sigma_2c) / sigma_ND 
</i><br/> 
here without a factor <i>1/2</i> (or else it would have been 
double-counted). This gives the correct 
<br/><i> 
(sigma_2b + sigma_2c) * sigma'_1 = sigma_1a * sigma_2b 
+ sigma_1a * sigma_2c + sigma_1c * sigma_2b + sigma_1c * sigma_2c/2 
</i><br/> 
The second interaction can be handled in exact analogy. 
 
<p/> 
For the considerations above it is assumed that the phase space cuts 
are the same for the two processes. It is possible to set the mass and 
transverse momentum cuts differently, however. This changes nothing 
for processes that already are different. For two collisions of the 
same type it is partly a matter of interpretation what is intended. 
If we consider the case of the same process in two non-overlapping 
phase space regions, most likely we want to consider them as 
separate processes, in the sense that we expect a factor 2 relative 
to Poissonian statistics from either of the two hardest processes 
populating either of the two phase space regions. In total we are 
therefore lead to adopt the same strategy as in case (3) above: 
only in the overlapping part of the two allowed phase space regions 
could two processes be identical and thus appear with a 1/2 factor, 
elsewhere the two processes are never identical and do not 
include the 1/2 factor. We reiterate, however, that the case of 
partly but not completely overlapping phase space regions for one and 
the same process is tricky, and not to be used without prior 
deliberation. 
 
<p/> 
The listing obtained with the <code>pythia.stat()</code> 
already contain these corrections factors, i.e. cross sections 
are for the occurrence of two interactions of the specified kinds. 
There is not a full tabulation of the matrix of all the possible 
combinations of a specific first process together with a specific 
second one (but the information is there for the user to do that, 
if desired). Instead <code>pythia.stat()</code> shows this 
matrix projected onto the set of processes and associated cross 
sections for the first and the second interaction, respectively. 
Up to statistical fluctuations, these two sections of the 
<code>pythia.stat()</code> listing both add up to the same 
total cross section for the event sample. 
 
<p/> 
There is a further special feature to be noted for this listing, 
and that is the difference between the number of "selected" events 
and the number of "accepted" ones. Here is how that comes about. 
Originally the first and second process are selected completely 
independently. The generation (in)efficiency is reflected in the 
different number of initially tried events for the first and second 
process, leading to the same number of selected events. While 
acceptable on their own, the combination of the two processes may 
be unacceptable, however. It may be that the two processes added 
together use more energy-momentum than kinematically allowed, or, 
even if not, are disfavoured when the PYTHIA approach to provide 
correlated parton densities is applied. Alternatively, referring 
to case (3) above, it may be because half of the events should 
be thrown for identical processes. Taken together, it is these 
effects that reduced the event number from "selected" to "accepted". 
(A further reduction may occur if a 
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>user hook</a> rejects some events.) 
 
<p/> 
It is allowed to use external Les Houches Accord input for the 
hardest process, and then pick an internal one for the second hardest. 
In this case PYTHIA does not have access to your thinking concerning 
the external process, and cannot know whether it overlaps with the 
internal or not. (External events <i>q qbar' &rarr; e+ nu_e</i> could 
agree with the internal <i>W</i> ones, or be a <i>W'</i> resonance 
in a BSM scenario, to give one example.) Therefore the combined cross 
section is always based on the scenario (1) above. Corrections for 
correlated parton densities are included also in this case, however. 
That is, an external event that takes a large fraction of the incoming 
beam momenta stands a fair chance of being rejected when it has to be 
combined with another hard process. For this reason the "selected" and 
"accepted" event numbers are likely to disagree. 
 
<p/> 
In the cross section calculation above, the <i>sigma'_1</i> 
cross sections are based on the number of accepted events, while 
the <i>f_1norm</i> factor is evaluated based on the cross sections 
for selected events. That way the suppression by correlations 
between the two processes does not get to be double-counted. 
 
<p/> 
The <code>pythia.stat()</code> listing contains two final 
lines, indicating the summed cross sections <i>sigma_1sum</i> and 
<i>sigma_2sum</i> for the first and second set of processes, at 
the "selected" stage above, plus information on the <i>sigma_ND</i> 
and <i>&lt;f_impact&gt;</i> used. The total cross section 
generated is related to this by 
<br/><i> 
&lt;f_impact&gt; * (sigma_1sum * sigma_2sum / sigma_ND) * 
(n_accepted / n_selected) 
</i><br/> 
 with an additional factor of <i>1/2</i> for case 2 above. 
 
<p/> 
The error quoted for the cross section of a process is a combination 
in quadrature of the error on this process alone with the error on 
the normalization factor, including the error on 
<i>&lt;f_impact&gt;</i>. As always it is a purely statistical one 
and of course hides considerably bigger systematic uncertainties. 
 
<br/> 
<br/><b>Warning:</b> the calculational machinery above has not (yet) 
been implemented for the case that the two interactions are to be 
associated with different impact-parameter profiles, as is the case 
for <code>MultipartonInteractions:bProfile = 4</code>, i.e. when the 
radius depends on the <i>x</i> value. Results for the double 
cross section therefore cannot be trusted in this case. 
 
<h3>Event information</h3> 
 
Normally the <code>process</code> event record only contains the 
hardest interaction, but in this case also the second hardest 
is stored there. If both of them are <i>2 &rarr; 2</i> ones, the 
first would be stored in lines 3 - 6 and the second in 7 - 10. 
For both, status codes 21 - 29 would be used, as for a hardest 
process. Any resonance decay chains would occur after the two 
main processes, to allow normal parsing. The beams in 1 and 2 
only appear in one copy. This structure is echoed in the 
full <code>event</code> event record. 
 
<p/> 
Most of the properties accessible by the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>pythia.info</a></code> 
methods refer to the first process, whether that happens to be the 
hardest or not. The code and <i>pT</i> scale of the second process 
are accessible by the <code>info.codeMPI(1)</code> and 
<code>info.pTMPI(1)</code>, however. 
 
<p/> 
The <code>sigmaGen()</code> and <code>sigmaErr()</code> methods provide 
the cross section and its error for the event sample as a whole, 
combining the information from the two hard processes as described 
above. In particular, the former should be used to give the 
weight of the generated event sample. The statistical error estimate 
is somewhat cruder and gives a larger value than the 
subprocess-by-subprocess one  employed in 
<code>pythia.stat()</code>, but this number is 
anyway less relevant, since systematical errors are likely to dominate. 
 
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

if($_POST["1"] != "off")
{
$data = "SecondHard:generate = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "SecondHard:TwoJets = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "SecondHard:PhotonAndJet = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "SecondHard:TwoPhotons = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "SecondHard:Charmonium = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "SecondHard:Bottomonium = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "SecondHard:SingleGmZ = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "SecondHard:SingleW = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "SecondHard:GmZAndJet = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "SecondHard:WAndJet = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "SecondHard:TopPair = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "SecondHard:SingleTop = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "SecondHard:TwoBJets = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
