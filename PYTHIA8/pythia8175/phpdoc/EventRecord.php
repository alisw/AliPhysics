<html>
<head>
<title>The Event Record</title>
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

<form method='post' action='EventRecord.php'>

<h2>The Event Record</h2>

A <code>Pythia</code> instance contains two members of the  
<code>Event</code> class. The one called <code>process</code> provides
a brief summary of the main steps of the hard process, while the 
one called <code>event</code> contains the full history. The
user would normally interact mainly with the second one, so 
we will exemplify primarily with that one.

<p/>
The <code>Event</code> class to first approximation is a vector of 
<code>Particle</code>s, so that it can expand to fit the current 
event size. The index operator is overloaded, so that e.g.
<code>event[i]</code> corresponds to the <i>i</i>'th particle 
of the object <code>event</code>. Thus <code>event[i].id()</code> 
returns the identity of the <i>i</i>'th particle, and so on. 
Therefore the methods of the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>Particle</a></code> class 
are at least as essential as those of the <code>Event</code> class 
itself. 

<p/>
As used inside PYTHIA, some conventions are imposed on the structure
of the event record. Entry 0 of the <code>vector&lt;Particle&gt;</code> 
is used to represent the event as a whole, with its total four-momentum 
and invariant mass, but does not form part of the event history. 
Lines 1 and 2 contains the two incoming beams, and only from here on 
history tracing works as could be expected. That way unassigned mother 
and daughter indices can be put 0 without ambiguity. Depending on the 
task at hand, a loop may therefore start at index 1 rather than 0 
without any loss. Specifically, for translation to other event record 
formats such as HepMC [<a href="Bibliography.php" target="page">Dob01</a>], where the first index is 1, the 
Pythia entry 0 definitely ought to be skipped in order to minimize the 
danger of indexing errors. 

<p/>
In the following we will list the methods available.
Only a few of them have a function to fill in normal user code.

<h3>Basic output methods</h3>

Some methods are available to read out information on the 
current event record:

<a name="method1"></a>
<p/><strong>Particle& Event::operator[](int i) &nbsp;</strong> <br/>
  
<strong>const Particle& Event::operator[](int i) &nbsp;</strong> <br/>
  
<strong>Particle& Event::at(int i) &nbsp;</strong> <br/>
returns a (<code>const</code>) reference to the <i>i</i>'th particle
in the event record, which can be used to get (or set) all the 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>properties</a> of this particle. 
  

<a name="method2"></a>
<p/><strong>int Event::size() &nbsp;</strong> <br/>
The event size, i.e. the size of the <code>vector&lt;Particle&gt;</code>.
Thus valid particles, to be accessed by the above indexing operator, 
are stored in the range <i>0 &lt;= i &lt; size()</i>. See comment 
above about the (ir)relevance of entry 0. 
  

<a name="method3"></a>
<p/><strong>void Event::list() &nbsp;</strong> <br/>
  
<strong>void Event::list(ostream& os) &nbsp;</strong> <br/>
  
<strong>void Event::list(bool showScaleAndVertex, bool showMothersAndDaughters = false) &nbsp;</strong> <br/>
  
<strong>void Event::list(bool showScaleAndVertex, bool showMothersAndDaughters, ostream& os) &nbsp;</strong> <br/>
Provide a listing of the whole event, i.e. of the 
<code>vector&lt;Particle&gt;</code>. The methods with fewer arguments 
call the final one with the respective default values, and are 
non-inlined so they can be used in a debugger. The basic identity 
code, status, mother, daughter, colour, four-momentum and mass data 
are always given, but the methods can also be called with a few 
optional arguments for further information:
<br/><code>argument</code><strong> showScaleAndVertex </strong> (<code>default = <strong>false</strong></code>) :  optionally give a 
second line for each particle, with the production scale (in GeV), 
the particle polarization (dimensionless), the production vertex 
(in mm or mm/c) and the invariant lifetime (also in mm/c).
  
<br/><code>argument</code><strong> showMothersAndDaughters </strong> (<code>default = <strong>false</strong></code>) : 
gives a list of all daughters and mothers of a particle, as defined by 
the <code>motherList(i)</code> and <code>daughterList(i)</code> methods 
described below. It is mainly intended for debug purposes. 
  
<br/><code>argument</code><strong> os </strong> (<code>default = <strong>cout</strong></code>) :  a reference to the <code>ostream</code>
object to which the event listing will be directed.
  

  

<p/>
Each <code>Particle</code> has two mother and two daughter indices.
These may be used to encode zero, one, two or more mothers/daughters,
depending on the combination of values and status code, according to 
well-defined <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>rules</a>. The 
two methods below can do this job easier for you.

<a name="method4"></a>
<p/><strong>vector&lt;int&gt; Event::motherList(int i) &nbsp;</strong> <br/>
returns a vector of all the mother indices of the particle at index 
<i>i</i>. This list is empty for entries 0, 1 and 2,
i.e. the "system" in line 0 is not counted as part of the history. 
Normally the list contains one or two mothers, but it can also be more, 
e.g. in string fragmentation the whole fragmenting system is counted 
as mothers to the primary hadrons. Many particles may have the same
<code>motherList</code>. Mothers are listed in ascending order.
  

<a name="method5"></a>
<p/><strong>vector&lt;int&gt; Event::daughterList(int i) &nbsp;</strong> <br/>
returns a vector of all the daughter indices of the particle at index 
<i>i</i>. This list is empty for a particle that did 
not decay (or, if the evolution is stopped early enough, a parton
that did not branch), while otherwise it can contain a list of 
varying length, from one to many. For the two incoming beam particles, 
all shower initiators and beam remnants are counted as daughters, 
with the one in slot 0 being the one leading up to the hardest 
interaction. The "system" in line 0 does not have any daughters, 
i.e. is not counted as part of the history. Many partons may have the 
same <code>daughterList</code>. Daughters are listed in ascending order.
  

<a name="method6"></a>
<p/><strong>int Event::statusHepMC(int i) &nbsp;</strong> <br/>
returns the status code according to the HepMC conventions agreed in
February 2009. This convention does not preserve the full information
provided by the internal PYTHIA status code, as obtained by
<code>Particle::status()</code>, but comes reasonably close. 
The allowed output values are:
<ul>
<li>0 : an empty entry, with no meaningful information and therefore
to be skipped unconditionally (should not occur in PYTHIA);</li> 
<li>1 : a final-state particle, i.e. a particle that is not decayed
further by the generator (may also include unstable particles that 
are to be decayed later, as part of the detector simulation);</li> 
<li>2 : a decayed Standard Model hadron or tau or mu lepton, excepting 
virtual intermediate states thereof (i.e. the particle must undergo
a normal decay, not e.g. a shower branching);</li>
<li>3 : a documentation entry (not used in PYTHIA);</li> 
<li>4 : an incoming beam particle;</li> 
<li>11 - 200 : an intermediate (decayed/branched/...) particle that does 
not fulfill the criteria of status code 2, with a generator-dependent 
classification of its nature; in PYTHIA the absolute value of the normal 
status code is used.</li> 
</ul>

  

<h3>Further output methods</h3>

The above methods are the main ones that a normal user would make
frequent use of. There are some further methods that also could come
in handy, in the exploration of the history of an event, but where 
the outcome is not always obvious if one is not familiar with the
detailed structure of an event record.

<a name="method7"></a>
<p/><strong>int Event::iTopCopy(int i) &nbsp;</strong> <br/>
  
<strong>int Event::iBotCopy(int i) &nbsp;</strong> <br/>
are used to trace carbon copies of the particle at index <i>i</i> up 
to its top mother or down to its bottom daughter. If there are no such 
carbon copies, <i>i</i> itself will be returned. A carbon copy is 
when the "same" particle appears several times in the event record, but 
with changed momentum owing to recoil effects. 
  

<a name="method8"></a>
<p/><strong>int Event::iTopCopyId(int i) &nbsp;</strong> <br/>
  
<strong>int Event::iBotCopyId(int i) &nbsp;</strong> <br/>
also trace top mother and bottom daughter, but do not require carbon 
copies, only that one can find an unbroken chain, of mothers or daughters, 
with the same flavour <code>id</code> code. When it encounters ambiguities,
say a <i>g -> g g</i> branching or a <i>u u -> u u</i> hard scattering,
it will stop the tracing and return the current position. It can be confused
by nontrivial flavour changes, e.g. a hard process <i>u d -> d u</i>  
by <i>W^+-</i> exchange will give the wrong answer. These methods
therefore are of limited use for common particles, in particular for the
gluon, but should work well for "rare" particles. 
  

<a name="method9"></a>
<p/><strong>vector&lt;int&gt; Event::sisterList(int i) &nbsp;</strong> <br/>
returns a vector of all the sister indices of the particle at index 
<i>i</i>, i.e. all the daughters of the first mother, except the 
particle itself. 
  

<a name="method10"></a>
<p/><strong>vector&lt;int&gt; Event::sisterListTopBot(int i,bool widenSearch = true) &nbsp;</strong> <br/>
returns a vector of all the sister indices of the particle at index 
<i>i</i>, tracking up and back down through carbon copies 
if required. That is, the particle is first traced up with 
<code>iTopCopy()</code> before its mother is found, and then all 
the particles in the <code>daughterList()</code> of this mother are 
traced down with <code>iBotCopy()</code>, omitting the original 
particle itself. Any non-final particles are removed from the list.
Should this make the list empty the search criterion is widened so that
all final daughters are allowed, not only carbon-copy ones. A second
argument <code>false</code> inhibits the second step, and increases 
the risk that an empty list is returned. A typical example of this
is for ISR cascades, e.g. <i>e -> e gamma</i> where the photon 
may not have any obvious sister in the final state if the bottom copy 
of the photon is an electron that annihilates and thus is not part of 
the final state.  
  

<a name="method11"></a>
<p/><strong>bool Event::isAncestor(int i, int iAncestor) &nbsp;</strong> <br/>
traces the particle <i>i</i> upwards through mother, grandmother, 
and so on, until either <i>iAncestor</i> is found or the top of 
the record is reached. Normally one unique mother is required,
as is the case e.g. in decay chains or in parton showers, so that
e.g. the tracing through a hard scattering would not work. For
hadronization, first-rank hadrons are identified with the respective 
string endpoint quark, which may be useful e.g. for <i>b</i> physics, 
while higher-rank hadrons give <code>false</code>. Currently also 
ministrings that collapsed to one single hadron and junction topologies 
give <code>false</code>.  
  

<p/>
One data member in an <code>Event</code> object is used to keep track 
of the largest <code>col()</code> or <code>acol()</code> colour tag set 
so far, so that new ones do not clash. 

<br/><br/><table><tr><td><strong>Event:startColTag  </td><td></td><td> <input type="text" name="1" value="100" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100</strong></code>; <code>minimum = 0</code>; <code>maximum = 1000</code>)</td></tr></table>
This sets the initial colour tag value used, so that the first one 
assigned is <code>startColTag + 1</code>, etc. The Les Houches accord 
[<a href="Bibliography.php" target="page">Boo01</a>] suggests this number to be 500, but 100 works equally 
well.
  

<a name="method12"></a>
<p/><strong>void Event::initColTag(int colTag = 0) &nbsp;</strong> <br/>
forces the current colour tag value to be the larger of the input
<code>colTag</code> and the above <code>Event:startColTag</code>
values. 
  

<a name="method13"></a>
<p/><strong>int Event::lastColTag() &nbsp;</strong> <br/>
returns the current maximum colour tag.
  

<a name="method14"></a>
<p/><strong>int Event::nextColTag() &nbsp;</strong> <br/>
increases the current maximum colour tag by one and returns this 
new value. This method is used whenever a new colour tag is needed.
   

<p/>
Many event properties are accessible via the <code>Info</code> class, 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>see here</a>. Since they are used 
directly in the event generation, a few are stored directly in the
<code>Event</code> class, however. 

<a name="method15"></a>
<p/><strong>void Event::scale( double scaleIn) &nbsp;</strong> <br/>
  
<strong>double Event::scale() &nbsp;</strong> <br/>
set or get the scale (in GeV) of the hardest process in the event.
Matches the function of the <code>scale</code> variable in the
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a>. 
  

<a name="method16"></a>
<p/><strong>void Event::scaleSecond( double scaleSecondIn) &nbsp;</strong> <br/>
  
<strong>double Event::scaleSecond() &nbsp;</strong> <br/>
set or get the scale (in GeV) of a second hard process in the event,
in those cases where such a one 
<?php $filepath = $_GET["filepath"];
echo "<a href='SecondHardProcess.php?filepath=".$filepath."' target='page'>";?>has been requested</a>.
  

<h3>Constructors and modifications of the event record</h3>
 
Although you would not normally need to create your own 
<code>Event</code> instance, there may be times where that
could be convenient. The typical example would be if you want to
create a new event record that is the sum of a few different ones,
e.g. if you want to simulate pileup events. There may also be cases
where you want to add one or a few particles to an existing event
record.  

<a name="method17"></a>
<p/><strong>Event::Event(int capacity = 100) &nbsp;</strong> <br/>
creates an empty event record, but with a reserved size
<i>capacity</i> for the <code>Particle</code> vector.  
  

<a name="method18"></a>
<p/><strong>Event& Event::operator=(const Event& oldEvent) &nbsp;</strong> <br/>
copies the input event record.
  

<a name="method19"></a>
<p/><strong>Event& Event::operator+=(const Event& addEvent) &nbsp;</strong> <br/>
appends an event to an existing one. For the appended particles 
mother, daughter and colour tags are shifted to make a consistent 
record. The zeroth particle of the appended event is not copied, 
but the zeroth particle of the combined event is updated to the 
full energy-momentum content.
  

<a name="method20"></a>
<p/><strong>void Event::init(string headerIn = &quot;&quot;, ParticleData* particleDataPtrIn = 0, int startColTagIn = 100) &nbsp;</strong> <br/>
initializes colour, the pointer to the particle database, and the 
header specification used for the event listing. We remind that a 
<code>Pythia</code> object contains two event records 
<code>process</code> and <code>event</code>. Thus one may e.g. 
call either  <code>pythia.process.list()</code> or 
<code>pythia.event.list()</code>. To distinguish those two rapidly 
at visual inspection, the <code>"Pythia Event Listing"</code> header 
is printed out differently, in one case adding 
<code>"(hard process)"</code> and in the other 
<code>"(complete event)"</code>. When <code>+=</code> is used to 
append an event, the modified event is printed with 
<code>"(combination of several events)"</code> as a reminder.
  

<a name="method21"></a>
<p/><strong>void Event::clear() &nbsp;</strong> <br/>
empties event record. Specifically the <code>Particle</code> vector 
size is reset to zero.
  

<a name="method22"></a>
<p/><strong>void Event::reset() &nbsp;</strong> <br/>
empties the event record, as <code>clear()</code> above, but then 
fills the zero entry of the <code>Particle</code> vector with the 
pseudoparticle used to represent the event as a whole. At this point
the pseudoparticle is not assigned any momentum or mass.
  

<a name="method23"></a>
<p/><strong>void Event::popBack(int n = 1) &nbsp;</strong> <br/>
removes the last <i>n</i> particle entries; must be a positive 
number.
  

<a name="method24"></a>
<p/><strong>bool Event::undoDecay(int i) &nbsp;</strong> <br/>
removes the decay chain of the particle <i>i</i> and thus restores 
it to its undecayed state. It is only intended for "normal" particle
decay chains, and will return false in other cases, notably if 
the particle is coloured. The procedure would not work if non-local
momentum shifts have been performed, such as with a Bose-Einstein 
shift procedure (or for a dipole shower recoiler).
  

<a name="method25"></a>
<p/><strong>int Event::append(Particle entryIn) &nbsp;</strong> <br/>
appends a particle to the bottom of the event record and 
returns the index of this position. 
  

<a name="method26"></a>
<p/><strong>int Event::append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz,  double e, double m = 0., double scale = 0., double pol = 9.) &nbsp;</strong> <br/>
appends a particle to the bottom of the event record and 
returns the index of this position; 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>see here</a> for the meaning
of the various particle properties.
  

<a name="method27"></a>
<p/><strong>int Event::append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Vec4 p, double m = 0., double scale = 0., double pol = 9.) &nbsp;</strong> <br/>
appends a particle to the bottom of the event record and 
returns the index of this position, as above but with four-momentum
as a <code>Vec4</code>.
  

<a name="method28"></a>
<p/><strong>int Event::append(int id, int status, int col, int acol, double px, double py, double pz, double e, double m = 0., double scale = 0., double pol = 9.) &nbsp;</strong> <br/>
  
<strong>int Event::append(int id, int status, int col, int acol, Vec4 p, double m = 0., double scale = 0., double pol = 9.) &nbsp;</strong> <br/>
appends a particle to the bottom of the event record and 
returns the index of this position, as above but with vanishing
(i.e. zero) mother and daughter indices.
  

<a name="method29"></a>
<p/><strong>int Event::setPDTPtr(int iSet = -1) &nbsp;</strong> <br/>
send in a pointer to the <code>ParticleData</code> database for 
particle <code>iSet</code>, by default the most recently appended 
particle. Also generates a pointer to the 
<code>ParticleDataEntry</code> object of the identity code
of the particle.
  

<a name="method30"></a>
<p/><strong>int Event::copy(int iCopy, int newStatus = 0) &nbsp;</strong> <br/>
copies the existing particle in entry <code>iCopy</code> to the
bottom of the event record and returns the index of this position.
By default, i.e. with <code>newStatus = 0</code>, everything is 
copied precisely as it is, which means that history information 
has to be modified further by hand to make sense. With a positive 
<code>newStatus</code>, the new copy is set up to be the daughter of 
the old, with status code <code>newStatus</code>, while the status 
code of <code>iCopy</code> is negated. With a negative 
<code>newStatus</code>, the new copy is instead set up to be the 
mother of <code>iCopy</code>. An attempt to copy an out-of-range 
entry will return -1.
  

<a name="method31"></a>
<p/><strong>Particle& Event::back() &nbsp;</strong> <br/>
returns a reference to the last particle in the event record.
  

<a name="method32"></a>
<p/><strong>void Event::restorePtrs() &nbsp;</strong> <br/>
each particle in the event record has a pointer to the whole database
and another to the particle species itself, used to find some particle
properties. The latter pointer is automatically set/changed whenever 
the particle identity is set/changed by one of the normal methods. 
(It is the "changed" part that prompts the inclusion of a pointer 
to the whole database.) Of course the pointer values are specific to 
the memory locations of the current run, and so it has no sense to 
save them if events are written to file. Should you use some
persistency scheme that bypasses the normal methods when the event is 
read back in, you can use <code>restorePtrs()</code> afterwards to set 
these pointers appropriately.
  

<p/>
A few methods exist to rotate and boost events. These derive from the
<?php $filepath = $_GET["filepath"];
echo "<a href='FourVectors.php?filepath=".$filepath."' target='page'>";?>Vec4</a> methods, and affect both the 
momentum and the vertex (position) components of all particles. 

<a name="method33"></a>
<p/><strong>void Event::rot(double theta, double phi) &nbsp;</strong> <br/>
rotate all particles in the event by this polar and azimuthal angle 
(expressed in radians). 
  

<a name="method34"></a>
<p/><strong>void Event::bst(double betaX, double betaY, double betaZ) &nbsp;</strong> <br/>
  
<strong>void Event::bst(double betaX, double betaY, double betaZ, double gamma) &nbsp;</strong> <br/>
  
<strong>void Event::bst(const Vec4& vec) &nbsp;</strong> <br/>
boost all particles in the event by this three-vector. 
Optionally you may provide the <i>gamma</i> value as a fourth argument, 
which may help avoid roundoff errors for big boosts. You may alternatively 
supply a <code>Vec4</code> four-vector, in which case the boost vector
becomes <i>beta = p/E</i>.
  

<a name="method35"></a>
<p/><strong>void Event::rotbst(const RotBstMatrix& M) &nbsp;</strong> <br/>
rotate and boost by the combined action encoded in the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='FourVectors.php?filepath=".$filepath."' target='page'>";?>RotBstMatrix</a> M</code>.
  

<h3>The Junction Class</h3>

The event record also contains a vector of junctions, which often
is empty or else contains only a very few per event. Methods are
available to add further junctions or query the current junction list.
This is only for the expert user, however, and is not discussed
further here, but only the main points.  

<p/>
A junction stores the properties associated with a baryon number that
is fully resolved, i.e. where three different colour indices are 
involved. There are two main applications,
<ol> 
<li>baryon beams, where at least two valence quarks are kicked out,
and so the motion of the baryon number is nontrivial;</li>
<li>baryon-number violating processes, e.g. in SUSY with broken
<i>R</i>-parity.</li>
</ol>
Information on junctions is set, partly in the process generation,
partly in the beam remnants machinery, and used by the fragmentation 
routines, but the normal user does not have to know the details. 

<p/>
For each junction, information is stored on the kind of junction, and 
on the three (anti)colour indices that are involved in the junction.
The possibilities foreseen are:
<ul>
<li><code>kind = 1</code> : incoming colourless particle to three 
outgoing colours (e.g. baryon beam remnant or 
<i>neutralino -> q q q</i>);</li>
<li><code>kind = 2</code> : incoming colourless particle to three 
outgoing anticolours;</li>
<li><code>kind = 3</code> : one incoming anticolour (stored first) 
and two outgoing  colours (e.g. antisquark decaying to two quarks, or
  gluino decay to three quarks);</li>
<li><code>kind = 4</code> : one incoming colour (stored first) and two 
outgoing anticolours (e.g. squark decaying to two antiquarks, or
  gluino decaying to three antiquarks);</li>
<li><code>kind = 5</code> : two incoming anticolours (stored first) 
and one outgoing colour (e.g. resonant squark production through RPV);</li>
<li><code>kind = 6</code> : two incoming colours (stored first) 
and one outgoing anticolour (e.g. resonant antisquark production
  through RPV);
</li>
</ul>
The odd (even) <code>kind</code> codes corresponds to a +1 (-1) change in 
baryon number across the junction.

<p/>
The kind and colour information in the list of junctions can be set 
or read with methods of the <code>Event</code> class, but are not of 
common interest and so not described here.

<p/>
A listing of current junctions can be obtained with the 
<code>listJunctions()</code> method.
 
<h3>Subsystems</h3>

Separate from the event record as such, but closely tied to it is the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='AdvancedUsage.php?filepath=".$filepath."' target='page'>";?>PartonSystems</a></code> class, 
which mainly stores the parton indices of incoming and outgoing partons, 
classified by collision subsystem. Such information is needed to 
interleave multiparton interactions, initial-state showers and final-state 
showers, and append beam remnants. It could also be used in other places. 
It is intended to be accessed only by experts, such as implementors of 
<?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>new showering models</a>.

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

if($_POST["1"] != "100")
{
$data = "Event:startColTag = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
