*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:36  cernlib
* Geant
*
*
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
>Name GBROWS
 
>Graphics
 
>Browse VOLU 'Volumes data structure' GXOBJ
 List       .  ' '
 Create     .  '-Svol; -Spos; +Editv'
 Position   .  '-Spos; +Editv'
 Divide     .  '-Sdvn; +Editv'
+
 'Save data structures in RZ file'  .  '-rz/fil'
 'Read data structures from RZ file'  .  '-rz/fil'
 
>Class Box 'Shape box volumes' big_Box sm_Box
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Trd1 'Shape trd1 volumes' big_Trd1 sm_Trd1
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Trd2 'Shape trd2 volumes' big_Trd2 sm_Trd2
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Trap 'Shape trap volumes' big_Trap sm_Trap
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Tube 'Shape tube volumes' big_Tube sm_Tube
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Tubs 'Shape tubs volumes' big_Tubs sm_Tubs
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Cone 'Shape cone volumes' big_Cone sm_Cone
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Cons 'Shape cons volumes' big_Cons sm_Cons
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Sphe 'Shape sphe volumes' big_Sphe sm_Sphe
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Para 'Shape para volumes' big_Para sm_Para
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Pgon 'Shape pgon volumes' big_Pgon sm_Pgon
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Pcon 'Shape pcon volumes' big_Pcon sm_Pcon
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Eltu 'Shape eltu volumes' big_Eltu sm_Eltu
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Hype 'Shape hype volumes' big_Hype sm_Hype
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Gtra 'Shape gtra volumes' big_Gtra sm_Gtra
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Ctub 'Shape ctub volumes' big_Ctub sm_Ctub
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class New 'New_shape' big_New sm_New
 Spec       .  '+Dspec [this]'
 Tree       .  'Dtree [this] 3 111'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Print      .  '+Pvolu [that]'
 Satt       .  'Satt [this]'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
+
 Draw       .  'Draw [this]'
 Cbox       .  '-Box'
 Ctub       .  '-Tube'
 Ccon       .  '-Cone'
 Csph       .  '-Sphe'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift'
 Satt       .  'Satt [this]'
 Move       .  '-Move'
 
>Class Pick 'Pick_volum' big_Pick sm_Pick
+
 Print   .  '+Pvolu [this1]'
 Spec    .  'Changewk;option nzfl; next; +Dspec [this]; Resetwk; option zfl1'
 Tree    .  'Changewk;option nzfl; next; Dtree [this] 3 111; Resetwk; option zfl1'
 Satt       .  'Satt [this]'
 Edit       .  '-Editv; +Editv'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
 Draw       .  '-Draw [this]'
 Cbox       .  '-Box [this]'
 Ctub       .  '-Tube [this]'
 Ccon       .  '-Cone [this]'
 Csph       .  '-Sphe [this]'
 Bomb       .  '-Bomb'
 Shif       .  '-Draw/Shift [this]'
 
>Class Tree 'Dtree' big_Tree sm_Tree
+
 Spec       .  'Changewk; option nzfl; next; +Dspec [this]; Resetwk; option zfl1'
 Satt       .  'Satt [this]'
 Spec3d     .  'box [this] 0 1000 0 1000 -1000 1000; +D3dspec [this]; _
-D3dspec [this]'
 Draw       .  '-Draw [this]'
 Edit       .  '-Editv; +Editv'
 Move3d     .  'next; +move3d [this]; -move3d [this]'
 
>Class Arrow 'Levels' big_Arrow sm_Arrow
+
 Tree       .  'Next; Dtree [this1] [this] 111'
 Spec       .  'Changewk; option nzfl; next; +Dspec [this1]; Resetwk; option zfl1'
 Satt       .  'Satt [this1]'
 
>Browse MATE 'Materials data structure' GXOBJ
 List       .  ' '
 Def_mat    .  '-Smate'
 Def_mix    .  '-Smixt'
 
>Class Elem 'Basic materials' big_Elem sm_Elem
 Edit       .  '-Smate [that]'
 Print      .  '+Pmate [that]'
 Plot_x-sec .  '-Drmat [that]'
 
>Class Mixt 'Mixtures and compounds' big_Mixt sm_Mixt
 Edit       .  '-Smixt [that]'
 Print      .  '+Pmate [that]'
 Plot_x-sec .  '-Drmat [that]'
 
>Browse TMED 'Tracking medium parameters' GXOBJ
 List       .  ' '
 Define     .  '-Stmed'
 
>Class Med 'Tracking media' big_Med sm_Med
 Edit_med   .  '-Stmed [that]'
 Ed_cut_mec .  '-Stpar [that]'
 Print      .  '+Ptmed [that]'
 
>Browse PART 'Particles data structure' GXOBJ
 List       .  ' '
 Define     .  '-Spart'
 
>Class Part 'Particles' big_Part sm_Part
 Edit       .  '-Spart [that]'
 Print      .  '+Ppart [that]'
 
>Browse KINE 'Kinematics data structure' GXOBJ
 List       .  ' '
 
>Class Kine 'Tracks' big_Track sm_Track
 Print      .  '+Prkine [that]'
+
 Print      .  '+Prkine [this]'
 
>Browse HITS 'Hits data structure' GXOBJ
 List       .  ' '
 
>Class /Hitset 'Sethit' big_Set sm_Set
 List       .  ' '
 Print      .  '+Phits [this]'
+
 Print      .  '+Phits [this]'
 
>Class Hitdet 'Dethit' big_Det sm_Det
 Print      .  '+Phits * [this]'
+
 Print      .  '+Phits [this1] [this] 0'
 
>Class Hitnum 'Numhit' big_Num sm_Num
+
 Print      .  '+Phits [this2] [this1] [this]'
 
>Browse ROTM 'Rotation matrix' GXOBJ
 List       .  ' '
 Create     .  '-Srotm'
 
>Class Rmatr 'Rotation matrix' big_Rmatr sm_Rmatr
 Edit       .  '-Srotm [that]'
 Print      .  '+Protm [that]'
 
>Browse VIEW 'View banks in memory' GXOBJ
 List       .  ' '
 Open       .  '-Dopen'
 Close      .  '+Dclose'
!Delete     .  '-G/del'
 
>Class VB 'View banks id' big_VB sm_VB
 Show       .  'Dshow [that]'
 Zoom       .  '-Zoom'
 Lens       .  '-Lens'
+
 Zoom       .  '-Zoom'
 Show       .  '+Dshow [that]'
 Lens       .  '-Lens'
 
>Icon_bitmaps
 
