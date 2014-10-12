##\mainpage  VMD extension functions
#
# A collection of TCL-VMD functions that support a number of
# structural transformations and shortcuts for TCL-VMD programmers.
# Features easy-to-remember semantics to
#
#  * Iterate a block of code over frames
#  * Iterate a block of code over trajectory files
#  * Compute the number and fraction of native contacts
#  * Compute distance matrices
#  * Compute root-mean-square alignments and related measures
#  * ...and more
# 
# Please refer to the <a href="modules.html">table of contents</a> for
# the full feature list.
#


##\page download Get the code
#
# Download the latest release from
# https://github.com/tonigi/vmd_extensions/releases/latest . 
#
# Right now there is no installation procedure - just issue
#
# \code source VMDextensions.tcl
# \endcode

##\page license License
#
#\copyright 2010-2014 CNR and UPF
#\author  toni.giorgino  isib cnr it.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 



##\defgroup forFrames Iterate a TCL block over all frames
# @{
# 
# Convenience functions to iterate a block of TCL code over
# a set of frames. See documentation and examples below. 
#
# **Note.** If the selection is subject to changes between frames,
# do perform a "$sel update" inside the block. If you are using
# more than one selection, pass only the first as an argument and, if
# necessary, update the others inside the block.
#
# Example
# ---------
#
# Iterate and compute secondary structure:
# \code
#     forFrames fn $kid {
#         animate goto $fn
#         mol ssrecalc [$kid molid]
#         puts "$fn [vecmean [$kid get alpha_helix]]" 
#     }
# \endcode

##  Iterators over frames in the current selection.  While in ''block'',
#  the first argument is set to the (integer) current frame
#  number. The selection is used to get the frame range, and is
#  automatically set with "frame" before executing the block.  Note
#  that, like '''for''', the variable name does not require the dollar
#  symbol.  Example:
# 
#        forFrames fno $sel { puts "$fno: [measure rgyr $sel ]"  }
#
proc forFrames {framenum sel block {step 1}} {
    upvar 1 $framenum pwvar
    set n [ molinfo [$sel molid] get numframes ]
    for { set i 0 } { $i < $n } { incr i $step } {
	set pwvar $i
	$sel frame $i
	uplevel $block
    }
}

## Similar to forFrames (in fact it is a plug-in replacement). If the
# block returns a value, also builds a list with the returned values.
proc forFramesMeasure {framenum sel block} {
    upvar 1 $framenum pwvar
    set n [ molinfo [$sel molid] get numframes ]
    set r {}
    for { set i 0 } { $i < $n } { incr i } {
	set pwvar $i
	$sel frame $i
	lappend r [uplevel $block]
    }
    return $r
}

## @}

# ----------------------------------------

##\defgroup forFiles Iterate a TCL block over a set of files
# @{
# Iterator over files matched by the given pattern. The matched files
# are loaded in sequence, sorted in [natural
# order](http://sourcefrog.net/projects/natsort/), and the block is
# executed. When executing the block, the first argument is set to the
# current file name (note that, like *for*, it does not require the
# dollar symbol). Will discard currently-loaded frames! 
#
# For example:
# \code
#       forFiles fn {*.dcd} {
#             puts "$fn has [ molinfo top get numframes ] frames"
#       } 
# \endcode

## Iterator over files matched by the given pattern sorted in [natural
# order](http://sourcefrog.net/projects/natsort/), and execute
# block. See detailed description above.
# 
#       forFiles fn {*.dcd} {
#             puts "$fn has [ molinfo top get numframes ] frames"
#       } 
proc forFiles {filename pattern block} {
    set contents [lsort -dictionary [glob  $pattern ]]
    upvar 1 $filename pwvar

    foreach item $contents {
	animate delete all
 	set pwvar $item
	mol addfile $item  first 0 last -1 step 1  waitfor  all
	uplevel $block
    }
    animate delete all
}


##\private Iterate over chains matched in the named pattern prefix, a la
# GPUGRID.  Will discard currently-loaded frames!
#
#       forChains chain {../results} {
#              puts "DCD files prefixed with $chain have [ molinfo top get numframes ] frames" 
#       }
proc forChains {chain dir block} {
    set allfiles [glob -tails -directory $dir *.dcd ]
    foreach fn $allfiles { lappend allnames [ lindex [ split $fn - ] 0 ] }
    set chainlist [ lsort -uniq -dictionary $allnames ]
    upvar 1 $chain pwvar

    foreach item $chainlist {
	animate delete all
 	set pwvar $item
	foreach fn [ lsort -dictionary [ glob -directory $dir $item-*.dcd] ] {
	    mol addfile $fn  first 0 last -1 step 1  waitfor  all
	}
	uplevel $block
    }
    animate delete all
}

## @}


# ----------------------------------------



##\defgroup native Compute the number of native contacts
# 
# Computation of the number of native contacts requires two steps:
#
#  1. creation of a "reference" list, i.e.  the native state; 
#  2. the actual measurement. 
#
# In the current version, both should come
# from the same structure. Contacts is based on VMD's *measure
# contacts* feature and therefore can be performed on one
# (intra-molecular) or two (inter-molecular) selections. 
# 
# **Note:** For native contacts, the "reference" structure must be
# *equivalent* to the one to be analyzed (same number of atoms)
#
# For a GUI to compute native contacts, see the [RMSD trajectory tool
# + NC extension](http://www.multiscalelab.org/utilities/RMSDTTNC).
# 
# Example usage:
# \code 
#      set chaina [ atomselect top "chain A" ];
#      set ref [ prepareNativeContacts 7 $chaina ];
#      # (go or select the frame of interest)
#      measureNativeContacts $ref 7 $chaina;
# \endcode
# 
# Another (deliberately verbose) example:
# \code
#       ## The atom selection
#       set chaina [ atomselect top "chain A and name CA" ]
# 
#       ## Use first trajectory frame as a reference
#       $chaina frame 0
#       set ref [ prepareNativeContacts 7 $chaina ]
# 
#       ## Get the number of native contacts (for computing their fraction)
#       set nnc [ llength $ref ]
#       puts "There are $nnc contacts in the native state"
# 
#                     # Now, for each frame,
#       forFrames fno $chaina {
#                     # compute number of native contacts,
#               set nc [measureNativeContacts $ref 7 $chaina]  
#                     # their fraction,
#               set qnc [ expr 100.0 * $nc / $nnc ]
#                     # and print both.
#               puts [ format "Frame %d: %f, %.3f%%" $fno $nc $qnc ]
#       }
# \endcode
#
#  @{

## Given the frame, compute whether each residue in query is making
# contact (< cutoff) with target. Return list like { { R1 0 } {R2 1 }
# ... } where R1, R2 etc belong to query
proc residueContacts {cutoff query target} {

    # Uniquified list of residues (will be the columns in the ouptut)
    set kidres  [ lsort -unique -integer [ $query get resid ] ]

    # Residues in kid which are in contact (uniquified)
    set llkid [lindex [measure contacts $cutoff $query $target] 0]

    # if no residues make contact, append a dummy resid so that the
    # following does not fail
    if { $llkid == {} } {
	set kidresc { -1 }
    } else {
	set tmp     [ atomselect top "index $llkid"]
	set kidresc [ $tmp get resid ]
	$tmp delete
    }

    # For each KID residue, have a column 0 or 1 depending on whether it is 
    # listed in kidresc
    set kidlist {}
    foreach i $kidres {
	if { [lsearch -integer $kidresc $i ] != -1 } {
	    lappend kidlist [ list $i 1 ]
	} else {
	    lappend kidlist [ list $i 0 ]
	}
    }

    return $kidlist
}




## Like \ref  measureNativeContacts, but return resid's rather than
# indices. Return a list of (unique-ified) residue contacts as { R1 R2
# } { R1 R3 } ...  where the first element of a pair is a resid in
# sel1 and the second is in sel2. All computations are taken in the
# selections' frames. 
proc residueContactPairs { cutoff sel1 sel2 } {
    set cl [ measure contacts $cutoff $sel1 $sel2 ]

    if { $cl == {} } {
	return {}
    }

    # index->resid maps for the two selections
    set tmp [$sel1 get {index resid}]
    array set s1ir [concat {*}$tmp]
    set tmp [$sel2 get {index resid}]
    array set s2ir [concat {*}$tmp]

    # List of (duplicated) residues in contact
    set l1 [ lindex $cl 0 ]
    set l2 [ lindex $cl 1 ]
    set n [ llength $l1 ]

    # Now uniquify the pairs 
    array set pairs {}
    foreach i1 $l1 i2 $l2 {
	set r1 $s1ir($i1)
	set r2 $s2ir($i2)
	# don't count self-contacts
	if { $r1 != $r2 } {
	    # count AB / BA contacts only once
	    if { ! [ info exists pairs($r2,$r1) ] } {
		set pairs($r1,$r2) 1
	    }
	}
    }

    set plist {}
    foreach p [ array names pairs ] {
	set pl [ split $p , ]
	lappend plist $pl
    }

    # And return the contact pairs
    return $plist
}


## Return a matrix (as a serialized array) with the time for first
# contact of each contact pair. Times may be missing if they never do
# a first contact! Example:
#       set fctm [firstContactTimeMatrix 5 $a $b]
proc firstContactTimeMatrix { cutoff sel1 sel2 } {
    set res {}
    forFrames fn $sel1 {
	$sel2 frame $fn
	set cpairs [residueContactPairs $cutoff $sel1 $sel2]
	set ctpairs {}
	foreach r1r2 $cpairs {
	    lassign $r1r2 r1 r2
	    if {![info exists a($r1,$r2)]} {
		set a($r1,$r2) $fn
		lappend rl $r1
		lappend cl $r2
	    }
	}
    }
    return [list [array get a] [lsort -uniq -dict $rl] [lsort -uniq -dict $cl]]
}





## Prepare the "reference" list of native contacts, e.g. from the
# crystal structure. Return value: a list of native contact pairs
# (only useful to be passed as an argument to measureNativeContacts,
# or to get its length).
#
# *Note.* The function is subject to change. Later can save sizes of
# atomselections (for safety checking) and/or cutoff. 
proc prepareNativeContacts { cutoff sel1 {sel2 0} } {
    if { $sel2 == 0 } { set sel2 $sel1 }
    return [ transpose [ measure contacts $cutoff $sel1 $sel2 ] ]
}

## TBD.
proc measureNativeContacts { nclist cutoff sel1 {sel2 0} } {
	set n 0
	if { $sel2 == 0 } { set sel2 $sel1 }
	# current contacts
	set cc [ transpose [ measure contacts $cutoff $sel1 $sel2 ] ]
	# check if current is in native set (intersect)
	# iterate over current pairs cp
	foreach cp $cc {
		if { 	[lsearch $nclist $cp] != -1 ||
			[lsearch $nclist [lreverse $cp]] != -1 } {incr n }
	}
	return $n
}

## @}


# ----------------------------------------


##\defgroup dm Compute distance matrices
# @{
# Usage:
# \code
#      set tt [distanceProfileFull $b $a]; true  
# #or:
#      set tt [distanceProfile $b $a]; true  
# \endcode
# Then:
# \code
#      set ttf [open tt1.dat w]; printTable $tt "" $ttf; close $ttf 
#  #or:
#      set ttf [open tt2.dat w]; printTable [reshapeToWide $tt] "" $ttf; close $ttf
# \endcode


## Given two atom selection, return the distance matrix (in
# longitudinal form: <tt>{ { r1 r2 dist} ... }</tt> ) between residues
# in given selections. All selected atoms will be considered for each
# residue, and the minimum distance per residue pair returned. This
# can be slow. If using only an atom per residue (eg CA), consider
# using \ref atomDistanceMatrix for speed. Residues are taken from the
# *resid* attribute.
proc residueDistanceMatrix { sel1 sel2 } {
    set r {}
    
    foreach r1 [lsort -uniq -integer [$sel1 get resid] ] {
	set as1  [atomselect top "resid $r1"] 
	foreach r2 [lsort -uniq -integer [$sel2 get resid] ] {
	    set as2  [atomselect top "resid $r2"]
	    set d [minDist $as1 $as2]
	    lappend r [ list $r1 $r2 $d ]
	    $as2 delete
	}
	$as1 delete
    }
    return $r
}


## Similar to \ref residueDistanceMatrix, but much faster
# implementation which returns a discrete distance (taken from the
# bins list).  For each residue pair, the distance will be the next
# higher value in {bins}. Runtime is O([llength $bins]). Returns a
# list of <tt>{ {rid1 rid2 dist} ... }</tt>
proc residueDistanceMatrixApprox { sel1 sel2 bins } {
    array set rp {}
    array set idx1res {}
    array set idx2res {}

    set rbins [lsort -real -decreasing $bins]
    set rl1  [lsort -uniq -integer [$sel1 get resid] ] 
    set rl2  [lsort -uniq -integer [$sel2 get resid] ] 

    # Prepare index-to-residue mappings
    foreach idx1 [ $sel1 get index ] res1 [ $sel1 get resid ] { set idx1res($idx1) $res1 }
    foreach idx2 [ $sel2 get index ] res2 [ $sel2 get resid ] { set idx2res($idx2) $res2 }

    # Initialize to highest value of the distance list
    set maxd [lindex $rbins 0]
    set rbins [ lreplace $rbins 0 0 ]
    foreach r1 $rl1 {
	foreach r2 $rl2 {
	    set rp($r1,$r2) $maxd
	}
    }

    foreach d $rbins {
	set cl [ measure contacts $d $sel1 $sel2 ]
	foreach i1 [ lindex $cl 0 ] i2 [ lindex $cl 1 ] {
	    set r1 $idx1res($i1)
	    set r2 $idx2res($i2)
	    set rp($r1,$r2) $d
	}
    }

    set r {}
    foreach r1 $rl1 {
	foreach r2 $rl2 {
	    lappend r [list $r1 $r2 $rp($r1,$r2)]
	}
    }
    return $r
}



## Return the atom distance matrix (in longitudinal form: { { r1 r2
# dist} ... } ) between residues in given selections.  There should be
# only one atom per residue in the selection. Residues are taken from
# the ''resid'' attribute.
proc atomDistanceMatrix { sel1 sel2 } {
    set r {}
    set al1 [ $sel1 get index ]
    set al2 [ $sel2 get index ]
    set rl1 [ $sel1 get resid ]
    set rl2 [ $sel2 get resid ]
    foreach a1 $al1 r1 $rl1 {
	foreach a2 $al2 r2 $rl2 {
	    set d [ measure bond [ list $a1 $a2 ] ]
	    lappend r [ list $r1 $r2 [format "%.2f" $d] ]
	}
    }
    return $r
}



## Compute the distance profile between a ligand and a protein, i.e.,
# by timestep and residue, the minimum distance between the two.
# Requires two atomselections (likely from the same molecule) and an
# optional list of distance bins. Returns a list of { { frame
# resid_lig resid_prot distance } ... } . Unlike \ref distanceProfile, this can be post-processed
# to compute all kind of profiles and istance matrices. 
proc distanceProfileFull {lig prot {step 1} {dbins "2 3 4 5 6 7 10" }} {
    # count the number of protein residues
    set nligresid [llength [lsort -uniq -integer [$lig get resid]]]; 
    set nprotresid [llength [lsort -uniq -integer [$prot get resid]]]; 

    set omat {};		# output matrix

    # for each frame, compute the discretized distance matrix
    forFrames fn $lig {
	$prot frame $fn
	set dm [residueDistanceMatrixApprox $lig $prot $dbins]
	# dm is now a "long" list of {{rid1 rid2 dist} ...}
	# prepend the frame no and paste into omat
	foreach lpd  $dm {
	    lassign $lpd rlig rprot dist
	    lappend omat [list $fn $rlig $rprot $dist]
	} 
    } $step

    return $omat
}



## Compute the distance profile between a ligand and a protein, i.e.,
# by timestep and residue, the minimum distance between the two.
# Requires two atomselections (likely from the same molecule) and an
# optional list of distance bins. Returns a list of { { frame resid
# distance } ... }. 
proc distanceProfile {lig prot {dbins "2 3 4 5 6 7 10" }} {
    # count the number of protein residues
    set nprotresid [llength [lsort -uniq -integer [$prot get resid]]]; 

    # treat the ligand as a single residue
    set oligresid [$lig get resid]
    $lig set resid 10000;

    set omat {};		# output matrix

    # for each frame, compute the discretized distance matrix
    forFrames fn $lig {
	$prot frame $fn
	set dm [residueDistanceMatrixApprox $lig $prot $dbins]
	# replace unused ligand resid with frame number, so that we
	# can reuse the reshape functions etc. 
	set dmt [transpose $dm]
	# $dm are actually two lists
	foreach rid [lindex $dmt 1] dist [lindex $dmt 2] {
	    lappend omat [list $fn $rid $dist]
	} 
    }

    $lig set resid $oligresid
    return $omat
}

## @}


# ----------------------------------------


##\defgroup qplot Quick plots
#@{

## Quick plot function. Takes either a list of y values, or two lists,
# with x and y values.  If two vectors are given, they are interpreted
# as ''x'' and ''y'' values respectively.
proc qplot { li { ly 0}  } {
    if {$ly==0} {
	return [ multiplot -y $li -plot -lines -marker point \
		     -radius 2 -fillcolor "#ff0000" -color "#ff0000" ]
    } else {
	return [ multiplot -x $li -y $ly -plot -lines -marker point \
		     -radius 2 -fillcolor "#ff0000" -color "#ff0000" ]
    }
}

## Compute and plot the histogram of the values in the ''list'' passed
# as a second argument, binned in ''nbins'' equal-sized bins. Returns
# a list of two lists, the former being lower boundaries for the bins,
# and the latter are then counts.
proc qhist { bins li } {
    set li [lsort -real $li]
    set lmin [lindex $li 0]
    set lmax [lindex [lreverse $li] 0]
    set n [llength $li]
    set h {};			# counts
    set x0 {};			# bin centers
    for {set i 0} {$i<$n} {incr i} {
	lappend h 0
	lappend x0 [expr $lmin+$i*($lmax-$lmin)/$n]
    }
    foreach x $li {
	set b [expr int($n*($x-$lmin)/($lmax-$lmin))]
	if {$b==$n} { incr b -1 }
	lset h $b [expr [lindex $h $b]+1]
    }
    multiplot -x $x0 -y $h -plot -marker point -radius 2 \
	-ylabel Count  -fillcolor "#ff0000" -color "#ff0000"
    return [list $x0 $h]
}

## @}


# ----------------------------------------

##\defgroup lf Load multiple files
# @{

## Usage:  loadFrames 53-*.coor
# Will load all files starting with 53 (coordinates, dcd, whatever) in natural sort
# e.g. 44-bla-2-100.coor < 44-bla11-100.coor
proc loadFrames {pattern} {
    animate delete all
    appendFrames $pattern
}



## Usage: loadFrames 53-*.coor Will load all files starting with 53
# (coordinates, dcd, whatever) in natural sort, e.g. `44-bla-2-100.coor`
# < `44-bla11-100.coor`. Error handling is tricky. The
# /tmp/appendFrames.[pid].log will be true in any case. The ... .ff
# file will contain infos on the last successfully loaded frame for
# each trajectory. Frames are 0-based.
proc appendFrames {pattern} {
    set dir {.}
    set contents [lsort -dictionary [glob  $pattern ]]
    set ifile [open /tmp/appendFrames.[pid].log w]
    set ffile [open /tmp/appendFrames.[pid].ff w]

    foreach item $contents {
	set ofr [molinfo top get frame]
	if [catch {mol addfile $item  first 0 last -1 step 1  waitfor  all} err] {
	    puts "$err"
	} else {
	    set nfr [molinfo top get frame]
	    puts $ffile "$nfr $item"
	}
	set nfr [molinfo top get frame]
	puts $ifile "$item from $ofr to $nfr"
    }
    close $ifile
    close $ffile
}


## Identify current Frame
proc identifyFrame {} {
    set r -1
    set ifile [open /tmp/appendFrames.[pid].log r]
    while { [gets $ifile line] != -1 } {
	if { [lindex $line end] == [molinfo top get frame] } {
	    set r [lindex $line 0]
	    break
	}
    }
    close $ifile
    return $r
}

## @}

# ----------------------------------------



##\defgroup ren Renumber residues
# @{

## Renumbers the residues in the atom selection so that they start from
# the given integer. Useful to re-match standard numbering. Note that
# if there are duplicate residue IDs (e.g. in the case of
# homo-multimers), the renumbered IDs will also be duplicate.
proc renumber { sel start } {
	if { [$sel num] == 0 } {
		puts "Error in renumber: empty selection!"
		return
	}
	set oresid [ $sel get resid ]
	set delta [ expr $start - [ lindex $oresid 0] ]
	set nresid { }
	foreach r $oresid {
		lappend nresid [ expr $r + $delta ]
	}
	$sel set resid $nresid
}

## Renumbers residues starting from 1, just as '''tleap''' does. Note
# that '''all''' residues will be renumbered irrespective of their
# original ''resid''. As a consequence, the identity of homo-multimers
# will be lost.
proc renumber_from_1 { sel } {
	if { [$sel num] == 0 } {
		puts "Error in renumber: empty selection!"
		return
	}
	set oresid [ $sel get residue ]
	set delta 1
	set nresid { }
	foreach r $oresid {
		lappend nresid [ expr $r + $delta ]
	}
	$sel set resid $nresid
}

## @}

# ----------------------------



##\defgroup geo Structural analysis and geometry
# @{


## Count the fraction of residues that have phi/psi in the canonical
# alpha region of the Ramachandran plot (&phi;,&psi;)=(-57,-47). Pass
# a selection of CA only.  The ''tolerance'' argument sets the allowed
# slack (default 40 degrees). The result is normalized to 1.0 for a
# fully alpha-helical segment.
proc helicity { sel { tol 40 } } {
    set n [$sel num]
    set nh 0.0
    set tol 40
    set phi0 -57
    set psi0 -47
    set phil [ expr $phi0 - $tol ]
    set phih [ expr $phi0 + $tol ]
    set psil [ expr $psi0 - $tol ]
    set psih [ expr $psi0 + $tol ]

    foreach phi [ $sel get phi ] psi [ $sel get psi ] {
	if {  $phil < $phi && $phi < $phih &&  $psil < $psi && $psi < $psih } {
	    set nh [ expr $nh + 1.0 ]
	}
    }
    return [ expr $nh / $n ]
}



## Count helicities as in: Kelley J Mol Biol. 2009 May 22;
# 388(5):919–927.  At least 3 residues within the canonical
# Ramachandran alpha region are required to count one. Pass a
# selection of CA only. The result is normalized to 1.0 for a fully
# alpha-helical segment.
proc helicity_3 { sel { tol 40 } { debug 0 } } {
    set n [$sel num]
    set nh 0.0
    set tol 40
    set phi0 -57
    set psi0 -47
    set phil [ expr $phi0 - $tol ]
    set phih [ expr $phi0 + $tol ]
    set psil [ expr $psi0 - $tol ]
    set psih [ expr $psi0 + $tol ]
    
    set hi 0
    set s ""

    foreach phi [ $sel get phi ] psi [ $sel get psi ] {
	if {  $phil < $phi && $phi < $phih &&  $psil < $psi && $psi < $psih } {
	    incr hi
	    set s "${s}H"
	} else {
	    set hi 0
	    set s "${s}_"
	}
	if { $hi >= 3 } {
	    set nh [ expr $nh + 1.0 ]
	}
    }
    if {$debug} { puts "DEBUG $s $nh over $n" }
    return [ expr $nh / ($n-2) ]
}

## Refine an atom selection: return a new atomselection consisting of
# the atoms in the old one, as long as they ALSO match the selection
# text txt
proc atomselectRefine {sel txt} {
	set mid [$sel molid]
	set otx [$sel text]
	set ns [atomselect $mid "($otx) and ($txt)"]
	$ns uplevel 1
	return $ns
}


## Return a zero-centered version of the input list
proc veccenter {l} {
    set m [vecmean $l]
    set N [llength $l]
    set mN [lrepeat $N $m]
    set r [vecsub $l $mN]
    return $r
}

## Returns the angle between two vectors
proc vecangle {d1 d2} {
	set cosangle [expr [vecdot $d1 $d2]/[veclength $d1]/[veclength $d2]]
	return [expr acos($cosangle)*180./3.141592653589793]
}


## Return the matrix which reorients the principal axes of inertia with
# x, y, z. The largest inertia axis will be aligned along z.
proc transinertia {sel} {
	set m [lindex [measure inertia $sel] 1]
	set v0 "[lindex $m 0] 0"
	set v1 "[lindex $m 1] 0"
	set v2 "[lindex $m 2] 0"
	set v3 "0 0 0 1"
	set m4 [list $v0 $v1 $v2 $v3]
	return $m4
}

## Return the center of the bounding box for the selection $sel.
proc boundingBoxCenter {sel} {
	set bb [measure minmax $sel]
	set a [lindex $bb 0] 
	set b [lindex $bb 1]
	# Box size
	set delta [vecsub $b $a]
	# Half box
	set deltah [vecscale .5 $delta]
	# Center
	set centroid [vecadd $a $deltah]
	return $centroid
}



## Compute the minimum distance between atoms in selections s1 and s2
proc minDist {s1 s2} {
    set md 1e6
    foreach i1 [$s1 get index] {
	foreach i2 [$s2 get index] {
	    if { $i1 == $i2 } {
		set d 0
	    } else {
		set d [ measure bond [ list $i1 $i2 ] ]
	    }
	    if { $d < $md } {
		set md $d 
	    }
	}
    }
    return $md
}

##\private TBD
proc doudouVolume {dx dy dz kx ky kz {kbt 0.59}} {
    set pi 3.14159265358979
    set dx1  [expr $dx + sqrt(2*$pi*$kbt/$kx) ]
    set dy1  [expr $dy + sqrt(2*$pi*$kbt/$ky) ]
    set dz1  [expr $dz + sqrt(2*$pi*$kbt/$kz) ]
    return [expr $dx1 * $dy1 * $dz1 ]
}

## @}

# --------------------------------

##\defgroup manip Structural manipulation
# @{

## Mutate first and last residue of a selection so that tleap will turn them
# into ACE and NME caps. "cap" may be "ACE", "NME", or "both" (default). 
proc addCaps { sel {cap both} } {
	set rl [lsort -integer -unique [$sel get residue]]
	set mid [$sel molid]
	if { $cap == "ACE" || $cap == "both" } {
		set ri [lindex $rl 0]
		set h2 [atomselect $mid "residue $ri and name H2"]
		$h2 set name C
		$h2 set resname ACE
		$h2 set resid [expr [$h2 get resid]-1]
		$h2 delete
	} 
	if { $cap == "NME" || $cap == "both" } {
		set rf [lindex $rl end]
		set oxt [atomselect $mid "residue $rf and name OXT"]
		$oxt set name N
		$oxt set resname NME
		$oxt set resid [expr [$oxt get resid]+1]
		$oxt delete
	} 
	puts "Structure mutated, now save PDB"
}


## Stronger version of \ref addCaps - replaces whatever first and last atoms
proc addCaps2 {sel {cap both}} {
	set sl [lsort -integer -unique [$sel get index]]
	set mid [$sel molid]
	if { $cap == "ACE" || $cap == "both" } {
		set ii [lindex $sl 0]
		set h2 [atomselect $mid "index $ii"]
		$h2 set name C
		$h2 set resname ACE
		$h2 set resid [expr [$h2 get resid]-1]
		$h2 delete
	} 
	if { $cap == "NME" || $cap == "both" } {
		set if [lindex $sl end]
		set oxt [atomselect $mid "index $if"]
		$oxt set name N
		$oxt set resname NME
		$oxt set resid [expr [$oxt get resid]+1]
		$oxt delete
	} 
	puts "Structure mutated, now save PDB"
}

## Move a selection so that its center of mass is placed at the given point.
# \param s1 The atomselection
# \param dest Destination point, as in `{x y z}`
proc moveSelTo { s1 dest } {
    set c [measure center $s1]
    $s1 moveby [vecsub $dest $c]
}

## Exchange the positions of two atomselections *s1* and *s2* (based
# on their centers).
proc swap { s1 s2 } {
    set mm [ trans center [ measure center $s1 ] offset [ measure center $s2 ] ]
    $s1 move $mm
    $s2 move [ measure inverse $mm ]
}

## @}


# -----------------------------------------

##\defgroup rms Root-mean square calculations 
# @{
#
# Example
# -------
#
# RMSD for all files (in a directory) and frames
#
# Variables \c $compare and \c $reference should be two atomselections
# in different molecules. The former should be the TOP molecule.
#
# \code
#     forFiles id {../filtered/ *.dcd} {
#         set outch [open $id.rmsd w]
#         forFrames fn $compare {
#             set trans_mat [measure fit $compare $reference]
#             $compare move $trans_mat
#             set rmsd [measure rmsd $compare $reference]
#             puts $outch "$rmsd $id $fn"
#         }
#         close $outch
#     }
# \endcode


## Compute RMSD, over all frames, of a selection of atoms *sel1* with
# respect to another *sel2* after aligning the set of atoms *ref1* to
# *ref2*. In ther words, for per each frame, align *ref1* to *ref2*,
# and measure RMSD of *sel1* with respect to *sel2*. *sel1* and *ref1*
# should belong to the same molecule (the trajectory under study,
# multiple frames).  *Sel2* and *ref2* should belong to the same
# molecule (the reference frame). Return a list of RMSD values (one
# per frame in *sel*). If *sel2* is the string `ROTATE`, RMSD is not
# computed, but *sel1* is rotated instead.
proc rmsdOf { sel1 sel2 ref1 ref2 } {
    set rmsdlist {}
    forFrames fn $sel1 {
	set oco [ $sel1 get { x y z } ]
	$ref1 frame $fn
	set xform [measure fit $ref1 $ref2]
	$sel1 move $xform
	if {$sel2 != "ROTATE"} {
	    lappend rmsdlist [measure rmsd $sel1 $sel2]
	    $sel1 set {x y z} $oco
	}
    }
    return $rmsdlist
}

## Compute average rmsf by sliding windows of width win. RMSF will be
# averaged by weight. See details in \ref rmsfTrajectoryColor.
# \param sel  atom selection
# \param win  window size (frames) 
# \param step stride 
proc rmsfTrajectory {sel {win 10} {step 1}} {
    set n [ molinfo [$sel molid] get numframes ]
    set m [ $sel get mass ]
    set mt [vecsum $m];		# total mass
    if { $mt==0 } {
	puts "Unknown masses, assuming all unitary"
	set m [ lrepeat 1 [llength $m]]
	set mt [vecsum $m]
    }
    set m [vecscale [expr 1./$mt] $m]
    
    set res {}
    for {set fr 0} {$fr<=[expr $n-$win]} {incr fr $step} {
	set rmsf [measure rmsf $sel first $fr last [expr $fr+$win-1]]
	set armsf [vecdot $m $rmsf]; # weighted average
	lappend res $armsf
    }

    # Pad the last frames with the last armsf value
    for {} {$fr < $n} {incr fr $step} {
	lappend res $armsf
    }

    return $res
}


## Compute average rmsf by sliding windows of width win. RMSF will be
# averaged by weight. Assign it to the "user" attribute at each frame, 
# which will hold the RMSF computed over the next *win* frames.
# \param sel  atom selection
# \param win  window size (frames)
proc rmsfTrajectoryColor {sel {win 10}} {
    set n [ molinfo [$sel molid] get numframes ]

    # Initial $win frames
    set rmsf [measure rmsf $sel first 0 last [expr $win-1]]
    for {set fr 0} {$fr<[expr $win/2]} {incr fr} {
	$sel frame $fr
	$sel set user $rmsf
    }

    # Final $win frames
    set rmsf [measure rmsf $sel first [expr $n-$win-1] last [expr $n-1]]
    for {set fr 0} {$fr<[expr $win/2]} {incr fr} {
	$sel frame [expr $n-$fr-1]
	$sel set user $rmsf
    }
    
    # Decent ones
    for {set fr 0} {$fr<=[expr $n-$win]} {incr fr} {
	set rmsf [measure rmsf $sel first $fr last [expr $fr+$win-1]]
	$sel frame [expr $fr+$win/2.]
	$sel set user $rmsf
    }
}

## @}

# ----------------------------


##\defgroup format Format conversions
# @{

## Write a crude PQR file using radiuses and masses in the topology
# (radii may not be appropriate for APBS calculations! use pdb2pqr
# instead!)  Note that this preserves the CHAIN id. VMD PQR loader
# misparses the file.
proc writePQR { sel filename } {
	set oldbeta [ $sel get beta ]
	set oldocc [ $sel get occupancy ]
	$sel set beta [ $sel get radius ]
	$sel set occupancy [ $sel get charge ]	
	$sel writepdb $filename
	$sel set beta $oldbeta
	$sel set occupancy $oldocc
}


## Write a PDB file using charges and masses in the topology
# (required by PLUMED's "driver" utility)
proc writePlumed { sel filename } {
	set oldbeta [ $sel get beta ]
	set oldocc [ $sel get occupancy ]
	$sel set beta [ $sel get charge ]
	$sel set occupancy [ $sel get mass ]	
	$sel writepdb $filename
	$sel set beta $oldbeta
	$sel set occupancy $oldocc
}

## Write a PDB file to be used as a reference in PATH cvs.  Assumes
# that occupancy and beta are set according to the align intentions.
proc writePlumedRef { sel filename } {
    set oldsegid [$sel get segid]
    set tmp "tmp_[pid]"

    set newsel [atomselectRefine $sel "beta != 0 or occupancy != 0"]
    $newsel set segid KEEP
    $sel writepdb $tmp
    $newsel delete

    set n 0
    set fdr [open $tmp r]
    set fdw [open $filename w]
    while { [gets $fdr line] != -1 } {
	if { [ regexp {KEEP} $line ] } {
	    # workaround plumed bug in PDB reader
	    incr n
	    set line [string replace $line 21 21 " "]
	    puts $fdw $line
	}
    }
    close $fdr; close $fdw
    file delete $tmp
    puts "Written $n lines"

    $sel set segid $oldsegid
}


## Write a PDB file containing TER cards to separate fragments. This is
# useful for software like ''tleap'' which requires them. ''$sel''
# must be an atom-selection function (as returned by atomselect). For
# a fuller implementation, see also [[utilities/PdbTer]].
proc writePDBTER { sel fname } {
    set tempf "$::env(TMPDIR)/writepdbter.[pid]"
    set oresid [$sel get resid]
    set ochain [$sel get chain]
    set osegid [$sel get segname]
    set oresidue [$sel get residue]
    
    # Need unique segids even in 4 chars
    foreach f  [$sel get fragment] {
    	lappend nsegid [expr $f % 10000]
    }
    
    # Workaround tleap's requirement of unique chain/resid combination
    foreach r $oresidue {
        lappend nresid [expr $r % 10000]
        lappend nchain [expr int($r/10000)]
    }
    $sel set resid $nresid
    $sel set chain $nchain
    $sel set segid $nsegid

    $sel writepdb $tempf

    $sel set resid $oresid
    $sel set chain $ochain
    $sel set segid $osegid

    set ch [ open $tempf r ]
    set och [ open $fname w ]
    set osid ""
    set orty ""
    while { -1 != [ gets $ch line ] } {
	set nsid [string range $line 72 75]
	set nrty [lindex $line 0]
	if { $nrty == "ATOM" && $orty == "ATOM" && $nsid != $osid } {
	    puts $och "TER"
	}
	puts $och $line
	set orty $nrty
	set osid $nsid
    }
    close $och
    close $ch
}


## Write null velocity file corresponding to the current structure
proc writeNullVelFile {as fname} {
    set oxyz [$as get {x y z}]
    $as set x 0;
    $as set y 0;
    $as set z 0;
    animate write namdbin $fname sel $as
    $as set {x y z} $oxyz
}

# proc writeZeroVelFile {n filename} {
# 	set fp [open $filename w]
# 	set m [expr $n*24]
# 	puts -nonewline $fp [binary format "i1x$m" $n]
# 	close $fp
# }



## Returns the 1-letter sequence of a selection (looking at CA only)
proc getFasta {osel} {
	array set atable {
		ALA	A		ARG	R
		ASN	N		ASP	D
		CYS	C		GLU	E
		GLN	Q		GLY	G
		HIE     H               HID     H
		HIS	H		ILE	I
		LEU	L		LYS	K
		MET	M		PHE	F
		PRO	P		SER	S
		THR	T		TRP	W
		TYR	Y		VAL	V
		SEC 	U		PYL	O
		CYX     C
	}

	set sel [atomselectRefine $osel "name CA"]
	set idl [$sel get resid]
	set nl [$sel get resname]
	set oid [expr [lindex $idl 0]-1]
	set seq ""
	foreach id $idl n $nl {
		if {[info exists atable($n)]} {
			set olc	$atable($n)
		} else {
			puts "Unknown resname $n"
			set olc "?"
		}
		set seq "${seq}$olc"
		if {$id != [expr $oid+1]} {
			puts "Discontinuity between RESIDs $oid and $id"
			set seq "${seq}/"
		}
		set oid $id
	}
        $sel delete
	return $seq
}

## @}




# ---------------------------------


##\defgroup graphics Miscellaneous graphics operations
# @{
#
# \ref labelAtoms shows labels with a chosen format on all atoms 
# matched by the given selection (string or atomselection object, default: "all").
# Existing labels are deleted. 
# For example: \code 
#   labelAtoms {%t} {element C}
# \endcode
proc labelAtoms {txt {sel all}} {
	if [ regexp {^atomselect} $sel ] {
		set as $sel
	} else {
		set as [atomselect top $sel]
	}

	label delete Atoms
	set i 0
	foreach a [$as list] {
		set tmp [format "%d/%d" [$as molid] $a]
		label add Atoms $tmp
		label textformat Atoms $i $txt
		incr i
	}
	if [ ! regexp {^atomselect} $sel ] { $as delete }
}


## @}


# ---------------------------------


##\defgroup save Save VMD state
# @{
#
# \ref saveFullState is a function to **almost completely** save the
# state of a VMD session. Saved state includes the view, materials,
# colors (as with the usual save state menu item), as well as
# *trajectories* and *topologies*.  Loading the saved state completely
# restores the state a VMD session was in.
#
# Usage
# -----
#
# State can then be saved with 
# \code saveFullState filename.vmd \endcode
#
# This will create:
#
# * a \b filename.vmd script, which can be reloaded to recover the state;
# * a \b filename.vmd.d directory, containing trajectory data and a snapshot image.
#
# If `filename` is not provided, a file selection dialog opens up.
#
# To recover a session, just `source` the script, or use it with the `vmd -e` option.
#
# Trajectory files are looked up with respect to the current directory. If they
# are not found, a second attempt is done with the absolute pathnames that were 
# valid at the time the state was saved.
#
#
# Limitations
# -----------
#
# Some data is known not to be restored, including:
#
# * volumetric data
# * time-varying variables, such as *user*
#
#
# To do
# -----
#  - preserve volume data
#  - preserve time-varying data (e.g. user)
#  - option to save only current frame
#  - preserve inferred topologies (parse filespec?)
#  - handle relative paths when loading states [DONE]
#
#
# Copyright
# ---------
#
# Most of this proc is copied from the \c save_state
# function built-in VMD, which is distributed under the terms of the
# UIUC Open Source License:
# http://www.ks.uiuc.edu/Research/vmd/current/LICENSE.html
#
#
# Modifications by Toni Giorgino <toni.giorgino isib.cnr.it>

## Save FULL visualization state, including trajectories Can be
# reloaded as usual.
proc saveFullState {{file EMPTYFILE}} {
  global representations
  global viewpoints
  save_viewpoint
  save_reps

  # If no file was given, get a filename.  Use the Tk file dialog if 
  # available, otherwise get it from stdin. 
  if {![string compare $file EMPTYFILE]} {
    set title "Enter filename to save current VMD state:"
    set filetypes [list {{VMD files} {.vmd}} {{All files} {*}}]
    if { [info commands tk_getSaveFile] != "" } {
      set file [tk_getSaveFile -defaultextension ".vmd"  -title $title -filetypes $filetypes]
    } else {
      puts "Enter filename to save current VMD state:"
      set file [gets stdin]
    }
  }
  if { ![string compare $file ""] } {
    return
  }

  set fildes [open $file w]
  puts $fildes "\#!/usr/local/bin/vmd"
  puts $fildes "\# VMD script written by save_full_state \$Revision: 1.44 $"

  set vmdversion [vmdinfo version]
  puts $fildes "\# VMD version: $vmdversion"

  puts $fildes "set viewplist {}"
  puts $fildes "set fixedlist {}"
  save_materials     $fildes
  save_atomselmacros $fildes
  save_display       $fildes

  # Data directory. Reusing file's path components is fine here,
  # i.e. when writing.  Not so in the state script, ie. when reading,
  # because CWD can be different. Proposed solution: use RELATIVE
  # pathnames in the state file, and make "load_state" resolve them
  # from the selected file (e.g. temporarily changing
  # directory). Alas, this would not solve other ways of loading
  # files, e.g. "play". Converting to absolute pathnames would work,
  # but implies that states can't be moved.
  set fildir $file.t
  set fildir_abs [file join [pwd] $fildir]
  file delete -force $fildir
  file mkdir $fildir
  puts "Saving trajectory data files in directory $fildir"

  render snapshot $fildir/preview.tga 

  # Does it depend on the top molecule?
  set current_frame [molinfo top get frame]

  foreach mol [molinfo list] {
      set mname [format "mol%04d" $mol]

      # in lack of a native format
      animate write psf $fildir/$mname.psf waitfor all $mol
      animate write pdb $fildir/$mname.pdb beg 0 end 0 waitfor all $mol
      animate write dcd $fildir/$mname.dcd waitfor all $mol

      # Try relative first, then absolute if it fails
      puts $fildes "if \[catch {
                                mol new $fildir/$mname.pdb waitfor all
                                animate delete all;
                                mol addfile $fildir/$mname.psf waitfor all
                                mol addfile $fildir/$mname.dcd waitfor all
                    } e \] {
                                puts \"Couldn't open original pathnames; trying $fildir_abs\";
                                mol new $fildir_abs/$mname.pdb waitfor all
                                animate delete all;
                                mol addfile $fildir_abs/$mname.psf waitfor all
                                mol addfile $fildir_abs/$mname.dcd waitfor all
                           } "
      
    # We load PDB for chain info, beta/okkupa, PSF for resid, topology, masses etc,
    # and DCD for coordinates. DCD could be replaced by a multi-frame
    # PDB, with the exception of (a) impossibility to save large
    # coordinates and (b) occupies more disk space.

    foreach g [graphics $mol list] {
      puts $fildes "graphics top [graphics $mol info $g]"
    }
    puts $fildes "mol delrep 0 top"
    if [info exists representations($mol)] {
      set i 0
      foreach rep $representations($mol) {
        foreach {r s c m pbc numpbc on selupd colupd colminmax smooth framespec cplist} $rep { break }
        puts $fildes "mol representation $r"
        puts $fildes "mol color $c"
        puts $fildes "mol selection {$s}"
        puts $fildes "mol material $m"
        puts $fildes "mol addrep top"
        if {[string length $pbc]} {
          puts $fildes "mol showperiodic top $i $pbc"
          puts $fildes "mol numperiodic top $i $numpbc"
        }
        puts $fildes "mol selupdate $i top $selupd"
        puts $fildes "mol colupdate $i top $colupd"
        puts $fildes "mol scaleminmax top $i $colminmax"
        puts $fildes "mol smoothrep top $i $smooth"
        puts $fildes "mol drawframes top $i {$framespec}"
        
        # restore per-representation clipping planes...
        set cpnum 0
        foreach cp $cplist {
          foreach { center color normal status } $cp { break }
          puts $fildes "mol clipplane center $cpnum $i top {$center}"
          puts $fildes "mol clipplane color  $cpnum $i top {$color }"
          puts $fildes "mol clipplane normal $cpnum $i top {$normal}"
          puts $fildes "mol clipplane status $cpnum $i top {$status}"
          incr cpnum
        }

        if { !$on } {
          puts $fildes "mol showrep top $i 0"
        }
        incr i
      } 
    }
    puts $fildes [list mol rename top [lindex [molinfo $mol get name] 0]]
    if {[molinfo $mol get drawn] == 0} {
      puts $fildes "molinfo top set drawn 0"
    }
    if {[molinfo $mol get active] == 0} {
      puts $fildes "molinfo top set active 0"
    }
    if {[molinfo $mol get fixed] == 1} {
      puts $fildes "lappend fixedlist \[molinfo top\]"
    }

    puts $fildes "set viewpoints(\[molinfo top\]) [list $viewpoints($mol)]"
    puts $fildes "lappend viewplist \[molinfo top\]"
    if {$mol == [molinfo top]} {
      puts $fildes "set topmol \[molinfo top\]"
    }
    puts $fildes "\# done with molecule $mol ----------------------------------"
  } 
  puts $fildes "foreach v \$viewplist \{"
  puts $fildes "  molinfo \$v set {center_matrix rotate_matrix scale_matrix global_matrix} \$viewpoints(\$v)"
  puts $fildes "\}"
  puts $fildes "foreach v \$fixedlist \{"
  puts $fildes "  molinfo \$v set fixed 1"
  puts $fildes "\}"
  puts $fildes "unset viewplist"
  puts $fildes "unset fixedlist"
  if {[llength [molinfo list]] > 0} {
    puts $fildes "mol top \$topmol"
    puts $fildes "unset topmol"
  }
  save_colors $fildes
  save_labels $fildes
  
  puts $fildes "animate goto $current_frame"
    
  close $fildes
}


## Create a list of "mol" commands that reproduce the current top
# molecule display.
proc dumpRepresentations {} {
 #    foreach fn [lindex [molinfo top get filename] 0] { 	puts "mol addfile $fn"     }
    puts "mol new [lindex [molinfo top get filename] 0 0] waitfor all"
    puts "mol addfile [lindex [molinfo top get filename] 0 end] waitfor all"
    puts "mol delrep 0 top"
    for {set i 0} {$i < [molinfo top get numreps]} {incr i} {
	lassign [molinfo top get "{rep $i} {selection $i} {color $i}"] a b c
	puts "  mol rep $a"
	puts "  mol selection $b"
	puts "  mol color $c"  
	puts "  mol addrep top"  
    }
    puts "  animate goto [molinfo top get frame]"
}


## @}





#\defgroup large Processing large trajectories in-memory
# @{

##\private Process in-memory a large trajectory that came from several different files.
# Something like
#
#        set gi [loadFrameGroup *.dcd 10];
#        forFrameGroup i g $gi { puts "$i $g ..." }
 proc loadFrameGroups {pattern {step 1}} {
    set flist [lsort -dictionary [glob $pattern]]
    set nframes {}
    foreach f $flist { 
	lappend nframes [molinfo top get numframes]
	mol addfile $f step $step waitfor all 
    }
    return [ list $flist $nframes ]
}

##\private Process in-memory a large trajectory that came from several different files.
# Something like
#
#        set gi [loadFrameGroup *.dcd 10];
#        forFrameGroup i g $gi { puts "$i $g ..." }
 proc forFrameGroups {framenum groupnum groupinfo  block} {
    upvar 1 $framenum  l_framenum
    upvar 1 $groupnum  l_groupnum

    set groupname_list [lindex $groupinfo 0]
    set groupstart_list  [lindex $groupinfo 1]

    set nmax [ molinfo top get numframes ]
    set ngroups [llength $groupname_list]
    lappend groupstart_list $nmax

    for {set gr 0} { $gr < $ngroups } { incr gr } {
	set gr_start [ lindex $groupstart_list $gr ]
	set gr_end [ lindex $groupstart_list [ expr $gr+1 ] ]
	set l_groupnum $gr
	for { set i $gr_start } { $i < $gr_end } { incr i } {
	    set l_framenum [expr $i-$gr_start]
	    uplevel $block
	}
    }
}
## @}

# ----------------------------------------


##\defgroup list List and matrix operations
# @{
# Miscellaneous TCL data-structure related operations.

## Transpose a 2D table (list of lists). From [[http://wiki.tcl.tk/2748]]
proc transpose {matrix} {
    set cmd list
    set i -1
    foreach col [lindex $matrix 0] {append cmd " \$[incr i]"}
    foreach row $matrix {
        set i -1
        foreach col $row {lappend [incr i] $col}
    }
    eval $cmd
}


## Sequence of integers from $from to $to, step by $step
proc lseq {from to {step 1}} {
    set ltmp {}
    for {set tmp $from} {$tmp<=$to} {incr tmp $step} {lappend ltmp $tmp}
    return $ltmp
}


## Given an array indexed by A,B pairs, return sorted row and column indices.
proc arrayIndices { arr } {
    array set a $arr
    set names [array names a]
    foreach i $names {
	lassign [split $i ,] r c
	lappend rl $r
	lappend cl $c
    }
    return [list [lsort -dict -uniq $rl] [lsort -dict -uniq $cl]]
}



## Pretty-print a 2D array (list of lists). Each line is prepended by
# the optional "extra" argument, a string, and output on the given
# channel (stdout if not given).
proc printTable { ll { prepend "" } { ch stdout } { fmt "%s " } } {
    foreach r $ll {
	puts -nonewline $ch "$prepend"
	foreach c $r {
	    puts -nonewline $ch [ format $fmt $c ]
	}
	puts  -nonewline $ch "\n"
    }
}



## Pretty-print an array, prepending rows and columns
proc printArray { tmp  { prepend "" } { ch stdout } { fmt "%s " }  } {
    lassign $tmp tmp2 rl cl
    array set a $tmp2
    puts $ch "     $cl";		# column headers
    foreach r [lsort -dict $rl] {
	puts -nonewline $ch "$r "
	foreach c [lsort -dict $cl] {
	    if [info exists a($r,$c)] {
		set v $a($r,$c)
	    } else {
		set v NA
	    }
	    puts -nonewline $ch "$v "
	}
	puts $ch " "
    }
}




## Reshape a matrix in "long" format into an array format.  No order is
# assumed. Returns a list containing: 1. the serialized array; 2. the
# list of rows (first index); 3. the list of columns (second
# index). This is a convenient matrix format for interchanging data.
proc reshapeLongToArray { table } {
    foreach tuple $table {
	lassign $tuple r c v
	lappend rl $r
	lappend cl $c
	set a($r,$c) $v
    }
    set rl [lsort -uniq -int $rl]
    set cl [lsort -uniq -int $cl]
    return [list [array get a] $rl $cl]
}

## Reshape array format into "wide".  No order is
# assumed. 
proc reshapeArrayToWide { tmp } {
    lassign $tmp tmp2 rl cl
    array set a $tmp2
    # prepare the table for returning
    set res {}
    foreach r [lsort -dict $rl] {
	set row {}
	foreach c [lsort -dict $cl] {
	    if [info exists a($r,$c)] {
		lappend row $a($r,$c)
	    } else {
		lappend row NA
	    }
	}
	lappend res $row
    }
    return $res
}


## Reshape a matrix in "long" format into a "wide" (rectangular) format
# The first column will be used as row index, the second as column.
# No order is assumed.
proc reshapeLongToWide { table } {
    set tmp [reshapeLongToArray $table]
    return [reshapeArrayToWide $tmp]
}


## Reverse of \ref reshapeLongToWide
proc reshapeArrayToLong { tmp } {
    lassign $tmp tmp2 rl cl
    array set a $tmp2
    foreach i [array names a] {
	lassign [split $i ,] r c
	set v $a($i)
	lappend res "$r $c $v"
    }
    return $res
}



## @}



# UTILITY --------------------------------------------------

##\defgroup asel_group Auto-cleanup atomselections
# @{
# The \ref asel and \ref asel_free functions provide a replacement for
# VMD's `atomselect` which behave like VMD's usual atomselect, but
# keeps a list of selections created in the scope of the calling
# function (it is called `_asel_list`, but you should not need
# it). Such atomselections can be deleted all at once when the
# function is exited via the \ref asel_free function. 
#
# \warning these functions are experimental and unnecessary for most
# people, because the behavior of recent VMD versions is to clean up
# atomselections created in a `proc` anyway (unless they are declared
# `global`).
#
# In short, just use \ref asel instead of `atomselect`, and call \ref
# asel_free at the end of your function.
# 
# The \ref ssel function allows one to write function which accept
# either selection strings or existing atomselections.


## Auto-cleanup version of `atomselect`. 
proc asel {args} {
    upvar _asel_list l
    set s [eval "atomselect $args"]
    $s global ;# uplevel
    lappend l $s
    return $s
}

## Auto-cleanup version of `atomselect`. If `str` is an atomselection,
# just return it as it is. Otherwise, make a selection of it with the
# top molecule, record it like \ref asel, and return it.
# \param str An atom selection string or an existing atomselection
proc ssel {str} {
    if [string is integer [$str num]] {
	return $str
    } else {
	upvar _asel_list l
	set s [atomselect $m $str]
	$s global ;# uplevel
	lappend l $s
	return $s
    }
}


## Cleanup the atom-selections performed by \ref asel and \ref ssel in the 
# current function. 
proc asel_free {} {
    upvar _asel_list l
    foreach s $l {
	puts "deleting $s"
	$s delete
    }
    set l {}
}


## @}

# http://wiki.tcl.tk/10876
# init.tcl
# Copyright 2001,2005 by Larry Smith
# Wild Open Source, Inc
# For license terms see "COPYING"

##\private Takes a list of variable-value pairs and creates and
# initializes the former with the latter in the calling
# context.  If the first parameter is "-using" it will
# same away the following argument and use it in a second
# pass to alter values set in the first pass.  In will not
# create or modify vars not mentioned in the initializing
# list.  Returns a list of unrecognized var val pairs.
# If a -var is followed by another -var or by the end of
# the list, the var is set to 1.

proc init { args } {
    if { [ llength $args ] == 0 } return
    if { [ llength $args ] == 1 } { eval set args $args }
    set arglist {}
    if { [ string equal -using [ lindex $args 0 ] ] } {
	set arglist [ lindex $args 1 ]; set args [ lrange $args 2 end ]
    }
    set varlist {}
    set rest {}
   # first pass - args on command line
    foreach { var val } $args {
	uplevel 1 set $var \{$val\}
     lappend varlist $var
    }
   # second pass - args passed in -using list
    set len [ llength $arglist ]
    for { set i 0 } { $i < $len } { incr i } {
	set var [ string range [ lindex $arglist $i ] 1 end ]
     incr i
	if { $i < [ llength $arglist ] } {
	    set val [ lindex $arglist $i ]
	} else {
       set val 1
	}
	if { [ string index $val 0 ] == "-" } {
       incr i -1
       set val 1
	}
	if { [ lsearch $varlist $var ] != -1 } {
	    uplevel 1 set $var \{$val\}
	} else {
       lappend rest -$var $val
	}
    }
   return $rest
}




##\private A tentative procedure to embed gnuplot plots.   Usage:
#      set where .c
#      canvas $where
#      pack $where -fill both -expand true
#      eval [gp {plot sin(x)}]
#      gnuplot $where
#      bind . <Configure> { gnuplot $where }
 proc gp { cli } {
    set in "set term tk
$cli
exit"
    set io [open "| gnuplot" r+]
    puts $io $in
    flush $io
    set msg [read $io]
    close $io
    return $msg
}

