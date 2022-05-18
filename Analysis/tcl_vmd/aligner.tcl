#package require Orient
#namespace import Orient::orient
puts "aligner loaded..."
puts "Available: geom_center {selection}\nalignGeomCenter {molID1 molID2 selText1 selText2}\nalign2XYZ {molID selText {frame 0}}\nalignSelections\ncenterByMol { mID centerMolselText outFilename patternFrameNum}\nprogressiveFit { mID outFilename {saveMatrix 0}} "
#geom_center -return vector of geometrical center of selection
#alignGeomCenter -translate selection from "selText1" by the difference of geometrical centerers of selections "selText2" - "selText1"
#centerSel {selection} - Move geometrical center of selection to {0 0 0}
#align2XYZ - #move all molecule to geometrical center of selection and align its selection axis with XYZ
#alignSelections - move selection "sel2" over "sel1" in sense of best least square fit
#centerByMol - allign selection ("centerMolselText") in all frames in trajectory to best overlap (in sense of best least square fit) \\
#               with selection from frame "patternFrameNum" and output the trajectory (dcd file type) to file "outFilename"
#progressiveFit - align all atoms in all frames in trajectory to best overlap (in sense of best least square fit) in progressive fashion \\
#                 allows to save transformation matrix per each step

proc geom_center {selection} {
        # set the geometrical center to 0
        set gc [veczero]
        # [$selection get {x y z}] returns a list of {x y z} 
        #    values (one per atoms) so get each term one by one
        foreach coord [$selection get {x y z}] {
           # sum up the coordinates
           set gc [vecadd $gc $coord]
        }
        # and scale by the inverse of the number of atoms
        return [vecscale [expr 1.0 /[$selection num]] $gc]
}

proc alignGeomCenter {molID1 molID2 selText1 selText2} {
  set selX [atomselect $molID1 $selText1]
  set selP [atomselect $molID2 $selText2]
  set allX [atomselect $molID1 all]
  set a [geom_center $selX]
  set b [geom_center $selP]
  set s [vecsub $b $a]
  set M [transoffset $s]
  $allX move $M
}

proc centerSel {selection} {
set s [geom_center $selection]
set M [transoffset $s]
$selection move $M
}

#proc align2XYZ {molID selText {frame 0}} {
#  #select molecules
#  set sel [atomselect $molID $selText frame $frame]
#  set allMol [atomselect $molID "all" frame $frame]
#  #move to {0 0 0}
#  set X [geom_center $sel]
#  $allMol moveby [vecinvert $X]
#  #align the axis
#  set I [draw principalaxes $sel]
#  set A [orient $sel [lindex $I 2] {0 0 1}]
#  $allMol move $A
#  set I [draw principalaxes $sel]
#  set A [orient $sel [lindex $I 1] {0 1 0}]
#  $allMol move $A
#  set I [draw principalaxes $sel]
#  draw delete all
#}

proc alignSelections {sel1 sel2} {
	set M [measure fit $sel2 $sel1]
	$sel2 move $M
}

proc centerByMol { mID centerMolselText outFilename patternFrameNum} {
#mol new $rFileName.pdb type {pdb} first 0 last -1 step 1 waitfor -1
#mol addfile $rFileName.dcd type {dcd} first 0 last -1 step 1 waitfor -1
#set molID [molinfo top get {id}]
set molID $mID
set sel1 [atomselect $molID $centerMolselText]
set sel2 [atomselect $molID $centerMolselText]
set sel3 [atomselect $molID "all"]
$sel1 frame $patternFrameNum
set nf [molinfo $molID get numframes]
for {set i 0} {$i < $nf} {incr i 1} {
	$sel2 frame $i
	set M [measure fit $sel2 $sel1]
    $sel3 frame $i
	$sel3 move $M
}

animate write dcd $outFilename.dcd beg 1 end [expr $nf - 1] skip 1 top
}

proc progressiveFit { mID {matrixFileHandle 0}} {
    set molID $mID
    set sel1 [atomselect $molID "all"]
    set sel2 [atomselect $molID "all"]
    set nf [molinfo $molID get numframes]
    set list_M {}
    for {set i 1} {$i < $nf} {incr i 1} {
        $sel1 frame [expr $i - 1]
	    $sel2 frame $i
	    set M [measure fit $sel2 $sel1]
	    if {$matrixFileHandle != 0} {lappend list_M $M}
	    $sel2 move $M
    }
    #animate write dcd aligned.dcd beg 1 end [expr $nf - 1] skip 1 $mID
    #puts "Aligned trajectory saved as aligned.dcd."
    if {$matrixFileHandle != 0} {
        foreach line $list_M {
	        puts $matrixFileHandle $line
        }
        close $matrixFileHandle
        puts "List of transformation matrixes saved to rot.mat."
    }
}
