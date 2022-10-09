package require pbctools
#Name:
#dipoleWriter
#Synopsis:
#	A Tcl script to write dipole moments of selection from GROMACS trajectory file(s) (*.xtc) to ASCI file.
#Note:
#	Charges are loaded from pdb file, where it must be writen as B-factor.
#Parameters:
#	rFileName - name of file with trajectory and pdb (like rFileName.pdb and rFileName.xtc)
#	aSelection - a VMD formad for residue selection
#	sFileName - name of file to write down the dipole moments
#	step - step (not write dipole moment for every frame)
#Example:
# 	dipoleWriter cys_*_md "water within 5 of resname CYS" cys_*_3ns_ws.dat
##################################################################################################################
#puts "Loaded BfactorAsCharge.tcl"
puts "avail:"
puts "dipoleWriter { rFileName aSelection sFileName {ps_per_step 1} {step 1}}"
puts "angleWriter { rFileName aSelection sFileName {step 1}}"
puts "dipoleWriterSel { selection sFileName {ps_per_step 1} {start 0} {stop -1} {step 1}}"

#Load gromacs trajectory and pdb. Assign charges from beta factor and write the dipole moment components and magnitude
# of selection givne by $aSelection string. Save data to $sFileName.
proc dipoleWriter { rFileName aSelection sFileName {ps_per_step 1} {step 1}} {
set molID [mol new $rFileName.pdb type {pdb} first 0 last -1 step 1 waitfor -1]
mol addfile $rFileName.dcd type {dcd} first 0 last -1 step 1 waitfor -1 $molID
animate delete beg 0 end 0 skip 0 $molID
#read charges from beta-factor column in pdb file (must be writen there before)
assignChargeFromBeta 
set molID [molinfo top get {id}]
set fp [open $sFileName w]
set sel [atomselect $molID $aSelection]
puts $fp "\# BoxVolume: [getVol]        A^3"
puts $fp "\# time    dip_x     dip_y      dip_z    |dip|"
set nf [molinfo $molID get numframes]
	for {set i 0} {$i < $nf} {incr i $step} {
            $sel frame $i
            $sel update
            if {! [catch {measure dipole $sel -debye -masscenter} vector]} {
                puts $fp "[expr $i*$ps_per_step]  [lindex $vector 0]  [lindex $vector 1]  [lindex $vector 2]  [veclength $vector]"
            }
        }
$sel delete
close $fp
mol delete $molID
}


#Load gromacs trajectory and pdb. Assign charges from beta factor and write the dipole moment components, magnitude and 
# angle between dipole moment vector of previous frame from selection givne by $aSelection string. Save data to $sFileName.
proc angleWriter { rFileName aSelection sFileName {step 1}} {

mol new $rFileName.pdb type {pdb} first 0 last -1 step 1 waitfor -1
mol addfile $rFileName.xtc type {xtc} first 0 last -1 step 1 waitfor -1
#read charges from beta-factor column in pdb file (must be writen there before)
assignChargeFromBeta 
set molID [molinfo top get {id}]
animate delete beg 0 end 0 skip 0 $molID
set fp [open $sFileName w]
set sel [atomselect $molID $aSelection]
puts $fp "\# BoxVolume: [getVol]        A^3"
puts $fp "\# frame    dip_x     dip_y      dip_z    |dip|     angle(deg)"
set oldVec {0 0 0}
set nf [molinfo $molID get numframes]
	for {set i 0} {$i < $nf} {incr i $step} {
            $sel frame $i
            $sel update
            if {! [catch {measure dipole $sel -debye -origcenter} vector]} {
				if { $i>0 } {
					set vecAngle [Rad2Deg [expr acos( [vecdot $oldVec $vector]/([veclength $oldVec]*[veclength $vector]) )]]
	                puts $fp "$i  [lindex $vector 0]  [lindex $vector 1]  [lindex $vector 2]  [veclength $vector]  $vecAngle"
				} else {
					puts $fp "$i  [lindex $vector 0]  [lindex $vector 1]  [lindex $vector 2]  [veclength $vector]  0.0 "
				}
            }
			set oldVec $vector
        }
$sel delete
close $fp
mol delete $molID
}

#Write the dipole moment of given $selection from frame $start to $stop with $step to file $sFileName 
proc dipoleWriterSel { selection sFileName {ps_per_step 1} {start 0} {stop -1} {step 1} } {
set molID [$selection molid]
set fp [open $sFileName w]

puts $fp "\# BoxVolume: [getVol]        A^3"
puts $fp "\# time    dip_x     dip_y      dip_z    |dip|"
if {$stop==-1} {
	set nf [molinfo $molID get numframes]
} else {
	set nf $stop
}
	for {set i $start} {$i < $nf} {incr i $step} {
            $selection frame $i
            $selection update
            if {! [catch {measure dipole $selection -debye -masscenter} vector]} {
                puts $fp "[format "%.1f" [expr $i * $ps_per_step]]  [lindex $vector 0]  [lindex $vector 1]  [lindex $vector 2]  [veclength $vector]"
            }
        }
close $fp
puts "Dipole moment of [$selection text] in molecule $molID from frame $start to frame $nf written to $sFileName."
}

proc Deg2Rad {degrees} {expr {$degrees*0.0174532925}}	
proc Rad2Deg {radians} {expr {$radians*57.29577951308}}	

proc cosine {degree} {expr {cos([Deg2Rad $degree])}}

proc unsqrt {alpha beta gama} {expr {sqrt(1+2*$alpha*$beta*$gama-$alpha*$alpha-$beta*$beta-$gama*$gama)}}

proc getVol {} {
   set dimensions [lindex [pbc get -last now] 0]
   set alpha [cosine [lindex $dimensions 3]]
   set beta [cosine [lindex $dimensions 4]]
   set gama [cosine [lindex $dimensions 5]]
   return [expr [lindex $dimensions 0]*[lindex $dimensions 1]*[lindex $dimensions 2]*[unsqrt $alpha $beta $gama]]
   }


