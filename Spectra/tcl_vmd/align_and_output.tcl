source BfactorAsCharge.tcl
source dipoleWriter.tcl

# this script align trajectories, write down transformation matrix, and
# outputs permanent dipole moment.

# Load twice reference structure pdb with charge in B-factor column 
# First as reference to align to
# Second as topology for trajectories
set ref [mol new reference.pdb]
set molID [mol new reference.pdb]
assignChargeFromBeta $molID
set nof_traj 100
set sel_text "protein and name CA"
# In loop load trajectories align and output data
set ref_ca [atomselect $ref $sel_text]
set sel [atomselect $molID all]
set calphas [atomselect $molID $sel_text]
for {set i 0} {$i < $nof_traj} {incr i} {
    animate delete beg 0 end -1 skip 0 $molID
	mol addfile "${i}/traj.dcd" type {dcd} first 0 last -1 step 1 waitfor -1 $molID
	puts "${i}"
    set fo [open "rot_${i}.mat" w]
    set nf [molinfo $molID get numframes]
    for {set n 0} {$n < $nf} {incr n 1} {
        $ref_ca frame $i
        $sel frame $i
        $calphas frame $i
        set M [measure fit $calphas $ref_ca]
        $sel move $M
        puts $fo $M
    }
    close $fo
    dipoleWriterSel $sel dip_perm_aligned_${i}.csv 0.6    		  
}
mol delete $molID
quit



