source aligner.tcl
source BfactorAsCharge.tcl

#this script write down transformation matrix to align trajectories with avaraged CA position
#Also it outputs permanent dipole moment for these.

proc decrlist {inlist} {
	foreach element $inlist {
	set new_element [expr $element - 1]
	lappend new_incrlist $new_element
	}
return $new_incrlist
}

set queen_ID [mol new queen_CA.pdb]

set molID [mol new jesslys_prot_charge.pdb]
assignChargeFromBeta $molID
set nof_traj 11
#select index for alignment- this list was selected based on C-alpha that least fluctuated beyond an arbitrary threshold.
set grom_idx {27   43   63   70   94  104  119  138  148  158  168  185  207  231  300  352 373  384  403  410  424  448  464  474  494  551  576  596  610  624 641  651  761  775  787  808  815  834  853 870  889  903  914  938  962 986  996 1010 1088 1123 1137 1156 1166 1180 1201 1213 1223 1234 1244 1263 1282 1293 1304 1316 1335 1349 1359 1370 1386 1400 1410 1420 1442 1464 1483 1859}
set vmd_idx [decrlist $grom_idx]
set queen [atomselect $queen_ID all]
set skeleton [atomselect $molID "index $vmd_idx"]
set sel [atomselect $molID all]
animate delete beg 0 end -1 skip 0 $molID
for {set i 10} {$i < $nof_traj} {incr i} {
	mol addfile "${i}/traj.dcd" type {dcd} first 0 last -1 step 1 waitfor -1 $molID
	puts "${i}"
    # First align first frame according to queen   
    $skeleton frame 0
    $sel frame 0
    set M [measure fit $skeleton $queen]
    $sel move $M
    # And save the transformation matrix
    set fo [open "rot_${i}.mat" w]
    puts $fo $M
    progressiveFit $molID $fo
    dipoleWriterSel $sel dip_perm_aligned_${i}.csv 0.6
	animate delete beg 0 end -1 skip 0 $molID      		  
}
mol delete $molID
mol delete $queen_ID


