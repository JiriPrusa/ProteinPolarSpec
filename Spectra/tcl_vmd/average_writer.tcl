###################################### PART 1 ###############################################################
#choose file.pdb to load
set molID [mol new jesslys_protein.pdb]
set nof_traj 100
set sel [atomselect $molID "all"]
animate delete beg 0 end -1 skip 0 $molID

# First write an average position for each C-Alpha subselection from each trajectory.dcd
for {set i 0} {$i < 100} {incr i} {
	mol addfile "/DATA/lysozyme_opm_2021/${i}/traj.dcd" type {dcd} first 0 last -1 step 1 waitfor -1 $molID
	puts "${i}"
	set ap [measure avpos $sel]
	$sel set {x y z} $ap
	$sel writepdb ap_all_${i}.pdb
	animate delete beg 0 end -1 skip 0 $molID      		  
}
mol delete $molID


######################################## PART 2 ##########################################################
# Load averaged C-Alpha position from each trajectory and make an avarage from it (The queen is born!)
set molID [mol new ap_0.pdb]
for {set i 1} {$i < 100} {incr i} {
	mol addfile "ap_${i}.pdb" type {pdb} waitfor -1 $molID 
}
set sel [atomselect $molID all]
set ap [measure avpos $sel]
$sel set {x y z} $ap
$sel writepdb queen_CA.pdb
mol delete $molID

quit
