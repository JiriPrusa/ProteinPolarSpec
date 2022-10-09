# tcl skript for assingning charge in VMD from B-factor value in pdb for molecule id "molID" (default top).
# USAGE: just load the pdb in VMD as top and source this script
puts "Loaded BfactorAsCharge.tcl"
puts "avail: assignChargeFromBeta {molID}"
proc assignChargeFromBeta {{molID top}} {
	set ev [atomselect top all]
    set indices [$ev get index]
   
    foreach ind $indices {
	    set at [atomselect $molID "index $ind"]
        $at set charge [$at get beta]
#		$at set radius [$at get occupancy]
		$at delete
    } 
$ev delete
}
