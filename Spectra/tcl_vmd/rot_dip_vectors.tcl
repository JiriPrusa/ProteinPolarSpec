set nof_traj 100
for {set i 0} {$i < $nof_traj} {incr i} {
    # read matrix in
    set fm [open rot_${i}.mat r]
    set rot_matrixes [split [read $fm] "\n"]
    close $fm
    set empty_line_index [expr [llength $rot_matrixes] - 1]
    set rot_matrixes [lreplace $rot_matrixes $empty_line_index $empty_line_index]
    # read dipole in
    set fd [open dip_${i}.csv r]
    set lines [split [read $fd] "\n"]   
    close $fd
    set dipoles [lreplace $lines 0 0]
    set dipoles [lreplace $dipoles $empty_line_index $empty_line_index]
    set fot [open rot_tot_${i}.dip w]
    set foi [open rot_ind_${i}.dip w]
    foreach M $rot_matrixes dip $dipoles {
        set tot_dip [lrange $dip 1 3]
        set ind_dip [lrange $dip 4 6]
        puts $fot [vectrans $M $tot_dip]
        puts $foi [vectrans $M $ind_dip]
    }
    close $fot
    close $foi

    
}
