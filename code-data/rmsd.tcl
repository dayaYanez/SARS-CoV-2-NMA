mol new  solvate.psf type psf
mol addfile minCov.dcd type dcd first 0 last -1 step 1 waitfor all
pbc wrap -centersel "protein" -center com -compound residue -all
set reference [atomselect top "protein and backbone" frame 0]
set compare [atomselect top "protein and backbone"]
set num_steps [molinfo top get numframes]
set outfile [open rmsd.dat w]

for {set frame 0} {$frame < $num_steps} {incr frame} {
		$compare frame $frame
		set trans_mat [measure fit $compare $reference]
		$compare move $trans_mat
		set rmsd [measure rmsd $compare $reference ]
   puts $outfile "$frame    $rmsd"
}
close $outfile



