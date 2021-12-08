## SASA 
set sel [atomselect top "protein"]
set n [molinfo top get numframes]
set output [open "SASA.dat" w]
# sasa calculation loop
for {set i 0} {$i < $n} {incr i} {
	molinfo top set frame $i
	set sasa [measure sasa 1.4 $sel -restrict $sel]
	puts "\t \t progress: $i/$n"
	puts $output "$i $sasa"
}
puts "\t \t progress: $n/$n"
puts "Done."	
puts "output file: SASA.dat"
close $output  
