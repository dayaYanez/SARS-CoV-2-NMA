set num [expr {$num_steps - 1}]

set outfile [open rmsf.dat w]
set sel [atomselect top " protein and name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
for {set i 0} {$i < [$sel num]} {incr i} {
  puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"
} 
close $outfile
