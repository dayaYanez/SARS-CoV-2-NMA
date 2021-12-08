proc gyr_radius {sel} {
  if {[$sel num] <= 0} {
    error "gyr_radius: must have at least one atom in selection"
  }
  set com [center_of_mass $sel]
  set sum 0
  foreach coord [$sel get {x y z}] {
    set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
  }
  return [expr sqrt($sum / ([$sel num] + 0.0))]
}
