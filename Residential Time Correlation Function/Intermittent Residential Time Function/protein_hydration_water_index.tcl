set out1 [open hydration_water_index.dat w]

set sys [mol new system_0ns.pdb]
mol addfile analysis.xtc waitfor all

set all [atomselect $sys all]
set p [atomselect $sys protein]

set nf [molinfo top get numframes]

for {set f 0} {$f < $nf} {incr f} {

	molinfo top set frame $f
	$all update
	$p update

	set w [atomselect top "water and within 5.5 of (protein and not hydrogen)"]


	set t [$w get index]

	puts $out1 "$t"
        puts $f
}

close $out1
