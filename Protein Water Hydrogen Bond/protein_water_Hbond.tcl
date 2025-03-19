set out1 [open donor.dat w]
set out2 [open acceptor.dat w]
set out3 [open donor_hydrogen.dat w]


set sys [mol new system_0ns.pdb]
mol addfile analysis.xtc waitfor all

set all [atomselect $sys all]
set p [atomselect $sys protein]

set nf [molinfo top get numframes]

for {set f 0} {$f < $nf} {incr f} {

	molinfo top set frame $f
	$all update

	set w [atomselect top "water and within 5.5 of (protein and not hydrogen)"]

	set c [measure hbonds 3.5 20 $p $w]
	set c1 [measure hbonds 3.5 20 $w $p]

	set t [lindex $c 0]
	set t1 [lindex $c1 0]
	set t2 [concat $t $t1]

	set t3 [lindex $c 1]
        set t4 [lindex $c1 1]
        set t5 [concat $t3 $t4]

	set t6 [lindex $c 2]
        set t7 [lindex $c1 2]
        set t8 [concat $t6 $t7]

	puts $out1 "$t2"
	puts $out2 "$t5"
	puts $out3 "$t8"
        puts "$f"

}

close $out1
close $out2
close $out3
exit
