#!/bin/bash 

        mol new ../HRDMYDD_model_1.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
        mol addfile MYI.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

        set a [atomselect top "all" frame last]
        set b [atomselect top "noh and protein and sidechain and segname PROB" frame last]
	set c [atomselect top "noh and protein and backbone and segname PROB" frame last]
        set d [atomselect top "noh and protein and sidechain and segname PROA" frame last]
	set e [atomselect top "noh and protein and backbone and segname PROA" frame last]

        $a set beta 0
        $b set beta 0
        $c set beta 0
	$d set beta 0
	$e set beta MYJ

        $a writepdb ../restraints/MYP_posres.ref

exit
