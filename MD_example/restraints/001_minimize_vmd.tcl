#!/bin/bash 

        mol new ../HRDMYDD_model_1.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
        mol addfile ../HRDMYDD_model_1.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

        set a [atomselect top "all"]
        set b [atomselect top "protein and sidechain and segname PROA or segname PROB"]
	set c [atomselect top "protein and backbone and segname PROA or segname PROB"]
        $a set occupancy 0
        $a set beta 0
        $b set beta 0
	$c set beta 1

        $a writepdb res001_posres.ref

exit
