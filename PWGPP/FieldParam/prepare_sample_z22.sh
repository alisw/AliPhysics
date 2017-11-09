#!/bin/bash
# cd field-parameterizer
mkdir -p dat_z22

# sample_field filenametag volumespec
function sample_field() {
	tag=$1
	spec=$2
	echo $tag $spec
	# SampleField.C(# of points, kGauss, r_min, p_min, z_min, r_max, p_max, z_max);
	[ -f dat_z22/$tag.sample.dat ] || aliroot -q -b -l "SampleField.C(100000, $spec)" | tail -n +13 > dat_z22/$tag.sample.dat
	[ -f dat_z22/$tag.test.dat   ] || aliroot -q -b -l "SampleField.C(100000, $spec)" | tail -n +13 > dat_z22/$tag.test.dat
}

function sample_field_all() {
	z_names=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
	z_range=(-550 -500 -450 -400 -350 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 350 400 450 500 550)

	r_names=("its" "tpc" "tof" "tofext" "cal")
	r_range=(0     80    250   400      423    500)

	p_range=("0" "TMath::PiOver2()" "TMath::Pi()" "TMath::Pi()*1.5" "TMath::TwoPi()")

	for mag in 2 5; do
		for ri in 0 1 2 3 4; do
			for zi in ${z_names[@]}; do
				for pi in 0 1 2 3; do
					tag=${r_names[ri]}${mag}k-z${z_names[zi]}-q$((pi+1))
					spec="${mag}, ${r_range[ri]}, ${p_range[pi]}, ${z_range[zi]}, ${r_range[ri+1]}, ${p_range[pi+1]}, ${z_range[zi+1]}"
					sample_field $tag "$spec"
				done
			done
		done
	done
}

sample_field_all
