#!/bin/bash

dirname='dataset'
pdbIDs=('1TOS' '1UAO' '1KUW' '1ID6' '1DNG' '1O53' '1DU1' '1DPK' '1HO7' '1CKZ' '1LFC' '1A11' '1HO0' '1MMC' '1D0R' '1ZWD' '1D1H' '1SPF' '1AML' '1BA4' '1C56')
chain='A'
datasetSpecifications='dIw_1_aIw_40'

methods=('iabp' 'ibp')

timeD='0'
timeH='01'
timeM='00'
timeS='00'
tolerance='0.01'
angularResolution='1.00'
numberOfSolutions='0'
sampleSize=('3' '5')

mkdir -p "run_files"

for ((k = 0; k < ${#sampleSize[@]}; k++)); do

	folder0="run_files/sample_size=${sampleSize[k]}"
	mkdir -p $folder0
	for ((i = 0; i < ${#pdbIDs[@]}; i++)); do
		pdbID=${pdbIDs[i]}
		instanceFilePath="$dirname/$pdbID/$datasetSpecifications/I_${pdbID}_model1_chain${chain}.dat"
		cliquesFilePath="$dirname/$pdbID/$datasetSpecifications/T_${pdbID}_model1_chain${chain}.dat"
		initialStructureFilePath="$dirname/$pdbID/$datasetSpecifications/X_${pdbID}_model1_chain${chain}.dat"
		for ((j = 0; j < 2; j++)); do
			method=${methods[j]}

			mkdir -p "$folder0/$method"
			mkdir -p "$folder0/$method/$pdbID"
			mkdir -p "$folder0/$method/$pdbID/$datasetSpecifications"

			fileIJK="$folder0/$method/$pdbID/$datasetSpecifications/${pdbID}_${chain}_inputfile.txt"

			echo "structure id: $pdbID" > $fileIJK
			echo "structure chain: $chain" >> $fileIJK
			echo "method: $method" >> $fileIJK
			echo "distance constraints file: $instanceFilePath" >> $fileIJK
			echo "cliques and given torsion angles file: $cliquesFilePath" >> $fileIJK
			echo "reference structure xyz file: $initialStructureFilePath" >> $fileIJK
			echo "time limit (days|hours|minutes|seconds): $timeD-$timeH:$timeM:$timeS" >> $fileIJK
			echo "tolerance (in Ångströms): $tolerance" >> $fileIJK
			echo "angular resolution (in degrees): $angularResolution" >> $fileIJK
			echo "number of solutions (set to 0 for all solutions): $numberOfSolutions" >> $fileIJK
			echo "sample size: ${sampleSize[k]}" >> $fileIJK
		done
	done
done