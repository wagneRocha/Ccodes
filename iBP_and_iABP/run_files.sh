#!/bin/bash

pdbIDs=('1TOS' '1UAO' '1KUW' '1ID6' '1DNG' '1O53' '1DU1' '1DPK' '1HO7' '1CKZ' '1LFC' '1A11' '1HO0' '1MMC' '1D0R' '1ZWD' '1D1H' '1SPF' '1AML' '1BA4' '1C56')
chain='A'
datasetSpecifications='dIw_1_aIw_40'
methods=('iabp' 'ibp')
#sampleSize=('3' '5')
sampleSize=('5')

mainPath='code'

mkdir -p "output"

for ((k = 0; k < ${#sampleSize[@]}; k++)); do
	
	inputFolder="run_files/sample_size=${sampleSize[k]}"
	outputFolder0="output/sample_size=${sampleSize[k]}"
	
	mkdir -p "$outputFolder0"
	
	for ((j = 0; j < 2; j++)); do
		method=${methods[j]}

		for ((i = 0; i < ${#pdbIDs[@]}; i++)); do
			pdbID=${pdbIDs[i]}
		
			mkdir -p "$outputFolder0/$method"
			mkdir -p "$outputFolder0/$method/$pdbID"
			outputFolder="$outputFolder0/$method/$pdbID/$datasetSpecifications"
			mkdir -p $outputFolder

			inputFile="$inputFolder/$method/$pdbID/$datasetSpecifications/${pdbID}_${chain}_inputfile.txt"
			
			echo "method: $method"
			echo "input file: $inputFile"
			echo "output folder: $outputFolder"
			echo "sample size = ${sampleSize[k]}"

			"$mainPath/./main" "$inputFile" "$outputFolder" > "$outputFolder/out.txt"

			echo 
			echo 
		done
	done
done