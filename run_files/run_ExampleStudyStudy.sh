for f in ../root/Input/CoverageStudy/*.root
do
	for Placement in ../PlacementFiles/*.txt
	do
		Size=$(echo $Placement | sed 's/[^0-9]*//g')
		File=$(echo $f | sed 's/\.\.\//\.\//g')
		Placement=$(echo $Placement | sed 's/\.\.\//\.\//g')
		run="run_"$(echo $f | sed 's/\.root//g' | sed 's/\.\.\/root\/Input\/CoverageStudy\///g' )"_$Size"
		OutputFile=$(echo $run | sed 's/run_/\.\/root\/Output\/CoverageStudy\/lightsim_/g' | sed 's/$/.root/g')
        	echo "cd /pc2014-data3/tdieminger && source setup_lAr.sh && cd /pc2014-data3/tdieminger/SoLAr_Env/qpixls && ./analyze_light $File $Placement $Size $OutputFile" > $run.sh
		chmod +x $run.sh
        	screen -S $run -m -d bash $run.sh
        	sleep 1
done
done
