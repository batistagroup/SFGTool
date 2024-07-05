####################################
### generate data
####################################

cd ppp/
	cd azimuthal_anisotropy/
		python ../../../sfg.py sfg.inp
	cd ..

	cd azimuthal_avg/
		python ../../../sfg.py sfg.inp
	cd ..
cd ..
cd ssp/
	python ../../sfg.py sfg.inp
cd ..
cd sps/
	python ../../sfg.py sfg.inp
cd ..
#cd psp/
#	python ../../sfg.py sfg.inp
#cd ..


