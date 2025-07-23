#!/bin/bash
#SBATCH -J name
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH -n 12

module load compiler/latest
module load compiler-rt/latest
module load mkl/latest
module load mpi/latest
module load tbb/latest


IMAX=10 

echo "{\"status\": \"started\"}" > status.json

I=0
mkdir run.0
cp POSCAR INCAR KPOINTS POTCAR run.0
cd run.0
mpirun vasp >& stdout
NSTEPS=$(cat stdout | grep E0 | wc -l)
cd ..

while [ $NSTEPS -gt 1 ] && [ $I -lt $IMAX ]
do
 I=$(($I+1))
 cp -r run.$(($I-1)) run.$I
 rm run.$(($I-1))/POTCAR
 cd run.$I
 cp CONTCAR POSCAR
 mpirun vasp >& stdout
 NSTEPS=$(cat stdout | grep E0 | wc -l)
 cd ..
done

I=$(($I+1))
cp -r run.$(($I-1)) run.final
rm run.$(($I-1))/POTCAR
cd run.final
cp CONTCAR POSCAR

sed -i "s/.*LREAL.*/LREAL = .FALSE./g" INCAR
sed -i "s/.*IBRION.*/IBRION = -1/g" INCAR
sed -i "s/.*NSW.*/NSW = 0/g" INCAR
sed -i "s/.*ISIF.*/ISIF = 2/g" INCAR
#sed -i "s/.*ISMEAR.*/ISMEAR = -5/g" INCAR
#sed -i "s/.*LCHARG.*/LCHARG = .TRUE./g" INCAR

mpirun vasp >& stdout
cd ..

echo "{\"status\": \"complete\"}" > status.json
