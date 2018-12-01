#eval `ssh-agent -s`
ssh-add
#mpirun -hosts cuda05,cuda06,cuda07,cuda08 -np 14 $PWD/run.sh
mpirun -hostfile hosts -np 14 $PWD/run.sh

