function writeClusterFile

for i = 1:66
    subfileName = ['sub_',num2str(i),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH --time 0-12:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP CMake GCCcore/8.3.0 Gurobi/7.0.1\n');
    
    
    
    fprintf(fptr,['a1=',num2str(i),'\n']);
    
    fprintf(fptr,'matlab -nodesktop -singleCompThread -r "FluxEstimationCluster($a1)"');
    
    fclose(fptr);
end
end