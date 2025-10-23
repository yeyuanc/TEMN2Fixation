//==============================================
//$id: run_script_gen.cpp, v1.0
// J.Tang, 2008/04/18
// X.Zhu, 2012/08/26: update due to change of PBS software used by RCAC/Purdue
// B.Zhao 2020/3/4: update due to Rice transition to Slurm
//==============================================
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>


using namespace std;

//geneate PBS shell script

//run_script_gen3 ofname #nparts #njb #ppn queue #walltime xtem tem4.para tem4.log
int main(int argc, char* argv[])
{
	int njb; // # of jobs to submit: 10
	int nps; // # of parts: 480
	int nrun, j; // nrun: # of parts per job: 480/10=48 (example for hansen with 48 cores per node)
	FILE* ofile[20000];
	char strtmp[200];
	int walltime;
	int ppn; // # of cores per node (example for hansen with 48 cores per node)
	char queue[200];

	cout << argc << endl;

	sscanf(argv[3], "%d", &njb);
	sscanf(argv[2], "%d", &nps);
	nrun = nps / njb; // 48
	ppn = atoi(argv[4]); // 48
	strcpy(queue, argv[5]);
	walltime = atoi(argv[6]);

	for(j = 0; j < njb; j ++)
	{
		sprintf(strtmp, "%s-%d",argv[1], j);
		ofile[j] = fopen(strtmp, "w");
		if(!ofile[j])
		{
			cout <<endl<<"Cannot open file "<<strtmp<<" for parameter copy"<<endl;
			exit(-1);
		}
	}

	for(j = 0; j < njb; j ++)
	{
		fprintf(ofile[j], "#!/bin/bash -l\n");
		fprintf(ofile[j], "#SBATCH -t %d:00:00\n", walltime);
		fprintf(ofile[j], "#SBATCH -A %s\n", queue);
		fprintf(ofile[j], "#SBATCH -N %d\n", nrun / ppn);
		//fprintf(ofile[j], "#SBATCH -N %d: -n %d,naccesspolicy=shared myjobsubmissionfile.sub\n", nrun / ppn, ppn); // hansen: 48/48=1, coates: 48/8=6
		fprintf(ofile[j], "#SBATCH -n %d\n", ppn);
		fprintf(ofile[j], "#SBATCH -o\n"); // join PBS output and error file
		fprintf(ofile[j], "#SBATCH --shared\n");
		//fprintf(ofile[j], "#SBATCH --mail-type=END\n");
		fprintf(ofile[j], ". /etc/profile\n");
		//fprintf(ofile[j], "module load mpich2/1.4.1p1_intel-12.1\n");
		//fprintf(ofile[j], "module load mpich2/1.4.1p1_intel-12.0.084\n");
		fprintf(ofile[j], "module load impi/5.1.2.150\n");
		fprintf(ofile[j], "cd $SLURM_SUBMIT_DIR\n");

		fprintf(ofile[j], "nd=`srun hostname | wc -l`\n"); // # of cores to require: nd=nodes*ppn=1*48=6*8=48
		//fprintf(ofile[j], "mpiexec.hydra -np $nd %s %s-%d %s\n", argv[7], argv[8], j, argv[9]);
		fprintf(ofile[j], "srun --mpi=pmi2 -n $nd %s %s-%d %s\n", argv[7], argv[8], j, argv[9]);
		fclose(ofile[j]);
	}
	return 0;

}

