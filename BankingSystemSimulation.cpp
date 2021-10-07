#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <gsl_randist.h>
#include <gsl_cdf.h>
#include "common.h"
#include <gsl_errno.h>
#include <gsl_matrix.h>
#include <gsl_rng.h>
#define NARGS 2
#define MAXSTRLEN 256
#define NCOLS 7
void Usage(char* s) {
	fprintf(stderr, "Usage: %s [-nointerbank | ­
		dumpbankloss | -dumpecontrib | -dumpecloss | ­
		randomseq] inputfile LargeLossLimit\n", s);
		exit(-1);
}
double R(double x) {
	return(0.12 * (1 - exp(-50 * x)) / (1 - exp(-50)) + 0.24 *
		(1 - (1 - exp(-50 * x)) / (1 - exp(-50))));
}
double ba(double x) {
	return(pow(0.11852 - 0.05478 * log(x), 2));
}
double regcap(double x, double y, double k) {
	double a;
	a = (k * gsl_cdf_ugaussian_P(
		pow(1 - R(x), -0.5) * gsl_cdf_ugaussian_Pinv(x) +
		pow(R(x) / (1 - R(x)), 0.5) * y) - x * k) *
		pow(1 - 1.5 * ba(x), -1) * 1.06;
	return(a);
}
void printDoubleArray(char* s, int n, double* p) {
	int i;
	if (s != NULL) {
		printf("-->ARRAY %s\n", s);
	}
	for (i = 0; i < n; i++) {
		printf("%f\n", p[i]);
	}
}
int main(int argc, char* argv[]) {
	double rho = 0;
	double LGD = 0.45;
	double BankLossCap = 0;
	char* filename = NULL;
	char* outfnameprefix1 = "Bankloss";
	char* outfnameprefix2 = "DIScharged";
	char* outfnameprefix3 = "Econtrib";
	char* outfnameprefix4 = "EcLoss";
	char* outfname1;
	char* outfname2;
	char* cont = "n";
	FILE* fp = NULL;
	FILE* ofp1 = NULL;
	FILE* ofp2 = NULL;
	FILE* ofp3 = NULL;
	FILE* ofp4 = NULL;
	char buffer[MAXSTRLEN];
	double** data = NULL;
	double** BankLoss = NULL;
	double* RContrib = NULL;
	double* EcLoss = NULL;
	double* DIScharged = NULL;
	double* ECcharged = NULL;
	double* market = NULL;
	double* corr = NULL;
	int* BankDefaulted = NULL;
	int* BankDefaultedOld = NULL;
	double InterbankUnitaryLoss = 0;
	double Baseloss = 0;
	int i, b, j, nDefaults;
	int nbanks = -1;
	int niter = -1;
	int hadDefault;
	int hadDefaultInterbank = 0;
	int numDefaultedBanks = 0;
	int allDefaulted;
	double InterBankLoss, InterBankDefault;
	int argnum = 1;
	int intsize1 = 0;
	int intsize2 = 0;
	int intsize3 = 0;
	int intsize4 = 0;
	int interbanksimulation = 1;
	int dumpbankloss = 0;
	int dumpecontrib = 0;
	int dumpecloss = 0;
	int nTotDefaults = 100000;
	int LargeLossLimit = 0;
	double tmp = 0;
	const gsl_rng_type* Ra;

	gsl_rng* r;

	gsl_rng_env_setup();

	Ra = gsl_rng_default;

	r = gsl_rng_alloc(Ra);

	/* checking number of args */
	if (argc < NARGS + 1) {

		Usage(argv[0]);

		exit(-1);

	}
	while (argnum < argc - NARGS) {
		if (!strcmp(argv[argnum], "-nointerbank")) {
			argnum++;
			interbanksimulation = 0;
			fprintf(stderr, "Interbank Simulation Disabled\n");
		}
		else if (!strcmp(argv[argnum], "-dumpbankloss")) {
			argnum++;
			dumpbankloss = 1;
			fprintf(stderr, "Dumping BankLoss Enabled\n");
		}
		else if (!strcmp(argv[argnum], "-dumpecloss")) {
			argnum++;

			dumpecloss = 1;

			fprintf(stderr, "Dumping Economic Losses

				Enabled\n");
		}
		else if (!strcmp(argv[argnum], "-dumpecontrib")) {
			argnum++;
			dumpecontrib = 1;
			fprintf(stderr, "Dumping Economic Risk
				Contributions Enabled\n");
		}
		else if (!strcmp(argv[argnum], "-randomseq")) {
			argnum++;
			srandom(time(NULL));
			fprintf(stderr, "Init random seed\n");
		}
		else {

			Usage(argv[0]);

		}

	}

	/* reading command line */

	filename = argv[argnum];

	if ((fp = fopen(filename, "r")) == NULL) {

		fprintf(stderr, "could not open file %s for reading
			\n", filename);
			exit(-1);
	}
	LargeLossLimit = atof(argv[argnum + 1]);
	if ((LargeLossLimit < 0)) {
		fprintf(stderr, "LargeLossLimit cannot be
			negative: % s\n", LargeLossLimit);
			exit(-1);
	}
	intsize1 = 10;

	intsize2 = floor(log10(nTotDefaults)) + 1;

	intsize3 = 15;

	outfnameprefix2 = filename;
	/* count number of banks (nrows in data matrix) */

	nbanks = 0;

	while (fgets(buffer, MAXSTRLEN, fp)) {

		nbanks++;

	}

	rewind(fp);

	gsl_matrix* m = gsl_matrix_alloc(1, nbanks);

	gsl_matrix* c = gsl_matrix_alloc(nbanks, nbanks);

	gsl_matrix* out = gsl_matrix_alloc(1, nbanks);

	/* allocate data matrix */

	if ((data = (double**)malloc(nbanks * sizeof

	(double*))) == NULL) {

		fprintf(stderr, "could not allocate data\n");

		exit(-1);

	}
	for (i = 0; i < nbanks; i++) {
		if ((data[i] = (double*)malloc(NCOLS * sizeof
		(double))) == NULL) {
			fprintf(stderr, "could not allocate data[%d]\n", i);
			exit(-1);
		}
		memset(data[i], 0, NCOLS * sizeof(double));
	}
	/* reading input file and initializing data matrix */
	for (i = 0; i < nbanks; i++) {
		212 Appendix: Software References and Tools
			if (fgets(buffer, MAXSTRLEN, fp) == NULL) {

				fprintf(stderr, "could not read line %d

					in input file\n", i);

					exit(-1);
			}
		if (sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf
			\n", &data[i][0], &data[i][1], &data[i][2], &data[i]
			[3], &data[i][4], &data[i][5], &data[i][6]) < NCOLS) {
			fprintf(stderr, "could not read the %d tokens expected
				\n", NCOLS);
				exit(-1);
		}

	}

	fclose(fp);
	/* allocating matrix BankLoss */

	if ((BankLoss = (double**)malloc

	(nTotDefaults * sizeof(double*))) == NULL) {
		fprintf(stderr, "could not allocate BankLoss\n");
		exit(-1);
	}

	for (i = 0; i < nTotDefaults; i++) {

		if ((BankLoss[i] = (double*)malloc

		(nbanks * sizeof(double))) == NULL) {

			fprintf(stderr, "could not allocate BankLoss

				[% d]\n", i);

				exit(-1);
		}
		memset(BankLoss[i], 0, nbanks * sizeof(double));
	}
	/* allocating RContrib */
	if ((RContrib = (double*)malloc(nbanks * sizeof
	(double))) == NULL) {
		fprintf(stderr, "could not allocate RContrib\n");
		exit(-1);
	}

	memset(RContrib, 0, nbanks * sizeof(double));

	/* allocating EcLoss */
	if ((EcLoss = (double*)malloc(nTotDefaults * sizeof
	(double))) == NULL) {
		fprintf(stderr, "could not allocate EcLoss\n");
		exit(-1);
	}

	memset(EcLoss, 0, nTotDefaults * sizeof(double));

	/* allocating DIScharged */

	if ((DIScharged = (double*)malloc

	(nTotDefaults * sizeof(double))) == NULL) {
		fprintf(stderr, "could not allocate DIScharged\n");
		exit(-1);
	}

	memset(DIScharged, 0, nTotDefaults * sizeof(double));

	/* allocating ECcharged */

	if ((ECcharged = (double*)malloc

	(nTotDefaults * sizeof(double))) == NULL) {
		fprintf(stderr, "could not allocate ECcharged\n");
		exit(-1);
	}

	memset(ECcharged, 0, nTotDefaults * sizeof(double));

	/* allocating market */

	if ((market = (double*)malloc(nbanks * sizeof

	(double))) == NULL) {
		fprintf(stderr, "could not allocate market\n");
		exit(-1);
	}

	memset(market, 0, nbanks * sizeof(double));

	/* allocating corr */
	if ((corr = (double*)malloc(nbanks * sizeof
	(double))) == NULL) {
		fprintf(stderr, "could not allocate corr\n");
		exit(-1);
	}
	memset(corr, 0, nbanks * sizeof(double));
	/* allocating BankDefaulted */
	if ((BankDefaulted = (int*)malloc(nbanks * sizeof
	(int))) == NULL) {
		fprintf(stderr, "could not allocate BankDefaulted\n");
		exit(-1);
	}
	memset(BankDefaulted, 0, nbanks * sizeof(int));
	/* allocating BankDefaultedOld */
	if ((BankDefaultedOld = (int*)malloc(nbanks * sizeof
	(int))) == NULL) {
		fprintf(stderr, "could not allocate BankDefaulted
			Old\n");
			exit(-1);
	}
	memset(BankDefaultedOld, 0, nbanks * sizeof(int));
	/* compute correlation matrix and its Cholesky decomp */
	for (i = 0; i < nbanks; i++) {
		corr[i] = pow((rho + data[i][6]), 0.5);
	}
	for (i = 0; i < nbanks; i++) {
		for (j = 0; j < nbanks; j++) {
			gsl_matrix_set(c, i, j, corr[i] * corr[j]);
			if (i == j) { gsl_matrix_set(c, i, j, 1); }
		}
	}

	gsl_linalg_cholesky_decomp(c);

	for (i = 0; i < nbanks; i++) {
		for (j = 0; j < nbanks; j++) {
			if (i < j) { gsl_matrix_set(c, i, j, 0); }
		}
	}
	/* loop until number of defaults equals nTotDefaults */

	nDefaults = 0;

	niter = -1;

	while (nDefaults < nTotDefaults) {

		niter++;

		hadDefault = 0;

		memset(BankDefaulted, 0, nbanks * sizeof(int));
		for (b = 0; b < nbanks; b++) {

			gsl_matrix_set(m, 0, b, gsl_ran_gaussian(r, 1));

		}
		for (j = 0; j < nbanks; j++) {

			tmp = 0;

			for (i = 0; i < nbanks; i++) {

				tmp += gsl_matrix_get(m, 0, i)

					* gsl_matrix_get(c, j, i);

			}
			gsl_matrix_set(out, 0, j, tmp);
		}
		for (i = 0; i < nbanks; i++) {

			market[i] = gsl_matrix_get(out, 0, i);

		}
		for (b = 0; b < nbanks; b++) {

			Baseloss = regcap(data[b][0], market[b], LGD);
			BankLoss[nDefaults][b] = data[b][3] * Baseloss;
			/*	 Check for defaults */
			if (BankLoss[nDefaults][b] > data[b][2]) {

				BankDefaulted[b] = 1;

				hadDefault = 1;

			}

		}

		for (j = 0; j < nbanks; j++) {
			BankDefaultedOld[j] = BankDefaulted[j];
		}
		if (interbanksimulation) {

			cont = "c";

			allDefaulted = 1;

			for (j = 0; j < nbanks; j++) {

				if (BankDefaultedOld[j] == 0) {

					allDefaulted = 0;

					break;

				}

			}

			hadDefaultInterbank = hadDefault;
			while ((allDefaulted == 0) && (hadDefaultInterbank)) {
				hadDefaultInterbank = 0;
				InterBankDefault = 0;
				InterBankLoss = 0;
				for (j = 0; j < nbanks; j++) {
					if (BankDefaulted[j]) {
						InterBankDefault += data
							[j][4];
					}
					InterBankLoss += data[j][5];
				}
				InterbankUnitaryLoss = (InterBankDefault /
					InterBankLoss);
				if ((InterBankDefault / InterBankLoss) > 1)
				{
					InterbankUnitaryLoss = 1;
				}
				for (j = 0; j < nbanks; j++) {
					BankDefaulted[j] = 0;
					BankLoss[nDefaults][j] += (InterbankUnitary
						Loss * data[j][5]);
					if ((BankLoss[nDefaults][j] > data[j][2]) &&
						(BankDefaultedOld[j] == 0)) {
						BankDefaulted[j] = 1;
						hadDefaultInterbank = 1;
						BankDefaultedOld[j] = 1;
					}
				}

				allDefaulted = 1;

				for (j = 0; j < nbanks; j++) {

					if (BankDefaultedOld[j] == 0) {

						allDefaulted = 0;

						break;

					}

				}

			}

		}

		/* end of interbank simulation */

/* update data structures before next iteration */
		for (j = 0; j < nbanks; j++) {
			if (BankDefaultedOld[j]) {
				numDefaultedBanks++;
				DIScharged[nDefaults] += data[j][1];
				ECcharged[nDefaults] += (BankLoss
					[nDefaults][j] - data[j][2]);
			}
		}
		if (hadDefault) {
			nDefaults++;
		}
		/* printf("Iterations=%d, nDefaults=%d/%d,
		nIBDefaults=%d\n", niter+1, nDefaults,
		numDefaultedBanks); */
}
intsize4 = floor(log10(numDefaultedBanks + 1)) + 1;
dump BankLoss : nbanks
dump EContrib : nbanks > LargeLossLimit
dump EcLoss nDefaults nbanks > LargeLossLimit
dump DIScharged nDefaults
/* dump BankLoss */
if (dumpbankloss) {
	outfname1 = malloc(strlen(outfnameprefix1)
		+ strlen("_") + strlen(cont) + strlen("_") + strlen
		(outfnameprefix2) +intsize1 + strlen("_") + intsize2 +
		strlen("-") + intsize4 + strlen(".dat") + 1);
	sprintf(outfname1, "%s_%s_%s_%d-%d-%d.dat",
		outfnameprefix1, outfnameprefix2, cont, LargeLossLimit,
		numDefaultedBanks, niter + 1);
	if ((ofp1 = fopen(outfname1, "w")) == NULL) {
		fprintf(stderr, "could not open %s for
			writing\n", outfname1);
			exit(-1);

	}

	for (i = 0; i < nTotDefaults; i++) {

		fprintf(ofp1, "%1.3lf", BankLoss[i][0]);
		for (b = 1; b < nbanks; b++) {
			fprintf(ofp1, " %1.3lf", BankLoss[i][b]);
		}
		fprintf(ofp1, "\n");
	}

	fclose(ofp1);

}

/* dump EContrib */

if (dumpecontrib) {

	outfname1 = malloc(strlen(outfnameprefix3)

		+ strlen("_") + strlen(outfnameprefix2) + strlen("_")
		+ strlen(cont) + strlen("_") + intsize1 + strlen
		("_") + intsize2 + strlen("-") + intsize4 + strlen
		(".dat") + 1);
	sprintf(outfname1, "%s_%s_%s_%d-%d-%d.dat",
		outfnameprefix3, outfnameprefix2, cont, LargeLossLimit,
		numDefaultedBanks, niter + 1);
	if ((ofp3 = fopen(outfname1, "w")) == NULL) {
		fprintf(stderr, "could not open %s for writing\n",
			outfname1);
		exit(-1);

	}

	for (i = 0; i < nDefaults; i++) {

		if (ECcharged[i] > LargeLossLimit) {
			for (b = 0; b < nbanks; b++) {
				if (BankLoss[i][b] > data[b][2]) {
					RContrib[b] += ((((ECcharged[i] -
						LargeLossLimit) / ECcharged[i]) * (BankLoss[i][b] - data
							[b][2])) / (niter + 1));

				}
			}

		}

		for (i = 0; i < nbanks; i++) {
			fprintf(ofp3, "%f\n", RContrib[i]);

		}

		fclose(ofp3);

	}

	/* dump EcLoss */

	if (dumpecloss) {

		outfname1 = malloc(strlen(outfnameprefix4)

			+ strlen("_") + strlen(outfnameprefix2) + strlen("_")
			+ strlen(cont) + strlen("_") + intsize1 + strlen
			("_") + intsize2 + strlen("-") + intsize4 + strlen(".
				dat") + 1);
				sprintf(outfname1, "%s_%s_%s_%d-%d-%d.dat",
					outfnameprefix4, outfnameprefix2, cont, LargeLossLimit,
					numDefaultedBanks, niter + 1);
		if ((ofp4 = fopen(outfname1, "w")) == NULL) {
			fprintf(stderr, "could not open %s for writing
				\n", outfname1);
				exit(-1);
		}
		for (i = 0; i < nDefaults; i++) {
			if (ECcharged[i] > LargeLossLimit) {
				for (b = 0; b < nbanks; b++) {
					BankLossCap = ((BankLoss[i][b]­
						data[b][2]) < (data[b][1] + data[b][4]) ? (BankLoss
							[i][b] - data[b][2]) : (data[b][1] + data[b][4]));
					if (BankLoss[i][b] > data[b][2]) {
						EcLoss[i] +=
							BankLossCap;
					}

				}

			}

		}

		for (i = 0; i < nDefaults; i++) {
			fprintf(ofp4, "%f\n", EcLoss[i]);

		}

		fclose(ofp4);

	}

	/* dump DIScharged */
	outfname2 = malloc(strlen(outfnameprefix2) + strlen
	("_") + strlen(cont) + strlen("_")
		+ intsize1 + strlen("_") + intsize2 + strlen("-")
		+ intsize4 + strlen(".dat") + 1);
	sprintf(outfname2, "%s_%s_%d-%d-%d.dat",

		outfnameprefix2, cont, LargeLossLimit,

		numDefaultedBanks, niter + 1);

	if ((ofp2 = fopen(outfname2, "w")) == NULL) {
		fprintf(stderr, "could not open %s for writing
			\n", outfname2);
			exit(-1);
	}
	for (i = 0; i < nTotDefaults; i++) {
		fprintf(ofp2, "%f\n", DIScharged[i]);

	}

	fclose(ofp2);

	return 0;
}
