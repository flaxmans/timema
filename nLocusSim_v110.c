/*
 *
 *  Created by Samuel Melvin Flaxman on 6/26/14.
 *  Copyright 2014 Samuel Melvin Flaxman. All rights reserved.
 *
 */

const char *version = "1.1.0";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#include "../MT/dSFMT.h"

// code for using Mersenne Twister RNG
dsfmt_t dsfmt;
#define seedRand(s) dsfmt_init_gen_rand(&dsfmt, s)
#define	randU() dsfmt_genrand_close_open(&dsfmt)
#define	randI() (unsigned)dsfmt_genrand_uint32(&dsfmt)


// parameters
#define MAX_GENERATIONS_DEFAULT 10000
#define ROWS_OF_PATCHES_DEFAULT 10
#define COLS_OF_PATCHES_DEFAULT 2
#define PROPORTION_ADENOSTOMA_DEFAULT 0.5
#define TOTAL_N_DEFAULT 2000
#define INIT_MEL_FREQ_CEANOTHUS_DEFAULT 0.5
#define INIT_MEL_FREQ_ADENOSTOMA_DEFAULT 0.5
#define nGENS_ALLOPATRY_DEFAULT 0
#define SD_MOVE_MELANISTIC_DEFAULT 0.02
#define SD_MOVE_GREEN_DEFAULT 0.016
#define RANDOM_LOCS_FOR_OFFSPRING_DEFAULT 1
// migration probability ends up being about 0.5 or 0.6 * SD_MOVE
#define COLOR_PATTERN_RECOMB_PROBABILITY_DEFAULT 0.5
#define DETERMINISTIC_DEFAULT 1
#define TIME_SAMPLING_INTERVAL_DEFAULT 50
#define PRINT_HEADERS_DEFAULT 1
#define S_COEFF_WRONG_STRIPE_DEFAULT 0.3
#define S_COEFF_MELANISTIC_IN_GREEN_DEFAULT 0.3
#define S_COEFF_GREEN_IN_DARK_DEFAULT 0.3
#define PROPORTION_HOST_DARK_NICHE_DEFAULT 0.1
#define MIGRATION_DEATH_PROBABILITY_DEFAULT 0.0
#define STRIPE_DOMINANCE_COEFFICIENT_DEFAULT 1.0 // dominance of UNSTRIPED_ALLELE over STRIPED ALLELE; 0.5 = co-dominance; > 0.5 = UNSTRIPE IS MORE DOMINANT
#define BASE_MATING_FITNESS_DEFAULT 0.3
#define CHC_MATCH_ADVANTAGE_DEFAULT 0.3
#define MELANISTIC_MATING_ADVANTAGE_M_DEFAULT 0.3
#define MELANISTIC_MATING_ADVANTAGE_F_DEFAULT 0.3
#define MULTIRUN_DEFAULT 0
#define LAYOUT_OF_HOSTS_RANDOM_DEFAULT 0
#define HOLD_MEL_ALLELE_FREQ_CONSTANT_DEFAULT 0
#define BIASED_NICHE_MATING_RATIO_DEFAULT 1.0
#define MAX_NUM_CHOICE_TRIALS_DEFAULT 5
#define ADDITIVE_MODEL_DEFAULT 0



// ALL OF THE FOLLOWING ARE NOT PARAMETERS BUT JUST FLAGS SO THAT SEEMINGLY ARBITRARY CONSTANTS DO NOT
// DOMINATE THE CODE
// locus type flags
// these are NOT the indexes of said loci, just flags
#define NEUTRAL_LOCUS 0
#define COLOR_LOCUS 1
#define PATTERN_LOCUS 2
#define SEX_LOCUS 3
#define ASSORTMENT_LOCUS 4

// allele flags
#define DIPLOID 2
// color
#define GREEN_ALLELE 1 // dominant
#define MELANISTIC_ALLELE 0 //recessive but when expressed, pattern does not matter
// pattern
#define UNSTRIPED_ALLELE 1 // dominant
#define STRIPED_ALLELE 0 // recessive
// neutral
#define NEUT_ADENO_ALLELE 0
#define NEUT_CEANO_ALLELE 1
// sex
#define FEMALE_ALLELE 0 // like an X chromosome
#define MALE_ALLELE 1 // like a Y chromosomes
// assortment
#define NO_ASSORT_ALLELE 0
#define ASSORT_ALLELE 1

// host flags
#define ADENOSTOMA_HOST 0
#define CEANOTHUS_HOST 1

// niche flags
#define NICHE_GREEN 0
#define NICHE_DARK 1

// phenotype flags in default model
#define PH_MELANISTIC 0
#define PH_STRIPED_GREEN_HOMO 1
#define PH_UNSTRIPED_GREEN_HOMO 2
#define PH_HETERO_STRIPE 3

// phenotype flags in purely additive model
#define PH_GU 0 // double homozygous green unstriped;
#define PH_GH 1 // homozygous green, heterozygous stripe;
#define PH_GS 2 // double homozygous green striped;
#define PH_HU 3 // heterozygous color; homozygous unstriped;
#define PH_HH 4 // heterozygous color; heterozygous stripe;
#define PH_HS 5 // heterozygous color; homozygous striped;
#define PH_MU 6 // double homozygous melanistic unstriped;
#define PH_MH 7 // homozygous melanistic, heterozygous stripe;
#define PH_MS 8 // double homozygous melanistic striped;

// flags for sex of individual
#define SEX_FEMALE 0
#define SEX_MALE 1



// function declarations
double boxMuller(double mu, double sd);
double calculateLD(double A, double B, int locus1, int locus2);
int calculateMetrics(double *FST);
double calculateProportionInEachNiche(void);
//void chooseParents1( int *momPt, int *dadPt, int *maleIndexes, int *femaleIndexes, int maleCount, int femaleCount );
void chooseParents2( int *momPt, int *dadPt, int *maleIndexes, int *femaleIndexes, int maleCount, int femaleCount );
//double computeFemaleFitness(int dadi, int momi);
int computePhenotype(int *ipt);
int computeSex(int *ipt);
void finalPrintsMultirun(FILE *mrfpt, int mrcount, int rngseed, double IHmatFreq, double mr, double *FST);
void finalPrintsAllRuns(int rngseed);
int findNiche(double x, double y);
void getOffspringAlleles(int parentIndex, int *offGtPt);
void initializeGenomeScheme(void);
void initializePatches(void);
void initializePopulation(void);
void makeViabilityScoeffMatrix(void);
double migration(void);
int multirunSetup(int seed);
void openDataRecordingFiles(void);
void printPopulation(char *fname);
double reproduction(void);
int RNGsetup(void);
void setMelanisticAlleleFrequency(void);
void sortPopulation(void);
void usage(char *s);
void viabilitySelection(void);


// global variables
int nPATCHES, nLOCI;
int ROWS_OF_PATCHES = ROWS_OF_PATCHES_DEFAULT;
int COLS_OF_PATCHES = COLS_OF_PATCHES_DEFAULT;
double PROPORTION_ADENOSTOMA = PROPORTION_ADENOSTOMA_DEFAULT;
int TOTAL_N = TOTAL_N_DEFAULT;
double SD_MOVE_MELANISTIC = SD_MOVE_MELANISTIC_DEFAULT;
double SD_MOVE_GREEN = SD_MOVE_GREEN_DEFAULT;
_Bool PRINT_HEADERS = PRINT_HEADERS_DEFAULT;
double INIT_MEL_FREQ_ADENOSTOMA = INIT_MEL_FREQ_ADENOSTOMA_DEFAULT;
double INIT_MEL_FREQ_CEANOTHUS = INIT_MEL_FREQ_CEANOTHUS_DEFAULT;
long int nGENS_ALLOPATRY = nGENS_ALLOPATRY_DEFAULT;
long int MAX_GENERATIONS = MAX_GENERATIONS_DEFAULT;
int DETERMINISTIC = DETERMINISTIC_DEFAULT;
int TIME_SAMPLING_INTERVAL = TIME_SAMPLING_INTERVAL_DEFAULT;
double S_COEFF_GREEN_IN_DARK_NICHE = S_COEFF_GREEN_IN_DARK_DEFAULT;
double S_COEFF_MELANISTIC_IN_GREEN_NICHE = S_COEFF_MELANISTIC_IN_GREEN_DEFAULT;
double S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE = S_COEFF_WRONG_STRIPE_DEFAULT;
double PROPORTION_HOST_DARK_NICHE = PROPORTION_HOST_DARK_NICHE_DEFAULT;
double darkNicheLowerBound, darkNichUpperBound;
double MIGRATION_DEATH_PROBABILITY = MIGRATION_DEATH_PROBABILITY_DEFAULT;
double STRIPE_DOMINANCE_COEFFICIENT = STRIPE_DOMINANCE_COEFFICIENT_DEFAULT;
int INDEX_OF_COLOR_LOCUS;
int INDEX_OF_PATTERN_LOCUS;
int INDEX_OF_SEX_LOCUS;
double BASE_MATING_FITNESS = BASE_MATING_FITNESS_DEFAULT;
double CHC_MATCH_ADVANTAGE = CHC_MATCH_ADVANTAGE_DEFAULT;
double MELANISTIC_MATING_ADVANTAGE_M = MELANISTIC_MATING_ADVANTAGE_M_DEFAULT;
double MELANISTIC_MATING_ADVANTAGE_F = MELANISTIC_MATING_ADVANTAGE_F_DEFAULT;
int RANDOM_LOCS_FOR_OFFSPRING = RANDOM_LOCS_FOR_OFFSPRING_DEFAULT;
const double TWOPI = 2.0 * M_PI;
long int t;
double COLOR_PATTERN_RECOMB_PROBABILITY = COLOR_PATTERN_RECOMB_PROBABILITY_DEFAULT;
//double MAX_FITNESS_INV;
int MULTIRUN = MULTIRUN_DEFAULT;
long int dadTries = 0, momTries = 0, reproTries = 0;
int LAYOUT_OF_HOSTS_RANDOM = LAYOUT_OF_HOSTS_RANDOM_DEFAULT;
int HOLD_MEL_ALLELE_FREQ_CONSTANT = HOLD_MEL_ALLELE_FREQ_CONSTANT_DEFAULT;
double BIASED_NICHE_MATING_RATIO = BIASED_NICHE_MATING_RATIO_DEFAULT;
int MAX_NUM_CHOICE_TRIALS = MAX_NUM_CHOICE_TRIALS_DEFAULT;
int ADDITIVE_MODEL = ADDITIVE_MODEL_DEFAULT;
double viabilityMatrix[9][2][2]; // for ADDITIVE_MODEL only

// global pointers for malloc calls
int *genotypes, *phenotypes, *patchLocations, *hostType, *nInEachPatch, *markedDead;
int *nicheWithinPatch, *locusTypes, *sexOfIndividual, *chcType;
double *coordinates, *recombinationProbabilities, *allele_frequencies;
int *locusStillVariable, *fixedAllele;
long int *fixationTimes;

// global data file pointers
FILE *migrationRecord, *patchAbundances, *FSTtimeSeries, *AlFreqByPatchTS;
FILE *AlFreqTS, *LDTS, *interHostMating;



int main(int argc, char *argv[])
{
    int ch, numFixed, i, seed, mrcount;
    char *progname = argv[0];
    double mr, crossMatingFreq;
    FILE *mrfpt;
    
    // read in optional command line arguments
    while ((ch = getopt(argc, argv, "a:A:b:B:c:C:d:D:f:F:g:G:h:H:L:m:M:N:o:p:P:r:R:s:S:t:T:w:z:Z:?")) != -1) {
        switch (ch) {
            case 'a':
                PROPORTION_ADENOSTOMA = strtod(optarg, (char **)NULL);
                break;
            case 'A':
                nGENS_ALLOPATRY = atoi(optarg);
                break;
            case 'b':
                HOLD_MEL_ALLELE_FREQ_CONSTANT = atoi(optarg);
                break;
            case 'B':
                BIASED_NICHE_MATING_RATIO = strtod(optarg, (char **)NULL);
                break;
            case 'c':
                MIGRATION_DEATH_PROBABILITY = strtod(optarg, (char **)NULL);
                break;
            case 'C':
                COLS_OF_PATCHES = atoi(optarg);
                break;
            case 'd':
                PROPORTION_HOST_DARK_NICHE = strtod(optarg, (char **)NULL);
                break;
            case 'D':
                DETERMINISTIC = atoi(optarg);
                break;
            case 'f':
                BASE_MATING_FITNESS = strtod(optarg, (char **)NULL);
                break;
            case 'F':
                CHC_MATCH_ADVANTAGE = strtod(optarg, (char **)NULL);
                break;
            case 'g':
                MELANISTIC_MATING_ADVANTAGE_M = strtod(optarg, (char **)NULL);
                break;
            case 'G':
                MELANISTIC_MATING_ADVANTAGE_F = strtod(optarg, (char **)NULL);
                break;
            case 'h':
                PRINT_HEADERS = atoi(optarg);
                break;
            case 'H':
                STRIPE_DOMINANCE_COEFFICIENT = strtod(optarg, (char **)NULL);
                break;
            case 'L':
                LAYOUT_OF_HOSTS_RANDOM = atoi(optarg);
                break;
            case 'm':
                SD_MOVE_MELANISTIC = strtod(optarg, (char **)NULL);
                break;
            case 'M':
                SD_MOVE_GREEN = strtod(optarg, (char **)NULL);
                break;
            case 'N':
                TOTAL_N = atoi(optarg);
                break;
            case 'o':
                RANDOM_LOCS_FOR_OFFSPRING = atoi(optarg);
                break;
            case 'p':
                INIT_MEL_FREQ_ADENOSTOMA = strtod(optarg, (char **)NULL);
                break;
            case 'P':
                INIT_MEL_FREQ_CEANOTHUS = strtod(optarg, (char **)NULL);
                break;
            case 'r':
                COLOR_PATTERN_RECOMB_PROBABILITY = strtod(optarg, (char **)NULL);
                break;
            case 'R':
                ROWS_OF_PATCHES = atoi(optarg);
                break;
            case 's':
                S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE = strtod(optarg, (char **)NULL);
                break;
            case 'S':
                S_COEFF_MELANISTIC_IN_GREEN_NICHE = strtod(optarg, (char **)NULL);
                break;
            case 't':
                S_COEFF_GREEN_IN_DARK_NICHE = strtod(optarg, (char **)NULL);
                break;
            case 'T':
                MAX_GENERATIONS = atoi(optarg);
                break;
            case 'w':
                MULTIRUN = atoi(optarg);
                break;
            case 'z':
                ADDITIVE_MODEL = atoi(optarg);
                break;
            case 'Z':
                TIME_SAMPLING_INTERVAL = atoi(optarg);
                break;
            case '?':
            default:
                usage(progname);
                exit(-1);
        }
    }
    
    // set up
    seed = RNGsetup();
    initializeGenomeScheme();
    double FST[nLOCI];
    initializePatches();
    initializePopulation();
    if ( ADDITIVE_MODEL )
        makeViabilityScoeffMatrix();
    if ( MULTIRUN ) {
        mrcount = multirunSetup(seed);
        mrfpt = fopen("MultirunEndpointData.csv","a");
    }
    else
        openDataRecordingFiles();
    
    // main operations
    t = 1;
    numFixed = 0;
    while ( numFixed < (nLOCI - 1) && t <= MAX_GENERATIONS ) { /* criterion is (nLOCI-1) b/c sex locus should NEVER fix */
        
        if ( t > nGENS_ALLOPATRY )
            mr = migration();
        
        viabilitySelection();
        
        crossMatingFreq = reproduction(); // includes differential reproduction
        
        numFixed = calculateMetrics( &FST[0] );
        
        t++;
    }
    
    if ( !MULTIRUN && ((t-1) % TIME_SAMPLING_INTERVAL != 0) ) {
        // grab the last set of data if they weren't written automatically elsewhere
        fprintf(migrationRecord, "%li,%f\n", t-1, mr);
        fprintf( patchAbundances, "%li", t-1);
        for ( i = 0; i < nPATCHES; i++ ) {
            fprintf(patchAbundances, ",%i", nInEachPatch[i]);
        }
        fprintf(patchAbundances, "\n");
        fprintf(interHostMating, "%li,%f\n", (t-1), crossMatingFreq);
    }
    
//    fprintf(stderr, "\nReproduction call counts (performance info):\ndadTries = %li, momTries = %li, reproTries = %li", dadTries, momTries, reproTries);
//    fprintf(stderr, "\nAverages: dad tries per repro = %f; mom tries per repro = %f\n\n", ((double)dadTries)/((double)reproTries), ((double)momTries)/((double)reproTries));
    
    if ( MULTIRUN )
        finalPrintsMultirun(mrfpt, mrcount, seed, crossMatingFreq, mr, &FST[0]);
    
    finalPrintsAllRuns(seed);
    
    // free blocks from malloc
    free(hostType);
    free(genotypes);
    free(phenotypes);
    free(patchLocations);
    free(nicheWithinPatch);
    free(markedDead);
    free(locusStillVariable);
    free(fixedAllele);
    free(fixationTimes);
    free(allele_frequencies);
    free(locusTypes);
    free(recombinationProbabilities);
    free(sexOfIndividual);
    
    // close files
    
    if ( MULTIRUN )
        fclose(mrfpt);
    else {
        fclose(migrationRecord);
        fclose(patchAbundances);
        fclose(FSTtimeSeries);
        fclose(AlFreqByPatchTS);
        fclose(AlFreqTS);
        fclose(LDTS);
        fclose(interHostMating);
    }
    return 0;
}


double boxMuller(double mu, double sd)
{
    /* boxmuller.c           Implements the Polar form of the Box-Muller
	 Transformation
	 
	 (c) Copyright 1994, Everett F. Carter Jr.
	 Permission is granted by the author to use
	 this software for any application provided this
	 copyright notice is preserved.
	 [ smf's note: accessed online 12.03.06 at
	 http://www.taygeta.com/random/boxmuller.html ]
	 
	 */
    
    double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;
	
	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ( randU() ) - 1.0;
			x2 = 2.0 * ( randU() ) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
		
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	
	return ( mu + y1 * sd );
}


double calculateLD(double A, double B, int locus1, int locus2)
{
	int i, j;
	double zz = 0.0, oo = 0.0, zo = 0.0, oz = 0.0;
	double Ninv = 0.5 / ((double) TOTAL_N);
	double exp11, DD, Dmax, Dprime, Delta;
	int *gpt1, *gpt2;
    
    
	for ( j = 0; j < DIPLOID; j++ ) {
		// accounts for both diploid sets, maternally and paternally inherited
		gpt1 = genotypes + j + (DIPLOID * locus1); // allele of first locus on chromosome
		gpt2 = genotypes + j + (DIPLOID * locus2); // allele on second locus on chromosome
		
		for ( i = 0; i < TOTAL_N; i++ ) {
			if ( *gpt1 ) {
				if ( *gpt2 ) {
					oo++;
				}
				else {
					oz++;
				}
			}
			else {
				if ( *gpt2 ) {
					zo++;
				}
				else {
					zz++;
				}
			}
			gpt1 += (DIPLOID * nLOCI);
			gpt2 += (DIPLOID * nLOCI);
		}
		
	}
	
	// haplotype frequencies
	oo *= Ninv; // observed 11
	oz *= Ninv; // observed 10
	zo *= Ninv; // observed 01
	zz *= Ninv; // observed 00
	
	DD = (( oo * zz ) - ( zo * oz )); // the textbook "D", ranges between -0.25 and 0.25
    Delta = DD / sqrt( A * B * (1.0 - A) * (1.0 - B) ); // the LD correlation coefficient, ranges between -1 and 1
	
    
    
    exp11 = A * B; // expected frequency of 11 haplotype
    
    // check on calcs
    if ( fabs(DD - (oo - exp11)) > 0.0001 )
    	fprintf(stderr, "Warning in calculateLD(): bad math, DD = %G, diff = %G\n", DD, (oo - exp11));
    if ( DD > 0.0 )
        Dmax = fmin( (A * (1.0 - B)), ((1.0 - A) * B) );
    else
        Dmax = fmin( (A * B), ((1.0 - A) * (1.0 - B)) );
    
    Dprime = DD / Dmax; // Dprime also ranges between -1 and 1
    
    if ( !MULTIRUN )
        fprintf(LDTS, ",%E", Delta);
    
    return Delta;
}



int calculateMetrics(double *FST)
{
    int i, j, patch, locus, numFixed;
	double totalAlleles[nLOCI], allelesByPatch[nLOCI][nPATCHES], Ninv, patchWeights[nPATCHES];
    double alleleFrequenciesByPatch[nLOCI][nPATCHES];
	double gtsum, HSsum, HT, p, Nsinv[nPATCHES], foo;
    int *spt, *stpt, ivar;
	
    Ninv = 1.0 / ((double) TOTAL_N);
    numFixed = 0;
    
    // initialize arrays
	for ( i = 0; i < nLOCI; i++ ) {
		totalAlleles[i] = 0.0;
		FST[i] = 0.0;
		for ( j = 0; j < nPATCHES; j++ ) {
			allelesByPatch[i][j] = 0.0;
		}
	}
	
	for ( i = 0; i < nPATCHES; i++ ) {
		patchWeights[i] = ((double) nInEachPatch[i]) * Ninv;
		Nsinv[i] = 0.5 / ((double) nInEachPatch[i]); // 0.5 since each has two alleles for each locus; reduces calc's below
	}
	
	for ( i = 0; i < TOTAL_N; i++ ) {
		patch = patchLocations[i];
		stpt = genotypes + ( i * nLOCI * DIPLOID );
		for ( locus = 0; locus < nLOCI; locus++ ) {
            spt = stpt + ( DIPLOID * locus );
            gtsum = ((double) ((*spt) + (*(spt+1))));
            totalAlleles[locus] = totalAlleles[locus] + gtsum;
            allelesByPatch[locus][patch] = allelesByPatch[locus][patch] + gtsum;
		}
	}
	
	
	// now calculate allele frequencies and FST
	for ( locus = 0; locus < nLOCI; locus++ ) {
		p = totalAlleles[locus] * Ninv * 0.5;
        allele_frequencies[locus] = p;
        if ( p <= 0.0 || p >= 1.0 ) { // inequalities here to account for possibility of imprecision in discrete math
            
            if ( locus == INDEX_OF_SEX_LOCUS ) {
                fprintf(stderr, "\nError! Sex locus fixed!\n");
                exit(-1);
            }
            
            FST[locus] = 0.0;
            numFixed++;
            if ( p >= 1.0 ) {
                for ( j = 0; j < nPATCHES; j++ ) {
                    alleleFrequenciesByPatch[locus][j] = 1.0;
                }
            }
            else {
                for ( j = 0; j < nPATCHES; j++ ) {
                    alleleFrequenciesByPatch[locus][j] = 0.0;
                }
            }

            if ( locusStillVariable[locus] ) {
                locusStillVariable[locus] = 0;
                fixationTimes[locus] = t;
                if ( p >= 1.0 )
                    fixedAllele[locus] = 1;
                else
                    fixedAllele[locus] = 0;
                fprintf(stderr, "\nLocus %i fixed for allele %i after %li generations\n", locus, fixedAllele[locus], fixationTimes[locus]);
            }
        }
        else {
            HT = 2.0 * p * (1.0 - p);
            HSsum = 0.0;
            for ( j = 0; j < nPATCHES; j++ ) {
                p = allelesByPatch[locus][j] * Nsinv[j];
                alleleFrequenciesByPatch[locus][j] = p;
                // /* // test check
                if ( p < 0.0 || p > 1.0 ) {
                    fprintf(stderr, "\nERROR!, p = %E\n",p);
                }
                // */ // end test check
                HSsum += ( 2.0 * p * (1.0 - p) ) * patchWeights[j];
            }
            FST[locus] = ( HT - HSsum ) / HT;
            // /* // test check
            if ( FST[locus] > 1.0 || FST[locus] < -1.0E-13 ) {
                fprintf(stderr, "\n\nFST ERROR!!\nFST[%i] = %E, HT = %E, HSsum = %E\n", locus, FST[locus],HT, HSsum);
                exit(1);
            }
            if ( FST[locus] < 0.0 )
                FST[locus] = 0.0;
            // */ // end test check
        }
    }
    
    if ( (t % TIME_SAMPLING_INTERVAL == 0 || numFixed == (nLOCI-1) ) && !MULTIRUN ) {
        fprintf(FSTtimeSeries, "%li", t);
        fprintf(AlFreqTS, "%li", t);
        fprintf(AlFreqByPatchTS, "%li", t);
        for ( i = 0; i < nLOCI; i++ ) {
            if ( locusStillVariable[i] )
                fprintf(FSTtimeSeries, ",%E", FST[i]);
            else
                fprintf(FSTtimeSeries, ",NA");
            fprintf(AlFreqTS, ",%E", allele_frequencies[i]);
            for ( j = 0; j < nPATCHES; j++ ) {
                fprintf(AlFreqByPatchTS, ",%E", alleleFrequenciesByPatch[i][j]);
            }
        }
        fprintf(FSTtimeSeries, "\n");
        fprintf(AlFreqTS, "\n");
        fprintf(AlFreqByPatchTS, "\n");
        
        fprintf(LDTS, "%li", t);
        for ( i = 0; i < (nLOCI-1); i++ ) {
            ivar = locusStillVariable[i];
            for ( j = i+1; j < nLOCI; j++ ) {
                if ( ivar && locusStillVariable[j] ) {
                    foo = calculateLD(allele_frequencies[i], allele_frequencies[j], i, j);
                }
                else {
                    fprintf(LDTS, ",NA");
                }
            }
        }
        fprintf(LDTS, "\n");
    }
    
    return numFixed;
}

double calculateProportionInEachNiche(void)
{
    int i, nDark = 0, nGreen = 0;
    double count = 0.0, p;
    
    for ( i = 0; i < TOTAL_N; i++ ) {
        if ( nicheWithinPatch[i] == NICHE_DARK ) {
            nDark++;
            count += 1.0;
        }
        else if ( nicheWithinPatch[i] == NICHE_GREEN )
            nGreen++;
    }
    
    if ( (nDark + nGreen) != TOTAL_N ) {
        fprintf(stderr, "\nError in calculateProportionInEachNiche():\n\tnDark = %i, nGreen = %i, sum = %i != TOTAL_N (%i)\n", nDark, nGreen, (nDark + nGreen), TOTAL_N);
    }
    else {
        p = count / ((double) TOTAL_N);
        fprintf(stderr, "\t%li\t%f\n", t, p);
        
    }
    
    return p;
}


//void chooseParents1( int *momPt, int *dadPt, int *maleIndexes, int *femaleIndexes, int maleCount, int femaleCount )
//{
//    // this choose function is for male choice; male phenotype does not influence his likelihood of mating (given that he has already survived viability selection)
//    // only female phenotype influences the chance of a mating occurring, via male choice
//    
//    int momi, dadi, dumi, dadWarningCount, dadNiche, momWarningCount;
//    _Bool seekingMom;
//    int triesThisTime;
//    double matingReduction, paccept;
//    
//    // choose dad
//    reproTries++;
//    triesThisTime = 0;
//    do {
//        dadTries++;
//        triesThisTime++;
//        dumi = randI() % maleCount;
//        dadi = maleIndexes[dumi];
//    } while ( markedDead[dadi] && triesThisTime < 100 );
//    
//    if ( markedDead[dadi] )
//        dadWarningCount++;
//    
//    dadNiche = nicheWithinPatch[dadi];
//    
//    // choose female
//    seekingMom = 1;
//    //            tries = 0;
//    //            displayOnce = 1;
//    
//    
//    triesThisTime = 0;
//    do {
//        momTries++;
//        triesThisTime++;
//        //                if ( tries > 10 && displayOnce ) {
//        //                    fprintf(stderr, "\nWarning in reproduction() finding momi:\n\ttries = %i, last paccept = %f!\n\tSimulation likely to be SLOW if this keeps happening!\n\tTry recoding with max number possible tries?\n", tries, paccept);
//        //                    displayOnce = 0;
//        //                }
//        
//        dumi = randI() % femaleCount;
//        momi = femaleIndexes[dumi];
//        
//        if ( nicheWithinPatch[momi] != dadNiche )
//            matingReduction = BIASED_NICHE_MATING_RATIO;
//        else
//            matingReduction = 1.0;
//        
//        if ( !markedDead[momi] ) {
//            paccept = (computeFemaleFitness(dadi, momi)) * matingReduction;
//            if ( randU() < paccept )
//                seekingMom = 0;
//        }
//    } while ( seekingMom && triesThisTime < 100 );
//    
//    if ( markedDead[momi] )
//        momWarningCount++;
//    
//    *momPt = momi;
//    *dadPt = dadi;
//
//}


void chooseParents2( int *momPt, int *dadPt, int *maleIndexes, int *femaleIndexes, int maleCount, int femaleCount )
{
    // in this incarnation,
    // there is a bit of what could be thought of as mutual mate choice because the
    // probability of a successful mating can be enhanced by the male being melanistic.
    // There is still male choice, but it is limited to a fixed number of trials
    // set by MAX_NUM_CHOICE_TRIALS
    
    int dumi, dadi, momi, nicheTries, MAXTRIES = 100;
    int choiceTrials, functionTries, triesHere;
    double maleMatingProb, paccept;
    _Bool seekingMom, proceed;
    
    
    
    functionTries = 0;
    do {
        functionTries++;
        reproTries++;
        
        triesHere = 0;
        do {
            
            // choose dad
            do {
                // just get a male who is alive
                dadTries++;
                triesHere++;
                dumi = randI() % maleCount;
                dadi = *(maleIndexes + dumi);
            } while ( ( *(markedDead + dadi) == 1 ) && (triesHere < MAXTRIES) );
            
            // give advantage to MEL. males
            if ( ADDITIVE_MODEL ) {
                if ( ((*(phenotypes + dadi)) == PH_MH) || ((*(phenotypes + dadi)) == PH_MS) || ((*(phenotypes + dadi)) == PH_MU) )
                    maleMatingProb = 1.0;
                else if ( ((*(phenotypes + dadi)) == PH_HH) || ((*(phenotypes + dadi)) == PH_HS) || ((*(phenotypes + dadi)) == PH_HU) )
                    maleMatingProb = (BASE_MATING_FITNESS + (0.5 * MELANISTIC_MATING_ADVANTAGE_M)) / (BASE_MATING_FITNESS + MELANISTIC_MATING_ADVANTAGE_M);
                else
                    maleMatingProb = BASE_MATING_FITNESS / (BASE_MATING_FITNESS + MELANISTIC_MATING_ADVANTAGE_M);
            }
            else {
                if ( *(phenotypes + dadi) == PH_MELANISTIC )
                    maleMatingProb = 1.0;
                else
                    maleMatingProb = BASE_MATING_FITNESS / (BASE_MATING_FITNESS + MELANISTIC_MATING_ADVANTAGE_M);
            }
            
        } while ( (randU() > maleMatingProb) && (triesHere < MAXTRIES) );
        
        seekingMom = 1;
        choiceTrials = 0;
        nicheTries = 0;
        do {
            // get a female who is alive
            triesHere = 0;
            do {
                momTries++;
                triesHere++;
                dumi =randI() % femaleCount;
                momi = *(femaleIndexes + dumi);
            } while ( ( *(markedDead + momi) == 1 ) && (triesHere < MAXTRIES) );
            
            // check niche match
            if ( *(nicheWithinPatch + dadi) != *(nicheWithinPatch + momi) ) {
                if ( randU() < BIASED_NICHE_MATING_RATIO )
                    proceed = 1;
                else
                    proceed = 0;
            }
            else
                proceed = 1;
            
            if ( proceed ) {
                // now do a choice trial
                
                paccept = BASE_MATING_FITNESS;
                
                if ( *(chcType + dadi) == *(chcType + momi) )
                    paccept += CHC_MATCH_ADVANTAGE;
                
                if ( ADDITIVE_MODEL ) {
                    if ( ((*(phenotypes + momi)) == PH_MH) || ((*(phenotypes + momi)) == PH_MU) || ((*(phenotypes + momi)) == PH_MS) )
                        paccept += MELANISTIC_MATING_ADVANTAGE_F;
                    else if ( ((*(phenotypes + momi)) == PH_HH) || ((*(phenotypes + momi)) == PH_HU) || ((*(phenotypes + momi)) == PH_HS) )
                        paccept += 0.5 * MELANISTIC_MATING_ADVANTAGE_F;
                }
                else {
                    if ( *(phenotypes + momi) == PH_MELANISTIC )
                        paccept += MELANISTIC_MATING_ADVANTAGE_F;
                }
                    
                paccept = paccept / (BASE_MATING_FITNESS + CHC_MATCH_ADVANTAGE + MELANISTIC_MATING_ADVANTAGE_F);
                
                if ( randU() <= paccept )
                    seekingMom = 0;
                else
                    choiceTrials++;
            }
            else
                nicheTries++;
            
        } while ( seekingMom && (choiceTrials <= MAX_NUM_CHOICE_TRIALS) && ( nicheTries < MAXTRIES ) );
        
    } while ( seekingMom && functionTries < MAXTRIES );
    
    if ( seekingMom ) {
        fprintf(stderr, "\nWarning in chooseParents2():\n\tfailed to find two parents!\n\tLast two randomly chosen parents being returned.\n");
    }
    
    
    *momPt = momi;
    *dadPt = dadi;
    
    //fprintf(stderr, "%i\t%i\t%i\t%i\t%i\t%f\n", dadi,momi,nicheTries,choiceTrials,functionTries, paccept);
}


//double computeFemaleFitness(int dadi, int momi)
//{
//    double f = BASE_MATING_FITNESS;
//    
//    if ( chcType[dadi] == chcType[momi] )
//        f += CHC_MATCH_ADVANTAGE;
//    
//    if ( ADDITIVE_MODEL ) {
//        if ( (phenotypes[momi] == PH_MH) || (phenotypes[momi] == PH_MS) || (phenotypes[momi] == PH_MU) )
//            f += MELANISTIC_MATING_ADVANTAGE_F;
//        else if ( (phenotypes[momi] == PH_HH) || (phenotypes[momi] == PH_HS) || (phenotypes[momi] == PH_HU) )
//            f += 0.5 * MELANISTIC_MATING_ADVANTAGE_F;
//    }
//    else {
//        if ( phenotypes[momi] == PH_MELANISTIC )
//            f += MELANISTIC_MATING_ADVANTAGE_F;
//    }
//    
//    f = f * MAX_FITNESS_INV; // normalization
//    
//    return f;
//}



int computePhenotype(int *ipt)
{
    int colorSum, patternSum, *colorpt, *patternpt;
    int homoS, homoU, hetUS, homoG, homoM, hetGM;
    
    
    homoS = STRIPED_ALLELE + STRIPED_ALLELE;
    homoU = UNSTRIPED_ALLELE + UNSTRIPED_ALLELE;
    hetUS = STRIPED_ALLELE + UNSTRIPED_ALLELE;
    homoG = GREEN_ALLELE + GREEN_ALLELE;
    homoM = MELANISTIC_ALLELE + MELANISTIC_ALLELE;
    hetGM = GREEN_ALLELE + MELANISTIC_ALLELE;
    
    // ipt is the pointer to first allele of this individual at its first locus
    colorpt = ipt + (DIPLOID * INDEX_OF_COLOR_LOCUS);
    patternpt = ipt + (DIPLOID * INDEX_OF_PATTERN_LOCUS);
    
    colorSum = (*colorpt) + (*(colorpt+1));
    patternSum = (*patternpt) + (*(patternpt+1));
    
    if ( ADDITIVE_MODEL ) {
        if ( colorSum == homoG ) {
            if ( patternSum == homoS )
                return PH_GS;
            else if ( patternSum == hetUS )
                return PH_GH;
            else if ( patternSum == homoU )
                return PH_GU;
            else {
                fprintf(stderr, "\nError in computePhenotype():\n\tpatternSum = %i not found\n", patternSum);
                exit(-1);
            }
        }
        else if ( colorSum == hetGM ) {
            if ( patternSum == homoS )
                return PH_HS;
            else if ( patternSum == hetUS )
                return PH_HH;
            else if ( patternSum == homoU )
                return PH_HU;
            else {
                fprintf(stderr, "\nError in computePhenotype():\n\tpatternSum = %i not found\n", patternSum);
                exit(-1);
            }
            
        }
        else if ( colorSum == homoM ) {
            if ( patternSum == homoS )
                return PH_MS;
            else if ( patternSum == hetUS )
                return PH_MH;
            else if ( patternSum == homoU )
                return PH_MU;
            else {
                fprintf(stderr, "\nError in computePhenotype():\n\tpatternSum = %i not found\n", patternSum);
                exit(-1);
            }
        }
        else {
            fprintf(stderr, "\nError in computePhenotype(): colorSum = %i not found!\n\tExiting...\n\n", colorSum);
            exit(-1);
        }
    }
    else {
        if ( colorSum < 1 ) {
            return PH_MELANISTIC; // recessive color expressed, causing pattern alleles to not matter
        }
        else {
            if ( patternSum == homoS )
                return PH_STRIPED_GREEN_HOMO; // recessive
            else if ( patternSum == hetUS )
                return PH_HETERO_STRIPE; // dominant
            else if ( patternSum == homoU )
                return PH_UNSTRIPED_GREEN_HOMO;
            else {
                fprintf(stderr, "\nError in computePhenotype():\n\tpatternSum = %i not found\n", patternSum);
                exit(-1);
            }
        }
    }
    
}


int computeSex(int *ipt)
{
    // ipt is the pointer to first allele of this individual at its first locus
    int a1, a2;
    
    a1 = *(ipt + (DIPLOID * INDEX_OF_SEX_LOCUS));
    a2 = *(ipt + 1 + (DIPLOID * INDEX_OF_SEX_LOCUS));
    
    if ( a1 == MALE_ALLELE && a2 == MALE_ALLELE ) {
        fprintf(stderr, "\nError in genotypes found in computeSex():\n\ttwo 'Y' chromosomes!\n");
        exit(-1);
    }
    else if ( a1 == MALE_ALLELE || a2 == MALE_ALLELE )
        return SEX_MALE;
    else if ( a1 == FEMALE_ALLELE && a2 == FEMALE_ALLELE )
        return SEX_FEMALE;
    else {
        fprintf(stderr, "\nError in computeSex():\n\tSex not found! a1 = %i, a2 = %i\n", a1, a2);
        exit(-1);
    }
}


void finalPrintsMultirun(FILE *mrfpt, int mrcount, int rngseed, double IHmatFreq, double mr, double *FST)
{
    FILE *fpt;
    int count, i, j;
    double propad, ldval;
    
    // the list of stuff and the order printed in the header on teh file
//
//
//
//    for ( i = 0; i < (nLOCI-1); i++ ) {
//        for ( j = (i+1); j < nLOCI; j++ ) {
//            fprintf(fpt2, ",LDloci%iand%i", i, j);
//        }
//    }
    
    propad = 0.0;
    for ( i = 0; i < nPATCHES; i++ ) {
        if ( hostType[i] == ADENOSTOMA_HOST ) {
            propad += 1.0;
        }
    }
    propad = propad / ((double) nPATCHES);
    
    fprintf(mrfpt,"%i,%i,'%s',%i,%li,%i,%i,%f,", mrcount,rngseed,version,nLOCI,MAX_GENERATIONS,ROWS_OF_PATCHES,COLS_OF_PATCHES,PROPORTION_HOST_DARK_NICHE);
    fprintf(mrfpt,"%f,%f,%i,%f,", PROPORTION_ADENOSTOMA,propad,TOTAL_N,INIT_MEL_FREQ_ADENOSTOMA);
    fprintf(mrfpt,"%f,%li,%f,%f,%i,", INIT_MEL_FREQ_CEANOTHUS,nGENS_ALLOPATRY,SD_MOVE_MELANISTIC,SD_MOVE_GREEN,RANDOM_LOCS_FOR_OFFSPRING);
    fprintf(mrfpt,"%f,%f,%f,", COLOR_PATTERN_RECOMB_PROBABILITY,S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE,S_COEFF_GREEN_IN_DARK_NICHE);
    fprintf(mrfpt,"%f,%f,%f,%f,", S_COEFF_MELANISTIC_IN_GREEN_NICHE,STRIPE_DOMINANCE_COEFFICIENT,BASE_MATING_FITNESS,CHC_MATCH_ADVANTAGE);
    fprintf(mrfpt,"%f,%f,%i,%i,%i,%f,%i,%i,%i", MELANISTIC_MATING_ADVANTAGE_M, MELANISTIC_MATING_ADVANTAGE_F,INDEX_OF_COLOR_LOCUS,INDEX_OF_PATTERN_LOCUS,INDEX_OF_SEX_LOCUS,BIASED_NICHE_MATING_RATIO,LAYOUT_OF_HOSTS_RANDOM,HOLD_MEL_ALLELE_FREQ_CONSTANT,ADDITIVE_MODEL);
    
    for ( i = 1; i < nLOCI; i++ )
        fprintf(mrfpt, ",%E", recombinationProbabilities[i]);
    fprintf(mrfpt, ",%li,%E,%E", (t-1),IHmatFreq,mr);
    
    for ( i = 0; i < nLOCI; i++ )
        fprintf(mrfpt, ",%E", allele_frequencies[i]);
    
    for ( i = 0; i < nLOCI; i++ ) {
        if ( locusStillVariable[i] )
            fprintf(mrfpt, ",%E", FST[i]);
        else
            fprintf(mrfpt, ",NA");
    }
    
    for ( i = 0; i < nLOCI; i++ ) {
        if ( locusStillVariable[i] )
            fprintf(mrfpt, ",NA");
        else
            fprintf(mrfpt, ",%li", fixationTimes[i]);
    }
    
    for ( i = 0; i < (nLOCI-1); i++ ) {
        for ( j = (i+1); j < nLOCI; j++ ) {
            if ( locusStillVariable[i] && locusStillVariable[j] ) {
                ldval = calculateLD( allele_frequencies[i], allele_frequencies[j], i, j );
                fprintf(mrfpt, ",%E", ldval);
            }
            else
                fprintf(mrfpt, ",NA");
        }
    }
    
    fprintf(mrfpt, "\n");
    
    
    
    fpt = fopen("MultirunCount.txt","r");
    fscanf(fpt, "%i", &count);
    fprintf(stderr, "\nRun %i finished\n\t\t**************\n\n", count);
    fclose(fpt);
}


void finalPrintsAllRuns(int rngseed)
{
    int i, j, *ipt, patchCount;
    FILE *params;
    double *dpt;
    
    //final population
    if ( !MULTIRUN )
        printPopulation( "FinalPopulation.csv" );
    
    // files for easily sourcing/reading in parameters and useful run information:
    if ( MULTIRUN )
        params = fopen("LastRunParameters.R","w");
    else
        params = fopen("parameters.R","w");
    
    fprintf(params, "codeVersion <- '%s';\n", version);
    fprintf(params, "rngseed <- %i; #random number seed for this run\n\n", rngseed);
    
    fprintf(params, "# parameter values:\n");
    fprintf(params, "nLOCI <- %i;\n", nLOCI);
    fprintf(params, "MAX_GENERATIONS <- %li;\n", MAX_GENERATIONS);
    fprintf(params, "ROWS_OF_PATCHES <- %i;\n", ROWS_OF_PATCHES);
    fprintf(params, "COLS_OF_PATCHES <- %i;\n", COLS_OF_PATCHES);
    fprintf(params, "PROPORTION_HOST_DARK_NICHE <- %f;\n", PROPORTION_HOST_DARK_NICHE);
    fprintf(params, "PROPORTION_ADENOSTOMA <- %f;\n", PROPORTION_ADENOSTOMA);
    fprintf(params, "TOTAL_N <- %i;\n", TOTAL_N);
    fprintf(params, "INIT_MEL_FREQ_ADENOSTOMA <- %f;\n", INIT_MEL_FREQ_ADENOSTOMA);
    fprintf(params, "INIT_MEL_FREQ_CEANOTHUS <- %f;\n", INIT_MEL_FREQ_CEANOTHUS);
    fprintf(params, "nGENS_ALLOPATRY <- %li;\n", nGENS_ALLOPATRY);
    fprintf(params, "SD_MOVE_MELANISTIC <- %f;\n", SD_MOVE_MELANISTIC);
    fprintf(params, "SD_MOVE_GREEN <- %f;\n", SD_MOVE_GREEN);
    fprintf(params, "RANDOM_LOCS_FOR_OFFSPRING <- %i;\n", RANDOM_LOCS_FOR_OFFSPRING);
    fprintf(params, "COLOR_PATTERN_RECOMB_PROBABILITY <- %f;\n", COLOR_PATTERN_RECOMB_PROBABILITY);
    fprintf(params, "TIME_SAMPLING_INTERVAL <- %i;\n", TIME_SAMPLING_INTERVAL);
    fprintf(params, "S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE <- %f;\n", S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE);
    fprintf(params, "S_COEFF_GREEN_IN_DARK_NICHE <- %f;\n", S_COEFF_GREEN_IN_DARK_NICHE);
    fprintf(params, "S_COEFF_MELANISTIC_IN_GREEN_NICHE <- %f;\n", S_COEFF_MELANISTIC_IN_GREEN_NICHE);
    fprintf(params, "STRIPE_DOMINANCE_COEFFICIENT <- %f;\n", STRIPE_DOMINANCE_COEFFICIENT);
    fprintf(params, "BASE_MATING_FITNESS <- %f;\n", BASE_MATING_FITNESS);
    fprintf(params, "CHC_MATCH_ADVANTAGE <- %f;\n", CHC_MATCH_ADVANTAGE);
    fprintf(params, "MELANISTIC_MATING_ADVANTAGE_M <- %f;\n", MELANISTIC_MATING_ADVANTAGE_M);
    fprintf(params, "MELANISTIC_MATING_ADVANTAGE_F <- %f;\n", MELANISTIC_MATING_ADVANTAGE_F);
    fprintf(params, "BIASED_NICHE_MATING_RATIO <- %f;\n", BIASED_NICHE_MATING_RATIO);
    fprintf(params, "LAYOUT_OF_HOSTS_RANDOM <- %i;\n", LAYOUT_OF_HOSTS_RANDOM);
    fprintf(params, "HOLD_MEL_ALLELE_FREQ_CONSTANT <- %i;\n", HOLD_MEL_ALLELE_FREQ_CONSTANT);
    fprintf(params, "ADDITIVE_MODEL <- %i;\n", ADDITIVE_MODEL);
    
    fprintf(params, "\n# allele, phenotype, and host coding flags; only make sense in appropriate context:\n");
    fprintf(params, "# color locus and its alleles:\n");
    fprintf(params, "INDEX_OF_COLOR_LOCUS <- %i;\n", INDEX_OF_COLOR_LOCUS);
    fprintf(params, "GREEN_ALLELE <- %i;\n", GREEN_ALLELE);
    fprintf(params, "MELANISTIC_ALLELE <- %i;\n", MELANISTIC_ALLELE);
    
    fprintf(params, "# pattern locus and its alleles:\n");
    fprintf(params, "INDEX_OF_PATTERN_LOCUS <- %i;\n", INDEX_OF_PATTERN_LOCUS);
    fprintf(params, "UNSTRIPED_ALLELE <- %i;\n", UNSTRIPED_ALLELE);
    fprintf(params, "STRIPED_ALLELE <- %i;\n", STRIPED_ALLELE);
    
    fprintf(params, "# sex locus and alleles; these 'alleles' work like XX X0 chromosomes:\n");
    fprintf(params, "INDEX_OF_SEX_LOCUS <- %i;\n", INDEX_OF_SEX_LOCUS);
    fprintf(params, "MALE_ALLELE <- %i;\n", MALE_ALLELE);
    fprintf(params, "FEMALE_ALLELE <- %i;\n", FEMALE_ALLELE);
    
    fprintf(params, "# neutral loci and alleles:\n");
    fprintf(params, "INDEXES_OF_NEUTRAL_LOCI <- c(");
    j = 0; // dummy counter variable for getting syntax correct for R in following for loop
    for ( i = 0; i < nLOCI; i++ ) {
        if ( locusTypes[i] == NEUTRAL_LOCUS ) {
            if ( j > 0 )
                fprintf(params, ",%i",i);
            else
                fprintf(params, "%i",i);
            j++;
        }
    }
    fprintf(params, ");\n");
    fprintf(params, "NEUT_ADENO_ALLELE <- %i;\n", NEUT_ADENO_ALLELE);
    fprintf(params, "NEUT_CEANO_ALLELE <- %i;\n", NEUT_CEANO_ALLELE);

    
    fprintf(params, "# color-pattern phenotype flags:\n");
    if ( !ADDITIVE_MODEL ) {
        fprintf(params, "PH_MELANISTIC <- %i;\n", PH_MELANISTIC);
        fprintf(params, "PH_STRIPED_GREEN_HOMO <- %i;\n", PH_STRIPED_GREEN_HOMO);
        fprintf(params, "PH_UNSTRIPED_GREEN_HOMO <- %i;\n", PH_UNSTRIPED_GREEN_HOMO);
        fprintf(params, "PH_HETERO_STRIPE <- %i;\n", PH_HETERO_STRIPE);
    }
    else {
        fprintf(params, "PH_GU <- %i; #double homozygous green unstriped\n", PH_GU);
        fprintf(params, "PH_GH <- %i; #homozygous green, heterozygous stripe\n", PH_GH);
        fprintf(params, "PH_GS <- %i; #double homozygous green striped\n", PH_GS);
        fprintf(params, "PH_HU <- %i; #heterozygous color, homozygous unstriped\n", PH_HU);
        fprintf(params, "PH_HH <- %i; #double heterozygote\n", PH_HH);
        fprintf(params, "PH_HS <- %i; #heterozygous color, homozygous striped\n", PH_HS);
        fprintf(params, "PH_MU <- %i; #double homozygous melanistic unstriped\n", PH_MU);
        fprintf(params, "PH_MH <- %i; #homozygous melanistic, heterozygous stripe\n", PH_MH);
        fprintf(params, "PH_MS <- %i; #double homozygous melanistic striped\n", PH_MS);
    }
    
    
    fprintf(params, "\n# Host plant spatial arrangement information:\n");
    fprintf(params, "# Host type flags:\n");
    fprintf(params, "ADENOSTOMA_HOST <- %i;\n", ADENOSTOMA_HOST);
    fprintf(params, "CEANOTHUS_HOST <- %i;\n", CEANOTHUS_HOST);
    fprintf(params, "NICHE_GREEN <- %i;\n", NICHE_GREEN);
    fprintf(params, "NICHE_DARK <- %i;\n", NICHE_DARK);
    fprintf(params, "# Host plant locations in this simulation run:\n");
    fprintf(params,"hostType <- matrix( c(");
    patchCount = 0;
    for ( i = 0; i < ROWS_OF_PATCHES; i++ ) {
        for ( j = 0; j < COLS_OF_PATCHES; j++ ) {
            fprintf(params, "%i", hostType[patchCount]);
            if ( patchCount < (nPATCHES-1) )
                fprintf(params, ",");
            patchCount++;
        }
    }
    fprintf(params,"), nrow=%i, ncol=%i, byrow=TRUE );\n", ROWS_OF_PATCHES, COLS_OF_PATCHES);
    fprintf(params,"hostPatchIndexes <- matrix( c(");
    patchCount = 0;
    for ( i = 0; i < ROWS_OF_PATCHES; i++ ) {
        for ( j = 0; j < COLS_OF_PATCHES; j++ ) {
            fprintf(params, "%i", patchCount);
            if ( patchCount < (nPATCHES-1) )
                fprintf(params, ",");
            patchCount++;
        }
    }
    fprintf(params,"), nrow=%i, ncol=%i, byrow=TRUE );\n", ROWS_OF_PATCHES, COLS_OF_PATCHES);
    fprintf(params,"# NOTE, for visualizing hostType I suggest the following commands, which I modified from\n# http://blog.snap.uaf.edu/2012/06/08/matrix-rotation-for-image-and-contour-plots-in-r/\n");
    fprintf(params,"myImage <- function(m) {newm <- t(m)[,nrow(m):1]; image(newm); grid(dim(newm)[1],dim(newm)[2], col = %cblack%c)};\n", 34, 34);
    fprintf(params,"# myImage(hostType)\n");
    
    fprintf(params, "\n# Potentially useful values of certain variables at end of simulation:\n");
    fprintf(params, "finalTimeAtEnd <- %li;\n", (t-1)); // it's t-1 because t is actually incremented at the END of the main work loop
    fprintf(params, "fixationTimesAtEachLocus <- c(");
    for ( i = 0; i < nLOCI; i++ ) {
        if ( i > 0 )
            fprintf(params, ",");
        if ( locusStillVariable[i] )
            fprintf(params, "NA");
        else
            fprintf(params, "%li", fixationTimes[i]);
    }
    fprintf(params, ");\n");
    
    fclose(params);
    
}


int findNiche(double x, double y)
{
    double offsetx, offsety, foo;
    
    // code for dark border and green interior
//    offsetx = modf( x , &foo );
//    offsety = modf( y , &foo );
//    if ( offsetx < darkNicheLowerBound || offsetx > darkNichUpperBound || offsety < darkNicheLowerBound || offsety > darkNichUpperBound ) {
//        return NICHE_DARK;
//    }
//    else
//        return NICHE_GREEN;
    
    // code for division simply into subrectangles
    offsetx = modf( x , &foo );
    if ( offsetx < PROPORTION_HOST_DARK_NICHE )
        return NICHE_DARK;
    else
        return NICHE_GREEN;
    
}


void getOffspringAlleles(int parentIndex, int *offGtPt)
{
    int i, j, *pgtpt;
    int haplotype;
    
    // parent genotype
    haplotype = 0;
    pgtpt = genotypes + ( DIPLOID * nLOCI * parentIndex );
    if ( randU() < 0.5 ) {
        // switch to other haplotype to start
        pgtpt++;
        haplotype = 1;
    }
    
    // first locus
    *offGtPt = *pgtpt;
    offGtPt += DIPLOID;
    
    // remaining loci
    for ( i = 1; i < nLOCI; i++ ) {
        if ( randU() < recombinationProbabilities[i] ) {
            // switch parental haplotype
            if ( haplotype == 1 ) {
                pgtpt++;
                haplotype = 0;
            }
            else {
                pgtpt += 3;
                haplotype = 1;
            }
        }
        else
            pgtpt += DIPLOID;
        
        *offGtPt = *pgtpt;
        
        offGtPt += DIPLOID;
    }
    
}


void initializeGenomeScheme(void)
{
    // this is a specific, custom scheme designed to answer a number of questions
    // about divergence as measured several ways
    nLOCI = 6;
    
    locusTypes = (int *) malloc( (sizeof(int) * nLOCI) );
    recombinationProbabilities = (double *) malloc( (sizeof(double) * nLOCI) );

    recombinationProbabilities[0] = 0.0; // this one is just a space filler
    
    
    locusTypes[0] = NEUTRAL_LOCUS;

    recombinationProbabilities[1] = 0.5; // b/t 0 and 1
    
    locusTypes[1] = NEUTRAL_LOCUS;
    
    recombinationProbabilities[2] = 0.001; // b/t 1 and 2
    
    locusTypes[2] = COLOR_LOCUS;
    INDEX_OF_COLOR_LOCUS = 2;
    
    recombinationProbabilities[3] = COLOR_PATTERN_RECOMB_PROBABILITY; // b/t 2 and 3, the color and pattern loci
    
    locusTypes[3] = PATTERN_LOCUS;
    INDEX_OF_PATTERN_LOCUS = 3;
    
    recombinationProbabilities[4] = 0.001;
    
    locusTypes[4] = NEUTRAL_LOCUS;

    recombinationProbabilities[5] = 0.5;
    
    locusTypes[5] = SEX_LOCUS; // like a sex chromosome
    INDEX_OF_SEX_LOCUS = 5;
    
    // record keeping
    FILE *gsf;
    int i;
    
    if ( ! MULTIRUN ) {
        gsf = fopen("GenomeScheme.csv","w");
        if ( PRINT_HEADERS )
            fprintf(gsf, "LocusIDnumber,LocusTypeCode,LocusType,RecombProbabilityWithPrevLocus\n");
        for ( i = 0; i < nLOCI; i++ ) {
            fprintf(gsf, "%i,%i,", i, locusTypes[i]);
            
            if ( locusTypes[i] == NEUTRAL_LOCUS )
                fprintf(gsf, "NEUTRAL_LOCUS,");
            else if ( locusTypes[i] == COLOR_LOCUS )
                fprintf(gsf, "COLOR_LOCUS,");
            else if ( locusTypes[i] == PATTERN_LOCUS )
                fprintf(gsf, "PATTERN_LOCUS,");
            else if ( locusTypes[i] == SEX_LOCUS )
                fprintf(gsf, "SEX_LOCUS,");
            //        else if ( locusTypes[i] == ASSORTMENT_LOCUS )
            //            fprintf(gsf, "ASSORTMENT_LOCUS,");
            else {
                fprintf(stderr, "\nError in initializeGenomeScheme():\n\ttype for locusTypes[i] (= %i) not recognized. Exiting\n", locusTypes[i]);
                exit(-1);
            }
            
            if ( i == 0 )
                fprintf(gsf, "NA\n");
            else
                fprintf(gsf, "%f\n", recombinationProbabilities[i]);
            
        }
        
        fclose(gsf);
    }
    
}


void initializePatches(void)
{
    int i, j, patchCount, aTarget, nAdeno, nCeano, halfwayIndex;
    double darkNicheWidth, foo;
    //FILE *hostLayout;
    
    nPATCHES = ROWS_OF_PATCHES * COLS_OF_PATCHES;
    
    
    aTarget = (int) ((((double) nPATCHES) * PROPORTION_ADENOSTOMA) + 0.5);
    
    hostType = (int *) malloc( (sizeof(int) * nPATCHES) );
    
    hostType[0] = ADENOSTOMA_HOST;
    hostType[(nPATCHES - 1)] = CEANOTHUS_HOST;
    nAdeno = 1;
    nCeano = 1;
    if ( nPATCHES > 2 ) {
        if ( LAYOUT_OF_HOSTS_RANDOM ) {
            for ( i = 1; i <= (nPATCHES - 2); i++ ) {
                if ( randU() < PROPORTION_ADENOSTOMA ) {
                    hostType[i] = ADENOSTOMA_HOST;
                    nAdeno++;
                }
                else {
                    hostType[i] = CEANOTHUS_HOST;
                    nCeano++;
                }
            }
        }
        else {
            halfwayIndex = aTarget - 1;
            
            for ( i = 1; i < (nPATCHES-1); i++ ) {
                if ( i <= halfwayIndex ) {
                    hostType[i] = ADENOSTOMA_HOST;
                    nAdeno++;
                }
                else {
                    hostType[i] = CEANOTHUS_HOST;
                    nCeano++;
                }
            }
            // check
            //fprintf(stderr, "\nhalfwayIndex = %i, nPATCHES = %i, nAdeno = %i, nCeano = %i\n", halfwayIndex, nPATCHES, nAdeno, nCeano);
            
            if ( ((nAdeno - aTarget) > 1) || ((aTarget - nAdeno) > 1) ) {
                fprintf(stderr, "\nError!  Patch number imbalance!\nnCeano = %i, nAdeno = %i\n", nCeano, nAdeno);
                exit(-1);
            }
        }
    }
    if ( (nPATCHES > 3) && (aTarget > 1) && (aTarget < (nPATCHES-1)) && LAYOUT_OF_HOSTS_RANDOM ) {
        // check test
        // fprintf(stderr, "\ninitially %i adeno, %i ceano, target = %i adeno. Host types:\n", nAdeno, nCeano, aTarget);
//        for ( i = 0; i < nPATCHES; i++ )
//            fprintf(stderr, "%i ", hostType[i]);
        if ( nAdeno < aTarget ) {
            do {
                // choose random patch index
                i = randI() % (nPATCHES - 2); // random integer between 0 and nPATCHES - 3
                i++; // random integer between 1 and nPATCHES - 2; 0 and nPATCHES-1 are off limits
                if ( i <= 0 || i >= (nPATCHES - 1) ) {
                    fprintf(stderr, "\nError! Patch type assignment scheme not working\n");
                    exit(-1);
                }
                // need to increase nAdeno
                if ( hostType[i] == CEANOTHUS_HOST ) {
                    hostType[i] = ADENOSTOMA_HOST;
                    nAdeno++;
                    nCeano--;
                    // check test
                    // fprintf(stderr, "\n\tpatch %i flipped C to A\n", i);
                }
            } while ( nAdeno < aTarget );
        }
        else if ( nAdeno > aTarget ) {
            do {
                // choose random patch index
                i = randI() % (nPATCHES - 2); // random integer between 0 and nPATCHES - 3
                i++; // random integer between 1 and nPATCHES - 2; 0 and nPATCHES-1 are off limits
                if ( i <= 0 || i >= (nPATCHES - 1) ) {
                    fprintf(stderr, "\nError! Patch type assignment scheme not working\n");
                    exit(-1);
                }
                // need to decrease nAdeno
                if ( hostType[i] == ADENOSTOMA_HOST ) {
                    hostType[i] = CEANOTHUS_HOST;
                    nAdeno--;
                    nCeano++;
                    // check test
                    // fprintf(stderr, "\n\tpatch %i flipped A to C\n", i);
                }
            } while ( nAdeno > aTarget );
        }
    }
    // check test
    nAdeno = 0;
    nCeano = 0;
    for ( i = 0; i < nPATCHES; i++ ) {
        if ( hostType[i] == ADENOSTOMA_HOST )
            nAdeno++;
        else if ( hostType[i] == CEANOTHUS_HOST )
            nCeano++;
    }
    foo = ((double) (nAdeno - aTarget));
    if ( (fabs( foo ) > 1.001) || ((nAdeno + nCeano) != nPATCHES) ) {
        fprintf(stderr, "\nWarning after checking patch numbers in initializePatches():\n\t%i adeno, %i ceano, %i total\n\n", nAdeno, nCeano, (nAdeno + nCeano));
    }
    
    // code for a setup with a "dark" border and a "green" interior
    darkNicheWidth = (1.0 - sqrt( 1.0 - PROPORTION_HOST_DARK_NICHE )) / 2.0;
    if ( darkNicheWidth < 0.0 || darkNicheWidth > 1.0 ) {
        fprintf(stderr, "\nError in initializePatches():\n\t\tdarkNicheWidth (= %f) computed badly\nExiting...\n", darkNicheWidth);
        exit(-1);
    }
    else {
        darkNicheLowerBound = darkNicheWidth;
        darkNichUpperBound = 1.0 - darkNicheWidth;
//        fprintf(stderr, "\nPROPORTION_HOST_DARK_NICHE = %f\n", PROPORTION_HOST_DARK_NICHE);
//        fprintf(stderr, "darkNicheWidth = %f\n", darkNicheWidth);
//        fprintf(stderr, "prop. green niche = %f\n", pow((1.0 - 2.0*darkNicheWidth),2));
//        fprintf(stderr, "prop. dark niche = %f\n",  (1.0 - pow((1.0 - 2.0*darkNicheWidth),2)));
//        fprintf(stderr, "upper and lower bounds:  %f, %f\n", darkNicheLowerBound, darkNichUpperBound);
    }
    
    
    
}


void initializePopulation(void)
{
    int n_per_patch, count, i, j, colorSum, patternSum, lt;
    int currentPatch, patchType, *ipt, *startpt, currentRow, currentColumn;
    double dum, minx, miny, *dpt;
    FILE *initPopn;
    
    n_per_patch = TOTAL_N / nPATCHES;
    
    // individual attributes
    genotypes = (int *) malloc( (sizeof(int) * TOTAL_N * DIPLOID * nLOCI) ); // diploid genotypes
    phenotypes = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    patchLocations = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    coordinates = (double *) malloc( (sizeof(double) * 2 * TOTAL_N) ); // 2-D euclidean coordinates
    nicheWithinPatch = (int *) malloc( (sizeof(int) * TOTAL_N) ); // niche currently occupied by individual
    markedDead = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    sexOfIndividual = (int *) malloc( (sizeof(int) * TOTAL_N) );
    chcType = (int *) malloc( (sizeof(int) * TOTAL_N) );
    
    // other data structures
    nInEachPatch = (int *) malloc( (sizeof(int) * nPATCHES) ); // count of individuals in each patch
    locusStillVariable = (int *) malloc( (sizeof(int) * nLOCI) ); // scalar designation
    fixedAllele = (int *) malloc( (sizeof(int) * nLOCI) );
    fixationTimes = (long int *) malloc( (sizeof(long int) * nLOCI) );
    allele_frequencies = (double *) malloc( (sizeof(double) * nLOCI) );
    
    ipt = genotypes;
    dpt = coordinates;
    
    for ( i = 0; i < nPATCHES; i++ )
        nInEachPatch[i] = 0;
    
    count = 0;
    currentPatch = 0;
    for ( i = 0; i < TOTAL_N; i++ ) {
        
        startpt = ipt;
        
        // set marked dead -- 1
        markedDead[i] = 0;
        
        if ( count >= n_per_patch ) {
            currentPatch++;
            count = 0;
        }
        if ( currentPatch == nPATCHES )
            currentPatch--; // a safety in case there's not a perfectly even number per patch
        patchType = hostType[currentPatch];
        
        // set chc type -- 2
        chcType[i] = patchType;
        
        // set up genotypes -- 3
        for ( j = 0; j < nLOCI; j++ ) {
            lt = locusTypes[j];
            
            // set sex alleles and sex -- 4
            if ( lt == SEX_LOCUS ) { // easiest case; host does not matter
                dum = randU();
                if ( dum < 0.5 ) {
                    // female
                    *ipt = FEMALE_ALLELE;
                    *(ipt + 1) = FEMALE_ALLELE;
                    sexOfIndividual[i] = SEX_FEMALE;
                }
                else {
                    // male
                    sexOfIndividual[i] = SEX_MALE;
                    if ( dum < 0.75 ) {  // make first copy male
                        *ipt = MALE_ALLELE;
                        *(ipt + 1) = FEMALE_ALLELE;
                    }
                    else { // make second copy male
                        *ipt = FEMALE_ALLELE;
                        *(ipt + 1) = MALE_ALLELE;
                    }
                }
            }
            else if ( patchType == ADENOSTOMA_HOST ) {
                if ( lt == NEUTRAL_LOCUS ) {
                    *ipt = NEUT_ADENO_ALLELE;
                    *(ipt + 1) = NEUT_ADENO_ALLELE;
                }
                else if ( lt == COLOR_LOCUS ) {
                    dum = randU();
                    if ( dum < (INIT_MEL_FREQ_ADENOSTOMA * INIT_MEL_FREQ_ADENOSTOMA) ) {
                        // homozygous
                        *ipt = MELANISTIC_ALLELE;
                        *(ipt + 1) = MELANISTIC_ALLELE;
                    }
                    else if ( dum < ((INIT_MEL_FREQ_ADENOSTOMA * INIT_MEL_FREQ_ADENOSTOMA) + (2.0 * INIT_MEL_FREQ_ADENOSTOMA * (1.0 - INIT_MEL_FREQ_ADENOSTOMA))) ) {
                        // heterozygous
                        if ( randU() < 0.5 ) {
                            *ipt = MELANISTIC_ALLELE;
                            *(ipt + 1) = GREEN_ALLELE;
                        }
                        else {
                            *ipt = GREEN_ALLELE;
                            *(ipt + 1) = MELANISTIC_ALLELE;
                        }
                    }
                    else {
                        // other homozygote
                        *ipt = GREEN_ALLELE;
                        *(ipt + 1) = GREEN_ALLELE;
                    }
                }
                else if ( lt == PATTERN_LOCUS ) {
                    *ipt = STRIPED_ALLELE;
                    *(ipt + 1) = STRIPED_ALLELE;
                }
                else {
                    fprintf(stderr, "\n Error in initializePopulation()!\n\tLocus type not recognized!\n");
                    exit(-1);
                }
            }
            else if ( patchType == CEANOTHUS_HOST ) {
                if ( lt == NEUTRAL_LOCUS ) {
                    *ipt = NEUT_CEANO_ALLELE;
                    *(ipt + 1) = NEUT_CEANO_ALLELE;
                }
                else if ( lt == COLOR_LOCUS ) {
                    // color locus
                    dum = randU();
                    if ( dum < (INIT_MEL_FREQ_CEANOTHUS * INIT_MEL_FREQ_CEANOTHUS) ) {
                        // homozygous
                        *ipt = MELANISTIC_ALLELE;
                        *(ipt + 1) = MELANISTIC_ALLELE;
                    }
                    else if ( dum < ((INIT_MEL_FREQ_CEANOTHUS * INIT_MEL_FREQ_CEANOTHUS) + (2.0 * INIT_MEL_FREQ_CEANOTHUS * (1.0 - INIT_MEL_FREQ_CEANOTHUS))) ) {
                        // heterozygous
                        if ( randU() < 0.5 ) {
                            *ipt = MELANISTIC_ALLELE;
                            *(ipt + 1) = GREEN_ALLELE;
                        }
                        else {
                            *ipt = GREEN_ALLELE;
                            *(ipt + 1) = MELANISTIC_ALLELE;
                        }
                    }
                    else {
                        // other homozygote
                        *ipt = GREEN_ALLELE;
                        *(ipt + 1) = GREEN_ALLELE;
                    }
                }
                else if ( lt == PATTERN_LOCUS ) {
                    *ipt = UNSTRIPED_ALLELE;
                    *(ipt + 1) = UNSTRIPED_ALLELE;
                }
                else {
                    fprintf(stderr, "\n Error in initializePopulation()!\n\tLocus type not recognized!\n");
                    exit(-1);
                }
            }
            else {
                fprintf(stderr, "\n Error in initializePopulation()!\n\tHost patch type not recognized!\n");
                exit(-1);
            }
            
            ipt += DIPLOID;

        }
        
        // phenotype -- 5
        phenotypes[i] = computePhenotype(startpt);
        
        // patchLocations -- 6
        patchLocations[i] = currentPatch;
        nInEachPatch[currentPatch] = nInEachPatch[currentPatch] + 1;
        
        // coordinates -- 7
        currentRow = currentPatch / COLS_OF_PATCHES;
        currentColumn = currentPatch % COLS_OF_PATCHES;
        miny = (double) currentRow;
        minx = (double) currentColumn;
        *dpt = minx + randU();
        *(dpt + 1) = miny + randU();
        // set niche -- 8
        nicheWithinPatch[i] = findNiche( (*dpt), (*(dpt+1)) );
        
        count++;
        dpt += 2; // 2 dimensional space
    }
    
    if ( !MULTIRUN )
        printPopulation( "InitialPopulation.csv" );
    
    for ( i = 0; i < nLOCI; i++ ) {
        locusStillVariable[i] = 1;
        fixationTimes[i] = -1;
        fixedAllele[i] = -1;
    }
    
    
    //MAX_FITNESS_INV = 1.0 / ( BASE_MATING_FITNESS + CHC_MATCH_ADVANTAGE + MELANISTIC_MATING_ADVANTAGE_M + MELANISTIC_MATING_ADVANTAGE_F );
    
}

void makeViabilityScoeffMatrix(void)
{
    // function specific to the ADDITIVE_MODEL
    
    // recall order of the phenotypes:
//    // phenotype flags in purely additive model
//    #define PH_GU 0 // double homozygous green unstriped;
//    #define PH_GH 1 // homozygous green, heterozygous stripe;
//    #define PH_GS 2 // double homozygous green striped;
//    #define PH_HU 3 // heterozygous color; homozygous unstriped;
//    #define PH_HH 4 // heterozygous color; heterozygous stripe;
//    #define PH_HS 5 // heterozygous color; homozygous striped;
//    #define PH_MU 6 // double homozygous melanistic unstriped;
//    #define PH_MH 7 // homozygous melanistic, heterozygous stripe;
//    #define PH_MS 8 // double homozygous melanistic striped;
    
//    // // host flags
//#define ADENOSTOMA_HOST 0
//#define CEANOTHUS_HOST 1
//    
//    // niche flags
//#define NICHE_GREEN 0
//#define NICHE_DARK 1
    
    // viabilityMatrix is a 3-D array:
    // first dimension, "rows", will index phenotype, from 0 to 8 (see above)
    // second dimension, "columns", will host type, 0 or 1  (see above)
    // third dimension, "depth", will index niche type, 0 or 1 (see above)

    // Homozygous green, unstriped bug GGUU
    viabilityMatrix[0][0][0] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE; // PH_GU, ADENO, NICHE_GREEN
    viabilityMatrix[0][0][1] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE + S_COEFF_GREEN_IN_DARK_NICHE; // PH_GU, ADENO, NICHE_DARK
    viabilityMatrix[0][1][0] = 0.0; // PH_GU, CEANO, NICHE_GREEN; good fit!
    viabilityMatrix[0][1][1] = S_COEFF_GREEN_IN_DARK_NICHE; // PH_GU, CEANO, NICHE_DARK
    
    // homozygous green, heterozygous stripe GGUs
    viabilityMatrix[1][0][0] = 0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE; // PH_GH, ADENO, NICHE_GREEN
    viabilityMatrix[1][0][1] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) + S_COEFF_GREEN_IN_DARK_NICHE; // PH_GH, ADENO, NICHE_DARK
    viabilityMatrix[1][1][0] = 0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE; // PH_GH, CEANO, NICHE_GREEN
    viabilityMatrix[1][1][1] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) + S_COEFF_GREEN_IN_DARK_NICHE; // PH_GH, CEANO, NICHE_DARK
    
    // homozygous green, homozygous striped GGss
    viabilityMatrix[2][0][0] = 0.0; // PH_GS, ADENO, NICHE_GREEN; good fit!
    viabilityMatrix[2][0][1] = S_COEFF_GREEN_IN_DARK_NICHE; // PH_GS, ADENO, NICHE_DARK
    viabilityMatrix[2][1][0] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE; // PH_GS, CEANO, NICHE_GREEN
    viabilityMatrix[2][1][1] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE + S_COEFF_GREEN_IN_DARK_NICHE; // PH_GS, CEANO, NICHE_DARK
    
    // heterozygous color, homozygous UNstriped GmUU
    viabilityMatrix[3][0][0] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE + (0.5 * S_COEFF_MELANISTIC_IN_GREEN_NICHE); // PH_HU, ADENO, NICHE_GREEN
    viabilityMatrix[3][0][1] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE + (0.5 * S_COEFF_GREEN_IN_DARK_NICHE); // PH_HU, ADENO, NICHE_DARK
    viabilityMatrix[3][1][0] = (0.5 * S_COEFF_MELANISTIC_IN_GREEN_NICHE); // PH_HU, CEANO, NICHE_GREEN
    viabilityMatrix[3][1][1] = (0.5 * S_COEFF_GREEN_IN_DARK_NICHE); // PH_HU, CEANO, NICHE_DARK
    
    // heterozygous color, heterozygous stripe GmUs
    viabilityMatrix[4][0][0] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) + (0.5 * S_COEFF_MELANISTIC_IN_GREEN_NICHE); // PH_HH, ADENO, NICHE_GREEN
    viabilityMatrix[4][0][1] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) + (0.5 * S_COEFF_GREEN_IN_DARK_NICHE); // PH_HH, ADENO, NICHE_DARK
    viabilityMatrix[4][1][0] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) + (0.5 * S_COEFF_MELANISTIC_IN_GREEN_NICHE); // PH_HH, CEANO, NICHE_GREEN
    viabilityMatrix[4][1][1] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) + (0.5 * S_COEFF_GREEN_IN_DARK_NICHE); // PH_HH, CEANO, NICHE_DARK
    
    // heterozygous color, homozygous striped Gmss
    viabilityMatrix[5][0][0] = (0.5 * S_COEFF_MELANISTIC_IN_GREEN_NICHE); // PH_HS, ADENO, NICHE_GREEN
    viabilityMatrix[5][0][1] = (0.5 * S_COEFF_GREEN_IN_DARK_NICHE); // PH_HS, ADENO, NICHE_DARK
    viabilityMatrix[5][1][0] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE + (0.5 * S_COEFF_MELANISTIC_IN_GREEN_NICHE); // PH_HS, CEANO, NICHE_GREEN
    viabilityMatrix[5][1][1] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE + (0.5 * S_COEFF_GREEN_IN_DARK_NICHE); // PH_HS, CEANO, NICHE_DARK
    
    // homozygous color, homozygous UNstriped mmUU
    viabilityMatrix[6][0][0] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE + S_COEFF_MELANISTIC_IN_GREEN_NICHE; // PH_MU, ADENO, NICHE_GREEN
    viabilityMatrix[6][0][1] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE; // PH_MU, ADENO, NICHE_DARK
    viabilityMatrix[6][1][0] = S_COEFF_MELANISTIC_IN_GREEN_NICHE; // PH_MU, CEANO, NICHE_GREEN
    viabilityMatrix[6][1][1] = 0.0; // PH_MU, CEANO, NICHE_DARK; good fit!
    
    // homozygous color, heterozygous stripe mmUs
    viabilityMatrix[7][0][0] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) + S_COEFF_MELANISTIC_IN_GREEN_NICHE; // PH_MH, ADENO, NICHE_GREEN
    viabilityMatrix[7][0][1] = 0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE; // PH_MH, ADENO, NICHE_DARK
    viabilityMatrix[7][1][0] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) + S_COEFF_MELANISTIC_IN_GREEN_NICHE; // PH_MH, CEANO, NICHE_GREEN
    viabilityMatrix[7][1][1] = (0.5 * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE); // PH_MH, CEANO, NICHE_DARK
    
    // homozygous color and stripe mmss
    viabilityMatrix[8][0][0] = S_COEFF_MELANISTIC_IN_GREEN_NICHE; // PH_MS, ADENO, NICHE_GREEN
    viabilityMatrix[8][0][1] = 0.0; // PH_MS, ADENO, NICHE_DARK; good fit!
    viabilityMatrix[8][1][0] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE + S_COEFF_MELANISTIC_IN_GREEN_NICHE; // PH_MS, CEANO, NICHE_GREEN
    viabilityMatrix[8][1][1] = S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE; // PH_MS, CEANO, NICHE_DARK
    
    
    int i, j, k;
    FILE *vimat;
    vimat = fopen("AdditiveModelViabilityMatrix.csv", "w");
    if ( PRINT_HEADERS )
        fprintf(vimat, "PhenotypeCode,HostTypeCode,NicheCode,Scoeff\n");
    for ( i = 0; i < 9; i++ ) {
        for ( j = 0; j < 2; j++ ) {
            for ( k = 0; k < 2; k++ ) {
                fprintf(vimat, "%i,%i,%i,%f\n", i, j, k, viabilityMatrix[i][j][k]);
            }
        }
    }
    fclose(vimat);
    
}


double migration(void)
{
    int i, j, oldPatch, newPatch, migrationCount = 0;
    double distance, direction, x, y, *xpt, *ypt;
    int xi, yi;
    double oldx, oldy, foo;
    _Bool once;
    
    
    xpt = coordinates;
    ypt = coordinates + 1;
    
    once = 1;
    
    for ( i = 0; i < TOTAL_N; i++ ) {
        direction = TWOPI * randU(); // direction in radians
        if ( ADDITIVE_MODEL ) {
            if ( phenotypes[i] == PH_MH || phenotypes[i] == PH_MU || phenotypes[i] == PH_MS )
                distance = fabs( boxMuller(0.0, SD_MOVE_MELANISTIC) );
            else if ( phenotypes[i] == PH_HH || phenotypes[i] == PH_HU || phenotypes[i] == PH_HS )
                distance = fabs( boxMuller(0.0, (0.5 * (SD_MOVE_GREEN + SD_MOVE_MELANISTIC))) );
            else
                distance = fabs( boxMuller(0.0, SD_MOVE_GREEN) );
        }
        else {
            if ( phenotypes[i] == PH_MELANISTIC )
                distance = fabs( boxMuller(0.0, SD_MOVE_MELANISTIC) );
            else
                distance = fabs( boxMuller(0.0, SD_MOVE_GREEN) );
        }
        
        x = (*xpt);
        y = (*ypt);
        
        oldx = x;
        oldy = y;
        
        x = x + (cos(direction) * distance);
        if ( x > ((double) COLS_OF_PATCHES) )
            *xpt = ((double) COLS_OF_PATCHES);
        else if ( x < 0.0 )
            *xpt = 0.0;
        else
            *xpt = x;
        
        y = y + (sin(direction) * distance);
        if ( y > ((double) ROWS_OF_PATCHES) )
            *ypt = ((double) ROWS_OF_PATCHES);
        else if ( y < 0.0 )
            *ypt = 0.0;
        else 
            *ypt = y;
        
        
        
        
        xi = (int) (*xpt); // column
        if ( xi == COLS_OF_PATCHES )
            xi--;
        yi = (int) (*ypt); // row
        if ( yi == ROWS_OF_PATCHES )
            yi--;
        
        newPatch = (yi * COLS_OF_PATCHES) + xi;
        oldPatch = patchLocations[i];
        
        if ( newPatch != oldPatch ) {
            migrationCount++;
            patchLocations[i] = newPatch;
            nInEachPatch[oldPatch] = nInEachPatch[oldPatch] - 1;
            nInEachPatch[newPatch] = nInEachPatch[newPatch] + 1;
            
            if ( randU() < MIGRATION_DEATH_PROBABILITY )
                markedDead[i] = 1;
//            if ( t <= 5 ) {
//                if ( once ) {
//                    fprintf(stderr, "\ni\tnewpat\toldpat\txi\tyi\toldx\t\toldy\t\tx\t\ty\t\tdistance\tdirection\n");
//                    once = 0;
//                }
//                fprintf(stderr, "%i\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\n", i, newPatch, oldPatch, xi, yi, oldx, oldy, *xpt, *ypt, distance, direction);
//            }
        }
        
        nicheWithinPatch[i] = findNiche( (*xpt), (*ypt) );
        
        xpt += 2; // 2-D space
        ypt += 2; // 2-D space
    }
    
    //fprintf(stderr, "\n");
    y = ((double) migrationCount) / ((double) TOTAL_N);
    
    if ( (t % TIME_SAMPLING_INTERVAL == 0) && !MULTIRUN ) {
        fprintf(migrationRecord, "%li,%f\n", t, y);
        fprintf( patchAbundances, "%li", t);
        for ( i = 0; i < nPATCHES; i++ ) {
            fprintf(patchAbundances, ",%i", nInEachPatch[i]);
        }
        fprintf(patchAbundances, "\n");
    }
    
    if ( migrationCount > 0 )
        sortPopulation();
    
    return y;
}


int multirunSetup(int seed)
{
    FILE *fpt, *fpt2;
    int rcount, count, i, j;
    
    // ensure new random number seed for next run
    fpt = fopen("RnumSeed.txt","w");
    fprintf(fpt,"%i\n",(seed+1));
    fclose(fpt);
    DETERMINISTIC = 1; // override to make sure seeds are all different in a multirun
    
    // index the multiple runs individually
    fpt = fopen("MultirunCount.txt","r");
    if (fpt == NULL) {
        fpt2 = fopen("MultirunCount.txt","w");
        fprintf(fpt2, "1\n");
        fclose(fpt2);
        count = 1;
    }
    else {
        rcount = fscanf(fpt, "%i", &count);
        //fprintf(stderr, "rcount = %i, count = %i\n", rcount, count);
        if ( rcount > 0 ) {
            count++;
            fclose(fpt);
        }
        else {
            fprintf(stderr, "\n\nWarning:\n\tnothing read from MultirunCount.txt file. Starting count at 1\n\n");
            count = 1;
        }
        fpt = fopen("MultirunCount.txt","w");
        fprintf(fpt, "%i\n", count);
        fclose(fpt);
    }
    fprintf(stderr, "\nRun %i started, seed = %i\n", count, seed);
    
    // data recording file
    fpt = fopen("MultirunEndpointData.csv","r");
    if (fpt == NULL) {
        fpt2 = fopen("MultirunEndpointData.csv","w");
        // add headers only when creating file anew
        // parameters
        fprintf(fpt2,"MultiRunIndex,RnumSeed,codeVersion,nLOCI,MAX_GENERATIONS,ROWS_OF_PATCHES,COLS_OF_PATCHES,PROPORTION_HOST_DARK_NICHE,");
        fprintf(fpt2,"PROPORTION_ADENOSTOMA,ACTUAL_PROPORTION_ADENOSTOMA,TOTAL_N,INIT_MEL_FREQ_ADENOSTOMA,");
        fprintf(fpt2,"INIT_MEL_FREQ_CEANOTHUS,nGENS_ALLOPATRY,SD_MOVE_MELANISTIC,SD_MOVE_GREEN,RANDOM_LOCS_FOR_OFFSPRING,");
        fprintf(fpt2,"COLOR_PATTERN_RECOMB_PROBABILITY,S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE,S_COEFF_GREEN_IN_DARK_NICHE,");
        fprintf(fpt2,"S_COEFF_MELANISTIC_IN_GREEN_NICHE,STRIPE_DOMINANCE_COEFFICIENT,BASE_MATING_FITNESS,CHC_MATCH_ADVANTAGE,");
        fprintf(fpt2,"MELANISTIC_MATING_ADVANTAGE_M,MELANISTIC_MATING_ADVANTAGE_F,INDEX_OF_COLOR_LOCUS,");
        fprintf(fpt2,"INDEX_OF_PATTERN_LOCUS,INDEX_OF_SEX_LOCUS,BIASED_NICHE_MATING_RATIO,LAYOUT_OF_HOSTS_RANDOM,");
        fprintf(fpt2,"HOLD_MEL_ALLELE_FREQ_CONSTANT,ADDITIVE_MODEL");
        
        for ( i = 1; i < nLOCI; i++ )
            fprintf(fpt2, ",RecombProbLoci%iand%i", (i-1), i );
        // results
        fprintf(fpt2, ",finalTimeAtEnd,InterHostMatingFreq,MigrationRate");
        for ( i = 0; i < nLOCI; i++ )
            fprintf(fpt2, ",AlleleFreqLocus%i", i);
        for ( i = 0; i < nLOCI; i++ )
            fprintf(fpt2, ",FSTlocus%i", i);
        for ( i = 0; i < nLOCI; i++ )
            fprintf(fpt2, ",FixationTimeLocus%i", i);
        for ( i = 0; i < (nLOCI-1); i++ ) {
            for ( j = (i+1); j < nLOCI; j++ ) {
                fprintf(fpt2, ",LDloci%iand%i", i, j);
            }
        }
        
        fprintf(fpt2,"\n");
        fclose(fpt2);
    }
    else
        fclose(fpt);
    
    return count;
}



void openDataRecordingFiles(void)
{
    int i, j;
    
    // data recording files
    migrationRecord = fopen("MigrationRecord.csv", "w");
    if ( PRINT_HEADERS )
        fprintf(migrationRecord, "Time,GrossMigrationRate\n"); // header line
    
    patchAbundances = fopen("PatchAbundances.csv", "w");
    if ( PRINT_HEADERS ) {
        fprintf(patchAbundances,"Time");
        for ( i = 0; i < nPATCHES; i++ ) {
            fprintf(patchAbundances, ",N%i", i);
        }
        fprintf( patchAbundances, "\n");
    }
    
    FSTtimeSeries = fopen("FSTtimeSeries.csv", "w");
    if ( PRINT_HEADERS ) {
        fprintf(FSTtimeSeries, "Time");
        for ( i = 0; i < nLOCI; i++ ) {
            fprintf(FSTtimeSeries, ",Locus%i", i);
        }
        fprintf( FSTtimeSeries, "\n");
    }
    
    AlFreqTS = fopen("AlleleFreqTimeSeries.csv", "w");
    if ( PRINT_HEADERS ) {
        fprintf( AlFreqTS, "Time");
        for ( i = 0; i < nLOCI; i++ ) {
            fprintf(AlFreqTS, ",Locus%i", i);
        }
        fprintf( AlFreqTS, "\n");
    }
    
    AlFreqByPatchTS = fopen("AlleleFreqByPatch.csv", "w");
    if ( PRINT_HEADERS ) {
        fprintf( AlFreqByPatchTS, "Time");
        for ( i = 0; i < nLOCI; i++ ) {
            for ( j = 0; j < nPATCHES; j++ ) {
                fprintf(AlFreqByPatchTS, ",Locus%iinPatch%i", i, j);
            }
        }
        fprintf(AlFreqByPatchTS, "\n");
    }
    
    LDTS = fopen("LDtimeSeries.csv", "w");
    if ( PRINT_HEADERS ) {
        fprintf(LDTS, "Time");
        for ( i = 0; i < (nLOCI-1); i++ ) {
            for ( j = i+1; j < nLOCI; j++ ) {
                fprintf(LDTS, ",Loci%iand%i", i, j);
            }
        }
        fprintf(LDTS, "\n");
    }
    
    interHostMating = fopen("InterHostMatingTS.csv","w");
    if ( PRINT_HEADERS )
        fprintf(interHostMating, "Time,interHostMatingFreq\n");
    
}


void printPopulation(char *fname)
{
    int i, j, *ipt;
    double *dpt;
    FILE *initPopn;
    char *phstr, *sexstr;
    
    initPopn = fopen( fname, "w" );
    if ( PRINT_HEADERS ) {
        fprintf(initPopn,"Individual,");
        for ( j = 0; j < nLOCI; j++ ) {
            
            
            if ( locusTypes[j] == NEUTRAL_LOCUS )
                fprintf(initPopn, "Locus%i_neutral_allele1,Locus%i_neutral_allele2,", j, j);
            else if ( locusTypes[j] == COLOR_LOCUS )
                fprintf(initPopn, "Locus%i_color_allele1,Locus%i_color_allele2,", j, j);
            else if ( locusTypes[j] == PATTERN_LOCUS )
                fprintf(initPopn, "Locus%i_pattern_allele1,Locus%i_pattern_allele2,", j, j);
            else if ( locusTypes[j] == SEX_LOCUS )
                fprintf(initPopn, "Locus%i_sexChrom_allele1,Locus%i_sexChrom_allele2,", j, j);
            //        else if ( locusTypes[i] == ASSORTMENT_LOCUS )
            //            fprintf(gsf, "ASSORTMENT_LOCUS,");
            else {
                fprintf(stderr, "\nError in initializeGenomeScheme():\n\ttype for locusTypes[i] (= %i) not recognized. Exiting\n", locusTypes[i]);
                exit(-1);
            }
            
        }
        fprintf(initPopn, "phenotypeCode,phenotype,sexCode,sex,patch,xcoordinate,ycoordinate,niche,hostType,chcType\n");
    }
    
    ipt = genotypes;
    dpt = coordinates;
    for ( i = 0; i < TOTAL_N; i++ ) {
        fprintf(initPopn, "%i", i);
        for ( j = 0; j < nLOCI; j++ ) {
            fprintf(initPopn, ",%i,%i", (*ipt), (*(ipt+1)));
            ipt += DIPLOID;
        }
        
        
        if ( ADDITIVE_MODEL ) {
            if ( phenotypes[i] == PH_GU )
                phstr = "PH_GU";
            else if ( phenotypes[i] == PH_GH )
                phstr = "PH_GH";
            else if ( phenotypes[i] == PH_GS )
                phstr = "PH_GS";
            else if ( phenotypes[i] == PH_HU )
                phstr = "PH_HU";
            else if ( phenotypes[i] == PH_HH )
                phstr = "PH_HH";
            else if ( phenotypes[i] == PH_HS )
                phstr = "PH_HS";
            else if ( phenotypes[i] == PH_MU )
                phstr = "PH_MU";
            else if ( phenotypes[i] == PH_MH )
                phstr = "PH_MH";
            else if ( phenotypes[i] == PH_MS )
                phstr = "PH_MS";
            else {
                fprintf(stderr, "\nError in printPopulation():\n\tphenotypes[%i] (= %i) code not found. Exiting\n", i, phenotypes[i]);
                exit(-1);
            }

        }
        else {
            if ( phenotypes[i] == PH_MELANISTIC )
                phstr = "PH_MELANISTIC";
            else if ( phenotypes[i] == PH_STRIPED_GREEN_HOMO )
                phstr = "PH_STRIPED_GREEN_HOMO";
            else if ( phenotypes[i] == PH_UNSTRIPED_GREEN_HOMO )
                phstr = "PH_UNSTRIPED_GREEN_HOMO";
            else if ( phenotypes[i] == PH_HETERO_STRIPE )
                phstr = "PH_HETERO_STRIPE";
            else {
                fprintf(stderr, "\nError in printPopulation():\n\tphenotypes[%i] (= %i) code not found. Exiting\n", i, phenotypes[i]);
                exit(-1);
            }
        }
        
        if ( sexOfIndividual[i] == SEX_FEMALE )
            sexstr = "SEX_FEMALE";
        else if ( sexOfIndividual[i] == SEX_MALE )
            sexstr = "SEX_MALE";
        else {
            fprintf(stderr, "\nError in printPopulation():\n\tsexOfIndividual[%i] (= %i) code not found. Exiting\n", i, sexOfIndividual[i]);
            exit(-1);
        }
        
        fprintf(initPopn, ",%i,%s,%i,%s,%i,%f,%f,%i,%i,%i\n", phenotypes[i], phstr, sexOfIndividual[i], sexstr, patchLocations[i], (*dpt), (*(dpt + 1)), nicheWithinPatch[i], hostType[(patchLocations[i])], chcType[i]);
        dpt += 2; // 2-D space
    }
    
    fclose(initPopn);
}



double reproduction(void)
{
    int i, j, startIndexes[nPATCHES], endIndexes[nPATCHES];
    int si, ei, n, pheno, patchType, momi, dadi, dumi;
    double fitnessArray[TOTAL_N], fitVal, fitSum, paccept, crossMatingFreq;
    int *offspringPhenotypes, *offPatchLocations, *offspringGenotypes;
    int *offspringNiches, *offspringSexes, *offspringChcType;
    double *offspringCoordinates, *dpt1, momf, dadf, *fapt, matingReduction;
    int *startpt, *iptGTs, *iptPHs, *iptPatLocs, *iptOffNiche, *iptsex, *chcpt;
    int row, col, totalCount, maleCount, femaleCount;
    int maleIndexes[TOTAL_N], femaleIndexes[TOTAL_N], crossMatingCount, dadNiche;
    int maleOffspringCount = 0, femaleOffspringCount = 0;
    
    offspringGenotypes = (int *) malloc( (sizeof(int) * TOTAL_N * DIPLOID * nLOCI) ); // diploid genotypes
    offspringPhenotypes = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    offPatchLocations = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    offspringCoordinates = (double *) malloc( (sizeof(double) * 2 * TOTAL_N) ); // 2-D euclidean coordinates
    offspringNiches = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    offspringSexes = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    offspringChcType = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    
    dpt1 = offspringCoordinates;
    iptGTs = offspringGenotypes;
    iptPHs = offspringPhenotypes;
    iptPatLocs = offPatchLocations;
    iptOffNiche = offspringNiches;
    iptsex = offspringSexes;
    chcpt = offspringChcType;
    // NOTE: markedDead not here because all marked dead will be reset to zero in this function
    
    fapt = &fitnessArray[0];
    
    crossMatingCount = 0;
    si = 0;
    ei = nInEachPatch[0] - 1;
    totalCount = 0;
    for ( i = 0; i < nPATCHES; i++ ) {
        
//        dadWarningCount = 0;
//        momWarningCount = 0;
        
        n = nInEachPatch[i];
        // first check
        if ( n < 2 ) {
            fprintf(stderr, "\nError in reproduction()!\n\tn = %i < 2 individuals in patch %i\n", n, i);
            exit(-1);
        }
        patchType = hostType[i];
        
        // build separate arrays of males and females indexes in this patch
        maleCount = 0;
        femaleCount = 0;
        for ( j = si; j <= ei; j++ ) {
            // a check on things
            if ( patchLocations[j] != i ) {
                fprintf(stderr, "\nError in reproduction()!\n\tpatchLocations[%i] = %i does not match expected patch %i.\n\tAlgorithm busted somewhere!\n", j, patchLocations[j], i);
                exit(-1);
            }
            // build male and female arrays
            if ( sexOfIndividual[j] == SEX_MALE ) {
                maleIndexes[maleCount] = j;
                maleCount++;
            }
            else if ( sexOfIndividual[j] == SEX_FEMALE ) {
                femaleIndexes[femaleCount] = j;
                femaleCount++;
            }
            else {
                fprintf(stderr, "\nError in reproduction():\n\tSex not found!\n");
                exit(-1);
            }
        }
        // another check on things
        if ( !maleCount || !femaleCount ) {
            if ( !maleCount )
                fprintf(stderr, "\nError in reproduction()!\n\tNo males in patch %i at time %li\n", i, t);
            if ( !femaleCount )
                fprintf(stderr, "\nError in reproduction()!\n\tNo females in patch %i at time %li\n", i, t);
            exit(-1);
        }
        if ( (maleCount + femaleCount) != n ) {
            fprintf(stderr, "\nError in reproduction()!\n\tmaleCount (%i) + femaleCount (%i) != n (%i)\n", maleCount, femaleCount, n);
            exit(-1);
        }
        
        // location in 2-D discrete patch space
        row = i / COLS_OF_PATCHES;
        col = i - ( row * COLS_OF_PATCHES );
        
        // make n offspring
        for ( j = 0; j < n; j++ ) {
            
            
            //chooseParents1( &momi, &dadi, &maleIndexes[0], &femaleIndexes[0], maleCount, femaleCount );
            
            chooseParents2( &momi, &dadi, &maleIndexes[0], &femaleIndexes[0], maleCount, femaleCount );
            if ( ( *(patchLocations + momi) != i ) || ( *(patchLocations + dadi) != i ) || momi == dadi ) {
                fprintf(stderr, "\nError in reproduction() --> chooseParents():\n\tyou chose bad parents!");
                fprintf(stderr, "\n\tmomi = %i, dadi = %i, patchLocations[momi] = %i, patchLocations[dadi] = %i, patch = %i\n", momi, dadi, *(patchLocations + momi), *(patchLocations + dadi), i);
                exit(-1);
            }
            
            
            
            if ( chcType[momi] != chcType[dadi] ) {
                crossMatingCount++;
            }

            // get gamete (alleles) from parent #1 -- first haplotype
            getOffspringAlleles(momi, iptGTs);
            
            // get gamete (alelles) from parent #2 -- second haplotype
            getOffspringAlleles(dadi, (iptGTs + 1));
            // set genotype is now done -- 1
            
            // set offspring coordinates -- 2
            if ( RANDOM_LOCS_FOR_OFFSPRING ) {
                *dpt1 = randU() + ((double) col); // "x" coordinate
                *(dpt1 + 1) = randU() + ((double) row); // "y" coordinate
            }
            else {
                *dpt1 = *(coordinates + (2 * momi)); // "x" coordinate
                *(dpt1 + 1) = *(coordinates + (2 * momi) + 1); // "y" coordinate
            }
            
            // set niche -- 3
            *iptOffNiche = findNiche( (*dpt1), (*(dpt1+1)) );
            
            // set offspring phenotype -- 4
            *iptPHs = computePhenotype( iptGTs );
            
            // offspring sex -- 5
            *iptsex = computeSex( iptGTs );
            if ( (*iptsex) == SEX_MALE )
                maleOffspringCount++;
            else
                femaleOffspringCount++;
            
            // chc type -- 6
            *chcpt = patchType;
            
            // record patch location -- 7
            *iptPatLocs = i;
            
            // increment pointers
            iptGTs += ( DIPLOID * nLOCI );
            iptPHs++;
            iptPatLocs++;
            iptOffNiche++;
            dpt1 += 2;
            iptsex++;
            chcpt++;
            
            totalCount++;
        }
        
//        if ( dadWarningCount > 0 ) {
//            fprintf(stderr, "\nWarning!  had to resurrect dead male(s) %i times to escape infinite loop!", dadWarningCount);
//            fprintf(stderr, "\nt = %li, maleCount = %i, nInEachPatch[%i] = %i\n", t, maleCount, i, nInEachPatch[i]);
//        }
//        if ( momWarningCount > 0 ) {
//            fprintf(stderr, "\nWarning!  had to resurrect dead female(s) %i times to escape infinite loop!", momWarningCount);
//            fprintf(stderr, "\nt = %li, femaleCount = %i, nInEachPatch[%i] = %i\n!", t, femaleCount, i, nInEachPatch[i]);
//        }
        
        
        if ( i < (nPATCHES-1) ) {
            si = ei + 1;
            ei += nInEachPatch[(i+1)];
            if ( (patchLocations[si] != patchLocations[ei]) || (patchLocations[si] != (i+1)) ) {
                fprintf(stderr, "\nError in reproduction(): Bad indexing!  Exiting!\n");
                fprintf(stderr, "si = %i, ei = %i, nInEachPatch[(i+1)] = %i\n", si, ei, nInEachPatch[(i+1)]);
                exit(-1);
            }
        }
        
    }
    
    // error checking
    if ( totalCount != TOTAL_N ) {
        fprintf(stderr, "\nError in reproduction()!  totalCount (%i) != TOTAL_N (%i)\n", totalCount, TOTAL_N);
        exit(-1);
    }
    if ( (maleOffspringCount < (TOTAL_N / 4)) || (femaleOffspringCount < (TOTAL_N / 4)) ) {
        fprintf(stderr, "\nWarning in reproduction():\n\tHighly biased sex ratio:\n\tnmaleOffspringCount = %i, femaleOffspringCount = %i\n", maleOffspringCount, femaleOffspringCount);
        exit(-1);
    }
    
    crossMatingFreq = ((double) crossMatingCount) / ((double) totalCount);
    if ( (t % TIME_SAMPLING_INTERVAL == 0) && !MULTIRUN )
        fprintf(interHostMating, "%li,%f\n", t, crossMatingFreq);
    
    // replace parents with offspring
    free(genotypes);
    genotypes = offspringGenotypes;
    
    free(phenotypes);
    phenotypes = offspringPhenotypes;
    
    free(patchLocations);
    patchLocations = offPatchLocations;
    
    free(coordinates);
    coordinates = offspringCoordinates;
    
    free(nicheWithinPatch);
    nicheWithinPatch = offspringNiches;
    
    free(sexOfIndividual);
    sexOfIndividual = offspringSexes;
    
    free(chcType);
    chcType = offspringChcType;
    
    // set markedDead -- 8
    for ( i = 0; i < TOTAL_N; i++ )
        markedDead[i] = 0;
    
    if ( HOLD_MEL_ALLELE_FREQ_CONSTANT )
        setMelanisticAlleleFrequency();
    
    return crossMatingFreq;
}


int RNGsetup(void)
{
    int seed, rcount, count;
    FILE *fpt, *fpt2;
    int stime;
    long ltime;
    FILE *rseed;
    long int i;

    
    if (DETERMINISTIC) {
        
        fpt = fopen("RnumSeed.txt","r");
        if (fpt == NULL) {
            perror("Can't open RnumSeed.txt");
            exit(-1);
        }
        
        rcount = fscanf(fpt,"%i",&seed);
        if ( rcount ) {
			seedRand(seed);	// fixed random number seed
			fclose(fpt);
		}
		else {
			fprintf(stderr, "\n\n\tError! nothing read from file! Exiting!\n\n");
			exit(-1);
		}
        //fprintf(stderr, "\n\nSeed = %i\n\n",seed);
    }
    else {
        /* use calendar time to seed random number generator.  Code adopted from Schildt's textbook */
        
        /* get the calendar time */
        ltime=time(NULL);
        stime=(unsigned) ltime/2;
        
        // generate and store random number seed
        rseed = fopen("RnumSeed.txt","w");
        fprintf(rseed,"%i\n",stime);
        fclose(rseed);
        seedRand(stime);	// get random number seed (system time)
        seed = stime;
    }
    
    // warm up
    for ( i = 0; i < 1000000; i++ )
        randU();
    
    return seed;
}


void setMelanisticAlleleFrequency(void)
{
    int i, p, *iptGTs, *iptPatchLoc, *iptIndivd;
    double dum, a0, ah, c0, ch;
    
    a0 = INIT_MEL_FREQ_ADENOSTOMA * INIT_MEL_FREQ_ADENOSTOMA; // p^2
    ah = 2.0 * INIT_MEL_FREQ_ADENOSTOMA * (1.0 - INIT_MEL_FREQ_ADENOSTOMA); // 2pq
    ah = ah + a0; // 2pq + p^2
    if ( ah > 1.0 || a0 < 0.0 || ah < a0 ) {
        fprintf(stderr, "\nError in setMelanisticAlleleFrequency():\n\ta0 = %f, ah = %f\n", a0, ah);
        exit(-1);
    }
    
    c0 = INIT_MEL_FREQ_CEANOTHUS * INIT_MEL_FREQ_CEANOTHUS;
    ch = 2.0 * INIT_MEL_FREQ_CEANOTHUS * ( 1.0 - INIT_MEL_FREQ_CEANOTHUS );
    ch = ch + c0;
    if ( ch > 1.0 || c0 < 0.0 || ch < c0 ) {
        fprintf(stderr, "\nError in setMelanisticAlleleFrequency():\n\tc0 = %f, ch = %f\n", c0, ch);
        exit(-1);
    }
    
    iptIndivd = genotypes;
    iptGTs = genotypes + (DIPLOID * INDEX_OF_COLOR_LOCUS);
    iptPatchLoc = patchLocations;
    
    for ( i = 0; i < TOTAL_N; i++ ) {
        p = *iptPatchLoc;
        dum = randU();
        
        if ( hostType[p] == ADENOSTOMA_HOST ) {
            if ( dum < a0 ) {
                *iptGTs = MELANISTIC_ALLELE;
                *(iptGTs + 1) = MELANISTIC_ALLELE;
            }
            else if ( dum < ch ) {
                if ( randU() < 0.5 ) {
                    *iptGTs = MELANISTIC_ALLELE;
                    *(iptGTs + 1) = GREEN_ALLELE;
                }
                else {
                    *iptGTs = GREEN_ALLELE;
                    *(iptGTs + 1) = MELANISTIC_ALLELE;
                }
            }
            else {
                *iptGTs = GREEN_ALLELE;
                *(iptGTs + 1) = GREEN_ALLELE;
            }
        }
        
        else if ( hostType[p] == CEANOTHUS_HOST ) {
            if ( dum < c0 ) {
                *iptGTs = MELANISTIC_ALLELE;
                *(iptGTs + 1) = MELANISTIC_ALLELE;
            }
            else if ( dum < ch ) {
                if ( randU() < 0.5 ) {
                    *iptGTs = MELANISTIC_ALLELE;
                    *(iptGTs + 1) = GREEN_ALLELE;
                }
                else {
                    *iptGTs = GREEN_ALLELE;
                    *(iptGTs + 1) = MELANISTIC_ALLELE;
                }
            }
            else {
                *iptGTs = GREEN_ALLELE;
                *(iptGTs + 1) = GREEN_ALLELE;
            }
        }
        
        else {
            fprintf(stderr, "\nError in setMelanisticAlleleFrequency():\n\tHost type not found!\n");
            exit(-1);
        }
        
        phenotypes[i] = computePhenotype( iptIndivd );
        
        iptGTs += ( DIPLOID * nLOCI );
        iptPatchLoc++;
        iptIndivd += ( DIPLOID * nLOCI );
    }
    
}


void sortPopulation(void)
{
    int i, j, count, currentPatch, index, *newMarkedDead, *newNicheWithinPatch;
    int *newGenotypes, *newPhenotypes, *newPatchLocations, *newChcType, *newSexOfIndividual;
    double *newCoordinates, *dptNew, *dptOld;
    int indexesToUse[nPATCHES], *iptNew, *ipgOld, *ipphOld, *ippaOld;
    int *ipMarkedDeadOld, *ipNicheWithinPatchOld, *ipChcOld, *ipSexOld;
    size_t copySize;
    
    copySize = ( sizeof(int) * DIPLOID * nLOCI ); // genotype copying size for memcpy() commands
    
    indexesToUse[0] = 0;
    for ( i = 1; i < nPATCHES; i++ ) {
        indexesToUse[i] = indexesToUse[(i-1)] + nInEachPatch[(i-1)];
    }
    
    newGenotypes = (int *) malloc( (sizeof(int) * TOTAL_N * 2 * nLOCI) ); // diploid genotypes
    newPhenotypes = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    newPatchLocations = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    newCoordinates = (double *) malloc( (sizeof(double) * 2 * TOTAL_N) ); // 2-D euclidean coordinates
    newMarkedDead = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    newNicheWithinPatch = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    newChcType = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    newSexOfIndividual = (int *) malloc( (sizeof(int) * TOTAL_N) ); // scalar designation
    
    ipgOld = genotypes;
    ipphOld = phenotypes;
    ippaOld = patchLocations;
    dptOld = coordinates;
    ipMarkedDeadOld = markedDead;
    ipNicheWithinPatchOld = nicheWithinPatch;
    ipChcOld = chcType;
    ipSexOld = sexOfIndividual;
    
    for ( i = 0; i < TOTAL_N; i++ ) {
        currentPatch = patchLocations[i];
        index = indexesToUse[currentPatch];
        
        // genotypes -- 1
        iptNew = newGenotypes + ( DIPLOID * nLOCI * index );
        memcpy( iptNew, ipgOld, copySize ); // dest, source, size
        ipgOld += (DIPLOID * nLOCI);
        
        // phenotypes -- 2
        iptNew = newPhenotypes + index;
        *iptNew = *ipphOld;
        ipphOld++;
        
        // niche -- 3
        iptNew = newNicheWithinPatch + index;
        *iptNew = *ipNicheWithinPatchOld;
        ipNicheWithinPatchOld++;
        
        // marked dead -- 4
        iptNew = newMarkedDead + index;
        *iptNew = *ipMarkedDeadOld;
        ipMarkedDeadOld++;
        
        // patch location -- 5
        iptNew = newPatchLocations + index;
        *iptNew = *ippaOld;
        ippaOld++;
        
        // coordinates -- 6
        dptNew = newCoordinates + (2 * index);
        *dptNew = *dptOld;
        *(dptNew+1) = *(dptOld+1);
        dptOld += 2;
        
        // chc -- 7
        iptNew = newChcType + index;
        *iptNew = *ipChcOld;
        ipChcOld++;
        
        // sex -- 8
        iptNew = newSexOfIndividual + index;
        *iptNew = *ipSexOld;
        ipSexOld++;
        
        indexesToUse[currentPatch] = indexesToUse[currentPatch] + 1;
    }
    
    // error checking
    count = 0;
    for ( i = 0; i < nPATCHES; i++ ) {
        count += nInEachPatch[i];
        if ( indexesToUse[i] != count ) {
            fprintf(stderr, "\nError!  Bad indexing in sort!  Exiting!\n\n");
            fprintf(stderr, "t = %li, i = %i, indexesToUse[i] = %i, nInEachPatch[i] = %i\n\n", t, i, indexesToUse[i], nInEachPatch[i]);
            exit(-1);
        }
    }

    free(genotypes);
    genotypes = newGenotypes;
    
    free(phenotypes);
    phenotypes = newPhenotypes;
    
    free(patchLocations);
    patchLocations = newPatchLocations;
    
    free(coordinates);
    coordinates = newCoordinates;
    
    free(markedDead);
    markedDead = newMarkedDead;
    
    free(nicheWithinPatch);
    nicheWithinPatch = newNicheWithinPatch;
    
    free(sexOfIndividual);
    sexOfIndividual = newSexOfIndividual;
    
    free(chcType);
    chcType = newChcType;
}



void usage(char *s)
{
    fprintf(stderr,  "\nUsage: %s [options]\n\n", s); // info print
	   
    
	fprintf(stderr,  "\tNOTE:  There is no error checking on user inputs.\n");
	fprintf(stderr,  "\tProgram behavior is unpredictable (not defined) for\n");
	fprintf(stderr,  "\tinputs that do not conform to guidelines below.\n\n");
    
    fprintf(stderr,  "\tN.B.:  Options are CaSe SeNsItIvE.\n\n");
    fprintf(stderr,  "\tN.B.:  Parameter names in ALL_CAPS denote actual names of these\n\t\tparameters as represented in the source code.\n\n");
    fprintf(stderr,  "\tN.B.:  You are going to do something #badass with this program.\n\n");
    fprintf(stderr,  "\tOK, enough with the #overkillMessages and #gratuitousHashTags!\n\tHere are the actual options you can call on the command line!\n\n");
	
	fprintf(stderr,  "\t[-a <proportion>]\n\t\tHost type proportions, a number between 0 and 1 that\n");
	fprintf(stderr,  "\t\trepresents the PROPORTION_ADENOSTOMA.\n");
	fprintf(stderr,  "\t\tThe default value for this proportion is %.3f\n", PROPORTION_ADENOSTOMA_DEFAULT);
    fprintf(stderr,  "\t\t(i.e., %.0f%% Adenostoma, %.0f%% Ceanothus).\n\n", (100.0*PROPORTION_ADENOSTOMA_DEFAULT), (100.0*(1.0 - PROPORTION_ADENOSTOMA_DEFAULT)));
    
    fprintf(stderr,  "\t[-A <generations>]\n\t\tnGENS_ALLOPATRY, the number of generations of zero migration.\n");
	fprintf(stderr,  "\t\tMust be an integer.  ");
	fprintf(stderr,  "The default value is %i\n\n", nGENS_ALLOPATRY_DEFAULT);
    
    fprintf(stderr,  "\t[-b <0,1>]\n\t\tHOLD_MEL_ALLELE_FREQ_CONSTANT or not.\n\t\t-b 0 causes the frequency of the melanistic\n\t\tallele to be able to change according to all usual\n\t\tevolutionary forces like other alleles.\n\t\t-b 1 causes the frequency to be held constant.\n");
    fprintf(stderr,  "\t\tOption must be 0 or 1.  ");
    fprintf(stderr,  "The default value is %i\n\n", HOLD_MEL_ALLELE_FREQ_CONSTANT_DEFAULT);
    
    fprintf(stderr,  "\t[-B <value>]\n\t\tBIASED_NICHE_MATING_RATIO, expressed as the ratio of mating\n\t\twith other niche to mating within your niche.\n\t\tA value of 1.0 means no bias; a value < 1.0 means bias toward\n\t\tmating within niche.\n\t\tWith the current code setup, the chosen value\n");
    fprintf(stderr,  "\t\tshould always lie between 0.0 and 1.0\n\t\t");
    fprintf(stderr,  "The default value is %f\n\n", BIASED_NICHE_MATING_RATIO_DEFAULT);
    
    fprintf(stderr,  "\t[-c <probability>]\n\t\tThe cost of migration, measured as MIGRATION_DEATH_PROBABILITY.\n");
	fprintf(stderr,  "\t\tMust be between 0 and 1.  ");
	fprintf(stderr,  "The default value is %f\n\n", MIGRATION_DEATH_PROBABILITY_DEFAULT);
    
    fprintf(stderr,  "\t[-C <integer>]\n\t\tCOLS_OF_PATCHES, the number of columns of discrete patches in\n\t\tthe entire habitat, where said habitat is viewed as a\n\t\trectangular grid of patches (see also -R below).\n");
	fprintf(stderr,  "\t\tMust be an integer >= 1.");
	fprintf(stderr,  "\n\t\tThe default value is %i\n\n", COLS_OF_PATCHES_DEFAULT);
    
    fprintf(stderr,  "\t[-d <proportion>]\n\t\tPROPORTION_HOST_DARK_NICHE, the proportion of a\n");
	fprintf(stderr,  "\t\thost plant that is a niche for the melanistic morph.");
	fprintf(stderr,  "\n\t\tValue should be between 0 and 1.  The default value is %.3f\n\n", PROPORTION_HOST_DARK_NICHE_DEFAULT);
    
    fprintf(stderr,  "\t[-D <0,1>]\n\t\tDETERMINISTIC mode.\n\t\t-D 0 causes the system to generate a new random number seed\n\t\tand store it in a file called RnumSeed.txt.\n\t\t-D 1 causes the program to read in a random number seed from an\n\t\texisting RnumSeed.txt file.\n");
	fprintf(stderr,  "\t\tMust be 0 or 1.  ");
	fprintf(stderr,  "The default value is %i\n\n", DETERMINISTIC_DEFAULT);

    
//case 'f':
//    BASE_MATING_FITNESS = strtod(optarg, (char **)NULL);
//    break;
//case 'F':
//    CHC_MATCH_ADVANTAGE = strtod(optarg, (char **)NULL);
//    break;
//case 'g':
//    MELANISTIC_MATING_ADVANTAGE = strtod(optarg, (char **)NULL);
    
    
    fprintf(stderr,  "\t[-f,-F,-g,-G <value>]\n\t\tRelative mating success values for phenotypes.\n");
    fprintf(stderr,  "\t\tNote that what really matters here is the RELATIVE value of\n");
    fprintf(stderr,  "\t\teach, because the probability of things working out is\n");
    fprintf(stderr,  "\t\tNORMALIZED by a SUM of relevant parameters to have a\n\t\tmaximum value of 1.\n");
    fprintf(stderr,  "\t\tOption\tCorresponding parameter\t\tDefault value\n");
    fprintf(stderr,  "\t\t-f\tBASE_MATING_FITNESS\t\t%.4f\n", BASE_MATING_FITNESS_DEFAULT);
    fprintf(stderr,  "\t\t-F\tCHC_MATCH_ADVANTAGE\t\t%.4f\n", CHC_MATCH_ADVANTAGE_DEFAULT);
	fprintf(stderr,  "\t\t-g\tMELANISTIC_MATING_ADVANTAGE_M\t%.4f\n", MELANISTIC_MATING_ADVANTAGE_M_DEFAULT);
    fprintf(stderr,  "\t\t-G\tMELANISTIC_MATING_ADVANTAGE_F\t%.4f\n", MELANISTIC_MATING_ADVANTAGE_F_DEFAULT);
    fprintf(stderr,  "\t\tMale mating fitness is normalized by the sum of f + g.\n");
    fprintf(stderr,  "\t\tFemale mating fitness is normalized by the sum of f + F + G.\n\n");
    
    fprintf(stderr,  "\t[-h <0,1>]\n\t\tPRINT_HEADERS in output .csv data files.\n\t\t-h 0 causes the output .csv data files to have NO headers.\n\t\t-h 1 causes the output data files to have headers.\n");
	fprintf(stderr,  "\t\tMust be 0 or 1.  ");
	fprintf(stderr,  "The default value is %i\n\n", PRINT_HEADERS_DEFAULT);
    
    fprintf(stderr,  "\t[-H <coeff>]\n\t\tSTRIPE_DOMINANCE_COEFFICIENT of the UNSTRIPED_ALLELE over the\n\t\tSTRIPED_ALLELE at the color locus.  ");
	fprintf(stderr,  "The value chosen should be\n\t\tbetween 0.5 and 1.  A value of 1 means that the UNSTRIPED_ALLELE\n\t\tis totally dominant; a value of 0.5 is perfect codominance.");
	fprintf(stderr,  "\n\t\tThe default value is %f\n\n", STRIPE_DOMINANCE_COEFFICIENT_DEFAULT);
    
    fprintf(stderr,  "\t[-L <0,1>]\n\t\tLAYOUT_OF_HOSTS_RANDOM or not. \n\t\t-L 0 causes patches of Adenostoma to be clumped with one\n\t\tanother and likewise for Ceanothus.\n\t\t-L 1 causes patches to be randomly located with respect to\n\t\thost type.  ");
    fprintf(stderr,  "Must be 0 or 1.  ");
    fprintf(stderr,  "\n\t\tThe default value is %i\n\n", LAYOUT_OF_HOSTS_RANDOM_DEFAULT);
    
    fprintf(stderr,  "\t[-m,-M <value>]\n\t\tStandard deviation of migration distance for a given phenotype:\n");
    fprintf(stderr,  "\t\tOption\tCorresponding parameter\t\tDefault value\n");
    fprintf(stderr,  "\t\t-m\tSD_MOVE_MELANISTIC\t\t%.3f\n", SD_MOVE_MELANISTIC_DEFAULT);
    fprintf(stderr,  "\t\t-M\tSD_MOVE_GREEN\t\t\t%.3f\n\n", SD_MOVE_GREEN_DEFAULT);

    fprintf(stderr,  "\t[-N <integer>]\n\t\tTOTAL_N, the population size of all patches (demes) combined.");
	fprintf(stderr,  "\n\t\tMust be an integer.  ");
	fprintf(stderr,  "The default value is %i\n\n", TOTAL_N_DEFAULT);
    
    fprintf(stderr,  "\t[-o <0,1>]\n\t\tRANDOM_LOCS_FOR_OFFSPRING.\n\t\t-o 1 causes offspring in a deme to have random locations within\n\t\tthat deme with respect to where their parents were.\n\t\t-o 0 causes an offspring to be 'born' at the location of one of\n\t\tits parents.\n\t\t\tThis matters because while patches (demes)\n\t\tare treated as discrete for the purposes of reproduction,\n\t\tspace within patches is indeed treated as continuous.\n");
	fprintf(stderr,  "\t\tMust be 0 or 1.  ");
	fprintf(stderr,  "The default value is %i\n\n", RANDOM_LOCS_FOR_OFFSPRING_DEFAULT);
    
    fprintf(stderr,  "\t[-p <proportion>]\n\t\tINIT_MEL_FREQ_ADENOSTOMA, the initial frequency of the\n");
	fprintf(stderr,  "\t\tmelanistic allele in Adenostoma patches (demes).");
	fprintf(stderr,  "\n\t\tValue should be between 0 and 1.  The default value is %.3f\n\n", INIT_MEL_FREQ_ADENOSTOMA_DEFAULT);
    
    fprintf(stderr,  "\t[-P <proportion>]\n\t\tINIT_MEL_FREQ_CEANOTHUS, the initial frequency of the\n");
	fprintf(stderr,  "\t\tmelanistic allele in Ceanothus patches (demes).");
	fprintf(stderr,  "\n\t\tValue should be between 0 and 1.  The default value is %.3f\n\n", INIT_MEL_FREQ_CEANOTHUS_DEFAULT);
    
    fprintf(stderr,  "\t[-r <proportion>]\n\t\tCOLOR_PATTERN_RECOMB_PROBABILITY, i.e., proportion of meioses in\n");
	fprintf(stderr,  "\t\twhich parentally inherited combinations of color and pattern\n\t\talleles at those two loci are broken up.  0.5 means they are");
	fprintf(stderr,  "\n\t\tunlinked.  The value should be between 0 and 0.5.\n\t\tThe default value is %.3f\n\n", COLOR_PATTERN_RECOMB_PROBABILITY_DEFAULT);
    
    fprintf(stderr,  "\t[-R <integer>]\n\t\tROWS_OF_PATCHES, the number of rows of discrete patches in\n\t\tthe entire habitat, where said habitat is viewed as a\n\t\trectangular grid of patches (see also -C above).\n");
	fprintf(stderr,  "\t\tMust be an integer > 1 (so that there are at least two patches).");
	fprintf(stderr,  "\n\t\tThe default value is %i\n\n", ROWS_OF_PATCHES_DEFAULT);
    
    fprintf(stderr,  "\t[-s,-S,-t <value>]\n\t\tSelection coefficients for a given phenotype on a given host\n\t\tand niche:\n");
    fprintf(stderr,  "\t\tOption\tCorresponding parameter\t\t\tDefault value\n");
    fprintf(stderr,  "\t\t-s\tS_COEFF_WRONG_STRIPE_IN_GREEN_NICHE\t%.2f\n", S_COEFF_WRONG_STRIPE_DEFAULT);
    fprintf(stderr,  "\t\t-S\tS_COEFF_MELANISTIC_IN_GREEN_NICHE\t%.2f\n", S_COEFF_MELANISTIC_IN_GREEN_DEFAULT);
    fprintf(stderr,  "\t\t-t\tS_COEFF_GREEN_IN_DARK_NICHE\t\t%.2f\n\n", S_COEFF_GREEN_IN_DARK_DEFAULT);
    
    fprintf(stderr,  "\t[-T <integer>]\n\t\tMAX_GENERATIONS, the maximum number of generations to run.");
	fprintf(stderr,  "\n\t\tMust be an integer.  ");
	fprintf(stderr,  "The default value is %i\n\n", MAX_GENERATIONS_DEFAULT);
    
    fprintf(stderr,  "\t[-w <0,1>]\n\t\tMULTIRUN mode.\n\t\t-w 0 causes the program to run like only a single\n\t\trun of the simulation will occur, writing full time series\n\t\tof data.\n\t\t-w 1 causes the program to run as if it is being called by a\n\t\twrapper that may script multiple runs (in sequence, not\n\t\tsimultaneously).  For this, only data on ending states will be\n\t\twritten, and those data will be appended to files such that a\n\t\tsingle data file may contain data on multiple independent\n\t\tinstances of the simulation.  Different random number seeds\n\t\twill be automatically generated with -w 1.\n");
	fprintf(stderr,  "\t\tMust be 0 or 1.  ");
	fprintf(stderr,  "The default value is %i\n\n", MULTIRUN_DEFAULT);
    
    fprintf(stderr,  "\t[-z <0,1>]\n\t\tuse a purely ADDITIVE_MODEL of fitness.");
    fprintf(stderr,  "\n\t\t0 causes the dominance and epistasis known from the real\n\t\tpopulation.  ");
    fprintf(stderr,  "\n\t\t1 causes a purely additive model of allele contributions to\n\t\tphenotypes and fitness.  ");
    fprintf(stderr,  "\n\t\tThe default value is %i\n\n", ADDITIVE_MODEL_DEFAULT);
    
    fprintf(stderr,  "\t[-Z <integer>]\n\t\tTIME_SAMPLING_INTERVAL, the frequency of data recording");
    fprintf(stderr,  "\n\t\texpressed in numbers of generations.  Only matters in regular\n\t\tmode (not multirun mode).  ");
    fprintf(stderr,  "Must be an integer > 0.\n\t\tThe default value is %i\n\n", TIME_SAMPLING_INTERVAL_DEFAULT);
    
    fprintf(stderr,  "\t[-?]\n\t\tPrint this list of options to screen and exit the program.\n\n");
    
}


void viabilitySelection(void)
{
    int i, j, phi, nichei, hosti;
    int *gtpt, patternGT, homoS, het, homoU;
    
    homoS = STRIPED_ALLELE + STRIPED_ALLELE;
    het = STRIPED_ALLELE + UNSTRIPED_ALLELE;
    homoU = UNSTRIPED_ALLELE + UNSTRIPED_ALLELE;
    
    gtpt = genotypes + (DIPLOID * INDEX_OF_PATTERN_LOCUS);
    
    for ( i = 0; i < TOTAL_N; i++ ) {
        
        phi = phenotypes[i];
        nichei = nicheWithinPatch[i];
        hosti = hostType[ (patchLocations[i]) ];
        
        
        if ( ADDITIVE_MODEL ) {
            if ( randU() < viabilityMatrix[phi][hosti][nichei] )
                markedDead[i] = 1;
        }
        else {
            if ( phi == PH_MELANISTIC ) { // melanistic
                if ( nichei == NICHE_GREEN ) {
                    if ( randU() < S_COEFF_MELANISTIC_IN_GREEN_NICHE ) {
                        markedDead[i] = 1;
                    }
                }
            }
            else if ( nichei == NICHE_DARK ) { // green of either phenotype in dark niche
                if ( randU() < S_COEFF_GREEN_IN_DARK_NICHE ) {
                    markedDead[i] = 1;
                }
            }
            else {
                // got a green bug in a green niche, but fitness now depends upon host type and bugs phenotype
                patternGT = (*gtpt) + (*(gtpt+1));
                if ( patternGT == homoS && hosti == CEANOTHUS_HOST ) {
                    // striped homozygote; bad on ceanothus
                    if ( randU() < S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE ) {
                        markedDead[i] = 1;
                    }
                }
                else if ( patternGT == homoU && hosti == ADENOSTOMA_HOST ) {
                    // unstriped homozygote; bad on adenostoma
                    if ( randU() < S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE ) {
                        markedDead[i] = 1;
                    }
                }
                else if ( patternGT == het ) {
                    if ( hosti == ADENOSTOMA_HOST ) {
                        // unstriped, mostly dominant, is bad
                        if ( randU() < (STRIPE_DOMINANCE_COEFFICIENT * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) ) {
                            markedDead[i] = 1;
                        }
                    }
                    else {
                        if ( randU() < ((1.0 - STRIPE_DOMINANCE_COEFFICIENT) * S_COEFF_WRONG_STRIPE_IN_GREEN_NICHE) ) {
                            markedDead[i] = 1;
                        }
                    }
                }
                else {
                    // check
                    if ( !(patternGT == homoS && hosti == ADENOSTOMA_HOST) && !(patternGT == homoU && hosti == CEANOTHUS_HOST) ) {
                        // at least one of those should be true if we didn't get inside any of the other conditionals
                        fprintf(stderr, "\nError in viabilitySelection():\n\tconditionals not set up the way you think\n");
                        fprintf(stderr, "\tpatternGT = %i, hosti = %i, t = %li, i = %i\n", patternGT, hosti, t, i);
                        exit(-1);
                    }
                }
            }
        }
        
        gtpt += (DIPLOID * nLOCI);
    }
    
}

