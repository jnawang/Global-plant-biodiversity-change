#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <dirent.h>
#include <regex.h>
#include <string.h>
int main(int argc, char *argv[]) {
	//
	setbuf(stdout, NULL);
	if (argc < 2) {
        printf("Please provide the filenum as a command-line argument.\n");
        return 1;
    }
    const char* filenum = argv[1];
	int ifilenum = atoi(filenum) - 1;   //convert string to integer
	printf("%i\n", ifilenum);
//	fflush(stdout);
	const char* scenario;
	//
	if (ifilenum / 10 == 0) {
		scenario = "present.dat";
	} else if (ifilenum / 10 == 1) {
		scenario = "ESM_126_2081-2100.dat";
	} else if (ifilenum / 10 == 2) {
		scenario = "Realized_ESM_126_2081-2100.dat";
	} else if (ifilenum / 10 == 3) {
		scenario = "ESM_245_2081-2100.dat";
	} else if (ifilenum / 10 == 4) {
		scenario = "Realized_ESM_245_2081-2100.dat";
	} else if (ifilenum / 10 == 5) {
		scenario = "ESM_370_2081-2100.dat";
	} else if (ifilenum / 10 == 6) {
		scenario = "Realized_ESM_370_2081-2100.dat";
	} else if (ifilenum / 10 == 7) {
		scenario = "ESM_585_2081-2100.dat";
	} else if (ifilenum / 10 == 8) {
		scenario = "Realized_ESM_585_2081-2100.dat";
	} else if (ifilenum / 10 == 9) {
		scenario = "Realized025_ESM_126_2081-2100.dat";
	} else if (ifilenum / 10 == 10) {
		scenario = "Realized975_ESM_126_2081-2100.dat";
	} else if (ifilenum / 10 == 11) {
		scenario = "Realized025_ESM_245_2081-2100.dat";
	} else if (ifilenum / 10 == 12) {
		scenario = "Realized975_ESM_245_2081-2100.dat";
	} else if (ifilenum / 10 == 13) {
		scenario = "Realized025_ESM_370_2081-2100.dat";
	} else if (ifilenum / 10 == 14) {
		scenario = "Realized975_ESM_370_2081-2100.dat";
	} else if (ifilenum / 10 == 15) {
		scenario = "Realized025_ESM_585_2081-2100.dat";
	} else if (ifilenum / 10 == 16) {
		scenario = "Realized975_ESM_585_2081-2100.dat";       // why divided by 10? because 10 thread is on for one dat file. 
	}
	printf("%s\n", scenario);
    // work on some chunks
	int chunks[] = {1, 54, 97, 146, 175, 188, 198, 206, 215, 230, 276}; 
	int ichunk   = ifilenum % 10;
    //
	printf("%i%i\n", chunks[ichunk], chunks[ichunk+1]);
    printf("start\n");
	//	
	long nsps = 80000;   // note: nsps is not actual number of species, which is COLS!
	int* matrix_pair_current = (int*)malloc(nsps * nsps * sizeof(int));          // Values in the matrix_pair vary [0, 1000000], to reprent [0,1] 
    int* matrix_pair_future  = (int*)malloc(nsps * nsps * sizeof(int));
//  initialize current species pair matrix
	for (int i = 0; i < nsps; i++) {
		for (int j = 0; j < nsps; j++) {
			matrix_pair_current[i*nsps + j] = 0;
		}
	}
//  read current species pair matrix
	char  filename_tmp[100] = "/home/jnawang/SDM_process/novel_community/species_pairs_merge/present.dat";
	printf("%s\n", filename_tmp);
//
	FILE* binfile = fopen(filename_tmp, "rb");
//
	int ROWS, COLS;
	if (binfile != NULL) {
		fread(&COLS, sizeof(int), 1, binfile);
		printf("%i\n", COLS);
		// only output the upper triangle elements
		for (int i = 0; i < COLS; i++) {
			int num;
			fread(&num, sizeof(int), 1, binfile);
			int* colptrloc = (int*)malloc(num * sizeof(int));    
			int* colptrval = (int*)malloc(num * sizeof(int));
			fread(colptrloc, sizeof(int), num, binfile);
			fread(colptrval, sizeof(int), num, binfile);
			for (int j = 0; j < num; j++) {
				matrix_pair_current[i * nsps + colptrloc[j] - 1] += colptrval[j];   // I minus 1 because colptrloc is normal column number.
			}
			free(colptrloc);
			free(colptrval);
		}
		//				
	}
	else {
			printf("Cannot open file %s\n", filename_tmp);
	}
//
	fclose(binfile);

//
//  initialize future species pair matrix
	for (int i = 0; i < nsps; i++) {
		for (int j = 0; j < nsps; j++) {
			matrix_pair_future[i*nsps + j] = 0;
		}
	}
//  read future species pair matrix
	sprintf(filename_tmp, "%s%s", "/home/jnawang/SDM_process/novel_community/species_pairs_merge/", scenario);
	printf("%s\n", filename_tmp);
//
	binfile = fopen(filename_tmp, "rb");
//

	if (binfile != NULL) {
		fread(&COLS, sizeof(int), 1, binfile);
		printf("%i\n", COLS);
		// only output the upper triangle elements
		for (int i = 0; i < COLS; i++) {
			int num;
			fread(&num, sizeof(int), 1, binfile);
			int* colptrloc = (int*)malloc(num * sizeof(int));    
			int* colptrval = (int*)malloc(num * sizeof(int));
			fread(colptrloc, sizeof(int), num, binfile);
			fread(colptrval, sizeof(int), num, binfile);
			for (int j = 0; j < num; j++) {
				matrix_pair_future[i * nsps + colptrloc[j] - 1] += colptrval[j];   // I minus 1 because colptrloc is normal column number.
			}
			free(colptrloc);
			free(colptrval);
		}
		//				
	}
	else {
			printf("Cannot open file %s\n", filename_tmp);
	}
	fclose(binfile);

// compute number of species pairs whose overlap change from <=0.25 to >=0.75
    int total_species_pair_current = 0;
	int total_species_pair_future = 0;
	int total_species_pair_more_overlap = 0;
	for (int i = 0; i < nsps; i++) {
		for (int j = 0; j < nsps; j++) {
			if (matrix_pair_current[i*nsps + j] > 0) {
				total_species_pair_current += 1;
			}
			if (matrix_pair_future[i*nsps + j] > 0) {
				total_species_pair_future += 1;
			}			
			if (matrix_pair_current[i*nsps + j] <= 250000 &&
			    matrix_pair_current[i*nsps + j] >  0      &&
                matrix_pair_future[i*nsps + j]  >= 500000 + matrix_pair_current[i*nsps + j]) {
				total_species_pair_more_overlap += 1;
			}
		}
	}
    printf("Current species pair number: %d\n", total_species_pair_current);
	printf("Future species pair number: %d\n", total_species_pair_future);
	printf("Species pair more overlap <=0.25 to increase by 0.5: %d\n", total_species_pair_more_overlap);
	
//
//  read in global species distribution data
	int global_gridID[1000000];
	float global_novel[1000000];
	int ngrid = 0;
	const char* directoryPath = "/home/jnawang/SDM_process/novel_community/map2bin/";  // Replace with the path to your desired directory
    for (int ifile = chunks[ichunk]; ifile < chunks[ichunk + 1]; ifile++) {
//      start to work on one file
		clock_t start_time = clock();
		char  filename_tmp[100];
		sprintf(filename_tmp, "%s%d%s", directoryPath, ifile, scenario);
		printf("%s\n", filename_tmp);
//
		FILE* binfile = fopen(filename_tmp, "rb");
		if (binfile != NULL) {
			fread(&ROWS, sizeof(int), 1, binfile);
			fread(&COLS, sizeof(int), 1, binfile);
			printf("%d, %d\n", ROWS, COLS);
//		loop through each row; cannot miss the row with num = 0
			for (int n = 0; n < ROWS; n++) {
				int nnovel = 0;
				int gridid, num;
				fread(&gridid, sizeof(int), 1, binfile);
				fread(&num, sizeof(int), 1, binfile);
				//printf("%d, %d\n", gridid, num);
				global_gridID[ngrid] = gridid;
				if (num <= 0) {
					global_novel[ngrid] = 0.0;
				} else {
					int* colptr = (int*)malloc(num * sizeof(int));
					fread(colptr, sizeof(int), num, binfile);
					//
					if (num == 1) {
						global_novel[ngrid] = 0.0;
					} else {
						for (int i = 0; i < num; i++) {
							for (int j = i + 1; j < num; j++) {
								// compare with current species pair matrix
								if (matrix_pair_current[(colptr[i] - 1) * nsps + colptr[j] - 1] >  0  &&
								    matrix_pair_current[(colptr[i] - 1) * nsps + colptr[j] - 1] <= 250000 &&
								    matrix_pair_future[(colptr[i] - 1)  * nsps + colptr[j] - 1] >= 500000 + matrix_pair_current[(colptr[i] - 1) * nsps + colptr[j] - 1]) {
									nnovel ++;
								}
							}
						}					
						global_novel[ngrid] = (float)nnovel * 2.0 / num / (num - 1);
					}
					//
					free(colptr);
				}
				ngrid ++;
			}
//
		}
		else {
				printf("Cannot open file %s\n", filename_tmp);
		}
		fclose(binfile);
//
		double elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		printf("Done in %f seconds\n", elapsed_time);
    }
	
// output global grid ID and global novelty matrix
	char  filename[100];
	sprintf(filename, "%s%d%s", "/home/jnawang/SDM_process/novel_community/species_pair_overlap_chunk/", ichunk + 1, scenario);
	puts(filename);
	FILE* outfile = fopen(filename, "wb"); // Open the file in binary write mode
	if (outfile == NULL) {
		printf("Failed to open output file.\n");
		return 1;
	}
	// Write the array to the file
	fwrite(&ngrid, sizeof(int), 1, outfile);
	//
	fwrite(global_gridID, sizeof(int), ngrid, outfile);
	fwrite(global_novel,  sizeof(float), ngrid, outfile);
    //	
	fclose(outfile);
	printf("global species pair overlap chunk has been written to the file.\n");
    //
    free(matrix_pair_current);
    free(matrix_pair_future);	
	return 0;
//
}
