#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <dirent.h>
#include <regex.h>
#include <string.h>
int main(int argc, char *argv[]) {
	if (argc < 2) {
        printf("Please provide the filenum as a command-line argument.\n");
        return 1;
    }
    const char* filenum = argv[1];
	int ifilenum = atoi(filenum) - 1;   //convert string to integer
	printf("%i\n", ifilenum);
	
	//
	const char* scenario;
	if (ifilenum == 0) {
		scenario = "present.dat";
	} else if (ifilenum == 1) {
		scenario = "ESM_126_2081-2100.dat";
	} else if (ifilenum == 2) {
		scenario = "Realized_ESM_126_2081-2100.dat";
	} else if (ifilenum == 3) {
		scenario = "ESM_245_2081-2100.dat";
	} else if (ifilenum == 4) {
		scenario = "Realized_ESM_245_2081-2100.dat";
	} else if (ifilenum == 5) {
		scenario = "ESM_370_2081-2100.dat";
	} else if (ifilenum == 6) {
		scenario = "Realized_ESM_370_2081-2100.dat";
	} else if (ifilenum == 7) {
		scenario = "ESM_585_2081-2100.dat";
	} else if (ifilenum == 8) {
		scenario = "Realized_ESM_585_2081-2100.dat";
	} else if (ifilenum == 9) {
		scenario = "Realized025_ESM_126_2081-2100.dat";
	} else if (ifilenum == 10) {
		scenario = "Realized975_ESM_126_2081-2100.dat";
	} else if (ifilenum == 11) {
		scenario = "Realized025_ESM_245_2081-2100.dat";
	} else if (ifilenum == 12) {
		scenario = "Realized975_ESM_245_2081-2100.dat";
	} else if (ifilenum == 13) {
		scenario = "Realized025_ESM_370_2081-2100.dat";
	} else if (ifilenum == 14) {
		scenario = "Realized975_ESM_370_2081-2100.dat";
	} else if (ifilenum == 15) {
		scenario = "Realized025_ESM_585_2081-2100.dat";
	} else if (ifilenum == 16) {
		scenario = "Realized975_ESM_585_2081-2100.dat";       // why divided by 10? because 10 thread is on for one dat file. 
	}
	printf("%s\n", scenario);
	char pattern0[50] = "[0-9]" ;
	strcat(pattern0, scenario);
//
	long nsps = 80000;   // need to update based on file data
	int* matrix_pair = (int*)malloc(nsps * nsps * sizeof(int));

	for (int i = 0; i < nsps; i++) {
		for (int j = 0; j < nsps; j++) {
			matrix_pair[i*nsps + j] = 0;
		}
	}
	puts("Finish allocate species pair matrix");
//
	// read species grid; this should have three columns;
	// define an array to put number of grids for each species
	int   species_grid[nsps];
	char  filename_tmp[100];
	sprintf(filename_tmp, "%s%d%s", "/home/jnawang/SDM_process/novel_community/species_pairs_merge/", ifilenum+1, "species_grid.bin");
//
	FILE* binfile = fopen(filename_tmp, "rb");
//
	if (binfile != NULL) {
		int num;
		fread(&num, sizeof(int), 1, binfile);
		printf("%i\n", num);
		//
		int* colptr = (int*)malloc(num * sizeof(int));
		fread(colptr, sizeof(int), num, binfile);
		for (int i = 0; i<num; i++){
			species_grid[i] = colptr[i];
		}
		for (int i = 0; i<10; i++){
			printf("%i\n", species_grid[i]);
		}
		free(colptr);
	}
	else {
		printf("Cannot open file %s\n", filename_tmp);
	}
//
	fclose(binfile);

// read files with a certain pattern
    const char* directoryPath = "/home/jnawang/SDM_process/novel_community/species_pairs_chunk/";  // Replace with the path to your desired directory
    const char* pattern = pattern0;       // Example: Match files ending with this
	DIR* directory;
	struct dirent* file;
	regex_t regex;
	//
    int regexCompilationResult = regcomp(&regex, pattern, REG_EXTENDED | REG_NOSUB);
    if (regexCompilationResult != 0) {
        printf("Failed to compile regular expression.\n");
        return 0;
    }
//
    directory = opendir(directoryPath);
    if (directory == NULL) {
        printf("Failed to open directory.\n");
        return 0;
    }
//
	int COLS;
    while ((file = readdir(directory)) != NULL) {
        if (regexec(&regex, file->d_name, 0, NULL, 0) == 0) {
            printf("%s\n", file->d_name);
//          start to work on one file
			clock_t start_time = clock();
			char  filename_tmp[100];
			strcpy(filename_tmp, directoryPath);
			strcat(filename_tmp, file->d_name);
			printf("%s\n", filename_tmp);
//
			FILE* binfile = fopen(filename_tmp, "rb");
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
						matrix_pair[i * nsps + colptrloc[j] - 1] += colptrval[j];   // I minus 1 because colptrloc is normal column number.
					}
					free(colptrloc);
					free(colptrval);
				}
				//				
			}
			else {
					printf("Cannot open file %s\n", file->d_name);
			}
//
			fclose(binfile);
//
			double elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
			printf("Done in %f seconds\n", elapsed_time);		
        }
    }
    closedir(directory);
    regfree(&regex);
// need to test species pair matrix before standarizing it.	I did the test!
	
// Standarize the species matrix based on species grid
	for (int i = 0; i < COLS; i++) {
		for (int j = i + 1; j < COLS; j++) {
			if (matrix_pair[i * nsps + j] != 0) {
				int min = (species_grid[i] < species_grid[j]) ? species_grid[i] : species_grid[j];
				matrix_pair[i * nsps + j] = matrix_pair[i * nsps + j] * (long)1000000 / min;    // multiply the result by 1 million
				// print some results to validate the calculation
				if (i == 0 || i == 1000 || i == COLS - 1) {
					printf("%i %i %i\n", i, j, matrix_pair[i * nsps + j]);
				}
			}
		}
	}
//  write down the species pair matrix
// I need to output large matrix before free it. 
	char  filename[100];
	sprintf(filename, "%s%s", "/home/jnawang/SDM_process/novel_community/species_pairs_merge/", scenario);
	puts(filename);
	FILE* outfile = fopen(filename, "wb"); // Open the file in binary write mode
	if (outfile == NULL) {
		printf("Failed to open output file.\n");
		return 1;
	}
//
	// Write the array to the file
	// write number of species: COLS
	fwrite(&COLS, sizeof(int), 1, outfile);
	int* colptrloc = (int*)malloc(COLS * sizeof(int));
	int* colptrval = (int*)malloc(COLS * sizeof(int));    
	// only output the upper triangle elements
	for (int i = 0; i < COLS; i++) {
		int num = 0;
		for (int j = i + 1; j < COLS; j++) {
			if (matrix_pair[i * nsps + j] != 0) {
				colptrloc[num] = j + 1;                         // this is the actual col to be consistent with R
				colptrval[num] = matrix_pair[i * nsps + j];     
				num++;
				// print some results to validate the model
				if (i == 0 || i == 1000 || i == COLS - 1) {
					printf("%i %i %i %i\n", i, j, num, matrix_pair[i * nsps + j]);
				}
			}
		}
		fwrite(&num, sizeof(int), 1, outfile);
		fwrite(colptrloc, sizeof(int), num, outfile);
		fwrite(colptrval, sizeof(int), num, outfile);
	}
	free(colptrloc);
	free(colptrval);
//
	fclose(outfile);
	printf("Species pair array has been written to the file.\n");
//
    free(matrix_pair);
	return 0;
//
}