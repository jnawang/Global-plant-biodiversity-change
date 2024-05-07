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
	int* matrix_pair = (int*)malloc(nsps * nsps * sizeof(int));
	printf("1\n");

	for (int i = 0; i < nsps; i++) {
		for (int j = 0; j < nsps; j++) {
			matrix_pair[i*nsps + j] = 0;
		}
	}
	printf("2\n");
//
    const char* directoryPath = "/home/jnawang/SDM_process/novel_community/map2bin/";  // Replace with the path to your desired directory  C:/Users/Junna/Box/Biodiversity/SDM/bin/
	int ROWS, COLS;
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
//
				for (int n = 0; n < ROWS; n++) {
						int gridid, num;
						fread(&gridid, sizeof(int), 1, binfile);
						fread(&num, sizeof(int), 1, binfile);
						if (num == 0)
								continue;
						int* colptr = (int*)malloc(num * sizeof(int));
						fread(colptr, sizeof(int), num, binfile);
						for (int i = 0; i < num; i++) {
								for (int j = i + 1; j < num; j++) {
										matrix_pair[(colptr[i] - 1) * nsps + colptr[j] - 1] ++;
								}
						}
						free(colptr);
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
// I need to output large matrix before free it. 
	char  filename[100];
	sprintf(filename, "%s%d%s", "/home/jnawang/SDM_process/novel_community/species_pairs_chunk/", ichunk + 1, scenario);
	puts(filename);
	FILE* outfile = fopen(filename, "wb"); // Open the file in binary write mode
	if (outfile == NULL) {
		printf("Failed to open output file.\n");
		return 1;
	}
	// Write the array to the file
	// write number of species: COLS
	fwrite(&COLS, sizeof(int), 1, outfile);
	//
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
}
