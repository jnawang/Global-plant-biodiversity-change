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
		scenario = "Realized975_ESM_585_2081-2100.dat";    
	}
	printf("%s\n", scenario);
	char pattern0[50] = "[0-9]" ;
	strcat(pattern0, scenario);	
	//
	int global_gridID[2000000];
	float global_novel[2000000];
	int ngrid = 0;
	// read files with a pattern
	const char* directoryPath = "/home/jnawang/SDM_process/novel_community/novel_species_pairs_chunk/";  // Replace with the path to your desired directory  C:/Users/Junna/Box/Biodiversity/SDM/bin/
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
				int num;
				fread(&num, sizeof(int), 1, binfile);
				printf("%i\n", num);
				int* colptrint     = (int*)malloc(num * sizeof(int));
				float* colptrfloat = (float*)malloc(num * sizeof(float));
				fread(colptrint, sizeof(int), num, binfile);
				fread(colptrfloat, sizeof(float), num, binfile);
				for (int i = 0; i < num; i++) {
					global_gridID[ngrid + i] = colptrint[i];
					global_novel[ngrid + i]  = colptrfloat[i];
				}
				// print out the first cell's result to check
				printf("%i %f\n", global_gridID[ngrid], global_novel[ngrid]);
				free(colptrint);
				free(colptrfloat);
				ngrid += num;
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
	// write out all the gridid and novelty
	char  filename[100];
	sprintf(filename, "%s%s", "/home/jnawang/SDM_process/novel_community/novel_species_pairs_merge/", scenario);
	puts(filename);
	FILE* outfile = fopen(filename, "wb"); // Open the file in binary write mode
	if (outfile == NULL) {
		printf("Failed to open output file.\n");
		return 1;
	}
	//
	fwrite(&ngrid, sizeof(int), 1, outfile);
	// need to ensure this works
	fwrite(global_gridID, sizeof(int), ngrid, outfile);
	fwrite(global_novel,  sizeof(float), ngrid, outfile);
    //	
	fclose(outfile);
	printf("global novelty chunk has been written to the file.\n");
    //
	return 0;
//
}
