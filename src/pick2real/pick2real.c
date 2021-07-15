#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#define maxlines 1000000
int main(int argc, char **argv) {
  int i, j;
  FILE *fp1, *fp3, *fp2;
  char inputfile_p[256], inputfile_s[256], output[256];
  char **inf_p, **inf_s, *inf0, *date, *station, *net;
  char *tt,*weight,*amplitude;
  int error = 0;
  char flag;
  char *year, *day, *month;
  char ymd_str[256];
  char *inf_p0, *inf_s0;
  for (i = 1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'P':
        sscanf(&argv[i][2], "%s", inputfile_p);
        break;
      case 'S':
        sscanf(&argv[i][2], "%s", inputfile_s);
        break;
      default:
        error = 1;
        break;
      }
    }
  }
  inf_p = (char **)malloc(sizeof(char *) * maxlines);
  inf_s = (char **)malloc(sizeof(char *) * maxlines);
  for (i = 0; i < maxlines; i++) {
    inf_p[i] = (char *)malloc(sizeof(char) * 256);
    inf_s[i] = (char *)malloc(sizeof(char) * 256);
  }
  if ((fp1 = fopen(inputfile_p, "r")) == NULL) {
    fprintf(stderr, "Unable to open file %s\n", inputfile_p);
    exit(-1);
  }
  if ((fp3 = fopen(inputfile_s, "r")) == NULL) {
    fprintf(stderr, "Unable to open file %s\n", inputfile_s);
    exit(-1);
  }
  for (i = 0; i < maxlines; i++) {
    fscanf(fp1, "%s", inf_p[i]);
    fscanf(fp3, "%s", inf_s[i]);
  }
  fclose(fp1);
  fclose(fp3);
  for (i = 0; i < maxlines; i++) {
    if (inf_p[i][0] == NULL)
      break;
    inf_p0 = inf_p[i];
    year = strtok(inf_p0, ",");
    month = strtok(NULL, ",");
    day = strtok(NULL, ",");
    net = strtok(NULL, ",");
    station = strtok(NULL, ",");
    flag = strtok(NULL, ",");
    tt = strtok(NULL, ",");
    weight = strtok(NULL, ",");
    amplitude = strtok(NULL, ",");
    sprintf(output, "%s%s%s/%s.%s.P.txt", year, month, day, net, station);
    sprintf(ymd_str, "%s%s%s", year, month, day);
    if (access(ymd_str, 0) == -1)
      mkdir(ymd_str, 0777);
    if ((fp2 = fopen(output, "at")) == NULL) {
      fprintf(stderr, "Unable to open file %s\n", output);
      exit(-1);
    }
    fprintf(fp2, "%s\t%s\t%s\n", tt, weight,amplitude);
    fclose(fp2);
  }
  for (i = 0; i < maxlines; i++) {
    if (inf_s[i][0] == NULL)
      break;
    inf_s0 = inf_s[i];
    year = strtok(inf_s0, ",");
    month = strtok(NULL, ",");
    day = strtok(NULL, ",");
    net = strtok(NULL, ",");
    station = strtok(NULL, ",");
    flag = strtok(NULL, ",");
    tt = strtok(NULL, ",");
    weight = strtok(NULL, ",");
    amplitude = strtok(NULL, ",");
    sprintf(output, "%s%s%s/%s.%s.S.txt", year, month, day, net, station);
    sprintf(ymd_str, "%s%s%s", year, month, day);
    if (access(ymd_str, 0) == -1)
      mkdir(ymd_str, 0777);
    if ((fp2 = fopen(output, "at")) == NULL) {
      fprintf(stderr, "Unable to open file %s\n", output);
      exit(-1);
    }
    fprintf(fp2, "%s\t%s\t%s\n", tt, weight,amplitude);
    fclose(fp2);
  }
  // free mem.
  for (i = 0; i < maxlines; i++) {
    free(inf_p[i]);
    free(inf_s[i]);
  }
  free(inf_p);
  free(inf_s);
  return 0;
}
