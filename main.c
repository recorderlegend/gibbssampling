
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define STRLEN 10
#define MOTIFLEN 7
#define NUMSTR 4
#define B 2
#define SC 0.5
#define SB 2

double randomgen(int x, int y) {
  // srand((unsigned)time(NULL) * getpid());

  float ran = (float)(rand());
  double val = x + (ran / RAND_MAX) * (y - x);

  return val;
}

int main() {
  char strings[NUMSTR][STRLEN] = {"ACCATGACAG", "GAGTATACCT", "CATGCTTACT",
                                  "CGGAATGCAT"};

  // array to store start index of random tuples
  int arridx[NUMSTR];
  // find a random index to select motif candidate from
  for (int i = 1; i < NUMSTR; i++) {
    int s = randomgen(0, STRLEN - MOTIFLEN);
    arridx[i] = s;
  }

  arridx[1] = 3;
  arridx[2] = 3;
  arridx[3] = 1;

  // alignment array
  char alignment[NUMSTR - 1][MOTIFLEN];
  // iterate for all of our strings except first one
  for (int i = 1; i < NUMSTR; i++) {

    // read motif candidate of size MOTIFLEN starting from position
    char tupmot[MOTIFLEN];
    for (int j = arridx[i]; j < arridx[i] + MOTIFLEN; j++) {
      tupmot[j - arridx[i]] = strings[i][j];
    }
    // set alignment array to the tuple
    for (int j = 0; j < MOTIFLEN; j++) {
      alignment[i - 1][j] = tupmot[j];
    }
  }
  for (int i = 0; i < NUMSTR - 1; i++) {
    for (int j = 0; j < MOTIFLEN; j++) {
      printf("%c", alignment[i][j]);
    }
    printf("\n");
  }

  // profile array
  // should only be 4 for this because only 4 nucleotides
  // 0:A 1:C 2:G 3:T
  int profile[4][MOTIFLEN + 1];
  memset(profile, 0, sizeof(profile));
  for (int i = 1; i < NUMSTR; i++) {
    // count the background frequency
    // copy string
    char currs[STRLEN];
    memcpy(currs, strings[i], STRLEN);
    for (int j = 0; j < STRLEN; j++) {
      if (j >= arridx[i] && j < arridx[i] + MOTIFLEN) {
        // printf("%c", currs[j]);
        if (currs[j] == 'A') {
          profile[0][j - arridx[i] + 1]++;
        }
        if (currs[j] == 'C') {
          profile[1][j - arridx[i] + 1]++;
        }
        if (currs[j] == 'G') {
          profile[2][j - arridx[i] + 1]++;
        }
        if (currs[j] == 'T') {
          profile[3][j - arridx[i] + 1]++;
        }
        continue;
      }
      // count up for background
      if (currs[j] == 'A') {
        profile[0][0]++;
      }
      if (currs[j] == 'C') {
        profile[1][0]++;
      }
      if (currs[j] == 'G') {
        profile[2][0]++;
      }
      if (currs[j] == 'T') {
        profile[3][0]++;
      }
    }
    // printf("\n");
  }

  // calculate probability profile array
  double profileprob[4][MOTIFLEN + 1];
  memset(profileprob, 0.0, sizeof(profileprob));

  // for background frequencies probablilty
  double bottommotif = ((NUMSTR - 1)) + SB;
  double bottombackground = ((NUMSTR - 1)*(STRLEN - MOTIFLEN)) + SB;
  //for background freq
  for (int i = 0; i < 4; i++) {
    double top = profile[i][0] + SC;
    profileprob[i][0] = top / bottombackground;
  }
  // for motif frequencies
  for (int i = 0; i < 4; i++) {
    for (int j = 1; j < MOTIFLEN + 1; j++) {
      double top = profile[i][j] + SC;
      profileprob[i][j] = top / bottommotif;
    }
  }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < MOTIFLEN + 1; j++) {
        printf("%f ", profileprob[i][j]);
        }
        printf("\n");
    }

  // calculate sampling step
  double samplingstep[STRLEN - MOTIFLEN + 1];
  memset(samplingstep, 0.0, sizeof(samplingstep));

  char nucleotides[4] = {'A', 'C', 'G', 'T'};
  // print profile array
  // for (int i = 0; i < 4; i++) {
  //   printf("%c: ", nucleotides[i]);
  //   for (int j = 0; j <= MOTIFLEN; j++) {
  //     printf("%f(%d) ", profileprob[i][j], profile[i][j]);
  //   }
  //   printf("\n");
  // }
}
