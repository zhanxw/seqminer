#include "BGenVariant.h"

#include "FileIO.h"

void BGenVariant::makeTable(int ploidy, int allele) const {
  assert(ploidy >= 0 && allele >= 0);
  if ((size_t)ploidy > table.size()) {
    table.resize(ploidy);
  }
  for (int i = 0; i < ploidy; ++i) {
    const int s = table[i].size();
    if (allele > s) {
      table[i].resize(allele);
    }
    for (int j = s; j < allele; ++j) {
      if (i == 0) {
        table[i][j] = 1;
        continue;
      }

      if (j == 0) {
        table[i][j] = 1;
      } else if (j == 1) {
        table[i][j] = i + 1;
      } else {
        table[i][j] = table[i][j - 1] * (i + j) / j;
      }
    }
  }
}

/**
 * Look up the genotype by @param index
 *
 * Essentailly, look up the first column and get the second column (into
 0-based)

   index  genotype    allele counts
    0     111         (3,0,0)
    1     112         (2,1,0)
    2     122         (1,2,0)
    3     222         (0,3,0)
    4     113         (2,0,1)
    5     123         (1,1,1)
    6     223         (0,2,1)
    7     133         (1,0,2)
    8     233         (0,1,2)
    9     333         (0,0,3)
 *
 */

void BGenVariant::findGenotype(int idx, int ploidy, int allele,
                               std::vector<int>* geno) const {
  std::vector<int>& g = *geno;
  g.resize(ploidy);

  makeTable(ploidy, allele);
  --ploidy;
  ++idx;
  while (ploidy >= 0) {
    int i = 0;
    int cumsum = 0;
    for (; i < allele; ++i) {
      if (cumsum + table[ploidy][i] >= idx) {
        idx -= cumsum;
        break;
      } else {
        cumsum += table[ploidy][i];
      }
    }

    g[ploidy] = i;
    ploidy--;
  }
}

void BGenVariant::printGT(int i, FileWriter* fp) const {
  if (isPhased) {
    if (missing[i]) {
      printGTMissingFromHaplotype(fp);
    } else {
      printGTFromHaplotype(i, fp);
    }
  } else {
    if (missing[i]) {
      printGTMissingFromGenotype(fp);
    } else {
      switch (K) {
        case 2:
          printGTAllele2FromGenotype(i, fp);
          break;
        case 1:
          printGTAllele1FromGenotype(i, fp);
          break;
        default:
          printGTAlleleGeneralFromGenotype(i, fp);
          break;
      }
    }
  }
}
void BGenVariant::printGTMissingFromHaplotype(FileWriter* fp) const {
  fp->write(".");
  for (int i = 1; i < ploidy[i]; ++i) {
    fp->write("|.");
  }
}
void BGenVariant::printGTMissingFromGenotype(FileWriter* fp) const {
  fp->write(".");
  for (int i = 1; i < ploidy[i]; ++i) {
    fp->write("/.");
  }
}
void BGenVariant::printGTAllele1FromGenotype(int i, FileWriter* fp) const {
  fp->write("0");
  for (int i = 1; i < ploidy[i]; ++i) {
    fp->write("/0");
  }
}
void BGenVariant::printGTAllele2FromGenotype(int i, FileWriter* fp) const {
  if (ploidy[i] == 2) {  // prob[index[i]] stores p(00), p(01), p(11)
    const float prob0 = prob[index[i]];
    const float prob1 = prob[index[i] + 1];
    const float prob2 = prob[index[i] + 2];

    if (prob0 > prob1 && prob0 > prob2) {
      fp->write("0/0");
    } else if (prob1 > prob0 && prob1 > prob2) {
      fp->write("0/1");
    } else {
      fp->write("1/1");
    }
  } else if (ploidy[i] == 1) {  // prob[index[i]] stores p(0), p(1)
    const float prob0 = prob[index[i]];
    const float prob1 = prob[index[i] + 1];

    if (prob0 > prob1) {
      fp->write("0");
    } else {
      fp->write("1");
    }
  } else {  // prob[index[i]] stores p(00..0), p(00...1), ... p(22..2)
    printGTAlleleGeneralFromGenotype(i, fp);
  }
}
void BGenVariant::printGTAlleleGeneralFromGenotype(int idx,
                                                   FileWriter* fp) const {
  // let Z = ploidity, K = alleles
  // total choose(N+K-1, K-1) genotypes
  int maxIdx = index[idx];
  float maxVal = prob[maxIdx];
  for (int i = maxIdx + 1; i < index[idx + 1]; ++i) {
    if (prob[i] > maxVal) {
      maxIdx = i;
      maxVal = prob[i];
    }
  }
  std::vector<int> geno;
  findGenotype(maxIdx - index[idx], ploidy[idx], K, &geno);
  for (size_t i = 0; i < geno.size(); ++i) {
    if (i) {
      fp->write("/");
    }
    fp->printf("%d", geno[i]);
  }
}
void BGenVariant::printGTFromHaplotype(int ii, FileWriter* fp) const {
  const int Z = ploidy[ii];
  //  K alleles
  int idx = index[ii];
  int maxIdx;
  float maxValue;
  for (int i = 0; i < Z; ++i) {
    maxIdx = 0;
    maxValue = prob[idx];
    idx++;
    for (int j = 1; j < K; ++j) {
      if (prob[idx] > maxValue) {
        maxValue = prob[idx];
        maxIdx = j;
      }
      idx++;
    }
    if (i) {
      fp->write("|");
    }
    fp->printf("%d", maxIdx);
  }
  assert(idx == index[ii + 1]);
}

/// Handle GP  //////////////////////////////////////////////////
/// output genotype probability
void BGenVariant::printGP(int i, FileWriter* fp) const {
  if (missing[i]) {
    printGPMissing(i, fp);
    return;
  }

  switch (K) {
    case 2:
      printGPAllele2(i, fp);
      break;
    case 1:
      printGPAllele1(i, fp);
      break;
    default:
      printGPAlleleGeneral(i, fp);
      break;
  }
}
void BGenVariant::printGPMissing(int idx, FileWriter* fp) const {
  for (int i = index[idx]; i < index[idx + 1]; ++i) {
    if (i != index[idx]) {
      fp->write(",");
    }
    fp->write(".");
  }
}
void BGenVariant::printGPAllele1(int i, FileWriter* fp) const {
  fp->write("1");
}
void BGenVariant::printGPAllele2(int i, FileWriter* fp) const {
  if (ploidy[i] == 2) {  // prob of 00, 01, 11
    const float prob0 = prob[index[i]];
    const float prob1 = prob[index[i] + 1];
    const float prob2 = prob[index[i] + 2];

    fp->printf("%g,%g,%g", prob0, prob1, prob2);
  } else if (ploidy[i] == 1) {  // prob of 0,1
    const float prob0 = prob[index[i]];
    const float prob1 = prob[index[i] + 1];

    fp->printf("%g,%g", prob0, prob1);
  } else {
    // let Z = ploidity, K = alleles
    // total choose(N+K-1, K-1) genotypes
    printGPAlleleGeneral(i, fp);
  }
}
void BGenVariant::printGPAlleleGeneral(int idx, FileWriter* fp) const {
  for (int i = index[idx]; i < index[idx + 1]; ++i) {
    if (i != index[idx]) {
      fp->write(",");
    }
    fp->printf("%g", prob[i]);
  }
}
/// Handle HP  //////////////////////////////////////////////////
// handle haplotype probability
void BGenVariant::printHP(int i, FileWriter* fp) const {
  if (missing[i]) {
    printHPMissing(i, fp);
    return;
  }
  printHPAlleleGeneral(i, fp);
}
void BGenVariant::printHPMissing(int idx, FileWriter* fp) const {
  for (int i = index[idx]; i < index[idx + 1]; ++i) {
    if (i != index[idx]) {
      fp->write(",");
    }
    fp->write(".");
  }
}
void BGenVariant::printHPAlleleGeneral(int idx, FileWriter* fp) const {
  for (int i = index[idx]; i < index[idx + 1]; ++i) {
    if (i != index[idx]) {
      fp->write(",");
    }
    fp->printf("%g", prob[i]);
  }
}

/// Handle dosage //////////////////////////////////////////////////
void BGenVariant::printDosage(int i, FileWriter* fp) const {
  const float dosage = computeDosage(i);
  if (dosage < 0) {
    fp->write(".");
  } else {
    fp->printf("%g", dosage);
  }
}

float BGenVariant::computeDosage(int i) const {
  if (missing[i]) {
    return -9;
  }
  if (ploidy[i] == 2 && K == 2) {
    if (isPhased) {
      const float prob1 = prob[index[i] + 1];
      const float prob2 = prob[index[i] + 3];

      return (prob1 + prob2);
    } else {
      // const float prob0 = prob[index[i]];
      const float prob1 = prob[index[i] + 1];
      const float prob2 = prob[index[i] + 2];

      return (prob1 + 2.0 * prob2);
    }
  } else {
    return -9;
  }
}
