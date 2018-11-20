#ifndef _BGENVARIANT_H_
#define _BGENVARIANT_H_

#include <stdint.h>  // uint32_t
#include <cassert>
#include <string>
#include <vector>

class FileWriter;

static std::vector<std::vector<int> >
    table;  // for quickly lookup genotypes in colex order

struct BGenVariant {
  BGenVariant(const uint32_t& n) : N(n){};
  const uint32_t& N;
  // variant id data
  std::string varid;
  std::string rsid;
  std::string chrom;
  uint32_t pos;
  uint16_t K;                        // max = 2^16 = 65536
  std::vector<std::string> alleles;  // K alleles

  // store parsed results
  uint8_t B;  // bits to represent probability
  std::vector<bool> missing;
  std::vector<uint8_t> ploidy;  // Z
  bool isPhased;
  // prob[ index[0] .. index[1]] is the probability for the first individual
  // prob [ index[N-1] .. index[N]] is the probability for the last individual
  // each individual has Z * K probabilites when phased, or
  // ( Z+K-1 ) choose ( K-1 ) probabilities when unphased
  // here: Z is the ploidy and K is the number of alleles
  std::vector<int> index;
  std::vector<float> prob;  // probability array

  /// Handle GT (genotye) //////////////////////////////////////////////////
  void printGT(int i, FileWriter* fp) const;
  // print missing genotypes when variant is phased
  void printGTMissingFromHaplotype(FileWriter* fp) const;
  // print missing genotypes when variant is unphased  
  void printGTMissingFromGenotype(FileWriter* fp) const;
  // print genotype when num of alleles = 1, and variant is unphased
  void printGTAllele1FromGenotype(int i, FileWriter* fp) const;
  // print genotype when num of alleles = 2, and variant is unphased  
  void printGTAllele2FromGenotype(int i, FileWriter* fp) const;
  // print genotype for any num of alleles, and variant is unphased    
  void printGTAlleleGeneralFromGenotype(int idx, FileWriter* fp) const;
  // print genotype for any num of alleles, and variant is phased
  void printGTFromHaplotype(int i, FileWriter* fp) const;

  /// Handle GP  //////////////////////////////////////////////////
  /// output genotype probability
  void printGP(int i, FileWriter* fp) const;
  void printGPMissing(int i, FileWriter* fp) const;
  void printGPAllele1(int i, FileWriter* fp) const;
  void printGPAllele2(int i, FileWriter* fp) const;
  void printGPAlleleGeneral(int i, FileWriter* fp) const;
  /// Handle HP  //////////////////////////////////////////////////
  // handle haplotype probability
  void printHP(int i, FileWriter* fp) const;
  void printHPMissing(int i, FileWriter* fp) const;
  void printHPAlleleGeneral(int i, FileWriter* fp) const;

  /// Handle dosage //////////////////////////////////////////////////
  void printDosage(int i, FileWriter* fp) const;
  float computeDosage(int i) const;
  /// auxillary functions
 private:
  void makeTable(int ploidy, int allele) const;
  void findGenotype(int idx, int ploidy, int allele,
                    std::vector<int>* geno) const;
  
};

#endif /* _BGENVARIANT_H_ */
