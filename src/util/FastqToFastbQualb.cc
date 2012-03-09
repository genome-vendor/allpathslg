///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FastqToFastbQualb.cc
 * \author tsharpe
 * \date Aug 12, 2009
 *
 * \brief Converts a Fastq file to a fastb and qualb.
 */
#include <cstdio>

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "system/ProcBuf.h"
#include <fstream>


static inline 
String Tag(String S = "FTFQ") { return Date() + " (" + S + "): "; } 


struct FastqEntry 
{
  char phred_min;
  char phred_max;
  BaseVec sequence;
  QualVec qualities;
};

class Swizzler
{
public:
    Swizzler( char offset ) : mOffset(offset), mMin('\x7f'), mMax('\x80') {}

    char operator()( char score )
    { if ( score < mMin ) mMin = score;
      if ( score > mMax ) mMax = score;
      return score > mOffset ? score-mOffset : 0; }

    char getMin() const { return mMin; }
    char getMax() const { return mMax; }

private:
    char mOffset;
    char mMin;
    char mMax;
};

bool NextFastqEntry( std::istream& is, FastqEntry &entry,
                            vec<size_t> * p_n_amb, const char phred_offset )
{
    bool result = false;
    std::string seqheader;
    std::string seq;
    std::string qualheader;
    std::string qual;
    if ( std::getline(is,seqheader) && seqheader.size() )
    {
        std::getline(is,seq);
        std::getline(is,qualheader);
        std::getline(is,qual);

        // ---- count number of ambiguous bases
        typedef std::string::iterator Itr;
        for ( Itr itr(seq.begin()), end(seq.end()); itr != end; ++itr )
            (*p_n_amb)[GeneralizedBase::fromChar(*itr).getAmbiguityCount()-1]++;

        entry.sequence = bvec(seq.begin(),seq.end(),GenCharToRandomBaseMapper());

        entry.qualities.clear().reserve(qual.size());
        Swizzler swiz(phred_offset);
        for ( Itr itr(qual.begin()), end(qual.end()); itr != end; ++itr )
            entry.qualities.push_back(swiz(*itr));
        entry.phred_min = swiz.getMin();
        entry.phred_max = swiz.getMax();
        result = true;
    }

    return result;
}









int main(int argc, char **argv) 
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(FASTQ);
  CommandArgument_String(OUT_HEAD);
  CommandArgument_UnsignedInt_OrDefault(READ_SIZE_THRESHOLD, 0);
  CommandArgument_String_OrDefault(FASTB_SUFFIX, ".fastb");
  CommandArgument_String_OrDefault(QUALB_SUFFIX, ".qualb");

  // Read trimming and reversing
  CommandArgument_Int_OrDefault(TRIM_START, 0);
  CommandArgument_Int_OrDefault_Doc(TRIM_END, 0, "If 0, don't trim");
  CommandArgument_Bool_OrDefault(REVERSE_READS, False);

  // Deal with +64 phred scores
  CommandArgument_Bool_OrDefault_Doc(PHRED_64, False, 
    "Phred score is encoded using +64, not +33 (default)");
  CommandArgument_Bool_OrDefault_Doc(FORCE_PHRED, False, 
    "True: accepts specified PHRED encoding. False: Aborts is incorrect PHRED detected.");
  

  CommandArgument_Bool_OrDefault(WRITE_FASTB, true);
  CommandArgument_Bool_OrDefault(WRITE_QUALB, true);
  EndCommandArguments;

  procbuf* pProcBuf = 0;
  istream* pIS = 0;
  if ( !FASTQ.EndsWith(".gz") )
      pIS = new std::ifstream(FASTQ.c_str());
  else
  {
      pProcBuf = new procbuf(("zcat "+FASTQ).c_str(),std::ios_base::in);
      pIS = new std::istream(pProcBuf);
  }

  // phred score can be encoded as score + 33 (default) or score + 64
  const int phred_offset = (PHRED_64 ? 64 : 33);

  BaseVecVec seqs;
  QualVecVec quals;

  seqs.reserve(10000);
  quals.reserve(10000);

  const bool trim_reads = (TRIM_START != 0 || TRIM_END != 0); 

  if (trim_reads)     cout << Tag() << "Trimming reads." << endl;
  if (REVERSE_READS)  cout << Tag() << "Reversing read pairs." << endl;
    
  FastqEntry e;
  e.phred_min = 127;
  e.phred_max = -127;
  BaseVec bv;
  QualVec qv;
  size_t n_reads_in = 0;
  size_t n_reads_out = 0;
  vec<size_t> n_amb(4,0);
  unsigned n_warns = 0;
  const unsigned n_warns_display = 10;
  while (NextFastqEntry(*pIS, e, &n_amb, phred_offset)) {
    if (e.sequence.size() > READ_SIZE_THRESHOLD) {

      const unsigned nb = e.sequence.size();
      const unsigned nq = e.qualities.size();
      const unsigned nb_missing = (nb < nq) ? nq - nb : 0;
      const unsigned nq_missing = (nb > nq) ? nb - nq : 0;
      
      if (nb != nq) {
        if (n_warns < n_warns_display)
          cout << Tag() << "WARNING: Read " << n_reads_in 
               << " has " << nb << " bases and " << nq << " quals." << endl;
        n_warns++;
      } 

      if (trim_reads) {
        const int end_base = (TRIM_END ? TRIM_END + 1 : e.sequence.size());
        bv.assign(e.sequence.begin() + TRIM_START, e.sequence.begin() + end_base);
        qv.assign(e.qualities.begin() + TRIM_START, e.qualities.begin() + end_base);
      }
      else {
        bv = e.sequence;
        qv = e.qualities;
      }
        
      if (REVERSE_READS) {
        bv.ReverseComplement();
        qv.ReverseMe();
      }

      seqs.push_back_reserve(bv);
      quals.push_back_reserve(qv);
      n_reads_out++;
    }
    n_reads_in++;
  }
  if (n_warns)
    cout << Tag() << "WARNING: There were " << n_warns 
         << " reads with inconsistent number of bases and quals." << endl;

  if (WRITE_FASTB) {  seqs.WriteAll(OUT_HEAD + FASTB_SUFFIX); }
  if (WRITE_QUALB) { quals.WriteAll(OUT_HEAD + QUALB_SUFFIX); }

  delete pIS;
  delete pProcBuf;

  cout << Tag() << setw(14) << n_reads_in << "  reads in fastq file." << endl;
  cout << Tag() << setw(14) << n_reads_out << "  reads in fastb/qualb files." << endl;

  {
    const double n_amb_tot = n_amb[0] + n_amb[1] + n_amb[2] + n_amb[3];
    const double frac_good = (double(n_amb[0]) + 
                              double(n_amb[1]) / 2.0 + 
                              double(n_amb[2]) / 3.0 +
                              double(n_amb[3]) / 4.0) / n_amb_tot;
    cout << Tag() << setw(14) << 100.0 * (1.0 - frac_good) 
         << "  % base errors due to ambiguities (all 'N' would yield 75%)." << endl;
  }
  
  cout << Tag() << setw(14) << (PHRED_64 ? "PHRED +64" : "PHRED +33") 
       << "  quality score encoding. Score set to Q0 if fastq score is negative." << endl;

  const int q_min = e.phred_min - phred_offset;  
  const int q_max = e.phred_max - phred_offset;

  cout << Tag() << setw(14) << q_min << "  lowest quality score." << endl;
  cout << Tag() << setw(14) << q_max << "  highest quality score." << endl;

  if (q_min < -5 || q_max > 61) {
    cout << Tag() << "!!!! QUALITY SCORES ARE NOT IN [-5:61] WITH SPECIFIED 'PHRED_64="
         << (PHRED_64 ? "True" : "False") << "'." 
         << endl;
    
    const int phred_offset_alt = (PHRED_64 ? 33 : 64);
    const int q_min_alt = e.phred_min - phred_offset_alt;  
    const int q_max_alt = e.phred_max - phred_offset_alt;
    cout << Tag() << "!!!! Set 'PHRED_64=" << (PHRED_64 ? "False" : "True") 
         << "' to get quality scores in [" << q_min_alt << "," << q_max_alt << "]";
    
    if ((q_min_alt >= -5 || q_max_alt <= 61)) cout << " which is correct." << endl;
    else                                      cout << " which is not correct either." << endl;


    if (!FORCE_PHRED) {
      cout << Tag() << "!!!! ABORTING because of irregular quality scores (use 'FORCE_PHRED=True' to override)." 
           << endl; 
      exit(1);
    }
  }
  cout << Tag() << "Done." << endl;

  return 0;
}
