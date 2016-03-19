#pragma once

#include "types.h"

namespace GLnexus {

/// Determine whether the given record is a gVCF reference confidence record
/// (or else a "normal" record with at least one specific ALT allele)
bool is_gvcf_ref_record(const bcf1_t* record);

/// Detect an idiosyncratic class of records from certain HaplotypeCaller
/// versions which have QUAL == 0.0 and 0/0 genotype calls...we treat these as
/// "pseudo" reference confidence records.
bool is_pseudo_ref_record(const bcf_hdr_t* hdr, bcf1_t* record);

// Helper class for keeping track of the per-allele depth of coverage info in
// a bcf1_t record. There are a couple different cases to handle, depending on
// whether we're looking at a gVCF reference confidence record or a "regular"
// VCF record.
class AlleleDepthHelper {
    size_t n_sample_ = 0, n_allele_ = 0;
    bool is_g_ = false;
    int32_t *v_ = nullptr;
    int vlen_ = 0;

public:

    // The AlleleDepthHelper is constructed into an undefined state. Load()
    // must be invoked, successfully, before it can be used.
    AlleleDepthHelper()
        {}

    ~AlleleDepthHelper() {
        free(v_);
    }

    // The helper can be reused for multiple records by calling Load()
    // repeatedly. This will be slightly more efficient than using a new
    // helper for each record.
    Status Load(const std::string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record);

    // The behavior of remaining methods is undefined until Load() has been
    // invoked successfully

    bool is_gvcf_ref() { return is_g_; }

    // get depth for sample i & allele j
    unsigned depth(unsigned sample, unsigned allele) {
        if (sample >= n_sample_ || allele >= n_allele_) return 0;
        if (is_g_) {
            // the MIN_DP array has just one integer per sample
            if (allele != 0) return 0;
            return v_[sample];
        } else {
            // the AD array has one integer per allele per sample
            return v_[sample*n_allele_+allele];
        }
    }
};

} // namespace GLnexus
