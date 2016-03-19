#include "BCFHelpers.h"
#include <sstream>
#include <set>

using namespace std;

namespace GLnexus {

static const set<string> gvcf_nonref_symbols = { "<NON_REF>", "<*>" };
bool is_gvcf_ref_record(const bcf1_t* record) {
    return record->n_allele == 2 && gvcf_nonref_symbols.find(record->d.allele[1]) != gvcf_nonref_symbols.end();
}

bool is_pseudo_ref_record(const bcf_hdr_t* hdr, bcf1_t* record) {
    if (record->qual != 0.0) {
        return false;
    }
    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(hdr, record, &gt, &gtsz);
    for (int i = 0; i < record->n_sample; i++) {
        assert(i*2+1<nGT);
        if (bcf_gt_is_missing(gt[i*2]) || bcf_gt_allele(gt[i*2]) != 0 ||
            bcf_gt_is_missing(gt[i*2+1]) || bcf_gt_allele(gt[i*2+1]) != 0) {
            free(gt);
            return false;
        }
    }
    free(gt);
    return true;
}

Status AlleleDepthHelper::Load(const string& dataset, const bcf_hdr_t* dataset_header, bcf1_t* record) {
    n_sample_ = record->n_sample;
    n_allele_ = record->n_allele;
    // is this a gVCF reference confidence record?
    is_g_ = is_gvcf_ref_record(record);

    if (is_g_) {
        // if so, look for the MIN_DP FORMAT field (or equivalent)
        int nv = bcf_get_format_int32(dataset_header, record, "MIN_DP", &v_, &vlen_);

        if (nv != record->n_sample) {
            ostringstream errmsg;
            errmsg << dataset << " " << range(record).str();
            return Status::Invalid("genotyper: gVCF reference MIN_DP field is missing or malformed", errmsg.str());
        }
    } else {
        // this is a regular VCF record, so look for the AD FORMAT field (or equivalent)
        int nv = bcf_get_format_int32(dataset_header, record, "AD", &v_, &vlen_);

        if (nv == -1 || nv == -3) {
            // We allow the AD field to not exist mainly as a (poor)
            // workaround for some of our test case gVCFs not having it...

            if (nv == -3) {
                // AD is declared in the header, but not present in the
                // FORMAT fields for this record. We'll tolerate this for
                // an unusual observed class of variant records which have
                // INFO DP=0 (gVCF test case DP0_noAD.yml)

                nv = bcf_get_info_int32(dataset_header, record, "DP", &v_, &vlen_);
                if (nv != 1 || v_[0] != 0) {
                    ostringstream errmsg;
                    errmsg << dataset << " " << range(record).str();
                    return Status::Invalid("genotyper: VCF allele DP field is missing", errmsg.str());
                }
            }

            nv = record->n_sample * record->n_allele;
            if (vlen_ < nv) {
                v_ = (int32_t*) realloc(v_, nv*sizeof(int32_t));
                vlen_ = nv;
            }
            memset(v_, 0, nv*sizeof(int32_t));
        } else if (nv != record->n_sample * record->n_allele) {
            ostringstream errmsg;
            errmsg << dataset << " " << range(record).str();
            return Status::Invalid("genotyper: VCF DP field is malformed", errmsg.str());
        }
    }
    return Status::OK();
}


}