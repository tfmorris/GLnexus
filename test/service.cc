#include <iostream>
#include <algorithm>
#include <vcf.h>
#include "service.h"
#include "unifier.h"
#include "genotyper.h"
#include "utils.cc"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

// Wrap BCFData to simulate I/O errors for quasi-random testing of errors
class SimFailBCFData : public BCFData {
    BCFData& inner_;
    const size_t fail_every_;
    mutable size_t i_ = 0;
    mutable bool failed_once_ = false;

    SimFailBCFData(BCFData& inner, size_t fail_every)
        : inner_(inner), fail_every_(fail_every) {}

public:
    static Status Open(BCFData& inner, size_t fail_every, unique_ptr<SimFailBCFData>& ans) {
        ans.reset(new SimFailBCFData(inner, fail_every));
        return Status::OK();
    }

    Status dataset_header(const string& dataset, shared_ptr<const bcf_hdr_t>& hdr) const override {
        if (i_++ % fail_every_ == 0) {
            failed_once_ = true;
            return Status::IOError("SIM");
        }
        return inner_.dataset_header(dataset, hdr);
    }

    Status dataset_range(const string& dataset, const bcf_hdr_t *hdr, const range& pos,
                         bcf_predicate predicate, vector<shared_ptr<bcf1_t> >& records) override {
        if (i_++ % fail_every_ == 0) {
            failed_once_ = true;
            return Status::IOError("SIM");
        }
        return inner_.dataset_range(dataset, hdr, pos, predicate, records);
    }

    bool failed_once() { return failed_once_; }
};

TEST_CASE("service::discover_alleles") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    SECTION("multiple ranges") {
        vector<range> ranges;
        ranges.push_back(range(0, 1000, 1001));
        ranges.push_back(range(0, 1001, 1002));
        ranges.push_back(range(0, 1010, 1013));
        ranges.push_back(range(1, 1000, 1001));
        ranges.push_back(range(1, 1010, 1012));
        ranges.push_back(range(1, 2000, 2100));
        vector<discovered_alleles> mals;
        s = svc->discover_alleles("<ALL>", ranges, mals);
        REQUIRE(s.ok());
        REQUIRE(mals.size() == ranges.size());

        REQUIRE(mals[0].size() == 2);
        REQUIRE(mals[0].find(allele(range(0, 1000, 1001), "A"))->second.copy_number == 6);
        REQUIRE(mals[0].find(allele(range(0, 1000, 1001), "G"))->second.copy_number == 6);

        REQUIRE(mals[1].size() == 4);
        REQUIRE(mals[1].find(allele(range(0, 1001, 1002), "A"))->second.copy_number == 6);
        REQUIRE(mals[1].find(allele(range(0, 1001, 1002), "C"))->second.copy_number == 2);
        REQUIRE(mals[1].find(allele(range(0, 1001, 1002), "G"))->second.copy_number == 2);
        REQUIRE(mals[1].find(allele(range(0, 1001, 1002), "T"))->second.copy_number == 2);

        REQUIRE(mals[2].size() == 4);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1012), "AG"))->second.copy_number == 3);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1012), "CC"))->second.copy_number == 3);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1013), "AGA"))->second.copy_number == 2);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1013), "CCC"))->second.copy_number == 4);

        REQUIRE(mals[3].size() == 2);
        REQUIRE(mals[3].find(allele(range(1, 1000, 1001), "A"))->second.copy_number == 2);
        REQUIRE(mals[3].find(allele(range(1, 1000, 1001), "AA"))->second.copy_number == 4);

        REQUIRE(mals[4].size() == 2);
        REQUIRE(mals[4].find(allele(range(1, 1010, 1012), "AG"))->second.copy_number == 3);
        REQUIRE(mals[4].find(allele(range(1, 1010, 1012), "CC"))->second.copy_number == 3);

        REQUIRE(mals[5].empty());
    }

    SECTION("simulate I/O errors - single") {
        unique_ptr<SimFailBCFData> faildata;
        bool worked = false;

        for (size_t fail_every = 1; fail_every < 100; fail_every++) {
            s = SimFailBCFData::Open(*data, fail_every, faildata);
            REQUIRE(s.ok());

            s = Service::Start(service_config(), *data, *faildata, svc);
            REQUIRE(s.ok());

            s = svc->discover_alleles("<ALL>", range(0, 0, 1099), als);
            if (faildata->failed_once()) {
                worked = true;
                REQUIRE(s == StatusCode::IO_ERROR);
                REQUIRE(s.str() == "IOError: SIM");
            } else {
                REQUIRE(s.ok());
            }
        }

        REQUIRE(worked);
    }

    SECTION("simulate I/O errors - multi") {
        unique_ptr<SimFailBCFData> faildata;
        bool worked = false;

        for (size_t fail_every = 1; fail_every < 100; fail_every++) {
            s = SimFailBCFData::Open(*data, fail_every, faildata);
            REQUIRE(s.ok());

            s = Service::Start(service_config(), *data, *faildata, svc);
            REQUIRE(s.ok());

            vector<range> ranges;
            ranges.push_back(range(0, 1000, 1001));
            ranges.push_back(range(0, 1001, 1002));
            ranges.push_back(range(0, 1010, 1013));
            ranges.push_back(range(1, 1000, 1001));
            ranges.push_back(range(1, 1010, 1012));
            vector<discovered_alleles> mals;
            s = svc->discover_alleles("<ALL>", ranges, mals);
            if (faildata->failed_once()) {
                worked = true;
                REQUIRE(s == StatusCode::IO_ERROR);
                REQUIRE(s.str() == "IOError: SIM");
            } else {
                REQUIRE(s.ok());
            }
        }

        REQUIRE(worked);
    }
}

TEST_CASE("service::discover_alleles gVCF") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"NA12878D_HiSeqX.21.10009462-10009469.gvcf", "NA12878D_HiSeqX.21.10009462-10009469.bogus.gvcf"}, data);
    cout << s.str() << endl;
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    SECTION("exclude symbolic alleles") {
        s = svc->discover_alleles("<ALL>", range(0, 10000000, 10010000), als);

        REQUIRE(als.size() == 2);
        auto p = als.find(allele(range(0, 10009463, 10009465), "TA"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.copy_number == 1.7152f);
        p = als.find(allele(range(0, 10009463, 10009465), "T"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.copy_number == 0.2847f);
    }

    SECTION("exclusion/detection of bogus alleles") {
        s = svc->discover_alleles("<ALL>", range(1, 10009463, 10009465), als);

        REQUIRE(als.size() == 2);
        auto p = als.find(allele(range(1, 10009463, 10009465), "TA"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.copy_number == 1.7152f);
        p = als.find(allele(range(1, 10009463, 10009465), "T"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.copy_number == 0.2847f);

        s = svc->discover_alleles("<ALL>", range(1, 10009465, 10009466), als);
        REQUIRE(s == StatusCode::INVALID);
    }
}

TEST_CASE("unified_sites") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    SECTION("trio1") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(unifier_config(), als, sites);
        REQUIRE(s.ok());

        vector<pair<string,size_t> > contigs;
        REQUIRE(data->contigs(contigs).ok());

        REQUIRE(sites.size() == 6);

        REQUIRE(sites[0].pos == range(0,1000,1001));
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[0] == "A");
        REQUIRE(sites[0].alleles[1] == "G");
        REQUIRE(sites[0].unification[allele(range(0,1000,1001),string("A"))] == 0);
        REQUIRE(sites[0].unification[allele(range(0,1000,1001),string("G"))] == 1);
        REQUIRE(sites[0].copy_number.size() == 2);
        REQUIRE(sites[0].copy_number[0] == 0);
        REQUIRE(sites[0].copy_number[1] == 4);

        REQUIRE(sites[1].pos == range(0,1001,1002));
        REQUIRE(sites[1].alleles.size() == 3);
        REQUIRE(sites[1].alleles[0] == "C");
        REQUIRE(sites[1].alleles[1] == "G");
        REQUIRE(sites[1].alleles[2] == "T");
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("C"))] == 0);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("G"))] == 1);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("T"))] == 2);
        REQUIRE(sites[1].copy_number.size() == 3);
        REQUIRE(sites[1].copy_number[0] == 0);
        REQUIRE(sites[1].copy_number[1] == 2);
        REQUIRE(sites[1].copy_number[2] == 2);

        REQUIRE(sites[2].pos == range(0,1010,1012));
        REQUIRE(sites[2].alleles.size() == 2);
        REQUIRE(sites[2].alleles[0] == "CC");
        REQUIRE(sites[2].alleles[1] == "AG");
        REQUIRE(sites[2].unification[allele(range(0,1010,1012),string("CC"))] == 0);
        REQUIRE(sites[2].unification[allele(range(0,1010,1012),string("AG"))] == 1);
        REQUIRE(sites[2].copy_number.size() == 2);
        REQUIRE(sites[2].copy_number[0] == 0);
        REQUIRE(sites[2].copy_number[1] == 3);

        REQUIRE(sites[3].pos == range(0,1100,1101));
        REQUIRE(sites[3].alleles.size() == 2);
        REQUIRE(sites[3].alleles[0] == "C");
        REQUIRE(sites[3].alleles[1] == "A");
        REQUIRE(sites[3].unification[allele(range(0,1100,1101),string("C"))] == 0);
        REQUIRE(sites[3].unification[allele(range(0,1100,1101),string("A"))] == 1);
        REQUIRE(sites[3].copy_number.size() == 2);
        REQUIRE(sites[3].copy_number[0] == 0);
        REQUIRE(sites[3].copy_number[1] == 3);

        REQUIRE(sites[4].pos == range(0,1102,1103));
        REQUIRE(sites[4].alleles.size() == 2);
        REQUIRE(sites[4].alleles[0] == "C");
        REQUIRE(sites[4].alleles[1] == "G");
        REQUIRE(sites[4].unification[allele(range(0,1102,1103),string("C"))] == 0);
        REQUIRE(sites[4].unification[allele(range(0,1102,1103),string("G"))] == 1);
        REQUIRE(sites[4].copy_number.size() == 2);
        REQUIRE(sites[4].copy_number[0] == 0);
        REQUIRE(sites[4].copy_number[1] == 3);

        REQUIRE(sites[5].pos == range(0,1200,1201));
        REQUIRE(sites[5].alleles.size() == 2);
        REQUIRE(sites[5].alleles[0] == "C");
        REQUIRE(sites[5].alleles[1] == "A");
        REQUIRE(sites[5].unification[allele(range(0,1200,1201),string("C"))] == 0);
        REQUIRE(sites[5].unification[allele(range(0,1200,1201),string("A"))] == 1);
        REQUIRE(sites[5].copy_number.size() == 2);
        REQUIRE(sites[5].copy_number[0] == 0);
        REQUIRE(sites[5].copy_number[1] == 3);
    }

    SECTION("2 trios") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(unifier_config(), als, sites);
        REQUIRE(s.ok());

        vector<pair<string,size_t> > contigs;
        REQUIRE(data->contigs(contigs).ok());

        REQUIRE(sites.size() == 6);

        REQUIRE(sites[0].pos == range(0,1000,1001));
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[0] == "A");
        REQUIRE(sites[0].alleles[1] == "G");
        REQUIRE(sites[0].unification[allele(range(0,1000,1001),string("A"))] == 0);
        REQUIRE(sites[0].unification[allele(range(0,1000,1001),string("G"))] == 1);
        REQUIRE(sites[0].copy_number.size() == 2);
        REQUIRE(sites[0].copy_number[0] == 0);
        REQUIRE(sites[0].copy_number[1] == 6);

        REQUIRE(sites[1].pos == range(0,1001,1002));
        REQUIRE(sites[1].alleles.size() == 4);
        REQUIRE(sites[1].alleles[0] == "C");
        REQUIRE(sites[1].alleles[1] == "A");
        REQUIRE(sites[1].alleles[2] == "G");
        REQUIRE(sites[1].alleles[3] == "T");
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("C"))] == 0);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("A"))] == 1);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("G"))] == 2);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("T"))] == 3);
        REQUIRE(sites[1].copy_number.size() == 4);
        REQUIRE(sites[1].copy_number[0] == 0);
        REQUIRE(sites[1].copy_number[1] == 6);
        REQUIRE(sites[1].copy_number[2] == 2);
        REQUIRE(sites[1].copy_number[3] == 2);

        REQUIRE(sites[2].pos == range(0,1010,1013));
        REQUIRE(sites[2].alleles.size() == 3);
        REQUIRE(sites[2].alleles[0] == "CCC");
        REQUIRE(sites[2].alleles[1] == "AGC");
        REQUIRE(sites[2].alleles[2] == "AGA");
        REQUIRE(sites[2].unification[allele(range(0,1010,1012),string("CC"))] == 0);
        REQUIRE(sites[2].unification[allele(range(0,1010,1013),string("CCC"))] == 0);
        REQUIRE(sites[2].unification[allele(range(0,1010,1012),string("AG"))] == 1);
        REQUIRE(sites[2].unification[allele(range(0,1010,1013),string("AGA"))] == 2);
        REQUIRE(sites[2].copy_number.size() == 3);
        REQUIRE(sites[2].copy_number[0] == 0);
        REQUIRE(sites[2].copy_number[1] == 3);
        REQUIRE(sites[2].copy_number[2] == 2);

        REQUIRE(sites[3].pos == range(0,1100,1101));
        REQUIRE(sites[3].alleles.size() == 2);
        REQUIRE(sites[3].alleles[0] == "C");
        REQUIRE(sites[3].alleles[1] == "A");
        REQUIRE(sites[3].unification[allele(range(0,1100,1101),string("C"))] == 0);
        REQUIRE(sites[3].unification[allele(range(0,1100,1101),string("A"))] == 1);
        REQUIRE(sites[3].copy_number.size() == 2);
        REQUIRE(sites[3].copy_number[0] == 0);
        REQUIRE(sites[3].copy_number[1] == 3);

        REQUIRE(sites[4].pos == range(0,1102,1103));
        REQUIRE(sites[4].alleles.size() == 2);
        REQUIRE(sites[4].alleles[0] == "C");
        REQUIRE(sites[4].alleles[1] == "G");
        REQUIRE(sites[4].unification[allele(range(0,1102,1103),string("C"))] == 0);
        REQUIRE(sites[4].unification[allele(range(0,1102,1103),string("G"))] == 1);
        REQUIRE(sites[4].copy_number.size() == 2);
        REQUIRE(sites[4].copy_number[0] == 0);
        REQUIRE(sites[4].copy_number[1] == 3);

        REQUIRE(sites[5].pos == range(0,1200,1201));
        REQUIRE(sites[5].alleles.size() == 2);
        REQUIRE(sites[5].alleles[0] == "C");
        REQUIRE(sites[5].alleles[1] == "A");
        REQUIRE(sites[5].unification[allele(range(0,1200,1201),string("C"))] == 0);
        REQUIRE(sites[5].unification[allele(range(0,1200,1201),string("A"))] == 1);
        REQUIRE(sites[5].copy_number.size() == 2);
        REQUIRE(sites[5].copy_number[0] == 0);
        REQUIRE(sites[5].copy_number[1] == 3);

        // An allele from trio2 that would collapse sites[3] and sites[4] was
        // pruned.

        /*
        for (const auto& site : sites) {
            cout << site.pos.str(contigs);
            for (unsigned i = 0; i < site.alleles.size(); i++) {
                cout << ' ' << site.alleles[i] << '[' << site.copy_number[i] << ']';
            }
            for (const auto& u : site.unification) {
                UNPAIR(u, p, i)
                UNPAIR(p, pos, suballele)
                cout << ' ' << suballele << '@' << pos << "->" << site.alleles[i];
            }
            cout << endl;
        }
        */
    }
}

TEST_CASE("genotyper placeholder") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    const string tfn("/tmp/GLnexus_unit_tests.bcf");

    SECTION("simulate I/O errors") {
        s = Service::Start(service_config(), *data, *data, svc);
        REQUIRE(s.ok());
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());
        vector<unified_site> sites;
        s = unified_sites(unifier_config(), als, sites);
        REQUIRE(s.ok());

        unique_ptr<SimFailBCFData> faildata;
        bool worked = false;

        for (size_t fail_every = 1; fail_every < 100; fail_every++) {
            s = SimFailBCFData::Open(*data, fail_every, faildata);
            REQUIRE(s.ok());

            s = Service::Start(service_config(), *data, *faildata, svc);
            REQUIRE(s.ok());

            s = svc->genotype_sites(genotyper_config(), string("<ALL>"), sites, tfn);
            if (faildata->failed_once()) {
                worked = true;
                REQUIRE(s == StatusCode::IO_ERROR);
                REQUIRE(s.str() == "IOError: SIM");
            } else {
                REQUIRE(s.ok());
            }
        }

        REQUIRE(worked);
    }

    SECTION("unification with multiple contigs") {
        discovered_alleles als0, als1;

        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als0);
        REQUIRE(s.ok());

        s = svc->discover_alleles("<ALL>", range(1, 0, 1000000), als1);
        REQUIRE(s.ok());

        als.clear();
        REQUIRE(merge_discovered_alleles(als0, als).ok());
        REQUIRE(merge_discovered_alleles(als1, als).ok());

        vector<unified_site> sites;
        s = unified_sites(unifier_config(), als, sites);
        cout << s.str() << endl;
        REQUIRE(s.ok());

        REQUIRE(is_sorted(sites.begin(), sites.end()));
        REQUIRE(sites[0].pos.rid == 0);
        REQUIRE(sites[sites.size()-1].pos.rid == 1);
    }
}

TEST_CASE("gVCF genotyper") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"NA12878D_HiSeqX.21.10009462-10009469.gvcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    const string tfn("/tmp/GLnexus_unit_tests.bcf");

    SECTION("no depth requirement") {
        vector<unified_site> sites;

        unified_site us(range(0,10009461,10009462));
        us.alleles.push_back("T");
        us.alleles.push_back("A");
        us.unification[allele(range(0,10009461,10009462),"T")] = 0;
        us.unification[allele(range(0,10009461,10009462),"A")] = 1;
        us.copy_number = { 1, 1 };
        sites.push_back(us);

        us.pos = range(0,10009462,10009463);
        us.alleles.clear();
        us.alleles.push_back("C");
        us.alleles.push_back("G");
        us.unification.clear();
        us.unification[allele(range(0,10009462,10009463),"C")] = 0;
        us.unification[allele(range(0,10009462,10009463),"G")] = 1;
        us.copy_number = { 1, 1 };
        sites.push_back(us);

        // this site spans two gVCF records and will not be called by the current algorithm.
        us.pos = range(0,10009465,10009467);
        us.alleles.clear();
        us.alleles.push_back("AA");
        us.alleles.push_back("GT");
        us.unification.clear();
        us.unification[allele(range(0,10009465,10009467),"AA")] = 0;
        us.unification[allele(range(0,10009465,10009467),"GT")] = 1;
        us.copy_number = { 1, 1 };
        sites.push_back(us);

        s = svc->genotype_sites(genotyper_config(), string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn);
        REQUIRE(s.ok());
    }

    SECTION("require depth > 12") {
        vector<unified_site> sites;

        unified_site us(range(0,10009463,10009465));
        us.alleles.push_back("TA");
        us.alleles.push_back("T");
        us.unification[allele(range(0,10009463,10009465),"TA")] = 0;
        us.unification[allele(range(0,10009463,10009465),"T")] = 1;
        us.copy_number = { 1, 1 };
        sites.push_back(us);

        us.pos = range(0,10009465,10009466);
        us.alleles.clear();
        us.alleles.push_back("A");
        us.alleles.push_back("G");
        us.unification.clear();
        us.unification[allele(range(0,10009465,10009466),"A")] = 0;
        us.unification[allele(range(0,10009465,10009466),"G")] = 1;
        us.copy_number = { 1, 1 };
        sites.push_back(us);

        us.pos = range(0,10009466,10009467);
        us.alleles.clear();
        us.alleles.push_back("A");
        us.alleles.push_back("C");
        us.unification.clear();
        us.unification[allele(range(0,10009466,10009467),"A")] = 0;
        us.unification[allele(range(0,10009466,10009467),"C")] = 1;
        us.copy_number = { 1, 1 };
        sites.push_back(us);

        genotyper_config cfg;
        cfg.required_dp = 13;

        s = svc->genotype_sites(cfg, string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn);
        REQUIRE(s.ok());
    }


    SECTION("require depth > 9") {
        vector<unified_site> sites;

        unified_site us(range(0,10009463,10009465));
        us.alleles.push_back("TA");
        us.alleles.push_back("T");
        us.unification[allele(range(0,10009463,10009465),"TA")] = 0;
        us.unification[allele(range(0,10009463,10009465),"T")] = 1;
        us.copy_number = { 1, 1 };
        sites.push_back(us);

        genotyper_config cfg;
        cfg.required_dp = 9;
        s = svc->genotype_sites(cfg, string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn);
        REQUIRE(s.ok());
    }
}

TEST_CASE("genotype residuals") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;
    s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
    REQUIRE(s.ok());

    vector<unified_site> sites;
    s = unified_sites(unifier_config(), als, sites);
    REQUIRE(s.ok());

    const string tfn("/tmp/GLnexus_unit_tests.bcf");
    genotyper_config cfg;
    cfg.output_residuals = true;
    s = svc->genotype_sites(cfg, string("<ALL>"), sites, tfn);
    REQUIRE(s.ok());

    // verify that the residuals file is in correct yaml format
    YAML::Node resFile = YAML::LoadFile("/tmp/GLnexus_unit_tests.residuals.yml");
    REQUIRE(resFile.size() == 3);

    // It seems that a list of documents split by --- symbols
    // are parsed as a yaml map.
    REQUIRE(resFile.IsMap());
}
