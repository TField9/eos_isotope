#!/bin/bash
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/he./he3_fda.root
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/he./he4_fda.root
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/li./li6_fda.root
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/li./li7_fda.root
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be7_fda.root
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be9_fda.root
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be10_fda.root
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/b./b10_fda.root
rm -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/b./b11_fda.root


hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/he./he3_fda.root  /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/he./he3_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/he./he4_fda.root  /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/he./he4_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/li./li6_fda.root  /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/li./li6_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/li./li7_fda.root  /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/li./li7_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be7_fda.root  /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be7_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be9_fda.root  /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be9_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be10_fda.root /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/be./be10_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/b./b10_fda.root   /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/b./b10_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/b./b11_fda.root   /eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda/b./b11_*.root