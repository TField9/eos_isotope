#!/bin/bash
rm -f /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be7_plot.root
rm -f /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be9_plot.root
rm -f /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be10_plot.root
    
hadd -f /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be7_plot.root /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be7_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be9_plot.root /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be9_*.root
hadd -f /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be10_plot.root /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be10_*.root