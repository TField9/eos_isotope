#!/bin/bash

#!/bin/bash
for ((i=1; i<=5000; i+=10)); do
    end=$((i+9))
    if [ $end -gt 5000 ]; then
        end=5000
    fi
    echo "sdst_be7_${i}_${end}.root"
done > /eos/ams/user/s/selu/mdst/tianye/sh/expected_files.txt

ls /eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_train|grep sdst_be7_ > /eos/ams/user/s/selu/mdst/tianye/sh/existing_files.txt

comm -23 <(sort \/eos/ams/user/s/selu/mdst/tianye/sh/expected_files.txt) <(sort /eos/ams/user/s/selu/mdst/tianye/sh/existing_files.txt) > /eos/ams/user/s/selu/mdst/tianye/sh/missing_files.txt