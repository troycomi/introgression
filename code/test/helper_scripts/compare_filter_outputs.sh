#! /bin/bash

actual=/tigress/tcomi/aclark4_temp/results/analysis_test/
expected=/tigress/tcomi/aclark4_temp/results/analysisp4e2/
echo starting comarison of $(basename $actual) to $(basename $expected)

for file in $(ls ${expected}*_filtered1.txt); do
    act=$(echo $file | sed 's/p4e2/_test/g')
    cmp <(sort $act) <(sort $file) \
        && echo $file passed! || echo $file failed #&& exit
done

for file in $(ls ${expected}*_filtered1intermediate.txt); do
    act=$(echo $file | sed 's/p4e2/_test/g')
    cmp <(sort $act | python intermediate_format.py) \
        <(sort $file | python intermediate_format.py) \
        && echo $file passed! || echo $file failed #&& exit
done
