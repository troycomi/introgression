#! /bin/bash

actual=/tigress/tcomi/aclark4_temp/results/analysis_test/
expected=/tigress/tcomi/aclark4_temp/results/analysisp4e2/
echo starting comarison of $(basename $actual) to $(basename $expected)

for file in $(ls ${expected}*_filtered1.txt); do
    act=$(echo $file | sed 's/p4e2_filtered1.txt/filter1.txt/g')
    cmp <(sort $act) <(sort $file) \
        && echo YAY! $file passed! || echo $file failed #&& exit
done

for file in $(ls ${expected}*_filtered1intermediate.txt); do
    act=$(echo $file | sed 's/p4e2_filtered1intermediate.txt/filter1inter.txt/g')
    cmp <(sort $act | python intermediate_format_1.py) \
        <(sort $file | python intermediate_format_1.py) \
        && echo YAY! $file passed! || echo $file failed #&& exit
done

for file in $(ls ${expected}*_filtered2.txt); do
    act=$(echo $file | sed 's/p4e2_filtered2.txt/filter2.txt/g')
    cmp <(sort $act) <(sort $file) \
        && echo YAY! $file passed! || echo $file failed #&& exit
done

for file in $(ls ${expected}*_filtered2intermediate.txt); do
    act=$(echo $file | sed 's/p4e2_filtered2intermediate.txt/filter2inter.txt/g')
    cmp <(sort $act | python intermediate_format_2.py) \
        <(sort $file | python intermediate_format_2.py) \
        && echo YAY! $file passed! || echo $file failed && exit
done

cmp <(sort ${expected}/filter_2_thresholds_p4e2.txt) \
    <(sort ${expected}/filter2_thresholds.txt) \
    && echo YAY! thresholds passed! || echo thresholds failed
