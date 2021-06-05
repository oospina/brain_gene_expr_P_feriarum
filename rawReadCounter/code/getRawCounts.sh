#!/bin/bash

## Argument 1 : Text file with I numbers, 1 per line
## Argument 2 : Text file with filepaths of BAM files, 1 per line

rm -f *.tmp

cat $1 | while read line; do
	inumber=$( fgrep $line $2 )
	samtools view $inumber | cut -f1,3 > ${line}_mappings.tmp
done

cat $1 | while read line; do
	cut -f1 ${line}_mappings.tmp > ${line}_mapreads.tmp
	cut -f2 ${line}_mappings.tmp > ${line}_mapcontigs.tmp
done

rm -f *_mappings.tmp

cat $1 | while read line; do
	paste ${line}_mapcontigs.tmp ${line}_mapreads.tmp > ${line}_BAMmappings.txt
done

rm -f *.tmp

cat $1 | while read line; do
        sed -r -i "s/\t/\t\t\t/g" ${line}_BAMmappings.txt
done

rm -f *.tmp

#cat $1 | while read line; do
#	echo "processing ${line}"
#	java -Xmx16g ExtractRawCountsFromCovBed2 ${line}_BAMmappings.txt ${line}_RawCountsFromBAM
#done > getRawCounts.log.txt
