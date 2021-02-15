#!/bin/bash
perl /home/pelmo/remapping_whole/scripts/remap_api.pl --mode asm-asm --annotation $i --from GCF_000003025.5 --dest GCF_000003025.6 --annot_out $i.annot --report_out $i.report

done
