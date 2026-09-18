# mLOY interval collection

`collect_mloy_intervals.py` extracts the minimum verified coverage information
needed to construct the male Y callable panel and calculate mLRR-Y. It reads a
completed staged SNV QC directory and never scans BAMs or reruns depth.

The collector independently verifies:

- the cohort-level `QC_COMPLETE.json` and expected membership;
- every per-sample `COMPLETE.json`, fingerprint, and receipt hash;
- every `coverage.intervals.tsv.gz` against its recorded SHA-256;
- a common Y parent-interval catalog across all samples;
- the frozen GRCh38 PAR/XTR resource and coordinate convention.

Because the retained summaries were segmented with an older XTR mask, a new
boundary can fall inside an aggregated segment. The collector never prorates or
interpolates that depth. A whole parent BED interval is excluded if any part
touches PAR or XTR, and the number of additionally excluded bases is reported.

The deterministic compressed output contains autosomal depth totals and Y-only
parent-interval depth summaries. It does not select callable intervals, define a
null distribution, estimate cell fraction, or call an event. Those operations
are gated on review of the collected male cohort.

```bash
python3 collect_mloy_intervals.py \
  --stages /path/to/full1329/stages \
  --mask resources/sex_compartments.GRCh38.v1.json \
  --output /private/path/mloy-y-intervals.PRIVATE.json.gz \
  --expected 1329 \
  --workers 4
```
