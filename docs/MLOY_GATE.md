# mLOY callable-panel gate

`mloy_gate.py` freezes an assay-specific callable-Y panel and calculates the
continuous mLRR-Y screening metric. It does not emit mLOY event calls.

## Callable panel

An eligible non-PAR/non-XTR Y parent interval is selected independently within
each capture family when at least 80% of its territory is at least 20x in at
least 95% of sex-gate-passing male samples. The sample metric is the median of
`log2(interval mean depth / autosomal mean depth)` over the selected intervals.

Samples with positive depth in less than 95% of the panel are held. Complete
absence of Y signal is a manual-review state, not an automatic mLOY event.

## Normalization and young null

The raw metric is centered on TP53-negative males:

- use the batch-specific control median when a batch has at least five eligible
  male controls;
- otherwise use the assay-wide control median.

The screening null is the youngest tertile of eligible TP53-negative males
within assay. The exploratory screen requires mLRR-Y below the young-null median
minus three robust standard deviations and an interpretive dosage loss of at
least 10%.

The screen is not an event call. `event_call_LOY` remains `NOT_CALIBRATED` until
BAM-level downsampling of Y non-PAR reads establishes recovery, bias and
sensitivity by assay and depth.

## FULL1329 gate result

- 273 male samples entered the module: 256 TP53-negative, 14 R337H and three
  non-R337H LFS;
- IDT panel: 160 intervals, 27,772 bases;
- Twist panel: 169 intervals, 28,019 bases;
- 271 samples produced eligible metrics;
- one sample had an incomplete panel and one had no Y signal;
- nine Twist samples screened positive before calibration: eight controls and
  one R337H, with ages spanning 9.2 to 79.7 years;
- no binary mLOY event was emitted.

The very young screening-positive control demonstrates that the statistical
screen alone includes technical variation and cannot be interpreted as mLOY
prevalence.
