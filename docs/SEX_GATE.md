# Coverage sex gate

`cnv_sex_gate.py` calibrates sex inference independently within each capture
family using the ratio of mean depth in X non-PAR territory to mean autosomal
depth. Reported sex is never overwritten; inferred sex and concordance are
written to separate fields.

## Frozen rule

For each capture family with known male and female controls:

- values at or below the maximum observed known-male X/autosome ratio are
  inferred as `Male`;
- values at or above the minimum observed known-female ratio are inferred as
  `Female`;
- values in the empty interval between those supports are `Ambiguous` and held.

The known labels are validated with leave-one-out midpoint classification before
the gate can be emitted. Overlapping known-sex supports, duplicate sample IDs,
missing assays, invalid ratios, or imperfect leave-one-out concordance fail
closed.

Y/autosome coverage is preserved for audit but is not used by this gate. The
original imported Y compartments used an approximate Y-XTR mask. The frozen
XTR intervals in `resources/sex_compartments.GRCh38.v1.json` therefore require
recomputation from retained interval summaries before any mLOY analysis.

## FULL1329 validation

The private validation receipt records:

- 1,095 known-sex samples and 100% leave-one-out concordance in both capture
  families;
- no reported/inferred discordance;
- 233 of 234 previously unknown samples assigned (193 female, 40 male);
- one sample in the empirical X/autosome gap, retained as `HOLD`.

No CNV caller or cohort analysis is launched by this step.

## Invocation

```bash
python3 cnv_sex_gate.py \
  --cohort /private/input/cohort.PRIVATE.json \
  --mask resources/sex_compartments.GRCh38.v1.json \
  --output-dir /private/output/sex_gate.PRIVATE \
  --expected 1329
```

The command writes a private cohort with separate inferred-sex fields, an
aggregate summary, and the assay-specific threshold/validation receipt.
