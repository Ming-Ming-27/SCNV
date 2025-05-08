import re

def assess_cnv(
    adata,
    pred_key='annotated_cnvs',
    label_key='simulated_cnvs'
):
    """
    Compare per-cell × per-chromosome predictions vs. labels,
    but if on the same chromosome any predicted segment is correct,
    count the chromosome as correctly predicted.

    Prints overall TP/FP/FN/TN, accuracy, precision, sensitivity.
    """
    # 1. get sorted list of chromosomes
    chroms = sorted(
        adata.var['chromosome'].unique(),
        key=lambda x: (float(x) if x.isdigit() else float('inf'), x)
    )
    # 2. parser: map chrom -> list of (start,end,cn)
    def parse_cnv_list(s):
        out = {}
        if not isinstance(s, str) or not s.strip():
            return out
        for part in s.split(','):
            part = part.strip()
            m = re.match(r'(.+?):(\d+)-(\d+)\s*\(CN\s*(\d+)\)', part)
            if not m:
                continue
            chrom, st, ed, cn = m.groups()
            out.setdefault(chrom, []).append((int(st), int(ed), int(cn)))
        return out

    TP = FP = FN = TN = 0

    for cell in adata.obs_names:
        pred_map  = parse_cnv_list(adata.obs.at[cell, pred_key])
        label_map = parse_cnv_list(adata.obs.at[cell, label_key])

        for chrom in chroms:
            preds = pred_map.get(chrom, [])
            labs  = label_map.get(chrom, [])

            if labs and preds:
                # both have segments: check if any pred matches any label
                match = False
                for (p0, p1, pcn) in preds:
                    for (l0, l1, lcn) in labs:
                        overlap = not (p1 < l0 or l1 < p0)
                        if overlap and pcn == lcn:
                            match = True
                            break
                    if match:
                        break
                if match:
                    TP += 1
                else:
                    FP += 1
            elif preds and not labs:
                # predicted variant but no true variant
                FP += 1
            elif not preds and labs:
                # missed a true variant
                FN += 1
            else:
                # no variant predicted and none truly present
                TN += 1

    total = TP + TN + FP + FN
    accuracy    = (TP + TN) / total if total else float('nan')
    precision   = TP / (TP + FP) if (TP + FP) else float('nan')
    sensitivity = TP / (TP + FN) if (TP + FN) else float('nan')

    print(f"TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
    print(f"Accuracy:    {accuracy:.3f}")
    print(f"Precision:   {precision:.3f}")
    print(f"Sensitivity: {sensitivity:.3f}")

    return {
        'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN,
        'accuracy': accuracy,
        'precision': precision,
        'sensitivity': sensitivity
    }