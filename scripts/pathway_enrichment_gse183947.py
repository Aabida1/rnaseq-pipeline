#!/usr/bin/env python3
"""Directional pathway enrichment for GSE183947 limma results.

Uses g:Profiler with all tested genes as the statistical background and the
FDR<0.05 up/down subsets as separate queries. Gene identifiers in this dataset
are gene symbols. Results are exploratory and depend on successful symbol mapping.
"""
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request

INPUT = "cancer-results/limma_tumor_vs_normal_all_genes.csv"
OUTDIR = "cancer-results"
API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "KEGG", "REAC", "WP"]


def read_results(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No differential-expression rows found in {path}")
    for row in rows:
        row["logFC"] = float(row["logFC"])
        row["adj.P.Val"] = float(row["adj.P.Val"])
    return rows


def extract_overlap_genes(term, query):
    # g:Profiler's "intersections" field contains per-query-gene evidence
    # arrays, not gene symbols. Map non-empty evidence entries back to the
    # exact ordered query list submitted to g:Profiler.
    intersections = term.get("intersections")
    if isinstance(intersections, list) and len(intersections) == len(query):
        genes = [gene for gene, evidence in zip(query, intersections) if evidence]
        if genes:
            return ";".join(genes)

    # Some API/client versions expose an explicit intersection list instead.
    intersection = term.get("intersection")
    if isinstance(intersection, list):
        genes = []
        for item in intersection:
            if isinstance(item, str):
                genes.append(item)
            elif isinstance(item, dict):
                value = item.get("name") or item.get("gene") or item.get("id")
                if value:
                    genes.append(str(value))
        if genes:
            return ";".join(dict.fromkeys(genes))
    if isinstance(intersection, str):
        return intersection
    return ""


def run_gprofiler(label, query, background):
    payload = {
        "organism": "hsapiens",
        "query": query,
        "background": background,
        "domain_scope": "custom",
        "sources": SOURCES,
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": False,
    }
    request = urllib.request.Request(
        API,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": "rnaseq-pipeline-GSE183947"},
        method="POST",
    )
    last_error = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                body = json.loads(response.read().decode("utf-8"))
            break
        except urllib.error.HTTPError as exc:
            try:
                detail = exc.read().decode("utf-8", errors="replace")
            except Exception:
                detail = str(exc)
            last_error = RuntimeError(f"HTTP {exc.code}: {detail[:1000]}")
            if exc.code < 500 or attempt == 2:
                break
            time.sleep(3 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(3 * (attempt + 1))
    else:
        body = None

    if last_error is not None:
        raise RuntimeError(f"g:Profiler request failed for {label}: {last_error}")
    if isinstance(body, dict):
        body = body.get("result", [])
    if not isinstance(body, list):
        raise RuntimeError(f"Unexpected g:Profiler response for {label}: {str(body)[:500]}")


    rows = []
    for term in body:
        overlap_genes = extract_overlap_genes(term, query)
        rows.append({
            "direction": label,
            "source": term.get("source", ""),
            "term_id": term.get("native", ""),
            "term_name": term.get("name", ""),
            "p_value": term.get("p_value", ""),
            "intersection_size": term.get("intersection_size", ""),
            "query_size": term.get("query_size", ""),
            "term_size": term.get("term_size", ""),
            "effective_domain_size": term.get("effective_domain_size", ""),
            "intersecting_genes": overlap_genes,
        })
    rows.sort(key=lambda x: float(x["p_value"]) if x["p_value"] != "" else 1.0)
    output = os.path.join(OUTDIR, f"pathway_enrichment_{label}.csv")
    fields = ["direction", "source", "term_id", "term_name", "p_value",
              "intersection_size", "query_size", "term_size",
              "effective_domain_size", "intersecting_genes"]
    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with_overlap = sum(bool(row["intersecting_genes"]) for row in rows)
    print(f"{label}: {len(rows)} enriched terms saved to {output}; "
          f"{with_overlap} terms include overlap gene names")
    return rows


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rows = read_results(INPUT)
    background = sorted({row["gene_id"].strip() for row in rows if row.get("gene_id", "").strip()})
    significant = [row for row in rows if row["adj.P.Val"] < 0.05]
    up = sorted({row["gene_id"].strip() for row in significant if row["logFC"] > 0})
    down = sorted({row["gene_id"].strip() for row in significant if row["logFC"] < 0})
    if not background or not up or not down:
        raise RuntimeError("Background or directional significant gene list is empty.")

    # Keep a local record of the exact input sets and parameters for reproducibility.
    with open(os.path.join(OUTDIR, "pathway_enrichment_inputs.json"), "w", encoding="utf-8") as handle:
        json.dump({
            "dataset": "GSE183947",
            "method": "g:Profiler g:GOSt over-representation analysis",
            "identifier_type": "gene symbols as supplied in GEO matrix",
            "background": "all 20,246 genes tested by limma",
            "significance_filter": "limma adjusted P value < 0.05",
            "direction": "up = positive tumour-vs-normal logFC; down = negative logFC",
            "sources": SOURCES,
            "correction": "g:SCS, threshold 0.05",
            "n_background": len(background),
            "background_note": "Unique gene symbols from the limma result table; custom domain scope",
            "n_significant": len(significant),
            "n_up": len(up),
            "n_down": len(down),
        }, handle, indent=2)

    up_terms = run_gprofiler("up", up, background)
    down_terms = run_gprofiler("down", down, background)
    with open(os.path.join(OUTDIR, "pathway_enrichment_summary.txt"), "w", encoding="utf-8") as handle:
        handle.write("GSE183947 directional pathway enrichment\n")
        handle.write("Method: g:Profiler over-representation; g:SCS adjusted threshold 0.05\n")
        handle.write("Input: gene symbols from limma FPKM-based analysis; background = all tested genes\n")
        handle.write("Interpretation is exploratory; review gene mapping, sample pairing and outliers first.\n\n")
        for label, terms in [("UP IN TUMOUR", up_terms), ("DOWN IN TUMOUR", down_terms)]:
            handle.write(f"{label}: {len(terms)} terms returned\n")
            for term in terms[:10]:
                handle.write(f"  {term['source']} | {term['term_name']} | p={term['p_value']} | "
                             f"overlap={term['intersection_size']} | genes={term['intersecting_genes']}\n")
            handle.write("\n")
    print("Directional pathway enrichment completed.")


if __name__ == "__main__":
    main()
