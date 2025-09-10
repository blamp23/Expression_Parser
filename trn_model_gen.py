#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

# run_parser invokes the Perl script for one rule and captures its output
def run_parser(parser_path: Path, rule: str) -> subprocess.CompletedProcess:
    """Invoke the Perl expression parser with a single rule string."""
    return subprocess.run(
        [str(parser_path), "--rule", rule],
        capture_output=True,
        text=True,
        check=False,
    )


def main(argv=None):
    p = argparse.ArgumentParser(
        description=(
            "Read boolean gene rules (one per line) from a text file, "
            "pass each rule to expressionParser_bl.pl, and output a TSV matrix."
        )
    )
    p.add_argument(
        "input",
        metavar="rules.txt",
        help="Text file with one boolean gene rule per line",
    )
    p.add_argument(
        "--parser",
        default="./expressionParser_bl.pl",
        help="Path to expression parser Perl script (default: ./expressionParser_bl.pl)",
    )
    p.add_argument(
        "--output",
        "-o",
        default="-",
        help="Write matrix TSV to file (default: stdout)",
    )
    args = p.parse_args(argv)

    parser_path = Path(args.parser)
    try:
        parser_path = parser_path.resolve()
    except Exception:
        pass
    if not parser_path.exists():
        print(f"Error: parser not found at {parser_path}", file=sys.stderr)
        return 2

    rules_path = Path(args.input)
    if not rules_path.exists():
        print(f"Error: input file not found: {rules_path}", file=sys.stderr)
        return 2

    out_stream = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8")
    try:
        with open(rules_path, "r", encoding="utf-8") as fh:
            import re

            columns = []  # list[(label, {var: coeff})]
            all_vars = []
            var_index = {}

            def add_var(name: str):
                if name not in var_index:
                    var_index[name] = len(all_vars)
                    all_vars.append(name)

            label_line_re = re.compile(r"^([^\s:]+):\s+(.*\S)\s*$")

            for idx, raw in enumerate(fh, start=1):
                line = raw.strip().rstrip("\r")
                if not line or line.startswith("#"):
                    continue

                # Support "GENE: (rule ...)" by converting first colon to a space
                if ":" in line:
                    head, tail = line.split(":", 1)
                    line = f"{head.strip()} {tail.strip()}"

                proc = run_parser(parser_path, line)
                if proc.returncode != 0:
                    msg = (
                        f"Parser failed (exit {proc.returncode}) for line {idx}: {line}"
                        f"\nSTDERR:\n{proc.stderr or ''}"
                    )
                    print(msg, file=sys.stderr)
                    return proc.returncode

                # Parse equations from parser stdout
                capturing = False
                for out_line in (proc.stdout or "").splitlines():
                    s = out_line.strip()
                    if s == "PRINTING RULES:":
                        capturing = True
                        continue
                    if not capturing:
                        continue
                    m = label_line_re.match(s)
                    if not m:
                        continue
                    label, eq = m.group(1), m.group(2)
                    pairs = re.findall(r"([+-]?\d+)\s+(\S+)", eq)
                    if not pairs:
                        continue
                    mapping = {}
                    for coeff_str, var in pairs:
                        coeff = int(coeff_str)
                        mapping[var] = mapping.get(var, 0) + coeff
                        add_var(var)
                    columns.append((label, mapping))

        # -------------------- NEW: Normalize protein rows --------------------
        # For every protein base seen as pr_<B>, NOT_pr_<B>, or NOT_<B> (and not *_AC),
        # ensure BOTH rows pr_<B> and NOT_pr_<B> exist (zeros if never referenced).
        def protein_base(v: str):
            if v.endswith("_AC"):
                return None
            if v.startswith("NOT_pr_"):
                return v[len("NOT_pr_") :]
            if v.startswith("pr_"):
                return v[len("pr_") :]
            if v.startswith("NOT_"):
                # allow NOT_SREBP1 style tokens to define the base SREBP1
                return v[len("NOT_") :]
            return None

        bases = set()
        for v in all_vars:
            b = protein_base(v)
            if b:
                bases.add(b)

        # Add missing pr_<base> and NOT_pr_<base> rows before building EX_/AV_
        for b in sorted(bases):
            if f"pr_{b}" not in var_index:
                add_var(f"pr_{b}")
            if f"NOT_pr_{b}" not in var_index:
                add_var(f"NOT_pr_{b}")
        # --------------------------------------------------------------------

        # ---------- Build EX_ and AV_ columns (after normalization) ----------
        # Identify gene bases (*_AC) and protein bases for AV_ construction
        gene_bases = set()
        protein_bases = set()
        for v in all_vars:
            if v.endswith("_AC"):
                gene_bases.add(v[:-3])  # strip "_AC"
            else:
                b = protein_base(v)
                if b:
                    protein_bases.add(b)

        # EX_<gene>: -1 on <gene>_AC row
        for g in sorted(gene_bases):
            ac_row = f"{g}_AC"
            if ac_row in var_index:
                columns.append((f"EX_{g}", {ac_row: -1}))

        # AV_<protein>: 1 on pr_<b>, NOT_pr_<b>, and NOT_<b> rows if present
        for b in sorted(protein_bases):
            mapping = {}
            for candidate in (f"pr_{b}", f"NOT_pr_{b}", f"NOT_{b}"):
                if candidate in var_index:
                    mapping[candidate] = 1
            if mapping:
                columns.append((f"AV_{b}", mapping))

        # Emit TSV matrix
        print("variable\t" + "\t".join(label for (label, _) in columns), file=out_stream)
        for var in sorted(all_vars):
            row = [var]
            for _, mapping in columns:
                row.append(str(mapping.get(var, 0)))
            print("\t".join(row), file=out_stream)

        return 0
    finally:
        if out_stream is not sys.stdout:
            out_stream.close()


if __name__ == "__main__":
    raise SystemExit(main())





