"""Literature column: a second adjudicator for each rule, independent of the data.

    python -m src.literature export results/run2-20bases results/run3-prompt-v2
    [retrieval -- see below]
    python -m src.literature fetch
    python -m src.literature judge
    python -m src.literature table

Everything lives in results/literature/ and is committed.

export   claims.json: id + claim text only. No data verdicts, so nothing
         downstream can be steered by them.
retrieval  Valency is a claude.ai connector, which headless `claude -p` cannot
         reach, so retrieval runs in an interactive Claude Code session:
         one semantic_search_papers call per claim, query = the claim text
         VERBATIM, limit 5, no filters. The returned ids/titles/urls and the
         request_id go into hits.json unedited. Retrieval involves no
         judgement, which is why it is safe to do from a session that has
         seen the data verdicts.
fetch    abstracts.json: abstracts pulled from PubMed / arXiv by id, cached.
judge    judgements.json: a blind `claude -p` with zero tools reads only the
         claim and its abstracts and returns a stance. Every quoted piece of
         evidence is checked to be a verbatim substring of its abstract;
         paper-level stances whose quote fails the check are discarded.
table    adjudication.md: data verdict next to literature stance, with
         disagreements flagged -- the most interesting rows.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

LIT = Path("results/literature")
FAMILIES = {"baseline.json": "textbook", "rules.json": "llm",
            "rules_prerevision.json": "llm_prerevision"}


def claim_id(text: str) -> str:
    return hashlib.sha1(text.encode()).hexdigest()[:10]


def _run_rows(run_dirs):
    """(run, family, row) for every scored rule in the given run folders."""
    for d in map(Path, run_dirs):
        for fname, family in FAMILIES.items():
            p = d / fname
            if p.exists():
                for row in json.loads(p.read_text()):
                    yield d.name, family, row


def export(run_dirs):
    LIT.mkdir(parents=True, exist_ok=True)
    claims = {}
    for _, _, row in _run_rows(run_dirs):
        claims.setdefault(claim_id(row["text"]), row["text"])
    out = [{"id": k, "text": v} for k, v in claims.items()]
    (LIT / "claims.json").write_text(json.dumps(out, indent=2))
    (LIT / "runs.json").write_text(json.dumps([str(d) for d in run_dirs]))
    print(f"{len(out)} unique claims -> {LIT / 'claims.json'}")


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "osi-falsifying-cf/0.1"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def _pubmed(ids):
    out = {}
    for i in range(0, len(ids), 100):
        xml = _get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                   f"?db=pubmed&retmode=xml&id={','.join(ids[i:i + 100])}")
        for art in ET.fromstring(xml).iter("PubmedArticle"):
            pmid = art.findtext(".//PMID")
            parts = ["".join(t.itertext()) for t in art.iter("AbstractText")]
            title = "".join(art.find(".//ArticleTitle").itertext()).strip()
            out[pmid] = (title, " ".join(p.strip() for p in parts if p.strip()))
    return out


def _arxiv(ids):
    out = {}
    ns = {"a": "http://www.w3.org/2005/Atom"}
    for i in range(0, len(ids), 50):
        xml = _get("http://export.arxiv.org/api/query?max_results=100"
                   f"&id_list={','.join(ids[i:i + 50])}")
        for e in ET.fromstring(xml).findall("a:entry", ns):
            aid = e.findtext("a:id", namespaces=ns).rsplit("/", 1)[-1]
            aid = re.sub(r"v\d+$", "", aid)
            out[aid] = (" ".join(e.findtext("a:title", "", ns).split()),
                        " ".join(e.findtext("a:summary", "", ns).split()))
    return out


def fetch():
    """Cache abstracts by id, and fill each hit's title from the same record
    (hits.json stores only what Valency returned: source, id, rank)."""
    hits_p = LIT / "hits.json"
    hits = json.loads(hits_p.read_text())
    cache_p = LIT / "abstracts.json"
    cache = json.loads(cache_p.read_text()) if cache_p.exists() else {}
    titles_p = LIT / "titles.json"
    titles = json.loads(titles_p.read_text()) if titles_p.exists() else {}
    want = {"pubmed": set(), "arxiv": set()}
    for h in hits.values():
        for p in h["papers"]:
            key = f"{p['source']}:{p['id']}"
            if key not in cache and p["source"] in want:
                want[p["source"]].add(p["id"])
    for src, fn in (("pubmed", _pubmed), ("arxiv", _arxiv)):
        if want[src]:
            for pid, (title, abstract) in fn(sorted(want[src])).items():
                cache[f"{src}:{pid}"] = abstract
                titles[f"{src}:{pid}"] = title
    for h in hits.values():
        for p in h["papers"]:
            p.setdefault("title", titles.get(f"{p['source']}:{p['id']}", ""))
            p.setdefault("url", f"https://pubmed.ncbi.nlm.nih.gov/{p['id']}/"
                         if p["source"] == "pubmed" else f"https://arxiv.org/abs/{p['id']}")
    cache_p.write_text(json.dumps(cache, indent=2))
    titles_p.write_text(json.dumps(titles, indent=2))
    hits_p.write_text(json.dumps(hits, indent=2))
    missing = sum(1 for h in hits.values() for p in h["papers"]
                  if not cache.get(f"{p['source']}:{p['id']}"))
    print(f"{len(cache)} abstracts cached; {missing} hits have no abstract (title only)")


JUDGE_SYSTEM = """You judge whether published literature supports a claim about \
the aqueous solubility of organic molecules. You see the claim and a few \
abstracts returned by a search engine. Judge ONLY from the text of the \
abstracts, not from your own knowledge. Many abstracts will be off-topic; say so.

Output JSON only, no prose, no markdown fences:
{"papers": [{"id": "<id as given>", "relevant": true|false,
             "stance": "supports"|"contradicts"|"neutral",
             "evidence": "<verbatim quote from that abstract, or empty string>"}],
 "overall": "supported"|"contradicted"|"mixed"|"no_evidence",
 "rationale": "<one sentence>"}

- "evidence" must be copied exactly from the abstract, character for character.
- "supports" / "contradicts" need a quote that bears on the claimed direction \
of the effect, not just on the topic.
- "overall" is "no_evidence" unless at least one relevant paper supports or \
contradicts the claim.
"""


def _judge_one(text: str, papers: list[dict], model: str, cid: str) -> dict:
    """One retry on unparseable JSON (verbatim quotes often contain unescaped
    double quotes). Raw responses are kept in raw/ for audit."""
    blocks = [f"[{p['source']}:{p['id']}] {p['title']}\n{p['abstract'] or '(no abstract)'}"
              for p in papers]
    prompt = f"Claim: {text}\n\nAbstracts:\n\n" + "\n\n".join(blocks)
    (LIT / "raw").mkdir(exist_ok=True)
    for attempt in (1, 2):
        txt = _call_judge(prompt, model)
        (LIT / "raw" / f"{cid}_{attempt}.txt").write_text(txt)
        try:
            return json.loads(txt)
        except json.JSONDecodeError as e:
            if attempt == 2:
                raise
            prompt += (f"\n\nYour previous reply was not valid JSON ({e}). Reply again "
                       "with valid JSON only; escape any double quotes inside strings.")


def _call_judge(prompt: str, model: str) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run(
            ["claude", "-p", "--tools", "", "--strict-mcp-config",
             "--mcp-config", '{"mcpServers":{}}', "--system-prompt", JUDGE_SYSTEM,
             "--model", model, "--no-session-persistence", "--output-format", "text"],
            input=prompt, capture_output=True, text=True, cwd=tmp, timeout=600,
        )
    if r.returncode != 0:
        raise RuntimeError(f"claude CLI failed: {r.stderr.strip()[:500]}")
    return r.stdout.strip().removeprefix("```json").removeprefix("```").removesuffix("```")


def _norm(s: str) -> str:
    return " ".join(s.split()).lower()


def judge(model: str = "claude-sonnet-4-6"):
    claims = json.loads((LIT / "claims.json").read_text())
    hits = json.loads((LIT / "hits.json").read_text())
    abstracts = json.loads((LIT / "abstracts.json").read_text())
    out_p = LIT / "judgements.json"
    out = json.loads(out_p.read_text()) if out_p.exists() else {}
    for c in claims:
        if c["id"] in out or c["id"] not in hits:
            continue
        papers = [dict(p, abstract=abstracts.get(f"{p['source']}:{p['id']}", ""))
                  for p in hits[c["id"]]["papers"]]
        j = _judge_one(c["text"], papers, model, c["id"])
        # Quote check: a stance only counts if its quote is really in the abstract.
        by_id = {f"{p['source']}:{p['id']}": p["abstract"] for p in papers}
        by_id.update({p["id"]: p["abstract"] for p in papers})
        for p in j.get("papers", []):
            ok = bool(p.get("evidence")) and _norm(p["evidence"]) in _norm(by_id.get(p["id"], ""))
            p["quote_verified"] = ok
            if p.get("stance") in ("supports", "contradicts") and not ok:
                p["stance_discarded"] = p["stance"]
                p["stance"] = "neutral"
        st = {p["stance"] for p in j.get("papers", []) if p.get("relevant")}
        j["overall_verified"] = ("mixed" if {"supports", "contradicts"} <= st else
                                 "supported" if "supports" in st else
                                 "contradicted" if "contradicts" in st else "no_evidence")
        j["model"] = model
        out[c["id"]] = j
        out_p.write_text(json.dumps(out, indent=2))
        print(f"  {j['overall_verified']:12} (model said {j.get('overall')}) {c['text'][:70]}")


# Claims that explicitly assert an effect beyond lipophilicity (prompt v2
# wording). For these, raw survival is not support: the data must show the
# effect at matched LogP too.
BEYOND_LOGP = re.compile(r"(matched|comparable) MolLogP|independently of lipophilicity|"
                         r"beyond what lipophilicity", re.I)


def _data_says(row) -> str:
    v = row.get("verdict", "")
    if v == "contradicted":
        return "contradicted"
    if BEYOND_LOGP.search(row["text"]):
        return {"survives_beyond_logp": "supported",
                "logp_only": "not supported beyond LogP"}.get(v, "unclear")
    if v in ("survives_beyond_logp", "logp_only", "untestable_beyond_logp"):
        return "supported"
    return "unclear"


def _agree(data: str, lit: str) -> str:
    if lit in ("no_evidence", "mixed", "not judged") or data == "unclear":
        return "n/a"
    if data == "not supported beyond LogP":
        return "DISAGREE" if lit == "supported" else "agree"
    return "agree" if lit == data else "DISAGREE"


def table():
    runs = json.loads((LIT / "runs.json").read_text())
    judg = json.loads((LIT / "judgements.json").read_text())
    hits = json.loads((LIT / "hits.json").read_text())
    lines = ["# Data vs literature", "",
             "Data: verdict on the held-out test split. Literature: blind judge over "
             "the top-5 Valency abstracts per claim, stances kept only with a "
             "verified verbatim quote. **DISAGREE** rows are the interesting ones.", "",
             "| run | family | data verdict | data says | raw | binned | literature | agree? "
             "| claim | papers |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    n = {"agree": 0, "DISAGREE": 0, "n/a": 0}
    for run, family, row in _run_rows(runs):
        cid = claim_id(row["text"])
        j = judg.get(cid)
        lit = j["overall_verified"] if j else "not judged"
        data = _data_says(row)
        agree = _agree(data, lit)
        n[agree] += 1
        cited = [p for p in (j or {}).get("papers", []) if p.get("stance") != "neutral"]
        urls = {p["id"]: p["url"] for p in hits.get(cid, {}).get("papers", [])}
        refs = " ".join(f"[{p['id']}]({urls.get(p['id'].split(':')[-1], '')}) {p['stance']}"
                        for p in cited)
        fmt = lambda x: "" if x is None else f"{x:+.2f}"
        lines.append(f"| {run} | {family} | {row.get('verdict')} | {data} | "
                     f"{fmt(row.get('effect'))} | "
                     f"{fmt(row.get('effect_strat'))} | {lit} | {agree} | "
                     f"{row['text'].replace('|', '/')} | {refs} |")
    lines += ["", f"agree {n['agree']} · DISAGREE {n['DISAGREE']} · n/a {n['n/a']}"]
    (LIT / "adjudication.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-1:]), f"-> {LIT / 'adjudication.md'}")


if __name__ == "__main__":
    cmd, *args = sys.argv[1:]
    {"export": lambda: export(args), "fetch": fetch, "judge": judge, "table": table}[cmd]()
